package Mods::StrainParts;

use strict;
use warnings;

use Exporter qw(import);
use File::Basename qw(dirname);
use File::Glob qw(bsd_glob);
use File::Path qw(make_path);
use Mods::GenoMetaAss qw(gzipopen gzipwrite);

our @EXPORT_OK = qw(
	balance_assembly_groups
	exact_worker_parts
	write_split_generation
	write_worker_completion
	split_generation_complete
	clear_split_generation
	resolve_fasta_artifact
	append_fasta_records_atomic
	sort_fasta_by_locus
);

sub balance_assembly_groups {
	my ($samples_by_group, $worker_count) = @_;
	die "assembly-group sample map must be a hash reference\n"
		unless ref($samples_by_group) eq 'HASH';
	die "worker count must be positive\n"
		unless defined($worker_count) && $worker_count =~ /^\d+$/ && $worker_count > 0;

	for my $group (keys %{$samples_by_group}) {
		die "assembly-group '$group' samples must be an array reference\n"
			unless ref($samples_by_group->{$group}) eq 'ARRAY';
	}
	my @worker_load = (0) x $worker_count;
	my %worker_for_group;
	for my $group (sort {
		scalar(@{$samples_by_group->{$b}}) <=> scalar(@{$samples_by_group->{$a}})
			|| $a cmp $b
	} keys %{$samples_by_group}) {
		my ($worker) = sort {
			$worker_load[$a] <=> $worker_load[$b] || $a <=> $b
		} 0 .. $worker_count - 1;
		$worker_for_group{$group} = $worker;
		# One unit represents the fixed group-level work; the sample count
		# represents sample-specific VCF/depth consensus work.
		$worker_load[$worker] += 1 + scalar(@{$samples_by_group->{$group}});
	}
	return (\%worker_for_group, \@worker_load);
}

sub exact_worker_parts {
	my ($prefix, $worker_count) = @_;
	die "worker-part prefix is required\n" unless defined($prefix) && length($prefix);
	die "worker count must be positive\n"
		unless defined($worker_count) && $worker_count =~ /^\d+$/ && $worker_count > 0;

	my @parts;
	for my $candidate (bsd_glob("$prefix.*")) {
		next unless $candidate =~ /^\Q$prefix\E\.(\d+)$/;
		my $worker = 0 + $1;
		next if $worker >= $worker_count;
		push @parts, [$worker, $candidate];
	}
	return map { $_->[1] }
		sort { $a->[0] <=> $b->[0] || $a->[1] cmp $b->[1] } @parts;
}

sub _atomic_write {
	my ($path, $contents) = @_;
	my $parent = dirname($path);
	make_path($parent) unless -d $parent;
	my $partial = "$path.part.$$";
	open my $fh, '>', $partial or die "Cannot create $partial: $!\n";
	print {$fh} $contents or die "Cannot write $partial: $!\n";
	close $fh or die "Cannot close $partial: $!\n";
	unless (rename $partial, $path) {
		my $rename_error = $!;
		if (-e $path) {
			unlink $path or die "Cannot replace $path: $!\n";
			rename $partial, $path or die "Cannot publish $path: $!\n";
		} else {
			die "Cannot publish $path: $rename_error\n";
		}
	}
	return $path;
}

sub write_split_generation {
	my ($manifest, $generation, $worker_count) = @_;
	die "split-generation manifest path is required\n"
		unless defined($manifest) && length($manifest);
	die "unsafe split-generation identifier\n"
		unless defined($generation) && $generation =~ /^[A-Za-z0-9_.:-]+$/;
	die "worker count must be positive\n"
		unless defined($worker_count) && $worker_count =~ /^\d+$/ && $worker_count > 0;
	return _atomic_write($manifest, "$generation\t$worker_count\n");
}

sub write_worker_completion {
	my ($stone, $generation) = @_;
	die "worker-completion path is required\n" unless defined($stone) && length($stone);
	die "unsafe split-generation identifier\n"
		unless defined($generation) && $generation =~ /^[A-Za-z0-9_.:-]+$/;
	return _atomic_write($stone, "$generation\n");
}

sub _read_first_line {
	my ($path) = @_;
	return unless -s $path;
	open my $fh, '<', $path or return;
	my $line = <$fh>;
	close $fh or return;
	return unless defined $line;
	chomp $line;
	return $line;
}

sub split_generation_complete {
	my ($manifest, $stone_prefix, $expected_workers) = @_;
	return 0 unless defined($manifest) && defined($stone_prefix);
	return 0 unless defined($expected_workers) && $expected_workers =~ /^\d+$/ && $expected_workers > 0;
	my $line = _read_first_line($manifest);
	return 0 unless defined($line) && $line =~ /^([A-Za-z0-9_.:-]+)\t(\d+)$/;
	my ($generation, $manifest_workers) = ($1, 0 + $2);
	return 0 unless $manifest_workers == $expected_workers;
	for my $worker (0 .. $expected_workers - 1) {
		my $completed_generation = _read_first_line("$stone_prefix.$worker.stone");
		return 0 unless defined($completed_generation) && $completed_generation eq $generation;
	}
	return 1;
}

sub clear_split_generation {
	my ($manifest, $stone_prefix) = @_;
	for my $path ($manifest, bsd_glob("$stone_prefix.*.stone")) {
		next unless defined($path) && length($path);
		next if $path ne $manifest && $path !~ /^\Q$stone_prefix\E\.\d+\.stone$/;
		unlink $path or die "Cannot remove split-generation state $path: $!\n"
			if -f $path || -l $path;
	}
}

sub resolve_fasta_artifact {
	my ($nominal_path) = @_;
	die "nominal FASTA path is required\n"
		unless defined($nominal_path) && length($nominal_path);
	my $plain_exists = -e $nominal_path;
	my $gzip_path = "$nominal_path.gz";
	my $gzip_exists = -e $gzip_path;
	die "Ambiguous FASTA sidecars exist at $nominal_path and $gzip_path; refusing to choose one\n"
		if $plain_exists && $gzip_exists;
	return $plain_exists ? $nominal_path : $gzip_exists ? $gzip_path : '';
}

sub append_fasta_records_atomic {
	my ($nominal_path, $records) = @_;
	$records = '' unless defined $records;
	return $nominal_path unless length $records;

	my $source = resolve_fasta_artifact($nominal_path);
	die "Cannot append to missing FASTA ${nominal_path}[.gz]\n" unless length $source;
	my $compressed = $source eq "$nominal_path.gz";
	my $partial = "$source.rewrite.$$" . ($compressed ? '.gz' : '');

	my ($in, $out);
	if ($compressed) {
		my $ok;
		($in, $ok) = gzipopen($source, 'gzip FASTA', 1, 0);
		die "Cannot read gzip FASTA $source\n" unless $ok && defined $in;
		$out = gzipwrite($partial, 'rewritten gzip FASTA', { threads => 1 });
	} else {
		open $in, '<', $source or die "Cannot read FASTA $source: $!\n";
		open $out, '>', $partial or die "Cannot create FASTA $partial: $!\n";
	}
	binmode $in;
	binmode $out;
	my $last = '';
	my $buffer;
	while (1) {
		my $read = read($in, $buffer, 1024 * 1024);
		die "Cannot read FASTA $source: $!\n" unless defined $read;
		last unless $read;
		print {$out} $buffer or die "Cannot write FASTA $partial: $!\n";
		$last = substr($buffer, -1, 1);
	}
	print {$out} "\n" if length($last) && $last ne "\n";
	print {$out} $records or die "Cannot append FASTA records to $partial: $!\n";
	close $in or die "Cannot close FASTA $source: $!\n";
	close $out or die "Cannot close FASTA $partial: $!\n";

	unless (rename $partial, $source) {
		my $rename_error = $!;
		unlink $partial if -e $partial;
		die "Cannot atomically replace FASTA $source: $rename_error\n";
	}
	return $source;
}

sub sort_fasta_by_locus {
	my ($path, $separator) = @_;
	die "FASTA path is required for locus sorting\n"
		unless defined($path) && length($path);
	die "FASTA locus separator is required\n"
		unless defined($separator) && length($separator);
	die "Cannot locus-sort missing or empty FASTA $path\n" unless -s $path;
	die "sort_fasta_by_locus expects a plain first-generation FASTA, not $path\n"
		if $path =~ /\.gz$/;

	open my $in, '<', $path or die "Cannot read FASTA $path for locus sorting: $!\n";
	binmode $in;
	my (@records, $record_start, $record_key);
	while (1) {
		my $line_start = tell($in);
		my $line = <$in>;
		last unless defined $line;
		next unless $line =~ /^>(\S+)/;
		push @records, [$record_key, $record_start, $line_start - $record_start]
			if defined $record_start;
		my $identifier = $1;
		my @parts = split /\Q$separator\E/, $identifier, -1;
		die "Cannot locus-sort FASTA header '$identifier': expected sample${separator}eggNOG${separator}gene_catalog_id\n"
			unless @parts >= 3 && length($parts[0]) && length($parts[1]) && length($parts[2]);
		my ($sample, $eggnog, $gene_id) = @parts[0, 1, 2];
		$record_key = join("\t", $eggnog, $gene_id, $sample, $identifier);
		$record_start = $line_start;
	}
	my $file_end = tell($in);
	push @records, [$record_key, $record_start, $file_end - $record_start]
		if defined $record_start;
	die "Cannot locus-sort FASTA without records: $path\n" unless @records;

	my $partial = "$path.sort.$$";
	open my $out, '>', $partial or die "Cannot create sorted FASTA $partial: $!\n";
	binmode $out;
	for my $record (sort {
		$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1]
	} @records) {
		seek($in, $record->[1], 0) or die "Cannot seek in FASTA $path: $!\n";
		my $remaining = $record->[2];
		while ($remaining > 0) {
			my $wanted = $remaining > 1024 * 1024 ? 1024 * 1024 : $remaining;
			my $read = read($in, my $buffer, $wanted);
			die "Cannot read FASTA record from $path: $!\n"
				unless defined($read) && $read > 0;
			print {$out} $buffer or die "Cannot write sorted FASTA $partial: $!\n";
			$remaining -= $read;
		}
	}
	close $in or die "Cannot close FASTA $path: $!\n";
	close $out or die "Cannot close sorted FASTA $partial: $!\n";
	unless (rename $partial, $path) {
		my $rename_error = $!;
		unlink $partial if -e $partial;
		die "Cannot atomically publish locus-sorted FASTA $path: $rename_error\n";
	}
	return scalar @records;
}

1;
