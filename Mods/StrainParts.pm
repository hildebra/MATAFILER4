package Mods::StrainParts;

use strict;
use warnings;

use Exporter qw(import);
use File::Basename qw(dirname);
use File::Glob qw(bsd_glob);
use File::Path qw(make_path);
use IO::Compress::Gzip qw($GzipError);
use IO::Uncompress::Gunzip qw($GunzipError);

our @EXPORT_OK = qw(
	exact_worker_parts
	write_split_generation
	write_worker_completion
	split_generation_complete
	clear_split_generation
	resolve_fasta_artifact
	append_fasta_records_atomic
);

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
	my $partial = "$source.rewrite.$$";

	my ($in, $out);
	if ($compressed) {
		$in = IO::Uncompress::Gunzip->new($source)
			or die "Cannot read gzip FASTA $source: $GunzipError\n";
		$out = IO::Compress::Gzip->new($partial)
			or die "Cannot create gzip FASTA $partial: $GzipError\n";
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

1;
