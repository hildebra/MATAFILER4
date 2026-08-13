#!/usr/bin/env perl
use strict;
use warnings;

use Digest::SHA;
use File::Basename qw(basename);
use File::Spec;
use Getopt::Long qw(GetOptions Configure);

my $VERSION = '1.00';
my ($staging, $manifest, $published_dir, $pigz, $mode, $help);
my $cores = 1;

Configure(qw(no_auto_abbrev no_ignore_case));
GetOptions(
	'staging=s' => \$staging,
	'manifest=s' => \$manifest,
	'publishedDir=s' => \$published_dir,
	'pigz=s' => \$pigz,
	'cores=i' => \$cores,
	'mode=s' => \$mode,
	'help|h' => \$help,
) or die usage();
if ($help) { print usage(); exit 0; }
die usage('unexpected positional arguments: '.join(' ', @ARGV)) if @ARGV;
die usage('-staging, -manifest, and -mode are required')
	unless defined($staging) && defined($manifest) && defined($mode);
die "Unsupported mode '$mode' (expected prepare or cleanup)\n"
	unless $mode eq 'prepare' || $mode eq 'cleanup';
die "Staged directory does not exist: $staging\n" unless -d $staging;
$staging = File::Spec->canonpath(File::Spec->rel2abs($staging));
$manifest = File::Spec->canonpath(File::Spec->rel2abs($manifest));
die "Shard manifest must be inside the staging directory: $manifest\n"
	unless dirname_is($manifest, $staging);

my $plan = read_manifest($manifest, $staging);
if ($mode eq 'cleanup') {
	die usage('-publishedDir is required in cleanup mode')
		unless defined($published_dir) && length($published_dir);
	cleanup_published_handoff($plan, $staging, $manifest, $published_dir);
	exit 0;
}
die usage('-pigz is required in prepare mode') unless defined($pigz) && length($pigz);
die "-cores must be positive\n" unless $cores =~ /^\d+$/ && $cores > 0;
prepare_handoff($plan, $staging, $manifest, $pigz, $cores, $published_dir);
exit 0;

sub prepare_handoff {
	my ($plan, $directory, $manifest_path, $pigz_path, $threads, $published) = @_;
	if (defined($published) && length($published)) {
		my $published_path = File::Spec->canonpath(
			File::Spec->rel2abs($published));
		if (cleanup_checkpoint_matches($plan, $directory)
				&& outputs_complete($plan, $published_path)) {
			cleanup_published_handoff(
				$plan, $directory, $manifest_path, $published_path);
			return;
		}
	}
	validate_controller_plan($plan, $directory);
	validate_part_inventory($plan, $directory);
	my @temporary;
	my $ok = eval {
		my ($fna_tmp, $fna_digest, $fna_records, $overlay_fna_digest,
			$overlay_fna_records) = prepare_fasta(
			$plan, $directory, 'fna', $pigz_path, $threads);
		push @temporary, $fna_tmp;
		my ($faa_tmp, $faa_digest, $faa_records, $overlay_faa_digest,
			$overlay_faa_records) = prepare_fasta(
			$plan, $directory, 'faa', $pigz_path, $threads);
		push @temporary, $faa_tmp;

		my ($category_tmp, $category_digest, $category_records,
			$category_loci, $category_samples, $category_sample_set,
			$overlay_category_digest, $overlay_category_records) = prepare_category(
			$plan, $directory, $pigz_path, $threads);
		push @temporary, $category_tmp;
		my ($qc_tmp, $qc_samples, $qc_sample_set) = prepare_qc(
			$plan, $directory, $pigz_path, $threads);
		push @temporary, $qc_tmp;
		my ($link_tmp, $link_digest, $link_records) = prepare_links(
			$plan, $directory, $pigz_path, $threads);
		push @temporary, $link_tmp;

		my $expected_records = $plan->{expected_records};
		for my $observed (
			['FNA', $fna_records], ['FAA', $faa_records],
			['category', $category_records], ['links', $link_records],
		) {
			die "$observed->[0] record count $observed->[1] does not match manifest count $expected_records\n"
				unless $observed->[1] == $expected_records;
		}
		for my $observed (
			['FAA', $faa_digest], ['category', $category_digest], ['links', $link_digest],
		) {
			die "Identifier order differs between FNA and $observed->[0] shards\n"
				unless $observed->[1] eq $fna_digest;
		}
		die "Outgroup FASTA overlays contain different identifiers\n"
			unless $overlay_fna_records == $overlay_faa_records
				&& $overlay_fna_digest eq $overlay_faa_digest;
		die "Outgroup category and FASTA overlays contain different identifiers\n"
			unless $overlay_category_records == $overlay_fna_records
				&& $overlay_category_digest eq $overlay_fna_digest;
		die "Final category contains $category_loci loci, expected $plan->{expected_loci}\n"
			unless $category_loci == $plan->{expected_loci};
		my $expected_final_samples = $plan->{expected_ingroup_samples}
			+ ($overlay_fna_records ? 1 : 0);
		die "Final category contains $category_samples samples, expected $expected_final_samples\n"
			unless $category_samples == $expected_final_samples;
		die "QC contains $qc_samples samples, expected $plan->{expected_ingroup_samples}\n"
			unless $qc_samples == $plan->{expected_ingroup_samples};
		for my $sample (keys %{$qc_sample_set}) {
			die "QC sample '$sample' is absent from the staged category\n"
				unless $category_sample_set->{$sample};
		}

		my $data_tmp = temporary_output($directory, $plan->{outputs}{data_log});
		write_gzip($data_tmp, "OG:$plan->{outgroup}\n", $pigz_path, $threads);
		push @temporary, $data_tmp;

		my @types = qw(fna faa category qc link data_log);
		for my $index (0 .. $#types) {
			my $final = File::Spec->catfile($directory,
				$plan->{outputs}{$types[$index]}.'.gz');
			atomic_replace($temporary[$index], $final);
			$temporary[$index] = '';
			my $plain = File::Spec->catfile($directory, $plan->{outputs}{$types[$index]});
			unlink $plain or die "Cannot remove stale staged aggregate $plain: $!\n"
				if -e $plain;
		}

		my $prepared = File::Spec->catfile($directory,
			'.strain_tree_input.prepared.tsv');
		atomic_write($prepared, join("\t", 'strain-staged-input-v1',
			$plan->{outgroup}, $category_loci, $category_samples, $plan->{mgs})."\n");
		1;
	};
	my $error = $@;
	for my $path (@temporary) {
		next unless defined($path) && length($path) && -e $path;
		unlink $path or warn "Cannot remove rejected staged output $path: $!\n";
	}
	die $error unless $ok;
	print "Fused and validated $plan->{expected_records} strain records from "
		.scalar(@{$plan->{workers}})." worker shard set(s) for $plan->{mgs}\n";
}

sub prepare_fasta {
	my ($plan, $directory, $type, $pigz_path, $threads) = @_;
	my (%records, %sort_key);
	my $digest = Digest::SHA->new(256);
	my $record_count = 0;
	for my $worker (@{$plan->{workers}}) {
		my $part = part_for($plan, $worker->{id}, $type);
		my ($count) = read_fasta(File::Spec->catfile($directory,
			$part->{basename}), sub {
				my ($identifier, $sequence) = @_;
				die "Duplicate FASTA identifier '$identifier' in $type shards\n"
					if exists $records{$identifier};
				$records{$identifier} = $sequence;
				$sort_key{$identifier} = fasta_sort_key(
					$identifier, $plan->{separator});
				$digest->add($identifier, "\n");
			});
		die "Worker $worker->{id} $type count $count does not match manifest count $worker->{records}\n"
			unless $count == $worker->{records};
		$record_count += $count;
	}

	my $overlay_path = File::Spec->catfile($directory,
		'.strain_tree_input.outgroup.'.$type);
	my ($overlay_count, $overlay_digest) = (0, Digest::SHA->new(256)->hexdigest);
	if (-s $overlay_path) {
		($overlay_count, $overlay_digest) = read_fasta($overlay_path, sub {
			my ($identifier, $sequence) = @_;
			die "Outgroup FASTA identifier '$identifier' duplicates an ingroup record\n"
				if exists $records{$identifier};
			$records{$identifier} = $sequence;
			$sort_key{$identifier} = fasta_sort_key(
				$identifier, $plan->{separator});
		});
	}
	my $temporary = temporary_output($directory, $plan->{outputs}{$type});
	write_gzip_stream($temporary, $pigz_path, $threads, sub {
		my ($output) = @_;
		for my $identifier (sort {
			$sort_key{$a} cmp $sort_key{$b} || $a cmp $b
		} keys %records) {
			print {$output} ">$identifier\n$records{$identifier}\n"
				or die "Cannot write sorted $type output $temporary: $!\n";
		}
	});
	return ($temporary, $digest->hexdigest, $record_count,
		$overlay_digest, $overlay_count);
}

sub read_fasta {
	my ($path, $consume) = @_;
	open my $input, '<', $path or die "Cannot read FASTA shard $path: $!\n";
	my ($identifier, $sequence) = ('', '');
	my $count = 0;
	my $digest = Digest::SHA->new(256);
	my $flush = sub {
		return unless length($identifier);
		die "FASTA record '$identifier' in $path has no sequence\n"
			unless length($sequence);
		$consume->($identifier, $sequence);
		$digest->add($identifier, "\n");
		$count++;
	};
	while (my $line = <$input>) {
		$line =~ s/[\r\n]+\z//;
		if ($line =~ /^>(\S+)\s*$/) {
			$flush->();
			($identifier, $sequence) = ($1, '');
			next;
		}
		die "Sequence data precedes the first FASTA header in $path\n"
			unless length($identifier);
		$line =~ s/\s+//g;
		$sequence .= $line;
	}
	$flush->();
	close $input or die "Cannot close FASTA shard $path: $!\n";
	return ($count, $digest->hexdigest);
}

sub prepare_category {
	my ($plan, $directory, $pigz_path, $threads) = @_;
	my (%loci, %samples, %ingroup_samples);
	my $digest = Digest::SHA->new(256);
	my $records = 0;
	for my $worker (@{$plan->{workers}}) {
		my $part = part_for($plan, $worker->{id}, 'category');
		open my $input, '<', File::Spec->catfile($directory, $part->{basename})
			or die "Cannot read category shard $part->{basename}: $!\n";
		my $worker_records = 0;
		while (my $line = <$input>) {
			$line =~ s/[\r\n]+\z//;
			next unless length($line);
			my @field = split /\t/, $line, -1;
			die "Malformed category shard row in $part->{basename}: $line\n"
				unless @field >= 4 && $field[0] eq $plan->{mgs}
					&& length($field[1]) && length($field[2]) && length($field[3]);
			if (exists($loci{$field[1]}{$field[2]})
					&& $loci{$field[1]}{$field[2]} ne $field[3]) {
				die "Conflicting category identifiers for $field[1]/$field[2]\n";
			}
			$loci{$field[1]}{$field[2]} = $field[3];
			$samples{$field[2]} = $ingroup_samples{$field[2]} = 1;
			$digest->add($field[3], "\n");
			$worker_records++;
		}
		close $input or die "Cannot close category shard $part->{basename}: $!\n";
		die "Worker $worker->{id} category count $worker_records does not match manifest count $worker->{records}\n"
			unless $worker_records == $worker->{records};
		$records += $worker_records;
	}
	die "Category contains ".scalar(keys %ingroup_samples)
		." ingroup samples, expected $plan->{expected_ingroup_samples}\n"
		unless scalar(keys %ingroup_samples) == $plan->{expected_ingroup_samples};

	my $overlay_digest = Digest::SHA->new(256);
	my $overlay_records = 0;
	my $overlay = File::Spec->catfile($directory,
		'.strain_tree_input.outgroup.cat.tsv');
	if (-s $overlay) {
		open my $input, '<', $overlay
			or die "Cannot read category overlay $overlay: $!\n";
		while (my $line = <$input>) {
			$line =~ s/[\r\n]+\z//;
			next unless length($line);
			my @field = split /\t/, $line, -1;
			die "Malformed category overlay row: $line\n"
				unless @field == 3 && !grep { !length($_) } @field
					&& $field[1] eq $plan->{outgroup};
			die "Outgroup overlay refers to absent locus '$field[0]'\n"
				unless exists $loci{$field[0]};
			die "Duplicate outgroup category entry for locus '$field[0]'\n"
				if exists $loci{$field[0]}{$field[1]};
			$loci{$field[0]}{$field[1]} = $field[2];
			$samples{$field[1]} = 1;
			$overlay_digest->add($field[2], "\n");
			$overlay_records++;
		}
		close $input or die "Cannot close category overlay $overlay: $!\n";
	}
	my $temporary = temporary_output($directory, $plan->{outputs}{category});
	write_gzip_stream($temporary, $pigz_path, $threads, sub {
		my ($output) = @_;
		for my $locus (sort keys %loci) {
			print {$output} join("\t", map { $loci{$locus}{$_} }
				sort keys %{$loci{$locus}}), "\n"
				or die "Cannot write finalized category $temporary: $!\n";
		}
	});
	return ($temporary, $digest->hexdigest, $records,
		scalar(keys %loci), scalar(keys %samples), \%samples,
		$overlay_digest->hexdigest, $overlay_records);
}

sub prepare_qc {
	my ($plan, $directory, $pigz_path, $threads) = @_;
	my %sample;
	for my $worker (@{$plan->{workers}}) {
		my $part = part_for($plan, $worker->{id}, 'qc');
		open my $input, '<', File::Spec->catfile($directory, $part->{basename})
			or die "Cannot read QC shard $part->{basename}: $!\n";
		my $rows = 0;
		while (my $line = <$input>) {
			$line =~ s/[\r\n]+\z//;
			next if $line eq '' || $line =~ /^#/ || $line =~ /^MGS\t/;
			my @field = split /\t/, $line, -1;
			die "Malformed sample-QC shard row in $part->{basename}: $line\n"
				unless @field >= 6 && $field[0] eq $plan->{mgs}
					&& length($field[1])
					&& !grep { $_ !~ /^\d+(?:\.\d+)?$/ } @field[3 .. 5];
			my ($mgs, $sample_id, $status, $ambiguous, $csp, $loci) = @field[0 .. 5];
			die "Duplicate QC sample '$sample_id' across worker shards\n"
				if exists $sample{$sample_id};
			$sample{$sample_id} = [$mgs, $status, 0 + $ambiguous, 0 + $csp, 0 + $loci];
			$rows++;
		}
		close $input or die "Cannot close QC shard $part->{basename}: $!\n";
		die "Worker $worker->{id} QC count $rows does not match manifest count $worker->{rows}\n"
			unless $rows == $worker->{rows};
	}
	my $temporary = temporary_output($directory, $plan->{outputs}{qc});
	write_gzip_stream($temporary, $pigz_path, $threads, sub {
		my ($output) = @_;
		print {$output} join("\t", qw(MGS sample status ambiguous_fraction csp_fraction validated_loci)), "\n";
		for my $sample_id (sort keys %sample) {
			print {$output} join("\t", $sample{$sample_id}[0], $sample_id,
				@{$sample{$sample_id}}[1 .. 4]), "\n"
				or die "Cannot write finalized QC $temporary: $!\n";
		}
	});
	return ($temporary, scalar(keys %sample), { map { $_ => 1 } keys %sample });
}

sub prepare_links {
	my ($plan, $directory, $pigz_path, $threads) = @_;
	my $digest = Digest::SHA->new(256);
	my $records = 0;
	my $temporary = temporary_output($directory, $plan->{outputs}{link});
	write_gzip_stream($temporary, $pigz_path, $threads, sub {
		my ($output) = @_;
		for my $worker (@{$plan->{workers}}) {
			my $part = part_for($plan, $worker->{id}, 'link');
			open my $input, '<', File::Spec->catfile($directory, $part->{basename})
				or die "Cannot read link shard $part->{basename}: $!\n";
			my $worker_records = 0;
			while (my $line = <$input>) {
				$line =~ s/[\r\n]+\z//;
				next unless length($line);
				my ($identifier) = split /\t/, $line, 2;
				die "Malformed link shard row in $part->{basename}: $line\n"
					unless defined($identifier) && length($identifier);
				$digest->add($identifier, "\n");
				print {$output} $line, "\n"
					or die "Cannot write finalized link file $temporary: $!\n";
				$worker_records++;
			}
			close $input or die "Cannot close link shard $part->{basename}: $!\n";
			die "Worker $worker->{id} link count $worker_records does not match manifest count $worker->{records}\n"
				unless $worker_records == $worker->{records};
			$records += $worker_records;
		}
	});
	return ($temporary, $digest->hexdigest, $records);
}

sub read_manifest {
	my ($path, $directory) = @_;
	open my $input, '<', $path or die "Cannot read shard manifest $path: $!\n";
	my $format = <$input> // '';
	$format =~ s/[\r\n]+\z//;
	die "Unsupported shard manifest $path\n" unless $format eq 'strain-shard-input-v1';
	my (%scalar, %output, %worker, %part);
	while (my $line = <$input>) {
		$line =~ s/[\r\n]+\z//;
		next unless length($line);
		my @field = split /\t/, $line, -1;
		my $kind = shift @field;
		if ($kind eq 'value') {
			die "Malformed or duplicate manifest value row: $line\n"
				unless @field == 2 && length($field[0]) && !exists $scalar{$field[0]};
			$scalar{$field[0]} = $field[1];
		} elsif ($kind eq 'output') {
			die "Malformed or duplicate manifest output row: $line\n"
				unless @field == 2 && $field[0] =~ /^(?:fna|faa|link|category|qc|data_log)$/
					&& safe_basename($field[1]) && !exists $output{$field[0]};
			$output{$field[0]} = $field[1];
		} elsif ($kind eq 'worker') {
			die "Malformed or duplicate manifest worker row: $line\n"
				unless @field == 3 && $field[0] =~ /^\d+$/ && $field[1] =~ /^\d+$/
					&& $field[2] =~ /^\d+$/ && !exists $worker{$field[0]};
			$worker{$field[0]} = { id => 0 + $field[0], rows => 0 + $field[1],
				records => 0 + $field[2] };
		} elsif ($kind eq 'part') {
			die "Malformed or duplicate manifest part row: $line\n"
				unless @field == 4 && $field[0] =~ /^\d+$/
					&& $field[1] =~ /^(?:fna|faa|link|category|qc)$/
					&& safe_basename($field[2]) && $field[3] =~ /^\d+$/
					&& !exists $part{$field[0]}{$field[1]};
			$part{$field[0]}{$field[1]} = {
				basename => $field[2], bytes => 0 + $field[3],
			};
		} else {
			die "Unknown shard manifest row type '$kind' in $path\n";
		}
	}
	close $input or die "Cannot close shard manifest $path: $!\n";
	for my $name (qw(mgs outgroup generation separator expected_records
		expected_ingroup_samples expected_loci)) {
		die "Shard manifest is missing value '$name'\n" unless exists $scalar{$name};
	}
	die "Unsafe MGS value in shard manifest\n"
		unless $scalar{mgs} =~ /^[A-Za-z0-9][A-Za-z0-9_.:+-]*$/;
	die "Unsafe outgroup value in shard manifest\n"
		unless $scalar{outgroup} eq ''
			|| $scalar{outgroup} =~ /^[A-Za-z0-9][A-Za-z0-9_.:+-]*$/;
	die "Unsafe generation value in shard manifest\n"
		unless $scalar{generation} =~ /^[A-Za-z0-9_.:-]+$/;
	die "Unsafe separator in shard manifest\n"
		unless length($scalar{separator}) && $scalar{separator} !~ /[\t\r\n]/;
	for my $name (qw(expected_records expected_ingroup_samples expected_loci)) {
		die "Non-numeric manifest value '$name'\n" unless $scalar{$name} =~ /^\d+$/;
		$scalar{$name} = 0 + $scalar{$name};
	}
	for my $type (qw(fna faa link category qc data_log)) {
		die "Shard manifest is missing output '$type'\n" unless exists $output{$type};
	}
	die "Shard manifest has no contributing workers\n" unless keys %worker;
	for my $id (keys %worker) {
		for my $type (qw(fna faa link category qc)) {
			die "Shard manifest is missing $type part for worker $id\n"
				unless exists $part{$id}{$type};
		}
	}
	my @workers = map {
		my $entry = $worker{$_};
		$entry->{parts} = $part{$_};
		$entry;
	} sort { $a <=> $b } keys %worker;
	return { %scalar, outputs => \%output, workers => \@workers };
}

sub validate_controller_plan {
	my ($manifest_plan, $directory) = @_;
	my $path = File::Spec->catfile($directory, '.strain_tree_input.plan.tsv');
	open my $input, '<', $path or die "Cannot read controller handoff plan $path: $!\n";
	my @lines = <$input>;
	close $input or die "Cannot close controller handoff plan $path: $!\n";
	s/[\r\n]+\z// for @lines;
	die "Unsupported controller handoff plan $path\n"
		unless @lines >= 3 && $lines[0] eq 'strain-staged-input-v1';
	my %value;
	for my $line (@lines[1 .. $#lines]) {
		my ($key, $value) = split /\t/, $line, 2;
		die "Malformed controller handoff plan row: $line\n"
			unless defined($key) && defined($value) && !exists $value{$key};
		$value{$key} = $value;
	}
	die "Controller plan MGS does not match shard manifest\n"
		unless ($value{mgs} // '') eq $manifest_plan->{mgs};
	die "Controller plan outgroup does not match shard manifest\n"
		unless ($value{outgroup} // '') eq $manifest_plan->{outgroup};
}

sub validate_part_inventory {
	my ($plan, $directory) = @_;
	for my $worker (@{$plan->{workers}}) {
		for my $part (values %{$worker->{parts}}) {
			my $path = File::Spec->catfile($directory, $part->{basename});
			my @stat = stat($path);
			die "Shard part is missing: $path\n"
				unless @stat && $stat[7] > 0;
			die "Shard part size changed for $path: observed $stat[7], manifest $part->{bytes}\n"
				unless $stat[7] == $part->{bytes};
		}
	}
	my %expected = map {
		map { $_->{basename} => 1 } values %{$_->{parts}}
	} @{$plan->{workers}};
	opendir my $handle, $directory or die "Cannot read staging directory $directory: $!\n";
	while (my $name = readdir $handle) {
		next unless $name =~ /\.\d+$/;
		next unless grep { $name =~ /^\Q$plan->{outputs}{$_}\E(?:\.tmp)?\.\d+$/ }
			qw(fna faa link category qc);
		die "Unexpected worker shard is not covered by the manifest: $name\n"
			unless $expected{$name};
	}
	closedir $handle or die "Cannot close staging directory $directory: $!\n";
}

sub cleanup_checkpoint_matches {
	my ($plan, $directory) = @_;
	my $path = File::Spec->catfile(
		$directory, '.strain_tree_input.cleanup.tsv');
	return 0 unless -s $path;
	open my $input, '<', $path or return 0;
	my $line = <$input> // '';
	my $closed = close $input;
	return 0 unless $closed;
	$line =~ s/[\r\n]+\z//;
	return $line eq join("\t", 'strain-shard-cleanup-v1',
		$plan->{generation}, $plan->{mgs}, $plan->{outgroup});
}

sub cleanup_published_handoff {
	my ($plan, $directory, $manifest_path, $published) = @_;
	$published = File::Spec->canonpath(File::Spec->rel2abs($published));
	die "Published output directory does not exist: $published\n" unless -d $published;
	die "Cannot clean shard handoff before every finalized output is published\n"
		unless outputs_complete($plan, $published);
	my $checkpoint = File::Spec->catfile(
		$directory, '.strain_tree_input.cleanup.tsv');
	atomic_write($checkpoint, join("\t", 'strain-shard-cleanup-v1',
		$plan->{generation}, $plan->{mgs}, $plan->{outgroup})."\n");
	my @remove = map {
		map { File::Spec->catfile($directory, $_->{basename}) } values %{$_->{parts}}
	} @{$plan->{workers}};
	push @remove, map { File::Spec->catfile($directory, $_) } qw(
		.strain_tree_input.outgroup.fna
		.strain_tree_input.outgroup.faa
		.strain_tree_input.outgroup.cat.tsv
		.strain_tree_input.prepared.tsv
		.strain_tree_input.plan.tsv
	);
	for my $path (@remove) {
		next unless -e $path;
		unlink $path or die "Cannot remove committed shard input $path: $!\n";
	}
	unlink $manifest_path
		or die "Cannot remove committed shard manifest $manifest_path: $!\n"
		if -e $manifest_path;
	unlink $checkpoint
		or die "Cannot remove committed cleanup checkpoint $checkpoint: $!\n"
		if -e $checkpoint;
	print "Committed and cleaned worker-shard handoff for $plan->{mgs}\n";
}

sub outputs_complete {
	my ($plan, $directory) = @_;
	return 0 unless defined($directory) && -d $directory;
	for my $type (qw(fna faa link category qc data_log)) {
		my $plain = File::Spec->catfile($directory, $plan->{outputs}{$type});
		my $gzip = "$plain.gz";
		return 0 unless (-s $plain || -s $gzip);
	}
	return 1;
}

sub part_for {
	my ($plan, $worker_id, $type) = @_;
	for my $worker (@{$plan->{workers}}) {
		return $worker->{parts}{$type} if $worker->{id} == $worker_id;
	}
	die "No $type part for worker $worker_id\n";
}

sub fasta_sort_key {
	my ($identifier, $separator) = @_;
	my @part = split /\Q$separator\E/, $identifier, -1;
	die "Malformed strain FASTA identifier '$identifier'\n"
		unless @part == 3 && !grep { !length($_) } @part;
	return join("\t", $part[1], $part[2], $part[0], $identifier);
}

sub temporary_output {
	my ($directory, $basename) = @_;
	return File::Spec->catfile($directory, "$basename.gz.prepare.$$");
}

sub write_gzip {
	my ($path, $contents, $pigz_path, $threads) = @_;
	write_gzip_stream($path, $pigz_path, $threads, sub {
		my ($output) = @_;
		print {$output} $contents or die "Cannot write compressed output $path: $!\n";
	});
}

sub write_gzip_stream {
	my ($path, $pigz_path, $threads, $writer) = @_;
	unlink $path or die "Cannot clear stale temporary output $path: $!\n" if -e $path;
	open my $destination, '>', $path or die "Cannot create compressed output $path: $!\n";
	binmode $destination;
	my $pid = open(my $pipe, '|-');
	die "Cannot fork compressor $pigz_path: $!\n" unless defined $pid;
	if ($pid == 0) {
		open STDOUT, '>&', $destination or die "Cannot connect compressor output: $!\n";
		exec {$pigz_path} $pigz_path, '-p', $threads, '-c', '--';
		die "Cannot execute compressor $pigz_path: $!\n";
	}
	close $destination or die "Cannot close compressed destination $path: $!\n";
	my $ok = eval { $writer->($pipe); 1 };
	my $error = $@;
	my $closed = close $pipe;
	if (!$ok || !$closed) {
		unlink $path if -e $path;
		die $error if !$ok;
		die "Compressor failed while writing $path (status $?)\n";
	}
	die "Compressor produced an empty output: $path\n" unless -s $path;
}

sub atomic_write {
	my ($path, $contents) = @_;
	my $temporary = "$path.write.$$";
	open my $output, '>', $temporary or die "Cannot create $temporary: $!\n";
	print {$output} $contents or die "Cannot write $temporary: $!\n";
	close $output or die "Cannot close $temporary: $!\n";
	atomic_replace($temporary, $path);
}

sub atomic_replace {
	my ($source, $destination) = @_;
	rename $source, $destination
		or die "Cannot atomically replace $destination with $source: $!\n";
}

sub safe_basename {
	my ($value) = @_;
	return defined($value) && length($value) && basename($value) eq $value
		&& $value !~ /[\t\r\n]/ && $value ne File::Spec->curdir
		&& $value ne File::Spec->updir;
}

sub dirname_is {
	my ($path, $directory) = @_;
	my ($volume, $directories) = File::Spec->splitpath($path);
	my $parent = File::Spec->canonpath(File::Spec->catpath($volume, $directories, ''));
	$parent =~ s{[\\/]+\z}{};
	my $expected = $directory;
	$expected =~ s{[\\/]+\z}{};
	return $parent eq $expected;
}

sub usage {
	my ($error) = @_;
	my $text = <<'USAGE';
Usage: finalize_strain_tree_inputs.pl -staging DIR -manifest FILE
       -mode prepare -pigz FILE [-cores INT] [-publishedDir DIR]
   or: finalize_strain_tree_inputs.pl -staging DIR -manifest FILE
       -mode cleanup -publishedDir DIR
USAGE
	return defined($error) ? "$error\n$text" : $text;
}
