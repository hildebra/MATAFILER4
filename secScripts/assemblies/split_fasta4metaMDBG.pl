#!/usr/bin/env perl
use strict;
use warnings;

use constant PI => 4 * atan2(1, 1);
use Getopt::Long qw(GetOptions);
use Mods::GenoMetaAss qw(gzipwrite gzipopen);

# Generate reproducible, assembly-derived synthetic long reads for metaMDBG.
# Read anchors are sampled in proportion to mapped depth. Breakpoints are
# supplied by breakpoints.pl, keeping detection independent of simulation.

my ($fasta_file, $coverage_file, $breakpoint_file, $output_fastq, $length_sd, $help);
my @length_templates;
my $length_sample_size = 100_000;
my $mean_length = 5_000;
my $max_synthetic_depth = 20;
my $seed = 1;
Getopt::Long::Configure(qw(no_auto_abbrev no_ignore_case));
GetOptions(
	'assembly=s' => \$fasta_file,
	'coverage=s' => \$coverage_file,
	'breakpoints=s' => \$breakpoint_file,
	'output=s' => \$output_fastq,
	'mean-read-length=i' => \$mean_length,
	'read-length-sd=f' => \$length_sd,
	'length-template=s@' => \@length_templates,
	'length-sample-size=i' => \$length_sample_size,
	'max-synthetic-depth=f' => \$max_synthetic_depth,
	'seed=i' => \$seed,
	'help|h' => \$help,
) or die usage();
if ($help) {
	print usage();
	exit 0;
}
die usage("unexpected positional arguments: @ARGV") if (@ARGV);
die usage("--assembly, --coverage, --breakpoints and --output are required")
	unless (defined($fasta_file) && defined($coverage_file)
		&& defined($breakpoint_file) && defined($output_fastq));
$length_sd = $mean_length * 0.20 unless (defined $length_sd);

die "Mean synthetic read length must be a positive integer\n"
	unless ($mean_length =~ /^\d+$/ && $mean_length > 0);
die "Synthetic read length SD must be non-negative\n"
	unless ($length_sd =~ /^(?:\d+(?:\.\d*)?|\.\d+)$/ && $length_sd >= 0);
die "Maximum synthetic depth must be positive\n"
	unless ($max_synthetic_depth =~ /^(?:\d+(?:\.\d*)?|\.\d+)$/ && $max_synthetic_depth > 0);
die "Random seed must be an integer\n" unless ($seed =~ /^-?\d+$/);
die "Length sample size must be positive\n" unless $length_sample_size > 0;
if (@length_templates) {
	# Use the observed long-read distribution when real ONT/PacBio reads are
	# available. Sampling bounds startup I/O for very large datasets.
	($mean_length, $length_sd) = empirical_read_lengths(\@length_templates, $length_sample_size);
}
srand($seed);

print "Creating coverage-weighted synthetic long reads for metaMDBG\n"
	."Assembly: $fasta_file\nCoverage: $coverage_file\nOutput: $output_fastq\n"
	."Breakpoints: $breakpoint_file\nMean length: $mean_length; SD: $length_sd; "
	."maximum synthetic depth: $max_synthetic_depth; seed: $seed\n";

my ($lengths, $order) = read_fasta_lengths($fasta_file);
my %intervals = map { $_ => [] } @{$order};
my ($coverage_fh) = gzipopen($coverage_file, 'mapping coverage', 1);
my ($coverage_lines, $matched_lines, $unknown_lines) = (0, 0, 0);
while (my $line = <$coverage_fh>) {
	$line =~ s/[\r\n]+$//;
	next if ($line eq '' || $line =~ /^(?:#|track\b|browser\b)/);
	$coverage_lines++;
	my @fields = split /\s+/, $line;
	if (@fields == 2) {
		my ($raw_id, $depth) = @fields;
		my $id = canonical_id($raw_id);
		validate_depth($depth, $coverage_lines, $line);
		if (!exists $lengths->{$id}) {
			$unknown_lines++;
			next;
		}
		$matched_lines++;
		push @{$intervals{$id}}, [0, $lengths->{$id}, 0 + $depth];
		next;
	}
	die "Malformed mapping coverage line $coverage_lines: $line\n"
		unless (@fields >= 4);
	my ($raw_id, $start, $end, $depth) = @fields[0 .. 3];
	my $id = canonical_id($raw_id);
	die "Invalid mapping coverage coordinates at line $coverage_lines: $line\n"
		unless ($start =~ /^\d+$/ && $end =~ /^\d+$/ && $end > $start);
	validate_depth($depth, $coverage_lines, $line);
	if (!exists $lengths->{$id}) {
		$unknown_lines++;
		next;
	}
	$matched_lines++;
	next if ($start >= $lengths->{$id});
	$end = $lengths->{$id} if ($end > $lengths->{$id});
	push @{$intervals{$id}}, [0 + $start, 0 + $end, 0 + $depth];
}
close $coverage_fh or die "Cannot close mapping coverage $coverage_file: $!\n";
die "Mapping coverage $coverage_file has no usable intervals for the assembly\n"
	unless ($matched_lines);
warn "Ignored $unknown_lines coverage interval(s) for contigs absent from the assembly\n"
	if ($unknown_lines);

# Expand sparse bedGraph data to include implicit zero-depth gaps, then turn
# the accepted breakpoint TSV into independent blocks that reads cannot cross.
my (%coverage_runs, %allowed_blocks);
my ($breakpoint_count, $breakpoint_bases, $breakpoint_contigs) = (0, 0, 0);
my $breakpoints = read_breakpoints($breakpoint_file, $lengths);
for my $id (@{$order}) {
	$coverage_runs{$id} = complete_coverage_runs($id, $lengths->{$id}, $intervals{$id});
	my @breaks = sort { $a->[0] <=> $b->[0] } @{$breakpoints->{$id} || []};
	my (@blocks, $cursor); $cursor = 0;
	for my $break (@breaks) {
		die "Overlapping breakpoint intervals for '$id'\n" if $break->[0] < $cursor;
		push @blocks, [$cursor, $break->[0]] if $break->[0] > $cursor;
		$cursor = $break->[1];
	}
	push @blocks, [$cursor, $lengths->{$id}] if $cursor < $lengths->{$id};
	$allowed_blocks{$id} = \@blocks;
	$breakpoint_count += @breaks;
	$breakpoint_contigs++ if @breaks;
	$breakpoint_bases += $_->[1] - $_->[0] for @breaks;
}

open my $fasta_fh, '<', $fasta_file
	or die "Cannot open assembly '$fasta_file': $!\n";
my $output_fh = gzipwrite($output_fastq, 'metaMDBG synthetic long reads');
my ($header, $sequence) = ('', '');
my ($written_reads, $target_bases, $written_bases, $supported_blocks) = (0, 0, 0, 0);
my ($minimum_read_length, $maximum_read_length);
while (my $line = <$fasta_fh>) {
	$line =~ s/[\r\n]+$//;
	if ($line =~ /^>(.*)$/) {
		if ($header ne '') {
			my ($reads, $target, $bases, $minimum, $maximum, $blocks_used) = simulate_contig(
				$header, $sequence, $coverage_runs{$header}, $allowed_blocks{$header},
				$output_fh, $mean_length, $length_sd, $max_synthetic_depth,
			);
			$written_reads += $reads; $target_bases += $target; $written_bases += $bases;
			$supported_blocks += $blocks_used;
			$minimum_read_length = $minimum if (defined($minimum)
				&& (!defined($minimum_read_length) || $minimum < $minimum_read_length));
			$maximum_read_length = $maximum if (defined($maximum)
				&& (!defined($maximum_read_length) || $maximum > $maximum_read_length));
		}
		$header = canonical_id($1);
		$sequence = '';
		next;
	}
	die "Sequence data found before the first FASTA header in $fasta_file\n"
		if ($header eq '' && $line ne '');
	$sequence .= $line;
}
if ($header ne '') {
	my ($reads, $target, $bases, $minimum, $maximum, $blocks_used) = simulate_contig(
		$header, $sequence, $coverage_runs{$header}, $allowed_blocks{$header},
		$output_fh, $mean_length, $length_sd, $max_synthetic_depth,
	);
	$written_reads += $reads; $target_bases += $target; $written_bases += $bases;
	$supported_blocks += $blocks_used;
	$minimum_read_length = $minimum if (defined($minimum)
		&& (!defined($minimum_read_length) || $minimum < $minimum_read_length));
	$maximum_read_length = $maximum if (defined($maximum)
		&& (!defined($maximum_read_length) || $maximum > $maximum_read_length));
}
close $fasta_fh or die "Cannot close assembly $fasta_file: $!\n";
close $output_fh or die "Cannot close output $output_fastq: $!\n";
die "No synthetic reads were written; mapping coverage was too low across every allowed block\n"
	unless ($written_reads);
my $observed_mean = $written_bases / $written_reads;
my $target_rounded = int($target_bases + 0.5);
my $target_percent = $target_bases > 0 ? 100 * $written_bases / $target_bases : 0;
print "\nSynthetic read simulation summary\n"
	."  Assembly contigs:          ".scalar(@{$order})."\n"
	."  Coverage-supported blocks: $supported_blocks\n"
	."  Breakpoints identified:    $breakpoint_count across $breakpoint_contigs contig(s)\n"
	."  Breakpoint bases excluded: $breakpoint_bases\n"
	."  Simulated reads:           $written_reads\n"
	."  Simulated read bases:      $written_bases\n"
	.sprintf("  Read length mean/min/max: %.1f / %d / %d bp\n",
		$observed_mean, $minimum_read_length, $maximum_read_length)
	."  Target coverage integral:  $target_rounded read-bases\n"
	.sprintf("  Achieved target:          %.1f%%\n", $target_percent)
	."  Output FASTQ:              $output_fastq\n\n";

sub usage {
	my ($error) = @_;
	my $message = '';
	$message .= "Error: $error\n\n" if (defined $error && $error ne '');
	$message .= <<'USAGE';
Usage:
  split_fasta4metaMDBG.pl --assembly FILE --coverage FILE --breakpoints FILE --output FILE [options]

Required:
  --assembly FILE                 Short-read assembly in FASTA format
  --coverage FILE                 Per-base bedGraph coverage, optionally gzip-compressed
  --breakpoints FILE              TSV produced by breakpoints.pl
  --output FILE                   Synthetic FASTQ output (written as gzip)

Simulation options:
  --mean-read-length INT          Requested mean read length [5000]
  --read-length-sd FLOAT          Read-length standard deviation [20% of mean]
  --length-template FILE          FASTQ used to estimate length distribution (repeatable)
  --length-sample-size INT        Maximum template reads inspected [100000]
  --max-synthetic-depth FLOAT     Cap coverage used for simulation [20]
  --seed INT                      Reproducible random seed [1]
  --help, -h                      Show this help
USAGE
	return $message;
}

sub canonical_id {
	my ($value) = @_;
	$value = '' unless (defined $value);
	$value =~ s/^>//;
	$value =~ s/^\s+|\s+$//g;
	$value =~ s/\s.*$//;
	return $value;
}

sub validate_depth {
	my ($depth, $line_number, $line) = @_;
	die "Invalid mapping depth at line $line_number: $line\n"
		unless ($depth =~ /^(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/ && $depth >= 0);
}

sub empirical_read_lengths {
	my ($paths, $limit) = @_;
	my ($count, $mean, $m2) = (0, 0, 0);
	PATH: for my $path (@$paths) {
		my ($fh) = gzipopen($path, 'long-read length template', 1);
		while ($count < $limit) {
			my $header = <$fh>; last unless defined $header;
			my $sequence = <$fh>; my $plus = <$fh>; my $quality = <$fh>;
			die "Truncated FASTQ length template '$path'\n"
				unless defined($sequence) && defined($plus) && defined($quality);
			die "Invalid FASTQ record in length template '$path'\n"
				unless $header =~ /^\@/ && $plus =~ /^\+/;
			$sequence =~ s/[\r\n]+$//; my $length = length($sequence); next unless $length;
			$count++; my $delta = $length - $mean; $mean += $delta / $count;
			$m2 += $delta * ($length - $mean);
		}
		close $fh or die "Cannot close length template '$path': $!\n";
		last PATH if $count >= $limit;
	}
	die "No FASTQ reads found in length templates\n" unless $count;
	my $sd = $count > 1 ? sqrt($m2 / ($count - 1)) : 0;
	return (int($mean + 0.5), $sd);
}

sub read_fasta_lengths {
	my ($path) = @_;
	open my $fh, '<', $path or die "Cannot open assembly '$path': $!\n";
	my (%lengths, @order);
	my $id = '';
	while (my $line = <$fh>) {
		$line =~ s/[\r\n]+$//;
		if ($line =~ /^>(.*)$/) {
			$id = canonical_id($1);
			die "Empty FASTA identifier in $path\n" if ($id eq '');
			die "Duplicate FASTA identifier '$id' in $path\n" if (exists $lengths{$id});
			$lengths{$id} = 0;
			push @order, $id;
			next;
		}
		die "Sequence data found before the first FASTA header in $path\n"
			if ($id eq '' && $line ne '');
		$lengths{$id} += length($line) if ($id ne '');
	}
	close $fh or die "Cannot close assembly $path: $!\n";
	die "No FASTA records found in $path\n" unless (@order);
	for my $contig (@order) {
		die "FASTA record '$contig' has an empty sequence\n" unless ($lengths{$contig});
	}
	return (\%lengths, \@order);
}

sub complete_coverage_runs {
	my ($id, $length, $raw_intervals) = @_;
	# bedGraph may omit positions with no mapped reads. Insert those omitted
	# spans explicitly at depth zero so they can become breakpoint candidates.
	my @sorted = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$raw_intervals};
	my (@runs, $cursor);
	$cursor = 0;
	for my $interval (@sorted) {
		my ($start, $end, $depth) = @{$interval};
		die "Overlapping or unsorted mapping coverage for '$id' at $start-$end\n"
			if ($start < $cursor);
		push @runs, [$cursor, $start, 0] if ($start > $cursor);
		push @runs, [$start, $end, $depth];
		$cursor = $end;
	}
	push @runs, [$cursor, $length, 0] if ($cursor < $length);
	return \@runs;
}

sub read_breakpoints {
	my ($path, $lengths) = @_;
	my ($fh) = gzipopen($path, 'breakpoint TSV', 1);
	my (%breaks, $line_number);
	while (my $line = <$fh>) {
		$line =~ s/[\r\n]+$//;
		next if $line eq '' || $line =~ /^#/;
		$line_number++;
		next if $line_number == 1 && $line =~ /^contig\tstart\tend(?:\t|$)/;
		my ($id, $start, $end) = split /\t/, $line;
		$id = canonical_id($id);
		die "Malformed breakpoint line $line_number: $line\n"
			unless exists($lengths->{$id}) && defined($start) && defined($end)
			&& $start =~ /^\d+$/ && $end =~ /^\d+$/ && $end > $start
			&& $end <= $lengths->{$id};
		push @{$breaks{$id}}, [0+$start, 0+$end];
	}
	close $fh or die "Cannot close breakpoint TSV '$path': $!\n";
	return \%breaks;
}

sub weighted_runs_in_block {
	my ($runs, $block_start, $block_end, $max_depth) = @_;
	# Store cumulative depth integrals. A uniform draw over the final integral
	# can then be translated back to a reference position in O(number of runs).
	my (@weighted, $total);
	$total = 0;
	for my $run (@{$runs}) {
		my ($start, $end, $depth) = @{$run};
		last if ($start >= $block_end);
		next if ($end <= $block_start || $depth <= 0);
		my $overlap_start = $start > $block_start ? $start : $block_start;
		my $overlap_end = $end < $block_end ? $end : $block_end;
		# Cap simulation weight so extreme short-read depth cannot dominate the
		# genuine long-read evidence passed to the hybrid assembler.
		my $weighted_depth = $depth > $max_depth ? $max_depth : $depth;
		my $weight = ($overlap_end - $overlap_start) * $weighted_depth;
		next unless ($weight > 0);
		$total += $weight;
		push @weighted, [$overlap_start, $overlap_end, $total];
	}
	return (\@weighted, $total);
}

sub sample_anchor {
	my ($weighted, $total) = @_;
	my $draw = rand($total);
	for my $run (@{$weighted}) {
		next if ($draw >= $run->[2]);
		return $run->[0] + int(rand($run->[1] - $run->[0]));
	}
	return $weighted->[-1][1] - 1;
}

sub normal_read_length {
	my ($mean, $sd, $block_length) = @_;
	return $block_length if ($block_length <= 1 || $mean >= $block_length);
	my $minimum = int($mean - 3 * $sd + 0.5);
	$minimum = 1 if ($minimum < 1);
	my $maximum = int($mean + 3 * $sd + 0.5);
	$maximum = $block_length if ($maximum > $block_length);
	$minimum = $maximum if ($minimum > $maximum);
	for (1 .. 100) {
		# Box-Muller transform: two uniform random values yield a normal deviate.
		my $u1 = rand(); $u1 = 1e-12 if ($u1 <= 0);
		my $u2 = rand();
		my $length = int($mean + $sd * sqrt(-2 * log($u1)) * cos(2 * PI * $u2) + 0.5);
		return $length if ($length >= $minimum && $length <= $maximum);
	}
	my $fallback = int($mean + 0.5);
	$fallback = $minimum if ($fallback < $minimum);
	$fallback = $maximum if ($fallback > $maximum);
	return $fallback;
}

sub simulate_contig {
	my ($id, $sequence, $runs, $blocks, $fh, $mean, $sd, $max_depth) = @_;
	my ($written, $target_bases, $written_bases, $serial, $blocks_used) = (0, 0, 0, 0, 0);
	my ($minimum_length, $maximum_length);
	for my $block (@{$blocks}) {
		my ($block_start, $block_end) = @{$block};
		my $block_length = $block_end - $block_start;
		next unless ($block_length > 0);
		my ($weighted, $coverage_bases) = weighted_runs_in_block(
			$runs, $block_start, $block_end, $max_depth,
		);
		next unless ($coverage_bases > 0);
		$blocks_used++;
		$target_bases += $coverage_bases;
		# Coverage depth is measured in read-bases per reference base, so its
		# integral divided by a typical read length gives the desired read count.
		my $nominal_length = $mean < $block_length ? $mean : $block_length;
		my $read_count = int($coverage_bases / $nominal_length + 0.5);
		next unless ($read_count > 0);
		for (1 .. $read_count) {
			my $length = normal_read_length($mean, $sd, $block_length);
			$minimum_length = $length if (!defined($minimum_length) || $length < $minimum_length);
			$maximum_length = $length if (!defined($maximum_length) || $length > $maximum_length);
			# Pick a coverage-weighted anchor, then a random offset of that anchor
			# within the read. Clamping to the current block prevents breakpoint crossing.
			my $anchor = sample_anchor($weighted, $coverage_bases);
			my $start = $anchor - int(rand($length));
			$start = $block_start if ($start < $block_start);
			$start = $block_end - $length if ($start + $length > $block_end);
			my $end = $start + $length;
			my $read = substr($sequence, $start, $length);
			my $quality = 'I' x $length;
			$serial++;
			my $read_id = sprintf('%s_SIM_%06d_START_%d_END_%d_ANCHOR_%d',
				$id, $serial, $start, $end, $anchor);
			print {$fh} "\@$read_id\n$read\n+\n$quality\n";
			$written++; $written_bases += $length;
		}
	}
	return ($written, $target_bases, $written_bases,
		$minimum_length, $maximum_length, $blocks_used);
}
