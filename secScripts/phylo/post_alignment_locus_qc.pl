#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use FindBin qw($Bin);
use lib "$Bin/../..";

use Mods::GenoMetaAss qw(gzipopen);

my $VERSION = '0.10';
my %DEFAULT = (
	min_sequences => 3,
	min_occupancy => 0.35,
	min_comparable_fraction => 0.25,
	min_comparable_sites_nt => 30,
	min_comparable_sites_aa => 10,
	max_median_divergence_nt => 0.18,
	max_median_divergence_aa => 0.35,
	max_p90_divergence_nt => 0.30,
	max_p90_divergence_aa => 0.50,
	relative_modified_z => 8.0,
	min_relative_median_divergence_nt => 0.06,
	min_relative_median_divergence_aa => 0.12,
	min_relative_p90_divergence_nt => 0.12,
	min_relative_p90_divergence_aa => 0.25,
	min_loci_for_relative => 8,
);

my ($manifest, $report, $keep, $help);
my $sequence_type = 'nt';
my $min_sequences = $DEFAULT{min_sequences};
my $min_occupancy = $DEFAULT{min_occupancy};
my $min_comparable_fraction = $DEFAULT{min_comparable_fraction};
my $min_comparable_sites;
my $max_median_divergence;
my $max_p90_divergence;
my $relative_modified_z = $DEFAULT{relative_modified_z};
my $min_relative_median_divergence;
my $min_relative_p90_divergence;
my $min_loci_for_relative = $DEFAULT{min_loci_for_relative};

Getopt::Long::Configure(qw(no_auto_abbrev no_ignore_case));
GetOptions(
	'manifest=s' => \$manifest,
	'report=s' => \$report,
	'keep=s' => \$keep,
	'sequenceType=s' => \$sequence_type,
	'minSequences=i' => \$min_sequences,
	'minOccupancy=f' => \$min_occupancy,
	'minComparableFraction=f' => \$min_comparable_fraction,
	'minComparableSites=i' => \$min_comparable_sites,
	'maxMedianDivergence=f' => \$max_median_divergence,
	'maxP90Divergence=f' => \$max_p90_divergence,
	'relativeModifiedZ=f' => \$relative_modified_z,
	'minRelativeMedianDivergence=f' => \$min_relative_median_divergence,
	'minRelativeP90Divergence=f' => \$min_relative_p90_divergence,
	'minLociForRelative=i' => \$min_loci_for_relative,
	'help|h' => \$help,
) or die usage();
if ($help) {
	print usage();
	exit 0;
}
die usage('unexpected positional arguments: '.join(' ', @ARGV)) if @ARGV;
die usage('-manifest, -report, and -keep are required')
	unless defined($manifest) && defined($report) && defined($keep);
$sequence_type = lc($sequence_type);
die "-sequenceType must be nt or aa\n"
	unless $sequence_type eq 'nt' || $sequence_type eq 'aa';

my $type_suffix = $sequence_type eq 'nt' ? 'nt' : 'aa';
$min_comparable_sites //= $DEFAULT{"min_comparable_sites_$type_suffix"};
$max_median_divergence //= $DEFAULT{"max_median_divergence_$type_suffix"};
$max_p90_divergence //= $DEFAULT{"max_p90_divergence_$type_suffix"};
$min_relative_median_divergence //=
	$DEFAULT{"min_relative_median_divergence_$type_suffix"};
$min_relative_p90_divergence //=
	$DEFAULT{"min_relative_p90_divergence_$type_suffix"};

die "Manifest is missing or empty: $manifest\n" unless -s $manifest;
die "-minSequences must be at least 2\n" unless $min_sequences >= 2;
die "-minComparableSites and -minLociForRelative must be positive\n"
	unless $min_comparable_sites > 0 && $min_loci_for_relative > 0;
for my $setting (
	['minOccupancy', $min_occupancy],
	['minComparableFraction', $min_comparable_fraction],
	['maxMedianDivergence', $max_median_divergence],
	['maxP90Divergence', $max_p90_divergence],
	['minRelativeMedianDivergence', $min_relative_median_divergence],
	['minRelativeP90Divergence', $min_relative_p90_divergence],
) {
	die "-$setting->[0] must be between zero and one\n"
		if $setting->[1] < 0 || $setting->[1] > 1;
}
die "-relativeModifiedZ must be non-negative\n" if $relative_modified_z < 0;

open my $manifest_fh, '<', $manifest or die "Cannot open $manifest: $!\n";
my (@paths, %seen);
while (my $line = <$manifest_fh>) {
	$line =~ s/[\r\n]+$//;
	next if $line eq '' || $line =~ /^#/;
	die "Duplicate alignment in manifest: $line\n" if $seen{$line}++;
	push @paths, $line;
}
close $manifest_fh or die "Cannot close $manifest: $!\n";
die "Manifest contains no alignment paths: $manifest\n" unless @paths;

my @loci = map {
	evaluate_alignment($_, $sequence_type, $min_comparable_sites,
		$min_comparable_fraction)
} @paths;

for my $locus (@loci) {
	push @{$locus->{reasons}}, 'too_few_sequences'
		if $locus->{sequences} < $min_sequences;
	push @{$locus->{reasons}}, 'low_occupancy'
		if $locus->{occupancy} < $min_occupancy;
	push @{$locus->{reasons}}, 'insufficient_effective_sequences'
		if $locus->{effective_sequences} < $min_sequences;
	next if @{$locus->{reasons}};
	push @{$locus->{reasons}}, 'high_median_consensus_divergence'
		if $locus->{median_divergence} > $max_median_divergence;
	push @{$locus->{reasons}}, 'high_p90_consensus_divergence'
		if $locus->{p90_divergence} > $max_p90_divergence;
}

my @structurally_valid = grep {
	!grep {
		$_ eq 'malformed_alignment'
			|| $_ eq 'too_few_sequences'
			|| $_ eq 'low_occupancy'
			|| $_ eq 'insufficient_effective_sequences'
	} @{$_->{reasons}}
} @loci;

if (@structurally_valid >= $min_loci_for_relative) {
	my @median_values = map { $_->{median_divergence} } @structurally_valid;
	my @p90_values = map { $_->{p90_divergence} } @structurally_valid;
	my ($median_center, $median_mad) = robust_center_spread(\@median_values);
	my ($p90_center, $p90_mad) = robust_center_spread(\@p90_values);
	for my $locus (@structurally_valid) {
		$locus->{median_modified_z} =
			modified_z($locus->{median_divergence}, $median_center, $median_mad);
		$locus->{p90_modified_z} =
			modified_z($locus->{p90_divergence}, $p90_center, $p90_mad);
		push @{$locus->{reasons}}, 'relative_median_divergence_outlier'
			if $locus->{median_divergence} >= $min_relative_median_divergence
				&& $locus->{median_modified_z} > $relative_modified_z;
		push @{$locus->{reasons}}, 'relative_p90_divergence_outlier'
			if $locus->{p90_divergence} >= $min_relative_p90_divergence
				&& $locus->{p90_modified_z} > $relative_modified_z;
	}
}

write_outputs($report, $keep, \@loci);
my $passed = grep { !@{$_->{reasons}} } @loci;
print "Post-alignment locus QC v$VERSION: retained $passed/".scalar(@loci)
	." loci; report=$report\n";
exit 0;

sub evaluate_alignment {
	my ($path, $type, $minimum_sites, $minimum_fraction) = @_;
	my $result = {
		path => $path,
		reasons => [],
		sequences => 0,
		length => 0,
		occupancy => 0,
		effective_sequences => 0,
		consensus_sites => 0,
		median_divergence => 0,
		p90_divergence => 0,
		max_divergence => 0,
		median_modified_z => 0,
		p90_modified_z => 0,
	};
	unless (-s $path || -s "$path.gz") {
		push @{$result->{reasons}}, 'malformed_alignment';
		return $result;
	}

	my ($fh, $ok) = gzipopen($path, 'post-alignment locus QC', 0, 0);
	unless ($ok) {
		push @{$result->{reasons}}, 'malformed_alignment';
		return $result;
	}
	my (@sequences, $sequence, $saw_header);
	while (my $line = <$fh>) {
		$line =~ s/[\r\n]+$//;
		next if $line eq '';
		if ($line =~ /^>/) {
			push @sequences, uc($sequence) if defined $sequence;
			$sequence = '';
			$saw_header = 1;
			next;
		}
		unless ($saw_header) {
			close $fh;
			push @{$result->{reasons}}, 'malformed_alignment';
			return $result;
		}
		$line =~ s/\s+//g;
		$sequence .= $line;
	}
	push @sequences, uc($sequence) if defined $sequence;
	close $fh;
	$result->{sequences} = scalar(@sequences);
	unless (@sequences && length($sequences[0])) {
		push @{$result->{reasons}}, 'malformed_alignment';
		return $result;
	}
	my $alignment_length = length($sequences[0]);
	if (grep { length($_) != $alignment_length } @sequences) {
		push @{$result->{reasons}}, 'malformed_alignment';
		return $result;
	}
	$result->{length} = $alignment_length;

	my $valid_pattern = $type eq 'nt' ? qr/[ACGTU]/ : qr/[ACDEFGHIKLMNPQRSTVWY]/;
	my (@consensus, $valid_cells);
	$valid_cells = 0;
	for my $position (0 .. $alignment_length - 1) {
		my %counts;
		for my $seq (@sequences) {
			my $character = substr($seq, $position, 1);
			next unless $character =~ $valid_pattern;
			$counts{$character}++;
			$valid_cells++;
		}
		next unless %counts;
		my @ordered = sort { $counts{$b} <=> $counts{$a} || $a cmp $b } keys %counts;
		next if @ordered > 1 && $counts{$ordered[0]} == $counts{$ordered[1]};
		$consensus[$position] = $ordered[0];
		$result->{consensus_sites}++;
	}
	$result->{occupancy} = $valid_cells / (@sequences * $alignment_length);

	my $required_sites = int($result->{consensus_sites} * $minimum_fraction + 0.999999);
	$required_sites = $minimum_sites if $required_sites < $minimum_sites;
	my @divergences;
	for my $seq (@sequences) {
		my ($compared, $differences) = (0, 0);
		for my $position (0 .. $#consensus) {
			next unless defined $consensus[$position];
			my $character = substr($seq, $position, 1);
			next unless $character =~ $valid_pattern;
			$compared++;
			$differences++ if $character ne $consensus[$position];
		}
		next if $compared < $required_sites;
		push @divergences, $differences / $compared;
	}
	$result->{effective_sequences} = scalar(@divergences);
	if (@divergences) {
		$result->{median_divergence} = percentile(\@divergences, 0.5);
		$result->{p90_divergence} = percentile(\@divergences, 0.9);
		$result->{max_divergence} = percentile(\@divergences, 1);
	}
	return $result;
}

sub percentile {
	my ($values, $fraction) = @_;
	return 0 unless @{$values};
	my @sorted = sort { $a <=> $b } @{$values};
	return $sorted[0] if @sorted == 1;
	my $index = $fraction * $#sorted;
	my $lower = int($index);
	my $upper = $lower == $#sorted ? $lower : $lower + 1;
	my $weight = $index - $lower;
	return $sorted[$lower] * (1 - $weight) + $sorted[$upper] * $weight;
}

sub robust_center_spread {
	my ($values) = @_;
	my $center = percentile($values, 0.5);
	my @deviations = map { abs($_ - $center) } @{$values};
	return ($center, percentile(\@deviations, 0.5));
}

sub modified_z {
	my ($value, $center, $mad) = @_;
	return 0 if $value <= $center;
	return 1_000_000 if $mad == 0;
	return 0.67448975 * ($value - $center) / $mad;
}

sub write_outputs {
	my ($report_path, $keep_path, $loci) = @_;
	my $report_tmp = "$report_path.tmp.$$";
	my $keep_tmp = "$keep_path.tmp.$$";
	open my $report_fh, '>', $report_tmp
		or die "Cannot create $report_tmp: $!\n";
	open my $keep_fh, '>', $keep_tmp
		or die "Cannot create $keep_tmp: $!\n";
	print {$report_fh} join("\t", qw(
		alignment status reasons sequences alignment_sites occupancy
		effective_sequences consensus_sites median_consensus_divergence
		p90_consensus_divergence maximum_consensus_divergence
		median_modified_z p90_modified_z
	)), "\n";
	for my $locus (@{$loci}) {
		my $status = @{$locus->{reasons}} ? 'REJECT' : 'PASS';
		print {$report_fh} join("\t",
			$locus->{path}, $status,
			@{$locus->{reasons}} ? join(',', @{$locus->{reasons}}) : '.',
			$locus->{sequences}, $locus->{length},
			sprintf('%.5f', $locus->{occupancy}),
			$locus->{effective_sequences}, $locus->{consensus_sites},
			map { sprintf('%.5f', $_) } @{$locus}{qw(
				median_divergence p90_divergence max_divergence
				median_modified_z p90_modified_z
			)}
		), "\n";
		print {$keep_fh} "$locus->{path}\n" if $status eq 'PASS';
	}
	close $report_fh or die "Cannot close $report_tmp: $!\n";
	close $keep_fh or die "Cannot close $keep_tmp: $!\n";
	rename $report_tmp, $report_path
		or die "Cannot install $report_path: $!\n";
	rename $keep_tmp, $keep_path
		or die "Cannot install $keep_path: $!\n";
}

sub usage {
	my ($error) = @_;
	my $prefix = defined($error) ? "Error: $error\n\n" : '';
	return $prefix.<<"USAGE";
Usage: post_alignment_locus_qc.pl -manifest FILE -report FILE -keep FILE [options]

Single-pass, post-alignment locus QC for metagenomic strain phylogenies. It
measures usable-cell occupancy and per-sequence distance from the unambiguous
column consensus, then applies permissive absolute and cross-locus robust
outlier filters. One divergent strain does not cause rejection unless it drives
the locus-wide 90th percentile beyond the permissive threshold.

  -sequenceType nt|aa             Alignment alphabet [nt]
  -minSequences INT               Comparable sequences [$DEFAULT{min_sequences}]
  -minOccupancy FLOAT             Valid alignment-cell fraction [$DEFAULT{min_occupancy}]
  -minComparableFraction FLOAT    Consensus sites required per sequence [$DEFAULT{min_comparable_fraction}]
  -minComparableSites INT         Absolute sites per sequence [NT $DEFAULT{min_comparable_sites_nt}; AA $DEFAULT{min_comparable_sites_aa}]
  -maxMedianDivergence FLOAT       Absolute consensus distance [NT $DEFAULT{max_median_divergence_nt}; AA $DEFAULT{max_median_divergence_aa}]
  -maxP90Divergence FLOAT          Absolute 90th-percentile distance [NT $DEFAULT{max_p90_divergence_nt}; AA $DEFAULT{max_p90_divergence_aa}]
  -relativeModifiedZ FLOAT         Cross-locus robust outlier threshold [$DEFAULT{relative_modified_z}]
  -minRelativeMedianDivergence FLOAT  Relative-filter floor [NT $DEFAULT{min_relative_median_divergence_nt}; AA $DEFAULT{min_relative_median_divergence_aa}]
  -minRelativeP90Divergence FLOAT  Relative-filter floor [NT $DEFAULT{min_relative_p90_divergence_nt}; AA $DEFAULT{min_relative_p90_divergence_aa}]
  -minLociForRelative INT          Loci required for cross-locus QC [$DEFAULT{min_loci_for_relative}]

The report records every locus and reason. FILE given to -keep contains only
accepted alignment paths and is intended for direct pipeline consumption.
USAGE
}
