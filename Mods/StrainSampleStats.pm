package Mods::StrainSampleStats;

use strict;
use warnings;

use Exporter qw(import);

our @EXPORT_OK = qw(
	sample_stat_columns
	sample_summary_columns
	aggregate_sample_rows
	encode_loci_histogram
	loci_histogram_rows
);

my @SAMPLE_STAT_COLUMNS = qw(
	sample worker assembly_group status selected_mgs candidate_mgs candidate_loci consensus_proteins
	used_mgs skipped_mgs unaccounted_mgs used_fraction
	min_genes_per_mgs presort_genes max_genes qc_enabled min_gene_depth min_bad_loci
	multi_gene_fraction_max csp_gene_fraction_max csp_locus_score_max breakpoint_gene_flank
	abundance_min_loci abundance_min_fold abundance_max_fold abundance_max_modified_z
	capped_mgs capped_loci skipped_within_2_loci_of_min
	retained_loci median_loci_per_used_mgs mean_loci_per_used_mgs
	used_mgs_loci_histogram
	pre_abundance_loci post_abundance_loci missing_consensus_loci low_depth_loci
	breakpoint_loci csp_rejected_loci ambiguous_loci abundance_filtered_loci
	invalid_protein_loci placement_flagged_mgs skip_no_selected_loci
	skip_no_usable_loci skip_too_few_after_abundance skip_too_few_valid_sequences
);

my @CONTROL_COLUMNS = qw(
	selected_mgs min_genes_per_mgs presort_genes max_genes qc_enabled
	min_gene_depth min_bad_loci multi_gene_fraction_max csp_gene_fraction_max
	csp_locus_score_max breakpoint_gene_flank abundance_min_loci
	abundance_min_fold abundance_max_fold abundance_max_modified_z
);

my @SUM_COLUMNS = qw(
	candidate_mgs used_mgs skipped_mgs unaccounted_mgs candidate_loci
	consensus_proteins retained_loci capped_mgs capped_loci
	skipped_within_2_loci_of_min pre_abundance_loci post_abundance_loci
	missing_consensus_loci low_depth_loci breakpoint_loci csp_rejected_loci
	ambiguous_loci abundance_filtered_loci invalid_protein_loci
	placement_flagged_mgs skip_no_selected_loci skip_no_usable_loci
	skip_too_few_after_abundance skip_too_few_valid_sequences
);

my @SUMMARY_COLUMNS = (
	qw(scope samples workers assembly_groups status_counts processed_samples),
	@CONTROL_COLUMNS,
	qw(
		candidate_mgs used_mgs skipped_mgs unaccounted_mgs used_fraction
		candidate_loci consensus_proteins retained_loci mean_loci_per_used_mgs
		capped_mgs capped_loci skipped_within_2_loci_of_min
		pre_abundance_loci post_abundance_loci missing_consensus_loci
		low_depth_loci breakpoint_loci csp_rejected_loci ambiguous_loci
		abundance_filtered_loci invalid_protein_loci placement_flagged_mgs
		skip_no_selected_loci skip_no_usable_loci
		skip_too_few_after_abundance skip_too_few_valid_sequences
		used_mgs_loci_histogram
	)
);

my @HISTOGRAM_UPPER_BOUNDS = (9, 19, 49, 99, 199, 499, 999, 1999);

sub sample_stat_columns {
	return @SAMPLE_STAT_COLUMNS;
}

sub sample_summary_columns {
	return @SUMMARY_COLUMNS;
}

sub _histogram_bins {
	my ($minimum) = @_;
	die "Histogram minimum must be a non-negative integer\n"
		unless defined($minimum) && $minimum =~ /\A\d+\z/;
	my $lower = 0 + $minimum;
	my @bins;
	for my $upper (@HISTOGRAM_UPPER_BOUNDS) {
		next if $upper < $lower;
		push @bins, [$lower, $upper, "$lower-$upper"];
		$lower = $upper + 1;
	}
	push @bins, [$lower, undef, "$lower+"];
	return @bins;
}

sub encode_loci_histogram {
	my ($values, $minimum) = @_;
	die "Histogram loci must be an array reference\n" unless ref($values) eq "ARRAY";
	my @bins = _histogram_bins($minimum);
	my %counts = map { $_->[2] => 0 } @bins;
	for my $value (@$values) {
		die "Histogram locus count must be a non-negative integer\n"
			unless defined($value) && $value =~ /\A\d+\z/;
		die "Used-MGS locus count $value is below the configured minimum $minimum\n"
			if $value < $minimum;
		for my $bin (@bins) {
			next if defined($bin->[1]) && $value > $bin->[1];
			$counts{$bin->[2]}++;
			last;
		}
	}
	return join(",", map { $_->[2]."=".$counts{$_->[2]} } @bins);
}

sub _parse_loci_histogram {
	my ($encoded, $minimum, $sample) = @_;
	my @bins = _histogram_bins($minimum);
	my %expected = map { $_->[2] => 1 } @bins;
	my %counts;
	$encoded = "" unless defined($encoded);
	for my $item (split /,/, $encoded, -1) {
		next unless length $item;
		die "Malformed used-MGS locus histogram '$item' for $sample\n"
			unless $item =~ /\A([^=]+)=(\d+)\z/;
		my ($label, $count) = ($1, $2);
		die "Unexpected used-MGS locus histogram bin '$label' for $sample\n"
			unless $expected{$label};
		die "Duplicate used-MGS locus histogram bin '$label' for $sample\n"
			if exists $counts{$label};
		$counts{$label} = 0 + $count;
	}
	for my $bin (@bins) {
		my $label = $bin->[2];
		die "Missing used-MGS locus histogram bin '$label' for $sample\n"
			unless exists $counts{$label};
	}
	return \%counts;
}

sub loci_histogram_rows {
	my ($encoded, $minimum) = @_;
	my $counts = _parse_loci_histogram($encoded, $minimum, "summary");
	return map { [$_->[2], $counts->{$_->[2]}] } _histogram_bins($minimum);
}

sub _numeric {
	my ($value, $column, $sample) = @_;
	$value = 0 unless defined($value) && length($value);
	die "Non-numeric sample statistic '$value' in $column for $sample\n"
		unless $value =~ /\A-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\z/;
	return 0 + $value;
}

sub aggregate_sample_rows {
	my ($rows, $scope) = @_;
	die "Sample-statistic rows must be an array reference\n" unless ref($rows) eq 'ARRAY';
	$scope = 'ALL' unless defined($scope) && length($scope);

	my (%workers, %groups, %statuses, %control_values, %histogram_counts, %summary);
	$summary{scope} = $scope;
	$summary{samples} = scalar(@{$rows});
	$summary{$_} = 0 for @SUM_COLUMNS;

	for my $row (@{$rows}) {
		die "Sample-statistic row must be a hash reference\n" unless ref($row) eq 'HASH';
		my $sample = $row->{sample} // '';
		die "Sample-statistic row has no sample name\n" unless length($sample);
		$workers{$row->{worker} // ''} = 1;
		$groups{$row->{assembly_group} // ''} = 1;
		$statuses{$row->{status} // 'unknown'}++;
		for my $column (@CONTROL_COLUMNS) {
			my $value = defined($row->{$column}) ? $row->{$column} : '';
			$control_values{$column}{$value} = 1;
		}
		for my $column (@SUM_COLUMNS) {
			$summary{$column} += _numeric($row->{$column}, $column, $sample);
		}
		my $minimum = _numeric($row->{min_genes_per_mgs}, "min_genes_per_mgs", $sample);
		die "Non-integer min_genes_per_mgs '$minimum' for $sample\n"
			unless $minimum == int($minimum) && $minimum >= 0;
		my $histogram = _parse_loci_histogram(
			$row->{used_mgs_loci_histogram}, $minimum, $sample
		);
		my $histogram_total = 0;
		for my $label (keys %$histogram) {
			$histogram_counts{$label} += $histogram->{$label};
			$histogram_total += $histogram->{$label};
		}
		my $used_mgs = _numeric($row->{used_mgs}, "used_mgs", $sample);
		die "Used-MGS histogram count $histogram_total does not match used_mgs $used_mgs for $sample\n"
			unless $histogram_total == $used_mgs;
	}

	$summary{workers} = scalar(keys %workers);
	$summary{assembly_groups} = scalar(keys %groups);
	$summary{status_counts} = join(',', map { "$_=$statuses{$_}" } sort keys %statuses);
	$summary{processed_samples} = $statuses{processed} // 0;
	for my $column (@CONTROL_COLUMNS) {
		my @values = sort keys %{$control_values{$column} || {}};
		die "Inconsistent $column values in sample-statistic scope $scope: "
			.join(',', @values)."\n" if @values > 1;
		$summary{$column} = @values ? $values[0] : '';
	}
	$summary{used_fraction} = $summary{candidate_mgs}
		? sprintf('%.6f', $summary{used_mgs} / $summary{candidate_mgs}) : 0;
	$summary{mean_loci_per_used_mgs} = $summary{used_mgs}
		? sprintf('%.3f', $summary{retained_loci} / $summary{used_mgs}) : 0;
	if (@$rows) {
		my @bins = _histogram_bins($summary{min_genes_per_mgs});
		$summary{used_mgs_loci_histogram} = join(",", map {
			$_->[2]."=".($histogram_counts{$_->[2]} // 0)
		} @bins);
	} else {
		$summary{used_mgs_loci_histogram} = "";
	}
	return \%summary;
}

1;
