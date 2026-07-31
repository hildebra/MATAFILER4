package Mods::StrainSampleStats;

use strict;
use warnings;

use Exporter qw(import);

our @EXPORT_OK = qw(
	sample_stat_columns
	sample_summary_columns
	aggregate_sample_rows
);

my @SAMPLE_STAT_COLUMNS = qw(
	sample worker assembly_driver status selected_mgs candidate_mgs candidate_loci consensus_proteins
	used_mgs skipped_mgs unaccounted_mgs used_fraction
	min_genes_per_mgs presort_genes max_genes qc_enabled min_gene_depth min_bad_loci
	multi_gene_fraction_max csp_gene_fraction_max csp_locus_score_max breakpoint_gene_flank
	abundance_min_loci abundance_min_fold abundance_max_fold abundance_max_modified_z
	capped_mgs capped_loci skipped_within_2_loci_of_min
	retained_loci median_loci_per_used_mgs mean_loci_per_used_mgs
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
	qw(scope samples workers assembly_drivers status_counts processed_samples),
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
	)
);

sub sample_stat_columns {
	return @SAMPLE_STAT_COLUMNS;
}

sub sample_summary_columns {
	return @SUMMARY_COLUMNS;
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

	my (%workers, %drivers, %statuses, %control_values, %summary);
	$summary{scope} = $scope;
	$summary{samples} = scalar(@{$rows});
	$summary{$_} = 0 for @SUM_COLUMNS;

	for my $row (@{$rows}) {
		die "Sample-statistic row must be a hash reference\n" unless ref($row) eq 'HASH';
		my $sample = $row->{sample} // '';
		die "Sample-statistic row has no sample name\n" unless length($sample);
		$workers{$row->{worker} // ''} = 1;
		$drivers{$row->{assembly_driver} // ''} = 1;
		$statuses{$row->{status} // 'unknown'}++;
		for my $column (@CONTROL_COLUMNS) {
			my $value = defined($row->{$column}) ? $row->{$column} : '';
			$control_values{$column}{$value} = 1;
		}
		for my $column (@SUM_COLUMNS) {
			$summary{$column} += _numeric($row->{$column}, $column, $sample);
		}
	}

	$summary{workers} = scalar(keys %workers);
	$summary{assembly_drivers} = scalar(keys %drivers);
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
	return \%summary;
}

1;
