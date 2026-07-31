use strict;
use warnings;

use FindBin qw($Bin);
use Test::More;

use lib "$Bin/..";
use Mods::StrainSampleStats qw(
	sample_stat_columns sample_summary_columns aggregate_sample_rows
);

my @columns = sample_stat_columns();
is($columns[0], 'sample', 'sample is the first per-sample field');
ok(grep($_ eq 'skip_too_few_valid_sequences', @columns),
	'per-sample fields include the final sequence-filter outcome');

sub sample_row {
	my (%override) = @_;
	my %row = map { $_ => 0 } @columns;
	@row{qw(
		selected_mgs min_genes_per_mgs presort_genes max_genes qc_enabled
		min_gene_depth min_bad_loci multi_gene_fraction_max csp_gene_fraction_max
		csp_locus_score_max breakpoint_gene_flank abundance_min_loci
		abundance_min_fold abundance_max_fold abundance_max_modified_z
	)} = (230, 8, 1200, 400, 1, 1, 3, 0.25, 0.05, 0.1, 50, 8, 1 / 3, 3, 3.5);
	@row{keys %override} = values %override;
	return \%row;
}

my $first = sample_row(
	sample => 'sample.1', worker => 0, assembly_driver => 'driver.1',
	status => 'processed', candidate_mgs => 30, used_mgs => 26,
	skipped_mgs => 4, candidate_loci => 22021, consensus_proteins => 24368,
	retained_loci => 9608, pre_abundance_loci => 21734,
	post_abundance_loci => 21684, missing_consensus_loci => 10,
	low_depth_loci => 1, ambiguous_loci => 285,
	abundance_filtered_loci => 50, skip_too_few_valid_sequences => 4,
);
my $second = sample_row(
	sample => 'sample.2', worker => 1, assembly_driver => 'driver.2',
	status => 'unavailable', candidate_mgs => 10, skipped_mgs => 10,
);

my $summary = aggregate_sample_rows([$first, $second], 'ALL');
is($summary->{samples}, 2, 'aggregate counts samples');
is($summary->{workers}, 2, 'aggregate counts contributing workers');
is($summary->{assembly_drivers}, 2, 'aggregate counts assembly drivers');
is($summary->{processed_samples}, 1, 'aggregate counts processed samples');
is($summary->{status_counts}, 'processed=1,unavailable=1',
	'aggregate reports every terminal sample status');
is($summary->{candidate_mgs}, 40, 'candidate MGS are summed');
is($summary->{used_mgs}, 26, 'used MGS are summed');
is($summary->{skipped_mgs}, 14, 'skipped MGS are summed');
is($summary->{used_fraction}, '0.650000',
	'overall used fraction is weighted by candidate MGS');
is($summary->{retained_loci}, 9608, 'retained loci are summed');
is($summary->{mean_loci_per_used_mgs}, '369.538',
	'mean retained loci is calculated across all used MGS');
is($summary->{max_genes}, 400, 'selected control parameters are retained');

my @summary_columns = sample_summary_columns();
ok(grep($_ eq 'status_counts', @summary_columns),
	'summary schema exposes terminal status counts');
ok(grep($_ eq 'abundance_filtered_loci', @summary_columns),
	'summary schema exposes abundance-filter losses');

my $mixed = sample_row(
	sample => 'sample.3', worker => 1, assembly_driver => 'driver.3',
	status => 'processed', max_genes => 500,
);
eval { aggregate_sample_rows([$first, $mixed], 'mixed') };
like($@, qr/Inconsistent max_genes values/,
	'aggregation rejects worker tables produced with different parameters');

my $malformed = sample_row(
	sample => 'sample.4', worker => 1, assembly_driver => 'driver.4',
	status => 'processed', candidate_mgs => 'not-a-number',
);
eval { aggregate_sample_rows([$first, $malformed], 'malformed') };
like($@, qr/Non-numeric sample statistic/,
	'aggregation rejects malformed numeric statistics');

done_testing();
