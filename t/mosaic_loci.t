use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::MosaicLoci qw(
	pair_key discover_mosaic_candidates read_paf_hits confirm_mosaic_candidates
	select_outgroup_panel read_mosaic_catalogue
);

sub hit {
	my ($query, $target, $identity, $matches) = @_;
	return {
		query => $query, target => $target, identity => $identity,
		matches => $matches, query_coverage => 0.95,
		target_coverage => 0.94, alignment_length => 1000,
	};
}

my @records = (
	{mgs => 'MGS1', cog => 'COG1', gene => '1', rank => 0},
	{mgs => 'MGS1', cog => 'COG1', gene => '2', rank => 1},
	{mgs => 'MGS1', cog => 'COG1', gene => '3', rank => 2},
);
my %members = (
	1 => 'S1__ctg_1,S2__ctg_1',
	2 => 'S3__ctg_1,S4__ctg_1',
	3 => 'S2__ctg_2,S5__ctg_1',
);
my %sequence = (1 => 'A' x 1000, 2 => 'A' x 950, 3 => 'A' x 980);
my $candidate = discover_mosaic_candidates(
	\@records, \%members, \%sequence, {minimum_length_ratio => 0.8},
);
is_deeply(
	[map { [$_->{left}, $_->{right}] } @{$candidate}],
	[['1', '2'], ['2', '3']],
	'candidate discovery retains non-cooccurring, length-compatible same-COG seeds',
);
my %rare_overlap_members = (
	1 => join(',', map { "S$_"."__ctg_1" } 1 .. 10),
	2 => join(',', 'S1__ctg_2', map { "T$_"."__ctg_1" } 2 .. 10),
);
my $rare_overlap = discover_mosaic_candidates(
	[@records[0, 1]], \%rare_overlap_members, \%sequence,
	{
		minimum_length_ratio => 0.8,
		maximum_overlap_samples => 1,
		maximum_sample_overlap_fraction => 0.15,
	},
);
is(scalar(@{$rare_overlap}), 1,
	'one rare co-occurrence remains a candidate so mixed-strain samples do not force locus loss');

my %unique_hits = (
	1 => [hit(1, 2, 0.94, 940), hit(1, 2, 0.93, 900), hit(1, 9, 0.91, 850)],
	2 => [hit(2, 1, 0.94, 930), hit(2, 8, 0.90, 820)],
);
my ($confirmed, $rejected) = confirm_mosaic_candidates(
	[$candidate->[0]], \%unique_hits, {minimum_score_margin => 0.02},
);
is(scalar(@{$confirmed}), 1,
	'reciprocal partners with a clear whole-catalogue lead are confirmed');
is(scalar(@{$rejected}), 0, 'the reciprocal unique pair is not rejected');

my %nonunique_hits = %unique_hits;
$nonunique_hits{1} = [hit(1, 2, 0.94, 940), hit(1, 9, 0.94, 935)];
($confirmed, $rejected) = confirm_mosaic_candidates(
	[$candidate->[0]], \%nonunique_hits, {minimum_score_margin => 0.02},
);
is(scalar(@{$confirmed}), 0,
	'a near-tied third catalogue homologue prevents mosaic merging');
is($rejected->[0]{reason}, 'nonunique_catalogue_hit',
	'rejection records why orthology is not uniquely interpretable');

my @panel_records = (
	{mgs => 'MGS1', cog => 'C1', gene => 'a1'},
	{mgs => 'MGS1', cog => 'C2', gene => 'a2'},
	{mgs => 'MGS1', cog => 'C3', gene => 'a3'},
	{mgs => 'MGS2', cog => 'C1', gene => 'b1'},
	{mgs => 'MGS2', cog => 'C2', gene => 'b2'},
	{mgs => 'MGS2', cog => 'C3', gene => 'b3'},
	{mgs => 'MGS3', cog => 'C1', gene => 'c1'},
	{mgs => 'MGS3', cog => 'C2', gene => 'c2'},
);
my %panel_hits = (
	a1 => [hit('a1', 'b1', 0.88, 880), hit('a1', 'c1', 0.90, 900)],
	a2 => [hit('a2', 'b2', 0.87, 870), hit('a2', 'c2', 0.90, 900)],
	a3 => [hit('a3', 'b3', 0.89, 890)],
);
my ($outgroup, $gene_map) = select_outgroup_panel(
	\@panel_records, \%panel_hits,
	{minimum_loci => 2, target_identity => 0.88},
);
is($outgroup->{MGS1}{target_mgs}, 'MGS2',
	'outgroup consolidation favours broad locus representation before closeness');
is(scalar(keys %{$gene_map->{MGS1}}), 3,
	'selected outgroup retains the best per-locus homologues');

my $tmp = tempdir(CLEANUP => 1);
my $paf = File::Spec->catfile($tmp, 'catalogue.paf');
open my $paf_fh, '>', $paf or die $!;
print {$paf_fh} join("\t",
	qw(1 1000 0 950 + 2 1000 0 940 900 960 60)), "\n";
close $paf_fh;
my $parsed_hits = read_paf_hits($paf);
cmp_ok(abs($parsed_hits->{1}[0]{identity} - 0.9375), '<', 1e-8,
	'PAF identity is calculated from aligned matches and block length');
is($parsed_hits->{1}[0]{query_coverage}, 0.95,
	'PAF query coverage is retained for minimum-alignment checks');

my $catalogue = File::Spec->catfile($tmp, 'mosaics.tsv');
open my $fh, '>', $catalogue or die $!;
print {$fh} "MOSAIC\tMGS1\tCOG1\t1\t2\t0.94\t0.95\t0.94\n";
print {$fh} "OUTGROUP\tMGS1\tMGS2\t3\t0.88\n";
print {$fh} "OUTGROUP_GENE\tMGS1\tMGS2\ta1\tb1\t0.88\t0.95\n";
close $fh;
my ($pairs, $outgroups, $outgroup_genes) = read_mosaic_catalogue($catalogue);
ok($pairs->{pair_key(1, 2)}, 'confirmed pair catalogue round-trips into an allowlist');
is($outgroups->{MGS1}, 'MGS2', 'preferred outgroup is parsed');
is($outgroup_genes->{MGS1}{a1}, 'b1', 'preferred per-locus outgroup gene is parsed');

done_testing();
