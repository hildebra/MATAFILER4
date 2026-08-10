use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::StrainPlacement qw(
	read_sample_qc split_strict_backbone
	read_epa_jplace filter_epa_placement_outliers write_epa_placed_tree
);

sub write_file {
	my ($path, $text) = @_;
	open my $fh, '>', $path or die $!;
	print {$fh} $text;
	close $fh or die $!;
}

sub slurp {
	my ($path) = @_;
	open my $fh, '<', $path or die $!;
	local $/;
	return <$fh>;
}

my $tmp = tempdir(CLEANUP => 1);
my $full = File::Spec->catfile($tmp, 'full.fna');
my $backbone = File::Spec->catfile($tmp, 'backbone.fna');
my $queries = File::Spec->catfile($tmp, 'placement.fna');
write_file($full, join('',
	">A\nAACCGGTTAACC\n",
	">B\nAACCGGTTAACA\n",
	">C\nAACCGGTTAAGG\n",
	">D\nAACCGG------\n",
	">E\nAACCGGTTAACC\n",
	">F\nAAC---------\n",
	">G\n------------\n",
));
my $split = split_strict_backbone(
	$full, $backbone, $queries, {E => 'placement'},
	{coverage_fraction => 0.35, minimum_backbone => 3},
);
is_deeply($split->{backbone}, [qw(A B C D E)],
	'strict backbone retains all but severe coverage outliers');
is_deeply($split->{placement}, [qw(F)],
	'only the severe low-coverage sample is deferred to placement');
is_deeply($split->{excluded}, [qw(G)],
	'all-gap samples are excluded before EPA-ng placement');
is($split->{reason}{G}, 'no_informative_alignment_sites',
	'all-gap placement exclusion records an actionable audit reason');
is($split->{reason}{E}, 'retained_after_locus_qc_masking',
	'locus-QC status alone does not remove a well-covered sample after masking');

my $eligible_backbone = File::Spec->catfile($tmp, 'eligible-backbone.fna');
my $eligible_queries = File::Spec->catfile($tmp, 'eligible-placement.fna');
my $eligible_split = split_strict_backbone(
	$full, $eligible_backbone, $eligible_queries, {E => 'placement'},
	{
		coverage_fraction => 0.35,
		minimum_backbone => 3,
		placement_eligible => {F => 0},
		placement_ineligible_reason => {F => 'below_placement_gene_fraction'},
	},
);
is_deeply($eligible_split->{backbone}, [qw(A B C D E)],
	'coverage-adequate samples still define the broad initial backbone');
is_deeply($eligible_split->{placement}, [],
	'an ineligible sparse sample is not appended as a placement');
is_deeply($eligible_split->{excluded}, [qw(F G)],
	'a low-coverage sample failing the restored gene threshold is explicitly excluded');
like($eligible_split->{reason}{F}, qr/below_placement_gene_fraction/,
	'placement exclusion preserves the threshold reason for audit');

my $policy_backbone = File::Spec->catfile($tmp, 'policy-backbone.fna');
my $policy_queries = File::Spec->catfile($tmp, 'policy-placement.fna');
my $policy_split = split_strict_backbone(
	$full, $policy_backbone, $policy_queries, {},
	{
		coverage_fraction => 0.35,
		minimum_backbone => 3,
		backbone_eligible => {D => 0},
		backbone_ineligible_reason => {D => 'below_backbone_gene_fraction'},
		placement_eligible => {D => 1},
	},
);
is_deeply($policy_split->{placement}, [qw(D F)],
	'backbone-specific coverage rejection defers an otherwise placeable sample');
like($policy_split->{reason}{D}, qr/below_backbone_gene_fraction/,
	'backbone deferral reason remains available for the final classification audit');

my $overlap_full = File::Spec->catfile($tmp, 'overlap-full.fna');
my $overlap_backbone = File::Spec->catfile($tmp, 'overlap-backbone.fna');
my $overlap_queries = File::Spec->catfile($tmp, 'overlap-placement.fna');
my $overlap_partition = File::Spec->catfile($tmp, 'overlap.partition');
write_file($overlap_full, join('',
	">A\nAACCGGTTAACC\n",
	">B\nAACCGGTTAACA\n",
	">C\nAACCGGTTAAGG\n",
	">Q\nAA----TG----\n",
	">R\nAA----------\n",
	">S\nAAAA--------\n",
));
write_file($overlap_partition, join('',
	"DNA, locus_1 = 1-6\n",
	"DNA, locus_2 = 7-12\n",
));
my $overlap_split = split_strict_backbone(
	$overlap_full, $overlap_backbone, $overlap_queries, {},
	{
		coverage_fraction => 0.5,
		minimum_backbone => 3,
		partition_file => $overlap_partition,
		minimum_backbone_overlap_nt => 4,
		minimum_backbone_overlap_loci => 2,
	},
);
is_deeply($overlap_split->{placement}, [qw(Q)],
	'a sparse query is retained only when it overlaps the backbone in two loci');
is_deeply($overlap_split->{excluded}, [qw(R S)],
	'queries failing either actual backbone-overlap gate are excluded');
like($overlap_split->{reason}{R}, qr/below_backbone_overlap_nt/,
	'backbone-overlap NT failure is explicit in the audit');
like($overlap_split->{reason}{S}, qr/below_backbone_overlap_loci/,
	'backbone-overlap locus failure is explicit in the audit');
is($overlap_split->{backbone_overlap}{Q}{backbone_overlap_nt}, 4,
	'backbone overlap counts only jointly supported coordinates');
is($overlap_split->{backbone_overlap}{Q}{backbone_overlap_loci}, 2,
	'backbone overlap reports the number of intersecting loci');
cmp_ok(abs($overlap_split->{backbone_overlap}{Q}{backbone_state_divergence} - 0.25),
	'<', 1e-12, 'backbone-state divergence is recorded without filtering the query');

my $fallback_backbone = File::Spec->catfile($tmp, 'fallback-backbone.fna');
my $fallback_queries = File::Spec->catfile($tmp, 'fallback-placement.fna');
my $fallback = split_strict_backbone(
	$full, $fallback_backbone, $fallback_queries,
	{E => 'placement', F => 'placement'},
	{coverage_fraction => 0.95, minimum_backbone => 5},
);
ok($fallback->{fallback}, 'underpowered strict backbones fall back explicitly');
ok(exists($fallback->{requested_reason}{F}),
	'fallback retains the QC reason for audit rather than silently relabelling the sample');

my $jplace = File::Spec->catfile($tmp, 'epa_result.jplace');
my $placed_tree = File::Spec->catfile($tmp, 'placed.treefile');
write_file($jplace, <<'JPLACE');
{"tree":"(A:0.1{1},(B:0.1{2},C:0.1{3}):0.2{4});","placements":[
 {"p":[[1,-10.0,0.90,0.04,0.01],[2,-12.0,0.10,0.03,0.02]],"n":["D"]},
 {"p":[[1,-11.0,0.80,0.06,0.02]],"n":["E"]}
]}
JPLACE
my $epa = read_epa_jplace($jplace, [qw(D E F)]);
is($epa->{placements}{D}{edge}, 1,
	'EPA-ng parser selects the highest-likelihood-weight edge');
is($epa->{placements}{D}{likelihood_weight_ratio}, 0.90,
	'EPA-ng parser retains likelihood-weight support');
is($epa->{placements}{D}{candidate_placements}, 2,
	'EPA-ng audit retains the number of candidate placements');
cmp_ok(abs($epa->{placements}{D}{edpl} - 0.1188), '<', 1e-12,
	'EDPL records spatial placement uncertainty on the reference tree');
is($epa->{placements}{E}{edpl}, 0,
	'a query with one candidate placement has zero EDPL');
is($epa->{placements}{F}{status}, 'not_reported',
	'queries absent from jplace remain explicit in the audit report');
write_epa_placed_tree($epa->{tree}, $placed_tree, $epa->{placements});
my $placed_text = slurp($placed_tree);
like($placed_text, qr/D:0\.01/, 'display tree contains the EPA-ng pendant branch');
like($placed_text, qr/E:0\.02/, 'multiple EPA-ng placements on one edge are retained');
like($placed_text, qr/A:0\.04/, 'EPA-ng distal branch length is applied from the edge child');
like($placed_text, qr/B:0\.1/, 'unaffected backbone topology is retained');
unlike($placed_text, qr/\{\d+\}/, 'EPA-ng edge labels are not leaked into the published tree');

my $outlier_tree = '(A:0.001{1},B:0.002{2},C:0.003{3},MGS.out:0.5{4});';
my %outlier_placements = (
	near => {status => 'placed', edge => 1, distal_length => 0.0005,
		pendant_length => 0.01},
	far => {status => 'placed', edge => 2, distal_length => 0.001,
		pendant_length => 0.05},
);
my $outlier_qc = filter_epa_placement_outliers(
	$outlier_tree, \%outlier_placements, {outgroup => 'MGS.out'});
cmp_ok(abs($outlier_qc->{threshold} - 0.02), '<', 1e-12,
	'the conservative absolute floor controls a very compact backbone cutoff');
is_deeply($outlier_qc->{retained}, ['near'],
	'a placement within the adaptive pendant cutoff is retained');
is_deeply($outlier_qc->{excluded}, ['far'],
	'a clearly separated pendant-branch outlier is excluded');
is($outlier_placements{far}{status}, 'excluded_outlier',
	'the excluded placement cannot be published into the final tree');
is($outlier_placements{far}{placement_filter_reason}, 'pendant_length_outlier',
	'the placement report receives a stable outlier reason');
my $filtered_tree = File::Spec->catfile($tmp, 'outlier-filtered.treefile');
write_epa_placed_tree($outlier_tree, $filtered_tree, \%outlier_placements);
my $filtered_text = slurp($filtered_tree);
like($filtered_text, qr/near:0\.01/, 'retained placement remains in the final tree');
unlike($filtered_text, qr/far:/, 'pendant outlier is absent from the final tree');

my $qc = File::Spec->catfile($tmp, 'sampleQC.tsv');
write_file($qc, join('',
	"MGS\tsample\tstatus\tambiguous_fraction\tcsp_fraction\tvalidated_loci\n",
	"MGS1\tA\tbackbone\t0\t0\t20\n",
	"MGS1\tA\tplacement\t0.4\t0\t12\n",
));
is(read_sample_qc($qc)->{A}, 'placement',
	'placement status dominates duplicate QC rows');

done_testing();
