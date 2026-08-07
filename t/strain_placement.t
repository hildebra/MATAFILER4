use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::StrainPlacement qw(
	read_sample_qc split_strict_backbone
	read_epa_jplace write_epa_placed_tree
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
));
my $split = split_strict_backbone(
	$full, $backbone, $queries, {E => 'placement'},
	{coverage_fraction => 0.35, minimum_backbone => 3},
);
is_deeply($split->{backbone}, [qw(A B C D E)],
	'strict backbone retains all but severe coverage outliers');
is_deeply($split->{placement}, [qw(F)],
	'only the severe low-coverage sample is deferred to placement');
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
is_deeply($eligible_split->{excluded}, [qw(F)],
	'a low-coverage sample failing the restored gene threshold is explicitly excluded');
like($eligible_split->{reason}{F}, qr/below_placement_gene_fraction/,
	'placement exclusion preserves the threshold reason for audit');

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
is($epa->{placements}{F}{status}, 'not_reported',
	'queries absent from jplace remain explicit in the audit report');
write_epa_placed_tree($epa->{tree}, $placed_tree, $epa->{placements});
my $placed_text = slurp($placed_tree);
like($placed_text, qr/D:0\.01/, 'display tree contains the EPA-ng pendant branch');
like($placed_text, qr/E:0\.02/, 'multiple EPA-ng placements on one edge are retained');
like($placed_text, qr/A:0\.04/, 'EPA-ng distal branch length is applied from the edge child');
like($placed_text, qr/B:0\.1/, 'unaffected backbone topology is retained');
unlike($placed_text, qr/\{\d+\}/, 'EPA-ng edge labels are not leaked into the published tree');

my $qc = File::Spec->catfile($tmp, 'sampleQC.tsv');
write_file($qc, join('',
	"MGS\tsample\tstatus\tambiguous_fraction\tcsp_fraction\tvalidated_loci\n",
	"MGS1\tA\tbackbone\t0\t0\t20\n",
	"MGS1\tA\tplacement\t0.4\t0\t12\n",
));
is(read_sample_qc($qc)->{A}, 'placement',
	'placement status dominates duplicate QC rows');

done_testing();
