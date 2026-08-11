use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

my $root = File::Spec->catdir($Bin, q{..});
my $script = File::Spec->catfile($root, q{secScripts}, q{phylo}, q{buildTree5.pl});

open my $script_handle, '<', $script or die "Cannot read $script: $!";
my $script_text = do { local $/; <$script_handle> };
close $script_handle or die "Cannot close $script: $!";

my ($epa_resource_helper) = $script_text =~
	/(sub epaResourcePlan \{.*?return \(\$threads, \$memoryMB\);\n\})/s;
my ($iqtree_explicit_helper, $iqtree_model_helper) = $script_text =~
	/(sub iqtreeExplicitEpaModel \{.*?\n\})\n\n(sub iqtreePlacementModel \{.*?\n\})\n\nsub epaModelArtifact/s;
BAIL_OUT('Cannot extract EPA-ng helper functions')
	unless defined($epa_resource_helper) && defined($iqtree_explicit_helper)
		&& defined($iqtree_model_helper);
my $epa_helpers = "$epa_resource_helper\n$iqtree_explicit_helper\n$iqtree_model_helper";
my $helpers_loaded = eval "package TestBuildTreeEpaHelpers; $epa_helpers; 1;";
ok($helpers_loaded, 'EPA-ng model and resource helpers load independently')
	or diag($@);
my ($classification_helper) = $script_text =~
	/(sub readStrictBackboneClassification \{.*?\n\})\n\nsub runEpaOnlyPlacement/s;
BAIL_OUT('Cannot extract strict-backbone classification reader')
	unless defined $classification_helper;
my $classification_helper_loaded = eval
	"package TestBuildTreeEpaClassification; $classification_helper; 1;";
ok($classification_helper_loaded,
	'EPA-only strict-backbone classification reader loads independently')
	or diag($@);

my $temporary = tempdir(CLEANUP => 1);
my $iqtree_prefix = File::Spec->catfile($temporary, 'IQtree_allsites.backbone');
sub write_test_file {
	my ($path, $contents) = @_;
	open my $handle, '>', $path or die "Cannot write $path: $!";
	print {$handle} $contents;
	close $handle or die "Cannot close $path: $!";
}
write_test_file("$iqtree_prefix.iqtree", <<'IQTREE');
Best-fit model according to BIC: GTR+F+G2
Rate parameter R:
  A-C: 0.9
  A-G: 3.1
  A-T: 1.2
  C-G: 0.8
  C-T: 4.2
  G-T: 1.0
State frequencies: (empirical counts from alignment)
  pi(A) = 0.29
  pi(C) = 0.21
  pi(G) = 0.23
  pi(T) = 0.27
Gamma shape alpha: 0.73
IQTREE
is(TestBuildTreeEpaHelpers::iqtreePlacementModel($iqtree_prefix),
	'GTR{0.9/3.1/1.2/0.8/4.2/1}+FU{0.29/0.21/0.23/0.27}+G2{0.73}',
	'IQ-TREE fitted GTR rates, base frequencies, and gamma shape are passed explicitly');

my $invariantReport = <<'IQTREE';
Rate parameter R:
A-C: 1
A-G: 2
A-T: 3
C-G: 4
C-T: 5
G-T: 1
pi(A) = 0.1
pi(C) = 0.2
pi(G) = 0.3
pi(T) = 0.4
Proportion of invariable sites: 0.15
Gamma shape alpha: 0.6
IQTREE
is(TestBuildTreeEpaHelpers::iqtreeExplicitEpaModel(
	'GTR+F+I+G4', $invariantReport),
	'GTR{1/2/3/4/5/1}+FU{0.1/0.2/0.3/0.4}+I{0.15}+G4{0.6}',
	'invariant-site proportion and gamma categories are retained in the explicit descriptor');

write_test_file("$iqtree_prefix.iqtree", "IQ-TREE report without a model label\n");
write_test_file("$iqtree_prefix.log",
	"Command: iqtree3 -s alignment.fna -m HKY+F+G4 -T 12\n");
is(TestBuildTreeEpaHelpers::iqtreePlacementModel($iqtree_prefix), 'HKY+F+G4',
	'IQ-TREE command-line model is the final parser fallback');
is_deeply(
	[TestBuildTreeEpaHelpers::epaResourcePlan(12, 12, -1, 4500)],
	[2, 2700],
	'EPA threads are capped to one per GB of the derived planning budget');
is_deeply(
	[TestBuildTreeEpaHelpers::epaResourcePlan(4, 12, 0, 4500)],
	[4, 0],
	'zero planning memory keeps the requested core-capped thread count');

my $classification = File::Spec->catfile($temporary, 'strict_backbone.samples.tsv');
write_test_file($classification, join("\n",
	join("\t", qw(sample tree_role reason informative_positions q90_informative
		backbone_overlap_nt backbone_overlap_loci backbone_state_divergence)),
	join("\t", qw(backbone1 backbone validated 1000 900 1000 8 0)),
	join("\t", qw(query1 placement sparse 450 900 425 4 0.03)),
	join("\t", qw(outlier1 excluded divergent 200 900 100 1 0.8)),
)."\n");
my $classification_state =
	TestBuildTreeEpaClassification::readStrictBackboneClassification($classification);
is_deeply($classification_state->{placement}, ['query1'],
	'EPA-only recovery reuses exactly the samples classified for placement');
is($classification_state->{backbone_overlap}{query1}{backbone_overlap_nt}, 425,
	'EPA-only placement reporting retains the original backbone-overlap metric');

like($script_text,
	qr/"epaOnly=i" => \\\$epaOnly.*?if \(\$epaOnly\).*?runEpaOnlyPlacement\(.*?exit\(0\)/s,
	'EPA-only mode exits through its dedicated placement path before ordinary MSA and inference work');
like($script_text,
	qr/sub runEpaOnlyPlacement.*?requires a validated IQ-TREE backbone.*?map_epa_placements_to_backbone\(.*?write_epa_placed_tree\(\$backboneTreeText, \$primaryTree.*?backbone retained=\$backboneTree/s,
	'EPA-only mode maps jplace edges and grafts placements onto the retained backbone');

my $coordinate_bounds_checks = () = $script_text =~ /next if \$position >= length\(\$sequence\);/g;
cmp_ok($coordinate_bounds_checks, '>=', 2,
	'taxon-aware raw and alignment coordinate scorers safely skip uneven sequence tails');
like($script_text,
	qr/if \(defined\(\$candidateSelection->\{terminal_reason\}\).*?writeOutcomeMarker\(\$terminalMarker, 'valid_no_tree', \$reason.*?exit\(0\)/s,
	'the candidate-selection terminal result is published as a persistent valid no-tree outcome');
like($script_text,
	qr/unless \(keys %metrics\).*?terminal_reason => 'taxon_aware_no_category_with_three_usable_samples'/s,
	'a taxon-aware candidate set with fewer than three usable samples returns a stable terminal reason');
my $placement_outlier_calls = () = $script_text =~ /filter_epa_placement_outliers\(/g;
cmp_ok($placement_outlier_calls, '>=', 2,
	'fresh and EPA-only placement publication both apply pendant-branch outlier QC');
my $backbone_mapping_calls =
	() = $script_text =~ /map_epa_placements_to_backbone\(/g;
cmp_ok($backbone_mapping_calls, '>=', 2,
	'normal and EPA-only publication both map jplace edges onto the backbone');
like($script_text, qr/strict_backbone\.epa_backbone_grafts\.tsv/,
	'EPA backbone grafting publishes a per-edge mapping report');
unlike($script_text, qr/write_epa_placed_tree\(\$epaResult->\{tree\}/,
	'the jplace Newick tree is never used as the publication template');
like($script_text,
	qr/pendant_outlier_limit placement_filter_reason reason/s,
	'EPA placement reports expose the applied cutoff and exclusion reason');
is(system($^X, q{-I}.$root, q{-c}, $script), 0, q{buildTree5.pl compiles});

done_testing;
