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

my ($epa_helpers) = $script_text =~
	/(sub epaMemoryPrefix \{.*?\n\}\n\nsub iqtreePlacementModel \{.*?\n\})\n\nsub epaModelArtifact/s;
BAIL_OUT('Cannot extract EPA-ng helper functions') unless defined $epa_helpers;
my $helpers_loaded = eval "package TestBuildTreeEpaHelpers; $epa_helpers; 1;";
ok($helpers_loaded, 'EPA-ng model and memory helpers load independently')
	or diag($@);

my $temporary = tempdir(CLEANUP => 1);
my $iqtree_prefix = File::Spec->catfile($temporary, 'IQtree_allsites.backbone');
sub write_test_file {
	my ($path, $contents) = @_;
	open my $handle, '>', $path or die "Cannot write $path: $!";
	print {$handle} $contents;
	close $handle or die "Cannot close $path: $!";
}
write_test_file("$iqtree_prefix.iqtree",
	"Best-fit model according to BIC: GTR+F+G2\n");
is(TestBuildTreeEpaHelpers::iqtreePlacementModel($iqtree_prefix), 'GTR+F+G2',
	'IQ-TREE model is parsed without the legacy Model of substitution label');

write_test_file("$iqtree_prefix.iqtree", "IQ-TREE report without a model label\n");
write_test_file("$iqtree_prefix.log",
	"Command: iqtree3 -s alignment.fna -m HKY+F+G4 -T 12\n");
is(TestBuildTreeEpaHelpers::iqtreePlacementModel($iqtree_prefix), 'HKY+F+G4',
	'IQ-TREE command-line model is the final parser fallback');
is(TestBuildTreeEpaHelpers::epaMemoryPrefix(2048),
	"ulimit -S -v 2097152;\n",
	'EPA memory MB are converted to a direct child-shell virtual-memory limit');

my $coordinate_bounds_checks = () = $script_text =~ /next if \$position >= length\(\$sequence\);/g;
cmp_ok($coordinate_bounds_checks, '>=', 2,
	'taxon-aware raw and alignment coordinate scorers safely skip uneven sequence tails');
is(system($^X, q{-I}.$root, q{-c}, $script), 0, q{buildTree5.pl compiles});

done_testing;
