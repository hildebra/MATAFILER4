use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::StrainPlacement qw(
	read_sample_qc split_strict_backbone
	nearest_backbone_placements write_placed_tree
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

my $placement_backbone = File::Spec->catfile($tmp, 'placement-backbone.fna');
my $placement_queries = File::Spec->catfile($tmp, 'placement-queries.fna');
write_file($placement_backbone, join('',
	">A\nAACCGGTTAACC\n",
	">B\nAACCGGTTAACA\n",
	">C\nAACCGGTTAAGG\n",
));
write_file($placement_queries, join('',
	">D\nAACCGG------\n",
	">E\nAACCGGTTAACC\n",
));
my $placements = nearest_backbone_placements(
	$placement_backbone, $placement_queries, 5,
);
is($placements->{D}{anchor}, 'A',
	'lower-coverage sample is placed using only its validated overlap');
is($placements->{D}{overlap}, 6, 'placement reports the validated overlap');
is($placements->{E}{anchor}, 'A', 'QC-marked complete sample has an exact anchor');

my $tree = File::Spec->catfile($tmp, 'backbone.treefile');
my $placed_tree = File::Spec->catfile($tmp, 'placed.treefile');
write_file($tree, "(A:0.1,(B:0.1,C:0.1):0.1);\n");
write_placed_tree($tree, $placed_tree, $placements);
my $placed_text = slurp($placed_tree);
like($placed_text, qr/D:/, 'display tree contains the low-coverage placement');
like($placed_text, qr/E:/, 'display tree contains the QC-deferred placement');
like($placed_text, qr/B:0\.1/, 'unaffected backbone topology is retained');

my $qc = File::Spec->catfile($tmp, 'sampleQC.tsv');
write_file($qc, join('',
	"MGS\tsample\tstatus\tambiguous_fraction\tcsp_fraction\tvalidated_loci\n",
	"MGS1\tA\tbackbone\t0\t0\t20\n",
	"MGS1\tA\tplacement\t0.4\t0\t12\n",
));
is(read_sample_qc($qc)->{A}, 'placement',
	'placement status dominates duplicate QC rows');

done_testing();
