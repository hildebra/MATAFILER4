use strict;
use warnings;

use File::Spec;
use FindBin qw($Bin);
use Test::More;

sub slurp {
	my ($path) = @_;
	open my $fh, '<', $path or die "Cannot read $path: $!";
	local $/;
	return <$fh>;
}

my $root = File::Spec->catdir($Bin, '..');
my $strain = slurp(File::Spec->catfile($root, 'secScripts', 'MGS', 'strain_within.pl'));
my $mgs = slurp(File::Spec->catfile($root, 'secScripts', 'MGS.pl'));
my $tree = slurp(File::Spec->catfile($root, 'secScripts', 'phylo', 'buildTree5.pl'));
my $mosaic = slurp(File::Spec->catfile(
	$root, 'secScripts', 'MGS', 'prepare_mosaic_loci.pl',
));

like($strain, qr/my \$noFilter = \$disableQC \? 1 : 0/,
	'QC disabling is controlled only by the explicit disableQC option');
like($strain, qr/\$takeAll = \$noGeneLimit/,
	'the no-gene-limit option controls only count capping');
like($strain, qr/next if !\$noGeneLimit && \$locCnt >= \$maxNGenes/,
	'unlimited mode does not accidentally reject every validated locus');
unlike($strain, qr/\$mode\s*=\s*["']MGSall["']/,
	'unlimited mode no longer enters the legacy QC-free MGSall mode');
like($strain, qr/my %FILTER_DEFAULT = \(.*?sub usage \{.*?maximum_genes_per_sample/s,
	'strain runtime defaults are reused by its help text');
like($strain,
	qr/prepare_mosaic_loci => 1.*?"prepareMosaicLoci=i".*?default \$default->\{prepare_mosaic_loci\}/s,
	'automatic mosaic preparation is enabled and documented from a shared default');

like($mgs,
	qr/MGS_mosaic_scr.*?prepare_mosaic|MGS_mosaic_scr/s,
	'MGS resolves the catalogue-wide mosaic preprocessing step');
like($mgs,
	qr/\$mosaicScr .*?-output \$mosaicCatalogue.*?\$strain1scr .*?-mosaicLoci \$mosaicCatalogue/s,
	'MGS confirms mosaics before passing their catalogue to strain filtering');
like($mosaic,
	qr/discover_mosaic_candidates.*?minimap2.*?confirm_mosaic_candidates.*?select_outgroup_panel/s,
	'mosaic preprocessing creates candidates, aligns catalogue-wide, confirms pairs, and consolidates outgroups');

like($tree, qr/use Mods::StrainPlacement.*?split_strict_backbone/s,
	'buildTree5 uses the tested strict-backbone splitter');
like($tree,
	qr/treeAtHeart\(\$tOhr\).*?nearest_backbone_placements.*?write_placed_tree/s,
	'buildTree5 places deferred samples only after backbone inference');
like($tree, qr/strict_backbone\.placements\.tsv/,
	'buildTree5 emits auditable overlap and distance placements');

done_testing();
