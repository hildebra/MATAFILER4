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

like($strain,
	qr/getProgPaths\("MGS_mosaic_scr"\).*?qsubSystem\(/s,
	'strain workflow resolves and submits catalogue-wide mosaic preprocessing');
like($strain,
	qr/qsubSystem\(.*?"MosaicMGS".*?qsubSystemJobAlive\(\[\$mosaicDependency\].*?unless -s \$mosaicLociFile/s,
	'strain workflow waits for and validates Mosaic before loading its catalogue');
like($mosaic,
	qr/discover_mosaic_candidates.*?my \$paf_path.*?system\(\@rtk_command\).*?read_rtk_mosaic_results.*?select_outgroup_panel/s,
	'mosaic preprocessing creates diagnostics, aligns catalogue-wide, lets rtk confirm pairs, and consolidates outgroups');
like($mosaic, qr/minimap_preset => 'asm20'.*?write_summary/s,
	'mosaic preprocessing uses an outgroup-sensitive alignment preset and records stage diagnostics');
like($mosaic,
	qr/select_interesting_records.*?Aligning .*?interesting genes before rtk2.*?'-c', '-D'.*?'-N', \$DEFAULT\{max_secondary_hits\}.*?Running rtk2 mosaic/s,
	'mosaic preprocessing bulk-aligns only informative genes before abundance-aware rtk confirmation');
like($mosaic,
	qr/open my \$minimap_fh, '-\|', \@command.*?rename \$temporary_paf, \$paf_path.*?system\(\@rtk_command\).*?read_paf_hits\(\$paf_path/s,
	'mosaic preprocessing materializes minimap PAF, runs rtk, then reuses that PAF for outgroups');
like($mosaic,
	qr/Mosaic preprocessing summary.*?Unique MGS-outgroup links:.*?Proposed outgroup gene links:.*?write_outgroup_table/s,
	'mosaic preprocessing prints useful final statistics and writes explicit outgroup proposals');
like($strain,
	qr/my \$mosaicDirectory = dirname\(\$mosaicLociFile\).*?prepare_mosaic_loci\.log.*?prepare_mosaic_loci\.sh/s,
	'strain stores Mosaic catalogues, logs, and submission scripts together');
like($strain,
	qr/my \$mosaicThreads = \$maxCores > 0 \? \$maxCores : \$numCores.*?qsubSystem\(.*?\$mosaicThreads.*?"\$\{mosaicMemGb\}G".*?"MosaicMGS"/s,
	'the submitted Mosaic job receives maximum strain cores and dedicated memory');
like($mgs,
	qr/"prepareMosaicLoci=i" => \\\$prepareMosaicLoci.*?\$ph2Cmd \.= "-mosaicLoci \$mosaicCatalogue " if \$prepareMosaicLoci/s,
	'MGS can omit mosaic preprocessing and the mosaic catalogue from strain analysis');
like($strain,
	qr/Mosaic checks disabled; same-COG catalogue clusters will remain separate and tree-based outgroups remain available/,
	'strain analysis reports its safe no-mosaic behavior explicitly');

like($tree, qr/use Mods::StrainPlacement.*?split_strict_backbone/s,
	'buildTree5 uses the tested strict-backbone splitter');
like($tree,
	qr/treeAtHeart\(\$tOhr\).*?nearest_backbone_placements.*?write_placed_tree/s,
	'buildTree5 places deferred samples only after backbone inference');
like($tree, qr/strict_backbone\.placements\.tsv/,
	'buildTree5 emits auditable overlap and distance placements');

done_testing();
