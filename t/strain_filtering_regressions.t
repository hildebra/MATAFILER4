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
	qr/my \$records = read_mgs_records\(\$mgs_file.*?my \$query_fasta.*?'-p', '0'.*?\$query_fasta, \$query_fasta.*?\$rtk, 'mosaic'.*?'-reference', \$mgs_file/s,
	'mosaic preprocessing self-aligns the complete raw MGS gene set before rtk confirmation');
unlike($mosaic, qr/'-groups'/,
	'rtk Mosaic discovery is not restricted to genes sharing a NOG');
like($mosaic,
	qr/my \$core_records.*?my \@outgroup_records.*?select_outgroup_panel\(\s*\\\@outgroup_records/s,
	'the core MGS table is reserved for stable outgroup selection');
like($mosaic, qr/minimap_preset => 'asm20'.*?write_summary/s,
	'mosaic preprocessing uses an outgroup-sensitive alignment preset and records stage diagnostics');
like($mosaic,
	qr/open my \$minimap_fh, '-\|', \@command.*?rename \$temporary_paf, \$paf_path.*?system\(\@rtk_command\).*?read_paf_hits\(\$paf_path/s,
	'mosaic preprocessing materializes minimap PAF, runs rtk, then reuses that PAF for outgroups');
like($mosaic,
	qr/my \$work_dir = tempdir.*?\$paf_path = "\$work_dir\/raw_mgs\.minimap2\.paf".*?my \$rtk_prefix = "\$work_dir\/rtk".*?copy_atomic\(\$rtk_report, \$candidate_output\).*?remove_legacy_intermediates\(\$output,/s,
	'generated FASTA, PAF, and native rtk outputs stay temporary while audit outputs are published');
unlike($mosaic, qr/write_(?:rejections|outgroup_table)\(/,
	'rejected and outgroup-only duplicate tables are no longer published');
like($strain,
	qr/my \$mosaicRunDirectory = tempdir\(.*?DIR => \$mosaicDirectory.*?\$mosaicRunDirectory, 'prepare_mosaic_loci\.log'.*?\$mosaicRunDirectory, 'prepare_mosaic_loci\.sh'.*?remove_tree\(\$mosaicRunDirectory\)/s,
	'successful scheduler scripts and logs are confined to and removed with a per-run workspace');
like($strain,
	qr/cleanupMosaicIntermediates\(\$mosaicLociFile\).*?sub cleanupMosaicIntermediates.*?\.minimap2\.paf.*?\.outgroups\.tsv/s,
	'strain reruns remove obsolete regenerable Mosaic artifacts once the catalogue exists');
like($strain,
	qr/strainRecovery\.tsv.*?mergeRecoveryLogs.*?writeRecoveryRow\(\$MGS, \$sd3, 'filtered'.*?writeRecoveryRow\(\$MGS, \$sd3, 'recovered'/s,
	'every evaluated sample-MGS is persisted as recovered or filtered with a reason');
like($strain,
	qr/sub writeStrainSummary.*?average_genes_per_recovered_MAG.*?recovered_MAGs\.genes_gt_.*?filtered_reason\./s,
	'the output-folder summary reports gene statistics, cumulative MAG thresholds, and filter reasons');
like($strain,
	qr/recovered_mosaic_loci.*?mosaic_outgroups_used.*?used_mosaic_outgroup/s,
	'the summary counts retained Mosaic loci and Mosaic-derived outgroup use, including resumed inputs');
like($strain,
	qr/my \$mosaicThreads = \$maxCores > 0 \? \$maxCores : \$numCores.*?qsubSystem\(.*?\$mosaicThreads.*?"\$\{mosaicMemGb\}G".*?"MosaicMGS"/s,
	'the submitted Mosaic job receives maximum strain cores and dedicated memory');
like($mgs,
	qr/"prepareMosaicLoci=i" => \\\$prepareMosaicLoci.*?-MGS \$finalClustersFilt -mosaicMGS \$finalClusters2.*?\$ph2Cmd \.= "-mosaicLoci \$mosaicCatalogue " if \$prepareMosaicLoci/s,
	'MGS passes raw clusters for Mosaic while retaining core clusters for strain analysis');
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
