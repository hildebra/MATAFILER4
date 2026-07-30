use strict;
use warnings;

use File::Spec;
use FindBin qw($Bin);
use Test::More;

my $root = File::Spec->catdir($Bin, '..');
my $script = File::Spec->catfile($root, 'secScripts', 'phylo', 'buildTree5.pl');

open my $fh, '<', $script or die "Cannot read $script: $!";
my $source = do { local $/; <$fh> };
close $fh;

my $compile_status = system($^X, '-I'.$root, '-c', $script);
is($compile_status, 0, 'buildTree5.pl compiles');

like($source, qr/my \$version = 5\.22;/,
	'Perl-owned tree job lifecycle increments the workflow version');
like($source,
	qr/"strainWithinPreset=i".*?if \(\$strainWithinPreset\) \{.*?\$useAA4tree = 0;.*?\$ntCntTotal = 400;.*?\$strictBackbone = 1;.*?\$continue = 1;.*?\$doDNDS = 0;.*?\$doTheta = 0;/s,
	"buildTree strain preset owns the fixed strain-tree settings");
like($source,
	qr/"stagedInputDir=s".*?"tmpSubdir=s".*?"completionMarker=s".*?publishStagedTreeInputs\(\$stagedInputDir.*?for my \$input_spec/s,
	"buildTree5 publishes staged inputs before validating its input paths");
like($source,
	qr/sub publishStagedTreeInputs.*?unless \(\@missing\).*?Using existing persistent tree inputs.*?opendir.*?sortFastaForCompression.*?move\(\$source, \$destination\).*?Tree inputs remain incomplete/s,
	"persistent inputs take precedence and staged publication is validated in Perl");
like($source,
	qr/length\(\$tmpSubdir\).*?\$ENV\{TMPDIR\}.*?File::Spec->catdir\(\$temporaryRoot.*?prepareTemporaryBase/s,
	"TMPDIR-relative work paths are resolved inside buildTree5");
like($source,
	qr/safeRemoveTree\(\$tmpD, \$tmpBase\).*?writeCompletionMarker\(\$completionMarker, \$\{\$trRetH\}\{nwk\}.*?sub writeCompletionMarker.*?nonempty primary tree.*?rename \$temporaryMarker, \$marker/s,
	"buildTree5 atomically publishes a completion marker only after validating its primary tree");
like($source,
	qr/post_alignment_locus_qc\.tsv.*?sub runPostAlignmentLocusQC.*?getProgPaths\("MSAfix"\).*?-manifest.*?-report.*?-keep/s,
	'buildTree invokes native MSAfix locus QC before concatenation and retains its report');
unlike($source, qr/postAlignmentLocusQC_scr/,
	'buildTree no longer invokes the Perl locus-QC script');
like($source, qr/post-alignment-loci-XXXXXX.*?UNLINK => 1.*?post-alignment-keep-XXXXXX.*?UNLINK => 1/s,
	'locus-QC manifest and keep-list temporaries are always scheduled for cleanup');
like($source, qr/my \@temporaryFiles = \(.*?bsd_glob\(quotemeta\(\$reportFile\)\."\.tmp\.\*"\).*?bsd_glob\(quotemeta\(\$keepFile\)\."\.tmp\.\*"\).*?unlink \$temporaryFile/s,
	'wrapper and partial native locus-QC files are explicitly deleted after every invocation');
like($source,
	qr/my %POST_ALIGNMENT_QC_DEFAULT = \(.*?enabled => 1.*?minimum_occupancy => 0\.35.*?relative_modified_z => 8\.0/s,
	'post-alignment QC is enabled with permissive metagenomic defaults');
like($source,
	qr/existing multi-locus alignment predates post-alignment.*?safeRemoveTree\(\$MsaD.*?safeRemoveTree\(\$treeD/s,
	'a legacy concatenated checkpoint is rebuilt once when its locus-QC audit is absent');
like($source,
	qr/\@MSAs = grep \{ \$keepPath\{\$_\} \} \@MSAs.*?\@MSAsSyn = grep \{ \$keepStem\{alignmentFileStem\(\$_\)\} \} \@MSAsSyn/s,
	'primary, synonymous, and nonsynonymous alignment sets stay locus-consistent');
like($source,
	qr/"iqMemMB=i" => \\\$iqMemMB.*?"iqPathogen=i" => \\\$iqPathogen.*?"iqLegacy=i" => \\\$iqLegacy/s,
	'buildTree exposes memory-capped pathogen and legacy IQ-TREE modes');
like($source, qr/-iqPathogen and -iqLegacy are mutually exclusive/,
	'buildTree rejects conflicting modern and legacy IQ-TREE modes');
like($source,
	qr/iqMemMB => \$iqMemMB.*?iqPathogen => \$iqPathogen.*?iqLegacy => \$iqLegacy/s,
	'buildTree forwards IQ-TREE execution controls to phyloTools');
like($source,
	qr/BuildTree pipeline v\$version.*?Inputs:.*?Paths:.*?Mode:.*?Alignment:.*?Filtering:.*?Trees:.*?Additional analyses:/s,
	'buildTree starts with a structured runtime configuration header');
like($source, qr/sub limitedWarn.*?No more '\$category' warning examples/s,
	'buildTree caps repetitive warning examples');
like($source, qr/END \{.*?Suppressed \$suppressed additional/s,
	'buildTree reports suppressed warning totals');
like($source, qr/Per-locus alignment summary:.*?Synonymous-site classification summary:/s,
	'buildTree reports aggregate alignment and site-classification progress');
like($source, qr/Alignment merge summary:.*?Overlap filtering summary:/s,
	'buildTree consolidates sequence and overlap filtering diagnostics');
unlike($source, qr/print \$cmd\."\\n"|print \$cmd\."\\n\\n"/,
	'buildTree does not echo routine execution commands');
like($source, qr/-outD is required/, 'an output directory is explicitly required');
like($source, qr/Refusing to use filesystem root/, 'filesystem roots are rejected as output directories');
like($source, qr/buildTree5_\$\{tmpTag\}_\$\$/, 'work is isolated in a process-owned temporary directory');
like($source,
	qr/prepareTemporaryBase\(\$tmpBase\).*?Requested temporary path is unusable:.*?falling back to \$fallbackTmpBase.*?prepareTemporaryBase\(\$tmpBase\)/s,
	'an unusable requested temporary path falls back to output-local workspace');
like($source,
	qr/sub prepareTemporaryBase .*?tempfile\(.*?DIR => \$path.*?print \{\$probeHandle\}.*?unlink \$probePath/s,
	'a temporary base must pass a create, write, close, and cleanup probe');
like($source,
	qr/my \$reusableAlignment = \$isAligned \|\| \(.*?fileGZe\(\$multAli\).*?if \(\$continue\).*?\$treesDone.*?\$reusableAlignment.*?no reusable alignment or complete tree checkpoint.*?safeRemoveTree\(\$MsaD.*?safeRemoveTree\(\$treeD/s,
	'continue mode retains only validated checkpoints and restarts incomplete alignment/tree stages');
like($source, qr/safeRemoveTree\(\$tmpD, \$tmpBase\)/, 'cleanup is limited to the owned temporary directory');

unlike($source, qr/touch \$IQtreef/, 'an empty IQ-TREE checkpoint is not manufactured');
unlike($source, qr/\$calcSyn\s*=\s*0\s*;\s*\$calcNonSyn\s*=\s*0\s*;\s*\n\s*if \(\$cogCats/, 
	'synonymous and nonsynonymous tree options are retained');
like($source, qr/my \$lengthInNt = .*\? \$totalNTs\{\$sp\} : \$totalNTs\{\$sp\} \* 3;/,
	'NTfiltCount compares a nucleotide-equivalent length');
like($source, qr/systemW\(\$cmd1\."\\n"\.\$cmd2\."\\n"\)/,
	'alignment and post-filter commands use checked execution');
like($source, qr/MSA command completed without producing/, 'alignment output is verified');
like($source, qr/runMSAFix\(\$tmpOutMSA, \$maxGapPerCol\)/,
	'per-locus nucleotide alignments use the guarded MSAfix path');
like($source, qr/runMSAFix\(\$multAli, \$maxGapPerCol\)/,
	'single-gene nucleotide alignments use the guarded MSAfix path');
like($source,
	qr/sub runMSAFix.*?\$tmpOutput = "\$alignment\.MSAfix\.\$\$\.fna".*?"-o", shellQuote\(\$tmpOutput\).*?if \(!-s \$tmpOutput\).*?rename \$tmpOutput, \$alignment/s,
	'MSAfix writes a nonempty sibling temporary file before atomically replacing its input');
like($source, qr/if \(!\$ok\).*?unlink \$tmpOutput if -e \$tmpOutput.*?die \$error/s,
	'a failed MSAfix attempt removes its partial output and preserves the original alignment');
like($source,
	qr/my \$ntAlignmentOK = eval \{.*?runMSAFix\(\$tmpOutMSA, \$maxGapPerCol\).*?if \(!\$ntAlignmentOK\).*?excluding locus \$gene from future calculations.*?next;/s,
	'a failed per-locus MSAfix or nucleotide conversion warns and excludes only that locus');
like($source,
	qr/my \$msaCommandOK = 1;.*?eval \{.*?systemW\(\$cmd1\."\\n"\.\$cmd2\."\\n"\).*?failed locus alignment.*?next;/s,
	'a failed per-locus aligner command does not terminate the multi-locus run');
like($source,
	qr/Per-locus alignment summary:.*?\$failedLoci failed and were excluded/s,
	'failed locus exclusions are included in the alignment summary');
like($source,
	qr/my \$distanceOK = eval \{.*?failed optional locus distance matrix.*?retaining locus \$gene/s,
	'an optional distance-matrix failure retains the successfully aligned locus');
like($source,
	qr/my \@unequal = grep.*?invalid locus MSA.*?excluding alignment \$MSAf during merge.*?next;/s,
	'a malformed unequal-length locus alignment is skipped before concatenation');
like($source,
	qr/\$excludedLoci\{\$gene\} = 1.*?\$excludedLoci\{\$geneF\}.*?skipping previously excluded locus/s,
	'a failed alignment locus is excluded from later fastGEAR processing');
like($source,
	qr/my \$phylipOK = eval \{.*?failed optional per-locus PHYLIP conversion.*?next;/s,
	'an optional per-locus PHYLIP conversion failure does not abort the run');
like($source,
	qr/my \$subtreeOK = eval \{.*?failed locus subtree.*?next;.*?No usable locus subtrees remain/s,
	'individual subtree failures are skipped while an impossible empty supertree remains fatal');
like($source,
	qr/my \$fastgearOK = eval \{.*?failed fastGEAR locus.*?next;/s,
	'an individual fastGEAR tool failure does not abort other loci');

like($source, qr/\$pigzBin -d .*\$partiF\.gz/, 'compressed partition restoration names the gzip file');
like($source,
	qr/if \(\$gzipInput\).*?basename\(\$inputFile\).*?sortFastaForCompression\(\$inputFile\).*?allFAAs\.faa.*?allFNAs\.fna.*?\$pigzBin -p \$ncore/s,
	'buildTree sorts only the named plain FNA/FAA inputs immediately before compressing them');
like($source,
	qr/sub fastaCompressionSortKey.*?parseSeqId\(\$identifier, "compression-sort FASTA header",1\).*?join\("\\t", \$gene, \$sample, \$identifier\).*?sub sortFastaForCompression.*?tempfile\(.*?DIR => dirname\(\$inputFile\).*?rename \$tmpFile, \$inputFile/s,
	'compression sorting orders FASTA records locus-first and replaces the input atomically');
unlike($source, qr/\$partiF\s*=\s*""\s+unless\s*\(-e \$partiF\)/,
	'a fresh multi-locus run does not discard its not-yet-created partition path');
like($source, qr/my \$partition = \$treeOpts\{partition\} \/\/ "";.*?\$treeOpts\{partition\} = "" unless \$partition ne "" && -s \$partition;/s,
	'the partition path is resolved after alignment concatenation, immediately before tree execution');
unlike($source, qr/\$continue && -e (?:\$treeOpts\{(?:fastTrOut|VfastTrOut|RAXNGtreeout|RAXtreeout)\}|"\$IQtree\.treefile")/,
	'resume gates do not accept empty tree outputs');
like($source, qr/\$continue && -s "\$IQtree\.treefile"/,
	'IQ-TREE resume requires a nonempty tree');
like($source, qr/unlink \$treeOpts\{RAXtreeout\}.*?if -e \$treeOpts\{RAXtreeout\} && !-s \$treeOpts\{RAXtreeout\}/s,
	'an empty legacy RAxML tree cannot suppress continuation recovery');
like($source,
	qr/filter_alignment_by_overlap\(\\%MFAA, \$isAA, \$minOverlapMSA\).*?push\(\@lengthsParts,\$len\)/s,
	'minimum taxon overlap is applied per locus before partition lengths are recorded');
like($source, qr/for my \$disM \(\@subfls\)/, 'all discovered distance matrices are merged');
like($source, qr/\$ffd\{\$k\} = 4/, 'fourfold degeneracy is classified by codon family');

for my $environment_name (qw(
	MF4_GUBBINS_BIN MF4_CLONALFRAMEML_BIN MF4_FASTGEAR_BIN
	MF4_FASTGEAR_MATLAB_BIN MF4_FASTGEAR_PARAM_FILE
)) {
	like($source, qr/\Q$environment_name\E/, "$environment_name documents dormant-tool reactivation");
}

done_testing();
