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

like($source, qr/my \$version = 5\.13;/,
	'buildTree per-locus fault isolation increments the workflow version');
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
