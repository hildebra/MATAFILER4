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

like($source, qr/my \$version = 5\.30;/,
	'default taxon-aware locus selection increments the workflow version');
like($source,
	qr/"strainWithinPreset=i".*?if \(\$strainWithinPreset\) \{.*?\$withinSpecies = 1;.*?\$useAA4tree = 0;.*?\$ntCntTotal = 400;.*?\$strictBackbone = 1;.*?\$continue = 1;.*?\$doDNDS = 0;.*?\$doTheta = 0;/s,
	"buildTree strain preset owns the fixed strain-tree settings and within-species mode");
like($source,
	qr/my \$withinSpecies = 0;.*?my \$minOverlapMSA;.*?"withinSpecies=i".*?\$minOverlapMSA = \$withinSpecies \? 2 : 0 unless defined \$minOverlapMSA;.*?\$postAlignmentLocusQC = \$withinSpecies.*?unless defined \$postAlignmentLocusQC;.*?\$postAlignmentDivergenceQC = \$withinSpecies \? 1 : 0/s,
	"between-species locus retention is the default and within-species filtering is explicit");
like($source,
	qr/my \@divergenceArguments = \$postAlignmentDivergenceQC.*?"-maxMedianDivergence", 1,.*?"-maxP90Divergence", 1,.*?"-relativeModifiedZ", 1_000_001/s,
	"explicit broad-tree QC disables strain-divergence rejection");
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
	qr/my %POST_ALIGNMENT_QC_DEFAULT = \(.*?between_species_enabled => 0.*?within_species_enabled => 1.*?minimum_occupancy => 0\.35.*?relative_modified_z => 8\.0.*?my \$postAlignmentLocusQC;/s,
	'broad trees retain all loci by default while within-species trees retain structural QC');
like($source,
	qr/post_alignment_locus_qc\.policy\.tsv.*?"schema=4".*?"enabled=\$postAlignmentLocusQC".*?"per_gene_length_fraction=\$ntFracGene".*?"minimum_category_q90_fraction=\$fracMaxGenes90pct".*?"minimum_gene_fraction_per_species=\$GeneFracPSpec".*?"taxon_aware=\$taxonAwareLocusSelection".*?\$legacyWithinSpeciesQCAudit = !\$taxonAwareLocusSelection && \$withinSpecies.*?!-e \$postAlignmentQCPolicyFile.*?\$postAlignmentQCAuditCurrent = \$postAlignmentQCPolicyMatches.*?!\$postAlignmentLocusQC.*?existing multi-locus alignment predates the current.*?safeRemoveTree\(\$MsaD.*?safeRemoveTree\(\$treeD/s,
	'changed locus-retention policies rebuild stale checkpoints while legacy within-species audits remain compatible');
like($source,
	qr/sub writePostAlignmentQCPolicy.*?post-alignment-policy-XXXXXX.*?UNLINK => 1.*?rename \$temporaryPolicy, \$policyFile.*?writePostAlignmentQCPolicy\(\$policyFile, \$policyText\)/s,
	'the exact locus-retention policy is published atomically');
like($source,
	qr/elsif \(\$cogCats ne ""\).*?Post-alignment locus QC disabled; retaining all \$candidateCount prepared loci.*?unlink \$postAlignmentQCReport.*?writePostAlignmentQCPolicy\(\$postAlignmentQCPolicyFile, \$postAlignmentQCPolicy\)/s,
	'disabled locus QC retains all prepared alignments and replaces stale QC audit state');
like($source,
	qr/my \$minimumCategorySequences = \$GenesQtl90 \* \$fracMaxGenes90pct;.*?= 1 if \$minimumCategorySequences < 1;.*?if \(\@spl >= \$minimumCategorySequences\).*?minimum category sequences=\$minimumCategorySequences/s,
	'a zero category-prevalence fraction retains every nonempty locus and reports the effective threshold');
like($source,
	qr/\@MSAs = grep \{ \$keepPath\{\$_\} \} \@MSAs.*?\@MSAsSyn = grep \{ \$keepStem\{alignmentFileStem\(\$_\)\} \} \@MSAsSyn/s,
	'primary, synonymous, and nonsynonymous alignment sets stay locus-consistent');
like($source,
	qr/my %TAXON_AWARE_DEFAULT = \(.*?enabled => 1.*?maximum_loci => 500.*?core_loci => 400.*?candidate_extra => 150.*?target_loci_per_sample => 25.*?target_nt_per_sample => 7500.*?"taxonAwareLocusSelection=i"/s,
	'taxon-aware locus selection is enabled by default and has bounded candidate/core defaults');
like($source,
	qr/if \(\$taxonAwareLocusSelection\) \{.*?selectTaxonAwareCandidateLoci\(.*?candidate_limit => \$taxonAwareMaxLoci \+ \$taxonAwareCandidateExtra.*?\@linesCats3 = \@\{\$candidateSelection->\{categories\}\}/s,
	'the taxon-aware pre-MSA pass chooses robust, rescue, and QC-backfill candidates in buildTree5');
like($source,
	qr/sub chooseTaxonAwareLoci.*?robust_core.*?\$coverageGain.*?qc_backfill.*?taxon_rescue/s,
	'locus choice combines a robust core with rarity-weighted taxon rescue and backfill');
like($source,
	qr/runPostAlignmentLocusQC\(.*?if \(\$taxonAwareLocusSelection && \$cogCats ne ""\).*?selectTaxonAwareFinalLoci\(.*?taxonAwareAlignmentMetrics.*?parsimony_informative_sites/s,
	'the final selector runs after MSAfix QC and scores actual occupancy and informative sites');
like($source,
	qr/sub classifyTaxonAwareSamples.*?below_minimum_anchor_nt.*?backbone_candidate.*?placement_candidate/s,
	'sparse but anchored taxa are retained for placement and final sample decisions are audited');
like($source, qr/taxon_aware_locus_selection\.tsv.*?taxon_aware_sample_selection\.tsv/s,
	'taxon-aware final locus and sample decisions have persistent audit tables');
like($source,
	qr/if \(\$taxonAwareLocusSelection && \$multAliF eq \$multAli\).*?\(\$num1 \* \$factor\) < \$minimumAnchorNT.*?else \{.*?\$qtl90NTcnts \* \$ntFrac/s,
	'the final primary merge honors the absolute taxon-aware anchor instead of rerunning relative sample filtering');
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
	qr/my \$primaryAlignmentReady = fileGZe\(\$multAli\).*?my \$siteAlignmentsReady =.*?my \$reusableAlignment = \$isAligned.*?if \(\$continue\).*?\$treesDone.*?\$reusableAlignment.*?no reusable alignment or complete tree checkpoint.*?safeRemoveTree\(\$MsaD.*?safeRemoveTree\(\$treeD.*?my \$calcMSA = !\$treesDone && !\$primaryAlignmentReady.*?\$doMSA = !\(/s,
	'continue mode derives reporting and MSA recovery from one validated alignment state');
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
like($source,
	qr/sub requestedTreeMethods.*?name => "IQ-TREE".*?iqtree => 1.*?sub treeMethodState.*?iqtreeOutputComplete\(\$hr->\{\$method->\{outputKey\}\}, \$hr->\{inMSA\}, \\\$validationReason\).*?checkpointComplete => \(\$continue && \$outputComplete/s,
	'IQ-TREE resume uses the shared successful-log and exact-taxon-parity checkpoint state');
like($source,
	qr/sub treePresent.*?treeMethodState\(\$_, \$hr\).*?sub treeAtHeart.*?my %treeState = map.*?treeMethodState\(\$_, \\%treeOpts\)/s,
	'tree recovery and execution derive completion from the same method-state helper');
like($source,
	qr/my %BACKBONE_DEFAULT = \(.*?coverage_fraction => 0\.35/s,
	'strict strain backbones use the relaxed severe-outlier coverage threshold');
like($source,
	qr/\$tOhr->\{IQtreeout\} \.= "\.backbone" if \$strictBackbone && \$doIQTree/,
	'strict IQ-TREE inference uses a dedicated backbone output prefix');
like($source,
	qr/my \$dedicatedBackbone = \$primaryTree =~ s\/\\\.backbone\\\.treefile\$\/\.treefile\/.*?write_placed_tree\(\$backboneTree, \$primaryTree, \$placements\).*?\$\{\$trRetH\}\{backbone_nwk\} = \$backboneTree;.*?\$\{\$trRetH\}\{nwk\} = \$primaryTree;/s,
	'the placed tree becomes the primary .treefile while the ML tree remains .backbone.treefile');
like($source, qr/unlink \$treeOpts\{RAXtreeout\}.*?if -e \$treeOpts\{RAXtreeout\} && !-s \$treeOpts\{RAXtreeout\}/s,
	'an empty legacy RAxML tree cannot suppress continuation recovery');
like($source,
	qr/filter_alignment_by_overlap\(\\%MFAA, \$isAA, \$minOverlapMSA\).*?push\(\@lengthsParts,\$len\)/s,
	'minimum taxon overlap is applied per locus before partition lengths are recorded');
like($source, qr/for my \$disM \(\@subfls\)/, 'all discovered distance matrices are merged');
like($source, qr/\$ffd\{\$k\} = 4/, 'fourfold degeneracy is classified by codon family');

for my $config_key (qw(
	gubbins clonalframeml fastgear fastgearSummary fastgearReorder
	fastgearMatlab fastgearParam
)) {
	like($source, qr/requireConfiguredTool\("\Q$config_key\E"/,
		"$config_key documents config-backed dormant-tool reactivation");
}

done_testing();
