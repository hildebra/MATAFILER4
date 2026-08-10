use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::GenoMetaAss qw(readClstrRev readFasta);

sub write_file {
	my ($path, $contents) = @_;
	open my $fh, '>', $path or die "Cannot write $path: $!";
	print {$fh} $contents or die "Cannot write $path: $!";
	close $fh or die "Cannot close $path: $!";
}

sub slurp {
	my ($path) = @_;
	open my $fh, '<', $path or die "Cannot read $path: $!";
	local $/;
	return <$fh>;
}

my $tmp = tempdir(CLEANUP => 1);
my $fasta = File::Spec->catfile($tmp, 'records.fa');
write_file($fasta, <<'FASTA');
>keep1 D=8 CSP=0.01
AAAA
>drop D=2 CSP=0.50
CCCC
>keep3 D=5 CSP=0.02
GGGG
FASTA

my %wanted = (keep1 => 1, keep3 => 1);
is_deeply(
	readFasta($fasta, 1, '\\s', \%wanted),
	{ keep1 => 'AAAA', keep3 => 'GGGG' },
	'FASTA subset selection applies independently to intermediate and final records',
);
is_deeply(
	readFasta($fasta, 0, '\\s', \%wanted),
	{
		'keep1 D=8 CSP=0.01' => 'AAAA',
		'keep3 D=5 CSP=0.02' => 'GGGG',
	},
	'FASTA subset lookup can use short IDs while retaining full headers',
);

my $glob_dir = File::Spec->catdir($tmp, 'glob');
mkdir $glob_dir or die "Cannot create $glob_dir: $!";
write_file(File::Spec->catfile($glob_dir, 'a_empty.fa'), '');
write_file(File::Spec->catfile($glob_dir, 'b_records.fa'), ">later\nACGT\n");
is_deeply(
	readFasta(File::Spec->catfile($glob_dir, '*.fa'), 1, '\\s'),
	{ later => 'ACGT' },
	'an empty member of a FASTA glob does not suppress later files',
);

my $cluster_index = File::Spec->catfile($tmp, 'cluster.idx');
write_file($cluster_index, "seed1\tsample1__gene1,sample2__gene2\n");
my (undef, $empty_cluster_subset) = readClstrRev($cluster_index, 0, {}, {});
is_deeply(
	$empty_cluster_subset,
	{},
	'an explicitly empty cluster-member subset does not fall back to the complete catalogue',
);

my $strain = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'MGS', 'strain_within.pl'));
like($strain, qr/sub consensusInputState .*?\$nt_ready && \$aa_ready.*?return 'regenerate' if \$vcf_ready/s,
	'consensus resume requires the paired NT and AA outputs and repairs from VCF');
like($strain, qr/readFasta\(\$fastaf,1,"\\\\s",\\%subG\).*?readFasta\(\$fastafAA,0,"\\\\s",\\%subG\)/s,
	'within-strain extraction reads only candidate consensus genes');
like($strain, qr/-completionMarker "?\.shellQuote\(\$treeStone\)/,
	"tree jobs delegate validated completion markers to buildTree5");
unlike($strain, qr/test -s "?\.shellQuote\(\$IQtreef\).*?touch "?\.shellQuote\(\$treeStone\)/s,
	"tree jobs no longer encode completion validation in shell");
like($strain, qr/my \@staleOutputs = \$record->\{epa_only\}.*?qw\(stone terminal\).*?qw\(stone tree terminal placement_pending\).*?retry_unlink/s,
	'an EPA-only submission retains its placement marker and backbone-derived primary path while ordinary retries clear stale outputs');
like($strain, qr/clear_split_generation\(\$splitManifest.*?write_split_generation\(\$splitManifest.*?printf '%s\\\\n'/s,
	'a new split generation clears stale state and tags each worker completion');
like($strain, qr/ConspecificMGS\.\$subJob\.log.*?sub mergeConspecificLogs/s,
	'split workers write isolated conspecific logs that are explicitly merged');
like($strain, qr/\$onlySubmit == 0 && !\$subJob/,
	'split children cannot recursively clean shared MGS output directories');
like($strain,
	qr/\$publishedInputsReady = !exists\(\$legacyLocusMGS\{\$MGS\}\).*?persistentMGSInputState\(\$MGS\) eq 'complete'.*?if \(\$publishedInputsReady && !\$mustRegenerateInputs\).*?combineMGSgenesDir\(\$MGS,\$tmpD,\$tmpD\)/s,
	'complete published inputs bypass missing scratch aggregates during tree recovery');
like($strain, qr/has neither complete published inputs nor complete combined worker input/,
	'incomplete worker input is reported only when published recovery inputs are also incomplete');
like($strain,
	qr/sub indexRecoveryRow .*?recoveryWorkersByMGS.*?recoveryRecordsByMGS.*?for my \$worker \(0 \.\. \$#parts\).*?indexRecoveryRow\(\$worker/s,
	"recovery rows index the exact workers and record cardinality expected for each MGS");
like($strain,
	qr/exact_worker_parts\(\$prefix, \$workerCount\).*?split_generation_complete\(\$splitManifest.*?\@missing.*?\@unexpected.*?Rejecting merge/s,
	"worker merge requires the complete recovery-ledger contributor set from a completed generation");
like($strain,
	qr/records=\$fnaRows expected=\$expectedRecords.*?QC=.*?expected=\$expectedRows.*?if \(\@validationErrors\).*?retry_unlink\(\$_.*?values %mergeFileByName.*?return \$aggregateComplete.*?retry_unlink\(\$mergeCheckpoint.*?retry_rename\(\$mergeFileByName.*?retry_unlink\(\$part/s,
	"merged record and sample cardinalities are validated before aggregate publication or part deletion");
like($strain,
	qr/my \@coreRequired = \(.*?\$LINKstdof.*?my \@filesets = \(.*?\$QCstdof\.tmp.*?my \@contributorNames = \(.*?\$LINKstdof.*?\$QCstdof\.tmp/s,
	'link and QC files are mandatory members of every worker contribution set');
like($strain,
	qr/sub stagedMGSInputsReady .*?coreRequiredNames.*?CATstdof\.tmp.*?QCstdof\.tmp.*?\|\|.*?CATstdof.*?QCstdof.*?mergeCheckpoint/s,
	'staged readiness requires core merge artifacts, a raw or final category/QC pair, and the last-written commit checkpoint');
like($strain,
	qr/unless \(\$recoveryContributionIndexReady\).*?Rejecting fresh merge.*?return \$aggregateComplete/s,
	'fresh merges without persisted recovery provenance are rejected without discarding a committed aggregate');
like($strain,
	qr/Digest::SHA->new\(256\).*?\$digest->add.*?identifier order differs: \$FNAstdof vs \$name/s,
	'FNA, FAA, category, and link streams must contain the same identifiers in the same order');
like($strain,
	qr/loadRecoveryContributionIndex\(\) unless \$recoveryContributionIndexReady.*?sub writeRecoveryContributionIndex.*?retry_rename\(\$temporary, \$path.*?sub loadRecoveryContributionIndex.*?\$recoveryContributionIndexReady = 1.*?sub mergeRecoveryLogs.*?writeRecoveryContributionIndex\(\);.*?retry_rename\(\$temporary, \$final.*?retry_unlink\(\$_.*?for \@parts/s,
	'worker provenance is persisted before disposable recovery logs are removed and is reloadable after restart');
like($strain,
	qr/\$aggregateComplete &&= -s \$mergeCheckpoint.*?my \@expectedWorkers.*?if \(\@validationErrors\).*?return \$aggregateComplete.*?retry_unlink\(\$mergeCheckpoint.*?retry_rename\(\$checkpointTemporary, \$mergeCheckpoint.*?retry_unlink\(\$part/s,
	'a last-written checkpoint protects committed aggregates and worker parts survive until validation and commit');
like($strain, qr/qsubSystemJobAlive\([^\n]+QSBoptHR[^\n]+if [^\n]+doSubmit/,
	'dry runs do not poll scheduler jobs that were never submitted');
like($strain, qr/\$nxtCmd \.= "-submit \$doSubmit ";.*?-qsubSystem/s,
	'postprocessing inherits submission state and the selected queue backend');
like($strain, qr/sub assertSafeWorkflowRemoval .*?resolved_default.*?Refusing to remove unowned custom output directory/s,
	'custom recursive output removal requires a workflow-owned directory');
like($strain, qr/sub limitedWarn .*?warningExampleLimit.*?Further '\$category' warnings are suppressed/s,
	'repetitive strain warnings retain examples and announce suppression');
like($strain, qr/Suppressed warning summary:.*?sort grep/s,
	'suppressed strain warnings receive a categorized exit summary');
unlike($strain, qr/print "\$cD\\n"/,
	'strain extraction no longer prints a raw working-directory path for every sample');
like($strain, qr/my \$version = 0\.99;/,
	'prepared Phase-II scratch reuse increments the workflow version');
like($strain,
	qr/my \@sampleStatColumns = sample_stat_columns\(\);.*?GetOptions\(.*?printEarlyRunHeader\(\)/s,
	'sample-statistics columns are initialized before the executable workflow begins');
unlike($strain, qr/print STDERR "\nAT SMPL::/,
	"sample progress does not emit a leading blank line per assembly group");
like($strain,
	qr/printEarlyRunHeader\(\);.*?read_mosaic_catalogue\(.*?prepRun\(\)/s,
	'the autoflushed basic header is emitted before Mosaic, map, and catalogue loading');
like($strain,
	qr/sub printEarlyRunHeader \{.*?Strain_within v\$version.*?Started:.*?Requested output:.*?Initializing paths, maps, and catalogues/s,
	'the immediate header identifies the run before expensive initialization starts');
like($strain,
	qr/print \{\$sampleStatsFH\} \$sampleStatsHeader.*?my %sampleStatsSeen;.*?local \*STDOUT;.*?open STDOUT.*?STDERR.*?foreach my \$sm \(\@srtdSmpls\).*?readGenesSample_Singl/s,
	'STDOUT is redirected once around the complete sample loop while the duplicated handle carries only TSV records');
like($strain,
	qr/sub writeSampleStats \{.*?without a sample name.*?duplicate row.*?Refusing to emit an empty.*?for my \$target \(\$fh, \$sampleStatsPartFH\).*?print \{\$target\} \$row, "\\n"/s,
	'every sample-statistics record has a sample name, is nonempty, and is emitted to stdout and its worker table at most once');
like($strain,
	qr/sub mergeSampleStats .*?Wrong sample-statistics field count.*?Duplicate sample-statistics row.*?aggregate_sample_rows.*?STEP 1 SAMPLE SUMMARY \(all workers\)/s,
	'all worker tables are validated, aggregated, saved, and reported at the end of Step 1');
like($strain,
	qr/STEP 1 SAMPLE SUMMARY \(all workers\).*?join\(" ", \@summaryPairs\).*?loci_histogram_rows.*?Used MGS retained-loci histogram/s,
	'the all-worker stdout summary uses key:value pairs and includes a retained-locus histogram');
like($strain,
	qr/mergeRecoveryLogs\(\) unless \$maxSubJob.*?mergeSampleStats\(\) unless \$maxSubJob.*?if \(\$maxSubJob && !\$subJob\).*?mergeRecoveryLogs\(\);.*?mergeSampleStats\(\);/s,
	'both single-worker and split-worker extraction produce the combined sample summary');
like($strain,
	qr/\$phase1SelfCmd -subjob \$sj &&\\n.*?write_worker_completion/s,
	'split workers publish completion only after sample statistics and extraction finish successfully');
like($strain,
	qr/my \$mosaicDirectory = File::Spec->catdir\(dirname\(\$mosaicMGSFile\), 'mosaic'\).*?basename\(\$mosaicMGSFile\)\."\.mosaic_loci\.\$clusterID\.confirmed\.tsv".*?prepare_mosaic_loci\.log/s,
	'a missing default Mosaic catalogue is named from the raw MGS table and uses a temporary job log');
like($strain,
	qr/qsubSystem\(.*?"MosaicMGS".*?qsubSystemJobAlive\(\[\$mosaicDependency\].*?Prerequisite Mosaic catalogue is ready/s,
	'Mosaic is submitted as a prerequisite and awaited before strain work continues');
like($strain,
	qr/unless \(\$doSubmit\).*?stopping before Mosaic-dependent strain extraction.*?exit 0/s,
	'a no-submission run generates the Mosaic script without consuming absent results');
like($strain,
	qr/Reusing existing confirmed Mosaic catalogue: \$mosaicLociFile.*?if \(length\(\$mosaicLociFile\) && !-s \$mosaicLociFile\).*?Raw MGS assignment file for Mosaic is missing or empty/s,
	'an existing confirmed Mosaic catalogue bypasses prerequisite generation and raw-MGS input requirements');

like($strain,
	qr/Mosaic outgroup \$source -> \$PreferredOutgroup\{\$source\}.*?Mosaic outgroup proposals loaded:.*?unique MGS-to-MGS connection.*?gene-to-gene link/s,
	'strain workflow reports loaded outgroup connections and proposed gene links');
like($strain,
	qr/Loading confirmed Mosaic catalogue for split worker.*?Using Mosaic catalogue:.*?next if \$subJob/s,
	'split workers load and summarize Mosaic data without connection-by-connection previews');
like($strain,
	qr/sub stepComplete .*?STEP COMPLETE: \$step/s,
	'step completion messages use one consistent formatter');
like($strain,
	qr/configuration and map initialization.*?assembly-group expansion.*?MGS and seed-locus selection.*?existing-output and resume audit.*?if \(\$runPartI\).*?consensus-input audit/s,
	'startup stages emit consistent completion messages, with the consensus audit limited to Phase I');
like($strain,
	qr/locus-model construction.*?catalogue_drivers=.*?resolved_loci=.*?consensus-gene extraction and publication.*?tree-input sizing/s,
	'major extraction and tree-preparation stages also report concise completion statistics');
like($strain,
	qr/historical exclusion loading.*?excluded_MGS=.*?outgroup-reference preparation.*?reference_NT=.*?MGS_with_outgroup_candidates=/s,
	'historical exclusions and outgroup-reference preparation report their final counts');
like($strain,
	qr/my %treeDisposition.*?\$treeDisposition\{\$epaOnlyRetry \? 'EPA-only retry job' : 'eligible tree job'\}\+\+.*?Tree submission accounting:.*?Tree submission pass complete:/s,
	'tree submission reports every eligible and skipped MGS disposition before waiting');
like($strain,
	qr/my \@pendingTreeJobs;.*?push \@pendingTreeJobs, \{.*?command => \$Tcmd\.\$outgS.*?tmp_space => \$QSBoptHR->\{tmpSpace\}.*?dispatchPendingTreeJobs\(.*?blocking => 0.*?Tree preparation pass complete:.*?dispatchPendingTreeJobs\(.*?blocking => 1.*?qsubSystemJobAlive\( \\\@jobs.*?writeTreeFailureAudit.*?without a valid output were quarantined/s,
	'eligible trees queue after conversion, drain opportunistically under capacity, then are tracked, awaited, and output-validated');
like($strain, qr/nonblockingMaxConcurrentJobs\} = 1 unless \$blocking/,
	'queued tree dispatch uses the non-blocking scheduler-capacity path');
like($strain, qr/deferredSubmissionDependency\(\).*?Phase II continues converting inputs/s,
	'queued tree dispatch retains deferred jobs while conversion continues');
like($strain, qr/\$options->\{tmpSpace\} = \$record->\{tmp_space\}/,
	'queued tree dispatch restores each job\'s stored temporary-space setting');
like($strain, qr/\$options->\{useLongQueue\} = \$record->\{use_long_queue\}/,
	'queued tree dispatch restores each job\'s stored queue setting');
like($strain,
	qr/\@treeJobAccounting.*?requested_mb => int\(\$totMem\).*?qsubSystemJobAlive.*?slurm_tree_memory_summary.*?format_slurm_tree_memory_summary/s,
	'completed Slurm tree jobs report MaxRSS against their requested memory');
like($strain,
	qr/append_fasta_records_atomic\(\$FNAtf.*?append_fasta_records_atomic\(\$FAAtf.*?if \(\$temporaryInput\).*?sort_fasta_by_locus\(\$FNAtf, \$SaSe\).*?sort_fasta_by_locus\(\$FAAtf, \$SaSe\)/s,
	'first-generation FNA and FAA inputs are locus-sorted after outgroup publication');
like($strain,
	qr/The following wait count reports jobs still present, not jobs omitted/,
	'the scheduler wait count is explicitly distinguished from submission coverage');
like($strain,
	qr/my \$treeJobOrdinal = \$cnt \+ 1;.*?"FT\$treeJobOrdinal"/s,
	'tree scheduler labels use one-based submission ordinals');
like($strain,
	qr/END \{.*?Suppressed warning summary:.*?Repeated status summary:.*?FATAL: strain_within\.pl terminated:.*?FINISH:/s,
	'shutdown summaries precede a final fatal or successful completion diagnostic');
like($strain,
	qr/\$completionMessage = "strain_within\.pl completed normally;.*?exit\(0\)/s,
	'the regular main-process exit records an explicit FINISH message');
like($strain,
	qr/my \$iqPathogen = 0.*?"iqPathogen=i"\s+=> \\\$iqPathogen.*?\$Tcmd .= "-iqPathogen 1 " if \$iqPathogen/s,
	'within-strain pathogen mode defaults off and is applied only by the parent tree command');
like($strain,
	qr/my \$iqMemMB = int\(\$totMem \* 0\.9\).*?if \(\$phyloProg == 1\)\{.*?"-iqMemMB \$iqMemMB ".*?"-iqPathogen 1 " if \$iqPathogen/s,
	'within-strain IQ-TREE always uses the standard resource-limited command and enables CMAPLE only by explicit request');
like($strain,
	qr/my \$placementRequested = \$strictBackbone \? 1 : 0;.*?\$baseMemMult = 150 if \$placementRequested.*?\$minimumMemMB = \(\$placementRequested \? 10240 : 5000\) \* \$memMulti;.*?\$minimumMemMB = 10240 if \$placementRequested.*?\$totMem = \$minimumMemMB if \$totMem < \$minimumMemMB/s,
	'within-strain gives EPA-ng placement jobs a 10 GiB floor and larger input-size estimate');
unlike($strain, qr/-iqLegacy\s+1|legacyMGTK/,
	'within-strain does not expose or submit the obsolete IQ-TREE legacy-kernel flag');
like($strain,
	qr/"recalcTrees=i"\s+=> \\\$recalcTrees.*?-recalcTrees must be 0 or 1.*?\$onlySubmit = 1 if \$recalcTrees/s,
	'tree recalculation is validated and forced into input-recovery-only mode');
like($strain,
	qr/-recalcTrees cannot be combined with -repairCAT, -deepRepair, or -redoSubmissionData.*?-recalcTrees must be launched by the main strainWithin process/s,
	'tree recalculation rejects input-regeneration modes and split-worker execution');
like($strain,
	qr/my \$runPartI = \(.*?\|\| \(\$recalcTrees && \$dirsNOTPrepped\).*?if \(\$runPartI\).*?Part I:: extracting relevant core MGS genes/s,
	'tree recalculation reruns extraction when required per-MGS inputs are absent');
like($strain,
	qr/my \$runPartI = .*?if \(\$runPartI\).*?preComputeConsSNP\(\).*?\} else \{.*?Skipping Part I.*?reportSavedSampleStats\(\)/s,
	'completed Phase I skips the extraction-only consensus audit and reports saved statistics');
like($strain,
	qr/sub reportSavedSampleStats .*?\$sampleStatsSummaryLogName.*?scope\} .*?eq 'ALL'.*?printSampleStatsSummary\(\$allSummary\)/s,
	'a Phase-I resume loads the persisted all-worker row before rendering its accounting');
like($strain,
	qr/sub printSampleStatsSummary .*?STEP 1 SAMPLE SUMMARY \(all workers\).*?Used MGS retained-loci histogram/s,
	'the same all-worker summary and retained-locus histogram are available after Phase I has completed');
like($strain,
	qr/sub recoverCompletedSplitPhaseI .*?split_generation_complete.*?incomplete recovery ledgers.*?incomplete sample-statistics ledgers.*?mergeConspecificLogs\(\).*?mergeRecoveryLogs\(\).*?mergeSampleStats\(\)/s,
	'a restart after every split worker finished merges validated ledgers before Phase II uses their staged inputs');
ok(index($strain, '$dirsNOTPrepped == 0') >= 0
	&& index($strain, 'tree-only resume skips obsolete Phase-I ledger validation') >= 0
	&& index($strain, 'every MGS input passed the completed audit') >= 0
	&& index($strain, 'continuing to Phase II') >= 0,
	'legacy runs with every MGS input already reusable do not rebuild Phase I solely because historical worker ledgers are absent');
like($strain,
	qr/next if \$recalcTrees && !\$MGSneedsExtraction\{\$MGS\}.*?\$MGSneedsExtraction\{\$MGS\} = 1/s,
	'tree recalculation limits its extraction model to MGS with missing inputs');
like($strain,
	qr/sub stagedMGSInputsReady .*?aggregateComplete.*?hasFreshParts.*?split_generation_complete.*?return 0 if grep.*?stagedMGSInputsReady\(\$MGS\)/s,
	'the resume audit accepts only a complete staged FNA/FAA/category set');
like($strain,
	qr/sub persistentMGSInputState .*?persistentMGSInputStateCache.*?\$FNAstdof, \$FAAstdof, \$CATstdof.*?my \$state = .*?'complete'.*?'incomplete'.*?sub scratchMGSInputState .*?scratchMGSInputStateCache.*?stagedMGSInputsReady.*?'complete'/s,
	'published reuse requires the complete FNA/FAA/category triplet while complete Stage-I staging remains reusable');
like($strain,
	qr/sub stagedMGSInputsReady .*?return 1 if \$aggregateComplete;.*?exact_worker_parts/s,
	'a committed staged aggregate avoids repeated worker-part directory scans');
ok(index($strain, 'my $preparedScratchInput') >= 0
	&& index($strain, 'merge.complete.tsv') >= 0
	&& index($strain, 'Stage-I input: reusing prepared scratch tree inputs') >= 0
	&& index($strain, 'return (scalar(keys %samples_seen), $genes_seen, $preparedOG, 1, 1);') >= 0,
	'fully prepared Phase-II scratch inputs remain resumable without repeated outgroup, sorting, or category work');
ok(index($strain, 'sub preparedOutgroupLog') >= 0
	&& index($strain, 'fileGZe($log_path)') >= 0
	&& index($strain, '$publishedPrepared') >= 0
	&& index($strain, '$scratchPrepared') >= 0
	&& index($strain, 'data.log.write.$$') >= 0
	&& index($strain, 'remove stale prepared outgroup log') >= 0,
	'outgroup preparation is committed atomically, reused only after that commit, and invalidated by fresh Stage-I input');
like($strain,
	qr/my \(%persistentMGSInputStateCache, %scratchMGSInputStateCache\).*?sub invalidateMGSInputState .*?delete \@persistentMGSInputStateCache.*?delete \@scratchMGSInputStateCache/s,
	'published and scratch triplet states are cached and explicitly invalidated after mutations');
like($strain,
	qr/my \$completedTree = .*?treeDone\.sto.*?fileGZs\(\$completedTree\).*?Avoid opening and decompressing its category sidecar again.*?next;.*?fileGZe\("\$SIdirs\{\$MGS\}\/\$CATstdof"\)/s,
	'a validated completed tree bypasses compressed category-sidecar inspection');
my ($quickWorkerValidation) = $strain =~
	/(sub validatePhase1WorkerLedger .*?)(?=sub phase1WorkersNeedingRetry)/s;
ok(defined($quickWorkerValidation),
	'Phase-I worker prevalidation is available for resume repair');
unlike($quickWorkerValidation, qr/while\s*\(/,
	'Phase-I worker prevalidation checks stones and headers without rescanning ledger rows');
like($strain,
	qr/sub mergeSampleStats .*?while \(my \$line = <\$in>\).*?Wrong sample-statistics field count.*?sub mergeRecoveryLogs .*?Unexpected MAG recovery header.*?while \(my \$line = <\$in>\) \{ indexRecoveryRow/s,
	'deep row and cardinality checks remain in the single merge pass');
like($strain,
	qr/sub resolveScratchDirectory .*?Reusing recorded scratch directory.*?sub persistScratchDirectory .*?retry_rename\(\$temporary, \$manifest.*?\.strain_within\.scratch\.tsv.*?resolveScratchDirectory\(\$derivedScratch.*?if \(\$subJob\).*?return;.*?persistScratchDirectory\(\$scratchManifest/s,
	'main and worker resumes restore an atomically persisted catalogue/output-bound scratch directory');
like($strain,
	qr/my \$publishedInputState = persistentMGSInputState\(\$MGS\).*?if \(\$publishedInputState ne 'complete'\).*?stagedMGSInputsReady\(\$MGS\).*?\$MGSneedsExtraction\{\$MGS\} = 1/s,
	'incomplete published or scratch triplets are marked for extraction without discarding a complete staged recovery set');
like($strain,
	qr/tree_input_sizing\.tsv.*?too_few_samples.*?incomplete_published.*?incomplete_scratch.*?empty_extraction/s,
	'tree-input sizing separates too-few, incomplete published, incomplete scratch, and empty extraction inputs');
like($strain,
	qr/sub recordValidatedEmptyExtractions.*?persistentMGSInputState\(\$MGS\) eq 'missing'.*?scratchMGSInputState\(\$MGS\) ne 'missing'.*?writeNoRecoverableLociMarker\(\$SIdirs\{\$MGS\}, 'empty_extraction'\).*?\$MGSnoTreeReason\{\$MGS\} = 'no_recoverable_loci'/s,
	'a completed Stage I persists validated no-recoverable-locus outcomes for future resumes');
like($strain,
	qr/\$multiSmpl > 2 && \$ngenes >= \$MGStoolowGsThr.*?too_few_usable_genes.*?writeTooFewMarker.*?sub validateTreeInputResolution.*?tree_input_resolution\.tsv.*?repair_required.*?tree_input_repair\.queue\.tsv.*?no catalogue-wide abort was triggered/s,
	'insufficient tree inputs are terminally marked while incomplete triplets enter a persistent repair queue');
like($strain,
	qr/last if \$cntShrCogs >= \$MGStoolowGsThr.*?if \(\$cntShrCogs < \$MGStoolowGsThr\)/s,
	'outgroup and tree eligibility share the configured eight-locus default');
like($strain,
	qr/my \$workerMGSSubset = \$recalcTrees.*?grep \{ \$MGSneedsExtraction\{\$_\} \} \@specis.*?'-MGSsubset', \$workerMGSSubset/s,
	'split extraction workers inherit the missing-input MGS subset');
like($strain,
	qr/Stage-I extraction scope: \$stageIScope.*?target_MGS=.*?Workers are balanced by assembly group/s,
	'split Stage I reports whether its MGS scope is explicit or recovery-driven');
like($strain,
	qr/my \$maxSubJob = -1;.*?choose_auto_worker_count\(.*?Automatic Stage-I splitting:.*?target \$\{targetGroupsPerWorker\} groups\/worker/s,
	'automatic Stage-I splitting is the default and reports its selected granularity');
like($strain,
	qr/'-submit', 0, '-onlySubmit', 1.*?'-MGSphylo', \$treeFile.*?'-flushEvery'.*?'-MGset', \$useGTDBmg/s,
	'extraction workers receive only extraction and outgroup inputs, not tree-submission behavior');
unlike($strain,
	qr/'-rateMergePartitions', \$rateMergePartitions.*?'-iqPathogen', \$iqPathogen.*?'-rmMSA', 0/s,
	'extraction-worker commands do not forward buildTree5-only model and MSA options');
like($strain,
	qr/my \$treeTmpGb = int\(.*?\$QSBoptHR->\{tmpSpace\} = \$nodeTmpConfigured \? \$treeTmpGb : 0.*?\? "-tmpSubdir ".*?strain_within\/\$MGS.*?: "-tmpD "/s,
	'tree jobs request and use node-local scratch when it is configured');
like($strain,
	qr/my \$publishedInputsReady = !exists\(\$legacyLocusMGS\{\$MGS\}\).*?persistentMGSInputState\(\$MGS\) eq 'complete'.*?if \(\$recalcTrees\).*?unless \(\$publishedInputsReady\).*?\$scratchInputsReady = combineMGSgenesDir\(\$MGS,\$tmpD,\$tmpD\).*?unless \(\$publishedInputsReady \|\| \$scratchInputsReady\).*?no recoverable inputs for recalculation.*?resetMGSTreeOutputs\(\$outD2, \$MGS\)/s,
	'tree outputs are reset only after complete published or recoverable staged per-MGS inputs are verified');
like($strain,
	qr/\$scratchInputsReady \|\|= combineMGSgenesDir\(\$MGS,\$tmpD,\$tmpD\).*?Stage-I input: using complete scratch FNA\/FAA\/category files.*?tree job will publish them to the MGS directory/s,
	'recalculated trees continue staged-input transformation and publication through the normal tree-job path');
like($strain,
	qr/staged input sets recovered for -recalcTrees: \$recalcScratchRecovered/,
	'tree submission accounting reports staged recalculation recovery separately from skipped dispositions');
like($strain,
	qr/sub resetMGSTreeOutputs .*?dirname\(\$resolvedMGS\) eq \$resolvedRoot.*?basename\(\$resolvedMGS\) eq \$MGS.*?remove_tree\(\$phyloDir, \{safe => 1\}\).*?retry_unlink\(\$treeStone/s,
	'tree-only reset is confined to the selected MGS phylo directory and completion checkpoint');
like($strain,
	qr/my \$locCl2G2 = \$cl2gene2\{\$sm\}.*?my \$COGprios1 = \$COGprios->\{\$MGS\}.*?\@candidates == 1.*?reason => 'unique'.*?\$LocusSeedProteins\{\$locus\} \|\|=.*?choose_locus_candidate/s,
	'within-strain extraction avoids hot-loop container copies and scoring unique candidates');
like($strain,
	qr/include_member_to_seed => 0.*?include_gene_to_locus => 0/s,
	'within-strain extraction omits unused locus indexes');
like($strain,
	qr/Only identifiers are needed.*?readFastaIDs\(\$resolvedFNA\).*?sub readFastaIDs/s,
	'within-strain outgroup handling scans only existing FASTA identifiers');
like($strain,
	qr/my \$cat_write = "\$CATtf\.write\.\$\$".*?print \{\$cat_out\}.*?rename \$cat_write, \$CATtf/s,
	'within-strain category publication streams through an atomic temporary file');
like($strain,
	qr/"flushEvery=i"\s+=> \\\$appendWriteTrigger.*?%outgroupGeneCache = \(\).*?'-flushEvery', \$appendWriteTrigger/s,
	'within-strain extraction exposes its buffer bound to workers and releases per-MGS outgroup caches');
like($strain,
	qr/contextMembersNeeded.*?contextLociNeeded.*?my %keptMemberContext.*?\$MemberContext = \\%keptMemberContext.*?my %keptLocusContext.*?\$LocusContext = \\%keptLocusContext/s,
	'within-strain extraction retains scoring contexts only for potentially ambiguous loci');
unlike($strain,
	qr/normalizeVCFHeaders\.pl/,
	'within-strain consensus regeneration does not invoke VCF normalization');
like($strain,
	qr/sub createAGlist.*?push \@\{\$AGlist\{\$cAssGrp\}\}, \$smpl.*?sub histoMGS/s,
	'within-strain assembly groups retain every sample for consensus extraction');
unlike($strain,
	qr/sub createAGlist.*?CntAimMap.*?sub histoMGS/s,
	'within-strain assembly groups are not collapsed to the last mapping-group sample');
like($strain,
	qr/\@subSds = \@\{\$AGlist\{\$cAssGrp\}\}.*?foreach my \$sd3 \(\@subSds\).*?createConsFastas\(\$cD, \$sd3/s,
	'each assembly-group sample receives sample-specific consensus regeneration');
like($strain,
	qr/Partition whole assembly groups.*?samplesByGroup.*?ownedGroup.*?\$mine\{\$alias\} = 1/s,
	'split extraction assigns complete assembly groups and their catalogue aliases to one worker');
like($strain,
	qr/balance_assembly_groups\(\\%samplesByGroup, \$maxSubJob\).*?\$workerForGroup->\{\$group\} == \$subJob.*?\$plannedSamples \+= scalar\(\@\{\$samplesByGroup\{\$_\}\}\).*?estimated sample load \$workerLoads->\[\$subJob\]\/\$totalWorkerLoad/s,
	'split extraction keeps assembly groups intact while balancing and reporting sample-level work');
like($strain,
	qr/pre-restricted to .*?sample driver\(s\) with target loci/s,
	'split-worker diagnostics distinguish post-index usable sample drivers from assembly groups');
like($strain,
	qr/readGenesSample_Singl\(\s*\$sm, \$writeLink, \$sttime, .*?\$appCnt, \$sampleStatsFH, .*?sampleStatsSeen.*?\$\{\$bufferedSamplesRef\}\+\+.*?appendWriteMGSgenes\(\$writeLink\)/s,
	'expanded assembly-group output is accounted once and flushed by sample to retain the RAM bound');
like($strain,
	qr/if \(\$mySamplesHR\).*?\$unrepresentedWorkerLoci\+\+.*?unless \$maxSubJob/s,
	'split-worker sparsity is summarized instead of reported as missing catalogue data');
like($strain,
	qr/-withinSpecies 1 -relativeNTFraction \$relativeNTFraction .*?-NTfiltPerGene \$GeneLengthMin -GenesPerSpecies \$GenesPerSpecies/s,
	'unfinished trees explicitly pass the named strain coverage filters to buildTree');
unlike($strain, qr/-NTfilt \$relativeNTFraction/,
	'strain workflow does not emit the retired ambiguous NTfilt option');
like($strain,
	qr/my \$GenesPerSpecies = 0\.2;.*?my \$GeneLengthMin = 0\.3;.*?my \$relativeNTFraction = 0\.1;.*?\$placementGenesPerSpecies = 0\.02; \$placementRelativeNTFraction = 0\.01;.*?my \$taxonAwareLocusSelection = 1;.*?"taxonAwareLocusSelection=i" => \\\$taxonAwareLocusSelection.*?-taxonAwareLocusSelection \$taxonAwareLocusSelection/s,
	'strainWithin uses stricter backbone defaults, lower explicit placement thresholds, and taxon-aware selection');
like($strain,
	qr/my \$rateMergePartitions = 1;.*?my \$rateMergeMaxBins = 8;.*?my \$rateMergeTargetSites = 30_000;.*?my \$rateMergeMinLoci = 20;.*?my \$rateMergeMinSites = 20_000;.*?"rateMergePartitions=i" => \\\$rateMergePartitions.*?-rateMergePartitions \$rateMergePartitions.*?-rateMergeMaxBins \$rateMergeMaxBins.*?-rateMergeTargetSites \$rateMergeTargetSites.*?-rateMergeMinLoci \$rateMergeMinLoci.*?-rateMergeMinSites \$rateMergeMinSites/s,
	'strainWithin enables deterministic rate merging and forwards all bin controls');
like($strain,
	qr/maximum_genes_per_sample => 600.*?maximum_tree_loci => 400.*?\$taxonAwareGeneBudget = \$treeLocusBudget < \$presortGenes.*?taxonAwareLocusBudgets\(\$taxonAwareGeneBudget\).*?-taxonAwareMaxLoci \$taxonAwareMaxLoci.*?-taxonAwareCoreLoci \$taxonAwareCoreLoci.*?-taxonAwareCandidateExtra \$taxonAwareCandidateExtra.*?sub taxonAwareLocusBudgets.*?\$maximumLoci \* 0\.8.*?\$maximumLoci \* 0\.3/s,
	'strainWithin scales 80% core, 20% rescue capacity, and 30% QC backfill to its effective gene budget');
like($strain,
	qr/my \$strictBackbone = 1;.*?my \$strictBackboneFraction = 0\.35;.*?my \$strictBackboneMinSamples = 3;.*?my \$placementMinOverlap = 400;.*?"strictBackbone=i"\s+=> \\\$strictBackbone.*?"strictBackboneFraction=f"\s+=> \\\$strictBackboneFraction.*?"strictBackboneMinSamples=i"\s+=> \\\$strictBackboneMinSamples.*?"placementMinOverlap=i"\s+=> \\\$placementMinOverlap/s,
	'strainWithin exposes default-active strict-backbone controls');
like($strain,
	qr/-strictBackbone \$strictBackbone .*?-strictBackboneFraction \$strictBackboneFraction .*?-strictBackboneMinSamples \$strictBackboneMinSamples .*?-placementMinOverlap \$placementMinOverlap/s,
	'strainWithin forwards all backbone controls to buildTree5');
like($strain,
	qr/-tmpSubdir .*?strain_within\/\$MGS.*?-stagedInputDir .*?\$tmpD.*?-completionMarker .*?\$treeStone/s,
	"tree jobs pass lifecycle paths to buildTree5 as ordinary options");
unlike($strain, qr/sub treeInputPrecopyCommand|staged_inputs=\(\)|mapfile -d|ready_test/,
	"strain_within no longer generates Bash input-publication logic");
unlike($strain, qr/\$\{TMPDIR\}\/strain_within|my \$postCmd|touch "?\.shellQuote\(\$treeStone\)/,
	"tree commands contain neither shell TMPDIR expansion nor shell checkpoints");

like($strain,
	qr/sub phase1WorkerCommand.*?Stage-I workers receive extraction\/consensus controls only.*?sub recoverCompletedSplitPhaseI.*?phase1WorkersNeedingRetry.*?Resubmitting invalid Phase-I worker.*?phase1_worker_repair\.queue\.tsv/s,
	'live and resumed Phase I share extraction-only worker commands and durable targeted repair');
like($strain,
	qr/No automatic full-tree resubmission was attempted.*?sub writeTreeFailureAudit.*?failed_missing_output.*?valid_no_tree.*?placement_pending/s,
	'tree outcomes are classified and quarantined without automatic tree resubmission');
like($strain,
	qr/\$totMem = int\(\$totMem \* 2\).*?\$numCoreL = 1.*?-epaOnly 1.*?\$outD2\/treeCmd\.epa_retry\.sh/s,
	'an EPA-only retry gets a one-core doubled-memory job and explicit BuildTree mode');
like($strain,
	qr/sub epaOnlyRetryReady.*?\$onlySubmit.*?IQtree_allsites\.backbone\.treefile.*?strict_backbone\.samples\.tsv.*?status\\tplacement_pending.*?explicit_pending/s,
	'a tree-only resume validates retained placement inputs and recognizes an explicit pending marker');
like($strain,
	qr/sub epaOnlyRetryReady.*?IQtree_allsites\.treefile.*?return '' if -s \$finalTree.*?legacy_missing_final.*?sub prepareEpaOnlyRetryState.*?clear stale completion missing final placed tree.*?create legacy placement-pending marker/s,
	'a legacy retained backbone without the final non-backbone tree is prepared for isolated EPA recovery');
like($strain,
	qr/my \$epaOnlyRetry = exists\(\$MGSepaOnlyRetry\{\$MGS\}\).*?if \(!\$epaOnlyRetry && exists \$MGSnoTree\{\$MGS\}\).*?if \(!\$epaOnlyRetry && exists\(\$ConspecificMGS\{\$MGS\}\)/s,
	'a validated EPA-only retry bypasses later historical no-tree and multicopy filters');
like($strain,
	qr/Placement-only recovery has already paid.*?exists\(\$MGSepaOnlyRetry\{\$specis\[\$b\]\}\).*?\$sizeOfDirs\[\$b\].*?Validated EPA-only recovery queue/s,
	'validated EPA-only retries are sorted ahead of ordinary full-tree retries');
like($strain,
	qr/strain_within\.heartbeat\.tsv.*?strain_within\.failure\.tsv.*?sub writeStrainWorkflowHeartbeat.*?sub writeStrainWorkflowFailure/s,
	'strain workflow persists stage heartbeats and fatal-stage diagnostics');
unlike($strain,
	qr/preflightStrainWorkflow|preflight_executable|preflight_directory|filesystem_capacity/,
	'strain workflow does not preflight environment-wrapped commands as local executables');

like($strain,
	qr/require_complete_linkage => 1.*?Mosaic complete-linkage protection rejected/s,
	'strain extraction requires pairwise confirmation throughout multi-seed Mosaic loci');
like($strain,
	qr/\$workflowHeartbeatPath = File::Spec->catfile\(\$LOGDIR.*?\$SNPconsLOGs = "\$LOGDIR\/SNPconsCalls.*?my \$final = "\$LOGDIR\/\$sampleStatsLogName".*?my \$summary = "\$LOGDIR\/\$sampleStatsSummaryLogName".*?my \$final = "\$LOGDIR\/\$recoveryLogName".*?my \$summary = "\$LOGDIR\/\$summaryLogName"/s,
	'worker heartbeats, SNP logs, sample statistics, recovery accounting, and summaries share LOGandSUB');
like($strain,
	qr/script => \$epaOnlyRetry.*?"\$outD2\/treeCmd\.epa_retry\.sh".*?"\$outD2\/treeCmd\.sh".*?qsubSystem\(\$LOGDIR\."strainAnalysis2\.sh"/s,
	'per-MGS tree scripts retain their compatibility paths while downstream strain-analysis scripts use LOGandSUB');
like($strain,
	qr/sub migrateLegacyOperationalLogs.*?strain_within.*?SNPconsCalls.*?strainSampleStats.*?strainRecovery.*?migrate legacy strain log.*?migrateLegacyOperationalLogs\(\)/s,
	'legacy top-level operational logs are safely migrated into LOGandSUB');
like($strain,
	qr/sub writeSelectionAttritionSummary.*?selection_attrition\.tsv.*?strainSelectionAttrition\.tsv/s,
	'strain summary aggregates completed BuildTree attrition reports');

my $build_tree = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'phylo', 'buildTree5.pl'));
like($build_tree, qr/if \(\$numSeq < 3\)/,
	'three-sample MGS accepted by the wrapper are retained for a minimal tree');

done_testing();
