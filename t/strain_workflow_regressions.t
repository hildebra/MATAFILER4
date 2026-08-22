use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::GenoMetaAss qw(
	readClstrRev readFasta writeClstrRevBinaryShards readClstrRevBinaryShard
	writeSequenceBinaryCache readSequenceBinaryCache
);

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

write_file($cluster_index, join("\n",
	"seed1\tsample1__gene1,sample2__gene2,alias1__gene3",
	"seed2\tsample2__gene4,sample3__gene5",
	"drop\tsample1__ignored",
).'\n');
my @binary_shards = map { File::Spec->catfile($tmp, "cluster.worker.$_.bin") }
	0 .. 2;
my $shard_fingerprint = 'a' x 64;
my $shard_metadata = writeClstrRevBinaryShards(
	$cluster_index,
	{ seed1 => 1, seed2 => 1 },
	{ sample1 => 0, alias1 => 0, sample2 => 1, sample3 => 1 },
	\@binary_shards,
	$shard_fingerprint,
);
is_deeply(
	[map { $_->{records} } @{$shard_metadata}],
	[1, 2, 0],
	'binary cluster-index publication records only selected clusters in each worker partition',
);
is_deeply(
	readClstrRevBinaryShard($binary_shards[0], $shard_fingerprint, 0, 3),
	{ seed1 => 'sample1__gene1,alias1__gene3' },
	'binary cluster-index shard preserves member order and catalogue aliases',
);
is_deeply(
	readClstrRevBinaryShard($binary_shards[1], $shard_fingerprint, 1, 3),
	{ seed1 => 'sample2__gene2', seed2 => 'sample2__gene4,sample3__gene5' },
	'binary cluster-index shard contains only the assigned worker members',
);
is_deeply(
	readClstrRevBinaryShard($binary_shards[2], $shard_fingerprint, 2, 3),
	{},
	'an empty worker receives a valid binary shard rather than falling back to the full index',
);
my $wrong_shard_provenance = eval {
	readClstrRevBinaryShard($binary_shards[0], 'b' x 64, 0, 3);
	1;
};
ok(!$wrong_shard_provenance && $@ =~ /Invalid binary cluster-index shard header/,
	'binary cluster-index reader rejects a shard from another provenance generation');
my $corrupt_shard = File::Spec->catfile($tmp, 'cluster.worker.corrupt.bin');
my $corrupt_contents = slurp($binary_shards[0]);
substr($corrupt_contents, -1, 1) = chr(ord(substr($corrupt_contents, -1, 1)) ^ 1);
write_file($corrupt_shard, $corrupt_contents);
my $corrupt_shard_accepted = eval {
	readClstrRevBinaryShard($corrupt_shard, $shard_fingerprint, 0, 3);
	1;
};
ok(!$corrupt_shard_accepted && $@ =~ /payload digest mismatch/,
	'binary cluster-index reader rejects payload or trailer corruption');

my $protein_cache = File::Spec->catfile($tmp, 'catalog.proteins.bin');
my $protein_fingerprint = 'c' x 64;
my $protein_metadata = writeSequenceBinaryCache(
	{ seed1 => 'MPEPTIDE', seed2 => 'MSECOND', empty => '' },
	$protein_cache,
	$protein_fingerprint,
);
is_deeply(
	$protein_metadata,
	{ records => 2, bytes => -s $protein_cache },
	'binary sequence cache records only nonempty reference proteins and reports its durable size',
);
is_deeply(
	readSequenceBinaryCache($protein_cache, $protein_fingerprint),
	{ seed1 => 'MPEPTIDE', seed2 => 'MSECOND' },
	'binary sequence cache round-trips the common catalogue-protein subset exactly',
);
my $wrong_protein_provenance = eval {
	readSequenceBinaryCache($protein_cache, 'd' x 64);
	1;
};
ok(!$wrong_protein_provenance && $@ =~ /Invalid binary sequence-cache header/,
	'binary sequence cache rejects a different selected-catalogue generation');
my $corrupt_protein_cache = File::Spec->catfile($tmp, 'catalog.proteins.corrupt.bin');
my $corrupt_protein_contents = slurp($protein_cache);
substr($corrupt_protein_contents, -1, 1) =
	chr(ord(substr($corrupt_protein_contents, -1, 1)) ^ 1);
write_file($corrupt_protein_cache, $corrupt_protein_contents);
my $corrupt_protein_accepted = eval {
	readSequenceBinaryCache($corrupt_protein_cache, $protein_fingerprint);
	1;
};
ok(!$corrupt_protein_accepted && $@ =~ /payload digest mismatch/,
	'binary sequence cache rejects same-size protein-cache corruption');

my $strain = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'MGS', 'strain_within.pl'));
my $strain2 = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'MGS', 'strain_within_2.2.pl'));
my $mgs = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'MGS.pl'));
my $build_tree = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'phylo', 'buildTree5.pl'));
my $internal_config = slurp(File::Spec->catfile($Bin, '..', 'Mods', 'config_internal.txt'));
my $site_config_template = slurp(File::Spec->catfile($Bin, '..', 'Mods', 'config.old'));
my $neighbor_tree_r = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'R_scripts', 'neighborTree.R'));

my ($durable_output_helper_source) = $strain =~
	/(sub strainOutputHasDurablePhaseIState \{.*?\n\})\n\nsub phase1PathStatComponent/s;
ok(defined($durable_output_helper_source),
	'existing-output evidence helper is available for isolated testing');
$durable_output_helper_source =~
	s/\Asub strainOutputHasDurablePhaseIState/sub/;
my $durable_output_helper = eval $durable_output_helper_source;
die "Cannot compile existing-output evidence helper: $@" if $@;

my $fresh_strain_output = File::Spec->catdir($tmp, 'fresh_strain_output');
mkdir $fresh_strain_output or die "Cannot create $fresh_strain_output: $!";
for my $operational (qw(LOGandSUB stones strainsScr1 .scratch)) {
	my $directory = File::Spec->catdir($fresh_strain_output, $operational);
	mkdir $directory or die "Cannot create $directory: $!";
}
ok(!$durable_output_helper->(
		$fresh_strain_output,
		File::Spec->catfile($fresh_strain_output, 'LOGandSUB', 'missing.summary'),
	),
	'a fresh output containing only operational directories permits a subset build');
my $existing_mgs_directory = File::Spec->catdir($fresh_strain_output, 'MGS.1');
mkdir $existing_mgs_directory
	or die "Cannot create $existing_mgs_directory: $!";
ok($durable_output_helper->($fresh_strain_output),
	'an existing per-MGS directory blocks a destructive subset rebuild');
my $summary_evidence = File::Spec->catfile($tmp, 'strainSampleStats.summary.tsv');
write_file($summary_evidence, "durable\n");
ok($durable_output_helper->(File::Spec->catdir($tmp, 'absent_output'), $summary_evidence),
	'a durable Phase-I summary blocks a destructive subset rebuild without a directory scan');

like($strain, qr/sub consensusInputState .*?\$nt_ready && \$aa_ready.*?return 'regenerate' if \$vcf_ready/s,
	'consensus resume requires the paired NT and AA outputs and repairs from VCF');
like($strain,
	qr/my \$onlyMSA = 0;.*?"onlyMSA=i"\s+=> \\\$onlyMSA.*?-onlyMSA must be 0 or 1.*?-onlyMSA 1 cannot be combined with -placeOnBackbone 1/s,
	'strainWithin exposes and validates alignment-only mode');
like($strain,
	qr/msaOnly\.complete\.tsv.*?MSA\/MSAli\.fna\.gz.*?-MSAprogram \$MSAprog -onlyMSA \$onlyMSA/s,
	'strainWithin forwards MSA-only mode and audits its durable alignment outcome');
like($strain,
	qr/if \(\$onlyMSA\) \{.*?strain_within_2\.2\.pl were not launched.*?exit\(0\);.*?my \$MGSabundance/s,
	'MSA-only controller runs stop before tree-dependent postprocessing');
like($build_tree,
	qr/"onlyMSA=i" => \\\$onlyMSA.*?if \(\$onlyMSA\) \{.*?finalizeMSAArtifacts.*?msa_complete.*?exit\(0\);.*?my \$trRetH/s,
	'BuildTree finalizes and records the concatenated MSA before tree inference');
like($strain, qr/readFasta\(\$fastaf,1,"\\\\s",\\%subG\).*?readFasta\(\$fastafAA,0,"\\\\s",\\%subG\)/s,
	'within-strain extraction reads only candidate consensus genes');
like($strain, qr/-completionMarker "?\.shellQuote\(\$treeStone\)/,
	"tree jobs delegate validated completion markers to buildTree5");
unlike($strain, qr/test -s "?\.shellQuote\(\$IQtreef\).*?touch "?\.shellQuote\(\$treeStone\)/s,
	"tree jobs no longer encode completion validation in shell");
like($strain, qr/my \@staleOutputs = \$record->\{epa_only\}.*?qw\(stone terminal\).*?\$record->\{msa_only\}.*?qw\(stone terminal placement_pending\).*?qw\(stone tree terminal placement_pending\).*?retry_unlink/s,
	'EPA-only and MSA-only submissions retain their reusable primary outputs while ordinary retries clear stale trees');
like($strain, qr/clear_split_generation\(\$splitManifest.*?write_split_generation\(\$splitManifest.*?printf '%s\\\\n'/s,
	'a new split generation clears stale state and tags each worker completion');
like($strain, qr/ConspecificMGS\.\$subJob\.log.*?sub mergeConspecificLogs/s,
	'split workers write isolated conspecific logs that are explicitly merged');
like($strain, qr/\$onlySubmit == 0 && !\$subJob/,
	'split children cannot recursively clean shared MGS output directories');
ok(index($strain, 'my $publishedInputsReady = !$epaOnlyRetry') >= 0
	&& index($strain, 'if (!$epaOnlyRetry && !($publishedInputsReady && !$mustRegenerateInputs))') >= 0
	&& index($strain, '$scratchInputsReady ||= prepareMGSInputSet($MGS,$tmpD);') >= 0,
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
	qr/loadRecoveryContributionIndex\(\)\s+unless \$recoveryContributionIndexReady \|\| \$leanOnlySubmitResume.*?sub writeRecoveryContributionIndex.*?retry_rename\(\$temporary, \$path.*?sub loadRecoveryContributionIndex.*?\$recoveryContributionIndexReady = 1.*?sub mergeRecoveryLogs.*?writeRecoveryContributionIndex\(\);.*?retry_rename\(\$temporary, \$final.*?retry_unlink\(\$_.*?for \@parts/s,
	'worker provenance is persisted before disposable recovery logs are removed and is reloadable after restart');
like($strain,
	qr/\$aggregateComplete &&= -s \$mergeCheckpoint.*?my \@expectedWorkers.*?if \(\@validationErrors\).*?return \$aggregateComplete.*?retry_unlink\(\$mergeCheckpoint.*?retry_rename\(\$checkpointTemporary, \$mergeCheckpoint.*?retry_unlink\(\$part/s,
	'a last-written checkpoint protects committed aggregates and worker parts survive until validation and commit');
like($strain, qr/qsubSystemJobAlive\([^\n]+QSBoptHR[^\n]+if [^\n]+doSubmit/,
	'dry runs do not poll scheduler jobs that were never submitted');
like($strain, qr/\$nxtCmd \.= "-submit \$doSubmit ";.*?-qsubSystem/s,
	'postprocessing inherits submission state and the selected queue backend');
like($strain, qr/\$nxtCmd \.= "-MGSphylo "\.shellQuote\(\$treeFile\).*?if \$treeFile ne ""/,
	'postprocessing receives the source MGS tree for outgroup recovery');
like($strain,
	qr{loadTreeOutgroupCandidates\(\$targetMGS\).*?sub loadTreeOutgroupCandidates .*?" --all".*?open my \$bulk.*?\$TreeOutgroupCandidatesBulkLoaded = 1.*?sub treeOutgroupCandidates .*?A failed bulk call.*?--preferred-tip}s,
	'Phase II imports all source-tree neighbour candidates in one call and retains a Mosaic-aware individual fallback');
like($strain,
	qr{tempfile\(.*?strain_mosaic_outgroups.*?TMPDIR => 1, UNLINK => 1.*?print \{\$preferredFh\} "\$MGS\\t\$PreferredOutgroup\{\$MGS\}\\n".*?" --preferred "}s,
	'the bulk R call receives all Mosaic preferences through one automatically removed temporary file');
like($strain,
	qr{my \(\$MGS, \$decision, \$preferred, \$preferredDistance, \$cutoff, \$candidateText\).*?split /\\t/, \$line, 6.*?Mosaic decisions:}s,
	'Perl imports the authoritative R ordering and summarizes Mosaic plausibility decisions');
like($strain,
	qr{if \(length\(\$treeFile\)\).*?Do not reinsert a rejected preference.*?push \@candidates, treeOutgroupCandidates\(\$MGS\).*?elsif \(exists\(\$PreferredOutgroup}s,
	'the reference preload follows R ordering when a phylogeny is present and uses Mosaic directly only without one');
like($strain,
	qr{if \(\$treeFile ne ""\).*?This order is authoritative.*?push \@candidates, treeOutgroupCandidates\(\$MGS\).*?elsif \(exists\(\$PreferredOutgroup}s,
	'the final outgroup chooser cannot reinsert a Mosaic proposal rejected by R');
like($neighbor_tree_r,
	qr{identical\(target, "--all"\).*?ape::cophenetic\.phylo\(tree\).*?for \(tip in tree\$tip\.label\).*?ranked\$decision.*?paste\(ranked\$candidates}s,
	'neighborTree bulk mode computes distances once and emits authoritative decision and candidate columns per tree tip');
like($neighbor_tree_r,
	qr{--preferred.*?--preferred-tip.*?ranked_neighbors <- function.*?stats::quantile.*?nearestDistance \* preferredNearestFactor.*?preferredDistance > cutoff.*?candidateNames\[candidateNames != preferred\].*?c\(preferred, candidateNames}s,
	'a plausible Mosaic outgroup is promoted while an extreme-distance proposal is excluded from the R result');
like($strain2, qr/"MGSphylo=s"\s*=>\s*\\\$MGSphylo.*?sub resolveOutgroup .*?data\.log.*?treeCmd\.sh.*?MGSphylo/s,
	'postprocessing preserves logged or saved outgroups and falls back to the source MGS tree');
like($strain, qr/sub assertSafeWorkflowRemoval .*?resolved_default.*?Refusing to remove unowned custom output directory/s,
	'custom recursive output removal requires a workflow-owned directory');
like($strain2,
	qr/\$waitForAnalysis->\('strainStats'\);.*?my \$shouldCombineStrainStats = \$forceStrainStats \|\| \$strainTaskCount > 0 \|\| !-s \$RsummaryTab;.*?if \(\$shouldCombineStrainStats\) \{.*?combineResults\(0\);.*?\} else \{.*?Reusing existing combined strainStats overview/s,
	'strainStats stores are combined only when missing or newly stale before the population phase is awaited');
like($strain2,
	qr/\$combineResultsR --path .*?shellQuote\(\$FMGpD\).*?--outDir .*?shellQuote\(\$FMGpD\)/s,
	'combineResults uses only supported path and output options; it auto-detects population-genetics stores');
unlike($strain2, qr/--include-popgen/,
	'combineResults is not passed the obsolete population-genetics inclusion option');
like($strain2,
	qr/if \(\$includePopGen\).*?did not produce the population overview table/s,
	'combineResults validates population outputs only for the population-inclusive pass');
unlike($strain2, qr/open my \$summary_fh.*?\$TXTreport/s,
	'postprocessing no longer concatenates per-MGS text reports into its overview');
unlike($strain2, qr/test -s "\.shellQuote\(\$analysisReport\)/,
	'the legacy text report is not required for RDS-based completion');
like($strain2,
	qr/"popGenStats=i".*?"popGenStrictOutgroup=i".*?"popGenGeneticCode=i".*?"popGenCodonStart=i".*?"popGenSeed=i".*?"popGenLegacyTextOutput=i".*?my \$popGenStatsR = \$doPopGenStats \? getProgPaths\("pogenStats"\).*?popGenStats\.output\.Rds.*?\$popGenStatsR .*?\$destBaseD, \$refMap, \$destD.*?--subsample .*?\$popGenSubsample.*?--ncore \$jobCores.*?--genetic-code \$popGenGeneticCode.*?--codon-start \$popGenCodonStart.*?--seed \$popGenSeed.*?--individual-column .*?\$individualVar.*?--outgroup .*?\$OG.*?--strict-outgroup.*?--legacy-text-output.*?\$strainFile = "\$destD\/IQtree_allsites\.strains\.txt".*?PopGenStats owns validation and fallback behavior.*?\$popGenCommand --strain-file .*?\$strainFile.*?test -s .*?\$popGenStore/s,
	'population genetics forwards parallel, outgroup, codon, seed, identity, and the strain-file path to its durable RDS workflow');
like($strain2,
	qr/\$popGenStatsReady = !\$doPopGenStats \|\| -s \$popGenStore.*?if \(\$doPopGenStats\) \{.*?\$waitForAnalysis->\('popGenStats'\);.*?my \$shouldCombinePopGenStats = \$forcePopGenStats \|\| \$popGenTaskCount > 0 \|\| !-s \$popGenSummaryTab;.*?if \(\$shouldCombinePopGenStats\) \{.*?combineResults\(1\);.*?\} else \{.*?Reusing existing combined PopGenStats overview.*?combineResults\.R did not produce the population overview table \$popGenSummaryTab/s,
	'existing population RDS stores and aggregate tables are reused, while missing or newly stale tables are combined after the population phase');
like($strain2,
	qr/\$popGenSubsampleSummaryTab = "\$FMGpD\/popGenStats\.subsamples\.tsv".*?Combined subsampled population-genetics overview/s,
	'subsampled population-genetics output is surfaced separately from the full population table');
unlike($strain2, qr/if \(0\)\{#rerun popgen stats\?\?/,
	'population genetics is no longer hidden behind a disabled legacy block');
like($strain2,
	qr/sub newickSampleCount \{.*?\$commas = \(\) = \$newick =~ \/,\/g;.*?\$tips = \$commas \+ 1;.*?return \$tips > 0 \? \$tips : 1;/s,
	'R-job batch cost uses an exact Newick tip/sample count rather than file size');
like($strain2,
	qr/my \@k2d = sort \{ \$treeSamples\{\$b\} <=> \$treeSamples\{\$a\}.*?my \$largestMGS = \$k2d\[0\];.*?\$batchSampleBudget = \$treeSamples\{\$largestMGS\};.*?for my \$analysisKind \(qw\(strainStats popGenStats\)\).*?\$batchSampleCost = \$treeSampleCount < \$batchSampleBudget.*?3 \* \$treeSampleCount.*?\$curBatchSamples \+ \$batchSampleCost > \$batchSampleBudget.*?\$curBatchSamples \+= \$batchSampleCost.*?\$curBatchSamples >= \$batchSampleBudget/s,
	'largest MGS defines the sample budget while each smaller MGS counts three times for batching overhead');
like($strain2,
	qr/my \$submitRAnalysisBatch = sub \{.*?qsubSystem\(.*?\$script, \$batchCmd, \$batchCores, \$batchMemory, \$batchLabel/s,
	'R-analysis batch submissions retain the standard core count and memory request through one reusable helper');
like($strain2,
	qr/for my \$analysisKind \(qw\(strainStats popGenStats\)\).*?\$jobCores = \$nCore;.*?\$submitRAnalysisBatch->\(.*?\$analysisKind, \$batchDestD\."\$analysisKind\.Ranalysis\.sh".*?\$batchCmd, \$batchCores, \$batchMemory, \$batchLabel, \$batchStores/s,
	'separate strainStats and PopGenStats batches use the requested standard core count');
like($strain2,
	qr/getProgPaths\("maxMF4mem", 0\).*?my \$rAnalysisOOMMaxMemoryGB = 512;.*?\$rAnalysisOOMMaxMemoryGB = \$1 \+ 0;.*?my \$rAnalysisMemoryBaseGB = 24;.*?my \$rAnalysisMemoryCoreThreshold = 4;.*?my \$rAnalysisRetryRounds = 2;.*?sub rAnalysisMemoryForCores \{.*?\$extraCores = \$cores - \$rAnalysisMemoryCoreThreshold;.*?return \$memoryGB\."G";/s,
	'R-analysis reads maxMF4mem, grows for extra cores, and has two OOM-aware retry rounds');
like($strain2,
	qr/if \(\$analysisKind eq 'strainStats'.*?\$cmd \.= "if ! test -s ".*?rm -f .*?\$analysisStore.*?\$cmd \.= "fi\\n";.*?\$batchStores->\{\$d\} = \$analysisStore.*?if \(\$analysisKind eq 'popGenStats'.*?\$cmd \.= "if ! test -s ".*?rm -f .*?\$popGenStore.*?\$cmd \.= "fi\\n";.*?\$batchStores->\{\$d\} = \$popGenStore/s,
	'replayed R-analysis batch scripts skip nonempty RDS outputs and retain the exact missing-output inventory');
like($strain2,
	qr/use Mods::SlurmAccounting qw\(slurm_tree_memory_summary next_oom_retry_memory_mb\);.*?my \$batchHasMissingStores = sub \{.*?my \$rAnalysisOOMBatches = sub \{.*?slurm_tree_memory_summary.*?while \(1\).*?last if \$round >= \$rAnalysisRetryRounds;.*?Retrying \$analysisKind round \$round\/\$rAnalysisRetryRounds.*?increaseRAnalysisMemory\(\$batch->\{memory\}\).*?\.retry\$round\.sh.*?sub increaseRAnalysisMemory \{.*?next_oom_retry_memory_mb.*?sub reportRAnalysisSchedulerFailures \{.*?warn "Slurm reported failed .*?incomplete result batches will be retried/s,
	'incomplete R-analysis batches are retried twice, using Mods Slurm accounting to double only OOM memory requests');
unlike($strain2, qr/die "Slurm reported failed \$analysisKind job\(s\)/,
	'R-analysis failures do not abort before their bounded retry rounds complete');
like($strain2,
qr/"redoStrainStats=i".*?"redoPopGenStats=i".*?my \$forceStrainStats = \$rewriteRanalysis \|\| \$redoStrainStats;.*?my \$forcePopGenStats = \$rewriteRanalysis \|\| \$redoPopGenStats;.*?\$analysisKind eq \x27strainStats\x27 && \$forceStrainStats.*?\$analysisKind eq \x27popGenStats\x27 && \$forcePopGenStats.*?rm -f .*?\$analysisStore.*?rm -f .*?\$popGenStore/s,
	'targeted statistic flags force only their selected result stores and combined overview');
my @analysisPhaseMarkers = (
	q{$waitForAnalysis->('strainStats');},
	q{combineResults(0);},
	q{my ($networkDep, $networkStone) = strainNetwork();},
	q{my ($treeWasDep, $treeWasStone) = treeWas();},
	q{my ($phyloFigureDep, $phyloFigureStone) = visualizeSignPhylos();},
	q{$waitForAnalysis->('popGenStats');},
	q{combineResults(1);},
);
my $previousPhaseMarker = -1;
my $analysisPhasesInOrder = 1;
for my $phaseMarker (@analysisPhaseMarkers) {
	my $phaseMarkerPosition = index($strain2, $phaseMarker);
	$analysisPhasesInOrder = 0
		if $phaseMarkerPosition < 0 || $phaseMarkerPosition <= $previousPhaseMarker;
	$previousPhaseMarker = $phaseMarkerPosition;
}
ok($analysisPhasesInOrder,
	'strain summaries and submitted network/treeWAS/phylogeny work start before the independent population phase is awaited');
like($strain2,
	qr/sub visualizeSignPhylos\{.*?\$phyloFigureStone = "\$FMGpD\/phyloFigures\.sto";.*?\$cmdPic \.= "touch ".*?if \(!-e \$phyloFigureStone\).*?qsubSystem\(\s*"\$FMGpD\/phyloFigures\.sh", \$cmdPic, 1, "24G", "phyloFigures".*?return \(\$dep, \$phyloFigureStone\);/s,
	'significant-phylogeny plotting is submitted with a durable checkpoint and adequate memory');
like($strain2,
	qr/\$phyloFigureCheckpoint = "\$FMGpD\/phyloFigures\.sto";.*?unlink \$phyloFigureCheckpoint.*?if -e \$phyloFigureCheckpoint/s,
	'explicit rewrites invalidate the submitted phylogeny-figure checkpoint');

like($strain2,
	qr/my \$cmdPrelude = "set -euo pipefail\\nulimit -s 20000\\n";/s,
	'generated R-analysis scripts explicitly preserve non-zero R exits, including Slurm OOM termination');
unlike($strain2, qr/MATAFILER_R_ANALYSIS_CORES/,
	'R and PopGenStats commands use the resolved allocated core count rather than an unexpanded shell variable');
unlike($strain2, qr/\$batchSize/,
	'R-job submission no longer uses a fixed phylogeny count per batch');
unlike($strain2, qr/if \(\$doSubmit && -d \$destD\)/,
	'partial result recovery does not erase an entire within directory outside an explicit rewrite');
like($strain2,
	qr/my \$networkDir = "\$FMGpD\/networks";.*?remove_tree\(\$networkDir\) if -d \$networkDir/s,
	'explicit rewrites clear the workflow-owned network cache');
like($strain2,
	qr/\$networkGraph = "\$netDir\/strain_graph\.Rds".*?\$networkStone && !-s \$networkGraph.*?Ignoring incomplete network checkpoint.*?unlink \$networkStone.*?test -s .*?\$networkGraph.*?touch .*?\$networkStone/s,
	'network completion requires a nonempty graph result as well as its checkpoint');
like($internal_config,
	qr/^combineResults_R\t\[Rscript\] \[MGSTKDir\]\/combineResults\.R\tenv:MGSTK$/m,
	'the combineResults command is configured through the MG-STK R environment');
like($strain, qr/sub limitedWarn .*?warningExampleLimit.*?Further '\$category' warnings are suppressed/s,
	'repetitive strain warnings retain examples and announce suppression');
like($strain, qr/Suppressed warning summary:.*?sort grep/s,
	'suppressed strain warnings receive a categorized exit summary');
unlike($strain, qr/print "\$cD\\n"/,
	'strain extraction no longer prints a raw working-directory path for every sample');
like($strain, qr/my \$version = 1\.36;/,
	'workflow behavior changes retain an explicit version marker');
like($strain,
	qr/my \$SNPcaller = "MPI";.*?"SNPcaller=s"\s*=> .*?SNPcaller.*?-SNPcaller must be MPI or FB.*?genes\.shrtHD\.SNPc\.\$\{SNPcaller\}\.fna\.gz.*?allSNP\.\$\{SNPcaller\}\.vcf\.gz/s,
	'strain Phase I selects MPI or FB caller-specific consensus and VCF filenames');
like($strain,
	qr/my \$legacyMPIContract = \$phase1ContractState eq 'missing' && \$SNPcaller eq 'MPI'.*?if \(\$onlySubmit && !\$subJob.*?'building_match'.*?Cannot reuse Phase-I outputs.*?if \(\$onlySubmit && \$subJob.*?'building_match'.*?!\$legacyMPIContract.*?Split worker caller contract/s,
	'caller provenance resumes compatible builds and rejects missing FB or incompatible state');
like($strain,
	qr/my \$phase1InputContractVersion = 3;.*?sub phase1PathStatComponent.*?different device number.*?\$contractVersion <= 2.*?\@metadata\[0, 1, 7, 9\].*?\@metadata\[1, 7, 9\]/s,
	'Phase-I provenance excludes node-local device identity from current cross-node contracts');
like($strain,
	qr/sub phase1GuideStatFingerprint.*?strain-phase1-guide-stat-v2.*?strain-phase1-guide-stat-v3.*?\$canonical\.srt.*?\$canonical\.srt\.gene2MGS.*?\$observation/s,
	'Phase-I provenance includes the original, sorted, indexed, and observation guide inputs');
like($strain,
	qr/sub phase1CatalogStatFingerprint.*?subset\.cats.*?compl\.incompl\.\$identity\.fna.*?fna\.clstr\.idx.*?prot\.faa.*?eggNOGmapper_NOG\.geneAss.*?split\(\/,\/, \$mapSpec/s,
	'Phase-I provenance includes each direct catalog, marker, annotation, and map input');
like($strain,
	qr/catalog_inputs_fingerprint.*?my \@values = \(\$phase1InputContractVersion, \$status, \$phase1CatalogIdentity,\s*\$phase1CatalogInputFingerprint, \$phase1MGSGuideFingerprint,\s*\$useGTDBmg, \$clusterID/s,
	'Phase-I contract persists the catalog and MGS-guide fingerprints');
like($strain,
	qr/\$recordedPhase1ContractVersion < \$phase1InputContractVersion.*?Cannot upgrade legacy Phase-I input contract.*?persistPhase1InputContract/s,
	'a resumed parent upgrades compatible legacy provenance before dispatching repair workers');
like($strain,
	qr/phase1CatalogStatFingerprint\(.*?\$recordedVersion\).*?phase1GuideStatFingerprint\(\$MGSfileOri, \$recordedVersion\)/s,
	'legacy provenance recalculation uses the preserved original guide rather than its sorted working path');
like($strain,
	qr/if \(!\$subJob && \$doPopGenStats && \$rmMSA\).*?Population genetics requires per-locus nucleotide MSAs/s,
	'population-genetics MSA warnings and overrides are restricted to the parent process');
ok(index($strain, q{persistPhase1InputContract($activePhase1InputContract, 'building')}) >= 0
	&& index($strain, q{persistPhase1InputContract(File::Spec->catfile($LOGDIR, $phase1InputContractName))}) >= 0
	&& index($strain, q{atomic_write_text($path, phase1InputContractContents($status)}) >= 0,
	'Phase-I caller provenance is atomically marked building and completed only after handoff');
like($strain,
	qr/-redo all cannot be combined with -MGSsubset.*?if \$redoMode eq 'all' && length\(\$subsMGSstr\).*?Tree for outgroup specified, but file is missing or empty:.*?!-s \$treeFile/s,
	'unsafe subset full-redo and empty outgroup-tree inputs fail before downstream work');
like($strain,
	qr/my \$unsafeSubsetRebuild = !\$onlySubmit && !\$subJob && length\(\$subsMGSstr\).*?strainOutputHasDurablePhaseIState.*?if \(\$unsafeSubsetRebuild\).*?shared non-subset results.*?would be cleared/s,
	'existing output state refuses ordinary destructive subset rebuilds while fresh roots remain eligible');
like($strain,
	qr/my \$rmMSA = 1;.*?my \$doPopGenStats = 1;.*?my \$popGenStrictOutgroup = 0;.*?my \$popGenGeneticCode = 1;.*?my \$popGenCodonStart = 1;.*?my \$popGenSeed = 1;.*?"popGenStats=i"\s*=> \\\$doPopGenStats.*?"popGenStrictOutgroup=i".*?"individualVar=s".*?if \(!\$subJob && \$doPopGenStats && \$rmMSA\).*?\$rmMSA = 0;.*?-rmMSA \$rmMSA.*?-popGenStats \$doPopGenStats.*?-popGenStrictOutgroup \$popGenStrictOutgroup.*?-popGenGeneticCode \$popGenGeneticCode.*?-popGenCodonStart \$popGenCodonStart.*?-popGenSeed \$popGenSeed.*?-popGenLegacyTextOutput \$popGenLegacyTextOutput/s,
	'population genetics retains MSAs and forwards reproducible configuration through strainwithin2');
like($strain, qr/Retain the Phase-I locus map.*?second catalogue-wide gene2tax scan.*?\$SIgenes and \$COGprios are reused/s,
	'Phase II reuses the Phase-I selected gene map rather than clearing and rebuilding it');
like($strain,
	qr/Preparing core-first exact outgroup-reference demands.*?my %broadCOG = map.*?\$cogTaxa\{\$_\} >= \$broadMinimumTaxa.*?exists\(\$preferredCoreGeneSet->\{\$gene\}\).*?readFasta\(\$refFAA, 1, "\\\\s", \\%requiredAA,.*?readFasta\(\$refFNA, 1, "\\\\s", \\%requiredNT,/s,
	'outgroup references use a core-first, broad-fallback demand manifest and stream only exact requested FNA/FAA records');
like($strain,
	qr/my \$outgroupCoreMinLoci = 0;.*?"outgroupCoreMinLoci=i".*?\$outgroupCoreMinLoci = int\(\$treeLocusBudget \* 0\.20 \+ 0\.999999\).*?if \$outgroupCoreMinLoci == 0;.*?\$minimumOutgroupLoci = \$outgroupDemandMinimum\{\$MGS\} \/\/ \$MGStoolowGsThr/s,
	'the outgroup floor defaults to 20% of the final-tree locus budget and is enforced per MGS');
like($strain,
	qr/my \$outgroupReferenceGeneCap = 2500;.*?\(\$candidateSIgenes, \$candidateGene2COG.*?readGene2tax\(.*?\$outgroupReferenceGeneCap.*?my \$addCandidate = sub.*?return if \$retainedForMGS >= \$outgroupReferenceGeneCap/s,
	'candidate reference maps use a generous 2,500-gene-per-outgroup-MGS cap while prioritizing the acceptance demand');
like($strain,
	qr/&& exists\(\$PreferredOutgroupGene.*?&& \(\$broadCOG/s,
	'an exact Mosaic link is usable only for a broadly available or preferred-core locus');
unlike($strain, qr/prepareSelectiveOutgroupReferenceCache|outgroupReferenceCacheActive|outgroup_reference_cache/,
	'the Phase II selective outgroup cache and its index lifecycle are absent');
like($strain, qr/sub outgroupRequirementLoci.*?preparedOutgroupLog.*?CATstdof\.tmp/s,
	'only raw staged inputs without a finalized outgroup overlay contribute outgroup requirements');
like($strain, qr/sub outgroupRequirementLoci.*?selected_gene_map.*?outgroup requirement category/s,
	'resume uses already-selected loci before considering a raw category scan');
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
	qr/my %sampleStatsSeen;.*?my \$nextSampleProgress = \$extractionStarted \+ 60;.*?readGenesSample_Singl\(.*?stepProgress\("consensus-gene extraction"/s,
	'STDOUT is redirected once around the complete sample loop while the duplicated handle carries only TSV records');
like($strain,
	qr/sub writeSampleStats \{.*?without a sample name.*?duplicate row.*?Refusing to emit an empty.*?for my \$target \(\$fh, \$sampleStatsPartFH\).*?print \{\$target\} \$row, "\\n"/s,
	'every sample-statistics record has a sample name, is nonempty, and is emitted to stdout and its worker table at most once');
like($strain,
	qr/sub mergeSampleStats .*?Wrong sample-statistics field count.*?Duplicate sample-statistics row.*?aggregate_sample_rows.*?STAGE I SAMPLE SUMMARY \(all workers\)/s,
	'all worker tables are validated, aggregated, saved, and reported at the end of Step 1');
like($strain,
	qr/STAGE I SAMPLE SUMMARY \(all workers\).*?join\("; ", \@summaryPairs\).*?loci_histogram_rows.*?Used MGS retained-loci histogram/s,
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
	qr/locus-model construction.*?catalogue_drivers=.*?resolved_loci=.*?consensus-gene extraction and publication.*?full-tree input sizing/s,
	'major extraction and tree-preparation stages also report concise completion statistics');
like($strain,
	qr/locus-model cluster-index scan.*?represented_seed_clusters=.*?locus-model catalogue-protein loading.*?loaded_proteins=.*?locus-group construction.*?ranked_clusters=.*?resolved_loci=.*?worker-member materialization.*?locus_sample_combinations=/s,
	'locus-model construction reports scan, protein, grouping, and worker-projection subphase timings');
like($strain,
	qr/sub phase1SelectedGeneFingerprint.*?phase1PathStatComponent\(\$gene2taxF\).*?\$presortGenes.*?for my \$mgs \(\@specis\).*?sub phase1IndexShardFingerprint.*?phase1PathStatComponent\(\$clusterIndex\).*?phase1SelectedGeneFingerprint.*?sort keys %\{\$workerForSample\}.*?sub publishPhase1IndexShards.*?writeClstrRevBinaryShards.*?retry_rename\(\$temporary\[\$worker\], \$final\[\$worker\].*?atomic_write_text\(File::Spec->catfile\(\$base, 'manifest\.tsv'\).*?sub loadPhase1ClusterIndex/s,
	'parent publishes provenance-bound binary worker shards atomically before the cache manifest');
like($strain,
	qr/phase1IndexShardCacheState.*?readClstrRevBinaryShard.*?binary_worker_shard.*?readClstrRev\(\$clusterIndex, 0, \$Gene2COG, \$mySamples\).*?full_index_fallback/s,
	'split workers validate their binary shard and retain the original full-index fallback');
like($strain,
	qr/sub phase1ProteinCacheFingerprint.*?phase1PathStatComponent\(\$proteinFile\).*?phase1SelectedGeneFingerprint.*?sub publishPhase1ProteinCache.*?writeSequenceBinaryCache.*?retry_rename\(\$temporary, \$final.*?atomic_write_text\(File::Spec->catfile\(\$base, 'manifest\.tsv'\).*?sub loadPhase1CatalogProteins.*?unless \(\$maxSubJob > 1\).*?catalogue_fasta.*?readSequenceBinaryCache.*?binary_common_subset.*?readFasta\(\$proteinFile, 1, "\\\\s", \$Gene2COG, \{ fai => 1 \}\).*?catalogue_fasta/s,
	'parent publishes one provenance-bound binary protein subset and workers retain FASTA fallback');
like($strain,
	qr/historical exclusion loading.*?excluded_MGS=.*?outgroup-reference preparation.*?reference_NT=.*?MGS_with_outgroup_candidates=/s,
	'historical exclusions and outgroup-reference preparation report their final counts');
like($strain,
	qr/"redoEPAfilter:i" => sub.*?if \(\$redoEPAfilter\).*?epa_result\.jplace.*?IQtree_allsites\.treefile.*?retry_unlink\(\$placedTree.*?treeDone\.sto.*?retry_unlink\(\$completion.*?placementPending\.sto.*?retry_unlink\(\$pending.*?Continuing through the normal controller workflow.*?if \(length\(\$MGSfile\) && !\$preparedMainBranchFastPath\)/s,
	'-redoEPAfilter clears derived lifecycle state and rejoins the normal controller workflow');
like($strain, qr/\$Tcmd \.= "-redoEPAfilter 1 " if \$redoEPAfilter/,
	'generated BuildTree commands pass an explicit numeric redo-EPA value');
unlike($strain, qr/epaFilterOnly|treeCmd\.epa_filter/,
	'the retired filter-only controller and special command path are absent');
like($strain,
	qr/my \$leanOnlySubmitRequested = \$onlySubmit && !\$subJob && !\$recalcTrees.*?!\$redoSubmissionData.*?!\$repairCAT.*?!\$deepRepair.*?!\$redoEPAfilter.*?!\$reSubmit.*?my \$resumePhaseISummary = File::Spec->catfile.*?my \$leanOnlySubmitResume = \$leanOnlySubmitRequested && -s \$resumePhaseISummary.*?running extraction where required/s,
	'lean dispatch requires durable Phase-I evidence while new outputs and strict modes retain full auditing');
like($strain,
	qr/if \(\$leanOnlySubmitResume\) \{.*?\$SIdirs\{\$_\} = "\$outD\/\$_\/" for \@specis.*?mode=deferred_per_MGS.*?global_metadata_scans=0.*?\} else \{.*?evalFileStatus\(\)/s,
	'ordinary only-submit startup performs no all-MGS filesystem audit and strict modes keep the existing audit');
like($strain,
	qr/my \$fullTreeInputsInitialized = \$leanOnlySubmitResume \? 1 : 0.*?scheduling hint only.*?tree_input_sizing\.tsv.*?missing estimates will be read only when their MGS reaches submission/s,
	'lean dispatch bypasses global input sizing and treats any old sizing table only as a resource hint');
like($strain,
	qr/sub prepareMGSInputSet \{.*?return 1 if \$leanOnlySubmitResume && -s "\$tmpD\/merge\.complete\.tsv".*?collectMGSShardHandoff/s,
	'an atomic aggregate checkpoint prevents worker-shard globbing before each lean submission');
like($strain,
	qr/loadRecoveryContributionIndex\(\)\s+unless \$recoveryContributionIndexReady \|\| \$leanOnlySubmitResume.*?sub prepareMGSInputSet.*?merge\.complete\.tsv.*?loadRecoveryContributionIndex\(\)\s+if \$leanOnlySubmitResume && !\$recoveryContributionIndexReady/s,
	'the recovery-contributor index is also loaded only if a just-in-time legacy shard fallback needs it');
like($strain,
	qr/\$leanOnlySubmitResume.*?merge\.complete\.tsv.*?fileGZe\(\$rawCategory\)/s,
	'lean staged-input handling trusts the commit marker but still validates the category it immediately consumes');
like($strain,
	qr/sub outgroupRequirementLoci \{.*?\$leanOnlySubmitResume.*?selected_gene_map_deferred_validation.*?scratchMGSInputState/s,
	'outgroup demand construction uses the retained locus map without probing every staged directory first');
like($strain,
	qr/MGS_SUBMISSION:.*?opendir\(my \$resumeDirectory, \$outD2\).*?treeDone\.sto.*?terminalMarkers.*?prepareMGSInputSet\(\$MGS,\$tmpD\).*?addOutgroup2MGS\(\$MGS,\$OG,\$tmpD\).*?dispatchPendingTreeJobs/s,
	'useful completion, input, and category checks occur just in time in the same loop that submits each MGS');
like($strain,
	qr/if \(!\$leanOnlySubmitResume && \$onlySubmit && !\$subJob.*?preparedMainBranchInputSet\(.*?\$preparedMainBranchFastPath = 1/s,
	'the exhaustive prepared-input preflight remains available only outside latency-sensitive lean dispatch');
like($strain,
	qr/sub resubmitExistingTreeCommands .*?treeCmd\.sh.*?placementPending\.sto.*?skipping Mosaic, map, and catalogue loading.*?qsubSystemWaitMaxJobs\(.*?qsubSystem2\(/s,
	'direct tree-command resubmission reuses saved scripts with scheduler-capacity throttling, including EPA recovery');
my ($directTreeResume) = $strain =~
	/(sub resubmitExistingTreeCommands .*?)(?=sub markStrainWorkflowDirectory)/s;
ok(defined($directTreeResume),
	'direct tree-command resume helper is available for isolated inspection');
unlike($directTreeResume, qr/\$guide|open my \$input/,
	'direct tree-command resume scans saved output scripts instead of reading the MGS guide');
like($directTreeResume,
	qr/bsd_glob.*?my \$treeDone.*?completionMarkerTree\(\$treeDone.*?next if !\$force && -s \$treeDone && length\(\$completedTree\).*?my \$publicationResume.*?epa_result\.jplace.*?treeCmd\.epa_retry\.sh.*?epa_only/s,
	'direct tree-command resume notices a removed placed tree despite treeDone and reuses saved EPA retry scripts');
like($directTreeResume, qr/\$redoEpa && !\$publicationResume.*?next;.*?elsif \(!\$publicationResume/s,
	'forced EPA filtering selects only retained-jplace publication resumes and never EPA-only recovery');
like($directTreeResume, qr/local \$ENV\{MATAFILER_REDO_EPA_FILTER\} = 1.*?qsubSystem2/s,
	'forced filtering propagates into older saved tree commands without rewriting them');
like($directTreeResume,
	qr/\$publicationResume.*?elsif \(!\$publicationResume\).*?FNAstdof.*?FAAstdof.*?CATstdof/s,
	'a retained-jplace publication resume skips unnecessary sequence-input checks');
like($strain,
	qr/my \$requiresOutgroupReference = \$runPartI \|\| \$CatNotPrepped \|\| \$repairCAT.*?my \$initializeOutgroupReferences = sub.*?unless \(\$requiresOutgroupReference.*?readFasta\(\$refFAA.*?readFasta\(\$refFNA/s,
	'tree-only resumes load reference FASTA catalogues only for input regeneration or repair');
like($strain,
	qr/sub addOutgroup2MGS\{.*?if \(\$outputReady.*?return \(.*?my \$preparedScratchInput.*?return \(.*?my \$stageReady.*?if \(\$requiresOutgroupReference && !\$outgroupReferenceInitialized\).*?\$initializeOutgroupReferences->\(\\\@fullTreeCandidates\).*?staged category scan for \$MGS/s,
	'already-overlaid MGS bypass the catalogue, while the first raw MGS streams one shared reference set before its overlay');
like($strain,
	qr/outgroup candidate discovery.*?outgroup protein FASTA streaming.*?outgroup nucleotide FASTA streaming/s,
	'direct outgroup lookup reports candidate and sequential FASTA-streaming progress');
unlike($strain, qr/nonEpaTreeAbsences/,
	'a missing final tree no longer makes reference catalogue loading mandatory');
like($strain,
	qr/my %treeDisposition.*?\$treeDisposition\{\$epaOnlyRetry \? 'EPA-only retry job'.*?\$onlyMSA \? 'eligible MSA-only job' : 'eligible tree job'\}\+\+.*?Tree submission accounting:.*?Tree submission pass complete:/s,
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
	qr/sub addOutgroup2MGS.*?\.strain_tree_input\.outgroup\.fna.*?\.strain_tree_input\.plan\.tsv/s,
	'the controller writes only compact outgroup overlays and a plan before handing final input construction to buildTree5');
unlike($strain, qr/sort_fasta_by_locus|append_fasta_records_atomic|readFastaIDs/,
	'the serial controller no longer rewrites, sorts, or fully scans staged FASTA inputs');
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
	qr/my \$iqMemMB = int\(\$totMem \* 0\.9\).*?if \(!\$onlyMSA && \$phyloProg == 1\)\{.*?"-iqMemMB \$iqMemMB ".*?"-iqPathogen 1 " if \$iqPathogen/s,
	'within-strain IQ-TREE always uses the standard resource-limited command and enables CMAPLE only by explicit request');
like($strain,
	qr/my \$placementRequested = \$strictBackbone \? 1 : 0;.*?\$baseMemMult = 150 if \$placementRequested.*?\$minimumMemMB = \(\$placementRequested \? 10240 : 5000\) \* \$memMulti;.*?\$minimumMemMB = 10240 if \$placementRequested.*?\$totMem = \$minimumMemMB if \$totMem < \$minimumMemMB/s,
	'within-strain gives EPA-ng placement jobs a 10 GiB floor and larger input-size estimate');
unlike($strain, qr/-iqLegacy\s+1|legacyMGTK/,
	'within-strain does not expose or submit the obsolete IQ-TREE legacy-kernel flag');
like($strain,
	qr/"redo=s"\s+=> \\\$redoMode.*?qw\(none tree input all\).*?if \(\$redoMode eq 'tree'\) \{ \$recalcTrees = 1; \$onlySubmit = 1; \}.*?elsif \(\$redoMode eq 'input'\) \{ \$repairCAT = 1; \$deepRepair = 1; \$onlySubmit = 1; \}.*?elsif \(\$redoMode eq 'all'\) \{ \$redoSubmissionData = 1; \$onlySubmit = 0; \}/s,
	'canonical redo modes resolve to the established tree, input, and full-rebuild paths');
like($strain,
	qr/"recalcTrees=i"\s+=> sub.*?'-redo tree'.*?"repairCAT=i"\s+=> sub.*?'-redo input'.*?"redoSubmissionData=i"\s+=> sub.*?'-redo all'.*?deprecated redo\/repair flags/s,
	'legacy recovery options remain temporary, warning compatibility aliases');
like($strain,
	qr/-redo cannot be combined with deprecated redo\/repair flags.*?-redo tree must be launched by the main strainWithin process/s,
	'redo modes reject mixed legacy input and split-worker tree execution');
like($mgs,
	qr/"redo=s"\s+=> \\\$strainRedo.*?-redo must be one of: none, tree, input, all.*?strainSampleStats\.summary\.tsv.*?my \$strainOnlySubmit = -s \$strainPhaseISummary \? 1 : 0.*?my \@strainArguments = \(.*?'-SNPcaller', \$SNPcaller.*?'-onlySubmit', \$strainOnlySubmit.*?'-redo', \$strainRedo.*?map \{ _shell_quote\(\$_\) \} \@strainArguments/s,
	'MGS.pl forwards caller and redo through a quoted argument array and enters only-submit after durable Phase I');
unlike($mgs,
	qr/\$strain1scr .*?-(?:reSubmit|redoSubmissionData|rmMSA)\b/s,
	'MGS.pl no longer supplies redundant or deprecated strain defaults');
like($strain,
	qr/my \$runPartI = \(.*?\|\| \(\$recalcTrees && \$dirsNOTPrepped\).*?if \(\$runPartI\).*?Stage I: consensus-gene extraction/s,
	'tree recalculation reruns extraction when required per-MGS inputs are absent');
like($strain,
	qr/my \$runPartI = .*?if \(\$runPartI\).*?preComputeConsSNP\(\).*?\} else \{.*?Skipping Part I.*?reportSavedSampleStats\(\)/s,
	'completed Phase I skips the extraction-only consensus audit and reports saved statistics');
like($strain,
	qr/sub reportSavedSampleStats .*?\$sampleStatsSummaryLogName.*?scope\} .*?eq 'ALL'.*?printSampleStatsSummary\(\$allSummary\)/s,
	'a Phase-I resume loads the persisted all-worker row before rendering its accounting');
like($strain,
	qr/sub reportSavedSampleStats .*?unless \(\$header eq \$expectedHeader\).*?older schema.*?continuing tree recovery.*?return 0/s,
	'an older reporting-only sample-summary schema cannot abort Phase-II tree recovery');
unlike($strain, qr/die "Unexpected saved sample-summary header/,
	'legacy saved sample-summary headers are no longer fatal');
like($strain,
	qr/sub printSampleStatsSummary .*?STAGE I SAMPLE SUMMARY \(all workers\).*?Used MGS retained-loci histogram/s,
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
	&& index($strain, 'return (scalar(keys %samplesSeen), $genesSeen, $preparedOG, 1, 1);') >= 0,
	'legacy fully prepared Phase-II scratch inputs remain resumable without redoing their controller-side work');
ok(index($strain, 'sub preparedOutgroupLog') >= 0
	&& index($strain, 'fileGZe($log_path)') >= 0
	&& index($strain, '$publishedPrepared') >= 0
	&& index($strain, '$scratchPrepared') >= 0
	&& index($strain, '.strain_tree_input.plan.tsv') >= 0
	&& index($strain, '.strain_tree_input.shards.tsv') >= 0
	&& index($strain, 'writeMGSShardManifest') >= 0
	&& index($strain, '; outgroup ') >= 0
	&& index($strain, '; $multiSmpl samples; $ngenes genes; $numCoreL cores; $totMem MB; $memoryProfile') >= 0,
	'new outgroup preparation writes a shard manifest and reports one compact per-MGS summary');
unlike($strain, qr/Controller staged-overlay preparation|Tree input hand-off: raw FNA|Tree input: \$multiSmpl samples|Tree input: using complete published|Stage-I input: reusing controller-prepared|Recovery state: validated backbone/,
	'per-MGS progress omits verbose staging and overlay lines');
like($strain,
	qr/my \(%persistentMGSInputStateCache, %scratchMGSInputStateCache\).*?sub invalidateMGSInputState .*?delete \@persistentMGSInputStateCache.*?delete \@scratchMGSInputStateCache/s,
	'published and scratch triplet states are cached and explicitly invalidated after mutations');
like($strain,
	qr/my \$completedTree = .*?treeDone\.sto.*?fileGZs\(\$completedTree\).*?BuildTree publishes treeDone\.sto atomically.*?next;.*?fileGZe\("\$SIdirs\{\$MGS\}\/\$CATstdof"\)/s,
	'a validated completed tree bypasses compressed category-sidecar inspection');
like($strain,
	qr/my \$completedTree = "\$outD2\/phylo\/\$treeFile";.*?my \$treeCompletion = "\$outD2\/treeDone\.sto";.*?\(\$onlySubmit != 0 \|\| \$subJob\).*?BuildTree publishes treeDone\.sto atomically.*?\$completedTreeFastPaths\+\+.*?next;.*?my \$tooFewMarker/s,
	'tree-only audits prioritize the durable completion marker and primary tree before deeper MGS probes');
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
	qr/\$minimumOutgroupLoci = \$outgroupDemandMinimum\{\$MGS\} \/\/ \$MGStoolowGsThr.*?last if \$represented >= \$minimumOutgroupLoci.*?if \(\$represented < \$minimumOutgroupLoci\)/s,
	'outgroup acceptance uses the per-MGS core/broad demand floor rather than the generic eight-locus minimum');
like($strain,
	qr/my \$workerMGSSubset = \$recalcTrees.*?grep \{ \$MGSneedsExtraction\{\$_\} \} \@specis.*?'-MGSsubset', \$workerMGSSubset/s,
	'split extraction workers inherit the missing-input MGS subset');
like($strain,
	qr/Stage-I extraction scope: \$stageIScope.*?target_MGS=.*?Workers are balanced by assembly group/s,
	'split Stage I reports whether its MGS scope is explicit or recovery-driven');
like($strain,
	qr/my \$maxSubJob = -1;.*?phase1SamplesByGroup\(\).*?effectiveGroupCount.*?choose_auto_worker_count\(.*?Automatic Stage-I splitting:.*?standalone.*?target \$\{targetGroupsPerWorker\} groups\/worker/s,
	'automatic Stage-I splitting counts standalone samples as effective schedulable groups');
like($strain,
	qr/sub phase1WorkerCommand.*?'-taxonAwareLocusSelection', \$taxonAwareLocusSelection.*?'-prepareMosaicLoci', \$prepareMosaicLoci.*?'-SNPcaller', \$SNPcaller/s,
	'split extraction workers inherit caller, Mosaic preparation, and Phase-I locus-selection controls');
my ($phase1WorkerCommandSource) = $strain =~
	/(sub phase1WorkerCommand \{.*?^\})/ms;
ok(defined($phase1WorkerCommandSource),
	'the Phase-I worker command builder is available for option-contract auditing');
my %phase1WorkerFlag = map { $_ => 1 }
	$phase1WorkerCommandSource =~ /'(-[A-Za-z][A-Za-z0-9]*)'/g;
my @requiredPhase1WorkerFlags = qw(
	-GCd -outD -MGS -clusterID -submit -onlySubmit -maxSubJob
	-MGSminGenesPSmpl -multiGeneSmplMax -conspGeneSmplMax
	-minBadLociPSmpl -presortGenes -maxGenes -treeLocusBudget
	-taxonAwareLocusSelection -disableQC -breakpointGeneFlank
	-abundanceMinLoci -abundanceMinFold -abundanceMaxFold
	-abundanceMaxModifiedZ -prepareMosaicLoci -flushEvery -MGset
	-SNPcaller -minSNPDepth -minSNPCallQual -forceSNPcalls
	-preCompConsSNP -skipIndels -SNPadaptiveQual -SNPdepthFilterScale
	-SNPindelRangeFilt -tmpD -mosaicLoci -MGSabundance -MGSsubset
);
is_deeply(
	[grep { !$phase1WorkerFlag{$_} } @requiredPhase1WorkerFlags],
	[],
	'all mandatory and conditional Phase-I execution flags are propagated to workers',
);
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
	qr/my \$publishedInputsReady = !\$epaOnlyRetry\s*&& !exists\(\$legacyLocusMGS\{\$MGS\}\).*?persistentMGSInputState\(\$MGS\) eq 'complete'.*?if \(\$recalcTrees\).*?unless \(\$publishedInputsReady\).*?\$scratchInputsReady = prepareMGSInputSet\(\$MGS,\$tmpD\).*?unless \(\$publishedInputsReady \|\| \$scratchInputsReady\).*?no recoverable inputs for recalculation.*?resetMGSTreeOutputs\(\$outD2, \$MGS\)/s,
	'tree outputs are reset only after complete published or recoverable staged per-MGS inputs are verified');
ok(index($strain, 'sub prepareMGSInputSet') >= 0
	&& index($strain, 'collectMGSShardHandoff($MGS, $tmpD)') >= 0
	&& index($strain, 'return combineMGSgenesDir($MGS, $tmpD);') >= 0
	&& index($strain, '$scratchInputsReady ||= prepareMGSInputSet($MGS,$tmpD);') >= 0
	&& index($strain, '"-stagedInputDir "') >= 0,
	'normal tree submission prefers worker-shard handoff and retains aggregate merging as a compatibility fallback');
like($strain,
	qr/staged input sets recovered for -redo tree: \$recalcScratchRecovered/,
	'tree submission accounting reports staged redo recovery separately from skipped dispositions');
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
	qr/my \@rawCategorySources = \$shardHandoff.*?parts\}\{category\}\{path\}.*?for my \$categorySource \(\@rawCategorySources\).*?gzipopen\(\$categorySource.*?\.strain_tree_input\.outgroup\.cat\.tsv/s,
	'within-strain outgroup handling scans category shards directly and emits only small overlays');
ok(index($strain, '.strain_tree_input.plan.tsv') >= 0
	&& index($strain, 'strain-staged-input-v1\noutgroup\t$OG\nmgs\t$MGS\n') >= 0
	&& index($strain, '.strain_tree_input.shards.tsv') >= 0
	&& index($strain, q{my @line = ('strain-shard-input-v1');}) >= 0,
	'within-strain records an explicit worker-shard finalization contract without publishing aggregates');
like($strain,
	qr/"flushEvery=i"\s+=> \\\$appendWriteTrigger.*?%outgroupGeneCache = \(\).*?'-flushEvery', \$appendWriteTrigger/s,
	'within-strain extraction exposes its buffer bound to workers and releases per-MGS outgroup caches');
like($strain,
	qr/phase1_flush_samples => 50.*?my \$appendWriteTrigger = \$FILTER_DEFAULT\{phase1_flush_samples\}/s,
	'within-strain bounds the default Stage-I buffer to fifty sample rows');
like($strain,
	qr/sub appendWriteMGSgenes.*?my \@pendingMGS = grep.*?Flushing buffered MGS records: \$pendingCount MGS to durable scratch.*?stepProgress\("buffered MGS publication"/s,
	'within-strain reports progress while durable worker shards are published');
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
	qr/sub phase1EstimatedInputBytes.*?fileGZs\(\$nominal\).*?sub phase1SampleWorkEstimate.*?phase1EstimatedInputBytes\(\$readyNT\).*?phase1EstimatedInputBytes\(\$vcf\).*?'regenerate'.*?sub phase1GroupWorkEstimates/s,
	'Phase I estimates FASTA scan size and penalizes consensus regeneration');
like($strain,
	qr/phase1GroupWorkEstimates\(\$samplesByGroup\).*?balance_assembly_groups\(\$samplesByGroup, \$maxSubJob, \$groupWork\).*?writePhase1WorkerPlan\("\$LOGDIR\/phase1_worker_plan\.tsv".*?my \$worker = \$workerForGroup->\{\$group\}.*?next unless \$worker == \$subJob.*?\$plannedSamples \+= scalar\(\@\{\$samplesByGroup->\{\$_\}\}\).*?estimated work/s,
	'split extraction keeps assembly groups intact while balancing estimated work and auditing its plan');
like($strain,
	qr/pre-restricted to .*?sample driver\(s\) with target loci/s,
	'split-worker diagnostics distinguish post-index usable sample drivers from assembly groups');
like($strain,
	qr/readGenesSample_Singl\(\s*\$sm, \$writeLink, \$sttime, .*?\$appCnt, undef, .*?sampleStatsSeen.*?\$\{\$bufferedSamplesRef\}\+\+.*?appendWriteMGSgenes\(\$writeLink\)/s,
	'expanded assembly-group output is accounted once and flushed by sample to retain the RAM bound');
like($strain,
	qr/if \(\$mySamplesHR\).*?\$unrepresentedWorkerLoci\+\+.*?unless \$maxSubJob/s,
	'split-worker sparsity is summarized instead of reported as missing catalogue data');
like($strain,
	qr/-withinSpecies 1 -relativeNTFraction \$relativeNTFraction .*?-NTfiltPerGene \$GeneLengthMin .*?-GeneLengthIncludeMin \$GeneLengthIncludeMin .*?-GenesPerSpecies \$GenesPerSpecies/s,
	'unfinished trees explicitly pass the named strain coverage filters to buildTree');
unlike($strain, qr/-NTfilt \$relativeNTFraction/,
	'strain workflow does not emit the retired ambiguous NTfilt option');
like($strain,
	qr/my \$GenesPerSpecies = 0\.2;.*?my \$GeneLengthMin = 0\.3;.*?my \$GeneLengthIncludeMin = 0\.03;.*?my \$relativeNTFraction = 0\.1;.*?\$placementGenesPerSpecies = 0.04; \$placementRelativeNTFraction = 0.03;.*?my \$taxonAwareLocusSelection = 1;.*?"GeneLengthIncludeMin=f" => \\\$GeneLengthIncludeMin.*?"taxonAwareLocusSelection=i" => \\\$taxonAwareLocusSelection.*?-GeneLengthIncludeMin \$GeneLengthIncludeMin.*?-taxonAwareLocusSelection \$taxonAwareLocusSelection/s,
	'strainWithin separates high-threshold QC from lower MSA inclusion while retaining balanced placement filters');
like($strain,
	qr/my \$taxonAwareRescueMinPrevalence = 0\.8;.*?"taxonAwareRescueMinPrevalence=f" => \\\$taxonAwareRescueMinPrevalence.*?-taxonAwareRescueMinPrevalence \$taxonAwareRescueMinPrevalence/s,
	'strainWithin exposes and forwards the broad-locus rescue prevalence guard');
like($strain,
	qr/my \$rateMergePartitions = 1;.*?my \$rateMergeMaxBins = 8;.*?my \$rateMergeTargetSites = 30_000;.*?my \$rateMergeMinLoci = 20;.*?my \$rateMergeMinSites = 20_000;.*?"rateMergePartitions=i" => \\\$rateMergePartitions.*?-rateMergePartitions \$rateMergePartitions.*?-rateMergeMaxBins \$rateMergeMaxBins.*?-rateMergeTargetSites \$rateMergeTargetSites.*?-rateMergeMinLoci \$rateMergeMinLoci.*?-rateMergeMinSites \$rateMergeMinSites/s,
	'strainWithin enables deterministic rate merging and forwards all bin controls');
like($strain,
	qr/maximum_genes_per_sample => 600.*?maximum_tree_loci => 400.*?\$taxonAwareGeneBudget = \$treeLocusBudget < \$presortGenes.*?taxonAwareLocusBudgets\(\$taxonAwareGeneBudget\).*?-taxonAwareMaxLoci \$taxonAwareMaxLoci.*?-taxonAwareCoreLoci \$taxonAwareCoreLoci.*?-taxonAwareCandidateExtra \$taxonAwareCandidateExtra.*?sub taxonAwareLocusBudgets.*?\$maximumLoci \* 0\.8.*?\$maximumLoci \* 0\.3/s,
	'strainWithin scales 80% core, 20% rescue capacity, and 30% QC backfill to its effective gene budget');
like($strain,
	qr/my \$strictBackbone = 0;.*?my \$strictBackboneFraction = 0\.35;.*?my \$strictBackboneMinSamples = 3;.*?my \$placementMinOverlap = 10_000;.*?"placeOnBackbone=i"\s+=> sub.*?\$placeOnBackboneSpecified = 1.*?"strictBackbone=i"\s+=> sub.*?'-placeOnBackbone'.*?"strictBackboneFraction=f"\s+=> \\\$strictBackboneFraction.*?"placementMinOverlap=i"\s+=> \\\$placementMinOverlap.*?-placeOnBackbone cannot be combined/s,
	'strainWithin exposes opt-in backbone placement and retains a deprecated switch alias');
like($strain,
	qr/-placeOnBackbone \$strictBackbone .*?if \(\$strictBackbone\) \{.*?-placementGenesPerSpecies.*?-strictBackboneFraction \$strictBackboneFraction .*?-placementMinOverlap \$placementMinOverlap/s,
	'strainWithin forwards placement-only controls only when backbone placement is enabled');
like($build_tree,
	qr/"placeOnBackbone=i" => sub.*?"strictBackbone=i" => sub.*?-placeOnBackbone cannot be combined.*?Option -strictBackbone is deprecated/s,
	'buildTree exposes the canonical placement switch with a guarded compatibility alias');
like($build_tree,
	qr/if \(\$strictBackbone\) \{\s*my \$backboneEligibility = classifyTaxonAwareCoverageEligibility.*?my \$placementEligibility = classifyTaxonAwareCoverageEligibility.*?\} else \{.*?placement disabled/s,
	'BuildTree skips backbone and placement eligibility filtering when placement is disabled');
like($strain,
	qr/sub writeGeneLengthSampleSummary .*?gene_length_filter\.samples\.tsv.*?"\$\{mgs\}:\$_".*?strainGeneLengthFilter\.samples\.tsv.*?gene_length_sample_audit\\t\$geneLengthSampleSummary/s,
	'strainWithin consolidates per-MGS length-gate decisions into a sample-wise run report');
like($strain,
	qr/my \$epaPendantOutlierFactor = 5;.*?my \$epaPendantMinThreshold = 0\.02;.*?"epaPendantOutlierFactor=f" => \\\$epaPendantOutlierFactor.*?"epaPendantMinThreshold=f" => \\\$epaPendantMinThreshold.*?-epaPendantOutlierFactor \$epaPendantOutlierFactor.*?-epaPendantMinThreshold \$epaPendantMinThreshold/s,
	'strainWithin enables and forwards adaptive EPA pendant-branch outlier QC');
like($strain,
	qr/-tmpSubdir .*?strain_within\/\$MGS.*?-stagedInputDir .*?\$tmpD.*?-completionMarker .*?\$treeStone/s,
	"tree jobs pass lifecycle paths to buildTree5 as ordinary options");
unlike($strain, qr/sub treeInputPrecopyCommand|staged_inputs=\(\)|mapfile -d|ready_test/,
	"strain_within no longer generates Bash input-publication logic");
unlike($strain, qr/\$\{TMPDIR\}\/strain_within|my \$postCmd|touch "?\.shellQuote\(\$treeStone\)/,
	"tree commands contain neither shell TMPDIR expansion nor shell checkpoints");

like($strain,
	qr/sub phase1WorkerCommand.*?Stage-I workers receive extraction, consensus, and extraction-relevant.*?sub recoverCompletedSplitPhaseI.*?phase1WorkersNeedingRetry.*?Resubmitting invalid Phase-I worker.*?phase1_worker_repair\.queue\.tsv/s,
	'live and resumed Phase I share extraction-only worker commands and durable targeted repair');
like($strain,
	qr/\$noGeneLimit = 1 if \$maxNGenes <= 0;.*?\$maxNGenes = 0 if \$noGeneLimit;.*?sub phase1WorkerCommand.*?'-presortGenes', \$presortGenes, '-maxGenes', \$maxNGenes,.*?'-disableQC', \$disableQC,/s,
	'unlimited extraction is canonicalized to maxGenes zero before worker commands are built');
unlike($strain,
	qr/sub phase1WorkerCommand \{.*?'-noGeneLimit'.*?^\}/ms,
	'split-worker commands do not forward the deprecated noGeneLimit alias');
like($strain,
	qr/No automatic full-tree resubmission was attempted.*?sub writeTreeFailureAudit.*?failed_missing_output.*?valid_no_tree.*?placement_pending/s,
	'tree outcomes are classified and quarantined without automatic tree resubmission');
like($strain,
	qr/my \$unresolvedInputs = validateTreeInputResolution\(\);.*?if \(\$unresolvedInputs\).*?tree_outcomes_quarantined=\$incompleteTreeOutcomes.*?exit\(0\);.*?if \(\$incompleteTreeOutcomes\).*?proceeding with downstream strain analysis for completed trees.*?qsubSystem\(\$LOGDIR\."strainAnalysis2\.sh"/s,
	'quarantined tree outcomes do not block step two once all tree inputs are resolved');
like($strain,
	qr/sub validateTreeInputResolution.*?my \(\$ready, \$terminal, \$excluded\) = \(0, 0, 0\).*?Tree-input resolution audit: ready=\$ready/s,
	'tree-input resolution summaries initialize zero-valued counters');
like($strain2,
	qr/my \@nonTreeOutcomeMarkers = qw\(.*?tooFewSamples\.sto.*?noRecoverableLoci\.sto.*?noTree\.sto.*?placementPending\.sto.*?\);.*?my \@outcomeMarkers = grep.*?if \(\@outcomeMarkers\).*?\$terminalTreeMGS\+\+;.*?next;/s,
	'step two explicitly skips MGS with valid no-tree or placement-pending markers');
like($strain2,
	qr/my \$treeCompletion = "\$FMGpD\/\$entry\/treeDone\.sto";.*?if \(-s \$treeCompletion\).*?completedTreeSize.*?completionMarkerFastPaths\+\+.*?next;.*?my \@outcomeMarkers/s,
	'step two uses the same durable completion-marker and primary-tree fast path as strain_within');
like($strain,
	qr/sub lifecycleMarkerReason.*?\^reason\\t/s,
	'BuildTree lifecycle-marker reasons have one reusable parser');
like($strain,
	qr/\$MGSnoTreeReason\{\$MGS\} = lifecycleMarkerReason\(\$buildTreeTerminalMarker/s,
	'future strain resumes retain the BuildTree terminal reason when skipping an MGS');
like($strain,
	qr/valid_no_tree_buildtree.*?lifecycleMarkerReason\("\$SIdirs\{\$MGS\}\/noTree\.sto"/s,
	'tree-input resolution audits include the specific BuildTree terminal reason');
like($strain,
	qr/\$totMem = int\(\$totMem \* 2\).*?\$numCoreL = 1.*?-epaOnly 1.*?\$outD2\/treeCmd\.epa_retry\.sh/s,
	'an EPA-only retry gets a one-core doubled-memory job and explicit BuildTree mode');
like($strain,
	qr/my \$treeOOMMaxMemGB = 512;.*?my \$treeOOMMaxMemGBSpecified = 0;.*?"treeOOMMaxMemGB=f" => sub \{.*?\$treeOOMMaxMemGBSpecified = 1;.*?getProgPaths\("maxMF4mem", 0\).*?-treeOOMMaxMemGB must be positive.*?maximum_rounds => \$treeOOMRetryRounds/s,
	'automatic tree OOM recovery uses maxMF4mem by default while retaining an explicit per-run override');
like($site_config_template, qr/^maxMF4mem\t512$/m,
	'the installed site-config template sets the shared strain OOM ceiling to 512 GiB');
unlike($internal_config, qr/^downloadQueue\t/m,
	'the internal program config leaves the archive queue to the site config');
like($strain,
	qr/sub retryOOMTreeJobs.*?for my \$round \(1 \.\. \$maximumRounds\).*?oom_jobs.*?next_oom_retry_memory_mb.*?epaOnlyRetryReady\(\$mgsDirectory, 1\).*?-epaThreads\\s\+\\d\+\/\$1-epaThreads 1.*?treeCmd\.epa_retry\.sh.*?qsubSystemJobAlive/s,
	'only accounting-confirmed OOM jobs are retried and EPA-stage retries use one thread');
like($strain,
	qr/sub dispatchPendingTreeJobs.*?submission_record => \{ %\{\$record\} \}.*?sub retryOOMTreeJobs/s,
	'tree submission accounting retains the exact command record needed for bounded OOM retries');
like($strain,
	qr/sub epaOnlyRetryReady.*?\$onlySubmit.*?IQtree_allsites\.backbone\.treefile.*?strict_backbone\.samples\.tsv.*?status\\tplacement_pending.*?explicit_pending/s,
	'a tree-only resume validates retained placement inputs and recognizes an explicit pending marker');
like($strain,
	qr/sub epaOnlyRetryReady.*?IQtree_allsites\.treefile.*?return '' if -s \$finalTree.*?legacy_missing_final.*?sub prepareEpaOnlyRetryState.*?clear stale completion missing final placed tree.*?create legacy placement-pending marker/s,
	'a legacy retained backbone without the final non-backbone tree is prepared for isolated EPA recovery');
like($strain,
	qr/my \$epaOnlyRetry = exists\(\$MGSepaOnlyRetry\{\$MGS\}\).*?my \$epaRecovery = \$epaOnlyRetry.*?if \(!\$epaRecovery && exists \$MGSnoTree\{\$MGS\}\).*?if \(!\$epaRecovery && exists\(\$ConspecificMGS\{\$MGS\}\)/s,
	'a validated EPA-only retry bypasses later historical no-tree and multicopy filters');
like($strain,
	qr/Placement-only recovery has already paid.*?my \@epaRecoveryMGS = grep.*?my \@fullTreeMGS = grep.*?\@specis = \(\@epaRecoveryMGS, \@fullTreeMGS\).*?my \$epaQueueBoundary = scalar\(\@epaRecoveryMGS\).*?Validated EPA-only recovery queue/s,
	'validated EPA-only retries are queued ahead of ordinary full-tree retries');
like($strain,
	qr/strain_within\.state\.tsv.*?strain_within\.heartbeat\.tsv.*?strain_within\.failure\.tsv.*?sub writeStrainWorkflowState.*?sub writeStrainWorkflowHeartbeat.*?sub writeStrainWorkflowFailure/s,
	'strain workflow stores running, completed, and failed status in one state record');
unlike($strain,
	qr/preflightStrainWorkflow|preflight_executable|preflight_directory|filesystem_capacity/,
	'strain workflow does not preflight environment-wrapped commands as local executables');

like($strain,
	qr/require_complete_linkage => 1.*?Mosaic complete-linkage protection rejected/s,
	'strain extraction requires pairwise confirmation throughout multi-seed Mosaic loci');
like($strain,
	qr/\$workflowStatePath = File::Spec->catfile\(\$LOGDIR.*?\$SNPconsLOGs = "\$LOGDIR\/SNPconsCalls.*?my \$final = "\$LOGDIR\/\$sampleStatsLogName".*?my \$summary = "\$LOGDIR\/\$sampleStatsSummaryLogName".*?my \$final = "\$LOGDIR\/\$recoveryLogName".*?my \$summary = "\$LOGDIR\/\$summaryLogName"/s,
	'worker state, SNP logs, sample statistics, recovery accounting, and summaries share LOGandSUB');
like($strain,
	qr/script => \$epaOnlyRetry \? "\$outD2\/treeCmd\.epa_retry\.sh".*?\$onlyMSA \? "\$outD2\/treeCmd\.msa_only\.sh" : "\$outD2\/treeCmd\.sh".*?qsubSystem\(\$LOGDIR\."strainAnalysis2\.sh"/s,
	'per-MGS normal, MSA-only, and EPA-retry scripts use distinct paths while downstream strain-analysis scripts use LOGandSUB');
like($strain,
	qr/sub migrateLegacyOperationalLogs.*?strain_within.*?SNPconsCalls.*?strainSampleStats.*?strainRecovery.*?migrate legacy strain log.*?migrateLegacyOperationalLogs\(\)/s,
	'legacy top-level operational logs are safely migrated into LOGandSUB');
like($strain,
	qr/sub writeSelectionAttritionSummary.*?selection_attrition\.tsv.*?strainSelectionAttrition\.tsv/s,
	'strain summary aggregates completed BuildTree attrition reports');
like($strain,
	qr/sub writeMGSSampleHistograms.*?backbone_samples.*?placement_samples.*?strict_backbone\.samples\.tsv.*?strainMGSSampleCounts\.tsv.*?role\\tlower\\tupper\\tbin\\tMGS_count\\tfraction.*?qw\(backbone placement\)/s,
	'across-MGS sample histograms report backbone and placement distributions separately');
like($strain,
	qr/writeMGSSampleHistograms\(\).*?MGS_sample_counts.*?MGS_sample_histogram.*?for my \$role \(qw\(backbone placement\)\).*?\$\{role\}_samples_per_MGS/s,
	'the run summary links exact sample counts and both role-specific histograms');

like($build_tree, qr/if \(\$numSeq < 3\)/,
	'three-sample MGS accepted by the wrapper are retained for a minimal tree');
like($build_tree,
	qr/enabled => 0,.*?my \$strictBackbone = \$BACKBONE_DEFAULT\{enabled\}/s,
	'buildTree5 keeps strict-backbone EPA placement disabled by default');
like($build_tree,
	qr/"GeneLengthIncludeMin=f" => \\\$ntFracGeneInclude.*?geneLengthIncludeByGene.*?qualification_sequences => \\%geneLengthQCSequence/s,
	'buildTree5 aligns lower-threshold recovery data while keeping final sample QC on high-threshold sequences');
like($build_tree,
	qr/-minGoodPosFrac", \(\$cogCats ne '' \? \$ntFracGeneInclude : 0\.6\)/,
	'category-based MSA cleaning does not impose the obsolete 60% coverage floor on recovered fragments');
like($build_tree,
	qr/my \$retainedJplace = File::Spec->catfile.*?if \(\$continue && \$dedicatedBackbone && !-s \$primaryTree.*?-s \$retainedJplace\).*?read_epa_jplace.*?reapplying placement filtering.*?else \{.*?runEpaNgPlacement.*?filter_epa_placement_outliers.*?write_epa_placed_tree/s,
	'normal BuildTree continuation reuses a retained jplace when its placed tree is missing');
like($build_tree,
	qr/"redoEPAfilter:i" => sub.*?if \(\$subsetSmpls >0\).*?if \(\$redoEPAfilter\).*?runRedoEpaFilter\(.*?exit\(0\).*?warn "MSAprobs.*?prepGenoDirs/s,
	'BuildTree executes forced EPA filtering before normal workflow startup');
unlike($build_tree, qr/epaFilterOnly|runEpaFilterOnly/,
	'BuildTree has no separate filter-only option or execution path');

done_testing();
