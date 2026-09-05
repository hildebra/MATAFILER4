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

# Exercise the Phase-I guide identity independently of the full controller. A
# split worker validates the contract before prepRun() assigns $outD, whereas
# the parent records it afterwards. Both moments must select the run-local
# sorted guide even when an older catalogue-side .srt also exists.
my ($phase1_fingerprint_helpers) = $strain =~
	/(sub phase1PathStatComponent \{.*?^\}\n\nsub phase1GuideStatFingerprint \{.*?^\})\n\nsub phase1CatalogStatFingerprint/ms;
BAIL_OUT('Cannot extract Phase-I guide fingerprint helpers')
	unless defined $phase1_fingerprint_helpers;
my $phase1_helpers_loaded = eval <<"PERL";
package TestPhase1GuideFingerprint;
use strict;
use warnings;
use Cwd qw(abs_path);
use Digest::SHA qw(sha256_hex);
use File::Basename qw(basename);
use File::Spec;
our \$phase1InputContractVersion = 3;
our \$outD = '';
sub resolveExistingFile {
	my (\$path) = \@_;
	return -f \$path ? \$path : undef;
}
$phase1_fingerprint_helpers
1;
PERL
ok($phase1_helpers_loaded, 'Phase-I guide fingerprint helpers load independently')
	or diag($@);
my $contract_catalogue = File::Spec->catdir($tmp, 'contract-catalogue');
my $contract_output = File::Spec->catdir($tmp, 'contract-output');
mkdir $contract_catalogue or die "Cannot create $contract_catalogue: $!";
mkdir $contract_output or die "Cannot create $contract_output: $!";
my $contract_guide = File::Spec->catfile(
	$contract_catalogue, 'SB.clusters.core');
write_file($contract_guide, "MGS.1\tgene1\n");
write_file(File::Spec->catfile($contract_catalogue, 'SB.clusters.obs'),
	"MGS.1\t1\n");
write_file("$contract_guide.srt", "MGS.catalogue\tgene-old\n");
write_file("$contract_guide.srt.gene2MGS", "gene-old\tMGS.catalogue\n");
my $staged_sorted = File::Spec->catfile(
	$contract_output, 'SB.clusters.core.srt');
write_file($staged_sorted, "MGS.1\tgene1\n");
write_file("$staged_sorted.gene2MGS", "gene1\tMGS.1\n");
$TestPhase1GuideFingerprint::outD = $contract_output;
my $parent_guide_fingerprint =
	TestPhase1GuideFingerprint::phase1GuideStatFingerprint($contract_guide);
$TestPhase1GuideFingerprint::outD = '';
my $worker_guide_fingerprint =
	TestPhase1GuideFingerprint::phase1GuideStatFingerprint(
		$contract_guide, undef, $contract_output);
is($worker_guide_fingerprint, $parent_guide_fingerprint,
	'pre-initialization worker and initialized parent fingerprint the same run-local guide');
isnt(
	TestPhase1GuideFingerprint::phase1GuideStatFingerprint($contract_guide),
	$parent_guide_fingerprint,
	'a catalogue-side sorted guide has a distinct identity and cannot be selected accidentally');
my $mgs = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'MGS.pl'));
my $build_tree = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'phylo', 'buildTree5.pl'));
my $build_tree_locus_module = slurp(File::Spec->catfile($Bin, q{..}, q{Mods}, q{MGSLocus.pm}));
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
	qr/msaOnly\.complete\.tsv.*?my \$msaOnlyOutput = "\$outD2\/MSA".*?-MSAprogram \$MSAprog -onlyMSA \$onlyMSA/s,
	'strainWithin forwards MSA-only mode and audits its durable alignment outcome');
like($strain,
	qr/sub msaOnlyArtifactsReady.*?status.*?msa_complete.*?\$name =~ \/\^MSAli\/.*?fna.*?fileGZs/s,
	'strainWithin validates MSA-only completion from non-merged per-locus artifacts');
like($strain,
	qr/if \(\$onlyMSA.*?strain_within_2\.2\.pl were not launched.*?exit\(0\);.*?my \$MGSabundance/s,
	'MSA-only controller runs stop before tree-dependent postprocessing');
like($build_tree,
	qr/"onlyMSA=i" => \\\$onlyMSA.*?filterMSA\(.*?runMSAFix\(\$tmpOutMSA.*?publishCompressedMSAArtifact\(\$finOutMSA.*?if \(\$onlyMSA.*?msa_complete.*?exit\(0\);.*?POST-ALIGNMENT WORKFLOW.*?mergeMSAs/s,
	'BuildTree exits after localized per-locus processing and before combined-MSA postprocessing or concatenation');
ok(index($strain, 'my $ensureLocusMSAs = !$onlyMSA && !$rmMSA') >= 0
	&& index($strain, q{$Tcmd .= "-ensureLocusMSAs 1 " if $ensureLocusMSAs;}) >= 0
	&& index($strain, 'msa_only => $msaOnlyJob') >= 0
	&& index($strain, 'treeCmd.msa_retain.sh') >= 0,
	'-rmMSA 0 queues a marker-backed retained-locus recovery without clearing the completed tree'
);
ok(index($build_tree, 'my $locusMSARecovery = $ensureLocusMSAs && $treesDone') >= 0
	&& index($build_tree, 'my $calcMSA = $locusMSARecovery ||') >= 0
	&& index($build_tree, 'unless $locusMSARecovery;') >= 0
	&& index($build_tree, 'localized per-locus alignments retained after successful tree construction') >= 0,
	'current BuildTree backfills only policy-matched locus MSAs, preserves tree completion, and records retained artifacts'
);

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
like($strain,
	qr/waitPhase1WorkersWithOOMScan\(.*?\) if \@jobsMain && \$doSubmit;.*?qsubSystemJobAlive\(\\\@pendingJobs, \$options\)\s+if \@pendingJobs && \$options->\{doSubmit\};/s,
	'dry runs do not poll scheduler jobs that were never submitted');
like($strain, qr/\$nxtCmd \.= "-submit \$doSubmit ";.*?-qsubSystem/s,
	'postprocessing inherits submission state and the selected queue backend');
like($strain, qr/\$nxtCmd \.= "-MGSphylo "\.shellQuote\(\$treeFile\).*?if \$treeFile ne ""/,
	'postprocessing receives the source MGS tree for outgroup recovery');
like($strain,
	qr{loadTreeOutgroupCandidates\(\$targetMGS\).*?sub loadTreeOutgroupCandidates .*?" --all --max-candidates 1".*?open my \$bulk.*?\$TreeOutgroupCandidatesBulkLoaded = 1.*?sub treeOutgroupCandidates .*?A failed bulk call.*?--preferred-tip .*?--max-candidates 1}s,
	'Phase II asks one bulk R call, and its Mosaic-aware individual fallback, for exactly one outgroup');
like($strain,
	qr{tempfile\(.*?strain_mosaic_outgroups.*?TMPDIR => 1, UNLINK => 1.*?print \{\$preferredFh\} "\$MGS\\t\$PreferredOutgroup\{\$MGS\}\\n".*?" --preferred "}s,
	'the bulk R call receives all Mosaic preferences through one automatically removed temporary file');
like($strain,
	qr{my \(\$MGS, \$decision, \$preferred, \$preferredDistance, \$cutoff, \$candidateText\).*?split /\\t/, \$line, 6.*?Mosaic decisions:}s,
	'Perl imports the authoritative R ordering and summarizes Mosaic plausibility decisions');
like($strain,
	qr{if \(length\(\$treeFile\)\).*?returned at most one authoritative outgroup.*?push \@candidates, treeOutgroupCandidates\(\$MGS\).*?\$SelectedOutgroup\{\$MGS\} = \$candidates\[0\]}s,
	'the reference preload stores only the one authoritative R-selected outgroup');
like($strain,
	qr{my \@requiredLoci = sort grep.*?\$OG = \$SelectedOutgroup\{\$MGS\} // ''.*?for my \$locus \(\@requiredLoci\).*?Predetermined outgroup \$OG supplies}s,
	'the final outgroup chooser validates only the controller-selected outgroup');
unlike($strain, qr/sub addOutgroup2MGS\{.*?for my \$candidate \(\@candidates\).*?my \(\$overlayFNA/s,
	'the per-MGS outgroup path has no candidate fallback loop');
like($neighbor_tree_r,
	qr{identical\(target, "--all"\).*?ape::cophenetic\.phylo\(tree\).*?for \(tip in tree\$tip\.label\).*?ranked\$decision.*?paste\(ranked\$candidates}s,
	'neighborTree bulk mode computes distances once and emits authoritative decision and candidate columns per tree tip');
like($neighbor_tree_r,
	qr{--preferred.*?--preferred-tip.*?ranked_neighbors <- function.*?stats::quantile.*?nearestDistance \* preferredNearestFactor.*?preferredDistance > cutoff.*?candidateNames\[candidateNames != preferred\].*?c\(preferred, candidateNames}s,
	'a plausible Mosaic outgroup is promoted while an extreme-distance proposal is excluded from the R result');
like($neighbor_tree_r,
	qr{--max-candidates.*?limit_candidates <- function.*?head\(candidates, maxCandidates\).*?candidates = limit_candidates\(candidateNames\)}s,
	'neighborTree can cap both ordinary and preferred decisions to one selected outgroup');
like($strain2, qr/"MGSphylo=s"\s*=>\s*\\\$MGSphylo.*?sub resolveOutgroup .*?data\.log.*?treeCmd\.sh.*?MGSphylo/s,
	'postprocessing preserves logged or saved outgroups and falls back to the source MGS tree');
like($strain, qr/sub assertSafeWorkflowRemoval .*?resolved_default.*?Refusing to remove unowned custom output directory/s,
	'custom recursive output removal requires a workflow-owned directory');
like($strain2,
	qr/\$waitForAnalysis->\('strainStats'\);.*?my \$shouldCombineStrainStats = \$forceStrainStats \|\| \$strainTaskCount > 0 \|\| !-s \$RsummaryTab;.*?if \(\$shouldCombineStrainStats\) \{.*?combineResults\(0\);.*?\} else \{.*?Reusing existing combined strainStats overview/s,
	'strainStats stores are combined only when missing or newly stale before the population phase is awaited');
like($strain2, qr/our \$version = 0\.50;/,
	'Results-directory behavior change has an explicit postprocessing version');
like($strain2,
	qr/my \$resultsDir = "\$FMGpD\/Results";.*?my \$RsummaryTab = "\$resultsDir\/strainStats\.tsv";.*?my \$popGenSummaryTab = "\$resultsDir\/popGenStats\.tsv";.*?my \$popGenSubsampleSummaryTab = "\$resultsDir\/popGenStats\.subsamples\.tsv";/s,
	'combined overview validation follows the MG-STK Results subdirectory contract');
like($strain2,
	qr/\$combineResultsR --path .*?shellQuote\(\$FMGpD\).*?--outDir .*?shellQuote\(\$FMGpD\)/s,
	'combineResults uses only supported path and output options; it auto-detects population-genetics stores');
unlike($strain2, qr/--outDir .*?shellQuote\(\$resultsDir\)/s,
	'the MG-STK output root is not changed to Results, which would create Results/Results');
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
	qr/"popGenStats=i".*?"popGenStrictOutgroup=i".*?"popGenGeneticCode=i".*?"popGenCodonStart=i".*?"popGenSeed=i".*?"popGenLegacyTextOutput=i".*?my \$popGenStatsR = \$doPopGenStats \? getProgPaths\("pogenStats"\).*?popGenStats\.output\.Rds.*?\$popGenStatsR .*?\$destBaseD, \$refMap, \$destD.*?--subsample .*?\$popGenSubsample.*?--ncore \$jobCores.*?--genetic-code \$popGenGeneticCode.*?--codon-start \$popGenCodonStart.*?--seed \$popGenSeed.*?--individual-column .*?\$individualVar.*?--category .*?shellQuote\(\$popGenCategory\) if length\(\$popGenCategory\).*?--outgroup .*?\$OG.*?--strict-outgroup.*?--legacy-text-output.*?\$strainFile = "\$destD\/IQtree_allsites\.strains\.txt".*?PopGenStats owns validation and fallback behavior.*?\$isolatedAnalysisBlock->\(\s*"\$popGenCommand --strain-file .*?\$strainFile.*?\$popGenStore/s,
	'population genetics forwards parallel, outgroup, codon, seed, identity, categories, and the strain-file path to its durable RDS workflow');
like($strain2,
	qr/\$popGenStatsReady = !\$doPopGenStats \|\| -s \$popGenStore.*?if \(\$doPopGenStats\) \{.*?\$waitForAnalysis->\('popGenStats'\);.*?my \$shouldCombinePopGenStats = \$forcePopGenStats \|\| \$popGenTaskCount > 0 \|\| !-s \$popGenSummaryTab;.*?if \(\$shouldCombinePopGenStats\) \{.*?combineResults\(1\);.*?\} else \{.*?Reusing existing combined PopGenStats overview.*?combineResults\.R did not produce the population overview table \$popGenSummaryTab/s,
	'existing population RDS stores and aggregate tables are reused, while missing or newly stale tables are combined after the population phase');
like($strain2,
	qr/\$popGenSubsampleSummaryTab = "\$resultsDir\/popGenStats\.subsamples\.tsv".*?Combined subsampled population-genetics overview/s,
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
	qr/sub visualizeSignPhylos\{.*?\$phyloFigureStone = "\$resultsDir\/phyloFigures\.sto";.*?\$cmdPic \.= "touch ".*?if \(!-e \$phyloFigureStone\).*?qsubSystem\(\s*"\$resultsDir\/phyloFigures\.sh", \$cmdPic, 1, "24G", "phyloFigures".*?return \(\$dep, \$phyloFigureStone\);/s,
	'significant-phylogeny plotting is submitted with a durable checkpoint and adequate memory');
like($strain2,
	qr/\$phyloFigureCheckpoint = "\$resultsDir\/phyloFigures\.sto";.*?unlink \$phyloFigureCheckpoint.*?if -e \$phyloFigureCheckpoint/s,
	'explicit rewrites invalidate the submitted phylogeny-figure checkpoint');

like($strain2,
	qr/my \$cmdPrelude = "set -eo pipefail\\nulimit -s 20000 2>\/dev\/null \|\| true\\n"\s*\."MF_ANALYSIS_FAILURES=0\\n";/s,
	'generated R-analysis scripts keep errexit but drop nounset and treat the stack limit as advisory');
unlike($strain2, qr/set -euo pipefail/,
	'nounset never reaches the generated scripts, whose conda activation sources hooks that are not nounset-clean');
like($strain2,
	qr/my \$isolatedAnalysisBlock = sub \{.*?"if \(\\n" \. \$body.*?"test -s " \. shellQuote\(\$store\).*?"\); then\\n".*?MF_ANALYSIS_FAILURES=\\\$\(\(MF_ANALYSIS_FAILURES \+ 1\)\)/s,
	'each MGS analysis runs in a tested subshell, so one failure cannot abort the rest of its batch');
like($strain2,
	qr/my \$cmdEpilogue = .*?MF_ANALYSIS_FAILURES.*?-gt 0.*?exit 1/s,
	'a batch that attempted every MGS still reports failure so Slurm OOM escalation keeps working');
like($strain2,
	qr/\$batchCmd \.= \$cmdEpilogue;\s*my \(\$dep, \$qcmd\) = qsubSystem\(.*?command => \$batchCmd,/s,
	'the retried command includes the same failure epilogue as the original submission');
my %isolatedAnalysisCommand = (
	strainStats => 'strainStatsR', popGenStats => 'popGenCommand',
);
for my $analysis (sort keys %isolatedAnalysisCommand) {
	my $variable = $isolatedAnalysisCommand{$analysis};
	like($strain2,
		qr/\$cmd \.= \$isolatedAnalysisBlock->\(\s*"\$\Q$variable\E\b/,
		"$analysis analyses are individually isolated within their batch");
}
unlike($strain2, qr/MATAFILER_R_ANALYSIS_CORES/,
	'R and PopGenStats commands use the resolved allocated core count rather than an unexpanded shell variable');
unlike($strain2, qr/\$batchSize/,
	'R-job submission no longer uses a fixed phylogeny count per batch');
unlike($strain2, qr/if \(\$doSubmit && -d \$destD\)/,
	'partial result recovery does not erase an entire within directory outside an explicit rewrite');
like($strain2,
	qr/my \$networkDir = "\$resultsDir\/networks";.*?remove_tree\(\$networkDir\) if -d \$networkDir/s,
	'explicit rewrites clear the workflow-owned network cache');
like($strain2,
	qr/sub strainNetwork.*?my \$netDir = "\$resultsDir\/networks\/";.*?-o .*?shellQuote\(\$netDir\)/s,
	'network tables and PDFs are published below Results');
like($strain2,
	qr/sub treeWas.*?my \$treewasOut = "\$resultsDir\/GeneEnrich\/";.*?-o .*?shellQuote\(\$treewasOut\)/s,
	'treeWAS tables are published below Results');
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
like($strain,
	qr/\$threadMemFactor = \$numCoreL \/ \$treeMemThreadDivisor;\s*\$threadMemFactor = 1 if \$threadMemFactor < 1;\s*\$totMem = int\(\$memoryPlanningInputMB \* \$baseMemMult \* \$memMulti \* \$threadMemFactor\)/s,
	'the initial tree memory request scales with the thread count IQ-TREE will use');
like($strain,
	qr/my \$treeOOMRetryRounds = 8;/,
	'the OOM retry round default is high enough to reach the configured memory ceiling');
like($strain,
	qr/die "OOM retry rounds must be between 0 and 12\\n"/,
	'the retry-round guard no longer caps recovery far below the memory ceiling');
like($strain,
	qr/my \$retryIqMemMB = int\(\$nextMB \* 0\.9\);\s*\$retry\{command\} =~ s\/\(\^\|\\s\)-iqMemMB\\s\+\\d\+\/\$1-iqMemMB \$retryIqMemMB\//s,
	'an escalated OOM retry updates the memory allowance carried in its saved command');
like($strain,
	qr/"treeOOMRetryRounds=i" => \\\$treeOOMRetryRounds,.*?"treeMemThreadDivisor=f" => \\\$treeMemThreadDivisor,/s,
	'both tree-memory controls are tunable from the command line');
like($build_tree_locus_module,
	qr/my \$sample_set_seeds = \$options->\{sample_set_seeds\};/,
	'only the sample sets may be restricted to the seeds that can merge');
like($build_tree_locus_module,
	qr/\$sample_set->\{\$gene\}\{\$sample\} = 1 if \$wantSampleSet;.*?\$gene_context->\{\$entries\[\$i\]\[1\]\}\{\$token\}\+\+ if \$include_gene_context;/s,
	'gene contexts retain the same directed neighbour counts for each focal seed');
like($build_tree_locus_module,
	qr/my \$context_seeds = \$options->\{context_seeds\}.*?my %relevant_contigs;.*?\$relevant_contigs\{\$sample\}\{\$contig\} = 1.*?next if defined\(\$context_seeds\).*?!\$context_seeds->\{\$entries\[\$i\]\[1\]\}/s,
	'pre-budgeting restricts focal contexts while retaining every ranked neighbour on their contigs');
like($strain,
	qr/# The scan always runs, even with no mosaic pair to merge.*?my \$scanScope = catalogueLocusContext\(/s,
	'the catalogue-wide scan is unconditional, since locus context is needed with or without a merge');
like($build_tree_locus_module,
	qr/s\/\^>\/\/ if substr\(\$_, 0, 1\) eq '>';\s*s\/\^\\s\+\|\\s\+\$\/\/g if \/\\s\/;/s,
	'per-member cleaning is guarded, since it runs once for every catalogue member');
like($strain,
	qr/sub fastRemoveTree \{.*?my \$victim = rename\(\$target, \$parked\) \? \$parked : \$target;.*?bsd_glob\("\$target\.deleting\.\*"\).*?system\('sh', '-c', 'rm -rf -- "\$1" &', 'sh', \$path\)/s,
	'large trees are freed by one rename and unlinked by a backgrounded rm, with leftovers swept');
like($strain,
	qr/sub fastRemoveTree \{.*?\$systemRemoveAvailable = \$\^O eq 'MSWin32' \? 0.*?my \$failed = !eval \{ remove_tree\(\$path\); 1 \};/s,
	'a platform without a usable rm still removes the tree in-process');
for my $bigTree ('$outD', '$scratchD', '$preConDir', '$outD2', '$scratch_mgs') {
	like($strain, qr/fastRemoveTree\(\Q$bigTree\E[,)]/,
		"the $bigTree tree is removed through the fast path");
}
unlike($strain, qr/remove_tree\(\$outD\)|remove_tree\(\$scratchD\)/,
	'initialization no longer walks the output or scratch trees from Perl');
like($strain, qr/remove_tree\(\$locSpace\) if -d \$locSpace;/,
	'small per-sample temporaries stay in-process, where forking rm would cost more than it saves');
like($strain, qr/my \$version = 1\.64;/,
	'workflow behavior changes retain an explicit version marker');
like($strain,
	qr/my \$resumeOutD = .*?my \$parentRunLock;.*?if \(!\$subJob\).*?\$parentRunLockPath = "\$lockBase\.strain_within\.lock".*?acquire_workflow_lock\(.*?prepRun\(\);/s,
	'the parent acquires a stable sibling lock before prepRun can remove output or scratch state');
like($strain,
	qr/my \$SNPcaller = "MPI";.*?"SNPcaller=s"\s*=> .*?SNPcaller.*?-SNPcaller must be MPI or FB.*?genes\.shrtHD\.SNPc\.\$\{SNPcaller\}\.fna\.gz.*?allSNP\.\$\{SNPcaller\}\.vcf\.gz/s,
	'strain Phase I selects MPI or FB caller-specific consensus and VCF filenames');
like($strain,
	qr/my \$legacyMPIContract = \$phase1ContractState eq 'missing' && \$SNPcaller eq 'MPI'.*?if \(\$onlySubmit && !\$subJob.*?'building_match'.*?Cannot reuse Phase-I outputs.*?if \(\$onlySubmit && \$subJob.*?'building_match'.*?!\$legacyMPIContract.*?Split worker Phase-I input contract/s,
	'caller provenance resumes compatible builds and rejects missing FB or incompatible state');
like($strain,
	qr/my \$resumeOutD = .*?\$phase1MGSGuideFingerprint =\s*phase1GuideStatFingerprint\(\$MGSfileOri, undef, \$resumeOutD\);.*?phase1InputContractState\(\$resumePhaseIInputContract\)/s,
	'workers fingerprint the run-local guide before validating the parent contract');
like($strain,
	qr/my \$phase1InputContractVersion = 3;.*?sub phase1PathStatComponent.*?different device number.*?\$contractVersion <= 2.*?\@metadata\[0, 1, 7, 9\].*?\@metadata\[1, 7, 9\]/s,
	'Phase-I provenance excludes node-local device identity from current cross-node contracts');
like($strain,
	qr/sub phase1GuideStatFingerprint.*?strain-phase1-guide-stat-v2.*?strain-phase1-guide-stat-v3.*?\$staged\.srt.*?\$sorted\.gene2MGS.*?\$observation/s,
	'Phase-I provenance tracks the original guide plus the run-local sorted, indexed and observation inputs');
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
	qr/my \$rmMSA = 1;.*?my \$doPopGenStats = 0;.*?my \$popGenStrictOutgroup = 0;.*?my \$popGenGeneticCode = 1;.*?my \$popGenCodonStart = 1;.*?my \$popGenSeed = 1;.*?"popGenStats=i"\s*=> \\\$doPopGenStats.*?"popGenStrictOutgroup=i".*?"popGenCategory=s"\s*=> \\\$popGenCategory.*?"individualVar=s".*?if \(!\$subJob && \$doPopGenStats && \$rmMSA\).*?\$rmMSA = 0;.*?-rmMSA \$rmMSA.*?-popGenStats \$doPopGenStats.*?-popGenStrictOutgroup \$popGenStrictOutgroup.*?-popGenGeneticCode \$popGenGeneticCode.*?-popGenCodonStart \$popGenCodonStart.*?-popGenSeed \$popGenSeed.*?-popGenLegacyTextOutput \$popGenLegacyTextOutput.*?-popGenCategory ".*?shellQuote\(\$popGenCategory\).*?if \$popGenCategory ne ""/s,
	'population genetics is opt-in, retains MSAs when enabled, and forwards reproducible configuration and categories through strainwithin2');
like($strain, qr/Retain the Phase-I locus map.*?second catalogue-wide gene2tax scan.*?\$SIgenes and \$COGprios are reused/s,
	'Phase II reuses the Phase-I selected gene map rather than clearing and rebuilding it');
like($strain,
	qr/my %broadCOG = map.*?\$cogTaxa\{\$_\} >= \$broadMinimumTaxa.*?exists\(\$preferredCoreGeneSet->\{\$gene\}\).*?unless \(%outgroupCatalogueMGS\).*?return;.*?Preparing core-first exact outgroup-reference demands.*?readFasta\(\$refFAA, 1, "\\\\s", \\%requiredAA,.*?readFasta\(\$refFNA, 1, "\\\\s", \\%requiredNT,/s,
	'outgroup viability is decided before candidate-map and exact FNA/FAA reference loading');
like($strain,
	qr/my \$outgroupCoreMinLoci = 0;.*?"outgroupCoreMinLoci=i".*?\$outgroupCoreMinLoci = int\(\$treeLocusBudget \* 0\.20 \+ 0\.999999\).*?if \$outgroupCoreMinLoci == 0;.*?\$minimumOutgroupLoci = \$outgroupDemandMinimum\{\$MGS\} \/\/ \$minLociPerMGS/s,
	'the outgroup floor defaults to 20% of the final-tree locus budget and is enforced per MGS');
like($strain,
	qr/my \$outgroupReferenceGeneCap = 2500;.*?readGene2tax\(.*?\$outgroupReferenceGeneCap.*?allowed_cogs_by_mgs => \\%eligibleCogsByOutgroup.*?my \$addCandidate = sub.*?return if \$retainedForMGS >= \$outgroupReferenceGeneCap/s,
	'candidate maps retain only eligible COGs for viable selected outgroups before applying the 2,500-gene cap');
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
	qr/locus-model catalogue-protein loading.*?loaded_proteins=.*?locus-model cluster-index scan.*?represented_seed_clusters=.*?locus-group construction.*?ranked_clusters=.*?resolved_loci=.*?worker-member materialization.*?locus_sample_combinations=/s,
	'locus-model construction reports protein, scan, grouping, and worker-projection subphase timings');
like($strain,
	qr/if \(\$maxSubJob > 1\).*?phase1LocusModelFingerprint.*?loadPhase1LocusModel.*?merge_candidate_seeds.*?catalogueLocusContext.*?publishPhase1LocusModel/s,
	'a split run reuses a published locus model and otherwise rebuilds it from a candidate-restricted catalogue scan');
like($strain,
	qr/\$digest->add\('strain-phase1-locus-model-v1', "\\0", \$version, "\\0",/,
	'a published locus model is invalidated by a workflow version change, not only by its inputs');
like($strain,
	qr/sub catalogueLocusContext \{.*?sample_set_seeds => \$mergeCandidates.*?readClstrRevBinaryShard\(\$phase1ShardPaths->\[\$worker\].*?accumulate_locus_context\(\$accumulator.*?\$clusters = \{\};.*?return "streamed_/s,
	'the catalogue-wide scan streams one published shard at a time instead of holding every sample');
like($strain,
	qr/sub catalogueLocusContext \{.*?my \$scanStarted = \$scanOptions->\{started\} \/\/ time.*?stepProgress\("catalogue-wide locus-model scan".*?\$scanStarted/s,
	'catalogue scan progress uses the scan start rather than global process start');
unlike($strain,
	qr/stepProgress\("catalogue-wide locus-model scan", \$scanned, \$maxSubJob,\s*\$\^T/,
	'catalogue scan progress no longer labels global runtime as step elapsed');
like($strain,
	qr/sub catalogueLocusContext \{.*?my \(undef, \$fullClusters\) = readClstrRev\(\$clusterIndex, 0, \$Gene2COG\);.*?return .whole_catalogue_read./s,
	'a missing or unreadable shard set still falls back to one whole-catalogue read');
like($strain,
	qr/loadPhase1ClusterIndex\(\s*\$cluster_index, \$workerForSampleHR, \$mySamplesHR\);.*?if \(\$maxSubJob > 1\) \{\s*\$modelFingerprint/s,
	'the worker slice is loaded before the model scan, so the published shards exist to stream');
like($strain,
	qr/preselect_locus_records\(\\\@records, \$treeLocusBudget.*?context_seeds => \\%modelContextSeeds.*?buildSelectedLocusGroups\(\$modelRecordsRef, \{\}, \{.*?precomputed_context => \\%catalogueContext.*?prebudget_excluded => \$prebudgetExcluded/s,
	'non-taxon-aware grouping pre-budgets focal contexts but retains merge backfill and streamed summaries');
like($strain,
	qr/sub buildSelectedLocusGroups.*?precomputed_context => \$options->\{precomputed_context\}.*?Catalogue-wide locus context was computed but grouping retained no context/s,
	'the workflow forwards precomputed context and refuses to publish a silently empty handoff');
like($build_tree_locus_module,
	qr/if \(my \$precomputed = \$options->\{precomputed_context\}\).*?\$sample_set = \{ %\{\$precomputed->\{sample_set\} \|\| \{\}\} \}/s,
	'precomputed summaries are shallow-copied so grouping cannot consume the caller index');
like($strain,
	qr/my %selectedContextSeeds.*?member_context_map\(\\\@records, \$cl2gene,\s*\{ context_seeds => \\%selectedContextSeeds \}\)/s,
	'split workers derive focal member contexts while retaining all ranked neighbours');
like($strain,
	qr/sub publishPhase1LocusModel.*?retry_open\('>', \$groupTemporary.*?print \{\$groupOut\}.*?retry_open\('>', \$contextTemporary.*?print \{\$contextOut\}.*?contextRows\+\+/s,
	'the corrected nonempty common model is streamed instead of duplicated in row arrays and strings');
like($strain,
	qr/stepComplete\("locus-group construction".*?model_source=\$modelSource/s,
	'the locus-model source is reported so divergent worker models are visible in the logs');
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
	qr/my \@pendingTreeJobs;.*?push \@pendingTreeJobs, \{.*?command => \$Tcmd\.\$outgS.*?tmp_space => \$QSBoptHR->\{tmpSpace\}.*?dispatchPendingTreeJobs\(.*?blocking => 0.*?Tree preparation pass complete:.*?dispatchPendingTreeJobs\(.*?blocking => 1.*?retryOOMTreeJobs\(\s*jobs => \\\@jobs,.*?writeTreeFailureAudit.*?without a valid output were quarantined/s,
	'eligible trees queue after conversion, then are globally submitted, tracked, awaited, and output-validated');
like($strain,
	qr/if \(!\$doSubmit \|\| \(\$epaOnlyRetry && time >= \$nextQueuedTreeSubmissionProbe\)\).*?Tree preparation pass complete:.*?blocking => 1/s,
	'ordinary Phase II jobs remain queued for global priority, while EPA-only recovery can still dispatch promptly');
like($strain,
	qr/\@\{\$queue\} = sort \{.*?epa_only.*?cores.*?workload_cells.*?sample_count.*?requested_mb.*?gene_count.*?priority_ordinal.*?\} \@\{\$queue\}/s,
	'Phase II submission retains EPA recovery priority, then orders full jobs by cores and the shared approximate workload');
like($strain,
	qr/push \@pendingTreeJobs, \{.*?cores => \$numCoreL.*?sample_count =>.*?workload_cells => \$workloadCells.*?gene_count =>.*?requested_mb => int\(\$totMem\).*?priority_ordinal => \$treeJobOrdinal/s,
	'each queued job retains the size signals used for cores, memory, and stable ordering');
like($strain, qr/nonblockingMaxConcurrentJobs\} = 1 unless \$blocking/,
	'queued tree dispatch uses the non-blocking scheduler-capacity path');
like($strain, qr/deferredSubmissionDependency\(\).*?retaining "\.scalar\(\@\{\$queue\}\).*?until live jobs fall below/s,
	'queued tree dispatch retains deferred jobs until scheduler capacity frees');
like($strain, qr/\$options->\{tmpSpace\} = \$record->\{tmp_space\}/,
	'queued tree dispatch restores each job\'s stored temporary-space setting');
like($strain, qr/\$options->\{useLongQueue\} = \$record->\{use_long_queue\}/,
	'queued tree dispatch restores each job\'s stored queue setting');
like($strain,
	qr/\@treeJobAccounting.*?requested_mb => int\(\$totMem\).*?qsubSystemJobAlive.*?slurm_oom_retry_plan.*?format_slurm_tree_memory_summary/s,
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
	&& index($strain, q{return (scalar(keys %samplesSeen), $genesSeen, $preparedOG, 1, 1, $ingroupSeen);}) >= 0,
	'legacy fully prepared Phase-II scratch inputs remain resumable without redoing their controller-side work');
ok(index($strain, 'sub preparedOutgroupLog') >= 0
	&& index($strain, 'fileGZe($log_path)') >= 0
	&& index($strain, '$publishedPrepared') >= 0
	&& index($strain, '$scratchPrepared') >= 0
	&& index($strain, '.strain_tree_input.plan.tsv') >= 0
	&& index($strain, '.strain_tree_input.shards.tsv') >= 0
	&& index($strain, 'writeMGSShardManifest') >= 0
	&& index($strain, '; outgroup ') >= 0
	&& index($strain, q{; $ingroupSmpl ingroup samples ($multiSmpl tips); $ngenes genes; }) >= 0,
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
	/(sub validatePhase1WorkerLedger .*?)(?=^sub )/ms;
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
	qr/\$ingroupSmpl > 2 && \$ngenes >= \$minLociPerMGS.*?too_few_usable_genes.*?writeTooFewMarker.*?sub validateTreeInputResolution.*?tree_input_resolution\.tsv.*?repair_required.*?tree_input_repair\.queue\.tsv.*?no catalogue-wide abort was triggered/s,
	'insufficient tree inputs are terminally marked while incomplete triplets enter a persistent repair queue');
like($strain,
	qr/\$minimumOutgroupLoci = \$outgroupDemandMinimum\{\$MGS\} \/\/ \$minLociPerMGS.*?my \@requiredLoci = sort grep.*?\$OG = \$SelectedOutgroup\{\$MGS\}.*?\@requiredLoci < \$minimumOutgroupLoci.*?\$represented < \$minimumOutgroupLoci/s,
	'the predetermined outgroup is checked against the per-MGS core/broad demand floor');
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
	qr/my \$treeTmpGb = int\(\(\$inputFNAsize \* 5 \+ 1023\) \/ 1024\);.*?\$treeTmpGb = 20 if \$treeTmpGb < 20/s,
	'strain BuildTree jobs reserve five times compressed input size with a 20 GiB scratch floor');
like($strain,
	qr/addOutgroup2MGS\(\$MGS,\$OG,\$tmpD\).*?choose_tree_core_count\(\$multiSmpl, \$maxCores\).*?\$Tcmd \.= "-cores \$numCoreL "/s,
	'BuildTree core guidance reuses the exact submitted-sample count collected during input preparation');
like($strain,
	qr/addOutgroup2MGS\(\$MGS,\$OG,\$tmpD\).*?\$workloadCells = \$multiSmpl \* \$ngenes.*?\$taxonLocusInputMB = \$workloadCells \/ 1024.*?\$memoryPlanningInputMB = \$taxonLocusInputMB > \$inputFNAsize.*?\$totMem = int\(\$memoryPlanningInputMB \* \$baseMemMult \* \$memMulti \* \$threadMemFactor\)/s,
	'Phase II memory mixes the already collected sample and usable-gene counts without a second sizing pass');
like($strain,
	qr/my \$treeFlag = \$onlyMSA.*?my \$Tcmd=.*?\$treeFlag /s,
	'BuildTree command retains the selected inference-engine flag while appending cores after input preparation');
unlike($strain, qr/largestFullTreeInput|sqrt\(\$inputFNAsize/,
	'BuildTree core planning no longer uses input byte size');
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
	qr/sub resetMGSTreeOutputs .*?dirname\(\$resolvedMGS\) eq \$resolvedRoot.*?basename\(\$resolvedMGS\) eq \$MGS.*?qw\(phylo MSA within\).*?next unless -d \$directory;.*?fastRemoveTree\(\$directory\).*?retry_unlink\(\$treeStone/s,
	'tree-only reset is confined to selected existing MGS tree-stage directories and clears them asynchronously');
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
	qr/my \$GenesPerSpecies = 0\.2;.*?my \$GeneLengthMin = 0\.4;.*?my \$GeneLengthIncludeMin = \$GeneLengthMin;.*?my \$geneLengthIncludeMinSpecified = 0;.*?my \$relativeNTFraction = 0\.1;.*?\$placementGenesPerSpecies = 0.04; \$placementRelativeNTFraction = 0.03;.*?my \$taxonAwareLocusSelection = 0;.*?"GeneLengthIncludeMin=f" => sub \{.*?\$GeneLengthIncludeMin = \$_\[1\];.*?\$geneLengthIncludeMinSpecified = 1;.*?"taxonAwareLocusSelection=i" => \\\$taxonAwareLocusSelection.*?\$GeneLengthIncludeMin = \$GeneLengthMin unless \$geneLengthIncludeMinSpecified;.*?-GeneLengthIncludeMin \$GeneLengthIncludeMin.*?-taxonAwareLocusSelection \$taxonAwareLocusSelection/s,
	'strainWithin defaults to one 40% gene-length gate and hard locus filtering while retaining explicit overrides');
like($strain,
	qr/my \$postAlignmentSequenceOutlierMask = 1;.*?"postAlignmentSequenceOutlierMask=i" => \\\$postAlignmentSequenceOutlierMask.*?-postAlignmentSequenceOutlierMask ".*?\$postAlignmentSequenceOutlierMask/s,
	'strainWithin enables and forwards the native within-locus sequence-outlier masker');
like($strain,
	qr/minimum_informative_nt_per_sample => 5000,.*?my \$NTfiltCount = \$FILTER_DEFAULT\{minimum_informative_nt_per_sample\};/s,
	'strainWithin retains the 5 kb absolute informative-position floor');
like($strain,
	qr/multi_gene_sample_max => 0\.10,.*?my \$multiGeneSmplMax = \$FILTER_DEFAULT\{multi_gene_sample_max\};/s,
	'strainWithin limits unresolved multigene loci to 10% by default');
like($strain,
	qr/"GenesPerSpecies=f" => sub.*?\$genesPerSpeciesSpecified = 1.*?"relativeNTFraction=f" => sub.*?\$relativeNTFractionSpecified = 1.*?"placementGenesPerSpecies=f" => sub.*?\$placementGenesPerSpeciesSpecified = 1.*?"placementRelativeNTFraction=f" => sub.*?\$placementRelativeNTFractionSpecified = 1.*?resolvePairedOptionDefault/s,
	'strainWithin resolves each coverage pair only after recording explicit options');
like($build_tree,
	qr/"relativeNTFraction=f" => sub.*?\$ntFracSpecified = 1.*?"GenesPerSpecies=f" => sub.*?\$geneFracPSpecSpecified = 1.*?"placementGenesPerSpecies=f" => sub.*?\$placementGeneFracPSpecSpecified = 1.*?"placementRelativeNTFraction=f" => sub.*?\$placementNTFracSpecified = 1.*?resolvePairedOptionDefault/s,
	'BuildTree applies the same explicit-option-aware coverage-pair defaults');
like($build_tree,
	qr/my \$ntFrac =0\.2; my \$ntFracGene = 0\.4;.*?my \$geneLengthIncludeMinSpecified = 0;.*?my \$ntFracGeneInclude = \$ntFracGene;.*?my \$fracMaxGenes90pct = 0\.3;.*?"GeneLengthIncludeMin=f" => sub \{.*?\$ntFracGeneInclude = \$_\[1\];.*?\$geneLengthIncludeMinSpecified = 1;.*?\$ntFracGeneInclude = \$ntFracGene unless \$geneLengthIncludeMinSpecified;/s,
	'BuildTree uses one 40% gene-length gate by default and restores the 30%-of-Q90 hard locus filter');
like($build_tree,
	qr/my %TAXON_AWARE_DEFAULT = \(\s*enabled => 0,.*?my \$taxonAwareHeterogeneityScoring = 1;/s,
	'BuildTree disables taxon-aware selection by default but restores its heterogeneity score when enabled');
like($build_tree,
	qr/between_species_sequence_outlier_mask => 0,.*?within_species_sequence_outlier_mask => 1,.*?"postAlignmentSequenceOutlierMask=i" => \\\$postAlignmentSequenceOutlierMask.*?-maskSequenceOutliers.*?-sequenceOutlierReport.*?-sequenceOutlierExemptPrefix/s,
	'BuildTree mode defaults and MSAfix invocation keep masking strain-specific and exempt the outgroup');
like($build_tree,
	qr/sub mergeMSAs.*?informativeSequenceLength\(\$bigMSAFAA\{\$kk\}, \$isAA\).*?classifyTaxonAwareCoverageEligibility\(.*?minimum_nt => \$ntCntTotal/s,
	'BuildTree rechecks sample coverage from the final overlap-filtered concatenation without another alignment scan');
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
	qr/my \$backboneEligibility = classifyTaxonAwareCoverageEligibility.*?if \(!\$strictBackbone && \$enforceSampleCoverage\) \{.*?delete \$samples\{\$_\} for \@coverageRemoved;.*?if \(\$strictBackbone\) \{\s*my \$placementEligibility = classifyTaxonAwareCoverageEligibility/s,
	'BuildTree removes samples below the coverage thresholds when placement is disabled, and routes them to placement when it is not');
like($build_tree,
	qr/my \$excludeFlaggedSamples = 0;.*?my \$enforceSampleCoverage = 0;/s,
	'both optional sample filters are generic mechanisms that stay off unless a caller sets a policy');
like($strain,
	qr/\$Tcmd \.= "-excludeFlaggedSamples \$excludeMixedStrainSamples ";\s*\$Tcmd \.= "-enforceSampleCoverage \$enforceSampleCoverage ";/s,
	'the strain workflow is the caller that turns both BuildTree sample filters on');
like($strain,
	qr/my \$minLociPerMGS = \$FILTER_DEFAULT\{minimum_loci_per_mgs\};/,
	'the per-MGS locus floor is separate from the per-sample extraction prefilter');
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
	qr/waitPhase1WorkersWithOOMScan\(\s*jobs => \\\@jobsMain.*?my \@failedWorkers = phase1WorkersNeedingRetry\(\$splitGeneration\).*?retryPhase1Workers\(.*?script_kind => "retry"/s,
	'fresh Phase-I validation supervises the live wave and routes failed workers through the standard retry path');
like($strain,
	qr/sub recoverCompletedSplitPhaseI.*?phase1WorkersNeedingRetry\(\$generation\).*?retryPhase1Workers\(.*?script_kind => "resume"/s,
	'resumed Phase-I repair routes failed workers through the same retry path');
like($strain,
	qr/sub retryPhase1Workers.*?slurm_oom_retry_plan.*?ceiling_reached.*?OOM escalation.*?memory_mb => \$retryMemoryMB\{\$worker\}/s,
	'the shared Phase-I retry path escalates only accounting-confirmed OOM memory and honors the ceiling');
like($strain,
	qr/sub submitPhase1Worker.*?\$memoryMB\."M", "Str1\.\$worker"/s,
	'one place resubmits a Phase-I worker for both the live OOM scan and the ledger audit');
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
	'automatic Phase-I and tree OOM recovery use maxMF4mem while retaining an explicit per-run override');
like($site_config_template, qr/^maxMF4mem\t512$/m,
	'the installed site-config template sets the shared strain OOM ceiling to 512 GiB');
unlike($internal_config, qr/^downloadQueue\t/m,
	'the internal program config leaves the archive queue to the site config');
like($strain,
	qr/sub retryOOMTreeJobs.*?qsubSystemJobAlive\(\s*\\\@pendingJobs, \$options, 0, -1, \$budget\).*?slurm_oom_retry_plan.*?by_job_id.*?next_mb.*?epaOnlyRetryReady\(\$mgsDirectory, 1\).*?-epaThreads.*?epaThreads 1.*?treeCmd\.epa_retry\.sh/s,
	'tree retries rescan the shared accounting-confirmed OOM plan on a bounded wait, and EPA-stage retries use one thread');
like($strain,
	qr/\(\$b->\{oom_retry\} \/\/ 0\) <=> \(\$a->\{oom_retry\} \/\/ 0\)\s*\|\|\s*\(\$b->\{epa_only\} \/\/ 0\) <=> \(\$a->\{epa_only\} \/\/ 0\)/s,
	'OOM retries take the leading dispatch tier, ahead of EPA recovery and bulk work');
like($strain,
	qr/\$retry\{oom_retry\} = 1;\s*\$retry\{job_nice\} = 0;/s,
	'an escalated tree retry drops the bulk priority handicap');
like($strain,
	qr/unshift \@\{\$pendingQueue\}, \@retryQueue;/,
	'escalations are injected ahead of the ordinary jobs still awaiting capacity');
like($strain,
	qr/my \$jobNice = 2500;.*?my \$maxQueuedJobs = 0;/s,
	'bulk submissions carry a priority handicap and an optional queue ceiling');
like($strain,
	qr/\$QSBoptHR->\{jobNice\} = \$jobNice;\s*\$QSBoptHR->\{maxConcurrentJobs\} = \$maxQueuedJobs;/s,
	'the scheduler options carry the run-wide priority and capacity settings');
like($strain,
	qr/my \$round = \(\$retriesByMGS\{\$original->\{mgs\}\} \|\| 0\) \+ 1;.*?if \(\$round > \$maximumRounds\).*?exhausted for/s,
	'the tree OOM retry budget is spent per MGS, not per submission wave');
like($strain,
	qr/my \$oomScanMinutes = 60;.*?my \$oomMinRetries = 3;/s,
	'OOM outcomes are rescanned hourly with a minimum per-job retry contract');
like($strain,
	qr/"oomScanMinutes=f" => \\\$oomScanMinutes,\s*"oomMinRetries=i" => \\\$oomMinRetries,/s,
	'both OOM rescan controls are tunable from the command line');
like($strain,
	qr/\$maximumRounds = \$minimumRounds if \$maximumRounds < \$minimumRounds;/,
	'the per-MGS tree retry budget never falls below the minimum OOM contract');
like($strain,
	qr/sub phase1RetryBudget.*?\$oomConfirmed && \$oomMinRetries > \$phase1WorkerRetries\s*\? \$oomMinRetries : \$phase1WorkerRetries/s,
	'a confirmed Phase-I OOM keeps the minimum retry contract, ordinary failures the smaller budget');
like($strain,
	qr/sub waitPhase1WorkersWithOOMScan.*?qsubSystemJobAlive\(\s*\\\@pendingJobs, \$QSBoptHR, 0, -1, \$scanSeconds\).*?escalatePhase1WorkerOOM/s,
	'Phase-I waves are rescanned for OOM while the remaining workers keep running');
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
	qr/script => \$epaOnlyRetry \? "\$outD2\/treeCmd\.epa_retry\.sh".*?\$onlyMSA \? "\$outD2\/treeCmd\.msa_only\.sh".*?treeCmd\.msa_retain\.sh" : "\$outD2\/treeCmd\.sh".*?qsubSystem\(\$LOGDIR\."strainAnalysis2\.sh"/s,
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
	qr/"GeneLengthIncludeMin=f" => sub \{.*?\$geneLengthIncludeMinSpecified = 1;.*?geneLengthIncludeByGene.*?qualification_sequences => \\%geneLengthQCSequence/s,
	'buildTree5 supports explicit lower-threshold recovery while keeping final sample QC on high-threshold sequences');
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


# Stage-I worker memory. The output buffers hold one string per MGS the worker
# has touched, so a sample count cannot bound them; and once the locus model is
# built a worker consults neither the ranked seed index nor most catalogue
# proteins again, because it exits before Phase II.
like($strain,
	qr/my \$flushOutputMB = 2048;\s*my \$flushOutputByteLimit = \$flushOutputMB \* 1024 \* 1024;\s*my \$bufferedOutputBytes = 0;/s,
	'Stage-I output buffering has an explicit byte budget as well as a sample count');
like($strain,
	qr/\$bufferedOutputBytes \+= length\(\$aaChunk\) \+ length\(\$ntChunk\)\s*\+ length\(\$linkChunk\) \+ length\(\$catChunk\);/s,
	'every buffered record is counted towards that budget');
like($strain,
	qr/if \(\$\{\$bufferedSamplesRef\} >= \$appendWriteTrigger\s*\|\| \$bufferedOutputBytes >= \$flushOutputByteLimit\)/s,
	'reaching the byte budget flushes early, so peak memory does not scale with the MGS count');
like($strain,
	qr/sub appendWriteMGSgenes.*?\$bufferedOutputBytes = 0;/s,
	'the byte budget is reset when the buffers are published');
like($strain,
	qr/if \(\$subJob\) \{.*?for my \$locus \(keys \%contextLociNeeded\).*?\$keptProteins\{\$seed\} = \$catalogProteins->\{\$seed\}.*?\$catalogProteins = \\%keptProteins;.*?\$SIgenes = \{\}; \$Gene2COG = \{\};/s,
	'a worker keeps catalogue proteins only for loci that can be ambiguous, and drops the seed index');
is(scalar(() = $strain =~ /\$catalogProteins = \\\%keptProteins;/g), 1,
	'the protein release happens in exactly one place, the split-worker branch');


# Two strain_within runs over one gene catalogue must not share, overwrite or
# delete each other's derived MGS guide. The sorter names its output after the
# guide it is handed, so the guide is staged inside the output directory and
# every product follows it there.
like($strain,
	qr/my \$stagedGuide = File::Spec->catfile\(\$outD, basename\(\$MGSfileOri\)\);/,
	'the MGS guide is staged inside the run output directory');
like($strain,
	qr/my \$sortedGuide = \$guide;.*?File::Spec->catfile\(\$outputDirectory, basename\(\$guide\)\.'\.srt'\)/s,
	'a resume looks for the sorted guide in its own output directory');
like($strain,
	qr/\$observedName->\(\$MGSfileOri\),\s*\$observedName->\(\$stagedGuide\)/s,
	'the optional occurrence table is staged under the name the sorter derives, so prevalence is unchanged');
like($strain,
	qr/\$sortMGSgenes \. " "\s*\. join\(" ", map \{ shellQuote\(\$_\) \} \(\$GCd, \$stagedGuide,/s,
	'the sorter is run on the staged guide, so it writes its output per-run');
unlike($strain, qr/glob\("\$MGSfile\.srt\*"\)/,
	'no run deletes the catalogue-level guide products another run may be reading');
like($strain,
	qr/symlink\(\$stagedGuide, \$sortedMGS\)/,
	'MGSall mode links its sorted guide inside the output directory too');

done_testing();
