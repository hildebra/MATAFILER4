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
like($strain, qr/test -s "?\.shellQuote\(\$IQtreef\).*?touch "?\.shellQuote\(\$treeStone\)/s,
	'a tree completion stone is conditional on a nonempty tree');
like($strain, qr/if \(\$doSubmit\) \{.*?unlink \$treeStone.*?if -e \$treeStone/s,
	'a submitted tree retry cannot pass through a stale completion stone');
like($strain, qr/unlink \$IQtreef.*?stale tree output/s,
	'a submitted tree retry must publish a fresh nonempty tree');
like($strain, qr/clear_split_generation\(\$splitManifest.*?write_split_generation\(\$splitManifest.*?printf '%s\\\\n'/s,
	'a new split generation clears stale state and tags each worker completion');
like($strain, qr/ConspecificMGS\.\$subJob\.log.*?sub mergeConspecificLogs/s,
	'split workers write isolated conspecific logs that are explicitly merged');
like($strain, qr/\$onlySubmit == 0 && !\$subJob/,
	'split children cannot recursively clean shared MGS output directories');
like($strain,
	qr/\$publishedInputsReady = !exists\(\$legacyLocusMGS\{\$MGS\}\).*?fileGZe\("\$outD2\/\$FNAstdof"\).*?fileGZe\("\$outD2\/\$FAAstdof"\).*?fileGZe\("\$outD2\/\$CATstdof"\).*?if \(\$publishedInputsReady && !\$mustRegenerateInputs\).*?combineMGSgenesDir\(\$MGS,\$tmpD,\$tmpD\)/s,
	'complete published inputs bypass missing scratch aggregates during tree recovery');
like($strain, qr/has neither complete published inputs nor complete combined worker input/,
	'incomplete worker input is reported only when published recovery inputs are also incomplete');
like($strain, qr/hasFreshParts.*?lacks required.*?return 0/s,
	'fresh worker parts replace stale combined inputs only when every required part exists');
like($strain, qr/exact_worker_parts\(\$prefix, \$workerCount\).*?split_generation_complete\(\$splitManifest.*?return \$aggregateComplete/s,
	'partial retries and merge scratch files cannot replace a complete aggregate');
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
like($strain, qr/my \$version = 0\.55;/,
	'within-strain mosaic and placement QC increments the workflow version');
like($strain,
	qr/sub stepComplete .*?STEP COMPLETE: \$step/s,
	'step completion messages use one consistent formatter');
like($strain,
	qr/configuration and map initialization.*?assembly-group expansion.*?consensus-input audit.*?MGS and seed-locus selection.*?existing-output and resume audit/s,
	'startup stages emit consistent completion messages with elapsed time and statistics');
like($strain,
	qr/locus-model construction.*?catalogue_drivers=.*?resolved_loci=.*?consensus-gene extraction and publication.*?tree-input sizing/s,
	'major extraction and tree-preparation stages also report concise completion statistics');
like($strain,
	qr/historical exclusion loading.*?excluded_MGS=.*?outgroup-reference preparation.*?reference_NT=.*?MGS_with_outgroup_candidates=/s,
	'historical exclusions and outgroup-reference preparation report their final counts');
like($strain,
	qr/my %treeDisposition.*?\$treeDisposition\{'eligible tree job'\}\+\+.*?Tree submission accounting:.*?Tree submission pass complete:/s,
	'tree submission reports every eligible and skipped MGS disposition before waiting');
like($strain,
	qr/qsubSystem\(\$outD2\."treeCmd\.sh".*?\$cnt \+\+.*?push \(\@jobs,\$dep\).*?\$expectedTreeOutputs\{\$MGS\}.*?qsubSystemJobAlive\( \\\@jobs.*?Tree jobs completed without valid tree outputs/s,
	'every eligible tree is submitted, tracked, awaited, and output-validated');
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
	qr/my \$iqPathogen = 0.*?"iqPathogen=i"\s+=> \\\$iqPathogen.*?-iqPathogen', \$iqPathogen/s,
	'within-strain pathogen mode defaults off and is propagated explicitly to split workers');
like($strain,
	qr/"legacyMGTK=i"\s+=> \\\$legacyMGTK.*?-iqPathogen and -legacyMGTK are mutually exclusive.*?-legacyMGTK', \$legacyMGTK/s,
	'within-strain legacy IQ-TREE selection remains exclusive and is propagated to split workers');
like($strain,
	qr/my \$iqMemMB = int\(\$totMem \* 0\.9\).*?if \(\$legacyMGTK\).*?"-iqLegacy 1 ".*?"-iqMemMB \$iqMemMB ".*?"-iqPathogen 1 " if \$iqPathogen/s,
	'within-strain IQ-TREE defaults to bounded standard mode and enables CMAPLE only by explicit request');
like($strain,
	qr/"recalcTrees=i"\s+=> \\\$recalcTrees.*?-recalcTrees must be 0 or 1.*?\$onlySubmit = 1 if \$recalcTrees/s,
	'tree recalculation is validated and forced into input-recovery-only mode');
like($strain,
	qr/-recalcTrees cannot be combined with -repairCAT, -deepRepair, or -redoSubmissionData.*?-recalcTrees must be launched by the main strainWithin process/s,
	'tree recalculation rejects input-regeneration modes and split-worker execution');
like($strain,
	qr/if \(!\$recalcTrees && \(.*?Part I:: extracting relevant core MGS genes/s,
	'tree recalculation bypasses consensus and per-MGS input regeneration');
like($strain,
	qr/my \$publishedInputsReady = !exists\(\$legacyLocusMGS\{\$MGS\}\).*?if \(\$recalcTrees\).*?unless \(\$publishedInputsReady\).*?\$scratchInputsReady = combineMGSgenesDir\(\$MGS,\$tmpD,\$tmpD\).*?unless \(\$publishedInputsReady \|\| \$scratchInputsReady\).*?no recoverable inputs for recalculation.*?resetMGSTreeOutputs\(\$outD2, \$MGS\)/s,
	'tree outputs are reset only after complete published or recoverable staged per-MGS inputs are verified');
like($strain,
	qr/\$scratchInputsReady \|\|= combineMGSgenesDir\(\$MGS,\$tmpD,\$tmpD\).*?using complete staged FNA\/FAA\/category files.*?tree job will publish them to the MGS directory/s,
	'recalculated trees continue staged-input transformation and publication through the normal tree-job path');
like($strain,
	qr/staged input sets recovered for -recalcTrees: \$recalcScratchRecovered/,
	'tree submission accounting reports staged recalculation recovery separately from skipped dispositions');
like($strain,
	qr/sub resetMGSTreeOutputs .*?dirname\(\$resolvedMGS\) eq \$resolvedRoot.*?basename\(\$resolvedMGS\) eq \$MGS.*?remove_tree\(\$phyloDir, \{safe => 1\}\).*?unlink \$treeStone/s,
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
	qr/"flushEvery=i"\s+=> \\\$appendWriteTrigger.*?'-flushEvery', \$appendWriteTrigger.*?%outgroupGeneCache = \(\)/s,
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
	qr/readGenesSample_Singl\(\$sm, \$writeLink,\$sttime,\\\$appCnt\).*?\$\{\$bufferedSamplesRef\}\+\+.*?appendWriteMGSgenes\(\$writeLink\)/s,
	'expanded assembly-group output is flushed by sample to retain the RAM bound');
like($strain,
	qr/if \(\$mySamplesHR\).*?\$unrepresentedWorkerLoci\+\+.*?unless \$maxSubJob/s,
	'split-worker sparsity is summarized instead of reported as missing catalogue data');
like($strain,
	qr/buildTree5 validates its persistent checkpoints.*?my \$contPhylo = 1;.*?-continue \$contPhylo/s,
	'unfinished trees delegate checkpoint recovery to buildTree continue mode');
like($strain,
	qr/sub treeInputPrecopyCommand .*?if \( \$ready_test \).*?Using existing persistent tree inputs.*?staged_inputs=\(\).*?if \[\[ -d \$staging_q \]\].*?if \(\( \$\{#staged_inputs\[\@\]\} \)\).*?No usable staged tree inputs found.*?if ! \( \$ready_test \)/s,
	'persistent tree inputs take precedence and recovery uses staging only when required');
like($strain, qr/test -s .*?test -s .*?\.gz.*?tree inputs are incomplete in both staging and persistent storage/s,
	'tree recovery accepts compressed persistent inputs and reports incomplete recovery data');
unlike($strain, qr/\$pigzBin -p \$numCoreL \$tmpD\/\*/,
	'generated tree jobs no longer pass an unmatched scratch glob to pigz');

my $build_tree = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'phylo', 'buildTree5.pl'));
like($build_tree, qr/if \(\$numSeq < 3\)/,
	'three-sample MGS accepted by the wrapper are retained for a minimal tree');

done_testing();
