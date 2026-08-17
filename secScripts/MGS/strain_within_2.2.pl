#!/usr/bin/perl
#submits R scripts to work through phylos
#args: [GeneCat dir] [intra_phylo_dir] [abundance MGS] [mapping file[ [cores]

use warnings;
use strict;

use Getopt::Long qw( GetOptions );

use Mods::GenoMetaAss qw( readClstrRev systemW readMapS readFasta gzipopen);
use Mods::Subm qw(qsubSystem emptyQsubOpt qsubSystemJobAlive qsubSystemWaitMaxJobs slurmJobFailureSummary);
use Mods::SlurmAccounting qw(slurm_tree_memory_summary next_oom_retry_memory_mb);
use Mods::IO_Tamoc_progs qw(getProgPaths );
use Mods::geneCat qw(readGene2tax createGene2MGS);
use Mods::TamocFunc qw ( getFileStr );
use Mods::CatalogPaths qw(resolve_catalog_maps);
use File::Path qw(make_path remove_tree);
use File::Basename qw(dirname);
use Text::ParseWords qw(shellwords);

sub sumSummaries;
sub strainNetwork;
sub treeWas;
sub visualizeSignPhylos;
sub combineResults;
sub newickSampleCount;
sub newickNodeCount;
sub rAnalysisMemoryForCores;
sub increaseRAnalysisMemory;
sub reportRAnalysisSchedulerFailures;
sub resolveOutgroup;
sub loggedOutgroup;
sub savedCommandOutgroup;
sub mgsTreeOutgroup;

my $MGSTKdir = getProgPaths("MGSTKDir");
my $configuredMaxMF4mem = getProgPaths("maxMF4mem", 0);
my $rAnalysisOOMMaxMemoryGB = 512;
if (defined($configuredMaxMF4mem) && $configuredMaxMF4mem =~ /^([0-9]+(?:\.[0-9]+)?)$/ && $1 > 0) {
	$rAnalysisOOMMaxMemoryGB = $1 + 0;
} elsif (defined($configuredMaxMF4mem) && length($configuredMaxMF4mem)) {
	warn "Ignoring invalid maxMF4mem setting '$configuredMaxMF4mem'; using 512 GiB\n";
}
#.21: added $DiscTests $ContTests
#.22: R interface updated
#.23: Perl interface updated
#.24: added familyVar and groupStabilityVars arguments for stability calculations 
#.25: 3.4.24: added Kasia's scoary scripts
#.26: 11.7.24: replaced scoary scripts with treewas
#.27: 6.1.26: added Anthony's network scripts, more qsub processes
#.28: 2.3.26: adapted for different phylo names
#.29: validate inputs, honour dry-run mode, and repair outgroup/checkpoint handling
#.30: make tree selection, queue handling and R-analysis completion checks robust
#.31: keep dry runs non-destructive and wait for validated downstream analyses
#.32: preserve summary headers and require validated network and treeWAS outputs
#.33: summarize per-tree progress and make runtime configuration explicit
#.34: resolve catalog maps through LOGandSUB/inmap.txt
#.35: exclude MGS with terminal no-tree or retained placement-pending outcomes
#.36: prioritize durable completed-tree evidence during postprocessing discovery
#.37: recover missing outgroups from saved tree metadata or the source MGS tree
#.38: combine versioned strainStats result stores instead of concatenating reports
#.39: run and combine versioned population-genetics stores alongside strain summaries
#.40: use durable RDS outputs directly as per-MGS completion evidence
#.41: scale R-analysis memory with explicitly requested cores and report scheduler failures
#.42: forward reproducible population-genetics settings and retain strict shell failures
#.43: phase strainStats and PopGenStats independently and weight small-MGS R batches
#.44: submit significant-phylogeny figures with the strain-only postprocessing phase
#.45: retry incomplete R-analysis batches twice and double memory after Slurm OOMs
my $version = 0.45;

my $rewriteRanalysis = 0; my $doSubmit = 1;
my $checkMaxNumJobs = 400;
my $doPopGenStats = 1;
my $popGenSubsample = "10,20,30,100,200,500";
my $popGenStrictOutgroup = 0;
my $popGenGeneticCode = 1;
my $popGenCodonStart = 1;
my $popGenSeed = 1;
my $popGenLegacyTextOutput = 0;

my $DiscTests =""; my $ContTests = "";
my $familyVar = ""; my $groupStabilityVars = "";
my $GCd = "";#$ARGV[0];
my $nCore = 4;#$ARGV[4];
my $nCoreHeavy = 12;
# strainStats parallelizes within a single R process. Its working memory grows
# with that process's worker count, so a fixed request was prone to Slurm OOMs.
my $rAnalysisMemoryBaseGB = 24;
my $rAnalysisMemoryCoreThreshold = 4;
my $rAnalysisMemoryPerExtraCoreGB = 1;
# Failed R-analysis batches are retried twice; an accounting-confirmed OOM
# doubles that batch's Slurm memory request for its next attempt.
my $rAnalysisRetryRounds = 2;
my $refMap = "";#$ARGV[3];#"/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/maps/drama4.map";
#my $FMGpD = "$GCd/MGS/phylo";
my $FMGpD = "";#$ARGV[1];# if (@ARGV > 1);
my $abMatrix = "";#$ARGV[2];
#discrete and continous tests done on each MGS' strains.. based on column header in map file, comma separated
#$DiscTests = $ARGV[5]; $ContTests = $ARGV[6]; 
my $individualVar = "AssmblGrps";
my $qsubSystem = "";

my $MGSphylo = "";
GetOptions(
	"GCd=s"          => \$GCd,
	"map=s"          => \$refMap,
	"submit=i"       => \$doSubmit,
	"reSubmit=i"     => \$rewriteRanalysis,
	"popGenStats=i"  => \$doPopGenStats,
	"popGenSubsample=s" => \$popGenSubsample,
	"popGenStrictOutgroup=i" => \$popGenStrictOutgroup,
	"popGenGeneticCode=i" => \$popGenGeneticCode,
	"popGenCodonStart=i" => \$popGenCodonStart,
	"popGenSeed=i" => \$popGenSeed,
	"popGenLegacyTextOutput=i" => \$popGenLegacyTextOutput,
	"cores=i"        => \$nCore,
	"Hcores=i"        => \$nCoreHeavy, #for heavy jobs (like scoary)
	"FMGdir=s"     => \$FMGpD,
	"MGSmatrix=s"      => \$abMatrix,
	"DiscTests=s" => \$DiscTests, #discrete categoried (from map) to be tested via perMANOVA 
	"ContTests=s"     => \$ContTests, #continous variables (from map) to be tested for strain phylo signal
	"familyVar=s"      => \$familyVar, #column name in metadata containing family id
	"groupStabilityVars=s"      => \$groupStabilityVars, #column names of categories used for calculation of resilience and persistence
	"individualVar=s"      => \$individualVar, #column name specifying individual IDs, AssmblGrps by default
	"qsubSystem=s"    => \$qsubSystem, #optional queue backend override; dry runs default to bash
	"MGSphylo=s"       => \$MGSphylo, #source MGS tree used by strain_within to choose fallback outgroups
) or die "Invalid strain_within_2.2.pl options\n";
die "Unexpected positional arguments: @ARGV\n" if @ARGV;
die "-GCd, -map, -FMGdir and -MGSmatrix are required\n"
	unless length($GCd) && length($refMap) && length($FMGpD) && length($abMatrix);
die "Gene-catalog and phylogeny directories must exist\n" unless -d $GCd && -d $FMGpD;
die "Map or MGS abundance matrix is missing\n" unless -s $refMap && -s $abMatrix;
die "Core requests must be positive and R-analysis retry rounds must be at least two\n"
	unless $nCore > 0 && $nCoreHeavy > 0 && $rAnalysisRetryRounds >= 2;
die "-submit and -reSubmit must be 0 or 1\n" unless $doSubmit =~ /^[01]$/ && $rewriteRanalysis =~ /^[01]$/;
die "-popGenStats must be 0 or 1\n" unless $doPopGenStats =~ /^[01]$/;
die "-popGenStrictOutgroup and -popGenLegacyTextOutput must be 0 or 1\n"
	unless $popGenStrictOutgroup =~ /^[01]$/ && $popGenLegacyTextOutput =~ /^[01]$/;
die "-popGenGeneticCode must be positive, -popGenCodonStart must be 1, 2, or 3, and -popGenSeed must be non-negative\n"
	unless $popGenGeneticCode > 0 && $popGenCodonStart >= 1 && $popGenCodonStart <= 3
		&& $popGenSeed >= 0;
$FMGpD =~ s{/+$}{} unless $FMGpD eq "/";
die "-MGSphylo does not exist or is empty: $MGSphylo\n"
	if length($MGSphylo) && !-s $MGSphylo;


my $cpDir = "";#$GCd/MGS/R_analysis/";
#my $outD
#my $defTreeFile = "IQtree_allsites.treefile";
my @defTreeFiles = ("IQtree_allsites.treefile","VERYFASTTREE_allsites.nwk","FASTTREE_allsites.nwk"); #multiple options to test for..
my $defTreeFileBase = "IQtree_allsites";
#die;
die "MGS phylo dir doesn't exist!\n$FMGpD\n" unless (-d $FMGpD);

my $queueMode = $qsubSystem;
$queueMode = "bash" if !$doSubmit && $queueMode eq "";
my $QSBoptHR = emptyQsubOpt($doSubmit,"",$queueMode);
$QSBoptHR->{tmpSpace} = 0;
my $bts = getProgPaths("buildTree_scr");
my $neighborTree = length($MGSphylo) ? getProgPaths("neighborTree") : "";
my $SaSe = "|";
my $SaSe2 = "\|"; #for regex
my $numCores = 40; #used for phylos..
my $RsummaryTab = "$FMGpD/strainStats.tsv";
my $popGenSummaryTab = "$FMGpD/popGenStats.tsv";
my $popGenSubsampleSummaryTab = "$FMGpD/popGenStats.subsamples.tsv";
my $legacyRsummaryTab = "$FMGpD/Rsummary.tab";


my %dirs;my %destDs; my %baseD;
my @nonTreeOutcomeMarkers = qw(
	tooFewSamples.sto noRecoverableLoci.sto noTree.sto placementPending.sto
);
my $terminalTreeMGS = 0;
my $completionMarkerFastPaths = 0;

print "=====================================================\n";
print "Strain postprocessing v$version\n";
print "Mode: " . ($doSubmit ? "submit" : "dry run") . "; scheduler: $QSBoptHR->{qmode}; "
	. "rewrite existing results: " . ($rewriteRanalysis ? "yes" : "no") . "\n";
print "Inputs: gene catalog=$GCd; phylogenies=$FMGpD; map=$refMap; abundance matrix=$abMatrix\n";
print "Fallback outgroup tree: ".(length($MGSphylo) ? $MGSphylo : "<none>")."\n";
print "Paths: strain summary=$RsummaryTab; population summary=$popGenSummaryTab; subsampled population summary=$popGenSubsampleSummaryTab; toolkit=$MGSTKdir\n";
print "Resources: standard cores=$nCore; heavy-analysis cores=$nCoreHeavy; R-analysis memory="
	."${rAnalysisMemoryBaseGB}G + ${rAnalysisMemoryPerExtraCoreGB}G per core above $rAnalysisMemoryCoreThreshold; "
	."retry rounds=$rAnalysisRetryRounds; OOM ceiling=${rAnalysisOOMMaxMemoryGB}G; queue throttle=$checkMaxNumJobs jobs\n";
print "Metadata: individual=" . ($individualVar || "<none>")
	. "; family=" . ($familyVar || "<none>")
	. "; stability groups=" . ($groupStabilityVars || "<none>") . "\n";
print "Association tests: discrete=" . ($DiscTests || "<none>")
	. "; continuous=" . ($ContTests || "<none>") . "\n";
print "Population genetics: " . ($doPopGenStats
	? "enabled (subsamples: ".($popGenSubsample || "<none>")."; seed=$popGenSeed; genetic code=$popGenGeneticCode; codon start=$popGenCodonStart; strict outgroup=$popGenStrictOutgroup)"
	: "disabled") . "\n";
print "=====================================================\n";



if ($rewriteRanalysis){ #faster to do once for all..
	warn "Rewriting existing strain postprocessing results because -reSubmit 1 was requested\n";
	if ($doSubmit) {
		#print "Removing old strain2 analysis..\n";
		for my $within (glob("$FMGpD/*/within")) { remove_tree($within) if -d $within; }
		for my $summary ($RsummaryTab, $popGenSummaryTab, $popGenSubsampleSummaryTab, $legacyRsummaryTab) {
			unlink $summary or die "Cannot remove $summary: $!\n" if -e $summary;
		}
		my $networkDir = "$FMGpD/networks";
		remove_tree($networkDir) if -d $networkDir;
		my $treeWasCheckpoint = "$FMGpD/GeneEnrich/treeWAS.sto";
		unlink $treeWasCheckpoint or die "Cannot remove stale checkpoint $treeWasCheckpoint: $!\n" if -e $treeWasCheckpoint;
		my $phyloFigureCheckpoint = "$FMGpD/phyloFigures.sto";
		unlink $phyloFigureCheckpoint or die "Cannot remove stale checkpoint $phyloFigureCheckpoint: $!\n" if -e $phyloFigureCheckpoint;
	} else {
		print "Dry run: existing strain2 results and checkpoints will be preserved.\n";
	}
}


print "Scanning for completed within-MGS phylogenies\n";
opendir DIR, $FMGpD or die "Cannot open $FMGpD: $!\n";
#loop over intra-phylo dir and check for file presence..
my %treeNodes;
my %treeSamples;
while ( my $entry = readdir DIR ) {
    next if $entry eq '.' or $entry eq '..';
    next unless -d $FMGpD . '/' . $entry;
	my $phyloDirectory = "$FMGpD/$entry/phylo/";
	next unless -d $phyloDirectory;
	my $treeCompletion = "$FMGpD/$entry/treeDone.sto";
	if (-s $treeCompletion) {
		my $completedTreeSize = 0;
		my $completedTreeIndex = 0;
		my $completedTreePath = "";
		while ($completedTreeSize == 0 && $completedTreeIndex < @defTreeFiles) {
			# -s returns undef for a missing candidate; keep the numeric loop
			# sentinel defined while trying the remaining tree formats.
			$completedTreePath = "$phyloDirectory$defTreeFiles[$completedTreeIndex]";
			$completedTreeSize = (-s $completedTreePath) // 0;
			$completedTreeIndex++;
		}
		if ($completedTreeSize) {
			# Match strain_within.pl: this atomic BuildTree completion marker is
			# written only after the final primary tree validates successfully.
			$dirs{$entry} = $phyloDirectory;
			$baseD{$entry} = "$FMGpD/$entry";
			$treeNodes{$entry} = newickNodeCount($completedTreePath);
			$treeSamples{$entry} = newickSampleCount($completedTreePath);
			$completionMarkerFastPaths++;
			next;
		}
	}
	my @outcomeMarkers = grep {
		-s "$FMGpD/$entry/$_"
	} @nonTreeOutcomeMarkers;
	if (@outcomeMarkers) {
		$terminalTreeMGS++;
		next;
	}
	#my $destD = "$FMGpD/$entry/within/";
	#system "cp $destD/$entry.nwk $FMGpD/$entry/phylo/IQtree.treefile " if (-e "$destD/$entry.nwk");
	my $legacyTreeSize = 0; my $x=0;
	my $legacyTreePath = "";
	while ($legacyTreeSize == 0 && $x < @defTreeFiles){
		# A missing file is not an error: try the next supported tree format.
		$legacyTreePath = "$phyloDirectory$defTreeFiles[$x]";
		$legacyTreeSize = (-s $legacyTreePath) // 0;
		$x++;
	}
	next unless ($legacyTreeSize);
	#genuine legacy MGS phylo dir without a current completion marker
	$dirs{$entry} = $phyloDirectory;
	$baseD{$entry} = "$FMGpD/$entry";
	$treeNodes{$entry} = newickNodeCount($legacyTreePath);
	$treeSamples{$entry} = newickSampleCount($legacyTreePath);
}

closedir DIR or die "Cannot close $FMGpD: $!\n";
print "Found ".scalar(keys %dirs)." MGS directories with a nonempty calculated tree";
print "; completion-marker fast paths=$completionMarkerFastPaths";
print "; skipped $terminalTreeMGS MGS with valid no-tree or placement-pending markers"
	if $terminalTreeMGS;
print "\n";
my $MGstats = "$GCd/metagStats.txt";
$MGstats = "-1" unless (-e $MGstats);
my $treeAbsent = 0;
my @k2d = sort { $treeSamples{$b} <=> $treeSamples{$a} || $a cmp $b } keys(%treeSamples);
if (!@k2d) {
	print "No nonempty phylogenies found; skipping strain postprocessing.\n";
	exit 0;
}
my $largestMGS = $k2d[0];
my $batchSampleBudget = $treeSamples{$largestMGS};

# qsubSystem also enables -e and pipefail, but repeat the complete strict mode
# in every generated R-analysis script so a Slurm OOM exit is never ignored.
my $cmdPrelude = "set -euo pipefail\nulimit -s 20000\n";
my %analysisAttempted;
my ($reusedStrainStats, $reusedPopGenStats) = (0, 0);
my $strainStatsR = getProgPaths("treeSubGrpsR");
my $popGenStatsR = $doPopGenStats ? getProgPaths("pogenStats") : "";
my $combineResultsR = getProgPaths("combineResults_R");
my %outgroupSources;
my %jobsByAnalysis;
my %submittedRAnalysisJobs;
my %submittedByAnalysis;
my %rAnalysisBatches;
my $wrHead=1;
my $recordRAnalysisJob = sub {
	my ($analysisKind, $dependency, $label, $cores, $memory, $jobRecords) = @_;
	return '' unless $doSubmit && ($QSBoptHR->{qmode} || '') eq 'slurm';
	return '' unless defined($dependency) && length($dependency);
	my $jobID = $dependency;
	my $runTag = $QSBoptHR->{rTag} || '';
	$jobID =~ s/^\Q$runTag\E// if length($runTag);
	return '' unless $jobID =~ /^\d+$/;
	$jobRecords ||= ($submittedRAnalysisJobs{$analysisKind} ||= {});
	$jobRecords->{$jobID} = {
		requested_name => "$analysisKind.$label",
		cores => $cores,
		memory => $memory,
	};
	return $jobID;
};
my $submitRAnalysisBatch = sub {
	my ($analysisKind, $script, $batchCmd, $batchCores, $batchMemory, $batchLabel, $stores) = @_;
	my ($dep, $qcmd) = qsubSystem(
		$script, $batchCmd, $batchCores, $batchMemory, $batchLabel, "", "", 1, [], $QSBoptHR);
	my $jobID = $recordRAnalysisJob->(
		$analysisKind, $dep, $batchLabel, $batchCores, $batchMemory);
	push(@{$jobsByAnalysis{$analysisKind}}, $dep);
	push(@{$rAnalysisBatches{$analysisKind}}, {
		script => $script,
		command => $batchCmd,
		cores => $batchCores,
		memory => $batchMemory,
		label => $batchLabel,
		stores => { %{$stores} },
		dependency => $dep,
		job_id => $jobID,
	});
	return $dep;
};

for my $analysisKind (qw(strainStats popGenStats)) {
	next if $analysisKind eq 'popGenStats' && !$doPopGenStats;
	my $cnt=-1; my $curBatch=0; my $curBatchSamples=0;
	my $cmd = $cmdPrelude; my $destD = "";
	my $batchDestD = ""; my $batchLabel = "";
	my $phaseSubmitted = 0;
	my $batchStores = {};
foreach my $d (@k2d){#loop over MGS intra-phylo dirs, submit R analysis
	$cnt++;
	$destD = $dirs{$d};
	die "Unexpected phylogeny directory: $destD\n" unless $destD =~ s{/phylo/?$}{/within/};
	my $destBaseD = $dirs{$d};
	die "Unexpected phylogeny directory: $destBaseD\n" unless $destBaseD =~ s{/phylo/?$}{/};
	$destDs{$d} = $destD;
	#my $locTree = "$destD ../phylo/$defTreeFile"; #two args in one..
	my $treePath = ""; my $x=0;
	my $defTree="";
	while (!-s $treePath && $x < @defTreeFiles){
		$treePath = "$dirs{$d}/$defTreeFiles[$x]";
		$defTree = $defTreeFiles[$x];
		$x++;
	}
	if (!-s $treePath){
		$treeAbsent++ if $analysisKind eq 'strainStats';
		next;
	}
	if ($rewriteRanalysis && $analysisKind eq 'strainStats' && $doSubmit && -d $destD) {
		unlink $_ or die "Cannot remove $_: $!\n" for grep { -f $_ || -l $_ } glob("$destD/*");
		remove_tree($_) for grep { -d $_ } glob("$destD/*");
	}
	# The versioned RDS stores are the authoritative completion artifacts.
	my $analysisLog = "$destD/$d.Ranalysis.log";
	my $analysisReport = "$destD/$d.analysis.txt";
	my $analysisStore = "$destD/strainStats.output.Rds";
	my $popGenStore = "$destD/popGenStats.output.Rds";
	my $strainStatsReady = -s $analysisStore;
	my $popGenStatsReady = !$doPopGenStats || -s $popGenStore;
	if (!$rewriteRanalysis) {
		if ($analysisKind eq 'strainStats' && $strainStatsReady) {
			$reusedStrainStats++;
			next;
		}
		if ($analysisKind eq 'popGenStats' && $popGenStatsReady) {
			$reusedPopGenStats++;
			next;
		}
	}
	make_path($destD);
	my (undef, $OG, $outgroupSource) = resolveOutgroup($d, $destBaseD);
	$outgroupSources{$outgroupSource}++ if $analysisKind eq 'strainStats';
	#system "cp $dirs{$d}/$defTreeFile $destD/$d.nwk";
	my $BinN = 1000;
	if ($d =~ m/MB2bin(\d+)/){$BinN = $1;}
	# Use the requested standard core count for every R invocation. PopGenStats
	# applies this deterministically as its per-MSA PSOCK worker count.
	my $jobCores = $nCore;
	my $treeSampleCount = $treeSamples{$d};
	my $batchSampleCost = $treeSampleCount < $batchSampleBudget
		? 3 * $treeSampleCount : $treeSampleCount;
	if ($curBatchSamples > 0 && $curBatchSamples + $batchSampleCost > $batchSampleBudget) {
		qsubSystemWaitMaxJobs($checkMaxNumJobs,0,$QSBoptHR) if $doSubmit;
		my $batchCores = $nCore;
		my $batchMemory = rAnalysisMemoryForCores($batchCores);
		my $batchCmd = $cmd;
		$submitRAnalysisBatch->(
			$analysisKind, $batchDestD."$analysisKind.Ranalysis.sh",
			$batchCmd, $batchCores, $batchMemory, $batchLabel, $batchStores);
		$phaseSubmitted++;
		$curBatch = 0; $curBatchSamples = 0; $cmd=$cmdPrelude;
		$batchStores = {};
		$batchDestD = ""; $batchLabel = "";
	}
	$cmd .= "echo ".shellQuote("At tree $d")."\n";
	my $OGstr = $OG ne "" ? "--outgroup ".shellQuote($OG)." " : "";
	if ($analysisKind eq 'strainStats' && (!$strainStatsReady || $rewriteRanalysis)) {
		$cmd .= "if ! test -s ".shellQuote($analysisStore)."; then\n";
		$cmd .= "rm -f ".join(" ", map { shellQuote($_) } ($analysisLog, $analysisReport, $analysisStore))."\n";
		$cmd .= "echo ".shellQuote("Starting strainStats for $d with $jobCores core(s)")."\n";
		$cmd .= "$strainStatsR --path ".shellQuote($destD)." --tree ".shellQuote("../phylo/$defTree")." --taxN ".shellQuote($d)." $OGstr --map ".shellQuote($refMap)." --metagStats ".shellQuote($MGstats)." --abMat ".shellQuote($abMatrix)." --ncore $jobCores --siteMode 1 --MFDir ".shellQuote($MGSTKdir)." --wrColNms $wrHead --discPermTests ".shellQuote($DiscTests)." --contPermTests ".shellQuote($ContTests)." --familyCol ".shellQuote($familyVar)." --groupStabilityVars ".shellQuote($groupStabilityVars)." > ".shellQuote($analysisLog)."\n";
		$cmd .= "test -s ".shellQuote($analysisStore)."\n";
		$cmd .= "echo ".shellQuote("Completed strainStats for $d")."\n";
		$cmd .= "fi\n";
		$batchStores->{$d} = $analysisStore;
		$analysisAttempted{$d}{strainStats} = { store => $analysisStore };
		$wrHead=0;
	}
	if ($analysisKind eq 'popGenStats' && (!$popGenStatsReady || $rewriteRanalysis)) {
		$cmd .= "if ! test -s ".shellQuote($popGenStore)."; then\n";
		$cmd .= "rm -f ".shellQuote($popGenStore)."\n";
		$cmd .= "echo ".shellQuote("Starting popGenStats for $d")."\n";
		my $popGenCommand = "$popGenStatsR ".join(" ", map { shellQuote($_) }
			($destBaseD, $refMap, $destD));
		$popGenCommand .= " --subsample ".shellQuote($popGenSubsample)
			if length($popGenSubsample);
		$popGenCommand .= " --ncore $jobCores"
			." --genetic-code $popGenGeneticCode"
			." --codon-start $popGenCodonStart"
			." --seed $popGenSeed"
			." --individual-column ".shellQuote($individualVar);
		$popGenCommand .= " --outgroup ".shellQuote($OG) if $OG ne "";
		$popGenCommand .= " --strict-outgroup" if $popGenStrictOutgroup;
		$popGenCommand .= " --legacy-text-output" if $popGenLegacyTextOutput;
		my $strainFile = "$destD/IQtree_allsites.strains.txt";
		# PopGenStats owns validation and fallback behavior for this optional input.
		$cmd .= "$popGenCommand --strain-file ".shellQuote($strainFile)."\n";
		$cmd .= "test -s ".shellQuote($popGenStore)."\n";
		$cmd .= "echo ".shellQuote("Completed popGenStats for $d")."\n";
		$cmd .= "fi\n";
		$batchStores->{$d} = $popGenStore;
		$analysisAttempted{$d}{popGenStats} = { store => $popGenStore };
	}
	
	#print $cmd;
	#next;
	#system $cmd."\n";
	#$QSBoptHR->{useLongQueue} = 1;
	$curBatch++;
	$curBatchSamples += $batchSampleCost;
	$batchDestD = $destD; $batchLabel = "$analysisKind.R$cnt";
	if ($curBatchSamples >= $batchSampleBudget){
		qsubSystemWaitMaxJobs($checkMaxNumJobs,0,$QSBoptHR) if $doSubmit;
		my $batchCores = $nCore;
		my $batchMemory = rAnalysisMemoryForCores($batchCores);
		my $batchCmd = $cmd;
		$submitRAnalysisBatch->(
			$analysisKind, $batchDestD."$analysisKind.Ranalysis.sh",
			$batchCmd, $batchCores, $batchMemory, $batchLabel, $batchStores);
		$phaseSubmitted++;
		$curBatch = 0; $curBatchSamples = 0; $cmd=$cmdPrelude;
		$batchStores = {};
		$batchDestD = ""; $batchLabel = "";
	}
	#die;
	#last if ($cnt > 5);
}
if ($curBatch > 0){
	my $batchCores = $nCore;
	my $batchMemory = rAnalysisMemoryForCores($batchCores);
	my $batchCmd = $cmd;
	$submitRAnalysisBatch->(
		$analysisKind, $batchDestD."$analysisKind.Ranalysis.sh",
		$batchCmd, $batchCores, $batchMemory, $batchLabel, $batchStores);
	$curBatch = 0; $curBatchSamples = 0; $cmd=$cmdPrelude;
	$batchStores = {};
	$phaseSubmitted++;
}
	$submittedByAnalysis{$analysisKind} = $phaseSubmitted;
}

my $strainTaskCount = scalar grep { exists($analysisAttempted{$_}{strainStats}) } keys %analysisAttempted;
my $popGenTaskCount = scalar grep { exists($analysisAttempted{$_}{popGenStats}) } keys %analysisAttempted;
print "R-analysis plan: $strainTaskCount strainStats MGS analyses in "
	. ($submittedByAnalysis{strainStats} || 0) . " job script(s); "
	. "$popGenTaskCount PopGenStats MGS analyses in "
	. ($submittedByAnalysis{popGenStats} || 0) . " job script(s); "
	. "reused strainStats=$reusedStrainStats, reused PopGenStats=$reusedPopGenStats; "
	. "sample budget=$batchSampleBudget (all smaller MGS count three times)";
print ", $treeAbsent tree(s) became unavailable while planning" if $treeAbsent;
print "\n";

print "Outgroup recovery: ".join(", ", map {
	"$_=$outgroupSources{$_}"
} sort keys %outgroupSources)."\n" if %outgroupSources;
if (!$doSubmit) {
	print "Dry run: R-analysis submission scripts were prepared; summaries and downstream analyses were not run.\n";
	exit 0;
}

my $verifyAnalysisStores = sub {
	my ($analysisKind) = @_;
	my @failed = grep {
		exists($analysisAttempted{$_}{$analysisKind})
			&& !-s $analysisAttempted{$_}{$analysisKind}{store};
	} sort keys %analysisAttempted;
	die "$analysisKind failed or produced incomplete RDS output for: ".join(", ", @failed)."\n"
		if @failed;
};
my $batchHasMissingStores = sub {
	my ($batch) = @_;
	return scalar grep { !-s $batch->{stores}{$_} }
		keys %{ $batch->{stores} || {} };
};
my $rAnalysisOOMBatches = sub {
	my ($batches) = @_;
	return {} unless $doSubmit && ($QSBoptHR->{qmode} || '') eq 'slurm';
	my @records;
	for my $batch (@{$batches}) {
		next unless defined($batch->{job_id}) && $batch->{job_id} =~ /^\d+$/;
		my $memory = $batch->{memory} || '';
		my $requestedMB = 0;
		$requestedMB = int($1 * 1024 + 0.5) if $memory =~ /^([\d.]+)G$/i;
		$requestedMB = int($1 + 0.5) if $memory =~ /^([\d.]+)M$/i;
		next unless $requestedMB > 0;
		push @records, {
			job_id => $batch->{job_id},
			requested_mb => $requestedMB,
			mgs => $batch->{label},
		};
	}
	return {} unless @records;
	my $summary = slurm_tree_memory_summary(\@records);
	if (!$summary->{available}) {
		warn "Could not query Slurm memory accounting for R-analysis OOM recovery: "
			.($summary->{error} || 'unknown error')."\n";
		return {};
	}
	return { map { $_->{job_id} => 1 } @{$summary->{oom_jobs} || []} };
};
my $waitForAnalysis = sub {
	my ($analysisKind) = @_;
	my @allBatches = @{$rAnalysisBatches{$analysisKind} || []};
	my @roundBatches = @allBatches;
	my $roundJobRecords = $submittedRAnalysisJobs{$analysisKind} || {};
	my $round = 0;
	while (1) {
		my @roundJobs = grep { defined($_) && length($_) }
			map { $_->{dependency} } @roundBatches;
		if (@roundJobs) {
			print "\n\nwaiting for $analysisKind ".($round ? "retry round $round" : "jobs")." to finish\n";
			qsubSystemJobAlive(\@roundJobs, $QSBoptHR);
		}
		reportRAnalysisSchedulerFailures($roundJobRecords, $QSBoptHR, $analysisKind);
		my @retryBatches = grep { $batchHasMissingStores->($_) } @allBatches;
		last unless @retryBatches;
		last if $round >= $rAnalysisRetryRounds;
		$round++;
		my $oomBatches = $rAnalysisOOMBatches->(\@roundBatches);
		my $oomCount = scalar grep { $oomBatches->{$_->{job_id}} } @retryBatches;
		print "Retrying $analysisKind round $round/$rAnalysisRetryRounds for ".scalar(@retryBatches)
			." incomplete batch(es)".($oomCount ? "; increasing memory for $oomCount OOM batch(es)" : '')."\n";
		my %retryJobRecords;
		my @nextRoundBatches;
		for my $batch (@retryBatches) {
			my $retryMemory = $batch->{memory};
			if ($oomBatches->{$batch->{job_id}}) {
				my $increasedMemory = increaseRAnalysisMemory($batch->{memory});
				if (defined($increasedMemory)) {
					$retryMemory = $increasedMemory;
				} else {
					warn "R-analysis OOM retry memory ceiling reached for $batch->{label}; retaining $retryMemory\n";
				}
			}
			qsubSystemWaitMaxJobs($checkMaxNumJobs, 0, $QSBoptHR) if $doSubmit;
			my $retryScript = $batch->{script};
			$retryScript =~ s/\.sh\z/.retry$round.sh/;
			$retryScript .= ".retry$round.sh" if $retryScript eq $batch->{script};
			my $retryLabel = "$batch->{label}.retry$round";
			my ($dep, $qcmd) = qsubSystem(
				$retryScript, $batch->{command}, $batch->{cores}, $retryMemory,
				$retryLabel, "", "", 1, [], $QSBoptHR);
			my $jobID = $recordRAnalysisJob->(
				$analysisKind, $dep, $retryLabel, $batch->{cores}, $retryMemory, \%retryJobRecords);
			$batch->{dependency} = $dep;
			$batch->{job_id} = $jobID;
			$batch->{memory} = $retryMemory;
			$batch->{label} = $retryLabel;
			push @nextRoundBatches, $batch;
		}
		@roundBatches = @nextRoundBatches;
		$roundJobRecords = \%retryJobRecords;
	}
	$verifyAnalysisStores->($analysisKind);
};

$waitForAnalysis->('strainStats');
my $shouldCombineStrainStats = $rewriteRanalysis || $strainTaskCount > 0 || !-s $RsummaryTab;
if ($shouldCombineStrainStats) {
	combineResults(0);
} else {
	print "Reusing existing combined strainStats overview: $RsummaryTab\n";
}

## Start strain-only downstream work while PopGenStats remains in the scheduler.
my ($networkDep, $networkStone) = strainNetwork();
my ($treeWasDep, $treeWasStone) = treeWas();
my ($phyloFigureDep, $phyloFigureStone) = visualizeSignPhylos();

if ($doPopGenStats) {
	$waitForAnalysis->('popGenStats');
	my $shouldCombinePopGenStats = $rewriteRanalysis || $popGenTaskCount > 0 || !-s $popGenSummaryTab;
	if ($shouldCombinePopGenStats) {
		combineResults(1);
	} else {
		print "Reusing existing combined PopGenStats overview: $popGenSummaryTab\n";
	}
}

my @postAnalysisJobs = grep { defined($_) && length($_) }
	($networkDep, $treeWasDep, $phyloFigureDep);
if (@postAnalysisJobs) {
	print "Waiting for network, treeWAS, and phylogeny-figure analyses to finish\n";
	qsubSystemJobAlive(\@postAnalysisJobs, $QSBoptHR);
}
my @missingPostAnalysisStones = grep { !-e $_ }
	($networkStone, $treeWasStone, $phyloFigureStone);
die "Downstream strain analyses failed or did not publish checkpoints: "
	.join(", ", @missingPostAnalysisStones)."\n"
	if @missingPostAnalysisStones;

if (0){#get within strain nuc div
	my $countStrains=0;
	open O,">$FMGpD/withinStrain.tab";
	foreach my $d (@k2d){
		my $destBaseD = $dirs{$d}; $destBaseD =~ s/(.*)\/phylo/$1\//; 
		my $strainFile = "$destBaseD/within/IQtree_allsites.strains.txt";
		next unless (-e $strainFile);
		print O "$d\tNA\n";
		my $wiStF="$destBaseD/codeml/WithinStrainDiv.txt";
		open I,"<$wiStF" or die "can't find file $wiStF\n";
		while (<I>){ print O $_;}
		close I;
		my %seensStrains;
		open I,"<$strainFile" or die "can;t open strain file $strainFile\n";
		while (my $lin = <I>){
			chomp $lin; my @spl = split /\t/,$lin;
			#print $lin."\n"; # if (@spl < 2 || !defined($spl[1]));
			$seensStrains{$spl[1]} = 1;# unless ($seensStrains{$spl[1]} eq "NA");
		}
		close I;
		$countStrains += (scalar(%seensStrains)-1);
	}
	close O;
}

my %TS; 
#--------------------------------------------------------------
#summary of dnds
if (0){
	#fubar summaries
	sumSummaries("hyphy.fubar.txt","$FMGpD/fubar.tab");
	sumSummaries("hyphy.fubar.s10.txt","$FMGpD/fubar.s10.tab");
	sumSummaries("hyphy.fubar.s20.txt","$FMGpD/fubar.s20.tab");
	sumSummaries("hyphy.fubar.s30.txt","$FMGpD/fubar.s30.tab");
	sumSummaries("hyphy.fubar.s100.txt","$FMGpD/fubar.s100.tab");
	sumSummaries("hyphy.fubar.s200.txt","$FMGpD/fubar.s200.tab");
	sumSummaries("hyphy.fubar.s500.txt","$FMGpD/fubar.s500.tab");
	sumSummaries("hyphy.fubar.unID.txt","$FMGpD/fubar.unID.tab");
	#popstats
	sumSummaries("PopStats.txt","$FMGpD/PopStats.tab");
	sumSummaries("PopStats.10.txt","$FMGpD/PopStats.10.tab");
	sumSummaries("PopStats.20.txt","$FMGpD/PopStats.20.tab");
	sumSummaries("PopStats.30.txt","$FMGpD/PopStats.30.tab");
	sumSummaries("PopStats.100.txt","$FMGpD/PopStats.100.tab");
	sumSummaries("PopStats.200.txt","$FMGpD/PopStats.200.tab");
	sumSummaries("PopStats.500.txt","$FMGpD/PopStats.500.tab");
	#sumSummaries("PopStats.500.txt","$FMGpD/PopStats.unID.tab");
	sumSummaries("hyphy.Theta.log","$FMGpD/Theta.tab");
}


print "\nFinished with strain postprocessing\n";
exit (0);

#--------------------------------------------------------------
#now starts the real work, take cluster files and make trees based on this subset.
#this part depends on pre created sets of fastas containing all genes for a given species
#then it will resort the cat file to only include samples deduced from R tree
foreach my $d (@k2d){
	my $clsts = "$destDs{$d}/${d}.cl_IndFam.txt";
	my %clusters;
	open I,"<$clsts" or die "can't open $clsts\n";
	while (<I>){
		chomp; my @spl = split /\t/;
		push(@{$clusters{$spl[1]}}, $spl[0]);
	}
	close I;
	
	my $fnFile = "$baseD{$d}/allFNAs.fna";
	my $faFile = "$baseD{$d}/allFAAs.faa";
	my $catFile = "$baseD{$d}/all.cat";
	if (!-e $fnFile || !-e $catFile){
		die "Can't find input files:\n$catFile\n$fnFile\n";
	}
	my %genePres;
	open I,"<$catFile" or die "Can't open the category file $catFile\n";
	while (<I>){
		chomp;
		my @spl = split /\t/;
		foreach my $gn (@spl){
			my ($sample, $locus) = split /\Q$SaSe\E/, $gn, 2;
			unless (defined($sample) && length($sample) && defined($locus) && length($locus)) {
				next;
			}
			$genePres{$locus}{$sample} = 1; #{COG|primaryGeneID}{SMPL}
		}
	}
	close I;
	my $cnt=0;
	#sort out which genes go into same tree..
	foreach my $clN (keys %clusters){
		my @smpls = @{$clusters{$clN}};
		my $nOD = "$baseD{$d}/Cl_$clN/";
		system "mkdir -p $nOD" unless (-d $nOD);
		my $clOF = "$nOD/cl$clN.cat";
		#TODO: add 3 smpls from other cluster as outliers..
		my $outg = "";
		open O,">$clOF";
		#each line one cog
		foreach my $cog (keys %genePres){
			my @cp;
			foreach my $s (@smpls){
				if (exists($genePres{$cog}{$s})){ #gene exists for sample, put in cat file
					push (@cp,"$s$SaSe$cog");
				}
			}
			if (@cp > 0 ){
				print O join("\t",@cp)."\n";
			}
		}
		close O;
		#cat file created, now submit job
		my $Tcmd= "$bts -fna $fnFile -aa $faFile -smplSep '\\$SaSe' -cats $clOF -outD $nOD -runIQtree 1 -iqFast 1 -runFastTree 0 -cores $numCores  ";
		$Tcmd .= "-withinSpecies 1 -AAtree 0 -bootstrap 000 -NTfiltCount 3000 -NTfilt 0.1 -NTfiltPerGene 0.6 -runRaxMLng 0 -minOverlapMSA 2 -MSAprogram 2 -AutoModel 0 \n";
		#die "$cmd\n" if ($cnt ==10);
		$QSBoptHR->{useLongQueue} = 1;
		my ($dep,$qcmd) = qsubSystem($nOD."subtree_$clN.sh",$Tcmd,$numCores,"1G","FT$cnt","","",1,[],$QSBoptHR);
		$cnt ++;
		if ($outg ne ""){
			open O,">$nOD/TODO";
			close O;
		}
	}
	
}






exit(0); #don't do FST for now..
#--------------------------------------------------------------
#calculate FST between countries / cities..
#depends on Alex's script for calculating FST
my $mapF = resolve_catalog_maps($GCd);
my ($hr1,$hr2) = readMapS($mapF,-1,"Country");
my %map = %{$hr1}; my %AsGrps = %{$hr2};
my $FSTbin = getProgPaths("FSTpy");
foreach my $d (@k2d){
	my $clsts = "$destDs{$d}/${d}.cl_IndFam.txt";
	my $cntrFcnt=0;
	my %cntr;
	open I,"<$clsts" or die "can't open $clsts\n";
	while (<I>){
		chomp; my @spl = split /\t/;
		if (exists($map{$spl[0]}{"Country"})){
			my $tcnt=$map{$spl[0]}{"Country"};
			push(@{$cntr{$tcnt}},$spl[0]) ;
			$cntrFcnt++;
		}
	}
	close I;
	#start FST script from alex
	my $FSTd = "$baseD{$d}/FST/";
	system "mkdir -p $FSTd" unless (-d $FSTd);
	my @presCntrs = sort keys %cntr;
	open O, ">$FSTd/FST.cntry.guide" or die $!;
	foreach my $cn (@presCntrs){
		print O "$cn\t".join(",",@{$cntr{$cn}}) . "\n";
	}
	close O;
	my $refpop = "USA";
	if (!exists($cntr{USA})){$refpop = $presCntrs[0];}
	my $Fcmd = "$FSTbin -i $baseD{$d}/MSA/MSAli.fna -b 100 -p $FSTd/FST.cntry.guide --refpop $refpop -o  $FSTd/FST.cntry.test\n";
	die "$Fcmd\n";
}





sub sumSummaries($ $){
	my ($inF,$outF) = @_;
	open O1,">$outF" or die "Can't open $outF\n";
	my $first=1;
	foreach my $d (@k2d){
		#next;
		my $destD = $dirs{$d}; $destD =~ s/(.*)\/phylo/$1\/codeml/; 
		if (-e "$destD/$inF"){
			open I,"<$destD/$inF";my $tmp = <I>; 
			print O1 "\t$tmp" if ($first);while (<I>){print O1 "$d\t$_";	}close I; 
		}
		$first=0;
	}
	close O1; 
}




sub strainNetwork{ #submits Anthony's script to build a network
	my $netDir = "$FMGpD/networks/";
	my $networkStone = "$netDir/networks.sto";
	my $networkGraph = "$netDir/strain_graph.Rds";
	my $dep = "";
	if (-e $networkStone && !-s $networkGraph) {
		warn "Ignoring incomplete network checkpoint $networkStone\n";
		unlink $networkStone or die "Cannot remove incomplete network checkpoint $networkStone: $!\n";
	}
	if (!-e $networkStone){
		my $networkScr = getProgPaths("runNetworks_R");#"Rscript $MGSTKdir/runNetworks.R";
		make_path($netDir);
		my $edgeTresh = 4;
		my $cmd = "$networkScr -i ".shellQuote($FMGpD)." -o ".shellQuote($netDir)." -m ".shellQuote($refMap)." -e $edgeTresh\n";
		$cmd .= "#consider the following options to change: -c [Column for clustering samples] -e [num shared strains for edges]\n";
		$cmd .= "test -s ".shellQuote($networkGraph)."\n";
		$cmd .= "touch ".shellQuote($networkStone)."\n";
		print "Submitting shared-strain network analysis; output: $netDir\n";
		#system $cmd;
		my $nCore = 1;
		my ($submittedDep,$qcmd) = qsubSystem($netDir."Network.sh",$cmd,$nCore,"20G","Network","","",1,[],$QSBoptHR);
		$dep = $submittedDep;

	}
	return ($dep, $networkStone);
}




sub treeWas{
	# functional enrichments of strains in conditions defined by user
	my $funCmd = "";
	my $treewasRun_R = getProgPaths("treewasRun_R");
	my $processTreewas_R = getProgPaths("processTreewas_R");
	my $treewasOut = "$FMGpD/GeneEnrich/";
	my $treewasOutfile = "$treewasOut/treeWAS_results.csv";
	my $summaryOutfile = "$treewasOut/treeWAS_results_functions.csv";
	my $treeWasStone = "$FMGpD/GeneEnrich/treeWAS.sto";
	my $MGSd = dirname($FMGpD);
	if (-e $treeWasStone && (!-s $treewasOutfile || !-s $summaryOutfile)) {
		warn "Ignoring incomplete treeWAS checkpoint $treeWasStone\n";
		unlink $treeWasStone or die "Cannot remove incomplete treeWAS checkpoint $treeWasStone: $!\n";
	}
	make_path($treewasOut);
	$funCmd .= "mkdir -p ".shellQuote($treewasOut)."\n";
	$funCmd .= "rm -f ".join(" ", map { shellQuote($_) }
		($treewasOutfile, $summaryOutfile, $treeWasStone))."\n";
	$funCmd .= "#1st command: run treewas job\n";
	$funCmd .= "$treewasRun_R --gene_cat_dir ".shellQuote($GCd)." --n_threads $nCoreHeavy --metadata_vars ".shellQuote($groupStabilityVars)." -o ".shellQuote($treewasOut)." --mgs_dir ".shellQuote($MGSd)." --metadata_file ".shellQuote($refMap)." -r ".shellQuote($MGSTKdir)." -i ".shellQuote($individualVar)."\n";
	$funCmd .= "#2nd command: process results\n";
	$funCmd .= "$processTreewas_R -i ".shellQuote($treewasOutfile)." --gene_cat_dir ".shellQuote($GCd)." --annot_files ".shellQuote("NOG,CZy,KGM")." --out_file ".shellQuote($summaryOutfile)." --n_threads $nCoreHeavy -r ".shellQuote($MGSTKdir)."\n";
	$funCmd .= "test -s ".shellQuote($treewasOutfile)."\n";
	$funCmd .= "test -s ".shellQuote($summaryOutfile)."\n";
	$funCmd .= "touch ".shellQuote($treeWasStone)."\n";

	my $dep = "";
	if (!-e $treeWasStone){
		my ($submittedDep,$qcmd) = qsubSystem($treewasOut."treeWAS.sh",$funCmd,$nCoreHeavy,"6G","treewas","","",1,[],$QSBoptHR);
		$dep = $submittedDep;
	}
	return ($dep, $treeWasStone);
}



sub visualizeSignPhylos{
	#my $taxFile = "$GCd/Anno/Tax/GTDBmg_MGS/specI.tax";
	my $vizPhylos = getProgPaths("vizPhylosSign_R");
	my $taxFile = "$FMGpD/../Annotation/MGS.GTDB.LCA.tax";
	my $phyloFigureStone = "$FMGpD/phyloFigures.sto";
	my $cmdPic = "$vizPhylos ".join(" ", map { shellQuote($_) }
		($RsummaryTab, $taxFile, $FMGpD, "phylo", $MGSTKdir, $refMap, "-1"))."\n";
	$cmdPic .= "touch ".shellQuote($phyloFigureStone)."\n";

	my $dep = "";
	if (!-e $phyloFigureStone) {
		print "Submitting figures of most significant phylogenies\nThis might take several hours..\n";
		my ($submittedDep, $qcmd) = qsubSystem(
			"$FMGpD/phyloFigures.sh", $cmdPic, 1, "24G", "phyloFigures", "", "", 1, [], $QSBoptHR);
		$dep = $submittedDep;
	}
	return ($dep, $phyloFigureStone);
}

sub resolveOutgroup {
	my ($mgs, $directory) = @_;
	my ($logged, $outgroup) = loggedOutgroup($directory);
	return (1, $outgroup, 'data.log') if $logged;

	my ($saved, $savedOutgroup) = savedCommandOutgroup($directory);
	return (1, $savedOutgroup, 'treeCmd.sh') if $saved;

	$outgroup = mgsTreeOutgroup($mgs);
	return (1, $outgroup, 'MGSphylo') if length($outgroup);
	return (0, '', 'none');
}

sub loggedOutgroup {
	my ($directory) = @_;
	my $logPath = "$directory/data.log";
	return (0, '') unless -e $logPath || -e "$logPath.gz";
	my ($log_fh) = gzipopen($logPath, "strain outgroup log");
	while (my $line = <$log_fh>) {
		$line =~ s/[\r\n]+\z//;
		if ($line =~ /^OG:(.*)\z/) {
			close $log_fh or die "Cannot close strain outgroup log for $directory: $!\n";
			return (1, $1);
		}
	}
	close $log_fh or die "Cannot close strain outgroup log for $directory: $!\n";
	return (0, '');
}

sub savedCommandOutgroup {
	my ($directory) = @_;
	my $commandPath = "$directory/treeCmd.sh";
	return (0, '') unless -s $commandPath;
	open my $command_fh, '<', $commandPath
		or die "Cannot read saved tree command $commandPath: $!\n";
	local $/;
	my $command = <$command_fh> // '';
	close $command_fh or die "Cannot close saved tree command $commandPath: $!\n";
	my @tokens = eval { shellwords($command) };
	return (0, '') if $@ || !@tokens;
	return (0, '') unless grep { $_ eq '-withinSpecies' } @tokens;
	for my $index (0 .. $#tokens - 1) {
		next unless $tokens[$index] eq '-outgroup' || $tokens[$index] eq '--outgroup';
		return (1, $tokens[$index + 1]);
	}
	# A saved within-species command without this option records an intentional
	# ingroup-only tree, so do not replace it with a later heuristic.
	return (1, '');
}

sub mgsTreeOutgroup {
	my ($mgs) = @_;
	return '' unless length($MGSphylo) && length($neighborTree);
	my $call = "$neighborTree ".shellQuote($MGSphylo)." ".shellQuote($mgs);
	my $candidates = `$call`;
	if ($? != 0) {
		warn "Cannot recover the outgroup for $mgs from $MGSphylo; command failed: $call\n";
		return '';
	}
	for my $candidate (split /\s+/, $candidates) {
		return $candidate
			if $candidate =~ /\A[A-Za-z0-9][A-Za-z0-9_.:+-]*\z/;
	}
	warn "No usable outgroup candidate for $mgs was returned from $MGSphylo\n";
	return '';
}

sub combineResults {
	die "combineResults_R is not configured\n" unless length($combineResultsR);
	my ($includePopGen) = @_;
	$includePopGen = 0 unless defined($includePopGen);
	my $command = "$combineResultsR --path ".shellQuote($FMGpD)
		." --outDir ".shellQuote($FMGpD)." --include-popgen $includePopGen\n";
	my $combinedKinds = $includePopGen ? 'strainStats and PopGenStats' : 'strainStats';
	print "Combining per-MGS $combinedKinds result stores into $RsummaryTab\n";
	systemW($command);
	die "combineResults.R did not produce the overview table $RsummaryTab\n"
		unless -s $RsummaryTab;
	print "Combined strain overview: $RsummaryTab\n";
	if ($includePopGen) {
		die "combineResults.R did not produce the population overview table $popGenSummaryTab\n"
			unless -s $popGenSummaryTab;
		print "Combined population-genetics overview: $popGenSummaryTab\n";
		if (length($popGenSubsample) && -s $popGenSubsampleSummaryTab) {
			print "Combined subsampled population-genetics overview: $popGenSubsampleSummaryTab\n";
		} elsif (length($popGenSubsample)) {
			warn "No subsampled population-genetics rows were available for $popGenSubsampleSummaryTab\n";
		}
	}
}

sub newickSampleCount {
	my ($treePath) = @_;
	open my $tree_fh, '<', $treePath or die "Cannot read tree $treePath: $!\n";
	local $/;
	my $newick = <$tree_fh> // '';
	close $tree_fh or die "Cannot close tree $treePath: $!\n";
	$newick =~ s/\[[^\]]*\]//gs;
	$newick =~ s/\s+//g;
	return 1 unless length($newick);
	my $commas = () = $newick =~ /,/g;
	my $tips = $commas + 1;
	return $tips > 0 ? $tips : 1;
}

sub newickNodeCount {
	my ($treePath) = @_;
	open my $tree_fh, '<', $treePath or die "Cannot read tree $treePath: $!\n";
	local $/;
	my $newick = <$tree_fh> // '';
	close $tree_fh or die "Cannot close tree $treePath: $!\n";
	$newick =~ s/\[[^\]]*\]//gs;
	$newick =~ s/\s+//g;
	return 1 unless length($newick);
	my $internal = () = $newick =~ /\(/g;
	my $commas = () = $newick =~ /,/g;
	my $tips = $commas + 1;
	my $nodes = $tips + $internal;
	return $nodes > 0 ? $nodes : 1;
}

sub rAnalysisMemoryForCores {
	my ($cores) = @_;
	die "R-analysis core count must be positive\n" unless defined($cores) && $cores > 0;
	my $extraCores = $cores - $rAnalysisMemoryCoreThreshold;
	$extraCores = 0 if $extraCores < 0;
	my $memoryGB = $rAnalysisMemoryBaseGB
		+ $extraCores * $rAnalysisMemoryPerExtraCoreGB;
	return $memoryGB."G";
}

sub increaseRAnalysisMemory {
	my ($memory) = @_;
	return undef unless defined($memory);
	my $currentMB;
	$currentMB = int($1 * 1024 + 0.5) if $memory =~ /^([\d.]+)G$/i;
	$currentMB = int($1 + 0.5) if $memory =~ /^([\d.]+)M$/i;
	return undef unless defined($currentMB) && $currentMB > 0;
	my $maximumMB = int($rAnalysisOOMMaxMemoryGB * 1024 + 0.5);
	my $nextMB = next_oom_retry_memory_mb($currentMB, $maximumMB);
	return defined($nextMB) ? $nextMB.'M' : undef;
}

sub reportRAnalysisSchedulerFailures {
	my ($jobs, $options, $analysisKind) = @_;
	$analysisKind ||= 'R-analysis';
	my $emptySummary = { failed => 0, categories => {} };
	return $emptySummary unless $options->{doSubmit} && ($options->{qmode} || '') eq 'slurm';
	return $emptySummary unless scalar(keys %{$jobs});
	my $summary = slurmJobFailureSummary($jobs, $options);
	if (!defined($summary)) {
		warn "Could not query Slurm accounting for $analysisKind jobs\n";
		return $emptySummary;
	}
	return $summary unless $summary->{failed};
	my @failures;
	for my $category (sort keys %{$summary->{categories}}) {
		my $counts = $summary->{categories}{$category}{failures} || {};
		push @failures, $category.": ".join(", ", map {
			"$_=$counts->{$_}"
		} sort keys %{$counts});
	}
	warn "Slurm reported failed $analysisKind job(s): ".join("; ", @failures)
		.". Inspect the matching $analysisKind.Ranalysis.sh.etxt and *.Ranalysis.log files; "
		."incomplete result batches will be retried and OOM batches receive additional memory.\n";
	return $summary;
}

sub shellQuote {
	my ($value) = @_;
	die "Cannot shell-quote an undefined value\n" unless defined $value;
	$value =~ s/'/'"'"'/g;
	return "'$value'";
}
