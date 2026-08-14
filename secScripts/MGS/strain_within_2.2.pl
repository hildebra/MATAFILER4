#!/usr/bin/perl
#submits R scripts to work through phylos
#args: [GeneCat dir] [intra_phylo_dir] [abundance MGS] [mapping file[ [cores]

use warnings;
use strict;

use Getopt::Long qw( GetOptions );

use Mods::GenoMetaAss qw( readClstrRev systemW readMapS readFasta gzipopen);
use Mods::Subm qw(qsubSystem emptyQsubOpt qsubSystemJobAlive qsubSystemWaitMaxJobs);
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
sub newickNodeCount;
sub resolveOutgroup;
sub loggedOutgroup;
sub savedCommandOutgroup;
sub mgsTreeOutgroup;

my $MGSTKdir = getProgPaths("MGSTKDir");

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
my $version = 0.39;

my $rewriteRanalysis = 0; my $doSubmit = 1;
my $checkMaxNumJobs = 400;
my $doPopGenStats = 1;
my $popGenSubsample = "10,20,30,100,200,500";

my $DiscTests =""; my $ContTests = "";
my $familyVar = ""; my $groupStabilityVars = "";
my $GCd = "";#$ARGV[0];
my $nCore = 4;#$ARGV[4];
my $nCoreHeavy = 12;
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
die "Core requests must be positive\n" unless $nCore > 0 && $nCoreHeavy > 0;
die "-submit and -reSubmit must be 0 or 1\n" unless $doSubmit =~ /^[01]$/ && $rewriteRanalysis =~ /^[01]$/;
die "-popGenStats must be 0 or 1\n" unless $doPopGenStats =~ /^[01]$/;
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
print "Resources: standard cores=$nCore; heavy-analysis cores=$nCoreHeavy; queue throttle=$checkMaxNumJobs jobs\n";
print "Metadata: individual=" . ($individualVar || "<none>")
	. "; family=" . ($familyVar || "<none>")
	. "; stability groups=" . ($groupStabilityVars || "<none>") . "\n";
print "Association tests: discrete=" . ($DiscTests || "<none>")
	. "; continuous=" . ($ContTests || "<none>") . "\n";
print "Population genetics: " . ($doPopGenStats ? "enabled (subsamples: ".($popGenSubsample || "<none>").")" : "disabled") . "\n";
print "=====================================================\n";



if ($rewriteRanalysis){ #faster to do once for all..
	warn "Rewriting existing strain postprocessing results because -reSubmit 1 was requested\n";
	if ($doSubmit) {
		#print "Removing old strain2 analysis..\n";
		for my $within (glob("$FMGpD/*/within")) { remove_tree($within) if -d $within; }
		for my $summary ($RsummaryTab, $popGenSummaryTab, $popGenSubsampleSummaryTab, $legacyRsummaryTab) {
			unlink $summary or die "Cannot remove $summary: $!\n" if -e $summary;
		}
		for my $checkpoint ("$FMGpD/networks/networks.sto", "$FMGpD/GeneEnrich/treeWAS.sto") {
			unlink $checkpoint or die "Cannot remove stale checkpoint $checkpoint: $!\n" if -e $checkpoint;
		}
	} else {
		print "Dry run: existing strain2 results and checkpoints will be preserved.\n";
	}
}


print "Scanning for completed within-MGS phylogenies\n";
opendir DIR, $FMGpD or die "Cannot open $FMGpD: $!\n";
#loop over intra-phylo dir and check for file presence..
my %treeNodes;
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
}

closedir DIR or die "Cannot close $FMGpD: $!\n";
print "Found ".scalar(keys %dirs)." MGS directories with a nonempty calculated tree";
print "; completion-marker fast paths=$completionMarkerFastPaths";
print "; skipped $terminalTreeMGS MGS with valid no-tree or placement-pending markers"
	if $terminalTreeMGS;
print "\n";
my $cnt=-1; my $curBatch=0; my $curBatchNodes=0; my $submitted=0;my @jobs;
my $MGstats = "$GCd/metagStats.txt";
$MGstats = "-1" unless (-e $MGstats);
my $treeAbsent = 0;
#my @k2d = sort keys %dirs;
my @k2d = sort { $treeNodes{$b} <=> $treeNodes{$a} || $a cmp $b } keys(%treeNodes);
if (!@k2d) {
	print "No nonempty phylogenies found; skipping strain postprocessing.\n";
	exit 0;
}
my $batchNodeBudget = $treeNodes{$k2d[0]};

my $cmdPrelude = "ulimit -s 20000\n";
my $cmd = $cmdPrelude;my $destD =""; my $wrHead=1;
my %analysisAttempted; my $legacyCompleted = 0; my $reusedAnalysis = 0;
my $strainStatsR = getProgPaths("treeSubGrpsR");
my $popGenStatsR = $doPopGenStats ? getProgPaths("pogenStats") : "";
my $combineResultsR = getProgPaths("combineResults_R");
my %outgroupSources;
my $batchDestD = ""; my $batchLabel = "";

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
		$treeAbsent++;
		next;
	}
	if ($rewriteRanalysis && $doSubmit && -d $destD) {
		unlink $_ or die "Cannot remove $_: $!\n" for grep { -f $_ || -l $_ } glob("$destD/*");
		remove_tree($_) for grep { -d $_ } glob("$destD/*");
	}
	# The versioned result store is the authoritative strainStats completion artifact.
	# New jobs also write an explicit success stone below.
	my $analysisLog = "$destD/$d.Ranalysis.log";
	my $analysisReport = "$destD/$d.analysis.txt";
	my $analysisStone = "$destD/$d.Ranalysis.sto";
	my $analysisStore = "$destD/strainStats.output.Rds";
	my $popGenStone = "$destD/$d.popGen.sto";
	my $popGenStore = "$destD/popGenStats.output.Rds";
	my $strainStatsReady = -s $analysisStore && (-e $analysisStone || -s $analysisLog);
	my $popGenStatsReady = !$doPopGenStats || -s $popGenStore;
	if (!$rewriteRanalysis && $strainStatsReady && $popGenStatsReady) {
		$legacyCompleted++ if !-e $analysisStone;
		$reusedAnalysis++;
		next;
	}
	if ($doSubmit && -d $destD) {
		unlink $_ or die "Cannot remove $_: $!\n" for grep { -f $_ || -l $_ } glob("$destD/*");
		remove_tree($_) for grep { -d $_ } glob("$destD/*");
	}
	make_path($destD);
	my (undef, $OG, $outgroupSource) = resolveOutgroup($d, $destBaseD);
	$outgroupSources{$outgroupSource}++;
	#system "cp $dirs{$d}/$defTreeFile $destD/$d.nwk";
	my $BinN = 1000;
	if ($d =~ m/MB2bin(\d+)/){$BinN = $1;}
	my $jobCores = $nCore;
	my $treeNodeCount = $treeNodes{$d};
	if ($curBatchNodes > 0 && $curBatchNodes + $treeNodeCount > $batchNodeBudget) {
		qsubSystemWaitMaxJobs($checkMaxNumJobs,0,$QSBoptHR) if $doSubmit;
		my ($dep,$qcmd) = qsubSystem($batchDestD."Ranalysis.sh",$cmd,$nCore,"20G",$batchLabel,"","",1,[],$QSBoptHR);
		push(@jobs,$dep);
		$submitted++;
		$curBatch = 0; $curBatchNodes = 0; $cmd=$cmdPrelude;
		$batchDestD = ""; $batchLabel = "";
	}
	$cmd .= "echo ".shellQuote("At tree $d")."\n";
	my $OGstr = $OG ne "" ? "--outgroup ".shellQuote($OG)." " : "";
	if (!$strainStatsReady || $rewriteRanalysis) {
		$cmd .= "rm -f ".join(" ", map { shellQuote($_) } ($analysisLog, $analysisReport, $analysisStone, $analysisStore))."\n";
		$cmd .= "$strainStatsR --path ".shellQuote($destD)." --tree ".shellQuote("../phylo/$defTree")." --taxN ".shellQuote($d)." $OGstr --map ".shellQuote($refMap)." --metagStats ".shellQuote($MGstats)." --abMat ".shellQuote($abMatrix)." --ncore $jobCores --siteMode 1 --MFDir ".shellQuote($MGSTKdir)." --wrColNms $wrHead --discPermTests ".shellQuote($DiscTests)." --contPermTests ".shellQuote($ContTests)." --familyCol ".shellQuote($familyVar)." --groupStabilityVars ".shellQuote($groupStabilityVars)." > ".shellQuote($analysisLog)."\n";
		$cmd .= "test -s ".shellQuote($analysisStore)."\n";
		$cmd .= "touch ".shellQuote($analysisStone)."\n";
		$analysisAttempted{$d}{strainStats} = { store => $analysisStore, stone => $analysisStone };
		$wrHead=0;
	}
	if ($doPopGenStats && (!$popGenStatsReady || $rewriteRanalysis)) {
		$cmd .= "rm -f ".join(" ", map { shellQuote($_) } ($popGenStone, $popGenStore))."\n";
		$cmd .= "$popGenStatsR ".join(" ", map { shellQuote($_) } ($destBaseD, $refMap, $destD));
		$cmd .= " --subsample ".shellQuote($popGenSubsample) if length($popGenSubsample);
		$cmd .= "\n";
		$cmd .= "test -s ".shellQuote($popGenStore)."\n";
		$cmd .= "touch ".shellQuote($popGenStone)."\n";
		$analysisAttempted{$d}{popGenStats} = { store => $popGenStore, stone => $popGenStone };
	}
	
	#print $cmd;
	#next;
	#system $cmd."\n";
	#$QSBoptHR->{useLongQueue} = 1;
	$curBatch++;
	$curBatchNodes += $treeNodeCount;
	$batchDestD = $destD; $batchLabel = "R$cnt";
	if ($curBatchNodes >= $batchNodeBudget){
		qsubSystemWaitMaxJobs($checkMaxNumJobs,0,$QSBoptHR) if $doSubmit;
		my ($dep,$qcmd) = qsubSystem($batchDestD."Ranalysis.sh",$cmd,$nCore,"20G",$batchLabel,"","",1,[],$QSBoptHR);
		push(@jobs,$dep);
		$submitted++;
		$curBatch = 0; $curBatchNodes = 0; $cmd=$cmdPrelude;
		$batchDestD = ""; $batchLabel = "";
	}
	#die;
	#last if ($cnt > 5);
}
if ($curBatch > 0){
	my ($dep,$qcmd) = qsubSystem($batchDestD."Ranalysis.sh",$cmd,$nCore,"20G",$batchLabel,"","",1,[],$QSBoptHR);
	$curBatch = 0; $curBatchNodes = 0; $cmd=$cmdPrelude;
	push(@jobs,$dep);
	$submitted++;
}

print "R-analysis plan: " . scalar(keys %analysisAttempted) . " MGS analyses prepared in "
	. "$submitted job script(s), $reusedAnalysis existing result(s) reused; size budget=$batchNodeBudget nodes";
print ", $treeAbsent tree(s) became unavailable while planning" if $treeAbsent;
print "\n";

print "Outgroup recovery: ".join(", ", map {
	"$_=$outgroupSources{$_}"
} sort keys %outgroupSources)."\n" if %outgroupSources;
if (!$doSubmit) {
	print "Dry run: R-analysis submission scripts were prepared; summaries and downstream analyses were not run.\n";
	exit 0;
}

if (@jobs){ #wait for all submitted R scripts, then continue in script
	
	#automate
	
	print "\n\nwaiting for R analysis to finish before subclustering step\n";
	qsubSystemJobAlive( \@jobs,$QSBoptHR );

}

print "Accepted $legacyCompleted legacy R-analysis completions without success stones.\n"
	if $legacyCompleted;
my @failedAnalysis = grep {
	my $artifacts = $analysisAttempted{$_};
	grep { !-s $artifacts->{$_}{store} || !-e $artifacts->{$_}{stone} } keys %$artifacts;
} sort keys %analysisAttempted;
die "R analysis failed or produced incomplete output for: ".join(", ", @failedAnalysis)."\n"
	if @failedAnalysis;

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

#die;
# Combine named result-store summaries rather than positional text reports.
combineResults();

## run network of similar samples
my ($networkDep, $networkStone) = strainNetwork();

#die;

# functional enrichments of strains in conditions defined by user
my ($treeWasDep, $treeWasStone) = treeWas();

my @postAnalysisJobs = grep { defined($_) && length($_) } ($networkDep, $treeWasDep);
if (@postAnalysisJobs) {
	print "Waiting for network and treeWAS analyses to finish\n";
	qsubSystemJobAlive(\@postAnalysisJobs, $QSBoptHR);
}
my @missingPostAnalysisStones = grep { !-e $_ } ($networkStone, $treeWasStone);
die "Downstream strain analyses failed or did not publish checkpoints: "
	.join(", ", @missingPostAnalysisStones)."\n"
	if @missingPostAnalysisStones;


visualizeSignPhylos();



#die "$FMGpD/strainStats.tsv";

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
	my $dep = "";
	if (!-e $networkStone){
		my $networkScr = getProgPaths("runNetworks_R");#"Rscript $MGSTKdir/runNetworks.R";
		make_path($netDir);
		my $edgeTresh = 4;
		my $cmd = "$networkScr -i ".shellQuote($FMGpD)." -o ".shellQuote($netDir)." -m ".shellQuote($refMap)." -e $edgeTresh\n";
		$cmd .= "#consider the following options to change: -c [Column for clustering samples] -e [num shared strains for edges]\n";
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
	my $cmdPic = "$vizPhylos ".join(" ", map { shellQuote($_) } ($RsummaryTab,$taxFile,$FMGpD,"phylo",$MGSTKdir,$refMap,"-1"))."\n";

	print "Printing figures of most significant phylogenies\nThis might take several hours..\n";
	systemW($cmdPic);
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
	my $command = "$combineResultsR --path ".shellQuote($FMGpD)
		." --outDir ".shellQuote($FMGpD)."\n";
	print "Combining per-MGS strainStats result stores into $RsummaryTab\n";
	systemW($command);
	die "combineResults.R did not produce the overview table $RsummaryTab\n"
		unless -s $RsummaryTab;
	print "Combined strain overview: $RsummaryTab\n";
	if ($doPopGenStats) {
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

sub shellQuote {
	my ($value) = @_;
	die "Cannot shell-quote an undefined value\n" unless defined $value;
	$value =~ s/'/'"'"'/g;
	return "'$value'";
}
