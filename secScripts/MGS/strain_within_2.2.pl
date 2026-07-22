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
use File::Path qw(make_path remove_tree);
use File::Copy qw(copy);
use File::Basename qw(basename dirname);

sub sumSummaries;
sub strainNetwork;
sub treeWas;
sub visualizeSignPhylos;


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
my $version = 0.32;

my $rewriteRanalysis = 0; my $doSubmit = 1;
my $checkMaxNumJobs = 400;

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

GetOptions(
	"GCd=s"          => \$GCd,
	"map=s"          => \$refMap,
	"submit=i"       => \$doSubmit,
	"reSubmit=i"     => \$rewriteRanalysis,
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
) or die "Invalid strain_within_2.2.pl options\n";
die "Unexpected positional arguments: @ARGV\n" if @ARGV;
die "-GCd, -map, -FMGdir and -MGSmatrix are required\n"
	unless length($GCd) && length($refMap) && length($FMGpD) && length($abMatrix);
die "Gene-catalog and phylogeny directories must exist\n" unless -d $GCd && -d $FMGpD;
die "Map or MGS abundance matrix is missing\n" unless -s $refMap && -s $abMatrix;
die "Core requests must be positive\n" unless $nCore > 0 && $nCoreHeavy > 0;
die "-submit and -reSubmit must be 0 or 1\n" unless $doSubmit =~ /^[01]$/ && $rewriteRanalysis =~ /^[01]$/;
$FMGpD =~ s{/+$}{} unless $FMGpD eq "/";


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
my $SaSe = "|";
my $SaSe2 = "\|"; #for regex
my $numCores = 40; #used for phylos..
my $RsummaryTab = "$FMGpD/Rsummary.tab";
my $existingSummaryHeader = '';
if (!$rewriteRanalysis && -s $RsummaryTab) {
	open my $old_summary, '<', $RsummaryTab or die "Cannot read $RsummaryTab: $!\n";
	$existingSummaryHeader = <$old_summary> // '';
	close $old_summary or die "Cannot close $RsummaryTab: $!\n";
}


my %dirs;my %destDs; my %baseD;

print "Strain analysis v $version\n";



if ($rewriteRanalysis){ #faster to do once for all..
	print "\nWARNING:: Rewriting strain2 results!\n" ;
	if ($doSubmit) {
		#print "Removing old strain2 analysis..\n";
		for my $within (glob("$FMGpD/*/within")) { remove_tree($within) if -d $within; }
		unlink $RsummaryTab or die "Cannot remove $RsummaryTab: $!\n" if -e $RsummaryTab;
		for my $checkpoint ("$FMGpD/networks/networks.sto", "$FMGpD/GeneEnrich/treeWAS.sto") {
			unlink $checkpoint or die "Cannot remove stale checkpoint $checkpoint: $!\n" if -e $checkpoint;
		}
	} else {
		print "Dry run: existing strain2 results and checkpoints will be preserved.\n";
	}
}


print "Reading dirs..\n";
opendir DIR, $FMGpD or die "Cannot open $FMGpD: $!\n";
#loop over intra-phylo dir and check for file presence..
my %sizTrees;
while ( my $entry = readdir DIR ) {
    next if $entry eq '.' or $entry eq '..';
    next unless -d $FMGpD . '/' . $entry;
	next unless (-d "$FMGpD/$entry/phylo/");
	#my $destD = "$FMGpD/$entry/within/";
	#system "cp $destD/$entry.nwk $FMGpD/$entry/phylo/IQtree.treefile " if (-e "$destD/$entry.nwk");
	my $sizTree = 0; my $x=0;
	while ($sizTree == 0 && $x < @defTreeFiles){
		$sizTree = -s "$FMGpD/$entry/phylo/$defTreeFiles[$x]" if (-e "$FMGpD/$entry/phylo/$defTreeFiles[$x]");
		$x++;
	}
	next unless ($sizTree);
	#genuine MGS phylo dir-> store in %dirs %baseD
	$dirs{$entry} = "$FMGpD/$entry/phylo/"; 
	$baseD{$entry} = "$FMGpD/$entry";
	$sizTrees{$entry} = $sizTree;
}

closedir DIR;
print "Found ".scalar(keys %dirs)." dirs with calculated tree\n";
my $cnt=-1; my $curBatch=0; my $batchSize = 0; my $submitted=0;my @jobs;
my $MGstats = "$GCd/metagStats.txt";
$MGstats = "-1" unless (-e $MGstats);
my $treeAbsent = 0;
#my @k2d = sort keys %dirs;
my @k2d = sort { $sizTrees{$b} <=> $sizTrees{$a} } keys(%sizTrees);
if (!@k2d) {
	print "No nonempty phylogenies found; skipping strain postprocessing.\n";
	exit 0;
}

my $cmdPrelude = "ulimit -s 20000\n";
my $cmd = $cmdPrelude;my $destD =""; my $wrHead=1; my $headerOwner='';
my %analysisAttempted; my $legacyCompleted = 0;
my $strainStatsR = getProgPaths("treeSubGrpsR");

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
	# Accept legacy log+report completions so existing successful runs are not
	# needlessly repeated. New jobs also write an explicit success stone below.
	my $analysisLog = "$destD/$d.Ranalysis.log";
	my $analysisReport = "$destD/$d.analysis.txt";
	my $analysisStone = "$destD/$d.Ranalysis.sto";
	if (!$rewriteRanalysis && -s $analysisReport && (-e $analysisStone || -s $analysisLog)) {
		$legacyCompleted++ if !-e $analysisStone;
		next;
	}
	if ($doSubmit && -d $destD) {
		unlink $_ or die "Cannot remove $_: $!\n" for grep { -f $_ || -l $_ } glob("$destD/*");
		remove_tree($_) for grep { -d $_ } glob("$destD/*");
	}
	make_path($destD);
	my $OG = "";
	if (-e "$destBaseD/data.log" || -e "$destBaseD/data.log.gz") {
		my ($log_fh) = gzipopen("$destBaseD/data.log", "strain outgroup log");
		while (my $line = <$log_fh>) {
			if ($line =~ /^OG:(.*)$/) { $OG = $1; chomp $OG; last; }
		}
		close $log_fh;
	}
	#system "cp $dirs{$d}/$defTreeFile $destD/$d.nwk";
	my $BinN = 1000;
	if ($d =~ m/MB2bin(\d+)/){$BinN = $1;}
	my $jobCores = $nCore;
	
	$cmd .= "rm -f ".join(" ", map { shellQuote($_) } ($analysisLog, $analysisReport, $analysisStone))."\n";
	$cmd .= "echo ".shellQuote("At tree $d")."\n";
	my $OGstr = $OG ne "" ? "--outgroup ".shellQuote($OG)." " : "";
	$cmd .= "$strainStatsR --path ".shellQuote($destD)." --tree ".shellQuote("../phylo/$defTree")." --taxN ".shellQuote($d)." $OGstr --map ".shellQuote($refMap)." --metagStats ".shellQuote($MGstats)." --abMat ".shellQuote($abMatrix)." --ncore $jobCores --siteMode 1 --MFDir ".shellQuote($MGSTKdir)." --wrColNms $wrHead --discPermTests ".shellQuote($DiscTests)." --contPermTests ".shellQuote($ContTests)." --familyCol ".shellQuote($familyVar)." --groupStabilityVars ".shellQuote($groupStabilityVars)." > ".shellQuote($analysisLog)."\n";
	$cmd .= "test -s ".shellQuote($analysisReport)."\n";
	$cmd .= "touch ".shellQuote($analysisStone)."\n";
	$analysisAttempted{$d} = { report => $analysisReport, stone => $analysisStone };
	$headerOwner = $d if $wrHead;
	$wrHead=0;
	if (0){#rerun popgen stats??
		my $RpogenS = getProgPaths("pogenStats");
		$cmd .= "$RpogenS $destBaseD $refMap $destBaseD/codeml/ $destBaseD/MSA/clnd/ 10,20,30,100,200,500\n";
	}
	
	#print $cmd;
	#next;
	#system $cmd."\n";
	#$QSBoptHR->{useLongQueue} = 1;
	print "$d: "; 
	$curBatch++;
	if ($curBatch > $batchSize){
		
		qsubSystemWaitMaxJobs($checkMaxNumJobs,0,$QSBoptHR) if $doSubmit;

		my ($dep,$qcmd) = qsubSystem($destD."Ranalysis.sh",$cmd,$jobCores,"20G","R$cnt","","",1,[],$QSBoptHR);
		#die " $destD\n";
		push(@jobs,$dep);
		$curBatch = 0; $cmd=$cmdPrelude;
		$submitted++;
	}
	#die;
	#last if ($cnt > 5);
}
if ($curBatch > 0){
	my ($dep,$qcmd) = qsubSystem($destD."Ranalysis.sh",$cmd,$nCore,"20G","R$cnt","","",1,[],$QSBoptHR);
	$curBatch = 0; $cmd=$cmdPrelude;
	push(@jobs,$dep);
}

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
	!-s $analysisAttempted{$_}{report} || !-e $analysisAttempted{$_}{stone}
} sort keys %analysisAttempted;
die "R analysis failed or produced incomplete output for: ".join(", ", @failedAnalysis)."\n"
	if @failedAnalysis;

print "$treeAbsent phylos absent\n";

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
	print "Total Strains seen: $countStrains\n";
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
#summary  of R stats
#create summary tables 
if (1 || !-e $RsummaryTab){
	my $summaryHeader = $existingSummaryHeader;
	if ($summaryHeader eq '' && $headerOwner ne '') {
		my $headerReport = "$destDs{$headerOwner}/${headerOwner}.analysis.txt";
		open my $header_fh, '<', $headerReport or die "Cannot read header report $headerReport: $!\n";
		$summaryHeader = <$header_fh> // '';
		close $header_fh or die "Cannot close header report $headerReport: $!\n";
	}
	open my $summary_fh, ">", $RsummaryTab or die "Cannot reset $RsummaryTab: $!\n";
	print {$summary_fh} $summaryHeader if $summaryHeader ne '';

	foreach my $d (@k2d){
		my $clsts = "$destDs{$d}/${d}.Ranalysis.log";
		if ($cpDir ne ""){
			make_path("$cpDir/$d");
			foreach my $source (glob("$destDs{$d}/${d}*")) {
				next unless -f $source;
				copy($source, "$cpDir/$d/".basename($source))
					or die "Cannot copy $source to $cpDir/$d: $!\n";
			}
		}
		my $SCtrig=0;
		
		my $TXTreport = "$destDs{$d}/${d}.analysis.txt";

		if (-e $TXTreport && -s $TXTreport){
			open my $report_fh, "<", $TXTreport or die "Cannot read $TXTreport: $!\n";
			while (my $line = <$report_fh>) {
				next if $summaryHeader ne '' && $line eq $summaryHeader;
				print {$summary_fh} $line;
			}
			close $report_fh or die "Cannot close $TXTreport: $!\n";
		}
	}
	close $summary_fh or die "Cannot close $RsummaryTab: $!\n";
}

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



#die "$FMGpD/Rsummary.tab";

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
				warn "Ignoring malformed tree sequence identifier '$gn' in $catFile\n";
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
		$Tcmd .= "-AAtree 0 -bootstrap 000 -NTfiltCount 3000 -NTfilt 0.1 -NTfiltPerGene 0.6 -runRaxMLng 0 -minOverlapMSA 2 -MSAprogram 2 -AutoModel 0 \n";
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
my $mapF = `cat $GCd/LOGandSUB/GCmaps.inf`;
#$mapF = $GCd."LOGandSUB/inmap.txt" if ($mapF eq "");
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
		} else {
			print "not found: $spl[0]\n";
		}
	}
	close I;
	print "$d : $cntrFcnt country\n";
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
		} else {
			print "missing: $destD/$inF\n";
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
		print "Running network of shared strains..\n$cmd\n";
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
	print $cmdPic;
	systemW($cmdPic);
}

sub shellQuote {
	my ($value) = @_;
	die "Cannot shell-quote an undefined value\n" unless defined $value;
	$value =~ s/'/'"'"'/g;
	return "'$value'";
}
