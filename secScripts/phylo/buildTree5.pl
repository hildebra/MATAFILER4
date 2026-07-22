#!/usr/bin/env perl

#ARGS: ./buildTree.pl -fna [FNA] -faa [FAA] -cat [categoryFile] -outD [outDir] -cores [CPUs] -useEte [1=ETE,0=this script] -NTfilt [filter]
#versions: ver 2 makes a link to nexus file formats, to be used in MrBayes and BEAST etc
#8.12.17: added mod3 from Mechthild
#2.1.20: rewrite of workflow to extend superTrees to superCheck
#version 5 added: hyphy fubar, R scripts for theta, guidance2
#2.1.25: v5.02: reduced threshold for including genes from MGS
#14.2.25: v5.03: added treshrink
#18.2.26: 5.04: MSA gz'ping
#2.3.26: 5.05: added veryFastTRee
#5.06: 15.4.26: added famse
#5.07: 16.7.26: validate paths/options and repair filtering, resume, optional-tree, and cleanup paths

use warnings;
use strict;
#use threads ('yield','stack_size' => 64*4096,'exit' => 'threads_only','stringify');
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::GenoMetaAss qw( fileGZe fileGZs gzipopen systemW readFasta readFastHD writeFasta quantile);
use Mods::phyloTools qw(convertMSA2NXS MSA filterMSA getTreeLeafs calcDisPos2 runRaxML runRaxMLng runQItree 
			runFasttree runVeryFasttree fixHDs4Phylo getGenoGenes getFMG readFMGdir );
			
			
use Getopt::Long qw( GetOptions );
#use Mods::ext::TreeIO;
#use Mods::IO::MaybeXS qw(encode_json decode_json);
use Mods::IO::PP qw (decode_json);
#use JSON qw( decode_json ); 
use Data::Dumper;
use Mods::math qw (medianArray avgArray meanArray);
use Cwd qw(abs_path);
use File::Basename qw(basename);
use File::Copy qw(copy move);
use File::Path qw(make_path remove_tree);
use File::Spec;


sub convertMultAli2NT;
sub mergeMSAs;
sub synPosOnly;
sub calcDisPos;#gets only the dissimilar positions of an MSA, as well as %id similarity
sub calcDisPos2;#de novo aligns pairwise via vsearch and calcs id (iddef 2)
sub calcDiffDNA;
sub selecAnalysis;
sub runFastgear;
sub mergePids; sub WattTheta;
sub singleGeneMSAprocess;
sub pruneTree;
sub prepGenoDirs;
sub createTreeOpt;
sub treePresent;
sub parseSeqId;
sub geneFileStem;
sub safeRemoveTree;
sub requireConfiguredTool;
sub shellQuote;

my $doPhym= 0;
my $version = 5.07;

my $pigzBin  = getProgPaths("pigz");
#my $trDist = getProgPaths("treeDistScr");


my $gubbinsBin = ""; # resolved lazily if the dormant Gubbins mode is reactivated
my $pamlBin = "";    # resolved lazily when codeml is requested
#my $evoConda = getProgPaths("evoEnv");#systemW "source $evoConda";



my $partiExt=".partition.RAXML";

#die "TODO $trimalBin\n";
#trimal -in /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5/T2/tesssst/MSA/COG0185.faa -out /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5/T2/tesssst/MSA/tst.fna -backtrans /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5/T2/tesssst/inMSA0.fna -keepheader -keepseqs -noallgaps -automated1 -ignorestopcodon
#some runtim options...
#my $ncore = 20;#RAXML cores
my $ntFrac =0.2; my $ntFracGene = 0.1; 
my $GeneFracPSpec = 0.1; #replacement for ntFracGene, as works also with supertrees
my $MSAprog = 2; #do MSA with clustal (1) or msaprobs (0), mafft(2), guidance2(3), MUSCLE5 (4)
my $calcDistMat = 0; #distmat of either AA or NT (depending on MSA)
my $calcDistMatExt = 0; #distmat of other AA or NT (depending on MSA), e.g. running two times an MSA
my $calcDistMatExtGo = 0;
my $treeAutoModel=1; #iqtree: choose model automatically (a bit slower)
my $fracMaxGenesFilter = 0.2;
my $fracMaxGenes90pct = 0.25; #gene cats to keep, e.g. 25% of 90th percentile


my $ntCntTotal =0; my $bootStrap=0; my $subsetSmpls = -1;
my ($fnFna, $aaFna,$cogCats,$outD,$ncore,$Ete, $filt,$smplDef,$smplSep,$calcSyn,$calcNonSyn,
			$useAA4tree,$calcDNAdiff,$tmpD ) = ("","","","",1,0,0.8,1,"_",0,0,0,0,"");

my ($continue,$isAligned) = (0,0);#overwrite already existing files?
my $outgroup="";
my $fixHeaders = 0;
my ($doGubbins,$doCFML,$doRAXML,$doFastTree,$doVeryFastTree, $doIQTree,$doRAXMLng) = (0,0,0,0,0, 0, 0);#fastree as default tree builder

#check length of fasta to avoid frameshift
my $doLengthCheck=1;
my $doDNDS=0;
my $doTheta=0;
my $mapF = "";
my $doGenesToPh=0;
my $doFastGear=0;
my $doFastGearSummary=0;
my $postFilter = "";
my $clusterName="";
my $MSAreq = 1;
my $iqFast=0;
my $minOverlapMSA = 0;
my $maxGapPerCol = 1 ;
my $minPcId = 0;
my $doSuperTree =0;
my $doSuperCheck=0;#check if tree's of single genes behave "strange"
my $gzipInput =0; my $removeMSA = 0;
my $useTreeShrink =0;


#EDBUGGIN
my $reparseHyphyJson = 1;


## parameters for pogen stats
my $selGene=0; 		#run dnds just on subset (given by genesForDNDS) of genes 
my @genesExtra;	#list with selected genes just for dnds
my @model= (0,1,2);  # codeml models: 0 -> neutral selection,  (1,2,7,8) -> test for postive selection
my @omegas = (0.3, 1.3); #try different omega values to check for convergence	
my $repeatCounts=2;	#set how often each model should be repeated to check for convergence 
my $codemlOutD=""; 
my $treeFile="";
my $genoindir = "";
my $wildcardflag = "";#"/*\.fna";
my $subsetPopgenStats = "10,20,30,100,200,500"; #maybe add later to options..
my $MSAsubsD =""; #only needed to create subsets of MSAs for hyphy etc calcs


die "no input args!\n" if (@ARGV == 0 );


GetOptions(
	"genoInD=s" => \$genoindir, #provide a dir with complete genomes, will extract FGMs and build tree between genomes (NT/AA flag)
	"wildcardflag=s" => \$wildcardflag,
	"fna=s" => \$fnFna,
	"aa=s"      => \$aaFna,
	"cats=s"      => \$cogCats,
	"outD=s"      => \$outD,
	"tmpD=s" => \$tmpD,
	"cores=i" => \$ncore,
	"superTree=i" => \$doSuperTree,
	"superCheck=i" => \$doSuperCheck,
	"fixHeaders=i" => \$fixHeaders, ## fix the fasta headers, if too long or containing not allowed symbols (nwk reserved)
	"useEte=i"      => \$Ete,
	"NTfilt=f"      => \$filt,
	"NTfiltPerGene=f"      => \$ntFracGene,
	"GenesPerSpecies=f" => \$GeneFracPSpec,
	"fracMaxGenes90pct=f" => \$fracMaxGenes90pct,
	"NTfiltCount=i" => \$ntCntTotal,
	"smplDef=i"	=> \$smplDef, #is the genome somehow quantified with a delimiter (_) ?
	"smplSep=s" => \$smplSep, #set the delimiter
	"outgroup=s"	=> \$outgroup,
	"AAtree=i" => \$useAA4tree,
	"MSAprogram=i" => \$MSAprog, #(0) MSAprobs, (1) clustalO, (2) mafft, (4) MUSCLE5, (5) FAMSA2 (only AA)
	"minOverlapMSA=i" => \$minOverlapMSA, #min overlap in MSA columns, in order to retain column
	"maxGapPerCol=f" =>\$maxGapPerCol, #same as minOverlapMSA, but for MSAfix and %of gaps allowed in a column
	"calcDistMat=i" => \$calcDistMat,
	"calcDistMatExt=i" => \$calcDistMatExt,
	"calcDiffDNA=i" => \$calcDNAdiff,
	"minPcId=f" => \$minPcId, #sequence is filtered from data, unless the average minPcId is >= $minPcId
	"SynTree=i"	=> \$calcSyn,
	"NonSynTree=i"	=> \$calcNonSyn,
	"continue=i" => \$continue,
	"bootstrap=i" => \$bootStrap,
	"subsetSmpls=i" => \$subsetSmpls,
	"postFilter=s" => \$postFilter, # "," sep list of zorro,guidance2,macse
	"rmMSA=i" => \$removeMSA, #to save diskspace
	"gzInput=i" => \$gzipInput, #to save diskspace
	"isAligned=i" => \$isAligned,
	"runRAxML=i" => \$doRAXML,
	"runRaxMLng=i" => \$doRAXMLng,
	"runFastTree=i" => \$doFastTree,
	"runVeryFastTree=i" => \$doVeryFastTree,
	"treeShrink=i" => \$useTreeShrink,
	"runIQtree=i" => \$doIQTree,
	"AutoModel=i" => \$treeAutoModel,
	"iqFast=i" => \$iqFast, #fast qiTree mode
	"runClonalFrameML=i" => \$doCFML,
	"runGubbins=i" => \$doGubbins,
	"runLengthCheck=i" => \$doLengthCheck,		#check that sequence length can be divided by 3
	"runDNDS=i" => \$doDNDS,			#run dNdS analysis
	"runTheta=i" => \$doTheta,
	"genesForDNDS=s{,}" => \@genesExtra,		#list with selected genes just for dnds
	"DNDSonSubset=i" => \$selGene,			#run dnds just on subset (given by genesForDNDS) of genes
	"codemlRepeats=i" => \$repeatCounts,
#	"treefileCodeml=s" => \$treeFile,
	"outDCodeml=s"=> \$codemlOutD,
	"genesToPhylip=i" => \$doGenesToPh,	
	"runFastgear=i" => \$doFastGear,
	"runFastGearPostProcessing=i" => \$doFastGearSummary,
	"map=s" =>\$mapF,
	"clustername=s" => \$clusterName,
) or die("Error in command line arguments\n");
die "Unexpected positional arguments: @ARGV\n" if @ARGV;

die "-cores must be a positive integer\n" if $ncore < 1;
die "-bootstrap must be zero or greater\n" if $bootStrap < 0;
die "-NTfiltCount must be zero or greater\n" if $ntCntTotal < 0;
die "-smplSep must not be empty\n" if $smplSep eq "";
eval { qr/$smplSep/ } or die "Invalid -smplSep regular expression '$smplSep': $@";
for my $fraction_name_value (
	["NTfilt", $filt], ["NTfiltPerGene", $ntFracGene],
	["GenesPerSpecies", $GeneFracPSpec], ["fracMaxGenes90pct", $fracMaxGenes90pct],
	["maxGapPerCol", $maxGapPerCol],
) {
	my ($name, $value) = @{$fraction_name_value};
	die "-$name must be between 0 and 1\n" if $value < 0 || $value > 1;
}
die "Unsupported -MSAprogram $MSAprog (expected 0, 1, 2, 4, or 5)\n"
	unless grep { $MSAprog == $_ } (0, 1, 2, 4, 5);


######### indir
if ($genoindir ne ""){
	my $genoindir2 = $genoindir;
	if (!-d $genoindir2){ $genoindir2 =~ s/[^\/]+$//;}
	if ($outD eq ""){$outD = $genoindir2."/phylo/";}
}
die "-outD is required (unless it can be derived from -genoInD)\n" if $outD eq "";
$outD = File::Spec->canonpath(File::Spec->rel2abs($outD));
my ($outVolume) = File::Spec->splitpath($outD, 1);
my $volumeRoot = File::Spec->canonpath(File::Spec->catpath($outVolume, File::Spec->rootdir, ""));
die "Refusing to use filesystem root '$outD' as -outD\n" if $outD eq $volumeRoot;
make_path($outD) unless -d $outD;
die "Output path is not a directory: $outD\n" unless -d $outD;

print "BuildTree 5 script v$version\nOutDir: $outD\n";

##### setup dirs
$codemlOutD = File::Spec->catdir($outD, "codeml") if ($codemlOutD eq "");

my $tmpBase = $tmpD eq "" ? File::Spec->catdir($outD, "tmp") : File::Spec->canonpath(File::Spec->rel2abs($tmpD));
make_path($tmpBase) unless -d $tmpBase;
die "Temporary path is not a directory: $tmpBase\n" unless -d $tmpBase;
my $tmpTag = $clusterName eq "" ? "default" : $clusterName;
$tmpTag =~ s/[^A-Za-z0-9_.-]+/_/g;
$tmpD = File::Spec->catdir($tmpBase, "buildTree5_${tmpTag}_$$");
make_path($tmpD);
my $treeD = File::Spec->catdir($outD, "phylo");#raxml, fasttree, phyml tree output dir

my $MsaD = File::Spec->catdir($outD, "MSA");
if ($removeMSA){
	$MsaD = File::Spec->catdir($tmpD, "MSA_$tmpTag");
}

$MSAsubsD = "$MsaD/clnd/";


if ($subsetSmpls >0){
	$MsaD =~ s/\/$/_S$subsetSmpls\//;
	$treeD =~ s/\/$/_S$subsetSmpls\//;
}
######

if ($MSAprog == 0){print "Warning:  MSAprobs with trimal gives warnings (ignore them)\n";}
if ($doCFML && !$doRAXML){die "Need RaxML alignment, if Clonal fram is to be run..\n";}

if ($aaFna eq "" || $useAA4tree){	$calcSyn=0;$calcNonSyn=0;}
if ($filt <1){$ntFrac=$filt; print "Using filter with $ntFrac fraction of nts\n";}
if ($outgroup ne ""){print "Using outgroup $outgroup\n";}
if ($bootStrap>0){print "Using bootstrapping in tree building\n";}
#if (($calcDistMat || $calcDistMatExt) && $isAligned || !$MSAprog){die"Can't calc distance mat, unless clustalO is being used for MSA\n";}
#else {$ntCntTotal = $filt;}

$MSAreq = 0 if (!$doFastTree && !$doVeryFastTree && !$doRAXML && !$doRAXMLng && !$doCFML && !$doGubbins && !$doIQTree);

make_path($tmpD) unless -d $tmpD;
my $cmd =""; my %usedGeneNms;


my $outD_clust = File::Spec->catdir($outD, "fastGear_work_$tmpTag");

#------------------------------------------
#sorting by COG, MSA & syn position extraction
if ($Ete){
	my $eteBin = getProgPaths("ete3");
	my $eteOut = File::Spec->catdir($outD, "tree");
	make_path($eteOut);
	$cmd = "$eteBin build -n $fnFna -a $aaFna -w clustalo_default-none-none-none  -m sptree_raxml_all --cpu $ncore -o $eteOut --clearall --nt-switch 0.0 --noimg";
	$cmd .= " --cogs $cogCats" unless ($cogCats eq "");
	print "Running tree analysis ..";
	print $cmd."\n";
	systemW($cmd . " > $eteOut/ETE.log 2>&1");
	print " Done.\n$eteOut\n";
	exit(0);
}


#general routine, starting with MSA
#followed by merge of MSA / NT -> AA conversion, MSA filter etc
#followed by tree building

if ($fixHeaders){
	if ($cogCats ne ""){die"implement fix hds for cats\naborting\n";}
	$aaFna=fixHDs4Phylo($aaFna);$fnFna = fixHDs4Phylo($fnFna); 
}


prepGenoDirs($genoindir);
for my $input_spec (["fna", $fnFna], ["aa", $aaFna], ["cats", $cogCats]) {
	my ($name, $path) = @{$input_spec};
	next if $path eq "";
	die "-$name input does not exist or is empty: $path\n" unless fileGZs($path);
}
die "A sequence input (-aa or -fna) is required\n" if $aaFna eq "" && $fnFna eq "";
die "Category-based alignments require -aa\n" if $cogCats ne "" && $aaFna eq "";
die "Nucleotide trees require -fna\n" if !$useAA4tree && $fnFna eq "";

if (!$continue){
	safeRemoveTree($treeD, $outD);
	safeRemoveTree($MsaD, $removeMSA ? $tmpD : $outD);
}
make_path($MsaD) unless -d $MsaD;
make_path($treeD) unless -d $treeD;
my $multAli = "$MsaD/MSAli.fna";
my $multAliSyn = $multAli.".syn.fna";
my $multAliNonSyn = $multAli.".nonsyn.fna";
my @theRealMSAs;
my $partiFile="";#partitioning for multi gene MSAs
my %specList; #list of species (without _COG00012 tag);
my %samples; 
my $MSAcat = "$MsaD/MSAcat.fna";

#prep tree Options
my $tOhr = createTreeOpt($multAli,"allsites","",0,"");
my %Tree1 = %{$tOhr};
my $tOhrNSun = createTreeOpt($multAliNonSyn,"nonsyn","",0,$Tree1{nwk});
my $tOhrSyn = createTreeOpt($multAliSyn,"syn","",0,$Tree1{nwk});



#DEBUG
#mergePids("$outD/MSA/",40, "NT") ;die;
#my $tmp = "/g/bork5/hildebra/results/TEC2/v5/T2dphylo/rDNA2/fullGenomes/ini16S.fna";
#calcDisPos2($tmp,"$outD/MSA/percID_syn.txt",1); die;


my @MSAs; my @MSA_AA; my @MSAsSyn; my @MSAsNonSyn;#full MSAs and MSAs with syn / nonsyn pos only
my @MSrm; 
my %FAA ; my %FNA ; my @geneList; my @geneListF;
my $doMSA = 1;
my $treesDone = treePresent($tOhr)
	&& (!$calcNonSyn || treePresent($tOhrNSun))
	&& (!$calcSyn || treePresent($tOhrSyn));
my $calcMSA = !$treesDone && !fileGZe($multAli);
#if (!$treesDone){#cleanup, avoid checkpoints..
#	system "rm -f $treeD/*";
#}
$doMSA =0 if ($isAligned || (
			$continue && (fileGZe($multAli) || !$calcMSA) && (fileGZe($multAliSyn) ||!$calcSyn)&& (fileGZe($multAliNonSyn) ||!$calcNonSyn)) );  ## checks if MSA already exists
			
#test if MSA is gzed
if (!-e $multAli && -e "${multAli}.gz"){
	my $gunCmd = "$pigzBin -p $ncore -d ${multAli}.gz\n";
	systemW($gunCmd);
}
#die "$doMSA $calcMSA $treesDone $continue $multAli\n";
#die "$doDNDS\n";
#my @xx = keys %FAA; die "$xx[0] $xx[1]\n$FAA{HM29_COG0185}\n";
if ($isAligned){
	my $alignedInput = $useAA4tree ? $aaFna : $fnFna;
	die "-isAligned currently requires an uncompressed input file: $alignedInput\n" unless -f $alignedInput;
	if (-e $alignedInput){
		unlink $multAli if -e $multAli || -l $multAli;
		my $source = abs_path($alignedInput) or die "Cannot resolve aligned input $alignedInput: $!\n";
		symlink($source, $multAli) || copy($source, $multAli)
			or die "Cannot link or copy aligned input $source to $multAli: $!\n";
		print "Using ".($useAA4tree ? "AA" : "NT")." sequences to build tree..\n\n";
	}
} elsif (!$doMSA && $cogCats ne ""){
	fillGeneList($cogCats);
} elsif ($doMSA && $cogCats ne ""){
	my $hr = readFasta($fnFna,1); my %FNA = %{$hr};
	my %FAA;
	if ($aaFna ne ""){
		$hr = readFasta($aaFna,1); %FAA = %{$hr};
	}
	print "Read Fasta(s)\n";

	############# test length of fna sequences can be divided by 3 ##############################
		
	if($doLengthCheck){
		my @notThrees;
		my $FNAseq; my $length; my $div;
		while (($FNAseq) = each (%FNA)){
			$length = length($FNA{$FNAseq});
			$div = $length/3;
			push(@notThrees, $FNAseq) if($div =~ /\D/);
		}
		print "FNA seq not divisible by 3 (N=". scalar(@notThrees) . "): @notThrees\n" if (@notThrees);
	}

	############# test if enough seq in Sample to add to tree (avoids confusion in MSA) ##############################

	my $cnt = -1;  # line count in cats file
	my $ogrpCnt=0;#my %genCats; 
	my %totalNTs;#overall nt counts (not N)
	my %charCnts;#per gene NT counts
	my %maxNtCnt;# per gene max NTs observed
	my %meanNTcnt; #probably more sensible to use this re overpredicting gene length etc
	my %qtl90NTcnt;
	my @genesPerCat;
	my $geneTooShort = 0; #count genes with too little NTs in gene..
	my $geneTooLong = 0;
	
	my ($xI,$ST)= gzipopen($cogCats,"CogCATs phylo");
	#open my $xI,"<$cogCats" or die "Can't open cogcats $cogCats\n";
	chomp(my @linesCats = <$xI>);
	close $xI;
	#first cleanup of cat file..
	my @linesCats2; my @linesCats3;
	foreach (@linesCats){ #check first some parameters..
		$cnt++; my @spl = split /\t/;
		if (@spl && $spl[0] =~ m/^#/){shift @spl;}
		@spl = grep !/^NA$/, @spl;#remove NAs
		die "No sequence identifiers in category line ".($cnt + 1)."\n" unless @spl;
		my @spl2;
		#$genesPerCat[$cnt] = scalar(@spl) ;
		my @geneLs;
		my ($sp, $gene) = parseSeqId($spl[0], "category line ".($cnt + 1));
		foreach my $seq (@spl){### $seq = genomeX_NOGY
			($sp) = parseSeqId($seq, "category line ".($cnt + 1));
			die "can't find AA seq $seq\n" if ($aaFna ne "" && !exists ($FAA{$seq}));
			die "can't find fna seq $seq\n" if (!exists ($FNA{$seq}) && !$useAA4tree);
			#print "$MFAA{$curK}\n";			#my $ss = $FAA{$seq}; 			#filter per sequence 
			my $geneL;
			if ($useAA4tree || ! keys(%FNA)){
				my $num1 = $FAA{$seq} =~ tr/[\-Xx]//;
				$geneL = (length( $FAA{$seq})-$num1);
			} else{
				my $num1 = ($FNA{$seq} =~ tr/[\-Nn]// ) ;
				$geneL = (length( $FNA{$seq})-$num1);
			}
			#AA length
			$charCnts{$sp}{$seq} = $geneL;
			push(@geneLs, $geneL);
		}
		my $qtl = quantile(0.9,@geneLs);#values(%{$charCnts{$sp}}));
		#print "Q$qtl $gene @geneLs\n";
		$qtl90NTcnt{$gene} = $qtl;#
		foreach my $seq (@spl){
			my ($sp) = parseSeqId($seq, "category line ".($cnt + 1));
			#quantile(0.8,values(%{$charCnts{$sp}}));
			if ( $charCnts{$sp}{$seq} >= ($qtl90NTcnt{$gene}  * $ntFracGene)){
				push(@spl2, $seq);
				$geneTooLong++;
			} else {
				$geneTooShort++;
			}
		}
		push(@linesCats2,\@spl2);
		#has to work with what is actually there, not what could have been..
		$genesPerCat[$cnt] = scalar(@spl2);
		#die;
	}
	#die;
	my $GenesQtl90 = quantile(0.9,@genesPerCat);
	my $GenesQtl50 = quantile(0.5,@genesPerCat);
	$cnt=-1;
	foreach my $aRef (@linesCats2){ #remove genes with just too few genes..
		$cnt++; my @spl = @{$aRef};
		if (@spl >= (($GenesQtl90 * $fracMaxGenes90pct) ) ){ #$GenesQtl50 || 
			push(@linesCats3,\@spl);
		}
		#print @spl . " ";
	}
	
	print "\n\n-----------------  Prefilter  ------------------\n";
	print "Remaining gene cats: ". scalar(@linesCats3) . "/" . scalar(@linesCats)."; Removed $geneTooShort/$geneTooLong genes < $ntFracGene qtl90 gene length\n";
	print "Warning:: Size linesCats3:: " . @linesCats3 . " linesCats2:: " .@linesCats2 ."\n$GenesQtl50 || $GenesQtl90 * $fracMaxGenes90pct\n" if (@linesCats3 < 20);
	@linesCats2 = (); #make space..
	$cnt=-1;
	foreach my $aRef (@linesCats3){
		$cnt++; my @spl = @{$aRef};
		my ($firstSample, $gene) = parseSeqId($spl[0], "filtered category line ".($cnt + 1));
		foreach my $seq (@spl){
			my ($sp, $seqGene) = parseSeqId($seq, "filtered category line ".($cnt + 1));
			die "Wrong gene in $seq, expected $gene!\n" if ($seqGene ne $gene);
			$specList{$sp} ++;
			#my $seq2 = $seq;
			if (!exists($maxNtCnt{$gene})){
				$maxNtCnt{$gene} = $charCnts{$sp}{$seq};
			} elsif (
				$charCnts{$sp}{$seq} > $maxNtCnt{$gene}){ $maxNtCnt{$gene} = $charCnts{$sp}{$seq};
			}
			#push(@geneLgt,$charCnts{$sp}{$seq});
			$meanNTcnt{$gene}+=$charCnts{$sp}{$seq};
		}
		$meanNTcnt{$gene} /= scalar(@spl);
		# sum( ( sort { $a <=> $b } @_ )[ int( $#_/2 ), ceil( $#_/2 ) ] )/2;
		#$qtl90NTcnt{$gene} = (sort { $a <=> $b }( @geneLgt)) [ int($#geneLgt*0.9) ];
		#print "DEB: $gene $meanNTcnt{$gene} $qtl90NTcnt{$gene}\n";
		#second round, do some prefiltering already, now that we have qtl and mean gene size
		foreach my $seq (@spl){
			my ($sp) = parseSeqId($seq, "filtered category line ".($cnt + 1));
			#first check if gene gets removed
			#next if ( ($charCnts{$sp}{$seq} < $maxNtCnt{$gene} ) * $ntFracGene);
			#next if ( $charCnts{$sp}{$seq} < ($qtl90NTcnt{$gene} * $ntFracGene));
			$totalNTs{$sp} += $charCnts{$sp}{$seq};
		}
	}
	#die "@genesPerCat\n$GenesQtl90\n".$fracMaxGenes90pct*$GenesQtl90."\n" ;
	my @specs = keys %specList;
	#print "specs:: @specs\n";
	die "No species left after filtering!!\n" if (@specs == 0);
	
	
	my $maxGenes=0; my $maxNtCntTotal=0; #my @allNTcnts;
	foreach my $sp (@specs){
		if ($specList{$sp}>$maxGenes){$maxGenes = $specList{$sp};}
		if (!exists($totalNTs{$sp})){$totalNTs{$sp}=0;}
		if ($maxNtCntTotal< $totalNTs{$sp}){$maxNtCntTotal = $totalNTs{$sp};}
		#push(@allNTcnts,$totalNTs{$sp});
	}
	my $qtl90NTcntAll =  quantile(0.9,values(%totalNTs));#@allNTcnts);#(sort { $a <=> $b }( @allNTcnts)) [ int($#allNTcnts*0.9) ];
	my $qtl95Genes = quantile(0.95,values(%specList));
	my $qtl90Genes = quantile(0.9,values(%specList));

	my %smplsRmvd; my $tooFewGenes=0;my $tooFewNTs=0;my $tooFewNTs2=0; my $specsRemain = 0;
	#print "Samples removed due to low gene presence:\n";
	my $OGfnd=0;
	foreach my $sp (@specs){
		my $isOG=0;  if ($outgroup ne "" && $outgroup eq $sp){$isOG = 1;$OGfnd++;}
		
		my $NTfilter = 0; $NTfilter =1 if ( $totalNTs{$sp} < ($qtl90NTcntAll * $ntFrac));
		my $lengthInNt = (!$useAA4tree && keys(%FNA)) ? $totalNTs{$sp} : $totalNTs{$sp} * 3;
		my $NTfilter2 =  0;$NTfilter2 = 1 if ($lengthInNt < $ntCntTotal);
		my $NTlengFilt = 0; $NTlengFilt =1 if ($specList{$sp} <  ($qtl90Genes * $GeneFracPSpec) );

		if (!$isOG && ($NTlengFilt || $NTfilter || $NTfilter2) ){
			$smplsRmvd{$sp}=1;
			$tooFewNTs++ if ($NTfilter);
			$tooFewNTs2++ if ($NTfilter2);
			$tooFewGenes++ if ($NTlengFilt);
			
			#print " $sp:$specList{$sp}:$totalNTs{$sp}; ";
		} else {
			$specsRemain ++;
		}
	}
	
	if ($OGfnd == 0 && $outgroup ne ""){
		die "could not find outgroups in sequence set!\n$outgroup\n$fnFna\n";
	}
	#############################################################################################

	
	print "------------------------------------------------\n";
	print "Per species: MaxGenes: $maxGenes, Qtl90Genes: $qtl90Genes, MaxAA: $maxNtCntTotal, Qtl90 NTs: $qtl90NTcntAll\n";
	print "Species/Smpls removed: <NTs($ntFrac,$ntCntTotal):$tooFewNTs,$tooFewNTs2 ; <genes/species($GeneFracPSpec):$tooFewGenes\n";
	print "Remaining Smpls/Strains: $specsRemain/".scalar(keys%specList)."\n";
	print "------------------------------------------------\n";
	#die "$maxGenes\n";
	@linesCats = (); #empty array


	#die;
	$cnt=-1; #line counter
	foreach my $aRef (@linesCats3){#go over each gene category, building MSA for each
	#----------------- main MSA loop ----------------------
		$cnt++; my @spl = @{$aRef};
		if (@spl ==0){print "No categories in cat file line $cnt\n";next;}
		if ($spl[0] =~ m/^#/){shift @spl;}
		my @spl2 = parseSeqId($spl[0], "category line ".($cnt + 1));
		my $gene = $spl2[1];
		my $gene_file_stem = geneFileStem($gene);
		#die "@spl\n";		
		my $ogrGenes = "";
		if ($outgroup ne ""){
			foreach my $seq (@spl){
				my ($seqSample) = parseSeqId($seq, "category line ".($cnt + 1));
				if ($seqSample eq $outgroup){$ogrGenes = $seq; $ogrpCnt ++ ;last;}
			}
		}
		
		#die "$spl2[0]\t$spl2[1]\n";
		

		die "Double gene name in tree build pre-concat: $spl2[0] $spl2[1]\n" if (exists($usedGeneNms{$spl2[1]}));
		$usedGeneNms{$spl2[1]} = 1;
		#die "@spl\n";
		my $tmpInMSA = "$tmpD/inMSA$cnt.faa";
		my $tmpInMSAnt = "$tmpD/inMSA$cnt.fna";
		my $tmpOutMSAaa = "$tmpD/$gene_file_stem.$cnt.faa";
		my $tmpOutMSA = "$tmpD/$gene_file_stem.$cnt.fna";
		my $finOutMSAaa = "$MsaD/$gene_file_stem.$cnt.faa";
		my $finOutMSA = "$MsaD/$gene_file_stem.$cnt.fna";
		
		my $endFileExists=0; $endFileExists =1 if (fileGZs($finOutMSAaa) && fileGZs($finOutMSA));
		
		open O,">$tmpInMSA" or die "Can;t open tmp faa file for MSA: $tmpInMSA\n";
		open O2,">$tmpInMSAnt" or die "Can;t open tmp fna file for MSA: $tmpInMSAnt\n";
		my $seqType = "AA";my $seqTypeOth = "NT";
		my $seqLength = 0; my $numSeq =0;
		
		#1st: collate sequences
		#do here already per gene length check .. probably better for alignment
		foreach my $seq (@spl){### $seq = genomeX_NOGY
			my ($sp) = parseSeqId($seq, "category line ".($cnt + 1));
			next if (exists($smplsRmvd{$sp}));
			if ($specList{$sp} <  ($qtl90Genes * $GeneFracPSpec) ){die "buildTree: GeneFracPSpec maxGenes shouldn't be here!\n";}
			my $seq2 = $seq;
			#just for this singular case applying..
			#next if ( $charCnts{$sp}{$seq} < ($qtl90NTcnt{$gene}  * $ntFracGene));  #maxNtCnt{$gene}
			if (!$smplDef){#create artificial head tag
				#TODO.. don't need it now for tec2, since no good NCBI taxid currently...
			}
			$samples{$sp} = 1; #$genCats{$spl2[1]} = 1; 
			$FAA{$seq} =~ s/\*$//g if ($MSAprog != 0);
			$FAA{$seq} =~ s/\x00//g;

			print O ">$seq2\n$FAA{$seq}\n";
			if (!$useAA4tree){
				$FNA{$seq} =~ s/\x00//g;
				$FNA{$seq} =~ s/-//g;
				print O2 ">$seq2\n$FNA{$seq}\n";
				$seqLength += length($FNA{$seq});
			} else {
				$seqLength += length($FAA{$seq});
			}
			$numSeq++;
		}
		close O;close O2;
		#done, samples are in O2
		if ($numSeq < 3){ #three tips are sufficient for the minimal resolved unrooted tree
			#system "rm -f $tmpInMSA $tmpInMSAnt";#
			unlink  $tmpInMSA; unlink $tmpInMSAnt;next;
		}

		$seqLength /= $numSeq;
		#print "$tmpInMSA,$tmpOutMSA2\n";
		my $cmd1 = MSA($tmpInMSA,$tmpOutMSAaa,$ncore,$MSAprog,$numSeq);
		#zorro/macse filter.. by default not used ("")
		my $cmd2 = filterMSA($tmpInMSA,$tmpOutMSAaa,$ncore,$postFilter,$useAA4tree);
		
		if ($endFileExists){$cmd1="";$cmd2="";}
		#if ($tmpOutMSA2 eq ""){ #filter failed...
		#	print "skipped protein, too many bad positions\n";
		#	system "rm -f $tmpOutMSAaa $tmpOutMSA $tmpInMSA $tmpInMSAnt";
		#	next;
		#}
		#$cmdGrand .= $cmd1."\n".$cmd2."\n";
		#print "$cmd1\n$cmd2\n";
		systemW($cmd1."\n".$cmd2."\n");
		die "MSA command completed without producing $tmpOutMSAaa\n" unless -s $tmpOutMSAaa;

		
		
		#move into its own for loop
		#dist mat related
		my $tmpDMatOth = "$MsaD/${seqTypeOth}_clustalo_percID_${cnt}_".int($seqLength).".txt";
		my $tmpDMat = "$MsaD/${seqType}_clustalo_percID_${cnt}_".int($seqLength).".txt";
		my $inFastaOth = $tmpInMSAnt;
		my $percIDhr; my $avgID; my $pIDsmplhr;
		if ($calcDistMat){ #for dmat: calc each gene spearately and merge scores later
			($avgID,$pIDsmplhr,$percIDhr)  = calcDisPos2($tmpInMSA,$tmpDMat,0,$ncore,$tmpD);
		#			if ($calcDistMatExt && -e $inFastaOth){
			if (-e $inFastaOth){ 
				($avgID,$pIDsmplhr,$percIDhr) = calcDisPos2($inFastaOth,$tmpDMatOth,1,$ncore,$tmpD);
				$calcDistMatExtGo = $avgID;
			} else { $calcDistMatExtGo = 0;}
		}
		
		
		
		#die "@MSAs\n";
		if (!$useAA4tree){
			#this part now is all concerned about NT level things..
			my ($tmpOutMSAsyn,$tmpOutMSAnonsyn) = ($tmpOutMSA, $tmpOutMSA);
			$tmpOutMSAnonsyn =~ s/\.fna/\.nonsyn\.fna/;$tmpOutMSAsyn =~ s/\.fna/\.syn\.fna/;

			if (!$endFileExists){
				convertMultAli2NT($tmpOutMSAaa,$tmpInMSAnt,$tmpOutMSA);
				($tmpOutMSAsyn,$tmpOutMSAnonsyn) = synPosOnly($tmpOutMSA,$tmpOutMSAaa,0,$ogrGenes,$calcSyn,$calcNonSyn);
				#this will not affect 4-fold only etc..
				my $msaFbin = getProgPaths("MSAfix");
				$cmd = "$msaFbin -i $tmpOutMSA  -maskLowID -maskBorderGap -rmGapColsGreater ".$maxGapPerCol." -minGoodPosFrac 0.6\n";
				systemW $cmd;
			}
			push (@MSAs,$finOutMSA);
			push (@MSAsSyn,$tmpOutMSAsyn) if ($tmpOutMSAsyn ne "");
			push (@MSAsNonSyn,$tmpOutMSAnonsyn) if ($tmpOutMSAnonsyn ne "");
			#die "@MSAs\n";
		} else {
			push (@MSA_AA,$finOutMSAaa);
		}
		#system "rm -f $tmpInMSA $tmpInMSAnt";# $tmpOutMSAaa";
		unlink  $tmpInMSA; unlink $tmpInMSAnt;
		push (@MSrm,$finOutMSAaa,$finOutMSA);
		#die "$MSrm[1]\n";
		print "$cnt "; 
		move($tmpOutMSAaa, $finOutMSAaa) or die "Cannot move $tmpOutMSAaa to $finOutMSAaa: $!\n"
			if (!fileGZs($finOutMSAaa) && -e $tmpOutMSAaa);
		move($tmpOutMSA, $finOutMSA) or die "Cannot move $tmpOutMSA to $finOutMSA: $!\n"
			if (!fileGZs($finOutMSA) && -e $tmpOutMSA);
	}
	
	my $mergPIDtag = "_merge";
	mergePids("$MsaD/",$cnt, "AA",$mergPIDtag) if ($calcDistMat); #merge different percIDs
	mergePids("$MsaD/",$cnt, "NT",$mergPIDtag) if ($calcDistMat); #merge different percIDs
	
	#could be used to filter genes further, but not for now
	#pogenStatsFilter();
	if ($outgroup ne ""){print "Found $ogrpCnt of $cnt outgroup sequences\n";}
} elsif ($doMSA) {#no marker way, single gene
	my $r1; my $r2;
	#,$r1,$r2)
	$multAli = singleGeneMSAprocess($multAli)#;,\@MSAs,\@MSA_AA);
	#@MSAs = @{$r1};	@MSA_AA = @{$r2};
}

#die "@MSA_AA\n\n";
if ($calcMSA && $cogCats ne "" && @MSAs == 0 && @MSA_AA == 0 && !fileGZs($multAli)){
	die "No usable MSAs were generated; no tree was created\n";
}
if ($calcMSA && $cogCats eq "" && !fileGZs($multAli)){
	die "Single-gene alignment was not generated: $multAli\n";
}

#prep final MSA file that is correct NT or AA and is merged
if (!$useAA4tree) {
	if ($cogCats eq ""){ #single gene case
		my ($hr,$OK) = readFasta($multAli,1); writeFasta($hr,$multAli);#complicated way to shorted headers of infile
	}
	mergeMSAs(\@MSAs,\%samples,$multAli,0,0);
	mergeMSAs(\@MSAsSyn,\%samples,$multAliSyn,1,0) if ($calcSyn);
	mergeMSAs(\@MSAsNonSyn,\%samples,$multAliNonSyn,1,0) if ($calcNonSyn);
	@theRealMSAs = @MSAs;

} else {#useAA4tree
	mergeMSAs(\@MSA_AA,\%samples,$multAli,0,1); #sames files as in @MSrm
	@theRealMSAs = @MSA_AA;
}

#phylip conversion??
if ( $doGenesToPh){ 
	my $phylipD = File::Spec->catdir($outD, "phylip");
	make_path($phylipD) unless -d $phylipD;
	my $fasta2phylip = getProgPaths("fasta2phylip_scr");

	foreach my $MSAfn (@MSAs){
		my $phylipOut = File::Spec->catfile($phylipD, basename($MSAfn).".ph");
		unlink $_ or die "Cannot remove stale PHYLIP output $_: $!\n" for glob("$phylipOut*");
		my $cmd2 = "$fasta2phylip -c 50 ".shellQuote($MSAfn)." > ".shellQuote($phylipOut)."\n";
		systemW $cmd2;
		die "PHYLIP conversion did not produce $phylipOut\n" unless -s $phylipOut;
		push(@geneList, $phylipOut);
	}
}


#die;

#-------------------------------------------
#Tree prep phase (MSA clean up, conversion, 4fold sites etc)
#convert fasta again




#-------------------------------------------
#Supertrees, gubbins etc
#-------------------------------------------

my $phyloTree = "";
if ($doSuperTree || $doSuperCheck){#can be for 2 reasons: 1) build actual super tree 2) quality control
	my @treeCol;
	for (my $i=0;$i<@theRealMSAs;$i++){
		print "===============>  Subtree $i  <===============\n";
		my $tOhrST = createTreeOpt($theRealMSAs[$i],"allsites",$i,1,"");
		my $trRetH = treeAtHeart($tOhrST);
		push(@treeCol,${$trRetH}{IQtreeout}.".treefile");
		if ($calcSyn){
		} 
		if ($calcNonSyn){
		}
	}
	if ($doSuperTree){
		my $outST = "$treeD/IQtree_allsites.treefile";
		my $specFile = "$treeD/IQtree_allsites.species";
		open OU,">$specFile" or die "Cannot write supertree species file $specFile: $!\n";
		print OU join "\n",sort keys (%specList);
		close OU;
		my $stBin = getProgPaths("supertree",0);
		my $cmd = "$stBin -s $specFile -F -o - @treeCol  | grep '\\[F01\\]' | cut -f2 -d' ' > $outST"; #-F
		#die "ST:\n$cmd\n";
		systemW $cmd;
		die "Supertree command did not produce $outST\n" unless -s $outST;
		$phyloTree = $outST;
	} elsif ($doSuperCheck) {
		
	}
}


if ($calcDNAdiff){
	calcDiffDNA($multAli,"$MsaD/percID_2.txt");
}


if ($doGubbins){
	$gubbinsBin = requireConfiguredTool("MF4_GUBBINS_BIN", "Gubbins") if $gubbinsBin eq "";
	my $gubbinsOutDir = File::Spec->catdir($outD, "gubbins");
	make_path($gubbinsOutDir) unless -d $gubbinsOutDir;
	my $outDG = File::Spec->catfile($gubbinsOutDir, "GD");
	if ($continue && -e $outDG.".final_tree.tre"&& -e $outDG.".summary_of_snp_distribution.vcf"){
		print "Gubbins result already exists in output folder, run will be skipped\n";
	} else {
		unlink $outDG if -e $outDG;
		my $cmdG = "$gubbinsBin --filter_percentage 50  --tree_builder hybrid --prefix $outDG --threads $ncore $multAli";
		if (0){$cmdG.=" --outgroup $outgroup";}
		systemW $cmdG;
		#die $cmdG."\n";
		print "Gubbins run finished\n";
	}
}


#distamce matrix, this is fast
#system "$trimalBin -in $multAli -gt 0.1 -cons 100 -out /dev/null -sident 2> /dev/null > $outD/MSA/percID.txt\n";
if (0 && !$useAA4tree){ #this is outdated
	calcDisPos($multAli,"$MsaD/percID.txt",1) unless(-e "$MsaD/percID.txt" && $continue);
	if ($calcSyn){
	#		system "$trimalBin -in $multAliSyn -gt 0.1 -cons 100 -out /dev/null -sident 2> /dev/null > $outD/MSA/percID_syn.txt\n";
	#	calcDisPos($multAliSyn,"$outD/MSA/percID_syn.txt",1);
	}
	if ($calcNonSyn){
		#system "$trimalBin -in $multAliNonSyn -gt 0.1 -cons 100 -out /dev/null -sident 2> /dev/null > $outD/MSA/percID_nonsyn.txt\n";
	#	calcDisPos($multAliNonSyn,"$outD/MSA/percID_nonsyn.txt",1);
	}
} elsif (0) {
	calcDisPos($multAli,"$MsaD/AA_percID.txt",1) unless(-e "$MsaD/percID.txt" && $continue);
}

#print "$theRealMSAs[0]\n\n\n";

#-------------------------------------------
#Tree building part with RaxML, IQtree, fasttree2, phyml
#-------------------------------------------
die "Expected a non-empty merged alignment before tree construction: $multAli\n"
	if $MSAreq && !fileGZs($multAli);

my $trRetH;
if ($doSuperTree){
	$Tree1{nwk} = $phyloTree;
	$trRetH = \%Tree1;
} else {
	$trRetH = treeAtHeart($tOhr);
	if ($calcSyn){ #tree at syn pos
		treeAtHeart($tOhrSyn);
	}
	if ($calcNonSyn){ #tree at non-syn pos
		treeAtHeart($tOhrNSun);
	}
}
#system "rm -f $multAli.ph $multAliSyn.ph $multAliNonSyn.ph";

if ($useTreeShrink){
	my $trShr = getProgPaths("treeshrink");
	my $inputTree = ${$trRetH}{nwk} // "";
	die "TreeShrink requested but no completed tree is available\n" unless $inputTree ne "" && -s $inputTree;
	my $cmd = "$trShr -i $outD -t ".shellQuote($inputTree)." -q 0.05  -O TS. -f";
	print $cmd."\n";
	systemW($cmd);
}

#die "$distTree_scr -d -a --dist-output $raxD/distance.syn.txt $raxD/RXML_sym.nwk\n";


#### post tree building #######
#### compute dNdS and other popgen values ######

pogenStatsFilter();

#might require $doGenesToPh for paml??
if($doDNDS){
	make_path($codemlOutD) unless -d $codemlOutD;
#die "@geneList";
	if($selGene){@geneList=@genesExtra;	}
	#my $tmpDir = $codemlOutD."/tmp/";;
	make_path($tmpD) unless -d $tmpD;
	$treeFile = $Tree1{nwk} if ($treeFile eq "");
	selecAnalysis(\@geneList, $treeFile, $codemlOutD, $tmpD);   

}
#if ($doTheta){ #Watterman estimator (Theta)
#	if($selGene){@geneList=@genesExtra;	}
#	WattTheta(\@geneList,$MsaD,$codemlOutD);
#}

FastGear();
if ($removeMSA){
	safeRemoveTree($MsaD, $tmpD);
} elsif ($gzipInput){
	for my $msaFile (grep { -f $_ && $_ !~ /\.gz$/ } glob(File::Spec->catfile($MsaD, "*"))){
		systemW("$pigzBin -p $ncore ".shellQuote($msaFile));
	}
}
if ($gzipInput){
	for my $inputFile ($aaFna, $fnFna, $cogCats){
		next if $inputFile eq "" || $inputFile =~ /\.gz$/ || !-f $inputFile;
		systemW("$pigzBin -p $ncore ".shellQuote($inputFile));
	}
}

safeRemoveTree($tmpD, $tmpBase);
	###################### ETE ######################3

print "All done: $outD \n\n";
exit(0);













##########################################################################################
##########################################################################################





sub treePresent{
	my ($hr) = @_;
	my %treeOpts = %{$hr};
	my $ret = 1;
	my $checked = 0;
	if ($doFastTree){
		$checked = 1;
		$ret=0 unless ($continue && -s $treeOpts{fastTrOut});
	}
	if ($doVeryFastTree){
		$checked = 1;
		$ret=0 unless ($continue && -s $treeOpts{VfastTrOut});
	}
	if ($doIQTree){
		$checked = 1;
		my $IQtree = "$treeOpts{IQtreeout}";
		$ret=0 unless ($continue && -s "$IQtree.treefile");
	}
	if ($doRAXMLng){
		$checked = 1;
		$ret=0 unless ($continue && -s $treeOpts{RAXNGtreeout});
	}
	if ($doRAXML){
		$checked = 1;
		$ret=0 unless ($continue && -s $treeOpts{RAXtreeout});
	}
	return $checked ? $ret : 0;
}



sub createTreeOpt{
	my ($multF,$siteTag,$tcnt,$silent,$consTree) = @_;
	#$siteTag="allsites";
	my $isSubTree = 0;
	my $outgroupL = $outgroup;
	$isSubTree = 1 if ($tcnt ne "");
	$outgroupL = "" if ($isSubTree);
	my $partiF=$multF.$partiExt;
	if (-e "$partiF.gz"){systemW("$pigzBin -d ".shellQuote("$partiF.gz"));}
	# Keep the expected path even on a fresh run: mergeMSAs creates this file
	# after the tree options are assembled.  Its existence is resolved only
	# immediately before a tree program is invoked.
	#object to transfer options to tree (and get them back..)
	my $BStag = ""; if ($bootStrap>0){$BStag="_BS$bootStrap";}
	my %treeOpts = (inMSA => $multF,
					IQtreeout => "$treeD/IQtree${tcnt}_${siteTag}",
					ncore => $ncore,
					tcnt => $tcnt,
					outgr => $outgroupL,
					bootStrap => $bootStrap,
					useAA => $useAA4tree,
					iqtreeFast => $iqFast,
					autoModel => $treeAutoModel,
					cont => $continue,
					silent => $silent,
					partition => $partiF,
					constraintTree => $consTree,
					isSubTree => $isSubTree,
					PhymTree => "$treeD/phyml${tcnt}_${siteTag}.nwk",
					fastTrOut => "$treeD/FASTTREE${tcnt}_${siteTag}.nwk",
					VfastTrOut => "$treeD/VERYFASTTREE${tcnt}_${siteTag}.nwk",
					RAXtreeout => "$treeD/RXML${tcnt}_${siteTag}$BStag.nwk",
					RAXNGtreeout => "$treeD/RXng${tcnt}_${siteTag}$BStag.nwk",
					runSafe => 0,
					);
	return \%treeOpts;
}

#core routine to calculte (start) phylo reconstruction
sub treeAtHeart{
	my ($hr) = @_;
	my %treeOpts = %{$hr};
	my $consTree = $treeOpts{constraintTree}; my $multF = $treeOpts{inMSA};
	my $silent = $treeOpts{silent}; my $tcnt = $treeOpts{tcnt};
	my $partition = $treeOpts{partition} // "";
	$treeOpts{partition} = "" unless $partition ne "" && -s $partition;
	
	if ($consTree ne ""){
		die "ref tree $consTree does not exist\n" unless (-e $consTree);
		
		my $hr = getTreeLeafs($consTree);
		my %nwLfs= %{$hr};
		
		my $cntMissTree=0;
		$hr = readFasta($multF,1,"input MSA to check for constraint tree");
		my %FNA = %{$hr};
		my $unpresent=0;
		foreach my $ge (keys %FNA){
			unless (exists($nwLfs{$ge})){$unpresent++; last;}
		}
		my @genomeList;
		if ($unpresent){
			my $cntMissTree=0;
			open O,">$multF.tmp" or die "can't open MSA out $multF.tmp\n";
			foreach my $genome (keys %FNA){
				#print "$genome\n";
				unless (exists($nwLfs{$genome})){$cntMissTree++; next;}
				my $seq = $FNA{$genome};
				print O ">$genome\n$seq\n";
				push (@genomeList, $genome);
			} 
			close O;
			if ($cntMissTree>0){print "Removed $cntMissTree extra sequences from MSA not in constraint tree\n";}
			move("$multF.tmp", $multF) or die "Cannot replace constrained MSA $multF: $!\n";
		} else {
			@genomeList = keys %FNA;
		}
		
		#my $ar = readFastHD($multF);
		my $nwkPrune = $multF."pruned.nwk";
		if (scalar(@genomeList) > 3){
			pruneTree($consTree,\@genomeList,$nwkPrune);
			$treeOpts{constraintTree} = $nwkPrune;
		} else {
			$treeOpts{constraintTree} = "";
		} 

	}
	#print "cons: $treeOpts{constraintTree}\n";
	#convert MSA to NEXUS
	#convertMSA2NXS($multAli,"$multAli.nxs");
	#format conversion for raxml..
	if ($doRAXML){
		my $f = $treeOpts{inMSA};
		my $fasta2phylip = getProgPaths("fasta2phylip_scr");
		my $tcmd = "rm -f $f.ph*; $fasta2phylip -c 50 $f > $f.ph\n";
		systemW $tcmd;#) {die "fasta2phylim failed:\n$tcmd\n";}
	}


	if ($doFastTree){
		unless ($continue && -s $treeOpts{fastTrOut}){
			runFasttree($treeOpts{inMSA},$treeOpts{fastTrOut},$treeOpts{useAA},$treeOpts{ncore});
		}
		$phyloTree = $treeOpts{fastTrOut};
	}
	if ($doVeryFastTree){
		unless ($continue && -s $treeOpts{VfastTrOut}){
			runVeryFasttree($treeOpts{inMSA},$treeOpts{VfastTrOut},$treeOpts{useAA},$treeOpts{ncore});
		}
		$phyloTree = $treeOpts{VfastTrOut};
	}
	if ($doIQTree){
		my $IQtree = "$treeOpts{IQtreeout}";
		unless ($continue && -s "$IQtree.treefile"){
			runQItree(\%treeOpts);
		}
		$phyloTree = "$IQtree.treefile";
	}
	if ($doRAXMLng){
		unless ($continue && -s $treeOpts{RAXNGtreeout}){
			runRaxMLng(\%treeOpts);
		}
		$phyloTree = $treeOpts{RAXNGtreeout};
	}
	if ($doRAXML){
		$treeOpts{inMSA} = "$multF.ph";
		if (!-s $treeOpts{inMSA}){ die "Can't find nonempty expected *.ph file: $multF.ph";}
		# runRaxML's continuation logic historically keys on path existence.  Do
		# not let a zero-byte published tree suppress recovery of partial work.
		unlink $treeOpts{RAXtreeout}
			or die "Cannot remove empty RAxML tree $treeOpts{RAXtreeout}: $!\n"
			if -e $treeOpts{RAXtreeout} && !-s $treeOpts{RAXtreeout};
		runRaxML(\%treeOpts);#"$multF.ph",$bootStrap,$outgroup,"$treeD/RXML_allsites$BStag.nwk",$ncore,$continue,!$useAA4tree);
		$phyloTree = $treeOpts{RAXtreeout};#"$treeD/RXML_$siteTag$BStag.nwk";
	}
	if ($doCFML && !$treeOpts{isSubTree}){
		my $outDG = File::Spec->catdir($outD, "clonalFrameML");
		make_path($outDG) unless -d $outDG;
		$outDG = File::Spec->catfile($outDG, "CFML");
		my $CFMLbin = requireConfiguredTool("MF4_CLONALFRAMEML_BIN", "ClonalFrameML");
		my $cmd = "$CFMLbin $phyloTree $multF $outDG\n";
		systemW($cmd);
		die "ClonalFrameML did not produce its expected labelled-tree output for $outDG\n"
			unless glob("${outDG}*labelled_tree.newick");
	}

	#phyml

	if ($doPhym){
		die "Phym is outdated, check code to reactivate..\n";
		my @thrs;
		my $tcmd = "";
		my $nwkFile = $treeOpts{PhymTree};#"$treeD/phyml${tcnt}_${siteTag}.nwk";
		my $phymlBin = getProgPaths("phyml");

		$tcmd = "$phymlBin --quiet -m GTR --no_memory_check -d nt -f m -v e -o tlr --nclasses 4 -b 2 -a e -i $multF.ph > $nwkFile\n";
		push(@thrs, threads->create(sub{system $tcmd;}));
		for (my $t=0;$t<@thrs;$t++){
			my $state = $thrs[$t]->join();
			if ($state){die "Thread $t exited with state $state\nSomething went wrong with RaxML\n";}
		}
	}
	$treeOpts{nwk} = $phyloTree;
	return (\%treeOpts);
}



sub singleGeneMSAprocess($){
#$multAli,\@MSAs,\@MSA_AA
	my ($multAli) = @_;#,$MSAsar,$MSA_AAar) = @_;
	#my @MSAs = @{$MSAsar};
	#my @MSA_AA=@{$MSA_AAar};
	print "No gene categories given, assumming 1 gene / species in input\n";
	my $tmpInMSA = $aaFna;
	#my $tmpInMSAnt = $fnFna;
	my $tmpOutMSAaa = "$tmpD/outMSA.faa";
	#my $tmpOutMSAsyn = $multAliSyn;#"$tmpD/outMSA.syn.fna";
	#my $tmpOutMSAnonsyn = $multAliNonSyn;
	
	#print "$tmpInMSAnt\n";
	my $numFas;
	if (-e $aaFna){ $numFas = `grep -c '^>' $aaFna`;#just faster to use aa file..
	} else {$numFas = `grep -c '^>' $fnFna`;
	}
	chomp $numFas;
	#die "$fnFna $numFas\n"; 
	my $seqType = "AA";
	my $seqTypeOth = "NT";
	my $inFasta = $aaFna;
	my $inFastaOth = $fnFna;
	my $mainTypeIsAA=1;
	if ($aaFna eq ""){ #always choose AA as default alignment, nt is only fallback
		$inFasta = $fnFna ;		$inFastaOth = $aaFna;
		$seqType = "NT";$seqTypeOth = "AA";	$mainTypeIsAA = 0;
	}
	die "Not enough sequences to build an alignment (found $numFas)\n" if $numFas <= 1;
	if ($isAligned){
		copy($inFasta, $tmpOutMSAaa) or die "Cannot copy aligned input $inFasta to $tmpOutMSAaa: $!\n";
	} else{
		#MSA calculation
		my $msaCmd = MSA($inFasta,$tmpOutMSAaa,$ncore,$MSAprog,$numFas);
		systemW($msaCmd);
		die "Single-gene MSA command completed without producing $tmpOutMSAaa\n" unless -s $tmpOutMSAaa;
	}

	if ($calcDistMat){
		if (-e $tmpInMSA){
			calcDisPos2($tmpInMSA,"$outD/MSA/${seqType}_vsearch_percID.txt",!$mainTypeIsAA,$ncore,$tmpD);
		}
		if (-e $inFastaOth){
			calcDisPos2($inFastaOth,"$outD/MSA/${seqTypeOth}_vsearch_percID.txt",$mainTypeIsAA,$ncore,$tmpD);
		}
	}

	if ($tmpInMSA ne "" && !$useAA4tree){
		convertMultAli2NT($tmpOutMSAaa,$fnFna,$multAli);
		my $msaFbin = getProgPaths("MSAfix");
		$cmd = "$msaFbin -i $multAli  -maskLowID -maskBorderGap -rmGapColsGreater ".$maxGapPerCol." -minGoodPosFrac 0.6\n";
		systemW $cmd;
		($multAliSyn, $multAliNonSyn) = synPosOnly($multAli,$tmpOutMSAaa,0,"",$calcSyn,$calcNonSyn);

		#system "rm $tmpInMSA $fnFna $tmpOutMSAaa";
		unlink $tmpOutMSAaa or die "Cannot remove temporary alignment $tmpOutMSAaa: $!\n" if -e $tmpOutMSAaa;
	} else {
		move($tmpOutMSAaa, $multAli) or die "Cannot move $tmpOutMSAaa to $multAli: $!\n";
	}
	#$multAli = $tmpOutMSA; $multAliSyn = $tmpOutMSAsyn;
	
	print "OUTPUT:: $multAli\n";
	
	return $multAli;#,\@MSAs,\@MSA_AA);

}

#### fastgear ##

sub FastGear{
	my ($fastgearSummaryBin, $fastgearReorderBin, $matlabBin);
	if ($doFastGearSummary){
		$fastgearSummaryBin = requireConfiguredTool("MF4_FASTGEAR_SUMMARY_BIN", "fastGEAR summary");
		$fastgearReorderBin = requireConfiguredTool("MF4_FASTGEAR_REORDER_BIN", "fastGEAR reorder");
		$matlabBin = requireConfiguredTool("MF4_FASTGEAR_MATLAB_BIN", "fastGEAR MATLAB runtime");
	}

	if($doFastGear){
#		open I,"<$cogCats" or die "Can't open cogcats $cogCats\n"; 
		my ($xI,$ST)= gzipopen($cogCats,"CogCATs phylo");

		while (<$xI>){
			my $cnt3 =0;
			chomp; my @splF = split /\t/;
			@splF = grep !/^NA$/, @splF;#remove NAs
			if (@splF ==0){print "No categories in cat file line $cnt3\n";next;}
			my @splF2 = parseSeqId($splF[0], "fastGEAR category line ".($cnt3 + 1));
			#die "@spl\n";
			push(@geneListF, $splF2[1]);
			$cnt3 ++;
		}
		close $xI;

		my $MsaDF1 = "$outD/MSA";
		my $MsaDF2 = "$outD_clust/MSA_FG";
		make_path($outD_clust);
		systemW("cp -r ".shellQuote($MsaDF1)." ".shellQuote($MsaDF2));
		

		foreach my $geneF (@geneListF){
			my $gene_file_stem = geneFileStem($geneF);
			my $outFG = "$outD/fastGear/fastGear_Results/$gene_file_stem";
			make_path($outFG) unless -d $outFG;
			my $outFileFG = "$outFG/${gene_file_stem}_res.mat";
			my @gene_msas = sort glob("$MsaDF2/$gene_file_stem.*.fna");
			die "No MSA files found for locus $geneF in $MsaDF2\n" unless @gene_msas;
			my $fastgear_input = "$MsaDF2/$gene_file_stem.fna";
			open my $fg_out, '>', $fastgear_input
				or die "Cannot create fastGEAR input $fastgear_input: $!\n";
			for my $msa (@gene_msas) {
				open my $fg_in, '<', $msa or die "Cannot read fastGEAR source MSA $msa: $!\n";
				while (my $line = <$fg_in>) {
					if ($line =~ /^>(\S+)/) {
						my ($sample) = parseSeqId($1, "fastGEAR MSA header in $msa");
						$line = ">$sample\n";
					}
					print {$fg_out} $line or die "Cannot write fastGEAR input $fastgear_input: $!\n";
				}
				close $fg_in or die "Cannot close fastGEAR source MSA $msa: $!\n";
			}
			close $fg_out or die "Cannot close fastGEAR input $fastgear_input: $!\n";
			my $FGparFile = requireConfiguredTool("MF4_FASTGEAR_PARAM_FILE", "fastGEAR parameter file");
			runFastgear($gene_file_stem, $outFileFG, $MsaDF2, $FGparFile);
		}
		safeRemoveTree($outD_clust, $outD);
		#die;
	}
		


	### postprocessing fastgear output ##
		
	if($doFastGearSummary){
		my $FGDataD = "$outD/fastGear/";
		my $summaryD = "$FGDataD/fastGear_Summaries";
		make_path($summaryD) unless -d $summaryD;
		my $resultD = "$FGDataD/fastGear_Results";
		die "no fastgear results found\n" unless(-d $resultD);
			
		#die "$treeFileFG\n";
		my $allNamesFile ="$summaryD/allNamesFromTop.txt";
		my $treeNamesFile ="$summaryD/subtreeNamesFromTop.txt";
		
		# get a list with all genomes in tree
		if(! -e $allNamesFile |! -e $treeNamesFile){
			my @genomeListFG;		
			my @FNAheader_all = `grep '^>' $multAli`;
			#die "@FNAheader_all\n";
			foreach my $genome (@FNAheader_all){
				my ($genome2) = $genome =~ m/>(.*)?/;
				push (@genomeListFG, $genome2);
				}
			open T1,">$allNamesFile";
			print T1 join("\n",@genomeListFG);
			close T1;	
			#die;
			open T2,">$treeNamesFile";
			print T2 join("\n",@genomeListFG);
			close T2;
			#die "@genomeListFG\n";
		}
		# reorder files -> required for further steps?
		my $reorder_cmd = "$fastgearReorderBin $matlabBin $FGDataD fastGear_ allNamesFromTop.txt both";
		systemW($reorder_cmd);
		
		## collect Recombination statistics #

		my $SRC_cmd = "$fastgearSummaryBin $matlabBin $FGDataD fastGear_";
		systemW($SRC_cmd);
		print "fastgear collect recombination statistics finished";
		my $FG_sumOut = "$summaryD/fastGear__recSummaries.txt";
		if(-e $FG_sumOut){safeRemoveTree($resultD, $FGDataD);}
		#die;
	}

}

sub mergePids($ $ $ $){
	my ($dir,$unusedMax,$seqType,$tag) = @_;
	#$outD/MSA/${seqType}_clustalo_percID_$cnt.txt
	#my $seqType = "AA";  
	opendir(my $dirHandle, $dir) or die "Cannot open distance-matrix directory $dir: $!\n";
	my @subfls = grep { /^\Q$seqType\E.*_percID_\d+_\d+\.txt$/ } readdir($dirHandle);
	closedir($dirHandle);
	@subfls = sort {
		my ($ai) = $a =~ /_percID_(\d+)_/;
		my ($bi) = $b =~ /_percID_(\d+)_/;
		$ai <=> $bi;
	} @subfls;
	return if (@subfls == 0);
	#die "@subfls\n$dir\n";
	my %bigMat; my %bigCnt;
	for my $disM (@subfls){
		my ($curL) = $disM =~ /_(\d+)\.txt$/;
		die "Cannot determine alignment length from distance matrix $disM\n" unless defined $curL;
		my $disPath = File::Spec->catfile($dir, $disM);
		die "can't find distance matrix $disPath\n" unless -e $disPath;
		my $cc=-2;
		my @IDS; my %cL;
		#print "$disM\n";
		open my $matrixIn, "<", $disPath or die "Cannot open distance matrix $disPath: $!\n";
		while (my $line = <$matrixIn>){
			$cc++;
			next if ($cc==-1);
			chomp $line;
			my @spl = split /\s+/,$line;
			my $id = shift @spl;
			push (@IDS,$id);
			$cL{$cc} = \@spl;
		}
		close $matrixIn;
		unlink $disPath or die "Cannot remove merged distance matrix $disPath: $!\n";
		#matrix in mem, now relate to actual dist matrix
		for (my $j=0;$j<@IDS;$j++){
			my ($id1) = parseSeqId($IDS[$j], "distance-matrix identifier");
			for (my $k=$j;$k<@IDS;$k++){
				my $curPID = $cL{$j}[$k];
				my ($id2) = parseSeqId($IDS[$k], "distance-matrix identifier");
				#print "$id1 $id2 $curPID\n";;
				$bigMat{$id1}{$id2} += $curPID * $curL;
				$bigCnt{$id1}{$id2} += $curL;
			}
		}
		#die;
	}
	my $oMat = "$dir/${seqType}${tag}_percID.txt";
	open O,">$oMat" or die "Can't open out dis mat $oMat\n";
	my @IDS = sort(keys %bigMat);
	print O "\t".join("\t",@IDS)."\n";
	for (my $j=0;$j<@IDS;$j++){
		my $id1 = $IDS[$j];
		print O $id1;
		for (my $k=0;$k<@IDS;$k++){
			my $id2 = $IDS[$k];
			if (exists($bigMat{$id1}{$id2} )){
				print O "\t".$bigMat{$id1}{$id2} / $bigCnt{$id1}{$id2};
			} elsif (exists($bigMat{$id2}{$id1} ) ) {
				print O "\t".$bigMat{$id2}{$id1} / $bigCnt{$id2}{$id1};
			} else {
				print O "\tNA"
				#die "$id1 $id2\n";
			}
		}
		print O "\n";
	}
	
	close O;
	#die $oMat;
}




#just get positions different between alignments, and their relative position
sub calcDiffDNA($ $){
	my ($MSA,$opID) = @_;
	my $kr = readFasta($MSA);
	my %MS = %{$kr};
	my $isNT = 1;
	my %diffArs;my %perID;
	print "Calculating distance matrix..\n";
	foreach my $k1 (keys %MS){
		$MS{$k1} = uc ($MS{$k1});
	}
	foreach my $k1 (keys %MS){
		my $ss1 = $MS{$k1};
		foreach my $k2 (keys %MS){
			next if ($k2 eq $k1);
			my $ss2 = $MS{$k2};
			my $mask = $ss1 ^ $ss2;
			my $diff=0;
			my$N2=($ss2 =~ tr/[-]//);
			my$N1=($ss1 =~ tr/[-]//);
			if ($isNT){
				$N1+=($ss1 =~ tr/[N]//);$N2+=($ss2 =~ tr/[N]//);
				while ($mask =~ /[^\0]/g) {
					my ($s1,$s2) = ( substr($ss1,$-[0],1),  substr($ss2,$-[0],1));#, ' ', $-[0], "\n";
					if ($s1 eq "-"){
						$N2++;next;
					}
					if ($s2 eq "-" ){#missing data, position doesn't matter
						$N1++;next;
					}
					$diffArs{$-[0]}=1;
					$diff++;
				}
			} else {
				while ($mask =~ /[^\0]/g) {
					my ($s1,$s2) = ( substr($ss1,$-[0],1),  substr($ss2,$-[0],1));#, ' ', $-[0], "\n";
					if ($s1 eq "N" ||$s1 eq "-"){
						$N2++;next;
					}
					if ($s2 eq "N" ||$s2 eq "-" ){#missing data, position doesn't matter
						$N1++;next;
					}
					$diffArs{$-[0]}=1;
					$diff++;
				}
			}
			my $nonDiff = ($mask =~ tr/[\0]//);
			$nonDiff -= $N1;
			$perID{$k1}{$k2}= $nonDiff/($diff+$nonDiff)*100;
		}
	}
	open O,">$opID" or die "Cant open out perc ID file $opID\n";
	my @smpls = keys %MS;
	print O "percID\t".join("\t",@smpls)."";
	foreach my $k1 (@smpls){
		print O "\n$k1\t";
		foreach my $k2 (@smpls){
			if ($k1 eq $k2){print O "\t100";
			} else {
				print O "\t".$perID{$k1}{$k2};
			}
		}
	}
	close O;
	my @NTdiffs = sort {$a <=> $b} (keys %diffArs);
	#die "Cant open out perc ID file $opID\n";
	my $MSAredF = $MSA; $MSAredF =~ s/\.[^\.]+$//;
	my $NSAposF = $MSAredF; $MSAredF.=".reduced.fna";
	$NSAposF .= ".reduced.pos";
	open O2,">$NSAposF" or die "Can't open reduced MSA position file $NSAposF\n";
	foreach my $i (@NTdiffs){
		print O2 "$i\n";
	} close O2;

	open O,">$MSAredF" or die "Can't open reduced MSA file $MSAredF\n";
	foreach my $k1 (@smpls){
		my $seq1 = $MS{$k1};
		my $seq = "";my $pos="";
		foreach my $i (@NTdiffs){
			$seq.=substr($seq1,$i,1);
			$pos .="$i\n";
		}
		print O ">$k1\n$seq\n";

	}
	close O;
	print "Dont calculating Distance matrix\n";
	#die "done\n";
}

sub calcDisPos($ $ $){
	my ($MSA,$opID, $isNT) = @_;
	#$cmd = $clustaloBin." -i $MSA -o $MSA.tmp --outfmt=fasta --percent-id --use-kimura --distmat-out $opID --threads=$ncore --force --full\n";
	#$cmd .= "rm -f $MSA.tmp\n";
	
	#too slow
	my $kr = readFasta($MSA);
	my %MS = %{$kr};
	my %diffArs;my %perID;
	print "Calculating distance matrix..\n";
	foreach my $k1 (keys %MS){
		$MS{$k1} = uc ($MS{$k1});
	}
	foreach my $k1 (keys %MS){
		my $ss1 = $MS{$k1};
		foreach my $k2 (keys %MS){
			next if ($k2 eq $k1);
			my $ss2 = $MS{$k2};
			my $mask = $ss1 ^ $ss2;
			my $diff=0;
			my$N2=($ss2 =~ tr/[-]//);
			my$N1=($ss1 =~ tr/[-]//);
			if ($isNT){
				$N1+=($ss1 =~ tr/[N]//);$N2+=($ss2 =~ tr/[N]//);
				while ($mask =~ /[^\0]/g) {
					my ($s1,$s2) = ( substr($ss1,$-[0],1),  substr($ss2,$-[0],1));#, ' ', $-[0], "\n";
					if ($s1 eq "-"){
						$N2++;next;
					}
					if ($s2 eq "-" ){#missing data, position doesn't matter
						$N1++;next;
					}
					$diffArs{$-[0]}=1;
					$diff++;
				}
			} else {
				while ($mask =~ /[^\0]/g) {
					my ($s1,$s2) = ( substr($ss1,$-[0],1),  substr($ss2,$-[0],1));#, ' ', $-[0], "\n";
					if ($s1 eq "X" ||$s1 eq "-"){
						$N2++;next;
					}
					if ($s2 eq "X" ||$s2 eq "-" ){#missing data, position doesn't matter
						$N1++;next;
					}
					$diffArs{$-[0]}=1;
					$diff++;
				}
			}
			my $nonDiff = ($mask =~ tr/[\0]//);
			$nonDiff -= $N1;
			$perID{$k1}{$k2}= $nonDiff/($diff+$nonDiff)*100;
		}
	}
	open O,">$opID" or die "Cant open out perc ID file $opID\n";
	my @smpls = keys %MS;
	print O "percID\t".join("\t",@smpls)."";
	foreach my $k1 (@smpls){
		print O "\n$k1\t";
		foreach my $k2 (@smpls){
			if ($k1 eq $k2){print O "\t100";
			} else {
				print O "\t".$perID{$k1}{$k2};
			}
		}
	}
	close O;
	my @NTdiffs = sort {$a <=> $b} (keys %diffArs);
	#die "Cant open out perc ID file $opID\n";
	my $MSAredF = $MSA; $MSAredF =~ s/\.[^\.]+$//;$MSAredF.=".reduced.fna";
	open O,">$MSAredF" or die "Can't open reduced MSA file $MSAredF\n";
	foreach my $k1 (@smpls){
		my $seq1 = $MS{$k1};
		my $seq = "";
		foreach my $i (@NTdiffs){
			$seq.=substr($seq1,$i,1);
		}
		print O ">$k1\n$seq\n";
	}
	close O;
	print "Dont calculating Distance matrix\n";
	#die "done\n";
}

sub mergeMSAs($ $ $ $){
	my ($MSAsAr,$samplesHr,$multAliF,$del,$isAA) = @_;
	my @MSAs = @{$MSAsAr}; my %samples = %{$samplesHr};
	my @smps = sort keys %samples;
	if (@smps == 0){#no cats file
		push(@MSAs ,$multAliF);
		return;
	}
	my %bigMSAFAAnxs;my %bigMSAFAA;foreach my $sm (@smps){$bigMSAFAA{$sm} ="";$bigMSAFAAnxs{$sm}="";}
	my @lengthsParts;
	foreach my $MSAf (@MSAs){
		#print $MSAf."\n"; 
		my $hit =0; my $miss =0;
		my $hr = readFasta($MSAf,1); my %MFAA = %{$hr};
		unlink $MSAf or die "Cannot remove merged MSA component $MSAf: $!\n" if $del && -e $MSAf;
		my @Mkeys = sort keys %MFAA;
		next if (@Mkeys == 0);
		my ($firstMsaSample, $gcat, $separator) = parseSeqId($Mkeys[0], "MSA header in $MSAf");
		my $len = length($MFAA{$Mkeys[0]});
		if ($len == 0){print STDERR "0 length sequence discovered in $MSAf\n";next ;}
		push(@lengthsParts,$len);
		foreach my $sm (@smps){
			my $curK = $sm.$separator.$gcat;
			if ( exists($MFAA{$curK}) && defined($MFAA{$curK})  ) {
				my $seq = $MFAA{$curK}; 
				die "Sequence lengths within MSA unequal: $len != ".length($seq)."\n" if (length($seq) != $len);
				$hit++;
				$bigMSAFAA{$sm} .= $seq;
				$seq =~ s/^(-+)/"?" x length($1)/e;
				$seq =~ s/(-+)$/"?" x length($1)/e;
				#die $seq;
				$bigMSAFAAnxs{$sm} .= $seq;
			} else {
				$bigMSAFAA{$sm} .= "-"x$len;
				$bigMSAFAAnxs{$sm} .= "?"x$len;
				$miss++;
			}
		}
		
		#die "$hit - $miss\n";
	}
	#filter part - count "-" in each seq
	my $factor = 1; $factor = 3 if ($isAA);
	my @ksMSAFAA = sort keys %bigMSAFAA;
	my $iniSeqNum = @ksMSAFAA; my $remSeqNum = 0;
	my %charCnts; my $maxNtCnt=0;
	#my @usedPos; #is now implemented in C++ program..
	#simply count gaps and N's
	foreach my $kk (@ksMSAFAA){
		#my $strCpy = $bigMSAFAA{$kk};
		my $num1 = 0;
		if ($isAA){
			$num1 = $bigMSAFAA{$kk} =~ tr/[\-Xx]//;
		} else {
			$num1 = $bigMSAFAA{$kk} =~ tr/[\-Nn]//;
		}
		
		$charCnts{$kk} = (length($bigMSAFAA{$kk})-$num1);
		#print "$kk GGGG  $charCnts{$kk} $num1\n";
		if ( $charCnts{$kk} > $maxNtCnt){
			$maxNtCnt = $charCnts{$kk};
		}
		next; #overlap implemented in C++
	}
	if ($maxNtCnt == 0){ #something really wrong
		die "No usable MSA positions remain after concatenation and filtering\n";
	}
	
	my $qtl90NTcnts = quantile(0.9,values(%charCnts));
	my $qtl50NTcnts = quantile(0.5,values(%charCnts));
	my $qtl25NTcnts = quantile(0.25,values(%charCnts));
	
	
	#final check on MSA's that enough data is present
	foreach my $kk (@ksMSAFAA){
		my $num1 = $charCnts{$kk};
		#print "$num1\n";
		if ( $maxNtCnt == 0 ||  ($num1 < ($qtl90NTcnts *$ntFrac) && $num1 < $qtl25NTcnts) || ($num1 < ($ntCntTotal/$factor) ) ){
			print "($num1 / $qtl90NTcnts ) < $ntFrac) && $num1 < $qtl25NTcnts || ($num1 < ($ntCntTotal/$factor) )\n";
			delete $bigMSAFAA{$kk}; delete $bigMSAFAAnxs{$kk}; $remSeqNum++; 
			if ($maxNtCnt == 0){
				print "$kk $num1  infinite $ntFrac \n";
			} else {
				print "$kk $num1  " . int($num1*1000 / ($maxNtCnt) )/1000 ." $ntFrac \n";
			}
		}
		#print "$num1  $kk \n";#$bigMSAFAA{$kk}\n\n"; last;
	}
	open O,">$multAliF" or die "Can't open MSA outfile $multAliF\n";
	open O2,">$multAliF.nxs" or die "Can't open MSA nexus outfile $multAliF.nxs\n";
	my @allKs = sort keys %bigMSAFAA;
	if (@allKs == 0){die "no genes for nexus output format.\nAborting\n";}
	print O2 "#NEXUS\nBegin data;\nDimensions ntax=".scalar(@allKs)." nchar=".length($bigMSAFAAnxs{$allKs[0]}).";\nFormat datatype=dna missing=? gap=-;\nMatrix\n";
	my $scnt=0;
	foreach my $kk (@allKs){
		print O ">$kk\n"; my $s1 = $bigMSAFAA{$kk};
		print O2 "\n$kk\t"; my $s2 = $bigMSAFAAnxs{$kk};
		#if (length($s1) !=  @usedPos){die "s1 and usedPos are not same length: ".length($s1)." : ".@usedPos."\n";}
		#if ($minOverlapMSA <=0){
		print O "$s1\n"; print O2 "$s2\n";
		#} else {
		#	my $usedP=0;
		#	for (my $i=0; $i<length($s1);$i++){
		#		if ($usedPos[$i] >= $minOverlapMSA){
		#			print O substr($s1,$i,1);
		#			print O2 substr($s2,$i,1);
		#			$usedP++;
		#		}
		#	}
		#	print O "\n"; print O2 "\n";
		#	print "Used $usedP positions of ". length($s1)  . " total positions.\n" if ($scnt==0);
		#	$scnt++;			
		#}
	}
	print O2 "\n;\nend;";

	close O;close O2;
	#die "$multAliF\n";
	
	#prepare partition file to record segment lengths..
	my $partiFile = $multAliF.$partiExt;
	open O,">$partiFile" or die "Can't open output partioning file $partiFile\n";
	my $lastP=0;
	my $TypeTag = "DNA";
	$TypeTag = "LG" if ($isAA); #this is the model to be used...
	for (my $i=0;$i<@lengthsParts;$i++){
		#DNA, part1 = 1-100
		print O "$TypeTag, part".($i+1) ." = ". ($lastP+1) ."-". ($lengthsParts[$i]+$lastP) ."\n";
		$lastP+=$lengthsParts[$i];
	}
	close O;
	
	print "Removed $remSeqNum of $iniSeqNum sequences\n";
}


sub convertMultAli2NT($ $ $){
	my ($inMSA,$NTs,$outMSA) = @_;
	my $tmpMSA=0;
	if ($inMSA eq $outMSA){$outMSA .= ".tmp"; $tmpMSA=1;}
	my $cmd = "";
	my $trimalBin = getProgPaths("trimal");
	#my $pal2nal = getProgPaths("pal2nal"); #"perl /g/bork3/home/hildebra/bin/pal2nal.v14/pal2nal.pl";

	#"$pal2nal $inMSA $NTs -output fasta -nostderr -codontable 11 > $outMSA\n";
	
	#$cmd = "$trimalBin -in $inMSA -out $outMSA -backtrans $NTs -keepheader -keepseqs -noallgaps -automated1 -ignorestopcodon\n";
	$cmd = "$trimalBin -in $inMSA -out $outMSA -backtrans $NTs -keepheader -ignorestopcodon  -gt 0.1 -cons 60 2>/dev/null\n";
	#die "$cmd\n$inMSA,$NTs,$outMSA\n";
	#my $hr1= readFasta($inMSA);
	#my %MSA = %{$hr1};
	#$hr1= readFasta($NTs);
	#my %NTs = %{$hr1};
	if ($tmpMSA){$cmd .= "rm -f $inMSA;mv $outMSA $inMSA;\n";}
	#print $cmd;
	#die "$cmd\n";
	systemW($cmd);
}

sub synPosOnlyAA($ $){#only leaves "constant" AA positions in MSA file.. 
#stupid, don't know if pal2nal can handle this.. prob not
	my ($inMSA,$outMSA) = @_;
	#print "Syn";
	my $hr = readFasta($inMSA,1); my %FNA = %{$hr};
	my @aSeq = keys %FNA;
	my $len = length ($FNA{$aSeq[0]});
	for (my $i=0; $i< $len; $i+=3){
		my $cod = substr $FNA{$aSeq[0]},$i,3;
		my $iniAA = "A";
		for (my $j=1;$j<@aSeq;$j++){
		}
	}
	#print " only\n";

}

sub synPosOnly{#now finished, version is cleaner
	my ($inMSA,$inAAMSA, $ffold, $outgroup, $doSyn, $doNSyn) = @_;
	
	my $outMSA = $inMSA; 	my $outMSAns = $inMSA;
	$outMSAns =~ s/\.fna/\.nonsyn\.fna/;
	$outMSA =~ s/\.fna/\.syn\.fna/;

	#print "Syn NT";
	my %convertor = (
    'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S',    # Serine
    'TTC' => 'F', 'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L', 'TTG' => 'L',    # Leucine
    'TAC' => 'Y',  'TAT' => 'Y',    # Tyrosine
    'TAA' => '*', 'TAG' => '*', 'TGA' => '*',    # Stop
    'TGC' => 'C', 'TGT' => 'C',    # Cysteine   
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L', 'CTT' => 'L',    # Leucine
    'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P',    # Proline
    'CAC' => 'H', 'CAT' => 'H',    # Histidine
    'CAA' => 'Q', 'CAG' => 'Q',    # Glutamine
    'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R',    # Arginine
    'ATA' => 'I', 'ATC' => 'I', 'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T',    # Threonine
    'AAC' => 'N','AAT' => 'N',    # Asparagine
    'AAA' => 'K', 'AAG' => 'K',    # Lysine
    'AGC' => 'S', 'AGT' => 'S',    # Serine
    'AGA' => 'R','AGG' => 'R',    # Arginine
    'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V',    # Valine
    'GCA' => 'A','GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A',    # Alanine
    'GAC' => 'D', 'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E', 'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G','GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G',    # Glycine
    );
	my %ffd;
	if ($ffold){ #calc 4fold deg codons in advance to real data
		foreach my $k (keys %convertor){
			my $subk = $k; my $iniAA = $convertor{$subk} ;
			my $cnt=0;
			foreach my $sNT ( ("A","T","G","C") ){
				
				substr ($subk,2,1) = $sNT;
				#print $subk ." " ;
				$cnt++ if ($convertor{$subk} eq $iniAA);
				
			}
#			if( $cnt ==4){ $ffd{$k} = 4;
#			} else {$ffd{$k} = 1;}
			if( $cnt ==4){ $ffd{$k} = 4;
			} else {$ffd{$k} = 1;}
			#die"\n$ffd{$iniAA}\n";
		}
	}

	#assumes correct 3 frame for all sequences in inMSA
	my $hr = readFasta($inMSA,1); my %FNA = %{$hr};
	#my %FAA;
	#if (0  || !$ffold){
	#	$hr = readFasta($inAAMSA); %FAA = %{$hr};
	#}
	#print "$inMSA\n$inAAMSA\n$outMSA\n";
	my @aSeq = sort keys %FNA;
	die "No sequences found in nucleotide MSA $inMSA\n" unless @aSeq;
	my %outFNA;#syn
	my %outFNAns;#non syn
	for (my $j=0;$j<@aSeq;$j++){$outFNA{$aSeq[$j]}="";}
	my $len = length ($FNA{$aSeq[0]});
	my $nsyn=0;my $syn=0;
	for (my $i=0; $i< $len; $i+=3){ #goes over every position
		my $j =0;
		my $iniAA = "-";
		my $iniCodon ;
		while (1){ #check for first informative position
			$iniCodon = substr $FNA{$aSeq[$j]},$i,3;
			if ($iniCodon =~ m/---/ || $iniCodon =~ m/N/i){$j++; last if ($j >= @aSeq); next;}
			die "error: $iniCodon\n" if ($iniCodon =~ m/-/); #should not happen
			die "codon doesn't exist $iniCodon \n" unless (exists($convertor{$iniCodon}));
			$iniAA = $convertor{$iniCodon};#substr $FAA{$aSeq[0]},$i,1; 
			last;
		}
		next if $j >= @aSeq && $iniAA eq "-";
		#die "$iniAA\n";
		my $isSame = 1;my $ntSame = 1;
		next unless (!$ffold || $ffd{$iniCodon} == 4);
	#print $i." $iniAA ";
		for (;$j<@aSeq;$j++){
			my ($seqSample) = parseSeqId($aSeq[$j], "synonymous-site MSA header", 1);
			if ($outgroup ne "" && $seqSample eq $outgroup){next;}
			my $newCodon = substr $FNA{$aSeq[$j]},$i,3;
			my $newAA = "-";
			if ($newCodon !~ m/-/ && $newCodon =~ m/[ACTG]{3}/i){
				die "Unkown AA $newCodon\n" unless (exists $convertor{$newCodon} );
				$newAA = $convertor{$newCodon} ; # substr $FAA{$aSeq[$j]},$i,1;
			} elsif ($newCodon =~ m/---/ || $newCodon =~ m/[NWYRSKMDVHB]/i){
			} else {
				die "newCodon wrong $newCodon\n" ;
			}
			if ($iniAA ne $newAA && $newAA ne "-"){
				$isSame =0; $ntSame =0; last;
			}
			if ($iniCodon ne $newCodon){
				$ntSame=0;
			}
		}
		$syn++ if !$ntSame && $isSame;
		$nsyn++ if !$isSame;
		for (my $j=0;$j<@aSeq;$j++){
			my $curCod = substr $FNA{$aSeq[$j]},$i,3;
			if ($ntSame){#add nts to file
				$outFNA{$aSeq[$j]} .= $curCod;
				$outFNAns{$aSeq[$j]} .= $curCod;
			} elsif ($isSame){#add nts to file
				if ($ffold){
					$outFNA{$aSeq[$j]} .= substr $FNA{$aSeq[$j]},($i)+2,1;
				} else {
					$outFNA{$aSeq[$j]} .= $curCod;
				}
				#print substr $FNA{$aSeq[$j]},$i*3,3 . " ";
			} else {
				$outFNAns{$aSeq[$j]} .= $curCod;
			}
		}
	}
	#die $inMSA."\n";
	if ($doSyn){
		if ($syn ==0){
			$outMSA = "";
		} else {
			open O ,">$outMSA" or die "Can't open outMSA $outMSA\n";
			for (my $j=0;$j<@aSeq;$j++){
				print O ">$aSeq[$j]\n$outFNA{$aSeq[$j]}\n";
			}
			close O;
		}
	}
	if ($doNSyn){
		if ($nsyn ==0){
			$outMSAns = "";
		} else {
			open O ,">$outMSAns" or die "Can't open outMSA $outMSAns\n";
			for (my $j=0;$j<@aSeq;$j++){
				print O ">$aSeq[$j]\n$outFNAns{$aSeq[$j]}\n";
			}
			close O;
		}
	}
	my ($reportedSample, $reportedGene) = parseSeqId($aSeq[0], "synonymous-site MSA header", 1);
	$reportedGene = "alignment" if $reportedGene eq "";
	#die "$outMSA\n";
	print "$reportedGene ($syn / $nsyn) ".@aSeq." seqs \n";
	#print " only\n";
	#print "\n";
	return ($outMSA,$outMSAns);
}

sub codeml{
	my ($MSAfile2,$codemlOutDTmp,$gene,$nwkFile_gene2,$repeatCounts) = @_;
	$pamlBin = getProgPaths("codeml") if $pamlBin eq "";
	my @omegaStart = @omegas;			
	my $codemlOutDFile = "$codemlOutD/${gene}_run2";
	system "mkdir -p  $codemlOutDFile" unless(-d $codemlOutDFile);

	chdir $codemlOutDTmp;
	my $modelName;
	for (my $mod=0; $mod < scalar(@model); $mod++){		
		$modelName = $model[$mod];
		my @repSel;
		for (my $rep=1; $rep <= $repeatCounts; $rep++) {

			open M0,">$codemlOutDTmp/${gene}_${rep}_${modelName}.c" or die "Can't open control file for codeml: $gene\n";
	
			print M0 "seqfile = $MSAfile2\n";
			print M0 "verbose = 2\n";
			print M0 "treefile = $nwkFile_gene2\n";
			print M0 "outfile = $codemlOutDTmp/codemlOut_${gene}_${rep}_${modelName}.txt\n";
			print M0 "aaDist = 0\n";
			print M0 "fix_blength = 0\n";
			print M0 "runmode = 0\n";
			print M0 "seqtype = 1\n";
			print M0 "CodonFreq = 2\n";
			print M0 "clock = 0\n";
			print M0 "model = 0\n";
			print M0 "NSsites = $modelName\n";
			print M0 "fix_omega = 0\n";
			print M0 "omega = $omegaStart[($rep-1)]\n";
			print M0 "cleandata = 0\n";
			print M0 "getSE = 0\n";	
			print M0 "icode = 0\n";
			print M0 "fix_kappa = 0\n";
			print M0 "kappa = 2\n";
			print M0 "Mgene = 0\n";
			print M0 "ncatG = 8\n";
			print M0 "RateAncestor = 0\n";
			print M0 "Small_Diff = 1e-6\n";
			print M0 "noisy = 0\n";
			

			close M0;

			## run codeml  
			$cmd = "$pamlBin $codemlOutDTmp/${gene}_${rep}_${modelName}.c\n";			
			systemW $cmd; 
			print "finished codeml Model $modelName run $rep on $gene\n";
			#die;
			
			#get lnL and push to array
			open(CM, "<$codemlOutDTmp/codemlOut_${gene}_${rep}_${modelName}.txt" ) or die "could not find $!";
			while (my $line = <CM>) {
					if ($line =~ /^lnL*/) {
							$line =~ m/^.*\):\s*([-+]?[0-9]*\.?[0-9]+)\s.*$/;
					push @repSel, $1;
					last;
					}				
			}
			close CM;
			#system "rm $codemlOutDTmp/rst1";

		} 
		
		my $idxMax = 0;
			$repSel[$idxMax] > $repSel[$_] or $idxMax = $_ for 1 .. $#repSel; 
		my $repSelected = $idxMax+1;
		#die "@repSel\n$repSelected\n";
		copy("$codemlOutDTmp/codemlOut_${gene}_${repSelected}_${modelName}.txt", "$codemlOutDFile/out_M$modelName.txt")
			or die "Cannot copy selected codeml result for $gene model $modelName: $!\n";
	}	
}


#starts my R script to calc popgenStats and filter out bad genes (in multi gene approach)
sub pogenStatsFilter{
	return [] unless ($doTheta);
	system "mkdir -p $MSAsubsD" unless (-d $MSAsubsD);
	system "mkdir -p $codemlOutD" unless (-d $codemlOutD);
	system "rm -f $codemlOutD/PopStats.*";
	my $condaAct = getProgPaths("CONDA");
	my $cmd = "";
	my $RpogenS = getProgPaths("pogenStats");

	#$cmd .= "$condaAct\nconda activate r_env\n"; #this was a workaround since the packages didn't run correctly on R cluster and different env was needed..
	$cmd .= "$RpogenS $outD $mapF $codemlOutD $MSAsubsD $subsetPopgenStats\n";
	#$cmd .= "conda deactivate\n";
	#die "$cmd\n";
	my $ret = `$cmd`;
	$ret =~ m/Outliers: (.*);;/;
	my @spl = split /,/,$ret;
	return \@spl;
} 

sub hyphy{
	my ($MSAfile2,$codemlOutDTmp,$gene,$nwkFile_gene2,$log) = @_;
	my $hyphyBin=getProgPaths("hyphy");
	my $cmd = "";#"source activate hyphy\n";
	$cmd .= "$hyphyBin CPU=$ncore fubar --alignment $MSAfile2 --tree $nwkFile_gene2 > $log\n";
	$cmd .= "gzip -c $MSAfile2.FUBAR.json > $log.json.gz\n";
	$cmd .= "rm -f $MSAfile2.FUBAR.*\n";
	print $log."\n";
	systemW $cmd ;#if (!-e $log);
	#die $cmd ;#if (!-e $log);
	return $log;
}

sub pruneTree($ $ $){
	my ($nwkFile,$aR,$nwkFile_gene) = @_;
	my @genomeList = @{$aR};
	my $eteBin = getProgPaths("ete3");
	my $cmd_prune = "$eteBin mod -t $nwkFile --prune @genomeList --unroot -o $nwkFile_gene";
	systemW $cmd_prune . " > $nwkFile_gene";
}

sub fubarXML($){
	my $inp = $_[0];
	return "" if (!-e $inp || (-e "$inp.json.gz" && !$reparseHyphyJson));
	
	if (0){#python
		return "";
	}
	my $rDNm; my $rDSm;
	my $txt = `zcat $inp.json.gz`;
	my $str=decode_json($txt);
	my %JS = %{$str};
	my @XX = @{$JS{MLE}{"content"}{"0"}};
	#print @XX."  BB  @XX\n";
	my $negSel=0; my $posSel=0; 
	my @ds = (); my @dn = ();
	for (my $i=0;$i<@XX;$i++){
		my @YY = @{$XX[$i]};
		push @ds,$YY[0];
		push @dn,$YY[1];
		next if (@YY<4 || !defined($YY[4]));
		$negSel++ if ($YY[3]>0.9);
		$posSel++ if ($YY[4]>0.9);
	}
	my $dnX=$2;my $dsX = $1;
	my $rDNmed = medianArray(@dn);my $rDSmed = medianArray(@ds);
	$rDNm = meanArray(\@dn); $rDSm = meanArray(\@ds);
#	$txt = `cat $inp`;
#	$txt =~ m/\* synonymous rate =  ([\d\.]+)\n.*\* non-synonymous rate =  ([\d\.]+)/;
	#print "hy report dnds: $dnX $dsX; $rDNm  $rDSm; $rDNmed $rDSmed\n"; #this is actually probably wrong
	return $JS{input}{"number of sequences"} ."\t". $JS{input}{"number of sites"} . "\t$negSel\t$posSel\t$rDNm\t$rDSm";
}

sub coreHyPhy{
	my ($MSADir,$gene,$xtra,$nwkFile,$codemlOutDTmp,$logF) = @_;
	return if (-e $logF && -e "$logF.json.gz");
	my $runCodeML = 0;
	my @genomeList;
	opendir(DIR, $MSADir);
	my @MSAfile = grep(/\A\Q$gene\E\.\d+.*$xtra\.fna\z/,readdir(DIR));
	closedir(DIR);
	if (@MSAfile == 0){
		#die "$gene.*$xtra\.fna";
		return;
	}
	#die "@MSAfile\n";
	
	my $MSAfile2 = "$MSADir/$MSAfile[0]";
	my $MSAfile3 = "$codemlOutDTmp/tmp.$MSAfile[0]";
	#get entries in nwk to match these up as well..
	my $hr = getTreeLeafs($nwkFile);
	my %nwLfs= %{$hr};
	
	#die "@nwLfs\n";
	my $cntMissTree=0;
	$hr = readFasta($MSAfile2,1,"input MSA for selection analysis");
	my %FNA = %{$hr};
	#		print "$MSAfile2\n";

	open O,">$MSAfile3" or die "can't open MSA out $MSAfile3\n";
	foreach my $genome (keys %FNA){
		#print "$genome\n";
		my ($genome2) = parseSeqId($genome, "selection-analysis MSA header");
		next if ($genome2 =~ m/$outgroup/);
		unless (exists($nwLfs{$genome2})){$cntMissTree++; next;}
		my $seq = $FNA{$genome};
		$seq =~ s/T[AG][GA]$//; #remove stop codon
		my $x=0;
		foreach my $sto (("TGA","TAG","TAA")){
			while ( 1 ){
				$x = index($seq,$sto,$x);
				last if ($x < 0);
				if ($x % 3 ==0){
					print "HIT\n";
					substr($seq,$x,3) = "NNN";
				}
				$x++;
			}
		}
		print O ">$genome2\n$seq\n";
		push (@genomeList, $genome2);
	} 
	close O;
	if ($cntMissTree>0){print "Removed from MSA $cntMissTree sequences due to not being present in tree\n";}
	next if(scalar @genomeList <3);
			
	my $nwkFile_gene = "$codemlOutDTmp/subtree_$gene.nwk";

	####################### Prepare tree ################################## 
	pruneTree($nwkFile,\@genomeList,$nwkFile_gene);
	#$thrs[$thrCnt]->join();
	if ($runCodeML){
		codeml($MSAfile3,$codemlOutDTmp,$gene,$nwkFile_gene,$repeatCounts) #$thrs[$thrCnt] = threads->create( );
	} else {
		hyphy($MSAfile3,$codemlOutDTmp,$gene,$nwkFile_gene,$logF);
	}
	system "rm -f $MSAfile3";
	system "rm -f $nwkFile_gene";
}

sub selecAnalysis($ $ $ $ $){
	my ($geneListRef, $nwkFile, $codemlOutD, $codemlOutDTmp) = @_;
	my @geneListFin = @{$geneListRef};
	#system "source activate hyphy" unless ($runCodeML);
	my $jsonExrScr = getProgPaths("fubarJson_scr");
	my $stdJSONheader = "Gene\tNseqs\tNsites\tnegSel\tposSel\tdn\tds\tdn2\tds2\n";
	#my @thrs;my $thrCnt=0;
	#for ($thrCnt=0;$thrCnt<$ncore;$thrCnt++){		$thrs[$thrCnt] = threads->create(sub{my $x=0;});	}
	#$thrCnt=0;
	#first go through $MsaD
	system "rm -f $codemlOutD/hyphy.fubar*" if ($reparseHyphyJson);
	my $logF1 =  "$codemlOutD/hyphy.fubar.txt";
	if (!-e $logF1 ){
		my %logs;
		foreach my $gene (@geneListFin){
			my $gene_file_stem = geneFileStem($gene);
			my $logF =  "$codemlOutD/$gene_file_stem.hyphy.fubar.log";
			$logs{$gene} = "$logF";
			coreHyPhy($MsaD,$gene_file_stem,"",$nwkFile,$codemlOutDTmp,$logF);
		}
		my $sumTxt=$stdJSONheader;
		print "reading jsons..";
		foreach my $gene (keys %logs){
			#my $summary = fubarXML("$logs{$gene}");
			next unless (-e "$logs{$gene}.json.gz");
			my $summary = `$jsonExrScr $logs{$gene}.json.gz`;
			next if ($summary eq "");
			$sumTxt .= "$gene\t$summary";
			#print "$sumTxt\n";
		}
		print "Done reading\n";
		open O,">$logF1";		print O $sumTxt;		close O;
	}
	#add by hand the unique seqs subset..
	$logF1 =  "$codemlOutD/hyphy.fubar.unID.txt";
	if (!-e $logF1 ){
		my %logs;
		foreach my $gene (@geneListFin){
			my $gene_file_stem = geneFileStem($gene);
			my $logF =  "$codemlOutD/$gene_file_stem.hyphy.fubar.unID.log";
			$logs{$gene} = $logF;
			coreHyPhy($MSAsubsD,$gene_file_stem,"\\.uInd",$nwkFile,$codemlOutDTmp,$logF);
		}
		my $sumTxt=$stdJSONheader;
		print "reading jsons unID..";
		foreach my $gene (keys %logs){
			#my $summary = fubarXML("$logs{$gene}");
			next unless (-e "$logs{$gene}.json.gz");
			my $summary = `$jsonExrScr $logs{$gene}.json.gz`;
			next if ($summary eq "");

			$sumTxt .= "$gene\t$summary";
			#system "rm -f $logs{$gene}*";
		}
		print "Done reading\n";
		open O,">$logF1";		print O $sumTxt;		close O;
	}
	#now go through subsets..
	
	foreach my $subs (split /,/,$subsetPopgenStats){
		my %logs;
		#print "$subs\n";
		my $logF1 =  "$codemlOutD/hyphy.fubar.s$subs.txt";
		next if (-e $logF1 );
		foreach my $gene (@geneListFin){
			my $gene_file_stem = geneFileStem($gene);
			my $logF =  "$codemlOutD/$gene_file_stem.hyphy.s$subs.fubar.log";
			$logs{$gene} = $logF;
			#COG0008.0.uInd.s20.fna
			coreHyPhy($MSAsubsD,$gene_file_stem,"uInd\\.s$subs",$nwkFile,$codemlOutDTmp,$logF);
		}
		my $sumTxt=$stdJSONheader;
		print "reading jsons s$subs..";
		foreach my $gene (keys %logs){
			#my $summary = fubarXML("$logs{$gene}");next if ($summary eq "");
			next unless (-e "$logs{$gene}.json.gz");
			my $summary = `$jsonExrScr $logs{$gene}.json.gz`;
			next if ($summary eq "");
			$sumTxt .= "$gene\t$summary";
			#system "rm -f $logs{$gene}*";
		}
		print "Done reading\n";
		open O,">$logF1";		print O $sumTxt;		close O;
	}
	system "rm -fr $codemlOutDTmp" if (-e $codemlOutDTmp);
	
}


#\@geneList,$MsaD,$codemlOutD
sub WattTheta{
	#hyphy /path/to/WattetrsonTheta.txt --alignment /path/to/alignment
	my ($geneListRef, $MSADir, $codemlOutD) = @_;
	my @geneListFin = @{$geneListRef};
	my $runCodeML = 0;
	#system "source activate hyphy" unless ($runCodeML);
	my $hyphyBin=getProgPaths("hyphy");
	#die "XX\n";
	my %logs ;
	my $logF1 =  "$codemlOutD/hyphy.Theta.log";
	return if (-e $logF1);
	my $otxt = "Gene\tNseqs\tNsites\tSegSites\tWattTheta\n";
	foreach my $gene (@geneListFin){
		my $cnt = 0;
		my $gene_file_stem = geneFileStem($gene);
		my $logF =  "$codemlOutD/$gene_file_stem.hyphy.Theta.log";
#		print "$logF\n";
		$logs{$gene} = $logF;
		#next if (-e $logF);
		my @genomeList;
				
		opendir(DIR, $MSADir);
		my @MSAfile = grep(/\A\Q$gene_file_stem\E\.\d+.*\.fna\z/,readdir(DIR));
		closedir(DIR);
		#print "@MSAfile\t$gene\t$MSADir\n";
		next if (@MSAfile == 0);
		my $MSAfile2 = "$MSADir/$MSAfile[0]";
		my $hyphyBin=getProgPaths("hyphy");
		my $cmd = "";#"source activate hyphy\n";
		$cmd .= "$hyphyBin CPU=$ncore /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/secScripts/phylo/WattetrsonTheta.hyphy --alignment $MSAfile2 ";#> $logF\n";
		my $txt = `$cmd`;
#		print $txt."\n";
		$txt =~ m/Sequences          = (\d+)\nSites              = (\d+)\nSegregating Sites  = (\d+)\n.*Watterson.s theta  = ([\d\.]+)/;
		$otxt.="$gene\t$1\t$2\t$3\t$4\n";
	}
	#foreach my $g (keys %logs){	my $txt = `cat $logs{$g}`;	}
	open O,">$logF1" or die "Can't open log out $logF1\n";
	print O $otxt;
	close O;
}
sub fillGeneList{
	my ($cogCats) = @_;
#	open my $xI,"<$cogCats" or die "Can't open cogcats $cogCats\n";
	my ($xI,$ST)= gzipopen($cogCats,"CogCATs phylo");
	chomp(my @linesCats = <$xI>);
	close $xI;
	foreach (@linesCats){
		chomp; my @spl = split /\t/;
		@spl = grep !/^NA$/, @spl;#remove NAs
		if ($spl[0] =~ m/^#/){shift @spl;}
		my @spl2 = parseSeqId($spl[0], "gene-list category");
		push(@geneList, $spl2[1]);				
	}
}


sub parseSeqId{
	my ($seqId, $context, $allowUndelimited) = @_;
	$context ||= "sequence identifier";
	if (defined($seqId)
		&& $seqId =~ /^(?<sample>.*?)(?<separator>$smplSep)(?<gene>.+)$/
		&& $+{sample} ne "" && $+{gene} ne ""){
		return ($+{sample}, $+{gene}, $+{separator});
	}
	return ($seqId, "", "") if $allowUndelimited && defined($seqId) && $seqId ne "";
	die "Cannot split $context '$seqId' with -smplSep '$smplSep'\n";
}


sub geneFileStem{
	my ($gene) = @_;
	die "Cannot create a filename for an empty gene/locus identifier\n"
		unless defined($gene) && length($gene);
	my $stem = $gene;
	$stem =~ s/_/__/g;
	$stem =~ s/([^A-Za-z0-9.-])/sprintf("_%02X", ord($1))/ge;
	return $stem;
}


sub shellQuote{
	my ($value) = @_;
	$value = "" unless defined $value;
	$value =~ s/'/'"'"'/g;
	return "'$value'";
}


sub safeRemoveTree{
	my ($path, $parent) = @_;
	return unless defined($path) && ($path ne "") && (-d $path || -l $path);
	my $absolutePath = File::Spec->canonpath(File::Spec->rel2abs($path));
	my $absoluteParent = File::Spec->canonpath(File::Spec->rel2abs($parent));
	my $relative = File::Spec->abs2rel($absolutePath, $absoluteParent);
	die "Refusing to remove $absolutePath outside $absoluteParent\n"
		if $relative eq File::Spec->curdir || $relative =~ /^\.\.(?:[\\\/]|$)/;
	my $errors;
	remove_tree($absolutePath, {error => \$errors});
	if ($errors && @{$errors}){
		my @messages;
		for my $entry (@{$errors}){
			my ($failedPath, $message) = %{$entry};
			push @messages, "$failedPath: $message";
		}
		die "Failed to remove directory tree $absolutePath: ".join("; ", @messages)."\n";
	}
}


sub requireConfiguredTool{
	my ($environmentName, $description) = @_;
	my $path = $ENV{$environmentName} // "";
	die "$description support is dormant and not configured. Set $environmentName to reactivate it.\n"
		if $path eq "";
	die "$description configured by $environmentName does not exist: $path\n" unless -e $path;
	return shellQuote(File::Spec->canonpath(File::Spec->rel2abs($path)));
}



### Fastgear -> test for recombination 
sub runFastgear($ $ $ $){
	my ($geneFG, $outFile, $inD, $parFile) = @_;
	my $fastgearBin = requireConfiguredTool("MF4_FASTGEAR_BIN", "fastGEAR");
	my $matlabBin = requireConfiguredTool("MF4_FASTGEAR_MATLAB_BIN", "fastGEAR MATLAB runtime");

	$cmd = "$fastgearBin $matlabBin ".shellQuote("$inD/$geneFG.fna")
		." ".shellQuote($outFile)." $parFile";
	systemW($cmd);
	die "fastGEAR did not produce $outFile for $geneFG\n" unless -s $outFile;
	print "fastgear on $geneFG finished";
	#die;
}




sub prepGenoDirs($){
	my ($genoindir) = @_;
	# -genoInD '/hpc-home/hildebra/grp/data/DB/Genomes/Ecoli_ref/*.fna'
	if ($genoindir eq ""){return;}
	$aaFna="$outD/autoFAA.faa";$fnFna="$outD/autoFNA.fna";$cogCats="$outD/autoCAT.cat";
	my $SaSe = "_"; my %catT;
	$smplSep = "_";
	
	#return(); #DEBUG
	if (-d $genoindir){
		if ($wildcardflag ne ""){
			$genoindir .= $wildcardflag;#"/*\.fna";
		} else {
			print "Warning: please give wildcard for genoIndir.\nAssumming *.fna now..\n";
			$genoindir .= "/*\.fna";
		}
	}
	#print "$genoindir\n\n";
	my @sfiles = glob($genoindir);
	if (scalar(@sfiles) == 0){die "$genoindir contains no files!\n";}
	
	#die "@sfiles\n";
	print "External Genomes: $genoindir\nProcessing ".scalar(@sfiles)." genomes\n";
	my %FNAfmg; my %FAAfmg;my %MGSFMG;
	
	my $cnt=0;
	foreach my $tarG(@sfiles){
		next if ($tarG  =~ /\*/ || $tarG  =~ /\.genes\./);
	#my $tarG = "$tarDir/$genoN";
		my $tag = $tarG;
		$tag =~ s/.*\///;
		$tag =~ s/\.[^\.]*$//;
		print "$tarG\n";
		my ($genes,$prots) = getGenoGenes($tarG,0,$ncore);
		my $FMGdir = getFMG("",$prots,$genes,$ncore);
		my ($hrN,$hrA,$hrC) = readFMGdir( "$FMGdir",$tag ,".");
		$cnt++;
		my %FAA=%{$hrA};
		my %FNA=%{$hrN};
		my %COGcat=%{$hrC};
		#add to categories..
		foreach my $k (keys %FAA){
			$FAAfmg{$k}=$FAA{$k};
			$FNAfmg{$k}=$FNA{$k};
			die "unkown key $k \n" unless (exists($COGcat{$k}));
			$MGSFMG{$tag}{$COGcat{$k}} = $k;
		}
		#print " $tag ";
		#die;
		#last if ($cnt>5);#DEBUG
	}
	
	
	open OA,">$aaFna"  or die "Can't open faa out file $aaFna\n"; 
	open ON,">$fnFna"  or die "Can't open faa out file $fnFna\n"; 
	foreach my $mg (keys %MGSFMG){
		foreach my $cog (keys %{$MGSFMG{$mg}}){
			my $ng = "$mg$SaSe$cog";
	#		print ON ">$ng\n$FNAfmg{$MGSFMG{$mg}{$cog}}\n";
			die "$MGSFMG{$mg}{$cog}\n" unless (exists( $FAAfmg{$MGSFMG{$mg}{$cog}} ));
			print OA ">$ng\n$FAAfmg{$MGSFMG{$mg}{$cog}}\n";
			print ON ">$ng\n$FNAfmg{$MGSFMG{$mg}{$cog}}\n";
			push(@{$catT{$cog}},"$ng");
		}
	}
	close OA; close ON;

	open OC,">$cogCats" or die "Can't open cat file $cogCats\n";
	foreach my $cg (keys %catT){
		print OC join("\t",@{$catT{$cg}})."\n";
	}
	close OC;
	
	
	#die "$cogCats\n";
	
}
