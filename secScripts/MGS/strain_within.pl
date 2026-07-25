#!/usr/bin/perl
#This script collects representative (consensus SNP called) DNA sequences from different metagenomic samples for each MGS, and submits buildTree5.pl for each to reconstruct a phylogeny
# check performance: /ei/projects/8/88e80936-2a5d-4f4a-afab-6f74b374c765/data/geneCats/famDrama7/Bin_SB/intra4_28Feb_01D2SV/MGS.10

use warnings;
use strict;

use Getopt::Long qw( GetOptions );
use List::Util qw/shuffle/;
use File::Path qw(make_path remove_tree);
use File::Glob qw(bsd_glob);
use File::Copy qw(copy);
use File::Basename qw(basename dirname);
use File::Spec;
use Cwd qw(abs_path getcwd);



use Mods::GenoMetaAss qw(gzipopen fileGZe fileGZs readClstrRev systemW median mean readMapS readFasta getAssemblPath getAssemblGFF getAssemblContigs checkSeqTech);
use Mods::Subm qw(qsubSystem emptyQsubOpt qsubSystem2 qsubSystemJobAlive qsubSystemWaitMaxJobs);
use Mods::IO_Tamoc_progs qw(getProgPaths truePath);
use Mods::TamocFunc qw(readTabbed getFileStr checkMF);
use Mods::geneCat qw(readGene2tax createGene2MGS);
use Mods::math qw(quantileArray);
use Mods::MGSLocus qw(build_locus_groups choose_locus_candidate protein_kmer_similarity robust_depth_mask);
use Mods::StrainParts qw(
	exact_worker_parts write_split_generation write_worker_completion
	split_generation_complete clear_split_generation resolve_fasta_artifact
	append_fasta_records_atomic
);

sub extractFNAFAA2genes;
sub histoMGS;
sub readGenesSample_Singl;
sub reportingsMGS;
sub prepRun;
sub prepGene2MGS;
sub createAGlist; sub preComputeConsSNP;
sub mergeConspecificLogs;
sub timeNice;
#sub combineMGSgenes;
sub combineMGSgenesDir; sub splitWorkerPartsRemain; sub getInputSize;
sub evalFileStatus;
sub addOutgroup2MGS;
sub writeTooFewMarker;
sub treeInputPrecopyCommand;
sub readFastaIDs;

sub limitedWarn;sub limitedNotice;


my %limitedWarningStats;
my %limitedNoticeStats;
my $warningExampleLimit = 5;
sub limitedWarn {
	my ($category, $message) = @_;
	my $entry = $limitedWarningStats{$category} ||= { total => 0, suppressed => 0 };
	$entry->{total}++;
	if ($entry->{total} <= $warningExampleLimit) {
		warn $message;
	} else {
		$entry->{suppressed}++;
		warn "Further '$category' warnings are suppressed; a total will be reported at exit.\n"
			if $entry->{total} == $warningExampleLimit + 1;
	}
}

sub limitedNotice {
	my ($category, $message) = @_;
	my $entry = $limitedNoticeStats{$category} ||= { total => 0, suppressed => 0 };
	$entry->{total}++;
	if ($entry->{total} <= $warningExampleLimit) {
		print $message;
	} else {
		$entry->{suppressed}++;
		print "Further '$category' messages are suppressed; a total will be reported at exit.\n"
			if $entry->{total} == $warningExampleLimit + 1;
	}
}

END {
	my @suppressed = sort grep {
		($limitedWarningStats{$_}{suppressed} || 0) > 0
	} keys %limitedWarningStats;
	if (@suppressed) {
		warn "\nSuppressed warning summary:\n";
		for my $category (@suppressed) {
			my $entry = $limitedWarningStats{$category};
			warn "  $category: $entry->{total} total; $entry->{suppressed} not shown\n";
		}
	}
	my @noticeSuppressed = sort grep {
		($limitedNoticeStats{$_}{suppressed} || 0) > 0
	} keys %limitedNoticeStats;
	if (@noticeSuppressed) {
		print "\nRepeated status summary:\n";
		for my $category (@noticeSuppressed) {
			my $entry = $limitedNoticeStats{$category};
			print "  $category: $entry->{total} total; $entry->{suppressed} not shown\n";
		}
	}
}


#v.14: reworked massively how many genes get included
#v.15: included lessons learned from MGS.pl v0.21 
#v.16 added familyVar and groupStabilityVars arguments for stability calculations
#v.17: considerations to improve speed of intial fna/faa extractions..
#v.18: 16.11.24: handling genes occurring >1 in a single sample/assembly
#v.19: 17.11.24: stricter filtering of genes, removing entire MGS if too many "bad genes" in them; added abundance based filtering of genes/sample
#v0.20: 22.11.24: fixed bug with v0.19 no longer accepting assmblGrps. code refactor that makes it a lot easier to understand
#v0.21: 2.1.25: v0.20 fix, to only select single gene instead of COG; further changed how genes are selected, to reomve potentially conspecific MGS per sample (instead of removing entire gene)
#v0.22: added per sample (not assmblGrp) MGS filtering based on multigenes
#v.23: removed MGS conspecific filter: was too harsh and didn't make sense to have a global filter: MGS are conspecific in a single sample, not all samples..
#.24: 31.10.25: on-the-fly creation of SNP consensus fastas, if correct vcf present
#.25: 22.12.25: precompute for vcf2fna added
#.26: 28.12.25: code refactor to later enable parallelization of main gene-collecting routine
#.27: 12.2.26: made code faster and more stable. changed default MSA aligner
#.28: 22.2.26: claude suggested code improvements
#.29: 24.2.26: switched to multi output file for subjobs (waits were to long/inconsistent performance and errors with file blocks)
#.30: 25.2.26: new code for combining files, subfiles written to scratch to improve speed further
#.31: 26.2.26: better integration new temp files, pick up from previous job, sorting jobs
#.32: 27.2.26: allows for subsets of MGS only to be calculated.. (good for testing)
#.33: 7.3.26: speed improvements across the board, more options for vcf2dna
#.34: 28.4.26: custom bin file
#.35: validate inputs and repair resume, outgroup, consensus, and temporary-directory handling
#.36: preserve locus-level same-COG genes, resolve paralogs, and make sample filters robust to sparse inputs
#.37: expose tree IDs as sample|COG|primaryGeneID while retaining MGS-qualified internal locus keys
#.38: validate paired consensus inputs, split-job logs, scheduler state, and destructive paths
#.39: make split retries generation-safe, merges atomic, and compressed outgroup updates reliable
#.40: bound repetitive data warnings, summarize suppressed diagnostics, and clarify progress output
#.41: make generated tree-input publication safe to rerun after scratch files have already moved
#.42: resubmit unfinished trees from published inputs without requiring scratch aggregates
#.43: avoid redundant candidate scoring and hot-loop container copies during extraction
#.44: reduce locus-model, FASTA scan, and category-publication peak memory
#.45: normalize repeated VCF headers and distinguish split-worker sparsity from missing catalogue data
my $version = 0.45;


my $cmdCall = join(" ", $0, @ARGV) . "\n";

my $pigzBin  = getProgPaths("pigz");

#input args..
my $GCd = "";#$ARGV[0];
my $MGSfile = "";#$ARGV[1];
my $clusterID = 95;
my $geneSelFile = "";
my $numCores = 4;#$ARGV[2];
my $subJob=0;#if 0, is main submitting job..
my $maxSubJob = 0;#into how many subjobs to split??
my $outDpre = "";
my $locTmpDir = ""; my $locTmpDir1 = "";
my $maxCores = -1;
my $onlySubmit =0;#extract genes anew?
my $reSubmit=0;
my $treeFile = "";
my $doSubmit=0;
my $subMode="";
my $multiGeneSmplMax = 0.25; #no higher than this rate in single samples conspec genes..
my $conspGeneSmplMax = 0.05; #no higher than this conspecific genes/MGS/sample
my $minDepthGene  = 1;
my $minBadLociForSampleSkip = 3; #avoid rejecting sparse MGS for a single bad locus
my $noIndels = 1;



my $repairCAT=0;

my $maxNGenes = 400;
my @subsetMGS=(); my $subsMGSstr="";
my $MSAprog = 2; ##(0) MSAprobs, (1) clustalO, (2) mafft, (4) MUSCLE5
my $phyloProg = 1; # #1=iqtree-fast, 2=veryfasttree, 3=fasttree
my $GenesPerSpecies = 0.2; #was previously 0.1.. maybe too low?
my $GeneLengthMin = 0.5;
my $presortGenes = 1200;
my $checkMaxNumJobs = 400;
my $useGTDBmg = "GTDB";
my $selfMemGb = 10;
my $redoSubmissionData = 0;
my $deepRepair = 0;
my $rmMSA = 1; #argument passed to buildTree5.pl 
my $contTests = ""; my $discTests = ""; #stat tests to be given to strain_within_2.2.pl
my $familyVar = ""; my $groupStabilityVars = "";

#SNP calling
my $minSNPDepth = 2; #changed to two: seems to give better results
my $minSNPCallQual = 20; #this is weak evidence in metag context
my $useAdaptiveQual = 0.0; #adaptive quality filtering in vcf2fna (based on depth)? Default: 0 (deactivated)
my $depthFilterScale =0.15; # if DP < mean contig depth *x, filter. Default: 0.15
my $indelRange = 5; #SNPs in range of X bp indels will be excluded
my $forceVCF2FNA = 0; #force the recalc of cons fasta from vcf..
my $SNPconsLOGs = ""; #logs for recalculating cons SNPs
my $preCompCons=0; #if >0, precompute in these blocks

my $takeAll = 0;
my $conspecificSpThr = 0.1; #higher fraction of genes being two copies in the same sample (abundance >0), and the whole MGS is removed from that sample
my $MGStoolowGsThr = 10; #less genes than this in a single sample -> rm MGS from sample for strains
my $mode = "MGS";
my $appendWriteTrigger = 200; #every Xth samples, genes are written (to manage memory); #limit this, perl seems to have some issues with too large strings..
my $startSubFromMGS = ""; #debug option: only start resubmitting tree building from this MGS (e.g. "MGS.1382" )
#define local files..
my $lSNPdir="SNP"; my $lMAPdir = "mapping";
my $lConsFNA = "genes.shrtHD.SNPc.MPI.fna.gz";
my $lConsCTG = "contig.SNPc.MPI.fna.gz";
my $lConsFAA = "proteins.shrtHD.SNPc.MPI.faa.gz";
my $SNPcaller = "MPI";
my $lConsVCF = "allSNP.${SNPcaller}.vcf.gz";
my $lConsVCFsup = "allSNP.${SNPcaller}-sup.vcf.gz";


#set up some base paths specific to pipeline..
my $FNAstdof = "allFNAs.fna"; my $FAAstdof = "allFAAs.faa";
my $LINKstdof = "link2GC.txt"; my $CATstdof = "all.cat";
my $abundF="/assemblies/metag/ContigStats/Coverage.pergene.gz";
my $bamDepthFsuffix = "-smd.bam.coverage.gz";
my $bamDepthFsuffixSup = ".sup-smd.bam.coverage.gz";
my $mapF2 = "";
my $memMulti = 1; #for buildTree script


checkMF();
#$treeFile = $ARGV[3] if (@ARGV > 3);$onlySubmit = $ARGV[4] if (@ARGV > 4);
#$doSubmit = $ARGV[5] if (@ARGV > 5);$subMode = $ARGV[6] if (@ARGV > 6);


GetOptions(
	"GCd=s"          => \$GCd,
	"clusterID=i"    => \$clusterID,
	"outD=s"         => \$outDpre,
	"MGS=s"          => \$MGSfile,
	#"geneSel=s"      => \$geneSelFile,
	"map2=s"         => \$mapF2, #to be given to strain2 script
	"nodeTmp|tmpD=s" => \$locTmpDir1, 
	"submit=i"       => \$doSubmit,
	"selfMemGb=i"    => \$selfMemGb,
	"onlySubmit=i"   => \$onlySubmit, #submit only jobs, or also recreate input fna/faa files? (can take days)
	"reSubmit=i"     => \$reSubmit, #for all MGS: resubmit tree phylo building
	"repairCAT=i"    => \$repairCAT,
	"deepRepair=i"   => \$deepRepair, #for missing MGS phylos: will resubmit phylo and rebuild fna/faa 
	"redoSubmissionData=i" => \$redoSubmissionData,  #for all MGS: will resubmit phylo and rebuild the fna/faa files..
	#workflow HPC usage
	"subjob=i"       => \$subJob,
	"maxSubJob=i"    => \$maxSubJob,
	"treeSubFromMGS=s" => \$startSubFromMGS, #debug option..
	#"cores=i"        => \$numCores, #not used any longer..
	"maxCores=i"     => \$maxCores, #superseedes -cores, will dynamically allocate num cores based on input file size, if defined
	"presortGenes=i" => \$presortGenes, #how many potential genes to include, of the original MGS (receovered will vary strongly  between samples)
	"maxGenes=i"     => \$maxNGenes, #how many genes to try to include? -> will be decided on each samples
	"flushEvery=i"   => \$appendWriteTrigger, #samples buffered before per-MGS records are flushed
	
	"forceSNPcalls=i"  => \$forceVCF2FNA,
	"preCompConsSNP=i"   => \$preCompCons,
	"MGSsubset=s"    => \$subsMGSstr,
	"submissionMode=s"      => \$subMode,
	"MGset=s"        => \$useGTDBmg,
	
	#used genes fine tuning..
	"MGSminGenesPSmpl=i" => \$MGStoolowGsThr, #less genes than this in a single sample -> rm MGS from sample for strains. default 10
	"multiGeneSmplMax=f" => \$multiGeneSmplMax, #default 0.15
	"conspGeneSmplMax=f" => \$conspGeneSmplMax, #default 0.05
	"minBadLociPSmpl=i" => \$minBadLociForSampleSkip,
	
	#transferred to buildTRee script..
	"GenesPerSpecies=f" => \$GenesPerSpecies,
	"GeneLengthMin=f" => \$GeneLengthMin,
	"MSAprog=i"      => \$MSAprog, #2=MAFFT, 4=muscle5
	"phyloProg=i"    => \$phyloProg, #1=iqtree-fast, 2=veryfasttree
	"rmMSA=i"        => \$rmMSA, #remove MSA, to save diskspace
	"phyloMemMulti=f" => \$memMulti, #mem used for buildtree. Default: 1.0
	
	"MGSphylo=s"     => \$treeFile,
	#transferred to MG-STK
	"ContTests=s"      => \$contTests, #continous stat tests to be handed to next step (just a passthrough)
	"DiscTests=s"      => \$discTests, #discrete stat tests to be handed to next step (just a passthrough)
	"familyVar=s"      => \$familyVar, #column name in metadata containing family id
	"groupStabilityVars=s"      => \$groupStabilityVars, #column names of categories used for calculation of resilience and persistence
	
	#SNP calling
	"minSNPDepth=i"  => \$minSNPDepth,
	"minSNPCallQual=i"  => \$minSNPCallQual,
	"skipIndels=i"     => \$noIndels,
	"SNPadaptiveQual=f" => \$useAdaptiveQual, #Default 0 (not active, recommended 0.15-0.5
	"SNPdepthFilterScale=f" => \$depthFilterScale, #Default 0.15
	"SNPindelRangeFilt=i" => \$indelRange,

) or die "Invalid strain_within.pl options\n";
die "Unexpected positional arguments: @ARGV\n" if @ARGV;
die "-GCd is required and must be a directory\n" unless length($GCd) && -d $GCd;
die "Either -MGS or -outD is required\n" unless length($MGSfile) || length($outDpre);
die "MGS file missing or empty: $MGSfile\n" if length($MGSfile) && !-s $MGSfile;
die "-MGset must be GTDB or FMG\n" unless $useGTDBmg eq "GTDB" || $useGTDBmg eq "FMG";
die "-clusterID must be between 1 and 100\n" unless $clusterID >= 1 && $clusterID <= 100;
die "Invalid subjob settings\n" if $maxSubJob < 0 || $subJob < 0 || ($maxSubJob && $subJob >= $maxSubJob);
die "Core, memory, and precompute settings must be non-negative\n"
	if $maxCores == 0 || $selfMemGb <= 0 || $preCompCons < 0;
die "-minBadLociPSmpl must be positive\n" unless $minBadLociForSampleSkip > 0;
die "-flushEvery must be positive\n" unless $appendWriteTrigger > 0;
die "Fractional filtering options must be between 0 and 1\n"
	if grep { $_ < 0 || $_ > 1 } ($multiGeneSmplMax, $conspGeneSmplMax, $GenesPerSpecies, $GeneLengthMin);
die "SNP depth, quality, adaptive filtering, and indel-range settings must be non-negative\n"
	if $minSNPDepth < 0 || $minSNPCallQual < 0 || $useAdaptiveQual < 0
		|| $depthFilterScale < 0 || $indelRange < 0;
die "-phyloMemMulti must be positive\n" unless $memMulti > 0;
die "-phyloProg must be 1 (IQ-TREE), 2 (VeryFastTree), or 3 (FastTree)\n"
	unless $phyloProg >= 1 && $phyloProg <= 3;
die "-MSAprog must be 0, 1, 2, or 4\n"
	unless grep { $MSAprog == $_ } (0, 1, 2, 4);

@subsetMGS = split /,/,$subsMGSstr if ($subsMGSstr ne "");
#print "SUBSMGS:: @subsetMGS\n";
#die timeNice(20) ." ".timeNice(12252)."\n"; #TEST

#define global vars
my $queueMode = $subMode;
$queueMode = "bash" if !$doSubmit && $queueMode eq "";
my $QSBoptHR = emptyQsubOpt($doSubmit,"",$queueMode);
my $MGSfileOri = $MGSfile; #save for later..

my $bindir;my $outD;my $scratchD;my $preConDir;my $LOGDIR;my $mapF;
my %map; my %AsGrps;my @samples;#map and assembly groups
my %ConspecificMGS; #list of conspecific MGS
my %MGSnoTree; #MGS known to have too few samples for a meaningful tree
my $legacyLocusOutputs = 0;
my %legacyLocusMGS;

my $gene2taxF; #where to find info what genes (gene cat)
my $sttime = time;	

prepRun();


my %AGlist; #list of assembly groups that need to be processed together;
createAGlist();
#foreach (sort keys %AGlist) {   print "$_ : @{$AGlist{$_}}\n";}die;

my %preCompSNPs;
my %unavailableSamples;
preComputeConsSNP();


my %replN; #my %genesWrite; #keep stats/track

#my %allFNA; my %allFAA; #big hash with all genes in @allGenes
#my %gene2genes; #no longer needed
#contains link from GCgene to fasta header assembly, cleaned up for multi copy already..
#structure: $cl2gene2{sample}{locus_id} = { member_gene => seed_catalogue_gene, ... }
#(this used to be a plain array of members plus a fully parallel %candidateSeed
#hash-of-hash-of-hash holding the same member names again just to carry the seed;
#folding the seed into the same hash removes that duplicate nesting.)
my %cl2gene2;
my $LocusByID = {}; my $MemberContext = {}; my $LocusContext = {};
my $catalogProteins = {};
my %LocusSeedProteins;
#my %SIcat;


my $SIgenes; my $Gene2COG; my $Gene2MGS; my $COGprios;
my %SIdirs; #unified storage of dirs per SI (SI==MGS)





#key step to determine with set of genes (representing MGS) is to be MSA'd for strain phylos
#these might be very limited number of genes here..
($SIgenes,$Gene2COG,$Gene2MGS,$COGprios) = readGene2tax($gene2taxF,$presortGenes,\@subsetMGS);#
#%SIgenes=%{$hr1};%Gene2COG=%{$hr2}; %Gene2MGS = %{$hr3}; %COGprios = %{$hr4};
my @specis = sort(keys(%{$SIgenes}));
$Gene2MGS = {}; #not consumed by the within-strain workflow
die "No MGS matched the selected input"
	. ($subsMGSstr ne "" ? " or -MGSsubset $subsMGSstr" : "") . "\n"
	unless @specis;
for my $MGS (@specis) {
	die "Unsafe MGS identifier '$MGS': use only letters, digits, dot, underscore, colon, plus, and hyphen\n"
		unless defined($MGS) && $MGS =~ /\A[A-Za-z0-9][A-Za-z0-9_.:+-]*\z/;
}
#sort specis by numbers, so start with MGS1, MGS2 etc
my %sis; foreach (@specis){if (m/(\d+)$/){ $sis{$_}=int($1);} else {$sis{$_}=1; print "Unknown code: $_";}}
@specis = sort {$sis{$a} <=> $sis{$b} || $a cmp $b } keys %sis;


#die "specis::\n@specis\n";
my $cnt=0; my $SaSe = "|"; 


my ($dirsNOTPrepped , $CatFileMiss , $CatNotPrepped , $treeAbsent, $doneDirs, $PhylosExist) 
			= evalFileStatus();
#DEBUG:getInputSize();


my %smplsPerMGS; #stats: MGS is represented in how many different samples?

#hashes of strings that keep results to be written for each species..
my %OCstrH ; my %OFstrH ; my %OAstrH ; my %OLstrH ;
my $splitGeneration = '';
my $splitManifest = "$LOGDIR/mainExtr.generation";
my $splitStonePrefix = "$LOGDIR/mainExtr";


if (($dirsNOTPrepped/@specis > 0.1) || $onlySubmit == 0
			|| $subJob || $redoSubmissionData
			|| $legacyLocusOutputs
			|| ($deepRepair && $dirsNOTPrepped)
			|| ($repairCAT && $CatFileMiss)){
	#$PhylosExist=0;
	
	print "\n\n----------------------------------------------------\nPart I:: extracting relevant core MGS genes (SNP consensus called) from original assemblies". "Elapsed time : ", timeNice(time - $sttime) . "\n----------------------------------------------------\n\n";
	
	prepGene2MGS();
	$Gene2COG = {}; #delete, no longer needed..
	
	reportingsMGS();
	%smplsPerMGS = (); #reporting-only sample/locus counts can be large
	$SIgenes = {}; #replaced locus selection is represented by $COGprios
	
	my @jobsMain;

	if ($maxSubJob && !$subJob){
		# A generation manifest prevents an isolated worker retry from being
		# mistaken for a complete replacement of previously merged inputs.
		clear_split_generation($splitManifest, $splitStonePrefix);
		$splitGeneration = join('.', time, $$, int(rand(1_000_000_000)));
		write_split_generation($splitManifest, $splitGeneration, $maxSubJob);
		#here needs to submit itself maxSubJob times
		my $strain1scr = getProgPaths("MGS_strain1_scr"); #self reference
		my @selfArgs = (
			'-GCd', $GCd, '-outD', $outD, '-MGS', $MGSfileOri,
			'-clusterID', $clusterID, '-submit', $doSubmit, '-onlySubmit', 0,
			'-reSubmit', 0, '-maxSubJob', $maxSubJob,
			'-MGSminGenesPSmpl', $MGStoolowGsThr,
			'-multiGeneSmplMax', $multiGeneSmplMax,
			'-conspGeneSmplMax', $conspGeneSmplMax,
			'-minBadLociPSmpl', $minBadLociForSampleSkip, '-MGSphylo', $treeFile,
			'-presortGenes', $presortGenes, '-maxGenes', $maxNGenes,
			'-flushEvery', $appendWriteTrigger,
			'-MGset', $useGTDBmg, '-redoSubmissionData', 0, '-deepRepair', 0,
			'-rmMSA', 0, '-minSNPDepth', $minSNPDepth,
			'-minSNPCallQual', $minSNPCallQual, '-forceSNPcalls', $forceVCF2FNA,
			'-preCompConsSNP', $preCompCons, '-skipIndels', $noIndels,
			'-SNPadaptiveQual', $useAdaptiveQual,
			'-SNPdepthFilterScale', $depthFilterScale,
			'-SNPindelRangeFilt', $indelRange,
		);
		push @selfArgs, ('-tmpD', $locTmpDir1) if $locTmpDir1 ne "";
		push @selfArgs, ('-MGSsubset', $subsMGSstr) if $subsMGSstr ne "";
		push @selfArgs, ('-submissionMode', $subMode) if $subMode ne "";
		my $selfCmd = $strain1scr . " " . join(" ", map { shellQuote($_) } @selfArgs);
		
		my $tmpHDD=$QSBoptHR->{tmpSpace} ; $QSBoptHR->{tmpSpace} =15; #request some basic amount
		
		#submission of self-subjobs..
		for (my $sj = 1; $sj < $maxSubJob; $sj ++){
			my $cmdX = "$selfCmd -subjob $sj;\n";
			my $checkF = "$LOGDIR/mainExtr.${sj}.stone";
			$cmdX .= "printf '%s\\n' ".shellQuote($splitGeneration)
				." > ".shellQuote($checkF)."\n";
			#die "$cmdX\n\n";
			print $LOGDIR."Strain1_B${sj}.sh\n";
			my ($dep,$qcmd) = qsubSystem($LOGDIR."Strain1_B${sj}.sh",$cmdX,1,"${selfMemGb}G","Str1.$sj","","",1,[],$QSBoptHR);
			push(@jobsMain,$dep);
		}
		$QSBoptHR->{tmpSpace} = $tmpHDD;
	}
	
	#and extract the corresponding fna/ faa from every other dir.. main single core work
	#this will also determine how many genes per MGS are now extracted..
	extractFNAFAA2genes();#@allGenes);
	%cl2gene2 = (); #no longer needed, delete
	$LocusByID = {};
	$MemberContext = {};
	$LocusContext = {};
	%LocusSeedProteins = ();
	$COGprios = {};
	#write logs to found genes etc.
	writeLogsStep1();
	write_worker_completion("$splitStonePrefix.0.stone", $splitGeneration)
		if $maxSubJob && !$subJob;
	
	if ($subJob){
		print "Finished subJob ${subJob}/$maxSubJob. Exiting..\n";
		exit(0);
	}

	if ($maxSubJob && !$subJob){ # second part for main worker: check that everything else is finished..
		if (@jobsMain && !$doSubmit) {
			print "Split-worker scripts were generated but not submitted; stopping before incomplete outputs are combined.\n";
			exit(0);
		}
		qsubSystemJobAlive( \@jobsMain,$QSBoptHR ) if @jobsMain && $doSubmit;
		die "Split extraction generation is incomplete; refusing to merge worker subsets\n"
			unless split_generation_complete($splitManifest, $splitStonePrefix, $maxSubJob);
		mergeConspecificLogs();
		
		#combineMGSgenes();
	}
	
	print "\nGene extraction & redistribution finished, ready to proceed to phylogeny jobs\n";

} else {
	print "Skipping Part I, outdir already prepared.\n";
}

#die;


#load some log files..
#if (scalar(keys(%genesWrite)) == 0) { #load genes found..
#	#read logs of found genes etc.
#	foreach my $MGS (@specis){
#		my $outD2 = $SIdirs{$MGS}; my $llogF="$outD2/geneFnd.log";
#		next unless (-e $llogF);
#		my $Lstr = `cat $llogF`; $Lstr =~ m/Total genes write (\S+): (\d+)/; 
#		$genesWrite{$1} = $2;
#		die "$llogF incorrect: $1 != $MGS\n" if ($1 ne $MGS);
#		$PhylosExist =0 if (!-d "$outD2/pjylo/");
#	}
#}
if (scalar(keys(%ConspecificMGS)) == 0){
	my $conlog = "$LOGDIR/ConspecificMGS.log";
	my $legacy_conlog = "$bindir/LOGandSUB/ConspecificMGS.log";
	$conlog = $legacy_conlog if !-s $conlog && -s $legacy_conlog;
	if (-s $conlog) {
		open I,"<$conlog" or die "Can't open conspecific $conlog\n";
		while (my $l = <I>){my @spl = split /\t/,$l;$ConspecificMGS{$spl[0]} = [split(/,/,$spl[1])];}
		close I;
	} else {
		warn "No prior conspecific-sample log found at $conlog; continuing without historical exclusions\n";
	}
}




my $FNAref = {}; my $FAAref = {};
my $SIgenes_OG = {}; my %OGgenesByCOG;
my %outgroupGeneCache;

my $geneCatLoaded=0;
#read in genecat to create outgroup fasta sequences..
if ($CatNotPrepped || $treeAbsent || $repairCAT || $deepRepair || $dirsNOTPrepped || $onlySubmit == 0 || $redoSubmissionData == 1){
	#also read reference gene seqs (for outgroup)
	my $refFNA = ""; my $refFAA = ""; my $refNameL = "unknw";
	if ($mode eq "MGS" || $mode eq "MGSall"){
		print "Reading reference genecat genes, to create outgroup sequences\n";
		$refFNA = "$GCd/compl.incompl.$clusterID.fna"; $refFAA = "$GCd/compl.incompl.$clusterID.prot.faa";
		$geneCatLoaded=1;$refNameL = "geneCat";
	} elsif ($mode eq "FMG"){
		print "reading FMG ref genes..";
		$refFNA = "$GCd/FMG/COG*.fna"; $refFAA = "$GCd/FMG/COG*.faa";
		$refNameL = "FMG ref";
	}
	
	# Outgroups can lie outside an explicitly requested target subset.
	my ($hr1,$Gene2COG_OG,$hr3,$hr4) = readGene2tax($gene2taxF,$presortGenes,[]);
	$SIgenes_OG = $hr1;
	for my $MGS (keys %{$hr4}) {
		for my $locus (@{$hr4->{$MGS}}) {
			my $gene = $hr1->{$MGS}{$locus};
			next unless defined($gene) && defined($Gene2COG_OG->{$gene});
			push @{$OGgenesByCOG{$MGS}{$Gene2COG_OG->{$gene}}}, $gene;
		}
	}
	#%SIgenes_OG=%{$hr1}; my %Gene2COG_OG=%{$hr2}; 
	$FAAref = readFasta($refFAA,1,"\\s",$Gene2COG_OG);
	$FNAref = readFasta($refFNA,1,"\\s",$Gene2COG_OG);
	print "read ". scalar(keys %{$FNAref})." genes from $refNameL\n";
	print "done\n";
}



print "\n\n----------------------------------------------------\n";
print "Part II:: resort .cat files, submit intraStrain phylogenies for " . scalar(@specis) . " MGS. ". "Elapsed time : ", timeNice(time - $sttime) ."\n----------------------------------------------------\n\n";


die "Tree for outgroup specified, but file not found:$treeFile\nAborting..\n" if  ($treeFile ne "" && !-e $treeFile);

#sort by largest dir first..
my @sizeOfDirs = getInputSize();
my @idx = sort { $sizeOfDirs[$b] <=> $sizeOfDirs[$a] } 0 .. $#sizeOfDirs;
@specis=@specis[@idx];@sizeOfDirs=@sizeOfDirs[@idx];
#print "SIZE2:: $sizeOfDirs[0] $sizeOfDirs[1] $specis[0] $specis[1]\n"; die;


#die;
#go through every SpecI;
$cnt=0; my $lcnt=-1; my @jobs; my %expectedTreeOutputs; my $Nspecis = @specis;
foreach my $MGS (@specis){ #loop creates per specI file structure to run buildTreeScript on..
	$lcnt++;
	if (!$reSubmit && !$repairCAT && !$redoSubmissionData && $CatFileMiss==0 && $CatNotPrepped==0 && $treeAbsent ==0){
		print "\nAll submission dirs prepared, nothing to do..\n";
		last;
	}
	if (exists $MGSnoTree{$MGS}) {
		limitedNotice('MGS skipped after too-few-samples extraction',
			"Skipping $MGS: previous extraction found too few samples for a tree.\n");
		next;
	}
	# previous condition was too lax: ( ($CatNotPrepped/$#specis) < 0.1)  , just check if we can resubmit anything here..
	if (exists($ConspecificMGS{$MGS}) && $ConspecificMGS{$MGS}->[0] =~ m/multicopy/){
		limitedNotice('MGS skipped as conspecific or multicopy',
			"Skipping $MGS due to inclusion in conspecific MGS list.\n");next;
	}
	if ($startSubFromMGS ne "" ){
		if ($MGS ne $startSubFromMGS){next;
		} else { $startSubFromMGS = "";} #deactivate now
	}
	my $outD2 = $SIdirs{$MGS};
	my $treeStone = "$outD2/treeDone.sto";
	my $IQtreef= "$outD2/phylo/IQtree_allsites.treefile";
	$IQtreef = "$outD2/phylo/VERYFASTTREE_allsites.nwk" if ($phyloProg == 2);
	$IQtreef = "$outD2/phylo/FASTTREE_allsites.nwk" if ($phyloProg == 3);
	
	if (!$reSubmit && !$repairCAT && !$redoSubmissionData && !exists($legacyLocusMGS{$MGS})
			&& -e $treeStone && -s $IQtreef ){
		limitedNotice('MGS skipped with existing trees',
			"Skipping $MGS: a valid tree already exists.\n");
		next;
	}
	
	print "Processing $MGS (".($lcnt + 1)."/$Nspecis); elapsed ".timeNice(time - $sttime)."\n";
	my $inputFNAsize = $sizeOfDirs[$lcnt];
	#PART I: create fasta files required by tree
	make_path($outD2) unless -d $outD2;
	my $tmpD  = "$scratchD/outs/$MGS/";
	if ($inputFNAsize ==0){
		limitedNotice('MGS skipped with empty input', "Skipping $MGS: input is empty.\n");
		next;
	} #empty input
	my $publishedInputsReady = fileGZe("$outD2/$FNAstdof")
		&& fileGZe("$outD2/$FAAstdof")
		&& fileGZe("$outD2/$CATstdof");
	my $mustRegenerateInputs = $repairCAT || $deepRepair || $redoSubmissionData
		|| exists($legacyLocusMGS{$MGS});
	if ($publishedInputsReady && !$mustRegenerateInputs) {
		print "  Recovery input: using complete published FNA/FAA/category files\n";
	} else {
		unless (combineMGSgenesDir($MGS,$tmpD,$tmpD)) {#$outD2); -> keep in tmpdir for now..
			limitedWarn('MGS with incomplete combined worker input',
				"$MGS has neither complete published inputs nor complete combined worker input; leaving it for an extraction repair run\n");
			next;
		}
	}
	
	#final locations (after copying etc)
	my $FNAtf = "$outD2/$FNAstdof"; my $FAAtf = "$outD2/$FAAstdof";
	my $CATtf = "$outD2/$CATstdof"; #my $Linkf = "$outD2/$LINKstdof";
	my $MSAdir = "$outD2/MSA/";
	
	
	my $outgS = "";my $OG = "";
	if (fileGZe("$outD2/data.log")) {
		my ($log_fh) = gzipopen("$outD2/data.log", "outgroup log");
		$OG = <$log_fh> // "";
		close $log_fh;
		chomp $OG;
		$OG =~ s/^OG://;
	}
	# buildTree5 validates its persistent checkpoints and restarts invalid stages.
	my $contPhylo = 1;
	
	#main command to build within species strain tree.. missing outgroup so far ($outgS)
	
	#fileGZs($FNAtf) / (1024 * 1024); #size in MB
	#$inputFNAsize*=5 if ($FNAtf =~ m/\.gz$/); #account for compressed input
	if ( 0&& ($MSAprog==4 && $inputFNAsize>700) ){ $QSBoptHR->{useLongQueue} = 1 ;	}
	my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0";
	my $baseMemMult = 150; $baseMemMult = 50 if ($phyloProg ==3 || $phyloProg ==2);
	my $totMem = int($inputFNAsize *$baseMemMult * $memMulti);
	$totMem = 10000*$memMulti if ($totMem < 10000);$totMem = 220000*$memMulti if ($totMem > 220000);
	my $numCoreL = $numCores;	
	if ($maxCores >0){ #scale cores according to used memory size
		$numCoreL = int($maxCores * sqrt($inputFNAsize/$sizeOfDirs[0]));
		$numCoreL = 4 if ($numCoreL < 4);		$numCoreL = $maxCores if ($numCoreL > $maxCores);
	}
	
	my $subsSmpl = -1;my $useSuperTree = 0;
	my $bts = getProgPaths("buildTree_scr");
	my $treeFlag = "-runIQtree 1 "; 
	if ($phyloProg == 2){$treeFlag = "-runVeryFastTree 1 ";}if ($phyloProg == 3){$treeFlag = "-runFastTree 1 ";}
	my $tree_sample_separator = quotemeta($SaSe);
	my $Tcmd= "$bts -fna ".shellQuote($FNAtf)." -aa ".shellQuote($FAAtf)." -smplSep ".shellQuote($tree_sample_separator)." -cats ".shellQuote($CATtf)." -outD ".shellQuote($outD2)." $treeFlag -cores $numCoreL  ";
	$Tcmd .= "-AAtree 0 -bootstrap 0 -NTfiltCount 400 -NTfilt 0.07 -NTfiltPerGene $GeneLengthMin -GenesPerSpecies $GenesPerSpecies -runRaxMLng 0 -minOverlapMSA 2 ";
	$Tcmd .= "-subsetSmpls $subsSmpl -fracMaxGenes90pct 0.7 "; #concentrate on almost complete gene groups.. can yield more samples overall and speeds up calc..
	$Tcmd .= "-rmMSA $rmMSA -gzInput 1 "; #save more diskspace..
	$Tcmd .= "-SynTree 0 -NonSynTree 0 -MSAprogram $MSAprog -continue $contPhylo -AutoModel 0 -iqFast 1 -superTree $useSuperTree ";
	$Tcmd .= "-runDNDS 0 -runTheta 0 -tmpD ".shellQuote("$scratchD/$MGS/")." -map ".shellQuote($mapF)." ";
	my $postCmd = "\n\ntest -s ".shellQuote($IQtreef)."\ntouch ".shellQuote($treeStone)."\n";
		#die "$cmd\n" if ($cnt ==10);
	
	#if (!fileGZe($FNAtf) || !fileGZe($FAAtf) ||  ( !fileGZe($CATtf) && !-e "$CATtf.tmp") ){
	#	print "Can't find required input files:\n$FNAtf\n$FAAtf\n$CATtf\n" ;
	#	die if ($cnt <= 1);
	#	next;
	#}
	
	qsubSystemWaitMaxJobs($checkMaxNumJobs,0,$QSBoptHR) if $doSubmit;
	#reformat .cat.tmp -> .cat and add outgroup fna seqs
	my $multiSmpl;my $ngenes; my $needsCopy = 0; my $inputReady = 0;
	($multiSmpl,$ngenes,$OG,$needsCopy,$inputReady)=
		addOutgroup2MGS($MGS,$OG,$tmpD); #$outD2 $tmpD
	# Locus names are MGS-qualified, so cached outgroup choices have no reuse
	# after this MGS and would otherwise accumulate for the entire submission.
	%outgroupGeneCache = ();
	unless ($inputReady) {
		$QSBoptHR->{tmpSpace} = $tmpSHDD;
		$QSBoptHR->{useLongQueue} = 0;
		limitedNotice('MGS awaiting input repair',
			"$MGS: input files are not ready; leaving it unmarked so a repair run can retry it.\n");
		next;
	}
	
	$outgS = " -outgroup ".shellQuote($OG)." "  if ($OG ne "");
	my $preCmd = "";
	if ($needsCopy){
		$preCmd = treeInputPrecopyCommand(
			$tmpD, $outD2, $numCoreL, $FNAtf, $FAAtf, $CATtf
		);
	}

	if ($multiSmpl>2){
		print "  Tree input: $multiSmpl samples, $ngenes potential genes; $numCoreL cores, $totMem memory\n";
	} else {
		limitedNotice('MGS with too few samples for tree statistics',
			"$MGS: too few samples ($multiSmpl) for tree statistics\n");
		writeTooFewMarker($outD2, $multiSmpl, $ngenes);
		remove_tree($tmpD) if $needsCopy && -d $tmpD;
		$QSBoptHR->{tmpSpace} = $tmpSHDD;
		$QSBoptHR->{useLongQueue} = 0;
		next;
	}
	unlink "$outD2/tooFewSamples.sto" if -e "$outD2/tooFewSamples.sto";
	
	#PART II: qsub tree build command
	
	#die "$cmd\n" if ($cnt ==10);
	if ($doSubmit) {
		unlink $treeStone or die "Cannot remove stale tree checkpoint $treeStone: $!\n"
			if -e $treeStone;
		unlink $IQtreef or die "Cannot remove stale tree output $IQtreef: $!\n"
			if -e $IQtreef;
	}
	my ($dep,$qcmd) = qsubSystem($outD2."treeCmd.sh",$preCmd.$Tcmd.$outgS.$postCmd,$numCoreL,int($totMem) ."M","FT$cnt","","",1,[],$QSBoptHR);
	$QSBoptHR->{tmpSpace} =$tmpSHDD;
	$QSBoptHR->{useLongQueue} = 0;
	$cnt ++;
	push (@jobs,$dep) if defined($dep) && length($dep);
	$expectedTreeOutputs{$MGS} = [$IQtreef, $treeStone];
	#die $outD2."treeCmd.sh\n";

}
if ($maxSubJob
		&& split_generation_complete($splitManifest, $splitStonePrefix, $maxSubJob)
		&& !splitWorkerPartsRemain()) {
	clear_split_generation($splitManifest, $splitStonePrefix);
}
#too many jobs to use as job dependency..
qsubSystemJobAlive( \@jobs,$QSBoptHR ) if @jobs && $doSubmit;
if ($doSubmit) {
	my @failed = grep {
		my ($tree, $stone) = @{$expectedTreeOutputs{$_}};
		!-s $tree || !-e $stone;
	} sort keys %expectedTreeOutputs;
	die "Tree jobs completed without valid tree outputs for: ".join(",", @failed)."\n"
		if @failed;
}
print "\nAll done for $cnt Bins\nRun strain_within_2.pl for summary stats:\n";

my $outDX =  $MGSfile;#"$GCd/$mode/intra_phylo/";
$outDX =~ s/[^\/]+$//;
my $MGSabundance = "$GCd/Anno/Tax/GTDBmg_MGS/specI.mat";
$MGSabundance = "$bindir/Annotation/Abundance/MGS.matL7.txt";

my $strain2Scr = getProgPaths("MGS_strain2_scr");

my $nxtCmd = "$strain2Scr -GCd ".shellQuote($GCd)." -FMGdir ".shellQuote($outD)." -MGSmatrix ".shellQuote($MGSabundance)." -cores 4 -reSubmit 0 -DiscTests ".shellQuote($discTests)." -ContTests ".shellQuote($contTests)." -familyVar ".shellQuote($familyVar)." -groupStabilityVars ".shellQuote($groupStabilityVars)." ";
$nxtCmd .= "-submit $doSubmit ";
$nxtCmd .= "-qsubSystem ".shellQuote($subMode)." " if $subMode ne "";
$nxtCmd .= "-Hcores $maxCores " if $maxCores > 0;
if ($mapF2 eq ""){$nxtCmd .= "-map ".shellQuote($mapF)." ";} else {$nxtCmd .= "-map ".shellQuote($mapF2)." ";}

$nxtCmd .= "\n";

#$GCd/MB2.clusters.ext.can.Rhcl.matL0.txt
	my ($dep,$qcmd) = qsubSystem($outD."/strainAnalysis2.sh",$nxtCmd,1,"60G","2StrainSub","","",1,[],$QSBoptHR);
print "\n". $nxtCmd."\n";


#cleanup
remove_tree($locTmpDir) if -d $locTmpDir;
remove_tree($preConDir) if ($preCompCons && -d $preConDir);


exit(0);

 

#########################################################################################
#########################################################################################


#	combineMGSgenesDir($MGS,$outD2);
sub combineMGSgenesDir{
	my ($MGS,$tmpD,$outD2) = @_;
	my @required = ("$outD2/$FNAstdof", "$outD2/$FAAstdof", "$outD2/$CATstdof.tmp");

	#my $outD3 = $tmpD; #work locally, copy later..
	my @filesets = (
		[$FNAstdof,      "$tmpD/$FNAstdof",      "$tmpD/$FNAstdof"],
		[$FAAstdof,      "$tmpD/$FAAstdof",      "$tmpD/$FAAstdof"],
		[$LINKstdof,     "$tmpD/$LINKstdof",     "$tmpD/$LINKstdof"],
		["$CATstdof.tmp","$tmpD/$CATstdof.tmp",  "$tmpD/$CATstdof.tmp"],
	);
	my %partsByName;
	my $workerCount = $maxSubJob || 1;
	for my $set (@filesets) {
		my ($name, $prefix) = @$set;
		# Exact suffix matching deliberately excludes abandoned
		# "$outfile.merge.PID" scratch files from worker input.
		$partsByName{$name} = [exact_worker_parts($prefix, $workerCount)];
	}
	my $hasFreshParts = grep { @{$partsByName{$_}} }
		($FNAstdof, $FAAstdof, $LINKstdof, "$CATstdof.tmp");
	my $aggregateComplete = !grep { !fileGZe($_) } @required;
	return $aggregateComplete unless $hasFreshParts;
	if ($maxSubJob && !split_generation_complete($splitManifest, $splitStonePrefix, $maxSubJob)) {
		limitedWarn('partial worker retries without a complete generation',
			"Ignoring partial worker retry for $MGS: no complete matching split-extraction generation is present\n");
		return $aggregateComplete;
	}
	for my $requiredName ($FNAstdof, $FAAstdof, "$CATstdof.tmp") {
		unless (@{$partsByName{$requiredName}}) {
			limitedWarn('worker extractions missing required parts',
				"Fresh worker extraction for $MGS lacks required $requiredName parts; retaining parts for repair\n");
			return 0;
		}
	}

	my @consumedParts;
	for my $set (@filesets) {

		my ($name, $prefix, $outfile) = @$set;
		my @parts = @{$partsByName{$name}};
		next unless @parts;

		my $mergeFile = "$outfile.merge.$$";
		open my $out, ">", $mergeFile or die "Cannot create $mergeFile: $!\n";
		binmode $out;

		for my $file (@parts) {
			open my $in, "<", $file or die "Cannot read $file: $!";
			binmode $in;
			while (my $line = <$in>) {
				print {$out} $line or die "Cannot write $mergeFile: $!\n";
			}
			close $in or die "Cannot close $file: $!\n";
		}

		close $out or die "Cannot close $mergeFile: $!\n";
		rename $mergeFile, $outfile or die "Cannot replace $outfile: $!\n";
		push @consumedParts, @parts;
	}
	my $complete = !grep { !fileGZe($_) } @required;
	if ($complete) {
		for my $part (@consumedParts) {
			unlink $part or warn "Cannot remove combined part $part: $!\n";
		}
	} elsif (@consumedParts) {
		limitedWarn('combined MGS inputs missing source parts',
			"Incomplete combined input for $MGS; retaining all source parts for repair\n");
	}
	if ($outD2 ne $tmpD) {
		make_path($outD2);
		for my $source (grep { -f $_ } glob("$tmpD/*")) {
			copy($source, "$outD2/" . basename($source))
				or die "Cannot copy $source to $outD2: $!\n";
		}
	}
	return $complete;
}

sub splitWorkerPartsRemain {
	return 0 unless $maxSubJob;
	for my $mgs_dir (grep { -d $_ } bsd_glob("$scratchD/outs/*")) {
		for my $name ($FNAstdof, $FAAstdof, $LINKstdof, "$CATstdof.tmp") {
			return 1 if exact_worker_parts("$mgs_dir/$name", $maxSubJob);
		}
	}
	return 0;
}


	
sub locusParts {
	my ($locus, $default_mgs) = @_;
	my @parts = split /\|/, ($locus // ''), -1;
	return @parts if @parts == 3;
	return ('', '', '') if @parts > 3;
	return ($default_mgs // '', $parts[0] // '', $parts[1] // '') if @parts == 2;
	return ($default_mgs // '', $parts[0] // '', '');
}

sub externalLocusName {
	my ($locus, $default_mgs) = @_;
	my (undef, $cog, $primary_gene) = locusParts($locus, $default_mgs);
	die "Cannot create an external name for malformed locus '$locus'\n"
		unless length($cog) && length($primary_gene);
	return join($SaSe, $cog, $primary_gene);
}

sub internalLocusName {
	my ($locus, $default_mgs) = @_;
	my ($mgs, $cog, $primary_gene) = locusParts($locus, $default_mgs);
	return '' unless length($mgs) && length($cog) && length($primary_gene);
	return '' if defined($default_mgs) && length($default_mgs) && $mgs ne $default_mgs;
	return join($SaSe, $mgs, $cog, $primary_gene);
}

sub outgroupGeneForLocus {
	my ($outgroup, $locus, $default_mgs) = @_;
	my $cache_key = join("\t", $outgroup, $locus);
	return $outgroupGeneCache{$cache_key} if exists $outgroupGeneCache{$cache_key};
	my (undef, $cog, $primary_gene) = locusParts($locus, $default_mgs);
	my @candidates = @{$OGgenesByCOG{$outgroup}{$cog} || []};
	return $outgroupGeneCache{$cache_key} = '' unless @candidates;
	my $target_sequence = $FAAref->{$primary_gene} // $catalogProteins->{$primary_gene};
	unless (defined($target_sequence) && length($target_sequence)) {
		return $outgroupGeneCache{$cache_key} = $candidates[0];
	}
	my ($best_gene, $best_score) = ('', -1);
	for my $candidate (@candidates) {
		next unless defined($FAAref->{$candidate}) && length($FAAref->{$candidate});
		my $length_ratio = length($target_sequence) < length($FAAref->{$candidate})
			? length($target_sequence) / length($FAAref->{$candidate})
			: length($FAAref->{$candidate}) / length($target_sequence);
		next if $length_ratio < 0.5;
		my $score = protein_kmer_similarity($target_sequence, $FAAref->{$candidate});
		if ($score > $best_score || ($score == $best_score && ($best_gene eq '' || $candidate cmp $best_gene) < 0)) {
			($best_gene, $best_score) = ($candidate, $score);
		}
	}
	return $outgroupGeneCache{$cache_key} = $best_gene;
}

sub addOutgroup2MGS{
	my ($MGS,$OG,$tmpD) = @_;
	my $outD2 = $SIdirs{$MGS};
	my $outD3 = $tmpD;
	my $outputReady = fileGZe("$outD2/$FNAstdof")
		&& fileGZe("$outD2/$FAAstdof") && fileGZe("$outD2/$CATstdof");
	if ($outputReady && !$repairCAT && !$deepRepair && !$redoSubmissionData
			&& !exists($legacyLocusMGS{$MGS})){
		my %samples_seen;
		my $genes_seen = 0;
		my ($cat_fh) = gzipopen("$outD2/$CATstdof", "existing category file");
		while (my $line = <$cat_fh>) {
			chomp $line;
			my @entries = split /\t/, $line;
			$genes_seen++;
			for my $entry (@entries) {
				my ($sample) = split /\Q$SaSe\E/, $entry, 2;
				$samples_seen{$sample} = 1 if defined $sample && length $sample;
			}
		}
		close $cat_fh or die "Cannot close existing category file for $MGS: $!\n";
		return(scalar(keys %samples_seen),$genes_seen,$OG,0,1);
	}
	my $temporaryInput = fileGZe("$tmpD/$FNAstdof") && fileGZe("$tmpD/$FAAstdof")
		&& (fileGZe("$tmpD/$CATstdof.tmp") || fileGZe("$tmpD/$CATstdof"));
	if (exists($legacyLocusMGS{$MGS}) && !$temporaryInput) {
		limitedWarn('MGS with stale identifiers and no regenerated input',
			"$MGS has stale sequence identifiers but no regenerated temporary input; skipping it until extraction can be rerun\n");
		return(0, 0, $OG, 0, 0);
	}
	$outD3 = $outD2 if !$temporaryInput && $outputReady;
	my $FNAtf = "$outD3/$FNAstdof"; my $FAAtf = "$outD3/$FAAstdof";
	my $CATtf = "$outD3/$CATstdof"; #my $Linkf = "$outD3/$LINKstdof";
	#my $IQtreef= "$outD3/phylo/IQtree_allsites.treefile";
	my $rmCatTmp=0;
	my $MSAdir = "$outD3/MSA/";
	die "Gene cat wasn't loaded, check program logic.\n!$deepRepair && $redoSubmissionData == 0 && $onlySubmit==1 && !$dirsNOTPrepped && !-e $CATtf.tmp \n" if (!$geneCatLoaded);
	unless (fileGZe($FNAtf) && fileGZe($FAAtf)) {
		limitedWarn('MGS missing NT or AA input', "Missing NT or AA input for $MGS in $outD3\n");
		return(0, 0, $OG, 0, 0);
	}
	my %SIcatLoc;
	my $malformedCatEntries = 0;
	if (fileGZe( "$CATtf.tmp") && !fileGZe( "$CATtf")){
		my ($ICT,$status) = gzipopen("$CATtf.tmp","Can't open cat file $CATtf.tmp\n",0);
		#open ICT,"<$CATtf.tmp" or die "Can't open cat file $CATtf.tmp\n";
		while (<$ICT>){
			chomp; my @spl = split /\t/;
			if (@spl < 4){
				limitedWarn('malformed temporary category rows',
					"Malformed category row in $CATtf: $_\n");
				$malformedCatEntries++;
				next;
			}
			die "$_\n$spl[0] not eq $MGS\n" unless ($spl[0] eq $MGS);
			$SIcatLoc {$spl[1]} {$spl[2]} = $spl[3];
		}
		close $ICT;
		$rmCatTmp = 1;
	} elsif (fileGZe( $CATtf)) {
		#print OC $SIcatLoc{$cog}{$smpl};				print OC "\t".$SIcatLoc{$cog}{$smpl};		my $ng = "$OG$SaSe$cog";
		#print "Reconstructing tmp cat file.";
		my ($ICT,$startus) = gzipopen($CATtf,"Can't open (precompiled) cat file $CATtf\nConsider deleting strain dir and rerunning strainMGS script\n",1);
		#open ICT,"<$CATtf" or die "Can't open (precompiled) cat file $CATtf\nConsider deleting strain dir and rerunning strainMGS script\n";
		my $catLines=0;my $cntItems=0;
		#  $repairCAT .. auto implemented..
		while (<$ICT>){
			chomp; my @spl = split /\t/;
			foreach my $tags (@spl){
				#my @spl2 = split (/\\$SaSe/,$tags);
				#$SIcatLoc {$spl2[1]}{$spl2[0]}  = $tags;print "$spl2[1] : $spl2[0]  = $tags\n";
				my ($sample, $external_locus) = split /\Q$SaSe\E/, $tags, 2;
				my $locus = defined($external_locus) ? internalLocusName($external_locus, $MGS) : '';
				if (defined($sample) && length($sample) && length($locus)){
					$SIcatLoc {$locus}{$sample}  = $tags;
				} else {
					$malformedCatEntries++;
				}
				$cntItems++;
				#print "$2 $1  = $tags\n";
			}
			$catLines++;
		}
		close $ICT;
		if ($catLines != keys (%SIcatLoc)){ #redo MSA
			remove_tree($MSAdir) if -d $MSAdir;
		}
		print "  Category input: $catLines lines, ".scalar(keys %SIcatLoc)
			." loci, $cntItems entries ($CATtf)\n";
	} else {
		limitedWarn('MGS missing category inputs',
			"$MGS has neither .cat nor .cat.tmp input in $outD3\n");
		return(0, 0, $OG, 0, 0);

	}
	
	
	#my @curCogs = sort keys %{$SIcat{$MGS}};
	my @curCogs = sort keys %SIcatLoc;
	if (scalar(@curCogs) < 10){
		if (!@curCogs && $malformedCatEntries) {
			limitedWarn('MGS with malformed category input',
				"$MGS category input is malformed; leaving it unmarked for repair\n");
			return(0, 0, $OG, 0, 0);
		}
		my %fewSamples;
		for my $cog (@curCogs) {
			$fewSamples{$_} = 1 for keys %{$SIcatLoc{$cog}};
		}
		limitedWarn('MGS with too few usable genes for tree construction',
			"$MGS has only ".scalar(@curCogs)." usable genes; skipping tree construction\n");
		return(scalar(keys %fewSamples), scalar(@curCogs), $OG, $outD3 ne $outD2, 1);
	}
	#print "COGs: $curCogs[0] $curCogs[1]\n";
	
	# --------------------------- OUTGROUP ----------------------------------------
	#include outgroup?
	if ($treeFile ne ""){
		my $neiTree = getProgPaths("neighborTree");
		my $call = "$neiTree ".shellQuote($treeFile)." ".shellQuote($MGS);
		#print "$call\n";
		my $OG1 = `$call`;
		if ($? != 0) {
			limitedWarn('outgroup lookup command failures',
				"Can't find outgroup from call $call; building an ingroup-only tree\n");
			$OG1 = "";
		}
		chomp $OG1;
		my @sspl = grep { length } split /\s+/,$OG1; $OG = "";
		limitedWarn('MGS without outgroup candidates',
			"No outgroup candidates returned for $MGS; building an ingroup-only tree\n")
			if @sspl == 0;
		my $cntShrCogs=0;
		for my $candidate (@sspl) {
			$cntShrCogs = 0;
			$OG = $candidate;
			if (!exists($SIgenes_OG->{$OG})){
				next;
			}
			#$cntShrCogs=0;
			#if (exists($SIgenes_OG{$OG})){
			foreach my $cog (@curCogs){
				my (undef, $annotation) = locusParts($cog, $MGS);
				next if $annotation =~ m/^uniq\d+$/;
				my $outgroup_gene = outgroupGeneForLocus($OG, $cog, $MGS);
				next unless length($outgroup_gene) && exists($FNAref->{$outgroup_gene});
				$cntShrCogs ++;
			}
			last if $cntShrCogs >= 10;
		}
		if ($cntShrCogs < 10){
			my @locus_preview = @curCogs[0 .. ($#curCogs < 9 ? $#curCogs : 9)];
			limitedWarn('MGS without a sufficiently represented outgroup',
				"Could not find a sufficiently represented outgroup for $MGS; candidates: @sspl; loci: @locus_preview\n");
			$OG = "";
		}
		if ($OG ne "" && !exists($SIgenes_OG->{$OG})){
			limitedWarn('selected outgroups absent from gene catalogue',
				"Selected outgroup $OG for $MGS is absent from the gene catalogue\n");
			$OG="";
		} 
	print "  Using outgroup $OG\n" if $OG ne '';
		#next;
	}
	
	#and fasta/faa/cat files..
	#open OL,">$Linkf" or die "Can't open link file $Linkf\n";
	#append to FNA/FAA ..for outgroups
	
	my %uniqSmpls;my $OGgenesUsed=0;
	my $tmpFAAog = ""; my $tmpFNAog = "";
	my $resolvedFNA = resolve_fasta_artifact($FNAtf);
	my $resolvedFAA = resolve_fasta_artifact($FAAtf);
	# Only identifiers are needed for duplicate detection.  Loading all
	# sequences here used to duplicate both complete per-MGS FASTA files.
	my $existingFNA = length($resolvedFNA) ? readFastaIDs($resolvedFNA) : {};
	my $existingFAA = length($resolvedFAA) ? readFastaIDs($resolvedFAA) : {};
	foreach my $cog (@curCogs){
		if ($OG ne ""){
			my (undef, $annotation) = locusParts($cog, $MGS);
			my $geneKey = outgroupGeneForLocus($OG, $cog, $MGS);

			next unless length($geneKey) && exists($FNAref->{$geneKey}) && exists($FAAref->{$geneKey});
			next if ($annotation =~ m/^uniq\d+$/);
			my $ng = "$OG$SaSe" . externalLocusName($cog, $MGS);
			$tmpFNAog .= ">$ng\n$FNAref->{$geneKey}\n" unless exists $existingFNA->{$ng};
			$tmpFAAog .= ">$ng\n$FAAref->{$geneKey}\n" unless exists $existingFAA->{$ng};
			#$SIcat{$MGS}{$cog}{$OG} = $ng;
			$SIcatLoc{$cog}{$OG} = $ng;
			$OGgenesUsed++;
			#if ($new){ print OC "$ng";$new=0;
			#} else { print OC "\t$ng"; }
		}
	}
	append_fasta_records_atomic($FNAtf, $tmpFNAog);
	append_fasta_records_atomic($FAAtf, $tmpFAAog);
	
	#print "used $OGgenesUsed genes  ";
	my $cat_write = "$CATtf.write.$$";
	open my $cat_out, '>', $cat_write or die "Can't open temporary cat file $cat_write: $!\n";
	foreach my $cog (@curCogs){
		my $cntL=0;
		foreach my $smpl (sort keys %{$SIcatLoc{$cog}}){
			print {$cat_out} ($cntL ? "\t" : ""), $SIcatLoc{$cog}{$smpl}
				or die "Can't write temporary cat file $cat_write: $!\n";
			$cntL++;
			$uniqSmpls{$smpl} = 1;
		}
		print {$cat_out} "\n" or die "Can't write temporary cat file $cat_write: $!\n";
	}
	close $cat_out or die "Can't close temporary cat file $cat_write: $!\n";
	rename $cat_write, $CATtf or die "Can't replace cat file $CATtf: $!\n";
	print "  Generated category file for $MGS\n";
	if ($OGgenesUsed ==0 && $OG ne ""){
		limitedWarn('MGS with no usable outgroup genes',
			"Couldn't include any outgroup genes for $OG; building $MGS without an outgroup\n");
		$OG = "";
	}

	#note done somewhere how many genes these actually are..
	my $multiSmpl=0;
	$multiSmpl = scalar(keys %uniqSmpls);
	#system "rm -f $CATtf.tmp*\n" if (fileGZe("$CATtf.tmp"));
	#system "echo \"OG:$OG\" > $outD3/data.log";
	for my $stale_cat (glob("$CATtf.tmp*")) {
		unlink $stale_cat or die "Cannot remove $stale_cat: $!\n";
	}
	open my $log, '>', "$outD3/data.log" or die "Cannot create $outD3/data.log: $!\n";
	print $log "OG:$OG\n";
	close $log or die "Cannot close $outD3/data.log: $!\n";

	
	
	#if ($outD3 ne $outD2){system "cp $outD3/* $outD2;";
	#local? -> no, give to slurm job..
	
	my $needsCopy = $outD3 ne $outD2 ? 1 : 0;

	return ($multiSmpl,scalar(@curCogs),$OG,$needsCopy,1);
	
# --------------------------- OUTGROUP DONE ----------------------------------------
}



sub writeLogsStep1{


	#print log file
	my $conlog = $maxSubJob
		? "$LOGDIR/ConspecificMGS.$subJob.log"
		: "$LOGDIR/ConspecificMGS.log";
	make_path($LOGDIR) unless -d $LOGDIR;
	open LO,">$conlog" or die "Can't open conspecific log file: $conlog\n";
	foreach my $MGS (sort keys %ConspecificMGS){
		my %seen;
		my @samples = sort grep { defined($_) && length($_) && !$seen{$_}++ }
			@{$ConspecificMGS{$MGS}};
		print LO $MGS . "\t" . join(",", @samples) . "\n";
	}
	close LO or die "Can't close conspecific log file: $conlog\n";
}

sub mergeConspecificLogs {
	return unless $maxSubJob;
	my %merged;
	for my $worker (0 .. $maxSubJob - 1) {
		my $part = "$LOGDIR/ConspecificMGS.$worker.log";
		die "Missing conspecific worker log: $part\n" unless -e $part;
		open my $in, '<', $part or die "Can't open conspecific worker log $part: $!\n";
		while (my $line = <$in>) {
			chomp $line;
			next unless length $line;
			my ($mgs, $sample_list) = split /\t/, $line, 2;
			die "Malformed conspecific worker log row in $part: $line\n"
				unless defined($mgs) && length($mgs) && defined($sample_list);
			$merged{$mgs}{$_} = 1 for grep { length } split /,/, $sample_list;
		}
		close $in or die "Can't close conspecific worker log $part: $!\n";
	}

	my $canonical = "$LOGDIR/ConspecificMGS.log";
	my $temporary = "$canonical.tmp.$$";
	open my $out, '>', $temporary or die "Can't write merged conspecific log $temporary: $!\n";
	for my $mgs (sort keys %merged) {
		print {$out} $mgs, "\t", join(',', sort keys %{$merged{$mgs}}), "\n"
			or die "Can't write merged conspecific log $temporary: $!\n";
	}
	close $out or die "Can't close merged conspecific log $temporary: $!\n";
	rename $temporary, $canonical
		or die "Can't install merged conspecific log $canonical: $!\n";
	%ConspecificMGS = map {
		$_ => [sort keys %{$merged{$_}}]
	} keys %merged;

	my @snp_parts = grep { -e $_ }
		map { "$outD/SNPconsCalls.$_.log" } 0 .. $maxSubJob - 1;
	if (@snp_parts) {
		my $snp_log = "$outD/SNPconsCalls.log";
		my $snp_tmp = "$snp_log.tmp.$$";
		open my $snp_out, '>', $snp_tmp or die "Can't write merged SNP consensus log $snp_tmp: $!\n";
		for my $part (@snp_parts) {
			open my $snp_in, '<', $part or die "Can't open SNP consensus worker log $part: $!\n";
			while (my $line = <$snp_in>) {
				print {$snp_out} $line or die "Can't write merged SNP consensus log $snp_tmp: $!\n";
			}
			close $snp_in or die "Can't close SNP consensus worker log $part: $!\n";
		}
		close $snp_out or die "Can't close merged SNP consensus log $snp_tmp: $!\n";
		rename $snp_tmp, $snp_log or die "Can't install merged SNP consensus log $snp_log: $!\n";
	}
}

sub writeTooFewMarker{
	my ($outD2, $sampleCount, $geneCount) = @_;
	make_path($outD2) unless -d $outD2;
	my $marker = "$outD2/tooFewSamples.sto";
	open my $out, '>', $marker or die "Cannot create $marker: $!\n";
	print {$out} "samples\t$sampleCount\ngenes\t$geneCount\n"
		or die "Cannot write $marker: $!\n";
	close $out or die "Cannot close $marker: $!\n";
}


sub prepGene2MGS{
	print "Preparing base strain alignments, per MGS\nThis might take a good while..\n";

	#If this run is split into subjobs, each worker only ever processes 1/maxSubJob of
	#the samples (see the identical stride logic later in extractFNAFAA2genes()). Previously
	#every worker still built the *complete* per-sample locus model (all samples, all MGS)
	#and only discarded the unneeded samples afterwards. Computing the worker's own sample
	#set up front lets us restrict the cluster-index parse itself, so the discarded data is
	#never materialized in this process at all.
	my $mySamplesHR = undef;
	if ($maxSubJob){
		my @srtdAllSmpls = sort @samples;
		my $Ndirs = scalar(@srtdAllSmpls);
		my %mine;
		for (my $i = $subJob; $i < $Ndirs; $i += $maxSubJob){
			$mine{$srtdAllSmpls[$i]} = 1;
		}
		$mySamplesHR = \%mine;
		print "Subjob ${subJob}/$maxSubJob: restricting locus-model construction to "
			. scalar(keys %mine) . " of $Ndirs samples\n";
	}

	my ($hr1,$cl2gene) = readClstrRev("$GCd/compl.incompl.$clusterID.fna.clstr.idx",0,$Gene2COG,$mySamplesHR);
	$hr1 = {};

	my $protein_file = "$GCd/compl.incompl.$clusterID.prot.faa";
	if (fileGZe($protein_file)) {
		$catalogProteins = readFasta($protein_file,1,"\\s",$Gene2COG);
	} else {
		warn "Catalogue protein file $protein_file is unavailable; keeping same-COG catalogue clusters separate\n";
	}

	my @records;
	for my $MGS (keys %{$COGprios}) {
		my $rank = 0;
		for my $seed_locus (@{$COGprios->{$MGS}}) {
			my $gene = $SIgenes->{$MGS}{$seed_locus};
			next unless defined $gene;
			push @records, {
				mgs => $MGS, cog => $Gene2COG->{$gene}, gene => $gene, rank => $rank++,
			};
		}
	}
	my $locus_model = build_locus_groups(
		\@records, $cl2gene, $catalogProteins,
		{
			# These indexes are useful to general callers but duplicate large
			# parts of the cluster model and are not consumed by this workflow.
			include_member_to_seed => 0,
			include_gene_to_locus => 0,
		},
	);
	my $ranked_record_count = scalar(@records);
	@records = ();
	$LocusByID = $locus_model->{locus_by_id};
	$MemberContext = $locus_model->{member_context};
	$LocusContext = $locus_model->{locus_context};

	my ($new_si_genes, $new_priorities) = ({}, {});
	for my $group (@{$locus_model->{groups}}) {
		$new_si_genes->{$group->{mgs}}{$group->{locus_id}} = $group->{primary_gene};
		push @{$new_priorities->{$group->{mgs}}}, $group->{locus_id};
	}
	$SIgenes = $new_si_genes;
	$COGprios = $new_priorities;

	my ($gene_sample_combinations, $ambiguous_seed_samples, $missing_clusters) = (0, 0, 0);
	my $unrepresentedWorkerLoci = 0;
	my (%contextMembersNeeded, %contextLociNeeded);
	for my $group (@{$locus_model->{groups}}) {
		my %per_sample;
		for my $seed (@{$group->{genes}}) {
			#delete (not just read) so the raw comma-joined membership string is freed the
			#moment it's consumed, rather than staying resident until a bulk clear at the
			#very end of this loop (which previously doubled peak memory: the fully-built
			#cl2gene2/candidateSeed structures existed alongside the still-intact $cl2gene).
			my $gene_string = delete $cl2gene->{$seed};
			unless (defined $gene_string) {
				if ($mySamplesHR) {
					# The cluster reader intentionally omits selected loci with
					# no member in this worker's sample partition.
					$unrepresentedWorkerLoci++;
					next;
				}
				limitedWarn('selected catalogue genes absent from cluster index',
					"Could not find selected catalogue gene $seed in the cluster index\n");
				$missing_clusters++;
				next;
			}
			for my $member (split /,/, $gene_string) {
				$member =~ s/^>//;
				next unless length $member;
				my ($sample) = split /__/, $member, 2;
				unless (defined($sample) && length($sample)) {
					limitedWarn('malformed catalogue cluster members',
						"Ignoring malformed catalogue member '$member' for seed $seed\n");
					next;
				}
				#belt-and-braces: readClstrRev already restricted members to this worker's
				#sample slice when $mySamplesHR was given, so this should normally be a no-op.
				next if $mySamplesHR && !exists $mySamplesHR->{$sample};
				$per_sample{$sample}{$member} = $seed;
			}
		}
		for my $sample (keys %per_sample) {
			#fold seed provenance directly into cl2gene2 (member => seed) instead of
			#keeping a fully parallel %candidateSeed hash-of-hash-of-hash with the same
			#member names duplicated again purely to carry the seed value.
			$cl2gene2{$sample}{$group->{locus_id}} = $per_sample{$sample};
			$smplsPerMGS{$group->{mgs}}{$sample}++;
			$gene_sample_combinations++;
			if (scalar(keys %{$per_sample{$sample}}) > 1) {
				$ambiguous_seed_samples++;
				$contextLociNeeded{$group->{locus_id}} = 1;
				$contextMembersNeeded{$_} = 1 for keys %{$per_sample{$sample}};
			}
		}
	}
	$cl2gene = {}; #any leftover (unconsumed) entries are dropped here
	# Context contributes only to multi-candidate resolution.  Unique candidates
	# bypass scoring, so retaining contexts for them only increases steady-state
	# extraction memory.
	my %keptMemberContext;
	for my $member (keys %contextMembersNeeded) {
		$keptMemberContext{$member} = $MemberContext->{$member}
			if exists $MemberContext->{$member};
	}
	$MemberContext = \%keptMemberContext;
	my %keptLocusContext;
	for my $locus (keys %contextLociNeeded) {
		$keptLocusContext{$locus} = $LocusContext->{$locus}
			if exists $LocusContext->{$locus};
	}
	$LocusContext = \%keptLocusContext;
	print "Prepared ".scalar(@{$locus_model->{groups}})." loci from $ranked_record_count"
		." ranked catalogue clusters; merged $locus_model->{merged_seeds} compatible same-COG seeds. "
		."$gene_sample_combinations locus-sample combinations, $ambiguous_seed_samples with multiple candidates"
		.($missing_clusters ? ", $missing_clusters missing cluster-index entries" : "")
		.($unrepresentedWorkerLoci ? ", $unrepresentedWorkerLoci loci outside this worker's sample slice" : "")
		.".\n";
}

sub prepRun{

	$mode = "FMG" if ($MGSfile eq "");
	if ($mode eq "FMG"){$takeAll = 0;}
	die "FMG mode does not support -maxSubJob; run it as a single extraction job\n"
		if $mode eq "FMG" && $maxSubJob;
	$takeAll = 1 if ($maxNGenes <= 0);
	if ($takeAll){$maxNGenes = -1;$mode="MGSall"; }


	$bindir = $MGSfile;$bindir =~ s/[^\/]+$//; 
	$bindir = $GCd if $bindir eq "";
	my $defaultOutD = $bindir."/intra_phylo/";
	$outD = $defaultOutD;#"$GCd/$mode/intra_phylo/";
	if ($outDpre ne ""){
		$outD = $outDpre ; 
		$outD .= "/" unless ($outD =~ m/\/$/);
		}
	my $safeDefaultOutD = $outDpre eq "" ? $defaultOutD : "";
	my $outputWasPresent = -d $outD ? 1 : 0;
	$LOGDIR = "$outD/LOGandSUB/";
	$SNPconsLOGs = "$outD/SNPconsCalls.$subJob.log" if ($SNPconsLOGs eq "");

	my $GCname = basename($GCd);
	my $outDname = basename($outD);
	die "Could not derive safe temporary-directory names\n" unless length($GCname) && length($outDname);
	$scratchD = getProgPaths("globalTmpDir",0);
	$scratchD = "$outD/.scratch" if $scratchD eq "";
	$scratchD .= "/strainsScr1/$GCname.$outDname/";
	#die "$scratchD  :$GCname :$GCd\n";
	if ($locTmpDir1 eq ""){
		my $locTmpN = getProgPaths("nodeTmpDir",0) ;
		my $suffix = ""; $suffix = "/SJ.${subJob}/" if ($subJob);
		if ($locTmpN eq ""){
			$locTmpDir =  "$outD/strainsScr1/$GCname.$outDname/$suffix" ; 
		} else {
			#my $tmp = `echo \$SLURM_LOCAL_SCRATCH`;
			#print STDERR "echo $locTmpN\n$tmp\n";
			#$locTmpN =~ s/\$/\\\$/;$locTmpN = `echo $locTmpN;`; #eval in sys
			$locTmpN=truePath($locTmpN,1);
			#die $locTmpN."\n";
			$locTmpDir =  "$locTmpN/strainsScr1/$GCname.$outDname/$suffix" ; 
		}
	} else {
		my $suffix = $subJob ? "/SJ.${subJob}" : "";
		$locTmpDir = "$locTmpDir1/strainsScr1/$GCname.$outDname$suffix/";
	}
	
	$preConDir = "$scratchD/preComp/";


	print "\n!! WARNING !!: RESUBMISSION mode selected (will resubmit MSA + phylos even for already completed MGS) !!\n" if ($reSubmit);
	print "\n!! WARNING !!: REDOSUBMISSIONDATA mode selected (will redo and resubmit MSA + phylos even for already completed MGS) !!\n" if ($redoSubmissionData);

	open my $map_info, '<', "$GCd/LOGandSUB/GCmaps.inf"
		or die "Cannot open $GCd/LOGandSUB/GCmaps.inf: $!\n";
	$mapF = <$map_info> // "";
	close $map_info or die "Cannot close map information file: $!\n";
	chomp $mapF;
	die "Mapping-file reference is empty or missing: $mapF\n" if ($mapF !~ /,/ && (!length($mapF) || !-e $mapF));
	
	#read info gene <-> taxonomy from this file, depends on config..
	$gene2taxF = "$GCd/FMG/gene2specI.txt";
	$gene2taxF = "$GCd/GTDBmg/gene2specI.txt" if ($useGTDBmg eq "GTDB");
	#die;

	#---------------
	#everything after is only for main submission job..
	if ($subJob){
		print "=============\n=============\nStrain_within v$version, subjob ${subJob}/$maxSubJob\n=============\n=============\n";
	} else {
			print "============= Strain_within v$version =============\n";
		print "Creating within species strains for ${mode}s in $GCd\n";
		print "Outdir: $outD\nTmpDir: $locTmpDir\nScratchDir: $scratchD\n";
		print "GC dir: $GCd\nIn Cluster: $MGSfile\nCores: $numCores (max: ${maxCores})\n";
		print "MAP: $mapF\n";
		#print "Ref tree: $treeFile\n";
		print "Using tree $treeFile to create automatically outgroups\n" if ($treeFile ne "");
		print "MGs: $useGTDBmg\nGene2Tax: $gene2taxF\n";
		print "Using $presortGenes genes from each MGS for location\n";
		print "Deep repariing remaining submission files\n" if ($deepRepair);
		print "Pre-creating ConsSNPs in $preConDir in $preCompCons runs\n" if ($preCompCons);
		print "-minSNPDepth $minSNPDepth, -minSNPCallQual $minSNPCallQual";
		print ", -SNPadaptiveQual $useAdaptiveQual, -SNPindelRangeFilt: $indelRange";
		if ($depthFilterScale){print ", depthFiltScale $depthFilterScale\n";}else {print "\n";}
		print "DiscTests=$discTests\n" unless ($discTests eq "");
		print "ContTests=$contTests\n" unless ($contTests eq "");
		print "familyVar=$familyVar\n" unless ($familyVar eq "");
		
		print "groupStabilityVars=$groupStabilityVars\n" unless ($groupStabilityVars eq "");
		print "MSAaligner: $MSAprog, GenesPerSpecies: $GenesPerSpecies, GeneLengthMin: $GeneLengthMin\n";
		
		
		if ($takeAll){print "**************** Take all genes MGS mode\n";}
		else {print "Using first $maxNGenes genes found per sample\n";}
		print "==============================================\n";
		if ($onlySubmit){print "Only submission mode\n";
		} elsif (!$subJob) {
			print "Creation of strain genes, old data might be deleted!\nDo you want to continue? (10s wait, use Ctrl-c to abort)\n"; sleep 10;
		}
	}
	
	

	#$mapF = $GCd."LOGandSUB/inmap.txt" if ($mapF eq "");
	my ($hr1,$hr2) = readMapS($mapF,-1);
	%map = %{$hr1}; %AsGrps = %{$hr2};
	#get all samples in assembly group, but only last in mapgroup
	@samples = @{$map{opt}{smpl_order}};
	my %sample_seen;
	for my $sample (@samples) {
		die "Unsafe sample identifier '$sample': use only letters, digits, dot, underscore, colon, plus, and hyphen\n"
			unless defined($sample) && $sample =~ /\A[A-Za-z0-9][A-Za-z0-9_.:+-]*\z/;
		die "Duplicate sample identifier in map: $sample\n" if $sample_seen{$sample}++;
	}


	if ($mode eq "MGS" || $mode eq "MGSall"){
		my $sortedMGS = "$MGSfile.srt";
		if ($subJob) {
			die "Sorted MGS guide is missing for subjob: $sortedMGS\n" unless -s $sortedMGS;
		} elsif ($mode eq "MGSall" && !-e $sortedMGS) {
			assertSafeWorkflowRemoval($outD, $safeDefaultOutD, $GCd, $MGSfileOri, $bindir, getcwd()) if -d $outD;
			remove_tree($outD) if -d $outD;
			remove_tree($scratchD) if -d $scratchD;
			unlink $_ or die "Cannot remove stale $_: $!\n"
				for grep { -f $_ || -l $_ } glob("$MGSfile.srt*");
			symlink($MGSfile, $sortedMGS)
				or die "Cannot link $sortedMGS to $MGSfile: $!\n";
		} elsif (!$onlySubmit || !-s $sortedMGS) {
			print "base files missing.. preparing complete resubmission and recalc of data\n";
			assertSafeWorkflowRemoval($outD, $safeDefaultOutD, $GCd, $MGSfileOri, $bindir, getcwd()) if -d $outD;
			remove_tree($outD) if -d $outD;
			remove_tree($scratchD) if -d $scratchD;
			unlink $_ or die "Cannot remove stale $_: $!\n"
				for grep { -f $_ || -l $_ } glob("$MGSfile.srt*");
			my $sortMGSgenes = getProgPaths("sortMGSGeneImport_scr");
			my $cmd = $sortMGSgenes . " "
				. join(" ", map { shellQuote($_) } ($GCd, $MGSfile, $useGTDBmg, $mode, $clusterID)) . "\n";
			print "$cmd\n";
			systemW $cmd;
			die "MGS sorting did not create $sortedMGS\n" unless -s $sortedMGS;
		} else {
			print "Continuing on prepared .srt files\n";
		}

		$MGSfile = $sortedMGS;
		$gene2taxF = createGene2MGS($MGSfile,$GCd);
		print "Using sorted MGS from $MGSfile, adding eggNOG in: $gene2taxF\n";
		print "\nnew MGS file: $MGSfile\n\n";
	} elsif ($subJob && $maxSubJob) {
		die "FMG mode does not support split MGS extraction jobs\n";
	}

	if ($subJob){
		return;
	}

	make_path($locTmpDir, $scratchD, $outD, $LOGDIR);
	my $outputBase = basename(File::Spec->canonpath($outD));
	my $owner = File::Spec->catfile($outD, '.matafiler-strain-workdir');
	markStrainWorkflowDirectory($outD)
		if !$outputWasPresent || !$onlySubmit || $outputBase eq 'intra_phylo'
			|| $outputBase eq 'within_phylo' || -e $owner;
	open FO, ">$LOGDIR/strainCmd.txt" or die "Cannot write $LOGDIR/strainCmd.txt: $!\n";
	print FO $cmdCall;
	close FO or die "Cannot close $LOGDIR/strainCmd.txt: $!\n";
	
	#DEBUG
	if ($preCompCons && !$subJob) {
		remove_tree($preConDir) if -d $preConDir;
		make_path($preConDir);
	}

	#STONES
	make_path("$outD/stones/");
	
	make_path($locTmpDir);
	open my $tmp_test, '>', "$locTmpDir/test.txt" or die "Couldn't create test file in local dir $locTmpDir: $!\n";
	close $tmp_test or die "Couldn't close test file in local dir $locTmpDir: $!\n";
	if ( ! -e "$locTmpDir/test.txt"){die "Couldn't create test file in local dir $locTmpDir\n";}
#die "passed $locTmpDir\n";

	return;
}


sub preComputeConsSNP{
	my $inputChk = "$outD/stones/0.fileChk.sto";
	my $fileAbsent = 0;
	my @missing_samples;
	my $submPreComp = 1;#DEBUG
	$submPreComp = 0 if ($subJob);

	
	my @accumVCFcmds; my $BatchCnt=0;my @jobsPre;
	foreach my $smpl (@samples){ # just check that files are there..
		# Always revalidate paired consensus files; an old checkpoint cannot prove
		# that both the nucleotide and protein cache still exist.
		unless (exists($map{$smpl}) && defined($map{$smpl}{wrdir}) && length($map{$smpl}{wrdir})) {
			limitedWarn('samples without working directories',
				"No working directory is configured for $smpl; sample will be skipped\n");
			$fileAbsent = 1;
			$unavailableSamples{$smpl} = "missing map working directory";
			push @missing_samples, $smpl;
			next;
		}
		my $cD = $map{$smpl}{wrdir}."/";
		if (-e "$cD/SMPL.empty") {
			$unavailableSamples{$smpl} = "sample is marked empty";
			next;
		}
		#my $tarF = $cD."/SNP/genes.shrtHD.SNPc.MPI.fna.gz";
		my $tarF = $cD."/$lSNPdir/$lConsFNA";
		my $tarF2 = $cD."/$lSNPdir/$lConsFAA";
		my $tarVCF = $cD."/$lSNPdir/$lConsVCF";
		my $input_state = consensusInputState(
			fileGZe($tarVCF), fileGZe($tarF), fileGZe($tarF2), $forceVCF2FNA
		);
		if ($input_state eq 'missing') {
			limitedWarn('samples without usable consensus inputs',
				"Can't find a complete consensus pair or a VCF to repair it for $smpl in $cD; sample will be skipped\n");
			$fileAbsent = 1;
			$unavailableSamples{$smpl} = "missing consensus pair and repair VCF";
			push @missing_samples, $smpl;
			next;
		}
		
		if ($preCompCons && $input_state eq 'regenerate'){
			#store these in scratch, uncompressed (much faster)
			my $fastaf = "$preConDir/$smpl.cons.genes.fna.gz";
			my $fastafAA = "$preConDir/$smpl.cons.prots.faa.gz";
			my $vcf2fnaCmd = createConsFastas($cD, $smpl, $fastaf, $fastafAA, 0, 1);
			$preCompSNPs{$smpl}{NT}=$fastaf;$preCompSNPs{$smpl}{AA}=$fastafAA;

			push(@accumVCFcmds,$vcf2fnaCmd);
			
			if (@accumVCFcmds >= $preCompCons){
				print "Precomp batch $BatchCnt " if ($submPreComp);
				my $cmdX = "\necho \"BATCH $BatchCnt\"\nmkdir -p ".shellQuote($preConDir).";\n\n" . join("\n",@accumVCFcmds);
				my $tmpSHDD=$QSBoptHR->{tmpSpace} ; $QSBoptHR->{tmpSpace} =0;
				my ($dep,$qcmd) = qsubSystem($LOGDIR."PreCompConsSNP_B${BatchCnt}.sh",$cmdX,1,"10G","ConsSNP$BatchCnt","","",$submPreComp,[],$QSBoptHR);
				$QSBoptHR->{tmpSpace} =$tmpSHDD;

				push(@jobsPre,$dep) if defined($dep) && length($dep);
				#reset counters etc
				$BatchCnt++; @accumVCFcmds=();

				#die;
			}
		}
	}
	#last batch of jobs..
	if (@accumVCFcmds){
		
		my $cmdX = "\necho \"BATCH $BatchCnt\"\nmkdir -p ".shellQuote($preConDir).";\n\n" . join("\n",@accumVCFcmds);
		my ($dep,$qcmd) = qsubSystem($LOGDIR."PreCompConsSNP_B${BatchCnt}.sh",$cmdX,1,"10G","ConsSNP$BatchCnt","","",$submPreComp,[],$QSBoptHR);
		push(@jobsPre,$dep) if defined($dep) && length($dep);

	}
	if (@jobsPre && $doSubmit){
		qsubSystemJobAlive( \@jobsPre,$QSBoptHR );
	}
	for my $smpl (keys %preCompSNPs) {
		my $nt = $preCompSNPs{$smpl}{NT};
		my $aa = $preCompSNPs{$smpl}{AA};
		next if fileGZe($nt) && fileGZe($aa);
		limitedWarn('incomplete precomputed consensus outputs',
			"Precomputed consensus output is incomplete for $smpl; falling back to on-the-fly generation\n");
		delete $preCompSNPs{$smpl};
	}
	if ($fileAbsent) {
		warn scalar(@missing_samples)." samples lack required SNP inputs and will be skipped: "
			.join(",", @missing_samples[0 .. ($#missing_samples < 9 ? $#missing_samples : 9)])
			.(@missing_samples > 10 ? ",..." : "")."\n";
		if (-e $inputChk) {
			unlink $inputChk or warn "Cannot remove stale input checkpoint $inputChk: $!\n";
		}
	}
	unless ($fileAbsent || -e "$inputChk"){
		print "All samples have SNP calls\n";
		open my $checkpoint, '>', $inputChk or die "Cannot create $inputChk: $!\n";
		close $checkpoint or die "Cannot close $inputChk: $!\n";
	}
}


sub createAGlist{
	foreach my $smpl (@samples){ #fill up AGlist
	#and fill %AGlist .. so always let run..
		next if ($map{$smpl}{AssGroup} eq "-1");
		my $cAssGrp = $map{$smpl}{AssGroup};
		my $cMapGrp = $map{$smpl}{MapGroup};
		die "Can't find mapping-group counters for $cMapGrp\n"
			unless exists($AsGrps{$cMapGrp}) && exists($AsGrps{$cMapGrp}{CntAimMap});
		#print "$smpl $cAssGrp $cMapGrp\n";
		$AsGrps{$cMapGrp}{CntMap} = 0 unless exists $AsGrps{$cMapGrp}{CntMap};
		$AsGrps{$cMapGrp}{CntMap} ++;
		next if ($AsGrps{$cMapGrp}{CntMap}  < $AsGrps{$cMapGrp}{CntAimMap} );
		push(@{$AGlist{$cAssGrp}} , $smpl);
		#if ($AsGrps{$cMapGrp}{CntMap}  < $AsGrps{$cMapGrp}{CntAimMap} ){			next;		}
	}
}

sub histoMGS{#specifically for MGS..
	my ($aref,$msg) = @_;
	my @cnts = @{$aref};
	my @binSiz = (10,20,30,50,70,100,200,300,500,700,1000,2000,5000,10000,1e6);
	my %binC; #my $prevC=0;
	foreach (@binSiz){$binC{$_} = 0;}
	foreach my $c(@cnts){
		#print "$c ";
		my $bs=0;
		while ($bs < @binSiz - 1 && $c > $binSiz[$bs]) {
			$bs++;
		}
		#print " X$c:${bs}X ";
		$binC{$binSiz[$bs]} ++;
	}
	#display bin counts..
	print $msg.": ";#"Bin size distribution: ";
	foreach (@binSiz){print " <$_:$binC{$_} " if ($binC{$_}>0);}
	print "\n";
	#DEBUG
	#print @cnts." : @cnts\n";
}

sub getInputSize{
	# fileGZs($FNAtf) / (1024 * 1024)
	my @out; my @missedMGS;
	foreach my $MGS (@specis){
		my $tmpD  = "$scratchD/outs/$MGS/";
		my $FNAtf = "$tmpD/$FNAstdof";
		my $inputFNAsize=0;
		my @multiM = glob("$FNAtf*");
		if (@multiM){
			foreach(@multiM){
			$inputFNAsize += fileGZs($_) / (1024 * 1024) ;
			}
		} elsif (fileGZe( "$SIdirs{$MGS}/$FNAstdof")) {
			#$FNAtf = "$SIdirs{$MGS}/$FNAstdof";
			#fileGZs
			if (-e "$SIdirs{$MGS}/$FNAstdof"){
				$inputFNAsize = fileGZs("$SIdirs{$MGS}/$FNAstdof") / (1024 * 1024) ;
			} elsif (-e "$SIdirs{$MGS}/$FNAstdof.gz"){
				$inputFNAsize = fileGZs("$SIdirs{$MGS}/$FNAstdof.gz") / (1024 * 1024)*40 ;
			}
		} else {
			push(@missedMGS,$MGS);
			$inputFNAsize = 0;
		}
		push(@out, $inputFNAsize); 
	}
	if (@missedMGS){print "\ngetInputSize:: could not find FNA for: @missedMGS\n";}
	#print "SIZE: @out\n";
	#die;
	return @out;
}


sub evalFileStatus{
	my $dirsNOTPrepped = 0; my $CatFileMiss = 0;my $CatNotPrepped = 0; my $treeAbsent=0;
	my $doneDirs=0;
	my $tooFewDirs=0;
	my $PhylosExist = 1;
	
	my $treeFile= "IQtree_allsites.treefile";
	if ($phyloProg == 2){$treeFile = "VERYFASTTREE_allsites.nwk";} elsif ($phyloProg == 3){$treeFile = "FASTTREE_allsites.nwk";}


	foreach my $MGS (@specis){ #loop creates per specI file structure to run buildTreeScript on..
		#PART I: create fasta files required by tree
		my $outD2 = "$outD/$MGS/";
		$SIdirs{$MGS} = $outD2;
		#print "$outD2\n";
		if (-d $outD2 && $onlySubmit == 0 && !$subJob){#only the parent may clean shared folders
			remove_tree($outD2);
			my $scratch_mgs = "$scratchD/outs/$MGS";
			remove_tree($scratch_mgs) if -d $scratch_mgs;
		}
		make_path($outD2) unless -d $outD2;
		my $tooFewMarker = "$outD2/tooFewSamples.sto";
		if (-s $tooFewMarker && !$deepRepair && !$redoSubmissionData && $onlySubmit != 0) {
			$MGSnoTree{$MGS} = 1;
			$tooFewDirs++;
			next;
		}
		unlink $tooFewMarker or die "Cannot remove stale $tooFewMarker: $!\n"
			if -e $tooFewMarker;
		
	#	if ( !-d $outD2 ||){ # first phase only has "all.cat.tmp" file..
	#		$dirsNOTPrepped ++;
	#	} els
		
		if (fileGZe("$SIdirs{$MGS}/$CATstdof")) {
			my ($category_fh) = gzipopen("$SIdirs{$MGS}/$CATstdof", "existing locus category file", 0);
			my $first_entry = '';
			if ($category_fh) {
				while (my $line = <$category_fh>) {
					chomp $line;
					next unless length $line;
					($first_entry) = split /\t/, $line, 2;
					last;
				}
				close $category_fh;
			}
			my @identifier_parts = split /\Q$SaSe\E/, $first_entry, -1;
			if (@identifier_parts != 3 || grep { !length } @identifier_parts) {
				limitedWarn('MGS with legacy sequence identifiers',
					"$MGS does not use the required sample|COG|primaryGeneID identifier format; scheduling input regeneration\n");
				$legacyLocusOutputs++;
				$legacyLocusMGS{$MGS} = 1;
				$dirsNOTPrepped++;
				$CatFileMiss++;
				next;
			}
		}
		if (!fileGZe("$SIdirs{$MGS}/$CATstdof")){
			$CatFileMiss ++ ; 
			my @multiM=glob("$scratchD/outs/$MGS/all.cat.tmp*");
			if ( -e "$outD2/all.cat.tmp" || @multiM ){
				$CatNotPrepped ++; #tmp file exists, needs Part II
			} else {
				#print "$scratchD/outs/$MGS\n";
				$dirsNOTPrepped ++; #needs Part I
			}
			#print "$SIdirs{$MGS}\n";
			#system "rm $SIdirs{$MGS}\n";
		}elsif(!fileGZs("$SIdirs{$MGS}/phylo/$treeFile")){
			$treeAbsent++;
			remove_tree("$scratchD/outs/$MGS") if -d "$scratchD/outs/$MGS";
		} elsif(fileGZe("$SIdirs{$MGS}/phylo/$treeFile")) {
			$doneDirs++;
			remove_tree("$scratchD/outs/$MGS") if -d "$scratchD/outs/$MGS";
		}
	}
	$PhylosExist = 0 if ($CatFileMiss/scalar(@specis) > 0.1); #only activate if more than 10% missing..

	print "Output dirs status: \nCatFileFinalMiss: $CatFileMiss, CatFileConvert: $CatNotPrepped, Dir not done: $dirsNOTPrepped, phylo absent: $treeAbsent, Dir done: $doneDirs, too few samples: $tooFewDirs, Phylo complete: $PhylosExist \n";
	#die;
	return($dirsNOTPrepped , $CatFileMiss , $CatNotPrepped , $treeAbsent, $doneDirs, $PhylosExist);
}


sub appendWriteMGSgenes {
    my ($writeLink) = @_;
	print "Flushing buffered MGS records\n";

    my $wrMGS = 0;
    my $suffix = ".$subJob";
    my $baseOut = "$scratchD/outs";

    foreach my $MGS (keys %OFstrH) {

        my $nt = $OFstrH{$MGS} or next;
		next if ($nt eq "");

        my $aa   = $OAstrH{$MGS};
        my $cat  = $OCstrH{$MGS};
        my $link = $OLstrH{$MGS};

        my $outD = "$baseOut/$MGS";
        make_path($outD) unless -d $outD;

        my $FNAtf = "$outD/$FNAstdof$suffix";
        my $FAAtf = "$outD/$FAAstdof$suffix";
        my $CATtf = "$outD/$CATstdof.tmp$suffix";

		open my $fh_nt, ">>", $FNAtf or die $!;
		print {$fh_nt} $nt or die "Cannot append $FNAtf: $!\n";
		close $fh_nt or die "Cannot close $FNAtf: $!\n";

		open my $fh_aa, ">>", $FAAtf or die $!;
		print {$fh_aa} $aa or die "Cannot append $FAAtf: $!\n";
		close $fh_aa or die "Cannot close $FAAtf: $!\n";

        if ($writeLink) {
			my $Linkf = "$outD/$LINKstdof$suffix";
			open my $fh_link, ">>", $Linkf or die $!;
			print {$fh_link} $link or die "Cannot append $Linkf: $!\n";
			close $fh_link or die "Cannot close $Linkf: $!\n";
		}

		open my $fh_cat, ">>", $CATtf or die $!;
		print {$fh_cat} $cat or die "Cannot append $CATtf: $!\n";
		close $fh_cat or die "Cannot close $CATtf: $!\n";

        $OFstrH{$MGS} = "";
        $OAstrH{$MGS} = "";
        $OCstrH{$MGS} = "";
        $OLstrH{$MGS} = "";

        $wrMGS++;
    }

    print "\nwrote for $wrMGS MGS data..\n";
}


sub appendWriteMGSgene_olds{
	#write genes to respective MGS intra phyla..
	my ($writeLink) = @_;
	print "Flushing buffered MGS records\n";
	my $wrMGS=0;
	my @SpecSet = keys(%OFstrH);
	my @specSetS = shuffle(@SpecSet); #shuffle to further reduce chance of multiple jobs writing consistently to the same files..
	my $FileSuff = ""; 
	$FileSuff = ".$subJob";# if ($subJob);
	foreach my $MGS (@specSetS){
		next if ($OFstrH{$MGS} eq "");#(!exists($OFstrH{$MGS}) || scalar(@{$OFstrH{$MGS}}) == 0 );
		my $hasSlept=0;
		#handle file paths..
		my $outD2 = $SIdirs{$MGS};
		$outD2 = "$scratchD/outs/$MGS/";
		system "mkdir -p $outD2" unless (-d $outD2);
		my $FNAtf = "$outD2/$FNAstdof$FileSuff"; my $FAAtf = "$outD2/$FAAstdof$FileSuff";my $Linkf = "$outD2/$LINKstdof$FileSuff";
		my $CATtf = "$outD2/$CATstdof.tmp$FileSuff";
		my $blockF = "$outD2/block.tmp";
		
		#block sample for other writes..
		#while (-e $blockF){sleep(5);$hasSlept=1;}while ($hasSlept && -e $blockF){sleep(8);}#second security layer..
		#system "touch $blockF";sleep(4) if ($hasSlept);#security that other process has finished writes completely
		#deactivate, go for unique file instead..

		#writing strings out..
		open OF,">>$FNAtf" or die "Can't append NT file $FNAtf\n";print OF $OFstrH{$MGS}; close OF;
		open OA,">>$FAAtf" or die "Can't append AA file $FAAtf\n";print OA $OAstrH{$MGS}; close OA;
		if ($writeLink){open OL,">>$Linkf" or die "Can't append link file $Linkf\n" ; print OL $OLstrH{$MGS}; close OL;}
		#this is only a temp file, that needs to be rewritten later..
		open OC,">>$CATtf" or die "Can't append to CAT file $CATtf\n";print OC $OCstrH{$MGS} ; close OC;
		
		#open OF,">>$FNAtf" or die "Can't append NT file $FNAtf\n";foreach(@{$OFstrH{$MGS}}){print OF $_;} close OF;
#		open OA,">>$FAAtf" or die "Can't append AA file $FAAtf\n";print OA join("",@{$OAstrH{$MGS}}); close OA;
		#open OA,">>$FAAtf" or die "Can't append AA file $FAAtf\n";foreach(@{$OAstrH{$MGS}}){print OA $_;} close OA;#print OA join("",@{$OAstrH{$MGS}}); close OA;
		#if ($writeLink){open OL,">>$Linkf" or die "Can't append link file $Linkf\n" ; foreach(@{$OLstrH{$MGS}}){print OL $_;}  close OL;}
		#this is only a temp file, that needs to be rewritten later..
		#open OC,">>$CATtf" or die "Can't append to CAT file $CATtf\n"; foreach(@{$OCstrH{$MGS}}){print OC $_;}   close OC;
		
		
		#$OCstrH{$MGS} = []; $OFstrH{$MGS} = []; $OAstrH{$MGS} = []; $OLstrH{$MGS} = [];
		$OCstrH{$MGS} = ""; $OFstrH{$MGS} = ""; $OAstrH{$MGS} = ""; $OLstrH{$MGS} = "";
		$wrMGS++;
		
		system "rm -f $blockF";
	}
	print "\nwrote for $wrMGS MGS data..\n";
}


sub reportingsMGS{
	#eval #sample/MGS
	my %smplPmgs;
	foreach my $MGS (keys %smplsPerMGS){
		foreach my $sm (keys %{$smplsPerMGS{$MGS}}){
			$smplPmgs{$MGS}++ if ($smplsPerMGS{$MGS}{$sm} > 10);
		}
	}
	my @smplNs = values(%smplPmgs);@smplNs = sort { $b <=> $a}  @smplNs;
	if (@smplNs) {
		my $qt50=quantileArray(0.5,@smplNs);my $qt90=quantileArray(0.90,@smplNs);
		my @top = @smplNs[0 .. ($#smplNs < 4 ? $#smplNs : 4)];
		print "Samples/MGS: QTL 50,90: $qt50 $qt90 . Top 5: ".join(" ",@top)."\n";
	} else {
		print "Samples/MGS: none have more than 10 candidate loci\n";
	}
	#die;
	
}

sub timeNice($){
	my ($tIN) = @_;
	$tIN = int($tIN);
	if ($tIN > (3600)){
		my $remMin = ($tIN%3600);
		return int($tIN/3600)."h".int($remMin/60)."m" . ($remMin%60) . "s";
	}
	if ($tIN > 60){
		return int($tIN/60)."m" . ($tIN%60) . "s";
	}
	return $tIN . "s";
}

sub shellQuote {
	my ($value) = @_;
	$value = "" unless defined $value;
	$value =~ s/'/'"'"'/g;
	return "'$value'";
}

sub treeInputPrecopyCommand {
	my ($staging_dir, $output_dir, $cores, @required_inputs) = @_;
	die "Cannot create a tree-input publication command without required inputs\n"
		unless @required_inputs;

	my $staging_q = shellQuote($staging_dir);
	my $output_q = shellQuote($output_dir);
	my $pigz_q = shellQuote($pigzBin);
	my $ready_test = join(" && ", map {
		"(test -s ".shellQuote($_)." || test -s ".shellQuote("$_.gz").")"
	} @required_inputs);
	my $required_array = join(" ", map { shellQuote($_) } @required_inputs);

	my $command = "\n# Persistent inputs take precedence during recovery; use staging only when needed.\n";
	$command .= "if ( $ready_test ); then\n";
	$command .= "  echo \"Using existing persistent tree inputs\"\n";
	$command .= "else\n";
	$command .= "  staged_inputs=()\n";
	$command .= "  if [[ -d $staging_q ]]; then\n";
	$command .= "    mapfile -d '' -t staged_inputs < <(find $staging_q -mindepth 1 -maxdepth 1 -type f -print0)\n";
	$command .= "  fi\n";
	$command .= '  if (( ${#staged_inputs[@]} )); then' . "\n";
	$command .= "    echo \"Publishing staged tree inputs to persistent storage\"\n";
	$command .= '    for staged in "${staged_inputs[@]}"; do' . "\n";
	$command .= '      if [[ "$staged" != *.gz ]]; then' . "\n";
	$command .= "        $pigz_q -p $cores -- " . '"$staged"' . "\n";
	$command .= "      fi\n";
	$command .= "    done\n";
	$command .= "    mapfile -d '' -t staged_inputs < <(find $staging_q -mindepth 1 -maxdepth 1 -type f -print0)\n";
	$command .= '    mv -- "${staged_inputs[@]}" ' . "$output_q\n";
	$command .= "  else\n";
	$command .= "    echo \"No usable staged tree inputs found\"\n";
	$command .= "  fi\n";
	$command .= "fi\n";
	$command .= "if ! ( $ready_test ); then\n";
	$command .= "  echo \"ERROR: tree inputs are incomplete in both staging and persistent storage\" >&2\n";
	$command .= "  required_inputs=($required_array)\n";
	$command .= '  for required in "${required_inputs[@]}"; do' . "\n";
	$command .= '    if [[ ! -s "$required" && ! -s "$required.gz" ]]; then' . "\n";
	$command .= q{      printf '  missing: %s[.gz]\n' "$required" >&2} . "\n";
	$command .= "    fi\n";
	$command .= "  done\n";
	$command .= "  exit 66\n";
	$command .= "fi\n";
	$command .= "echo \"Tree inputs ready in persistent storage\"\n\n";
	return $command;
}

sub readFastaIDs {
	my ($path) = @_;
	my %ids;
	return \%ids unless defined($path) && length($path) && fileGZe($path);
	my ($fh) = gzipopen($path, "FASTA identifier scan", 0);
	while (my $line = <$fh>) {
		$ids{$1} = 1 if $line =~ /^>(\S+)/;
	}
	close $fh or die "Cannot close FASTA identifier input $path: $!\n";
	return \%ids;
}

sub consensusInputState {
	my ($vcf_ready, $nt_ready, $aa_ready, $force_regeneration) = @_;
	return 'missing' if $force_regeneration && !$vcf_ready;
	return 'ready' if !$force_regeneration && $nt_ready && $aa_ready;
	return 'regenerate' if $vcf_ready;
	return 'missing';
}

sub assertSafeWorkflowRemoval {
	my ($target, $default_target, @protected) = @_;
	return unless -d $target;
	my $resolved = abs_path($target)
		or die "Cannot resolve workflow output directory before removal: $target\n";
	$resolved = File::Spec->canonpath($resolved);
	my ($volume) = File::Spec->splitpath($resolved, 1);
	my $root = File::Spec->canonpath(File::Spec->catpath($volume, File::Spec->rootdir(), ''));
	my $compare_target = $^O eq 'MSWin32' ? lc($resolved) : $resolved;
	my $compare_root = $^O eq 'MSWin32' ? lc($root) : $root;
	die "Refusing to remove filesystem root as a strain workflow directory: $resolved\n"
		if $compare_target eq $compare_root;

	my $prefix = $compare_target;
	$prefix .= File::Spec->catfile('', '') unless $prefix =~ m{[\\/]$};
	for my $protected (@protected) {
		next unless defined($protected) && length($protected) && -e $protected;
		my $resolved_protected = abs_path($protected) or next;
		$resolved_protected = File::Spec->canonpath($resolved_protected);
		my $compare_protected = $^O eq 'MSWin32' ? lc($resolved_protected) : $resolved_protected;
		die "Refusing to remove $resolved because it contains protected path $resolved_protected\n"
			if $compare_protected eq $compare_target || index($compare_protected, $prefix) == 0;
	}

	my $owner = File::Spec->catfile($resolved, '.matafiler-strain-workdir');
	my $is_default = 0;
	if (defined($default_target) && length($default_target) && -d $default_target) {
		my $resolved_default = abs_path($default_target);
		if (defined $resolved_default) {
			$resolved_default = File::Spec->canonpath($resolved_default);
			$resolved_default = lc($resolved_default) if $^O eq 'MSWin32';
			$is_default = 1 if $resolved_default eq $compare_target;
		}
	}
	die "Refusing to remove unowned custom output directory $resolved; expected $owner\n"
		unless $is_default || -f $owner;
}

sub markStrainWorkflowDirectory {
	my ($target) = @_;
	make_path($target) unless -d $target;
	my $owner = File::Spec->catfile($target, '.matafiler-strain-workdir');
	return if -e $owner;
	open my $fh, '>', $owner or die "Cannot create strain workflow ownership marker $owner: $!\n";
	print {$fh} "strain_within\t$version\n"
		or die "Cannot write strain workflow ownership marker $owner: $!\n";
	close $fh or die "Cannot close strain workflow ownership marker $owner: $!\n";
}




#this routine hast to get genes out of each sample, that are needed
#and save them to be later written per specI
sub extractFNAFAA2genes{
	# Each worker owns one numeric suffix.  A retry must replace, not append to,
	# that worker's previous partial extraction.
	for my $pattern (
		"$scratchD/outs/*/$FNAstdof.$subJob",
		"$scratchD/outs/*/$FAAstdof.$subJob",
		"$scratchD/outs/*/$LINKstdof.$subJob",
		"$scratchD/outs/*/$CATstdof.tmp.$subJob",
	) {
		for my $part (bsd_glob($pattern)) {
			unlink $part or die "Cannot remove stale worker part $part: $!\n"
				if -f $part || -l $part;
		}
	}
	my %perMGScnts;
	my %representedLocus;
	my $gnCnt=0;
	#my %totGnes;
	#create gene to genes list
	foreach my $sm (keys %cl2gene2){
		#my @locGenes;
		#print "$sm ";
		my $MGSgeneCnt=0;
		foreach my $gn (keys %{$cl2gene2{$sm}}){
			#$totGnes{$gn} = 1;
			$gnCnt++;
			if (exists($LocusByID->{$gn}) && !exists($representedLocus{$gn})){
				$representedLocus{$gn} = 1;
				$perMGScnts{$LocusByID->{$gn}{mgs}}++;
				#print "1";
				$MGSgeneCnt++;
			}
		}
		#print "$sm  $gnCnt $MGSgeneCnt \n";
	}
	my @histoMGScnts ;#= values %perMGScnts;
	my $lowCandidateMGS = 0;
	foreach my $MGS (keys %perMGScnts){
		my $perMGSgenes = $perMGScnts{$MGS};
		push(@histoMGScnts,  $perMGSgenes);
		if ($perMGSgenes < 10){
			$lowCandidateMGS++;
			limitedWarn('MGS with fewer than 10 candidate loci',
				"Only $perMGSgenes genes/COGs for MGS $MGS; MGS genes might be multi-copy\n")
				unless $maxSubJob;
		}
	}
	print "$lowCandidateMGS MGS have fewer than 10 candidate loci in this worker's sample slice; "
		."this is expected for sparse split-worker partitions\n"
		if $maxSubJob && $lowCandidateMGS;
	#DBUG
	my $represented_mgs = scalar(keys(%perMGScnts));
	my $average_loci = $represented_mgs ? int(0.5 + $gnCnt / $represented_mgs) : 0;
	print "Loci per MGS (prefiltering, N= ". $gnCnt  ." loci, $represented_mgs MGS, avg $average_loci loci/MGS):\n";
	histoMGS(\@histoMGScnts,"Theorectical best Bin sizes: ");
	#some stats on genes/MGS
	my @srtdSmpls = sort (keys %cl2gene2);
	
	#subjob? samples are already restricted to this worker's slice: prepGene2MGS() now
	#builds %cl2gene2 directly from a pre-filtered cluster-index parse (see the
	#$mySamplesHR restriction there), so @srtdSmpls (= sort keys %cl2gene2) is already
	#exactly this worker's share. Re-applying the same stride split here on the already-
	#reduced key set would incorrectly select only every Nth *remaining* sample and
	#silently drop the rest, so we no longer do that -- just report what we got.
	if ($maxSubJob){
		my $Ndirs = scalar(@samples);
		my $Nsmpls=0;
		foreach my $sd(keys %AGlist){
			$Nsmpls += scalar (@{$AGlist{$sd}}); #@{$AGlist{$cAssGrp}}
		}
		print "total samples: $Nsmpls , total in map: $Ndirs\n";
		print "\nSUBJOB ${subJob}/$maxSubJob: pre-restricted to " . scalar(@srtdSmpls)
			. " samples: @srtdSmpls\n\n";
	}
	
	
	
	print "Extracting GC genes from " . scalar(@srtdSmpls). " (of " . scalar(keys(%cl2gene2)) . ") ASsembly Groups\n";

	
	#different way to go over genes..
	 my $smCnt=1;
	 #storage hash for raw fasta/faa/link files, needs to be written separately
	#goes over every assembly group to extract SNP corrected genes that fall into each MGS
	my $writeLink = 1; my $appCnt=0;
		#DEBUG	@srtdSmpls = ("PDB3.F");
	
	
	foreach my $sm (@srtdSmpls){

		print "\nAT SMPL:: $smCnt/" . scalar(@srtdSmpls) ." $sm - ". "Elapsed time : ", timeNice(time - $sttime) . "\n";
		#readGenesSample_Singl($sm, $OFstrHR, $OAstrHR, $OCstrHR, $OLstrHR, $writeLink,$sttime);
		
		readGenesSample_Singl($sm, $writeLink,$sttime);
		$smCnt++; $appCnt++;
		
		if ($appCnt >= $appendWriteTrigger){
			appendWriteMGSgenes( $writeLink);
			$appCnt=0;
		}
	}
	
	
	appendWriteMGSgenes($writeLink);
	print "Done writing all genes to subdirs, elapsed time: " . timeNice(time - $sttime)  . "\n";
	$appCnt=0;
	#done at the point with gene extractions
	return;
}

sub createConsFastas{
	my ($cD,$sm, $oFNA, $oFAA,$append2LOG,$returnCmd) = @_;
	my $vcf2fnaBin = getProgPaths("vcf2fna");
	#my $normalizeVCF = abs_path(File::Spec->catfile(dirname(abs_path($0)), '..', 'SNP', 'normalizeVCFHeaders.pl'));
	#die "Cannot locate VCF header normalizer\n" unless defined($normalizeVCF) && -f $normalizeVCF;
	my $vcf2fnaOpt = "";
	#my $seqPlatf = "hiSeq"; #-> get this from .map ..
	my $refFA = getAssemblContigs($cD); my $refGFF = getAssemblGFF($cD);
	my $depthFile = "$cD$lMAPdir/$sm$bamDepthFsuffix";
	my $ofasCons = "$cD/$lSNPdir/$lConsCTG";
	my $vcfFile = "$cD/$lSNPdir/$lConsVCF";
	my $normalizedVCF = $vcfFile;#"$oFNA.input.vcf";
	#my @normalizeCommands = ("perl ".shellQuote($normalizeVCF)." -input ".shellQuote($vcfFile)." -output ".shellQuote($normalizedVCF)	);
	
	#DEBUG
	
	my $secSeqTechS = "";#secondary reads..
	my $support_reads = defined($map{$sm}{"SupportReads"}) ? $map{$sm}{"SupportReads"} : "";
	if ($support_reads =~ m/PB:/){$secSeqTechS = "PB" ;
	} elsif ($support_reads =~ m/ONT:/) {$secSeqTechS = "ONT" ;}
	my $seqPlatf = defined($map{$sm}{SeqTech}) ? $map{$sm}{SeqTech} : ""; #primary reads

	my $cmd ="";
	if ($seqPlatf eq ""){$seqPlatf = "hiSeq";} #if empty, assume hiSeq
	my $skipTerm = $noIndels ? " -skipINDELs" : "";
	my $commonOpt = "-t 1$skipTerm -minCallDepth $minSNPDepth -minCallQual $minSNPCallQual"
		. " -minCallQualAdaptive $useAdaptiveQual"
		. " -depthFilterScale $depthFilterScale -indelRange $indelRange";
	if ($secSeqTechS eq ""){
		#in case of only illumina:
		
		#checkSeqTech($seqPlatf);
		$vcf2fnaOpt = "-seqPlatform ".shellQuote($seqPlatf)." $commonOpt";
		$cmd = "$vcf2fnaBin $vcf2fnaOpt -ref ".shellQuote($refFA)
			." -inVCF ".shellQuote($normalizedVCF)." -depthF ".shellQuote($depthFile)."  ";
	} else {
		#die;
		#in case of both PacBio and illumina:
		#$vcf2fnaOpt = "-seqPlatform $SNPIHR->{SeqTech},$SNPIHR->{SeqTechSuppl} -t 1 -minCallDepth $minDepth,$minDepth -minCallQual $minCallQual ";
		#$cmd = "$vcf2fnaBin $vcf2fnaOpt -ref $refFA -inVCF $vcfFile,$vcfFileS -depthF $depthFile,$depthFileS ";# -oCtg $ofasCons.gz " ;
		my $vcfFileS = "$cD/$lSNPdir/$lConsVCFsup";
		my $normalizedVCFS = $vcfFileS;#"$oFNA.input.sup.vcf";
		#push @normalizeCommands,"perl ".shellQuote($normalizeVCF)." -input ".shellQuote($vcfFileS)." -output ".shellQuote($normalizedVCFS);
		my $depthFileS = "$cD$lMAPdir/$sm$bamDepthFsuffixSup";
		$vcf2fnaOpt = "-seqPlatform ".shellQuote("$seqPlatf,$secSeqTechS")." $commonOpt";
		$cmd = "$vcf2fnaBin $vcf2fnaOpt -ref ".shellQuote($refFA)
			." -inVCF ".shellQuote("$normalizedVCF,$normalizedVCFS")
			." -depthF ".shellQuote("$depthFile,$depthFileS")." -oCtg /dev/null ";
	}

	$cmd .= "-gff ".shellQuote($refGFF)." -oGeneNT ".shellQuote($oFNA)." -oGeneAA ".shellQuote($oFAA);
	#$cmd = join("\n", @normalizeCommands, $cmd);
	if ($append2LOG){$cmd.=" >> ".shellQuote($SNPconsLOGs)."\n";
	} else {$cmd .= "\n";}
	if ($returnCmd){ #don't excecute
		return $cmd;
	}
	
	#local excecution.. probably takes forever..
	#print "$cmd\n";
	#system "echo \$SLURM_LOCAL_SCRATCH";
	systemW $cmd;
}

sub readGenesSample_Singl{
	#go into curSpl dir and extract all marked gene reps.. 
	#write to correct format so they can be used in phylo later
	my ($sm, $writeLink,$sttime) = @_;
	#my %subG = %{$subGHR};#$_[0]};
	
	my %subG; my %locMGScnt;
	# This structure is read-only here; retaining the reference avoids copying
	# every locus entry at the start of each sample.
	my $locCl2G2 = $cl2gene2{$sm};

	my $noFilter =0;
	$noFilter = 1 if ($mode eq "MGSall");
	
	foreach my $gn (keys %{$locCl2G2}){
		#put genes into hash to avoid duplicates.. (locCl2G2{$gn} is now {member=>seed})
		foreach(keys %{$locCl2G2->{$gn}}){$subG{$_} = 1;}
		
		my $MGS = exists($LocusByID->{$gn}) ? $LocusByID->{$gn}{mgs} : undef;
		#stats collection on MGS usage
		if (defined $MGS){#exists($Gene2MGS->{$gn})){
			$locMGScnt{$MGS}++;
		}
	}
	print scalar(keys(%subG))." genes, " . scalar(keys(%locMGScnt)). " MGS\n";
	my @histoMGScnts = values %locMGScnt;
	histoMGS(\@histoMGScnts, "Possible Bins in sample");

	
	
	
	my $sd = $sm; #this is current sample
	my $sd2 = $sd;
#	my $writeLink = 1;
	if (exists(  $map{altNms}{$sd}  )){
		$sd2 = $map{altNms}{$sd}; $replN{$sd} = $sd2;
	}
	#print "SMMM: $sd $sd2 $replN{$sd}\n";
	#check if sample in map
		#print "map s: " .scalar(keys%map)."\n";

	unless (exists ($map{$sd2}) ) {
		limitedWarn('assembly groups absent from the map',
			"Can't find map entry for $sd; assembly group will be skipped\n");
		return;
	}
	my @subGKs = keys %subG;
	unless (@subGKs) {
		limitedWarn('samples without candidate genes',
			"No candidate genes found for sample $sm; sample will be skipped\n");
		return;
	}
	unless ($subGKs[0] =~ m/^(.*)__/) {
		limitedWarn('unparseable catalogue members',
			"Cannot parse catalogue member '$subGKs[0]' for sample $sm; sample will be skipped\n");
		return;
	}
	#find out if other samples are in the same assmblGrp..
	my @subSds = ($sd2);
	my $cAssGrp = $map{$sd2}{AssGroup};
	if (exists($AGlist{$cAssGrp})){
		@subSds = @{$AGlist{$cAssGrp}};
	}
	

	
	#print "$map{$sm}{SeqTech}\t2:$map{$sd2}{SeqTech}\t3:$map{$subSds[0]}{SeqTech}\n";
	
	#print "YY @subSds : $sd2 $sd\n";#die;
	#go into each sample ($sd3) from assembly group ($sd), that an assembly might be associated to (across multiple assemblies in assmblGrp)
	foreach my $sd3 (@subSds){
		if (exists $unavailableSamples{$sd3}) {
			limitedWarn('unavailable samples', "Skipping $sd3: $unavailableSamples{$sd3}\n");
			next;
		}
		unless (exists($map{$sd3}) && defined($map{$sd3}{wrdir}) && length($map{$sd3}{wrdir})) {
			limitedWarn('samples missing map entries or working directories',
				"Skipping $sd3: missing map entry or working directory\n");
			next;
		}
		#print "Time A: " . timeNice(time - $sttime)  . "\n";
		my $locSpace = "$locTmpDir/$sd3.cons/"; 
		#my $locSpace = "$preConDir/$sd3.cons/"; 
		
		my %locFAA; my %locFNA;my%locCSP;
		my %locMGSgenes; #keep track of genes written for each MGS..
		my $cD = $map{$sd3}{wrdir}."/";
		if (-e "$cD/SMPL.empty"){
			print ".. Empty->skip ";
			next;
		}
		my $rename = 0;
		$rename = 1 if ($sd2 ne $sd3);
		#print "r:$rename $sd3  (from $sd2) ";
		my $metaGD = getAssemblPath($cD,"",0);
		if ($metaGD eq ""){
			limitedWarn('samples without assemblies',
				"Assembly not available for $sd3 in $cD; sample will be skipped\n");
			next;
		}
		#get NT's
		#my $tar = $metaGD."genePred/genes.shrtHD.fna";
		
		#pre-calculated, as in old MGTK versions (pre 0.69):
		my $fastaf = "$cD/$lSNPdir/$lConsFNA";
		my $fastafAA = "$cD/$lSNPdir/$lConsFAA";
		my $fastafVCF = "$cD/$lSNPdir/$lConsVCF";
		my $locForceVCF2FNA=$forceVCF2FNA;
		
		if (exists($preCompSNPs{$sd3})
				&& fileGZe($preCompSNPs{$sd3}{NT}) && fileGZe($preCompSNPs{$sd3}{AA})){
			# The phase summary already reports precomputed consensus usage; avoid
			# printing one full filesystem path per sample here.
			$fastaf=$preCompSNPs{$sd3}{NT};
			$fastafAA=$preCompSNPs{$sd3}{AA};
			$locForceVCF2FNA=0;
		} elsif (exists($preCompSNPs{$sd3})) {
			limitedWarn('incomplete precomputed consensus files',
				"Ignoring incomplete precomputed consensus files for $sd3\n");
			delete $preCompSNPs{$sd3};
		}
		
		my $input_state = consensusInputState(
			fileGZe($fastafVCF), fileGZe($fastaf), fileGZe($fastafAA), $locForceVCF2FNA
		);
		if ($input_state eq 'missing') {
			limitedWarn('samples without repairable consensus files',
				"Skipping $sd3: consensus NT/AA files are incomplete and no repair VCF is available\n");
			next;
		}
		# Rebuild both members of the pair whenever either is missing.  Writing
		# into sample-local scratch avoids appending a new sidecar beside a .gz cache.
		if ($input_state eq 'regenerate'){
			make_path($locSpace) unless -d $locSpace;
			print "Recreating consensus fasta files on the fly.. ";
			#store these in scratch, uncompressed (much faster)
			$fastaf = "$locSpace/$sd3.cons.genes.fna";
			$fastafAA = "$locSpace/$sd3.cons.prots.faa";
			createConsFastas($cD, $sd3, $fastaf, $fastafAA, 1, 0);
		}
		#print "$fastaf\n";
		unless (fileGZe($fastaf) && fileGZe($fastafAA)){
			print "\n=====================================\nIncomplete consensus pair $fastaf / $fastafAA -> skip sample\n=====================================\n";
			#die;
			next;
		}
		#print "Time A1: " . timeNice(time - $sttime)  . "\n";
		#print "$fastaf\n";
		#read the assemble nt and AA genes from the sample
		my $FNA = readFasta($fastaf,1,"\\s",\%subG);
		#my %FNA = %{$hr};
		my $FAA2 = readFasta($fastafAA,0,"\\s",\%subG);# retain full headers for depth/CSP parsing
		my %FAA ;#= {};
		my %depths;
		#my $abunHR = readTabbed($cD.$abundF);
		#print "Time B: " . timeNice(time - $sttime)  . "\n";

		#my %FAA = %{$hr};
		#convert FAA hd
		my %conspSc;#read conspecific strain score from SNP consensus call..
		foreach my $k(keys %{$FAA2}){
			# Transfer, rather than copy, each sequence while normalizing its
			# header so the full-header hash shrinks throughout conversion.
			my $protein_sequence = delete $FAA2->{$k};
			#$k =~ m/^(\S+)\s.*CSP=([0-9\.]+)/;
			#requires vcf2fn v 0.25
			unless ($k =~ m/^(\S+)\sD=([0-9.]+)\s.*CSP=([0-9.]+)/) {
				limitedWarn('malformed consensus protein headers',
					"Malformed consensus protein header, skipping: $k\n");
				next;
			}
			my ($tmp, $depth, $csp) = ($1, $2, $3);
			$conspSc{$tmp} = $csp;
			$depths{$tmp} = $depth;
			$FAA{$tmp} = $protein_sequence;
		}
		$FAA2 = {};
		#print "Time C: " . timeNice(time - $sttime)  . "\n";

		#some stats on gene extractions..
		my $missGene=0; my $foundGene=0; my $SInum=0; my $conspGen=0;my $SNPresFail=0;
		my $doubleGenes=0; my $MGStoolowGskip=0;my $missAbundance=0;
		#stats on different ways to filter genes
		my $geneLost=0; my $conSpecFail=0; my $abundFail=0; my $doubleGsFail=0;
		
		
		#3rd part: genes were read and renamed.. now write them out already here to save mem overall
		#currently takes too long in large GCs..
		my $COGpriosZero=0;
		my $MGScnt = scalar((keys %locMGScnt));
		foreach my $MGS (keys %locMGScnt) {
			# The priority list is immutable during extraction, so do not copy
			# as many as $presortGenes entries for every sample/MGS pair.
			my $COGprios1 = $COGprios->{$MGS};
			if (!$COGprios1 || !@{$COGprios1}){
				$COGpriosZero++;
				next;
			}

			my $locConSpecGen=0; my $accAbu=0; my $LmissG=0; my $doubleCntL=0;
			my $LmuissAbu=0; my $evaluableLoci=0;
			my @genes2 = (); #stores semi-final list of genes
			my @abunGs = (); #abundance vector of genes
			my %curLocus;
			my %linkStr; #temp storage for links to gene cat etc of catalogues genes
			foreach my $locus (@{$COGprios1}){
				next unless length($locus) && exists($locCl2G2->{$locus});
				my $membersHR = $locCl2G2->{$locus}; #{member => seed}
				my @genes = keys %{$membersHR};
				my @candidates;
				my $had_evaluable = 0;
				my $csp_rejected = 0;
				for my $gX (@genes) {
					next unless length $gX;
					if (!exists($FAA{$gX}) || !exists($FNA->{$gX})) {
						$LmissG++;
						next;
					}
					my $depth = $depths{$gX};
					if (!$noFilter && (!defined($depth) || $depth < $minDepthGene)) {
						$LmuissAbu++;
						next;
					}
					$depth = 1 unless defined($depth) && $depth > 0;
					$had_evaluable = 1;
					if (!$noFilter && defined($conspSc{$gX}) && $conspSc{$gX} > $conspecificSpThr) {
						$csp_rejected++;
						next;
					}
					push @candidates, {
						id => $gX, protein => $FAA{$gX}, depth => $depth,
						seed => $membersHR->{$gX},
						context => $MemberContext->{$gX} || {},
					};
				}
				$evaluableLoci++ if $had_evaluable;
				if ($had_evaluable && !@candidates && $csp_rejected) {
					$locConSpecGen++;
					next;
				}
				next unless @candidates;

				my $group = $LocusByID->{$locus};
				# A unique viable candidate is always selected by
				# choose_locus_candidate.  Most loci take this path, so skip
				# the otherwise-unused protein k-mer and context scoring.
				my $selection;
				if (@candidates == 1) {
					$selection = {
						status => 'selected', candidate => $candidates[0], reason => 'unique',
					};
				} else {
					# Cache this invariant map lazily: ambiguous loci can recur
					# across samples, while unique loci need no extra storage.
					$LocusSeedProteins{$locus} ||= {
						map {
							defined($catalogProteins->{$_}) ? ($_ => $catalogProteins->{$_}) : ()
						} @{$group->{genes}}
					};
					$selection = choose_locus_candidate(
						\@candidates,
						$LocusSeedProteins{$locus},
						$LocusContext->{$locus},
					);
				}
				if ($selection->{status} ne 'selected') {
					$doubleCntL++;
					next;
				}
				my $curG = $selection->{candidate}{id};
				my $bestAB = $selection->{candidate}{depth};
				push @genes2, $curG;
				$curLocus{$curG} = $locus;
				push @abunGs, $bestAB;
				$accAbu += $bestAB;
				my $laterHd = "$sd3$SaSe" . externalLocusName($locus, $MGS);
				$linkStr{$curG} = "$laterHd\t$locus\t$group->{primary_gene}\t".scalar(@genes)
					."\t".join(",",@genes)."\t$selection->{reason}\n";
			}

			$doubleGenes += $doubleCntL;
			$conspGen+=$locConSpecGen;
			$missGene += $LmissG;
			$missAbundance += $LmuissAbu;
			my $double_failure = !$noFilter && $evaluableLoci > 0 && $doubleCntL >= $minBadLociForSampleSkip
				&& ($doubleCntL / $evaluableLoci) > $multiGeneSmplMax;
			my $csp_failure = !$noFilter && $evaluableLoci > 0 && $locConSpecGen >= $minBadLociForSampleSkip
				&& ($locConSpecGen / $evaluableLoci) > $conspGeneSmplMax;
			if ($double_failure || $csp_failure){
				push(@{$ConspecificMGS{$MGS}}, "$sd3" ); 
				$doubleGsFail++ if $double_failure;
				$conSpecFail++ if $csp_failure;
				next;
			}
			next unless @genes2;

			my $depth_mask = $noFilter ? [(1) x scalar(@abunGs)] : robust_depth_mask(\@abunGs);
			my @genes3=();
			for (my $i=0;$i<scalar(@abunGs);$i++){
				unless ($depth_mask->[$i]) {
					$abundFail++; next;
				}
				push (@genes3, $genes2[$i]);
			}
			
			if (scalar(@genes3)< $MGStoolowGsThr){
				$MGStoolowGskip++;next;
			}
				
			#now write MGS into local temp storage for later tree building..
			my $locCnt=0;
			my @OCstr; my @OFstr; my @OAstr; my @OLstr ;
			foreach my $gX (  @genes3 ){
				unless (exists($FAA{$gX}) && exists($FNA->{$gX})){
					limitedWarn('catalogue genes absent from consensus sequences',
						"Could not find '$gX' gene in consensus sequences\n");
					next;
				}
				my $strCpy = ""; $strCpy = $FAA{$gX};# if (exists($locFAA{$gX}));
				my $AAlen = 0; $AAlen = int(length($strCpy)) if (defined($strCpy));
				if ($AAlen == 0){$SNPresFail++; next;}
				my $num1 = $strCpy =~ tr/\-Xx//;
				if ($num1 >= ($AAlen-1)){ $SNPresFail++; next;} #all X, exclude..
				if (!$noFilter){
				if ($locCnt >= $maxNGenes){ next;}
				}
				
				$locCnt++;
				#write gene out
				my $ng = "$sd3$SaSe" . externalLocusName($curLocus{$gX}, $MGS);
				# Tree-facing identifier: sample|COG|primary catalogue gene.
				#die;
				push(@OFstr , ">$ng\n$FNA->{$gX}\n"); #FNA
				push(@OAstr ,">$ng\n$strCpy\n"); #FAA
				#add to category for later..
				push(@OCstr , "$MGS\t$curLocus{$gX}\t$sd3\t$ng\n");
				#$SIcat{$MGS}{$cog}{$sd3} = $ng;
				#$genesWrite{$MGS}++;
				
				if ($writeLink){
					push(@OLstr, $linkStr{$gX});
				}
			}
			
			
			if (scalar(@OFstr) == 0 || $locCnt < $MGStoolowGsThr){ #5 genes is really too little to be considered valid as good strain rep..
				$MGStoolowGskip++;
				#delete $locMGSgenes{$MGS};
				next;
			}
			$locMGSgenes{$MGS} = $locCnt;
			
			if (!exists($OAstrH{$MGS})){#set up base strings
				$OAstrH{$MGS} = "";$OFstrH{$MGS} = "";$OLstrH{$MGS} = "";$OCstrH{$MGS} = "";
			}
			if ($locCnt>0){
				#save in tmp hash (faster than opening bunch of files..
				$OAstrH{$MGS} .= join("",@OAstr);$OFstrH{$MGS} .= join("",@OFstr);
				$OLstrH{$MGS} .= join("",@OLstr);$OCstrH{$MGS} .= join("",@OCstr);
				#push(@{$OAstrH{$MGS}},join("",@OAstr));push(@{$OFstrH{$MGS}},join("",@OFstr));
				#push(@{$OLstrH{$MGS}}, join("",@OLstr));push(@{$OCstrH{$MGS}},join("",@OCstr));
				$SInum ++ ;
				$foundGene+=$locCnt;
			}
			#clenup tmp
		} #loop over MGS
		#print "Time D: " . timeNice(time - $sttime)  . "\n";
		remove_tree($locSpace) if -d $locSpace;

		my @genesPmgs = values %locMGSgenes; 	@genesPmgs = sort { $a <=> $b}  @genesPmgs;
		histoMGS(\@genesPmgs,"Detected Bin Genes:");
		
		print "$sd3 - Missed/MissAbund/lost/abundFilterFail/SNPresFail Gs: ${missGene}/${missAbundance}/${geneLost}/${abundFail}/$SNPresFail\tConspecGs/consMGS/doublGs/failcMGS: ${conspGen}/$conSpecFail/${doubleGenes}/$doubleGsFail\tFoundGs: $foundGene/". scalar(keys %FAA) . "\tused MGS/skipped MGS: ${SInum}/$MGStoolowGskip\t";
		if (@genesPmgs) {
			print "GperMGS (median,mean): " . median(@genesPmgs) . "/". int(mean(@genesPmgs)+0.5);
		} else {
			print "GperMGS (median,mean): 0/0";
		}
		if ($COGpriosZero>=$MGScnt*0.95){print " $COGpriosZero / $MGScnt no COGprio list! ";}
		print "\n";
	}
}
