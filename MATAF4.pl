#!/usr/bin/env perl
#The Metagenomic Assembly, Genomic Recovery and Assembly Independent Mapping Tool (MATAFILER)
#main MATAFILER routine
# (c) Falk Hildebrand, 2016-2025
#examples
#./MATAF4.pl map2tar test/refCtg.fasta,test/refCtg.fasta test1,test2
#./MATAF4.pl map2tar test/TEC2/v5/TEC2.MM4.BEE.GF.rn.fa TEC2
#./MATAF4.pl -map dir/map

use warnings;
use strict;
use File::Basename;
use File::Find ();
use File::Path qw(make_path);
use File::Spec;
use Cwd 'abs_path';
use POSIX;
use Getopt::Long qw( GetOptions );
use List::Util qw(max sum);
use Text::Wrap qw(wrap);
use Text::ParseWords qw(shellwords);

use vars qw($CONFIG_FILE);


#load MF specific modules
use Mods::GenoMetaAss qw(readMap readMapS getDirsPerAssmblGrp checkAssmblGrp lcp readFastHD prefixFAhd prefix_find 
			gzipopen fileGZe fileGZs contig_stats_coverage_complete
			readFasta writeFasta systemW getAssemblPath  filsizeMB resetAsGrps
			iniCleanSeqSetHR checkSeqTech is3rdGenSeqTech hasSuppRds 
			addFileLocs2AssmGrp getRawLibrariesAssmGrp getCleanLibrariesAssmGrp
			parseSupportReads discoverReadFiles);
use Mods::ReadLibrary qw(
	newReadLibrary cloneReadLibraries readLibrariesFromArrays
	ensureSeqSetLibraries ensureCleanSeqSetLibraries
	syncSeqSetLegacy syncCleanSeqSetLegacy replaceScopeLibraries
	readLibrariesByScope libraryFiles libraryPairs libraryTechnology
);
use Mods::IO_Tamoc_progs qw(getProgPaths setConfigFile jgi_depth_cmd inputFmtSpadesLibraries inputFmtMegahitRuntimeLibraries createGapFillopt
			buildMapperIdx mapperDBbuilt decideMapper  checkMapsDoneSH greaterComputeSpace);
use Mods::SNP qw(SNPconsensus_vcf SVcall_vcf);
use Mods::TamocFunc qw (cram2bsam getSpecificDBpaths getFileStr displayPOTUS bam2cram checkMF checkMFFInstall);
use Mods::phyloTools qw(fixHDs4Phylo);
use Mods::Binning qw (getBinSubdirName binningOutputsComplete );
use Mods::Subm qw (qsubSystemWaitMaxJobs qsubSystem emptyQsubOpt findQsubSys qsubSystemJobAlive MFnext add2SampleDeps numUserJobs numLiveUserJobs numActiveUserJobs recordSampleLockJobs sampleLockActiveJobs primeSampleLockJobSnapshot slurmJobFailureSummary submitSlurmWithDependencyRecovery deferredSubmissionDependency submissionDependencyDeferred handleSubmissionFailure);
use Mods::WorkflowState qw(inspect_workflow_state encode_state_report);
use Mods::WorkflowPlan qw(build_workflow_plan encode_workflow_plan);
use Mods::WorkflowRunner qw(run_workflow_preflight);
use Mods::WorkflowControl qw(
	advance_loop_window overlap_loop_window rolling_completed_frontier priority_outputs_complete parse_loop_spec should_rerun_locked_window assembly_cores_for_input assembly_group_output_dirs parse_ignored_samples
	balanced_parallel_batches
	hybrid_group_ready hybrid_package_complete hybrid_package_sample_id missing_input_files source_input_files
	hybrid_local_scratch_gb
	sample_base_output_dir sample_is_ignored workflow_members_match
	normalise_job_dependencies append_job_dependencies deferred_command_dependencies augment_deferred_submission
	commands_are_lightweight_filesystem cleanup_stage_barrier
);


#some useful HPC commands..
#squ | grep 'r' | cut -f11 -d' ' | xargs -t -i scontrol update TimeLimit=84:00:00 jobid={}
#squ | grep 'dencyNev' | cut -f11 -d' ' | xargs  -t -i scancel {}



#local subs
sub announce_MF4;
sub smplStats; sub checkDrives; 
sub isLastSampleInAssembly;
sub cleanupCompletionRequirements; sub finishedCleanupArguments; sub runFinishedCleanup; sub submitFinishedCleanup;
sub uploadRawFilePrep; sub unploadRawFilePostprocess;
sub seedUnzip2tmp; sub cleanInput; #unzipping reads; removing these at later stages ; remove tmp dirs

sub sdmClean; sub sdmOptSet;  #qual filter reads
sub mergeReads; #merge reads via flash
sub removeHostSeqs; sub krakenTaxEst;sub prepKraken;
sub loop2C_check;
sub primeLoopSchedulerSnapshot;
sub deferLoopProducerWave;

sub metagAssemblyRun;
sub buildAssemblyMapIdx;
sub createPsAssLongReads; #pseudo assembler
sub prepPreAssmbl; #hybrid pacbio/ill assemblies
sub genePredictions; sub run_prodigal; #gene prediction
sub mapReadsToRef;  sub scndMap2Genos;
sub contigStatsOutputsComplete; sub runContigStats;#sub bam2cram;
sub mocat_reorder; sub postSubmQsub; 
sub detectRibo;  sub riboSummary;

sub checkRawProgsFin; sub prepareDiamondRerun; sub publishKrakenResults;

sub runOrthoPlacement;
sub runDiamond; sub DiaPostProcess;
sub nopareil; sub calcCoverage2nd;sub d2metaDist; 
sub metphlanMapping; sub mergeMP2Table;
sub mOTU2Mapping; sub mergeMotu2Table; sub prepMOTU2;
sub genoSize; sub check_map_done; sub check_depth_done;
sub mapping_reference_matches;
sub postprocess;
sub reportSlurmJobFailures;
sub setDefaultMFconfig;
sub getCmdLineOptions;
sub setupHPC;
sub runStateInspection;
sub runAutomaticWorkflowPreflight;
sub workflowStateOptions;
sub sampleReadSet;
sub discoverSampleInputs; sub populateInputSizesFast; sub spaceInAssGrp;

sub createConsSNPandSVs;


#------- version history MATAFILER / MATAFILER --------
#.75: 4.3.26: ini MATAFILER4 version
#4.01: 13.3.26: removing bugs from hybrid assembly detection, switching to 4.x versioning
#4.02: 15.4.26: updated internal logic for passing read paths, enabled hybrid mode in complex assembly groups
#4.03: 18.4.26: further fixes to hybrid assembly logic. separateCongigs.pl tech hardened.
#4.04: 24.4.26: hybrid assembly logic, adapting GC to different binners
#4.10: 16.5.26: CahtGPT 5.6 sol: completely rewrite of calling logic, multiple bugs fixed across various scripts, better metagStats reporting, hybrid assemblies strengthned.
#4.11: 21.7.26: loopTillComplete overlaps light passes and merges the next block on each final pass.
#       Retained-lock windows get a bounded retry when at least one but fewer than
#       1% of their samples or at most the configured number of jobs remain active.
#       Loop specifications and incompatible rewrite modes fail early.
#4.12: 22.7.26: reconcile aged Slurm dependencies through accounting after MinJobAge,
#       omitting successfully completed jobs while reporting failed or unknown jobs.
#       Track fresh submission times to avoid routine accounting queries, and remove
#       unnecessary cross-sample ContigStats-to-ContigStats dependencies.
#4.13: 24.7.26: loopTillComplete distinguishes executing jobs from queued
#       dependencies and starts its next pass at a configurable active-job threshold.
#4.14: 26.7.26: cache one validated input discovery per sample for sizing and staging.
#4.15: 26.7.26: group run state and unify raw/clean read sets around library records.
#4.16: 26.7.26: consolidate runtime options and checkpoint paths into named hashes.
#4.17: 26.7.26: replace per-sample RMLOCK scheduler jobs with job-ID lock ledgers
#       that MATAF4 safely releases after their recorded jobs leave the scheduler.
#4.18: 26.7.26: summarize Slurm failures, including OOM and timeout outcomes,
#       as an occurrence matrix grouped by MATAFILER job category.
#4.19: 26.7.26: defer sample scratch and LOGandSUB creation until unfinished
#       work is confirmed, keeping completed-sample passes free of staging churn.
#4.20: 26.7.26: centralize completed-sample filesystem publication and cleanup
#       policy in cleanup_finished_sample.pl.
#4.21: 26.7.26: report primary and supplementary input sizes independently
#       and reject physical input files assigned to both read scopes. Complete
#       the user-facing transition from the former toolkit name to MATAFILER.
#4.22: 28.7.26: complete assembly-independent workflows without assembly
#       checkpoints, release their sample scratch safely, and reject assembly-only
#       binning or variant options when assembly is disabled.
#4.23: 30.7.26: use a rolling completed-sample frontier, a final full-range
#       verification pass, shared batched Slurm lock/accounting snapshots, and
#       running-plus-pending admission control. Repeated completed-sample visits
#       use an ordered priority-output probe before deeper filesystem checks.
#4.24: 30.7.26: prevent loopTillComplete from blocking indefinitely at the
#       live-job cap; defer only submissions, retain sample cleanup checks, and
#       continue past failed deferred dependency chains.
#4.25: 30.7.26: submit loopTillComplete producers in readiness waves, so input
#       staging, quality/host filtering, assembly, mapping, and contig statistics
#       complete before their consumers are submitted.
#4.26: 31.7.26: after the first final full-range verification, wait for every
#       job submitted by this invocation to leave the scheduler before starting
#       another full-range pass.
my $MATFILER_ver = 4.26;

#----------------- defaults ----------------- 

# Runtime controls set directly by the command line.  Keeping these together
# distinguishes invocation state from persistent pipeline configuration.
my %runOptions = (
	operationMode => @ARGV > 0 ? $ARGV[0] : "",
	sharedTmpDir => "",
	nodeTmpDir => "",
	baseID => "",
	from => 0,
	to => 999999999999,
	submit => 1,
	loopCount => "0",
	loopInitialCount => 0,
	loopWindowSize => 0,
);

#runtime paths
my $logDir = ""; #this is the local logdir
my $baseOut = "";


#counter on where MF is in map
my $JNUM=0;

#config is more overall configuration for MATAFILER
my %MFconfig; 

#MFcontstants: object to store essential paths/file endings
my %MFcontstants;

#MFopt: global object with options for MATAFILER. Added in MF v0.5, slowly rebuild MF around this system
my %MFopt; 

# Keep machine-readable inspection and plan output clean. When either read-only
# mode was requested, ordinary diagnostics use STDERR and JSON is written
# explicitly to STDOUT.
my $readOnlyStateRequested = 0;
for (my $argI = 0; $argI < @ARGV; $argI++) {
	if ($ARGV[$argI] =~ m{^--?(?:inspectState|planState)=(\d+)$}) {
		$readOnlyStateRequested = 1 if ($1 != 0);
	} elsif ($ARGV[$argI] =~ m{^--?(?:inspectState|planState)$}) {
		$readOnlyStateRequested = 1
			unless ($argI + 1 < @ARGV && $ARGV[$argI + 1] eq '0');
	}
}
select STDERR if ($readOnlyStateRequested);


#keep track of DBs that the metagenome will be filtered against..
my @filterHostDB = ();
#track secondary mapping and ref DBs
my %map2ndTogRefDB;my %make2ndMapDecoy;
my @bwt2outD =(); my @DBbtRefX = (); my @DBbtRefGFF=(); my @bwt2ndMapNmds;
my @scaffTarExternalOLib1; my @scaffTarExternalOLib2;

my @EBIjobs = (); #keeps track of $MFconfig{uploadRawRds} jobs, submits postprocessing (md5)
#----------- map all reads to a specific reference - options ---------

#progStats: object to track progress of programs/submissions
my %progStats;#count up progress of submitted jobs in current run

#HDDspace: object to handle HDD usage: Always format as "XXG" XX = space requirements in Gb. Excecption: "-1"
my %HDDspace;

my %locStats; #keeps statistics of samples in hash (from already finished samples
#keeps globally used vars (mostly dirs/paths)
my %MFglobal;

#fixed checkpoint names
my %checkpointNames = (
	preAssemblyDone => "preassmblDone.sto",
	assemblyDone => "ass.done.sto",
);
#preAssebmly (hybrid assemblies etc) is done, main assembly

#fixed dirs for specific set of samples
my %preDIRs = (dir_ContigStats => "/assemblies/metag/ContigStats/",dir2MePhl => "pseudoGC/Phylo/MP2/",dir2RiboF => "pseudoGC/Phylo/RiboFind/");

# Show help before initialising defaults or checking the installation.  Apart
# from keeping the output clean, this makes the option reference available on
# systems that have not been configured for a pipeline run yet.
my $helpRequested = grep { /^(?:--?help|-h|-\?)$/ } @ARGV;

#say hello to user 
announce_MF4();
help() if ($helpRequested);
setDefaultMFconfig();
getCmdLineOptions;
die "-loopTillCompleteActiveJobs requires a non-negative integer\n"
	if ($MFconfig{loopTillCompleteActiveJobs} < 0);
die "-schedulerPollSeconds requires a positive integer\n"
	if ($MFconfig{schedulerPollSeconds} < 1);
die "-assemblCores requires a non-negative integer (0 enables automatic scaling)\n"
	if ($MFopt{AssemblyCores} < 0);
die "-minBinnerAssemblyMB requires a non-negative number\n"
	if ($MFopt{minBinnerAssemblyMB} < 0);
$MFconfig{inspectState} = 1 if ($MFconfig{planState});
if (!$MFconfig{inspectState} && $runOptions{loopCount} && (
		$MFopt{redoAssMapping} || $MFopt{BinnerRedoAll} || $MFopt{redoAssembly}
		|| $MFopt{redoSNPcons} || $MFopt{redoSNPgene}
		|| $MFopt{rewriteAllIfAnyDiamond} || $MFopt{rewriteDiamond}
		|| $MFopt{RedoRiboFind} || $MFopt{RedoRiboAssign}
	)) {
	die "MATAFILER rewrite options cannot be combined with -loopTillComplete; " .
		"disable rewrite options before starting a looped run.\n";
}
checkMF(1) unless ($MFconfig{inspectState});



#programs of global (pun) importance --------------------------
my $smtBin = getProgPaths("samtools");#
my $pigzBin  = getProgPaths("pigz");
my $avx2Constr =  getProgPaths("avx2_constraint",0);

 
#set up link to submission system on cluster. Inspection mode deliberately
#does not initialise or query the scheduler.
my $QSBoptHR;
$QSBoptHR = setupHPC() unless ($MFconfig{inspectState});

#----------- map all reads to a specific reference, preparation ---------
my $map2ndMpde=0;#0=map2tar;2=map2DB;3=map2GC
#-----------   scaffolding external contigs parameters
my $scaffTarExternal = "";my $scaffTarExternalName = ""; 
my $scaffTarExtLibTar = ""; my $bwt2ndMapDep = ""; 

# the map and some base parameters (base ID, in path, out path) can be (re)set
my %map; my %AsGrps; my %DOs;#DOs only required for metabat, to use all mappings within an assembly group
my $workflowIteration = 0;
my $ignoredSamplesHR = parse_ignored_samples($MFconfig{ignoreSmpl});
prepareMap();
my @samples = @{$map{opt}{smpl_order}}; 
populateInputSizesFast($_) for @samples;

if ($MFconfig{inspectState}) {
	runStateInspection();
	exit(0);
}
runAutomaticWorkflowPreflight($workflowIteration) if ($MFconfig{autoStatePlan});


# Central reporting state. Repeated workflow passes replace a sample's
# snapshot without losing samples collected in earlier loop windows.
my %runReport = (
	samples => {},
	order => [],
	seen => {},
	context => {},
	empty_samples => {},
	present_assemblies => 0,
);
#the sample currently worked on, important to keep as a global variable
my $curSmpl = "";
#compare to previous dir
my $baseoutPrev = "";
# Inputs accumulated for cross-sample d2 distance calculation.
my %d2Inputs = (
	samples => {},
	filtered_read1 => [],
	filtered_read2 => [],
	dependencies => "",
);

my @unzipjobs; 
my @grandDeps; #used for loop2completion , collects dependencies 

#profiler / human filtering preps
my $krakDeps = prepKraken();
prepMetaphlan();


if ($runOptions{to} > @samples){
	print "Reset range of samples to ". @samples."\n"; $runOptions{to} = @samples;
}
my $from = $runOptions{from}; my $to = $runOptions{to};
my ($selectedFrom, $selectedTo) = ($runOptions{from}, $runOptions{to});
#die "\"@samples\"\n";
if ($runOptions{loopWindowSize} > 0){
	$to = $from + $runOptions{loopWindowSize}; $to = $runOptions{to} if ($to > $runOptions{to});
}
my $loopIterationSubmissionStart = $QSBoptHR->{submittedJobs} || 0;
my $loopIterationExtended = 0;
my $loopFinalLockRetryUsed = 0;
my %loopSubmittedJobIds;
my %loopSampleCompleted;
my $loopFinalVerification = 0;
my $loopSawActiveLocks = 0;
$QSBoptHR->{nonblockingMaxConcurrentJobs} = $runOptions{loopCount} ? 1 : 0;
primeLoopSchedulerSnapshot($from, $to) if $runOptions{loopCount};


#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------

#die $to."\n";
#for loop that goes over every single sample in the .map
for ($JNUM=$from; $JNUM<$to;$JNUM++){
	
	unless ($runOptions{loopCount}) {
		qsubSystemWaitMaxJobs(
			$MFconfig{checkMaxNumJobs}, $MFconfig{killDepNever}, $QSBoptHR,
		);
	}
	$curSmpl = $samples[$JNUM];
	
	
	#set up initial local paths for a given sample
	my $dir2rd=""; my $curDir = "";
	my $curOutDir = $map{$curSmpl}{wrdir};	
	#key IDs for sample
	my $cAssGrp = $curSmpl;my $cMapGrp = $map{$curSmpl}{MapGroup};
	if ($map{$curSmpl}{AssGroup} ne "-1"){ $cAssGrp = $map{$curSmpl}{AssGroup};}
	my @sampleDeps = (); #catalogues all job dependencies created in this loop
	

	#local flow control
	#print "SMPL::$curSmpl\n";
	$dir2rd = $map{$curSmpl}{dir};  $dir2rd = $map{$curSmpl}{prefix} if ($dir2rd eq "");
	my $SmplName = $map{$curSmpl}{SmplID};
	if ($dir2rd eq "" ){#very specific read dir..
		if ($map{$curSmpl}{SupportReads} ne ""){
			$curDir = "";	#$curOutDir = $map{$curSmpl}{wrdir};	#$curOutDir = "$baseOut$SmplName/";	
		} else { die "Can;t find valid path for $SmplName\n";
		}
		$dir2rd = $SmplName;
	} else {
		$curDir = $map{$curSmpl}{rddir};	
	}
	
	$baseOut = sample_base_output_dir($curOutDir, $curSmpl);
	if (!-d $baseOut || $baseoutPrev ne $baseOut){
		$MFglobal{globalLogDir} = $baseOut."LOGandSUB/"; #this is the gloabl logdir (across all samples in current run)
		system("mkdir -p $MFglobal{globalLogDir}/sdm") unless (-d "$MFglobal{globalLogDir}/sdm");
		open $QSBoptHR->{LOG},">",$MFglobal{globalLogDir}."qsub.log";# unless ($runOptions{submit} == 0);
		$MFglobal{collectFinished} = $baseOut."runFinished.log";
		foreach (split /,/,$MFconfig{mapFile}) {system "cp $_ $MFglobal{globalLogDir}/";}
		print $MFglobal{globalLogDir}."qsub.log\n";
		$baseoutPrev = $baseOut;
	}

	$logDir = "$curOutDir/LOGandSUB/";

	#ignore samples .. for various reasons ------------------------------------------------------------------------------------
	if ($MFconfig{ignoreSmpl} ne ""){
		if (sample_is_ignored($ignoredSamplesHR, $SmplName)){print "\n ======= Ignoring sample $SmplName =======\n";
		loop2C_check($cAssGrp,\@sampleDeps);next;}
	}
	my $smplLockF = "$logDir/$MFcontstants{DefaultSampleLock}";
	if (-e $smplLockF){
		if ($MFconfig{rmSmplLocks}){ system "rm -f $smplLockF";
		} else {
			my $activeLockJobs = sampleLockActiveJobs($smplLockF, $QSBoptHR);
			if (defined($activeLockJobs) && $activeLockJobs == 0) {
				unlink $smplLockF
					or die "Cannot release completed sample lock $smplLockF: $!\n";
				print "Released completed sample lock for $SmplName\n"
					unless $MFconfig{silent};
			} else {
				$loopSawActiveLocks = 1
					if defined($activeLockJobs) && $activeLockJobs > 0;
				my $lockReason = defined($activeLockJobs)
					? "$activeLockJobs recorded job(s) still queued or running"
					: "legacy/unreadable lock has no verifiable job ledger";
				print "\n    >>>>>>>>>> Sample $SmplName is locked: $lockReason <<<<<<<<<<  \n";
				loop2C_check($cAssGrp,\@sampleDeps);next;
			}
		}
	} 
	$QSBoptHR->{LOCKfile} = $smplLockF; #set lockfile to be created if any job is submitted
	

	%locStats = ();
	#$locStats{hasPaired} = 0;	$locStats{hasSingle} = 0; 
	$d2Inputs{samples}{$SmplName}=1;
	
	print "\n======= $SmplName - $JNUM - $dir2rd =======\n" unless($MFconfig{silent});
	
	#set up dirs ------------------------------------------------------------------------------------
	my $smplTmpDir = "$MFglobal{runTmpDirGlobal}$SmplName/"; #curTmpDir
	my $nodeSpTmpD = "$runOptions{nodeTmpDir}/$SmplName";

	my $finalCommAssDir = "$curOutDir/assemblies/metag/";
	my $finalCommAssDirSingle = $finalCommAssDir; #this is only used for checking..
	my $finalMapDir = "$curOutDir/mapping/";
	#$DBpath="$curOutDir/readDB/";
	
	
	$AsGrps{$cAssGrp}{AssemblSmplDirs} .= $curOutDir."\n";
	my $AssemblyGo=0; my $MappingGo=0; #controls if assemblies / mappings are done in respective groups
	
	#complicated flow control for multi sample assemblies
	die "cAssGrp eq  \"$cAssGrp\" ".$runOptions{baseID}."\n" if ($cAssGrp eq "");
	my $assmGrpTag = "AssmblGrp_$cAssGrp";
	#die "AssemblGrp exists already!!: \"$assmGrpTag\"\n" if (exists($assmblGrpLog{$assmGrpTag}));
	#$assmblGrpLog{$assmGrpTag} = 1;

	$finalCommAssDir = "$baseOut/$assmGrpTag/metag/" if ( !exists($AsGrps{$cAssGrp}{CntAimAss}) || $AsGrps{$cAssGrp}{CntAimAss}>1);
	my $asmDir = $finalCommAssDir;
	$asmDir =~ s/metag\/$//;
	my $metaGpreAssmblDir = "$curOutDir/assemblies/pre$assmGrpTag/"; #sample-specific durable hybrid-assembly handoff

	#assign job name (dependency) only ONCE
	if ( !exists($AsGrps{$cAssGrp}{CntAss}) || $AsGrps{$cAssGrp}{CntAss} == 0){
		$AsGrps{$cAssGrp}{AssemblJobName} = "";
	}
	
	$AsGrps{$cAssGrp}{CntAss} ++;
	print "AssmblGrp: " . $AsGrps{$cAssGrp}{CntAss} .":".$AsGrps{$cAssGrp}{CntAimAss}. ";\n"
		if ($AsGrps{$cAssGrp}{CntAimAss} > 1 && !$MFconfig{silent});
	if ($AsGrps{$cAssGrp}{CntAss}  >= $AsGrps{$cAssGrp}{CntAimAss} ){
		#print "running assembluy\n";
		$AssemblyGo = 1;
	}
	
	
	#mapping groups?
	$AsGrps{$cMapGrp}{CntMap} ++;
	#my $hasPrimaryRds= 1;$hasPrimaryRds = 0 if ($map{$curSmpl}{prefix} eq "" && $map{$curSmpl}{dir} eq "");
	

	print "MapGroup: ".$AsGrps{$cMapGrp}{CntMap} .":".$AsGrps{$cMapGrp}{CntAimMap}.";\n"
		if ($AsGrps{$cMapGrp}{CntAimMap} > 1 && !$MFconfig{silent});
	if (!exists($AsGrps{$cMapGrp}{CntMap})){ die "Can;t find CntMap for $cMapGrp";}
	if ($AsGrps{$cMapGrp}{CntMap}  >= $AsGrps{$cMapGrp}{CntAimMap} 
			&& $map{$curSmpl}{hasPrimaryRds} ) {  #ensure primary reads are present and registered
		#print "running mapping";
		$MappingGo = 1;
	}
	#die "MG: $MappingGo\n";
	if ($MFopt{DoAssembly} ==0 ){$AssemblyGo=0;}
		#mapping related
	my $cramthebam = 1;
	my $bamcramMap = "bam"; if ($cramthebam){$bamcramMap = "cram";}

	
	#----------------------------  2  ----------------------------------
	#set up dirs for this sample
	my $metagAssDir = $asmDir."metag/";
	my $geneDir = $metagAssDir."genePred/";
	#system("mkdir -p $metagAssDir $geneDir");
	my $metaGassembly=$metagAssDir."scaffolds.fasta.filt"; 
	my $finalCommScaffDir = "$finalCommAssDir/scaffolds/";
	my $metaGscaffDir = "$metagAssDir/scaffolds/";
	my $STOfinScaff = "$finalCommScaffDir/scaffDone.sto";
	my $pseudoAssFile = "$metagAssDir/longReads.fasta.filt";
	my $pseudoAssFileFinal = "$finalCommAssDir/longReads.fasta.filt";
	my $finAssLoc = "$finalCommAssDir/scaffolds.fasta.filt";
	my $ContigStatsDir  = "$curOutDir/$preDIRs{dir_ContigStats}/";
	my $coveragePerCtg = "$ContigStatsDir/Coverage.percontig.gz";
	my $markerGenesPerCtg = "$finalCommAssDir/ContigStats/GTDBmg/marker_genes_meta.tsv"; #GTDB marker genes
	my $suppCoveragePerCtg = "$ContigStatsDir/Cov.sup.percontig.gz";
	my $nonParDir = $curOutDir."nonpareil/";
	#SNP calling on assembly related files
	my $SNPdir = "$curOutDir/SNP/";
	my $contigsSNP = "$SNPdir/contig.SNPc.$MFopt{SNPcallerFlag}.fna"; #flag $MFopt{saveConsFastas} controls this.. #keep without gz, although will be gz'd
	my $genePredSNP = "$SNPdir/genes.shrtHD.SNPc.$MFopt{SNPcallerFlag}.fna.gz";
	my $genePredAASNP = "$SNPdir/proteins.shrtHD.SNPc.$MFopt{SNPcallerFlag}.faa.gz";
	my $vcfSNP = "$SNPdir/allSNP.$MFopt{SNPcallerFlag}.vcf";	$vcfSNP = "" if (!$MFopt{saveVCF});
	my $vcfSNPsupp = "$SNPdir/allSNP.$MFopt{SNPcallerFlag}-sup.vcf";	$vcfSNPsupp = "" if (!$MFopt{saveVCF});
	my $SVdir = "$curOutDir/SV/";
	my $vcfSV = "$SVdir/allSV.$MFopt{SVcallerFlag}.bcf"; my $vscSVsupp  = "$SVdir/allSV-sup.$MFopt{SVcallerFlag}.bcf";
	my $CRAMmap = "$finalMapDir/$SmplName-smd.cram";
	my $SupCRAMmap = "$finalMapDir/$SmplName.sup-smd.cram";
	my $inputRawFile = "$curOutDir/input_raw.txt";
	my $binningDir = "$finalCommAssDir/Binning/";
	my @smplIDs = ("");#smpls in current assembly group
	if ($MFopt{DoMetaBat2}){ 
		@smplIDs = @{$DOs{$cAssGrp}{SmplID}};
	}
	$binningDir .= (getBinSubdirName($MFopt{DoMetaBat2})) . "/";
	#if ($MFopt{DoMetaBat2} == 3){$binningDir .= "MD/";} elsif ($MFopt{DoMetaBat2} == 2){$binningDir .= "SB/";} elsif ($MFopt{DoMetaBat2} == 1){$binningDir .= "MB2/";} else {die "Unkown binning option $MFopt{DoMetaBat2}\n";}
	my $BinningOut = "$binningDir/$smplIDs[-1]";#	$BinningOut = "$binningDir/SB/$smplIDs[-1]" if ($MFopt{DoMetaBat2} == 2);	$BinningOut = "$binningDir/MD/$smplIDs[-1]" if ($MFopt{DoMetaBat2} == 3);
	
	
	
	# Resolved checkpoints for the current sample.
	my %sampleCheckpoints = (
		primaryMapping => "$CRAMmap.sto",
		supportMapping => "$SupCRAMmap.sto",
		mappingComplete => "$finalMapDir/done.sto",
		primaryConsensus => "$contigsSNP.SNP.cons.stone",
		supportConsensus => "$contigsSNP.SNP.supp.cons.stone",
	);
	if ($MFopt{normSNPindels}) {
		$sampleCheckpoints{primaryConsensus} =~ s/\.stone/\.norm\.stone/;
		$sampleCheckpoints{supportConsensus} =~ s/\.stone/\.norm\.stone/;
	}

	# Once a sample has passed the full completion and cleanup path in this run,
	# revisit only its most informative requested outputs. Checks are deliberately
	# ordered so a missing early-stage product avoids lower-priority filesystem IO.
	if ($runOptions{loopCount} && !$loopFinalVerification
			&& $loopSampleCompleted{$JNUM}) {
		my $supportMappingRequested = $MFopt{mapSupport2Assembly}
			&& ($map{$curSmpl}{SupportReads} || '') ne '';
		my @priorityStages = (
			{
				name => 'mapping',
				required => $MFopt{map2Assembly} && $map{$curSmpl}{hasPrimaryRds},
				kind => 'exists', all => [$sampleCheckpoints{mappingComplete}],
			},
			{
				name => 'support mapping', required => $supportMappingRequested,
				kind => 'exists', all => [$sampleCheckpoints{supportMapping}],
			},
			{
				name => 'depth',
				required => $MFopt{DoAssembly} && $MFopt{map2Assembly}
					&& $map{$curSmpl}{hasPrimaryRds},
				any => [
					$coveragePerCtg, "$finalMapDir/$SmplName-smd.bam.coverage.gz",
					"$finalMapDir/$SmplName-smd.cram.coverage.gz",
				],
			},
			{
				name => 'assembly checkpoint', required => $MFopt{DoAssembly},
				kind => 'exists',
				all => ["$finalCommAssDir/$checkpointNames{assemblyDone}"],
			},
			{
				name => 'assembly outputs', required => $MFopt{DoAssembly},
				all => [$finAssLoc],
				any => [
					"$finalCommAssDir/genePred/proteins.shrtHD.faa",
					"$finalCommAssDir/genePred/proteins.shrtHD.faa.gz",
					"$finalCommAssDir/genePred/proteins.bac.shrtHD.faa",
					"$finalCommAssDir/genePred/proteins.bac.shrtHD.faa.gz",
				],
			},
			{
				name => 'binning assignment', required => $MFopt{DoMetaBat2},
				kind => 'exists', all => [$BinningOut],
			},
			{
				name => 'binning statistics', required => $MFopt{DoMetaBat2},
				all => ["$BinningOut.assStat"],
			},
			{
				name => 'CheckM', required => $MFopt{DoMetaBat2} && $MFopt{useCheckM1},
				kind => 'exists', all => ["$BinningOut.cm"],
			},
			{
				name => 'CheckM2', required => $MFopt{DoMetaBat2} && $MFopt{useCheckM2},
				all => ["$BinningOut.cm2"],
			},
			{
				name => 'SNP consensus',
				required => $MFopt{DoConsSNP} && $map{$curSmpl}{hasPrimaryRds},
				kind => 'exists', all => [$sampleCheckpoints{primaryConsensus}],
			},
			{
				name => 'SNP VCF',
				required => $MFopt{DoConsSNP} && $map{$curSmpl}{hasPrimaryRds}
					&& $MFopt{saveVCF},
				any => [$vcfSNP, "$vcfSNP.gz"],
			},
			{
				name => 'support SNP consensus',
				required => $MFopt{DoSuppConsSNP} && $supportMappingRequested,
				kind => 'exists', all => [$sampleCheckpoints{supportConsensus}],
			},
			{
				name => 'support SNP VCF',
				required => $MFopt{DoSuppConsSNP} && $supportMappingRequested
					&& $MFopt{saveVCF},
				any => [$vcfSNPsupp, "$vcfSNPsupp.gz"],
			},
			{
				name => 'structural variants',
				required => $MFopt{callSVs} && $map{$curSmpl}{hasPrimaryRds},
				all => [$vcfSV],
			},
			{
				name => 'support structural variants',
				required => $MFopt{callSVsSupp} && $supportMappingRequested,
				all => [$vscSVsupp],
			},
		);
		my $priority = priority_outputs_complete(\@priorityStages);
		if ($priority->{complete}
				&& (!$MFconfig{rmScratchTmp} || !-d $smplTmpDir)) {
			print "Sample remains complete after priority output check\n"
				unless $MFconfig{silent};
			loop2C_check($cAssGrp, \@sampleDeps);
			next;
		}
		delete $loopSampleCompleted{$JNUM};
	}
	
	# collect stats on seq qual, assembly etc
	if ($MFconfig{alwaysDoStats}){
		push @{$runReport{order}}, $SmplName unless $runReport{seen}{$SmplName}++;
		$runReport{context}{$SmplName} = {
			DIR => $dir2rd,
			input_dir => $curOutDir,
			assembly_dir => $asmDir,
		};
	}
	
	$map{$curSmpl}{inputFilesEmpty} = 0; 
	if (-e "$curOutDir/SMPL.empty"){$map{$curSmpl}{inputFilesEmpty} = 1;}

	#detect what already exists..
	my $efinAssLoc = 0; $efinAssLoc = 1  if (-s $finAssLoc && -e "$finalCommAssDir/$checkpointNames{assemblyDone}");
	#die "$efinAssLoc\n$finalCommAssDir/$checkpointNames{assemblyDone}\n";
	# Only canonical MATAFILER4 outputs participate in workflow recovery. Scratch
	# trees left by older releases are deliberately ignored.
	#activate if two assemblies for single sample required, e.g. hybrid assemblies
	#my $doPreAssmFlag = 0; my $postPreAssmblGo =0 ;
	
	
	my ($ePreAssmbly,$doPreAssmFlag,$postPreAssmblGo,$ePreAssmblPck) = prepPreAssmbl($finalCommAssDir,$metaGpreAssmblDir,$finalMapDir, "$nodeSpTmpD/preAssmblData/",
				$ContigStatsDir, $cAssGrp, $finAssLoc,$finalCommAssDir);#moves files to new locations
	
	
	#print "$ePreAssmbly,$doPreAssmFlag,$postPreAssmblGo,$ePreAssmblPck\n$finalCommAssDir/$checkpointNames{preAssemblyDone}\n$metaGpreAssmblDir/moved.sto\n";
	
	my $eCovAsssembly = 1; $eCovAsssembly = 0 if (!fileGZe($coveragePerCtg) );
	$eCovAsssembly = 0 if (!fileGZe( $markerGenesPerCtg) && $AssemblyGo);
	#die "COV: $eCovAsssembly\n";
	my $eSuppCovAsssembly = 0; $eSuppCovAsssembly = 1 if (fileGZe($suppCoveragePerCtg));
	#die "$eCovAsssembly  $eSuppCovAsssembly  $suppCoveragePerCtg\n";
	#will be created in contigstats step (not related to bowtie & sortbam)
	my $eFinMapCovGZ = 0; $eFinMapCovGZ = 1 if (-e $sampleCheckpoints{primaryMapping} && -e "$finalMapDir/$SmplName-smd.bam.coverage.gz");#"$finalMapDir/$SmplName-smd.bam.coverage.gz";
	#$MFopt{mapSupport2Assembly}
	my $locMapSup2Assembly =0; $locMapSup2Assembly =1 if ($MFopt{mapSupport2Assembly} && $map{$curSmpl}{"SupportReads"} ne "");
	my $eFinSupMapCovGZ = 0; $eFinSupMapCovGZ = 1 if ($locMapSup2Assembly && -e $sampleCheckpoints{supportMapping} && fileGZe("$finalMapDir/$SmplName.sup-smd.bam.coverage"));
	#die "$locMapSup2Assembly $eFinSupMapCovGZ  $MFopt{mapSupport2Assembly}  $finalMapDir\n$sampleCheckpoints{supportMapping}\n";
	my $dfinalCommAssDir = 0 ; $dfinalCommAssDir = 1 if (-d $finalCommAssDir);
	my $eFinalMapDir = 0; $eFinalMapDir = 1 if (-s $sampleCheckpoints{mappingComplete});
	#upload2EBI 
	my $DoUploadRawReads = 0; $DoUploadRawReads = 1 if ($MFconfig{uploadRawRds} ne ""); 

	
	
	#die "$eFinMapCovGZ $sampleCheckpoints{primaryMapping} && -e $finalMapDir/$SmplName-smd.bam.coverage.gz\n";
	#die "$metagAssDir\n";
	

	
	
	my $locRedoSNPcalling =0; 
	my $locRedoSVs = 0;
	my $locRedoAssMapping = $MFopt{redoAssMapping};
	# A hybrid preassembly and the final hybrid assembly intentionally use the
	# same canonical mapping directory.  Stones and coverage files alone cannot
	# tell which assembly an existing CRAM targets.  New mappings carry the size
	# and mtime of their reference; once the final assembly exists, reject a
	# missing (legacy) or non-matching stamp so preassembly CRAMs are remapped.
	if ($MFopt{DoAssembly} == 5 && -s $finAssLoc
			&& -e "$finalCommAssDir/$checkpointNames{assemblyDone}") {
		my @mappingReferenceChecks;
		push @mappingReferenceChecks, [
			"$finalMapDir/$SmplName-smd.reference.stat", 'primary'
		] if ($eFinMapCovGZ);
		push @mappingReferenceChecks, [
			"$finalMapDir/$SmplName.sup-smd.reference.stat", 'supplementary'
		] if ($eFinSupMapCovGZ);
		for my $referenceCheck (@mappingReferenceChecks) {
			my ($stamp, $kind) = @{$referenceCheck};
			next if mapping_reference_matches($stamp, $finAssLoc);
			print "Hybrid $kind mapping does not identify the current final assembly; remapping $SmplName\n";
			$locRedoAssMapping = 1;
			last;
		}
	}

	#check if current assembly group is the same as before!
	my $locRewrite = 0; my $locRedoAssembl = 0;
	if ($efinAssLoc && -e "$finalCommAssDir/smpls_used.txt"){
		my @actualMembers;
		open my $memberFH, '<', "$finalCommAssDir/smpls_used.txt"
			or die "Could not read assembly-group membership: $!\n";
		while (my $member = <$memberFH>){
			chomp $member;
			push @actualMembers, $member unless ($member =~ m/^\s*$/);
		}
		close $memberFH;
		my $expectedMembers = assembly_group_output_dirs(\%map, $cAssGrp);
		if (!workflow_members_match($expectedMembers, \@actualMembers)){
			print "$cAssGrp assembly-group membership has changed! (previous: "
				.scalar(@actualMembers).", current: ".scalar(@{$expectedMembers}).")\n"
				."$finalCommAssDir\nRemoving assembly and all processed reads\n";
			unless ($MFconfig{OKtoRWassGrps}) {print "Stopping MATAFILER, human intervention needed.. use the flag \"-OKtoRWassGrps 1\" to allow MATAFILER to delete files\n"; die;}
			$locRewrite=1;
		}
	}
	
	if ($eFinMapCovGZ && (!$ePreAssmblPck && !$ePreAssmbly && !$efinAssLoc ) ){#impossible, so reason must be severe!
		print "Mapping exists, but no assembly, removing mapping..\n$finalMapDir\n$finAssLoc\n$metaGassembly\n";
		die unless ($MFconfig{OKtoRWassGrps});
		$locRewrite = 0; $locRedoAssembl = 0;
	}
	
	if ($efinAssLoc && $finAssLoc ne "$finalCommAssDirSingle/scaffolds.fasta.filt" && -s "$finalCommAssDirSingle/scaffolds.fasta.filt"){
		print "Something wrong.. assembly group assembly and single assembly present:\n$finAssLoc\n$finalCommAssDirSingle/scaffolds.fasta.filt\n";
		$locRedoAssMapping=1;$locRedoSNPcalling=1;$locRedoSVs=1;
		system "rm -fr $finalCommAssDirSingle; mkdir -p $finalCommAssDirSingle;\n";
		$eCovAsssembly=0;$eFinMapCovGZ=0;$eFinalMapDir=0;$eFinSupMapCovGZ=0;$eSuppCovAsssembly=0;$eSuppCovAsssembly=0;
	}
	
	# Full sample statistics are collected after the submission loop. They are
	# reporting data and must not force every sample to open all historical logs
	# before the fast completion path can run.
	if ( ($MFconfig{skipWrongPairedSmpls} || $MFconfig{OKtoRWassGrps}) && -e "$logDir/sdmReadCleaner.sh.etxt" && `tail -n 70 $logDir/sdmReadCleaner.sh.etxt | grep 'invalid paired read' ` ne ""){
		print "$logDir/sdmReadCleaner.sh.etxt problems! Delete outdir\n";
		if ($MFconfig{OKtoRWassGrps}){
			$locRewrite=1 ;
		}elsif ($MFconfig{skipWrongPairedSmpls}){
			loop2C_check($cAssGrp,\@sampleDeps);next;
		}
	}
	
	
	#--------------------------  DELETION SECTION  -----------------------------------
	#DELETION SECTION
	#redo run - or parts thereof	
		
	if ($MFopt{rewriteGenePred}){
		print "Deleting gene Predictions and dependent files..\n";
		system "rm -rf $finalCommAssDir/genePred" if (-d "$finalCommAssDir/genePred");
		system "rm -rf $finalCommAssDir/ContigStats" if (-d "$finalCommAssDir/ContigStats");
		system "rm -rf $ContigStatsDir/*pergene* $ContigStatsDir/GeneStats.tx* $ContigStatsDir/GTDBmg $ContigStatsDir/FMG";
	}
	if ($MFconfig{OKtoRWassGrps} && $locRewrite){
		print "Deleting previous results.. rerun MATAFILER for sample\n";
		system ("rm -r -f $asmDir $finalCommAssDir");
		system("rm -f -r $curOutDir $smplTmpDir $MFglobal{collectFinished} ");
		#next; #too deep, needs a complete new round over dir..
		#$efinAssLoc = 0;	$eFinMapCovGZ = 0;	
		$efinAssLoc =0 ;  $dfinalCommAssDir =0;
		$eCovAsssembly = 0; $eSuppCovAsssembly=0; $eFinSupMapCovGZ=0; $eFinMapCovGZ = 0;$eFinalMapDir = 0;$locRewrite = 0; $locRedoAssembl = 0;$eSuppCovAsssembly=0;
	} 
	if ($eCovAsssembly && -e "$ContigStatsDir/Coverage.percontig" && -e "$ContigStatsDir/Coverage.percontig.gz"){
		print "redoing covverage calculations..\n";
		$eCovAsssembly = 0; system "rm $ContigStatsDir/Coverage*";
	}

	
	#automatically delete mapping, if assembly no longer exists..
	#print "locRedoAssMapping : $locRedoAssMapping\n";
	if ($MFopt{map2Assembly} ){
		if ($eFinMapCovGZ && !$eFinalMapDir){$locRedoAssMapping = 1 ; print "R0 ";}
		my $mappingArtifactsPresent = $eFinMapCovGZ || $eFinalMapDir || -e $sampleCheckpoints{primaryMapping}
			|| -e $CRAMmap || fileGZe("$finalMapDir/$SmplName-smd.bam.coverage.gz");
		if (!$efinAssLoc && !$ePreAssmbly && $mappingArtifactsPresent){$locRedoAssMapping = 1 ;}
		if (!$MappingGo && $map{$curSmpl}{hasPrimaryRds} && $eFinalMapDir){$locRedoAssMapping = 1 ;print "R2 ";}
		if (-e $sampleCheckpoints{primaryMapping} && (!fileGZe( "$finalMapDir/$SmplName-smd.bam.coverage.gz") || !-e $CRAMmap) ){$locRedoAssMapping = 1 ;print "R3 ";}
		#if ($eFinMapCovGZ && (exists($locStats{uniqAlign}) && $locStats{uniqAlign} > 20) && -s $CRAMmap <300){$locRedoAssMapping = 1 ;print "R4";}
		#print "$CRAMmap :: $locRedoAssMapping\n";
		if ($locRedoAssMapping){# && -e $CRAMmap){
			print "redo assem mapping!" . " -s $CRAMmap \n" ;
			#die;
		}
		
		#die "$sampleCheckpoints{primaryMapping} && !-e $finalMapDir/$SmplName-smd.bam.coverage.gz";
		
	}

	#die "locRedoAssMapping : $locRedoAssMapping $finalMapDir   : !$efinAssLoc \n" if ($locRedoAssMapping);
	#my $sizemap = -s $CRAMmap;#print "size map: " . $sizemap . "\n";
	#delete assembly
	if ($MFopt{redoAssembly} || $locRedoAssembl){
		if ($AsGrps{$cAssGrp}{CntAimAss} > 1 && !$MFconfig{OKtoRWassGrps}){
			die "Refusing to rebuild shared assembly group $cAssGrp without -OKtoRWassGrps 1\n";
		}
		print "Removing assembly ... \n" if (-e $metaGassembly );
		system "rm -fr $finalCommAssDir";
		$efinAssLoc = 0;	
		$locRedoAssMapping=1;
	}
	#delete mapping to assembly
	if ($locRedoAssMapping){
		#die "locDel\n";
		my @cramStat = stat($CRAMmap);
		my $fileSiz = @cramStat && $cramStat[7] > 0 ? $cramStat[7] : "NA";
		print "Deleting previous assembly mapping, size map: ". $fileSiz . " ; $CRAMmap\n" if ($MappingGo && $eFinalMapDir);
		#die "$finAssLoc && !-e $metaGassembly\n";
		system "rm -fr $finalMapDir $ContigStatsDir/Coverage.* $ContigStatsDir/Cov.sup.*";
		$eFinMapCovGZ = 0;	
		$eCovAsssembly = 0; $eSuppCovAsssembly=0; $eFinSupMapCovGZ=0; $eFinalMapDir = 0;
		#are there SNPs called? remove as well..
		$locRedoSNPcalling=1; $locRedoSVs=1;
	}
	#Case: primary assembly mapping was done, support reads were not yet mapped.. need to redo binning 
	#print "$locMapSup2Assembly && !$eFinSupMapCovGZ) && ($MFopt{map2Assembly} && $eFinMapCovGZ \n";
	my @binningBaseStat = stat($BinningOut);
	my $binningArtifactsPresent = @binningBaseStat || -e "$BinningOut.cm"
		|| -e "$BinningOut.cm2" || -e "$BinningOut.assStat";
	if ( ($locMapSup2Assembly && !$eFinSupMapCovGZ) && ($MFopt{map2Assembly} && $eFinMapCovGZ )
			&& $binningArtifactsPresent ){
		print "redoing binning due to support mapping not included..\n";
		system("rm -rf  $binningDir/");
	}
	#debug case: binning was empty
	if ($MFopt{DoMetaBat2} && ( $MFopt{BinnerRedoAll}
			|| ($MFopt{BinnerRedoEmpty} && @binningBaseStat && $binningBaseStat[7] == 0) ) ){
		print "redoing binning due to empty bins (flag -redoEmptyBins 1) ..\n";
		system "rm -rf $binningDir";
	}
	if ( (!$doPreAssmFlag || !$ePreAssmblPck) && !$map{$curSmpl}{inputFilesEmpty} &&
				((!$eFinMapCovGZ && $eCovAsssembly) || ($eSuppCovAsssembly && !$eFinSupMapCovGZ) ) #redo only contigstats related to coverage..
				 || $MFconfig{redoCS}){
		#print "redoing contig stats global..\n";
		#die "(!$doPreAssmFlag || !$ePreAssmblPck) && ((!$eFinMapCovGZ && $eCovAsssembly) || ($eSuppCovAsssembly && !$eFinSupMapCovGZ) )\n";
		system("rm -rf $finalCommAssDir/ContigStats/ $ContigStatsDir $binningDir/");
		$eCovAsssembly = 0; $eSuppCovAsssembly=0; #contigstats needs redoing..
	}
	my $KrakenOD = $curOutDir."Tax/kraken/$MFopt{globalKraTaxkDB}/";
	if ($MFopt{RedoKraken} && -d $KrakenOD) {system "rm -r $KrakenOD" ;}
	if ($MFopt{RedoRiboFind}){system "rm -rf $curOutDir/ribos";}
	if ($MFopt{RedoRiboAssign}){system "rm -rf $curOutDir/ribos//ltsLCA";}
	if ($MFopt{DoRibofind} && -e "$curOutDir/LOGandSUB/RiboLCA.sh.etxt"){
		#my $LCAetxt = `cat $curOutDir/LOGandSUB/RiboLCA.sh.etxt`;
		#if ($LCAetxt =~ m/ParseError thrown: Unexpected character .\@. found/){system "rm -rf $curOutDir/ribos";}
	}
	if ($locRedoSNPcalling){system "rm -fr $SNPdir";}
	
	if ($locRedoSVs){system "rm -fr $SVdir";}
	if ($MFopt{redoSNPcons}){		system "rm -rf $SNPdir $genePredSNP* $contigsSNP* $genePredAASNP* $logDir/SNP";
	} elsif ($MFopt{redoSNPgene}){		system "rm -rf $genePredSNP* $genePredAASNP* ";
	}
	my $boolGenePredOK=0;
	if ($MFopt{DoEukGenePred}){
		$boolGenePredOK = 1 if ( fileGZe("$finalCommAssDir/genePred/proteins.bac.shrtHD.faa") || ($MFopt{pseudoAssembly} && fileGZe("$finalCommAssDir/genePred/proteins.bac.shrtHD.faa")));
	} else {
		$boolGenePredOK = 1 if (fileGZe("$finalCommAssDir/genePred/proteins.shrtHD.faa") || ($MFopt{pseudoAssembly} && fileGZe("$finalCommAssDir/genePred/proteins.shrtHD.faa")) );
	}
	#DEBUG to gzip outputs..
	if ($MFopt{genePredGZenforce} && $boolGenePredOK && -e "$finalCommAssDir/genePred/genes.gff"){$boolGenePredOK =0;}
	#die "$boolGenePredOK\n$finalCommAssDir\n";

	
	#central flag that decides if an assembly is done
	my $boolAssemblyOK=0;
	$boolAssemblyOK=1 if ($boolGenePredOK && $efinAssLoc );#&& (!$MFopt{map2Assembly} || $eFinMapCovGZ ) );
	my $assemblyOutputsRequired = $MFopt{DoAssembly} ? 1 : 0;
	my $assemblyWorkflowComplete = !$assemblyOutputsRequired || $boolAssemblyOK;
	#die "$boolGenePredOK && $efinAssLoc && (!$MFopt{map2Assembly} || $eFinMapCovGZ ) $locRedoAssMapping\n";
			#&& (-s "$finalMapDir/$SmplName-smd.bam" || -s "$finalMapDir/$SmplName-smd.cram")
	#die "$boolAssemblyOK\n$finalCommAssDir/genePred/proteins.shrtHD.faa\n$finalMapDir/$SmplName-smd.bam.coverage.gz\n";
	

	if ( ( !$boolAssemblyOK && $MFconfig{unfiniRew}==1 ) ){
		print "Deleting unfinished previous results; this sample will be rebuilt on the next pipeline pass.\n";
		system('rm', '-rf', '--', $curOutDir, $smplTmpDir, $MFglobal{collectFinished}) == 0
			or die "Failed to remove unfinished results for $curSmpl\n";
		loop2C_check($cAssGrp,\@sampleDeps);next;
	}

	
#--------------------- secondary map deletions & flags --------------------------
	my $boolScndMappingOK = 1; my $iix =0;
	my $boolScndCoverageOK = 1;
	if ($MFopt{MapRewrite2nd}){ 
		print "rewriting secondary map\n";
		foreach my $bwt2outDTT (@bwt2outD){
			my $expectedMapCovGZ = "$bwt2outDTT/$bwt2ndMapNmds[$iix]"."_".$SmplName."-0-smd.bam.coverage.gz"; #$bamcramMap : 2nd map only has .bam output
			system "rm -f $expectedMapCovGZ*";
			$eFinMapCovGZ = 0;	
		}
	}
	#die "eFinMapCovGZ $eFinMapCovGZ\n";
	foreach my $bwt2outDTT (@bwt2outD){
		my $expectedMapCovGZ = "$bwt2outDTT/$bwt2ndMapNmds[$iix]"."_".$SmplName."-0-smd.bam.coverage.gz";
		my $expectedMapBam = "$bwt2outDTT/$bwt2ndMapNmds[$iix]"."_".$SmplName."-0-smd.bam";
		$iix++;
		
		#print $expectedMapCovGZ."\n";
		if ( -e "$expectedMapCovGZ" && -e $expectedMapBam && $MappingGo  ){
			$boolScndMappingOK=1;
		}else{
			$boolScndMappingOK=0; $boolScndCoverageOK=0;
			#be clean
			system "rm -f $expectedMapBam*";
			last;
		}
		if ($MFopt{mapModeCovDo} && (!-e $expectedMapCovGZ.".median.percontig" || !-e $expectedMapCovGZ.".percontig"|| !-e $expectedMapCovGZ.".pergene")){
			$boolScndCoverageOK=0;
			system "rm -f $expectedMapCovGZ.*";
		}
	}
	if (@bwt2outD == 0 ){$boolScndMappingOK = 1 ; $boolScndCoverageOK=1;}#|| !$MappingGo);	
	#print $boolScndMappingOK."\n$boolAssemblyOK\n";

#--------------------- other flags --------------------------
	#contamination flag: redo/do contaminant removal to eg make sure human contamination is correctly logged
	my $calcContamination = 0;
	if ($efinAssLoc && $MFopt{completeContaStats}) {
		my $contamination = getContamination(
			"$curOutDir/LOGandSUB/KrakHS.sh.etxt",
			"$curOutDir/LOGandSUB/KrakHS.sh.otxt", '',
		);
		$calcContamination = 1
			if (($contamination->{FilteredContaRdsPerc} || '') eq "?\t");
	}
	#print "Conta: $calcContamination   \"$locStats{contamination}\"\n";

	
	
#	#-----------------------  FLAGS  ------------------------  
	
	
	#check on processes not dependent on assemblies
	prepareDiamondRerun($curOutDir);
	my ($calcKraken,$calcDiamond,$calcDiaParse,$calcRibofind,$calcRiboAssign,$calcGenoSize,
			$calcMetaPhlan, $calcMOTU2,$calcTaxaTar) = checkRawProgsFin($curOutDir,$SmplName);
	publishKrakenResults($curOutDir,$SmplName) if ($MFopt{DoKraken} && !$calcKraken);
	#not complete yet? Then delete..
	if ($MFconfig{redoFails} && ($calcRibofind||$calcDiamond || $calcDiaParse ||$calcMOTU2 || $calcMetaPhlan || $calcTaxaTar)){
		print "Removing failed results for $curSmpl; this sample will be rebuilt on the next pipeline pass.\n";
		my @failedSampleTargets = ($curOutDir, $smplTmpDir, $MFglobal{collectFinished});
		if ($AsGrps{$cAssGrp}{CntAimAss} <= 1){
			push @failedSampleTargets, $asmDir, $finalCommAssDir;
		} else {
			print "Retaining shared assembly-group outputs; use -redoAssembly with "
				."-OKtoRWassGrps 1 for an explicitly authorized group rebuild.\n";
		}
		system('rm', '-rf', '--', @failedSampleTargets) == 0
			or die "Failed to remove failed results for $curSmpl\n";
		loop2C_check($cAssGrp,\@sampleDeps);next;
	}



	# Support coverage belongs to the final hybrid assembly. During preassembly it
	# is intentionally absent and must not keep resubmitting no-op ContigStats jobs.
	my $supportCoverageRequired = $locMapSup2Assembly && $efinAssLoc
		&& !$doPreAssmFlag && !$ePreAssmblPck;
	my $allMapDone =0;#used for SNP calling and Binning - but binning requires info if all maps are finished from all samples
	$allMapDone = 1 if ( 
				(!$map{$curSmpl}{hasPrimaryRds} || ($eFinMapCovGZ && -e "$finalMapDir/$SmplName-smd.$bamcramMap" && $eCovAsssembly )) #primary
				&& ($eSuppCovAsssembly || !$supportCoverageRequired)  #secondary
				&& $AsGrps{$cAssGrp}{MapDeps} !~ m/[^;]/ );
	#die "$allMapDone\n-e $finalMapDir/$SmplName-smd.$bamcramMap && $eCovAsssembly && !$ePreAssmbly && ($eSuppCovAsssembly || !$locMapSup2Assembly) && $AsGrps{$cAssGrp}{MapDeps} !~ m/[^;]/\n";
	
	#coverage done?
	#my $allCovDone = 0; $allCovDone = 1 if ( ($eSuppCovAsssembly || !$locMapSup2Assembly) && ($eCovAsssembly || !$map{$curSmpl}{hasPrimaryRds}) );
	my $calcCoverage = 0; $calcCoverage =1 if ((($map{$curSmpl}{hasPrimaryRds} && !$eCovAsssembly) || (!$eSuppCovAsssembly && $supportCoverageRequired) ) && $MFopt{map2Assembly});
	#print "$calcCoverage = 1 if (($eSuppCovAsssembly || !$locMapSup2Assembly) && ($eCovAsssembly || !$MappingGo) );\n";
	
	#binning done?
	my $calcBinning = 0;
	my $binningComplete = binningOutputsComplete(
		$BinningOut, $MFopt{useCheckM1}, $MFopt{useCheckM2},
	);
	# A downstream job need not wait for its inputs to exist at submission time.
	# For a normal (non-hybrid-preassembly) group, AssemblyGo means that this
	# invocation will schedule the final assembly and all of its mappings.  The
	# binner is submitted later with those scheduler dependencies.
	my $supportMappingPublished = !$locMapSup2Assembly || $eFinSupMapCovGZ;
	if ($MFopt{DoMetaBat2} && !$doPreAssmFlag && !$ePreAssmblPck && $AssemblyGo
			&& $supportMappingPublished
			&& !$binningComplete) {
		$calcBinning=$MFopt{DoMetaBat2};
		#die "$MFopt{DoMetaBat2} && $boolAssemblyOK && $AssemblyGo && $AsGrps{$cAssGrp}{MapDeps} !~ m/[^;]/ &&  (!-e $BinningOut.cm || !-s $BinningOut.cm2\n";
	}
	

	#die "$metaGassembly\n$finAssLoc\n$nodeSpTmpD\n$eCovAsssembly\n";

#	#and some more flags for subprocesses
	my $nonPareilFlag = !-s "$nonParDir/$SmplName.npo" && $MFopt{DoNonPareil} ;
	my $scaffoldFlag = 0; if ( !-e $STOfinScaff && $map{$curSmpl}{"SupportReads"} =~ m/mate/i ){$scaffoldFlag = 1 ;}# print "SUPP:: $map{$curSmpl}{SupportReads}\n";}
	my $assemblyFlag = 0; $assemblyFlag = 1 if ( $MFopt{DoAssembly} && !$boolAssemblyOK && !$efinAssLoc );

	#die "$assemblyFlag = 1 if ( $MFopt{DoAssembly} && !$boolAssemblyOK && !$efinAssLoc && !-e $metaGassembly\n $doPreAssmFlag\n";
	my $calcReadMerge = 0;
	$calcReadMerge = 1 if ($MFopt{doReadMerge} && ($MFopt{calcOrthoPlacement} || $calcDiamond || $calcGenoSize));
	my $mapAssFlag = 0; $mapAssFlag = 1 if ($map{$curSmpl}{hasPrimaryRds} && $MFopt{map2Assembly} && !$eFinMapCovGZ  );
	#only for support reads (from hybrid assemblies)
	my $mapSuppAssFlag =0;$mapSuppAssFlag = 1 if ($supportCoverageRequired && !$eFinSupMapCovGZ);#hasSuppRds(\%AsGrps,$cAssGrp,$curSmpl ) );
	my $calcSuppCoverage = 0; $calcSuppCoverage =1 if ($MFopt{mapSupport2Assembly} && !$eSuppCovAsssembly && $map{$curSmpl}{"SupportReads"} ne "" && $mapSuppAssFlag && !$doPreAssmFlag);
	
	#die "$calcCoverage\n$mapAssFlag , $mapSuppAssFlag :: $ePreAssmbly && $doPreAssmFlag XX $postPreAssmblGo,$ePreAssmblPck\n";
	
	#die "$mapSuppAssFlag = 1 if ($locMapSup2Assembly && !$eFinSupMapCovGZ && $efinAssLoc \n";

	#requires only bam/cram && assembly
	my $calcConsSNP=0; 
	if ($MFopt{DoConsSNP} && $map{$curSmpl}{hasPrimaryRds} && !$doPreAssmFlag && !$ePreAssmblPck){
		my $exSNPf= fileGZe($vcfSNP);
		if ($exSNPf && fileGZs($vcfSNP) == 0){system "rm -f $vcfSNP*";$exSNPf=0;} #some old versions produced an empty vcf file..
		$calcConsSNP=0; $calcConsSNP =1 if ( !-e $sampleCheckpoints{primaryConsensus} || ($MFopt{saveConsFastas} &&  fileGZe($genePredSNP)==0  ) || ($MFopt{saveVCF} &&  !$exSNPf ) ) ;
	}
	my $calcSuppConsSNP=0; $calcSuppConsSNP =1 if (!$doPreAssmFlag && !$ePreAssmblPck && $locMapSup2Assembly && $MFopt{DoSuppConsSNP} && (!-e  $sampleCheckpoints{supportConsensus}  ));
	
	
	
	#structural variants calcs
	my $calcSVs = 0; $calcSVs = 1 if ( $map{$curSmpl}{hasPrimaryRds} && $MFopt{callSVs} && !-e $vcfSV && !$doPreAssmFlag && !$ePreAssmblPck);
	my $calcSVsSupp = 0; $calcSVsSupp = 1 if ($locMapSup2Assembly && !$doPreAssmFlag && !$ePreAssmblPck && $MFopt{callSVsSupp} && !-e $vscSVsupp );
	if (!$mapAssFlag && ($calcConsSNP || $calcSVsSupp || $calcSVs || $calcSuppConsSNP) && $eFinMapCovGZ && !$allMapDone ){
		#die $mapAssFlag;
		$mapAssFlag = 1; #reactivate mapping of assembly, to allow SNP consensus calling..
		$eFinMapCovGZ = 0;
	}


	
	my $calc2ndMapSNP = 0; $calc2ndMapSNP = 1 if ($MFopt{Do2ndMapSNP});
	my $pseudAssFlag = 0; $pseudAssFlag = 1 if ($MFopt{pseudoAssembly} && $map{$curSmpl}{ExcludeAssem} eq "0" && (!-e $pseudoAssFileFinal.".sto" || !$boolGenePredOK));
	my $dowstreamAnalysisFlag = 0; 
	#$unpackZip simulates dowstreamAnalysisFlag, just to get sdm running..
	$dowstreamAnalysisFlag=1 if ( $MFconfig{unpackZip} || $calcContamination || $MFopt{calcOrthoPlacement} || $scaffTarExternal ne "" || $assemblyFlag  
		|| $pseudAssFlag || $scaffoldFlag  || $nonPareilFlag || $calcGenoSize || $calcDiamond || $calcDiaParse
		|| $MFopt{DoCalcD2s} || $calcKraken || $calcRibofind || $calcRiboAssign || $calcMOTU2 ||  $calcMetaPhlan || $calcTaxaTar );
	my $requireRawReadsFlag = 0;#only list modules that really need raw reads
	$requireRawReadsFlag = 1 if ( !$boolScndMappingOK || $calcContamination);
	
	my $primaryCleanPending = !-e "$smplTmpDir/seqClean/filterDone.stone";
	my $supportCleanPending = ($map{$curSmpl}{SupportReads} || '') ne ''
		&& !-e "$smplTmpDir/seqClean/filterSupplDone.stone";
	my $seqCleanFlag = $dowstreamAnalysisFlag && ($primaryCleanPending || $supportCleanPending) ? 1 : 0;
	my $porechopFlag = 0;
	$porechopFlag = 1 if ($MFopt{usePorechop} && $dowstreamAnalysisFlag && !-e "$smplTmpDir/rawRds/poreChopped.stone");
	#die "$assemblyFlag\t$seqCleanFlag\t$boolScndMappingOK\n";
	my $calcUnzip=0;
	$calcUnzip=1 if ($calcDiamond || $porechopFlag || $seqCleanFlag  || $mapAssFlag || $mapSuppAssFlag || (!$MFopt{useUnmapped} && !$boolScndMappingOK) || $MFconfig{uploadRawRds} ne ""); 
	#print "chk1 $mapSuppAssFlag $calcSuppCoverage $eSuppCovAsssembly\n" ;

	# Cleanup is destructive to read-cleaning stones, mapping indexes and staged
	# reads. A generic list of jobs submitted in this pass is not proof that the
	# terminal assembly analyses are complete; track each producer and terminal
	# analysis explicitly, including same-pass ConsSNP dependencies.
	my $cleanupContigSubparts = $MFconfig{defaultContigSubs}."gFG";
	$cleanupContigSubparts .= "m" if ($MFopt{DoBinning});
	$cleanupContigSubparts .= "k" if ($MFopt{kmerAssembly});
	$cleanupContigSubparts .= "4" if ($MFopt{kmerPerGene});
	my $cleanupContigStatsComplete = contigStatsOutputsComplete(
		$curOutDir, $finalCommAssDir, $cleanupContigSubparts, 1,
		$curSmpl, $supportCoverageRequired,
	);
	my $cleanupPrimaryConsensusRequired = $MFopt{DoConsSNP} && $map{$curSmpl}{hasPrimaryRds};
	my $cleanupSupportConsensusRequired = $MFopt{DoSuppConsSNP} && $locMapSup2Assembly;
	my $cleanupConsensusRequired = $cleanupPrimaryConsensusRequired || $cleanupSupportConsensusRequired;
	my $cleanupVariantRequired =
		$cleanupConsensusRequired
		|| ($MFopt{callSVs} && $map{$curSmpl}{hasPrimaryRds})
		|| ($MFopt{callSVsSupp} && $locMapSup2Assembly);
	my $cleanupVariantsComplete = !$cleanupVariantRequired
		|| (!$calcConsSNP && !$calcSuppConsSNP && !$calcSVs && !$calcSVsSupp);
	my $terminalOutputsComplete = !$assemblyOutputsRequired || (
		!$doPreAssmFlag && !$ePreAssmblPck
		&& $efinAssLoc && $cleanupContigStatsComplete
		&& (!$MFopt{DoMetaBat2} || $binningComplete)
		&& $cleanupVariantsComplete
	);
	my $cleanupRequirements = cleanupCompletionRequirements(
		contig_dir => $ContigStatsDir,
		assembly_dir => $finalCommAssDir,
		contig_subparts => $assemblyOutputsRequired ? $cleanupContigSubparts : '',
		primary_coverage_required => $assemblyOutputsRequired && $map{$curSmpl}{hasPrimaryRds},
		support_coverage_required => $assemblyOutputsRequired && $supportCoverageRequired,
		binning_base => ($assemblyOutputsRequired && $MFopt{DoMetaBat2}) ? $BinningOut : '',
		primary_snp_stone => ($assemblyOutputsRequired && $cleanupPrimaryConsensusRequired) ? $sampleCheckpoints{primaryConsensus} : '',
		support_snp_stone => ($assemblyOutputsRequired && $cleanupSupportConsensusRequired) ? $sampleCheckpoints{supportConsensus} : '',
		consensus_contigs => ($assemblyOutputsRequired && $cleanupConsensusRequired && $MFopt{saveConsFastas}) ? $contigsSNP : '',
		consensus_genes => ($assemblyOutputsRequired && $cleanupConsensusRequired && $MFopt{saveConsFastas}) ? [$genePredSNP, $genePredAASNP] : [],
		primary_vcf => ($assemblyOutputsRequired && $cleanupPrimaryConsensusRequired && $MFopt{saveVCF}) ? $vcfSNP : '',
		support_vcf => ($assemblyOutputsRequired && $cleanupSupportConsensusRequired && $MFopt{saveVCF}) ? $vcfSNPsupp : '',
		primary_sv => ($assemblyOutputsRequired && $MFopt{callSVs} && $map{$curSmpl}{hasPrimaryRds}) ? $vcfSV : '',
		support_sv => ($assemblyOutputsRequired && $MFopt{callSVsSupp} && $locMapSup2Assembly) ? $vscSVsupp : '',
	);

	if ($scaffTarExternal ne "" &&  $map{$curSmpl}{"SupportReads"} !~ /mate/i && $scaffTarExtLibTar ne $curSmpl ){print"scNxt\n";loop2C_check($cAssGrp,\@sampleDeps);next;}
	

#	#-----------------------  END FLAGS  ------------------------  

	#some more flow control..
	if ( !$DoUploadRawReads && $boolScndMappingOK && !$MFopt{DoCalcD2s} &&
		!$calcConsSNP && !$calcSuppConsSNP && !$calcSVs && !$calcSVsSupp &&
		!$calcBinning && !$calc2ndMapSNP && $assemblyWorkflowComplete && $boolScndCoverageOK
		 && !$calcCoverage && !$calcSuppCoverage && !$dowstreamAnalysisFlag
		 && $terminalOutputsComplete
		#&& !$calcRibofind && !$calcRiboAssign && !$MFopt{calcOrthoPlacement} && !$calcGenoSize && !$calcDiamond && !$calcDiaParse && 
		#!$calcMetaPhlan && !$calcTaxaTar && !$calcMOTU2 && !$calcKraken && $scaffTarExternal eq ""
	){
		if ( ($boolAssemblyOK || ($doPreAssmFlag && $ePreAssmbly && !$ePreAssmblPck)) && !$locRedoAssMapping ){ #causes a lot of overhead but mainly to avoid unpacking reads again..
			$runReport{present_assemblies}++;#= $AsGrps{$cAssGrp}{CntAimAss};
		}
		my $cleanupComplete = runFinishedCleanup(finishedCleanupArguments(
			$curSmpl, $SmplName, $finalCommAssDir, $finalMapDir,
			$smplTmpDir, $finAssLoc, $logDir, $cleanupRequirements,
			$assemblyOutputsRequired,
		));
		if ($cleanupComplete && (!$MFconfig{rmScratchTmp} || !-d $smplTmpDir)) {
			$loopSampleCompleted{$JNUM} = 1;
		} else {
			delete $loopSampleCompleted{$JNUM};
		}
		print "Sample already complete; no jobs submitted\n" unless $MFconfig{silent};
		MFnext($smplLockF,\@sampleDeps,$JNUM ,$QSBoptHR); 
		loop2C_check($cAssGrp,\@sampleDeps);
		next;
	}

	# Scratch and submission directories are needed only beyond the completed
	# sample shortcut. Do not create and immediately remove them on no-op passes.
	make_path($logDir) unless -d $logDir;
	make_path($smplTmpDir) unless -d $smplTmpDir;
	my @checkLocs = ($smplTmpDir);
	push(@checkLocs,$curDir) if ($curDir ne "");
	$QSBoptHR->{LocationCheckStrg} = checkDrives(\@checkLocs);
	


	#report for debugging:
	#print " !$DoUploadRawReads && $boolScndMappingOK && !$MFopt{DoCalcD2s} &&
	#	!$calcConsSNP && !$calcSuppConsSNP && !$calcSVs && !$calcSVsSupp &&
	#	!$calcBinning && !$calc2ndMapSNP && $boolAssemblyOK && $boolScndCoverageOK 
	#	 && !$calcCoverage && !$calcSuppCoverage && !$dowstreamAnalysisFlag\n";
	 #die;
	








#-----------------------------------------------------------------------------------------
#                   actual job submissions starts from here
#-----------------------------------------------------------------------------------------

#unzipping & sorting of files into scratch dir..
#also does either trimomatic or porechop..
#die "calcUnzip $calcUnzip\n" ;#if ($calcUnzip);
#	my ($jdep,$cfp1ar,$cfp2ar,$cfpsar,$WT,$rawFiles, $mmpuNum, $libInfoRef, $inputRawSize) = 
	#my %seqSet = (pa1 => \@pa1, pa2 => \@pa2, pas => \@pas, paX1 => \@paX1, paX2 => \@paX2, paXs => \@paXs, libInfo => \@libInfo, libInfoX => \@libInfoX,
	#		totalInputSizeMB => $totalInputSizeMB,rawReads => $rawReads,mmpu => $mmpu, WT => $WT);
	#unzip and change ifastap & cfp1/cfp2
	my $curUnzipDep = ""; 
	if ($MFconfig{maxUnzpJobs} > 0 && @unzipjobs >= $MFconfig{maxUnzpJobs}){#only run X jobs in parallel, lest the cluster IO breaks down..
		$curUnzipDep = $unzipjobs[-($MFconfig{maxUnzpJobs})];#join(";",@last_n);
	}
	
	if ($map{$curSmpl}{SRA_download} ne "" || $map{$curSmpl}{ENA_download} ne "") {
		 #-> do downloads to tmp dir
		 #function that downloads to tmp dir
	}

	
	my ($jdep) =  #,$hrefSeqSet
			seedUnzip2tmp($curDir,$curSmpl,$curUnzipDep,$nodeSpTmpD,
			$smplTmpDir,$calcUnzip,$finalMapDir,
			$porechopFlag,$inputRawFile);
	#my %seqSet = %{$hrefSeqSet};
	#print "$seqSet{pa1}   $seqSet{seqTech}   $seqSet{seqTechX}\n";
	push (@unzipjobs,$jdep) unless ($jdep eq "");
	
	if ($map{$curSmpl}{inputFilesEmpty} && exists($map{$curSmpl}{inputFileSizeMB}) && $map{$curSmpl}{inputFileSizeMB} >= $MFconfig{skipSmallSmplsMB}){
			system "rm -f $curOutDir/SMPL.empty"; $map{$curSmpl}{inputFilesEmpty} = 0;
	}
	#print "small skip: $MFconfig{skipSmallSmplsMB} $map{$curSmpl}{inputFileSizeMB} $map{$curSmpl}{inputXFileSizeMB} $AssemblyGo \n";
	if (-e "$curOutDir/SMPL.empty" || $jdep eq "EMPTY_DO_NEXT" || (($map{$curSmpl}{inputFileSizeMB}  + $map{$curSmpl}{inputXFileSizeMB}) < $MFconfig{skipSmallSmplsMB} ) && (!$AssemblyGo || $AsGrps{$cAssGrp}{CntAimAss} <= 1 ) ){
		
		if ( ($map{$curSmpl}{inputFileSizeMB} + $map{$curSmpl}{inputXFileSizeMB}) < $MFconfig{skipSmallSmplsMB} ){
			print "Skipping sample $curSmpl due to $map{$curSmpl}{inputFileSizeMB} < $MFconfig{skipSmallSmplsMB} MB\n";
			#still create essentials
			#	my $coveragePerCtg = "$ContigStatsDir/Coverage.percontig.gz";
			#!$eFinMapCovGZ && $eCovAsssembly) || ($eSuppCovAsssembly && !$eFinSupMapCovGZ)
			my $geneCovTmpFile = $coveragePerCtg; $geneCovTmpFile =~ s/percontig.gz$/count_pergene/;
			my $geneCntTmpFile = $coveragePerCtg; $geneCntTmpFile =~ s/percontig.gz$/pergene/;
			my $geneMedTmpFile = $coveragePerCtg; $geneMedTmpFile =~ s/percontig.gz$/median.pergene/;
			my $tmpMapCovGZ = "$finalMapDir/$SmplName-smd.bam.coverage.gz";
			system "mkdir -p $ContigStatsDir" unless (-d $ContigStatsDir);
			system "mkdir -p $finalMapDir" unless (-d $finalMapDir);
			system "mkdir -p $curOutDir" unless (-d $curOutDir); $coveragePerCtg =~ s/\.gz$//;
			system "touch $tmpMapCovGZ $coveragePerCtg $geneCovTmpFile $geneCntTmpFile $geneMedTmpFile $curOutDir/SMPL.empty";
		} else {
			#print "Sample empty.. next\n";
		}
		$runReport{empty_samples}{$curSmpl} =
			$map{$curSmpl}{inputFileSizeMB} + $map{$curSmpl}{inputXFileSizeMB};
		reduceProgStats(); #reduce counters for riboFind etc.. no find in this sample!
		MFnext($smplLockF,\@sampleDeps,$JNUM ,$QSBoptHR); 
		loop2C_check($cAssGrp,\@sampleDeps);next;
	}
	append_job_dependencies(\$AsGrps{$cAssGrp}{SeqClnDeps}, $jdep) if ($assemblyFlag);
	if (deferLoopProducerWave(
			'input staging', $jdep, $smplLockF, $cAssGrp, \@sampleDeps,
	)) { next; }
	
	
	#$mmpuOutTab .= $dir2rd."\t".$seqSet{"mmpu"}."\n";
	append_job_dependencies(\$AsGrps{$cMapGrp}{SeqUnZDeps}, $jdep);
	append_job_dependencies(\$AsGrps{$cAssGrp}{UnzpDeps}, $jdep);
	$AsGrps{$cAssGrp}{readDeps} = $AsGrps{$cAssGrp}{UnzpDeps};
	my $UZdep = $jdep;
	
	
	#empty links for assembler and nonpareil
	#my($arp1,$arp2,$singAr,$matAr,$sdmjN) = ([],[],[],[],"");
	my $sdmjN = ""; #job on main  reads
	#upload2EBI?? -> this has to happen before sdm etc, ensuring raw reads are being processed, in the same files as initially used...
	my $uplJob = uploadRawFilePrep($smplTmpDir."uploadPrep/",$curSmpl,$jdep,0);
	my $uplJobX = uploadRawFilePrep($smplTmpDir."uploadPrep/",$curSmpl,$jdep,1);
	
	push (@EBIjobs, $uplJob,$uplJobX);


	#empty links and objects for merging of reads
	my($mergJbN) = ("");
	# punsh the whole thing through sdm.. 
	if ( (!$boolAssemblyOK||$calcContamination) && $MFopt{useSDM}!=0 ){#&& !$is3rdGen) {
		$sdmjN  = sdmClean($curOutDir, $smplTmpDir."seqClean/",$jdep,$dowstreamAnalysisFlag,0) ;
		# check for support reads as well..
		#job on support  reads
		my $sdmjN2 = sdmClean($curOutDir, $smplTmpDir."seqClean/",$jdep,$dowstreamAnalysisFlag,1) ;
		$sdmjN .= ";$sdmjN2" if ($sdmjN2 ne "");
	}  
	append_job_dependencies(\$AsGrps{$cAssGrp}{SeqClnDeps}, $sdmjN) if ($assemblyFlag);
	if (deferLoopProducerWave(
			'quality filtering', $sdmjN, $smplLockF, $cAssGrp, \@sampleDeps,
	)) { next; }
	#adds raw and cleaned read file location to the whole assembly group
	
	my $cleanSeqSetHR = sampleReadSet($curSmpl, "clean");
	my $assGrpHR = addFileLocs2AssmGrp(
		\%AsGrps, $cAssGrp, $SmplName, $cleanSeqSetHR,
		sampleReadSet($curSmpl, "raw"),
	);
	%AsGrps = %{$assGrpHR};
	
	
	
	#make sure that no mapping is started when postAssmbl not yet done..
	#needs to happen after reads registered (for assembly of multiple samples..
	#$MappingGo=0 
	if ( !$efinAssLoc ){ #there should be no support mapping, unless the hybrid assembly has finished!
		$mapSuppAssFlag=0;$eSuppCovAsssembly=0;
	}
	if ($ePreAssmblPck && !$efinAssLoc && !$postPreAssmblGo){$mapAssFlag =0; $mapSuppAssFlag=0;}
	#print "mapCHK $mapSuppAssFlag $mapAssFlag $ePreAssmblPck && ! $efinAssLoc && ! $postPreAssmblGo\n";
	if ( !$mapAssFlag && $boolGenePredOK && ( ( ($eCovAsssembly || $ePreAssmblPck) && !$postPreAssmblGo && !$efinAssLoc ) #nothing to do until $doPreAssmFlag releases
			|| (!$ePreAssmblPck && $eCovAsssembly && $ePreAssmbly && !$postPreAssmblGo && !$efinAssLoc )  )
	){ #last sample (assembly) should not map while other maps are still running..
		print "next due to waiting for preassemblies..";
		print " flags: !$ePreAssmblPck && $ePreAssmbly && !$postPreAssmblGo && !$efinAssLoc \n";
		MFnext($smplLockF,\@sampleDeps,$JNUM ,$QSBoptHR); loop2C_check($cAssGrp,\@sampleDeps); next;
	}


	#filter human or other hosts..
	$sdmjN = removeHostSeqs($nodeSpTmpD,$sdmjN,1) if ($dowstreamAnalysisFlag && (!$boolAssemblyOK || $calcContamination));
	append_job_dependencies(\$AsGrps{$cAssGrp}{SeqClnDeps}, $sdmjN) if ($assemblyFlag);
	if (deferLoopProducerWave(
			'host filtering', $sdmjN, $smplLockF, $cAssGrp, \@sampleDeps,
	)) { next; }
	#merge reads?
	($mergJbN) = mergeReads($sdmjN,$smplTmpDir."merge_clean/",$calcReadMerge,$dowstreamAnalysisFlag);
	if (deferLoopProducerWave(
			'read merging', $mergJbN, $smplLockF, $cAssGrp, \@sampleDeps,
	)) { next; }
	
	#raw files only required for mapping reads to assemblies, so delete o/w
	#$cfp1ar,$cfp2ar,
	if (!$MFopt{DoAssembly} && $MFconfig{importMocat}==0 && $MFconfig{removeInputAgain} && !$requireRawReadsFlag){ $sdmjN = cleanInput($sdmjN,$smplTmpDir);}
	append_job_dependencies(\$AsGrps{$cAssGrp}{readDeps}, $mergJbN);


	#keeps track of all sdm jobs
	$d2Inputs{dependencies} .= ";".$sdmjN;
	
	#TODO
	my $primaryCleanLibraries = readLibrariesByScope($cleanSeqSetHR, 'primary', 1, $curSmpl);
	push(@{$d2Inputs{filtered_read1}}, @{libraryFiles($primaryCleanLibraries, 'r1')});
	push(@{$d2Inputs{filtered_read2}}, @{libraryFiles($primaryCleanLibraries, 'r2')});
	#die;
	if ($MFconfig{unpackZip}){
		print "next due to onlyFilterZip 1\n";MFnext($smplLockF,\@sampleDeps,$JNUM ,$QSBoptHR); 
		loop2C_check($cAssGrp,\@sampleDeps);next;
	}
	
	
	
	
	
	#-----------------------------------------------------------------
	#---------functions only dependent on UN/FILTERED reads-----------
	#-----------------------------------------------------------------
	#take long reads and filter for very long reads (454 etc might be long enough)
	#then do gene predictions *instead* of gene predictions on assemblies
	#die $pseudAssFlag."\n";
	
	#create new dependency flag for both unzip and sdm dep:
	my $primaryDep = "$UZdep;$sdmjN";
	
	if ($pseudAssFlag ){ #don't do assembly directly, use pre existing files instead..
		my ($psAssDep, $psFile, $metagDir) = createPsAssLongReads($mergJbN.";".$primaryDep, $pseudoAssFile, $finalCommAssDir, $SmplName);#pseudoAssFileFinal
		if ($psAssDep ne ""){
			$AsGrps{$cAssGrp}{pseudoAssmblDep} = $psAssDep;
		}
		#predict genes on assembly
		my $prodRun = genePredictions($psFile,$metagDir."/genePred/",$psAssDep,$finalCommAssDir,"","$nodeSpTmpD/genePred/",1);
		if ($prodRun ne ""){
			$AsGrps{$cAssGrp}{pseudoAssmblDep}  = $prodRun;
		}
		append_job_dependencies(\$AsGrps{$cAssGrp}{readDeps}, $prodRun);
	}
	
	if ($MFopt{calcOrthoPlacement}){
		runOrthoPlacement($curOutDir."orthos/",$nodeSpTmpD."/orthos/",
					$mergJbN.";".$primaryDep);
	}
	if ($calcDiamond || $calcDiaParse){
		my ($djname,$djCln) = runDiamond($curOutDir."diamond/",$baseOut."DB/DiamDB/",$nodeSpTmpD."/diaRefDB/",
					$mergJbN.";".$primaryDep,$MFopt{reqDiaDB}); #GlbTmpPath
		$AsGrps{$cAssGrp}{DiamDeps} = normalise_job_dependencies($djname);
		append_job_dependencies(\$AsGrps{global}{DiamDeps}, $djname);
		$AsGrps{global}{DiamCln} = $djCln unless($djCln eq "");
		append_job_dependencies(\$AsGrps{$cAssGrp}{readDeps}, $djname);
	}
	
	#DEBUG
	#my ($arp1,$arp2,$singAr,$mergRdsHr) = ("","","","");
	#non pareil (estimate community size etc)
	if ($nonPareilFlag){
		my $globalNPD = $baseOut."NonPareil/";
		nopareil(libraryFiles($primaryCleanLibraries, 'r1'),$nonParDir, $globalNPD, $SmplName,$primaryDep);
		MFnext($smplLockF,\@sampleDeps,$JNUM ,$QSBoptHR); 
		loop2C_check($cAssGrp,\@sampleDeps);next;
	}
	if ($calcGenoSize){#use microbeCensus to get avg genome size
		my $gsJdep = genoSize($curOutDir."MicroCens/",$mergJbN.";".$primaryDep);
		append_job_dependencies(\$AsGrps{$cAssGrp}{readDeps}, $gsJdep);
	}
	
	#kraken (estimate taxa abundance
	if ($calcKraken){
		my $krJdep = krakenTaxEst($KrakenOD, $nodeSpTmpD."krak/",$SmplName,$primaryDep);
		append_job_dependencies(\$AsGrps{$cAssGrp}{readDeps}, $krJdep);
		
	}
	if ($calcRibofind || $calcRiboAssign){
#		die "STOP ribo\n";
		my $ITSrun = detectRibo($nodeSpTmpD."ITS/",$curOutDir."ribos/",$primaryDep,$SmplName,$baseOut."DB/");
		#$AsGrps{$cAssGrp}{ITSDeps} .= $ITSrun.";";
		append_job_dependencies(\$AsGrps{$cAssGrp}{readDeps}, $ITSrun);
	}
	
	#metaphlan2 - taxa abudnance estimates
	if ($calcMetaPhlan){
		my $dir_MP2 = $baseOut.$preDIRs{dir2MePhl};#"pseudoGC/Phylo/MP2/"; #metaphlan 2 dir

		my $MP2jname = metphlanMapping($nodeSpTmpD."MP2/",$dir_MP2,$SmplName,$MFopt{MapperCores},$primaryDep); #\@cfp1,\@cfp2
		append_job_dependencies(\$AsGrps{$cAssGrp}{readDeps}, $MP2jname);
	}
	
	if ($calcTaxaTar){
		my $dir_TaxTar = $baseOut."pseudoGC/Phylo/TaxaTarget/"; #taxaTar dir
		my $taxTarjname = TaxaTarget($nodeSpTmpD."TaxTar/",$dir_TaxTar,$SmplName,$MFopt{MapperCores},$primaryDep); 
		append_job_dependencies(\$AsGrps{$cAssGrp}{readDeps}, $taxTarjname);
	}
	
	#mOTU2  - taxa abundance estimates
	if ($calcMOTU2){
		my $dir_mOTU2 = $baseOut."pseudoGC/Phylo/mOTU2/"; #mOUT 2 dir
		my $MP2jname = mOTU2Mapping($nodeSpTmpD."Motu2/",$dir_mOTU2,$SmplName,$MFopt{MapperCores},$primaryDep); #\@cfp1,\@cfp2
		append_job_dependencies(\$AsGrps{$cAssGrp}{readDeps}, $MP2jname);
	}
	
	my $SmplNameX = $SmplName;
	if ($AsGrps{$cAssGrp}{CntAss} > 1){$SmplNameX .= "M".$AsGrps{$cAssGrp}{CntAss};}
	#my @tmp = @{$arp1};die ("ASSflag: ".$assemblyFlag."\n@tmp\n");
	#die "$assemblyFlag\n $AssemblyGo\n";
	
	
	#add dependencies to wait on for loopUntil
	add2SampleDeps(\@sampleDeps, [$AsGrps{$cAssGrp}{readDeps},$AsGrps{$cAssGrp}{DiamDeps},$AsGrps{$cMapGrp}{SeqUnZDeps},$mergJbN,$sdmjN]);
	






	#-----------------------------------------------------------------
	#------------------------ ASSEMBLY -------------------------------
	#-----------------------------------------------------------------
	append_job_dependencies(\$AsGrps{$cAssGrp}{SeqClnDeps}, $sdmjN) if ($assemblyFlag);
	if ($assemblyFlag && deferLoopProducerWave(
			'input preparation', $sdmjN, $smplLockF, $cAssGrp, \@sampleDeps,
	)) { next; }
	if ($AssemblyGo && deferLoopProducerWave(
			'assembly-group input preparation', $AsGrps{$cAssGrp}{SeqClnDeps},
			$smplLockF, $cAssGrp, \@sampleDeps,
	)) { next; }

	if ( ($assemblyFlag || $scaffoldFlag || $scaffTarExternal ne "") && $AssemblyGo){ #assembly does not exist
		die "Can't do assembly and pseudoassembly on the same sample!\n" if ($pseudAssFlag || $MFopt{pseudoAssembly});
		#print "preAsmChk: $ePreAssmbly, $ePreAssmblPck, $doPreAssmFlag, $postPreAssmblGo\n";
		#die;
		metagAssemblyRun( $cAssGrp,"$nodeSpTmpD/ass",$metagAssDir ,$geneDir,  $SmplNameX,$scaffoldFlag,$metaGscaffDir,
					$assemblyFlag,$AssemblyGo,$ePreAssmbly, $doPreAssmFlag, $postPreAssmblGo,$finalCommAssDir);
		if (deferLoopProducerWave(
				'assembly', $AsGrps{$cAssGrp}{AssemblJobName},
				$smplLockF, $cAssGrp, \@sampleDeps,
		)) { next; }
		my $producedAssemblyDir = ($MFopt{DoAssembly} == 5 && $doPreAssmFlag)
			? $metagAssDir : $finalCommAssDir;
		$metaGassembly = "$producedAssemblyDir/scaffolds.fasta.filt";
		$geneDir = "$producedAssemblyDir/genePred/";
		#call genes, depends on assembly
		$AsGrps{$cAssGrp}{prodRun} = genePredictions($metaGassembly,$geneDir,$AsGrps{$cAssGrp}{AssemblJobName},$finalCommAssDir,"","$nodeSpTmpD/genePred/",1);
	} elsif ($boolAssemblyOK && !$locRedoAssMapping) {
		$runReport{present_assemblies}++;
	}
	
	#die "$assemblyFlag || ($ePreAssmbly && $doPreAssmFlag) \n";
	if (!$assemblyFlag || ($ePreAssmbly && $doPreAssmFlag) ){   # gene predictions on assembly, assemblies already do exist
		$metaGassembly = $finAssLoc; #print "No Assembly routines required\n" if ($MFopt{DoAssembly}==0);
		#print "GP:: !$boolGenePredOK && $AssemblyGo \n";
		if (!$boolGenePredOK && $AssemblyGo ){
			$geneDir = $finalCommAssDir."/genePred/";
			$AsGrps{$cAssGrp}{prodRun} = genePredictions($metaGassembly,$geneDir,$AsGrps{$cAssGrp}{AssemblJobName},$finalCommAssDir,"","$nodeSpTmpD/genePred/",1);
		}
	}
	if (deferLoopProducerWave(
			'gene prediction', $AsGrps{$cAssGrp}{prodRun},
			$smplLockF, $cAssGrp, \@sampleDeps,
	)) { next; }
	

	
	my $currentMappingDeps = '';
	my $finalAssemblyScheduled = $efinAssLoc || $MFopt{DoAssembly} != 5 || $postPreAssmblGo;
	if ($AssemblyGo && $finalAssemblyScheduled && $AsGrps{$cAssGrp}{PostAssemblCmd} ne "") {
		print "Submitting deferred assembly-group mapping jobs\n";
		my $deferredDeps = postSubmQsub(
			"$logDir/MultiMapper.sh", $AsGrps{$cAssGrp}{PostAssemblCmd},
			$AsGrps{$cAssGrp}{AssemblJobName},
		);
		append_job_dependencies(\$AsGrps{$cAssGrp}{MapDeps}, $deferredDeps);
		append_job_dependencies(\$AsGrps{$cAssGrp}{BinDeps}, $deferredDeps);
		append_job_dependencies(\$currentMappingDeps, $deferredDeps);
		add2SampleDeps(\@sampleDeps, [$deferredDeps]);
		$AsGrps{$cAssGrp}{PostAssemblCmd} = "";
	}
	
	#-----------------------------------------------------------------
	#------------------------  MAPPING -------------------------------
	#-----------------------------------------------------------------

	
	#-----------------   mapping to other tars -----------------
	#2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd 2nd
	# maps on given reference genome(s), not on the assembly itself!!
	if (@bwt2outD>0 && $MappingGo && (!$boolScndMappingOK || !$boolScndCoverageOK || $calc2ndMapSNP) ){#map reads to specific tar
		scndMap2Genos($SmplName,$cMapGrp,$cAssGrp,$curOutDir,$nodeSpTmpD,
			\@sampleDeps,$calc2ndMapSNP,$boolScndCoverageOK);
	} elsif ($MFopt{mapModeActive}) {
		MFnext($smplLockF,\@sampleDeps,$JNUM ,$QSBoptHR); 
		loop2C_check($cAssGrp,\@sampleDeps); next;
	}
	
	
	#%%%%%%%%%%%%%%%%   functions dependent on assembly -> submit post-assembly   #%%%%%%%%%%%%%%%%
	#  ---------  mapping 2  assembly  -----------------------   map reads to Assembly      ------------------------------

	#need index for mapper?
	my $doMapping = $MappingGo && !$eFinMapCovGZ && $MFopt{map2Assembly} && ($MFopt{DoAssembly} || $mapAssFlag);
	
	my $assemblyBuildIndexFlag=0; 	$assemblyBuildIndexFlag=1 if (-s $finAssLoc && ($doMapping || $mapSuppAssFlag));
	#die "$assemblyBuildIndexFlag  $MFopt{DoAssembly}  && -s $finAssLoc  && $MFopt{map2Assembly} && ($mapAssFlag || $mapSuppAssFlag ) \n". mapperDBbuilt($finAssLoc,$MFopt{MapperProg})  ."\n";
	#print "build $assemblyBuildIndexFlag   $MFopt{DoAssembly} && !$assemblyFlag && $MFopt{map2Assembly} && $mapAssFlag && $MFopt{MapperProg}\n";
	if ($assemblyBuildIndexFlag && $AsGrps{$cAssGrp}{AssemblJobName} eq ""){ #in this case asembly was done, but index was never built
		buildAssemblyMapIdx($finAssLoc, $cAssGrp, $mapAssFlag,$mapSuppAssFlag,$SmplName);
		if (deferLoopProducerWave(
				'mapping index', $AsGrps{$cAssGrp}{AssemblJobName},
				$smplLockF, $cAssGrp, \@sampleDeps,
		)) { next; }
	}

	
	my $mappingDeferred = 0;
	if ($doMapping){ #mapping to the assembly (can be multi-sample assembly as well)
		my $mapNow =0;
		$mapNow = 1 if ($AssemblyGo ||  $efinAssLoc || ($ePreAssmbly && $doPreAssmFlag) );#controls if several samples are assembled together and this needs to be waited for
		#print "main map\n";
		my $unAlDir = "$finalMapDir/unaligned/";$unAlDir = "" if (!$MFopt{SaveUnalignedReads}); #in most cases we don't need unaligned reads..
		my %dirset = 	(nodeTmp=>$nodeSpTmpD,outDir => "$finalMapDir/", unalDir => $unAlDir,
						sbj => $metaGassembly, assGrp => $cAssGrp,  smplName => $SmplName,mappingStarted =>1,
						glbTmp => $nodeSpTmpD."_mapWork/",glbMapDir => $finalMapDir,mapSupport => 0,
						libraries => getRawLibrariesAssmGrp(\%AsGrps,$cAssGrp,0,$SmplName),
						readTec => '', submit => 1,submNow => $mapNow,
						sortCores => $MFopt{bamSortCores}, mapCores => $MFopt{MapperCores}, cramAlig => $cramthebam,
						# Primary short- and long-read assemblies always receive a breakpoint report.
						breakpointOutput => "$finalMapDir/$SmplName-smd.bam.breakpoints.tsv.gz",
						strictHybridCoverage => ($MFopt{DoAssembly} == 5 && $ePreAssmbly && $doPreAssmFlag ? 1 : 0));
		# primary mapping (onto de novo assembly)
		my ($map2Ctgs,$delaySubmCmd,$mapOptHr) = mapReadsToRef(\%dirset, normalise_job_dependencies($AsGrps{$cAssGrp}{AssemblJobName}, $jdep));#\@libsCFP);
		my ($map2Ctgs_2,$delaySubmCmd_2,$mapStat)  = bamDepth(\%dirset,$map2Ctgs,$mapOptHr);
		$delaySubmCmd .= "\n".$delaySubmCmd_2;
		append_job_dependencies(\$AsGrps{$cAssGrp}{MapDeps}, $map2Ctgs_2);
		append_job_dependencies(\$currentMappingDeps, $map2Ctgs_2);
		append_job_dependencies(\$AsGrps{$cAssGrp}{BinDeps}, $map2Ctgs_2);
		if (!${$mapOptHr}{immediateSubm} && $delaySubmCmd =~ /\S/ ){
			$mappingDeferred = 1;
			#store command for later..
			$AsGrps{$cAssGrp}{PostAssemblCmd} .= $delaySubmCmd;
		}
	}
	
	
	if ($mapSuppAssFlag){ #supplementary mappings (eg PacBio in hybrid assemblies)
		print "mapping support reads\n";
		#print "DEBUG:: !$doPreAssmFlag && !$ePreAssmblPck\n";
		my $mapNow = 1;
		my $unAlDir = "$finalMapDir/unaligned_supp/";$unAlDir = "" if (!$MFopt{SaveUnalignedReads});
		my %dirset = 	(nodeTmp=>$nodeSpTmpD,outDir => "$finalMapDir/", unalDir => $unAlDir,
						sbj => $metaGassembly, assGrp => $cAssGrp,  smplName => $SmplName,
						glbTmp => $nodeSpTmpD."_mapWorkSupp/",glbMapDir => $finalMapDir, mapSupport => 1,
						libraries => getRawLibrariesAssmGrp(\%AsGrps,$cAssGrp,1,$SmplName),
						readTec => "", submit => 1,submNow => $mapNow,mappingStarted=>1,
						sortCores => $MFopt{bamSortCores}, mapCores => $MFopt{MapperCores}, cramAlig => $cramthebam);
		# primary mapping of support reads (onto de novo assembly)
		my ($mapSup2Ctgs,$delaySubmCmd,$mapOptHr) = mapReadsToRef(\%dirset, normalise_job_dependencies($AsGrps{$cAssGrp}{AssemblJobName}, $jdep));#\@libsCFP);
		my ($mapSup2Ctgs_2,$delaySubmCmd_2,$mapStat)  = bamDepth(\%dirset,$mapSup2Ctgs,$mapOptHr);
			$delaySubmCmd .= "\n".$delaySubmCmd_2;
		append_job_dependencies(\$AsGrps{$cAssGrp}{MapDeps}, $mapSup2Ctgs_2);
		append_job_dependencies(\$currentMappingDeps, $mapSup2Ctgs_2);
		append_job_dependencies(\$AsGrps{$cAssGrp}{BinDeps}, $mapSup2Ctgs_2);
		if ($delaySubmCmd =~ /\S/ && !${$mapOptHr}{immediateSubm}) {
			$mappingDeferred = 1;
			$AsGrps{$cAssGrp}{PostAssemblCmd} .= $delaySubmCmd;
		}
#		die;
	}
	
	my $mappingWaveDeps = normalise_job_dependencies(
		$currentMappingDeps, $AsGrps{$cAssGrp}{MapDeps},
	);
	if (deferLoopProducerWave(
			'assembly mapping', $mappingWaveDeps,
			$smplLockF, $cAssGrp, \@sampleDeps,
	)) { next; }
#		die "$MappingGo && !$eFinMapCovGZ && $MFopt{map2Assembly} && ($MFopt{DoAssembly} || $mapAssFlag)";

	
	#-----------------------------------------------------------------
	#---------------- producer barriers, downstream analysis ---------------
	#-----------------------------------------------------------------
	my ($fullContigStatsDep, $binningJobDep, $variantJobDep) = ("", "", "");
	my $currentContigStatsDeps = '';
	
	
	# Completed assemblies and mappings are published by their producer jobs.
	# This is the explicit publication barrier for downstream work; there is no
	# separate result-moving cleanup job.
	my $publicationDeps = normalise_job_dependencies(
		$jdep, $AsGrps{$cAssGrp}{MapDeps}, $AsGrps{$cAssGrp}{prodRun},
		$AsGrps{$cAssGrp}{AssemblJobName}, $uplJob,
	);
	append_job_dependencies(\$AsGrps{$cAssGrp}{BinDeps}, $publicationDeps);
	add2SampleDeps(\@sampleDeps, [$jdep , $AsGrps{$cAssGrp}{MapDeps} , $AsGrps{$cAssGrp}{scndMapping},$AsGrps{$cAssGrp}{prodRun} ]);

	#die;
	
	#this flag is essentially $allCovDone 
	#my $needsContigStats =0; #flag to activate contigStats later..
	#$needsContigStats = 1 if (($MappingGo && !$eCovAsssembly) || ($mapSuppAssFlag && !$eSuppCovAsssembly));	
	# calc statsitics concercing readqual, mappings, genes & contigs
	#print "COV: $calcCoverage \n";
	if ($pseudAssFlag || ($AssemblyGo && $MFopt{DoAssembly}) || 
				isLastSampleInAssembly($finalCommAssDir,$curOutDir) ) {
		my $subprts = $MFconfig{defaultContigSubs}."gFG"; $subprts .= "m" if ($MFopt{DoBinning});
		$subprts .= "k" if ($MFopt{kmerAssembly} );$subprts .= "4" if ($MFopt{kmerPerGene});

		if ($MFopt{DoMetaBat2} == 4 && $subprts !~ m/F/) {
			$subprts .= "F";
			print "GenomeFace requires FetchMG - adding it as a step in contigStats\n";
		}

		my ($contRun, $tmpCDd) = (runContigStats(
			$curOutDir, normalise_job_dependencies($publicationDeps,$AsGrps{$cAssGrp}{prodRun}),
			$finalCommAssDir, $subprts, 1, $nodeSpTmpD, 1, 6, $curSmpl,
			$supportCoverageRequired,
		))[0, 2];

		#run contig stats
		my $deferredContigDeps = postSubmQsub(
			"$logDir/MultiContigStats.sh", $AsGrps{$cAssGrp}{PostClnCmd}, $contRun,
		);
		$AsGrps{$cAssGrp}{PostClnCmd} = "";
		$jdep = normalise_job_dependencies($contRun, $deferredContigDeps);
		append_job_dependencies(\$currentContigStatsDeps, $jdep);
		$fullContigStatsDep = $jdep if ($tmpCDd && $jdep ne "");
		append_job_dependencies(\$AsGrps{$cAssGrp}{BinDeps}, $deferredContigDeps);

		if ($contRun ne "") {
			append_job_dependencies(\$AsGrps{$cAssGrp}{BinDeps}, $contRun);
			print "Added main contig stats as a downstream dependency\n"
				if ($MFopt{DoMetaBat2} == 4);
		}

	} elsif (((exists($AsGrps{$cAssGrp}{MapDeps}) && $AsGrps{$cAssGrp}{MapDeps} =~ m/[^;\s]/ ) || $calcCoverage) ) {
		#die "test23  $AsGrps{$cAssGrp}{MapDeps}\n";
		#calculate solely abundance / gene after producer publication and assembly contig stats
		my $submitContigNow = $mappingDeferred ? 0 : 1;
		my ($jn,$delaySubmCmd2,$tmpCDd) = runContigStats($curOutDir,$publicationDeps,$finalCommAssDir,$MFconfig{defaultContigSubs},$submitContigNow,$nodeSpTmpD,$AssemblyGo,1, $curSmpl,$supportCoverageRequired);
		$AsGrps{$cAssGrp}{PostClnCmd} .= $delaySubmCmd2;
		$jdep = $jn;
		append_job_dependencies(\$currentContigStatsDeps, $jdep);
		append_job_dependencies(\$AsGrps{$cAssGrp}{BinDeps}, $jdep) if ($jdep ne "");
	}
#	die;
	add2SampleDeps(\@sampleDeps, [$publicationDeps,$jdep]);
	my $contigStatsWaveDeps = normalise_job_dependencies(
		$currentContigStatsDeps, $AsGrps{$cAssGrp}{BinDeps},
	);
	if (deferLoopProducerWave(
			'contig statistics', $contigStatsWaveDeps,
			$smplLockF, $cAssGrp, \@sampleDeps,
	)) { next; }
	#Binning, SNP calling: only after copying files from tmp and running contig stats
	if ( $calcBinning && $AssemblyGo ){  #$allMapDone rm: this is checked now via $AsGrps{$cAssGrp}{MapDeps}
		my $binnerTmp = $nodeSpTmpD;
		$binnerTmp = $smplTmpDir if ($MFopt{useBinnerScratch});
		$binningJobDep = submitGenomeBinner($binnerTmp,$finAssLoc,$BinningOut, $cAssGrp,$smplIDs[-1]);
		add2SampleDeps(\@sampleDeps, [$binningJobDep]);
	}
	
	#die "AT CONS SNP\n$allMapDone\n";
	my $primaryVariantRequested = $calcConsSNP || $calcSVs;
	my $supportVariantRequested = $calcSuppConsSNP || $calcSVsSupp;
	my $variantWorkRequested = $primaryVariantRequested || $supportVariantRequested;
	my $geneConsensusRequested = $MFopt{saveConsFastas}
		&& ($calcConsSNP || $calcSuppConsSNP);
	# Assembly, annotation, mapping, and ContigStats outputs can all be valid
	# future outputs of jobs submitted in this pass. SNP region planning now runs
	# inside the Cons allocation, so scheduler dependencies can safely replace
	# launch-time file-existence checks for those pending products.
	my $assemblyDownstreamScheduled = $AssemblyGo && !$doPreAssmFlag && !$ePreAssmblPck
		&& ($efinAssLoc || $MFopt{DoAssembly} != 5 || $postPreAssmblGo)
		&& (!$map{$curSmpl}{hasPrimaryRds} || ($MFopt{map2Assembly} && $MappingGo));
	my $assemblyDownstreamDeferred = $mappingDeferred
		&& !$doPreAssmFlag && !$ePreAssmblPck;
	my $variantCommonInputsPublished = $efinAssLoc
		&& (!$geneConsensusRequested || ($boolGenePredOK
			&& fileGZe("$finalCommAssDir/genePred/genes.gff")));
	my $primaryVariantInputsPublished = !$primaryVariantRequested || (
		-s "$finalMapDir/$SmplName-smd.$bamcramMap"
		&& (!$calcConsSNP || ($eFinMapCovGZ
			&& fileGZe("$finalMapDir/$SmplName-smd.bam.coverage") && $eCovAsssembly))
	);
	my $supportVariantInputsPublished = !$supportVariantRequested || (
		-s "$finalMapDir/$SmplName.sup-smd.$bamcramMap"
		&& (!$calcSuppConsSNP || ($eFinSupMapCovGZ
			&& fileGZe("$finalMapDir/$SmplName.sup-smd.bam.coverage") && $eSuppCovAsssembly))
	);
	my $variantInputsMayBePending = (!$variantCommonInputsPublished
			|| !$primaryVariantInputsPublished || !$supportVariantInputsPublished)
		&& ($assemblyDownstreamScheduled || $assemblyDownstreamDeferred);
	my $variantCommonInputsReady = $variantCommonInputsPublished
		|| $variantInputsMayBePending;
	my $primaryVariantInputsReady = $primaryVariantInputsPublished
		|| ($variantInputsMayBePending && $doMapping);
	my $supportVariantInputsReady = $supportVariantInputsPublished
		|| ($variantInputsMayBePending && $mapSuppAssFlag);
	my $variantSubmissionDeferred = $variantInputsMayBePending
		&& $assemblyDownstreamDeferred;
	if ($variantWorkRequested && $variantCommonInputsReady
			&& $primaryVariantInputsReady && $supportVariantInputsReady){
		my $rawReadSet = sampleReadSet($curSmpl, "raw");
		my $variantPrimaryTechnology = libraryTechnology(
			readLibrariesByScope($rawReadSet, "primary", 0, $curSmpl),
			"primary variant reads for $curSmpl", 0,
		);
		my $variantSupportTechnology = libraryTechnology(
			readLibrariesByScope($rawReadSet, "support", 0, $curSmpl),
			"support variant reads for $curSmpl", 0,
		);
		#die "conssnp:: $calcConsSNP $allMapDone $finalMapDir\n";
		#my $ofas = "$curOutDir/SNP/genePred/genes.shrtHD.SNPc.fna";
	
		my $consensusGff = "$finalCommAssDir/genePred/genes.gff";
		$consensusGff .= ".gz" if ($MFopt{GenePredGZ});
		my %SNPinfo = (gff => "$finalCommAssDir/genePred/genes.gff",
						assembly => $finAssLoc,
						mapD => "$finalMapDir",normIndels => $MFopt{normSNPindels},
						SNPcaller => $MFopt{SNPcallerFlag},hasPrimaryRds => $map{$curSmpl}{hasPrimaryRds},
						createFastas => $MFopt{saveConsFastas},
						saveVCF => $MFopt{saveVCF},
						ofas => $contigsSNP, #primary file of contigs
						genefna => $genePredSNP,genefaa => $genePredAASNP,
						vcfFile => $vcfSNP,vcfFileSupp => $vcfSNPsupp, gffFile => $consensusGff,
						nodeTmpD => $nodeSpTmpD,scratch => "$nodeSpTmpD/SNP/",
						smpl => $SmplName,bamcram => $bamcramMap,minDepth => $MFopt{consSNPminDepth},
						depthF => $coveragePerCtg,firstInSample => 1, #($i == 0 ? 1 : 0)
						bpSplit => 1e6,runLocal => 1,
						SeqTech => $variantPrimaryTechnology,
						SeqTechSuppl => $variantSupportTechnology,
						cmdFileTag => "ConsAssem",maxCores => $MFopt{maxSNPcores},#memReq => $MFopt{memSNPcall},
						jdeps => $AsGrps{$cAssGrp}{BinDeps},split_jobs => $MFopt{SNPconsJobsPsmpl},
						inputSizeMB => ($map{$curSmpl}{inputFileSizeMB} || 0)
							+ ($map{$curSmpl}{inputXFileSizeMB} || 0),
						allowPendingInputs => ($variantInputsMayBePending ? 1 : 0),
						immediateSubm => ($variantSubmissionDeferred ? 0 : 1),
						overwrite => $MFopt{redoSNPcons}, memPJob => $MFopt{memPJob},
						STOconSNP => $sampleCheckpoints{primaryConsensus}, STOconSNPsupp => "",
						minCallQual => $MFopt{SNPminCallQual},
						callConsSNP => $calcConsSNP, callConsSNPSupp => $calcSuppConsSNP,
						#struct vars
						callSVs => $MFopt{callSVs}, vcfSVfile => $vcfSV, vcfSVfileS => $vscSVsupp, callSVsSupp => $MFopt{callSVsSupp},
					);
		if ($calcSuppConsSNP){
			$SNPinfo{STOconSNPsupp} = $sampleCheckpoints{supportConsensus}; #trigger for also looking at cons SNP for support reads
		}
		
		my $variantSubmissionCommands;
		($variantJobDep,$variantSubmissionCommands) = createConsSNPandSVs(\%SNPinfo); #SNP calls on assembly
		$AsGrps{$cAssGrp}{PostConsCmd} .= $variantSubmissionCommands
			if ($variantSubmissionDeferred && defined($variantSubmissionCommands));
		add2SampleDeps(\@sampleDeps, [$variantJobDep]);
		#push(@sampleDeps, $consSNPdep) if (defined $consSNPdep && $consSNPdep ne "");
	}
	my $deferredVariantDeps = "";
	if ($AssemblyGo && ($AsGrps{$cAssGrp}{PostConsCmd} || "") =~ /\S/) {
		print "Submitting deferred assembly-group Cons jobs\n";
		$deferredVariantDeps = postSubmQsub(
			"$logDir/MultiConsensus.sh", $AsGrps{$cAssGrp}{PostConsCmd},
			normalise_job_dependencies($AsGrps{$cAssGrp}{BinDeps}, $publicationDeps, $jdep),
		);
		$AsGrps{$cAssGrp}{PostConsCmd} = "" if ($runOptions{submit});
		$variantJobDep = normalise_job_dependencies($variantJobDep, $deferredVariantDeps);
		add2SampleDeps(\@sampleDeps, [$deferredVariantDeps]);
	}
	my $cleanupBarrier = cleanup_stage_barrier(
		{
			name => 'final assembly publication', required => $assemblyOutputsRequired,
			complete => (!$doPreAssmFlag && !$ePreAssmblPck && $efinAssLoc),
			dependencies => (!$doPreAssmFlag && !$ePreAssmblPck && $AssemblyGo)
				? $publicationDeps : '',
		},
		{
			name => 'contig stats', required => $assemblyOutputsRequired,
			complete => $cleanupContigStatsComplete,
			dependencies => $fullContigStatsDep,
		},
		{
			name => 'binning', required => ($assemblyOutputsRequired && $MFopt{DoMetaBat2}) ? 1 : 0,
			complete => $binningComplete,
			dependencies => $binningJobDep,
		},
		{
			name => 'consSNP/variant analysis', required => ($assemblyOutputsRequired && $cleanupVariantRequired) ? 1 : 0,
			complete => $cleanupVariantsComplete,
			dependencies => $variantJobDep,
		},
	);
	if ($MFconfig{rmScratchTmp} && !$MFopt{DoCalcD2s} && @sampleDeps) {
		if ($cleanupBarrier->{ready}) {
			my @cleanupDependencies = split /;/, normalise_job_dependencies(
				\@sampleDeps, $cleanupBarrier->{dependencies},
			);
			my $cleanupJob = submitFinishedCleanup(
				"$logDir/FinishedCleanup.sh", "CLN$JNUM", \@cleanupDependencies,
				finishedCleanupArguments(
					$curSmpl, $SmplName, $finalCommAssDir, $finalMapDir,
					$smplTmpDir, $finAssLoc, $logDir, $cleanupRequirements,
					$assemblyOutputsRequired,
				),
			);
			add2SampleDeps(\@sampleDeps, [$cleanupJob]) if $cleanupJob ne '';
		} else {
			print "Deferring finished cleanup for $SmplName; waiting for "
				.join(', ', @{$cleanupBarrier->{blocked}})."\n";
		}
	}
	# Persist every sample owner, including cleanup. A later MATAF4 pass releases
	# the lock only after all recorded scheduler jobs have terminated.
	MFnext($smplLockF,\@sampleDeps,$JNUM ,$QSBoptHR);
	### loop2complete functionality
	loop2C_check($cAssGrp,\@sampleDeps);

	#print "END2\n@sampleDeps\n".@sampleDeps."\n";
}






print "\n\n###################################\n".$baseOut."\nFINISHED MATAFILER submission loop\n";

postprocess();


print "###################################\n\n";






close $QSBoptHR->{LOG};
exit(0);







































#--------------------------------------------------------------
#get ref genomes
#print("/g/bork5/hildebra/bin/cdbfasta/cdbfasta $PaulRefGenomes");
#if (!-f $thisRefSeq){system("/g/bork5/hildebra/bin/cdbfasta/cdbyank -a $GID $PaulRefGenomes.cidx > $thisRefSeq");}

#--------------------------------------------------------------
#and the assembly
#$cmd = "/g/bork5/hildebra/dev/Perl/assemblies/./multAss.pl $DBpath2 $assDir2 $assDir2/tmp/ $Ref 4 $CutSeq";

#####################################################
#
#####################################################

sub cleanupCompletionRequirements {
	my (%args) = @_;
	my (@exists, @nonempty, @files, @nonemptyFiles);
	my $contigDir = $args{contig_dir};
	$contigDir =~ s{/+$}{};
	if ($args{primary_coverage_required}) {
		push @exists, "$contigDir/Coverage.stone";
		push @files, map { "$contigDir/Coverage.$_" }
			qw(percontig median.percontig pergene count_pergene);
	}
	if ($args{support_coverage_required}) {
		push @exists, "$contigDir/Cov.sup.stone";
		push @files, map { "$contigDir/Cov.sup.$_" }
			qw(percontig median.percontig pergene count_pergene);
	}
	my $subparts = $args{contig_subparts} || '';
	push @nonemptyFiles, "$args{assembly_dir}/ContigStats/FMG/FMGids.txt"
		if $subparts =~ /F/;
	push @nonemptyFiles, "$args{assembly_dir}/ContigStats/GTDBmg/marker_genes_meta.tsv"
		if $subparts =~ /G/;
	push @nonemptyFiles, "$contigDir/scaff.pergene.4kmer.pm5"
		if $MFopt{kmerPerGene} && $subparts =~ /4/;

	if ($args{binning_base}) {
		push @exists, $args{binning_base};
		push @nonempty, "$args{binning_base}.assStat";
		push @exists, "$args{binning_base}.cm" if $MFopt{useCheckM1};
		push @nonempty, "$args{binning_base}.cm2" if $MFopt{useCheckM2};
	}
	push @exists, grep { defined($_) && $_ ne '' }
		@args{qw(primary_snp_stone support_snp_stone)};
	push @nonemptyFiles, grep { defined($_) && $_ ne '' }
		@args{qw(consensus_contigs primary_vcf support_vcf)};
	push @nonemptyFiles, @{$args{consensus_genes} || []};
	push @nonempty, grep { defined($_) && $_ ne '' }
		@args{qw(primary_sv support_sv)};
	return {
		exists => \@exists,
		nonempty => \@nonempty,
		files => \@files,
		nonempty_files => \@nonemptyFiles,
	};
}

sub finishedCleanupArguments {
	my ($sampleKey, $sampleName, $assemblyDir, $mappingDir,
		$sampleTemp, $assembly, $sampleLogDir, $requirements,
		$assemblyRequired) = @_;
	my @memberKeys = exists($map{$sampleKey}{AG_members})
		? @{$map{$sampleKey}{AG_members}} : ($sampleKey);
	my @memberArgs = map { ('--member', $map{$_}{SmplID}) } @memberKeys;
	my @memberLockArgs = map {
		('--member-lock', "$map{$_}{wrdir}/LOGandSUB/$MFcontstants{DefaultSampleLock}")
	} @memberKeys;
	my @arguments = (
		'--sample', $sampleName, @memberArgs, @memberLockArgs,
		'--state-dir', "$assemblyDir/.cleanup-indexes",
		'--allowed-root', $baseOut,
		'--mapping-dir', $mappingDir,
		'--sample-temp', $sampleTemp,
		'--scratch-root', $MFglobal{runTmpDirGlobal},
		($MFconfig{rmScratchTmp} ? '--remove-temporary' : '--no-remove-temporary'),
		'--snp-log-dir', "$sampleLogDir/SNP",
	);
	push @arguments,
		'--assembly', $assembly,
		'--assembly-path-file', "$map{$sampleKey}{wrdir}/assemblies/metag/assembly.txt",
		'--assembly-dir', $assemblyDir
		if $assemblyRequired;
	push @arguments, ('--remove-alignment', "$mappingDir/$sampleName-smd.cram")
		if $MFopt{map2Assembly} && !$MFopt{mapSaveCram} && $MFopt{DoMetaBat2};
	$requirements ||= {};
	my %requirementFlags = (
		exists => '--require-exists',
		nonempty => '--require-nonempty',
		files => '--require-file',
		nonempty_files => '--require-nonempty-file',
	);
	for my $kind (qw(exists nonempty files nonempty_files)) {
		push @arguments, map { ($requirementFlags{$kind}, $_) }
			@{$requirements->{$kind} || []};
	}
	return @arguments;
}

sub runFinishedCleanup {
	my @arguments = @_;
	my @cleaner = shellwords(getProgPaths("finished_sample_cleanup"));
	if (!@cleaner) {
		warn "Finished-sample cleanup command is not configured\n";
		return 0;
	}
	my $status = system(@cleaner, @arguments);
	if ($status != 0) {
		my $exitCode = $status == -1 ? -1 : $status >> 8;
		warn "Finished-sample cleanup failed with exit code $exitCode; retained temporary files\n";
		return 0;
	}
	return 1;
}

sub submitFinishedCleanup {
	my ($scriptFile, $jobName, $dependencies, @arguments) = @_;
	return '' unless ref($dependencies) eq 'ARRAY' && @{$dependencies};
	my @cleaner = shellwords(getProgPaths("finished_sample_cleanup"));
	if (!@cleaner) {
		warn "Finished-sample cleanup command is not configured\n";
		return '';
	}
	my $command = _shell_command(@cleaner, @arguments)."\n";
	my $dependencyString = normalise_job_dependencies($dependencies);
	return '' if $dependencyString eq '';

	my ($oldTmp, $oldShort, $oldAfterAny, $oldSubmissionConfig) =
		@{$QSBoptHR}{qw(tmpSpace useShortQueue afterAny submissionConfig)};
	$QSBoptHR->{tmpSpace} = 0;
	$QSBoptHR->{useShortQueue} = 1;
	$QSBoptHR->{afterAny} = 0;
	if (($QSBoptHR->{qmode} || '') eq 'slurm'
			&& ($QSBoptHR->{submissionConfig} || '') !~ /(?:^|;)--kill-on-invalid-dep=yes(?:;|$)/) {
		$QSBoptHR->{submissionConfig} = join(';',
			grep { defined($_) && $_ ne '' }
				($QSBoptHR->{submissionConfig}, '--kill-on-invalid-dep=yes'));
	}
	my ($cleanupJob, $submissionError);
	eval {
		($cleanupJob) = qsubSystem(
			$scriptFile, $command, 1, "1G", $jobName,
			$dependencyString, "", 1, [], $QSBoptHR,
		);
		1;
	} or $submissionError = $@ || 'unknown cleanup submission failure';
	@{$QSBoptHR}{qw(tmpSpace useShortQueue afterAny submissionConfig)} =
		($oldTmp, $oldShort, $oldAfterAny, $oldSubmissionConfig);
	die $submissionError if defined $submissionError;
	return $cleanupJob;
}

sub primeLoopSchedulerSnapshot {
	my ($start, $stop) = @_;
	return unless $runOptions{loopCount};
	$start = $selectedFrom if !defined($start) || $start < $selectedFrom;
	$stop = $selectedTo if !defined($stop) || $stop > $selectedTo;
	return if $start >= $stop;
	my @lockFiles = map {
		my $sampleKey = $samples[$_];
		"$map{$sampleKey}{wrdir}/LOGandSUB/$MFcontstants{DefaultSampleLock}";
	} $start .. $stop - 1;
	primeSampleLockJobSnapshot(\@lockFiles, $QSBoptHR);
}

sub deferLoopProducerWave {
	my ($reason, $dependencies, $lockFile, $assemblyGroup, $sampleDeps) = @_;
	return 0 unless $runOptions{loopCount};
	my $barrier = normalise_job_dependencies($dependencies);
	return 0 if $barrier eq '';
	add2SampleDeps($sampleDeps, [$barrier]);
	print "Producer wave '$reason' is pending for $curSmpl; "
		."deferring downstream jobs until a later loop pass.\n"
		unless $MFconfig{silent};
	MFnext($lockFile, $sampleDeps, $JNUM, $QSBoptHR);
	loop2C_check($assemblyGroup, $sampleDeps);
	return 1;
}


sub loop2C_check(){
	my ($cAssGrp,$sampleDeps_AR) = @_;
	if ($runOptions{loopCount} ){
		push (@grandDeps, @{$sampleDeps_AR});
		if ($JNUM == ($to-1)){
			my $submittedThisIteration =
				($QSBoptHR->{submittedJobs} || 0) - $loopIterationSubmissionStart;
			my $capacityDeferred = delete($QSBoptHR->{capacityDeferred}) ? 1 : 0;
			delete $QSBoptHR->{capacityDeferralAnnounced};

			my $activeJobs = 0;
			my $rerunLockedWindow = 0;
			# A pass that submitted work already reruns after its dependency wait.
			# Check the scheduler only for a no-op pass that may have skipped locks.
			if (!$loopFinalVerification && !$capacityDeferred && $runOptions{submit}
					&& !$MFconfig{rmSmplLocks} && $submittedThisIteration == 0) {
				$activeJobs = numActiveUserJobs(
					$QSBoptHR, 1, [keys %loopSubmittedJobIds],
				);
				$rerunLockedWindow = should_rerun_locked_window(
					active_jobs => $activeJobs,
					sample_count => $to - $from,
					active_job_threshold => $MFconfig{loopTillCompleteActiveJobs},
					remove_locks => $MFconfig{rmSmplLocks},
				);
			}
			if ($rerunLockedWindow
					&& ($runOptions{loopCount} > 1 || !$loopFinalLockRetryUsed)) {
				if ($runOptions{loopCount} > 1) {
					$runOptions{loopCount}--;
				} else {
					# Grant one final extra scan, but never spin indefinitely while a
					# small set of unrelated or long-running jobs remains active.
					$loopFinalLockRetryUsed = 1;
				}
				print "Retained sample locks with $activeJobs active job(s); " .
					"rerunning samples $from till $to before advancing the window.\n";
				@grandDeps = ();
				resetAsGrps(\%AsGrps);
				$loopIterationExtended = 0;
				$JNUM = $from - 1;
				primeLoopSchedulerSnapshot($from, $to);
				return;
			}
			my $lastWindowPass = $runOptions{loopCount} == 1 ? 1 : 0;
			my $overlapWindow = !$capacityDeferred && $runOptions{loopWindowSize} > 0
				? overlap_loop_window(
					to => $to, upper => $runOptions{to}, window_size => $runOptions{loopWindowSize},
					submitted_jobs => $submittedThisIteration,
					already_extended => $loopIterationExtended,
					last_pass => $lastWindowPass,
				)
				: { extended => 0 };
			if ($overlapWindow->{extended}) {
				my $previousTo = $to;
				$to = $overlapWindow->{to};
				$loopIterationExtended = 1;
				$loopFinalLockRetryUsed = 0;
				# The newly admitted block must receive a full retry allowance. The
				# current pass counts as its first iteration at the extended boundary.
				# Keep this pass's dependency list and submission snapshot open so the
				# eventual wait covers jobs from both blocks.
				$runOptions{loopCount} = $runOptions{loopInitialCount};
				my $overlapReason = $lastWindowPass
					? "Final loop pass"
					: "Light loop pass ($submittedThisIteration submitted job(s), threshold " .
						"$overlapWindow->{job_limit})";
				print "$overlapReason; extending sample window from " .
					"$from -> $previousTo to $from -> $to before waiting.\n";
				primeLoopSchedulerSnapshot($from, $to);
				return;
			}
			$loopIterationSubmissionStart = $QSBoptHR->{submittedJobs} || 0;
			$loopIterationExtended = 0;
			$loopFinalLockRetryUsed = 0;
			my $continueCurrentWindow = $submittedThisIteration > 0 || $capacityDeferred;
			if ($capacityDeferred) {
				# Capacity is temporary scheduler state, not a consumed workflow pass.
				$runOptions{loopCount} = $runOptions{loopInitialCount};
			} else {
				$runOptions{loopCount} = 0 unless $continueCurrentWindow;
				$runOptions{loopCount}-- if $continueCurrentWindow;
			}
			print "\n\n-------------------------------------------\n-------------------------------------------\n";
			if ($continueCurrentWindow) {
				if ($capacityDeferred) {
					print "Loop pass reached the Slurm job limit; the same range will be revisited.\n";
				} else {
					print "Completed loop iteration with $submittedThisIteration submitted job(s); " .
						($runOptions{loopCount} ? "repeating with $runOptions{loopCount} iteration(s) remaining.\n"
						                  : "iteration limit reached.\n");
				}
			} else {
				print "No jobs were submitted in the current iteration; ending this loop early.\n";
			}

			#print "L2C:: $runOptions{loopCount}  @{$sampleDeps_AR}\n";
			if ($continueCurrentWindow) {
				my $loopJobTag = $QSBoptHR->{rTag} || '';
				for my $loopJobId (
					split /;/, normalise_job_dependencies(\@grandDeps)
				) {
					$loopJobId =~ s/^\Q$loopJobTag\E//;
					$loopSubmittedJobIds{$loopJobId} = 1
						if ($loopJobId =~ /^\d+$/);
				}
				qsubSystemJobAlive(
					\@grandDeps, $QSBoptHR, 1,
					$MFconfig{loopTillCompleteActiveJobs},
				);
				if ($capacityDeferred) {
					sleep($MFconfig{schedulerPollSeconds});
				}
				# Reinspect after every completed submission pass. This lets hybrid
				# preassembly packages and assembly-group outputs become dependencies
				# for the next pass without requiring a separate user command.
				unless ($loopFinalVerification) {
					$workflowIteration++;
					runAutomaticWorkflowPreflight($workflowIteration) if ($MFconfig{autoStatePlan});
				}
			}
			# Reset pass-local state before either another iteration or a new window.
			@grandDeps = ();
			resetAsGrps(\%AsGrps);

			my $hadActiveLocks = $loopSawActiveLocks;
			$loopSawActiveLocks = 0;
			unless ($loopFinalVerification) {
				my $frontier = rolling_completed_frontier(
					from => $from, upper => $selectedTo,
					window_size => $runOptions{loopWindowSize},
					completed => \%loopSampleCompleted,
				);
				if ($frontier->{advanced}) {
					my $oldFrom = $from;
					$from = $frontier->{from};
					if ($runOptions{loopWindowSize} > 0) {
						$to += $frontier->{advanced};
						$to = $selectedTo if $to > $selectedTo;
					}
					print "Advanced completed-sample frontier from $oldFrom to $from; "
						."active scan is $from -> $to.\n";
				}
			}

			if ($loopFinalVerification) {
				my @invocationJobIds = sort { $a <=> $b }
					grep { /^\d+$/ }
					keys %{$QSBoptHR->{submittedJobRecords} || {}};
				my $liveInvocationJobs = @invocationJobIds
					? numLiveUserJobs($QSBoptHR, 1, \@invocationJobIds) : 0;
				if ($continueCurrentWindow || $hadActiveLocks || $liveInvocationJobs > 0) {
					my $reason = $continueCurrentWindow
						? 'newly submitted work'
						: $hadActiveLocks ? 'active sample locks'
						: "$liveInvocationJobs pending or running invocation job(s)";
					print "Final verification found $reason; waiting for "
						."$liveInvocationJobs active job(s) submitted by this invocation "
						."before another full-range pass.\n";
					if (@invocationJobIds) {
						# No active-job threshold here: pending and running jobs from this
						# invocation must all leave the scheduler before the next full scan.
						qsubSystemJobAlive(\@invocationJobIds, $QSBoptHR, 1);
					} elsif ($hadActiveLocks) {
						# Locks can predate this invocation. Avoid a tight rescan when no
						# locally submitted job ID is available to wait on.
						sleep($MFconfig{schedulerPollSeconds});
					}
					$workflowIteration++;
					runAutomaticWorkflowPreflight($workflowIteration)
						if ($MFconfig{autoStatePlan});
					$from = $selectedFrom;
					$to = $selectedTo;
					$runOptions{loopCount} = 1;
					$JNUM = $from - 1;
					print "Submitted jobs have finished; starting another full-range verification pass "
						."$from -> $to.\n";
				} else {
					$runOptions{loopCount} = 0;
					print "Final full-range verification completed; sample statistics may now be collected.\n";
				}
			} elsif ($runOptions{loopCount}) {
				$JNUM = $from - 1;
				print "Reanalyzing rolling sample range $from till $to\n";
			} elsif ($runOptions{loopWindowSize} > 0 && $to < $selectedTo) {
				my $previousTo = $to;
				$to += $runOptions{loopWindowSize};
				$to = $selectedTo if $to > $selectedTo;
				$runOptions{loopCount} = $runOptions{loopInitialCount};
				$JNUM = $from - 1;
				%loopSubmittedJobIds = ();
				print "Expanded rolling sample range $from -> $previousTo to $from -> $to.\n";
			} else {
				$loopFinalVerification = 1;
				$from = $selectedFrom;
				$to = $selectedTo;
				$runOptions{loopCount} = 1;
				$JNUM = $from - 1;
				%loopSubmittedJobIds = ();
				print "Starting final full-range verification pass $from -> $to before statistics.\n";
			}
			primeLoopSchedulerSnapshot($from, $to) if $runOptions{loopCount};
			print "-------------------------------------------\n-------------------------------------------\n";
			
		}
	}
}



sub postprocess{
	#print "\n\n###################################\nMain Loop done\n######################################\n";
	#global clean up cmds (like DB removals from scratch)
	print "Postprocessing:\n";
	if ($MFconfig{alwaysDoStats}) {
		for my $sampleName (@{$runReport{order}}) {
			my $context = $runReport{context}{$sampleName} || next;
			%locStats = ();
			$runReport{samples}{$sampleName} = {
				DIR => $context->{DIR},
				values => smplStats(
					$context->{input_dir}, $context->{assembly_dir}, $sampleName,
				),
			};
		}
	}
	if (@{$QSBoptHR->{submissionErrors} || []}) {
		print STDERR "MATAFILER continued after ".scalar(@{$QSBoptHR->{submissionErrors}})
			." scheduler submission problem(s). Existing sample locks were retained; failed work can be retried on a later run.\n";
	}

	#print input files, sorted by samples
	##transfer first..
	my %inputRawFQs;
	foreach my $cS (@samples){
		my $rawReadSet = sampleReadSet($cS, "raw");
		if (ref($rawReadSet) eq "HASH" && exists($rawReadSet->{rawReads})){
			$inputRawFQs{$cS} = $rawReadSet->{rawReads};
		} else {
			$inputRawFQs{$cS} = "";
		}
	}

	if (scalar(keys(%inputRawFQs)) == scalar(keys %{$d2Inputs{samples}})){
		open O,">$baseOut/Input_raw.txt";
		foreach my $sampleName (keys %{$d2Inputs{samples}}){
			print O "$sampleName\t$inputRawFQs{$sampleName}\n"
				if (exists($inputRawFQs{$sampleName}));
		}
		close O;
	}
	my $MGSfile = "$baseOut/metagStats.txt";
	my $MGShtml = "$baseOut/metagStatsReport.html";
	my $prevRep=0; #how many samples reported on?
	if (-e $MGSfile) {
		open my $previous, '<', $MGSfile;
		$prevRep++ while <$previous>;
		close $previous;
		$prevRep-- if $prevRep;
	}

	if (%{$runReport{samples}}) {
		my $statsText = _metag_stats_text($runReport{samples}, $runReport{order});
		my $temporary = "$MGSfile.tmp.$$";
		open my $statsFH, '>', $temporary or die "Cannot write temporary metagStats '$temporary': $!\n";
		print {$statsFH} $statsText or die "Cannot write metagStats data: $!\n";
		close $statsFH or die "Cannot close temporary metagStats '$temporary': $!\n";
		rename $temporary, $MGSfile or die "Cannot promote metagStats '$temporary' to '$MGSfile': $!\n";
	}
	print "Stats in $MGSfile \n";
	if ($MFopt{writeStats} && (@{$runReport{order}} > $prevRep || !-e $MGShtml ) && -s $MGSfile ){
		my $qcMakeHTMLReport = getProgPaths("qcMakeHTMLReport");
		my $Rpath = getProgPaths("Rpath");
		my $call = "$qcMakeHTMLReport $Rpath $MGSfile $MGShtml 1> /dev/null 2>&1;\n";
		#print $call."\n";
		#my $QCRes = `$call`; chomp $QCRes;
		system $call;
		#print $QCRes
		print "HTML Report in $MGShtml\n";
	}
	my $dir_MP2 = $baseOut.$preDIRs{dir2MePhl}; #metaphlan 2 dir
	#my $dir_RibFind = $baseOut.$preDIRs{dir2RiboF}; #metaphlan 2 dir
	#my $dir_KrakFind = $baseOut."pseudoGC/Phylo/KrakenTax/"; #metaphlan 2 dir
	my $dir_KrakFind = $baseOut."pseudoGC/Phylo/KrakenTax/$MFopt{globalKraTaxkDB}/"; #kraken dir
	my $dir_mOTU2 = $baseOut."pseudoGC/Phylo/mOTU2/"; #mOUT 2 dir


	#merging of per sample assignments..
	riboSummary();
	mergeMP2Table($dir_MP2);
	mergeMotu2Table($dir_mOTU2);
	unploadRawFilePostprocess();
	DiaPostProcess("",$baseOut);
	d2metaDist(
		$d2Inputs{samples}, $d2Inputs{filtered_read1}, $d2Inputs{filtered_read2},
		$d2Inputs{dependencies}, $baseOut."/d2StarComp/",
	);

	if ($MFopt{DoKraken} ){
		if( $progStats{KrakTaxFailCnts}){
			print "$progStats{KrakTaxFailCnts} samples with incomplete KrakenTax\n";} 
		else {print "All samples have assigned Kraken Taxonomy\n";
			my $mergeTblScript = getProgPaths("metPhl2Merge");#"/g/bork3/home/hildebra/bin/metaphlan2/utils/merge_metaphlan_tables.py";
			opendir D, $dir_KrakFind; my @krkF = grep {/0.*/ && -d $dir_KrakFind.$_} readdir(D); closedir D;
			my $mrgCmd = "";
			foreach my $kf (@krkF){
				#$kf =~ m/krak\.(.*)\.cnt\.tax/; my $thr = $1;# die $thr."  $kf\n";
				$mrgCmd .= "$mergeTblScript $dir_KrakFind/$kf/*.krak.txt > $dir_KrakFind/Krak.$kf.mat\n";
			}
			#die $mrgCmd."\n$dir_KrakFind\n";
			system "$mrgCmd";
		}
	}

	#open O,">$baseOut/MMPU.txt";print O $mmpuOutTab;close O;

	#final message with instructions on how to build a gene cat..
	if ($MFopt{DoAssembly}){
		my $warnMsg = "";
		$warnMsg = "(may be inaccurate due to loop2complete)" if ($runOptions{loopInitialCount});
		print "Found ". ($runReport{present_assemblies} + scalar(keys %{$runReport{empty_samples}}))." of ".
			(scalar keys %{$d2Inputs{samples}}) ." samples assembled (or ignored) and all tasks done.\n";
		print "$warnMsg\n" if ($warnMsg ne "");
		
	}
	my $totalScratchUse=0;
	foreach my $smpl (@samples){
		#input size for 1) raw 2) filtered 3) cram + assembl + consensus
		$totalScratchUse += ($map{$smpl}{inputFileSizeMB}) *3 if (exists($map{$smpl}{inputFileSizeMB}));
		$totalScratchUse += ($map{$smpl}{inputXFileSizeMB}) *3 if (exists($map{$smpl}{inputXFileSizeMB}));
	}
	print "Estimated scratch use: " . int($totalScratchUse/1024)."G\n";
	if (keys %{$runReport{empty_samples}} > 0){
		print "Found Empty/too small (<". $MFconfig{skipSmallSmplsMB} ."MB) samples (N=".
			scalar(keys %{$runReport{empty_samples}}) . "):\n".
			join(",", keys %{$runReport{empty_samples}}) ."\n\n";
	}


	if ($MFopt{DoAssembly} && $runReport{present_assemblies} > 0
		&& ($runReport{present_assemblies} + scalar(keys %{$runReport{empty_samples}}))
			== scalar(keys %{$d2Inputs{samples}})){
		my $gcScr = getProgPaths("geneCat_scr");
		my $GCsub = $baseOut."/GeneCat_pre.sh";
		my $gcmd = "";
		$gcmd .= "#creates gene catalog in the specified outdir with specified cores, attempting to reuse existing dirs (in case catalog creation failed):\n";
		my $sugGCmem = 100; 
		$sugGCmem = 200 if ($runReport{present_assemblies} > 100);
		$sugGCmem = 700 if ($runReport{present_assemblies} > 500);
		$sugGCmem = 1200 if ($runReport{present_assemblies} > 1000);
		$sugGCmem = 2500 if ($runReport{present_assemblies} > 5000);
		my $sugGCcores = 12;
		$sugGCcores = 24 if ($runReport{present_assemblies} > 100);
		$sugGCcores = 32 if ($runReport{present_assemblies} > 500);
		$sugGCcores = 48 if ($runReport{present_assemblies} > 1000);
		$sugGCcores = 72 if ($runReport{present_assemblies} > 5000);

		$gcmd .= "$gcScr -map $MFconfig{mapFile} -GCd [insert outdir] -mem $sugGCmem -cores $sugGCcores -clusterID 95 -doStrains $MFopt{DoConsSNP} -continue 1 -binSpeciesMG $MFopt{DoMetaBat2} -useCheckM2 $MFopt{useCheckM2} -useCheckM1 $MFopt{useCheckM1} -MGset GTDB \n";
		print "\n\nNext step, create a genecatalog with call to (but modify .sh first!): \nsbatch $GCsub\n";
		$QSBoptHR->{doSubmit} = 0;
		my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0; 
		my ($jobN, $tmpCmd) = qsubSystem($GCsub,$gcmd,1,"80G","GeCat","","",1,[],$QSBoptHR) ;
		$QSBoptHR->{tmpSpace} =$tmpSHDD;
		$QSBoptHR->{doSubmit} = 1;
	}
	reportSlurmJobFailures();

}

sub reportSlurmJobFailures {
	return unless (($QSBoptHR->{qmode} || '') eq 'slurm');
	my %jobs = %{$QSBoptHR->{submittedJobRecords} || {}};
	for my $job_id (keys %{$QSBoptHR->{accountingJobIds} || {}}) {
		$jobs{$job_id} ||= { requested_name => '' };
	}
	return unless %jobs;

	my $summary = slurmJobFailureSummary(\%jobs, $QSBoptHR);
	if (!defined $summary) {
		warn "Could not query Slurm accounting for the end-of-run failure summary\n";
		return;
	}
	return unless $summary->{failed};

	my %observedFailure;
	for my $entry (values %{$summary->{categories}}) {
		$observedFailure{$_} = 1 for keys %{$entry->{failures}};
	}
	my @failureColumns = grep { $observedFailure{$_} }
		qw(OOM TIMEOUT DEPENDENCY NODE_FAIL PREEMPTED CANCELLED FAILED);

	print STDERR "\nSlurm terminal failures by MATAFILER job category "
		."(N=$summary->{failed}):\n";
	print STDERR join("\t", "Job_category", "Total", @failureColumns), "\n";
	for my $category (sort keys %{$summary->{categories}}) {
		my $entry = $summary->{categories}{$category};
		print STDERR join("\t",
			$category,
			$entry->{total},
			map { $entry->{failures}{$_} || 0 } @failureColumns,
		), "\n";
	}
}





sub sampleReadSet {
	my ($sample, $phase, $replacement) = @_;
	die "Unknown read-set phase '$phase'\n"
		unless ($phase eq "raw" || $phase eq "clean");
	if (@_ >= 3) {
		$map{$sample}{reads} ||= {};
		$map{$sample}{reads}{$phase} = $replacement;
	}
	return undef unless ref($map{$sample}{reads}) eq "HASH";
	return $map{$sample}{reads}{$phase};
}

sub discoverSampleInputs {
	my ($sample, $primaryDir) = @_;
	$primaryDir = $map{$sample}{rddir} || "" unless defined($primaryDir);
	my $prefix = $map{$sample}{prefix} || "";
	my $supportSpec = $map{$sample}{SupportReads} || "";
	my $cached = $map{$sample}{inputDiscovery};
	return $cached
		if (ref($cached) eq "HASH"
			&& $cached->{primary_dir} eq $primaryDir
			&& $cached->{support_spec} eq $supportSpec);

	my $result = {
		primary_dir => $primaryDir,
		support_spec => $supportSpec,
		support_technology => "",
		primary => { read1 => [], read2 => [], single => [], bam => [], file_sizes => {} },
		support => { read1 => [], read2 => [], single => [], bam => [], file_sizes => {} },
		primary_bytes => 0,
		support_bytes => 0,
		primary_error => "",
		primary_missing_dir => 0,
		support_error => "",
	};

	if ($primaryDir ne "") {
		if (!-d $primaryDir) {
			$result->{primary_error} = "Infile dir not existing: $primaryDir\n";
			$result->{primary_missing_dir} = 1;
		} else {
			my $found = eval {
				discoverReadFiles($primaryDir, $prefix, {
					read1 => $MFconfig{readsRpairs} != 0 ? $MFconfig{rawFileSrchStr1} : "",
					read2 => $MFconfig{readsRpairs} != 0 ? $MFconfig{rawFileSrchStr2} : "",
					single => $MFconfig{rawFileSrchStrSingl},
					bam => $MFconfig{rawFileBamSrchSing},
					prefer_single => $MFconfig{prefSinglFQgreps},
				});
			};
			if (!$found) {
				$result->{primary_error} = $@ || "Could not discover primary inputs in $primaryDir\n";
			} else {
				$result->{primary} = $found;
				if ($MFconfig{doDateFileCheck}) {
					my (@datedRead1, @datedRead2);
					for (my $i = 0; $i < @{$found->{read1}}; $i++) {
						my @fileStat = stat("$primaryDir/$found->{read1}[$i]");
						next unless @fileStat;
						my $month = POSIX::strftime("%m", localtime($fileStat[9]));
						if ($month < 3) {
							push @datedRead1, $found->{read1}[$i];
							push @datedRead2, $found->{read2}[$i];
						}
					}
					$found->{read1} = \@datedRead1;
					$found->{read2} = \@datedRead2;
					$result->{primary_error} =
						"still too many file: $primaryDir\n @datedRead1\n@datedRead2\n"
						if (@datedRead1 != 1);
				}
				my %selected = map { $_ => 1 }
					(@{$found->{read1}}, @{$found->{read2}},
					 @{$found->{single}}, @{$found->{bam}});
				$result->{primary_bytes} =
					sum(0, map { $found->{file_sizes}{$_} || 0 } keys %selected);
			}
		}
	}

	if ($supportSpec ne "") {
		my ($technology, $supportPaths) = parseSupportReads($supportSpec);
		$result->{support_technology} = $technology;
		if (grep { -d $_ } @{$supportPaths}) {
			if (@{$supportPaths} != 1 || !-d $supportPaths->[0]) {
				$result->{support_error} =
					"SupportReads may specify one directory or a list of files, not a mixture: @{$supportPaths}\n";
			} else {
				my $supportDir = $supportPaths->[0];
				my $found = eval {
					discoverReadFiles($supportDir, "", {
						read1 => $MFconfig{rawFileSrchStrXtra1},
						read2 => $MFconfig{rawFileSrchStrXtra2},
						single => '\\.(?:f(?:ast)?q|f(?:ast)?a)(?:\\.(?:gz|bz2))?$',
						bam => '\\.bam$', prefer_single => 0,
					});
				};
				if (!$found) {
					$result->{support_error} =
						$@ || "Could not discover support inputs in $supportDir\n";
				} else {
					foreach my $type (qw(read1 read2 single bam)) {
						$result->{support}{$type} =
							[map { "$supportDir/$_" } @{$found->{$type}}];
					}
					$result->{support}{file_sizes} = {
						map { ("$supportDir/$_" => ($found->{file_sizes}{$_} || 0)) }
							(@{$found->{read1}}, @{$found->{read2}},
							 @{$found->{single}}, @{$found->{bam}})
					};
					$result->{support_error} =
						"Can't find supported FASTQ, FASTA, or BAM inputs in support directory $supportDir\n"
						if (!@{$found->{read1}} && !@{$found->{single}} && !@{$found->{bam}});
				}
			}
		} else {
			foreach my $supportFile (@{$supportPaths}) {
				if (!-f $supportFile) {
					$result->{support_error} = "Can't find support-read file $supportFile\n";
					last;
				}
				if ($supportFile !~ /\.(?:bam|f(?:ast)?q|f(?:ast)?a)(?:\.(?:gz|bz2))?$/i) {
					$result->{support_error} =
						"Unsupported SupportReads format for $supportFile; expected BAM, FASTQ, or FASTA, optionally gz/bz2 compressed.\n";
					last;
				}
				my $type = $supportFile =~ /\.bam$/i ? "bam" : "single";
				push @{$result->{support}{$type}}, $supportFile;
				my @fileStat = stat($supportFile);
				$result->{support}{file_sizes}{$supportFile} = $fileStat[7];
			}
		}

		if ($result->{support_error} eq "") {
			my %supportBasenames;
			foreach my $supportFile (
				@{$result->{support}{read1}}, @{$result->{support}{read2}},
				@{$result->{support}{single}}, @{$result->{support}{bam}}
			) {
				my $name = basename($supportFile);
				if ($supportBasenames{$name}++) {
					$result->{support_error} =
						"SupportReads contains more than one source named '$name'; staged filenames must be unique.\n";
					last;
				}
			}
		}
		my %selected = map { $_ => 1 }
			(@{$result->{support}{read1}}, @{$result->{support}{read2}},
			 @{$result->{support}{single}}, @{$result->{support}{bam}});
		$result->{support_bytes} =
			sum(0, map { $result->{support}{file_sizes}{$_} || 0 } keys %selected);
	}

	# Primary discovery stores names relative to its input directory, whereas
	# supplementary discovery stores resolved source paths. Compare canonical
	# identities so a symlink or repeated path cannot be counted in both scopes.
	if ($result->{primary_error} eq "" && $result->{support_error} eq "") {
		my %primaryIdentity;
		for my $file (
			@{$result->{primary}{read1}}, @{$result->{primary}{read2}},
			@{$result->{primary}{single}}, @{$result->{primary}{bam}}
		) {
			my $path = File::Spec->file_name_is_absolute($file)
				? $file : "$primaryDir/$file";
			my $identity = abs_path($path) || File::Spec->canonpath($path);
			$primaryIdentity{$identity} = $path;
		}
		my @overlap;
		for my $file (
			@{$result->{support}{read1}}, @{$result->{support}{read2}},
			@{$result->{support}{single}}, @{$result->{support}{bam}}
		) {
			my $identity = abs_path($file) || File::Spec->canonpath($file);
			push @overlap, $file if exists $primaryIdentity{$identity};
		}
		$result->{support_error} =
			"Primary and supplementary inputs resolve to the same file(s): "
			.join(", ", @overlap)."\n"
			if @overlap;
	}

	$map{$sample}{inputDiscovery} = $result;
	return $result;
}

sub populateInputSizesFast {
	my ($sample) = @_;
	my $inputs = discoverSampleInputs($sample);
	if ($MFconfig{abortOnEmptyInput}) {
		die $inputs->{primary_error} if $inputs->{primary_error} ne "";
		die $inputs->{support_error} if $inputs->{support_error} ne "";
	}
	$map{$sample}{inputFileSizeMB} = $inputs->{primary_bytes} / (1024 * 1024);
	$map{$sample}{inputXFileSizeMB} = $inputs->{support_bytes} / (1024 * 1024);
}

sub spaceInAssGrp{
#determine how much space is used in total assembly group..
	my ($curSmplX,$includeSupport) = @_;
	$includeSupport ||= 0;
	my $primaryInputSize = 0;
	my $supportInputSize = 0;
	#print "map: $map{$curSmplX}";
	my @curMems = (); 
	@curMems = @{$map{$curSmplX}{AG_members}} if (defined $map{$curSmplX}{AG_members});
	my $missedSmpls=0;my $totalSmpls = scalar(@curMems);
	foreach my $memsSmpls (@curMems){
		#print "$memsSmpls $map{$curSmpl}{inputFileSizeMB}{$memsSmpls}\n";
		if (exists($map{$memsSmpls}{inputFileSizeMB})){
			$primaryInputSize += $map{$memsSmpls}{inputFileSizeMB} ;
		} else {
			#print "Warning: file size estimation not registered for $memsSmpls\n";
			$missedSmpls++;
		}
		$supportInputSize += $map{$memsSmpls}{inputXFileSizeMB}
			if ($includeSupport && exists($map{$memsSmpls}{inputXFileSizeMB}));
	}
	if (!@curMems){
		$primaryInputSize += $map{$curSmplX}{inputFileSizeMB};
		$supportInputSize += $map{$curSmplX}{inputXFileSizeMB}
			if ($includeSupport && exists($map{$curSmplX}{inputXFileSizeMB}));
	}
	#estimate to account for missed samples..
	if ($missedSmpls){
		my $knownSmpls = $totalSmpls - $missedSmpls;
		die "Cannot estimate assembly-group input size: no member has a registered input size\n"
			if ($knownSmpls == 0);
		$primaryInputSize *= $totalSmpls / $knownSmpls;
	}
	return $primaryInputSize + $supportInputSize;
}

sub mapping_reference_matches {
	my ($stamp, $reference) = @_;
	my @stampStat = stat($stamp);
	my @referenceStat = stat($reference);
	return 0 unless (@stampStat && $stampStat[7] > 0
		&& @referenceStat && $referenceStat[7] > 0);
	open my $stampFH, '<', $stamp or return 0;
	my $line = <$stampFH>;
	close $stampFH;
	# GNU stat does not interpret \t in --format consistently across deployed
	# versions.  Accept the short-lived literal "\\t" format as well as normal
	# whitespace so mappings produced by the affected release remain usable.
	return 0 unless defined($line) && $line =~ /^(\d+)(?:\s+|\\t)(\d+)\s*$/;
	my ($recordedSize, $recordedMtime) = ($1, $2);
	return $recordedSize == $referenceStat[7]
		&& $recordedMtime == $referenceStat[9];
}

sub rmEmptySmpls{
	my @dirs = @_;
	my @outs;
	foreach my $D(@dirs){
		next if (-e "$D/SMPL.empty");
		push(@outs,$D);
	}
	return @outs;
}

sub submitGenomeBinner{
	my ($nodeSpTmpD,$metaGassembly,$MetaBat2out,$cAssGrp,$smplIDs1) = @_; #$finalCommAssDir,
	#$MFopt{DoMetaBat2} = 1: metabat2, 2: SemiBin 3: MetaDecoder, 4: GenomeFace, 5: SCGBinner
	#MetaBat2out = file with contigs per bin
	
	#die;
	#if (!-s "$MetaBat2out.cm2"){system "rm $MetaBat2out.cm2";}
	my $BinDir = $MetaBat2out;$BinDir =~ s/[^\/]+$//; #\/[^\/]+\/
	#my $finalCommAssDir = $MetaBat2out;$finalCommAssDir =~ s/\/Binning[^\/]+$//;
	my $MB2coresL = $MFopt{BinnerCores};
	my $HDDspL = $HDDspace{metabat2};
	my $inputSizeloc = spaceInAssGrp($curSmpl);

	my $nodeSpTmpD2= $nodeSpTmpD."/Bin$JNUM$cAssGrp$MFopt{DoMetaBat2}/";
	
	#if ($MFconfig{rmBinFailAssmbly}){print"Warning submitGenomeBinner::\n\nremoving $finalCommAssDir\n\n";system"rm -r $finalCommAssDir";return;}
	my @paths = @{$DOs{$cAssGrp}{wrdir}};#split /,/,$allPaths;
	@paths = rmEmptySmpls(@paths);
	die "No non-empty sample directories are available for binning assembly group $cAssGrp\n"
		unless @paths;
	my $smplIncl = scalar(@paths);
	
	my $CM1done = 0; my $CM2done = 0; my $eBinAssStat=0;
	$CM1done = 1 if (-e "$MetaBat2out.cm" );
	$eBinAssStat =1 if (-e "$MetaBat2out.assStat");
	$CM2done = 1 if (-s "$MetaBat2out.cm2" );
	
	if (binningOutputsComplete($MetaBat2out, $MFopt{useCheckM1}, $MFopt{useCheckM2})){
		return;
	}
	#die "$MetaBat2out\n";
	my $totMem = 80;#80;
	$totMem = 90 if ($CM1done || !$MFopt{useCheckM1});
	#die "$HDDspL $totMem\n";
	if ($smplIncl >20 || $inputSizeloc > 1e5 ){$HDDspL = greaterComputeSpace($HDDspL,480);$totMem = greaterComputeSpace($totMem,400);
	}elsif ($smplIncl >10 || $inputSizeloc > 5e4 ){$HDDspL = greaterComputeSpace($HDDspL,280);$totMem = greaterComputeSpace($totMem,260);
	}elsif  ($smplIncl >5 || $inputSizeloc > 2e4 ){$HDDspL = greaterComputeSpace($HDDspL,180);$totMem = greaterComputeSpace($totMem,180);
	}elsif  ($smplIncl >3 || $inputSizeloc > 2e4 ){$HDDspL = greaterComputeSpace($HDDspL,120);$totMem = greaterComputeSpace($totMem,140);
	}elsif  ($smplIncl >1 || $inputSizeloc > 1e4 ){$HDDspL = greaterComputeSpace($HDDspL,100);$totMem = greaterComputeSpace($totMem,80);
	}
	$totMem = greaterComputeSpace($totMem,220) if ($MFopt{useCheckM1});
	$totMem = greaterComputeSpace($totMem,$MFopt{BinnerMem}) if ($MFopt{BinnerMem} > 0);
	
	
	#problems with NRP cluster: some jobs require more mem and no nodes available for these..
	if ($MFopt{useBinnerScratch} ){$HDDspL=0;}
	
	#die "$totMem\n";
	my $MBcmd = "mkdir -p $nodeSpTmpD2\n";
	
	$MBcmd .= "\n\n#checking that all required mappings are done\n". checkMapsDoneSH(\@paths) ."#checks done\n\n";
	#die "$MBcmd\n\n binner\n";
	
	my @binLibraries = (
		@{getRawLibrariesAssmGrp(\%AsGrps,$cAssGrp,0)},
		@{getRawLibrariesAssmGrp(\%AsGrps,$cAssGrp,1)},
	);
	my @longLibraries = grep { $_->{is_long} } @binLibraries;
	my @shortLibraries = grep { !$_->{is_long} } @binLibraries;
	my %longTechnologies = map { ($_->{technology} || '') => 1 } @longLibraries;
	delete $longTechnologies{''};
	my $seqTec = "hiSeq";
	if (@longLibraries) {
		$seqTec = @shortLibraries || keys(%longTechnologies) != 1
			? "hybrid"
			: (keys %longTechnologies)[0];
	}
	$seqTec = "hybrid" if ($MFopt{DoAssembly} == 5);
	#die "$seqTec\n";
	

	#execute calcs later in perl script..
	my $BinnerScr = getProgPaths("Binner_scr");
	$MBcmd .= "$BinnerScr -binner $MFopt{DoMetaBat2} -binD $BinDir -smplID \"$smplIDs1\" -tmpD \"$nodeSpTmpD2\" -assmbl $metaGassembly -assmblGrp $cAssGrp -cores $MB2coresL -smplDirs " . join(",",@paths) . " -seqTec \"$seqTec\" -logDir \"$paths[-1]LOGandSUB\" ";
	$MBcmd .= "-minAssemblySizeMB $MFopt{minBinnerAssemblyMB} ";
	$MBcmd .= "-SB_env $MFopt{SB_env} " if ($MFopt{SB_env} ne "");
	$MBcmd .= ";\n";


	my $BinnerName = getBinSubdirName($MFopt{DoMetaBat2});

	# An empty assignment is a valid "no bins found" result.  Reuse any
	# published assignment for quality/statistics repair; -redoEmptyBins is the
	# explicit opt-in path that removes and recomputes an empty result.
	$MBcmd = "" if (-e $MetaBat2out);
	
	#die $MBcmd;
	
	#my $MetaBat2out = "$finalCommAssDir/Binning/MB2/$smplIDs[-1]";
	my $postCmd = "";
	if ( ( $MFopt{useCheckM1} && !$CM1done) || (!$CM2done && $MFopt{useCheckM2})  || !$eBinAssStat){
		my $mb2Qual = getProgPaths("mb2qualCheck_scr");
		$postCmd = "\n\nrm -rf $nodeSpTmpD2; mkdir -p $nodeSpTmpD2;\n";
		$postCmd .= "$mb2Qual -asm $metaGassembly -binF $MetaBat2out -tmpD $nodeSpTmpD2 -ncore $MB2coresL -checkM2 $MFopt{useCheckM2} -checkM1 $MFopt{useCheckM1} -binner $MFopt{DoMetaBat2} " ;
		$postCmd .= "\n";
	} 
	$postCmd .= "\nrm -rf $nodeSpTmpD2\n";
	#die "$MBcmd$postCmd";
	my $jobName = "Bin$JNUM";
	my $preHDDspace=$QSBoptHR->{tmpSpace};$QSBoptHR->{tmpSpace} = $HDDspL;
	#die "$AsGrps{$cAssGrp}{BinDeps}  , $cAssGrp\n";
	my $deps = ""; $deps = $AsGrps{$cAssGrp}{BinDeps} if (defined($AsGrps{$cAssGrp}{BinDeps}));
	my ($jobName2, $tmpCmd) = qsubSystem($paths[-1]."LOGandSUB/Binner$BinnerName.sh", $MBcmd.$postCmd,
			$MB2coresL, int($totMem)."G" , $jobName, $deps,"",1,[],$QSBoptHR);
	$QSBoptHR->{tmpSpace} = $preHDDspace;$QSBoptHR->{useGPUQueue} = 0;
#	die $jobName2;
	return $jobName2;
}


sub createConsSNPandSVs{
	my ($SNPinfohr) = @_;
	my %SNPinfo = %{$SNPinfohr};
	#print "SNP\nX$SNPinfo{assembly}X\n";
	my $preHDDspace=${$QSBoptHR}{tmpSpace};
	$SNPinfo{qsubDir} = "$logDir/SNP/" unless (exists($SNPinfo{qsubDir}));
	$SNPinfo{JNUM} = $JNUM;
	my $ASFS = 0; 
	if (exists($SNPinfo{assembly})){
		$ASFS = filsizeMB($SNPinfo{assembly});
	} else {die "createConsSNPandSVs:: object SNPinfo{assembly} missing.\n";}
#	${$QSBoptHR}{tmpSpace}  = int($map{$curSmpl}{inputFileSizeMB}*15/1024)+15  ."G"; #*20 for SNPconsensus_vcf #= $HDDspace{SNPcall};
	${$QSBoptHR}{tmpSpace}  = int($ASFS*400/1024)+15  ."G"; #*20 for SNPconsensus_vcf #= $HDDspace{SNPcall};

	#die "${$QSBoptHR}{tmpSpace}\n";
	$SNPinfo{QSHR} = $QSBoptHR;

	# $SNPinfo{jdeps} not used
	#system "rm -f $SNPinfo{mapD}/multi.vcf*" if ($SNPinfo{overwrite});
	unless (exists($SNPinfo{MAR})){
		my @mapping = ($SNPinfo{mapD}."/".$SNPinfo{smpl}."-smd.".$SNPinfo{bamcram});
		$SNPinfo{MAR} = \@mapping;
	}
	
	#check for supp:	my $locMapSup2Assembly =0; $locMapSup2Assembly =1 if ($MFopt{mapSupport2Assembly} && $map{$curSmpl}{"SupportReads"} ne "");
	#die "$SNPinfo{STOconSNPsupp}\n";
	if ($SNPinfo{STOconSNPsupp} ne "" && !exists($SNPinfo{MARsupp})){
		my @mapping = ($SNPinfo{mapD}."/".$SNPinfo{smpl}.".sup-smd.".$SNPinfo{bamcram});
		$SNPinfo{MARsupp} = \@mapping;
		#$finalMapDir/$SmplName.sup-smd.bam.coverage.gz
		
	} 
#	my ($ovcf,$jdep) = SNPconsensus_vcf2(\%SNPinfo);
	my ($jdep,$submissionCommands) = ("", "");
	my $runConsensus = exists($SNPinfo{callConsSNP}) || exists($SNPinfo{callConsSNPSupp})
		? (($SNPinfo{callConsSNP} || 0) || ($SNPinfo{callConsSNPSupp} || 0))
		: 1;
	($jdep,$submissionCommands) = SNPconsensus_vcf(\%SNPinfo) if $runConsensus;

	${$QSBoptHR}{tmpSpace} = $preHDDspace;
	#SNPconsensus_fasta($ovcf,\%SNPinfo,$jdep,$QSBoptHR);
	
	
	
	#2nd part: call SVs
	if (($SNPinfo{callSVs} || 0) || ($SNPinfo{callSVsSupp} || 0)){
		my ($jdep2,$qcmd2) = SVcall_vcf(\%SNPinfo);
		$jdep .= ";$jdep2" if ($jdep2 ne "");
		$submissionCommands .= $qcmd2 if (defined($qcmd2) && $qcmd2 ne "");
	}

	
	return wantarray ? ($jdep,$submissionCommands) : $jdep;
}


sub DiaPostProcess(){
	my ($cmd,$baseOut) = @_;
	#also do higher level summaries of diamond to DB mappings
	
	return if (0 || !$MFopt{DoDiamond} ); #temporary deactivated
	
	#die;
	my $mrgDiScr = getProgPaths("mrgDia_scr");
	my @DBS = split/,/,$MFopt{reqDiaDB};
	my $mapFiles = $MFconfig{mapFile};
	foreach my $DB (@DBS){
		$progStats{$DB}{DiaDBSearchCompl} =0 unless (exists($progStats{$DB}{DiaDBSearchCompl}));
		$progStats{$DB}{DiaDBSearchIncomplete}=0 unless (exists($progStats{$DB}{DiaDBSearchIncomplete}));
		my $countFile = "$baseOut/pseudoGC/FUNCT/$DB.compl";
		my $refDone = 0; $refDone = int(getFileStr($countFile)) if (-e $countFile);
		next if (!exists($progStats{$DB}{DiaDBSearchCompl}) || $progStats{$DB}{DiaDBSearchCompl} <= $refDone );
		print "$DB :: $progStats{$DB}{DiaDBSearchCompl} ($refDone previously done)\n";
		$cmd .= "$mrgDiScr $baseOut $DB $mapFiles\necho $progStats{$DB}{DiaDBSearchCompl} > $countFile\n" if ($progStats{$DB}{DiaDBSearchCompl}>= 1 );
		#`echo $progStats{$DB}{DiaDBSearchCompl} > $countFile`;
		#print "$DB: complete $progStats{$DB}{DiaDBSearchCompl} >= incomplete $progStats{$DB}{DiaDBSearchIncomplete}\n"
	}
	$cmd .= $AsGrps{global}{DiamCln};
	print "\nMerg diamond:$cmd\n\n";
	systemW $cmd;

}
sub mergeMotu2Table($){
	my ($inD) = @_;
	return unless ($MFopt{DoMOTU2});
	if ($progStats{mOTU2FailCnts}){print "$progStats{mOTU2FailCnts} / ". ($progStats{mOTU2ComplCnts}+$progStats{mOTU2FailCnts}) ." samples with incomplete mOTU2 assignments\n"; return;}
	print "\nAll samples ($progStats{mOTU2ComplCnts}) have mOTU2 assignments.\n";
	my $outD = $inD;
	$outD =~ s/[^\/]+\/?$//;

	my $m2mrgSto = "$outD/m2.Smpl.cnts.stone";
	my $oldM2Cnt = 0;
	if (-e $m2mrgSto){$oldM2Cnt = `cat $m2mrgSto`; chomp $oldM2Cnt;}
	#print "$m2mrgSto $oldM2Cnt $progStats{mOTU2ComplCnts}\n";
	if ($oldM2Cnt >= $progStats{mOTU2ComplCnts} && -e "$outD/m2.class.txt"){return;}

	
	#if (-e "$outD/m2.$taxLvlN[0].txt"){return;}
	
	print "Merging motu2 files into tax tables..\n";
	my $mrgMOTU2 = getProgPaths("mrgMOUT2_scr");
	my $mrgCmd = "$mrgMOTU2 $inD $progStats{mOTU2ComplCnts}\n";
	#die "$mrgCmd\n";
	my ($jobN, $tmpCmd) = qsubSystem($baseOut."/LOGandSUB/MOTU2merg.sh",$mrgCmd,1,"80G","MOTU2mrg","","",1,[],$QSBoptHR) ;
}



#calculates the hmm based freq estimates and divergence from these -> used for dist matrix
sub d2metaDist{
	my ($hrSmpls,$arPaths1,undef,$deps,$outPath) = @_;
	if(!$MFopt{DoCalcD2s}){return;}
	my @paths = @{$arPaths1}; my @Smpls = keys (%{$hrSmpls});
	if (@paths < 1){print "Not enough samples for d2s!\n";return;}
	my $d2metaBin = getProgPaths("d2meta");#"/g/bork3/home/hildebra/bin/d2Meta/d2Meta/d2Meta.out";
	print "Calculating kmer distances for ".@paths." samples\n";
	system "mkdir -p $outPath/LOGandSUB";
	my $d2MK = 6;
	my $smplFile = "$outPath/mapd2s.txt";
	open O,">$smplFile";
	for (my $i=0;$i<@paths;$i++){
		print O "$paths[$i] $Smpls[$i]\n";
	}
	close O;
	my $cmd = "";
	$cmd .= "cd $outPath\n";
	$cmd .= "$d2metaBin $d2MK $smplFile Q\n"; #or Q for fastQ
	$cmd .= "touch $outPath/d2meta.stone\n";
	my $jobN = ""; my $tmpCmd;
	unless (-e "$outPath/d2meta.stone"){
		$jobN = "_d2met";
		($jobN, $tmpCmd) = qsubSystem($outPath."/LOGandSUB/d2Met.sh",$cmd,1,"80G",$jobN,"$deps","",1, $QSBoptHR->{General_Hosts},$QSBoptHR) ;
	}
	return $jobN;
}


sub postSubmQsub {
	my ($outf, $commands, $dependencies, $options) = @_;
	return "" if ($commands eq "" || !$runOptions{submit});
	$options ||= {};
	die 'postSubmQsub options must be a hash reference'
		unless ref($options) eq 'HASH';
	my @submitted;
	my @augmented_commands;
	for my $command (grep { /\S/ } split /\r?\n/, $commands) {
		my ($script_path) = $command =~ /(?:^|\s|<)(\S+\.sh)(?:\s|$)/;
		die "Deferred submission does not identify a job script: $command\n"
			unless (defined $script_path && -e $script_path);
		open my $script_fh, '<', $script_path
			or die "Cannot read deferred job script $script_path: $!\n";
		my $script = do { local $/; <$script_fh> };
		close $script_fh;
		# Deferred batches are independent by default and share only their producer
		# barrier. A genuinely multi-stage batch must request predecessor chaining.
		my $command_dependencies = deferred_command_dependencies(
			dependencies => $dependencies, submitted => \@submitted,
			chain_previous => $options->{chain_previous},
		);
		if (submissionDependencyDeferred($command_dependencies)) {
			push @submitted, deferredSubmissionDependency();
			last;
		}
		my $augmented = augment_deferred_submission(
			qmode => $QSBoptHR->{qmode}, command => $command, script => $script,
			dependencies => $command_dependencies, run_tag => $QSBoptHR->{rTag},
		);
		if ($augmented->{script} ne $script) {
			open my $script_out, '>', $script_path
				or die "Cannot update deferred job script $script_path: $!\n";
			print {$script_out} $augmented->{script};
			close $script_out or die "Cannot close deferred job script $script_path: $!\n";
		}
		push @augmented_commands, $augmented->{command};
		my $submitted_before = scalar @submitted;
		my $scheduler_job_id = "";
		my $capacityAvailable = qsubSystemWaitMaxJobs(
			$MFconfig{checkMaxNumJobs}, $MFconfig{killDepNever}, $QSBoptHR,
		);
		unless ($capacityAvailable) {
			push @submitted, deferredSubmissionDependency();
			last;
		}
		my ($output, $status) = $QSBoptHR->{qmode} eq 'slurm'
			? submitSlurmWithDependencyRecovery(
				$augmented->{command}, $script_path, $QSBoptHR,
			)
			: do {
				my $scheduler_output = `$augmented->{command} 2>&1`;
				($scheduler_output, $?);
			};
		if ($status != 0) {
			my $message = "Deferred job submission failed: "
				."$augmented->{command}\n$output";
			push @submitted, handleSubmissionFailure($QSBoptHR, $message);
			last;
		}
		if ($QSBoptHR->{qmode} eq 'slurm') {
			if ($output =~ /^Submitted batch job (\d+)\s*$/m) {
				$scheduler_job_id = $1;
				push @submitted, $QSBoptHR->{rTag}.$scheduler_job_id;
			}
		} elsif ($QSBoptHR->{qmode} eq 'sge') {
			if ($output =~ /\bYour job(?:-array)?\s+(\d+)\b/) {
				$scheduler_job_id = $1;
				push @submitted, $QSBoptHR->{rTag}.$scheduler_job_id;
			}
		} elsif ($QSBoptHR->{qmode} eq 'lsf') {
			if ($output =~ /\bJob <(\d+)>/) {
				$scheduler_job_id = $1;
				push @submitted, $QSBoptHR->{rTag}.$scheduler_job_id;
			}
		}
		die "Could not parse deferred scheduler job id: $output\n"
			if ($QSBoptHR->{qmode} ne 'bash' && @submitted == $submitted_before);
		$QSBoptHR->{submittedJobs} = 0 unless defined $QSBoptHR->{submittedJobs};
		$QSBoptHR->{submittedJobs}++;
		if ($QSBoptHR->{qmode} eq 'slurm' && $scheduler_job_id ne "") {
			$QSBoptHR->{slurmDependencySubmittedAt} ||= {};
			$QSBoptHR->{slurmDependencySubmittedAt}{$scheduler_job_id} = time;
			my ($slurmJobName) = $augmented->{script} =~ /^#SBATCH\s+-J\s+(\S+)/m;
			$QSBoptHR->{submittedJobRecords}{$scheduler_job_id} = {
				requested_name => $slurmJobName || '',
			};
		}
		recordSampleLockJobs(
			$QSBoptHR->{LOCKfile}, [$scheduler_job_id], $QSBoptHR,
		) if $scheduler_job_id ne "";
	}
	open my $audit_fh, '>', $outf or die "Cannot write deferred submission audit $outf: $!\n";
	print {$audit_fh} join("\n", @augmented_commands), "\n";
	close $audit_fh or die "Cannot close deferred submission audit $outf: $!\n";
	return normalise_job_dependencies(\@submitted);
}





sub RiboMeta($ $ $ $){
	my ($calcRibofind,$calcRiboAssign,$curOutDir,$SmplName) = @_;
	my $dir_RibFind = $baseOut.$preDIRs{dir2RiboF}; #ribofind dir

	if ($calcRibofind || $calcRiboAssign){
		if ($MFopt{RedoRiboThatFailed} ){
			system "rm -r $curOutDir/ribos/";
			$calcRibofind = 1;
		}
		$progStats{riboFindFailCnts} ++ ;
	} elsif ($MFopt{DoRibofind} && !$calcRiboAssign) { #all done, copy files to central dir for postprocessing..
		my @RFtags = ("SSU","LSU");#"ITS",
		foreach my $RFtag (@RFtags){
			system "mkdir -p $dir_RibFind/$RFtag/" unless (-d "$dir_RibFind/$RFtag/"); #system "mkdir -p $dir_RibFind/SSU/" unless (-d "$dir_RibFind/SSU/"); system "mkdir -p $dir_RibFind/LSU/" unless (-d "$dir_RibFind/LSU/");
			my $fromCp = "$curOutDir/ribos/ltsLCA/${RFtag}riboRun_bl.hiera.txt"; my $toCpy = "$dir_RibFind/$RFtag/$SmplName.$RFtag.hiera.txt";
			my @sourceStat = stat($fromCp);
			my @sourceGzipStat = stat("$fromCp.gz");
			my @destinationStat = stat($toCpy);
			my @destinationGzipStat = stat("$toCpy.gz");
			if ($MFopt{checkRiboNonEmpty}){
				#pretty hard check
				my $numLines=0;
				if (-e "$fromCp.gz"){$numLines = `zcat $fromCp.gz | wc -l`;
				} else {$numLines = `wc -l $fromCp`;} $numLines =~ /(\d+)/; $numLines=$1;
				#die $numLines."\n";
				if ($numLines<=1){$calcRiboAssign=1;$calcRibofind=1;
					system "rm -r $curOutDir/ribos//ltsLCA $curOutDir/ribos/*.sto ";last;
				}
			}
			my $sourceSize = @sourceStat ? $sourceStat[7]
				: (@sourceGzipStat ? $sourceGzipStat[7] : undef);
			my $destinationSize = @destinationStat ? $destinationStat[7]
				: (@destinationGzipStat ? $destinationGzipStat[7] : undef);
			if (!defined($destinationSize)
					|| (defined($sourceSize) && $sourceSize != $destinationSize)){
				unlink "$toCpy" if @destinationStat;
				if (@sourceGzipStat){
					#system "zcat $fromCp.gz > $toCpy" ;
					system "rm -f $toCpy.gz;ln -s $fromCp.gz $toCpy.gz";
				} elsif (@sourceStat) {
					system "gzip $fromCp";
					system "rm -f $toCpy.gz;ln -s $fromCp.gz $toCpy.gz";
				} elsif($RFtag eq "SSU") {#just redo.. SSU is only essential thing
					system "rm -rf $curOutDir/ribos\n"; 
					$calcRibofind = 1; $calcRiboAssign=1;
					last;
				}
			}
			#system "gzip $fromCp" unless (-e "$fromCp.gz");
			system "rm -f $curOutDir/ribos/ltsLCA/inter${RFtag}riboRun_bl.fna" if (-e "$curOutDir/ribos/ltsLCA/inter${RFtag}riboRun_bl.fna");
		}
		$progStats{riboFindComplCnts} ++; #completed already
	} 
}


sub prepareDiamondRerun($){
	my ($curOutDir) = @_;
	return unless ($MFopt{DoDiamond});
	my @alldbs = split /,/,$MFopt{reqDiaDB};
	my $all_requested = @alldbs == $MFopt{maxReqDiaDB};

	if ($MFopt{rewriteDiamond} && $all_requested) {
		if (-d "$curOutDir/diamond/") {
			system('rm', '-rf', '--', "$curOutDir/diamond/") == 0
				or die "Failed to remove $curOutDir/diamond/ for Diamond rebuild\n";
		}
		return;
	}
	if ($MFopt{redoDiamondParse} && $all_requested) {
		for my $target (glob("$curOutDir/diamond/CNT*")) {
			system('rm', '-rf', '--', $target) == 0
				or die "Failed to remove Diamond parse target $target\n";
		}
	}
	return unless ($MFopt{redoDiamondParse});
	my $secCogBin = getProgPaths("secCogBin_scr");
	foreach my $term (@alldbs){
		unlink glob("$curOutDir/diamond/dia.$term.blast.*.stone");
		# Preserve the legacy parser input until its XX layout is replaced by a
		# declared stage artifact.
		system("$secCogBin -i $curOutDir/diamond/XX -DB $term -eval $MFopt{diaEVal} -mode 4") == 0
			or die "Diamond reparsing failed for database $term\n";
		if ($term eq "ABR" && -d "$curOutDir/diamond/ABR/") {
			system('rm', '-rf', '--', "$curOutDir/diamond/ABR/") == 0
				or die "Failed to remove stale Diamond ABR output\n";
		}
	}
}

sub IsDiaRunFinished($){
	my ($curOutDir) = @_;
	my @alldbs = split /,/,$MFopt{reqDiaDB};
	if (!$MFopt{DoDiamond}){return (0,0);}
	my $cD = 0; my $pD = 0; #dia_calc, dia_parse
	if ($MFopt{rewriteDiamond}){$MFopt{redoDiamondParse} = 1;}
	foreach my $term (@alldbs){
		#print $term."   $cD, $pD\n";
		if ($MFopt{rewriteDiamond} ){$pD=1; $cD=1;}
		#die "$MFopt{rewriteDiamond}\n";
#print "$curOutDir/diamond/dia.$term.blast.gz\n";
		#die "$curOutDir/diamond/dia.$term.blast.gz\n$curOutDir/diamond/dia.$term.blast.srt.gz";
		#|| !-e "$curOutDir/diamond/dia.$term.blast.srt.gz"
		if (!$cD && (!-e "$curOutDir/diamond/dia.$term.blast.gz" && !-e "$curOutDir/diamond/dia.$term.blast.srt.gz" )){$cD = 1; }#system "rm $curOutDir/diamond/dia.$term.blas*.gz";}
		$pD = 1 if (!-e "$curOutDir/diamond/dia.$term.blast.srt.gz.stone");#  <- last version always requires .srt.gz
		#print "$cD, $pD  $curOutDir/diamond/dia.$term.blast.gz\n";
	}
	#die "$cD, $pD\n";
	$cD = 0 if ($pD==0);
	
	if ($MFopt{rewriteAllIfAnyDiamond} && ($cD  || $pD)){ #just delete everything..
		system "rm -r $curOutDir/diamond/" if (-d "$curOutDir/diamond/"); $cD=1; $pD=1;
	}
	if (!$cD && !$pD){
		foreach my $curDB (@alldbs){$progStats{$curDB}{DiaDBSearchCompl}++;}
	}
	
	return ($cD,$pD);
}


#called in case sample "is empty", reduce some counters
sub reduceProgStats{
	if ($MFopt{DoMetaPhlan}){
		$progStats{metaPhl2FailCnts}--;
	}
	if ($MFopt{DoMOTU2}){
		$progStats{mOTU2FailCnts}--;
	}
	if ($MFopt{DoTaxaTarget}){
		$progStats{taxTarFailCnts}--;
	}
	if ($MFopt{DoRibofind}){
		$progStats{riboFindFailCnts} -- ;
	}
}

#check if programes have finished, that rely only on raw reads
sub checkRawProgsFin{
	my ($curOutDir,$SmplName) = @_;
	my ($calcDiamond,$calcDiaParse) = IsDiaRunFinished($curOutDir);
	
	my $calcRibofind = 0; my $calcRiboAssign = 0;my $calcGenoSize=0; 
	my $calcMetaPhlan=0;	my $calcMOTU2=0;	my $calcTaxaTar = 0;
	my $calcKraken =0;
	
	
	#Kraken .. totally outdated way of doing things..
	
	my $KrakenOD = $curOutDir."Tax/kraken/$MFopt{globalKraTaxkDB}/";
	$calcKraken = 1 if ($MFopt{DoKraken} && (!-d $KrakenOD || !-e "$KrakenOD/krakDone.sto"));
	$progStats{KrakTaxFailCnts}++ if ($calcKraken && $MFopt{DoKraken});

	
	$calcGenoSize=1 if ($MFopt{DoGenoSizeEst} && 	!-e "$curOutDir/MicroCens/MC.0.result");
	$calcRibofind = 1 if ($MFopt{DoRibofind} && (!-e "$curOutDir/ribos//SSU_pull.sto"|| !-e "$curOutDir/ribos//LSU_pull.sto" || ($MFopt{doRiboAssembl} && !-e "$curOutDir/ribos/Ass/allAss.sto" ))); #!-e "$curOutDir/ribos//ITS_pull.sto"|| 
	$calcRiboAssign = 1 if ($MFopt{DoRibofind} && ( !-e "$curOutDir/ribos//ltsLCA/Assigned.sto"  || !-e "$curOutDir/ribos//ltsLCA/SSU_ass.sto") );		#!-e "$curOutDir/ribos//ltsLCA/ITS_ass.sto"||  #ITS no longer required.. unreliable imo
	RiboMeta($calcRibofind,$calcRiboAssign,$curOutDir,$SmplName);

	#die $dir_MP2."$SmplName.MP2.sto";

	if ($MFopt{DoMetaPhlan}){
		my $dir_MP2 = $baseOut.$preDIRs{dir2MePhl};#"pseudoGC/Phylo/MP2/"; #metaphlan 2 dir
		if (!-e $dir_MP2."$SmplName.MP2.sto"){
			$calcMetaPhlan=1 ;
			$progStats{metaPhl2FailCnts}++;
		} else { $progStats{metaPhl2ComplCnts}++;}
	}
	if ($MFopt{DoMOTU2}){
		my $dir_mOTU2 = $baseOut."pseudoGC/Phylo/mOTU2/"; #mOUT 2 dir
		if (!-e $dir_mOTU2."$SmplName.Motu2.sto"){
			$calcMOTU2=1;
			$progStats{mOTU2FailCnts}++;
		} else { $progStats{mOTU2ComplCnts}++;}
	}
	if ($MFopt{DoTaxaTarget}){
		my $dir_TaxTar = $baseOut."pseudoGC/Phylo/TaxaTarget/"; #taxaTar dir
		if (!-e $dir_TaxTar."$SmplName.TaxTar.sto"){
			$calcTaxaTar=1 ;
			$progStats{taxTarFailCnts}++;
		} else { $progStats{taxTarComplCnts}++;}
	}
	return ($calcKraken,$calcDiamond,$calcDiaParse,$calcRibofind,$calcRiboAssign,$calcGenoSize,
			$calcMetaPhlan, $calcMOTU2,$calcTaxaTar) ;
}

sub riboSummary{
	return if (!$MFopt{DoRibofind} );
	if( $progStats{riboFindFailCnts}>0){
		print "$progStats{riboFindFailCnts} / $progStats{riboFindComplCnts} samples with incomplete RiboFind\n";
		return;
	} 
	my $dir_RibFind = $baseOut.$preDIRs{dir2RiboF};#"pseudoGC/Phylo/RiboFind/"; #ribofinder dir
	my $prevItems = 0;
	if ( -e "$dir_RibFind/SSU.cnt.stone"){
		my $tmp = `cat $dir_RibFind/SSU.cnt.stone`;
		chomp $tmp;  $prevItems = $tmp+0;
	}
	my $mergeMiTagScript = getProgPaths("mrgMiTag_scr");#"/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/secScripts/miTagTaxTable.pl";
	print "All samples ($progStats{riboFindComplCnts}) have assigned RiboFinds\n";
	my $ITSpres=0;
	my @lvls = ("domain","phylum","class","order","family","genus","species","hit2db");
	#unless (!-d "$dir_RibFind/ITS/"){$mrgCmd .= "$mergeMiTagScript ".join(",",@lvls)." $dir_RibFind/ITS.miTag $dir_RibFind/ITS/ \n"; $ITSpres=1;}
	my $mrgCmd = "$mergeMiTagScript ".join(",",@lvls)." $dir_RibFind/SSU.miTag $dir_RibFind/SSU/ \n" unless (!-d "$dir_RibFind/SSU/");
	$mrgCmd .= "echo \"$progStats{riboFindComplCnts}\" > $dir_RibFind/SSU.cnt.stone\n";
	my $mrgCmd2 = "$mergeMiTagScript ".join(",",@lvls)." $dir_RibFind/LSU.miTag $dir_RibFind/LSU/ \n" unless (!-d "$dir_RibFind/LSU/");
	$mrgCmd2 .= "echo \"$progStats{riboFindComplCnts}\" > $dir_RibFind/LSU.cnt.stone\n";
	#die $mrgCmd."\n";
	my $of_exist = 1;
	foreach my $lvl (@lvls){ 
		$of_exist=0 unless ((!$ITSpres || -e "$dir_RibFind/ITS.miTag.$lvl.txt") && -e "$dir_RibFind/LSU.miTag.$lvl.txt" && -e "$dir_RibFind/SSU.miTag.$lvl.txt"); 
		#die "$dir_RibFind/ITS.miTag.$lvl.txt\n$dir_RibFind/LSU.miTag.$lvl.txt\n$dir_RibFind/SSU.miTag.$lvl.txt\n" if (!$of_exist);
	}
	if ($prevItems < $progStats{riboFindComplCnts}){
		print "Redoing ribo tables, as more samples currently available ($progStats{riboFindComplCnts}, prev:$prevItems)\n";
		$of_exist = 0;
	}
	#die;
	#system $mrgCmd."\n" ; die;
	if ($of_exist == 0){
		my ($jobN, $tmpCmd) = qsubSystem($baseOut."/LOGandSUB/SSUmerge.sh",$mrgCmd,1,"80G","SSUmrg","","",1,[],$QSBoptHR) ;
		($jobN, $tmpCmd) = qsubSystem($baseOut."/LOGandSUB/LSUmerge.sh",$mrgCmd2,1,"80G","LSUmrg","","",1,[],$QSBoptHR) ;
	}
	#system "$mrgCmd";
}

# I've removed sortmerna (and ITS-related sortmerna) from this, since 
# the index is now only read once into memory. If this doesn't work, we will need to
# add it back in. 
sub detectRibo(){
	my ( $tmpP,$outP,$jobd,$SMPN,$glbTmpDDB) = @_;
	my $cleanSeqSetHR = sampleReadSet($curSmpl, "clean");
	my $libraries = readLibrariesByScope($cleanSeqSetHR, 'primary', 1, $curSmpl);
	#print "DEP: $jobd\n";
	my $numCore = 12;
	my $numCore2 = 12;
	my $cLSUSSUscript = getProgPaths("cLSUSSU_scr");#"perl /g/bork3/home/hildebra/dev/Perl/16Stools/catchLSUSSU.pl";
	#my $lambdaIdxBin = getProgPaths("lambdaIdx");
	my $lambdaBin = getProgPaths("lambda");#"/g/bork3/home/hildebra/dev/lotus//bin//lambda/lambda";
	
	
	
	my $pairs = libraryPairs($libraries);
	my @re1 = map { $_->{files}{r1} } @{$pairs};
	my @re2 = map { $_->{files}{r2} } @{$pairs};
	my @singl = @{libraryFiles($libraries, 'single')};
	#print "ri"; 
	if (@re1 > 1 || @re2 > 1){
		#print"\nWARNING::\nOnly the first read file will be searched for ribosomes\n";
		#first create new tmp file
	}
	
	#copy DB to server
	my $DBrna = "$glbTmpDDB/rnaDB/";
	my $DBrna2 = "$glbTmpDDB/LCADB/";
	if ($MFopt{globalRiboDependence}->{DBcp} eq "" ){

		my $DBcmd = "";
		my $DBcores = 1;
		$MFopt{globalRiboDependence}->{DBcp}="alreadyCopied";
		#my $ITSfilePref = $1;


		#die @DBs."@DBs\n";

		###  Uncomment this section if you do actually end up needing sortmerna index ###
		#my @DBs = split(/,/,getProgPaths("srtMRNA_DBs"));
		#my @DBsIdx = @DBs;my @DBsTestIdx = @DBs; my $filesCopied = 1;
		#my $srtMRNA_path = getProgPaths("srtMRNA_path");
		#if ( !-d $DBrna ){
		#	$DBcmd .= "mkdir -p $DBrna\n";
		#}
		#for (my $ii=0;$ii<@DBsIdx;$ii++){
		#	$DBsIdx[$ii] =~ s/\.fasta$/\.idx/;
		#	$DBsTestIdx[$ii] =~ s/\.fasta$/\.idx\.kmer_0\.dat/;
		#	if ( -e "$DBrna//$DBs[$ii]"  && -e "$DBrna//$DBsTestIdx[$ii]"  ){
		#		next;
		#	}
		#	die "\nCould not find expected sortmerna file:\n$srtMRNA_path/rRNA_databases/$DBs[$ii]\n" if ( !-e "$srtMRNA_path/rRNA_databases/$DBs[$ii]"  );
		#	if ( !-e "$srtMRNA_path/rRNA_databases/$DBsTestIdx[$ii]"  ){
		#		$DBcmd .= "$srtMRNA_path./indexdb_rna --ref $srtMRNA_path/rRNA_databases/$DBs[$ii],$srtMRNA_path/rRNA_databases/$DBsIdx[$ii]\n";
		#	}
		#	$DBcmd .= "\ncp $srtMRNA_path/rRNA_databases/$DBs[$ii] $srtMRNA_path/rRNA_databases/$DBsIdx[$ii]* $DBrna\n";
		#}
		
		#my $ITSDBfa = getProgPaths("ITSdbFA",0);
		#if ($ITSDBfa ne ""){ #only do if not empty.. otherwise ignore (not required)
		#	my $ITSDBpref = $ITSDBfa;$ITSDBpref =~ s/\.fa.*$//;
		#	my $ITSDBidx = $ITSDBfa; $ITSDBidx =~ s/\.fa.*$/\.idx/;
		#	$ITSDBpref =~ m/\/([^\/]+)$/;
		#	$ITSDBpref=~ m/(^.*)\/[^\/]+/;
			
		#	if (!-e "$ITSDBpref.idx.kmer_0.dat"){ #ITS DBs
				
		#		if (!-e "$ITSDBfa"){
		#			print "Missing $ITSDBfa  ITS DB file!\n"; exit(32);
		#		}
		#		if (!-e "$ITSDBidx.kmer_0.dat"){
		#			$DBcmd .= "\n$srtMRNA_path./indexdb_rna --ref $ITSDBfa,$ITSDBidx\n";
		#		}
		#		$DBcmd .= "\ncp ${ITSDBpref}* $DBrna\n";
				#has to be noted that this doesn't need to happen again
				#print "ribo DB already present\n";
				
				#$DBcmd = "";
		#	} 
		#}
		### End of sortmerna index section ###


		#and get flash DBs as well over to that dir
		my @DBn = ("LSUdbFA","LSUtax","SSUdbFA","SSUtax");#,"PR2dbFA","PR2tax"); #"ITSdbFA","ITStax",
		my $LCAar = getProgPaths(\@DBn,0);
		my @LCAdbs = @{$LCAar}; 
#		die "@LCAdbs\n".@LCAdbs."\n";
		#check first if the lambda DB was already built
		my $doCopyDBtoScratch = 0;
		$DBcmd .= "\nmkdir -p $DBrna2\n" if (!-d $DBrna2);
		for (my $kk=0;$kk<@LCAdbs; $kk+=2){
			my $DB = $LCAdbs[$kk];
			#first test if already copied to scratch..
			$LCAdbs[$kk] =~ m/\/([^\/]+)$/;
			my $SLVtestNme = $1;
			if (!-e $DBrna2."/$SLVtestNme.lambda/index.lf.drp"|| !-e $DBrna2."/$SLVtestNme" ){
				$doCopyDBtoScratch = 1;
			} else {next;}
			#print "$DB\n";
			unless (-e $DB){print "$DBn[$kk] not found, won't use it\n";next;}
			die "wrong DB checked for index built:$DB\n" if ($DB =~ m/\.tax$/);
#			if (!-f $DB.".dna5.fm.sa.val"  ) { #old lambda 1.0x style
			if (!-f "$DB.lba.gz" ){#new lambda3 style  # !-f $DB.".lambda/index.lf.drp"  ) { #new 1.9x style
				print "Building LAMBDA index anew (may take up to an hour)..\n";
				$DBcores = 12;
				#$DBcmd .= "\n$lambdaIdxBin -p blastn -t $DBcores -d $DB\n";
				$DBcmd .= "$lambdaBin mkindexn -t $DBcores -d $DB;\n";
				$DBcmd .= "$pigzBin -p $DBcores $DB.lba;\n";

				#die "$DBcmd\n";
			}
			if ($doCopyDBtoScratch){
				#my $jstr = join("* ",@LCAdbs);$jstr =~ s/\s\*//g;	$jstr =~ s/^\*//g;
				$DBcmd.= "\ncp -r $LCAdbs[$kk]* ".$LCAdbs[$kk+1]."* $DBrna2\n";
			}
		}
		#die "$DBcmd\n$DBrna2\n";
		#die "$SLVlsuNme\n";
		#print "$DBrna2/SLV_128_LSU.tax || !-e $DBrna2/gb203_pr2_all_10_28_99p.fasta  $DBrna2/ITS_comb.fa.dna5.fm.sa.val\n";
		
		my $jN = "_RRDB$JNUM"; my $tmpCmd;
		#die $DBcmd."$DBrna/ITS_comb.idx.kmer_0.dat";
		if ($DBcmd ne ""){
			#die $DBcmd."\n";
			my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0; 
			($jN, $tmpCmd) = qsubSystem($logDir."RiboDBprep.sh",$DBcmd,$DBcores,"10G",$jN,"","",1,$QSBoptHR->{General_Hosts},$QSBoptHR);
			$QSBoptHR->{tmpSpace} =$tmpSHDD;
			$MFopt{globalRiboDependence}->{DBcp} = $jN;
		}
	}
	#die;
	
	#first detect LSU/SSU/ITS in metag
	my $tmpDY = "$tmpP/SMRNA/";
	my $cmd = "";
	#die "$MFconfig{readsRpairs}\n";
	my $readConfig = 1;
	if (@re1 > 0){
		$cmd .= "\n$cLSUSSUscript -R1 '".join(",",@re1)."' -R2 '". join(",",@re2)."' ";
		if (@singl>0){
			$cmd .= " -RS '".join(",",@singl) . "' "; #tmpP < scratch too slow
		} else {
			$cmd .=" -RS '-1' "; #tmpP < scratch too slow
		}
		$cmd .= "-tmpDir $tmpDY -alignDir $outP -cores $numCore -smplID $SMPN -assmblRibos $MFopt{doRiboAssembl} \n";#"$tmpDY $outP $numCore $SMPN $MFopt{doRiboAssembl} $DBrna\n\n";
	} else {
		$readConfig = 0; 
		$cmd .= "\n$cLSUSSUscript -R1 '-1' -R2 '-1' -RS '".join(",",@singl)."' -tmpDir $tmpDY -alignDir $outP -cores $numCore -smplID $SMPN -assmblRibos $MFopt{doRiboAssembl} \n\n"; #tmpP < scratch too slow
	}
	my $sto1 = "$outP/RibFnd.sto";my $stoLCAL = "$outP//ltsLCA/LSU_ass.sto";	my $stoLCAS = "$outP//ltsLCA/SSU_ass.sto";
	$cmd .= "touch $sto1\n" unless ($cmd eq "");
	#this part assigns tax
	my $tmpDX = "$tmpP/LCA/"; 
	my $numCoreL2 = int($numCore2) ;  
	#my $readConfig =$MFconfig{readsRpairs};
	$readConfig=2 if (@singl>0 && $readConfig);
	#ver 2
	#my $cmd2 = "$lotusLCA_cLSU $outP $SMPN $numCoreL2 $DBrna2 $tmpDX $readConfig\n\n";
	#ver 3
	my $cfgstr = "";
	$cfgstr = "-config $MFconfig{configFile} " if ($MFconfig{configFile} ne "" );
	my $lotusLCA_cLSU = getProgPaths("lotusLCA_cLSU_scr");#"perl /g/bork3/home/hildebra/dev/Perl/16Stools/lotus_LCA_blast2.pl";
	my $cmd2 = "";
	$cmd2 .= "rm -rf $tmpDX/*;\nmkdir -p $tmpDX\n\n";
	$cmd2 .= "$lotusLCA_cLSU -dir $outP -smplID $SMPN -cores $numCoreL2 -DBdir $DBrna2 -tmpD $tmpDX -keepReads $MFopt{riboStoreRds} -pairedRds $readConfig $cfgstr -maxReadNum $MFopt{riboLCAmaxRds} -simMode 2 \n\n";
	
	
	#die $cmd2."\n";
	my $jobName="";
	my $Scmd = "";
	my $allLCAstones = 0; 
	$allLCAstones = 1 if ( -e $stoLCAL && -e $stoLCAS );#&& -e "$outP//ltsLCA/ITS_ass.sto");

	if (-d $outP  && -e "$outP/SSU_pull.sto" && -e "$outP/LSU_pull.sto" &&  #&& -e "$outP/ITS_pull.sto"
		-s "$outP//ltsLCA/LSUriboRun_bl.hiera.txt" && $allLCAstones &&
		($MFopt{doRiboAssembl} && -e $outP."/Ass/allAss.sto" ) ){
		#really everything done
		$jobName = $jobd;
	} else {
		$jobd.=";".$MFopt{globalRiboDependence}->{DBcp} unless ($MFopt{globalRiboDependence}->{DBcp} eq "alreadyCopied");
		my $tmpCmd; my $mem = "3G";
		#better to double check calcRiboFind
		my $calcRiboFind=0;$calcRiboFind = 1 if( !-e "$outP/SSU_pull.sto"|| !-e "$outP/LSU_pull.sto" || ($MFopt{doRiboAssembl} && !-e $outP."/Ass/allAss.sto"));
		if (  $calcRiboFind ){ #!-e "$outP/ITS_pull.sto"||
			$jobName = "_RF$JNUM"; 
			#die "RIBOFIND\n$outP/SSU_pull.sto\n"; 
			my $tmpSHDD = $QSBoptHR->{tmpSpace};
			my $curSHFF = int($map{$SMPN}{inputFileSizeMB}/1024*17)+5;
			my $predefSHDD = $HDDspace{Ribos}; $predefSHDD =~ s/G$//;
			# catchLSUSSU writes extracted reads and SortMeRNA work files below
			# the node-local directory. Preserve the configured floor for small
			# samples and scale above it for large inputs.
			$curSHFF = $predefSHDD if ($curSHFF < $predefSHDD);
			$QSBoptHR->{tmpSpace}= $curSHFF . "G";
			
			($jobName, $tmpCmd) = qsubSystem($logDir."RiboFinder.sh",$cmd,$numCore,$mem,$jobName,$jobd,"",1,[],$QSBoptHR);
			$QSBoptHR->{tmpSpace} = $tmpSHDD; 
		} else {
			$jobName = $jobd;
		}
		if (!-e "$outP//ltsLCA/LSUriboRun_bl.hiera.txt" || !-e "$outP//ltsLCA/SSUriboRun_bl.hiera.txt" || !$allLCAstones ){ #|| !-e "$outP//ltsLCA/ITSriboRun_bl.hiera.txt" 
			$jobd=$jobName; $mem="4G";
			$jobd .= ";".$MFopt{globalRiboDependence}->{DBcp} unless ($MFopt{globalRiboDependence}->{DBcp} eq "alreadyCopied");
			$QSBoptHR->{useLongQueue} = 0;
			my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = $HDDspace{Ribos};
			($jobName, $tmpCmd) = qsubSystem($logDir."RiboLCA.sh",$cmd2,$numCore2,$mem,"_RA$JNUM",$jobName,"",1,[],$QSBoptHR);
			$QSBoptHR->{tmpSpace} =$tmpSHDD;
			$QSBoptHR->{useLongQueue} = 0;
		}
	}
	return $jobName;
}




sub prepDiamondDB($ $ $ $){#takes care of copying the respective DB over to scratch
	my ($curDB,$CLrefDBD,$ncore, $searchMode) = @_;
	#searchMode 1: diamond, 2: mmseqs2
	my $diaBin = getProgPaths("diamond");
	my $mmseqs2Bin = getProgPaths("mmseqs2");
	my $dbSuffix = ".db.dmnd";

	my ($DBpath ,$refDB ,$shrtDB) = getSpecificDBpaths($curDB,0);
	system "mkdir -p $CLrefDBD" unless (-d $CLrefDBD);
	my $clnCmd = "";
	if ($MFopt{globalDiamondDependence}->{$curDB} eq "" ){
		my $DBcmd = ""; 
		# check for pr.db.dmnd
		my $ncoreDB=1;
		my $epoTS = 0; my $refDBnew=0;
		my $epoTSdia = 99999999999999;
		
		if ($searchMode == 1 ){
			$epoTSdia = ( stat "$diaBin" )[9] if (-f $diaBin);
		} else {
			$epoTSdia = ( stat "$mmseqs2Bin" )[9];
			$dbSuffix = ".db.mms2";
		}
		$epoTS = ( stat "$DBpath$refDB${dbSuffix}" )[9] if (-e "$DBpath$refDB${dbSuffix}");
		#$timestamp = POSIX::strftime( "%d%m%y", localtime( $epoTS));
		#print "time: $epoTSdia $epoTS \n";
		if (!-e "$DBpath$refDB${dbSuffix}"){
		# age check completely deactivated..
		#if ($epoTS < $epoTSdia){#checks age of binary vs DB creation.. 
			$ncoreDB = $ncore;$refDBnew=1; 
			system "rm -f $DBpath$refDB${dbSuffix}*";
			if ($searchMode == 1){
				$DBcmd .= "$diaBin makedb --in $DBpath$refDB -d $DBpath$refDB${dbSuffix} -p $ncoreDB\n";
			} else {
				$DBcmd .= "$mmseqs2Bin createdb $DBpath$refDB $DBpath$refDB${dbSuffix} --compressed 1\n";
			}
		}
		unless (-e "$DBpath$refDB.length"){
			my $genelengthScript = getProgPaths("genelength_scr");#
			$DBcmd .= "$genelengthScript $DBpath$refDB $DBpath$refDB.length\n";
		}
		#$clnCmd .= "rm -rf $CLrefDBD;" if (length($CLrefDBD)>6);
		#idea here is to copy to central hdd (like /scratch)
		system "rm -f $CLrefDBD/$refDB${dbSuffix}*"  if (($refDBnew || $MFopt{rewriteDiamond} )&& -d $CLrefDBD);
		#print " !-e $CLrefDBD/$refDB.db.dmnd && !-e $CLrefDBD/$refDB.length\n ";
		if ( -e "$CLrefDBD/$refDB${dbSuffix}" && -e "$CLrefDBD/$refDB.length" 
			#&& ($curDB eq "NOG" && !-e "$CLrefDBD/NOG.members.tsv") &&
			#(($curDB ne "KGB" && $curDB ne "KGM" && $curDB ne "KGE")|| -s "$CLrefDBD/genes_ko.list")
			){
			#has to be noted that this doesn't need to happen again
			$MFopt{globalDiamondDependence}->{$curDB}="$shrtDB-1";
		} else {
			$DBcmd .= "mkdir -p $CLrefDBD\n";
			$DBcmd .= "cp $DBpath$refDB${dbSuffix}*  $DBpath$refDB.length  $CLrefDBD\n";
		}
		
		#specialized file copy
		if ($curDB eq "NOG" && !-s "$CLrefDBD/NOG.members.tsv" && !-s "$CLrefDBD/NOG.annotations.tsv"){
			system "rm -f $CLrefDBD/NOG* $CLrefDBD/all_species_data.txt";
			$DBcmd .= "cp $DBpath/all_species_data.txt $DBpath/NOG.members.tsv $DBpath/NOG.annotations.tsv $CLrefDBD\n";
		}
		if ($curDB eq "CZy" && !-s "$CLrefDBD/MohCzy.tax"){
			system "rm -f $CLrefDBD/MohCzy.tax $CLrefDBD/cazy_substrate_info.txt";
			$DBcmd .= "cp $DBpath/MohCzy.tax $DBpath/cazy_substrate_info.txt $CLrefDBD\n";
		}
		if (($curDB eq "KGB" || $curDB eq "KGM" || $curDB eq "KGE") && !-s "$CLrefDBD/genes_ko.list"){
			system "rm -f $CLrefDBD/genes_ko.list $CLrefDBD/kegg.tax.list";
			$DBcmd .= "cp $DBpath/genes_ko.list $DBpath/kegg.tax.list $CLrefDBD\n";
		}
		if ($curDB eq "ABRc" && !-s "$CLrefDBD/card.parsed.f11.tab.map"){ 
			system "rm -f $CLrefDBD/card*";
			$DBcmd .= "cp $DBpath/card*.txt $DBpath/card*.map $CLrefDBD\n";
		}
		if ($curDB eq "PTV" && !-s "$CLrefDBD/PATRIC_VF2.tab"){ 
			system "rm -f $CLrefDBD/PATRIC_VF2.tab";
			$DBcmd .= "cp $DBpath/PATRIC_VF2.tab $CLrefDBD\n";
		}
		if ($curDB eq "VDB" && !-s "$CLrefDBD/VF.tab"){ 
			system "rm -f $CLrefDBD/VF.tab";
			$DBcmd .= "cp $DBpath/VF.tab $CLrefDBD\n";
		}
		if ($curDB eq "PAB" && !-s "$CLrefDBD/all_species_data.txt"){
			#copy NOG taxonomy
			my ($DBpathN) = getSpecificDBpaths("NOG",0);
			$DBcmd .= "cp $DBpathN/all_species_data.txt $CLrefDBD\n" unless (-e "$CLrefDBD/all_species_data.txt");
		}
		if ($curDB eq "TCDB" && !-s "$CLrefDBD/hir.txt"){ 
			$DBcmd .= "cp $DBpath/TCDBhir.txt  $CLrefDBD\n";
		}

		my $jN = "_DIDB$shrtDB$JNUM"; my $tmpCmd;
#		die "$DBcmd";
		if ($DBcmd ne ""){
			my @preConstr = @{$QSBoptHR->{constraint}};
			#push(@{$QSBoptHR->{constraint}}, "intel");
			my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0; 
			($jN, $tmpCmd) = qsubSystem($logDir."DiamondDBprep$shrtDB.sh",$DBcmd,$ncoreDB,int(100)."G",$jN,"","",1,$QSBoptHR->{General_Hosts},$QSBoptHR);
			$QSBoptHR->{tmpSpace} =$tmpSHDD;
			$MFopt{globalDiamondDependence}->{$curDB} = $jN;
			@{$QSBoptHR->{constraint}} = @preConstr;
			#die "$MFopt{globalDiamondDependence}->{$curDB}";
		}
		#die $DBcmd."\n";
	}
	return ($refDB,$shrtDB,$clnCmd);
}


sub getRdLibraries {
	my ($libraries, $mergedLibrary) = @_;
	my $pairs = libraryPairs($libraries);
	my @reads1 = map { $_->{files}{r1} } @{$pairs};
	my @reads2 = map { $_->{files}{r2} } @{$pairs};
	my @single = @{libraryFiles($libraries, 'single')};
	my $mergeMode = ref($mergedLibrary) eq 'HASH' && ($mergedLibrary->{files}{single} || '') ne '';
	my %result;
	$result{0} = \@single if @single;
	my $usedRead2 = $mergeMode ? [$mergedLibrary->{files}{r2}] : \@reads2;
	my $usedRead1 = $mergeMode ? [$mergedLibrary->{files}{r1}] : \@reads1;
	$result{1} = $usedRead2 if grep { defined($_) && $_ ne '' } @{$usedRead2};
	$result{2} = $usedRead1 if grep { defined($_) && $_ ne '' } @{$usedRead1};
	$result{3} = [$mergedLibrary->{files}{single}] if ($mergeMode);
	return %result;
}

sub isLastSampleInAssembly{
	my ($assD,$curD) = @_;
	return 0 if (!-d $assD);
	return 0 if (!-e "$assD/smpls_used.txt");
	open I,"<$assD/smpls_used.txt" or die $!;
	my $lastEntr="";
	while (<I>){
		chomp $_; 
		$lastEntr = $_ if (length($_)>2);
	}
	close I;
	#print "$lastEntr  $curD \n";
	return 1 if ($lastEntr eq $curD);
	return 0;
}


sub workflowStateOptions {
	return {
		assembly_mode => $MFopt{DoAssembly},
		map_to_assembly => $MFopt{map2Assembly},
		map_support_to_assembly => $MFopt{mapSupport2Assembly},
		run_tmp_dir => $MFglobal{runTmpDirGlobal},
	};
}

sub publishKrakenResults{
	my ($curOutDir,$SmplName) = @_;
	my $KrakenOD = $curOutDir."Tax/kraken/$MFopt{globalKraTaxkDB}/";
	return unless (-d $KrakenOD && -e "$KrakenOD/krakDone.sto");
	my $dir_KrakFind = $baseOut."pseudoGC/Phylo/KrakenTax/$MFopt{globalKraTaxkDB}/";
	opendir(my $dh, $KrakenOD) or die "Can't read completed Kraken directory $KrakenOD\n";
	my @krkF = grep { /^krak\.(.+)\.cnt\.tax$/ && -s "$KrakenOD/$_" } readdir($dh);
	closedir($dh);
	foreach my $kf (@krkF){
		$kf =~ /^krak\.(.+)\.cnt\.tax$/;
		my $thr = $1;
		my $dest = "$dir_KrakFind/$thr";
		system('mkdir', '-p', $dest) == 0 or die "Can't create Kraken publication directory $dest\n";
		system('cp', "$KrakenOD/$kf", "$dest/$SmplName.$thr.krak.txt") == 0
			or die "Can't publish Kraken result $KrakenOD/$kf\n";
	}
}


sub runStateInspection {
	my $report = inspect_workflow_state(
		map => \%map,
		groups => \%AsGrps,
		options => workflowStateOptions(),
	);
	if ($MFconfig{stateReport} ne '') {
		my $json = encode_state_report($report);
		open my $fh, '>', $MFconfig{stateReport}
			or die "Cannot write state report $MFconfig{stateReport}: $!\n";
		print {$fh} $json;
		close $fh;
		print "Wrote read-only state report to $MFconfig{stateReport}\n";
	}

	if ($MFconfig{planState}) {
		my $plan = build_workflow_plan($report);
		my $json = encode_workflow_plan($plan);
		if ($MFconfig{planReport} ne '') {
			open my $fh, '>', $MFconfig{planReport}
				or die "Cannot write workflow plan $MFconfig{planReport}: $!\n";
			print {$fh} $json;
			close $fh;
			print "Wrote read-only repair/submission plan to $MFconfig{planReport}\n";
		} else {
			print STDOUT $json;
		}
	} elsif ($MFconfig{stateReport} eq '') {
		print STDOUT encode_state_report($report);
	}
}


sub runAutomaticWorkflowPreflight {
	my ($iteration) = @_;
	my $applyRepairs = $runOptions{submit} && $MFconfig{autoRepairState};
	my $result = run_workflow_preflight(
		map => \%map,
		groups => \%AsGrps,
		options => workflowStateOptions(),
		apply_repairs => $applyRepairs,
		allow_group_rewrite => $MFconfig{OKtoRWassGrps},
		iteration => $iteration,
	);

	my $auditDir = $baseOut ne '' ? "$baseOut/LOGandSUB/workflow" : '';
	if ($auditDir ne '') {
		make_path($auditDir) unless (-d $auditDir);
		my $suffix = sprintf('%03d', $iteration);
		my $statePath = "$auditDir/state.iteration-$suffix.json";
		my $planPath = "$auditDir/plan.iteration-$suffix.json";
		open my $stateFH, '>', $statePath
			or die "Cannot write automatic state report $statePath: $!\n";
		print {$stateFH} encode_state_report($result->{state});
		close $stateFH;
		open my $planFH, '>', $planPath
			or die "Cannot write automatic workflow plan $planPath: $!\n";
		print {$planFH} encode_workflow_plan($result->{plan});
		close $planFH;
	}

	my $repairSummary = $result->{repairs};
	my $planSummary = $result->{plan}{summary};
	my $repairWord = $applyRepairs ? 'removed' : 'would remove';
	my $repairCount = $applyRepairs
		? $repairSummary->{removed_targets} : $repairSummary->{would_remove_targets};
	print "Workflow preflight iteration $iteration: $planSummary->{submissions} pending submissions; "
		."$repairWord $repairCount safe partial targets; "
		."$repairSummary->{blocked_repairs} protected repairs require explicit authorization.\n";
	return $result;
}


#preparation of secondary mapping, including wildcard resolution, gene calling, index building
sub prepareMap{
	

	my ($hr,$hr2) = readMapS($MFconfig{mapFile},0,\%map,\%AsGrps,$MFconfig{oldStylFolders});
	checkAssmblGrp($hr);

	%AsGrps = %{$hr2}; %map = %{$hr};
	if ($MFopt{DoMetaBat2}){ #do any binning?
		my ($hrD) = getDirsPerAssmblGrp(\%map,\%AsGrps);
		%DOs = %{$hrD};
	}



	#dirs from config file--------------------------
	#can be overwritten by $map{opt}{GlbTmpD} $map{opt}{NodeTmpD}
	$runOptions{sharedTmpDir} = getProgPaths("globalTmpDir",0) unless ($runOptions{sharedTmpDir} ne "");
	$runOptions{nodeTmpDir} = getProgPaths("nodeTmpDir",0) unless ($runOptions{nodeTmpDir} ne "");
	#useless, because baseout can change..
	$baseOut = $map{opt}{outDir} if (exists($map{opt}{outDir} ) && $map{opt}{outDir}  ne "" && $map{opt}{outDir} !~ m/,/);
	$runOptions{baseID} = $map{opt}{baseID} if (exists($map{opt}{baseID} ) && $map{opt}{baseID}  ne "");
	#die "baseout dir has \",\": $baseOut\nDid you supply multiple maps? (not supported)\n";
	#overwrite tmp dirs??
	if ($map{opt}{GlbTmpD} ne ""){print "Taking Global temp dir from map: $map{opt}{GlbTmpD} \n";$runOptions{sharedTmpDir} = $map{opt}{GlbTmpD} ;}
	if ($map{opt}{NodeTmpD} ne ""){print "Taking Node temp dir from map: $map{opt}{NodeTmpD} \n";$runOptions{nodeTmpDir} = $map{opt}{NodeTmpD} ;}




	$MFglobal{runTmpDirGlobal} = "$runOptions{sharedTmpDir}/$runOptions{baseID}/";
	# Inspection is a read-only planning path: do not create scratch directories,
	# prepare databases, or enqueue any prerequisite jobs.
	return if ($MFconfig{inspectState});


	
	unless ($runOptions{operationMode} eq "scaffold" || $runOptions{operationMode} eq "map2tar" || $runOptions{operationMode} eq "map2DB" || $runOptions{operationMode} eq "map2GC"){
		#not asked for 2nd map? ok, deactivate all related parameters
		$MFopt{mapModeTogether} = 0; $MFopt{DoMapModeDecoy} = 0; $MFopt{mapModeActive} = 0;
		return;
	}
	if ($runOptions{operationMode} eq "scaffold"){
		die"update scaffold\n";
		$scaffTarExternal = $ARGV[1];
		if (!-f $scaffTarExternal){
			die "Could not find scaffold file:\n$scaffTarExternal\n";
		}
		$scaffTarExternalName = $ARGV[2];
		if (@ARGV>3){
			$scaffTarExtLibTar = $ARGV[3];
		}
		return;
			

	}

#in this case primary focus is on mapping and not on assemblies
	if ($runOptions{operationMode} eq "map2DB" || $runOptions{operationMode} eq "map2GC"){
		$MFopt{mapModeCovDo}=0;$MFopt{DoMapModeDecoy}=0;$MFopt{mapModeTogether}=0;$map2ndMpde=2;
	}
	if ($MFopt{map2Assembly} || $MFopt{DoAssembly}){
		print "Mapping mode: reference ";
		print " Decoy" if ($MFopt{DoMapModeDecoy});
		print " Competitive" if ($MFopt{mapModeTogether}== 1);
		print " combined map, seperate reporting" if ($MFopt{mapModeTogether}== 2);
		print " combined map, combined reporting" if ($MFopt{mapModeTogether}== -1);
		print "\nDeactivating assembly and dependent modules.\n";
		$MFopt{map2Assembly}=0; $MFopt{DoAssembly} =0;
	}
	my $DBsubmCnt=0; my $GENEsubmCnt=0;
	my @refDB1 = split(/,/,$MFopt{refDBall});
	#die "@refDB1\n";
	#$MFopt{mapModeTogether} = 0 if (@refDB1 == 1);#could still be a wildcard, wrong assumption
	my @bwt2Name1;
	if (defined $MFopt{bwt2NameAll} && $MFopt{bwt2NameAll} ne "" ){
		@bwt2Name1 = split(/,/,$MFopt{bwt2NameAll}) ;
	} elsif ($runOptions{operationMode} eq "map2DB"  ){
		@bwt2Name1 = ("refDB");
	} elsif ($runOptions{operationMode} eq "map2GC"){
		@bwt2Name1 = ("GC");
		$map2ndMpde=3;
	} else {
		@bwt2Name1 = ("auto") x scalar(@refDB1);
	}
	#die "@bwt2Name1\n";
	my @refDB;my @bwt2Name ;
	my %FNrefDB2ndmap;
	
	for (my $i=0;$i<@refDB1;$i++){
		my @sfiles;
		if ($runOptions{operationMode} eq "map2GC"){
			die "can only have one ref to GC: @refDB1\n" unless (@refDB1 == 1);
			@sfiles = ($refDB1[0]."/compl.incompl.95.fna");
			die "Could not find reference in GC dir. Expected: $sfiles[0]\n" unless (-e $sfiles[0]);
		} else {
			@sfiles = glob($refDB1[$i]);
		}
		#die "@sfiles\n$refDB1[$i]\n".@sfiles."\n";
		my $iniBwtNm = $bwt2Name1[$i];
		#die "@sfiles\n$refDB1[$i]\n";
		if (@sfiles>1){
			for (my $j=0;$j<@sfiles;$j++){
				push(@refDB,$sfiles[$j]);
				if ($iniBwtNm eq "auto"){
					#$sfiles[$j] =~ m/ssemblyfind_list_(.*)\.txt\/(.*)\.contigs_/; my $nmnew = $1.$2; $nmnew =~ s/#/_/g;
					$sfiles[$j] =~ m/\/([^\/]+)\.f.*a$/;
					my $nmnew = $1;
					#die "$nmnew\n";
					push(@bwt2Name,$nmnew);

				} else {
					push(@bwt2Name,$bwt2Name1[$i].$j);
				}
			}
		} elsif (@sfiles == 1) {
			push(@refDB,$sfiles[0]);
			if ($iniBwtNm eq "auto"){
				$sfiles[0] =~ m/.*\/([^\/]+)$/;
				$iniBwtNm = $1;
				$iniBwtNm =~ s/\.[^\.]+$//;
			}
			push(@bwt2Name,$iniBwtNm);
		} else {
			die "Could not find file for entry $refDB1[$i]\n";
		}
	}
	#die "@refDB\n@bwt2Name\n";
	my $shrtMapNm = "";	$shrtMapNm = "Comb_" if (@bwt2Name > 1);
	for ( my $i=0;$i< @bwt2Name; $i++){
		my $substrl = length($bwt2Name[$i]);	if (@bwt2Name > 3){$substrl = 8;}	
		if (@bwt2Name > 5){$substrl = 6;}	if (@bwt2Name > 7){$substrl = 3;}
		if ($i==0){
			$shrtMapNm .= substr($bwt2Name[$i],0,$substrl);
		} else {
			$shrtMapNm .= ".". substr($bwt2Name[$i],0,$substrl);
		}
		if ($i>5){$shrtMapNm .= ".X".(@bwt2Name - $i)."X"; last;}
	}
	$map2ndTogRefDB{DB} = "$baseOut/GlbMap/$shrtMapNm/$shrtMapNm.fa";
	#die"$map2ndTogRefDB{DB}\n";
	#decoy mapping setup (only required in map2tar
	$make2ndMapDecoy{Lib} = "";
	#die "decoy mapping not ready for multi fastas\n" if (@refDB > 1);
	if ($MFopt{DoMapModeDecoy} || $MFopt{mapModeTogether}){
		if ($MFopt{mapModeTogether} && $MFopt{MapRewrite2nd} ){
			#system("rm -r -f $map2ndTogRefDB{DB}*");#;mkdir -p $bwt2outDl
		}
		for (my $i=0;$i<@refDB; $i++){ #take care of ref DB decoy prep
			#last if (-e $map2ndTogRefDB{DB});
			#print "$refDB[$i]\n";
			my $aref;
			if ($MFopt{mapModeTogether}){#read into mem & combine
				my $hr = readFasta($refDB[$i],1); 
				$hr = prefixFAhd($hr,$bwt2Name[$i]);#rename header to fasta files...
				my %FN = %{$hr};
				%FNrefDB2ndmap = (%FNrefDB2ndmap , %FN);
				$aref = [keys %FN];
			} else {
				$aref = readFastHD($refDB[$i]);
			}
			if ($MFopt{mapModeTogether}>0 || $MFopt{DoMapModeDecoy}){
				push(@{$make2ndMapDecoy{regions}}, join(" ",@{$aref}) );
				push(@{$make2ndMapDecoy{region_lcs}},  lcp(@{$aref}) );#prefix_find($aref)  );
			}
		}
	}
	#		die "@{$make2ndMapDecoy{region_lcs}}\n";
	
	#build of combined refDB, if competitive mapping..
	my $bwtDBcore = 10; 
	
	if ($MFopt{mapModeTogether}){#the fasta's were already combined into %FNrefDB2ndmap, just built idx now..
		system "mkdir -p $baseOut/GlbMap/LOGandSUB/" unless (-d "$baseOut/GlbMap/LOGandSUB/");
		system "mkdir -p $baseOut/GlbMap/$shrtMapNm" unless (-d "$baseOut/GlbMap/$shrtMapNm");
		#system "mkdir -p $map2ndTogRefDB{DB}" unless (-d $map2ndTogRefDB{DB});
		print "$map2ndTogRefDB{DB}\n";
		#die;
		my @combinedReferenceStat = stat($map2ndTogRefDB{DB});
		writeFasta(\%FNrefDB2ndmap,"$map2ndTogRefDB{DB}")
			unless (@combinedReferenceStat && $combinedReferenceStat[7] > 0);
		if ($MFopt{mapModeTogether}==-1){
			@refDB = ($map2ndTogRefDB{DB});
			@bwt2Name = ($shrtMapNm);
			$MFopt{mapModeTogether} = 0;#deactivate, as from now will be treated as singular ref
		} else {
			my ($cmd,$DBbtRef, $chkFile) = buildMapperIdx($map2ndTogRefDB{DB},$bwtDBcore,$MFopt{largeMapperDB},$MFopt{MapperProg}) ;
			if (!-e $chkFile ){
				my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0; 
				($bwt2ndMapDep,$cmd) = qsubSystem("$baseOut/GlbMap/LOGandSUB/builBwtIdx_comp.sh",$cmd,$bwtDBcore,(int(25)+1) ."G","BWI_compe","","",1,[],$QSBoptHR) ;
				$QSBoptHR->{tmpSpace} =$tmpSHDD;
			}
		}
	}
	
	print "\n=======================\nmap to $MFopt{refDBall}\n with map mode $MFopt{mapModeTogether}\n=======================\n\n";
	$MFopt{mapModeActive} =1; my $cmdBIG = "";my $cmdBIGgene = "";
	#die "@refDB\n";
	
	#build index for each fasta, and predict genes on these
	die "MATAFILER map2nd check failed: @bwt2Name != @refDB" if (@bwt2Name != @refDB);
	#die "X:: @refDB  @bwt2Name\n";
	for (my $i=0;$i<@refDB; $i++){
		#$refDB[$i] =~ m/(.*\/)[^\/]+/;
		#my $refDir = $1;
		my $bwt2outDl = "$baseOut/GlbMap/$bwt2Name[$i]/";
		if ($map2ndMpde == 3){#outdir should be in the GC dir
			$bwt2outDl = "$refDB1[0]/unmappedMap/"
		}
		push @bwt2ndMapNmds , $bwt2Name[$i];
		push(@bwt2outD,$bwt2outDl);
		system "mkdir -p $bwt2outDl" unless (-d $bwt2outDl);
		if ($MFopt{mapModeTogether} >= 0 && $MFopt{MapRewrite2nd} ){
			print "Rebuilding requested mapping DBs..\n" if ($i==0);
		}
		
		$refDB[$i] =~ m/.*\/([^\/]+)$/;
		system "cp $refDB[$i] $bwt2outDl" if ($MFopt{mapModeCovDo}  && !-e "$bwt2outDl/$1");
		#print "\n$refDB[$i]\n";
		#die "$bwt2outDl/$1\n";
		#system "mkdir -p $bwt2outDl/LOGandSUB" unless (-d "$bwt2outDl/LOGandSUB");
		
		
		
		my ($cmd,$DBbtRef,$chkFile) = buildMapperIdx($refDB[$i],$bwtDBcore,$MFopt{largeMapperDB},$MFopt{MapperProg}) ;
		# The FASTA and its adjacent mapper index are durable inputs.  Keeping the
		# canonical reference here also ensures CRAM and FASTA-based mappers do not
		# receive a path into global scratch.
		$DBbtRef = $refDB[$i];
		my $idxNFini = !mapperDBbuilt($refDB[$i],$MFopt{MapperProg});
		#print $cmd."\n";
		if (!$MFopt{mapModeTogether} && $idxNFini && $cmd ne ""){ #not required for these map modi
			#system $cmd 
			#	my ($bwt2ndMapDep2,$cmd2) = qsubSystem($bwt2outDl."/LOGandSUB/builBwtIdx$i.sh",$cmdBIG,$bwtDBcore,(int(20/$bwtDBcore)+1) ."G","BWI".$i,"","",1,[],$QSBoptHR) ;$bwt2ndMapDep .= ";$bwt2ndMapDep2";
			$cmdBIG .= "\n\n#====== $i =======\n".$cmd;
			$DBsubmCnt++;
		}
		
		#$DBbtRefX = $DBbtRef;
		#die "$DBbtRef\n";
		push(@DBbtRefX,$DBbtRef);
		if($MFopt{mapModeCovDo} && $map2ndMpde != 3){ #get the coverage per gene etc; for this I need a gene prediction
												#but not for GC mapping (these are genes already)
			my $gDir = $bwt2outDl."";
			my $nativeGFF = $refDB[$i];$nativeGFF =~ s/\.[^\.]+$/\.gff/;
			my $gffF = "genes.$bwt2Name[$i].gff";
			#die "$nativeGFF\n";
			if (-e "$gDir/$gffF"){
				;
			}elsif (-e $nativeGFF){
				system "cp $nativeGFF $gDir/$gffF";
			} else {
				system "mkdir -p $gDir";
				$logDir = $gDir;
				my $dEGP = $MFopt{DoEukGenePred}; $MFopt{DoEukGenePred} = 0;
				my $tmpDep1 = genePredictions($refDB[$i],$gDir,"",$gDir,"iGP$i","",0);
				$MFopt{DoEukGenePred} = $dEGP;
				$cmdBIGgene .= "#====== $i =======\n".$tmpDep1."cp $gDir/genes.gff $nativeGFF; mv $gDir/genes.gff $gDir/$gffF\n\n\n";
				$GENEsubmCnt++;
				#my ($tmpDep,$tmpCmd) = qsubSystem( $bwt2outDl."/LOGandSUB/cpGenes.sh",  "cp $gDir/genes.gff $nativeGFF; mv $gDir/genes.gff $gDir/$gffF\n",
				#1,"1G","genecop".$i,$tmpDep1,"",1,[],$QSBoptHR);
				#$bwt2ndMapDep .= ";".$tmpDep;
			}
			push(@DBbtRefGFF,$gDir."/$gffF");#"genePred/genes.gff"
		}
		$logDir="";
	} #end for loop building bwtIdx
	#now submit all together as single call..
	my $bwt2outDl = "$baseOut/GlbMap/LOGandSUB/"; system "mkdir -p $bwt2outDl" unless (-d $bwt2outDl);
	#submit mapping index build
	#die "$cmdBIG\n\n";
	if ($DBsubmCnt>0){
		my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0; 
		my $mapperMemDB= 20; $mapperMemDB = 40 if ($MFopt{largeMapperDB});
		my ($bwt2ndMapDep2,$cmd2) = qsubSystem($bwt2outDl."/builBwtIdxBIG.sh",$cmdBIG,$bwtDBcore,(int($mapperMemDB)+1) ."G","BWIbig","","",1,[],$QSBoptHR) ;
		$QSBoptHR->{tmpSpace} =$tmpSHDD;
		$bwt2ndMapDep .= ";$bwt2ndMapDep2";
	}
	#submit gene predictions
	if ($GENEsubmCnt>0){
		my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0; 
		my ($tmpDep,$tmpCmd) = qsubSystem( $bwt2outDl."/GenesPredBIG.sh",  $cmdBIGgene,
			1,"1G","genePred","","",1,[],$QSBoptHR) ;
		$QSBoptHR->{tmpSpace} =$tmpSHDD;
		$bwt2ndMapDep .= ";".$tmpDep;
	}
	
	
	if ($map2ndMpde == 3 && !$MFopt{useUnmapped}){
		die "Did you mean to activate \"-mapUnmapped\"?\n";
	}

	#die "@bwt2outD\n";
}


sub runOrthoPlacement(){
	my ($outD,$tmpP,$jdep) = @_;
	die "runOrthoPlacement no longer active\n";
	my $cleanSeqSetHR = sampleReadSet($curSmpl, "clean");
	my $libraries = readLibrariesByScope($cleanSeqSetHR, 'primary', 1, $curSmpl);
	my $mergedLibrary = $cleanSeqSetHR->{merged_library};
	
	my $sdmBin = getProgPaths("sdm");#"/g/bork3/home/hildebra/dev/C++/sdm/./sdm";
	my $fna2faaBin = getProgPaths("fna2faa");
	system "mkdir -p $outD $tmpP" unless (-d $outD && -d $tmpP);
	my %RdLibs = getRdLibraries($libraries,$mergedLibrary);
	my $numCore=22;
	my $scrP = "";my $cmd="";
	$cmd .= "mkdir -p $tmpP\n";
	#my $tmpF = ${$RdLibs{$kk}}[0];
	system "rm -f $tmpP/6fr.fna";
	foreach my $kk (sort {$a <=> $b} keys %RdLibs){
		my @rds = @{$RdLibs{$kk}};
		$scrP = $rds[0] if ($scrP eq "");
		foreach my $fnaI (@rds){
			$cmd .= "$sdmBin -i $fnaI -o_fna - | $fna2faaBin -q - >> $tmpP/6fr.fna\n";
		}
	}
	$scrP =~ s/\/[^\/]+$//;
	#die $scrP."\n";
	my $hmmscr = "python /g/bork3/home/hildebra/dev/Perl/SoilHelpers/extract_domains.py";
	my $hmmD = "/g/bork3/home/hildebra/DB/HMMs/FungiAB/";
	my @hmmsDB=("Condensation.hmm","dmat.hmm","AMP-binding.hmm","PKS_KS.hmm","Terpene_synth_C.hmm");
	my @hmmNms = ("Condensation","dmat","AMP","PKS","Terpene");
	my $onceExtr=0;
	my $redo=0;
	for (my $i=0; $i<@hmmsDB; $i++){	
		my $hmmM  = "$hmmD/$hmmsDB[$i]";
		my $hmmName = $hmmsDB[$i]; $hmmName=~s/\.hmm//;
		if (!-e "$outD/$hmmName.fna" || $redo){
			$cmd .= "$hmmscr $tmpP/6fr.fna $hmmM $outD/$hmmName.fna $numCore\n" ;
			$onceExtr=1;
		}
	}
	if($onceExtr ==0 && !$redo){
		$cmd = "";
	}
	$cmd .= "mkdir -p $tmpP\n";
		
	#second new part .. alignment to existing tree
	$onceExtr=0;
	for (my $i=0; $i<@hmmsDB; $i++){	
		my $hmmName = $hmmsDB[$i]; $hmmName=~s/\.hmm//;
		my $hmmN2 = $hmmNms[$i];
		my $tarFile = "$outD/$hmmName.fna";
		my @targetStat = stat($tarFile);
		if (@targetStat && $targetStat[7] == 0){next;}
		#$tarFile = fixHDs4Phylo($tarFile);
		#fix fasta headers
#		$cmd .= "cat $tarFile | perl -p -e 's/[:|#+]/_/g' > $tarFile\n";

		my $clustaloBin = getProgPaths("clustalo");
		my $raxMLbin = getProgPaths("raxml");
		my $treeDBdir = "/g/bork3/home/hildebra/data/SoilABdoms/Jaime/JaimeT2/";
#		my $refTREE = "$treeDBdir/$hmmN2/final_alg/phylo/FASTTREE_allsites.nwk";
		my $refTREE = "$treeDBdir/$hmmN2/final_alg/phylo/IQtree_fast_allsites.treefile";
		my $refMSA = "$treeDBdir/$hmmN2/final_alg/outMSA.faa";
		
		if (!-e "$outD/RAxML_classification.${hmmName}_place"){
			$onceExtr=1;
			$cmd .= "$clustaloBin --p1 $refMSA -i $tarFile --threads $numCore > $tmpP/$hmmName.msa.faa\n";
			$cmd .= "sed -i 's/[:|#+]/_/g' $tmpP/$hmmName.msa.faa\n";
			$cmd .= "rm -f $tmpP/RAxML*\n";
			$cmd .= "$raxMLbin -f v -s $tmpP/$hmmName.msa.faa  -w $tmpP -t $refTREE -m PROTGAMMALG -n ${hmmName}_place -T$numCore\n";
			#$cmd .= "mv $tmpP/RAxML_fastTreeSH_Support.${hmmName}_place $outD/${hmmName}_place.nwk\n";
			$cmd .= "mv $tmpP/*${hmmName}_place $outD/\n";
		}
	}
	
		#die $cmd;
	
	
	$cmd .= "rm -r $tmpP\n";
	#die $cmd;
	my $jobName = "_ORTH$JNUM"; 
	if ($onceExtr){
		my ($jobName2, $tmpCmd) = qsubSystem($logDir."PABpred.sh",$cmd,$numCore,"4G",$jobName,$jdep,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR);
		return $jobName2;
	} else { return "";}
}




sub runDiamond(){
	my ($outD,$CLrefDBD,$tmpP,$jdep,$curDB_o) = @_;
	my $cleanSeqSetHR = sampleReadSet($curSmpl, "clean");
	my $libraries = readLibrariesByScope($cleanSeqSetHR, 'primary', 1, $curSmpl);
	my $mergedLibrary = $cleanSeqSetHR->{merged_library};
	
	my $secCogBin = getProgPaths("secCogBin_scr");
	my $diaBin = getProgPaths("diamond");
	my $mmseqs2Bin = getProgPaths("mmseqs2");
	my $searchMode = 1;
	
	my %RdLibs = getRdLibraries($libraries,$mergedLibrary);
	
	die "No reads found for sample $outD in runDiamond sub\n" if (scalar(keys(%RdLibs)) == 0);
	
	my $clnCmd="";my $jobN2="";
	#die;
	system "mkdir -p $outD" unless (-d $outD && -d $tmpP);
	#my $diaOfmt = "tab";#old
	#$Query,$Subject,$id,$AlLen,$mistmatches,$gapOpe,$qstart,$qend,$sstart,$send,$eval,$bitSc
	my $diaOfmt = "-f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore "; #diamond
	
	my $mmsOfmt = "--format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"; #mmseqs
	my $ncore = $MFopt{diaCores}; my $sensBlast = "";
	$sensBlast = " --sensitive " if ($MFopt{diaRunSensitive});
	$sensBlast .= " -F $MFopt{diaFrameshift} " if ($MFopt{diaFrameshift});
	foreach my $curDB (split /,/,$curDB_o){
		$progStats{$curDB}{DiaDBSearchCompl} = 0 unless (defined($progStats{$curDB}{DiaDBSearchCompl}));
		$progStats{$curDB}{DiaDBSearchIncomplete} = 0 unless (defined($progStats{$curDB}{DiaDBSearchIncomplete}));
		#print "$curDB";
		my ($refDB,$shrtDB,$clnCmd) = prepDiamondDB($curDB,$CLrefDBD,$ncore,$searchMode);
		my $doInterpret = 1;
		#print "$refDB ,$shrtDB,$clnCmd\n";
		#$doInterpret = 0 if ($shrtDB eq "ABR");
		my $getQSeq = 0;
		$getQSeq = 0 if ($curDB eq "PAB");
		my $diaOfmt2 = $diaOfmt;
		$diaOfmt2 .= " qseq" if ($getQSeq); #in case, I want to get the query sequence (matching)
		#run actual diamond
		my @collect = (); my @collectSingl=();
		my $cmd ="mkdir -p $tmpP\n";
		foreach my $kk (sort {$a <=> $b} keys %RdLibs){
			my @rds = @{$RdLibs{$kk}};
			for (my $ii=0;$ii<@rds;$ii++){
				my $query = $rds[$ii];	
				my $outF = "$tmpP/DiaAssignment.sub.$shrtDB.$kk.$ii";
				#my $tmpcnt = `grep -c '^>' $query`; chomp $tmpcnt
				#--comp-based-stats 0
				if ($searchMode==1){
				$cmd .= "$diaBin blastx $diaOfmt2 --masking 0 --comp-based-stats 0 --compress 1 --quiet -t $tmpP --min-orf 25 -d $CLrefDBD$refDB.db -q $query -k 5 -e 1e-4 -o $outF $sensBlast -p $ncore\n"; #
				} elsif ($searchMode==2){
				$cmd .= "$mmseqs2Bin easy-search $query $CLrefDBD$refDB.db.mms2 $outF.gz $tmpP --threads $ncore --max-accept 500 --compressed 1 -s 4 $mmsOfmt \n";
				}
				
				#$cmd .= "$diaBin view -a $outF.tmp -o $outF -f tab\nrm $outF.tmp.daa\n";
				if ($kk==0 || $kk == 3){#single or ext fragments, doesn't need to be sorted
					push(@collectSingl,$outF.".gz");
				} else {
					push(@collect,$outF.".gz");
				}
			}
		}
		
		#die "$cmd\n";
		
		my $out = $outD."dia.$shrtDB.blast";
		my $outgz = "$out.srt.gz";
		#die "$outgz\n";
		#unzip, sort, zip
		if (@collect) {
			$cmd .= "zcat ".join( " ",@collect) ." | sort -t\$'\\t' -k1 -T $tmpP | $pigzBin --stdout -p $ncore > $outgz \n";
			$cmd .= "rm -f ". join( " ",@collect) . "\n";
		} else {
			$cmd .= "$pigzBin -c </dev/null > $outgz\n";
		}
		if (@collectSingl >= 1){
			#append on gzip, can be done with gzip
			$cmd .= "cat ".join( " ",@collectSingl) ." >> $outgz\nrm -f " . join( " ",@collectSingl) ."\n"; #$out.srt
		}
		$cmd.= "rm -r $tmpP\n";
		#die $cmd."\n";
		my $cmd2 = "$secCogBin -i $outgz -DB $shrtDB -eval $MFopt{diaEVal} -percID $MFopt{DiaPercID} -minAlignLen $MFopt{DiaMinAlignLen} -minFractQueryCov $MFopt{DiaMinFracQueryCov} -mode 0 -LF $CLrefDBD/$refDB.length -reportDomains $getQSeq -DButil $CLrefDBD -tmp $tmpP";
		#$cmd2 .= " " if ($getQSeq);
		if ($curDB eq "ABR"){
			my $KrisABR = getProgPaths("KrisABR_scr");#"perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/secScripts/ABRblastFilter.pl";
			$cmd2 = "$KrisABR $outgz $outD/ABR/ABR.genes.txt $outD/ABR/ABR.cats.txt $CLrefDBD\n";
		} elsif ($curDB eq "PAB" && $MFopt{PABtaxChk}){ #NOG assignments
			$cmd2 .= " -NOGtaxChk $outD/dia.NOG.blast.srt ";
		}
		$cmd2 .= "\n";
		$cmd2 .= "rm -f $outgz\n" if ($MFopt{DiaRmRawHits});
		
		#check if the secondary routines for parsing blast out still need to be run
		if ($doInterpret){
			if (!(-e  "$out.gz.stone" || -e  "$out.srt.gz.stone" || -e  "$out.stone")){$doInterpret=1;
			} else {$doInterpret=0;
			}
		}
		my $jobName = $jdep;
		my $globDep = $MFopt{globalDiamondDependence}->{$curDB};
		$globDep = "" if ($MFopt{globalDiamondDependence}->{$curDB} eq "$shrtDB-1");
		
		my $memu = $MFopt{diamondMem} . "G"; my $tmpCmd;
		if (!-d $outD || !(-e "$out" || -e "$out.gz"|| -e "$out.srt.gz") ){ #diamond alignments
			$jobName = "_D$shrtDB$JNUM"; 
			my @preConstr = @{$QSBoptHR->{constraint}};
			push(@{$QSBoptHR->{constraint}}, $avx2Constr);
			my $tmpSHDD = $QSBoptHR->{tmpSpace};
			$QSBoptHR->{tmpSpace} = $HDDspace{diamond}; #set option how much tmp space is required, and reset afterwards
			($jobName, $tmpCmd) = qsubSystem($logDir."Diamo$shrtDB.sh",$cmd,$ncore,$memu,$jobName,$jdep.";".$globDep,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR);
			@{$QSBoptHR->{constraint}} = @preConstr;
			$QSBoptHR->{tmpSpace} = $tmpSHDD;
			$doInterpret=1;#run interpret step in any case
			
		} else {
			$jobName = $globDep;
		}
		if ($doInterpret){ #parsing of dia output
			my $jobName2 =  "_DP$shrtDB$JNUM";
			$memu = "30G";
			($jobName, $tmpCmd) = qsubSystem($logDir."Diamo_parse$shrtDB.sh",$cmd2,1,$memu,$jobName2,$jobName,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR) ;
			#die "bo";
			$progStats{$curDB}{DiaDBSearchIncomplete}++;#this needs to be done first..
		} else {
			$progStats{$curDB}{DiaDBSearchCompl}++;
		}
		if ($jobN2 eq ""){ $jobN2 = $jobName; } else {$jobN2 .= ";".$jobName;}
	}
	return ($jobN2,$clnCmd);
}
sub nopareil(){
	my ($ar1,$outD,$Gdir,$name,$jobd) = @_;
	my $numCore = 4;
	my @re1 = @{$ar1};
	
	my $npBin = getProgPaths("nonpareil");#"/g/bork5/hildebra/bin/nonpareil/nonpareil";
	my $globalNPD = $baseOut."NonPareil/";
	system "mkdir -p $Gdir" unless (-d "$globalNPD");

	my $sumOut = "$name.npo";
	my $cmd = "mkdir -p $outD\n";
	$cmd .= "$npBin -s $re1[0] -f fastq -t 20 -m 40000 -b $outD$name -n 10240 -i 0.1 -m 0.2\n";#-t $numCore  -o $sumOut
	$cmd .= "cp $outD/$sumOut $Gdir";
	#R part
	#source('/g/bork5/hildebra/bin/nonpareil/utils/Nonpareil.R');
	#Nonpareil.curve('$outD/$sumOut');
	my $jobName = "_NP$JNUM"; my $tmpCmd;
	if (!-d $outD || !-e "$outD/$sumOut"){
		my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0; 
		($jobName,$tmpCmd) = qsubSystem($logDir."NonPar.sh",$cmd,1,"42G",$jobName,$jobd,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR);
		$QSBoptHR->{tmpSpace} =$tmpSHDD;
	}
	return $jobName;
}
sub contigStatsOutputsComplete {
	my ($path,$assD,$subprts,$AssemblyGo,$smpl,$requireSupportCoverage) = @_;
	my $ContigStatsDir  = "$path/$preDIRs{dir_ContigStats}";
	my $CSfilesComplete = 1;
	my $primaryCoverageComplete = contig_stats_coverage_complete($ContigStatsDir, "Coverage");
	$CSfilesComplete = 0 if (!$primaryCoverageComplete && $map{$smpl}{hasPrimaryRds});
	$CSfilesComplete = 0  if ($MFopt{kmerPerGene} && $AssemblyGo && ! fileGZs("$ContigStatsDir/scaff.pergene.4kmer.pm5" ));
	my $supportCoverageComplete = contig_stats_coverage_complete($ContigStatsDir, "Cov.sup");
	$CSfilesComplete = 0 if ($requireSupportCoverage && !$supportCoverageComplete);
	$CSfilesComplete = 0  if ($subprts =~ m/F/ && ! fileGZs( "$assD/ContigStats//FMG/FMGids.txt" ));
	$CSfilesComplete = 0  if ($subprts =~ m/G/ && ! fileGZs( "$assD/ContigStats/GTDBmg/marker_genes_meta.tsv" ));
	return $CSfilesComplete;
}

sub runContigStats{
	my ($path,$jobd,$assD,$subprts,$immSubm,  $tmpD,$AssemblyGo,$Nthr,$smpl,$requireSupportCoverage) = @_;
	my $sepCtsScript = getProgPaths("sepCts_scr");#
	my $CSfilesComplete = contigStatsOutputsComplete(
		$path, $assD, $subprts, $AssemblyGo, $smpl, $requireSupportCoverage,
	);
	my $rawReadSet = sampleReadSet($curSmpl, "raw");
	my $readL = $rawReadSet->{samplReadLength};
	my $readLX = $rawReadSet->{samplReadLengthX};

	#die;
	#die "$CSfilesComplete";
	return ("","",0) if ($CSfilesComplete);
	
	print "Deferring Contig Stats until its assembly-group mapping is submitted\n"
		unless ($immSubm);
	$Nthr =1 unless ($subprts =~m/[FGEm]/);
	
	
	my $jobName = "_CS$JNUM"; $QSBoptHR->{LocationCheckStrg}=""; my $tmpCmd="";
	my $jobDep = "";
	#system "mkdir -p $cwd" unless ($cwd eq "" || -d $cwd);
	my $cmd = "";
	#$cmd .= "mkdir -p $cwd\n" unless ($cwd eq "" );
	$cmd .= "$sepCtsScript -inD $path -assD $assD -subparts $subprts -readLength $readL -readLengthSup $readLX -tmpD $tmpD -threads $Nthr -smplID $smpl";
	#$jobName = "_CS$JNUM"; 
	($jobDep,$tmpCmd) = qsubSystem($logDir."ContigStats.sh",$cmd,$Nthr,int(50)."G",$jobName,$jobd,"",$immSubm,[],$QSBoptHR);
	$tmpCmd = "" if ($immSubm);
	#die "$cmd\n$logDir.ContigStats.sh,$cmd,$Nthr,int(50/$Nthr).G,$jobName,$jobd,$cwd,$immSubm\n";
	return ($jobDep,$tmpCmd, 1);
}
sub calcCoverage2nd
{
	my ($cov,$gff,$RL,$cstNme,$jobd,$dirsHr) = @_;
	my $readCov_Bin =getProgPaths("readCov");
	my $jobName = ""; $QSBoptHR->{LocationCheckStrg}=""; my $tmpCmd="";
	my $cmd = "$readCov_Bin $cov $gff $RL";
	my $qdir = $logDir; $qdir = ${$dirsHr}{qsubDir} if (exists( ${$dirsHr}{qsubDir} ));
	if (!-s $cov.".pergene" || !-s $cov.".percontig" || !-s $cov.".median.percontig" ){
		$jobName = "_COV$JNUM";
		$jobName = "$cstNme"."_$JNUM" if ($cstNme ne "");
		#die "$cmd\n";
		if (${$dirsHr}{submit}){
			my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0; 
			($jobName,$tmpCmd) = qsubSystem($qdir."COV$cstNme.sh",$cmd,1,"20G",$jobName,$jobd,"",1,[],$QSBoptHR);
			$QSBoptHR->{tmpSpace} =$tmpSHDD;
		} else {
			$tmpCmd = $cmd;
		}
	} else {
		$jobName = $jobd;
	}
	return ($jobName,$tmpCmd);
}

sub checkDrives{
	my ($aref) = @_;
	my @locs = @{$aref};
	my $retStr = "\n####### BEGIN file location check ######\n";
	foreach my $llo (@locs){$retStr.="ls -l $llo > /dev/null \n";}
	$retStr.=" \nsleep 3\n"; my $cnt=0;
	foreach my $llo (@locs){ $cnt++;
		$retStr.="if [ ! -d \"$llo\" ]; then echo \'Location $cnt does not exist\'; exit 5; fi\n";
	}
	$retStr .= "####### END file location check ######\n";
	return $retStr;
}

# $sdmjN = cleanInput($cfp1ar,$cfp2ar,$sdmjN,$smplTmpDir);
sub cleanInput( $ $ $){
	my ($sdmjN,$saveD) = @_;
	my $libraries = ensureSeqSetLibraries(sampleReadSet($curSmpl, "raw"), $curSmpl);
	my $cmd = "";
	foreach my $library (@{$libraries}) {
		my @files = grep { defined($_) && $_ ne '' } @{$library->{files}}{qw(r1 r2 single)};
		my @temporary = grep { /\Q$saveD\E/ && -e $_ } @files;
		$cmd .= "rm -f ".join(" ", @temporary)."\n" if (@temporary);
	}

	#die $cmd."\n";
	my $jobName = $sdmjN;
	if ($cmd ne ""){
		print "Removing raw input fastqs..\n";
		# Must wait for $sdmjN (the SDM job consuming these files) to finish before
		# deleting them: SDM is submitted asynchronously, so an immediate `system`
		# call here raced the queued job and deleted its inputs before it ran.
		my $tmpCmd;
		my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0";
		($jobName, $tmpCmd) = qsubSystem($logDir."ClnUnzip.sh",$cmd,1,"1G","_PC$JNUM",$sdmjN,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR);
		$QSBoptHR->{tmpSpace} =$tmpSHDD;
	}
	return $jobName;
 }
 
 
#Gap Filler to refine scaffolds
sub GapFillCtgs{
	my($libraries,$scaffolds,$GFdir_a,$dep,$xtrTag) = @_; #.= "_GFI1";
	my $GFbin = getProgPaths("gapfiller");#"perl /g/bork5/hildebra/bin/GapFiller/GapFiller_n.pl";
	my $pairs = libraryPairs($libraries);
	my @inserts;
	my $prefi = "GF";
	my $numCore = 16;
	mkdir($GFdir_a.$prefi);
	my $log = $GFdir_a.$prefi."/GapFiller.log";
	#system("mkdir -p $GFdir_a$prefi");
	my $GFlib = ($GFdir_a."GFlib.opt");
	my @libFiles;
	for (my $i=0;$i<@{$pairs};$i++){
		push(@libFiles, ($pairs->[$i]{files}{r1}.",".$pairs->[$i]{files}{r2}));
		my $metadata = $pairs->[$i]{metadata} || {};
		push(@inserts, $metadata->{insert_size} || 450);
	}
	createGapFillopt($GFlib,\@libFiles,\@inserts);
#GF round 1
	my $cmd = "";
	$cmd .= "mkdir -p $GFdir_a$prefi\n";
	$cmd .= $GFbin . " -l $GFlib -s $scaffolds -m 75 -o 2 -r 0.7 -d 70 -t 10 -g 1 -T $numCore -b ".$prefi." -D $GFdir_a > $log \n";
	$cmd .= "rm -r $GFdir_a$prefi/alignoutput $GFdir_a$prefi/intermediate_results $GFdir_a$prefi/reads\n";
	#die $cmd."\n";
	if ($cmd ne ""){
		my $tmpCmd;
		my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0; 
		($dep, $tmpCmd) = qsubSystem($logDir."GapFill_ext_$xtrTag.sh",$cmd,$numCore,int(80)."G","_GFE$JNUM",$dep,"",1,[],$QSBoptHR);
		$QSBoptHR->{tmpSpace} =$tmpSHDD;
	}
	return $dep;
}

#scaffolding via mate pairs
sub scaffoldCtgs{
	my ($AsgHR,$ASG, $externalLibraries, $refCtgs, $tmpD1,$outD,$dep,$Ncore,$smplName,$spadesRef,$xtraTag) = @_;
	my $libraries = getRawLibrariesAssmGrp($AsgHR,$ASG,0);
	my $pairs = libraryPairs($libraries);
	my $bwt2Bin = getProgPaths("bwt2");#"/g/bork5/hildebra/bin/bowtie2-2.2.9/bowtie2";
	my $besstBin = getProgPaths("BESST");#"/g/bork3/home/hildebra/bin/BESST/./runBESST";
	my $externalPairs = libraryPairs($externalLibraries || []);
	return ("",$dep) if (@{$pairs} == 0 );
	my $tmpD = $tmpD1."/scaff/";
	my $bwtIdx = $refCtgs.$MFcontstants{bwt2IdxFileSuffix}; my $cmdDB="";my $chkFile = "";
	my $clnCmd = ""; my $spadesDir = ""; my $spadFakeDir = "";
	if ($spadesRef ){#&& -e ("$refCtgs/scaffolds.fasta")){#this is supposed to be the spades dir
		my $oldDir = "spades_ori";
		$spadFakeDir = $refCtgs;
		$clnCmd .= "mkdir -p $refCtgs/../$oldDir/; mv $refCtgs/* $refCtgs/../$oldDir/; mv $refCtgs/../$oldDir $refCtgs\n";
		$spadesDir = "$refCtgs$oldDir"; 
		$clnCmd .= "cp $spadesDir/smpls_used.txt $refCtgs\n";
		$refCtgs .= "$oldDir/scaffolds.fasta";
		$clnCmd .="\ngzip $spadesDir/*";
		$clnCmd .="\ngunzip $refCtgs";
		($cmdDB,$bwtIdx,$chkFile) = buildMapperIdx("$refCtgs",$Ncore,0,$MFopt{MapperProg});#$Ncore);
	} elsif (!-e $bwtIdx){
		($cmdDB,$bwtIdx,$chkFile) = buildMapperIdx("$refCtgs",$Ncore,0,$MFopt{MapperProg});
	}
	
	
	#my @rd1; my @rd2;
	my $algCmd = "$clnCmd\n$cmdDB\n";
	$algCmd .= "mkdir -p $tmpD\n";
	my @bams; my $cnt=0; my @insSiz; my @orientations;
	for (my $i=0;$i<@{$pairs};$i++){
		next unless (($pairs->[$i]{label} || '') =~ m/mate/i);
		#push @rd1,$rds[$i];push @rd2,$rds[$i+1];
		my $tmpOut = "$tmpD/tmpMateAlign$cnt.bam";
		my $tmpBAM = "$tmpD/tmpMateAlign$cnt.srt.bam";
		$algCmd .= "$bwt2Bin --no-unal --end-to-end -p $Ncore -x $bwtIdx -X $MFconfig{mateInsertLength} -1 $pairs->[$i]{files}{r1} -2 $pairs->[$i]{files}{r2} | $smtBin view -b -F 4 - > $tmpOut\n";
		$algCmd .= "$smtBin sort -@ $Ncore -T kk -O bam -o $tmpBAM $tmpOut; $smtBin index $tmpBAM\n";
		#$algCmd .= "$novosrtBin --ram 50G -o $tmpBAM -i $tmpOut \n";
		$algCmd .= "rm $tmpOut\n\n";
		my $metadata = $pairs->[$i]{metadata} || {};
		push(@bams,$tmpBAM);
		push(@insSiz,$metadata->{insert_size} || 10000);
		push(@orientations,$metadata->{orientation} || "fr");
		$cnt++;
		#print "sc  mat\n";
	}
	for (my $i=0;$i<@{$externalPairs};$i++){
		#push @rd1,$rds[$i];push @rd2,$rds[$i+1];
		my $tmpOut = "$tmpD/tmpMateAlign$cnt.bam";
		my $tmpBAM = "$tmpD/tmpMateAlign$cnt.srt.bam";
		$algCmd .= "$bwt2Bin --no-unal --end-to-end -p $Ncore -x $bwtIdx  -1 $externalPairs->[$i]{files}{r1} -2 $externalPairs->[$i]{files}{r2} | $smtBin view -b -F 4 - > $tmpOut\n";
		$algCmd .= "$smtBin sort -@ $Ncore -T kk -O bam -o $tmpBAM $tmpOut; $smtBin index $tmpBAM\n";
		#$algCmd .= "$novosrtBin --ram 50G -o $tmpBAM -i $tmpOut \n";
		$algCmd .= "rm $tmpOut\n\n";
		my $metadata = $externalPairs->[$i]{metadata} || {};
		push(@bams,$tmpBAM);
		push(@insSiz,$metadata->{insert_size} || 500);
		push(@orientations,$metadata->{orientation} || "fr");
		$cnt++;
		#print "sc  mat\n";
	}
	
	#die "\n\n\n$algCmd\n\n\n";
	return ("",$dep) if (@bams ==0 );
	#create bowtie2 mapping

	my $zcmd = "-z 5000";  #-z 10000";
	 my $ori = "--orientation ".join(" ",@orientations);
	#for (my $i=1;$i<@bams;$i++){ $ori .= " fr"}#$zcmd .= " 10000";
	my $cmd = "";
	my $bams = join(" ",@bams);
	system "mkdir -p $outD" unless (-d $outD);
	$cmd .= "mkdir -p $outD\n";
	$cmd .= "$besstBin $zcmd $ori -f $bams -o $outD -c $refCtgs\n"; #-q $Ncore <- unstable?
	my $jobName = $dep;
	$cmd .= "rm -r $tmpD\n";
	#cleanup2, fake spades result folder
	my $renameCtgScr = getProgPaths("renameCtg_scr");#"perl renameCtgs.pl";
	if ($spadesDir ne ""){
		my $assStatScr = getProgPaths("assStat_scr");#"perl assemblathon_stats.pl";
		my $sizFiltScr = getProgPaths("sizFilt_scr");#"perl sizeFilterFas.pl";
		$cmd .= "mv $outD/BESST_output/pass1/Scaffolds_pass1.fa $spadFakeDir/scaffolds.fasta\n";
		$cmd .= "$renameCtgScr $spadFakeDir/scaffolds.fasta $smplName\n";
		$cmd .= "$sizFiltScr $spadFakeDir/scaffolds.fasta $MFopt{scaffoldMinSize} 200\n";
		$cmd .= "$assStatScr -scaff_size $MFopt{scaffoldMinSize} $spadFakeDir/scaffolds.fasta > $spadFakeDir/AssemblyStats.500.txt\n";
		$cmd .= "$assStatScr $spadFakeDir/scaffolds.fasta > $spadFakeDir/AssemblyStats.ini.txt\n";
		$cmd .= "$assStatScr $spadFakeDir/scaffolds.fasta.filt > $spadFakeDir/AssemblyStats.txt\n\n";
		my ($cmdX) = buildMapperIdx("$spadFakeDir/scaffolds.fasta.filt",$Ncore,0,$MFopt{MapperProg});
		$cmd .= $cmdX."\n";

	}
	my $newScaffFNA = "$outD/BESST_output/pass2/Scaffolds_pass2.fa";
	$cmd .= "$renameCtgScr $newScaffFNA $smplName\n";
	$cmd .= "\ntouch $outD/scaffDone.sto\n" ;
	$cmd = "" if (-e "$outD/scaffDone.sto");
	
#die "scaff cmd $cmd\n";
	if ($cmd ne ""){
		my $tmpCmd;
		($dep, $tmpCmd) = qsubSystem($logDir."BesstScaff_ext_$xtraTag.sh",$algCmd.$cmd,$Ncore,int(50)."G","_BBE$JNUM$xtraTag",$dep,"",1,[],$QSBoptHR);
	}
	return ($newScaffFNA,$dep);
}

#preprocess mate pairs (use nxtrim on them and communicate results)
 sub check_mates($ $ $ $ $){
	my ($ar,$ifastasPre,$mateD,$doMateCln,$dep) = @_;
	my @mat = @{$ar};
	my $cmd = "";
	my @sarPre; my $mateC=0;
	my @mates;
	my $nxtrimBin = getProgPaths("nxtrim");#"/g/bork3/home/hildebra/bin/NxTrim/./nxtrim";

	foreach my $matp (@mat){
		system "mkdir -p $mateD" unless (-d $mateD);
		my @mateX = split /,/,$matp;
		push(@sarPre,"$mateD/mate.${mateC}.se.fastq.gz");
		push(@mates,"$mateD/mate.${mateC}_R1.unknown.fastq.gz","$mateD/mate.${mateC}_R2.unknown.fastq.gz","$mateD/mate.${mateC}_R1.mp.fastq.gz","$mateD/mate.${mateC}_R2.mp.fastq.gz");
		next if ( -e "$mateD/matesDone.sto");
		#--rf keeps reads in rf; --joinreads joins pe
		$cmd .= "$nxtrimBin --ignorePF --separate -1 $mateX[0] -2 $mateX[1] -O $mateD/mate.$mateC\n";
		#$nxtrimBin --stdout-mp -1 $rd1 -2 $rd2 | $bwaBin mem $refCtgs -p - > out.sam
		$ifastasPre .= ";$mateD/mate.${mateC}_R1.pe.fastq.gz,$mateD/mate.${mateC}_R2.pe.fastq.gz";
		$mateC++;
	}
	#die "@mates SDS\n";
	$cmd .= "touch $mateD/matesDone.sto\n";
	#die $cmd;
	$cmd = "" if ($doMateCln == 0 || $mateC==0);
	my $jobName = $dep;
	if ($cmd ne ""){
		my $tmpCmd;
		($jobName, $tmpCmd) = qsubSystem($logDir."mateClean.sh",$cmd,1,"1G","_NX$JNUM",$dep,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR);
	}
	return ($ifastasPre,\@sarPre,\@mates,$jobName) ;
 }

 
 #preprocess mate pairs (use nxtrim on them and communicate results)
 sub check_matesL($ $ $ $){
	my ($pa1,$pa2,$mateD,$doMateCln) = @_;
	my %ret;
	my $cmd = "";
	my $nxtrimBin = getProgPaths("nxtrim");#"/g/bork3/home/hildebra/bin/NxTrim/./nxtrim";

	$cmd .= "rm -rf $mateD; mkdir -p $mateD\n";
	 my $mateC=0;
	#foreach my $matp (@mat){
	if ($mateC > 0){die "mate pairs only supports single library\n";}
	system "mkdir -p $mateD" unless (-d $mateD);
	#my @mateX = split /,/,$matp;
	#--rf keeps reads in rf; --joinreads joins pe
	$cmd .= "$nxtrimBin --ignorePF --separate -1 $pa1 -2 $pa2 -O $mateD/mate.$mateC\n";
	$cmd .= "rm -f $pa1 $pa2\n";
	#$nxtrimBin --stdout-mp -1 $rd1 -2 $rd2 | $bwaBin mem $refCtgs -p - > out.sam
	$ret{pe1} = "$mateD/mate.${mateC}_R1.pe.fastq.gz"; $ret{pe2} = "$mateD/mate.${mateC}_R2.pe.fastq.gz";
	$ret{se} =  "$mateD/mate.${mateC}.se.fastq.gz";
	$ret{un1} = "$mateD/mate.${mateC}_R1.unknown.fastq.gz"; $ret{un2} = "$mateD/mate.${mateC}_R2.unknown.fastq.gz";
	$ret{mp1} = "$mateD/mate.${mateC}_R1.mp.fastq.gz"; $ret{mp2} = "$mateD/mate.${mateC}_R2.mp.fastq.gz";
	$mateC++;
	#}
	#die "@mates SDS\n";
	my $locStone = "$mateD/matesDone.sto";
	$cmd .= "touch $locStone\n";
	#die $cmd;
	$cmd = "" if ($doMateCln == 0 || $mateC==0);
	
	return (\%ret,$cmd,$locStone) ;
 }

 
 #determines what read types are present and starts cleaning of mates, if required
 sub get_ifa_mifa($ $ $ $ $ $ $){
	my ($ifastas,$ifastasS,$libInfoAr, $mateD, $doMateCln, $jdep, $singlReadMode) = @_;
	my @allFastas = split /;/,$ifastas;
	my @libInfo = @{$libInfoAr}; my @miSeqFastas = (); my $mcnt=0;
	my @mates;
	#die "@libInfo\n";
	foreach (@libInfo){
		if ($_ =~ m/.*miseq.*/i) {
			push (@miSeqFastas, $allFastas[$mcnt]);
			splice(@allFastas, $mcnt, 1);
		} elsif ($_ =~ m/.*mate.*/i){ #remove from process
			push (@mates, $allFastas[$mcnt]);
			splice(@allFastas,$mcnt,1);
		}
		$mcnt++;
	}
	$ifastas = join(";",@allFastas);
	my $mi_ifastas = join(";",@miSeqFastas);
	my $singleIfas = $ifastasS; my $matRef = [];
	#die "$ifastas\n$singleIfas\n\n";
	#moved to seedUnzip2tmp
	#($ifastas,$singleAddAr,$matRef,$jdep) = check_mates(\@mates,$ifastas,$mateD,$doMateCln,$jdep);

	return ($ifastas, $mi_ifastas, $singleIfas, $matRef, $jdep);
 }
 sub get_sdm_outf($ $ $ $ $ $){
	my ($ifastas, $mi_ifastas,$finD,$singlReadMode,$pairedReadMode, $useXtras) = @_;
	my @ret1;my @ret2;my @sret;
	my $fEnd = "fq";
 	if ($MFopt{gzipSDMOut}){
		$fEnd = "fq.gz";
	}
	my $baseFname = "filtered";
	$baseFname = "filtered.suppl" if ($useXtras);
	
	if ( $pairedReadMode){
		push(@ret1,$finD."$baseFname.1.$fEnd");push( @ret2, ($finD."$baseFname.2.$fEnd"));push(@sret, ($finD."$baseFname.singl.$fEnd"));
	} 
	if ($singlReadMode && !$pairedReadMode){ #only single reads avaialble, different file ending..
		push(@sret, $finD."$baseFname.s.$fEnd");
	}
	if ($mi_ifastas ne ""){
		push(@ret1,$finD."${baseFname}_mi.1.$fEnd");push( @ret2, ($finD."${baseFname}_mi.2.$fEnd"));push(@sret, ($finD."${baseFname}_mi.singl.$fEnd"));
	}

	return (\@ret1,\@ret2,\@sret);
}
#($mergRdsHsh,$mergJbN) = mergeReads($arp1,$arp2,$sdmjN,$smplTmpDir."merge_clean/");
sub mergeReads(){
	my ($jdep,$outdir,$doMerge,$runThis) = @_;

	my $cleanSeqSetHR = sampleReadSet($curSmpl, "clean");
	my $pairs = libraryPairs(readLibrariesByScope($cleanSeqSetHR, 'primary', 1, $curSmpl));

	my $flashBin = getProgPaths("flash");
	my $numCores = 8;
	my $outT = "sdmCln";
	my %ret = (mrg => "", pair1 => "", pair2 => "");
	#print "XSADS\n!$doMerge || !$MFconfig{readsRpairs} || !$runThis\n";
	if (!$doMerge || !$MFconfig{readsRpairs} || @{$pairs} == 0 || !$runThis){return ("");}
	#print "XSADS1\n";
	#if (@{$arp1} > 1 ){die "Array with reads provided to merging routine is too large!\n@{$arp1}\n";}
	my $mergCmd  = "";
	for (my $i=0; $i<@{$pairs};$i++){
		my $outTL = $outT;
		$outTL .= ".$i" if ($i > 0);
		$mergCmd .= "$flashBin -M 250 -z -o $outTL -d $outdir -t $numCores $pairs->[$i]{files}{r1} $pairs->[$i]{files}{r2}\n";
		if ($i > 0){
			$mergCmd .= "cat $outTL.extendedFrags.fastq.gz >> $outT.extendedFrags.fastq.gz;cat $outTL.notCombined_2.fastq.gz >> $outT.notCombined_2.fastq.gz; ";
			$mergCmd .= "cat $outTL.notCombined_1.fastq.gz >> $outT.notCombined_1.fastq.gz;\n";
		}
	}
	my $stone = "$outdir/$outT.sto";
	$mergCmd .="touch $stone\n ";
	my $jobName = "";
	if (-e $stone && -e "$outdir/$outT.extendedFrags.fastq" && !-e "$outdir/$outT.extendedFrags.fastq.gz"){
		#zip
		$mergCmd = "$pigzBin -f -p $numCores $outdir/$outT.extendedFrags.fastq $outdir/$outT.notCombined_1.fastq $outdir/$outT.notCombined_2.fastq\n";
		system "rm $stone";
	}
	if (!-e $stone){
		$jobName = "_FL$JNUM"; my $tmpCmd;
		($jobName, $tmpCmd) = qsubSystem($logDir."flashMrg.sh",$mergCmd,$numCores,"3G",$jobName,$jdep,"",1,[],$QSBoptHR);
	}
	
	#die "$logDir/flashMrg.sh\n$mergCmd\n";
	$ret{mrg} = "$outdir/$outT.extendedFrags.fastq.gz";
	$ret{pair1} = "$outdir/$outT.notCombined_1.fastq.gz";
	$ret{pair2} = "$outdir/$outT.notCombined_2.fastq.gz";
	
	$cleanSeqSetHR->{merged_library} = newReadLibrary(
		id => "$curSmpl:primary:merged", sample => $curSmpl, scope => 'primary',
		technology => libraryTechnology($pairs, "merged reads for $curSmpl", 1),
		is_long => 0, label => 'merged', phase => 'merged',
		files => {r1 => $ret{pair1}, r2 => $ret{pair2}, single => $ret{mrg}, bam => ''},
	);
	syncCleanSeqSetLegacy($cleanSeqSetHR);

	
	return ($jobName);
}
 
sub _shell_quote {
	my ($value) = @_;
	die "Cannot quote an undefined shell argument\n" unless defined $value;
	die "Shell argument contains a NUL or newline\n" if $value =~ /[\0\r\n]/;
	$value =~ s/'/'"'"'/g;
	return "'$value'";
}

sub _shell_command {
	return join(' ', map { _shell_quote($_) } @_);
}

sub _staged_read_files_present {
	my ($root, @markers) = @_;
	return 0 unless defined($root) && -d $root;
	my %marker = map { $_ => 1 } grep { defined($_) && $_ ne '' } @markers;
	my $present = 0;
	File::Find::find({
		no_chdir => 1,
		wanted => sub {
			return if $present || !-f $File::Find::name || $marker{$File::Find::name};
			$present = 1;
		},
	}, $root);
	return $present;
}

sub _validate_sdm_integer_setting {
	my ($name, $value, $minimum) = @_;
	$value = 0 unless defined $value;
	die "sdmClean::Invalid integer for $name: '$value'\n"
		unless $value =~ /^-?\d+$/ && $value >= $minimum;
	return int($value);
}

sub _set_sdm_option {
	my ($text, $name, $value) = @_;
	my $displayName = defined($name) ? $name : '<undefined>';
	die "Invalid SDM option name '$displayName'\n"
		unless defined($name) && $name =~ /^[A-Za-z][A-Za-z0-9_]*$/;
	die "Invalid value for SDM option '$name'\n"
		unless defined($value) && $value !~ /[\0\r\n\t]/;
	die "Cannot update '$name' in an undefined SDM option file\n" unless defined $text;
	my $replacement = "$name\t$value\n";
	my $matches = ($text =~ s/^\Q$name\E\t[^\r\n]*(?:\r?\n|\z)/$replacement/mg);
	$matches ||= 0;
	die "Expected one active '$name' entry in the SDM option file, found $matches\n"
		unless $matches == 1;
	return $text;
}

sub adaptSDMopt{
	my ($baseSF, $oDir, $readLength, $technology, $variant) = @_;
	$variant ||= '';
	die "Invalid SDM read length '$readLength'\n"
		unless defined($readLength) && $readLength =~ /^\d+$/;
	die "Invalid SDM option-file tag '$technology'/'$variant'\n"
		unless defined($technology) && "$technology$variant" =~ /^[A-Za-z0-9_.-]*$/;
	my $tag = join('_', grep { $_ ne '' } ($readLength, $technology, $variant));
	my $newSDMf = "$oDir/sdmo_$tag.txt";

	open my $inputFH, '<', $baseSF or die "Can't open SDM options '$baseSF': $!\n";
	my $str = do { local $/; <$inputFH> };
	close $inputFH or die "Can't close SDM options '$baseSF': $!\n";

	if ($technology ne 'proto' && $technology ne 'PB' && $technology ne 'ONT' && $readLength != 0){
		my $maxError = $readLength < 50 ? '0.5' : $readLength < 90 ? '1.2' : $readLength < 200 ? '2.5' : undef;
		$str = _set_sdm_option($str, 'maxAccumulatedError', $maxError) if defined $maxError;
		my ($windowWidth, $windowThreshold, $maxAmbiguous) = $readLength < 90
			? (8, 16, 1)
			: (18, 20, 2);
		$str = _set_sdm_option($str, 'TrimWindowWidth', $windowWidth);
		$str = _set_sdm_option($str, 'TrimWindowThreshhold', $windowThreshold);
		$str = _set_sdm_option($str, 'maxAmbiguousNT', $maxAmbiguous);
	}
	$str = _set_sdm_option($str, 'BinErrorModelAlpha', -1)
		unless $MFopt{sdmProbabilisticFilter};
	foreach my $option (sort keys %{$MFopt{sdm_opt}}){
		$str = _set_sdm_option($str, $option, $MFopt{sdm_opt}{$option});
	}

	my $temporary = "$newSDMf.$$";
	open my $outputFH, '>', $temporary or die "Can't write SDM options '$temporary': $!\n";
	print {$outputFH} $str or die "Can't write SDM options '$temporary': $!\n";
	close $outputFH or die "Can't close SDM options '$temporary': $!\n";
	rename $temporary, $newSDMf or die "Can't publish SDM options '$newSDMf': $!\n";
	return $newSDMf;
}



sub sdmOptSet{
	my ($samplReadLength, $technology) = @_;
	die "Invalid SDM read length '$samplReadLength'\n"
		unless defined($samplReadLength) && $samplReadLength =~ /^\d+$/;
	die "Cannot select SDM options without a sequencing technology\n"
		unless defined($technology) && $technology ne '';
	if ($MFopt{sdmOpt} ne ""){
		die "-customSDMopt must point to a file (currently: $MFopt{sdmOpt})\n"
			unless -f $MFopt{sdmOpt};
		die "-customSDMopt must point to a non-empty file (currently: $MFopt{sdmOpt})\n"
			unless -s $MFopt{sdmOpt};
		return ($MFopt{sdmOpt},$MFopt{sdmOpt});
	}

	my $pairOpt = $MFopt{baseSDMopt};
	$pairOpt = getProgPaths('baseSDMoptPacBio') if $technology eq 'PB';
	$pairOpt = getProgPaths('baseSDMoptONT') if $technology eq 'ONT';
	$pairOpt = $MFopt{baseSDMoptMiSeq} if $technology eq 'miSeq';
	$pairOpt = getProgPaths('baseSDMoptAVITI') if $technology eq 'AVITI';
	$pairOpt = getProgPaths('baseSDMoptProto') if $technology eq 'proto';
	my $singleOpt = $technology eq '454' ? getProgPaths('baseSDMopt454') : $pairOpt;

	# Always materialise an adapted file: global overrides and probabilistic-filter
	# settings also apply when read length is unknown or the library is long-read.
	if ($singleOpt eq $pairOpt){
		$pairOpt = adaptSDMopt($pairOpt, $MFglobal{globalLogDir}, $samplReadLength, $technology);
		$singleOpt = $pairOpt;
	} else {
		$pairOpt = adaptSDMopt($pairOpt, $MFglobal{globalLogDir}, $samplReadLength, $technology, 'pair');
		$singleOpt = adaptSDMopt($singleOpt, $MFglobal{globalLogDir}, $samplReadLength, $technology, 'single');
	}
	if (($technology eq 'PB' || $technology eq 'ONT') && $samplReadLength != 0 && $samplReadLength < 1000){
		print "WARNING: $technology reads have an unusually short configured length ($samplReadLength). "
			."Check -inputReadLength or -inputReadLengthSuppl.\n";
	}
	return ($pairOpt, $singleOpt);
}


sub sdmClean(){
	my ($curOutDir,$finD,$jobd,$runThis, $useXtras ) = @_;
	my $seqSet = sampleReadSet($curSmpl, "raw");
	my $cleanSeqSetHR = sampleReadSet($curSmpl, "clean");
	die "sdmClean::Could not find raw read set for $curSmpl\n"
		unless ref($seqSet) eq 'HASH';
	die "sdmClean::Could not find clean read set for $curSmpl\n"
		unless ref($cleanSeqSetHR) eq 'HASH';
	my $scope = $useXtras ? 'support' : 'primary';
	my $libraries = readLibrariesByScope($seqSet, $scope, 0, $curSmpl);
	return '' unless @{$libraries};
	$finD .= '/' unless $finD =~ m{/$};

	my %integerSetting = (
		sdmCores => _validate_sdm_integer_setting('sdmCores', $MFopt{sdmCores}, 1),
		# -1 is the established disabled sentinel; zero also emits no SDM flag.
		XfirstReads => _validate_sdm_integer_setting('XfirstReads', $MFconfig{XfirstReads}, -1),
		cut5pR1 => _validate_sdm_integer_setting('cut5pR1', $map{$curSmpl}{cut5pR1}, 0),
		cut5pR2 => _validate_sdm_integer_setting('cut5pR2', $map{$curSmpl}{cut5pR2}, 0),
		firstXrdsRd => _validate_sdm_integer_setting('firstXrdsRd', $map{$curSmpl}{firstXrdsRd}, 0),
		firstXrdsWr => _validate_sdm_integer_setting('firstXrdsWr', $map{$curSmpl}{firstXrdsWr}, 0),
	);

	my $samplReadLength = 0;
	$samplReadLength = $seqSet->{samplReadLength} if (defined($seqSet->{samplReadLength}));
	if ($useXtras){
		$samplReadLength = $seqSet->{samplReadLengthX} || 0;
	}

	my $sdmBin = getProgPaths("sdm");# sdm program from LotuS2 pipeline
	my $fEnd = $MFopt{gzipSDMOut} ? 'fq.gz' : 'fq';
	my $baseFname = $useXtras ? 'filtered.suppl' : 'filtered';
	my $logStem = $useXtras ? 'filterSuppl' : 'filter';
	my $stone = $finD.($useXtras ? 'filterSupplDone.stone' : 'filterDone.stone');
	my $sdmLogDir = "$logDir/sdm";
	my $cmd = _shell_command('mkdir', '-p', '--', $finD)."\n";
	$cmd .= _shell_command('rm', '-f', '--', $stone)."\n";
	$cmd .= _shell_command('mkdir', '-p', '--', $sdmLogDir)."\n";
	my @sdmCommon = (
		'-ignore_IO_errors', 1, '-i_qual_offset', 'auto',
		'-binomialFilterBothPairs', 1, '-threads', $integerSetting{sdmCores},
	);
	push @sdmCommon, ('-XfirstReads', $integerSetting{XfirstReads})
		if $integerSetting{XfirstReads} > 0;
	my @sdmExtra = $MFopt{SDMlogQualvsLen} ? ('-logLvsQ', 1) : ();
	my @sdmCut;
	push @sdmCut, ('-5PR1cut', $integerSetting{cut5pR1}) if $integerSetting{cut5pR1} > 0;
	push @sdmCut, ('-5PR2cut', $integerSetting{cut5pR2}) if $integerSetting{cut5pR2} > 0;
	push @sdmCut, ('-XfirstReadsRead', $integerSetting{firstXrdsRd}) if $integerSetting{firstXrdsRd} > 0;
	push @sdmCut, ('-XfirstReadsWritten', $integerSetting{firstXrdsWr}) if $integerSetting{firstXrdsWr} > 0;

	my @cleanLibraries;
	my @requiredOutputs;
	for (my $i = 0; $i < @{$libraries}; $i++) {
		my $library = $libraries->[$i];
		my $technology = $library->{technology} || "";
		checkSeqTech($technology, "MATAF4.pl::sdmClean library $library->{id}");
		die "MATAF4.pl::sdmClean library $library->{id} has no sequencing technology\n"
			if $technology eq '';
		my $isLong = $library->{is_long} || is3rdGenSeqTech($technology);
		my ($sdmPairOpt, $sdmSingleOpt) = sdmOptSet($samplReadLength, $technology);
		my @libraryArgs = @sdmCommon;
		push @libraryArgs, ('-illuminaClip', 1) if ($MFopt{trimAdapters} && !$isLong);
		my $suffix = $i == 0 ? '' : ".lib$i";
		my $prefix = "$finD$baseFname$suffix";
		my $hasPair = ($library->{files}{r1} || '') ne '';
		my $hasSingle = ($library->{files}{single} || '') ne '';
		my ($outR1, $outR2, $outSingle) = ('', '', '');
		my $logSuffix = $i == 0 ? '' : ".$i";

		if ($hasPair) {
			$outR1 = "$prefix.1.$fEnd";
			$outR2 = "$prefix.2.$fEnd";
			$outSingle = "$prefix.singl.$fEnd";
			$cmd .= _shell_command('rm', '-f', '--', $outR1, $outR2, $outSingle)."\n";
			$cmd .= _shell_command(
				$sdmBin, '-i', "$library->{files}{r1},$library->{files}{r2}",
				'-o_fastq', "$outR1,$outR2", '-options', $sdmPairOpt,
				'-paired', 2, @sdmExtra, '-log', "$sdmLogDir/$logStem$logSuffix.log",
				@libraryArgs, @sdmCut,
			)."\n";
			my $recoveredPattern = "$prefix.*.singl.$fEnd";
			$cmd .= 'mapfile -t recovered < <(compgen -G '._shell_quote($recoveredPattern)." || true)\n";
			$cmd .= 'if (( ${#recovered[@]} )); then cat -- "${recovered[@]}" > '
				._shell_quote($outSingle).'; rm -f -- "${recovered[@]}"; else ';
			$cmd .= $MFopt{gzipSDMOut}
				? _shell_command($pigzBin, '-c').' </dev/null > '._shell_quote($outSingle)
				: ': > '._shell_quote($outSingle);
			$cmd .= "; fi\n";
			push @requiredOutputs, $outR1, $outR2, $outSingle;
		}

		if ($hasSingle) {
			$outSingle = "$prefix.s.$fEnd" unless $hasPair;
			my $tmpSingle = "$prefix.input-single.fq";
			$cmd .= _shell_command('rm', '-f', '--', $tmpSingle)."\n";
			$cmd .= _shell_command(
				$sdmBin, '-i', $library->{files}{single}, '-o_fastq', $tmpSingle,
				'-options', $sdmSingleOpt, @sdmExtra, '-paired', 1,
				'-log', "$sdmLogDir/$logStem.S$logSuffix.log", @libraryArgs, @sdmCut,
			)."\n";
			my $redirect = $hasPair ? '>>' : '>';
			if ($MFopt{gzipSDMOut}) {
				$cmd .= _shell_command($pigzBin, '-p', $integerSetting{sdmCores}, '-c', $tmpSingle)
					." $redirect "._shell_quote($outSingle)."\n";
			} else {
				$cmd .= _shell_command('cat', '--', $tmpSingle)." $redirect "._shell_quote($outSingle)."\n";
			}
			$cmd .= _shell_command('rm', '-f', '--', $tmpSingle)."\n";
			push @requiredOutputs, $outSingle unless $hasPair;
		}

		push @cleanLibraries, newReadLibrary(
			id => $library->{id}, sample => $library->{sample} || $curSmpl,
			scope => $scope, technology => $technology, is_long => $isLong,
			label => $library->{label}, phase => 'clean',
			files => {r1 => $outR1, r2 => $outR2, single => $outSingle, bam => ''},
			source_files => $library->{files},
			metadata => {source_library_id => $library->{id}},
		);
	}

	if (!-e "$curOutDir/input_fil.txt" && !$useXtras) {
		open my $inputFH, '>', "$curOutDir/input_fil.txt" or die "Cannot write $curOutDir/input_fil.txt: $!\n";
		print {$inputFH} join(';', map {
			$_->{files}{r1} ? "$_->{files}{r1},$_->{files}{r2}" : $_->{files}{single}
		} @{$libraries});
		close $inputFH or die "Cannot close $curOutDir/input_fil.txt: $!\n";
	}

	$cmd .= _shell_command('touch', '--', $stone)."\n";
	my $jobName = "";
	my $presence = -e $stone ? 1 : 0;
	$presence = 0 if grep { !-e $_ } @requiredOutputs;
	my $qsubFile = $logDir."sdmReadCleaner.sh";
	replaceScopeLibraries($cleanSeqSetHR, $scope, \@cleanLibraries, 1, $curSmpl);
	$qsubFile = $logDir."sdmReadCleanerSuppl.sh" if ($useXtras);

	if (!$presence && $runThis){
		print "sdm'ing support reads..\n" if ($useXtras);
		$jobName = "_SDM${useXtras}_$JNUM"; my $tmpCmd;
		local $QSBoptHR->{tmpSpace} = 0;
		($jobName, $tmpCmd) = qsubSystem($qsubFile,$cmd,$integerSetting{sdmCores},$MFopt{sdmMem},$jobName,$jobd,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR);
	}
	return $jobName;
}


sub mocat_reorder(){
	my($ar1,$ar2,$singlAr,$inJob) = @_;
	
	#my @pairs = split(";",$ifastas);
	my @ret1 = @{$ar1}; my @ret2 = @{$ar2}; 
	#print "@ret1\n@ret2\n";
	#foreach (@pairs){		my @spl = split /,/;		push(@ret1,$spl[0]);push(@ret2,$spl[1]);	}
	#my $jobName = $inJob;
	return (\@ret1,\@ret2,$singlAr,$inJob);
}
sub SEEECER(){
	my ($p1ar,$p2ar,$tmpD) = @_;
	die("SEECER deactived\n");
	my @p1 = @{$p1ar}; my @p2 = @{$p2ar};
	if ( $p1[0] =~ m/.*\.gz/){
		system("gunzip -c ".join(" ",@p1)." > $tmpD/pair.1.fastq");	system("gunzip -c ".join(" ",@p2)." > $tmpD/pair.2.fastq");
		@p1 = ("$tmpD/pair.1.fastq"); @p2 = ("$tmpD/pair.2.fastq");
	} elsif ( @p1 > 1 ){
		system("cat ".join(" ",@p1)." > $tmpD/pair.1.fastq");	system("cat ".join(" ",@p2)." > $tmpD/pair.2.fastq");
		@p1 = ("$tmpD/pair.1.fastq"); @p2 = ("$tmpD/pair.2.fastq");
	}
	my $SEEbin = "bash /g/bork5/hildebra/bin/SEECER-0.1.3/SEECER/bin/run_seecer.sh";
	system("mkdir -p $tmpD/tmpS");
	my $cmd = $SEEbin . " -t $tmpD/tmpS $p1[0] $p2[0]";
	#qsubSystem($logDir."SEECERCleaner.sh",$cmd,1,"30G",1);
	my @ret1= ($tmpD."pair.1.fastq_corrected.fa");my @ret2= ($tmpD."pair.2.fastq_corrected.fa");
	return (\@ret1,\@ret2);
}


#calculates md5 sums needed to upload filtered raw fastq's to EBI/SRA
sub unploadRawFilePostprocess{
	#die "XXn\n";
	if (@EBIjobs == 0){return;}
	if ($MFconfig{uploadRawRds} eq ""){return ;}
	if (-s "$MFconfig{uploadRawRds}/R1.md5"){return;}
	print "Postprocess upload postprocess\n";
	#md5sum --version
	
	my $cmd = "md5sum $MFconfig{uploadRawRds}/*.R2.fq.gz > $MFconfig{uploadRawRds}/R2.md5 ";
	$cmd .= "& \n md5sum $MFconfig{uploadRawRds}/*.R1.fq.gz > $MFconfig{uploadRawRds}/R1.md5 ";
	unless (my @files = glob("\Q$MFconfig{uploadRawRds}\E/*.Rsingl.fq.gz")) {
		$cmd .= "###\n### in case needed:\n### md5sum $MFconfig{uploadRawRds}/*.Rsingl.fq.gz > $MFconfig{uploadRawRds}/Rsingl.md5\n";
	} else {
		$cmd .= "& \nmd5sum $MFconfig{uploadRawRds}/*.Rsingl.fq.gz > $MFconfig{uploadRawRds}/Rsingl.md5\n";
	}
	
	my ($jobN, $tmpCmd) = qsubSystem("$MFglobal{globalLogDir}/postEBI.sh",$cmd,3,"20G","_PP",join(";",@EBIjobs),"",1,$QSBoptHR->{General_Hosts},$QSBoptHR) ;
}

#cleans raw fastqs fastq's (human DNA) to prepare upload to EBI/SRA
sub uploadRawFilePrep{
	if ($MFconfig{uploadRawRds} eq ""){return "" ;}

	my ($tmpD,$smplID, $jdep, $useXtras) = @_;
	my $totalInputSizeMB = $useXtras
		? ($map{$smplID}{inputXFileSizeMB} || 0)
		: ($map{$smplID}{inputFileSizeMB} || 0);

	my $seqSet = sampleReadSet($smplID, "raw");
	my $scope = $useXtras ? 'support' : 'primary';
	my $libraries = readLibrariesByScope($seqSet, $scope, 0, $smplID);
	return "" unless @{$libraries};
	my $tag = ""; $tag = "X." if ($useXtras);
	if ($useXtras){
		print "preparing xtra raw fastq upload for ENA/SRA (hostfilter $MFopt{humanFilter}).. ";
	} else {
		print "preparing raw fastq upload for ENA/SRA (hostfilter $MFopt{humanFilter}).. ";
	}
	
	my $fastqhdsChk = getProgPaths("fastHdChkENA");

	#prepare databases
	my $DBdir = $MFglobal{krakenDBDirGlobal}."/";
	my @DBname = ("hum1stTry"); 
	if (@filterHostDB > 0){
		@DBname = (); my $cnt=0;
		foreach my $fhdb (@filterHostDB){
			if (!-d $fhdb){die "-filterHostKrakDB is not a dir: \"$fhdb\"\n";}
			$DBname[$cnt] = $fhdb; $DBdir=""; $cnt++;
		}
	}


	my $outD  ="$MFconfig{uploadRawRds}"; 
	my $numThr = 4;
	system "mkdir -p $outD/tmp/ " unless (-d "$outD/tmp");
	$tmpD = "$tmpD/tmp$tag/";
	my $cmd = "";#"rm -rf $outD/tmp/;mkdir -p $outD/tmp/\n";
	for (my $idx = 0; $idx < @{$libraries}; $idx++) {
		my $library = $libraries->[$idx];
		my $technology = $library->{technology} || '';
		checkSeqTech($technology, "MATAF4.pl::uploadRawFilePrep library $library->{id}") if $technology ne '';
		my $isLong = $library->{is_long} || is3rdGenSeqTech($technology);
		my $xtra = $tag."$idx.";
		$xtra = "mate.$xtra" if (($library->{label} || '') =~ /mate/i);
		$xtra = "miSeq.$xtra" if (($library->{label} || '') =~ /miseq/i);
		$xtra = "PB.$xtra" if ($technology eq 'PB');
		$xtra = "ONT.$xtra" if ($technology eq 'ONT');

		if ($library->{files}{r1}) {
			my ($r1, $r2) = @{$library->{files}}{qw(r1 r2)};
			die "Upload library $library->{id} mixes compressed and uncompressed mates\n"
				if (($r1 =~ /\.gz$/) != ($r2 =~ /\.gz$/));
			my $gz = $r1 =~ /\.gz$/ ? '.gz' : '';
			my $of1 = "$tmpD/$smplID.${xtra}R1.fq$gz";
			my $of2 = "$tmpD/$smplID.${xtra}R2.fq$gz";
			my $tmp1 = "$tmpD/tmp/$smplID.${xtra}TEMP.R1.fq$gz";
			my $tmp2 = "$tmpD/tmp/$smplID.${xtra}TEMP.R2.fq$gz";
			my $final1 = "$outD/".basename($of1);
			my $final2 = "$outD/".basename($of2);
			unless (-e $final1 && -e $final2) {
				$cmd .= "rm -fr $tmpD;\nmkdir -p $tmpD/tmp/ $outD\n";
				$cmd .= "ln -s $r1 $tmp1\nln -s $r2 $tmp2\n";
				$cmd .= hostRmBase($tmp1,$tmp2,$MFopt{humanFilter},$isLong,$numThr,$tmpD,"$DBdir$DBname[0]");
				$cmd .= "$fastqhdsChk $tmp1 1; $fastqhdsChk $tmp2 2;\n";
				$cmd .= "mv $tmp1 $of1\nmv $tmp2 $of2\nmv $of1 $of2 $outD\n\n";
			}
		}

		if ($library->{files}{single}) {
			my $single = $library->{files}{single};
			my $gz = $single =~ /\.gz$/ ? '.gz' : '';
			my $of = "$tmpD/$smplID.${xtra}Rsingle.fq$gz";
			my $tmp = "$tmpD/tmp/$smplID.${xtra}TEMP.Rsingle.fq$gz";
			my $final = "$outD/".basename($of);
			unless (-e $final) {
				$cmd .= "rm -fr $tmpD;\nmkdir -p $tmpD/tmp/ $outD\n";
				$cmd .= "ln -s $single $tmp\n";
				$cmd .= hostRmBase($tmp,"",$MFopt{humanFilter},$isLong,$numThr,$tmpD,"$DBdir$DBname[0]");
				$cmd .= "$fastqhdsChk $tmp 3;\n";
				$cmd .= "mv $tmp $of\nmv $of $outD\n\n";
			}
		}
	}
	$cmd .= "rm -rf $tmpD\n" if ($cmd ne "");
	
	#$cmd = "rm -rf /local/hildebra/MF/\n";
	
	#die "$cmd\n";
	my $retJob = "";
	if ($cmd ne ""){
		#systemW $cmd;
		my $preHDDspace=$QSBoptHR->{tmpSpace};
		$QSBoptHR->{tmpSpace} = int($totalInputSizeMB/1024*6)+30  ."G";
#		$QSBoptHR->{tmpSpace} = $HDDspace{prepPub}; #increase local space..  # use $totalInputSizeMB ??
		#print "$QSBoptHR->{tmpSpace}\n";
		my ($jobN, $tmpCmd) = qsubSystem("$logDir/prepEBI$tag.sh",$cmd,$numThr,int(50) . "G","EBI$tag$JNUM","$jdep;$krakDeps","",1,$QSBoptHR->{General_Hosts},$QSBoptHR) ;
		$QSBoptHR->{tmpSpace} = $preHDDspace;
		$retJob = $jobN;
	} else {
		print "all done\n";
	}
	return $retJob;
}
#function to copy files into mocat compatible format .. due to mocat inflexibility for input *sic*
sub mocatFileCpy($$$$$){
	if ($MFconfig{mocatLinkDir} eq ""){return 0 ;}
	my ($inD,$cfp1ar,$cfp2ar,$cfpsar,$smplID) = @_;
	my @pa1 = @{$cfp1ar}; my @pa2 = @{$cfp2ar}; my @sa = @{$cfpsar};
	my $outD  ="$MFconfig{mocatLinkDir}/$smplID"; 
	my $cmd = "";
	$cmd = "mkdir -p $outD\n" unless (-d $outD);
	my $idx=0;
	foreach my $f (@pa1){
		my $of = "$outD/raw.$idx.1.fq"; $of .= ".gz" if ($f =~ m/\.gz$/);
		$cmd .= "ln -s $inD$f $of\n" unless (-e $of);
		$idx++;
	}
	$idx=0;
	foreach my $f (@pa2){
		my $of = "$outD/raw.$idx.2.fq"; $of .= ".gz" if ($f =~ m/\.gz$/);
		
		$cmd .= "ln -s $inD$f $of\n"  unless (-e $of);
		$idx++;
	}
	$idx=0;
	foreach my $f (@sa){
		my $of  = "$outD/raw.$idx.single.fq"; $of .= ".gz" if ($f =~ m/\.gz$/);
		$cmd .= "ln -s $inD$f $of\n"  unless (-e $of);
		$idx++;
	}
	
	#die "$cmd\n";
	if ($cmd ne ""){
		systemW $cmd;
		#my ($jobN, $tmpCmd) = qsubSystem("$logDir/cp2mocat.sh",$cmd,1,"1G","_CPM$JNUM",$jdep,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR) ;
	}
	return 1;
}



sub complexGunzCpMv($ $ $ $ $ $){
	my ($fastap,$in,$scrathD,$finDest,$ncore, $allowLinks) = @_;
	my $unzipcmd = "";
	my $lowEffort=1;
	
	my $out = basename($in);
	if (0&&$in =~ m/\.gz$/){ #deactivated.. takes too much space on file sys
		$unzipcmd .= "$pigzBin  -f -p $ncore -d -c $fastap$in";
		$out =~ s/\.gz$//;
		$unzipcmd .= " > $finDest$out \n";
		$lowEffort=0;
	}elsif ($in =~ m/\.bz2/){
		my $bzip2B  = getProgPaths("bzip2");
		
		$out =~ s/\.bz2$/\.gz/;
		$unzipcmd .= "$bzip2B --decompress --stdout $fastap$in | $pigzBin  -f -p $ncore -c -f > $finDest$out";
		$lowEffort=0;
	}elsif ($allowLinks){
		$unzipcmd .= "ln -s $fastap$in $finDest/$out;\n";
		
	} else {
		$unzipcmd .= "cp $fastap$in $finDest;\n";
		$lowEffort=0;
	}
	if (0 && $scrathD ne $finDest){
		$unzipcmd .= "mv -f $scrathD$in $finDest\n";
		$lowEffort=0;
		die "$scrathD ne $finDest";
	}
	my $newRDf = $finDest."$out";
	
	return ($unzipcmd,$newRDf,$lowEffort);
}

sub outfiles_trimall($ $){
	my ($opath,$fil) = @_;
	my $OFp1 = "$opath/$fil"; 
	my $OFu1 = "$opath/$fil";	
	if ($OFu1 =~ m/\.gz$/){
		$OFu1 =~ s/(\.[^\.]+\.gz)$/\.sing$1/ ;#$OFp1 =~ s/\.gz$// 
	} elsif ($OFu1 =~ m/\.f[^\.]*q$/){
		$OFu1 =~ s/(\.f[^\.]*q)$/\.sing$1/ ;
	} else {
		die "Unknown file ending: $fil\n";
	}
	#die "$OFp1\n$OFu1\n";
	return ($OFp1,$OFu1);
}
sub outfiles_Bam($ $){
	my ($opath,$fil) = @_;
	my $OFu1 = "$opath/$fil"; 
	if ($OFu1 =~ m/\.bam$/){
		$OFu1 =~ s/(\.bam)$/\.unbam\.fq\.gz/ ;
	} else {
		die "Unknown file ending: $fil\n";
	}
	#die "$fil\n$OFu1\n";
	return $OFu1;
}


sub valid_files{
	my ($path, $paIR) = @_;
	my @paI = @{$paIR};
	my @paO;
	foreach my $ff (@paI){
		if (-e "$path/$ff" && !-d "$path/$ff"){
			push (@paO, $ff) ;
			print "$path/$ff\n";
		}
	}
	return \@paO;
}



sub seedUnzip2tmp{
	my ($fastp,$curSmpl,$jDepe,$tmpPath,$finDest, 
		$calcUnzp,$finalMapDir,$porechopFlag,$inputRawFile) = @_;
	my $himipeSeqAd = getProgPaths("illuminaTS3pe"); #for trimomatic
	my $trimJar = getProgPaths("trimomatic");

	my $smplPrefix = $map{$curSmpl}{prefix};
	my $doMateCln = 1; #nxtrim  mate pairs ?
	my $numCore=$MFopt{unzipCores};
	my $rawReads=""; my $mmpu = "";
	my $inputDiscovery = discoverSampleInputs($curSmpl, $fastp);
	my $xtraRdsTech = $inputDiscovery->{support_technology};
	my $totalInputSizeMB=10000; #default to something sensible
	my $totalXInputSizeMB=0; #normally no suppl present..

	if ($fastp eq "") { print "No primary dir.. \n"; }
	#main source for input files..
	my @pa1 = (); my @pa2 = (); my @pas = (); 
	#support reads.. 
	my @paX1= (); my @paXs= (); my @paX2= ();
	#information on libraries (name of library)
	my @libInfo= (); my @libInfoX= ();
	#information on seq tech
	my $seqTech = "hiSeq"; my $seqTechX = ""; #assume by defauly illumina as read tec
	#store bam formated input..
	my @paBam; my @paBamX;
	
	if (exists($map{$curSmpl}{"SeqTech"})){$seqTech = $map{$curSmpl}{"SeqTech"};}
	
	
		#die "Mapping $pa1 to ref\n";
	if ($xtraRdsTech ne ""){
		checkSeqTech($xtraRdsTech,"MATAF4.pl::SupportReads");
		$seqTechX = $xtraRdsTech;
	}
	my $is3rdGen = is3rdGenSeqTech($seqTech);
	my $is3rdGenX = is3rdGenSeqTech($seqTechX);
	my $singleSeqTech = $map{$curSmpl}{SeqTechSingl} || $seqTech;
	checkSeqTech($singleSeqTech,"MATAF4.pl::singleton reads");
	
	#die "$seqTech\n$is3rdGen\n";
	
	#create empty return object
	my %seqSet = (libraries => [], pa1 => \@pa1, pa2 => \@pa2, pas => \@pas, seqTech => $seqTech, is3rdGen => $is3rdGen,
			paX1 => \@paX1, paX2 => \@paX2, paXs => \@paXs, seqTechX => $seqTechX, is3rdGenX =>  $is3rdGenX,
			libInfo => \@libInfo, libInfoX => \@libInfoX,
			totalInputSizeMB => -1, totalXInputSizeMB => -1,
			rawReads => $rawReads,
			mmpu => $mmpu, 
			);

	
	#check if unmapped reads requested..
	if ($MFopt{useUnmapped}){
		die "useUnmapped needs checking before use.. Exiting\n";
		$seqSet{libraries} = [newReadLibrary(
			id => "$curSmpl:primary:unmapped", sample => $curSmpl, scope => 'primary',
			technology => $seqTech, is_long => $is3rdGen, label => 'unmapped', phase => 'staged',
			files => {r1 => "$finalMapDir/unaligned/unal.1.fq.gz", r2 => "$finalMapDir/unaligned/unal.2.fq.gz",
				single => "$finalMapDir/unaligned/unal.fq.gz", bam => ''},
		)];
		syncSeqSetLegacy(\%seqSet);
		$seqSet{totalInputSizeMB} = filsizeMB((@pa1,@pa2,@pas));
#		return ("",\@pa1,\@pa2, \@pas, 0, "", "",\@libInfo, "",$totalInputSizeMB);
		return ("", \%seqSet);
	}
	

	if ($inputDiscovery->{primary_error} ne "") {
		if ($inputDiscovery->{primary_missing_dir} && !$MFconfig{abortOnEmptyInput}) {
			print $inputDiscovery->{primary_error};
			$seqSet{totalInputSizeMB}=0;
			return ("EMPTY_DO_NEXT", \%seqSet);
		}
		die $inputDiscovery->{primary_error};
	}

	die $inputDiscovery->{support_error}
		if ($inputDiscovery->{support_error} ne "");
	@pa1 = @{$inputDiscovery->{primary}{read1}};
	@pa2 = @{$inputDiscovery->{primary}{read2}};
	@pas = @{$inputDiscovery->{primary}{single}};
	@paBam = @{$inputDiscovery->{primary}{bam}};
	@paX1 = @{$inputDiscovery->{support}{read1}};
	@paX2 = @{$inputDiscovery->{support}{read2}};
	@paXs = @{$inputDiscovery->{support}{single}};
	@paBamX = @{$inputDiscovery->{support}{bam}};
	
	#create libinfo array 
	for (my $i = 0; $i < max((scalar(@pas),scalar(@pa1))); $i++) {
		$libInfo[$i] = "lib$i";
	}
	if(@paBam > 0){ #assumes as input a mix of fq's and bams
		my $firstBamLibrary = scalar(@libInfo);
		my $lastBamLibrary = $firstBamLibrary + scalar(@paBam);
		for (my $i = $firstBamLibrary; $i < $lastBamLibrary; $i++) {
			$libInfo[$i] = "lib$i";
		}
	}
	#$locStats{hasPaired} = 1 if (@pa1 > 0);$locStats{hasSingle} = 1 if (@pas > 0);
	if ($xtraRdsTech ne ""){
		@libInfoX = ($xtraRdsTech) x max(scalar(@paX1), scalar(@paXs) + scalar(@paBamX));
	}


	$totalInputSizeMB = $inputDiscovery->{primary_bytes} / (1024 * 1024);
	$map{$curSmpl}{inputFileSizeMB} = $totalInputSizeMB;
	#and file size for suppl files..
	$totalXInputSizeMB = $inputDiscovery->{support_bytes} / (1024 * 1024);
	$map{$curSmpl}{inputXFileSizeMB} = $totalXInputSizeMB;


	#die "@pa1";
	#die "For dir $fastp, unual fastq files exist:P1\n".join("\n",@pa1)."\nP2:\n".join("\n",@pa2)."\n";
	
	#check if raw file is a symlink and if this is valid & create raw read link (HD times):
	for (my $i = 0; $i<@pa1; $i++){
		my $pp = $fastp;
		my $read1Path = "$pp/$pa1[$i]";
		my $read2Path = "$pp/$pa2[$i]";
		my $realP = $read1Path;
		if (-l $read1Path){
			$realP = abs_path($read1Path) || "";
			if ($realP eq '' || !-f $realP){die "File $pa1[$i] is not file.\n";}
		}
		if (-l $read2Path){
			my $realMate = abs_path($read2Path) || "";
			if ($realMate eq '' || !-f $realMate){die "File $pa2[$i] does not exist.\n";}
		}
#		print "$pp$pa1[$i]\n";
		if (defined $realP && $realP =~ m/\/(MMPU[^\/]+)\//){
			$mmpu = $1;
		}
	}
	#die $rawReads."\n";
	#DEBUG fix to reduce file sizes
	#@fastap2 = ($fastap2[0]); @fastap1 = ($fastap1[0]);
	
	#could not find any files.. report this to user
	if (@pa1 == 0 && @pas ==0 && @paBam == 0 && @paX1 == 0 && @paXs ==0 && @paBamX ==0 ){
		my $msg = "Can;t find files in $fastp\nUsing search pattern: $smplPrefix$MFconfig{rawFileSrchStr1}  $smplPrefix$MFconfig{rawFileSrchStr2}\n$smplPrefix$MFconfig{rawFileSrchStrSingl}\n";
		die $msg if ($MFconfig{abortOnEmptyInput});
		print $msg;$totalInputSizeMB=0;$map{$curSmpl}{inputFileSizeMB} =0;$map{$curSmpl}{inputXFileSizeMB} =0;
		#return ("EMPTY_DO_NEXT",\@pa1,\@pa2, \@pas, 0, "", "",\@libInfo, "",$totalInputSizeMB);
		$seqSet{libraries} = readLibrariesFromArrays(
			sample => $curSmpl, scope => 'primary', phase => 'raw', technology => $seqTech,
			is_long => $is3rdGen, r1 => \@pa1, r2 => \@pa2, single => \@pas, labels => \@libInfo,
		);
		syncSeqSetLegacy(\%seqSet);
		$map{$curSmpl}{inputFilesEmpty} = 1;
		return ("EMPTY_DO_NEXT", \%seqSet);

	}
	
	#set up stone to mark end of run
	my $finishStone = "$finDest/rawRds/done.sto";
	my $trimoStone = "";
	my $porechStone = "";
	if ($porechopFlag){
		$porechStone = "$finDest/rawRds/poreChopped.stone";
	}
	if (!$porechopFlag ){
		$MFconfig{filterFromSource} = 0;#too dangerous..
	}
	
	#report on what was found so far..
	if ($map{$curSmpl}{inputFileSizeMB} > 0 || $map{$curSmpl}{inputXFileSizeMB} > 0){
		printf("Raw primary input size: %.1f Mb", $map{$curSmpl}{inputFileSizeMB});     
		printf("    Suppl: %.1f Mb", $map{$curSmpl}{inputXFileSizeMB})
			if (@paX1 || @paXs || @paBamX);
		#print "Input size raw (Mb): " . int($map{$curSmpl}{inputFileSizeMB}) ;
		#print "; Suppl: " . int($map{$curSmpl}{inputXFileSizeMB} );
		if (@pa1 || @pas){
			print " Fastq pairs: " . scalar(@pa1) if (@pa1);
			print " Fastq Singls: " . scalar (@pas) if (@pas);
			print " tech: $seqTech"; if($is3rdGen){print ", 3rd gen " ;} else {print " ";}
		}
		print " Fastq Supports Singls: " . scalar (@paXs) if (@paXs);
		print " Bam Singls: " . scalar (@paBam) if (@paBam);
		print " Bam Supports Singls: " . scalar (@paBamX) if (@paBamX);
		
		
		print "\n" ;
	}
	
	
	#relinking in mocat file structure, if requested ## not used any longer
	#my $mocatFCDone = mocatFileCpy($fastp,\@pa1,\@pa2,\@pas,$curSmpl);
	
	#prepare files to be uploaded to EBI etc, if requested
	#my $uplDone = uploadRawFilePrep($fastp,$tmpPath,\@pa1,\@pa2,\@pas,$curSmpl,\@libInfo,$totalInputSizeMB);
	#push (@EBIjobs, $uplDone);

	#no longer used.. just process files as if real files, even if upload2EBI flag used..
#	if (0 ||                              $mocatFCDone || $uplDone ne ""){ 
#		$totalInputSizeMB=0;$map{$curSmpl}{inputFileSizeMB}=0;
		#return ("EMPTY_DO_NEXT",\@pa1,\@pa2, \@pas, 0, "", "",\@libInfo, "",$totalInputSizeMB);
#		$seqSet{pa1} = \@pa1;$seqSet{pa2} = \@pa2;$seqSet{pas} = \@pas;
#		$seqSet{libInfo} = \@libInfo;
#		return ("EMPTY_DO_NEXT", \%seqSet);
#	}
	# Preserve the authoritative read locations before the arrays are rewritten
	# to their generated rawRds destinations. A retry must validate the sources,
	# because missing rawRds files are exactly what this stage recreates.
	my @sourceInputs = @{source_input_files($fastp, @pa1, @pa2, @pas, @paBam)};
	push @sourceInputs, @{source_input_files('', @paX1, @paX2, @paXs, @paBamX)};
	$rawReads = join(";", @sourceInputs);
	die "tmpPath empty: $tmpPath" if ($tmpPath eq "");
	$tmpPath.="/rawRds/";
	my $unzipcmd = "";
	#$unzipcmd .= "set -e\n"; #sleep $WT\n
	$unzipcmd .= "rm -rf $finishStone $trimoStone $porechStone $finDest/rawRds/\nmkdir -p $finDest/rawRds/;\n";
	my $unzipcmdTMP = "rm -r -f $tmpPath;\nmkdir -p $tmpPath;\n";
	#make sure input is unzipped <- deprecated, in newer MF versions input is .gz
	my $testf1 = "";my $testf2 = "";
	my $lowEffort =-1; my $allowLinks=0;
	$allowLinks = 1 if (@pa1 == 1);
	$allowLinks = 1 if (@pa1 == 0 && @pas == 1 && !$porechopFlag);
	my $illCLip = $himipeSeqAd;
	if ($map{$curSmpl}{clip} ne ""){$illCLip = $map{$curSmpl}{clip};}
	#die "Can't find illumina trimming file: $illCLip\n" if ($useTrimomatic && !-e $illCLip);
	
	# A support-only staging directory stores its reads below rawRds/Support.
	# Inspect recursively and invalidate only stale markers; never delete valid
	# staged inputs while another cleaner may be consuming them.
	if (-d "$finDest/rawRds"
		&& !_staged_read_files_present("$finDest/rawRds", $finishStone, $trimoStone, $porechStone)){
		unlink grep { defined($_) && $_ ne '' && -e $_ } ($finishStone, $trimoStone, $porechStone);
	}

	
	
	#first conversion of bam to fastqs::
	#this is currently only working with unpaired reads!
	my $primarySingleSourceCount = scalar(@pas);
	my $supportSingleSourceCount = scalar(@paXs);
	for (my $i=0; $i<@paBam; $i++){
		my $pp = $fastp;
		my $smtBin = getProgPaths("samtools");
		my $bamFastq = outfiles_Bam("$finDest/rawRds/",basename($paBam[$i]));
		push @pas, $bamFastq;
		$unzipcmd .= "\necho \"Converting bam $i to fastq\"\n";
		$unzipcmd .= "$smtBin fastq -@ $numCore -t $pp/$paBam[$i] -0 $bamFastq;\n"; #| $pigzBin -p $numCore -c >
		$lowEffort = 0;
	}
	#and also take care of support reads in bam format
	for (my $i=0; $i<@paBamX; $i++){
		my $smtBin = getProgPaths("samtools");
		my $BamF = basename($paBamX[$i]);
		my $supportDir = "$finDest/rawRds/Support/";
		system "mkdir -p $supportDir" if ($i==0 && !-d $supportDir);
		my $bamFastq = outfiles_Bam($supportDir,$BamF);
		push @paXs, $bamFastq;
		$unzipcmd .= "echo \"Converting support bam $i to fastq\"\n";
		$unzipcmd .= "mkdir -p $supportDir;\n" if ($i==0);
		$unzipcmd .= "$smtBin fastq -@ $numCore -t $paBamX[$i] -0 $bamFastq;\n"; # | $pigzBin -p $numCore -c >
		$lowEffort = 0;
	}
	
	for (my $i=0; $i<@pa1; $i++){
		#print $pa1[$i]."\n";
		my $pp = $fastp."/";
		if ($MFconfig{filterFromSource}){
			$lowEffort =1 if ($lowEffort != 0);
			$pa1[$i] = $pp.$pa1[$i]; $pa2[$i] = $pp.$pa2[$i];
		} elsif ($porechopFlag){
#			$unzipcmd .= "$porechBin \n";
			die "porechop is not implemented for read pairs!\n";
		} elsif (0){#$useTrimomatic){
			die "trimomatic no longer supported.. function is now implemented in sdm!\n";
			#$pp = $fastp2 if ($libInfo[$i] eq $xtraRdsTech);
			#trimomatic instead of unzip
			my ($OFp1,$OFu1) = outfiles_trimall("$finDest/rawRds/",$pa1[$i]);
			my ($OFp2,$OFu2) = outfiles_trimall("$finDest/rawRds/",$pa2[$i]);
			$unzipcmd .= "java -jar $trimJar PE -threads $numCore $pp/$pa1[$i] $pp/$pa2[$i] $OFp1 $OFu1 $OFp2 $OFu2 ILLUMINACLIP:$illCLip:2:30:10\n";
			#for now: discard of singletons
			$unzipcmd .= "rm -f $OFu1 $OFu2\n";
			$pa1[$i] = $OFp1; $pa2[$i] = $OFp2;
			$lowEffort=0;
		} else {
		#old style
			my ($tmpCmd,$newF,$LEloc) = complexGunzCpMv($pp,$pa1[$i],$tmpPath,$finDest."/rawRds/",$numCore,$allowLinks);
			$unzipcmd .= $tmpCmd."\n";
			$pa1[$i] = $newF;
			$lowEffort = 0 if ($LEloc==0);
			($tmpCmd,$newF,$LEloc) = complexGunzCpMv($pp,$pa2[$i],$tmpPath,$finDest."/rawRds/",$numCore,$allowLinks);
			$unzipcmd .= $tmpCmd."\n";
			$pa2[$i] = $newF;
			
			$lowEffort = 0 if ($LEloc==0);
			#$lowEffort = 1 if ($allowLinks && $lowEffort == -1);
		}
	}
	
	
	#$unzipcmd .= "chmod +w ".join(" ",@pa1)."\n" if (!$MFconfig{filterFromSource} && @pa1>0);
	#$unzipcmd .= "chmod +w ".join(" ",@pa2)."\n" if (!$MFconfig{filterFromSource} && @pa2>0);
	#die "$unzipcmd\n";
	#for porechop this might be tons of files..
	
	for (my $i=0; $i<$primarySingleSourceCount; $i++){
		#next if (scalar(@paBam));
		my $porechopped = "$finDest/rawRds/$pas[$i]"; $porechopped .= ".gz" unless ($porechopped =~ m/\.gz$/);
		my $pp = $fastp;
		if ($i==0 && ($porechopFlag && $is3rdGen) && !$allowLinks){$unzipcmd .=  "\nrm -f $porechopped\ntouch $porechopped\n\n";}
		#print "$libInfo[$i] eq $xtraRdsTech\n";
		#$pp = $fastp2 if ($libInfo[$i] eq $xtraRdsTech);
		if ($MFconfig{filterFromSource}){
			$pas[$i] = $pp.$pas[$i]; 
		} elsif ($porechopFlag){
			#porechop is running really slow and instable, probably better to get fast5 and use modern basecaller, that will do this automatically..
			my $porechBin = getProgPaths("porechop");
			$unzipcmd .= "$porechBin -i $pp/$pas[$i] -t $numCore  --adapter_threshold 90 |gzip -c >> $porechopped\n";
			if (@pa1 > 0 ){die "no paired end reads can be given together with porechopped long reads!\n";}
		} else {
			my ($tmpCmd,$newF) = complexGunzCpMv($pp,$pas[$i],$tmpPath,$finDest."/rawRds/",$numCore,$allowLinks);
			$unzipcmd .= $tmpCmd."\n";
			$pas[$i] = $newF;
			if ($MFconfig{splitFastaInput} != 0){ #in case of input assemblies (MG-RAST.. arghh!!)
				my $sizSplitScr = getProgPaths("sizSplit_scr");
				$unzipcmd .= "\n$sizSplitScr $pas[$i] $MFconfig{splitFastaInput}\n";
			}
			$lowEffort = 1 if ($allowLinks &&  $lowEffort != 0);
		}
		if ($porechopFlag && $is3rdGen){
			$unzipcmd .= "touch $porechStone\n" if ($porechopFlag);
			$pas[$i] = ($porechopped);
		}
	}

	# Supplementary reads have their own resolved source paths and destination.
	# Keeping them separate prevents a support directory from being interpreted
	# relative to the primary input directory.
	my $supportDest = "$finDest/rawRds/Support/";
	$unzipcmd .= "mkdir -p $supportDest;\n" if (@paX1 || $supportSingleSourceCount);
	for (my $i=0; $i<@paX1; $i++) {
		my ($tmpCmd,$newF,$LEloc) = complexGunzCpMv("",$paX1[$i],$tmpPath,$supportDest,$numCore,0);
		$unzipcmd .= $tmpCmd."\n"; $paX1[$i] = $newF;
		$lowEffort = 0 if ($LEloc==0);
		($tmpCmd,$newF,$LEloc) = complexGunzCpMv("",$paX2[$i],$tmpPath,$supportDest,$numCore,0);
		$unzipcmd .= $tmpCmd."\n"; $paX2[$i] = $newF;
		$lowEffort = 0 if ($LEloc==0);
	}
	for (my $i=0; $i<$supportSingleSourceCount; $i++) {
		my ($tmpCmd,$newF,$LEloc) = complexGunzCpMv("",$paXs[$i],$tmpPath,$supportDest,$numCore,0);
		$unzipcmd .= $tmpCmd."\n"; $paXs[$i] = $newF;
		$lowEffort = 0 if ($LEloc==0);
	}
		#die "@pa1\n@pas\n";

	$lowEffort = 0 if ($lowEffort == -1); #FALLBACK option
	#die "$lowEffort\n";

	$unzipcmd .= "touch $finishStone\n";
	#if ($useTrimomatic && !$porechopFlag){$unzipcmd .= "touch $trimoStone\n" ;}
	my $jobN = "";
	#die;
	#die "unipss $unzipcmd\n";
	#print "  HH ".-s $testf2 < -s $testf1." FF \n";
	my $tmpCmd;
	#die "$unzipcmd\n$calcUnzp\n";
	if ($calcUnzp && !-e $finishStone && !$MFconfig{filterFromSource}){ #submit & check for files
		my @missingInputs = @{missing_input_files(@sourceInputs)};
		die "Missing or empty input files before unzip for $curSmpl:\n"
			.join("\n", @missingInputs)."\n" if (@missingInputs);
		if (!-e $finishStone){#|| ($useTrimomatic && !-e $trimoStone) ){
			my $lightweightLocal = commands_are_lightweight_filesystem($unzipcmd)
				&& normalise_job_dependencies($jDepe) eq '';
			if ($lightweightLocal){
				print "Executing lightweight UZ setup locally for $curSmpl\n";
				systemW $unzipcmd;
			} else {
				$jobN = "_UZ$JNUM"; 
				$unzipcmd = $unzipcmdTMP . $unzipcmd ;

				$unzipcmd = "" if (-e $finishStone && -e $trimoStone);
				my $tmpSHDD = $QSBoptHR->{tmpSpace};
				$QSBoptHR->{tmpSpace} = $HDDspace{kraken}; #set option how much tmp space is required, and reset afterwards
				($jobN, $tmpCmd) = qsubSystem($logDir."UNZP.sh",$unzipcmd,$numCore,"20G",$jobN,$jDepe,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR) ;
				$QSBoptHR->{tmpSpace} = $tmpSHDD;
			}
			#print " FDFS ";
		}
	#	die($jobN);
		#### 2 : remove human contamination
		#DB just needs to be loaded way too often.. do after sdm to have single files
		#$jobN = removeHostSeqs(\@pa1,\@pa2, \@pas,$tmpPath,$jobN);
	}
	#die "$unzipcmd\n";
	
		#check already here for mate pair support reads, deactivate fastp2
	my $mateCmd = "";
#	print "@libInfo\n";
	my $ii=0; my $mateLibraryIndex = 0; my @matePrps=();
	while($ii<@libInfo){
		#print $ii." \n";
		#next;
		unless ($libInfo[$ii] =~ m/.*mate.*/i){$ii++;next;} #remove from process
		my $sourceLabel = $libInfo[$ii];
		my $mateDir = $finDest."mateCln/lib$mateLibraryIndex/";
		my ($href,$libraryMateCmd,$mateSto) = check_matesL($pa1[$ii],$pa2[$ii],$mateDir,$doMateCln);
		$mateCmd .= $libraryMateCmd unless -e $mateSto;
		#remove this ori file from raw reads
		splice(@pa1,$ii,1);splice(@pa2,$ii,1);splice(@libInfo,$ii,1);
		push(@matePrps,{files => $href, label => $sourceLabel, index => $mateLibraryIndex});
		$mateLibraryIndex++;
	}
	if ($mateCmd ne "" && $calcUnzp){
		($jobN, $tmpCmd) = qsubSystem($logDir."MATE.sh",$mateCmd,1,"20G","_MT$JNUM",$jobN,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR) ;
	}
	
	#principally done.. now catalog and save for later
	my @inputRawStat = stat($inputRawFile);
	if (!@inputRawStat || $inputRawStat[7] == 0){
		open(my $rawInputFH, ">", $inputRawFile) or die "Cannot write $inputRawFile: $!\n";
		print {$rawInputFH} $rawReads;
		close($rawInputFH) or die "Cannot close $inputRawFile: $!\n";
	}

	
#set global var
	$map{$curSmpl}{inputFileSizeMB} = $totalInputSizeMB;
	$map{$curSmpl}{inputXFileSizeMB} = $totalXInputSizeMB;
	$map{$curSmpl}{inputFilesEmpty} = 0;

	#die "HJASD:@pa1\n@pas\n";
	my $primaryLibraries = readLibrariesFromArrays(
		sample => $curSmpl, scope => 'primary', phase => 'staged',
		technology => $seqTech, pair_technology => $seqTech,
		single_technology => $singleSeqTech, is_long => $is3rdGen, separate_roles => 1,
		r1 => \@pa1, r2 => \@pa2, single => \@pas, labels => \@libInfo,
	);
	foreach my $mateLibrary (@matePrps) {
		my $files = $mateLibrary->{files};
		my $mateIndex = $mateLibrary->{index};
		my $sourceLabel = $mateLibrary->{label} || "mate$mateIndex";
		push @{$primaryLibraries}, newReadLibrary(
			id => "$curSmpl:primary:mate:$mateIndex:pe", sample => $curSmpl,
			scope => 'primary', technology => $seqTech, is_long => $is3rdGen,
			label => "$sourceLabel.pe", phase => 'staged',
			files => {r1 => $files->{pe1}, r2 => $files->{pe2}, single => $files->{se}, bam => ''},
		);
		push @{$primaryLibraries}, newReadLibrary(
			id => "$curSmpl:primary:mate:$mateIndex:mp", sample => $curSmpl,
			scope => 'primary', technology => $seqTech, is_long => $is3rdGen,
			label => "$sourceLabel.mate", phase => 'staged',
			files => {r1 => $files->{mp1}, r2 => $files->{mp2}, single => '', bam => ''},
		);
		push @{$primaryLibraries}, newReadLibrary(
			id => "$curSmpl:primary:mate:$mateIndex:unknown", sample => $curSmpl,
			scope => 'primary', technology => $seqTech, is_long => $is3rdGen,
			label => "$sourceLabel.mate_unknown", phase => 'staged',
			files => {r1 => $files->{un1}, r2 => $files->{un2}, single => '', bam => ''},
		);
	}
	my $supportLibraries = readLibrariesFromArrays(
		sample => $curSmpl, scope => 'support', phase => 'staged',
		technology => $seqTechX, pair_technology => $seqTechX,
		single_technology => $seqTechX, is_long => $is3rdGenX, separate_roles => 1,
		r1 => \@paX1, r2 => \@paX2, single => \@paXs, labels => \@libInfoX,
	);
	%seqSet = (libraries => [@{$primaryLibraries}, @{$supportLibraries}],
			totalInputSizeMB => $totalInputSizeMB, inputXFileSizeMB => $totalXInputSizeMB,
			rawReads => $rawReads,
			mmpu => $mmpu, 
			samplReadLength => $MFconfig{defaultReadLength}, #some default value for typically short paired reads..
			samplReadLengthX => $MFconfig{defaultReadLengthX}, #for any supplementary reads (eg PacBio)
			);
	syncSeqSetLegacy(\%seqSet);
	if (exists $map{$curSmpl}{readLength} && $map{$curSmpl}{readLength} != 0){
		$seqSet{samplReadLength} = $map{$curSmpl}{readLength};
	}
	if (exists $map{$curSmpl}{readLengthX} && $map{$curSmpl}{readLengthX} != 0){
		$seqSet{samplReadLengthX} = $map{$curSmpl}{readLengthX};
	}
	
	
	my $cleanSeqSetHR = iniCleanSeqSetHR(\%seqSet);
	
	#crucial step that attaches read info to map, to make globally accessible
	sampleReadSet($curSmpl, "raw", \%seqSet);
	sampleReadSet($curSmpl, "clean", $cleanSeqSetHR);
	
	#print "$seqSet{pa1}   $seqSet{seqTech}   $seqSet{seqTechX}\n";

	#return ($jobN,\@pa1,\@pa2, \@pas, $WT, $rawReads, $mmpu,\@libInfo, $totalInputSizeMB);
	return ($jobN);
}


sub mOTU2Mapping{
	my ($tmpD,$finOutD,$smp,$Ncore,$deps) = @_;
	my $cleanSeqSetHR = sampleReadSet($curSmpl, "clean");
	my $libraries = readLibrariesByScope($cleanSeqSetHR, 'primary', 1, $curSmpl);

	
	
	#die "yes\n\n";
	my $DBdir = getProgPaths("motus2_DB",0);
	my $m2Bin = getProgPaths("motus2");#"python "."$m2Glb/$m2sub/motus";

	my $stone = $finOutD."$smp.Motu2.sto";
	my $DBstr = ""; $DBstr = "-db $DBdir " unless ($DBdir eq "");
	if (!$MFopt{DoMOTU2} ||  -e $stone){return;}
	my $pairs = libraryPairs($libraries);
	my @car1 = map { $_->{files}{r1} } @{$pairs};
	my @car2 = map { $_->{files}{r2} } @{$pairs};
	my @sar = @{libraryFiles($libraries, 'single')};
	my $inF1 = join(",",@car1); my $inF2 = join(",",@car2); my $inFS = join(",",@sar); 
	system "mkdir -p $finOutD\n" unless (-d $finOutD);
	my $cmd = "mkdir -p $tmpD\n";
	#$cmd .= "which bwa\n";
	$cmd .= "$m2Bin profile ";
	$cmd .= "-f $inF1 -r $inF2 " if (@car1 > 0 );
	$cmd .= "-s $inFS " if (@sar > 0);
	$cmd .= " -n $smp $DBstr -q -c -t $Ncore -q | gzip -c > $finOutD/$smp.motu2.tab.gz\n";
	$cmd .= "if [ -s $finOutD/$smp.motu2.tab.gz ] ; then touch  $stone ; fi\n";
	#$cmd .= "touch  $stone\n";
	my $jobN = "mOT$JNUM";
	# mOTUs uses the job working directory for sizeable alignment intermediates.
	# Scale the scheduler request with the compressed input instead of inheriting
	# the generic per-job scratch default.
	my $previousTmpSpace = $QSBoptHR->{tmpSpace};
	$QSBoptHR->{tmpSpace} =
		int(($map{$curSmpl}{inputFileSizeMB} * 4) / 1024) + 15 ."G";
	my ($jobN2,$tmpCmd) = qsubSystem($logDir."mOTU2_prof.sh",$cmd,$Ncore,"3G",$jobN,$deps,"",1,[],$QSBoptHR);
	$QSBoptHR->{tmpSpace} = $previousTmpSpace;
	$jobN  = $jobN2;
	return $jobN;
}


sub prepKraken(){
	#my ($DBdir) = @_;
	$MFglobal{krakenDBDirGlobal} = "";
	my $usesDefaultHostDB = ($MFopt{humanFilter} > 0 && $MFopt{humanFilter} < 3
		&& $MFopt{filterHostDB1} eq "");
	return "" unless ($MFopt{DoKraken} || $MFopt{DoEukGenePred} || $usesDefaultHostDB);
	if ($MFopt{DoKraken} && $MFopt{globalKraTaxkDB} eq ""){
		die "Kraken profiling requested, but no database was selected with -krakenDB\n";
	}

	# Kraken2_path_DB is the current configuration key. Keep the old name as a
	# non-required fallback for installations with a legacy local config.
	my $oriKrakDir = getProgPaths("Kraken2_path_DB",0);
	$oriKrakDir = getProgPaths("Kraken_path_DB",0) if ($oriKrakDir eq "");
	die "Kraken is enabled, but neither Kraken2_path_DB nor legacy Kraken_path_DB is configured\n"
		if ($oriKrakDir eq "");
	$oriKrakDir .= "/" unless ($oriKrakDir =~ m{/$});
	$MFglobal{krakenDBDirGlobal} = $oriKrakDir;

	my %DBname;
	if ($usesDefaultHostDB){
		$DBname{"hum1stTry"} = 1;
	}
	
	
	#die "$MFopt{humanFilter} && $MFopt{filterHostDB1}\n";
	$DBname{"minikraken_2015"} = 1 if ($MFopt{DoEukGenePred});
	$DBname{$MFopt{globalKraTaxkDB}} = 1 if ($MFopt{globalKraTaxkDB} ne "");
	return "" if (keys %DBname == 0);
	foreach my $kk (keys %DBname ){
		if (!-d "$oriKrakDir$kk"){die "can't find kraken db $oriKrakDir$kk\n";}
	}
	return "";
}


sub hostRmBase{
	my ($r1,$r2,$hostRMVer,$seqGen,$numThr,$tmpD,$krRefDB) = @_;
	my $cmd = "";
	if ($numThr eq ""){die "hostRmBase:: thread arg not correct\n";}
	return $cmd if ($r1 eq "");
	my $tmpF ="$tmpD/krak.tmp.fq";
	my $gzFlag = "";$gzFlag =  "--gzip-compressed" if ($r1 =~ m/\.gz$/);
	my $gzEnd = ""; $gzEnd = ".gz" if ($gzFlag ne "");
	if ($hostRMVer==2){
		die "kraken 1 no longer suported for host removal\n";
		my $unsplBin = getProgPaths("unsplitKrak_scr");
		my $krkBin = getProgPaths("kraken");
		$cmd .= "$krkBin --preload --threads $numThr --fastq-input $gzFlag --unclassified-out $tmpF --db $krRefDB  $r1 $r2 > /dev/null\n";
		#overwrites input files
		$r1 =~ s/\.gz$//; $r2 =~ s/\.gz$//;
		$cmd .= "$unsplBin $tmpF $r1 $r2\nrm -f $tmpF*\n";
		if ($gzFlag ne ""){
			$cmd .= "$pigzBin -f -p $numThr $r1 $r2\n";
		}
	} elsif($hostRMVer==1) { #kraken2
		my $krk2Bin = getProgPaths("kraken2");my $tmpF1 ;my $tmpF2 ; my $kr2flags= "";
		if ($r2 eq ""){
			$tmpF2=""; $tmpF1 = $tmpF;
		} else {
			$tmpF ="$tmpD/krak.tmp#.fq";  $tmpF1 = "$tmpD/krak.tmp_1.fq"; $tmpF2 = "$tmpD/krak.tmp_2.fq"; $kr2flags = "--paired ";
		}
		my $kr2QuiMod = ""; $kr2QuiMod = $MFopt{filterHostKr2QuickMode}{$seqGen} if (exists($MFopt{filterHostKr2QuickMode}{$seqGen}));
		$cmd .= "$krk2Bin --threads $numThr $gzFlag $kr2flags --unclassified-out $tmpF --db $krRefDB --output - $kr2QuiMod --confidence $MFopt{krakHostConf} $r1 $r2 \n";
		$cmd .= "$pigzBin -f -p $numThr $tmpF1 $tmpF2\n" unless ($gzFlag eq "");
		$cmd .= "rm -f $r1 $r2; \n";
		$cmd .= "mv $tmpF1$gzEnd $r1;\n ";
		$cmd .= "mv $tmpF2$gzEnd $r2;\n" if ($r2 ne "");
	} elsif($hostRMVer==3) { #hostile
		my $hostileBin = getProgPaths("hostile");
		my $hostileDB = getProgPaths("hostileDB");
		#	$MFopt{hostileIndex} = "human-t2t-hla";
		system "rm $hostileDB/$MFopt{hostileIndex}.mmi" if (-z "$hostileDB/$MFopt{hostileIndex}.mmi");
		$cmd .= "export HOSTILE_CACHE_DIR=$hostileDB\n";
		my $input = "--fastq1 $r1 "; $input .= "--fastq2 $r2 " if ($r2 ne "");
		$cmd .= "$hostileBin clean $input --index $hostileDB/$MFopt{hostileIndex} --aligner auto --output $tmpD/ --threads $numThr --airplane --force\n";
		if ($r2 ne ""){
			my $newR1 = $r1; $newR1 =~ s/.*\///; $newR1 =~ s/\.fq([\.gz])?/\.clean_1\.fastq$1/;#remove path
			$cmd .= "rm -f $r1 \nmv $tmpD/$newR1 $r1;\n";
			my $newR2 = $r2; $newR2 =~ s/.*\///; $newR2 =~ s/\.fq([\.gz])?/\.clean_2\.fastq$1/;
			$cmd .= " rm -f $r2; mv $tmpD/$newR2 $r2\n";
		} else {
			my $newR1 = $r1; $newR1 =~ s/.*\///;$newR1 =~ s/\.fq([\.gz])?/\.clean\.fastq$1/;
			$cmd .= "rm -f $r1 \nmv $tmpD/$newR1 $r1;\n";
		}
	}
	return $cmd;
}


sub removeHostSeqs($ $ $){
	my ($tmpD,$jDep,$checkIfExists) = @_;
	my $cleanSeqSetHR = sampleReadSet($curSmpl, "clean");

	return $jDep unless ($MFopt{humanFilter});
	my $libraries = ensureCleanSeqSetLibraries($cleanSeqSetHR, $curSmpl);
	
	my $outputExists=1;
	my $fileDir = "";

	if ($fileDir eq ""){
		my ($firstFile) = map { $_->{files}{r1} || $_->{files}{single} }
			grep { $_->{files}{r1} || $_->{files}{single} } @{$libraries};
		$firstFile =~ m/(.*\/)[^\/]+$/ and $fileDir = $1 if defined($firstFile);
	}
	#return $jDep if ($fileDir eq "");
	if ($fileDir eq ""){
		die "Not able to find fileDir (\"$fileDir\") for sample $curSmpl in removeHostSeqs step.. Exit\n";
	}
	#my $hostRMVer=$MFopt{humanFilter}; #0: no, 1:kraken2, 2: kraken1, 3:hostile
	#files don't need to be checked.. it's in any case just sdm files..
	$outputExists = 0 if (!-e "$fileDir/krak.stone");
	#die "$fileDir\n$outputExists && $checkIfExists\n";
	return $jDep if ($outputExists && $checkIfExists);
	

	#my ($DBdir,$DBname) = @_;
	
	my $DBdir = $MFglobal{krakenDBDirGlobal}."/";
	
	my @DBname = ("hum1stTry"); my $numThr = 4;
	$numThr = 10 if ($MFopt{humanFilter} == 3); #hostile a bit slower, use more cores
	if (@filterHostDB > 0){
		@DBname = (); my $cnt=0;
		foreach my $fhdb (@filterHostDB){
			if (!-d $fhdb){die "-filterHostKrakDB is not a dir: \"$fhdb\"\n";}
			$DBname[$cnt] = $fhdb; $DBdir=""; $cnt++;
		}
	}
	my $cmd = "\n\nmkdir -p $tmpD\n\n";  
	
	#loop around different DBs..
	
	for (my $j=0;$j<@DBname ; $j++){
		foreach my $library (@{$libraries}) {
			$cmd .= hostRmBase($library->{files}{r1},$library->{files}{r2},$MFopt{humanFilter},$library->{is_long},$numThr,$tmpD,"$DBdir$DBname[$j]")
				if ($library->{files}{r1});
			$cmd .= hostRmBase($library->{files}{single},"",$MFopt{humanFilter},$library->{is_long},$numThr,$tmpD,"$DBdir$DBname[$j]")
				if ($library->{files}{single});
		}
	}
	$cmd .= "\n\n";
	$cmd .= "touch $fileDir/krak.stone\n";
	#die "$cmd\n";
	my $jobN = ""; my $tmpCmd="";
	unless (-e "$fileDir/krak.stone"){
		my $tmpSHDD = $QSBoptHR->{tmpSpace};
		my $reqSpace = $HDDspace{kraken};
		if ($map{$curSmpl}{inputFileSizeMB} + $map{$curSmpl}{inputXFileSizeMB} > 10000){
			#convert input to a) Gb and b) *6 to account for extra space needed by kraken
			$reqSpace = int($map{$curSmpl}{inputFileSizeMB}*6/1024)+30  ."G";
		} 
		#die "krakHS : $reqSpace $map{$curSmpl}{inputFileSizeMB} \n\n";
		$QSBoptHR->{tmpSpace} = $reqSpace; #set option how much tmp space is required, and reset afterwards
		#### 1 : UNZIP
		$jobN = "_KR$JNUM";
		($jobN, $tmpCmd) = qsubSystem($logDir."KrakHS.sh",$cmd,$numThr,"16G",$jobN,$jDep.";$krakDeps","",1,$QSBoptHR->{General_Hosts},$QSBoptHR) ;
			$QSBoptHR->{tmpSpace} = $tmpSHDD;
	}
	return $jobN;
}

sub genoSize(){
	my ($oD,$jdep) = @_;
	my $cleanSeqSetHR = sampleReadSet($curSmpl, "clean");

	my $pairs = libraryPairs(readLibrariesByScope($cleanSeqSetHR, 'primary', 1, $curSmpl));
	my $cmd = "mkdir -p $oD\n";
	my $microCensBin = getProgPaths("microCens");
	for (my $i=0;$i<@{$pairs};$i++){
		$cmd .= "$microCensBin -t 2 $pairs->[$i]{files}{r1},$pairs->[$i]{files}{r2} $oD/MC.$i.result\n";
	}
#	die $cmd;
	my $jobName = "_GS$JNUM"; my $tmpCmd;
	
	($jobName,$tmpCmd) = qsubSystem($logDir."MicroCens.sh",$cmd,2,"40G",$jobName,$jdep,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR);
	return $jobName;
}
sub krakenTaxEst(){
	my ($outD, $tmpD,$name,$jobd) = @_;

	my $cleanSeqSetHR = sampleReadSet($curSmpl, "clean");
	my $libraries = readLibrariesByScope($cleanSeqSetHR, 'primary', 1, $curSmpl);
	my $pairs = libraryPairs($libraries);
	my @pas = @{libraryFiles($libraries, 'single')};
	#$outD.= "$MFopt{globalKraTaxkDB}/";
	my $krakStone = "$outD/krakDone.sto";
	return $jobd if (-d $outD && -e $krakStone);
	my $numCore = $MFopt{krakenCores};
	my @thrs = (0.01,0.02,0.04,0.06,0.1,0.2,0.3);
	my $cmd = "mkdir -p $outD\nrm -rf $tmpD\nmkdir -p $tmpD\n";
	my $it =0;
	my $curDB = "$MFglobal{krakenDBDirGlobal}/$MFopt{globalKraTaxkDB}";

	my $krkBin = "";my $krak1 = 1;
	$krkBin = getProgPaths("kraken") if ($krak1);

	#die $curDB."\n";
	#paired read tax assign
	for (my $i=0;$i<@{$pairs};$i++){
		my $r1 = $pairs->[$i]{files}{r1}; my $r2 = $pairs->[$i]{files}{r2};
		$cmd .= "$krkBin --paired --preload --threads $numCore --fastq-input  --db $curDB  $r1 $r2 >$tmpD/rawKrak.$it.out\n";
		for (my $j=0;$j< @thrs;$j++){
			$cmd .= "$krkBin-filter --db $curDB  --threshold $thrs[$j] $tmpD/rawKrak.$it.out | $krkBin-translate --mpa-format --db $curDB > $tmpD/krak_$thrs[$j]"."_$it.out\n";
		}
		$it++;
	}
	$cmd .= "\n\n";
	#single read tax assign
	for (my $i=0;$i<@pas;$i++){
		my $rs = $pas[$i]; 
		$cmd .= "$krkBin --preload --threads $numCore --fastq-input  --db $curDB  $rs >$tmpD/rawKrak.$it.out\n";
		#$cmd .= " | tee ";#$krkBin-filter --db $curDB  --threshold 0.01 | $krkBin-translate --mpa-format --db $curDB > $tmpD/krak$it.out\n";
		for (my $j=0;$j< @thrs;$j++){
			$cmd .= "$krkBin-filter --db $curDB  --threshold $thrs[$j] $tmpD/rawKrak.$it.out | $krkBin-translate --mpa-format --db $curDB > $tmpD/krak_$thrs[$j]"."_$it.out\n";
		}
		#overwrites input files
		$it++;
	}
	
	#TODO: 1: make table; 2: copy to outD
	$cmd .= "\n\n";
	my $krakCnts1 = getProgPaths("krakCnts_scr");
	for (my $j=0;$j< @thrs;$j++){
		$cmd .= "cat $tmpD/krak_$thrs[$j]"."_*.out > $tmpD/allkrak$thrs[$j].out\n";
		$cmd .= "$krakCnts1 $tmpD/allkrak$thrs[$j].out $outD/krak.$thrs[$j].cnt.tax\n";
		$cmd .= "[ -s $outD/krak.$thrs[$j].cnt.tax ] || exit 4\n";
	}

	$cmd .= "touch $krakStone\n";
	$cmd .= "rm -fr $tmpD";

	my $jobName = "_KT$JNUM"; my $tmpCmd;
	if (!-d $outD || !-e $krakStone){
		my $previousTmpSpace = $QSBoptHR->{tmpSpace};
		my $requestedTmpSpace = $HDDspace{kraken};
		if ($map{$curSmpl}{inputFileSizeMB} > 10000) {
			# raw classifications plus the translated threshold series coexist in
			# $tmpD until the final count tables have been validated.
			$requestedTmpSpace =
				int(($map{$curSmpl}{inputFileSizeMB} * 6) / 1024) + 30 ."G";
		}
		$QSBoptHR->{tmpSpace} = $requestedTmpSpace;
		($jobName,$tmpCmd) = qsubSystem($logDir."KrkTax.sh",$cmd,$numCore,"20G",$jobName,$jobd.";$krakDeps","",1,$QSBoptHR->{General_Hosts},$QSBoptHR);
		$QSBoptHR->{tmpSpace} = $previousTmpSpace;
	}
	return $jobName;
}

sub check_map_done{
	my ($doCram, $finalD, $baseN) = @_;
	# Only canonical outputs count as complete. Older scratch-only results are
	# ignored so partial outputs cannot cross MATAFILER version boundaries.
	if ($doCram){
		return 1 if (-s "$finalD/$baseN-smd.cram" && -e "$finalD/$baseN-smd.cram.sto");
	} else {
		return 1 if (-s "$finalD/$baseN-smd.bam");
	}
	return 0;
}
sub check_depth_done{
	my (undef, $finalD, $baseN) = @_;
	return 1 if (fileGZs("$finalD/$baseN-smd.bam.coverage"));
	return 0;
}

#create string for mapper to register libraries
sub getRgStr{ 
	my ($smpl,$libsOri,$libsOriX,$usePairs,$mapper,$readTechnology) = @_;
	my $rgStr ="noReg";
	my $platform = 'ILLUMINA';
	$platform = 'PACBIO' if (($readTechnology || '') eq 'PB');
	$platform = 'ONT' if (($readTechnology || '') eq 'ONT');
	if ($mapper > 1 || $mapper == -2){ #bwa/minimap2 have same format..
		$rgStr = '\'@RG\\tID:$smpl\\tSM:'.$smpl.'\\tPL:'.$platform;
		$rgStr .= '\\tLB:'.$libsOri.'\'';
	}
	if ($mapper==1 || $mapper ==5){ #bowtie2 & strobealign
		my $sep=" "; $sep = "=" if ($mapper ==5);
		$rgStr = "--rg-id${sep}$smpl --rg${sep}SM:$smpl --rg${sep}PL:$platform "; #PU:lib1
		if ($usePairs){
			$rgStr .= "--rg${sep}LB:$libsOri ";
			$rgStr .= " -X $MFconfig{mateInsertLength} " if ($libsOri =~ m/mate/ && $mapper==1);
		} else {
			$rgStr .= "--rg${sep}LB:$libsOriX ";
			$rgStr .= " -X $MFconfig{mateInsertLength} " if ($libsOriX =~ m/mate/ && $mapper==1);
		}
	}
	return $rgStr;
}



sub getAlgnCmdBase{
	my ($MapperProg,$NcoreL, $readTec,$mapModeTogether,$unaligned) = @_;
	my $algCmdBase = "";
	if ($MapperProg==1){ #bowtie2
		my $bwt2Bin = getProgPaths("bwt2");
		$algCmdBase = "$bwt2Bin  --end-to-end -p $NcoreL  "; #--no-unal
		if ($mapModeTogether == 2 || $mapModeTogether == -1){$algCmdBase .= "-a ";}
	} elsif ($MapperProg==2){
		my $bwaBin = getProgPaths("bwa");
		$algCmdBase = "$bwaBin mem -t $NcoreL ";
	} elsif ($MapperProg==3){#minimap2 
		my $mini2Bin = getProgPaths("minimap2");
		$algCmdBase = "$mini2Bin -2 -t $NcoreL --secondary=no ";
		if ($readTec eq "ONT"){
			$algCmdBase .= " -x map-ont" ; #use nanopore optimzed for now..
		} elsif ($readTec eq "PB"){
			$algCmdBase .= " -x map-pb" ; #use nanopore optimzed for now..
		} else {
			print"Warning: Minimap2 used for short reads; not recommended\n";
			$algCmdBase .= " -x sr" ;
		}
		#$algCmdBase .= " --sam-hit-only " unless ($unaligned); #deactivated as important for counting
	} elsif ($MapperProg==4){ #kma 
		my $kmaBin = getProgPaths("kma");
		my $consID = 0.95;my $minPhred = 15; my $minMapQ = 20; my $minQueryCov = 0.2; #-bc $minPhred -tmp $nodeTmp/${baseN}.kmatmp/
		$algCmdBase = "$kmaBin -nc -na -nf -sam 4 -apm p -mrc $minQueryCov -mq $minMapQ -bcd 1 -ID $consID -ref_fsa -t $NcoreL   ";
		if ($readTec eq "ONT"){$algCmdBase .= " -bcNano -ont " ;
		} elsif ($readTec eq "PB"){$algCmdBase .= " -mint3  " ;
		} else {$algCmdBase .= " -mint2 " ;
		}
	} elsif ($MapperProg==5){ #strobealign
		my $stroBin = getProgPaths("strobealign");
		$algCmdBase = "$stroBin -t $NcoreL --no-progress ";
		#$algCmdBase .= " -U " unless ($unaligned);

	} else {
		die "getAlgnCmdBase:: unknown mapper \"$MapperProg\"\nAborting..\n";
	}
	return $algCmdBase;
}



sub announce_MF4{
	
	
#	print "888b     d888        d8888 88888888888     d8888 8888888888 8888888 888      8888888888 8888888b.  \n";#
#	print "8888b   d8888       d88888     888        d88888 888          888   888      888        888   Y88b \n";
#	print "88888b.d88888      d88P888     888       d88P888 888          888   888      888        888    888 \n";
#	print "888Y88888P888     d88P 888     888      d88P 888 8888888      888   888      8888888    888   d88P \n";
#	print "888 Y888P 888    d88P  888     888     d88P  888 888          888   888      888        8888888P\"  \n";
#	print "888  Y8P  888   d88P   888     888    d88P   888 888          888   888      888        888 T88b   \n";
#	print "888   \"   888  d8888888888     888   d8888888888 888          888   888      888        888  T88b  \n";
#	print "888       888 d88P     888     888  d88P     888 888        8888888 88888888 8888888888 888   T88b \n";

	print "/------------------------------------------------------------------------\\\n";
	print "|  _______ _______ _______ _______ _______ _____        _______  ______  |\n";
	print "|  |  |  | |_____|    |    |_____| |______   |   |      |______ |_____/  |\n";
	print "|  |  |  | |     |    |    |     | |       __|__ |_____ |______ |    \\_  |\n";
	print "|                                                                        |\n";
	print "\\------------------------------------------------------------------------/\n";
#	print "/--------------------------------------------\\\n";
#	print "|  ███╗   ███╗ ██████╗    ████████╗██╗  ██╗  |\n";
#	print "|  ████╗ ████║██╔════╝    ╚══██╔══╝██║ ██╔╝  |\n";
#	print "|  ██╔████╔██║██║  ███╗█████╗██║   █████╔╝   |\n";
#	print "|  ██║╚██╔╝██║██║   ██║╚════╝██║   ██╔═██╗   |\n";
#	print "|  ██║ ╚═╝ ██║╚██████╔╝      ██║   ██║  ██╗  |\n";
#	print "|  ╚═╝     ╚═╝ ╚═════╝       ╚═╝   ╚═╝  ╚═╝  |\n";
#	print "\\--------------------------------------------/\n";

	print "This is MATAFILER4 v$MATFILER_ver\n";
}




sub getMapProgNm{
	my ($MapperProg) = @_;
	my $mapProgNm = "undefined mapper";
	if ($MapperProg==1){$mapProgNm = "bowtie2";
	}elsif ($MapperProg==2){$mapProgNm = "bwa";
	}elsif ($MapperProg==3){$mapProgNm = "minimap2";
	}elsif ($MapperProg==4){$mapProgNm = "kma";
	}elsif ($MapperProg==5){$mapProgNm = "strobealign";
	} else {print "coult not recognize mapper $MapperProg!!\n";
	}
	return $mapProgNm
}

#			my %postTreat = (MapperProg=>$MapperProg, readTec=>$readTec, NcoreL => $NcoreL,nodeTmp =>$nodeTmp,tmpOut21=>$tmpOut21, 
#										subBamsAR => \@subBams,unaligned =>$unaligned,baseN => $baseN,xtraSamSteps1 => $xtraSamSteps1,
#										decoyModeActive => $decoyModeActive, map2ndTogether => $map2ndTogether,regsAR => \@regs, 
#										doCram => $doCram, finalDSar => \@finalDS, outNmsAR => \@outNms, mappDirAR => \@mappDir, 
#										reg_lcsAR => \@reg_lcs)
sub alignPostTreat{
	my ($postTreatHR, $i, $kk) = @_;
	my %postTreat = %{$postTreatHR};
	my $bamHdFilt_scr = getProgPaths("bamHdFilt_scr");
	
	my $algCmd = "";
	my ($MapperProg, $readTec, $NcoreL,$nodeTmp,$tmpOut21, $subBamsAR,$unaligned,$baseN) = ($postTreat{MapperProg}, $postTreat{readTec}, $postTreat{NcoreL},$postTreat{nodeTmp}, $postTreat{tmpOut21}, $postTreat{subBamsAR},$postTreat{unaligned},$postTreat{baseN});
	my $xtraSamSteps1 = $postTreat{xtraSamSteps1}; 
	
	my $bamfilter = getProgPaths("bamFilter_scr");
		#my $filterStep="";
	if ( ($MapperProg==3 || $MapperProg==5 ) && ($readTec eq "ONT" || $readTec eq "PB")){ #low id long reads...
		if ($readTec eq "PB" ){
			$algCmd .= " | $bamfilter $MFopt{bamfilterPB}  ";
		} else {
			$algCmd .= " | $bamfilter $MFopt{bamfilterONT}  "; #defaults to ONT (safer parameters)
		}
	} else { #illumina parameters
		my $ill_filter = $postTreat{strictHybridCoverage}
			? $MFopt{bamfilterHybridIll} : $MFopt{bamfilterIll};
		$algCmd .= " | $bamfilter $ill_filter ";
	}


	
	my $iTO = "$tmpOut21.$i.$kk";
	if ($unaligned ne "" ){ #basic sort of unmapped reads via samtools (general purpose step)
		$algCmd .= " | $smtBin view -b1 -@ $NcoreL > $iTO.t\n";
		$algCmd .= "$smtBin view -u -h $iTO.t | $xtraSamSteps1 $smtBin view -b1 -@ $NcoreL -F 4 -  > $iTO\n";
		#sort out unaligned reads
		$algCmd .= "$smtBin view -u -h -@ $NcoreL -f 4 $iTO.t | $smtBin fastq -1 $nodeTmp/$baseN.$i.1.fq.gz -2 $nodeTmp/$baseN.$i.2.fq.gz -s $nodeTmp/$baseN.$i.s.fq.gz - \n";
		#and copy them already to final destination.. no reason to keep them around..
		$algCmd .= "cat $nodeTmp/$baseN.$i.1.fq.gz >> $unaligned/unal.1.fq.gz;\ncat $nodeTmp/$baseN.$i.2.fq.gz >> $unaligned/unal.2.fq.gz;\n cat $nodeTmp/$baseN.$i.s.fq.gz >> $unaligned/unal.fq.gz;\n";
		#and remove all the temp files..
		$algCmd .= "rm -f $nodeTmp/$baseN.$i.*fq.gz $iTO.t\n";
	} else {
		$algCmd .= " | $xtraSamSteps1  $smtBin view -b1 -@ $NcoreL -F 4 - > $iTO\n";
	}
	
	if ($postTreat{decoyModeActive} || $postTreat{map2ndTogether}>0){#remove unnecessary reads & filter specific reads for each ref genome into separate bam
		$algCmd .= "$smtBin index $iTO\n";
		for (my $k=0;$k<@{$postTreat{regsAR}};$k++){
			#		print "$k\t$finalDS[$k]\n";
			if(check_map_done(${$postTreat{doCram}}, ${$postTreat{finalDSar}}[$k], ${$postTreat{outNmsAR}}[$k])){$$subBamsAR[$k]="";next;}
			$algCmd .= "\n\nset +e \n" if ($k==0);
			$algCmd .= "#  %%%%%%%%%%%%%%%% $k %%%%%%%%%%%%%%%% \n";
			$algCmd .= "$bamHdFilt_scr $iTO ${$postTreat{reg_lcsAR}}[$k] 0 > $iTO.decoy.sam.$k\n";
			$algCmd .= "sleep 0.05\n";
			$algCmd.= "  $smtBin view -@ $NcoreL $iTO ${$postTreat{regsAR}}[$k] >> $iTO.decoy.sam.$k\n";
			$algCmd.= "  $smtBin view -b -h -@ $NcoreL $iTO.decoy.sam.$k > $iTO.decoy.$k\nrm -f $iTO.decoy.sam.$k \n";
			$$subBamsAR[$k] .= " $iTO.decoy.$k";  
		}
		
		$algCmd.= "rm -f $iTO\n";# mv $tmpOut21.decoy.$i $tmpOut21.$i\n";
		$algCmd .= "\n\nset -e \n";

	} else {
		$$subBamsAR[0] .= " $iTO"; 
	}
	
	#die "$algCmd\n";
	return  ($algCmd,$subBamsAR);
}


sub mapReadsToRef{
	my ($dirsHr, $jDepe) = @_;
	#my ($par1,$par2,$parS,$liar,$rear) = getRawSeqsAssmGrp($AsgHR,$ASG,$supportRds);
	#wrong: needs to be reads of the sample, not assembly group!
	my $outName = $dirsHr->{smplName};my $ASG = $dirsHr->{assGrp};
	my $is2ndMap = $dirsHr->{is2ndMap};
	my $doCram = $dirsHr->{cramAlig};
	my $immediateSubm = $dirsHr->{submNow};
	my $finalUnaligned = $dirsHr->{unalDir}; #final unaligned-read directory; "" disables retention
	my $Ncore = $dirsHr->{mapCores};
	my $supportRds = $dirsHr->{mapSupport};
	my $supTag = ""; if ($supportRds){$supTag = ".sup";}
	#die "$supTag $supportRds $dirsHr->{mapSupport} \n";

	my $REF = $dirsHr->{sbj}; #target to map onto, can by ","-spearated list
	#print "$outName\n";
	my $libraries = ref($dirsHr->{libraries}) eq 'ARRAY'
		? $dirsHr->{libraries}
		: getRawLibrariesAssmGrp(\%AsGrps,$ASG,$supportRds,$outName);
	my $pairs = libraryPairs($libraries);
	my @pa1 = map { $_->{files}{r1} } @{$pairs};
	my @pa2 = map { $_->{files}{r2} } @{$pairs};
	my @paS = @{libraryFiles($libraries, 'single')};
	my @singleLibraries = grep { $_->{files}{single} } @{$libraries};
	my @mappingLibraries = (@{$pairs}, @singleLibraries);
	my @libsOri = map { $_->{label} || $_->{id} || 'library' } @mappingLibraries;
	my $recordReadTechnology = libraryTechnology($libraries, "mapping sample $outName", 1);
	#simple rule for mapper program: for now set to bowtie2

	#die "$REF\n";
			# ($dirsHr,$outName, $is2ndMap, $par1,$par2,$Ncore,#hm,m,x,x,x,x,x
	#	$REF,$jDepe, $unaligned,$immediateSubm,#m,x,x,x,m
	#	$doCram,$smpl,$libAR) = @_;#x,x,x
	#get mapper progs..
	#my $bwaBin = "";
	# if ($MapperProg==2);

		#get essential dirs...
	my $mappDirPre = ${$dirsHr}{glbMapDir};	my $nodeTmp = ${$dirsHr}{nodeTmp}; #m,s
	$nodeTmp.="_map${supTag}/";
	my $qdir = $logDir; $qdir = ${$dirsHr}{qsubDir} if (exists( ${$dirsHr}{qsubDir} ));
	my $declaredReadTechnology = ${$dirsHr}{readTec} || '';
	die "Mapping technology '$declaredReadTechnology' disagrees with library records '$recordReadTechnology' for $outName\n"
		if ($declaredReadTechnology ne '' && $recordReadTechnology ne ''
			&& $declaredReadTechnology ne $recordReadTechnology);
	my $readTec = $recordReadTechnology || $declaredReadTechnology;
	my $tmpOut = ${$dirsHr}{glbTmp};	my $finalD = ${$dirsHr}{outDir}; #node-local work, final mapping
	my $unaligned = $finalUnaligned eq "" ? "" : $tmpOut."/unaligned/";
	my $mapperProgLoc = decideMapper($MFopt{MapperProg},$readTec);
	my $mappingInputSizeMB = $supportRds
		? ($map{$curSmpl}{inputXFileSizeMB} || 0)
		: ($map{$curSmpl}{inputFileSizeMB} || 0);
	#die "$finalD\n";
	my @finalDS = split /,/,$finalD;
	my @mappDir = split /,/,$mappDirPre;

	#my $bwtIdx = $REF."$MFcontstants{bwt2IdxFileSuffix}";
	my @bwtIdxs;
	#already exists
	#print "IS${unaligned}SI\n";
	my $baseN = "$outName$supTag";#$RNAME."_".$QNAME;
	
	my @outNms = split /,/,$outName;
	if ($outName ne "" ){
		$baseN = $outNms[0].$supTag;
	} else {
		@outNms = ($baseN);
	}
	if ($supportRds){for (my $ii=0;$ii<@outNms; $ii++){$outNms[$ii] = $outNms[$ii].$supTag;}}
	#die "$baseN @outNms\n";
	my @tmpOut22 = ($tmpOut."/$baseN.iniAlignment.bam");
	my @tmpOutxtra = ($nodeTmp."/$baseN.iniAlignment.xtra");
	#global value overwrites local value
	if ($doCram){$doCram = $MFopt{doBam2Cram};}
	#calculate total input size (to get handle on req disk space
	

	
	my %params = (
		mappingCommand => "",
		mappingDependencies => $jDepe,
		mappingStarted => 1,
	);
	my $bamFresh = 0; #is the bam newly being created?
	my $decoyModeActive=0; #decoy mapping
	my $map2ndTogether = $MFopt{mapModeTogether}; #map competetively among all reference genomes provided might change mapping result, if other genome set is used)
	$decoyModeActive=1 if ( $MFopt{DoMapModeDecoy} && exists($make2ndMapDecoy{Lib}) && -e $make2ndMapDecoy{Lib});
	my $isSorted = 0;		$isSorted=1 if (@pa1==1 && $MFopt{DoMapModeDecoy} && $decoyModeActive);
#	die $REF;
	if (!$decoyModeActive){
		@bwtIdxs = split /,/,$REF;
		for (my $kk=0;$kk<@bwtIdxs;$kk++){
			$bwtIdxs[$kk] .= $MFcontstants{bwt2IdxFileSuffix};
		}
	}
	my $referencePreparationCommand = "";
	if ($REF =~ /\.gz$/) {
		my $compressedReference = $REF;
		my $stagedReference = "$nodeTmp/".basename($compressedReference);
		$stagedReference =~ s/\.gz$//;
		$referencePreparationCommand .= "mkdir -p $nodeTmp\n";
		$referencePreparationCommand .= "$pigzBin -dc $compressedReference > $stagedReference\n";
		$referencePreparationCommand .= "test -s $stagedReference\n";
		if ($mapperProgLoc == 2 || $mapperProgLoc == 4) {
			$referencePreparationCommand .= "for f in $compressedReference.*; do [ -e \"\$f\" ] || continue; "
				."suffix=\${f#$compressedReference}; cp \"\$f\" \"$stagedReference\$suffix\"; done\n";
		}
		$REF = $stagedReference;
	}
	$params{mappingReference} = $REF;
	$params{referencePreparationCommand} = $referencePreparationCommand;
	$params{mappingInputSizeMB} = $mappingInputSizeMB;
	#die "$isSorted\n";
	my $anyUsedPairs= 1; $anyUsedPairs = 0 if (scalar @pa1 == 0);
	$params{sortedbam}=$isSorted; $params{bamIsNew} = $bamFresh; $params{is2ndMap} = $is2ndMap;
	$params{immediateSubm} =  $immediateSubm; $params{usePairs} = $anyUsedPairs;

	my $outputExistsNEx = 0;
	for (my $k=0;$k<@outNms;$k++){
		$tmpOut22[$k] = $tmpOut."/$outNms[$k].iniAlignment.bam"; 
		$tmpOutxtra[$k] = $tmpOut."/$outNms[$k].iniAlignment.xtra"; 
		#print "$tmpOut22[$k]\n";
	}
	for (my $k=0;$k<@outNms;$k++){
		if ($is2ndMap && $MFopt{MapRewrite2nd}){$outputExistsNEx++; system "rm -fr $tmpOut22[$k] $finalDS[$k]/$outNms[$k]-smd*"; next;}
		my $outstat = check_map_done($doCram, $finalDS[$k], $outNms[$k]);
		next if ($outstat);#-e "$finalDS[$k]/$outNms[$k]-smd.bam.coverage.gz" );#|| -e "$mappDir[$k]/$outNms[$k]-smd.bam.coverage.gz");
		#print "-e $finalDS[$k]/$outNms[$k]-smd.bam.coverage.gz || -e $mappDir/$outNms[$k]-smd.bam.coverage.gz";
		my @temporaryOutputStat = stat($tmpOut22[$k]);
		if (@temporaryOutputStat && $temporaryOutputStat[7] < 100){
			unlink $tmpOut22[$k];
			@temporaryOutputStat = stat($tmpOut22[$k]);
		}
		if (!@temporaryOutputStat){ $outputExistsNEx++; }
	}
	#die "@outNms \n$outputExistsNEx\n";
	if ($outputExistsNEx == 0){
		#die;
		return ("","",\%params);
	} 
	
	my $NcoreL = int($Ncore);#correction for threads
	#my $nxtCRAM = "$tmpOut/$baseN-smd.cram";
	my $tmpOut21 = $nodeTmp."/$baseN.iniAlignment.bam";
	#my $sortTMP = $nodeTmp."/$baseN.srt";supportRds
	#print "$nxtBAM\n";
	#die("too far\n");
	
	my $retCmds="";my $tmpCmd=""; my $xtraSamSteps1="";
	#my $nxtCRAM = "$mappDir/$baseN-smd.cram";

	# All working paths are created by the combined mapping job on its assigned
	# node.  Submission-time mkdir calls are incorrect for node-local paths.
	#move ref DB & unzip
	#system("rm -r $tmpOut\n mkdir -p $tmpOut\n cp $REF $tmpOut");
	#$REF = $tmpOut.basename($REF);
	my $jobN = "";
	my $tmpUna = $tmpOut."/unalTMP/";
	#Error: No EOF block on /g/scb/bork/hildebra/SNP/MeHiAss/MH0411//tmp//mapping/Alignment.bam, possibly truncated file.
	my $bashN = "";	if ($is2ndMap){$bashN = "$outName"; $bashN =~ s/,/./g; if (length($bashN)>50){$bashN=substr($bashN,0,40)."_etc";}}
	if (0&& -e $qdir.$bashN."bwtMap2.sh.etxt"){open I,"<$qdir/".$bashN."bwtMap2.sh.etxt" or die "Can't open old bowtie_2 logfile $qdir\n"; my $str = join("", <I>); close I;
		if ($str =~ /Error: No EOF block on (.*), possibly truncated file\./){system("rm -f $1 ".join (" ",@tmpOut22) );}	
		close I;	
	}
	#depending on setup, mappDir == tmpOut
	my $algCmd = "";
	$algCmd .= "\nrm -rf $tmpOut $nodeTmp\nmkdir -p $tmpOut\nmkdir -p $nodeTmp\n"; 
	$algCmd .= $referencePreparationCommand;
	#die "$algCmd\n@mappDir\n";
	#my $algCmd = "rm -rf ".join(" ",split(/,/,$mappDir))."\nmkdir -p ".join(" ",split(/,/,$mappDir))." $tmpOut $nodeTmp\n"; 
	$algCmd .= "mkdir -p $qdir/mapStats/\n";
	
	my @regs;#subset of DB seqs to filter for 
	my @reg_lcs;
	$REF =~ m/([^\/]+)$/; my $REFnm = $1; $REFnm =~ s/\.f.*a$//; 
	#die "$REF\n$REFnm\n";
	if ($decoyModeActive || $map2ndTogether){
		#first create a new ref DB, including the targets and the assemblies from this dir
		my $bwtIdx = "$nodeTmp/$REFnm.decoyDB.fna";
		if ($map2ndTogether){
			$bwtIdx = "$map2ndTogRefDB{DB}";
		} else {
			my $decoyDBscr = getProgPaths("decoyDB_scr"); 
			$algCmd.= "\n$decoyDBscr $REF $make2ndMapDecoy{Lib} $bwtIdx $NcoreL $outName $finalD\n";
		}
		#die "$algCmd\n";
		$bwtIdx .= $MFcontstants{bwt2IdxFileSuffix};
		
		$xtraSamSteps1 = "$smtBin sort -@ $NcoreL -T $nodeTmp./$baseN.srt - |";
		@regs = @{$make2ndMapDecoy{regions}};
		@reg_lcs =  @{$make2ndMapDecoy{region_lcs}};
		@bwtIdxs = ($bwtIdx);
	} 
	my $bwt2DBsuf = "bt2"; $bwt2DBsuf = "bt2l" if ($MFopt{largeMapperDB}); 
	$algCmd .= "if [ ! -e $bwtIdxs[0].1.$bwt2DBsuf ] ;then\n	echo \"Could not find assembly bowtie2 index: $bwtIdxs[0].1.$bwt2DBsuf\";\n	exit 23;\nfi\n" if ($mapperProgLoc == 1); #needs to exit on error

	
	my $cntAli=0; my $totlRefs = scalar(@bwtIdxs);
	my $numLib = scalar @pa1 + scalar @paS; 
	#if there's too many refgenomes, copy reads onto tmp dir
	if ($totlRefs > 5){
		$algCmd .= "\n\n#copying read files to local tmp\n";
		for (my $ii=0;$ii<@pa1;$ii++){
			$algCmd .= "cp $pa1[$ii] $nodeTmp\n";$pa1[$ii]=~ m/\/([^\/]+$)/;$pa1[$ii] = $nodeTmp."/$1";
		}
		for (my $ii=0;$ii<@pa2;$ii++){
			$algCmd .= "cp $pa2[$ii] $nodeTmp\n";$pa2[$ii]=~ m/\/([^\/]+$)/;$pa2[$ii] = $nodeTmp."/$1";
		}
		for (my $ii=0;$ii<@paS;$ii++){
			$algCmd .= "cp $paS[$ii] $nodeTmp\n";$paS[$ii]=~ m/\/([^\/]+$)/;$paS[$ii] = $nodeTmp."/$1";
		}
	}
	
	my $actualMappings=0;
#	if ($numLib==0){
#		$numLib = scalar @paS; 
#	}
	#die "@bwtIdxs\n"; 
	#organize mapping against 1 or more refs
	for (my $kk=0;$kk<@bwtIdxs;$kk++){  #iterator over different genomes
		#die  "$finalDS[$kk]/$outNms[$kk]-smd.bam\nXXX\n$mappDir[$kk]\n";
		if ($decoyModeActive || $map2ndTogether>0){
			my $allDone=1;
			for (my $k=0;$k<@outNms;$k++){
				#print "$finalDS[$k] $outNms[$k] $mappDir[$k]\n";
				$allDone =0 unless(check_map_done($doCram, $finalDS[$k], $outNms[$k]));
			}
			#die;
			last if ($allDone);
		} elsif (-e $tmpOut22[$kk] || check_map_done($doCram, $finalDS[$kk], $outNms[$kk])){ #this check is for non-decoy mode
			#(-e "$finalDS[$kk]/$outNms[$kk]-smd.bam" && -e "$finalDS[$kk]/$outNms[$kk]-smd.bam.coverage.gz") ){
			next;
		} #|| -e "$mappDir[$kk]/$outNms[$kk]-smd.bam.coverage.gz"
		$cntAli++;
		$bamFresh=1;
		#print "ali\n";
		$algCmd .= "#--------------- $kk ---------------\nmkdir -p $mappDir[$kk];\n";
		$algCmd .= "rm -rf $unaligned $tmpUna;\nmkdir -p $unaligned $tmpUna;\n" if ($unaligned ne "");
		
		my $algCmdBase = getAlgnCmdBase($mapperProgLoc,$NcoreL, $readTec,$MFopt{mapModeTogether},$unaligned);
#die "algCmdBase::$algCmdBase\n";
		my @subBams; 
		my @accR1=(); my @accR2=(); my @accRS=();
		for (my $i=0; $i< $numLib; $i++){
			$actualMappings++;
			my $usePairs=1;
			if ($i >= scalar @pa1){$usePairs=0;}
			$anyUsedPairs =1 if ($usePairs);
			my $iS = $i - scalar @pa1;
			
			#test if more reads can be accummulated in the next round...
			#print "@libsOri  : $i $numLib ". scalar @pa1 ."\n";
			if ($mapperProgLoc==1 ){
				#accumulate reads..
				if ($usePairs){ push(@accR1,$pa1[$i]); push(@accR2, $pa2[$i]); 
				} else { push(@accRS,$paS[$iS]); }
				if ( ($i+1) < $numLib && ((($i+1) >= scalar @pa1 && $usePairs==0) || (($i+1) < scalar(@pa1) && $usePairs==1))
						&& $libsOri[$i] eq $libsOri[$i+1]){ #same reads, same lib in next round, all set!
					next;
				}
			} else { #other mappers can't use multiple input files..
				if ($usePairs){ @accR1 = ($pa1[$i]); @accR2 = ($pa2[$i]); 
				} else { @accRS  = ($paS[$iS]); }
			}
			
			#$pa1[$i] =~ m/\/([^\/]+)\.f.*q$/;
			#my $rgID = "$outName";
			my $rgStr = getRgStr($outName,$libsOri[$i],$libsOri[$i],$usePairs,$mapperProgLoc,$readTec);
			#die "$rgStr\n";
			if ($mapperProgLoc==1){ #bowtie2
				if ($usePairs){
					$algCmd .= "$algCmdBase -x $bwtIdxs[$kk] -1 ".join(",",@accR1) ." -2 ".join(",",@accR2);
				} else {
					$algCmd .="$algCmdBase -x $bwtIdxs[$kk] -U ".join(",",@accRS);#$paS[$iS];
				}
				$algCmd .= " $rgStr "; #--rg-id $rgID -> this is handled by getRgStr()
			} elsif ($mapperProgLoc==2){ #bwa
				die "single end mapping not implemented for bwa\n" if (!$usePairs);
				$algCmd .= $algCmdBase." -R $rgStr $REF " . join(",",@accR1). " " . join(",",@accR2); ##$pa1[$i]." ".$pa2[$i];
			} elsif ($mapperProgLoc==3){ #minimap2
				$algCmd .= $algCmdBase." -R $rgStr -a $REF ";#$MFcontstants{mini2IdxFileSuffix} ";
				if ($usePairs){ $algCmd .= join(",",@accR1) . " " . join(",",@accR2); #$pa1[$i]." ".$pa2[$i]
				} else { $algCmd .= join(" ",@accRS)
				}
			}elsif ($mapperProgLoc==4){ #kma
				$algCmd .= "$algCmdBase -t_db $REF$MFcontstants{kmaIdxFileSuffix}  ";
				if ($usePairs) {$algCmd .= " -ipe " . join(",",@accR1). " " . join(",",@accR2) ;}# $pa1[$i] $pa2[$i] " 
				else {$algCmd .= " -i " . join(" ",@accRS) ;}
				$algCmd .= " -o $tmpOutxtra[$i] ";
			}elsif ($mapperProgLoc==5){ #strobealign
				$algCmd .= "$algCmdBase $rgStr $REF  ";
				if ($usePairs) {$algCmd .=  join(",",@accR1). " " . join(",",@accR2) ;}# $pa1[$i] $pa2[$i] " 
				else {$algCmd .= " " . join(" ",@accRS) ;}
			}
			@accR1=(); @accR2=(); @accRS=(); #and empty what was already mapped

			#samtools filter and other postfiltering..
			my %postTreat = (MapperProg=>$mapperProgLoc, readTec=>$readTec, NcoreL => $NcoreL,nodeTmp =>$nodeTmp,tmpOut21=>$tmpOut21, 
										subBamsAR => \@subBams,unaligned =>$unaligned,baseN => $baseN,xtraSamSteps1 => $xtraSamSteps1,
										decoyModeActive => $decoyModeActive, map2ndTogether => $map2ndTogether,regsAR => \@regs, 
										doCram => $doCram, finalDSar => \@finalDS, outNmsAR => \@outNms, mappDirAR => \@mappDir, 
										reg_lcsAR => \@reg_lcs,
										strictHybridCoverage => ($dirsHr->{strictHybridCoverage} || 0));
			
			my ($algCmdPost,$subBamAR) = alignPostTreat(\%postTreat, $i, $kk);
			$algCmd .= $algCmdPost;
			@subBams = @{$subBamAR};
			
		}

		my $filterStep="";

		#bottleneck step: bam filtering
		#my $bamfilter = getProgPaths("bamFilter_scr");
		#my $filterStep="";
		#if ( ($mapperProgLoc==3 || $mapperProgLoc==5 ) && ($readTec eq "ONT" || $readTec eq "PB")){ #low id long reads...
		#	$filterStep = " | $smtBin view -@ $NcoreL - | $bamfilter $MFopt{bamfilterPB} | $smtBin view -@ $NcoreL -b1 - ";
		#} else { #illumina parameters
	#		$filterStep = " | $smtBin view -@ $NcoreL - | $bamfilter $MFopt{bamfilterIll} | $smtBin view -@ $NcoreL -b1 - ";
	#	}

		for (my $k=0;$k<@subBams;$k++){
			#$tmpOut22 = $tmpOut."/$outNms[$k].iniAlignment.bam";
			#print $tmpOut22."\n";
			#if (-e $tmpOut22){ next;}
			$actualMappings++;
			my $tarBam = $tmpOut22[$kk];
			$tarBam = $tmpOut22[$k] if ($decoyModeActive || $map2ndTogether>0);
			next if ($subBams[$k] eq "");#case that decoy map has already created parts of the mappings..
			my @bamParts = grep { $_ ne "" } split /\s+/, $subBams[$k];
			if (@bamParts > 1){
				$algCmd .= "\n$smtBin cat ".$subBams[$k]." $filterStep > $tarBam\n"; #$k here, because this refers to @regs
				$algCmd .= "\nrm -f ".$subBams[$k]."\n";
			} elsif (@bamParts == 1) {
				$algCmd .=  "\nmv $bamParts[0] $tarBam\n";
			}
		}
	}
	#die "$algCmd\n";
	my $mapProgNm = getMapProgNm($mapperProgLoc);

	foreach my $mpdSS (@mappDir) {
		if (!-e "$mpdSS/map.sto"){$bamFresh=1;}
	}
	$params{sortedbam}=$isSorted; $params{bamIsNew} = $bamFresh; $params{is2ndMap} = $is2ndMap;
	$params{immediateSubm} = $immediateSubm; $params{usePairs} = $anyUsedPairs;
	#@paS + @pa1 + @paBam
	#die "$actualMappings\n";
	if ((@paS + @pa1) == 0){
		print "No mapping needed, no primary reads..";
		$params{mappingStarted} = 0;
		return("","",\%params);
	}
	$params{mappingStarted} = 1;
	#print "mappingStarted :: $params{mappingStarted}\n";

	#print "$cntAli new alignments\n" if ($cntAli > 0);
	if ($bamFresh) {
		# Mapping and post-processing must execute in one scheduler allocation so
		# the intermediate alignment can remain on node-local storage.
		$params{mappingCommand} = $algCmd;
		$params{mappingDependencies} = $jDepe;
		$params{mappingJobName} = "_MAP$JNUM$supTag.$outNms[0]";
		$params{mappingScript} = $qdir.$bashN."map$supTag.sh";
		$params{mappingCores} = $Ncore;
		$params{mappingMemoryGB} = int($MFopt{MapperMemory}+1);
		$params{mappingWorkDir} = $tmpOut;
		$params{mapperNodeDir} = $nodeTmp;
		$params{unalignedWorkDir} = $unaligned;
		$params{unalignedFinalDir} = $finalUnaligned;
	} else {
		$params{mappingCommand} = "";
	}

	return("","",\%params);
}	
	










sub bamDepth{
	my ($dirsHr, $jDep,$mapparhr) = @_;
	#die "bamdep\n";
	if (!$dirsHr->{mappingStarted}
			|| (exists($mapparhr->{mappingStarted}) && !$mapparhr->{mappingStarted})){
		#die;
		return("","","");
	}
	my $readCov_Bin = getProgPaths("readCov");
	my $outName = $dirsHr->{smplName};my $ASG = $dirsHr->{assGrp};
	my $doCram =  $dirsHr->{cramAlig};
	my $mappDir = ${$dirsHr}{glbMapDir};	my $nodeTmp = ${$dirsHr}{nodeTmp}."_bamDep/$outName/";
	my $tmpOut = ${$dirsHr}{glbTmp};	my $finalD = ${$dirsHr}{outDir};
	my $qdir = $logDir; $qdir = ${$dirsHr}{qsubDir} if (exists( ${$dirsHr}{qsubDir} ));
	my $supportRds = $dirsHr->{mapSupport};
	my $supTag = ""; if ($supportRds){$supTag = ".sup";}
	my $REF = $dirsHr->{sbj}; #target to map onto, can by ","-spearated list	
	#my $bedCovBin = getProgPaths("bedCov");#"/g/bork5/hildebra/bin/bedtools2-2.21.0/bin/genomeCoverageBed";
	#my $mosDepBin = getProgPaths("mosdepth");
	my $allowDeleteMap = 1;
	#my ($par1,$par2,$parS,$liar,$rear) = getRawSeqsAssmGrp($AsgHR,$ASG,$supportRds,$outName);
	#my @libsOri = @{$liar};
	
	#readTec

	my %params = %{$mapparhr};
	$REF = $params{mappingReference} if (($params{mappingReference} || "") ne "");
	my ($isSorted , $bamFresh,$is2ndMap, $usePairs) =($params{sortedbam},$params{bamIsNew},$params{is2ndMap},$params{usePairs});
	my $mappingCommand = delete($mapparhr->{mappingCommand}) || "";
	my $mappingDependencies = $params{mappingDependencies} || $jDep;
	my $referencePreparationCommand = $params{referencePreparationCommand} || "";
	my $mappingInputSizeMB = $params{mappingInputSizeMB} || 0;
	#recreate base pars from map2tar sub
	my $locDoRmDup = $MFopt{MapperRmDup};
	$locDoRmDup = 0 if (!$usePairs);
	my $immediateSubm = $params{immediateSubm} ;
	#print "pairs used:$usePairs\n";

	my $baseN = "$outName$supTag";
	my $bashN = "";	if ($is2ndMap){$bashN = "$outName"; $bashN =~ s/,/./;}
	my $mappingRes = $tmpOut."/$baseN.iniAlignment.bam"; #this is the input, result of previous mapping steps
	my $cramSTO = "$finalD/$baseN-smd.cram.sto";
	my $nxtBAM = "$nodeTmp/$baseN-smd.bam";
	my $finalBam = "$finalD/$baseN-smd.bam";
	my $sortTMP = $nodeTmp."/$baseN.srt";
	my $sortTMP2 = $nodeTmp."/$baseN.2.srt";
 
	#check if already done
	my $outstat = check_map_done($doCram, $finalD, $baseN);
	my $outstat2 = check_depth_done($doCram, $finalD, $baseN);
	my $breakpointOut = $dirsHr->{breakpointOutput} || "";
	my $breakpointDone = $breakpointOut eq "" || -s $breakpointOut;
	my $breakpointWork = ($breakpointOut ne "" && !$breakpointDone)
		? "$nodeTmp/".basename($breakpointOut) : "";
	#die "$outstat $outstat2 $mappingRes $nxtBAM\n$mappDir\n$finalD\n";
	if ($outstat2 && $outstat && $breakpointDone){return ("","",$outstat);}

	#die "$mappDir/$baseN-smd.bam.coverage.gz\n$finalD/$baseN-smd.bam.coverage.gz\n" ;#if (-e "$finalD/$baseN-smd.bam.coverage.gz");
		
	my $numCore = ${$dirsHr}{sortCores};#new functionality with sambamba
	$numCore = 1 if (!$numCore || $numCore < 1);
	my $locSrtMem = $MFopt{mapSortMemGb};
	my $baseMem=20;
	if ($locSrtMem < 0){
		$locSrtMem = $baseMem + (2 * $mappingInputSizeMB/1024);
		$locSrtMem += $baseMem if ($MFopt{largeMapperDB});
	}
	# Duplicate removal normally streams name-sort -> fixmate -> coordinate-sort.
	# Both sorts are then alive at once and -m applies per thread.  For large
	# inputs, materialize the name-sorted BAM so only one sort is resident at a
	# time.  The mapping scratch request already budgets eight times compressed
	# input size, which leaves room for this intermediate on local SSD.
	my $largeSortInputMB = 4 * 1024;
	my $serialiseDuplicateSorts = $locDoRmDup
		&& $mappingInputSizeMB >= $largeSortInputMB;
	my $sortProcessCount = ($locDoRmDup && !$serialiseDuplicateSorts) ? 2 : 1;
	my $sortMemoryMB = int((($locSrtMem * 1024) - 2048)
		/ ($numCore * $sortProcessCount));
	$sortMemoryMB = 256 if ($sortMemoryMB < 256);
	# Large sorts create substantial buffers outside the nominal -m arena.
	# Keep their per-thread arena at samtools' conservative default even when an
	# automatic or user-provided total budget would permit a much larger value.
	my $sortMemoryCapMB = $serialiseDuplicateSorts ? 768 : 2048;
	$sortMemoryMB = $sortMemoryCapMB if ($sortMemoryMB > $sortMemoryCapMB);
	# A completed mapping may still predate breakpoint output.  Repair that
	# single derivative from the canonical coverage without trying to sort a
	# node-local alignment that no longer exists.
	if ($outstat && $outstat2 && !$breakpointDone && $mappingCommand eq "") {
		my $breakpointScr = getProgPaths("breakpoints_scr");
		my $breakpointDir = dirname($breakpointOut);
		my $breakpointName = basename($breakpointOut);
		my $breakpointStage = "$breakpointDir/.$breakpointName.stage";
		my $breakpointCoverage = -s "$finalBam.coverage.gz"
			? "$finalBam.coverage.gz" : "$finalBam.coverage";
		my $breakpointCmd = "rm -rf $nodeTmp $breakpointStage\n"
			."mkdir -p $nodeTmp $breakpointDir $breakpointStage\n"
			.$referencePreparationCommand
			."test -s $breakpointCoverage\n";
		$breakpointCmd .= "$pigzBin -t $breakpointCoverage\n" if ($breakpointCoverage =~ /\.gz$/);
		$breakpointCmd .= "$breakpointScr --assembly $REF --coverage $breakpointCoverage --output $breakpointWork "
			."--breakpoint-depth $MFopt{breakpointDepth} --min-breakpoint-length $MFopt{breakpointMinLength} "
			."--smooth-gap $MFopt{breakpointSmoothGap} --flank-length $MFopt{breakpointFlankLength} "
			."--min-flank-depth $MFopt{breakpointMinFlankDepth} --max-flank-fraction $MFopt{breakpointMaxFlankFraction}\n"
			."test -s $breakpointWork\n"
			."mv $breakpointWork $breakpointStage/$breakpointName\n"
			."mv -f $breakpointStage/$breakpointName $breakpointOut\n"
			."rmdir $breakpointStage\nrm -rf $nodeTmp";
		$breakpointCmd .= " $params{mapperNodeDir}" if (($params{mapperNodeDir} || "") ne "");
		$breakpointCmd .= "\n";
		my ($breakpointJob, $breakpointSubmission) = ($jDep, "");
		if (${$dirsHr}{submit}) {
			($breakpointJob, $breakpointSubmission) = qsubSystem(
				$qdir.$bashN."breakpoint$supTag.sh", $breakpointCmd,
				$numCore, "20G", "_MAP$JNUM$supTag.$outName",
				$mappingDependencies, "", $immediateSubm,
				$QSBoptHR->{General_Hosts}, $QSBoptHR);
		} else {
			$breakpointSubmission = $breakpointCmd;
		}
		return ($breakpointJob, $breakpointSubmission, $outstat);
	}
	# A complete canonical BAM/CRAM is authoritative. If only coverage (and
	# possibly its breakpoint derivative) is missing, derive those files directly
	# without sorting, duplicate removal, or CRAM conversion a second time.
	if ($outstat && !$outstat2 && $mappingCommand eq "") {
		my $canonicalAlignment = $doCram
			? "$finalD/$baseN-smd.cram" : "$finalD/$baseN-smd.bam";
		my $repairPrefix = "$nodeTmp/$baseN-smd.bam";
		my $repairCoverage = "$repairPrefix.coverage.gz";
		my $repairStage = "$finalD/.$baseN.coverage-stage";
		my $sam2bed = getProgPaths("samcov2bed");
		my $depthFilters = $dirsHr->{strictHybridCoverage}
			? "-Q $MFopt{hybridMinMapQ} -q $MFopt{hybridMinBaseQ}" : "";
		my $referenceOption = $doCram ? "--reference $REF " : "";
		my $repairCmd = "rm -rf $nodeTmp $repairStage\nmkdir -p $nodeTmp $finalD $repairStage\n";
		$repairCmd .= $referencePreparationCommand;
		$repairCmd .= "test -s $canonicalAlignment\n$smtBin quickcheck $canonicalAlignment\n";
		$repairCmd .= "$smtBin depth ${referenceOption}-aa $depthFilters -@ $numCore $canonicalAlignment "
			."| $sam2bed | $pigzBin -p $numCore -c > $repairCoverage\n";
		$repairCmd .= "test -s $repairCoverage\n$pigzBin -t $repairCoverage\n";
		$repairCmd .= jgi_depth_cmd([$canonicalAlignment],$repairPrefix,95,$numCore,$REF)
			if ($MFopt{DoJGIcoverage});
		if (!$breakpointDone) {
			my $breakpointScr = getProgPaths("breakpoints_scr");
			my $breakpointDir = dirname($breakpointOut);
			my $breakpointStage = "$breakpointOut.stage";
			$repairCmd .= "mkdir -p $breakpointDir\n";
			$repairCmd .= "$breakpointScr --assembly $REF --coverage $repairCoverage --output $breakpointWork "
				."--breakpoint-depth $MFopt{breakpointDepth} --min-breakpoint-length $MFopt{breakpointMinLength} "
				."--smooth-gap $MFopt{breakpointSmoothGap} --flank-length $MFopt{breakpointFlankLength} "
				."--min-flank-depth $MFopt{breakpointMinFlankDepth} --max-flank-fraction $MFopt{breakpointMaxFlankFraction}\n";
			$repairCmd .= "test -s $breakpointWork\nrm -f $breakpointStage\n"
				."mv $breakpointWork $breakpointStage\nmv -f $breakpointStage $breakpointOut\n";
		}
		$repairCmd .= "for f in $repairPrefix*; do [ -e \"\$f\" ] || continue; mv \"\$f\" $repairStage/; done\n";
		$repairCmd .= "test -s $repairStage/$baseN-smd.bam.coverage.gz\n"
			."for f in $repairStage/*; do mv -f \"\$f\" $finalD/; done\nrmdir $repairStage\n";
		$repairCmd .= "rm -rf $nodeTmp";
		$repairCmd .= " $params{mapperNodeDir}" if (($params{mapperNodeDir} || "") ne "");
		$repairCmd .= "\n";
		my ($repairJob, $repairSubmission) = ($jDep, "");
		if (${$dirsHr}{submit}) {
			my $repairScratch = $QSBoptHR->{tmpSpace};
			my $baseMapHDD = $HDDspace{mapping}; $baseMapHDD =~ s/G$//;
			$QSBoptHR->{tmpSpace} = int((2.0 * $mappingInputSizeMB * $baseMapHDD) / 1024) + 10 ."G";
			($repairJob,$repairSubmission) = qsubSystem(
				$qdir.$bashN."coverageRepair$supTag.sh",$repairCmd,$numCore,
				(int($locSrtMem)+1)."G","_COV$JNUM$supTag.$outName",
				$mappingDependencies,"",$immediateSubm,$QSBoptHR->{General_Hosts},$QSBoptHR);
			$QSBoptHR->{tmpSpace} = $repairScratch;
		} else {
			$repairSubmission = $repairCmd;
		}
		return ($repairJob,$repairSubmission,$outstat);
	}
	#my $biobambamBin = "bamsormadup";
	#depth profile
	#my $cmd = "$smtBin faidx $REF\n $smtBin view -btS $REF.fai $tmpOut/$baseN.sam > $tmpOut/$baseN.bam\n";
	my $cmd = "sleep 1\n";
	#$cmd .= "mkdir -p $nodeTmp\n";
	$cmd .= "mkdir -p $mappDir $tmpOut $nodeTmp\n";
	#sort & remove duplicates
	if ($outstat && !$outstat2 ){
		if ($doCram){
			$cmd .= "#uncramming already stored results..\n" . cram2bsam("$finalD/$baseN-smd.cram",$REF,$mappingRes,1,$numCore) ."\n" ;
		} elsif (-e "$finalD/$baseN-smd.bam") {
			$mappingRes = "$finalD/$baseN-smd.bam";
			$allowDeleteMap = 0;
		} else {
			die "bamDepth:: canonical mapping output is missing..\n";
		}
	}

	if ($locDoRmDup){
	#bit complicated in samtools now: first sort by name, then fixmate (samtools 1.8)
	#".int($MFopt{mapSortMemGb}/$numCore)."G
		$cmd .= "echo \"samtools sort budget: ${sortMemoryMB}M x $numCore threads x $sortProcessCount concurrent sorts\"\n";
		$cmd .= "echo \"Sorting .bam and fixing mates ...\"\n";
		#$cmd .= "$smtBin sort -n -m 768M -T $sortTMP -@ $numCore $mappingRes | $smtBin fixmate -m -@ $numCore - - | $smtBin sort -T $sortTMP2 -m 768M -@ $numCore - | $smtBin markdup -s -r -@ $numCore - $nxtBAM\n";
		if ($serialiseDuplicateSorts) {
			my $nameSortedBam = "$nodeTmp/$baseN.name-sorted.bam";
			$cmd .= "echo \"Large mapping input: serializing name and coordinate sorts\"\n";
			$cmd .= "$smtBin sort -n -m ${sortMemoryMB}M -u -T $sortTMP -@ $numCore -o $nameSortedBam $mappingRes;\n";
			$cmd .= "$smtBin fixmate -m -@ $numCore -u $nameSortedBam - | $smtBin sort -T $sortTMP2 -m ${sortMemoryMB}M -@ $numCore -o $nxtBAM.nf;\n";
			$cmd .= "rm -f $nameSortedBam;\n";
		} else {
			$cmd .= "$smtBin sort -n -m ${sortMemoryMB}M -u -T $sortTMP -@ $numCore $mappingRes | $smtBin fixmate -m -@ $numCore -u - - | $smtBin sort -T $sortTMP2 -m ${sortMemoryMB}M -@ $numCore -o $nxtBAM.nf;\n";
		}
		$cmd .= "echo \"marking duplicates ...\"\n";
		$cmd .= "$smtBin markdup --use-read-groups --no-multi-dup -d 1000 -T $sortTMP -s -r -O BAM -@ $numCore $nxtBAM.nf $nxtBAM;\n";
		$cmd .= "rm -f $nxtBAM.nf;\n";
	} else {
		if ($isSorted ){
			$cmd .= "mv $mappingRes $nxtBAM\n";
		} else {
			$cmd .= "echo \"Sorting .bam ...\"\n";
			$cmd .= "$smtBin sort -@ $numCore -m ${sortMemoryMB}M -O BAM -T $sortTMP2 $mappingRes > $nxtBAM\n";
		}
	}
	$cmd .= "echo \"Building .bam index...\"\n";
	$cmd .= "$smtBin index -@ $numCore $nxtBAM\n";
	$cmd .= "[ -s $nxtBAM ] || exit 3\n";#check if output is empty (samtools crashed?)
	$cmd .= "rm -f $mappingRes \n" if (!$isSorted && $allowDeleteMap);

	#$cmd .= "$smtBin view -bS $tmpOut21 > $tmpOut/$baseN.bam\n";
	#my $mdJar = "/g/bork5/hildebra/bin/picard-tools-1.119/MarkDuplicates.jar";
	#$cmd .= "java -Xms1g -Xmx24g -XX:ParallelGCThreads=2 -XX:MaxPermSize=1g -XX:+CMSClassUnloadingEnabled ";		$cmd .= "-jar $mdJar INPUT=$mappDir/$baseN-s.bam OUTPUT=$mappDir/$baseN-smd.bam METRICS_FILE=$mappDir/$baseN-smd.metrics ";		
	#$cmd .= "AS=TRUE VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 REMOVE_DUPLICATES=TRUE\n\n";
	
	#$cmd .= "$smtBin index $nxtBAM\n";#$mappDir/$baseN-smd.bam\n";
	#takes too much space, calc on the fly!
	#$cmd .= "$smtBin depth $mappDir/$baseN-s.bam > $mappDir/$baseN.depth.bam.txt\n";#depth file (different estimator)
	my $covCmd = "";
	if (!-e "$nxtBAM.coverage.gz"){
		$covCmd .= "echo \"Creating coverage of reads...\"\n";
		#$covCmd .= "$mosDepBin -t 4 -x $nxtBAM.coverage $nxtBAM\n";
		#$covCmd .= "mv $nxtBAM.coverage.per-base.bed.gz $nxtBAM.coverage.gz\n";
		#$covCmd = "$bedCovBin -ibam $nxtBAM -bg > $nxtBAM.coverage\n";
		#$covCmd .= "rm -f $nxtBAM.coverage.gz\n$pigzBin -f -p $numCore $nxtBAM.coverage\n";
		my $sam2bed = getProgPaths("samcov2bed");
		my $depthFilters = $dirsHr->{strictHybridCoverage}
			? "-Q $MFopt{hybridMinMapQ} -q $MFopt{hybridMinBaseQ}" : "";
		$covCmd .= "$smtBin depth -aa $depthFilters -@ $numCore $nxtBAM | $sam2bed | $pigzBin -p $numCore -c > $nxtBAM.coverage.gz\n";
		#$covCmd .= "rm -f $nxtBAM.coverage.gz\n$pigzBin -f -p $numCore $nxtBAM.coverage\n";
		if ($map2ndMpde == 3){ #map unmapped to genecat
			$covCmd .= "$readCov_Bin $nxtBAM.coverage.gz - 100\n"; #read length doesn't matter, since no win is given out and no gene is present
		}
	}
	if ($breakpointOut ne "" && !-s $breakpointOut) {
		my $breakpointScr = getProgPaths("breakpoints_scr");
		my $breakpointDir = dirname($breakpointWork);
		my $breakpointCoverage = -s "$finalBam.coverage.gz" ? "$finalBam.coverage.gz"
			: (-s "$finalBam.coverage" ? "$finalBam.coverage" : "$nxtBAM.coverage.gz");
		$covCmd .= "mkdir -p $breakpointDir\n";
		$covCmd .= "$breakpointScr --assembly $REF --coverage $breakpointCoverage --output $breakpointWork "
			."--breakpoint-depth $MFopt{breakpointDepth} --min-breakpoint-length $MFopt{breakpointMinLength} "
			."--smooth-gap $MFopt{breakpointSmoothGap} --flank-length $MFopt{breakpointFlankLength} "
			."--min-flank-depth $MFopt{breakpointMinFlankDepth} --max-flank-fraction $MFopt{breakpointMaxFlankFraction}\n";
		$covCmd .= "[ -s $breakpointWork ] || exit 34\n";
	}
	#die "$covCmd\n";
	
	#jgi depth profile - not required, different call to metabat
	$covCmd .= jgi_depth_cmd([$nxtBAM],$nxtBAM,95) if ($MFopt{DoJGIcoverage});
	#$covCmd .= "/g/bork5/hildebra/bin/bedtools2-2.21.0/bin/genomeCoverageBed -ibam $nxtBAM -bg  | ".'awk \'BEGIN {pc=""} {	c=$1;	if (c == pc) {		cov=cov+$2*$5;	} else {		print pc,cov;		cov=$2*$5;	pc=c}';
	#$covCmd .= "} END {print pc,cov}\' $nxtBAM.coverage | tail -n +2 > $nxtBAM.coverage.percontig";
	#$covCmd .= "rm -f $nxtBAM.bai;\n";
	my ($CRAMcmd,$CRAMf) = bam2cram($nxtBAM,$REF,1,$doCram,"", $numCore);
	$CRAMcmd = "echo \"Building .cram ...\"\n$CRAMcmd" if ($CRAMcmd ne "");
	my $publishStage = "$finalD/.$baseN.mapping-stage";
	$CRAMcmd .= "\n# Publish the complete mapping from node-local work into its canonical directory.\n";
	$CRAMcmd .= "rm -rf $publishStage\nmkdir -p $publishStage $finalD\n";
	$CRAMcmd .= "for f in $CRAMf*; do [ -e \"\$f\" ] || continue; mv \"\$f\" $publishStage/; done\n"
		unless ($CRAMf eq "");
	if ($breakpointWork ne "") {
		# breakpointWork shares the $nxtBAM prefix. Exclude it from the generic
		# sidecar loop so the required, validated move below does not publish it twice.
		$CRAMcmd .= "for f in $nxtBAM*; do [ -e \"\$f\" ] || continue; [ \"\$f\" = \"$breakpointWork\" ] && continue; mv \"\$f\" $publishStage/; done\n";
		$CRAMcmd .= "test -s $breakpointWork\nmv $breakpointWork $publishStage/\n";
	} else {
		$CRAMcmd .= "for f in $nxtBAM*; do [ -e \"\$f\" ] || continue; mv \"\$f\" $publishStage/; done\n";
	}
	if ($doCram) {
		$CRAMcmd .= "test -s $publishStage/$baseN-smd.cram\n$smtBin quickcheck $publishStage/$baseN-smd.cram\n";
	} else {
		$CRAMcmd .= "test -s $publishStage/$baseN-smd.bam\n$smtBin quickcheck $publishStage/$baseN-smd.bam\n";
	}
	$CRAMcmd .= "test -s $publishStage/$baseN-smd.bam.coverage.gz\n$pigzBin -t $publishStage/$baseN-smd.bam.coverage.gz\n";
	# Record which concrete assembly file this alignment targets.  This is a
	# cheap durable discriminator between a hybrid preassembly and the final
	# assembly even though both occupy the same path at different workflow
	# stages.  Secondary multi-reference mappings are not part of that state
	# transition and may have comma-separated references.
	if (!$is2ndMap && $REF !~ /,/) {
		$CRAMcmd .= "test -s \"$REF\"\n";
		$CRAMcmd .= "stat -c '%s %Y' \"$REF\" > $publishStage/$baseN-smd.reference.stat\n";
		$CRAMcmd .= "test -s $publishStage/$baseN-smd.reference.stat\n";
	}
	$CRAMcmd .= "for f in $publishStage/*; do mv -f \"\$f\" $finalD/; done\nrmdir $publishStage\n";
	if (($params{unalignedFinalDir} || "") ne "") {
		my $unalignedStage = "$finalD/.$baseN.unaligned-stage";
		$CRAMcmd .= "if [ -d $params{unalignedWorkDir} ]; then rm -rf $unalignedStage; mkdir -p $unalignedStage; "
			."cp -a $params{unalignedWorkDir}/. $unalignedStage/; rm -rf $params{unalignedFinalDir}; "
			."mv $unalignedStage $params{unalignedFinalDir}; fi\n";
	}
	# Completion markers are the publication commit record and must be last.
	$CRAMcmd .= "touch $cramSTO\n" if ($doCram);
	$CRAMcmd .= "echo \"".basename($nxtBAM)."\" > $finalD/done.sto\n" if(!$supportRds || !$map{$curSmpl}{hasPrimaryRds});


	my $newJobN = $params{mappingJobName} || "_MAP$JNUM$supTag.$outName";
	#subsequent jobs not dependent on this one
	my $jobN2 = $jDep; my $retCmds="";
	$covCmd = "" if (fileGZs("$nxtBAM.coverage") && $breakpointDone);
	
	#my $cleaner = "mv $tmpD $fin
	if (0){
		print "B1 " if ( $doCram && !-e $cramSTO); print "B2 " if (  !-s "$nxtBAM.coverage"); print "B3 " if ( $bamFresh);
		print " ".$cramSTO."\n";
	}
	my $nodeCln = "\nrm -rf $nodeTmp";
	# Multi-reference mapping reuses the mapper's node-local alignments across
	# several bamDepth command blocks. Its caller removes shared work only after
	# the last reference has been published.
	unless ($dirsHr->{deferMappingCleanup}) {
		$nodeCln .= " $params{mappingWorkDir}" if (($params{mappingWorkDir} || "") ne "");
		$nodeCln .= " $params{mapperNodeDir}" if (($params{mapperNodeDir} || "") ne "");
	}
	$nodeCln .= ";\necho \"DONE mapping\"\n";
	#die "$cmd\n$covCmd\n$CRAMcmd\n";
	if ( $mappingCommand ne "" || ($doCram && !-e $cramSTO) || (!fileGZs("$finalBam.coverage")) || !$breakpointDone ){#|| $bamFresh){
		my $preHDDspace=$QSBoptHR->{tmpSpace};		my $baseMapHDD = $HDDspace{mapping} ;  $baseMapHDD =~ s/G$//;
		$QSBoptHR->{tmpSpace} = int((2.0 * $mappingInputSizeMB*$baseMapHDD) /1024)+30  ."G";		if (${$dirsHr}{submit}){
		#die "map2:: $qdir\n$cramSTO\n$nxtBAM.coverage\n";
		my $combinedCores = $numCore > ($params{mappingCores} || 0) ? $numCore : ($params{mappingCores} || $numCore);
		my $combinedMem = int($locSrtMem)+1;
		my $controlledSortMB = $sortMemoryMB * $numCore * $sortProcessCount;
		# -m is approximate and excludes compression threads, fixmate, pipes, Perl,
		# libraries and allocator fragmentation.  Add 25% plus 4 GiB rather than
		# requesting only the nominal samtools arenas.
		my $sortRequiredMem = int((($controlledSortMB * 1.25) + 4096 + 1023) / 1024);
		$combinedMem = $sortRequiredMem if ($sortRequiredMem > $combinedMem);
		$combinedMem = $params{mappingMemoryGB} if (($params{mappingMemoryGB} || 0) > $combinedMem);
		my $combinedScript = $params{mappingScript} || $qdir.$bashN."map$supTag.sh";
		($jobN2,$retCmds) = qsubSystem($combinedScript,
				$mappingCommand."\n".$cmd."\n".$covCmd."\n".$CRAMcmd."\n$nodeCln\n"
				,$combinedCores,  $combinedMem."G",$newJobN,$mappingDependencies,"",$immediateSubm,$QSBoptHR->{General_Hosts},$QSBoptHR);
		$QSBoptHR->{tmpSpace}=$preHDDspace;
		} else {
			$cmd =~ s/sleep \d+//;
			$retCmds = $mappingCommand."\n".$cmd."\n".$covCmd."\n".$CRAMcmd."\n$nodeCln\n";
		}
	} 
	#die();
	return($jobN2,$retCmds,$outstat);
}

sub mergeMP2Table($){
	my ($dir_MP2) = @_;
	return if (!$MFopt{DoMetaPhlan});
	if ($progStats{metaPhl2FailCnts}){
		print "$progStats{metaPhl2FailCnts} / ". ($progStats{metaPhl2ComplCnts}+$progStats{metaPhl2FailCnts}) ." samples with incomplete Metaphlan assignments\n";
		return;
	}
	my $outD = $dir_MP2;
	$outD =~ s/[^\/]+\/?$//;

	my $mergeTblScript = getProgPaths("metPhl2Merge");#"/g/bork3/home/hildebra/bin/metaphlan2/utils/merge_metaphlan_tables.py";
	my $MP2Mstone = "$outD/MP2.cnt.stone";
	my $prevCnts = 0;
	if (-e $MP2Mstone){
		$prevCnts = `cat $MP2Mstone`; chomp $prevCnts; 
	}
	my $redoMP2Tables = 0;
	print "\nAll samples ($progStats{metaPhl2ComplCnts}) have metaphlan assignments.\n";
	if ($prevCnts < $progStats{metaPhl2ComplCnts}){
		print "Redoing metaphlan2 table merge ,due to higher number of samples detected ($progStats{metaPhl2ComplCnts}, prev: $prevCnts)\n";
		$redoMP2Tables = 1;
	} else {
		return;
	}
	my $getHDerCmd = "head -n1 ";
	if ($MFopt{DoMetaPhlan} >= 3){
		$getHDerCmd = "head -n2 ";
	}
	my @slvl = ("k__","p__","c__","o__","f__","g__","s__"); my @llvl = ("kingdom","phylum","class","order","family","genus","species");
	my $mrgCmd = ""; 
	my $tmprawF = " $outD/MePh.Raw.mat";
	my $premrgCmd = "$mergeTblScript $dir_MP2/*MP2.txt > $tmprawF\n" ;
	for (my $i=0;$i<@slvl;$i++){
		my $outF = "$outD/MePh.all.$llvl[$i].mat";
		if (!-e "$outF" || $redoMP2Tables){
			unless ($premrgCmd eq ""){$mrgCmd.=$premrgCmd;$premrgCmd="";	}
			$mrgCmd .= "$getHDerCmd  $tmprawF | tail -n1 > $outF\n"; #needed for version 3
			$mrgCmd .= "cat  $tmprawF | grep \"".$slvl[$i]."[^\\|]*\\s\" | sed 's/|/;/g' >> $outF\n" ;
		}
	}
	
	$mrgCmd .= "echo \"$progStats{metaPhl2ComplCnts}\" > $MP2Mstone";

	#systemW "$mrgCmd";
	my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0; 
	my ($jobN, $tmpCmd) = qsubSystem($baseOut."/LOGandSUB/MP2merg.sh",$mrgCmd,1,"80G","MP2mrg","","",1,[],$QSBoptHR) ;
	$QSBoptHR->{tmpSpace} =$tmpSHDD;
	return;
	
	
	#this part is no longer used since MATAFILER v 0.23
	#die "$mrgCmd\n";

	$premrgCmd = "$mergeTblScript $dir_MP2/*MP2.noV.noB.txt > $dir_MP2/All.MP2.noV.noB.mat\n";
	for (my $i=0;$i<@slvl;$i++){
		my $outF = "$dir_MP2/MP2.noV.noB.$llvl[$i].mat";
		if (!-e "$outF" || $redoMP2Tables){
			unless ($premrgCmd eq ""){$mrgCmd.=$premrgCmd;$premrgCmd="";	}
			$mrgCmd .= "$getHDerCmd  $dir_MP2/tmp.All.MP2.mat > $outF\n";
			$mrgCmd .= "cat $dir_MP2/All.MP2.noV.noB.mat | grep \"".$slvl[$i]."[^\\|]*\\s\" | sed 's/|/;/g' >> $outF\n"  unless (-e "$outF"&& !$redoMP2Tables);
		}
	}
	$mrgCmd .= "rm -f $dir_MP2/All.MP2.noV.noB.mat\n\n";
	$premrgCmd = "$mergeTblScript $dir_MP2/*MP2.noV.txt > $dir_MP2/All.MP2.noV.mat\n";
	for (my $i=0;$i<@slvl;$i++){
		my $outF = "$dir_MP2/MP2.noV.$llvl[$i].mat";
		if (!-e "$outF" || $redoMP2Tables){
			unless ($premrgCmd eq ""){$mrgCmd.=$premrgCmd;$premrgCmd="";	}
			$mrgCmd .= "$$getHDerCmd  $dir_MP2/tmp.All.MP2.mat > $outF\n";
			$mrgCmd .= "cat $dir_MP2/All.MP2.noV.mat | grep \"".$slvl[$i]."[^\\|]*\\s\" | sed 's/|/;/g'  >> $outF\n"  unless (-e "$outF"&& !$redoMP2Tables);
		}
	}
	$mrgCmd .= "rm -f $dir_MP2/All.MP2.noV.mat\n\n";
	#viruses
	$premrgCmd = "$mergeTblScript $dir_MP2/*MP2.VirusOnly.txt > $dir_MP2/All.MP2.VirusOnly.mat\n";
	for (my $i=0;$i<@slvl;$i++){
		my $outF = "$dir_MP2/MP2.VirusOnly.$llvl[$i].mat";
		if (!-e "$outF" || $redoMP2Tables){
			unless ($premrgCmd eq ""){$mrgCmd.=$premrgCmd;$premrgCmd="";	}
			$mrgCmd .= "$$getHDerCmd  $dir_MP2/tmp.All.MP2.mat > $outF\n";
			$mrgCmd .= "cat $dir_MP2/All.MP2.VirusOnly.mat | grep \"".$slvl[$i]."[^\\|]*\\s\" | sed 's/|/;/g'  >> $outF\n"  unless (-e "$outF"&& !$redoMP2Tables);
		}
	}
	$mrgCmd .= "rm -f $dir_MP2/tmp.All.MP2.mat $dir_MP2/All.MP2.VirusOnly.mat\n\n";

	#die $mrgCmd."\n";
}

sub prepMetaphlan{
	#die;
	return unless ($MFopt{DoMetaPhlan});
	#return; #takes too much time..
	my $metaPhlBin = getProgPaths("metPhl2");
	print "Checking metaphlan version .. ";
	my $vstr = "";
	$vstr = `$metaPhlBin --version 2>/dev/null`;
	$vstr =~ m/version ([\.\d]+)/;
	#print "$1 ";
	my $MPver = $1; my $MPverL = int(substr($MPver,0,1));
	if ($MFopt{DoMetaPhlan}){
		$MFopt{DoMetaPhlan} = 1;
		if ( $MPverL < 3.0){ print "Metaphlan version below 3, reinstall updated metaphlan3\n"; exit(33);}
		print "will use MetaPhlan version $MPver\n";
		$MFopt{DoMetaPhlan} = $MPverL;
	}
}


sub metphlanMapping{
	my ($tmpD,$finOutD,$smp,$Ncore,$deps) = @_;
	
	my $cleanSeqSetHR = sampleReadSet($curSmpl, "clean");

	my $libraries = readLibrariesByScope($cleanSeqSetHR, 'primary', 1, $curSmpl);

	
	my $bwt2Bin = getProgPaths("bwt2");#"/g/bork5/hildebra/bin/bowtie2-2.2.9/bowtie2";
	my $metPhlaBin = getProgPaths("metPhl2");#"/g/bork3/home/hildebra/bin/metaphlan2/metaphlan2.py";
	#path to metaphlan DB
	my $mpDB = getProgPaths("metPhl2_db",0);#metaphlan2/db_v20/mpa_v20_m200
	my $metPhl2Merge = getProgPaths("metPhl2Merge");#"/g/bork3/home/hildebra/bin/metaphlan2/utils/merge_metaphlan_tables.py";
	
	#split metaphlan DB path for v3
	$mpDB =~ m/^(.*)\/([^\/]+)$/;
	my $mpDB1 = $1."/";my $mpDB2 = $2;
	
	my $stone = $finOutD."$smp.MP2.sto";
	if (!$MFopt{DoMetaPhlan} ||  -e $stone){return;}
	my $pairs = libraryPairs($libraries);
	my @car1 = map { $_->{files}{r1} } @{$pairs};
	my @car2 = map { $_->{files}{r2} } @{$pairs};
	my @sar = @{libraryFiles($libraries, 'single')};
	my $inF1 = join(",",@car1); my $inF2 = join(",",@car2); my $inFS = join(",",@sar); 
	system "mkdir -p $finOutD\n" unless (-d $finOutD);
	my $finOut = $finOutD."$smp.MP2.txt";
	my $finOut_noV = $finOutD."$smp.MP2.noV.txt";
	my $finOut_noVB = $finOutD."$smp.MP2.noV.noB.txt";
	my $finOut_Vo = $finOutD."$smp.MP2.VirusOnly.txt";
	my $sam = "$tmpD/metph2.sam";
	my $vXparams = "";#my $v2params = ""; 
	if ($MFopt{DoMetaPhlan} >=4){ #version 4+
		$vXparams = " --sample_id $smp --nproc $Ncore --unclassified_estimation --add_viruses --nreads \$readN -o " ; #if ($MFopt{DoMetaPhlan3}); #--unknown_estimation -> requires --nreads
	} elsif ($MFopt{DoMetaPhlan} >=3){ #version 3+
		$vXparams = " --sample_id $smp --nproc $Ncore --unknown_estimation --nreads \$readN -o " ; #if ($MFopt{DoMetaPhlan3}); #--unknown_estimation -> requires --nreads
	} else { #version 2
		$vXparams = "--ignore_viruses >";# if (!$MFopt{DoMetaPhlan3});
	}
	my $taxinfo = "--mpa_pkl $mpDB.pkl ";
	$taxinfo = "--bowtie2db $mpDB1 -x $mpDB2 " if ($MFopt{DoMetaPhlan} >=3);
	my $qsubFile = $logDir."metaPhl$MFopt{DoMetaPhlan}.sh";
	
	#$ bowtie2 --sam-no-hd --sam-no-sq --no-unal --very-sensitive -S metagenome.sam -x metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103  -U metagenome.fastq
	my $cmd = "mkdir -p $tmpD\n$bwt2Bin --sam-no-hd --sam-no-sq --no-unal --very-sensitive -S $sam -p $Ncore -x $mpDB ";
	$cmd .= "-1 $inF1 -2 $inF2 " if (@car1 > 0 );;
	$cmd .= "-U $inFS " if (@sar > 0);
	$cmd .= "\n\n";
	$cmd .=  "sleep 1\nreadN=\$(grep -v 'Warning:' $qsubFile.etxt |  head -n1 | cut -f1 -d' ')\necho \$readN\n";
	$cmd .= "$metPhlaBin $sam --input_type sam $taxinfo $vXparams $finOut\n";
	#$cmd .= "$metPhlaBin $sam --input_type sam $taxinfo $v2params $v3params $finOut_noV\n" if (!$MFopt{DoMetaPhlan3});
	#$cmd .= "$metPhlaBin $sam --input_type sam --ignore_bacteria --ignore_archaea $taxinfo $v2params $v3params $finOut_noVB\n";
	#$cmd .= "$metPhlaBin $sam --input_type sam --ignore_bacteria --ignore_eukaryotes --ignore_archaea $taxinfo $v2params $v3params $finOut_Vo\n";
	$cmd .= "rm -f $sam\n";
	my $mergeStr = "$metPhl2Merge *.MP2.txt > $finOutD/comb.MP2.txt";
	$cmd .= "[ -s $finOut ] || exit 4\n";
	$cmd .= "echo \' $mergeStr \' > $stone\n";
	my $jobN = "MP$MFopt{DoMetaPhlan}$JNUM";

	# Bowtie2 materialises an uncompressed SAM in node-local storage. Its size
	# can substantially exceed the compressed FASTQ inputs, so do not leave this
	# submission on the generic per-job scratch default.
	my $previousTmpSpace = $QSBoptHR->{tmpSpace};
	$QSBoptHR->{tmpSpace} =
		int(($map{$curSmpl}{inputFileSizeMB} * 6) / 1024) + 15 ."G";
	my ($jobN2,$tmpCmd) = qsubSystem($qsubFile,
			$cmd,$Ncore,"3G",$jobN,$deps,"",1,[],$QSBoptHR);
	$QSBoptHR->{tmpSpace} = $previousTmpSpace;
	$jobN  = $jobN2;
	return $jobN;
}

sub TaxaTarget{
	my ($tmpD,$finOutD,$smp,$Ncore,$deps) = @_;
	
	my $cleanSeqSetHR = sampleReadSet($curSmpl, "clean");
	my $libraries = readLibrariesByScope($cleanSeqSetHR, 'primary', 1, $curSmpl);

	my $stone = $finOutD."$smp.TaxTar.sto";
	if (!$MFopt{DoTaxaTarget} ||  -e $stone){return;}

	my $pairs = libraryPairs($libraries);
	my @car1 = map { $_->{files}{r1} } @{$pairs};
	my @car2 = map { $_->{files}{r2} } @{$pairs};
	my @sar = @{libraryFiles($libraries, 'single')};
	die "TaxaTarget requires at least one paired-read library for sample $smp\n" if (!@car1);
	die "TaxaTarget does not support singleton reads for sample $smp; disable it or provide paired reads only\n" if (@sar);
	my $taxTBin = getProgPaths("TaxaTarget");#"/g/bork5/hildebra/bin/bowtie2-2.2.9/bowtie2";
	my $taxTarDir = $taxTBin;
	$taxTarDir =~ s/run_pipeline_scripts\/run_protist_pipeline_fda.py//;
	$taxTarDir =~ s/python //;
	my $sampleOut = "$finOutD/$smp";
	my $cmd = "mkdir -p $sampleOut $tmpD\n";
	my @library_outputs;
	for (my $i=0; $i<@car1; $i++) {
		my $libraryOut = "$sampleOut/lib$i";
		push @library_outputs, $libraryOut;
		$cmd .= "mkdir -p $libraryOut\n";
		$cmd .= "$taxTBin -r $car1[$i] -r2 $car2[$i] -e $taxTarDir/run_pipeline_scripts/environment.txt --tmp -t $Ncore -o $libraryOut\n";
		$cmd .= "find $libraryOut -type f -size +0c -print -quit | grep -q .\n";
	}
	$cmd .= "printf '%s\\n' '".join("' '", @library_outputs)."' > $stone\n";
	my $jobN = "TT$JNUM";

	my ($jobN2,$tmpCmd) = qsubSystem($logDir."taxtar.sh",
			$cmd,$Ncore,"3G",$jobN,$deps,"",1,[],$QSBoptHR);
			
	#die $logDir;
	return $jobN2;
}



sub remComma($){
	my ($in) = @_;
	$in =~ s/,//g;
	return $in;
}


sub _smpl_stats_columns {
	my @sdm = qw(totRds Rejected1 Rejected2 Accepted1 Accepted2 Singl1 Singl2 AvgSeqLen MaxSeqLength AvgSeqQual accErr);
	my @binners;
	for my $mode (1 .. 5) {
		my $name = getBinSubdirName($mode);
		push @binners, "HQ_bins_$name", "MQ_bins_$name", "${name}_other_bins", "${name}_total_bins";
	}
	return (
		# Input discovery and raw-upload preparation.
		qw(RawInputSize RawInputSizeSub InputIsPaired InputIsSingle
		FilteredContaRdsPerc_EBI FilteredContaRds_EBI FilteredNonContaRds_EBI),
		# SDM primary/support cleaning, followed by host filtering and FLASH.
		@sdm, qw(SDMAcceptedPercent), (map { "${_}_Sup" } @sdm), qw(SDMAcceptedPercent_Sup),
		qw(FilteredContaRdsPerc FilteredContaRds FilteredNonContaRds),
		qw(Merged NotMerged AvgGenomeSizeEst TotalGenomesEst),
		# Assembly and the optional hybrid comparison are produced together.
		qw(ContigN50 NScaff400 NScaffG1k NScaffG10k NScaffG100k NScaffG1M
		ScaffN50 ScaffMaxSize ScaffSize CircCtgs CircCtgG1M
		HybridPreassemblyCount
		HybridPreassemblyContigs HybridFinalContigs HybridContigsDelta HybridFinalToPreContigsRatio
		HybridPreassemblyBases HybridFinalBases HybridBasesDelta HybridFinalToPreBasesRatio
		HybridPreassemblyN50 HybridFinalN50 HybridN50Delta HybridFinalToPreN50Ratio
		HybridPreassemblyN90 HybridFinalN90 HybridN90Delta HybridFinalToPreN90Ratio
		HybridPreassemblyLongest HybridFinalLongest HybridLongestDelta HybridFinalToPreLongestRatio
		HybridPreassemblyGCPercent HybridFinalGCPercent HybridGCPercentDelta HybridFinalToPreGCPercentRatio),
		# Assembly mapping, duplicate handling, and breakpoint detection.
		qw(ReadsPaired AlignedReads OverallAlignment UniqueAlgned MultAlign DisconcAlign SingleUniqAlign SingleMultiAlign
		OpticalDuplicates PCRduplicates PassedMD EstLibSize
		BreakpointCount BreakpointContigs BreakpointBases BreakpointMeanLength BreakpointMaxLength AssemblyBreakpointPercent),
		# ContigStats gene summary.
		qw(GeneNumber AvgGeneLength AvgComplGeneLength BpGenes BpNotGenes GeneCodingPercent Gcomplete G5pComplete G3pComplete Gincomplete),
		# Binning is submitted after ContigStats.
		@binners,
		# Consensus SNP/INDEL calling follows binning submission.
		qw(SNP_TotalResolvedBp SNP_fastaEntries SNP_Num SNP_Passed SNP_resolved SNPsPerMbp INDEL_Num INDEL_Passed INDELsPerMbp),
	);
}

sub _metag_stats_text {
	my ($stats, $order) = @_;
	my @preferred = _smpl_stats_columns();
	my %known;
	my @duplicates;
	for my $column (@preferred) {
		push @duplicates, $column if $known{$column};
		$known{$column} = 1;
	}
	die "Duplicate field(s) in sample-stat schema: ".join(', ', @duplicates)."\n" if @duplicates;
	my %observed;
	for my $sample (@$order) {
		my $values = $stats->{$sample}{values} || {};
		for my $name (keys %$values) {
			$observed{$name} = 1 if defined($values->{$name}) && $values->{$name} ne '';
		}
	}
	my @unknown = sort grep { !$known{$_} } keys %observed;
	warn "Appending unknown sample-stat field(s): ".join(', ', @unknown)."\n" if @unknown;
	my @columns = (grep { $observed{$_} } @preferred, @unknown);
	my @lines = (join("\t", 'SMPLID', 'DIR', @columns));
	for my $sample (@$order) {
		next unless exists $stats->{$sample};
		my $values = $stats->{$sample}{values} || {};
		my @row = ($sample, $stats->{$sample}{DIR});
		for my $name (@columns) {
			my $value = defined($values->{$name}) ? $values->{$name} : '';
			$value =~ s/[\t\r\n]+/ /g;
			push @row, $value;
		}
		push @lines, join("\t", @row);
	}
	return join("\n", @lines)."\n";
}

sub _sdm_version {
	my ($text) = @_;
	return '' unless defined($text);
	return $1 if $text =~ m/\bsdm(?:\s+\([^\n)]*\))?\s+(\d+(?:\.\d+)+)\b/i;
	return '';
}

sub _sdm_version_at_least {
	my ($version, $required_major, $required_minor) = @_;
	return 0 unless defined($version) && $version =~ m/^(\d+)\.(\d+)/;
	return $1 > $required_major || ($1 == $required_major && $2 >= $required_minor);
}

sub _parse_sdm_stats_text {
	my ($filStats, $MaxLengthHistBased, $suffix) = @_;
	$filStats ||= '';
	$MaxLengthHistBased ||= 0;
	$suffix ||= '';
	my $sdmVersion = _sdm_version($filStats);
	my ($totRds,$Rejected1,$Rejected2,$Accepted1,$Accepted2,$Singl1,$Singl2,$AvgLen,$MaxLength,$AvgQual,$accErr) =
		("0","0","0","0","0","0","0","0","0","0","0");
	my $parsed = 0;

	# SDM 3.40 introduced a labelled, column-aligned summary and reports
	# high- and mid-quality accepted reads separately.  Accepted1 remains the
	# total accepted single-end reads in metagStats, so combine both classes.
	if (_sdm_version_at_least($sdmVersion, 3, 40)
			&& $filStats =~ m/^\s*Reads processed:\s*([0-9,]+)/m) {
		$parsed = 1;
		if ($filStats =~ m/^\s*Reads processed:\s*([0-9,]+)\s+([0-9,]+)\s*$/m) {
			# Paired 3.40+ summaries use aligned Read 1 / Read 2 columns.
			$totRds = remComma($1) + remComma($2);
			if ($filStats =~ m/^\s*Rejected:\s*([0-9,]+)\s+\([^)]*\)\s+([0-9,]+)/m) {
				$Rejected1 = remComma($1); $Rejected2 = remComma($2);
			}
			my ($high1, $high2, $mid1, $mid2) = (0, 0, 0, 0);
			if ($filStats =~ m/^\s*Accepted \(high quality\):\s*([0-9,]+)\s+\([^)]*\)\s+([0-9,]+)/mi) {
				$high1 = remComma($1); $high2 = remComma($2);
			}
			if ($filStats =~ m/^\s*Accepted \(mid quality\):\s*([0-9,]+)\s+\([^)]*\)\s+([0-9,]+)/mi) {
				$mid1 = remComma($1); $mid2 = remComma($2);
			}
			$Accepted1 = $high1 + $mid1; $Accepted2 = $high2 + $mid2;
			if ($filStats =~ m/^\s*Recovered singleton reads:\s*([0-9,]+)\s+([0-9,]+)\s*$/mi) {
				$Singl1 = remComma($1); $Singl2 = remComma($2);
			}
		} else {
			$filStats =~ m/^\s*Reads processed:\s*([0-9,]+)/m;
			$totRds = remComma($1);
			if ($filStats =~ m/^\s*Rejected:\s*([0-9,]+)/m) {
				$Rejected1 = remComma($1);
			}
			my ($acceptedHigh, $acceptedMid) = (0, 0);
			if ($filStats =~ m/^\s*Accepted \(high quality\):\s*([0-9,]+)/mi) {
				$acceptedHigh = remComma($1);
			}
			if ($filStats =~ m/^\s*Accepted \(mid quality\):\s*([0-9,]+)/mi) {
				$acceptedMid = remComma($1);
			}
			$Accepted1 = $acceptedHigh + $acceptedMid;
			$Singl1 = $Accepted1;
		}
		if ($filStats =~ m/^\s*-\s*[Ss]equence Length\s*:\s*[\d.]+\s*\/\s*([\d.]+)\s*\/\s*([\d.]+)/m) {
			$AvgLen = $1; $MaxLength = $2;
		}
		if ($filStats =~ m/^\s*-\s*Quality\s*:\s*[\d.]+\s*\/\s*([\d.]+)\s*\/\s*[\d.]+/m) {
			$AvgQual = $1;
		}
		if ($filStats =~ m/^\s*-\s*Accum\. Error\s+([\d.]+)/m) {
			$accErr = $1;
		}
	} elsif ($filStats =~ m/Reads processed: ([0-9,]+); ([0-9,]+) \(pa/) {
		$parsed = 1;
		$totRds = remComma($1) + remComma($2);
		if ($filStats =~ m/Rejected: ([0-9,]+); ([0-9,]+)\n/) {
			$Rejected1 = remComma($1); $Rejected2 = remComma($2);
		}
		if ($filStats =~ m/Accepted \(High qual\): ([0-9,]+); ([0-9,]+)\s/) {
			$Accepted1 = remComma($1); $Accepted2 = remComma($2);
		}
		if ($filStats =~ m/Singletons among these: ([0-9,]+); ([0-9,]+)\n/) {
			$Singl1 = remComma($1); $Singl2 = remComma($2);
		}
		if ($filStats =~ m/- [Ss]equence Length :\s*\d+.*\/([^\/]+)\/([^\/]+)\n/) {
			$AvgLen = $1; $MaxLength = $2;
		}
		if ($filStats =~ m/- Quality :\s*\d+.*\/(\d+.*)\/\d+.*\n/) {
			$AvgQual = $1;
		}
		if ($filStats =~ m/- Accum\. Error ([\d\.]+)/) {
			$accErr = $1;
		}
	} elsif ($filStats =~ m/Reads processed: ([0-9,]+)/) {
		$parsed = 1;
		$totRds = remComma($1);
		if ($filStats =~ m/Rejected: ([0-9,]+)\n/) {
			$Rejected1 = remComma($1);
		}
		if ($filStats =~ m/Accepted \(High qual\): ([0-9,]+)\s/) {
			$Accepted1 = remComma($1); $Singl1 = $Accepted1;
		}
		if ($filStats =~ m/- [Ss]equence Length :\s*\d+.*\/([^\/]+)\/([^\/]+)\n/) {
			$AvgLen = $1; $MaxLength = $2;
		}
		if ($filStats =~ m/.*- Quality :\s*\d+.*\/(\d+.*)\/\d+.*\n/) {
			$AvgQual = $1;
		}
		if ($filStats =~ m/.*- Accum\. Error ([\d\.]+)/) {
			$accErr = $1;
		}
	}
	if ($MaxLengthHistBased > $MaxLength) {$MaxLength = $MaxLengthHistBased;}
	my @names = qw(totRds Rejected1 Rejected2 Accepted1 Accepted2 Singl1 Singl2 AvgSeqLen MaxSeqLength AvgSeqQual accErr);
	my @values = ($totRds,$Rejected1,$Rejected2,$Accepted1,$Accepted2,$Singl1,$Singl2,$AvgLen,$MaxLength,$AvgQual,$accErr);
	my %result;
	my @output_values = $parsed ? @values : map { '' } @names;
	@result{map { "$_$suffix" } @names} = @output_values;
	return \%result;
}

sub _sdm_histogram_max_length {
	my ($inD) = @_;
	my $MaxLengthHistBased = 0;
	my $length_histogram = getFileStr("$inD/LOGandSUB/sdm/filter_lenHist.txt",0);
	if ($length_histogram ne '') {
		my @lines = split(/\n/, $length_histogram);
		$MaxLengthHistBased = $1 if @lines && $lines[-1] =~ m/^(\d+)\s/;
	}
	return $MaxLengthHistBased;
}

sub sdmStats {
	my ($inF,$inD,$suffix) = @_;
	my $MaxLengthHistBased = _sdm_histogram_max_length($inD);
	my $filStats = getFileStr($inF,0,70);
	return _parse_sdm_stats_text($filStats, $MaxLengthHistBased, $suffix);
}

sub sdmStatsMany {
	my ($files, $inD, $suffix) = @_;
	my @countFields = qw(totRds Rejected1 Rejected2 Accepted1 Accepted2 Singl1 Singl2);
	my @averageFields = qw(AvgSeqLen AvgSeqQual accErr);
	my %combined = map { $_ => 0 } (@countFields, @averageFields, 'MaxSeqLength');
	my %averageWeights;
	my $MaxLengthHistBased = _sdm_histogram_max_length($inD);
	foreach my $file (@{$files || []}) {
		next unless defined($file) && $file ne '';
		my $filStats = getFileStr($file,0,70);
		next if $filStats eq '';
		my $stats = _parse_sdm_stats_text($filStats, $MaxLengthHistBased, '');
		my $weight = ($stats->{totRds} || 0) + 0;
		$combined{$_} += ($stats->{$_} || 0) for @countFields;
		foreach my $field (@averageFields) {
			next unless defined($stats->{$field}) && $stats->{$field} ne '';
			$combined{$field} += $stats->{$field} * $weight;
			$averageWeights{$field} += $weight;
		}
		$combined{MaxSeqLength} = $stats->{MaxSeqLength}
			if (($stats->{MaxSeqLength} || 0) > $combined{MaxSeqLength});
	}
	$combined{$_} = $averageWeights{$_}
		? $combined{$_} / $averageWeights{$_} : '' for @averageFields;
	return {map { ($_.$suffix) => $combined{$_} }
		(@countFields, @averageFields, 'MaxSeqLength')};
}


sub bwtLogRd($$$){
	my ($splAr,$idx,$rhr) = @_;
	my @spl = @{$splAr};
	my %ret = %{$rhr};
	#DEFAULTS:
	$ret{uniqAlign}=-1;$ret{multAlign}=0;$ret{DisconcAlign}=0;
	return \%ret if $idx < 0 || $idx + 1 >= @spl;
	my $single_end = $spl[$idx+1] =~ m/\(100.00%\) were unpaired; of these:/ ? 1 : 0;
	my $last_required = $single_end ? $idx + 5 : $idx + 14;
	if ($last_required >= @spl) {
		warn "Incomplete bowtie2 statistics block; expected through line index $last_required\n";
		return \%ret;
	}
	if ($single_end){#single end read mapping!
		#die "SE\n";
		if ($spl[$idx+2] =~ m/(\d+) \(.+\) aligned 0 times/){ $ret{notAlign} = $1;} else { $ret{notAlign}=0; print "bwtOut wrg1.1\n";}
		 $ret{uniqAlign}=0;
		 $ret{multAlign}=0;
		$ret{DisconcAlign}=0;
		if ($spl[$idx+3] =~ m/(\d+) \(.+\) aligned exactly 1 time/ ){ $ret{SinglAlign} = $1;} else { $ret{SinglAlign}=0; warn "Could not parse single-read unique alignments\n";}
		if ($spl[$idx+4] =~ m/(\d+) \(.+\) aligned >1 times/ ){ $ret{SinglAlignMult} = $1;} else { $ret{SinglAlignMult}=0; warn "Could not parse single-read multiple alignments\n";}
		if ( $spl[$idx+5] =~ m/(\d+\.\d+)\% overall alignment rate/ ){ $ret{AlignmRate} = $1;} else {$ret{AlignmRate} = -1; warn "Could not parse overall alignment rate\n";}
	} else {
		#die "Pair\n";
		if ($spl[$idx+2] =~ m/(\d+) \(.+\) aligned concordantly 0 times/){ $ret{notAlign} = $1;} else { $ret{notAlign}=0; print "bwtOut wrg1\n";}
		if ($spl[$idx+3] =~ m/(\d+) \(.+\) aligned concordantly exactly 1 time/ ){ $ret{uniqAlign} = $1;} else { $ret{uniqAlign}=0;print "bwtOut wrg2\n";}
		if ($spl[$idx+4] =~ m/(\d+) \(.+\) aligned concordantly >1 times/ ){ $ret{multAlign} = $1;} else { $ret{multAlign}=0;print "bwtOut wrg3  $spl[6]\n";}
		if ($spl[$idx+7] =~ m/(\d+) \(.+\) aligned discordantly 1 time/ ){ $ret{DisconcAlign} = $1;} else { $ret{DisconcAlign}=0;print "bwtOut wrg4 $spl[9]\n";}
		if ($spl[$idx+12] =~ m/(\d+) \(.+\) aligned exactly 1 time/ ){ $ret{SinglAlign} = $1;} else { $ret{SinglAlign}=0;print "bwtOut wrg5 $spl[14]\n";}
		if ($spl[$idx+13] =~ m/(\d+) \(.+\) aligned >1 times/ ){ $ret{SinglAlignMult} = $1;} else { $ret{SinglAlignMult}=0; warn "Could not parse paired-read multiple alignments\n";}
		if ($spl[$idx+14] =~ m/(\d+\.\d+)\% overall alignment rate/ ){ $ret{AlignmRate} = $1;} else {$ret{AlignmRate} = -1; warn "Could not parse overall alignment rate\n";}
	}
	return (\%ret);
}


sub getMapStats{
	my ($inP) = @_;
	my $inFi = "$inP/map.sh.etxt"; 
	$inFi = "$inP/bwtMap.sh.etxt" if (!-e $inFi && -e "$inP/bwtMap.sh.etxt"); #old MF file names..
	#my $outStrDesc = "";
	my @spl = ();
	my $alignStats = getFileStr($inFi,0);
	if (-s $inFi){
		@spl = split(/\n/,$alignStats);
	}
	#die "$alignStats\n";
	my $idx =-1;my $dobwtStat=0;
	if ($alignStats =~ m/This is strobealign/){
		$dobwtStat=3;
	}elsif($alignStats =~ m/\[M::worker_pipeline/){
		$dobwtStat=2;
	}elsif (@spl > 12 && $alignStats =~ m/reads; of these:/){
		$dobwtStat=1;
		for my $candidate (0 .. $#spl) {
			if ($spl[$candidate] =~ m/\d+ reads; of these:/) { $idx = $candidate; last; }
		}
		if ($idx>=0 && $spl[$idx+0] =~ m/(\d+) reads; of these:/){
			$locStats{totReadPairs} = $1;
		} else {
			$dobwtStat=0;#die "wrong bwtOut: $inD/LOGandSUB/bwtMap.sh.etxt X $spl[2]";
		}
	}
	#die "$spl[$idx]\n$locStats{totReadPairs}\n";
	my @columns = qw(ReadsPaired AlignedReads OverallAlignment UniqueAlgned MultAlign DisconcAlign SingleUniqAlign SingleMultiAlign);
	my %result = map { $_ => '' } @columns;
	
	my $incoming =0;my $retained =0;my $removed = 0;
	if ($alignStats =~ m/^Inentries: (\d+)/m){
		my @matc = $alignStats =~ m/^Inentries: (\d+)/mg;		foreach (@matc) {$incoming += int($_);}
		 @matc =$alignStats =~ m/^TotalRetained: (\d+)/mg;	foreach (@matc) {$retained += int($_);}
		 @matc = $alignStats =~ m/^TotalRm: (\d+)/mg;	foreach (@matc) {$removed += int($_);}
		 $locStats{totReadPairs} = $incoming;
		 $locStats{uniqAlign} = $retained;
	} else {
		$locStats{totReadPairs} = -1;
		$locStats{uniqAlign} = -1;
	}
	
	if (!$dobwtStat){
		#$outStr .= "\t" x 7;
		#die "X!\n";
	} elsif ($dobwtStat == 1){#$alignStats =~ m/reads; of these:/){
		my ($rhr) = bwtLogRd(\@spl,$idx,\%locStats);
		%locStats = %{$rhr};
		$result{ReadsPaired} = $locStats{totReadPairs};
		$result{AlignedReads} = $retained;
		$result{OverallAlignment} = $locStats{AlignmRate};
		if (($locStats{totReadPairs} || 0) > 0) {
			$result{UniqueAlgned} = $locStats{uniqAlign}/$locStats{totReadPairs}*100;
			$result{MultAlign} = $locStats{multAlign}/$locStats{totReadPairs}*100;
			$result{DisconcAlign} = $locStats{DisconcAlign}/$locStats{totReadPairs}*100;
			$result{SingleUniqAlign} = $locStats{SinglAlign}/$locStats{totReadPairs}*100;
			$result{SingleMultiAlign} = $locStats{SinglAlignMult}/$locStats{totReadPairs}*100;
		}
	} elsif ($dobwtStat == 2){ #minimap2
		my @matches = ($alignStats =~ m/\[M::worker_pipeline::.*\] mapped (\d+) sequences/g);
		#print "@matches\n";
		my $sum=0; $sum += $_ foreach (@matches);
		if ($sum>0){
			my $frac = (1 - ($removed/$sum)) * 100;
			@result{qw(ReadsPaired AlignedReads OverallAlignment)} = ($sum, $retained, $frac);
		}
		#die "minimap!!$sum\n";
	} elsif ($dobwtStat == 3){ #strobealign
		if ($incoming > 0){
			#$locStats{totReadPairs} = -1 if (!exists($locStats{totReadPairs}));
			my $frac = 0; $frac = $incoming/$locStats{totReadPairs} if( $locStats{totReadPairs}>0);
			$locStats{AlignmRate}=$frac;
			@result{qw(ReadsPaired AlignedReads OverallAlignment UniqueAlgned)} =
				($locStats{totReadPairs}, $retained, $frac,
				 $locStats{uniqAlign}/$locStats{totReadPairs}*100);
		}
	}
	return \%result;
}
sub optiDups{
	my ($inP) = @_;
	my $inFi = "$inP/map2.sh.etxt"; 
	$inFi = "$inP/bwtMap2.sh.etxt" if (!-e $inFi); #old MF file names..
	#my $outStrDesc = "";
	$locStats{duplOptic}=0; $locStats{duplPCR}=0;$locStats{duplPass}=0; $locStats{EstLibSize} = 0;
	my $doDup = 0;my $alignStats2 = getFileStr($inFi,0);
	if ($alignStats2 ne "" ){if (defined ($alignStats2)){$doDup=1;}}
	if ($doDup){
		
		if ($alignStats2 =~ m/samtools markdup/){#samtools dedup
			my $pairDup =0; my $singlDup=0; my $optiPairDup=0; my $optiSinglDup=0;
			my $Npairs=0; my $Nsingle=0;
			if ($alignStats2 =~ m/PAIRED: (\d+)/){$Npairs=$1;}
			if ($alignStats2 =~ m/SINGLE: (\d+)/){$Nsingle=$1;}
			if ($alignStats2 =~ m/WRITTEN: (\d+)/){$locStats{duplPass}=$1;}
			if ($alignStats2 =~ m/DUPLICATE PAIR: (\d+)/){$pairDup=$1;}
			if ($alignStats2 =~ m/DUPLICATE SINGLE: (\d+)/){$singlDup=$1;}
			if ($alignStats2 =~ m/DUPLICATE PAIR OPTICAL: (\d+)/){$optiPairDup=$1;}
			if ($alignStats2 =~ m/DUPLICATE SINGLE OPTICAL: (\d+)/){$optiSinglDup=$1;}
			if ($alignStats2 =~ m/ESTIMATED_LIBRARY_SIZE: (\d+)/){$locStats{EstLibSize}=$1;}
			$locStats{duplPCR} = ($pairDup+$singlDup);#/($Npairs+$Nsingle);
			$locStats{duplOptic} = ($optiPairDup+$optiSinglDup);#/($Npairs+$Nsingle);
			
		} else {
			#print "$alignStats2 YUYS\n";
			if ($alignStats2 =~ m/Proper Pairs\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/){$locStats{EstLibSize} = $4; $locStats{duplOptic} = $3; $locStats{duplPCR} = $2; $locStats{duplPass} = $1;} else {$doDup=0;}#die "Can't find dupl stats: $alignStats2\n";}
			if ($alignStats2 =~ m/Improper Pairs\s+(\d+)\s+(\d+)\s+(\d+)/){$locStats{duplOptic} += $3; $locStats{duplPCR} += $2; $locStats{duplPass} += $1;} else {$doDup=0;}#{die "Can't find impr. dupl stats: $alignStats2\n";}
		}
	}
	my %result = map { $_ => '' } qw(OpticalDuplicates PCRduplicates PassedMD EstLibSize);
	@result{qw(OpticalDuplicates PCRduplicates PassedMD EstLibSize)} =
		@locStats{qw(duplOptic duplPCR duplPass EstLibSize)} if $doDup;
	return \%result;
}

sub getContamination{
	my ($inFi,$inFi2,$idx) = @_;
	my $suffix = $idx eq '' ? '' : "_$idx";
	my %result = (
		"FilteredContaRdsPerc$suffix" => '',
		"FilteredContaRds$suffix" => '',
		"FilteredNonContaRds$suffix" => '',
	);
	my $filStats = getFileStr($inFi,0);#`cat $inD/LOGandSUB/KrakHS.sh.etxt`; chomp $filStats;
	#if ($filStats eq "" ){$outStr .= "?\t?\t?\t";	return ($outStr,$outStrDesc);}
	
	my @matches = ($filStats =~ m/\d+ sequences classified \((\d+\.?\d*)%\)/g);
	my @hits = ($filStats =~ m/(\d+) sequences classified \(\d+\.?\d*%\)/g);
	my @nonhits = ($filStats =~ m/(\d+) sequences unclassified \(\d+\.?\d*%\)/g);
	
	if (@nonhits == 0){#check if this was done via hostile..
		$filStats = getFileStr($inFi2,0); #paired reads should be counted as two.. but are counted as one in hostile..
		@hits = ($filStats =~ m/"reads_removed": (\d*),/g); #"reads_removed": 202,
		@nonhits = ($filStats =~ m/"reads_out": (\d*),/g); 
		@matches = ($filStats =~ m/"reads_removed_proportion": (\d*),/g); 
	}
	
	
	if (@nonhits > 0){
		my $totHits=0; my $totNH=0;
		my $pair_count = @hits < @nonhits ? scalar(@hits) : scalar(@nonhits);
		for (my $i=0;$i<$pair_count;$i++){
			$totHits += $hits[$i];$totNH += $nonhits[$i];
		}
		$result{"FilteredContaRdsPerc$suffix"} = join(";",@matches);
		$result{"FilteredContaRds$suffix"} = $totHits;
		$result{"FilteredNonContaRds$suffix"} = $totNH;
	} 
	return \%result;
}

sub getBinnerStats{
	my ($tmpassD,$SmplN) = @_;
	my %result;
	for my $i (1 .. 5) {
		my $name = getBinSubdirName($i);
		@result{"HQ_bins_$name", "MQ_bins_$name", "${name}_other_bins", "${name}_total_bins"} = ('', '', '', '');
	}
	return \%result unless defined($tmpassD) && $tmpassD ne '' && -d $tmpassD;
	#all binners..
	for (my $i=1; $i < 6; $i++){
		my $SCdir = getBinSubdirName($i);
		my $SBbinCM2 = "$tmpassD/Binning/$SCdir/$SmplN.cm2";
		if (-s $SBbinCM2){
			my $HQbinCnt = 0; my $MQbinCnt = 0;my $totBins=0;
			open my $bin_fh, '<', $SBbinCM2 or do { warn "Cannot read bin statistics '$SBbinCM2': $!\n"; next; };
			while (<$bin_fh>){my @spl = split/\t/; next unless @spl >= 3; next if ($spl[1] eq "Completeness");
				if ($spl[1] >= 90 && $spl[2] <= 5){$HQbinCnt ++ ;
				} elsif ($spl[1] >= 80 && $spl[2] <= 5){$MQbinCnt ++ ;}
				$totBins ++;
			}
			close $bin_fh or warn "Cannot close bin statistics '$SBbinCM2': $!\n";
			@result{"HQ_bins_$SCdir", "MQ_bins_$SCdir", "${SCdir}_other_bins", "${SCdir}_total_bins"} =
				($HQbinCnt, $MQbinCnt, $totBins-$HQbinCnt-$MQbinCnt, $totBins);
		} else {
			@result{"HQ_bins_$SCdir", "MQ_bins_$SCdir", "${SCdir}_other_bins", "${SCdir}_total_bins"} = ('', '', '', '');
		}
	}
	return \%result;
}

sub getSNPStats{
	my ($inFi) = @_;
	my $geneStats = getFileStr("${inFi}.etxt",0);
	my @columns = qw(SNP_TotalResolvedBp SNP_fastaEntries SNP_Num SNP_Passed SNP_resolved SNPsPerMbp INDEL_Num INDEL_Passed INDELsPerMbp);
	my %empty = map { $_ => '' } @columns;
	my $has_stats = 0;
	#Total bp written: 34950480 (0 not resolved) on 26727 entries
	my $bps = 0; my $entrs=0;my $confl=0; my $resol=0;my $snpNum=0; my $indelNuml=0;
	my $SNPpassed=0; my $INDpassed=0;
	if ($geneStats =~ m/Total bp written: (\d+) \(\d+ not resolved\) on (\d+) entries/){
		$has_stats = 1;
		$bps = $1; $entrs=$2;
		if ($geneStats =~ m/Conflicting calls: (\d+) Resolved with second line: (\d+)/) {
			$confl=$1; $resol=$2;
		}
		$snpNum=$1 if $geneStats =~ m/Total SNPs detected: (\d+)/;
		#$outStr = "$bps\t$entrs\t$confl\t$resol\t$snpNum\t";
		
	} else { #try in .otxt
		$geneStats = getFileStr("${inFi}.otxt",0);
		if ($geneStats =~ m/Total bp that can be determined: (\d+) in (\d+) entries./){
			$has_stats = 1;
			$bps = $1; $entrs=$2;
			#  - Found 11 SNPs and 5 INDELS.
			if ($geneStats =~ m/  - Found (\d+) SNPs and (\d+) INDELS./) {
				$snpNum=$1; $indelNuml=$2;
			}
			if ($geneStats =~ m/  - Passed (\d+);(\d+) SNPs and INDELS. Conflicts resolved: 0/) {
				$SNPpassed=$1; $INDpassed=$2;
			}
		}
	}
	return \%empty unless $has_stats;
	my $per_mbp = $bps > 0 ? 1_000_000 / $bps : 0;
	return { SNP_TotalResolvedBp=>$bps, SNP_fastaEntries=>$entrs, SNP_Num=>$snpNum,
		SNP_Passed=>$SNPpassed, SNP_resolved=>$resol,
		SNPsPerMbp=>sprintf('%.3f', $snpNum * $per_mbp), INDEL_Num=>$indelNuml,
		INDEL_Passed=>$INDpassed, INDELsPerMbp=>sprintf('%.3f', $indelNuml * $per_mbp) };
}

sub getGeneStats{
	my ($inFi) = @_;
	my @columns = qw(GeneNumber AvgGeneLength AvgComplGeneLength BpGenes BpNotGenes Gcomplete G5pComplete G3pComplete Gincomplete);
	my %result = map { $_ => '' } @columns;
	my $geneStats = getFileStr("$inFi",0);
	my @spl1 = split("\n", $geneStats);

	#if (!-e "$inFi" || -s "$inFi" == 0){	
	#	return ($outStr, $outStrDesc);
	#}
	#open I,"<$inFi" or die $!; 
	#GeneNumber	AvgGeneLength	AvgComplGeneLength	BpGenes	BpNotGenesGcomplete	G5pComplete	G3pComplete	Gincomplete
	#my $tmp ="";$tmp = <I>;#if (defined($tmp)){chomp($tmp);$outStrDesc .= $tmp."\t";}$tmp = <I>;
	#while (my $tmp = <I>){
	foreach my $tmp (@spl1){
		chomp($tmp);#$outStr = $tmp."\t"; $outStr5 .= $tmp."\t" if ($do500Stat); 
		next if ($tmp =~ m/^GeneNumber/);
		my @spl = split /\t/,$tmp; 
		#print @spl . " @spl\n";
		if (@spl >=9){
			@result{@columns} = @spl[0..8];
		}
	}
	return \%result;
}
sub getASsemblyStats{
	my ($tmpassD,$assemblStatsFile,$doCirc) = @_;
	my @columns = qw(ContigN50 NScaff400 NScaffG1k NScaffG10k NScaffG100k NScaffG1M ScaffN50 ScaffMaxSize ScaffSize CircCtgs CircCtgG1M);
	my %result = map { $_ => '' } @columns;
	return \%result unless defined($tmpassD) && $tmpassD ne '';
	$tmpassD =~ s{/$}{};
	my $stats_path = "$tmpassD/$assemblStatsFile";
	my $assStats = getFileStr($stats_path,0);
	if ($assStats ne ""){
		my %patterns = (
			ContigN50 => qr/N50 contig length\s+(\d+)/,
			NScaff400 => qr/Number of scaffolds\s+(\d+)/,
			ScaffSize => qr/Total size of scaffolds\s+(\d+)/,
			ScaffMaxSize => qr/Longest scaffold\s+(\d+)/,
			ScaffN50 => qr/N50 scaffold length\s+(\d+)/,
			NScaffG1k => qr/Number of scaffolds > 1K nt\s+(\d+)/,
			NScaffG10k => qr/Number of scaffolds > 10K nt\s+(\d+)/,
			NScaffG100k => qr/Number of scaffolds > 100K nt\s+(\d+)/,
			NScaffG1M => qr/Number of scaffolds > 1M nt\s+(\d+)/,
		);
		my @missing;
		for my $name (keys %patterns) {
			if ($assStats =~ $patterns{$name}) { $result{$name} = $1; }
			else { push @missing, $name; }
		}
		warn "Incomplete assembly statistics '$stats_path'; missing: ".join(', ', sort @missing)."\n" if @missing;
		$locStats{CtgN50} = $result{ContigN50} if $result{ContigN50} ne '';
		$locStats{NScaff} = $result{NScaff400} if $result{NScaff400} ne '';
		$locStats{ScaffSize} = $result{ScaffSize} if $result{ScaffSize} ne '';
	}
	if ($assStats ne '') { $result{CircCtgs} = 0; $result{CircCtgG1M} = 0; }
	if ($doCirc && -e "$tmpassD/scaffolds.fasta.circ"){
		$result{CircCtgs} = 0; $result{CircCtgG1M} = 0;
		$assStats = getFileStr("$tmpassD/scaffolds.fasta.circ",0);
		$result{CircCtgs} = $assStats =~ tr/>//;
		my @matchs = ($assStats =~ m/.*_L=(\d+)=/g);
		foreach(@matchs){$result{CircCtgG1M}++ if ($_ > 1000000);}
	}
	return \%result;
}

sub getBreakpointStats {
	my ($path) = @_;
	my @columns = qw(BreakpointCount BreakpointContigs BreakpointBases BreakpointMeanLength BreakpointMaxLength);
	my %result = map { $_ => '' } @columns;
	return \%result unless fileGZe($path);
	my $fh;
	eval { ($fh) = gzipopen($path, 'breakpoint statistics', 1); 1 }
		or do { warn "Cannot read breakpoint report '$path': $@"; return \%result; };
	my ($count, $bases, $maximum) = (0, 0, 0);
	my (%contigs, $malformed);
	while (my $line = <$fh>) {
		$line =~ s/[\r\n]+$//;
		next if $line eq '' || $line =~ /^(?:#|contig\tstart\tend(?:\t|$))/;
		my @fields = split /\t/, $line;
		unless (@fields >= 3 && $fields[1] =~ /^\d+$/
			&& $fields[2] =~ /^\d+$/ && $fields[2] > $fields[1]) { $malformed++; next; }
		my $length = (@fields >= 4 && $fields[3] =~ /^\d+$/)
			? $fields[3] : $fields[2] - $fields[1];
		$count++; $bases += $length; $maximum = $length if $length > $maximum;
		$contigs{$fields[0]} = 1;
	}
	close $fh or warn "Cannot close breakpoint report '$path': $!\n";
	warn "Skipped $malformed malformed breakpoint row(s) in '$path'\n" if $malformed;
	my $mean = $count ? sprintf('%.1f', $bases / $count) : 0;
	@result{@columns} = ($count, scalar(keys %contigs), $bases, $mean, $maximum);
	return \%result;
}

sub getHybridAssemblyStats {
	my ($path) = @_;
	my %result;
	return \%result unless defined($path) && -s $path;
	open my $fh, '<', $path or do { warn "Cannot read hybrid assembly report '$path': $!\n"; return \%result; };
	my %prefix = (
		contigs => 'Contigs', total_bp => 'Bases', N50 => 'N50', N90 => 'N90',
		longest => 'Longest', GC_percent => 'GCPercent',
	);
	my $malformed = 0;
	while (my $line = <$fh>) {
		$line =~ s/[\r\n]+$//;
		next if $line eq '' || $line =~ /^metric\t/;
		my @fields = split /\t/, $line, -1;
		if ($fields[0] eq 'source_preassembly_count') {
			$result{HybridPreassemblyCount} = $fields[1] if defined($fields[1]) && $fields[1] ne '';
			next;
		}
		my $name = $prefix{$fields[0]};
		unless (defined($name) && @fields >= 5) { $malformed++; next; }
		@result{"HybridPreassembly$name", "HybridFinal$name", "Hybrid${name}Delta", "HybridFinalToPre${name}Ratio"}
			= @fields[1..4];
	}
	close $fh or warn "Cannot close hybrid assembly report '$path': $!\n";
	warn "Skipped $malformed malformed hybrid comparison row(s) in '$path'\n" if $malformed;
	return \%result;
}

sub smplStats {
	my ($inD,$assDir,$SmplN) = @_;
	my %values;
	my $merge = sub {
		my ($named) = @_;
		return unless ref($named) eq 'HASH';
		@values{keys %$named} = values %$named;
	};
	my $rawReadSet = sampleReadSet($SmplN, "raw");
	my $seq_set = ref($rawReadSet) eq 'HASH' ? $rawReadSet : {};
	my $input_libraries = readLibrariesByScope($seq_set, 'primary', 0, $SmplN);
	$values{InputIsPaired} = @{libraryPairs($input_libraries)} ? 1 : 0;
	$values{InputIsSingle} = @{libraryFiles($input_libraries, 'single')} ? 1 : 0;
	$values{RawInputSize} = exists($map{$SmplN}{inputFileSizeMB})
		? sprintf('%.3fG', $map{$SmplN}{inputFileSizeMB}/1024) : -1;
	$values{RawInputSizeSub} = exists($map{$SmplN}{inputXFileSizeMB})
		? sprintf('%.3fG', $map{$SmplN}{inputXFileSizeMB}/1024) : -1;

	# Raw-upload preparation is submitted before read cleaning.
	$merge->(getContamination("$inD/LOGandSUB/prepEBI.sh.etxt", "$inD/LOGandSUB/prepEBI.sh.otxt", 'EBI'));

	my @primary_logs = grep { $_ !~ /filterSuppl/ } glob("$inD/LOGandSUB/sdm/filter*.log");
	@primary_logs = ("$inD/LOGandSUB/sdmReadCleaner.sh.etxt")
		if (!@primary_logs && -s "$inD/LOGandSUB/sdmReadCleaner.sh.etxt");
	if (@primary_logs) {
		my $stats = sdmStatsMany(\@primary_logs, $inD, '');
		$merge->($stats);
		my $accepted = ($stats->{Accepted1} || 0) + ($stats->{Accepted2} || 0);
		$values{SDMAcceptedPercent} = sprintf('%.3f', 100 * $accepted / $stats->{totRds})
			if (($stats->{totRds} || 0) > 0);
		$locStats{$_} = $stats->{$_} for qw(totRds Rejected1 Rejected2 Accepted1 Accepted2 Singl1 Singl2);
	}
	my @support_logs = glob("$inD/LOGandSUB/sdm/filterSuppl*.log");
	@support_logs = ("$inD/LOGandSUB/sdmReadCleanerSuppl.sh.etxt")
		if (!@support_logs && -s "$inD/LOGandSUB/sdmReadCleanerSuppl.sh.etxt");
	if (@support_logs) {
		my $stats = sdmStatsMany(\@support_logs, $inD, '_Sup');
		$merge->($stats);
		my $accepted = ($stats->{Accepted1_Sup} || 0) + ($stats->{Accepted2_Sup} || 0);
		$values{SDMAcceptedPercent_Sup} = sprintf('%.3f', 100 * $accepted / $stats->{totRds_Sup})
			if (($stats->{totRds_Sup} || 0) > 0);
	}

	my $contamination = getContamination("$inD/LOGandSUB/KrakHS.sh.etxt", "$inD/LOGandSUB/KrakHS.sh.otxt", '');
	$merge->($contamination);
	$locStats{contamination} = $contamination->{FilteredContaRdsPerc};

	my $text = getFileStr("$inD/LOGandSUB/flashMrg.sh.otxt",0);
	my @mergedCounts = $text =~ /\[FLASH\]\s+Combined pairs:\s+(\d+)/g;
	my @unmergedCounts = $text =~ /\[FLASH\]\s+Uncombined pairs:\s+(\d+)/g;
	$values{Merged} = sum(@mergedCounts) if @mergedCounts;
	$values{NotMerged} = sum(@unmergedCounts) if @unmergedCounts;
	$text = getFileStr("$inD/MicroCens/MC.0.result",0);
	$values{AvgGenomeSizeEst} = $1 if $text =~ /average_genome_size:\s*([\d.]+)/;
	$values{TotalGenomesEst} = $1 if $text =~ /genome_equivalents:\s*([\d.]+)/;

	my $assembly_dir = '';
	if (-e "$inD/assemblies/metag/assembly.txt") {
		my $pointer = getFileStr("$inD/assemblies/metag/assembly.txt",0);
		($assembly_dir) = grep { $_ ne '' } map { s/^\s+|\s+$//gr } split /[\r\n]+/, $pointer;
		$assembly_dir ||= '';
	} elsif (-e "$inD/assemblies/metag/AssemblyStats.txt") {
		$assembly_dir = "$inD/assemblies/metag";
	} elsif (-e "$assDir/metag/AssemblyStats.txt") {
		$assembly_dir = "$assDir/metag";
	}
	if ($assembly_dir ne '' && !-d $assembly_dir) {
		warn "Assembly-statistics pointer is not a directory: '$assembly_dir'\n";
		$assembly_dir = '';
	}
	$merge->(getASsemblyStats($assembly_dir, 'AssemblyStats.txt', 1));
	$merge->(getHybridAssemblyStats("$assembly_dir/HybridAssemblyComparison.tsv")) if $assembly_dir ne '';
	$merge->(getMapStats("$inD/LOGandSUB/"));
	$merge->(optiDups("$inD/LOGandSUB/"));
	$merge->(getBreakpointStats("$inD/mapping/$SmplN-smd.bam.breakpoints.tsv.gz"));
	$values{AssemblyBreakpointPercent} = sprintf('%.6f', 100 * $values{BreakpointBases} / $values{ScaffSize})
		if (($values{ScaffSize} || 0) > 0 && defined($values{BreakpointBases}) && $values{BreakpointBases} ne '');
	$merge->(getGeneStats("$inD/$preDIRs{dir_ContigStats}/GeneStats.txt"));
	my $annotated_bases = ($values{BpGenes} || 0) + ($values{BpNotGenes} || 0);
	$values{GeneCodingPercent} = sprintf('%.3f', 100 * $values{BpGenes} / $annotated_bases)
		if ($annotated_bases > 0);
	$merge->(getBinnerStats($assembly_dir,$SmplN));
	$merge->(getSNPStats("$inD/LOGandSUB/SNP/ConsAssem.oSNPc.sh"));
	return \%values;
}

# smplStats is implemented above using the fixed named-value schema.

sub spadesHosts{
	#figure out if only certain node subset has enough HDD space
	if (0 && $HDDspace{spades} > 100){
		my $locHosts = `bhosts | grep ok | cut -d" " -f 1 | grep compute | tr "\\n" ","`;
		my $tmpStr = `pdsh -w $locHosts -u 4 "df -l" `;
		my $srchTerm = '/$';
		if (`hostname` =~ m/submaster$/){
			$srchTerm = '/tmp';
		}
		#print "psdh done\n";
		#die $tmpStr;
		#54325072 = 52G
		foreach my $l (split(/\n/,$tmpStr)){
			next if ($l !~ m/compute\S+:/);
			next unless ($l =~ m/$srchTerm/);
			#print $l."\n";
			my @hosts = split(/:/,$l);
			my @spl = split(/\s+/,$hosts[1]);
			#die $spl[1]." XX ".$spl[2]." XX ".$spl[3]." XX ".$spl[4]." \n ";
			#print $spl[4] / 1024 / 1024 . "\n";
			if (($spl[4] / 1024 / 1024) > $HDDspace{spades}){
				push(@{$QSBoptHR->{Spades_Hosts}},$hosts[0]);
			}
			if (($spl[4] / 1024 / 1024) > 40){
				push(@{$QSBoptHR->{General_Hosts}},$hosts[0]);
			}
		}
		print "Found ".scalar @{$QSBoptHR->{Spades_Hosts}}." host machines with > $HDDspace{spades} G space\n";
		print "Found ".scalar @{$QSBoptHR->{General_Hosts}}." host machines with > 40 G space\n";
		if (scalar @{$QSBoptHR->{Spades_Hosts}} ==0 ){die "Not enough hosts for spades temp space found\n";}
		sleep (2);
	}
}


sub readG2M($){
	my ($inf) = @_;
	open K,"<",$inf;
	my %COGs2motus;
	my %motus;
	while (my $li = <K>){
		chomp($li);
		my @spl = split("\t",$li);
		$motus{$spl[0]} = $spl[3];
		$COGs2motus{$spl[0]} = $spl[2];
	}
	close K;
	print "Read gene to map file\n";
}

sub scndMap2Genos{
	# ----------------- 2nd mapping (map to ref genomes supplied by user) ---------------------
	my ($SmplName,$cMapGrp,$cAssGrp,$curOutDir,$nodeSpTmpD,
			$sampleDepsAR,$calc2ndMapSNP,$boolScndCoverageOK) = @_;
	my @mapOutXS;my @bamBaseNameS;
	my @pendingSecondMapSNP;
	my $cleanSeqSetHR = sampleReadSet($curSmpl, "clean");
	my $samplReadLength = sampleReadSet($curSmpl, "raw")->{samplReadLength};

	for (my $i=0;$i<@bwt2outD; $i++){
		push @mapOutXS, $bwt2outD[$i];
		my $fname = $bwt2ndMapNmds[$i]."_".$SmplName."-0";
		#if (length($fname ) > 20){$bwt2ndMapNmds[$i]."_".$SmplName."-0";
		push @bamBaseNameS, $fname;
	}
#		system "rm -rf ".join(" ",@mapOutXS) if ($redo2ndMapping);

	$make2ndMapDecoy{Lib} = $curOutDir if ($MFopt{DoMapModeDecoy});
	my $cramthebam=0;
	my $secondMapLibraries = getRawLibrariesAssmGrp(\%AsGrps,$cAssGrp,0,$SmplName);
	my $secondMapTechnology = libraryTechnology(
		$secondMapLibraries, "secondary mapping reads for $SmplName", 1,
	);
	#map to all refs at once		
	my %dirset = 	(nodeTmp=>$nodeSpTmpD,outDir => join(",",@bwt2outD),unalDir=>"",
					sbj => join(",",@DBbtRefX),assGrp => $cAssGrp,
					smplName => $SmplName,#join(",",@bamBaseNameS),
					glbTmp => $nodeSpTmpD."_xtraMapWork/", is2ndMap => 1,
					qsubDir => "$logDir/map2nd/",mapSupport => 0,
					glbMapDir => join(",",@mapOutXS),mappingStarted =>1,
					libraries => $secondMapLibraries,
					readTec => '',
					submit => 1, submNow => 1, cramAlig => $cramthebam,
					sortCores => $MFopt{bamSortCores}, mapCores => $MFopt{MapperCores},
					deferMappingCleanup => 1);
	#die "XX @DBbtRefX\n";
	# ----------------- 2nd mapping (map to ref genomes supplied by user) ---------------------
	my ($map2CtgsX, undef, $mapOptHr) =
		mapReadsToRef(\%dirset,$AsGrps{$cMapGrp}{SeqUnZDeps}.";".$bwt2ndMapDep );#$localAssembly);
	#------------  and calc coverage for each separate
	#different strategy, a bit hacky: deactivate submissions and collect for one call
	$dirset{submit} = 0; my $bigMap = "";my $bigCov = "";
	for (my $i=0;$i<@bwt2outD; $i++){
		 my $mapOutX = $mapOutXS[$i]; 
		$dirset{outDir} = $bwt2outD[$i];$dirset{glbMapDir} =$mapOutXS[$i];
		$dirset{smplName} = $bamBaseNameS[$i];$dirset{sbj} =$DBbtRefX[$i];
		my ($map2CtgsY,$delaySubmCmdY,$mapStat)  = bamDepth(\%dirset,$map2CtgsX,$mapOptHr);
		$bigMap .= "\n\n#----------$i -------------\n$delaySubmCmdY\n" if ($delaySubmCmdY ne "");
		$delaySubmCmdY = "";
		
		#die "$map2CtgsY\n$delaySubmCmdY\n$mapStat\n";
		# Derive per-reference abundance from the canonical coverage. When the
		# mapping is new this command follows its publisher in the same job.
		if ($MFopt{mapModeCovDo} && !$boolScndCoverageOK){
			($map2CtgsY,$delaySubmCmdY) = calcCoverage2nd(
				"$mapOutX/$bamBaseNameS[$i]-smd.bam.coverage.gz", $DBbtRefGFF[$i],
				$samplReadLength, $SmplName.".$i",
				$map2CtgsY.";".$bwt2ndMapDep, \%dirset);
			#die "ASD\n";
		}
		# Per-reference jobs are collected into the combined sort/coverage submission below.
		$bigCov .= "\n\n#---------- $i -------------\n$delaySubmCmdY\n" if ($delaySubmCmdY ne "");
		
		
		
		if ($calc2ndMapSNP){ #2nd mapping SNP calling (consensus) (not assembly!!)
			my %SNPinfo = (
				assembly => "$bwt2outD[$i]/$bwt2ndMapNmds[$i].fa",#$DBbtRefX[$i],
				MAR => ["$bwt2outD[$i]/$bamBaseNameS[$i]-smd.bam"],
				SNPcaller => $MFopt{SNPcallerFlag},bamcram=>"bam",normIndels => $MFopt{normSNPindels},
				#doesn't work, if contig name and length is not given..
				#depthF => "$bwt2outD[$i]/$bamBaseNameS[$i]-smd.bam.coverage.gz.percontig",
				ofas => "$bwt2outD[$i]/$bamBaseNameS[$i].SNPc.$MFopt{SNPcallerFlag}.fna", #only output needed for this.. unless I later want to add also a gene calling.. (not needed for TEC2 reb)
				#vcfFile => "$bwt2outD[$i]/$bamBaseNameS[$i].$MFopt{SNPcallerFlag}.vcf",
				firstInSample => ($i == 0 ? 1 : 0), 
				SeqTech => $secondMapTechnology, SeqTechSuppl => "",
				nodeTmpD => $dirset{nodeTmp}, #$nodeSpTmpD,
				scratch => $dirset{glbTmp},
				qsubDir => $dirset{qsubDir}, jdeps => $map2CtgsY.";$bwt2ndMapDep",
				cmdFileTag => $bwt2ndMapNmds[$i], minDepth => $MFopt{consSNPminDepth},
				callSVs => $MFopt{callSVs}, vcfSVfile => "$bwt2outD[$i]/$bamBaseNameS[$i]-smd.SV.vcf", vcfSVfileS => "", callSVsSupp => 0,
				smpl => $bamBaseNameS[$i], maxCores => $MFopt{maxSNPcores}, #memReq => $MFopt{memSNPcall},
				bpSplit => 4e5,	runLocal => 1, split_jobs => $MFopt{SNPconsJobsPsmpl}, overwrite => $MFopt{redoSNPcons},
				minCallQual => $MFopt{SNPminCallQual},
				);
				

			push @pendingSecondMapSNP, \%SNPinfo;
			#push(@sampleDeps, $consSNPdep) if (defined $consSNPdep && $consSNPdep ne "");
			#die if ($i==2);
		}
	
	}
	
	# One scheduler allocation maps, sorts, calculates coverage, validates, and
	# publishes every requested secondary-reference mapping.
	my $secondMapCmd = $bigMap."\n".$bigCov;
	if ($secondMapCmd =~ /\S/) {
		my @sharedMapWork = grep { defined($_) && $_ ne "" }
			($mapOptHr->{mappingWorkDir}, $mapOptHr->{mapperNodeDir});
		$secondMapCmd .= "\nrm -rf ".join(" ", @sharedMapWork)."\n" if (@sharedMapWork);
	}
	my ($sortJD,$tmpCmd) = ("", "");
	if ($secondMapCmd =~ /\S/) {
		my $secondMapCores = $MFopt{bamSortCores} > ($mapOptHr->{mappingCores} || 0)
			? $MFopt{bamSortCores} : ($mapOptHr->{mappingCores} || $MFopt{bamSortCores});
		my $secondMapSortMem = $MFopt{mapSortMemGb};
		if ($secondMapSortMem < 0) {
			$secondMapSortMem = 20 + (2 * $map{$curSmpl}{inputFileSizeMB}/1024);
			$secondMapSortMem += 20 if ($MFopt{largeMapperDB});
		}
		my $secondMapMem = int($secondMapSortMem) + 1;
		$secondMapMem = $mapOptHr->{mappingMemoryGB}
			if (($mapOptHr->{mappingMemoryGB} || 0) > $secondMapMem);
		my $previousTmpSpace = $QSBoptHR->{tmpSpace};
		my $baseMapHDD = $HDDspace{mapping};
		$baseMapHDD =~ s/G$//;
		$QSBoptHR->{tmpSpace} = int((2.0 * $map{$curSmpl}{inputFileSizeMB} * $baseMapHDD) / 1024) + 30 ."G";
		($sortJD,$tmpCmd) = qsubSystem($dirset{qsubDir}."MAP$SmplName.sh",$secondMapCmd,
			$secondMapCores,$secondMapMem."G",$SmplName."MAP2nd",
			($mapOptHr->{mappingDependencies} || ""),"",1,[],$QSBoptHR);
		$QSBoptHR->{tmpSpace} = $previousTmpSpace;
	}

	append_job_dependencies(\$AsGrps{$cMapGrp}{MapDeps}, $sortJD);
	append_job_dependencies(\$AsGrps{$cAssGrp}{scndMapping}, $sortJD);
	for my $SNPinfo (@pendingSecondMapSNP) {
		my $secondReference = $SNPinfo->{assembly};
		my $secondMapping = $SNPinfo->{MAR}->[0];
		# Do not construct a consensus workflow around future scheduler outputs.
		# A later pipeline pass will return here after both files are published.
		next unless -s $secondReference && -s $secondMapping;
		$SNPinfo->{jdeps} = normalise_job_dependencies($SNPinfo->{jdeps}, $sortJD);
		my $consSNPdep = createConsSNPandSVs($SNPinfo);
		add2SampleDeps($sampleDepsAR, [$consSNPdep]);
	}
	#@cleans = {}; #just deactivate clean up for sec mapping..
}

sub ReadsFromMapping{
	my ($tmpFile,$linkFile) = @_;
	my $newOfile = ""; my $GID = "";
	open TO,">",$linkFile;
	open I,"<",$tmpFile;
	while (my $line = <I>){
		chomp($line);
		my @spl = split(/\t/,$line);
		my $readID = $spl[0];
		my $geneMtch = $spl[2];
		#get read and genomeMatch
		#unless (exists($motus{$geneMtch})){print "could not find gene ID $geneMtch\n"; next;}
		#my $newHd = "@".$COGs2motus{$geneMtch}."___".$motus{$geneMtch};
		$newOfile = $GID.".rdM";
		#die ($newOfile);
		
		$readID=~m/(.*)#/;
		print TO $1."\t"."\t".$newOfile."\n";#$newHd."\n".$spl[9]."\n+\n".$spl[10]."\n";
		#die;
	}
	close TO;
	close I;
}

sub RayAssembly(){
 "mpiexec -n 1 /g/bork5/hildebra/bin/Ray-2.3.1/ray-build/Ray -o test -p test/test_1.fastq test/test_2.fastq -k 31"
 }
 
 
sub buildAssemblyMapIdx{
	my ($finAssLoc,$cAssGrp, $mainRds, $suppRds, $smpl) = @_;
	my %requiredMappers;
	for my $request ([$mainRds, 0], [$suppRds, 1]) {
		next unless $request->[0];
		my $libraries = getRawLibrariesAssmGrp(\%AsGrps,$cAssGrp,$request->[1],$smpl);
		my $readTechnology = libraryTechnology($libraries,
			"assembly-group $cAssGrp ".($request->[1] ? 'support' : 'primary')." mapper index", 1);
		my $mapper = decideMapper($MFopt{MapperProg},$readTechnology);
		# minimap2 and strobealign consume the FASTA directly in mapReadsToRef.
		next if ($mapper == 3 || $mapper == 5);
		$requiredMappers{$mapper} = 1;
	}
	return unless keys %requiredMappers;
	my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0;
	my $submitted = 0;
	for my $mapper (sort {$a <=> $b} keys %requiredMappers) {
		my ($cmdDB) = buildMapperIdx($finAssLoc,$MFopt{MapperCores},$MFopt{largeMapperDB},$mapper);
		next if ($cmdDB eq "");
		my ($jname) = qsubSystem($logDir."mapperIdx.$mapper.sh",$cmdDB,
			int($MFopt{MapperCores}),(int($MFopt{bwtIdxAssMem})+1)."G",
			"mapIdx${mapper}_$JNUM","","",1,[],$QSBoptHR);
		append_job_dependencies(\$AsGrps{$cAssGrp}{AssemblJobName}, $jname);
		$submitted++ if ($jname ne "");
	}
	$QSBoptHR->{tmpSpace} = $tmpSHDD;
	print "Building $submitted mapper index job(s) for assembly $finAssLoc\n" if $submitted;
	
}

 
 sub createPsAssLongReads(){
	my ($jdep, $pseudoAssFile, $Fdir, $smplName) = @_;
	
	my $cleanSeqSetHR = sampleReadSet($curSmpl, "clean");
	my $libraries = readLibrariesByScope($cleanSeqSetHR, 'primary', 1, $curSmpl);
	die "Pseudoassembly for $curSmpl requires singleton long-read libraries\n"
		if @{libraryPairs($libraries)};
	my @nonLong = grep { !$_->{is_long} } @{$libraries};
	die "Pseudoassembly for $curSmpl received non-long libraries: ".join(', ', map { $_->{id} } @nonLong)."\n"
		if @nonLong;
	my @allRds = @{libraryFiles($libraries, 'single')};
	die "Pseudoassembly for $curSmpl has no singleton reads\n" unless @allRds;
	$Fdir =~ s{/+$}{};
	my $psFinal = "$Fdir/".basename($pseudoAssFile);
	my $psStage = "$psFinal.stage";
	my $pseudoAssFileFlag = $psFinal.".sto";
	my $psFile = $psFinal;

	#die "@allRds\n";
	my $renameCtgScr = getProgPaths("renameCtg_scr");#"perl renameCtgs.pl";
	my $sizFiltScr = getProgPaths("sizFilt_scr");#"perl sizeFilterFas.pl";
	my $cmd = "";
	$cmd .= "mkdir -p $Fdir\n";
	$cmd .= "rm -f $psStage\n";
	$cmd .= "$sizFiltScr ".join(",",@allRds)." $MFopt{scaffoldMinSize} -1 $psStage\n";
	$cmd .= "$renameCtgScr $psStage $smplName\n";
	$cmd .= "test -s $psStage\nmv -f $psStage $psFinal\n";
	$cmd .= "touch $pseudoAssFileFlag\n";
	#die $cmd;
	my ($jname,$tmpCmd) = ("","");
	if (!-e $pseudoAssFileFlag){
		$jname = "_PA$JNUM";;
		my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0; 
		($jname,$tmpCmd) = qsubSystem($logDir."pseudoAssembly.sh",$cmd,1,"5G",$jname,$jdep,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR) ;
		$QSBoptHR->{tmpSpace} =$tmpSHDD;

	} else {$cmd="";}
	my $megDir = $psFile; $megDir =~ s/[^\/]+$//;
	#die "$megDir\n";
	return ($jname,$psFile, $megDir);
 }

#spadesAssembly( \%AsGrps,$cAssGrp,"$nodeSpTmpD/ass",$metagAssDir,$MFopt{spadesBayHam} ,$MFglobal{shortAssembly}, $SmplNameX,$hostFilter,$scaffoldFlag) ;
sub spadesAssembly{
	my ($asHr,$cAsGrp,$nodeTmp,$finalOut,$doClean,$helpAssembl,$smplName,$hostFilter,$mateFlag) = @_;
	$finalOut =~ s{/+$}{};
	my $stageOut = "$finalOut.assembly-stage";
	my $backupOut = "$finalOut.previous";

	
	#my $p1ar = $AsGrps{$cAsGrp}{FilterSeq1};my $p2ar = $AsGrps{$cAsGrp}{FilterSeq2};my $singlAr = $AsGrps{$cAsGrp}{FilterSeqS};my $cReadTecAr = $AsGrps{$cAsGrp}{ReadTec};
	my $libraries = getCleanLibrariesAssmGrp($asHr, $cAsGrp, 0);
	my $jDepe = $AsGrps{$cAsGrp}{SeqClnDeps};
	
	


	my $spadesBin = getProgPaths("spades");
	my $isCloudSpades = 0;
	foreach my $library (@{$libraries}){if (($library->{technology} || '') =~ m/SLR/){$isCloudSpades = 1;}}
	if ($isCloudSpades){
		$spadesBin = getProgPaths("cloudspades") ;
		print "Using CloudSpades\n";
	}
	

	#print all samples used 
	my $nCores = $MFopt{AssemblyCores};#6
	my $noTmpOnNode = 0; #prevent usage of tmp space on node
	if ($noTmpOnNode){
		$nodeTmp = $finalOut;
	}
	my $cmd = "rm -rf $nodeTmp $stageOut\nmkdir -p $nodeTmp\n\n";
	$cmd .= "\necho '". $AsGrps{$cAsGrp}{AssemblSmplDirs}. "' > $nodeTmp/smpls_used.txt\n\n";
	$cmd .= "echo \"Starting Spades assembly\"\n";
	my $defTotMem = $MFopt{AssemblyMemory};#60;
	if ($defTotMem == -1){ #auto set mem
		$defTotMem = (spaceInAssGrp($curSmpl)*4+1e5)/1024;
	}

	$cmd .= $spadesBin;
	my $K = $MFopt{AssemblyKmers} ;
	#insert single reads
	my $errStep = "";
	$errStep = "--only-assembler " if ($doClean == 0);
	my $sprds = inputFmtSpadesLibraries($libraries,$logDir);
	$cmd .= " --meta ";
	$cmd .= " $K $sprds -t $nCores $errStep -m $defTotMem ";#--mismatch-correction "; # --meta  --sc "; #> $log #--meta :buggy in 3.6
	$cmd .= " --mismatch-correction " if ($MFopt{spadesMisMatCor});
	if ($helpAssembl ne ""){
		$cmd .= "--untrusted-contigs $helpAssembl ";
	}
	
	$cmd .= "-o $nodeTmp\n";
	#from here could as well be separate 1 core job
	#cleanup assembly
	$cmd .= "\nrm -f -r $nodeTmp/K* $nodeTmp/tmp $nodeTmp/mismatch_corrector/*\n";
	#dual size filter
	$cmd .= "echo \"Assembly finished, renaming contigs, size sorting contigs, collating assembly statistics\"\n";
	my $renameCtgScr = getProgPaths("renameCtg_scr");#"perl renameCtgs.pl";
	my $sizFiltScr = getProgPaths("sizFilt_scr");#"perl sizeFilterFas.pl";
	$cmd .= "$renameCtgScr $nodeTmp/scaffolds.fasta $smplName\n";
	$cmd .= "$sizFiltScr $nodeTmp/scaffolds.fasta $MFopt{scaffoldMinSize} 200\n";
	
	if ($JNUM > 1 && $helpAssembl ne ""){
		#cluster short reads using CD-HIT, replaces scaffolds.fasta.filt2 file
		my $secAss = $nodeTmp."secondary_shorts/";
		$cmd .= "mkdir -p $secAss\n";
		$cmd .= "cat $nodeTmp/scaffolds.fasta >> $nodeTmp/scaffolds.fasta2\n";
		$cmd .= $spadesBin." -s $nodeTmp/scaffolds.fasta2 --only-assembler -m 1000 -t $nCores $K -o $secAss\n";
		#cleanup
		$cmd .= "cp  $secAss/contigs.fasta $nodeTmp/scaffolds.fasta2\nrm -f -r $secAss\n";
	}
	#die "$cmd\n\n";
	my $assStatScr = getProgPaths("assStat_scr");#"perl /g/bork3/home/hildebra/dev/Perl/assemblies/assemblathon_stats.pl";
	$cmd .= "$assStatScr -scaff_size $MFopt{scaffoldMinSize} $nodeTmp/scaffolds.fasta > $nodeTmp/AssemblyStats.500.txt\n";
	$cmd .= "$assStatScr $nodeTmp/scaffolds.fasta > $nodeTmp/AssemblyStats.ini.txt\n";
	$cmd .= "$assStatScr $nodeTmp/scaffolds.fasta.filt > $nodeTmp/AssemblyStats.txt\n";

	my ($cmdDB,$bwtIdx,$chkFile) = buildMapperIdx("$nodeTmp/scaffolds.fasta.filt",$nCores,$MFopt{largeMapperDB},$MFopt{MapperProg});#$nCores);
	$cmd .= $cmdDB unless($mateFlag || !$MFopt{map2Assembly}); #doesn't need bowtie index
	
	#clean up
	$cmd .= "echo \"Zipping non-essential files\"\n";
	$cmd .= "\n $pigzBin -f -p $nCores -r $nodeTmp/scaffolds.fasta $nodeTmp/misc/\n $pigzBin -f -p $nCores -r $nodeTmp/contigs.paths  $nodeTmp/*contigs.fa* \n";
	$cmd .=  "rm -rf $nodeTmp/assembly_graph*.gfa $nodeTmp/corrected $nodeTmp/*.fastg $nodeTmp/before_rr*\n";
	unless ($noTmpOnNode){
		$cmd .=  "mkdir -p $stageOut\ncp -r $nodeTmp/* $stageOut\n" ;
		$cmd .= "rm -rf $nodeTmp\n";
	}
	my $jname = "";
	
	if (-e $logDir."spaderun.sh.otxt"){	#check for out of mem
		open I,"<$logDir/spaderun.sh.otxt" or die "Can't open old assembly logfile $logDir\n"; my $str = join("", <I>); close I;
		if ($str =~ / Error in malloc\(\): out of memory/ ||$str =~ m/TERM_MEMLIMIT: job killed after reaching LSF memory usage limit/){ #memory error for real
			my $replMem  = "";
			if ($str =~ /\n    Max Memory :     (\d+) MB\n/){	$replMem = int(($1/1024)*1.7+0.5);
			} elsif ($str =~ /\nMAX MEM ([\d.]+)G\n/){	$replMem = int($1*1.7+0.5);}
			unless ($replMem eq ""){
				$replMem = 50 if ($replMem < 50);
				$defTotMem = $replMem;
				print $defTotMem."G: new total SPAdes memory\n";
			}
		}
	}
	# Keep SPAdes' own limit aligned with the scheduler request after an OOM retry.
	$cmd =~ s/ -m [\d.]+ / -m $defTotMem /;
	$cmd .= "echo \"MAX MEM ".$defTotMem."G\"\n";
	$cmd .= "echo \"SPADES\" > $stageOut/$checkpointNames{assemblyDone}\n";
	$cmd .= "test -s $stageOut/scaffolds.fasta.filt && test -s $stageOut/AssemblyStats.txt && test -s $stageOut/$checkpointNames{assemblyDone} || exit 37\n";
	$cmd .= "rm -rf $backupOut\nif [ -e $finalOut ]; then mv $finalOut $backupOut; fi\n";
	$cmd .= "if ! mv $stageOut $finalOut; then if [ -e $backupOut ]; then mv $backupOut $finalOut; fi; exit 38; fi\n";
	$cmd .= "rm -rf $backupOut\n";
	#print "in Assembly\n$jDepe\n";
	#print "$finalOut/scaffolds.fasta.filt\n";
	my $locDiskSpace = $HDDspace{assembler};
	if ($locDiskSpace eq "-1" || $locDiskSpace eq "-1G"){
		$locDiskSpace = $HDDspace{spades};
	}
	#die "$locDiskSpace\n";
	unless (-e "$finalOut/scaffolds.fasta.filt" && !-z "$finalOut/scaffolds.fasta.filt" && !-z "$finalOut/AssemblyStats.txt"){
		#my $size_in_mb = (-s $fh) / (1024 * 1024);
		my $tmpCmd="";
		$jname = "_A$JNUM";#$givenJName;
		#$QSBoptHR->{useLongQueue} = 1;
		my $tmpSHDD = $QSBoptHR->{tmpSpace};
		# Every SPAdes branch writes its complete work tree to $nodeTmp.
		# Host selection must not decide whether the matching scratch request is
		# attached to the submitted job.
		$QSBoptHR->{tmpSpace} = $locDiskSpace;
		if ($hostFilter || $MFopt{SpadesAlwaysHDDnode}){
			$QSBoptHR->{useLongQueue} = $MFopt{SpadesLongtime};
			#$QSBoptHR->{tmpSpace} = $HDDspace{spades}; #set option how much tmp space is required, and reset afterwards
			($jname,$tmpCmd) = qsubSystem($logDir."spaderun.sh",$cmd,int($nCores),int($defTotMem)."G",$jname,$jDepe,"",1,$QSBoptHR->{Spades_Hosts},$QSBoptHR) ;
			$QSBoptHR->{useLongQueue} = 0;
		} else {
			($jname,$tmpCmd) = qsubSystem($logDir."spaderun.sh",$cmd,int($nCores),int($defTotMem)."G",$jname,$jDepe,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR) ;
		}
		$QSBoptHR->{tmpSpace} = $tmpSHDD;
		#$QSBoptHR->{useLongQueue} = 0;
	} else {
		print "Spades: Assembly already present in final location\n";
	}
	#die("SPADE\n");
	return ($jname);
}



sub movePreAssmData{
	my ($metagD, $mvD,$mapD, $CSdir, $smplID, $breakpointTsv, $groupID ) = @_; #$tmpD ,
	# Staging and backup directories must be siblings of the final package.
	# Callers commonly include a trailing slash; strip it before suffixing $mvD,
	# otherwise the rotation below tries to move the package into itself.
	$mvD =~ s{/+$}{};
	#very thorough checks that everything is correctly prepped
	my $mvSTO = "$mvD/moved.sto";
#	system "rm -fr $metagD;\n" if ($AssemblyGo && -e $mvSTO);
	if ( -e $mvSTO){
		print "Moved assmbl already.. $mvSTO \n";
		return;
	}
	print "Preparing preassembly package..";
	#die;
	die "movePreAssmData::Coverage does not exist in $CSdir .. can't move preAssembly\n" if (!fileGZe("$CSdir/Coverage.percontig"));
	die "movePreAssmData::Coverage.median.percontig not in $CSdir .. aborting\n" if (!fileGZe("$CSdir/Coverage.median.percontig"));
	die "movePreAssmData::Not a preassebmly?: $metagD\n" unless (-e "$metagD/$checkpointNames{preAssemblyDone}");
	die "movePreAssmData::Couldn't find ContigStats in $metagD\n" unless (-d $CSdir);
	die "movePreAssmData::Couldn't find Assembly $metagD/scaffolds.fasta.filt\n" unless (-e "$metagD/scaffolds.fasta.filt");
	my $mappingCov = "$mapD/$smplID-smd.bam.coverage.gz";
	die "movePreAssmData::Mapping coverage not found: $mappingCov\n" unless (-s $mappingCov);
	die "movePreAssmData::Breakpoint report not found: $breakpointTsv\n" unless (-s $breakpointTsv);
	my $cmd = "";
	# Build a complete sibling directory first. An interrupted copy can never be
	# mistaken for a usable package, and an older partial package is retained for diagnosis.
	my $packageTag = time().".$$";
	my $stage = "$mvD.stage.$packageTag";
	$cmd .= "rm -rf $stage; mkdir -p $stage;\n";
	$cmd .= "cp -rf $metagD/scaffolds.fasta.filt $metagD/$checkpointNames{preAssemblyDone} $CSdir/Coverage* $stage;\n";
	$cmd .= "cp $mappingCov $stage/mapping.coverage.gz; cp $breakpointTsv $stage/breakpoints.tsv.gz;\n";
	$cmd .= "printf 'key\\tvalue\\nschema_version\\t2\\nsample_id\\t$smplID\\nassembly_group\\t$groupID\\nsimulator\\tcoverage_weighted_v2\\nmax_synthetic_depth\\t$MFopt{hybridSyntheticMaxDepth}\\n' > $stage/package.manifest.tsv;\n";
	$cmd .= "test -s $stage/scaffolds.fasta.filt && test -s $stage/mapping.coverage.gz && test -s $stage/breakpoints.tsv.gz && test -s $stage/package.manifest.tsv || exit 35;\n";
	if (-e "$metagD/AssemblyStats.txt"){
		$cmd .= "mv $metagD/AssemblyStats.txt $logDir/preAssmStat.txt;\n";
	} elsif (-e "$metagD/AssemblyStats.ini.txt"){
		$cmd .= "mv $metagD/AssemblyStats.ini.txt $logDir/preAssmStat.txt;\n";
	}
	# Mapper indexes remain with the source assembly; packaging is non-destructive.
	#if ($AssemblyGo){$cmd .= "rm -fr $metagD;\n" ; print "removed preASsmbl dir";}
	#$cmd .= "cp $tmpD/$checkpointNames{preAssemblyDone} $metagD;\n";
	#$cmd .= "cp -rf $tmpD/* $mvD/;\n";
	#$cmd .= "rm -fr $tmpD\n";
	$cmd .= "touch $stage/moved.sto\n";
	$cmd .= "if [ -e $mvD ]; then mv $mvD $mvD.incomplete.$packageTag; fi; mv $stage $mvD\n";
	#die $cmd;
	systemW $cmd;
	print " Done \n";
	#return $newCovFile;
	return;
}

#used in hybrid assemblies
sub prepPreAssmbl{
	my ($metagD, $mvD,$mapD, $tmpD , $CSdir,  $cAssGrp, $finAssLoc,$finalCommAssDir) = @_;
	#print "$mvD\n";
	
	$AsGrps{$cAssGrp}{CntPreAss} = 0 unless (exists($AsGrps{$cAssGrp}{CntPreAss}));
	$AsGrps{$cAssGrp}{CntPreAssMiss} = 0 unless (exists($AsGrps{$cAssGrp}{CntPreAssMiss}));
	$AsGrps{$cAssGrp}{CntPreAssNoPrim} = 0 unless (exists($AsGrps{$cAssGrp}{CntPreAssNoPrim}));
	
	my $hasPrimary = $map{$curSmpl}{hasPrimaryRds};
	
	#my $hasPrimaryRds= 1;$hasPrimaryRds = 0 if ($map{$curSmpl}{prefix} eq "" && $map{$curSmpl}{dir} eq "");
	
	my $hybridAssemblyRequested = $MFopt{DoAssembly} == 5
		&& ($map{$curSmpl}{SupportReads} || "") =~ m/(?:PB|ONT):/;
	print "precnt: $AsGrps{$cAssGrp}{CntPreAss}\n"
		if ($hybridAssemblyRequested && !$MFconfig{silent});

	my $ePreAssmbly = 0; $ePreAssmbly = 1 if (-s $finAssLoc && -e "$finalCommAssDir/$checkpointNames{preAssemblyDone}");
	#die "$ePreAssmbly\n";
	my $ePreAssmblPck = 0; 
	if (hybrid_package_complete($mvD)) {
		$ePreAssmbly=1;$ePreAssmblPck = 1;
	}
	my $doPreAssmFlag = 0;

	if (-s $finAssLoc && -e "$finalCommAssDir/$checkpointNames{assemblyDone}"){ #indication that hybrid assembly is already done
		#die;
		# The package remains on disk, but it must no longer gate mappings and
		# downstream analyses of the final assembly.
		return ($ePreAssmbly,$doPreAssmFlag,0,0);
	}
	
	if (!$hasPrimary){#should not be included at all: nothing to assemble within preassembly..
		$AsGrps{$cAssGrp}{CntPreAssNoPrim}++ ;
		my $postAssemblyGo = hybrid_group_ready(
			$AsGrps{$cAssGrp}{CntPreAss},
			$AsGrps{$cAssGrp}{CntPreAssNoPrim},
			$AsGrps{$cAssGrp}{CntAimAss},
		);
		return ($ePreAssmbly,0,$postAssemblyGo,$ePreAssmblPck);
	}
	my $eCOV = 0; $eCOV = 1 if (fileGZe( "$CSdir/Coverage.percontig"));
	my $eCOVmv = 0; $eCOVmv = 1 if (fileGZe("$mvD/Coverage.percontig"));
	#print "$eCOV $eCOVmv   $CSdir/Coverage.percontig    $mvD/Coverage.percontig\n";
	if ($MFopt{DoAssembly} == 5 && $AsGrps{$cAssGrp}{SupportReads} =~ m/(?:PB|ONT):/){#$map{$curSmpl}{"SupportReads"} =~ m/PB:/ ){
		#condition: right assembly mode and actually secondary support reads
		$doPreAssmFlag = 1 ;
		#print "XAS\n";
		if ((!$eCOVmv && !$eCOV) || !$ePreAssmbly || $map{$curSmpl}{inputFilesEmpty}){
			#print "preAssmbl: nothing done yet.. \n$mvD\n$metagD\n";
			#die;
			if ($map{$curSmpl}{inputFilesEmpty}){
				$AsGrps{$cAssGrp}{CntPreAssNoPrim} ++ ; #needs to be counted as "existing"
			}
			return ($ePreAssmbly,$doPreAssmFlag, 0, $ePreAssmblPck );
		}
	} else {
		#die;
		#print "out3\n";
		return (0,0,0,0); 
	}
	
	#my $eCOV = 0; $eCOV =1 if ( -e "$CSdir/Coverage.percontig");
	
	#this moves the contig stats (and removes in ori location)
	if ((!$ePreAssmblPck || (($eCOV || !$hasPrimary) && ($eCOVmv || $map{$curSmpl}{"SupportReads"} eq "" ))) && -e "$metagD/$checkpointNames{preAssemblyDone}" ){
		#die "preAssmX: $PostAssemblyGo $doPreAssmFlag     $AsGrps{$cAssGrp}{CntPreAss} >= $AsGrps{$cAssGrp}{CntAimAss}\n";
		$doPreAssmFlag = 0 ;#no prep needed any longer.. files will/are saved already!
		#all ready for second assembly step!
		movePreAssmData($metagD, $mvD,$mapD,  $CSdir, $curSmpl,
			"$mapD/$curSmpl-smd.bam.breakpoints.tsv.gz", $cAssGrp) ;#$tmpD ,
		$ePreAssmblPck = 1;
	}  
	#die "XAS\n";
	if ($ePreAssmblPck){
		#print "$ePreAssmblPck || ($eCOV || !$hasPrimary) && ($eCOVmv || $map{$curSmpl}{SupportReads} eq \"\" )\n";
		if ($ePreAssmblPck || () ){
			$AsGrps{$cAssGrp}{CntPreAss} ++ ;
			push(@{$AsGrps{$cAssGrp}{preAsmblDir}}, $mvD);
		}
		$doPreAssmFlag = 0 ; #everyone else needs to keep
		
	}
	my $PostAssemblyGo = 0;
	#print "UUU $doPreAssmFlag\n";
	my $finJobs = ($AsGrps{$cAssGrp}{CntPreAss}+$AsGrps{$cAssGrp}{CntPreAssNoPrim} ); #+ $AsGrps{$cAssGrp}{CntPreAssMiss}
	print "FIN: $finJobs ($AsGrps{$cAssGrp}{CntPreAss} packages + "
		."$AsGrps{$cAssGrp}{CntPreAssNoPrim} without-primary) >= $AsGrps{$cAssGrp}{CntAimAss}\n";
	$PostAssemblyGo = 1 if (!$doPreAssmFlag && ( $finJobs >= $AsGrps{$cAssGrp}{CntAimAss}) ); #has already seen enough complete preAssmblies
	$doPreAssmFlag = 1 if (!$PostAssemblyGo); 
	#print "-e $CSdir/Coverage.percontig   $metagD/$checkpointNames{preAssemblyDone}\n" ;
	#print "preAssm:  $doPreAssmFlag     $AsGrps{$cAssGrp}{CntPreAss} +$AsGrps{$cAssGrp}{CntPreAssNoPrim} >= $AsGrps{$cAssGrp}{CntAimAss} :: $ePreAssmblPck $PostAssemblyGo\n";
	#die "$doPreAssmFlag\n";
	return ($ePreAssmbly,$doPreAssmFlag,$PostAssemblyGo,$ePreAssmblPck);
}


sub longRdAssembly{
	my ($asHr,$cAsGrp,$nodeTmp,$finalOut,$helpAssembl,$smplName, $useSupportRds,$LassP) = @_;
	# $finalOut is normally passed as ".../metag/".  Staging and backup
	# directories must be siblings of that directory; suffixing a trailing-slash
	# path creates children such as "metag/.hybrid-stage", after which the
	# publication rotation attempts to move metag into itself and the next loop
	# resubmits metaMDBG.  Canonicalise before constructing the atomic paths.
	$finalOut =~ s{/+$}{};
	
	my $libraries = getCleanLibrariesAssmGrp($asHr, $cAsGrp, $useSupportRds);
	my $pairs = libraryPairs($libraries);
	if (@{$pairs}){print "Paired reads are defined, but long-read assemblies rely on singleton reads!\nAborting\n";die;}
	my $singlAr = libraryFiles($libraries, 'single');
	my $numInLibs = scalar @{$singlAr};
	my %long_read_tech = map { ($_->{technology} || '') => 1 }
		grep { ($_->{technology} || '') =~ /^(?:PB|ONT)$/ } @{$libraries};
	die "Assembly group $cAsGrp mixes ONT and PacBio reads; split these technologies into separate assembly groups\n"
		if keys(%long_read_tech) > 1;
	my ($long_read_tech) = keys %long_read_tech;
	#print "$numInLibs libs\n";
		
	if ($LassP == 5){#hybrid mode.. check that all required files are present or stop here
		#if (${$cReadTecAr}[0] ne "PB"){print"Expected PacBio (\"PB\") readTech for metaMDBG, found \"${$cReadTecAr}[0]\"\n";die;}
		my $long_reads_detected = keys(%long_read_tech) ? 1 : 0;
		if (!$long_reads_detected){
			print "Hybrid Assembly.. expected \"PB\" or \"ONT\" reads for support reads!";
			print "(found technologies \"".join(',', map { $_->{technology} || '' } @{$libraries})."\")\nAborting..\n";
			die;
		}
		
		return "" unless (-e $helpAssembl || $helpAssembl eq "hybridmMDBG");
	}
	#my ($p1arX,$p2arX,$singlArX,$cReadTecArX) = getCleanSeqsAssmGrp($asHr, $cAsGrp, 1);

	#my $singlAr = $AsGrps{$cAsGrp}{FilterSeqS};#my $cReadTec = $AsGrps{$cAsGrp}{ReadTec};
	my $jDepe = $AsGrps{$cAsGrp}{SeqClnDeps};
	
	my $nameProg= "flye"; $nameProg="mMDBG"if($MFopt{DoAssembly}==4 || $MFopt{DoAssembly}==5);
	if (grep { ($_->{technology} || '') =~ m/SLR/i } @{$libraries}){die "Can't use synthetic long reads (SLR) with $nameProg\n";}
	my $nCores = $MFopt{AssemblyCores};#6
	
	my $nodeTmp2 = "$nodeTmp/tmpRawRds/";
	my $stageOut = "$finalOut.hybrid-stage";
	my $backupOut = "$finalOut.previous";
	my $cmd = "rm -rf $nodeTmp $stageOut\nmkdir -p $nodeTmp $nodeTmp2\n\n"; #\n  mkdir -p $nodeTmp/tmp\n
	#input reads for assembly .. expected unpaired, long reads
	my @inRds;# = @{$singlAr};
	my @hybridPreassemblies;
	
	#metaMDBG "hack" to impute illumina assemblies:
	my $cmdPre = "";
	if ($helpAssembl eq "hybridmMDBG" && $MFopt{DoAssembly} == 5){
		#$cmd .= "echo \"Hybrid Assembly prep\"\n";
		my $spl4m = getProgPaths("split_fasta4metaMDBG_scr");
		my $runPar=1;
		#my $illPathS = ;
		my @illDirs = @{$AsGrps{$cAsGrp}{preAsmblDir}}; #split /,/,$illPathS;
		@hybridPreassemblies = map { "$_/scaffolds.fasta.filt" } @illDirs;

		$cmdPre .= "echo \"presplitting helper assembly\"\n";
		my @simulatorCommands;
		my ($lengthTemplate) = grep { defined($_) && $_ ne "" } @{$singlAr};
		my $lengthTemplateArg = defined($lengthTemplate) ? "--length-template $lengthTemplate" : "";
		# Preassembly packages and long-read libraries are independent collections:
		# generate one synthetic input per package, then append every PB/ONT input.
		for (my $i=0;$i<@illDirs;$i++){
			my $illD = $illDirs[$i];
			die "longRdAssembly:: $illD is not a dir!" unless (-d $illD);
			my $contigCov = "$illD/mapping.coverage.gz";
			die "Hybrid package $illD lacks mapping.coverage.gz\n" unless fileGZe($contigCov);
			my $breakpointTsv = "$illD/breakpoints.tsv.gz";
			die "Hybrid package $illD lacks breakpoints.tsv.gz\n" unless -s $breakpointTsv;
			my $packageSample = hybrid_package_sample_id($illD);
			die "Hybrid package $illD lacks a sample_id in package.manifest.tsv\n"
				if ($packageSample eq "");
			$packageSample =~ s/[^A-Za-z0-9_.-]+/_/g;
			my $dupiAssmbl = "$nodeTmp2/$packageSample.synthetic.fastq.gz";
			my $preAssmbl = "$illD/scaffolds.fasta.filt";
			push @simulatorCommands,
				"( $spl4m --assembly $preAssmbl --coverage $contigCov --breakpoints $breakpointTsv --output $dupiAssmbl "
				."--mean-read-length $MFconfig{defaultReadLengthX} $lengthTemplateArg --max-synthetic-depth $MFopt{hybridSyntheticMaxDepth} )";
			$inRds[$i] = $dupiAssmbl;
		}
		#instead of combining ill + PacBio fastas in $cmdLater, just add PB as extra libs (makes more sense in any case..)
		for (my $i=0;$i<@{$singlAr};$i++){
			push(@inRds,$singlAr->[$i]);
		}
		#transfer commands to main #cmd stream..
		if ($runPar) {
			my $batchSizes = balanced_parallel_batches(scalar(@simulatorCommands), 20);
			my $simulatorIndex = 0;
			for (my $batch = 0; $batch < @{$batchSizes}; $batch++) {
				my $batchEnd = $simulatorIndex + $batchSizes->[$batch];
				$cmdPre .= "\necho \"Starting synthetic-read simulator batch ".($batch + 1)
					."/".scalar(@{$batchSizes})." (".$batchSizes->[$batch]." jobs)\"\n";
				for (my $i = $simulatorIndex; $i < $batchEnd; $i++) {
					$cmdPre .= $simulatorCommands[$i]." &\nsim_pid_$i=\$!\n";
				}
				$cmdPre .= "sim_failed=0\n";
				for (my $i = $simulatorIndex; $i < $batchEnd; $i++) {
					$cmdPre .= "if ! wait \$sim_pid_$i; then echo \"synthetic-read simulator $i failed\" >&2; sim_failed=1; fi\n";
				}
				$cmdPre .= "[ \$sim_failed -eq 0 ] || exit 33\n";
				$simulatorIndex = $batchEnd;
			}
			$cmdPre .= "\n";
		}

		$cmd .= $cmdPre."\n\n";
		for (my $i=0;$i<@illDirs;$i++){
			next if ($inRds[$i] eq "");
			$cmd .= "if [ ! -s $inRds[$i] ]; then echo \"$inRds[$i] not present.. abort\"; exit 33; fi \n";
			$cmd .= "$pigzBin -t $inRds[$i] || exit 33\n";
		}
		$cmd .= "echo \"presplitting done\"\n\n";
		#die "$cmd\n";
	} else {
		@inRds = @{$singlAr};
	}
	
	#die "@inRds\n";
	#print $cmd;#die;
	
	$cmd .= "echo \"Starting $nameProg assembly\"\n";
	
	my $contigRecovery = "";
	if ($MFopt{DoAssembly}==3){#FLYE
		my $technology = $long_read_tech || '';
		if ($technology ne "ONT"){print"Expected Oxford Nanopore (\"ONT\") readTech for flye, found \"$technology\"\n";die;}
		my $flyeBin = getProgPaths("flye");
		$cmd .= $flyeBin;
		$cmd .= " --nano-raw ".join(" ", @inRds)." -t $nCores --meta -g 3g ";
		if ($helpAssembl ne ""){
			$cmd .= "--subassemblies $helpAssembl ";
		}
		$cmd .= "-o $nodeTmp\n";  #--tmp-dir $nodeTmp/tmp/
		#$outAssemblyF = "$nodeTmp/assembly.fasta";
		$contigRecovery .= "\nrm -fr $nodeTmp/00-assembly/ $nodeTmp/10-consensus/ $nodeTmp/20-repeat/ $nodeTmp/30-contigger/ $nodeTmp/40-polishing/ $nodeTmp/assembly_graph.gv\n";
		$contigRecovery .= "\nmv $nodeTmp/assembly.fasta $nodeTmp/scaffolds.fasta\n\n";
	} elsif($MFopt{DoAssembly}==4 || $MFopt{DoAssembly}==5){ #metaMDBG
		my $inFileFlag = "--in-hifi";
		$inFileFlag = "--in-ont" if (defined($long_read_tech) && $long_read_tech eq "ONT");
		my $mMDBG = getProgPaths("metaMDBG");
		$cmd .= "$mMDBG asm --threads $nCores --out-dir $nodeTmp $inFileFlag " . join(" ",@inRds) . "\n";
		$cmd .= "rm -rf $nodeTmp/tmp/;\n";
		# metaMDBG publishes gzip-compressed contigs, while the shared assembly
		# cleanup below requires an uncompressed scaffolds.fasta. Convert once,
		# validate it, then discard the redundant compressed copy.
		$contigRecovery .= "test -s $nodeTmp/contigs.fasta.gz || exit 34\n";
		$contigRecovery .= "$pigzBin -dc -p $nCores $nodeTmp/contigs.fasta.gz > $nodeTmp/scaffolds.fasta || exit 34\n";
		$contigRecovery .= "test -s $nodeTmp/scaffolds.fasta || exit 34\n";
		$contigRecovery .= "rm -f $nodeTmp/contigs.fasta.gz\n\n";
		
		#contigs.fasta.gz
	} else {
		die "Wrong assembler selected\n";
	}
	$cmd .= "\nrm -rf $nodeTmp2\n";
	
	#$cmd .= getProgPaths("activateBase")."\n"; #not needed any longer.. all assemblers are in base env
	
	
	#from here could as well be separate 1 core job
	#cleanup assembly
	#$cmd .= "\necho '". chomp($AsGrps{$cAsGrp}{AssemblSmplDirs}). "' > $nodeTmp/smpls_used.txt\n\n";
	$cmd .= "\necho '". $AsGrps{$cAsGrp}{AssemblSmplDirs}. "' > $nodeTmp/smpls_used.txt\n\n";

	$cmd .= $contigRecovery;#
	#dual size filter
	$cmd .= "echo \"Assembly finished, renaming contigs, size sorting contigs, collating assembly statistics\"\n";
	my $renameCtgScr = getProgPaths("renameCtg_scr");#"perl renameCtgs.pl";
	my $sizFiltScr = getProgPaths("sizFilt_scr");#"perl sizeFilterFas.pl";
	$cmd .= "$renameCtgScr $nodeTmp/scaffolds.fasta $smplName\n";
	$cmd .= "$sizFiltScr $nodeTmp/scaffolds.fasta $MFopt{scaffoldMinSize} 200\n";
	
	my $assStatScr = getProgPaths("assStat_scr");#"perl /g/bork3/home/hildebra/dev/Perl/assemblies/assemblathon_stats.pl";
	$cmd .= "$assStatScr -scaff_size $MFopt{scaffoldMinSize} $nodeTmp/scaffolds.fasta > $nodeTmp/AssemblyStats.500.txt\n";
	$cmd .= "$assStatScr $nodeTmp/scaffolds.fasta > $nodeTmp/AssemblyStats.ini.txt\n";
	$cmd .= "$assStatScr $nodeTmp/scaffolds.fasta.filt > $nodeTmp/AssemblyStats.txt\n";
	if (@hybridPreassemblies) {
		my $compareScr = getProgPaths("compareHybridAssemblies_scr");
		my $comparisonPreassembly = $hybridPreassemblies[0];
		$cmd .= "$compareScr --final $nodeTmp/scaffolds.fasta.filt --preassembly $comparisonPreassembly --output $nodeTmp/HybridAssemblyComparison.tsv\n";
		$cmd .= "[ -s $nodeTmp/HybridAssemblyComparison.tsv ] || exit 36\n";
	}
	
	my ($cmdDB,$bwtIdx,$chkFile) = buildMapperIdx("$nodeTmp/scaffolds.fasta.filt",$nCores,$MFopt{largeMapperDB},$MFopt{MapperProg});#$nCores);
	$cmd .= $cmdDB unless( !$MFopt{map2Assembly}); #doesn't need bowtie index
	
	#clean up
	$cmd .= "echo \"Zipping non-essential files\"\n";
	$cmd .= "\n $pigzBin -f -p $nCores $nodeTmp/scaffolds.fasta\n";
	$cmd .=  "mkdir -p $stageOut\ncp -r $nodeTmp/* $stageOut\n" ;
	$cmd .= "rm -rf $nodeTmp\n";
	my $jname = "";

	my $defTotMem = $MFopt{AssemblyMemory};#60;
	if ($defTotMem == -1){ #auto set mem
		$defTotMem = (spaceInAssGrp($curSmpl,$useSupportRds)*8+1e4)/1024;
	}

	my $defMem = ($defTotMem);

	$cmd .= "echo \"MAX MEM ".$defTotMem."G\"\n";
	$cmd .= "echo \"$nameProg\" > $stageOut/$checkpointNames{assemblyDone}\n";
	$cmd .= "test -s $stageOut/scaffolds.fasta.filt && test -s $stageOut/AssemblyStats.txt && test -s $stageOut/$checkpointNames{assemblyDone} || exit 37\n";
	$cmd .= "test -s $stageOut/HybridAssemblyComparison.tsv || exit 37\n" if @hybridPreassemblies;
	$cmd .= "rm -rf $backupOut\n";
	$cmd .= "if [ -e $finalOut ]; then mv $finalOut $backupOut; fi\n";
	$cmd .= "if ! mv $stageOut $finalOut; then if [ -e $backupOut ]; then mv $backupOut $finalOut; fi; exit 38; fi\n";
	$cmd .= "rm -rf $backupOut\n";
	#die "$cmd\n\n";
	if (!-e "$finalOut/scaffolds.fasta.filt"  || !-e "$finalOut/AssemblyStats.txt"){
		#my $size_in_mb = (-s $fh) / (1024 * 1024);
		my $tmpCmd="";
		$jname = "$nameProg$JNUM";#$givenJName;
		$QSBoptHR->{useLongQueue} = 0;#super fast, doesn't need long queue
		my $tmpSHDD = $QSBoptHR->{tmpSpace};
		my $assemblerScratchGB = $HDDspace{assembler};
		$assemblerScratchGB =~ s/G$//;
		if ($assemblerScratchGB < 0) {
			$assemblerScratchGB = $nameProg eq "mMDBG"
				? $HDDspace{metaMDBG}
				: int(25 + (4 * spaceInAssGrp($curSmpl,$useSupportRds) / 1024));
			$assemblerScratchGB =~ s/G$//;
		}
		if ($nameProg eq "mMDBG") {
			my $preassemblyBytes = 0;
			$preassemblyBytes += -s $_ for grep { defined($_) && -s $_ } @hybridPreassemblies;
			$assemblerScratchGB = hybrid_local_scratch_gb(
				assembler_gb => $assemblerScratchGB,
				preassembly_bytes => $preassemblyBytes,
				max_synthetic_depth => $MFopt{hybridSyntheticMaxDepth},
			);
		}
		# Flye and metaMDBG both create their complete assembly work trees below
		# $nodeTmp, so both require an explicit scheduler scratch request.
		$QSBoptHR->{tmpSpace} = $assemblerScratchGB."G";
		if ( $MFopt{SpadesAlwaysHDDnode}){
			($jname,$tmpCmd) = qsubSystem($logDir."$nameProg.sh",$cmd,(int($nCores)),int($defMem)."G",$jname,$jDepe,"",1,$QSBoptHR->{Spades_Hosts},$QSBoptHR) ;
		} else {
			($jname,$tmpCmd) = qsubSystem($logDir."$nameProg.sh",$cmd,(int($nCores)),int($defMem)."G",$jname,$jDepe,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR) ;
		}
		$QSBoptHR->{tmpSpace} = $tmpSHDD;
		$QSBoptHR->{useLongQueue} = 0;
	} else {
		print "longReadAssm: Assembly already present in final location\n";
	}
	return ($jname);

	
	
}
#( \%AsGrps,$cAssGrp,"$nodeSpTmpD/ass",$metagAssDir ,$MFglobal{shortAssembly}, $SmplNameX,$hostFilter,$scaffoldFlag)
sub megahitAssembly{
	my ($asHr,$cAsGrp,$nodeTmp,$finalOut,$helpAssembl,$smplName,$hostFilter,$mateFlag) = @_;
	$finalOut =~ s{/+$}{};
	my $stageOut = "$finalOut.assembly-stage";
	my $backupOut = "$finalOut.previous";
	my $megahitBin = getProgPaths("megahit");
	my $stoneAssmbl = "$checkpointNames{assemblyDone}";
	
#	my $p1ar = $AsGrps{$cAsGrp}{FilterSeq1};my $p2ar = $AsGrps{$cAsGrp}{FilterSeq2};my $singlAr = $AsGrps{$cAsGrp}{FilterSeqS};
	my $jDepe = $AsGrps{$cAsGrp}{SeqClnDeps};
#	my $cReadTec = $AsGrps{$cAsGrp}{ReadTec};
	my $libraries = getCleanLibrariesAssmGrp($asHr, $cAsGrp, 0);

	if (grep { ($_->{technology} || '') =~ m/SLR/i } @{$libraries}){die "Can't use synthetic long reads (SLR) with megahit\n";}
	#print all samples used 
	my $nCores = $MFopt{AssemblyCores};#6
	my $noTmpOnNode = 0; #prevent usage of tmp space on node
	if ($noTmpOnNode){
		$nodeTmp = $finalOut;
	}
	my $nodePreD = $nodeTmp;$nodePreD=~s/\/[^\/]+\/*$/\//;
	#die "$nodePreD\n$nodeTmp\n";
	
	my $cmd = "rm -rf $nodeTmp $stageOut\nmkdir -p $nodePreD\n\n"; #\n  mkdir -p $nodeTmp/tmp\n
	$cmd .= "echo \"Starting megaHit assembly\"\n";
	my $inputSizeloc = spaceInAssGrp($curSmpl);
	#die "DS$inputSizeloc\n";
	my $defTotMem = $MFopt{AssemblyMemory};#60;
	if ($defTotMem == -1){ #auto set mem
		#die "ts:$inputSloc\n";
		$defTotMem = ($inputSizeloc*1.85 + 5e4)/1024;
	}
	my $defMem = int($defTotMem);
	#$defTotMem = $defMem; #total really available mem (in GB)
	
	my $locDiskSpace = $HDDspace{assembler};
	if ($locDiskSpace eq "-1G" || $locDiskSpace eq "-1"){
		$locDiskSpace = int(25+(4*$inputSizeloc/1024))."G";
		#$locDiskSpace = $HDDspace{megaHit};
	}
	#print "diskspace: $locDiskSpace\n";die;


	my $K = $MFopt{AssemblyKmers} ;
	$K =~ s/-k //;
	#check that k's are closer than 28
	my @spl  = sort {$a <=> $b}(split /,/,$K);
	for (my $i=1;$i<@spl;$i++){
		if ($spl[$i] - $spl[$i-1] > 28){die "kmers for megahit: $K: steps must be <28\n";}
	}
	$K = join(",",@spl);
	#insert single reads
	my $numInLibs = scalar @{$libraries};
	my ($megahitInputSetup, $sprds) = inputFmtMegahitRuntimeLibraries($libraries, 'megahit_inputs');
	$cmd .= $megahitInputSetup;
	$cmd .= $megahitBin;
	$cmd .= " --k-list $K $sprds -t $nCores -m ". int($defTotMem*1024*1024*1024*0.8) ." --out-prefix megaAss ";
	if ($helpAssembl ne ""){
		if ($helpAssembl eq "preAssmbl"){ #check for keywords
			$stoneAssmbl = $checkpointNames{preAssemblyDone};
			$helpAssembl = "";
		} elsif (-f $helpAssembl) {
			die "MEGAHIT does not support a file-valued helper assembly ('$helpAssembl'); use SPAdes or the hybrid preassembly workflow\n";
		} else { die "Can't decipher helpAssmbl input to megahitAssembly: $helpAssembl\n";}
	}
	$cmd .= "-o $nodeTmp  \n";  #--tmp-dir $nodeTmp/tmp/
	#from here could as well be separate 1 core job
	#cleanup assembly
	$cmd .= "echo \"Assembly finished, renaming contigs, size sorting contigs, collating assembly statistics\"\n";
	$cmd .= "\necho '". $AsGrps{$cAsGrp}{AssemblSmplDirs}. "' > $nodeTmp/smpls_used.txt\n\n";
	$cmd .= "\nrm -fr $nodeTmp/tmp/ $nodeTmp/intermediate_contigs/ \nmv $nodeTmp/megaAss.contigs.fa $nodeTmp/scaffolds.fasta\n\n";
	#dual size filter
	my $renameCtgScr = getProgPaths("renameCtg_scr");#"perl renameCtgs.pl";
	my $sizFiltScr = getProgPaths("sizFilt_scr");#"perl sizeFilterFas.pl";
	$cmd .= "$renameCtgScr $nodeTmp/scaffolds.fasta $smplName\n";
	$cmd .= "$sizFiltScr $nodeTmp/scaffolds.fasta $MFopt{scaffoldMinSize} 200\n";
	
	if ($JNUM > 1 && $helpAssembl ne ""){
		die("helper assemblies not supported with megahit: $helpAssembl\n");
		#cluster short reads using CD-HIT, replaces scaffolds.fasta.filt2 file
		my $secAss = $nodeTmp."secondary_shorts/";
		#$cmd .= "mkdir -p $secAss\n";
		#$cmd .= "cat $nodeTmp/scaffolds.fasta >> $nodeTmp/scaffolds.fasta2\n";
		#$cmd .= $spadesBin." -s $nodeTmp/scaffolds.fasta2 --only-assembler -m 1000 -t $nCores $K -o $secAss\n";
		#cleanup
		#$cmd .= "cp  $secAss/contigs.fasta $nodeTmp/scaffolds.fasta2\nrm -f -r $secAss\n";
	}
	
	
	my $assStatScr = getProgPaths("assStat_scr");#"perl /g/bork3/home/hildebra/dev/Perl/assemblies/assemblathon_stats.pl";
	$cmd .= "$assStatScr -scaff_size $MFopt{scaffoldMinSize} $nodeTmp/scaffolds.fasta > $nodeTmp/AssemblyStats.500.txt\n";
	$cmd .= "$assStatScr $nodeTmp/scaffolds.fasta > $nodeTmp/AssemblyStats.ini.txt\n";
	$cmd .= "$assStatScr $nodeTmp/scaffolds.fasta.filt > $nodeTmp/AssemblyStats.txt\n";
	
	my ($cmdDB,$bwtIdx,$chkFile) = buildMapperIdx("$nodeTmp/scaffolds.fasta.filt",$nCores,$MFopt{largeMapperDB},$MFopt{MapperProg});#$nCores);
	
	
	#DEBUG
	#print "\n\n$cmdDB ==  || $mateFlag || !$MFopt{map2Assembly} \n";
	
	unless($cmdDB eq "" || $mateFlag || !$MFopt{map2Assembly}){
		$cmd .= "echo \"Building mapper index on assembly\"\n";
		$cmd .= $cmdDB ; #doesn't need bowtie index
	}
	#print "DB:: $cmdDB\nMATE::$mateFlag\n";
	#clean up
	$cmd .= "echo \"Zipping non-essential files\"\n";
	$cmd .= "\n $pigzBin -f -p $nCores $nodeTmp/scaffolds.fasta.filt2 $nodeTmp/scaffolds.fasta.lnk\n";
	$cmd .= "rm -f $nodeTmp/scaffolds.fasta\n";
	unless ($noTmpOnNode){
		$cmd .=  "mkdir -p $stageOut\ncp -r $nodeTmp/* $stageOut\n" ;
		$cmd .= "rm -rf $nodeTmp\n";
	}
	my $jname = "";
	
#	if (-e $logDir."megahitrun.sh.otxt"){	#check for out of mem
#	}

	$cmd .= "echo \"MAX MEM ".$defTotMem."G\"\n";
	$cmd .= "echo \"MEGAHIT\" > $stageOut/$stoneAssmbl\n";
	$cmd .= "test -s $stageOut/scaffolds.fasta.filt && test -s $stageOut/AssemblyStats.txt && test -s $stageOut/$stoneAssmbl || exit 37\n";
	$cmd .= "rm -rf $backupOut\nif [ -e $finalOut ]; then mv $finalOut $backupOut; fi\n";
	$cmd .= "if ! mv $stageOut $finalOut; then if [ -e $backupOut ]; then mv $backupOut $finalOut; fi; exit 38; fi\n";
	$cmd .= "rm -rf $backupOut\n";
	#die "$cmd\n\n";
	if (!-e "$finalOut/scaffolds.fasta.filt"  || !-e "$finalOut/AssemblyStats.txt"){
		#my $size_in_mb = (-s $fh) / (1024 * 1024);
		my $tmpCmd="";
		$jname = "mA$JNUM";#$givenJName;
		my $tmpSHDD = $QSBoptHR->{tmpSpace};
		$QSBoptHR->{tmpSpace} = $locDiskSpace;#$HDDspace{megaHit};
		($jname,$tmpCmd) = qsubSystem($logDir."megahitrun.sh",$cmd,$nCores,$defMem."G",$jname,$jDepe,"",1,$QSBoptHR->{Spades_Hosts},$QSBoptHR) ;
		$QSBoptHR->{tmpSpace} = $tmpSHDD;
	} else {
		print "MegaHit: Assembly already present in final location\n";
	}
	#die("SPADE\n");
	return ($jname);
}


sub metagAssemblyRun{
	my ( $cAssGrp,$nodeTmp,$metagAssDir, $geneDir, $SmplNameX,$scaffoldFlag,$metaGscaffDir,
				$assemblyFlag,$AssemblyGo,$ePreAssmbly, $doPreAssmFlag, $postAssmblGo,$finalCommAssDir) = @_;
				#metagAssembly( $cAssGrp,"$nodeSpTmpD/ass",$metagAssDir ,$MFglobal{shortAssembly}, $SmplNameX,$scaffoldFlag,
				#	$assemblyFlag,$AssemblyGo,$ePreAssmbly, $doPreAssmFlag, $postPreAssmblGo);
	my $assemblyInputMB = spaceInAssGrp($curSmpl, 1);
	my $assemblyCores = assembly_cores_for_input(
		input_mb => $assemblyInputMB,
		configured_cores => $MFopt{AssemblyCores},
	);
	local $MFopt{AssemblyCores} = $assemblyCores;
	printf "Assembly step using %d core(s) for %.1f MiB of assembly-group input. ",
		$assemblyCores, $assemblyInputMB;
	my $cleanSeqSetHR = sampleReadSet($curSmpl, "clean");
	my $hostFilter = 0;$hostFilter = 1 if ($AsGrps{$cAssGrp}{CntAimAss} > 3);#reset required HDD space
	my $tmpN ="";
	my $LasseP = $MFopt{DoAssembly};
	my ($primaryLibraries,$supportLibraries) = ([], []);
	
	if ($LasseP == 5){
		$primaryLibraries = getCleanLibrariesAssmGrp(\%AsGrps, $cAssGrp, 0);
		$supportLibraries = getCleanLibrariesAssmGrp(\%AsGrps, $cAssGrp, 1);
	}
	if ($LasseP == 5 && $AsGrps{$cAssGrp}{SupportReads} !~ /(?:PB|ONT):/){#$map{$curSmpl}{"SupportReads"} eq ""){#! scalar(@{$singlArX}) ){#decide on single tech
		$LasseP = 2; #go for megahit by default..
		my $technology = libraryTechnology($primaryLibraries,
			"assembly-group $cAssGrp primary assembly selection", 1);
		$LasseP  = 4 if ($technology eq "PB");
		$LasseP  = 3 if ($technology eq "ONT");
		#die "${$cReadTecAr}[0]\nXXXZ\n$LasseP\n";
	}
	# Preassemblies remain package-local intermediates. Complete ordinary,
	# long-read, and final hybrid assemblies publish to the canonical directory.
	my $assemblyOutDir = ($LasseP == 5 && $doPreAssmFlag)
		? $metagAssDir : $finalCommAssDir;
	my $publishedScaffDir = "$assemblyOutDir/scaffolds/";
	
	if ($LasseP == 5){ #hybrid assembly: first megahit(2), then metaMDBG(4)
		if ($doPreAssmFlag == 1){
			print "Preassembly step: ";
			$tmpN = megahitAssembly( \%AsGrps,$cAssGrp,$nodeTmp,$assemblyOutDir ,
				"preAssmbl", $SmplNameX,$hostFilter,$scaffoldFlag) if (!$ePreAssmbly);
		} elsif ($doPreAssmFlag == 0 && $postAssmblGo) {
			print "Final combining long assembly step: ";
			#$finalMapDir
			#die;
			$tmpN = longRdAssembly( \%AsGrps,$cAssGrp,"$nodeTmp",$assemblyOutDir,
				"hybridmMDBG",$SmplNameX,1,$LasseP) ; #$metaGpreAssmblDir, 
		} elsif ($doPreAssmFlag == 0 ){
			print "Preassembly: waiting for all samples\n";
		} else {
			die "Unknown preassembly state: $doPreAssmFlag\n";
		}
#					die;
	}elsif($LasseP == 1){
		$tmpN = spadesAssembly( \%AsGrps,$cAssGrp,"$nodeTmp",$assemblyOutDir,$MFopt{spadesBayHam} ,
			$MFglobal{shortAssembly}, $SmplNameX,$hostFilter,$scaffoldFlag) ;
	}elsif($LasseP == 2){
		$tmpN = megahitAssembly( \%AsGrps,$cAssGrp,"$nodeTmp",$assemblyOutDir ,
			$MFglobal{shortAssembly}, $SmplNameX,$hostFilter,$scaffoldFlag) ;
	} elsif( ($LasseP == 3 ||  $LasseP == 4)
		&& grep { $_->{is_long} } @{readLibrariesByScope($cleanSeqSetHR, 'primary', 1, $curSmpl)} ){
		$tmpN = longRdAssembly( \%AsGrps,$cAssGrp,"$nodeTmp",$assemblyOutDir,
		$MFglobal{shortAssembly}, $SmplNameX,0,$LasseP) ;
	}
	#die ;
	#if mates available, do them here
	append_job_dependencies(\$AsGrps{$cAssGrp}{AssemblJobName}, $tmpN); #always add in dep on read extraction
	
	
	# 2nd assembly step: scaffolding; maybe move later further down?
	if ($scaffoldFlag){
	#$finalCommScaffDir "$finalCommScaffDir/scaffDone.sto" $finAssLoc $metaGassembly
		#my $curAssLoc = $metaGassembly;
		#$curAssLoc = $finAssLoc if ($efinAssLoc);
		my ($newScaff,$sdep) = scaffoldCtgs(\%AsGrps,$cAssGrp,
				[],
				$assemblyOutDir,$nodeTmp."/scaff/",$publishedScaffDir,$AsGrps{$cAssGrp}{AssemblJobName},$MFopt{MapperCores}, $SmplNameX,1,"");
		unless ($sdep eq ""){
			append_job_dependencies(\$AsGrps{$cAssGrp}{AssemblJobName}, $sdep);
		}
	}
	
	#external contigs to be scaffolded (e.g. TEC2 extracts)
	if ($scaffTarExternal ne ""){
	#die "inscaff\n";
	#die "@scaffTarExternalOLib1\n";
		my $metaGscaffDirExt = "$metaGscaffDir/$scaffTarExternalName/";
		my ($newScaff,$sdep) = scaffoldCtgs(\%AsGrps,$cAssGrp,
				[],
				$scaffTarExternal,$nodeTmp."/SCFEX$scaffTarExternalName/",$metaGscaffDirExt,$AsGrps{$cAssGrp}{AssemblJobName},$MFopt{MapperCores}, 
				$SmplNameX,0,$scaffTarExternalName);
	#my($ar1,$ar2,$scaffolds,$GFdir_a) = @_; #.= "_GFI1";
		if (@scaffTarExternalOLib1 > 0 ){
			#die "in gapfill\n$newScaff\n";
			my $gapFillLibraries = readLibrariesFromArrays(
				sample => $SmplNameX, scope => 'primary', phase => 'external',
				technology => '', r1 => \@scaffTarExternalOLib1, r2 => \@scaffTarExternalOLib2,
			);
			GapFillCtgs($gapFillLibraries,$newScaff,$metaGscaffDirExt."GapFill/",$sdep,$scaffTarExternalName);
		}

		#last;
	}

	
	
	return;
} 


sub genePredictions($ $ $ $ $) {
	my ($inputScaff, $outDir, $jobDepend,$finDir,$specJname,$scrathD,$locSubm) = @_;
	my $tmpGene = $outDir."/tmpCalls/";
	my $numThr = 8; #defaul 4 cores..
	
	my $gzOut = ""; 
	
	$gzOut = ".gz" if ($MFopt{GenePredGZ});
	
	#system("mkdir -p $outDir");
	my $output_format_prodigal = "gff"; 
	my $expectedD = "$finDir/genePred/";
	my $augCmd  = ""; my $cmpCmd = "";
	my $tmpCmd=""; #placeholder for qsub return
	my $splitDep = $jobDepend;
	my $inputBac = $inputScaff; my $inputEuk = "$scrathD/euk.kraken.fasta";
	my $bacmark="";
	if ($MFopt{DoEukGenePred} ){ $bacmark = ".bac";	}
#	die;

	
	my $pprodigalBin = getProgPaths("pprodigal");
	#my $prodigalBin = getProgPaths("prodigal");
	
	#print "$outDir/proteins$bacmark.shrtHD.faa\n";
#	if ( (-s "$expectedD/proteins$bacmark.shrtHD.faa" && -s "$expectedD/genes$bacmark.gff") ||
#			(-s "$outDir/proteins$bacmark.shrtHD.faa" && -s "$outDir/genes$bacmark.gff") ){
	if ( ( fileGZs("$expectedD/proteins$bacmark.shrtHD.faa") ) 
			||(fileGZs("$outDir/proteins$bacmark.shrtHD.faa") ) 
			){
		#check if protein/genes were already gzip'd
		if ( $MFopt{GenePredGZ} && -e "$expectedD/genes$bacmark.gff" ){
			my $cmd = "";
			$cmd .= "rm -f $expectedD/proteins$bacmark.shrtHD.faa.gz ; $pigzBin -p $numThr $expectedD/proteins$bacmark.shrtHD.faa\n" if (-s "$expectedD/proteins$bacmark.shrtHD.faa");
			$cmd .= "rm -f  $expectedD/genes$bacmark.shrtHD.fna.gz; $pigzBin -p $numThr $expectedD/genes$bacmark.shrtHD.fna\n" if (-s "$expectedD/genes$bacmark.shrtHD.fna") ;
			$cmd .= "rm -f  $expectedD/genes$bacmark.per.ctg.gz; $pigzBin -p $numThr $expectedD/genes$bacmark.per.ctg\n" if (-s "$expectedD/genes$bacmark.per.ctg");
			$cmd .= "rm -f  $expectedD/genes$bacmark.gff.gz; $pigzBin -p $numThr $expectedD/genes$bacmark.gff\n" if (-s "$expectedD/genes$bacmark.gff");
			#$cmd .= "$pigzBin -p $numThr $expectedD/proteins$bacmark.shrtHD.faa $expectedD/genes$bacmark.shrtHD.fna $expectedD/genes$bacmark.per.ctg $expectedD/genes$bacmark.gff;"; 
			#my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0; 
			#my ($jname,$tmpCmd) = qsubSystem($logDir."prodigalGZ.sh",$cmd,$numThr,"2G","GZAA$JNUM","","",1,$QSBoptHR->{General_Hosts},$QSBoptHR); 
			#$QSBoptHR->{tmpSpace} =$tmpSHDD;
			#print $cmd;
			print "gzipping genePred/ files\n";
			system "$cmd"; #faster to just run locally..
		}
		return "";
	}
	#use kraken to classify contigs..
	if ($MFopt{DoEukGenePred} ){ #stupid way of doing things..
		#first split euk vs bac contigs :: 22.6.23:: deprecated: use whokaryote that requires prodigal annotations first
		die;
		my $splCores = 10;
		my ($refDB,$shrtDB,$clnCmd) = prepDiamondDB("NOG",$baseOut."DB/DiamDB/",$splCores,1);
		my $scmd = "mkdir -p $outDir\n";
		my $splitKgdContig = getProgPaths("contigKgdSplit_scr");
		$scmd .= "$splitKgdContig $inputScaff $scrathD $MFglobal{krakenDBDirGlobal} ".$baseOut."DB/DiamDB/ $splCores\n";
#		my $splitKgdContig2 = getProgPaths("contigKgdSplit_scr2");
#		$scmd .= "$splitKgdContig2 $inputScaff $scrathD $gff $splCores\n";
		$scmd .= "touch $scrathD/bac.euk.split.sto\n";
		#die $scmd."\n";
		$inputBac = "$scrathD/bact.kraken.fasta";
		#die "$inputEuk || !-e $inputBac\n";
		if (!-e $inputEuk || !-e $inputBac || !-e "$scrathD/bac.euk.split.sto"){
			($splitDep,$tmpCmd) = qsubSystem($logDir."splitAssembly.sh",$scmd,$splCores,int(50)."G","_SC$JNUM","$jobDepend;$krakDeps","",1,$QSBoptHR->{General_Hosts},$QSBoptHR); 
		}
	}

#die "toofar";
#run prodigal on bacterial contigs..
	my $gzCMD = "";
	if ($gzOut eq ".gz"){ $gzCMD = " | $pigzBin --stdout -p $numThr ";}
	my $prodigal_cmd = "";
	$prodigal_cmd .= "rm -rf $outDir;\nmkdir -p $tmpGene $outDir;\n";
	
	$prodigal_cmd.="if [ -e $inputBac ] && [ ! -s $inputBac ] ; then\n";
	$prodigal_cmd.="touch $tmpGene/proteins$bacmark.shrtHD.faa $tmpGene/genes$bacmark.per.ctg $tmpGene/genes$bacmark.shrtHD.fna $tmpGene/genes$bacmark.gff\n";
	$prodigal_cmd.="else\n";
	
	$prodigal_cmd .= ("$pprodigalBin -i $inputBac -o $tmpGene/genes$bacmark.gff -a $tmpGene/proteins.faa -d $tmpGene/genes.fna -f $output_format_prodigal -p meta -T $numThr \n");
	$prodigal_cmd .= "sleep 1;cut -f1 -d \" \" $tmpGene/proteins.faa $gzCMD > $tmpGene/proteins$bacmark.shrtHD.faa$gzOut\n" ;
	$prodigal_cmd .= "cut -f1 -d \" \" $tmpGene/genes.fna $gzCMD > $tmpGene/genes$bacmark.shrtHD.fna$gzOut\n" ;
	$prodigal_cmd .= "cut -f1 $tmpGene/genes$bacmark.gff | sort | uniq -c | grep -v '#' |  awk -v OFS='\\t' {'print \$2, \$1'} $gzCMD > $tmpGene/genes$bacmark.per.ctg$gzOut\n";
	$prodigal_cmd .= "$pigzBin -p $numThr $tmpGene/genes$bacmark.gff;\n" if ($gzOut eq ".gz");
	$prodigal_cmd .= "rm -f $tmpGene/proteins.faa $tmpGene/genes.fna\n";	
	#copy everything to the right place (for now, might have to change this later to take care of Eukarya)
	$prodigal_cmd .= "mkdir -p $expectedD\n";
		
	$prodigal_cmd.="fi\n";
	
	
	
	if ($MFopt{DoEukGenePred} ){ #inconvenient way of doing things..
		#1st set mark that prodigal genes are on bacteria only
		die "\$MFopt{DoEukGenePred} option needs some updates!\n";
		my $scrDir = "";#"/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/helpers/euk_gene_caller/bin";
		my $augustusBin = getProgPaths("augustus");#"/g/bork3/home/hildebra/bin/augustus-3.2.1/bin/augustus";
		$prodigal_cmd .= "touch $tmpGene/genes.prodigal.bact.only\n";
		$augCmd .= "mkdir -p $tmpGene\n";
		my $eukmark = ".euk";

		$augCmd.="if [ ! -s $inputEuk ];then\n";
		$augCmd.="touch $tmpGene/genes$eukmark.gff $tmpGene/genes$eukmark.codingseq $tmpGene/genes$eukmark.fna\n";
		$augCmd.="else\n";
		$augCmd .= "$augustusBin --species=ustilago_maydis $inputEuk --protein=off --codingseq=on > $tmpGene/augustus.gff;\n";
		$augCmd .= "perl $scrDir/getAnnoFast.pl --seqfile=$inputEuk $tmpGene/augustus.gff;";
		$augCmd.="mv $tmpGene/augustus.gff $tmpGene/genes$eukmark.gff\n mv $tmpGene/augustus.codingseq $tmpGene/genes$eukmark.codingseq\nmv $tmpGene/augustus.cdsexons $tmpGene/genes$eukmark.fna\n";
		$augCmd .= "cut -f1 -d \" \" $tmpGene/genes$eukmark.fna > $tmpGene/genes$eukmark.shrtHD.fna\n" ; 
		#TODO: gz needs to be added here..
		$augCmd .= "rm -f $tmpGene/genes$eukmark.fna\n";
		$augCmd.="fi\n";

		
		my $axl = "$tmpGene/augustus.list"; my $pxl = "$tmpGene/genes.list";
		my $cmpCmd = "grep -v \"^#\" $tmpGene/genes.gff > $pxl;\n";
		$cmpCmd .= "grep -w \"transcript\" $tmpGene/augustus.gff > $axl;\n";
		$cmpCmd .= "grep \">\" $inputEuk > $tmpGene.scaff.list;\n";
		$cmpCmd .= "awk \'!/^>/ { printf \"%s\", $0; n = \"\n\" } /^>/ { print n $0; n = \"\" } END { printf \"%s\", n }\' $tmpGene/augustus.codingseq > $axl.nolines.fa;\n";
		$cmpCmd .= "awk \'!/^>/ { printf \"%s\", $0; n = \"\n\" } /^>/ { print n $0; n = \"\" } END { printf \"%s\", n }\' $tmpGene/genes.shrtHD.fna > $pxl.nolines.fa;\n";
		$cmpCmd .= "perl $scrDir/process_nodes.pl $pxl $axl;\n";
		$cmpCmd .= "perl $scrDir/match_nodes_euk.pl $axl.output.txt $axl.nolines.fa;\n";
		$cmpCmd .= "perl $scrDir/match_nodes_bac.pl $pxl.output.txt $pxl.nolines.fa;\n";
		$cmpCmd .= "perl $scrDir/overlap_nodes.pl $tmpGene.scaff.list $axl.output.txt $pxl.list.output.txt ;\n"; #$f.blast.matched.txt
		$cmpCmd .= "perl $scrDir/finalgeneparsing.pl $tmpGene.scaff.list.ALLmatched.txt\n";
		#now, let's pull the appropriated gene. for each contig, need to retrieve ALL genes for that contig.
		$cmpCmd .= "cut -f 4 $tmpGene.scaff.list.ALLmatched.txt.decision.txt | sort |uniq> $tmpGene/formatching.txt\n";
		#$cmpCmd .= "perl ./bin/getgenes.pl $tmpGene/formatching.txt $f.aug.nolines.fa.renamed.fa \n";
		#$cmpCmd .= "perl ./bin/getgenes.pl $tmpGene/formatching.txt $f.mgm.nolines.fa.renamed.fa\n";
		#$cmpCmd .= "cat $f.aug.nolines.fa.renamed.fa.output.txt $f.mgm.nolines.fa.renamed.fa.output.txt > $f.genecalled.fa\n";
		#die $augCmd;
	}
	
	#finish up file transfer etc
	$prodigal_cmd .= "mv $tmpGene/* $outDir\nrm -fr $tmpGene\n";

	#if (!-e "$inputScaff") {die "Input scaffold $inputScaff does not exists.\n";}
	#die $prodigal_cmd."\n";
	#if (system $prodigal_cmd){die "Finished prodigal with errors";}
	my $jname = "";
	#print "$expectedD/proteins.shrtHD.faa\n$expectedD/genes.gff\n";
#	if ( (!-s "$expectedD/proteins$bacmark.shrtHD.faa" || !-s "$expectedD/genes$bacmark.gff") &&
#			(!-s "$outDir/proteins$bacmark.shrtHD.faa" || !-s "$outDir/genes$bacmark.gff") ){
	$jname = "_GP$JNUM"; 
	$jname = $specJname.$JNUM if ($specJname ne "");
	#print "genesubm\n";
	if ($locSubm){
		my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = 0; 
		($jname,$tmpCmd) = qsubSystem($logDir."prodigalrun.sh",$augCmd.$prodigal_cmd,$numThr,"2G",$jname,$splitDep,"",1,$QSBoptHR->{General_Hosts},$QSBoptHR); 
		$QSBoptHR->{tmpSpace} =$tmpSHDD;
	} else { #just return run command, can be exectuted from outside
		$jname = $augCmd.$prodigal_cmd."\n";
	}

#	}
	return $jname;
}




sub setupHPC{
	my $QSBoptHR1 = emptyQsubOpt($runOptions{submit},"",$MFconfig{submSytem});
	# Rejecting one job blocks its dependent chain without aborting unrelated
	# samples in this invocation.
	$QSBoptHR1->{continueOnSubmitError} = 1;
	$QSBoptHR1->{submissionErrors} = [];
	my $currentJobs = numUserJobs($QSBoptHR1,1);
	print "Found $currentJobs jobs registered to user.";
	#could be only the submitting job is active? is this within a submission?
	if ($currentJobs == 0 && $MFconfig{rmSmplLocks} ==0){print " Auto-removing sample locks (to override set \"-rmSmplLocks -1\").\n";$MFconfig{rmSmplLocks}=1;}
	elsif ($MFconfig{rmSmplLocks} == -1 ){$MFconfig{rmSmplLocks}=0;} #clearly a debug option.. not documented
	else {print "\n";}
	#die;
	#my %QSBopt = %{$QSBoptHR1};
	$QSBoptHR1->{tmpSpace} = $MFconfig{nodeHDDspace};
	$QSBoptHR1->{wcKeysForJob} = $MFconfig{wcKeysForJob};
	$QSBoptHR1->{excludeNodes} = $MFconfig{excludeNodes};
	$QSBoptHR1->{jobPollSeconds} = $MFconfig{schedulerPollSeconds};
	$QSBoptHR1->{maxConcurrentJobs} = $MFconfig{checkMaxNumJobs};
	$QSBoptHR1->{killDependencyNever} = $MFconfig{killDepNever};
	# Queue-size checks used to run once per sample. Cache a successful query,
	# while counting every locally submitted job conservatively against it.
	$QSBoptHR1->{pendingJobCheckInterval} =
		$MFconfig{schedulerPollSeconds} * 3;
	$QSBoptHR1->{pendingJobCheckInterval} = 30
		if $QSBoptHR1->{pendingJobCheckInterval} < 30;
	$QSBoptHR1->{Spades_Hosts} = []; $QSBoptHR1->{General_Hosts} = [];
	spadesHosts();
	#my $LocationCheckStrg=""; #command that is put in front of every qsub, to check if drives are connected, sub checkDrives
	#queing capability
	return $QSBoptHR1
}



sub setDefaultMFconfig{
	
	print "Setting default parameters..  ";
	
	#MFglobal
	$MFglobal{shortAssembly} = ""; $MFglobal{prevAssembly} = ""; 
	
	#progStats
	#progStats: object to track progress of programs/submissions
	$progStats{riboFindFailCnts}=0; $progStats{riboFindComplCnts} = 0; 
	$progStats{taxTarFailCnts}=0;  $progStats{taxTarComplCnts}=0; 
	$progStats{mOTU2FailCnts}=0;  $progStats{mOTU2ComplCnts}=0; 
	$progStats{metaPhl2FailCnts}=0; $progStats{metaPhl2ComplCnts}=0;
	$progStats{KrakTaxFailCnts} =0;


	#HDDspace: object to handle HDD usage: Always format as "XXG" XX = space requirements in Gb. Excecption: "-1"
	$HDDspace{kraken} = "120G"; 
	$HDDspace{assembler} = "-1";
	$HDDspace{spades} = "250G"; 
	$HDDspace{megaHit} = "140G";
	$HDDspace{riboFind} = "100G";
	$HDDspace{metaMDBG} = "150G";
	$HDDspace{SNPcall} = "60G";
	$HDDspace{metabat2} = "120G";
	$HDDspace{mapping} = "4G"; #scales linear per input 1 GB filesize
	$HDDspace{diamond} = "80G";
	$HDDspace{prepPub} = "80G";
	$HDDspace{Ribos} = "80G";


	#MFcontstants: object to store essential paths/file endings
	$MFcontstants{bwt2IdxFileSuffix} = ".bw2";
	$MFcontstants{mini2IdxFileSuffix} = ".mmi";
	$MFcontstants{kmaIdxFileSuffix} = ".kma";
	#locking samples during processing (so no other jobs are started in these)
	$MFcontstants{DefaultSampleLock} = "MGTK.locked";



	#MFopt: global object with options for MATAFILER. Added in MF v0.5, slowly rebuild MF around this system

	#non-asembly based tax + functional assignments #DoMetaPhlan3 merges into $DoMetaPhlan , after version check
	#my $DoMetaPhlan = 0;   my $DoMetaPhlan3 = 0; my $DoMOTU2 = 0; my $DoTaxaTarget = 0; 
	$MFopt{DoMetaPhlan}=0;#$MFopt{DoMetaPhlan3}=0;
	$MFopt{DoMOTU2}=0;$MFopt{DoTaxaTarget}=0;
	$MFopt{PABtaxChk} =0;
	$MFconfig{skipSmallSmplsMB} = 1; #skips samples with less MB in input  files


	#my $DoRibofind = 0; my $doRiboAssembl = 0; my $RedoRiboFind = 0; #ITS/SSU/LSU detection
	$MFopt{DoRibofind}=0; $MFopt{doRiboAssembl}=0; $MFopt{RedoRiboFind}=0;#ITS/SSU/LSU detection
	#my $riboLCAmaxRds = 250000; #set to 50k by default, larger gets slower too fast; MF0.45: 250k due to usage of lambda3
	$MFopt{riboLCAmaxRds} = 250000;
	$MFopt{riboStoreRds}=0;$MFopt{RedoRiboAssign}=0;$MFopt{RedoRiboThatFailed}=0;$MFopt{checkRiboNonEmpty}=0;
	$MFopt{globalRiboDependence} = {DBcp => ""};

	#Binning related options
	$MFopt{useBinnerScratch} = 0;$MFopt{BinnerMem} = 0;$MFopt{useCheckM2} = 1;
	$MFopt{DoBinning} = 0;$MFopt{useCheckM1} = 0;$MFopt{BinnerCores} = 9;
	$MFopt{DoMetaBat2} = 0; $MFopt{BinnerRedoEmpty} = 0;  $MFopt{SB_env} = "";
	$MFopt{BinnerRedoAll} =0; $MFopt{minBinnerAssemblyMB} = 2;

	#read preprocessing
	$MFopt{unzipCores} = 3; 
	$MFopt{trimAdapters} =1;
	$MFopt{usePorechop} = 0;
	$MFopt{SDMlogQualvsLen} = 0; #sdm log of qual per read vs length (eg for PacBio qual checks..)
	$MFopt{sdmCores} = 6; #sdm specific cores  #currently set to 1, sdm multi thread instability
	$MFopt{sdmMem} = "15G"; #total mem for sdm job in Gb, default 15
	$MFopt{sdm_opt} = {}; #empty object that can be used to modify default sdm parameters
	$MFopt{tmpSdmminSL} =0; $MFopt{tmpSdmmaxSL}=0;
	$MFopt{gzipSDMOut} = 1;#zip sdm filtered files
	$MFopt{sdmProbabilisticFilter} =1;

	#which read filtering option to use in sdm?
	$MFopt{useSDM} = 2;	
	$MFopt{sdmOpt} = "";
	$MFopt{baseSDMopt} = getProgPaths("baseSDMopt_rel"); 
	if ($MFopt{useSDM} ==2 ){$MFopt{baseSDMopt} = getProgPaths("baseSDMopt");}
	$MFopt{baseSDMoptMiSeq} = getProgPaths("baseSDMoptMiSeq_rel");
	if ($MFopt{useSDM} ==2 ){$MFopt{baseSDMoptMiSeq} = getProgPaths("baseSDMoptMiSeq");}


	$MFopt{writeStats} = 0;


	#Assembly related options
	$MFopt{doReadMerge} = 0;
	$MFopt{DoAssembly} = 0;  #1=Spades, 2=MegaHIT, 3= flye, 4=metaMDBG, 5=hybrid ill-PB (megahit, metaMDBG), 0=deactivated
	$MFopt{SpadesAlwaysHDDnode} = 1;$MFopt{spadesBayHam} = 0; 
	$MFopt{spadesMisMatCor} = 0; $MFopt{redoAssembly} =0 ; $MFopt{SpadesLongtime} = 0;
	$MFopt{pseudoAssembly} = 0; #in case no assembly is possible (soil single reads), just filter for reads X long 
	$MFopt{AssemblyCores} = 0; $MFopt{AssemblyMemory} = -1; #0 auto-scales from 8 to 48; memory in GB
	$MFopt{AssemblyKmers} = "27,43,67,87,101,127" ; #"27,33,55,71";
	$MFopt{kmerPerGene} = 0; #calculate kmer frequencies for each gene instead of per scaffold
	$MFopt{kmerAssembly} = 0; #calculate kmer frequencies for each scaffold
	$MFopt{scaffoldMinSize} = 500; #all scaffolds/contigs below this will be dropped

	#mapping related options
	$MFopt{MapperProg} = -1;#1=bowtie2, 2=bwa, 3=minimap2, 4=kma, 5=strobealign, -1=auto (bowtie2 short, minimap2 long reads)
	$MFopt{map2Assembly} = 1; $MFopt{mapSortMemGb} = -1; #in Gb
	$MFopt{SaveUnalignedReads} =0; $MFopt{useUnmapped} = 0;
	$MFopt{mapSupport2Assembly} = 1;
	$MFopt{bwtIdxAssMem} = 40; #total mem in GB, not core adjusted, for building index from assembly
	$MFopt{doBam2Cram} = 1; $MFopt{redoAssMapping} =0;
	$MFopt{DoJGIcoverage} = 0; #only required for metabat binning, not required any longer..
	$MFopt{bamfilterIll} = "0.05 0.75 20 3"; $MFopt{bamfilterPB} = "0.05 0.5 30 0"; $MFopt{bamfilterONT} = "0.15 0.5 10 0";
	# Hybrid preassembly coverage is deliberately conservative: alignments must
	# cover more of the read and have substantially higher mapping confidence.
	$MFopt{bamfilterHybridIll} = "0.03 0.90 40 5";
	$MFopt{hybridMinMapQ} = 40; $MFopt{hybridMinBaseQ} = 20;
	$MFopt{breakpointDepth} = 0.10; $MFopt{breakpointMinLength} = 100;
	$MFopt{breakpointSmoothGap} = 100; $MFopt{breakpointFlankLength} = 500;
	$MFopt{breakpointMinFlankDepth} = 1; $MFopt{breakpointMaxFlankFraction} = 0.10;
	$MFopt{hybridSyntheticMaxDepth} = 20;
	$MFopt{mapSaveCram} = 1; #by default, keep the back-mapping bams/crams 
	$MFopt{MapperCores} = 8;  $MFopt{MapperRmDup} = 1; #mapping cores; ??? ; remove Dups (can be costly if many ref seqs present)
	$MFopt{bamSortCores} = -1;
	$MFopt{MapperMemory} = -1; #total mem for mapping job in Gb, Default -1
	$MFopt{redoMapping} = 0; #rewrite mapping?
	$MFopt{mapModeActive}=0; $MFopt{mapModeCovDo}=1;#get the coverage per gene etc
	$MFopt{mapModeTogether} = -1; #-1: comb map & report, 0: separate bwt runs, 1:competitive, 2: non-competitive mapping, sep reports per file
	$MFopt{largeMapperDB} = 0; #use flags in mapper index built for large ref DBs?

	$MFopt{DoMapModeDecoy} =1;  $MFopt{MapRewrite2nd} = 0;
	$MFopt{Do2ndMapSNP} = 0;
	$MFopt{refDBall} =""; $MFopt{bwt2NameAll}="" ; #user options for refDBs

	#extra programs to be used?
	$MFopt{DoNonPareil}=0; #non-pareil
	$MFopt{DoCalcD2s} = 0;  #D2s distance .. bit outdated, no longer used..
	$MFopt{calcOrthoPlacement} =0; #Jaime's 6 frame translation and hmmersearch (AB production)
	$MFopt{DoGenoSizeEst} =0; #genome size estimate via microcensus

	#Kraken related config
	$MFopt{DoKraken} = 0; $MFopt{RedoKraken} = 0; 
	$MFopt{krakenCores} = 9;
	$MFopt{completeContaStats} = 1;
	$MFopt{humanFilter} = 0; ##0: no, 1:kraken2, 2: kraken1, 3:hostile
	$MFopt{filterHostDB1} = "";  #customize host org to filter (e.g. human, chicken ..)
	$MFopt{krakHostConf} = 0.01; #confidence needed to assign read to host ref
	$MFopt{filterHostKr2QuickMode}{0} = "";$MFopt{filterHostKr2QuickMode}{1} = "";# 0/1 for is3rdGen? "--quick "; deactivated for now..
	$MFopt{hostileIndex} = "human-t2t-hla";
	$MFopt{globalKraTaxkDB} = "";
	$MFopt{globalDiamondDependence} = {CZy=>"",MOH2 => "", MOH=>"",NOG=>"",ABR=>"",ABRc=>"",KGB=>"",KGE=>"",ACL=>"",KGM=>"", PTV=>"", PAB => "", URE=>"", URacc=>"", AMI=>""};
	



	#SNPs
	$MFopt{DoConsSNP}=0; $MFopt{DoSuppConsSNP}=0; $MFopt{redoSNPcons} = 0; $MFopt{redoSNPgene} =0; $MFopt{SNPconsJobsPsmpl} = 0;
	$MFopt{SNPminCallQual} = 20; $MFopt{memPJob} = 0; #set to 0 to indicate default estimation
	$MFopt{saveVCF} = 1; $MFopt{saveConsFastas} = 0;
    #$MFopt{memSNPcall} = 23; -> no longer used
	$MFopt{maxSNPcores} = 10;  $MFopt{consSNPminDepth} = 0; $MFopt{normSNPindels} = 1;
	$MFopt{SNPcallerFlag} = "MPI"; #"MPI" mpileup or ".FB" for freebayes
	$MFopt{callSVs} = 0; #0=not, 1=delly, 2=gridss
	$MFopt{callSVsSupp} = 0; #same as "callSVs" but for supplemental reads
	$MFopt{SVcallerFlag} = "DL"; #DL for delly, or GR for gridss


	#Func annotation
	$MFopt{DoDiamond} = 0; $MFopt{rewriteDiamond} =0; $MFopt{redoDiamondParse} = 0; #redoes matching of reads; redoes interpretation
	$MFopt{rewriteAllIfAnyDiamond}=0;
	$MFopt{maxReqDiaDB} = 6; #max number of databases supported by MATAFILER
	$MFopt{reqDiaDB} = "";#,NOG,MOH,ABR,ABRc,ACL,KGM,PTV,PAB";#,ACL,KGM,ABRc,CZy";#"NOG,CZy"; #"NOG,MOH,CZy,ABR,ABRc,ACL,KGM"   #old KGE,KGB
	$MFopt{diaEVal} = "1e-7"; $MFopt{diaCores} = 12; ; $MFopt{DiaRmRawHits} = 0; $MFopt{diaRunSensitive} = 0;
	$MFopt{diaFrameshift} = 0; #diamond -F/--frameshift penalty; 0 disables frameshift alignment mode (recommended for long, error-prone reads e.g. ONT/PacBio)
	$MFopt{diamondMem} = 7; #GB memory to request for diamond jobs (qsub --mem)
	$MFopt{DiaMinAlignLen} = 20; $MFopt{DiaMinFracQueryCov} = 0.1; $MFopt{DiaPercID} =40;

	#gene prediction related
	$MFopt{DoEukGenePred} = 0;
	$MFopt{GenePredGZ} = 1; #gzip output to reduce storage usage
	$MFopt{genePredGZenforce}=1; #recheck and rezip if not already done?
	$MFopt{rewriteGenePred} = 0; #rewrite gene predictions and the few stats related?
	
	
	#MFconfig configuration with defaults
	$MFconfig{mapFile} = "";
	$MFconfig{inspectState} = 0;
	$MFconfig{planState} = 0;
	$MFconfig{autoStatePlan} = 0;
	$MFconfig{autoRepairState} = 1;
	$MFconfig{stateReport} = "";
	$MFconfig{planReport} = "";
	$MFconfig{configFile} = "";
	$MFconfig{nodeHDDspace} = 30; #30 Gb; default value
	$MFconfig{maxUnzpJobs} = 20; #how many unzip jobs to run in parallel (not to overload HPC IO)

	$MFconfig{rawFileSrchStr1} = '.*1\.f[^\.]*q\.gz$';
	$MFconfig{rawFileSrchStr2} = '.*2\.f[^\.]*q\.gz$';
	$MFconfig{rawFileBamSrchSing} = ""; #inputBAMregex 
	$MFconfig{rawFileSrchStrSingl} = "";
	$MFconfig{rawFileSrchStrXtra1} = '.*1_sequence\.f[^\.]*q\.gz$';
	$MFconfig{rawFileSrchStrXtra2} = '.*2_sequence\.f[^\.]*q\.gz$';
	$MFconfig{submSytem} = ""; #user supplied submission system flag.. default is empty and autodetect


	$MFconfig{mateInsertLength} = 20000; #controls expected mate insert size , import for bowtie2 mappings

	#more specific control: unfiniRew=rewrite unfinished sample dir; $redoCS = redo ContigStats completely; 
	#removeInputAgain=remove unzipped files from scratch, after sdm; remove_reads_tmpDir = leave cleaned reads on scratch after everything finishes
	$MFconfig{unfiniRew}=0; $MFconfig{redoCS}=0; $MFconfig{removeInputAgain}=1; $MFconfig{remove_reads_tmpDir}=1; 
	$MFconfig{OKtoRWassGrps} = 0; $MFconfig{rmBinFailAssmbly} = 0;
	$MFconfig{skipWrongPairedSmpls} =1; #check in sdm output if wrong read pairs present
	$MFconfig{defaultContigSubs} = "as"; #default subprts for contigstats..

	$MFconfig{silent} = 0;
	$MFconfig{redoFails} = 0;$MFconfig{XfirstReads} = -1;
	$MFconfig{killDepNever} = 0;  $MFconfig{checkMaxNumJobs} = 0;  #slurm related.. $killDepNever=1 kills jobs in state "DependencyNeverFinished" (happens a lot), while $checkMaxNumJobs=X halts the pipeline if more than X jobs are already queued up
	$MFconfig{loopTillCompleteActiveJobs} = 3;
	$MFconfig{schedulerPollSeconds} = 20;
	$MFconfig{excludeNodes} = ""; #excluding certain nodes..
	$MFconfig{readsRpairs} =-1; #are reads given in pairs? default: -1 = no clue
	#my $useTrimomatic=0; #trimmomatic step now replaced by sdm solution -> $MFopt{trimAdapters}

	$MFconfig{abortOnEmptyInput} = 0; 
	#my $relaxedSmplNames = 0;#don't abort when SMPLID in map contains strange chars  -> this is now handled by adding to map "#RelaxSMPLID	TRUE"
	$MFconfig{ignoreSmpl} = "";
	$MFconfig{splitFastaInput} = 0; #assembly as input..
	$MFconfig{importMocat} = 0; $MFconfig{mocatFiltPath} = "reads.screened.screened.adapter.on.hg19.solexaqa/"; 
	$MFconfig{alwaysDoStats} = 1; 
	$MFconfig{rmScratchTmp}=1;#Default; extremely important option as this adds a lot of overhead and scratch usage space, but reduces later overhead a lot and makes IO more stable
	$MFconfig{unpackZip} = 0; #only goes through with read filtering, needed to get files for luis
	$MFconfig{filterFromSource}=1; #powerful option that skips the unzip step.. use careful
	$MFconfig{doDateFileCheck} = 0; #very specific option for Moh's reads that were of different dates..
	$MFconfig{DoFreeGlbTmp} = 0; 
	$MFconfig{defaultReadLength} = 150; $MFconfig{defaultReadLengthX} = 5000;
	$MFconfig{oldStylFolders} =0; #0=smpl name as out folder; 1=inputdir as out foler (legacy)
	$MFconfig{mocatLinkDir} = "";
	$MFconfig{wcKeysForJob} = ""; #EI specific system to register jobs under certain flag
	$MFconfig{prefSinglFQgreps} = 0; #if grep of files (rawSrchString) has multi assignments, which grep to trust more?
	$MFconfig{rmSmplLocks} = 0;
	$MFconfig{uploadRawRds} = ""; #prepare raw input fastq's for upload 2 EBI? Clean reads will be stored in this dir


#statistics collections etc, things that can count/change during run

	
	print "Done. ";

}



sub help {
	my $scriptPath = abs_path(__FILE__) || __FILE__;
	my $reference = dirname($scriptPath)."/docs/flag_reference.md";
	my @sections;
	my $currentSection;

	if (open(my $referenceFH, '<', $reference)) {
		my $inMataf4Reference = 0;
		while (my $line = <$referenceFH>) {
			if ($line =~ /^##\s+MATAF4\.pl\s*$/i) {
				$inMataf4Reference = 1;
				next;
			}
			last if ($inMataf4Reference
				&& $line =~ /^##\s+Flag comparison against previous manual\.md\s*$/i);
			next unless ($inMataf4Reference);

			if ($line =~ /^##\s+(.+?)\s*$/) {
				$currentSection = { title => $1, options => [] };
				push @sections, $currentSection;
				next;
			}
			next unless ($line =~ /^\|\s*`-/);

			my @cells = split(/\|/, $line, -1);
			next unless (@cells >= 7 && defined($currentSection));
			my ($aliases, $type, $default, $status, $description)
				= map { _plainHelpCell($_) } @cells[1..5];
			push @{$currentSection->{options}}, {
				aliases => $aliases,
				type => $type,
				default => $default,
				status => $status,
				description => $description,
			};
		}
		close($referenceFH);
	}

	print "\nMATAFILER4 command-line help (version $MATFILER_ver)\n";
	print "Usage: MATAF4.pl -map <mapping-file> [options]\n";
	print "       MATAF4.pl -help | -h | -?\n\n";
	print "Options may be written with one or two leading dashes.\n";

	if (@sections && grep { @{$_->{options}} } @sections) {
		print "The complete option list below is read from docs/flag_reference.md.\n";
		foreach my $section (@sections) {
			next unless (@{$section->{options}});
			print "\n$section->{title}:\n";
			foreach my $option (@{$section->{options}}) {
				_printHelpOption($option);
			}
		}
		print "\nFull reference: $reference\n";
	} else {
		print "\nUnable to read the MATAF4.pl option tables from:\n  $reference\n";
		print "See that file for the complete flag reference.\n";
	}
	print "\n";
	exit(0);
}

sub _plainHelpCell {
	my ($text) = @_;
	$text = "" unless (defined($text));
	$text =~ s/^\s+|\s+$//g;
	$text =~ s/`//g;
	$text =~ s/\*\*//g;
	$text =~ s/\[([^\]]+)\]\([^\)]+\)/$1/g;
	return $text;
}

sub _printHelpOption {
	my ($option) = @_;
	my $signature = $option->{aliases};
	$signature .= " <$option->{type}>"
		if ($option->{type} ne "" && lc($option->{type}) ne "flag");

	my $details = $option->{description};
	my @metadata;
	push @metadata, "default: $option->{default}" if ($option->{default} ne "");
	push @metadata, "status: $option->{status}"
		if ($option->{status} ne "" && lc($option->{status}) ne "stable");
	$details .= " " if ($details ne "" && @metadata);
	$details .= "(".join("; ", @metadata).")" if (@metadata);

	local $Text::Wrap::columns = 100;
	local $Text::Wrap::huge = 'overflow';
	print "  $signature\n";
	print wrap("      ", "      ", $details)."\n" if ($details ne "");
}



sub getCmdLineOptions{
	
	print "Reading command line options.. ";
	
	GetOptions(

	#base options
		"help|?|h" => \&help,
		"checkInstall" => sub { checkMFFInstall("",1) },
		"map=s"      => \$MFconfig{mapFile},
		"config=s" => \$MFconfig{configFile},
		"inspectState=i" => \$MFconfig{inspectState},
		"planState=i" => \$MFconfig{planState},
		"stateReport=s" => \$MFconfig{stateReport},
		"planReport=s" => \$MFconfig{planReport},
		"autoStatePlan=i" => \$MFconfig{autoStatePlan},
		"autoRepairState=i" => \$MFconfig{autoRepairState},

	#flow related
		"redoFails=i" =>\$MFconfig{redoFails}, #if any step of requested analysis failed, just redo everything (extraction etc)
		"redoContigStats=i" => \$MFconfig{redoCS}, #runContigStats (coverage per gene, kmers, GC content) will be deleted & started again
		"submSystem=s" => \$MFconfig{submSytem},  #qsub,SGE,bsub,LSF.. by default will try to autodetect
		"submit=i" => \$runOptions{submit},  #submit any jobs at all? (0= no submission, just for trying if everything is correctly set up)
		"from=i" => \$runOptions{from},  #start at which samples from map file?
		"to=i" => \$runOptions{to},   #stop at which samples from map file?
		"loopTillComplete=s" => \$runOptions{loopCount}, #rolling completion loop followed by one final full-range verification pass
		#use syntax "X:Y" where X is the pass budget and Y is the active rolling window size
		"loopTillCompleteActiveJobs=i" => \$MFconfig{loopTillCompleteActiveJobs}, #start the next pass once no more than this many submitted jobs are running
		"schedulerPollSeconds=i" => \$MFconfig{schedulerPollSeconds}, #seconds between loopTillComplete scheduler queries
		"excludeNodes=s" => \$MFconfig{excludeNodes}, #exclude certain nodes?
		"maxConcurrentJobs=i" => \$MFconfig{checkMaxNumJobs}, #max running plus pending user jobs, enforced before every Slurm submission
		"killDepNever=i" => \$MFconfig{killDepNever}, #kill jobs in "Dependency never finished" state? 
		"requireInput=i" => \$MFconfig{abortOnEmptyInput},  #in case input reads are no longer present, 0 will continue pipeline, 1 will abort
		"ignoreSmpls=s" => \$MFconfig{ignoreSmpl},  #skip a certain sample (sample id)
		"rmSmplLocks=i" => \$MFconfig{rmSmplLocks}, #remove existing sample locks (useful if jobs have crashed, leaving abondened sample locks)
		#"rmRawRds=i" => \$MFconfig{DoFreeGlbTmp}, #rm raw sequences, once all jobs have finished <-- redundant with reduceScratchUse
		"silent" => \$MFconfig{silent},
		"maxUnzpJobs=i" => \$MFconfig{maxUnzpJobs}, #how many unzip jobs to run in parallel (not to overload HPC IO). Default:20
		"skipSmallSmplsMB=i" => \$MFconfig{skipSmallSmplsMB},  #skip samples with a combined input smaller than this in MB (raw file size, independent of compressed or raw)
		"forceWriteStats=i" => \$MFopt{writeStats}, # force (re)writing of the metagStats report and text file
	
	#input FQ related
	#file strucuture
		#"relaxedSmplNames=i" => \$relaxedSmplNames, #add instead to map: #RelaxSMPLID	TRUE
		"rm_tmpdir_reads=i" => \$MFconfig{remove_reads_tmpDir}, #Default 1, remove tmpdir with reads
		"rm_tmpInput=i" => \$MFconfig{removeInputAgain},#remove raw, human / adaptor filtered reads, if sdm clean created? (and not needed any longer)
		"reduceScratchUse=i" => \$MFconfig{rmScratchTmp}, #should always be 1, unless debugging..
		"globalTmpDir=s" => \$runOptions{sharedTmpDir},#absolute path to global shared tmp dir (like a scratch dir)
		"nodeTmpDir=s" => \$runOptions{nodeTmpDir},#absolute path to tmp dir on local HDD of each executing node
		"nodeHDDspace=s" => \$MFconfig{nodeHDDspace},#HDD tmp space to be requested for each node (in Gb). Some systems don't support this
		"legacyFolders=i" => \$MFconfig{oldStylFolders}, #legacy option, controls if output folders will use the read dir as name (1) or the name in the mapping file (0). Default=0

	#preprocessing (cleaning reads etc)
	#	"useTrimomatic=i" => \$useTrimomatic, #rm adapter seq from input reads
		"usePorechop=i" => \$MFopt{usePorechop}, #adapter rm for Nanopore.. should actually be automatically with newer sdm (not implemented)
		"inputFQregex1=s" => \$MFconfig{rawFileSrchStr1}, #regex for detecting read pair 1 in input fastq files
		"inputFQregex2=s" => \$MFconfig{rawFileSrchStr2}, #regex for detecting read pair 2 in input fastq files
		"inputFQregexSingle=s" => \$MFconfig{rawFileSrchStrSingl}, #regex for detecting single end reads in input fastq files
		"inputFQregexTrustSingle=i" => \$MFconfig{prefSinglFQgreps} , #if grep of files (rawSrchString) has multi assignments, which grep to trust more?
		"inputBAMregex=s" => \$MFconfig{rawFileBamSrchSing}, #bams that will be converted to fastq
		"splitFastaInput=i" => \$MFconfig{splitFastaInput},
		"mergeReads=i" => \$MFopt{doReadMerge},  #merge read pair 1+2 before assembly etc? (usually doesn't help assembly, but useful for mapping to ref database in some rare instances)
		"ProbRdFilter=i" => \$MFopt{sdmProbabilisticFilter},
		"pairedReadInput=i" => \$MFconfig{readsRpairs}, #determines if read pairs are expected in each in dir
		"inputReadLengthSuppl=i" => \$MFconfig{defaultReadLengthX},
		"filterHostRds|filterHumanRds=i" => \$MFopt{humanFilter}, #0: no, 1:kraken2, 2: kraken1, 3:hostile
		"filterHostKrak2DB=s" => \ $MFopt{filterHostDB1}, #customize host org to filter (e.g. human, chicken ..)
		"filterHostKr2Conf=s" => \$MFopt{krakHostConf},
		"filterHostKr2Quick=s" => \$MFopt{filterHostKr2QuickMode}{0},
		"hostileIndex=s" => \$MFopt{hostileIndex},
		"onlyFilterZip=i" => \$MFconfig{unpackZip},
		"mocatFiltered=i" => \$MFconfig{importMocat},
		"logQualvsLen=i" => \$MFopt{SDMlogQualvsLen}, #sdm log file.. can be quite large; logs qual of read vs read length

	#sdm related
		"inputReadLength=i" => \$MFconfig{defaultReadLength},
		"gzipSDMout=i" => \$MFopt{gzipSDMOut},
		"XfirstReads=i" => \$MFconfig{XfirstReads},
		"minReadLength=i" => \$MFopt{tmpSdmminSL},
		"maxReadLength=i" => \$MFopt{tmpSdmmaxSL},
		"filterAdapters=i" => \$MFopt{trimAdapters},
		"customSDMopt=s"  => \$MFopt{sdmOpt},
		"sdmMem=s" => \$MFopt{sdmMem}, #total mem for sdm job in Gb, default 15

	#assembly related
		"spadesCores|assemblCores=i" => \$MFopt{AssemblyCores},
		"spadesMemory|assemblMemory=i" => \$MFopt{AssemblyMemory}, #in GB
		"spadesKmers|assemblyKmers=s" => \$MFopt{AssemblyKmers}, #comma delimited list
		"reAssembleMG=i" => \$MFopt{redoAssembly},
		"asssemblyHddSpace=i" => \$HDDspace{assembler},
		"assembleMG=i" => \$MFopt{DoAssembly}, #1=Spades, 2=MegaHIT, 3= flye, 4=metaMDBG, 5=hybrid ill-PB (megahit, metaMDBG)
		"assemblyLongTime=i" => \$MFopt{SpadesLongtime},
		"assemblyScaffMinSize=i" => \$MFopt{scaffoldMinSize},
	#binning
		"Binner|MetaBat2|binSpeciesMG=i" => \$MFopt{DoMetaBat2}, #0=no, 1=metaBat2, 2=SemiBin, 3: MetaDecoder, 4: GenomeFace, 5: SCGBinner
		"BinnerCores=i" => \$MFopt{BinnerCores}, #cores used for Binning process (and checkM)
		"BinnerMem=i" => \$MFopt{BinnerMem}, # define binning memory, Gb, 0=auto
		"minBinnerAssemblyMB=f" => \$MFopt{minBinnerAssemblyMB}, #skip binning assemblies below this many million sequence bases; 0 disables
		"checkM2=i" => \$MFopt{useCheckM2},
		"checkM1=i" => \$MFopt{useCheckM1},
		"BinnerScratchTmp=i" => \$MFopt{useBinnerScratch}, #very specific (undocumented) use of scratch instead of nodetmp dir
		#"binSpeciesMG=i" => \$MFopt{DoBinning},#deactivated, replaced by MetaBat2
		"redoEmptyBins=i" => \$MFopt{BinnerRedoEmpty}, #debug option; redo bins that are empty (no bin detected). Note: this can sometimes happen for metagenomes
		"redoBinning=i" => \$MFopt{BinnerRedoAll}, 
		"SB_env=s"  => \$MFopt{SB_env}, #semiBin environment; if given, will avoid re-training de novo binning model. Default: "" (autotrain). should be #human_gut/dog_gut/ocean/soil/cat_gut/human_oral/mouse_gut/pig_gut/built_environment/wastewater/chicken_caecum/global
	#gene prediction on assembly
		"predictEukGenes=i" => \$MFopt{DoEukGenePred},#severely limits total predicted gene amount (~25% of total genes)
		"kmerPerGene=i" => \$MFopt{kmerPerGene}, #calculate kmer frequencies for each gene instead of per scaffold
		"genePredGZenforce=i" => \$MFopt{genePredGZenforce},
		"rewriteGenePred=i"   => \$MFopt{rewriteGenePred},
	#mapping
		"mapper=i" => \$MFopt{MapperProg}, ##1=bowtie2, 2=bwa, 3=minimap2, 4=kma, 5=strobealign -1=auto (bowtie2 short, minimap2 long reads), -2=auto(strobealign short, minimap2 long) 
		"mapUnmapped=i" => \$MFopt{useUnmapped},
		"mappingCoverage=i" => \$MFopt{mapModeCovDo},
		"mappingMem=i" => \$MFopt{MapperMemory}, #total mem for mini2/kma/bwa/bwt2 in GB
		"mapSortMem=i" => \$MFopt{mapSortMemGb}, #total mem for samtools sort in GB
		"rmDuplicates=i" => \$MFopt{MapperRmDup},
		"mappingCores=i" => \$MFopt{MapperCores},
		"mapperFilterIll=s" => \$MFopt{bamfilterIll}, # max NM edit rate, min query coverage, min mapping quality, min clip at both ends (0 disables); default: "0.05 0.75 20 3"
		"mapperFilterHybridIll=s" => \$MFopt{bamfilterHybridIll},
		"hybridMinMapQ=i" => \$MFopt{hybridMinMapQ},
		"hybridMinBaseQ=i" => \$MFopt{hybridMinBaseQ},
		"breakpointDepth=f" => \$MFopt{breakpointDepth},
		"breakpointMinLength=i" => \$MFopt{breakpointMinLength},
		"breakpointSmoothGap=i" => \$MFopt{breakpointSmoothGap},
		"breakpointFlankLength=i" => \$MFopt{breakpointFlankLength},
		"breakpointMinFlankDepth=f" => \$MFopt{breakpointMinFlankDepth},
		"breakpointMaxFlankFraction=f" => \$MFopt{breakpointMaxFlankFraction},
		"hybridSyntheticMaxDepth=f" => \$MFopt{hybridSyntheticMaxDepth},
		"mapperFilterPB=s" => \$MFopt{bamfilterPB},
		"mapperFilterONT=s" => \$MFopt{bamfilterONT},
		"mapSaveCRAM=i" => \$MFopt{mapSaveCram},
		#"redoMapping=i" =>\$MFopt{redoMapping},
	#mapping related (2) (assembly)
		"remap2assembly|redoMap2assembly|redoMapping=i" => \$MFopt{redoAssMapping},
		"JGIdepths=i" => \$MFopt{DoJGIcoverage},
		"mapReadsOntoAssembly=i" => \$MFopt{map2Assembly} ,  #map original reads back on assembly, to estimate abundance etc
		"mapSupportReadsOntoAssembly=i" => \$MFopt{mapSupport2Assembly}, # (1) map "SupportReads" onto assembly. Default: 0
		"saveReadsNotMap2Assembly=i" => \$MFopt{SaveUnalignedReads},
	#map2tar / map2DB / map2GC
		"decoyMapping=i" => \$MFopt{DoMapModeDecoy},	#1: "Decoy mapping": map against reference genome AND against assembly of metagenome (drawing obvious better hits to metagenome, the "decoy")
		"competitive2ndmap=i" => \$MFopt{mapModeTogether}, #1: Competitive, 2: combined but report separately per input genome, -1: combined and report all together 
		"ref=s" => \$MFopt{refDBall},
		"mapperLargeRef=i" => \$MFopt{largeMapperDB}, #use flags in mapper index built for large ref DBs?
		"mapnms=s" => \$MFopt{bwt2NameAll}, #name for this final files
		"redo2ndmap=i" => \$MFopt{MapRewrite2nd},
	#SNPs 
		"get2ndMappingConsSNP=i" => \$MFopt{Do2ndMapSNP},#SNPs (onto mapping)
		"getAssemblConsSNP=i" => \$MFopt{DoConsSNP},  #SNPs (onto self assembly) #calculates consensus SNP of assembly (useful for checking assembly gets consensus and Assmbl_grps)
		"getAssemblConsSNPsuppRds=i" => \$MFopt{DoSuppConsSNP}, #same as getAssemblConsSNP, but SNP calling for support reads
		"redoAssmblConsSNP=i" => \$MFopt{redoSNPcons},
		"SNPmem=i" => \$MFopt{memPJob}, #memory per assigned core, in GB
		"redoGeneExtrSNP=i" => \$MFopt{redoSNPgene},
		"SNPjobSsplit=i" => \$MFopt{SNPconsJobsPsmpl}, #parallel jobs per sample; 0 estimates from alignment size
		"SNPminCallQual=i" => \$MFopt{SNPminCallQual},
		"SNPsaveVCF=i" => \$MFopt{saveVCF}, #save vcf of SNP calles? DEfault : 1
		"SNPsaveConsFasta=i" => \$MFopt{saveConsFastas}, #Save consensus fasta from vcf calls? Default: 0 -> too large, can be quickly recreated..
		"SNPcaller=s" => \$MFopt{SNPcallerFlag},
		"SNPcores=i" => \$MFopt{maxSNPcores},
		#"SNPmem=i" => \$MFopt{memSNPcall}, #memory for consensus SNP job in Gb
		"SNPconsMinDepth=i" => \$MFopt{consSNPminDepth}, #how many reads coverage to include position for consensus call?
		"SNPnormINDEL=i" => \$MFopt{normSNPindels}, #using bcftools norm to left-align indels
		"SVcaller=i"  => \$MFopt{callSVs}, #calling structural variants: 1=delly, 2=gridss. Default (0).

	#functional profiling (diamond)
		"profileFunct=i"=> \$MFopt{DoDiamond},
		"reParseFunct=i" => \$MFopt{redoDiamondParse},
		"reProfileFunct=i" => \$MFopt{rewriteDiamond},
		"reProfileFuncTogether=i" => \$MFopt{rewriteAllIfAnyDiamond}, #if any func database needs to be redone, than redo all indicated databases (useful if number of reads used changes..)
		"DiaCores=i" => \$MFopt{diaCores},
		"DiaMem=i" => \$MFopt{diamondMem}, # memory in GB for diamond alignment jobs
		"DiaParseEvals=s" => \$MFopt{diaEVal}, #evalues at which to accept hits to func database
		"DiaSensitiveMode=i" => \$MFopt{diaRunSensitive},
		"DiaFrameshift=i" => \$MFopt{diaFrameshift}, #diamond -F/--frameshift penalty (e.g. 15); enables frameshift-aware alignment for long, error-prone reads. 0 = off (default)
		"rmRawDiamondHits=i" => \$MFopt{DiaRmRawHits},
		"DiaMinAlignLen=i" => \$MFopt{DiaMinAlignLen},
		"DiaMinFracQueryCov=f" =>  \$MFopt{DiaMinFracQueryCov},
		"DiaPercID=i" => \$MFopt{DiaPercID},
		"DiaDBs=s" => \$MFopt{reqDiaDB},#NOG,MOH,ABR,ABRc,ACL,KGM,CZy,PTV,PAB,MOH2,URE,URacc,AMI
	#functional profiling (Jaime tree)
		"orthoExtract=i" => \$MFopt{calcOrthoPlacement},
	#ribo profiling (miTag)
		"profileRibosome=i" => \$MFopt{DoRibofind},
		"riobsomalAssembly=i"  => \$MFopt{doRiboAssembl},
		"reProfileRibosome=i" => \$MFopt{RedoRiboFind} ,  
		"reRibosomeLCA=i"=> \$MFopt{RedoRiboAssign},
		"riboMaxRds=i" => \$MFopt{riboLCAmaxRds},
		"saveRiboRds=i" => \$MFopt{riboStoreRds},
		"thoroughCheckRiboFinish=i" => \$MFopt{checkRiboNonEmpty},
	#other tax profilers..
		"profileMetaphlan=i"=> \$MFopt{DoMetaPhlan},
		#"profileMetaphlan3=i"=> \$MFopt{DoMetaPhlan3},
		"profileMOTU2=i" => \$MFopt{DoMOTU2},
		"profileKraken=i"=> \$MFopt{DoKraken},
		"profileTaxaTarget=i" => \$MFopt{DoTaxaTarget},
		"estGenoSize=i" => \$MFopt{DoGenoSizeEst}, #estimate average size of genomes in data
		"krakenDB=s"=> \$MFopt{globalKraTaxkDB}, #"virusDB";#= "minikraken_2015/";
	#D2s distance
		"calcInterMGdistance=i" => \$MFopt{DoCalcD2s},
	#IO for specific uses
		"newFileStructure=s" => \$MFconfig{mocatLinkDir},#just relink raw files for use in mocat
		"upload2EBI=s" => \$MFconfig{uploadRawRds}, #copy human read removed raw files to this dir, named after sample
	#institute specific: EI
		"wcKeyJobs=s" => \$MFconfig{wcKeysForJob},
	#DEBUG
		"OKtoRWassGrps=i" => \$MFconfig{OKtoRWassGrps}, # can delete assemblies, if suspects error in them
	);
	
	
	# ------------------------------------------ options post processing ------------------------------------------
	setConfigFile($MFconfig{configFile});

	die "ERROR:: No mapping file provided (-map)\n" if ($MFconfig{mapFile} eq "");
	if (!$MFopt{DoAssembly}){
		$MFopt{mapSupport2Assembly}=0;$MFopt{map2Assembly}=0;
	}
	die "ERROR:: \"-mappingMem\" argument contains characters: $MFopt{MapperMemory}" if ($MFopt{MapperMemory} !~ m/[\d-]+/);
	die "ERROR:: \"-mapSortMem\" argument contains characters: $MFopt{mapSortMemGb}" if ($MFopt{mapSortMemGb} !~ m/[\d-]+/);
	die "ERROR:: \"-assemblMemory\" argument contains characters: $MFopt{AssemblyMemory}" if ($MFopt{AssemblyMemory} !~ m/[\d-]+/);
	die "ERROR:: -Binner must be one of 0..5\n"
		unless $MFopt{DoMetaBat2} >= 0 && $MFopt{DoMetaBat2} <= 5;
	die "ERROR:: binning requires -assembleMG to select an assembly mode\n"
		if !$MFopt{DoAssembly} && $MFopt{DoMetaBat2};
	die "ERROR:: -BinnerCores must be a positive integer\n"
		unless $MFopt{BinnerCores} > 0;
	die "ERROR:: -BinnerMem must be zero or a positive integer\n"
		unless $MFopt{BinnerMem} >= 0;
	die "ERROR:: binning requires -useCheckM1 1 and/or -useCheckM2 1\n"
		if $MFopt{DoMetaBat2} && !$MFopt{useCheckM1} && !$MFopt{useCheckM2};
	die "ERROR:: -SB_env is only valid with -Binner 2 (SemiBin)\n"
		if $MFopt{SB_env} ne '' && $MFopt{DoMetaBat2} != 2;
	#die "ERROR:: \"-SNPmem\" argument contains characters: $MFopt{memSNPcall}" if ($MFopt{memSNPcall} !~ m/[\d-]+/);
	die "ERROR:: \"-diamondMem\" argument contains characters: $MFopt{diamondMem}" if ($MFopt{diamondMem} !~ m/[\d-]+/);
	$MFopt{diamondMem} = 7 if ($MFopt{diamondMem} <= 0);
	if ($MFopt{MapperMemory} == -1 ){
		if($MFopt{MapperProg} >2){$MFopt{MapperMemory} = 35 ;
		} else {$MFopt{MapperMemory} = 20 ;}	
	}

	if ($MFopt{DoDiamond} && $MFopt{reqDiaDB} eq ""){die "ERROR:: Functional profiling was requested (-profileFunct 1), but no DB to map against was defined (-diamondDBs)\n";}
	$MFopt{AssemblyKmers} = "-k $MFopt{AssemblyKmers}" unless ($MFopt{AssemblyKmers} =~ m/^-k/);
	$MFopt{sdm_opt}->{minSeqLength}=$MFopt{tmpSdmminSL} if ($MFopt{tmpSdmminSL} > 0);
	$MFopt{sdm_opt}->{maxSeqLength}=$MFopt{tmpSdmmaxSL} if ($MFopt{tmpSdmmaxSL} > 0);
	$MFconfig{remove_reads_tmpDir} = 1 if ($MFconfig{DoFreeGlbTmp} || $MFconfig{rmScratchTmp});
	if ($MFopt{bamSortCores} == -1){$MFopt{bamSortCores} = $MFopt{MapperCores};}
	$MFconfig{filterFromSource}=1 if ($MFconfig{unpackZip} );
	@filterHostDB = split /,/,$MFopt{filterHostDB1};
	
	die "ERROR:: SNPcaller argument invalid, has to be \"MPI\" or \"FB\"\n" if ($MFopt{SNPcallerFlag} ne "MPI" && $MFopt{SNPcallerFlag} ne "FB");
	if ($MFopt{DoConsSNP} && !$MFopt{saveConsFastas} && !$MFopt{saveVCF}){die "ERROR:: Can't use -SNPsaveVCF 0 and -SNPsaveConsFasta 0 -> SNP calling would not be saved in any way..\n";}
	die "ERROR:: assembly consensus SNP calling requires -assembleMG to select an assembly mode\n"
		if !$MFopt{DoAssembly} && ($MFopt{DoConsSNP} || $MFopt{DoSuppConsSNP});

	
	#structural variants..
	if ($MFopt{callSVs} == 0){ ;
	}elsif ($MFopt{callSVs} == 1){	$MFopt{SVcallerFlag} = "DL"; #delly
	}elsif ($MFopt{callSVs} == 2){	$MFopt{SVcallerFlag} = "GY"; # gridss
	}else {die"ERROR:: Invalid callSVs option: $MFopt{callSVs}\n";}
	die "ERROR:: structural-variant calling requires -assembleMG to select an assembly mode\n"
		if !$MFopt{DoAssembly} && $MFopt{callSVs};
	
	if ($MFopt{SB_env} ne ""){
		if ($MFopt{SB_env} ne "human_gut" && $MFopt{SB_env} ne "dog_gut" && $MFopt{SB_env} ne "ocean" && $MFopt{SB_env} ne "soil" && 
		$MFopt{SB_env} ne "cat_gut" && $MFopt{SB_env} ne "human_oral" && $MFopt{SB_env} ne "mouse_gut" && $MFopt{SB_env} ne "pig_gut" && 
		$MFopt{SB_env} ne "built_environment" && $MFopt{SB_env} ne "wastewater" && $MFopt{SB_env} ne "chicken_caecum" && $MFopt{SB_env} ne "global"){
			die "-SB_env must be either: human_gut/dog_gut/ocean/soil/cat_gut/human_oral/mouse_gut/pig_gut/built_environment/wastewater/chicken_caecum/global\n";
		}
	}
	if ($MFconfig{defaultReadLength} < 10){print "Warning: extremely short read length set(-inputReadLength): $MFconfig{defaultReadLength}\n";
	}
	if ($MFconfig{defaultReadLengthX} < 10){print "Warning: extremely short read length set(-inputReadLengthSuppl): $MFconfig{defaultReadLengthX}\n";
	}

	if ($MFopt{DoMetaBat2} == 4) {
		if (!$MFopt{DoAssembly}) {
			die "ERROR:: GenomeFace binner requires an assembly!\n";
		}
		print "GenomeFace binner selected:\n";
		print " - FetchMG (F) will be automatically included in contig stats and run first\n";
		print " - GPU queue must be used\n";
	}
	if ($MFopt{DoMetaBat2} == 5) {
		if (!$MFopt{DoAssembly}) {
			die "ERROR:: SCGBinner requires an assembly!\n";
		}
		print "SCGBinner selected: runs on CPU medium queue (~36h, ~32G)\n";
	}

	$MFopt{memPJob} = int($MFopt{memPJob});
	
	#check HDDspace format
	foreach my $k (keys (%HDDspace)){
		$HDDspace{$k} .= "G" unless ($HDDspace{$k} =~ m/G$/);
	}
	
	if ($MFopt{DoAssembly} == 5 && !$MFopt{mapSaveCram}){
		print "deactivating \"-mapSaveCram\", not supported for hybrid assemblies\n";
		$MFopt{mapSaveCram} = 1;
	}
	
	$MFopt{callSVsSupp} = $MFopt{callSVs}; #for now it just enforces doing both suppl and main SVs, if SVs requested at all..
	
	#set up further dependencies for MF
	my $loopSpec = parse_loop_spec($runOptions{loopCount});
	$runOptions{loopCount} = $loopSpec->{loop_count};
	$runOptions{loopInitialCount} = $loopSpec->{loop_count};
	$runOptions{loopWindowSize} = $loopSpec->{window_size};
	print "Loop2completion=$runOptions{loopCount}; Window size=$runOptions{loopWindowSize}\n"
		if ($runOptions{loopCount});
	
	print "Done. ";

}
