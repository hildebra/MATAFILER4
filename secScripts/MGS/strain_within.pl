#!/usr/bin/perl
#This script collects representative (consensus SNP called) DNA sequences from different metagenomic samples for each MGS, and submits buildTree5.pl for each to reconstruct a phylogeny
# check performance: /ei/projects/8/88e80936-2a5d-4f4a-afab-6f74b374c765/data/geneCats/famDrama7/Bin_SB/intra4_28Feb_01D2SV/MGS.10

use warnings;
use strict;

use Getopt::Long qw( GetOptions );
use File::Path qw(make_path remove_tree);
use File::Glob qw(bsd_glob);
use File::Basename qw(basename dirname);
use File::Spec;
use File::Temp qw(tempdir);
use Cwd qw(abs_path getcwd);
use Digest::SHA qw(sha256_hex);



use Mods::GenoMetaAss qw(gzipopen fileGZe fileGZs resolveExistingFile readClstrRev systemW median mean readMapS readFasta getAssemblPath getAssemblGFF getAssemblContigs);
use Mods::Subm qw(qsubSystem emptyQsubOpt qsubSystemJobAlive qsubSystemWaitMaxJobs
	deferredSubmissionDependency);
use Mods::IO_Tamoc_progs qw(getProgPaths truePath);
use Mods::TamocFunc qw(checkMF);
use Mods::geneCat qw(readGene2tax createGene2MGS);
use Mods::math qw(quantileArray);
use Mods::MGSLocus qw(build_locus_groups choose_locus_candidate protein_kmer_similarity robust_depth_mask);
use Mods::MosaicLoci qw(read_mosaic_catalogue);
use Mods::StrainQC qw(breakpoint_gene_mask abundance_pattern_mask);
use Mods::StrainParts qw(
	balance_assembly_groups choose_auto_worker_count exact_worker_parts write_split_generation write_worker_completion
	split_generation_complete clear_split_generation
);
use Mods::SlurmAccounting qw(
	slurm_tree_memory_summary format_slurm_tree_memory_summary
	next_oom_retry_memory_mb
);
use Mods::WorkflowResilience qw(
	retry_operation retry_unlink retry_rename retry_open retry_close
	write_workflow_record
);
use Mods::CatalogPaths qw(catalog_identity resolve_catalog_maps);
use Mods::StrainSampleStats qw(
	sample_stat_columns sample_summary_columns aggregate_sample_rows
	encode_loci_histogram loci_histogram_rows
);

sub extractFNAFAA2genes;
sub histoMGS;
sub readGenesSample_Singl;
sub reportingsMGS;
sub prepRun;
sub resolveScratchDirectory;
sub persistScratchDirectory;
sub migrateLegacyOperationalLogs;
sub printEarlyRunHeader;
sub prepGene2MGS;
sub createAGlist; sub preComputeConsSNP;
sub mergeConspecificLogs;
sub timeNice;
#sub combineMGSgenes;
sub combineMGSgenesDir; sub prepareMGSInputSet; sub collectMGSShardHandoff;
sub writeMGSShardManifest; sub readSplitGeneration;
sub splitWorkerPartsRemain; sub getInputSize;
sub persistentMGSInputState; sub scratchMGSInputState;
sub invalidateMGSInputState;
sub stagedMGSInputsReady; sub evalFileStatus;
sub preparedOutgroupLog {
	my ($directory) = @_;
	my $log_path = "$directory/data.log";
	return (0, '') unless fileGZe($log_path);
	my ($log) = gzipopen($log_path, "prepared outgroup log", 1);
	my $line = <$log> // '';
	close $log or return (0, '');
	$line =~ s/[\r\n]+\z//;
	return (0, '') unless $line =~ /^OG:(.*)\z/;
	return (1, $1);
}


sub addOutgroup2MGS;
sub writeTooFewMarker;
sub treeOutgroupCandidates;
sub writeNoRecoverableLociMarker;
sub recordValidatedEmptyExtractions;
sub validateTreeInputResolution;
sub epaOnlyRetryReady;
sub prepareEpaOnlyRetryState;
sub resetMGSTreeOutputs;
sub stepComplete;
sub stepProgress;
sub preparedMainBranchInputSet;
sub outgroupRequirementLoci;
sub writeSortedIdentifierSet;
sub prepareSelectiveOutgroupReferenceCache;
sub dispatchPendingTreeJobs;
sub retryOOMTreeJobs;
sub usage;
sub writeRecoveryRow;
sub mergeRecoveryLogs;
sub indexRecoveryRow;
sub writeRecoveryContributionIndex;
sub loadRecoveryContributionIndex;
sub writeStrainSummary;
sub writeSelectionAttritionSummary;
sub writeMGSSampleHistograms;
sub mergeSampleStats;
sub reportSavedSampleStats;
sub printSampleStatsSummary;
sub recoverCompletedSplitPhaseI;
sub taxonAwareLocusBudgets;
sub phase1WorkersNeedingRetry;
sub phase1WorkerCommand;
sub writePhase1RepairQueue;
sub validatePhase1WorkerLedger;
sub writeStrainWorkflowState;
sub cleanupLegacyStrainWorkflowStateFiles;
sub writeStrainWorkflowHeartbeat;
sub writeStrainWorkflowFailure;
sub writeTreeFailureAudit;
sub lifecycleMarkerReason;

sub limitedWarn;sub limitedNotice;


my %limitedWarningStats;
my %limitedNoticeStats;
my $warningExampleLimit = 5;
my ($workflowStage, $workflowStatePath, $workflowStatus, $workflowReason) =
	('startup', '', 'running', '');
my ($legacyWorkflowHeartbeatPath, $legacyWorkflowFailurePath) = ('', '');

my $completionMessage = "";


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
#.45: distinguish split-worker sparsity from missing catalogue data
#.46: restore every sample in shared assembly groups
#.47: use bounded-memory IQ-TREE 3 pathogen mode with an exact legacy-tree switch
#.48: add tree-only reset and resubmission from published per-MGS inputs
#.49: make IQ-TREE pathogen/CMAPLE mode explicitly opt-in
#.50: account for every tree-submission decision and make termination status explicit
#.51: summarize completion of major initialization and workflow steps
#.52: let -recalcTrees recover complete staged inputs before resubmitting trees
#.53: report per-tree SLURM MaxRSS, OOM events, and requested-memory headroom
#.54: locus-sort first-generation per-MGS FNA/FAA files before compression
#.55: separate the gene cap from QC and add mosaic, breakpoint, abundance, and placement QC
#.56: prepare a confirmed mosaic catalogue on demand for legacy direct invocations
#.57: use persistent catalog identity, map manifests, and explicit abundance paths
#.58: store automatic mosaic catalogues and logs in the binner-local mosaic directory
#.59: bulk-align only genes with mosaic or outgroup comparison potential
#.60: validate and report unique MGS-outgroup connections and their gene links
#.61: submit missing mosaic preprocessing as a prerequisite job and wait for it
#.62: discover mosaics across the raw MGS gene set and merge confirmed chains transitively
#.63: keep only rerun and audit outputs while cleaning Mosaic intermediates
#.64: persist run-wide recovered/filtered MAG, Mosaic, and outgroup statistics
#.65: reuse an existing confirmed Mosaic catalogue without requiring raw inputs
#.66: let tree recalculation extract missing inputs and use job-local tree scratch
#.67: parameter-driven gene caps, lower minimum, and per-sample TSV statistics
#.68: remove blank per-sample progress separators from captured output
#.69: delegate tree input publication and completion checkpoints to buildTree5
#.70: require and cardinality-check every recovered split-worker contribution
#.71: persist merge provenance, isolate per-sample TSV output, and print startup immediately
#.72: initialize sample-statistics columns before the executable workflow
#.73: persist, merge, and summarize sample statistics across extraction workers
#.74: print key:value summaries, use assembly-group terminology, and aggregate retained-locus histograms
#.75: balance Step 1 workers by assembly-group count and samples per group
#.76: relax strain-tree defaults and expose buildTree5 taxon-aware locus selection
#.77: enable taxon-aware selection and scale its hierarchy to the strain gene budget
#.78: use deterministic rate/GC partition merging for strain trees by default
#.79: target deterministic rate/GC partitions by effective called sites
#.80: restore legacy IQ-TREE strain inference by default and gate sparse placements
#.81: enable EPA-ng strict-backbone placement by default and expose its controls
#.82: make split Stage-I scope explicit and keep extraction workers tree-option free
#.83: choose split-worker count automatically from assembly-group and sample load
#.84: balance indivisible assembly groups by their sample-level Phase-I workload
#.85: queue prepared Phase-II tree jobs while scheduler capacity is full
#.86: skip extraction-only consensus audits on Phase-I resume and report saved worker statistics
#.87: finish a completed split Phase-I ledger after a main-worker restart
#.88: persist validated no-tree outcomes and reject unclassified tree inputs
#.89: retry Phase-I workers, quarantine terminal MGS outcomes, and harden filesystem publication
#.90: cache catalogue-wide input states and avoid duplicate full-ledger validation scans
#.91: persist the exact shared scratch directory for reliable cross-run resume
#.96: use the authoritative Phase-I input audit for legacy ledger-free resumes
#.97: bound EPA-ng placement memory and worker threads independently of tree inference
#1.06: make buildTree5 finalize staged category/QC/outgroup overlays and input sorting
#1.07: hand validated worker shards to buildTree5 and compact per-MGS submission output
#1.08: skip catalogue and Mosaic initialization when every final tree input is published
#1.09: cache only exact Mosaic and same-COG outgroup references through successful tree validation
#1.01: republish existing EPA placements through final outlier filtering only
#1.02: invalidate EPA-derived completion state before ordinary saved-command resume
#1.03: accept bare and explicit numeric redo-EPA flags
#1.04: propagate forced EPA redo into existing saved tree commands
#1.05: keep EPA redo in the normal controller path through downstream analysis
#1.10: consolidate controller heartbeat and failure records into one state file
#1.11: begin strain postprocessing from completed trees while retaining quarantined tree outcomes
#1.12: prioritize durable completed-tree evidence during tree-only resume audits
#1.13: pass the source MGS tree to postprocessing for outgroup recovery
#1.14: require broadly prevalent loci for taxon-aware rescue selection
#1.15: retain per-locus nucleotide MSAs whenever population genetics is enabled
#1.16: prefer universal-core guide loci and consolidate final taxon-aware diagnostics
#1.17: require stronger multi-locus/backbone overlap before sparse-sample placement
my $version = 1.17;


my $cmdCall = join(" ", $0, @ARGV) . "\n";


#input args..
my $GCd = "";#$ARGV[0];
my $MGSfile = "";#$ARGV[1];
my $clusterID = 95;
my $numCores = 4;#$ARGV[2];
my $subJob=0;#if 0, is main submitting job..
my $maxSubJob = -1;#-1 auto; 0 disables splitting; positive values are explicit worker counts
my $outDpre = "";
my $locTmpDir = ""; my $locTmpDir1 = "";
my $maxCores = -1;
my $onlySubmit =0;#extract genes anew?
my $reSubmit=0;
my $recalcTrees=0; #remove tree-stage outputs and rebuild from published or complete staged per-MGS inputs
my $treeFile = "";
my $doSubmit=0;
my $subMode="";
my %FILTER_DEFAULT = (
	multi_gene_sample_max => 0.25,
	conspecific_gene_sample_max => 0.05,
	minimum_gene_depth => 1,
	minimum_bad_loci_for_sample_skip => 3,
	minimum_mgs_genes_per_sample => 8,
	maximum_genes_per_sample => 600,
	maximum_tree_loci => 400,
	breakpoint_gene_flank => 50,
	abundance_minimum_loci => 8,
	abundance_minimum_fold => 1 / 3,
	abundance_maximum_fold => 3,
	abundance_maximum_modified_z => 3.5,
	prepare_mosaic_loci => 1,
);
my $multiGeneSmplMax = $FILTER_DEFAULT{multi_gene_sample_max};
my $conspGeneSmplMax = $FILTER_DEFAULT{conspecific_gene_sample_max};
my $minDepthGene  = $FILTER_DEFAULT{minimum_gene_depth};
my $minBadLociForSampleSkip = $FILTER_DEFAULT{minimum_bad_loci_for_sample_skip};
my $noIndels = 1;



my $repairCAT=0;

my $maxNGenes = $FILTER_DEFAULT{maximum_genes_per_sample};
my $treeLocusBudget = $FILTER_DEFAULT{maximum_tree_loci};
my $noGeneLimit = 0;
my $disableQC = 0;
my ($mosaicLociFile, $mosaicMGSFile) = ("", "");
my $MGSabundanceOverride = "";
my $prepareMosaicLoci = $FILTER_DEFAULT{prepare_mosaic_loci};
my $breakpointGeneFlank = $FILTER_DEFAULT{breakpoint_gene_flank};
my $abundanceMinimumLoci = $FILTER_DEFAULT{abundance_minimum_loci};
my $abundanceMinimumFold = $FILTER_DEFAULT{abundance_minimum_fold};
my $abundanceMaximumFold = $FILTER_DEFAULT{abundance_maximum_fold};
my $abundanceMaximumModifiedZ = $FILTER_DEFAULT{abundance_maximum_modified_z};
my @subsetMGS=(); my $subsMGSstr="";
my $MSAprog = 2; ##(0) MSAprobs, (1) clustalO, (2) mafft, (4) MUSCLE5
my $phyloProg = 1; #1=IQ-TREE, 2=VeryFastTree, 3=FastTree
my $iqPathogen = 0; #opt in to IQ-TREE 3 pathogen/CMAPLE mode
my $GenesPerSpecies = 0.2;
my $GeneLengthMin = 0.3;
my $relativeNTFraction = 0.1;
my $NTfiltCount = 0;
my ($placementGenesPerSpecies, $placementRelativeNTFraction, $placementNTfiltCount);
$placementGenesPerSpecies = 0.04; $placementRelativeNTFraction = 0.03;
my $taxonAwareLocusSelection = 1;
my $taxonAwareRescueMinPrevalence = 0.8;
my $preferredCoreGenes = "";
my $compactTaxonAwareDiagnostics = 1;
my $rateMergePartitions = 1;
my $rateMergeMaxBins = 8;
my $rateMergeTargetSites = 30_000;
my $rateMergeMinLoci = 20;
my $rateMergeMinSites = 20_000;
my $strictBackbone = 1;
my $strictBackboneFraction = 0.35;
my $strictBackboneMinSamples = 3;
my $placementMinOverlap = 10_000;
my $epaThreads = 2;
my $epaMaxMemMB = -1; # derive from the per-tree IQ-TREE allowance in buildTree5
my $epaPendantOutlierFactor = 5;
my $epaPendantMinThreshold = 0.02;
my $redoEPAfilter = 0;
my $presortGenes = 1200;
my $useGTDBmg = "GTDB";
my $selfMemGb = 10;
my $mosaicMemGb = 150;
my $phase1WorkerRetries = 2;
my $treeOOMMaxMemGB = 1500;
my $treeOOMRetryRounds = 3;
my $redoSubmissionData = 0;
my $deepRepair = 0;
my $rmMSA = 1; #remove per-locus MSAs unless a downstream analysis requires them
my $doPopGenStats = 1;
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

my $conspecificSpThr = 0.1; #higher fraction of genes being two copies in the same sample (abundance >0), and the whole MGS is removed from that sample
my $MGStoolowGsThr = $FILTER_DEFAULT{minimum_mgs_genes_per_sample};
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
my $QCstdof = "sampleQC.tsv";
my $abundF="/assemblies/metag/ContigStats/Coverage.pergene.gz";
my $recoveryLogName = "strainRecovery.tsv";
my $summaryLogName = "strain_within.summary.log";
my $recoveryLogFH;
my $sampleStatsLogName = "strainSampleStats.tsv";
my $sampleStatsSummaryLogName = "strainSampleStats.summary.tsv";
my $sampleStatsPartFH;
my $bamDepthFsuffix = "-smd.bam.coverage.gz";
my $bamDepthFsuffixSup = ".sup-smd.bam.coverage.gz";
my $mapF2 = "";
my $memMulti = 1; #for buildTree script
my $help = 0;


my @sampleStatColumns = sample_stat_columns();

#$treeFile = $ARGV[3] if (@ARGV > 3);$onlySubmit = $ARGV[4] if (@ARGV > 4);
#$doSubmit = $ARGV[5] if (@ARGV > 5);$subMode = $ARGV[6] if (@ARGV > 6);


GetOptions(
	"GCd=s"          => \$GCd,
	"clusterID=i"    => \$clusterID,
	"outD=s"         => \$outDpre,
	"MGS=s"          => \$MGSfile,
	"map2=s"         => \$mapF2, #to be given to strain2 script
	"nodeTmp|tmpD=s" => \$locTmpDir1, 
	"submit=i"       => \$doSubmit,
	"selfMemGb=i"    => \$selfMemGb,
	"mosaicMemGb=i"  => \$mosaicMemGb,
	"phase1WorkerRetries=i" => \$phase1WorkerRetries,
	"treeOOMMaxMemGB=f" => \$treeOOMMaxMemGB,
	"onlySubmit=i"   => \$onlySubmit, #submit only jobs, or also recreate input fna/faa files? (can take days)
	"reSubmit=i"     => \$reSubmit, #for all MGS: resubmit tree phylo building
	"recalcTrees=i"  => \$recalcTrees, #delete tree outputs and rebuild from existing per-MGS inputs
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
	"maxGenes=i"     => \$maxNGenes, #maximum validated genes retained for each MGS/sample
	"treeLocusBudget=i" => \$treeLocusBudget, #bounded final loci passed to BuildTree selection
	"noGeneLimit=i"  => \$noGeneLimit, #remove only the gene-count cap; QC remains enabled
	"disableQC=i"    => \$disableQC, #expert/debug option: disable biological QC independently of the gene cap
	"mosaicLoci=s"   => \$mosaicLociFile, #catalogue-wide confirmed mosaic/outgroup table
	"mosaicMGS=s"    => \$mosaicMGSFile, #raw SB.clusters used for comprehensive Mosaic discovery
	"MGSabundance=s" => \$MGSabundanceOverride, #explicit MGS abundance matrix for nonstandard guide locations
	"prepareMosaicLoci=i" => \$prepareMosaicLoci, #create the default catalogue if absent
	"flushEvery=i"   => \$appendWriteTrigger, #samples buffered before per-MGS records are flushed
	
	"forceSNPcalls=i"  => \$forceVCF2FNA,
	"preCompConsSNP=i"   => \$preCompCons,
	"MGSsubset=s"    => \$subsMGSstr,
	"submissionMode=s"      => \$subMode,
	"MGset=s"        => \$useGTDBmg,
	
	#used genes fine tuning..
	"MGSminGenesPSmpl=i" => \$MGStoolowGsThr, #less genes than this in a single sample -> rm MGS from sample for strains. default 8
	"multiGeneSmplMax=f" => \$multiGeneSmplMax, #default 0.15
	"conspGeneSmplMax=f" => \$conspGeneSmplMax, #default 0.05
	"minBadLociPSmpl=i" => \$minBadLociForSampleSkip,
	"breakpointGeneFlank=i" => \$breakpointGeneFlank,
	"abundanceMinLoci=i" => \$abundanceMinimumLoci,
	"abundanceMinFold=f" => \$abundanceMinimumFold,
	"abundanceMaxFold=f" => \$abundanceMaximumFold,
	"abundanceMaxModifiedZ=f" => \$abundanceMaximumModifiedZ,
	
	#transferred to buildTRee script..
	"GenesPerSpecies=f" => \$GenesPerSpecies,
	"GeneLengthMin=f" => \$GeneLengthMin,
	"relativeNTFraction=f" => \$relativeNTFraction,
	"NTfiltCount=i" => \$NTfiltCount,
	"placementGenesPerSpecies=f" => \$placementGenesPerSpecies,
	"placementRelativeNTFraction=f" => \$placementRelativeNTFraction,
	"placementNTfiltCount=i" => \$placementNTfiltCount,
	"preferredCoreGenes=s" => \$preferredCoreGenes,
	"compactTaxonAwareDiagnostics=i" => \$compactTaxonAwareDiagnostics,
	"taxonAwareLocusSelection=i" => \$taxonAwareLocusSelection,
	"taxonAwareRescueMinPrevalence=f" => \$taxonAwareRescueMinPrevalence,
	"rateMergePartitions=i" => \$rateMergePartitions,
	"rateMergeMaxBins=i" => \$rateMergeMaxBins,
	"rateMergeTargetSites=i" => \$rateMergeTargetSites,
	"rateMergeMinLoci=i" => \$rateMergeMinLoci,
	"rateMergeMinSites=i" => \$rateMergeMinSites,
	"strictBackbone=i" => \$strictBackbone,
	"strictBackboneFraction=f" => \$strictBackboneFraction,
	"strictBackboneMinSamples=i" => \$strictBackboneMinSamples,
	"placementMinOverlap=i" => \$placementMinOverlap,
	"epaThreads=i" => \$epaThreads,
	"epaMaxMemMB=i" => \$epaMaxMemMB,
	"epaPendantOutlierFactor=f" => \$epaPendantOutlierFactor,
	"epaPendantMinThreshold=f" => \$epaPendantMinThreshold,
	"redoEPAfilter:i" => sub { $redoEPAfilter = $_[1] || 1; },
	"MSAprog=i"      => \$MSAprog, #2=MAFFT, 4=muscle5
	"phyloProg=i"    => \$phyloProg, #1=IQ-TREE, 2=VeryFastTree, 3=FastTree
	"iqPathogen=i"   => \$iqPathogen, #explicitly enable IQ-TREE 3 pathogen/CMAPLE mode
	"rmMSA=i"        => \$rmMSA, #remove MSA, to save diskspace
	"popGenStats=i"  => \$doPopGenStats, #requires retained per-locus nucleotide MSAs
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
	"help|h" => \$help,

) or die "Invalid strain_within.pl options\n";
if ($help) {
	print usage(\%FILTER_DEFAULT);
	exit 0;
}
die "Unexpected positional arguments: @ARGV\n" if @ARGV;
die "-redoEPAfilter must be 0 or 1\n"
	unless $redoEPAfilter == 0 || $redoEPAfilter == 1;
checkMF();
die "-GCd is required and must be a directory\n" unless length($GCd) && -d $GCd;
die "Either -MGS or -outD is required\n" unless length($MGSfile) || length($outDpre);
die "MGS file missing or empty: $MGSfile\n" if length($MGSfile) && !-s $MGSfile;
die "-MGset must be GTDB or FMG\n" unless $useGTDBmg eq "GTDB" || $useGTDBmg eq "FMG";
die "-clusterID must be between 1 and 100\n" unless $clusterID >= 1 && $clusterID <= 100;
die "Invalid subjob settings: -maxSubJob must be -1 (auto), 0 (disabled), or positive\n"
	if $maxSubJob < -1 || $subJob < 0 || ($maxSubJob > 0 && $subJob >= $maxSubJob);
die "Core, memory, and precompute settings must be non-negative\n"
	if $maxCores == 0 || $selfMemGb <= 0 || $mosaicMemGb <= 0 || $preCompCons < 0;
die "-minBadLociPSmpl must be positive\n" unless $minBadLociForSampleSkip > 0;
die "-MGSminGenesPSmpl and -presortGenes must be positive\n"
	unless $MGStoolowGsThr > 0 && $presortGenes > 0;
die "-treeLocusBudget must be positive\n" unless $treeLocusBudget > 0;
die "-flushEvery must be positive\n" unless $appendWriteTrigger > 0;
die "-phase1WorkerRetries must be between 0 and 10\n"
	unless $phase1WorkerRetries >= 0 && $phase1WorkerRetries <= 10;
die "-treeOOMMaxMemGB must be positive\n" unless $treeOOMMaxMemGB > 0;
die "Fractional filtering options must be between 0 and 1\n"
	if grep { $_ < 0 || $_ > 1 } ($multiGeneSmplMax, $conspGeneSmplMax,
		$GenesPerSpecies, $GeneLengthMin, $relativeNTFraction,
		$taxonAwareRescueMinPrevalence,
		grep { defined } ($placementGenesPerSpecies, $placementRelativeNTFraction));
die "-NTfiltCount and -placementNTfiltCount must be non-negative\n"
	if $NTfiltCount < 0 || (defined($placementNTfiltCount) && $placementNTfiltCount < 0);
die "-compactTaxonAwareDiagnostics must be 0 or 1\n"
	unless $compactTaxonAwareDiagnostics == 0 || $compactTaxonAwareDiagnostics == 1;
die "-taxonAwareLocusSelection must be 0 or 1\n"
	unless $taxonAwareLocusSelection == 0 || $taxonAwareLocusSelection == 1;
die "-rateMergePartitions must be 0 or 1\n"
	unless $rateMergePartitions == 0 || $rateMergePartitions == 1;
die "-rateMergeMaxBins, -rateMergeTargetSites, -rateMergeMinLoci, and -rateMergeMinSites must be positive\n"
	if grep { $_ < 1 } ($rateMergeMaxBins, $rateMergeTargetSites, $rateMergeMinLoci, $rateMergeMinSites);
die "-strictBackbone must be 0 or 1\n"
	unless $strictBackbone == 0 || $strictBackbone == 1;
die "-strictBackboneFraction must be between 0 and 1\n"
	if $strictBackboneFraction < 0 || $strictBackboneFraction > 1;
die "-strictBackboneMinSamples must be at least 3\n"
	if $strictBackboneMinSamples < 3;
die "-placementMinOverlap must be non-negative\n"
	if $placementMinOverlap < 0;
die "-epaThreads must be positive\n" if $epaThreads < 1;
die "-epaMaxMemMB must be -1 (derived), 0 (no memory-based scaling), or positive\n"
	if $epaMaxMemMB < -1;
die "-epaPendantOutlierFactor and -epaPendantMinThreshold must be non-negative\n"
	if $epaPendantOutlierFactor < 0 || $epaPendantMinThreshold < 0;
if ($redoEPAfilter) {
	die "-redoEPAfilter requires -strictBackbone 1 and -phyloProg 1\n"
		unless $strictBackbone && $phyloProg == 1;
	die "-redoEPAfilter cannot be combined with tree/input regeneration modes\n"
		if $recalcTrees || $reSubmit || $repairCAT || $deepRepair || $redoSubmissionData;
	die "-redoEPAfilter is a parent-only resume mode\n" if $subJob;
	$onlySubmit = 1;
}
my ($taxonAwareGeneBudget, $taxonAwareMaxLoci,
	$taxonAwareCoreLoci, $taxonAwareCandidateExtra) = (0, 0, 0, 0);
if ($taxonAwareLocusSelection) {
	$taxonAwareGeneBudget = $treeLocusBudget < $presortGenes
		? $treeLocusBudget : $presortGenes;
	($taxonAwareMaxLoci, $taxonAwareCoreLoci, $taxonAwareCandidateExtra) =
		taxonAwareLocusBudgets($taxonAwareGeneBudget);
}
die "SNP depth, quality, adaptive filtering, and indel-range settings must be non-negative\n"
	if $minSNPDepth < 0 || $minSNPCallQual < 0 || $useAdaptiveQual < 0
		|| $depthFilterScale < 0 || $indelRange < 0;
die "-phyloMemMulti must be positive\n" unless $memMulti > 0;
die "-phyloProg must be 1 (IQ-TREE), 2 (VeryFastTree), or 3 (FastTree)\n"
	unless $phyloProg >= 1 && $phyloProg <= 3;
die "-iqPathogen must be 0 or 1\n" unless $iqPathogen == 0 || $iqPathogen == 1;
die "-iqPathogen requires -phyloProg 1\n" if $iqPathogen && $phyloProg != 1;
die "-recalcTrees must be 0 or 1\n" unless $recalcTrees == 0 || $recalcTrees == 1;
die "-noGeneLimit and -disableQC must each be 0 or 1\n"
	unless ($noGeneLimit == 0 || $noGeneLimit == 1)
		&& ($disableQC == 0 || $disableQC == 1);
die "-prepareMosaicLoci must be 0 or 1 "
	."(default $FILTER_DEFAULT{prepare_mosaic_loci})\n"
	unless $prepareMosaicLoci == 0 || $prepareMosaicLoci == 1;
die "-breakpointGeneFlank must be non-negative\n" unless $breakpointGeneFlank >= 0;
die "Abundance-pattern settings are invalid\n"
	unless $abundanceMinimumLoci > 0 && $abundanceMinimumFold > 0
		&& $abundanceMaximumFold >= $abundanceMinimumFold
		&& $abundanceMaximumModifiedZ >= 0;
die "-recalcTrees cannot be combined with -repairCAT, -deepRepair, or -redoSubmissionData\n"
	if $recalcTrees && ($repairCAT || $deepRepair || $redoSubmissionData);
die "-recalcTrees must be launched by the main strainWithin process, not a split worker\n"
	if $recalcTrees && $subJob;
die "-MSAprog must be 0, 1, 2, or 4\n"
	unless grep { $MSAprog == $_ } (0, 1, 2, 4);
die "-rmMSA must be 0 or 1\n" unless $rmMSA == 0 || $rmMSA == 1;
die "-popGenStats must be 0 or 1\n"
	unless $doPopGenStats == 0 || $doPopGenStats == 1;
if ($doPopGenStats && $rmMSA) {
	warn "Population genetics requires per-locus nucleotide MSAs; overriding -rmMSA 1 to -rmMSA 0\n";
	$rmMSA = 0;
}

$GCd = abs_path($GCd);
$GCd .= "/" unless $GCd =~ m{/$};
$MGSfile = abs_path($MGSfile) if length $MGSfile;
$mosaicMGSFile = File::Spec->rel2abs($mosaicMGSFile) if length $mosaicMGSFile;
if (!length($preferredCoreGenes) && length($MGSfile)) {
	my $companionCoreGuide = $MGSfile =~ /\.core\z/
		? $MGSfile
		: "$MGSfile.core";
	$preferredCoreGenes = $companionCoreGuide if -s $companionCoreGuide;
}
if (length($preferredCoreGenes)) {
	$preferredCoreGenes = abs_path($preferredCoreGenes)
		or die "Cannot resolve -preferredCoreGenes: $preferredCoreGenes\n";
	die "-preferredCoreGenes is missing or empty: $preferredCoreGenes\n"
		unless -s $preferredCoreGenes;
}
$outDpre = File::Spec->rel2abs($outDpre) if length $outDpre;
$mosaicLociFile = File::Spec->rel2abs($mosaicLociFile) if length $mosaicLociFile;
$MGSabundanceOverride = File::Spec->rel2abs($MGSabundanceOverride)
	if length $MGSabundanceOverride;

if (!length($mosaicMGSFile) && length($MGSfile)) {
	$mosaicMGSFile = $MGSfile;
	$mosaicMGSFile =~ s/\.core\z//;
}

$noGeneLimit = 1 if $maxNGenes <= 0; #backward-compatible no-cap spelling; QC is unchanged
die "-maxGenes must be at least -MGSminGenesPSmpl unless -noGeneLimit 1 is used\n"
	if !$noGeneLimit && $maxNGenes < $MGStoolowGsThr;
$maxNGenes = -1 if $noGeneLimit;

$onlySubmit = 1 if $recalcTrees; #tree-only recovery reuses published or complete staged inputs
printEarlyRunHeader();

@subsetMGS = split /,/,$subsMGSstr if ($subsMGSstr ne "");
#print "SUBSMGS:: @subsetMGS\n";
#die timeNice(20) ." ".timeNice(12252)."\n"; #TEST

#define global vars
my $queueMode = $subMode;
$queueMode = "bash" if !$doSubmit && $queueMode eq "";
my $QSBoptHR = emptyQsubOpt($doSubmit,"",$queueMode);
my $MGSfileOri = $MGSfile; #save for later..


my $resumeBindir = $MGSfile;
$resumeBindir =~ s/[^\/]+$//;
$resumeBindir = $GCd if $resumeBindir eq "";
my $resumeOutD = length($outDpre) ? $outDpre : "$resumeBindir/intra_phylo/";

my ($preparedMainBranchFastPath, @preparedMainBranchMGS) = (0);
my %preparedMainBranchCategoryValidated;
if ($onlySubmit && !$subJob && length($MGSfile)
		&& !$redoSubmissionData && !$repairCAT && !$deepRepair) {
	my ($ready, $mgs, $reason) = preparedMainBranchInputSet(
		$MGSfile, $resumeOutD, \@subsetMGS,
	);
	if ($ready) {
		$preparedMainBranchFastPath = 1;
		@preparedMainBranchMGS = @{$mgs};
		print "Prepared-input recovery: ".scalar(@preparedMainBranchMGS)
			." MGS have final FASTA/category/outgroup state; skipping Mosaic "
			."and gene-catalogue initialization while retaining the normal workflow. "
			."Global elapsed ".timeNice(time - $^T)."\n";
	} else {
		print STDERR "Prepared-input recovery unavailable ($reason); using full initialization.\n";
	}
}

# Redoing EPA filtering is a local recovery step: invalidate only the
# published EPA-derived artifacts, then rejoin the shared controller path.
# placementPending.sto remains reserved for the separate EPA-only recovery path.
if ($redoEPAfilter) {
	die "-redoEPAfilter output directory does not exist: $resumeOutD\n"
		unless -d $resumeOutD;
	my %subset = map { $_ => 1 } @subsetMGS;
	my ($retained, $removed, $completionRemoved, $pendingRemoved,
		$alreadyMissing) = (0, 0, 0, 0, 0);
	for my $jplace (bsd_glob(File::Spec->catfile(
			$resumeOutD, '*', 'phylo', 'epa-ng', 'epa_result.jplace'))) {
		next unless -s $jplace;
		my $mgsDir = dirname(dirname(dirname($jplace)));
		my $mgs = basename($mgsDir);
		next if %subset && !$subset{$mgs};
		next if -s File::Spec->catfile($mgsDir, 'noTree.sto');
		$retained++;
		my $placedTree = File::Spec->catfile(
			$mgsDir, 'phylo', 'IQtree_allsites.treefile');
		if (-e $placedTree) {
			if ($doSubmit) {
				retry_unlink($placedTree,
					label => "invalidate EPA-filtered tree for $mgs");
				$removed++;
			} else {
				print "Would remove EPA-filtered tree $placedTree\n";
			}
		} else {
			$alreadyMissing++;
		}
		my $completion = File::Spec->catfile($mgsDir, 'treeDone.sto');
		if (-e $completion) {
			if ($doSubmit) {
				retry_unlink($completion,
					label => "invalidate EPA-filtered completion marker for $mgs");
				$completionRemoved++;
			} else {
				print "Would remove EPA-filtered completion marker $completion\n";
			}
		}
		my $pending = File::Spec->catfile($mgsDir, 'placementPending.sto');
		if (-e $pending) {
			if ($doSubmit) {
				retry_unlink($pending,
					label => "clear stale EPA-only marker for normal resume of $mgs");
				$pendingRemoved++;
			} else {
				print "Would remove stale EPA-only marker $pending\n";
			}
		}
	}
	print "Redo EPA filter resume: retained_jplace=$retained, "
		."placed_trees_removed=$removed, "
		."completion_markers_removed=$completionRemoved, "
		."pending_markers_removed=$pendingRemoved, "
		."already_missing=$alreadyMissing. "
		."Continuing through the normal controller workflow.\n";
}
if (length($MGSfile) && !$preparedMainBranchFastPath) {
	my $explicitMosaicCatalogue = length($mosaicLociFile);
	if (!$explicitMosaicCatalogue && $prepareMosaicLoci) {
		my $mosaicDirectory = File::Spec->catdir(dirname($mosaicMGSFile), 'mosaic');
		make_path($mosaicDirectory);
		$mosaicLociFile = File::Spec->catfile(
			$mosaicDirectory,
			basename($mosaicMGSFile).".mosaic_loci.$clusterID.confirmed.tsv",
		);
	}
	if (length($mosaicLociFile) && -s $mosaicLociFile) {
		print STDERR ($subJob
			? "Loading confirmed Mosaic catalogue for split worker: $mosaicLociFile\n"
			: "Reusing existing confirmed Mosaic catalogue: $mosaicLociFile\n");
	}
	if (length($mosaicLociFile) && !-s $mosaicLociFile) {
		if (!$prepareMosaicLoci) {
			die "Confirmed mosaic catalogue is missing or empty: $mosaicLociFile\n"
				if $explicitMosaicCatalogue;
			$mosaicLociFile = "";
		} elsif ($subJob) {
			die "Confirmed mosaic catalogue is unavailable to split worker $subJob: "
				."$mosaicLociFile. Run the main strain_within.pl process first.\n";
		} else {
			my $mosaicDirectory = dirname($mosaicLociFile);
			die "Raw MGS assignment file for Mosaic is missing or empty: $mosaicMGSFile\n"
				unless -s $mosaicMGSFile;
			make_path($mosaicDirectory);
			my $mosaicRunDirectory = tempdir(
				'prepare-mosaic-XXXXXX', DIR => $mosaicDirectory,
				CLEANUP => 0,
			);
			my $mosaicScript = getProgPaths("MGS_mosaic_scr");
			my $mosaicThreads = $maxCores > 0 ? $maxCores : $numCores;
			$mosaicThreads = 1 if $mosaicThreads < 1;
			my $mosaicCommand = join(" ",
				$mosaicScript,
				"-GCd", shellQuote($GCd),
				"-MGS", shellQuote($mosaicMGSFile),
				"-coreMGS", shellQuote($MGSfile),
				"-clusterID", $clusterID,
				"-threads", $mosaicThreads,
				"-output", shellQuote($mosaicLociFile),
			);
			$mosaicCommand .= " -tmpD ".shellQuote($locTmpDir1)
				if length($locTmpDir1) && -d $locTmpDir1;
			my $mosaicLog = File::Spec->catfile(
				$mosaicRunDirectory, 'prepare_mosaic_loci.log',
			);
			$mosaicCommand .= " > ".shellQuote($mosaicLog)." 2>&1\n";
			$mosaicCommand .= "test -s ".shellQuote($mosaicLociFile)."\n";
			my $mosaicJobScript = File::Spec->catfile(
				$mosaicRunDirectory, 'prepare_mosaic_loci.sh',
			);
			print "Confirmed mosaic catalogue is absent; submitting prerequisite "
				."Mosaic job with $mosaicThreads cores and ${mosaicMemGb}G memory\n";
			my ($mosaicDependency, $mosaicSubmissionCommand) = qsubSystem(
				$mosaicJobScript, $mosaicCommand, $mosaicThreads,
				"${mosaicMemGb}G", "MosaicMGS", "", "", 1, [], $QSBoptHR,
			);
			unless ($doSubmit) {
				$completionMessage = "Mosaic prerequisite script was generated at "
					."$mosaicJobScript; submission was disabled.";
				print "Mosaic prerequisite was not submitted because -submit 0; "
					."stopping before Mosaic-dependent strain extraction.\n";
				exit 0;
			}
			print "Waiting for prerequisite Mosaic job to finish before loading "
				."loci/outgroups or starting extraction workers\n";
			qsubSystemJobAlive([$mosaicDependency], $QSBoptHR);
			die "Mosaic preprocessing completed without a nonempty catalogue: "
				."$mosaicLociFile. Inspect $mosaicLog and "
				."$mosaicJobScript.etxt\n" unless -s $mosaicLociFile;
			print "Prerequisite Mosaic catalogue is ready: $mosaicLociFile\n";
			remove_tree($mosaicRunDirectory);
			print "Removed successful Mosaic job workspace $mosaicRunDirectory\n";
		}
	}
}
cleanupMosaicIntermediates($mosaicLociFile) if !$preparedMainBranchFastPath && length($mosaicLociFile) && -s $mosaicLociFile;

my $bindir;my $outD;my $scratchD;my $preConDir;my $LOGDIR;my $mapF;
my %map; my %AsGrps;my @samples;#map and assembly groups
my %ConspecificMGS; #list of conspecific MGS
my %MGSnoTree; #MGS known to have a persistent valid no-tree outcome
my %MGSnoTreeReason;
my %MGSepaOnlyRetry;
my %MGSsubmissionComplete;
my %deferredScratchCleanup;
my $legacyLocusOutputs = 0;
my %legacyLocusMGS;
my (%ConfirmedMosaicPairs, %PreferredOutgroup, %PreferredOutgroupGene);
if (length($mosaicLociFile) && !$preparedMainBranchFastPath) {
	my $mosaicStarted = time;
	my $nextMosaicProgress = time + 60;
	my ($pairs, $outgroups, $outgroup_genes) = read_mosaic_catalogue(
		$mosaicLociFile,
		sub {
			my ($status) = @_;
			return if time < $nextMosaicProgress;
			stepProgress("Mosaic catalogue loading", $status->{rows_scanned}, undef,
				$mosaicStarted, "file=".basename($status->{file}));
			$nextMosaicProgress = time + 60;
		},
	);
	%ConfirmedMosaicPairs = %{$pairs};
	%PreferredOutgroup = %{$outgroups};
	%PreferredOutgroupGene = %{$outgroup_genes};
	if ($subJob) {
		print "Using Mosaic catalogue: ".scalar(keys %ConfirmedMosaicPairs)
			." confirmed pair(s), ".scalar(keys %PreferredOutgroup)
			." outgroup choice(s)\n";
	} else {
		print "Loaded ".scalar(keys %ConfirmedMosaicPairs)." confirmed mosaic pair(s) and "
			.scalar(keys %PreferredOutgroup)." consolidated outgroup choice(s) from $mosaicLociFile\n";
	}
	my $connection_number = 0;
	my $gene_link_number = 0;
	for my $source (sort keys %PreferredOutgroup) {
		my @queries = sort keys %{$PreferredOutgroupGene{$source} || {}};
		$gene_link_number += scalar(@queries);
		$connection_number++;
		next if $subJob;
		next if $connection_number > 10;
		my @preview_queries = @queries
			? @queries[0 .. ($#queries < 5 ? $#queries : 5)]
			: ();
		my @preview = map { $_.'->'.$PreferredOutgroupGene{$source}{$_} }
			@preview_queries;
		print STDERR "  Mosaic outgroup $source -> $PreferredOutgroup{$source}: "
			.scalar(@queries)." proposed gene link(s)"
			.(@preview ? " (".join(', ', @preview)
				.(scalar(@queries) > @preview ? ", ..." : "").")" : "")."\n";
	}
	if (!$subJob) {
		print "Mosaic outgroup proposals loaded: ".scalar(keys %PreferredOutgroup)
			." unique MGS-to-MGS connection(s), $gene_link_number gene-to-gene link(s)";
		print "; detailed previews limited to 10 connections" if $connection_number > 10;
		print "\n";
	}
} elsif (!$preparedMainBranchFastPath) {
	warn($prepareMosaicLoci
		? "No confirmed mosaic catalogue supplied; same-COG catalogue clusters will remain separate\n"
		: "Mosaic checks disabled; same-COG catalogue clusters will remain separate and tree-based outgroups remain available\n");
}

my $gene2taxF; #where to find info what genes (gene cat)
my $sttime = $^T;
my $stepStarted = time;
prepRun();
$workflowStatePath = File::Spec->catfile($LOGDIR,
	$subJob ? "strain_within.worker.$subJob.state.tsv" : 'strain_within.state.tsv');
$legacyWorkflowHeartbeatPath = File::Spec->catfile($LOGDIR,
	$subJob ? "strain_within.worker.$subJob.heartbeat.tsv" : 'strain_within.heartbeat.tsv');
$legacyWorkflowFailurePath = File::Spec->catfile($LOGDIR,
	$subJob ? "strain_within.worker.$subJob.failure.tsv" : 'strain_within.failure.tsv');
writeStrainWorkflowHeartbeat('configuration');
stepComplete("configuration and map initialization", $stepStarted,
	"samples=".scalar(@samples), "mode=$mode", "output=$outD");


my %AGlist; #list of assembly groups that need to be processed together;
$stepStarted = time;
if ($preparedMainBranchFastPath) {
	$maxSubJob = 0 if $maxSubJob == -1;
	stepComplete("assembly-group expansion", $stepStarted,
		"status=skipped_prepared_inputs", "groups=not_loaded",
		"samples=".scalar(@samples));
} else {
createAGlist();
my $groupedSampleCount = 0;
$groupedSampleCount += scalar(@{$AGlist{$_}}) for keys %AGlist;
if ($maxSubJob == -1) {
	my ($automaticWorkers, $targetGroupsPerWorker) = choose_auto_worker_count(
		scalar(keys %AGlist), scalar(@samples),
	);
	$maxSubJob = $automaticWorkers;
	print "Automatic Stage-I splitting: ".scalar(keys %AGlist)." assembly groups, "
		.scalar(@samples)." samples, target ${targetGroupsPerWorker} groups/worker; "
		.($maxSubJob ? "using $maxSubJob workers" : "using the main process only")."\n";
}
stepComplete("assembly-group expansion", $stepStarted,
	"groups=".scalar(keys %AGlist), "grouped_samples=$groupedSampleCount",
	"standalone_samples=".(scalar(@samples) - $groupedSampleCount));
}
#foreach (sort keys %AGlist) {   print "$_ : @{$AGlist{$_}}\n";}die;

my %preCompSNPs;
my %unavailableSamples;


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
my %MGSneedsExtraction; #selected MGS with neither published nor staged tree inputs
my (%persistentMGSInputStateCache, %scratchMGSInputStateCache);






#key step to determine with set of genes (representing MGS) is to be MSA'd for strain phylos
#these might be very limited number of genes here..
$stepStarted = time;
my @specis;
if ($preparedMainBranchFastPath) {
	($SIgenes, $Gene2COG, $Gene2MGS, $COGprios) = ({}, {}, {}, {});
	@specis = @preparedMainBranchMGS;
	stepComplete("MGS and seed-locus selection", $stepStarted,
		"selected_MGS=".scalar(@specis), "status=skipped_catalogue_index",
		"seed_loci=not_loaded", "catalogue_genes=not_loaded");
} else {
my $seedIndexStarted = time;
my $nextSeedIndexProgress = time + 60;
($SIgenes,$Gene2COG,$Gene2MGS,$COGprios) = readGene2tax(
	$gene2taxF, $presortGenes, \@subsetMGS,
	sub {
		my ($status) = @_;
		return if time < $nextSeedIndexProgress;
		stepProgress("MGS and seed-locus selection", $status->{rows_scanned}, undef,
			$seedIndexStarted, "included_genes=$status->{included_genes}");
		$nextSeedIndexProgress = time + 60;
	},
);#
#%SIgenes=%{$hr1};%Gene2COG=%{$hr2}; %Gene2MGS = %{$hr3}; %COGprios = %{$hr4};
@specis = sort(keys(%{$SIgenes}));
$Gene2MGS = {}; #not consumed by the within-strain workflow
die "No MGS matched the selected input"
	. ($subsMGSstr ne "" ? " or -MGSsubset $subsMGSstr" : "") . "\n"
	unless @specis;
for my $MGS (@specis) {
	die "Unsafe MGS identifier '$MGS': use only letters, digits, dot, underscore, colon, plus, and hyphen\n"
		unless defined($MGS) && $MGS =~ /\A[A-Za-z0-9][A-Za-z0-9_.:+-]*\z/;
}
my $selectedSeedLoci = 0;
$selectedSeedLoci += scalar(keys %{$SIgenes->{$_}}) for @specis;
stepComplete("MGS and seed-locus selection", $stepStarted,
	"selected_MGS=".scalar(@specis), "seed_loci=$selectedSeedLoci",
	"catalogue_genes=".scalar(keys %{$Gene2COG}));
}
#sort specis by numbers, so start with MGS1, MGS2 etc
my %sis; foreach (@specis){if (m/(\d+)$/){ $sis{$_}=int($1);} else {$sis{$_}=1; print "Unknown code: $_";}}
@specis = sort {$sis{$a} <=> $sis{$b} || $a cmp $b } keys %sis;


#die "specis::\n@specis\n";
my $cnt=0; my $SaSe = "|"; 


$stepStarted = time;
my ($dirsNOTPrepped , $CatFileMiss , $CatNotPrepped , $treeAbsent, $doneDirs, $PhylosExist,
	$noRecoverableLociDirs, $completedTreeFastPaths)
			= evalFileStatus();
my $epaOnlyRetryCount = scalar(keys %MGSepaOnlyRetry);
my $legacyEpaRetryCount = scalar(grep {
	($MGSepaOnlyRetry{$_} // '') eq 'legacy_missing_final'
} keys %MGSepaOnlyRetry);
my $fullTreeRetryCount = $treeAbsent - $epaOnlyRetryCount;
$fullTreeRetryCount = 0 if $fullTreeRetryCount < 0;
stepComplete("existing-output and resume audit", $stepStarted,
	"prepared_trees=$doneDirs", "completion_marker_fast_paths=$completedTreeFastPaths",
	"missing_trees=$treeAbsent", "incomplete_tree_inputs=$CatFileMiss",
	"directories_needing_extraction=$dirsNOTPrepped",
	"validated_no_locus=$noRecoverableLociDirs",
	"epa_only_retries=$epaOnlyRetryCount", "legacy_epa_retries=$legacyEpaRetryCount",
	"full_tree_retries=$fullTreeRetryCount");
#DEBUG:getInputSize();


my %smplsPerMGS; #stats: MGS is represented in how many different samples?

#hashes of strings that keep results to be written for each species..
my %OCstrH ; my %OFstrH ; my %OAstrH ; my %OLstrH ; my %OQstrH ;
my $splitGeneration = '';
my $splitManifest = "$LOGDIR/mainExtr.generation";
my $splitStonePrefix = "$LOGDIR/mainExtr";
my (%recoveryWorkersByMGS, %recoveryRecordsByMGS, %recoveryRowsByMGS);
my (%recoveryWorkerRecordsByMGS, %recoveryWorkerRowsByMGS);
my (%recoverySamplesByMGS, %recoveryUniqueSamplesByMGS);
my %stagedShardHandoff;
my $recoveryContributionIndexReady = 0;


my $runPartI = (($dirsNOTPrepped/@specis > 0.1) || $onlySubmit == 0
			|| $subJob || $redoSubmissionData
			|| $legacyLocusOutputs
			|| ($deepRepair && $dirsNOTPrepped)
			|| ($repairCAT && $CatFileMiss)
			|| ($recalcTrees && $dirsNOTPrepped));
if ($runPartI){
	#$PhylosExist=0;
	
	print "\n\n----------------------------------------------------\nPart I:: extracting relevant core MGS genes (SNP consensus called) from original assemblies". "Elapsed time : ", timeNice(time - $sttime) . "\n----------------------------------------------------\n\n";
	# Checking every sample's consensus/VCF inputs is extraction-only work.  In
	# a tree-only resume it used to dominate startup despite no sample data being
	# read afterwards, so keep it strictly within the Phase-I path.
	$stepStarted = time;
	if (!$subJob) {
		preComputeConsSNP();
		stepComplete("consensus-input audit", $stepStarted,
			"usable_samples=".(scalar(@samples) - scalar(keys %unavailableSamples)),
			"unavailable_samples=".scalar(keys %unavailableSamples),
			"precomputed_consensus=".scalar(keys %preCompSNPs));
	} else {
		# The parent validates and, when requested, precomputes every sample before
		# it submits workers. Repeating that full scan in every worker multiplied
		# filesystem metadata traffic without changing extraction decisions.
		print "Split worker $subJob: reusing parent consensus preflight; validating only owned samples during extraction.\n";
		stepComplete("consensus-input audit", $stepStarted,
			"scope=parent_preflight", "worker=$subJob");
	}
	
	$stepStarted = time;
	my @stageIExtractionMGS = $recalcTrees
		? grep { $MGSneedsExtraction{$_} } @specis
		: @specis;
	my $stageIScope = $recalcTrees
		? 'recalcTrees: incomplete input triplets only'
		: length($subsMGSstr) ? 'explicit -MGSsubset'
		: 'all selected MGS';
	print "Stage-I extraction scope: $stageIScope; target_MGS="
		.scalar(@stageIExtractionMGS)."; split_workers=$maxSubJob. "
		."Workers are balanced by assembly group, so a worker with no locus for the target MGS writes no FASTA records.\n"
		if $maxSubJob && !$subJob;
	prepGene2MGS();
	stepComplete("locus-model construction", $stepStarted,
		"catalogue_drivers=".scalar(keys %cl2gene2),
		"resolved_loci=".scalar(keys %{$LocusByID}),
		"selected_MGS=".scalar(keys %{$COGprios}));
	$Gene2COG = {}; #delete, no longer needed..
	
	reportingsMGS();
	%smplsPerMGS = (); #reporting-only sample/locus counts can be large
	$SIgenes = {}; #replaced locus selection is represented by $COGprios
	
	my @jobsMain;
	my $phase1SelfCmd = '';

	if ($maxSubJob && !$subJob){
		# A generation manifest prevents an isolated worker retry from being
		# mistaken for a complete replacement of previously merged inputs.
		clear_split_generation($splitManifest, $splitStonePrefix);
		$splitGeneration = join('.', time, $$, int(rand(1_000_000_000)));
		write_split_generation($splitManifest, $splitGeneration, $maxSubJob);
		#here needs to submit itself maxSubJob times
		$phase1SelfCmd = phase1WorkerCommand();
		
		my $tmpHDD=$QSBoptHR->{tmpSpace} ; $QSBoptHR->{tmpSpace} =15; #request some basic amount
		
		#submission of self-subjobs..
		for (my $sj = 1; $sj < $maxSubJob; $sj ++){
			my $cmdX = "$phase1SelfCmd -subjob $sj &&\n";
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
	my $extractionDriverCount = scalar(keys %cl2gene2);
	my $extractionLocusCount = scalar(keys %{$LocusByID});
	$stepStarted = time;
	extractFNAFAA2genes();#@allGenes);
	%cl2gene2 = (); #no longer needed, delete
	$LocusByID = {};
	$MemberContext = {};
	$LocusContext = {};
	%LocusSeedProteins = ();
	$COGprios = {};
	#write logs to found genes etc.
	writeLogsStep1();
	mergeRecoveryLogs() unless $maxSubJob;
	mergeSampleStats() unless $maxSubJob;
	stepComplete("consensus-gene extraction and publication", $stepStarted,
		"catalogue_drivers=$extractionDriverCount",
		"resolved_loci=$extractionLocusCount",
		"worker=$subJob");
	write_worker_completion("$splitStonePrefix.0.stone", $splitGeneration)
		if $maxSubJob && !$subJob;
	
	if ($subJob){
		$completionMessage = "strain_within.pl subjob ${subJob}/$maxSubJob completed normally.";
		print "Finished subJob ${subJob}/$maxSubJob. Exiting..\n";
		exit(0);
	}

	if ($maxSubJob && !$subJob){ # second part for main worker: check that everything else is finished..
		if (@jobsMain && !$doSubmit) {
			$completionMessage = "split-worker scripts were generated successfully; submission was disabled.";
			print "Split-worker scripts were generated but not submitted; stopping before incomplete outputs are combined.\n";
			exit(0);
		}
		qsubSystemJobAlive( \@jobsMain,$QSBoptHR ) if @jobsMain && $doSubmit;
		my $workerRetryRound = 0;
		while (1) {
			my @failedWorkers = phase1WorkersNeedingRetry($splitGeneration);
			last unless @failedWorkers;
			if ((grep { $_ == 0 } @failedWorkers)
					|| $workerRetryRound >= $phase1WorkerRetries) {
				my $queue = writePhase1RepairQueue($splitGeneration, \@failedWorkers,
					'live Phase-I worker validation failed');
				$completionMessage = "Phase I requires worker repair before Phase II; no tree jobs were submitted.";
				print "Phase-I processing paused safely; repair queue: $queue. Invalid workers: "
					.join(',', @failedWorkers)."\n";
				exit(0);
			}
			$workerRetryRound++;
			print "Retrying Phase-I worker(s) ".join(',', @failedWorkers)
				." (round $workerRetryRound/$phase1WorkerRetries)\n";
			my @retryJobs;
			my $savedTmp = $QSBoptHR->{tmpSpace}; $QSBoptHR->{tmpSpace} = 15;
			for my $sj (@failedWorkers) {
				my $stone = "$splitStonePrefix.$sj.stone";
				retry_unlink($stone, fatal => 0, label => "clear worker $sj completion");
				my $cmdX = "$phase1SelfCmd -subjob $sj &&\n"
					."printf '%s\\n' ".shellQuote($splitGeneration)
					." > ".shellQuote($stone)."\n";
				my ($dependency) = qsubSystem("$LOGDIR/Strain1_B${sj}.retry${workerRetryRound}.sh",
					$cmdX,1,"${selfMemGb}G","Str1.$sj","","",1,[],$QSBoptHR);
				push @retryJobs, $dependency;
			}
			$QSBoptHR->{tmpSpace} = $savedTmp;
			qsubSystemJobAlive(\@retryJobs, $QSBoptHR);
		}
		my $generationComplete = retry_operation(
			label => 'validate completed split-extraction generation', fatal => 0,
			code => sub { split_generation_complete(
				$splitManifest, $splitStonePrefix, $maxSubJob) },
		);
		unless ($generationComplete) {
			my @workers = 0 .. $maxSubJob - 1;
			my $queue = writePhase1RepairQueue($splitGeneration, \@workers,
				'split generation remained incomplete after bounded filesystem retries');
			$completionMessage = "Phase I generation validation is incomplete; no tree jobs were submitted.";
			print "Phase-I processing paused safely; repair queue: $queue\n";
			exit(0);
		}
		mergeConspecificLogs();
		mergeRecoveryLogs();
		mergeSampleStats();
		
		#combineMGSgenes();
	}
	invalidateMGSInputState(@stageIExtractionMGS);
	my $emptyMGS = recordValidatedEmptyExtractions(\@stageIExtractionMGS);
	print "Stage-I terminal-input classification: recorded $emptyMGS MGS with no recoverable loci.\n"
		if $emptyMGS;
	
	print "\nGene extraction & redistribution finished, ready to proceed to phylogeny jobs\n";

} else {
	print "Skipping Part I, all required per-MGS inputs are already prepared.\n";
	my $mergedSampleStats = recoverCompletedSplitPhaseI();
	if ($mergedSampleStats) {
		invalidateMGSInputState(@specis);
	} else {
		reportSavedSampleStats();
	}
}
loadRecoveryContributionIndex() unless $recoveryContributionIndexReady;

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
$stepStarted = time;
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
stepComplete("historical exclusion loading", $stepStarted,
	"excluded_MGS=".scalar(keys %ConspecificMGS));




my $FNAref = {}; my $FAAref = {};
my $SIgenes_OG = {}; my %OGgenesByCOG;
my %outgroupGeneCache;
my %TreeOutgroupCandidates;
my %outgroupCategoryPreflight;

# Reference data is initialized lazily after EPA-only recovery jobs have been
# submitted. Only MGS with raw staged input and no finalized outgroup overlay
# need catalogue sequences; complete published inputs already contain them.
my $requiresOutgroupReference = $runPartI || $CatNotPrepped || $repairCAT
	|| $deepRepair || $redoSubmissionData;
my %outgroupCatalogueMGS;
my $outgroupReferenceInitialized = 0;
my $outgroupReferenceCacheDir = "$scratchD/outgroup_reference_cache";
my $outgroupReferenceCacheActive = -s "$outgroupReferenceCacheDir/complete.sto" ? 1 : 0;
my $initializeOutgroupReferences = sub {
	my ($targetMGS) = @_;
	return if $outgroupReferenceInitialized;
	$outgroupReferenceInitialized = 1;
	my $referenceStarted = time;
	unless ($requiresOutgroupReference && @{$targetMGS}) {
		stepComplete("outgroup-reference preparation", $referenceStarted,
			"status=not_required", "reference_NT=0",
			"required=".($requiresOutgroupReference ? 1 : 0),
			"reference_AA=0", "MGS_with_outgroup_candidates=0");
		return;
	}

	my (%targetLoci, %candidateOutgroupsByMGS);
	my $candidateStarted = time;
	my $candidateMGS = 0;
	my $nextCandidateProgress = time + 60;
	for my $MGS (@{$targetMGS}) {
		my ($loci, $locusSource, $sampleCount) = outgroupRequirementLoci($MGS);
		next unless @{$loci};
		$targetLoci{$MGS} = $loci;
		if (defined($sampleCount)) {
			$outgroupCategoryPreflight{$MGS} = {
				loci => $loci, sample_count => $sampleCount,
			};
		}
		my @candidates;
		push @candidates, $PreferredOutgroup{$MGS}
			if exists($PreferredOutgroup{$MGS}) && length($PreferredOutgroup{$MGS});
		push @candidates, treeOutgroupCandidates($MGS) if length($treeFile);
		my %seen;
		@candidates = grep {
			defined($_) && /\A[A-Za-z0-9][A-Za-z0-9_.:+-]*\z/ && !$seen{$_}++
		} @candidates;
		$candidateOutgroupsByMGS{$MGS} = \@candidates;
		$outgroupCatalogueMGS{$_} = 1 for @candidates;
		$candidateMGS++;
		if (time >= $nextCandidateProgress) {
			stepProgress("outgroup requirement discovery", $candidateMGS,
				scalar(@{$targetMGS}), $candidateStarted,
				"candidate_MGS=".scalar(keys %outgroupCatalogueMGS),
				"current_source=$locusSource");
			$nextCandidateProgress = time + 60;
		}
	}
	unless (%outgroupCatalogueMGS && %targetLoci) {
		stepComplete("outgroup-reference preparation", $referenceStarted,
			"status=not_required", "reference_NT=0", "required=1",
			"reference_AA=0", "MGS_with_outgroup_candidates=0");
		return;
	}

	my ($refFNA, $refFAA, $refNameL) = ("", "", "unknown");
	if ($mode eq "MGS" || $mode eq "MGSall") {
		$refFNA = "$GCd/compl.incompl.$clusterID.fna";
		$refFAA = "$GCd/compl.incompl.$clusterID.prot.faa";
		$refNameL = "geneCat";
		$refFNA = resolveExistingFile($refFNA) // $refFNA;
		$refFAA = resolveExistingFile($refFAA) // $refFAA;
	} elsif ($mode eq "FMG") {
		$refFNA = "$GCd/FMG/COG*.fna";
		$refFAA = "$GCd/FMG/COG*.faa";
		$refNameL = "FMG ref";
	}
	print "Resolving selective $refNameL outgroup references for "
		.scalar(keys %targetLoci)." staged MGS and "
		.scalar(keys %outgroupCatalogueMGS)." candidate outgroup MGS; global elapsed "
		.timeNice(time - $^T)."\n";

	my @outgroupMGS = sort keys %outgroupCatalogueMGS;
	my %mapMGS = map { $_ => 1 } (@outgroupMGS, keys %targetLoci);
	my @mapMGS = sort keys %mapMGS;
	my $mapStarted = time;
	my $nextMapProgress = time + 60;
	my ($hr1, $Gene2COG_OG, $hr4);
	if (!@subsetMGS && keys(%{$SIgenes}) && keys(%{$Gene2COG}) && keys(%{$COGprios})) {
		# The ordinary all-MGS run already loaded this presorted map. Reuse it;
		# the sequence cache still retains only Mosaic/direct and same-COG references.
		($hr1, $Gene2COG_OG, $hr4) = ($SIgenes, $Gene2COG, $COGprios);
		print "Reusing the selected-MGS gene map for selective outgroup lookup; global elapsed "
			.timeNice(time - $^T)."\n";
	} else {
		my $unusedGene2MGS;
		($hr1, $Gene2COG_OG, $unusedGene2MGS, $hr4) = readGene2tax(
			$gene2taxF, $presortGenes, \@mapMGS,
			sub {
				my ($status) = @_;
				return if time < $nextMapProgress;
				stepProgress("outgroup gene-map loading", $status->{rows_scanned}, undef,
					$mapStarted, "included_genes=$status->{included_genes}");
				$nextMapProgress = time + 60;
			},
		);
	}
	$SIgenes_OG = $hr1;
	my @mappedMGS = grep { exists($outgroupCatalogueMGS{$_}) } keys %{$hr4};
	my $mappedCount = 0;
	my $nextMappingProgress = time + 60;
	for my $MGS (@mappedMGS) {
		for my $locus (@{$hr4->{$MGS}}) {
			my $gene = $hr1->{$MGS}{$locus};
			next unless defined($gene) && defined($Gene2COG_OG->{$gene});
			push @{$OGgenesByCOG{$MGS}{$Gene2COG_OG->{$gene}}}, $gene;
		}
		$mappedCount++;
		if (time >= $nextMappingProgress) {
			stepProgress("outgroup locus-index construction", $mappedCount,
				scalar(@mappedMGS), $mapStarted);
			$nextMappingProgress = time + 60;
		}
	}

	# Mosaic gives an exact target gene for many loci. Only loci without such a
	# mapping need the same-COG candidate panel and primary-gene protein used by
	# the existing similarity tie-breaker.
	my (%requiredNT, %requiredAA);
	my ($directMosaicLoci, $fallbackLoci) = (0, 0);
	for my $MGS (keys %targetLoci) {
		for my $locus (@{$targetLoci{$MGS}}) {
			my (undef, $cog, $primaryGene) = locusParts($locus, $MGS);
			next unless length($cog) && length($primaryGene);
			for my $outgroup (@{$candidateOutgroupsByMGS{$MGS} || []}) {
				if (($PreferredOutgroup{$MGS} // '') eq $outgroup
						&& exists($PreferredOutgroupGene{$MGS}{$primaryGene})) {
					my $gene = $PreferredOutgroupGene{$MGS}{$primaryGene};
					if (defined($gene) && length($gene)) {
						$requiredNT{$gene} = 1;
						$requiredAA{$gene} = 1;
						$directMosaicLoci++;
					}
					next;
				}
				my @fallback = @{$OGgenesByCOG{$outgroup}{$cog} || []};
				next unless @fallback;
				$requiredNT{$_} = 1 for @fallback;
				$requiredAA{$_} = 1 for @fallback;
				$requiredAA{$primaryGene} = 1;
				$fallbackLoci++;
			}
		}
	}
	unless (%requiredNT && %requiredAA) {
		stepComplete("outgroup-reference preparation", $referenceStarted,
			"status=no_resolvable_loci", "reference_NT=0", "required=1",
			"reference_AA=0",
			"MGS_with_outgroup_candidates=".scalar(keys %outgroupCatalogueMGS));
		return;
	}

	my ($cacheFNA, $cacheFAA) = prepareSelectiveOutgroupReferenceCache(
		$refFNA, $refFAA, \%requiredNT, \%requiredAA,
	);
	my $referenceMode = 'selective_indexed_cache';
	if (length($cacheFNA) && length($cacheFAA)) {
		$FNAref = readFasta($cacheFNA, 1, "\\s");
		$FAAref = readFasta($cacheFAA, 1, "\\s");
		$outgroupReferenceCacheActive = 1;
	} else {
		$referenceMode = 'shared_fasta_read_fallback';
		warn "Selective outgroup-reference cache was unavailable; using the shared indexed subset reader with streaming fallback\n";
		$FAAref = readFasta($refFAA, 1, "\\s", $Gene2COG_OG, { fai => 1 });
		$FNAref = readFasta($refFNA, 1, "\\s", $Gene2COG_OG, { fai => 1 });
	}
	print "Loaded ".scalar(keys %{$FNAref})." nucleotide and "
		.scalar(keys %{$FAAref})." protein outgroup reference genes from $refNameL "
		."($referenceMode; requested NT=".scalar(keys %requiredNT)
		.", AA=".scalar(keys %requiredAA).", Mosaic-direct=$directMosaicLoci, "
		."fallback-comparisons=$fallbackLoci)\n";
	stepComplete("outgroup-reference preparation", $referenceStarted,
		"status=loaded", "mode=$referenceMode",
		"reference_NT=".scalar(keys %{$FNAref}), "required=1",
		"reference_AA=".scalar(keys %{$FAAref}),
		"MGS_with_outgroup_candidates=".scalar(keys %outgroupCatalogueMGS));
};



print "\n\n----------------------------------------------------\n";
print "Part II:: submit intraStrain phylogenies for " . scalar(@specis) . " MGS. ". "Elapsed time : ", timeNice(time - $sttime) ."\n----------------------------------------------------\n\n";


die "Tree for outgroup specified, but file not found:$treeFile\nAborting..\n" if  ($treeFile ne "" && !-e $treeFile);

# Placement-only recovery has already paid for alignment and backbone
# inference. It is prepared and submitted before any catalogue-wide full-tree
# sizing or outgroup-reference initialization.
my @epaRecoveryMGS = grep { exists($MGSepaOnlyRetry{$_}) } @specis;
my @fullTreeMGS = grep { !exists($MGSepaOnlyRetry{$_}) } @specis;
my @fullTreeCandidates = grep {
	!exists($MGSsubmissionComplete{$_})
		&& !exists($MGSnoTree{$_})
		&& !(exists($ConspecificMGS{$_})
			&& ($ConspecificMGS{$_}->[0] // "") =~ m/multicopy/)
} @fullTreeMGS;
my %inputSizeByMGS;
for my $MGS (@epaRecoveryMGS) {
	my $retainedMSA = "$SIdirs{$MGS}/MSA/MSAli.fna";
	my $retainedMSASize = -s $retainedMSA;
	$retainedMSASize ||= -s "$retainedMSA.gz";
	$inputSizeByMGS{$MGS} = int(($retainedMSASize || 0) / 1024) || 1;
}
@specis = (@epaRecoveryMGS, @fullTreeMGS);
my $epaQueueBoundary = scalar(@epaRecoveryMGS);
my $fullTreeInputsInitialized = 0;
my $largestFullTreeInput = 1;
print "Validated EPA-only recovery queue: $epaOnlyRetryCount MGS; "
	."placement-only jobs will be submitted before full-tree input initialization; "
	."global elapsed ".timeNice(time - $^T)."\n"
	if $epaOnlyRetryCount;


#die;
#go through every SpecI;
$cnt=0; my $lcnt; my @jobs; my @treeJobAccounting; my %expectedTreeOutputs; my $Nspecis = @specis;
my @pendingTreeJobs;
my $submittedTreeJobs = 0;
my $nextQueuedTreeSubmissionProbe = 0;
my $treeSubmissionProbeSeconds = 20;
my %mosaicOutgroupsUsed;
my $treeMGSVisited = 0;
my %treeDisposition;
my $recalcScratchRecovered = 0;
for ($lcnt = 0; $lcnt < @specis; $lcnt++) {
	if (!$fullTreeInputsInitialized && $lcnt == $epaQueueBoundary) {
		if ($epaQueueBoundary) {
			my $epaSubmissionStarted = time;
			print "EPA-only preparation complete; draining "
				.scalar(@pendingTreeJobs)." queued placement job(s) before full-tree initialization; "
				."global elapsed ".timeNice(time - $^T)."\n";
			my $drain = dispatchPendingTreeJobs(
				queue => \@pendingTreeJobs, options => $QSBoptHR,
				jobs => \@jobs, accounting => \@treeJobAccounting,
				blocking => 1,
			);
			$submittedTreeJobs += $drain->{submitted};
			stepComplete("EPA-only recovery submission", $epaSubmissionStarted,
				"eligible=$epaQueueBoundary", "submitted=$submittedTreeJobs",
				"pending=".scalar(@pendingTreeJobs));
		}

		my $sizingStarted = time;
		print "Starting targeted full-tree input sizing for "
			.scalar(@fullTreeCandidates)." actionable MGS; global elapsed "
			.timeNice(time - $^T)."\n";
		my ($fullSizeRef, $treeInputAudit) = getInputSize(\@fullTreeCandidates);
		@inputSizeByMGS{@fullTreeCandidates} = @{$fullSizeRef};
		my $nonemptyTreeInputs = scalar(grep { $_ > 0 } @{$fullSizeRef});
		my $treeInputMB = 0;
		$treeInputMB += $_ for @{$fullSizeRef};
		$largestFullTreeInput = $_ > $largestFullTreeInput
			? $_ : $largestFullTreeInput for @{$fullSizeRef};
		my %fullTreeOrder;
		@fullTreeOrder{@fullTreeMGS} = 0 .. $#fullTreeMGS;
		my @sortedFullTreeMGS = sort {
			($inputSizeByMGS{$b} // 0) <=> ($inputSizeByMGS{$a} // 0)
				|| $fullTreeOrder{$a} <=> $fullTreeOrder{$b}
		} @fullTreeMGS;
		splice @specis, $lcnt, scalar(@fullTreeMGS), @sortedFullTreeMGS;
		stepComplete("full-tree input sizing", $sizingStarted,
			"target_MGS=".scalar(@fullTreeCandidates),
			"nonempty_inputs=$nonemptyTreeInputs",
			"estimated_uncompressed_MB=".int($treeInputMB + 0.5),
			"tooFewSamples=$treeInputAudit->{too_few_samples}",
			"noRecoverableLoci=$treeInputAudit->{no_recoverable_loci}",
			"incomplete_published=$treeInputAudit->{incomplete_published}",
			"incomplete_scratch=$treeInputAudit->{incomplete_scratch}",
			"empty_extraction=$treeInputAudit->{empty_extraction}",
			"audit=$treeInputAudit->{audit_file}");
		$initializeOutgroupReferences->(\@fullTreeCandidates);
		$fullTreeInputsInitialized = 1;
	}
	my $MGS = $specis[$lcnt]; # one per-MGS tree preparation/submission decision
	my $epaOnlyRetry = exists($MGSepaOnlyRetry{$MGS}) ? 1 : 0;
	my $epaRecovery = $epaOnlyRetry;
	if (!$recalcTrees && !$reSubmit && !$repairCAT && !$redoSubmissionData && $CatFileMiss==0 && $CatNotPrepped==0 && $treeAbsent ==0){
		$treeDisposition{'submission pass unnecessary'} += $Nspecis - $treeMGSVisited;
		print "\nAll submission dirs prepared, nothing to do..\n";
		last;
	}
	$treeMGSVisited++;
	if (!$epaRecovery && exists $MGSnoTree{$MGS}) {
		my $reason = $MGSnoTreeReason{$MGS} // 'too_few_samples';
		$treeDisposition{"valid no-tree: $reason"}++;
		limitedNotice('MGS skipped after valid no-tree classification',
			"Skipping $MGS: previous extraction recorded terminal no-tree state '$reason'.\n");
		next;
	}
	# previous condition was too lax: ( ($CatNotPrepped/$#specis) < 0.1)  , just check if we can resubmit anything here..
	if (!$epaRecovery && exists($ConspecificMGS{$MGS}) && $ConspecificMGS{$MGS}->[0] =~ m/multicopy/){
		$treeDisposition{'conspecific or multicopy'}++;
		limitedNotice('MGS skipped as conspecific or multicopy',
			"Skipping $MGS due to inclusion in conspecific MGS list.\n");next;
	}
	if ($startSubFromMGS ne "" ){
		if ($MGS ne $startSubFromMGS){
			$treeDisposition{'before requested submission start'}++;
			next;
		} else { $startSubFromMGS = "";} #deactivate now
	}
	my $outD2 = $SIdirs{$MGS};
	my $tmpD  = "$scratchD/outs/$MGS/";
	my $treeStone = "$outD2/treeDone.sto";
	my $terminalTreeMarker = "$outD2/noTree.sto";
	my $placementPendingMarker = "$outD2/placementPending.sto";
	my $IQtreef= "$outD2/phylo/IQtree_allsites.treefile";
	$IQtreef = "$outD2/phylo/VERYFASTTREE_allsites.nwk" if ($phyloProg == 2);
	$IQtreef = "$outD2/phylo/FASTTREE_allsites.nwk" if ($phyloProg == 3);
	my $publishedInputsReady = !$epaOnlyRetry
		&& !exists($legacyLocusMGS{$MGS})
		&& persistentMGSInputState($MGS) eq 'complete';
	if ($epaOnlyRetry && !prepareEpaOnlyRetryState(
			$outD2, $MGSepaOnlyRetry{$MGS})) {
		$treeDisposition{'placement completed during resume audit'}++;
		limitedNotice('EPA retry no longer required',
			"Skipping $MGS: the final non-backbone tree appeared during resume preparation.\n");
		next;
	}
	my $scratchInputsReady = 0;
	if ($recalcTrees) {
		# The sizing pass already recognizes staged inputs. Recover and combine
		# those inputs here as well, before deciding that this MGS is ineligible.
		# This also permits a complete new-format staging set to replace legacy
		# published identifiers without rerunning consensus/extraction.
		unless ($publishedInputsReady) {
			$scratchInputsReady = prepareMGSInputSet($MGS,$tmpD);
			$recalcScratchRecovered++ if $scratchInputsReady;
		}
		unless ($publishedInputsReady || $scratchInputsReady) {
			$treeDisposition{'no recoverable inputs for recalculation'}++;
			limitedWarn('MGS missing recoverable inputs for tree recalculation',
				"Skipping $MGS: -recalcTrees found neither complete published inputs nor a complete staged FNA/FAA/category set.\n");
			next;
		}
		resetMGSTreeOutputs($outD2, $MGS);
	}
	
	if (!$recalcTrees && !$reSubmit && !$repairCAT && !$redoSubmissionData && !exists($legacyLocusMGS{$MGS})
			&& -e $treeStone && -s $IQtreef ){
		$treeDisposition{'valid tree already present'}++;
		limitedNotice('MGS skipped with existing trees',
			"Skipping $MGS: a valid tree already exists.\n");
		next;
	}
	
	my $inputFNAsize = $inputSizeByMGS{$MGS} // 0;
	if ($epaRecovery && !$inputFNAsize) {
		my $retainedMSA = "$outD2/MSA/MSAli.fna";
		my $retainedMSASize = -s $retainedMSA;
		$retainedMSASize ||= -s "$retainedMSA.gz";
		$inputFNAsize = int($retainedMSASize / 1024) || 1;
	}
	#PART I: create fasta files required by tree
	make_path($outD2) unless -d $outD2;
	if ($inputFNAsize ==0 && !$epaRecovery){
		$treeDisposition{'empty input'}++;
		limitedNotice('MGS skipped with empty input', "Skipping $MGS: input is empty.\n");
		next;
	} #empty input
	my $mustRegenerateInputs = $repairCAT || $deepRepair || $redoSubmissionData
		|| exists($legacyLocusMGS{$MGS});
	if (!$epaOnlyRetry && !($publishedInputsReady && !$mustRegenerateInputs)) {
		$scratchInputsReady ||= prepareMGSInputSet($MGS,$tmpD);
		unless ($scratchInputsReady) {#$outD2); -> keep in tmpdir for now..
			$treeDisposition{'incomplete published and worker inputs'}++;
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
	
	#main command to build within species strain tree.. missing outgroup so far ($outgS)
	
	#fileGZs($FNAtf) / (1024 * 1024); #size in MB
	#$inputFNAsize*=5 if ($FNAtf =~ m/\.gz$/); #account for compressed input
	if ( 0&& ($MSAprog==4 && $inputFNAsize>700) ){ $QSBoptHR->{useLongQueue} = 1 ;	}
	my $tmpSHDD = $QSBoptHR->{tmpSpace};
	my $nodeTmpConfigured = getProgPaths("nodeTmpDir",0) ne "";
	my $treeTmpGb = int(($inputFNAsize * 4 + 1023) / 1024);
	$treeTmpGb = 15 if $treeTmpGb < 15;
	$QSBoptHR->{tmpSpace} = $nodeTmpConfigured ? $treeTmpGb : 0;
	# Placement retains likelihood vectors across the reference tree and can use
	# substantially more memory than tree inference for long concatenated MSAs.
	# Reserve a distinct scheduler profile whenever EPA-ng placement is requested.
	my $placementRequested = $strictBackbone ? 1 : 0;
	my $baseMemMult = 75; $baseMemMult = 15 if ($phyloProg ==3 || $phyloProg ==2);
	$baseMemMult = 150 if $placementRequested && $baseMemMult < 150;
	my $memoryProfile = $placementRequested ? 'EPA-ng placement' : 'tree-only';
	my $minimumMemMB = ($placementRequested ? 10240 : 5000) * $memMulti;
	$minimumMemMB = 10240 if $placementRequested && $minimumMemMB < 10240;
	my $maximumMemMB = 110000 * $memMulti;
	$maximumMemMB = $minimumMemMB if $maximumMemMB < $minimumMemMB;
	my $totMem = int($inputFNAsize * $baseMemMult * $memMulti);
	$totMem = $minimumMemMB if $totMem < $minimumMemMB;
	$totMem = $maximumMemMB if $totMem > $maximumMemMB;
	my $numCoreL = $numCores;	
	if ($maxCores >0){ #scale cores according to used memory size
		my $scaleReference = $epaRecovery
			? ($inputFNAsize || 1) : $largestFullTreeInput;
		$numCoreL = int($maxCores * sqrt($inputFNAsize/$scaleReference));
		$numCoreL = 4 if ($numCoreL < 4);		$numCoreL = $maxCores if ($numCoreL > $maxCores);
	}
	if ($epaOnlyRetry) {
		# EPA-ng retry is deliberately single-threaded. The doubled request is a
		# scheduler/cgroup allowance, not an EPA-ng --maxmem argument.
		$totMem = int($totMem * 2);
		$totMem = 20480 if $totMem < 20480;
		my $retryMaximumMemMB = int(220000 * $memMulti);
		$totMem = $retryMaximumMemMB if $totMem > $retryMaximumMemMB;
		$numCoreL = 1;
		$memoryProfile = 'EPA-ng placement-only retry';
	}
	my $iqMemMB = int($totMem * 0.9); #also supplies EPA planning-memory reporting
	

	my $bts = getProgPaths("buildTree_scr");
	my $treeFlag = "-runIQtree 1 "; 
	if ($phyloProg == 2){$treeFlag = "-runVeryFastTree 1 ";}if ($phyloProg == 3){$treeFlag = "-runFastTree 1 ";}
	my $tree_sample_separator = quotemeta($SaSe);
	my $Tcmd= "$bts -fna ".shellQuote($FNAtf)." -aa ".shellQuote($FAAtf)." -smplSep ".shellQuote($tree_sample_separator)." -cats ".shellQuote($CATtf)." -outD ".shellQuote($outD2)." $treeFlag -cores $numCoreL  ";
	$Tcmd .= "-withinSpecies 1 -relativeNTFraction $relativeNTFraction "
		."-NTfiltPerGene $GeneLengthMin -GenesPerSpecies $GenesPerSpecies "
		."-NTfiltCount $NTfiltCount -iqFast 1 ";
	$Tcmd .= "-placementGenesPerSpecies $placementGenesPerSpecies "
		if defined $placementGenesPerSpecies;
	$Tcmd .= "-placementRelativeNTFraction $placementRelativeNTFraction "
		if defined $placementRelativeNTFraction;
	$Tcmd .= "-placementNTfiltCount $placementNTfiltCount "
		if defined $placementNTfiltCount;
	$Tcmd .= "-taxonAwareLocusSelection $taxonAwareLocusSelection ";
	if ($taxonAwareLocusSelection) {
		$Tcmd .= "-taxonAwareMaxLoci $taxonAwareMaxLoci "
			."-taxonAwareCoreLoci $taxonAwareCoreLoci "
			."-taxonAwareCandidateExtra $taxonAwareCandidateExtra "
			."-taxonAwareRescueMinPrevalence $taxonAwareRescueMinPrevalence ";
		$Tcmd .= "-preferredCoreGenes ".shellQuote($preferredCoreGenes)." "
			if length($preferredCoreGenes) && !$epaOnlyRetry;
	}
	$Tcmd .= "-compactTaxonAwareDiagnostics $compactTaxonAwareDiagnostics ";
	$Tcmd .= "-rateMergePartitions $rateMergePartitions "
		."-rateMergeMaxBins $rateMergeMaxBins "
		."-rateMergeTargetSites $rateMergeTargetSites "
		."-rateMergeMinLoci $rateMergeMinLoci "
		."-rateMergeMinSites $rateMergeMinSites ";
	$Tcmd .= "-rmMSA $rmMSA -MSAprogram $MSAprog ";
	$Tcmd .= "-strictBackbone $strictBackbone "
		."-strictBackboneFraction $strictBackboneFraction "
		."-strictBackboneMinSamples $strictBackboneMinSamples "
		."-placementMinOverlap $placementMinOverlap "
		."-epaThreads ".($epaOnlyRetry ? 1 : $epaThreads)
		." -epaMaxMemMB $epaMaxMemMB "
		."-epaPendantOutlierFactor $epaPendantOutlierFactor "
		."-epaPendantMinThreshold $epaPendantMinThreshold ";
	if ($phyloProg == 1){
		$Tcmd .= "-iqMemMB $iqMemMB ";
		$Tcmd .= "-iqPathogen 1 " if $iqPathogen;
	}
	my $treeTmpOption = $nodeTmpConfigured
		? "-tmpSubdir ".shellQuote("strain_within/$MGS")
		: "-tmpD ".shellQuote("$scratchD/$MGS/");
	$Tcmd .= "$treeTmpOption -map ".shellQuote($mapF)." ";
		#die "$cmd\n" if ($cnt ==10);
	
	#if (!fileGZe($FNAtf) || !fileGZe($FAAtf) ||  ( !fileGZe($CATtf) && !-e "$CATtf.tmp") ){
	#	print "Can't find required input files:\n$FNAtf\n$FAAtf\n$CATtf\n" ;
	#	die if ($cnt <= 1);
	#	next;
	#}
	
	# Keep global catalogue-dependent outgroup selection in this controller, but
	# hand all large FASTA rewrites, category regrouping, QC finalization, sorting,
	# compression, and publication to the independently scheduled tree job.
	my $multiSmpl;my $ngenes; my $needsCopy = 0; my $inputReady = 0;
	if ($epaOnlyRetry) {
		$inputReady = 1;
	} else {
		($multiSmpl,$ngenes,$OG,$needsCopy,$inputReady)=
			addOutgroup2MGS($MGS,$OG,$tmpD);
	}
	# Locus names are MGS-qualified, so cached outgroup choices have no reuse
	# after this MGS and would otherwise accumulate for the entire submission.
	%outgroupGeneCache = ();
	$mosaicOutgroupsUsed{$MGS} = 1 if !$epaRecovery && $inputReady && exists($PreferredOutgroup{$MGS})
		&& length($OG) && $OG eq $PreferredOutgroup{$MGS};
	invalidateMGSInputState($MGS) if $inputReady;
	unless ($inputReady) {
		$treeDisposition{'input awaiting repair'}++;
		$QSBoptHR->{tmpSpace} = $tmpSHDD;
		$QSBoptHR->{useLongQueue} = 0;
		limitedNotice('MGS awaiting input repair',
			"$MGS: input files are not ready; leaving it unmarked so a repair run can retry it.\n");
		next;
	}
	
	$outgS = " -outgroup ".shellQuote($OG)." "  if ($OG ne "");
	$Tcmd .= "-sampleQC ".shellQuote("$outD2/$QCstdof")." "
		if !$epaRecovery && (fileGZe("$outD2/$QCstdof") || fileGZe("$tmpD/$QCstdof")
			|| fileGZe("$tmpD/$QCstdof.tmp") || $stagedShardHandoff{$MGS});
	$Tcmd .= "-stagedInputDir ".shellQuote($tmpD)." " if !$epaRecovery && $needsCopy;
	$Tcmd .= "-redoEPAfilter 1 " if $redoEPAfilter
		&& -s "$outD2/phylo/epa-ng/epa_result.jplace";
	$Tcmd .= "-epaOnly 1 " if $epaOnlyRetry;
	$Tcmd .= "-continue 1 -completionMarker ".shellQuote($treeStone)." "
		."-terminalMarker ".shellQuote($terminalTreeMarker)." "
		."-placementPendingMarker ".shellQuote($placementPendingMarker)." ";

	if ($epaOnlyRetry) {
		print "$MGS (".($lcnt + 1)."/$Nspecis); elapsed ".timeNice(time - $sttime)
			."; outgroup ".(length($OG) ? $OG : 'none')
			."; samples n/a; genes n/a; 1 core; $totMem MB; EPA-ng placement-only retry\n";
	} elsif ($multiSmpl > 2 && $ngenes >= $MGStoolowGsThr){
		print "$MGS (".($lcnt + 1)."/$Nspecis); elapsed ".timeNice(time - $sttime)
			."; outgroup ".(length($OG) ? $OG : 'none')
			."; $multiSmpl samples; $ngenes genes; $numCoreL cores; $totMem MB; $memoryProfile\n";
	} else {
		my $reason = $multiSmpl <= 2 ? 'too_few_samples' : 'too_few_usable_genes';
		$treeDisposition{"valid no-tree: $reason"}++;
		limitedNotice('MGS with insufficient tree input',
			"$MGS: $reason (samples=$multiSmpl, usable_genes=$ngenes); skipping tree construction\n");
		writeTooFewMarker($outD2, $multiSmpl, $ngenes, $reason);
		remove_tree($tmpD) if $needsCopy && -d $tmpD;
		$QSBoptHR->{tmpSpace} = $tmpSHDD;
		$QSBoptHR->{useLongQueue} = 0;
		next;
	}
	unlink "$outD2/tooFewSamples.sto" if -e "$outD2/tooFewSamples.sto";
	
	# PART II: retain this completely prepared tree job even if the scheduler
	# is full.  Submission is drained opportunistically below, rather than
	# blocking further category conversion and outgroup preparation.
	my $treeJobOrdinal = $cnt + 1;
	push @pendingTreeJobs, {
		mgs => $MGS,
		script => $epaOnlyRetry ? "$outD2/treeCmd.epa_retry.sh" : "$outD2/treeCmd.sh",
		command => $Tcmd.$outgS."\n",
		cores => $numCoreL,
		memory => int($totMem)."M",
		requested_mb => int($totMem),
		job_name => $epaOnlyRetry ? "EPA$treeJobOrdinal" : "FT$treeJobOrdinal",
		epa_only => $epaOnlyRetry,
		terminal => $terminalTreeMarker,
		placement_pending => $placementPendingMarker,
		tree => $IQtreef,
		stone => $treeStone,
		tmp_space => $QSBoptHR->{tmpSpace},
		use_long_queue => $QSBoptHR->{useLongQueue},
	};
	$QSBoptHR->{tmpSpace} =$tmpSHDD;
	$QSBoptHR->{useLongQueue} = 0;
	$cnt ++;
	$treeDisposition{$epaOnlyRetry ? 'EPA-only retry job' : 'eligible tree job'}++;
	$expectedTreeOutputs{$MGS} = [$IQtreef, $treeStone,
		$terminalTreeMarker, $placementPendingMarker];
	if (!$doSubmit || time >= $nextQueuedTreeSubmissionProbe) {
		my $drain = dispatchPendingTreeJobs(
			queue => \@pendingTreeJobs, options => $QSBoptHR,
			jobs => \@jobs, accounting => \@treeJobAccounting,
			blocking => 0,
		);
		$submittedTreeJobs += $drain->{submitted};
		$nextQueuedTreeSubmissionProbe = time + $treeSubmissionProbeSeconds
			if @pendingTreeJobs && $doSubmit;
	}
	#die $outD2."treeCmd.sh\n";

}
my $treeAccounted = 0;
$treeAccounted += $_ for values %treeDisposition;
print "\nTree submission accounting: $treeAccounted/$Nspecis selected MGS accounted for; "
	. "$treeMGSVisited visited by the submission loop.\n";
for my $reason (sort keys %treeDisposition) {
	print "  $reason: $treeDisposition{$reason}\n";
}
print "  staged input sets recovered for -recalcTrees: $recalcScratchRecovered\n"
	if $recalcTrees;
if ($doSubmit) {
	print "Tree preparation pass complete: $cnt eligible tree job(s), "
		."$submittedTreeJobs submitted so far, ".scalar(@pendingTreeJobs)
		." awaiting scheduler capacity.\n";
	my $drain = dispatchPendingTreeJobs(
		queue => \@pendingTreeJobs, options => $QSBoptHR,
		jobs => \@jobs, accounting => \@treeJobAccounting,
		blocking => 1,
	);
	$submittedTreeJobs += $drain->{submitted};
	die "Internal error: tree submission queue was not drained\n" if @pendingTreeJobs;
	print "Tree submission pass complete: $submittedTreeJobs eligible tree job(s) submitted; "
		.scalar(@jobs)." scheduler job ID(s) tracked. "
		."The following wait count reports jobs still present, not jobs omitted.\n";
} else {
	my $drain = dispatchPendingTreeJobs(
		queue => \@pendingTreeJobs, options => $QSBoptHR,
		jobs => \@jobs, accounting => \@treeJobAccounting,
		blocking => 1,
	);
	print "Tree submission pass complete: $cnt eligible tree script(s) generated; scheduler submission disabled.\n";
}
if (%deferredScratchCleanup) {
	my $cleanupStarted = time;
	my @cleanupPaths = sort keys %deferredScratchCleanup;
	my $cleaned = 0;
	my $nextCleanupProgress = time + 60;
	print "Starting deferred cleanup of ".scalar(@cleanupPaths)
		." completed/published scratch directories after tree submission; global elapsed "
		.timeNice(time - $^T)."\n";
	for my $path (@cleanupPaths) {
		remove_tree($path) if -d $path;
		$cleaned++;
		if (time >= $nextCleanupProgress) {
			stepProgress("deferred completed-tree scratch cleanup", $cleaned,
				scalar(@cleanupPaths), $cleanupStarted);
			$nextCleanupProgress = time + 60;
		}
	}
	stepComplete("deferred completed-tree scratch cleanup", $cleanupStarted,
		"directories=$cleaned");
}
if ($maxSubJob
		&& split_generation_complete($splitManifest, $splitStonePrefix, $maxSubJob)
		&& !splitWorkerPartsRemain()) {
	clear_split_generation($splitManifest, $splitStonePrefix);
}
#too many jobs to use as job dependency..
qsubSystemJobAlive( \@jobs,$QSBoptHR ) if @jobs && $doSubmit;
if (@treeJobAccounting) {
	retryOOMTreeJobs(
		accounting => \@treeJobAccounting,
		options => $QSBoptHR,
		maximum_mb => int($treeOOMMaxMemGB * 1024 + 0.5),
		maximum_rounds => $treeOOMRetryRounds,
	);
}
my $incompleteTreeOutcomes = 0;
if ($doSubmit) {
	my ($failed, $pending, $terminal) =
		writeTreeFailureAudit(\%expectedTreeOutputs);
	$incompleteTreeOutcomes = @{$failed} + @{$pending};
	warn "Tree jobs without a valid output were quarantined: ".join(',', @{$failed})."\n"
		if @{$failed};
	warn "Tree jobs with a retained backbone and pending placement: "
		.join(',', @{$pending})."\n" if @{$pending};
	print "Valid terminal no-tree outcomes: ".join(',', @{$terminal})."\n"
		if @{$terminal};
}
writeStrainSummary(\%treeDisposition, \%mosaicOutgroupsUsed);
my $unresolvedInputs = validateTreeInputResolution();
if ($unresolvedInputs) {
	$completionMessage = "strain_within.pl preserved completed work but stopped before downstream "
		."strain analysis; tree_outcomes_quarantined=$incompleteTreeOutcomes, "
		."tree_inputs_pending=$unresolvedInputs.";
	print "Workflow is partially complete; consult tree_job_outcomes.tsv and "
		."tree_input_resolution.tsv. No automatic full-tree resubmission was attempted.\n";
	exit(0);
}
if ($incompleteTreeOutcomes) {
	print "Tree-job outcomes remain quarantined in tree_job_outcomes.tsv, but all tree "
		."inputs are resolved; proceeding with downstream strain analysis for completed trees.\n";
}
if ($doSubmit && !$incompleteTreeOutcomes && $outgroupReferenceCacheActive
		&& length($outgroupReferenceCacheDir) && -d $outgroupReferenceCacheDir) {
	my $cacheCleanupStarted = time;
	remove_tree($outgroupReferenceCacheDir);
	stepComplete("outgroup-reference cache cleanup", $cacheCleanupStarted,
		"status=all_tree_outcomes_validated");
} elsif ($outgroupReferenceCacheActive && length($outgroupReferenceCacheDir)) {
	print "Retaining selective outgroup-reference cache until all submitted phylogenies "
		."validate successfully: $outgroupReferenceCacheDir\n";
}
print "\nAll done for $cnt Bins\nRun strain_within_2.pl for summary stats:\n";

my $MGSabundance = $MGSabundanceOverride ne ""
	? $MGSabundanceOverride
	: "$bindir/Annotation/Abundance/MGS.matL7.txt";
die "MGS abundance matrix is missing or empty: $MGSabundance\n" unless -s $MGSabundance;

my $strain2Scr = getProgPaths("MGS_strain2_scr");

my $nxtCmd = "$strain2Scr -GCd ".shellQuote($GCd)." -FMGdir ".shellQuote($outD)." -MGSmatrix ".shellQuote($MGSabundance)." -cores 4 -reSubmit 0 -DiscTests ".shellQuote($discTests)." -ContTests ".shellQuote($contTests)." -familyVar ".shellQuote($familyVar)." -groupStabilityVars ".shellQuote($groupStabilityVars)." ";
$nxtCmd .= "-MGSphylo ".shellQuote($treeFile)." " if $treeFile ne "";
$nxtCmd .= "-popGenStats $doPopGenStats ";
$nxtCmd .= "-submit $doSubmit ";
$nxtCmd .= "-qsubSystem ".shellQuote($subMode)." " if $subMode ne "";
$nxtCmd .= "-Hcores $maxCores " if $maxCores > 0;
if ($mapF2 eq ""){$nxtCmd .= "-map ".shellQuote($mapF)." ";} else {$nxtCmd .= "-map ".shellQuote($mapF2)." ";}

$nxtCmd .= "\n";

#$GCd/MB2.clusters.ext.can.Rhcl.matL0.txt
	my ($dep,$qcmd) = qsubSystem($LOGDIR."strainAnalysis2.sh",$nxtCmd,1,"60G","2StrainSub","","",1,[],$QSBoptHR);
print "\n". $nxtCmd."\n";


#cleanup
remove_tree($locTmpDir) if -d $locTmpDir;
remove_tree($preConDir) if ($preCompCons && -d $preConDir);

writeStrainWorkflowHeartbeat('complete');
if ($doSubmit && $incompleteTreeOutcomes) {
	$completionMessage = "strain_within.pl started downstream strain analysis from completed trees; "
		."$incompleteTreeOutcomes tree outcome(s) remain quarantined for inspection.";
} else {
	$completionMessage = "strain_within.pl completed normally; $cnt eligible tree job(s) were "
		. ($doSubmit ? "submitted and validated." : "generated without scheduler submission.");
}
exit(0);

 

#########################################################################################
#########################################################################################


sub readSplitGeneration {
	return 'unsplit' unless $maxSubJob;
	open my $input, '<', $splitManifest
		or die "Cannot read completed split generation $splitManifest: $!\n";
	my $line = <$input> // '';
	close $input or die "Cannot close completed split generation $splitManifest: $!\n";
	$line =~ s/[\r\n]+\z//;
	die "Malformed completed split generation $splitManifest\n"
		unless $line =~ /^([A-Za-z0-9_.:-]+)\t(\d+)$/ && $2 == $maxSubJob;
	my $generation = $1;
	die "Split generation is no longer complete: $splitManifest\n"
		unless split_generation_complete($splitManifest, $splitStonePrefix, $maxSubJob);
	return $generation;
}

sub collectMGSShardHandoff {
	my ($MGS, $tmpD) = @_;
	my @types = (
		['fna', $FNAstdof], ['faa', $FAAstdof], ['link', $LINKstdof],
		['category', "$CATstdof.tmp"], ['qc', "$QCstdof.tmp"],
	);
	my $workerCount = $maxSubJob || 1;
	my (%parts, %workerPart);
	for my $type (@types) {
		my ($label, $name) = @{$type};
		my @found = exact_worker_parts("$tmpD/$name", $workerCount);
		$parts{$label} = \@found;
		for my $path (@found) {
			die "Cannot determine worker suffix for $path\n" unless $path =~ /\.(\d+)\z/;
			$workerPart{$1}{$label} = $path;
		}
	}
	return unless grep { @{$parts{$_}} } map { $_->[0] } @types;
	return unless !$maxSubJob
		|| split_generation_complete($splitManifest, $splitStonePrefix, $maxSubJob);
	return unless $recoveryContributionIndexReady;
	my @expectedWorkers = sort { $a <=> $b }
		keys %{$recoveryWorkersByMGS{$MGS} || {}};
	return unless @expectedWorkers;
	my %expected = map { $_ => 1 } @expectedWorkers;
	for my $type (@types) {
		my ($label) = @{$type};
		my @actual = sort { $a <=> $b } map {
			/\.(\d+)\z/ ? 0 + $1 : ()
		} @{$parts{$label}};
		my @missing = grep { !exists $workerPart{$_}{$label} } @expectedWorkers;
		my @unexpected = grep { !$expected{$_} } @actual;
		if (@missing || @unexpected) {
			limitedWarn('incomplete worker shard handoff',
				"Cannot hand $MGS/$label to buildTree5: missing="
				.(@missing ? join(',', @missing) : 'none')." unexpected="
				.(@unexpected ? join(',', @unexpected) : 'none')."\n");
			return;
		}
	}
	my @workers;
	my ($workerRows, $workerRecords) = (0, 0);
	for my $worker (@expectedWorkers) {
		my %workerParts;
		for my $type (@types) {
			my ($label) = @{$type};
			my $path = $workerPart{$worker}{$label};
			my @stat = stat($path);
			return unless @stat && $stat[7] > 0;
			$workerParts{$label} = {
				path => $path, basename => basename($path), bytes => $stat[7],
			};
		}
		my $rows = $recoveryWorkerRowsByMGS{$MGS}{$worker} // 0;
		my $records = $recoveryWorkerRecordsByMGS{$MGS}{$worker} // 0;
		$workerRows += $rows;
		$workerRecords += $records;
		push @workers, {
			id => $worker,
			rows => $rows,
			records => $records,
			parts => \%workerParts,
		};
	}
	my $expectedRecords = $recoveryRecordsByMGS{$MGS} // 0;
	my $expectedSamples = $recoveryUniqueSamplesByMGS{$MGS}
		// scalar(keys %{$recoverySamplesByMGS{$MGS} || {}});
	die "Recovery contribution totals changed while preparing $MGS: "
		."worker rows=$workerRows expected samples=$expectedSamples; "
		."worker records=$workerRecords expected records=$expectedRecords\n"
		unless $workerRows == $expectedSamples && $workerRecords == $expectedRecords;
	return {
		generation => readSplitGeneration(), workers => \@workers,
		expected_records => $expectedRecords,
		expected_ingroup_samples => $expectedSamples,
	};
}

sub prepareMGSInputSet {
	my ($MGS, $tmpD) = @_;
	if (my $handoff = collectMGSShardHandoff($MGS, $tmpD)) {
		$stagedShardHandoff{$MGS} = $handoff;
		return 1;
	}
	delete $stagedShardHandoff{$MGS};
	return combineMGSgenesDir($MGS, $tmpD);
}

sub writeMGSShardManifest {
	my ($path, $handoff, $MGS, $OG, $loci, $ingroupSamples, $writer) = @_;
	die "Cannot write a shard manifest without a validated handoff for $MGS\n"
		unless $handoff && @{$handoff->{workers} || []};
	die "Shard/category sample count mismatch for $MGS: category=$ingroupSamples manifest=$handoff->{expected_ingroup_samples}\n"
		unless $ingroupSamples == $handoff->{expected_ingroup_samples};
	my @line = ('strain-shard-input-v1');
	push @line,
		join("\t", 'value', 'mgs', $MGS),
		join("\t", 'value', 'outgroup', $OG),
		join("\t", 'value', 'generation', $handoff->{generation}),
		join("\t", 'value', 'separator', $SaSe),
		join("\t", 'value', 'expected_records', $handoff->{expected_records}),
		join("\t", 'value', 'expected_ingroup_samples', $ingroupSamples),
		join("\t", 'value', 'expected_loci', $loci),
		join("\t", 'output', 'fna', $FNAstdof),
		join("\t", 'output', 'faa', $FAAstdof),
		join("\t", 'output', 'link', $LINKstdof),
		join("\t", 'output', 'category', $CATstdof),
		join("\t", 'output', 'qc', $QCstdof),
		join("\t", 'output', 'data_log', 'data.log');
	for my $worker (@{$handoff->{workers}}) {
		push @line, join("\t", 'worker', $worker->{id},
			$worker->{rows}, $worker->{records});
		for my $type (qw(fna faa link category qc)) {
			my $part = $worker->{parts}{$type};
			push @line, join("\t", 'part', $worker->{id}, $type,
				$part->{basename}, $part->{bytes});
		}
	}
	$writer->($path, join("\n", @line)."\n", 'staged worker-shard manifest');
}

sub combineMGSgenesDir{
	my ($MGS,$tmpD) = @_;
	my @coreRequired = (
		"$tmpD/$FNAstdof", "$tmpD/$FAAstdof", "$tmpD/$LINKstdof",
	);
	my $mergeCheckpoint = "$tmpD/merge.complete.tsv";
	my $aggregateComplete = !grep { !fileGZe($_) } @coreRequired;
	$aggregateComplete &&= (
		(fileGZe("$tmpD/$CATstdof.tmp") && fileGZe("$tmpD/$QCstdof.tmp"))
		|| (fileGZe("$tmpD/$CATstdof") && fileGZe("$tmpD/$QCstdof"))
	);
	$aggregateComplete &&= -s $mergeCheckpoint;
	my @filesets = (
		[$FNAstdof,       "$tmpD/$FNAstdof",       "$tmpD/$FNAstdof"],
		[$FAAstdof,       "$tmpD/$FAAstdof",       "$tmpD/$FAAstdof"],
		[$LINKstdof,      "$tmpD/$LINKstdof",      "$tmpD/$LINKstdof"],
		["$CATstdof.tmp", "$tmpD/$CATstdof.tmp", "$tmpD/$CATstdof.tmp"],
		["$QCstdof.tmp",  "$tmpD/$QCstdof.tmp",  "$tmpD/$QCstdof.tmp"],
	);
	my (%partsByName, %partByWorker);
	my $workerCount = $maxSubJob || 1;
	for my $set (@filesets) {
		my ($name, $prefix) = @$set;
		my @parts = exact_worker_parts($prefix, $workerCount);
		$partsByName{$name} = \@parts;
		for my $part (@parts) {
			die "Cannot determine worker suffix for $part\n" unless $part =~ /\.(\d+)\z/;
			$partByWorker{$name}{$1} = $part;
		}
	}
	my $hasFreshParts = grep { @{$partsByName{$_}} }
		($FNAstdof, $FAAstdof, $LINKstdof, "$CATstdof.tmp", "$QCstdof.tmp");
	return $aggregateComplete unless $hasFreshParts;
	if ($maxSubJob && !split_generation_complete($splitManifest, $splitStonePrefix, $maxSubJob)) {
		limitedWarn('partial worker retries without a complete generation',
			"Ignoring partial worker retry for $MGS: no complete matching split-extraction generation is present\n");
		return $aggregateComplete;
	}

	unless ($recoveryContributionIndexReady) {
		limitedWarn('missing recovery contribution index',
			"Rejecting fresh merge for $MGS: no recovery contribution index is available; retaining worker parts\n");
		return $aggregateComplete;
	}
	my @expectedWorkers = sort { $a <=> $b }
		keys %{$recoveryWorkersByMGS{$MGS} || {}};
	unless (@expectedWorkers) {
		limitedWarn('worker parts without recovery rows',
			"Rejecting worker parts for $MGS: the recovery ledger has no recovered samples for this MGS\n");
		return $aggregateComplete;
	}
	my %expected = map { $_ => 1 } @expectedWorkers;
	my @contributorNames = (
		$FNAstdof, $FAAstdof, $LINKstdof, "$CATstdof.tmp", "$QCstdof.tmp",
	);
	for my $name (@contributorNames) {
		my @actual = sort { $a <=> $b } keys %{$partByWorker{$name} || {}};
		my @missing = grep { !exists $partByWorker{$name}{$_} } @expectedWorkers;
		my @unexpected = grep { !$expected{$_} } @actual;
		if (@missing || @unexpected) {
			my $details = 'missing='.(@missing ? join(',', @missing) : 'none')
				.' unexpected='.(@unexpected ? join(',', @unexpected) : 'none');
			limitedWarn('incomplete worker contribution set',
				"Rejecting merge for $MGS/$name: $details; retaining worker parts\n");
			return $aggregateComplete;
		}
	}

	my (%mergeFileByName, %observedRows, %identifierDigest);
	my @consumedParts;
	for my $set (@filesets) {
		my ($name, undef, $outfile) = @$set;
		my @parts = @{$partsByName{$name}};
		next unless @parts;
		my $mergeFile = "$outfile.merge.$$";
		open my $out, '>', $mergeFile or die "Cannot create $mergeFile: $!\n";
		binmode $out;
		my $rows = 0;
		my $digest = $name eq "$QCstdof.tmp" ? undef : Digest::SHA->new(256);
		for my $file (@parts) {
			open my $part_fh, '<', $file or die "Cannot read $file: $!\n";
			binmode $part_fh;
			while (my $line = <$part_fh>) {
				print {$out} $line or die "Cannot write $mergeFile: $!\n";
				my $identifier;
				if ($name eq $FNAstdof || $name eq $FAAstdof) {
					if ($line =~ /^>(\S+)/) {
						$rows++;
						$identifier = $1;
					}
				} elsif ($line =~ /\S/) {
					$rows++;
					my @field = split /\t/, $line, -1;
					if ($name eq "$CATstdof.tmp") {
						die "Malformed category worker row in $file\n" unless @field >= 4 && length $field[3];
						$identifier = $field[3];
					} elsif ($name eq $LINKstdof) {
						die "Malformed link worker row in $file\n" unless @field && length $field[0];
						$identifier = $field[0];
					}
				}
				$identifier =~ s/[\r\n]+\z// if defined $identifier;
				$digest->add($identifier, "\n") if $digest && defined $identifier;
			}
			close $part_fh or die "Cannot close $file: $!\n";
		}
		close $out or die "Cannot close $mergeFile: $!\n";
		$mergeFileByName{$name} = $mergeFile;
		$observedRows{$name} = $rows;
		$identifierDigest{$name} = $digest->hexdigest if $digest;
		push @consumedParts, @parts;
	}

	my @validationErrors;
	my $fnaRows = $observedRows{$FNAstdof} // -1;
	my $faaRows = $observedRows{$FAAstdof} // -1;
	my $catRows = $observedRows{"$CATstdof.tmp"} // -1;
	push @validationErrors, "FNA=$fnaRows FAA=$faaRows category=$catRows"
		unless $fnaRows >= 0 && $fnaRows == $faaRows && $fnaRows == $catRows;
	my $linkRows = $observedRows{$LINKstdof} // -1;
	push @validationErrors, "links=$linkRows FNA=$fnaRows" unless $linkRows == $fnaRows;
	for my $name ($FAAstdof, "$CATstdof.tmp", $LINKstdof) {
		push @validationErrors, "identifier order differs: $FNAstdof vs $name"
			unless ($identifierDigest{$FNAstdof} // '') eq ($identifierDigest{$name} // '');
	}
	if ($recoveryContributionIndexReady) {
		my $expectedRecords = $recoveryRecordsByMGS{$MGS} // 0;
		my $expectedRows = $recoveryRowsByMGS{$MGS} // 0;
		my $uniqueSamples = $recoveryUniqueSamplesByMGS{$MGS}
			// scalar keys %{$recoverySamplesByMGS{$MGS} || {}};
		push @validationErrors, "records=$fnaRows expected=$expectedRecords"
			unless $fnaRows == $expectedRecords;
		push @validationErrors, "QC=".($observedRows{"$QCstdof.tmp"} // -1)." expected=$expectedRows"
			unless ($observedRows{"$QCstdof.tmp"} // -1) == $expectedRows;
		push @validationErrors, "recovery_rows=$expectedRows unique_samples=$uniqueSamples"
			unless $expectedRows == $uniqueSamples;
		if (exists $mergeFileByName{$LINKstdof}) {
			push @validationErrors, "links=$observedRows{$LINKstdof} expected=$expectedRecords"
				unless $observedRows{$LINKstdof} == $expectedRecords;
		}
	}
	if (@validationErrors) {
		retry_unlink($_, fatal => 0, label => "clean rejected worker merge") for values %mergeFileByName;
		limitedWarn('worker merge count mismatch',
			"Rejecting merge for $MGS: ".join('; ', @validationErrors)."; retaining worker parts\n");
		return $aggregateComplete;
	}

	# Invalidate the previous commit only after the complete retry has passed all
	# contributor, cardinality, and identifier-order checks. A validation error or
	# crash while constructing merge files therefore leaves the old aggregate usable.
	retry_unlink($mergeCheckpoint, label => "invalidate stale merge checkpoint");
	# A completed earlier Phase II may have finalised scratch category/QC sidecars.
	# A fresh worker generation supersedes them, so do not reuse stale final inputs.
	for my $preparedSidecar ("$tmpD/$CATstdof", "$tmpD/$QCstdof") {
		retry_unlink($preparedSidecar, fatal => 0,
			label => "remove stale prepared scratch sidecar");
	}
	retry_unlink("$tmpD/data.log", fatal => 0,
		label => "remove stale prepared outgroup log");
	retry_unlink("$tmpD/data.log.gz", fatal => 0,
		label => "remove stale compressed prepared outgroup log");
	for my $set (@filesets) {
		my ($name, undef, $outfile) = @$set;
		next unless exists $mergeFileByName{$name};
		retry_rename($mergeFileByName{$name}, $outfile,
			label => "publish merged MGS input $outfile");
	}
	my $checkpointTemporary = "$mergeCheckpoint.write.$$";
	my $checkpointFH = retry_open('>', $checkpointTemporary,
		label => 'create MGS merge checkpoint');
	print {$checkpointFH} join("\t", "strain-merge-v1", $MGS,
		$recoveryRowsByMGS{$MGS} // 0, $fnaRows, $identifierDigest{$FNAstdof} // ''), "\n"
		or die "Cannot write $checkpointTemporary: $!\n";
	retry_close($checkpointFH, 'close MGS merge checkpoint');
	retry_rename($checkpointTemporary, $mergeCheckpoint,
		label => 'publish MGS merge checkpoint');
	my $complete = !grep { !fileGZe($_) } @coreRequired;
	$complete &&= (
		fileGZe("$tmpD/$CATstdof.tmp")
		&& fileGZe("$tmpD/$QCstdof.tmp")
	);
	$complete &&= -s $mergeCheckpoint;
	if ($complete) {
		for my $part (@consumedParts) {
			retry_unlink($part, fatal => 0, label => "clean combined worker part");
		}
		if ($recoveryContributionIndexReady) {
			warn "Merged $MGS: workers=".join(',', @expectedWorkers)
				.", samples=$recoveryRowsByMGS{$MGS}, records=$recoveryRecordsByMGS{$MGS}\n";
		}
	} else {
		limitedWarn('combined MGS inputs missing source parts',
			"Incomplete combined input for $MGS; retaining all source parts for repair\n");
	}
	return $complete;
}

sub splitWorkerPartsRemain {
	return 0 unless $maxSubJob;
	for my $mgs_dir (grep { -d $_ } bsd_glob("$scratchD/outs/*")) {
		for my $name ($FNAstdof, $FAAstdof, $LINKstdof, "$CATstdof.tmp", "$QCstdof.tmp") {
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
	if (defined($default_mgs)
		&& ($PreferredOutgroup{$default_mgs} // '') eq $outgroup
		&& exists($PreferredOutgroupGene{$default_mgs}{$primary_gene})) {
		return $outgroupGeneCache{$cache_key}
			= $PreferredOutgroupGene{$default_mgs}{$primary_gene};
	}
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

sub treeOutgroupCandidates {
	my ($MGS) = @_;
	return @{$TreeOutgroupCandidates{$MGS}}
		if exists($TreeOutgroupCandidates{$MGS});
	my @candidates;
	if (defined($treeFile) && length($treeFile) && -e $treeFile) {
		my $neiTree = getProgPaths("neighborTree");
		my $call = "$neiTree ".shellQuote($treeFile)." ".shellQuote($MGS);
		my $outgroup_text = `$call`;
		if ($? != 0) {
			limitedWarn('outgroup lookup command failures',
				"Can't find outgroup from call $call; trying catalogue-derived candidates\n");
		} else {
			@candidates = split /\s+/, $outgroup_text;
		}
	}
	my %seen;
	@candidates = grep {
		/\A[A-Za-z0-9][A-Za-z0-9_.:+-]*\z/ && !$seen{$_}++
	} @candidates;
	$TreeOutgroupCandidates{$MGS} = \@candidates;
	return @candidates;
}

sub addOutgroup2MGS{
	my ($MGS,$OG,$tmpD) = @_;
	my $outD2 = $SIdirs{$MGS};
	my $shardHandoff = $stagedShardHandoff{$MGS};
	my $outputReady = fileGZe("$outD2/$FNAstdof")
		&& fileGZe("$outD2/$FAAstdof") && fileGZe("$outD2/$CATstdof");
	my ($publishedPrepared, $publishedOG) = preparedOutgroupLog($outD2);
	if ($outputReady && $publishedPrepared && !$repairCAT && !$deepRepair && !$redoSubmissionData
			&& !exists($legacyLocusMGS{$MGS})){
		my (%samplesSeen, $genesSeen);
		my ($catFh) = gzipopen("$outD2/$CATstdof", "existing category file");
		while (my $line = <$catFh>) {
			chomp $line;
			next unless length($line);
			$genesSeen++;
			for my $entry (split /\t/, $line) {
				my ($sample) = split /\Q$SaSe\E/, $entry, 2;
				$samplesSeen{$sample} = 1 if defined($sample) && length($sample);
			}
		}
		close $catFh or die "Cannot close existing category file for $MGS: $!\n";
		return (scalar(keys %samplesSeen), $genesSeen, $publishedOG, 0, 1);
	}

	# Compatibility for controller runs that had already completed the old
	# controller-side Phase II before this version was installed.
	my $preparedScratchInput = fileGZe("$tmpD/$FNAstdof")
		&& fileGZe("$tmpD/$FAAstdof") && fileGZe("$tmpD/$LINKstdof")
		&& fileGZe("$tmpD/$CATstdof") && fileGZe("$tmpD/$QCstdof")
		&& -s "$tmpD/merge.complete.tsv";
	my ($scratchPrepared, $preparedOG) = preparedOutgroupLog($tmpD);
	if (!$shardHandoff && $preparedScratchInput && $scratchPrepared && !$repairCAT && !$deepRepair && !$redoSubmissionData
			&& !exists($legacyLocusMGS{$MGS})) {
		my (%samplesSeen, $genesSeen);
		my ($catFh) = gzipopen("$tmpD/$CATstdof", "prepared scratch category file");
		while (my $line = <$catFh>) {
			chomp $line;
			next unless length($line);
			$genesSeen++;
			for my $entry (split /\t/, $line) {
				my ($sample) = split /\Q$SaSe\E/, $entry, 2;
				$samplesSeen{$sample} = 1 if defined($sample) && length($sample);
			}
		}
		close $catFh or die "Cannot close prepared scratch category file for $MGS: $!\n";
		return (scalar(keys %samplesSeen), $genesSeen, $preparedOG, 1, 1);
	}

	my $rawCategory = "$tmpD/$CATstdof.tmp";
	my @rawCategorySources = $shardHandoff
		? map { $_->{parts}{category}{path} } @{$shardHandoff->{workers}}
		: ($rawCategory);
	my $stageReady = $shardHandoff ? scalar(@rawCategorySources)
		: fileGZe("$tmpD/$FNAstdof") && fileGZe("$tmpD/$FAAstdof")
			&& fileGZe($rawCategory) && fileGZe("$tmpD/$QCstdof.tmp")
			&& -s "$tmpD/merge.complete.tsv";
	if (!$stageReady) {
		limitedWarn('MGS missing raw staged tree input',
			"$MGS has no complete raw staged FNA/FAA/category/QC input in $tmpD; leaving it for repair\n");
		return (0, 0, $OG, 0, 0);
	}

	# Reuse a raw-category requirement prepass when one was necessary. Normal
	# resumes derive requirements from the selected gene map and retain this
	# just-in-time streaming scan for exact locus and sample counts.
	my (%locusSeen, %sampleSeen);
	my $ingroupSampleCount;
	if (my $preflight = delete $outgroupCategoryPreflight{$MGS}) {
		$locusSeen{$_} = 1 for @{$preflight->{loci} || []};
		$ingroupSampleCount = $preflight->{sample_count} // 0;
	} else {
		my $categoryScanStarted = time;
		my $rawCategoryRows = 0;
		my $nextCategoryProgress = time + 60;
		for my $categorySource (@rawCategorySources) {
			my ($rawFh) = gzipopen($categorySource, "raw staged category input", 1);
			while (my $line = <$rawFh>) {
				$line =~ s/[\r\n]+\z//;
				next unless length($line);
				$rawCategoryRows++;
				if (time >= $nextCategoryProgress) {
					stepProgress("staged category scan for $MGS", $rawCategoryRows, undef,
						$categoryScanStarted, "loci=".scalar(keys %locusSeen),
						"samples=".scalar(keys %sampleSeen));
					$nextCategoryProgress = time + 60;
				}
				my @fields = split /\t/, $line, -1;
				die "Malformed raw staged category row for $MGS: $line\n"
					unless @fields >= 4 && $fields[0] eq $MGS
						&& length($fields[1]) && length($fields[2]) && length($fields[3]);
				$locusSeen{$fields[1]} = 1;
				$sampleSeen{$fields[2]} = 1;
			}
			close $rawFh or die "Cannot close raw staged category input $categorySource: $!\n";
		}
		$ingroupSampleCount = scalar keys %sampleSeen;
	}
	my @curCogs = sort keys %locusSeen;
	if (@curCogs < $MGStoolowGsThr) {
		limitedWarn('MGS with too few usable genes for tree construction',
			"$MGS has only ".scalar(@curCogs)." usable genes; skipping tree construction\n");
		return ($ingroupSampleCount, scalar(@curCogs), $OG, 1, 1);
	}

	if ($treeFile ne "" || exists($PreferredOutgroup{$MGS})) {
		my @candidates;
		push @candidates, $PreferredOutgroup{$MGS}
			if exists($PreferredOutgroup{$MGS}) && length($PreferredOutgroup{$MGS});
		push @candidates, treeOutgroupCandidates($MGS) if $treeFile ne "";
		my %seenCandidate;
		@candidates = grep { !$seenCandidate{$_}++ } @candidates;
		$OG = "";
		limitedWarn('MGS without outgroup candidates',
			"No outgroup candidates returned for $MGS; building an ingroup-only tree\n")
			if @candidates == 0;
		my $represented = 0;
		for my $candidate (@candidates) {
			$represented = 0;
			$OG = $candidate;
			next unless exists($SIgenes_OG->{$OG});
			for my $locus (@curCogs) {
				my (undef, $annotation) = locusParts($locus, $MGS);
				next if $annotation =~ m/^uniq\d+$/;
				my $outgroupGene = outgroupGeneForLocus($OG, $locus, $MGS);
				next unless length($outgroupGene) && exists($FNAref->{$outgroupGene});
				$represented++;
			}
			last if $represented >= $MGStoolowGsThr;
		}
		if ($represented < $MGStoolowGsThr) {
			my @preview = @curCogs[0 .. ($#curCogs < 9 ? $#curCogs : 9)];
			limitedWarn('MGS without a sufficiently represented outgroup',
				"Could not find a sufficiently represented outgroup for $MGS; candidates: @candidates; loci: @preview\n");
			$OG = "";
		}
		if ($OG ne "" && !exists($SIgenes_OG->{$OG})) {
			limitedWarn('selected outgroups absent from gene catalogue',
				"Selected outgroup $OG for $MGS is absent from the gene catalogue\n");
			$OG = "";
		}
	}

	my ($overlayFNA, $overlayFAA, $overlayCategory, $outgroupGenes) = ('', '', '', 0);
	if ($OG ne '') {
		for my $locus (@curCogs) {
			my (undef, $annotation) = locusParts($locus, $MGS);
			next if $annotation =~ m/^uniq\d+$/;
			my $gene = outgroupGeneForLocus($OG, $locus, $MGS);
			next unless length($gene) && exists($FNAref->{$gene}) && exists($FAAref->{$gene});
			my $identifier = "$OG$SaSe" . externalLocusName($locus, $MGS);
			$overlayFNA .= ">$identifier\n$FNAref->{$gene}\n";
			$overlayFAA .= ">$identifier\n$FAAref->{$gene}\n";
			$overlayCategory .= join("\t", $locus, $OG, $identifier)."\n";
			$outgroupGenes++;
		}
		if (!$outgroupGenes) {
			limitedWarn('MGS with no usable outgroup genes',
				"Couldn't stage any outgroup genes for $OG; building $MGS without an outgroup\n");
			$OG = '';
		}
	}
	my $treeSampleCount = $ingroupSampleCount + ($OG ne '' && $outgroupGenes ? 1 : 0);

	# The plan and small overlays are the only Phase-II outputs written by the
	# controller.  They contain already selected reference sequences, so the tree
	# job needs neither gene2tax nor any other large catalogue index.
	for my $stale (bsd_glob("$tmpD/.strain_tree_input.*")) {
		retry_unlink($stale, fatal => 0, label => "clear stale staged tree overlay $stale");
	}
	my $writeOverlay = sub {
		my ($path, $contents, $label) = @_;
		my $temporary = "$path.write.$$";
		open my $output, '>', $temporary or die "Cannot create $label $temporary: $!\n";
		print {$output} $contents or die "Cannot write $label $temporary: $!\n";
		close $output or die "Cannot close $label $temporary: $!\n";
		retry_rename($temporary, $path, label => "publish $label $path");
	};
	$writeOverlay->("$tmpD/.strain_tree_input.outgroup.fna", $overlayFNA,
		'staged outgroup nucleotide overlay') if length($overlayFNA);
	$writeOverlay->("$tmpD/.strain_tree_input.outgroup.faa", $overlayFAA,
		'staged outgroup protein overlay') if length($overlayFAA);
	$writeOverlay->("$tmpD/.strain_tree_input.outgroup.cat.tsv", $overlayCategory,
		'staged outgroup category overlay') if length($overlayCategory);
	$writeOverlay->("$tmpD/.strain_tree_input.plan.tsv",
		"strain-staged-input-v1\noutgroup\t$OG\nmgs\t$MGS\n", 'staged tree-input plan');
	writeMGSShardManifest("$tmpD/.strain_tree_input.shards.tsv",
		$shardHandoff, $MGS, $OG, scalar(@curCogs), $ingroupSampleCount, $writeOverlay)
		if $shardHandoff;
	return ($treeSampleCount, scalar(@curCogs), $OG, 1, 1);
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

	my @snp_parts = grep { -e $_ } map {
		my $current = "$LOGDIR/SNPconsCalls.$_.log";
		my $legacy = "$outD/SNPconsCalls.$_.log";
		-e $current ? $current : $legacy;
	} 0 .. $maxSubJob - 1;
	if (@snp_parts) {
		my $snp_log = "$LOGDIR/SNPconsCalls.log";
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
	my ($outD2, $sampleCount, $geneCount, $reason) = @_;
	make_path($outD2) unless -d $outD2;
	my $marker = "$outD2/tooFewSamples.sto";
	open my $out, '>', $marker or die "Cannot create $marker: $!\n";
	print {$out} "reason\t".($reason // 'too_few_samples')."\nsamples\t$sampleCount\ngenes\t$geneCount\n"
		or die "Cannot write $marker: $!\n";
	close $out or die "Cannot close $marker: $!\n";
}

sub writeNoRecoverableLociMarker {
	my ($outD2, $reason) = @_;
	make_path($outD2) unless -d $outD2;
	my $marker = "$outD2/noRecoverableLoci.sto";
	my $temporary = "$marker.write.$$";
	open my $out, '>', $temporary or die "Cannot create $temporary: $!\n";
	print {$out} "reason\t".($reason // 'empty_extraction')."\n"
		or die "Cannot write $temporary: $!\n";
	close $out or die "Cannot close $temporary: $!\n";
	rename $temporary, $marker or die "Cannot publish $marker: $!\n";
}

sub recordValidatedEmptyExtractions {
	my ($targets) = @_;
	return 0 unless ref($targets) eq 'ARRAY';
	my $recorded = 0;
	for my $MGS (@{$targets}) {
		next if exists $MGSnoTree{$MGS};
		next unless persistentMGSInputState($MGS) eq 'missing';
		next if scratchMGSInputState($MGS) ne 'missing';
		# This is called only after every worker in the current Stage-I generation
		# has completed and its ledgers have been merged.  With no partial files
		# and no aggregate triplet left, no sample supplied a recoverable locus.
		writeNoRecoverableLociMarker($SIdirs{$MGS}, 'empty_extraction');
		$MGSnoTree{$MGS} = 1;
		$MGSnoTreeReason{$MGS} = 'no_recoverable_loci';
		$recorded++;
	}
	return $recorded;
}

sub validateTreeInputResolution {
	my $audit = "$LOGDIR/tree_input_resolution.tsv";
	my $temporary = "$audit.write.$$";
	make_path($LOGDIR) unless -d $LOGDIR;
	my $out = retry_open('>', $temporary, label => 'create tree-input resolution audit');
	print {$out} join("\t", qw(MGS resolution persistent_state scratch_state reason)), "\n"
		or die "Cannot write $temporary: $!\n";
	my (@repairRequired, %repairState, $ready, $terminal, $excluded);
	for my $MGS (@specis) {
		my $persistent = persistentMGSInputState($MGS);
		my $scratch = scratchMGSInputState($MGS);
		my ($resolution, $reason) = ('', '');
		if (-s "$SIdirs{$MGS}/tooFewSamples.sto") {
			$resolution = 'valid_no_tree_too_few_samples';
			$terminal++;
		} elsif (-s "$SIdirs{$MGS}/noRecoverableLoci.sto") {
			$resolution = 'valid_no_tree_no_recoverable_loci';
			$terminal++;
		} elsif (-s "$SIdirs{$MGS}/noTree.sto") {
			$resolution = 'valid_no_tree_buildtree';
			$reason = lifecycleMarkerReason("$SIdirs{$MGS}/noTree.sto",
				'buildtree_no_usable_alignment');
			$terminal++;
		} elsif (exists($ConspecificMGS{$MGS}) && $ConspecificMGS{$MGS}->[0] =~ /multicopy/) {
			$resolution = 'excluded_conspecific_or_multicopy';
			$excluded++;
		} elsif ($persistent eq 'complete' || $scratch eq 'complete') {
			$resolution = 'tree_input_ready';
			$ready++;
		} else {
			$resolution = 'repair_required';
			push @repairRequired, $MGS;
			$repairState{$MGS} = [$persistent, $scratch];
		}
		print {$out} join("\t", $MGS, $resolution, $persistent, $scratch, $reason), "\n"
			or die "Cannot write $temporary: $!\n";
	}
	retry_close($out, 'close tree-input resolution audit');
	retry_rename($temporary, $audit, label => 'publish tree-input resolution audit');
	print "Tree-input resolution audit: ready=$ready, valid_no_tree=$terminal, "
		."excluded=$excluded, repair_required=".scalar(@repairRequired)."; $audit\n";
	my $repairQueue = "$LOGDIR/tree_input_repair.queue.tsv";
	if (@repairRequired) {
		my $queueTemporary = "$repairQueue.write.$$";
		my $queue = retry_open('>', $queueTemporary, label => 'create tree-input repair queue');
		print {$queue} "MGS\tpersistent_state\tscratch_state\n";
		for my $MGS (@repairRequired) {
			print {$queue} join("\t", $MGS, @{$repairState{$MGS}}), "\n"
				or die "Cannot write tree-input repair queue: $!\n";
		}
		retry_close($queue, 'close tree-input repair queue');
		retry_rename($queueTemporary, $repairQueue,
			label => 'publish tree-input repair queue');
		warn "Tree-input repair remains queued for ".scalar(@repairRequired)
			." MGS; completed trees were retained and no catalogue-wide abort was triggered\n";
	} else {
		retry_unlink($repairQueue, fatal => 0, label => 'clear obsolete tree-input repair queue');
	}
	return scalar(@repairRequired);
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
		# Partition whole assembly groups, never individual samples.  The
		# catalogue has one shared-reference driver, while every member still
		# needs its sample-specific VCF/depth consensus below.
		my (%samplesByGroup, %groupForSample);
		for my $sample (@samples) {
			my $group = defined($map{$sample}{AssGroup}) && $map{$sample}{AssGroup} ne '-1'
				? $map{$sample}{AssGroup} : "__standalone__${sample}";
			push @{$samplesByGroup{$group}}, $sample;
			$groupForSample{$sample} = $group;
		}
		my @groups = sort keys %samplesByGroup;
		my ($workerForGroup, $workerLoads) =
			balance_assembly_groups(\%samplesByGroup, $maxSubJob);
		my (%mine, %ownedGroup);
		for my $group (@groups) {
			next unless $workerForGroup->{$group} == $subJob;
			$ownedGroup{$group} = 1;
			$mine{$_} = 1 for @{$samplesByGroup{$group}};
		}
		# Assembly catalogues can use generated aliases such as sampleM2.
		for my $alias (keys %{$map{altNms} || {}}) {
			my $sample = $map{altNms}{$alias};
			my $group = $groupForSample{$sample};
			$mine{$alias} = 1 if defined($group) && $ownedGroup{$group};
		}
		$mySamplesHR = \%mine;
		my $totalWorkerLoad = 0;
		$totalWorkerLoad += $_ for @{$workerLoads};
		my $plannedSamples = 0;
		$plannedSamples += scalar(@{$samplesByGroup{$_}}) for keys %ownedGroup;
		print "Subjob ${subJob}/$maxSubJob: restricting locus-model construction to "
			. scalar(keys %ownedGroup)." of ".scalar(@groups)
			." assembly groups ($plannedSamples planned sample(s), "
			.scalar(keys %mine)." sample/alias identifiers; "
			."estimated sample load $workerLoads->[$subJob]/$totalWorkerLoad)\n";
	}

	my ($hr1,$cl2gene) = readClstrRev("$GCd/compl.incompl.$clusterID.fna.clstr.idx",0,$Gene2COG,$mySamplesHR);
	$hr1 = {};

	my $protein_file = "$GCd/compl.incompl.$clusterID.prot.faa";
	if (fileGZe($protein_file)) {
		$catalogProteins = readFasta($protein_file,1,"\\s",$Gene2COG, { fai => 1 });
	} else {
		warn "Catalogue protein file $protein_file is unavailable; keeping same-COG catalogue clusters separate\n";
	}

	my @records;
	for my $MGS (keys %{$COGprios}) {
		# In tree-recalculation mode, published or complete staged inputs are
		# reused. Build an extraction model only for MGS that have no such input.
		next if $recalcTrees && !$MGSneedsExtraction{$MGS};
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
			# Every member pair in a multi-seed Mosaic locus must be independently
			# confirmed.  This prevents an A-B-C chain from silently merging A and C.
			allowed_merge_pairs => \%ConfirmedMosaicPairs,
			require_complete_linkage => 1,
			allow_confirmed_cooccurrence => 1,
		},
	);
	my $ranked_record_count = scalar(@records);
	my $linkage_rejections = $locus_model->{incomplete_linkage_rejections} || 0;
	print "Mosaic complete-linkage protection rejected $linkage_rejections "
		."transitive component merge(s)\n" if $linkage_rejections;
	@records = ();
	my @selected_locus_groups = @{$locus_model->{groups}};
	my $locus_budget_excluded = 0;
	if (!$taxonAwareLocusSelection) {
		my %selected_loci_by_mgs;
		@selected_locus_groups = grep {
			if (($selected_loci_by_mgs{$_->{mgs}} // 0) >= $treeLocusBudget) {
				$locus_budget_excluded++;
				0;
			} else {
				$selected_loci_by_mgs{$_->{mgs}}++;
				1;
			}
		} @selected_locus_groups;
	}
	$LocusByID = {
		map { $_->{locus_id} => $_ } @selected_locus_groups
	};
	$MemberContext = $locus_model->{member_context};
	$LocusContext = $locus_model->{locus_context};

	my ($new_si_genes, $new_priorities) = ({}, {});
	for my $group (@selected_locus_groups) {
		$new_si_genes->{$group->{mgs}}{$group->{locus_id}} = $group->{primary_gene};
		push @{$new_priorities->{$group->{mgs}}}, $group->{locus_id};
	}
	$SIgenes = $new_si_genes;
	$COGprios = $new_priorities;

	my ($gene_sample_combinations, $ambiguous_seed_samples, $missing_clusters) = (0, 0, 0);
	my $unrepresentedWorkerLoci = 0;
	my (%contextMembersNeeded, %contextLociNeeded);
	for my $group (@selected_locus_groups) {
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
	print "Prepared ".scalar(@selected_locus_groups)." loci from $ranked_record_count"
		." ranked catalogue clusters; merged $locus_model->{merged_seeds} compatible same-COG seeds. "
		."$gene_sample_combinations locus-sample combinations, $ambiguous_seed_samples with multiple candidates"
		.($missing_clusters ? ", $missing_clusters missing cluster-index entries" : "")
		.($unrepresentedWorkerLoci ? ", $unrepresentedWorkerLoci loci outside this worker's sample slice" : "")
		.($locus_budget_excluded ? ", $locus_budget_excluded lower-ranked loci excluded by -treeLocusBudget" : "")
		.".\n";
}

sub resolveScratchDirectory {
	my ($derived, $manifest, $catalogIdentity, $canonicalOutD, $outDname) = @_;
	return $derived unless -s $manifest;

	open my $input, '<', $manifest or do {
		limitedWarn('scratch manifest',
			"Cannot read scratch-directory manifest $manifest: $!; using $derived\n");
		return $derived;
	};
	my $header = <$input> // '';
	my $row = <$input> // '';
	my $extra = '';
	while (my $line = <$input>) {
		$extra .= $line if $line =~ /\S/;
	}
	close $input;
	$header =~ s/[\r\n]+\z//;
	$row =~ s/[\r\n]+\z//;
	my @fields = split /\t/, $row, -1;
	my $valid = $header eq join("\t",
		qw(version catalog_identity output_path scratch_directory))
		&& @fields == 4 && $fields[0] eq '1'
		&& $fields[1] eq $catalogIdentity
		&& $fields[2] eq $canonicalOutD
		&& !$extra && File::Spec->file_name_is_absolute($fields[3])
		&& $fields[3] !~ /[\t\r\n]/;
	if ($valid) {
		my $normalized = File::Spec->canonpath($fields[3]);
		my @parts = File::Spec->splitdir($normalized);
		pop @parts while @parts && $parts[-1] eq '';
		my $leaf = pop(@parts) // '';
		my $catalogPart = pop(@parts) // '';
		my $namespace = pop(@parts) // '';
		$valid = $namespace eq 'strainsScr1'
			&& $catalogPart eq $catalogIdentity
			&& $leaf =~ /^\Q$outDname\E\.[0-9a-f]{16}\z/;
		if ($valid) {
			$normalized .= '/' unless $normalized =~ m{/\z};
			print "Reusing recorded scratch directory: $normalized\n";
			return $normalized;
		}
	}
	limitedWarn('scratch manifest',
		"Invalid or stale scratch-directory manifest $manifest; using $derived\n");
	return $derived;
}

sub persistScratchDirectory {
	my ($manifest, $catalogIdentity, $canonicalOutD, $directory) = @_;
	die "Unsafe scratch-directory manifest fields\n"
		if grep { !defined($_) || /[\t\r\n]/ }
			($catalogIdentity, $canonicalOutD, $directory);
	my $contents = join("\t",
		qw(version catalog_identity output_path scratch_directory))."\n"
		.join("\t", 1, $catalogIdentity, $canonicalOutD, $directory)."\n";
	if (-s $manifest && open(my $existing, '<', $manifest)) {
		local $/;
		my $current = <$existing> // '';
		close $existing;
		return if $current eq $contents;
	}
	my $temporary = "$manifest.tmp.$$";
	retry_unlink($temporary, fatal => 0, label => 'clean scratch manifest temporary');
	my $output = retry_open('>', $temporary,
		label => 'create scratch-directory manifest');
	print {$output} $contents
		or die "Cannot write scratch-directory manifest $temporary: $!\n";
	retry_close($output, 'close scratch-directory manifest');
	retry_rename($temporary, $manifest,
		label => 'publish scratch-directory manifest');
}

sub migrateLegacyOperationalLogs {
	return 0 unless -d $outD && -d $LOGDIR;
	opendir my $directory, $outD
		or die "Cannot inspect strain output directory $outD: $!\n";
	my @legacy = grep {
		/\A(?:
			strain_within(?:\.worker\.\d+)?\.(?:heartbeat|failure)\.tsv
			|SNPconsCalls(?:\.\d+)?\.log
			|strainSampleStats(?:\.summary)?\.tsv(?:\.\d+)?
			|strainRecovery\.tsv(?:\.\d+|\.contributors\.tsv)?
			|strain_within\.summary\.log
			|strainSelectionAttrition\.tsv
			|strainCmd\.txt
			|strainAnalysis2\.sh(?:\.(?:o|e)txt)?
			|ConspecificMGS(?:\.\d+)?\.log
		)\z/x
	} readdir $directory;
	closedir $directory
		or die "Cannot close strain output directory $outD: $!\n";
	my $migrated = 0;
	for my $name (sort @legacy) {
		my $source = File::Spec->catfile($outD, $name);
		my $destination = File::Spec->catfile($LOGDIR, $name);
		next unless -f $source;
		if (-e $destination) {
			my $suffix = 1;
			my $legacyDestination;
			do {
				$legacyDestination = "$destination.legacy.$suffix";
				$suffix++;
			} while (-e $legacyDestination);
			limitedNotice('legacy operational log migration collision',
				"Preserving legacy $source as $legacyDestination\n");
			$destination = $legacyDestination;
		}
		retry_rename($source, $destination,
			label => "migrate legacy strain log $name");
		$migrated++;
	}
	print "Migrated $migrated legacy strain workflow file(s) into $LOGDIR\n"
		if $migrated;
	return $migrated;
}

sub prepRun{

	$mode = "FMG" if ($MGSfile eq "");
	$maxSubJob = 0 if $mode eq "FMG" && $maxSubJob == -1;
	die "FMG mode does not support positive -maxSubJob; run it as a single extraction job\n"
		if $mode eq "FMG" && $maxSubJob > 0;


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
	$SNPconsLOGs = "$LOGDIR/SNPconsCalls.$subJob.log" if ($SNPconsLOGs eq "");

	my $GCname = basename($GCd);
	my $outDname = basename($outD);
	die "Could not derive safe temporary-directory names\n" unless length($GCname) && length($outDname);
	my $catalogIdentity = catalog_identity($GCd);
	my $canonicalOutD = abs_path($outD) || File::Spec->rel2abs($outD);
	my $manifestOutputPath = File::Spec->canonpath(File::Spec->rel2abs($outD));
	my $outputIdentity = substr(sha256_hex($canonicalOutD), 0, 16);
	my $scratchRoot = getProgPaths("globalTmpDir",0);
	$scratchRoot = "$outD/.scratch" if $scratchRoot eq "";
	my $derivedScratch = "$scratchRoot/strainsScr1/$catalogIdentity/$outDname.$outputIdentity/";
	$derivedScratch = File::Spec->canonpath(File::Spec->rel2abs($derivedScratch));
	$derivedScratch .= '/' unless $derivedScratch =~ m{/\z};
	my $scratchManifest = "$outD/.strain_within.scratch.tsv";
	$scratchD = resolveScratchDirectory($derivedScratch, $scratchManifest,
		$catalogIdentity, $manifestOutputPath, $outDname);
	#die "$scratchD  :$GCname :$GCd\n";
	if ($locTmpDir1 eq ""){
		my $locTmpN = getProgPaths("nodeTmpDir",0) ;
		my $suffix = ""; $suffix = "/SJ.${subJob}/" if ($subJob);
		if ($locTmpN eq ""){
			$locTmpDir = "$outD/strainsScr1/$catalogIdentity/$outDname.$outputIdentity/$suffix";
		} else {
			#my $tmp = `echo \$SLURM_LOCAL_SCRATCH`;
			#print STDERR "echo $locTmpN\n$tmp\n";
			#$locTmpN =~ s/\$/\\\$/;$locTmpN = `echo $locTmpN;`; #eval in sys
			$locTmpN=truePath($locTmpN,1);
			#die $locTmpN."\n";
			$locTmpDir = "$locTmpN/strainsScr1/$catalogIdentity/$outDname.$outputIdentity/$suffix";
		}
	} else {
		my $suffix = $subJob ? "/SJ.${subJob}" : "";
		$locTmpDir = "$locTmpDir1/strainsScr1/$catalogIdentity/$outDname.$outputIdentity$suffix/";
	}
	
	$preConDir = "$scratchD/preComp/";


	print "\n!! WARNING !!: RESUBMISSION mode selected (will resubmit MSA + phylos even for already completed MGS) !!\n" if ($reSubmit);
	print "\n!! WARNING !!: RECALCTREES mode selected (will delete per-MGS tree outputs and rebuild from published or complete staged inputs) !!\n" if ($recalcTrees);
	print "\n!! WARNING !!: REDOSUBMISSIONDATA mode selected (will redo and resubmit MSA + phylos even for already completed MGS) !!\n" if ($redoSubmissionData);

	$mapF = resolve_catalog_maps($GCd);
	
	#read info gene <-> taxonomy from this file, depends on config..
	$gene2taxF = "$GCd/FMG/gene2specI.txt";
	$gene2taxF = "$GCd/GTDBmg/gene2specI.txt" if ($useGTDBmg eq "GTDB");
	#die;

	#---------------
	#everything after is only for main submission job..
	if ($subJob){
		print "----- Resolved run configuration: worker ${subJob}/$maxSubJob -----\n";
	} else {
		print "----- Resolved run configuration -----\n";
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
		print "MSAaligner: $MSAprog, backbone GenesPerSpecies: $GenesPerSpecies, "
			."GeneLengthMin: $GeneLengthMin, backbone relativeNTFraction: $relativeNTFraction, "
			."placement GenesPerSpecies: $placementGenesPerSpecies, "
			."placement relativeNTFraction: $placementRelativeNTFraction, "
			."taxonAwareLocusSelection: $taxonAwareLocusSelection\n";
		print "Rate/GC partition merging: enabled=$rateMergePartitions, "
			."maximumBins=$rateMergeMaxBins, targetBin=$rateMergeTargetSites effective sites, "
			."minimumBin=$rateMergeMinLoci loci/"
			."$rateMergeMinSites sites\n";
		print "Taxon-aware locus hierarchy: extractionPool="
			.($noGeneLimit ? 'unlimited' : $maxNGenes)
			.", finalTreeBudget=$taxonAwareGeneBudget, "
			."robustCore=$taxonAwareCoreLoci, taxonRescue="
			.($taxonAwareMaxLoci - $taxonAwareCoreLoci)
			.", qcBackfill=$taxonAwareCandidateExtra, rescueMinPrevalence="
			."$taxonAwareRescueMinPrevalence, preferredCore="
			.(length($preferredCoreGenes) ? $preferredCoreGenes : '<none>')
			.", compactDiagnostics=$compactTaxonAwareDiagnostics\n"
			if $taxonAwareLocusSelection;
		
		
		if ($noGeneLimit){print "No per-sample gene limit; biological QC remains "
			.($disableQC ? "disabled by explicit request\n" : "enabled\n");}
		else {print "Using at most $maxNGenes genes per sample after QC\n";}
		print "Filtering defaults: MGSminGenesPSmpl=$MGStoolowGsThr, "
			."multiGeneSmplMax=$multiGeneSmplMax, conspGeneSmplMax=$conspGeneSmplMax, "
			."breakpointGeneFlank=$breakpointGeneFlank, abundanceMinLoci=$abundanceMinimumLoci, "
			."abundanceFold=$abundanceMinimumFold-$abundanceMaximumFold, "
			."abundanceMaxModifiedZ=$abundanceMaximumModifiedZ\n";
		print "==============================================\n";
		if ($onlySubmit){print "Only submission mode\n";
		} elsif (!$subJob) {
			print "Creation of strain genes, old data might be deleted!\nDo you want to continue? (10s wait, use Ctrl-c to abort)\n"; sleep 10;
		}
	}
	
	

	#$mapF = $GCd."LOGandSUB/inmap.txt" if ($mapF eq "");
	my ($hr1,$hr2) = readMapS($mapF,-1);
	%map = %{$hr1}; %AsGrps = %{$hr2};
	#get every sample: assembly references can be shared, but VCFs, depths, and
	#consensus sequences remain sample-specific
	@samples = @{$map{opt}{smpl_order}};
	my %sample_seen;
	for my $sample (@samples) {
		die "Unsafe sample identifier '$sample': use only letters, digits, dot, underscore, colon, plus, and hyphen\n"
			unless defined($sample) && $sample =~ /\A[A-Za-z0-9][A-Za-z0-9_.:+-]*\z/;
		die "Duplicate sample identifier in map: $sample\n" if $sample_seen{$sample}++;
	}


	if ($mode eq "MGS" || $mode eq "MGSall"){
		my $sortedMGS = "$MGSfile.srt";
		if ($preparedMainBranchFastPath) {
			$MGSfile = -s $sortedMGS ? $sortedMGS : $MGSfile;
			$gene2taxF = "";
			print "Prepared-input recovery: skipped MGS sorting and gene-to-MGS index creation.\n";
		} else {
		if ($subJob) {
			die "Sorted MGS guide is missing for subjob: $sortedMGS\n" unless -s $sortedMGS;
		} elsif ($recalcTrees && !-s $sortedMGS) {
			die "-recalcTrees requires the existing sorted MGS guide: $sortedMGS\n";
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
		}
	} elsif ($subJob && $maxSubJob) {
		die "FMG mode does not support split MGS extraction jobs\n";
	}

	if ($subJob){
		return;
	}

	make_path($locTmpDir, $scratchD, $outD, $LOGDIR);
	migrateLegacyOperationalLogs();
	my $outputBase = basename(File::Spec->canonpath($outD));
	my $owner = File::Spec->catfile($outD, '.matafiler-strain-workdir');
	persistScratchDirectory($scratchManifest, $catalogIdentity,
		$manifestOutputPath, $scratchD);
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
	my $submPreComp = 1;

	
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
	%AGlist = ();
	my %seen;
	foreach my $smpl (@samples){
		die "Sample $smpl has no assembly-group assignment\n"
			unless exists($map{$smpl}) && defined($map{$smpl}{AssGroup});
		my $cAssGrp = $map{$smpl}{AssGroup};
		next if $cAssGrp eq "-1"; #the caller handles a standalone sample directly
		next if $seen{$cAssGrp}{$smpl}++;
		# Mapping-group members may share the reference assembly, but each
		# sample has its own VCF/depth and must contribute a consensus.
		push @{$AGlist{$cAssGrp}}, $smpl;
	}
}

sub histoMGS{#specifically for MGS..
	my ($aref,$msg) = @_;
	my @cnts = @{$aref};
	my @binSiz = (10,20,30,50,70,100,200,300,700,1000,2000,5000,10000,1e6);
	push @binSiz, $maxNGenes if !$noGeneLimit && $maxNGenes > 0;
	my %seenBin;
	@binSiz = sort { $a <=> $b } grep { !$seenBin{$_}++ } @binSiz;
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
	# Emit the complete diagnostic in one write.  A bare newline write can
	# appear as an empty stdout record in scheduler stream collectors.
	my @displayBins = map { " <=$_:$binC{$_} " }
		grep { $binC{$_} > 0 } @binSiz;
	print $msg.": ".join("", @displayBins)."\n";
	#DEBUG
	#print @cnts." : @cnts\n";
}

sub outgroupRequirementLoci {
	my ($MGS) = @_;
	my $published = $SIdirs{$MGS} // "$outD/$MGS/";
	my $scratch = "$scratchD/outs/$MGS";
	my $stagedReady = scratchMGSInputState($MGS) eq 'complete';
	my $reusePrepared = !$repairCAT && !$deepRepair && !$redoSubmissionData
		&& !exists($legacyLocusMGS{$MGS});
	if ($reusePrepared && persistentMGSInputState($MGS) eq 'complete') {
		my ($prepared) = preparedOutgroupLog($published);
		return ([], 'published_overlay', 0) if $prepared;
	}
	if ($reusePrepared && $stagedReady) {
		my ($prepared) = preparedOutgroupLog($scratch);
		return ([], 'staged_overlay', 0) if $prepared;
	}

	return ([], 'no_complete_staging', 0) unless $stagedReady;
	if (exists($COGprios->{$MGS}) && @{$COGprios->{$MGS}}) {
		my %seen;
		my @selected = grep { defined($_) && length($_) && !$seen{$_}++ }
			@{$COGprios->{$MGS}};
		return (\@selected, 'selected_gene_map', undef) if @selected;
	}
	my @sources;
	my $aggregate = "$scratch/$CATstdof.tmp";
	if (fileGZe($aggregate)) {
		@sources = ($aggregate);
	} else {
		@sources = exact_worker_parts($aggregate, $maxSubJob || 1);
	}
	return ([], 'no_raw_category', 0) unless @sources;
	my (%seen, %samples, @loci);
	my $rows = 0;
	my $scanStarted = time;
	my $nextProgress = time + 60;
	for my $source (@sources) {
		my ($input) = gzipopen($source, "outgroup requirement category", 1);
		while (my $line = <$input>) {
			$rows++;
			if (time >= $nextProgress) {
				stepProgress("outgroup requirement category scan for $MGS",
					$rows, undef, $scanStarted, "loci=".scalar(keys %seen),
					"samples=".scalar(keys %samples));
				$nextProgress = time + 60;
			}
			$line =~ s/[\r\n]+\z//;
			next unless length($line);
			my @field = split /\t/, $line, -1;
			die "Malformed raw staged category row for $MGS: $line\n"
				unless @field >= 4 && $field[0] eq $MGS
					&& length($field[1]) && length($field[2]) && length($field[3]);
			$samples{$field[2]} = 1 if defined($field[2]) && length($field[2]);
			my ($locusMGS, $cog, $primaryGene) = locusParts($field[1], $MGS);
			next unless $locusMGS eq $MGS && length($cog) && length($primaryGene);
			my $locus = $field[1];
			push @loci, $locus unless $seen{$locus}++;
		}
		close $input or die "Cannot close outgroup requirement category $source: $!\n";
	}
	return (\@loci, @sources > 1 ? 'worker_categories' : 'aggregate_category',
		scalar(keys %samples));
}

sub writeSortedIdentifierSet {
	my ($path, $identifiers) = @_;
	my $temporary = "$path.tmp.$$";
	my $output = retry_open('>', $temporary,
		label => "create outgroup-reference identifier list");
	for my $identifier (sort keys %{$identifiers}) {
		die "Unsafe outgroup-reference identifier '$identifier'\n"
			unless defined($identifier) && length($identifier)
				&& $identifier !~ /[\s\x00-\x1f\x7f]/;
		print {$output} "$identifier\n"
			or die "Cannot write outgroup-reference identifier list $temporary: $!\n";
	}
	retry_close($output, 'close outgroup-reference identifier list');
	retry_rename($temporary, $path,
		label => "publish outgroup-reference identifier list");
}

sub prepareSelectiveOutgroupReferenceCache {
	my ($sourceFNA, $sourceFAA, $requiredNT, $requiredAA) = @_;
	return ('', '') if $mode eq 'FMG' || $sourceFNA =~ /[*?\[]/
		|| $sourceFAA =~ /[*?\[]/;
	$outgroupReferenceCacheDir = "$scratchD/outgroup_reference_cache";
	make_path($outgroupReferenceCacheDir) unless -d $outgroupReferenceCacheDir;
	my $ntIDs = "$outgroupReferenceCacheDir/required.nt.ids";
	my $aaIDs = "$outgroupReferenceCacheDir/required.aa.ids";
	writeSortedIdentifierSet($ntIDs, $requiredNT);
	writeSortedIdentifierSet($aaIDs, $requiredAA);
	my $helper = getProgPaths('MGS_outgroup_ref_scr');
	my $command = join(' ', $helper,
		'-nt', shellQuote($sourceFNA), '-aa', shellQuote($sourceFAA),
		'-ntIDs', shellQuote($ntIDs), '-aaIDs', shellQuote($aaIDs),
		'-outD', shellQuote($outgroupReferenceCacheDir),
	);
	$command .= ' -mosaic '.shellQuote($mosaicLociFile)
		if length($mosaicLociFile) && -s $mosaicLociFile;
	print "Preparing indexed selective outgroup-reference cache: NT="
		.scalar(keys %{$requiredNT}).", AA=".scalar(keys %{$requiredAA})
		."; cache=$outgroupReferenceCacheDir; global elapsed "
		.timeNice(time - $^T)."\n";
	my $status = system('bash', '-o', 'pipefail', '-c', $command);
	if ($status != 0) {
		my $detail = $status == -1 ? $! : 'exit '.($status >> 8);
		limitedWarn('selective outgroup-reference cache failures',
			"Selective outgroup-reference helper failed ($detail); retaining its cache workspace for diagnosis\n");
		return ('', '');
	}
	my $cachedFNA = "$outgroupReferenceCacheDir/references.fna";
	my $cachedFAA = "$outgroupReferenceCacheDir/references.faa";
	unless (-s $cachedFNA && -s "$cachedFNA.fai"
			&& -s $cachedFAA && -s "$cachedFAA.fai"
			&& -s "$outgroupReferenceCacheDir/complete.sto") {
		limitedWarn('incomplete selective outgroup-reference cache',
			"Selective outgroup-reference helper returned success without a complete indexed cache\n");
		return ('', '');
	}
	return ($cachedFNA, $cachedFAA);
}

sub preparedMainBranchInputSet {
	my ($guide, $outputDirectory, $subset) = @_;
	return (0, [], "output directory is absent")
		unless defined($outputDirectory) && -d $outputDirectory;
	return (0, [], "MGS guide is unspecified")
		unless defined($guide) && length($guide);
	my $sortedGuide = $guide =~ /\.srt\z/ ? $guide : "$guide.srt";
	$guide = $sortedGuide if -s $sortedGuide;
	return (0, [], "MGS guide is absent") unless -s $guide;
	open my $input, "<", $guide
		or return (0, [], "MGS guide cannot be opened: $!");
	my %requested = map { $_ => 1 } @{$subset || []};
	my (%seen, @mgs);
	my $lineNumber = 0;
	my $nextProgress = time + 60;
	while (my $line = <$input>) {
		$lineNumber++;
		next if $line =~ /^\s*\z/;
		my $tab = index($line, "\t");
		unless ($tab > 0) {
			close $input;
			return (0, [], "MGS guide row $lineNumber is not tab-delimited");
		}
		my $mgs = substr($line, 0, $tab);
		next if %requested && !$requested{$mgs};
		unless ($mgs =~ /\A[A-Za-z0-9][A-Za-z0-9_.:+-]*\z/) {
			close $input;
			return (0, [], "unsafe MGS identifier on guide row $lineNumber");
		}
		push @mgs, $mgs unless $seen{$mgs}++;
		if (time >= $nextProgress) {
			print STDERR "Prepared-input preflight: scanned $lineNumber guide rows; "
				.scalar(@mgs)." selected MGS; global elapsed "
				.timeNice(time - $^T)."\n";
			$nextProgress = time + 60;
		}
	}
	close $input
		or return (0, [], "MGS guide cannot be closed: $!");
	if (%requested) {
		my @missing = grep { !$seen{$_} } sort keys %requested;
		return (0, [], "selected MGS absent from guide: ".join(",", @missing))
			if @missing;
	}
	return (0, [], "MGS guide selected no MGS") unless @mgs;
	my $checkedMGS = 0;
	for my $mgs (@mgs) {
		$checkedMGS++;
		if (time >= $nextProgress) {
			print STDERR "Prepared-input preflight: checked $checkedMGS/".scalar(@mgs)
				." MGS input sets; current_MGS=$mgs; global elapsed "
				.timeNice(time - $^T)."\n";
			$nextProgress = time + 60;
		}
		my $directory = File::Spec->catdir($outputDirectory, $mgs);
		my $terminal = grep { -s File::Spec->catfile($directory, $_) }
			qw(tooFewSamples.sto noRecoverableLoci.sto noTree.sto);
		next if $terminal;
		my @missing = grep {
			!fileGZe(File::Spec->catfile($directory, $_))
		} ($FNAstdof, $FAAstdof, $CATstdof);
		return (0, [], "$mgs lacks final ".join("/", @missing)) if @missing;
		my ($category) = gzipopen(
			File::Spec->catfile($directory, $CATstdof),
			"prepared-input category preflight", 1,
		);
		return (0, [], "$mgs category cannot be opened") unless $category;
		my $firstIdentifier = "";
		while (my $line = <$category>) {
			next unless $line =~ /\S/;
			($firstIdentifier) = split /\t/, $line, 2;
			last;
		}
		close $category
			or return (0, [], "$mgs category cannot be closed");
		$firstIdentifier =~ s/[\r\n]+\z//;
		my @identifierParts = split /\Q|\E/, $firstIdentifier, -1;
		return (0, [], "$mgs has a legacy or empty category file")
			unless @identifierParts == 3 && !grep { !length($_) } @identifierParts;
		$preparedMainBranchCategoryValidated{$mgs} = 1;
		my ($outgroupPrepared) = preparedOutgroupLog($directory);
		return (0, [], "$mgs lacks a valid final outgroup record")
			unless $outgroupPrepared;
	}
	return (1, \@mgs, "ready");
}

sub persistentMGSInputState {
	my ($MGS, $refresh) = @_;
	return $persistentMGSInputStateCache{$MGS}
		if !$refresh && exists $persistentMGSInputStateCache{$MGS};
	my $mgsDir = $SIdirs{$MGS} // "$outD/$MGS";
	my @required = ($FNAstdof, $FAAstdof, $CATstdof);
	my @present = grep { fileGZe("$mgsDir/$_") } @required;
	my $state = @present == @required ? 'complete'
		: @present ? 'incomplete' : 'missing';
	return $persistentMGSInputStateCache{$MGS} = $state;
}

sub scratchMGSInputState {
	my ($MGS, $refresh) = @_;
	return $scratchMGSInputStateCache{$MGS}
		if !$refresh && exists $scratchMGSInputStateCache{$MGS};
	if (stagedMGSInputsReady($MGS)) {
		return $scratchMGSInputStateCache{$MGS} = 'complete';
	}
	my $mgsDir = "$scratchD/outs/$MGS";
	my @required = ($FNAstdof, $FAAstdof, $LINKstdof, "$CATstdof.tmp", "$QCstdof.tmp", $CATstdof, $QCstdof);
	my $hasAny = grep {
		fileGZe("$mgsDir/$_") || bsd_glob("$mgsDir/$_.*")
	} @required;
	return $scratchMGSInputStateCache{$MGS} =
		$hasAny ? 'incomplete' : 'missing';
}

sub invalidateMGSInputState {
	my (@MGS) = @_;
	delete @persistentMGSInputStateCache{@MGS};
	delete @scratchMGSInputStateCache{@MGS};
}

sub getInputSize{
	my ($targets) = @_;
	my @targetMGS = defined($targets) ? @{$targets} : @specis;
	my (@out, %counts);
	my $sizingStarted = time;
	my $sizedMGS = 0;
	my $nextSizingProgress = time + 60;
	my $auditFile = "$LOGDIR/tree_input_sizing.tsv";
	my $auditTemporary = "$auditFile.tmp.$$";
	make_path($LOGDIR) unless -d $LOGDIR;
	retry_unlink($auditTemporary, fatal => 0, label => "clean stale tree-input audit");
	open my $audit, '>', $auditTemporary
		or die "Cannot create tree-input audit $auditTemporary: $!\n";
	print {$audit} join("\t", qw(
		MGS selected_state source estimated_uncompressed_MB persistent_state scratch_state
	)), "\n";
	foreach my $MGS (@targetMGS){
		my $tmpD = "$scratchD/outs/$MGS";
		my $publishedState = persistentMGSInputState($MGS);
		# A complete published triplet is authoritative for ordinary tree recovery.
		# Avoid worker-part globs and merge-checkpoint probes in scratch for it.
		my $scratchState = $publishedState eq 'complete'
			? 'not_checked' : scratchMGSInputState($MGS);
		my $tooFew = -s "$SIdirs{$MGS}/tooFewSamples.sto" ? 1 : 0;
		my $noRecoverableLoci = -s "$SIdirs{$MGS}/noRecoverableLoci.sto" ? 1 : 0;
		my ($state, $source, $inputFNAsize) = ('empty_extraction', 'none', 0);
		if ($tooFew) {
			($state, $source) = ('too_few_samples', 'marker');
			$counts{too_few_samples}++;
		} elsif ($noRecoverableLoci) {
			($state, $source) = ('no_recoverable_loci', 'marker');
			$counts{no_recoverable_loci}++;
		} elsif ($scratchState eq 'complete') {
			($state, $source) = ('ready', 'scratch');
			$counts{ready}++;
			my @sizeSources;
			if (-e "$tmpD/$FNAstdof") {
				@sizeSources = ("$tmpD/$FNAstdof");
			} elsif (-e "$tmpD/$FNAstdof.gz") {
				@sizeSources = ("$tmpD/$FNAstdof.gz");
			} else {
				@sizeSources = exact_worker_parts(
					"$tmpD/$FNAstdof", $maxSubJob || 1);
			}
			for my $path (@sizeSources) {
				$inputFNAsize += fileGZs($path) / (1024 * 1024);
			}
		} elsif ($publishedState eq 'complete') {
			($state, $source) = ('ready', 'published');
			$counts{ready}++;
			my $fna = "$SIdirs{$MGS}/$FNAstdof";
			if (-e $fna) {
				$inputFNAsize = fileGZs($fna) / (1024 * 1024);
			} elsif (-e "$fna.gz") {
				$inputFNAsize = fileGZs("$fna.gz") / (1024 * 1024) * 40;
			}
		} elsif ($scratchState eq 'incomplete' || $publishedState eq 'incomplete') {
			$state = $scratchState eq 'incomplete' && $publishedState eq 'incomplete'
				? 'incomplete_scratch_and_published'
				: $scratchState eq 'incomplete' ? 'incomplete_scratch' : 'incomplete_published';
			$source = 'none';
			$counts{incomplete_scratch}++ if $scratchState eq 'incomplete';
			$counts{incomplete_published}++ if $publishedState eq 'incomplete';
		} else {
			$counts{empty_extraction}++;
		}
		print {$audit} join("\t", $MGS, $state, $source,
			sprintf('%.6f', $inputFNAsize), $publishedState, $scratchState), "\n";
		push @out, $inputFNAsize;
		$sizedMGS++;
		if (time >= $nextSizingProgress) {
			stepProgress("full-tree input sizing", $sizedMGS, scalar(@targetMGS),
				$sizingStarted, "current_MGS=$MGS", "ready=".($counts{ready} // 0),
				"incomplete=".(($counts{incomplete_scratch} // 0)
					+ ($counts{incomplete_published} // 0)));
			$nextSizingProgress = time + 60;
		}
	}
	close $audit or die "Cannot close tree-input audit $auditTemporary: $!\n";
	rename $auditTemporary, $auditFile
		or die "Cannot publish tree-input audit $auditTemporary as $auditFile: $!\n";
	return (\@out, {
		too_few_samples => $counts{too_few_samples} // 0,
		no_recoverable_loci => $counts{no_recoverable_loci} // 0,
		incomplete_published => $counts{incomplete_published} // 0,
		incomplete_scratch => $counts{incomplete_scratch} // 0,
		empty_extraction => $counts{empty_extraction} // 0,
		audit_file => $auditFile,
	});
}


sub stagedMGSInputsReady {
	my ($MGS) = @_;
	my $mgsDir = "$scratchD/outs/$MGS";
	my @coreRequiredNames = ($FNAstdof, $FAAstdof, $LINKstdof);
	my $mergeCheckpoint = "$mgsDir/merge.complete.tsv";
	my $aggregateComplete = !grep { !fileGZe("$mgsDir/$_") } @coreRequiredNames;
	$aggregateComplete &&= (
		(fileGZe("$mgsDir/$CATstdof.tmp") && fileGZe("$mgsDir/$QCstdof.tmp"))
		|| (fileGZe("$mgsDir/$CATstdof") && fileGZe("$mgsDir/$QCstdof"))
	);
	$aggregateComplete &&= -s $mergeCheckpoint;
	# The checkpoint commits the complete aggregate. Fresh worker parts are
	# considered by combineMGSgenesDir when a merge is requested; readiness and
	# sizing do not need five directory globs for an already reusable aggregate.
	return 1 if $aggregateComplete;
	my $workerCount = $maxSubJob || 1;
	my @requiredNames = (
		$FNAstdof, $FAAstdof, $LINKstdof, "$CATstdof.tmp", "$QCstdof.tmp",
	);
	my %parts = map {
		$_ => [exact_worker_parts("$mgsDir/$_", $workerCount)]
	} @requiredNames;
	my $hasFreshParts = grep { @{$parts{$_}} } @requiredNames;
	return $aggregateComplete unless $hasFreshParts;
	if ($maxSubJob && !split_generation_complete(
			"$LOGDIR/mainExtr.generation", "$LOGDIR/mainExtr", $maxSubJob)) {
		return $aggregateComplete;
	}
	return 0 if grep { !@{$parts{$_}} } @requiredNames;
	return 1;
}

sub evalFileStatus{
	my $dirsNOTPrepped = 0; my $CatFileMiss = 0;my $CatNotPrepped = 0; my $treeAbsent=0;
	my $doneDirs=0;
	my $completedTreeFastPaths=0;
	my $tooFewDirs=0;
	my $noRecoverableLociDirs=0;
	my $PhylosExist = 1;
	my $auditStarted = time;
	my $auditedMGS = 0;
	my $nextAuditProgress = time + 60;
	
	my $treeFile= "IQtree_allsites.treefile";
	if ($phyloProg == 2){$treeFile = "VERYFASTTREE_allsites.nwk";} elsif ($phyloProg == 3){$treeFile = "FASTTREE_allsites.nwk";}


	foreach my $MGS (@specis){ #loop creates per specI file structure to run buildTreeScript on..
		#PART I: create fasta files required by tree
		my $outD2 = "$outD/$MGS/";
		$SIdirs{$MGS} = $outD2;
		my $completedTree = "$outD2/phylo/$treeFile";
		my $treeCompletion = "$outD2/treeDone.sto";
		if (!$recalcTrees && !$reSubmit && !$repairCAT && !$deepRepair
				&& !$redoSubmissionData && ($onlySubmit != 0 || $subJob)
				&& -s $treeCompletion && fileGZs($completedTree)) {
			# BuildTree publishes treeDone.sto atomically only after validating the
			# primary tree and clearing terminal lifecycle markers.  On a tree-only
			# resume this pair is authoritative, so avoid per-MGS directory creation,
			# terminal-marker checks, and input-sidecar probes.
			$doneDirs++;
			$completedTreeFastPaths++;
			$MGSsubmissionComplete{$MGS} = 1;
			$deferredScratchCleanup{"$scratchD/outs/$MGS"} = 1
				if -d "$scratchD/outs/$MGS";
			next;
		}
		#print "$outD2\n";
		if (-d $outD2 && $onlySubmit == 0 && !$subJob && !$recalcTrees){#only the parent may clean shared folders
			remove_tree($outD2);
			my $scratch_mgs = "$scratchD/outs/$MGS";
			remove_tree($scratch_mgs) if -d $scratch_mgs;
		}
		make_path($outD2) unless -d $outD2;
		my $tooFewMarker = "$outD2/tooFewSamples.sto";
		my $noRecoverableLociMarker = "$outD2/noRecoverableLoci.sto";
		my $buildTreeTerminalMarker = "$outD2/noTree.sto";
		if (-s $tooFewMarker && !$deepRepair && !$redoSubmissionData && ($onlySubmit != 0 || $recalcTrees)) {
			$MGSnoTree{$MGS} = 1;
			$MGSnoTreeReason{$MGS} = 'insufficient_tree_input';
			$tooFewDirs++;
			next;
		}
		if (-s $noRecoverableLociMarker && !$deepRepair && !$redoSubmissionData
				&& ($onlySubmit != 0 || $recalcTrees)) {
			$MGSnoTree{$MGS} = 1;
			$MGSnoTreeReason{$MGS} = 'no_recoverable_loci';
			$noRecoverableLociDirs++;
			next;
		}
		if (-s $buildTreeTerminalMarker && !$deepRepair && !$redoSubmissionData
				&& ($onlySubmit != 0 || $recalcTrees)) {
			$MGSnoTree{$MGS} = 1;
			$MGSnoTreeReason{$MGS} = lifecycleMarkerReason($buildTreeTerminalMarker,
				'buildtree_no_usable_alignment');
			$noRecoverableLociDirs++;
			next;
		}
		retry_unlink($tooFewMarker, label => "clear stale too-few marker");
		retry_unlink($noRecoverableLociMarker, label => "clear stale no-locus marker");
		retry_unlink($buildTreeTerminalMarker, label => "clear stale BuildTree terminal marker");
		
		# Placement recovery is independent of FNA/FAA/category publication.
		# Recognize it before probing a compressed category sidecar.
		if (my $epaRetryState = epaOnlyRetryReady($SIdirs{$MGS})) {
			$MGSepaOnlyRetry{$MGS} = $epaRetryState;
			$treeAbsent++;
			$deferredScratchCleanup{"$scratchD/outs/$MGS"} = 1
				if -d "$scratchD/outs/$MGS";
			next;
		}

	#	if ( !-d $outD2 ||){ # first phase only has "all.cat.tmp" file..
	#		$dirsNOTPrepped ++;
	#	} els
		
		if (fileGZe("$SIdirs{$MGS}/$CATstdof")) {
			my $first_entry = '';
			my $category_fh;
			if ($preparedMainBranchFastPath
					&& $preparedMainBranchCategoryValidated{$MGS}) {
				$first_entry = "prepared|category|$MGS";
			} else {
				($category_fh) = gzipopen("$SIdirs{$MGS}/$CATstdof",
					"existing locus category file", 0);
			}
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
				$MGSneedsExtraction{$MGS} = 1;
				$dirsNOTPrepped++;
				$CatFileMiss++;
				next;
			}
		}
		my $publishedInputState = persistentMGSInputState($MGS);
		if ($publishedInputState ne 'complete') {
			$CatFileMiss ++ ;
			if (stagedMGSInputsReady($MGS)) {
				# Preserve a complete Stage-I staging set. A later run can combine
				# and publish it in Part II, even when its published copy is partial.
				$CatNotPrepped ++;
			} else {
				$MGSneedsExtraction{$MGS} = 1;
				$dirsNOTPrepped ++; #missing or incomplete triplet requires Part I
			}
			#print "$SIdirs{$MGS}\n";
			#system "rm $SIdirs{$MGS}\n";
		} elsif(!fileGZs("$SIdirs{$MGS}/phylo/$treeFile")){
			$treeAbsent++;
			$deferredScratchCleanup{"$scratchD/outs/$MGS"} = 1
				if -d "$scratchD/outs/$MGS";
		} elsif(fileGZe("$SIdirs{$MGS}/phylo/$treeFile")) {
			$doneDirs++;
			$deferredScratchCleanup{"$scratchD/outs/$MGS"} = 1
				if -d "$scratchD/outs/$MGS";
		}
	} continue {
		$auditedMGS++;
		if (time >= $nextAuditProgress) {
			stepProgress("existing-output and resume audit", $auditedMGS,
				scalar(@specis), $auditStarted,
				"complete=$doneDirs", "EPA_only=".scalar(keys %MGSepaOnlyRetry),
				"staged=$CatNotPrepped", "extraction_needed=$dirsNOTPrepped");
			$nextAuditProgress = time + 60;
		}
	}
	$PhylosExist = 0 if ($CatFileMiss/scalar(@specis) > 0.1); #only activate if more than 10% missing..

	print "Output dirs status: \nIncomplete tree inputs: $CatFileMiss, complete staged inputs: $CatNotPrepped, Dir not done: $dirsNOTPrepped, phylo absent: $treeAbsent, Dir done: $doneDirs, completion-marker fast paths: $completedTreeFastPaths, too few samples: $tooFewDirs, no recoverable loci: $noRecoverableLociDirs, Phylo complete: $PhylosExist \n";
	#die;
	return($dirsNOTPrepped , $CatFileMiss , $CatNotPrepped , $treeAbsent, $doneDirs, $PhylosExist,
		$noRecoverableLociDirs, $completedTreeFastPaths);
}

sub epaOnlyRetryReady {
	my ($mgsDirectory, $activeOOMRetry) = @_;
	my $retryModeAllowed = $activeOOMRetry || (
		$onlySubmit && !$recalcTrees && !$reSubmit
			&& !$repairCAT && !$deepRepair && !$redoSubmissionData
	);
	return '' unless $retryModeAllowed
		&& $strictBackbone && $phyloProg == 1;
	return '' unless defined($mgsDirectory) && -d $mgsDirectory;
	my $pending = File::Spec->catfile($mgsDirectory, 'placementPending.sto');
	my $terminal = File::Spec->catfile($mgsDirectory, 'noTree.sto');
	my $finalTree = File::Spec->catfile(
		$mgsDirectory, 'phylo', 'IQtree_allsites.treefile');
	my $jplace = File::Spec->catfile(
		$mgsDirectory, 'phylo', 'epa-ng', 'epa_result.jplace');
	# A retained jplace with no final tree needs only normal publication/filter
	# continuation. EPA-only recovery is reserved for placement that never
	# produced a reusable jplace.
	return '' if -s $finalTree || -s $terminal || -s $jplace;
	my @required = (
		File::Spec->catfile($mgsDirectory, 'MSA', 'MSAli.fna'),
		File::Spec->catfile($mgsDirectory, 'MSA', 'MSAli.placement.fna'),
		File::Spec->catfile($mgsDirectory, 'phylo', 'IQtree_allsites.backbone.treefile'),
		File::Spec->catfile($mgsDirectory, 'phylo', 'IQtree_allsites.backbone.log'),
		File::Spec->catfile($mgsDirectory, 'phylo', 'strict_backbone.samples.tsv'),
	);
	return '' unless fileGZe($required[0]) && fileGZe($required[1]);
	return '' if grep { !-s $_ } @required[2 .. $#required];
	return 'legacy_missing_final' unless -s $pending;
	open my $marker, '<', $pending or return 'legacy_missing_final';
	my $placementPending = 0;
	while (my $line = <$marker>) {
		if ($line =~ /^status\tplacement_pending\s*$/) {
			$placementPending = 1;
			last;
		}
	}
	close $marker or return 'legacy_missing_final';
	return $placementPending ? 'explicit_pending' : 'legacy_missing_final';
}

sub prepareEpaOnlyRetryState {
	my ($mgsDirectory, $state) = @_;
	my $finalTree = File::Spec->catfile(
		$mgsDirectory, 'phylo', 'IQtree_allsites.treefile');
	return 0 if -s $finalTree;
	my $completion = File::Spec->catfile($mgsDirectory, 'treeDone.sto');
	# A completion stone without the final non-backbone tree is stale and would
	# otherwise make BuildTree refuse the isolated recovery run.
	retry_unlink($completion, label => 'clear stale completion missing final placed tree')
		if -e $completion;
	return 1 if ($state // '') eq 'explicit_pending';
	my $pending = File::Spec->catfile($mgsDirectory, 'placementPending.sto');
	my $temporary = "$pending.tmp.$$";
	retry_unlink($temporary, fatal => 0,
		label => 'clear legacy placement marker temporary');
	my $marker = retry_open('>', $temporary,
		label => 'create legacy placement-pending marker');
	print {$marker} join("\n",
		"status\tplacement_pending",
		"reason\tlegacy run retained a strict backbone but has no final non-backbone tree",
		"retry_mode\tepa_only",
	), "\n" or die "Cannot write legacy placement marker $temporary: $!\n";
	retry_close($marker, 'close legacy placement-pending marker');
	retry_rename($temporary, $pending,
		label => 'publish legacy placement-pending marker');
	print "  Legacy placement recovery: final IQtree_allsites.treefile is absent; "
		."prepared an isolated BuildTree EPA retry.\n";
	return 1;
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
        my $qc   = $OQstrH{$MGS} // "";

        my $outD = "$baseOut/$MGS";
        make_path($outD) unless -d $outD;

        my $FNAtf = "$outD/$FNAstdof$suffix";
        my $FAAtf = "$outD/$FAAstdof$suffix";
        my $CATtf = "$outD/$CATstdof.tmp$suffix";
        my $QCtf = "$outD/$QCstdof.tmp$suffix";

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
		if (length $qc) {
			open my $fh_qc, ">>", $QCtf or die $!;
			print {$fh_qc} $qc or die "Cannot append $QCtf: $!\n";
			close $fh_qc or die "Cannot close $QCtf: $!\n";
		}

        $OFstrH{$MGS} = "";
        $OAstrH{$MGS} = "";
        $OCstrH{$MGS} = "";
        $OLstrH{$MGS} = "";
        $OQstrH{$MGS} = "";

        $wrMGS++;
    }

    print "wrote for $wrMGS MGS data..\n";
}


sub reportingsMGS{
	#eval #sample/MGS
	my %smplPmgs;
	foreach my $MGS (keys %smplsPerMGS){
		foreach my $sm (keys %{$smplsPerMGS{$MGS}}){
			$smplPmgs{$MGS}++ if ($smplsPerMGS{$MGS}{$sm} >= $MGStoolowGsThr);
		}
	}
	my @smplNs = values(%smplPmgs);@smplNs = sort { $b <=> $a}  @smplNs;
	if (@smplNs) {
		my $qt50=quantileArray(0.5,@smplNs);my $qt90=quantileArray(0.90,@smplNs);
		my @top = @smplNs[0 .. ($#smplNs < 4 ? $#smplNs : 4)];
		print "Samples/MGS: QTL 50,90: $qt50 $qt90 . Top 5: ".join(" ",@top)."\n";
	} else {
		print "Samples/MGS: none have at least $MGStoolowGsThr candidate loci\n";
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

sub stepComplete {
	my ($step, $started, @statistics) = @_;
	writeStrainWorkflowHeartbeat($step);
	my $elapsed = timeNice(time - $started);
	my $globalElapsed = timeNice(time - $^T);
	my $details = @statistics ? "; ".join(", ", @statistics) : "";
	print "STEP COMPLETE: $step (${elapsed}; global elapsed $globalElapsed)$details\n";
}
sub stepProgress {
	my ($step, $completed, $total, $started, @statistics) = @_;
	writeStrainWorkflowHeartbeat($step);
	my $scope = defined($total) && $total > 0
		? "$completed/$total" : "$completed";
	my $details = @statistics ? "; ".join(", ", @statistics) : "";
	print "STEP PROGRESS: $step $scope; step elapsed "
		.timeNice(time - $started)."; global elapsed ".timeNice(time - $^T)
		.$details."\n";
}
sub writeStrainWorkflowState {
	return 0 unless length($workflowStatePath);
	my $written = write_workflow_record($workflowStatePath,
		status => $workflowStatus, stage => $workflowStage,
		reason => $workflowReason);
	cleanupLegacyStrainWorkflowStateFiles() if $written;
	return $written;
}

sub cleanupLegacyStrainWorkflowStateFiles {
	for my $legacy ($legacyWorkflowHeartbeatPath, $legacyWorkflowFailurePath) {
		next unless length($legacy) && (-e $legacy || -l $legacy);
		retry_unlink($legacy, fatal => 0,
			label => "remove legacy strain workflow state $legacy");
	}
}

sub writeStrainWorkflowHeartbeat {
	my ($stage) = @_;
	$workflowStage = $stage if defined($stage) && length($stage);
	$workflowStatus = $workflowStage eq 'complete' ? 'completed' : 'running';
	$workflowReason = '';
	writeStrainWorkflowState();
}

sub writeStrainWorkflowFailure {
	my ($error) = @_;
	$error //= 'unknown failure';
	$workflowStatus = 'failed';
	$workflowReason = $error;
	writeStrainWorkflowState();
}

sub phase1WorkerCommand {
	my $strain1scr = getProgPaths("MGS_strain1_scr");
	my @selfArgs = (
		# Stage-I workers receive extraction/consensus controls only; BuildTree
		# model, alignment, submission, and placement flags stay in the parent.
		'-GCd', $GCd, '-outD', $outD, '-MGS', $MGSfileOri,
		'-clusterID', $clusterID, '-submit', 0, '-onlySubmit', 1,
		'-maxSubJob', $maxSubJob,
		'-MGSminGenesPSmpl', $MGStoolowGsThr,
		'-multiGeneSmplMax', $multiGeneSmplMax,
		'-conspGeneSmplMax', $conspGeneSmplMax,
		'-minBadLociPSmpl', $minBadLociForSampleSkip, '-MGSphylo', $treeFile,
		'-presortGenes', $presortGenes, '-maxGenes', $maxNGenes,
		'-treeLocusBudget', $treeLocusBudget,
		'-noGeneLimit', $noGeneLimit, '-disableQC', $disableQC,
		'-breakpointGeneFlank', $breakpointGeneFlank,
		'-abundanceMinLoci', $abundanceMinimumLoci,
		'-abundanceMinFold', $abundanceMinimumFold,
		'-abundanceMaxFold', $abundanceMaximumFold,
		'-abundanceMaxModifiedZ', $abundanceMaximumModifiedZ,
		'-flushEvery', $appendWriteTrigger,
		'-MGset', $useGTDBmg, '-minSNPDepth', $minSNPDepth,
		'-minSNPCallQual', $minSNPCallQual, '-forceSNPcalls', $forceVCF2FNA,
		'-preCompConsSNP', $preCompCons, '-skipIndels', $noIndels,
		'-SNPadaptiveQual', $useAdaptiveQual,
		'-SNPdepthFilterScale', $depthFilterScale,
		'-SNPindelRangeFilt', $indelRange,
	);
	push @selfArgs, ('-tmpD', $locTmpDir1) if $locTmpDir1 ne "";
	push @selfArgs, ('-mosaicLoci', $mosaicLociFile) if $mosaicLociFile ne "";
	push @selfArgs, ('-MGSabundance', $MGSabundanceOverride)
		if $MGSabundanceOverride ne "";
	my $workerMGSSubset = $recalcTrees
		? join(",", grep { $MGSneedsExtraction{$_} } @specis)
		: $subsMGSstr;
	push @selfArgs, ('-MGSsubset', $workerMGSSubset) if $workerMGSSubset ne "";
	return $strain1scr . " " . join(" ", map { shellQuote($_) } @selfArgs);
}

sub writePhase1RepairQueue {
	my ($generation, $workers, $reason) = @_;
	my $path = "$LOGDIR/phase1_worker_repair.queue.tsv";
	my $handle = retry_open('>', "$path.tmp.$$", label => 'create Phase-I repair queue');
	print {$handle} join("\t", qw(generation worker reason)), "\n";
	for my $worker (@{$workers || []}) {
		my (undef, $workerReason) = validatePhase1WorkerLedger($worker, $generation);
		$workerReason ||= $reason || 'validation failed';
		$workerReason =~ s/[\t\r\n]+/ /g;
		print {$handle} join("\t", $generation, $worker, $workerReason), "\n";
	}
	retry_close($handle, 'close Phase-I repair queue');
	retry_rename("$path.tmp.$$", $path, label => 'publish Phase-I repair queue');
	return $path;
}

sub validatePhase1WorkerLedger {
	my ($worker, $generation) = @_;
	my $stone = "$splitStonePrefix.$worker.stone";
	my $stoneGeneration = '';
	if (-s $stone && open(my $stoneFH, '<', $stone)) {
		$stoneGeneration = <$stoneFH> // '';
		close $stoneFH;
		chomp $stoneGeneration;
	}
	return (0, 'missing or generation-mismatched completion stone')
		unless $stoneGeneration eq $generation;
	my @ledgers = (
		["$LOGDIR/$recoveryLogName.$worker",
		 join("\t", qw(MGS sample outcome reason retained_genes qc_status
			ambiguous_failure conspecific_failure recovered_mosaic_loci)), 'recovery'],
		["$LOGDIR/$sampleStatsLogName.$worker", join("\t", @sampleStatColumns), 'sample'],
	);
	for my $ledger (@ledgers) {
		my ($path, $expectedHeader, $kind) = @{$ledger};
		return (0, "missing $kind ledger $path") unless -s $path;
		open my $input, '<', $path or return (0, "unreadable $kind ledger $path: $!");
		my $header = <$input> // ''; chomp $header;
		unless ($header eq $expectedHeader) {
			close $input;
			return (0, "unexpected $kind ledger header in $path");
		}
		# The worker completion stone proves that its producer closed these ledgers.
		# Full row/cardinality validation is performed while merging each ledger;
		# rescanning every row here doubled Phase-I recovery I/O.
		close $input or return (0, "cannot close $kind ledger $path: $!");
	}
	my $conspecific = "$LOGDIR/ConspecificMGS.$worker.log";
	return (0, "missing conspecific ledger $conspecific") unless -e $conspecific;
	return (1, '');
}

sub phase1WorkersNeedingRetry {
	my ($generation) = @_;
	my @failed;
	for my $worker (0 .. $maxSubJob - 1) {
		my ($valid, $reason) = validatePhase1WorkerLedger($worker, $generation);
		next if $valid;
		limitedWarn('invalid Phase-I worker output',
			"Phase-I worker $worker requires retry: $reason\n");
		push @failed, $worker;
	}
	return @failed;
}

sub writeRecoveryRow {
	my (@fields) = @_;
	die "MAG recovery log is not open\n" unless $recoveryLogFH;
	print {$recoveryLogFH} join("\t", @fields), "\n"
		or die "Cannot write MAG recovery statistics: $!\n";
}

sub indexRecoveryRow {
	my ($worker, $line, $source) = @_;
	my $copy = $line;
	$copy =~ s/[\r\n]+\z//;
	return unless length $copy;
	my @field = split /\t/, $copy, -1;
	die "Malformed MAG recovery row in $source: expected at least 9 tab-delimited fields\n"
		unless @field >= 9;
	my ($mgs, $sample, $outcome, undef, $retained_genes) = @field[0 .. 4];
	return unless $outcome eq "recovered";
	die "Malformed retained-gene count '$retained_genes' for $mgs/$sample in $source\n"
		unless $retained_genes =~ /^\d+\z/;
	$recoveryWorkersByMGS{$mgs}{$worker} = 1;
	$recoveryRecordsByMGS{$mgs} += $retained_genes;
	$recoveryRowsByMGS{$mgs}++;
	$recoveryWorkerRecordsByMGS{$mgs}{$worker} += $retained_genes;
	$recoveryWorkerRowsByMGS{$mgs}{$worker}++;
	$recoverySamplesByMGS{$mgs}{$sample} = 1;
}
sub recoverCompletedSplitPhaseI {
	# A previous main worker can end after every extraction worker has published
	# its completion stone, but before it merges their ledgers.  The aggregate
	# FNA/FAA/category files then look reusable to the resume audit while the
	# contribution index required to validate their worker parts is absent.
	# Recover only a proven-complete generation: never merge partial retries.
	return 0 unless $maxSubJob && !$subJob;
	if ($dirsNOTPrepped == 0) {
		limitedNotice('tree-only resume skips obsolete Phase-I ledger validation',
			"Tree-only resume: every MGS input passed the completed audit; skipping obsolete Phase-I worker-ledger validation and continuing to Phase II.\n");
		retry_unlink("$LOGDIR/phase1_worker_repair.queue.tsv", fatal => 0,
			label => "clear obsolete Phase-I repair queue for tree-only resume");
		return 0;
	}
	return 0 unless split_generation_complete(
		$splitManifest, $splitStonePrefix, $maxSubJob,
	);
	my $manifestHandle = retry_open('<', $splitManifest,
		label => 'read completed Phase-I generation');
	my $manifestLine = <$manifestHandle> // '';
	retry_close($manifestHandle, 'close completed Phase-I generation');
	chomp $manifestLine;
	die "Malformed completed Phase-I generation manifest: $splitManifest\n"
		unless $manifestLine =~ /^([A-Za-z0-9_.:-]+)\t\Q$maxSubJob\E$/;
	my $generation = $1;
	my @recoveryParts = map { "$LOGDIR/$recoveryLogName.$_" } 0 .. $maxSubJob - 1;
	my @sampleStatsParts = map { "$LOGDIR/$sampleStatsLogName.$_" } 0 .. $maxSubJob - 1;
	my $hasRecoveryParts = grep { -e $_ } @recoveryParts;
	my $hasSampleStatsParts = grep { -e $_ } @sampleStatsParts;
	my @failedWorkers = phase1WorkersNeedingRetry($generation);
	my $workerRetryRound = 0;
	while (@failedWorkers) {
		my $cannotRepair = grep { $_ == 0 } @failedWorkers;
		$cannotRepair ||= !$doSubmit || $workerRetryRound >= $phase1WorkerRetries;
		if ($cannotRepair) {
			my $queue = writePhase1RepairQueue($generation, \@failedWorkers,
				'resumed Phase-I worker validation failed');
			$completionMessage = "Phase I requires worker repair before Phase II; no tree jobs were submitted.";
			print "Phase-I recovery paused safely; repair queue: $queue. Invalid workers: ".join(',', @failedWorkers)."\n";
			exit(0);
		}
		$workerRetryRound++;
		print "Resubmitting invalid Phase-I worker(s) from completed generation: "
			.join(',', @failedWorkers)." (round $workerRetryRound/$phase1WorkerRetries)\n";
		my @retryJobs;
		my $savedTmp = $QSBoptHR->{tmpSpace};
		$QSBoptHR->{tmpSpace} = 15;
		my $workerCommand = phase1WorkerCommand();
		for my $worker (@failedWorkers) {
			my $stone = "$splitStonePrefix.$worker.stone";
			retry_unlink($stone, fatal => 0, label => "clear worker $worker completion");
			my $cmdX = "$workerCommand -subjob $worker &&\n"
				."printf '%s\\n' ".shellQuote($generation)." > ".shellQuote($stone)."\n";
			my ($dependency) = qsubSystem("$LOGDIR/Strain1_B${worker}.resume${workerRetryRound}.sh",
				$cmdX,1,"${selfMemGb}G","Str1.$worker","","",1,[],$QSBoptHR);
			push @retryJobs, $dependency;
		}
		$QSBoptHR->{tmpSpace} = $savedTmp;
		qsubSystemJobAlive(\@retryJobs, $QSBoptHR);
		@failedWorkers = phase1WorkersNeedingRetry($generation);

	}
	retry_unlink("$LOGDIR/phase1_worker_repair.queue.tsv", fatal => 0,
		label => 'clear obsolete Phase-I repair queue');
	$hasRecoveryParts = grep { -e $_ } @recoveryParts;
	$hasSampleStatsParts = grep { -e $_ } @sampleStatsParts;

	return 0 unless $hasRecoveryParts || $hasSampleStatsParts;

	my @missingRecovery = grep { !-s $_ } @recoveryParts;
	my @missingStats = grep { !-s $_ } @sampleStatsParts;
	die "Completed split Phase I has incomplete recovery ledgers: ".join(',', @missingRecovery)."\n"
		if @missingRecovery;
	die "Completed split Phase I has incomplete sample-statistics ledgers: ".join(',', @missingStats)."\n"
		if @missingStats;

	print "Recovering completed split Phase I: merging $maxSubJob worker ledgers before tree submission.\n";
	mergeConspecificLogs();
	mergeRecoveryLogs() unless $recoveryContributionIndexReady;
	mergeSampleStats();
	return 1;
}

sub mergeSampleStats {
	make_path($outD) unless -d $outD;
	my @parts = $maxSubJob
		? map { "$LOGDIR/$sampleStatsLogName.$_" } 0 .. $maxSubJob - 1
		: ("$LOGDIR/$sampleStatsLogName.0");
	my @missing = grep { !-s $_ } @parts;
	die "Missing per-worker sample statistics: ".join(',', @missing)."\n" if @missing;

	my $expectedHeader = join("\t", @sampleStatColumns);
	my $final = "$LOGDIR/$sampleStatsLogName";
	my $finalTemporary = "$final.write.$$";
	my $merged = retry_open('>', $finalTemporary,
		label => 'create merged sample statistics');
	print {$merged} $expectedHeader, "\n"
		or die "Cannot write $finalTemporary: $!\n";
	my (@allRows, %rowsByWorker, %seenSample);
	for my $worker (0 .. $#parts) {
		my $part = $parts[$worker];
		open my $in, '<', $part or die "Cannot read $part: $!\n";
		my $header = <$in>;
		die "Per-worker sample statistics have no header: $part\n" unless defined $header;
		$header =~ s/[\r\n]+\z//;
		die "Unexpected per-worker sample-statistics header in $part\n"
			unless $header eq $expectedHeader;
		my $lineNumber = 1;
		while (my $line = <$in>) {
			$lineNumber++;
			$line =~ s/[\r\n]+\z//;
			die "Empty sample-statistics row in $part line $lineNumber\n" unless length $line;
			my @values = split /\t/, $line, -1;
			die "Wrong sample-statistics field count in $part line $lineNumber: got "
				.scalar(@values).", expected ".scalar(@sampleStatColumns)."\n"
				unless @values == @sampleStatColumns;
			my %row;
			@row{@sampleStatColumns} = @values;
			my $sample = $row{sample} // '';
			die "Sample-statistics row has no sample in $part line $lineNumber\n"
				unless length $sample;
			die "Duplicate sample-statistics row for $sample across workers\n"
				if $seenSample{$sample}++;
			die "Sample-statistics worker mismatch for $sample: row=$row{worker}, file=$worker\n"
				unless defined($row{worker}) && $row{worker} =~ /^\d+\z/
					&& $row{worker} == $worker;
			print {$merged} $line, "\n" or die "Cannot write $finalTemporary: $!\n";
			push @allRows, \%row;
			push @{$rowsByWorker{$worker}}, \%row;
		}
		close $in or die "Cannot close $part: $!\n";
	}
	retry_close($merged, 'close merged sample statistics');
	die "No sample-statistics rows were recovered from worker tables\n" unless @allRows;

	my @summaryColumns = sample_summary_columns();
	my @summaryRows = map {
		aggregate_sample_rows($rowsByWorker{$_} || [], "worker.$_")
	} 0 .. $#parts;
	$summaryRows[$_]->{workers} = 1 for 0 .. $#parts;
	my $allSummary = aggregate_sample_rows(\@allRows, 'ALL');
	$allSummary->{workers} = scalar(@parts);
	push @summaryRows, $allSummary;
	my $summary = "$LOGDIR/$sampleStatsSummaryLogName";
	my $summaryTemporary = "$summary.write.$$";
	my $summaryFH = retry_open('>', $summaryTemporary,
		label => 'create sample-statistics summary');
	print {$summaryFH} join("\t", @summaryColumns), "\n"
		or die "Cannot write $summaryTemporary: $!\n";
	for my $row (@summaryRows) {
		my @values = map {
			my $value = defined($row->{$_}) ? $row->{$_} : '';
			$value =~ s/[\t\r\n]+/ /g;
			$value;
		} @summaryColumns;
		print {$summaryFH} join("\t", @values), "\n"
			or die "Cannot write $summaryTemporary: $!\n";
	}
	retry_close($summaryFH, 'close sample-statistics summary');
	retry_rename($finalTemporary, $final, label => 'publish merged sample statistics');
	retry_rename($summaryTemporary, $summary, label => 'publish sample-statistics summary');

	warn "Per-sample statistics cover ".scalar(@allRows)." of ".scalar(@samples)
		." configured samples\n" if @samples && @allRows != @samples;
	printSampleStatsSummary($allSummary);
	print "Merged per-sample statistics: $final\n";
	print "Per-worker and all-worker summary: $summary\n";
	return $allSummary;
}

sub reportSavedSampleStats {
	my $summary = "$LOGDIR/$sampleStatsSummaryLogName";
	my $legacySummary = "$outD/$sampleStatsSummaryLogName";
	$summary = $legacySummary if !-s $summary && -s $legacySummary;
	unless (-s $summary) {
		warn "Phase I statistics summary is unavailable at $summary; continuing without the saved sample histogram\n";
		return 0;
	}
	my @summaryColumns = sample_summary_columns();
	my $expectedHeader = join("\t", @summaryColumns);
	open my $in, '<', $summary or die "Cannot read saved sample summary $summary: $!\n";
	my $header = <$in> // '';
	$header =~ s/[\r\n]+\z//;
	unless ($header eq $expectedHeader) {
		close $in or warn "Cannot close incompatible saved sample summary $summary: $!\n";
		warn "Saved Phase-I sample summary uses an older schema at $summary; "
			."continuing tree recovery without replaying its optional histogram\n";
		return 0;
	}
	my $allSummary;
	while (my $line = <$in>) {
		$line =~ s/[\r\n]+\z//;
		next unless length($line);
		my @values = split /\t/, $line, -1;
		die "Wrong saved sample-summary field count in $summary\n"
			unless @values == @summaryColumns;
		my %row;
		@row{@summaryColumns} = @values;
		next unless ($row{scope} // '') eq 'ALL';
		die "Duplicate ALL row in saved sample summary $summary\n" if $allSummary;
		$allSummary = \%row;
	}
	close $in or die "Cannot close saved sample summary $summary: $!\n";
	die "Saved sample summary has no ALL row: $summary\n" unless $allSummary;
	print "Reusing completed Phase I sample accounting: $summary\n";
	printSampleStatsSummary($allSummary);
	return 1;
}

sub printSampleStatsSummary {
	my ($allSummary) = @_;
	die "Sample summary must be a hash reference\n" unless ref($allSummary) eq 'HASH';
	my @humanColumns = grep { $_ ne "used_mgs_loci_histogram" } sample_summary_columns();
	my @summaryPairs = map {
		my $value = defined($allSummary->{$_}) ? $allSummary->{$_} : "";
		$value =~ s/\s+/_/g;
		"$_:$value";
	} @humanColumns;
	print "STEP 1 SAMPLE SUMMARY (all workers)\n";
	print join(" ", @summaryPairs), "\n";
	my @histogramRows = loci_histogram_rows(
		$allSummary->{used_mgs_loci_histogram}, $allSummary->{min_genes_per_mgs}
	);
	my $largestBin = 0;
	for my $row (@histogramRows) {
		$largestBin = $row->[1] if $row->[1] > $largestBin;
	}
	print "Used MGS retained-loci histogram (MGS-sample observations):\n";
	for my $row (@histogramRows) {
		my ($label, $count) = @$row;
		my $fraction = $allSummary->{used_mgs}
			? 100 * $count / $allSummary->{used_mgs} : 0;
		my $barWidth = $largestBin ? int(30 * $count / $largestBin + 0.5) : 0;
		$barWidth = 1 if $count && !$barWidth;
		printf "  %-10s %8d %6.2f%% %s\n", $label, $count, $fraction, "#" x $barWidth;
	}
}

sub writeRecoveryContributionIndex {
	my $path = "$LOGDIR/$recoveryLogName.contributors.tsv";
	my $temporary = "$path.write.$$";
	my $out = retry_open('>', $temporary,
		label => 'create recovery contribution index');
	print {$out} join("\t", qw(MGS worker recovered_rows retained_records unique_samples)), "\n"
		or die "Cannot write $temporary: $!\n";
	for my $mgs (sort keys %recoveryWorkersByMGS) {
		my $unique = scalar keys %{$recoverySamplesByMGS{$mgs} || {}};
		$recoveryUniqueSamplesByMGS{$mgs} = $unique;
		die "Recovery ledger contains duplicate recovered samples for $mgs\n"
			unless $unique == ($recoveryRowsByMGS{$mgs} // 0);
		for my $worker (sort { $a <=> $b } keys %{$recoveryWorkersByMGS{$mgs}}) {
			print {$out} join("\t", $mgs, $worker,
				$recoveryWorkerRowsByMGS{$mgs}{$worker} // 0,
				$recoveryWorkerRecordsByMGS{$mgs}{$worker} // 0, $unique), "\n"
				or die "Cannot write $temporary: $!\n";
		}
	}
	retry_close($out, 'close recovery contribution index');
	retry_rename($temporary, $path, label => 'publish recovery contribution index');
}

sub loadRecoveryContributionIndex {
	my $path = "$LOGDIR/$recoveryLogName.contributors.tsv";
	return 0 unless -s $path;
	open my $in, '<', $path or die "Cannot read $path: $!\n";
	my $header = <$in> // '';
	chomp $header;
	die "Unsupported recovery contribution index header in $path\n"
		unless $header eq join("\t", qw(MGS worker recovered_rows retained_records unique_samples));
	while (my $line = <$in>) {
		chomp $line;
		next unless length $line;
		my ($mgs, $worker, $rows, $records, $unique) = split /\t/, $line, -1;
		die "Malformed recovery contribution index row in $path: $line\n"
			unless defined($mgs) && length($mgs) && defined($unique)
				&& $worker =~ /^\d+\z/ && $rows =~ /^\d+\z/
				&& $records =~ /^\d+\z/ && $unique =~ /^\d+\z/;
		die "Duplicate recovery contribution index row for $mgs worker $worker\n"
			if $recoveryWorkersByMGS{$mgs}{$worker};
		$recoveryWorkersByMGS{$mgs}{$worker} = 1;
		$recoveryWorkerRowsByMGS{$mgs}{$worker} = $rows;
		$recoveryWorkerRecordsByMGS{$mgs}{$worker} = $records;
		$recoveryRowsByMGS{$mgs} += $rows;
		$recoveryRecordsByMGS{$mgs} += $records;
		if (exists $recoveryUniqueSamplesByMGS{$mgs}) {
			die "Inconsistent unique-sample count for $mgs in $path\n"
				unless $recoveryUniqueSamplesByMGS{$mgs} == $unique;
		} else {
			$recoveryUniqueSamplesByMGS{$mgs} = $unique;
		}
	}
	close $in or die "Cannot close $path: $!\n";
	for my $mgs (keys %recoveryRowsByMGS) {
		die "Recovery contribution index has duplicate samples for $mgs\n"
			unless $recoveryRowsByMGS{$mgs} == $recoveryUniqueSamplesByMGS{$mgs};
	}
	$recoveryContributionIndexReady = 1;
	warn "Loaded recovery contribution index: $path\n";
	return 1;
}

sub mergeRecoveryLogs {
	make_path($LOGDIR) unless -d $LOGDIR;
	my @parts = $maxSubJob
		? map { "$LOGDIR/$recoveryLogName.$_" } 0 .. $maxSubJob - 1
		: ("$LOGDIR/$recoveryLogName.0");
	return unless grep { -e $_ } @parts;
	my @missing = grep { !-s $_ } @parts;
	die "Missing MAG recovery worker log(s): ".join(',', @missing)."\n" if @missing;
	my $final = "$LOGDIR/$recoveryLogName";
	my $temporary = "$final.write.$$";
	my $out = retry_open('>', $temporary,
		label => 'create merged MAG recovery ledger');
	my $header_written = 0;
	for my $worker (0 .. $#parts) {
		my $part = $parts[$worker];
		open my $in, '<', $part or die "Cannot read $part: $!\n";
		my $header = <$in>;
		die "MAG recovery worker log has no header: $part\n" unless defined $header;
		$header =~ s/[\r\n]+\z//;
		my $expectedHeader = join("\t", qw(MGS sample outcome reason retained_genes
			qc_status ambiguous_failure conspecific_failure recovered_mosaic_loci));
		die "Unexpected MAG recovery header in $part\n"
			unless $header eq $expectedHeader;
		$header .= "\n";
		print {$out} $header unless $header_written++;
		while (my $line = <$in>) { indexRecoveryRow($worker, $line, $part); print {$out} $line or die "Cannot write $temporary: $!\n"; }
		close $in or die "Cannot close $part: $!\n";
	}
	retry_close($out, 'close merged MAG recovery ledger');
	writeRecoveryContributionIndex();
	retry_rename($temporary, $final, label => 'publish merged MAG recovery ledger');
	retry_unlink($_, fatal => 0, label => "clean merged recovery ledger") for @parts;
	$recoveryContributionIndexReady = 1;
	print "MAG recovery accounting: $final\n";
}

sub writeSelectionAttritionSummary {
	my ($recoveryMetrics, $filterReasons) = @_;
	$recoveryMetrics ||= {};
	$filterReasons ||= {};
	my @rows;
	push @rows, ['recovery', $_, $recoveryMetrics->{$_}]
		for sort keys %{$recoveryMetrics};
	push @rows, ['recovery', "filtered_reason.$_", $filterReasons->{$_}]
		for sort keys %{$filterReasons};

	my $sampleSummary = "$LOGDIR/$sampleStatsSummaryLogName";
	my $legacySampleSummary = "$outD/$sampleStatsSummaryLogName";
	$sampleSummary = $legacySampleSummary
		if !-s $sampleSummary && -s $legacySampleSummary;
	if (-s $sampleSummary) {
		open my $sampleInput, '<', $sampleSummary
			or die "Cannot read sample attrition summary $sampleSummary: $!\n";
		my $header = <$sampleInput> // '';
		$header =~ s/[\r\n]+\z//;
		my @columns = split /\t/, $header, -1;
		my $allRow;
		while (my $line = <$sampleInput>) {
			$line =~ s/[\r\n]+\z//;
			next unless length $line;
			my @values = split /\t/, $line, -1;
			next unless @values == @columns;
			my %row;
			@row{@columns} = @values;
			if (($row{scope} // '') eq 'ALL') {
				$allRow = \%row;
				last;
			}
		}
		close $sampleInput
			or die "Cannot close sample attrition summary $sampleSummary: $!\n";
		if ($allRow) {
			for my $metric (qw(
				candidate_mgs candidate_loci pre_abundance_loci post_abundance_loci
				missing_consensus_loci low_depth_loci breakpoint_loci csp_rejected_loci
				ambiguous_loci abundance_filtered_loci invalid_protein_loci retained_loci
				capped_mgs capped_loci skipped_within_2_loci_of_min
				skip_no_selected_loci skip_no_usable_loci
				skip_too_few_after_abundance skip_too_few_valid_sequences
			)) {
				push @rows, ['extraction', $metric, $allRow->{$metric}]
					if exists $allRow->{$metric};
			}
		}
	}

	my (%treeTotals, %treeMetricReports);
	my $treeReports = 0;
	for my $mgs (@specis) {
		next unless defined($SIdirs{$mgs}) && length($SIdirs{$mgs});
		my $report = File::Spec->catfile(
			$SIdirs{$mgs}, 'phylo', 'selection_attrition.tsv');
		next unless -s $report;
		open my $treeInput, '<', $report
			or die "Cannot read tree selection attrition $report: $!\n";
		my $header = <$treeInput> // '';
		$header =~ s/[\r\n]+\z//;
		die "Unexpected tree selection attrition header in $report\n"
			unless $header eq "metric\tvalue";
		$treeReports++;
		while (my $line = <$treeInput>) {
			$line =~ s/[\r\n]+\z//;
			next unless length $line;
			my ($metric, $value) = split /\t/, $line, 2;
			next if !defined($metric) || $metric eq 'schema';
			next unless defined($value)
				&& $value =~ /\A-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\z/;
			$treeTotals{$metric} += 0 + $value;
			$treeMetricReports{$metric}++;
		}
		close $treeInput
			or die "Cannot close tree selection attrition $report: $!\n";
	}
	push @rows, ['tree', 'reports_expected', scalar(@specis)];
	push @rows, ['tree', 'reports_available', $treeReports];
	for my $metric (sort keys %treeTotals) {
		push @rows, ['tree', $metric, $treeTotals{$metric}];
		push @rows, ['tree', "$metric.reports", $treeMetricReports{$metric}];
	}

	my $path = "$LOGDIR/strainSelectionAttrition.tsv";
	my $temporary = "$path.write.$$";
	open my $output, '>', $temporary
		or die "Cannot create strain selection attrition $temporary: $!\n";
	print {$output} "scope\tmetric\tvalue\n"
		or die "Cannot write strain selection attrition header: $!\n";
	for my $row (@rows) {
		my @value = map { defined($_) ? $_ : '' } @{$row};
		print {$output} join("\t", @value), "\n"
			or die "Cannot write strain selection attrition row: $!\n";
	}
	close $output or die "Cannot close strain selection attrition $temporary: $!\n";
	rename $temporary, $path
		or die "Cannot publish strain selection attrition $path: $!\n";
	return $path;
}

sub writeMGSSampleHistograms {
	my @records;
	for my $mgs (@specis) {
		next unless defined($SIdirs{$mgs}) && length($SIdirs{$mgs});
		my ($backbone, $placement, $excluded, $source);
		my $attrition = File::Spec->catfile(
			$SIdirs{$mgs}, 'phylo', 'selection_attrition.tsv');
		if (-s $attrition) {
			open my $input, '<', $attrition
				or die "Cannot read MGS sample attrition $attrition: $!\n";
			my $header = <$input> // '';
			$header =~ s/[\r\n]+\z//;
			my %metric;
			if ($header eq "metric\tvalue") {
				while (my $line = <$input>) {
					$line =~ s/[\r\n]+\z//;
					my ($name, $value) = split /\t/, $line, 2;
					next unless defined($name) && defined($value)
						&& $value =~ /\A\d+\z/;
					$metric{$name} = 0 + $value;
				}
			}
			close $input or die "Cannot close MGS sample attrition $attrition: $!\n";
			if (exists($metric{backbone_samples}) && exists($metric{placement_samples})) {
				($backbone, $placement, $excluded, $source) = (
					$metric{backbone_samples}, $metric{placement_samples},
					$metric{excluded_samples} // 0, 'selection_attrition',
				);
			} elsif (exists($metric{final_samples})) {
				($backbone, $placement, $excluded, $source) = (
					$metric{final_samples}, 0, 0, 'selection_attrition.final_samples',
				);
			}
		}
		if (!defined($backbone) || !defined($placement)) {
			my $classification = File::Spec->catfile(
				$SIdirs{$mgs}, 'phylo', 'strict_backbone.samples.tsv');
			if (-s $classification) {
				open my $input, '<', $classification
					or die "Cannot read MGS sample classification $classification: $!\n";
				my $header = <$input> // '';
				$header =~ s/[\r\n]+\z//;
				my @columns = split /\t/, $header, -1;
				my %column = map { $columns[$_] => $_ } 0 .. $#columns;
				if (exists($column{sample}) && exists($column{tree_role})) {
					my (%seen, %roleCount);
					while (my $line = <$input>) {
						$line =~ s/[\r\n]+\z//;
						next unless length($line);
						my @value = split /\t/, $line, -1;
						my $sample = $value[$column{sample}] // '';
						my $role = $value[$column{tree_role}] // '';
						next unless length($sample) && !$seen{$sample}++;
						$roleCount{$role}++ if $role =~ /\A(?:backbone|placement|excluded)\z/;
					}
					($backbone, $placement, $excluded, $source) = (
						$roleCount{backbone} // 0, $roleCount{placement} // 0,
						$roleCount{excluded} // 0, 'strict_backbone.samples.tsv',
					);
				}
				close $input
					or die "Cannot close MGS sample classification $classification: $!\n";
			}
		}
		next unless defined($backbone) && defined($placement);
		my $finalTree = File::Spec->catfile(
			$SIdirs{$mgs}, 'phylo', 'IQtree_allsites.treefile');
		my $treeStatus = -s $finalTree ? 'complete'
			: -s File::Spec->catfile($SIdirs{$mgs}, 'placementPending.sto')
				? 'placement_pending' : 'tree_missing';
		push @records, {
			mgs => $mgs, backbone => $backbone, placement => $placement,
			excluded => $excluded // 0, included => $backbone + $placement,
			source => $source, tree_status => $treeStatus,
		};
	}

	my $detailPath = "$LOGDIR/strainMGSSampleCounts.tsv";
	my $detailTemporary = "$detailPath.write.$$";
	my $detail = retry_open('>', $detailTemporary,
		label => 'create per-MGS included-sample counts');
	print {$detail} join("\t", qw(
		MGS backbone_samples placement_samples included_samples excluded_samples
		tree_status source
	)), "\n" or die "Cannot write $detailTemporary: $!\n";
	for my $record (sort { $a->{mgs} cmp $b->{mgs} } @records) {
		print {$detail} join("\t", @{$record}{qw(
			mgs backbone placement included excluded tree_status source
		)}), "\n" or die "Cannot write $detailTemporary: $!\n";
	}
	retry_close($detail, 'close per-MGS included-sample counts');
	retry_rename($detailTemporary, $detailPath,
		label => 'publish per-MGS included-sample counts');

	my @bins = (
		[0, 0, '0'], [1, 2, '1-2'], [3, 4, '3-4'], [5, 9, '5-9'],
		[10, 19, '10-19'], [20, 49, '20-49'], [50, 99, '50-99'],
		[100, 199, '100-199'], [200, 499, '200-499'],
		[500, 999, '500-999'], [1000, 1999, '1000-1999'],
		[2000, undef, '2000+'],
	);
	my $histogramPath = "$LOGDIR/strainMGSSampleHistogram.tsv";
	my $histogramTemporary = "$histogramPath.write.$$";
	my $histogram = retry_open('>', $histogramTemporary,
		label => 'create across-MGS sample histogram');
	print {$histogram} "role\tlower\tupper\tbin\tMGS_count\tfraction\n"
		or die "Cannot write $histogramTemporary: $!\n";
	my %statistics;
	for my $role (qw(backbone placement)) {
		my @values = map { $_->{$role} } @records;
		my @counts = (0) x scalar(@bins);
		for my $value (@values) {
			for my $index (0 .. $#bins) {
				my ($lower, $upper) = @{$bins[$index]};
				next if $value < $lower;
				next if defined($upper) && $value > $upper;
				$counts[$index]++;
				last;
			}
		}
		my $maximumBin = @counts ? (sort { $b <=> $a } @counts)[0] : 0;
		print ucfirst($role), " samples per MGS (", scalar(@values), " MGS):\n";
		for my $index (0 .. $#bins) {
			my ($lower, $upper, $label) = @{$bins[$index]};
			my $fraction = @values ? $counts[$index] / @values : 0;
			print {$histogram} join("\t", $role, $lower,
				defined($upper) ? $upper : '', $label, $counts[$index],
				sprintf('%.6f', $fraction)), "\n"
				or die "Cannot write $histogramTemporary: $!\n";
			next unless $counts[$index];
			my $barWidth = $maximumBin
				? int(30 * $counts[$index] / $maximumBin + 0.5) : 0;
			$barWidth = 1 if !$barWidth;
			printf "  %-10s %7d %6.2f%% %s\n", $label, $counts[$index],
				100 * $fraction, '#' x $barWidth;
		}
		$statistics{$role} = {
			minimum => @values ? (sort { $a <=> $b } @values)[0] : 0,
			maximum => @values ? (sort { $b <=> $a } @values)[0] : 0,
			median => @values ? median(@values) : 0,
			mean => @values ? mean(@values) : 0,
		};
	}
	retry_close($histogram, 'close across-MGS sample histogram');
	retry_rename($histogramTemporary, $histogramPath,
		label => 'publish across-MGS sample histogram');
	print "Across-MGS sample histogram: $histogramPath\n";
	print "Per-MGS sample counts: $detailPath\n";
	return {
		histogram => $histogramPath, details => $detailPath,
		mgs_count => scalar(@records), statistics => \%statistics,
	};
}

sub writeStrainSummary {
	my ($tree_disposition, $mosaic_outgroups_used) = @_;
	my %used_mosaic_outgroup = %{$mosaic_outgroups_used || {}};
	for my $mgs (@specis) {
		next unless exists($PreferredOutgroup{$mgs});
		for my $path ("$SIdirs{$mgs}/data.log", "$scratchD/outs/$mgs/data.log") {
			next unless fileGZe($path);
			my ($fh) = gzipopen($path, "Mosaic outgroup usage log");
			my $line = <$fh> // '';
			close $fh;
			chomp $line;
			$line =~ s/^OG://;
			$used_mosaic_outgroup{$mgs} = 1 if length($line) && $line eq $PreferredOutgroup{$mgs};
			last if exists($used_mosaic_outgroup{$mgs});
		}
	}
	my $recovery = "$LOGDIR/$recoveryLogName";
	my $legacyRecovery = "$outD/$recoveryLogName";
	$recovery = $legacyRecovery if !-s $recovery && -s $legacyRecovery;
	my ($evaluated, $recovered, $filtered, $gene_sum, $mosaic_loci) = (0, 0, 0, 0, 0);
	my (@recovered_genes, %filter_reason, %recovered_status, %represented_samples, %represented_mgs);
	if (-s $recovery) {
		open my $fh, '<', $recovery or die "Cannot read $recovery: $!\n";
		my $header = <$fh>;
		while (my $line = <$fh>) {
			chomp $line;
			next unless length $line;
			my ($mgs, $sample, $outcome, $reason, $genes, $status, $ambiguous, $conspecific, $row_mosaic_loci) = split /\t/, $line, -1;
			$evaluated++;
			$represented_samples{$sample} = 1;
			$represented_mgs{$mgs} = 1;
			if ($outcome eq 'recovered') {
				$recovered++;
				$recovered_status{$status || 'unspecified'}++;
				$gene_sum += $genes;
				$mosaic_loci += $row_mosaic_loci || 0;
				push @recovered_genes, 0 + $genes;
			} else {
				$filtered++;
				$filter_reason{$reason || 'unspecified'}++;
			}
		}
		close $fh or die "Cannot close $recovery: $!\n";
	}
	my $average = $recovered ? $gene_sum / $recovered : 0;
	my $median_genes = @recovered_genes ? median(@recovered_genes) : 0;
	my @thresholds = (10, 50, 100, 200, 500, 1000, 2000);
	my %above;
	for my $threshold (@thresholds) { $above{$threshold} = scalar(grep { $_ > $threshold } @recovered_genes); }
	my $selectionAttrition = writeSelectionAttritionSummary({
		evaluated_sample_mgs => $evaluated,
		recovered_mags => $recovered,
		filtered_mags => $filtered,
		recovered_loci => $gene_sum,
		average_loci_per_recovered_mag => sprintf('%.6f', $average),
		median_loci_per_recovered_mag => $median_genes,
		recovered_mosaic_loci => $mosaic_loci,
	}, \%filter_reason);
	my $sampleHistograms = writeMGSSampleHistograms();
	my @lines = (
		"Strain-within recovery summary (v$version)",
		"output_directory\t$outD",
		"recovery_accounting\t".(-s $recovery ? $recovery : 'not_available'),
		"selection_attrition\t$selectionAttrition",
		"MGS_sample_counts\t$sampleHistograms->{details}",
		"MGS_sample_histogram\t$sampleHistograms->{histogram}",
		"MGS_with_sample_counts\t$sampleHistograms->{mgs_count}",
		"input_samples\t".scalar(@samples),
		"usable_samples\t".(scalar(@samples) - scalar(keys %unavailableSamples)),
		"unavailable_samples\t".scalar(keys %unavailableSamples),
		"selected_MGS\t".scalar(@specis),
		"evaluated_sample_MGS\t$evaluated",
		"recovered_MAGs\t$recovered",
		"filtered_MAGs\t$filtered",
		sprintf("average_genes_per_recovered_MAG\t%.2f", $average),
		"median_genes_per_recovered_MAG\t$median_genes",
		"recovered_mosaic_loci\t$mosaic_loci",
		"mosaic_outgroups_used\t".scalar(keys %used_mosaic_outgroup),
		"samples_with_evaluated_MAGs\t".scalar(keys %represented_samples),
		"MGS_with_evaluated_samples\t".scalar(keys %represented_mgs),
	);
	for my $role (qw(backbone placement)) {
		my $stats = $sampleHistograms->{statistics}{$role};
		push @lines,
			"${role}_samples_per_MGS.minimum\t$stats->{minimum}",
			"${role}_samples_per_MGS.median\t$stats->{median}",
			sprintf("${role}_samples_per_MGS.mean\t%.2f", $stats->{mean}),
			"${role}_samples_per_MGS.maximum\t$stats->{maximum}";
	}
	push @lines, map { "recovered_MAGs.genes_gt_$_\t$above{$_}" } @thresholds;
	push @lines, map { "recovered_status.$_\t$recovered_status{$_}" } sort keys %recovered_status;
	push @lines, map { "filtered_reason.$_\t$filter_reason{$_}" } sort keys %filter_reason;
	push @lines, map { "tree_disposition.$_\t$tree_disposition->{$_}" } sort keys %{$tree_disposition};
	my $summary = "$LOGDIR/$summaryLogName";
	my $temporary = "$summary.write.$$";
	open my $out, '>', $temporary or die "Cannot create $temporary: $!\n";
	print {$out} join("\n", @lines), "\n" or die "Cannot write $temporary: $!\n";
	close $out or die "Cannot close $temporary: $!\n";
	rename $temporary, $summary or die "Cannot install $summary: $!\n";
	print "\n", join("\n", @lines), "\nStatistics log\t$summary\n";
}

sub printEarlyRunHeader {
	my $earlyMode = length($MGSfile) ? 'MGS' : 'FMG';
	my $inputDirectory = length($MGSfile) ? dirname($MGSfile) : $GCd;
	my $requestedOutput = length($outDpre)
		? $outDpre : File::Spec->catdir($inputDirectory, 'intra_phylo');
	my $selected = select(STDOUT);
	$| = 1;
	select($selected);
	print "============= Strain_within v$version =============\n";
	print "Started: ".scalar(localtime())."\n";
	print "Mode: $earlyMode".($subJob ? "; split worker $subJob/$maxSubJob" : '')."\n";
	print "GC dir: $GCd\n";
	print "MGS input: ".(length($MGSfile) ? $MGSfile : '(FMG mode)')."\n";
	print "Requested output: $requestedOutput\n";
	print "Cores: $numCores (max: $maxCores); submit=$doSubmit; "
		."onlySubmit=$onlySubmit; recalcTrees=$recalcTrees; redoSubmissionData=$redoSubmissionData; "
		."redoEPAfilter=$redoEPAfilter\n";
	print "Tree OOM recovery: rounds=$treeOOMRetryRounds; maximum memory=${treeOOMMaxMemGB}GB\n";
	print "Initializing paths, maps, and catalogues...\n";
	print "==============================================\n";
}

sub cleanupMosaicIntermediates {
	my ($prefix) = @_;
	my @paths = (
		(map { $prefix.$_ } qw(
			.minimap2.paf
			.rtk.mosaic.tsv
			.rtk.mosaic.summary.tsv
			.rtk.concat.list
			.rejected.tsv
			.outgroups.tsv
		)),
		(map { File::Spec->catfile(dirname($prefix), $_) } qw(
			prepare_mosaic_loci.log
			prepare_mosaic_loci.sh
			prepare_mosaic_loci.sh.otxt
			prepare_mosaic_loci.sh.etxt
		)),
	);
	my $removed = 0;
	for my $path (@paths) {
		next unless -e $path;
		if (unlink $path) {
			$removed++;
		} else {
			warn "Cannot remove obsolete Mosaic intermediate $path: $!\n";
		}
	}
	print "Removed $removed obsolete Mosaic intermediate(s)\n" if $removed;
}

sub shellQuote {
	my ($value) = @_;
	$value = "" unless defined $value;
	$value =~ s/'/'"'"'/g;
	return "'$value'";
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


sub lifecycleMarkerReason {
	my ($marker, $fallback) = @_;
	my $reason = defined($fallback) ? $fallback : '';
	return $reason unless defined($marker) && -s $marker;
	if (open(my $markerFH, '<', $marker)) {
		while (my $line = <$markerFH>) {
			if ($line =~ /^reason\t(.*)/) {
				$reason = $1;
				last;
			}
		}
		close $markerFH;
	}
	$reason =~ s/[\t\r\n]+/ /g;
	return length($reason) ? $reason : (defined($fallback) ? $fallback : '');
}

sub writeTreeFailureAudit {
	my ($expected) = @_;
	my $path = "$LOGDIR/tree_job_outcomes.tsv";
	my $temporary = "$path.write.$$";
	my $output = retry_open('>', $temporary, label => 'create tree-job outcome audit');
	print {$output} "MGS\tstatus\ttree\tcompletion_marker\treason\n"
		or die "Cannot write tree-job outcome audit $temporary: $!\n";
	my (@failed, @pending, @terminal);
	for my $mgs (sort keys %{$expected}) {
		my ($tree, $stone, $terminalMarker, $pendingMarker) = @{$expected->{$mgs}};
		my ($status, $reason) = ('failed_missing_output', '');
		if (-s $tree && -s $stone) {
			$status = 'complete';
		} elsif (-s $terminalMarker) {
			$status = 'valid_no_tree';
			push @terminal, $mgs;
		} elsif (-s $pendingMarker) {
			$status = 'placement_pending';
			push @pending, $mgs;
		} else {
			push @failed, $mgs;
		}
		my $marker = $status eq 'valid_no_tree' ? $terminalMarker
			: $status eq 'placement_pending' ? $pendingMarker : '';
		$reason = lifecycleMarkerReason($marker, '') if $marker;
		print {$output} join("\t", $mgs, $status, $tree, $stone, $reason), "\n"
			or die "Cannot write tree-job outcome audit $temporary: $!\n";
	}
	retry_close($output, 'close tree-job outcome audit');
	retry_rename($temporary, $path, label => 'publish tree-job outcome audit');
	print "Tree-job outcome audit: $path\n";
	return (\@failed, \@pending, \@terminal);
}
sub dispatchPendingTreeJobs {
	my %args = @_;
	my $queue = $args{queue} || [];
	my $options = $args{options} || {};
	my $jobs = $args{jobs} || [];
	my $accounting = $args{accounting} || [];
	my $blocking = $args{blocking} ? 1 : 0;
	my $submitted = 0;

	return { submitted => 0, pending => 0 } unless @{$queue};
	my $saved_nonblocking = $options->{nonblockingMaxConcurrentJobs};
	my $saved_tmp_space = $options->{tmpSpace};
	my $saved_long_queue = $options->{useLongQueue};
	delete $options->{capacityDeferred};
	delete $options->{capacityDeferralAnnounced};
	$options->{nonblockingMaxConcurrentJobs} = 1 unless $blocking;

	while (@{$queue}) {
		my $record = $queue->[0];
		for my $required (qw(mgs script command cores memory job_name tree stone terminal placement_pending)) {
			die "Queued tree job is missing '$required'\n"
				unless defined($record->{$required}) && length($record->{$required});
		}
		if ($options->{doSubmit}) {
			my @staleOutputs = $record->{epa_only}
				? qw(stone terminal) : qw(stone tree terminal placement_pending);
			retry_unlink($record->{$_}, label => "clear stale tree-job $_")
				for @staleOutputs;
		}
		$options->{tmpSpace} = $record->{tmp_space};
		$options->{useLongQueue} = $record->{use_long_queue};
		my ($dependency) = qsubSystem(
			$record->{script}, $record->{command}, $record->{cores},
			$record->{memory}, $record->{job_name}, "", "", 1, [], $options,
		);
		if (defined($dependency) && $dependency eq deferredSubmissionDependency()) {
			print "Scheduler capacity is full; retaining ".scalar(@{$queue})
				." prepared tree job(s) while Phase II continues converting inputs.\n";
			last;
		}
		shift @{$queue};
		$submitted++ if $options->{doSubmit};
		push @{$jobs}, $dependency if defined($dependency) && length($dependency);
		if ($options->{doSubmit} && ($options->{qmode} || '') eq 'slurm'
				&& defined($dependency)) {
			my $schedulerJobID = $dependency;
			$schedulerJobID =~ s/^\Q$options->{rTag}\E//;
			push @{$accounting}, {
				job_id => $schedulerJobID, mgs => $record->{mgs},
				requested_mb => $record->{requested_mb},
				retry_round => $record->{retry_round} // 0,
				submission_record => { %{$record} },
			} if $schedulerJobID =~ /^\d+$/;
		}
	}
	$options->{nonblockingMaxConcurrentJobs} = $saved_nonblocking;
	$options->{tmpSpace} = $saved_tmp_space;
	$options->{useLongQueue} = $saved_long_queue;
	return { submitted => $submitted, pending => scalar(@{$queue}) };
}

sub retryOOMTreeJobs {
	my %args = @_;
	my $accounting = $args{accounting} || [];
	my $options = $args{options} || {};
	my $maximumMB = $args{maximum_mb};
	my $maximumRounds = $args{maximum_rounds} // 3;
	return 0 unless @{$accounting};
	return 0 unless $options->{doSubmit} && ($options->{qmode} || '') eq 'slurm';
	die "OOM retry rounds must be between 0 and 3\n"
		unless $maximumRounds >= 0 && $maximumRounds <= 3;

	my @roundAccounting = @{$accounting};
	my $summary = slurm_tree_memory_summary(\@roundAccounting);
	print format_slurm_tree_memory_summary($summary);
	my $retried = 0;
	for my $round (1 .. $maximumRounds) {
		last unless $summary->{available};
		my @oom = @{$summary->{oom_jobs} || []};
		last unless @oom;
		my @retryQueue;
		for my $oom (@oom) {
			my $original = $oom->{submission_record};
			unless (ref($original) eq 'HASH') {
				warn "Cannot retry OOM tree job $oom->{job_id}: submission record is unavailable\n";
				next;
			}
			my $nextMB = next_oom_retry_memory_mb(
				$original->{requested_mb}, $maximumMB);
			unless (defined($nextMB)) {
				warn "OOM retry ceiling reached for $original->{mgs}: "
					."$original->{requested_mb} MB already meets -treeOOMMaxMemGB "
					."$treeOOMMaxMemGB\n";
				next;
			}
			my %retry = %{$original};
			my $mgsDirectory = dirname($retry{terminal});
			my $epaState = epaOnlyRetryReady($mgsDirectory, 1);
			my $epaStage = $retry{epa_only} || length($epaState);
			if ($epaStage) {
				unless (length($epaState)) {
					warn "Cannot isolate EPA OOM retry for $retry{mgs}: retained placement state is incomplete\n";
					next;
				}
				unless (prepareEpaOnlyRetryState($mgsDirectory, $epaState)) {
					print "Skipping OOM retry for $retry{mgs}: final placed tree is now complete.\n";
					next;
				}
				$retry{epa_only} = 1;
				$retry{cores} = 1;
				$retry{command} =~ s/(^|\s)-epaThreads\s+\d+/$1-epaThreads 1/;
				$retry{command} =~ s/(^|\s)-cores\s+\d+/$1-cores 1/;
				unless ($retry{command} =~ /(?:^|\s)-epaOnly\s+1(?:\s|$)/) {
					$retry{command} =~ s/\s+\z//;
					$retry{command} .= " -epaOnly 1\n";
				}
				$retry{script} = File::Spec->catfile(
					$mgsDirectory, 'treeCmd.epa_retry.sh');
			}
			$retry{retry_round} = $round;
			$retry{requested_mb} = $nextMB;
			$retry{memory} = $nextMB.'M';
			$retry{job_name} = 'OOM'.$round.'.'.$retry{mgs};
			$retry{script} = File::Spec->catfile(
				$mgsDirectory, "treeCmd.oom_retry.$round.sh") unless $epaStage;
			push @retryQueue, \%retry;
			print "OOM retry round $round/$maximumRounds for $retry{mgs}: "
				."$original->{requested_mb} MB -> $nextMB MB; "
				.($epaStage ? 'EPA-only with 1 thread' : "$retry{cores} core full-tree resume")
				."\n";
		}
		last unless @retryQueue;
		my (@retryJobs, @retryAccounting);
		my $drain = dispatchPendingTreeJobs(
			queue => \@retryQueue, options => $options,
			jobs => \@retryJobs, accounting => \@retryAccounting,
			blocking => 1,
		);
		die "Internal error: OOM retry queue was not drained\n" if @retryQueue;
		last unless $drain->{submitted} && @retryAccounting;
		qsubSystemJobAlive(\@retryJobs, $options) if @retryJobs;
		push @{$accounting}, @retryAccounting;
		$retried += scalar(@retryAccounting);
		@roundAccounting = @retryAccounting;
		$summary = slurm_tree_memory_summary(\@roundAccounting);
		print "OOM retry round $round accounting:\n";
		print format_slurm_tree_memory_summary($summary);
	}
	if ($summary->{available} && @{$summary->{oom_jobs} || []}) {
		warn "Tree OOM recovery stopped after at most $maximumRounds retry round(s); "
			."remaining OOM outcomes stay quarantined for inspection\n";
	}
	return $retried;
}

sub resetMGSTreeOutputs {
	my ($mgsDir, $MGS) = @_;
	die "Cannot reset tree outputs for an absent MGS directory: $mgsDir\n"
		unless -d $mgsDir;
	my $resolvedRoot = abs_path($outD)
		or die "Cannot resolve within-strain output directory $outD: $!\n";
	my $resolvedMGS = abs_path($mgsDir)
		or die "Cannot resolve MGS output directory $mgsDir: $!\n";
	die "Refusing to reset tree outputs outside the selected within-strain directory: $resolvedMGS\n"
		unless dirname($resolvedMGS) eq $resolvedRoot
			&& basename($resolvedMGS) eq $MGS;

	my $phyloDir = File::Spec->catdir($resolvedMGS, "phylo");
	my $treeStone = File::Spec->catfile($resolvedMGS, "treeDone.sto");
	my $terminalMarker = File::Spec->catfile($resolvedMGS, "noTree.sto");
	my $placementMarker = File::Spec->catfile($resolvedMGS, "placementPending.sto");
	if (-d $phyloDir) {
		remove_tree($phyloDir, {safe => 1});
		die "Cannot completely remove tree output directory $phyloDir\n" if -e $phyloDir;
	}
	retry_unlink($treeStone, label => "remove tree completion checkpoint");
	print "  Reset tree outputs: removed $phyloDir and tree completion checkpoint\n";
	retry_unlink($terminalMarker, label => "remove terminal no-tree marker");
	retry_unlink($placementMarker, label => "remove placement-pending marker");
}

sub completionMarkerTree {
	my ($marker, $output_directory) = @_;
	return '' unless defined($marker) && -s $marker
		&& defined($output_directory) && -d $output_directory;
	open my $input, '<', $marker or return '';
	my $line = <$input> // '';
	close $input or return '';
	$line =~ s/[\r\n]+\z//;
	my ($producer, $marker_version, $tree_path) = split /\t/, $line, 3;
	return '' unless defined($producer) && $producer eq 'buildTree5'
		&& defined($marker_version) && $marker_version =~ /^\d+(?:\.\d+)?\z/
		&& $marker_version >= 5.40
		&& defined($tree_path) && length($tree_path);
	my $output = File::Spec->canonpath(File::Spec->rel2abs($output_directory));
	my $tree = File::Spec->canonpath(File::Spec->rel2abs($tree_path, $output));
	my $relative = File::Spec->abs2rel($tree, $output);
	return '' if $relative eq File::Spec->curdir
		|| $relative =~ /^\.\.(?:[\\\/]|\z)/;
	return -s $tree ? $tree : '';
}

sub directResumeStagedInputsReady {
	my ($scratch_directory, $mgs, $script) = @_;
	return 0 unless defined($scratch_directory) && length($scratch_directory)
		&& defined($mgs) && length($mgs) && defined($script) && -s $script;
	open my $input, '<', $script or return 0;
	my $script_text = do { local $/; <$input> // '' };
	close $input or return 0;
	return 0 unless $script_text =~ /(?:^|\s)-stagedInputDir(?:\s|=)/;
	my $staged_directory = File::Spec->catdir($scratch_directory, 'outs', $mgs);
	return !grep {
		!fileGZe(File::Spec->catfile($staged_directory, $_))
	} ($FNAstdof, $FAAstdof, $CATstdof);
}

sub resubmitExistingTreeCommands {
	my %args = @_;
	my $outdir = $args{outdir} // '';
	my $force = $args{force} ? 1 : 0;
	my $redoEpa = $args{redo_epa} ? 1 : 0;
	my $scratch_directory = $args{scratch_directory} // '';
	my $subset = $args{subset} || [];
	my $options = $args{options} || {};
	return (0, 0) unless -d $outdir && $options->{doSubmit};
	print "Resubmitting only phylogenies\n";

	my %requested;
	for my $mgs (@{$subset}) {
		return (0, 0) unless defined($mgs)
			&& $mgs =~ /\A[A-Za-z0-9][A-Za-z0-9_.:+-]*\z/;
		$requested{$mgs} = 1;
	}
	if (!%requested) {
		for my $script (bsd_glob(File::Spec->catfile($outdir, '*', 'treeCmd.sh'))) {
			my $mgs = basename(dirname($script));
			next unless $mgs =~ /\A[A-Za-z0-9][A-Za-z0-9_.:+-]*\z/;
			$requested{$mgs} = 1;
		}
	}
	return (0, 0) unless %requested;

	my @scripts;
	for my $mgs (sort keys %requested) {
		my $mgs_dir = File::Spec->catdir($outdir, $mgs);
		unless (-d $mgs_dir) {
			limitedWarn('direct resume missing MGS directory',
				"Skipping $mgs: saved-command output directory is absent: $mgs_dir\n");
			next;
		}
		next if -s File::Spec->catfile($mgs_dir, 'noTree.sto');
		my $finalTree = File::Spec->catfile(
			$mgs_dir, 'phylo', 'IQtree_allsites.treefile');
		my $treeDone = File::Spec->catfile($mgs_dir, 'treeDone.sto');
		my $completedTree = completionMarkerTree($treeDone, $mgs_dir);
		if (!length($completedTree)) {
			for my $candidate (map {
				File::Spec->catfile($mgs_dir, 'phylo', $_)
			} qw(IQtree_allsites.treefile VERYFASTTREE_allsites.nwk FASTTREE_allsites.nwk)) {
				if (-s $candidate) {
					$completedTree = $candidate;
					last;
				}
			}
		}
		next if !$force && -s $treeDone && length($completedTree);
		my $pending = File::Spec->catfile($mgs_dir, 'placementPending.sto');
		my $publicationResume = !$force && !-s $finalTree
			&& -s File::Spec->catfile(
				$mgs_dir, 'phylo', 'IQtree_allsites.backbone.treefile')
			&& -s File::Spec->catfile(
				$mgs_dir, 'phylo', 'epa-ng', 'epa_result.jplace');
		my ($script, $mode) = (
			File::Spec->catfile($mgs_dir, 'treeCmd.sh'),
			$redoEpa ? 'redo_epa' : 'full',
		);
		my $stagedInputsReady = directResumeStagedInputsReady(
			$scratch_directory, $mgs, $script);
		if ($redoEpa && !$publicationResume) {
			limitedWarn('redo EPA filter missing retained publication state',
				"Skipping $mgs: -redoEPAfilter requires its retained backbone and jplace\n");
			next;
		} elsif (!$publicationResume && !$force && -s $pending) {
			my $retry_script = File::Spec->catfile($mgs_dir, 'treeCmd.epa_retry.sh');
			$script = $retry_script if -s $retry_script;
			$mode = 'epa_only';
		} elsif (!$publicationResume) {
			my @missing = grep {
				!fileGZe(File::Spec->catfile($mgs_dir, $_))
			} ($FNAstdof, $FAAstdof, $CATstdof);
			if (@missing && !$stagedInputsReady) {
				limitedWarn('direct resume missing tree input',
					"Skipping $mgs: saved full-tree command lacks "
					.join(', ', @missing)." and has no complete staged-input triplet\n");
				next;
			}
			print "Direct resume $mgs: reusing validated staged tree inputs\n" if @missing;
		}
		unless (-s $script) {
			limitedWarn('direct resume missing saved command',
				"Skipping $mgs: saved tree command is absent or empty: $script\n");
			next;
		}
		push @scripts, [$mgs, $script, $mode];
	}

	print "Direct tree-command resume: ".scalar(@scripts)
		." saved treeCmd.sh job(s); skipping Mosaic, map, and catalogue loading.\n";
	for my $record (@scripts) {
		# Saved full-tree commands need stale final outputs cleared; EPA-only
		# retries retain their pending marker as the BuildTree recovery contract.
		my $mgs_dir = dirname($record->[1]);
		my @stale = (File::Spec->catfile($mgs_dir, 'treeDone.sto'));
		if ($record->[2] eq 'epa_only') {
			push @stale, File::Spec->catfile($mgs_dir, 'noTree.sto');
		} else {
			push @stale, File::Spec->catfile($mgs_dir, 'placementPending.sto'),
				map { File::Spec->catfile($mgs_dir, 'phylo', $_) }
					qw(IQtree_allsites.treefile VERYFASTTREE_allsites.nwk FASTTREE_allsites.nwk);
		}
		for my $stale (@stale) {
			retry_unlink($stale, label => 'clear stale direct-resume tree output') if -e $stale;
		}
		qsubSystemWaitMaxJobs(
			$options->{maxConcurrentJobs} || 0,
			$options->{killDependencyNever} || 0, $options,
		);
		my $submission;
		if ($record->[2] eq 'redo_epa') {
			local $ENV{MATAFILER_REDO_EPA_FILTER} = 1;
			$submission = qsubSystem2($record->[1], $options);
		} else {
			$submission = qsubSystem2($record->[1], $options);
		}
		print "  Resubmitted $record->[0] from $record->[1]: $submission";
	}
	return (1, scalar(@scripts));
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





sub newSampleStats {
	my ($sample, $status, $candidateMGS, $candidateLoci, $assemblyGroup) = @_;
	my %stats = map { $_ => 0 } @sampleStatColumns;
	$stats{sample} = $sample;
	$stats{worker} = $subJob;
	$stats{assembly_group} = defined($assemblyGroup) ? $assemblyGroup : $sample;
	$stats{status} = $status;
	$stats{selected_mgs} = scalar(@specis);
	$stats{candidate_mgs} = $candidateMGS;
	$stats{candidate_loci} = $candidateLoci;
	$stats{min_genes_per_mgs} = $MGStoolowGsThr;
	$stats{presort_genes} = $presortGenes;
	$stats{max_genes} = $noGeneLimit ? q{unlimited} : $maxNGenes;
	$stats{tree_locus_budget} = $treeLocusBudget;
	$stats{qc_enabled} = $disableQC ? 0 : 1;
	$stats{min_gene_depth} = $minDepthGene;
	$stats{min_bad_loci} = $minBadLociForSampleSkip;
	$stats{multi_gene_fraction_max} = $multiGeneSmplMax;
	$stats{csp_gene_fraction_max} = $conspGeneSmplMax;
	$stats{csp_locus_score_max} = $conspecificSpThr;
	$stats{breakpoint_gene_flank} = $breakpointGeneFlank;
	$stats{abundance_min_loci} = $abundanceMinimumLoci;
	$stats{abundance_min_fold} = $abundanceMinimumFold;
	$stats{abundance_max_fold} = $abundanceMaximumFold;
	$stats{abundance_max_modified_z} = $abundanceMaximumModifiedZ;
	$stats{used_mgs_loci_histogram} = encode_loci_histogram([], $MGStoolowGsThr);
	return \%stats;
}

sub writeSampleStats {
	my ($fh, $stats, $seen) = @_;
	my $sample = $stats->{sample} // '';
	die "Cannot emit per-sample statistics without a sample name\n" unless length $sample;
	die "Per-sample statistics would contain a duplicate row for $sample\n"
		if $seen && $seen->{$sample}++;
	my @values = map {
		my $value = defined($stats->{$_}) ? $stats->{$_} : q{};
		$value =~ s/[\t\r\n]+/ /g;
		$value;
	} @sampleStatColumns;
	my $row = join("\t", @values);
	die "Refusing to emit an empty per-sample statistics row for $sample\n"
		unless $row =~ /\S/;
	for my $target ($fh, $sampleStatsPartFH) {
		next unless $target;
		print {$target} $row, "\n"
			or die "Cannot write per-sample statistics: $!\n";
	}
}

#this routine hast to get genes out of each sample, that are needed
#and save them to be later written per specI
sub extractFNAFAA2genes{
	my $recovery_part = $maxSubJob
		? "$LOGDIR/$recoveryLogName.$subJob"
		: "$LOGDIR/$recoveryLogName.0";
	open $recoveryLogFH, '>', $recovery_part
		or die "Cannot create MAG recovery log $recovery_part: $!\n";
	print {$recoveryLogFH} join("\t", qw(
		MGS sample outcome reason retained_genes qc_status
		ambiguous_failure conspecific_failure recovered_mosaic_loci
	)), "\n";
	my $sampleStatsHeader = join("\t", @sampleStatColumns);
	my $sample_stats_part = "$LOGDIR/$sampleStatsLogName.$subJob";
	open $sampleStatsPartFH, '>', $sample_stats_part
		or die "Cannot create per-worker sample statistics $sample_stats_part: $!\n";
	print {$sampleStatsPartFH} $sampleStatsHeader, "\n"
		or die "Cannot write per-worker sample-statistics header: $!\n";

	# Each worker owns one numeric suffix.  A retry must replace, not append to,
	# that worker's previous partial extraction.
	for my $pattern (
		"$scratchD/outs/*/$FNAstdof.$subJob",
		"$scratchD/outs/*/$FAAstdof.$subJob",
		"$scratchD/outs/*/$LINKstdof.$subJob",
		"$scratchD/outs/*/$CATstdof.tmp.$subJob",
		"$scratchD/outs/*/$QCstdof.tmp.$subJob",
	) {
		for my $part (bsd_glob($pattern)) {
			retry_unlink($part, label => "remove stale worker part");
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
		if ($perMGSgenes < $MGStoolowGsThr){
			$lowCandidateMGS++;
			limitedWarn("MGS with fewer than $MGStoolowGsThr candidate loci",
				"Only $perMGSgenes genes/COGs for MGS $MGS; MGS genes might be multi-copy\n")
				unless $maxSubJob;
		}
	}
	print "$lowCandidateMGS MGS have fewer than $MGStoolowGsThr candidate loci in this worker's sample slice; "
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
		my @preview = @srtdSmpls > 10 ? @srtdSmpls[0 .. 9] : @srtdSmpls;
		print "\nSUBJOB ${subJob}/$maxSubJob: pre-restricted to " . scalar(@srtdSmpls)
			. " sample driver(s) with target loci"
			. (@preview ? ": ".join(' ', @preview) : '')
			. (scalar(@srtdSmpls) > @preview ? " ..." : "") . "\n\n";
	}
	
	
	
	print "Extracting GC genes from " . scalar(@srtdSmpls). " (of " . scalar(keys(%cl2gene2)) . ") ASsembly Groups\n";

	
	#different way to go over genes..
	 my $smCnt=1;
	 #storage hash for raw fasta/faa/link files, needs to be written separately
	#goes over every assembly group to extract SNP corrected genes that fall into each MGS
	my $writeLink = 1; my $appCnt=0;
		#DEBUG	@srtdSmpls = ("PDB3.F");
	
	
	my $previousFH = select(STDOUT);
	$| = 1;
	select($previousFH);
	open my $sampleStatsFH, q{>&}, \*STDOUT
		or die "Cannot duplicate STDOUT for per-sample statistics: $!\n";
	$previousFH = select($sampleStatsFH);
	$| = 1;
	select($previousFH);
	print {$sampleStatsFH} $sampleStatsHeader, "\n"
		or die "Cannot write per-sample statistics header: $!\n";
	my %sampleStatsSeen;
	{
		# Redirect diagnostics once for the complete accumulation loop. Reopening
		# STDOUT for every assembly group produced empty scheduler records on some
		# systems. The duplicated handle above remains the TSV-only stream.
		local *STDOUT;
		open STDOUT, q{>&}, \*STDERR
			or die "Cannot redirect sample diagnostics to STDERR: $!\n";
		foreach my $sm (@srtdSmpls){
			print STDERR "AT SMPL:: $smCnt/" . scalar(@srtdSmpls) ." $sm - ". "Elapsed time : ", timeNice(time - $sttime) . "\n";
			readGenesSample_Singl(
				$sm, $writeLink, $sttime, \$appCnt, $sampleStatsFH, \%sampleStatsSeen,
			);
			$smCnt++;
		}
	}
	close $sampleStatsFH or die "Cannot close per-sample statistics stream: $!\n";
	close $sampleStatsPartFH
		or die "Cannot close per-worker sample statistics $sample_stats_part: $!\n";
	undef $sampleStatsPartFH;
	warn "Per-sample statistics emitted: ".scalar(keys %sampleStatsSeen)." nonempty row(s)\n";
	
	
	appendWriteMGSgenes($writeLink);
	close $recoveryLogFH or die "Cannot close MAG recovery log $recovery_part: $!\n";
	undef $recoveryLogFH;
	print "Done writing all genes to subdirs, elapsed time: " . timeNice(time - $sttime)  . "\n";
	$appCnt=0;
	#done at the point with gene extractions
	return;
}

sub reduceSeqTech{
	my ($inST) = @_;
	#		if (platform != "ill" && platform != "PB" && platform != "ONT" && platform != "unspecified") {

	return $inST if ($inST eq "ill" || $inST eq "ONT" || $inST eq "PB" || $inST eq "unspecified");
	
	return "ill" if ($inST eq "hiSeq" || $inST eq "miSeq");
	return "unspecified";
	
}

sub createConsFastas{
	my ($cD,$sm, $oFNA, $oFAA,$append2LOG,$returnCmd) = @_;
	my $vcf2fnaBin = getProgPaths("vcf2fna");
	# VCF normalization is intentionally not part of this workflow.
	my $vcf2fnaOpt = "";
	#my $seqPlatf = "hiSeq"; #-> get this from .map ..
	my $refFA = getAssemblContigs($cD); my $refGFF = getAssemblGFF($cD);
	my $depthFile = "$cD$lMAPdir/$sm$bamDepthFsuffix";
	my $ofasCons = "$cD/$lSNPdir/$lConsCTG";
	my $vcfFile = "$cD/$lSNPdir/$lConsVCF";
	my $inputVCF = $vcfFile;
	
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
		
		$seqPlatf = reduceSeqTech($seqPlatf);
		$vcf2fnaOpt = "-seqPlatform ".shellQuote($seqPlatf)." $commonOpt";
		$cmd = "$vcf2fnaBin $vcf2fnaOpt -ref ".shellQuote($refFA)
			." -inVCF ".shellQuote($inputVCF)." -depthF ".shellQuote($depthFile)."  ";
	} else {
		#die;
		#in case of both PacBio and illumina:
		#$vcf2fnaOpt = "-seqPlatform $SNPIHR->{SeqTech},$SNPIHR->{SeqTechSuppl} -t 1 -minCallDepth $minDepth,$minDepth -minCallQual $minCallQual ";
		#$cmd = "$vcf2fnaBin $vcf2fnaOpt -ref $refFA -inVCF $vcfFile,$vcfFileS -depthF $depthFile,$depthFileS ";# -oCtg $ofasCons.gz " ;
		my $vcfFileS = "$cD/$lSNPdir/$lConsVCFsup";
		my $inputVCFS = $vcfFileS;
		$seqPlatf = reduceSeqTech($seqPlatf);
		$secSeqTechS = reduceSeqTech($secSeqTechS);
		my $depthFileS = "$cD$lMAPdir/$sm$bamDepthFsuffixSup";
		$vcf2fnaOpt = "-seqPlatform ".shellQuote("$seqPlatf,$secSeqTechS")." $commonOpt";
		$cmd = "$vcf2fnaBin $vcf2fnaOpt -ref ".shellQuote($refFA)
			." -inVCF ".shellQuote("$inputVCF,$inputVCFS")
			." -depthF ".shellQuote("$depthFile,$depthFileS")." -oCtg /dev/null ";
	}

	$cmd .= "-gff ".shellQuote($refGFF)." -oGeneNT ".shellQuote($oFNA)." -oGeneAA ".shellQuote($oFAA);
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
	my ($sm, $writeLink,$sttime,$bufferedSamplesRef,$sampleStatsFH,$sampleStatsSeen) = @_;
	#my %subG = %{$subGHR};#$_[0]};
	
	my %subG; my %locMGScnt;
	# This structure is read-only here; retaining the reference avoids copying
	# every locus entry at the start of each sample.
	my $locCl2G2 = $cl2gene2{$sm};

	my $noFilter = $disableQC ? 1 : 0;
	
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
	my $candidateLoci = 0;
	$candidateLoci += $_ for values %locMGScnt;
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
		my $stats = newSampleStats($sm, q{driver_map_missing}, scalar(keys %locMGScnt), $candidateLoci);
		$stats->{skipped_mgs} = $stats->{candidate_mgs};
		writeSampleStats($sampleStatsFH, $stats, $sampleStatsSeen);
		return;
	}
	my @subGKs = keys %subG;
	unless (@subGKs) {
		limitedWarn('samples without candidate genes',
			"No candidate genes found for sample $sm; sample will be skipped\n");
		my $stats = newSampleStats($sm, q{no_candidate_genes}, scalar(keys %locMGScnt), $candidateLoci);
		$stats->{skipped_mgs} = $stats->{candidate_mgs};
		writeSampleStats($sampleStatsFH, $stats, $sampleStatsSeen);
		return;
	}
	unless ($subGKs[0] =~ m/^(.*)__/) {
		limitedWarn('unparseable catalogue members',
			"Cannot parse catalogue member '$subGKs[0]' for sample $sm; sample will be skipped\n");
		my $stats = newSampleStats($sm, q{unparseable_catalogue_member}, scalar(keys %locMGScnt), $candidateLoci);
		$stats->{skipped_mgs} = $stats->{candidate_mgs};
		writeSampleStats($sampleStatsFH, $stats, $sampleStatsSeen);
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
		my $sampleStats = newSampleStats($sd3, q{processing}, scalar(keys %locMGScnt), $candidateLoci, $sm);
		if (exists $unavailableSamples{$sd3}) {
			limitedWarn('unavailable samples', "Skipping $sd3: $unavailableSamples{$sd3}\n");
			$sampleStats->{status} = q{unavailable};
			$sampleStats->{skipped_mgs} = $sampleStats->{candidate_mgs};
			writeSampleStats($sampleStatsFH, $sampleStats, $sampleStatsSeen);
			next;
		}
		unless (exists($map{$sd3}) && defined($map{$sd3}{wrdir}) && length($map{$sd3}{wrdir})) {
			limitedWarn('samples missing map entries or working directories',
				"Skipping $sd3: missing map entry or working directory\n");
			$sampleStats->{status} = q{map_or_workdir_missing};
			$sampleStats->{skipped_mgs} = $sampleStats->{candidate_mgs};
			writeSampleStats($sampleStatsFH, $sampleStats, $sampleStatsSeen);
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
			$sampleStats->{status} = q{empty_sample};
			$sampleStats->{skipped_mgs} = $sampleStats->{candidate_mgs};
			writeSampleStats($sampleStatsFH, $sampleStats, $sampleStatsSeen);
			next;
		}
		my $rename = 0;
		$rename = 1 if ($sd2 ne $sd3);
		#print "r:$rename $sd3  (from $sd2) ";
		my $metaGD = getAssemblPath($cD,"",0);
		if ($metaGD eq ""){
			limitedWarn('samples without assemblies',
				"Assembly not available for $sd3 in $cD; sample will be skipped\n");
			$sampleStats->{status} = q{assembly_missing};
			$sampleStats->{skipped_mgs} = $sampleStats->{candidate_mgs};
			writeSampleStats($sampleStatsFH, $sampleStats, $sampleStatsSeen);
			next;
		}
		#get NT's
		#my $tar = $metaGD."genePred/genes.shrtHD.fna";
		
		#pre-calculated, as in older MATAFILER versions (pre 0.69):
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
		} elsif ($preCompCons) {
				# Split workers do not rebuild the parent's full precompute map. The
				# persistent scratch paths are deterministic, so inspect only this
				# sample's completed pair when it is actually extracted.
				my $precomputedNT = "$preConDir/$sd3.cons.genes.fna.gz";
				my $precomputedAA = "$preConDir/$sd3.cons.prots.faa.gz";
				if (fileGZe($precomputedNT) && fileGZe($precomputedAA)) {
					($fastaf, $fastafAA, $locForceVCF2FNA) =
						($precomputedNT, $precomputedAA, 0);
				}
			}
		
		my $input_state = consensusInputState(
			fileGZe($fastafVCF), fileGZe($fastaf), fileGZe($fastafAA), $locForceVCF2FNA
		);
		if ($input_state eq 'missing') {
			limitedWarn('samples without repairable consensus files',
				"Skipping $sd3: consensus NT/AA files are incomplete and no repair VCF is available\n");
			$sampleStats->{status} = q{consensus_missing};
			$sampleStats->{skipped_mgs} = $sampleStats->{candidate_mgs};
			writeSampleStats($sampleStatsFH, $sampleStats, $sampleStatsSeen);
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
			$sampleStats->{status} = q{consensus_incomplete};
			$sampleStats->{skipped_mgs} = $sampleStats->{candidate_mgs};
			writeSampleStats($sampleStatsFH, $sampleStats, $sampleStatsSeen);
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
		my $breakpointMask = {};
		unless ($noFilter) {
			my $breakpoint_file = "$cD/$lMAPdir/$sd3-smd.bam.breakpoints.tsv.gz";
			if (fileGZe($breakpoint_file)) {
				my $gff_file = getAssemblGFF($cD);
				my $breakpoint_ok = eval {
					$breakpointMask = breakpoint_gene_mask(
						$gff_file, resolveExistingFile($breakpoint_file),
						\%subG, $breakpointGeneFlank,
					);
					1;
				};
				unless ($breakpoint_ok) {
					my $error = $@ || "unknown breakpoint/GFF parsing failure";
					$error =~ s/\s+$//;
					limitedWarn('unusable breakpoint reports',
						"Could not apply breakpoint masking for $sd3: $error\n");
					$breakpointMask = {};
				}
			}
		}
		#print "Time C: " . timeNice(time - $sttime)  . "\n";

		#some stats on gene extractions..
		my $missGene=0; my $foundGene=0; my $SInum=0; my $conspGen=0;my $SNPresFail=0;
		my $doubleGenes=0; my $MGStoolowGskip=0;my $missAbundance=0;
		#stats on different ways to filter genes
		my $abundFail=0;
		my $breakpointFail=0;
		my ($cappedMGS, $cappedLoci, $nearThresholdMGS) = (0, 0, 0);
		my ($preAbundanceLoci, $postAbundanceLoci) = (0, 0);
		my ($skipNoSelected, $skipNoUsable, $skipAfterAbundance, $skipAfterSequence) = (0, 0, 0, 0);
		my $placementFlaggedMGS = 0;
		
		
		#3rd part: genes were read and renamed.. now write them out already here to save mem overall
		#currently takes too long in large GCs..
		my $MGScnt = scalar((keys %locMGScnt));
		foreach my $MGS (keys %locMGScnt) {
			# The priority list is immutable during extraction, so do not copy
			# as many as $presortGenes entries for every sample/MGS pair.
			my $COGprios1 = $COGprios->{$MGS};
			if (!$COGprios1 || !@{$COGprios1}){
				$MGStoolowGskip++;
				$skipNoSelected++;
				writeRecoveryRow($MGS, $sd3, 'filtered', 'no_selected_loci', 0, '', 0, 0, 0);
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
					if (!$noFilter && exists($breakpointMask->{$gX})) {
						$breakpointFail++;
						next;
					}
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
			my $sampleQCStatus = ($double_failure || $csp_failure) ? 'placement' : 'backbone';
			if ($double_failure || $csp_failure){
				$placementFlaggedMGS++;
				push(@{$ConspecificMGS{$MGS}}, "$sd3" ); 
				# Questionable loci were already masked above. Retain the
				# remaining validated loci, but keep this sample out of the
				# strict backbone and place it after backbone inference.
			}
			unless (@genes2) {
				$MGStoolowGskip++;
				$skipNoUsable++;
				writeRecoveryRow($MGS, $sd3, 'filtered', 'no_usable_loci', 0,
					$sampleQCStatus, $double_failure, $csp_failure, 0);
				next;
			}

			$preAbundanceLoci += scalar(@genes2);
			my $depth_mask = $noFilter ? [(1) x scalar(@abunGs)]
				: abundance_pattern_mask(\@abunGs, {
					minimum_count => $abundanceMinimumLoci,
					minimum_fold => $abundanceMinimumFold,
					maximum_fold => $abundanceMaximumFold,
					maximum_modified_z => $abundanceMaximumModifiedZ,
				});
			my @genes3=();
			for (my $i=0;$i<scalar(@abunGs);$i++){
				unless ($depth_mask->[$i]) {
					$abundFail++; next;
				}
				push (@genes3, $genes2[$i]);
			}
			
			$postAbundanceLoci += scalar(@genes3);
			if (scalar(@genes3)< $MGStoolowGsThr){
				$MGStoolowGskip++;
				$skipAfterAbundance++;
				$nearThresholdMGS++ if scalar(@genes3) >= ($MGStoolowGsThr > 2 ? $MGStoolowGsThr - 2 : 1);
				writeRecoveryRow($MGS, $sd3, 'filtered', 'too_few_after_abundance',
					scalar(@genes3), $sampleQCStatus, $double_failure, $csp_failure, 0);
				next;
			}
				
			#now write MGS into local temp storage for later tree building..
			my $locCnt=0; my $locMosaicCnt=0;
			my $wasCapped=0; my $locCappedLoci=0;
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
				if (!$noGeneLimit && $locCnt >= $maxNGenes) {
					$wasCapped = 1;
					$locCappedLoci++;
					next;
				}
				
				$locCnt++;
				#write gene out
				my $retained_group = $LocusByID->{$curLocus{$gX}};
				$locMosaicCnt++ if $retained_group && @{$retained_group->{genes} || []} > 1;
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
			
			
			if ($wasCapped) {
				$cappedMGS++;
				$cappedLoci += $locCappedLoci;
			}
			if (scalar(@OFstr) == 0 || $locCnt < $MGStoolowGsThr){ #5 genes is really too little to be considered valid as good strain rep..
				$MGStoolowGskip++;
				$skipAfterSequence++;
				$nearThresholdMGS++ if $locCnt >= ($MGStoolowGsThr > 2 ? $MGStoolowGsThr - 2 : 1);
				#delete $locMGSgenes{$MGS};
				writeRecoveryRow($MGS, $sd3, 'filtered', 'too_few_valid_sequences',
					$locCnt, $sampleQCStatus, $double_failure, $csp_failure, $locMosaicCnt);
				next;
			}
			$locMGSgenes{$MGS} = $locCnt;
			
			if (!exists($OAstrH{$MGS})){#set up base strings
				$OAstrH{$MGS} = "";$OFstrH{$MGS} = "";$OLstrH{$MGS} = "";$OCstrH{$MGS} = "";
				$OQstrH{$MGS} = "";
			}
			if ($locCnt>0){
				#save in tmp hash (faster than opening bunch of files..
				$OAstrH{$MGS} .= join("",@OAstr);$OFstrH{$MGS} .= join("",@OFstr);
				my $recovery_reason = $double_failure && $csp_failure
					? 'placement_ambiguous_and_conspecific'
					: $double_failure ? 'placement_ambiguous'
					: $csp_failure ? 'placement_conspecific' : 'passed_qc';
				writeRecoveryRow($MGS, $sd3, 'recovered', $recovery_reason,
					$locCnt, $sampleQCStatus, $double_failure, $csp_failure,
					$locMosaicCnt);
				$OLstrH{$MGS} .= join("",@OLstr);$OCstrH{$MGS} .= join("",@OCstr);
				my $ambiguous_fraction = $evaluableLoci ? $doubleCntL / $evaluableLoci : 0;
				my $csp_fraction = $evaluableLoci ? $locConSpecGen / $evaluableLoci : 0;
				$OQstrH{$MGS} .= join("\t", $MGS, $sd3, $sampleQCStatus,
					sprintf('%.6f', $ambiguous_fraction),
					sprintf('%.6f', $csp_fraction), $locCnt)."\n";
				#push(@{$OAstrH{$MGS}},join("",@OAstr));push(@{$OFstrH{$MGS}},join("",@OFstr));
				#push(@{$OLstrH{$MGS}}, join("",@OLstr));push(@{$OCstrH{$MGS}},join("",@OCstr));
				$SInum ++ ;
				$foundGene+=$locCnt;
			}
			#clenup tmp
		} #loop over MGS
		#print "Time D: " . timeNice(time - $sttime)  . "\n";
		remove_tree($locSpace) if -d $locSpace;

		my @genesPmgs = sort { $a <=> $b } values %locMGSgenes;
		histoMGS(\@genesPmgs,"Detected Bin Genes");
		my $unaccountedMGS = $MGScnt - $SInum - $MGStoolowGskip;
		limitedWarn(q{unaccounted per-sample MGS},
			"$sd3 has $unaccountedMGS candidate MGS without a terminal outcome\n")
			if $unaccountedMGS;
		$sampleStats->{status} = q{processed};
		$sampleStats->{consensus_proteins} = scalar(keys %FAA);
		$sampleStats->{used_mgs} = $SInum;
		$sampleStats->{skipped_mgs} = $MGStoolowGskip;
		$sampleStats->{unaccounted_mgs} = $unaccountedMGS;
		$sampleStats->{used_fraction} = $MGScnt ? sprintf(q{%.6f}, $SInum / $MGScnt) : 0;
		$sampleStats->{capped_mgs} = $cappedMGS;
		$sampleStats->{capped_loci} = $cappedLoci;
		$sampleStats->{skipped_within_2_loci_of_min} = $nearThresholdMGS;
		$sampleStats->{retained_loci} = $foundGene;
		$sampleStats->{median_loci_per_used_mgs} = @genesPmgs ? median(@genesPmgs) : 0;
		$sampleStats->{mean_loci_per_used_mgs} = @genesPmgs ? sprintf(q{%.3f}, mean(@genesPmgs)) : 0;
		$sampleStats->{used_mgs_loci_histogram} = encode_loci_histogram(
			\@genesPmgs, $MGStoolowGsThr
		);
		$sampleStats->{pre_abundance_loci} = $preAbundanceLoci;
		$sampleStats->{post_abundance_loci} = $postAbundanceLoci;
		$sampleStats->{missing_consensus_loci} = $missGene;
		$sampleStats->{low_depth_loci} = $missAbundance;
		$sampleStats->{breakpoint_loci} = $breakpointFail;
		$sampleStats->{csp_rejected_loci} = $conspGen;
		$sampleStats->{ambiguous_loci} = $doubleGenes;
		$sampleStats->{abundance_filtered_loci} = $abundFail;
		$sampleStats->{invalid_protein_loci} = $SNPresFail;
		$sampleStats->{placement_flagged_mgs} = $placementFlaggedMGS;
		$sampleStats->{skip_no_selected_loci} = $skipNoSelected;
		$sampleStats->{skip_no_usable_loci} = $skipNoUsable;
		$sampleStats->{skip_too_few_after_abundance} = $skipAfterAbundance;
		$sampleStats->{skip_too_few_valid_sequences} = $skipAfterSequence;
		writeSampleStats($sampleStatsFH, $sampleStats, $sampleStatsSeen);
		if ($bufferedSamplesRef) {
			${$bufferedSamplesRef}++;
			if (${$bufferedSamplesRef} >= $appendWriteTrigger) {
				appendWriteMGSgenes($writeLink);
				${$bufferedSamplesRef} = 0;
			}
		}
	}
}

sub usage {
	my ($default) = @_;
	return <<"USAGE";
Usage: strain_within.pl -GCd DIR -MGS FILE [options]

Workflow splitting:
  -maxSubJob INT                Stage-I worker count: -1 selects automatically
                                 (default; 50-150 assembly groups/worker), 0
                                 disables splitting, positive values are explicit
  -subjob INT                   Internal split-worker index; supplied only by
                                 the parent process
  -treeOOMMaxMemGB FLOAT        Maximum memory for automatic tree OOM retries;
                                 each OOM round doubles the previous request and
                                 at most three rounds are attempted [default 1500]

Gene selection and biological QC:
  -maxGenes INT                  Maximum validated loci per MGS/sample
                                 [default $default->{maximum_genes_per_sample}]
  -noGeneLimit 0|1              Remove the locus cap; QC remains active [default 0]
  -treeLocusBudget INT           Maximum final loci selected for each tree;
                                 candidate alignment remains bounded to this plus
                                 the QC-backfill allowance
                                 [default $default->{maximum_tree_loci}]
  -disableQC 0|1                Disable biological QC (expert/debug only) [default 0]
  -MGSminGenesPSmpl INT         Minimum validated loci retained per MGS/sample
                                 [default $default->{minimum_mgs_genes_per_sample}]
  -multiGeneSmplMax FLOAT       Maximum ambiguous-locus fraction
                                 [default $default->{multi_gene_sample_max}]
  -conspGeneSmplMax FLOAT       Maximum conspecific-signal locus fraction
                                 [default $default->{conspecific_gene_sample_max}]
  -minBadLociPSmpl INT          Minimum bad loci before deferring a sample
                                 [default $default->{minimum_bad_loci_for_sample_skip}]
  -mosaicLoci FILE              Confirmed catalogue-wide mosaic/outgroup table
  -mosaicMGS FILE               Raw SB.clusters assignment table used to discover
                                 mosaics; inferred by removing .core from -MGS
  -MGSabundance FILE            Explicit MGS abundance matrix; recommended when
                                 the MGS guide is outside its Bin_* directory
  -prepareMosaicLoci 0|1        If the catalogue is absent, create it beside the
                                 MGS file in a submitted prerequisite job, wait,
                                 and validate it before extraction
                                 [default $default->{prepare_mosaic_loci}]
  -mosaicMemGb INT              Total memory requested for that prerequisite job
                                 [default 150]
  -breakpointGeneFlank INT      Mask genes this many bases around mapping breakpoints
                                 [default $default->{breakpoint_gene_flank}]
  -abundanceMinLoci INT         Loci required for robust abundance-pattern filtering
                                 [default $default->{abundance_minimum_loci}]
  -abundanceMinFold FLOAT       Lowest accepted locus/median depth ratio
                                 [default $default->{abundance_minimum_fold}]
  -abundanceMaxFold FLOAT       Highest accepted locus/median depth ratio
                                 [default $default->{abundance_maximum_fold}]
  -abundanceMaxModifiedZ FLOAT  Modified-Z threshold for depth outliers
                                 [default $default->{abundance_maximum_modified_z}]

Questionable sample-locus observations are masked first. Samples with excessive
ambiguity among the remaining observations are retained for post-tree placement
instead of being used to infer the strict backbone.

Tree locus filtering:
  -GeneLengthMin FLOAT          Minimum fraction of a locus length-Q90 retained
                                 [default 0.3]
  -GenesPerSpecies FLOAT        Backbone minimum relative locus coverage per sample
                                 [default 0.2]
  -relativeNTFraction FLOAT     Backbone minimum relative informative-NT coverage
                                 [default 0.1]
  -NTfiltCount INT              Backbone minimum informative NT after final MSA
                                 [default 0]
  -placementGenesPerSpecies FLOAT  Placement gene fraction [default 0.04]
  -placementRelativeNTFraction FLOAT  Placement NT fraction [default 0.03]
  -placementNTfiltCount INT     Placement minimum informative NT; defaults to
                                 -NTfiltCount when omitted
  -taxonAwareLocusSelection 0|1 Align a robust-plus-backfill candidate set, then
                                 select robust/core and taxon-rescue loci after MSA QC
                                 [default 1]
  -taxonAwareRescueMinPrevalence FLOAT  Minimum fraction of usable taxa carrying
                                 a locus before taxon rescue/QC backfill may select it
                                 [default 0.8]
  -rateMergePartitions 0|1      Merge final loci into deterministic rate/GC bins
                                 before IQ-TREE [default 1]
  -preferredCoreGenes FILE       Prefer universal-core seed loci listed in this
                                 raw .core guide. When omitted, use -MGS itself
                                 if it ends in .core, otherwise a readable
                                 sibling -MGS.core file [default auto]
  -compactTaxonAwareDiagnostics 0|1  Merge final taxon-aware/rate audit TSVs
                                 into phylo/taxon_aware_diagnostics.tsv
                                 [default 1]
  -rateMergeMaxBins INT         Maximum deterministic partition bins [default 8]
  -rateMergeTargetSites INT     Target effective called sites per initial bin
                                 [default 30000]
  -rateMergeMinLoci INT         Minimum loci per bin before nearest-bin merging
                                 [default 20]
  -rateMergeMinSites INT        Minimum alignment sites per bin before merging
                                 [default 20000]
  -strictBackbone 0|1           Infer a broad ML backbone and place only deferred
                                 sparse samples with EPA-ng [default 1]
  -strictBackboneFraction FLOAT Defer a sample only below this fraction of the
                                 informative-site Q90 [default 0.35]
  -strictBackboneMinSamples INT  Minimum retained backbone samples before using
                                 the complete alignment as fallback [default 3]
  -placementMinOverlap INT      Minimum informative positions shared with the
                                 inferred backbone [default 10000]
  -epaThreads INT                Requested EPA-ng threads; BuildTree caps these by
                                 cores and 1 thread/GB planning memory [default 2]
  -epaMaxMemMB INT               EPA-ng thread-planning budget; -1 derives 60% of
                                 each IQ-TREE allowance, 0 disables memory scaling
                                 [default -1]
  -epaPendantOutlierFactor FLOAT Exclude placements whose pendant branch exceeds
                                 this multiple of backbone terminal-branch Q95;
                                 zero disables the filter [default 5]
  -epaPendantMinThreshold FLOAT  Minimum pendant-branch cutoff, substitutions/site
                                 [default 0.02]
  -redoEPAfilter                Rebuild each final EPA-placed tree from its retained
                                 jplace and backbone, then continue through the normal
                                 controller validation and downstream strain analysis

A tree-only resume (-onlySubmit 1) follows the regular controller flow: it
reuses complete published or staged inputs, submits unfinished tree jobs, waits
for their outcomes, validates the result, and then submits strain_within_2.2.pl.

On a tree-only resume (-onlySubmit 1), a placementPending.sto accompanied by a
validated retained IQ-TREE backbone, MSA, query alignment, and sample
classification is retried automatically in EPA-only mode. The retry uses one
core and a doubled scheduler-memory request (minimum 20 GB), without rerunning
alignment or tree inference. Legacy runs are handled as well: when these
retained placement inputs exist but phylo/IQtree_allsites.treefile does not,
strain_within reconstructs the pending marker and restarts BuildTree in the
same isolated EPA-only mode.
USAGE
}

sub taxonAwareLocusBudgets {
	my ($geneBudget) = @_;
	die "Taxon-aware locus selection requires a positive strain gene budget\n"
		unless defined($geneBudget) && $geneBudget > 0;
	my $maximumLoci = int($geneBudget);
	my $coreLoci = int($maximumLoci * 0.8 + 0.5);
	$coreLoci = 1 if $coreLoci < 1;
	$coreLoci = $maximumLoci if $coreLoci > $maximumLoci;
	my $candidateExtra = int($maximumLoci * 0.3 + 0.5);
	return ($maximumLoci, $coreLoci, $candidateExtra);
}

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
	my $exitStatus = $?;
	if ($exitStatus != 0) {
		writeStrainWorkflowFailure($@ || 'non-zero process exit');
	} elsif (length($workflowStatePath) && $workflowStatus eq 'running') {
		$workflowStatus = 'completed';
		writeStrainWorkflowState();
	}
	my $fatalError = $exitStatus != 0 ? $@ : "";
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
	if (defined($fatalError) && length($fatalError)) {
		$fatalError =~ s/\s+$//;
		print STDERR "\nFATAL: strain_within.pl terminated: $fatalError\n";
		# The explicit final diagnostic above follows every shutdown summary.
		# Clear the active exception so Perl does not print it again out of order.
		$@ = "";
	} elsif (length($completionMessage)) {
		print "\nFINISH: $completionMessage\n";
	}
	$? = $exitStatus;
}
