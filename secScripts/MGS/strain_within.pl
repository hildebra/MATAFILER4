#!/usr/bin/perl
#This script collects representative (consensus SNP called) DNA sequences from different metagenomic samples for each MGS, and submits buildTree5.pl for each to reconstruct a phylogeny
# check performance: /ei/projects/8/88e80936-2a5d-4f4a-afab-6f74b374c765/data/geneCats/famDrama7/Bin_SB/intra4_28Feb_01D2SV/MGS.10

use warnings;
use strict;

use Getopt::Long qw( GetOptions );
use File::Path qw(make_path remove_tree);
use File::Glob qw(bsd_glob);
use File::Basename qw(basename dirname);
use File::Copy qw(copy);
use File::Spec;
use File::Temp qw(tempdir tempfile);
use Cwd qw(abs_path getcwd);
use Digest::SHA qw(sha256_hex);



use Mods::GenoMetaAss qw(gzipopen fileGZe fileGZs resolveExistingFile readClstrRev
	writeClstrRevBinaryShards readClstrRevBinaryShard
	writeSequenceBinaryCache readSequenceBinaryCache
	systemW median mean readMapS readFasta getAssemblPath getAssemblGFF getAssemblContigs);
use Mods::Subm qw(qsubSystem emptyQsubOpt qsubSystemJobAlive qsubSystemWaitMaxJobs
	deferredSubmissionDependency);
use Mods::IO_Tamoc_progs qw(getProgPaths truePath);
use Mods::FlagReference qw(printFlagHelp resolvePairedOptionDefault);
use Mods::TamocFunc qw(checkMF);
use Mods::geneCat qw(readGene2tax createGene2MGS);
use Mods::math qw(quantileArray);
use Mods::MGSLocus qw(build_locus_groups choose_locus_candidate member_context_map
	accumulate_locus_context merge_candidate_seeds preselect_locus_records
	protein_kmer_similarity robust_depth_mask);
use Mods::MosaicLoci qw(read_mosaic_catalogue);
use Mods::StrainQC qw(breakpoint_gene_mask abundance_pattern_mask);
use Mods::StrainParts qw(
	balance_assembly_groups choose_auto_worker_count choose_tree_core_count exact_worker_parts write_split_generation write_worker_completion
	split_generation_complete clear_split_generation
);
use Mods::SlurmAccounting qw(
	slurm_tree_memory_summary format_slurm_tree_memory_summary
	slurm_oom_retry_plan slurm_job_id_from_dependency
);
use Mods::WorkflowResilience qw(
	retry_operation retry_unlink retry_rename retry_open retry_close
	atomic_write_text write_workflow_record acquire_workflow_lock
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
sub phase1SamplesByGroup; sub phase1EstimatedInputBytes; sub phase1SampleWorkEstimate;
sub phase1GroupWorkEstimates; sub writePhase1WorkerPlan;
sub phase1SelectedGeneFingerprint;
sub phase1IndexShardFingerprint; sub phase1IndexShardCacheState;
sub publishPhase1IndexShards; sub loadPhase1ClusterIndex;
sub phase1ProteinCacheFingerprint; sub phase1ProteinCacheState;
sub publishPhase1ProteinCache; sub loadPhase1CatalogProteins;
sub phase1LocusModelFingerprint; sub phase1LocusModelState;
sub publishPhase1LocusModel; sub loadPhase1LocusModel; sub buildSelectedLocusGroups;
sub catalogueLocusContext;
sub mergeConspecificLogs;
sub timeNice;
sub stageStart;
#sub combineMGSgenes;
sub combineMGSgenesDir; sub prepareMGSInputSet; sub collectMGSShardHandoff;
sub writeMGSShardManifest; sub readSplitGeneration;
sub splitWorkerPartsRemain; sub getInputSize;
sub persistentMGSInputState; sub scratchMGSInputState;
sub invalidateMGSInputState;
sub msaOnlyArtifactsReady;
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
sub loadTreeOutgroupCandidates;
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
sub readPreferredCoreGeneSet;
sub treePreferredCoreGuide;
sub dispatchPendingTreeJobs;
sub retryOOMTreeJobs;
sub oomScanSeconds;
sub phase1RetryBudget;
sub submitPhase1Worker;
sub escalatePhase1WorkerOOM;
sub waitPhase1WorkersWithOOMScan;
sub phase1AutoMemoryMB;
sub phase1DefaultWorkerMemoryMB;
sub writeRecoveryRow;
sub mergeRecoveryLogs;
sub indexRecoveryRow;
sub writeRecoveryContributionIndex;
sub loadRecoveryContributionIndex;
sub writeStrainSummary;
sub writeSelectionAttritionSummary;
sub writeGeneLengthSampleSummary;
sub writeMGSSampleHistograms;
sub mergeSampleStats;
sub reportSavedSampleStats;
sub printSampleStatsSummary;
sub recoverCompletedSplitPhaseI;
sub taxonAwareLocusBudgets;
sub phase1WorkersNeedingRetry;
sub phase1WorkerCommand;
sub retryPhase1Workers;
sub strainOutputHasDurablePhaseIState;
sub phase1PathStatComponent;
sub phase1GuideStatFingerprint;
sub phase1CatalogStatFingerprint;
sub phase1InputContractContents;
sub phase1InputContractState;
sub persistPhase1InputContract;
sub writePhase1RepairQueue;
sub validatePhase1WorkerLedger;
sub writeStrainWorkflowState;
sub cleanupLegacyStrainWorkflowStateFiles;
sub writeStrainWorkflowHeartbeat;
sub writeStrainWorkflowFailure;
sub writeTreeFailureAudit;
sub lifecycleMarkerReason;
sub fastRemoveTree;
#Cached once: without a usable rm the recursive unlink must run in-process.
my $systemRemoveAvailable;

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
#.81: expose EPA-ng strict-backbone placement controls
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
#.98: disable EPA-ng strain placement by default; retain explicit opt-in
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
#1.18: balance Phase-I workers by consensus input size and regeneration cost
#1.20: stream Phase-II outgroup references directly from the catalogue
#1.21: load only core-first exact outgroup-reference demands during Phase II
#1.22: rank Mosaic outgroup proposals authoritatively against the source phylogeny
#1.23: stream parent only-submit resumes without a controller-wide file audit
#1.24: expose reproducible PopGenStats configuration to the postprocessing launcher
#1.25: recover partial loci after high-threshold sample QC and report both length gates
#1.26: consolidate destructive recovery under -redo and deprecate redundant flags
#1.27: expose EPA placement as -placeOnBackbone and fully gate placement-only filters
#1.28: prevent new outputs from entering lean tree-only resume without Phase-I evidence
#1.29: support alignment-only BuildTree runs without tree-dependent postprocessing
#1.30: bind Phase-I outputs and split workers to the selected SNP caller
#1.31: bound Stage-I output buffering and report durable shard publication
#1.32: report locus-model scan, protein, grouping, and worker projection timings
#1.33: fan the catalogue cluster index into atomic binary split-worker shards
#1.34: publish one validated selected-catalogue protein cache for all Phase-I workers
#1.35: canonicalize unlimited split-worker extraction as -maxGenes 0
#1.36: make Phase-I contracts portable across compute-node device namespaces
#1.37: escalate accounting-confirmed Phase-I OOM retries and increase tree scratch headroom
#1.38: guide BuildTree core requests by the submitted sample count
#1.39: preselect one viable outgroup per MGS before loading exact reference loci
#1.40: unify Phase-II core, memory, and submission priority around prepared job size
#1.41: validate MSA-only completion from per-locus alignments instead of a concatenation
#1.42: build the Phase-I locus model once in the parent and publish it to split workers
#1.43: exclude mixed-strain samples from the de novo tree, independent of placement
#1.44: make the coverage thresholds the primary sample filter and split the per-MGS locus floor
#1.45: stream the catalogue-wide locus-model scan so the parent no longer holds every sample
#1.46: scale tree memory with thread count and let the OOM ceiling, not the round count, stop retries
#1.47: restrict the catalogue-wide locus-model scan to contigs that carry a merge candidate
#1.48: park large trees with one rename and unlink them with a backgrounded rm
#1.49: bound Stage-I output buffers by bytes and release post-model state in workers
#1.50: give sample retention an absolute informative-NT floor, hand buildTree5 the
#	ranked extraction guide, and rename the sample QC verdict to name the finding
#	(single_strain/mixed_strain) rather than a downstream disposition
#1.51: rescan Slurm accounting on a fixed interval so OOM jobs are resubmitted
#	while the rest of their wave still runs, with the retry budget counted per
#	job instead of per wave
#1.52: let OOM retries overtake the queued bulk wave, by handicapping ordinary
#	submissions with --nice, dispatching retries from a leading tier, and
#	optionally capping how much of the wave is queued at once
#1.53: make population genetics opt-in and forward -popGenCategory, the map
#	columns that group its per-sample statistics
#1.54: keep every run-derived MGS guide product inside the output directory, so
#	concurrent runs over one catalogue cannot share or delete each other's files
#1.55: publish the catalogue context actually computed by the parent, pre-budget
#	its focal loci without dropping lower-ranked neighbours, and report true scan time
#1.56: serialize parent controllers per output so an accidental concurrent redo
#	cannot park/delete a live run's staged guide, worker scratch, or published state
#1.57: let either explicitly supplied gene-count/NT-fraction coverage threshold
#	provide the otherwise implicit paired threshold, including the placement pair
#1.58: restore hard locus-prevalence selection and use one 40% locus-length
#	threshold for both locus QC and MSA inclusion by default
#1.59: lower the default unresolved multigene-locus fraction to 0.10
#1.60: enable MSAfix isolated within-locus sequence-outlier masking by default
my $version = 1.60;


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
my $redoMode = "none";
my %deprecatedOptionSeen;
my $treeFile = "";
my $doSubmit=0;
my $subMode="";
my %FILTER_DEFAULT = (
	multi_gene_sample_max => 0.10,
	conspecific_gene_sample_max => 0.05,
	minimum_gene_depth => 1,
	minimum_bad_loci_for_sample_skip => 3,
	minimum_mgs_genes_per_sample => 8,
	minimum_loci_per_mgs => 8,
	maximum_genes_per_sample => 600,
	maximum_tree_loci => 400,
	#Absolute floor under the -relativeNTFraction/-GenesPerSpecies gate, which
	#is otherwise purely relative to the cohort's own 0.9 quantile: without it a
	#uniformly shallow cohort keeps every sample while the same sample is dropped
	#from a cohort that happens to contain a few deep ones. 5000 informative NT is
	#roughly five well-covered loci, below -MGSminGenesPSmpl's 8-locus extraction
	#minimum, so it only removes samples that carry too little sequence to place
	#at all and never overrides the relative gate.
	minimum_informative_nt_per_sample => 5000,
	breakpoint_gene_flank => 50,
	abundance_minimum_loci => 8,
	abundance_minimum_fold => 1 / 3,
	abundance_maximum_fold => 3,
	abundance_maximum_modified_z => 3.5,
	prepare_mosaic_loci => 1,
	phase1_flush_samples => 50,
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
my $onlyMSA = 0; #retain localized per-locus alignments; skip combined-MSA work and trees
my $phyloProg = 1; #1=IQ-TREE, 2=VeryFastTree, 3=FastTree
my $iqPathogen = 0; #opt in to IQ-TREE 3 pathogen/CMAPLE mode
my $GenesPerSpecies = 0.2;
my $GeneLengthMin = 0.4;
my $GeneLengthIncludeMin = $GeneLengthMin;
my $geneLengthIncludeMinSpecified = 0;
my $relativeNTFraction = 0.1;
my ($genesPerSpeciesSpecified, $relativeNTFractionSpecified) = (0, 0);
my $NTfiltCount = $FILTER_DEFAULT{minimum_informative_nt_per_sample};
my ($placementGenesPerSpecies, $placementRelativeNTFraction, $placementNTfiltCount);
$placementGenesPerSpecies = 0.04; $placementRelativeNTFraction = 0.03;
my ($placementGenesPerSpeciesSpecified,
	$placementRelativeNTFractionSpecified) = (0, 0);
my @pairedDefaultInheritances;
my $taxonAwareLocusSelection = 0;
my $taxonAwareRescueMinPrevalence = 0.8;
my $outgroupCoreMinLoci = 0; # derive as 20% of -treeLocusBudget unless overridden
my $preferredCoreGenes = "";
my $compactTaxonAwareDiagnostics = 1;
my $rateMergePartitions = 1;
my $outgroupReferenceGeneCap = 2500; # candidate reference genes retained per outgroup MGS
my $rateMergeMaxBins = 8;
my $rateMergeTargetSites = 30_000;
my $rateMergeMinLoci = 20;
my $rateMergeMinSites = 20_000;
my $postAlignmentSequenceOutlierMask = 1;
my $strictBackbone = 0;
#Extraction QC flags a sample as mixed when too many of its loci hold
#unresolvable same-COG paralogs (-multiGeneSmplMax) or conspecific consensus
#calls (-conspGeneSmplMax). Those samples describe a strain mixture, not one
#strain, so they are dropped from the tree. Sparse-but-clean samples are a
#separate decision and stay governed by the tree-inclusion filters.
my $excludeMixedStrainSamples = 1;
#-GenesPerSpecies/-relativeNTFraction/-NTfiltCount are this workflow's primary
#sample-inclusion policy. buildTree5 offers applying them as a removal only as a
#generic, default-off mechanism, because it also builds trees for callers with
#no such policy; the strain workflow is the caller that switches it on.
my $enforceSampleCoverage = 1;
my ($placeOnBackboneSpecified, $legacyStrictBackboneSpecified) = (0, 0);
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
my $selfMemAuto = 0; #-selfMemGb auto (or -1): size Phase-I workers from the input
#Phase-I worker memory model, consulted only in auto mode. A split worker holds
#its own shard of the catalogue cluster index, the shared locus model published
#by the parent, and a per-sample working set. The coefficients are deliberately
#generous: an underestimate is recovered by the Slurm OOM escalation already in
#retryPhase1Workers, but only after a wasted run, whereas an overestimate merely
#queues a larger allocation. Tune these if the reported estimate is off.
my $phase1MemBaseMB = 2048;	#interpreter, code and fixed structures
my $phase1MemIndexFactor = 3.0;	#in-memory expansion over raw cluster-index bytes
my $phase1MemLocusKB = 2;	#per resolved locus; every worker loads the whole model
my $phase1MemSampleMB = 8;	#per owned assembly group
my $phase1MemFloorMB = 4096;	#never request less than this per worker
my $phase1AutoMemoryCacheMB = 0;
my $mosaicMemGb = 150;
my $phase1WorkerRetries = 2;
# Default is read from maxMF4mem in MATAFILERcfg.txt after option parsing.
# Retain a built-in fallback for existing installations with an older config.
my $treeOOMMaxMemGB = 512;
my $treeOOMMaxMemGBSpecified = 0;
#Doubling from the initial request must be able to reach -treeOOMMaxMemGB:
#with the 5 GB floor that is seven rounds to 512 GB. A lower bound made the
#round count, rather than the configured ceiling, decide when to give up.
my $treeOOMRetryRounds = 8;
#Both phases submit their biggest, most OOM-prone jobs first and leave a long
#tail of short ones behind them. Waiting for that whole wave before the first
#accounting scan delayed every escalation until the tail had drained, so instead
#poll Slurm accounting on this cadence and resubmit whatever already died while
#the rest of the wave keeps running.
my $oomScanMinutes = 60;
#The retry budget is therefore counted per job, not per wave: a shared round
#counter would let an early scan spend the retries a job that only fails much
#later still needs. Every accounting-confirmed OOM job keeps at least this many
#escalations of its own, whichever scan first sees it fail.
my $oomMinRetries = 3;
#Rescanning alone cannot help a retry that is queued behind thousands of jobs
#the same run submitted earlier: Slurm orders pending jobs by priority, and an
#OOM retry is always the youngest. Slurm subtracts --nice from priority and lets
#any user raise it, so ordinary jobs carry this handicap while every OOM retry
#submits at nice 0 and therefore outranks the whole backlog. Set 0 to disable.
my $jobNice = 5000;
#Optional ceiling on this user's live (running + pending) scheduler jobs. With a
#cap the wave is submitted in batches, so an OOM retry only has to overtake the
#jobs already queued, not every job the run will ever submit. 0 keeps the
#historical behaviour of handing the whole wave to the scheduler at once.
my $maxQueuedJobs = 0;
#While a retained queue is still draining, top up scheduler capacity on this
#cadence. The accounting rescan keeps its own, much slower -oomScanMinutes.
my $treeQueueDrainProbeSeconds = 60;
my $nextCapacityNotice = 0;
#IQ-TREE allocates per-thread likelihood buffers, so a 60-thread job needs far
#more than the same alignment at 4 threads. The base multiplier is calibrated at
#four cores; this divisor turns the requested core count into that factor. Raise
#it to weaken the scaling, or set it above -maxCores to disable it.
my $treeMemThreadDivisor = 4;
my $redoSubmissionData = 0;
my $deepRepair = 0;
my $rmMSA = 1; #remove per-locus MSAs unless a downstream analysis requires them
my $doPopGenStats = 0; #opt-in: PopGenStats is much slower than strainStats
my $popGenStrictOutgroup = 0;
my $popGenGeneticCode = 1;
my $popGenCodonStart = 1;
my $popGenSeed = 1;
my $popGenLegacyTextOutput = 0;
my $popGenCategory = ""; #map columns grouping the population statistics
my $contTests = ""; my $discTests = ""; #stat tests to be given to strain_within_2.2.pl
my $familyVar = ""; my $groupStabilityVars = "";
my $individualVar = "AssmblGrps";

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
#Cheap extraction-time prefilter: a sample's loci for one MGS.  It is not the
#sample-inclusion policy - that is -GenesPerSpecies/-relativeNTFraction/
#-NTfiltCount at tree time, which measure selected, aligned, informative loci
#relative to the cohort.  This floor only avoids writing and aligning records
#that cannot clear any of them.
my $MGStoolowGsThr = $FILTER_DEFAULT{minimum_mgs_genes_per_sample};
#Distinct loci an MGS needs before a tree is worth building at all.  This is a
#property of the MGS, not of any sample, and was previously conflated with the
#per-sample floor above.
my $minLociPerMGS = $FILTER_DEFAULT{minimum_loci_per_mgs};
my $mode = "MGS";
my $appendWriteTrigger = $FILTER_DEFAULT{phase1_flush_samples};
#A sample count cannot bound the Stage-I output buffers, which hold one string
#per MGS this worker has touched. Flush on accumulated bytes as well so peak
#memory stays flat no matter how many MGS or loci a sample contributes.
my $flushOutputMB = 2048;
my $flushOutputByteLimit = $flushOutputMB * 1024 * 1024;
my $bufferedOutputBytes = 0;
# Publish Stage-I records often enough to bound Perl's large string buffers.
# The worker shards themselves remain on durable shared scratch because the
# parent consumes them after the worker's node-local storage has been removed.
my $startSubFromMGS = ""; #debug option: only start resubmitting tree building from this MGS (e.g. "MGS.1382" )
#define local files..
my $lSNPdir="SNP"; my $lMAPdir = "mapping";
my $SNPcaller = "MPI";
my ($lConsFNA, $lConsCTG, $lConsFAA, $lConsVCF, $lConsVCFsup,
	$phase1CatalogIdentity, $phase1CatalogInputFingerprint,
	$phase1MGSGuideFingerprint, $phase1MapSpec) = ("") x 9;


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
my $phase1InputContractName = "strainPhaseI.input.tsv";
my $phase1InputContractVersion = 3;
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
	"tmpD=s"         => \$locTmpDir1,
	"nodeTmp=s"      => sub { $locTmpDir1 = $_[1]; $deprecatedOptionSeen{nodeTmp} = '-tmpD'; },
	"submit=i"       => \$doSubmit,
	"selfMemGb=s"    => sub {
		my (undef, $value) = @_;
		$value = '' unless defined $value;
		if (lc($value) eq 'auto' || $value eq '-1') {
			$selfMemAuto = 1;
		} elsif ($value =~ /^\d+(?:\.\d+)?$/ && $value > 0) {
			$selfMemAuto = 0;
			$selfMemGb = $value;
		} else {
			die "-selfMemGb must be a positive number of GB, or \"auto\"\n";
		}
	},
	"mosaicMemGb=i"  => \$mosaicMemGb,
	"phase1WorkerRetries=i" => \$phase1WorkerRetries,
	"treeOOMMaxMemGB=f" => sub {
		$treeOOMMaxMemGB = $_[1];
		$treeOOMMaxMemGBSpecified = 1;
	},
	"treeOOMRetryRounds=i" => \$treeOOMRetryRounds,
	"oomScanMinutes=f" => \$oomScanMinutes,
	"oomMinRetries=i" => \$oomMinRetries,
	"jobNice=i" => \$jobNice,
	"maxQueuedJobs=i" => \$maxQueuedJobs,
	"treeMemThreadDivisor=f" => \$treeMemThreadDivisor,
	"onlySubmit=i"   => \$onlySubmit, #submit only jobs, or also recreate input fna/faa files? (can take days)
	"redo=s"         => \$redoMode,
	"reSubmit=i"     => sub { $reSubmit = $_[1]; $deprecatedOptionSeen{reSubmit} = '-redo tree'; },
	"recalcTrees=i"  => sub { $recalcTrees = $_[1]; $deprecatedOptionSeen{recalcTrees} = '-redo tree'; },
	"repairCAT=i"    => sub { $repairCAT = $_[1]; $deprecatedOptionSeen{repairCAT} = '-redo input'; },
	"deepRepair=i"   => sub { $deepRepair = $_[1]; $deprecatedOptionSeen{deepRepair} = '-redo input'; },
	"redoSubmissionData=i" => sub { $redoSubmissionData = $_[1]; $deprecatedOptionSeen{redoSubmissionData} = '-redo all'; },
	#workflow HPC usage
	"subjob=i"       => \$subJob,
	"maxSubJob=i"    => \$maxSubJob,
	"treeSubFromMGS=s" => sub { $startSubFromMGS = $_[1]; $deprecatedOptionSeen{treeSubFromMGS} = '-MGSsubset'; },
	#"cores=i"        => \$numCores, #not used any longer..
	"maxCores=i"     => \$maxCores, # supersedes -cores; when positive, scale tree cores from the submitted sample count
	"presortGenes=i" => \$presortGenes, #how many potential genes to include, of the original MGS (receovered will vary strongly  between samples)
	"maxGenes=i"     => \$maxNGenes, #maximum validated genes retained for each MGS/sample
	"treeLocusBudget=i" => \$treeLocusBudget, #bounded final loci passed to BuildTree selection
	"noGeneLimit=i"  => sub { $noGeneLimit = $_[1]; $deprecatedOptionSeen{noGeneLimit} = '-maxGenes 0'; },
	"disableQC=i"    => \$disableQC, #expert/debug option: disable biological QC independently of the gene cap
	"mosaicLoci=s"   => \$mosaicLociFile, #catalogue-wide confirmed mosaic/outgroup table
	"mosaicMGS=s"    => \$mosaicMGSFile, #raw SB.clusters used for comprehensive Mosaic discovery
	"MGSabundance=s" => \$MGSabundanceOverride, #explicit MGS abundance matrix for nonstandard guide locations
	"prepareMosaicLoci=i" => \$prepareMosaicLoci, #create the default catalogue if absent
	"flushEvery=i"   => \$appendWriteTrigger, #samples buffered before per-MGS records are flushed
	"flushMemMB=i"   => \$flushOutputMB, #buffered Stage-I output before an early flush
	
	"forceSNPcalls=i"  => \$forceVCF2FNA,
	"preCompConsSNP=i"   => \$preCompCons,
	"MGSsubset=s"    => \$subsMGSstr,
	"submissionMode=s"      => \$subMode,
	"MGset=s"        => \$useGTDBmg,
	
	#used genes fine tuning..
	"MGSminGenesPSmpl=i" => \$MGStoolowGsThr, #extraction prefilter: loci a sample needs for one MGS
	"minLociPerMGS=i" => \$minLociPerMGS, #distinct loci an MGS needs before a tree is attempted
	"multiGeneSmplMax=f" => \$multiGeneSmplMax, #default 0.10
	"conspGeneSmplMax=f" => \$conspGeneSmplMax, #default 0.05
	"minBadLociPSmpl=i" => \$minBadLociForSampleSkip,
	"breakpointGeneFlank=i" => \$breakpointGeneFlank,
	"abundanceMinLoci=i" => \$abundanceMinimumLoci,
	"abundanceMinFold=f" => \$abundanceMinimumFold,
	"abundanceMaxFold=f" => \$abundanceMaximumFold,
	"abundanceMaxModifiedZ=f" => \$abundanceMaximumModifiedZ,
	
	#transferred to buildTRee script..
	"GenesPerSpecies=f" => sub {
		$GenesPerSpecies = $_[1];
		$genesPerSpeciesSpecified = 1;
	},
	"GeneLengthMin=f" => \$GeneLengthMin,
	"GeneLengthIncludeMin=f" => sub {
		$GeneLengthIncludeMin = $_[1];
		$geneLengthIncludeMinSpecified = 1;
	},
	"relativeNTFraction=f" => sub {
		$relativeNTFraction = $_[1];
		$relativeNTFractionSpecified = 1;
	},
	"NTfiltCount=i" => \$NTfiltCount,
	"placementGenesPerSpecies=f" => sub {
		$placementGenesPerSpecies = $_[1];
		$placementGenesPerSpeciesSpecified = 1;
	},
	"placementRelativeNTFraction=f" => sub {
		$placementRelativeNTFraction = $_[1];
		$placementRelativeNTFractionSpecified = 1;
	},
	"placementNTfiltCount=i" => \$placementNTfiltCount,
	"preferredCoreGenes=s" => \$preferredCoreGenes,
	"compactTaxonAwareDiagnostics=i" => \$compactTaxonAwareDiagnostics,
	"taxonAwareLocusSelection=i" => \$taxonAwareLocusSelection,
	"taxonAwareRescueMinPrevalence=f" => \$taxonAwareRescueMinPrevalence,
	"outgroupCoreMinLoci=i" => \$outgroupCoreMinLoci,
	"rateMergePartitions=i" => \$rateMergePartitions,
	"rateMergeMaxBins=i" => \$rateMergeMaxBins,
	"outgroupReferenceGeneCap=i" => \$outgroupReferenceGeneCap,
	"rateMergeTargetSites=i" => \$rateMergeTargetSites,
	"rateMergeMinLoci=i" => \$rateMergeMinLoci,
	"rateMergeMinSites=i" => \$rateMergeMinSites,
	"postAlignmentSequenceOutlierMask=i" => \$postAlignmentSequenceOutlierMask,
	"placeOnBackbone=i" => sub { $strictBackbone = $_[1]; $placeOnBackboneSpecified = 1; },
	"excludeMixedStrainSamples=i" => \$excludeMixedStrainSamples,
	"enforceSampleCoverage=i" => \$enforceSampleCoverage,
	"strictBackbone=i" => sub {
		$strictBackbone = $_[1];
		$legacyStrictBackboneSpecified = 1;
		$deprecatedOptionSeen{strictBackbone} = '-placeOnBackbone';
	},
	"strictBackboneFraction=f" => \$strictBackboneFraction,
	"strictBackboneMinSamples=i" => \$strictBackboneMinSamples,
	"placementMinOverlap=i" => \$placementMinOverlap,
	"epaThreads=i" => \$epaThreads,
	"epaMaxMemMB=i" => \$epaMaxMemMB,
	"epaPendantOutlierFactor=f" => \$epaPendantOutlierFactor,
	"epaPendantMinThreshold=f" => \$epaPendantMinThreshold,
	"redoEPAfilter:i" => sub { $redoEPAfilter = $_[1] || 1; },
	"MSAprog=i"      => \$MSAprog, #2=MAFFT, 4=muscle5
	"onlyMSA=i"      => \$onlyMSA,
	"phyloProg=i"    => \$phyloProg, #1=IQ-TREE, 2=VeryFastTree, 3=FastTree
	"iqPathogen=i"   => \$iqPathogen, #explicitly enable IQ-TREE 3 pathogen/CMAPLE mode
	"rmMSA=i"        => \$rmMSA, #remove MSA, to save diskspace
	"popGenStats=i"  => \$doPopGenStats, #requires retained per-locus nucleotide MSAs
	"popGenStrictOutgroup=i" => \$popGenStrictOutgroup,
	"popGenGeneticCode=i" => \$popGenGeneticCode,
	"popGenCodonStart=i" => \$popGenCodonStart,
	"popGenSeed=i" => \$popGenSeed,
	"popGenLegacyTextOutput=i" => \$popGenLegacyTextOutput,
	"popGenCategory=s" => \$popGenCategory, #forwarded to popGenStats.R --category
	"phyloMemMulti=f" => \$memMulti, #mem used for buildtree. Default: 1.0
	
	"MGSphylo=s"     => \$treeFile,
	#transferred to MG-STK
	"ContTests=s"      => \$contTests, #continous stat tests to be handed to next step (just a passthrough)
	"DiscTests=s"      => \$discTests, #discrete stat tests to be handed to next step (just a passthrough)
	"familyVar=s"      => \$familyVar, #column name in metadata containing family id
	"groupStabilityVars=s"      => \$groupStabilityVars, #column names of categories used for calculation of resilience and persistence
	"individualVar=s" => \$individualVar, #sample ID column used by population-genetics and treeWAS analyses
	
	#SNP calling
	"SNPcaller=s"    => \$SNPcaller,
	"minSNPDepth=i"  => \$minSNPDepth,
	"minSNPCallQual=i"  => \$minSNPCallQual,
	"skipIndels=i"     => \$noIndels,
	"SNPadaptiveQual=f" => \$useAdaptiveQual, #Default 0 (not active, recommended 0.15-0.5
	"SNPdepthFilterScale=f" => \$depthFilterScale, #Default 0.15
	"SNPindelRangeFilt=i" => \$indelRange,
	"help|h" => \$help,

) or die "Invalid strain_within.pl options\n";
if ($help) {
	#option tables and their defaults come from docs/flag_reference.md
	printFlagHelp(
		script  => "strain_within.pl",
		version => $version,
		usage   => ["strain_within.pl -GCd DIR -MGS FILE [options]",
			"strain_within.pl -help | -h"],
		summary => "Within-MGS locus extraction, quality control, tree preparation/submission "
			."and hand-off to strain_within_2.2.pl. Normally launched by MGS.pl.",
		exit    => 1,
	);
}
$GeneLengthIncludeMin = $GeneLengthMin unless $geneLengthIncludeMinSpecified;
for my $pair (
	[
		{ name => 'GenesPerSpecies', value_ref => \$GenesPerSpecies,
			specified => $genesPerSpeciesSpecified },
		{ name => 'relativeNTFraction', value_ref => \$relativeNTFraction,
			specified => $relativeNTFractionSpecified },
	],
	[
		{ name => 'placementGenesPerSpecies', value_ref => \$placementGenesPerSpecies,
			specified => $placementGenesPerSpeciesSpecified },
		{ name => 'placementRelativeNTFraction',
			value_ref => \$placementRelativeNTFraction,
			specified => $placementRelativeNTFractionSpecified },
	],
) {
	my $inheritance = resolvePairedOptionDefault(
		first => $pair->[0], second => $pair->[1]);
	push @pairedDefaultInheritances, $inheritance if $inheritance;
}
die "-placeOnBackbone cannot be combined with deprecated -strictBackbone\n"
	if $placeOnBackboneSpecified && $legacyStrictBackboneSpecified;
for my $name (sort keys %deprecatedOptionSeen) {
	warn "Option -$name is deprecated; use $deprecatedOptionSeen{$name} instead. "
		."Compatibility support will be removed in a future release.\n";
}

my %validRedoMode = map { $_ => 1 } qw(none tree input all);
die "-redo must be one of: none, tree, input, all\n"
	unless $validRedoMode{$redoMode};
my %legacyRedoModes;
$legacyRedoModes{tree} = 1 if $reSubmit || $recalcTrees;
$legacyRedoModes{input} = 1 if $repairCAT || $deepRepair;
$legacyRedoModes{all} = 1 if $redoSubmissionData;
my @legacyRedoModes = sort keys %legacyRedoModes;
die "-redo cannot be combined with deprecated redo/repair flags\n"
	if $redoMode ne 'none' && @legacyRedoModes;
die "Deprecated redo/repair flags request conflicting modes: @legacyRedoModes\n"
	if @legacyRedoModes > 1;
$redoMode = $legacyRedoModes[0] if $redoMode eq 'none' && @legacyRedoModes;
$reSubmit = $recalcTrees = $repairCAT = $deepRepair = $redoSubmissionData = 0;
if ($redoMode eq 'tree') { $recalcTrees = 1; $onlySubmit = 1; }
elsif ($redoMode eq 'input') { $repairCAT = 1; $deepRepair = 1; $onlySubmit = 1; }
elsif ($redoMode eq 'all') { $redoSubmissionData = 1; $onlySubmit = 0; }
die "-redo all cannot be combined with -MGSsubset because the shared output "
	."would be cleared before only the subset was rebuilt; use -redo tree/input "
	."for a subset, or omit -MGSsubset for a full rebuild\n"
	if $redoMode eq 'all' && length($subsMGSstr);

if (!$treeOOMMaxMemGBSpecified) {
	my $configuredMaxMF4mem = getProgPaths("maxMF4mem", 0);
	if (defined($configuredMaxMF4mem) && $configuredMaxMF4mem =~ /^([0-9]+(?:\.[0-9]+)?)$/ && $1 > 0) {
		$treeOOMMaxMemGB = $1 + 0;
	} elsif (defined($configuredMaxMF4mem) && length($configuredMaxMF4mem)) {
		warn "Ignoring invalid maxMF4mem setting '$configuredMaxMF4mem'; using 512 GiB\n";
	}
}
die "Unexpected positional arguments: @ARGV\n" if @ARGV;
die "-SNPcaller must be MPI or FB\n"
	unless $SNPcaller eq 'MPI' || $SNPcaller eq 'FB';
$lConsFNA = "genes.shrtHD.SNPc.${SNPcaller}.fna.gz";
$lConsCTG = "contig.SNPc.${SNPcaller}.fna.gz";
$lConsFAA = "proteins.shrtHD.SNPc.${SNPcaller}.faa.gz";
$lConsVCF = "allSNP.${SNPcaller}.vcf.gz";
$lConsVCFsup = "allSNP.${SNPcaller}-sup.vcf.gz";
die "Tree for outgroup specified, but file is missing or empty: $treeFile\nAborting..\n"
	if $treeFile ne "" && !-s $treeFile;
die "-redoEPAfilter must be 0 or 1\n"
	unless $redoEPAfilter == 0 || $redoEPAfilter == 1;
die "-onlyMSA must be 0 or 1\n"
	unless $onlyMSA == 0 || $onlyMSA == 1;
die "-onlyMSA 1 cannot be combined with -placeOnBackbone 1\n"
	if $onlyMSA && $strictBackbone;
die "-onlyMSA 1 cannot be combined with -redo tree or -redoEPAfilter\n"
	if $onlyMSA && ($redoMode eq 'tree' || $redoEPAfilter);
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
die "-minLociPerMGS must be positive\n" unless $minLociPerMGS > 0;
die "-treeLocusBudget must be positive\n" unless $treeLocusBudget > 0;
$outgroupCoreMinLoci = int($treeLocusBudget * 0.20 + 0.999999)
	if $outgroupCoreMinLoci == 0;
die "-outgroupCoreMinLoci must be positive\n" unless $outgroupCoreMinLoci > 0;
die "-flushEvery must be positive\n" unless $appendWriteTrigger > 0;
die "-flushMemMB must be positive\n" unless $flushOutputMB > 0;
$flushOutputByteLimit = $flushOutputMB * 1024 * 1024;
die "-phase1WorkerRetries must be between 0 and 10\n"
	unless $phase1WorkerRetries >= 0 && $phase1WorkerRetries <= 10;
die "-outgroupReferenceGeneCap must be positive\n" unless $outgroupReferenceGeneCap > 0;
die "-treeOOMMaxMemGB must be positive\n" unless $treeOOMMaxMemGB > 0;
die "-treeOOMRetryRounds must be between 0 and 12
"
	unless $treeOOMRetryRounds >= 0 && $treeOOMRetryRounds <= 12;
die "-treeMemThreadDivisor must be positive
" unless $treeMemThreadDivisor > 0;
die "-oomScanMinutes must be positive
" unless $oomScanMinutes > 0;
die "-oomMinRetries must be between 0 and 12
"
	unless $oomMinRetries >= 0 && $oomMinRetries <= 12;
die "-jobNice must not be negative; Slurm only lets an operator raise priority
" unless $jobNice >= 0;
die "-maxQueuedJobs must not be negative
" unless $maxQueuedJobs >= 0;
die "Fractional filtering options must be between 0 and 1\n"
	if grep { $_ < 0 || $_ > 1 } ($multiGeneSmplMax, $conspGeneSmplMax,
		$GenesPerSpecies, $GeneLengthMin, $GeneLengthIncludeMin,
		$relativeNTFraction,
		$taxonAwareRescueMinPrevalence,
		$strictBackbone
			? grep { defined } ($placementGenesPerSpecies, $placementRelativeNTFraction)
			: ());
die "-GeneLengthIncludeMin cannot exceed -GeneLengthMin because the inclusion "
	."threshold must not be stricter than the QC threshold\n"
	if $GeneLengthIncludeMin > $GeneLengthMin;
die "-NTfiltCount must be non-negative\n" if $NTfiltCount < 0;
die "-placementNTfiltCount must be non-negative\n"
	if $strictBackbone && defined($placementNTfiltCount) && $placementNTfiltCount < 0;
die "-compactTaxonAwareDiagnostics must be 0 or 1\n"
	unless $compactTaxonAwareDiagnostics == 0 || $compactTaxonAwareDiagnostics == 1;
die "-taxonAwareLocusSelection must be 0 or 1\n"
	unless $taxonAwareLocusSelection == 0 || $taxonAwareLocusSelection == 1;
die "-rateMergePartitions must be 0 or 1\n"
	unless $rateMergePartitions == 0 || $rateMergePartitions == 1;
die "-postAlignmentSequenceOutlierMask must be 0 or 1\n"
	unless $postAlignmentSequenceOutlierMask == 0
		|| $postAlignmentSequenceOutlierMask == 1;
die "-rateMergeMaxBins, -rateMergeTargetSites, -rateMergeMinLoci, and -rateMergeMinSites must be positive\n"
	if grep { $_ < 1 } ($rateMergeMaxBins, $rateMergeTargetSites, $rateMergeMinLoci, $rateMergeMinSites);
die "-placeOnBackbone must be 0 or 1\n"
	unless $strictBackbone == 0 || $strictBackbone == 1;
if ($strictBackbone) {
	die "-strictBackboneFraction must be between 0 and 1\n"
		if $strictBackboneFraction < 0 || $strictBackboneFraction > 1;
	die "-strictBackboneMinSamples must be at least 3\n"
		if $strictBackboneMinSamples < 3;
	die "-placementMinOverlap must be non-negative\n" if $placementMinOverlap < 0;
	die "-epaThreads must be positive\n" if $epaThreads < 1;
	die "-epaMaxMemMB must be -1 (derived), 0 (no memory-based scaling), or positive\n"
		if $epaMaxMemMB < -1;
	die "-epaPendantOutlierFactor and -epaPendantMinThreshold must be non-negative\n"
		if $epaPendantOutlierFactor < 0 || $epaPendantMinThreshold < 0;
}
if ($redoEPAfilter) {
	die "-redoEPAfilter requires -placeOnBackbone 1 and -phyloProg 1\n"
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
die "-redo tree must be launched by the main strainWithin process, not a split worker\n"
	if $recalcTrees && $subJob;
die "-MSAprog must be 0, 1, 2, or 4\n"
	unless grep { $MSAprog == $_ } (0, 1, 2, 4);
die "-rmMSA must be 0 or 1\n" unless $rmMSA == 0 || $rmMSA == 1;
die "-excludeMixedStrainSamples must be 0 or 1\n"
	unless $excludeMixedStrainSamples == 0 || $excludeMixedStrainSamples == 1;
die "-enforceSampleCoverage must be 0 or 1\n"
	unless $enforceSampleCoverage == 0 || $enforceSampleCoverage == 1;
die "-popGenStats must be 0 or 1\n"
	unless $doPopGenStats == 0 || $doPopGenStats == 1;
die "-popGenStrictOutgroup and -popGenLegacyTextOutput must be 0 or 1\n"
	unless ($popGenStrictOutgroup == 0 || $popGenStrictOutgroup == 1)
		&& ($popGenLegacyTextOutput == 0 || $popGenLegacyTextOutput == 1);
die "-popGenGeneticCode must be positive, -popGenCodonStart must be 1, 2, or 3, and -popGenSeed must be non-negative\n"
	unless $popGenGeneticCode > 0 && $popGenCodonStart >= 1 && $popGenCodonStart <= 3
		&& $popGenSeed >= 0;
die "-popGenCategory requires -popGenStats 1\n" if length($popGenCategory) && !$doPopGenStats;
die "-individualVar must not be empty\n" unless length($individualVar);
if (!$subJob && $doPopGenStats && $rmMSA) {
	warn "Population genetics requires per-locus nucleotide MSAs; overriding -rmMSA 1 to -rmMSA 0\n";
	$rmMSA = 0;
}

$GCd = abs_path($GCd);
$GCd .= "/" unless $GCd =~ m{/$};
$MGSfile = abs_path($MGSfile) if length $MGSfile;
$phase1CatalogIdentity = catalog_identity($GCd);
$phase1MapSpec = resolve_catalog_maps($GCd);
$phase1CatalogInputFingerprint = phase1CatalogStatFingerprint(
	$GCd, $clusterID, $useGTDBmg, $phase1MapSpec);
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
my $preferredCoreGeneSet = length($preferredCoreGenes)
	? readPreferredCoreGeneSet($preferredCoreGenes) : {};
$outDpre = File::Spec->rel2abs($outDpre) if length $outDpre;
$mosaicLociFile = File::Spec->rel2abs($mosaicLociFile) if length $mosaicLociFile;
$MGSabundanceOverride = File::Spec->rel2abs($MGSabundanceOverride)
	if length $MGSabundanceOverride;

if (!length($mosaicMGSFile) && length($MGSfile)) {
	$mosaicMGSFile = $MGSfile;
	$mosaicMGSFile =~ s/\.core\z//;
}

$noGeneLimit = 1 if $maxNGenes <= 0; #backward-compatible no-cap spelling; QC is unchanged
die "-maxGenes must be at least -MGSminGenesPSmpl, or <=0 to remove the cap\n"
	if !$noGeneLimit && $maxNGenes < $MGStoolowGsThr;
# Every sample walks the same shared priority order and is truncated at
# -maxGenes, so a cap below -presortGenes is a prefix cut that bites only the
# well-covered samples: they stop before the tail of the pool, while shallow
# samples never reach the cap and keep sporadic low-priority loci. Loci ranked
# beyond the cap therefore show artificially low prevalence among exactly the
# samples that recovered the most, which is the opposite of what the tree-time
# selector reads prevalence to mean.
# The bias exists whenever the cap is below the pool, but it only reaches the
# tree when the locus budget can select into the affected range. Keep the
# always-true structural fact as a reported number (below, with the locus
# hierarchy) and warn only about the configuration that actually breaks:
# a budget no single sample can satisfy, which leaves every sample short of
# full coverage and makes the Q90-relative retention gate unreachable.
warn "-treeLocusBudget $treeLocusBudget exceeds -maxGenes $maxNGenes: no sample "
	."can contribute more than $maxNGenes loci, so no sample can reach the "
	."selected-locus count the tree budget implies, and loci ranked beyond "
	."position $maxNGenes carry prevalence biased downward against exactly the "
	."best-covered samples. Raise -maxGenes to at least -treeLocusBudget, or "
	."lower -treeLocusBudget.\n"
	if !$noGeneLimit && $maxNGenes > 0 && $treeLocusBudget > $maxNGenes;
# Keep one canonical unlimited value after parsing. In particular, split-worker
# commands must not carry both the deprecated -noGeneLimit alias and a negative
# -maxGenes value: -maxGenes 0 is the documented no-cap contract.
$maxNGenes = 0 if $noGeneLimit;

$onlySubmit = 1 if $recalcTrees; #tree-only recovery reuses published or complete staged inputs
# Ordinary parent -onlySubmit runs are latency-sensitive cluster dispatchers,
# but only after a durable Phase-I summary proves this is an actual resume.
my $leanOnlySubmitRequested = $onlySubmit && !$subJob && !$recalcTrees
	&& !$redoSubmissionData && !$repairCAT && !$deepRepair
	&& !$redoEPAfilter && !$reSubmit;
printEarlyRunHeader();

@subsetMGS = split /,/,$subsMGSstr if ($subsMGSstr ne "");
#print "SUBSMGS:: @subsetMGS\n";
#die timeNice(20) ." ".timeNice(12252)."\n"; #TEST

#define global vars
my $queueMode = $subMode;
$queueMode = "bash" if !$doSubmit && $queueMode eq "";
my $QSBoptHR = emptyQsubOpt($doSubmit,"",$queueMode);
#Ordinary submissions carry the priority handicap and respect the queue ceiling;
#OOM recovery overrides the handicap per job when it is dispatched.
$QSBoptHR->{jobNice} = $jobNice;
$QSBoptHR->{maxConcurrentJobs} = $maxQueuedJobs;
my $MGSfileOri = $MGSfile; #save for later..


my $resumeBindir = $MGSfile;
$resumeBindir =~ s/[^\/]+$//;
$resumeBindir = $GCd if $resumeBindir eq "";
my $resumeOutD = length($outDpre) ? $outDpre : "$resumeBindir/intra_phylo/";
# The output itself is deliberately replaceable during a full rebuild, so its
# lock cannot live underneath that tree: fastRemoveTree would rename it away and
# a second controller could lock a new inode at the recreated path. Keep one
# stable sibling lock open for the lifetime of every parent controller. Split
# workers are part of that controller's generation and must not contend for it.
my $parentRunLock;
if (!$subJob) {
	my $lockBase = File::Spec->canonpath(File::Spec->rel2abs($resumeOutD));
	$lockBase =~ s{[\\/]+\z}{};
	my $parentRunLockPath = "$lockBase.strain_within.lock";
	my $host = $ENV{HOSTNAME} // $ENV{HOST} // 'unknown';
	my $job = $ENV{SLURM_JOB_ID} // $ENV{JOB_ID} // $ENV{LSB_JOBID} // 'none';
	$parentRunLock = acquire_workflow_lock(
		$parentRunLockPath,
		label => "strain_within parent for $resumeOutD",
		owner => join(' ', "pid=$$", "host=$host", "job=$job",
			'started='.time, "redo=$redoMode", "onlySubmit=$onlySubmit"),
	);
	print "Acquired exclusive strain parent lock: $parentRunLockPath\n";
}
# Contract validation happens before prepRun() assigns $outD. Fingerprint the
# same run-local sorted guide that prepRun() and the parent use after staging;
# otherwise workers compare the catalogue-side .srt against the parent's staged
# .srt and reject an otherwise identical SNP-input contract.
$phase1MGSGuideFingerprint =
	phase1GuideStatFingerprint($MGSfileOri, undef, $resumeOutD);
my $resumePhaseISummary = File::Spec->catfile(
	$resumeOutD, 'LOGandSUB', $sampleStatsSummaryLogName);
my $resumePhaseIInputContract = File::Spec->catfile(
	$resumeOutD, 'LOGandSUB', $phase1InputContractName);
my $unsafeSubsetRebuild = !$onlySubmit && !$subJob && length($subsMGSstr)
	&& strainOutputHasDurablePhaseIState(
		$resumeOutD, $resumePhaseISummary, $resumePhaseIInputContract,
		File::Spec->catfile($resumeOutD, $sampleStatsSummaryLogName),
	);
if ($unsafeSubsetRebuild) {
	die "Refusing an input-rebuilding run (-onlySubmit 0) with -MGSsubset in "
		."existing strain output $resumeOutD because shared non-subset results "
		."would be cleared. Use a fresh -outD, a non-destructive resume mode, "
		."or rebuild without -MGSsubset.\n";
}
my ($phase1ContractState, $phase1ContractReason,
	$recordedPhase1ContractVersion, $recordedPhase1ContractStatus) =
	phase1InputContractState($resumePhaseIInputContract);
my $legacyMPIContract = $phase1ContractState eq 'missing' && $SNPcaller eq 'MPI';
if ($onlySubmit && !$subJob
		&& $phase1ContractState ne 'match'
		&& $phase1ContractState ne 'building_match'
		&& !$legacyMPIContract) {
	die "Cannot reuse Phase-I outputs with -SNPcaller $SNPcaller: "
		."$phase1ContractReason. Run a complete rebuild with -redo all "
		."(without -MGSsubset) before resuming.\n";
}
if ($onlySubmit && $subJob
		&& $phase1ContractState ne 'match'
		&& $phase1ContractState ne 'building_match'
		&& !$legacyMPIContract) {
	die "Split worker Phase-I input contract is incompatible with this run "
		."(-SNPcaller $SNPcaller): "
		."$phase1ContractReason\n";
}
if ($onlySubmit && !$subJob && $phase1ContractState eq 'building_match') {
	print STDERR "Resuming compatible incomplete Phase-I state: "
		."$phase1ContractReason\n";
}
if ($onlySubmit && !$subJob && $legacyMPIContract) {
	print STDERR "No Phase-I SNP-input contract found at $resumePhaseIInputContract; "
		."treating existing outputs as legacy MPI state. A successful resume will "
		."record the explicit contract.\n";
}
my $leanOnlySubmitResume = $leanOnlySubmitRequested && -s $resumePhaseISummary;
if ($leanOnlySubmitRequested && !$leanOnlySubmitResume) {
	print STDERR "Lean only-submit resume unavailable: no durable Phase-I "
		."completion summary at $resumePhaseISummary; auditing inputs and "
		."running extraction where required.\n";
}

my ($preparedMainBranchFastPath, @preparedMainBranchMGS) = (0);
my %preparedMainBranchCategoryValidated;
if (!$leanOnlySubmitResume && $onlySubmit && !$subJob && length($MGSfile)
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
			fastRemoveTree($mosaicRunDirectory);
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
$phase1MGSGuideFingerprint = phase1GuideStatFingerprint($MGSfileOri);
my $activePhase1InputContract =
	File::Spec->catfile($LOGDIR, $phase1InputContractName);
if (!$subJob && !$onlySubmit) {
	persistPhase1InputContract($activePhase1InputContract, 'building');
} elsif (!$subJob && defined($recordedPhase1ContractVersion)
		&& $recordedPhase1ContractVersion < $phase1InputContractVersion
		&& ($phase1ContractState eq 'match'
			|| $phase1ContractState eq 'building_match')) {
	my ($preparedState, $preparedReason) =
		phase1InputContractState($activePhase1InputContract);
	die "Cannot upgrade legacy Phase-I input contract before worker dispatch: "
		."$preparedReason\n"
		unless $preparedState eq 'match'
			|| $preparedState eq 'building_match';
	persistPhase1InputContract(
		$activePhase1InputContract, $recordedPhase1ContractStatus);
}
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
my ($phase1Groups) = phase1SamplesByGroup();
my $groupedSampleCount = 0;
$groupedSampleCount += scalar(@{$AGlist{$_}}) for keys %AGlist;
my $effectiveGroupCount = scalar(keys %{$phase1Groups});
my $standaloneSampleCount = scalar(@samples) - $groupedSampleCount;
if ($maxSubJob == -1) {
	my ($automaticWorkers, $targetGroupsPerWorker) = choose_auto_worker_count(
		$effectiveGroupCount, scalar(@samples),
	);
	$maxSubJob = $automaticWorkers;
	print "Automatic Stage-I splitting: $effectiveGroupCount effective groups "
		."(".scalar(keys %AGlist)." shared, $standaloneSampleCount standalone), "
		.scalar(@samples)." samples, target ${targetGroupsPerWorker} groups/worker; "
		.($maxSubJob ? "using $maxSubJob workers" : "using the main process only")."\n";
}
stepComplete("assembly-group expansion", $stepStarted,
	"groups=".scalar(keys %AGlist), "effective_groups=$effectiveGroupCount",
	"grouped_samples=$groupedSampleCount", "standalone_samples=$standaloneSampleCount");
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
#Published cluster-index shards, recorded so the catalogue-wide locus-model scan
#can stream them one at a time instead of holding every sample at once.
my $phase1ShardPaths; my $phase1ShardFingerprint = "";
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
my ($dirsNOTPrepped, $CatFileMiss, $CatNotPrepped, $treeAbsent, $doneDirs,
	$PhylosExist, $noRecoverableLociDirs, $completedTreeFastPaths);
if ($leanOnlySubmitResume) {
	# -onlySubmit is an explicit request to dispatch trees from an existing Phase-I
	# handoff. Populate paths in memory and defer the few useful checks to the MGS
	# immediately being prepared, so cluster work can begin without an all-MGS
	# metadata barrier.
	$SIdirs{$_} = "$outD/$_/" for @specis;
	($dirsNOTPrepped, $CatFileMiss, $CatNotPrepped, $treeAbsent, $doneDirs,
		$PhylosExist, $noRecoverableLociDirs, $completedTreeFastPaths)
		= (0, scalar(@specis), scalar(@specis), scalar(@specis), 0, 0, 0, 0);
	stepComplete("existing-output and resume audit", $stepStarted,
		"mode=deferred_per_MGS", "selected_MGS=".scalar(@specis),
		"global_metadata_scans=0");
	print "Lean only-submit resume: Phase-I input, terminal-marker, and tree checks "
		."will run once for each MGS immediately before its submission decision.\n";
} else {
	($dirsNOTPrepped, $CatFileMiss, $CatNotPrepped, $treeAbsent, $doneDirs,
		$PhylosExist, $noRecoverableLociDirs, $completedTreeFastPaths)
		= evalFileStatus();
	my $auditEpaOnly = scalar(keys %MGSepaOnlyRetry);
	my $auditLegacyEpa = scalar(grep {
		($MGSepaOnlyRetry{$_} // '') eq 'legacy_missing_final'
	} keys %MGSepaOnlyRetry);
	my $auditFullTree = $treeAbsent - $auditEpaOnly;
	$auditFullTree = 0 if $auditFullTree < 0;
	stepComplete("existing-output and resume audit", $stepStarted,
		"prepared_trees=$doneDirs", "completion_marker_fast_paths=$completedTreeFastPaths",
		"missing_trees=$treeAbsent", "incomplete_tree_inputs=$CatFileMiss",
		"directories_needing_extraction=$dirsNOTPrepped",
		"validated_no_locus=$noRecoverableLociDirs",
		"epa_only_retries=$auditEpaOnly", "legacy_epa_retries=$auditLegacyEpa",
		"full_tree_retries=$auditFullTree");
}
my $epaOnlyRetryCount = scalar(keys %MGSepaOnlyRetry);
my $legacyEpaRetryCount = scalar(grep {
	($MGSepaOnlyRetry{$_} // '') eq 'legacy_missing_final'
} keys %MGSepaOnlyRetry);
my $fullTreeRetryCount = $treeAbsent - $epaOnlyRetryCount;
$fullTreeRetryCount = 0 if $fullTreeRetryCount < 0;
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
	
	stageStart(q{Stage I: consensus-gene extraction},
		q{Extracting relevant core MGS loci from SNP-consensus assemblies});
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
	# Retain the Phase-I locus map until Phase II has selected its outgroups.
	# It avoids a second catalogue-wide gene2tax scan before tree submission.
	
	reportingsMGS();
	%smplsPerMGS = (); #reporting-only sample/locus counts can be large
	# $SIgenes and $COGprios are reused by direct Phase-II outgroup lookup.
	
	my @jobsMain;
	my (%phase1JobByWorker, %phase1MemoryByWorker);
	#Retry bookkeeping is per worker and outlives a single wave, so the live OOM
	#supervisor and the later ledger audit share one budget per worker.
	my (%phase1RetriesByWorker, %phase1HandledJobs, %phase1TerminalWorkers);
	my $phase1InitialMemoryMB = 0;
	my $phase1SelfCmd = '';

	if ($maxSubJob && !$subJob){
		# Sized here because this is the first point at which the locus model and
		# the sample assignment are both known.
		$phase1InitialMemoryMB = phase1DefaultWorkerMemoryMB(
			workers => $maxSubJob,
			loci => scalar(keys %{$LocusByID}),
			samples => scalar(keys %cl2gene2),
		);
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
			my ($dep,$qcmd) = qsubSystem($LOGDIR."Strain1_B${sj}.sh",$cmdX,1,"${phase1InitialMemoryMB}M","Str1.$sj","","",1,[],$QSBoptHR);
			push(@jobsMain,$dep);
			my $jobID = slurm_job_id_from_dependency($dep, $QSBoptHR->{rTag});
			$phase1JobByWorker{$sj} = $jobID if defined($jobID);
			$phase1MemoryByWorker{$sj} = $phase1InitialMemoryMB;
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
		# Keep the selected locus/protein map until Phase II has staged all
		# outgroup overlays. This prevents a second full gene2tax read and keeps
		# same-COG outgroup matching sequence-aware.
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
		#Workers are balanced by assembly group, so a wave mixes one very large
		#worker with many short ones. Rescan Slurm accounting on the -oomScanMinutes
		#cadence and resubmit an OOM-killed worker at once, rather than letting the
		#tail of the wave hold its escalation back for hours.
		waitPhase1WorkersWithOOMScan(
			jobs => \@jobsMain, generation => $splitGeneration,
			worker_command => $phase1SelfCmd, script_kind => "retry",
			job_by_worker => \%phase1JobByWorker,
			memory_by_worker => \%phase1MemoryByWorker,
			retries_by_worker => \%phase1RetriesByWorker,
			handled_jobs => \%phase1HandledJobs,
			terminal_workers => \%phase1TerminalWorkers,
			default_mb => $phase1InitialMemoryMB,
		) if @jobsMain && $doSubmit;
		my @failedWorkers = phase1WorkersNeedingRetry($splitGeneration);
		if (@failedWorkers) {
			my $remaining = retryPhase1Workers(
				generation => $splitGeneration,
				workers => \@failedWorkers,
				worker_command => $phase1SelfCmd,
				script_kind => "retry",
				job_by_worker => \%phase1JobByWorker,
				memory_by_worker => \%phase1MemoryByWorker,
				retries_by_worker => \%phase1RetriesByWorker,
				handled_jobs => \%phase1HandledJobs,
				terminal_workers => \%phase1TerminalWorkers,
			);
			if (@{$remaining}) {
				my $queue = writePhase1RepairQueue($splitGeneration, $remaining,
					"live Phase-I worker validation failed");
				$completionMessage = "Phase I requires worker repair before Phase II; no tree jobs were submitted.";
				print "Phase-I processing paused safely; repair queue: $queue. Invalid workers: "
					.join(",", @{$remaining})."\n";
				exit(0);
			}
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
	
	print "Stage I hand-off: consensus-gene inputs are published and ready for phylogeny preparation. "
		."Global elapsed ".timeNice(time - $^T)."\n";

} else {
	print "Skipping Part I, all required per-MGS inputs are already prepared.\n";
	my $mergedSampleStats = recoverCompletedSplitPhaseI();
	if ($mergedSampleStats) {
		invalidateMGSInputState(@specis);
	} else {
		reportSavedSampleStats();
	}
}
persistPhase1InputContract(File::Spec->catfile($LOGDIR, $phase1InputContractName))
	unless $subJob;
loadRecoveryContributionIndex()
	unless $recoveryContributionIndexReady || $leanOnlySubmitResume;

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
my %OGgenesByCOG;
my %outgroupGeneCache;
my %TreeOutgroupCandidates;
my %SelectedOutgroup;
my %outgroupCategoryPreflight;

# Reference data is initialized lazily after EPA-only recovery jobs have been
# submitted. It reuses the Phase-I locus map and streams the catalogue FASTAs
# once, with no .fai lookup or on-disk outgroup cache.
my $requiresOutgroupReference = $runPartI || $CatNotPrepped || $repairCAT
	|| $deepRepair || $redoSubmissionData;
my %outgroupEligibleLoci;
my %outgroupCatalogueMGS;
my $TreeOutgroupCandidatesBulkLoaded = 0;
my $outgroupReferenceInitialized = 0;
my %outgroupDemandLoci;
my %outgroupDemandMinimum;
my $initializeOutgroupReferences = sub {
	my ($targetMGS) = @_;
	return if $outgroupReferenceInitialized;
	$outgroupReferenceInitialized = 1;
	my $referenceStarted = time;
	unless ($requiresOutgroupReference && @{$targetMGS || []}) {
		stepComplete("outgroup-reference preparation", $referenceStarted,
			"status=not_required", "reference_NT=0",
			"required=".($requiresOutgroupReference ? 1 : 0),
			"reference_AA=0", "MGS_with_outgroup_candidates=0");
		return;
	}
	my (%demandCogsByOutgroup, %eligibleCogsByOutgroup);
	my $candidateStarted = time;
	loadTreeOutgroupCandidates($targetMGS) if length($treeFile);
	my $candidateMGS = 0;
	my $nextCandidateProgress = time + 60;
	for my $MGS (@{$targetMGS}) {
		my @candidates;
		if (length($treeFile)) {
			# R has already resolved the Mosaic proposal against the phylogeny
			# and returned at most one authoritative outgroup.
			push @candidates, treeOutgroupCandidates($MGS);
		} elsif (exists($PreferredOutgroup{$MGS}) && length($PreferredOutgroup{$MGS})) {
			push @candidates, $PreferredOutgroup{$MGS};
		}
		my %seenCandidate;
		@candidates = grep {
			defined($_) && /\A[A-Za-z0-9][A-Za-z0-9_.:+-]*\z/ && !$seenCandidate{$_}++
		} @candidates;
		$SelectedOutgroup{$MGS} = $candidates[0] if @candidates;
		$candidateMGS++;
		if (time >= $nextCandidateProgress) {
			stepProgress("outgroup candidate discovery", $candidateMGS,
				scalar(@{$targetMGS}), $candidateStarted,
				"selected_outgroup_MGS=".scalar(keys %SelectedOutgroup));
			$nextCandidateProgress = time + 60;
		}
	}
	unless (%SelectedOutgroup) {
		stepComplete("outgroup-reference preparation", $referenceStarted,
			"status=no_selected_outgroups", "reference_NT=0", "required=1",
			"reference_AA=0", "MGS_with_outgroup_candidates=0");
		return;
	}

	my ($refFNA, $refFAA, $refNameL) = ("", "", "unknown");
	if ($mode eq "MGS" || $mode eq "MGSall") {
		$refFNA = resolveExistingFile("$GCd/compl.incompl.$clusterID.fna")
			// "$GCd/compl.incompl.$clusterID.fna";
		$refFAA = resolveExistingFile("$GCd/compl.incompl.$clusterID.prot.faa")
			// "$GCd/compl.incompl.$clusterID.prot.faa";
		$refNameL = "geneCat";
	} elsif ($mode eq "FMG") {
		$refFNA = "$GCd/FMG/COG*.fna";
		$refFAA = "$GCd/FMG/COG*.faa";
		$refNameL = "FMG ref";
	}

	my $mapStarted = time;
	my $nextMapProgress = time + 60;
	my ($hr1, $Gene2COG_OG, $hr4);
	if (keys(%{$SIgenes}) && keys(%{$Gene2COG}) && keys(%{$COGprios})) {
		($hr1, $Gene2COG_OG, $hr4) = ($SIgenes, $Gene2COG, $COGprios);
		print "Reusing the Phase-I selected-MGS gene map for core-first target-locus lookup; global elapsed "
			.timeNice(time - $^T)."\n";
	} else {
		my $unusedGene2MGS;
		($hr1, $Gene2COG_OG, $unusedGene2MGS, $hr4) = readGene2tax(
			$gene2taxF, $presortGenes, $targetMGS,
			sub {
				my ($status) = @_;
				return if time < $nextMapProgress;
				stepProgress("outgroup target gene-map loading", $status->{rows_scanned}, undef,
					$mapStarted, "included_genes=$status->{included_genes}");
				$nextMapProgress = time + 60;
			},
		);
	}

	my (%cogTaxa, $targetMapCount);
	for my $MGS (@{$targetMGS}) {
		next unless exists($hr4->{$MGS});
		$targetMapCount++;
		my %seenCOG;
		for my $locus (@{$hr4->{$MGS}}) {
			my (undef, $cog, $gene) = locusParts($locus, $MGS);
			next unless length($cog) && length($gene);
			$seenCOG{$cog} = 1;
		}
		$cogTaxa{$_}++ for keys %seenCOG;
	}
	my $broadMinimumTaxa = int($targetMapCount * $taxonAwareRescueMinPrevalence + 0.999999);
	$broadMinimumTaxa = 1 if $broadMinimumTaxa < 1;
	my %broadCOG = map { $_ => 1 } grep {
		$cogTaxa{$_} >= $broadMinimumTaxa
	} keys %cogTaxa;

	my (%requiredNT, %requiredAA, %directReferenceByMGS);
	my ($coreSatisfiedMGS, $broadFallbackMGS, $unmetDemandMGS, $demandLoci,
		$noSelectedMGS) = (0, 0, 0, 0, 0);
	for my $MGS (@{$targetMGS}) {
		my @targetLoci = @{$hr4->{$MGS} || []};
		my $selectedOutgroup = $SelectedOutgroup{$MGS} // '';
		my (%demandCOG, %seenCoreCOG, %seenBroadCOG, %seenDemandLocus);
		my @demand;
		my $addDemand = sub {
			my ($locus) = @_;
			return if !defined($locus) || !length($locus) || $seenDemandLocus{$locus}++;
			my (undef, $cog, $gene) = locusParts($locus, $MGS);
			return unless length($cog) && length($gene);
			$demandCOG{$cog} = 1;
			push @demand, $locus;
		};
		for my $locus (@targetLoci) {
			my (undef, $cog, $gene) = locusParts($locus, $MGS);
			next unless length($cog) && length($gene)
				&& exists($preferredCoreGeneSet->{$gene}) && !$seenCoreCOG{$cog}++;
			$addDemand->($locus);
		}
		for my $locus (@targetLoci) {
			my (undef, $cog, $gene) = locusParts($locus, $MGS);
			next unless length($cog) && length($gene)
				&& length($selectedOutgroup)
				&& $selectedOutgroup eq ($PreferredOutgroup{$MGS} // '')
				&& exists($PreferredOutgroupGene{$MGS}{$gene})
				&& ($broadCOG{$cog} || exists($preferredCoreGeneSet->{$gene}));
			$addDemand->($locus);
			my $directGene = $PreferredOutgroupGene{$MGS}{$gene};
			$directReferenceByMGS{$MGS}{$directGene} = 1
				if defined($directGene) && length($directGene);
		}
		my $coreAndDirectCount = scalar(@demand);
		for my $locus (@targetLoci) {
			last if @demand >= $outgroupCoreMinLoci;
			my (undef, $cog, $gene) = locusParts($locus, $MGS);
			next unless length($cog) && length($gene) && $broadCOG{$cog}
				&& !$seenBroadCOG{$cog}++;
			$addDemand->($locus);
		}
		$coreSatisfiedMGS++ if length($selectedOutgroup)
			&& $coreAndDirectCount >= $outgroupCoreMinLoci;
		$broadFallbackMGS++ if length($selectedOutgroup)
			&& $coreAndDirectCount < $outgroupCoreMinLoci
			&& @demand >= $outgroupCoreMinLoci;
		$unmetDemandMGS++ if length($selectedOutgroup)
			&& @demand < $outgroupCoreMinLoci;
		$noSelectedMGS++ unless length($selectedOutgroup);
		$outgroupDemandLoci{$MGS}{$_} = 1 for @demand;
		$outgroupDemandMinimum{$MGS} = $outgroupCoreMinLoci;
		$demandLoci += scalar(@demand);
		unless (length($selectedOutgroup) && @demand >= $outgroupCoreMinLoci) {
			delete $SelectedOutgroup{$MGS};
			next;
		}
		my %eligibleCOG = %demandCOG;
		for my $locus (@targetLoci) {
			my (undef, $cog, $gene) = locusParts($locus, $MGS);
			next unless length($cog) && length($gene)
				&& ($broadCOG{$cog} || exists($preferredCoreGeneSet->{$gene}));
			$outgroupEligibleLoci{$MGS}{$locus} = 1;
			$eligibleCOG{$cog} = 1;
		}
		$outgroupCatalogueMGS{$selectedOutgroup} = 1;
		$demandCogsByOutgroup{$selectedOutgroup}{$_} = 1 for keys %demandCOG;
		$eligibleCogsByOutgroup{$selectedOutgroup}{$_} = 1 for keys %eligibleCOG;
		$requiredNT{$_} = 1 for keys %{$directReferenceByMGS{$MGS} || {}};
		$requiredAA{$_} = 1 for keys %{$directReferenceByMGS{$MGS} || {}};
	}

	print "Outgroup core-demand manifest: target_MGS=".scalar(@{$targetMGS})
		."; core_or_direct=$coreSatisfiedMGS; broad_fallback=$broadFallbackMGS"
		."; below_minimum=$unmetDemandMGS; no_selection=$noSelectedMGS"
		."; demanded_loci=$demandLoci; viable_selected_MGS="
		.scalar(keys %SelectedOutgroup)
		."; candidate_gene_cap=$outgroupReferenceGeneCap\n";
	unless (%outgroupCatalogueMGS) {
		my $status = $unmetDemandMGS && !$noSelectedMGS
			? "all_below_minimum"
			: "no_viable_outgroups";
		stepComplete("outgroup-reference preparation", $referenceStarted,
			"status=$status", "reference_NT=0", "required=1",
			"reference_AA=0", "MGS_with_outgroup_candidates=0");
		return;
	}

	print "Preparing core-first exact outgroup-reference demands for "
		.scalar(keys %outgroupCatalogueMGS)
		." selected MGS; no FASTA index or outgroup cache; global elapsed "
		.timeNice(time - $^T)."\n";

	my @outgroupMGS = sort keys %outgroupCatalogueMGS;
	my ($candidateSIgenes, $candidateGene2COG, $candidateCOGprios);
	my $candidateMapStarted = time;
	my $nextCandidateMapProgress = time + 60;
	my $unusedCandidateGene2MGS;
	($candidateSIgenes, $candidateGene2COG, $unusedCandidateGene2MGS, $candidateCOGprios) = readGene2tax(
		$gene2taxF, $outgroupReferenceGeneCap, \@outgroupMGS,
		sub {
			my ($status) = @_;
			return if time < $nextCandidateMapProgress;
			stepProgress("outgroup candidate gene-map loading", $status->{rows_scanned}, undef,
				$candidateMapStarted, "included_genes=$status->{included_genes}");
			$nextCandidateMapProgress = time + 60;
		},
		{ allowed_cogs_by_mgs => \%eligibleCogsByOutgroup },
	);

	my @mappedMGS = sort grep { exists($outgroupCatalogueMGS{$_}) } keys %{$candidateCOGprios};
	my ($mappedCount, $mappedCandidateGenes) = (0, 0);
	my $nextMappingProgress = time + 60;
	for my $MGS (@mappedMGS) {
		my $eligibleCOG = $eligibleCogsByOutgroup{$MGS} || {};
		my $demandCOG = $demandCogsByOutgroup{$MGS} || {};
		next unless keys %{$eligibleCOG};
		my (@candidateRecords, %seenCandidateGene);
		for my $locus (@{$candidateCOGprios->{$MGS} || []}) {
			my $gene = $candidateSIgenes->{$MGS}{$locus};
			next unless defined($gene) && defined($candidateGene2COG->{$gene})
				&& !$seenCandidateGene{$gene}++;
			my $cog = $candidateGene2COG->{$gene};
			next unless $eligibleCOG->{$cog};
			push @candidateRecords, [$cog, $gene, $demandCOG->{$cog} ? 1 : 0];
		}
		my $retainedForMGS = 0;
		my (%retainedCandidateGene, %retainedDemandCOG);
		my $addCandidate = sub {
			my ($cog, $gene) = @_;
			return if $retainedForMGS >= $outgroupReferenceGeneCap
				|| $retainedCandidateGene{$gene}++;
			push @{$OGgenesByCOG{$MGS}{$cog}}, $gene;
			$requiredNT{$gene} = 1;
			$requiredAA{$gene} = 1;
			$retainedForMGS++;
		};
		# Keep one candidate per demanded COG before using the remaining pool for
		# copy choice, so a small cap cannot crowd out the acceptance evidence.
		for my $record (@candidateRecords) {
			next unless $record->[2];
			next if $retainedDemandCOG{$record->[0]}++;
			$addCandidate->($record->[0], $record->[1]);
		}
		for my $demandPass (1, 0) {
			for my $record (@candidateRecords) {
				next unless $record->[2] == $demandPass;
				$addCandidate->($record->[0], $record->[1]);
			}
		}
		$mappedCandidateGenes += $retainedForMGS;
		$mappedCount++;
		if (time >= $nextMappingProgress) {
			stepProgress("outgroup core-demand lookup", $mappedCount,
				scalar(@mappedMGS), $mapStarted,
				"requested_reference_genes=".scalar(keys %requiredNT),
				"candidate_genes=$mappedCandidateGenes");
			$nextMappingProgress = time + 60;
		}
	}
	unless (%requiredNT && %requiredAA) {
		stepComplete("outgroup-reference preparation", $referenceStarted,
			"status=no_resolvable_loci", "reference_NT=0", "required=1",
			"reference_AA=0", "MGS_with_outgroup_candidates=".scalar(keys %outgroupCatalogueMGS));
		return;
	}
	print "Outgroup exact-reference request: viable_outgroup_MGS="
		.scalar(keys %outgroupCatalogueMGS)."; candidate_genes=$mappedCandidateGenes"
		."; requested_reference_genes=".scalar(keys %requiredNT)."\n";

	my $faaStarted = time;
	my $nextFaaProgress = time + 60;
	$FAAref = readFasta($refFAA, 1, "\\s", \%requiredAA,
		sub {
			my ($status) = @_;
			return if time < $nextFaaProgress;
			stepProgress("outgroup protein FASTA streaming",
				$status->{records_scanned}, undef, $faaStarted,
				"retained=$status->{records_retained}",
				"requested=".scalar(keys %requiredAA),
				"file=".basename($status->{file}));
			$nextFaaProgress = time + 60;
		});
	my $fnaStarted = time;
	my $nextFnaProgress = time + 60;
	$FNAref = readFasta($refFNA, 1, "\\s", \%requiredNT,
		sub {
			my ($status) = @_;
			return if time < $nextFnaProgress;
			stepProgress("outgroup nucleotide FASTA streaming",
				$status->{records_scanned}, undef, $fnaStarted,
				"retained=$status->{records_retained}",
				"requested=".scalar(keys %requiredNT),
				"file=".basename($status->{file}));
			$nextFnaProgress = time + 60;
		});
	print "Loaded ".scalar(keys %{$FNAref})." nucleotide and "
		.scalar(keys %{$FAAref})." protein exact outgroup reference genes from $refNameL by sequential streaming\n";
	stepComplete("outgroup-reference preparation", $referenceStarted,
		"status=loaded", "mode=core_first_streaming",
		"demanded_loci=$demandLoci", "outgroup_core_min_loci=$outgroupCoreMinLoci",
		"candidate_gene_cap=$outgroupReferenceGeneCap",
		"reference_NT=".scalar(keys %{$FNAref}),
		"reference_AA=".scalar(keys %{$FAAref}),
		"MGS_with_outgroup_candidates=".scalar(keys %outgroupCatalogueMGS));
};


stageStart(q{Stage II: phylogeny preparation and submission},
	q{Preparing intra-strain phylogenies for }.scalar(@specis).q{ MGS});


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
my $fullTreeInputsInitialized = $leanOnlySubmitResume ? 1 : 0;
if ($leanOnlySubmitResume) {
	# A prior sizing table is a scheduling hint only. Reusing it avoids a fresh
	# all-MGS metadata pass; missing entries are sized just in time below.
	my $cachedSizingFile = "$LOGDIR/tree_input_sizing.tsv";
	my %wantedSize = map { $_ => 1 } @fullTreeCandidates;
	my $cachedSizes = 0;
	if (-s $cachedSizingFile && open(my $cachedSizing, '<', $cachedSizingFile)) {
		my $header = <$cachedSizing> // '';
		if ($header =~ /^MGS\tselected_state\tsource\testimated_uncompressed_MB\b/) {
			while (my $line = <$cachedSizing>) {
				$line =~ s/[\r\n]+\z//;
				my ($MGS, undef, undef, $megabytes) = split /\t/, $line, 5;
				next unless defined($MGS) && $wantedSize{$MGS}
					&& defined($megabytes) && $megabytes =~ /^\d+(?:\.\d+)?\z/
					&& $megabytes > 0;
				$inputSizeByMGS{$MGS} = 0 + $megabytes;
				$cachedSizes++;
			}
		}
		close $cachedSizing;
	}
	print "Lean only-submit sizing: reused $cachedSizes cached resource estimate(s); "
		."missing estimates will be read only when their MGS reaches submission.\n";
}
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
MGS_SUBMISSION:
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
		next MGS_SUBMISSION;
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
	my $msaOnlyStone = "$outD2/msaOnly.complete.tsv";
	my $msaOnlyOutput = "$outD2/MSA";
	my $terminalTreeMarker = "$outD2/noTree.sto";
	my $placementPendingMarker = "$outD2/placementPending.sto";
	my $IQtreef= "$outD2/phylo/IQtree_allsites.treefile";
	$IQtreef = "$outD2/phylo/VERYFASTTREE_allsites.nwk" if ($phyloProg == 2);
	$IQtreef = "$outD2/phylo/FASTTREE_allsites.nwk" if ($phyloProg == 3);
	my %resumeEntry;
	if ($leanOnlySubmitResume) {
		if (opendir(my $resumeDirectory, $outD2)) {
			$resumeEntry{$_} = 1 for readdir($resumeDirectory);
			closedir($resumeDirectory);
		}
		if ($onlyMSA && msaOnlyArtifactsReady($outD2)) {
			$treeDisposition{'valid MSA already present'}++;
			limitedNotice('MGS skipped with existing MSA-only result',
				"Skipping $MGS: a completed MSA-only result is already present.\n");
			next;
		}
		if (!$onlyMSA && $resumeEntry{'treeDone.sto'} && -s $IQtreef) {
			$treeDisposition{'valid tree already present'}++;
			limitedNotice('MGS skipped with existing trees',
				"Skipping $MGS: a completed tree is already present.\n");
			next;
		}
		my @terminalMarkers = (
			['tooFewSamples.sto', 'insufficient_tree_input'],
			['noRecoverableLoci.sto', 'no_recoverable_loci'],
			['noTree.sto', 'buildtree_no_usable_alignment'],
		);
		for my $terminal (@terminalMarkers) {
			next unless $resumeEntry{$terminal->[0]};
			my $marker = "$outD2/$terminal->[0]";
			my $reason = lifecycleMarkerReason($marker, $terminal->[1]);
			$treeDisposition{"valid no-tree: $reason"}++;
			limitedNotice('MGS skipped after valid no-tree classification',
				"Skipping $MGS: terminal no-tree state '$reason'.\n");
			next MGS_SUBMISSION;
		}
		if ($resumeEntry{'placementPending.sto'}) {
			my $epaState = epaOnlyRetryReady($outD2);
			if (length($epaState)) {
				$MGSepaOnlyRetry{$MGS} = $epaState;
				$epaOnlyRetry = 1;
				$epaRecovery = 1;
				$epaOnlyRetryCount++;
			}
		}
	}
	my $publishedInputsReady = !$epaOnlyRetry
		&& !exists($legacyLocusMGS{$MGS})
		&& (!$leanOnlySubmitResume || $resumeEntry{'data.log'}
			|| $resumeEntry{'data.log.gz'})
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
				"Skipping $MGS: -redo tree found neither complete published inputs nor a complete staged FNA/FAA/category set.\n");
			next;
		}
		resetMGSTreeOutputs($outD2, $MGS);
	}
	
	if (!$leanOnlySubmitResume && !$recalcTrees && !$reSubmit && !$repairCAT
			&& !$redoSubmissionData && !exists($legacyLocusMGS{$MGS})
			&& ($onlyMSA
				? msaOnlyArtifactsReady($outD2)
				: (-e $treeStone && -s $IQtreef))) {
		$treeDisposition{$onlyMSA ? 'valid MSA already present' : 'valid tree already present'}++;
		limitedNotice($onlyMSA ? 'MGS skipped with existing MSA-only result'
				: 'MGS skipped with existing trees',
			"Skipping $MGS: a valid ".($onlyMSA ? 'MSA-only result' : 'tree')." already exists.\n");
		next;
	}
	
	my $inputFNAsize = $inputSizeByMGS{$MGS} // 0;
	if ($leanOnlySubmitResume && !$epaRecovery
			&& !exists($inputSizeByMGS{$MGS})) {
		my $inputBytes = fileGZs("$outD2/$FNAstdof");
		$inputBytes ||= fileGZs("$tmpD/$FNAstdof");
		$inputFNAsize = $inputBytes > 0 ? $inputBytes / (1024 * 1024) : 1;
		$inputSizeByMGS{$MGS} = $inputFNAsize;
	}
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
	# Allow headroom for decompressed alignments, engine temporaries, and atomic publication.
	my $treeTmpGb = int(($inputFNAsize * 5 + 1023) / 1024);
	$treeTmpGb = 20 if $treeTmpGb < 20;
	$QSBoptHR->{tmpSpace} = $nodeTmpConfigured ? $treeTmpGb : 0;
	# Placement retains likelihood vectors across the reference tree and can use
	# substantially more memory than tree inference for long concatenated MSAs.
	# Reserve a distinct scheduler profile whenever EPA-ng placement is requested.
	my $placementRequested = $strictBackbone ? 1 : 0;
	my $baseMemMult = 75; $baseMemMult = 15 if ($phyloProg ==3 || $phyloProg ==2);
	$baseMemMult = 150 if $placementRequested && $baseMemMult < 150;
	my $memoryProfile = $onlyMSA ? 'MSA-only'
		: $placementRequested ? 'EPA-ng placement' : 'tree-only';
	my $minimumMemMB = ($placementRequested ? 10240 : 5000) * $memMulti;
	$minimumMemMB = 10240 if $placementRequested && $minimumMemMB < 10240;
	my $maximumMemMB = 110000 * $memMulti;
	$maximumMemMB = $minimumMemMB if $maximumMemMB < $minimumMemMB;
	my $totMem = int($inputFNAsize * $baseMemMult * $memMulti);
	$totMem = $minimumMemMB if $totMem < $minimumMemMB;
	$totMem = $maximumMemMB if $totMem > $maximumMemMB;
	my $taxonLocusInputMB = 0;
	my $workloadCells = 0;
	my $threadMemFactor = 1;
	my $memoryPlanningInputMB = $inputFNAsize;
	my $numCoreL = $numCores;
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
	

	my $bts = getProgPaths("buildTree_scr");
	my $treeFlag = $onlyMSA ? "" : "-runIQtree 1 ";
	if (!$onlyMSA && $phyloProg == 2){$treeFlag = "-runVeryFastTree 1 ";}
	if (!$onlyMSA && $phyloProg == 3){$treeFlag = "-runFastTree 1 ";}
	my $tree_sample_separator = quotemeta($SaSe);
	my $Tcmd= "$bts -fna ".shellQuote($FNAtf)." -aa ".shellQuote($FAAtf)." -smplSep ".shellQuote($tree_sample_separator)." -cats ".shellQuote($CATtf)." -outD ".shellQuote($outD2)." $treeFlag ";
	$Tcmd .= "-withinSpecies 1 -relativeNTFraction $relativeNTFraction "
		."-NTfiltPerGene $GeneLengthMin "
		."-GeneLengthIncludeMin $GeneLengthIncludeMin "
		."-GenesPerSpecies $GenesPerSpecies "
		."-NTfiltCount $NTfiltCount -iqFast 1 ";
	$Tcmd .= "-taxonAwareLocusSelection $taxonAwareLocusSelection ";
	if ($taxonAwareLocusSelection) {
		$Tcmd .= "-taxonAwareMaxLoci $taxonAwareMaxLoci "
			."-taxonAwareCoreLoci $taxonAwareCoreLoci "
			."-taxonAwareCandidateExtra $taxonAwareCandidateExtra "
			."-taxonAwareRescueMinPrevalence $taxonAwareRescueMinPrevalence ";
		# Hand buildTree5 the guide that actually drove extraction priority. Its
		# row order is the presorter's importance ranking (markers first, then
		# chimera/paralogy/prevalence-window evidence and expected informative
		# yield), and buildTree5 now scores loci on that order. $preferredCoreGenes
		# is resolved before sorting, so it still names the unranked .core table;
		# using it here would feed filterMB2core's row order in as a ranking.
		my $treeCoreGuide = treePreferredCoreGuide($preferredCoreGenes);
		$Tcmd .= "-preferredCoreGenes ".shellQuote($treeCoreGuide)." "
			if length($treeCoreGuide) && !$epaOnlyRetry;
	}
	$Tcmd .= "-compactTaxonAwareDiagnostics $compactTaxonAwareDiagnostics ";
	$Tcmd .= "-rateMergePartitions $rateMergePartitions "
		."-rateMergeMaxBins $rateMergeMaxBins "
		."-rateMergeTargetSites $rateMergeTargetSites "
		."-rateMergeMinLoci $rateMergeMinLoci "
		."-rateMergeMinSites $rateMergeMinSites ";
	$Tcmd .= "-postAlignmentSequenceOutlierMask "
		."$postAlignmentSequenceOutlierMask ";
	$Tcmd .= "-rmMSA $rmMSA -MSAprogram $MSAprog -onlyMSA $onlyMSA ";
	$Tcmd .= "-placeOnBackbone $strictBackbone ";
	# buildTree5 provides both as generic, default-off mechanisms; the strain
	# workflow is the caller that sets a policy for them.
	$Tcmd .= "-excludeFlaggedSamples $excludeMixedStrainSamples ";
	$Tcmd .= "-enforceSampleCoverage $enforceSampleCoverage ";
	if ($strictBackbone) {
		$Tcmd .= "-placementGenesPerSpecies $placementGenesPerSpecies "
			if defined $placementGenesPerSpecies;
		$Tcmd .= "-placementRelativeNTFraction $placementRelativeNTFraction "
			if defined $placementRelativeNTFraction;
		$Tcmd .= "-placementNTfiltCount $placementNTfiltCount "
			if defined $placementNTfiltCount;
		$Tcmd .= "-strictBackboneFraction $strictBackboneFraction "
			."-strictBackboneMinSamples $strictBackboneMinSamples "
			."-placementMinOverlap $placementMinOverlap "
			."-epaThreads ".($epaOnlyRetry ? 1 : $epaThreads)
			." -epaMaxMemMB $epaMaxMemMB "
			."-epaPendantOutlierFactor $epaPendantOutlierFactor "
			."-epaPendantMinThreshold $epaPendantMinThreshold ";
	}
	my $treeTmpOption = $nodeTmpConfigured
		? "-tmpSubdir ".shellQuote("strain_within/$MGS")
		: "-tmpD ".shellQuote("$scratchD/$MGS/");
	$Tcmd .= "$treeTmpOption -map ".shellQuote($mapF)." ";

	my $multiSmpl;my $ngenes; my $needsCopy = 0; my $inputReady = 0;
	#$multiSmpl counts tree tips and therefore includes the staged outgroup; it
	#stays the basis for core and memory planning. $ingroupSmpl excludes the
	#outgroup and is what decides whether a tree can carry any signal at all.
	my $ingroupSmpl;
	if ($epaOnlyRetry) {
		$inputReady = 1;
	} else {
		($multiSmpl,$ngenes,$OG,$needsCopy,$inputReady,$ingroupSmpl)=
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
	
	if (!$epaOnlyRetry && $maxCores > 0) {
		$numCoreL = choose_tree_core_count($multiSmpl, $maxCores);
	}
	unless ($epaOnlyRetry) {
		# addOutgroup2MGS already counted both dimensions while finalizing this
		# input. Treat each sample-by-gene cell as roughly one kilobase of
		# likelihood-matrix work; no extra file scan or preparation pass is needed.
		$workloadCells = $multiSmpl * $ngenes;
		$taxonLocusInputMB = $workloadCells / 1024;
		$memoryPlanningInputMB = $taxonLocusInputMB > $inputFNAsize
			? $taxonLocusInputMB : $inputFNAsize;
		# IQ-TREE keeps per-thread partial-likelihood buffers, so the same
		# alignment costs far more at 60 threads than at the four the base
		# multiplier is calibrated for. Ignoring this made wide, sample-rich MGS
		# start an order of magnitude under what they need.
		$threadMemFactor = $numCoreL / $treeMemThreadDivisor;
		$threadMemFactor = 1 if $threadMemFactor < 1;
		$totMem = int($memoryPlanningInputMB * $baseMemMult * $memMulti * $threadMemFactor);
		$totMem = $minimumMemMB if $totMem < $minimumMemMB;
		$totMem = $maximumMemMB if $totMem > $maximumMemMB;
	}
	my $iqMemMB = int($totMem * 0.9); #also supplies EPA planning-memory reporting
	$Tcmd .= "-cores $numCoreL ";
	if (!$onlyMSA && $phyloProg == 1){
		$Tcmd .= "-iqMemMB $iqMemMB ";
		$Tcmd .= "-iqPathogen 1 " if $iqPathogen;
	}

	$outgS = " -outgroup ".shellQuote($OG)." "  if ($OG ne "");
	$Tcmd .= "-sampleQC ".shellQuote("$outD2/$QCstdof")." "
		if !$epaRecovery && (fileGZe("$outD2/$QCstdof") || fileGZe("$tmpD/$QCstdof")
			|| fileGZe("$tmpD/$QCstdof.tmp") || $stagedShardHandoff{$MGS});
	$Tcmd .= "-stagedInputDir ".shellQuote($tmpD)." " if !$epaRecovery && $needsCopy;
	$Tcmd .= "-redoEPAfilter 1 " if $redoEPAfilter
		&& -s "$outD2/phylo/epa-ng/epa_result.jplace";
	$Tcmd .= "-epaOnly 1 " if $epaOnlyRetry;
	$Tcmd .= "-continue 1 ";
	$Tcmd .= "-completionMarker ".shellQuote($treeStone)." " unless $onlyMSA;
	$Tcmd .= "-terminalMarker ".shellQuote($terminalTreeMarker)." "
		."-placementPendingMarker ".shellQuote($placementPendingMarker)." ";

	if ($epaOnlyRetry) {
		print "$MGS (".($lcnt + 1)."/$Nspecis); elapsed ".timeNice(time - $sttime)
			."; outgroup ".(length($OG) ? $OG : 'none')
			."; samples n/a; genes n/a; 1 core; $totMem MB; EPA-ng placement-only retry\n";
	} elsif ($ingroupSmpl > 2 && $ngenes >= $minLociPerMGS){
		print "$MGS (".($lcnt + 1)."/$Nspecis); elapsed ".timeNice(time - $sttime)
			."; outgroup ".(length($OG) ? $OG : 'none')
			."; $ingroupSmpl ingroup samples ($multiSmpl tips); $ngenes genes; "
			."$numCoreL cores; $totMem MB; $memoryProfile\n";
	} else {
		# Gate on the ingroup alone: with two ingroup samples the only unrooted
		# topology over the resulting three tips is fixed in advance, so the tree
		# job can spend a full allocation and still terminate in no_usable_loci.
		my $reason = $ingroupSmpl <= 2 ? 'too_few_samples' : 'too_few_usable_genes';
		$treeDisposition{"valid no-tree: $reason"}++;
		limitedNotice('MGS with insufficient tree input',
			"$MGS: $reason (ingroup_samples=$ingroupSmpl, tree_tips=$multiSmpl, "
			."usable_genes=$ngenes); skipping tree construction\n");
		writeTooFewMarker($outD2, $ingroupSmpl, $ngenes, $reason);
		remove_tree($tmpD) if $needsCopy && -d $tmpD;
		$QSBoptHR->{tmpSpace} = $tmpSHDD;
		$QSBoptHR->{useLongQueue} = 0;
		next;
	}
	unlink "$outD2/tooFewSamples.sto" if -e "$outD2/tooFewSamples.sto";
	
	# PART II: retain each completely prepared job. Full trees are submitted
	# together after preparation so the current core selector can define a global
	# largest-first order; EPA-only recovery remains latency-prioritized.
	my $treeJobOrdinal = $cnt + 1;
	push @pendingTreeJobs, {
		mgs => $MGS,
		script => $epaOnlyRetry ? "$outD2/treeCmd.epa_retry.sh"
			: $onlyMSA ? "$outD2/treeCmd.msa_only.sh" : "$outD2/treeCmd.sh",
		command => $Tcmd.$outgS."\n",
		cores => $numCoreL,
		sample_count => defined($multiSmpl) ? $multiSmpl : 0,
		workload_cells => $workloadCells,
		gene_count => defined($ngenes) ? $ngenes : 0,
		memory => int($totMem)."M",
		requested_mb => int($totMem),
		priority_ordinal => $treeJobOrdinal,
		memory_planning_input_mb => $memoryPlanningInputMB,
		job_name => $epaOnlyRetry ? "EPA$treeJobOrdinal"
			: $onlyMSA ? "MSA$treeJobOrdinal" : "FT$treeJobOrdinal",
		epa_only => $epaOnlyRetry,
		msa_only => $onlyMSA,
		terminal => $terminalTreeMarker,
		placement_pending => $placementPendingMarker,
		tree => $onlyMSA ? $msaOnlyOutput : $IQtreef,
		stone => $onlyMSA ? $msaOnlyStone : $treeStone,
		tmp_space => $QSBoptHR->{tmpSpace},
		use_long_queue => $QSBoptHR->{useLongQueue},
		job_nice => $jobNice,
	};
	$QSBoptHR->{tmpSpace} =$tmpSHDD;
	$QSBoptHR->{useLongQueue} = 0;
	$cnt ++;
	$treeDisposition{$epaOnlyRetry ? 'EPA-only retry job'
		: $onlyMSA ? 'eligible MSA-only job' : 'eligible tree job'}++;
	$expectedTreeOutputs{$MGS} = [$onlyMSA ? $msaOnlyOutput : $IQtreef,
		$onlyMSA ? $msaOnlyStone : $treeStone,
		$terminalTreeMarker, $placementPendingMarker, $onlyMSA];
	if (!$doSubmit || ($epaOnlyRetry && time >= $nextQueuedTreeSubmissionProbe)) {
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
print "  staged input sets recovered for -redo tree: $recalcScratchRecovered\n"
	if $recalcTrees;
if ($doSubmit) {
	print "Tree preparation pass complete: $cnt eligible tree job(s), "
		."$submittedTreeJobs submitted so far, ".scalar(@pendingTreeJobs)
		." awaiting scheduler capacity. "
		."The following wait count reports jobs still present, not jobs omitted.\n";
	#The rest of the queue is drained by the OOM supervisor below. Draining it
	#here would put every escalation behind the whole remaining wave, which is
	#exactly the ordering problem the supervisor exists to avoid.
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
		fastRemoveTree($path);
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
#The dispatch order puts the largest, most OOM-prone trees first and leaves a
#long tail of short jobs behind them, so waiting for the whole wave before the
#first accounting scan delayed every escalation until that tail had drained.
#retryOOMTreeJobs owns the wait instead and rescans every -oomScanMinutes.
retryOOMTreeJobs(
	jobs => \@jobs,
	pending_queue => \@pendingTreeJobs,
	submitted_ref => \$submittedTreeJobs,
	accounting => \@treeJobAccounting,
	options => $QSBoptHR,
	maximum_mb => int($treeOOMMaxMemGB * 1024 + 0.5),
	maximum_rounds => $treeOOMRetryRounds,
);
die "Internal error: tree submission queue was not drained\n"
	if @pendingTreeJobs && $doSubmit;
print "Tree submission pass complete: $submittedTreeJobs eligible tree job(s) submitted; "
	.scalar(@jobs)." scheduler job ID(s) tracked.\n" if $doSubmit;
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
if ($onlyMSA) {
	$completionMessage = "strain_within.pl completed MSA-only processing; "
		."MSA_outcomes_quarantined=$incompleteTreeOutcomes. Tree inference and "
		."tree-dependent strain postprocessing were intentionally skipped.";
	print "MSA-only workflow complete. Tree inference, EPA-ng placement, and "
		."strain_within_2.2.pl were not launched.\n";
	exit(0);
}
if ($incompleteTreeOutcomes) {
	print "Tree-job outcomes remain quarantined in tree_job_outcomes.tsv, but all tree "
		."inputs are resolved; proceeding with downstream strain analysis for completed trees.\n";
}
print "\nAll done for $cnt Bins\nRun strain_within_2.pl for summary stats:\n";

my $MGSabundance = $MGSabundanceOverride ne ""
	? $MGSabundanceOverride
	: "$bindir/Annotation/Abundance/MGS.matL7.txt";
die "MGS abundance matrix is missing or empty: $MGSabundance\n" unless -s $MGSabundance;

my $strain2Scr = getProgPaths("MGS_strain2_scr");

my $nxtCmd = "$strain2Scr -GCd ".shellQuote($GCd)." -FMGdir ".shellQuote($outD)." -MGSmatrix ".shellQuote($MGSabundance)." -cores 4 -DiscTests ".shellQuote($discTests)." -ContTests ".shellQuote($contTests)." -familyVar ".shellQuote($familyVar)." -groupStabilityVars ".shellQuote($groupStabilityVars)." -individualVar ".shellQuote($individualVar)." ";
$nxtCmd .= "-MGSphylo ".shellQuote($treeFile)." " if $treeFile ne "";
$nxtCmd .= "-popGenStats $doPopGenStats -popGenStrictOutgroup $popGenStrictOutgroup "
	."-popGenGeneticCode $popGenGeneticCode -popGenCodonStart $popGenCodonStart "
	."-popGenSeed $popGenSeed -popGenLegacyTextOutput $popGenLegacyTextOutput ";
$nxtCmd .= "-popGenCategory ".shellQuote($popGenCategory)." " if $popGenCategory ne "";
$nxtCmd .= "-submit $doSubmit ";
$nxtCmd .= "-qsubSystem ".shellQuote($subMode)." " if $subMode ne "";
$nxtCmd .= "-Hcores $maxCores " if $maxCores > 0;
if ($mapF2 eq ""){$nxtCmd .= "-map ".shellQuote($mapF)." ";} else {$nxtCmd .= "-map ".shellQuote($mapF2)." ";}

$nxtCmd .= "\n";

#$GCd/MB2.clusters.ext.can.Rhcl.matL0.txt
	my ($dep,$qcmd) = qsubSystem($LOGDIR."strainAnalysis2.sh",$nxtCmd,1,"60G","2StrainSub","","",1,[],$QSBoptHR);
print "\n". $nxtCmd."\n";


#cleanup
fastRemoveTree($locTmpDir);
fastRemoveTree($preConDir) if $preCompCons;

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
	# merge.complete.tsv is published only after the complete aggregate passed
	# contributor/cardinality/order validation. In lean dispatch mode, trust that
	# commit and let the immediate category read detect any external corruption.
	return 1 if $leanOnlySubmitResume && -s "$tmpD/merge.complete.tsv";
	# Legacy shard-only runs need the compact contributor index, but load it only
	# when the first checkpoint-less MGS actually reaches submission.
	loadRecoveryContributionIndex()
		if $leanOnlySubmitResume && !$recoveryContributionIndexReady;
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

sub loadTreeOutgroupCandidates {
	my ($targetMGS) = @_;
	return 1 if $TreeOutgroupCandidatesBulkLoaded;
	return 1 unless defined($treeFile) && length($treeFile) && -e $treeFile;
	my %wanted = map { $_ => 1 } @{$targetMGS || []};
	return 1 unless keys %wanted;
	my $started = time;
	my $neiTree = getProgPaths("neighborTree");
	my $call = "$neiTree ".shellQuote($treeFile)." --all --max-candidates 1";

	# Give the one bulk R process the consolidated Mosaic choices. File::Temp
	# removes this transient run-level input; no per-MGS artefacts are created.
	my ($preferredFh, $preferredPath);
	my $preferredCount = 0;
	for my $MGS (sort keys %wanted) {
		next unless exists($PreferredOutgroup{$MGS})
			&& length($PreferredOutgroup{$MGS});
		unless ($MGS =~ /\A[A-Za-z0-9][A-Za-z0-9_.:+-]*\z/
			&& $PreferredOutgroup{$MGS} =~ /\A[A-Za-z0-9][A-Za-z0-9_.:+-]*\z/) {
			limitedWarn('invalid Mosaic outgroup identifiers',
				"Skipping invalid Mosaic outgroup preference $MGS -> $PreferredOutgroup{$MGS}\n");
			next;
		}
		if (!$preferredFh) {
			($preferredFh, $preferredPath) = tempfile(
				'strain_mosaic_outgroups.XXXXXX', TMPDIR => 1, UNLINK => 1);
		}
		print {$preferredFh} "$MGS\t$PreferredOutgroup{$MGS}\n"
			or die "Cannot write temporary Mosaic outgroup preferences $preferredPath: $!\n";
		$preferredCount++;
	}
	if ($preferredFh) {
		close $preferredFh
			or die "Cannot close temporary Mosaic outgroup preferences $preferredPath: $!\n";
		$call .= " --preferred ".shellQuote($preferredPath);
	}

	print "Discovering tree-neighbour candidates for ".scalar(keys %wanted)
		." actionable MGS in one R call, with $preferredCount Mosaic preference(s); "
		."global elapsed ".timeNice(time - $^T)."\n";
	open my $bulk, "$call |"
		or do {
			limitedWarn('bulk outgroup lookup command failures',
				"Can't start bulk tree-neighbour lookup $call; falling back to individual lookups\n");
			return 0;
		};
	my ($lines, $loaded, $malformed) = (0, 0, 0);
	my %preferenceDecisionCount;
	while (my $line = <$bulk>) {
		$lines++;
		$line =~ s/[\r\n]+\z//;
		my ($MGS, $decision, $preferred, $preferredDistance, $cutoff, $candidateText)
			= split /\t/, $line, 6;
		next unless defined($MGS) && exists($wanted{$MGS});
		unless (defined($decision) && $decision =~ /\A(?:none|accepted|same_as_target|absent_from_tree|non_finite_distance|too_close|no_eligible_neighbors|too_distant)\z/) {
			$malformed++;
			limitedWarn('malformed bulk outgroup rows',
				"Ignoring malformed tree-neighbour row for $MGS from $treeFile\n");
			next;
		}
		$preferenceDecisionCount{$decision}++;
		my %seen;
		my @candidates = grep {
			/\A[A-Za-z0-9][A-Za-z0-9_.:+-]*\z/ && !$seen{$_}++
		} split /\s+/, ($candidateText // '');
		@candidates = ($candidates[0]) if @candidates > 1;
		$TreeOutgroupCandidates{$MGS} = \@candidates;
		$loaded++;
	}
	my $bulkOk = close $bulk;
	unless ($bulkOk && !$malformed) {
		delete $TreeOutgroupCandidates{$_} for keys %wanted;
		limitedWarn('bulk outgroup lookup command failures',
			"Bulk tree-neighbour lookup failed or returned incompatible rows for $treeFile; falling back to individual lookups\n");
		return 0;
	}
	$TreeOutgroupCandidates{$_} = [] for grep {
		!exists($TreeOutgroupCandidates{$_})
	} keys %wanted;
	$TreeOutgroupCandidatesBulkLoaded = 1;
	my $decisionSummary = join(',', map {
		$_.'='.$preferenceDecisionCount{$_}
	} sort keys %preferenceDecisionCount);
	print "Loaded bulk tree-neighbour candidates for $loaded/".scalar(keys %wanted)
		." actionable MGS from $lines tree rows in ".timeNice(time - $started)
		."; Mosaic decisions: ".($decisionSummary || 'none')."\n";
	return 1;
}

sub treeOutgroupCandidates {
	my ($MGS) = @_;
	return @{$TreeOutgroupCandidates{$MGS}}
		if exists($TreeOutgroupCandidates{$MGS});
	my @candidates;
	# A failed bulk call retains the previous single-MGS lookup as a resilience
	# fallback. Successful bulk loading never starts one process per MGS. The
	# same Mosaic plausibility decision remains authoritative in fallback mode.
	if (defined($treeFile) && length($treeFile) && -e $treeFile) {
		my $neiTree = getProgPaths("neighborTree");
		my $call = "$neiTree ".shellQuote($treeFile)." ".shellQuote($MGS);
		$call .= " --preferred-tip ".shellQuote($PreferredOutgroup{$MGS})
			if exists($PreferredOutgroup{$MGS}) && length($PreferredOutgroup{$MGS});
		$call .= " --max-candidates 1";
		my $outgroup_text = `$call`;
		if ($? != 0) {
			limitedWarn('outgroup lookup command failures',
				"Can't find an authoritative outgroup order from call $call; leaving its candidate list empty\n");
		} else {
			@candidates = split /\s+/, $outgroup_text;
		}
	}
	my %seen;
	@candidates = grep {
		/\A[A-Za-z0-9][A-Za-z0-9_.:+-]*\z/ && !$seen{$_}++
	} @candidates;
	@candidates = ($candidates[0]) if @candidates > 1;
	$TreeOutgroupCandidates{$MGS} = \@candidates;
	return @candidates;
}

sub addOutgroup2MGS{
	my ($MGS,$OG,$tmpD) = @_;
	my $outD2 = $SIdirs{$MGS};
	my $shardHandoff = $stagedShardHandoff{$MGS};
	my ($publishedPrepared, $publishedOG) = preparedOutgroupLog($outD2);
	my $outputReady = $publishedPrepared
		&& fileGZe("$outD2/$FNAstdof")
		&& fileGZe("$outD2/$FAAstdof") && fileGZe("$outD2/$CATstdof");
	if ($outputReady && !$repairCAT && !$deepRepair && !$redoSubmissionData
			&& !exists($legacyLocusMGS{$MGS})) {
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
		#The published category already carries the outgroup overlay rows, so the
		#sample count read back from it is one above the ingroup count.
		my $ingroupSeen = scalar(keys %samplesSeen) - (defined($publishedOG)
			&& length($publishedOG) && $samplesSeen{$publishedOG} ? 1 : 0);
		return (scalar(keys %samplesSeen), $genesSeen, $publishedOG, 0, 1, $ingroupSeen);
	}

	# Compatibility for controller runs that had already completed the old
	# controller-side Phase II before this version was installed.
	my ($scratchPrepared, $preparedOG) = preparedOutgroupLog($tmpD);
	my $preparedScratchInput = $scratchPrepared && (
		$leanOnlySubmitResume
			? (-s "$tmpD/merge.complete.tsv" && fileGZe("$tmpD/$CATstdof"))
			: fileGZe("$tmpD/$FNAstdof") && fileGZe("$tmpD/$FAAstdof")
				&& fileGZe("$tmpD/$LINKstdof") && fileGZe("$tmpD/$CATstdof")
				&& fileGZe("$tmpD/$QCstdof") && -s "$tmpD/merge.complete.tsv"
	);
	if (!$shardHandoff && $preparedScratchInput && !$repairCAT && !$deepRepair
			&& !$redoSubmissionData && !exists($legacyLocusMGS{$MGS})) {
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
		my $ingroupSeen = scalar(keys %samplesSeen) - (defined($preparedOG)
			&& length($preparedOG) && $samplesSeen{$preparedOG} ? 1 : 0);
		return (scalar(keys %samplesSeen), $genesSeen, $preparedOG, 1, 1, $ingroupSeen);
	}

	my $rawCategory = "$tmpD/$CATstdof.tmp";
	my @rawCategorySources = $shardHandoff
		? map { $_->{parts}{category}{path} } @{$shardHandoff->{workers}}
		: ($rawCategory);
	my $stageReady = $shardHandoff ? scalar(@rawCategorySources)
		: $leanOnlySubmitResume
			? (-s "$tmpD/merge.complete.tsv" && fileGZe($rawCategory))
			: fileGZe("$tmpD/$FNAstdof") && fileGZe("$tmpD/$FAAstdof")
				&& fileGZe($rawCategory) && fileGZe("$tmpD/$QCstdof.tmp")
				&& -s "$tmpD/merge.complete.tsv";
	if (!$stageReady) {
		limitedWarn('MGS missing raw staged tree input',
			"$MGS has no complete raw staged FNA/FAA/category/QC input in $tmpD; leaving it for repair\n");
		return (0, 0, $OG, 0, 0, 0);
	}

	# Published or already-overlaid scratch inputs returned above without touching
	# the catalogue. The first genuinely raw MGS initializes the shared core-first
	# reference set once; subsequent raw MGS reuse it.
	if ($requiresOutgroupReference && !$outgroupReferenceInitialized) {
		$initializeOutgroupReferences->(\@fullTreeCandidates);
		print "Completed core-first sequential outgroup-reference loading; resuming per-MGS outgroup addition and tree submission; global elapsed "
			.timeNice(time - $^T)."\n";
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
	if (@curCogs < $minLociPerMGS) {
		limitedWarn('MGS with too few usable genes for tree construction',
			"$MGS has only ".scalar(@curCogs)." usable genes; skipping tree construction\n");
		return ($ingroupSampleCount, scalar(@curCogs), $OG, 1, 1, $ingroupSampleCount);
	}

	if ($treeFile ne "" || exists($PreferredOutgroup{$MGS})) {
		my $minimumOutgroupLoci = $outgroupDemandMinimum{$MGS} // $minLociPerMGS;
		my @requiredLoci = sort grep {
			$locusSeen{$_}
		} keys %{$outgroupDemandLoci{$MGS} || {}};
		$OG = $SelectedOutgroup{$MGS} // '';
		if (!length($OG)) {
			limitedWarn('MGS without a selected outgroup',
				"No viable predetermined outgroup remains for $MGS; building an ingroup-only tree\n");
		} elsif (@requiredLoci < $minimumOutgroupLoci) {
			limitedWarn('MGS with insufficient staged outgroup demand',
				"$MGS retains only ".scalar(@requiredLoci)." of $minimumOutgroupLoci required outgroup loci; building an ingroup-only tree\n");
			$OG = "";
		} else {
			my $represented = 0;
			for my $locus (@requiredLoci) {
				my (undef, $annotation) = locusParts($locus, $MGS);
				next if $annotation =~ m/^uniq\d+$/;
				my $outgroupGene = outgroupGeneForLocus($OG, $locus, $MGS);
				next unless length($outgroupGene) && exists($FNAref->{$outgroupGene});
				$represented++;
			}
			if ($represented < $minimumOutgroupLoci) {
				my @preview = @requiredLoci[0 .. ($#requiredLoci < 9 ? $#requiredLoci : 9)];
				limitedWarn('MGS without a sufficiently represented predetermined outgroup',
					"Predetermined outgroup $OG supplies $represented/$minimumOutgroupLoci required loci for $MGS; loci: @preview\n");
				$OG = "";
			}
		}
	}

	my ($overlayFNA, $overlayFAA, $overlayCategory, $outgroupGenes) = ('', '', '', 0);
	if ($OG ne '') {
		for my $locus (@curCogs) {
			next unless exists($outgroupEligibleLoci{$MGS}{$locus});
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
	return ($treeSampleCount, scalar(@curCogs), $OG, 1, 1, $ingroupSampleCount);
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
	my (@repairRequired, %repairState);
	my ($ready, $terminal, $excluded) = (0, 0, 0);
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


sub phase1SamplesByGroup {
	my (%samplesByGroup, %groupForSample);
	for my $sample (@samples) {
		die "Sample $sample has no assembly-group assignment\n"
			unless exists($map{$sample}) && defined($map{$sample}{AssGroup});
		my $group = $map{$sample}{AssGroup} ne '-1'
			? $map{$sample}{AssGroup} : "__standalone__".$sample;
		push @{$samplesByGroup{$group}}, $sample;
		$groupForSample{$sample} = $group;
	}
	return (\%samplesByGroup, \%groupForSample);
}

sub phase1EstimatedInputBytes {
	my ($path) = @_;
	my $nominal = $path;
	$nominal =~ s/\.gz\z//;
	return fileGZs($nominal);
}

sub phase1SampleWorkEstimate {
	my ($sample) = @_;
	return (1, 'skipped')
		unless exists($map{$sample}) && defined($map{$sample}{wrdir})
			&& length($map{$sample}{wrdir});
	my $cD = $map{$sample}{wrdir};
	return (1, 'skipped') if -e "$cD/SMPL.empty";

	my $localNT = "$cD/$lSNPdir/$lConsFNA";
	my $localAA = "$cD/$lSNPdir/$lConsFAA";
	my $vcf = "$cD/$lSNPdir/$lConsVCF";
	my ($readyNT, $readyAA);
	if (exists($preCompSNPs{$sample})
			&& fileGZe($preCompSNPs{$sample}{NT})
			&& fileGZe($preCompSNPs{$sample}{AA})) {
		($readyNT, $readyAA) = @{$preCompSNPs{$sample}}{qw(NT AA)};
	} elsif ($preCompCons && defined($preConDir) && length($preConDir)) {
		my $precomputedNT = "$preConDir/$sample.cons.genes.fna.gz";
		my $precomputedAA = "$preConDir/$sample.cons.prots.faa.gz";
		($readyNT, $readyAA) = ($precomputedNT, $precomputedAA)
			if fileGZe($precomputedNT) && fileGZe($precomputedAA);
	}
	if (!defined($readyNT) && !$forceVCF2FNA
			&& fileGZe($localNT) && fileGZe($localAA)) {
		($readyNT, $readyAA) = ($localNT, $localAA);
	}
	if (defined $readyNT) {
		my $bytes = phase1EstimatedInputBytes($readyNT) + phase1EstimatedInputBytes($readyAA);
		$bytes = 0 if $bytes < 0;
		# The fixed component retains per-sample consensus/depth overhead; fileGZs
		# estimates uncompressed bytes for gzip inputs, matching the FASTA scan cost.
		return (16 + $bytes / (1024 * 1024), 'ready');
	}
	return (1, 'skipped') unless fileGZe($vcf);
	my $vcfBytes = phase1EstimatedInputBytes($vcf);
	$vcfBytes = 1024 * 1024 if $vcfBytes < 1024 * 1024;
	# Regeneration also reads reference/depth/GFF inputs. A fixed penalty plus
	# twice the estimated VCF input prevents these samples from clustering merely
	# because their pre-existing consensus FASTAs are absent.
	return (64 + 2 * $vcfBytes / (1024 * 1024), 'regenerate');
}

sub phase1GroupWorkEstimates {
	my ($samplesByGroup) = @_;
	my (%groupWork, %groupMeta);
	for my $group (keys %{$samplesByGroup}) {
		my $work = 8; # indivisible assembly-driver/locus-model overhead
		my %counts = (ready => 0, regenerate => 0, skipped => 0);
		for my $sample (@{$samplesByGroup->{$group}}) {
			my ($sampleWork, $state) = phase1SampleWorkEstimate($sample);
			$work += $sampleWork;
			$counts{$state}++;
		}
		$groupWork{$group} = $work;
		$groupMeta{$group} = {
			samples => scalar(@{$samplesByGroup->{$group}}), %counts,
		};
	}
	return (\%groupWork, \%groupMeta);
}

sub writePhase1WorkerPlan {
	my ($path, $samplesByGroup, $workerForGroup, $groupWork,
		$groupMeta, $workerLoads) = @_;
	my $temporary = "$path.write.$$";
	my $out = retry_open('>', $temporary, label => 'create Phase-I worker plan');
	print {$out} "record\tworker\tassembly_group\tsamples\testimated_work\tready\tregenerate\tskipped\n";
	for my $group (sort {
		$workerForGroup->{$a} <=> $workerForGroup->{$b} || $a cmp $b
	} keys %{$samplesByGroup}) {
		my $safeGroup = $group;
		$safeGroup =~ s/[\t\r\n]+/_/g;
		my $meta = $groupMeta->{$group};
		print {$out} join("\t", 'group', $workerForGroup->{$group}, $safeGroup,
			$meta->{samples}, sprintf('%.2f', $groupWork->{$group}),
			@{$meta}{qw(ready regenerate skipped)}), "\n";
	}
	for my $worker (0 .. $#{$workerLoads}) {
		print {$out} join("\t", 'worker_total', $worker, '', '',
			sprintf('%.2f', $workerLoads->[$worker]), '', '', ''), "\n";
	}
	retry_close($out, 'close Phase-I worker plan');
	retry_rename($temporary, $path, label => 'publish Phase-I worker plan');
}

sub phase1SelectedGeneFingerprint {
	my $digest = Digest::SHA->new(256);
	$digest->add('strain-phase1-selected-genes-v1', "\0",
		phase1PathStatComponent($gene2taxF), "\0", $presortGenes, "\0");
	# The selected seed set is deterministic from the stat-bound gene2tax input,
	# its rank limit, and these selected MGS identifiers. Binding those inputs
	# avoids re-hashing millions of Gene2COG keys in every split worker.
	for my $mgs (@specis) {
		$digest->add('mgs', "\0", $mgs, "\0");
	}
	return $digest->hexdigest;
}

sub phase1IndexShardFingerprint {
	my ($clusterIndex, $workerForSample) = @_;
	return unless ref($workerForSample) eq 'HASH' && $maxSubJob > 0;
	my $digest = Digest::SHA->new(256);
	$digest->add('strain-phase1-cluster-index-shards-v1', "\0",
		phase1PathStatComponent($clusterIndex), "\0",
		phase1SelectedGeneFingerprint(), "\0", $maxSubJob, "\0");
	for my $sample (sort keys %{$workerForSample}) {
		$digest->add('sample', "\0", $sample, "\0",
			$workerForSample->{$sample}, "\0");
	}
	return $digest->hexdigest;
}

sub phase1IndexShardCacheState {
	my ($base, $fingerprint) = @_;
	my $manifest = File::Spec->catfile($base, 'manifest.tsv');
	return unless -s $manifest;
	open my $input, '<', $manifest or return;
	my @lines = <$input>;
	return unless close $input;
	s/\r?\n\z// for @lines;
	return unless @lines == $maxSubJob + 2
		&& $lines[0] eq join("\t", qw(version fingerprint worker_count))
		&& $lines[1] eq join("\t", 1, $fingerprint, $maxSubJob);
	my (@paths, @records);
	for my $worker (0 .. $maxSubJob - 1) {
		my @fields = split /\t/, $lines[$worker + 2], -1;
		my $expectedName = "worker.$worker.bin";
		return unless @fields == 5 && $fields[0] eq 'worker'
			&& $fields[1] eq "$worker" && $fields[2] =~ /\A\d+\z/
			&& $fields[3] =~ /\A\d+\z/ && $fields[4] eq $expectedName;
		my $path = File::Spec->catfile($base, $expectedName);
		my $bytes = -s $path;
		return unless -f $path && defined($bytes) && $bytes == $fields[3];
		push @paths, $path;
		push @records, 0 + $fields[2];
	}
	return (\@paths, \@records, $manifest);
}

sub publishPhase1IndexShards {
	my ($clusterIndex, $workerForSample, $base, $fingerprint) = @_;
	make_path($base) unless -d $base;
	my @final = map { File::Spec->catfile($base, "worker.$_.bin") }
		0 .. $maxSubJob - 1;
	my @temporary = map { "$_.tmp.$$" } @final;
	my $metadata;
	my $ok = eval {
		$metadata = writeClstrRevBinaryShards($clusterIndex, $Gene2COG,
			$workerForSample, \@temporary, $fingerprint);
		for my $worker (0 .. $maxSubJob - 1) {
			retry_rename($temporary[$worker], $final[$worker],
				label => "publish Phase-I cluster-index shard $worker");
		}
		my @manifest = (
			join("\t", qw(version fingerprint worker_count)),
			join("\t", 1, $fingerprint, $maxSubJob),
		);
		for my $record (@{$metadata}) {
			push @manifest, join("\t", 'worker', $record->{worker},
				$record->{records}, $record->{bytes},
				"worker.$record->{worker}.bin");
		}
		atomic_write_text(File::Spec->catfile($base, 'manifest.tsv'),
			join("\n", @manifest)."\n",
			label => 'publish Phase-I cluster-index shard manifest');
		1;
	};
	my $error = $@;
	for my $path (@temporary) {
		retry_unlink($path, fatal => 0,
			label => 'clean incomplete Phase-I cluster-index shard')
			if -e $path || -l $path;
	}
	die $error unless $ok;
	my $totalRecords = 0;
	my $totalBytes = 0;
	$totalRecords += $_->{records} for @{$metadata};
	$totalBytes += $_->{bytes} for @{$metadata};
	print "Published $maxSubJob atomic binary cluster-index shard(s): "
		."$totalRecords worker-cluster records, $totalBytes bytes; cache=$base\n";
	return 1;
}

sub loadPhase1ClusterIndex {
	my ($clusterIndex, $workerForSample, $mySamples) = @_;
	if ($maxSubJob > 0 && ref($workerForSample) eq 'HASH') {
		my $fingerprint = phase1IndexShardFingerprint($clusterIndex, $workerForSample);
		my $base = File::Spec->catdir($scratchD, 'phase1_cluster_index', $fingerprint);
		my ($paths, $records) = phase1IndexShardCacheState($base, $fingerprint);
		if (!$paths && !$subJob) {
			my $published = eval {
				publishPhase1IndexShards($clusterIndex, $workerForSample,
					$base, $fingerprint);
				1;
			};
			unless ($published) {
				my $error = $@ || 'unknown publication error';
				$error =~ s/\s+\z//;
				limitedWarn('Phase-I cluster-index shard publication',
					"Could not publish Phase-I cluster-index shards; using the full index as a compatibility fallback: $error\n");
			}
			($paths, $records) = phase1IndexShardCacheState($base, $fingerprint);
		}
		if ($paths) {
			# Retained so the catalogue-wide locus-model scan can stream every
			# worker's slice instead of materializing the whole catalogue.
			($phase1ShardPaths, $phase1ShardFingerprint) = ($paths, $fingerprint);
			my $clusters = eval {
				readClstrRevBinaryShard($paths->[$subJob], $fingerprint,
					$subJob, $maxSubJob);
			};
			if ($clusters) {
				print "Loaded Phase-I binary cluster-index shard $subJob/$maxSubJob: "
					."$records->[$subJob] cluster record(s), source=$paths->[$subJob]\n";
				return ($clusters, 'binary_worker_shard');
			}
			my $error = $@ || 'unknown shard validation error';
			$error =~ s/\s+\z//;
			limitedWarn('Phase-I cluster-index shard validation',
				"Phase-I cluster-index shard $subJob is invalid; using the full index as a compatibility fallback: $error\n");
		}
	} elsif ($maxSubJob > 0) {
		limitedWarn('Phase-I cluster-index shard assignment',
			"Worker/sample aliases cannot be assigned uniquely; using the full cluster index for this split generation\n");
	}
	my (undef, $clusters) = readClstrRev($clusterIndex, 0, $Gene2COG, $mySamples);
	return ($clusters, 'full_index_fallback');
}

sub phase1ProteinCacheFingerprint {
	my ($proteinFile) = @_;
	return sha256_hex(join("\0", 'strain-phase1-protein-cache-v1',
		phase1PathStatComponent($proteinFile),
		phase1SelectedGeneFingerprint()));
}

sub phase1ProteinCacheState {
	my ($base, $fingerprint) = @_;
	my $manifest = File::Spec->catfile($base, 'manifest.tsv');
	return unless -s $manifest;
	open my $input, '<', $manifest or return;
	my @lines = <$input>;
	return unless close $input;
	s/\r?\n\z// for @lines;
	return unless @lines == 2
		&& $lines[0] eq join("\t", qw(version fingerprint records bytes file));
	my @fields = split /\t/, $lines[1], -1;
	return unless @fields == 5 && $fields[0] eq '1'
		&& $fields[1] eq $fingerprint && $fields[2] =~ /\A\d+\z/
		&& $fields[3] =~ /\A\d+\z/ && $fields[4] eq 'catalog.proteins.bin';
	my $path = File::Spec->catfile($base, $fields[4]);
	my $bytes = -s $path;
	return unless -f $path && defined($bytes) && $bytes == $fields[3];
	return ($path, 0 + $fields[2], $manifest);
}

sub publishPhase1ProteinCache {
	my ($proteins, $base, $fingerprint) = @_;
	make_path($base) unless -d $base;
	my $final = File::Spec->catfile($base, 'catalog.proteins.bin');
	my $temporary = "$final.tmp.$$";
	my $metadata;
	my $ok = eval {
		$metadata = writeSequenceBinaryCache($proteins, $temporary, $fingerprint);
		retry_rename($temporary, $final,
			label => 'publish common Phase-I catalogue-protein cache');
		my $contents = join("\t", qw(version fingerprint records bytes file))."\n"
			.join("\t", 1, $fingerprint, $metadata->{records},
				$metadata->{bytes}, 'catalog.proteins.bin')."\n";
		atomic_write_text(File::Spec->catfile($base, 'manifest.tsv'), $contents,
			label => 'publish common Phase-I catalogue-protein manifest');
		1;
	};
	my $error = $@;
	retry_unlink($temporary, fatal => 0,
		label => 'clean incomplete Phase-I catalogue-protein cache')
		if -e $temporary || -l $temporary;
	die $error unless $ok;
	print "Published common Phase-I catalogue-protein cache: "
		."$metadata->{records} sequence(s), $metadata->{bytes} bytes; cache=$final\n";
	return 1;
}

sub loadPhase1CatalogProteins {
	my ($proteinFile) = @_;
	my $available = fileGZe($proteinFile) ? 1 : 0;
	return ({}, 0, 'unavailable') unless $available;
	# A common cache only amortizes its publication when Phase I is split. Avoid
	# adding a large scratch write to the ordinary single-process path.
	unless ($maxSubJob > 1) {
		return (readFasta($proteinFile, 1, "\\s", $Gene2COG, { fai => 1 }),
			1, 'catalogue_fasta');
	}
	my $fingerprint = phase1ProteinCacheFingerprint($proteinFile);
	my $base = File::Spec->catdir($scratchD, 'phase1_catalogue_proteins', $fingerprint);
	my ($cachePath, $expectedRecords) = phase1ProteinCacheState($base, $fingerprint);
	if ($cachePath) {
		my $proteins = eval { readSequenceBinaryCache($cachePath, $fingerprint) };
		if ($proteins && scalar(keys %{$proteins}) == $expectedRecords) {
			print "Loaded common Phase-I catalogue-protein cache: $expectedRecords sequence(s), source=$cachePath\n";
			return ($proteins, 1, 'binary_common_subset');
		}
		my $error = $@ || 'sequence count does not match the manifest';
		$error =~ s/\s+\z//;
		limitedWarn('Phase-I catalogue-protein cache validation',
			"Common Phase-I catalogue-protein cache is invalid; using the source FASTA as a compatibility fallback: $error\n");
	}
	my $proteins = readFasta($proteinFile, 1, "\\s", $Gene2COG, { fai => 1 });
	if (!$subJob) {
		my $published = eval {
			publishPhase1ProteinCache($proteins, $base, $fingerprint);
			1;
		};
		unless ($published) {
			my $error = $@ || 'unknown publication error';
			$error =~ s/\s+\z//;
			limitedWarn('Phase-I catalogue-protein cache publication',
				"Could not publish the common Phase-I catalogue-protein cache; continuing with the source FASTA result: $error\n");
		}
	}
	return ($proteins, 1, 'catalogue_fasta');
}

sub phase1LocusModelFingerprint {
	my ($clusterIndex, $proteinFile, $modelMGS) = @_;
	my $digest = Digest::SHA->new(256);
	# Everything that can change a locus boundary: the ranked seed input, the
	# cluster membership it is grouped over, the proteins compared during
	# merging, the confirmed Mosaic allowlist, and the post-grouping budget.
	# The workflow version is part of the identity: the model's contents depend on
	# how it was derived, not only on its inputs, so a model published by an
	# earlier build must never be reused silently after that derivation changes.
	$digest->add('strain-phase1-locus-model-v1', "\0", $version, "\0",
		phase1PathStatComponent($gene2taxF), "\0", $presortGenes, "\0",
		phase1PathStatComponent($clusterIndex), "\0",
		phase1PathStatComponent($proteinFile), "\0",
		$taxonAwareLocusSelection, "\0", $treeLocusBudget, "\0");
	$digest->add('mgs', "\0", $_, "\0") for @{$modelMGS};
	$digest->add('mosaic', "\0", $_, "\0") for sort keys %ConfirmedMosaicPairs;
	return $digest->hexdigest;
}

sub phase1LocusModelState {
	my ($base, $fingerprint) = @_;
	my $manifest = File::Spec->catfile($base, 'manifest.tsv');
	return unless -s $manifest;
	open my $input, '<', $manifest or return;
	my @lines = <$input>;
	return unless close $input;
	s/\r?\n\z// for @lines;
	my @columns = qw(version fingerprint groups group_bytes context_rows
		context_bytes ranked_records merged_seeds linkage_rejections budget_excluded);
	return unless @lines == 2 && $lines[0] eq join("\t", @columns);
	my @fields = split /\t/, $lines[1], -1;
	return unless @fields == @columns && $fields[0] eq '1'
		&& $fields[1] eq $fingerprint;
	return if grep { $fields[$_] !~ /\A\d+\z/ } 2 .. $#columns;
	my %state;
	@state{@columns[2 .. $#columns]} = map { 0 + $_ } @fields[2 .. $#columns];
	for my $part (['groups', 'locus_groups.tsv', 'group_bytes'],
			['context', 'locus_context.tsv', 'context_bytes']) {
		my ($label, $name, $sizeKey) = @{$part};
		my $path = File::Spec->catfile($base, $name);
		my $bytes = -f $path ? -s $path : undef;
		return unless defined($bytes) && $bytes == $state{$sizeKey};
		$state{"${label}_path"} = $path;
	}
	return \%state;
}

sub publishPhase1LocusModel {
	my ($groups, $locusContext, $base, $fingerprint, $counters) = @_;
	make_path($base) unless -d $base;
	my $groupPath = File::Spec->catfile($base, 'locus_groups.tsv');
	my $contextPath = File::Spec->catfile($base, 'locus_context.tsv');
	my $groupTemporary = "$groupPath.tmp.$$";
	my $contextTemporary = "$contextPath.tmp.$$";
	my ($groupRows, $contextRows) = (0, 0);
	my $ok = eval {
		my $groupOut = retry_open('>', $groupTemporary,
			label => 'create common Phase-I locus groups');
		for my $group (@{$groups}) {
			print {$groupOut} join("\t", $group->{mgs}, $group->{cog},
				$group->{primary_gene}, $group->{rank},
				join(",", @{$group->{genes}})), "\n"
				or die "Cannot write $groupTemporary: $!\n";
			$groupRows++;
		}
		retry_close($groupOut, 'close common Phase-I locus groups');
		retry_rename($groupTemporary, $groupPath,
			label => 'publish common Phase-I locus groups');

		my $contextOut = retry_open('>', $contextTemporary,
			label => 'create common Phase-I locus contexts');
		for my $group (@{$groups}) {
			my $context = $locusContext->{$group->{locus_id}} || {};
			for my $token (sort keys %{$context}) {
				print {$contextOut} join("\t", $group->{locus_id}, $token,
					$context->{$token}), "\n"
					or die "Cannot write $contextTemporary: $!\n";
				$contextRows++;
			}
		}
		retry_close($contextOut, 'close common Phase-I locus contexts');
		retry_rename($contextTemporary, $contextPath,
			label => 'publish common Phase-I locus contexts');
		1;
	};
	my $error = $@;
	for my $temporary ($groupTemporary, $contextTemporary) {
		retry_unlink($temporary, fatal => 0,
			label => 'clean incomplete common Phase-I locus model')
			if -e $temporary || -l $temporary;
	}
	die $error unless $ok;
	my $groupBytes = -s $groupPath;
	my $contextBytes = -s $contextPath;
	die "Cannot size published common Phase-I locus model\n"
		unless defined($groupBytes) && defined($contextBytes);
	my @columns = qw(version fingerprint groups group_bytes context_rows
		context_bytes ranked_records merged_seeds linkage_rejections budget_excluded);
	atomic_write_text(File::Spec->catfile($base, 'manifest.tsv'),
		join("\t", @columns)."\n"
		.join("\t", 1, $fingerprint, $groupRows, $groupBytes,
			$contextRows, $contextBytes,
			$counters->{ranked_records}, $counters->{merged_seeds},
			$counters->{linkage_rejections}, $counters->{budget_excluded})."\n",
		label => 'publish common Phase-I locus-model manifest');
	print "Published common Phase-I locus model: $groupRows locus/loci, "
		."$contextRows context record(s); cache=$base\n";
	return 1;
}

sub loadPhase1LocusModel {
	my ($base, $fingerprint, $expectedRecords) = @_;
	my $state = phase1LocusModelState($base, $fingerprint);
	return unless $state;
	# A model built for a different ranked-seed set cannot describe this
	# process's loci, so reject it rather than silently mixing generations.
	return if defined($expectedRecords) && $state->{ranked_records} != $expectedRecords;
	open my $groupsIn, '<', $state->{groups_path} or return;
	my (@groups, %seen);
	while (my $line = <$groupsIn>) {
		$line =~ s/[\r\n]+\z//;
		next unless length $line;
		my @fields = split /\t/, $line, -1;
		unless (@fields == 5 && length($fields[0]) && length($fields[1])
				&& length($fields[2]) && $fields[3] =~ /\A\d+\z/ && length($fields[4])) {
			close $groupsIn;
			return;
		}
		my ($mgs, $cog, $primary, $rank, $geneList) = @fields;
		my @genes = grep { length } split /,/, $geneList;
		my $locus = join('|', $mgs, $cog, $primary);
		if (!@genes || $genes[0] ne $primary || $seen{$locus}++) {
			close $groupsIn;
			return;
		}
		push @groups, {
			mgs => $mgs, cog => $cog, primary_gene => $primary,
			genes => \@genes, rank => 0 + $rank, locus_id => $locus,
		};
	}
	return unless close $groupsIn;
	return unless scalar(@groups) == $state->{groups};

	open my $contextIn, '<', $state->{context_path} or return;
	my (%locusContext, $contextRows);
	$contextRows = 0;
	while (my $line = <$contextIn>) {
		$line =~ s/[\r\n]+\z//;
		next unless length $line;
		my @fields = split /\t/, $line, -1;
		unless (@fields == 3 && length($fields[0]) && length($fields[1])
				&& $fields[2] =~ /\A\d+\z/ && exists($seen{$fields[0]})) {
			close $contextIn;
			return;
		}
		$locusContext{$fields[0]}{$fields[1]} = 0 + $fields[2];
		$contextRows++;
	}
	return unless close $contextIn;
	return unless $contextRows == $state->{context_rows};
	return (\@groups, \%locusContext, $state);
}

#Derive the catalogue-wide seed summaries the merge decisions need, without ever
#holding the whole catalogue's members. Synteny context never crosses a sample,
#and both summaries are additive over disjoint sample sets, so the already
#published per-worker shards can be streamed one at a time. Only when no shard
#set exists does this fall back to a single whole-catalogue read.
sub catalogueLocusContext {
	my ($accumulator, $records, $mergeCandidates, $clusterIndex, $scanOptions) = @_;
	$scanOptions ||= {};
	my $scanStarted = $scanOptions->{started} // time;
	# Sample sets are needed only for possible merges. Context is focal only for
	# records that can survive the final budget, while every ranked record remains
	# available as a neighbour on those focal contigs.
	my %options = (
		sample_set_seeds => $mergeCandidates,
		context_seeds => $scanOptions->{context_seeds},
		consume_cluster_members => 1,
	);
	if ($maxSubJob > 1 && ref($phase1ShardPaths) eq 'ARRAY'
			&& @{$phase1ShardPaths} == $maxSubJob) {
		my $scanned = 0;
		for my $worker (0 .. $maxSubJob - 1) {
			my $clusters = eval {
				readClstrRevBinaryShard($phase1ShardPaths->[$worker],
					$phase1ShardFingerprint, $worker, $maxSubJob);
			};
			unless ($clusters) {
				my $error = $@ || 'unknown shard read failure';
				$error =~ s/\s+\z//;
				limitedWarn('Phase-I locus-model shard scan',
					"Cannot stream cluster-index shard $worker; falling back to one "
					."whole-catalogue read: $error\n");
				$scanned = 0;
				last;
			}
			accumulate_locus_context($accumulator, $records, $clusters, \%options);
			$clusters = {};
			$scanned++;
			stepProgress("catalogue-wide locus-model scan", $scanned, $maxSubJob,
				$scanStarted,
				"context_seeds=".scalar(keys %{$accumulator->{gene_context} || {}}));
		}
		return "streamed_$scanned/${maxSubJob}_shards" if $scanned == $maxSubJob;
		%{$accumulator} = ();
	}
	my (undef, $fullClusters) = readClstrRev($clusterIndex, 0, $Gene2COG);
	accumulate_locus_context($accumulator, $records, $fullClusters, \%options);
	$fullClusters = {};
	return 'whole_catalogue_read';
}

#Group the ranked catalogue seeds into loci and apply the post-grouping budget.
#Shared by the parent's catalogue-wide build and the compatibility fallback.
sub buildSelectedLocusGroups {
	my ($records, $clusterMembers, $options) = @_;
	$options ||= {};
	my $locus_model = build_locus_groups(
		$records, $clusterMembers, $catalogProteins,
		{
			# These indexes are useful to general callers but duplicate large
			# parts of the cluster model and are not consumed by this workflow.
			include_member_to_seed => 0,
			include_gene_to_locus => 0,
			# The parent publishes member contexts per worker slice instead of
			# retaining a catalogue-wide copy it would immediately discard.
			include_member_context => $options->{member_context} ? 1 : 0,
			# A streamed catalogue scan already derived these seed summaries. Keep
			# the handoff explicit: silently dropping it publishes an empty context
			# model after spending hours computing the data.
			precomputed_context => $options->{precomputed_context},
			# A catalogue-wide membership map is discarded straight after the
			# build, so release its member strings as they are summarized.
			consume_cluster_members => $options->{consume} ? 1 : 0,
			# Every member pair in a multi-seed Mosaic locus must be independently
			# confirmed.  This prevents an A-B-C chain from silently merging A and C.
			allowed_merge_pairs => \%ConfirmedMosaicPairs,
			require_complete_linkage => 1,
			allow_confirmed_cooccurrence => 1,
		},
	);
	if ($options->{precomputed_context}
			&& keys %{$options->{precomputed_context}{gene_context} || {}}) {
		my $hasContext = 0;
		for my $context (values %{$locus_model->{locus_context} || {}}) {
			if (keys %{$context}) {
				$hasContext = 1;
				last;
			}
		}
		die "Catalogue-wide locus context was computed but grouping retained no context; "
			."refusing to publish an empty Phase-I context model\n"
			unless $hasContext;
	}
	my @groups = @{$locus_model->{groups}};
	my $budgetExcluded = $options->{prebudget_excluded} || 0;
	if (!$taxonAwareLocusSelection) {
		my %selected_loci_by_mgs;
		@groups = grep {
			if (($selected_loci_by_mgs{$_->{mgs}} // 0) >= $treeLocusBudget) {
				$budgetExcluded++;
				0;
			} else {
				$selected_loci_by_mgs{$_->{mgs}}++;
				1;
			}
		} @groups;
	}
	return (\@groups, $locus_model->{locus_context}, $locus_model->{member_context},
		{
			ranked_records => $options->{ranked_records_total} // scalar(@{$records}),
			merged_seeds => $locus_model->{merged_seeds} || 0,
			linkage_rejections => $locus_model->{incomplete_linkage_rejections} || 0,
			budget_excluded => $budgetExcluded,
		});
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
	my $workerForSampleHR = undef;
	if ($maxSubJob){
		# Partition whole assembly groups, never individual samples. The
		# catalogue has one shared-reference driver, while every member still
		# needs its sample-specific VCF/depth consensus below. Input size and
		# regeneration make equal sample counts poor proxies for elapsed time.
		my ($samplesByGroup, $groupForSample) = phase1SamplesByGroup();
		my ($groupWork, $groupMeta) = phase1GroupWorkEstimates($samplesByGroup);
		my @groups = sort keys %{$samplesByGroup};
		my ($workerForGroup, $workerLoads) =
			balance_assembly_groups($samplesByGroup, $maxSubJob, $groupWork);
		writePhase1WorkerPlan("$LOGDIR/phase1_worker_plan.tsv",
			$samplesByGroup, $workerForGroup, $groupWork, $groupMeta, $workerLoads)
			unless $subJob;
		my (%mine, %ownedGroup, %workerForSample);
		my $ambiguousAliasAssignment = 0;
		for my $group (@groups) {
			my $worker = $workerForGroup->{$group};
			$workerForSample{$_} = $worker for @{$samplesByGroup->{$group}};
			next unless $worker == $subJob;
			$ownedGroup{$group} = 1;
			$mine{$_} = 1 for @{$samplesByGroup->{$group}};
		}
		# Assembly catalogues can use generated aliases such as sampleM2.
		for my $alias (keys %{$map{altNms} || {}}) {
			my $sample = $map{altNms}{$alias};
			my $group = $groupForSample->{$sample};
			next unless defined $group;
			my $worker = $workerForGroup->{$group};
			if (exists($workerForSample{$alias})
					&& $workerForSample{$alias} != $worker) {
				$ambiguousAliasAssignment = 1;
			} else {
				$workerForSample{$alias} = $worker;
			}
			$mine{$alias} = 1 if $ownedGroup{$group};
		}
		$mySamplesHR = \%mine;
		$workerForSampleHR = \%workerForSample unless $ambiguousAliasAssignment;
		my $totalWorkerLoad = 0;
		$totalWorkerLoad += $_ for @{$workerLoads};
		my $plannedSamples = 0;
		$plannedSamples += scalar(@{$samplesByGroup->{$_}}) for keys %ownedGroup;
		print "Subjob ${subJob}/$maxSubJob: restricting locus-model construction to "
			. scalar(keys %ownedGroup)." of ".scalar(@groups)
			." assembly groups ($plannedSamples planned sample(s), "
			.scalar(keys %mine)." sample/alias identifiers; "
			."estimated work ".sprintf('%.2f', $workerLoads->[$subJob])."/"
			.sprintf('%.2f', $totalWorkerLoad)." units)\n";
	}

	my $cluster_index = "$GCd/compl.incompl.$clusterID.fna.clstr.idx";
	my $modelSubstepStarted = time;
	my $protein_file = "$GCd/compl.incompl.$clusterID.prot.faa";
	my ($loadedCatalogProteins, $protein_file_available, $proteinSource) =
		loadPhase1CatalogProteins($protein_file);
	$catalogProteins = $loadedCatalogProteins;
	if (!$protein_file_available) {
		warn "Catalogue protein file $protein_file is unavailable; keeping same-COG catalogue clusters separate\n";
	}
	stepComplete("locus-model catalogue-protein loading", $modelSubstepStarted,
		"worker=$subJob", "source_available=$protein_file_available",
		"loaded_proteins=".scalar(keys %{$catalogProteins}),
		"load_mode=$proteinSource",
		"catalogue=$protein_file");

	$modelSubstepStarted = time;
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
	my $ranked_record_count = scalar(@records);
	my %modelMGSseen;
	$modelMGSseen{$_->{mgs}} = 1 for @records;
	my @modelMGS = sort keys %modelMGSseen;
	my ($modelRecordsRef, $prebudgetExcluded, $prebudgetMergeCandidates) =
		(\@records, 0, undef);
	if ($maxSubJob > 1 && !$taxonAwareLocusSelection) {
		($modelRecordsRef, $prebudgetExcluded, $prebudgetMergeCandidates) =
			preselect_locus_records(\@records, $treeLocusBudget,
				\%ConfirmedMosaicPairs);
	}
	my $model_record_count = scalar(@{$modelRecordsRef});
	my %modelContextSeeds = map { $_->{gene} => 1 } @{$modelRecordsRef};
	print "Phase-I context pre-budget: $model_record_count/$ranked_record_count "
		."ranked seed(s) remain focal"
		.($prebudgetExcluded
			? "; $prebudgetExcluded lower-ranked seed(s) retained only as neighbours"
			: q{})
		."\n" if $maxSubJob > 1;

	# Same-COG seeds merge into one locus using catalogue-wide co-occurrence and
	# synteny context.  A sample-restricted shard cannot reconstruct either, so
	# workers that each build their own model can disagree about locus identity
	# and emit one biological locus under two names.  The parent therefore builds
	# the model once over the complete cluster index and publishes it; every
	# worker reuses that model and only derives the member contexts of its own
	# samples, which depend on one sample's contigs alone.
	my ($selectedGroupsRef, $modelLocusContext, $modelMemberContext, $modelCounters);
	my $modelSource = 'local_build';
	my ($modelBase, $modelFingerprint) = ('', '');
	# The parent's own sample slice is loaded first so that the catalogue-wide
	# scan below can stream the published shards one at a time. Holding the whole
	# catalogue's members at once costs this process roughly maxSubJob times the
	# memory of a worker and was enough to exhaust a 30 GB allocation.
	my $sliceStarted = time;
	my ($cl2gene, $clusterIndexSource) = loadPhase1ClusterIndex(
		$cluster_index, $workerForSampleHR, $mySamplesHR);
	stepComplete("locus-model cluster-index scan", $sliceStarted,
		"worker=$subJob",
		"sample_filter=".($mySamplesHR ? scalar(keys %{$mySamplesHR}) : 'all'),
		"represented_seed_clusters=".scalar(keys %{$cl2gene}),
		"source=$clusterIndexSource",
		"index=$cluster_index");
	if ($maxSubJob > 1) {
		$modelFingerprint =
			phase1LocusModelFingerprint($cluster_index, $protein_file, \@modelMGS);
		$modelBase = File::Spec->catdir($scratchD, 'phase1_locus_model', $modelFingerprint);
		my ($loadedGroups, $loadedContext, $loadedState) =
			loadPhase1LocusModel($modelBase, $modelFingerprint, $ranked_record_count);
		if ($loadedGroups) {
			($selectedGroupsRef, $modelLocusContext) = ($loadedGroups, $loadedContext);
			$modelCounters = {
				ranked_records => $loadedState->{ranked_records},
				merged_seeds => $loadedState->{merged_seeds},
				linkage_rejections => $loadedState->{linkage_rejections},
				budget_excluded => $loadedState->{budget_excluded},
			};
			$modelSource = 'published_common_model';
			print "Loaded common Phase-I locus model: ".scalar(@{$selectedGroupsRef})
				." locus/loci, source=$modelBase\n";
		} else {
			# Deriving the model from this process's own shard would give it locus
			# identities the other workers do not share, so read the complete
			# membership even when only the published copy is missing.
			limitedWarn('Phase-I locus-model cache',
				"No common Phase-I locus model at $modelBase; rebuilding it from the "
				."complete cluster index so this worker's locus identities still match "
				."the other workers\n") if $subJob;
			my $fullScanStarted = time;
			my %catalogueContext;
			my $mergeCandidates = $prebudgetMergeCandidates
				|| merge_candidate_seeds($modelRecordsRef, \%ConfirmedMosaicPairs);
			# The scan always runs, even with no mosaic pair to merge: every locus
			# needs its synteny context so that a sample offering two candidates
			# can still be resolved instead of losing the locus as ambiguous.
			my $scanScope = catalogueLocusContext(\%catalogueContext, \@records,
				$mergeCandidates, $cluster_index,
				{ started => $fullScanStarted, context_seeds => \%modelContextSeeds });
			print "Catalogue-wide locus-model scan: "
				.scalar(keys %{$catalogueContext{gene_context} || {}})
				." seed(s) with synteny context, "
				.scalar(keys %{$mergeCandidates})." merge candidate(s) needing a "
				."sample set; source=$scanScope; elapsed "
				.timeNice(time - $fullScanStarted)."\n";
			($selectedGroupsRef, $modelLocusContext, undef, $modelCounters) =
				buildSelectedLocusGroups($modelRecordsRef, {}, {
					member_context => 0,
					precomputed_context => \%catalogueContext,
					prebudget_excluded => $prebudgetExcluded,
					ranked_records_total => $ranked_record_count,
				});
			$modelSource = $subJob ? 'catalogue_wide_rebuild' : 'catalogue_wide_build';
			if (!$subJob) {
				# Publication must succeed before workers are submitted: without it
				# each of them would rebuild the model from the complete index, which
				# is correct but costs one full index scan per worker.
				my $published = eval {
					publishPhase1LocusModel($selectedGroupsRef, $modelLocusContext,
						$modelBase, $modelFingerprint, $modelCounters);
					1;
				};
				unless ($published) {
					my $error = $@ || 'unknown publication error';
					$error =~ s/\s+\z//;
					limitedWarn('Phase-I locus-model publication',
						"Could not publish the common Phase-I locus model to $modelBase; "
						."each split worker will rebuild it from the complete cluster "
						."index instead: $error\n");
				}
			}
		}
	}

	if ($selectedGroupsRef) {
		my %selectedContextSeeds;
		for my $group (@{$selectedGroupsRef}) {
			$selectedContextSeeds{$_} = 1 for @{$group->{genes} || []};
		}
		$modelMemberContext = member_context_map(\@records, $cl2gene,
			{ context_seeds => \%selectedContextSeeds });
	} else {
		($selectedGroupsRef, $modelLocusContext, $modelMemberContext, $modelCounters) =
			buildSelectedLocusGroups(\@records, $cl2gene,
				{ member_context => 1, consume => 0 });
	}
	my $linkage_rejections = $modelCounters->{linkage_rejections};
	print "Mosaic complete-linkage protection rejected $linkage_rejections "
		."transitive component merge(s)\n" if $linkage_rejections;
	@records = ();
	my @selected_locus_groups = @{$selectedGroupsRef};
	my $locus_budget_excluded = $modelCounters->{budget_excluded};
	$LocusByID = {
		map { $_->{locus_id} => $_ } @selected_locus_groups
	};
	$MemberContext = $modelMemberContext;
	$LocusContext = $modelLocusContext;

	my ($new_si_genes, $new_priorities) = ({}, {});
	for my $group (@selected_locus_groups) {
		$new_si_genes->{$group->{mgs}}{$group->{locus_id}} = $group->{primary_gene};
		push @{$new_priorities->{$group->{mgs}}}, $group->{locus_id};
	}
	$SIgenes = $new_si_genes;
	$COGprios = $new_priorities;
	stepComplete("locus-group construction", $modelSubstepStarted,
		"worker=$subJob", "ranked_clusters=$ranked_record_count",
		"context_focal_clusters=$model_record_count",
		"prebudget_excluded=$prebudgetExcluded",
		"resolved_loci=".scalar(@selected_locus_groups),
		"merged_seeds=$modelCounters->{merged_seeds}",
		"linkage_rejections=$linkage_rejections",
		"model_source=$modelSource");

	$modelSubstepStarted = time;
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
	if ($subJob) {
		# A worker exits after extraction, so from here it only ever consults
		# catalogue proteins through $LocusSeedProteins, and only for a locus that
		# offered more than one candidate in some sample. Holding a protein for
		# every selected gene costs gigabytes that nothing will read again.
		my %keptProteins;
		for my $locus (keys %contextLociNeeded) {
			my $group = $LocusByID->{$locus} or next;
			for my $seed (@{$group->{genes} || []}) {
				$keptProteins{$seed} = $catalogProteins->{$seed}
					if exists $catalogProteins->{$seed};
			}
		}
		my $droppedProteins = scalar(keys %{$catalogProteins}) - scalar(keys %keptProteins);
		$catalogProteins = \%keptProteins;
		# The ranked seed index and the gene/COG map were inputs to the locus
		# model; extraction reads neither, and Phase II never runs here.
		my $droppedSeedLoci = 0;
		$droppedSeedLoci += scalar(keys %{$SIgenes->{$_}}) for keys %{$SIgenes};
		my $droppedGeneCOG = scalar(keys %{$Gene2COG});
		$SIgenes = {}; $Gene2COG = {};
		print "Worker $subJob released post-model state: $droppedProteins catalogue "
			."protein(s), $droppedSeedLoci ranked seed locus/loci, $droppedGeneCOG "
			."gene/COG entries; retained ".scalar(keys %keptProteins)
			." protein(s) for ".scalar(keys %contextLociNeeded)." ambiguous locus/loci\n";
	}
	stepComplete("worker-member materialization", $modelSubstepStarted,
		"worker=$subJob", "catalogue_drivers=".scalar(keys %cl2gene2),
		"locus_sample_combinations=$gene_sample_combinations",
		"ambiguous_combinations=$ambiguous_seed_samples",
		"unrepresented_worker_loci=$unrepresentedWorkerLoci");
	print "Prepared ".scalar(@selected_locus_groups)." loci from $ranked_record_count"
		." ranked catalogue clusters; merged $modelCounters->{merged_seeds} compatible same-COG seeds. "
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


	print "\n!! WARNING !!: REDO mode '$redoMode' selected; matching existing results may be invalidated and rebuilt. !!\n"
		if $redoMode ne 'none';

	$mapF = $phase1MapSpec;
	
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
		print "Redoing incomplete strain inputs\n" if $redoMode eq 'input';
		print "Pre-creating ConsSNPs in $preConDir in $preCompCons runs\n" if ($preCompCons);
		print "SNP caller: $SNPcaller; consensus NT=$lConsFNA; AA=$lConsFAA; "
			."contig=$lConsCTG; VCF=$lConsVCF\n";
		print "-minSNPDepth $minSNPDepth, -minSNPCallQual $minSNPCallQual";
		print ", -SNPadaptiveQual $useAdaptiveQual, -SNPindelRangeFilt: $indelRange";
		if ($depthFilterScale){print ", depthFiltScale $depthFilterScale\n";}else {print "\n";}
		print "DiscTests=$discTests\n" unless ($discTests eq "");
		print "ContTests=$contTests\n" unless ($contTests eq "");
		print "familyVar=$familyVar\n" unless ($familyVar eq "");
		
		print "groupStabilityVars=$groupStabilityVars\n" unless ($groupStabilityVars eq "");
		print "MSAaligner: $MSAprog, tree GenesPerSpecies: $GenesPerSpecies, "
			."GeneLengthMin: $GeneLengthMin, GeneLengthIncludeMin: $GeneLengthIncludeMin, "
			."tree relativeNTFraction: $relativeNTFraction, "
			."taxonAwareLocusSelection: $taxonAwareLocusSelection\n";
		if ($strictBackbone) {
			print "Backbone placement enabled: GenesPerSpecies=$placementGenesPerSpecies, "
				."relativeNTFraction=$placementRelativeNTFraction, "
				."strictBackboneFraction=$strictBackboneFraction\n";
		} else {
			print "Backbone placement disabled; backbone- and placement-only filters are inactive\n";
		}
		print "Rate/GC partition merging: enabled=$rateMergePartitions, "
			."maximumBins=$rateMergeMaxBins, targetBin=$rateMergeTargetSites effective sites, "
			."minimumBin=$rateMergeMinLoci loci/"
			."$rateMergeMinSites sites\n";
		print "Taxon-aware locus hierarchy: candidatePool=$presortGenes"
			.", perSampleCap=".($noGeneLimit ? 'unlimited' : $maxNGenes)
			.((!$noGeneLimit && $maxNGenes > 0 && $maxNGenes < $presortGenes)
				? " (loci ranked beyond $maxNGenes are reachable only by samples"
					." that never fill the cap, so their prevalence reads low)" : "")
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
		print "Filtering defaults: MGSminGenesPSmpl=$MGStoolowGsThr, minLociPerMGS=$minLociPerMGS, "
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
		# Everything this run derives from the MGS guide belongs to the run, not to
		# the shared catalogue directory. Two strain_within runs over one catalogue
		# used to write - and delete - the same .srt/.gene2MGS pair beside the input.
		# The sorter names its output after the guide it is given, so staging the
		# guide inside the output directory keeps every derived file per-run without
		# touching anything the other run may be reading.
		my $stagedGuide = File::Spec->catfile($outD, basename($MGSfileOri));
		# The sorter also reads an optional occurrence table beside the guide, using
		# this exact derivation; stage it under the matching name or the prevalence
		# estimate silently changes.
		my $observedName = sub {
			my ($path) = @_;
			$path =~ s/\.core$//;
			return "$path.obs";
		};
		my $stageGuideInput = sub {
			my ($source, $target) = @_;
			return 0 unless -e $source;
			return 1 if -e $target || -l $target;
			my $resolved = abs_path($source);
			return 1 if defined($resolved) && symlink($resolved, $target);
			# A filesystem without symlinks still gets an exact copy.
			copy($source, $target)
				or die "Cannot stage MGS input $source as $target: $!\n";
			return 1;
		};
		my $stageGuide = sub {
			make_path($outD) unless -d $outD;
			$stageGuideInput->($MGSfileOri, $stagedGuide);
			$stageGuideInput->($observedName->($MGSfileOri),
				$observedName->($stagedGuide));
		};
		my $sortedMGS = "$stagedGuide.srt";
		if ($preparedMainBranchFastPath) {
			$MGSfile = -s $sortedMGS ? $sortedMGS : $MGSfile;
			$gene2taxF = "";
			print "Prepared-input recovery: skipped MGS sorting and gene-to-MGS index creation.\n";
		} else {
		if ($subJob) {
			die "Sorted MGS guide is missing for subjob: $sortedMGS\n" unless -s $sortedMGS;
		} elsif ($recalcTrees && !-s $sortedMGS) {
			die "-redo tree requires the existing sorted MGS guide: $sortedMGS\n";
		} elsif ($mode eq "MGSall" && !-e $sortedMGS) {
			assertSafeWorkflowRemoval($outD, $safeDefaultOutD, $GCd, $MGSfileOri, $bindir, getcwd()) if -d $outD;
			fastRemoveTree($outD);
			fastRemoveTree($scratchD);
			$stageGuide->();
			symlink($stagedGuide, $sortedMGS)
				or die "Cannot link $sortedMGS to $stagedGuide: $!\n";
		} elsif (!$onlySubmit || !-s $sortedMGS) {
			print "base files missing.. preparing complete resubmission and recalc of data\n";
			assertSafeWorkflowRemoval($outD, $safeDefaultOutD, $GCd, $MGSfileOri, $bindir, getcwd()) if -d $outD;
			fastRemoveTree($outD);
			fastRemoveTree($scratchD);
			$stageGuide->();
			my $sortMGSgenes = getProgPaths("sortMGSGeneImport_scr");
			my $cmd = $sortMGSgenes . " "
				. join(" ", map { shellQuote($_) } ($GCd, $stagedGuide, $useGTDBmg, $mode, $clusterID)) . "\n";
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
		fastRemoveTree($preConDir);
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
	if ($leanOnlySubmitResume && exists($COGprios->{$MGS})
			&& @{$COGprios->{$MGS}}) {
		my %seen;
		my @selected = grep { defined($_) && length($_) && !$seen{$_}++ }
			@{$COGprios->{$MGS}};
		return (\@selected, 'selected_gene_map_deferred_validation', undef)
			if @selected;
	}
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

sub treePreferredCoreGuide {
	# $MGSfile has been repointed at the sorted guide by the time trees are
	# submitted, so it carries the extraction priority order: the presorter's
	# ranking in MGS mode, and (via the symlink) the input's own order in MGSall
	# mode, which is still exactly the order extraction prioritised. Prefer it
	# over the pre-sort .core fallback, which is a membership table only.
	my ($fallbackGuide) = @_;
	$fallbackGuide = '' unless defined $fallbackGuide;
	return $MGSfile if length($MGSfile) && $MGSfile =~ /\.srt\z/ && -s $MGSfile;
	return $fallbackGuide;
}

sub readPreferredCoreGeneSet {
	my ($path) = @_;
	return {} unless defined($path) && length($path);
	open my $input, '<', $path
		or die "Cannot open preferred-core guide $path: $!\n";
	my (%genes, $ignored);
	while (my $line = <$input>) {
		$line =~ s/[\r\n]+\z//;
		next if $line eq '' || $line =~ /^\s*#/;
		my @field = split /\t/, $line, -1;
		if (@field < 2 || !length($field[1])) {
			$ignored++;
			next;
		}
		for my $gene (split /,/, $field[1]) {
			$gene =~ s/^\s+|\s+$//g;
			next unless length($gene);
			next if $gene =~ /\A(?:gene|gene_id|seed)\z/i;
			$genes{$gene} = 1;
		}
	}
	close $input or die "Cannot close preferred-core guide $path: $!\n";
	die "Preferred-core guide $path has no usable seed genes\n" unless keys %genes;
	print "Loaded ".scalar(keys %genes)." preferred universal-core seed gene(s) from $path for Phase II outgroup selection"
		.($ignored ? "; ignored $ignored malformed row(s)" : '')."\n";
	return \%genes;
}

sub preparedMainBranchInputSet {
	my ($guide, $outputDirectory, $subset) = @_;
	return (0, [], "output directory is absent")
		unless defined($outputDirectory) && -d $outputDirectory;
	return (0, [], "MGS guide is unspecified")
		unless defined($guide) && length($guide);
	# The sorted guide is a per-run product and now lives in the output directory.
	# Runs made before that still have it beside the input, so accept either.
	my $sortedGuide = $guide;
	unless ($guide =~ /\.srt\z/) {
		$sortedGuide = File::Spec->catfile($outputDirectory, basename($guide).'.srt');
		$sortedGuide = "$guide.srt" if !-s $sortedGuide && -s "$guide.srt";
	}
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
	return 1 if $leanOnlySubmitResume && -s $mergeCheckpoint;
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

sub msaOnlyArtifactsReady {
	my ($outputDirectory) = @_;
	return 0 unless defined($outputDirectory) && length($outputDirectory);
	my $marker = File::Spec->catfile($outputDirectory, 'msaOnly.complete.tsv');
	return 0 unless -s $marker;
	open my $markerHandle, '<', $marker or return 0;
	my $status = '';
	while (my $line = <$markerHandle>) {
		if ($line =~ /^status\t([^\r\n]+)/) {
			$status = $1;
			last;
		}
	}
	close $markerHandle or return 0;
	return 0 unless $status eq 'msa_complete';
	my $msaDirectory = File::Spec->catdir($outputDirectory, 'MSA');
	return 0 unless -d $msaDirectory;
	opendir my $msaHandle, $msaDirectory or return 0;
	my $ready = 0;
	while (my $name = readdir $msaHandle) {
		next if $name =~ /^MSAli/ || $name !~ /\.fna\.gz\z/;
		my $path = File::Spec->catfile($msaDirectory, $name);
		if (fileGZs($path)) {
			$ready = 1;
			last;
		}
	}
	closedir $msaHandle or return 0;
	return $ready;
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
				&& ($onlyMSA
					? msaOnlyArtifactsReady($outD2)
					: (-s $treeCompletion && fileGZs($completedTree)))) {
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
			fastRemoveTree($outD2);
			my $scratch_mgs = "$scratchD/outs/$MGS";
			fastRemoveTree($scratch_mgs);
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
		} elsif($onlyMSA
				? !msaOnlyArtifactsReady($outD2)
				: !fileGZs("$SIdirs{$MGS}/phylo/$treeFile")){
			$treeAbsent++;
			$deferredScratchCleanup{"$scratchD/outs/$MGS"} = 1
				if -d "$scratchD/outs/$MGS";
		} elsif($onlyMSA
				? msaOnlyArtifactsReady($outD2)
				: fileGZe("$SIdirs{$MGS}/phylo/$treeFile")) {
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

    my $wrMGS = 0;
    my $suffix = ".$subJob";
    my $baseOut = "$scratchD/outs";
	my @pendingMGS = grep {
		defined($OFstrH{$_}) && length($OFstrH{$_})
	} keys %OFstrH;
	my $pendingCount = scalar(@pendingMGS);
	my $flushStarted = time;
	my $nextFlushProgress = $flushStarted + 60;
	print "Flushing buffered MGS records: $pendingCount MGS to durable scratch $baseOut\n";

    foreach my $MGS (@pendingMGS) {

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
		if (time >= $nextFlushProgress) {
			stepProgress("buffered MGS publication", $wrMGS, $pendingCount,
				$flushStarted, "worker=$subJob");
			$nextFlushProgress = time + 60;
		}
    }

    $bufferedOutputBytes = 0;
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

sub stageStart {
	my ($stage, $description) = @_;
	print "\n========== $stage ==========\n";
	print "$description\n" if defined($description) && length($description);
	print "Global elapsed: ".timeNice(time - $^T)."\n";
	print "========================================\n";
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

sub strainOutputHasDurablePhaseIState {
	my ($outputDirectory, @durableFiles) = @_;
	return 1 if grep {
		defined($_) && length($_) && (-s $_ || -l $_)
	} @durableFiles;
	return 0 unless defined($outputDirectory) && -d $outputDirectory;
	opendir my $directory, $outputDirectory
		or die "Cannot inspect existing strain output $outputDirectory: $!\n";
	my %operationalDirectory = map { $_ => 1 }
		qw(LOGandSUB stones strainsScr1 .scratch);
	while (defined(my $entry = readdir $directory)) {
		next if $entry eq '.' || $entry eq '..'
			|| $operationalDirectory{$entry};
		my $candidate = File::Spec->catdir($outputDirectory, $entry);
		next unless -d $candidate;
		closedir $directory
			or die "Cannot close existing strain output $outputDirectory: $!\n";
		return 1;
	}
	closedir $directory
		or die "Cannot close existing strain output $outputDirectory: $!\n";
	return 0;
}

sub phase1PathStatComponent {
	my ($path, $contractVersion) = @_;
	$contractVersion //= $phase1InputContractVersion;
	return join("\0", '<unspecified>', 'missing')
		unless defined($path) && length($path);
	my $resolved = resolveExistingFile($path);
	return join("\0", $path, 'missing') unless defined($resolved);
	my $canonical = abs_path($resolved);
	die "Cannot resolve Phase-I input $resolved\n"
		unless defined($canonical) && length($canonical);
	my @metadata = stat($canonical);
	die "Cannot stat Phase-I input $canonical: $!\n"
		unless @metadata;
	# Keep this constant-cost. Contract v2 included st_dev, but a shared file can
	# have a different device number in different compute-node mount namespaces.
	# V3 therefore uses canonical path/inode/size/mtime. Omitting ctime prevents a
	# chmod alone from forcing an expensive Phase-I rebuild.
	my @stableMetadata = $contractVersion <= 2
		? @metadata[0, 1, 7, 9]
		: @metadata[1, 7, 9];
	return join("\0", $path, $canonical, @stableMetadata);
}

sub phase1GuideStatFingerprint {
	my ($path, $contractVersion, $outputDirectory) = @_;
	$contractVersion //= $phase1InputContractVersion;
	my $fingerprintSchema = $contractVersion <= 2
		? 'strain-phase1-guide-stat-v2'
		: 'strain-phase1-guide-stat-v3';
	return sha256_hex(join("\0", $fingerprintSchema, 'FMG'))
		unless defined($path) && length($path);
	my $canonical = abs_path($path);
	die "Cannot resolve Phase-I MGS/core guide $path\n"
		unless defined($canonical) && length($canonical);
	my $observation = $canonical;
	$observation =~ s/\.core\z//;
	$observation .= '.obs';
	# The sorted guide and its gene-to-MGS index are per-run products in the output
	# directory. Track them there, still accepting the pre-relocation layout beside
	# the input so an older run's contract keeps describing the same files.
	$outputDirectory = $outD
		if (!defined($outputDirectory) || !length($outputDirectory))
			&& defined($outD) && length($outD);
	my $staged = defined($outputDirectory) && length($outputDirectory)
		? File::Spec->catfile($outputDirectory, basename($canonical)) : q{};
	my $sorted = length($staged) && -s "$staged.srt" ? "$staged.srt" : "$canonical.srt";
	my @inputs = ($canonical, $sorted, "$sorted.gene2MGS", $observation);
	return sha256_hex(join("\0",
		$fingerprintSchema,
		map { phase1PathStatComponent($_, $contractVersion) } @inputs,
	));
}

sub phase1CatalogStatFingerprint {
	my ($catalog, $identity, $markerSet, $mapSpec, $contractVersion) = @_;
	$contractVersion //= $phase1InputContractVersion;
	my $markerFile = $markerSet eq 'GTDB'
		? "$catalog/GTDBmg.subset.cats" : "$catalog/FMG.subset.cats";
	my @inputs = (
		"$catalog/compl.incompl.$identity.fna",
		"$catalog/compl.incompl.$identity.fna.clstr.idx",
		"$catalog/compl.incompl.$identity.prot.faa",
		$markerFile,
		"$catalog/Anno/Func/emapper/eggNOGmapper_NOG.geneAss",
		grep { length } split(/,/, $mapSpec // ''),
	);
	my $fingerprintSchema = $contractVersion <= 2
		? 'strain-phase1-catalog-stat-v1'
		: 'strain-phase1-catalog-stat-v2';
	return sha256_hex(join("\0",
		$fingerprintSchema,
		map { phase1PathStatComponent($_, $contractVersion) } @inputs,
	));
}

sub phase1InputContractContents {
	my ($status) = @_;
	$status //= 'complete';
	die "Internal error: invalid Phase-I SNP-input contract status '$status'\n"
		unless $status eq 'building' || $status eq 'complete';
	my @columns = qw(version status catalog_identity catalog_inputs_fingerprint
		mgs_guide_fingerprint marker_set cluster_id snp_caller consensus_nt
		consensus_aa consensus_contig vcf vcf_support);
	my @values = ($phase1InputContractVersion, $status, $phase1CatalogIdentity,
		$phase1CatalogInputFingerprint, $phase1MGSGuideFingerprint,
		$useGTDBmg, $clusterID, $SNPcaller, $lConsFNA, $lConsFAA,
		$lConsCTG, $lConsVCF, $lConsVCFsup);
	return join("\t", @columns)."\n".join("\t", @values)."\n";
}

sub phase1InputContractState {
	my ($path) = @_;
	return ('missing', "no Phase-I SNP-input contract exists at $path")
		unless -e $path || -l $path;
	return ('invalid', "Phase-I SNP-input contract is not a nonempty file at $path")
		unless -f $path && -s $path;
	my $input = retry_open('<', $path, label => 'read Phase-I SNP-input contract');
	local $/;
	my $contents = <$input> // '';
	retry_close($input, 'close Phase-I SNP-input contract');
	my @lines = split /\n/, $contents, -1;
	pop @lines while @lines && $lines[-1] eq '';
	s/\r\z// for @lines;
	my $expectedHeader = join("\t", qw(version status catalog_identity
		catalog_inputs_fingerprint mgs_guide_fingerprint marker_set cluster_id
		snp_caller consensus_nt consensus_aa consensus_contig vcf vcf_support));
	return ('invalid', "Phase-I SNP-input contract has an invalid schema at $path")
		unless @lines == 2 && $lines[0] eq $expectedHeader;
	my @fields = split /\t/, $lines[1], -1;
	return ('invalid', "Phase-I SNP-input contract has an invalid record at $path")
		unless @fields == 13 && ($fields[0] eq '2'
				|| $fields[0] eq "$phase1InputContractVersion")
			&& ($fields[1] eq 'building' || $fields[1] eq 'complete');
	my $recordedVersion = 0 + $fields[0];
	my ($expectedCatalogFingerprint, $expectedGuideFingerprint) =
		$recordedVersion == $phase1InputContractVersion
		? ($phase1CatalogInputFingerprint, $phase1MGSGuideFingerprint)
		: (
			phase1CatalogStatFingerprint(
				$GCd, $clusterID, $useGTDBmg, $phase1MapSpec,
				$recordedVersion),
			phase1GuideStatFingerprint($MGSfileOri, $recordedVersion),
		);
	my @expectedInputs = ($phase1CatalogIdentity,
		$expectedCatalogFingerprint, $expectedGuideFingerprint,
		$useGTDBmg, $clusterID, $SNPcaller, $lConsFNA, $lConsFAA,
		$lConsCTG, $lConsVCF, $lConsVCFsup);
	my $inputsMatch = join("\0", @fields[2 .. 12])
		eq join("\0", @expectedInputs);
	if ($inputsMatch && $fields[1] eq 'complete') {
		return ('match', "Phase-I input contract matches the catalog, MGS/core guide, "
			."and -SNPcaller $SNPcaller", $recordedVersion, $fields[1]);
	}
	if ($inputsMatch) {
		return ('building_match', "Phase I for -SNPcaller $SNPcaller is marked incomplete at $path",
			$recordedVersion, $fields[1]);
	}
	my @mismatch;
	push @mismatch, 'catalog identity'
		if $fields[2] ne $phase1CatalogIdentity;
	push @mismatch, 'catalog or map input identity'
		if $fields[3] ne $expectedCatalogFingerprint;
	push @mismatch, 'MGS/core guide identity'
		if $fields[4] ne $expectedGuideFingerprint;
	push @mismatch, 'marker set or cluster identity'
		if $fields[5] ne $useGTDBmg || $fields[6] ne "$clusterID";
	push @mismatch, 'SNP caller or filenames'
		if join("\0", @fields[7 .. 12])
			ne join("\0", $SNPcaller, $lConsFNA, $lConsFAA,
				$lConsCTG, $lConsVCF, $lConsVCFsup);
	my $recordedCaller = $fields[7];
	$recordedCaller =~ s/[^A-Za-z0-9_.:+-]/?/g;
	my $legacyAdvice = $recordedVersion < $phase1InputContractVersion
		? "; legacy contract v$recordedVersion includes a node-local filesystem "
			."device identity; restart the main controller to upgrade it before "
			."dispatching workers"
		: "";
	return ('mismatch', "Phase-I SNP-input contract at $path records "
		."v$recordedVersion $fields[1] caller '$recordedCaller' with incompatible "
		.join(', ', @mismatch).$legacyAdvice,
		$recordedVersion, $fields[1]);
}

sub persistPhase1InputContract {
	my ($path, $status) = @_;
	$status //= 'complete';
	my ($state, undef, $recordedVersion) = phase1InputContractState($path);
	return 1 if $recordedVersion
		&& $recordedVersion == $phase1InputContractVersion
		&& $status eq 'complete' && $state eq 'match';
	return 1 if $recordedVersion
		&& $recordedVersion == $phase1InputContractVersion
		&& $status eq 'building' && $state eq 'building_match';
	make_path(dirname($path)) unless -d dirname($path);
	atomic_write_text($path, phase1InputContractContents($status),
		label => 'publish Phase-I SNP-input contract');
	print "Phase-I SNP-input contract: version=$phase1InputContractVersion; "
		."status=$status; caller=$SNPcaller; file=$path\n";
	return 1;
}

sub phase1WorkerCommand {
	my $strain1scr = getProgPaths("MGS_strain1_scr");
	my @selfArgs = (
		# Stage-I workers receive extraction, consensus, and extraction-relevant
		# locus controls; BuildTree model/submission flags stay in the parent.
		'-GCd', $GCd, '-outD', $outD, '-MGS', $MGSfileOri,
		'-clusterID', $clusterID, '-submit', 0, '-onlySubmit', 1,
		'-maxSubJob', $maxSubJob,
		'-MGSminGenesPSmpl', $MGStoolowGsThr,
		'-minLociPerMGS', $minLociPerMGS,
		'-multiGeneSmplMax', $multiGeneSmplMax,
		'-conspGeneSmplMax', $conspGeneSmplMax,
		'-minBadLociPSmpl', $minBadLociForSampleSkip, '-MGSphylo', $treeFile,
		'-presortGenes', $presortGenes, '-maxGenes', $maxNGenes,
		'-outgroupCoreMinLoci', $outgroupCoreMinLoci,
		'-outgroupReferenceGeneCap', $outgroupReferenceGeneCap,
		'-treeLocusBudget', $treeLocusBudget,
		'-taxonAwareLocusSelection', $taxonAwareLocusSelection,
		'-disableQC', $disableQC,
		'-breakpointGeneFlank', $breakpointGeneFlank,
		'-abundanceMinLoci', $abundanceMinimumLoci,
		'-abundanceMinFold', $abundanceMinimumFold,
		'-abundanceMaxFold', $abundanceMaximumFold,
		'-abundanceMaxModifiedZ', $abundanceMaximumModifiedZ,
		'-prepareMosaicLoci', $prepareMosaicLoci,
		'-flushEvery', $appendWriteTrigger,
		'-flushMemMB', $flushOutputMB,
		'-MGset', $useGTDBmg, '-SNPcaller', $SNPcaller,
		'-minSNPDepth', $minSNPDepth,
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

#Estimate what one Phase-I split worker needs, from the inputs it will actually
#hold. Every term is measured or counted rather than assumed, so a small run no
#longer queues the same fixed allocation as a catalogue-scale one, and a large
#run no longer starts below its working set and pays for an OOM round trip.
sub phase1AutoMemoryMB {
	my (%args) = @_;
	my $workers = $args{workers} || 1;
	$workers = 1 if $workers < 1;
	my $loci = $args{loci} || 0;
	my $samples = $args{samples} || 0;

	#the cluster index is the dominant per-worker structure; each worker loads
	#only the shard covering its own samples
	my $index = "$GCd/compl.incompl.$clusterID.fna.clstr.idx";
	my $indexBytes = -s $index;
	my $indexNote = "index";
	if (!$indexBytes) {
		my $compressed = -s "$index.gz";
		if ($compressed) {
			$indexBytes = $compressed * 4; #typical text ratio, decompressed size
			$indexNote = "index.gz x4";
		} else {
			$indexBytes = 0;
			$indexNote = "index unreadable";
		}
	}

	my $shardMB = ($indexBytes / $workers) * $phase1MemIndexFactor / (1024 * 1024);
	my $modelMB = $loci * $phase1MemLocusKB / 1024;
	my $sampleMB = ($samples / $workers) * $phase1MemSampleMB;
	my $estimate = $phase1MemBaseMB + $shardMB + $modelMB + $sampleMB;

	my $ceilingMB = int($treeOOMMaxMemGB * 1024 + 0.5);
	my $chosen = int($estimate + 0.999);
	my $clamped = '';
	if ($chosen < $phase1MemFloorMB) {
		$chosen = $phase1MemFloorMB; $clamped = ", raised to floor";
	} elsif ($chosen > $ceilingMB) {
		$chosen = $ceilingMB; $clamped = ", capped at -treeOOMMaxMemGB";
	}

	printf "Phase-I auto memory: %d MB per worker%s\n"
		."  base %d MB + index shard %.0f MB + locus model %.0f MB + samples %.0f MB\n"
		."  from %.2f GB %s / %d workers, %d loci, %d assembly groups\n",
		$chosen, $clamped, $phase1MemBaseMB, $shardMB, $modelMB, $sampleMB,
		$indexBytes / (1024 ** 3), $indexNote, $workers, $loci, $samples;
	warn "Phase-I auto memory could not read $index; the estimate excludes the "
		."cluster-index shard and may be low\n" if $indexNote eq "index unreadable";
	return $chosen;
}

#One place that answers "how much memory does a Phase-I worker start with",
#for both the initial submission and any later retry or resume.
sub phase1DefaultWorkerMemoryMB {
	my (%args) = @_;
	return int($selfMemGb * 1024 + 0.5) unless $selfMemAuto;
	return $phase1AutoMemoryCacheMB if $phase1AutoMemoryCacheMB;
	$phase1AutoMemoryCacheMB = phase1AutoMemoryMB(%args);
	return $phase1AutoMemoryCacheMB;
}

#Interval between accounting scans, expressed in seconds for the scheduler wait.
sub oomScanSeconds {
	return int($oomScanMinutes * 60 + 0.5);
}

#The OOM contract is per job: an accounting-confirmed OOM failure always keeps
#at least -oomMinRetries escalations of its own, no matter how many scans have
#already run or how many sibling jobs failed before it. Ordinary (non-OOM)
#failures keep the smaller -phase1WorkerRetries budget.
sub phase1RetryBudget {
	my ($oomConfirmed) = @_;
	return $oomConfirmed && $oomMinRetries > $phase1WorkerRetries
		? $oomMinRetries : $phase1WorkerRetries;
}

#One place that resubmits a single Phase-I worker, shared by the live OOM
#supervisor and by the ledger-driven retry loop so both keep the same job and
#memory bookkeeping.
sub submitPhase1Worker {
	my %args = @_;
	my $worker = $args{worker};
	my $memoryMB = $args{memory_mb};
	my $generation = $args{generation};
	my $workerCommand = $args{worker_command};
	my $label = $args{label};
	#Every Phase-I resubmission is recovery that the whole run waits on, so it
	#drops the bulk priority handicap unless the caller asks otherwise.
	my $nice = defined($args{job_nice}) ? $args{job_nice} : 0;
	my $jobByWorker = $args{job_by_worker} || {};
	my $memoryByWorker = $args{memory_by_worker} || {};
	my $stone = "$splitStonePrefix.$worker.stone";
	retry_unlink($stone, fatal => 0,
		label => "clear worker $worker completion");
	my $cmdX = "$workerCommand -subjob $worker &&\n"
		."printf '%s\\n' ".shellQuote($generation)
		." > ".shellQuote($stone)."\n";
	my $savedTmp = $QSBoptHR->{tmpSpace};
	my $savedNice = $QSBoptHR->{jobNice};
	$QSBoptHR->{tmpSpace} = 15;
	$QSBoptHR->{jobNice} = $nice;
	my ($dependency) = qsubSystem(
		"$LOGDIR/Strain1_B${worker}.${label}.sh",
		$cmdX, 1, $memoryMB."M", "Str1.$worker",
		"", "", 1, [], $QSBoptHR,
	);
	$QSBoptHR->{tmpSpace} = $savedTmp;
	$QSBoptHR->{jobNice} = $savedNice;
	$memoryByWorker->{$worker} = $memoryMB;
	my $jobID = slurm_job_id_from_dependency($dependency, $QSBoptHR->{rTag});
	if (defined($jobID)) {
		$jobByWorker->{$worker} = $jobID;
	} else {
		delete $jobByWorker->{$worker};
	}
	return $dependency;
}

#One accounting pass over the tracked Phase-I worker jobs. A worker Slurm has
#already killed for memory is resubmitted straight away with a doubled request,
#instead of waiting for the siblings that still have hours of work left.
#Returns the new dependencies, so the caller keeps waiting for them too.
sub escalatePhase1WorkerOOM {
	my %args = @_;
	return () unless $doSubmit && ($QSBoptHR->{qmode} || "") eq "slurm";
	my $jobByWorker = $args{job_by_worker} || {};
	my $memoryByWorker = $args{memory_by_worker} || {};
	my $retriesByWorker = $args{retries_by_worker} || {};
	my $handledJobs = $args{handled_jobs} || {};
	my $terminalOOMWorker = $args{terminal_workers} || {};
	my $defaultMemoryMB = $args{default_mb}
		|| phase1DefaultWorkerMemoryMB(workers => $maxSubJob);
	my @records;
	for my $worker (sort { $a <=> $b } keys %{$jobByWorker}) {
		my $jobID = $jobByWorker->{$worker};
		next unless defined($jobID) && $jobID =~ /^\d+$/;
		next if $handledJobs->{$jobID};
		push @records, {
			job_id => $jobID,
			worker => $worker,
			requested_mb => $memoryByWorker->{$worker} || $defaultMemoryMB,
		};
	}
	return () unless @records;
	my $plan = slurm_oom_retry_plan(\@records,
		int($treeOOMMaxMemGB * 1024 + 0.5));
	unless ($plan->{summary}{available}) {
		warn "Phase-I Slurm OOM accounting unavailable; deferring escalation to the "
			."next scan: ".($plan->{summary}{error} || "unknown error")."\n";
		return ();
	}
	my %workerByJob = map { $_->{job_id} => $_->{worker} } @records;
	my @dependencies;
	for my $jobID (sort { $a <=> $b } keys %{$plan->{by_job_id}}) {
		my $oom = $plan->{by_job_id}{$jobID};
		my $worker = $workerByJob{$jobID};
		next unless defined $worker;
		#Whatever the verdict, this outcome has now been ruled on: a later scan
		#must not rediscover it and escalate or warn about it a second time.
		$handledJobs->{$jobID} = 1;
		if ($oom->{ceiling_reached}) {
			warn "Cannot retry OOM Phase-I worker $worker: $oom->{requested_mb} MB "
				."already meets -treeOOMMaxMemGB $treeOOMMaxMemGB\n";
			$terminalOOMWorker->{$worker} = 1;
			next;
		}
		my $budget = phase1RetryBudget(1);
		my $spent = $retriesByWorker->{$worker} || 0;
		if ($spent >= $budget) {
			warn "Phase-I worker $worker exhausted its OOM retry budget "
				."($spent of $budget escalation(s) spent); its outcome stays "
				."quarantined for the ledger audit\n";
			next;
		}
		my $round = ++$retriesByWorker->{$worker};
		print "Phase-I worker $worker OOM escalation (attempt $round/$budget): "
			."$oom->{requested_mb} MB -> $oom->{next_mb} MB\n";
		push @dependencies, submitPhase1Worker(
			worker => $worker, memory_mb => $oom->{next_mb},
			generation => $args{generation},
			worker_command => $args{worker_command},
			label => ($args{script_kind} || "retry")."OOM${round}",
			job_by_worker => $jobByWorker,
			memory_by_worker => $memoryByWorker,
		);
	}
	return @dependencies;
}

#Wait for a wave of Phase-I workers, returning to the accounting scan every
#-oomScanMinutes so OOM failures are escalated while the wave is still running.
sub waitPhase1WorkersWithOOMScan {
	my %args = @_;
	my @pendingJobs = grep { defined($_) && length($_) } @{$args{jobs} || []};
	return unless @pendingJobs && $doSubmit;
	unless (($QSBoptHR->{qmode} || "") eq "slurm") {
		qsubSystemJobAlive(\@pendingJobs, $QSBoptHR);
		return;
	}
	my $scanSeconds = oomScanSeconds();
	my $scan = 0;
	while (@pendingJobs) {
		my $remaining = qsubSystemJobAlive(
			\@pendingJobs, $QSBoptHR, 0, -1, $scanSeconds);
		@pendingJobs = @{$remaining || []};
		$scan++;
		my @escalated = escalatePhase1WorkerOOM(%args);
		next unless @escalated;
		print "Phase-I OOM scan $scan resubmitted ".scalar(@escalated)
			." worker(s); ".scalar(@pendingJobs)." earlier job(s) still queued.\n";
		push @pendingJobs, @escalated;
	}
	return;
}

sub retryPhase1Workers {
	my %args = @_;
	my $generation = $args{generation};
	my @failedWorkers = @{$args{workers} || []};
	my $workerCommand = $args{worker_command};
	my $scriptKind = $args{script_kind} || "retry";
	my $jobByWorker = $args{job_by_worker} || {};
	my $memoryByWorker = $args{memory_by_worker} || {};
	#Budgets, ruled-on outcomes and ceiling verdicts are shared with the live OOM
	#supervisor: a worker escalated while its wave was still running must not get
	#a second, independent allowance once the ledger audit reaches it.
	my $retriesByWorker = $args{retries_by_worker} || {};
	my $handledJobs = $args{handled_jobs} || {};
	my $terminalOOMWorker = $args{terminal_workers} || {};
	die "Phase-I retry requires a generation and worker command\n"
		unless defined($generation) && length($generation)
			&& defined($workerCommand) && length($workerCommand);
	#On the resume path no memory ledger survives, so auto mode re-derives the
	#same baseline rather than silently dropping back to the fixed default.
	my $defaultMemoryMB = phase1DefaultWorkerMemoryMB(workers => $maxSubJob);
	my $maximumMemoryMB = int($treeOOMMaxMemGB * 1024 + 0.5);

	while (@failedWorkers) {
		last if grep { $_ == 0 } @failedWorkers;
		last unless $doSubmit;
		my @eligibleWorkers = grep { !$terminalOOMWorker->{$_} } @failedWorkers;
		last unless @eligibleWorkers;

		my %oomByJob;
		if (($QSBoptHR->{qmode} || "") eq "slurm") {
			my @records;
			for my $worker (@eligibleWorkers) {
				my $jobID = $jobByWorker->{$worker};
				next unless defined($jobID) && $jobID =~ /^\d+$/;
				next if $handledJobs->{$jobID};
				push @records, {
					job_id => $jobID,
					worker => $worker,
					requested_mb => $memoryByWorker->{$worker} || $defaultMemoryMB,
				};
			}
			if (@records) {
				my $plan = slurm_oom_retry_plan(\@records, $maximumMemoryMB);
				if (!$plan->{summary}{available}) {
					warn "Phase-I Slurm OOM accounting unavailable; retrying with the previous memory request: "
						.($plan->{summary}{error} || "unknown error")."\n";
				} else {
					%oomByJob = %{$plan->{by_job_id}};
				}
			}
		}

		my (%retryMemoryMB, %retryBudget);
		for my $worker (@eligibleWorkers) {
			my $currentMB = $memoryByWorker->{$worker} || $defaultMemoryMB;
			my $jobID = $jobByWorker->{$worker};
			my $oom = defined($jobID) ? $oomByJob{$jobID} : undef;
			if ($oom) {
				$handledJobs->{$jobID} = 1;
				if ($oom->{ceiling_reached}) {
					warn "Cannot retry OOM Phase-I worker $worker: $currentMB MB already meets "
						."-treeOOMMaxMemGB $treeOOMMaxMemGB\n";
					$terminalOOMWorker->{$worker} = 1;
					next;
				}
				$currentMB = $oom->{next_mb};
				print "Phase-I worker $worker OOM escalation: $oom->{requested_mb} MB -> "
					."$currentMB MB\n";
			}
			my $budget = phase1RetryBudget($oom ? 1 : 0);
			my $spent = $retriesByWorker->{$worker} || 0;
			if ($spent >= $budget) {
				warn "Phase-I worker $worker exhausted its retry budget "
					."($spent attempt(s) spent, budget $budget)\n";
				next;
			}
			$retryBudget{$worker} = $budget;
			$retryMemoryMB{$worker} = $currentMB;
		}
		my @retryableWorkers = grep { exists($retryMemoryMB{$_}) } @eligibleWorkers;
		last unless @retryableWorkers;

		my @retryJobs;
		for my $worker (@retryableWorkers) {
			my $round = ++$retriesByWorker->{$worker};
			print "Retrying Phase-I worker $worker "
				."(attempt $round/$retryBudget{$worker})\n";
			push @retryJobs, submitPhase1Worker(
				worker => $worker, memory_mb => $retryMemoryMB{$worker},
				generation => $generation, worker_command => $workerCommand,
				label => "${scriptKind}${round}",
				job_by_worker => $jobByWorker,
				memory_by_worker => $memoryByWorker,
			);
		}
		waitPhase1WorkersWithOOMScan(
			jobs => \@retryJobs, generation => $generation,
			worker_command => $workerCommand, script_kind => $scriptKind,
			job_by_worker => $jobByWorker, memory_by_worker => $memoryByWorker,
			retries_by_worker => $retriesByWorker, handled_jobs => $handledJobs,
			terminal_workers => $terminalOOMWorker,
			default_mb => $defaultMemoryMB,
		);
		@failedWorkers = phase1WorkersNeedingRetry($generation);
	}
	return \@failedWorkers;
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
		my $message = $leanOnlySubmitResume
			? "Lean tree-only resume: deferring Phase-I input checks to each MGS and skipping obsolete worker-ledger validation.\n"
			: "Tree-only resume: every MGS input passed the completed audit; skipping obsolete Phase-I worker-ledger validation and continuing to Phase II.\n";
		limitedNotice('tree-only resume skips obsolete Phase-I ledger validation', $message);
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
	if (@failedWorkers) {
		my $remaining = retryPhase1Workers(
			generation => $generation,
			workers => \@failedWorkers,
			worker_command => phase1WorkerCommand(),
			script_kind => "resume",
		);
		if (@{$remaining}) {
			my $queue = writePhase1RepairQueue($generation, $remaining,
				"resumed Phase-I worker validation failed");
			$completionMessage = "Phase I requires worker repair before Phase II; no tree jobs were submitted.";
			print "Phase-I recovery paused safely; repair queue: $queue. Invalid workers: "
				.join(",", @{$remaining})."\n";
			exit(0);
		}
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
	my @summaryPairs = (
		"samples=".($allSummary->{samples} // 0),
		"processed=".($allSummary->{processed_samples} // 0),
		"used_MGS=".($allSummary->{used_mgs} // 0)."/".($allSummary->{candidate_mgs} // 0),
		"retained_loci=".($allSummary->{retained_loci} // 0),
		"mean_loci_per_used_MGS=".($allSummary->{mean_loci_per_used_mgs} // 0),
		"skipped_MGS=".($allSummary->{skipped_mgs} // 0),
		"status=".($allSummary->{status_counts} // q{}),
	);
	print "STAGE I SAMPLE SUMMARY (all workers)\n";
	print join("; ", @summaryPairs), "\n";
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

sub writeGeneLengthSampleSummary {
	my @sourceColumns = qw(
		sample gene_length_min gene_length_include_min sample_prefilter_status
		input_loci gene_length_min_pass_loci gene_length_min_dropped_loci
		qc_pass_loci qc_dropped_loci gene_length_include_min_pass_loci
		gene_length_include_min_dropped_loci recovery_candidate_loci
		recovered_for_msa_loci gene_length_min_dropped_genes qc_dropped_genes
		gene_length_include_min_dropped_genes recovery_candidate_genes
		recovered_for_msa_genes
	);
	my $expectedHeader = join("\t", @sourceColumns);
	my @countColumns = qw(
		input_loci gene_length_min_pass_loci gene_length_min_dropped_loci
		qc_pass_loci qc_dropped_loci gene_length_include_min_pass_loci
		gene_length_include_min_dropped_loci recovery_candidate_loci
		recovered_for_msa_loci
	);
	my @geneColumns = qw(
		gene_length_min_dropped_genes qc_dropped_genes
		gene_length_include_min_dropped_genes recovery_candidate_genes
		recovered_for_msa_genes
	);
	my %sample;
	my $reports = 0;
	for my $mgs (@specis) {
		next unless defined($SIdirs{$mgs}) && length($SIdirs{$mgs});
		my $report = File::Spec->catfile(
			$SIdirs{$mgs}, 'phylo', 'gene_length_filter.samples.tsv');
		next unless -s $report;
		open my $input, '<', $report
			or die "Cannot read gene-length sample audit $report: $!\n";
		my $header = <$input> // '';
		$header =~ s/[\r\n]+\z//;
		die "Unexpected gene-length sample-audit header in $report\n"
			unless $header eq $expectedHeader;
		$reports++;
		while (my $line = <$input>) {
			$line =~ s/[\r\n]+\z//;
			next unless length $line;
			my @values = split /\t/, $line, -1;
			die "Wrong field count in gene-length sample audit $report\n"
				unless @values == @sourceColumns;
			my %row;
			@row{@sourceColumns} = @values;
			my $sampleId = $row{sample};
			die "Gene-length sample audit has no sample in $report\n"
				unless defined($sampleId) && length($sampleId);
			my $aggregate = $sample{$sampleId} ||= {
				gene_length_min => $row{gene_length_min},
				gene_length_include_min => $row{gene_length_include_min},
			};
			die "Inconsistent gene-length thresholds for sample $sampleId\n"
				unless $aggregate->{gene_length_min} eq $row{gene_length_min}
					&& $aggregate->{gene_length_include_min}
						eq $row{gene_length_include_min};
			$aggregate->{mgs_reports}++;
			if (($row{sample_prefilter_status} // '')
					eq 'removed_by_high_threshold_qc') {
				$aggregate->{mgs_removed_by_high_threshold_qc}++;
			} else {
				$aggregate->{mgs_retained_or_pending}++;
			}
			for my $column (@countColumns) {
				die "Non-numeric $column for $sampleId in $report\n"
					unless ($row{$column} // '') =~ /\A\d+\z/;
				$aggregate->{$column} += $row{$column};
			}
			for my $column (@geneColumns) {
				push @{$aggregate->{$column}}, map { "${mgs}:$_" }
					grep { length } split /,/, ($row{$column} // '');
			}
		}
		close $input or die "Cannot close gene-length sample audit $report: $!\n";
	}
	return 'not_available' unless $reports;

	my @outputColumns = (
		qw(sample gene_length_min gene_length_include_min mgs_reports
			mgs_retained_or_pending mgs_removed_by_high_threshold_qc),
		@countColumns, @geneColumns,
	);
	my $path = "$LOGDIR/strainGeneLengthFilter.samples.tsv";
	my $temporary = "$path.write.$$";
	my $output = retry_open('>', $temporary,
		label => 'create strain gene-length sample summary');
	print {$output} join("\t", @outputColumns), "\n"
		or die "Cannot write strain gene-length sample-summary header: $!\n";
	for my $sampleId (sort keys %sample) {
		my $row = $sample{$sampleId};
		$row->{sample} = $sampleId;
		$row->{$_} //= 0 for qw(
			mgs_reports mgs_retained_or_pending mgs_removed_by_high_threshold_qc
		), @countColumns;
		$row->{$_} = join(',', sort @{$row->{$_} || []}) for @geneColumns;
		print {$output} join("\t", map { $row->{$_} // '' } @outputColumns), "\n"
			or die "Cannot write strain gene-length summary for $sampleId: $!\n";
	}
	retry_close($output, 'close strain gene-length sample summary');
	retry_rename($temporary, $path,
		label => 'publish strain gene-length sample summary');
	print "Strain-wide gene-length sample audit: $path ($reports MGS reports)\n";
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
	my $geneLengthSampleSummary = writeGeneLengthSampleSummary();
	my $sampleHistograms = writeMGSSampleHistograms();
	my @lines = (
		"Strain-within recovery summary (v$version)",
		"output_directory\t$outD",
		"recovery_accounting\t".(-s $recovery ? $recovery : 'not_available'),
		"selection_attrition\t$selectionAttrition",
		"gene_length_sample_audit\t$geneLengthSampleSummary",
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
	print "SNP caller: $SNPcaller ($lConsFNA; $lConsFAA; $lConsVCF)\n";
	for my $inheritance (@pairedDefaultInheritances) {
		print "Paired option default: -$inheritance->{target}=$inheritance->{value} "
			."inherited from explicit -$inheritance->{source}\n";
	}
	print "Cores: $numCores (max: $maxCores); submit=$doSubmit; "
		."onlySubmit=$onlySubmit; redo=$redoMode; redoEPAfilter=$redoEPAfilter\n";
	print "Tree OOM recovery: rounds=$treeOOMRetryRounds; maximum memory=${treeOOMMaxMemGB}GB; "
		."per-thread memory scaling=cores/$treeMemThreadDivisor\n";
	print "OOM rescan: every $oomScanMinutes min while jobs are still running; at "
		."least $oomMinRetries escalation(s) per OOM job\n";
	print "Submission priority: ordinary jobs at --nice=$jobNice, OOM retries at "
		."--nice=0; live-job ceiling="
		.($maxQueuedJobs ? $maxQueuedJobs : 'none')."\n";
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


#File::Path::remove_tree lstats every entry from Perl; coreutils rm walks the
#same tree in C and is several times faster on the multi-GB scratch and output
#trees this workflow produces. Renaming first also makes the removal atomic for
#the caller: the old path is gone after one metadata operation, so a slow or
#interrupted unlink can never leave a half-emptied directory that still looks
#like usable state. Callers must have already run assertSafeWorkflowRemoval.
sub fastRemoveTree {
	my ($target, %options) = @_;
	return 0 unless defined($target) && length($target) && -d $target;
	my $parked = join('.', $target, 'deleting', time, $$, int(rand(1_000_000_000)));
	# The parked name is a sibling, so the rename stays on one filesystem.
	my $victim = rename($target, $parked) ? $parked : $target;
	my @victims = ($victim);
	# Sweep anything an interrupted earlier run parked but never finished. This is
	# also what recovers the space if a backgrounded unlink is killed with its job.
	push @victims, grep { -d $_ && $_ ne $victim }
		bsd_glob("$target.deleting.*") unless $target =~ /[*?\[\]{}]/;
	# Without a usable rm the unlink has to be done in-process, and then it must
	# also be synchronous. Probe once rather than per call.
	unless (defined $systemRemoveAvailable) {
		$systemRemoveAvailable = $^O eq 'MSWin32' ? 0
			: system('sh', '-c', 'command -v rm >/dev/null 2>&1') == 0 ? 1 : 0;
	}
	my $background = $systemRemoveAvailable && !$options{wait};
	my $removed = 0;
	for my $path (@victims) {
		my $status = -1;
		if ($systemRemoveAvailable) {
			# The shell starts rm and exits, so system() returns immediately and
			# the unlink is reparented instead of being left as a child to reap.
			# Only reclaiming the space is deferred: the caller's path is already
			# gone, so nothing downstream can observe the directory.
			$status = $background
				? system('sh', '-c', 'rm -rf -- "$1" &', 'sh', $path)
				: system('rm', '-rf', '--', $path);
		}
		if ($status != 0) {
			my $failed = !eval { remove_tree($path); 1 };
			if ($failed || -d $path) {
				warn "Could not remove $path: "
					.($@ || 'directory still present')."\n";
				next;
			}
		}
		$removed++;
	}
	return $removed;
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
		my ($tree, $stone, $terminalMarker, $pendingMarker, $msaOnly) =
			@{$expected->{$mgs}};
		my ($status, $reason) = ('failed_missing_output', '');
		my $outputReady = $msaOnly ? msaOnlyArtifactsReady(dirname($stone)) : -s $tree;
		if ($outputReady && -s $stone) {
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
	my $saved_nice = $options->{jobNice};
	delete $options->{capacityDeferred};
	delete $options->{capacityDeferralAnnounced};
	$options->{nonblockingMaxConcurrentJobs} = 1 unless $blocking;

	# OOM recovery takes the first tier and EPA recovery the second: an escalation
	# injected mid-wave must reach the scheduler before the ordinary jobs that are
	# still waiting for capacity, or it inherits their whole backlog. Within each
	# tier, submit jobs requesting the most cores first, then use the already
	# collected sample count-by-gene workload proxy as a descending tie-breaker.
	# These are the same values that selected cores and memory above, rather than
	# a second sizing scan.
	@{$queue} = sort {
		($b->{oom_retry} // 0) <=> ($a->{oom_retry} // 0)
			|| ($b->{epa_only} // 0) <=> ($a->{epa_only} // 0)
			|| ($b->{cores} // 0) <=> ($a->{cores} // 0)
			|| ($b->{workload_cells} // 0) <=> ($a->{workload_cells} // 0)
			|| ($b->{sample_count} // 0) <=> ($a->{sample_count} // 0)
			|| ($b->{requested_mb} // 0) <=> ($a->{requested_mb} // 0)
			|| ($b->{gene_count} // 0) <=> ($a->{gene_count} // 0)
			|| ($a->{priority_ordinal} // 0) <=> ($b->{priority_ordinal} // 0)
			|| ($a->{mgs} // '') cmp ($b->{mgs} // '')
	} @{$queue};
	while (@{$queue}) {
		my $record = $queue->[0];
		for my $required (qw(mgs script command cores sample_count gene_count workload_cells memory requested_mb priority_ordinal job_name tree stone terminal placement_pending)) {
			die "Queued tree job is missing '$required'\n"
				unless defined($record->{$required}) && length($record->{$required});
		}
		if ($options->{doSubmit}) {
			my @staleOutputs = $record->{epa_only}
				? qw(stone terminal)
				: $record->{msa_only}
					? qw(stone terminal placement_pending)
					: qw(stone tree terminal placement_pending);
			retry_unlink($record->{$_}, label => "clear stale tree-job $_")
				for @staleOutputs;
		}
		$options->{tmpSpace} = $record->{tmp_space};
		$options->{useLongQueue} = $record->{use_long_queue};
		#An OOM retry records nice 0, so it outranks the queued ordinary wave.
		$options->{jobNice} = defined($record->{job_nice})
			? $record->{job_nice} : $saved_nice;
		my ($dependency) = qsubSystem(
			$record->{script}, $record->{command}, $record->{cores},
			$record->{memory}, $record->{job_name}, "", "", 1, [], $options,
		);
		if (defined($dependency) && $dependency eq deferredSubmissionDependency()) {
			#The supervisor retries this every minute for as long as the wave
			#lasts, so report it on a slow cadence instead of once per attempt.
			if (time >= $nextCapacityNotice) {
				print "Scheduler capacity is full; retaining ".scalar(@{$queue})
					." prepared tree job(s) until live jobs fall below "
					.($options->{maxConcurrentJobs} || 0)."\n";
				$nextCapacityNotice = time + 600;
			}
			last;
		}
		shift @{$queue};
		$submitted++ if $options->{doSubmit};
		push @{$jobs}, $dependency if defined($dependency) && length($dependency);
		if ($options->{doSubmit} && ($options->{qmode} || '') eq 'slurm'
				&& defined($dependency)) {
			my $schedulerJobID = slurm_job_id_from_dependency(
				$dependency, $options->{rTag});
			push @{$accounting}, {
				job_id => $schedulerJobID, mgs => $record->{mgs},
				requested_mb => $record->{requested_mb},
				retry_round => $record->{retry_round} // 0,
				submission_record => { %{$record} },
			} if defined($schedulerJobID);
		}
	}
	$options->{nonblockingMaxConcurrentJobs} = $saved_nonblocking;
	$options->{tmpSpace} = $saved_tmp_space;
	$options->{useLongQueue} = $saved_long_queue;
	$options->{jobNice} = $saved_nice;
	return { submitted => $submitted, pending => scalar(@{$queue}) };
}

sub retryOOMTreeJobs {
	my %args = @_;
	my $accounting = $args{accounting} || [];
	my $options = $args{options} || {};
	my $maximumMB = $args{maximum_mb};
	my $maximumRounds = $args{maximum_rounds} // 8;
	my $minimumRounds = $args{minimum_rounds} // $oomMinRetries;
	my $scanSeconds = $args{scan_seconds} // oomScanSeconds();
	my $pendingQueue = $args{pending_queue} || [];
	my $submittedRef = $args{submitted_ref};
	my @pendingJobs = grep { defined($_) && length($_) } @{$args{jobs} || []};
	# The configured memory ceiling, not the round count, should decide when to
	# stop: doubling out of the 5 GB floor needs seven rounds to reach 512 GB.
	# The count is spent per MGS rather than per wave, because a periodic scan
	# escalates whichever tree died since the last pass: a shared counter would
	# let the first failures consume the retries a late failure still needs.
	$maximumRounds = $minimumRounds if $maximumRounds < $minimumRounds;
	die "OOM retry rounds must be between 0 and 12\n"
		unless $maximumRounds >= 0 && $maximumRounds <= 12;
	#Without Slurm accounting no escalation is possible, so fall back to the
	#historical behaviour: drain the whole queue, then wait for it.
	unless ($options->{doSubmit} && ($options->{qmode} || '') eq 'slurm') {
		if (@{$pendingQueue}) {
			my (@queuedJobs, @queuedAccounting);
			my $drain = dispatchPendingTreeJobs(
				queue => $pendingQueue, options => $options,
				jobs => \@queuedJobs, accounting => \@queuedAccounting,
				blocking => 1,
			);
			${$submittedRef} += $drain->{submitted} if ref($submittedRef) eq 'SCALAR';
			push @pendingJobs, grep { defined($_) && length($_) } @queuedJobs;
			push @{$accounting}, @queuedAccounting;
		}
		qsubSystemJobAlive(\@pendingJobs, $options)
			if @pendingJobs && $options->{doSubmit};
		return 0;
	}

	my (%retriesByMGS, %handledJob);
	my $retried = 0;
	my $scan = 0;
	my $summary = { available => 0, error => 'no accounting scan was reached' };
	my $nextScan = 0;
	while (1) {
		# Top up scheduler capacity from the retained queue first. Escalations
		# were unshifted onto its front and sort into their own leading tier, so
		# they reach the scheduler before the ordinary jobs still waiting here.
		if (@{$pendingQueue}) {
			my (@queuedJobs, @queuedAccounting);
			my $drain = dispatchPendingTreeJobs(
				queue => $pendingQueue, options => $options,
				jobs => \@queuedJobs, accounting => \@queuedAccounting,
				blocking => 0,
			);
			${$submittedRef} += $drain->{submitted} if ref($submittedRef) eq 'SCALAR';
			push @pendingJobs, grep { defined($_) && length($_) } @queuedJobs;
			push @{$accounting}, @queuedAccounting;
		}
		unless (time >= $nextScan) {
			last unless @pendingJobs || @{$pendingQueue};
			# A queue that is still draining needs a much shorter wait than the
			# accounting cadence, so freed capacity is refilled promptly.
			my $budget = @{$pendingQueue}
				? $treeQueueDrainProbeSeconds : ($nextScan - time);
			$budget = 1 if $budget < 1;
			if (@pendingJobs) {
				my $remaining = qsubSystemJobAlive(
					\@pendingJobs, $options, 0, -1, $budget);
				@pendingJobs = @{$remaining || []};
			} else {
				sleep($budget);
			}
			next;
		}
		$nextScan = time + $scanSeconds;
		$scan++;
		# Outcomes already ruled on are dropped from the query, so a retried or
		# refused job is never rediscovered by a later scan.
		my @candidates = grep { !$handledJob{$_->{job_id}} } @{$accounting};
		my $oomPlan = slurm_oom_retry_plan(\@candidates, $maximumMB);
		$summary = $oomPlan->{summary};
		print "\nTree OOM scan $scan: ".scalar(@pendingJobs)." job(s) with the "
			."scheduler, ".scalar(@{$pendingQueue})." awaiting capacity; "
			.scalar(@candidates)." outcome(s) inspected.\n";
		print format_slurm_tree_memory_summary($summary);
		my @oom = $summary->{available}
			? map { $oomPlan->{by_job_id}{$_} }
				sort { $a <=> $b } keys %{$oomPlan->{by_job_id}}
			: ();
		my @retryQueue;
		for my $oom (@oom) {
			$handledJob{$oom->{job_id}} = 1;
			my $original = $oom->{submission_record};
			unless (ref($original) eq 'HASH') {
				warn "Cannot retry OOM tree job $oom->{job_id}: submission record is unavailable\n";
				next;
			}
			my $nextMB = $oom->{next_mb};
			unless (defined($nextMB)) {
				warn "OOM retry ceiling reached for $original->{mgs}: "
					."$original->{requested_mb} MB already meets -treeOOMMaxMemGB "
					."$treeOOMMaxMemGB\n";
				next;
			}
			my $round = ($retriesByMGS{$original->{mgs}} || 0) + 1;
			if ($round > $maximumRounds) {
				warn "OOM retry budget of $maximumRounds round(s) is exhausted for "
					."$original->{mgs}; its outcome stays quarantined\n";
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
			$retriesByMGS{$original->{mgs}} = $round;
			#Recovery, not bulk work: take the leading dispatch tier locally and
			#drop the nice handicap so Slurm ranks it above the queued wave too.
			$retry{oom_retry} = 1;
			$retry{job_nice} = 0;
			$retry{retry_round} = $round;
			$retry{requested_mb} = $nextMB;
			$retry{memory} = $nextMB.'M';
			# The saved command still carries the original allowance, which is what
			# BuildTree reports and uses for EPA planning. Track the new request.
			my $retryIqMemMB = int($nextMB * 0.9);
			$retry{command} =~ s/(^|\s)-iqMemMB\s+\d+/$1-iqMemMB $retryIqMemMB/;
			$retry{job_name} = 'OOM'.$round.'.'.$retry{mgs};
			$retry{script} = File::Spec->catfile(
				$mgsDirectory, "treeCmd.oom_retry.$round.sh") unless $epaStage;
			push @retryQueue, \%retry;
			print "OOM retry round $round/$maximumRounds for $retry{mgs}: "
				."$original->{requested_mb} MB -> $nextMB MB; "
				.($epaStage ? 'EPA-only with 1 thread' : "$retry{cores} core full-tree resume")
				."\n";
		}
		if (@retryQueue) {
			#Ahead of the ordinary jobs still waiting for capacity, not behind.
			unshift @{$pendingQueue}, @retryQueue;
			$retried += scalar(@retryQueue);
			print "Injected ".scalar(@retryQueue)." OOM retry job(s) ahead of "
				.(scalar(@{$pendingQueue}) - scalar(@retryQueue))
				." queued ordinary tree job(s).\n";
		}
		last unless @pendingJobs || @{$pendingQueue};
	}
	print "Tree OOM recovery complete: $retried escalated retry job(s) across "
		.scalar(keys %retriesByMGS)." MGS in $scan accounting scan(s); "
		."budget was $maximumRounds round(s) per MGS.\n";
	if ($summary->{available} && @{$summary->{oom_jobs} || []}) {
		warn "Tree OOM recovery ended with ".scalar(@{$summary->{oom_jobs}})
			." OOM outcome(s) it could not escalate further; they stay quarantined "
			."for inspection\n";
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
		if ($perMGSgenes < $minLociPerMGS){
			$lowCandidateMGS++;
			limitedWarn("MGS with fewer than $minLociPerMGS candidate loci",
				"Only $perMGSgenes genes/COGs for MGS $MGS; MGS genes might be multi-copy\n")
				unless $maxSubJob;
		}
	}
	print "$lowCandidateMGS MGS have fewer than $minLociPerMGS candidate loci in this worker's sample slice; "
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
	
	
	my %sampleStatsSeen;
	my $extractionStarted = time;
	my $nextSampleProgress = $extractionStarted + 60;
	for my $sm (@srtdSmpls) {
		readGenesSample_Singl(
			$sm, $writeLink, $sttime, \$appCnt, undef, \%sampleStatsSeen,
		);
		if (time >= $nextSampleProgress) {
			stepProgress("consensus-gene extraction", $smCnt, scalar(@srtdSmpls),
				$extractionStarted, "worker=$subJob",
				"sample_rows=".scalar(keys %sampleStatsSeen));
			$nextSampleProgress = time + 60;
		}
		$smCnt++;
	}
	close $sampleStatsPartFH
		or die "Cannot close per-worker sample statistics $sample_stats_part: $!\n";
	undef $sampleStatsPartFH;
	print "Stage I sample accounting: ".scalar(keys %sampleStatsSeen)." sample row(s) saved to $sample_stats_part.\n";
	
	
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
	my $candidateLoci = 0;
	$candidateLoci += $_ for values %locMGScnt;

	
	
	
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
			limitedWarn(q{empty samples}, "Skipping $sd3: sample is marked empty\n");
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
			limitedNotice(q{on-the-fly consensus regeneration}, "Regenerating consensus FASTA for $sd3\n");
			#store these in scratch, uncompressed (much faster)
			$fastaf = "$locSpace/$sd3.cons.genes.fna";
			$fastafAA = "$locSpace/$sd3.cons.prots.faa";
			createConsFastas($cD, $sd3, $fastaf, $fastafAA, 1, 0);
		}
		#print "$fastaf\n";
		unless (fileGZe($fastaf) && fileGZe($fastafAA)){
			limitedWarn(q{incomplete consensus pairs}, "Skipping $sd3: regenerated consensus FASTA pair is incomplete\n");
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
			# Name the finding, not a disposition: whether a mixed-strain sample is
			# deleted, deferred to placement, or kept is decided downstream by
			# -excludeFlaggedSamples and -placeOnBackbone.
			my $sampleQCStatus = ($double_failure || $csp_failure)
				? 'mixed_strain' : 'single_strain';
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
				my $aaChunk = join("",@OAstr); my $ntChunk = join("",@OFstr);
				my $linkChunk = join("",@OLstr); my $catChunk = join("",@OCstr);
				$bufferedOutputBytes += length($aaChunk) + length($ntChunk)
					+ length($linkChunk) + length($catChunk);
				$OAstrH{$MGS} .= $aaChunk; $OFstrH{$MGS} .= $ntChunk;
				my $recovery_reason = $double_failure && $csp_failure
					? 'placement_ambiguous_and_conspecific'
					: $double_failure ? 'placement_ambiguous'
					: $csp_failure ? 'placement_conspecific' : 'passed_qc';
				writeRecoveryRow($MGS, $sd3, 'recovered', $recovery_reason,
					$locCnt, $sampleQCStatus, $double_failure, $csp_failure,
					$locMosaicCnt);
				$OLstrH{$MGS} .= $linkChunk; $OCstrH{$MGS} .= $catChunk;
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
			# A sample count alone does not bound these buffers: they hold every
			# MGS this worker has touched, so their size depends on how many loci
			# each sample contributes, not on how many samples have been seen.
			# Flushing on accumulated bytes keeps peak memory flat instead of
			# scaling with the number of MGS in the run.
			if (${$bufferedSamplesRef} >= $appendWriteTrigger
					|| $bufferedOutputBytes >= $flushOutputByteLimit) {
				appendWriteMGSgenes($writeLink);
				${$bufferedSamplesRef} = 0;
			}
		}
	}
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
