#!/usr/bin/env perl

#ARGS: ./buildTree.pl -fna [FNA] -faa [FAA] -cat [categoryFile] -outD [outDir] -cores [CPUs] -useEte [1=ETE,0=this script] -relativeNTFraction [filter]
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
#5.08: 22.7.26: restore partitions, require nonempty resume outputs, and activate per-locus overlap filtering
#5.09: run MSAfix through a validated temporary output before replacing an alignment
#5.10: consolidate runtime configuration, progress, and repetitive diagnostics
#5.11: recover to output-local work space when the requested temporary path is unusable
#5.12: validate persistent continuation checkpoints and restart incomplete stages
#5.13: isolate recoverable alignment failures to their individual loci
#5.14: add bounded-memory/auto-thread IQ-TREE mode, pathogen support, and legacy compatibility
#5.15: fall back from pathogen mode when alignments exceed CMAPLE's compiled length limit
#5.16: omit IQ-TREE -mem for incompatible partition-model runs
#5.17: sort allFNA/allFAA records by locus before buildTree compresses them
#5.18: use the allocated IQ-TREE thread count directly; avoid costly per-tree AUTO benchmarking
#5.19: infer a strict validated backbone and place sparse samples afterwards
#5.20: filter anomalous per-locus alignments before concatenation
#5.21: use native MSAfix locus QC and clean its temporary files on every exit
#5.22: own staged-input publication, node-local temp selection, and completion markers
#5.23: default to broad/inter-species locus filtering; make strain-level filtering explicit
#5.24: retain all prepared loci by default for broad phylogenies; preserve opt-in locus QC
#5.25: fingerprint input filters and allow rare nonempty categories in broad marker trees
#5.26: validate complete IQ-TREE outputs, retry numerical underflow safely, and relax strain backbones
#5.27: publish the placed strain tree by default and retain ML inference as .backbone.treefile
#5.28: centralize reusable alignment and per-method tree checkpoint state
#5.29: add opt-in taxon-aware two-stage locus selection around the existing MSA workflow
#5.30: enable taxon-aware locus selection by default
#5.31: keep IQ-TREE inference unrooted; root published output downstream
#5.32: reject stronger locus-level divergence outliers and enable partition merging by default for strain trees
#5.33: restore the fast fixed IQ-TREE model as the strain-tree default; retain AutoModel as an opt-in
#5.34: deterministically merge strain loci into rate/GC partition bins before IQ-TREE
#5.35: size deterministic rate/GC partitions by effective called sites, not locus count
#5.36: retain broad backbones while applying coverage criteria to sparse-sample placement
#5.37: use EPA-ng maximum-likelihood placement for sparse strict-backbone samples
#5.38: separate backbone admission from placement eligibility and retain sparse taxa through MSA
#5.39: consume overlap-filtered rate/GC metrics from MSAfix rather than rescanning loci in Perl
#5.40: add structured post-alignment reporting, fixed-model defaults, and resumable lifecycle outcomes
#5.41: add bounded filesystem retries, preflight validation, durable locus checkpoints, and isolated EPA recovery
#5.42: reuse durable completion state and MSAfix statistics; cache unchanged IQ-TREE validation
#5.43: summarize repetitive per-locus MSAfix cleaning diagnostics and remove temporary captured logs
#5.44: bound EPA-ng placement memory and worker threads for strain-backbone placement
#5.45: execute configured EPA-ng environment wrappers as shell code rather than quoted paths
#5.46: run every MSAfix alignment and locus-QC invocation with the configured MSA core count
#5.47: apply MSAfix v2.15 coding-NT technical-offset repair after protein-guided alignment
#5.48: treat uneven raw or partial alignment tails as missing in taxon-aware coordinate scoring
#5.49: pass parsed IQ-TREE models directly to EPA-ng
#5.50: replace EPA ulimit with memory-aware threads and bounded query chunks
#5.53: persist taxon-aware candidate exhaustion as a valid terminal no-tree outcome
#5.54: exclude EPA placements with pendant branches far outside the backbone distribution
#5.55: reuse a retained EPA jplace when normal continuation is missing only its placed tree
#5.56: pass fitted IQ-TREE GTR parameters directly to EPA-ng
#5.57: restore authoritative backbone branch lengths after EPA-ng placement
#5.58: graft EPA placements directly onto the persisted backbone tree
#5.59: force retained-jplace filtering through the ordinary continuation path
#5.60: accept bare and explicit numeric redo-EPA flags
#5.61: redo retained EPA filtering before alignment and inference startup
#5.62: inherit forced redo state when strain_within resubmits an older tree command
#5.63: forward fitted IQ-TREE parameters when available, warn on fallback, and report the EPA-ng command
#5.64: parse IQ-TREE 3 compact model tables without silently collapsing partitioned EPA models
#5.65: state that EPA-ng does not optimize missing symbolic-model parameters
#5.66: refit one unpartitioned GTR model on fixed partitioned-backbone topology for EPA-ng
#5.67: remove per-locus MSA artifacts and retain compressed concatenated alignments

#5.68: finalize strain staged category/QC/outgroup overlays in the tree job
#5.69: delegate fused strain worker-shard finalization to a standalone helper
#5.70: distinguish MSA-selection checkpoints from downstream tree-stage settings
#5.71: consolidate internal workflow policy and lifecycle checkpoints into one state file
#5.72: validate retained concatenated alignment checkpoints before resuming tree inference
#5.73: require broadly prevalent loci for taxon rescue and QC backfill
#5.74: retain per-locus nucleotide MSAs when -rmMSA 0 is requested
#5.75: prefer universal-core guide loci and consolidate final taxon-aware diagnostics
#5.76: require 10k shared backbone sites for sparse-sample placement by default
#5.77: recover partial sample loci after high-threshold QC and audit both length gates
#5.80: pass explicit IQ-TREE sequence types and report retained-MSA resume counts
#5.81: keep plain MSAs and partitions on scratch and publish only compressed checkpoints
#5.82: exit -onlyMSA after localized per-locus processing and before combined-MSA work
#5.83: repair -subsetSmpls output naming, syn/nonsyn outgroup exclusion and codon
#      tolerance, the all-sites constraint tree, -rmMSA 0 protein-MSA retention,
#      and compile the sample/gene separator pattern once
use warnings;
use strict;
#use threads ('yield','stack_size' => 64*4096,'exit' => 'threads_only','stringify');
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::FlagReference qw(printFlagHelp helpRequested);
use Mods::GenoMetaAss qw( fileGZe fileGZs gzipopen systemW readFasta readFastHD writeFasta quantile);
use Mods::phyloTools qw(convertMSA2NXS MSA filterMSA getTreeLeafs calcDisPos2 runRaxML runRaxMLng runQItree 
			runFasttree runVeryFasttree iqtreeOutputComplete cleanupIQTreeTransients
			fixHDs4Phylo getGenoGenes getFMG readFMGdir );
use Mods::PhyloAlignment qw(filter_alignment_by_overlap);
use Mods::StrainPlacement qw(
	read_sample_qc canonical_sample_qc_status split_strict_backbone
	read_epa_jplace filter_epa_placement_outliers map_epa_placements_to_backbone write_epa_placed_tree
);
			
			
use Getopt::Long qw( GetOptions Configure );
#use Mods::ext::TreeIO;
#use Mods::IO::MaybeXS qw(encode_json decode_json);
use Mods::IO::PP qw (decode_json);
#use JSON qw( decode_json ); 
use Data::Dumper;
use Mods::math qw (medianArray avgArray meanArray);
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);
use File::Copy qw(copy move);
use File::Glob qw(bsd_glob);
use File::Path qw(make_path remove_tree);
use File::Spec;
use File::Temp qw(tempfile);
use FindBin qw($Bin);
use Mods::WorkflowResilience qw(
	retry_operation retry_unlink retry_rename retry_open retry_close
	preflight_directory filesystem_capacity
);
use Mods::StrainParts qw(append_fasta_records_atomic);


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
sub requestedTreeMethods;
sub treeMethodState;
sub parseSeqId;
sub compileSampleSeparator;
sub geneFileStem;
sub safeRemoveTree;
sub cachedIQTreeOutputComplete;
sub requireConfiguredTool;
sub shellQuote;
sub runMSAFix;
sub summarizeMSAFixLog;
sub limitedWarn;
sub prepareTemporaryBase;
sub sortFastaForCompression;
sub fastaCompressionSortKey;
sub runPostAlignmentLocusQC;
sub selectTaxonAwareCandidateLoci;
sub selectTaxonAwareFinalLoci;
sub chooseTaxonAwareLoci;
sub classifyTaxonAwareSamples;
sub classifyTaxonAwareCoverageEligibility;
sub taxonAwareAlignmentMetrics;
sub informativeSequenceLength;
sub bestGeneSequencesBySample;
sub writeTaxonAwareLocusAudit;
sub writeTaxonAwareSampleAudit;
sub preferredCoreGeneSet;
sub catalogueGeneFromLocus;
sub presortRankFromLocus;
sub assignPresortScores;
sub compactTaxonAwareDiagnostics;
sub writeSelectionAttritionAudit;
sub writeGeneLengthSampleAudit;
sub legacyPolicyFileMatches;
sub alignmentFileStem;
sub readPostAlignmentRateMetrics;
sub deterministicRatePartitions;
sub writeRatePartitionAudit;
sub publishStagedTreeInputs;
sub runStagedStrainShardHelper;
sub stagedTreeInputFiles;
sub prepareStagedStrainInputs;
sub finalizeStagedStrainCategory;
sub finalizeStagedSampleQC;
sub stagedOverlayRecordsPresent;
sub stagedOverlayCompletion;
sub writeStagedPreparationMarker;
sub writeCompletionMarker;
sub reusableCompletionTree;
sub writeOutcomeMarker;
sub completeTaxonAwareOutgroupAnchorTerminal;
sub clearLifecycleMarker;
sub preflightBuildTree;
sub treeAlignmentCheckpointStatus;
sub restoreCompressedMSAArtifact;
sub finalizeMSAArtifacts;
sub msaOnlyArtifacts;
sub inputFingerprint;
sub epaModelArtifact;
sub epaRefitIqtreeModel;
sub iqtreeGtrPartitionCount;
sub runEpaNgPlacement;
sub readStrictBackboneClassification;
sub runEpaOnlyPlacement;
sub runRedoEpaFilter;
sub writeEpaBackboneGraftAudit;
sub printEpaBackboneGraftSummary;
sub writeEpaPlacementFilterSummary;
sub printEpaPlacementFilterSummary;
sub epaFilterMetricValue;
sub readEpaFilterBackboneTree;
sub epaResourcePlan;
sub iqtreePlacementModel;
sub iqtreeExplicitEpaModel;
sub postAlignmentStep;
sub elapsedTimeText;
sub alignmentCollectionStats;
sub alignmentCollectionStatsFromReport;
sub partitionLocusRangeCount;
sub rawCoordinateInformation;

sub readBuildTreeState;
sub buildTreeStatePolicyMatches;
sub writeBuildTreeState;
sub cleanupLegacyBuildTreeStateFiles;
sub writeWorkflowHeartbeat;
sub writeWorkflowFailure;
my $doPhym= 0;
my $version = "5.83";
my %iqtreeValidationCache;
my %limitedWarningCounts;
my %limitedWarningLimits;
my $synSummaryCount = 0;
my $synSiteTotal = 0;
my ($msaFixCleanedLoci, $msaFixBorderMasked, $msaFixLowIdMasked, $msaFixMinGoodRemoved) = (0, 0, 0, 0);
my $nonSynSiteTotal = 0;
my ($workflowStage, $workflowStateFile, $workflowStatus, $workflowReason) =
	(q{initialization}, q{}, q{running}, q{});
my ($workflowMsaSelectionPolicy, $workflowTreeStagePolicy) = (q{}, q{});

END {
	writeWorkflowFailure($@ || 'non-zero process exit') if $? != 0;
	for my $category (sort keys %limitedWarningCounts) {
		my $limit = $limitedWarningLimits{$category} // 5;
		my $suppressed = $limitedWarningCounts{$category} - $limit;
		warn "Suppressed $suppressed additional '$category' warning(s)\n"
			if $suppressed > 0;
	}
	if ($msaFixCleanedLoci) {
		warn "MSAfix cleaning summary: loci=$msaFixCleanedLoci, border-gap-masked=$msaFixBorderMasked, low-ID-masked=$msaFixLowIdMasked, below-minimum-good-positions=$msaFixMinGoodRemoved\n";
	}
}

#answered before the first getProgPaths() call, so -help needs no site config
if (helpRequested(@ARGV)) {
	printFlagHelp(
		script  => "buildTree5.pl",
		version => $version,
		usage   => ["buildTree5.pl -fna FILE -aa FILE -cats FILE -outD DIR [options]",
			"buildTree5.pl -genoInD DIR -outD DIR [options]",
			"buildTree5.pl -help | -h | -?"],
		summary => "Phylogenetic tree construction and related MSA/population-genetic analyses. "
			."Between-species trees are the default; use -withinSpecies 1 for strain trees.",
		exit    => 1,
	);
}

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
my $ntFracGeneInclude = 0.03;
my $GeneFracPSpec = 0.1; #replacement for ntFracGene, as works also with supertrees
my $MSAprog = 2; #do MSA with clustal (1) or msaprobs (0), mafft(2), guidance2(3), MUSCLE5 (4)
my $calcDistMat = 0; #distmat of either AA or NT (depending on MSA)
my $calcDistMatExt = 0; #distmat of other AA or NT (depending on MSA), e.g. running two times an MSA
my $calcDistMatExtGo = 0;
my $treeAutoModel=0; #fixed GTR+F+G2 by default; automatic model selection is opt-in
my $treeAutoModelExplicit=0;
my $fracMaxGenesFilter = 0.2;
my $fracMaxGenes90pct = 0.25; #gene cats to keep, e.g. 25% of 90th percentile


my $ntCntTotal =0; my $bootStrap=0; my $subsetSmpls = -1;
my ($fnFna, $aaFna,$cogCats,$outD,$ncore,$Ete)= ("","","","",1,0);
my ($smplDef,$smplSep,$calcSyn,$calcNonSyn,$useAA4tree,$calcDNAdiff,$tmpD ) = (1,"_",0,0,0,0,"");
my $sampleGeneRegex; #compiled once from $smplSep by compileSampleSeparator
my ($stagedInputDir, $tmpSubdir, $completionMarker) = ("", "", "");
my $onlyMSA = 0;
my ($terminalMarker, $placementPendingMarker) = ("", "");
my $withinSpecies = 0;
my $strainWithinPreset = 0;

my ($continue,$isAligned) = (0,0);#overwrite already existing files?
my $epaOnly = 0;
my $redoEPAfilter =
	($ENV{MATAFILER_REDO_EPA_FILTER} // '') eq '1' ? 1 : 0;
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
my $iqFast=1;
my $iqMemMB=0;
my $iqPathogen=0;
my $iqLegacy=0;
my $iqLegacyExplicit=0;
my $minOverlapMSA;
my $msaFixRecoverTechnicalOffsets = 1; #coding-NT repair; applies only to MSAfix single-alignment mode
my $msaFixCodingFrame = 1;
my $msaFixGeneticCode = 11; #bacterial/archaeal/plastid genetic code
my $msaFixRecoveryBand = 3;
my $maxGapPerCol = 1 ;
my $minPcId = 0;
my $doSuperTree =0;
my $doSuperCheck=0;#check if tree's of single genes behave "strange"
my $gzipInput =0; my $removeMSA = 1; # retain per-locus nucleotide and protein MSAs only when explicitly requested
my $useTreeShrink =0;
my %BACKBONE_DEFAULT = (
	enabled => 0,
	coverage_fraction => 0.35,
	minimum_overlap => 10_000,
	minimum_samples => 3,
);
my $strictBackbone = $BACKBONE_DEFAULT{enabled};
my ($placeOnBackboneSpecified, $legacyStrictBackboneSpecified) = (0, 0);
my $strictBackboneFraction = $BACKBONE_DEFAULT{coverage_fraction};
my $placementMinOverlap = $BACKBONE_DEFAULT{minimum_overlap};
my $strictBackboneMinSamples = $BACKBONE_DEFAULT{minimum_samples};
my ($placementGeneFracPSpec, $placementNTFrac, $placementNTCntTotal);
my %EPA_NG_DEFAULT = (threads => 4, memory_fraction => 0.60,
	memory_per_thread_mb => 1024, chunk_size => 16,
	pendant_outlier_factor => 5, pendant_minimum_threshold => 0.02);
my $epaThreads = $EPA_NG_DEFAULT{threads};
my $epaMaxMemMB = -1; # -1 derives a thread-planning budget; 0 disables memory-based thread scaling
my $epaPendantOutlierFactor = $EPA_NG_DEFAULT{pendant_outlier_factor};
my $epaPendantMinThreshold = $EPA_NG_DEFAULT{pendant_minimum_threshold};
my $sampleQCFile = "";
#Two optional sample filters. Both are mechanisms only: what counts as a flagged
#sample, and which coverage thresholds apply, are decided by the caller. Both
#default to off so that callers which do not ask for them keep the behaviour
#they had before, and only a caller that sets a deliberate policy - such as
#strain_within.pl - turns them on.
#
#(1) Drop samples that -sampleQC marks unfit for the tree. Without -sampleQC
#this is inert.
my $excludeFlaggedSamples = 0;
#(2) Apply -GenesPerSpecies/-relativeNTFraction/-NTfiltCount as a removal rather
#than as backbone/placement routing. The taxon-aware selector otherwise keeps
#every sample holding a single informative site, and with -placeOnBackbone 0
#there is no later stage that would reconsider it. The non-taxon-aware prefilter
#has always removed on these same thresholds, so this makes the two consistent.
my $enforceSampleCoverage = 0;
my %sampleQCStatus;
my %flaggedExcluded;
my %POST_ALIGNMENT_QC_DEFAULT = (
	between_species_enabled => 0,
	within_species_enabled => 1,
	minimum_sequences => 3,
	minimum_occupancy => 0.35,
	relative_modified_z => 5.0,
	minimum_loci_for_relative => 8,
);
my $postAlignmentLocusQC;
my $postAlignmentMinSequences = $POST_ALIGNMENT_QC_DEFAULT{minimum_sequences};
my $postAlignmentMinOccupancy = $POST_ALIGNMENT_QC_DEFAULT{minimum_occupancy};
my $postAlignmentRelativeZ = $POST_ALIGNMENT_QC_DEFAULT{relative_modified_z};
my $postAlignmentMinLociRelative =
	$POST_ALIGNMENT_QC_DEFAULT{minimum_loci_for_relative};
my $postAlignmentDivergenceQC;
my %RATE_MERGE_DEFAULT = (
	enabled => 0,
	maximum_bins => 8,
	target_sites_per_bin => 30_000,
	minimum_loci_per_bin => 20,
	minimum_sites_per_bin => 20_000,
);
my $rateMergePartitions = $RATE_MERGE_DEFAULT{enabled};
my $rateMergePartitionsExplicit = 0;
my $rateMergeMaxBins = $RATE_MERGE_DEFAULT{maximum_bins};
my $rateMergeTargetSites = $RATE_MERGE_DEFAULT{target_sites_per_bin};
my $rateMergeMinLoci = $RATE_MERGE_DEFAULT{minimum_loci_per_bin};
my $rateMergeMinSites = $RATE_MERGE_DEFAULT{minimum_sites_per_bin};
my %TAXON_AWARE_DEFAULT = (
	enabled => 1,
	maximum_loci => 500,
	core_loci => 400,
	candidate_extra => 150,
	minimum_sequence_nt => 60,
	target_loci_per_sample => 25,
	target_nt_per_sample => 7500,
	rescue_minimum_prevalence => 0.8,
);
my $taxonAwareLocusSelection = $TAXON_AWARE_DEFAULT{enabled};
my $taxonAwareMaxLoci = $TAXON_AWARE_DEFAULT{maximum_loci};
my $taxonAwareCoreLoci = $TAXON_AWARE_DEFAULT{core_loci};
my $taxonAwareCandidateExtra = $TAXON_AWARE_DEFAULT{candidate_extra};
my $taxonAwareMinSequenceNT = $TAXON_AWARE_DEFAULT{minimum_sequence_nt};
my $taxonAwareTargetLoci = $TAXON_AWARE_DEFAULT{target_loci_per_sample};
my $taxonAwareTargetNT = $TAXON_AWARE_DEFAULT{target_nt_per_sample};
my $taxonAwareRescueMinPrevalence =
	$TAXON_AWARE_DEFAULT{rescue_minimum_prevalence};
#Scale-dependent shape of the final locus score. informationScore saturates at
#the saturation density; excessVariationPenalty ramps from onset over span. The
#between-species values are the historical ones. Within a species the two must
#move together: leaving saturation at 2% and onset at 20% opens a ten-fold band
#in which a locus collects full information credit and no penalty at all, which
#is exactly where cross-species read recruitment, unresolved paralogy and
#recombination sit while genuine core loci sit below the credit threshold.
my %TAXON_AWARE_SCORE_DEFAULT = (
	within_information_saturation => 0.005,
	between_information_saturation => 0.02,
	within_excess_variation_onset => 0.05,
	between_excess_variation_onset => 0.20,
	within_excess_variation_span => 0.10,
	between_excess_variation_span => 0.30,
);
my ($taxonAwareInformationSaturation, $taxonAwareExcessVariationOnset,
	$taxonAwareExcessVariationSpan);
#Rescue eligibility is measured on recovered prevalence, which is capped by
#sequencing depth: a locus in every genome cannot reach an absolute 0.8 once a
#fifth of the samples are shallow, and the coverage phase then selects nothing.
#'relative' reads -taxonAwareRescueMinPrevalence as a fraction of the highest
#prevalence any locus actually attains, so the gate tracks the achievable
#ceiling; 'absolute' restores the previous fixed threshold.
my $taxonAwareRescuePrevalenceMode = 'relative';
#Derive the per-sample coverage targets from the same Q90 basis that decides
#sample retention, so the greedy phase optimises toward the threshold that
#actually governs inclusion rather than a much lower fixed floor.
my $taxonAwareTargetsFromGate = 1;
my $preferredCoreGenes = "";
#The -preferredCoreGenes guide arrives in the presorter's own rank order, which
#encodes marker status, chimera/paralogy rejection, the MGS prevalence window
#and expected informative yield - evidence this script cannot re-derive from
#sequence lengths alone. Carry that order into both scoring stages as a weighted
#term rather than discarding it. 0 restores the previous rank-blind behaviour.
my $taxonAwarePresortWeight = 0.15;
my $compactTaxonAwareDiagnostics = 1;
my $preferredCoreGeneSet = {};
my $ntFiltExplicit = 0;


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


# Do not accept abbreviated switches: in particular, the retired -NTfilt
# spelling must not ambiguously match -NTfiltCount/-NTfiltPerGene.
Configure('no_auto_abbrev');
GetOptions(
	"genoInD=s" => \$genoindir, #provide a dir with complete genomes, will extract FGMs and build tree between genomes (NT/AA flag)
	"wildcardflag=s" => \$wildcardflag,
	"fna=s" => \$fnFna,
	"aa=s"      => \$aaFna,
	"cats=s"      => \$cogCats,
	"outD=s"      => \$outD,
	"tmpD=s" => \$tmpD,
	"stagedInputDir=s" => \$stagedInputDir,
	"tmpSubdir=s" => \$tmpSubdir,
	"completionMarker=s" => \$completionMarker,
	"onlyMSA=i" => \$onlyMSA,
	"terminalMarker=s" => \$terminalMarker,
	"placementPendingMarker=s" => \$placementPendingMarker,
	"withinSpecies=i" => \$withinSpecies,
	"strainWithinPreset=i" => \$strainWithinPreset,
	"cores=i" => \$ncore,
	"superTree=i" => \$doSuperTree,
	"superCheck=i" => \$doSuperCheck,
	"fixHeaders=i" => \$fixHeaders, ## fix the fasta headers, if too long or containing not allowed symbols (nwk reserved)
	"useEte=i"      => \$Ete,
	"relativeNTFraction=f" => \$ntFrac,
	"NTfiltPerGene=f"      => \$ntFracGene,
	"GeneLengthIncludeMin=f" => \$ntFracGeneInclude,
	"GenesPerSpecies=f" => \$GeneFracPSpec,
	"placementGenesPerSpecies=f" => \$placementGeneFracPSpec,
	"placementRelativeNTFraction=f" => \$placementNTFrac,
	"placementNTfiltCount=i" => \$placementNTCntTotal,
	"fracMaxGenes90pct=f" => \$fracMaxGenes90pct,
	"NTfiltCount=i" => \$ntCntTotal,
	"smplDef=i"	=> \$smplDef, #is the genome somehow quantified with a delimiter (_) ?
	"smplSep=s" => \$smplSep, #set the delimiter
	"outgroup=s"	=> \$outgroup,
	"AAtree=i" => \$useAA4tree,
	"MSAprogram=i" => \$MSAprog, #(0) MSAprobs, (1) clustalO, (2) mafft, (4) MUSCLE5, (5) FAMSA2 (only AA)
	"MSAfixRecoverTechnicalOffsets=i" => \$msaFixRecoverTechnicalOffsets, #repair coding-NT gap offsets after back-translation
	"MSAfixCodingFrame=i" => \$msaFixCodingFrame,
	"MSAfixGeneticCode=i" => \$msaFixGeneticCode,
	"MSAfixRecoveryBand=i" => \$msaFixRecoveryBand,
	"minOverlapMSA=f" => \$minOverlapMSA, #minimum called-sequence fraction per retained MSA column
	"maxGapPerCol=f" =>\$maxGapPerCol, #same as minOverlapMSA, but for MSAfix and %of gaps allowed in a column
	"calcDistMat=i" => \$calcDistMat,
	"calcDistMatExt=i" => \$calcDistMatExt,
	"calcDiffDNA=i" => \$calcDNAdiff,
	"minPcId=f" => \$minPcId, #sequence is filtered from data, unless the average minPcId is >= $minPcId
	"SynTree=i"	=> \$calcSyn,
	"NonSynTree=i"	=> \$calcNonSyn,
	"continue=i" => \$continue,
	"epaOnly=i" => \$epaOnly,
	"redoEPAfilter:i" => sub { $redoEPAfilter = $_[1] || 1; },
	"bootstrap=i" => \$bootStrap,
	"subsetSmpls=i" => \$subsetSmpls,
	"postFilter=s" => \$postFilter, # "," sep list of zorro,guidance2,macse
	"rmMSA=i" => \$removeMSA, #1 removes per-locus MSAs to save diskspace, 0 retains them for resume
	"gzInput=i" => \$gzipInput, #to save diskspace
	"isAligned=i" => \$isAligned,
	"runRAxML=i" => \$doRAXML,
	"runRaxMLng=i" => \$doRAXMLng,
	"runFastTree=i" => \$doFastTree,
	"runVeryFastTree=i" => \$doVeryFastTree,
	"treeShrink=i" => \$useTreeShrink,
	"sampleQC=s" => \$sampleQCFile,
	"excludeFlaggedSamples=i" => \$excludeFlaggedSamples,
	"enforceSampleCoverage=i" => \$enforceSampleCoverage,
	"placeOnBackbone=i" => sub { $strictBackbone = $_[1]; $placeOnBackboneSpecified = 1; },
	"strictBackbone=i" => sub { $strictBackbone = $_[1]; $legacyStrictBackboneSpecified = 1; },
	"strictBackboneFraction=f" => \$strictBackboneFraction,
	"placementMinOverlap=i" => \$placementMinOverlap,
	"strictBackboneMinSamples=i" => \$strictBackboneMinSamples,
	"epaThreads=i" => \$epaThreads,
	"epaMaxMemMB=i" => \$epaMaxMemMB,
	"epaPendantOutlierFactor=f" => \$epaPendantOutlierFactor,
	"epaPendantMinThreshold=f" => \$epaPendantMinThreshold,
	"postAlignmentLocusQC=i" => \$postAlignmentLocusQC,
	"postAlignmentMinSequences=i" => \$postAlignmentMinSequences,
	"postAlignmentMinOccupancy=f" => \$postAlignmentMinOccupancy,
	"postAlignmentDivergenceQC=i" => \$postAlignmentDivergenceQC,
	"postAlignmentRelativeZ=f" => \$postAlignmentRelativeZ,
	"postAlignmentMinLociRelative=i" => \$postAlignmentMinLociRelative,
	"rateMergePartitions=i" => sub {
		$rateMergePartitions = $_[1];
		$rateMergePartitionsExplicit = 1;
	},
	"rateMergeMaxBins=i" => \$rateMergeMaxBins,
	"rateMergeTargetSites=i" => \$rateMergeTargetSites,
	"rateMergeMinLoci=i" => \$rateMergeMinLoci,
	"rateMergeMinSites=i" => \$rateMergeMinSites,
	"taxonAwareLocusSelection=i" => \$taxonAwareLocusSelection,
	"taxonAwareMaxLoci=i" => \$taxonAwareMaxLoci,
	"taxonAwareCoreLoci=i" => \$taxonAwareCoreLoci,
	"taxonAwareCandidateExtra=i" => \$taxonAwareCandidateExtra,
	"taxonAwareMinSequenceNT=i" => \$taxonAwareMinSequenceNT,
	"taxonAwareTargetLoci=i" => \$taxonAwareTargetLoci,
	"taxonAwareTargetNT=i" => \$taxonAwareTargetNT,
	"preferredCoreGenes=s" => \$preferredCoreGenes,
	"taxonAwarePresortWeight=f" => \$taxonAwarePresortWeight,
	"compactTaxonAwareDiagnostics=i" => \$compactTaxonAwareDiagnostics,
	"taxonAwareRescueMinPrevalence=f" => \$taxonAwareRescueMinPrevalence,
	"taxonAwareRescuePrevalenceMode=s" => \$taxonAwareRescuePrevalenceMode,
	"taxonAwareTargetsFromGate=i" => \$taxonAwareTargetsFromGate,
	"taxonAwareInformationSaturation=f" => \$taxonAwareInformationSaturation,
	"taxonAwareExcessVariationOnset=f" => \$taxonAwareExcessVariationOnset,
	"taxonAwareExcessVariationSpan=f" => \$taxonAwareExcessVariationSpan,
	"runIQtree=i" => \$doIQTree,
	"AutoModel=i" => sub {
		$treeAutoModel = $_[1];
		$treeAutoModelExplicit = 1;
	},
	"iqFast=i" => \$iqFast, #fast qiTree mode
	"iqMemMB=i" => \$iqMemMB, #IQ-TREE RAM cap in MB; 0 leaves IQ-TREE uncapped
	"iqPathogen=i" => \$iqPathogen, #IQ-TREE 3 CMAPLE/native low-divergence selection
	"iqLegacy=i" => sub {
		$iqLegacy = $_[1];
		$iqLegacyExplicit = 1;
	}, #restore the pre-5.14 IQ-TREE command
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
die "-placeOnBackbone cannot be combined with deprecated -strictBackbone\n"
	if $placeOnBackboneSpecified && $legacyStrictBackboneSpecified;
warn "Option -strictBackbone is deprecated; use -placeOnBackbone instead. "
	."Compatibility support will be removed in a future release.\n"
	if $legacyStrictBackboneSpecified;
die "-withinSpecies must be 0 or 1\n"
	unless $withinSpecies == 0 || $withinSpecies == 1;
die "-strainWithinPreset must be 0 or 1\n"
	unless $strainWithinPreset == 0 || $strainWithinPreset == 1;
die "-epaOnly must be 0 or 1\n" unless $epaOnly == 0 || $epaOnly == 1;
die "-AutoModel must be 0 or 1\n"
	unless $treeAutoModel == 0 || $treeAutoModel == 1;
die "-redoEPAfilter must be 0 or 1\n"
	unless $redoEPAfilter == 0 || $redoEPAfilter == 1;
die "-onlyMSA must be 0 or 1\n" unless $onlyMSA == 0 || $onlyMSA == 1;
die "-onlyMSA 1 cannot be combined with -placeOnBackbone 1\n"
	if $onlyMSA && $strictBackbone;

$minOverlapMSA = $withinSpecies ? 0.35 : 0 unless defined $minOverlapMSA;
$postAlignmentLocusQC = $withinSpecies
	? $POST_ALIGNMENT_QC_DEFAULT{within_species_enabled}
	: $POST_ALIGNMENT_QC_DEFAULT{between_species_enabled}
	unless defined $postAlignmentLocusQC;
$postAlignmentDivergenceQC = $withinSpecies ? 1 : 0
	unless defined $postAlignmentDivergenceQC;
$taxonAwareInformationSaturation = $TAXON_AWARE_SCORE_DEFAULT{
	($withinSpecies ? 'within' : 'between').'_information_saturation'}
	unless defined $taxonAwareInformationSaturation;
$taxonAwareExcessVariationOnset = $TAXON_AWARE_SCORE_DEFAULT{
	($withinSpecies ? 'within' : 'between').'_excess_variation_onset'}
	unless defined $taxonAwareExcessVariationOnset;
$taxonAwareExcessVariationSpan = $TAXON_AWARE_SCORE_DEFAULT{
	($withinSpecies ? 'within' : 'between').'_excess_variation_span'}
	unless defined $taxonAwareExcessVariationSpan;
die "Unexpected positional arguments: @ARGV\n" if @ARGV;

die "-cores must be a positive integer\n" if $ncore < 1;
die "-rmMSA must be 0 or 1\n" unless $removeMSA == 0 || $removeMSA == 1;
die "-bootstrap must be zero or greater\n" if $bootStrap < 0;
die "-NTfiltCount must be zero or greater\n" if $ntCntTotal < 0;
$placementGeneFracPSpec = $GeneFracPSpec unless defined $placementGeneFracPSpec;
$placementNTFrac = $ntFrac unless defined $placementNTFrac;
$placementNTCntTotal = $ntCntTotal unless defined $placementNTCntTotal;
die "-placementNTfiltCount must be zero or greater\n"
	if $strictBackbone && $placementNTCntTotal < 0;
die "-minOverlapMSA must be between zero and one\n"
	if $minOverlapMSA < 0 || $minOverlapMSA > 1;
die "-MSAfixRecoverTechnicalOffsets must be 0 or 1\n"
	unless $msaFixRecoverTechnicalOffsets == 0 || $msaFixRecoverTechnicalOffsets == 1;
die "-MSAfixCodingFrame must be 1, 2, or 3\n"
	unless $msaFixCodingFrame >= 1 && $msaFixCodingFrame <= 3;
die "-MSAfixGeneticCode must be a positive NCBI translation-table ID\n"
	if $msaFixGeneticCode < 1;
die "-MSAfixRecoveryBand must be zero or greater\n"
	if $msaFixRecoveryBand < 0;
die "-iqMemMB must be zero or greater\n" if $iqMemMB < 0;
die "-iqPathogen must be 0 or 1\n" unless $iqPathogen == 0 || $iqPathogen == 1;
die "-iqLegacy must be 0 or 1\n" unless $iqLegacy == 0 || $iqLegacy == 1;
die "-iqPathogen and -iqLegacy are mutually exclusive\n" if $iqPathogen && $iqLegacy;
die "-placeOnBackbone must be 0 or 1 (default $BACKBONE_DEFAULT{enabled})\n"
	unless $strictBackbone == 0 || $strictBackbone == 1;
if ($strictBackbone) {
	die "-strictBackboneFraction must be between 0 and 1 "
		."(default $BACKBONE_DEFAULT{coverage_fraction})\n"
		if $strictBackboneFraction < 0 || $strictBackboneFraction > 1;
	die "-placementMinOverlap must be non-negative "
		."(default $BACKBONE_DEFAULT{minimum_overlap})\n"
		if $placementMinOverlap < 0;
	die "-strictBackboneMinSamples must be at least 3 "
		."(default $BACKBONE_DEFAULT{minimum_samples})\n"
		if $strictBackboneMinSamples < 3;
	die "-epaThreads must be a positive integer (default $EPA_NG_DEFAULT{threads})\n"
		if $epaThreads < 1;
	die "-epaMaxMemMB must be -1 (derived), 0 (no memory-based scaling), or a positive MB value\n"
		if $epaMaxMemMB < -1;
	die "-epaPendantOutlierFactor and -epaPendantMinThreshold must be non-negative\n"
		if $epaPendantOutlierFactor < 0 || $epaPendantMinThreshold < 0;
}
if ($epaOnly) {
	die "-epaOnly requires -continue 1\n" unless $continue;
	die "-epaOnly requires -placeOnBackbone 1\n" unless $strictBackbone;
	die "-epaOnly currently requires exactly -runIQtree 1\n"
		unless $doIQTree && !$doRAXML && !$doRAXMLng
			&& !$doFastTree && !$doVeryFastTree;
	die "-epaOnly requires -completionMarker and -placementPendingMarker\n"
		unless length($completionMarker) && length($placementPendingMarker);
}
if ($redoEPAfilter) {
	die "-redoEPAfilter requires -continue 1\n" unless $continue;
	die "-redoEPAfilter requires -placeOnBackbone 1\n" unless $strictBackbone;
	die "-redoEPAfilter currently requires exactly -runIQtree 1\n"
		unless $doIQTree && !$doRAXML && !$doRAXMLng
			&& !$doFastTree && !$doVeryFastTree;
	die "-redoEPAfilter cannot be combined with -epaOnly\n"
		if $epaOnly;
	die "-redoEPAfilter requires -completionMarker\n"
		unless length($completionMarker);
}
die "-postAlignmentLocusQC must be 0 or 1 "
	."(default: 0 between species, 1 within species)\n"
	unless $postAlignmentLocusQC == 0 || $postAlignmentLocusQC == 1;
die "-postAlignmentMinSequences must be at least 2 "
	."(default $POST_ALIGNMENT_QC_DEFAULT{minimum_sequences})\n"
	if $postAlignmentMinSequences < 2;
die "-postAlignmentMinOccupancy must be between 0 and 1 "
	."(default $POST_ALIGNMENT_QC_DEFAULT{minimum_occupancy})\n"
	if $postAlignmentMinOccupancy < 0 || $postAlignmentMinOccupancy > 1;
die "-postAlignmentDivergenceQC must be 0 or 1 "
	."(default: 0 between species, 1 within species)\n"
	unless $postAlignmentDivergenceQC == 0 || $postAlignmentDivergenceQC == 1;
die "-postAlignmentRelativeZ must be non-negative "
	."(default $POST_ALIGNMENT_QC_DEFAULT{relative_modified_z})\n"
	if $postAlignmentRelativeZ < 0;
die "-postAlignmentMinLociRelative must be positive "
	."(default $POST_ALIGNMENT_QC_DEFAULT{minimum_loci_for_relative})\n"
	if $postAlignmentMinLociRelative < 1;
die "-rateMergePartitions must be 0 or 1\n"
	unless $rateMergePartitions == 0 || $rateMergePartitions == 1;
die "-rateMergeMaxBins, -rateMergeTargetSites, -rateMergeMinLoci, and -rateMergeMinSites must be positive\n"
	if grep { $_ < 1 } ($rateMergeMaxBins, $rateMergeTargetSites, $rateMergeMinLoci, $rateMergeMinSites);
die "-rateMergePartitions requires nucleotide trees (-AAtree 0)\n"
	if $rateMergePartitions && $useAA4tree;
die "-rateMergePartitions requires -postAlignmentLocusQC 1 because MSAfix "
	."supplies the overlap-filtered rate/GC metrics\n"
	if $rateMergePartitions && !$postAlignmentLocusQC;
die "-taxonAwareLocusSelection must be 0 or 1\n"
	unless $taxonAwareLocusSelection == 0 || $taxonAwareLocusSelection == 1;
die "-taxonAwareMaxLoci, -taxonAwareCoreLoci, -taxonAwareMinSequenceNT, "
	."-taxonAwareTargetLoci, and -taxonAwareTargetNT must be positive\n"
	if grep { $_ < 1 } ($taxonAwareMaxLoci, $taxonAwareCoreLoci,
		$taxonAwareMinSequenceNT, $taxonAwareTargetLoci, $taxonAwareTargetNT);
die "-taxonAwareCandidateExtra must be zero or greater\n"
	if $taxonAwareCandidateExtra < 0;
die "-taxonAwareCoreLoci cannot exceed -taxonAwareMaxLoci\n"
	if $taxonAwareCoreLoci > $taxonAwareMaxLoci;
die "-taxonAwareRescueMinPrevalence must be between 0 and 1 "
	."(default $TAXON_AWARE_DEFAULT{rescue_minimum_prevalence})\n"
	if $taxonAwareRescueMinPrevalence < 0
		|| $taxonAwareRescueMinPrevalence > 1;
die "-taxonAwareRescuePrevalenceMode must be 'relative' or 'absolute'\n"
	unless $taxonAwareRescuePrevalenceMode eq 'relative'
		|| $taxonAwareRescuePrevalenceMode eq 'absolute';
die "-taxonAwareTargetsFromGate must be 0 or 1\n"
	unless $taxonAwareTargetsFromGate == 0 || $taxonAwareTargetsFromGate == 1;
die "-taxonAwareInformationSaturation must be greater than 0 and at most 1\n"
	if $taxonAwareInformationSaturation <= 0 || $taxonAwareInformationSaturation > 1;
die "-taxonAwareExcessVariationOnset must be greater than 0 and at most 1\n"
	if $taxonAwareExcessVariationOnset <= 0 || $taxonAwareExcessVariationOnset > 1;
die "-taxonAwareExcessVariationSpan must be greater than 0\n"
	if $taxonAwareExcessVariationSpan <= 0;
die "-taxonAwareExcessVariationOnset must not sit below "
	."-taxonAwareInformationSaturation, or a locus would be penalised before it "
	."is fully credited\n"
	if $taxonAwareExcessVariationOnset < $taxonAwareInformationSaturation;
die "-taxonAwarePresortWeight must be between 0 and 1\n"
	if $taxonAwarePresortWeight < 0 || $taxonAwarePresortWeight > 1;
die "-compactTaxonAwareDiagnostics must be 0 or 1\n"
	unless $compactTaxonAwareDiagnostics == 0 || $compactTaxonAwareDiagnostics == 1;
if (length($preferredCoreGenes)) {
	$preferredCoreGenes = abs_path($preferredCoreGenes)
		or die "Cannot resolve -preferredCoreGenes: $preferredCoreGenes\n";
	die "-preferredCoreGenes is missing or empty: $preferredCoreGenes\n"
		unless -s $preferredCoreGenes;
	$preferredCoreGeneSet = preferredCoreGeneSet($preferredCoreGenes);
}
die "-tmpD and -tmpSubdir are mutually exclusive\n" if length($tmpD) && length($tmpSubdir);
if (length($tmpSubdir)) {
	die "-tmpSubdir must be a safe relative path\n"
		if File::Spec->file_name_is_absolute($tmpSubdir)
			|| grep { $_ eq File::Spec->updir } File::Spec->splitdir($tmpSubdir);
}
die "-smplSep must not be empty\n" if $smplSep eq "";
eval { qr/$smplSep/ } or die "Invalid -smplSep regular expression '$smplSep': $@";
compileSampleSeparator();
my @fractionNameValues = (
	["relativeNTFraction", $ntFrac],
	["NTfiltPerGene", $ntFracGene],
	["GeneLengthIncludeMin", $ntFracGeneInclude],
	["GenesPerSpecies", $GeneFracPSpec],
	["fracMaxGenes90pct", $fracMaxGenes90pct],
	["maxGapPerCol", $maxGapPerCol],
);
push @fractionNameValues,
	["placementRelativeNTFraction", $placementNTFrac],
	["placementGenesPerSpecies", $placementGeneFracPSpec]
	if $strictBackbone;
for my $fraction_name_value (@fractionNameValues) {
	my ($name, $value) = @{$fraction_name_value};
	die "-$name must be between 0 and 1\n" if $value < 0 || $value > 1;
}
die "-GeneLengthIncludeMin cannot exceed -NTfiltPerGene because the inclusion "
	."threshold must not be stricter than the QC threshold\n"
	if $ntFracGeneInclude > $ntFracGene;
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
$terminalMarker ||= File::Spec->catfile($outD, 'noTree.sto');
$placementPendingMarker ||= File::Spec->catfile($outD, 'placementPending.sto');
$workflowStateFile = File::Spec->catfile($outD, q{buildTree.state.tsv});
writeWorkflowHeartbeat('preflight');

if (length($stagedInputDir)) {
	my @requiredInputs = grep { defined($_) && length($_) } ($fnFna, $aaFna, $cogCats);
	publishStagedTreeInputs($stagedInputDir, $outD, $ncore, \@requiredInputs,
		$sampleQCFile, $outgroup);
}
die "-sampleQC does not exist or is empty: $sampleQCFile\n"
	if length($sampleQCFile) && !fileGZs($sampleQCFile);
die "-excludeFlaggedSamples must be 0 or 1\n"
	unless $excludeFlaggedSamples == 0 || $excludeFlaggedSamples == 1;
die "-enforceSampleCoverage must be 0 or 1\n"
	unless $enforceSampleCoverage == 0 || $enforceSampleCoverage == 1;
if (length($sampleQCFile)) {
	%sampleQCStatus = %{read_sample_qc($sampleQCFile)};
	if ($excludeFlaggedSamples) {
		# The outgroup is a catalogue reference rather than an extracted sample, so
		# it carries no QC verdict; never let a stray one remove the root either.
		%flaggedExcluded = map { $_ => 1 }
			grep { !length($outgroup) || $_ ne $outgroup }
			grep { ($sampleQCStatus{$_} // '') eq 'mixed_strain' } keys %sampleQCStatus;
	}
}

##### setup dirs
$codemlOutD = File::Spec->catdir($outD, "codeml") if ($codemlOutD eq "");

my $requestedTmpBase;
if (length($tmpSubdir)) {
	my $environmentTmp = $ENV{TMPDIR} // "";
	my $temporaryRoot = length($environmentTmp)
		? $environmentTmp
		: File::Spec->catdir($outD, "tmp");
	$requestedTmpBase = File::Spec->canonpath(
		File::Spec->catdir($temporaryRoot, File::Spec->splitdir($tmpSubdir))
	);
} else {
	$requestedTmpBase = $tmpD eq ""
		? File::Spec->catdir($outD, "tmp")
		: File::Spec->canonpath(File::Spec->rel2abs($tmpD));
}
my $tmpBase = $requestedTmpBase;
my ($tmpReady, $tmpError) = prepareTemporaryBase($tmpBase);
if (!$tmpReady && (length($tmpD) || length($tmpSubdir))) {
	my $fallbackTmpBase = File::Spec->catdir($outD, "tmp");
	warn "Requested temporary path is unusable: $tmpBase ($tmpError); "
		. "falling back to $fallbackTmpBase\n";
	$tmpBase = $fallbackTmpBase;
	($tmpReady, $tmpError) = prepareTemporaryBase($tmpBase);
}
die "No usable temporary path is available: $tmpBase ($tmpError)\n"
	unless $tmpReady;
my $tmpTag = $clusterName eq "" ? "default" : $clusterName;
$tmpTag =~ s/[^A-Za-z0-9_.-]+/_/g;
$tmpD = File::Spec->catdir($tmpBase, "buildTree5_${tmpTag}_$$");
make_path($tmpD);
my $treeD = File::Spec->catdir($outD, "phylo");#raxml, fasttree, phyml tree output dir

my $MsaD = File::Spec->catdir($outD, "MSA");
# Per-locus checkpoints remain available until terminal workflow finalization.
if ($subsetSmpls >0){
	# catdir never returns a trailing separator, so these directories have to be
	# renamed by appending the subset tag; a trailing-slash substitution would be
	# a silent no-op and let a subset run overwrite the full run's outputs.
	$MsaD .= "_S$subsetSmpls";
	$treeD .= "_S$subsetSmpls";
}
my $MsaWorkD = File::Spec->catdir($tmpD, basename($MsaD));
make_path($MsaWorkD) unless -d $MsaWorkD;
$MSAsubsD = File::Spec->catdir($MsaWorkD, 'clnd');

my $multAliArtifact = File::Spec->catfile($MsaD, 'MSAli.fna');
my $multAli = File::Spec->catfile($MsaWorkD, 'MSAli.fna');
my $multAliSynArtifact = $multAliArtifact.'.syn.fna';
my $multAliSyn = $multAli.'.syn.fna';
my $multAliNonSynArtifact = $multAliArtifact.'.nonsyn.fna';
my $multAliNonSyn = $multAli.'.nonsyn.fna';
my $placementAlignmentArtifact = File::Spec->catfile($MsaD, 'MSAli.placement.fna');
my $placementAlignment = File::Spec->catfile($MsaWorkD, 'MSAli.placement.fna');
if ($redoEPAfilter) {
	runRedoEpaFilter(
		$treeD, File::Spec->catfile($treeD, 'strict_backbone.samples.tsv'));
	exit(0);
}

######

warn "MSAprobs may emit non-fatal trimming warnings\n" if $MSAprog == 0;
if ($doCFML && !$doRAXML){die "Need RaxML alignment, if Clonal fram is to be run..\n";}

if ($aaFna eq "" || $useAA4tree){	$calcSyn=0;$calcNonSyn=0;}

$MSAreq = 0 if (!$doFastTree && !$doVeryFastTree && !$doRAXML && !$doRAXMLng && !$doCFML && !$doGubbins && !$doIQTree);

make_path($tmpD) unless -d $tmpD;
my %msaProgramNames = (
	0 => "MSAprobs",
	1 => "Clustal Omega",
	2 => "MAFFT",
	4 => "MUSCLE5",
	5 => "FAMSA2",
);
my @treeMethods;
push @treeMethods, "IQ-TREE" if $doIQTree;
push @treeMethods, "RAxML" if $doRAXML;
push @treeMethods, "RAxML-NG" if $doRAXMLng;
push @treeMethods, "FastTree" if $doFastTree;
push @treeMethods, "VeryFastTree" if $doVeryFastTree;
push @treeMethods, "ClonalFrameML" if $doCFML;
push @treeMethods, "Gubbins" if $doGubbins;
my @inputDescriptions;
push @inputDescriptions, "NT=$fnFna" if $fnFna ne "";
push @inputDescriptions, "AA=$aaFna" if $aaFna ne "";
push @inputDescriptions, "categories=$cogCats" if $cogCats ne "";
push @inputDescriptions, "genomes=$genoindir" if $genoindir ne "";
push @inputDescriptions, "preferred-core=$preferredCoreGenes" if $preferredCoreGenes ne "";
print "=====================================================\n";
print "BuildTree pipeline v$version\n";
print "Inputs: " . join("; ", @inputDescriptions) . "\n";
print "Paths: output=$outD; temporary=$tmpD; MSA work=$MsaWorkD; MSA checkpoints=$MsaD; trees=$treeD\n";
print "Mode: " . ($cogCats ne "" ? "multi-locus" : "single-locus")
	. "; scope=" . ($withinSpecies ? "within-species" : "between-species/broad")
	. "; sequence=" . ($useAA4tree ? "amino acid" : "nucleotide")
	. "; input aligned=" . ($isAligned ? "yes" : "no")
	. "; continue=" . ($continue ? "yes" : "no") . "\n";
print "Alignment: $msaProgramNames{$MSAprog}; cores=$ncore; post-filter="
	. ($postFilter || "<none>")
	. ($removeMSA ? "; per-locus MSA cleanup=enabled" : "; per-locus nucleotide and protein MSAs retained (resumable)")
	. "; MSAli compression=always\n";
print "MSAfix coding-NT repair: " . ($msaFixRecoverTechnicalOffsets ? "enabled; frame=$msaFixCodingFrame; genetic code=$msaFixGeneticCode; band=$msaFixRecoveryBand" : "disabled") . "\n"
	unless $useAA4tree;
print "Filtering: per-gene QC length fraction=$ntFracGene; post-QC inclusion fraction=$ntFracGeneInclude; category Q90 fraction=$fracMaxGenes90pct; "
	. "backbone NT fraction=$ntFrac; backbone gene fraction=$GeneFracPSpec; "
	. "backbone minimum NT=$ntCntTotal; minimum overlap=$minOverlapMSA; maximum gap fraction=$maxGapPerCol\n";
print "Flagged-sample exclusion: "
	. (!length($sampleQCFile) ? "inactive (no -sampleQC table supplied)"
		: !$excludeFlaggedSamples ? "disabled by -excludeFlaggedSamples 0"
		: "enabled; ".scalar(keys %flaggedExcluded)." of "
			.scalar(keys %sampleQCStatus)." sample(s) marked unfit by the caller")
	. "\n";
print "Sample coverage filter: "
	. ($enforceSampleCoverage
		? "enabled; samples below the thresholds above are removed"
		: "disabled by -enforceSampleCoverage 0; every aligned sample is retained")
	. ($strictBackbone ? " (superseded by -placeOnBackbone routing)" : "") . "\n";
print "Post-alignment locus QC: enabled="
	. ($postAlignmentLocusQC ? "yes" : "no")
	. "; divergence QC=" . ($postAlignmentDivergenceQC ? "yes" : "no")
	. "; minimum sequences=$postAlignmentMinSequences"
	. "; minimum occupancy=$postAlignmentMinOccupancy"
	. "; relative modified-Z="
	. ($postAlignmentDivergenceQC ? $postAlignmentRelativeZ : "<disabled>")
	. "; minimum loci for relative QC=$postAlignmentMinLociRelative\n";
print "Partition merging: enabled=" . ($rateMergePartitions ? "yes" : "no")
	. "; maximum bins=$rateMergeMaxBins"
	. "; target size=$rateMergeTargetSites effective sites/bin"
	. "; minimum bin size=$rateMergeMinLoci loci/$rateMergeMinSites sites\n";
print "Taxon-aware locus selection: enabled="
	. ($taxonAwareLocusSelection ? "yes" : "no")
	. "; final loci=$taxonAwareMaxLoci; robust core=$taxonAwareCoreLoci"
	. "; alignment backfill=$taxonAwareCandidateExtra"
	. "; minimum sequence NT=$taxonAwareMinSequenceNT"
	. "; rescue minimum prevalence=$taxonAwareRescueMinPrevalence"
	. " ($taxonAwareRescuePrevalenceMode)"
	. "; targets from retention gate=" . ($taxonAwareTargetsFromGate ? "yes" : "no")
	. "; presort-rank weight=$taxonAwarePresortWeight"
	. "; information saturation=$taxonAwareInformationSaturation"
	. "; excess-variation penalty from $taxonAwareExcessVariationOnset over $taxonAwareExcessVariationSpan"
	. "; preferred universal core="
	. (keys(%{$preferredCoreGeneSet}) ? scalar(keys %{$preferredCoreGeneSet}) : "<none>")
	. "; sample target=$taxonAwareTargetLoci loci/$taxonAwareTargetNT NT\n";
if ($strictBackbone) {
	print "Backbone placement: enabled; sample QC=" . ($sampleQCFile || "<none>")
		. "; coverage fraction=$strictBackboneFraction"
		. "; placement gene fraction=$placementGeneFracPSpec"
		. "; placement NT fraction=$placementNTFrac"
		. "; placement minimum NT=$placementNTCntTotal"
		. "; minimum placement overlap=$placementMinOverlap"
		. "; minimum backbone samples=$strictBackboneMinSamples\n";
	my ($epaReportedThreads, $epaReportedMaxMemMB) =
		epaResourcePlan($epaThreads, $ncore, $epaMaxMemMB, $iqMemMB);
	print "EPA-ng placement resources: threads=$epaReportedThreads"
		." (requested=$epaThreads); planning memory="
		.($epaReportedMaxMemMB ? "${epaReportedMaxMemMB}MB at $EPA_NG_DEFAULT{memory_per_thread_mb}MB/thread" : "disabled")
		."; hard memory limit=scheduler/cgroup; query chunk=$EPA_NG_DEFAULT{chunk_size}\n";
	print "EPA-ng placement outlier QC: "
		.($epaPendantOutlierFactor > 0
			? "maximum pendant branch=max($epaPendantMinThreshold, "
				."$epaPendantOutlierFactor x backbone terminal-branch Q95)"
			: "disabled") . "\n";
} else {
	print "Backbone placement: disabled; backbone- and placement-only filters are inactive\n";
}
print "Trees: " . (@treeMethods ? join(", ", @treeMethods) : "<none>")
	. "; bootstrap=$bootStrap; outgroup=" . ($outgroup || "<none>")
	. "; supertree=" . ($doSuperTree ? "yes" : "no")
	. "; IQ-TREE mode=" . ($iqLegacy ? "legacy" : $iqPathogen ? "pathogen" : "standard")
	. "; IQ-TREE model=" . ($treeAutoModel ? "AutoModel" : $useAA4tree ? "LG+F+G" : "GTR+F+G2")
	. "; IQ-TREE memory=" . ($iqMemMB ? "${iqMemMB}MB" : "auto") . "\n";
print "Additional analyses: synonymous=" . ($calcSyn ? "yes" : "no")
	. "; nonsynonymous=" . ($calcNonSyn ? "yes" : "no")
	. "; distance matrix=" . ($calcDistMat ? "yes" : "no")
	. "; dN/dS=" . ($doDNDS ? "yes" : "no")
	. "; TreeShrink=" . ($useTreeShrink ? "yes" : "no") . "\n";
print "=====================================================\n";
my $cmd =""; my %usedGeneNms; my %excludedLoci;


my $outD_clust = File::Spec->catdir($tmpD, "fastGear_work_$tmpTag");


#------------------------------------------
#sorting by COG, MSA & syn position extraction
if ($Ete){
	my $eteBin = getProgPaths("ete3");
	my $eteOut = File::Spec->catdir($outD, "tree");
	make_path($eteOut);
	$cmd = "$eteBin build -n $fnFna -a $aaFna -w clustalo_default-none-none-none  -m sptree_raxml_all --cpu $ncore -o $eteOut --clearall --nt-switch 0.0 --noimg";
	$cmd .= " --cogs $cogCats" unless ($cogCats eq "");
	print "Running ETE tree analysis; detailed output: $eteOut/ETE.log\n";
	systemW($cmd . " > $eteOut/ETE.log 2>&1");
	print "ETE tree analysis completed; output: $eteOut\n";
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
preflightBuildTree($outD, $tmpD);

if ($epaOnly) {
	restoreCompressedMSAArtifact($multAliArtifact, $multAli);
	restoreCompressedMSAArtifact($placementAlignmentArtifact, $placementAlignment);
	my $epaTreeOptions = createTreeOpt($multAli, 'allsites', '', 0, '');
	$epaTreeOptions->{IQtreeout} .= '.backbone';
	my $epaOnlyComplete = runEpaOnlyPlacement(
		$epaTreeOptions, $multAli, $placementAlignment, $treeD,
		File::Spec->catfile($treeD, 'strict_backbone.samples.tsv'),
	);
	finalizeMSAArtifacts($MsaD, $MsaWorkD) unless $epaOnlyComplete;
	exit(0);
}

if (!$continue){
	safeRemoveTree($treeD, $outD);
	safeRemoveTree($MsaD, $outD);
}
make_path($MsaD) unless -d $MsaD;
make_path($MsaWorkD) unless -d $MsaWorkD;
make_path($treeD) unless -d $treeD;
my @theRealMSAs;
my $partiFile="";#partitioning for multi gene MSAs
my %specList; #list of species (without _COG00012 tag);
my %samples; 
my $MSAcat = File::Spec->catfile($MsaWorkD, 'MSAcat.fna');

#prep tree Options
my $tOhr = createTreeOpt($multAli,"allsites","",0,"");
$tOhr->{IQtreeout} .= ".backbone" if $strictBackbone && $doIQTree;
my %Tree1 = %{$tOhr};
# The all-sites tree does not exist yet: createTreeOpt only sets 'nwk' once
# treeAtHeart has run. The constraint is attached immediately before the
# syn/nonsyn trees are inferred, not captured from this pre-inference snapshot.
my $tOhrNSun = createTreeOpt($multAliNonSyn,"nonsyn","",0,"");
my $tOhrSyn = createTreeOpt($multAliSyn,"syn","",0,"");



#DEBUG
#mergePids("$outD/MSA/",40, "NT") ;die;
#my $tmp = "/g/bork5/hildebra/results/TEC2/v5/T2dphylo/rDNA2/fullGenomes/ini16S.fna";
#calcDisPos2($tmp,"$outD/MSA/percID_syn.txt",1); die;


my @MSAs; my @MSA_AA; my @MSAsSyn; my @MSAsNonSyn;#full MSAs and MSAs with syn / nonsyn pos only
my @MSrm; 
my %FAA ; my %FNA ; my @geneList; my @geneListF;
my (%primaryAlignmentGene, %taxonAwarePreMetrics, %taxonAwareUniverseSamples);
my (%partitionRateProxy, %partitionSelectionPhase, %taxonAwareFinalMetricByPath);
my (%taxonAwareBackboneEligibility, %taxonAwareBackboneIneligibleReason);
my (%taxonAwarePlacementEligibility, %taxonAwarePlacementIneligibleReason);
my $strictSplit;
my $strictPlacementMinimumNT = $placementMinOverlap;
my $strictPlacementMinimumLoci = 2;
my $postAlignmentQCReport = "$treeD/post_alignment_locus_qc.tsv";
my $selectionAttritionReport = "$treeD/selection_attrition.tsv";
my $geneLengthSampleReport = "$treeD/gene_length_filter.samples.tsv";
my (%geneLengthSampleAudit, %geneLengthQCSequence, %geneLengthIncludeByGene);
my %geneLengthRecoveredForMSA;
my %selectionAttrition = map { $_ => 'NA' } qw(
	input_loci input_sequences input_samples
	qc_excluded_samples qc_excluded_sequences
	qc_emptied_loci coverage_excluded_samples length_retained_sequences
	length_filtered_sequences length_include_retained_sequences
	length_include_filtered_sequences length_recovery_candidate_sequences
	length_recovered_msa_sequences gene_length_min_dropped_loci
	gene_length_include_min_dropped_loci gene_length_recovery_candidate_loci
	gene_length_recovered_msa_loci eligible_loci candidate_loci candidate_samples
	aligned_loci alignment_failed_loci post_qc_loci final_loci final_samples
	backbone_samples placement_samples excluded_samples
);
#A run that never reaches the coverage filter removed nothing through it, which
#is a measurement rather than an absent one.
$selectionAttrition{coverage_excluded_samples} = 0;
my $legacyPostAlignmentQCPolicyFile = "$treeD/post_alignment_locus_qc.policy.tsv";
my $legacyAlignmentWorkPolicyFile = "$MsaD/alignment_work.policy.tsv";
my $legacyPostAlignmentPolicyFile = "$treeD/post_alignment.policy.tsv";
my $postAlignmentQCPolicy = join("\t",
	"schema=16",
	# Flagged-sample exclusion changes which sequences reach the alignment, so a
	# flipped setting must invalidate cached per-locus alignments and QC.
	"qc_sample_exclusion=".($excludeFlaggedSamples ? 1 : 0),
	"qc_excluded=".scalar(keys %flaggedExcluded),
	"msa_program=$MSAprog",
	"post_filter=$postFilter",
	"sample_definition=$smplDef",
	"sample_separator=$smplSep",
	"preferred_core_input=".inputFingerprint($preferredCoreGenes),
	"outgroup=$outgroup",
	"synonymous_sites=$calcSyn",
	"nonsynonymous_sites=$calcNonSyn",
	"fna_input=".inputFingerprint($fnFna),
	"faa_input=".inputFingerprint($aaFna),
	"category_input=".inputFingerprint($cogCats),
	"enabled=$postAlignmentLocusQC",
	"scope=".($withinSpecies ? "within" : "between"),
	"sequence=".($useAA4tree ? "aa" : "nt"),
	"per_gene_length_fraction=$ntFracGene",
	"per_gene_inclusion_fraction=$ntFracGeneInclude",
	"minimum_category_q90_fraction=$fracMaxGenes90pct",
	"backbone_nt_fraction=$ntFrac",
	"backbone_gene_fraction=$GeneFracPSpec",
	"backbone_minimum_nt=$ntCntTotal",
	"minimum_overlap=$minOverlapMSA",
	"maximum_gap_fraction=$maxGapPerCol",
	"msafix_recover_technical_offsets=$msaFixRecoverTechnicalOffsets",
	"msafix_coding_frame=$msaFixCodingFrame",
	"msafix_genetic_code=$msaFixGeneticCode",
	"msafix_recovery_band=$msaFixRecoveryBand",
	"minimum_sequences=$postAlignmentMinSequences",
	"minimum_occupancy=$postAlignmentMinOccupancy",
	"divergence_qc=$postAlignmentDivergenceQC",
	"relative_modified_z=".($postAlignmentDivergenceQC
		? $postAlignmentRelativeZ : "disabled"),
	"minimum_loci_relative=$postAlignmentMinLociRelative",
	"rate_partition_merge=$rateMergePartitions",
	"rate_partition_maximum_bins=$rateMergeMaxBins",
	"rate_partition_target_sites=$rateMergeTargetSites",
	"rate_partition_minimum_loci=$rateMergeMinLoci",
	"rate_partition_minimum_sites=$rateMergeMinSites",
	"taxon_aware=$taxonAwareLocusSelection",
	"taxon_aware_maximum_loci=$taxonAwareMaxLoci",
	"taxon_aware_core_loci=$taxonAwareCoreLoci",
	"taxon_aware_candidate_extra=$taxonAwareCandidateExtra",
	"taxon_aware_minimum_sequence_nt=$taxonAwareMinSequenceNT",
	"taxon_aware_target_loci=$taxonAwareTargetLoci",
	"taxon_aware_target_nt=$taxonAwareTargetNT",
	"taxon_aware_rescue_minimum_prevalence=$taxonAwareRescueMinPrevalence",
	"taxon_aware_rescue_prevalence_mode=$taxonAwareRescuePrevalenceMode",
	"taxon_aware_targets_from_gate=$taxonAwareTargetsFromGate",
	"taxon_aware_presort_weight=$taxonAwarePresortWeight",
	"taxon_aware_information_saturation=$taxonAwareInformationSaturation",
	"taxon_aware_excess_variation_onset=$taxonAwareExcessVariationOnset",
	"taxon_aware_excess_variation_span=$taxonAwareExcessVariationSpan",
)."\n";
my $postAlignmentPolicy = join("\t",
	"schema=1",
	"only_msa=$onlyMSA",
	"tree_methods=iqtree:$doIQTree,raxml:$doRAXML,raxmlng:$doRAXMLng,fasttree:$doFastTree,veryfasttree:$doVeryFastTree",
	"bootstrap=$bootStrap",
	"iqtree_auto_model=$treeAutoModel",
	"iqtree_fast=$iqFast",
	"iqtree_memory_mb=$iqMemMB",
	"iqtree_pathogen=$iqPathogen",
	"iqtree_legacy=$iqLegacy",
	"epa_threads=".($strictBackbone ? $epaThreads : "disabled"),
	"epa_memory_mb=".($strictBackbone ? $epaMaxMemMB : "disabled"),
	"epa_pendant_outlier_factor=".($strictBackbone ? $epaPendantOutlierFactor : "disabled"),
	"epa_pendant_minimum_threshold=".($strictBackbone ? $epaPendantMinThreshold : "disabled"),
	"tree_shrink=$useTreeShrink",
	"clonal_frame=$doCFML",
	"gubbins=$doGubbins",
	"dna_distance=$calcDNAdiff",
	"strict_backbone=$strictBackbone",
	"strict_backbone_fraction=".($strictBackbone ? $strictBackboneFraction : "disabled"),
	"strict_backbone_minimum_samples=".($strictBackbone ? $strictBackboneMinSamples : "disabled"),
	"placement_nt_fraction=".($strictBackbone ? $placementNTFrac : "disabled"),
	"placement_gene_fraction=".($strictBackbone ? $placementGeneFracPSpec : "disabled"),
	"placement_minimum_nt=".($strictBackbone ? $placementNTCntTotal : "disabled"),
	"placement_minimum_overlap=".($strictBackbone ? $placementMinOverlap : "disabled"),
)."\n";
$workflowMsaSelectionPolicy = $postAlignmentQCPolicy;
$workflowTreeStagePolicy = $postAlignmentPolicy;
my $buildTreeState = readBuildTreeState($workflowStateFile);
my $stateHasPolicies = defined($buildTreeState->{msa_selection_policy})
	&& length($buildTreeState->{msa_selection_policy})
	&& defined($buildTreeState->{tree_stage_policy})
	&& length($buildTreeState->{tree_stage_policy});
my $legacyMsaSelectionPolicyMatches = legacyPolicyFileMatches(
	$legacyAlignmentWorkPolicyFile, $postAlignmentQCPolicy,
	"legacy alignment-work policy")
	|| legacyPolicyFileMatches($legacyPostAlignmentQCPolicyFile,
		$postAlignmentQCPolicy, "legacy locus-QC policy");
my $legacyTreeStagePolicyMatches = legacyPolicyFileMatches(
	$legacyPostAlignmentPolicyFile, $postAlignmentPolicy,
	"legacy post-alignment tree-stage policy");
my $postAlignmentQCPolicyMatches = $stateHasPolicies
	? buildTreeStatePolicyMatches($buildTreeState,
		'msa_selection_policy', $postAlignmentQCPolicy)
	: $legacyMsaSelectionPolicyMatches;
my $alignmentWorkPolicyMatches = $postAlignmentQCPolicyMatches;
my $postAlignmentPolicyMatches = $stateHasPolicies
	? buildTreeStatePolicyMatches($buildTreeState,
		'tree_stage_policy', $postAlignmentPolicy)
	: $legacyTreeStagePolicyMatches;
my $legacyWithinSpeciesQCAudit = !$stateHasPolicies
	&& !$taxonAwareLocusSelection && !$rateMergePartitions && $withinSpecies
	&& -s $postAlignmentQCReport && !-e $legacyPostAlignmentQCPolicyFile
	&& !-e $legacyAlignmentWorkPolicyFile;
my $postAlignmentQCAuditCurrent = $postAlignmentQCPolicyMatches
	&& (!$postAlignmentLocusQC || -s $postAlignmentQCReport);
$postAlignmentQCAuditCurrent = 1 if $legacyWithinSpeciesQCAudit;
my $durableCompletionTree = reusableCompletionTree($completionMarker, $outD);
my $requestedPrimaryMethods = $doIQTree + $doRAXML + $doRAXMLng
	+ $doFastTree + $doVeryFastTree;
my $completedTreeName = length($durableCompletionTree)
	? basename($durableCompletionTree) : '';
my $completionMatchesMethod = $requestedPrimaryMethods == 1 && (
	($doIQTree && $completedTreeName =~ /^IQtree.*\.treefile\z/)
	|| ($doRAXML && $completedTreeName =~ /^RXML.*\.nwk\z/)
	|| ($doRAXMLng && $completedTreeName =~ /^RXng.*\.nwk\z/)
	|| ($doFastTree && $completedTreeName =~ /^FASTTREE.*\.nwk\z/)
	|| ($doVeryFastTree && $completedTreeName =~ /^VERYFASTTREE.*\.nwk\z/)
);
my $hasAdditionalAnalysis = $Ete || $calcDistMat || $calcDNAdiff
	|| $doGenesToPh || $doSuperTree || $doSuperCheck || $doGubbins
	|| $doCFML || $useTreeShrink || $doDNDS || $doTheta
	|| $doFastGear || $doFastGearSummary || $gzipInput;
if (length($durableCompletionTree) && $completionMatchesMethod
		&& !$hasAdditionalAnalysis
		&& ($cogCats eq '' || ($alignmentWorkPolicyMatches
		&& $postAlignmentQCAuditCurrent && $postAlignmentPolicyMatches))) {
	# The marker is published only after tree validation and all requested standard
	# stages finish. A matching policy therefore avoids reopening every locus and
	# rescanning the concatenated alignment on a duplicate/resumed invocation.
	finalizeMSAArtifacts($MsaD, $MsaWorkD);
	safeRemoveTree($tmpD, $tmpBase);
	clearLifecycleMarker($terminalMarker, 'clear obsolete terminal no-tree marker');
	clearLifecycleMarker($placementPendingMarker, 'clear completed placement-pending marker');
	writeBuildTreeState();
	cleanupLegacyBuildTreeStateFiles();
	compactTaxonAwareDiagnostics();
	writeWorkflowHeartbeat('complete');
	print "Recovery state: durable completion marker and current policy match; "
		."skipping alignment/QC/tree revalidation ($durableCompletionTree)\n";
	exit(0);
}
my $doMSA = 1;
my $treesDone = treePresent($tOhr)
	&& (!$calcNonSyn || treePresent($tOhrNSun))
	&& (!$calcSyn || treePresent($tOhrSyn));
if ($strictBackbone && $treesDone
		&& (!-s "$treeD/strict_backbone.samples.tsv"
			|| !-s "$treeD/strict_backbone.epa_placements.tsv")) {
	if (-s $placementPendingMarker && fileGZe($multAliArtifact)) {
		print "Recovery state: validated backbone has pending EPA-ng placement; "
			."retaining inference and retrying placement only\n";
	} else {
		print "Recovery state: existing tree predates strict-backbone EPA-ng placement; "
			."rebuilding tree outputs from the retained alignment\n";
		safeRemoveTree($treeD, $outD);
		make_path($treeD);
		$treesDone = 0;
	}
}
if ($cogCats ne "" && $continue && !$alignmentWorkPolicyMatches) {
	print "Recovery state: MSA-selection policy changed; rebuilding per-locus alignments and tree outputs\n";
	safeRemoveTree($MsaD, $outD);
	safeRemoveTree($treeD, $outD);
	make_path($MsaD);
	make_path($treeD);
	$treesDone = 0;
} elsif ($cogCats ne "" && $continue && !$postAlignmentQCAuditCurrent) {
	print "Recovery state: post-alignment QC checkpoint is unavailable; rebuilding per-locus alignments and tree outputs\n";
	safeRemoveTree($MsaD, $outD);
	safeRemoveTree($treeD, $outD);
	make_path($MsaD);
	make_path($treeD);
	$treesDone = 0;
} elsif ($cogCats ne "" && $continue && !$postAlignmentPolicyMatches
		&& ($treesDone || fileGZe($multAliArtifact))) {
	my $postAlignmentQCBackup = "";
	if ($postAlignmentLocusQC && -s $postAlignmentQCReport) {
		my ($backupHandle, $backupPath) = tempfile(
			"post-alignment-locus-qc-XXXXXX", DIR => $tmpD, UNLINK => 1);
		retry_close($backupHandle, "close post-alignment QC backup");
		copy($postAlignmentQCReport, $backupPath)
			or die "Cannot preserve post-alignment QC report $postAlignmentQCReport: $!\n";
		$postAlignmentQCBackup = $backupPath;
	}
	print "Recovery state: downstream tree-stage policy changed; retaining the selected MSA and rebuilding tree outputs\n";
	clearLifecycleMarker($completionMarker, "clear completion before tree-stage rebuild");
	safeRemoveTree($treeD, $outD);
	make_path($treeD);
	if ($postAlignmentQCBackup ne "") {
		copy($postAlignmentQCBackup, $postAlignmentQCReport)
			or die "Cannot restore post-alignment QC report $postAlignmentQCReport: $!\n";
		retry_unlink($postAlignmentQCBackup, label => "remove post-alignment QC backup");
	}
	$treesDone = 0;
}
my ($primaryAlignmentReady, $primaryAlignmentReason, $primaryAlignmentMetadata) =
	treeAlignmentCheckpointStatus($multAliArtifact, $useAA4tree);
if ($continue && fileGZe($multAliArtifact) && !$primaryAlignmentReady) {
	warn "Recovery state: ignoring unusable retained alignment checkpoint "
		."$multAliArtifact ($primaryAlignmentReason); rebuilding it from the locus inputs\n";
}
my $siteAlignmentsReady = (!$calcSyn || fileGZe($multAliSynArtifact))
	&& (!$calcNonSyn || fileGZe($multAliNonSynArtifact));
my $reusableAlignment = $isAligned
	|| ($primaryAlignmentReady && $siteAlignmentsReady);
if ($continue) {
	if ($treesDone) {
		print "Recovery state: complete nonempty tree output found; retaining completed tree stages\n";
	} elsif ($reusableAlignment) {
		print "Recovery state: reusable nonempty alignment checkpoint found; rebuilding missing tree stages\n";
	} else {
		if ($cogCats ne '' && $alignmentWorkPolicyMatches) {
			print "Recovery state: retaining policy-matched completed per-locus alignments; "
				."rebuilding unfinished alignment and tree stages\n";
		} else {
			print "Recovery state: no reusable alignment or policy-matched locus checkpoint; "
				."restarting alignment and tree stages from input FASTA/category files\n";
			safeRemoveTree($MsaD, $outD);
		}
		safeRemoveTree($treeD, $outD);
		make_path($MsaD) unless -d $MsaD;
		make_path($treeD) unless -d $treeD;
	}
}
writeBuildTreeState();
cleanupLegacyBuildTreeStateFiles();
my $calcMSA = !$treesDone && !$primaryAlignmentReady;
my $retainedConcatenatedCheckpoint = !$calcMSA && $primaryAlignmentReady
	&& $cogCats ne '';
my $retainedCheckpointSampleCount = $retainedConcatenatedCheckpoint
	? ($primaryAlignmentMetadata->{sequences} // 0) : 0;
my $retainedPartitionLocusCount = $retainedConcatenatedCheckpoint
	? partitionLocusRangeCount($multAliArtifact.$partiExt) : 0;
$doMSA = !(
	$isAligned
	|| ($continue && ($treesDone || $primaryAlignmentReady) && $siteAlignmentsReady)
);
			
# Reopen retained compressed alignments in scratch without consuming the durable checkpoints.
restoreMSAArtifactSet($multAliArtifact, $multAli);
restoreMSAArtifactSet($multAliSynArtifact, $multAliSyn) if $calcSyn;
restoreMSAArtifactSet($multAliNonSynArtifact, $multAliNonSyn) if $calcNonSyn;
restoreMSAArtifactSet($placementAlignmentArtifact, $placementAlignment) if $strictBackbone;
my ($alignedLoci, $failedLoci, $candidateLoci) = (0, 0, 0);
if ($isAligned){
	my $alignedInput = $useAA4tree ? $aaFna : $fnFna;
	die "-isAligned currently requires an uncompressed input file: $alignedInput\n" unless -f $alignedInput;
	if (-e $alignedInput){
		unlink $multAli if -e $multAli || -l $multAli;
		my $source = abs_path($alignedInput) or die "Cannot resolve aligned input $alignedInput: $!\n";
		symlink($source, $multAli) || copy($source, $multAli)
			or die "Cannot link or copy aligned input $source to $multAli: $!\n";
		print "Pre-aligned input staged as $multAli\n";
	}
} elsif (!$doMSA && $cogCats ne ""){
	fillGeneList($cogCats);
} elsif ($doMSA && $cogCats ne ""){
	my $hr = readFasta($fnFna,1); my %FNA = %{$hr};
	my %FAA;
	if ($aaFna ne ""){
		$hr = readFasta($aaFna,1); %FAA = %{$hr};
	}
	print "Loaded sequence inputs: " . scalar(keys %FNA) . " nucleotide and "
		. scalar(keys %FAA) . " protein records\n";

	############# test length of fna sequences can be divided by 3 ##############################
		
	if($doLengthCheck){
		my @notThrees;
		my $FNAseq; my $length; my $div;
		while (($FNAseq) = each (%FNA)){
			$length = length($FNA{$FNAseq});
			$div = $length/3;
			push(@notThrees, $FNAseq) if($div =~ /\D/);
		}
		if (@notThrees) {
			my @examples = @notThrees > 10 ? @notThrees[0 .. 9] : @notThrees;
			warn "Nucleotide sequences not divisible by 3: " . scalar(@notThrees)
				. "; examples: " . join(", ", @examples)
				. (@notThrees > 10 ? " (+".(@notThrees - 10)." more)" : "") . "\n";
		}
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
	my ($qcExcludedSequences, $qcEmptiedLoci) = (0, 0);
	my ($geneLengthMinDroppedSequences, $geneLengthIncludeRetainedSequences,
		$geneLengthIncludeDroppedSequences, $geneLengthRecoveryCandidateSequences,
		$geneLengthRecoveredMSASequences) = (0, 0, 0, 0, 0);
	
	my ($xI,$ST)= gzipopen($cogCats,"CogCATs phylo");
	#open my $xI,"<$cogCats" or die "Can't open cogcats $cogCats\n";
	chomp(my @linesCats = <$xI>);
	close $xI;
	#first cleanup of cat file..
	my @linesCats2; my @linesCats3;
	my %inputSamples;
	my $inputSequences = 0;
	foreach (@linesCats){ #check first some parameters..
		$cnt++; my @spl = split /\t/;
		if (@spl && $spl[0] =~ m/^#/){shift @spl;}
		@spl = grep !/^NA$/, @spl;#remove NAs
		die "No sequence identifiers in category line ".($cnt + 1)."\n" unless @spl;
		if (%flaggedExcluded) {
			# Drop flagged samples before any length or prevalence statistic is
			# taken, so they cannot shift the per-locus Q90 length reference or the
			# locus occupancy that the retained samples are then judged against.
			my @retained = grep {
				my ($sp) = parseSeqId($_, "category line ".($cnt + 1));
				!$flaggedExcluded{$sp};
			} @spl;
			$qcExcludedSequences += scalar(@spl) - scalar(@retained);
			@spl = @retained;
			unless (@spl) {
				# Every observation of this locus came from an excluded sample.
				$qcEmptiedLoci++;
				push(@linesCats2, []);
				$genesPerCat[$cnt] = 0;
				next;
			}
		}
		my (@spl2, @splInclude);
		#$genesPerCat[$cnt] = scalar(@spl) ;
		my @geneLs;
		my ($sp, $gene) = parseSeqId($spl[0], "category line ".($cnt + 1));
		foreach my $seq (@spl){### $seq = genomeX_NOGY
			($sp) = parseSeqId($seq, "category line ".($cnt + 1));
			$inputSamples{$sp} = 1;
			$inputSequences++;
			$geneLengthSampleAudit{$sp}{$gene}{input} = 1;
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
			my $ntEquivalentLength = $charCnts{$sp}{$seq} * ($useAA4tree ? 3 : 1);
			my $passesGeneLengthMin = $charCnts{$sp}{$seq}
				>= ($qtl90NTcnt{$gene} * $ntFracGene);
			my $passesAbsoluteMinimum = !$taxonAwareLocusSelection
				|| $ntEquivalentLength >= $taxonAwareMinSequenceNT;
			my $passesQC = $passesGeneLengthMin && $passesAbsoluteMinimum;
			my $passesInclude = $charCnts{$sp}{$seq}
				>= ($qtl90NTcnt{$gene} * $ntFracGeneInclude);
			$geneLengthSampleAudit{$sp}{$gene}{gene_length_min} = 1
				if $passesGeneLengthMin;
			$geneLengthSampleAudit{$sp}{$gene}{absolute_minimum} = 1
				if $passesAbsoluteMinimum;
			$geneLengthSampleAudit{$sp}{$gene}{qc} = 1 if $passesQC;
			$geneLengthSampleAudit{$sp}{$gene}{include} = 1 if $passesInclude;
			if ($passesQC) {
				push(@spl2, $seq);
				$geneLengthQCSequence{$seq} = 1;
				$geneTooLong++;
			} else {
				$geneTooShort++;
			}
			$geneLengthMinDroppedSequences++ unless $passesGeneLengthMin;
			if ($passesInclude) {
				push @splInclude, $seq;
				$geneLengthIncludeRetainedSequences++;
				$geneLengthRecoveryCandidateSequences++ unless $passesQC;
			} else {
				$geneLengthIncludeDroppedSequences++;
			}
		}
		$geneLengthIncludeByGene{$gene} = \@splInclude;
		push(@linesCats2,\@spl2);
		#has to work with what is actually there, not what could have been..
		$genesPerCat[$cnt] = scalar(@spl2);
		#die;
	}
	#die;
	my $GenesQtl90 = quantile(0.9,@genesPerCat);
	$selectionAttrition{input_loci} = scalar(@linesCats);
	$selectionAttrition{input_sequences} = $inputSequences;
	$selectionAttrition{input_samples} = scalar(keys %inputSamples);
	$selectionAttrition{qc_excluded_samples} = scalar(keys %flaggedExcluded);
	$selectionAttrition{qc_excluded_sequences} = $qcExcludedSequences;
	$selectionAttrition{qc_emptied_loci} = $qcEmptiedLoci;
	if (%flaggedExcluded) {
		print "Flagged-sample exclusion: removed ".scalar(keys %flaggedExcluded)
			." sample(s) and $qcExcludedSequences sequence(s) marked unfit by $sampleQCFile"
			.($qcEmptiedLoci ? "; $qcEmptiedLoci locus/loci lost every observation" : "")
			."; per-sample verdicts and fractions are in $sampleQCFile\n";
	}
	$selectionAttrition{length_include_retained_sequences} =
		$geneLengthIncludeRetainedSequences;
	$selectionAttrition{length_include_filtered_sequences} =
		$geneLengthIncludeDroppedSequences;
	$selectionAttrition{length_recovery_candidate_sequences} =
		$geneLengthRecoveryCandidateSequences;
	my $initialGeneLengthAudit = writeGeneLengthSampleAudit(
		$geneLengthSampleReport, \%geneLengthSampleAudit, {}, {},
		$ntFracGene, $ntFracGeneInclude,
	);
	$selectionAttrition{gene_length_min_dropped_loci} =
		$initialGeneLengthAudit->{gene_length_min_dropped_loci};
	$selectionAttrition{gene_length_include_min_dropped_loci} =
		$initialGeneLengthAudit->{gene_length_include_min_dropped_loci};
	$selectionAttrition{gene_length_recovery_candidate_loci} =
		$initialGeneLengthAudit->{recovery_candidate_loci};
	my $GenesQtl50 = quantile(0.5,@genesPerCat);
	my $minimumCategorySequences = $GenesQtl90 * $fracMaxGenes90pct;
	$minimumCategorySequences = 1 if $minimumCategorySequences < 1;
	if ($taxonAwareLocusSelection) {
		my $candidateSelection = selectTaxonAwareCandidateLoci(
			categories => \@linesCats2,
			char_counts => \%charCnts,
			sequences => $useAA4tree ? \%FAA : \%FNA,
			candidate_limit => $taxonAwareMaxLoci + $taxonAwareCandidateExtra,
			final_limit => $taxonAwareMaxLoci,
			core_limit => $taxonAwareCoreLoci,
			target_loci => $taxonAwareTargetLoci,
			target_nt => $taxonAwareTargetNT,
			rescue_min_prevalence => $taxonAwareRescueMinPrevalence,
			rescue_prevalence_mode => $taxonAwareRescuePrevalenceMode,
			targets_from_gate => $taxonAwareTargetsFromGate,
			gate_gene_fraction => $GeneFracPSpec,
			gate_nt_fraction => $ntFrac,
			gate_minimum_nt => $ntCntTotal,
			presort_weight => $taxonAwarePresortWeight,
			preferred_core_genes => $preferredCoreGeneSet,
			use_aa => $useAA4tree,
			report => "$treeD/taxon_aware_locus_candidates.tsv",
		);
		if (defined($candidateSelection->{terminal_reason})
				&& length($candidateSelection->{terminal_reason})) {
			my $reason = $candidateSelection->{terminal_reason};
			$selectionAttrition{eligible_loci} = 0;
			$selectionAttrition{candidate_loci} = 0;
			$selectionAttrition{candidate_samples} = 0;
			$selectionAttrition{length_retained_sequences} = $geneTooLong;
			$selectionAttrition{length_filtered_sequences} = $geneTooShort;
			$selectionAttrition{aligned_loci} = 0;
			$selectionAttrition{alignment_failed_loci} = 0;
			$selectionAttrition{post_qc_loci} = 0;
			$selectionAttrition{final_loci} = 0;
			$selectionAttrition{final_samples} = 0;
			$selectionAttrition{backbone_samples} = 0;
			$selectionAttrition{placement_samples} = 0;
			$selectionAttrition{excluded_samples} = scalar(keys %inputSamples);
			clearLifecycleMarker($completionMarker, 'clear stale tree completion');
			clearLifecycleMarker($placementPendingMarker,
				'clear stale placement-pending marker');
			cleanupLegacyBuildTreeStateFiles();
			writeSelectionAttritionAudit($selectionAttritionReport, \%selectionAttrition);
			writeOutcomeMarker($terminalMarker, 'valid_no_tree', $reason, {
				input_loci => scalar(@linesCats), input_sequences => $inputSequences,
				input_samples => scalar(keys %inputSamples),
				length_retained_sequences => $geneTooLong,
				length_filtered_sequences => $geneTooShort,
			}, $outD);
			finalizeMSAArtifacts($MsaD, $MsaWorkD);
			safeRemoveTree($tmpD, $tmpBase);
			compactTaxonAwareDiagnostics();
			writeWorkflowHeartbeat('complete');
			print "BuildTree completed with a valid terminal no-tree outcome: $reason\n";
			exit(0);
		}
		@linesCats3 = @{$candidateSelection->{categories}};
		%taxonAwarePreMetrics = %{$candidateSelection->{metrics}};
		%taxonAwareUniverseSamples = %{$candidateSelection->{samples}};
		$selectionAttrition{eligible_loci} = scalar(keys %taxonAwarePreMetrics);
		$selectionAttrition{candidate_samples} = scalar(keys %taxonAwareUniverseSamples);
		print "Taxon-aware gene-category prefilter: retained "
			. scalar(@linesCats3) . "/" . scalar(@linesCats)
			. " alignment candidates (up to $taxonAwareCandidateExtra are QC backfill); "
			. "removed $geneTooShort of " . ($geneTooShort + $geneTooLong)
			. " sequence(s) below the $ntFracGene gene-length Q90/minimum-NT rule\n";
	} else {
		$cnt=-1;
		foreach my $aRef (@linesCats2){ #remove genes with just too few genes..
			$cnt++; my @spl = @{$aRef};
			if (@spl >= $minimumCategorySequences){ #$GenesQtl50 ||
				push(@linesCats3,\@spl);
			}
			#print @spl . " ";
		}
		print "Gene-category prefilter: retained " . scalar(@linesCats3) . "/"
			. scalar(@linesCats) . " categories; removed $geneTooShort of "
			. ($geneTooShort + $geneTooLong) . " sequence(s) below $ntFracGene of their gene-length Q90; "
			. "category-size Q90=$GenesQtl90, minimum category sequences=$minimumCategorySequences\n";
	}
	$selectionAttrition{candidate_loci} = scalar(@linesCats3);
	$selectionAttrition{length_retained_sequences} = $geneTooLong;
	$selectionAttrition{length_filtered_sequences} = $geneTooShort;
	warn "Only " . scalar(@linesCats3) . " gene categories remain after prefiltering "
		. "(Q50=$GenesQtl50, Q90=$GenesQtl90, category threshold="
		. "$minimumCategorySequences)\n"
		if @linesCats3 < 20;
	@linesCats2 = (); #make space..
	$cnt=-1;
	foreach my $aRef (@linesCats3){
		$cnt++; my @spl = @{$aRef};
		my ($firstSample, $gene) = parseSeqId($spl[0], "filtered category line ".($cnt + 1));
		foreach my $seq (@spl){
			my ($sp, $seqGene) = parseSeqId($seq, "filtered category line ".($cnt + 1));
			die "Wrong gene in $seq, expected $gene!\n" if ($seqGene ne $gene);
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
		#The MSA loop below feeds each retained category from the
		#GeneLengthIncludeMin set, not from the NTfiltPerGene-passing subset, so
		#the species prefilter has to score the same sequences the aligner will
		#actually use. Scoring only the QC-passing subset removed samples whose
		#sequences are mostly short before their own inclusion-based recovery
		#could apply, leaving GeneLengthIncludeMin without effect for exactly the
		#samples it exists to rescue.
		my $contributing = $geneLengthIncludeByGene{$gene} // \@spl;
		foreach my $seq (@{$contributing}){
			my ($sp) = parseSeqId($seq, "filtered category line ".($cnt + 1));
			#first check if gene gets removed
			#next if ( ($charCnts{$sp}{$seq} < $maxNtCnt{$gene} ) * $ntFracGene);
			#next if ( $charCnts{$sp}{$seq} < ($qtl90NTcnt{$gene} * $ntFracGene));
			$specList{$sp} ++;
			$totalNTs{$sp} += $charCnts{$sp}{$seq};
		}
	}
	if ($taxonAwareLocusSelection) {
		$specList{$_} //= 0 for keys %taxonAwareUniverseSamples;
		$totalNTs{$_} //= 0 for keys %taxonAwareUniverseSamples;
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

	my (%smplsRmvd, %samplesPassedHighThresholdQC);
	my $tooFewGenes=0;my $tooFewNTs=0;my $tooFewNTs2=0; my $specsRemain = 0;
	#print "Samples removed due to low gene presence:\n";
	my $OGfnd=0;
	if ($taxonAwareLocusSelection) {
		my ($sampleSelection, $selectionError);
		{
			local $@;
			$sampleSelection = eval {
				classifyTaxonAwareSamples(
					metrics => \%taxonAwarePreMetrics,
					samples => \%taxonAwareUniverseSamples,
					target_loci => $taxonAwareTargetLoci,
					target_nt => $taxonAwareTargetNT,
					minimum_anchor_nt => 1,
					selected_only => 1,
					outgroup => $outgroup,
				);
			};
			$selectionError = $@;
		}
		if (length($selectionError)) {
			completeTaxonAwareOutgroupAnchorTerminal('candidate_selection', $selectionError)
				if $selectionError =~ /^Taxon-aware selection could not retain an aligned anchor for outgroup /;
			die $selectionError;
		}
		writeTaxonAwareSampleAudit(
			"$treeD/taxon_aware_sample_candidates.tsv", $sampleSelection);
		for my $sp (@specs) {
			$OGfnd++ if $outgroup ne "" && $outgroup eq $sp;
			if (($sampleSelection->{$sp}{role} // "remove") eq "remove") {
				$smplsRmvd{$sp} = 1;
			} else {
				$samplesPassedHighThresholdQC{$sp} = 1;
				$specsRemain++;
			}
		}
		print "Taxon-aware species prefilter: retained $specsRemain/"
			. scalar(keys %specList) . "; only samples without any usable NT "
			. "in the selected candidates were removed before MSA; audit: "
			. "$treeD/taxon_aware_sample_candidates.tsv\n";
	} else {
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
			} else {
				$samplesPassedHighThresholdQC{$sp} = 1;
				$specsRemain ++;
			}
		}
	}
	$smplsRmvd{$_} = 1 for grep {
		!$samplesPassedHighThresholdQC{$_}
	} keys %inputSamples;
	
	if ($OGfnd == 0 && $outgroup ne ""){
		warn "Configured outgroup '$outgroup' is absent from the retained sequence set; "
			."continuing with an unrooted phylogeny\n";
		$outgroup = "";
	}
	#############################################################################################

	
	print "Species prefilter: retained $specsRemain/" . scalar(keys %specList)
		. "; removed for relative NT=$tooFewNTs, minimum NT=$tooFewNTs2, minimum genes=$tooFewGenes\n"
		unless $taxonAwareLocusSelection;
	print "Species input statistics: maximum genes=$maxGenes; gene-count Q90=$qtl90Genes; "
		. "maximum informative NT=$maxNtCntTotal; informative-NT Q90=$qtl90NTcntAll; "
		. "scored on the GeneLengthIncludeMin=$ntFracGeneInclude sequence set\n";
	#die "$maxGenes\n";
	@linesCats = (); #empty array


	#die;
	$cnt=-1; #line counter
	$alignedLoci = 0;
	$failedLoci = 0;
	$candidateLoci = scalar @linesCats3;
	foreach my $aRef (@linesCats3){#go over each gene category, building MSA for each
	#----------------- main MSA loop ----------------------
		$cnt++; my @spl = @{$aRef};
		if (@spl ==0){
			limitedWarn("empty gene category", "Ignoring empty filtered gene-category entry at index $cnt\n");
			next;
		}
		if ($spl[0] =~ m/^#/){shift @spl;}
		my @spl2 = parseSeqId($spl[0], "category line ".($cnt + 1));
		my $gene = $spl2[1];
		if (exists $geneLengthIncludeByGene{$gene}) {
			if ($taxonAwareLocusSelection) {
				my ($bestSequence) = bestGeneSequencesBySample(
					$geneLengthIncludeByGene{$gene}, \%charCnts,
					"recovery category ".($cnt + 1), $gene,
				);
				@spl = map { $bestSequence->{$_} } sort keys %{$bestSequence};
			} else {
				@spl = @{$geneLengthIncludeByGene{$gene}};
			}
		}
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
		my $finOutMSAaa = File::Spec->catfile($MsaWorkD, "$gene_file_stem.$cnt.faa");
		my $finOutMSA = File::Spec->catfile($MsaWorkD, "$gene_file_stem.$cnt.fna");
		my $publishedOutMSAaa = File::Spec->catfile($MsaD, "$gene_file_stem.$cnt.faa");
		my $publishedOutMSA = File::Spec->catfile($MsaD, "$gene_file_stem.$cnt.fna");
		my $endFileExists = $useAA4tree
			? fileGZs($publishedOutMSAaa)
			: fileGZs($publishedOutMSAaa) && fileGZs($publishedOutMSA);
		if ($endFileExists) {
			restoreCompressedMSAArtifact($publishedOutMSAaa, $finOutMSAaa);
			restoreCompressedMSAArtifact($publishedOutMSA, $finOutMSA)
				unless $useAA4tree;
		}
		
		open O,">$tmpInMSA" or die "Can;t open tmp faa file for MSA: $tmpInMSA\n";
		open O2,">$tmpInMSAnt" or die "Can;t open tmp fna file for MSA: $tmpInMSAnt\n";
		my $seqType = "AA";my $seqTypeOth = "NT";
		my $seqLength = 0; my $numSeq =0;
		my %currentRecoveredSampleGene;
		my $currentRecoveredSequences = 0;
		
		#1st: collate sequences
		#do here already per gene length check .. probably better for alignment
		foreach my $seq (@spl){### $seq = genomeX_NOGY
			my ($sp) = parseSeqId($seq, "category line ".($cnt + 1));
			next unless $samplesPassedHighThresholdQC{$sp};
			if (!$geneLengthSampleAudit{$sp}{$gene}{qc}) {
				$currentRecoveredSampleGene{$sp}{$gene} = 1;
				$currentRecoveredSequences++;
			}
			if (!$taxonAwareLocusSelection
					&& $specList{$sp} < ($qtl90Genes * $GeneFracPSpec)) {
				die "buildTree: GeneFracPSpec maxGenes shouldn't be here!\n";
			}
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
		for my $sample (keys %currentRecoveredSampleGene) {
			$geneLengthRecoveredForMSA{$sample}{$_} = 1
				for keys %{$currentRecoveredSampleGene{$sample}};
		}
		$geneLengthRecoveredMSASequences += $currentRecoveredSequences;

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
		my $msaCommandOK = 1;
		if (!$endFileExists) {
			$msaCommandOK = eval {
				systemW($cmd1."\n".$cmd2."\n");
				1;
			};
		}
		if (!$msaCommandOK || (!$endFileExists && !-s $tmpOutMSAaa)) {
			my $error = $@ || "MSA command completed without producing a nonempty output";
			$error =~ s/\s+$//;
			$failedLoci++;
			$excludedLoci{$gene} = 1;
			limitedWarn("failed locus alignment",
				"Warning: excluding locus $gene from future calculations: $error\n");
			unlink $_ for grep { defined($_) && -e $_ }
				($tmpInMSA, $tmpInMSAnt, $tmpOutMSAaa, $tmpOutMSA);
			next;
		}

		
		
		#move into its own for loop
		#dist mat related
		my $tmpDMatOth = "$MsaD/${seqTypeOth}_clustalo_percID_${cnt}_".int($seqLength).".txt";
		my $tmpDMat = "$MsaD/${seqType}_clustalo_percID_${cnt}_".int($seqLength).".txt";
		my $inFastaOth = $tmpInMSAnt;
		my $percIDhr; my $avgID; my $pIDsmplhr;
		if ($calcDistMat){ #for dmat: calc each gene spearately and merge scores later
			my $distanceOK = eval {
				($avgID,$pIDsmplhr,$percIDhr) =
					calcDisPos2($tmpInMSA,$tmpDMat,0,$ncore,$tmpD);
				if (-e $inFastaOth){
					($avgID,$pIDsmplhr,$percIDhr) =
						calcDisPos2($inFastaOth,$tmpDMatOth,1,$ncore,$tmpD);
					$calcDistMatExtGo = $avgID;
				} else {
					$calcDistMatExtGo = 0;
				}
				1;
			};
			if (!$distanceOK) {
				my $error = $@ || "unknown distance-matrix failure";
				$error =~ s/\s+$//;
				limitedWarn("failed optional locus distance matrix",
					"Warning: retaining locus $gene but omitting its distance matrix: $error\n");
				unlink $_ for grep { -e $_ } ($tmpDMat, $tmpDMatOth);
			}
		}
		
		
		
		#die "@MSAs\n";
		if (!$useAA4tree){
			#this part now is all concerned about NT level things..

			if (!$endFileExists){
				my $ntAlignmentOK = eval {
					convertMultAli2NT($tmpOutMSAaa,$tmpInMSAnt,$tmpOutMSA);
					die "AA-to-NT conversion completed without producing a nonempty output\n"
						unless -s $tmpOutMSA;
					# MSAfix remains immediately after creation and before publication.
					runMSAFix($tmpOutMSA, $maxGapPerCol);
					1;
				};
				if (!$ntAlignmentOK) {
					my $error = $@ || "unknown nucleotide-alignment failure";
					$error =~ s/\s+$//;
					$failedLoci++;
					$excludedLoci{$gene} = 1;
					limitedWarn("failed locus alignment",
						"Warning: excluding locus $gene from future calculations: $error\n");
					unlink $_ for grep { defined($_) && -e $_ }
						($tmpInMSA, $tmpInMSAnt, $tmpOutMSAaa, $tmpOutMSA,
							$finOutMSAaa, $finOutMSA,
							$tmpDMat, $tmpDMatOth);
					next;
				}
			}
		}
		unlink $tmpInMSA; unlink $tmpInMSAnt;
		if (!$endFileExists) {
			move($tmpOutMSAaa, $finOutMSAaa)
				or die "Cannot move $tmpOutMSAaa to scratch MSA $finOutMSAaa: $!\n";
			move($tmpOutMSA, $finOutMSA)
				or die "Cannot move $tmpOutMSA to scratch MSA $finOutMSA: $!\n"
				unless $useAA4tree;
			publishCompressedMSAArtifact($finOutMSAaa, $publishedOutMSAaa);
			publishCompressedMSAArtifact($finOutMSA, $publishedOutMSA)
				unless $useAA4tree;
		}

		if (!$useAA4tree) {
			my ($scratchMSAsyn, $scratchMSAnonsyn);
			my $siteSubsetOK = eval {
				($scratchMSAsyn,$scratchMSAnonsyn) =
					synPosOnly($finOutMSA,$finOutMSAaa,0,
						($ogrGenes ne "" ? $outgroup : ""),$calcSyn,$calcNonSyn);
				1;
			};
			if (!$siteSubsetOK) {
				my $error = $@ || "unknown synonymous-site derivation failure";
				$error =~ s/\s+$//;
				$failedLoci++;
				$excludedLoci{$gene} = 1;
				limitedWarn("failed locus site subsets",
					"Warning: excluding locus $gene from future calculations: $error\n");
				next;
			}
			push (@MSAs,$finOutMSA);
			$primaryAlignmentGene{$finOutMSA} = $gene;
			push (@MSAsSyn,$scratchMSAsyn)
				if defined($scratchMSAsyn) && $scratchMSAsyn ne "" && fileGZs($scratchMSAsyn);
			push (@MSAsNonSyn,$scratchMSAnonsyn)
				if defined($scratchMSAnonsyn) && $scratchMSAnonsyn ne "" && fileGZs($scratchMSAnonsyn);
		} else {
			push (@MSA_AA,$finOutMSAaa);
			$primaryAlignmentGene{$finOutMSAaa} = $gene;
		}
		push (@MSrm,$finOutMSAaa,$finOutMSA);
		$alignedLoci++;
		print "Prepared $alignedLoci/$candidateLoci locus alignments\n"
			if $alignedLoci == 1 || $alignedLoci % 25 == 0;
	}
	print "Per-locus alignment summary: $alignedLoci/$candidateLoci candidate loci prepared"
		. ($failedLoci ? "; $failedLoci failed and were excluded" : "") . "\n";
	$selectionAttrition{aligned_loci} = $alignedLoci;
	$selectionAttrition{alignment_failed_loci} = $failedLoci;
	$selectionAttrition{length_recovered_msa_sequences} =
		$geneLengthRecoveredMSASequences;
	my $finalGeneLengthAudit = writeGeneLengthSampleAudit(
		$geneLengthSampleReport, \%geneLengthSampleAudit, \%smplsRmvd,
		\%geneLengthRecoveredForMSA, $ntFracGene, $ntFracGeneInclude,
	);
	$selectionAttrition{gene_length_recovered_msa_loci} =
		$finalGeneLengthAudit->{recovered_for_msa_loci};
	print "Gene-length recovery summary: QC retained $geneTooLong sequence(s); "
		."$geneLengthRecoveryCandidateSequences passed inclusion-only; "
		."$geneLengthRecoveredMSASequences entered candidate MSA input; "
		."$geneLengthIncludeDroppedSequences failed GeneLengthIncludeMin; "
		."sample audit=$geneLengthSampleReport\n";
	
	my $mergPIDtag = "_merge";
	mergePids("$MsaD/",$cnt, "AA",$mergPIDtag) if ($calcDistMat); #merge different percIDs
	mergePids("$MsaD/",$cnt, "NT",$mergPIDtag) if ($calcDistMat); #merge different percIDs
	
	#could be used to filter genes further, but not for now
	#pogenStatsFilter();
	print "Outgroup coverage: $ogrpCnt/$candidateLoci candidate loci contained '$outgroup'\n"
		if $outgroup ne "";
} elsif ($doMSA) {#no marker way, single gene
	my $r1; my $r2;
	#,$r1,$r2)
	$multAli = singleGeneMSAprocess($multAli)#;,\@MSAs,\@MSA_AA);
	#@MSAs = @{$r1};	@MSA_AA = @{$r2};
}
if ($synSummaryCount) {
	print "Synonymous-site classification summary: $synSummaryCount alignment(s), "
		. "$synSiteTotal synonymous-variable and $nonSynSiteTotal nonsynonymous-variable codon(s)\n";
}

my $msaOnlyCompletionMarker = File::Spec->catfile($outD, 'msaOnly.complete.tsv');
if ($onlyMSA) {
	if ($cogCats eq '') {
		my $extension = $useAA4tree ? 'faa' : 'fna';
		my $singleAlignment = File::Spec->catfile($MsaD, "single_locus.$extension");
		publishCompressedMSAArtifact($multAli, $singleAlignment);
	}
	my $artifacts = msaOnlyArtifacts($MsaD);
	die "MSA-only mode did not produce a non-empty per-locus alignment artifact\n"
		unless @{$artifacts};
	clearLifecycleMarker($completionMarker, 'clear tree completion in MSA-only mode');
	clearLifecycleMarker($terminalMarker, 'clear obsolete terminal no-tree marker');
	clearLifecycleMarker($placementPendingMarker,
		'clear obsolete placement-pending marker');
	writeOutcomeMarker($msaOnlyCompletionMarker, 'msa_complete',
		'localized per-locus alignments completed; combined-MSA postprocessing, concatenation, and phylogeny intentionally skipped',
		{ alignment_directory => $MsaD, artifacts => scalar(@{$artifacts}) }, $outD);
	writeBuildTreeState();
	cleanupLegacyBuildTreeStateFiles();
	writeWorkflowHeartbeat('complete');
	safeRemoveTree($tmpD, $tmpBase);
	print "BuildTree MSA-only mode completed successfully; per-locus alignments="
		.scalar(@{$artifacts})."; directory=$MsaD\n";
	exit(0);
}
clearLifecycleMarker($msaOnlyCompletionMarker,
	'clear obsolete MSA-only completion marker');

print "\n---------------- POST-ALIGNMENT WORKFLOW ----------------\n";
my $postAlignmentStepStarted = time;

if ($postAlignmentLocusQC && $cogCats ne "") {
	my $primaryAlignments = $useAA4tree ? \@MSA_AA : \@MSAs;
	if (@{$primaryAlignments}) {
		my $kept = runPostAlignmentLocusQC(
			$primaryAlignments,
			$useAA4tree ? 'aa' : 'nt',
			$postAlignmentQCReport,
		);
		my %keepPath = map { $_ => 1 } @{$kept};
		my %keepStem = map { alignmentFileStem($_) => 1 } @{$kept};
		if ($useAA4tree) {
			@MSA_AA = grep { $keepPath{$_} } @MSA_AA;
		} else {
			@MSAs = grep { $keepPath{$_} } @MSAs;
			@MSAsSyn = grep { $keepStem{alignmentFileStem($_)} } @MSAsSyn;
			@MSAsNonSyn = grep { $keepStem{alignmentFileStem($_)} } @MSAsNonSyn;
		}
	}
} elsif ($cogCats ne "") {
	my $primaryAlignments = $useAA4tree ? \@MSA_AA : \@MSAs;
	my $candidateCount = scalar @{$primaryAlignments};
	print "Post-alignment locus QC disabled; retaining all $candidateCount prepared loci\n"
		if $candidateCount;
	if (-e $postAlignmentQCReport) {
		unlink $postAlignmentQCReport
			or die "Cannot remove stale locus-QC report $postAlignmentQCReport: $!\n";
	}
}
my $postQCPrimary = $useAA4tree ? \@MSA_AA : \@MSAs;
my $postAlignmentStats = $postAlignmentLocusQC && $cogCats ne ""
	&& -s $postAlignmentQCReport
	? alignmentCollectionStatsFromReport($postAlignmentQCReport)
	: alignmentCollectionStats($postQCPrimary);
if ($retainedConcatenatedCheckpoint && !$postAlignmentStats->{loci}) {
	$postAlignmentStats->{loci} = $retainedPartitionLocusCount;
}
my $reportedPostQCLoci = $retainedConcatenatedCheckpoint
	? $postAlignmentStats->{loci} : scalar(@{$postQCPrimary});
$selectionAttrition{post_qc_loci} = $reportedPostQCLoci;
postAlignmentStep("locus QC", $postAlignmentStepStarted,
	"enabled=".($postAlignmentLocusQC ? 1 : 0),
	"retained_loci=$reportedPostQCLoci",
	"source=".($retainedConcatenatedCheckpoint ? 'retained_checkpoint' : 'current_run'),
	"report=$postAlignmentQCReport");
$postAlignmentStepStarted = time;

postAlignmentStep("alignment inventory", $postAlignmentStepStarted,
	"loci=$postAlignmentStats->{loci}",
	"mean_sequences_per_locus=$postAlignmentStats->{mean_sequences}",
	"mean_alignment_length=$postAlignmentStats->{mean_length}",
	"total_alignment_sites=$postAlignmentStats->{total_sites}",
	"sequence_range=$postAlignmentStats->{minimum_sequences}-$postAlignmentStats->{maximum_sequences}",
	"length_range=$postAlignmentStats->{minimum_length}-$postAlignmentStats->{maximum_length}");
$postAlignmentStepStarted = time;

if ($taxonAwareLocusSelection && $cogCats ne "") {
	my $primaryAlignments = $useAA4tree ? \@MSA_AA : \@MSAs;
	if (@{$primaryAlignments}) {
		my $postQCAlignmentCount = scalar @{$primaryAlignments};
		my ($finalSelection, $selectionError);
		{
			local $@;
			$finalSelection = eval {
				selectTaxonAwareFinalLoci(
					alignments => $primaryAlignments,
					path_gene => \%primaryAlignmentGene,
					pre_metrics => \%taxonAwarePreMetrics,
					samples => \%taxonAwareUniverseSamples,
					maximum_loci => $taxonAwareMaxLoci,
					core_loci => $taxonAwareCoreLoci,
					target_loci => $taxonAwareTargetLoci,
					target_nt => $taxonAwareTargetNT,
					rescue_min_prevalence => $taxonAwareRescueMinPrevalence,
					rescue_prevalence_mode => $taxonAwareRescuePrevalenceMode,
					targets_from_gate => $taxonAwareTargetsFromGate,
					gate_gene_fraction => $GeneFracPSpec,
					gate_nt_fraction => $ntFrac,
					gate_minimum_nt => $ntCntTotal,
					presort_weight => $taxonAwarePresortWeight,
					information_saturation => $taxonAwareInformationSaturation,
					excess_variation_onset => $taxonAwareExcessVariationOnset,
					excess_variation_span => $taxonAwareExcessVariationSpan,
					preferred_core_genes => $preferredCoreGeneSet,
					minimum_anchor_nt => 1,
					use_aa => $useAA4tree,
					qualification_sequences => \%geneLengthQCSequence,
					outgroup => $outgroup,
					locus_report => "$treeD/taxon_aware_locus_selection.tsv",
					sample_report => "$treeD/taxon_aware_sample_selection.tsv",
				);
			};
			$selectionError = $@;
		}
		if (length($selectionError)) {
			completeTaxonAwareOutgroupAnchorTerminal('final_selection', $selectionError)
				if $selectionError =~ /^Taxon-aware selection could not retain an aligned anchor for outgroup /;
			die $selectionError;
		}
		my %keepPath = map { $_ => 1 } @{$finalSelection->{alignments}};
		my %keepStem = map { alignmentFileStem($_) => 1 }
			@{$finalSelection->{alignments}};
		if ($useAA4tree) {
			@MSA_AA = grep { $keepPath{$_} } @MSA_AA;
		} else {
			@MSAs = grep { $keepPath{$_} } @MSAs;
			@MSAsSyn = grep { $keepStem{alignmentFileStem($_)} } @MSAsSyn;
			@MSAsNonSyn = grep { $keepStem{alignmentFileStem($_)} } @MSAsNonSyn;
		}
		for my $sample (keys %samples) {
			delete $samples{$sample}
				if ($finalSelection->{sample_metrics}{$sample}{role} // "remove") eq "remove";
		}
		%taxonAwareFinalMetricByPath = map {
			my $metric = $finalSelection->{locus_metrics}{$_};
			$metric->{path} => $metric
		} grep {
			$finalSelection->{locus_metrics}{$_}{selected}
		} keys %{$finalSelection->{locus_metrics}};
		# Backbone coverage is the primary sample filter. With -placeOnBackbone 1
		# a sample failing it is deferred to placement; with placement off there is
		# nothing behind it, so it is removed here instead of entering the tree on
		# the strength of a single informative nucleotide.
		my $backboneEligibility = classifyTaxonAwareCoverageEligibility(
			sample_metrics => $finalSelection->{sample_metrics},
			gene_fraction => $GeneFracPSpec,
			nt_fraction => $ntFrac,
			minimum_nt => $ntCntTotal,
			minimum_loci_floor => 1,
			role => 'backbone', outgroup => $outgroup,
		);
		%taxonAwareBackboneEligibility = map {
			$_ => $backboneEligibility->{samples}{$_}{eligible}
		} keys %{$backboneEligibility->{samples}};
		%taxonAwareBackboneIneligibleReason = map {
			$_ => $backboneEligibility->{samples}{$_}{reason}
		} grep { !$backboneEligibility->{samples}{$_}{eligible} }
			keys %{$backboneEligibility->{samples}};
		my $backboneAudit = "$treeD/taxon_aware_backbone_eligibility.tsv";
		open my $backboneOutput, '>', $backboneAudit
			or die "Cannot write taxon-aware backbone eligibility $backboneAudit: $!\n";
		print {$backboneOutput} "sample\tselected_loci\tselected_nt\tbackbone_eligible\treason\tminimum_loci\tminimum_nt\n";
		for my $sample (sort keys %{$backboneEligibility->{samples}}) {
			my $entry = $backboneEligibility->{samples}{$sample};
			print {$backboneOutput} join("\t", $sample, $entry->{selected_loci},
				$entry->{selected_nt}, $entry->{eligible} ? 1 : 0, $entry->{reason},
				$backboneEligibility->{minimum_loci}, $backboneEligibility->{minimum_nt}), "\n";
		}
		close $backboneOutput
			or die "Cannot close taxon-aware backbone eligibility $backboneAudit: $!\n";
		if (!$strictBackbone && $enforceSampleCoverage) {
			my @coverageRemoved = sort grep {
				exists($samples{$_}) && !$taxonAwareBackboneEligibility{$_}
			} keys %taxonAwareBackboneEligibility;
			delete $samples{$_} for @coverageRemoved;
			$selectionAttrition{coverage_excluded_samples} = scalar(@coverageRemoved);
			my %removedReason;
			$removedReason{$taxonAwareBackboneIneligibleReason{$_} // 'backbone_coverage_not_met'}++
				for @coverageRemoved;
			print "Sample coverage filter: removed ".scalar(@coverageRemoved)." of "
				.scalar(keys %taxonAwareBackboneEligibility)." sample(s) below "
				."$backboneEligibility->{minimum_loci} selected locus/loci or "
				."$backboneEligibility->{minimum_nt} informative NT "
				."(-GenesPerSpecies=$GeneFracPSpec, -relativeNTFraction=$ntFrac, "
				."-NTfiltCount=$ntCntTotal)"
				.(%removedReason ? "; reasons: ".join(", ",
					map { "$_=$removedReason{$_}" } sort keys %removedReason) : "")
				."; audit=$backboneAudit\n";
		} else {
			$selectionAttrition{coverage_excluded_samples} = 0;
			print "Sample coverage filter: "
				.($strictBackbone
					? "samples below the gate are deferred to EPA-ng placement"
					: "disabled by -enforceSampleCoverage 0; every aligned sample is retained")
				."; audit=$backboneAudit\n";
		}
		if ($strictBackbone) {
		my $placementEligibility = classifyTaxonAwareCoverageEligibility(
			sample_metrics => $finalSelection->{sample_metrics},
			gene_fraction => $placementGeneFracPSpec,
			nt_fraction => $placementNTFrac,
			minimum_nt => $placementNTCntTotal,
			minimum_overlap => $placementMinOverlap,
			minimum_loci_floor => 2,
			role => 'placement',
			outgroup => $outgroup,
		);
		$strictPlacementMinimumNT = $placementEligibility->{minimum_nt};
		%taxonAwarePlacementEligibility = map {
			$_ => $placementEligibility->{samples}{$_}{eligible}
		} keys %{$placementEligibility->{samples}};
		%taxonAwarePlacementIneligibleReason = map {
			$_ => $placementEligibility->{samples}{$_}{reason}
		} grep { !$placementEligibility->{samples}{$_}{eligible} }
			keys %{$placementEligibility->{samples}};
		my $placementAudit = "$treeD/taxon_aware_placement_eligibility.tsv";
		open my $placementOutput, '>', $placementAudit
			or die "Cannot write taxon-aware placement eligibility $placementAudit: $!\n";
		print {$placementOutput} join("\t", qw(
			sample selected_loci selected_nt placement_eligible reason
			minimum_loci minimum_nt
		)), "\n";
		for my $sample (sort keys %{$placementEligibility->{samples}}) {
			my $entry = $placementEligibility->{samples}{$sample};
			print {$placementOutput} join("\t",
				$sample, $entry->{selected_loci}, $entry->{selected_nt},
				$entry->{eligible} ? 1 : 0, $entry->{reason},
				$placementEligibility->{minimum_loci},
				$placementEligibility->{minimum_nt},
			), "\n";
		}
		close $placementOutput
			or die "Cannot close taxon-aware placement eligibility $placementAudit: $!\n";
		print "Taxon-aware final selection retained "
			. scalar(@{$finalSelection->{alignments}}) . "/"
			. $postQCAlignmentCount
			. " post-QC loci; reports: $treeD/taxon_aware_locus_selection.tsv, "
			. "$treeD/taxon_aware_sample_selection.tsv, $backboneAudit, $placementAudit\n";
		} else {
			# Backbone eligibility is also the generic sample-coverage audit when
			# placement is disabled. Keep it for compactTaxonAwareDiagnostics() to
			# publish instead of silently deleting the evidence behind removals.
			retry_unlink("$treeD/taxon_aware_placement_eligibility.tsv");
			print "Taxon-aware final selection retained "
				. scalar(@{$finalSelection->{alignments}}) . "/"
				. $postQCAlignmentCount
				. " post-QC loci; placement disabled, reports: "
				. "$treeD/taxon_aware_locus_selection.tsv, "
				. "$treeD/taxon_aware_sample_selection.tsv, $backboneAudit\n";
		}
	}
}
my $postSelectionPrimary = $useAA4tree ? \@MSA_AA : \@MSAs;
my $reportedSelectedLoci = scalar(@{$postSelectionPrimary});
my $reportedSelectedSamples = scalar(keys %samples);
if ($retainedConcatenatedCheckpoint) {
	$reportedSelectedLoci = $retainedPartitionLocusCount;
	$reportedSelectedLoci = $reportedPostQCLoci unless $reportedSelectedLoci;
	$reportedSelectedSamples = $retainedCheckpointSampleCount;
}
$selectionAttrition{final_loci} = $reportedSelectedLoci;
$selectionAttrition{final_samples} = $reportedSelectedSamples;
postAlignmentStep("taxon-aware locus selection", $postAlignmentStepStarted,
	"enabled=".($taxonAwareLocusSelection ? 1 : 0),
	"selected_loci=$reportedSelectedLoci",
	"samples=$reportedSelectedSamples",
	"source=".($retainedConcatenatedCheckpoint ? 'retained_checkpoint' : 'current_run'));
$postAlignmentStepStarted = time;

if ($rateMergePartitions && $cogCats ne "") {
	%partitionRateProxy = %{readPostAlignmentRateMetrics($postAlignmentQCReport)};
	for my $alignment (keys %taxonAwareFinalMetricByPath) {
		$partitionSelectionPhase{$alignment} =
			$taxonAwareFinalMetricByPath{$alignment}{selection_phase} // '';
	}
}
postAlignmentStep("rate/GC partition preparation", $postAlignmentStepStarted,
	"enabled=".($rateMergePartitions ? 1 : 0),
	"locus_metrics=".scalar(keys %partitionRateProxy));
$postAlignmentStepStarted = time;

#die "@MSA_AA\n\n";
if ($calcMSA
		&& (($cogCats ne '' && @MSAs == 0 && @MSA_AA == 0) || $cogCats eq '')) {
	my $reason = $cogCats ne '' ? 'no_usable_loci' : 'single_gene_alignment_failed';
	clearLifecycleMarker($completionMarker, 'clear stale tree completion');
	clearLifecycleMarker($placementPendingMarker, 'clear stale placement-pending marker');
	cleanupLegacyBuildTreeStateFiles();
	$selectionAttrition{backbone_samples} = 0;
	$selectionAttrition{placement_samples} = 0;
	$selectionAttrition{excluded_samples} = $selectionAttrition{final_samples} eq 'NA' ? 0 : $selectionAttrition{final_samples};
	writeSelectionAttritionAudit($selectionAttritionReport, \%selectionAttrition);
	writeOutcomeMarker($terminalMarker, 'valid_no_tree', $reason, {
		candidate_loci => $candidateLoci, aligned_loci => $alignedLoci,
		failed_loci => $failedLoci, samples => scalar(keys %samples),
	}, $outD);
			finalizeMSAArtifacts($MsaD, $MsaWorkD);
			safeRemoveTree($tmpD, $tmpBase);
			compactTaxonAwareDiagnostics();
			writeWorkflowHeartbeat('complete');
	print "BuildTree completed with a valid terminal no-tree outcome: $reason\n";
	exit(0);
}

#prep final MSA file that is correct NT or AA and is merged. A validated
#concatenated checkpoint is already the selected/merged result, so do not call
#mergeMSAs with intentionally empty per-locus arrays during a tree-only resume.
if ($retainedConcatenatedCheckpoint) {
	# Retain the existing alignment and partition pair without rewriting either.
} elsif (!$useAA4tree) {
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
unless ($retainedConcatenatedCheckpoint) {
	publishMSAArtifactSet($multAli, $multAliArtifact);
	publishMSAArtifactSet($multAliSyn, $multAliSynArtifact) if $calcSyn;
	publishMSAArtifactSet($multAliNonSyn, $multAliNonSynArtifact) if $calcNonSyn;
}

postAlignmentStep("concatenation", $postAlignmentStepStarted,
	"loci=$reportedSelectedLoci", "samples=$reportedSelectedSamples",
	"source=".($retainedConcatenatedCheckpoint ? 'retained_checkpoint' : 'current_run'),
	"alignment=$multAli");
$postAlignmentStepStarted = time;

if ($strictBackbone) {
	my $fullAlignmentArtifact = File::Spec->catfile($MsaD, 'MSAli.full.fna');
	my $fullAlignment = File::Spec->catfile($MsaWorkD, 'MSAli.full.fna');
	restoreMSAArtifactSet($fullAlignmentArtifact, $fullAlignment)
		unless $calcMSA;
	if (!$calcMSA && -s $fullAlignment) {
		my ($fullAlignmentReady, $fullAlignmentReason) =
			treeAlignmentCheckpointStatus($fullAlignment, $useAA4tree);
		if (!$fullAlignmentReady) {
			warn "Strict-backbone: replacing unusable retained full alignment "
				."$fullAlignment ($fullAlignmentReason) from $multAli\n";
			copy($multAli, $fullAlignment)
				or die "Cannot replace unusable full alignment $fullAlignment: $!\n";
		}
	}
	if ($calcMSA || !-s $fullAlignment) {
		copy($multAli, $fullAlignment)
			or die "Cannot preserve full alignment as $fullAlignment: $!\n";
	}
	# Already read once during preflight; flagged samples were removed from
	# the alignment before this point, so what remains here are coverage decisions.
	my $sampleStatus = \%sampleQCStatus;
	$strictSplit = split_strict_backbone(
		$fullAlignment, $multAli, $placementAlignment, $sampleStatus,
		{
			is_aa => $useAA4tree,
			coverage_fraction => $strictBackboneFraction,
			minimum_backbone => $strictBackboneMinSamples,
			backbone_eligible => \%taxonAwareBackboneEligibility,
			backbone_ineligible_reason => \%taxonAwareBackboneIneligibleReason,
			placement_eligible => \%taxonAwarePlacementEligibility,
			placement_ineligible_reason => \%taxonAwarePlacementIneligibleReason,
			outgroup => $outgroup,
			partition_file => $multAli.$partiExt,
			minimum_backbone_overlap_nt => $strictPlacementMinimumNT,
			minimum_backbone_overlap_loci => $strictPlacementMinimumLoci,
		},
	);
	my ($strictAlignmentReady, $strictAlignmentReason) =
		treeAlignmentCheckpointStatus($multAli, $useAA4tree);
	die "Strict-backbone produced an unusable inference alignment $multAli: "
		."$strictAlignmentReason\n" unless $strictAlignmentReady;
	publishMSAArtifactSet($fullAlignment, $fullAlignmentArtifact);
	publishMSAArtifactSet($multAli, $multAliArtifact);
	publishMSAArtifactSet($placementAlignment, $placementAlignmentArtifact)
		if -s $placementAlignment;

	my $classificationFile = "$treeD/strict_backbone.samples.tsv";
	open my $classification, '>', $classificationFile
		or die "Cannot write $classificationFile: $!\n";
	print {$classification} join("\t",
		qw(sample tree_role reason informative_positions q90_informative
			backbone_overlap_nt backbone_overlap_loci backbone_state_divergence)), "\n";
	my %isPlacement = map { $_ => 1 } @{$strictSplit->{placement}};
	my %isExcluded = map { $_ => 1 } @{$strictSplit->{excluded} // []};
	for my $sample (sort(
		@{$strictSplit->{backbone}}, @{$strictSplit->{placement}}, @{$strictSplit->{excluded} // []}
	)) {
		my $reason = $strictSplit->{reason}{$sample} // 'validated_backbone';
		$reason = "backbone_fallback:".$strictSplit->{requested_reason}{$sample}
			if $strictSplit->{fallback}
				&& exists($strictSplit->{requested_reason}{$sample});
		my $overlapMetric = $strictSplit->{backbone_overlap}{$sample} || {};
		print {$classification} join("\t",
			$sample,
			$isExcluded{$sample} ? 'excluded' : $isPlacement{$sample} ? 'placement' : 'backbone',
			$reason,
			$strictSplit->{informative}{$sample},
			sprintf('%.2f', $strictSplit->{q90_informative}),
			defined($overlapMetric->{backbone_overlap_nt})
				? $overlapMetric->{backbone_overlap_nt} : 'NA',
			defined($overlapMetric->{backbone_overlap_loci})
				? $overlapMetric->{backbone_overlap_loci} : 'NA',
			defined($overlapMetric->{backbone_state_divergence})
				? sprintf('%.8g', $overlapMetric->{backbone_state_divergence}) : 'NA',
		), "\n";
	}
	close $classification or die "Cannot close $classificationFile: $!\n";
	print "Strict-backbone split: ".scalar(@{$strictSplit->{backbone}})
		." backbone and ".scalar(@{$strictSplit->{placement}})
		." placement sample(s), ".scalar(@{$strictSplit->{excluded} // []})
		." excluded from placement; required backbone overlap="
		."$strictPlacementMinimumNT NT across $strictPlacementMinimumLoci loci; "
		."full alignment retained at $fullAlignment\n";
	warn "Strict-backbone fallback: fewer than $strictBackboneMinSamples validated "
		."backbone samples remained, so all samples were used for inference; "
		."see $classificationFile\n"
		if $strictSplit->{fallback};
}
$selectionAttrition{backbone_samples} = $strictSplit
	? scalar(@{$strictSplit->{backbone}}) : $reportedSelectedSamples;
$selectionAttrition{placement_samples} = $strictSplit
	? scalar(@{$strictSplit->{placement}}) : 0;
$selectionAttrition{excluded_samples} = $strictSplit
	? scalar(@{$strictSplit->{excluded} || []}) : 0;
writeSelectionAttritionAudit($selectionAttritionReport, \%selectionAttrition);

postAlignmentStep("strict-backbone preparation", $postAlignmentStepStarted,
	"enabled=".($strictBackbone ? 1 : 0),
	"backbone_samples=".($strictSplit ? scalar(@{$strictSplit->{backbone}}) : $reportedSelectedSamples),
	"placement_samples=".($strictSplit ? scalar(@{$strictSplit->{placement}}) : 0));
$postAlignmentStepStarted = time;

#phylip conversion??
if ( $doGenesToPh){ 
	my $phylipD = File::Spec->catdir($outD, "phylip");
	make_path($phylipD) unless -d $phylipD;
	my $fasta2phylip = getProgPaths("fasta2phylip_scr");

	foreach my $MSAfn (@MSAs){
		my $phylipOut = File::Spec->catfile($phylipD, basename($MSAfn).".ph");
		retry_unlink($_, fatal => 0, label => "clean optional PHYLIP output")
			for glob("$phylipOut*");
		my $cmd2 = "$fasta2phylip -c 50 ".shellQuote($MSAfn)." > ".shellQuote($phylipOut)."\n";
		my $phylipOK = eval {
			systemW $cmd2;
			die "conversion did not produce a nonempty output\n" unless -s $phylipOut;
			1;
		};
		if (!$phylipOK) {
			my $error = $@ || "unknown PHYLIP conversion failure";
			$error =~ s/\s+$//;
			limitedWarn("failed optional per-locus PHYLIP conversion",
				"Warning: omitting PHYLIP output for $MSAfn: $error\n");
			unlink $phylipOut if -e $phylipOut;
			next;
		}
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
		print "Building subtree " . ($i + 1) . "/" . scalar(@theRealMSAs) . "\n"
			if $i == 0 || ($i + 1) % 25 == 0;
		my $tOhrST = createTreeOpt($theRealMSAs[$i],"allsites",$i,1,"");
		my $subtreeFile;
		my $subtreeOK = eval {
			my $trRetH = treeAtHeart($tOhrST);
			$subtreeFile = ${$trRetH}{nwk} // "";
			die "subtree method produced no nonempty tree\n" unless -s $subtreeFile;
			1;
		};
		if (!$subtreeOK) {
			my $error = $@ || "unknown subtree failure";
			$error =~ s/\s+$//;
			limitedWarn("failed locus subtree",
				"Warning: excluding subtree for $theRealMSAs[$i]: $error\n");
			next;
		}
		push(@treeCol,$subtreeFile);
		if ($calcSyn){
		} 
		if ($calcNonSyn){
		}
	}
	print "Subtree construction summary: " . scalar(@treeCol) . " subtree(s) prepared\n";
	if ($doSuperTree){
		die "No usable locus subtrees remain; no supertree can be created\n" unless @treeCol;
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
	$gubbinsBin = requireConfiguredTool("gubbins", "Gubbins") if $gubbinsBin eq "";
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
my $inferenceStarted = time;
if ($doSuperTree){
	$Tree1{nwk} = $phyloTree;
	$trRetH = \%Tree1;
} else {
	$trRetH = treeAtHeart($tOhr);
	# Constrain the site-subset topologies to the completed all-sites tree.
	my $allSitesTree = ${$trRetH}{nwk} // "";
	$allSitesTree = "" unless $allSitesTree ne "" && -s $allSitesTree;
	if ($calcSyn){ #tree at syn pos
		$tOhrSyn->{constraintTree} = $allSitesTree;
		treeAtHeart($tOhrSyn);
	}
	if ($calcNonSyn){ #tree at non-syn pos
		$tOhrNSun->{constraintTree} = $allSitesTree;
		treeAtHeart($tOhrNSun);
	}
}
postAlignmentStep("phylogeny inference", $inferenceStarted,
	"methods=".(@treeMethods ? join('+', @treeMethods) : 'none'),
	"model=".($treeAutoModel ? 'AutoModel' : $useAA4tree ? 'LG+F+G' : 'GTR+F+G2'),
	"primary_tree=".(${$trRetH}{nwk} // '<none>'));

my $placementStarted = time;
if ($strictSplit) {
	my $backboneTree = ${$trRetH}{nwk} // "";
	if ($backboneTree ne "" && -s $backboneTree) {
		my $primaryTree = $backboneTree;
		my $dedicatedBackbone = $primaryTree =~ s/\.backbone\.treefile$/.treefile/;
		my $report = "$treeD/strict_backbone.epa_placements.tsv";
		my @placementReportColumns = qw(
			sample status backbone_overlap_nt backbone_overlap_loci
			backbone_state_divergence edge likelihood likelihood_weight_ratio
			edpl candidate_placements distal_length backbone_distal_length pendant_length
			pendant_outlier_limit placement_filter_reason reason
		);
		if (@{$strictSplit->{placement}}) {
			my ($epaResult, $modelArtifact, $jplaceFile);
			my $retainedJplace = File::Spec->catfile(
				$treeD, 'epa-ng', 'epa_result.jplace');
			my $placementOK = eval {
				if ($continue && $dedicatedBackbone && !-s $primaryTree
						&& -s $retainedJplace) {
					$jplaceFile = $retainedJplace;
					$modelArtifact = 'retained EPA-ng result';
					$epaResult = read_epa_jplace(
						$jplaceFile, $strictSplit->{placement});
					print "Recovery state: final EPA-placed tree is missing; "
						."reusing $jplaceFile and reapplying placement filtering\n";
				} else {
					($epaResult, $modelArtifact, $jplaceFile) = runEpaNgPlacement(
						$tOhr, $backboneTree, $multAli, $placementAlignment,
						$strictSplit->{placement}, $treeD,
					);
				}
				1;
			};
			if (!$placementOK) {
				my $error = $@ || 'unknown EPA-ng placement failure';
				$error =~ s/\s+\z//;
				clearLifecycleMarker($completionMarker, 'clear stale tree completion');
				writeOutcomeMarker($placementPendingMarker, 'placement_pending', $error, {
					backbone_tree => $backboneTree,
					query_samples => scalar(@{$strictSplit->{placement}}),
				}, $outD);
				warn "EPA-ng placement deferred; the validated backbone and compressed MSA were retained: $error\n";
				finalizeMSAArtifacts($MsaD, $MsaWorkD);
				safeRemoveTree($tmpD, $tmpBase);
				print "BuildTree completed with placement pending; rerun with -continue 1 to retry placement only\n";
				exit(0);
			}
			my $placements = $epaResult->{placements};
			my $backboneTreeText = readEpaFilterBackboneTree($backboneTree);
			my $backboneGraftQC = map_epa_placements_to_backbone(
				$epaResult->{tree}, $backboneTreeText, $placements);
			my $backboneGraftReport =
				"$treeD/strict_backbone.epa_backbone_grafts.tsv";
			writeEpaBackboneGraftAudit($backboneGraftQC, $backboneGraftReport);
			printEpaBackboneGraftSummary($backboneGraftQC, $backboneGraftReport);
			my $placementQC = filter_epa_placement_outliers(
				$backboneTreeText, $placements,
				{
					pendant_outlier_factor => $epaPendantOutlierFactor,
					pendant_minimum_threshold => $epaPendantMinThreshold,
					outgroup => $outgroup,
				},
			);
			my $filterSummary = "$treeD/strict_backbone.epa_filter_summary.tsv";
			writeEpaPlacementFilterSummary($placementQC, $filterSummary);
			printEpaPlacementFilterSummary($placementQC, $filterSummary, $report);
			my $reportFh = retry_open('>', $report, label => "write EPA-ng placement report");
			print {$reportFh} join("\t", @placementReportColumns), "\n";
			for my $sample (sort keys %{$placements}) {
				my $entry = $placements->{$sample};
				my $overlapMetric =
					$strictSplit->{backbone_overlap}{$sample} || {};
				my @overlapValues = map {
					defined($overlapMetric->{$_})
						? sprintf('%.12g', $overlapMetric->{$_}) : 'NA'
				} qw(backbone_overlap_nt backbone_overlap_loci
					backbone_state_divergence);
				print {$reportFh} join("\t",
					$sample, $entry->{status}, @overlapValues,
					map({ defined($entry->{$_}) ? sprintf('%.12g', $entry->{$_}) : 'NA' }
						qw(edge likelihood likelihood_weight_ratio edpl candidate_placements distal_length backbone_distal_length pendant_length)),
					defined($entry->{pendant_outlier_limit})
						? sprintf('%.12g', $entry->{pendant_outlier_limit}) : 'NA',
					$entry->{placement_filter_reason} // '',
					$strictSplit->{reason}{$sample} // '',
				), "\n";
			}
			retry_close($reportFh, "close EPA-ng placement report");
			if (!$dedicatedBackbone) {
				$primaryTree =~ s/\.treefile$/.placed.treefile/;
				$primaryTree .= ".placed.treefile" if $primaryTree eq $backboneTree;
			}
			my $publicationOK = eval {
				write_epa_placed_tree($backboneTreeText, $primaryTree, $placements);
				1;
			};
			if (!$publicationOK) {
				my $error = $@ || 'unknown EPA-ng tree-publication failure';
				$error =~ s/\s+\z//;
				clearLifecycleMarker($completionMarker, 'clear stale tree completion');
				writeOutcomeMarker($placementPendingMarker, 'placement_pending', $error, {
					backbone_tree => $backboneTree,
					query_samples => scalar(@{$strictSplit->{placement}}),
					jplace => $jplaceFile,
				}, $outD);
				warn "EPA-ng placement publication deferred; the validated backbone, jplace, and compressed MSA were retained: $error\n";
				finalizeMSAArtifacts($MsaD, $MsaWorkD);
				safeRemoveTree($tmpD, $tmpBase);
				print "BuildTree completed with placement pending; rerun with -continue 1 to retry placement publication\n";
				exit(0);
			}
			print "EPA-ng ML placements: $report; jplace: $jplaceFile; model: $modelArtifact; "
				."primary tree: $primaryTree; backbone tree: $backboneTree\n";
		} else {
			my $reportFh = retry_open('>', $report, label => "write empty EPA-ng placement report");
			print {$reportFh} join("\t", @placementReportColumns), "\n";
			retry_close($reportFh, "close empty EPA-ng placement report");
			if ($dedicatedBackbone) {
				my $temporaryPrimary = "$primaryTree.tmp.$$";
				retry_unlink($temporaryPrimary, label => "clear primary-tree temporary");
				retry_operation(
					label => "copy backbone tree to primary-tree temporary",
					code => sub { copy($backboneTree, $temporaryPrimary) && -s $temporaryPrimary },
				);
				retry_rename($temporaryPrimary, $primaryTree,
					label => "publish primary tree $primaryTree");
			}
			print "No samples required EPA-ng placement; primary tree: $primaryTree; "
				."backbone tree: $backboneTree\n";
		}
		if ($primaryTree ne $backboneTree) {
			${$trRetH}{backbone_nwk} = $backboneTree;
			${$trRetH}{nwk} = $primaryTree;
			$phyloTree = $primaryTree;
		}
	} else {
		warn "Strict-backbone samples were separated, but no completed backbone tree "
			."was available for post-inference placement\n";
	}
}
postAlignmentStep("EPA-ng placement and tree publication", $placementStarted,
	"enabled=".($strictSplit ? 1 : 0),
	"placement_samples=".($strictSplit ? scalar(@{$strictSplit->{placement}}) : 0),
	"primary_tree=".(${$trRetH}{nwk} // '<none>'));
#system "rm -f $multAli.ph $multAliSyn.ph $multAliNonSyn.ph";

if ($useTreeShrink){
	my $trShr = getProgPaths("treeshrink");
	my $inputTree = ${$trRetH}{nwk} // "";
	die "TreeShrink requested but no completed tree is available\n" unless $inputTree ne "" && -s $inputTree;
	my $cmd = "$trShr -i $outD -t ".shellQuote($inputTree)." -q 0.05  -O TS. -f";
	print "Running TreeShrink on $inputTree\n";
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
finalizeMSAArtifacts($MsaD, $MsaWorkD);
if ($gzipInput){
	# Release input caches before a compression-time ordered rewrite.  Only
	# buildTree-owned plain-to-gzip conversions are sorted; existing .gz inputs
	# and FASTA files with other names are left untouched.
	%FNA = ();
	%FAA = ();
	for my $inputFile ($aaFna, $fnFna, $cogCats){
		next if $inputFile eq "" || $inputFile =~ /\.gz$/ || !-f $inputFile;
		my $inputBasename = basename($inputFile);
		sortFastaForCompression($inputFile)
			if $inputBasename eq "allFAAs.faa" || $inputBasename eq "allFNAs.fna";
		systemW("$pigzBin -p $ncore ".shellQuote($inputFile));
	}
}

safeRemoveTree($tmpD, $tmpBase);
clearLifecycleMarker($terminalMarker, 'clear obsolete terminal no-tree marker');
clearLifecycleMarker($placementPendingMarker, 'clear completed placement-pending marker');
writeCompletionMarker($completionMarker, ${$trRetH}{nwk}, $outD)
	if length($completionMarker);
cleanupLegacyBuildTreeStateFiles();
compactTaxonAwareDiagnostics();
writeWorkflowHeartbeat('complete');
	###################### ETE ######################3

print "BuildTree completed successfully\n";
print "Outputs: alignments=$MsaD; trees=$treeD; run directory=$outD\n";
exit(0);













##########################################################################################
##########################################################################################





sub requestedTreeMethods{
	return grep { $_->{enabled} } (
		{name => "FastTree",     enabled => $doFastTree,     outputKey => "fastTrOut"},
		{name => "VeryFastTree", enabled => $doVeryFastTree, outputKey => "VfastTrOut"},
		{name => "IQ-TREE",      enabled => $doIQTree,       outputKey => "IQtreeout", iqtree => 1},
		{name => "RAxML-NG",     enabled => $doRAXMLng,      outputKey => "RAXNGtreeout"},
		{name => "RAxML",        enabled => $doRAXML,        outputKey => "RAXtreeout"},
	);
}

sub cachedIQTreeOutputComplete {
	my ($prefix, $alignment, $reasonRef) = @_;
	my $fingerprint = join("\t", $prefix, $alignment,
		inputFingerprint("$prefix.treefile"),
		inputFingerprint("$prefix.log"), inputFingerprint($alignment));
	if (exists $iqtreeValidationCache{$fingerprint}) {
		my ($complete, $reason) = @{$iqtreeValidationCache{$fingerprint}};
		${$reasonRef} = $reason if ref($reasonRef) eq 'SCALAR';
		return $complete;
	}
	my $reason = '';
	my $complete = iqtreeOutputComplete($prefix, $alignment, \$reason);
	$iqtreeValidationCache{$fingerprint} = [$complete, $reason];
	${$reasonRef} = $reason if ref($reasonRef) eq 'SCALAR';
	return $complete;
}

sub treeMethodState{
	my ($method, $hr) = @_;
	my $output = $method->{iqtree}
		? "$hr->{$method->{outputKey}}.treefile"
		: $hr->{$method->{outputKey}};
	my $validationReason = '';
	my $outputComplete = $method->{iqtree}
		? cachedIQTreeOutputComplete($hr->{$method->{outputKey}}, $hr->{inMSA}, \$validationReason)
		: (-s $output ? 1 : 0);
	return {
		%{$method},
		output => $output,
		outputComplete => $outputComplete,
		checkpointComplete => ($continue && $outputComplete ? 1 : 0),
		validationReason => $validationReason,
	};
}

sub treePresent{
	my ($hr) = @_;
	my @states = map { treeMethodState($_, $hr) } requestedTreeMethods();
	return @states && !grep { !$_->{checkpointComplete} } @states;
}



sub createTreeOpt{
	my ($multF,$siteTag,$tcnt,$silent,$consTree) = @_;
	#$siteTag="allsites";
	my $isSubTree = 0;
	my $outgroupL = $outgroup;
	$isSubTree = 1 if ($tcnt ne "");
	$outgroupL = "" if ($isSubTree);
	my $partiF=$multF.$partiExt;
	# mergeMSAs creates the tree input and its partition sidecar on scratch.
	# Its existence is resolved only immediately before a tree program is invoked.
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
					iqMemMB => $iqMemMB,
					iqPathogen => $iqPathogen,
					iqLegacy => $iqLegacy,
					cont => 0,
					restartIncomplete => $continue,
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

sub epaResourcePlan {
	my ($requestedThreads, $availableCores, $configuredMemoryMB, $treeMemoryMB,
		$memoryFraction, $memoryPerThreadMB) = @_;
	$memoryFraction = 0.60 unless defined $memoryFraction;
	$memoryPerThreadMB = 1024 unless defined $memoryPerThreadMB;
	my $threads = $requestedThreads < $availableCores
		? $requestedThreads : $availableCores;
	my $memoryMB = $configuredMemoryMB;
	$memoryMB = int($treeMemoryMB * $memoryFraction)
		if $memoryMB < 0 && $treeMemoryMB > 0;
	$memoryMB = 0 if $memoryMB < 0;
	if ($memoryMB > 0) {
		my $memoryThreads = int($memoryMB / $memoryPerThreadMB);
		$memoryThreads = 1 if $memoryThreads < 1;
		$threads = $memoryThreads if $memoryThreads < $threads;
	}
	return ($threads, $memoryMB);
}

sub readStrictBackboneClassification {
	my ($path) = @_;
	die "EPA recovery requires a non-empty strict-backbone classification: $path\n"
		unless -s $path;
	open my $input, '<', $path
		or die "Cannot read strict-backbone classification $path: $!\n";
	my $header = <$input> // '';
	$header =~ s/[\r\n]+\z//;
	my @columns = split /\t/, $header, -1;
	my %column = map { $columns[$_] => $_ } 0 .. $#columns;
	for my $required (qw(sample tree_role reason backbone_overlap_nt
		backbone_overlap_loci backbone_state_divergence)) {
		die "Strict-backbone classification $path lacks column '$required'\n"
			unless exists $column{$required};
	}
	my $split = {
		backbone => [], placement => [], excluded => [],
		reason => {}, backbone_overlap => {},
	};
	my %seen;
	while (my $line = <$input>) {
		$line =~ s/[\r\n]+\z//;
		next unless length $line;
		my @field = split /\t/, $line, -1;
		my $sample = $field[$column{sample}] // '';
		my $role = $field[$column{tree_role}] // '';
		die "Malformed strict-backbone classification row in $path\n"
			unless length($sample) && grep { $role eq $_ } qw(backbone placement excluded);
		die "Duplicate sample '$sample' in strict-backbone classification $path\n"
			if $seen{$sample}++;
		push @{$split->{$role}}, $sample;
		$split->{reason}{$sample} = $field[$column{reason}] // '';
		for my $metric (qw(backbone_overlap_nt backbone_overlap_loci
			backbone_state_divergence)) {
			my $value = $field[$column{$metric}] // 'NA';
			$split->{backbone_overlap}{$sample}{$metric} = 0 + $value
				if $value ne 'NA' && $value =~ /\A-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\z/;
		}
	}
	close $input or die "Cannot close strict-backbone classification $path: $!\n";
	return $split;
}

sub runEpaOnlyPlacement {
	my ($treeOpts, $backboneAlignment, $queryAlignment, $treeDirectory,
		$classificationFile) = @_;
	die "EPA-only recovery refuses to run without placement-pending state: $placementPendingMarker\n"
		unless -s $placementPendingMarker;
	open my $marker, '<', $placementPendingMarker
		or die "Cannot read placement-pending marker $placementPendingMarker: $!\n";
	my %markerState;
	while (my $line = <$marker>) {
		$line =~ s/[\r\n]+\z//;
		my ($key, $value) = split /\t/, $line, 2;
		$markerState{$key} = $value if defined($key) && defined($value);
	}
	close $marker or die "Cannot close placement-pending marker $placementPendingMarker: $!\n";
	die "EPA-only recovery requires status=placement_pending in $placementPendingMarker\n"
		unless ($markerState{status} // '') eq 'placement_pending';
	die "EPA-only recovery refuses to replace a completed result: $completionMarker\n"
		if -s $completionMarker;
	die "EPA-only recovery requires retained backbone alignment $backboneAlignment\n"
		unless -s $backboneAlignment;
	die "EPA-only recovery requires retained placement alignment $queryAlignment\n"
		unless -s $queryAlignment;

	my $backboneTree = "$treeOpts->{IQtreeout}.treefile";
	my $validationReason = '';
	die "EPA-only recovery requires a validated IQ-TREE backbone $backboneTree: $validationReason\n"
		unless cachedIQTreeOutputComplete(
			$treeOpts->{IQtreeout}, $backboneAlignment, \$validationReason);
	my $split = readStrictBackboneClassification($classificationFile);
	die "EPA-only recovery found no placement samples in $classificationFile\n"
		unless @{$split->{placement}};
	my $primaryTree = $backboneTree;
	die "EPA-only recovery expected a dedicated .backbone.treefile: $backboneTree\n"
		unless $primaryTree =~ s/\.backbone\.treefile\z/.treefile/;
	my $report = "$treeDirectory/strict_backbone.epa_placements.tsv";
	my @reportColumns = qw(
		sample status backbone_overlap_nt backbone_overlap_loci
		backbone_state_divergence edge likelihood likelihood_weight_ratio
		edpl candidate_placements distal_length backbone_distal_length pendant_length
		pendant_outlier_limit placement_filter_reason reason
	);

	writeWorkflowHeartbeat('EPA-only placement');
	my ($epaResult, $modelArtifact, $jplaceFile);
	my $placementOK = eval {
		($epaResult, $modelArtifact, $jplaceFile) = runEpaNgPlacement(
			$treeOpts, $backboneTree, $backboneAlignment, $queryAlignment,
			$split->{placement}, $treeDirectory,
		);
		1;
	};
	if (!$placementOK) {
		my $error = $@ || 'unknown EPA-ng placement failure';
		$error =~ s/\s+\z//;
		writeOutcomeMarker($placementPendingMarker, 'placement_pending', $error, {
			backbone_tree => $backboneTree,
			query_samples => scalar(@{$split->{placement}}),
			retry_mode => 'epa_only',
		}, $outD);
		warn "EPA-only placement remains pending; the backbone and MSA were not modified: $error\n";
		safeRemoveTree($tmpD, $tmpBase);
		writeWorkflowHeartbeat('placement_pending');
		return 0;
	}

	my $placements = $epaResult->{placements};
	my $backboneTreeText = readEpaFilterBackboneTree($backboneTree);
	my $backboneGraftQC = map_epa_placements_to_backbone(
		$epaResult->{tree}, $backboneTreeText, $placements);
	my $backboneGraftReport =
		"$treeDirectory/strict_backbone.epa_backbone_grafts.tsv";
	writeEpaBackboneGraftAudit($backboneGraftQC, $backboneGraftReport);
	printEpaBackboneGraftSummary($backboneGraftQC, $backboneGraftReport);
	my $placementQC = filter_epa_placement_outliers(
		$backboneTreeText, $placements,
		{
			pendant_outlier_factor => $epaPendantOutlierFactor,
			pendant_minimum_threshold => $epaPendantMinThreshold,
			outgroup => $outgroup,
		},
	);
	my $filterSummary = "$treeDirectory/strict_backbone.epa_filter_summary.tsv";
	writeEpaPlacementFilterSummary($placementQC, $filterSummary);
	printEpaPlacementFilterSummary($placementQC, $filterSummary, $report);
	my $reportHandle = retry_open('>', $report,
		label => 'write EPA-only placement report');
	print {$reportHandle} join("\t", @reportColumns), "\n";
	for my $sample (sort keys %{$placements}) {
		my $entry = $placements->{$sample};
		my $overlap = $split->{backbone_overlap}{$sample} || {};
		print {$reportHandle} join("\t",
			$sample, $entry->{status},
			map({ defined($overlap->{$_}) ? sprintf('%.12g', $overlap->{$_}) : 'NA' }
				qw(backbone_overlap_nt backbone_overlap_loci backbone_state_divergence)),
			map({ defined($entry->{$_}) ? sprintf('%.12g', $entry->{$_}) : 'NA' }
				qw(edge likelihood likelihood_weight_ratio edpl candidate_placements
					distal_length backbone_distal_length pendant_length)),
			defined($entry->{pendant_outlier_limit})
				? sprintf('%.12g', $entry->{pendant_outlier_limit}) : 'NA',
			$entry->{placement_filter_reason} // '',
			$split->{reason}{$sample} // '',
		), "\n";
	}
	retry_close($reportHandle, 'close EPA-only placement report');
	write_epa_placed_tree($backboneTreeText, $primaryTree, $placements);
	die "EPA-only placement did not publish its primary tree: $primaryTree\n"
		unless -s $primaryTree;

	finalizeMSAArtifacts($MsaD, $MsaWorkD);
	writeCompletionMarker($completionMarker, $primaryTree, $outD);
	clearLifecycleMarker($terminalMarker, 'clear obsolete terminal no-tree marker');
	clearLifecycleMarker($placementPendingMarker, 'clear completed placement-pending marker');
	cleanupLegacyBuildTreeStateFiles();
	safeRemoveTree($tmpD, $tmpBase);
	compactTaxonAwareDiagnostics();
	writeWorkflowHeartbeat('complete');
	print "EPA-only recovery completed; primary tree=$primaryTree; "
		."backbone retained=$backboneTree; jplace=$jplaceFile; model=$modelArtifact\n";
	return 1;
}

sub runRedoEpaFilter {
	my ($treeDirectory, $classificationFile) = @_;
	die "Forced EPA filter redo requires an existing tree directory: $treeDirectory\n"
		unless defined($treeDirectory) && -d $treeDirectory;
	my $backboneTree = File::Spec->catfile(
		$treeDirectory, 'IQtree_allsites.backbone.treefile');
	my $primaryTree = File::Spec->catfile(
		$treeDirectory, 'IQtree_allsites.treefile');
	my $jplaceFile = File::Spec->catfile(
		$treeDirectory, 'epa-ng', 'epa_result.jplace');
	die "Forced EPA filter redo requires retained backbone $backboneTree\n"
		unless -s $backboneTree;
	die "Forced EPA filter redo requires retained jplace $jplaceFile\n"
		unless -s $jplaceFile;

	my $split = readStrictBackboneClassification($classificationFile);
	die "Forced EPA filter redo found no placement samples in $classificationFile\n"
		unless @{$split->{placement}};
	writeWorkflowHeartbeat('redo EPA filtering');
	clearLifecycleMarker($completionMarker,
		'clear completion before forced EPA filter redo');
	clearLifecycleMarker($placementPendingMarker,
		'clear stale EPA-only state before forced EPA filter redo');

	my $epaResult = read_epa_jplace($jplaceFile, $split->{placement});
	my $placements = $epaResult->{placements};
	my $backboneTreeText = readEpaFilterBackboneTree($backboneTree);
	my $backboneGraftQC = map_epa_placements_to_backbone(
		$epaResult->{tree}, $backboneTreeText, $placements);
	my $backboneGraftReport =
		"$treeDirectory/strict_backbone.epa_backbone_grafts.tsv";
	writeEpaBackboneGraftAudit($backboneGraftQC, $backboneGraftReport);
	printEpaBackboneGraftSummary($backboneGraftQC, $backboneGraftReport);

	my $placementQC = filter_epa_placement_outliers(
		$backboneTreeText, $placements,
		{
			pendant_outlier_factor => $epaPendantOutlierFactor,
			pendant_minimum_threshold => $epaPendantMinThreshold,
			outgroup => $outgroup,
		},
	);
	my $report = "$treeDirectory/strict_backbone.epa_placements.tsv";
	my $filterSummary =
		"$treeDirectory/strict_backbone.epa_filter_summary.tsv";
	writeEpaPlacementFilterSummary($placementQC, $filterSummary);
	printEpaPlacementFilterSummary($placementQC, $filterSummary, $report);

	my @reportColumns = qw(
		sample status backbone_overlap_nt backbone_overlap_loci
		backbone_state_divergence edge likelihood likelihood_weight_ratio
		edpl candidate_placements distal_length backbone_distal_length pendant_length
		pendant_outlier_limit placement_filter_reason reason
	);
	my $reportHandle = retry_open('>', $report,
		label => 'write forced EPA placement report');
	print {$reportHandle} join("\t", @reportColumns), "\n";
	for my $sample (sort keys %{$placements}) {
		my $entry = $placements->{$sample};
		my $overlap = $split->{backbone_overlap}{$sample} || {};
		print {$reportHandle} join("\t",
			$sample, $entry->{status},
			map({ defined($overlap->{$_})
				? sprintf('%.12g', $overlap->{$_}) : 'NA' }
				qw(backbone_overlap_nt backbone_overlap_loci
					backbone_state_divergence)),
			map({ defined($entry->{$_})
				? sprintf('%.12g', $entry->{$_}) : 'NA' }
				qw(edge likelihood likelihood_weight_ratio edpl
					candidate_placements distal_length
					backbone_distal_length pendant_length)),
			defined($entry->{pendant_outlier_limit})
				? sprintf('%.12g', $entry->{pendant_outlier_limit}) : 'NA',
			$entry->{placement_filter_reason} // '',
			$split->{reason}{$sample} // '',
		), "\n";
	}
	retry_close($reportHandle, 'close forced EPA placement report');

	retry_unlink($primaryTree,
		label => 'remove superseded EPA-placed tree before publication')
		if -e $primaryTree;
	write_epa_placed_tree($backboneTreeText, $primaryTree, $placements);
	die "Forced EPA filter redo did not publish its primary tree: $primaryTree\n"
		unless -s $primaryTree;
	finalizeMSAArtifacts($MsaD, $MsaWorkD);
	writeCompletionMarker($completionMarker, $primaryTree, $outD);
	clearLifecycleMarker($terminalMarker, 'clear obsolete terminal no-tree marker');
	cleanupLegacyBuildTreeStateFiles();
	safeRemoveTree($tmpD, $tmpBase);
	compactTaxonAwareDiagnostics();
	writeWorkflowHeartbeat('complete');
	print "Forced EPA filter redo completed without alignment or inference; "
		."primary tree=$primaryTree; backbone retained=$backboneTree; "
		."jplace=$jplaceFile\n";
	return 1;
}

sub readEpaFilterBackboneTree {
	my ($backboneTree) = @_;
	die "EPA placement filtering requires a persisted backbone tree\n"
		unless defined($backboneTree) && -s $backboneTree;
	my $backboneHandle = retry_open('<', $backboneTree,
		label => 'read EPA placement-filter backbone tree');
	local $/;
	my $backboneText = <$backboneHandle>;
	retry_close($backboneHandle, 'close EPA placement-filter backbone tree');
	die "EPA placement-filter backbone tree is empty: $backboneTree\n"
		unless defined($backboneText) && $backboneText =~ /\S/;
	return $backboneText;
}

sub writeEpaBackboneGraftAudit {
	my ($graftQC, $auditFile) = @_;
	die "EPA backbone-graft audit requires mapping metrics\n"
		unless ref($graftQC) eq 'HASH'
			&& ref($graftQC->{rows}) eq 'ARRAY';
	die "EPA backbone-graft audit requires an output path\n"
		unless defined($auditFile) && length($auditFile);
	my @columns = qw(
		edge edge_type terminal split_size jplace_length backbone_length
		difference changed placement_count clamped_placement_count
	);
	my $temporary = "$auditFile.tmp.$$";
	retry_unlink($temporary, fatal => 0,
		label => 'clear EPA backbone-graft audit temporary');
	my $auditHandle = retry_open('>', $temporary,
		label => 'write EPA backbone-graft audit');
	print {$auditHandle} join("\t", @columns), "\n";
	for my $row (@{$graftQC->{rows}}) {
		print {$auditHandle} join("\t",
			map {
				my $value = $row->{$_};
				$_ eq 'edge_type' || $_ eq 'terminal'
					? (defined($value) ? $value : '')
					: epaFilterMetricValue($value)
			} @columns
		), "\n" or die "Cannot write EPA backbone-graft audit $temporary: $!\n";
	}
	retry_close($auditHandle, 'close EPA backbone-graft audit');
	retry_rename($temporary, $auditFile,
		label => 'publish EPA backbone-graft audit');
	return $auditFile;
}

sub printEpaBackboneGraftSummary {
	my ($graftQC, $auditFile) = @_;
	print "EPA-ng backbone graft mapping: matched "
		.epaFilterMetricValue($graftQC->{compared_edge_count})
		." jplace edges to the original backbone; differing jplace lengths="
		.epaFilterMetricValue($graftQC->{different_length_count})
		.", missing jplace lengths="
		.epaFilterMetricValue($graftQC->{jplace_length_missing_count})
		.", zero-length backbone edges="
		.epaFilterMetricValue($graftQC->{zero_backbone_edge_count})
		."; mapped placements="
		.epaFilterMetricValue($graftQC->{mapped_placement_count})
		.", clamped attachment coordinates="
		.epaFilterMetricValue($graftQC->{clamped_placement_count})
		."; max absolute difference="
		.epaFilterMetricValue($graftQC->{max_absolute_difference})
		."; publication template=original backbone; details=$auditFile\n";
}

sub epaFilterMetricValue {
	my ($value) = @_;
	return 'NA' unless defined $value;
	return sprintf('%.12g', $value)
		if $value =~ /^[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?$/;
	return $value;
}

sub writeEpaPlacementFilterSummary {
	my ($placementQC, $summaryFile) = @_;
	die "EPA placement-filter summary requires filter metrics\n"
		unless ref($placementQC) eq 'HASH';
	die "EPA placement-filter summary requires an output path\n"
		unless defined($summaryFile) && length($summaryFile);
	my $temporary = "$summaryFile.tmp.$$";
	retry_unlink($temporary, fatal => 0,
		label => 'clear EPA placement-filter summary temporary');
	my $summaryHandle = retry_open('>', $temporary,
		label => 'write EPA placement-filter summary');
	print {$summaryHandle} "metric\tvalue\n";
	my @metrics = (
		['filter_enabled', $placementQC->{enabled} ? 1 : 0],
		['placed_queries_evaluated', $placementQC->{placed_query_count}],
		['retained_queries', scalar(@{$placementQC->{retained} // []})],
		['excluded_queries', scalar(@{$placementQC->{excluded} // []})],
		['backbone_terminal_branch_count', $placementQC->{backbone_terminal_count}],
		['backbone_terminal_branch_q95', $placementQC->{backbone_q95}],
		['pendant_outlier_factor', $placementQC->{factor}],
		['pendant_minimum_threshold', $placementQC->{minimum_threshold}],
		['scaled_backbone_q95_threshold', $placementQC->{scaled_threshold}],
		['applied_pendant_threshold', $placementQC->{threshold}],
		['threshold_source', $placementQC->{threshold_source}],
		['query_pendant_length_count', $placementQC->{query_pendant_count}],
		['query_pendant_length_missing_count', $placementQC->{query_pendant_missing_count}],
		['query_pendant_length_min', $placementQC->{query_pendant_min}],
		['query_pendant_length_median', $placementQC->{query_pendant_median}],
		['query_pendant_length_q95', $placementQC->{query_pendant_q95}],
		['query_pendant_length_max', $placementQC->{query_pendant_max}],
		['excluded_samples', join(',', @{$placementQC->{excluded} // []}) || 'none'],
	);
	for my $metric (@metrics) {
		print {$summaryHandle} $metric->[0], "\t", epaFilterMetricValue($metric->[1]), "\n"
			or die "Cannot write EPA placement-filter summary $temporary: $!\n";
	}
	retry_close($summaryHandle, 'close EPA placement-filter summary');
	retry_rename($temporary, $summaryFile,
		label => 'publish EPA placement-filter summary');
	return $summaryFile;
}

sub printEpaPlacementFilterSummary {
	my ($placementQC, $summaryFile, $detailReport) = @_;
	if (!$placementQC->{enabled}) {
		print "EPA-ng pendant-branch QC: disabled; summary=$summaryFile; details=$detailReport\n";
		return;
	}
	print "EPA-ng pendant-branch QC: evaluated "
		.epaFilterMetricValue($placementQC->{placed_query_count})
		." placement queries; retained ".scalar(@{$placementQC->{retained}})
		.", excluded ".scalar(@{$placementQC->{excluded}})
		."; backbone terminals=".epaFilterMetricValue($placementQC->{backbone_terminal_count})
		.", Q95=".epaFilterMetricValue($placementQC->{backbone_q95})
		."; query pendant min/median/Q95/max="
		.join('/', map { epaFilterMetricValue($placementQC->{$_}) }
			qw(query_pendant_min query_pendant_median query_pendant_q95 query_pendant_max))
		."; cutoff=max(floor ".epaFilterMetricValue($placementQC->{minimum_threshold})
		.", factor ".epaFilterMetricValue($placementQC->{factor})
		." x Q95=".epaFilterMetricValue($placementQC->{scaled_threshold})
		.")=".epaFilterMetricValue($placementQC->{threshold})
		." [".epaFilterMetricValue($placementQC->{threshold_source})."]"
		."; summary=$summaryFile; details=$detailReport\n";
	if (@{$placementQC->{excluded}}) {
		my @excluded = @{$placementQC->{excluded}};
		my $preview_count = @excluded > 20 ? 20 : scalar(@excluded);
		my $preview = join(',', @excluded[0 .. $preview_count - 1]);
		$preview .= " (+".(@excluded - $preview_count)." more)"
			if @excluded > $preview_count;
		print "EPA-ng pendant-branch QC excluded: "
			.$preview."\n";
	}
}

sub iqtreeExplicitEpaModel {
	my ($model, $text) = @_;
	return '' unless defined($model) && defined($text) && $model =~ /^GTR(?:\+|\z)/i;
	my $modifiers = uc($model);
	$modifiers =~ s/^GTR//;
	$modifiers =~ s/\+(?:FQ|FO|F)//g;
	$modifiers =~ s/\+I//g;
	$modifiers =~ s/\+G\d+//g;
	return '' if length $modifiers;
	my $number = qr/[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?/;
	my @iqtree3Compact;
	while ($text =~ /^\s*\d+\s+(\S+)\s+\S+\s+GTR\{([^}]*)\}\+F\{([^}]*)\}(?:\+I\{([^}]*)\})?(?:\+G(\d+)\{([^}]*)\})?\s*$/mig) {
		next unless uc($1) eq uc($model);
		push @iqtree3Compact, {
			rates => $2, frequencies => $3, invariant => $4,
			gamma_categories => $5, gamma_alpha => $6,
		};
	}
	# IQ-TREE 3 writes compact fitted parameters in its substitution-process
	# table. EPA-ng has a single model partition, so never select one row from
	# a multi-partition table here.
	return '' if @iqtree3Compact > 1;
	if (@iqtree3Compact) {
		my $entry = $iqtree3Compact[0];
		my @rates = $entry->{rates} =~ /($number)/g;
		return '' unless @rates == 5 || @rates == 6;
		push @rates, 1 if @rates == 5; # IQ-TREE omits its normalized G-T rate.
		return '' if grep { $_ <= 0 } @rates;
		my @frequencies = $entry->{frequencies} =~ /($number)/g;
		return '' unless @frequencies == 4;
		return '' if grep { $_ <= 0 } @frequencies;
		my $frequencyTotal = 0;
		$frequencyTotal += $_ for @frequencies;
		return '' unless $frequencyTotal > 0;
		my $descriptor = 'GTR{'.join('/', map { sprintf('%.12g', $_) } @rates)
			.'}+FU{'.join('/', map { sprintf('%.12g', $_ / $frequencyTotal) } @frequencies).'}';
		if ($model =~ /\+I(?:\+|\z)/i) {
			return '' unless defined($entry->{invariant}) && $entry->{invariant} =~ /^$number\z/;
			my $invariant = 0 + $entry->{invariant};
			return '' unless $invariant >= 0 && $invariant < 1;
			$descriptor .= '+I{'.sprintf('%.12g', $invariant).'}';
		}
		if ($model =~ /\+G(\d+)(?:\+|\z)/i) {
			my $categories = 0 + $1;
			return '' unless $categories > 0
				&& defined($entry->{gamma_categories})
				&& $entry->{gamma_categories} == $categories
				&& defined($entry->{gamma_alpha}) && $entry->{gamma_alpha} =~ /^$number\z/;
			my $alpha = 0 + $entry->{gamma_alpha};
			return '' unless $alpha > 0;
			$descriptor .= '+G'.$categories.'{'.sprintf('%.12g', $alpha).'}';
		}
		return $descriptor;
	}
	my %rate;
	while ($text =~ /^\s*([ACGT])\s*-\s*([ACGT])\s*:\s*($number)\s*$/mg) {
		my ($left, $right, $value) = (uc($1), uc($2), 0 + $3);
		my $pair = index('ACGT', $left) < index('ACGT', $right)
			? "$left$right" : "$right$left";
		$rate{$pair} = $value if $value > 0;
	}
	my @rateOrder = qw(AC AG AT CG CT GT);
	if (grep { !exists $rate{$_} } @rateOrder) {
		while ($text =~ /^\s*(?:Rate parameters?|Substitution rates?)(?:\s*\([^)]*\))?\s*:\s*([^\r\n]+)$/mig) {
			my @values = $1 =~ /($number)/g;
			next unless @values == @rateOrder;
			@rate{@rateOrder} = map { 0 + $_ } @values;
			last;
		}
	}
	return '' if grep { !exists $rate{$_} } @rateOrder;
	my %frequency;
	while ($text =~ /pi\s*\(\s*([ACGT])\s*\)\s*=\s*($number)/ig) {
		my ($state, $value) = (uc($1), 0 + $2);
		$frequency{$state} = $value if $value > 0;
	}
	my @stateOrder = qw(A C G T);
	if (grep { !exists $frequency{$_} } @stateOrder) {
		while ($text =~ /^\s*(?:Base|State) frequencies(?:\s*\([^)]*\))?\s*:\s*([^\r\n]+)$/mig) {
			my @values = $1 =~ /($number)/g;
			next unless @values == @stateOrder;
			@frequency{@stateOrder} = map { 0 + $_ } @values;
			last;
		}
	}
	return '' if grep { !exists $frequency{$_} } @stateOrder;
	my $frequencyTotal = 0;
	$frequencyTotal += $frequency{$_} for @stateOrder;
	return '' unless $frequencyTotal > 0;
	my $descriptor = 'GTR{'
		.join('/', map { sprintf('%.12g', $rate{$_}) } @rateOrder).'}+FU{'
		.join('/', map { sprintf('%.12g', $frequency{$_} / $frequencyTotal) } @stateOrder).'}';
	if ($model =~ /\+I(?:\+|\z)/i) {
		return '' unless $text =~ /^\s*Proportion of invariable sites:\s*($number)\s*$/mi;
		my $invariant = 0 + $1;
		return '' unless $invariant >= 0 && $invariant < 1;
		$descriptor .= '+I{'.sprintf('%.12g', $invariant).'}';
	}
	if ($model =~ /\+G(\d+)(?:\+|\z)/i) {
		my $categories = 0 + $1;
		return '' unless $categories > 0
			&& $text =~ /^\s*Gamma shape alpha:\s*($number)\s*$/mi;
		my $alpha = 0 + $1;
		return '' unless $alpha > 0;
		$descriptor .= '+G'.$categories.'{'.sprintf('%.12g', $alpha).'}';
	}
	return $descriptor;
}

sub iqtreePlacementModel {
	my ($prefix) = @_;
	my @reports;
	for my $file ("$prefix.iqtree", "$prefix.log") {
		next unless -s $file;
		open my $handle, '<', $file
			or die "Cannot read IQ-TREE model output $file: $!\n";
		my $text = do { local $/; <$handle> };
		close $handle or die "Cannot close IQ-TREE model output $file: $!\n";
		push @reports, { file => $file, text => $text };
	}
	my $combinedText = join("\n", map { $_->{text} } @reports);
	for my $pattern (
		qr/^\s*Model of substitution:\s*([^\s,;]+)/mi,
		qr/^\s*Best-fit model(?: according to [^:]+)?:\s*([^\s,;]+)/mi,
		qr/^\s*(?:Substitution model|Model):\s*([^\s,;]+)/mi,
		qr/(?:^|\s)-m\s+['"]?([A-Za-z0-9_.+{}=-]+)['"]?/m,
		qr/^\s*\d+\s+([A-Za-z][A-Za-z0-9]*(?:\+[A-Za-z0-9]+)*)\s+\S+\s+[A-Za-z][A-Za-z0-9]*\{/mi,
	) {
		for my $report (@reports) {
			next unless $report->{text} =~ $pattern;
			my $model = $1;
			next if $model =~ /^(?:TEST|AUTO|MFP(?:\+MERGE)?)$/i;
			my $explicit = iqtreeExplicitEpaModel($model, $combinedText);
			return $explicit if length $explicit;
			my $partitionRows = () = $combinedText =~ /^\s*\d+\s+GTR(?:\+[A-Za-z0-9]+)*\s+\S+\s+GTR\{/mig;
			my $warning = $model =~ /^GTR(?:\+|\z)/i && $partitionRows > 1
				? "Warning: IQ-TREE selected $model and reported $partitionRows fitted GTR parameter sets "
					."for separate partitions in $prefix.iqtree and $prefix.log. EPA-ng accepts "
					."one model for a concatenated alignment, so BuildTree will not silently use "
					."one partition's rates; using the generic symbolic descriptor instead. EPA-ng "
					."does not refit its parameters (apart from empirical base frequencies for +F).\n"
				: $model =~ /^GTR(?:\+|\z)/i
				? "Warning: IQ-TREE selected $model, but BuildTree could not parse its "
					."complete fitted GTR rates, base frequencies, and rate-heterogeneity "
					."parameters from $prefix.iqtree and $prefix.log; using the generic symbolic "
					."descriptor instead. EPA-ng does not refit missing parameters (apart from "
					."empirical base frequencies for +F).\n"
				: "Warning: BuildTree could not serialize fitted IQ-TREE parameters for "
					."$model from $prefix.iqtree and $prefix.log; using the generic symbolic "
					."descriptor instead. EPA-ng does not refit missing model parameters.\n";
			warn $warning;
			return $model;
		}
	}
	die "Cannot determine the fitted IQ-TREE model from $prefix.iqtree or "
		."$prefix.log for EPA-ng placement\n";
}

sub iqtreeGtrPartitionCount {
	my ($prefix) = @_;
	my $report = "$prefix.iqtree";
	return 0 unless -s $report;
	open my $handle, '<', $report
		or die "Cannot read IQ-TREE substitution-process report $report: $!\n";
	my $text = do { local $/; <$handle> };
	close $handle or die "Cannot close IQ-TREE substitution-process report $report: $!\n";
	return () = $text =~ /^\s*\d+\s+GTR(?:\+[A-Za-z0-9]+)*\s+\S+\s+GTR\{/mig;
}

sub epaModelArtifact {
	my ($treeOpts, $backboneTree, $backboneAlignment) = @_;
	my $iqtree = "$treeOpts->{IQtreeout}.treefile";
	if ($backboneTree eq $iqtree) {
		return epaRefitIqtreeModel($treeOpts, $backboneTree, $backboneAlignment)
			if iqtreeGtrPartitionCount($treeOpts->{IQtreeout}) > 1;
		return iqtreePlacementModel($treeOpts->{IQtreeout});
	}
	my $raxmlng = $treeOpts->{RAXNGtreeout};
	if ($backboneTree eq $raxmlng) {
		my $model = $raxmlng;
		$model =~ s/\.[^.]+$//;
		$model .= '.bestModel';
		return $model if -s $model;
		die "EPA-ng placement requires the completed RAxML-NG model file $model to "
			."reuse its fitted model parameters\n";
	}
	my $raxml = $treeOpts->{RAXtreeout};
	if ($backboneTree eq $raxml) {
		my $report = $raxml;
		$report =~ s/\.[^.]+$/.raxml.info/;
		return $report if -s $report;
		die "EPA-ng placement requires the completed RAxML fitted-model report $report "
			."to reuse its fitted model parameters\n";
	}
	die "EPA-ng strict-backbone placement supports a matching IQ-TREE, RAxML-NG, or RAxML "
		."backbone model; no reusable model input is available for $backboneTree\n";
}

sub epaRefitIqtreeModel {
	my ($treeOpts, $backboneTree, $backboneAlignment) = @_;
	die "EPA-ng's one-model IQ-TREE refit requires a non-empty backbone alignment: "
		."$backboneAlignment\n" unless -s $backboneAlignment;
	die "EPA-ng's one-model IQ-TREE refit requires a non-empty backbone tree: "
		."$backboneTree\n" unless -s $backboneTree;
	die "EPA-ng's one-model IQ-TREE refit is only defined for nucleotide backbones\n"
		if $treeOpts->{useAA};
	my $refitPrefix = "$treeOpts->{IQtreeout}.epa_model";
	my $marker = "$refitPrefix.inputs";
	my $fingerprint = join("\t", 'epa-ng-iqtree-single-model-refit-v1', 'GTR+F+G2',
		inputFingerprint($backboneAlignment), inputFingerprint($backboneTree));
	my $savedFingerprint = '';
	if (-s $marker) {
		my $markerHandle = retry_open('<', $marker,
			label => "read IQ-TREE EPA-ng model-refit marker");
		$savedFingerprint = <$markerHandle> // '';
		retry_close($markerHandle, "close IQ-TREE EPA-ng model-refit marker");
		chomp $savedFingerprint;
	}
	my $validationReason = '';
	my $reusable = $savedFingerprint eq $fingerprint && -s "$refitPrefix.iqtree"
		&& cachedIQTreeOutputComplete($refitPrefix, $backboneAlignment, \$validationReason);
	if ($reusable) {
		my $model = eval { iqtreePlacementModel($refitPrefix) };
		if (!$@ && $model =~ /^GTR\{/) {
			cleanupIQTreeTransients($refitPrefix);
			print "Reusing fixed-topology, unpartitioned IQ-TREE EPA-ng model refit: "
				."$refitPrefix.iqtree\n";
			return $model;
		}
		warn "Cached IQ-TREE EPA-ng model refit cannot provide a complete explicit GTR "
			."descriptor; rebuilding it.\n";
	}
	print "IQ-TREE EPA-ng model refit: estimating one unpartitioned GTR+F+G2 model "
		."on the fixed retained backbone topology\n";
	my %refitOpts = %{$treeOpts};
	$refitOpts{inMSA} = $backboneAlignment;
	$refitOpts{IQtreeout} = $refitPrefix;
	$refitOpts{partition} = '';
	$refitOpts{constraintTree} = '';
	$refitOpts{fixedTree} = $backboneTree;
	$refitOpts{bootStrap} = 0;
	$refitOpts{useAA} = 0;
	$refitOpts{iqtreeFast} = 0;
	$refitOpts{autoModel} = 0;
	$refitOpts{iqPathogen} = 0;
	$refitOpts{iqLegacy} = 0;
	$refitOpts{runSafe} = 0;
	runQItree(\%refitOpts);
	$validationReason = '';
	die "IQ-TREE completed the EPA-ng model refit, but its output is incomplete: "
		."$validationReason\n"
		unless -s "$refitPrefix.iqtree"
			&& cachedIQTreeOutputComplete($refitPrefix, $backboneAlignment, \$validationReason);
	my $model = iqtreePlacementModel($refitPrefix);
	die "IQ-TREE EPA-ng model refit did not provide a complete explicit GTR descriptor\n"
		unless $model =~ /^GTR\{/;
	my $temporaryMarker = "$marker.tmp.$$";
	my $markerHandle = retry_open('>', $temporaryMarker,
		label => "write IQ-TREE EPA-ng model-refit marker");
	print {$markerHandle} "$fingerprint\n"
		or die "Cannot write IQ-TREE EPA-ng model-refit marker $temporaryMarker: $!\n";
	retry_close($markerHandle, "close IQ-TREE EPA-ng model-refit marker");
	retry_rename($temporaryMarker, $marker,
		label => "publish IQ-TREE EPA-ng model-refit marker");
	print "EPA-ng uses fixed-topology, unpartitioned IQ-TREE model refit: "
		."$refitPrefix.iqtree\n";
	return $model;
}

sub runEpaNgPlacement {
	my ($treeOpts, $backboneTree, $backboneAlignment, $queryAlignment,
		$queries, $treeDirectory) = @_;
	die "EPA-ng placement requires a non-empty backbone alignment: $backboneAlignment\n"
		unless -s $backboneAlignment;
	die "EPA-ng placement requires a non-empty query alignment: $queryAlignment\n"
		unless -s $queryAlignment;
	my $epaNg = getProgPaths("epa-ng", 0);
	die "Strict-backbone placement requested, but epa-ng is not configured. "
		."Set epa-ng in the selected MATAFILER configuration.\n"
		unless defined($epaNg) && length($epaNg);
	my $placementModel = epaModelArtifact($treeOpts, $backboneTree, $backboneAlignment);
	my $epaDirectory = "$treeDirectory/epa-ng";
	safeRemoveTree($epaDirectory, $treeDirectory) if -d $epaDirectory || -l $epaDirectory;
	make_path($epaDirectory);
	my ($placementThreads, $placementMemoryBudgetMB) =
		epaResourcePlan($epaThreads, $ncore, $epaMaxMemMB, $iqMemMB,
			$EPA_NG_DEFAULT{memory_fraction}, $EPA_NG_DEFAULT{memory_per_thread_mb});
	my @command = (
		'--ref-msa', shellQuote($backboneAlignment),
		'--tree', shellQuote($backboneTree),
		'--query', shellQuote($queryAlignment),
		'--outdir', shellQuote($epaDirectory),
		'-m', shellQuote($placementModel),
		'--threads', $placementThreads,
		'--chunk-size', $EPA_NG_DEFAULT{chunk_size},
	);
	my $command = $epaNg;
	$command =~ s/\s+\z//;
	$command .= " ".join(' ', @command) . "\n";
	print STDERR "EPA-ng command: $command";
	print "Running EPA-ng ML placement with model $placementModel; "
		."threads=$placementThreads; planning memory="
		.($placementMemoryBudgetMB ? "${placementMemoryBudgetMB}MB" : "disabled")
		."; hard limit=scheduler/cgroup; query chunk=$EPA_NG_DEFAULT{chunk_size}\n";
	systemW($command);
	my $jplace = "$epaDirectory/epa_result.jplace";
	die "EPA-ng completed without its expected placement file $jplace\n" unless -s $jplace;
	my $result = read_epa_jplace($jplace, $queries);
	return ($result, $placementModel, $jplace);
}

#core routine to calculte (start) phylo reconstruction
sub treeAtHeart{
	my ($hr) = @_;
	my %treeOpts = %{$hr};
	my $consTree = $treeOpts{constraintTree} // ""; my $multF = $treeOpts{inMSA};
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
	my @treeMethods = requestedTreeMethods();
	my %treeState = map {
		my $state = treeMethodState($_, \%treeOpts);
		($state->{name} => $state);
	} @treeMethods;

	if ($doFastTree){
		unless ($treeState{"FastTree"}{checkpointComplete}){
			runFasttree($treeOpts{inMSA},$treeOpts{fastTrOut},$treeOpts{useAA},$treeOpts{ncore});
		}
		$phyloTree = $treeState{"FastTree"}{output};
	}
	if ($doVeryFastTree){
		unless ($treeState{"VeryFastTree"}{checkpointComplete}){
			runVeryFasttree($treeOpts{inMSA},$treeOpts{VfastTrOut},$treeOpts{useAA},$treeOpts{ncore});
		}
		$phyloTree = $treeState{"VeryFastTree"}{output};
	}
	if ($doIQTree){
		my $state = $treeState{"IQ-TREE"};
		my $IQtree = $treeOpts{IQtreeout};
		unless ($state->{checkpointComplete}){
			print "IQ-TREE checkpoint will be rebuilt/resumed: $state->{validationReason}\n"
				if $continue && (-e "$IQtree.treefile" || -e "$IQtree.log");
			runQItree(\%treeOpts);
		} else {
			cleanupIQTreeTransients($IQtree);
		}
		my $postRunState = treeMethodState($state, \%treeOpts);
		die "IQ-TREE output failed post-run validation: $postRunState->{validationReason}\n"
			unless $postRunState->{outputComplete};
		$phyloTree = $postRunState->{output};
	}
	if ($doRAXMLng){
		unless ($treeState{"RAxML-NG"}{checkpointComplete}){
			runRaxMLng(\%treeOpts);
		}
		$phyloTree = $treeState{"RAxML-NG"}{output};
	}
	if ($doRAXML){
		unless ($treeState{"RAxML"}{checkpointComplete}){
			my $f = $treeOpts{inMSA};
			my $fasta2phylip = getProgPaths("fasta2phylip_scr");
			my $tcmd = "rm -f $f.ph*; $fasta2phylip -c 50 $f > $f.ph\n";
			systemW $tcmd;
			$treeOpts{inMSA} = "$multF.ph";
			die "Can't find nonempty expected *.ph file: $multF.ph"
				unless -s $treeOpts{inMSA};
			# runRaxML's continuation logic historically keys on path existence. Do
			# not let a zero-byte published tree suppress recovery of partial work.
			unlink $treeOpts{RAXtreeout}
				or die "Cannot remove empty RAxML tree $treeOpts{RAXtreeout}: $!\n"
				if -e $treeOpts{RAXtreeout} && !-s $treeOpts{RAXtreeout};
			runRaxML(\%treeOpts);
		}
		$phyloTree = $treeState{"RAxML"}{output};
	}
	if ($doCFML && !$treeOpts{isSubTree}){
		my $outDG = File::Spec->catdir($outD, "clonalFrameML");
		make_path($outDG) unless -d $outDG;
		$outDG = File::Spec->catfile($outDG, "CFML");
		my $CFMLbin = requireConfiguredTool("clonalframeml", "ClonalFrameML");
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
	print "Single-locus mode: no gene-category file supplied\n";
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
		runMSAFix($multAli, $maxGapPerCol);
		($multAliSyn, $multAliNonSyn) = synPosOnly($multAli,$tmpOutMSAaa,0,"",$calcSyn,$calcNonSyn);

		#system "rm $tmpInMSA $fnFna $tmpOutMSAaa";
		unlink $tmpOutMSAaa or die "Cannot remove temporary alignment $tmpOutMSAaa: $!\n" if -e $tmpOutMSAaa;
	} else {
		move($tmpOutMSAaa, $multAli) or die "Cannot move $tmpOutMSAaa to $multAli: $!\n";
	}
	#$multAli = $tmpOutMSA; $multAliSyn = $tmpOutMSAsyn;
	
	print "Single-locus alignment ready: $multAli\n";
	
	return $multAli;#,\@MSAs,\@MSA_AA);

}

#### fastgear ##

sub FastGear{
	my ($fastgearSummaryBin, $fastgearReorderBin, $matlabBin);
	if ($doFastGearSummary){
		$fastgearSummaryBin = requireConfiguredTool("fastgearSummary", "fastGEAR summary");
		$fastgearReorderBin = requireConfiguredTool("fastgearReorder", "fastGEAR reorder");
		$matlabBin = requireConfiguredTool("fastgearMatlab", "fastGEAR MATLAB runtime");
	}

	if($doFastGear){
#		open I,"<$cogCats" or die "Can't open cogcats $cogCats\n"; 
		my ($xI,$ST)= gzipopen($cogCats,"CogCATs phylo");

		my $cnt3 = 0;
		while (<$xI>){
			chomp; my @splF = split /\t/;
			@splF = grep !/^NA$/, @splF;#remove NAs
			if (@splF ==0){
				limitedWarn("empty fastGEAR category",
					"Ignoring empty fastGEAR category at line " . ($cnt3 + 1) . "\n");
				$cnt3++;
				next;
			}
			my @splF2 = parseSeqId($splF[0], "fastGEAR category line ".($cnt3 + 1));
			#die "@spl\n";
			push(@geneListF, $splF2[1]);
			$cnt3 ++;
		}
		close $xI;

		my $MsaDF1 = $MsaWorkD;
		my $MsaDF2 = "$outD_clust/MSA_FG";
		make_path($outD_clust);
		systemW("cp -r ".shellQuote($MsaDF1)." ".shellQuote($MsaDF2));
		

		my $fastgearDone = 0;
		foreach my $geneF (@geneListF){
			if ($excludedLoci{$geneF}) {
				limitedWarn("excluded fastGEAR locus",
					"Warning: skipping previously excluded locus $geneF in fastGEAR\n");
				next;
			}
			my $gene_file_stem = geneFileStem($geneF);
			my $outFG = "$outD/fastGear/fastGear_Results/$gene_file_stem";
			make_path($outFG) unless -d $outFG;
			my $outFileFG = "$outFG/${gene_file_stem}_res.mat";
			my @gene_msas = sort glob("$MsaDF2/$gene_file_stem.*.fna");
			unless (@gene_msas) {
				limitedWarn("missing fastGEAR locus alignment",
					"Warning: skipping fastGEAR locus $geneF because no MSA files were found in $MsaDF2\n");
				next;
			}
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
			my $FGparFile = requireConfiguredTool("fastgearParam", "fastGEAR parameter file");
			my $fastgearOK = eval {
				runFastgear($gene_file_stem, $outFileFG, $MsaDF2, $FGparFile);
				1;
			};
			if (!$fastgearOK) {
				my $error = $@ || "unknown fastGEAR failure";
				$error =~ s/\s+$//;
				limitedWarn("failed fastGEAR locus",
					"Warning: skipping failed fastGEAR locus $geneF: $error\n");
				next;
			}
			$fastgearDone++;
			print "fastGEAR progress: $fastgearDone/" . scalar(@geneListF) . " loci\n"
				if $fastgearDone == 1 || $fastgearDone % 25 == 0;
		}
		print "fastGEAR summary: $fastgearDone/" . scalar(@geneListF) . " loci completed\n";
		safeRemoveTree($outD_clust, $tmpD);
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
		print "fastGEAR recombination summary completed: $summaryD\n";
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
	print "Distance matrix written: $opID\n";
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
	print "Distance matrix written: $opID\n";
	#die "done\n";
}

sub readPostAlignmentRateMetrics {
	my ($reportFile) = @_;
	return {} unless defined($reportFile) && -s $reportFile;
	open my $report, '<', $reportFile
		or die "Cannot read post-alignment locus-QC report $reportFile: $!\n";
	my $header = <$report>;
	die "Post-alignment locus-QC report is empty: $reportFile\n"
		unless defined $header;
	$header =~ s/[\r\n]+$//;
	my @columns = split /\t/, $header, -1;
	my %columnIndex = map { $columns[$_] => $_ } 0 .. $#columns;
	for my $required (qw(
		alignment status p90_consensus_divergence
		called_cells gc_cells gc_fraction effective_sites
	)) {
		die "Post-alignment locus-QC report $reportFile lacks column '$required'; "
			."deterministic rate/GC merging requires MSAfix v2.14 or later\n"
			unless exists $columnIndex{$required};
	}
	my %metrics;
	while (my $line = <$report>) {
		$line =~ s/[\r\n]+$//;
		next unless length($line);
		my @fields = split /\t/, $line, -1;
		next unless ($fields[$columnIndex{status}] // '') eq 'PASS';
		my $path = $fields[$columnIndex{alignment}] // '';
		my ($rate, $called, $gc, $effective) = @fields[@columnIndex{qw(
			p90_consensus_divergence called_cells gc_fraction effective_sites
		)}];
		next unless length($path)
			&& defined($rate) && defined($called) && defined($gc) && defined($effective)
			&& $rate =~ /\A(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?\z/
			&& $called =~ /\A\d+(?:\.\d*)?\z/
			&& $gc =~ /\A(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?\z/
			&& $effective =~ /\A(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?\z/;
		$metrics{$path} = {
			value => 0 + $rate,
			source => 'p90_consensus_divergence',
			called_cells => 0 + $called,
			gc_fraction => 0 + $gc,
			effective_sites => 0 + $effective,
		};
	}
	close $report
		or die "Cannot close post-alignment locus-QC report $reportFile: $!\n";
	return \%metrics;
}

sub partitionLocusRangeCount {
	my ($partitionFile) = @_;
	return 0 unless defined($partitionFile) && length($partitionFile);
	my ($partition, $opened) =
		gzipopen($partitionFile, 'retained alignment partition', 0, 0);
	return 0 unless $opened && $partition;
	my $ranges = 0;
	while (my $line = <$partition>) {
		$line =~ s/[\r\n]+\z//;
		next unless $line =~ /=\s*(.+)\z/;
		my $coordinates = $1;
		$ranges++ while $coordinates =~ /(?:\A|,\s*)\d+\s*-\s*\d+(?=\s*(?:,|\z))/g;
	}
	close $partition
		or die "Cannot close retained alignment partition $partitionFile: $!\n";
	return $ranges;
}

sub deterministicRatePartitions {
	my ($loci) = @_;
	die "Deterministic partition merging received no loci\n" unless @{$loci};
	my @rescueLoci = grep { ($_->{selection_phase} // '') eq 'taxon_rescue' } @{$loci};
	my @binningLoci = grep { ($_->{selection_phase} // '') ne 'taxon_rescue' } @{$loci};
	unless (@binningLoci) {
		@binningLoci = @{$loci};
		@rescueLoci = ();
	}
	my @rates = map { $_->{rate_proxy} } @{$loci};
	my @gcFractions = map { $_->{gc_fraction} } @{$loci};

	my $summarize = sub {
		my ($members) = @_;
		my ($alignmentSites, $effectiveSites, $rate, $gc) = (0, 0, 0, 0);
		for my $member (@{$members}) {
			$alignmentSites += $member->{length};
			$effectiveSites += $member->{effective_sites};
			$rate += $member->{rate_proxy};
			$gc += $member->{gc_fraction};
		}
		my $count = scalar @{$members};
		return {
			loci => $count,
			sites => $effectiveSites,
			alignment_sites => $alignmentSites,
			mean_rate => $count ? $rate / $count : 0,
			mean_gc => $count ? $gc / $count : 0,
		};
	};
	my ($minimumRate, $maximumRate) = (sort { $a <=> $b } @rates)[0, -1];
	my ($minimumGC, $maximumGC) = (sort { $a <=> $b } @gcFractions)[0, -1];
	my $rateRange = $maximumRate - $minimumRate;
	my $gcRange = $maximumGC - $minimumGC;
	my $totalEffectiveSites = 0;
	$totalEffectiveSites += $_->{effective_sites} for @{$loci};
	my $desiredBins = int(($totalEffectiveSites + $rateMergeTargetSites - 1)
		/ $rateMergeTargetSites);
	$desiredBins = 1 if $desiredBins < 1;
	$desiredBins = $rateMergeMaxBins if $desiredBins > $rateMergeMaxBins;
	$desiredBins = scalar(@binningLoci) if $desiredBins > @binningLoci;

	# Refine the largest current bin until the site-driven target is met.  Each
	# split is placed near its effective-site median and uses whichever of P90
	# divergence or GC still has the stronger local normalized spread.
	my $splitMetric = sub {
		my ($members, $field, $range) = @_;
		return undef if @{$members} < 2 || !$range;
		my @ordered = sort {
			$a->{$field} <=> $b->{$field}
				|| $a->{start} <=> $b->{start}
		} @{$members};
		my $total = 0;
		$total += $_->{effective_sites} for @ordered;
		my ($leftSites, $bestIndex, $bestBalance);
		for my $index (0 .. $#ordered - 1) {
			$leftSites += $ordered[$index]{effective_sites};
			next if $ordered[$index]{$field} == $ordered[$index + 1]{$field};
			my $balance = abs($total - 2 * $leftSites);
			if (!defined($bestBalance) || $balance < $bestBalance) {
				($bestIndex, $bestBalance) = ($index, $balance);
			}
		}
		return undef unless defined $bestIndex;
		my @left = @ordered[0 .. $bestIndex];
		my @right = @ordered[$bestIndex + 1 .. $#ordered];
		my ($minimum, $maximum) = ($ordered[0]{$field}, $ordered[-1]{$field});
		return {
			field => $field,
			left => \@left,
			right => \@right,
			balance => $bestBalance,
			span => ($maximum - $minimum) / $range,
		};
	};
	my $bestSplit = sub {
		my ($members) = @_;
		my @candidates = grep { defined } (
			$splitMetric->($members, 'rate_proxy', $rateRange),
			$splitMetric->($members, 'gc_fraction', $gcRange),
		);
		return undef unless @candidates;
		return (sort {
			$b->{span} <=> $a->{span}
				|| $a->{balance} <=> $b->{balance}
				|| $a->{field} cmp $b->{field}
		} @candidates)[0];
	};
	my %bins = (root => \@binningLoci);
	my $splitSerial = 0;
	while (keys(%bins) < $desiredBins) {
		my @splitCandidates;
		for my $key (keys %bins) {
			my $split = $bestSplit->($bins{$key});
			push @splitCandidates, { key => $key, split => $split,
				summary => $summarize->($bins{$key}) } if $split;
		}
		last unless @splitCandidates;
		my $candidate = (sort {
			$b->{summary}{sites} <=> $a->{summary}{sites}
				|| $b->{split}{span} <=> $a->{split}{span}
				|| $a->{key} cmp $b->{key}
		} @splitCandidates)[0];
		my $axis = $candidate->{split}{field} eq 'rate_proxy' ? 'p90' : 'gc';
		my $leftKey = $candidate->{key}."_${axis}".(++$splitSerial).'L';
		my $rightKey = $candidate->{key}."_${axis}".$splitSerial.'H';
		$_->{initial_bin} = $leftKey for @{$candidate->{split}{left}};
		$_->{initial_bin} = $rightKey for @{$candidate->{split}{right}};
		delete $bins{$candidate->{key}};
		$bins{$leftKey} = $candidate->{split}{left};
		$bins{$rightKey} = $candidate->{split}{right};
	}
	$_->{initial_bin} //= 'root' for @binningLoci;
	for my $locus (@rescueLoci) {
		my %summary = map { $_ => $summarize->($bins{$_}) } keys %bins;
		my ($target) = sort {
			my $distanceA = ($rateRange
				? abs($locus->{rate_proxy} - $summary{$a}{mean_rate}) / $rateRange : 0)
				+ ($gcRange
					? abs($locus->{gc_fraction} - $summary{$a}{mean_gc}) / $gcRange : 0);
			my $distanceB = ($rateRange
				? abs($locus->{rate_proxy} - $summary{$b}{mean_rate}) / $rateRange : 0)
				+ ($gcRange
					? abs($locus->{gc_fraction} - $summary{$b}{mean_gc}) / $gcRange : 0);
			$distanceA <=> $distanceB || $a cmp $b;
		} keys %bins;
		$locus->{initial_bin} = 'taxon_rescue_to_'.$target;
		push @{$bins{$target}}, $locus;
	}

	while (keys(%bins) > 1) {
		my %summary = map { $_ => $summarize->($bins{$_}) } keys %bins;
		my @undersized = sort {
			$summary{$a}{loci} <=> $summary{$b}{loci}
				|| $summary{$a}{sites} <=> $summary{$b}{sites}
				|| $a cmp $b
		} grep {
			$summary{$_}{loci} < $rateMergeMinLoci
				|| $summary{$_}{sites} < $rateMergeMinSites
		} keys %bins;
		last unless @undersized;
		my $source = $undersized[0];
		my @targets = grep { $_ ne $source } keys %bins;
		my ($target) = sort {
			my $distanceA = ($rateRange
				? abs($summary{$source}{mean_rate} - $summary{$a}{mean_rate}) / $rateRange : 0)
				+ ($gcRange
					? abs($summary{$source}{mean_gc} - $summary{$a}{mean_gc}) / $gcRange : 0);
			my $distanceB = ($rateRange
				? abs($summary{$source}{mean_rate} - $summary{$b}{mean_rate}) / $rateRange : 0)
				+ ($gcRange
					? abs($summary{$source}{mean_gc} - $summary{$b}{mean_gc}) / $gcRange : 0);
			$distanceA <=> $distanceB || $a cmp $b;
		} @targets;
		push @{$bins{$target}}, @{$bins{$source}};
		delete $bins{$source};
	}

	my %summary = map { $_ => $summarize->($bins{$_}) } keys %bins;
	my @orderedKeys = sort {
		$summary{$a}{mean_rate} <=> $summary{$b}{mean_rate}
			|| $summary{$a}{mean_gc} <=> $summary{$b}{mean_gc}
			|| $a cmp $b
	} keys %bins;
	my @partitions;
	for my $index (0 .. $#orderedKeys) {
		my $key = $orderedKeys[$index];
		my $name = sprintf('rateGC%02d', $index + 1);
		my @members = sort { $a->{start} <=> $b->{start} } @{$bins{$key}};
		$_->{partition} = $name for @members;
		push @partitions, {
			name => $name,
			members => \@members,
			%{$summary{$key}},
		};
	}
	return \@partitions;
}

sub writeRatePartitionAudit {
	my ($path, $loci, $partitions) = @_;
	make_path(dirname($path)) unless -d dirname($path);
	my ($audit, $temporary) = tempfile(
		'rate-merged-partitions-XXXXXX',
		DIR => dirname($path),
		UNLINK => 1,
	);
	print {$audit} join("\t", qw(
		locus alignment selection_phase start end alignment_sites effective_called_sites
		rate_proxy rate_proxy_source gc_fraction called_cells initial_bin partition
		partition_loci partition_effective_sites partition_alignment_sites
	)), "\n";
	my %partitionByName = map { $_->{name} => $_ } @{$partitions};
	for my $locus (sort { $a->{start} <=> $b->{start} } @{$loci}) {
		my $partition = $partitionByName{$locus->{partition}};
		print {$audit} join("\t",
			$locus->{locus}, $locus->{alignment}, $locus->{selection_phase},
			@{$locus}{qw(start end length)},
			sprintf('%.8g', $locus->{effective_sites}),
			sprintf('%.8g', $locus->{rate_proxy}), $locus->{rate_proxy_source},
			sprintf('%.8g', $locus->{gc_fraction}), $locus->{called_cells},
			$locus->{initial_bin}, $locus->{partition},
			$partition->{loci}, sprintf('%.8g', $partition->{sites}),
			$partition->{alignment_sites},
		), "\n";
	}
	close $audit or die "Cannot close rate-partition audit $temporary: $!\n";
	rename $temporary, $path
		or die "Cannot publish rate-partition audit $path: $!\n";
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
	my @partitionLoci;
	my $partitionCoordinate = 0;
	my $applyRateMerging = $rateMergePartitions && !$isAA && $multAliF eq $multAli;
	my $overlapFilteredLoci = 0;
	my $overlapColumnsRemoved = 0;
	foreach my $MSAf (@MSAs){
		#print $MSAf."\n"; 
		my $hit =0; my $miss =0;
		my $hr;
		my $readOK = eval {
			$hr = readFasta($MSAf,1);
			1;
		};
		if (!$readOK) {
			my $error = $@ || "unreadable alignment";
			$error =~ s/\s+$//;
			limitedWarn("invalid locus MSA",
				"Warning: excluding alignment $MSAf during merge: $error\n");
			next;
		}
		my %MFAA = %{$hr};
		unlink $MSAf or die "Cannot remove merged MSA component $MSAf: $!\n" if $del && -e $MSAf;
		my @Mkeys = sort keys %MFAA;
		next if (@Mkeys == 0);
		my ($firstMsaSample, $gcat, $separator);
		my $headerOK = eval {
			($firstMsaSample, $gcat, $separator) =
				parseSeqId($Mkeys[0], "MSA header in $MSAf");
			1;
		};
		if (!$headerOK) {
			my $error = $@ || "unparseable alignment header";
			$error =~ s/\s+$//;
			limitedWarn("invalid locus MSA",
				"Warning: excluding alignment $MSAf during merge: $error\n");
			next;
		}
		my $len = length($MFAA{$Mkeys[0]});
		if ($len == 0){
			limitedWarn("zero-length MSA", "Ignoring zero-length alignment $MSAf\n");
			$excludedLoci{$gcat} = 1;
			next;
		}
		my $originalLen = $len;
		my ($filtered, $retainedLen, $removedColumns);
		my $minimumOverlapCount = int(scalar(@Mkeys) * $minOverlapMSA + 0.999999);
		my $overlapOK = eval {
			($filtered, $retainedLen, $removedColumns) =
				filter_alignment_by_overlap(\%MFAA, $isAA, $minimumOverlapCount);
			1;
		};
		if (!$overlapOK) {
			my $error = $@ || "overlap filtering failed";
			$error =~ s/\s+$//;
			limitedWarn("invalid locus MSA",
				"Warning: excluding alignment $MSAf during merge: $error\n");
			$excludedLoci{$gcat} = 1;
			next;
		}
		%MFAA = %{$filtered};
		if ($retainedLen == 0) {
			limitedWarn("empty overlap-filtered MSA",
				"Overlap filtering removed every column from $MSAf; skipping this locus\n");
			$excludedLoci{$gcat} = 1;
			next;
		}
		$len = $retainedLen;
		my @unequal = grep { length($MFAA{$_}) != $len } @Mkeys;
		if (@unequal) {
			limitedWarn("invalid locus MSA",
				"Warning: excluding alignment $MSAf during merge because "
				.scalar(@unequal)." sequence(s) have unequal lengths\n");
			$excludedLoci{$gcat} = 1;
			next;
		}
		if ($removedColumns) {
			$overlapFilteredLoci++;
			$overlapColumnsRemoved += $removedColumns;
		}
		push(@lengthsParts,$len);
		if ($applyRateMerging) {
			my $rateMetric = $partitionRateProxy{$MSAf};
			die "No MSAfix v2.14 rate/GC metrics for retained alignment $MSAf\n"
				unless $rateMetric;
			push @partitionLoci, {
				locus => $gcat,
				alignment => $MSAf,
				start => $partitionCoordinate + 1,
				end => $partitionCoordinate + $len,
				length => $len,
				rate_proxy => $rateMetric->{value},
				rate_proxy_source => $rateMetric->{source},
				selection_phase => $partitionSelectionPhase{$MSAf} // '',
				gc_fraction => $rateMetric->{gc_fraction},
				called_cells => $rateMetric->{called_cells},
				effective_sites => $rateMetric->{effective_sites},
			};
		}
		$partitionCoordinate += $len;
		foreach my $sm (@smps){
			my $curK = $sm.$separator.$gcat;
			if ( exists($MFAA{$curK}) && defined($MFAA{$curK})  ) {
				my $seq = $MFAA{$curK}; 
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
	my @removedSeqExamples;
	my %charCnts; my $maxNtCnt=0;
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
		my $removeSample;
		if ($taxonAwareLocusSelection && $multAliF eq $multAli) {
			my $minimumAnchorNT = $ntCntTotal > $placementMinOverlap
				? $ntCntTotal : $placementMinOverlap;
			$removeSample = $maxNtCnt == 0 || ($num1 * $factor) < $minimumAnchorNT;
		} else {
			$removeSample = $maxNtCnt == 0; #|| ($num1 < ($qtl90NTcnts * $ntFrac) && $num1 < $qtl25NTcnts) || ($num1 < ($ntCntTotal / $factor));
		}
		if ($removeSample){
			delete $bigMSAFAA{$kk}; delete $bigMSAFAAnxs{$kk}; $remSeqNum++; 
			push @removedSeqExamples, "$kk($num1 informative positions)"
				if @removedSeqExamples < 5;
		}
		#print "$num1  $kk \n";#$bigMSAFAA{$kk}\n\n"; last;
	}
	open O,">$multAliF" or die "Can't open MSA outfile $multAliF\n";
	open O2,">$multAliF.nxs" or die "Can't open MSA nexus outfile $multAliF.nxs\n";
	my @allKs = sort keys %bigMSAFAA;
	if (@allKs == 0){die "no genes for nexus output format.\nAborting\n";}
	print O2 "#NEXUS\nBegin data;\nDimensions ntax=".scalar(@allKs)." nchar=".length($bigMSAFAAnxs{$allKs[0]}).";\nFormat datatype=dna missing=? gap=-;\nMatrix\n";
	foreach my $kk (@allKs){
		print O ">$kk\n"; my $s1 = $bigMSAFAA{$kk};
		print O2 "\n$kk\t"; my $s2 = $bigMSAFAAnxs{$kk};
		print O "$s1\n"; print O2 "$s2\n";
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
	if ($applyRateMerging) {
		my $partitions = deterministicRatePartitions(\@partitionLoci);
		for my $partition (@{$partitions}) {
			my @ranges = map { $_->{start}."-".$_->{end} } @{$partition->{members}};
			print O "$TypeTag, $partition->{name} = ".join(", ", @ranges)."\n";
		}
		my $auditPath = "$treeD/rate_merged_partitions.tsv";
		writeRatePartitionAudit($auditPath, \@partitionLoci, $partitions);
		print "Deterministic rate/GC partition merging: ".scalar(@partitionLoci)
			." loci -> ".scalar(@{$partitions})." partition(s); audit=$auditPath\n";
	} else {
		for (my $i=0;$i<@lengthsParts;$i++){
			#DNA, part1 = 1-100
			print O "$TypeTag, part".($i+1) ." = ". ($lastP+1) ."-". ($lengthsParts[$i]+$lastP) ."\n";
			$lastP+=$lengthsParts[$i];
		}
		my $auditPath = "$treeD/rate_merged_partitions.tsv";
		unlink $auditPath
			or die "Cannot remove stale rate-partition audit $auditPath: $!\n"
			if $multAliF eq $multAli && -e $auditPath;
	}
	close O;
	
	print "Alignment merge summary: retained " . ($iniSeqNum - $remSeqNum) . "/$iniSeqNum sequences";
	print "; removed examples: " . join(", ", @removedSeqExamples) if @removedSeqExamples;
	print "\n";
	print "Overlap filtering summary: $overlapFilteredLoci locus/loci changed; "
		. "$overlapColumnsRemoved columns removed\n"
		if $overlapFilteredLoci;
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
	#$outgroupSample is a sample name, not a sequence identifier: the MSA headers
	#are split with parseSeqId before they are compared against it.
	my ($inMSA,$inAAMSA, $ffold, $outgroupSample, $doSyn, $doNSyn) = @_;

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
	#Resolve each header to its sample once: the outgroup test below runs per
	#codon per sequence and must not re-split identifiers that many times.
	my @sampleOfSeq = map {
		(parseSeqId($_, "synonymous-site MSA header", 1))[0]
	} @aSeq;
	my %outFNA;#syn
	my %outFNAns;#non syn
	for (my $j=0;$j<@aSeq;$j++){$outFNA{$aSeq[$j]}="";}
	my $len = length ($FNA{$aSeq[0]});
	die "Nucleotide MSA $inMSA has length $len, which is not a multiple of 3; "
		."synonymous-site classification requires codon-aligned input\n" if $len % 3;
	my $nsyn=0;my $syn=0;
	for (my $i=0; $i< $len; $i+=3){ #goes over every position
		my $j =0;
		my $iniAA = "-";
		my $iniCodon ;
		#Anchor on the first callable codon, applying exactly the tolerance the
		#comparison loop below uses: gaps and any IUPAC ambiguity code are missing
		#data, and the outgroup never defines the reference amino acid.
		while ($j < @aSeq){ #check for first informative position
			if ($outgroupSample ne "" && $sampleOfSeq[$j] eq $outgroupSample){$j++; next;}
			$iniCodon = substr $FNA{$aSeq[$j]},$i,3;
			my $iniCodonUC = uc($iniCodon);
			if ($iniCodonUC !~ m/-/ && $iniCodonUC =~ m/^[ACTG]{3}$/){
				die "codon doesn't exist $iniCodon \n" unless (exists($convertor{$iniCodonUC}));
				$iniAA = $convertor{$iniCodonUC};#substr $FAA{$aSeq[0]},$i,1;
				last;
			} elsif ($iniCodonUC =~ m/^-{3}$/ || $iniCodonUC =~ m/[NWYRSKMDVHB]/){
				$j++; next;
			} else {
				die "iniCodon wrong $iniCodon\n";
			}
		}
		next if $j >= @aSeq;
		#die "$iniAA\n";
		my $isSame = 1;my $ntSame = 1;
		next unless (!$ffold || $ffd{uc($iniCodon)} == 4);
	#print $i." $iniAA ";
		for (;$j<@aSeq;$j++){
			next if $outgroupSample ne "" && $sampleOfSeq[$j] eq $outgroupSample;
			my $newCodon = uc(substr $FNA{$aSeq[$j]},$i,3);
			my $newAA = "-";
			if ($newCodon !~ m/-/ && $newCodon =~ m/^[ACTG]{3}$/){
				die "Unkown AA $newCodon\n" unless (exists $convertor{$newCodon} );
				$newAA = $convertor{$newCodon} ; # substr $FAA{$aSeq[$j]},$i,1;
			} elsif ($newCodon =~ m/^-{3}$/ || $newCodon =~ m/[NWYRSKMDVHB]/){
			} else {
				die "newCodon wrong $newCodon\n" ;
			}
			if ($iniAA ne $newAA && $newAA ne "-"){
				$isSame =0; $ntSame =0; last;
			}
			if (uc($iniCodon) ne $newCodon){
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
	$synSummaryCount++;
	$synSiteTotal += $syn;
	$nonSynSiteTotal += $nsyn;
	print "Site classification example: $reportedGene; synonymous=$syn; nonsynonymous=$nsyn; "
		. scalar(@aSeq) . " sequences\n"
		if $synSummaryCount <= 5;
	print "Further per-locus site-classification messages are suppressed; an aggregate follows alignment generation\n"
		if $synSummaryCount == 6;
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
		print "codeml summary: gene=$gene; model=$modelName; repeats=$repeatCounts; selected repeat=$repSelected\n";
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
	my $maskedInternalStops = 0;
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
					substr($seq,$x,3) = "NNN";
					$maskedInternalStops++;
				}
				$x++;
			}
		}
		print O ">$genome2\n$seq\n";
		push (@genomeList, $genome2);
	} 
	close O;
	limitedWarn("internal stop codons",
		"Masked $maskedInternalStops in-frame stop codon(s) while preparing $gene selection analysis\n")
		if $maskedInternalStops;
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
			coreHyPhy($MsaWorkD,$gene_file_stem,"",$nwkFile,$codemlOutDTmp,$logF);
		}
		my $sumTxt=$stdJSONheader;
		foreach my $gene (keys %logs){
			#my $summary = fubarXML("$logs{$gene}");
			next unless (-e "$logs{$gene}.json.gz");
			my $summary = `$jsonExrScr $logs{$gene}.json.gz`;
			next if ($summary eq "");
			$sumTxt .= "$gene\t$summary";
			#print "$sumTxt\n";
		}
		open O,">$logF1";		print O $sumTxt;		close O;
		print "Selection summary written: $logF1\n";
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
		foreach my $gene (keys %logs){
			#my $summary = fubarXML("$logs{$gene}");
			next unless (-e "$logs{$gene}.json.gz");
			my $summary = `$jsonExrScr $logs{$gene}.json.gz`;
			next if ($summary eq "");

			$sumTxt .= "$gene\t$summary";
			#system "rm -f $logs{$gene}*";
		}
		open O,">$logF1";		print O $sumTxt;		close O;
		print "Selection summary written: $logF1\n";
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
		foreach my $gene (keys %logs){
			#my $summary = fubarXML("$logs{$gene}");next if ($summary eq "");
			next unless (-e "$logs{$gene}.json.gz");
			my $summary = `$jsonExrScr $logs{$gene}.json.gz`;
			next if ($summary eq "");
			$sumTxt .= "$gene\t$summary";
			#system "rm -f $logs{$gene}*";
		}
		open O,">$logF1";		print O $sumTxt;		close O;
		print "Selection summary written: $logF1\n";
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
		$cmd .= "$hyphyBin CPU=$ncore ".getProgPaths("wattersonTheta_scr")." --alignment $MSAfile2 ";#> $logF\n";
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


sub compileSampleSeparator{
	# parseSeqId runs once per sequence per pass, so the pattern is compiled once
	# rather than re-interpolated on every call. Groups 1 and 2 precede the
	# interpolated separator and therefore keep stable numbers even when
	# -smplSep contains capture groups of its own; the gene is read from the
	# offset just past group 2 instead of a group number that would shift.
	$sampleGeneRegex = qr/^(.*?)($smplSep)(.+)$/;
}

sub parseSeqId{
	my ($seqId, $context, $allowUndelimited) = @_;
	$context ||= "sequence identifier";
	if (defined($seqId) && $seqId =~ $sampleGeneRegex){
		my ($sample, $separator) = ($1, $2);
		my $gene = substr($seqId, $+[2]);
		return ($sample, $gene, $separator) if $sample ne "" && $gene ne "";
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

sub fastaCompressionSortKey{
	my ($header) = @_;
	my ($identifier) = split /\s+/, $header, 2;
	my ($sample,$gene) = parseSeqId($identifier, "compression-sort FASTA header",1);
	return length($gene)
		? join("\t", $gene, $sample, $identifier)
		: $identifier;
}

sub sortFastaForCompression{
	my ($inputFile) = @_;
	my $records = readFasta($inputFile,0);
	die "Cannot sort empty FASTA input before compression: $inputFile\n"
		unless keys %{$records};
	my %sortKeys = map { $_ => fastaCompressionSortKey($_) } keys %{$records};
	my @headers = sort {
		$sortKeys{$a} cmp $sortKeys{$b}
			|| $a cmp $b
	} keys %{$records};
	my ($out,$tmpFile) = tempfile(
		basename($inputFile).".sort.XXXXXX",
		DIR => dirname($inputFile),
		UNLINK => 0,
	);
	for my $header (@headers){
		print {$out} ">$header\n$records->{$header}\n"
			or die "Cannot write sorted FASTA temporary file $tmpFile: $!\n";
	}
	close $out or die "Cannot close sorted FASTA temporary file $tmpFile: $!\n";
	rename $tmpFile, $inputFile
		or die "Cannot replace $inputFile with locus-sorted FASTA: $!\n";
	print "Sorted ".scalar(@headers)." records by gene/locus before compressing $inputFile\n";
}


sub shellQuote{
	my ($value) = @_;
	$value = "" unless defined $value;
	$value =~ s/'/'"'"'/g;
	return "'$value'";
}

sub informativeSequenceLength {
	my ($sequence, $useAA) = @_;
	$sequence = uc($sequence // "");
	if ($useAA) {
		$sequence =~ s/[^ACDEFGHIKLMNPQRSTVWY]//g;
	} else {
		$sequence =~ s/[^ACGTU]//g;
	}
	return length($sequence);
}

sub bestGeneSequencesBySample {
	my ($category, $charCounts, $context, $expectedGene) = @_;
	die "Best-sequence selection requires a category array\n"
		unless ref($category) eq 'ARRAY';
	die "Best-sequence selection requires per-sample character counts\n"
		unless ref($charCounts) eq 'HASH';
	$context ||= 'gene category';
	my (%bestSequence, %bestSites);
	for my $sequenceId (@{$category}) {
		my ($sample, $gene) = parseSeqId($sequenceId, $context);
		die "Wrong gene in $sequenceId, expected $expectedGene!\n"
			if defined($expectedGene) && length($expectedGene) && $gene ne $expectedGene;
		my $sites = $charCounts->{$sample}{$sequenceId} // 0;
		if (!exists($bestSites{$sample}) || $sites > $bestSites{$sample}) {
			$bestSites{$sample} = $sites;
			$bestSequence{$sample} = $sequenceId;
		}
	}
	return (\%bestSequence, \%bestSites);
}

sub selectTaxonAwareCandidateLoci {
	my %args = @_;
	my $categories = $args{categories};
	my $charCounts = $args{char_counts};
	my (%metrics, %samples);
	my $categoryIndex = -1;
	for my $category (@{$categories}) {
		$categoryIndex++;
		next unless @{$category};
		my (undef, $gene) = parseSeqId(
			$category->[0], "taxon-aware category ".($categoryIndex + 1));
		die "Duplicate locus '$gene' in taxon-aware category input\n"
			if exists $metrics{$gene};
		my ($bestSequence, $bestSites) = bestGeneSequencesBySample(
			$category, $charCounts,
			"taxon-aware category ".($categoryIndex + 1), $gene,
		);
		if ($args{use_aa}) {
			$_ *= 3 for values %{$bestSites};
		}
		next if keys(%{$bestSequence}) < 3;
		my @siteCounts = values %{$bestSites};
		my $q90 = quantile(0.9, @siteCounts);
		my $medianSites = quantile(0.5, @siteCounts);
		my @absoluteDeviations = map { abs($_ - $medianSites) } @siteCounts;
		my $mad = @absoluteDeviations ? quantile(0.5, @absoluteDeviations) : 0;
		my @completeness = map { $q90 > 0 ? ($_ / $q90 > 1 ? 1 : $_ / $q90) : 0 }
			@siteCounts;
		my $medianCompleteness = @completeness ? quantile(0.5, @completeness) : 0;
		my $lengthStability = $medianSites > 0
			? 1 - ($mad / $medianSites > 1 ? 1 : $mad / $medianSites)
			: 0;
		my @sequenceIds = map { $bestSequence->{$_} } sort keys %{$bestSequence};
		# Raw coordinates are not homologous when indels are present. Defer all
		# information scoring until after alignment and avoid this extra scan.
		my $potentialInformation = { variable_sites => 0, parsimony_informative_sites => 0 };
		$metrics{$gene} = {
			gene => $gene,
			preferred_core => catalogueGeneFromLocus($gene, $args{preferred_core_genes}) ? 1 : 0,
			category => \@sequenceIds,
			sample_sites => $bestSites,
			sample_count => scalar(keys %{$bestSequence}),
			q90_nt => $q90,
			median_completeness => $medianCompleteness,
			length_stability => $lengthStability,
			potential_variable_sites => $potentialInformation->{variable_sites},
			potential_parsimony_informative_sites => $potentialInformation->{parsimony_informative_sites},
			potential_informative_nt => $potentialInformation->{parsimony_informative_sites}
				* ($args{use_aa} ? 3 : 1),
			quality_score => 0,
		};
		for my $sample (keys %{$bestSites}) {
			$samples{$sample}{available_loci}++;
			$samples{$sample}{available_nt} += $bestSites->{$sample};
		}
	}
	unless (keys %metrics) {
		writeTaxonAwareLocusAudit($args{report}, "candidate", \%metrics);
		return {
			categories => [], metrics => {}, samples => {},
			terminal_reason => 'taxon_aware_no_category_with_three_usable_samples',
		};
	}
	my $universeSampleCount = scalar(keys %samples) || 1;
	my $presortRanked = assignPresortScores(\%metrics, $args{preferred_core_genes});
	my $presortWeight = $presortRanked ? ($args{presort_weight} // 0) : 0;
	my $measuredWeight = 1 - $presortWeight;
	for my $metric (values %metrics) {
		my $prevalence = $metric->{sample_count} / $universeSampleCount;
		$metric->{prevalence} = $prevalence;
		$metric->{potential_information_score} = 0;
		$metric->{robust_score} = $measuredWeight * (
				0.55 * $prevalence
				+ 0.30 * $metric->{median_completeness}
				+ 0.15 * $metric->{length_stability})
			+ $presortWeight * $metric->{presort_score};
		$metric->{quality_score} = $metric->{robust_score};
	}
	printf "Taxon-aware candidate scoring: %d/%d locus/loci carry a presorter rank; "
		."presort weight=%.2f\n", $presortRanked, scalar(keys %metrics), $presortWeight
		if $presortRanked;
	my $selectedGenes = chooseTaxonAwareLoci(
		metrics => \%metrics,
		limit => $args{candidate_limit},
		core_limit => $args{core_limit},
		final_limit => $args{final_limit},
		target_loci => $args{target_loci},
		target_nt => $args{target_nt},
		rescue_min_prevalence => $args{rescue_min_prevalence},
		rescue_prevalence_mode => $args{rescue_prevalence_mode},
		targets_from_gate => $args{targets_from_gate},
		gate_gene_fraction => $args{gate_gene_fraction},
		gate_nt_fraction => $args{gate_nt_fraction},
		gate_minimum_nt => $args{gate_minimum_nt},
		stage => "candidate",
	);
	my @selectedCategories = map { $metrics{$_}{category} } @{$selectedGenes};
	writeTaxonAwareLocusAudit($args{report}, "candidate", \%metrics);
	return {
		categories => \@selectedCategories,
		metrics => \%metrics,
		samples => \%samples,
	};
}

sub rawCoordinateInformation {
	my ($sequenceIds, $sequences, $useAA) = @_;
	return { variable_sites => 0, parsimony_informative_sites => 0 }
		unless ref($sequenceIds) eq 'ARRAY' && ref($sequences) eq 'HASH';
	my @sequence = map { uc($sequences->{$_} // '') } @{$sequenceIds};
	my $maximumLength = 0;
	for my $sequence (@sequence) {
		$maximumLength = length($sequence) if length($sequence) > $maximumLength;
	}
	my ($variable, $informative) = (0, 0);
	for my $position (0 .. $maximumLength - 1) {
		my %states;
		for my $sequence (@sequence) {
			# Variable raw gene lengths are expected; an absent tail is missing data.
			next if $position >= length($sequence);
			my $state = substr($sequence, $position, 1);
			next unless $useAA
				? $state =~ /^[ACDEFGHIKLMNPQRSTVWY]$/
				: $state =~ /^[ACGTU]$/;
			$states{$state}++;
		}
		$variable++ if keys(%states) >= 2;
		my $repeatedStates = grep { $_ >= 2 } values(%states);
		$informative++ if $repeatedStates >= 2;
	}
	return { variable_sites => $variable, parsimony_informative_sites => $informative };
}

sub chooseTaxonAwareLoci {
	my %args = @_;
	my $metrics = $args{metrics};
	my $rescueMinimumPrevalence = $args{rescue_min_prevalence}
		// $TAXON_AWARE_DEFAULT{rescue_minimum_prevalence};
	my $preferredCoreRescueBonus = 0.05; # secondary to actual coverage gain
	for my $metric (values %{$metrics}) {
		delete @{$metric}{qw(
			selected selection_rank selection_phase selection_objective
			coverage_rescue_eligible coverage_rescue_reason
		)};
	}
	my @eligible = sort {
		($metrics->{$b}{preferred_core} // 0) <=> ($metrics->{$a}{preferred_core} // 0)
			|| ($metrics->{$b}{quality_score} // 0) <=> ($metrics->{$a}{quality_score} // 0)
			|| ($metrics->{$b}{presort_score} // 0) <=> ($metrics->{$a}{presort_score} // 0)
			|| $a cmp $b
	} grep { ($metrics->{$_}{sample_count} // 0) >= 3 } keys %{$metrics};
	# Rescue eligibility is measured on recovered prevalence, which sequencing
	# depth caps well below 1: in a cohort with a low-coverage tail no locus can
	# reach an absolute 0.8 and the coverage phase then selects nothing at all.
	# Reading the setting against the highest prevalence actually attained keeps
	# the gate on the same scale as the data it is applied to.
	my $maximumPrevalence = 0;
	for my $gene (@eligible) {
		my $prevalence = $metrics->{$gene}{prevalence} // 0;
		$maximumPrevalence = $prevalence if $prevalence > $maximumPrevalence;
	}
	my $prevalenceMode = $args{rescue_prevalence_mode} // 'absolute';
	my $effectiveRescuePrevalence = $prevalenceMode eq 'relative'
		? $rescueMinimumPrevalence * $maximumPrevalence
		: $rescueMinimumPrevalence;
	my $rescueEligibleCount = 0;
	for my $gene (@eligible) {
		my $prevalence = $metrics->{$gene}{prevalence} // 0;
		my $broadlyAvailable = $prevalence + 1e-12 >= $effectiveRescuePrevalence;
		$metrics->{$gene}{coverage_rescue_eligible} = $broadlyAvailable ? 1 : 0;
		$metrics->{$gene}{coverage_rescue_reason} = $broadlyAvailable
			? "broad_prevalence_met" : "below_rescue_minimum_prevalence";
		$rescueEligibleCount++ if $broadlyAvailable;
	}
	my $limit = $args{limit} < @eligible ? $args{limit} : scalar(@eligible);
	my $coreLimit = $args{core_limit} < $limit ? $args{core_limit} : $limit;
	my (%selected, %sampleLoci, %sampleSites, %availability);
	for my $gene (@eligible) {
		next unless $metrics->{$gene}{coverage_rescue_eligible};
		$availability{$_}++ for keys %{$metrics->{$gene}{sample_sites}};
	}
	my $maximumAvailability = 1;
	for my $available (values %availability) {
		$maximumAvailability = $available if $available > $maximumAvailability;
	}
	my @chosen;
	for my $gene (@eligible[0 .. $coreLimit - 1]) {
		push @chosen, $gene;
		$selected{$gene} = 1;
		$metrics->{$gene}{selection_phase} = "robust_core";
		$metrics->{$gene}{selection_objective} = $metrics->{$gene}{quality_score};
		for my $sample (keys %{$metrics->{$gene}{sample_sites}}) {
			$sampleLoci{$sample}++;
			$sampleSites{$sample} += $metrics->{$gene}{sample_sites}{$sample};
		}
	}
	# Sample retention is decided on 0.9-quantile fractions of the selected
	# coverage, while the greedy phase optimises an absolute per-sample floor.
	# Left apart, the greedy phase declares a sample satisfied far below the
	# level at which it is actually kept, so project the final Q90 from the
	# completed core phase and raise the targets to the retention basis. The
	# configured targets stay a floor: this only ever tightens them.
	my ($targetLoci, $targetNT) = ($args{target_loci}, $args{target_nt});
	my ($gateTargetLoci, $gateTargetNT) = (0, 0);
	if ($args{targets_from_gate} && $coreLimit > 0 && @chosen) {
		my $projection = $limit / $coreLimit;
		my $projectedLoci = quantile(0.9, values %sampleLoci) * $projection;
		my $projectedNT = quantile(0.9, values %sampleSites) * $projection;
		$gateTargetLoci = int(($args{gate_gene_fraction} // 0) * $projectedLoci + 0.999999);
		$gateTargetNT = int(($args{gate_nt_fraction} // 0) * $projectedNT + 0.999999);
		$gateTargetNT = $args{gate_minimum_nt}
			if ($args{gate_minimum_nt} // 0) > $gateTargetNT;
		$targetLoci = $gateTargetLoci if $gateTargetLoci > $targetLoci;
		$targetNT = $gateTargetNT if $gateTargetNT > $targetNT;
	}
	while (@chosen < $limit) {
		my ($bestGene, $bestObjective);
		for my $gene (@eligible) {
			next if $selected{$gene};
			next unless $metrics->{$gene}{coverage_rescue_eligible};
			my $coverageGain = 0;
			for my $sample (keys %{$metrics->{$gene}{sample_sites}}) {
				my $rarityWeight = sqrt(
					$maximumAvailability / ($availability{$sample} || 1));
				$rarityWeight = 6 if $rarityWeight > 6;
				$coverageGain += $rarityWeight / $targetLoci
					if ($sampleLoci{$sample} // 0) < $targetLoci;
				my $siteDeficit = $targetNT - ($sampleSites{$sample} // 0);
				if ($siteDeficit > 0) {
					my $siteGain = $metrics->{$gene}{sample_sites}{$sample};
					$siteGain = $siteDeficit if $siteGain > $siteDeficit;
					$coverageGain += $rarityWeight * $siteGain / $targetNT;
				}
			}
			my $objective = $coverageGain
				+ 0.25 * ($metrics->{$gene}{quality_score} // 0)
				+ $preferredCoreRescueBonus * ($metrics->{$gene}{preferred_core} // 0);
			if (!defined($bestObjective) || $objective > $bestObjective
					|| ($objective == $bestObjective
						&& (($metrics->{$gene}{presort_score} // 0)
								<=> ($metrics->{$bestGene}{presort_score} // 0)
							|| $bestGene cmp $gene) > 0)) {
				($bestGene, $bestObjective) = ($gene, $objective);
			}
		}
		last unless defined $bestGene;
		push @chosen, $bestGene;
		$selected{$bestGene} = 1;
		my $nextRank = scalar(@chosen);
		$metrics->{$bestGene}{selection_phase} =
			$args{stage} eq "candidate"
				&& $nextRank > ($args{final_limit} // $args{limit})
			? "qc_backfill" : "taxon_rescue";
		$metrics->{$bestGene}{selection_objective} = $bestObjective;
		for my $sample (keys %{$metrics->{$bestGene}{sample_sites}}) {
			$sampleLoci{$sample}++;
			$sampleSites{$sample} += $metrics->{$bestGene}{sample_sites}{$sample};
		}
	}
	for my $index (0 .. $#chosen) {
		$metrics->{$chosen[$index]}{selected} = 1;
		$metrics->{$chosen[$index]}{selection_rank} = $index + 1;
	}
	my $stage = $args{stage} // 'selection';
	if (@chosen < $limit) {
		# The coverage phase ran out of admissible loci before the budget was
		# spent. This used to end the loop silently and forfeit the difference.
		warn "Taxon-aware $stage filled only ".scalar(@chosen)."/$limit locus "
			."slot(s): the coverage phase had $rescueEligibleCount eligible "
			."locus/loci at prevalence >= "
			.sprintf('%.4f', $effectiveRescuePrevalence)
			." ($prevalenceMode mode, setting="
			.sprintf('%.2f', $rescueMinimumPrevalence)
			.", highest observed prevalence=".sprintf('%.4f', $maximumPrevalence)
			."). Lower -taxonAwareRescueMinPrevalence, or use "
			."-taxonAwareRescuePrevalenceMode relative, to spend the whole budget.\n";
	}
	printf "Taxon-aware %s: %d/%d locus slot(s) filled; robust core=%d, "
		."coverage phase=%d (eligible=%d at prevalence >= %.4f, %s mode); "
		."per-sample targets=%d loci/%d NT%s\n",
		$stage, scalar(@chosen), $limit, ($coreLimit < 0 ? 0 : $coreLimit),
		scalar(@chosen) - ($coreLimit < 0 ? 0 : $coreLimit), $rescueEligibleCount,
		$effectiveRescuePrevalence, $prevalenceMode, $targetLoci, $targetNT,
		(($gateTargetLoci || $gateTargetNT)
			? sprintf(" (raised to the retention basis from %d/%d)",
				$args{target_loci}, $args{target_nt})
			: "");
	return \@chosen;
}

sub taxonAwareAlignmentMetrics {
	my ($alignment, $useAA, $gene, $qualificationSequences) = @_;
	my $records = readFasta($alignment, 1);
	die "Taxon-aware selector cannot read alignment $alignment\n"
		unless ref($records) eq 'HASH' && keys %{$records};
	my @sequenceIds = sort keys %{$records};
	my $alignmentLength = 0;
	for my $sequenceId (@sequenceIds) {
		my $length = length($records->{$sequenceId} // "");
		$alignmentLength = $length if $length > $alignmentLength;
	}
	my (%sampleSites, @sequences);
	my $validCells = 0;
	for my $sequenceId (@sequenceIds) {
		my $sequence = uc($records->{$sequenceId} // "");
		push @sequences, $sequence;
		my ($sample) = parseSeqId($sequenceId, "taxon-aware alignment $gene");
		my $sites = informativeSequenceLength($sequence, $useAA);
		$sites *= 3 if $useAA;
		my $qualifiesForSampleQC = !defined($qualificationSequences)
			|| $qualificationSequences->{$sequenceId};
		$sampleSites{$sample} = $sites
			if $qualifiesForSampleQC
				&& (!exists($sampleSites{$sample}) || $sites > $sampleSites{$sample});
		$validCells += $useAA ? int($sites / 3) : $sites;
	}
	my ($variableSites, $parsimonyInformativeSites) = (0, 0);
	for my $position (0 .. $alignmentLength - 1) {
		my %stateCount;
		for my $sequence (@sequences) {
			# Retain robustness for a partial external alignment: absent tails are missing data.
				next if $position >= length($sequence);
			my $state = substr($sequence, $position, 1);
			next unless $useAA
				? $state =~ /^[ACDEFGHIKLMNPQRSTVWY]$/
				: $state =~ /^[ACGTU]$/;
			$stateCount{$state}++;
		}
		$variableSites++ if keys(%stateCount) >= 2;
		my $repeatedStates = grep { $_ >= 2 } values %stateCount;
		$parsimonyInformativeSites++ if $repeatedStates >= 2;
	}
	my $occupancyDenominator = scalar(@sequenceIds) * $alignmentLength;
	return {
		gene => $gene,
		path => $alignment,
		sample_sites => \%sampleSites,
		sample_count => scalar(keys %sampleSites),
		alignment_length_nt => $alignmentLength * ($useAA ? 3 : 1),
		occupancy => $occupancyDenominator
			? $validCells / $occupancyDenominator : 0,
		variable_sites => $variableSites,
		parsimony_informative_sites => $parsimonyInformativeSites,
		information_density => $alignmentLength ? $parsimonyInformativeSites / $alignmentLength : 0,
		variable_density => $alignmentLength ? $variableSites / $alignmentLength : 0,
		called_cells => $validCells,
	};
}

sub classifyTaxonAwareSamples {
	my %args = @_;
	my $metrics = $args{metrics};
	my %result;
	for my $sample (sort keys %{$args{samples}}) {
		my ($selectedLoci, $selectedNT) = (0, 0);
		for my $metric (values %{$metrics}) {
			next if $args{selected_only} && !$metric->{selected};
			next unless ($metric->{sample_sites}{$sample} // 0) > 0;
			$selectedLoci++;
			$selectedNT += $metric->{sample_sites}{$sample};
		}
		my $available = $args{samples}{$sample};
		my $availableLoci = ref($available) eq 'HASH'
			? ($available->{available_loci} // 0) : 0;
		my $availableNT = ref($available) eq 'HASH'
			? ($available->{available_nt} // 0) : 0;
		my ($role, $reason);
		if ($selectedNT <= 0) {
			($role, $reason) = ("remove", "no_selected_anchor");
		} elsif ($selectedNT < $args{minimum_anchor_nt}) {
			($role, $reason) = ("remove", "below_minimum_anchor_nt");
		} elsif ($selectedLoci >= $args{target_loci}
				&& $selectedNT >= $args{target_nt}) {
			($role, $reason) = ("backbone_candidate", "coverage_target_met");
		} else {
			($role, $reason) = ("placement_candidate", "usable_sparse_anchor");
		}
		$result{$sample} = {
			available_loci => $availableLoci,
			available_nt => $availableNT,
			selected_loci => $selectedLoci,
			selected_nt => $selectedNT,
			role => $role,
			reason => $reason,
		};
	}
	if (($args{outgroup} // "") ne ""
			&& (!exists($result{$args{outgroup}})
				|| $result{$args{outgroup}}{role} eq "remove")) {
		die "Taxon-aware selection could not retain an aligned anchor for outgroup "
			."'$args{outgroup}'\n";
	}
	return \%result;
}

sub classifyTaxonAwareCoverageEligibility {
	my %args = @_;
	my $metrics = $args{sample_metrics};
	my $role = $args{role} // 'placement';
	die ucfirst($role)." eligibility requires final taxon-aware sample metrics\n"
		unless ref($metrics) eq 'HASH' && keys %{$metrics};
	my @selectedLoci = map { $_->{selected_loci} // 0 } values %{$metrics};
	my @selectedNT = map { $_->{selected_nt} // 0 } values %{$metrics};
	my $minimumLoci = int(quantile(0.9, @selectedLoci)
		* ($args{gene_fraction} // 0) + 0.999999);
	my $minimumLociFloor = $args{minimum_loci_floor} // 1;
	$minimumLoci = $minimumLociFloor if $minimumLoci < $minimumLociFloor;
	my $relativeNT = int(quantile(0.9, @selectedNT)
		* ($args{nt_fraction} // 0) + 0.999999);
	my $minimumNT = $args{minimum_nt} // 0;
	$minimumNT = $relativeNT if $relativeNT > $minimumNT;
	$minimumNT = $args{minimum_overlap}
		if ($args{minimum_overlap} // 0) > $minimumNT;
	my %result;
	for my $sample (sort keys %{$metrics}) {
		my $metric = $metrics->{$sample};
		my ($eligible, $reason) = (1, "${role}_coverage_met");
		if (($args{outgroup} // '') ne '' && $sample eq $args{outgroup}) {
			($eligible, $reason) = (1, 'outgroup_retained');
		} elsif (($metric->{role} // 'remove') eq 'remove') {
			($eligible, $reason) = (0, 'not_retained_after_locus_selection');
		} elsif (($metric->{selected_loci} // 0) < $minimumLoci) {
			($eligible, $reason) = (0, "below_${role}_gene_fraction");
		} elsif (($metric->{selected_nt} // 0) < $minimumNT) {
			($eligible, $reason) = (0, "below_${role}_nt_fraction");
		}
		$result{$sample} = {
			eligible => $eligible,
			reason => $reason,
			selected_loci => $metric->{selected_loci} // 0,
			selected_nt => $metric->{selected_nt} // 0,
		};
	}
	return {
		minimum_loci => $minimumLoci,
		minimum_nt => $minimumNT,
		samples => \%result,
	};
}

sub selectTaxonAwareFinalLoci {
	my %args = @_;
	my (%metrics, %seenGene);
	my $universeSampleCount = scalar(keys %{$args{samples} || {}}) || 1;
	for my $alignment (@{$args{alignments}}) {
		my $gene = $args{path_gene}{$alignment};
		die "Taxon-aware selector has no locus mapping for alignment $alignment\n"
			unless defined($gene) && length($gene);
		die "Taxon-aware selector received duplicate alignment for locus $gene\n"
			if $seenGene{$gene}++;
		my $metric = taxonAwareAlignmentMetrics(
			$alignment, $args{use_aa}, $gene, $args{qualification_sequences});
		my $preScore = $args{pre_metrics}{$gene}{robust_score} // 0;
		my $prevalence = $metric->{sample_count} / $universeSampleCount;
		my $informationDensity = $metric->{information_density} // 0;
		my $informationSupport = ($metric->{parsimony_informative_sites} // 0) / 5;
		$informationSupport = 1 if $informationSupport > 1;
		my $informationScore =
			($informationDensity / $args{information_saturation}) * $informationSupport;
		$informationScore = 1 if $informationScore > 1;
		my $variableDensity = $metric->{variable_density} // 0;
		my $excessVariationPenalty = $variableDensity > $args{excess_variation_onset}
			? ($variableDensity - $args{excess_variation_onset})
				/ $args{excess_variation_span}
			: 0;
		$excessVariationPenalty = 1 if $excessVariationPenalty > 1;
		$metric->{robust_score} = $preScore;
		$metric->{preferred_core} = ($args{pre_metrics}{$gene}{preferred_core}
			// catalogueGeneFromLocus($gene, $args{preferred_core_genes})) ? 1 : 0;
		$metric->{prevalence} = $prevalence;
		$metric->{information_score} = $informationScore;
		$metric->{information_density} = $informationDensity;
		$metric->{variable_density} = $variableDensity;
		$metric->{excess_variation_penalty} = $excessVariationPenalty;
		$metrics{$gene} = $metric;
	}
	# The presort term is normalised across the surviving loci, so scoring has to
	# wait until every metric exists.
	my $presortRanked = assignPresortScores(\%metrics, $args{preferred_core_genes});
	my $presortWeight = $presortRanked ? ($args{presort_weight} // 0) : 0;
	my $measuredWeight = 1 - $presortWeight;
	for my $metric (values %metrics) {
		$metric->{quality_score} = $measuredWeight * (
				0.30 * $metric->{robust_score}
				+ 0.25 * $metric->{occupancy}
				+ 0.25 * $metric->{prevalence}
				+ 0.20 * $metric->{information_score})
			+ $presortWeight * $metric->{presort_score}
			- 0.10 * $metric->{excess_variation_penalty};
		$metric->{quality_score} = 0 if $metric->{quality_score} < 0;
	}
	printf "Taxon-aware final scoring: %d/%d locus/loci carry a presorter rank; "
		."presort weight=%.2f; information saturates at %.4f variable-site density; "
		."excess-variation penalty from %.4f over %.4f\n",
		$presortRanked, scalar(keys %metrics), $presortWeight,
		$args{information_saturation}, $args{excess_variation_onset},
		$args{excess_variation_span};
	my $selectedGenes = chooseTaxonAwareLoci(
		metrics => \%metrics,
		limit => $args{maximum_loci},
		core_limit => $args{core_loci},
		final_limit => $args{maximum_loci},
		target_loci => $args{target_loci},
		target_nt => $args{target_nt},
		rescue_min_prevalence => $args{rescue_min_prevalence},
		rescue_prevalence_mode => $args{rescue_prevalence_mode},
		targets_from_gate => $args{targets_from_gate},
		gate_gene_fraction => $args{gate_gene_fraction},
		gate_nt_fraction => $args{gate_nt_fraction},
		gate_minimum_nt => $args{gate_minimum_nt},
		stage => "final",
	);
	die "Taxon-aware final selection found no usable locus\n"
		unless @{$selectedGenes};
	my @selectedAlignments = map { $metrics{$_}{path} } @{$selectedGenes};
	my $sampleMetrics = classifyTaxonAwareSamples(
		metrics => \%metrics,
		samples => $args{samples},
		target_loci => $args{target_loci},
		target_nt => $args{target_nt},
		minimum_anchor_nt => $args{minimum_anchor_nt},
		selected_only => 1,
		outgroup => $args{outgroup},
	);
	writeTaxonAwareLocusAudit($args{locus_report}, "final", \%metrics);
	writeTaxonAwareSampleAudit($args{sample_report}, $sampleMetrics);
	return {
		alignments => \@selectedAlignments,
		sample_metrics => $sampleMetrics,
		locus_metrics => \%metrics,
	};
}

sub writeSelectionAttritionAudit {
	my ($path, $stats) = @_;
	die "Selection attrition statistics must be a hash reference\n"
		unless ref($stats) eq 'HASH';
	my @order = qw(
		input_loci input_sequences input_samples
		qc_excluded_samples qc_excluded_sequences
		qc_emptied_loci coverage_excluded_samples length_retained_sequences
		length_filtered_sequences length_include_retained_sequences
		length_include_filtered_sequences length_recovery_candidate_sequences
		length_recovered_msa_sequences gene_length_min_dropped_loci
		gene_length_include_min_dropped_loci gene_length_recovery_candidate_loci
		gene_length_recovered_msa_loci eligible_loci candidate_loci candidate_samples
		aligned_loci alignment_failed_loci post_qc_loci final_loci final_samples
		backbone_samples placement_samples excluded_samples
	);
	my %previous;
	if (-s $path) {
		open my $existing, '<', $path
			or die "Cannot read existing selection attrition audit $path: $!\n";
		my $header = <$existing> // '';
		$header =~ s/[\r\n]+\z//;
		if ($header eq "metric\tvalue") {
			while (my $line = <$existing>) {
				$line =~ s/[\r\n]+\z//;
				my ($metric, $value) = split /\t/, $line, 2;
				next unless defined($metric) && defined($value)
					&& $metric ne 'schema' && $value ne 'NA';
				$previous{$metric} = $value;
			}
		}
		close $existing
			or die "Cannot close existing selection attrition audit $path: $!\n";
	}
	make_path(dirname($path)) unless -d dirname($path);
	my $temporary = "$path.tmp.$$";
	open my $output, '>', $temporary
		or die "Cannot write selection attrition audit $temporary: $!\n";
	print {$output} "metric\tvalue\n"
		or die "Cannot write selection attrition header $temporary: $!\n";
	print {$output} "schema\t3\n"
		or die "Cannot write selection attrition schema $temporary: $!\n";
	for my $metric (@order) {
		my $value = defined($stats->{$metric}) ? $stats->{$metric} : 'NA';
		$value = $previous{$metric}
			if $value eq 'NA' && exists($previous{$metric});
		print {$output} "$metric\t$value\n"
			or die "Cannot write selection attrition metric $metric: $!\n";
	}
	close $output or die "Cannot close selection attrition audit $temporary: $!\n";
	rename $temporary, $path
		or die "Cannot publish selection attrition audit $path: $!\n";
}

sub writeGeneLengthSampleAudit {
	my ($path, $audit, $removedSamples, $recoveredForMSA,
		$qcThreshold, $includeThreshold) = @_;
	die "Gene-length sample audit requires a sample hash\n"
		unless ref($audit) eq 'HASH';
	$removedSamples ||= {};
	$recoveredForMSA ||= {};
	my @columns = qw(
		sample gene_length_min gene_length_include_min sample_prefilter_status
		input_loci gene_length_min_pass_loci gene_length_min_dropped_loci
		qc_pass_loci qc_dropped_loci gene_length_include_min_pass_loci
		gene_length_include_min_dropped_loci recovery_candidate_loci
		recovered_for_msa_loci gene_length_min_dropped_genes qc_dropped_genes
		gene_length_include_min_dropped_genes recovery_candidate_genes
		recovered_for_msa_genes
	);
	my %summary = map { $_ => 0 } qw(
		input_loci gene_length_min_dropped_loci qc_dropped_loci
		gene_length_include_min_dropped_loci recovery_candidate_loci
		recovered_for_msa_loci
	);
	make_path(dirname($path)) unless -d dirname($path);
	my $temporary = "$path.tmp.$$";
	my $output = retry_open('>', $temporary,
		label => 'write gene-length sample audit');
	print {$output} join("\t", @columns), "\n"
		or die "Cannot write gene-length sample-audit header $temporary: $!\n";
	for my $sample (sort keys %{$audit}) {
		my @genes = sort keys %{$audit->{$sample}};
		my @lengthMinPass = grep {
			$audit->{$sample}{$_}{gene_length_min}
		} @genes;
		my @lengthMinDropped = grep {
			!$audit->{$sample}{$_}{gene_length_min}
		} @genes;
		my @qcPass = grep { $audit->{$sample}{$_}{qc} } @genes;
		my @qcDropped = grep { !$audit->{$sample}{$_}{qc} } @genes;
		my @includePass = grep { $audit->{$sample}{$_}{include} } @genes;
		my @includeDropped = grep { !$audit->{$sample}{$_}{include} } @genes;
		my @recoveryCandidate = grep {
			!$audit->{$sample}{$_}{qc} && $audit->{$sample}{$_}{include}
		} @genes;
		my @recovered = sort keys %{$recoveredForMSA->{$sample} || {}};
		my %row = (
			sample => $sample,
			gene_length_min => $qcThreshold,
			gene_length_include_min => $includeThreshold,
			sample_prefilter_status => $removedSamples->{$sample}
				? 'removed_by_high_threshold_qc' : 'retained_or_pending',
			input_loci => scalar(@genes),
			gene_length_min_pass_loci => scalar(@lengthMinPass),
			gene_length_min_dropped_loci => scalar(@lengthMinDropped),
			qc_pass_loci => scalar(@qcPass),
			qc_dropped_loci => scalar(@qcDropped),
			gene_length_include_min_pass_loci => scalar(@includePass),
			gene_length_include_min_dropped_loci => scalar(@includeDropped),
			recovery_candidate_loci => scalar(@recoveryCandidate),
			recovered_for_msa_loci => scalar(@recovered),
			gene_length_min_dropped_genes => join(',', @lengthMinDropped),
			qc_dropped_genes => join(',', @qcDropped),
			gene_length_include_min_dropped_genes => join(',', @includeDropped),
			recovery_candidate_genes => join(',', @recoveryCandidate),
			recovered_for_msa_genes => join(',', @recovered),
		);
		$summary{input_loci} += scalar(@genes);
		$summary{gene_length_min_dropped_loci} += scalar(@lengthMinDropped);
		$summary{qc_dropped_loci} += scalar(@qcDropped);
		$summary{gene_length_include_min_dropped_loci} += scalar(@includeDropped);
		$summary{recovery_candidate_loci} += scalar(@recoveryCandidate);
		$summary{recovered_for_msa_loci} += scalar(@recovered);
		print {$output} join("\t", map { $row{$_} // '' } @columns), "\n"
			or die "Cannot write gene-length sample-audit row for $sample: $!\n";
	}
	retry_close($output, 'close gene-length sample audit');
	retry_rename($temporary, $path, label => 'publish gene-length sample audit');
	return \%summary;
}

sub writeTaxonAwareLocusAudit {
	my ($path, $stage, $metrics) = @_;
	make_path(dirname($path)) unless -d dirname($path);
	open my $output, '>', $path
		or die "Cannot write taxon-aware locus audit $path: $!\n";
	print {$output} join("\t", qw(
		stage gene selected rank phase quality_score robust_score occupancy prevalence
		information_score information_density variable_density excess_variation_penalty
		preferred_core presort_rank presort_score sample_count q90_nt alignment_length_nt variable_sites
		parsimony_informative_sites potential_variable_sites
		potential_parsimony_informative_sites potential_informative_nt
		potential_information_score
		median_completeness length_stability
		selection_objective coverage_rescue_eligible coverage_rescue_reason alignment
	))."\n";
	for my $gene (sort {
		($metrics->{$a}{selected} ? 0 : 1) <=> ($metrics->{$b}{selected} ? 0 : 1)
			|| ($metrics->{$a}{selection_rank} // 1_000_000)
				<=> ($metrics->{$b}{selection_rank} // 1_000_000)
			|| $a cmp $b
	} keys %{$metrics}) {
		my $metric = $metrics->{$gene};
		my @values = (
			$stage, $gene, $metric->{selected} ? 1 : 0,
			$metric->{selection_rank} // "", $metric->{selection_phase} // "",
			map({ defined($metric->{$_}) ? sprintf("%.6f", $metric->{$_}) : "" }
				qw(quality_score robust_score occupancy prevalence information_score
					information_density variable_density excess_variation_penalty)),
			$metric->{preferred_core} ? 1 : 0,
			$metric->{presort_rank} // 0,
			defined($metric->{presort_score})
				? sprintf("%.6f", $metric->{presort_score}) : "",
			$metric->{sample_count} // "", $metric->{q90_nt} // "",
			$metric->{alignment_length_nt} // "", $metric->{variable_sites} // "",
			$metric->{parsimony_informative_sites} // "",
			$metric->{potential_variable_sites} // "",
			$metric->{potential_parsimony_informative_sites} // "",
			$metric->{potential_informative_nt} // "",
			defined($metric->{potential_information_score})
				? sprintf("%.6f", $metric->{potential_information_score}) : "",
			map({ defined($metric->{$_}) ? sprintf("%.6f", $metric->{$_}) : "" }
				qw(median_completeness length_stability selection_objective)),
			defined($metric->{coverage_rescue_eligible})
				? $metric->{coverage_rescue_eligible} : "",
			$metric->{coverage_rescue_reason} // "",
			$metric->{path} // "",
		);
		print {$output} join("\t", @values)."\n";
	}
	close $output or die "Cannot close taxon-aware locus audit $path: $!\n";
}

sub writeTaxonAwareSampleAudit {
	my ($path, $samples) = @_;
	make_path(dirname($path)) unless -d dirname($path);
	open my $output, '>', $path
		or die "Cannot write taxon-aware sample audit $path: $!\n";
	print {$output} "sample\tavailable_loci\tavailable_nt\tselected_loci\tselected_nt\trole\treason\n";
	for my $sample (sort keys %{$samples}) {
		my $metric = $samples->{$sample};
		print {$output} join("\t", $sample,
			map({ $metric->{$_} // "" } qw(
				available_loci available_nt selected_loci selected_nt role reason
			)))."\n";
	}
	close $output or die "Cannot close taxon-aware sample audit $path: $!\n";
}

sub preferredCoreGeneSet {
	my ($path) = @_;
	return {} unless defined($path) && length($path);
	open my $input, '<', $path
		or die "Cannot open preferred-core guide $path: $!\n";
	my (%genes, $ignored);
	my $rank = 0;
	while (my $line = <$input>) {
		$line =~ s/[\r\n]+\z//;
		next if $line eq '' || $line =~ /^\s*#/;
		my @field = split /\t/, $line, -1;
		if (@field < 2 || !length($field[1])) {
			$ignored++;
			next;
		}
		# Raw .core guides have one seed per row; sorter output uses the
		# same second column for a comma-separated ranked seed list. Either way
		# the file order is the presorter's importance order, so record the
		# 1-based position and keep the first (best) one seen per gene.
		for my $gene (split /,/, $field[1]) {
			$gene =~ s/^\s+|\s+$//g;
			next unless length($gene);
			next if $gene =~ /\A(?:gene|gene_id|seed)\z/i;
			$genes{$gene} = ++$rank unless exists $genes{$gene};
		}
	}
	close $input or die "Cannot close preferred-core guide $path: $!\n";
	die "Preferred-core guide $path has no usable seed genes\n" unless keys %genes;
	print "Loaded ".scalar(keys %genes)." ranked preferred universal-core seed gene(s) from $path"
		.($ignored ? "; ignored $ignored malformed row(s)" : '')."\n";
	return \%genes;
}

sub presortRankFromLocus {
	my ($locus, $preferred) = @_;
	return 0 unless defined($locus) && length($locus)
		&& ref($preferred) eq 'HASH' && keys %{$preferred};
	my @part = split /\|/, $locus, -1;
	my $catalogueGene = pop @part;
	return $preferred->{$catalogueGene} // 0;
}

sub assignPresortScores {
	# Normalise the guide's ranks over the loci actually present in this stage,
	# never over the guide's own length: a guide pooled across every MGS would
	# otherwise map this MGS's loci into a narrow band of near-identical scores
	# and the term would carry no information. Renormalising per stage also lets
	# the final stage re-spread the survivors across the full range.
	my ($metrics, $preferred) = @_;
	my @ranked;
	for my $gene (keys %{$metrics}) {
		my $rank = presortRankFromLocus($gene, $preferred);
		$metrics->{$gene}{presort_rank} = $rank;
		$metrics->{$gene}{presort_score} = 0;
		push @ranked, $gene if $rank > 0;
	}
	return 0 unless @ranked;
	@ranked = sort {
		$metrics->{$a}{presort_rank} <=> $metrics->{$b}{presort_rank} || $a cmp $b
	} @ranked;
	my $span = $#ranked;
	for my $position (0 .. $#ranked) {
		$metrics->{$ranked[$position]}{presort_score} =
			$span > 0 ? 1 - $position / $span : 1;
	}
	return scalar(@ranked);
}

sub catalogueGeneFromLocus {
	# Membership only; presortRankFromLocus carries the guide's ordering.
	my ($locus, $preferred) = @_;
	return presortRankFromLocus($locus, $preferred) ? 1 : 0;
}

sub compactTaxonAwareDiagnostics {
	return 0 unless $taxonAwareLocusSelection && $compactTaxonAwareDiagnostics;
	my @names = qw(
		taxon_aware_locus_candidates.tsv
		taxon_aware_sample_candidates.tsv
		taxon_aware_locus_selection.tsv
		taxon_aware_sample_selection.tsv
		taxon_aware_backbone_eligibility.tsv
		taxon_aware_placement_eligibility.tsv
		rate_merged_partitions.tsv
	);
	my @sources = map { File::Spec->catfile($treeD, $_) }
		grep { -s File::Spec->catfile($treeD, $_) } @names;
	return 0 unless @sources;
	my $output = File::Spec->catfile($treeD, 'taxon_aware_diagnostics.tsv');
	my $temporary = "$output.write.$$";
	open my $merged, '>', $temporary
		or die "Cannot create consolidated taxon-aware diagnostics $temporary: $!\n";
	print {$merged} "# MATAFILER taxon-aware diagnostics v1\n"
		or die "Cannot write consolidated taxon-aware diagnostics header: $!\n";
	for my $source (@sources) {
		open my $input, '<', $source
			or die "Cannot read taxon-aware diagnostic $source: $!\n";
		print {$merged} "## ".basename($source)."\n"
			or die "Cannot write diagnostics section header: $!\n";
		while (my $line = <$input>) {
			print {$merged} $line
				or die "Cannot write consolidated diagnostic from $source: $!\n";
		}
		close $input or die "Cannot close taxon-aware diagnostic $source: $!\n";
		print {$merged} "\n" or die "Cannot separate consolidated diagnostics: $!\n";
	}
	close $merged or die "Cannot close consolidated taxon-aware diagnostics $temporary: $!\n";
	retry_rename($temporary, $output,
		label => "publish consolidated taxon-aware diagnostics $output");
	for my $source (@sources) {
		retry_unlink($source, fatal => 0,
			label => "remove consolidated taxon-aware diagnostic $source");
	}
	print "Consolidated ".scalar(@sources)
		." taxon-aware/rate diagnostic file(s) into $output\n";
	return scalar(@sources);
}

sub legacyPolicyFileMatches {
	my ($policyFile, $policyText, $description) = @_;
	return 0 unless -s $policyFile;
	my $policyRead = retry_open(q{<}, $policyFile,
		label => "read ".($description || "legacy workflow policy"));
	my $existingPolicy = do { local $/; <$policyRead> };
	retry_close($policyRead, "close ".($description || "legacy workflow policy"));
	return $existingPolicy eq $policyText;
}

sub alignmentFileStem {
	my ($path) = @_;
	my $stem = basename($path);
	$stem =~ s/\.gz$//;
	$stem =~ s/\.(?:syn|nonsyn)\.fna$//;
	$stem =~ s/\.(?:fna|faa)$//;
	return $stem;
}

sub runPostAlignmentLocusQC {
	my ($alignments, $sequenceType, $reportFile) = @_;
	die "Post-alignment locus QC requires at least one alignment\n"
		unless @{$alignments};
	make_path(dirname($reportFile)) unless -d dirname($reportFile);
	my ($manifestFH, $manifestFile) = tempfile(
		"post-alignment-loci-XXXXXX",
		DIR => $tmpD,
		UNLINK => 1,
	);
	print {$manifestFH} "$_\n" for @{$alignments};
	close $manifestFH or die "Cannot close locus-QC manifest $manifestFile: $!\n";
	my ($keepFH, $keepFile) = tempfile(
		"post-alignment-keep-XXXXXX",
		DIR => $tmpD,
		UNLINK => 1,
	);
	close $keepFH or die "Cannot close locus-QC keep file $keepFile: $!\n";

	my $msaFix = getProgPaths("MSAfix");
	my @divergenceArguments = $postAlignmentDivergenceQC
		? ("-relativeModifiedZ", $postAlignmentRelativeZ)
		: (
			"-maxMedianDivergence", 1,
			"-maxP90Divergence", 1,
			"-relativeModifiedZ", 1_000_001,
		);
	my $command = join(" ",
		shellQuote($msaFix),
		"-manifest", shellQuote($manifestFile),
		"-report", shellQuote($reportFile),
		"-keep", shellQuote($keepFile),
		"-threads", $ncore,
		"-sequenceType", $sequenceType,
		"-minSequences", $postAlignmentMinSequences,
		"-minOccupancy", $postAlignmentMinOccupancy,
		"-minOverlapMSA", $minOverlapMSA,
		@divergenceArguments,
		"-minLociForRelative", $postAlignmentMinLociRelative,
	) . "\n";
	my @kept;
	my $ok = eval {
		systemW($command);
		die "Native MSAfix locus QC did not produce its report: $reportFile\n"
			unless -s $reportFile;

		open my $keepRead, '<', $keepFile
			or die "Cannot open locus-QC keep file $keepFile: $!\n";
		while (my $line = <$keepRead>) {
			$line =~ s/[\r\n]+$//;
			push @kept, $line if length($line);
		}
		close $keepRead or die "Cannot close locus-QC keep file $keepFile: $!\n";
		1;
	};
	my $error = $@;
	my @temporaryFiles = (
		$manifestFile,
		$keepFile,
		bsd_glob(quotemeta($reportFile).".tmp.*"),
		bsd_glob(quotemeta($keepFile).".tmp.*"),
	);
	my (%seenTemporary, @cleanupErrors);
	for my $temporaryFile (@temporaryFiles) {
		next unless defined($temporaryFile) && -e $temporaryFile;
		next if $seenTemporary{$temporaryFile}++;
		push @cleanupErrors, "$temporaryFile: $!"
			unless unlink $temporaryFile;
	}
	if (@cleanupErrors) {
		$error .= "Native MSAfix locus-QC temporary cleanup failed:\n"
			.join("\n", @cleanupErrors)."\n";
	}
	if (!$ok || @cleanupErrors) {
		$error ||= "Native MSAfix locus QC failed\n";
		die $error;
	}
	die "Post-alignment locus QC rejected all ".scalar(@{$alignments})
		." loci; see $reportFile\n" unless @kept;
	print "Post-alignment locus QC retained ".scalar(@kept)."/"
		.scalar(@{$alignments})." loci; report: $reportFile\n";
	return \@kept;
}



sub inputFingerprint {
	my ($path) = @_;
	return '' unless defined($path) && length($path);
	my $actual = -e $path ? $path : -e "$path.gz" ? "$path.gz" : $path;
	my @stat = stat($actual);
	return "missing:$actual" unless @stat;
	my $absolute = File::Spec->canonpath(File::Spec->rel2abs($actual));
	return join(':', $absolute, $stat[7], $stat[9]);
}


sub preflightBuildTree {
	my ($outputDirectory, $temporaryDirectory) = @_;
	preflight_directory($outputDirectory, 'BuildTree output directory');
	preflight_directory($temporaryDirectory, 'BuildTree temporary directory');
	for my $entry (['output', $outputDirectory], ['temporary', $temporaryDirectory]) {
		my $capacity = filesystem_capacity($entry->[1]);
		warn "Preflight warning: $entry->[0] filesystem has less than 2 GiB available\n"
			if defined($capacity->{available_kb})
				&& $capacity->{available_kb} < 2 * 1024 * 1024;
		warn "Preflight warning: $entry->[0] filesystem has fewer than 10,000 inodes available\n"
			if defined($capacity->{available_inodes})
				&& $capacity->{available_inodes} < 10_000;
	}
	print "Preflight complete: writable paths, disk space, and inodes checked\n";
}
sub prepareTemporaryBase {
	my ($path) = @_;
	my $created = eval {
		make_path($path) unless -d $path;
		1;
	};
	if (!$created || !-d $path) {
		my $error = $@ || "path is not a directory";
		$error =~ s/\s+$//;
		return (0, $error);
	}

	my ($probeHandle, $probePath);
	my $writeable = eval {
		($probeHandle, $probePath) = tempfile(
			".buildTree5-writecheck-XXXXXX",
			DIR => $path,
			UNLINK => 0,
		);
		print {$probeHandle} "buildTree5 temporary-path check\n"
			or die "write failed: $!";
		close $probeHandle or die "close failed: $!";
		undef $probeHandle;
		1;
	};
	if (!$writeable) {
		my $error = $@ || "write test failed";
		$error =~ s/\s+$//;
		close $probeHandle if defined($probeHandle) && fileno($probeHandle);
		unlink $probePath if defined($probePath) && -e $probePath;
		return (0, $error);
	}
	unless (unlink $probePath) {
		return (0, "cannot remove temporary-path write test $probePath: $!");
	}
	return (1, "");
}

sub elapsedTimeText {
	my ($seconds) = @_;
	$seconds = 0 unless defined($seconds) && $seconds >= 0;
	return sprintf('%.1fs', $seconds) if $seconds < 60;
	my $minutes = int($seconds / 60);
	my $remaining = $seconds - $minutes * 60;
	return sprintf('%dm%.1fs', $minutes, $remaining) if $minutes < 60;
	my $hours = int($minutes / 60);
	$minutes %= 60;
	return sprintf('%dh%dm%.1fs', $hours, $minutes, $remaining);
}

sub postAlignmentStep {
	my ($name, $started, @details) = @_;
	writeWorkflowHeartbeat($name);
	my $elapsed = elapsedTimeText(time - $started);
	my @clean = grep { defined($_) && length($_) } @details;
	print 'POST-ALIGNMENT STEP: '.$name.' ('.$elapsed.')'
		.(@clean ? '; '.join(', ', @clean) : '')."\n";
}

sub readBuildTreeState {
	my ($stateFile) = @_;
	my %state;
	return \%state unless defined($stateFile) && -s $stateFile;
	my $stateHandle = retry_open('<', $stateFile,
		label => 'read BuildTree state');
	while (my $line = <$stateHandle>) {
		$line =~ s/[\r\n]+\z//;
		my ($key, $value) = split /\t/, $line, 2;
		next unless defined($key) && length($key) && defined($value);
		$state{$key} = $value;
	}
	retry_close($stateHandle, 'close BuildTree state');
	return \%state;
}

sub buildTreeStatePolicyMatches {
	my ($state, $key, $policyText) = @_;
	return 0 unless ref($state) eq 'HASH' && exists($state->{$key});
	$policyText //= '';
	$policyText =~ s/[\r\n]+\z//;
	return $state->{$key} eq $policyText;
}

sub writeBuildTreeState {
	return 0 unless length($workflowStateFile) && length($outD) && -d $outD;
	my $reason = $workflowReason // '';
	$reason =~ s/[\t\r\n]+/ /g;
	my $msaSelectionPolicy = $workflowMsaSelectionPolicy // '';
	my $treeStagePolicy = $workflowTreeStagePolicy // '';
	$msaSelectionPolicy =~ s/[\r\n]+\z//;
	$treeStagePolicy =~ s/[\r\n]+\z//;
	my $previousState = readBuildTreeState($workflowStateFile);
	$msaSelectionPolicy = $previousState->{msa_selection_policy}
		if !length($msaSelectionPolicy)
			&& exists($previousState->{msa_selection_policy});
	$treeStagePolicy = $previousState->{tree_stage_policy}
		if !length($treeStagePolicy)
			&& exists($previousState->{tree_stage_policy});
	my ($stateHandle, $temporaryState) = tempfile(
		'buildTree-state-XXXXXX', DIR => $outD, UNLINK => 1,
	);
	print {$stateHandle} join("\n",
		join("\t", schema => 1),
		join("\t", version => $version),
		join("\t", status => ($workflowStatus || 'running')),
		join("\t", stage => ($workflowStage || 'initialization')),
		join("\t", timestamp => time),
		join("\t", pid => $$),
		join("\t", reason => $reason),
		join("\t", msa_selection_policy => $msaSelectionPolicy),
		join("\t", tree_stage_policy => $treeStagePolicy),
	)."\n" or die "Cannot write BuildTree state $temporaryState: $!\n";
	retry_close($stateHandle, "close BuildTree state $temporaryState");
	retry_rename($temporaryState, $workflowStateFile,
		label => "publish BuildTree state $workflowStateFile");
	return 1;
}

sub cleanupLegacyBuildTreeStateFiles {
	return unless length($workflowStateFile) && -s $workflowStateFile;
	my @legacyFiles = (
		File::Spec->catfile($outD, 'buildTree.heartbeat.tsv'),
		File::Spec->catfile($outD, 'buildTree.failure.tsv'),
		File::Spec->catfile($MsaD, 'alignment_work.policy.tsv'),
		File::Spec->catfile($treeD, 'post_alignment.policy.tsv'),
		File::Spec->catfile($treeD, 'post_alignment_locus_qc.policy.tsv'),
	);
	for my $legacyFile (@legacyFiles) {
		next unless -e $legacyFile || -l $legacyFile;
		retry_unlink($legacyFile, fatal => 0,
			label => "remove legacy BuildTree state $legacyFile");
	}
}

sub writeWorkflowHeartbeat {
	my ($stage) = @_;
	$workflowStage = $stage if defined($stage) && length($stage);
	$workflowStatus = $workflowStage eq 'complete' ? 'complete'
		: $workflowStage eq 'placement_pending' ? 'placement_pending' : 'running';
	$workflowReason = '' unless $workflowStatus eq 'placement_pending';
	eval { writeBuildTreeState(); 1; };
}

sub writeWorkflowFailure {
	my ($error) = @_;
	$error //= 'unknown failure';
	$workflowStatus = 'failed';
	$workflowReason = $error;
	eval { writeBuildTreeState(); 1; };
}
sub alignmentCollectionStatsFromReport {
	my ($reportFile) = @_;
	open my $report, '<', $reportFile
		or die "Cannot read alignment-statistics report $reportFile: $!\n";
	my $header = <$report> // '';
	$header =~ s/[\r\n]+\z//;
	my @columns = split /\t/, $header, -1;
	my %columnIndex = map { $columns[$_] => $_ } 0 .. $#columns;
	for my $required (qw(status sequences alignment_sites)) {
		die "Alignment-statistics report $reportFile lacks column '$required'\n"
			unless exists $columnIndex{$required};
	}
	my ($loci, $totalSequences, $totalSites) = (0, 0, 0);
	my ($minimumSequences, $maximumSequences, $minimumLength, $maximumLength);
	while (my $line = <$report>) {
		$line =~ s/[\r\n]+\z//;
		next unless length($line);
		my @fields = split /\t/, $line, -1;
		my $status = $fields[$columnIndex{status}] // '';
		next unless $status eq 'PASS';
		my $sequenceCount = $fields[$columnIndex{sequences}] // '';
		my $alignmentLength = $fields[$columnIndex{alignment_sites}] // '';
		next unless $sequenceCount =~ /^\d+\z/ && $alignmentLength =~ /^\d+\z/;
		$loci++;
		$totalSequences += $sequenceCount;
		$totalSites += $alignmentLength;
		$minimumSequences = $sequenceCount
			if !defined($minimumSequences) || $sequenceCount < $minimumSequences;
		$maximumSequences = $sequenceCount
			if !defined($maximumSequences) || $sequenceCount > $maximumSequences;
		$minimumLength = $alignmentLength
			if !defined($minimumLength) || $alignmentLength < $minimumLength;
		$maximumLength = $alignmentLength
			if !defined($maximumLength) || $alignmentLength > $maximumLength;
	}
	close $report or die "Cannot close alignment-statistics report $reportFile: $!\n";
	return {
		loci => $loci,
		mean_sequences => $loci ? sprintf('%.1f', $totalSequences / $loci) : 0,
		mean_length => $loci ? sprintf('%.1f', $totalSites / $loci) : 0,
		total_sites => $totalSites,
		minimum_sequences => $minimumSequences // 0,
		maximum_sequences => $maximumSequences // 0,
		minimum_length => $minimumLength // 0,
		maximum_length => $maximumLength // 0,
	};
}

sub alignmentCollectionStats {
	my ($alignments) = @_;
	die "Alignment collection must be an array reference\n"
		unless ref($alignments) eq 'ARRAY';
	my ($loci, $totalSequences, $totalSites) = (0, 0, 0);
	my ($minimumSequences, $maximumSequences, $minimumLength, $maximumLength);
	for my $alignment (@{$alignments}) {
		next unless defined($alignment) && fileGZs($alignment);
		my ($input, $ok) = gzipopen($alignment, 'post-alignment statistics', 1);
		die "Cannot read alignment $alignment for post-alignment statistics\n"
			unless $ok && $input;
		my ($sequenceCount, $alignmentLength, $currentLength) = (0, undef, 0);
		while (my $line = <$input>) {
			if ($line =~ /^>/) {
				if ($sequenceCount) {
					die "Unequal sequence lengths in alignment $alignment\n"
						if defined($alignmentLength) && $currentLength != $alignmentLength;
					$alignmentLength = $currentLength unless defined $alignmentLength;
				}
				$sequenceCount++;
				$currentLength = 0;
				next;
			}
			$line =~ s/\s+//g;
			$currentLength += length($line);
		}
		if ($sequenceCount) {
			die "Unequal sequence lengths in alignment $alignment\n"
				if defined($alignmentLength) && $currentLength != $alignmentLength;
			$alignmentLength = $currentLength unless defined $alignmentLength;
		}
		close $input or die "Cannot close alignment $alignment: $!\n";
		next unless $sequenceCount && defined $alignmentLength;
		$loci++;
		$totalSequences += $sequenceCount;
		$totalSites += $alignmentLength;
		$minimumSequences = $sequenceCount if !defined($minimumSequences) || $sequenceCount < $minimumSequences;
		$maximumSequences = $sequenceCount if !defined($maximumSequences) || $sequenceCount > $maximumSequences;
		$minimumLength = $alignmentLength if !defined($minimumLength) || $alignmentLength < $minimumLength;
		$maximumLength = $alignmentLength if !defined($maximumLength) || $alignmentLength > $maximumLength;
	}
	return {
		loci => $loci,
		mean_sequences => $loci ? sprintf('%.1f', $totalSequences / $loci) : 0,
		mean_length => $loci ? sprintf('%.1f', $totalSites / $loci) : 0,
		total_sites => $totalSites,
		minimum_sequences => $minimumSequences // 0,
		maximum_sequences => $maximumSequences // 0,
		minimum_length => $minimumLength // 0,
		maximum_length => $maximumLength // 0,
	};
}

sub summarizeMSAFixLog {
	my ($logFile, $alignment, $failed) = @_;
	return '' unless -s $logFile;
	open my $input, '<', $logFile
		or return "Cannot read captured MSAfix output $logFile: $!\n";
	my ($borderMasked, $lowIdMasked, $minGoodRemoved) = (0, 0, 0);
	my @allLines;
	while (my $line = <$input>) {
		$line =~ s/[\r\n]+\z//;
		next unless length $line;
		push @allLines, $line;
		if ($line =~ /^Border gap check: Masked (\d+) sequences$/) {
			$borderMasked += $1;
		} elsif ($line =~ /^Low ID check: Masked (\d+) sequences$/) {
			$lowIdMasked += $1;
		} elsif ($line =~ /^Removed (\d+) sequences due to less than [0-9.eE+-]+ good positions$/) {
			$minGoodRemoved += $1;
		} else {
			limitedWarn('MSAfix unrecognized diagnostic',
				"MSAfix ".alignmentFileStem($alignment).": $line\n", 5);
		}
	}
	close $input;
	my $total = $borderMasked + $lowIdMasked + $minGoodRemoved;
	if ($total) {
		$msaFixCleanedLoci++;
		$msaFixBorderMasked += $borderMasked;
		$msaFixLowIdMasked += $lowIdMasked;
		$msaFixMinGoodRemoved += $minGoodRemoved;
		limitedWarn('MSAfix per-locus cleaning',
			"MSAfix ".alignmentFileStem($alignment).": border-gap=$borderMasked; low-ID=$lowIdMasked; below-minimum-good-positions=$minGoodRemoved\n", 5);
	}
	return '' unless $failed;
	splice @allLines, 20 if @allLines > 20;
	return @allLines
		? "Captured MSAfix output for $alignment (first ".scalar(@allLines)." line(s)):\n".join("\n", @allLines)."\n"
		: '';
}

sub limitedWarn {
	my ($category, $message, $limit) = @_;
	$limit = 5 unless defined $limit;
	$limitedWarningLimits{$category} = $limit;
	my $count = ++$limitedWarningCounts{$category};
	warn $message if $count <= $limit;
	warn "No more '$category' warning examples will be shown\n"
		if $count == $limit + 1;
}

sub runMSAFix {
	my ($alignment, $maxGapFraction) = @_;
	die "MSAfix input is missing or empty: $alignment\n" unless -s $alignment;

	my $msaFbin = getProgPaths("MSAfix");
	my @codingRecoveryArguments = $msaFixRecoverTechnicalOffsets
		? ("-recoverTechnicalOffsets", "-codingFrame", $msaFixCodingFrame,
			"-geneticCode", $msaFixGeneticCode, "-recoveryBand", $msaFixRecoveryBand)
		: ();
	my $tmpOutput = "$alignment.MSAfix.$$.fna";
	retry_unlink($tmpOutput, label => "clear stale MSAfix output");
	my $tmpLog = "$alignment.MSAfix.$$.log";
	retry_unlink($tmpLog, fatal => 0, label => "clear stale MSAfix log");

	my $cmd = join(" ",
		shellQuote($msaFbin),
		"-i", shellQuote($alignment),
		"-o", shellQuote($tmpOutput),
		"-maskLowID",
		"-threads", $ncore,
		@codingRecoveryArguments,
		"-maskBorderGap",
		"-rmGapColsGreater", $maxGapFraction,
		"-minGoodPosFrac", ($cogCats ne '' ? $ntFracGeneInclude : 0.6),
	);

	my $ok = eval {
		systemW($cmd." > ".shellQuote($tmpLog)." 2>&1");
		1;
	};
	my $outputReady = -s $tmpOutput;
	my $capturedOutput = summarizeMSAFixLog($tmpLog, $alignment, !$ok || !$outputReady);
	retry_unlink($tmpLog, fatal => 0, label => "clean captured MSAfix log");
	if (!$ok) {
		my $error = $@ || "MSAfix failed for $alignment\n";
		retry_unlink($tmpOutput, fatal => 0, label => "clean failed MSAfix output");
		die $error.$capturedOutput;
	}
	if (!$outputReady) {
		retry_unlink($tmpOutput, fatal => 0, label => "clean empty MSAfix output");
		die "MSAfix completed without producing a nonempty output for $alignment\n".$capturedOutput;
	}

	retry_rename($tmpOutput, $alignment, label => "publish validated MSAfix alignment");
}


sub runStagedStrainShardHelper {
	my ($mode, $staging, $manifest, $output, $cores) = @_;
	my $helper = File::Spec->canonpath(File::Spec->catfile(
		$Bin, File::Spec->updir, 'MGS', 'finalize_strain_tree_inputs.pl'));
	die "Staged strain shard helper is missing or empty: $helper\n" unless -s $helper;
	my @command = ($^X, $helper,
		'-staging', $staging, '-manifest', $manifest, '-mode', $mode);
	if ($mode eq 'prepare') {
		push @command, ('-pigz', $pigzBin, '-cores', $cores,
			'-publishedDir', $output);
	} elsif ($mode eq 'cleanup') {
		push @command, ('-publishedDir', $output);
	} else {
		die "Unsupported staged strain shard helper mode '$mode'\n";
	}
	my $status = system {$^X} @command;
	die "Could not execute staged strain shard helper $helper: $!\n" if $status == -1;
	die "Staged strain shard helper failed in $mode mode with status $status\n"
		if $status != 0;
}

sub prepareStagedStrainInputs {
	my ($staging, $requiredInputs, $sampleQCPath, $configuredOutgroup) = @_;
	my $plan = File::Spec->catfile($staging, '.strain_tree_input.plan.tsv');
	return unless -s $plan;
	open my $planIn, '<', $plan or die "Cannot read staged strain input plan $plan: $!\n";
	my $format = <$planIn> // '';
	$format =~ s/[\r\n]+\z//;
	die "Unsupported staged strain input plan $plan\n"
		unless $format eq 'strain-staged-input-v1';
	my $outgroupLine = <$planIn> // '';
	$outgroupLine =~ s/[\r\n]+\z//;
	my ($field, $plannedOutgroup) = split /\t/, $outgroupLine, 2;
	die "Malformed staged strain input plan $plan\n"
		unless defined($field) && $field eq 'outgroup' && defined($plannedOutgroup)
			&& $plannedOutgroup !~ /[\r\n]/;
	my $plannedMGS = '';
	my $mgsLine = <$planIn>;
	if (defined($mgsLine)) {
		$mgsLine =~ s/[\r\n]+\z//;
		if (length($mgsLine)) {
			my ($mgsField, $value) = split /\t/, $mgsLine, 2;
			die "Malformed staged strain MGS plan row in $plan\n"
				unless defined($mgsField) && $mgsField eq 'mgs'
					&& defined($value) && length($value) && $value !~ /[\r\n]/;
			$plannedMGS = $value;
		}
	}
	while (my $extra = <$planIn>) {
		die "Unexpected extra row in staged strain input plan $plan\n"
			if $extra =~ /\S/;
	}
	close $planIn or die "Cannot close staged strain input plan $plan: $!\n";
	my $requestedOutgroup = defined($configuredOutgroup) ? $configuredOutgroup : '';
	die "Staged strain input outgroup '$plannedOutgroup' does not match buildTree outgroup '$requestedOutgroup'\n"
		if $plannedOutgroup ne $requestedOutgroup;

	my $categoryPath = $requiredInputs->[2] // '';
	return unless length($categoryPath);
	my $categoryName = basename($categoryPath);
	my $rawCategory = File::Spec->catfile($staging, "$categoryName.tmp");
	my $finalCategory = File::Spec->catfile($staging, $categoryName);
	my $sampleQCName = length($sampleQCPath) ? basename($sampleQCPath) : '';
	my $rawSampleQC = length($sampleQCName)
		? File::Spec->catfile($staging, "$sampleQCName.tmp") : '';
	my $finalSampleQC = length($sampleQCName)
		? File::Spec->catfile($staging, $sampleQCName) : '';
	my $prepared = File::Spec->catfile($staging, '.strain_tree_input.prepared.tsv');
	if (-s $prepared && fileGZs($finalCategory)
		&& (!length($sampleQCName) || fileGZs($finalSampleQC))) {
		open my $preparedIn, '<', $prepared
			or die "Cannot read staged strain preparation marker $prepared: $!\n";
		my $preparedLine = <$preparedIn> // '';
		close $preparedIn
			or die "Cannot close staged strain preparation marker $prepared: $!\n";
		$preparedLine =~ s/[\r\n]+\z//;
		my @preparedFields = split /\t/, $preparedLine, -1;
		die "Malformed staged strain preparation marker $prepared\n"
			unless @preparedFields >= 4 && $preparedFields[0] eq 'strain-staged-input-v1'
				&& $preparedFields[1] eq $plannedOutgroup;
		die "Staged strain preparation marker $prepared does not match MGS $plannedMGS\n"
			if length($plannedMGS)
				&& (!defined($preparedFields[4]) || $preparedFields[4] ne $plannedMGS);
		return;
	}
	die "Staged strain category input is missing: $rawCategory\n" unless fileGZs($rawCategory);

	my $overlayCategory = File::Spec->catfile($staging, '.strain_tree_input.outgroup.cat.tsv');
	my ($loci, $samples) = finalizeStagedStrainCategory(
		$rawCategory, $overlayCategory, $finalCategory, $plannedMGS);
	if (length($sampleQCName)) {
		die "Staged sample QC input is missing: $rawSampleQC\n" unless fileGZs($rawSampleQC);
		finalizeStagedSampleQC($rawSampleQC, $finalSampleQC, $plannedMGS);
	}
	my @sequenceInputs = @{$requiredInputs}[0, 1];
	my @overlayNames = (
		'.strain_tree_input.outgroup.fna',
		'.strain_tree_input.outgroup.faa',
	);
	for my $index (0 .. $#sequenceInputs) {
		my $inputPath = $sequenceInputs[$index] // '';
		next unless length($inputPath);
		my $source = File::Spec->catfile($staging, basename($inputPath));
		my $overlay = File::Spec->catfile($staging, $overlayNames[$index]);
		stagedOverlayCompletion($source, $overlay,
			File::Spec->catfile($staging, ".strain_tree_input.overlay.$index.attempt"),
			File::Spec->catfile($staging, ".strain_tree_input.overlay.$index.done"));
	}
	my $outgroupLog = File::Spec->catfile($staging, 'data.log');
	writeStagedPreparationMarker($outgroupLog, "OG:$plannedOutgroup\n");
	writeStagedPreparationMarker($prepared, join("\t",
		'strain-staged-input-v1', $plannedOutgroup, $loci, $samples, $plannedMGS)."\n");
	for my $temporary ($rawCategory, $rawSampleQC, $overlayCategory,
		map { File::Spec->catfile($staging, $_) } @overlayNames) {
		next unless defined($temporary) && length($temporary);
		retry_unlink($temporary, fatal => 0,
			label => "clean finalized staged strain input $temporary");
	}
	print "Finalized staged strain input: $loci loci, $samples samples; outgroup="
		.(length($plannedOutgroup) ? $plannedOutgroup : '<none>')."\n";
}

sub finalizeStagedStrainCategory {
	my ($rawCategory, $overlayCategory, $finalCategory, $expectedMGS) = @_;
	my (%loci, %samples);
	my $observedMGS = '';
	my ($input) = gzipopen($rawCategory, 'staged strain category input', 1);
	while (my $line = <$input>) {
		$line =~ s/[\r\n]+\z//;
		next unless length($line);
		my @fields = split /\t/, $line, -1;
		die "Malformed staged strain category row in $rawCategory: $line\n"
			unless @fields >= 4 && length($fields[0]) && length($fields[1])
				&& length($fields[2]) && length($fields[3]);
		$observedMGS = $fields[0] unless length($observedMGS);
		die "Staged strain category input mixes MGS $observedMGS and $fields[0]\n"
			unless $fields[0] eq $observedMGS;
		die "Staged strain category input belongs to $fields[0], expected $expectedMGS\n"
			if defined($expectedMGS) && length($expectedMGS) && $fields[0] ne $expectedMGS;
		$loci{$fields[1]}{$fields[2]} = $fields[3];
		$samples{$fields[2]} = 1;
	}
	close $input or die "Cannot close staged strain category input $rawCategory: $!\n";
	if (fileGZs($overlayCategory)) {
		my ($overlay) = gzipopen($overlayCategory, 'staged outgroup category overlay', 1);
		while (my $line = <$overlay>) {
			$line =~ s/[\r\n]+\z//;
			next unless length($line);
			my @fields = split /\t/, $line, -1;
			die "Malformed staged outgroup category overlay row in $overlayCategory: $line\n"
				unless @fields == 3 && !grep { !length($_) } @fields;
			$loci{$fields[0]}{$fields[1]} = $fields[2];
			$samples{$fields[1]} = 1;
		}
		close $overlay or die "Cannot close staged outgroup category overlay $overlayCategory: $!\n";
	}
	my $temporary = "$finalCategory.write.$$";
	open my $output, '>', $temporary
		or die "Cannot create finalized staged category $temporary: $!\n";
	for my $locus (sort keys %loci) {
		print {$output} join("\t", map { $loci{$locus}{$_} } sort keys %{$loci{$locus}}), "\n"
			or die "Cannot write finalized staged category $temporary: $!\n";
	}
	close $output or die "Cannot close finalized staged category $temporary: $!\n";
	retry_rename($temporary, $finalCategory,
		label => "publish finalized staged category $finalCategory");
	return (scalar(keys %loci), scalar(keys %samples));
}

sub finalizeStagedSampleQC {
	my ($rawSampleQC, $finalSampleQC, $expectedMGS) = @_;
	my %sample;
	my $observedMGS = '';
	my ($input) = gzipopen($rawSampleQC, 'staged strain sample QC', 1);
	while (my $line = <$input>) {
		$line =~ s/[\r\n]+\z//;
		next if $line eq '' || $line =~ /^#/ || $line =~ /^MGS\t/;
		my @fields = split /\t/, $line, -1;
		die "Malformed staged sample QC row in $rawSampleQC: $line\n"
			unless @fields >= 6 && length($fields[0]) && length($fields[1]);
		my ($mgs, $sampleId, $status, $ambiguous, $csp, $loci) = @fields[0 .. 5];
		$observedMGS = $mgs unless length($observedMGS);
		die "Staged sample QC input mixes MGS $observedMGS and $mgs\n"
			unless $mgs eq $observedMGS;
		die "Staged sample QC input belongs to $mgs, expected $expectedMGS\n"
			if defined($expectedMGS) && length($expectedMGS) && $mgs ne $expectedMGS;
		my $entry = $sample{$sampleId};
		if (!$entry) {
			$sample{$sampleId} = [$mgs, $status, 0 + $ambiguous, 0 + $csp, 0 + $loci];
		} else {
			die "Staged sample QC sample '$sampleId' belongs to both $entry->[0] and $mgs\n"
				unless $entry->[0] eq $mgs;
			# Raw staged rows may still carry the pre-1.08 spelling.
			$entry->[1] = $status
				if canonical_sample_qc_status($status) eq 'mixed_strain';
			$entry->[2] = $ambiguous if $ambiguous > $entry->[2];
			$entry->[3] = $csp if $csp > $entry->[3];
			$entry->[4] = $loci if $loci > $entry->[4];
		}
	}
	close $input or die "Cannot close staged sample QC $rawSampleQC: $!\n";
	my $temporary = "$finalSampleQC.write.$$";
	open my $output, '>', $temporary
		or die "Cannot create finalized staged sample QC $temporary: $!\n";
	print {$output} join("\t", qw(MGS sample status ambiguous_fraction csp_fraction validated_loci)), "\n";
	for my $sampleId (sort keys %sample) {
		print {$output} join("\t", $sample{$sampleId}[0], $sampleId, @{$sample{$sampleId}}[1 .. 4]), "\n"
			or die "Cannot write finalized staged sample QC $temporary: $!\n";
	}
	close $output or die "Cannot close finalized staged sample QC $temporary: $!\n";
	retry_rename($temporary, $finalSampleQC,
		label => "publish finalized staged sample QC $finalSampleQC");
	return scalar(keys %sample);
}

sub stagedOverlayRecordsPresent {
	my ($source, $overlay) = @_;
	return 1 unless fileGZs($overlay);
	return 0 unless fileGZs($source);
	my %needed;
	my ($overlayIn) = gzipopen($overlay, 'staged outgroup FASTA overlay', 1);
	while (my $line = <$overlayIn>) {
		$needed{$1} = 1 if $line =~ /^>(\S+)/;
	}
	close $overlayIn or die "Cannot close staged outgroup FASTA overlay $overlay: $!\n";
	return 1 unless keys %needed;
	my ($sourceIn) = gzipopen($source, 'staged FASTA overlay recovery scan', 1);
	while (my $line = <$sourceIn>) {
		delete $needed{$1} if $line =~ /^>(\S+)/;
		last unless keys %needed;
	}
	close $sourceIn or die "Cannot close staged FASTA overlay recovery scan $source: $!\n";
	return !keys %needed;
}

sub stagedOverlayCompletion {
	my ($source, $overlay, $attempt, $complete) = @_;
	return unless fileGZs($overlay);
	return if -s $complete;
	die "Staged FASTA input is missing: $source\n" unless fileGZs($source);
	if (-e $attempt) {
		append_fasta_records_atomic($source, do {
			open my $overlayIn, '<', $overlay
				or die "Cannot read staged outgroup FASTA overlay $overlay: $!\n";
			local $/;
			my $records = <$overlayIn> // '';
			close $overlayIn or die "Cannot close staged outgroup FASTA overlay $overlay: $!\n";
			stagedOverlayRecordsPresent($source, $overlay) ? '' : $records;
		});
	} else {
		writeStagedPreparationMarker($attempt, "attempt\n");
		open my $overlayIn, '<', $overlay
			or die "Cannot read staged outgroup FASTA overlay $overlay: $!\n";
		local $/;
		my $records = <$overlayIn> // '';
		close $overlayIn or die "Cannot close staged outgroup FASTA overlay $overlay: $!\n";
		append_fasta_records_atomic($source, $records);
	}
	writeStagedPreparationMarker($complete, "complete\n");
}

sub writeStagedPreparationMarker {
	my ($path, $contents) = @_;
	my $temporary = "$path.write.$$";
	open my $output, '>', $temporary or die "Cannot create staged preparation file $temporary: $!\n";
	print {$output} $contents or die "Cannot write staged preparation file $temporary: $!\n";
	close $output or die "Cannot close staged preparation file $temporary: $!\n";
	retry_rename($temporary, $path, label => "publish staged preparation file $path");
}

sub stagedTreeInputFiles {
	my ($staging) = @_;
	return () unless -d $staging;
	opendir my $directoryHandle, $staging
		or die "Cannot read staged tree-input directory $staging: $!\n";
	my @stagedFiles = sort map { File::Spec->catfile($staging, $_) }
		grep {
			$_ ne File::Spec->curdir && $_ ne File::Spec->updir
				&& $_ !~ /^\./ && $_ !~ /\.(?:tmp|write|rewrite|sort|merge)(?:\.|\z)/
				&& $_ !~ /\.\d+\z/ && $_ ne 'merge.complete.tsv'
				&& -f File::Spec->catfile($staging, $_)
				&& -s File::Spec->catfile($staging, $_)
		} readdir $directoryHandle;
	closedir $directoryHandle
		or die "Cannot close staged tree-input directory $staging: $!\n";
	return @stagedFiles;
}

sub publishStagedTreeInputs {
	my ($stagingDirectory, $outputDirectory, $cores, $requiredInputs,
		$sampleQCPath, $configuredOutgroup) = @_;
	my @missing = grep { !fileGZs($_) } @{$requiredInputs};
	my $staging = File::Spec->canonpath(File::Spec->rel2abs($stagingDirectory));
	my $output = File::Spec->canonpath(File::Spec->rel2abs($outputDirectory));
	my $stagedPlan = File::Spec->catfile($staging, ".strain_tree_input.plan.tsv");
	my $shardManifest = File::Spec->catfile($staging, ".strain_tree_input.shards.tsv");
	my $hasShardHandoff = -s $shardManifest ? 1 : 0;
	my $stagedPrimaryInput = -s $stagedPlan && -d $staging && grep {
		fileGZs(File::Spec->catfile($staging, basename($_)))
	} @{$requiredInputs};
	my @stagedFiles = stagedTreeInputFiles($staging);
	my $stagedResidualInput = -s $stagedPlan && (@stagedFiles || $hasShardHandoff);
	unless (@missing || $stagedPrimaryInput || $stagedResidualInput) {
		print "Using existing persistent tree inputs\n";
		return;
	}

	die "Staged tree-input directory does not exist: $staging\n" unless -d $staging;
	die "Staged tree-input directory must differ from output directory: $staging\n"
		if $staging eq $output;
	if ($hasShardHandoff) {
		runStagedStrainShardHelper(
			'prepare', $staging, $shardManifest, $output, $cores);
		$hasShardHandoff = -s $shardManifest ? 1 : 0;
	}
	prepareStagedStrainInputs($staging, $requiredInputs, $sampleQCPath, $configuredOutgroup);
	@missing = grep { !fileGZs($_) } @{$requiredInputs};

	@stagedFiles = stagedTreeInputFiles($staging);
	if (!@stagedFiles) {
		die "No usable staged tree inputs found in $staging\n" if @missing;
		runStagedStrainShardHelper('cleanup', $staging, $shardManifest, $output, $cores)
			if $hasShardHandoff;
		print "Tree inputs ready in persistent storage\n";
		return;
	}
	my %stagedBasename = map { basename($_) => 1 } @stagedFiles;
	my @unavailable = grep {
		my $requiredBasename = basename($_);
		!$stagedBasename{$requiredBasename} && !$stagedBasename{"$requiredBasename.gz"};
	} @missing;
	die "Required tree inputs are absent from both persistent and staged storage: "
		.join(", ", @unavailable)."\n" if @unavailable;

	print "Publishing ".scalar(@stagedFiles)." staged tree-input file(s) to $output\n";
	for my $source (@stagedFiles) {
		if ($source !~ /\.gz$/) {
			my $sourceBasename = basename($source);
			sortFastaForCompression($source)
				if $sourceBasename eq "allFAAs.faa" || $sourceBasename eq "allFNAs.fna";
			my $compressionStatus = system {$pigzBin} $pigzBin, "-p", $cores, "--", $source;
			die "Could not execute pigz for staged tree input $source: $!\n"
				if $compressionStatus == -1;
			die "pigz failed for staged tree input $source with status $compressionStatus\n"
				if $compressionStatus != 0;
			$source .= ".gz";
			die "Compression did not produce staged tree input $source\n" unless -s $source;
		}
		my $destination = File::Spec->catfile($output, basename($source));
		retry_operation(
			label => "publish staged tree input $destination",
			code => sub { move($source, $destination) && -s $destination },
		);
	}

	@missing = grep { !fileGZs($_) } @{$requiredInputs};
	if (@missing) {
		die "Tree inputs remain incomplete after staged publication; missing: "
			.join(", ", map { $_ . "[.gz]" } @missing)."\n";
	}
	runStagedStrainShardHelper('cleanup', $staging, $shardManifest, $output, $cores)
		if $hasShardHandoff;
	print "Tree inputs ready in persistent storage\n";
}

sub reusableCompletionTree {
	my ($markerPath, $outputDirectory) = @_;
	return '' unless $continue && defined($markerPath) && -s $markerPath;
	open my $marker, '<', $markerPath or return '';
	my $line = <$marker> // '';
	close $marker or return '';
	$line =~ s/[\r\n]+\z//;
	my ($producer, $markerVersion, $treePath) = split /\t/, $line, 3;
	return '' unless defined($producer) && $producer eq 'buildTree5'
		&& defined($markerVersion)
		&& $markerVersion =~ /^\d+(?:\.\d+)?\z/ && $markerVersion >= 5.40
		&& defined($treePath) && length($treePath);

	my $output = File::Spec->canonpath(File::Spec->rel2abs($outputDirectory));
	my $tree = File::Spec->canonpath(File::Spec->rel2abs($treePath, $output));
	my $relative = File::Spec->abs2rel($tree, $output);
	return '' if $relative eq File::Spec->curdir
		|| $relative =~ /^\.\.(?:[\\\/]|$)/;
	return -s $tree ? $tree : '';
}

sub writeCompletionMarker {
	my ($markerPath, $treePath, $outputDirectory) = @_;
	die "Cannot create completion marker without a nonempty primary tree: $treePath\n"
		unless defined($treePath) && length($treePath) && -s $treePath;

	my $marker = File::Spec->canonpath(File::Spec->rel2abs($markerPath));
	my $output = File::Spec->canonpath(File::Spec->rel2abs($outputDirectory));
	my $relative = File::Spec->abs2rel($marker, $output);
	my @relativeParts = File::Spec->splitdir($relative);
	die "Completion marker must be inside the output directory: $marker\n"
		if $relative eq File::Spec->curdir || grep { $_ eq File::Spec->updir } @relativeParts;

	my $temporaryMarker = "$marker.tmp.$$";
	retry_unlink($temporaryMarker, label => "clear stale completion marker");
	my $markerHandle = retry_open('>', $temporaryMarker, label => "create completion marker");
	print {$markerHandle} "buildTree5\t$version\t$treePath\n"
		or die "Cannot write completion marker $temporaryMarker: $!\n";
	retry_close($markerHandle, "close completion marker $temporaryMarker");
	retry_rename($temporaryMarker, $marker, label => "publish completion marker $marker");
	print "Validated primary tree and published completion marker: $marker\n";
}



sub writeOutcomeMarker {
	my ($markerPath, $status, $reason, $details, $outputDirectory) = @_;
	return unless defined($markerPath) && length($markerPath);
	die "Unsafe lifecycle status '$status'\n" unless $status =~ /^[a-z0-9_]+$/;
	my $marker = File::Spec->canonpath(File::Spec->rel2abs($markerPath));
	my $output = File::Spec->canonpath(File::Spec->rel2abs($outputDirectory));
	my $relative = File::Spec->abs2rel($marker, $output);
	die "Lifecycle marker must be inside the output directory: $marker\n"
		if $relative eq File::Spec->curdir || $relative =~ /^\.\.(?:[\\\/]|$)/;
	my $temporary = "$marker.tmp.$$";
	retry_unlink($temporary, label => "clear stale lifecycle marker");
	my $handle = retry_open('>', $temporary, label => "create lifecycle marker");
	$reason //= ''; $reason =~ s/[\t\r\n]+/ /g;
	print {$handle} "status\t$status\nreason\t$reason\nversion\t$version\n"
		or die "Cannot write lifecycle marker $temporary: $!\n";
	for my $key (sort keys %{$details || {}}) {
		my $value = defined($details->{$key}) ? $details->{$key} : '';
		$value =~ s/[\t\r\n]+/ /g;
		print {$handle} "$key\t$value\n"
			or die "Cannot write lifecycle marker $temporary: $!\n";
	}
	retry_close($handle, "close lifecycle marker $temporary");
	retry_rename($temporary, $marker, label => "publish lifecycle marker $marker");
	print "Published lifecycle outcome: $status ($marker)\n";
}

sub completeTaxonAwareOutgroupAnchorTerminal {
	my ($stage, $error) = @_;
	my $reason = 'taxon_aware_outgroup_no_selected_anchor';
	$error //= '';
	$error =~ s/[\r\n]+\z//;
	clearLifecycleMarker($completionMarker, 'clear stale tree completion');
	clearLifecycleMarker($placementPendingMarker,
		'clear stale placement-pending marker');
	cleanupLegacyBuildTreeStateFiles();
	$selectionAttrition{final_loci} = 0;
	$selectionAttrition{final_samples} = 0;
	$selectionAttrition{backbone_samples} = 0;
	$selectionAttrition{placement_samples} = 0;
	$selectionAttrition{excluded_samples} = scalar(keys %samples);
	writeSelectionAttritionAudit($selectionAttritionReport, \%selectionAttrition);
	writeOutcomeMarker($terminalMarker, 'valid_no_tree', $reason, {
		stage => $stage,
		outgroup => $outgroup,
		candidate_loci => $selectionAttrition{candidate_loci} // 0,
		aligned_loci => $selectionAttrition{aligned_loci} // 0,
		post_qc_loci => $selectionAttrition{post_qc_loci} // 0,
		error => $error,
	}, $outD);
	finalizeMSAArtifacts($MsaD, $MsaWorkD);
	safeRemoveTree($tmpD, $tmpBase);
	compactTaxonAwareDiagnostics();
	writeWorkflowHeartbeat('complete');
	print "BuildTree completed with a valid terminal no-tree outcome: $reason\n";
	exit(0);
}

sub clearLifecycleMarker {
	my ($marker, $label) = @_;
	return 1 unless defined($marker) && length($marker);
	return retry_unlink($marker, fatal => 0, label => $label || "clear lifecycle marker");
}
sub safeRemoveTree{
	my ($path, $parent) = @_;
	return unless defined($path) && ($path ne "") && (-d $path || -l $path);
	my $absolutePath = File::Spec->canonpath(File::Spec->rel2abs($path));
	my $absoluteParent = File::Spec->canonpath(File::Spec->rel2abs($parent));
	my $relative = File::Spec->abs2rel($absolutePath, $absoluteParent);
	die "Refusing to remove $absolutePath outside $absoluteParent\n"
		if $relative eq File::Spec->curdir || $relative =~ /^\.\.(?:[\\\/]|$)/;
	retry_operation(
		label => "remove directory tree $absolutePath",
		code => sub {
			my $errors;
			remove_tree($absolutePath, {error => \$errors});
			if ($errors && @{$errors}) {
				my @messages = map { my ($p, $m) = %{$_}; "$p: $m" } @{$errors};
				die join('; ', @messages)."\n";
			}
			return !-e $absolutePath && !-l $absolutePath;
		},
	);
}
sub treeAlignmentCheckpointStatus {
	my ($path, $isAA) = @_;
	return (0, 'file is missing or empty') unless fileGZe($path);
	my ($input, $opened) = gzipopen($path, 'tree alignment checkpoint', 1);
	return (0, 'file cannot be read') unless $opened && $input;
	my (%seen, $identifier, $sequenceCount, $alignmentLength);
	my $sequence = '';
	my $finishRecord = sub {
		return (1, '') unless defined $identifier;
		return (0, "empty sequence for '$identifier'") unless length($sequence);
		return (0, "unequal sequence length for '$identifier'")
			if defined($alignmentLength) && length($sequence) != $alignmentLength;
		$alignmentLength = length($sequence) unless defined $alignmentLength;
		my $alphabet = $isAA
			? qr/^[ACDEFGHIKLMNPQRSTVWYBXZJUO?*-]+$/i
			: qr/^[ACGTURYSWKMBDHVN?-]+$/i;
		return (0, "unsupported sequence character for '$identifier'")
			unless $sequence =~ $alphabet;
		my $informative = $sequence;
		$informative =~ s/[^ACDEFGHIKLMNPQRSTVWY]//gi if $isAA;
		$informative =~ s/[^ACGTU]//gi unless $isAA;
		return (0, "no informative residues for '$identifier'") unless length($informative);
		$sequenceCount++;
		return (1, '');
	};
	while (my $line = <$input>) {
		$line =~ s/[\r\n]+\z//;
		next unless length($line);
		if ($line =~ /^>(\S+)/) {
			my ($ok, $reason) = $finishRecord->();
			if (!$ok) {
				close $input;
				return (0, $reason);
			}
			$identifier = $1;
			if ($seen{$identifier}++) {
				close $input;
				return (0, "duplicate FASTA identifier '$identifier'");
			}
			$sequence = '';
			next;
		}
		if (!defined $identifier) {
			close $input;
			return (0, 'sequence occurs before its first FASTA header');
		}
		$line =~ s/\s+//g;
		$sequence .= $line;
	}
	my ($ok, $reason) = $finishRecord->();
	my $closed = close $input;
	return (0, 'compressed alignment cannot be read completely') unless $closed;
	return (0, $reason) unless $ok;
	return (0, 'no FASTA sequences') unless $sequenceCount;
	return (0, 'fewer than two FASTA sequences') if $sequenceCount < 2;
	return (0, 'zero alignment length') unless $alignmentLength;
	return (1, '', {
		sequences => $sequenceCount,
		alignment_length => $alignmentLength,
	});
}
sub copyStreamToPathAtomically {
	my ($input, $destination, $description) = @_;
	$description ||= "publish $destination";
	make_path(dirname($destination)) unless -d dirname($destination);
	my $temporary = "$destination.write.$$";
	retry_unlink($temporary, label => "clear temporary $temporary");
	my $output = retry_open('>', $temporary, label => "open temporary $temporary");
	my $ok = eval {
		my $buffer;
		while (1) {
			my $read = read($input, $buffer, 1024 * 1024);
			die "Cannot read stream for $description: $!\n" unless defined $read;
			last unless $read;
			print {$output} $buffer
				or die "Cannot write temporary $temporary: $!\n";
		}
		retry_close($output, "close temporary $temporary");
		die "Input stream failed while attempting to $description\n"
			unless close($input);
		die "Temporary artifact is empty after attempting to $description: $temporary\n"
			unless -s $temporary;
		retry_rename($temporary, $destination, label => $description);
		1;
	};
	if (!$ok) {
		my $error = $@ || "Unknown stream-copy failure while attempting to $description\n";
		close($output);
		close($input);
		retry_unlink($temporary, fatal => 0,
			label => "remove failed temporary $temporary");
		die $error;
	}
	return $destination;
}


sub publishCompressedMSAArtifact {
	my ($working, $published) = @_;
	die "Cannot publish missing or empty scratch MSA $working\n" unless -s $working;
	make_path(dirname($published)) unless -d dirname($published);
	my $compressed = "$published.gz";
	open my $compressor, '-|', $pigzBin, '-p', $ncore, '-c', '--', $working
		or die "Cannot start MSA compression for $working: $!\n";
	copyStreamToPathAtomically(
		$compressor, $compressed, "publish compressed MSA checkpoint $compressed");
	retry_unlink($published, label => "remove persistent plain MSA $published")
		if -e $published || -l $published;
	return $compressed;
}


sub restoreCompressedMSAArtifact {
	my ($published, $working) = @_;
	die "restoreCompressedMSAArtifact requires persistent and scratch paths\n"
		unless defined($published) && length($published)
			&& defined($working) && length($working);
	return $working if -s $working;
	my $source = -s "$published.gz" ? "$published.gz" : $published;
	my ($input, $opened, $resolved) =
		gzipopen($source, 'retained MSA artifact', 0, 0);
	return $working unless $opened && $input;
	make_path(dirname($working)) unless -d dirname($working);
	copyStreamToPathAtomically(
		$input, $working, "publish restored scratch MSA $working");
	if (defined($resolved) && File::Spec->canonpath($resolved)
			eq File::Spec->canonpath($published)) {
		# Migrate legacy uncompressed checkpoints only after the scratch copy is safe.
		publishCompressedMSAArtifact($working, $published);
	}
	retry_unlink($published, label => "remove stale persistent plain MSA $published")
		if -e $published || -l $published;
	return $working;
}

sub publishMSAArtifactSet {
	my ($working, $published) = @_;
	return unless -s $working;
	publishCompressedMSAArtifact($working, $published);
	publishCompressedMSAArtifact("$working.nxs", "$published.nxs")
		if -s "$working.nxs";
	publishCompressedMSAArtifact($working.$partiExt, $published.$partiExt)
		if -s $working.$partiExt;
	return "$published.gz";
}

sub restoreMSAArtifactSet {
	my ($published, $working) = @_;
	restoreCompressedMSAArtifact($published, $working);
	restoreCompressedMSAArtifact($published.$partiExt, $working.$partiExt);
	return $working;
}

sub finalizeMSAArtifacts {
	my ($directory, $workingDirectory) = @_;
	return unless defined($directory) && -d $directory;
	my $publishedCount = 0;
	if (defined($workingDirectory) && -d $workingDirectory) {
		opendir my $workingHandle, $workingDirectory
			or die "Cannot inspect scratch MSA directory $workingDirectory: $!\n";
		for my $name (sort readdir $workingHandle) {
			next unless $name =~ /^MSAli.*\.fna\z/;
			my $working = File::Spec->catfile($workingDirectory, $name);
			next unless -s $working;
			my $published = File::Spec->catfile($directory, $name);
			publishMSAArtifactSet($working, $published);
			$publishedCount++;
		}
		closedir $workingHandle
			or die "Cannot close scratch MSA directory $workingDirectory: $!\n";
	}

	opendir my $directoryHandle, $directory
		or die "Cannot inspect MSA directory $directory for finalization: $!\n";
	my (@singleLocusNucleotide, @singleLocusProtein, @retainedPlain);
	my $cleanedSubdirectory = '';
	for my $name (readdir $directoryHandle) {
		next if $name eq File::Spec->curdir || $name eq File::Spec->updir;
		my $path = File::Spec->catfile($directory, $name);
		if (-d $path) {
			$cleanedSubdirectory = $path if $name eq 'clnd';
			next;
		}
		next unless -f $path || -l $path;
		if ($name =~ /^MSAli.*(?:\.fna(?:\.nxs)?|\Q$partiExt\E)(?:\.gz)?\z/) {
			push @retainedPlain, $path unless $name =~ /\.gz\z/;
			next;
		}
		push @singleLocusNucleotide, $path if $name =~ /\.fna(?:\.gz)?\z/;
		push @singleLocusProtein, $path if $name =~ /\.faa(?:\.gz)?\z/;
	}
	closedir $directoryHandle
		or die "Cannot close MSA directory $directory after finalization scan: $!\n";

	for my $path (sort @retainedPlain) {
		if (-l $path || !-s $path) {
			retry_unlink($path, label => "remove completed scratch MSA link $path");
			next;
		}
		publishCompressedMSAArtifact($path, $path);
		$publishedCount++;
	}

	# Per-locus resume keys on the protein alignment in both sequence modes, so
	# -rmMSA 0 has to retain the .faa checkpoints alongside the .fna ones;
	# discarding either half would make the retained half unusable for resume.
	my $removedSingleLocus = 0;
	my $retainedSingleLocus = 0;
	for my $path (@singleLocusNucleotide, @singleLocusProtein) {
		if ($removeMSA) {
			retry_unlink($path, label => "remove completed locus MSA $path");
			$removedSingleLocus++;
			next;
		}
		if ($path =~ /\.gz\z/) {
			$retainedSingleLocus++;
		} elsif (-l $path || !-s $path) {
			retry_unlink($path, label => "remove completed scratch locus-MSA link $path");
		} else {
			publishCompressedMSAArtifact($path, $path);
			$retainedSingleLocus++;
		}
	}
	safeRemoveTree($cleanedSubdirectory, $directory) if $cleanedSubdirectory ne '';
	print "MSA finalization: removed $removedSingleLocus single-locus alignment file(s)"
		.($retainedSingleLocus ? "; retained $retainedSingleLocus compressed single-locus MSA(s)" : '')
		.($cleanedSubdirectory ne '' ? '; removed legacy MSA/clnd' : '')
		."; published $publishedCount compressed MSAli checkpoint(s)\n";
	return $publishedCount;
}

sub msaOnlyArtifacts {
	my ($directory) = @_;
	return [] unless defined($directory) && -d $directory;
	opendir my $handle, $directory
		or die "Cannot inspect MSA-only artifact directory $directory: $!\n";
	my @artifacts = map { File::Spec->catfile($directory, $_) }
		sort grep {
			$_ !~ /^MSAli/ && /\.(?:faa|fna)\.gz\z/
				&& -f File::Spec->catfile($directory, $_)
				&& -s File::Spec->catfile($directory, $_)
		} readdir $handle;
	closedir $handle
		or die "Cannot close MSA-only artifact directory $directory: $!\n";
	return \@artifacts;
}

sub requireConfiguredTool{
	my ($configKey, $description) = @_;
	my $configured = getProgPaths($configKey, 0);
	die "$description support is dormant and not configured. Set $configKey in the selected MATAFILER config to reactivate it.\n"
		if $configured eq "";
	return $configured;
}



### Fastgear -> test for recombination 
sub runFastgear($ $ $ $){
	my ($geneFG, $outFile, $inD, $parFile) = @_;
	my $fastgearBin = requireConfiguredTool("fastgear", "fastGEAR");
	my $matlabBin = requireConfiguredTool("fastgearMatlab", "fastGEAR MATLAB runtime");

	$cmd = "$fastgearBin $matlabBin ".shellQuote("$inD/$geneFG.fna")
		." ".shellQuote($outFile)." $parFile";
	systemW($cmd);
	die "fastGEAR did not produce $outFile for $geneFG\n" unless -s $outFile;
	#die;
}




sub prepGenoDirs($){
	my ($genoindir) = @_;
	# -genoInD '/hpc-home/hildebra/grp/data/DB/Genomes/Ecoli_ref/*.fna'
	if ($genoindir eq ""){return;}
	$aaFna="$outD/autoFAA.faa";$fnFna="$outD/autoFNA.fna";$cogCats="$outD/autoCAT.cat";
	my $SaSe = "_"; my %catT;
	$smplSep = "_"; compileSampleSeparator(); #keep parseSeqId in step with the forced separator
	
	#return(); #DEBUG
	if (-d $genoindir){
		if ($wildcardflag ne ""){
			$genoindir .= $wildcardflag;#"/*\.fna";
		} else {
			warn "No -wildcardflag supplied for -genoInD; defaulting to *.fna\n";
			$genoindir .= "/*\.fna";
		}
	}
	#print "$genoindir\n\n";
	my @sfiles = glob($genoindir);
	if (scalar(@sfiles) == 0){die "$genoindir contains no files!\n";}
	
	#die "@sfiles\n";
	print "External-genome preparation: " . scalar(@sfiles) . " input genome(s) matching $genoindir\n";
	my %FNAfmg; my %FAAfmg;my %MGSFMG;
	
	my $cnt=0;
	foreach my $tarG(@sfiles){
		next if ($tarG  =~ /\*/ || $tarG  =~ /\.genes\./);
	#my $tarG = "$tarDir/$genoN";
		my $tag = $tarG;
		$tag =~ s/.*\///;
		$tag =~ s/\.[^\.]*$//;
		my ($genes,$prots) = getGenoGenes($tarG,0,$ncore);
		my $FMGdir = getFMG("",$prots,$genes,$ncore);
		my ($hrN,$hrA,$hrC) = readFMGdir( "$FMGdir",$tag ,".");
		$cnt++;
		print "External-genome progress: $cnt/" . scalar(@sfiles) . " genomes\n"
			if $cnt == 1 || $cnt % 25 == 0;
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
	print "External-genome summary: $cnt genome(s) processed\n";
	
	
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
