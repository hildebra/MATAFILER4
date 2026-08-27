<!-- Documentation navigation -->
[Home](../README.md) | [Quick start](quickstart.md) | [Workflows](common_workflows.md) | [Strain-within guide](strainwithin.md) | [Outputs](outputs.md) | [Flag reference](flag_reference.md)

---

# Strain-within workflow guide

Quick links: [Before starting](#before-starting) | [Tutorial](#tutorial) | [`strain_within.pl` flags](#strain_withinpl-flag-groups) | [`strain_within_2.2.pl` flags](#strain_within_22pl-flags) | [Outputs to review](#outputs-to-review)

This guide covers the two strain workflow scripts (the wrapper is sometimes called `strainwithin.sh` in run notes; its repository entry point is `strain_within.pl`):

- `secScripts/MGS/strain_within.pl` prepares strain-resolved MGS inputs, filters loci and samples, submits tree jobs, and then launches postprocessing.
- `secScripts/MGS/strain_within_2.2.pl` combines completed within-MGS trees with metadata, runs `strainStats` and optional `PopGenStats`, and produces cohort-level summaries.

Normally `MGS.pl` launches `strain_within.pl` with the catalogue-specific paths already filled in. Run either script directly only when those inputs are known and a restart or targeted postprocessing operation is intended.

## Before starting

You need a completed gene catalogue (`-GCd`), a nonempty MGS guide (`-MGS`), a sample metadata map, and completed per-MGS phylogenies or the inputs from which they can be built. The wrapper derives its default output structure from the MGS guide; use `-outD` when resuming an existing strain-workflow directory.

Start with `-submit 0`. It validates paths and writes planned job commands without submitting them. Re-run the same command with `-submit 1` after reviewing the generated commands and scheduler configuration.

## Tutorial

### 1. Run the tree-building wrapper

The following is an illustrative direct invocation. Replace every path with paths from the current gene-catalogue run.

```bash
perl "$MF4DIR/secScripts/MGS/strain_within.pl" \
  -GCd /project/gene_catalogue \
  -MGS /project/gene_catalogue/Binning/Bin_SB/SB.clusters.core \
  -outD /project/gene_catalogue/within_phylo \
  -map2 /project/samples.map \
  -MGSabundance /project/gene_catalogue/Annotation/Abundance/MGS.matL7.txt \
  -submit 0
```

The wrapper extracts and QC-filters strain loci (Phase I), prepares each MGS tree (Phase II), waits for submitted jobs, and launches `strain_within_2.2.pl` for statistics. Use the same command with `-submit 1` to submit jobs.

By default, the wrapper builds one tree from the complete retained alignment. EPA-ng placement is **disabled**. To infer a well-covered backbone and place eligible sparse strains afterward, explicitly add `-placeOnBackbone 1`.

Split Phase-I workers share one locus model. Same-COG catalogue seeds merge into a single locus using catalogue-wide co-occurrence and synteny context, which no worker can reconstruct from its own sample slice, so the parent builds the model once over the complete cluster index and publishes it to shared scratch before submitting workers. Each worker loads that model and derives only the member contexts of its own samples. The `locus-group construction` step reports `model_source=`: `published_common_model` for the normal path, `catalogue_wide_build` for the parent that built and published it, `catalogue_wide_rebuild` when a worker had to rebuild it because the published copy was missing, and `local_build` for unsplit runs.

On Slurm, split Phase-I workers that fail validation use one shared bounded retry path. When scheduler accounting confirms `OUT_OF_MEMORY`, the affected worker's memory request doubles on its next retry, up to `-treeOOMMaxMemGB`; non-OOM failures retain the previous request. Strain BuildTree jobs request node-local scratch at five times the compressed input size, with a 20 GiB minimum, when node-local scratch is configured. When `-maxCores` is positive, a full BuildTree job requests `ceil(sqrt(submitted samples))` cores, with the existing four-core floor and `-maxCores` cap; the submitted-sample count is reused from input preparation rather than inferred from file size. EPA-only recovery remains single-threaded.

To generate only the per-locus alignments, add `-onlyMSA 1`. Each locus follows the normal localized pipeline—alignment filtering, NT backtranslation, MSAfix, and compressed checkpoint publication—and successful MGS retain the resulting non-merged `MSA/*.fna.gz` files plus `msaOnly.complete.tsv`. BuildTree exits before combined-MSA postprocessing and concatenation, while the wrapper skips phylogeny, EPA-ng placement, and `strain_within_2.2.pl`. Resume an interrupted alignment-only run with the same option plus `-onlySubmit 1`. Because no backbone exists in this mode, it cannot be combined with `-placeOnBackbone 1`.

### 2. Two independent reasons a sample leaves a tree

These are separate mechanisms and it is worth keeping them apart when reading the reports.

**Mixed strain** — extraction QC judges whether a sample describes *one* strain. A sample is flagged when at least `-minBadLociPSmpl` (3) of its loci are bad in one of two ways, and those exceed a fraction of its evaluable loci:

- unresolvable same-COG paralogs, above `-multiGeneSmplMax` (0.25);
- conspecific SNP-consensus calls, above `-conspGeneSmplMax` (0.05).

With `-excludeMixedStrainSamples 1` (the default) such samples are dropped from that MGS's tree, before any length or prevalence statistic is taken, so they cannot shift the per-locus Q90 length reference or the locus occupancy that the retained samples are judged against. Set it to `0` to keep them. **This never looks at how many loci a sample has.** Per-sample verdicts and their exact fractions are in `<MGS>/sampleQC.tsv`; run-wide counts are in `LOGandSUB/strainRecovery.tsv` (`qc_status`, `ambiguous_failure`, `conspecific_failure`) and in `selection_attrition.tsv` (`qc_excluded_samples`, `qc_excluded_sequences`, `qc_emptied_loci`).

**Too little data** — a clean but sparse sample is judged only by coverage, and it is never flagged as mixed strain. Two gates act in sequence, and it is worth knowing which one is the policy:

- `-MGSminGenesPSmpl` (8) during extraction is a cheap absolute prefilter. It only avoids writing and aligning records that could not clear any inclusion threshold; it cannot tell three well-covered loci from eight fragments, because at that point nothing is aligned yet.
- `-GenesPerSpecies`, `-relativeNTFraction` and `-NTfiltCount` at tree time are the **primary inclusion policy**. They are measured on selected, aligned, informative loci relative to the cohort Q90. With `-enforceSampleCoverage 1` (the default) a sample failing them is removed; with `-placeOnBackbone 1` it is deferred to EPA-ng placement instead.

`-minLociPerMGS` (8) is a third, unrelated threshold: the distinct loci an MGS needs before a tree is attempted at all. It used to share `-MGSminGenesPSmpl`'s value even though it asks about the MGS rather than about a sample.

Before lowering `-MGSminGenesPSmpl`, inspect the retained-locus distribution — and remember that recovered samples still have to clear the tree-time thresholds, so lowering the prefilter alone changes nothing unless `-GenesPerSpecies` is lowered with it:

```bash
awk -F'\t' '$4=="too_few_after_abundance"{n[$5]++} END{c=0; for(k=7;k>=1;k--){c+=n[k]+0; printf "%d loci: %8d rows | -MGSminGenesPSmpl %d recovers %d\n", k, n[k]+0, k, c}}' LOGandSUB/strainRecovery.tsv
```

### 3. Choose sample-inclusion stringency

These values control inclusion after locus QC. `-GenesPerSpecies` is the relative retained-locus floor, `-relativeNTFraction` is the relative informative-nucleotide floor, and `-NTfiltCount` is an absolute informative-nucleotide floor. Larger values retain fewer, better-covered strains.

Per-locus coverage now has two distinct gates. `-GeneLengthMin` (default `0.3`) supplies the high-confidence observations used to select loci and decide whether a sample passes the existing backbone or placement coverage filters. After that decision, `-GeneLengthIncludeMin` (default `0.03`) may recover partial observations into the MSA and phylogeny. Recovered-only observations do not increase the QC locus/NT totals that admitted the sample. They do remain subject to MSAfix masking, column-overlap filtering, and post-alignment locus occupancy/divergence QC.

#### Lenient inclusion

Use this only when sparse sampling is expected and you will inspect the resulting coverage/attrition reports carefully.

```bash
  -MGSminGenesPSmpl 5 \
  -GeneLengthMin 0.20 \
  -GeneLengthIncludeMin 0.03 \
  -GenesPerSpecies 0.10 \
  -relativeNTFraction 0.05 \
  -NTfiltCount 0 \
  -treeLocusBudget 400
```

This lowers the per-sample locus and relative-NT gates but leaves locus QC active. It does not enable EPA-ng placement; add `-placeOnBackbone 1` only if sparse strains should be placed against a backbone.

#### Conservative inclusion

Use this when the primary tree should contain only well-supported strains.

```bash
  -MGSminGenesPSmpl 15 \
  -GeneLengthMin 0.50 \
  -GeneLengthIncludeMin 0.03 \
  -GenesPerSpecies 0.35 \
  -relativeNTFraction 0.20 \
  -NTfiltCount 10000 \
  -treeLocusBudget 400
```

These are example presets, not universal biological cutoffs. Inspect `phylo/taxon_aware_diagnostics.tsv` and the sample-attrition reports before choosing a preset for a new cohort.

### 4. Enable optional sparse-strain placement

Placement is separate from basic sample inclusion. It requires an IQ-TREE backbone and enough overlap with that backbone.

```bash
  -placeOnBackbone 1 \
  -strictBackboneFraction 0.35 \
  -placementGenesPerSpecies 0.04 \
  -placementRelativeNTFraction 0.03 \
  -placementMinOverlap 10000 \
  -epaThreads 2
```

`-placeOnBackbone 1` is opt-in. Eligible sparse strains are placed with EPA-ng; samples that fail the placement coverage or overlap checks stay excluded from the published tree. With `-placeOnBackbone 0` (the default), `-strictBackboneFraction`, `-strictBackboneMinSamples`, and every `-placement*`/EPA-specific filter are inactive; the ordinary tree-inclusion and locus-QC settings still apply. `-redoEPAfilter 1` is a resume-only operation for refreshing retained EPA placement filtering and requires `-placeOnBackbone 1`.

### 5. Redo statistics only

Run the second script directly when trees already exist and only R-based summaries need recomputation:

```bash
perl "$MF4DIR/secScripts/MGS/strain_within_2.2.pl" \
  -GCd /project/gene_catalogue \
  -FMGdir /project/gene_catalogue/within_phylo \
  -map /project/samples.map \
  -MGSmatrix /project/gene_catalogue/Annotation/Abundance/MGS.matL7.txt \
  -submit 1 \
  -redoStrainStats 1
```

Use `-redoPopGenStats 1` instead to redo only population-genetic statistics. It requires `-popGenStats 1`. `-reSubmit 1` remains the broad reset: it clears all within-MGS postprocessing state and both statistics types.

## `strain_within.pl` flag groups

This is a task-oriented summary. See the [complete `strain_within.pl` flag reference](flag_reference.md#strain_withinpl) for every parsed option, type and source default.

| Group | Important flags | Purpose |
|---|---|---|
| Inputs and execution | `-GCd`, `-MGS`, `-outD`, `-map2`, `-MGSabundance`, `-submit`, `-submissionMode` | Define catalogue inputs, output location, metadata, abundance matrix, and submission mode. |
| Resume and redo | `-onlySubmit`, `-onlyMSA`, `-redo` | Resume completed work normally, stop after localized per-locus MSAs, or explicitly rebuild tree outputs, incomplete inputs, or all selected results. |
| Extraction QC | `-maxGenes`, `-treeLocusBudget`, `-MGSminGenesPSmpl`, `-GeneLengthMin`, `-GeneLengthIncludeMin`, `-multiGeneSmplMax`, `-conspGeneSmplMax`, `-excludeMixedStrainSamples`, `-minLociPerMGS` | Use high-coverage loci for QC, optionally recover partial loci after sample admission, and control other strain/locus filters. Set `-maxGenes 0` to remove only the gene-count cap; QC remains active. |
| Tree inclusion | `-GenesPerSpecies`, `-relativeNTFraction`, `-NTfiltCount`, `-enforceSampleCoverage`, `-taxonAwareLocusSelection`, `-taxonAwareRescueMinPrevalence` | The primary sample-inclusion policy, plus taxon-aware locus selection. |
| Sparse placement | `-placeOnBackbone`, `-strictBackboneFraction`, `-placementGenesPerSpecies`, `-placementRelativeNTFraction`, `-placementNTfiltCount`, `-placementMinOverlap`, `-epaThreads` | Optional EPA-ng placement controls. Off by default; enable with `-placeOnBackbone 1`. When disabled, all other options in this group are inactive. |
| Mosaic/outgroups | `-mosaicLoci`, `-mosaicMGS`, `-prepareMosaicLoci`, `-outgroupCoreMinLoci`, `-preferredCoreGenes` | Manage Mosaic discovery and choose broadly supported outgroup loci. |
| Tree resources | `-maxCores`, `-selfMemGb`, `-mosaicMemGb`, `-treeOOMMaxMemGB`, `-rateMergePartitions` | Set scheduler and tree-inference resources. |
| Downstream analyses | `-popGenStats`, `-popGenStrictOutgroup`, `-popGenGeneticCode`, `-popGenCodonStart`, `-popGenSeed`, `-individualVar`, `-DiscTests`, `-ContTests` | Forward population-genetic and association-analysis configuration to postprocessing. |

The `-redo` values are deliberately short: `-redo tree` deletes and rebuilds tree-stage outputs while reusing complete inputs, `-redo input` rebuilds missing or incomplete strain inputs and their dependent trees, and `-redo all` deletes and rebuilds all strain inputs and trees for the selected MGS. The default is `-redo none`. Combine a redo with `-MGSsubset` to restrict it to explicit MGS identifiers.

## `strain_within_2.2.pl` flags

This table covers the options normally used for restarts and analysis. See the [complete `strain_within_2.2.pl` flag reference](flag_reference.md#strain_within_22pl) for every parsed option.

| Flag | Default | Purpose |
|---|---:|---|
| `-GCd`, `-FMGdir`, `-map`, `-MGSmatrix` | required | Gene catalogue, within-phylogeny directory, sample map, and MGS abundance matrix. |
| `-submit` | `1` | Submit analyses; use `0` to generate a dry-run plan. |
| `-cores`, `-Hcores` | `4`, `12` | Standard R-analysis cores and heavy downstream-analysis cores. |
| `-popGenStats` | `1` | Enable per-MGS population-genetic analysis. |
| `-popGenSubsample` | `10,20,30,100,200,500` | Subsample sizes used by PopGenStats. |
| `-popGenStrictOutgroup` | `0` | Require the requested outgroup for population-genetic analysis. |
| `-popGenGeneticCode`, `-popGenCodonStart`, `-popGenSeed` | `1`, `1`, `1` | Genetic-code, codon-start, and reproducibility settings. |
| `-reSubmit` | `0` | Broadly clear and redo within-MGS postprocessing and both summary types. |
| `-redoStrainStats` | `0` | Redo only `strainStats`, its combined overview, and dependent network/treeWAS/figure checkpoints. |
| `-redoPopGenStats` | `0` | Redo only `PopGenStats` stores and population summaries; requires `-popGenStats 1`. |
| `-individualVar`, `-familyVar`, `-groupStabilityVars`, `-DiscTests`, `-ContTests` | empty except `individualVar=AssmblGrps` | Metadata columns used by population genetics, stability, and association analyses. |

## Outputs to review

- `<FMGdir>/strainStats.tsv`: combined strain-statistics overview.
- `<FMGdir>/popGenStats.tsv` and `popGenStats.subsamples.tsv`: population-genetic overviews when enabled.
- `<MGS>/within/strainStats.output.Rds` and `popGenStats.output.Rds`: durable per-MGS result stores.
- `<MGS>/phylo/gene_length_filter.samples.tsv`: per-sample counts and gene lists for both length gates and MSA recovery.
- `<MGS>/phylo/taxon_aware_diagnostics.tsv` and `selection_attrition.tsv`: retained-locus and sample-inclusion decisions, including `qc_excluded_samples` (mixed strain) and `coverage_excluded_samples` (below the inclusion thresholds).
- `<MGS>/phylo/taxon_aware_backbone_eligibility.tsv`: per-sample selected loci, informative NT, and the thresholds they were judged against.
- `<MGS>/sampleQC.tsv`: per-sample mixed-strain verdict with its ambiguous and conspecific locus fractions.
- `<FMGdir>/LOGandSUB/strainGeneLengthFilter.samples.tsv`: run-wide sample report with MGS-qualified dropped/recovered gene lists.
- `<FMGdir>/networks/`, `<FMGdir>/GeneEnrich/`, and `phyloFigures.sto`: downstream strain-only outputs/checkpoints.

For the full tree-engine option reference, see [BuildTree5 flags](flag_reference.md#buildtree5pl). For generic output locations, see [Outputs](outputs.md).
