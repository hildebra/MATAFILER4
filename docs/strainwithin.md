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

### 2. Choose sample-inclusion stringency

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

### 3. Enable optional sparse-strain placement

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

### 4. Redo statistics only

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
| Resume and redo | `-onlySubmit`, `-redo` | Resume completed work normally or explicitly rebuild tree outputs, incomplete inputs, or all selected results. |
| Extraction QC | `-maxGenes`, `-treeLocusBudget`, `-MGSminGenesPSmpl`, `-GeneLengthMin`, `-GeneLengthIncludeMin`, `-multiGeneSmplMax`, `-conspGeneSmplMax` | Use high-coverage loci for QC, optionally recover partial loci after sample admission, and control other strain/locus filters. Set `-maxGenes 0` to remove only the gene-count cap; QC remains active. |
| Tree inclusion | `-GenesPerSpecies`, `-relativeNTFraction`, `-NTfiltCount`, `-taxonAwareLocusSelection`, `-taxonAwareRescueMinPrevalence` | Set backbone inclusion thresholds and taxon-aware locus selection. |
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
- `<MGS>/phylo/taxon_aware_diagnostics.tsv` and `selection_attrition.tsv`: retained-locus, sample-inclusion, and aggregate length-recovery decisions.
- `<FMGdir>/LOGandSUB/strainGeneLengthFilter.samples.tsv`: run-wide sample report with MGS-qualified dropped/recovered gene lists.
- `<FMGdir>/networks/`, `<FMGdir>/GeneEnrich/`, and `phyloFigures.sto`: downstream strain-only outputs/checkpoints.

For the full tree-engine option reference, see [BuildTree5 flags](flag_reference.md#buildtree5pl). For generic output locations, see [Outputs](outputs.md).
