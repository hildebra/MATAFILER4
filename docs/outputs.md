<!-- Documentation navigation -->
[Home](../README.md) | [Quick start](quickstart.md) | [Installation](install.md) | [Configuration](configuration.md) | [Mapping files](mapping_files.md) | [Workflows](common_workflows.md) | [Outputs](outputs.md) | [Flag reference](flag_reference.md) | [FAQ](FAQ.md) | [Glossary](glossary.md)

---

# MATAFILER4 outputs

MATAFILER4 creates many intermediate and final files. For most downstream analyses, start with the following outputs.

## Start here: primary files

| Output | Typical location | Use first for |
|---|---|---|
| `metagStatsReport.html` | `#OutPath/#RunID/` | Browser-based run overview and QC. |
| `metagStats.txt` | `#OutPath/#RunID/` | Text summary of run completion and sample statistics. |
| `Matrix.mat.scaled.gz` | Gene catalog directory passed to `geneCat.pl -GCd` | Normalized gene abundance matrix. |
| `Matrix.mat.gz` | Gene catalog directory passed to `geneCat.pl -GCd` | Raw gene count matrix. |
| Functional abundance matrices | Gene catalog annotation folders or `pseudoGC/`, depending on workflow | KEGG, eggNOG, CAZy, ARG and other functional profiles. |
| MGS abundance and taxonomy tables | Gene catalog MGS/MAG output folders | Genome-resolved and species-level interpretation. |
| Sample-level `LOGandSUB/` | `#OutPath/#RunID/<SmplID>/LOGandSUB/` | Debugging failed or incomplete jobs. |

## Completion checklists

### Assembly-dependent MATAF4.pl run

```text
#OutPath/#RunID/metagStats.txt
#OutPath/#RunID/metagStatsReport.html
#OutPath/#RunID/<SmplID>/assemblies/
#OutPath/#RunID/<SmplID>/assemblies/metag/genePred/
#OutPath/#RunID/<SmplID>/assemblies/metag/ContigStats/
#OutPath/#RunID/<SmplID>/mapping/
```

### Gene catalog run

```text
<GCd>/compl.incompl.95.fna
<GCd>/compl.incompl.95.prot.faa
<GCd>/Matrix.mat.gz
<GCd>/Matrix.mat.scaled.gz
<GCd>/Anno/
<GCd>/Binning/ or <GCd>/Bin_SB/  # when MAG/MGS processing was enabled
```


This page summarizes where MATAFILER4 writes its major outputs and which files are most useful for downstream analyses. MATAFILER4 can run in assembly-dependent and assembly-independent modes, so not every file listed here is produced in every run.

## Quick output map

| Output class | Typical location | What it contains | Usually used for |
|---|---|---|---|
| Run-level output | `#OutPath/#RunID/` | Sample folders, logs, assembly-independent profiles, and run summaries | Checking run completion; finding sample-level outputs |
| Sample assemblies | `#OutPath/#RunID/<SmplID>/assemblies/` | Assembly FASTA files, assembly statistics, filtered scaffolds | Assembly QC; input for gene prediction and binning |
| Gene predictions | `#OutPath/#RunID/<SmplID>/assemblies/metag/genePred/` | Nucleotide and amino-acid genes, GFF files, per-contig gene tables | Building the gene catalog; linking genes back to contigs |
| Contig and gene coverage | `#OutPath/#RunID/<SmplID>/assemblies/metag/ContigStats/` | Coverage, GC, k-mer, marker-gene and per-gene statistics | Abundance estimation; binning; QC |
| Read mappings | `#OutPath/#RunID/<SmplID>/mapping/` | Read-to-assembly mappings, normally CRAM/SAM/BAM-derived files | Coverage estimation; troubleshooting mapping |
| Assembly-independent profiles | `#OutPath/#RunID/pseudoGC/` and related run folders | miTAG, mOTUs, MetaPhlAn, Kraken2 and functional read-based profiles, depending on flags | Taxonomic/functional profiling without assembly |
| Gene catalog | Path given to `geneCat.pl -GCd` | Non-redundant genes, gene abundance matrices, annotations and MGS outputs | Main integrated gene-level analysis |
| Functional annotations | `<gene_catalog>/Anno/Func/` | KEGG, eggNOG, CAZy, ARG and other functional abundance matrices | Functional interpretation |
| Taxonomic annotations | `<gene_catalog>/Anno/Tax/` | Taxonomic abundance matrices from gene-level or MGS-level assignments | Taxonomic interpretation |
| MAG/MGS output | `<gene_catalog>/Bin_SB/` or `<gene_catalog>/Binning/` | MAG/MGS clusters, quality estimates, taxonomies, phylogenies, MGS abundance matrices | Genome-resolved analysis and strain/species tracking |

## General conventions

Most final result tables are tab-separated feature abundance matrices. Rows are features, columns are samples, and values are abundances or counts. Features can be genes, functions, taxa, MAGs, MGS, or marker genes.

Common conventions:

- `*.mat*` usually denotes a matrix.
- `*.gz` files are gzip-compressed.
- `L0`, `L1`, `L2`, etc. denote hierarchical levels. `L0` is usually the most specific level and higher numbers are broader summaries, unless a tool-specific file states otherwise.
- `?` usually denotes an unknown or unresolved assignment at a given level.
- `-1` usually denotes unassigned genes or reads.
- `*.scaled*` matrices are normalized, commonly so that each sample column sums to 1.
- Raw or intermediate assignment files are useful for debugging, but summarized matrices are usually preferable for downstream analysis.

## Run-level output directory

The main run directory is defined by the mapping file:

```text
#OutPath /path/to/output_base
#RunID   run_name
```

The run-level output is therefore:

```text
#OutPath/#RunID/
```

This directory typically contains one subdirectory per sample, global log and submission folders, and outputs from assembly-independent profiling. Each sample directory is named after the `#SmplID` value from the mapping file, so the choice of sample IDs affects all downstream file names and matrix column names.

## Sample-level assembly outputs

For an assembly-dependent run, each sample or assembly group receives its own output folder:

```text
#OutPath/#RunID/<SmplID>/
```

Common subdirectories are:

```text
assemblies/
input_fil.txt
input_raw.txt
LOGandSUB/
mapping/
```

The exact files depend on the run flags. A minimal completed assembly usually contains the assembly FASTA files, assembly statistics, contig statistics, gene predictions and mapping-derived coverage information.

### `assemblies/`

The `assemblies/` directory contains the metagenomic assembly and assembly statistics. Common files include:

| File | Meaning |
|---|---|
| `scaffolds.fasta.gz` | Final scaffold FASTA, often compressed. |
| `scaffolds.fasta.filt.gz` | Scaffold FASTA after length filtering. |
| `AssemblyStats.txt` | Assembly statistics for the full assembly. |
| `AssemblyStats.500.txt` | Assembly statistics after filtering to scaffolds above the configured length threshold, often 500 bp. |
| `Ass.done.sto` | Stone/checkpoint file indicating that the assembly step finished. |
| `metag/` | Main assembly post-processing directory. |

### `assemblies/metag/ContigStats/`

`ContigStats/` stores abundance and sequence-composition files used by later steps, especially binning and gene catalog construction.

Common files include:

| File | Meaning |
|---|---|
| `Coverage.pergene` | Per-gene coverage estimated from mapped reads. |
| `Coverage.count_pergene` | Per-gene read counts. |
| `Coverage.median.pergene` | Median per-base coverage per gene. |
| `Coverage.percontig` | Per-contig coverage. |
| `Coverage.median.percontig` | Median per-base coverage per contig. |
| `Coverage.window` | Coverage in windows across contigs. |
| `GeneStats.txt` | Summary statistics for predicted genes. |
| `FMG` | Forty universal marker genes used in several MATAFILER routines. |
| `ess100genes` | Essential marker-gene set. |
| `scaff.GC.gz`, `scaff.pergene.GC3.gz` | GC and codon-position GC statistics. |
| `scaff.4mer.gz`, `scaff.pergene.4mer.gz` | k-mer composition statistics. |

### `assemblies/metag/genePred/`

`genePred/` contains the predicted genes for the sample assembly.

| File | Meaning |
|---|---|
| `genes.shrtHD.fna` | Nucleotide sequences of predicted genes with MATAFILER-shortened headers. |
| `proteins.shrtHD.faa` | Amino-acid sequences of predicted proteins with MATAFILER-shortened headers. |
| `genes.gff` | GFF file linking predicted genes to their contigs and coordinates. |
| `genes.per.ctg` | Number and identity of genes per contig. |

The shortened gene headers are important because they encode the sample, contig and gene position. A typical header has the form:

```text
>DLF0001__C1263_L=114017=_35
```

This identifies the originating sample (`DLF0001`), the contig (`C1263`), the contig length (`L=114017`) and the gene number on that contig (`35`). These identifiers allow MATAFILER to trace genes back from the gene catalog to individual samples and contigs.

### `mapping/`

The `mapping/` directory stores sample-named read mappings to the assembly, usually in compressed alignment formats such as CRAM. If `-saveReadsNotMap2Assembly 1` was used, unmapped reads are saved in a subdirectory such as:

```text
mapping/unaligned/
```

A sample is ready to be incorporated into a gene catalog when assembly, gene prediction, contig statistics and mapping-derived coverage outputs are complete.

## Gene catalog outputs

The gene catalog output directory is not automatically the same as `#OutPath/#RunID/`. It is set when running `geneCat.pl`, for example:

```bash
perl $MF4DIR/secScripts/geneCat.pl \
  -map /path/to/mapping_file.map \
  -GCd /path/to/results/genecat \
  -mem 200 -cores 24 -clusterID 95
```

The directory passed to `-GCd` stores the integrated gene catalog and its derived summaries.

### Core gene catalog files

Every current catalog also stores two durable workflow descriptors under
`LOGandSUB/`:

| File | Meaning |
| --- | --- |
| `inmap.txt` | Map manifest. It contains one relocatable catalog-local `map.N.txt` path per line, so catalogs constructed from multiple mapping files retain the complete ordered set. Older single-map catalogs whose `inmap.txt` is itself a mapping table remain readable. |
| `catalog.sha256` | Persistent catalog identity. The same SHA-256 value is reused by MGS and strain workflows, including after the complete catalog directory is moved. |

`GCmaps.inf` is still emitted as a write-only, comma-separated compatibility and
logging artifact, but current code does not read it. Existing catalogs that have
`map.N.txt` copies but no `inmap.txt` are migrated automatically when their maps
are next resolved.

| File | Meaning |
|---|---|
| `compl.incompl.95.fna` | Nucleotide representation of non-redundant genes clustered at 95% identity. |
| `compl.incompl.95.prot.faa` | Amino-acid representation of the same clustered genes. |
| `*.GC` | GC content for gene catalog entries. |
| `*.kmer` | k-mer frequencies for gene catalog entries. |
| `*.length` | Gene lengths. |
| `Matrix.mat.gz` | Main gene abundance matrix: read counts per gene per sample. |
| `Matrix.mat.scaled.gz` | Normalized gene abundance matrix, typically with each sample column summing to 1. |
| `Matrix.mat.sum` | Column sums of `Matrix.mat.gz`, useful for normalization checks. |
| `Matrix.FMG.mat` | Subset of the gene matrix containing universal bacterial marker genes. |
| `Mat.cov.mat.gz` | Gene coverage matrix, typically reads mapped per gene length. |
| `Mat.med.mat.gz` | Median per-base coverage matrix for genes. |

For routine biological interpretation, the summarized taxonomic, functional and MGS-level matrices are usually more convenient than the raw full gene matrix.

## Functional annotation outputs

Functional annotations are usually stored in:

```text
<gene_catalog>/Anno/Func/
```

MATAFILER can summarize functional annotations from several databases. Common prefixes include:

| Prefix | Database / category |
|---|---|
| `KGM` or `KEGG` | KEGG orthologs or modules |
| `NOG` | eggNOG |
| `CZy` | CAZy carbohydrate-active enzymes |
| `ABR` / `ABRc` | Antibiotic resistance annotations |
| `TCDB` | Transporter Classification Database |
| `PAB`, `PTV`, `ACL`, `MOH` | Other optional functional databases, depending on installation and flags |

For hierarchical databases, MATAFILER writes one matrix per level, for example:

```text
CZyL0.txt
CZyL1.txt
NOG.L0.txt
NOG.L1.txt
```

`L0` is normally the most specific annotation level. Higher levels are broader summaries.

### KEGG and module outputs

Module-style outputs can be found under:

```text
<gene_catalog>/Anno/Func/modules/
```

Common files include:

| File | Meaning |
|---|---|
| `KEGG.mat` | Module abundance matrix. |
| `KEGG.descr` | Module hierarchy and descriptions. |
| `KEGG.KOused` | KEGG orthologs used to infer each module in each sample. |
| `KEGG.MODscore` | Module completeness score per sample. |
| `KGML0.txt` | KEGG ortholog abundance matrix. |

### Gene-to-function assignments

To inspect which genes were assigned to which functions, use assignment files such as:

```text
Anno/Func/DIAass_<db>.srt.gzgeneAss.gz
```

where `<db>` is the database of interest. These files list gene catalog IDs and their assigned lowest-level functional annotations. A single gene can have more than one annotation.

## Taxonomic annotation outputs

Taxonomic summaries are usually stored in:

```text
<gene_catalog>/Anno/Tax/
```

These are feature abundance matrices derived from gene abundances and taxonomic annotation. Feature names are usually full or partial lineages separated by semicolons. Unresolved ranks are denoted with `?`.

Common taxonomic outputs include:

| Output | Meaning |
|---|---|
| `GTDBmg` or `GTDBmg_MGS` | GTDB-marker-gene-derived taxonomy outputs; some older files may be deprecated depending on the current workflow. |
| `krak*` | Kraken/Kraken2-derived taxonomic summaries. |
| `specI*` | specI-style marker-gene taxonomic summaries, historically used in MATAFILER. |
| `MGS.matL*.txt` | MGS-based abundance matrices by taxonomic level, in current MGS workflows. |

For current MGS-based taxonomic abundance, prefer the matrices generated from the MGS workflow, commonly under:

```text
<gene_catalog>/Bin_SB/Annotation/Abundance/
```

or the equivalent current MGS output directory configured with `MGS.pl -outD`.

## MAG and MGS outputs

MATAFILER reconstructs MAGs from single-sample assemblies and clusters them into MGS, usually at species-level resolution. Depending on the MATAFILER version and selected workflow, these outputs are commonly under one of:

```text
<gene_catalog>/Bin_SB/
<gene_catalog>/Binning/
```

The output directory can also be changed with:

```bash
MGS.pl -outD <output_directory>
```

### Important MGS files

| File | Meaning |
|---|---|
| `MB2.clusters.ext.can.Rhcl.matL0.txt` | MGS abundance across samples. |
| `MB2.clusters.ext.can.Rhcl` | Assignment of gene-catalog genes to MGS, combining MetaBAT2, canopy clustering and MATAFILER post-filtering. |
| `MB2.clusters.ext.can.Rhcl.cm` | CheckM quality and contamination estimates for MGS. |
| `GTDBTK.tax` | GTDB-Tk taxonomic assignment; generally the preferred modern taxonomy output when available. |
| `kraken2.LCA` | Kraken2 lowest-common-ancestor taxonomic assignment. |
| `between_phylo/phylo/IQtree_allsites.treefile` | De novo phylogeny of MGS. |
| `between_phylo/phylo/IQtree_allsites.pdf` | Abundance-annotated rendering, generated after the MGS abundance matrix is available. |
| `RhclClust/*.faa` | Proteins found in each MGS, stored separately per MGS. |
| `MAG.MB2.assStat.summary` | Summary of MetaBAT2 bins found in each sample. |
| `Bin_SB/MAGvsGC.txt.gz` | Links MAGs, MGS, marker genes and other gene-catalog genes. |
| `Bin_SB/Annotation/GTDBmg_MGS/` | MGS-specific SpecI annotation and abundance intermediates. These are kept inside the selected binner output so different MGS definitions cannot reuse one another's results. |

Sparse runs can finish successfully without manufacturing cluster or tree data:

- `Canopy_AC/SKIPPED.txt` (or `Canopy/SKIPPED.txt`) explains that Canopy was skipped because the abundance matrix had too few sample columns, or that Canopy completed without finding clusters.
- `Bin_<binner>/NO_MGS.txt` explains why no usable MGS could be reconstructed (for example, no minimally usable MAGs or no retained core genes).
- `Bin_<binner>/between_phylo/SKIPPED.txt` records that fewer than three marker-bearing MGS were available for a meaningful between-MGS phylogeny.

These files represent completed, expected low-cardinality outcomes. Missing or partial paired outputs (for example, only one of the Canopy cluster and profile files) remain errors.

The between-MGS tree starts immediately after Stage I publishes the MGS core set and runs concurrently with taxonomy and abundance. It is inferred from concatenated amino-acid FMG alignments with per-locus partitions, partition-aware ModelFinder, 1,000 ultrafast bootstrap replicates, and minimum marker/column occupancy filters. The topology is unrooted unless the standalone launcher is given suitable reference genomes and an explicit downstream rooting strategy.

### MGS abundance matrices

Current MGS workflows commonly write abundance matrices by taxonomic level:

```text
Bin_SB/Annotation/Abundance/MGS.matL0.txt
Bin_SB/Annotation/Abundance/MGS.matL1.txt
...
Bin_SB/Annotation/Abundance/MGS.matL7.txt
```

Typical interpretation:

| File | Interpretation |
|---|---|
| `MGS.matL0.txt` | Broadest rank, often domain-level. |
| `MGS.matL1.txt` - `MGS.matL6.txt` | Intermediate taxonomic ranks. |
| `MGS.matL7.txt` | Individual MGS-level abundance. |

At MGS level, features ending in `?` are taxa inferred through marker-gene LCA assignments for which no reconstructed MGS was available.

### MAG/MGS gene content

`Bin_SB/MAGvsGC.txt.gz` links MAGs and MGS to their gene content. Important columns include:

| Column | Meaning |
|---|---|
| `MAG` | MAG identifier. |
| `MGS` | MGS identifier to which the MAG belongs. |
| `Representative4MGS` | Marks the representative MAG for an MGS, often with `*`. |
| Marker-gene columns | Gene IDs matching individual marker genes. |
| `other_genes` | Comma-separated list of other genes in the MAG. Double commas indicate the start of a new contig. |

This file is useful for linking MGS or MAGs to functional annotations from `Anno/Func/`.

## Assembly-independent outputs

Assembly-independent profiling uses raw or filtered reads directly and is useful when assemblies are not feasible, for example in highly complex soil metagenomes.

Depending on flags, MATAFILER can write outputs from:

| Method | Typical purpose |
|---|---|
| Kraken2 | Read-based taxonomic profiling. |
| mOTUs | Marker-gene-based taxonomic profiling. |
| MetaPhlAn | Marker-gene-based taxonomic profiling. |
| miTAG / ribosomal profiling | Ribosomal small- or large-subunit profiling. |
| DIAMOND functional profiling | Read-based functional profiling against selected databases. |

These outputs are commonly summarized in or near:

```text
#OutPath/#RunID/pseudoGC/
```

and related run-level folders. The precise file names depend on flags such as `-profileFunct`, `-DiaDBs`, `-profileRibosome`, `-profileMOTU2`, `-profileMetaphlan2` or `-profileKraken`.

## Recommended files for common downstream questions

| Question | Start with these files |
|---|---|
| Did my samples assemble and map correctly? | `metagStats.txt`, `metagStatsReport.html`, sample `LOGandSUB/`, `AssemblyStats*.txt`, `Coverage.*` files. |
| What genes are present and abundant? | `<gene_catalog>/Matrix.mat.gz` and `<gene_catalog>/Matrix.mat.scaled.gz`. |
| What functions are present? | `<gene_catalog>/Anno/Func/*L0.txt` and broader `L*.txt` summaries. |
| What taxa are present? | Current MGS abundance matrices under `Bin_SB/Annotation/Abundance/`, or selected `Anno/Tax/` matrices. |
| Which MGS are present? | `MGS.matL7.txt` and MGS taxonomy files such as `GTDBTK.tax`. |
| Which genes belong to a MAG or MGS? | `Bin_SB/MAGvsGC.txt.gz` and `MB2.clusters.ext.can.Rhcl`. |
| What is the taxonomy of MGS? | `GTDBTK.tax`, `kraken2.LCA`, and MGS abundance feature names. |
| What are the proteins in an MGS? | `Binning/RhclClust/*.faa` or equivalent `Bin_SB` MGS protein folders. |

## Deprecated or lower-priority outputs

Some historical outputs remain in the directory structure for compatibility. Prefer current GTDB-Tk and MGS workflow outputs where available.

- `Anno/Tax/GTDBmg_MGS`: older gene-catalog/specI-style MGS abundance summaries; often superseded by `Bin_SB/Annotation/Abundance/`.
- `Anno/Tax/GTDBmg`: older marker-gene-only taxonomic summaries.
- `EstCoverage`: historical k-mer-based coverage estimation; generally deprecated.
- Intermediate DIAMOND folders such as `CNT_*` and `DIAass_*`: useful for debugging, but final `Anno/Func/*` matrices are usually easier to use.

## Minimal output checklist

Before moving to downstream analyses, check that the expected outputs exist for the workflow you ran.

For assembly-dependent gene catalog analysis:

1. Sample folders exist under `#OutPath/#RunID/<SmplID>/`.
2. Each included sample has an assembly and gene prediction outputs.
3. Coverage files exist under `assemblies/metag/ContigStats/`.
4. The gene catalog directory contains `Matrix.mat.gz` and annotation folders.
5. Functional summaries exist in `Anno/Func/` if functional annotation was requested.
6. Taxonomic or MGS summaries exist in `Anno/Tax/`, `Binning/`, or `Bin_SB/`, depending on the workflow.
7. Run logs do not indicate failed jobs that were not rerun.

For assembly-independent analysis:

1. Requested read-based profilers completed.
2. Corresponding abundance matrices exist in `pseudoGC/` or the relevant run-level output folder.
3. `metagStats.txt` and `metagStatsReport.html` are present if post-processing completed.
