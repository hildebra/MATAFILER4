<!-- Documentation navigation -->
[Home](../README.md) | [Quick start](quickstart.md) | [Installation](install.md) | [Configuration](configuration.md) | [Mapping files](mapping_files.md) | [Workflows](common_workflows.md) | [Examples](examples.md) | [Profiling tutorial](profiling_tutorial.md) | [Strain-within guide](strainwithin.md) | [Outputs](outputs.md) | [Flag reference](flag_reference.md) | [FAQ](FAQ.md) | [Glossary](glossary.md)

---

# Glossary

| Term | Meaning |
|---|---|
| Assembly-dependent workflow | Workflow where reads are assembled into contigs/scaffolds, genes are predicted, reads are mapped back to assemblies, and gene/MAG/MGS outputs are derived. |
| Assembly-independent workflow | Workflow that profiles reads directly against taxonomic or functional references without relying on metagenomic assembly. |
| `AssmblGrps` | Mapping-file column that groups samples into a shared assembly. Useful for longitudinal or related samples. |
| `MapGrps` | Mapping-file column that groups samples for reference mapping. Usually needed only in specialised workflows. |
| MAG | Metagenome-assembled genome; a genome bin reconstructed from metagenomic assemblies. |
| MGS | Metagenomic species; a dereplicated or clustered species-level genome/bin unit used for abundance, taxonomy and phylogeny. |
| Gene catalog | Non-redundant set of genes clustered across samples, usually at 95% identity in MATAFILER workflows. |
| `pseudoGC` | Run-level folder that stores several assembly-independent or read-based profiling outputs in accessible tabular formats. |
| `ContigStats` | Directory containing coverage, GC, k-mer and marker-gene summaries for contigs/genes. Used by binning and downstream analysis. |
| FMG | A set of marker genes used in several MATAFILER routines. |
| `shrtHD` | MATAFILER-shortened FASTA header format that preserves sample, contig and gene origin. Required for tracking genes across pipeline stages. |
| Stone / `.sto` file | Completion marker file used by MATAFILER4 to decide whether a step has finished and can be skipped on rerun. |
| `LOGandSUB` | Directory containing submitted job scripts and stdout/stderr logs. This is usually the first place to inspect failed runs. |
| `map2tar` | MATAFILER mode for mapping reads to a supplied target reference FASTA. |
| `map2DB` / `map2GC` | Reference/database mapping modes related to `map2tar`, used for database or gene-catalog mapping workflows. |
| Support reads | Additional reads from another technology, such as PacBio, ONT or mate-pair reads, used to support assembly or mapping. |
| Scaled matrix | Normalized abundance matrix, commonly with sample columns summing to 1. |


## Script entry points

- **`MATAF4.pl`**: main sample-level MATAFILER4 pipeline entry point.
- **`geneCat.pl`**: gene catalog construction and gene catalog annotation entry point.
- **`MGS.pl`**: MAG/MGS clustering, abundance and taxonomy entry point.
- **`buildTree5.pl`**: phylogenetic tree and MSA workflow entry point.
