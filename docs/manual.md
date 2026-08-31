<!-- Documentation navigation -->
[Home](../README.md) | [Quick start](quickstart.md) | [Installation](install.md) | [Configuration](configuration.md) | [Mapping files](mapping_files.md) | [Workflows](common_workflows.md) | [Examples](examples.md) | [Profiling tutorial](profiling_tutorial.md) | [Strain-within guide](strainwithin.md) | [Outputs](outputs.md) | [Flag reference](flag_reference.md) | [FAQ](FAQ.md) | [Glossary](glossary.md)

---

# MATAFILER4 manual

This page is a compact manual and navigation point. The previous long manual was retained as [`manual_legacy.md`](manual_legacy.md) for comparison. Current command-line options are validated against the repository Perl sources in [Flag reference](flag_reference.md).

## Core pipeline phases

MATAFILER4 commonly runs in four script-level phases:

1. **Sample-level processing with `MATAF4.pl`**: read detection, preprocessing, host filtering, assembly, read mapping, contig statistics, binning and optional read-based taxonomic/functional profiling.
2. **Gene catalog construction with `secScripts/geneCat.pl`**: clustering predicted genes across samples, creating gene abundance matrices, assigning functions and preparing genome/species-resolved outputs.
3. **MGS/MAG processing with `secScripts/MGS.pl`**: dereplication, immediate between-MGS marker-protein tree inference, MGS abundance/taxonomy and optional strain-oriented processing.
4. **Phylogenetic analysis with `secScripts/phylo/buildTree5.pl`**: MSA generation, tree inference and optional population-genetic analyses.

## Recommended reading order

1. [Quick start](quickstart.md)
2. [Configuration](configuration.md)
3. [Mapping files](mapping_files.md)
4. [Common workflows](common_workflows.md)
5. [Outputs](outputs.md)
6. [Flag reference](flag_reference.md)
7. [FAQ](FAQ.md)

## Source-validated flag references

The reference covers every command-line entry point, and each script renders its own
`-help` directly from it:

- [`MATAF4.pl` flags](flag_reference.md#mataf4pl)
- [`geneCat.pl` flags](flag_reference.md#genecatpl)
- [`MGS.pl` flags](flag_reference.md#mgspl)
- [`strain_within.pl` flags](flag_reference.md#strain_withinpl)
- [`strain_within_2.2.pl` flags](flag_reference.md#strain_within_22pl)
- [`buildTree5.pl` flags](flag_reference.md#buildtree5pl)

## Current vs legacy option names

Prefer the clearer current names in new `MATAF4.pl` documentation and examples:

| Prefer | Legacy alias still accepted by `MATAF4.pl` |
|---|---|
| `-assemblCores` | `-spadesCores` |
| `-assemblMemory` | `-spadesMemory` |
| `-assemblyKmers` | `-spadesKmers` |
| `-Binner` | `-MetaBat2`, `-binSpeciesMG` |
| `-filterHostRds` | `-filterHumanRds` |
| `-remap2assembly` | `-redoMap2assembly`, `-redoMapping` |

Do not carry these aliases across scripts unless the receiving script parses them. For example, `geneCat.pl` accepts `-binSpeciesMG`, but not `-Binner`; `buildTree5.pl` accepts `-aa`, but not `-faa`.

## Output documentation

The output description has been moved to [Outputs](outputs.md), including primary files, completion checklists, sample-level assembly outputs, gene catalog files, functional annotations and MGS/MAG outputs.
