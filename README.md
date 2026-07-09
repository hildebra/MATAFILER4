<p align="center"><img src="/../main/helpers/images/matafiler-logo.png?raw=true" alt="MATAFILER4" width=500 align="center"></img></p>

<!-- Documentation navigation -->
[Home](README.md) | [Quick start](docs/quickstart.md) | [Installation](docs/install.md) | [Configuration](docs/configuration.md) | [Mapping files](docs/mapping_files.md) | [Workflows](docs/common_workflows.md) | [Outputs](docs/outputs.md) | [Flag reference](docs/flag_reference.md) | [FAQ](docs/FAQ.md) | [Glossary](docs/glossary.md)

---

# MATAFILER4 documentation

MATAFILER4 is a Linux/HPC-oriented metagenomic processing pipeline for raw shotgun metagenomic reads. It supports assembly-dependent workflows for communities that assemble well, such as many host-associated microbiomes, and assembly-independent profiling workflows for highly complex communities such as soil.

This documentation snapshot was updated against `MATAF4.pl` version `4.04`.

## Which mode should I use?

| Use case | Recommended starting point |
|---|---|
| Gut or other host-associated shotgun metagenomes where assemblies are expected to recover substantial read content | [Assembly-dependent workflow](docs/common_workflows.md#assembly-dependent-metagenomic-assembly--gene-catalog) |
| Highly complex communities where assemblies are unlikely to be informative | [Assembly-independent profiling](docs/common_workflows.md#assembly-independent-profiling) |
| Short-read plus ONT/PacBio support data | [Hybrid assemblies](docs/common_workflows.md#hybrid-assemblies) |
| Mapping reads to a defined reference FASTA or database | [map2tar / map2DB / map2GC](docs/common_workflows.md#map2tar-map2db-and-map2gc-reference-mapping) |
| You already have output and need to know what files matter | [Outputs](docs/outputs.md) |
| You need the exact current command-line options | [Flag reference](docs/flag_reference.md) |

## Minimal installation

```bash
git clone https://github.com/hildebra/MATAFILER4.git
cd MATAFILER4
bash helpers/install/installer.sh
micromamba activate MF4
./MATAF4.pl -checkInstall
```

## Minimal test pattern

Create or adapt a mapping file, then run a dry-run first:

```bash
MAP=/path/to/mapping_file.map
perl $MF4DIR/MATAF4.pl -map "$MAP" -assembleMG 2 -submit 0 -from 0 -to 1
```

Only submit to the scheduler after the dry-run has validated paths and configuration:

```bash
perl $MF4DIR/MATAF4.pl -map "$MAP" -assembleMG 2 -submit 1 -from 0 -to 1
```

## Documentation map

1. [Installation](docs/install.md) — software, environments, databases and HPC setup.
2. [Quick start](docs/quickstart.md) — shortest path from install to a test run.
3. [Configuration](docs/configuration.md) — `config.txt`, temporary directories and cluster settings.
4. [Mapping files](docs/mapping_files.md) — required map structure and sample metadata fields.
5. [Common workflows](docs/common_workflows.md) — assembly-dependent, assembly-independent, hybrid and reference mapping examples.
6. [Outputs](docs/outputs.md) — file-by-file description of final and intermediate results.
7. [Flag reference](docs/flag_reference.md) — current options parsed from `MATAF4.pl`, `geneCat.pl`, `MGS.pl` and `buildTree5.pl`.
8. [FAQ](docs/FAQ.md) — troubleshooting and common failure modes.
9. [Glossary](docs/glossary.md) — terms used throughout the pipeline.


## Validated flag references

The command-line reference is generated from the uploaded Perl sources for `MATAF4.pl`, `geneCat.pl`, `MGS.pl` and `buildTree5.pl`. See [Flag reference](docs/flag_reference.md).
