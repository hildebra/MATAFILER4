<p align="center"><img src="/../main/helpers/images/matafiler-logo.png?raw=true" alt="MATAFILER4" width=500 align="center"></img></p>

<p align="center">
  <a href="https://github.com/hildebra/MATAFILER4"><img alt="GitHub repository" src="https://img.shields.io/badge/GitHub-MATAFILER4-181717?logo=github"></a>
  <a href="docs/install.md"><img alt="Installation" src="https://img.shields.io/badge/install-micromamba%20%2B%20installer.sh-blue"></a>
  <a href="docs/quickstart.md"><img alt="Quick start" src="https://img.shields.io/badge/quickstart-dry--run%20first-success"></a>
  <a href="docs/flag_reference.md"><img alt="Flags validated from source" src="https://img.shields.io/badge/flags-validated%20from%20Perl%20source-informational"></a>
  <img alt="Platform" src="https://img.shields.io/badge/platform-Linux%20%2F%20HPC-lightgrey">
  <img alt="Language" src="https://img.shields.io/badge/language-Perl-39457E?logo=perl">
  <img alt="MATAFILER4 version" src="https://img.shields.io/badge/MATAFILER4-v4.04-brightgreen">
</p>

<p align="center">
  <a href="https://github.com/hildebra/MATAFILER4/releases"><img alt="GitHub release" src="https://img.shields.io/github/v/release/hildebra/MATAFILER4?include_prereleases&label=release"></a>
  <a href="https://github.com/hildebra/MATAFILER4/commits/main"><img alt="Last commit" src="https://img.shields.io/github/last-commit/hildebra/MATAFILER4?label=last%20commit"></a>
  <a href="https://github.com/hildebra/MATAFILER4/issues"><img alt="Open issues" src="https://img.shields.io/github/issues/hildebra/MATAFILER4?label=issues"></a>
  <a href="https://github.com/hildebra/MATAFILER4/blob/main/LICENSE"><img alt="License" src="https://img.shields.io/github/license/hildebra/MATAFILER4"></a>
  <a href="https://github.com/hildebra/MATAFILER4/releases"><img alt="Downloads" src="https://img.shields.io/github/downloads/hildebra/MATAFILER4/total?label=downloads"></a>
</p>

<!--
Optional CI badges to enable once matching GitHub Actions workflow files exist.
Replace the workflow filenames if your repository uses different names.

[![Installer smoke test](https://github.com/hildebra/MATAFILER4/actions/workflows/installer-smoke-test.yml/badge.svg)](https://github.com/hildebra/MATAFILER4/actions/workflows/installer-smoke-test.yml)
[![Documentation links](https://github.com/hildebra/MATAFILER4/actions/workflows/docs.yml/badge.svg)](https://github.com/hildebra/MATAFILER4/actions/workflows/docs.yml)
[![Perl syntax checks](https://github.com/hildebra/MATAFILER4/actions/workflows/perl-syntax.yml/badge.svg)](https://github.com/hildebra/MATAFILER4/actions/workflows/perl-syntax.yml)
-->


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
