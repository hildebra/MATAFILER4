<!-- Documentation navigation -->
[Home](documentation/README.md) | [Quick start](documentation/quickstart.md) | [Installation](documentation/install.md) | [Configuration](documentation/configuration.md) | [Mapping files](documentation/mapping_files.md) | [Workflows](documentation/common_workflows.md) | [Outputs](documentation/outputs.md) | [Flag reference](documentation/flag_reference.md) | [FAQ](documentation/FAQ.md) | [Glossary](documentation/glossary.md)

---

<p align="center"><img src="/../main/helpers/images/matafiler-logo.png?raw=true" alt="MATAFILER4" width=500 align="center"></img></p>


# MATAFILER4 documentation

MATAFILER4 is a Linux/HPC-oriented metagenomic processing pipeline for raw shotgun metagenomic reads. It supports assembly-dependent workflows for communities that assemble well, such as many host-associated microbiomes, and assembly-independent profiling workflows for highly complex communities such as soil.

This documentation snapshot was updated against `MATAF4.pl` version `4.04`.

## Which mode should I use?

| Use case | Recommended starting point |
|---|---|
| Gut or other host-associated shotgun metagenomes where assemblies are expected to recover substantial read content | [Assembly-dependent workflow](common_workflows.md#assembly-dependent-metagenomic-assembly--gene-catalog) |
| Highly complex communities where assemblies are unlikely to be informative | [Assembly-independent profiling](common_workflows.md#assembly-independent-profiling) |
| Short-read plus ONT/PacBio support data | [Hybrid assemblies](common_workflows.md#hybrid-assemblies) |
| Mapping reads to a defined reference FASTA or database | [map2tar / map2DB / map2GC](common_workflows.md#map2tar-map2db-and-map2gc-reference-mapping) |
| You already have output and need to know what files matter | [Outputs](outputs.md) |
| You need the exact current command-line options | [Flag reference](flag_reference.md) |

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

1. [Installation](documentation/install.md) — software, environments, databases and HPC setup.
2. [Quick start](documentation/quickstart.md) — shortest path from install to a test run.
3. [Configuration](documentation/configuration.md) — `config.txt`, temporary directories and cluster settings.
4. [Mapping files](documentation/mapping_files.md) — required map structure and sample metadata fields.
5. [Common workflows](documentation/common_workflows.md) — assembly-dependent, assembly-independent, hybrid and reference mapping examples.
6. [Outputs](documentation/outputs.md) — file-by-file description of final and intermediate results.
7. [Flag reference](documentation/flag_reference.md) — current options parsed from `MATAF4.pl`, `geneCat.pl`, `MGS.pl` and `buildTree5.pl`.
8. [FAQ](documentation/FAQ.md) — troubleshooting and common failure modes.
9. [Glossary](documentation/glossary.md) — terms used throughout the pipeline.


## Validated flag references

The command-line reference is generated from the uploaded Perl sources for `MATAF4.pl`, `geneCat.pl`, `MGS.pl` and `buildTree5.pl`. See [Flag reference](documentation/flag_reference.md).


## Citing MATAFILER4

**Please cite MATAFILER4 with:**
- Assembly mode: Hildebrand, F. et al. Antibiotics-induced monodominance of a novel gut bacterial order. Gut 68, 1781–1790 (2019).
- Strain mode: Hildebrand, F. et al. Dispersal strategies shape persistence and evolution of human gut bacteria. Cell Host & Microbe 29, 1167-1176.e9 (2021).
- Assembly-independent mode: Bahram, M. et al. Metagenomic assessment of the global diversity and distribution of bacteria and fungi. Environmental Microbiology 23, 316–326 (2021).
- sdm, LCA: Özkurt, E. et al. Microbiome (2022).

Falk Hildebrand <Falk.Hildebrand@gmail.com>

## License

 Copyright (c) 2017-2026 Falk Hildebrand

 MATAFILER4 is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 2 of the License, or
 (at your option) any later version.

 MATAFILER4 is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 See the file LICENSE for more details.

 You should have received a copy of the GNU General Public License
 along with the source code.  If not, see <http://www.gnu.org/licenses/>.

</details>