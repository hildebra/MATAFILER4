<!-- Documentation navigation -->
[Home](README.md) | [Quick start](quickstart.md) | [Installation](install.md) | [Configuration](configuration.md) | [Mapping files](mapping_files.md) | [Workflows](common_workflows.md) | [Outputs](outputs.md) | [Flag reference](flag_reference.md) | [FAQ](FAQ.md) | [Glossary](glossary.md)

---

# Mapping files

The mapping file is the central description of samples, input paths and output locations. It must be a tab-delimited text file with Linux line endings.

## Required structure

The first row must start with `#SmplID`. Other columns can appear in any order.

```text
#SmplID	Path	SmplPrefix	AssmblGrps
#OutPath	/path/to/output_base
#RunID	my_run
#DirPath	/path/to/raw_reads
Sample01	Sample01_dir	Sample01_	S01
Sample02	Sample02_dir	Sample02_	S02
```

The run-level output directory is:

```text
#OutPath/#RunID/
```

Each sample-level output folder is named after `#SmplID`.

## Required map tags

| Tag | Meaning |
|---|---|
| `#OutPath` | Base output directory. Ensure there is enough space. |
| `#RunID` | Name of the run directory under `#OutPath`; also a global run identifier. |
| `#DirPath` | Base directory for subsequent sample `Path` entries. Can be used multiple times if raw reads are split across locations. |

## Common sample columns

| Column | Meaning | Recommendation |
|---|---|---|
| `#SmplID` | Unique sample identifier used throughout outputs and matrix columns | Use short alphanumeric IDs; avoid special characters. |
| `Path` | Relative folder path under the current `#DirPath` | Use when each sample has its own directory. |
| `SmplPrefix` | File prefix for sample read files | Use when several samples share a directory. |
| `AssmblGrps` | Samples with the same value are assembled together | Useful for longitudinal samples from the same subject. |
| `MapGrps` | Samples with the same value are combined for reference mapping | Use rarely; often easier to place technical read files together. |
| `SeqTech` | Sequencing technology: `ONT`, `PB`, `ill`, `miSeq`, `hiSeq`, `GAII`, `GAII_solexa`, `proto`, `454`, `AVITI`, `SLR` | Helps select filtering and assembly behaviour. |
| `ReadLength` | Expected read length | Usually autodetected; set manually only when necessary. |
| `SupportReads` | Additional reads such as `PB:/path/to/file`, `ONT:/path/to/file` or `mate:/path/to/file` | Used for hybrid/support-read workflows. |
| `ExcludeAssembly` | `1` excludes the sample from assembly | Useful for samples that should still contribute to mapping/profiling. |
| `cut5PR1`, `cut5PR2` | Trim bases from the 5′ end of read 1 or read 2 | Advanced per-sample trimming. |
| `firstXreadsRd`, `firstXreadsWr` | Limit reads read/written for a sample | Useful for testing or subsampling. |

## Optional run-level tags

| Tag | Meaning |
|---|---|
| `#NodeTmpDir` | Overrides `nodeTmpDir` from `config.txt`. |
| `#GlobalTmpDir` | Overrides `globalTmpDir` from `config.txt`. |
| `#mocatFiltPath` | Path to pre-filtered MOCAT reads when raw reads cannot be used. |
| `#RelaxSMPLID` | Use with caution to relax sample-ID format checks. |
| `#WARNING` | Use with caution; can disable stopping on some mapping-file warnings. |

## Practical recommendations

- Create the mapping file in a spreadsheet only if useful, then copy/paste to a tab-delimited text file.
- Do not load/save the mapping file through tools that may remove comment-like `#DirPath` rows or reorder rows.
- Keep samples belonging to the same `AssmblGrps` block together where possible.
- Validate with `-submit 0 -from 0 -to 1` before submitting a full run.
