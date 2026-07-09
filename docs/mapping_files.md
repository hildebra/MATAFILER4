<!-- Documentation navigation -->
[Home](../README.md) | [Quick start](quickstart.md) | [Installation](install.md) | [Configuration](configuration.md) | [Mapping files](mapping_files.md) | [Workflows](common_workflows.md) | [Outputs](outputs.md) | [Flag reference](flag_reference.md) | [FAQ](FAQ.md) | [Glossary](glossary.md)

---

# Mapping files

The MATAFILER4 mapping file is the central description of samples, input paths, grouping structure, and run-level output locations. It is parsed by `readMap()` in `Mods/GenoMetaAss.pm` and should be treated as a strict, tab-delimited control file.

Use Linux line endings where possible. The parser removes general newline characters, but many external tools can silently damage the map structure if they reorder rows, drop lines beginning with `#`, or convert tabs to spaces.

## Minimal valid map

The first row must be the header and must contain `#SmplID`. The second column must be either `Path` or `SmplPrefix`. Other recognised columns can appear in any order.

```text
#SmplID	Path	AssmblGrps	SeqTech
#OutPath	/path/to/output_base
#RunID	my_run
#DirPath	/path/to/raw_reads
Sample01	Sample01_dir	Subject01	hiSeq
Sample02	Sample02_dir	Subject01	hiSeq
```

The run-level output directory is built as:

```text
#OutPath/#RunID/
```

By default, each sample-level output directory is named after `#SmplID`:

```text
#OutPath/#RunID/<SmplID>/
```

Older MATAFILER folder layouts can still be inferred or forced by legacy options, but new analyses should use the sample-ID based folder layout.

## Header rules

`readMap()` applies several hard checks to the header:

| Rule | Consequence |
|---|---|
| The first data row is interpreted as the header. | Do not put comments or metadata above the header. |
| The header must contain `#SmplID`. | The parser stops if it cannot find this column. |
| The header must contain at least one of `Path` or `SmplPrefix`. | The parser stops if neither is present. |
| Either `Path` or `SmplPrefix` must be the second column. | The parser stops otherwise. |
| Other recognised columns can be elsewhere. | Their position is detected by column name. |

Recommended header patterns are therefore:

```text
#SmplID	Path	AssmblGrps	SeqTech
```

or:

```text
#SmplID	SmplPrefix	AssmblGrps	SeqTech
```

Do not set both `Path` and `SmplPrefix` for the same sample row. If both fields are non-empty, MATAFILER4 treats this as unsafe because it may use the same reads twice.

## Required run-level tags

The following tags are parsed from lines beginning with `#` after the header. They must appear before sample rows that depend on them.

| Tag | Required | Meaning |
|---|---:|---|
| `#OutPath` | yes | Base output directory. MATAFILER4 appends `#RunID` unless it is already present. |
| `#RunID` | yes | Run identifier and output subdirectory name below `#OutPath`. |
| `#DirPath` | yes | Base input directory for following sample rows. Can be repeated; each occurrence applies to subsequent samples until the next `#DirPath`. |

Example with two input locations:

```text
#SmplID	Path	AssmblGrps
#OutPath	/path/to/results
#RunID	study_A
#DirPath	/data/run_001
Sample01	Sample01_dir	Subject01
Sample02	Sample02_dir	Subject01
#DirPath	/data/run_002
Sample03	Sample03_dir	Subject02
Sample04	Sample04_dir	Subject02
```

## Optional run-level tags parsed by `readMap()`

| Tag | Value | Meaning |
|---|---|---|
| `#NodeTmpDir` | path | Overrides the node-local temporary directory from `config.txt`. |
| `#GlobalTmpDir` | path | Overrides the globally visible temporary directory from `config.txt`. |
| `#mocatFiltPath` | path | Path to pre-filtered MOCAT reads when raw reads cannot be used. Prefer raw reads unless there is a specific reason. |
| `#illuminaClip` | path | Per-map Illumina adapter/clip file path used by trimming-related routines. |
| `#RelaxSMPLID` | `TRUE` | Relaxes some sample-ID checks. Use only for legacy maps; new maps should not require this. |
| `#WARNING` | `OFF` | Disables selected map safety warnings, including duplicate read-location checks and numeric assembly-group warnings. Use only after manual review. |

`#RelaxSMPLID TRUE` does not make arbitrary sample IDs safe. It only relaxes selected checks; downstream file names, FASTA headers, and matrices are still easier to handle when sample IDs are simple.

## Recognised sample columns

| Column | Required | Meaning | Notes |
|---|---:|---|---|
| `#SmplID` | yes | Unique sample identifier used throughout output directories, sequence headers, and matrices. | Must not be empty or duplicated. |
| `Path` | conditional | Relative sample directory under the current `#DirPath`. | Use when each sample has its own directory. |
| `SmplPrefix` | conditional | File prefix under the current `#DirPath`. | Use when several samples share one directory. |
| `SeqTech` | no | Sequencing technology for primary reads. | Recognised values: `ONT`, `PB`, `ill`, `hiSeq`, `454`, `SLR`, `proto`, `miSeq`, `GAII`, `GAII_solexa`, or empty. |
| `SeqTechSingl` | no | Sequencing technology for single-read input. | Same recognised values as `SeqTech`. |
| `ReadLength` | no | Expected read length. | Usually autodetected; set only if required. |
| `AssmblGrps` | no | Samples with the same value are assembled together. | Use descriptive non-numeric group IDs such as `Subject01`. |
| `MapGrps` | no | Samples with the same value are combined for reference/secondary mapping logic. | Internally prefixed with `M_`. Use sparingly. |
| `FamilyGrps` | no | Family/group label for family-aware downstream routines. | Used by some MAG/MGS-related post-processing. |
| `SupportReads` | no | Additional support reads for hybrid assemblies or special workflows. | Use typed entries such as `PB:/path/to/file.bam`, `ONT:/path/to/file.fastq.gz`, or `mate:/path/to/dir/`. |
| `ExcludeAssembly` | no | `1` excludes the sample from assembly/gene-catalog gene collection. | Defaults to `0`. |
| `EstCoverage` | no | Legacy coverage-estimation switch. | Deprecated; defaults to `0`. |
| `cut5PR1` | no | Number of bases to trim from the 5′ end of read 1. | Defaults to `0`. |
| `cut5PR2` | no | Number of bases to trim from the 5′ end of read 2. | Defaults to `0`. |
| `firstXreadsRd` | no | Read only the first X reads for this sample. | Defaults to `0`, meaning no per-sample limit. |
| `firstXreadsWr` | no | Write only the first X reads for this sample after processing. | Defaults to `0`, meaning no per-sample limit. |
| `ENAdownload` | no | ENA accession or download identifier. | Parsed into the map object; used by download-aware workflows. |
| `SRAdownload` | no | SRA accession or download identifier. | Parsed into the map object; used by download-aware workflows. |

`Path` and `SmplPrefix` are alternative ways to locate primary reads. At least one of these columns must exist in the header, but each sample row should fill only one of them.

## Sample-ID rules

By default, `#SmplID` is deliberately strict because it is reused in output paths, sequence headers, and abundance matrices.

Avoid:

- empty sample IDs
- duplicate sample IDs
- reserved names: `opt`, `altNms`
- `$`
- `_`
- `,`
- `-` unless `#RelaxSMPLID TRUE` is set
- starting with a number unless `#RelaxSMPLID TRUE` is set
- characters outside `a-z`, `A-Z`, `0-9`, and `.` unless `#RelaxSMPLID TRUE` is set

Recommended examples:

```text
P001.T0
P001.T1
Mouse11T0
SoilA01
```

Not recommended:

```text
001_sample-A
sample_01
sample,01
patient 01
```

## Input-location rules

### Using `Path`

Use `Path` when each sample has its own directory below `#DirPath`:

```text
#SmplID	Path	AssmblGrps
#OutPath	/path/to/results
#RunID	path_example
#DirPath	/data/fastq
Sample01	Sample01	Subject01
Sample02	Sample02	Subject01
```

This points MATAFILER4 to:

```text
/data/fastq/Sample01/
/data/fastq/Sample02/
```

### Using `SmplPrefix`

Use `SmplPrefix` when files for several samples are in the same directory and can be distinguished by prefix:

```text
#SmplID	SmplPrefix	AssmblGrps
#OutPath	/path/to/results
#RunID	prefix_example
#DirPath	/data/fastq
Sample01	Sample01_	Subject01
Sample02	Sample02_	Subject01
```

The exact R1/R2/single-read files are then selected by the command-line regular expressions such as `-inputFQregex1`, `-inputFQregex2`, and `-inputFQregexSingle`.

### Duplicate read-location checks

With warnings enabled, MATAFILER4 stops if the same `Path` or `SmplPrefix` appears more than once, because this is likely to process the same reads twice. You can disable this check with `#WARNING OFF`, but this should be reserved for unusual, manually verified maps.

## Grouping rules

### `AssmblGrps`

Samples sharing an `AssmblGrps` value are assembled together. This is useful for longitudinal samples from the same individual or related samples where co-assembly is intended.

```text
#SmplID	Path	AssmblGrps
P001.T0	P001_T0	P001
P001.T1	P001_T1	P001
P002.T0	P002_T0	P002
```

Important details:

- If `AssmblGrps` is empty or absent, MATAFILER4 creates one assembly group per sample.
- Explicit assembly-group IDs should not be purely numeric when warnings are enabled.
- When combining several map files through `readMapS()`, explicit assembly-group IDs must be unique across the maps.
- Keep samples from the same assembly group adjacent where possible. The parser can track distributed groups, but adjacent rows make job submission and troubleshooting easier.

### `MapGrps`

Samples sharing a `MapGrps` value are treated as a mapping group for workflows that combine reads for reference/secondary mapping. Internally, MATAFILER4 prefixes these group IDs with `M_` to distinguish them from assembly groups.

Use this only when you explicitly want technical samples to contribute to a combined mapping unit. For most biological samples, keep `MapGrps` empty.

### `FamilyGrps`

`FamilyGrps` is parsed as an additional grouping label and is mainly useful for downstream routines that need a family or cohort grouping, for example extracting bins or summaries per family.

## Support reads and hybrid inputs

`SupportReads` stores additional read sources, usually for hybrid assemblies. The preferred format is:

```text
<tag>:/absolute/or/resolvable/path
```

Common tags are:

```text
PB:/path/to/pacbio.bam
ONT:/path/to/nanopore.fastq.gz
mate:/path/to/matepair_reads/
```

For multiple support-read entries, use semicolons between typed entries:

```text
PB:/path/to/pb1.bam;PB:/path/to/pb2.bam
```

MATAFILER4 resolves environment variables in paths such as `$PROJECT` and normalises repeated slashes.

For hybrid workflows, primary reads should normally be second-generation reads, with third-generation reads supplied through `SupportReads`.

## Complete example: longitudinal co-assembly

```text
#SmplID	Path	AssmblGrps	SeqTech	ReadLength
#OutPath	/path/to/results
#RunID	longitudinal_run
#DirPath	/data/study_A/fastq
P001.T0	P001_T0	P001	hiSeq	150
P001.T1	P001_T1	P001	hiSeq	150
P002.T0	P002_T0	P002	hiSeq	150
P002.T1	P002_T1	P002	hiSeq	150
```

## Complete example: hybrid Illumina plus PacBio

```text
#SmplID	Path	AssmblGrps	SeqTech	SupportReads
#OutPath	/path/to/results
#RunID	hybrid_run
#DirPath	/data/illumina_fastq
SampleA	SampleA_illumina	SampleA	hiSeq	PB:/data/pacbio/SampleA.bam
SampleB	SampleB_illumina	SampleB	hiSeq	PB:/data/pacbio/SampleB.bam
```

## Complete example: shared directory with prefixes

```text
#SmplID	SmplPrefix	AssmblGrps	SeqTech
#OutPath	/path/to/results
#RunID	prefix_run
#DirPath	/data/all_fastqs
Mouse11T0	Mouse11T0_	Mouse11	hiSeq
Mouse11T1	Mouse11T1_	Mouse11	hiSeq
Mouse12T0	Mouse12T0_	Mouse12	hiSeq
```

## Practical validation checklist

Before launching a full run, check that:

1. The file is tab-delimited, not comma-delimited or space-delimited.
2. The first row is the header and contains `#SmplID`.
3. The second column is either `Path` or `SmplPrefix`.
4. `#OutPath`, `#RunID`, and `#DirPath` are present before sample rows.
5. Each sample row has exactly one primary read locator: either `Path` or `SmplPrefix`.
6. Sample IDs are simple and unique.
7. Explicit `AssmblGrps` are descriptive and non-numeric.
8. Repeated `#DirPath` lines appear before the sample rows they apply to.
9. `SeqTech` values are one of the values recognised by `readMap()`.
10. The map passes a dry run, for example:

```bash
perl $MF4DIR/MATAF4.pl -map mapping_file.map -submit 0 -from 0 -to 1
```

## Common problems detected by `readMap()`

| Error/problem | Likely cause | Fix |
|---|---|---|
| `Could not find "#SmplID" in input map` | Header missing or not first row. | Move the header to the first row and include `#SmplID`. |
| `Expected to find at least "SmplPrefix" or "Path"` | No input-location column in header. | Add `Path` or `SmplPrefix`. |
| `Either "SmplPrefix" or "Path" has to be second column` | Input-location column is present but not second. | Reorder the header. |
| `Provide tag "#OutPath" in map` | Missing run-level output base. | Add `#OutPath`. |
| `Provide tag "#RunID" in map` | Missing run ID. | Add `#RunID`. |
| `Provide tag "#DirPath" in map` | Missing input base directory. | Add `#DirPath` before sample rows. |
| `Double sample ID` | Duplicate `#SmplID`. | Rename one sample. |
| `Found the sample path ... more than once` | Duplicate `Path` or `SmplPrefix`. | Check for unintended duplicated read use. |
| `AssmblGrps are not allowed to be purely numeric` | Assembly group is e.g. `1` or `23`. | Use a descriptive value such as `AG1`, `Subject23`, or `Mouse23`. |

## Notes on editing maps

A mapping file is order-dependent because `#DirPath`, `#OutPath`, and `#RunID` affect subsequent sample rows. Avoid editing workflows that may reorder rows or discard lines beginning with `#`. If using a spreadsheet, copy/paste into a plain-text editor and save as tab-delimited text.
