<!-- Documentation navigation -->
[Home](../README.md) | [Quick start](quickstart.md) | [Installation](install.md) | [Configuration](configuration.md) | [Mapping files](mapping_files.md) | [Workflows](common_workflows.md) | [Outputs](outputs.md) | [Flag reference](flag_reference.md) | [FAQ](FAQ.md) | [Glossary](glossary.md)

---

# Configuration

MATAFILER4 stores most site-specific settings in `config.txt`. The pipeline can also override some paths from the mapping file via `#GlobalTmpDir` and `#NodeTmpDir`.

## Required path settings

| Setting | Meaning | Notes |
|---|---|---|
| `MFLRDir` | MATAFILER4 installation directory | Should point to the directory containing `MATAF4.pl`. |
| `DBDir` | Main database directory | Used by several profiling, annotation and classification steps. |
| `globalTmpDir` | Shared scratch visible from all compute nodes | Critical for exchanging intermediate files between jobs. |
| `nodeTmpDir` | Node-local temporary directory | Used for fast local I/O where supported by the cluster. |

## Temporary directory strategy

MATAFILER4 creates many large temporary files. A stable configuration should use:

- a globally visible scratch area for files shared between jobs;
- a node-local temporary directory for high-I/O steps;
- enough local disk space for assembly, mapping, binning and Diamond jobs;
- automatic scratch cleanup enabled unless debugging.

Relevant flags include `-globalTmpDir`, `-nodeTmpDir`, `-nodeHDDspace`, `-reduceScratchUse`, `-rm_tmpdir_reads` and `-rm_tmpInput`.

With `-reduceScratchUse 1`, cleanup is released only after final assembly
publication, ContigStats, configured binning, and configured ConsSNP/variant
work are either already complete or represented by successful terminal scheduler
dependencies. The cleanup worker verifies their published outputs again before
removing sample-owned BAM/CRAM indexes, stale SNP BED files, or the sample
scratch directory. Mapper and FASTA indexes adjacent to a generated assembly
are retained until every member of that assembly group is complete. Indexes
adjacent to external references are never removed.

For assembly-independent runs (`-assembleMG 0`), assembly and ContigStats
barriers are not required. Cleanup waits for the configured profiling and
reference-mapping jobs, then removes only sample-owned temporary files and
indexes. Binning, assembly consensus SNP, and structural-variant options require
an enabled assembly mode and are rejected during option validation.

## Scheduler settings

MATAFILER4 can autodetect supported schedulers, but this can be overridden with:

```bash
-submSystem slurm
```

Supported systems in the documentation and source are Slurm, SGE/qsub and LSF/bsub. For large runs, use:

```bash
-maxConcurrentJobs 600
-killDepNever 1
```

`-killDepNever 1` can remove Slurm jobs stuck in dependency states, but use it only when this matches your local scheduler policy.

Sites where compute nodes do not normally have outbound network access can set
`netQueue` in `config.txt` to a network-enabled queue or Slurm partition. Code
that sets `useNetQueue` on its submission options will use that queue and the
configured `longTime` wall time for the next job. If `netQueue` is empty or
unset, MATAFILER4 falls back to `mediumQueue`.

## Database setup

GTDB and GTDB-Tk databases are used for MAG and MGS-related classification. The installer provides:

```bash
helpers/install/get_gtdb.py
```

Example:

```bash
./get_gtdb.py all -v 226 -t /path/to/download/to -d /path/to/extract/to --tk split
```

Delete the download directory only after confirming that the extracted database paths are configured correctly.

## Versioning

Record the following with every analysis:

```text
MATAFILER4 version or git commit
config.txt used
mapping file used
GTDB / GTDB-Tk version
functional database versions
scheduler and cluster partition used
```

This is especially important because MATAFILER4 is under active development and some option names are retained as legacy aliases.
