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

When a submitted job requests node-local scratch and `nodeTmpDir` is configured,
the submission wrapper creates a job-specific directory below `nodeTmpDir`,
exports it as `TMPDIR`, and runs the job from that directory. Generated submission
scripts and scheduler logs remain in their persistent `LOGandSUB` locations.

For Slurm, every job first starts in the persistent directory containing its generated
`LOGandSUB` script and receives `TMPDIR=/tmp`. This prevents dependent jobs from
inheriting a vanished parent-job node-scratch path. When node-local scratch is
requested, the wrapper then creates and enters `nodeTmpDir/matafiler4.$SLURM_JOB_ID`
before running the payload.

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
-schedulerCapacityCheckJobs 10
-killDepNever 1
```

`-maxConcurrentJobs` counts all of the user's running and pending Slurm jobs,
not only dependency-pending work. MATAFILER4 adds each accepted submission to
the cached count and refreshes the exact scheduler count after
`-schedulerCapacityCheckJobs` submissions (default 10), or sooner when the
conservative count reaches the cap. In `-loopTillComplete` mode, a full
queue defers only new submissions; sample completion and cleanup inspection
continues, followed by a bounded scheduler polling delay and a retry of the same
range. `-killDepNever 1` can remove Slurm jobs
stuck in dependency states, but use it only when this matches your local scheduler policy.

ENA and SRA acquisition jobs use the dedicated `downloadQueue` Slurm partition.
Its default is set in the local `Mods/MATAFILERcfg.txt`, which the installer
creates from `Mods/config.old`; the supplied value is `nbi-download`. To
override it, edit the selected site `config.txt`:

```text
downloadQueue	my-network-partition
```

An explicitly empty `downloadQueue` value falls back to `mediumQueue`. This
lets installations without a separate network partition retain the ordinary
default. The queue choice is one-shot and applies only to the archive acquisition
job; dependent filtering, assembly, and mapping jobs use their normal queues.

`netQueue` remains available for other code paths that explicitly request a
general network-enabled queue. It likewise falls back to `mediumQueue` when
empty and uses the configured `longTime` wall time.

`maxMF4mem` sets the maximum GiB requested by automatic OOM retries in the
strain workflows. The supplied default is `512`; `-treeOOMMaxMemGB` remains an
explicit per-run override for tree construction.

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
