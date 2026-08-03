<!-- Documentation navigation -->
[Home](../README.md) | [Quick start](quickstart.md) | [Installation](install.md) | [Configuration](configuration.md) | [Mapping files](mapping_files.md) | [Workflows](common_workflows.md) | [Outputs](outputs.md) | [Flag reference](flag_reference.md) | [FAQ](FAQ.md) | [Glossary](glossary.md)

---

# Quick start

This page gives the shortest practical route to a first MATAFILER4 run. For detailed installation and HPC configuration, see [Installation](install.md) and [Configuration](configuration.md).

## 1. Install and activate

```bash
git clone https://github.com/hildebra/MATAFILER4.git
cd MATAFILER4
bash helpers/install/installer.sh
micromamba activate MF4
./MATAF4.pl -checkInstall
```

## 2. Prepare `config.txt`

Check at least these entries:

```text
MFLRDir        /path/to/MATAFILER4/
DBDir          /path/to/database_dir/
globalTmpDir   /path/to/shared/scratch/
nodeTmpDir     /path/to/node/local/tmp/
```

`globalTmpDir` must be visible to all compute nodes. `nodeTmpDir` should point to local node storage if your cluster provides it.

## 3. Prepare a mapping file

A minimal assembly mapping file looks like this:

```text
#SmplID	Path	AssmblGrps
#OutPath	/path/to/output_base
#RunID	my_matafiler_run
#DirPath	/path/to/raw_reads
Sample01	Sample01_dir	Sample01
Sample02	Sample02_dir	Sample02
```

The file must be tab-delimited. Keep `#SmplID` values short, unique and free of special characters. Use either `Path` or `SmplPrefix` to locate a sample's primary reads; do not populate both for the same sample. See [Mapping files](mapping_files.md), including its section on primary BAM input, for the full format.

## 4. Dry-run one sample

```bash
MAP=/path/to/mapping_file.map
perl $MF4DIR/MATAF4.pl   -map "$MAP"   -assembleMG 2   -assemblCores 12   -assemblMemory 100   -mapReadsOntoAssembly 1   -requireInput 1   -submit 0   -from 0 -to 1
```

Use `-submit 0` first to check paths, scheduler setup and command construction without submitting the full analysis. `-requireInput 1` makes this initial run stop if a sample directory is missing or no input files match; the default is `0`, which allows continuation when original reads have intentionally been removed after earlier processing.

## 5. Submit

```bash
perl $MF4DIR/MATAF4.pl   -map "$MAP"   -assembleMG 2   -assemblCores 12   -assemblMemory 100   -mapReadsOntoAssembly 1   -requireInput 1   -submit 1   -from 0 -to 1
```

After submitted jobs finish, rerun the same command. MATAFILER4 uses completion marker files and should pick up unfinished or missing steps rather than restarting completed work.

Normal runs use the sample completion markers directly and do not perform the
additional full-workflow inspection. Add `-autoStatePlan 1` to inspect and plan
current state, repair narrowly safe partial outputs, and save JSON audit files
under `#OutPath/#RunID/LOGandSUB/workflow/`. With `-loopTillComplete`, an enabled
preflight repeats after each wait and before the next pass.

Input directories are checked lazily as unfinished samples reach input staging,
so startup no longer scans every mapped sample. Use `-precheckInputDirs 1` only
when an eager all-sample scan is desired; combine it with `-requireInput 1` to
abort on missing or invalid inputs.

Samples are revisited in a rolling window. Its start advances only across the
contiguous prefix whose
requested outputs are complete and whose temporary directories have been
cleaned, keeping unfinished samples in rotation while admitting new ones.
A sample that passes every requested output and cleanup check is closed with
`<SmplID>/MATAFILER.sample.complete.json`. Matching later visits read only that
record and report the sample complete. A changed workflow request, redo option,
or downstream ContigStats, SNP, or binner work removes the sentinel and restores
full checking. The loop still performs one final full-range verification pass
before publishing statistics from the sentinels. That first full-range pass
starts immediately; another full-range pass cannot start until all pending and
running jobs submitted by this invocation have finished. Each visit submits only the next ready producer wave, so a failed
input filter is retried before assembly, mapping, binning, or consensus jobs are
created. Scheduler state is collected once per pass for all sample locks,
and `-maxConcurrentJobs` counts both pending and running user jobs before every
Slurm submission. A full queue defers new submissions without skipping sample
completion or cleanup checks; the same rolling range is retried after the
configured scheduler polling interval.

## 6. Check completion

A successful assembly-dependent run should normally create:

```text
#OutPath/#RunID/metagStats.txt
#OutPath/#RunID/metagStatsReport.html
#OutPath/#RunID/<SmplID>/MATAFILER.sample.complete.json
#OutPath/#RunID/<SmplID>/assemblies/
#OutPath/#RunID/<SmplID>/assemblies/metag/genePred/
#OutPath/#RunID/<SmplID>/assemblies/metag/ContigStats/
#OutPath/#RunID/<SmplID>/mapping/
```

For downstream interpretation, continue with [Outputs](outputs.md). For a gene catalog, continue with [Common workflows](common_workflows.md#gene-catalog-creation).


## Source-validated option names

Use [Flag reference](flag_reference.md) when adapting commands. In particular, `MATAF4.pl` accepts `-Binner`, but `geneCat.pl` uses `-binSpeciesMG`; `buildTree5.pl` uses `-aa` for amino-acid FASTA input.
