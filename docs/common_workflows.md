<!-- Documentation navigation -->
[Home](../README.md) | [Quick start](quickstart.md) | [Installation](install.md) | [Configuration](configuration.md) | [Mapping files](mapping_files.md) | [Workflows](common_workflows.md) | [Examples](examples.md) | [Profiling tutorial](profiling_tutorial.md) | [Strain-within guide](strainwithin.md) | [Outputs](outputs.md) | [Flag reference](flag_reference.md) | [FAQ](FAQ.md) | [Glossary](glossary.md)

---

# Common workflows

This page collects practical command patterns. For exact current option names and aliases, see [Flag reference](flag_reference.md).

## Assembly-dependent metagenomic assembly + gene catalog

This mode is appropriate for metagenomes where assemblies are expected to capture a substantial fraction of reads.

### 1. Create a mapping file

See [Mapping files](mapping_files.md).

### 2. Create a run script

```bash
#!/bin/bash
#SBATCH -J SUB_MF
#SBATCH -N 1 --cpus-per-task=1 --mem=10024 --export=ALL
#SBATCH -o /path/to/run_mf4_mhit.mfc.otxt
#SBATCH -e /path/to/run_mf4_mhit.mfc.etxt
#SBATCH -p "ei-long,qib-long"

set -e
ulimit -c 0

MAP=/path/to/your/mapping/file/FILE.map

perl "$MF4DIR/MATAF4.pl"   -map "$MAP"   -assembleMG 2   -assemblCores 12   -assemblyKmers "25,43,67,87,111,131"   -assemblMemory 100   -mapReadsOntoAssembly 1   -kmerPerGene 0   -filterHostRds 1   -filterHostKrak2DB /path/to/kraken2/host_db/   -mappingMem 5   -profileMOTU2 0   -profileMetaphlan 1   -Binner 2   -maxConcurrentJobs 600   -from 0 -to 1   -submit 0   -getAssemblConsSNP 0
```

Notes:

- Start with `-submit 0`; switch to `-submit 1` only after checking the dry-run.
- `-spadesCores`, `-spadesKmers`, `-spadesMemory`, `-MetaBat2` and `-binSpeciesMG` are still accepted as aliases, but the clearer names are `-assemblCores`, `-assemblyKmers`, `-assemblMemory` and `-Binner`.
- Prefer `-profileMetaphlan`; `-profileMetaphlan3` remains accepted as a compatibility alias.

### 3. Run and rerun

```bash
sbatch run_mf4_mhit.mfc
```

Wait until submitted jobs finish, then rerun the same script. MATAFILER4 uses marker files to identify completed steps and continues missing work.

### 4. Expected output after MATAF4.pl completion

```text
#OutPath/#RunID/metagStats.txt
#OutPath/#RunID/metagStatsReport.html
#OutPath/#RunID/<SmplID>/assemblies/
#OutPath/#RunID/<SmplID>/assemblies/metag/genePred/
#OutPath/#RunID/<SmplID>/assemblies/metag/ContigStats/
#OutPath/#RunID/<SmplID>/mapping/
```

## Gene catalog creation

After sample assemblies and mapping have completed, run `secScripts/geneCat.pl` using the generated or adapted `GeneCat.sh` script.

```bash
perl "$MF4DIR/secScripts/geneCat.pl"   -map /path/to/mapping_file.map   -GCd /path/to/results/genecat   -mem 200   -cores 24   -clusterID 95   -doStrains 0   -continue 1   -binSpeciesMG 2   -useCheckM1 0   -useCheckM2 1   -MGset GTDB
```

Note: `geneCat.pl` uses `-binSpeciesMG` for binning/MGS selection. `-Binner` is only a `MATAF4.pl` alias and is not parsed by `geneCat.pl`.

Expected outputs include:

```text
<GCd>/compl.incompl.95.fna
<GCd>/compl.incompl.95.prot.faa
<GCd>/Matrix.mat.gz
<GCd>/Matrix.mat.scaled.gz
<GCd>/Anno/
<GCd>/Binning/ or <GCd>/Bin_SB/  # if MAG/MGS processing was enabled
```

## Explicit MGS.pl run

`geneCat.pl` can orchestrate MGS/MAG processing, but `MGS.pl` can also be run explicitly when needed.

```bash
perl "$MF4DIR/secScripts/MGS.pl"   -GCd /path/to/results/genecat   -tmp /path/to/scratch/MGS   -smallCores 8   -bottleneckCores 24   -mem 150   -binSpeciesMG 2   -useCheckM2 1   -useCheckM1 0   -MGset GTDB   -submit 1
```

Use `-outD` only if you do not want the default output directory below the gene catalog directory. `-MGset` must be either `GTDB` or `FMG`.


## Assembly-independent profiling

For database setup and individual RiboFind, functional, MetaPhlAn 4, mOTUs 4,
and Protal commands, see the [read-based profiling tutorial](profiling_tutorial.md).

This mode is useful for highly complex communities where assembly is not expected to be productive.

```bash
MAP=/path/to/dir/mapping_file.map

perl "$MF4DIR/MATAF4.pl"   -map "$MAP"   -inputFQregex1 '.*_1\.f[^\.]*q\.gz$'   -inputFQregex2 '.*_2\.f[^\.]*q\.gz$'   -mergeReads 0   -profileFunct 1   -reParseFunct 0   -reProfileFunct 0   -DiaDBs KGM,NOG,CZy   -DiaCores 8   -maxConcurrentJobs 300   -profileRibosome 1   -reProfileRibosome 0   -profileMOTU2 0   -profileMetaphlan 1   -filterHostRds 0   -inputReadLength 150   -assembleMG 0   -submit 1   -from 0 -to 40
```

Expected outputs are primarily in run-level profiling folders such as `pseudoGC/`, plus summary files such as `metagStats.txt` and `metagStatsReport.html`.

### Protal profiling from raw reads

For a profile per sample and the merged abundance table:

```bash
perl "$MF4DIR/MATAF4.pl" \
  -map "$MAP" \
  -assembleMG 0 \
  -profileProtal 1 \
  -ProtalCores 4 \
  -ProtalMem 100 \
  -submit 1
```

`-profileProtal 2` does the same, also keeps each sample's SAM and builds the strain
MSAs of the cohort from them:

```bash
perl "$MF4DIR/MATAF4.pl" \
  -map "$MAP" \
  -assembleMG 0 \
  -profileProtal 2 \
  -ProtalCores 4 \
  -ProtalMem 100 \
  -submit 1
```

Protal reads each sample's source files (yours, or the downloaded ones) rather than
staged or SDM-cleaned reads, and keeps the sample's profile
(`<sample>/Tax/Protal/profiles/<sample>.profile`). Its SAM is kept
(`<sample>/Tax/Protal/<sample>.sam.zst`) only with `-profileProtal 2`, which needs it
for the strain MSAs; with `1` it is deleted after profiling. Switching from `1` to `2`
therefore aligns the samples again. How a sample is aligned depends on where its
reads come from:

- **Downloaded samples** (`ENAdownload`/`SRAdownload` in the map) get a Protal job of
  their own as soon as their download is ready, also when they belong to an assembly
  group; the sample's scratch, and with it the download, is cleaned up after that
  job, so a cohort is never held on disk at once. The runs of a sample are profiled
  together.
- **Local samples** are aligned together, in one batch map job per pass, which loads
  Protal's index once. `-protalBatchSize N` splits a pass's local samples into jobs of
  at most `N` samples. A local sample whose reads need joining (several libraries) or
  decompressing (bzip2, xz) gets a job of its own.

Each sample is profiled with the reads of one kind: paired-end, else single-end,
else PacBio or ONT reads, using every library of that kind. With the default
`-protalIgnoreErrors 1`, other kinds of reads and BAM input are left out with a
warning, and a sample with nothing usable is skipped (`<sample>.Protal.skip`).
`-protalIgnoreErrors 0` stops instead. A failed Protal job is never hidden: in a
batch, the samples Protal finished are kept and the job fails for the others, which a
later pass aligns again.

The profiles of the selected `-from`/`-to` range are merged into
`pseudoGC/protal/Protal.abundance.tsv`; a sample with neither a profile nor a
submitted job defers the merge, preventing a partial table. In mode `2`, a run over
the complete mapped cohort (no `-from`/`-to`) then has Protal profile the kept SAMs
again, without aligning or loading its index, and write an MSA per species found in
two or more samples to `pseudoGC/protal/strains/<cohort_signature>/`
(`strains/current` points at the latest cohort). A run over part of the cohort
builds no MSAs. The same can be done by hand on any set of kept SAMs, for example
after adding samples: `protal --map` with a `SAM` column naming them, and no
`--no_strains`.

## Hybrid assemblies

Use hybrid workflows when short reads are supported by ONT, PacBio or mate-pair reads. In the mapping file, provide support reads through `SupportReads`, for example:

```text
Sample01	Sample01_dir	Sample01_	S01	PB:/path/to/pacbio_reads.bam
```

Relevant flags include:

```bash
-assembleMG 5
-mapSupportReadsOntoAssembly 1
-mapSaveCRAM 1
```

`-mapSaveCRAM` is kept enabled for hybrid assemblies because downstream support-read steps may need the mappings.

## map2tar, map2DB and map2GC reference mapping

The `map2tar` mode maps reads to a supplied reference FASTA.

```bash
perl "$MF4DIR/MATAF4.pl" map2tar   -map "$MAP"   -ref /path/to/reference.fa   -mapnms my_reference   -competitive2ndmap 1   -decoyMapping 1   -submit 1
```

Important options are `-ref`, `-mapnms`, `-mapUnmapped`, `-competitive2ndmap`, `-decoyMapping`, `-redo2ndmap` and `-mappingCoverage`.
