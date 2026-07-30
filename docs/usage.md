<!-- Documentation navigation -->
[Home](../README.md) | [Quick start](quickstart.md) | [Installation](install.md) | [Configuration](configuration.md) | [Mapping files](mapping_files.md) | [Workflows](common_workflows.md) | [Outputs](outputs.md) | [Flag reference](flag_reference.md) | [FAQ](FAQ.md) | [Glossary](glossary.md)

---

# Running MATAFILER4


## Optional state inspection workflow

For a normal run, no extra planning command or manual approval cycle is needed:

```bash
perl MATAF4.pl -map project.map [normal workflow flags] -submit 1
```

Normal runs now enter the ordinary pipeline logic directly. This avoids a
full-workflow metadata scan before the per-sample completion checks. To opt in
to the state planner and automatic safe repairs, add `-autoStatePlan 1`.

When enabled, the preflight:

1. Inspect files, completion markers, samples, and assembly groups.
2. Build a dependency-aware repair/submission plan.
3. Automatically invalidate narrowly scoped partial mapping, coverage, and
   hybrid preassembly-package outputs.
4. Reinspect the repaired state, then let the existing submission engine pick
   up unfinished work.

Users can still rerun the same command and MATAFILER4 resumes incomplete
samples from its ordinary completion markers. With `-autoStatePlan 1`,
`-submit 0` previews safe repairs without deleting their targets, and
`-autoRepairState 0` keeps inspection and planning but disables repairs.

Each preflight writes an audit snapshot beneath
`#OutPath/#RunID/LOGandSUB/workflow/`. Files are numbered by iteration, for
example `state.iteration-000.json` and `plan.iteration-000.json`.

When both `-loopTillComplete` and `-autoStatePlan 1` are enabled, the first
preflight runs before submission. Normally, at each loop boundary MATAFILER4
waits for the jobs submitted by the current pass, reinspects completed hybrid
packages and assembly-group outputs, applies safe repairs, and only then starts
the next pass. A sparse pass is handled differently: when it submits between
one and `floor(window-size/4)` jobs (with a minimum threshold of one), the next
sample block is added before the wait. The next block is also always added on
the final allowed pass of the current block, even when that pass is busy. Only
one block can be added per pass. This overlaps long-running tail work with new
work while retaining bounded admission. The expanded range gets a fresh pass
budget.
Completed members of a hybrid assembly group are retained while missing members
are resubmitted; final group assembly remains dependent on all required
preassembly packages.
As fully complete samples have their temporary data removed, the range start
advances across that continuous prefix and later samples replenish the rolling
range. A completed sample is revisited through an ordered fast probe (mapping,
depth, assembly, binning and requested variant outputs); a missing priority
output restores the full inspection path. When rolling processing reaches its
normal end, MATAFILER4 makes one final full-range verification pass. Statistics
are collected only after that pass finds neither newly submitted work nor
active sample locks.

With the default `-rmSmplLocks 0`, a pass may submit nothing because active
samples are still locked. MATAFILER4 checks the scheduler in that case and
reruns the current expanded range when at least one job remains active and the
count is either below 3 or strictly below 1% of the range's sample count. These
rescans consume the configured pass allowance; one additional final rescan is
permitted, after which normal overlap or window advancement resumes. This guard
does not run for dry runs, lock-removal mode, or when no scheduler jobs remain.

At the beginning of each loop pass, sample lock job IDs are collected first.
MATAFILER4 uses one `squeue` snapshot and bounded, multi-job `sacct` calls for
tracked IDs no longer present in `squeue`; all samples and dependency checks in
that pass reuse the result. `-maxConcurrentJobs` counts both running and pending
jobs and is enforced immediately before each scheduler submission. During
`-loopTillComplete`, reaching the cap defers submissions rather than blocking
inside one sample: output and cleanup checks continue for the remaining samples,
then the same rolling range is retried after `-schedulerPollSeconds`. Outside
loop mode, the cap retains its blocking behaviour.

Group-wide invalidation is intentionally not classified as an automatic safe
repair. An exact assembly-group membership change remains blocked unless the
existing `-OKtoRWassGrps 1` authorization is supplied.

### Optional expert diagnostics

The stabilization model still separates three concerns for debugging:

1. **Inspect:** read files and completion markers and report their current state.
2. **Plan:** convert that snapshot into reviewable repair and submission actions
   with explicit dependencies.
3. **Execute:** safely repair selected state internally, then use MATAFILER4's
   established submission functions.

The explicit inspect and plan commands are optional expert tools. They never
initialize the scheduler, submit jobs, create scratch directories, or repair
outputs. Plan documents are diagnostic rather than executable shell scripts;
`execution_supported` therefore remains `0` in plan JSON.

### Inspect state

Inspect existing sample, mapping, assembly, coverage, and assembly-group state
without creating scratch directories, deleting outputs, repairing files, or
submitting jobs:

```bash
perl MATAF4.pl -map project.map -inspectState 1
```

The report is JSON on standard output. To write it to a selected file instead:

```bash
perl MATAF4.pl -map project.map -inspectState 1 -stateReport state.json
```

Incomplete combinations, such as a completion stone without its expected
artifact or an assembly-group membership mismatch, are reported as issues. The
inspection command never applies repairs.

### Build a repair/submission plan

Generate an explicit dependency-ordered repair and submission plan from the same
inspection snapshot:

```bash
perl MATAF4.pl -map project.map -planState 1
```

Pass the same workflow flags that would be used for execution, such as
`-assembleMG 5` for hybrid mode. The plan only proposes stages that were
requested by those flags.

The plan is JSON on standard output. It embeds the source inspection report and
lists repairs, confirmation gates, submissions, expected outputs, and
`depends_on` action IDs. Hybrid mode (`-assembleMG 5`) is represented as
per-sample preassembly packages followed by the final assembly-group job;
downstream mappings and contig statistics depend on that final job.

To persist both documents explicitly:

```bash
perl MATAF4.pl -map project.map -planState 1 \
  -stateReport state.json -planReport plan.json
```

Plan generation remains read-only. It does not interpret the plan as shell
commands, delete any target, initialize the scheduler, or submit jobs.
Group-wide invalidation actions retain the `OKtoRWassGrps` authorization gate.

Each action has a stable `id`, `kind`, `operation`, `scope`, `reason_codes`, and
`depends_on` list. Repair actions additionally identify their targets,
`automatic_targets`, automatic policy, and required authorization; submission
actions list their expected outputs. The action graph is validated for missing
dependencies and cycles before it is emitted.


This page describes running behaviour and core concepts. New users should usually read [Quick start](quickstart.md), [Mapping files](mapping_files.md) and [Common workflows](common_workflows.md) before using the full [Flag reference](flag_reference.md).

## Running MATAFILER4

MATAFILER4 is designed for Linux HPC environments and for processing thousands of metagenomes. It therefore relies on job schedulers (slurm, SGE and LSF are supported) and multiple safeguards to resume failed jobs. Please see examples below for specific runs.

### Temporary and output files

The output path for non-temporary files, such as assemblies, bins, gene predictions and abundance tables, is defined in each mapping file separately, composed of the arguments "#OutPath" and "#RunID". Final run-level output is stored in "#OutPath/#RunID/", here each sample will have its own folder, and within this folder assemblies, gene predictions (assembly dir), mapping reads to the assemblies (mapping dir) and a detailed log of the steps run (LOGandSUB dir), will be stored.

Because the pipeline is expected to run on a compute cluster, temporary directories are critical for performance and for file exchange between compute nodes that are usually physically separated clusters.
The pipeline expects a path to storage that is globally available on all nodes and a temporary directory that is locally available on each node. These are configured with `globalTmpDir` and `nodeTmpDir` in `config.txt`. For a file-by-file description of final outputs, see [MATAFILER4 outputs](outputs.md). 

### Mapping file for MATAFILER4


The most important input is a mapping file that describes your samples and raw-read locations. See 'examples' dir for some map examples (also explaining how to do compound assemblies, compound mapping). These column names (headers) are reserved key words in the mapping file (other columns can be eg. metadata per sample etc):
- **#SmplID** [STRING] MATAFILER4 maps always need to have the first column names *#SmplID*. The string in this column will be used in all subsequent analyses, intermediate files, sequence heads etc to uniquely identify samples, therefore choose with extreme care! Good practice would be to include some basic information about the sample in the SMPLID, but should be as short and descriptive as possible. *DO NOT USE SPECIAL CHARACTERS IN THE SMPLID, keep it basic*!  
- **Path** [STRING] - is the relative path to primary read files for each sample (see #DirPath, this needs to be set to the absolute path). FASTQ files are selected with the `-inputFQregex*` options. Unpaired primary BAM reads, such as PacBio reads, can instead be selected with `-inputBAMregex '.*\.bam$'` and are converted to FASTQ with `samtools fastq`; see [Mapping files](mapping_files.md#using-bam-files-as-primary-input). 1 or 2 in paired FASTQ file names indicates the first or second read. E.g. al0-0_12s005629-2-1_lane3.2.fq.gz is the second read, here the pipeline expects to have al0-0_12s005629-2-1_lane3.1.fq.gz in the same dir.
Further, you can add the following specifics for each single sample:   
- **AssmblGrps** [STRING] - set this to a number or string. all samples with the same tag will be assembled together (e.g. samples from the same patient at different time points).  
- **MapGrps** [STRING] - set a tag here as in AssmblGrps. All reads from these samples will be thrown together, when mapping against target sequences (only works with option "map2tar" and "map2DB").
- **SeqTech** ['ONT', 'PB', 'ill', 'miSeq', 'hiSeq', 'GAII', 'GAII_solexa', 'proto', '454', 'AVITI', 'SLR'] - Sequencing technology used in sample: 3rd gen: Oxford Nanopore (`ONT`) and PacBio (`PB`); 2nd gen: Illumina short reads (`ill` and its subtechnologies `miSeq`, `hiSeq`, `GAII`, `GAII_solexa`), Ion Proton (`proto`), 454 (`454`), Elements AVITI (`AVITI`), or synthetic long reads (`SLR`). `AVITI` selects the AVITI-specific SDM filtering configuration.
- **ReadLength** - Expected read length in sample. Is usually automatically determined, use with caution!
- **EstCoverage** [0/1] - (Deprecated!!) Used to indicate if the avg coverage of genomes should be estimated in sample.
- **SupportReads** [tag:path] - Additional reads created with a different sequencing technology. E.g. miSeq (`miSeq:/path/to/file`), mate-pair (`mate:/path/to/file`) or PacBio (`PB:/path/to/file.bam`). For several files of the same technology, use one tag and comma-separated paths: `PB:/path/to/pb1.bam,/path/to/pb2.bam`.
- **ExcludeAssembly** [0/1] - Exclude sample from assemblies?
- **cut5PR1** [INT] - remove the first nts (from 5') on read 1
- **cut5PR2** [INT] - remove the first nts (from 5') on read 2
- **firstXreadsRd** [INT] - only read the first X reads (paired reads count as 2) for that sample.
- **firstXreadsWr** [INT] - only write the first X reads (paired reads count as 2) for that sample.
 
 The following tags can be added to a new line (ie row) in the map. Tag is followed by tab delimiter and specific input.

#### Required map tags
- **#OutPath**	[Path] Where to write the output (can be massive, make sure you have enough space)
- **#RunID**	[string] The directory below OutPath, where results are stored. Also serves as global identifier for this run
- **#DirPath**	[Path] Base directory where subdir with the fastqs can be found. You can insert this on several lines, if the base path changes for all samples afterwards.

#### Optional map tags

- **#NodeTmpDir**	[Path] temporary dir only accessible within each compute cluster node, overrides **nodeTmpDir** definition in config file
- **#GlobalTmpDir**	[Path] temporary dir (scratch) accessible from all compute nodes, overrides **globalTmpDir** definition in config file
- **#mocatFiltPath**	If for some reason you are forced to use mocat filtered fastqs and not the original, unfiltered files (strongly recommended), than you can indicate in which subdir these mocat files can be found
- **#RelaxSMPLID**	[TRUE/FALSE] 	Use FALSE to deactivate basic checks if the #SmplID adheres to MATAFILER4 formats. Caution: use on your own risk!
- **#WARNING**	[OFF/ON]	If **OFF** MATAFILER4 won't stop when an error is encountered in the map. Caution: use on your own risk!

After this follow the sample IDs and the relative path, where to find the input fastqs.  
See _examples/example_map_assemblies.map_ for a very complicated mapping file with several source dirs.

#### Example mapping file

```{sh}
#SmplID	Path	SmplPrefix	AssmblGrps
#OutPath	/hpc-home/path/to/your/results/folder
#RunID	NAMEofresultsFOLDER (#MATAFILER4 will make this folder with this name by itself)
#DESCRIPTION (#not important, but you can mark what is the run is about)
#DirPath	/path/to/folder/with/raw/reads
Mouse11t0		PID_C11T0_	M11
Mouse11t1		PID_C11T1_	M11
Mouse12t0		PID_C12T0_	M12
Mouse12t1		PID_C12T1_	M12
Mouse14t0		PID_C14T0_	M14
Mouse14t1		PID_C14T1_	M14
Mouse15t0		PID_PD11T0_	M15
#DirPath	/path/to/another/folder/with/more/raw/reads
Mouse15t1	SubDir1		M15
Mouse16t0	SubDir2		M16
Mouse16t1	SubDir3		M16
```		

#### Tips and recommendations for creating mapping files

- It is recommended to create the mapping file in **Excel** and copy-paste it in a **.map** text file afterwards (will be tab-delimited by default, the expected MATAFILER4 format). You can use functions like "=VLOOKUP()" to match sample IDs across different tables. 

- The **#SmplID** column determines the name of a sample all the way throughout the pipeline! Be very careful what ID you choose, as this will impact the sample names you'll have to deal with later, choose something a) short and b) descriptive. Avoid c) special characters (_|$%~\`\*& etc) in the SmplID!

- **AssmblGrps**: Assembly groups are useful for assembling samples from e.g. a time series together, giving a better assembly usually. Choose the name of an assembly group a) unique b) short and descriptive and c) avoid special chars (\[\]{}_|$%~\`\*& etc)!

- If using **assembly groups**, try to keep samples from the same assembly group as a block. MATAFILER4 can also deal with these assembly groups distributed across the map, but in terms of job submission strategy it's best to have these samples next to each other in the map (and also for you organizing your experiment).

- <ins>**Loading and saving a mapping file into R will likely lead to problems!**</ins> This is because the #DirPath tag sets the path for all samples underneath. Loading this into R will often skip the #DirPath line or reorder the samples, so saving this again will lead to wrong paths being set!







## Additional usage scenarios


MATAFILER4 can be used for a bulk of tasks not directly related to initial assembly, profiling or MAGs, but often related to postprocessing these. Two usage scenarios (map2tar and building phylogenies) are listed below.


### map2tar mode

- This mode maps the reads to `reference` genomes, e.g. from a mock community. This allows a fast profiling for specific purposes. The mode switches off assembly-based processes.

1. create mapping file `path/to/mapping_file.map`

```{sh}
#SmplID	SmplPrefix	AssmblGrps
#OutPath	/path/to/results/profilingMF
#RunID	name_of_the_run_dir
#DESCRIPTION
#DirPath	/path/to/raw/reads
BERG100	BERG100	BERG100
BERG100w	BERG100w	BERG100w
BERG10	BERG10	BERG10
BERG10w	BERG10w	BERG10w
BERGmock	BERGmock	BERGmock
```

2. make script with first command `run_profiling.sh`

```{sh}
MAP="/path/to/mapping_file.map"
perl $MF4DIR/MATAF4.pl map2tar \
	-map $MAP  -ref '/path/to/reference/mock_community/*.fasta' -filterHumanRds 0 -mappingCores 12 -mapperFilterIll '0.02 0.75 00'  -redo2ndmap 0 -mappingMem 15 -submit 1 -competitive2ndmap -1 -decoyMapping 0
```

- Explanation: ref are the .fasta formatted reference genomes you want to map your metagenomic reads to, metagenomic reads are defined in the map, as in other runs. -mapperFilterIll defines how the mapped reads will be quality filtered. -competitive2ndmap defines if reads will be mapped against all references at once (competitive) or separately against each single reference. -decoyMapping determines if an already created read assembly will be used to "decoy" map reads against (useful if you suspect that most reads aligning to your reference would be false positive hits).

- run with `bash run_profiling.sh`

- output files will contain coverage per window, contig, etc. which can be used for plotting.

### Building phylogenetic trees with MATAFILER4

1. create a file `phyloScript.pl`

```{sh}
#!/usr/bin/perl
use strict; use warnings;

my $bts = "/path/to/MATAFILER4/secScripts/phylo/buildTree5.pl";
my $inD = "/path/to/input/dir/phylo/";
my $outD = $inD."/bts/";
my $tempD = "/path/to/scratch/dir/treetest/";
my $numCore = 8;

my $cmd = "$bts -genoInD $inD -outD  $outD -tmpD $tempD -runIQtree 1 -iqFast 1 -AAtree 1 -cores $numCore -wildcardflag '/*.f*a' -continue 1 \n";

print $cmd;
system $cmd;
```

Explanation: $inD is an input dir with complete genomes, the script will extract FGMs and build tree between genomes. `-AAtree 1` tells the script to use AA MSAs to build the phylogeny via iqTree. `-wildcardflag '/*.f*a'` tells the script how to look for reference genomes in $inD. 

2. run the script `phyloScript.pl`

- Run with `perl /path/to/phyloScript.pl` together with submission script on cluster.

3. you can do many additional phylogeny / popgen related analysis with the `buildTree5.pl` script ([see `buildTree5.pl` flags](flag_reference.md#buildtree5pl))
