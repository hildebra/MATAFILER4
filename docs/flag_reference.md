<!-- Documentation navigation -->
[Home](../README.md) | [Quick start](quickstart.md) | [Installation](install.md) | [Configuration](configuration.md) | [Mapping files](mapping_files.md) | [Workflows](common_workflows.md) | [Examples](examples.md) | [Profiling tutorial](profiling_tutorial.md) | [Strain-within guide](strainwithin.md) | [Outputs](outputs.md) | [Flag reference](flag_reference.md) | [FAQ](FAQ.md) | [Glossary](glossary.md)

---

# Flag reference

Quick links: [`MATAF4.pl`](#mataf4pl) | [`geneCat.pl`](#genecatpl) | [`MGS.pl`](#mgspl) | [`strain_within.pl`](#strain_withinpl) | [`strain_within_2.2.pl`](#strain_within_22pl) | [`buildTree5.pl`](#buildtree5pl)

This page is validated against the repository Perl source files for `MATAF4.pl`, `geneCat.pl`, `MGS.pl`, `strain_within.pl`, `strain_within_2.2.pl` and `buildTree5.pl`.

| Script | Version in referenced source | Role |
|---|---:|---|
| `MATAF4.pl` | `4.46` | Main sample-level pipeline: read detection, preprocessing, host filtering, assembly, mapping, binning, SNP/SV calling and read-based profiling. |
| `geneCat.pl` | `0.58` | Gene catalog construction and downstream gene-catalog annotation/MGS orchestration. |
| `MGS.pl` | `0.55` | MGS/MAG dereplication, abundance/taxonomy and optional strain workflow orchestration. |
| `strain_within.pl` | `1.53` | Within-MGS locus extraction, quality control, tree preparation/submission and downstream hand-off. |
| `strain_within_2.2.pl` | `0.49` | Within-MGS tree postprocessing, strain statistics and optional population-genetic analysis. |
| `buildTree5.pl` | `5.83` | Phylogenetic tree construction and related MSA/population-genetic analyses. |

## How to read the tables

This page is the single source of the option lists. Every entry point renders its own `-help`
from the section below that carries its name, through `Mods::FlagReference`, so a table edited
here changes the command-line help as well. Editing a table therefore needs no matching change
in the scripts — but adding or renaming a flag in a script does need a matching row here, or it
becomes invisible to users.

- Each script accepts `-help`, `-h` or `-?` (one or two leading dashes), answered before any site configuration is loaded.
- **Aliases** are equivalent command-line names. Use any listed alias with a leading dash.
- **Type** is the argument type parsed by `Getopt::Long`: integer, float, string, string list or flag.
- **Default** is inferred from source-code assignments where possible. Empty cells mean no simple default could be inferred automatically.
- **Status** marks options that look stable, advanced/internal, deprecated/legacy or experimental/unsupported based on source comments and option names.

## Important validation notes

- `geneCat.pl` does **not** accept `-Binner`; use `-binSpeciesMG` for gene-catalog/MGS binning selection.
- `MATAF4.pl` accepts `-Binner`, `-MetaBat2` and `-binSpeciesMG` as aliases for the sample-level binning option.
- `-profileMetaphlan3` remains a compatibility alias for `-profileMetaphlan`.
- `-profileProtal 1` and `2` both take the first compatible short-read pair and skip samples without one; `-protalIgnoreErrors 0` instead requires exactly one paired-end library and no singleton/BAM input. Mode `1` merges per-sample profiles under `pseudoGC/protal_singular/` (no strain analysis); mode `2` runs one cohort-wide map and publishes the merged table plus strain MSAs under `pseudoGC/protal/`.
- The strain workflow entry points are named `strain_within.pl` and `strain_within_2.2.pl` in the repository; spellings without underscores are not repository filenames.
- `buildTree5.pl` accepts `-aa` for the amino-acid FASTA input; the older comment spelling `-faa` is not a parsed flag.
- `geneCat.pl -MGset` and `MGS.pl -MGset` are constrained in source to `GTDB` or `FMG`.

## MATAF4.pl

Main sample-level pipeline. `MATAF4.pl -help` prints the option tables below, so keep every description to one or two short sentences.

## Base options

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-help`, `-?`, `-h` | flag | `` | stable | Show help. |
| `-checkInstall` | flag | `` | stable | Check that core MATAFILER4 programs/environments are installed. |
| `-map` | string | `""` | stable | Mapping file describing samples and input/output paths. |
| `-config` | string | `""` | stable | Alternative configuration file. |
| `-precheckInputDirs` | integer | `0` | advanced | Pre-scan all mapped sample inputs and cache their sizes. Use `-requireInput 1` to fail on missing input. |
| `-inspectState` | integer | `0` | stable | Emit a read-only JSON snapshot of workflow artifacts and markers. |
| `-planState` | integer | `0` | stable | Emit a read-only, dependency-ordered repair/submission plan from the inspection snapshot. |
| `-stateReport` | string | `""` | stable | Write the inspection JSON to this explicit path. |
| `-planReport` | string | `""` | stable | Write the repair/submission plan JSON to this explicit path. |
| `-autoStatePlan` | integer | `0` | advanced | Run the state preflight before execution and at each loop boundary. Disabled by default. |
| `-autoRepairState` | integer | `1` | stable | Apply only preflight repairs classified as automatically safe when submission is enabled. |

## Flow related

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-redoFails` | integer | `0` | stable | Delete failed sample-local results and rebuild on the next pass; shared assembly-group output is kept. |
| `-redoContigStats` | integer | `0` | stable | runContigStats (coverage per gene, kmers, GC content) will be deleted & started again |
| `-submSystem` | string | `""` | stable | qsub,SGE,bsub,LSF.. by default will try to autodetect |
| `-submit` | integer | `1` | stable | submit any jobs at all? (0= no submission, just for trying if everything is correctly set up) |
| `-from` | integer | `0` | stable | start at which samples from map file? |
| `-to` | integer | `999999999999` | stable | stop at which samples from map file? |
| `-loopTillComplete` | string | `"0"` | advanced | Repeat samples in rolling windows until complete, then run one full verification pass. Syntax `X:Y` sets the pass budget and window size; statistics are collected only after a clean full pass. |
| `-loopTillCompleteActiveJobs` | integer | `3` | advanced | Start the next rolling pass once at most this many window jobs are still running. |
| `-schedulerPollSeconds` | integer | `20` | advanced | Seconds between scheduler queries while `loopTillComplete` waits. Values must be positive. |
| `-schedulerCapacityCheckJobs` | integer | `10` | advanced | Refresh the exact Slurm job count after this many submissions, or sooner at `-maxConcurrentJobs`. |
| `-excludeNodes` | string | `""` | stable | exclude certain nodes? |
| `-maxConcurrentJobs` | integer | `0` | stable | Maximum running plus pending user jobs; 0 = unlimited. Loop runs defer when full, non-loop runs wait. |
| `-killDepNever` | integer | `0` | stable | kill jobs in "Dependency never finished" state? |
| `-requireInput` | integer | `0` | stable | in case input reads are no longer present, 0 will continue pipeline, 1 will abort |
| `-ignoreSmpls` | string | `""` | stable | Comma-separated exact sample IDs to skip (no regex, no prefix match). |
| `-rmSmplLocks` | integer | `0` | stable | Remove existing sample lock files. |
| `-silent` | flag | `0` | stable | Suppress routine progress messages. |
| `-maxUnzpJobs` | integer | `20` | stable | how many unzip jobs to run in parallel (not to overload HPC IO). Default:20 |
| `-skipSmallSmplsMB` | integer | `1` | stable | skip samples with a combined input smaller than this in MB (raw file size, independent of compressed or raw) |
| `-forceWriteStats` | integer | `0` | stable | force (re)writing of the metagStats report and text file |

## File structure

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-rm_tmpdir_reads` | integer | `1` | stable | Default 1, remove tmpdir with reads |
| `-rm_tmpInput` | integer | `1` | stable | remove raw, human / adaptor filtered reads, if sdm clean created? (and not needed any longer) |
| `-reduceScratchUse` | integer | `1` | internal/advanced | remove sample scratch and rebuildable indexes once terminal outputs are published; 0 for debugging |
| `-globalTmpDir` | string | `""` | stable | absolute path to global shared tmp dir (like a scratch dir) |
| `-nodeTmpDir` | string | `""` | stable | absolute path to tmp dir on local HDD of each executing node |
| `-nodeHDDspace` | string | `30` | stable | HDD tmp space to be requested for each node (in Gb). Some systems don't support this |
| `-legacyFolders` | integer | `0` | stable | legacy option, controls if output folders will use the read dir as name (1) or the name in the mapping file (0). Default=0 |

## Preprocessing (cleaning reads etc)

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-usePorechop` | integer | `0` | stable | adapter rm for Nanopore.. should actually be automatically with newer sdm (not implemented) |
| `-inputFQregex1` | string | `'.*1\.f[^\.]*q\.gz$'` | stable | regex for detecting read pair 1 in input fastq files |
| `-inputFQregex2` | string | `'.*2\.f[^\.]*q\.gz$'` | stable | regex for detecting read pair 2 in input fastq files |
| `-inputFQregexSingle` | string | `""` | stable | regex for detecting single end reads in input fastq files |
| `-inputFQregexTrustSingle` | integer | `0` | stable | if grep of files (rawSrchString) has multi assignments, which grep to trust more? |
| `-inputBAMregex` | string | `""` | stable | regex for detecting BAM read files (e.g. `'.*\.bam$'`); matches are treated as unpaired and converted with `samtools fastq`. Empty disables BAM input |
| `-splitFastaInput` | integer | `0` | stable | Enable FASTA input splitting during read staging. |
| `-mergeReads` | integer | `0` | stable | merge read pair 1+2 before assembly etc? (usually doesn't help assembly, but useful for mapping to ref database in some rare instances) |
| `-ProbRdFilter` | integer | `1` | stable | Enable probabilistic SDM read filtering. |
| `-pairedReadInput` | integer | `-1` | stable | determines if read pairs are expected in each in dir |
| `-inputReadLengthSuppl` | integer | `5000` | stable | Default supplementary-read length used for planning. |
| `-filterHostRds`, `-filterHumanRds` | integer | `0` | stable | 0: no, 1: kraken2, 3: hostile. 2 (kraken1) is no longer supported and aborts. |
| `-filterHostKrak2DB` | string | `""` | stable | customize host org to filter (e.g. human, chicken ..) |
| `-filterHostKr2Conf` | string | `0.01` | stable | Kraken2 confidence threshold for host assignment. |
| `-filterHostKr2Quick` | string | `` | stable | Kraken2 quick-mode option for host filtering. |
| `-hostileIndex` | string | `"human-t2t-hla"` | stable | Hostile reference index used for host filtering. |
| `-onlyFilterZip` | integer | `0` | stable | Stage and filter reads, then stop before downstream analyses. |
| `-mocatFiltered` | integer | `0` | stable | Import reads from the configured MOCAT filtered-output layout. |
| `-logQualvsLen` | integer | `0` | stable | sdm log file.. can be quite large; logs qual of read vs read length |

## Sdm related

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-inputReadLength` | integer | `150` | stable | Default primary-read length used for planning. |
| `-gzipSDMout` | integer | `1` | stable | Compress SDM-filtered FASTQ output. |
| `-XfirstReads` | integer | `-1` | stable | Process only the first N reads; -1 processes all reads. |
| `-minReadLength` | integer | `0` | stable | Override the minimum read length accepted by SDM; 0 keeps its configured value. |
| `-maxReadLength` | integer | `0` | stable | Override the maximum read length accepted by SDM; 0 keeps its configured value. |
| `-filterAdapters` | integer | `1` | stable | Enable adapter trimming during SDM filtering. |
| `-customSDMopt` | string | `""` | stable | Use the specified non-empty SDM options file. |
| `-sdmMem` | string | `"15G"` | stable | total mem for sdm job in Gb, default 15 |

## Assembly related

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-spadesCores`, `-assemblCores` | integer | `0` | stable | 0 = autoscale with assembly-group input, 8 cores at <=500 MiB to 48 cores at >=10 GiB. Positive values override |
| `-spadesMemory`, `-assemblMemory` | integer | `-1` | stable | in GB |
| `-spadesKmers`, `-assemblyKmers` | string | `"27,43,67,87,101,127"` | stable | comma delimited list |
| `-reAssembleMG` | integer | `0` | stable | Rebuild an assembly; shared assembly groups additionally require `-OKtoRWassGrps 1`. |
| `-asssemblyHddSpace` | integer | `"-1"` | stable | Assembler scratch request in GB; -1 uses the assembler-specific default. |
| `-assembleMG` | integer | `0` | stable | 1=Spades, 2=MegaHIT, 3= flye, 4=metaMDBG, 5=hybrid ill-PB (megahit, metaMDBG) |
| `-assemblyLongTime` | integer | `0` | stable | Submit assembly jobs with the long-runtime setting. |
| `-assemblyScaffMinSize` | integer | `500` | stable | Minimum retained scaffold length in bases. |

## Binning

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-Binner`, `-MetaBat2`, `-binSpeciesMG` | integer | `0` | stable | 0=no, 1=metaBat2, 2=SemiBin, 3: MetaDecoder, 4: GenomeFace, 5: SCGBinner |
| `-BinnerCores` | integer | `9` | stable | cores used for Binning process (and checkM) |
| `-BinnerMem` | integer | `0` | stable | define binning memory, Gb, 0=auto |
| `-minBinnerAssemblyMB` | float | `2` | stable | Skip binning when the assembly holds fewer million sequence bases; empty binner outputs are published instead. `0` disables the cutoff. |
| `-checkM2` | integer | `1` | stable | Assess recovered bins with CheckM2. |
| `-checkM1` | integer | `0` | stable | Assess recovered bins with CheckM1. |
| `-BinnerScratchTmp` | integer | `0` | internal/advanced | very specific (undocumented) use of scratch instead of nodetmp dir |
| `-redoEmptyBins` | integer | `0` | internal/advanced | debug option; redo bins that are empty (no bin detected). Note: this can sometimes happen for metagenomes |
| `-redoBinning` | integer | `0` | stable | Remove and rebuild requested binning results. |
| `-SB_env` | string | `""` | stable | semiBin environment; if given, will avoid re-training de novo binning model. Default: "" (autotrain). should be #human_gut/dog_gut/ocean/soil/cat_gut/human_oral/mouse_gut/pig_gut/built_environment/wastewater/chicken_caecum/global |

## Gene prediction on assembly

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-predictEukGenes` | integer | `0` | stable | severely limits total predicted gene amount (~25% of total genes) |
| `-kmerPerGene` | integer | `0` | stable | calculate kmer frequencies for each gene instead of per scaffold |
| `-genePredGZenforce` | integer | `1` | stable | Require the canonical compressed gene-prediction outputs. |
| `-rewriteGenePred` | integer | `0` | stable | Remove and rebuild gene predictions and dependent statistics. |

## Mapping

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-mapper` | integer | `-1` | stable | 1=bowtie2, 2=bwa, 3=minimap2, 4=kma, 5=strobealign -1=auto (bowtie2 short, minimap2 long reads), -2=auto(strobealign short, minimap2 long) |
| `-mapUnmapped` | integer | `0` | stable | Pass reads left unmapped by one reference to the next mapping. |
| `-mappingCoverage` | integer | `1` | stable | Calculate per-contig and per-gene coverage for reference mappings. |
| `-mappingMem` | integer | `-1` | stable | total mem for mini2/kma/bwa/bwt2 in GB |
| `-mapSortMem` | integer | `-1` | stable | total mem for samtools sort in GB |
| `-rmDuplicates` | integer | `1` | stable | Remove duplicate alignments before coverage calculation. |
| `-mappingCores` | integer | `8` | stable | CPU cores requested for each mapping job. |
| `-mapperFilterIll` | string | `"0.05 0.75 20 3"` | stable | max NM edit rate, min query coverage, min mapping quality, min clip at both ends (0 disables clipping) |
| `-mapperFilterHybridIll` | string | `"0.03 0.90 40 5"` | advanced | as -mapperFilterIll, but for coverage of hybrid preassemblies (stricter) |
| `-hybridMinMapQ` | integer | `40` | advanced | Minimum mapping quality passed to `samtools depth` for hybrid-preassembly coverage. |
| `-hybridMinBaseQ` | integer | `20` | advanced | Minimum base quality passed to `samtools depth` for hybrid-preassembly coverage. |
| `-breakpointDepth` | float | `0.10` | advanced | Relative coverage threshold used to identify low-depth assembly breakpoints for hybrid read simulation. |
| `-breakpointMinLength` | integer | `100` | advanced | Minimum length of a low-depth region reported as a hybrid-assembly breakpoint. |
| `-breakpointSmoothGap` | integer | `100` | advanced | Maximum gap joined while smoothing adjacent low-depth breakpoint regions. |
| `-breakpointFlankLength` | integer | `500` | advanced | Number of bases inspected on each side of a candidate breakpoint. |
| `-breakpointMinFlankDepth` | float | `1` | advanced | Minimum flank depth required when accepting a candidate breakpoint. |
| `-breakpointMaxFlankFraction` | float | `0.10` | advanced | Maximum low-depth fraction allowed within breakpoint flanks. |
| `-hybridSyntheticMaxDepth` | float | `20` | advanced | Cap on synthetic read depth generated from each hybrid-preassembly package. |
| `-mapperFilterPB` | string | `"0.05 0.5 30 0"` | stable | PacBio alignment filter: maximum edit rate, minimum query coverage, mapping quality and clipping. |
| `-mapperFilterONT` | string | `"0.15 0.5 10 0"` | stable | Nanopore alignment filter: maximum edit rate, minimum query coverage, mapping quality and clipping. |
| `-mapSaveCRAM` | integer | `1` | stable | Retain assembly back-mapping alignments as CRAM files. |

## Mapping related (2) (assembly)

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-remap2assembly`, `-redoMap2assembly`, `-redoMapping` | integer | `0` | stable | Remove and rebuild read mappings to the assembly. |
| `-JGIdepths` | integer | `0` | stable | Generate JGI depth output for assembly mappings. |
| `-mapReadsOntoAssembly` | integer | `1` | stable | map original reads back on assembly, to estimate abundance etc |
| `-mapSupportReadsOntoAssembly` | integer | `1` | stable | Map `SupportReads` onto the assembly and calculate their coverage separately. |
| `-saveReadsNotMap2Assembly` | integer | `0` | stable | Save reads that do not map to the assembly. |

## Map2tar / map2db / map2gc

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-decoyMapping` | integer | `1` | stable | 1: "Decoy mapping": map against reference genome AND against assembly of metagenome (drawing obvious better hits to metagenome, the "decoy") |
| `-competitive2ndmap` | integer | `-1` | stable | -1: combined map, combined reporting; 0: a separate mapping run per reference; 1: competitive; 2: combined map, separate reporting per input genome |
| `-ref` | string | `""` | stable | Comma-separated reference FASTA files for secondary mapping. |
| `-mapperLargeRef` | integer | `0` | stable | use flags in mapper index built for large ref DBs? |
| `-mapnms` | string | `""` | stable | name for this final files |
| `-redo2ndmap` | integer | `0` | stable | Remove and rebuild secondary-reference mappings. |

## Snps

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-get2ndMappingConsSNP` | integer | `0` | stable | SNPs (onto mapping) |
| `-getAssemblConsSNP` | integer | `0` | stable | SNPs (onto self assembly) #calculates consensus SNP of assembly (useful for checking assembly gets consensus and Assmbl_grps) |
| `-getAssemblConsSNPsuppRds` | integer | `0` | stable | same as getAssemblConsSNP, but SNP calling for support reads |
| `-redoAssmblConsSNP` | integer | `0` | stable | Remove and rebuild assembly consensus SNP outputs. |
| `-SNPmem` | integer | `0` | stable | memory per assigned core, in GB |
| `-redoGeneExtrSNP` | integer | `0` | stable | Rebuild gene and protein sequences derived from consensus SNPs. |
| `-SNPjobSsplit` | integer | `0` | stable | parallel jobs per sample; 0 estimates 1-10 jobs from alignment size (capped by -SNPcores) |
| `-SNPminCallQual` | integer | `20` | stable | Minimum variant quality used for consensus calls. |
| `-SNPsaveVCF` | integer | `1` | stable | save vcf of SNP calles? DEfault : 1 |
| `-SNPsaveConsFasta` | integer | `0` | stable | Save consensus fasta from vcf calls? Default: 0 -> too large, can be quickly recreated.. |
| `-SNPcaller` | string | `"MPI"` | stable | Consensus variant caller: MPI for mpileup or FB for FreeBayes. |
| `-SNPcores` | integer | `10` | stable | Maximum cores used by consensus SNP calling. |
| `-SNPconsMinDepth` | integer | `0` | stable | how many reads coverage to include position for consensus call? |
| `-SNPnormINDEL` | integer | `1` | stable | using bcftools norm to left-align indels |
| `-SVcaller` | integer | `0` | stable | calling structural variants: 1=delly, 2=gridss. Default (0). |

## Functional profiling (diamond)

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-profileFunct` | integer | `0` | stable | Run DIAMOND-based functional profiling. |
| `-reParseFunct` | integer | `0` | stable | Reparse existing DIAMOND hits. |
| `-reProfileFunct` | integer | `0` | stable | Remove and rebuild DIAMOND alignments and parsed profiles. |
| `-reProfileFuncTogether` | integer | `0` | stable | if any func database needs to be redone, than redo all indicated databases (useful if number of reads used changes..) |
| `-DiaCores` | integer | `12` | stable | CPU cores requested for DIAMOND jobs. |
| `-DiaMem` | integer | `7` | stable | memory in GB for diamond alignment jobs |
| `-DiaParseEvals` | string | `"1e-7"` | stable | evalues at which to accept hits to func database |
| `-DiaSensitiveMode` | integer | `0` | stable | Enable DIAMOND sensitive mode. |
| `-DiaFrameshift` | integer | `0` | stable | diamond -F frameshift penalty for long, error-prone reads; 0 disables frameshift-aware alignment |
| `-rmRawDiamondHits` | integer | `0` | stable | Delete raw DIAMOND hits after successful parsing. |
| `-DiaMinAlignLen` | integer | `20` | stable | Minimum accepted DIAMOND alignment length. |
| `-DiaMinFracQueryCov` | float | `0.1` | stable | Minimum accepted fraction of the query aligned. |
| `-DiaPercID` | integer | `40` | stable | Minimum accepted DIAMOND percent identity. |
| `-DiaDBs` | string | `""` | stable | Comma-separated functional databases: NOG,MOH,MOH2,ABR,ABRc,ACL,KGM,KGB,KGE,CZy,PTV,PAB,URE,URacc,AMI. See the [profiling tutorial](profiling_tutorial.md) for the config key each one needs. |

## Functional profiling (jaime tree)

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-orthoExtract` | integer | `0` | stable | Run translated-read orthologue placement. |

## Ribo profiling (mitag)

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-profileRibosome` | integer | `0` | stable | Run SSU/LSU read extraction and taxonomic assignment. |
| `-riobsomalAssembly` | integer | `0` | stable | Assemble extracted ribosomal reads. |
| `-reProfileRibosome` | integer | `0` | stable | delete RiboFind extraction, assignments and merged profiles, then rerun; implies -profileRibosome 1 |
| `-reRibosomeLCA` | integer | `0` | stable | delete RiboFind assignments and merged results, then rerun LCA; implies -profileRibosome 1 |
| `-riboMaxRds` | integer | `250000` | stable | Maximum extracted reads assigned per ribosomal marker. |
| `-saveRiboRds` | integer | `0` | stable | Retain reads used during ribosomal assignment. |
| `-thoroughCheckRiboFinish` | integer | `0` | stable | Require non-empty RiboFind assignment output before accepting completion. |

## Other tax profilers..

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-profileMetaphlan`, `-profileMetaphlan3` | integer | `0` | stable | Run MetaPhlAn taxonomic profiling. |
| `-profileProtal` | integer | `0` | stable | 1: profile each sample separately and merge; 2: one map, Protal across the complete cohort |
| `-ProtalCores` | integer | `4` | stable | CPU cores requested for a singular or combined Protal job. |
| `-ProtalMem` | integer | `100` | stable | Total memory in GB requested for a singular or combined Protal job. |
| `-protalIgnoreErrors` | integer | `1` | stable | 1: use the first compatible short-read pair and skip samples without one; 0: strict input validation |
| `-profileMOTU2` | integer | `0` | stable | Run mOTUs taxonomic profiling. |
| `-profileKraken` | integer | `0` | stable | Run Kraken taxonomic profiling. |
| `-profileTaxaTarget` | integer | `0` | stable | Run target-taxon profiling. |
| `-estGenoSize` | integer | `0` | stable | estimate average size of genomes in data |
| `-krakenDB` | string | `""` | stable | "virusDB";#= "minikraken_2015/"; |

## D2s distance

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-calcInterMGdistance` | integer | `0` | stable | Calculate the legacy inter-metagenome D2 distance matrix. |

## Io for specific uses

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-newFileStructure` | string | `""` | internal/advanced | just relink raw files for use in mocat |
| `-upload2EBI` | string | `""` | internal/advanced | copy human read removed raw files to this dir, named after sample |

## Institute specific: ei

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-wcKeyJobs` | string | `""` | internal/advanced | EI specific: work-category/accounting key added to submitted jobs |

## Debug

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-OKtoRWassGrps` | integer | `0` | internal/advanced | can delete assemblies, if suspects error in them |

## Legacy flag names

`MATAF4.pl` parses 171 options under 182 accepted names (aliases included); every one is listed above. Option-like strings found in [`manual_legacy.md`](manual_legacy.md) that `MATAF4.pl` does **not** accept belong to `geneCat.pl`, `MGS.pl`, `strain_within.pl` or `buildTree5.pl`, or have been removed. Look them up in the sections below before using them.

Removed with no replacement: `-useTrimomatic` (superseded by `-filterAdapters`), `-rmRawRds` (superseded by `-reduceScratchUse`), `-profileMetaphlan2` (use `-profileMetaphlan`), `-binSpeciesMG` as a separate binner switch (now an alias of `-Binner`).

## geneCat.pl

Gene-catalog construction and downstream gene-catalog annotation/MGS orchestration.

### Directories/files

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-o`, `-GCd` | string |  | stable | main save location for gene catalog and supporting files |
| `-tmp` | string | `-glbTmp` value | stable | working tmp dir; defaults to `-glbTmp` |
| `-glbTmp` | string | `<globalTmpDir>/GC/` | stable | global tmp dir, from config key `globalTmpDir` |
| `-map` | string | `?` | stable | mapping file(s); comma-separate several .map files to combine datasets (e.g. `-map file1.map,f2.map`) |

### Run modes

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-m`, `-mode` | string | `geneCat` | stable | Operation mode; see the list below. |

- `geneCat` — build or resume the complete gene-catalog workflow (default).
- `subprepSmpls` — internal per-batch sample collation, driven by `-SmplStart`/`-SmplStop`/`-SmplBatch`.
- `mergeCLs` — merge the clustering stages into the final catalog.
- `protExtract` — extract representative proteins; needs `-map` or a stored map.
- `FuncAssign` — functional annotation via `-functDB` and `-functAligner`.
- `FuncEMAP` — eggNOG-mapper annotation.
- `FOAM`, `ABR` — HMM-based FOAM and antibiotic-resistance annotation.
- `kraken`, `kaiju`, `specI` — taxonomic annotation of catalog genes.
- `FMG_extr` — extract the marker genes into a separate folder.
- `CANOPY` — canopy clustering of the gene abundance matrix.
- `ntMatchGC` — map catalog genes onto `-refDB` and write `-out`.

### Cluster options

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-clusterID` | integer | `95` | stable | percent identity at which the gene catalog is clustered |
| `-minGeneL` | integer | `100` | stable | minimal gene length for gene to be included in gene catalog, default: 100 |
| `-extraGenesNT` | string |  | stable | add genes (nt) from external sources, e.g. from complete genomes |
| `-extraGenesAA` | string |  | stable | add genes (AA) from external sources, e.g. from complete genomes |
| `-mmseqC` | integer | `1` | stable | 1: use mmseqs2 instead of CD-HIT for gene clustering |
| `-decluterMatrix` | integer | `0` | stable | declutering of gene matrix? by default deactivated, was more useful for canopy based MGS, can intro unwanted biases as long as gene/prots don't get removed |

### Flow control

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-1stepClust` | integer | `1` | stable | cluster incomplete genes separate? |
| `-submitLocal` | integer | `1` | stable | pretty important run mode switch, to submit jobs while geneCat is runnning single core |
| `-submSystem` | string | `""` | stable | qsub,SGE,bsub,LSF..; empty autodetects |
| `-continue`, `-justCDhit` | integer | `1` | stable | flow control, 1: continue with found files 0: delete existing (partial) gene cat, start again |
| `-c`, `-cores` | integer | `20` | stable | cores for the main gene-catalog jobs |
| `-c0`, `-cores0` | integer | `-1` | stable | specifcally cores only for the big main clustering job.. |
| `-c3`, `-cores3` | integer | `-1` | stable | for small jobs that really don't require that much power.. |
| `-mem` | integer | `200` | stable | max mem |
| `-stone` | string | `""` | advanced/internal | checkpoint file to touch when the selected `-mode` finishes |
| `-mem3` | integer | `-1` | stable | max mem for smaller jobs |

### Sample processing related

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-calcSupplCovSmpls` | integer | `1` | stable | if suppl reads were mapped, report gene abundances in these as separate samples (columns?) |
| `-oldStyleFolders` | integer | `-1` | deprecated/legacy | deprecated. only used for results calculated with an older MATAFILER version |
| `-requireAllAssemblies` | integer | `1` | advanced/internal | normally not exposed, continues even if some assemblies not present.. |
| `-sampleBatches` | integer | `-1` | stable | how many batches to use for initial accumulation of genes? (200-500 samples per batch recommended) |

### Binning/MGS related

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-binSpeciesMG` | integer | `2` | stable | use MAGs to create MGS? 1= metaBat2, 2=SemiBin, 3=metaDecoder, 4=GF, 5=SC |
| `-useCheckM2` | integer | `1` | stable | 1: use checkM2 completeness predictions, Default: 1 |
| `-useCheckM1` | integer | `0` | stable | 1: use checkM completeness predictions, Default: 0 |
| `-doStrains` | integer | `0` | advanced/internal | 1: calculate intraSpecific phylogenies on each MGS |
| `-SNPcaller` | string | `MPI` | stable | Consensus caller whose compressed MATAF4 outputs are used by the downstream strain workflow; accepted values are `MPI` and `FB`. |
| `-doMags` | integer | `1` | stable | 1: start canopy clustering, metabat2 & subsequent merging into MGS |
| `-canopyAutoCorr` | float | `0.15` | stable | canopy clustering parameter to filter autocorrelated genes prior to canopy clustering |

### Marker Genes/ taxonomy

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-MGset` | string | `GTDB` | stable | use either FMG or GTDB marker genes to compare and merge MAGs and calculate their abundance |

### Flags for specific modes

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-out` | string |  | stable | output dir, only used in modes protExtract ntMatchGC |
| `-functDB` | string | `KGM,TCDB,CZy,ABRc` | stable | for FuncAssign mode: functional DBs to annotate gene cat to |
| `-refDB` | string |  | stable | for ntMatchGC mode: reference fasta DB |
| `-fastaSplit` | string | `500M` | stable | For FuncAssign mode: split gene catalog into chunks to parallelise jobs. Default: 500M. |
| `-functAligner` | string | `diamond` | stable | either "diamond" or "foldseek" |
| `-SmplStart` | integer | `-1` | stable | for subprepSmpls |
| `-SmplStop` | integer | `-1` | stable | for subprepSmpls |
| `-SmplBatch` | integer | `-1` | stable | for subprepSmpls |

### Flags for functional assignment

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-FuncMinBitSc` | float | `45` | stable | minimum bit score for a functional assignment |
| `-FuncMinAlLeng` | integer | `30` | stable | minimum alignment length (AA) for a functional assignment |
| `-FuncMinPercSbjCov` | float | `0.5` | stable | Minimum fraction of subject coverage for functional assignment. |
| `-FuncMinPerID` | float | `25` | stable | minimum percent identity for a functional assignment |
| `-FuncMinEVal` | float | `1e-8` | stable | maximum e-value for a functional assignment |

## MGS.pl

MGS/MAG dereplication, abundance/taxonomy and optional strain workflow orchestration.

### General options

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-GCd` | string |  | stable | gene catalog dir |
| `-clusterID` | integer | `95` | stable | gene-catalog clustering identity percentage |
| `-outD` | string | `<GCd>/Bin_SB/` | stable | output dir |
| `-tmp` | string |  | stable | temp dir |
| `-submit` | integer | `1` | stable | 1:submit jobs, 0: dry run. Default: 1 |
| `-canopies` | string |  | stable | location of canopy clustering output file (clusters.txt) |
| `-smallCores` | integer | `4` | stable | cores used for normal jobs (not intensive) |
| `-bottleneckCores` | integer | `12` | stable | cores for compute intensive jobs |
| `-redoCluster` | integer | `0` | stable | delete and redo the clusterMAGs dereplication |
| `-redoTax` | integer | `0` | stable | rewrite tax annotations |
| `-MGset` | string | `GTDB` | stable | marker genes used for MAG merging/abundance: `GTDB` or `FMG` |
| `-wait4stone` | string |  | stable | wait for these files to be created, refers currently exclusively to eggNOG annotations that are needed later |
| `-wait4stoneTimeout` | integer | `86400` | stable | maximum wait in seconds; `0` waits indefinitely |
| `-mem` | integer | `150` | stable | memory used for intensive jobs |
| `-strains` | integer | `0` | stable | 1: calc instra species strain phylogenies. Default: 0 |
| `-redo` | string | `none` | stable | strain-workflow redo mode forwarded downstream: `none`, `tree`, `input` or `all` |
| `-prepareMosaicLoci` | integer | `1` | stable | 1: confirm Mosaic loci/outgroups before strain analysis; 0: keep same-NOG seed clusters separate |
| `-SNPcaller` | string | `MPI` | stable | consensus caller whose inputs are forwarded to `strain_within.pl`: `MPI` or `FB` |
| `-useCheckM2` | integer | `0` | stable | CheckM2 default qual checking of MAGs/MGS |
| `-useCheckM1` | integer | `1` | stable | CheckM default qual checking of MAGs/MGS |
| `-binSpeciesMG` | integer | `2` | stable | 0=no, 1=metaBat2, 2=SemiBin, 3: MetaDecoder, 4 ,5 |
| `-ignoreIncompleteMAGs` | integer | `1` | stable | 1: assemblies without MAG calculations are ignored. Default: 1 |
| `-legacy` | integer | `0` | deprecated/legacy | pre-Dec-2022 clustering; needs `-MGset FMG`. No longer supported |
| `-perlClusterMAGs` | flag | off | advanced/internal | use the Perl clusterMAGs path instead of the clusterMAGs binary (compatibility/debug) |
| `-genomesPerFamily` | integer | `0` | advanced | extract bins per family (or per assembly group/sample when the family is missing) |

## strain_within.pl

Within-MGS locus extraction, quality control, tree orchestration and downstream hand-off. Normally `MGS.pl` supplies the catalogue paths; see the [strain-within workflow guide](strainwithin.md) before invoking this script directly.

For split Phase I, the parent scans the catalogue cluster index once and atomically publishes a provenance-bound binary membership shard for each worker under shared strain scratch. The manifest is published only after every shard is complete. Workers validate the generation, size, header, record count and payload digest before use; an absent, stale or corrupt cache falls back to the original full-index parser.

The parent also publishes one common binary subset of the selected catalogue proteins. All Phase-I workers use that subset for sequence-aware locus grouping and ambiguous-copy selection instead of independently scanning the full catalogue FAA. Catalogue nucleotide sequences are deliberately not copied at this point: Phase-I worker nucleotide records come from sample consensus FASTAs, while the parent determines and streams the exact catalogue FNA outgroup set only after Phase I. Creating an early catalogue FNA subset would add a large read and would still miss candidate genes when `-outgroupReferenceGeneCap` exceeds `-presortGenes`.

Phase-I input contracts use cross-node-stable file identity: canonical path, inode, size and modification time, plus the persistent catalogue identity. Filesystem device numbers are excluded because shared files may have different device numbers in separate compute-node mount namespaces. Version-2 contracts remain readable and a resumed parent upgrades a compatible legacy contract before dispatching workers.

### Inputs, execution and workflow control

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-GCd` | string | required | stable | Gene-catalogue directory. |
| `-MGS` | string | required unless `-outD` is given | stable | Nonempty MGS guide. |
| `-outD` | string | derived | stable | Output directory; also permits a targeted resume without `-MGS`. |
| `-map2` | string | `""` | stable | Sample metadata map forwarded to `strain_within_2.2.pl`. |
| `-MGSabundance` | string | auto | stable | Explicit MGS abundance matrix, recommended for a nonstandard guide location. |
| `-clusterID` | integer | `95` | stable | Gene-catalogue clustering identity; must be 1–100. |
| `-MGset` | string | `GTDB` | stable | Marker-gene set; accepted values are `GTDB` and `FMG`. |
| `-SNPcaller` | string | `MPI` | stable | Select caller-specific compressed consensus FASTA/VCF inputs; accepted values are `MPI` and `FB`. |
| `-MGSphylo` | string | `""` | stable | Source MGS tree forwarded for fallback outgroup selection. |
| `-tmpD` | string | `""` | stable | Node-local temporary directory. |
| `-submit` | integer | `0` | stable | Submit generated work (`1`) or perform a dry run (`0`). |
| `-submissionMode` | string | auto | stable | Scheduler/backend override; dry runs default to `bash`. |
| `-maxSubJob` | integer | `-1` | advanced | Phase-I worker count: `-1` auto-selects (50-150 assembly groups per worker), `0` disables splitting, positive is explicit. |
| `-MGSsubset` | string | `""` | advanced | Comma-separated MGS subset. |
| `-help`, `-h` | flag | `0` | stable | Show the built-in usage text. |

### Resume, redo and resources

On a tree-only resume (`-onlySubmit 1`), an MGS whose placement never finished is retried automatically in EPA-only mode:

- A retained `placementPending.sto` plus a validated backbone, MSA, query alignment and sample classification is resubmitted with `-epaOnly 1`, reusing the existing tree instead of rerunning alignment or inference.
- The retry runs on one core with `-epaThreads 1`, and doubles the ordinary memory estimate, then clamps it to a fixed 20 GB floor and a `-phyloMemMulti`-scaled 220 GB ceiling. The doubled figure is a scheduler/cgroup allowance, not an EPA-ng `--maxmem` argument.
- Legacy runs are recovered too: when those placement inputs exist but `phylo/IQtree_allsites.treefile` does not, the pending marker is reconstructed and BuildTree restarts in the same isolated EPA-only mode.

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-onlySubmit` | integer | `0` | stable | Reuse completed Phase-I preparation and submit only missing tree work. |
| `-onlyMSA` | integer | `0` | stable | Stop after the per-locus nucleotide alignments, before concatenation, phylogeny, placement and `strain_within_2.2.pl`. Records `msaOnly.complete.tsv`; resume with `-onlySubmit 1`. Incompatible with `-placeOnBackbone 1`, `-redo tree` and `-redoEPAfilter`. |
| `-redo` | string | `none` | stable | Destructive recovery: `none` resumes; `tree` rebuilds trees only; `input` rebuilds incomplete inputs and dependent trees; `all` rebuilds everything. Combine with `-MGSsubset` (not with `all`). |
| `-redoEPAfilter` | optional integer | `0` | advanced | Rebuild each final EPA-placed tree from its retained jplace and backbone, then continue through normal controller validation and downstream strain analysis. Bare flag implies `1`. |
| `-maxCores` | integer | `-1` | stable | Cap BuildTree jobs at this many cores; the request is `ceil(sqrt(samples))` with a four-core floor. `-1` keeps the fixed request. |
| `-selfMemGb` | string | `10` | stable | Memory in GiB for the controller and Phase-I workers. `auto` (or `-1`) sizes each worker from the cluster-index shard it loads, the resolved locus model and its assembly-group count, and prints the chosen value and its terms. Confirmed worker OOM retries double the request. |
| `-mosaicMemGb` | integer | `150` | stable | Total memory in GiB for the Mosaic prerequisite. |
| `-phase1WorkerRetries` | integer | `2` | advanced | Retry count for invalid Phase-I workers, 0-10. Confirmed OOM retries request more memory and use the larger `-oomMinRetries` budget. |
| `-treeOOMMaxMemGB` | float | config `maxMF4mem`, else `512` | advanced | Memory ceiling for automatic Phase-I worker and tree OOM retries. |
| `-treeOOMRetryRounds` | integer | `8` | advanced | Maximum OOM retry rounds per MGS, 0-12. High enough that `-treeOOMMaxMemGB` decides when to stop. |
| `-oomScanMinutes` | float | `60` | advanced | Rescan Slurm accounting this often while Phase-I workers or tree jobs still run, resubmitting whatever was killed for memory. Both phases submit their largest jobs first, so escalation no longer waits for the tail of short jobs. |
| `-oomMinRetries` | integer | `3` | advanced | Minimum OOM escalations guaranteed per failed job, 0-12; `0` removes the floor. Budgets are counted per job, so an early rescan cannot spend the retries a later failure needs. |
| `-jobNice` | integer | `5000` | advanced | Slurm `--nice` handicap on ordinary tree and Phase-I jobs; OOM retries always submit at `--nice=0` so they outrank a backlog this run already queued. Rescanning alone cannot reorder jobs Slurm has already accepted. `0` disables it. |
| `-maxQueuedJobs` | integer | `0` | advanced | Ceiling on this user's live (running + pending) scheduler jobs, so the tree wave is submitted in batches and an OOM retry only has to overtake what is already queued. `0` submits the whole wave at once. |
| `-treeMemThreadDivisor` | float | `4` | advanced | Scale the initial tree memory request by `cores/DIVISOR` (never below 1); raise it to request less on wide jobs. |
| `-phyloMemMulti` | float | `1` | advanced | Multiplier applied to BuildTree memory planning. |
| `-flushMemMB` | integer | `2048` | advanced | Also flush buffered Phase-I output once this much is held; lower it if split workers run out of memory. |

### Extraction and biological QC

Questionable sample-locus observations are masked first. A sample still carrying excessive ambiguity among the remaining observations is kept for post-tree placement rather than used to infer the backbone.

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-presortGenes` | integer | `1200` | stable | Potential loci considered before final tree selection. |
| `-maxGenes` | integer | `600` | stable | Maximum validated loci retained per MGS/sample; values `<=0` remove this cap without disabling QC. |
| `-treeLocusBudget` | integer | `400` | stable | Maximum final loci selected for each tree; the candidate alignment stays bounded to this plus the QC-backfill allowance. |
| `-MGSminGenesPSmpl` | integer | `8` | stable | Extraction prefilter: validated loci a sample needs for one MGS before it is written. Tree inclusion is decided later by `-GenesPerSpecies`/`-relativeNTFraction`/`-NTfiltCount`. |
| `-multiGeneSmplMax` | float | `0.25` | stable | Maximum ambiguous/multigene-locus fraction per sample. |
| `-conspGeneSmplMax` | float | `0.05` | stable | Maximum conspecific-signal locus fraction per sample. |
| `-minBadLociPSmpl` | integer | `3` | stable | Minimum bad loci before a sample is flagged as mixed strain. |
| `-excludeMixedStrainSamples` | integer | `1` | stable | Drop samples breaching `-multiGeneSmplMax` or `-conspGeneSmplMax` from their MGS tree (passes `-excludeFlaggedSamples 1` to BuildTree). Never considers locus counts. |
| `-enforceSampleCoverage` | integer | `1` | stable | Make `-GenesPerSpecies`/`-relativeNTFraction`/`-NTfiltCount` the sample-inclusion filter. `0` retains every aligned sample. |
| `-minLociPerMGS` | integer | `8` | stable | Distinct loci an MGS needs before a tree is attempted; a property of the MGS, not of a sample. |
| `-breakpointGeneFlank` | integer | `50` | stable | Bases around mapping breakpoints in which genes are masked. |
| `-abundanceMinLoci` | integer | `8` | stable | Minimum loci for robust abundance-pattern filtering. |
| `-abundanceMinFold` | float | `0.333333…` | stable | Lowest accepted locus/median depth ratio. |
| `-abundanceMaxFold` | float | `3` | stable | Highest accepted locus/median depth ratio. |
| `-abundanceMaxModifiedZ` | float | `3.5` | stable | Modified-Z threshold for depth outliers. |
| `-forceSNPcalls` | integer | `0` | advanced | Force regeneration of consensus FASTA from VCF. |
| `-preCompConsSNP` | integer | `0` | advanced | Precompute consensus SNPs in this many blocks; `0` disables it. |
| `-minSNPDepth` | integer | `2` | stable | Minimum SNP depth. |
| `-minSNPCallQual` | integer | `20` | stable | Minimum SNP call quality. |
| `-skipIndels` | integer | `1` | stable | Exclude indels from consensus processing. |
| `-SNPadaptiveQual` | float | `0` | advanced | Adaptive depth-based quality filtering; `0` disables it. |
| `-SNPdepthFilterScale` | float | `0.15` | stable | Filter when depth is below mean contig depth multiplied by this value. |
| `-SNPindelRangeFilt` | integer | `5` | stable | Exclude SNPs within this many bases of indels. |

### Internal and deprecated compatibility options

These options remain parsed so existing generated commands and worker scripts do not break, but they are not part of the normal user-facing interface.

| Option | Status | Replacement or purpose |
|---|---|---|
| `-nodeTmp` | deprecated | Use `-tmpD`. |
| `-reSubmit`, `-recalcTrees` | deprecated | Use `-redo tree`. |
| `-repairCAT`, `-deepRepair` | deprecated | Use `-redo input`. |
| `-redoSubmissionData` | deprecated | Use `-redo all`. |
| `-strictBackbone` | deprecated | Use `-placeOnBackbone`. |
| `-treeSubFromMGS` | deprecated | Use `-MGSsubset`. |
| `-noGeneLimit` | deprecated | Use `-maxGenes 0`. |
| `-subjob` | internal | Split-worker index supplied by the parent process. |
| `-flushEvery` | internal | Publish buffered Stage-I records after this many sample rows; default `50`. Lower values cut peak memory but reopen shard files more often. |
| `-disableQC` | internal | Developer override that disables biological QC. |

The parent canonicalizes all unlimited extraction spellings to `-maxGenes 0` when it creates split-worker commands. It does not forward the deprecated `-noGeneLimit` alias. Therefore the presence of `-maxGenes 0` in a worker script explicitly means that no per-sample locus-count cap is applied; biological QC and the separate `-presortGenes` and `-treeLocusBudget` bounds remain active.

### Mosaic and outgroup preparation

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-mosaicLoci` | string | auto | stable | Confirmed catalogue-wide mosaic/outgroup table. |
| `-mosaicMGS` | string | derived from `-MGS` | stable | Raw `SB.clusters` assignment used for comprehensive Mosaic discovery; inferred by removing `.core` from `-MGS`. |
| `-prepareMosaicLoci` | integer | `1` | stable | Create, wait for and validate a missing Mosaic catalogue. |
| `-outgroupCoreMinLoci` | integer | auto | advanced | Candidate outgroup overlap floor; `0` derives 20% of `-treeLocusBudget`. |
| `-outgroupReferenceGeneCap` | integer | `2500` | advanced | Candidate reference genes retained per outgroup MGS. |
| `-preferredCoreGenes` | string | auto | advanced | Universal-core guide used to prioritise eligible loci. |

### Tree selection and placement

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-GeneLengthMin` | float | `0.3` | stable | Minimum fraction of the locus length-Q90 used for sample/locus QC. |
| `-GeneLengthIncludeMin` | float | `0.03` | stable | Lower post-QC fraction allowed into backbone and placement MSA input; it cannot exceed `-GeneLengthMin`. |
| `-GenesPerSpecies` | float | `0.2` | stable | Backbone minimum relative locus coverage per sample. |
| `-relativeNTFraction` | float | `0.1` | stable | Backbone minimum relative informative-NT coverage. |
| `-NTfiltCount` | integer | `5000` | stable | Absolute backbone floor on informative NT after final MSA; roughly five well-covered loci. `0` leaves only the relative gates. |
| `-placementGenesPerSpecies` | float | `0.04` | advanced | Placement minimum relative locus coverage. |
| `-placementRelativeNTFraction` | float | `0.03` | advanced | Placement minimum relative informative-NT coverage. |
| `-placementNTfiltCount` | integer | mirrors `-NTfiltCount` | advanced | Placement absolute informative-NT floor. |
| `-taxonAwareLocusSelection` | integer | `1` | stable | Enable robust-core plus taxon-rescue locus selection after MSA QC. |
| `-taxonAwareRescueMinPrevalence` | float | `0.8` | advanced | Minimum usable-taxon prevalence for rescue/backfill loci. |
| `-compactTaxonAwareDiagnostics` | integer | `1` | stable | Merge final taxon-aware/rate audit tables into one diagnostic file. |
| `-rateMergePartitions` | integer | `1` | stable | Merge loci into deterministic divergence/GC partition bins. |
| `-rateMergeMaxBins` | integer | `8` | advanced | Maximum deterministic partition bins. |
| `-rateMergeTargetSites` | integer | `30000` | advanced | Target effective called sites per initial bin. |
| `-rateMergeMinLoci` | integer | `20` | advanced | Minimum loci per bin before nearest-bin merging. |
| `-rateMergeMinSites` | integer | `20000` | advanced | Minimum effective sites per bin before merging. |
| `-placeOnBackbone` | integer | `0` | advanced | Infer a well-covered backbone and place eligible sparse samples with EPA-ng. When `0`, every backbone- and placement-only filter is inactive. |
| `-strictBackboneFraction` | float | `0.35` | advanced | Severe aligned-coverage deferral threshold relative to sample Q90. |
| `-strictBackboneMinSamples` | integer | `3` | advanced | Minimum backbone samples before falling back to the complete alignment. |
| `-placementMinOverlap` | integer | `10000` | advanced | Minimum informative positions shared with the inferred backbone. |
| `-epaThreads` | integer | `2` | advanced | Requested EPA-ng threads; BuildTree caps them by available cores and by one thread per GB of planning memory. |
| `-epaMaxMemMB` | integer | `-1` | advanced | EPA-ng thread-planning budget; `-1` derives 60% of each IQ-TREE allowance, `0` disables memory-based scaling. |
| `-epaPendantOutlierFactor` | float | `5` | advanced | Pendant-branch cutoff multiplier; `0` disables this filter. |
| `-epaPendantMinThreshold` | float | `0.02` | advanced | Minimum pendant-branch cutoff in substitutions/site. |
| `-MSAprog` | integer | `2` | stable | Alignment program: 0 MSAProbs, 1 Clustal Omega, 2 MAFFT, or 4 MUSCLE5. |
| `-phyloProg` | integer | `1` | stable | Tree program: 1 IQ-TREE, 2 VeryFastTree, or 3 FastTree. |
| `-iqPathogen` | integer | `0` | advanced | Opt in to the IQ-TREE 3 pathogen/CMAPLE path. |
| `-rmMSA` | integer | `1` | stable | Remove per-locus MSAs unless downstream analysis requires them. |

### Downstream analysis and metadata

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-popGenStats` | integer | `0` | stable | Enable population-genetic analysis and retain required nucleotide MSAs. Off by default because it is much slower than the strain statistics; setting it to `1` also forces `-rmMSA 0`. |
| `-popGenCategory` | string | `""` | stable | Comma-separated metadata columns whose levels each receive their own population-genetic analysis; forwarded to `popGenStats.R --category`. Requires `-popGenStats 1`. |
| `-popGenStrictOutgroup` | integer | `0` | stable | Require the requested outgroup for population-genetic analysis. |
| `-popGenGeneticCode` | integer | `1` | stable | Genetic code forwarded to population-genetic analysis. |
| `-popGenCodonStart` | integer | `1` | stable | Codon frame start (1, 2 or 3). |
| `-popGenSeed` | integer | `1` | stable | Non-negative reproducibility seed. |
| `-popGenLegacyTextOutput` | integer | `0` | deprecated/legacy | Also emit legacy population-genetics text output. |
| `-ContTests` | string | `""` | stable | Comma-separated continuous metadata tests forwarded downstream. |
| `-DiscTests` | string | `""` | stable | Comma-separated discrete metadata tests forwarded downstream. |
| `-familyVar` | string | `""` | stable | Metadata column containing family identifiers. |
| `-groupStabilityVars` | string | `""` | stable | Metadata columns used for resilience/persistence calculations. |
| `-individualVar` | string | `AssmblGrps` | stable | Metadata column containing individual identifiers. |

## strain_within_2.2.pl

Within-MGS postprocessing, strain statistics and population-genetic analysis. It is normally launched by `strain_within.pl`; direct use is primarily for targeted postprocessing restarts.

### General options

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-GCd` | string | required | stable | Gene-catalogue directory. |
| `-FMGdir` | string | required | stable | Within-phylogeny directory containing the per-MGS results. |
| `-map` | string | required | stable | Nonempty sample metadata map. |
| `-MGSmatrix` | string | required | stable | Nonempty MGS abundance matrix. |
| `-MGSphylo` | string | `""` | stable | Source MGS tree used to recover fallback outgroups. |
| `-submit` | integer | `1` | stable | Submit analyses (`1`) or generate a dry-run plan (`0`). |
| `-qsubSystem` | string | auto | advanced | Queue backend override; dry runs default to `bash`. |
| `-cores` | integer | `4` | stable | Cores for standard R analyses. |
| `-Hcores` | integer | `12` | stable | Cores for heavier downstream analyses. |
| `-reSubmit` | integer | `0` | advanced | Clear and redo all within-MGS postprocessing and both statistic types. |
| `-redoStrainStats` | integer | `0` | stable | Redo only strain statistics and dependent postprocessing checkpoints. |
| `-redoPopGenStats` | integer | `0` | stable | Redo only population-genetic stores and summaries; requires `-popGenStats 1`. |

### Population genetics and metadata

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-popGenStats` | integer | `0` | stable | Enable per-MGS population-genetic analysis. Off by default because it is much slower than the strain statistics. |
| `-popGenSubsample` | string | `10,20,30,100,200,500` | stable | Comma-separated population-genetic subsample sizes. |
| `-popGenCategory` | string | `""` | stable | Comma-separated metadata columns passed to `popGenStats.R --category`. Each level gets its own full analysis (dN/dS, Tajima's D, nucleotide diversity, …) over the samples in that level, alongside the ungrouped result. Join columns with `BLOCK` to analyse each combination of their values. Requires `-popGenStats 1`. Independent of `-DiscTests`, where `BLOCK` instead means a permutation-test blocking factor. |
| `-popGenStrictOutgroup` | integer | `0` | stable | Require the requested outgroup. |
| `-popGenGeneticCode` | integer | `1` | stable | Positive genetic-code identifier. |
| `-popGenCodonStart` | integer | `1` | stable | Codon frame start (1, 2 or 3). |
| `-popGenSeed` | integer | `1` | stable | Non-negative reproducibility seed. |
| `-popGenLegacyTextOutput` | integer | `0` | deprecated/legacy | Also emit legacy population-genetics text output. |
| `-DiscTests` | string | `""` | stable | Discrete metadata variables tested for phylogenetic signal. |
| `-ContTests` | string | `""` | stable | Continuous metadata variables tested for phylogenetic signal. |
| `-familyVar` | string | `""` | stable | Metadata column containing family identifiers. |
| `-groupStabilityVars` | string | `""` | stable | Metadata columns used for resilience/persistence calculations. |
| `-individualVar` | string | `AssmblGrps` | stable | Metadata column containing individual identifiers. |

## buildTree5.pl

Phylogenetic tree construction and related MSA/population-genetic analyses.

### General options

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-genoInD` | string |  | stable | provide a dir with complete genomes, will extract FGMs and build tree between genomes (NT/AA flag) |
| `-wildcardflag` | string |  | stable | Glob pattern selecting the per-locus input FASTAs when `-fna`/`-aa` are not given. |
| `-fna` | string |  | stable | Nucleotide FASTA input. |
| `-aa` | string |  | stable | Amino-acid FASTA input. Note: source flag is -aa, not -faa. |
| `-cats` | string |  | stable | Gene-category (locus) assignment file linking sequences to loci. |
| `-outD` | string |  | stable | Output directory. |
| `-tmpD` | string |  | stable | Scratch directory holding all plain MSA, concatenated and tree files. |
| `-stagedInputDir` | string | `""` | advanced/internal | Directory holding input already staged by the caller; set by `strain_within.pl`. |
| `-tmpSubdir` | string | `""` | advanced/internal | Job-specific subdirectory created under `-tmpD`. |
| `-completionMarker` | string | `""` | advanced/internal | File written once the run finishes successfully. |
| `-terminalMarker` | string | `""` | advanced/internal | File written when the run ends in a valid terminal state that must not be retried. |
| `-placementPendingMarker` | string | `""` | advanced/internal | File marking a backbone that still needs EPA-ng placement, so `-epaOnly` can resume it. |
| `-withinSpecies` | integer | `0` | stable | Enable within-species filtering: require two called sequences per retained column by default and enable divergence-based locus QC. |
| `-strainWithinPreset` | integer | `0` | advanced/internal | Apply the fixed MATAFILER strain-tree preset; implies `-withinSpecies 1`. |
| `-cores` | integer | `1` | stable | CPU cores used by MSA and tree programs. |
| `-superTree` | integer | `0` | stable | Build a supertree from the per-locus trees instead of one concatenated tree. |
| `-superCheck` | integer | `0` | stable | Only check whether the supertree inputs are complete, then exit. |
| `-fixHeaders` | integer | `0` | stable | fix the fasta headers, if too long or containing not allowed symbols (nwk reserved) |
| `-useEte` | integer | `0` | stable | Use the ete3 toolkit workflow instead of the internal MSA/tree steps. |
| `-relativeNTFraction` | float | `0.2` direct; `0.1` strain workflow | stable | Backbone minimum informative nucleotide content as a fraction of the final selected-sample Q90. This explicit name replaces the retired ambiguous `-NTfilt` switch. |
| `-NTfiltPerGene` | float | `0.1` | stable | Per-locus minimum informative-NT fraction; sequences below it are dropped from that locus. |
| `-GeneLengthIncludeMin` | float | `0.03` | stable | Lower per-sample locus-length-Q90 fraction admitted to MSA only after QC at `-NTfiltPerGene`; recovered-only observations do not count toward backbone/placement eligibility. |
| `-GenesPerSpecies` | float | `0.1` | stable | Backbone minimum retained-locus count as a fraction of the final selected-sample Q90. |
| `-fracMaxGenes90pct` | float | `0.25` | stable | Minimum locus size as a fraction of the category-size Q90. Set to `0` to retain every category containing at least one length-filtered sequence. |
| `-NTfiltCount` | integer | `0` | stable | Backbone absolute minimum informative nucleotide count after final locus selection. |
| `-placementRelativeNTFraction` | float | mirrors direct; `0.03` strain workflow | advanced | Placement equivalent of the relative NT coverage filter. The strain default retains samples with at least 3% of the selected-sample NT Q90, subject to the absolute shared-backbone overlap floor. |
| `-placementGenesPerSpecies` | float | mirrors direct; `0.04` strain workflow | advanced | Placement equivalent of the relative locus-count filter. The strain default requires 4% of the selected-sample locus Q90, and a minimum of two loci is always enforced. |
| `-placementNTfiltCount` | integer | mirrors `-NTfiltCount` | advanced | Placement equivalent of the absolute informative-NT filter, before `-placementMinOverlap` is applied. |
| `-smplDef` | integer | `1` | stable | is the genome somehow quantified with a delimiter (_) ? |
| `-smplSep` | string | `_` | stable | set the delimiter |
| `-outgroup` | string |  | stable | Retained as the anchor for filtering, downstream backbone rooting, and placement; not passed to IQ-TREE, where `-o` is cosmetic under reversible models and can fail on an internal sparse tree. |
| `-AAtree` | integer | `0` | stable | Build the tree from amino acids instead of nucleotides. |
| `-MSAprogram` | integer | `2` | stable | (0) MSAprobs, (1) clustalO, (2) mafft, (4) MUSCLE5, (5) FAMSA2 (only AA) |
| `-MSAfixRecoverTechnicalOffsets` | integer | `1` | advanced | Repair coding-NT gap offsets after back-translation (MSAfix single-alignment mode only). |
| `-MSAfixCodingFrame` | integer | `1` | advanced | Reading frame assumed when repairing coding nucleotide alignments. |
| `-MSAfixGeneticCode` | integer | `11` | advanced | Genetic code used by MSAfix; `11` is bacterial/archaeal/plastid. |
| `-MSAfixRecoveryBand` | integer | `3` | advanced | Codon window searched on each side when recovering a technical offset. |
| `-minOverlapMSA` | float | `0` between species; `0.35` within species | stable | Minimum called-sequence fraction required to retain an MSA column. This is passed to MSAfix v2.14+; BuildTree converts it with `ceil(sequences × fraction)` for its final concatenation, giving the same retained columns. An explicit value overrides the `-withinSpecies` default. |
| `-maxGapPerCol` | float | `1` | stable | same as minOverlapMSA, but for MSAfix and %of gaps allowed in a column |
| `-calcDistMat` | integer | `0` | stable | Write a distance matrix for the alignment type used by the tree. |
| `-calcDistMatExt` | integer | `0` | stable | Also write a distance matrix for the other alignment type (runs a second MSA). |
| `-calcDiffDNA` | integer | `0` | stable | Report pairwise nucleotide differences between samples. |
| `-minPcId` | float | `0` | stable | sequence is filtered from data, unless the average minPcId is >= $minPcId |
| `-SynTree` | integer | `0` | stable | Additionally build a tree from synonymous sites only. |
| `-NonSynTree` | integer | `0` | stable | Additionally build a tree from non-synonymous sites only. |
| `-continue` | integer | `0` | stable | Reuse existing intermediate files instead of overwriting them. |
| `-onlyMSA` | integer | `0` | stable | Finish the existing per-locus alignment pipeline, including filtering, NT backtranslation, localized MSAfix, and `.gz` checkpoint publication, then exit before combined-MSA postprocessing and `mergeMSAs`. Writes `msaOnly.complete.tsv`, retains non-merged per-locus `MSA/*.fna.gz` files, and skips partitions, phylogeny, and EPA-ng. Incompatible with `-placeOnBackbone 1`. |
| `-epaOnly` | integer | `0` | advanced | Resume a placement-pending run: reuse the retained backbone and rerun only EPA-ng placement and publication. |
| `-redoEPAfilter` | optional integer | `0` | advanced | With `-continue 1`, rerun only placement filtering, its audit report and final-tree publication from retained backbone/jplace artifacts. Bare flag implies `1`. |
| `-bootstrap` | integer | `0` | stable | Bootstrap replicates; `0` disables bootstrapping. |
| `-subsetSmpls` | integer | `-1` | stable | Randomly subsample to this many samples; `-1` keeps all. |
| `-postFilter` | string |  | stable | "," sep list of zorro,guidance2,macse |
| `-rmMSA` | integer | `1` | stable | Remove per-locus checkpoints at normal terminal completion. With `0`, retain per-locus nucleotide alignments as `.fna.gz`; protein checkpoints are still removed. `-onlyMSA 1` exits before finalization and therefore retains its per-locus checkpoints regardless. Active plain MSAs are scratch files in either mode. |
| `-gzInput` | integer | `0` | stable | Compress buildTree-owned input files after a successful run; retained `MSAli*.fna` output is always stored as `.gz`. |
| `-isAligned` | integer | `0` | stable | Treat the input as already aligned and skip the MSA step. |
| `-runRAxML` | integer | `0` | stable | Build the tree with RAxML. |
| `-runRaxMLng` | integer | `0` | stable | Build the tree with RAxML-NG. |
| `-runFastTree` | integer | `0` | stable | Build the tree with FastTree. |
| `-runVeryFastTree` | integer | `0` | stable | Build the tree with VeryFastTree. |
| `-treeShrink` | integer | `0` | stable | Run TreeShrink to remove long-branch outlier tips. |
| `-placeOnBackbone` | integer | `0` | advanced | Infer an ML backbone from well-covered samples and place severe coverage outliers with EPA-ng. When `0`, all other backbone- and placement-only options are ignored. |
| `-strictBackbone` | integer | `0` | deprecated | Compatibility alias for `-placeOnBackbone`. |
| `-strictBackboneFraction` | float | `0.35` | advanced | Defer a sample only when its called alignment sites are below this fraction of the sample Q90. Locus-QC status alone does not defer a sample after questionable loci have been masked. |
| `-strictBackboneMinSamples` | integer | `3` | advanced | Fall back to the complete alignment if coverage filtering would leave fewer backbone samples. |
| `-placementMinOverlap` | integer | `10000` | advanced | Absolute minimum called alignment positions required both before placement and in actual overlap with the inferred backbone. This complements the relative locus/NT gates and is the primary sparse-placement reliability control. |
| `-epaPendantOutlierFactor` | float | `5` | advanced | Exclude EPA-ng placements whose pendant branch exceeds this multiple of the non-outgroup backbone terminal-branch Q95. Set to `0` to disable. |
| `-epaPendantMinThreshold` | float | `0.02` | advanced | Minimum adaptive pendant-branch cutoff in substitutions/site, preventing over-filtering of very compact backbones. |
| `-epaThreads` | integer | `4` | advanced | Requested EPA-ng threads, capped by available cores and by one thread per GB of planning memory. |
| `-epaMaxMemMB` | integer | `-1` | advanced | EPA-ng thread-planning budget; `-1` derives 60% of each IQ-TREE allowance, `0` disables memory-based scaling. |
| `-sampleQC` | string |  | advanced | Optional per-sample QC table. Rows marked `placement` name samples the caller judged unfit for the tree; `-excludeFlaggedSamples` decides whether they are dropped. Written by `strain_within.pl`; other callers can omit it. |
| `-excludeFlaggedSamples` | integer | `0` | stable | Remove samples that the `-sampleQC` table marks unfit, before any length or prevalence statistic is taken, so they cannot shift the per-locus Q90 length reference or locus occupancy. A generic mechanism: what counts as unfit is decided by the caller that writes the table. No-op without `-sampleQC`. Never considers locus counts or coverage. |
| `-enforceSampleCoverage` | integer | `0` | stable | Apply `-GenesPerSpecies`, `-relativeNTFraction` and `-NTfiltCount` as a sample removal in the taxon-aware branch, which otherwise retains every sample holding one informative site. Matches what the `-taxonAwareLocusSelection 0` prefilter has always done. Inactive under `-placeOnBackbone 1`, where those samples are routed to placement instead. Off by default so callers that set no inclusion policy keep their previous behaviour. |
| `-postAlignmentLocusQC` | integer | `0` between species; `1` within species | stable | Run native MSAfix locus-comparability QC before multi-locus concatenation. Broad trees retain every prepared locus unless QC is explicitly enabled. |
| `-postAlignmentMinSequences` | integer | `3` | stable | Minimum comparable sequences required for an aligned locus. |
| `-postAlignmentMinOccupancy` | float | `0.35` | stable | Minimum fraction of unambiguous alignment cells; permissive for metagenomic loci. |
| `-postAlignmentDivergenceQC` | integer | `0` between species; `1` within species | stable | Reject absolute and cross-locus divergence outliers. When locus QC is enabled, its structural checks remain active even if divergence QC is disabled. |
| `-postAlignmentRelativeZ` | float | `5.0` | stable | Modified-Z threshold for cross-locus consensus-divergence outliers when divergence QC is enabled. The stricter within-species default rejects anomalously fast loci, not samples. |
| `-postAlignmentMinLociRelative` | integer | `8` | stable | Minimum locus count before applying cross-locus robust outlier QC. |
| `-rateMergePartitions` | integer | `0` direct; `1` strain preset | stable | Replace one partition per retained locus with deterministic divergence/GC bins for the primary nucleotide tree. `strain_within.pl` enables this by default. Its audit is a section of `phylo/taxon_aware_diagnostics.tsv` unless `-compactTaxonAwareDiagnostics 0` is set. |
| `-rateMergeMaxBins` | integer | `8` | advanced | Upper bound on deterministic bins. |
| `-rateMergeTargetSites` | integer | `30000` | advanced | Target effective called sites per initial partition. The initial target is `ceil(total effective called sites / target)`, then capped by `-rateMergeMaxBins`. |
| `-rateMergeMinLoci` | integer | `20` | advanced | Merge a bin with its nearest normalized divergence/GC neighbour while it contains fewer loci than this threshold. |
| `-rateMergeMinSites` | integer | `20000` | advanced | Merge a bin with its nearest neighbour while its effective called sites are below this threshold. A locus contributes its mean number of called bases across retained taxa, so missing data reduce support. |
| `-taxonAwareLocusSelection` | integer | `1` | stable | Enable two-stage locus selection. A permissive robust/core-plus-rescue candidate set is aligned first; the final set is chosen after MSA QC using occupancy and parsimony-informative sites. Set to `0` for the legacy filters. |
| `-taxonAwareMaxLoci` | integer | `500` | advanced | Maximum loci retained in the final concatenated alignment. |
| `-taxonAwareCoreLoci` | integer | `400` | advanced | Highest-scoring robust loci selected before greedy taxon-coverage rescue. Must not exceed `-taxonAwareMaxLoci`. |
| `-taxonAwareCandidateExtra` | integer | `150` | advanced | Extra pre-MSA candidate loci available to backfill alignment failures and MSA-QC rejections. |
| `-taxonAwareMinSequenceNT` | integer | `60` | advanced | Minimum unambiguous nucleotide-equivalent sequence length for an occurrence in the taxon-aware candidate pass. |
| `-taxonAwareTargetLoci` | integer | `25` | advanced | Per-sample locus target used by greedy coverage rescue and backbone-candidate reporting. |
| `-taxonAwareTargetNT` | integer | `7500` | advanced | Per-sample informative-NT target used by greedy coverage rescue and backbone-candidate reporting. |
| `-taxonAwareRescueMinPrevalence` | float | `0.8` | advanced | Minimum fraction of usable taxa a locus must occur in before taxon rescue or QC backfill may select it. `0` restores an unrestricted rescue pool; robust-core ranking is unaffected. |
| `-taxonAwareRescuePrevalenceMode` | string | `relative` | advanced | Interpret `-taxonAwareRescueMinPrevalence` relative to taxa with any usable locus (`relative`) or to all taxa (`absolute`). |
| `-taxonAwarePresortWeight` | float | `0.15` | advanced | Weight of the pre-MSA presort score in the final locus ranking. |
| `-taxonAwareTargetsFromGate` | integer | `1` | advanced | Derive the per-sample locus/NT targets from the active inclusion gate instead of `-taxonAwareTargetLoci`/`-taxonAwareTargetNT`. |
| `-taxonAwareInformationSaturation` | float | `0.005` within; `0.02` between | advanced | Variable-site fraction at which a locus is scored as fully informative. |
| `-taxonAwareExcessVariationOnset` | float | `0.05` within; `0.20` between | advanced | Variable-site fraction above which a locus starts to be penalised as hypervariable. |
| `-taxonAwareExcessVariationSpan` | float | `0.10` within; `0.30` between | advanced | Width of the penalty ramp above `-taxonAwareExcessVariationOnset`. |
| `-preferredCoreGenes` | string | empty direct; auto companion `.core` in `strain_within.pl` | advanced | Prefer loci whose catalogue seed occurs in this universal-core guide. The guide may use raw `.core` rows or sorted `MGS<TAB>gene,gene` rows. It does not make a locus eligible when sequence/QC rules reject it. |
| `-compactTaxonAwareDiagnostics` | integer | `1` | stable | After a successful or valid terminal run, merge taxon-aware and rate-partition audit TSVs into `phylo/taxon_aware_diagnostics.tsv`; set `0` to retain the individual source files. |
| `-runIQtree` | integer | `0` | stable | Build the tree with IQ-TREE. |
| `-AutoModel` | integer | `0` | stable | IQ-TREE model selection (`MFP+MERGE` for partitioned alignments). Off by default, which uses fixed `GTR+F+G2` / `LG+F+G`. |
| `-iqFast` | integer | `1` | stable | fast IQ-TREE mode |
| `-iqMemMB` | integer | `0` | advanced | IQ-TREE RAM cap in MB; `0` leaves IQ-TREE uncapped. |
| `-iqPathogen` | integer | `0` | advanced | Use the IQ-TREE 3 pathogen/CMAPLE path for low-divergence alignments; disables `-iqLegacy`. |
| `-iqLegacy` | integer | `0` | advanced | Restore the pre-5.14 IQ-TREE command. `strain_within.pl` enables this by default; pass `0` for the modern command. |
| `-runClonalFrameML` | integer | `0` | stable | Run ClonalFrameML recombination detection on the finished tree. |
| `-runGubbins` | integer | `0` | stable | Run Gubbins recombination detection (dormant; the binary is resolved lazily). |
| `-runLengthCheck` | integer | `1` | stable | check that sequence length can be divided by 3 |
| `-runDNDS` | integer | `0` | stable | run dNdS analysis |
| `-runTheta` | integer | `0` | stable | Calculate population-genetic theta estimates. |
| `-genesForDNDS` | string |  | stable | list with selected genes just for dnds |
| `-DNDSonSubset` | integer | `0` | stable | run dnds just on subset (given by genesForDNDS) of genes |
| `-codemlRepeats` | integer | `2` | stable | Independent codeml restarts per dN/dS model fit. |
| `-outDCodeml` | string |  | stable | Output directory for codeml/dN/dS results; empty writes under `-outD`. |
| `-genesToPhylip` | integer | `0` | stable | Also write the per-locus alignments in PHYLIP format. |
| `-runFastgear` | integer | `0` | stable | Run fastGEAR recombination analysis. |
| `-runFastGearPostProcessing` | integer | `0` | stable | Summarise existing fastGEAR output. |
| `-map` | string |  | stable | MATAFILER mapping file, used to relabel tips with sample metadata. |
| `-clustername` | string |  | stable | Label used in output file names and reports. |

Completed IQ-TREE runs are accepted only when the log contains its completion signature and the inferred backbone contains exactly the backbone-alignment taxa. IQ-TREE's standard identical-sequence handling is retained; `-keep-ident` is not used. The safe likelihood kernel is enabled pre-emptively for alignments with at least 750 taxa or over 700 MB; a numerical-underflow diagnostic also triggers one clean retry with `-safe`. Successful runs remove IQ-TREE transient files such as `.uniqueseq.phy`, while retaining the final tree, report, and log.

In strict-backbone IQ-TREE mode, `IQtree_allsites.backbone.treefile` contains the ML-inferred backbone and `IQtree_allsites.treefile` is the default primary tree containing accepted EPA-ng maximum-likelihood placements. For single-partition GTR models, BuildTree parses IQ-TREE's fitted exchangeabilities, base frequencies, invariant-site proportion when present, gamma-category count, and gamma shape, then supplies them directly as a RAxML-NG-style EPA-ng model descriptor. IQ-TREE 3 compact substitution-process tables are supported. EPA-ng accepts only one model for a concatenated placement alignment and does not refit symbolic-model parameters; with `+F` it calculates empirical base frequencies, while unfitted exchangeabilities and gamma shape remain generic defaults. Therefore, when the selected IQ-TREE backbone has multiple fitted GTR partition rows, BuildTree reuses IQ-TREE to estimate one unpartitioned `GTR+F+G2` model on the fixed retained backbone topology. It caches the validated refit as `IQtree_allsites.backbone.epa_model.*`, parses its explicit fitted descriptor for EPA-ng, and keeps the original backbone tree for placement and publication. RAxML-NG and legacy RAxML backbones retain their `.bestModel` and `.raxml.info` fitted-model reports. The raw `phylo/epa-ng/epa_result.jplace` output is retained. Before tree publication, placement pendant lengths are compared with the larger of the configured minimum (`-epaPendantMinThreshold`) and the non-outgroup backbone terminal-branch Q95 multiplied by `-epaPendantOutlierFactor`; separated long-branch placements are marked `excluded_outlier` and omitted from the primary tree. `phylo/strict_backbone.epa_placements.tsv` records the selected edge, likelihood, likelihood-weight ratio, branch lengths, applied cutoff, and `pendant_length_outlier` reason. This placement path requires an IQ-TREE or RAxML backbone with its model artifact; it deliberately does not fall back to an approximate nearest-tip graft. Taxon-aware runs retain every sample with any usable selected sequence through MSA. The established `-GenesPerSpecies`, `-relativeNTFraction`, and `-NTfiltCount` then decide backbone eligibility; `-strictBackboneFraction` remains a separate severe aligned-coverage deferral threshold. Mirrored `-placementGenesPerSpecies`, `-placementRelativeNTFraction`, and `-placementNTfiltCount` decide whether a deferred sample may be placed, together with the two-locus floor and `-placementMinOverlap`. Omitting a mirrored placement flag preserves backward compatibility by inheriting its backbone counterpart. `phylo/taxon_aware_backbone_eligibility.tsv`, `phylo/taxon_aware_placement_eligibility.tsv`, and `phylo/strict_backbone.samples.tsv` record each decision and exclusion reason.

Broad or between-species phylogeny is the default. In this mode `buildTree5.pl` does not remove columns using a fixed taxon-count overlap threshold, and post-alignment QC checks locus structure/occupancy without rejecting deep AA divergence. Use `-withinSpecies 1` for strain or other within-species trees. `-minOverlapMSA` and `-postAlignmentDivergenceQC` can override the individual mode defaults. Existing `-strainWithinPreset 1` calls remain compatible and imply within-species mode. With `-continue 1`, a stored QC-policy marker prevents reuse of a concatenated alignment built under different broad-mode settings; legacy within-species audits remain compatible.

All plain per-locus, concatenated, synonymous/nonsynonymous, full-backbone, placement, NEXUS, and partition files exist only in the job-specific scratch directory selected by `-tmp`/`TMPDIR`; tree programs consume those scratch paths directly. A newly created coding nucleotide locus runs MSAfix there immediately after back-translation and before its first gzip checkpoint is published. Subsequent mutable processing also remains on scratch. BuildTree atomically publishes only continuation-safe `.gz` artifacts under `MSA/`; a continued job expands them into fresh scratch without consuming them. BuildTree neither creates an uncompressed alignment in `MSA/` nor exposes scratch through a symlink there; legacy plain artifacts are migrated to gzip and removed.

Taxon-aware selection is enabled by default and stays inside `buildTree5.pl` so it can reuse the normal alignment and MSAfix stages. For direct calls, the pre-MSA pass ranks length-stable, complete, prevalent loci and retains a 400-locus robust core. If `-preferredCoreGenes` is supplied, a locus whose catalogue seed is listed in that GTDB/core guide is ranked ahead of other eligible loci and receives a modest secondary preference during greedy taxon rescue; it still must pass normal sequence and QC filters, and it does not relax the hard broad-availability gate for taxon rescue or QC backfill. Greedy taxon rescue and backfill remain restricted to loci found in at least 80% of taxa with any usable candidate locus, preventing a sparse accessory or mobile-element locus from entering only because it carries an underrepresented taxon. After MSAfix, the final pass recalculates prevalence and reapplies the same guard. Sparse samples remain available only when a broadly shared selected locus anchors them. Candidate and final rows expose `preferred_core`, `coverage_rescue_eligible`, and `coverage_rescue_reason`. By default, the final taxon-aware and rate-partition reports are sectioned into `phylo/taxon_aware_diagnostics.tsv`; `phylo/selection_attrition.tsv` remains separate because the strain controller aggregates it. Set `-compactTaxonAwareDiagnostics 0` for the individual files.

`strain_within.pl` enables the selector but leaves backbone placement disabled by default; set `-placeOnBackbone 1` to opt in. It forwards the 0.8 rescue-prevalence guard, and automatically uses the supplied `-MGS` guide when it ends in `.core` (otherwise an available sibling `-MGS.core`) as `-preferredCoreGenes`. Give `-preferredCoreGenes` explicitly to override that choice. It scales the hierarchy to its effective selected-gene budget: 80% robust core, 20% taxon-rescue capacity, and an additional 30% QC-backfill candidate pool. The final-tree budget is `min(treeLocusBudget, presortGenes)`; `maxGenes` separately limits how many validated loci are extracted per MGS/sample. Thus 500 selected genes produce 400 core + 100 rescue + 150 backfill, while the default 400-gene cap produces 320 + 80 + 120. Its filtering defaults are `-GeneLengthMin 0.3` for QC, `-GeneLengthIncludeMin 0.03` for post-QC MSA recovery, ordinary-tree `-GenesPerSpecies 0.2` and `-relativeNTFraction 0.1`, and—only when placement is enabled—thresholds of `0.04` loci and `0.03` NT plus a 10,000-site shared-backbone floor. Within-species locus QC rejects cross-locus consensus-divergence outliers above modified-Z 5.0. The strain preset now uses the legacy IQ-TREE command by default; pass `-iqLegacy 0` to use the standard modern command, or `-iqPathogen 1` to select the IQ-TREE 3 pathogen path (which automatically disables the legacy default).

For split Phase I (`-maxSubJob > 0`), assembly groups remain indivisible, but assignment is now weighted by estimated work rather than raw sample count: each sample contributes a fixed extraction cost plus its estimated uncompressed consensus FNA/FAA size, while samples requiring consensus regeneration receive an additional VCF-size penalty. Automatic splitting (`-maxSubJob -1`) counts standalone samples as effective groups. The consolidated run-level audit `LOGandSUB/phase1_worker_plan.tsv` records every group assignment, input state, estimated work, and worker total; a dominant group visible there can still create an unavoidable straggler because it cannot safely be split across workers.

Phase II first resolves one authoritative outgroup per actionable MGS in `strain_within.pl`, then derives the core-first demand manifest before loading candidate mappings or reference FASTAs. When `-MGSphylo` is supplied, it writes the consolidated Mosaic target-to-outgroup proposals to one temporary file and asks `neighborTree.R --all --preferred --max-candidates 1` for all tree-tip decisions in a single process, rather than starting R and rereading the source tree for each MGS. R accepts a plausible Mosaic proposal as the selection; if that proposal is absent, effectively identical, non-finite, or beyond the configured distance cutoff, the nearest eligible tree neighbour becomes the single selection. Perl does not reinsert the rejected proposal or try another candidate when the selected outgroup lacks coverage. Direct `neighborTree.R` calls can use `--preferred-tip TIP`, adjust `--preferred-quantile` and `--preferred-nearest-factor`, and optionally set `--max-candidates INT`; its default of `0` preserves an unlimited result for other callers. The temporary preference file is automatically removed and creates no per-MGS output clutter. The demand manifest uses GTDB/core-guide seed loci first, permits an exact Mosaic homolog only when its source is the selected outgroup and the locus is seed-listed or broadly available, and fills any shortfall only from COGs present in at least `-taxonAwareRescueMinPrevalence` of actionable targets. The default floor is `ceil(0.20 * -treeLocusBudget)` (80 with the default budget of 400); set `-outgroupCoreMinLoci` to override it. Targets below this floor are removed before the outgroup gene map is loaded; if none remain, Phase II opens neither candidate reference mappings nor FNA/FAA catalogues. For viable targets, `readGene2tax` retains only approved COGs for the selected outgroup before applying the per-outgroup `-outgroupReferenceGeneCap` (2,500 by default), and the sequential FASTA reader retains only those exact requested genes. Alternative same-COG copies are considered only inside that predetermined outgroup. `buildTree5.pl` receives the already finalized overlay and selected `-outgroup`; it performs no outgroup preselection.

Phase II resource planning adds no sample/gene pre-scan. The required `addOutgroup2MGS` input-finalization result already supplies submitted samples and usable genes. Ordinary core requests continue to use `ceil(sqrt(samples))` with the four-core floor and `-maxCores` cap. The controller uses `samples * usable_genes` as a cheap workload proxy: memory takes the larger of the existing input-size estimate and `samples * usable_genes / 1024` MB before applying the existing engine multiplier, the `-treeMemThreadDivisor` thread factor, `-phyloMemMulti`, and clamps. Once all ordinary jobs are prepared, submission is ordered by descending requested cores and then descending workload; EPA-only recovery retains its explicit priority tier.

The strain wrapper forwards backbone and EPA placement controls only when `-placeOnBackbone 1` is enabled. In that mode, well-covered samples infer the ML backbone; only severe coverage outliers or samples failing backbone admission are considered for EPA-ng placement, samples below the separate placement coverage/overlap gates remain excluded, and placements forming a clearly separated long-branch tail are omitted from the published tree. The default `-placeOnBackbone 0` infers one tree from the complete retained alignment and deactivates all backbone- and placement-only filters.

An ordinary parent `-onlySubmit 1` run uses a lean streaming resume path. It trusts the atomic `merge.complete.tsv` Phase-I commit, avoids the controller-wide per-MGS state audit and input-size/sort pass, and checks each MGS's completion/terminal markers, required category input, and missing resource estimate only when that MGS reaches the submission loop. Already-overlaid published or scratch inputs also bypass reference-catalogue initialization; the one shared reference stream begins only at the first raw MGS. The resulting tree job is queued immediately before the controller advances to the next MGS. A prior `LOGandSUB/tree_input_sizing.tsv` is optional and is read only as a resource hint; a legacy worker-shard contributor index is likewise deferred unless a checkpoint-less MGS actually needs that fallback. Explicit `-redo tree`, `-redo input`, `-redo all`, and `-redoEPAfilter` requests retain the exhaustive audit because those modes must classify or mutate existing state. In strict modes, published tree input is reusable only when its FNA, FAA, and category files are all present and nonempty; incomplete input is queued for recovery rather than treated as zero-size.

Deterministic merging runs after MSA QC and final taxon-aware locus selection, so it never changes the selected loci or taxa. Its initial partition count is `ceil(total effective called sites / -rateMergeTargetSites)` (30,000 by default), capped by `-rateMergeMaxBins`. It then repeatedly splits the largest current robust/backfill bin near its effective-site midpoint, using the locally more heterogeneous of the native QC report's P90 consensus divergence and final gap-free GC fraction. Thus additional partitions progressively resolve whichever P90 or GC signal remains, without allowing one very long locus to dominate a bin. Taxon-rescue loci join the nearest existing robust bin and cannot create sparse rescue-only partitions. If locus QC is disabled, the taxon-aware variable-site fraction supplies the rate proxy. Bins are collapsed into the nearest normalized rate/GC centroid until every remaining bin meets both minimum-size rules, or only one bin remains. The site rule uses effective called sites—the sum of each locus's mean called bases across retained taxa—so sparse missing-data loci contribute less. The retained `MSA/MSAli.fna.partition.RAXML.gz` contains comma-separated, non-contiguous locus ranges and is restored to scratch for continuation; `phylo/rate_merged_partitions.tsv` records every locus, selection phase, coordinate, missing-data-aware site count, metric source, initial bin, and final partition. These settings are included in the continuation-policy fingerprint, so changing them rebuilds stale alignments and trees.
