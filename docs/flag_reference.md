<!-- Documentation navigation -->
[Home](../README.md) | [Quick start](quickstart.md) | [Installation](install.md) | [Configuration](configuration.md) | [Mapping files](mapping_files.md) | [Workflows](common_workflows.md) | [Outputs](outputs.md) | [Flag reference](flag_reference.md) | [FAQ](FAQ.md) | [Glossary](glossary.md)

---

# Flag reference

This page is validated against the uploaded Perl source files for `MATAF4.pl`, `geneCat.pl`, `MGS.pl` and `buildTree5.pl`.

| Script | Version in uploaded source | Role |
|---|---:|---|
| `MATAF4.pl` | `4.11` | Main sample-level pipeline: read detection, preprocessing, host filtering, assembly, mapping, binning, SNP/SV calling and read-based profiling. |
| `geneCat.pl` | `0.51` | Gene catalog construction and downstream gene-catalog annotation/MGS orchestration. |
| `MGS.pl` | `0.45` | MGS/MAG dereplication, abundance/taxonomy and optional strain workflow orchestration. |
| `buildTree5.pl` | `5.20` | Phylogenetic tree construction and related MSA/population-genetic analyses. |

## How to read the tables

- **Aliases** are equivalent command-line names. Use any listed alias with a leading dash.
- **Type** is the argument type parsed by `Getopt::Long`: integer, float, string, string list or flag.
- **Default** is inferred from source-code assignments where possible. Empty cells mean no simple default could be inferred automatically.
- **Status** marks options that look stable, advanced/internal, deprecated/legacy or experimental/unsupported based on source comments and option names.

## Important validation notes

- `geneCat.pl` does **not** accept `-Binner`; use `-binSpeciesMG` for gene-catalog/MGS binning selection.
- `MATAF4.pl` accepts `-Binner`, `-MetaBat2` and `-binSpeciesMG` as aliases for the sample-level binning option.
- `MATAF4.pl` accepts `-profileMetaphlan`, not `-profileMetaphlan3`.
- `buildTree5.pl` accepts `-aa` for the amino-acid FASTA input; the older comment spelling `-faa` is not a parsed flag.
- `geneCat.pl -MGset` and `MGS.pl -MGset` are constrained in source to `GTDB` or `FMG`.

## MATAF4.pl

Main sample-level pipeline. This section preserves the more complete MATAF4.pl default annotations from the previous source-derived reference and updates the overall reference to include all uploaded Perl entry points.

## Base options

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-help`, `-?`, `-h` | flag | `` | stable | Show help. |
| `-checkInstall` | flag | `` | stable | Check that core MATAFILER4 programs/environments are installed. |
| `-map` | string | `""` | stable | Mapping file describing samples and input/output paths. |
| `-config` | string | `""` | stable | Alternative configuration file. |
| `-inspectState` | integer | `0` | stable | Emit a read-only JSON snapshot of workflow artifacts and markers. |
| `-planState` | integer | `0` | stable | Emit a read-only, dependency-ordered repair/submission plan from the inspection snapshot. |
| `-stateReport` | string | `""` | stable | Write the inspection JSON to this explicit path. |
| `-planReport` | string | `""` | stable | Write the repair/submission plan JSON to this explicit path. |
| `-autoStatePlan` | integer | `0` | advanced | Opt in to the internal inspect/plan preflight before normal execution and at each `loopTillComplete` boundary. Disabled by default to avoid a full metadata scan of every sample. |
| `-autoRepairState` | integer | `1` | stable | Apply only preflight repairs classified as automatically safe when submission is enabled. |

## Flow related

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-redoFails` | integer | `0` | stable | Remove failed sample-local results and rebuild the sample on the next pass; shared assembly-group outputs are retained. |
| `-redoContigStats` | integer | `0` | stable | runContigStats (coverage per gene, kmers, GC content) will be deleted & started again |
| `-submSystem` | string | `""` | stable | qsub,SGE,bsub,LSF.. by default will try to autodetect |
| `-submit` | integer | `1` | stable | submit any jobs at all? (0= no submission, just for trying if everything is correctly set up) |
| `-from` | integer | `0` | stable | start at which samples from map file? |
| `-to` | integer | `999999999999` | stable | stop at which samples from map file? |
| `-loopTillComplete` | string | `"0"` | advanced | Loop over selected samples; `X:Y` processes blocks of at most `Y` samples for up to `X` passes. One following block is admitted before waiting when a pass submits at most `floor(Y/4)` jobs (minimum 1), and always on the final allowed pass of a block. At most one block is added per pass and the enlarged range receives a fresh pass budget. Preflight is repeated after each wait only when `-autoStatePlan 1` is enabled. |
| `-loopTillCompleteActiveJobs` | integer | `3` | advanced | Start the next loop pass once no more than this many dependencies submitted for the current loop window are actually executing. Queued dependency-pending jobs and unrelated user jobs do not inflate the active count. |
| `-schedulerPollSeconds` | integer | `20` | advanced | Seconds between scheduler queries while `loopTillComplete` waits. Values must be positive. |
| `-excludeNodes` | string | `""` | stable | exclude certain nodes? |
| `-maxConcurrentJobs` | integer | `0` | stable | max jobs in queue, useful for large samples sets, currently only works on slurm |
| `-killDepNever` | integer | `0` | stable | kill jobs in "Dependency never finished" state? |
| `-requireInput` | integer | `0` | stable | in case input reads are no longer present, 0 will continue pipeline, 1 will abort |
| `-ignoreSmpls` | string | `""` | stable | Comma-separated exact sample IDs to skip; values are not regular expressions or prefix matches. |
| `-rmSmplLocks` | integer | `0` | stable | Remove existing sample locks. With the default `0`, a no-submission `loopTillComplete` pass reruns its current range when at least one user job remains active and the count is at most `-loopTillCompleteActiveJobs` or strictly below 1% of the range's samples. Retries use the normal pass budget plus at most one final extra scan. |
| `-silent` | flag | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-maxUnzpJobs` | integer | `20` | stable | how many unzip jobs to run in parallel (not to overload HPC IO). Default:20 |
| `-skipSmallSmplsMB` | integer | `1` | stable | skip samples with a combined input smaller than this in MB (raw file size, independent of compressed or raw) |
| `-forceWriteStats` | integer | `0` | stable | force (re)writing of the metagStats report and text file |

## File structure

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-rm_tmpdir_reads` | integer | `1` | stable | Default 1, remove tmpdir with reads |
| `-rm_tmpInput` | integer | `1` | stable | remove raw, human / adaptor filtered reads, if sdm clean created? (and not needed any longer) |
| `-reduceScratchUse` | integer | `1` | internal/advanced | remove sample scratch and rebuildable indexes only after ContigStats and configured binning/ConsSNP terminal outputs are published; set to 0 for debugging |
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
| `-inputBAMregex` | string | `""` | stable | Regex for detecting primary BAM read files under the location selected by the map's `Path` or `SmplPrefix`. Matching BAMs are treated as unpaired reads and converted with `samtools fastq`; for example, use `'.*\.bam$'`. Empty disables BAM discovery. |
| `-splitFastaInput` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-mergeReads` | integer | `0` | stable | merge read pair 1+2 before assembly etc? (usually doesn't help assembly, but useful for mapping to ref database in some rare instances) |
| `-ProbRdFilter` | integer | `1` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-pairedReadInput` | integer | `-1` | stable | determines if read pairs are expected in each in dir |
| `-inputReadLengthSuppl` | integer | `5000` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-filterHostRds`, `-filterHumanRds` | integer | `0` | stable | 0: no, 1:kraken2, 2: kraken1, 3:hostile |
| `-filterHostKrak2DB` | string | `""` | stable | customize host org to filter (e.g. human, chicken ..) |
| `-filterHostKr2Conf` | string | `0.01` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-filterHostKr2Quick` | string | `` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-hostileIndex` | string | `"human-t2t-hla"` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-onlyFilterZip` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-mocatFiltered` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-logQualvsLen` | integer | `0` | stable | sdm log file.. can be quite large; logs qual of read vs read length |

## Sdm related

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-inputReadLength` | integer | `150` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-gzipSDMout` | integer | `1` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-XfirstReads` | integer | `-1` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-minReadLength` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-maxReadLength` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-filterAdapters` | integer | `1` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-customSDMopt` | string | `""` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-sdmMem` | string | `"15G"` | stable | total mem for sdm job in Gb, default 15 |

## Assembly related

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-spadesCores`, `-assemblCores` | integer | `0` | stable | `0` automatically scales assembly jobs from 8 cores at up to 500 MiB of assembly-group input to 48 cores at 10 GiB or more. A positive value is an explicit override. |
| `-spadesMemory`, `-assemblMemory` | integer | `-1` | stable | in GB |
| `-spadesKmers`, `-assemblyKmers` | string | `"27,43,67,87,101,127"` | stable | comma delimited list |
| `-reAssembleMG` | integer | `0` | stable | Rebuild an assembly; shared assembly groups additionally require `-OKtoRWassGrps 1`. |
| `-asssemblyHddSpace` | integer | `"-1"` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-assembleMG` | integer | `0` | stable | 1=Spades, 2=MegaHIT, 3= flye, 4=metaMDBG, 5=hybrid ill-PB (megahit, metaMDBG) |
| `-assemblyLongTime` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-assemblyScaffMinSize` | integer | `500` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |

## Binning

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-Binner`, `-MetaBat2`, `-binSpeciesMG` | integer | `0` | stable | 0=no, 1=metaBat2, 2=SemiBin, 3: MetaDecoder, 4: GenomeFace, 5: SCGBinner |
| `-BinnerCores` | integer | `9` | stable | cores used for Binning process (and checkM) |
| `-BinnerMem` | integer | `0` | stable | define binning memory, Gb, 0=auto |
| `-minBinnerAssemblyMB` | float | `5` | stable | Do not run the binner when the assembly contains fewer than this many million sequence bases. Publish the standard empty binner outputs instead. Set to `0` to disable the cutoff. |
| `-checkM2` | integer | `1` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-checkM1` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-BinnerScratchTmp` | integer | `0` | internal/advanced | very specific (undocumented) use of scratch instead of nodetmp dir |
| `-redoEmptyBins` | integer | `0` | internal/advanced | debug option; redo bins that are empty (no bin detected). Note: this can sometimes happen for metagenomes |
| `-redoBinning` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-SB_env` | string | `""` | stable | semiBin environment; if given, will avoid re-training de novo binning model. Default: "" (autotrain). should be #human_gut/dog_gut/ocean/soil/cat_gut/human_oral/mouse_gut/pig_gut/built_environment/wastewater/chicken_caecum/global |

## Gene prediction on assembly

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-predictEukGenes` | integer | `0` | stable | severely limits total predicted gene amount (~25% of total genes) |
| `-kmerPerGene` | integer | `0` | stable | calculate kmer frequencies for each gene instead of per scaffold |
| `-genePredGZenforce` | integer | `1` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-rewriteGenePred` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |

## Mapping

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-mapper` | integer | `-1` | stable | #1=bowtie2, 2=bwa, 3=minimap2, 4=kma, 5=strobealign -1=auto (bowtie2 short, minimap2 long reads), -2=auto(strobealign short, minimap2 long) |
| `-mapUnmapped` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-mappingCoverage` | integer | `1` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-mappingMem` | integer | `-1` | stable | total mem for mini2/kma/bwa/bwt2 in GB |
| `-mapSortMem` | integer | `-1` | stable | total mem for samtools sort in GB |
| `-rmDuplicates` | integer | `1` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-mappingCores` | integer | `8` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-mapperFilterIll` | string | `"0.05 0.75 20 3"` | stable | Maximum NM edit rate, minimum query coverage, minimum mapping quality, and minimum clipping at both ends (0 disables clipping filter). |
| `-mapperFilterHybridIll` | string | `"0.03 0.90 40 5"` | advanced | Stricter Illumina BAM filter used while deriving coverage for hybrid preassemblies: maximum edit rate, minimum query coverage, minimum mapping quality and minimum clipping at both ends. |
| `-hybridMinMapQ` | integer | `40` | advanced | Minimum mapping quality passed to `samtools depth` for hybrid-preassembly coverage. |
| `-hybridMinBaseQ` | integer | `20` | advanced | Minimum base quality passed to `samtools depth` for hybrid-preassembly coverage. |
| `-breakpointDepth` | float | `0.10` | advanced | Relative coverage threshold used to identify low-depth assembly breakpoints for hybrid read simulation. |
| `-breakpointMinLength` | integer | `100` | advanced | Minimum length of a low-depth region reported as a hybrid-assembly breakpoint. |
| `-breakpointSmoothGap` | integer | `100` | advanced | Maximum gap joined while smoothing adjacent low-depth breakpoint regions. |
| `-breakpointFlankLength` | integer | `500` | advanced | Number of bases inspected on each side of a candidate breakpoint. |
| `-breakpointMinFlankDepth` | float | `1` | advanced | Minimum flank depth required when accepting a candidate breakpoint. |
| `-breakpointMaxFlankFraction` | float | `0.10` | advanced | Maximum low-depth fraction allowed within breakpoint flanks. |
| `-hybridSyntheticMaxDepth` | float | `20` | advanced | Cap on synthetic read depth generated from each hybrid-preassembly package. |
| `-mapperFilterPB` | string | `"0.05 0.5 30 0"` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-mapperFilterONT` | string | `"0.15 0.5 10 0"` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-mapSaveCRAM` | integer | `1` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |

## Mapping related (2) (assembly)

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-remap2assembly`, `-redoMap2assembly`, `-redoMapping` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-JGIdepths` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-mapReadsOntoAssembly` | integer | `1` | stable | map original reads back on assembly, to estimate abundance etc |
| `-mapSupportReadsOntoAssembly` | integer | `1` | stable | Map `SupportReads` onto the assembly and calculate their coverage separately. |
| `-saveReadsNotMap2Assembly` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |

## Map2tar / map2db / map2gc

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-decoyMapping` | integer | `1` | stable | 1: "Decoy mapping": map against reference genome AND against assembly of metagenome (drawing obvious better hits to metagenome, the "decoy") |
| `-competitive2ndmap` | integer | `-1` | stable | 1: Competitive, 2: combined but report separately per input genome, -1: combined and report all together |
| `-ref` | string | `""` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-mapperLargeRef` | integer | `0` | stable | use flags in mapper index built for large ref DBs? |
| `-mapnms` | string | `""` | stable | name for this final files |
| `-redo2ndmap` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |

## Snps

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-get2ndMappingConsSNP` | integer | `0` | stable | SNPs (onto mapping) |
| `-getAssemblConsSNP` | integer | `0` | stable | SNPs (onto self assembly) #calculates consensus SNP of assembly (useful for checking assembly gets consensus and Assmbl_grps) |
| `-getAssemblConsSNPsuppRds` | integer | `0` | stable | same as getAssemblConsSNP, but SNP calling for support reads |
| `-redoAssmblConsSNP` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-SNPmem` | integer | `0` | stable | memory per assigned core, in GB |
| `-redoGeneExtrSNP` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-SNPjobSsplit` | integer | `0` | stable | parallel jobs per sample; 0 uses tiered alignment-size estimates from 2 cores at 300 MB to 10 cores at 10 GB |
| `-SNPminCallQual` | integer | `20` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-SNPsaveVCF` | integer | `1` | stable | save vcf of SNP calles? DEfault : 1 |
| `-SNPsaveConsFasta` | integer | `0` | stable | Save consensus fasta from vcf calls? Default: 0 -> too large, can be quickly recreated.. |
| `-SNPcaller` | string | `"MPI"` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-SNPcores` | integer | `10` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-SNPconsMinDepth` | integer | `0` | stable | how many reads coverage to include position for consensus call? |
| `-SNPnormINDEL` | integer | `1` | stable | using bcftools norm to left-align indels |
| `-SVcaller` | integer | `0` | stable | calling structural variants: 1=delly, 2=gridss. Default (0). |

## Functional profiling (diamond)

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-profileFunct` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-reParseFunct` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-reProfileFunct` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-reProfileFuncTogether` | integer | `0` | stable | if any func database needs to be redone, than redo all indicated databases (useful if number of reads used changes..) |
| `-diamondCores` | integer | `12` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-diamondMem` | integer | `7` | stable | memory in GB for diamond alignment jobs |
| `-DiaParseEvals` | string | `"1e-7"` | stable | evalues at which to accept hits to func database |
| `-DiaSensitiveMode` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-rmRawDiamondHits` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-DiaMinAlignLen` | integer | `20` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-DiaMinFracQueryCov` | float | `0.1` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-DiaPercID` | integer | `40` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-diamondDBs` | string | `""` | stable | NOG,MOH,ABR,ABRc,ACL,KGM,CZy,PTV,PAB,MOH2,URE,URacc,AMI |

## Functional profiling (jaime tree)

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-orthoExtract` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |

## Ribo profiling (mitag)

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-profileRibosome` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-riobsomalAssembly` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-reProfileRibosome` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-reRibosomeLCA` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-riboMaxRds` | integer | `250000` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-saveRiboRds` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-thoroughCheckRiboFinish` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |

## Other tax profilers..

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-profileMetaphlan` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-profileMOTU2` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-profileKraken` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-profileTaxaTarget` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |
| `-estGenoSize` | integer | `0` | stable | estimate average size of genomes in data |
| `-krakenDB` | string | `""` | stable | "virusDB";#= "minikraken_2015/"; |

## D2s distance

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-calcInterMGdistance` | integer | `0` | stable | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |

## Io for specific uses

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-newFileStructure` | string | `""` | internal/advanced | just relink raw files for use in mocat |
| `-upload2EBI` | string | `""` | internal/advanced | copy human read removed raw files to this dir, named after sample |

## Institute specific: ei

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-wcKeyJobs` | string | `""` | internal/advanced | Accepted by MATAF4.pl; inspect source or help output for detailed behaviour. |

## Debug

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-OKtoRWassGrps` | integer | `0` | internal/advanced | can delete assemblies, if suspects error in them |

## Flag comparison against previous manual.md

This comparison was produced by matching dash-prefixed names from the previous long-form manual against the flags parsed from `MATAF4.pl`. It is intentionally conservative because prose can contain option-like strings.

- Flags accepted by `MATAF4.pl`: **155**
- Option-like names previously seen in the legacy manual: **181**
- Accepted by code but not previously documented: **26**
- Mentioned in the legacy manual but not accepted directly by `MATAF4.pl`: **52**

### Accepted by `MATAF4.pl` but not previously documented

`-?`, `-BinnerScratchTmp`, `-MetaBat2`, `-SNPnormINDEL`, `-assemblyLongTime`, `-assemblyScaffMinSize`, `-binSpeciesMG`, `-checkInstall`, `-customSDMopt`, `-filterAdapters`, `-filterHumanRds`, `-getAssemblConsSNPsuppRds`, `-h`, `-help`, `-inputReadLengthSuppl`, `-logQualvsLen`, `-mapperFilterONT`, `-mapperFilterPB`, `-orthoExtract`, `-profileMetaphlan`, `-profileTaxaTarget`, `-reAssembleMG`, `-redoMap2assembly`, `-redoMapping`, `-rewriteGenePred`, `-spadesKmers`

### Previously mentioned but not accepted as direct `MATAF4.pl` flags

These may be legacy prose artefacts, GeneCat/MGS options, flags for other scripts, or removed options. Do not use them with `MATAF4.pl` unless they are documented elsewhere.

`-AutoModel`, `-DNDSonSubset`, `-GCd`, `-GenesPerSpecies`, `-MGset`, `-MSAprogram`, `-NTfilt`, `-NTfiltCount`, `-NTfiltPerGene`, `-NonSynTree`, `-SNPmemPerJob`, `-SynTree`, `-aa`, `-bootstrap`, `-calcDiffDNA`, `-calcDistMat`, `-calcDistMatExt`, `-cats`, `-clustername`, `-codemlRepeats`, `-fixHeaders`, `-fna`, `-fracMaxGenes90pct`, `-genesForDNDS`, `-genesToPhylip`, `-gzInput`, `-isAligned`, `-maxGapPerCol`, `-minOverlapMSA`, `-minPcId`, `-outDCodeml`, `-outgroup`, `-postFilter`, `-profileMetaphlan2`, `-quick`, `-rmMSA`, `-runClonalFrameML`, `-runDNDS`, `-runFastGearPostProcessing`, `-runFastTree`, `-runFastgear`, `-runGubbins`, `-runLengthCheck`, `-runRAxML`, `-runRaxMLng`, `-runTheta`, `-subsetSmpls`, `-superCheck`, `-superTree`, `-tmp`, `-useEte`, `-useTrimomatic`

## geneCat.pl

Gene-catalog construction and downstream gene-catalog annotation/MGS orchestration. The uploaded source reports version `0.51`.

### Directories/files

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-o`, `-GCd` | string |  | stable | main save location for gene catalog and supporting files |
| `-tmp` | string | `$GLBtmp` | stable | tmp dir, global availalbe |
| `-glbTmp` | string |  | stable | global tmp dir, same as -tmp usually |
| `-map` | string | `?` | stable | mapping file, can be a combination of several .map files to combine different datasets (e.g. -map file1.map,f2.map) |

### Run modes

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-m`, `-mode` | string | `geneCat` | stable | possible modes: mergeCLs CANOPY specI kraken kaiju FMG_extr FOAM ABR FuncAssign protExtract ntMatchGC geneCat subprepSmpls |

### Cluster options

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-clusterID` | integer | `95` | stable | identity at which to cluster gene catalog, default: 0.95 |
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
| `-submSystem` | string |  | stable | Accepted by geneCat.pl; see source/help output for detailed behaviour. |
| `-continue`, `-justCDhit` | integer | `1` | stable | flow control, 1: continue with found files 0: delete existing (partial) gene cat, start again |
| `-c`, `-cores` | integer | `20` | stable | Accepted by geneCat.pl; see source/help output for detailed behaviour. |
| `-c0`, `-cores0` | integer | `-1` | stable | specifcally cores only for the big main clustering job.. |
| `-c3`, `-cores3` | integer | `-1` | stable | for small jobs that really don't require that much power.. |
| `-mem` | integer | `200` | stable | max mem |
| `-stone` | string |  | stable | Accepted by geneCat.pl; see source/help output for detailed behaviour. |
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
| `-FuncMinBitSc` | float | `45` | stable | Accepted by geneCat.pl; see source/help output for detailed behaviour. |
| `-FuncMinAlLeng` | integer | `30` | stable | Accepted by geneCat.pl; see source/help output for detailed behaviour. |
| `-FuncMinPercSbjCov` | float | `0.5` | stable | Minimum fraction of subject coverage for functional assignment. |
| `-FuncMinPerID` | float | `25` | stable | Accepted by geneCat.pl; see source/help output for detailed behaviour. |
| `-FuncMinEVal` | float | `1e-8` | stable | Accepted by geneCat.pl; see source/help output for detailed behaviour. |

## MGS.pl

MGS/MAG dereplication, abundance/taxonomy and optional strain workflow orchestration. The uploaded source reports version `0.32`.

### General options

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-GCd` | string |  | stable | gene catalog dir |
| `-clusterID` | integer | `95` | stable | gene-catalog clustering identity percentage |
| `-outD` | string |  | stable | defaults to $inD/Bin_SB/ |
| `-tmp` | string |  | stable | temp dir |
| `-submit` | integer | `1` | stable | 1:submit jobs, 0: dry run. Default: 1 |
| `-canopies` | string |  | stable | location of canopy clustering output file (clusters.txt) |
| `-smallCores` | integer | `4` | stable | cores used for normal jobs (not intensive) |
| `-bottleneckCores` | integer | `12` | stable | cores for compute intensive jobs |
| `-redoCluster` | integer | `0` | stable | Accepted by MGS.pl; see source/help output for detailed behaviour. |
| `-redoTax` | integer | `0` | stable | rewrite tax annotations |
| `-MGset` | string | `FMG` | stable | GTDB or FMG, which marker genes are used? Default: GTDB |
| `-wait4stone` | string |  | stable | wait for these files to be created, refers currently exclusively to eggNOG annotations that are needed later |
| `-wait4stoneTimeout` | integer | `86400` | stable | maximum wait in seconds; zero waits indefinitely |
| `-mem` | integer | `150` | stable | memory used for intensive jobs |
| `-strains` | integer | `0` | stable | 1: calc instra species strain phylogenies. Default: 0 |
| `-prepareMosaicLoci` | integer | `1` | stable | 1: have strain_within submit, await, and validate Mosaic as a prerequisite before strain analysis; 0: skip Mosaic preprocessing and keep same-NOG seed clusters separate |
| `-useCheckM2` | integer | `0` | stable | CheckM2 default qual checking of MAGs/MGS |
| `-useCheckM1` | integer | `1` | stable | CheckM default qual checking of MAGs/MGS |
| `-binSpeciesMG` | integer | `2` | stable | 0=no, 1=metaBat2, 2=SemiBin, 3: MetaDecoder, 4 ,5 |
| `-ignoreIncompleteMAGs` | integer | `1` | stable | 1: assemblies without MAG calculations are ignored. Default: 1 |
| `-legacy` | integer | `0` | deprecated/legacy | 1: use legacy code as pre Dec `22 (clustering is a bit more muddy, reported abundances slightly different, remember to use -MGset FMG). No longer supported. Default: 0 |
| `-genomesPerFamily` | integer | `0` | stable | Accepted by MGS.pl; see source/help output for detailed behaviour. |

## buildTree5.pl

Phylogenetic tree construction and related MSA/population-genetic analyses. The source reports version `5.20`.

### General options

| Aliases | Type | Default | Status | Description |
|---|---:|---|---|---|
| `-genoInD` | string |  | stable | provide a dir with complete genomes, will extract FGMs and build tree between genomes (NT/AA flag) |
| `-wildcardflag` | string |  | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-fna` | string |  | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-aa` | string |  | stable | Amino-acid FASTA input. Note: source flag is -aa, not -faa. |
| `-cats` | string |  | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-outD` | string |  | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-tmpD` | string |  | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-cores` | integer | `1` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-superTree` | integer | `0` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-superCheck` | integer | `0` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-fixHeaders` | integer | `0` | stable | # fix the fasta headers, if too long or containing not allowed symbols (nwk reserved) |
| `-useEte` | integer | `0` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-NTfilt` | float | `0.8` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-NTfiltPerGene` | float | `0.1` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-GenesPerSpecies` | float | `0.1` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-fracMaxGenes90pct` | float | `0.25` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-NTfiltCount` | integer | `0` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-smplDef` | integer | `1` | stable | is the genome somehow quantified with a delimiter (_) ? |
| `-smplSep` | string | `_` | stable | set the delimiter |
| `-outgroup` | string |  | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-AAtree` | integer | `0` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-MSAprogram` | integer | `2` | stable | (0) MSAprobs, (1) clustalO, (2) mafft, (4) MUSCLE5, (5) FAMSA2 (only AA) |
| `-minOverlapMSA` | integer | `0` | stable | min overlap in MSA columns, in order to retain column |
| `-maxGapPerCol` | float | `1` | stable | same as minOverlapMSA, but for MSAfix and %of gaps allowed in a column |
| `-calcDistMat` | integer | `0` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-calcDistMatExt` | integer | `0` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-calcDiffDNA` | integer | `0` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-minPcId` | float | `0` | stable | sequence is filtered from data, unless the average minPcId is >= $minPcId |
| `-SynTree` | integer | `0` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-NonSynTree` | integer | `0` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-continue` | integer | `0` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-bootstrap` | integer | `0` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-subsetSmpls` | integer | `-1` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-postFilter` | string |  | stable | "," sep list of zorro,guidance2,macse |
| `-rmMSA` | integer | `0` | stable | to save diskspace |
| `-gzInput` | integer | `0` | stable | to save diskspace |
| `-isAligned` | integer | `0` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-runRAxML` | integer | `0` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-runRaxMLng` | integer | `0` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-runFastTree` | integer | `0` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-runVeryFastTree` | integer | `0` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-treeShrink` | integer | `0` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-postAlignmentLocusQC` | integer | `1` | stable | Run native MSAfix locus-comparability QC before multi-locus concatenation. |
| `-postAlignmentMinSequences` | integer | `3` | stable | Minimum comparable sequences required for an aligned locus. |
| `-postAlignmentMinOccupancy` | float | `0.35` | stable | Minimum fraction of unambiguous alignment cells; permissive for metagenomic loci. |
| `-postAlignmentRelativeZ` | float | `8.0` | stable | Modified-Z threshold for cross-locus consensus-divergence outliers. |
| `-postAlignmentMinLociRelative` | integer | `8` | stable | Minimum locus count before applying cross-locus robust outlier QC. |
| `-runIQtree` | integer | `0` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-AutoModel` | integer | `1` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-iqFast` | integer | `0` | stable | fast qiTree mode |
| `-runClonalFrameML` | integer | `0` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-runGubbins` | integer | `0` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-runLengthCheck` | integer | `1` | stable | check that sequence length can be divided by 3 |
| `-runDNDS` | integer | `0` | stable | run dNdS analysis |
| `-runTheta` | integer | `0` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-genesForDNDS` | string |  | stable | list with selected genes just for dnds |
| `-DNDSonSubset` | integer | `0` | stable | run dnds just on subset (given by genesForDNDS) of genes |
| `-codemlRepeats` | integer | `2` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-outDCodeml` | string |  | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-genesToPhylip` | integer | `0` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-runFastgear` | integer | `0` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-runFastGearPostProcessing` | integer | `0` | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-map` | string |  | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
| `-clustername` | string |  | stable | Accepted by buildTree5.pl; see source/help output for detailed behaviour. |
