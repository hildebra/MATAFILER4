# MATAFILER4 pipeline audit — 24 September 2026

## Scope and method

The whole active pipeline was reviewed: `MATAF4.pl` (in four slices), the shared `Mods/`, the gene catalogue and functional annotation (`geneCat.pl`, `secScripts/functions/`), the MGS stage, `strain_within.pl`, `buildTree5.pl`, and the helper scripts the pipeline invokes. The two most recent, previously unaudited commits (`0d6cc12` SAM/BAM/CRAM streaming, `bcb6e84` gene-catalogue coverage clean-up) got particular attention. Findings already listed as fixed or deliberate in the earlier audits were excluded.

Every fix below was checked against the code and its callers before it was made. Most were reproduced beforehand with small Perl probes that run the real subroutine or script. A fix was applied only when the intended behaviour was clear from the code, its callers or the docs, and the change could not break existing projects. Everything else is listed under "Found, not fixed".

The pipeline was **not** run end to end: there was no Slurm cluster and none of the external tools (aligners, assemblers, CheckM2, …) were available. Verification consisted of the unit-test suite (run under WSL Ubuntu, perl 5.38), `perl -c` on every changed script, and targeted probes.

## Fixed

### Installation and configuration

| Location | Defect | Fix |
|---|---|---|
| `Mods/config.old` | The template the installer copies into `config.txt` used `MFLRDir $MGTKDIR` and `CONDAbaseEnv MGTK`. A fresh install died on the first config read (the installer only exports `MF4DIR`); with `MGTKDIR` set, `checkMF` exited 23 because the documented env is `MF4`. `DBDir` was the literal `empty`. | `$MF4DIR`, `MF4`, and `DBDir $MF4DIR/data/DBs/`, where the installer puts databases. Existing configs are not touched. |

### `MATAF4.pl`

| Location | Defect | Fix |
|---|---|---|
| `alignPostTreat` | `${$postTreat{doCram}}` dereferenced a plain number under `strict`: decoy/competitive secondary mapping (map2tar default) died while building every job. | pass the value |
| `mapReadsToRef` (KMA) | `-o $tmpOutxtra[$i]` indexed a per-reference array with the library index: a bare `-o` for every library after the first. | one prefix per reference and library |
| `alignPostTreat` | `samtools fastq` without `-0`: unaligned single-end reads went to the job log instead of `unal.fq.gz`. | add `-0` output |
| `seedUnzip2tmp` | Since `0d6cc12` `rawRds/` is no longer wiped, but only the first porechop output was truncated before `>>`: a requeue duplicated the others. | truncate every output |
| `Mods/WorkflowControl.pm` | `0d6cc12` quoted the prepended `rm`, so local lightweight staging was never recognised and became a scheduler job plus an extra loop pass. | accept a single-quoted program word (test added) |
| `unploadRawFilePostprocess` | Looked for `*.Rsingl.fq.gz`; uploads are written as `*.Rsingle.fq.gz`. Background md5sums were not awaited. | correct name, `wait` |
| `getRgStr` | minimap2/bwa read-group ID was the literal `$smpl`. | real sample name |
| `mergeReads` | Multi-library merging `cat` used paths relative to the job's working directory, not FLASH's `-d` directory: the job always failed. | absolute paths |
| `detectRibo` | Scratch check looked for a lambda 1.9 index that lambda3 never writes: SILVA was recopied on every run. | check `.lba.gz` |
| `prepDiamondDB` | TCDB check tested `hir.txt`, the file copied is `TCDBhir.txt`: a DB-prep job on every run. | correct name |
| `prepareDiamondRerun` | Regression from `579c34c`: `-reProfileFunct 1` only re-aligned when all six databases were requested. | delete the requested databases' hits again |
| main loop / `loop2C_check` | Closed samples only counted towards RiboFind, so mOTU/MetaPhlAn/Diamond merges went stale when samples were added. Counters also accumulated across `-loopTillComplete` passes, so failures from early passes blocked the final merges. | count closed samples from the completion sentinel's component evidence; reset counters when a new pass starts (the final full pass keeps them) |
| `reduceProgStats` | Kraken failure count was not undone for empty samples. | undo it (guarded) |
| option guard | `-rewriteGenePred`, `-redoContigStats`, `-redo2ndmap`, `-reParseFunct`, `-redoKraken` and `-redoFails` delete outputs on every pass, so `-loopTillComplete` could not converge (`-redoFails` deletes samples whose jobs are still running). | rejected with `-loopTillComplete`, like the other rewrite options; noted in `flag_reference.md` |
| R3 remap rule | With `-mapSaveCRAM 0` and binning, cleanup removes the CRAM on purpose, but R3 treated its absence as a broken mapping and wiped mapping, coverage and SNP output of every reopened sample. | CRAM required only when it is kept; SNP/SV work that needs it still remaps |
| unzip throttle | `EMPTY_DO_NEXT` and failed/deferred markers were used as scheduler dependencies of later staging jobs. | only real job IDs |
| `contigStatsOutputsComplete`, cleanup | Primary coverage was required even with `-mapReadsOntoAssembly 0`: ContigStats resubmitted forever, samples never closed. | same rule as `sampleCompletionComponents` |
| k-mer completion checks | `scaff.pergene.4kmer.pm5` is written to the assembly-group ContigStats but was looked for per sample: co-assembly members never completed. | assembly-group path (three checks) |
| `prepareMap` | Secondary-reference gene prediction always gzips `genes.gff`; the job copied the plain file and failed, blocking secondary mapping. | decompress (or move the plain file from empty input) |
| `spadesAssembly` | Auto memory gave SPAdes a fractional `-m`. | integer |
| `optiDups` | Read `map2.sh.etxt`; markdup has run inside `map.sh` since the jobs were merged. | first log containing duplicate stats |
| `sdmStatsMany` | One hard-coded `filter_lenHist.txt` was used for every sdm log. | each log's own histogram |
| `getSNPStats` | Regex did not match current vcf2fna output; `SNP_Passed`/`INDEL_Passed` were always 0. | current and old wording |
| `getContamination` | hostile's fractional `reads_removed_proportion` never matched. | parsed and reported as a percentage |
| `getGeneStats` | Only read plain `GeneStats.txt`; it is published gzipped, so gene columns were always blank. | read `.gz` too |
| option check | Error message named geneCat's `-useCheckM1/2`. | `-checkM1/2` |

### Shared modules

| Location | Defect | Fix |
|---|---|---|
| `Mods/SNP.pm` | Restarting consensus calling after a partial run: finished chunks kept the BED file the planner rewrites, and the final check exited 33 every time. | finished chunks remove their BED |
| `Mods/IO_Tamoc_progs.pm` `jgi_depth_cmd` | A BAM that the same job will create was treated as a sample directory: `-JGIdepths 1` crashed job construction. | only an existing empty file is an error |
| `Mods/Subm.pm` | SGE `h_rss` got a unitless MB number, i.e. bytes. | `M` suffix |
| `Mods/FuncTools.pm` | foldseek's default third column is a 0–1 fraction, compared with `-percID` as a percentage: no annotations. | BLAST-like `--format-output` |
| `Mods/GenoMetaAss.pm` | Exported a `runDiamond` it does not define. | removed |

### Gene catalogue and functional annotation

| Location | Defect | Fix |
|---|---|---|
| `secScripts/GC/eggNOG_split.sh` | `awk '$2~!/-/'` kept rows only if column 2 contained a `1`/`0`: most CAZy/KEGG/EC/GO/BiGG/PFAM annotations were lost. | `$2 != "-"`, header skipped |
| `geneCat.pl` collation | `int($n/$b*($i+1))` rounds down for many cohort sizes (e.g. 2008 samples, 11 batches): the last sample was never collated. | multiply before dividing |
| `geneCat.pl` FuncAssign | `redo=>1` for the first database deleted its finished DIAMOND results at every rerun. | no forced redo |
| `combine_DIA.pl` | Missing COG and category files were recorded in each other's failure sets; `J\t5` style rows (3 characters) were skipped. | sets swapped back; rows need a tab |

### MGS, strains and phylogeny

| Location | Defect | Fix |
|---|---|---|
| `MGS.pl` | Rebuilding MGS clusters left the Mosaic catalogue of the old clusters, which strain_within reuses verbatim. | invalidated with the other derivatives |
| `clusterMAGs.pl` (`-perlClusterMAGs`) | Contig separators were counted as a gene (MGS.pl then died); thresholds shifted a tier when a tier was empty or its last MAG was ambiguous. | skip separators; use each MAG's own tier |
| `strain_within.pl` | Resuming an output from before v1.54 (sorted guide beside the input guide) wiped the output directory. | fall back to the legacy guide in resume modes |
| `strain_within.pl` | The tree OOM supervisor exited without a final accounting scan when all jobs finished between scans. | one final scan |
| `strain_within.pl`, `buildTree5.pl` | `-redoEPAfilter 0` enabled the redo. | `:1` option |
| `buildTree5.pl` | Resuming an interrupted run deleted all finished per-locus alignments (the QC reports it checks only exist after the locus loop). | only with merged output present |
| `buildTree5.pl` | A sparse outgroup kept by the prefilter hit an internal assertion. | outgroup exempt |

### Other invoked scripts

| Location | Defect | Fix |
|---|---|---|
| `runBinners.pl` | A failed GPU (GenomeFace) job was accepted as "0 bins" permanently. | require the job's stone |
| `deployMapDB.pl` | Built the decoy index with mapper 0, which is rejected: decoy mapping always died. | bowtie2, as the caller expects |
| `ENASRAdl.pl` | In a three-file paired ENA run the unpaired file became read 1. | positional fallback only for two files |
| `MG_LCA.pl` | A leftover `.m8` from an interrupted job blocked all retries of that COG. | removed when its LCA output is missing |
| `checkFQhds4ENA.pl` | Rewrote uncompressed FASTQ as gzip under the `.fq` name. | keep the input's compression |

## Found, not fixed

These need a decision, cluster-level testing, or touch behaviour existing projects rely on.

- **CheckM2 database path.** `config_DBs.txt` says `[DBDir]/checkm2/…`; the installer and `docs/install.md` use `data/DBs/CM2/…`, and `runCheckM2` exports the config value as `CHECKM2DB`. Aligning either side can break sites set up the other way, and changing the config value invalidates MGS canopy checkpoints.
- **`-completeContaStats` is inert** since `ba6d4db` (compares with the old `"?\t"` return). Reviving it would re-run host filtering on existing projects.
- **Duplicate between-MGS tree job**: rerunning `MGS.pl` while the tree job is queued submits a second buildTree5 into the same directory (needs job-ID tracking).
- **strain_within split Phase I**: recovery is gated on `$dirsNOTPrepped` instead of unmerged ledgers; `-redo input` schedules the freshly regenerated scratch inputs for deletion; `selected_mgs` differs between parent and workers under `-redo tree`. Each needs Slurm-level testing.
- **geneCat FuncAssign**: the first database's collection job deletes the split FASTAs the other databases' DIAMOND jobs may still read.
- **`extrAllE100GC.pl`**: the marker stone is written before the matrix job finishes; `Mattrix.FMG.mat` typos disable its skip checks (fixing them enables code paths that never ran).
- **`kmerPerGene.pl`** dies in geneCat's `-submitLocal 0` mode when k-mers were not computed (geneCat adds the step unconditionally).
- **`-thoroughCheckRiboFinish 1`** deletes valid empty ribosomal hierarchies and loops (opt-in).
- **buildTree5 opt-in paths**: `-genoInD` with `_` in genome names, `-runDNDS` analysing no loci, single-locus `-SynTree/-NonSynTree` file names, `-calcDiffDNA` identity maths (`tr/[-]//` counts `[`..`]`), an absent outgroup still passed to RAxML, gap codons in `synPosOnly`, the ClonalFrameML glob, the `-gzInput` fingerprint, single-locus `.gz` input, NEXUS `datatype=dna` for proteins.
- **SAM/BAM/CRAM input (`0d6cc12`)**: READ1/READ2-flagged records go to `/dev/null`, so a paired unaligned BAM maps as empty while sdm cleans it fully. Documented as singleton-only, not enforced. The new `inputCramReference*` signature keys invalidate every existing completion sentinel once.
- **Python files** (not edited): `get_ranks.py` expects NCBI rank `superkingdom` (renamed `domain` in 2025); `get_gtdb.py` still reads `MGTKDIR`.
- **Smaller items**: legacy Kraken1 `krak_count_tax.pl` assigns ranks by position; `createPsAssLongReads` passes `-1` (dormant); `MATAF4.pl:~1662` terminal-empty removal regex never matches (making it match would break group-membership checks); accession samples re-download when their sentinel is rejected; `bin/distv9.pl` misreads scientific-notation branch lengths (not called by live code).

## Retired scripts removed

The 2026-09-23 changelog retired a set of scripts and removed their config keys, but the files were never deleted from git. They caused every failure of `t/config_integrity.t` and `t/script_compile.t`. They have now been removed (48 files). Each was checked to be unreferenced by live code, tests and docs.

- Named in that changelog: `secScripts/phylo/buildTree4.pl`, `secScripts/MGS/strain_within_3.pl`, the vcf2cons SNP chain (`secScripts/SNP/` except the live `plan_consensus_regions.pl`, plus `secScripts/others/SNPcalls.pl` and `SnpSimus.sh`), `annotateMGwSpecIs.pl`/`2.pl`, `parseBlastFunct.pl`, `secScripts/unused/`, `helpers/deprecated/`, `helpers/documentation_old/`, and the two unused bundled modules `Mods/ext/TreeIO.pm` and `Mods/IO/Uncompress/AnyUncompress.pm`.
- Config key removed in `bcb6e84`: `sepReadLength.pl`, `helpers/growthRate.pl`, `FMGgenes2tree.pl`, `NormMiTag.R`, `filterDeepCan.pl`, `reformatProGenomes.pl`, `prepare_strain_outgroup_refs.pl`, `helpers/merge-paired-reads.sh`, `unmerge-paired-reads.sh`, `WattetrsonTheta.hyphy`.
- Broken and unreferenced, failing the compile/config tests: `helpers/getFMGsFromMB_R.pl`, `helpers/motu_genome_integrate.pl`, `secScripts/GC/compareGeneCats.pl`.

Kept because live code uses them: `Mods/IO/PP.pm` (buildTree5), `Mods/ext/FAlite.pm` (assemblathon_stats) and `plan_consensus_regions.pl` (Mods/SNP.pm). Still present but unreferenced: `secScripts/MGS/extractTECx.pl`, which only drove removed scripts, `secScripts/others/distv8.pl`, `secScripts/SNP/createGenoCmp.pl` and `prepFullGenoPhylo.pl`.

## Tests

- Baseline before the audit: 77 files, 2,513 tests, 10 failures, all in `config_integrity.t`/`script_compile.t` and all caused by the retired-but-tracked scripts.
- After the fixes and removals: 78 files, 2,499 tests, all passing. The count is lower because `script_compile.t` no longer compiles the removed scripts.
- New `t/audit_2026_09_24.t` (17 checks) covers the stats parsers, read-group string, gzipped gene statistics, duplicate-log selection, per-log sdm histograms, ENA roles, geneCat batch bounds and the eggNOG filter. 14 of its checks fail against the pre-audit code.
- `t/workflow_control.t` gained two checks for quoted lightweight commands.
