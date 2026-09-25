# MATAFILER4 pipeline audit, second pass — 25 September 2026

## Scope and method

This pass followed the [24 September audit](../2026-09-24/report.md), whose fixes were still uncommitted. It had two parts:

1. **Review of the 24 September changes.** Every edited hunk was read against its callers (counter keys, directory conventions, module signatures, the SNP chunk counter, the KMA output prefix, the combine_DIA matrix hashes). No regression was found; the baseline suite (78 files, 2,499 tests) passes when run as `prove -I. t/`.
2. **A fresh audit of the whole active pipeline** in eleven slices: `MATAF4.pl` (five slices), the shared `Mods/`, the workflow-state modules, the gene catalogue and functional annotation, the MGS stage, `strain_within.pl`, `buildTree5.pl`, and the remaining secondary scripts. Reviewers were told to exclude everything already listed in the earlier reports, to verify every claim against the producer and consumer of each file, and to substantiate with a probe of the real subroutine where possible. The installer, the database configuration table and the helper scripts were checked separately.

Every fix below was verified in the source before it was made. The pipeline was **not** run end to end (no scheduler, no external tools). Verification is `perl -c` on every changed file, the full unit-test suite under WSL (perl 5.38), and the new `t/audit_2026_09_25.t`.

## Fixed

### `MATAF4.pl`

| Location | Defect | Fix |
|---|---|---|
| `prepareMap` | Secondary-reference gene prediction called `genePredictions` with the mapping output directory `GlbMap/<name>/` as its output directory. `genePredictions` starts its job with `rm -rf <outDir>`, so the job deleted the reference copy, the index being built alongside, and, on an existing project that later enabled coverage or SNP calling, every published per-sample secondary BAM and coverage. With `-mapModeTogether -1` it deleted the combined reference itself, so the prediction failed and blocked all secondary mapping. | predict in `GlbMap/<name>/genePred/`; the GFF is copied out of it |
| `scndMap2Genos` | `mapReadsToRef` received one sample name while its output directory, reference and mapping-directory arguments were per-reference comma lists. The mapper wrote `<sample>.iniAlignment.bam` (and a destination-less `mv` for every further reference), whereas the per-reference sort/coverage step read `<ref>_<sample>-0.iniAlignment.bam`; the completion check tested a file name the pipeline never creates. map2tar/map2DB could not finish and was resubmitted every pass. | pass one name per reference (the expression that had been left commented out); the read-group string uses the first name |
| `prepareMap` | With `-mapRefSNP` the consensus caller reads the reference as `GlbMap/<name>/<name>.fa`, but the only copy made kept the source file's basename (`Ecoli.fna`), and none was made without `-mapCov`. Every secondary SNP job died with "reference is missing". | link the reference as `<name>.fa` when SNP calling is requested |
| `sdmClean` | Singleton recovery globbed `filtered.*.singl.fq.gz`, which also matches the other libraries' (`filtered.lib1.singl…`) and the support scope's (`filtered.suppl.singl…`) final singleton files in the same directory. On a requeued cleaning job the first library appended and deleted them, so those reads ended up twice in the clean set. | glob only sdm's two mate files, `filtered.[12].singl.*` |
| `seedUnzip2tmp` | With `-usePorechop 1`, porechop ran for every single-end sample, but the staged path was only recorded for long-read samples: short-read singletons got a relative path, the staging stone was rejected and the job resubmitted every pass. | porechop only for long-read samples (same condition as the path update) |
| `mapReadsToRef` | A paired plus a singleton library in decoy mode were marked sorted although the parts are joined with `samtools cat`; with `-MapperRmDup 0` the index step failed. | sorted only for a single library |
| `reduceProgStats` | Decremented profiler failure counters unconditionally, so a sample that had already finished mOTUs/MetaPhlAn and was later finalised as empty drove the counter to -1, which the merge gates treat as "failures present": the cohort tables were never rebuilt. | never below zero |
| `postprocess` | The HTML report and the Kraken table merge ran `env:`-prefixed commands through `/bin/sh`; the activation prologue is bash syntax and fails on dash systems. Return codes were ignored. | run through bash, warn on failure |
| `riboSummary` | The skip check for the SSU/LSU cohort merges tested plain `SSU.miTag.<rank>.txt`, but `miTagTaxTable.pl` gzips every table, so both 80 GB merge jobs were resubmitted on every invocation and every loop pass. | gzip-aware check; tutorial updated |
| `submitGenomeBinner` | The binner was skipped whenever the assignment file existed, even without `Binning.stone`, which `runBinners.pl` writes last. A binner job that created the file and then failed (e.g. GenomeFace without `bins.tsv`) was accepted as "0 bins" on the next pass and the sample closed. | reuse only with the stone (present since the first MATAFILER4 release) |

### `strain_within.pl`

| Location | Defect | Fix |
|---|---|---|
| tree reset, `prepRun`, `evalFileStatus` | `-submit 0` is documented as a dry run, but with `-redo tree` it deleted every existing tree, MSA and `within/` result of the selected MGS, and with `-redo all` the whole output directory, before writing the planned commands. The `-redoEPAfilter` path already guarded its removals. | removals only with `-submit 1`; the dry run prints what it would remove |

### Workflow state (`-autoStatePlan`, `-inspectState`)

| Location | Defect | Fix |
|---|---|---|
| `Mods/WorkflowState.pm` | The mapping stage required the CRAM. With `-mapSaveCRAM 0` and a binner the finished-sample cleanup removes it on purpose, so the automatic preflight classified every completed sample as partial and, with repairs enabled (the default), deleted its mapping stone and coverage at every loop boundary; the sample remapped, the cleanup removed the CRAM again, and so on. | new option `mapping_cram_kept`, set by `MATAF4.pl` from the same rule as the cleanup; the CRAM is only required when kept |
| `Mods/WorkflowState.pm` | The hybrid preassembly package was inspected under the run scratch directory, but since `3d3a181` it lives under the sample output (`assemblies/preAssmblGrp_<group>/`). Every complete package was reported missing and re-planned. | inspect the durable location; three tests updated |

### Shared modules

| Location | Defect | Fix |
|---|---|---|
| `Mods/GenoMetaAss.pm` `splitFastas` | An existing split was reused when its first and last chunk existed. A split interrupted mid-way, or one made from an earlier catalogue with the same file name, was accepted as complete, so eggNOG-mapper and FOAM annotated a truncated or stale protein set. | a stone written last records the input's size/mtime and every chunk's size; anything else is re-split, and stale chunks are removed first |
| `Mods/FuncTools.pm` | The `.length` table that `parseBlastFunct2.pl` needs was only produced inside the DIAMOND index job. With `-functAligner foldseek`, or a DIAMOND index that already existed, no job produced it and every collection job died. | the length step is submitted on its own when the table is missing |

### Gene catalogue, MGS and phylogeny

| Location | Defect | Fix |
|---|---|---|
| `phylo_MGS_between.pl` | The launcher rewrote `all.faa`/`all.cats` on every invocation. `buildTree5 -continue` fingerprints its inputs by size and modification time, so every relaunch of the between-MGS tree (e.g. after a wall-time kill) discarded all finished per-locus alignments. Its own reuse check compared the buildTree5 policy schema with the literal `13`, which buildTree5 stopped writing several releases ago. | write to temporaries and publish only when the content differs; compare the policy values, not the schema number |
| `annotateMGwSpecIs3.pl` | GTDB taxonomy padding only replaced the first and the trailing empty rank, so an MGS novel at family level carried empty rank names into `specI.tax`/`.genus`/`.species`. | every empty rank becomes `?` |
| `extrAllE100GC.pl` | Died for any sample without `assemblies/metag/assembly.txt`, although geneCat's collation was told to tolerate it (`-requireAllAssemblies 0`), blocking MG_LCA and MGS. | skip such samples like the collation does |
| `parseBlastFunct2.pl` | `-percID` was integer-only while geneCat's `-FuncMinPerID` is a float: a value such as `27.5` killed every collection job after the alignments had run. | float |
| `decluterGC.pl` | Hard-coded `sse4` Slurm constraint. | the configured `avx2_constraint`, like the other jobs |

## Found, not fixed

- **`-decluterMatrix 1` cannot be resumed.** `decluterGC.pl` replaces `Matrix.mat.gz` and deletes `compl.incompl.<id>.fna.clstr*`, both of which are outputs recorded in the matrix and publish checkpoints written moments earlier. A rerun re-enters clustering or dies "GC Files not found". Off by default; fixing it means recording the decluttered products in the checkpoints.
- **buildTree5 alignment policy ignores the sample-QC content.** Only the count of excluded samples is in the policy; a changed set with the same count reuses per-locus alignments built from the old sample set. Adding the fingerprint invalidates every existing alignment checkpoint once.
- **strain_within split Phase I resume can dead-end in the worker-repair queue.** `mergeRecoveryLogs` deletes the per-worker recovery ledgers, but a later non-lean resume (`-redoEPAfilter`, or a missing sample-stats summary) re-validates exactly those ledgers, declares every worker invalid and exits with "Phase I requires worker repair" although Phase I was fully merged. Narrow trigger; the fix is to treat a merged generation as merged.
- **`-normSNPindels` on a project finished with `-SNPsaveVCF 0`** submits one failing normalisation job per sample before the full rerun (opt-in path).
- **Perl `-perlClusterMAGs` report schema.** The naming is fixed (see below), but the report still has six metadata columns where `markersPerMGS.pl` expects 22 (MGS audit of 11 September, item 5).
- **Dead code noted by reviewers** (not changed): the SPAdes OOM-retry regexes only match LSF logs; `spadesHosts` is disabled; `getGenesInSmpl.pl` calls `getProgPaths` without importing it and has no live caller; `post_alignment_locus_qc.pl` is not on any production path; `FuncTools.pm`'s "random" `diaFunc/10/` directory is never used.

## Decisions applied (follow-up, 25 September)

After reviewing the open list, the maintainer decided the following; each was implemented and tested the same day.

| Item | Decision | Change |
|---|---|---|
| Broken `refs/codex/…` references | delete | The five checkpoint references of another tool, whose paths exceeded the Windows limit, were removed. Automatic repack and `git fetch` work again. |
| CheckM2 database path | follow the installer | `checkm2DB` is `[DBDir]/CM2/CheckM2_database/uniref100.KO.1.dmnd`, where `installer.sh` downloads it. Sites with their own `checkm2DB` are unaffected; others rebuild MGS canopy checkpoints once. |
| Shared assembly from a member subset | reopen every member | A group rewrite now removes the completion sentinels of all members. In addition, a group assembly never starts while members of the group took the "already complete" fast path in the same pass (they registered no reads): their sentinels are removed and the assembly starts on the next pass (`deferGroupAssemblyForClosedMembers`, both assembly start points). |
| Remap keeps old bins | leave as is | No change. |
| Paired unaligned BAM input | implement for sdm | The layout of every SAM/BAM/CRAM input is detected once per file (first primary record FLAG `0x1`, header `SO:`). sdm reads a paired file as true pairs (`-paired 2`); a file without `SO:queryname` is name-sorted into the cleaning directory first and the copy removed afterwards. The cleaned library is an ordinary R1/R2/singleton set. Raw-read paths (single-pass raw mapping and the raw FASTQ cache) write every mate to one stream instead of dropping them, so those mappers treat mates as unpaired reads; old caches are rebuilt once. Upload preparation stops for paired alignment input. |
| MGS abundance checkpoints | fix | Both abundance stages are checked with the resume validator and written without the empty-sample list. |
| extrAllE100GC.pl marker stone | fix | The script waits for its matrix-subset jobs and verifies their outputs before exiting, so geneCat's marker stone follows real products. The misspelled skip test now checks the complete output set. |
| kmerPerGene.pl without k-mers | fix | The script checks all assemblies first and skips cleanly (no table) when MATAF4 did not compute per-gene k-mers; geneCat compresses the table only when written. |
| hostile and ProstT5 paths | fix | `[DBDir]/hostile/` and `[DBDir]/PtostT5_W`, like the other installer-managed databases. |
| `-perlClusterMAGs` naming | fix | MAG ids are `sample__bin` and `Cano__x`, as the Binning readers and the native clusterMAGs expect. |
| `get_ranks.py` | Perl port | `secScripts/phylo/taxid2ranks.pl` reads the lineages from the Kraken2 database that assigned the IDs (`taxo.k2d`, layout checked against the Kraken2 source), or from an NCBI taxdump in `<db>/taxonomy/`. It accepts NCBI's 2025 `domain` rank as well as `superkingdom`, and reports unknown IDs with empty ranks instead of dropping them. geneCat calls it once instead of in 500-ID batches. No ete3 is needed. |
| `get_gtdb.py` | Perl port | `helpers/install/get_gtdb.pl` reproduces every subcommand (`download`, `extract`, `configure`, `all`) and version class, streams the large archives through `tar` instead of joining them on disk, and reads `MF4DIR` (then `MGTKDIR`, then its own location). It fixes several defects of the Python (crash without `MGTKDIR`, `--tk skip` crash, `--tk full` overwriting its archive, the test stub mismatch, repeated mash-line comments, a GTDB-Tk version table that would downgrade the installed environment) and writes the tables with LF line endings, since the CR of the Python output ended up in accession IDs. `t/get_gtdb.t` (149 checks) runs it end to end offline on the bundled stubs. The internal layout of the real multi-GB archives and releases older than R214 could not be checked offline. |

## Tests

- `t/decisions_2026_09_25.t`, `t/taxid2ranks.t` and new checks in `t/alignment_input.t` and `t/main_submission_graph.t` cover the decisions above. The real bundled sdm cleans a paired SAM fixture into complete `filtered.1`/`filtered.2` files.
- `t/audit_2026_09_25.t` (23 checks): the splitter's stone, reuse and re-split rules; the workflow-state CRAM rule and package location; the counter reset; GTDB rank padding; and source guards for the command-construction fixes that cannot run here.
- `t/workflow_plan.t` and `t/workflow_runner.t` now place the hybrid package under the sample output.
- Full suite after this pass: 79 files, 2,522 tests; after the follow-up decisions: 82 files, 2,722 tests. All passing (`prove -I. t/` under WSL, perl 5.38).
