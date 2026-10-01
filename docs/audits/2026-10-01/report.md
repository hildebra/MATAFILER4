# Functional-annotation audit — 1 October 2026

## Scope and method

This audit was triggered by KEGG DIAMOND chunk jobs (`DKGM.<i>.sh`) running out of node-local scratch. It covered the whole functional-assignment path, in five slices:

1. the geneCat orchestration (`FuncAssign`, `FuncEMAP`, `FOAM`/`ABR` modes, stage submission and checkpoints);
2. `Mods/FuncTools.pm` with its callers, `getSpecificDBpaths` and the shared catalogue split;
3. the hit parser `parseBlastFunct2.pl`, with `prepVFDB.pl` and `geneLengthFasta.pl`;
4. eggNOG-mapper, the KEGG module step and the HMM path;
5. read-based functional profiling in `MATAF4.pl` (`-profileFunct`) and the cohort merge `combine_DIA.pl`.

Reviewers excluded everything listed in the earlier reports. They checked every claim against both the producer and the consumer of each file, and confirmed most findings by running the real subroutine or script on synthetic inputs, with the scheduler stubbed. Every fix below was verified in the source before it was made.

The pipeline was **not** run end to end on a cluster. Verification:

- `perl -c` on every changed file;
- the full unit-test suite under WSL (perl 5.38): 85 files, all passing;
- the new `t/audit_2026_10_01.t` (13 subtests). Run against the unmodified code, it fails on the first defects it reaches.

Three suite failures that predated this audit were fixed separately (see the end of this report).

## Scratch space (requested change)

DIAMOND writes intermediate seed hits and per-block alignment results to `-t` (node-local scratch). For KEGG (a large, highly redundant reference), a 500 MB catalogue chunk ran out of the 320 GB scratch request.

| Jobs | Before | Now |
|---|---|---|
| FuncAssign chunk jobs for `KGM`, `KGE`, `KGB`, `NOG` | 320 G (KGM/NOG), 40 G (KGE/KGB) | 500 G |
| FuncAssign chunk jobs for every other database | 40 G | 250 G |
| eggNOG-mapper chunk jobs (`FuncEMAP`) | 77 G | 500 G |

If 500 G is still not enough, the other levers are:

- a smaller `-fastaSplit` (temp use scales with the query letters per job);
- DIAMOND `--hit-membuf` (seed hits in RAM);
- a KEGG reference restricted to genes with a KO (the parser only ever assigns KO-bearing subjects).

## Fixed

### Gene catalogue: `secScripts/geneCat.pl`

| Location | Defect | Fix |
|---|---|---|
| `geneCatFunc`, stage submission | Changing `-functAligner` or any `-FuncMin*` cutoff invalidated `10.func.stone`, but FuncAssign then found every per-database marker present and submitted nothing. Its done job re-stamped the stone with the new parameters over results computed with the old ones. `-redoFunc`, the documented remedy, was not forwarded from the main run. | Per-database record `Anno/Func/.<db>.params` of the aligner and the effective cutoffs (after `%funcDBcutoffs`). A changed aligner, or a less strict e-value than the alignments were made with, realigns that database. Any other change (including a stricter e-value) deletes only its per-gene assignments and matrix marker, so the existing alignments are re-interpreted; the alignment e-value is kept and also used for chunks still to be aligned. KGM module tables are removed so they are rebuilt. Products without a record (older runs) are kept. `-redoFunc` is forwarded and forces the stage. |
| FuncAssign/FuncEMAP modes, main flow | Re-running geneCat while a stage's jobs were queued submitted the whole graph again: duplicate chunk jobs, and the first run's final job deleted the catalog split under the second run's jobs. | In-progress marker `Anno/Func/.<stage>.inflight`: created exclusively before submitting, holds the final job, removed by it. While that job is queued or running neither the main flow nor the mode submits the stage. Stale markers (final job finished, Slurm `DependencyNeverSatisfied`, or a dead submitter) are removed. |
| stage submission | A legacy empty `10.func.stone` / `10.emap.stone` (`touch` era) counted as valid for any database list, aligner or cutoffs (`Checkpoint.pm` accepts empty stones). | Both stages require a non-empty stone. A finished legacy catalogue reruns each stage once; nothing is realigned. |
| `geneCatFunc` | Marker `.<db>.matrix.done` present but `DIAass_<db>*` deleted (documented as intermediates): FuncAssign submitted nothing and the stone job died "Cannot checkpoint missing output" on every pass. | The database is recomputed when its per-gene file is missing. |
| FuncAssign done job | The stone recorded only the per-gene files. With SGE `-hold_jid` or `afterAny`, a failed matrix job was followed by a valid stone. | The matrix markers are recorded as stone outputs. |
| main flow, `-submitLocal 0` | `CDHITexe.sh` ends with `rm -rf $tmpDir`. The FuncAssign/FuncEMAP stage commands run inside it only to submit jobs, and their split and chunk outputs live under `$tmpDir`. With `-doMags 0` every queued chunk job lost its input, and the stage could never finish. | The final cleanup keeps `funcSplit_*`, `GCanno_*` and `eggNOGmapper_*`; their own final jobs remove them. |
| stage submission | `-submSystem`, `-tmp` and `-glbTmp` were not passed to the FuncAssign/FuncEMAP subprocesses. With `-submSystem bash` on a host without a scheduler, FuncAssign died "No queueing system found" on every pass. | Forwarded. The resulting paths equal the former defaults. |
| `geneCatFunc_emapper` | eggNOG-mapper writes `<-o>.emapper.annotations` in place, and `--no_file_comments` drops its only end marker. A chunk killed mid-run left a file the resume test accepted, so genes went missing silently or CombineEMAP failed forever. | Runs with `-o <chunk>.part` and renames on success. |
| CombineEMAP | `grep … \|\| true` also swallowed exit 2 (missing chunk file), so a failed chunk was dropped from the merge. | Only grep exit 1 (no rows) is tolerated. |
| EM.* matrices | `-extHiera -hieraSrtDown` made rtk split `KOs;modules` on `,` before the hierarchy. Level 1 of `EM.KGM`/`EM.KGP` held `KO;module` path names, so modules and pathways were never totalled across KOs. | Flags removed. Single-level tables and level 0 are unchanged. |
| option check | `-fastaSplit` accepted `K` and lower-case suffixes that `splitFastas` treats as chunk counts (`500m` = 500 jobs). | Only a count, or a size ending in `M` or `G`. |
| FuncAssign pre-check | `mp3` passed the pre-check (no DB files) and then died in FuncTools after earlier databases were already submitted, on every pass. | The pre-check requires the `mp3` program path. FuncTools no longer requires DB files for mp3. |
| `FOAMassign` (modes FOAM/ABR) | The best-hit helper needed python2, which no installer environment provides. It also kept the wrong hit for every gene after the first (`maxscore=0` on a new gene) and crashed on a chunk without hits. FOAM and ABR shared one split, and each chunk job deleted its input chunk. | Perl port `secScripts/functions/hmmBestHit.pl`: highest score per gene, no sorted input needed. The Python helper was removed. The split is per mode, and chunks are removed by the collect job. |

### `Mods/FuncTools.pm`

| Location | Defect | Fix |
|---|---|---|
| chunk loop | Outputs are named by chunk index only. After a re-split with another `-fastaSplit`, leftovers were concatenated with new chunks: genes were duplicated, or missing with fewer chunks. | A chunk output older than its catalogue chunk is discarded (same rule as eggNOG-mapper). |
| chunk commands | Every duplicate submission of a chunk wrote the same `.tmp.gz`. | Job-unique temporary name (`.tmp.$$.gz`). |
| `assignFuncPerGene` | Per-gene file present and `DIAass` absent: every chunk was realigned with no collection job, and the done job then deleted the split under them. | No realignment when the per-gene file exists. |
| collection job | Did not wait for the DB-length / `VF.tab` jobs when no chunk job was submitted, so the parser raced `-LF`. | Depends on them. |
| DB build | `diamond makedb` wrote the final `.dmnd`, which is only tested with `-e`, so a killed build was accepted forever. All databases shared one `DiamondDBprep.sh`/`DBlength.sh` script name, so their logs clobbered each other. | Built under a temporary name and moved; per-database script names. |
| `VF.tab` | Built only when missing. An updated VFDB, or a table in the pre-0.61 layout, killed the collection after all alignments had run. | `vfTabStale`: rebuild when missing, older than a VFDB FASTA, or with fewer than 7 columns. |
| `calc_modules` | The module definitions (`Module_path_DB`) are not installed by the installer. Without them `rtk module` failed. Since c43dee1 this also blocked `10.emap.stone`, and with `-strains 1` MGS.pl then waited 24 h and died. | A module set with a missing `.list`/`.descr`/`_hiera.txt` is skipped with a warning. |
| node scratch | `getProgPaths("nodeTmpDir")` was required here but optional everywhere else. | Optional, falling back to the shared tmp dir. |

### Hit parser and database helpers

| Location | Defect | Fix |
|---|---|---|
| `parseBlastFunct2.pl` (TCDB, NOG, PTV per-gene files) | Free-text descriptions are written as hierarchy levels. geneCat turns tabs into `;` and rtk splits on `;` `,` `\|`, so a TCDB class such as "Porters (uniporters, symporters, antiporters)" became three features carrying the full count (TCDB L1/L3, NOG L2). | `,` `;` `\|` in descriptions become `_`. |
| `parseBlastFunct2.pl` (CAZy) | Multi-family subjects were written `GH13;CBM48`, so the second family became the substrate level and the substrates a new level. | Families joined with `,` (AND); duplicate substrates removed. |
| `parseBlastFunct2.pl` | A truncated `.gz` (trailer missing, cut on a line boundary) was accepted and the stone written. | `close` is checked. |
| `prepVFDB.pl` | A `[…]` in the protein description (e.g. `[2Fe-2S]`) started the VF bracket, so one factor became two L1 features. | VF and category names may not contain brackets. |
| `geneLengthFasta.pl` | Wrote the `.length` table in place; callers only test existence. | Temporary file and rename. |
| `getSpecificDBpaths` | Callers build `$DBpath$refDB`, so a configured path without a trailing `/` passed the checks and then failed. | Slash appended when missing. |

### Read-based profiling: `MATAF4.pl` and `combine_DIA.pl`

| Location | Defect | Fix |
|---|---|---|
| `runDiamond` (ABR) | Since 17608f1 the ABR filter gets the run-local DB copy as its database dir. That copy has no ARDB tables, so every ABR parse died and the sample never closed. | The configured ABR database dir is passed. |
| `prepDiamondDB` (VFA/VFB/VDB) | Only copied `VF.tab`, which only geneCat built. On a fresh install the DB-prep job failed, and the stage could never complete. | Builds `VF.tab` with `prepVFDB.pl` when stale or missing. |
| `runDiamond` | With `-rmRawDiamondHits 1`, a parsed database without raw hits was realigned and reparsed whenever another database of the sample was pending. | Not realigned when its parse marker exists. |
| `prepareDiamondRerun` | With exactly `maxReqDiaDB` (6) databases requested, `-reParseFunct` removed every `CNT*` table while keeping other databases' markers (permanently missing tables), and `-reProfileFunct` removed all of `diamond/`. | Always per database. |
| `prepareDiamondRerun` | The mode-4 cleanup omitted `-percID`, so it targeted `CNT_<eval>_20` instead of `CNT_<eval>_<DiaPercID>`. | `-percID` passed. |
| `runDiamond` | All databases' jobs of a sample shared one temp dir and each removed it at the end. | Per-database temp dir. |
| `IsDiaRunFinished` | On a parse-only pass (e.g. `-reParseFunct 1`) the merged library is not registered, so old merged-read hits without `MF4:read_count` were reparsed instead of regenerated. Merged pairs were counted once (contrary to the 11 September note on cached searches). | Such hits schedule the search stage when read merging is on. |
| `runDiamond` | Read hits were sorted with the locale's collation; the parser needs each read's lines adjacent. | `LC_ALL=C sort` (hardening, not observed). |
| `combine_DIA.pl` | Took the table layout from the first map sample only. If that sample had no `diamond/` (empty sample, `-from N`), the merge died inside the controller on every pass, and the later postprocess steps were skipped. | First sample with tables is used; a database without any is skipped. |
| `combine_DIA.pl` | TCDB and VFA/VFB/VDB category matrices were always header-only: it read `gene.cnts`, but the parser writes `.CATcnts` for them. | Added to the `.CATcnts` branch. |

### Documentation

`docs/outputs.md`:

- real module paths (`Anno/Func/modules/<set>/`);
- the per-database parameter record and re-run rules;
- the corrected `EM.KGM`/`EM.KGP` level 1;
- `VF.tab` rebuild.

`docs/flag_reference.md`: `-fastaSplit`, `-redoFunc`, `-rmRawDiamondHits`.

`docs/profiling_tutorial.md`: `VF.tab`.

## Found, not fixed

- **Pre-0.61 catalogues** have no parameter record, so their results are kept as computed. That is before the query-coverage alternative existed. To apply the current rules, run once with `-redoFunc 1`. This realigns everything; a re-interpretation-only switch would be cheaper if wanted.
- **Duplicate DB builders.** `VDB` and `VFB` map to the same FASTA, and MATAF4's `prepDiamondDB` builds the same `.dmnd`/`.length` in place. Concurrent first use can race; FuncTools' own build is now atomic.
- **`combine_DIA.pl` self-heal** removes `dia.<DB>.blast.gz.stone`, a name no longer used, so it never triggers. It was left inert on purpose: activating it would loop reparsing whenever samples differ in their `CNT_*` directories.
- **`lambdaBl`** (MG_LCA marker-gene taxonomy, outside this path) renames the shared reference FASTA while building its index. An interrupted build leaves it renamed and every retry fails.
- **`readTabbed3`** accepts rows one column short of the requested column (undef value). It was not changed, because other callers may rely on it; the `VF.tab` layout is checked instead.
- **Dead or dormant code:** `phyloTools::prepNOGSETgenomes` (calls `assignFuncPerGene` with 5 arguments, no caller); `passBlast` (no caller); foldseek in local mode (`buildFSdb` without options would die); the parser's legacy KGB+KGE concatenation (unreachable). DB index and length jobs are still submitted before the early return when nothing will be aligned.
- **Query-coverage-only cutoffs on 12-column input** (no `qlen`) reject every hit. No current caller combines the two.

## Notes for running projects

- Finished chunk outputs, `DIAass_*` and per-gene files are reused. On the next FuncAssign run each finished database gets a `.<db>.params` record for the current settings; nothing is recomputed.
- Catalogues with a legacy empty `10.emap.stone` rerun FuncEMAP once from the existing merged annotations: category split, matrices, modules. The `EM.KGM`/`EM.KGP` level-1 tables then change meaning, as described above.
- FOAM/ABR results change: the best hit per gene is now the highest-scoring one.

## Pre-existing test failures (fixed the same day)

Three test files failed on the unmodified `HEAD` as well:

| Test | Cause | Fix |
|---|---|---|
| `t/script_compile.t` 44 | `helpers/combineMarkers.pl` has had a syntax error since the first commit (`".@samples" samples"`, missing `.`). Machines without `List::MoreUtils` skipped the test, so it went unnoticed. | `scalar(@samples)` with the operator. |
| `t/ena_sra_download.t` | PR #12 (ENA run selection) requests `library_strategy` from ENA and requires the column; the test's ENA report fixture lacked it. With that fixed, the SRA case failed too: PR #11 passes the resolved `.sra` object to `vdb-validate` and `fasterq-dump` (verified on real downloads), but the test's fake tools still expected the prefetch directory. | Fixture has `library_strategy`; the fake tools expect the resolved SRA object. |
| `t/audit_2026_09_24.t` 13–14 | aff3b10 moved the sample-batch arithmetic into `_prep_batch_range`; the test still searched for the old inline expressions. | The test checks `_prep_batch_range` with the same cases. |
