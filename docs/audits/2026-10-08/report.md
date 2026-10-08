# RiboFind audit — 8 October 2026

## Scope and method

This audit covered the whole `-profileRibosome` path (SSU/LSU miTag profiling):

1. option handling and redo flags (`normalise_ribosome_request`, `prepare_ribosome_rerun`);
2. completion evidence (`ribosome_completion_evidence`, the sample sentinel), `checkRawProgsFin` and `RiboMeta`;
3. `detectRibo`: preparing the LCA reference copy and submitting jobs;
4. read extraction, `secScripts/miTag/catchLSUSSU.pl` (SortMeRNA);
5. assignment, `secScripts/miTag/lotus_LCA_blast3.pl` (FLASH, LAMBDA/VSEARCH/SortMeRNA, LCA);
6. the cohort merge (`riboSummary`, `secScripts/miTag/miTagTaxTable.pl`).

Issues already listed in the 2026-09-11, 2026-09-24 and 2026-09-25 reports were left out. A first pass reported its findings, and decisions on four of them followed the same day (see "Decisions" below).

Each defect was reproduced before it was fixed, by running the real scripts on synthetic inputs under WSL Ubuntu (perl 5.38). Most probes used stand-in programs. The SortMeRNA check used the real SortMeRNA 7.0.0 and 4.3.6 (installed into a scratch prefix) and the bundled `bin/LCA`. LAMBDA, VSEARCH and FLASH were not available. The pipeline was **not** run on a cluster.

Verification:

- `perl -c` on every changed file.
- The full unit-test suite under WSL. Before the changes: 86 files, 2,751 tests, all passing. After: 87 files, 2,810 tests, all passing.
- New `t/audit_2026_10_08.t` (59 tests). Against the old scripts and `MATAF4.pl`, every check it reaches fails (10), and it then stops at the sample-list merge, which the old merger does not support.

## Fixed

| Location | Defect | Fix |
|---|---|---|
| `catchLSUSSU.pl` (completion check) | **Endless resubmission.** MATAF4 requires all three published read files per marker (`reads_<tag>.r1/.r2/.fq.gz`). `catchLSUSSU` only checked the roles of the current input. Before 0.6, a pair-only sample got no singleton file (or a zero-byte placeholder). For such a sample `catchLSUSSU` exited "complete" without writing anything, and MATAF4 found the profile incomplete again: a new RiboFinder job on every pass, and the SSU/LSU tables were never merged. | Same contract as MATAF4 (stone plus all three files). A role the input does not have gets the empty gzip container that a fresh run writes, without a new SortMeRNA search. A role the input does have is never filled in this way: that marker is rerun. `catchLSUSSU` 0.7. |
| `lotus_LCA_blast3.pl` `materializeReadInput`; `catchLSUSSU.pl` | **Old reads assigned after re-extraction.** The first `lotus_LCA` versions ran `gunzip` on the published reads inside `ribos/`, and a plain `reads_<tag>*.fq` was preferred over the `.gz`. Since 4.39 such samples are extracted again (the `.gz` files were missing), but the old plain reads were the ones assigned. | The published `.gz` is preferred. When a marker is rerun, `catchLSUSSU` also removes the old plain read files. |
| `lotus_LCA_blast3.pl` | **Missing reads checkpointed as complete.** Without the extraction checkpoint (`<tag>_pull.sto`), missing reads looked like a sample without ribosomal reads: empty hierarchies, `<tag>_ass.sto` and `Assigned.sto` were written (manual runs, `afterany` dependencies). | A marker is assigned only after its extraction checkpoint exists. Otherwise the job fails and names the checkpoint. |
| `lotus_LCA_blast3.pl` `-simMode 3` | **SortMeRNA search mode ignored `-cover`.** `LCA` reads column 11 as the query length (checked with `bin/LCA`: query length 1000 for 150 aligned bases assigns 0 of 100 reads). SortMeRNA's `--blast 1` table has the e-value there, so every hit passed the coverage filter. MATAF4 itself uses LAMBDA (`-simMode 2`). | SortMeRNA's table is rewritten to the 11-column layout LAMBDA and VSEARCH produce, with each read's length from the query FASTA. |
| `detectRibo` (LCA reference copy) | **PR2 broke every LCA job.** `config_DBs.txt` offers `PR2dbFA`/`PR2tax` (commented out). `lotus_LCA_blast3` uses them for SSU, but `detectRibo` copied only LSU/SSU into `DB/LCADB/`, so every LCA job died with "Missing SSU reference database". | A configured PR2 reference is copied too. |
| `detectRibo` (LCA reference copy) | **Truncated copies kept for ever.** A copy counted as present if the FASTA and `.lba.gz` existed, so the truncated files of an interrupted copy job made every later LCA job fail. The taxonomy file was not checked. An unset taxonomy key would have turned `cp $tax*` into every file in the job directory. | `lca_reference_copy_current` compares FASTA, index and taxonomy with their sources by size. A missing reference file, or a missing or unset taxonomy, stops MATAF4 at startup instead of failing every LCA job. |
| `detectRibo` (LCA reference copy) | **No database copy for a second output folder.** "Already copied" was one global flag, but the copy goes to `<OutPath>/DB/LCADB/`. With several `#OutPath` folders, the samples after the first folder got no copy job and their LCA jobs failed. | The copy is tracked per database directory. |
| option handling | **`-riobsomalAssembly 1` failed every sample, every pass.** `catchLSUSSU` refuses assembly, and nothing writes the `Ass/allAss.sto` that completion then requires. The flag was documented as stable. | A nonzero value stops MATAF4 at startup. The flag reference marks it deprecated/legacy. |
| `RiboMeta` | The hierarchy was compressed with `gzip` without `-f`: after an interrupted gzip, the partial `.gz` was kept and published. | `compress_ribosome_hierarchies` uses `gzip -f`. A hierarchy lost after the completion check now repeats only the assignment (the extracted reads are intact); before, an SSU loss deleted all of `ribos/` and an LSU loss was ignored. |
| `riboSummary` | `my $mrgCmd = … unless (…)` is undefined behaviour in Perl. | Replaced (the merge was rewritten, see below). |
| `detectRibo` | RiboFinder's scratch size was computed from `$map{<SmplID>}`, which is keyed by map key: 0 MB, so always the floor, whenever the two differ. | Uses the map key. |
| `docs/profiling_tutorial.md` | It said the SortMeRNA keys (`SSUdbFAsrt`, `LSUdbFAsrt`, `SSUidx`, `LSUidx`) are "not read by the current `detectRibo()` workflow". `catchLSUSSU` requires both references. | Corrected, with how to build an index directory and the optional PR2 reference. |

## Decisions (follow-up, same day)

### LCA reruns do not clean reads again

Decision: rerunning the LCA step must not trigger read cleaning.

Before, `calcRiboAssign` was part of `$stagedReadsAnalysisFlag`. By default (`-reduceScratchUse 1`) a finished sample's scratch is removed, so `-reRibosomeLCA 1`, or any failed LCA job, staged the raw reads and ran sdm again on every affected sample. The LCA job waited for that, although it only reads `<sample>/ribos/reads_*.fq.gz`.

Now an assignment that is pending without extraction is handled like Protal. The LCA job is submitted right after the empty-sample check, with no staging or cleaning dependency, and the sample stays open until it finishes. `detectRibo` reads the cleaned-read set only when it submits an extraction. The LCA job takes the read layout from what `catchLSUSSU` published (`-pairedRds 2`), so it no longer needs the cleaned reads either. It waits for the LCA database copy directly. Before, it got that dependency only through the extraction job.

### SSU/LSU tables: exactly the samples of the map

Decision: merge only the samples of the current map, across several output folders or map files, and always the whole map, not the `-from`/`-to` part.

Before, the tables were merged from per-sample links in `<baseOut>/pseudoGC/Phylo/RiboFind/{SSU,LSU}/`, where `<baseOut>` is the folder of the sample being processed. With several `#OutPath` folders or map files, links and tables were split across folders. Samples dropped from the map stayed in the tables. The merge was triggered when the count of complete samples in the current pass grew. Closed samples never re-created a lost link (fixed in the first pass with a link-restoring step, now replaced).

Now `riboSummary` checks every sample of every `-map` file, independent of `-from`/`-to`, from its own output folder. It leaves out `-ignoreSmpl` samples, samples marked `SMPL.empty`, and samples whose completion sentinel records a `skipped_*` outcome. All other samples must have a complete profile and assignment (`ribosome_completion_evidence`). Otherwise nothing is merged, and the first missing samples are named.

The merge gets an explicit sample list (`<marker>.miTag.samples.tsv`: column name, hierarchy path in the sample's own folder). `miTagTaxTable.pl` accepts such a list, fails on a listed file that does not exist, and still accepts a directory, where it now names dangling links. The tables go to `pseudoGC/Phylo/RiboFind/` beside the first map sample's output, like the Protal cohort output.

`<marker>.miTag.sto` records a signature over the merged columns, paths, sizes and modification times. A new merge runs when it changes: samples added or removed, another map, a re-assigned sample. The per-sample link folders are no longer written, and the RiboFind progress counters are gone.

### FLASH merging in the LCA step: unchanged

Decision: keep it, because SortMeRNA handles paired reads itself and the read-merge path may be retired as a whole. For reference: with the bundled sdm 3.61 beta, `-merge_pairs_filter 1` reports merged pairs but writes an empty `.merg.fq` and keeps every pair in the pair files.

### SortMeRNA: version and options

`helpers/install/MF4.yml` installs `conda-forge::sortmerna=7.0.0`, the latest release (8 July 2026; 6.x replaced CMPH with BBHash and SSW with Parasail). Checked against the 7.0.0 source and by running it:

- Every option `catchLSUSSU` passes still exists with the assumed meaning: `--ref`, `--reads` (twice for pairs), `--idx-dir`, `--index 0` (use an existing index, never build one), `--workdir`, `--aligned`, `--fastx`, `--paired_in`, `--out2`, `--no-best` with `--num_alignments 1` (stop at the first alignment that passes `-e`), `--zip-out 1`, `--threads`, `-e`.
- Output names are `<prefix>_fwd.fq.gz`/`_rev.fq.gz` (paired) and `<prefix>.fq.gz` (single), which is what `catchLSUSSU` expects. On synthetic data (200 rRNA and 200 random pairs, 100 + 100 singletons), `catchLSUSSU` with 7.0.0 published exactly the 200 rRNA pairs and the 100 rRNA singletons, both with and without a configured index directory.
- Index directories built with 4.x still work: 7.0.0 names the four index files the same way and only checks that they exist. An index built by 4.3.6 had the same file sizes (internal IDs differ), and 7.0.0 returned exactly the reads it returns with its own index. The configured `SSUidx`/`LSUidx` directories need no rebuild. (The bioconda 4.3.7 build crashes with "Illegal instruction" on the test machine, so 4.3.6 was used.)
- The reference sets are unchanged: 7.0.0 ships the same `smr_v4.3_*` databases (same sizes) as 4.3.7.
- 7.0.0 can resume an interrupted run from its work directory. `catchLSUSSU` gives every SortMeRNA call a new work directory and removes it afterwards, so no stale state is resumed.

### What changes for existing projects

- Pair-only samples whose profile was missing the singleton file now close without re-extraction.
- Samples with old plain read files are assigned from their new extraction.
- The first pass checks `DB/LCADB/` by size. Complete copies are kept.
- The SSU/LSU tables are merged once more under the new rule. They then contain exactly the map's samples, written beside the first map sample's output. The old `RiboFind/SSU/`, `RiboFind/LSU/` link folders and `*.cnt.stone` files are no longer used. They are not deleted, except a sample's own links by `-reRibosomeLCA`/`-reProfileRibosome`.
- `-reRibosomeLCA 1` no longer stages or cleans reads.
- RiboFind runs stop at startup if a configured LCA reference or its taxonomy is missing, or if `-riobsomalAssembly` is nonzero.

## Found, not fixed

- **`LCA -readInput` is not used.** `bin/LCA` documents `-readInput` ("Are the inputs miTags? Default: off (assumes OTUs)"). `lotus_LCA_blast3.pl` assigns miTag reads but does not pass it. What the switch changes was not checked.
- **Unmerged read pairs are probably not assigned.** `mergePairs` joins R1 and reverse-complemented R2 with no spacer, and LAMBDA aligns locally. A pair with an insert longer than both reads together aligns over about one mate, which fails `-cover 0.85` (column 11 is the joined length). So effectively only FLASH-merged pairs and singletons are profiled. Kept as decided above; not confirmed with LAMBDA.
- **`-riboMaxRds` keeps the first N candidates,** not a random subsample: merged pairs, then unmerged pairs, then singletons, earlier libraries first. `-riboMaxRds` and `-saveRiboRds` are not part of the completion signature, so changing them leaves closed samples as they are.
- **PR2 normalisation in `miTagTaxTable.pl` handles Opisthokonta only.** Other PR2 eukaryotes keep their supergroup in the phylum column. Relevant now that a configured PR2 reference is used. The KSGP taxonomy was not checked.
- **MetaPhlAn and mOTU merges** still take every profile in their central folder, as the RiboFind merge did before.
- **`RedoRiboThatFailed`** (in `RiboMeta`) is set by no command-line option.
- **`-thoroughCheckRiboFinish 1`** still deletes valid empty hierarchies and loops (opt-in; listed on 2026-09-24).

## Tests

- `t/audit_2026_10_08.t` covers:
  - the extraction contract: a legacy pair-only profile, a missing input role not filled in, old plain files removed;
  - assignment: published reads preferred, refusal without an extraction checkpoint, the SortMeRNA search table;
  - the reference-copy check and hierarchy compression;
  - the merge from a sample list, and dangling links in directory mode;
  - the real `riboSummary`, run over a map spanning two output folders with empty, skipped, ignored and unprofiled samples: the merge input, the tables, no merge when nothing changed, a new merge after a map change or a re-assignment, no partial merge;
  - `-riobsomalAssembly`;
  - source checks of the main-loop wiring: an LCA-only rerun without staging, per-folder database copies.
- On Windows-native Perl, the three directory-merge checks that need POSIX symbolic links are skipped.
- `t/mitag_reference_validation.t`: the malformed-FASTQ fixture now has the extraction checkpoint a real extraction always writes.
- `t/audit_2026_09_25.t`: its source check of the merge's gzip-aware skip test follows the new code.
