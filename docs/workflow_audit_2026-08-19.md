# Workflow audit: MATAF4 to within-strain analysis

Date: 2026-08-19

## Scope and standard used

This audit follows the active path through:

1. `MATAF4.pl`
2. `secScripts/geneCat.pl`
3. `secScripts/MGS.pl`
4. `secScripts/MGS/strain_within.pl`

The review concentrated on concrete control-flow, resume, checkpoint, scheduler,
filename, and producer/consumer bugs. A change was made only where the intended
contract could be established from the current code, its callers, tests, or
repository history. Potential fixes that would silently change the biological
cohort, discard usable products, or add a large input scan/recalculation were
left as findings.

The review included static tracing, focused synthetic tests, repository-history
comparison, and the supplied MGS resume log. It did not execute the external
bioinformatics programs or a live HPC scheduler end to end.

## Workflow contracts after the audit

| Boundary | Enforced contract |
|---|---|
| MATAF4 -> geneCat | Binning off is represented as `doMags=0`; strains require MGS; binner is `0..5`; the downstream quality selection is exactly one checker when MGS is enabled; `SNPcaller` is preserved. |
| geneCat asynchronous stages -> MGS | Catalog publication covers both FASTA and cluster files; Canopy and functional stones are written by their actual completion jobs; an invalid eggNOG stone is removed before replacement work is exposed to MGS. |
| geneCat -> MGS | Saved `MGS.sh -outD` is reused only for the same canonical catalog path when resolvable, cluster identity, and binner. Controlled path arguments are shell-quoted. |
| MGS resume -> GTDB/BinExtr | Declared catalog-derived stage parameters and recorded output stats, rather than an unrelated controller-version fingerprint, decide reuse. Compatibility adoption is restricted to structured MGS 0.54/0.55 manifests. |
| MGS -> strain_within | The selected `MPI`/`FB` caller and exact generated arguments are forwarded. Phase-I reuse is bound to a durable input contract; split workers inherit the caller plus `taxonAwareLocusSelection` and `prepareMosaicLoci`. |

## Supplied MGS/GTDB resume log

The supplied log shows the characteristic failure sequence:

- MGS starts as version 0.55.
- It reports a Stage-I checkpoint mismatch but retains the existing MGS clusters.
- It recreates 894,144 gene sequences for 537 MGS.
- It submits a new `GTDB_MGS` job even though the taxonomy products existed.

Repository history identifies the strongest likely cause. Commit `458fc6a`
changed only the downstream strain-resume launcher while bumping MGS 0.54 to
0.55. The old global checkpoint parameter set included `pipeline_version`, so,
if the supplied run's checkpoints were written by 0.54, this field alone was
sufficient to invalidate BinExtr, GTDB, and abundance checkpoints even though
their biological inputs had not changed. The old checkpoint JSON was not in the
supplied log, so another parameter or output-stat mismatch cannot be excluded.

MGS now uses stage-specific contracts for representative-genome extraction and
GTDB taxonomy. New manifests record the declared catalog-derived parameters and
relevant outputs. Existing global-schema JSON manifests recorded by 0.54 or
0.55 are eligible for a one-way compatibility adoption only while the current
controller is 0.55 and the parameters available in that schema plus recorded
output size/mtime still match. Historical BinExtr manifests may contain
`GTDBtk.tar.gz` because GTDB was submitted before extraction outputs were
collected; migration ignores exactly that downstream archive while validating
the recorded core, gene FASTAs, representative contigs, and optional family
outputs. Old BinExtr manifests did not record MAG-report stats, so compatibility
uses the narrower temporal check that the current nonempty report is not newer
than the manifest. Empty stones, malformed manifests, older releases, missing
or changed recorded outputs, incompatible catalog-derived parameters, or
`-redoTax` still cause a rebuild.

The `/projects/8/...` command path becoming `/projects/1/...` in the log is
consistent with `/projects/8/...` being a symlink or alias whose `abs_path`
target is `/projects/1/...`; this was not independently verified. While both
aliases resolve, that does not explain the version mismatch. It is nevertheless
a separate path-resume risk because MGS canonicalizes `-GCd` but preserves the
spelling of an explicit `-outD`; see remaining flaw 9 below.

## Robust fixes applied

### MATAF4.pl

- Unknown command-line options and negative/reversed ranges now fail instead of
  continuing with an unintended configuration. The range is revalidated after
  `-to` is clamped to the map size.
- A partial/ranged/locked invocation can no longer replace the canonical
  `metagStats.txt` or emit a misleading full-catalog `GeneCat_pre.sh`; both
  require terminal records for the complete mapped cohort. An existing summary
  is retained on an incomplete pass.
- The generated geneCat handoff now preserves `SNPcaller`, explicitly separates
  `doMags` from the numeric binner, disables strains when no MGS can exist, and
  translates MATAF4's broader CheckM settings to MGS's exact-one-checker
  contract.

### secScripts/geneCat.pl

- Added and validated `SNPcaller=MPI|FB`, and forwarded it to MGS. MGS-only
  binner/checker constraints are enforced only when `doMags=1`; `doStrains=1`
  now requires MGS.
- Removed accumulated-command duplication in both clustering paths. Deferred
  script construction no longer tests outputs before the generated work has
  run, and marker extraction is no longer omitted in that mode.
- Canopy completion is owned and validated by the Canopy child. It preserves
  `canopyAutoCorr`; MGS no longer starts behind a controller-created premature
  stone.
- Functional annotation now creates its stone through a convergence job that
  depends on every database matrix. The recursive worker also preserves the
  requested FASTA split size and every `FuncMin*` threshold instead of silently
  restoring defaults.
- Invalid stale eggNOG stones are removed before resubmission, closing the race
  with MGS's existence-based wait. Chunk resume checks the output actually used
  by the combiner (`*.emapper.annotations`).
- Scheduler return values are used in their documented `(job_id, command)`
  order for functional/FOAM dependencies.
- Catalog publication validates and checkpoints both the nucleotide FASTA and
  cluster file, and restores backups if either final product is absent. A matrix
  affects Canopy sample-count decisions only with a valid matrix checkpoint.
- Multi-step children inherit the selected clustering engine and memory.
  `requireAllAssemblies=0` now also covers a completely absent optional sample
  path.
- A saved MGS output directory is accepted only for the current catalog
  (canonicalizing resolvable aliases), cluster identity, and binner. The MGS
  handoff quotes controlled paths while leaving the configured launcher command
  raw.

### secScripts/MGS.pl

- Added the stage-specific BinExtr/GTDB contracts and strictly bounded 0.54 to
  0.55 compatible-manifest reuse described above. New BinExtr manifests exclude
  downstream GTDB archives.
- A missing cheap `.core` derivative no longer forces expensive Stage-I
  reclustering when an active, preserved-unweighted, or weighted assignment is
  recoverable. Interrupted weighted handoffs are resolved before clustering,
  with dependent checkpoints/trees invalidated before any rename.
- A recovered assignment with mismatched provenance is no longer relabelled as
  current. New Stage-I checkpoints are written only for clustering performed in
  the current invocation and record both active assignments and `.core`.
  Every true rebuild removes stale derivatives first.
- Stage-I recovery now requires at least one parseable assignment rather than a
  merely nonempty file. The probe stops after its first assignment, so a
  header-only/truncated file cannot become a false valid-no-MGS outcome without
  adding a second full-file scan.
- Restored binner-local SpecI output under
  `Bin_<binner>/Annotation/<marker>_MGS` and passes that output explicitly to
  the worker, preventing the main products of different binners from colliding.
- Added/validated `SNPcaller`, forwards it to strain analysis, and shell-quotes
  each generated argument without rewriting the configured launcher command.

### secScripts/MGS/strain_within.pl

- Added strict `SNPcaller= MPI|FB` support and derives the exact compressed
  MATAF4 consensus NT, AA, contig, primary-VCF, and support-VCF filenames after
  option parsing. Caller selection is forwarded to split workers.
- Added an atomic Phase-I input contract with `building` and `complete` states.
  In version-2 state it binds catalog identity; requested/canonical path,
  device, inode, size, and mtime metadata for the catalog/index/protein/marker/
  eggNOG/map inputs and MGS-guide companions; marker/cluster identity; caller;
  and filenames. Ordinary replacement or rewriting of those bound files is
  detected without reading their contents, and a matching interrupted build
  can resume. Missing legacy state remains an unverified compatibility exception
  only for the historically hard-coded MPI caller.
- Split workers now inherit `taxonAwareLocusSelection` and
  `prepareMosaicLoci`; split and unsplit extraction no longer diverge by falling
  back to worker defaults.
- Refuses `-redo all -MGSsubset`, and refuses any input-rebuilding subset run
  against an output root with durable non-subset state. This prevents `prepRun`
  from clearing the shared output and rebuilding only a subset. A genuinely
  fresh subset output remains allowed.
- A supplied but empty outgroup tree now fails during input validation rather
  than much later in the workflow.

## Major remaining logical flaws

These are real failure paths, but an automatic change was not judged safe within
this audit.

### 1. MATAF4 partial ranges can split a shared assembly/map group

`CntAimAss` and `CntAimMap` are calculated from every member in the full map,
while the progress counters increment only after ignored and locked samples have
already been skipped. A `-from/-to` selection that contains only part of a
shared group, or a group with a permanently ignored member, may never reach its
target and therefore never submit the shared assembly/mapping work. A safe
repair needs an explicit policy: reject split groups, expand the requested
range, or deliberately redefine biological group membership. Choosing one in
code would change the user's cohort.

### 2. MATAF4 global post-processing has no single aggregate-root/cohort model

Comma-separated maps can have different output roots. The controller audit root
comes from the first sample, `baseOut` changes per sample, and global reports are
ultimately written under the last sample's root. In addition, several global
mergers (RiboFind, MetaPhlAn, mOTU, Diamond, d2, and related products) still run
after an intentional partial range. This audit protected the canonical
metagStats/GeneCat publications, but fixing every other aggregate requires a
declared common root or independent per-root completeness contracts.

### 3. geneCat checkpoints do not encode enough option provenance

Most geneCat stones are validated primarily by `cluster_id`. With the same
identity, changes to gene length/external genes, clustering engine, marker set,
functional databases/thresholds, or Canopy autocorrelation can reuse
incompatible results. Correct repair requires stage-specific manifests,
downstream invalidation, and a migration policy for legacy empty stones; it can
legitimately trigger expensive reclustering or annotation.

### 4. Parallel gene-catalog aggregation is not crash-idempotent

Each batch appends its gzip member to a shared aggregate, synchronizes it, and
then deletes the source. A crash after the append but before source deletion
causes the same batch to be appended again on retry. A correct fix needs a
durable batch ledger or staged aggregate publication, not another existence
check.

### 5. MGS deliberately retains Stage-I assignments after a genuine mismatch

When primary assignments exist, MGS warns about a checkpoint/input mismatch but
retains them unless the user requests `-redoCluster 1`. The audit stopped it
from laundering those assignments into a current checkpoint, but downstream
work can still combine old membership with changed catalog inputs. Automatically
discarding/reclustering would be expensive and would reverse the established
recovery policy.

### 6. `-ignoreIncompleteMAGs 0` cannot perform its documented recovery path

The per-assembly binner-generation loop begins with an unconditional `last`, so
turning off `ignoreIncompleteMAGs` does not generate missing bins and the later
validation can fail. Enabling this branch would launch substantial binner and
quality work and currently contains binner-specific unsupported branches; it
needs a separate functional/HPC test before activation.

### 7. Strain Phase-I provenance still does not cover every semantic input

The version-2 contract detects ordinary metadata changes to its catalog, map,
and MGS-guide inputs, but it does not hash their contents and does not stat-bind
the per-sample consensus FASTA/VCF/depth sources, Mosaic catalog/raw assignments,
or every Phase-I QC, locus, and consensus tuning flag. A same-path/inode/size/
mtime rewrite, a changed per-sample source, or changed tuning with the same bound
metadata can therefore reuse earlier extraction. Expanding this safely requires
a compatibility policy; hashing every large source on each lean resume would
also violate the runtime constraint.

### 8. GTDB taxonomy checkpoints omit tool/reference-release provenance

The stage contract binds catalog-derived inputs and output stats, but not the
configured GTDB-Tk command, `GTDBtk_DB` release/path, or skani reference. A tool
or reference-database upgrade can therefore reuse taxonomy produced by the old
release. The legacy manifests contain no trustworthy release token, so adding a
strict token now would either adopt unverifiable history or force the expensive
GTDB rerun that this audit was asked to avoid.

### 9. Catalog aliases can split input and output checkpoint namespaces

MGS canonicalizes `-GCd` with `abs_path` but converts an explicit `-outD` only
to an absolute spelling. This is visible in the supplied log: catalog inputs
are announced under `/projects/1/...`, while outputs and manifests remain under
`/projects/8/...`. Recorded output paths can become unresolvable if that alias
changes or disappears, causing avoidable rebuilds even when the data remains.
Canonicalizing explicit output roots automatically could also redirect a user's
requested publication location, so it was left as a documented path contract.

### 10. Strain input rebuilds mutate shared guide derivatives

Every `-onlySubmit 0` strain run can delete and recreate `$MGSfile.srt*` beside
the shared MGS/core guide, including runs using a separate `-outD` or an MGS
subset. A concurrent strain controller or split worker using the same guide can
therefore see `.srt` or `.srt.gene2MGS` disappear between checks. Output-root
guards cannot protect files outside that root. A safe correction needs shared
locking/freshness rules or output-local derivatives, not another existence
check.

## Other concrete remaining risks

- MGS Canopy preparation accepts cached `.filt`/quality products largely by
  presence; a changed Canopy source can therefore reuse stale filtering/QC.
- MGS `-submit 0` is not a pure dry run: direct Stage-I commands may still
  execute. Do not use it as a no-write inspection guarantee.
- In geneCat fire mode, the outer script can remove its temporary directory
  before exceptionally slow nested functional jobs converge. Serial waiting
  would lengthen the critical path; safe repair needs separate temp ownership.
- A strain split worker can fail late when primary consensus inputs exist but a
  required PB/ONT support VCF/depth file is absent. Whether to skip that sample
  or regenerate without support changes biological semantics.
- Strain `-redo input` currently regenerates inputs for every selected MGS, not
  only incomplete ones, and can do more work than its name suggests.
- Missing-contract legacy MPI strain output cannot prove which old catalog/guide
  produced it. It remains accepted to avoid forcing all historical runs through
  a full rebuild; the first successful accepted resume promotes that unverified
  state into a version-2 contract recording the adopted current metadata, not
  proof of its original provenance.
- The fast fingerprints use filesystem identity/size/timestamps rather than
  hashing very large catalogs on every resume. Carefully engineered same-stat
  in-place rewrites can evade them.
- MGS's generic `wait4stone` checks existence rather than producer parameters.
  The active eggNOG path is protected by producer-side invalidation, but the
  generic option cannot validate an arbitrary external producer.
- The optional Perl `clusterMAGs` compatibility output appears inconsistent
  with `markersPerMGS.pl`'s expected report columns. The binary path is the
  default; the compatibility schema needs a coordinated producer/consumer test.
- Although binner-local SpecI main products are now isolated, the SpecI worker
  still emits a catalog-global legacy `gene2specI.txt` sidecar, so simultaneous
  binner runs can race that compatibility file.
- Many older generated commands still interpolate paths without quoting. A
  global rewrite is unsafe because some configured executable strings
  intentionally include environment setup or shell composition.

## Runtime impact

The added resume checks use manifest parsing and constant-cost stats per bound
file, plus validation of output records already enumerated by checkpoints; the
total check is proportional to the number of bound files and manifest records.
No whole-catalog hashing or new biological analysis was added. Eligible
0.54/0.55 MGS runs can now avoid the expensive BinExtr/GTDB repetition seen in
the supplied log. Some newly enforced dependency waits can make the controller
wait for work it previously (incorrectly) treated as complete, but they do not
add a new analysis stage.

## Verification record

- Baseline before edits: full `prove -lr t` passed 52 files / 1,887 tests.
- All four audited Perl entry points compile with `perl -I. -c`.
- Integrated checkpoint, path, CLI, geneCat, MGS, strain, stability, and
  workflow-control suites pass: 10 files / 633 tests.
- Final full `prove -lr t` passes: 53 files / 1,955 tests.
- `git diff --check` passes.
