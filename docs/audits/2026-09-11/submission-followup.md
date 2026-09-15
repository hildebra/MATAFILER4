# Requested follow-up to the submission dependency audit

The follow-up implements the user's conservative fallback for finding 1 and fixes findings 2 and 3. Scaffolding code remains available for backward compatibility; no legacy helper was removed.

## Finding 1: prevent unsafe binning

The terminal-empty assembly-group release branch no longer calls `submitGenomeBinner`. It can still release assembly, gene prediction, and deferred mapping jobs, and cleanup retains those dependencies.

A robust full ContigStats chain is not a small extension of this branch. `runContigStats` obtains raw-read metadata from the current sample and expects per-sample coverage products; the release member is empty. Final hybrid support mappings may also remain unpublished. Selecting a different member's context would additionally require coordinated deferred statistics/variant release and completion-sentinel recovery. The authorized fallback avoids that larger change and prevents binning from racing or omitting its inputs.

This deliberately **does not claim to repair automatic ContigStats or variant release for terminal-empty groups**. Such groups can remain incomplete; the unsafe binner is withheld. The ordinary group's existing statistics → binning dependency chain is unchanged.

## Finding 2: preserve hybrid exclusion accounting across completed-sample early returns

After accepting a completion record, the main loop now counts empty-input, too-small, and cleaned-empty terminal members in `CntPreAssNoPrim` when assembly mode 5 is selected. This happens before the completed-sample early return, which otherwise bypasses `prepPreAssmbl`.

The existing group reset clears this counter between passes. Open samples continue to use `prepPreAssmbl`; the closed and open paths are mutually exclusive within a visit. This adds no packaging, file inspection, or scheduler query.

The [updated probe](submission-release.after-followup.json) executes the actual completed-sample branch, then the real preassembly helper with a complete package. It now reports:

- Two visited members and two target members.
- One complete preassembly package and one excluded empty member.
- Final hybrid readiness enabled.

The regression also checks all three accepted empty outcomes, a successful nonempty outcome, an SDM-warning outcome, and ordinary assembly mode.

## Finding 3: retain legacy scaffolding with explicit dependencies

`scaffoldCtgs` now combines its supplied dependency with the existing group raw-read staging barrier, `UnzpDeps`. This applies to both ordinary and external scaffolding, including the case where the assembly already exists and no assembly producer job is pending.

`metagAssemblyRun` returns the external scaffold and gap-filling consumer IDs. Both main-loop callers add these to sample dependencies immediately, before a producer-wave early exit. Ordinary scaffolding remains part of the assembly publication barrier. External scaffolding stays separate from that barrier, preserving concurrency for unrelated assembly mapping.

The legacy fixtures exposed two related preparation cases handled in this follow-up:

- No eligible scaffold result means no downstream gap-filling job is submitted.
- A fresh GapFiller output directory is created with the existing `make_path` helper before writing its options file; the former single-level `mkdir` failed when the parent directory did not exist.

Generated commands still use the legacy mate-pair libraries, BESST output, and GapFiller options. The existing `operationMode=scaffold` entry-point guard (`die "update scaffold"`) remains unchanged; this work retains and tests the underlying compatibility code without enabling a guarded operation mode.

## Validation and runtime

- `t/main_submission_graph.t`: 53 checks, including the actual terminal-empty branch, completed-sample early return, and hybrid helper.
- `t/legacy_scaffolding_dependencies.t`: 22 checks, exercising actual controller/helper bodies and generated command arguments with recording submission stubs.
- Focused dependency suite: **4 files, 276 checks passed**.
- Main controller syntax and scoped whitespace checks passed.
- Full suite: **68 files, 2,779 checks executed; 3 failed**, all in `t/program_and_file_checks.t` (checks 28, 33, and 36). Every other test file passed, including the changed dependency tests and the prior algorithm/provenance tests.
- Those three failures compare literal IQ-TREE source strings against `Mods/phyloTools.pm`, which already had unrelated working-tree edits when this follow-up started and was not changed here. The current code quotes alignment/partition paths with `_shellQuote` and selects the model through a variable, so the old literal-source patterns no longer match. These assertions and that unrelated module were left unchanged. The suite output is retained locally in `/tmp/mataf4-group-followup-full-suite.txt`.

These changes add no scheduler queries or broad filesystem scans. The ordinary path gains a small return-value check; the completed-empty path gains an in-memory counter update. Added staging-dependency work and directory creation run only inside legacy scaffolding helpers. The earlier coverage-check and DIAMOND collection optimizations are retained.

No production scheduler jobs or external biological tools were run. Concurrent unrelated working-tree changes were preserved.
