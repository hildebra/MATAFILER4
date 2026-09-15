# geneCat audit follow-up — 11 September 2026

Scope: `secScripts/geneCat.pl`, the directly affected catalogue helpers and their command registrations. `secScripts/MGS.pl` is excluded. Changes to other pipeline controllers already present in the shared checkout are outside this audit.

## Original findings

| Finding | Result |
|---|---|
| 1. Backups omit marker representatives | Pre-merge clusters use separate `unmerged.*.gz` backups. Final backups are written after marker merging, and the merged checkpoint validates them. A published catalogue remains resumable after temporary backup cleanup. |
| 2. Interrupted collation duplicates batches | A full collation retry clears the owned B0 aggregates, temporary marker batches, stale append locks and batch reports before replaying all batches. |
| 3. Incorrect SAM supplementary-alignment mask | No independent flag-filter change was made. The parser had only the removed multi-step caller and was deleted with that workflow. The terminology clarification below explains its scope. |
| 4. Repeated cluster conversion fails | Conversion accepts both original and already-numbered headers, regenerates the member index, and publishes replacement files only after validation. |
| 5. Multi-step clustering expects uncompressed inputs | Eliminated by removing multi-step clustering. The single-step workflow consumes the compressed collation outputs. |
| 6. Temporary path requires a trailing slash | Marker merging joins its paths with `File::Spec`; the controller normalizes the temporary directory. |
| 7. External representatives cannot be numbered/extracted | Nucleotide rewriting handles the existing `xtraSmpls` index. External proteins reuse `attachProteins3`, including strict validation that requested proteins exist. |
| 8. Marker extraction assumes identity 95 | The selected identity is passed to marker extraction, k-mer statistics and decluttering helpers. Optional positional identity arguments retain 95 as the default for existing standalone helper invocations. |
| 9. No Kraken classifications fails under pipefail | A configured `awk` filter accepts zero matching records. Successful empty output is published; a failed Kraken process is still rejected. |
| 10. Configured eggNOG command is ignored | Annotation workers execute the complete configured `emapper` command, including wrappers. |

## Supplementary sequencing inputs versus SAM alignment flags

MATAFILER supplementary data are additional sequencing files associated with a sample, such as PacBio data supplied alongside Illumina data. Their primary alignments remain primary alignments. Input-library provenance does not, by itself, set the SAM supplementary-alignment flag.

SAM's supplementary-alignment flag instead identifies an additional alignment segment of a split read. Its value is decimal 2048, written as hexadecimal `0x800`. The removed parser used `0x2048`, a different bit mask. The original audit finding concerned this parser's intermediate gene-to-gene minimap2 SAM files, not the primary alignments of MATAFILER supplementary sequencing inputs.

The parser has been removed only because its multi-step clustering caller was removed. The `-calcSupplCovSmpls` option and its forwarding to abundance-matrix generation remain in place. Regression fixtures exercise supplementary coverage both enabled and disabled.

## Consolidation and configuration

- A shared command builder produces the three abundance-matrix variants, and one output-aware completion contract is used in local and deferred execution.
- The multi-step clustering implementation, its unused supporting routines and `secScripts/unused/mergeCls.pl` are removed. Existing commands should omit the removed `-1stepClust` option. Marker-specific clustering remains part of the single-step catalogue workflow.
- Direct external commands in geneCat resolve through `getProgPaths`; shell builtins and native Perl file operations need no executable registration. Checkpoint writing and recursive geneCat calls also use configured script commands.
- `Mods/config_internal.txt` registers the shell utilities still used by generated jobs, `writeCheckpoint_scr`, `cdhit_est`, `kaiju`, and `kaiju_addTaxonNames`. Override the new executable keys in site configuration when programs live outside PATH; changing only the old `cdhit` or `kaijuDir` entry no longer selects these executables.
- Perl helper registrations explicitly select Perl. `hmmBestHit_scr` selects Python 2 because the existing legacy helper uses Python 2 syntax; it has not been ported in this change.
- The Kaiju name-expansion default is `kaiju-addTaxonNames`, as documented by the [Kaiju project](https://github.com/bioinformatics-centre/kaiju#adding-taxa-names-to-output-file). A legacy installation can override `kaiju_addTaxonNames` with its older executable.

## Validation

The focused recovery suite runs actual production helper bodies and generated shell commands with small fixtures. Clustering/scheduler and annotation executables are mocked; marker merging, gzip backups/restoration, file publication, external protein selection and marker FASTA extraction use the actual implementations.

Covered cases include local and deferred clustering; 40 marker representatives surviving backup restoration; replay without duplicate collation records; repeat cluster conversion and missing member-index recovery; external FASTA descriptions; missing required proteins; catalogue identity 97; configured command wrappers; successful empty Kraken output; and failed Kraken output remaining unpublished.

The compatibility suite covers matrix generation, workflow invariants, configuration and catalogue paths, checkpoints, compressed I/O, program/file checks, command-line documentation and shared-helper stability. Perl syntax and scoped whitespace checks pass.

No real biological catalogue, external clustering/annotation installation, or HPC scheduler was exercised. These fixtures verify controller behavior and command construction, not biological performance of external tools.
