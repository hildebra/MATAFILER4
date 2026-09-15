# BuildTree5 audit and consolidation — 2026-09-11

The implementation changes advance BuildTree from 5.92 to 5.93. The review
followed input staging, alignment preparation, checkpoint reuse, tree-engine
wrappers, EPA publication, and optional selection analysis. It also examined
`fasta2phylip.pl`, `fubarDNDS.py`, the staged strain-input finalizer,
`Mods/PhyloAlignment.pm`, and the locus-QC interface. Native MSAfix remains the
production locus-QC implementation; the older Perl QC script is not substituted
into the production workflow.

## Corrected defects

| Area | Defect and correction |
|---|---|
| Single-locus trees | The no-tree guard rejected every newly prepared single-locus alignment. Valid alignments now reach inference, with actual sample/locus counts. Normalization writes a new scratch file and renames it, avoiding modification of an input staged through a symlink. |
| Single-locus resume | Durable completion previously bypassed input and downstream-policy comparisons. Single-locus runs now use the existing policies, without requiring category-only QC reports. Unchanged completion still takes the cheap reuse path. |
| FastTree / VeryFastTree | Engines wrote directly to the final tree path; a failed run could leave a nonempty partial tree. Both now use one wrapper that requires successful execution and a nonempty temporary output before renaming it into place. |
| EPA publication | Fresh placement, EPA-only recovery, and filter-only redo now share the same mapping, filtering, report, and publication implementation. Redo no longer deletes the previous published tree before a replacement is ready. |
| EPA grafting | A tolerated endpoint overshoot could be clamped outside the grouping tolerance without consuming any placement, indefinitely repeating the loop. The loop now always consumes its first validated placement. Unexpected/duplicate query names, reference-tip collisions, and negative pendant lengths are rejected. |
| Tree labels | `getTreeLeafs` now reuses the existing Newick tip parser. Internal support/root labels are no longer mistaken for samples, and quoted tip names are handled. |
| PHYLIP conversion | The previous header omitted the final sequence's length and mixed 10- and 50-character name widths when padding sequences. Conversion now preserves aligned sequence lengths, validates all rows before stdout publication, and rejects identifier collisions. The RAxML caller publishes through a temporary file. |
| Selection analysis | An empty outgroup matched every sample; matching is now exact and conditional on a nonempty name. Terminal TGG was incorrectly removed, and some stop codons were stripped instead of masked. In-frame stops are masked without shortening alignments. Small loci return normally instead of jumping out of the caller's loop. The inferred primary tree is used when no separate tree is supplied. |
| Selection results and cleanup | Full, unique, and subsampled summaries share one loop. The older parser's leading tab no longer shifts report columns. Failed parser commands stop publication. HyPhy cleanup removes its own temporary directory, preserving the alignment workspace needed by later stages. |
| Staged categories | The legacy category finalizer rejects contradictory locus/sample identifiers and invalid outgroup overlays, matching checks already present in the shard finalizer. |
| Failure state | END processing preserves the process exit status. Placement failures record their reason before publishing pending workflow state. |

## Consolidation and compatibility

- BuildTree itself is 340 lines shorter (8,337 to 7,997); the four changed
  implementation files together are 396 lines shorter. Tests are additional.
- Three EPA publication implementations, three selection-summary loops, repeated
  QC-report backup/restore code, FastTree wrappers, IQ-TREE model selection, and
  tree-tip parsing were consolidated.
- Unreachable CodeML, XML FUBAR, Watterson-theta, synonymous-AA, PhyML, and
  Guidance implementation blocks were removed. These blocks were not reachable
  through the current supported options; exported module functions were retained.
- Existing artifact names, completion/state schemas, policy keys, and EPA report
  columns are preserved. No blanket version-based checkpoint invalidation was
  introduced. A single-locus run with a changed input/policy correctly rebuilds;
  an unchanged completed run does not rerun inference.
- The corrected PHYLIP output and FUBAR column alignment intentionally differ
  from malformed legacy output. Invalid unequal-length alignments and ambiguous
  truncated PHYLIP identifiers now fail explicitly.
- No new all-locus alignment scan, input-content hash, inference call, or
  scheduler job was added to the normal multi-locus path. The single-locus path
  reuses its normalization read for validation and sample counts and now enters
  the existing inventory reporting. Tree publication adds temporary output and
  rename operations. Runtime on a full HPC dataset was not benchmarked.

## Outstanding legacy findings

These are findings, not repaired or validated features in this change. They do
not run with the default empty `-postFilter` and the IQ-TREE settings in the
reported strain command. Their repair requires separate treatment of optional
biological processing or backend support, beyond the consolidation above.

1. **Zorro executes before the alignment exists.** BuildTree calls `filterMSA`
   while constructing the command, before executing the aligner. Unlike the
   MACSE branch, its Zorro branch immediately reads/modifies the output alignment.
   `zorroFilter` also calls `writeFasta` without importing it. Its missing-data
   masking and back-translation contract need to be settled together with the
   execution-order repair (`buildTree5.pl`, alignment loop;
   `Mods/phyloTools.pm`, `zorroFilter` and `filterMSA`).
2. **MACSE's refined alignment is not consumed.** The helper writes an `.2`
   output while later stages continue using the original alignment. The caller
   supplies an amino-acid alignment, although `refineAlignment` expects a
   nucleotide alignment, and the output-type flag is reversed relative to
   `useAA4tree`. Repair needs synchronized NT/AA outputs and an explicit decision
   about frameshift handling. See the authors' [refineAlignment
   documentation](https://www.agap-ge2pop.org/refinealignment/) and
   `Mods/phyloTools.pm::filterMSA`.
3. **Legacy RAxML capabilities do not fully match the shared options.**
   `runRaxMLng` reads the bootstrap option but never adds a bootstrap/support
   command; requesting bootstraps therefore still yields its ordinary best
   tree. Both RAxML wrappers eagerly resolve unused auxiliary tools, and the
   v8 bootstrap-recovery path retains unchecked shell operations. These paths
   need backend-specific integration coverage before stronger guarantees can
   be made (`Mods/phyloTools.pm::runRaxMLng` and `runRaxML`).

## Validation

The following suite passed: **10 files, 790 top-level tests**. The new audit
file also contains nested assertions exercising the defects above.

```sh
PERL5LIB=.:t/lib PERL5OPT=-MMFTestConfig prove \
  t/buildtree_audit.t t/build_tree5_regressions.t \
  t/build_tree5_taxon_aware_smoke.t t/iqtree_output_validation.t \
  t/strain_placement.t t/finalize_strain_tree_inputs.t \
  t/phylo_alignment.t t/post_alignment_locus_qc.t \
  t/strain_workflow_regressions.t t/strain_parts.t
```

Perl syntax checks passed for all four edited implementation files; the FUBAR
Python parser also compiled. Focused whitespace checks passed. Tests cover
resume/rebuild decisions, retained artifacts, failed publication, placement
filter redo, query validation and grafting, selection-input preparation,
PHYLIP conversion, and filenames with spaces/quotes in the changed wrappers.
External inference/selection programs are controlled fixtures where needed;
this is not an end-to-end validation with production IQ-TREE, EPA-ng, HyPhy,
MACSE, Zorro, or RAxML binaries and an HPC scheduler.
