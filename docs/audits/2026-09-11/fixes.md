# Follow-up: all 11 local algorithm defects fixed

All 11 findings in [the audit](report.md) have been corrected. The original numerical fixtures are retained in [results.json](results.json); the same reproducer now records the corrected behavior in [results.after-fixes.json](results.after-fixes.json).

## Changes by finding

| Finding | Implemented correction |
| --- | --- |
| 1. Secondary coverage filenames | `coverage_derivative_paths` and `coverage_derivatives_complete` in `Mods::GenoMetaAss` now give scheduling, completion checks, and sample evidence the same naming convention. Actual `rdCover` outputs are accepted; historical `.coverage.gz.<suffix>` names and compressed derivatives remain supported. Forced remapping also clears the actual derivative filenames. |
| 2. mOTUs across samples | Cached taxonomy parsing is separate from abundance accumulation. Every sample contributes to every rank, including mOTUs encountered earlier. |
| 3. Functional subject lengths | Every candidate uses its own validated positive integer subject length for coverage and CARD thresholds. GLN uses the selected subject's length. Missing or invalid lengths fail with a subject-specific diagnostic. |
| 4. Mate-2-only hits | Either mate survives alone. If mates hit different subjects, their singleton candidates compete in the existing best-hit selection; they are no longer discarded by an accidental intersection. |
| 5. Empty ribosomal profiles | The merger accepts the producer's zero-byte completed-empty representation and header-only hierarchies. Columns come from the complete input sample roster. An all-empty cohort emits a header containing every sample and no taxon rows. Malformed nonempty headers still fail. |
| 6. Mapping statistics | `Mods::StatsLogReader::parse_bam_filter_counters` reads current and legacy counter blocks and sums complete process blocks. Retained/input percentage includes already-unmapped and malformed records in the denominator. Bowtie's input total is preserved. Strobealign no longer reports a fraction of one or invented uniqueness. |
| 7. Functional multiplicity | Ordinary reads and genes count once; a combined pair counts twice. Raw-read searches append `MF4:read_count=2` to merged-library hits before mixing files; ordinary hits use the parser’s count-one default, avoiding redundant compression passes. Catalog annotation passes `-queryType genes`. The parser also supports explicit homogeneous `-queryType merged` input. |
| 8. Paired overlap | Both functional parsers use `Mods::FuncTools::mergeBlastPair`. Inclusive intersections handle containment, identical intervals, shared endpoints, disjoint intervals, and reversed coordinates. Original hit records are left intact. |
| 9. Gene k-mer windows | One bounded-buffer emitter handles steady-state output, contig transitions, and EOF. Each gene appears exactly once with a centered, boundary-truncated window. It reuses `avgArray` and `roundAr`. |
| 10. Uninformative final k-mer | Header boundaries and EOF use the same `emit_kmers` routine. Eligible sequences without a valid k-mer are omitted consistently, preserving the former non-final-record policy. |
| 11. GC / GC3 rows | One row writer always supplies the newline and the `-1` undefined value. Base counting normalizes lowercase sequence. Contig GC3 uses exact contig keys, so names such as `ctg` and `ctg2` stay separate. |

## Interpretation and compatibility

- **Gene windows:** the existing `pm5` output naming is interpreted as the gene plus up to five neighboring genes on each side: at most 11 genes. At an edge, average only available genes from the same contig. The second argument is a radius; zero reports the gene's own vector. The old routine mixed asymmetric windows with duplicated records, so its numerical output is intentionally not preserved.
- **Paired alignments:** ungapped overlap length is exact. BLAST's 12 summary fields cannot reconstruct shared alignment columns for gapped hits. The shared helper estimates overlapping columns from the smaller alignment-column density over the subject interval, rounds to an integer, and caps the subtraction at each hit's length. Mean identity, minimum E-value, and proportional combined-bit-score scaling retain their existing interpretation.
- **Counts versus GLN:** count multiplicity represents the number of reads supporting the selected candidate. GLN remains aligned length divided by subject length, with paired overlap removed; it is not multiplied again for a merged pair.
- **Subject-coverage option:** `-minFractQueryCov` remains a compatibility alias for `-minPercSbjCov`. The main caller now uses the accurate parser option name. MATAF4's public `-DiaMinFracQueryCov` flag retains its historical subject-coverage behavior, now described correctly in the current flag documentation. True query coverage would require a separate interface and query-length data.
- **Mapping table units:** historical column names are retained. `AlignedReads` is the retained mapped **SAM record** count, and `OverallAlignment` is its percentage of all input SAM records when filter counters exist. Secondary/supplementary alignments can make record counts differ from unique reads. `ReadsPaired` preserves Bowtie's input units (pairs for paired input; reads for single-end input), uses minimap2's processed sequence count when available, and otherwise falls back to filter input records. Bowtie's category rates remain aligner statistics. Missing or incomplete filter blocks do not become fabricated zero counts or 100% rates.
- **Cached functional searches:** new raw search files have a `.read-counts-v1.stone` marker. Reparsing preserves it, avoiding redundant alignment. If the functional caller encounters an old merged-library search without provenance, it regenerates those hits before interpretation because the missing multiplicity cannot be recovered from a bare query identifier.

The code changes do not retroactively rewrite completed biological result tables. Recompute affected stages to update existing results. For functional profiles, the existing `-reParseFunct 1` option requests reinterpretation; merged-library searches lacking provenance will need alignment again. The broader `-reProfileFunct 1` option requests fresh functional searches. No production dataset was rerun as part of this change.

## Validation

- Added `t/local_algorithm_regressions.t` and `t/functional_provenance.t`: **136 checks passed**. These exercise actual Perl scripts, extracted MATAF4 caller bodies, the bundled `rdCover` binary on Linux, and generated search/interpretation commands with a fixture aligner.
- Additional cases cover CARD length-scaled thresholds, invalid subject lengths, mixed library provenance, reparse cache reuse, reversed/contained/endpoint/gapped overlaps, both empty hierarchy representations, an all-empty cohort, multiple window radii and contig sizes, lowercase GC/GC3, and current/legacy/multiple filter blocks.
- The initial six affected existing test files passed **311 checks**.
- Main-script syntax and scoped whitespace validation passed.
- Full current-checkout regression suite: **64 files, 2,689 checks passed** with `prove -I. -It/lib t/*.t`.

External aligner/assembler algorithms, HPC execution, and biological accuracy on real datasets were not tested. Other catalog/MGS/strain changes were already present or appeared concurrently in the shared checkout and were not part of these 11 fixes.
