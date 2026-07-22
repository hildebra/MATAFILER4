# Changelog

## 2026-07-22 — MGS and strain workflow hardening

This release hardens the Perl orchestration around species clustering, taxonomic abundance estimation, and between- and within-species phylogeny construction. The `clusterMAGs` C++ clustering algorithm itself was not reviewed in this pass.

### Version updates

- `secScripts/MGS.pl`: 0.34 → 0.35
- `secScripts/MGS/clusterMAGs.pl`: 0.21 → 0.22
- `secScripts/MGS/strain_within.pl`: 0.38 → 0.39
- `secScripts/MGS/strain_within_2.2.pl`: 0.31 → 0.32
- Earlier changes in this review already raised `buildTree5.pl` to 5.08, `annotateMGwSpecIs2.pl` to 0.14, and `resortMGSgenes4importance.pl` to 0.15.

Scripts without an established runtime version variable retain that convention and now carry a dated source comment describing their changes.

### Changed and fixed

- Restored GTDB markers as the documented default and consistently propagated marker choice, core count, quality-checker mode, and catalog identity.
- Added checkpoint manifests with workflow parameters and input size/mtime fingerprints. Changed inputs now invalidate stale clustering and phylogeny products.
- Hardened empty, singleton, and sparse MGS outcomes, weighted resumes, representative-bin extraction, resource requests, and downstream output validation.
- Made the `clusterMAGs.pl` wrapper validate inputs, honour its fallback mapping input, create log paths safely, and clean temporary outputs using checked operations.
- Deduplicated marker-to-MGS links and excluded markers shared by multiple MGS from abundance and taxonomy correction, avoiding input-order-dependent assignment.
- Made taxonomy parsing reject malformed or duplicate records, retain zero-marker MGS in statistics, normalize empty ranks, merge bacterial and archaeal GTDB-Tk summaries under one header, and isolate each GTDB-Tk attempt.
- Made between-MGS phylogeny construction deterministic, skip sets with fewer than three marker-bearing taxa, and exclude ambiguous paralog cells instead of selecting an arbitrary first copy.
- Made within-MGS split jobs generation-aware, require exact worker-part completion, merge parts atomically, isolate logs, validate paired nucleotide/amino-acid consensus files, and protect output paths from unsafe deletion.
- Added safe in-place handling for compressed outgroup FASTA files, rejected ambiguous sidecars, validated tree outputs, and allowed the biologically valid three-tip minimum.
- Made downstream strain summaries header-safe and dry runs non-destructive, with dependency and output checks for network and treeWAS analyses.
- Restored partition-aware tree inference, activated per-locus overlap filtering with synchronized partitions, and required non-empty resume artifacts in `buildTree5.pl`.
- Fixed final-record subset FASTA parsing, representative-family EOF handling, deterministic gene ranking, and marker identifier trimming.

### Biological interpretation and follow-up

- Within-species trees represent one dominant haploid consensus per sample and do not resolve multiple coexisting strains in one sample.
- Recombination masking is not yet a default stage, so branches may reflect recombination as well as vertical ancestry.
- A between-MGS pseudo-tip can combine markers from different MAGs assigned to the same MGS; representative-genome or explicit locus-support policies should be evaluated.
- Arithmetic mean marker abundance remains sensitive to marker dropout and outliers. Median or trimmed estimators plus marker-support reporting would improve interpretation across uneven coverage.
- Fully side-effect-free MGS dry runs, workflow locking, and complete shell-argument quoting remain follow-up engineering work.

### Validation

- Focused regressions cover checkpoint invalidation, taxonomy merging, shared-marker exclusion, sparse phylogenies, split-generation recovery, compressed FASTA handling, and tree partition/overlap behaviour.
- At the end of the implementation pass, 349 focused assertions and syntax checks for 16 Perl scripts/modules passed.
