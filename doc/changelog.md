# Changelog

## 2026-08-19 — Post-QC recovery of partial strain loci

- Updated `buildTree5.pl` to 5.77 and `strain_within.pl` to 1.25. `-GeneLengthMin` remains the high-confidence per-sample locus gate used for candidate selection and backbone/placement eligibility; the new `-GeneLengthIncludeMin` (default `0.03`) recovers lower-coverage observations only after that QC gate.
- Recovered observations enter both the complete/backbone MSA and, when strict placement is explicitly enabled, the placement-query MSA. Final taxon-aware sample coverage metrics ignore recovered-only observations, while MSAfix masking and post-alignment occupancy/divergence QC evaluate the expanded alignments.
- Added per-MGS `phylo/gene_length_filter.samples.tsv`, run-wide `LOGandSUB/strainGeneLengthFilter.samples.tsv`, and separate selection-attrition totals for observations dropped by either threshold or recovered into MSA input.
- Strict-backbone EPA placement remains disabled by default in both scripts (`-strictBackbone 0`) and is still enabled only by an explicit user flag.

## 2026-08-06 — Legacy strain inference and placement eligibility

- Updated `buildTree5.pl` to 5.36 and `strain_within.pl` to 0.80. The historical IQ-TREE command is again the strain default; `-legacyMGTK 0` selects the standard modern command, while `-iqPathogen 1` selects pathogen mode and disables the default automatically.
- The relaxed 35%-of-Q90 backbone continues to maximise samples used for initial inference. A sample deferred for low coverage is now grafted only if it passes the existing `GenesPerSpecies` and relative-NT thresholds calculated from the final selected alignment, the configured NT/overlap floors, and an absolute minimum of two loci. This prevents single-locus re-additions.
- Added `phylo/taxon_aware_placement_eligibility.tsv`; `strict_backbone.samples.tsv` now records samples excluded from placement and their reason. The continuation fingerprint records the legacy-IQ-TREE choice.

## 2026-08-06 — Site-driven strain partition merging

- Updated `buildTree5.pl` to 5.35 and `strain_within.pl` to 0.79. Initial deterministic rate/GC partition count is now driven by effective called sites rather than retained locus count: `ceil(total effective sites / 30000)`, capped at eight bins by default.
- The splitter repeatedly refines the largest current bin at an effective-site-balanced boundary using the more heterogeneous local signal of post-alignment P90 divergence or GC fraction. Existing 20-locus and 20,000-effective-site minimums still merge weak bins into their nearest rate/GC neighbour.
- Added `-rateMergeTargetSites INT` (default `30000`) to both `buildTree5.pl` and `strain_within.pl`; it is recorded in the continuation-policy fingerprint.

## 2026-08-05 — Deterministic strain partition merging

- Updated `buildTree5.pl` to 5.34 and `strain_within.pl` to 0.78. Strain trees now deterministically pool final nucleotide loci into bounded rate/GC partitions while retaining fixed `GTR+F+G2`; direct `buildTree5.pl` calls can opt in with `-rateMergePartitions 1`.
- Initial partition capacity scales to 4, 6, or 8 bins for up to 100, 250, or more loci. Bins below 20 loci or 20,000 effective called sites merge into their nearest normalized divergence/GC neighbour; taxon-rescue loci join the nearest robust bin.
- The rate proxy reuses post-alignment P90 consensus divergence, with taxon-aware variable-site fraction as a fallback. `phylo/rate_merged_partitions.tsv` audits all metrics, coordinates, and assignments, and the continuation fingerprint includes every merge setting.

## 2026-08-05 — Fast fixed strain-tree model

- Updated `buildTree5.pl` to 5.33. The strain preset again uses fixed `GTR+F+G2` by default because full `MFP+MERGE` model and partition selection is too expensive for hundreds of loci. `-AutoModel 1` remains available for an explicit diagnostic run.

## 2026-08-05 — Strain locus outlier QC and partition merging

- Updated `buildTree5.pl` to 5.32. Within-species MSAfix locus QC now uses a modified-Z cutoff of 5.0 (previously 8.0) for consensus-divergence outliers, rejecting anomalous loci while leaving sample retention unchanged.
- The strain preset now defaults to IQ-TREE AutoModel. Partitioned nucleotide trees use `MFP+MERGE`, allowing IQ-TREE to merge compatible locus partitions instead of estimating an unstable independent rate for every short locus. `-AutoModel 0` remains an explicit fixed-model opt-out.
- The continuation policy records the AutoModel choice, so a cached result inferred under a different tree-model policy is rebuilt rather than reused.

## 2026-08-05 — Taxon-aware strain locus selection

- Updated `buildTree5.pl` to 5.31 and kept the selector localized in the tree workflow rather than splitting MSA orchestration into another script.
- Keep IQ-TREE inference unrooted: `-outgroup` remains available for locus anchoring and downstream backbone/placed-tree rooting, but is no longer passed as IQ-TREE's cosmetic `-o` flag.
- Added taxon-aware locus selection and subsequently made `-taxonAwareLocusSelection 1` the default. Direct calls align a bounded 400-core/100-rescue set, carry 150 extra loci as QC backfill, and choose the final 500 loci after MSAfix from observed occupancy and parsimony-informative sites.
- Replaced relative species rejection in the taxon-aware path with an anchor rule: well-covered samples are backbone candidates, sparse samples remain placement candidates, and only samples lacking usable placement overlap are removed.
- Added candidate/final locus and sample TSV audits under the tree output directory, and fingerprinted selector settings in continuation policy state.
- Updated `strain_within.pl` to 0.77. Taxon-aware selection is enabled by default and its 80% core, 20% rescue, and 30% additional backfill hierarchy scales to the effective `maxGenes`/`presortGenes` budget. The other defaults remain gene-length fraction 0.3, per-species gene fraction 0.05, and relative informative-NT fraction 0.02.

## 2026-08-04 — MATAFILER4 4.38 orchestration consolidation

- Route normal, SDM-warning, and empty/too-small completion through one sentinel-publication path with shared metagStats and input-size data.
- Reuse one parsed scheduler queue snapshot for live counts, sample locks, and immediate dependency waits; share bounded Slurm accounting batches between dependency reconciliation and failure summaries.
- Generate primary and supplementary SNP products from one scope model covering validation, region planning, pileup, concatenation, normalization, vcf2fna inputs, and completion stones.
- Centralize BuildTree alignment recovery and per-engine checkpoint state so completion tests and execution decisions cannot diverge.
- Remove the obsolete pending-only scheduler counter and unused SNP locals, and regenerate structural regression assertions for the consolidated paths.
- Limit SortMeRNA reference, optional-index, and downstream ribosomal-database preflight validation to path existence, leaving access validation to the tools on the allocated node.
- Build RiboFinder SortMeRNA commands with per-invocation work directories and place Boolean flags before value-bearing options, avoiding ambiguous trailing-flag parsing and shared default work state.

## 2026-08-04 — MATAFILER4 4.37 empty-state recovery

- Revalidate a cached `SMPL.empty` marker against the current map-resolved primary and supplementary input sizes before it can suppress requested work.
- Reopen a stored empty/too-small completion sentinel when its current input size no longer qualifies for that terminal outcome, even when the sentinel recorded the same nonempty size.
- Reject internally contradictory empty-sample classification instead of silently closing a valid sample.
- Accept `-profileMetaphlan3` as a compatibility alias for `-profileMetaphlan`.

## 2026-08-04 — IQ-TREE completion and strain-tip recovery

- Updated `buildTree5.pl` from 5.25 to 5.27.
- Retained IQ-TREE's standard identical-sequence handling without `-keep-ident`, while requiring exact alignment/tree taxon parity before accepting the inferred backbone as complete.
- Rejected partial tree files whose logs lack the IQ-TREE completion marker or end in an error. Numerical-underflow failures now restart once with the safe likelihood kernel; large strain alignments use that kernel pre-emptively.
- Removed completed-run IQ-TREE transients, including `.uniqueseq.phy`, checkpoints, distance matrices, and variable-site intermediates.
- Relaxed strict-backbone filtering from 70% to 35% of Q90 called-site coverage. A locus-QC placement flag alone no longer removes a well-covered sample after the affected loci have already been masked.
- Published the approximate nearest-backbone grafted tree as the default `IQtree_allsites.treefile` and retained the ML inference as `IQtree_allsites.backbone.treefile`; the grafts do not alter or improve the ML backbone topology.

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
