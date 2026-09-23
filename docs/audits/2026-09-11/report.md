# MATAF4 local algorithm and caller audit — 11 September 2026

The original audit found **11 reproducible defects**, plus a median defect corrected while reducing duplicated code. **All 11 have now been fixed in the follow-up implementation.** The descriptions below preserve the original failure evidence and audit locations; see [fixes and validation](fixes.md) for the resulting behavior and [results after the fixes](results.after-fixes.json) for rerun fixtures. [Original results](results.json) remain available for comparison.

P1 means a high-priority failure or silent corruption in an enabled workflow; P2 means a narrower input case or an incorrect derivative/statistic. These are code and fixture findings, not estimates of how much any existing biological dataset was affected.

## Original findings — all fixed

### 1. [P1] Secondary coverage checks use filenames the bundled calculator never writes

**Location:** `MATAF4.pl:5801`; also `1272` and `2599`. Caller: `scndMap2Genos`, around `10069`.

`calcCoverage2nd` receives `...bam.coverage.gz` and checks `$cov.pergene`, `$cov.percontig`, and `$cov.median.percontig`. The bundled `bin/rdCover` removes `.gz` before adding these suffixes. A real invocation on a 1,000-base contig with depth 10 produced `sample-smd.bam.coverage.pergene`, never `sample-smd.bam.coverage.gz.pergene`. All three names expected by the caller were absent, despite successful execution and correct mean depth.

The same incorrect names occur in initial completion checks and the sample-completion evidence. Enabling secondary per-gene coverage therefore keeps successfully processed samples incomplete and schedules repeated coverage work. The primary `separateContigs.pl::geneAbundance` path already uses the calculator's actual convention.

**Correction:** centralize the derivative stem calculation and reuse it in command construction, completion checks, and evidence records. Preserve compatibility with older outputs deliberately. The reproduction records the tested executable's SHA-256 because a configured replacement binary could differ.

### 2. [P1] mOTUs rank aggregation discards a taxon's abundance in later samples

**Location:** `secScripts/composition/mrgMotu2.pl:38–48`. Caller: `MATAF4.pl::mergeMotu2Table`, around `4331`.

`next if (exists $tax{$spl[0]})` skips both taxonomy parsing and the abundance additions for an mOTU encountered previously. Taxonomy needs parsing once; abundance must be accumulated for every sample.

**Reproduction:** the same mOTU has counts 12 in sample A and 3 in sample B. `m2.motu.txt` correctly contains `12,3`, while `m2.kingdom.txt` contains `12,0`. The omission affects every aggregated rank through species and depends on the alphabetical sample-processing order.

**Correction:** cache parsed lineages, then add each sample's abundance outside the cache-initialization branch. The existing compatibility test exercises only one sample, which explains why it does not detect this case.

### 3. [P1] Functional filtering and normalization reuse the first subject's length

**Location:** `secScripts/functions/parseBlastFunct2.pl:690–691`, used at `727`, `729`, and `905`. Callers: `MATAF4.pl::runDiamond`, around `5675`, and `Mods::FuncTools::assignFuncPerGene` through `secScripts/geneCat.pl`.

`$SbjLen` is set from the first candidate and is never updated while evaluating or selecting another subject. It controls the coverage threshold, the CARD length-scaled score threshold, and gene-length-normalized abundance.

**Reproduction:** candidate A is 1,000 aa; the stronger candidate B is 100 aa; alignment length is 50 aa. B is selected, but its GLN contribution is **0.05 instead of 0.5**. With a 0.4 coverage threshold, B is rejected even though it covers half its subject.

**Correction:** validate and use each candidate's own positive length during filtering, and the selected subject's length during normalization. Separately, `-minFractQueryCov` and `-minPercSbjCov` currently alias the same subject-length calculation: the former is not a query-coverage test despite its name. A true query-coverage option needs query-length data and translated-coordinate handling.

### 4. [P1] Functional profiling loses mate-2-only hits

**Location:** `secScripts/functions/parseBlastFunct2.pl:1224` and `1241–1244`. Caller: its read-group loop at `381`, invoked by `MATAF4.pl::runDiamond`.

`combineBlasts` returns the first mate unchanged when the second is empty, but has no symmetric return when the first is empty. The subsequent intersection-only loop drops every mate-2 hit. The singleton-preservation code below `next` is unreachable.

**Reproduction:** a single valid `read/1` hit produces B with count 1 and GLN 0.5; the otherwise identical `read/2` hit produces empty count and GLN tables. The script still exits successfully and writes its completion marker.

**Correction:** preserve both singleton orientations. Explicitly specify whether two mates hitting different subjects should compete independently or be rejected; that policy should not arise accidentally from dead code.

### 5. [P1] Valid empty ribosomal results cannot survive profile merging

**Location:** `secScripts/miTag/miTagTaxTable.pl:45–46` and `94`. Producer: `secScripts/miTag/lotus_LCA_blast3.pl:176–180`, with completion acceptance at `207–216`. Caller/publication: `MATAF4.pl::RiboMeta` and `riboSummary`.

When no marker reads exist, the assignment script explicitly creates an empty hierarchy and accepts it as complete. `RiboMeta` publishes it, but the merger requires a header and aborts. Merely adding a header does not repair the whole contract: sample columns come from `$sites`, which is populated only while reading data rows.

**Reproduction:** one populated hierarchy plus one empty hierarchy exits 255. Replacing the empty file with a header-only hierarchy exits 0, but the table contains only the populated sample. The zero-result sample disappears.

**Correction:** define one completed-empty hierarchy representation and preserve the complete input sample roster independently of observed taxa. SSU and LSU matrices should retain explicit zero columns for completed samples with no assignments.

### 6. [P2] Mapping statistics parse obsolete filter counters and report incorrect rates

**Location:** `MATAF4.pl:9557–9566` and `9585–9601`. Producer: `secScripts/assemblies/bamFilter.pl:191–202`.

The caller looks for `Inentries`, `TotalRetained`, and `TotalRm`. The filter now emits `Input records`, `Retained mapped records`, `Newly filtered records`, and separate already-unmapped/malformed counts. The caller consequently leaves retained/removed counts at zero and overwrites a parsed Bowtie2 total with -1. The Strobealign branch cannot populate its statistics because `incoming` stays zero.

**Reproduction:** actual filter output for two alignments, one passing and one failing mapping quality, combined with a two-sequence minimap2 log yields **AlignedReads=0 and OverallAlignment=100**, rather than 1 and 50.

**Correction:** put counter interpretation in a shared parser accepting current and legacy logs. Preserve the distinction between records, reads, and pairs, and include already-unmapped records in the appropriate denominator. Updating label regexes alone also leaves the Strobealign formula `$incoming / $locStats{totReadPairs}` equal to 1, rather than a percentage based on retained alignments.

### 7. [P2] Functional count weights infer a merged pair from the absence of a mate suffix

**Location:** `secScripts/functions/parseBlastFunct2.pl:822–826`. Callers: raw-read profiling and catalog annotation through `Mods::FuncTools`.

Every query without `/1` or `/2` is counted as two reads. Single-end reads and predicted gene identifiers commonly have neither suffix, and the caller does not pass their provenance into the parser. The same representation is also used for actual merged pairs.

**Reproduction:** changing the input identifier from `read/1` to the singleton identifier `read` changes count 1 to count 2 with identical alignment evidence. GLN remains 0.5. Mixed singleton/paired profiles therefore receive inconsistent count weights.

**Correction:** carry explicit read/gene/merged-pair provenance, and apply multiplicity from that provenance instead of guessing it from a missing suffix.

### 8. [P2] Paired-hit overlap calculation is wrong for containment and inclusive endpoints

**Location:** `secScripts/functions/parseBlastFunct2.pl:1271–1278`; duplicated calculation in `secScripts/functions/ABRblastFilter2.pl::combineBlasts`.

After ordering subject starts, overlap is computed as `end1 - start2`. It omits the inclusive endpoint and does not cap the overlap at the second hit's end. The result directly changes combined alignment length, bit score, threshold acceptance, and GLN.

**Reproduction:** ungapped subject intervals `[1,200]` and `[50,100]` have overlap 51 and union length 200. The function reports combined length **101**, also reducing the synthetic bit score to 101. Identical or one-base-overlapping intervals have the endpoint error as well.

**Correction:** share a validated inclusive-interval intersection calculation between both parsers. For gapped hits, document how reference-coordinate overlap is converted into alignment-column contributions; a bounding interval alone cannot reconstruct exact overlapping aligned columns.

### 9. [P2] Rolling gene k-mer windows repeat identifiers and assign different means to the same gene

**Location:** `secScripts/composition/kmer_Ngenes.pl:38–43`, `50–56`, and `63–68`. Caller: `separateContigs.pl` with window argument 5, enabled by `MATAF4.pl -kmerPerGene`.

The flush loop selects `$rSize - $numG` while simultaneously removing the first record. For several iterations this selects the same original gene. It also re-emits the last gene already emitted by the streaming branch. The same faulty flush logic is copied at contig transitions and EOF.

**Reproduction:** eight genes with a scalar k-mer feature 1 through 8 generate twelve rows; `ctg_4` appears five times with means 4.5, 4.5, 5, 5.5, and 6. A downstream keyed reader can overwrite earlier values, and the artifact does not contain one row per gene.

**Correction:** define the intended centered-window radius and track the next unreported gene independently of buffered records. Use one flush implementation for contig changes and EOF, reusing `avgArray`/`roundAr` or a shared rolling-sum helper.

### 10. [P2] The final k-mer record can crash on zero informative bases

**Location:** `secScripts/composition/calc.kmerfreq.pl:213–228`. Caller: `separateContigs.pl` for contig or gene k-mers.

The normal record-processing branch guards `totalkmers > 0`; its duplicated final-record branch divides without that guard. An eligible-length sequence containing only Ns, or with no complete valid k-mer, is skipped in one position but aborts in the last position.

**Reproduction:** a valid 120-base sequence followed by 120 Ns exits 255 with `Illegal division by zero` at line 228.

**Correction:** reuse a single record-emission function at header transitions and EOF, including an explicit policy for sequences without informative k-mers.

### 11. [P2] Undefined GC values concatenate adjacent output records

**Location:** `secScripts/composition/calcGC.pl:78–79`; analogous GC3 branches at `91–92` and `100–101`. Caller: `separateContigs.pl` with GC statistics enabled.

The zero-informative-base branch writes `ID\t-1` without a newline. The next sequence's result is appended to that same row, corrupting both identifiers and values. The contig GC3 empty branch also writes the gene tag instead of the contig key.

**Reproduction:** `ambiguous=NNNN`, followed by `valid=GCGC`, produces `ambiguous\t-1valid\t100.000` on one line, while returning success. Lowercase sequence also takes the uninformative branch because the counts only match uppercase bases.

**Correction:** emit every result through one newline-terminated row writer, keep contig and gene identifiers distinct, and normalize sequence case before base counting.

## Redundancy reduction completed

- `Mods::GenoMetaAss::median`, `Mods::MGSLocus::_median`, and `Mods::MosaicLoci::_median` now delegate to the existing `Mods::math::medianArray`, retaining each wrapper's empty-input result of zero. The former `GenoMetaAss` algorithm used `n-1` for parity/index selection: it returned 1 for `[1,3]` and 2 for `[1,3,9]`. They now return 2 and 3. This fixes marker-taxonomy median reporting, including the current `annotateMGwSpecIs3.pl` consumer, as part of the consolidation. Other quantile methods were left separate because they have different interpolation conventions.
- `Mods::Binning::MB2assigns` and `MB2assignedBinIds` share one parser and quality validation. `MB2assigns` accepts an omitted quality file for the quality-preparation stage; its existing two-argument return contract remains unchanged. `checkBinQual.pl` imports this helper instead of maintaining another parser. The ID-only path still stores only bin IDs, avoiding an unnecessary full contig-membership allocation.
- Focused tests cover odd/even/singleton/empty median inputs and agreement between pre-quality and quality-validated bin membership. Existing locus, mosaic, binning, and quality-processing tests also pass.

Other strain-controller and strain-documentation edits appeared concurrently in the shared checkout. They were not made or reverted by this audit.

## Evidence and scope

The original Python demonstrator (`reproduce.py`) has been retired; the project keeps its tests in Perl only. The fixtures it contained are covered by regression tests in [`t/local_algorithm_regressions.t`](../../../t/local_algorithm_regressions.t), which runs with the rest of the suite (`perl helpers/runTests.pl`). [Recorded results](results.json) and [results after the fixes](results.after-fixes.json) keep the numerical examples and the calculator hash from the original run.

Validation:

- Initial algorithm-focused baseline: **9 files, 440 checks passed**.
- Refactor-focused run: **7 files, 275 checks passed**.
- Full current-checkout regression suite: **59 files, 2,486 checks passed** (`prove -I. -It/lib t/*.t`). Existing tests passing does not negate the reproduced defects above.
- Refactor whitespace checks and main-script syntax validation were performed separately.

The pass traced active call boundaries for input libraries/cleaning, mapping filters, coverage and secondary mapping, contig statistics, bin assignment and quality checks, hybrid breakpoints/simulation, SNP region planning and consensus commands, functional profiling/catalog annotation, and mOTUs/ribosomal profile merging. Source inspection was most detailed in the algorithms and consumers cited above. Existing tests also exercised bedGraph coordinate conversion, breakpoint flanks/smoothing, coverage-weighted simulation, BAM CIGAR filtering, and consensus region coverage.

This is not a claim that every transitive script or branch is defect-free. The optional catalog/MGS/phylogeny workflow is much larger than the immediate MATAF4 call graph and was followed selectively for shared algorithm consumers. Dormant/deprecated branches and HPC scheduler behavior were not the focus. No real metagenomic dataset, aligner, assembler, or taxonomic database was run. Local C++ executables are present without their source; `rdCover` was checked with a uniform-depth fixture and `vcf2fna`'s supported interface was inspected, but their internal algorithms were not source-audited.
