# Mapper and search-tool options audit — 8 October 2026

## Scope and method

This audit checked every call of a read mapper or sequence search tool for options that are wrong, ignored, or never passed:

- bowtie2, bwa, minimap2, strobealign and kma (read mapping, indexes, decoy and competitive mapping, scaffolding);
- DIAMOND (read-based functional profiling, geneCat annotation, eggNOG-mapper);
- MMseqs2, LAMBDA, HMMER and vsearch;
- the profilers MetaPhlAn, mOTUs, Kraken and TaxaTarget.

RiboFind/miTag was audited separately today (`report.md`). Items in the earlier reports were left out.

Options were checked against the source code and manuals of the versions the installer pins (`helpers/install/*.yml`):

| Tool | Version |
|---|---|
| DIAMOND | 2.2.8 |
| bowtie2 | 2.5.4 |
| bwa | 0.7.19 |
| minimap2 | 2.31 |
| strobealign | 0.17.0 |
| MMseqs2 | 18-8cc5c |
| LAMBDA | 3.1.0 |
| eggNOG-mapper | 3.0.0-beta6 |
| MetaPhlAn | 4.2.6 |
| mOTUs | 4.1.0 |
| Kraken2 | 2.17.1 |
| samtools, bcftools | 1.23 |

None of the tools was run. The generated commands were run with stand-in tools under WSL (perl 5.38).

Verification:

- `perl -c` on every changed file.
- The full unit-test suite (`prove -I. t/`). Before the changes (bba777e): 87 files, 2,810 tests, one failing. That test still expected the bwa `.pac` check replaced in bba777e. After: 88 files, 2,837 tests, all passing.
- New `t/audit_2026_10_08_mappers.t` (44 tests). New cases in `t/alignment_input.t`, which runs the real `mapReadsToRef` commands with stand-in mappers.
- Against the code before these changes, the new tests fail at every fixed point.

## Fixed

| Location | Defect | Fix |
|---|---|---|
| `mapReadsToRef` | **Secondary mapping failed for every mapper except bowtie2.** bwa, minimap2, strobealign (and kma) were given `$REF`, which `scndMap2Genos` sets to all references joined by commas. More than one reference (`-competitive2ndmap 0`) or competitive modes 1/2 failed. With the default `-mapper -1`, long-read samples use minimap2, so long-read map2tar never worked. Staging a `.gz` reference also treated the whole list as one file. | Each reference is mapped against its own FASTA. Competitive modes use the combined DB and decoy mode the decoy FASTA. `.gz` references are staged one by one; two with the same file name get separate copies. |
| `mapReadsToRef`, `deployMapDB.pl` | **Decoy mapping silently skipped** for the same mappers: they mapped to the reference, not the decoy DB. `deployMapDB.pl` always built a bowtie2 index. | The mapper uses the decoy FASTA. `deployMapDB.pl` takes the mapper (new optional 7th argument) and builds a bwa index when needed. |
| `mapReadsToRef` | **`-mapper 2` (bwa) killed the controller** for any sample with singletons ("single end mapping not implemented for bwa"). | bwa maps single-end files, and alignment input (BAM/CRAM) from stdin like the other mappers. |
| `mapReadsToRef` | The bowtie2 index check accepted only the suffix implied by `-mapperLargeRef`. bowtie2-build writes `.bt2l` by itself above ~4 Gbp, and the decoy DB was always built as `.bt2`. Large references with `-mapperLargeRef 0`, or decoy mode with `-mapperLargeRef 1`, exited 23 on every pass. Only the first of several indexes was checked. | `.bt2` or `.bt2l`, for every reference. |
| `getAlgnCmdBase` | PacBio reads were mapped with minimap2 `-x map-pb` (CLR preset). Everywhere else MATAFILER treats `PB` as HiFi (metaMDBG `--in-hifi`, bcftools `-X pacbio-ccs`, the PB bamFilter cutoffs). | `-x map-hifi`. |
| kma (`-mapper 4`) | kma's SAM output has no `NM:i` tag, so `bamFilter.pl` dropped every record: an empty BAM without an error. Its presets also overrode `-mq/-bcd/-mrc`, and `-ID 0.95` meant 0.95 %. kma was not installed. | **kma removed** (code, config, docs, installer). `-mapper 4` stops MATAF4 at startup with a message. `cleanup_finished_sample.pl` still removes old `.kma*` index files. |
| `buildMapperIdx` (also `geneCat -ntMatchGC`) | minimap2 got a prebuilt `.mmi` (`-H`, default k/w). minimap2 takes `-H/-k/-w` from the index and ignores the preset of the mapping call, so `-ntMatchGC` lost its `-x asm20` seeding. MATAF4 never read the `.mmi`, yet `-mapper 3` submitted an index job for it. | minimap2 maps the FASTA, like strobealign. (Committed earlier today in bba777e.) |
| `buildMapperIdx`, `mapperDBbuilt` | A bwa index counted as complete once `.pac` existed. bwa writes `.pac` first and `.sa` last, so an interrupted index job was accepted. | Checked by `.sa`. (bba777e) |
| `scaffoldCtgs` | Mate libraries are always mapped with bowtie2, but the index was built for `-mapper`. With bwa/minimap2/strobealign, BESST scaffolding failed. `samtools sort -T kk` wrote into the job's working directory. | bowtie2 index; sort temp files in the scaffolding work directory. |
| `runDiamond` | **A read pair counted 1 instead of 2 when the read names had no `/1` `/2` suffix.** R1 and R2 were searched separately; the parser recognises mates only by that suffix. sdm adds it only to Illumina 1.8 headers, so SRA/ENA downloads (`fasterq-dump`) were counted at up to half weight (`cnt` and `GLN`) next to other samples. The `sort -k1` that brought mates together also reordered each read's hits by subject name, and the best-hit choice depends on order for near-ties. | **Both mates of a pair go through one DIAMOND run.** DIAMOND has no paired mode, so the new `secScripts/functions/interleaveMates.pl` interleaves R1/R2 on stdin, names the mates `<read>/1` and `<read>/2`, and stops on mates out of order or files of different length. DIAMOND writes hits in query order, so each pair's hits stay together in DIAMOND's score order, and no sort is needed. |
| `runDiamond` | `--min-orf 25` also applied in frameshift mode (`-DiaFrameshift`), where DIAMOND itself uses no ORF filter for error-prone reads. | No `--min-orf` with `-F`. |
| `runDiamond` | The search e-value was fixed at 1e-4, so any `-DiaParseEvals` value above it was silently capped. | The search uses the most permissive `-DiaParseEvals` value (at least 1e-4). |
| `metphlanMapping` | **`--nreads` counted pairs.** Mates were mapped as pairs, and the read total came from bowtie2's summary, which counts pairs. MetaPhlAn counts every SAM record and scales all abundances by mapped/total reads. For paired samples, the mapped fraction was up to 2× too high, UNCLASSIFIED too low and all taxa inflated. | All mates are mapped unpaired (`-U R1,R2,S`), as MetaPhlAn does itself. |
| `metphlanMapping` | The read total was read from the scheduler's `.etxt`, which does not exist with `-submSystem bash` (the job failed) and is fragile under LSF. | bowtie2's stderr goes to a log in the work directory, which is parsed. |
| `metphlanMapping` | The job requested 3 GB in total. bowtie2 loads the whole MetaPhlAn index (about 33 GB for vJan25), so the job is OOM-killed wherever memory is enforced. | Index size + 6 GB (40G if the index files cannot be seen at submission). |
| `prepMetaphlan`, `taxPerMGS_gtdb.pl` | The version probes ran the `env:` prologue (`[[ … ]]`) through /bin/sh. On dash hosts the environment was not activated and the probe failed. | Run through bash. |
| geneCat marker clustering | mmseqs ran with `--split-memory-limit 42G` and the main job's thread count inside `cogCluster.sh`, which requests 25G and `cores/2`. | Threads and memory of that job (geneCat 0.65). |
| `FuncTools::lambdaBl` (`MG_LCA.pl`) | LAMBDA 3.1 `mkindexn` accepts `-t` 2–1000, so a 1-core call failed. | At least 2. |
| `ABRblastFilter2.pl` | Any query ID ending in "2" counted as mate 2 (e.g. `SRR1.12`). | Only IDs ending in `/2`. |

### What changes for existing projects

- **Read-based DIAMOND.**
  - Existing searches (`dia.<DB>.blast.srt.gz`) are kept; only new searches use the joint mate run, the new e-value and the new ORF setting.
  - Samples whose read names had no `/1` `/2` suffix (SRA/ENA downloads) were counted with each pair as one read. To correct them, rerun with `-reProfileFunct 1`.
- **MetaPhlAn.** Existing profiles are kept. To recompute one with the corrected `--nreads`, delete its `<sample>.MP2.sto`.
- **Mapping.**
  - Existing BAMs and coverage are kept.
  - New PacBio mappings use `map-hifi`.
  - map2tar with bwa, minimap2 or strobealign and several references, decoy mapping with these mappers, and bwa on samples with singletons now run.
- **`-mapper 4`** stops at startup.
- **Unused `.mmi` files** next to assemblies are no longer built. Existing ones are ignored, and the finished-sample cleanup removes them.

## Found, not fixed

These need a decision or testing on a cluster.

- **DIAMOND read jobs request 7 GB.** `-DiaMem` defaults to 7, and no `-b/-c/--memory-limit` is passed. DIAMOND 2.2.8 blastx defaults to `-b 2.0`; the wiki says memory is "roughly six times" that number in GB, and the align stage budgets `--memory-limit` 16G. Suggest 16G, or pass `--memory-limit`.
- **Scratch for read DIAMOND is 80G for every database.** For KGM/NOG, geneCat needed 500G with 100M-letter chunks; a read query block is up to 2 G letters. Suggest the `bigFuncDB` sizing.
- **Long reads in `-profileFunct`.** They get the short-read search (`-k 5`, no `--range-culling`/`--long-reads`), and the parser keeps one function per read, so a 10 kb read counts one gene.
- **Changing DIAMOND search or parse options never invalidates results.** This covers `-DiaParseEvals`, `-DiaPercID`, `-DiaMinAlignLen`, `-DiaMinFracQueryCov`, `-DiaSensitiveMode` and `-DiaFrameshift`. The stone records no parameters; geneCat got a `.params` record on 10-01.
- **`--min-orf 25` for short reads.** DIAMOND's own default is no filter for frames under 30 aa and 20 aa below 100 aa; 25 masks reads under 75 nt completely.
- **geneCat FuncAssign sensitivity.** It runs DIAMOND in default mode (designed for >60 % identity) while the parser accepts 25 %; the eggNOG route runs `--sensitive`.
- **KGM/NOG chunk memory.** These jobs request 160G, but DIAMOND runs with its defaults (~12–20 GB).
- **`-competitive2ndmap 2` adds bowtie2 `-a`.** Secondary alignments get MAPQ 255, the tied primary gets MAPQ 0/1 and is removed by bamFilter, and `samtools depth` skips secondaries. Coverage therefore equals mode 1, at extra cost.
- **Decoy with competitive modes 1/2.** No decoy DB is built, although `docs/common_workflows.md` shows `-competitive2ndmap 1 -decoyMapping 1`.
- **`-profileKraken` is hard-wired to Kraken 1** (`kraken --preload --fastq-input`, `kraken-filter`, `kraken-translate`). MF4 installs Kraken2 and a `.k2d` database, so every job fails on a standard install.
- **TaxaTarget.** "No reads mapped" exits 1, so samples without protist hits are resubmitted on every pass.
- **Mosaic loci.** `-x asm20` implies `-s200`. At 90 % identity that needs about 400 aligned bp, and 80 % never passes, so the 0.80–0.95 outgroup window is largely unreachable. Consider `-s 40` after `-x`.
- **bcftools (`Mods/SNP.pm`).** `-X ont/pacbio-ccs` comes after `--min-BQ 30` and resets it to 5 for long reads.
- **mmseqs linclust** (catalogue) matches both strands; CD-HIT is sense-only.
- **FOAM hmmsearch** has no `-Z`, so `-E 1e-5` is applied per chunk.
- **mOTUs and TaxaTarget** request 3G; measure before relying on it.
- **Lower priority:**
  - `getMapStats` parses only the first bowtie2 summary.
  - minimap2 above 8 Gbases builds a multi-part index.
  - `-mapUnmapped` dies in `seedUnzip2tmp`.
  - The read group ID is the sample name for every library.
  - The mapper and samtools both run `N` threads in an `N`-core job.
  - The CD-HIT `-M` exceeds its job.
  - `decluterGC.pl` always rebuilds its mmseqs DB.
  - The MMseqs2 branch of `runDiamond` was unreachable and broken (`fident` is 0–1; `--compressed` does not gzip `.m8`). It was removed with the mate change.
