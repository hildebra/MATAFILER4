# Downloading and analysing public metagenomes

MATAFILER4 downloads public reads from ENA or SRA, validates them, analyses them like local samples and deletes the downloaded copy once the sample is finished. You list one sample per row in the mapping file with its run accession(s) in an `ENAdownload` or `SRAdownload` column; MATAF4 does the rest.

One step stays outside MATAF4: turning a study accession into one row per sample and dropping 16S amplicon runs. A study accession in a single row would download every run of the study into one sample, amplicons included. Step 2 below covers this.

```text
Inspect the study on ENA (Step 1)
  -> one map row per sample (Step 2)
  -> MATAF4 submits a download job per sample
  -> download job: fetch, check bytes, MD5, FASTQ records
       all files valid?  no  -> the job stops with an error; the next MATAF4 pass resubmits it
                         yes -> the next MATAF4 pass analyses the sample
  -> cleanup deletes the downloaded reads only
```

## Validation on PRJNA776530 (October 2026)

PRJNA776530 is a human gut microbiome diet-intervention study ([ENA record](https://www.ebi.ac.uk/ena/browser/view/PRJNA776530)). 21 healthy adults ate the same 7-day homogeneous average American diet, which tested whether a standardised diet reduces between-person variation in microbiome-dependent metabolites (phenylacetylglutamine, p-cresol sulfate and indoxyl sulfate as primary outcome). The record holds 105 BioSamples from these participants.

The study holds 210 sequencing runs from 105 BioSamples. Each BioSample was sequenced twice: once as 16S rRNA amplicon (Illumina MiSeq, `library_selection` PCR, about 2 x 300 bp) and once as shotgun metagenome (Illumina NovaSeq 6000, `library_selection` RANDOM, about 2 x 125 to 150 bp). Half the runs are therefore 16S and only the 105 shotgun runs, 85.8 GB, are wanted.

A map with the 105 shotgun runs produced 105 download jobs, all finished in one pass without intervention. An independent `md5sum` check against ENA confirmed 210 of 210 FASTQ files byte-identical, with no amplicon file present. Deleting one downloaded file and running MATAF4 again re-submitted the download for that sample only, and the new file matched ENA. Four further samples were run through download, analysis and cleanup, three of them with `-keepRawReads 1`.

## One-time setup

Install MATAFILER4 with its installer, then check these entries.

| Item | Where | What to set |
| --- | --- | --- |
| `globalTmpDir` | `Mods/MATAFILERcfg.txt` | any path you choose, usually your scratch; downloads go to `<globalTmpDir>/<RunID>/<SmplID>/accession_download/` |
| `#OutPath` | the mapping file | any path you choose; results and `input_accession.json` go to `<OutPath>/<RunID>/<SmplID>/` |
| `downloadQueue` | `Mods/MATAFILERcfg.txt` | a partition whose compute nodes can reach ENA and NCBI over HTTPS |
| `downloadMem` | `Mods/MATAFILERcfg.txt`, optional | leave unset (default 4G); a download job used about 400 MB in tests |
| `wget` | MF4 environment | needed for `ENAdownload` |
| `prefetch`, `vdb-validate`, `fasterq-dump` | MF4 environment | needed for `SRAdownload` only: `micromamba install -n MF4 -c bioconda -c conda-forge sra-tools` |

MATAF4 runs on the same nodes as usual; it sends the download jobs to `downloadQueue` itself. Only Step 1 and `makeAccessionMap.sh` with a study accession contact ENA directly, so run those where you have internet access.

If you keep a second MATAFILER4 install beside an older one, point `PERL5LIB` at the install you run (`export PERL5LIB=$MF4DIR`). Several helper scripts load `Mods::` through `PERL5LIB`, and an old path there loads the old modules without an error message.

## Step 1: inspect the study on ENA

Check what the study contains before downloading anything; this takes seconds and downloads nothing. It queries ENA, so run it where you have internet access:

```bash
ACC=PRJNA776530
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$ACC&result=read_run&fields=run_accession,sample_accession,library_strategy,library_selection,instrument_model,read_count,base_count&format=tsv" > $ACC.ena.tsv
tail -n+2 $ACC.ena.tsv | wc -l          # number of runs
tail -n+2 $ACC.ena.tsv | awk -F'\t' '{l=($6>0)?int($7/$6):0; print $3"\t"$4"\t"$5"\tavgLen~"int(l/50)*50}' | sort | uniq -c
```

The second command groups runs by strategy, selection, instrument and bases per spot. For PRJNA776530 it prints:

```
    105 AMPLICON  PCR     Illumina MiSeq         avgLen~600
    105 WGS       RANDOM  Illumina NovaSeq 6000  avgLen~250
```

What to look for:

- **Mixed studies are common.** Four of 13 studies in one diet meta-analysis mix 16S and shotgun under one accession.
- **`library_strategy` is typed by the submitter and is often wrong.** Shotgun runs appear as `OTHER` (PRJEB26842) or `WGA` (SRP045211), so never keep only `WGS`; exclude `AMPLICON` instead.
- **A group labelled `WGS` or `OTHER` but with `library_selection` `PCR`, or MiSeq reads near 2 x 300 bp, is probably amplicon.** Decide these by hand.
- **Runs with 0 reads exist.** ENA keeps adding empty runs to old studies (PRJNA703757 gained one in November 2025).

## Step 2: build the mapping file

Write one map row per sample, with that sample's run accession(s) in `ENAdownload`. `secScripts/fileManage/makeAccessionMap.sh` does this for a whole study. Given a study accession it queries ENA, so it also needs internet access; given the TSV from Step 1 it works offline:

```bash
bash $MF4DIR/secScripts/fileManage/makeAccessionMap.sh $ACC $ACC /path/to/matafiler_output/ > $ACC.map
grep -vc '^#' $ACC.map        # number of samples
```

What it writes:

```text
#SmplID	ENAdownload	AssmblGrps
#RunID	PRJNAxxxxxx
#OutPath	/path/to/matafiler_output/
SAMNaaaaaaa	SRRaaaaaaa	SAMNaaaaaaa
SAMNbbbbbbb	SRRbbbbbb1,SRRbbbbbb2	SAMNbbbbbbb
```

- **One row per BioSample.** Several runs of one sample go comma separated into one cell and are downloaded into the same sample.
- **Excluded runs:** `library_strategy` `AMPLICON`, `library_selection` `PCR`, and runs with 0 reads. Each excluded run is printed on stderr with the reason, followed by a summary line.
- **`SmplID` is the BioSample accession.** Replace the first and third columns if you want other sample names or co-assembly groups (`AssmblGrps`).
- **Optional columns** such as `SeqTech` can be added; see `docs/mapping_files.md`.
- The script cannot recognise an amplicon run labelled neither `AMPLICON` nor `PCR`. Step 1 is what catches those.

The accessions above are placeholders. The second row shows a sample sequenced in two runs.

## Step 3: run MATAF4

Run MATAF4 as you normally do, on the usual nodes, inside `tmux` or `screen`:

```bash
export MF4DIR=/path/to/MATAFILER4
export PERL5LIB=$MF4DIR
micromamba activate MF4
perl $MF4DIR/MATAF4.pl -map $ACC.map -submit 1 -maxConcurrentJobs 400 <your analysis options>
```

What each pass does:

1. **First pass:** one download job per sample goes to `downloadQueue`, and every downloaded sample is deferred (`sample processing is deferred until validation completes`). No analysis job is submitted yet.
2. **Download jobs:** each fetches its files, checks ENA byte counts and MD5, decompresses every FASTQ to check its records, then writes `input_accession.json`. Any failure ends the job with an error. The download partition may cap how many of your jobs run at once; the rest wait in the queue.
3. **Next passes:** run the same command again once the download jobs have finished. Samples whose download is complete start their analysis like local samples. Samples whose download failed or whose files are missing get a new download job.
4. **Repeat** after each round of jobs until no new job is submitted.

`-maxConcurrentJobs` is raised because MATAF4 otherwise pauses at "running + pending below 100" when you already have other jobs queued.

`-loopTillComplete` should repeat these passes for you. It has not yet been tested together with downloads; until it has, run the passes by hand.

## Step 4: check the downloads

Every downloaded sample has `input_accession.json` in its output folder, so a count against the map shows whether all downloads finished:

```bash
OUT=/path/to/matafiler_output/$ACC
ls $OUT/*/input_accession.json | wc -l          # equals the number of map rows
sacct -u $USER -S <date> -X --format=JobName%20,State,ExitCode -n | grep _DL | awk '{print $2, $3}' | sort | uniq -c
```

| Where | What is there |
| --- | --- |
| `<globalTmpDir>/<RunID>/<SmplID>/accession_download/` | the downloaded FASTQs, in scratch |
| `<OutPath>/<RunID>/<SmplID>/input_accession.json` | accession, file names, bytes, MD5, record count, download time; stays after scratch is cleaned |
| `<OutPath>/<RunID>/accession_manifest.tsv` | one row per downloaded file of the run, with its current path and whether it is still present |
| `<OutPath>/<RunID>/<SmplID>/LOGandSUB/accession_download.sh.otxt` / `.etxt` | the download job's log |

The download job has already checked every file against ENA. For an extra check that does not rely on MATAFILER, compare `md5sum` of each file against the `fastq_md5` column of the ENA filereport; this is how PRJNA776530 was validated.

## Step 5: after the analysis

When all required steps of a sample have finished, MATAF4's cleanup removes that sample's `accession_download/` folder. This happens by default (`-reduceScratchUse 1`) and frees the scratch the download used.

- **Only downloaded reads are deleted.** Cleanup is given the download folder only for samples with `ENAdownload` or `SRAdownload`, and refuses any path outside that sample's scratch folder. Local FASTQs named in `Path` or `SmplPrefix` are never touched.
- **Downloaded reads are kept until every step that reads them is done**, including a Protal job if `-profileProtal` is on.
- **Expect one extra download per sample at the end.** The cleanup job deletes the reads before the sample is marked complete. On the next pass MATAF4 finds the reads missing, downloads them again, and only on the pass after that marks the sample complete and deletes them once more. The cost is one more download of each sample; the results are not affected.
- **A rerun after completion does not download again.** Once a sample is marked complete (`MF4.sentinel.<SmplID>.json` in its output folder), MATAF4 skips it.

### Keeping the downloaded reads: `-keepRawReads 1`

With `-keepRawReads 1`, downloads go to `<globalTmpDir>/<RunID>/accession_reads/<SmplID>/` instead of the sample's scratch folder. Cleanup still removes the sample's scratch folder but never the kept reads, and later passes reuse them instead of downloading again. Use it when your quota allows and you expect to re-run with other parameters. Scratch purge policies still apply to kept reads.

If a later download of the same accession differs from the recorded `input_accession.json` (ENA MD5, or SRA record count), the download job stops and names the changed files, so results stay reproducible. Delete that `input_accession.json` to accept the new archive content.

## Known limitations

| Limitation | Effect | What to do |
| --- | --- | --- |
| A study accession in one map row | every run of the study, amplicons included, lands in one sample | one row per sample (Step 2) |
| `--exclude-strategy` cannot be set from the map | `ENASRAdl.pl` filters only when called by hand | filter while building the map (Step 2) |
| Cleanup deletes the reads one pass before the sample is marked complete | each sample is downloaded a second time before it closes | accept the extra download, or use `-keepRawReads 1` if scratch allows |
| Retries inside a download job are 3 immediate attempts | an ENA outage of minutes fails the job | run MATAF4 again later; only failed samples are resubmitted |
| Readiness check is "file exists and is not empty" | a file damaged on scratch after download, same size, is not noticed | rare; delete the file to force a new download |
| SRA downloads record no layout or platform | SeqTech defaults to Illumina | set `SeqTech` in the map for non-Illumina SRA data |
| Cleanup treats an empty `FMGids.txt` as unfinished | small samples with no marker gene keep their downloaded reads in scratch | delete by hand once results are checked |

Not yet tested: a real download failure followed by a rerun, and `-loopTillComplete` together with downloads.
