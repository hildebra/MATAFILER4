<!-- Documentation navigation -->
[Home](../README.md) | [Quick start](quickstart.md) | [Installation](install.md) | [Configuration](configuration.md) | [Mapping files](mapping_files.md) | [Workflows](common_workflows.md) | [Examples](examples.md) | [Profiling tutorial](profiling_tutorial.md) | [Strain-within guide](strainwithin.md) | [Outputs](outputs.md) | [Flag reference](flag_reference.md) | [FAQ](FAQ.md) | [Glossary](glossary.md)

---

# Read-based profiling tutorial

This tutorial covers RiboFind, DIAMOND functional profiling, MetaPhlAn 4,
mOTUs 4, and both Protal modes in `MATAF4.pl`. These workflows can be run with
`-assembleMG 0`; they do not require a metagenome assembly.

## What each profiler reads

| Function | Main flag | Input used by MATAFILER | Job layout |
|---|---|---|---|
| RiboFind (SSU/LSU) | `-profileRibosome 1` | Primary SDM-cleaned short reads | Per sample, followed by run-level SSU/LSU merges |
| Functional profiles | `-profileFunct 1 -DiaDBs ...` | Primary SDM-cleaned reads | Per sample and database, followed by run-level matrices |
| MetaPhlAn | `-profileMetaphlan 1` | Primary SDM-cleaned reads | Per sample, followed by run-level matrices |
| mOTUs | `-profileMOTU2 1` | Primary SDM-cleaned reads | Per sample, followed by run-level matrices |
| Protal singular | `-profileProtal 1` | One compatible primary raw short-read pair per sample | Per sample, followed by one profile merge |
| Protal combined | `-profileProtal 2` | One compatible primary raw short-read pair per sample | One generated map and one Protal job for the complete cohort |

Protal is deliberately the exception: it profiles staged raw reads. The other
functions use the cleaned read set produced by the normal MATAFILER read stage.

## 1. Install or update the profilers

Run the installer before configuring the database paths:

```bash
bash helpers/install/installer.sh
```

The versions checked against the upstream interfaces on 24 August 2026 are:

| Program | Bundled version | Conda environment | Database used by this tutorial |
|---|---:|---|---|
| MetaPhlAn | 4.2.6 | `MF4checkm2` | `mpa_vJan25_CHOCOPhlAnSGB_202503` |
| mOTUs | 4.1.0 | `MF4motus` | mOTUs marker-gene DB 4.1 |
| Protal | 0.6.0a | `MF4` | A compatible Protal full or mini database |

MetaPhlAn 4.2.6 and mOTUs 4.1.0 were the current upstream releases at that
audit date. See the [MetaPhlAn releases](https://github.com/biobakery/MetaPhlAn/releases)
and [mOTUs releases](https://github.com/motu-tool/mOTUs/releases) when updating
the pins later.

The installer downloads the MetaPhlAn and mOTUs databases on a best-effort
basis. A network or permission failure leaves the version marker absent, so a
later installer run retries it. To force both downloads again:

```bash
bash helpers/install/installer.sh --refresh-databases
```

Equivalent manual database commands are:

```bash
micromamba run -n MF4checkm2 metaphlan --install \
  --db_dir /path/to/DBs/MP4 \
  --index mpa_vJan25_CHOCOPhlAnSGB_202503

micromamba run -n MF4motus motus downloadMGDB \
  -db /path/to/DBs/mOTUs
```

For Protal, download a database using the
[official Protal instructions](https://protal.earlham.ac.uk/main.php?site=documentation#download-the-database).
The MATAFILER installer installs the program, not that external database.

Verify the programs after installation:

```bash
micromamba run -n MF4checkm2 metaphlan --version
micromamba run -n MF4motus motus --help
micromamba run -n MF4 protal --help
perl MATAF4.pl -checkInstall
```

## 2. Link the databases in `config_DB`

The bundled database-path file is `Mods/config_DBs.txt`. A site or user config
may override the same keys. All of these paths are resolved through
`getProgPaths()`, so keep the established key names even where they contain an
older tool-version suffix.

### Taxonomic profilers

```text
# RiboFind: reference FASTA and matching LCA taxonomy file
LSUdbFA      [DBDir]/MarkerG/SILVA/138.1/SLV_138.1_LSU.fasta
LSUtax       [DBDir]/MarkerG/SILVA/138.1/SLV_138.1_LSU.tax
SSUdbFA      [DBDir]/MarkerG/KSGP/v4.0/KSGPv4.0.fasta
SSUtax       [DBDir]/MarkerG/KSGP/v4.0/KSGP_plus2.tax

# MetaPhlAn: full database prefix, without .pkl or Bowtie2 suffixes
metPhl2_db   [DBDir]/MP4/mpa_vJan25_CHOCOPhlAnSGB_202503

# mOTUs: parent directory which contains db_mOTU/
motus2_DB    [DBDir]/mOTUs

# Protal: database directory, or leave empty and export PROTAL_DB_PATH
protal_db    /path/to/protal_database
```

Only `LSUdbFA`, `LSUtax`, `SSUdbFA`, and `SSUtax` are required by the current
RiboFind implementation. The `SSUdbFAsrt`, `LSUdbFAsrt`, `SSUidx`, and `LSUidx`
entries are retained for the older SortMeRNA path but are not read by the
current `detectRibo()` workflow. If a LAMBDA nucleotide index is missing,
MATAFILER builds it beside the configured FASTA before copying it to scratch;
the reference directory therefore needs to be writable for a first run.

`metPhl2_db` is a compatibility key name. For MetaPhlAn 4.2 it must point to
the index prefix. MATAFILER splits that value into `--db_dir` and `--index`.
Likewise, `motus2_DB` is retained for compatibility but now points to the mOTUs
4 parent directory rather than the `db_mOTU` child itself.

If `protal_db` is empty, export the database path in both the controller and
compute-job environment:

```bash
export PROTAL_DB_PATH=/path/to/protal_database
```

### Functional databases

Select one to six aliases with `-DiaDBs`. Each alias resolves to the following
`config_DB` directory key and required protein FASTA:

| `-DiaDBs` alias | `config_DB` key | Required reference file in that directory |
|---|---|---|
| `NOG` | `eggNOG40_path_DB` | `eggnog4.proteins.all.fa` |
| `MOH` | `Moh_path_DB` | `Extra_functions.faa` |
| `MOH2` | `Moh_path_DB` | `Nitrogen_cycl_genes.faa` |
| `CZy` | `CAZy_path_DB` | `Cazys_2019.fasta` |
| `ABR` | `ABRfors_path_DB` | `ardb_and_reforghits.fa` |
| `ABRc` | `ABRcard_path_DB` | `card.parsed.f11.faa` |
| `KGE` | `KEGG_path_DB` | `genus_eukaryotes.pep` |
| `KGB` | `KEGG_path_DB` | `species_prokaryotes.pep` |
| `KGM` | `KEGG_path_DB` | `euk_pro.pep` |
| `ACL` | `ACL_path_DB` | `aclame_proteins_all_0.4.fasta` |
| `TCDB` | `TCDB_path_DB` | `tcdb.faa` |
| `PTV` | `PATRIC_VIR_path_DB` | `PATRIC_VF.faa` |
| `PAB` | `ABprod_path_DB` | `dedup_best_prod_predictions.faa` |
| `VDB` | `VirDB_path_DB` | `VFDB_setB_pro.fas` |
| `URE` | `URE_path_DB` | `ualpha_gtdb_proteins.faa` |
| `URacc` | `URE_path_DB` | `urease_accessory_gtdb_proteins.faa` |
| `AMI` | `URE_path_DB` | `amidohydrolase_gtdb_proteins.faa` |

For the commonly used `KGM,NOG,CZy` set, a minimal path block is:

```text
KEGG_path_DB       [DBDir]/Funct/KEGG/
eggNOG40_path_DB   [DBDir]/Funct/eggNOG10/
CAZy_path_DB       [DBDir]/Funct/CAZy/
```

Some parsers also require annotation sidecars in the same directory:

- `NOG`: `all_species_data.txt`, `NOG.members.tsv`, and `NOG.annotations.tsv`.
- `CZy`: `MohCzy.tax` and `cazy_substrate_info.txt`.
- `KGE`, `KGB`, or `KGM`: `genes_ko.list` and `kegg.tax.list`.
- `ABRc`: the matching `card*.txt` and `card*.map` files.
- `PTV`: `PATRIC_VF2.tab`.
- `VDB`: `VF.tab`.
- `TCDB`: `TCDBhir.txt`.
- `PAB` with its taxonomy check enabled: `all_species_data.txt` from the
  configured `NOG` directory.

MATAFILER creates the DIAMOND `*.db.dmnd` and `*.length` files beside the
reference when absent, then copies the needed files to scratch. The database
directory must therefore be writable for initial index creation and must have
enough free space for the FASTA plus its index.

## 3. Start with a dry run

The following checks one sample without submitting jobs:

```bash
MAP=/path/to/cohort.map

perl "$MF4DIR/MATAF4.pl" \
  -map "$MAP" \
  -assembleMG 0 \
  -requireInput 1 \
  -profileRibosome 1 \
  -profileFunct 1 -DiaDBs KGM,NOG,CZy -DiaCores 8 \
  -profileMetaphlan 1 \
  -profileMOTU2 1 \
  -profileProtal 1 \
  -submit 0 \
  -from 0 -to 1
```

Inspect the printed database paths, program versions, generated commands, and
sample read discovery. Then submit the same range with `-submit 1`.

## 4. Run individual functions

### RiboFind

```bash
perl "$MF4DIR/MATAF4.pl" -map "$MAP" -assembleMG 0 \
  -profileRibosome 1 -submit 1
```

Use `-reProfileRibosome 1` to redo extraction and assignment, or
`-reRibosomeLCA 1` to retain extraction and redo only the LCA assignment.

### Functional profiling

```bash
perl "$MF4DIR/MATAF4.pl" -map "$MAP" -assembleMG 0 \
  -profileFunct 1 -DiaDBs KGM,NOG,CZy -DiaCores 8 -submit 1
```

`-reProfileFunct 1` repeats the read-to-database search. `-reParseFunct 1`
reparses retained search results where available.

### MetaPhlAn 4

```bash
perl "$MF4DIR/MATAF4.pl" -map "$MAP" -assembleMG 0 \
  -profileMetaphlan 1 -submit 1
```

The current interface maps cleaned reads with Bowtie2 and passes the resulting
SAM to MetaPhlAn 4.2 using `--input_type sam`, `--nreads`, `--db_dir`, and
`--index`. Unclassified-abundance estimation is already the MetaPhlAn 4.2
default. The removed `--add_viruses` and `--unclassified_estimation` options are
not sent to 4.2. MATAFILER retains its older branch for MetaPhlAn 3 and 4.0/4.1,
but the installer and this tutorial target 4.2.6.

### mOTUs 4

```bash
perl "$MF4DIR/MATAF4.pl" -map "$MAP" -assembleMG 0 \
  -profileMOTU2 1 -submit 1
```

MATAFILER passes multiple read files as separate arguments to mOTUs 4 and uses
its required `-o` output option. The temporary BAM, MGC, inserts, and relative-
abundance files remain in job scratch; only the compressed count profile is
published. Paired files must still contain corresponding reads in the same
order, as required by mOTUs.

### Protal per sample

```bash
perl "$MF4DIR/MATAF4.pl" -map "$MAP" -assembleMG 0 \
  -profileProtal 1 -ProtalCores 4 -ProtalMem 100 -submit 1
```

`-protalIgnoreErrors 1` is the default. It chooses the first compatible primary
raw short-read pair, ignores additional pairs and singleton/BAM streams with a
warning, and skips samples with no compatible pair. Set it to `0` for strict
input validation. Mode 1 disables strain analysis and retains only the durable
per-sample profiles required for the final merge.

### Protal combined cohort

```bash
perl "$MF4DIR/MATAF4.pl" -map "$MAP" -assembleMG 0 \
  -profileProtal 2 -ProtalCores 4 -ProtalMem 100 -submit 1
```

Mode 2 requires the complete mapped cohort, including all supplied mapping
files, so omit `-from` and `-to`. MATAFILER generates one Protal map, holds the
staged raw directories until the combined job finishes, retains strain MSAs in
the final directory, and then runs its validated scratch-cleanup script.

## 5. Complete the run and merge profiles

Wait for the submitted jobs, then rerun the identical controller command.
MATAFILER reuses completion markers and submits or runs the run-level mergers
only when the required sample profiles are complete. This rerun is important;
do not treat the end of the first controller invocation as the end of all batch
work.

Typical durable outputs are:

| Function | Per-sample or source output | Run-level output |
|---|---|---|
| RiboFind | `<sample>/ribos/` | `pseudoGC/Phylo/RiboFind/SSU.miTag.<rank>.txt` and `LSU.miTag.<rank>.txt` |
| Functional | `<sample>/diamond/` | `pseudoGC/FUNCT/<DB_alias>/` |
| MetaPhlAn | `pseudoGC/Phylo/MP2/<sample>.MP2.txt` | `pseudoGC/Phylo/MePh.all.<rank>.mat` |
| mOTUs | `pseudoGC/Phylo/mOTU2/<sample>.motu2.tab.gz` | `pseudoGC/Phylo/m2.motu.txt` and `m2.<rank>.txt` |
| Protal mode 1 | `<sample>/Tax/Protal/profiles/<sample>.profile` | `pseudoGC/protal_singular/Protal.abundance.tsv` |
| Protal mode 2 | `pseudoGC/protal/profiles/` | `pseudoGC/protal/Protal.abundance.tsv` and retained `strains/` MSAs |

## Troubleshooting checks

- If MetaPhlAn reports an unknown option, confirm version 4.2.6 and inspect the
  generated `LOGandSUB/metaPhl4.sh` for `--db_dir`, not `--bowtie2db`.
- If mOTUs reports that `-o` is required, an older generated job script is being
  reused; regenerate/resubmit the profiling job with the current `MATAF4.pl`.
- If mOTUs cannot load its database, `motus2_DB` must be the directory directly
  above `db_mOTU/`, and that child should contain `mOTUsv4.1.versions`.
- If RiboFind tries to build an index and fails, check write access beside the
  configured SSU/LSU FASTA files.
- If functional profiling stops in `getSpecificDBpaths()`, check both the
  directory key for the selected alias and the exact reference filename in the
  table above.
- If Protal is skipped unexpectedly, inspect the per-sample `.Protal.skip` file
  and rerun with `-protalIgnoreErrors 0` to obtain strict validation details.
