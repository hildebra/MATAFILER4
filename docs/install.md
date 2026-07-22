<!-- Documentation navigation -->
[Home](../README.md) | [Quick start](quickstart.md) | [Installation](install.md) | [Configuration](configuration.md) | [Mapping files](mapping_files.md) | [Workflows](common_workflows.md) | [Outputs](outputs.md) | [Flag reference](flag_reference.md) | [FAQ](FAQ.md) | [Glossary](glossary.md)

---

# Installing MATAFILER4

## Overview

MATAFILER4 is installed with the bundled installer script:

```bash
bash helpers/install/installer.sh
```

The installer uses `micromamba` to create and update the conda environments required by MATAFILER4. It can be rerun after major MATAFILER4 updates; existing environments are updated rather than fully reinstalled.

## Requirements

Before running the installer, make sure you have:

- a Linux/HPC environment
- `git`
- `bash`
- `micromamba` available in your `PATH`
- internet access during installation, unless environments and databases are managed manually by your local system administrators

The installer will abort if `micromamba` cannot be found.

## Basic installation

Clone the repository:

```bash
git clone https://github.com/hildebra/MATAFILER4.git
cd MATAFILER4
```

Run the installer:

```bash
bash helpers/install/installer.sh
```

The first installation can take a long time because several environments and databases are created or downloaded.

Legacy `MGTK*` environments are reported but left untouched by default. Remove them
explicitly only after reviewing their contents:

```bash
bash helpers/install/installer.sh --remove-legacy-envs
```

To force a fresh CheckM2 and MetaPhlAn database download:

```bash
bash helpers/install/installer.sh --refresh-databases
```

After the installer finishes, reload your shell configuration or start a new shell:

```bash
source ~/.bashrc
```

Activate the main MATAFILER4 environment:

```bash
micromamba activate MF4
```

Check the installation:

```bash
./MATAF4.pl -checkInstall
```

## What the installer does

The installer performs the following main actions:

1. Checks that `micromamba` is available.
2. Creates `config.txt` if it does not already exist, based on the bundled MATAFILER configuration template.
3. Adds the following entries to `~/.bashrc`, unless they are already present:

```bash
export MF4DIR=/path/to/MATAFILER4/
export PERL5LIB="$PERL5LIB:/path/to/MATAFILER4/"
```

4. Reports old legacy `MGTK*` environments. They are removed only when
   `--remove-legacy-envs` is supplied.
5. Creates or updates the MATAFILER4 conda environments, including:

| Environment | Main purpose |
|---|---|
| `MF4` | Main MATAFILER4 environment |
| `MF4gtdbtk` | GTDB-Tk-related tools |
| `MF4semibin` | SemiBin |
| `MF4binners` | Additional binning tools |
| `MF4genomeface` | Optional GenomeFace functionality; uses an external NERSC package channel |
| `MF4scgbinner` | SCG-based binning |
| `MF4checkm2` | CheckM2 and MetaPhlAn dependencies |
| `MF4phylo` | Phylogenetic tools |
| `MF4_R` | R-based helper scripts |

`MF4genomeface` is best-effort because its package metadata is hosted outside
conda-forge and Bioconda. If that external channel is unavailable, the installer
prints a warning and continues; rerun it later to install or update GenomeFace.

6. Clones a pinned `extract_gtdb_mg` revision into `gits/XGTDB/` if missing and
   verifies existing checkouts before use.
7. Downloads selected databases where possible, including:
   - the `hostile` human reference index `human-t2t-hla`, if `hostile` is available
   - the CheckM2 database
   - the MetaPhlAn database

Successful CheckM2 and MetaPhlAn downloads receive tool-version marker files.
Missing or mismatched markers cause the installer to try the database setup again;
`--refresh-databases` forces a fresh download. Database downloads are best-effort:
network or permission failures produce warnings but do not abort installation of the
software environments. If a database directory belongs to the current user but is
missing its owner write/search bits, the installer repairs those bits before retrying.

## Configuration after installation

After installation, check `config.txt` in the MATAFILER4 directory. At minimum, confirm these paths:

```text
MFLRDir      /path/to/MATAFILER4/
DBDir        /path/to/MATAFILER4/data/DBs/
globalTmpDir /path/to/global/scratch/
nodeTmpDir   /path/to/node/local/tmp/
```

For HPC use, `globalTmpDir` should be visible from all compute nodes, while `nodeTmpDir` should point to node-local temporary storage if available. Correct temporary-directory configuration is important for performance and for avoiding excessive I/O on shared storage.

## GTDB and GTDB-Tk databases

GTDB and GTDB-Tk databases are required for MAG classification and GTDB marker-gene workflows. The installer prints the recommended download command at the end of installation.

Example:

```bash
cd helpers/install
./get_gtdb.py all -v 226 -t /path/to/download -d /path/to/extract/to --tk split
```

The temporary download directory passed with `-t` can be removed after the database has been extracted and configured. See:

```bash
./get_gtdb.py -h
```

for additional options, including workflows where download and extraction need to be run separately.

## Updating MATAFILER4

To update MATAFILER4 itself:

```bash
cd /path/to/MATAFILER4
git pull
```

Then rerun the installer to update environments if dependency definitions changed:

```bash
bash helpers/install/installer.sh
```

This is expected usage and is usually much faster than the initial installation.

## Troubleshooting installation

### `micromamba could not be found`

Install micromamba and ensure it is available in your `PATH`, then rerun the installer.

### Old `MG-TK` entries in `.bashrc`

The installer checks for old `MG-TK` shell configuration blocks. If it finds them, remove the lines after the marker:

```text
##------------> MG-TK ADDED
```

Then rerun the installer.

### Conda dependency conflicts

The installer already uses flexible channel priority for micromamba environment creation. If conflicts persist, update micromamba and rerun the installer. On managed HPC systems, it may be preferable to ask local support to inspect the failing environment YAML file in `helpers/install/`.

### CheckM2 or MetaPhlAn database download fails

The database downloads are optional during software installation. A failure leaves
the corresponding version marker absent and the installer continues, so rerunning it
will try again. For a MetaPhlAn permission warning, inspect both Unix mode bits and
any HPC filesystem ACL:

```bash
ls -ld /path/to/MATAFILER4/data/DBs/MP4
getfacl /path/to/MATAFILER4/data/DBs/MP4  # if getfacl is available
```

If you own the directory, `chmod u+rwx /path/to/MATAFILER4/data/DBs/MP4` is normally
sufficient. Otherwise use a database location you own or ask the directory owner or
cluster administrator to grant access.

Activate the relevant environment and retry manually:

```bash
micromamba activate MF4checkm2
checkm2 database --download --path /path/to/MATAFILER4/data/DBs/CM2/
metaphlan --install --bowtie2db /path/to/MATAFILER4/data/DBs/MP4/
```

### Installation was interrupted

Rerun:

```bash
bash helpers/install/installer.sh
```

The installer is designed to create missing environments and update existing ones.

## After installation

A minimal sanity check is:

```bash
source ~/.bashrc
micromamba activate MF4
./MATAF4.pl -checkInstall
```

Then continue with the [quick start](quickstart.md) or the [common workflows](common_workflows.md).
