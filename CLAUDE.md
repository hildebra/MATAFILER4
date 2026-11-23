# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

MATAFILER4 is a Perl pipeline for shotgun metagenomics on Linux HPC clusters (Slurm mainly; SGE and LSF are also supported). A full run needs the micromamba environments from `helpers/install/`, large databases and a scheduler, so a development machine cannot run it end to end. Code changes are checked with the unit tests in `t/` and with dry runs (`-submit 0`). User documentation is in `docs/` (start at `README.md`).

## Language rule

Write Perl, not Python: no new Python scripts, tests, helpers or mock tools. Tests use `Test::More`; stand-in tools in tests are Perl or shell scripts.

## Commands

```bash
perl helpers/runTests.pl                       # all tests in t/, 4 files in parallel
perl helpers/runTests.pl strain mgs            # only test files whose name matches
perl helpers/runTests.pl -jobs 1 -v sdm_cleaner  # one file, full TAP output
perl t/checkpoint.t                            # one file directly, from the repository root
perl helpers/runTests.pl -list                 # list test files
```

- The runner needs `TAP::Harness`. The Perl in Git Bash on Windows lacks it, so run the suite in WSL or on Linux. Running a single file with `perl t/<name>.t` only needs `Test::More`.
- There is no build step and no linter. `t/script_compile.t` compiles every `.pl`/`.pm` in the repository.
- Tests that need a compiled helper from `bin/` (`sdm`, `rtk2`, `rdCover`, `MSAfix`, `vcf2fna`, `LCA`, `clusterMAGs`, `cc.bin`: Linux binaries) skip themselves when it cannot run. The runner's summary lists them.
- A pipeline dry run checks paths and configuration without submitting jobs:
  `perl MATAF4.pl -map run.map -assembleMG 2 -requireInput 1 -submit 0 -from 0 -to 1`.
  `-inspectState 1` and `-planState 1` print read-only JSON reports of what is complete and what would be submitted.
- `perl MATAF4.pl -checkInstall` checks a cluster installation.

## Architecture

### Entry points and pipeline order

1. `MATAF4.pl` (about 12,500 lines) processes each sample in the map file: read download and cleaning (`sdm`), host filtering, assembly (single, co-assembly groups, hybrid with long reads), gene prediction, mapping, contig stats, binning, and profiling (Kraken2, MetaPhlAn, mOTUs, Protal, RiboFind, DIAMOND).
2. `secScripts/geneCat.pl` builds a gene catalogue across samples, with abundance matrices and functional annotation. It starts `MGS.pl`.
3. `secScripts/MGS.pl` builds MAGs/MGS. It starts `secScripts/MGS/strain_within.pl` (strains within an MGS) and `secScripts/MGS/phylo_MGS_between.pl` (tree between MGS); both build their trees with `secScripts/phylo/buildTree5.pl`.

`secScripts/` also holds many small worker scripts that these entry points write into job scripts. `helpers/` holds the installer and stand-alone utilities.

### Controller model: rerun until complete

`MATAF4.pl` does almost no heavy work itself. It loops over the samples in the map file (`for ($JNUM=$from; $JNUM<$to; ...)`). For each step it checks completion markers, writes a job script into the sample's `LOGandSUB/` directory, and submits it through `Mods::Subm::qsubSystem` with scheduler dependencies on earlier jobs. Users rerun the same command after the jobs finish, or use `-loopTillComplete`, until every sample is done. This affects every change to a step:

- A step counts as done when its marker exists. Markers are mostly `*.sto` files (e.g. `ass.done.sto`, `<cram>.sto`), plus `Mods::Checkpoint` manifests and the per-sample `MF4.sentinel.*` records in `Mods::SampleCompletion`. A job must write its marker last, and only when it succeeded, otherwise reruns skip broken output.
- A sample with a lock file in `LOGandSUB/` is skipped while its recorded jobs are still queued or running.
- Samples in the same assembly group (`AssmblGrps` column) are co-assembled when the loop reaches the group's last member. The output goes to `<OutPath>/<RunID>/AssmblGrp_<group>/`.
- `Mods::JobGraph` is a supervised alternative, used mainly by geneCat. Jobs are recorded instead of submitted. Each job is submitted only after its prerequisites completed, and its outcome is read from a marker it writes. A job killed for lack of memory is resubmitted with 1.5x memory.
- `Mods/Workflow*.pm` (State, Plan, Apply, Runner, Control, Resilience) hold logic taken out of the controller: state inspection, planning, loop windows, dependency handling, atomic writes, retries and preflight checks. The unit tests exercise these modules, because the controller loop itself cannot run without a cluster.

### Configuration and tool paths

Do not hard-code tool or database paths. Look them up with `getProgPaths("key")` from `Mods/IO_Tamoc_progs.pm`. It reads three tab-separated files, and the first definition of a key wins:

1. The site config `Mods/MATAFILERcfg.txt`, or the file given with `-config`. It is gitignored; the installer creates it from `Mods/config.old` and links it as `config.txt` at the repository root.
2. `Mods/config_internal.txt`: the repository's own scripts and binaries (e.g. `geneCat_scr`, `MAGpipe`, `sdm`).
3. `Mods/config_DBs.txt`: database locations.

Values can contain `[MFLRDir]`, `[DBDir]`, `[BINDir]`, `[Rscript]` and similar placeholders, and `$ENV_VAR`. An entry with `env:<name>` comes back as a micromamba activation line followed by the command. Such a value is meant to be pasted into a job script, not used as a file path.

### Adding tools, environments and databases

Reuse what MATAFILER4 already installs before adding anything new:

- **Check first.** Look for an existing key in `Mods/config_internal.txt` or `Mods/config_DBs.txt`, a binary in `bin/`, or a package in one of the micromamba environments defined in `helpers/install/*.yml` (`MF4`, `MF4_R`, `MF4phylo`, `MF4binners`, ...). On an installed system, `perl MATAF4.pl -checkInstall` reports what is present.
- **Read tools from the config.** Call every tool through `getProgPaths("key")`. A new tool gets a key in `Mods/config_internal.txt` (a database: `Mods/config_DBs.txt`), with `env:<environment>` when it runs in an MF4 environment. Add its package to an existing environment's `.yml` where it fits; a new environment also needs an `ensure_environment` line in `helpers/install/installer.sh`. Add essential tools to the `-checkInstall` list (`checkMFFInstall` in `Mods/TamocFunc.pm`). A new or renamed repository script also gets a key; `t/config_integrity.t` checks that every key the code looks up resolves and that every repository path in the configs exists.
- **Link an external toolkit through one directory key**, as is done for MG-STK: the site config sets `MGSTKDir` once (`MGSTKDir	$MF4DIR` in `Mods/config.old`), and `config_internal.txt` refers to its scripts as `[MGSTKDir]/strainStats.R` with `env:MGSTK`. A new directory key of this kind must also be added to `loadConfigs` in `Mods/IO_Tamoc_progs.pm` (the foundation-key list and the placeholder substitution).
- **Find the MATAFILER4 directory with `$MF4DIR`.** The installer exports it in `~/.bashrc`. Config templates use it (`MFLRDir	$MF4DIR` and `DBDir	$MF4DIR/data/DBs/` in `Mods/config.old`), as do the command examples in `docs/` (`perl $MF4DIR/MATAF4.pl ...`). Inside config values, refer to the directory as `[MFLRDir]`.

### Command-line flags

`docs/flag_reference.md` is the single source of `-help` text for `MATAF4.pl`, `geneCat.pl`, `MGS.pl`, `strain_within.pl`, `strain_within_2.2.pl` and `buildTree5.pl` (rendered by `Mods::FlagReference`). When you add, rename or remove an option in a script's `GetOptions` block (for `MATAF4.pl`, in `sub getCmdLineOptions`), update that script's section in the reference. `t/command_line_help.t` fails if the two disagree.

### Map files

A run is described by a tab-separated map file. Its header row starts with `#SmplID`, and lines such as `#OutPath`, `#RunID` and `#DirPath` set run-wide values. Sample output goes to `<OutPath>/<RunID>/<SmplID>/`. `Mods::GenoMetaAss::readMap` parses it; the format is in `docs/mapping_files.md`.

## Conventions

- **Tests**: put new tests in `t/<topic>.t` with `Test::More`. A test should run a function or script on a small fixture and check the result, not match source text. `t/lib/MFTestConfig.pm` points the configuration at `t/MATAFILERcfg.txt` and puts the shims in `t/bin/` (e.g. a fake `pigz`) on `PATH`. Load the repository's modules with `use lib "$Bin/.."` (and `"$Bin/lib"`), so tests do not depend on the caller's include path.
- **No site-specific paths**: `t/private_paths.t` fails if institute paths (cluster home or project directories, personal temporary or user directories) appear in `MATAF4.pl`, `Mods/`, `helpers/`, `secScripts/`, `docs/` or `examples/`. Site notes belong in `docs/local_hpc_notes.md` or the local config.
- **Line endings are LF** (`.gitattributes` sets `eol=lf`). Generated job scripts run on Linux, so do not add CRLF.
- **Versions and changelog**: each entry point has its own version: `$MATFILER_ver` plus the version-history comment block near the top of `MATAF4.pl`, `our $version` in `geneCat.pl`, `$MGSpipelineVersion` in `MGS.pl`, `$version` in `buildTree5.pl`. For a release-worthy change, bump the version, add an entry to `doc/changelog.md` (problem, fix, tests), and name the version in the commit subject (e.g. `geneCat 0.69: ...`).
- **Docs**: when behaviour, outputs or options change, update the matching page in `docs/` in the same change.
- **Branches**: `main` is the default branch; current development is on `dev2`.
