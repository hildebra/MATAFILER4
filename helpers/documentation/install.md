
## Installing MATAFILER4

<details>
  <summary> Expand Installation section </summary>


### Requirements

MATAFILER4 requires a perl installation and sdm requires a fairly recent C++ compiler (like gcc or clang) that supports C++11; these will be automatically installed in the install script.
MATAFILER4 currently only works under linux, and is expected to run on a computer cluster. Since the pipeline includes a lot of external sofware, you will need fully installed Micromamba ([https://mamba.readthedocs.io/en/latest/installation.html](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)).


### Installation


MATAFILER4 can be downloaded directly from Github, using:

```bash
git clone https://github.com/hildebra/MATAFILER4.git
```

MATAFILER4 comes with an installation script, that uses micromamba. Ensure you have micromamba installed for your account on a linux HPC. Then run:

```bash
bash helpers/install/installer.sh
```

This will guide you through the installation (should run completely automatic) and requires internet access. Since a lot of packages will be installed, this can take an hour or longer. All required software will be downloaded and installed in the Conda/Mamba directories.

If you are having issues with package conflicts when `installer.sh` is creating environments, trying setting your channel priority to flexible: `micromamba config set channel_priority flexible`

Last, you can run 

```bash
# Activate MF4 environment every time you need to run MATAFILER4
micromamba activate MF4

./MATAF4.pl -checkInstall
```

to check that some essential programs have been correctly installed and are available in the expected environments.
Note that this is only a subset of programs, but should cover most use cases of MATAFILER4. (This will also automatically run after each installation of MATAFILER4)


### Updating MATAFILER4

MATAFILER4 will be frequently updated. To get the latest version, go to your MATAFILER4 directory and run
```
git pull
```
Sometimes new packages will be included or program versions modified. To obtain these changes, run the install script again (this will update existing environments - no worries, this is not a complete reinstall):
```
bash helpers/install/installer.sh
```

### Preparing MATAFILER4

- follow installation process (essentially `git clone https://github.com/hildebra/MATAFILER4.git` & run `bash helpers/install/installer.sh` )

- After the installation is complete, you will find the file named: "config.txt" inside of the MATAFILER4 directory. This is the main file where you have defined all the paths for directories and slurm configuation. Always check in order to ensure that all directories are correct: 

    - MFLRDir	`/path/to/your/MATAFILER4/installation/`
    - DBDir	`/path/to/your/database_dir/`

- change tmp dir (scratch space) with project scratch folder:

    - globalTmpDir	`/path/to/your/scratch/`
	- nodeTmpDir	`/path/on/node/to/tmp` -> on slurm systems this could be a variable, e.g. `$SLURM_LOCAL_SCRATCH/MATAF4/`

- follow either the assembly-dependent or assembly-independent example runs (Examples section below)

### Download databases

#### GTDB & GTDBtk

GTDB and the GTDBtk database are used for MAG classification in the gene catalog step.

After installing MATAFILER4, there is a utility script `helpers/install/get_gtdb.py` provided to download these databases and format them appropriately for MATAFILER4.

To download and extract the r226 databases, and configure MATAFILER4 to use them, run

```
./get_gtdb.py all -v 226 -t /path/to/download/to -d /path/to/extract/to --tk split
```

You can delete the download directory (`-t`) after everything is correctly
configured.

If the system you run MATAFILER4 on does not have internet access, the download process can be done separately.
See `./get_gtdb.py -h` for help on running these steps separately.

### Useful configurations to track and check on MATAFILER4 jobs

The most common reason why MATAFILER4 jobs fail are related to node configurations (available ram, hdd space, CPUs). There are several aliases that are useful in checking on slurm jobs that are running on your local HPC, understanding how MATAFILER4 processes your samples and fixing errors. Thus following up jobs and checking their error logs is essential in understanding limitations in your current environment and get your metagenomes processed effectively, as listed below:

These aliases can be directly added to your ~/.bashrc (just make sure the .bashrc is loaded):

```{sh}
#list running jobs with more relevant info
alias sq='squeue -u $USER -o "%8i %.4P %.14j %.2t %8M %.3C %.15R %20E"'
#check where job bash, std output, error output is stored, dependencies etc
alias si='scontrol show job'
#delete jobs that have DependencyNever status
alias scDN="squeue -u $USER | grep dencyNev | cut -f11 -d' ' | xargs  -t -i scancel {}"
#show the number of jobs currently running for different users on your cluster; useful for estimating how busy the HPC currently is
alias busy="squeue | sed -E 's/ +/\t/g' | cut -f5 | sort | uniq -c | sed -E 's/ +//' | sort -k1 -n -t' '"
#show output log of job
sio() {
JID=$1
if test "$#" -eq 0; then
JID=$(squeue -u hildebra | grep $USER | grep -v 'interact' | awk '{$1=$1};1' | cut -f1 -d' ' | head -1)
fi
cat $(scontrol show job $JID | grep 'StdOut' | sed 's/.*=//g')
}
#show error log of job
sie() {
JID=$1
if test "$#" -eq 0; then
JID=$(squeue -u $USER | grep -v 'interact' | awk '{$1=$1};1' | cut -f1 -d' ' | head -1)
fi
cat $(scontrol show job $JID | grep 'StdErr' | sed 's/.*=//g')
}
#show bash script (commands) of job
sis() {
JID=$1
if test "$#" -eq 0; then
JID=$(squeue -u $USER | grep -v 'interact' | awk '{$1=$1};1' | cut -f1 -d' ' | head -1)
fi
cat $(scontrol show job $JID | grep 'Command' | sed 's/.*=//g')
}
```




## I/O considerations

### I/O (Input/Output): important considerations and design decisions


Analysing a shotgun metagenomic experiment can be a computationally extremely demanding task, as in some experiments several TB of data can be accumulating. MATAFILER4 was designed with the latter case in mind, but can of course also handle smaller experiment.  
In order to be able to cope with these data amounts, a lot of 'file juggling' is happening behind the scenes. A lot of temporary files are being created that don't need to be saved on long term storage solution that are backed up and generally also slower. For this purpose big servers usually have a 'scratch' dir that is the global temporary storage on which all nodes in a cluster can write, but that is not backed up and might be cleaned infrequently. Further, usually each node has a local temp dir, to which only that specific node has access. Using these temporary solutions does make the whole cluster more stable and also enable other users to use a cluster more efficiently. To give you an example: if you have an IO heave process like searching with diamond through a lot of reads, you will use up the bandwith provided by your permanent storage very quickly. This could lead to situations where 500 cores on the cluster are busy with running in parallel diamond searches, but since the IO is so severely limited, only a small fraction of data trickles through to these jobs, effectively maybe giving 16 cores work. In this case the cluster would be unnecessarily blocked and the 500 core job would also take much longer than needed. That is the reason why file juggling is so important and why so much development effort went into optimizing this for MATAFILER4.  
To take advantage of this, I strongly recommend to ask your sysadmin where the local and global temp storage on your cluster are and set in the MATAFILER4 config the variables 'globalTmpDir' and 'nodeTmpDir' variables correspondingly.




## Known issues

This is a beta release of MATAFILER4. Some parts of the pipeline will currently not run, because we have not started yet linking in the various databases being used. Known DBs missing:
- LSU/SSU DBs ((needed for miTag approaches, flag -profileRibosome )
- all functional annotation databases (needed in gene catalog step or flag -profileFunct )
</details>


