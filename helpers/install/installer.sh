#!/usr/bin/bash
#installer script for MG-TK
#stay in helpers/install/ dir while executing

#some basic housekeeping
#set -e
ulimit -c 0;
#set

echo "MATAFILER4 installer script"
echo "This script will install several conda environments with most dependencies for MATAFILER4. It will use micromamba to run the installations, please ensure micromamba is installed natively for your account.";
echo "You can (and probably should) rerun this script every time you pull a major MF update. Changes to dependencies will be automatically updated and rerunning the script will be significantly faster than running it the first time.";

if ! command -v micromamba &> /dev/null
then
    echo "micromamba could not be found"
	echo "Make sure micromamba is in your \$PATH"
	echo "Aborting"
    exit
fi

if [ -z "${MAMBA_EXE}" ] ; then
MAMBA_E=$MAMBA_EXE
else
MAMBA_E=micromamba
fi
#which micromamba


eval "$($MAMBA_E shell hook --shell=bash)"

find_in_mamba_env(){
	$MAMBA_E env list | grep "${@}" >/dev/null 2>/dev/null
}

find_in_bashrc(){
	grep "${@}" ~/.bashrc >/dev/null 2>/dev/null
}


echo "Using micromamba version:"
$MAMBA_E --version
#exit


SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
MFdir=$(realpath -s $SCRIPT_DIR/../..)

#remove stone that declares programs are checked by MG-TK
rm -f $MFdir/helpers/install/progsChecked.sto
touch $MFdir/helpers/install/runningInstall.sto

mkdir -p $MFdir/gits/ 
INSTdir=$MFdir/helpers/install/
DBdir=$MFdir/data/DBs/


#MFLRDir	$MF3DIR
if [ ! -f $MFdir/config.txt ] ; then
	cp -f $MFdir/Mods/config.old $MFdir/Mods/MATAFILERcfg.txt
	#should be defaulted to $MF3DIR now
	#sed -i "s+MFLRDir.*+$\t$MFdir+" $MFdir/Mods/MATAFILERcfg.txt
	ln -s $MFdir/Mods/MATAFILERcfg.txt $MFdir/config.txt
	echo "Rewrote config.txt. Please modify as needed to local paths"
fi


if find_in_bashrc "##------------> MG-TK ADDED" ; then
	echo "It seems your ~/.bashrc still contains an old MG-TK install (lines after "
	echo "\"##------------> MG-TK ADDED\" )."
	echo "Please delete all MG-TK added lines from your bashrc, as this is incompatible with MATAFILER"
	echo "Rerun MATAFILER installer after this is done"
	echo "      Exiting.."
	exit 
fi


if ! find_in_bashrc "##------------> MF4 ADDED" ; then
	printf "\n\n##------------> MF4 ADDED <----------##\nexport MF4DIR=$MFdir/\nexport PERL5LIB=\"\$PERL5LIB:$MFdir/\"\n##------------> MF4 ADDED <----------##\n\n" >> ~/.bashrc
	echo "Added MATAFILER modules to .bashrc"
fi



#potential old MGTK versions .. uninstall environments

envsToDEL='MGTK MGTK_R MGTKbinners MGTKcheckm2 MGTKgtdbtk MGTKphylo MGTKsemibin MGTKwhokar'
for etd in $envsToDEL; do
	if find_in_mamba_env $etd ; then
		$MAMBA_E  env remove -y -n  $etd 
	fi 
done


#exit;


# For all micromamba installs, we use --channel-priority 1. This sets to 
# "flexible" priority - from conda docs "the solver may reach into lower 
# priority channels to fulfill dependencies, rather than raising an 
# unsatisfiable error". This should help avoid a lot of dependency problems
# during install, as micromamba defaults to "strict".
CHNLprio=1



if ! find_in_mamba_env "MF4\s" ; then
	echo "Creating base MATAFILER conda environment.. This might take awhile"
	#first install spades that seems to require a lot of mem..
	#$MAMBA_E create --channel-priority $CHNLprio -q -y -n MFF spades
#just in case it crashes, this often recovers it..
	$MAMBA_E create --channel-priority $CHNLprio -q -y -f $INSTdir/MF4.yml #-q -y
	$MAMBA_E activate MF4
	#echo "Installing R packages in MGTK environment";	Rscript $INSTdir/reqPackages.R
	#pip install biopython
else 
	echo "Updating base MATAFILER conda environment.. Please be patient"
					#	$MAMBA_E activate MGTK
	$MAMBA_E install --channel-priority $CHNLprio -q -y -f $INSTdir/MF4.yml #-q -y
	
	#echo "Updating R packages in MGTK environment"
	#{ Rscript $INSTdir/reqPackages.R
	#} || { 		echo "Rscript install failed.. trying direct excecution";		$INSTdir/./reqPackages.R;	}
fi

#exit

$MAMBA_E activate MF4

#took prostT5 out for now.. not used by default
#if command -v foldseek &> /dev/null ; then
#	echo "Preparing prostT5 for foldseek.."
#	#prepare foldseek search
#	if [ ! -d $DBdir/PtostT5_W ];then
#		foldseek databases ProstT5 $DBdir/PtostT5_W $DBdir/tmp;
#	fi
#fi

if command -v hostile &> /dev/null ; then
	echo "Installing hostile human reference database"
	export HOSTILE_CACHE_DIR=$DBdir/hostile/
	hostile index fetch --name human-t2t-hla
	#run a second time in cases it crashes..
	hostile index fetch --name human-t2t-hla
	
	
fi
echo ""
echo "Installing/updating further dependencies in additional conda environments.."
echo ""
echo "" 


#if ! find_in_mamba_env "checkm2" ; then
#	#git clone https://github.com/chklovski/CheckM2.git $MFdir/gits/checkm2/
#	echo "Installing checkm2"
#	$MAMBA_E create -y -q -f $INSTdir/checkm2.yml -n checkm2
#	$MAMBA_E activate checkm2
#	pip3 install --upgrade pip

#	pip install CheckM2 packaging
#	$MAMBA_E activate MFF
#	echo "checkm2 installed. you can verify this environment by running \"micromamba activate checkm2\" and \"checkm2\""
#else
#	$MAMBA_E activate checkm2
#	if [ ! command -v checkm2 > /dev/null 2>&1 ]; then
#		echo "Could not find checkm2. Please install via \"micromamba activate checkm2;pip install CheckM2 packaging\" and restart installer"
#		exit 5
#	fi
#	$MAMBA_E activate MFF
#fi

if [ ! -f "$SCRIPT_DIR/../../gits/XGTDB/extract_gtdb_mg.py" ]; then
	echo "Installing extractGTDB into $MFdir/gits/XGTDB/"
	git clone https://github.com/4less/extract_gtdb_mg.git $MFdir/gits/XGTDB/
fi

#additional dependencies not in the main yml..
if ! find_in_mamba_env "MF4gtdbtk" ; then
	echo "Creating MF4gtdbtk environment"
	$MAMBA_E create --channel-priority $CHNLprio -q -y -f $INSTdir/GTDBTK.yml 
else 
	echo "Updating MF4gtdbtk environment"
	$MAMBA_E install --channel-priority $CHNLprio -q -y -f $INSTdir/GTDBTK.yml 
fi

if ! find_in_mamba_env "MF4semibin" ; then
	echo "Creating MF4semibin environment"
	$MAMBA_E create --channel-priority $CHNLprio -q -y -f $INSTdir/SemiBin.yml 
else 
	echo "Updating MF4semibin environment"
	$MAMBA_E install --channel-priority $CHNLprio -q -y -f $INSTdir/SemiBin.yml 
fi


if ! find_in_mamba_env "MF4binners" ; then
	echo "Creating MF4binners environment"
	$MAMBA_E create --channel-priority $CHNLprio -q -y -f $INSTdir/Binners.yml
else 
	echo "Updating MF4binners environment"
	$MAMBA_E install --channel-priority $CHNLprio -q -y -f $INSTdir/Binners.yml
fi

if ! find_in_mamba_env "MF4genomeface" ; then
	echo "Creating MF4genomeface environment"
	$MAMBA_E create --channel-priority $CHNLprio -q -y -f $INSTdir/MF4genomeface.yml
else 
	echo "Updating MF4genomeface environment"
	$MAMBA_E install --channel-priority $CHNLprio -q -y -f $INSTdir/MF4genomeface.yml
fi

if ! find_in_mamba_env "MF4scgbinner" ; then
	echo "Creating MF4scgbinner environment"
	PIP_USER=false $MAMBA_E create --channel-priority $CHNLprio -q -y -f $INSTdir/SCGBinner.yml
else
	echo "Updating MF4scgbinner environment"
	PIP_USER=false $MAMBA_E install --channel-priority $CHNLprio -q -y -f $INSTdir/SCGBinner.yml
fi

CM2DB=$DBdir/CM2/
MP4DB=$DBdir/MP4/

if ! find_in_mamba_env "MF4checkm2" ; then
	echo "Creating MF4checkm2 environment"
	$MAMBA_E create --channel-priority $CHNLprio -q -y -f $INSTdir/checkm2.yml 
	$MAMBA_E activate MF4checkm2
	#install checkM2 DB
	checkm2 database --download --path $CM2DB
	#install metaphlan4 DB
	#currently there is a bug: can't install into a given dir via conda
	metaphlan --install --bowtie2db $MP4DB
	$MAMBA_E deactivate
else 
	echo "Updating checkm2 environment"
	$MAMBA_E install -y -q -f $INSTdir/checkm2.yml 
	if [ ! -d $CM2DB ]; then
		$MAMBA_E activate MF4checkm2
		checkm2 database --download --path $CM2DB
		$MAMBA_E deactivate
	fi
	if [ ! -d $MP4DB ]; then
		$MAMBA_E activate MF4checkm2
		metaphlan --install --bowtie2db $MP4DB
		$MAMBA_E deactivate
	fi
	
fi


#if ! find_in_mamba_env "comeBin" ; then
#	echo "Installing comeBin environment"
#	$MAMBA_E create -y -q -f $INSTdir/comeBin.yml
#else 
#	echo "Updating gtdbtk environment"
#	$MAMBA_E update -y -q -f $INSTdir/comeBin.yml 
#fi

#additional dependencies not in the main yml..
#if ! find_in_mamba_env "MGTKmetaMDBG" ; then
#	echo "Installing MGTKmetaMDBG environment"
#	$MAMBA_E create --channel-priority 0 -q -y -f $INSTdir/metaMDBG.yml
#else 
#	echo "Updating MGTKmetaMDBG environment"
#	$MAMBA_E activate MGTKmetaMDBG
#	$MAMBA_E update --channel-priority 0 -q -y -f $INSTdir/metaMDBG.yml 
#	$MAMBA_E deactivate 
#fi

#if ! find_in_mamba_env "motus" ; then
#	echo "Installing motus environment"
#	$MAMBA_E create -y -q -f $INSTdir/motus.yml
#else 
#	echo "Updating motus environment"
#	$MAMBA_E update -y -q -f $INSTdir/motus.yml 
#fi

if ! find_in_mamba_env "MF4phylo" ; then
	echo "Installing MF4phylo environment"
	$MAMBA_E create --channel-priority $CHNLprio -q -y -f $INSTdir/phylo.yml 
else 
	echo "Updating MF4phylo environment"
	$MAMBA_E install --channel-priority $CHNLprio -q -y -f $INSTdir/phylo.yml 
fi

#if ! find_in_mamba_env "MGTKwhokar" ; then
#	echo "Installing MGTKwhokar environment"
#	$MAMBA_E create --channel-priority $CHNLprio -q -y -f $INSTdir/whokaryote.yml 
#else 
#	echo "Updating MGTKwhokar environment"
#	$MAMBA_E install --channel-priority $CHNLprio -q -y -f $INSTdir/whokaryote.yml 
#fi

if ! find_in_mamba_env "MF4_R" ; then
	echo "Installing MF4_R environment"
	$MAMBA_E create --channel-priority $CHNLprio -q -y -f $INSTdir/MGTK_R.yml 
else 
	echo "Updating MF4_R environment"
	$MAMBA_E install --channel-priority $CHNLprio -q -y -f $INSTdir/MGTK_R.yml 
fi




#if ! find_in_mamba_env "Rbase" ; then
#	$MAMBA_E create -y -f $INSTdir/Rbase.yml
#	$MAMBA_E activate Rbase
#	Rscript reqPackages.R
#else 
#	$MAMBA_E update -y -f $INSTdir/Rbase.yml --prune --allow-uninstall
#fi

#later toadd..
#git clone https://github.com/GaetanBenoitDev/metaMDBG.git;mima create -y -f conda_env.yml;activate metaMDBG

rm -f $MFdir/helpers/install/runningInstall.sto





echo ""
echo ""
echo "How to download GTDB & GTDBtk databases"
echo ""
echo "The database for GTDB and GTDBtk are needed for MAG classification."
echo "These can be downloaded using the script 'helpers/install/get_gtdb.py' i.e."
echo "    ./get_gtdb.py all -v 226 -t /path/to/download -d /path/to/extract/to --tk split"
echo "will download and extract these databases, and configure MG-TK to use them."
echo "See script help (./get_gtdb.py -h) for more information on usage."
echo ""
echo ""

echo "Finished MG-TK install"
echo ""
echo "To run MG-TK, make sure you are in the MGTK environment (micromamba activate MGTK)."
echo "You can rerun the installer.sh anytime, to ensure package were installed or are being updated."
echo "Run \"MG-TK.pl -checkInstall\" to ensure that the installation was successful."

exit 
