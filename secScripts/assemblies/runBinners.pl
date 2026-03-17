#!/usr/bin/env perl

use warnings;
use strict;
use File::Basename;
use Getopt::Long qw( GetOptions );

use Mods::GenoMetaAss qw(systemW);
use Mods::Binning qw (runMetaBat runCheckM runSemiBin runMetaDecoder getBinSubdirName runGenomeFace runSCGBinner );
use Mods::IO_Tamoc_progs qw(getProgPaths jgi_depth_cmd);
use Mods::Subm qw(qsubSystem emptyQsubOpt qsubSystemJobAlive );

#v0.1: 2.3.24: ini version FH
#v0.11: 23.4.24: added support for hybrid assemblies (i.e. 2 crams / sample)
#v0.12: 17.1.26: added Genome Face support
#v0.13: 12.3.26: GPU job submission for GenomeFace
my $version = 0.13;


my $DoMetaBat2 = "";
my $BinDir = "";
my $smplIDs1 = "";
my $nodeSpTmpD2 = "";
my $cAssGrp = "";
my $metaGassembly = "";
my $MB2coresL = 1;
my $pathsPre = "";
my $seqTec = "ill";
my $GPUused = 0;
my $giveSBenv = ""; #human_gut/dog_gut/ocean/soil/cat_gut/human_oral/mouse_gut/pig_gut/built_environment/wastewater/chicken_caecum/global
my $logDir = "";


#"$BinnerScr -binner $DoMetaBat2 -binD $Bindir -smplID $smplIDs1 -tmpD $nodeSpTmpD2 -assmbl $metaGassembly -assmblGrp $cAssGrp -cores $MB2coresL -smplDirs " . join(",",@paths) . " -seqTec $seqTec

GetOptions(
#Directories/files
	"binner=i"  => \$DoMetaBat2, #1: metabat2 2: semibin 3: metadecoder, 4: GenomeFace
	#"tmp=s" => \$tmpDir,
	"logDir=s"  => \$logDir,
	"binD=s" => \$BinDir, #dir where bin result will be stored
	"smplID=s" => \$smplIDs1, #id of current (assembled) sample
	"tmpD=s" => \$nodeSpTmpD2, #node specific (local) tmp dir
	"assmbl=s" => \$metaGassembly, #file path to assembly, fna
	"assmblGrp=s" => \$cAssGrp, #assembly group of current sample
	"cores=i" => \$MB2coresL , #num cores to be used (locally in this case)
	"smplDirs=s" => \$pathsPre, #dirs to MF output dirs, where the map etc will be looked up
	"seqTec=s" => \$seqTec,  #PB, ONT, ill etc
	"SB_env=s" => \$giveSBenv,
	"useGPU=i" => \$GPUused,
);




my @paths = split /,/,$pathsPre;
my $BinnerName = getBinSubdirName($DoMetaBat2);
if ($logDir eq ""){
	$logDir = $BinDir;
}

print "======================================================================\n";
print "   runBinners.pl v $version\n";
print "     found " . scalar @paths . " sample dirs for sample \"$smplIDs1\"\n";
print "     using assembly $metaGassembly\n";
print "     using $MB2coresL cores, binner \"$DoMetaBat2\" to outdir $BinDir\n";
print "     using $giveSBenv environment\n" if ($giveSBenv ne "");
print "======================================================================\n";



my $preCmd = "";#"mkdir -p $nodeSpTmpD2";
my $BinCmd = "";
system "rm -rf $nodeSpTmpD2; mkdir -p $nodeSpTmpD2" unless (-d $nodeSpTmpD2);
system "rm -rf $BinDir; mkdir -p $BinDir;";


#actual work happens here:

if ($DoMetaBat2 == 1){
	$preCmd .= jgi_depth_cmd(\@paths,$nodeSpTmpD2."/depth$cAssGrp",95,$MB2coresL,$metaGassembly);# unless (-e );
	my $tmp = runMetaBat("$nodeSpTmpD2/depth$cAssGrp.jgi.depth.txt",$BinDir,$smplIDs1,$metaGassembly,$MB2coresL);
	if ($tmp eq "" ){$preCmd="";
	} else {$BinCmd .= $tmp;}
	#"MB2";
} elsif ($DoMetaBat2 == 2){#SemiBin
	#could work with jgidepth, but only for single sample.. not for co-assembly :(
	#"$nodeSpTmpD2/depth$cAssGrp.jgi.depth.txt"
	$BinCmd .= runSemiBin("",$BinDir,$nodeSpTmpD2, $smplIDs1,$metaGassembly,$MB2coresL,\@paths, $seqTec, $giveSBenv);
}elsif ($DoMetaBat2 == 3){#MetaDecoder
	#could work with jgidepth, but only for single sample.. not for co-assembly :(
	#"$nodeSpTmpD2/depth$cAssGrp.jgi.depth.txt"
	$BinCmd .= runMetaDecoder("",$BinDir,$nodeSpTmpD2, $smplIDs1,$metaGassembly,$MB2coresL,\@paths);
} elsif ($DoMetaBat2 == 4){#GenomeFace
	# Use $BinDir for depth file and tmpDir: GPU job runs on a different node with no access to node-local SSD
	$preCmd .= jgi_depth_cmd(\@paths,"$BinDir/depth$cAssGrp",95,$MB2coresL,$metaGassembly);
	$BinCmd .= runGenomeFace("$BinDir/depth$cAssGrp.jgi.depth.txt",$BinDir,$smplIDs1,$metaGassembly,$MB2coresL,$BinDir);
	$GPUused = 1;
} elsif ($DoMetaBat2 == 5){#SCGBinner
	$BinCmd .= runSCGBinner("",$BinDir,$nodeSpTmpD2,$smplIDs1,$metaGassembly,$MB2coresL,\@paths);
}

my $stone = "$BinDir/Binning.stone";
$BinCmd .= "\n echo $smplIDs1 > $stone\n";


print "running: $preCmd.$BinCmd\n";
#die ;
systemW $preCmd;

#submit this via GPU queue if GenomeFace, otherwise run directly
if ($GPUused){
	my $submitJob = 1;
	my $QSBoptHR = emptyQsubOpt(1,"");
	$QSBoptHR->{useGPUQueue} = 1;
	$QSBoptHR->{gpuCount} = 1;
	my @jobDeps = ();
	my ($jID, $blank) = qsubSystem("$logDir/$BinnerName.GPU.sh", $BinCmd, $submitJob,"25G","$BinnerName.$smplIDs1","","",1,[],$QSBoptHR);
	push(@jobDeps,$jID);
	qsubSystemJobAlive(\@jobDeps, $QSBoptHR, 1);
	$GPUused = 0;
} else {
	systemW $BinCmd;
}


print "Done executing binner $BinnerName\n";
