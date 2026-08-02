#!/usr/bin/env perl
#uses a multi sample assembly to bin contigs with metabat
# perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec//MGS.pl /g/scb/bork/hildebra/SNP/GCs/alienGC2/ /local/bork/hildebra/MB2test/ /g/scb/bork/hildebra/SNP/GCs/alienGC2/Canopy2/clusters.txt
# perl /g/bork3/home/hildebra/dev/Perl/reAssemble2Spec//helpers/MGS/compoundBinning.pl /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/ /local/bork/hildebra/MB2test/ /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/Canopy4_AC/clusters.txt
#perl /hpc-home/hildebra/dev/Perl/MATAF3//MGS.pl -GCd /g/bork3/home/hildebra/data/SNP/GCs/alienGC2/ -tmp /scratch/hildebra//GC/GC_Chicken//MAGs/ -nc 24 -canopies /g/bork3/home/hildebra/data/SNP/GCs/alienGC2/Canopy2/clusters.txt

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long qw( GetOptions );
use File::Path qw(make_path remove_tree);
use File::Basename qw(dirname);
use File::Spec;
use Cwd qw(abs_path);

#.12: checkm2, mem optimizations
#.13: semiBin included
#.14: proGenomes3 added
#.15: GTDBmg added
#.17: process streamlining, HDD usage updates
#.18: updated to "jelly core" MGS clustering, advanced stats on MGS
#.19: set GTDBmg as default, strain resolution updated
#.20: added .LCA MG info to proritize MAGs ("tax clean" MAGs)
#.22: -outD flag. complete rework of clusterMAGs.pl script (single processing)
#.23: clusterMAG binary
#.24: no-canopy fix
#.25: 9.12.23: added extraction of representative MGS genome in contigs (highest qual MAG)
#.26: 13.8.24: added -genomesPerFamily flag & function
#.27: 12.11.24: removed necessity for -canopies flag
#.28: 28.12.25: multi scaling flags for strainScr1 added
#.31: propagate catalog identity and store validated checkpoint manifests
#.32: remove deprecated hierarchical/deep-correlation post-clustering path
#.33: harden sparse Canopy/MAG/MGS outcomes, singleton output, weighted resumes, and phylogeny skips
#.34: restore the documented GTDB default and fingerprint workflow-defining checkpoint options
#.35: invalidate stale input-derived products and validate sparse/downstream outputs before resuming
#.36: consolidate repetitive MAG diagnostics and make stage progress easier to scan
#.37: launch the between-MGS tree immediately after Stage I and defer only its
#     abundance-dependent visualization.
#.38: exclude samples marked SMPL.empty before resolving assembly-group paths
#.39: stream large MGS inputs/outputs and record per-job wall time and peak RSS
#.40: announce startup before configuration and input metadata loading
#.41: use the clusterMAGs binary directly unless the Perl compatibility path is requested
#.42: consume the clusterMAGs binary's compressed MAG report without recompressing it
#.43: precompute catalogue-validated mosaic loci and consolidated outgroups
#.44: remove inline timing wrappers; rely on scheduler sacct accounting
#.45: standardize MAGvsGC.txt.gz in the Bin_<binner> directory
#.46: use the catalog map manifest/identity and isolate binner-specific annotations
#.47: build the between-MGS tree from the selected FMG or predefined GTDB markers
#.48: keep mosaic catalogues, diagnostics, and logs in the binner-local mosaic directory
#.49: bulk-align only genes with mosaic or outgroup comparison potential
#.50: allow strain analysis to run without mosaic preprocessing

my $MGSpipelineVersion = 0.51;
my $clusterID = 95;
my %checkpointParameters;

use Mods::IO_Tamoc_progs qw(getProgPaths jgi_depth_cmd);
use Mods::GenoMetaAss qw(readMap getDirsPerAssmblGrp unzipFileARezip getAssemblPath systemW gzipopen);
use Mods::Subm qw(qsubSystem emptyQsubOpt qsubSystemJobAlive);
use Mods::TamocFunc qw(checkMF);
use Mods::geneCat qw(readMG_LCA);
use Mods::Binning qw (getBinSubdirName createBin2 createBinCtgs runMetaBat runCheckM runCheckM2 createBinFAA readMGS MB2assignedBinIds);
use Mods::Checkpoint qw(write_checkpoint checkpoint_valid);
use Mods::CatalogPaths qw(catalog_identity resolve_catalog_maps);

sub getGoodMBstats;
sub printL;
sub CanopyPrep;
sub invertIndex;
sub _representative_contig_outputs_valid;

$| = 1;
print "Starting MGS pipeline v$MGSpipelineVersion; parsing configuration before loading inputs.\n";

#my $metab2Bin = getProgPaths("metabat2");

my $rareBin = getProgPaths("rare");


#add this? https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4748697/figure/fig-1/
my $inD = "";#$ARGV[0];
my $outD = "";

#$inD.="/" unless($inD =~ m/\/$/);
my $doSubmit = 1;
my $numCore = 4;
my $canCore = 12;
my $memG = 150;#used only for binner
my $legacyV = 0;#legacy (pre Dec `22) parameters
my $perlClusterMAGs = 0; #compatibility/debug implementation; binary is the default
my $rewrTAX = 0;
my $rewrClusterMAGs = 0; #redo clusterMAGs analysis
my $doStrains = 0;
my $prepareMosaicLoci = 1;
my $tmpD = ""; 
my $canopyF = "";
my $nodeTmpD = getProgPaths("nodeTmpDir");
my $useCheckM2 = 0; my $useCheckM1 = 1;
my $binSpeciesMG = 2; #defaults to semibin
my $ignoIncomplMAGs = 1; #if 1, will ignore if no Bin in MF3
my $useGTDBmg = "GTDB";
my $wait4stone = "";
my $wait4stoneTimeout = 86400;
my $useWeightedMGSscores = 1;
my $doBinCtgsPerFam =0 ; #extract bins for each family (or assembly_grp/sample if family missing)
my $customBinFile = ""; #overrides default per GC global behaviour


my $stopAfterCluster = 0; #DEBUG flag!!

#options to pipeline..
GetOptions(
	"GCd=s"      => \$inD, 				#gene catalog dir
	"clusterID=i" => \$clusterID,			#gene-catalog clustering identity percentage
	"outD=s"     => \$outD,				#defaults to $inD/Bin_SB/
	"tmp=s" => \$tmpD,					#temp dir
	"submit=i" => \$doSubmit,			#1:submit jobs, 0: dry run. Default: 1
	"canopies=s" => \$canopyF,			#location of canopy clustering output file (clusters.txt)
	"smallCores=i" => \$numCore,		#cores used for normal jobs (not intensive)
	"bottleneckCores=i" => \$canCore,	#cores for compute intensive jobs
	"redoCluster=i" => \$rewrClusterMAGs,
	"redoTax=i" => \$rewrTAX,			#rewrite tax annotations
	"MGset=s" => \$useGTDBmg,			#GTDB or FMG, which marker genes are used? Default: GTDB
	"wait4stone=s" => \$wait4stone,     #wait for these files to be created, refers currently exclusively to eggNOG annotations that are needed later
	"wait4stoneTimeout=i" => \$wait4stoneTimeout, #maximum wait in seconds; 0 retains unlimited waiting
	"mem=i" => \$memG,					#memory used for intensive jobs
	"strains=i" => \$doStrains,			#1: calc instra species strain phylogenies. Default: 0
	"prepareMosaicLoci=i" => \$prepareMosaicLoci, #1: confirm mosaic loci/outgroups before strain analysis; 0: keep seed clusters separate
	"useCheckM2=i" => \$useCheckM2,		#CheckM2 default qual checking of MAGs/MGS
	"useCheckM1=i" => \$useCheckM1,		#CheckM default qual checking of MAGs/MGS
	"binSpeciesMG=i" => \$binSpeciesMG,	#0=no, 1=metaBat2, 2=SemiBin, 3: MetaDecoder, 4 ,5
	"ignoreIncompleteMAGs=i" => \$ignoIncomplMAGs,	#1: assemblies without MAG calculations are ignored. Default: 1
	"legacy=i" => \$legacyV,			#1: use legacy code as pre Dec `22 (clustering is a bit more muddy, reported abundances slightly different, remember to use -MGset FMG). No longer supported. Default: 0
	"perlClusterMAGs!" => \$perlClusterMAGs,	#explicit compatibility/debug path; default uses the clusterMAGs binary directly
	"genomesPerFamily=i" => \$doBinCtgsPerFam,
	#"customBinFile=s" => \$customBinFile, #take all genes per bin from this file -> need to go to withinStr.pl
) or die "Invalid MGS.pl options\n";
die "Unexpected positional arguments: @ARGV\n" if @ARGV;

#check all correct
checkMF();
die "Select exactly one quality checker: -useCheckM1 1 -useCheckM2 0, or vice versa\n"
	unless ($useCheckM1 ? 1 : 0) + ($useCheckM2 ? 1 : 0) == 1;
die "-wait4stoneTimeout must be zero or a positive number of seconds\n" if $wait4stoneTimeout < 0;
die "Core and memory requests must be positive\n" unless $numCore > 0 && $canCore > 0 && $memG > 0;
die "-binSpeciesMG must be one of 1..5\n" unless $binSpeciesMG >= 1 && $binSpeciesMG <= 5;
die "-clusterID must be between 1 and 100\n" unless $clusterID >= 1 && $clusterID <= 100;
die "-prepareMosaicLoci must be 0 or 1\n"
	unless $prepareMosaicLoci == 0 || $prepareMosaicLoci == 1;

#die "$useCheckM2 $useCheckM1\n";

die "-MGset option has to be \"GTDB\" or \"FMG\"\n" unless ($useGTDBmg eq "GTDB" || $useGTDBmg eq "FMG");

%checkpointParameters = (
	cluster_id             => $clusterID,
	pipeline_version       => $MGSpipelineVersion,
	marker_set             => $useGTDBmg,
	binner                 => $binSpeciesMG,
	quality_checker        => $useCheckM2 ? 'checkm2' : 'checkm1',
	legacy_mode            => $legacyV ? 1 : 0,
	perl_cluster_mags      => $perlClusterMAGs ? 1 : 0,
	weighted_mgs_scores    => $useWeightedMGSscores ? 1 : 0,
	ignore_incomplete_mags => $ignoIncomplMAGs ? 1 : 0,
	genomes_per_family     => $doBinCtgsPerFam ? 1 : 0,
	canopy_file            => $canopyF,
);

die "Needs input dir arg (-GCd)!" if ($inD eq "");
die "Gene catalog directory does not exist: $inD\n" unless -d $inD;
die "-legacy is no longer supported\n" if $legacyV;
$inD = abs_path($inD);
$inD .= "/" unless $inD =~ m{/$};
my $catalogIdentity = catalog_identity($inD);
#die "$doStrains\n";
#set up basic structures
$tmpD = $inD."/tmp/" if ($tmpD eq "");
$tmpD = File::Spec->rel2abs($tmpD);
$tmpD .= "/" unless $tmpD =~ m{/$};
$tmpD .= "$catalogIdentity/";
my $QSBoptHR = emptyQsubOpt($doSubmit,"");
my %QSBopt = %{$QSBoptHR};

my $singleSample = 0;

die "No MGS set: -binSpeciesMG 0 \n" if ($binSpeciesMG == 0);
my $BinnerShrt = getBinSubdirName($binSpeciesMG);

#my $BinnerShrt = "MB2";
#if ($binSpeciesMG == 2){$BinnerShrt = "SB";}#SemiBin
#if ($binSpeciesMG == 3){$BinnerShrt = "MD";}
my $COGdir = "FMG";
if ($useGTDBmg eq "GTDB"){ 
	$COGdir = "GTDBmg";
}

$outD = $inD."/Bin_$BinnerShrt/" if ($outD eq "");
$outD = File::Spec->rel2abs($outD);
$outD .= "/" unless ($outD =~ m/\/$/);
my $logDir = $outD."LOGandSUB/";
my $annoDir = $outD."Annotation/";
#checkpoints
my $chkpDir = "$logDir/checkpoints/"; 
my $checkpointWriter = dirname(abs_path($0))."/../helpers/writeCheckpoint.pl";
my $ABmgsSton = "$chkpDir/abund.mgs.stone";
my $ABmgsSton2 = "$chkpDir/abund.mgs_core.stone";
my $st1ston = "$chkpDir/Stage1.stone";
my $iniMB2sto = "$chkpDir/$BinnerShrt.cm.stone";
my $GTDBtaxSto = "$chkpDir/GTDBTK.stone";
my $BinExtrSto = "$chkpDir/BinExtr.stone";
my $noMGSSto = "$chkpDir/no-usable-mgs.stone";
my $noMGSReport = "$outD/NO_MGS.txt";

#main guide files for MGS
my $finalClusters2 = "$outD/$BinnerShrt.clusters";
my $finalClustersW = "$outD/$BinnerShrt.Wclusters"; #unweighted verssion..
my $finalClustersFilt = $finalClusters2.".core";


if (-e "$inD/LOGandSUB/inmap.txt"){ #this is the outdir of a whole MATAFILER run, or geneCat, doesn't matter
	$singleSample = 0;
	print "Compound Assembly MetaBatting..\n";
} 

make_path($tmpD, $outD, $logDir, $annoDir, $chkpDir);



my $GCd = $inD;
my $mapF = resolve_catalog_maps($GCd);
$checkpointParameters{catalog_identity} = $catalogIdentity;

open LOG, '>', "$logDir/pipeline.log" or die "Cannot open $logDir/pipeline.log: $!\n";
printL "=====================================================\n";
printL "MGS pipeline v$MGSpipelineVersion\n";
printL "Mode: " . ($doSubmit ? "submit" : "dry run") . "; scheduler: $QSBopt{qmode}\n";
printL "Legacy parameter mode: " . ($legacyV ? "yes" : "no") . "\n";
printL "Inputs: gene catalog=$inD; map=$mapF\n";
printL "Paths: output=$outD; temporary=$tmpD; log=$logDir/pipeline.log\n";
printL "Clustering: binner=$BinnerShrt; marker set=$useGTDBmg; identity=$clusterID%; weighted scores="
	. ($useWeightedMGSscores ? "yes" : "no") . "\n";
printL "MAG clustering engine: "
	. ($perlClusterMAGs ? "Perl compatibility path (explicitly enabled)" : "clusterMAGs binary") . "\n";
printL "Quality: " . ($useCheckM2 ? "CheckM2" : "CheckM1")
	. "; ignore incomplete MAGs=" . ($ignoIncomplMAGs ? "yes" : "no") . "\n";
printL "Resources: standard cores=$numCore; bottleneck cores=$canCore; binner memory=${memG}G\n";
printL "Optional analyses: strains=" . ($doStrains ? "yes" : "no")
	. "; mosaic checks=" . ($prepareMosaicLoci ? "yes" : "no")
	. "; family genomes=" . ($doBinCtgsPerFam ? "yes" : "no") . "\n";
printL "Requested rebuilds: clustering=" . ($rewrClusterMAGs ? "yes" : "no")
	. "; taxonomy=" . ($rewrTAX ? "yes" : "no") . "\n";
printL "Requested Canopy assignments: $canopyF\n" if $canopyF ne "";
printL "Configuration accepted; loading mapping and catalogue metadata...\n";

# Fast provenance fingerprints for the primary biological inputs.  Size plus
# mtime avoids hashing very large catalogues on every resume while detecting
# normal replacements and rewrites.
my @checkpointInputs = (
	['matrix',      "$GCd/Matrix.mat.gz"],
	['catalog_fna', "$GCd/compl.incompl.$clusterID.fna"],
	['catalog_faa', "$GCd/compl.incompl.$clusterID.prot.faa"],
);
my $mapIndex = 0;
push @checkpointInputs, ["map_".$mapIndex++, $_] for grep { length } split /,/, $mapF;
push @checkpointInputs, ['canopy', $canopyF] if length($canopyF) && -e $canopyF;
for my $input (@checkpointInputs) {
	my ($label, $path) = @$input;
	my @stat = stat($path);
	die "Cannot fingerprint required MGS input $path\n" unless @stat;
	$checkpointParameters{"${label}_size"} = 0 + $stat[7];
	$checkpointParameters{"${label}_mtime"} = 0 + $stat[9];
}

#die "$mapF\n";
#figure out which compound assemblies there are..

#infer Assembly dirs & corrsponding bams with several Samples (compound assemblies)
my ($hrD,$hrM) = getDirsPerAssmblGrp($mapF);
my %map = %{$hrM};
my %DOs = %{$hrD};
my $rawNumSamples = scalar(@{$map{opt}{smpl_order}});
my @emptySamples = _exclude_empty_samples(\%DOs, \%map);
$hrM = \%map;
my @DoosD = sort keys %DOs; #dirs of assembly groups



my $numSamples = scalar(  @{$map{opt}{smpl_order}}  );#@DoosD;
die "No non-empty samples were found in the mapping input: $mapF\n" unless $numSamples;
$checkpointParameters{empty_samples} = join(',', @emptySamples);
my $profileSamples = _matrix_sample_count("$GCd/Matrix.mat.gz");
my $useCanopies=1;
if ($profileSamples<10 || $canopyF eq ""){$useCanopies=0;}
if ($useCanopies && !-s $canopyF) {
	warn "Requested Canopy assignments are missing or empty; continuing without Canopy MGS: $canopyF\n";
	$useCanopies = 0;
}

my @existingClusterProducts = grep { -e $_ } (
	glob("$outD/$BinnerShrt.clusters*"),
	glob("$outD/$BinnerShrt.Wclusters*"),
);
my $stage1ProvenanceInvalid =
	@existingClusterProducts && !_checkpoint_valid($st1ston);
warn "Existing MGS clustering does not match the current inputs/options; invalidating it before reclustering\n"
	if $stage1ProvenanceInvalid && !$rewrClusterMAGs;

printL "Input metadata loaded successfully.\n";
printL "Canopy assignments: $canopyF\n" if ($canopyF ne "" && $useCanopies);
if (!$useCanopies){
	my $reason = $profileSamples < 10 ? "N<10 matrix samples (N=$profileSamples)" : "no usable Canopy assignment file";
	printL "No Canopies used: $reason\n";
}
if (@emptySamples) {
	printL "Excluded " . scalar(@emptySamples)
		. " sample(s) marked SMPL.empty: " . join(", ", @emptySamples) . "\n";
}
printL "Samples in map: $rawNumSamples; eligible non-empty samples: $numSamples; "
	. "abundance profiles in matrix: $profileSamples\n";
printL "=====================================================\n";
my $cmSuffix = ".cm"; $cmSuffix = ".cm2" if ($useCheckM2); 

#clean up
if ($rewrClusterMAGs || $stage1ProvenanceInvalid) {
	for my $file (glob("$outD/$BinnerShrt.clusters*"), glob("$outD/$BinnerShrt.Wclusters*")) {
		unlink $file or die "Cannot remove $file: $!\n" if -f $file;
	}
	for my $checkpoint ($st1ston, $GTDBtaxSto, $BinExtrSto, $ABmgsSton, $ABmgsSton2, $noMGSSto) {
		unlink $checkpoint or die "Cannot invalidate downstream checkpoint $checkpoint: $!\n"
			if -e $checkpoint;
	}
	unlink $noMGSReport or die "Cannot remove stale $noMGSReport: $!\n" if -e $noMGSReport;
	for my $phyloDir ("$outD/between_phylo", "$outD/within_phylo") {
		remove_tree($phyloDir) if -d $phyloDir;
	}
}
my $ph1flag =
	(-s "$outD/$BinnerShrt.clusters.obs" && -s $finalClusters2 && _checkpoint_valid($st1ston))
	? 0 : 1;
#my $FMGsubs = `wc -l $GCd/Matrix.$COGdir.mat | cut -f1 -d' '`; chomp $FMGsubs; $FMGsubs = int($FMGsubs);
# a whole lot faster.. but imprecise!
my @marker_lca_files = glob("$GCd/$COGdir/*.LCA");
die "No marker-gene LCA files found in $GCd/$COGdir\n" unless @marker_lca_files;
my $FMGsubs = _count_lines_up_to(20, @marker_lca_files);
if ($FMGsubs < 20) {
	_finish_without_mgs("only $FMGsubs marker-gene LCA assignments were available for $profileSamples abundance profile(s)", $noMGSReport, $noMGSSto)
		if $profileSamples < 10;
	die "$GCd/$COGdir/*.LCA suspiciously small (N=$FMGsubs)\nPlease ensure correctness\n";
} #$GCd/Matrix.$COGdir.mat
#die;
my $usableCanopyCount = 0;
if ($useCanopies) {
	my $CanoDir = $canopyF;$CanoDir=~s/\/[^\/]+$/\//; $CanoDir .= "Bins/";
	#die "$CanoDir";
	$usableCanopyCount = CanopyPrep($canopyF,$CanoDir);
	if ($usableCanopyCount) {
		$canopyF .= ".filt";
	} else {
		warn "No Canopy MGS passed size/quality preparation; continuing with MAGs only\n";
		$useCanopies = 0;
	}
}

#run metabat on each assembly group
my $cnt=0; my @jobs;
printL "Found ".scalar(@DoosD) ." assembly groups, ";
if ($ph1flag){
	printL "clustering available binnings\n";
} else {
	printL "reusing existing $BinnerShrt MGS clustering\n";
}
foreach my $Doo (@DoosD){ #this loops ensures Binner predictions exist for each assembly
	#print "$Doo\n";
	last; #should be done in MATAFILER.. deactivate here..
	last if (!$ph1flag && _checkpoint_valid($iniMB2sto));
	my $bef = "";
	my $tmpD2 = "$tmpD$Doo/";
	my $nodeTmpD2 = "$nodeTmpD/checkM/C$Doo/";
	#print "$nodeTmpD2\n";
	$bef .= "mkdir -p $tmpD2\n";# unless (-d $tmpD2);
	#my $allPaths = $DOs{$Doo}{wrdir};
	#my $smplIDtmp = $DOs{$Doo}{SmplID};
	my @smplIDs = @{$DOs{$Doo}{SmplID}};#split /,/,$smplIDtmp;
	my @paths = @{$DOs{$Doo}{wrdir}};#split /,/,$allPaths;
	#next if (@paths <=1);
	my $metaGD = getAssemblPath($paths[-1]);
	my $refFA = $metaGD."/scaffolds.fasta.filt";
	
	
	my $MBout = "$metaGD/Binning/$BinnerShrt/$smplIDs[-1]";
	my $postCmd = "";

	my $CM1done = 0; my $CM2done = 0; my $eBinAssStat=0;
	$CM1done = 1 if (-e "$MBout.cm" );$CM2done = 1 if (-e "$MBout.cm2" );
	$eBinAssStat =1 if (-e "$MBout.assStat");
	#print $MBout." $useCheckM1 $CM1done $useCheckM2 $CM2done $eBinAssStat\n";
	next if ($eBinAssStat && ( ( $useCheckM1 && $CM1done) || ($useCheckM2 && $CM2done ) ) );
	next if ($ignoIncomplMAGs);

	
	#die "$refFA\n";
	#my $refFA = "$inD/$Doo/metag/scaffolds.fasta.filt";
	my $MBcmd = "";
	if ($binSpeciesMG == 1){
		$bef .= jgi_depth_cmd(\@paths,$tmpD2."/depth",95,$numCore,$refFA);# unless (-e );
		$MBcmd = runMetaBat("$tmpD2/depth.jgi.depth.txt",$metaGD."/Binning/$BinnerShrt/",$smplIDs[-1],$refFA);
	} elsif ($binSpeciesMG == 2){
		die "MGS.pl::SemiBin not implemented\n";
	} elsif ($binSpeciesMG == 3){
		die "MGS.pl::MetaDecoder not implemented\n";
	} else {
		die "Binning option $binSpeciesMG not implemented!\n";
	}
	#print $bef.$MBcmd;
	$bef = "" if ($MBcmd eq "");
	my $jobName = "Bin$cnt";
	my $mb2Qual = getProgPaths("mb2qualCheck_scr");
	$postCmd = "\n\nrm -rf $nodeTmpD2; mkdir -p $nodeTmpD2;\n$mb2Qual $refFA $MBout $nodeTmpD2 $numCore $useCheckM2\n" unless (-e "$MBout$cmSuffix"  && -e "$MBout.assStat");
	if ($MBcmd eq "" && $postCmd eq "") {next;}#print "next "; next;}
	#next;
	#die "$postCmd\n";
	$postCmd .= "rm -rf $tmpD2\n";
	#print "$MBout\n";
	#die "$bef$MBcmd$postCmd";
	my ($jobName2, $tmpCmd) = qsubSystem(
		$paths[-1]."LOGandSUB/${BinnerShrt}_bin.sh",
		$bef.$MBcmd.$postCmd,
		$numCore,int($memG)."G",$jobName,"","",1,[],\%QSBopt,
	);
	$cnt++;
	push (@jobs, $jobName2);
	#die $paths[-1]."LOGandSUB/MB2_bin.sh";
}
qsubSystemJobAlive( \@jobs,\%QSBopt );


#check that really all cm 's are there
$cnt=0; my @missedMAGs=(); my $usableMAGcount=0;
printL "Checking $BinnerShrt MAG availability across " . scalar(@DoosD) . " assembly groups\n";
foreach my $Doo (@DoosD){
	my @paths = @{$DOs{$Doo}{wrdir}};#split /,/,$allPaths;
	my @smplIDs = @{$DOs{$Doo}{SmplID}};#split /,/,$smplIDtmp;
	my $metaGD = getAssemblPath($paths[-1]);
	my $MBout = "$metaGD/Binning/$BinnerShrt/$smplIDs[-1]";
	if ($ignoIncomplMAGs && (!-e $MBout || -s $MBout ==0 ) ){
		push (@missedMAGs, $smplIDs[-1]);
		next; 
	}
	#check if maybe emtpy Bin?
	if (!-e "$MBout$cmSuffix" ){
		_touch_empty_file("$MBout$cmSuffix") if _bin_assignments_are_empty($MBout);
	}
	die "Bin $MBout seems incomplete\n" unless (-e "$MBout$cmSuffix" && -e "$MBout.assStat");
	my ($bin_assignments, $bin_quality) = MB2assignedBinIds($MBout, "$MBout$cmSuffix");
	for my $bin (keys %{$bin_assignments}) {
		next unless exists $bin_quality->{$bin};
		$usableMAGcount++
			if $bin_quality->{$bin}{compl} >= 60 && $bin_quality->{$bin}{conta} <= 10;
	}
	$cnt++;
}
printL "MAG availability summary: $cnt usable assembly group(s), "
	. scalar(@missedMAGs) . " missing/empty, $usableMAGcount bin(s) passed the 60% completeness/10% contamination screen\n";
if (@missedMAGs) {
	my @examples = @missedMAGs > 5 ? @missedMAGs[0 .. 4] : @missedMAGs;
	printL "Missing/empty MAG examples: " . join(", ", @examples) . "\n";
	printL "No more missing/empty MAG examples are shown here; the complete list is retained in pipeline.log\n"
		if @missedMAGs > 5;
	print LOG "All missing/empty MAG assembly groups: " . join(", ", @missedMAGs) . "\n";
}
unlink $iniMB2sto or die "Cannot invalidate stale checkpoint $iniMB2sto: $!\n"
	if @missedMAGs && -e $iniMB2sto;
_touch_checkpoint($iniMB2sto, 'per-sample-mag-quality') unless _checkpoint_valid($iniMB2sto) || @missedMAGs;

_finish_without_mgs("no assigned MAG passed the minimum 60% completeness/10% contamination screen and no usable Canopy MGS was available", $noMGSReport, $noMGSSto)
	unless $usableMAGcount || $usableCanopyCount;

#if ($cnt){	print "Waiting for jobs to finish.. restart when done\n";	exit(0);}

#----------------------------------------------------------------------------------------------------
#from here merging of MAGs into MGS

#just writes MAG.$BinnerShrt.assStat.summary, important for  reading in %valMBs
getGoodMBstats() if (!-e $finalClusters2 );#die;


#cluster MAGs based on shared genes between them
@jobs = ();
if ($ph1flag  || !-e "$outD/$BinnerShrt.clusters" ){
	my $cmd;
	if ($perlClusterMAGs) {
		my $clusscr = getProgPaths("clusterMGS_scr");
		my $canoIncl = $useCanopies ? "-canopies $canopyF" : "";
		$cmd = "$clusscr -GCd $GCd -BinDir $outD -logDir $logDir -binSpeciesMG $binSpeciesMG -MGset $useGTDBmg -clusterID $clusterID -useCheckM1 $useCheckM1 -useCheckM2 $useCheckM2 -legacy $legacyV -perlClusterMAGs -cores 1 $canoIncl ";
		$cmd .= "1>&2 > $logDir/clusterMGS_scr.log\n";
		warn "Using the explicitly requested Perl clusterMAGs compatibility path; this is slower and intended only for debugging or result comparison\n";
		printL "Clustering MAGs into MGS with the Perl compatibility implementation; detailed output: $logDir/clusterMGS_scr.log\n";
	} else {
		my $clusterBinary = getProgPaths("clusterMAGs");
		my $clusteringMapF = _maps_without_empty_samples(
			$mapF, \@emptySamples, "$logDir/nonempty_maps",
		);
		my $canoIncl = $useCanopies ? "-canopyDir $canopyF" : "";
		$cmd = "$clusterBinary -CMsuffix $cmSuffix -path2Bins Binning/$BinnerShrt/ -FILEtag $BinnerShrt -MGStag MGS. -geneCatIdx $GCd/compl.incompl.$clusterID.fna.clstr.idx -log $logDir/clusterMGS_scr.log -LCAdir $GCd/$COGdir -outDir $outD -map $clusteringMapF $canoIncl -MGfile $GCd/$COGdir.subset.cats \n";
		$cmd .= "test -s $outD/MAGvsGC.txt.gz\n";
		printL "Clustering MAGs into MGS directly with the clusterMAGs binary; detailed output: $logDir/clusterMGS_scr.log\n";
	}
	systemW $cmd;
}

die if ($stopAfterCluster); #DEBUGing only!!
qsubSystemJobAlive( \@jobs,\%QSBopt ) if (@jobs);

#decide between weighted and unweighted scores for binning
my $activeMGSCount = _mgs_count($finalClusters2);
my $weightedMGSCount = _mgs_count($finalClustersW);
my $preservedMGSCount = _mgs_count("${finalClusters2}UW");
my $activatedOnlyWeighted = 0;
if (!$activeMGSCount && $preservedMGSCount) {
	warn "Restoring preserved unweighted MGS assignments after an incomplete weighted handoff\n";
	rename "${finalClusters2}UW", $finalClusters2
		or die "Cannot restore ${finalClusters2}UW as $finalClusters2: $!\n";
	$activeMGSCount = $preservedMGSCount;
}
if (!$activeMGSCount && !$preservedMGSCount && $weightedMGSCount) {
	warn "Activating the only available weighted MGS assignments\n";
	rename $finalClustersW, $finalClusters2
		or die "Cannot activate $finalClustersW as $finalClusters2: $!\n";
	$activeMGSCount = $weightedMGSCount;
	$weightedMGSCount = 0;
	$activatedOnlyWeighted = 1;
}
if ($useWeightedMGSscores && !$preservedMGSCount && !$activatedOnlyWeighted){
	unlink "${finalClusters2}UW" or die "Cannot remove empty ${finalClusters2}UW: $!\n"
		if -e "${finalClusters2}UW";
	if ($activeMGSCount && $weightedMGSCount) {
		rename $finalClusters2, "${finalClusters2}UW" or die "Cannot preserve $finalClusters2: $!\n";
		rename $finalClustersW, $finalClusters2 or die "Cannot activate weighted clusters $finalClustersW: $!\n";
		$activeMGSCount = $weightedMGSCount;
	} elsif ($activeMGSCount) {
		warn "Weighted MGS assignments were not produced; retaining the valid unweighted assignments\n";
	}
	#system "touch $finalClustersW.mov";
}

_finish_without_mgs("MAG/Canopy clustering produced no MGS assignments", $noMGSReport, $noMGSSto)
	unless $activeMGSCount;

my $observation_file = "$outD/$BinnerShrt.clusters.obs";
if (!-s $observation_file) {
	if ($activeMGSCount == 1) {
		_write_single_mgs_observations($finalClusters2, $observation_file);
		printL "Synthesized the missing observation table for a single MGS\n";
	} else {
		die "MGS observation file is missing for $activeMGSCount MGS: $observation_file\n";
	}
}

#alt motulizer?
#system "motulizer"
#$hr = readmotulizertable();

#create a core of metabat2 clusters, based on gene occurrence (between different samples)
my $RfilterMB2 = getProgPaths("filterMB2core");
my $postCmd = "$RfilterMB2 $finalClusters2\n" ; #creates $outD/MB2.clusters.core & $outD/MB2.clusters.ext
die "MGS cluster file is missing or empty: $finalClusters2\n" unless -s $finalClusters2;
systemW $postCmd if !-s $finalClustersFilt;
my $coreMGSCount = _mgs_count($finalClustersFilt);
_finish_without_mgs("no MGS retained any core genes after post-filtering", $noMGSReport, $noMGSSto)
	unless $coreMGSCount;
unlink $noMGSReport or die "Cannot remove stale $noMGSReport: $!\n" if -e $noMGSReport;
unlink $noMGSSto or die "Cannot remove stale $noMGSSto: $!\n" if -e $noMGSSto;

#die;

printL "---------------------------------------------------------------------\n";
printL "Stage I clustering done, MGS calculated.\nProgressing to Stage II: annotations, phylogenies and abundances\n";
printL "Using $finalClustersFilt as MGS rep\n";
printL "---------------------------------------------------------------------\n";
_touch_checkpoint($st1ston, 'stage-1', $finalClustersFilt) unless _checkpoint_valid($st1ston);

# Start the between-MGS tree as soon as the newly published MGS core set exists.
# Tree inference depends only on the selected marker proteins and MGS membership.  Its
# abundance-annotated visualization is submitted later, after MGS abundance is
# available, so taxonomy and abundance can run concurrently with the tree.
my $treeMem = "120";
if ($numSamples > 2000){ #scale with the number of assembly groups
	$treeMem = "200";
}
my $phyloBetween = getProgPaths("MGSPhyloBetween_scr");
my $baseTreeCmd = "$phyloBetween -GCd $GCd -MGS $finalClustersFilt -MGset $useGTDBmg -mem $treeMem -c $canCore -MSAprogram 4 -fast 0 ";
my $wait4tree = 2;
my $outDphylo = "$outD/between_phylo/";
my $iniTree = "$outDphylo/phylo/IQtree_allsites.treefile";
my $treePdf = "$outDphylo/phylo/IQtree_allsites.pdf";
my $ph1Cmd = "$baseTreeCmd -outD $outDphylo -wait2finish $wait4tree -visualize 0 ";
my $treedep = "";
my $betweenTreeSkipped = 0;

if ($coreMGSCount < 3) {
	$betweenTreeSkipped = 1;
	printL "Skipping between-MGS phylogeny: at least 3 MGS are required, but only $coreMGSCount were retained\n";
} elsif (!-s $iniTree) {
	printL "Preparing between-MGS phylogeny immediately after MGS creation in $outDphylo\n";

	if ($useGTDBmg eq "FMG") {
		my $refTreeMsg = "\n################\n# If you want to include custom reference genomes, use\n# $baseTreeCmd -outD $outD/customRefs/ -refGenos [refs]\n################\n";
		print $refTreeMsg;
		$refTreeMsg = "#If you want to include custom reference genomes, use $baseTreeCmd -outD $outD/customRefs/ -refGenos [refs]";
		$ph1Cmd .= " -xtraMsg \"$refTreeMsg\";";
	} else {
		printL "Custom reference genomes are unavailable for predefined GTDB-marker trees; use -MGset FMG for FMG extraction from references\n";
	}

	if (!$doSubmit) {
		print "Dry run: between-MGS launcher was not executed.\n";
	} else {
		my $ph1OUT = `$ph1Cmd`;
		my $ph1Status = $?;
		die "Between-MGS phylogeny command failed (exit " . ($ph1Status >> 8) . "):\n$ph1Cmd\n$ph1OUT\n"
			if $ph1Status != 0;
		if ($ph1OUT =~ m/WAITID=(\d+)/) {
			$treedep = $1;
			printL "Between-MGS phylogeny submitted as job $treedep\n";
		} elsif ($ph1OUT =~ m/^SKIPPED=(.+)$/m) {
			$betweenTreeSkipped = 1;
			printL "Between-MGS phylogeny was skipped by its launcher: $1\n";
		} else {
			die "Between-MGS phylogeny command reported neither WAITID nor SKIPPED:\n$ph1OUT\n";
		}
	}
} else {
	printL "Reusing existing between-MGS phylogeny: $iniTree\n";
}


#get checkM quality for new Bins
my $binD = "$outD/Genomes/MGS_GC/";
my $binDctg = "$outD/Genomes/MGS_ctg/";
my $binDctgFam = "$outD/Genomes/MGS_ctg_fam/";
make_path($binD, $binDctg, $binDctgFam);
my $binExtractionValid = _checkpoint_valid($BinExtrSto);
$binExtractionValid &&= _representative_contig_outputs_valid($binDctg);
$binExtractionValid &&= _representative_contig_outputs_valid($binDctgFam)
	if $doBinCtgsPerFam;
unless ($binExtractionValid) {
	# Prevent removed/renamed MGS from surviving as stale genomes after a
	# clustering or catalogue change.
	remove_tree($_) for grep { -d $_ } ($binD, $binDctg, $binDctgFam);
	make_path($binD, $binDctg, $binDctgFam);
}

#gget repr genomes for each MGS
if (!$binExtractionValid){
	print "Creating reference genome fasta's for all MGS based on gene cat genes in\n$binD\n";
	createBin2($binD,"$finalClustersFilt","$GCd/compl.incompl.$clusterID.prot.faa","faa");
	createBin2($binD,"$finalClustersFilt","$GCd/compl.incompl.$clusterID.fna","fna");
}



#clean up a bit..
if ($rewrTAX) {
	for my $path (
		glob("$annoDir/GTDB*"), glob("$annoDir/kraken2*"), glob("$annoDir/specI*"),
		"$annoDir/${COGdir}_MGS",
		$ABmgsSton, $ABmgsSton2
	) {
		if (-d $path) { remove_tree($path); }
		elsif (-f $path) { unlink $path or die "Cannot remove $path: $!\n"; }
	}
}
my @jobs2wait=();

#basic quality checks are done at this point
#now get 1) taxonomy 2) phylogeny #) abundance matrix of MGS 4)abundance matrix MGS + specI (to capture unbinned species)

#GTDB tax & kraken2 tax
my $GTDBtaxF = "$annoDir/GTDBTK.tax";
if (!-e $GTDBtaxF || !-e"$annoDir/gtdbtk.summary.tsv" || !_checkpoint_valid($GTDBtaxSto)){
	my $GTDBtax = getProgPaths("taxPerMGSgtdb_scr");
	my $memGTDB = 230; 
	$memGTDB = 300;#high mem situation..
	my $cmd = "$GTDBtax $binD $numCore $nodeTmpD/GTDBmgs/ $outD\n";
	$cmd .= "test -s $outD/GTDBTK.tax\ntest -s $outD/gtdbtk.summary.tsv\n";
	$cmd .= "mv $outD/GTDBTK.tax $outD/gtdbtk.summary.tsv $annoDir\n";
	$cmd .= "test -s $GTDBtaxF\ntest -s $annoDir/gtdbtk.summary.tsv\n";
	$cmd .= _checkpoint_command($checkpointWriter, $GTDBtaxSto, 'gtdb-taxonomy',
		$GTDBtaxF, "$annoDir/gtdbtk.summary.tsv", $finalClustersFilt);
	#changed mem from 370 to 100 with GTDB-TK 2.1.0
	my $tmpSHDD = $QSBopt{tmpSpace};	$QSBopt{tmpSpace} = "150G"; 
	my ($jobName2, $tmpCmd) = qsubSystem($logDir."/GTDB.sh",
		$cmd,
		$numCore,int($memGTDB)."G","GTDB_MGS","","",1,[],\%QSBopt);
	$QSBopt{tmpSpace} =$tmpSHDD;
	push(@jobs2wait,$jobName2);
}

if (!$binExtractionValid){
	print "\n\nCreating reference genome fasta's for all MGS based on contigs in\n$binDctg\n\n";
	if ($doBinCtgsPerFam){
	#do I really need per family genomes??
		print "Also creating family-wise ref genomes\n";
		createBinCtgs($binDctgFam,$hrM,"$outD/MAGvsGC.txt.gz",1,$BinnerShrt);
		#die;
	}

	createBinCtgs($binDctg,$hrM,"$outD/MAGvsGC.txt.gz",0,$BinnerShrt);
	my @geneBinFiles = grep { -f $_ } glob("$binD/*");
	my @contigBinFiles = grep { -f $_ } glob("$binDctg/*");
	my @familyBinFiles = $doBinCtgsPerFam ? grep { -f $_ } glob("$binDctgFam/*") : ();
	die "Bin extraction produced no gene-based MGS genomes in $binD\n" unless @geneBinFiles;
	die "Bin extraction produced no representative contig genomes in $binDctg\n" unless @contigBinFiles;
	_touch_checkpoint($BinExtrSto, 'extract-bin-contigs',
		$finalClustersFilt, @geneBinFiles, @contigBinFiles, @familyBinFiles);
}

#wait for checkm/GTDB
qsubSystemJobAlive( \@jobs2wait,\%QSBopt ) if $doSubmit;
if ($doSubmit) {
	die "GTDB taxonomy stage incomplete\n$GTDBtaxF\n$annoDir/gtdbtk.summary.tsv\n$GTDBtaxSto\n"
		unless -s $GTDBtaxF && -s "$annoDir/gtdbtk.summary.tsv" && _checkpoint_valid($GTDBtaxSto);
}
#if (!-e "$finalClusters2$cmSuffix" && -e "$finalClusters3$cmSuffix"){system "mv $finalClusters3$cmSuffix $finalClusters2$cmSuffix";} #needs to be moved
die "Bin-quality scores missing\n$finalClusters2$cmSuffix\n" if (!-e "$finalClusters2$cmSuffix" );
#die "$cmSuffix\n";

#get only MGS at >$complThre compl, <5 contamination (middle qual), to be used in abundance, strains etc
#create filtered down final liste

#generate taxonomy from kraken2 assignments
my @kraken_jobs;
my $krakenInput = "$GCd/Anno/Tax/krak2.txt";
my $krakenSkipped = 0;
if (!-s "$annoDir/kraken2.LCA" || !-s "$annoDir/kraken2.tax"){
	if (!-s $krakenInput) {
		$krakenSkipped = 1;
		warn "Optional Kraken input is missing or empty; skipping MGS Kraken taxonomy:\n$krakenInput\n";
	} else {
		my $kr2taxScr = getProgPaths("taxPerMGS_scr");
		my $cmd =  "$kr2taxScr $finalClustersFilt $GCd $annoDir/kraken2\n";# unless (-e "$finalClusters2.LCA");
		my $tmpSHDD = $QSBopt{tmpSpace};	$QSBopt{tmpSpace} = "0";
		my ($jobName2, $tmpCmd) = qsubSystem($logDir."/krak2MGS.sh",
			$cmd,
			1,int(200/1)."G","KR2_MGS","","",1,[],\%QSBopt) ;
		$QSBopt{tmpSpace} =$tmpSHDD;
		push @kraken_jobs, $jobName2 if $jobName2;
	}
}


#redo tree and abundance:
#system "rm -r $outD/between_phylo/ $GCd//Anno/Tax/SpecI_MGS/ $outD/specI.tax";

#generate abundances per MGS
if (0 && !-e "$finalClusters2.matL0.txt"){ #deprecated, use specI based annotations instead..
	invertIndex($finalClusters2,"$finalClusters2.rev") unless (-e "$finalClusters2.rev");
	my $cmd = "$rareBin sumMat -i $GCd/Matrix.mat.gz -o $finalClusters2.mat -refD $finalClusters2.rev -t $numCore\n";
	$cmd .= "rm $finalClusters2.rev\n";
#	systemW $cmd;
	my $tmpSHDD = $QSBopt{tmpSpace};	$QSBopt{tmpSpace} = "0"; 
	my ($jobName2, $tmpCmd) = qsubSystem($logDir."/MGSabund.sh",
		$cmd,
		1,int(100)."G","AB1_MGS","","",1,[],\%QSBopt) ;
	$QSBopt{tmpSpace} =$tmpSHDD;
}

#die;
#once all tax annotations are done, infer consensus tax for MAGs
#annotate specI's with MAGs added..


my $specIoutDir = "$annoDir/${COGdir}_MGS";
my $specIabundance = "$specIoutDir/specI.mat";
my @annotation_jobs;
unless (_checkpoint_valid($ABmgsSton) && -s $specIabundance && -s "$annoDir/specI.tax"){
	my $specIabu = getProgPaths("specIGC_scr");
	my $cmdSI = "$specIabu -GCd $GCd -cores $canCore -MGS $finalClustersFilt -MGStax $GTDBtaxF -MGset $useGTDBmg -outD $specIoutDir\n";
	if ($legacyV){
		$specIabu = getProgPaths("specIGC_scr_v0");
		$cmdSI = "$specIabu $GCd $canCore $finalClustersFilt $GTDBtaxF\n";
	}
	$cmdSI .= "cp $specIoutDir/MGS2speci.txt $annoDir/specI.tax\n";
	$cmdSI .= "test -s $specIabundance\n";
	$cmdSI .= "test -s $annoDir/specI.tax\n";
	$cmdSI .= _checkpoint_command($checkpointWriter, $ABmgsSton, 'mgs-abundance',
		$specIabundance, "$annoDir/specI.tax", $finalClustersFilt, $GTDBtaxF);
	printL "Merge with SpecI & get abundance \n";
	#print "$cmdSI\n";

	# The selected marker-set LCA inputs were validated before clustering.  Always
	# run this expensive aggregation through the configured submission backend;
	# the previous FMG/tax/*.tmp.m8 shortcut was unrelated to GTDB marker runs and
	# could execute a large job directly on the launcher node.
	my $tmpSHDD = $QSBopt{tmpSpace};	$QSBopt{tmpSpace} = "0";
	my ($jobName2, $tmpCmd) = qsubSystem($logDir."/abundMGS.sh",
		$cmdSI,
		1,int(200/1)."G","AB2_MGS","","",1,[],\%QSBopt) ;
	$QSBopt{tmpSpace} =$tmpSHDD;
	push @annotation_jobs, $jobName2 if $jobName2;
}

qsubSystemJobAlive( \@annotation_jobs,\%QSBopt ) if $doSubmit && @annotation_jobs;
if ($doSubmit) {
	die "MGS/specI annotation stage incomplete\n$specIabundance\n$annoDir/specI.tax\n$ABmgsSton\n"
		unless -s $specIabundance && -s "$annoDir/specI.tax" && _checkpoint_valid($ABmgsSton);
}

my @marker_jobs;
unless (_checkpoint_valid($ABmgsSton2) && -s "$outD/Annotation/Abundance/MGS.matL0.txt" && -s "$outD/Annotation/Abundance/MGS.matL7.txt"){
	my $MMLscr = getProgPaths("MAGMGSLCA_scr");
	my $cmdSI2 = "$MMLscr -GCd $GCd -cores $numCore -MGset $useGTDBmg -Binner $BinnerShrt -binD $outD\n";
	$cmdSI2 .= "test -s $outD/Annotation/Abundance/MGS.matL0.txt\n";
	$cmdSI2 .= "test -s $outD/Annotation/Abundance/MGS.matL7.txt\n";
	$cmdSI2 .= _checkpoint_command($checkpointWriter, $ABmgsSton2, 'marker-mgs-abundance',
		"$outD/Annotation/Abundance/MGS.matL0.txt", "$outD/Annotation/Abundance/MGS.matL7.txt",
		$finalClustersFilt, $GTDBtaxF);
	my $tmpSHDD = $QSBopt{tmpSpace};	$QSBopt{tmpSpace} = "0"; 
	# qsubSystem's memory argument is emitted as total memory by the Slurm
	# backend; keep the intended 100 GiB request independent of thread count.
	my ($jobName2, $tmpCmd) = qsubSystem($logDir."/abundMGS_core.sh",
		$cmdSI2,
		$numCore,"100G","AB_MGS_core","","",1,[],\%QSBopt) ;
	$QSBopt{tmpSpace} =$tmpSHDD;
	push @marker_jobs, $jobName2 if $jobName2;
}

qsubSystemJobAlive( \@kraken_jobs,\%QSBopt ) if $doSubmit && @kraken_jobs;
qsubSystemJobAlive( \@marker_jobs,\%QSBopt ) if $doSubmit && @marker_jobs;
if ($doSubmit) {
	warn "Optional Kraken MGS taxonomy stage incomplete; continuing without Kraken-derived MGS taxonomy\n"
		unless $krakenSkipped || (-s "$annoDir/kraken2.LCA" && -s "$annoDir/kraken2.tax");
	die "Marker-based MGS abundance stage incomplete\n" unless _checkpoint_valid($ABmgsSton2)
		&& -s "$outD/Annotation/Abundance/MGS.matL0.txt"
		&& -s "$outD/Annotation/Abundance/MGS.matL7.txt";
}

# Visualization needs the abundance matrix, unlike tree inference itself.
# Submit it only now, and depend on a newly launched tree when necessary.
if (!$betweenTreeSkipped && !-s $treePdf) {
	my $treeAbundance = "$outD/Annotation/Abundance/MGS.matL7.txt";
	if (!$doSubmit) {
		print "Dry run: between-MGS tree visualization was not submitted.\n";
	} else {
		die "Cannot visualize between-MGS tree without abundance matrix: $treeAbundance\n"
			unless -s $treeAbundance;
		my $vizTree = getProgPaths("vizBtwPhylo_R");
		my $vizCmd = "$vizTree $treeAbundance $iniTree $treePdf\n";
		$vizCmd .= "test -s $treePdf\n";
		my ($vizDep, $vizSubmitCmd) = qsubSystem(
			$logDir."/interMGSphyloViz.sh",
			$vizCmd, 1, "20G",
			"MGSphyloViz", $treedep, "", 1, [], \%QSBopt,
		);
		printL "Between-MGS visualization submitted as job $vizDep"
			. ($treedep ne "" ? " after tree job $treedep" : "") . "\n";
	}
}



#die;



#need to rewrite to new format, also check for passed MGS to include in intra-strain analysis
#needs MGSselection.txt
#from here on relies on files from R script, but this could be also auto generated (>80 compl, <5 cota) <- task done (may 20)
#reformat_4phylo($finalClustersFilt) unless (-e "$finalClustersFilt.mgs");


#process ends here unless strains need to be calculated
if ($doStrains == 0){
	printL "\n\nCompound Binning script finished.. no strain analysis set\n";
	close LOG;
	exit(0);
}
print "\n\n########################\nStarting strain delineation MGS\n########################\n";

if ($wait4stone ne ""){
	my $cntWaits=0;
	my $waitStart = time;
	while ( !-e $wait4stone){
		if ($cntWaits == 0){print "\n\nWaiting for process generating $wait4stone\nIf the process aborted with error, please exit this routine as well\n\n";}
		if ($wait4stoneTimeout > 0 && time - $waitStart >= $wait4stoneTimeout) {
			die "Timed out after $wait4stoneTimeout seconds waiting for $wait4stone. Check the producing job before rerunning.\n";
		}
		$cntWaits++;
		sleep(10);
	}
}
#die "XX\n";

my $strain1scr = getProgPaths("MGS_strain1_scr");
my $memUsage = 30; #in Gb
my $NsubJobs = 0 ; #split job up?
my $preCompCons = 0;

if ($numSamples > 1500){#scale with the number of assembly groups
	$memUsage = 15; $NsubJobs = 30; $preCompCons = 0; 
} elsif ($numSamples > 700){
	$memUsage = 12; $NsubJobs = 10; $preCompCons = 0;#at a certain size it makes more sense to handle  preCompCons this in-job
} elsif ($numSamples > 400){
	$memUsage = 7; $NsubJobs = 4; $preCompCons = 10;
} elsif ($numSamples > 150){ #
	$memUsage = 5; $NsubJobs = 0; $preCompCons = 5;
}
#my $prunTree = "$outD/between_phylo/prunned.nwk";
#
#my $ph2Cmd = "$strain1scr $GCd $finalClustersFilt.mgs $canCore $iniTree 0 1\n";#$outD/between_phylo/phylo/IQtree.treefile\n";
my $mosaicDir = "$outD/mosaic/";
my $mosaicCatalogue = "$mosaicDir/$BinnerShrt.clusters.mosaic_loci.$clusterID.confirmed.tsv";
my $ph2Cmd = "mkdir -p "._shell_quote("$outD/within_phylo/")." || exit 65\n";
$ph2Cmd .= "$strain1scr -GCd $GCd -MGS $finalClustersFilt -mosaicMGS $finalClusters2 -MGSabundance $outD/Annotation/Abundance/MGS.matL7.txt -MGset $useGTDBmg -clusterID $clusterID -maxCores $canCore -rmMSA 1 -preCompConsSNP $preCompCons -selfMemGb $memUsage -mosaicMemGb $memG -onlySubmit 1 -submit $doSubmit -reSubmit 0 -maxSubJob $NsubJobs -redoSubmissionData 0 -outD $outD/within_phylo/ -prepareMosaicLoci $prepareMosaicLoci ";
$ph2Cmd .= "-mosaicLoci $mosaicCatalogue " if $prepareMosaicLoci;
$ph2Cmd .= "-MGSphylo $iniTree " if -s $iniTree || $treedep ne "";
$ph2Cmd .= "\n";

$ph2Cmd .= "#consider adapting further options: \n#-rmMSA 0 -presortGenes 1700 -maxGenes 500 -MGSminGenesPSmpl 5 -multiGeneSmplMax 0.15 -conspGeneSmplMax 0.05 -nodeTmp [path]\n";#$outD/between_phylo/phylo/IQtree.treefile\n";
$ph2Cmd .= "#-minSNPCallQual 20 -GenesPerSpecies 0.1 -GeneLengthMin 0.5 -skipIndels 0 -minSNPDepth 2 -SNPdepthFilterScale 0.1 -SNPindelRangeFilt 5 -SNPadaptiveQual 0.0";
#systemW $ph2Cmd;
my $launcherCores = 1;
my $launcherMemory = $memUsage;
printL "Preparing within-MGS strain analysis in $outD/within_phylo/ "
	. "($launcherCores cores, ${launcherMemory}G launcher memory, "
	. ($prepareMosaicLoci
		? "mosaic checks delegated to a ${canCore}-core prerequisite job"
		: "mosaic checks disabled")
	. ", $NsubJobs subjob partition(s))\n";
my $tmpSHDD = $QSBopt{tmpSpace};	$QSBopt{tmpSpace} = "20"; #needs some tmp space for on the fly creations.. 
my ($jobName2, $tmpCmd) = qsubSystem($logDir."/strainMGS.sh",
	$ph2Cmd,
	$launcherCores,int($launcherMemory)."G","strainKickoff",$treedep,"",1,[],\%QSBopt) ;
$QSBopt{tmpSpace} =$tmpSHDD;

#get phylogenies intra-species.. this requires a lot of power and best called from big cluster..
printL "Compound Binning script finished; within-MGS strain jobs were prepared in $logDir/strainMGS.sh\n";
close LOG;
exit(0);




















#####################################################################
#####################################################################
# Subroutines

sub _read_one_line {
	my ($file) = @_;
	open my $fh, '<', $file or die "Cannot open $file: $!\n";
	my $line = <$fh>;
	close $fh or die "Cannot close $file: $!\n";
	die "Expected content in $file, but it is empty\n" unless defined $line;
	chomp $line;
	return $line;
}

sub _count_lines_up_to {
	my ($limit, @files) = @_;
	die "Line-count limit must be positive\n" unless defined($limit) && $limit > 0;
	my $count = 0;
	for my $file (@files) {
		open my $fh, '<', $file or die "Cannot open $file: $!\n";
		while (<$fh>) {
			$count++;
			last if $count >= $limit;
		}
		close $fh or die "Cannot close $file: $!\n";
		last if $count >= $limit;
	}
	return $count;
}

sub _matrix_sample_count {
	my ($file) = @_;
	my ($fh, $ok) = gzipopen($file, 'gene abundance matrix', 1, 0);
	die "Cannot open gene abundance matrix $file\n" unless $ok && defined $fh;
	my $header = <$fh>;
	# This is deliberately a header-only read. Closing the pipe before pigz has
	# emitted the multi-GB matrix body gives pigz SIGPIPE and therefore a false
	# close() status; that is expected here and must not abort a resumed run.
	close $fh;
	die "Gene abundance matrix is empty: $file\n" unless defined $header;
	chomp $header;
	my @fields = split /\t/, $header, -1;
	pop @fields while @fields && $fields[-1] eq '';
	shift @fields if @fields;
	die "Gene abundance matrix has no sample columns: $file\n" unless @fields;
	return scalar @fields;
}

sub _exclude_empty_samples {
	my ($groups, $map) = @_;
	my %empty;

	for my $sample (@{$map->{opt}{smpl_order} || []}) {
		my $work_dir = $map->{$sample}{wrdir};
		next unless defined($work_dir) && length($work_dir);
		$empty{$sample} = 1 if -e "$work_dir/SMPL.empty";
	}

	for my $group (keys %{$groups}) {
		my @sample_ids = @{$groups->{$group}{SmplID} || []};
		my @work_dirs = @{$groups->{$group}{wrdir} || []};
		die "Assembly group $group has mismatched sample and working-directory lists\n"
			unless @sample_ids == @work_dirs;

		my (@eligible_ids, @eligible_dirs);
		for my $index (0 .. $#sample_ids) {
			if (-e "$work_dirs[$index]/SMPL.empty") {
				$empty{$sample_ids[$index]} = 1;
				next;
			}
			push @eligible_ids, $sample_ids[$index];
			push @eligible_dirs, $work_dirs[$index];
		}

		if (@eligible_ids) {
			$groups->{$group}{SmplID} = \@eligible_ids;
			$groups->{$group}{wrdir} = \@eligible_dirs;
		} else {
			delete $groups->{$group};
		}
	}

	@{$map->{opt}{smpl_order}} =
		grep { !$empty{$_} } @{$map->{opt}{smpl_order} || []};
	delete $map->{$_} for keys %empty;
	return sort keys %empty;
}

sub _maps_without_empty_samples {
	my ($map_files, $empty_samples, $target_dir) = @_;
	die "Mapping files are required for MAG clustering\n"
		unless defined($map_files) && length($map_files);
	die "Empty-sample list must be an array reference\n"
		unless ref($empty_samples) eq 'ARRAY';
	return $map_files unless @{$empty_samples};

	my %empty = map { $_ => 1 } @{$empty_samples};
	make_path($target_dir);
	my @filtered_maps;
	my $map_index = 0;
	for my $input_map (split /,/, $map_files) {
		my $output_map = "$target_dir/map.$map_index.txt";
		my $temporary = "$output_map.tmp.$$";
		open my $input, '<', $input_map or die "Cannot open map $input_map: $!\n";
		open my $output, '>', $temporary or die "Cannot write filtered map $temporary: $!\n";
		while (my $line = <$input>) {
			my ($sample) = split /\t/, $line, 2;
			next if $empty{$sample};
			print {$output} $line or die "Cannot write filtered map $temporary: $!\n";
		}
		close $input or die "Cannot close map $input_map: $!\n";
		close $output or die "Cannot close filtered map $temporary: $!\n";
		unlink $output_map or die "Cannot replace filtered map $output_map: $!\n"
			if -e $output_map;
		rename $temporary, $output_map
			or die "Cannot publish filtered map $output_map: $!\n";
		push @filtered_maps, $output_map;
		$map_index++;
	}
	return join(',', @filtered_maps);
}

sub _touch_checkpoint {
	my ($file, $stage, @outputs) = @_;
	write_checkpoint($file,
		parameters => { %checkpointParameters, stage => $stage },
		outputs => \@outputs,
	);
}

sub _touch_empty_file {
	my ($file) = @_;
	open my $fh, '>', $file or die "Cannot create checkpoint $file: $!\n";
	close $fh or die "Cannot close checkpoint $file: $!\n";
}

sub _checkpoint_valid {
	my ($file) = @_;
	# Empty legacy stones cannot encode marker-set, binner, or QC provenance.
	# Rebuild them in this workflow instead of silently accepting stale state.
	return 0 unless defined($file) && -s $file;
	return checkpoint_valid($file, parameters => \%checkpointParameters);
}

sub _representative_contig_outputs_valid {
	my ($directory) = @_;
	return 0 unless defined($directory) && -d $directory;
	my @outputs = grep { -f $_ } glob("$directory/*");
	return 0 unless @outputs;
	return !grep {
		$_ !~ /\.(?:fa|fna|fasta)\.gz\z/i
	} @outputs;
}

sub _shell_quote {
	my ($value) = @_;
	$value = '' unless defined $value;
	$value =~ s/'/'"'"'/g;
	return "'$value'";
}

sub _checkpoint_command {
	my ($writer, $stone, $stage, @outputs) = @_;
	my %parameters = (%checkpointParameters, stage => $stage);
	my @args = ('perl', $writer, '--stone', $stone);
	for my $key (sort keys %parameters) {
		push @args, ('--param', "$key=$parameters{$key}");
	}
	push @args, map { ('--output', $_) } @outputs;
	return join(' ', map { _shell_quote($_) } @args) . "\n";
}

sub _bin_assignments_are_empty {
	my ($file) = @_;
	open my $fh, '<', $file or die "Cannot open bin assignments $file: $!\n";
	while (my $line = <$fh>) {
		chomp $line;
		next unless length $line;
		my @fields = split /\t/, $line;
		next if @fields >= 2 && $fields[0] eq 'Sequence ID';
		if (defined $fields[1] && $fields[1] ne '0') {
			close $fh;
			return 0;
		}
	}
	close $fh or die "Cannot close $file: $!\n";
	return 1;
}

sub _mgs_ids {
	my ($file) = @_;
	return () unless defined($file) && -f $file;
	open my $fh, '<', $file or die "Cannot open MGS assignments $file: $!\n";
	my %ids;
	while (my $line = <$fh>) {
		chomp $line;
		next if $line =~ /^\s*$/;
		next if $line =~ /^#/;
		my @fields = split /\t/, $line, -1;
		next if $fields[0] eq 'Bin';
		die "Malformed MGS assignment in $file at line $.\n"
			unless @fields >= 2 && length($fields[0]) && length($fields[1]);
		$ids{$fields[0]} = 1;
	}
	close $fh or die "Cannot close MGS assignments $file: $!\n";
	return sort keys %ids;
}

sub _mgs_count {
	my ($file) = @_;
	my @ids = _mgs_ids($file);
	return scalar @ids;
}

sub _write_single_mgs_observations {
	my ($cluster_file, $observation_file) = @_;
	my @ids = _mgs_ids($cluster_file);
	die "Cannot synthesize observations for " . scalar(@ids) . " MGS in $cluster_file\n"
		unless @ids == 1;
	open my $fh, '>', $observation_file
		or die "Cannot write single-MGS observations $observation_file: $!\n";
	print {$fh} "Bin\tObservations\tQualTier\tMembers\n$ids[0]\t1\t1\t\n"
		or die "Cannot write single-MGS observations $observation_file: $!\n";
	close $fh or die "Cannot close $observation_file: $!\n";
}

sub _finish_without_mgs {
	my ($reason, $report_file, $checkpoint_file) = @_;
	open my $fh, '>', $report_file or die "Cannot write $report_file: $!\n";
	print {$fh} "MGS reconstruction completed without a usable MGS.\nReason: $reason\n"
		or die "Cannot write $report_file: $!\n";
	close $fh or die "Cannot close $report_file: $!\n";
	_touch_checkpoint($checkpoint_file, 'no-usable-mgs', $report_file);
	printL "MGS reconstruction completed without a usable MGS: $reason\n";
	close LOG or die "Cannot close MGS pipeline log: $!\n";
	exit 0;
}

# Legacy helper retained below the main routing.
sub reformat_4phylo{
	my ($FCF) =@_; #, $clusSelHR
	#my @allMGS=keys(%{$clusSelHR});
	#open I,"<$outD/MGSselection.txt" or die "Can't open $outD/MGSselection.txt\n"; 
	#while(<I>){chomp;push(@allMGS,$_);} close I;
	my $hr = readMGS($FCF); my %MGS = %{$hr};
	print "Selected ".scalar(keys %MGS) . " Bins for intrastrain phylo\n";
	open O,">$FCF.mgs" or die $!;
	my $MGcnt=0;
	foreach my $MG (keys %MGS){
		my $MG2 = $MG; $MG2 =~ s/_/:/;#for the MG2dram bins..
		$MG2 = $MG unless (exists($MGS{$MG2}));
		next unless (exists($MGS{$MG2}));
		$MG =~ s/:/_/; #make sure stupid : is completely gone...
		#print $MG."\n";
		$MGcnt++;
		print O $MG."\t".join(",",@{$MGS{$MG2}})."\n";
	}
	close O;
	print "Found $MGcnt MGs for withinstrain analysis\n";
}



sub invertIndex{
	my ($in,$out) = @_;
	open I,"<$in" or die "can;t open $in\n";
	open O,">$out" or die "can;t open $out\n";
	while (<I>){
		chomp; my @spl=split /\t/;
		$spl[1] =~ s/:/_/g;
		print O "$spl[1]\t$spl[0]\n";
	}
	close I; close O;
}
sub CanopyPrep{
	my ($inFc,$binCanDir) = @_;
	if ($canopyF eq ""){print"Canopy not requested, will skip step\n";return 0;}
	#read genes in canopies..
	my $ChkMevalF = "$inFc.filt$cmSuffix";
	return _mgs_count("$inFc.filt") if (-s $ChkMevalF);
	my %canCnts;
	printL "Prepping Canopy MGS (format, Bin quality)..\n";
	printL "$inFc\n";
	open my $canopy_input, '<', $inFc or die "can't open canopy file $inFc\n";
	while (my $line = <$canopy_input>){
		chomp $line;
		next if $line =~ /^\s*$/;
		my @spl = split /\t/, $line;
		die "Malformed Canopy assignment in $inFc at line $.\n" unless @spl >= 2 && length($spl[0]) && length($spl[1]);
		$canCnts{$spl[0]} ++;
	}
	close $canopy_input or die "Cannot close canopy file $inFc: $!\n";

	my %kept_canopies = map { $_ => 1 } grep { $canCnts{$_} >= 700 } keys %canCnts;
	my $canCNT = scalar keys %kept_canopies;
	my $filtered_tmp = "$inFc.filt.tmp.$$";
	open $canopy_input, '<', $inFc or die "can't reopen canopy file $inFc\n";
	open my $filtered_output, '>', $filtered_tmp
		or die "Cannot write filtered Canopy assignments $filtered_tmp: $!\n";
	while (my $line = <$canopy_input>) {
		next if $line =~ /^\s*$/;
		my ($canopy) = split /\t/, $line, 2;
		next unless $kept_canopies{$canopy};
		print {$filtered_output} $line
			or die "Cannot write filtered Canopy assignments $filtered_tmp: $!\n";
	}
	close $canopy_input or die "Cannot close canopy file $inFc: $!\n";
	close $filtered_output
		or die "Cannot close filtered Canopy assignments $filtered_tmp: $!\n";
	unlink "$inFc.filt" or die "Cannot replace $inFc.filt: $!\n" if -e "$inFc.filt";
	rename $filtered_tmp, "$inFc.filt"
		or die "Cannot publish filtered Canopy assignments $inFc.filt: $!\n";
	printL "Kept $canCNT/". scalar(keys(%canCnts)) . " Canopy MGS\n";
	return 0 unless $canCNT;
	createBinFAA($binCanDir,"$inFc.filt","$GCd/compl.incompl.$clusterID.prot.faa","faa");
	if ($useCheckM1 && !-s $ChkMevalF) {
		printL "running checkM on new Canopy MGS..\n";
		my $req_CMmem = 200;
		my $cmC = runCheckM($binCanDir,$ChkMevalF,"$nodeTmpD/cmCANO/",$numCore,0);
		my ($jobName2, $tmpCmd) = qsubSystem($logDir."/checkM.cano0.sh",
			$cmC,
			$numCore,int($req_CMmem)."G","ChMcano","","",1,[],\%QSBopt);
		push(@jobs2wait,$jobName2);
	}
	#checkM2
	if ( $useCheckM2 && !-s $ChkMevalF ){
		printL "running checkM2 on new Canopy MGS..\n";
		my $req_CMmem = 50;	my $cmC = "";
		$cmC .= runCheckM2($binCanDir,$ChkMevalF,"$nodeTmpD/cmCANO/",$canCore,0) ;	
		my ($jobName2, $tmpCmd) = qsubSystem($logDir."/checkM2.cano0.sh",
			$cmC,
			$canCore,int($req_CMmem)."G","ChM2cano","","",1,[],\%QSBopt);
		push(@jobs2wait,$jobName2);
	}

	#wait for checkM & checkM2
	qsubSystemJobAlive( \@jobs2wait,\%QSBopt );@jobs2wait = ();
	die "Canopy quality checking completed without producing $ChkMevalF\n" unless -s $ChkMevalF;
	return $canCNT;
}


sub getGoodMBstats{
	my $logfile = "$outD/MAG.$BinnerShrt.assStat.summary";
	printL "\n\nPer sample MAG stats are in $logfile\n\n";
	return if (-e $logfile && -s $logfile>0);# && -e $finalClusters);
	printL "recording per MAG assembly stats..";
	open L,">$logfile.1" or die $!;
	my $headWritten=0;
	foreach my $Doo (sort keys %DOs){
		my @smplIDs = @{$DOs{$Doo}{SmplID}}; my @paths = @{$DOs{$Doo}{wrdir}};
		my $metaGD = getAssemblPath($paths[-1]);
		my $MBf = $metaGD."/Binning/$BinnerShrt/$smplIDs[-1]";	my $MBfQual = $MBf.$cmSuffix;
		next unless (-e $MBfQual);
		my ($hr1,$hr2) = MB2assignedBinIds($MBf,$MBfQual);
		my %MB = %{$hr1};my %MBQ = %{$hr2}; my %valMBs;
		next if (scalar(keys %MB) == 0); #empty MAG file..
		foreach my $bin (keys %MB){
			if ($MBQ{$bin}{compl}< 80 || $MBQ{$bin}{conta}> 5 ){next;}
			$valMBs{$bin} = "$MBQ{$bin}{compl}\t$MBQ{$bin}{conta}";
		}
		open I,"<$MBf.assStat" or die "Can't open assStat $MBf.assStat\n";
		if (!$headWritten){
			$headWritten=1;
			my $hd = <I>; chomp $hd;
			print L "Sample\t$hd\tCM_Compl\tCM_Conta\n";
		}
		while (my $li = <I>){
			chomp $li;
			$li =~ m/^(\S+)\t/; next unless (exists($valMBs{$1}));
			print L "$smplIDs[-1]\t$li\t$valMBs{$1}\n";
		}
		close I;
	}
	close L;
	rename "$logfile.1", $logfile or die "Cannot replace $logfile with $logfile.1: $!\n";
	print " Done\n";
}




sub printL{
	my ($msg) = @_;
	print $msg;
	print LOG $msg;
}
