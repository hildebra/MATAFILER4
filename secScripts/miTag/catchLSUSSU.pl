#!/usr/bin/env perl
#./catchLSUSSU.pl /g/scb/bork/hildebra/data2/Soil_finland/Project_150108/Sample_S141_170/S141_170_TCCGCGAA-ATAGAGGC_L002_R1_001.fastq.gz /g/scb/bork/hildebra/data2/Soil_finland/Project_150108/Sample_S141_170/S141_170_TCCGCGAA-ATAGAGGC_L002_R2_001.fastq.gz /g/scb/bork/hildebra/data2/Soil_finland/tmp_16s/ /g/scb/bork/hildebra/data2/Soil_finland/tmp_16s/aligned/ 10 1
# ./catchLSUSSU.pl '/scratch/bork/hildebra/SimuB/simulated_metaG3SA_0//tmp/seqClean/filtered.1.fastq' '/scratch/bork/hildebra/SimuB/simulated_metaG3SA_0//tmp/seqClean/filtered.2.fastq' /tmp/hildebra/simulated_metaG3SA_0ITS//SMRNA/ /g/scb/bork/hildebra/SNP///SimuB/simulated_metaG3SA_0/ribos/ 12 Sb1 0 /scratch/bork/hildebra/SimuB//rnaDB/

use strict;
use warnings;
use Getopt::Long qw( GetOptions );
use Mods::GenoMetaAss qw(renameFastqCnts systemW);
use Mods::IO_Tamoc_progs qw(inputFmtSpades getProgPaths);
use File::Copy qw(copy);
use File::Path qw(make_path remove_tree);
use File::Spec;

sub smrnaRunCmd;
sub copySortmernaOutputs;
sub touchFile;
sub shellQuote;


#18.5.26: added versioning to 0.1
my $cLSUSSUver = 0.2;

#my $tmpP = "/g/scb/bork/hildebra/data2/Soil_finland/tmp_16s/";
my $tmpP = "";#$ARGV[3];
my $alignPath = "";#$ARGV[4];
my $threads = 1;#$ARGV[5];
my $smpN = "";#$ARGV[6];
my $doRiboAssembl = 0;#$ARGV[7];
my $path2DB = "";#$ARGV[8];
my 	$read1 ="";
my $read2 = "";
my $readS = "";


GetOptions(
	"R1=s" => \$read1, 
	"R2=s" => \$read2, 
	"RS=s" => \$readS,

	"alignDir=s" => \$alignPath,
	"smplID=s"      => \$smpN,
	"tmpDir=s"      => \$tmpP,
	"cores=i" => \$threads,
	#"DBdir=s" => \$path2DB,
	"assmblRibos=i" => \$doRiboAssembl,
) or die("Error in command line arguments\n");
die "Unexpected positional arguments: @ARGV\n" if @ARGV;
die "-alignDir is required\n" if $alignPath eq "";
die "-tmpDir is required\n" if $tmpP eq "";
die "-smplID is required\n" if $smpN eq "";
die "-cores must be a positive integer\n" if $threads < 1;
die "Ribosomal assembly is no longer supported\n" if $doRiboAssembl;

my @r1i = grep { $_ ne "" && $_ ne "-1" } split(",",$read1);
my @r2i = grep { $_ ne "" && $_ ne "-1" } split(",",$read2);
my @rSi = grep { $_ ne "" && $_ ne "-1" } split(",",$readS);
die "At least one -R1 or -RS input is required\n" unless @r1i || @rSi;
die "-R1 and -R2 must contain the same number of files\n" if @r2i && @r1i != @r2i;
die "-R2 was supplied without -R1\n" if @r2i && !@r1i;
for my $readFile (@r1i, @r2i, @rSi){
	die "Read input does not exist: $readFile\n" unless -e $readFile;
}
my $rS=join(",",@rSi);
my $singlMode = @r2i ? 0 : 1;

$alignPath = File::Spec->canonpath(File::Spec->rel2abs($alignPath));
my $tmpRoot = File::Spec->canonpath(File::Spec->rel2abs($tmpP));
make_path($alignPath) unless -d $alignPath;
make_path($tmpRoot) unless -d $tmpRoot;
die "-alignDir is not a directory: $alignPath\n" unless -d $alignPath;
die "-tmpDir is not a directory: $tmpRoot\n" unless -d $tmpRoot;
my $safeSample = $smpN;
$safeSample =~ s/[^A-Za-z0-9_.-]+/_/g;
$tmpP = File::Spec->catdir($tmpRoot, "catchLSUSSU_${safeSample}_$$");
make_path($tmpP);

if (-e "$alignPath/SSU_pull.sto" && -e "$alignPath/LSU_pull.sto"){
	print "All riboFind SortMeRNA targets are complete\n";
	remove_tree($tmpP);
	exit(0);
}

#if ($path2DB eq ""){die "database not defined (-DBdir) ! \n";}
#if (@ARGV<8){die "Not enough input arguments!!\n";}


#my $smrPath = getProgPaths("srtMRNA_path");
#my $path2DB = "$smrPath/rRNA_databases/";
#my $smrnaBin = "$smrPath/./sortmerna";

my $smrnaBin = getProgPaths("sortmerna");
announce($smrnaBin);
#my $mergeScript = getProgPaths("mergeRdScr");
#my $unmergeScript = getProgPaths("unmergeRdScr");
#my $spadesBin = getProgPaths("spades");
my $ltslcaP = "$alignPath/ltsLCA/";
if (!-e "$alignPath/SSU_pull.sto" || !-e "$alignPath/LSU_pull.sto") { #unpack reads
	#preparation of reads
	#my @r1i = split(";",$ARGV[0]); my $r1="$tmpP/read1.tmp.fq";
	#my @r2i = split(";",$ARGV[1]); my $r2="$tmpP/read2.tmp.fq";
	#my @rSi = split(";",$ARGV[2]); my $rS="$tmpP/readSingl.tmp.fq";
	#die $singlMode."\n";
	#die "$r1i[0]\n";
	
	my $read1Str = join(",", @r1i);
	my $read2Str = join(",", @r2i);
	
	
	
	#my $ITSDBfa = getProgPaths("ITSdbFA"); $ITSDBfa =~ m/([^\/]+)$/; $ITSDBfa = $1; my $ITSDBidx = $ITSDBfa; $ITSDBidx =~ s/\.fa*$/\.idx/;

	#my $refDBits = getProgPaths("ITSdbFAsrt");
	my $refDBssu = getProgPaths("SSUdbFAsrt");
	my $refDBlsu = getProgPaths("LSUdbFAsrt");
	my $idxSSU   = getProgPaths("SSUidx", 0);
	my $idxLSU   = getProgPaths("LSUidx", 0);

	#my $curStone = "$alignPath/ITS_pull.sto";
	#unless (-e $curStone){ #ITS seems to be in general unreliable (too diverse?)
	#	if ($refDBits ne ""){
	#		my $ltag = "reads_ITS";
	#		my $runner = smrnaRunCmd($tmpP."/$ltag",$refDBits,$read1Str,$read2Str,$alignPath,$singlMode,$kvdbITS);
	#		$runner .= "\n\n". smrnaRunCmd($tmpP."/$ltag",$refDBits,$rS,"",$alignPath,1,$kvdbITS) if ($rS ne "");
	#		systemW $runner;
	#		outfileCpy($tmpP."/$ltag",$alignPath);
	#		systemW "rm -f $ltslcaP/ITS_ass.sto" if (-e " $ltslcaP/ITS_ass.sto");
	#		systemW "touch $curStone" unless (`wc -l $alignPath/$ltag.r1.fq | cut -f1 -d' '` != `wc -l $alignPath/$ltag.r1.fq | cut -f1 -d' '` );
	#	} else {
	#		print "Skipping ITS since DB could not be found\n"
	#	}
	#}
	print "Skipping ITS\n";

	my $curStone = "$alignPath/SSU_pull.sto";
	unless (-e $curStone){
		my $ltag = "reads_SSU";
		my $runner = smrnaRunCmd($tmpP."/$ltag",$refDBssu,$read1Str,$read2Str,$alignPath,$singlMode,$idxSSU);
		$runner .= "\n\n". smrnaRunCmd($tmpP."/$ltag",$refDBssu,$rS,"",$alignPath,1,$idxSSU) if ($rS ne "");
		systemW $runner;
		copySortmernaOutputs($tmpP, $alignPath, $ltag, $singlMode || $rS ne "");
		touchFile($curStone);
		
		#renameFastqCnts($alignPath."/reads_SSU.r1.fq",$smpN."__SSU"); renameFastqCnts($alignPath."/reads_SSU.r2.fq",$smpN."__SSU");
		#system "touch $curStone";
	}
	$curStone = "$alignPath/LSU_pull.sto";
	unless (-e $curStone){
		my $ltag = "reads_LSU";
		my $runner = smrnaRunCmd($tmpP."/$ltag",$refDBlsu,$read1Str,$read2Str,$alignPath,$singlMode,$idxLSU);
		$runner .= "\n\n". smrnaRunCmd($tmpP."/$ltag",$refDBlsu,$rS,"",$alignPath,1,$idxLSU) if ($rS ne "");
		#print $runner."\n";
		unlink "$ltslcaP/LSU_ass.sto" if -e "$ltslcaP/LSU_ass.sto";
		systemW $runner;
		#if (system $runner) {print "Error in $runner\n"; exit(89);}#unless (system $runner) {die "Failed\n$runner\n";}
		#outfileCpy($tmpP."/$ltag",$alignPath);
		#renameFastqCnts($alignPath."/reads_LSU.r1.fq",$smpN."__LSU"); renameFastqCnts($alignPath."/reads_LSU.r2.fq",$smpN."__LSU");
		#system "touch $curStone";
		copySortmernaOutputs($tmpP, $alignPath, $ltag, $singlMode || $rS ne "");
		touchFile($curStone);

	}


}

remove_tree($tmpP) if -d $tmpP;

#ribo assemblies.. difficult as high chance for chimeras.. maybe switch assembler later?
my $outP = $alignPath;
if ($doRiboAssembl && (!-d $outP."/Ass_ITS" || !-e $outP."/Ass_ITS/scaffolds.fasta" || !-e "$outP/Ass/allAss.sto") ){#assembly of ITS, LSU, SSU seperately
	die "catchLSUSSU.pl::riboassembly no longer supported\n";
	# my $K =  "-k 27,33,55,71,125";
	# my $logDir = "$outP/Ass/logs/";
	# my $Scmd = "mkdir -p $logDir\n\n";
	# my $nodeTmp = $outP."Ass/SSU/";
	# if (!-z "$outP/reads_SSU.r1.fq"){ #check if input is empty
		# $Scmd .= "\nmkdir -p $nodeTmp\n$spadesBin $K ".inputFmtSpades([$outP."/reads_SSU.r1.fq"],[$outP."/reads_SSU.r2.fq"],[],$logDir)." -t $threads --meta -m 60 ";#SSU --mismatch-correction 
		# $Scmd .= "-o $nodeTmp\nrm -f -r $nodeTmp/K* $nodeTmp/tmp $nodeTmp/mismatch_corrector $nodeTmp/corrected $nodeTmp/misc \n";
	# }
	# $nodeTmp = $outP."Ass/LSU/";
	# if (!-z "$outP/reads_LSU.r1.fq"){
		# $Scmd .= "\nmkdir -p $nodeTmp\n$spadesBin $K ".inputFmtSpades([$outP."/reads_LSU.r1.fq"],[$outP."/reads_LSU.r2.fq"],[],$logDir)." -t $threads --meta -m 60 ";#LSU
		# $Scmd .= "-o $nodeTmp\nrm -f -r $nodeTmp/K* $nodeTmp/tmp $nodeTmp/mismatch_corrector $nodeTmp/corrected $nodeTmp/misc \n";
	# }
	# $nodeTmp = $outP."Ass/ITS/";
	# if (!-z "$outP/reads_ITS.r1.fq"){
		# $Scmd .= "\nmkdir -p $nodeTmp\n$spadesBin $K ".inputFmtSpades([$outP."/reads_ITS.r1.fq"],[$outP."/reads_ITS.r2.fq"],[],$logDir)." -t $threads --meta -m 60 ";#ITS
		# $Scmd .= "-o $nodeTmp\nrm -f -r $nodeTmp/K* $nodeTmp/tmp $nodeTmp/mismatch_corrector $nodeTmp/corrected $nodeTmp/misc \n";
	# }
	# $Scmd .= "touch $outP/Ass/allAss.sto\n";
	# systemW $Scmd;
}

 print "Finished pull out \n";
 print "(and eventually assembly)\n" if ($doRiboAssembl);
 exit(0);





sub smrnaRunCmd( $ $ $ $ $ $ $){
	my ($outFile, $refDB, $R1, $R2, $finD, $isSingl, $idxDir) = @_;
	$idxDir = "" unless defined $idxDir;

	return "" if ($R1 eq "");

	# Build --ref args from colon-separated FASTA paths
	my $refStr = join(" ", map { "--ref ".shellQuote($_) } split(":", $refDB));
	my $idxStr = ($idxDir ne "") ? "--idx-dir ".shellQuote($idxDir)." --index 0" : "";

	my $cmd = "";
	if ($isSingl){
		my @singleReads = grep { $_ ne "" } split(",", $R1);
		for (my $i = 0; $i < @singleReads; $i++){
			my $pfx = @singleReads > 1 ? "${outFile}.tmp${i}" : $outFile;
			$cmd .= "$smrnaBin --reads ".shellQuote($singleReads[$i])." $refStr $idxStr";
			$cmd .= " --kvdb ".shellQuote("${outFile}.kvdb${i}")." --readb ".shellQuote("${outFile}.readb${i}")." --aligned ".shellQuote($pfx);
			$cmd .= " --fastx --threads $threads -e 1e-12 --num_alignments 1 --no-best \n";
			$cmd .= "rm -rf ".shellQuote("${outFile}.kvdb${i}")." ".shellQuote("${outFile}.readb${i}")."\n";
		}
		if (@singleReads > 1){
			my $outputs = join(" ", map { shellQuote("${outFile}.tmp${_}.fq.gz") } 0..$#singleReads);
			$cmd .= "cat $outputs > ".shellQuote("${outFile}.fq.gz")." && rm -f $outputs\n";
		}
	} else {
		return "" if ($R2 eq "");
		my @r1s = split(",", $R1);
		my @r2s = split(",", $R2);
		die "Internal error: unequal paired-read lists\n" unless @r1s == @r2s;
		for (my $i = 0; $i < @r1s; $i++){
			my $pfx = (@r1s > 1) ? "${outFile}.tmp${i}" : $outFile;
			$cmd .= "$smrnaBin --reads ".shellQuote($r1s[$i])." --reads ".shellQuote($r2s[$i])." $refStr $idxStr";
			$cmd .= " --kvdb ".shellQuote("${outFile}.kvdb${i}")." --readb ".shellQuote("${outFile}.readb${i}")." --aligned ".shellQuote($pfx);
			$cmd .= " --fastx --threads $threads -e 1e-12 --num_alignments 1 --no-best --paired_in --out2 \n";
			$cmd .= "rm -rf ".shellQuote("${outFile}.kvdb${i}")." ".shellQuote("${outFile}.readb${i}")."\n";
		}
		if (@r1s > 1){
			# multiple input pairs: cat per-pair outputs into final files
			my $fwdF = join(" ", map { "'${outFile}.tmp${_}_fwd.fq.gz'" } 0..$#r1s);
			my $revF = join(" ", map { "'${outFile}.tmp${_}_rev.fq.gz'" } 0..$#r1s);
			$cmd .= "cat $fwdF > '${outFile}.r1.fq.gz' && rm -f $fwdF\n";
			$cmd .= "cat $revF > '${outFile}.r2.fq.gz' && rm -f $revF\n";
		} else {
			# v4 --out2 writes _fwd.fq / _rev.fq; rename to .r1.fq / .r2.fq
			# NOTE: verify these suffixes match your sortmerna 4.3.4 build if output is missing
			$cmd .= "mv -f '${outFile}_fwd.fq.gz' '${outFile}.r1.fq.gz'\n";
			$cmd .= "mv -f '${outFile}_rev.fq.gz' '${outFile}.r2.fq.gz'\n";
		}
	}
	return $cmd;
}


sub make_interleave{
	die "deprecated make_interleave\n";
}

sub copySortmernaOutputs{
	my ($sourceDir, $destinationDir, $tag, $includeSingle) = @_;
	for my $oldFile (glob(File::Spec->catfile($destinationDir, "${tag}.r*.fq.gz")),
		glob(File::Spec->catfile($destinationDir, "${tag}.fq.gz"))){
		unlink $oldFile or die "Cannot remove stale SortMeRNA output $oldFile: $!\n";
	}
	my @expected = ("${tag}.r1.fq.gz", "${tag}.r2.fq.gz");
	push @expected, "${tag}.fq.gz" if $includeSingle;
	for my $fileName (@expected){
		my $source = File::Spec->catfile($sourceDir, $fileName);
		my $destination = File::Spec->catfile($destinationDir, $fileName);
		if (-e $source){
			copy($source, $destination) or die "Cannot copy $source to $destination: $!\n";
		} else {
			touchFile($destination); # successful search with no hits
		}
	}
}


sub touchFile{
	my ($path) = @_;
	open my $touchHandle, ">", $path or die "Cannot create checkpoint/output $path: $!\n";
	close $touchHandle or die "Cannot close checkpoint/output $path: $!\n";
}


sub shellQuote{
	my ($value) = @_;
	$value = "" unless defined $value;
	$value =~ s/'/'"'"'/g;
	return "'$value'";
}


sub announce{
	my ($configuredSortmerna) = @_;
	print "catchLSUSSU v $cLSUSSUver\n";
	my $status = systemW("$configuredSortmerna --version", 0);
	warn "Could not query the configured SortMeRNA version\n" if $status;
}
