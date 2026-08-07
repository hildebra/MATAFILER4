#!/usr/bin/perl
#script that evaluates Bin qual within single assembly via checkM or checkM2, gets assembly stats per Bin and makes assembly-bins ready for collection via compoundBinning.pl
#perl extractMeBat2.pl /g/bork3/home/hildebra/data/SNP/GNMass3_a/alien-11-883-0/assemblies/metag/scaffolds.fasta.filt /g/bork3/home/hildebra/data/SNP/GNMass3_a/alien-11-883-0/assemblies/metag/Binning/MB2/MM20 /tmp/XX
use warnings;
use strict;
use Getopt::Long qw( GetOptions );

use Mods::GenoMetaAss qw(  systemW readFasta);
use Mods::Binning qw(runCheckM runCheckM2 MB2N50);
use Mods::math qw(medianArray);
sub MB2assigns; #sub MB2N50;
sub createBinFAA;

#0.11: 1.4.26: added getOpt interface, Pilea integration
my $version= 0.11;

my $refFA = "";# $ARGV[0];
my $MB2 = ""; #$ARGV[1];
my $tmpD = ""; #$ARGV[2];
my $ncore = 1;
#$ncore = $ARGV[3] if (@ARGV > 3);
my $usCheckM2 = 0;
#$usCheckM2 = $ARGV[4] if (@ARGV > 4);
my $usCheckM1 = 1;
#$usCheckM1 = $ARGV[5] if (@ARGV > 5);
my $BinnerChoice = 0;
#$BinnerChoice = $ARGV[6] if (@ARGV > 6);
my $read1raw = "";my $read2raw = ""; my $readSraw = ""; #comma delimted list of reads.. only needed for pilea..

if (@ARGV && $ARGV[0] !~ /^-/) {
	die "Legacy usage: $0 assembly bin-assignments tmp-dir [cores [checkM2 [checkM1 [binner]]]]\n"
		if @ARGV < 3 || @ARGV > 7;
	($refFA, $MB2, $tmpD) = splice(@ARGV, 0, 3);
	$ncore = shift @ARGV if @ARGV;
	$usCheckM2 = shift @ARGV if @ARGV;
	$usCheckM1 = shift @ARGV if @ARGV;
	$BinnerChoice = shift @ARGV if @ARGV;
}

GetOptions(
	"asm=s"  => \$refFA,
	"tmpD=s" => \$tmpD,
	"binF=s"  => \$MB2,
	"ncore=i" => \$ncore,
	"checkM2=i" => \$usCheckM2,
	"checkM1=i" => \$usCheckM1,
	"binner=i" => \$BinnerChoice,
	"read1=s" => \$read1raw,
	"read2=s" => \$read2raw,
	"readS=s" => \$readSraw,
) or die "Invalid checkBinQual.pl options\n";
die "Unexpected positional arguments: @ARGV\n" if @ARGV;
die "-asm, -tmpD and -binF are required\n"
	unless length($refFA) && length($tmpD) && length($MB2);
die "Assembly is missing or empty: $refFA\n" unless -s $refFA;
die "-ncore must be a positive integer\n" unless $ncore > 0;
die "-binner must be one of 0..5\n"
	unless $BinnerChoice =~ /^\d+$/ && $BinnerChoice >= 0 && $BinnerChoice <= 5;
die "Select at least one of -checkM1 or -checkM2\n" unless $usCheckM1 || $usCheckM2;

my $MB2Dir = $MB2; $MB2Dir =~ s/[^\/]+$//;

#die "$ncore\n";
#die "$tmpD\n";
print "Bin postprocessing v$version\n";
print "Extracting Bins, Qual check, reformatting into $MB2\n";

system "mkdir -p $tmpD" unless (-d $tmpD);
my $binD = "$tmpD/bins/";
system "mkdir -p $binD" unless (-d $binD);



#$isSemiBin = 1 if (-d "$MB2Dir/output_recluster_bins"); #!-e $MB2 &&  #likely SemiBin outdir 

my $emptyBin=0; #anything to do here, or no Bins found?
my %MB;

if ($BinnerChoice == 2 && !-s $MB2){
	#needs to create metabat like file..
	print "Assuming SemiBat output dir..\n" if ($BinnerChoice == 2);
	#first prepare to delete all unused files..
	system "mv $MB2Dir/* $tmpD/";
	system "cp $tmpD/*.stone $MB2Dir";
	if (-d "$tmpD/output_recluster_bins"){
		system "mv $tmpD/output_recluster_bins/* $binD" ;
	} elsif (-d "$tmpD/output_bins"){
		system "mv $tmpD/output_bins/* $binD" ;
	}
	
	my $repStr = "";
	opendir(DIR, $binD) or die "Could not open $binD\n";
	while (my $filename = readdir(DIR)) {
		#print "$binD/$filename\n";
		
		my $binN = $filename; $binN =~ s/\.fa$//;
		next unless (-f "$binD/$filename" && -s "$binD/$filename");
		my $FR = readFasta("$binD/$filename");
		my %FAS = %{$FR};
		foreach my $c (keys %FAS){
			$repStr .= "$c\t$binN\n";
		}
	}
	closedir(DIR);
	#write actual output, if something to report..
	if ($repStr ne ""){
		open O,">$MB2" or die $!;
		print O $repStr;
		close O;
	}elsif (!-e $MB2){
		system "touch $MB2";
	}
} elsif ($BinnerChoice) {#MetaBat2 processing..
	print "Detected standardized binner assignment output..\n";
} else {
	print "Asssuming generic binner..\n";
}

#read in contig to bin assignments
system "touch $MB2" unless (-e $MB2);
my $hr = MB2assigns($MB2);
%MB = %{$hr};

#standardized path to recreate bin groups, creates contigs per bin, 1 file each bin
createBinFAA($binD,$refFA);

my $outFile = $MB2.".cm";
my $outFile2 = $MB2.".cm2";

print "found ". int(keys %MB) ." Bins\n";

if ($usCheckM1){
	my $tmpD2 = $tmpD."/CM/";
	if ($emptyBin){#no bins in metag..
		system "touch $outFile";
	} else {
		runCheckM($binD,$outFile,$tmpD2,$ncore,1,"fna") unless (-e $outFile);
	}
}
if ($usCheckM2){
	my $tmpD2 = $tmpD."/CM2/";
	if ($emptyBin){
		#system "touch $outFile2";
		open O,">$outFile2";
		print O "Name\tCompleteness\tContamination\tCompleteness_Model_Used\tTranslation_Table_Used\tAdditional_Notes\n";
		close O;

	} else {
		runCheckM2($binD,$outFile2,$tmpD2,$ncore,1,"fna") unless (-e $outFile2);
	}
}
#and get N50 etc vals for each contig
if (!-e "$MB2.assStat"){
	my $hr = MB2N50(\%MB);
	my %asS = %{$hr};
	open O,">$MB2.assStat" or die $!;
	print O "MB2\ttotalL\tmeanL\tctgN\tN20\tN50\tN80\tG1k\tG10k\tG100k\tG1M\n";
	my @itKeys = qw(tL meanL cN N20 N50 N80 1K 10K 100K 1M);
	foreach my $ak (keys %asS){
		#print "$ak\n"; print "$it\n";
		print O "$ak";
		foreach my $it(@itKeys){print O "\t$asS{$ak}{$it}";} #
		print O "\n";
	}
	close O;
}



#DONE





sub createBinFAA($$){ #old version, no longer used
	my ($binD, $refFA) = @_;
	$hr = readFasta($refFA);
	system "mkdir -p $binD" unless (-d $binD);
	my %FAS = %{$hr};

	$emptyBin = 1 unless keys %MB;

	foreach my $bin (keys %MB){
		my @ctgs = @{$MB{$bin}};
		open O,">$binD/$bin.fna" or die $!;
		foreach my $ctg (@ctgs){
			die "can't find contig $ctg\n" unless (exists $FAS{$ctg});
			print O ">$ctg\n$FAS{$ctg}\n";
		}
		close O;
	}
	undef %FAS ;
}


sub MB2assigns($){
	my ($inF) = @_;
	my %ret;
	open I,"<$inF" or die "Can't open maxbin2 output $inF\n";
	while (<I>){
		chomp;
		next if /^\s*$/;
		my @spl  = split /\t/, $_, -1;
		die "Malformed bin assignment in $inF: $_\n"
			unless @spl >= 2 && length($spl[0]) && length($spl[1]);
		next if ($spl[1] eq "0");
		next if ($spl[0] eq "Sequence ID");
		push(@{$ret{$spl[1]}}, $spl[0]);
	}
	close I;
	return \%ret;
}
