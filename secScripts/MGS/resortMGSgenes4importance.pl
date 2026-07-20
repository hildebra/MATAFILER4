#!/usr/bin/perl
#script that takes a selection of MGS genes (canopy format) and sorts them based on a) marker genes b) copy unmber c) overall occurence
#./resortMGSgenes4importance.pl /g/bork3/home/hildebra/data/SNP/GCs/DramaGCv5/ /g/bork3/home/hildebra/data/SNP/GCs/DramaGCv5/Binning/MetaBat/MB2.clusters.ext.can.Rhcl.mgs
use warnings;
use strict;
use Mods::geneCat qw(readGene2tax createGene2MGS);
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::GenoMetaAss qw(gzipopen readFasta writeFasta systemW);
use Mods::TamocFunc qw(readTabbed);
use Mods::math qw(meanArray medianArray);

sub evalCurMGS;


#v0.1: adopt .core MGS files to get additional info for sorting genes by importance
#v0.11: 9.2.24: retain more genes/MGS
#v0.12: 11.2.24: adopted to weighted multiBin scores; more subs to make script more modifiable
#v0.13: flush the final MGS and handle marker-free groups safely
#v0.14: count distinct samples, use multicopy evidence, and make ranking deterministic
my $version = 0.14;


# set up some base variables
die "Usage: $0 GC-dir MGS-file GTDB|FMG mode [cluster-ID]\n" unless @ARGV == 4 || @ARGV == 5;
my $rareBin = getProgPaths("rare");
my $GCd = $ARGV[0];
my $MGSfile = $ARGV[1];
my $useGTDBmg = $ARGV[2];
my $mode = $ARGV[3];
my $clusterID = @ARGV == 5 ? $ARGV[4] : 95;
die "cluster-ID must be between 1 and 100\n"
	unless $clusterID =~ /^\d+$/ && $clusterID >= 1 && $clusterID <= 100;
my $obsFile = $MGSfile; #$ARGV[3];# if (@ARGV > 3);
$obsFile =~ s/\.core$//; $obsFile.=".obs";
die "ARG 2 option has to be \"GTDB\" or \"FMG\"\n" unless ($useGTDBmg eq "GTDB" || $useGTDBmg eq "FMG");

#main output file
my $finout = "$MGSfile.srt";
if (-e $finout && -s $finout){
	#print $finout."\n"; exit(0);
	print "Overwriting $finout\n";
}


print "\n--------------------------------------------------\nResorting MGS genes for importance in strain phylo ver $version\n--------------------------------------------------\n";
#my @FMG40 = ("COG0012","COG0016","COG0018","COG0048","COG0049","COG0052","COG0080","COG0081","COG0085","COG0087","COG0088","COG0090","COG0091","COG0092","COG0093","COG0094","COG0096","COG0097","COG0098","COG0099","COG0100","COG0102","COG0103","COG0124","COG0172","COG0184","COG0185","COG0186","COG0197","COG0200","COG0201","COG0202","COG0215","COG0256","COG0495","COG0522","COG0525","COG0533","COG0541","COG0552");
#my %FMG40 = map { $_ => 1 } @FMG40;


die "Can't find main infile $MGSfile\n" unless (-s $MGSfile);

#read MGS occurrence to understand distribution; a missing observation table
#must not prevent small/partial MGS sets from being ranked.
my %MGSobs;
if (-s $obsFile) {
	my $hr = readTabbed($obsFile);
	%MGSobs = %{$hr};
} else {
	warn "MGS occurrence table $obsFile is unavailable; estimating prevalence from genes\n";
}

#load GTDB/FMG genes directly..
my $inMGFile="$GCd/FMG.subset.cats";
if ($useGTDBmg eq "GTDB"){
	$inMGFile="$GCd/GTDBmg.subset.cats";
}

#load list of reference marker genes (to mark these later as important genes)
my %MGset=();
open I,"<$inMGFile" or die "resortMGSgenes4importance.pl: Couldn't open $inMGFile\n"; 
my $totMGSgenes=0;
while (<I>){
	my @spl1 = split /\t/;
	my @spl2 = split /,/,$spl1[2];
	foreach my $gene (@spl2){$MGset{$gene} = 1;$totMGSgenes++;}
}
close I;
print STDERR "Loaded $totMGSgenes $useGTDBmg marker genes from $inMGFile\n";



#alt: go with compl.incompl.95.fna.clstr.idx to calc gene occurrences..
my %geneOcc;
my ($I,$ST) = gzipopen("$GCd/compl.incompl.$clusterID.fna.clstr.idx","gene cat index file");
#my $maxOcc = 0;
while (<$I>){
	chomp; my @spl= split /\t/;
	if (@spl<2){next; }#die $1." no tab char\n";}
	my %samples;
	for my $member (split /,/, $spl[1]) {
		$member =~ s/^>//;
		my ($sample) = split /__/, $member, 2;
		$samples{$sample} = 1 if defined($sample) && length($sample);
	}
	$geneOcc{$spl[0]} = scalar(keys %samples);
	#$maxOcc = $spl[1] if ($maxOcc < $spl[1]);
}
close $I;
print STDERR "Counted gene occurrence for " . scalar(keys(%geneOcc)) . " genes\n";

#my $gene2taxF = createGene2MGS($MGSfile,$GCd);
#my %gen2Bin ; my %gene2COG;
#open I,"<$gene2taxF" or die "Cant' open $gene2taxF\n";
#while (<I>){
#	chomp; my @spl = split /\t/;
#	push (@{$gen2Bin{$spl[1]}}, $spl[0]);
#	$gene2COG{$spl[0]} = $spl[2] if (defined ($spl[2]));
#}
#close I;


#my $sortedOutFile = "$MGSfile.srt";
open O,">$finout" or die "can't open outfile $finout\n";
open I,"<$MGSfile" or die "cant open infil $MGSfile\n";
my $cn=0; my $MGScnt = 0; my $geneCnt=0;
my $curMGS=""; 
my %occ; my %multiCp; my %markers; my %multiBin; 
#foreach my $mg (keys %gen2Bin){
while (my $line = <I>){
	chomp $line;
	next if ($line =~ m/^#/); #commented line
	my @spl = split /\t/,$line;
	my $MGS = $spl[0];  
	$curMGS = $MGS if ($curMGS eq ""); 
	if ($MGS eq ""){die "Undefine MGS on line $line\n";next;}

	if ($MGS ne $curMGS){
		my $retS = evalCurMGS($MGS);
		print O $retS;
	}
	die "Malformed MGS row: $line\n" unless @spl >= 6 && defined $spl[1] && length $spl[1];
	my $gene = $spl[1];
	#push @genes,$gene; 
	if (exists($multiBin{$gene})) {
		warn "Ignoring duplicate gene $gene in MGS $MGS\n";
		next;
	}
	if (exists($MGset{$gene})){
		$markers{$gene}=$spl[2];
	} else {
		$occ{$gene} = $spl[2]; 
	}		
	$multiCp{$gene} = $spl[3]; $multiBin{$gene} = $spl[4];
	$cn ++;
	

}
print O evalCurMGS("") if $curMGS ne "";
close O;
close I;
#report that all went fine
print  "Finished \nProcessed $cn genes, used $geneCnt genes in $MGScnt MGS\nSaved in $finout\n";


exit(0);



sub evalCurMGS{
	#this routine decides which MGS genes (already pre-filtered for core genes) will be handed on to strain phylo construction.. should be "certain" cutoffs for removing genes (intra phylo will do another round of filtering)
	my ($MGS) = @_;
	my %finalList;
	my $mrkCnt=0;
	my @all_genes = (keys(%markers), keys(%occ));
	my @known_occ = map { $geneOcc{$_} } grep { defined($geneOcc{$_}) && $geneOcc{$_} > 0 } @all_genes;
	my @marker_occ = map { $geneOcc{$_} } grep { defined($geneOcc{$_}) && $geneOcc{$_} > 0 } keys %markers;
	my $expected_occ = @marker_occ >= 3 ? medianArray(@marker_occ)
		: @known_occ ? medianArray(@known_occ) : 1;
	my $MGSob = exists($MGSobs{$curMGS}) ? $MGSobs{$curMGS} : $expected_occ;

	my $copy_fraction = sub {
		my ($gene) = @_;
		my $observations = exists($markers{$gene}) ? $markers{$gene} : $occ{$gene};
		return $observations && $observations > 0 ? ($multiCp{$gene} || 0) / $observations : 0;
	};
	my $reject_gene = sub {
		my ($gene) = @_;
		my $observations = exists($markers{$gene}) ? $markers{$gene} : $occ{$gene};
		return 1 if defined($multiBin{$gene}) && $multiBin{$gene} >= 3;
		# A single duplicate observation is weak evidence in undersampled MGS.
		return 1 if $observations >= 5 && ($multiCp{$gene} || 0) >= 2
			&& $copy_fraction->($gene) > 0.2;
		# Only reject gross prevalence mismatches when enough samples exist.
		if ($expected_occ >= 4 && defined($geneOcc{$gene})) {
			return 1 if $geneOcc{$gene} < 0.25 * $expected_occ
				|| $geneOcc{$gene} > 4 * $expected_occ;
		}
		return 0;
	};
	my $rank_genes = sub {
		my ($left, $right, $occ_hr) = @_;
		return $copy_fraction->($left) <=> $copy_fraction->($right)
			|| ($multiBin{$left} // 0) <=> ($multiBin{$right} // 0)
			|| ($occ_hr->{$right} // 0) <=> ($occ_hr->{$left} // 0)
			|| abs(($geneOcc{$left} // $expected_occ) - $expected_occ)
				<=> abs(($geneOcc{$right} // $expected_occ) - $expected_occ)
			|| $left cmp $right;
	};

	my @mrks = sort { $rank_genes->($a, $b, \%markers) } keys %markers;
	for my $gn (@mrks) {
		next if $reject_gene->($gn);
		$finalList{$gn} = scalar(keys %finalList);
		$mrkCnt++;
	}
	my @srtedGenes = sort { $rank_genes->($a, $b, \%occ) } keys %occ;
	for my $gn (@srtedGenes) {
		next if $reject_gene->($gn);
		$finalList{$gn} = scalar(keys %finalList);
	}

	my $medMBi = @all_genes ? medianArray(map { $multiBin{$_} } @all_genes) : 0;
	my $avgMBi = @all_genes ? meanArray([map { $multiBin{$_} } @all_genes]) : 0;
	print "${curMGS} (".scalar(keys %multiBin)."):: " ;
	print scalar(keys %finalList) ." genes, $mrkCnt markerGs used, expected prevalence: "
		. int($expected_occ*100)/100 . ", MAG occ: $MGSob, median/avg multiBin $medMBi/"
		. int($avgMBi*100)/100 ."\n";
	
	$geneCnt+= scalar(keys %finalList);
	
	my @finalGs = sort { $finalList{$a} <=> $finalList{$b} } keys(%finalList);
	#my @finalGs = @finalList{@keys};
	#die "\n@finalGs\n";

	my $retStr= $curMGS."\t".join(",", @finalGs)."\n";
	
	
	$curMGS = $MGS; 
	%occ = (); %multiCp= (); %markers= (); %multiBin= ();
	$MGScnt++;
	return $retStr;
}





