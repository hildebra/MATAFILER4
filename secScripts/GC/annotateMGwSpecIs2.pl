#!/usr/bin/env perl
#annotates specI's in the dataset
#also creates specI abundance tables
#will create tax abundance table 
#dont forget to add custom genomes to specI database!! -> this was now replaced by 3rd arg being a MGS list
#  perl annotateMGwSpecIs.pl /g/bork3/home/hildebra/data/SNP/GCs/T2_HM3_GNM3_ABR 12 /g/bork3/home/hildebra/data/SNP/GCs/T2_HM3_GNM3_ABR/Canopy3/clusters.txt
# perl annotateMGwSpecIs.pl /g/bork3/home/hildebra/data/SNP/GCs/alienGC2 12 /g/bork3/home/hildebra/data/SNP/GCs/alienGC2/Binning/MetaBat/MB2.clusters.core /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/Binning/MetaBat/extended_Tax.txt
# perl annotateMGwSpecIs.pl /g/bork3/home/hildebra/data/SNP/GCs/DramaGCv5/ 12 /g/bork3/home/hildebra/data/SNP/GCs/DramaGCv5/Binning/MetaBat/MB2.clusters.ext.can.Rhcl.mgs.srt /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/Binning/MetaBat/extended_Tax.txt


use warnings;
use strict;

use Mods::GenoMetaAss qw( readClstrRev systemW median mean readFasta gzipopen);
use Mods::Subm qw(qsubSystem emptyQsubOpt );
use Mods::geneCat qw(calculate_spearman_correlation read_matrix correlation);
use Mods::FuncTools qw(passBlast lambdaBl);
use Mods::TamocFunc qw( readTable);
use Mods::math qw(nonZero);


use Mods::IO_Tamoc_progs qw(getProgPaths);
use List::Util;
use Getopt::Long qw( GetOptions );

#.11: reverted some of the drastic gene filtering
#.12 adapted for Bin_SB folder potentially being different..
#.13 validate inputs and repair MGS assignment/output handling
#.14 exclude marker genes shared by multiple MGS instead of assigning by input order
#.15 resolve assignments through one canonical gene/specI/COG index
my $version = 0.15;

my $matrixStorage;
my $FMGmatrix;
sub loadMarkerMatrix($){
	my ($matrixFile) = @_;
	return read_matrix($matrixFile) if $matrixStorage eq "dense";
	my ($nonzero,$cells) = $matrixStorage eq "auto" ? matrixDensity($matrixFile) : (undef,undef);
	my $useSparse = $matrixStorage eq "sparse" || ($cells && $nonzero / $cells <= 0.15);
	if ($matrixStorage eq "auto"){
		my $density = $cells ? 100 * $nonzero / $cells : 0;
		print sprintf("Matrix density %.2f%%; using %s storage\n",$density,$useSparse ? "sparse" : "dense");
	}
	return $useSparse ? readSparseMatrix($matrixFile) : read_matrix($matrixFile);
}

sub matrixDensity($){
	my ($matrixFile) = @_;
	my ($I,$status) = gzipopen($matrixFile,"Matrix File",1);
	die "Could not open matrix $matrixFile\n" unless ($status && defined $I);
	my ($lineCount,$nonzero,$cells) = (0,0,0);
	while (my $line = <$I>){
		chomp $line;
		$lineCount++;
		next if $lineCount == 1;
		my @values = split /\t/,$line,-1;
		shift @values;
		foreach my $value (@values){
			next if !defined($value) || $value eq "";
			$cells++;
			$nonzero++ if $value != 0;
		}
	}
	close $I;
	die "not enough lines in matrix $matrixFile\n" if $lineCount <= 2;
	return ($nonzero,$cells);
}

sub readSparseMatrix($){
	my ($matrixFile) = @_;
	print "Reading sparse matrix $matrixFile";
	my ($I,$status) = gzipopen($matrixFile,"Matrix File",1);
	die "Could not open matrix $matrixFile\n" unless ($status && defined $I);
	my (%rows,@header);
	my ($lineCount,$nonzero,$cells) = (0,0,0);
	while (my $line = <$I>){
		chomp $line;
		$lineCount++;
		my @values = split /\t/,$line,-1;
		my $gene = shift @values;
		if ($lineCount == 1){
			@header = @values;
			next;
		}
		my (@indices,@rowValues);
		for (my $index=0; $index<@values; $index++){
			my $value = $values[$index];
			next if !defined($value) || $value eq "" || $value == 0;
			push(@indices,$index);
			push(@rowValues,$value + 0);
			$nonzero++;
		}
		$cells += scalar(@values);
		$rows{$gene} = [\@indices,\@rowValues];
	}
	close $I;
	die "not enough lines in matrix $matrixFile\n" if $lineCount <= 2;
	my $density = $cells ? 100 * $nonzero / $cells : 0;
	print sprintf(" .. Done (%.2f%% non-zero)\n",$density);
	return {header=>\@header, rows=>\%rows, storage=>"sparse"};
}

sub matrixIsSparse(){
	return exists($FMGmatrix->{storage}) && $FMGmatrix->{storage} eq "sparse";
}

sub matrixRowExists($){
	my ($gene) = @_;
	return exists($FMGmatrix->{rows}{$gene}) if matrixIsSparse();
	return exists($FMGmatrix->{$gene});
}

sub matrixHeader(){
	return $FMGmatrix->{header};
}

sub matrixWidth(){
	return scalar(@{$FMGmatrix->{header}});
}

sub matrixGenes(){
	return keys %{$FMGmatrix->{rows}} if matrixIsSparse();
	return grep { $_ ne "header" } keys %{$FMGmatrix};
}

sub addMatrixRow($$){
	my ($target,$gene) = @_;
	die "Matrix entry missing: $gene\n" unless matrixRowExists($gene);
	my $width = matrixWidth();
	push(@{$target},(0) x ($width - scalar(@{$target}))) if scalar(@{$target}) < $width;
	if (matrixIsSparse()){
		my ($indices,$values) = @{$FMGmatrix->{rows}{$gene}};
		for (my $i=0; $i<@{$indices}; $i++){
			$target->[$indices->[$i]] += $values->[$i];
		}
		return;
	}
	my $row = $FMGmatrix->{$gene};
	for (my $i=0; $i<@{$row}; $i++){
		$target->[$i] += $row->[$i];
	}
}

sub matrixCorrelation($$){
	my ($dense,$gene) = @_;
	die "Matrix entry missing: $gene\n" unless matrixRowExists($gene);
	return correlation($dense,$FMGmatrix->{$gene}) unless matrixIsSparse();
	my ($indices,$values) = @{$FMGmatrix->{rows}{$gene}};
	my $count = scalar(@{$dense});
	my ($sumY,$sumYY) = (0,0);
	foreach my $value (@{$dense}){
		$sumY += $value;
		$sumYY += $value * $value;
	}
	my ($sumX,$sumXX,$sumXY) = (0,0,0);
	for (my $i=0; $i<@{$indices}; $i++){
		my $value = $values->[$i];
		$sumX += $value;
		$sumXX += $value * $value;
		$sumXY += $value * $dense->[$indices->[$i]];
	}
	my $ssxx = $sumXX - $sumX * $sumX / $count;
	my $ssyy = $sumYY - $sumY * $sumY / $count;
	my $ssxy = $sumXY - $sumX * $sumY / $count;
	return sprintf("%.4f",0) if $ssxx == 0 || $ssyy == 0 || $ssxy == 0;
	return sprintf("%.4f",$ssxy / sqrt($ssxx * $ssyy));
}
sub readMotuTax; sub fixGTDBtax;
sub readGene2mlinkage; sub readSpecIids;
sub readNCBI;sub read_speci_tax;
sub MGSassign; sub transferSI2MGS;
sub readMGS;
sub rebase0;

#sub calculate_spearman_correlation;
#sub read_matrix; 

sub getCorrs;
sub add2geneList;
sub resolveAssignments;
sub specImatrix;
#sub createAreadSpecItax; 
sub writeSpecItax;

#my $SpecID="/g/bork3/home/hildebra/DB/MarkerG/specI/"; my $freeze11=1;
my $freeze11=0;
my $doGTDBtax = 1; #Firmicutes_A etc
#progenomes.specIv2_2
my $globalCorrThreshold = 0.6; # determines cutoff, when still to accept correlating genes into specI
my $reblast=0;#do blast again?

my $rarBin = getProgPaths("rare");#"/g/bork5/hildebra/dev/C++/rare/rare";
my $samBin = getProgPaths("samtools");#"/g/bork5/hildebra/bin/samtools-1.2/samtools";
my $bts = getProgPaths("buildTree_scr");



if (@ARGV == 0){
	die "Not enough input args: use ./annotateMGwMotus.pl [path to GC] [# Cores]\n";
}

my $GCd ="";#$ARGV[0]."/";
my $mode = "specI";
my $BlastCores = 1;#$ARGV[1];
my $MGSfile = ""; #MGS annotations; will take precedence of blast specI annotations..
#$MGSfile = $ARGV[2] if (@ARGV > 2);
my $MGStax = ""; #MGS taxonomy (fitting to MGSfile).. only if given will try to calc abundance table
#$MGStax = $ARGV[3] if (@ARGV > 3);
my $useGTDBmg = "FMG"; #GTDB or FMG
my $tmpD = "";
my $minGenes = 10; #that many MGs are required, to include a species..
my $outD = "";
$matrixStorage = "auto"; # auto selects sparse rows only below a conservative density threshold
my $hr1; #general purpose hash ref pointer..




#options to pipeline..
GetOptions(
	"GCd=s"      => \$GCd,
	"outD=s"     => \$outD,
	"tmp=s"      => \$tmpD,
	"MGS=s"      => \$MGSfile,
	"cores=i"    => \$BlastCores,
	"MGStax=s"   => \$MGStax,
	"MGset=s"    => \$useGTDBmg,#GTDB or FMG
	"minGenes=i" => \$minGenes,
	"matrixStorage=s" => \$matrixStorage,
) or die "Invalid annotateMGwSpecIs2.pl options\n";
die "Unexpected positional arguments: @ARGV\n" if @ARGV;
die "Needs option -GCd $GCd\n" if ($GCd eq "");
die "Gene-catalog directory not found: $GCd\n" unless -d $GCd;
die "-cores and -minGenes must be positive\n" unless $BlastCores > 0 && $minGenes > 0;
die "-matrixStorage must be auto, dense, or sparse\n" unless $matrixStorage =~ /^(auto|dense|sparse)$/;


my $speciesLink = "specI_lnks"; my $speciesCutoff = "specI_cutoff";
my $speciesGTDB = "specI_GTDB"; my $speciesDir = "specIPath";
my $MGterm = "FMG";
#GTDB marker genes instead of FMG?? 
die "-MGset option has to be \"GTDB\" or \"FMG\"\n" unless ($useGTDBmg eq "GTDB" || $useGTDBmg eq "FMG");
if ($useGTDBmg eq "GTDB"){ 
	$MGterm = "GTDBmg";
	$speciesLink = "GTDB_lnks"; $speciesCutoff = "GTDB_cutoff";
	$speciesGTDB = "GTDB_GTDB"; $speciesDir = "GTDBPath";
}
my $MGdir = "$GCd/$MGterm/";
if ($outD eq "") {
	$outD = $MGSfile ne "" ? "$GCd/Anno/Tax/${MGterm}_MGS/" : "$GCd/Anno/Tax/$MGterm/";
}
my $inSImap = getProgPaths($speciesLink);
my $GTDBspecI = getProgPaths($speciesGTDB);
my $SpecID=getProgPaths($speciesDir);#directoy with all 40 SpecI marker genes


print "-------------------------------------------------------------------------\nannotateMG script v $version\n-------------------------------------------------------------------------\n";
print "Using gene cat in $GCd\nMGS: $MGSfile\nMGStax: $MGStax\ncores: $BlastCores MGset: $useGTDBmg minGenes: $minGenes matrixStorage: $matrixStorage\n";


#ystem "mkdir -p $MGdir/tax" unless (-d "$MGdir/tax");
system "mkdir -p $outD" unless (-d "$outD");

#die "Writing to $outD, using $MGdir\n";

#die $outD;
my $motuDir = "";#"/g/bork3/home/hildebra/DB/MarkerG/mOTU";

#annotate against DB using lambda
#my $hr1 = readTable($GTDBspecI,"\t",";" );
#my %specItax = %{$hr1};


#first check if taxpergene exists (can be used in LCA algo directly)
if (0 ){ #not needed: there should be *.LCA file now!
	my $taxPerGene = "$SpecID/specI.pergene.tax";
	my $prepPG = getProgPaths("progenomes_prep_scr");
	my $cmd = "$prepPG $SpecID $inSImap $GTDBspecI $taxPerGene\n";
	system $cmd if (!-e $taxPerGene);
}


#specItax is essentially only used to register what key is present..
#my $specIfullTaxHR = createAreadSpecItax(\%specItax,"$SpecID/specI.tax3",$GTDBspecI); #specI.tax3 doesn't need to exist any longer..
#
#tax per specI - new way relying on precomputed LCAs
my %specIfullTax = %{readTable($GTDBspecI,"\t" )};# %{$specIfullTaxHR};
my $xtrLab= "";$xtrLab= ".rep" if ($freeze11);

#assign each COG separately
system "mkdir -p $MGdir" unless (-d $MGdir);

my %SpecIgenes; my %gene2COG;

my %Q2S; #per-gene assignment map, compacted to counts before correlation

my %COGs;
#my %FMGlist;
my $COGgenes=0;
open IC,"<$GCd/${MGterm}.subset.cats" or die "Can't open $GCd/${MGterm}.subset.cats\n";
while (<IC>){
	chomp;	my @spl  = split /\t/;
	#$cats{$spl[0]} = $spl[2];
	my @genes = split(/,/,$spl[2]);
	my $curCOG = $spl[0];
	$COGs{$curCOG} = 1;
	foreach (@genes){$gene2COG{$_} = $curCOG; $COGgenes++;}
}
close IC;
print "Read $COGgenes COG 2 MGs genes\n";




#MGS related containers
my %Gene2MGS; my %speci2MGS; 
my %MGSlist; #$MGSlist{$curMGS} = 1;

readMGS($MGSfile); #only used if in MGS mode..
#undef %FMGlist; #no longer needed from here on
#die;
my $allOK=1;

#new routine: will read .LCA to get links..
print "Reading LCA's\n";
foreach my $COG (keys %COGs){
	# Keep only one LCA table in memory at a time. Copying it into %tax
	# temporarily doubled the largest per-COG allocation.
	my $tax = readTable("$MGdir/$COG.LCA","\t",";",1);
	foreach my $gen (keys %{$tax}){
		my $speci = $tax->{$gen}; #complete tax string.. well should be ok
		if (!exists($specIfullTax{$speci})){
			my @tmp = split /;/,$speci;
			$specIfullTax{$speci} = \@tmp;
		}
		# Just link the specI to MGS.
		if (exists($Gene2MGS{$gen})){
			my $curMGS = $Gene2MGS{$gen};
			$speci2MGS{$curMGS}{$speci}++;
		}
		# Preserve all initial assignments so getCorrs can exclude ambiguous genes.
		if (!exists($Q2S{$gen}{$speci})){
			$Q2S{$gen}{$speci} = scalar(keys(%{$Q2S{$gen}}));
		}
		push(@{$SpecIgenes{$speci}{$COG}},$gen);
	}
}

#die;

if ($allOK==0){die"One or more Blast files were not ok.. restart procedure\n";}
# Write out SpecI -> MGS assignments (and tax), then release the vote table.
my $SI2MGS = MGSassign();
# Transfer gene assignments to MGS, transfer SI tax -> MGS tax.
transferSI2MGS($SI2MGS);
undef %speci2MGS;
# After transfer, downstream code only needs to know that a gene belongs to an MGS.
# Release the MGS identifier strings while retaining that membership flag.
$Gene2MGS{$_} = 1 for keys %Gene2MGS;
# getCorrs only needs the number of distinct assignments per gene.
# Replace each nested per-species hash before the large matrix is resident.
for my $gene (keys %Q2S) {
	$Q2S{$gene} = scalar(keys %{$Q2S{$gene}});
}

print "Reading MG marker gene matrix..\n";
$FMGmatrix = loadMarkerMatrix("$GCd/Matrix.$MGterm.mat");

my %gene2specI; my %specIprofiles;
my %SpecIgenes2; #collects list of genes that could be associated to specI
# Sort out best multi hit by correlation analysis.
getCorrs();
# Subsequent passes use only the selected assignments in %gene2specI.
undef %SpecIgenes;
undef %Q2S;

# Build temporary profiles, then replace the mutable multi-assignment state with
# one deterministic owner per gene and one gene per owner/COG slot.
rebase0();
resolveAssignments();
rebase0(1);
print "resolving remainder genes..\n";

# getCorrs no longer needs the raw LCA assignment index after it returns.
my $xtraEntry=0;
foreach my $COG (keys %COGs){
	# Stream one LCA table at a time; do not create a duplicate hash.
	my $tax = readTable("$MGdir/$COG.LCA","\t",";",1);
	foreach my $gid (keys %{$tax}){
		next if exists($Gene2MGS{$gid});
		my $speci = $tax->{$gid}; #complete tax string.. well should be ok
		$speci = $SI2MGS->{$speci} if exists($SI2MGS->{$speci});
		# This MG has already been assigned in the high-confidence initial pass.
		next if exists($SpecIgenes2{$speci}{$COG});

		# Correlate to the species core before accepting the remainder assignment.
		next unless exists($specIprofiles{$speci});
		my $corr = matrixCorrelation($specIprofiles{$speci},$gid);
		next if $corr < $globalCorrThreshold;

		$xtraEntry++;
		# Block this COG slot for the species and record the final assignment.
		add2geneList($speci,$COG,$gid);
	}
}
undef $SI2MGS;
undef %SpecIgenes2;



print "Entries not in MGS: $xtraEntry\n";


#write specI assignments for markerG
open O,">$MGdir/gene2specI.txt";
foreach my $k (keys %gene2specI){
	if (!exists($gene2COG{$k})){
		print "No COG assignment $k!  ";
	} else {
		print O "$k\t". join(",",keys %{$gene2specI{$k}}) . "\t$gene2COG{$k}\n";
	}
}
close O;

#create abundance profile
rebase0(1);#create final specIprofile..
specImatrix("$outD/specI.mat",\%specIfullTax);


print "Finished SpecI annotations & matrix: $outD\n";
exit(0);




#####################################################
#####################################################
#####################################################
#####################################################



sub readMGS{
	my ($MGSfile) = @_;
	return if ($MGSfile eq "");
	print "Reading reference MGS: $MGSfile\nTax: $MGStax\n";
	my %tmpTax;
	unless ($MGStax eq ""){
		open IT,"<$MGStax" or die "can't open MGS tax file $MGStax\n";
		while (<IT>){
			chomp;my @spl=split/\t/;
			next if ($spl[0] eq "domain" || $spl[0] eq "user_genome");
			my $id = shift @spl;
			if (@spl == 1){
				my $tmp=$spl[0]; $tmp =~ s/;;/;\?;/;$tmp =~ s/;$/;\?/;
				@spl = split /;/,$tmp;
				push(@spl,"?") while (@spl < 7);
			}
			#@spl = fixGTDBtax(@spl) if ($doGTDBtax);
			$tmpTax{$id} = \@spl;
		}
		close IT;
	}
	my %markerOwners;
	open IM,"<$MGSfile" or die "Can't open MGS $MGSfile\n";
	while (<IM>){ 
		chomp;
		next if /^\s*$/ || /^#/;
		my @spl = split /\t/, $_, -1;
		next if $spl[0] eq "Bin";
		die "Malformed MGS row in $MGSfile at line $.: $_\n"
			unless @spl >= 2 && length($spl[0]) && length($spl[1]);
		my @spl2 = split /,/,$spl[1];
		my $curMGS = $spl[0];
		$MGSlist{$curMGS} = 1;
		foreach my $gen (@spl2){
			next unless (exists($gene2COG{$gen}));
			$markerOwners{$gen}{$curMGS} = 1;
		}
	}
	close IM;

	my $ambiguousMarkers = 0;
	foreach my $gen (sort keys %markerOwners) {
		my @owners = sort keys %{$markerOwners{$gen}};
		if (@owners != 1) {
			$ambiguousMarkers++;
			next;
		}
		my $curMGS = $owners[0];
		$Gene2MGS{$gen} = $curMGS;
		push(@{$SpecIgenes{$curMGS}{$gene2COG{$gen}}},$gen);
	}
	print "Excluded $ambiguousMarkers marker genes assigned to multiple MGS\n"
		if $ambiguousMarkers;

	my $taxF=0; my $taxN=0;
	foreach my $curMGS (sort keys %MGSlist) {
		if (exists($tmpTax{$curMGS})){
			$specIfullTax{$curMGS} = $tmpTax{$curMGS};
			$taxF++;
		} else {
			$specIfullTax{$curMGS} = ["Bins","?","?","?","?","?","?"]
				unless exists($specIfullTax{$curMGS});
			$taxN++;
		}
	}
	my $meanS=0; my @sizes; my $morethan1=0; my $only1=0;
	foreach my $MGS (sort keys %MGSlist){
		my $curS =0;
		foreach my $cat (sort keys %{$SpecIgenes{$MGS} || {}}){
			 my $lcurS = @{$SpecIgenes{$MGS}{$cat}};
			 $curS += $lcurS;
			 if ($lcurS > 1){$morethan1 ++ ;} else {$only1++;}
		}
		$meanS += $curS;
		push(@sizes, $curS);
	}
	@sizes = sort { $a <=> $b } @sizes;
	my $tax_total = $taxF + $taxN;
	my $copy_total = $only1 + $morethan1;
	my $mean_markers = $tax_total ? $meanS / $tax_total : 0;
	my $median_markers = @sizes ? median(@sizes) : 0;
	print "Found ". scalar(keys %MGSlist)." MGS with $taxF/$tax_total taxonomies. Mean MGs/MGS: $mean_markers; median: $median_markers. $morethan1/$copy_total with >1 copy.\n";
#die "@sizes\n";
	
	#print"T::@{$specIfullTax{MGS0287}}\n";

}


sub MGSassign{
	my %Si2MGS;
	return \%Si2MGS if ($MGSfile eq "");
	my $logfile = "$outD/MGS2speci.txt";
	my $collisionLog = "$outD/MGS2speci.collisions.tsv";
	open OM,">$logfile" or die "Can't open $logfile\n";
	open OC,">$collisionLog" or die "Can't open $collisionLog\n";
	print OC "specI\tMGSs\taction\n";

	my (@records, @maxAssi);
	my %ownersForSI;
	my $notAssigned=0;
	foreach my $MGS (sort keys %speci2MGS){
		my @SIsAdd;
		my $valSI=""; my $valStr=0;
		my $totalAssi =0;
		my $maxassFrac = 0;
		foreach my $si (keys %{$speci2MGS{$MGS}}){$totalAssi += $speci2MGS{$MGS}{$si}; }
		my @ranked_si = sort {
			$speci2MGS{$MGS}{$b} <=> $speci2MGS{$MGS}{$a} || $a cmp $b
		} keys %{$speci2MGS{$MGS}};
		if (@ranked_si && $speci2MGS{$MGS}{$ranked_si[0]} > 5) {
			$valSI = $ranked_si[0];
			$valStr = $speci2MGS{$MGS}{$valSI};
			$valSI = "" if @ranked_si > 1 && $speci2MGS{$MGS}{$ranked_si[1]} > 0.5 * $valStr;
		}
		foreach my $si (@ranked_si){
			push(@SIsAdd,"$si:$speci2MGS{$MGS}{$si}");
			my $assFrac = $totalAssi ? $speci2MGS{$MGS}{$si}/$totalAssi : 0;
			$maxassFrac = $assFrac if ($assFrac >= $maxassFrac);
		}
		push(@maxAssi,$maxassFrac);
		$ownersForSI{$valSI}{$MGS} = 1 if $valSI ne "";
		my $tax = $valSI ne "" && exists($specIfullTax{$valSI})
			? join(";",@{$specIfullTax{$valSI}}) : "?;?;?;?;?;?;?";
		push(@records, {mgs=>$MGS, si=>$valSI, tax=>$tax, votes=>join(",",@SIsAdd)});
	}

	my %collidingSI;
	foreach my $si (sort keys %ownersForSI){
		my @owners = sort keys %{$ownersForSI{$si}};
		next if @owners == 1;
		$collidingSI{$si} = 1;
		print OC "$si\t".join(",",@owners)."\twithheld: multiple MGS selected the same specI\n";
	}
	close OC;

	my ($assignedMGS, $withheldMGS) = (0, 0);
	foreach my $record (@records){
		print OM "$record->{mgs}\t$record->{tax}\t$record->{votes}\n";
		if ($record->{si} eq "") {
			$notAssigned++;
		} elsif ($collidingSI{$record->{si}}) {
			$withheldMGS++;
		} else {
			$Si2MGS{$record->{si}} = $record->{mgs};
			$assignedMGS++;
		}
	}
	close OM;

	if ($MGStax eq ""){
		warn "No MGStax given; continuing with placeholder MGS taxonomy where necessary.\n";
	}

	my $chimera=0;
	foreach (@maxAssi){
		$chimera ++ if ($_ <0.6);
	}
	my $totMGS = $notAssigned+$assignedMGS+$withheldMGS;
	print "Matched ${assignedMGS} of ". ($totMGS) ." MGS to taxa/specIs";
	print "; withheld $withheldMGS colliding MGS assignments" if $withheldMGS;
	print ": $logfile\n";
	my $mean_assignment = @maxAssi ? mean(@maxAssi) : 0;
	my $median_assignment = @maxAssi ? median(@maxAssi) : 0;
	print "Mean max assignments MGS: $mean_assignment , median: $median_assignment; pot. chimeric (species level and above): $chimera \n";
	return \%Si2MGS;
}


sub transferSI2MGS{
	my ($hr) = @_;
	my $Si2MGS = $hr;
	print "Comparing MGS to SpecIs..\n" if scalar(keys %{$Si2MGS});
	foreach my $valSI (keys %{$Si2MGS}){
		my $MGS = $Si2MGS->{$valSI};
		#print "$valSI  $MGS\n";
		#get tax transferred..
		$specIfullTax{$MGS} = $specIfullTax{$valSI} unless (exists( $specIfullTax{$MGS} ));
		#actually replace SI with this MGS
		foreach my $COG (keys %{$SpecIgenes{$valSI}}){
			#transfer Q2S, delete old entry for it...
			foreach my $gen (@{$SpecIgenes{$valSI}{$COG}}){
				$Q2S{$gen}{$MGS} = $Q2S{$gen}{$valSI};
				delete $Q2S{$gen}{$valSI};
			}
			my @adds= @{$SpecIgenes{$valSI}{$COG}};
			if (exists($SpecIgenes{$MGS}{$COG})){ 
				push(@adds,@{$SpecIgenes{$MGS}{$COG}});
				@{$SpecIgenes{$MGS}{$COG}} = do { my %seen; grep { !$seen{$_}++ } @adds };
			} else { @{$SpecIgenes{$MGS}{$COG}} = @adds;}
			delete $SpecIgenes{$valSI}{$COG};
		}
		delete $SpecIgenes{$valSI};
	}
}


sub rebase($){ # calculates the profile for each SI from its marker-gene rows
	my ($hr1)=@_;
	my $specIs = $hr1;
	my $skippedSpecies=0;
	foreach my $sp (keys %{$specIs}){
		my @tar;
		my $MGn=0;
		foreach my $gid (@{$specIs->{$sp}}){
			addMatrixRow(\@tar,$gid);
			$MGn++;
		}
		if ($MGn == 0){
			$skippedSpecies++;
			next;
		}
		for (my $j=0;$j<@tar;$j++){$tar[$j] /= $MGn;}
		$specIprofiles{$sp} = \@tar;
	}
	print "Rebbase: skipped $skippedSpecies species/MGS\n";
}
sub rebase0($){
	my ($replaceProfiles) = @_;
	my %specIs;
	foreach my $k (keys %gene2specI){
		my @spl = sort keys %{$gene2specI{$k}};
		die "gene not in gene2specI $k\n" if (!exists($gene2specI{$k}) || @spl == 0);
		push(@{$specIs{$spl[0]}},$k);
	}
	%specIprofiles = () if $replaceProfiles;
	rebase(\%specIs);
}

sub readMotuTax($){
	my ($inF) = @_;
	my %gene2motu;
	my %motu2tax;
	open I,"<$inF";
	while (my $l = <I>){
		chomp $l; my @spl=split/\t/,$l;
		$gene2motu{$spl[0]} = $spl[8];
		$motu2tax{$spl[8]} = join(";",@spl[1,2,3,4,5,6,7]) if (!exists $motu2tax{$spl[8]} );
	}
	close I;
	return (\%gene2motu,\%motu2tax);
}

sub readGene2mlinkage($){
	my ($inF) = @_;
	my %gene2LG;
	open I,"<$inF";
	while (my $l = <I>){
		chomp $l; my @spl=split/\t/,$l;
		$gene2LG{$spl[0]} = $spl[2];
	}
	close I;
	return (\%gene2LG);
}
sub read_speci_tax($){
	my ($inF) = @_;
	my %gene2LG;
	open I,"<$inF";
	while (my $l = <I>){
		chomp $l; my @spl=split/\t/,$l;
		$gene2LG{$spl[0]} = $spl[2];
	}
	close I;
	return (\%gene2LG);
}

sub readNCBI($){
	my ($inF) = @_;
	my %gene2LG;
	open I,"<$inF";
	while (my $l = <I>){
		chomp $l; my @spl=split/\t/,$l;
		$gene2LG{$spl[0]} = $spl[2];
	}
	close I;
	return (\%gene2LG);
}



sub specImatrix($ $){
	my ($oF,$hr) = @_;
	print "Creating specI matrix..\n";
	my $sTax = $hr;
	#print "@{$sTax{specI_v2_Cluster34}}\n";
	
	my %specIcnts; #just use for histgram
	foreach my $k (keys %gene2specI){$specIcnts{ (keys %{$gene2specI{$k}})[0] }++;}
	
	my $rmSpecs=0;my %delSIs;
	foreach my $si (sort keys %specIprofiles){
		if (exists($specIcnts{$si}) && $specIcnts{$si} < $minGenes){
			delete $specIprofiles{$si};
			$delSIs{$si} = 1;
			$rmSpecs++;
		}
	}
	print "Removed $rmSpecs specI's from final matrix due to having <$minGenes marker genes\n";
	
	#create background count of SpecI genes not assigned
	my %bkgrnd;
	foreach my $gid (matrixGenes()){
		if (exists ($gene2specI{$gid})  ){
			if ( exists( $delSIs{  (keys %{$gene2specI{$gid}})[0] }  )    ){
				delete $gene2specI{$gid};
			} else {
				next;
			}
		}
		next if exists($Gene2MGS{$gid});
		addMatrixRow($bkgrnd{$gene2COG{$gid}},$gid);
	}
	#split up background by COG/sample; the former dblCh vector was never used.
	my @bkgrnd1;
	for (my $j=0;$j<matrixWidth();$j++){
		my @meanVal;
		foreach my $cog (keys %bkgrnd){
			if ($bkgrnd{$cog}[$j] > 0){
				push(@meanVal, $bkgrnd{$cog}[$j]);
			}
		}
		$bkgrnd1[$j] = mean(@meanVal);
	}
	
	#print "@{$specIprofiles{specI_v2_Cluster34}}\n";
	
	open Ox,">$oF" or die "Can't open out mat $oF\n";
	print Ox "SpecI\t".join ("\t",@{matrixHeader()})."\n";
	#print O "SUM\t\t".join ("\t",@dblCh)."\n";
	print Ox "?\t".join ("\t",@bkgrnd1)."\n";
	foreach my $si (sort keys %specIprofiles){
		print Ox "$si\t". join ("\t",@{$specIprofiles{ $si }})."\n";
	}
	close Ox;
	my $oFx = $oF; $oFx =~ s/\.[^\.]*$//;
	

		

	#and print SI tax for later reference..
	open Ot,">$oFx.tax" or die "Can;t open tax file $oFx.tax\n";
	foreach my $si (sort keys %specIprofiles){
		if (!exists($sTax->{$si})){
			if (exists($MGSlist{$si})){
			} else {die "doesnt exist: $si\n";	}
		}
		print Ot "$si\t".join ("\t",@{$sTax->{$si}})."\n";
	}
	close Ot;
	
	#
	#
	#calculating higher level abundance matrix
	my @taxLs = ("superkingdom","phylum","class","order","family","genus","species");
	#print "@{$specIprofiles{specI_v2_Cluster34}}\n";
	for (my $t=0;$t<@taxLs;$t++){
		#sum up to hi lvl
		my %thisMap;
		foreach my $si (sort keys %specIprofiles){
			
			if (!exists($sTax->{$si})){
				if (exists($MGSlist{$si})){
				} else {
					die "doesnt exist: $si\n";
				}
			}
			print "ERR: $si    @{$sTax->{$si}}\n" if (@{$sTax->{$si}} <= $t);
			my $clvl = join (";",@{$sTax->{$si}}[0 .. $t]);
			
			#if ($clvl eq "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia albertii"){
				#print "$si\n@{$specIprofiles{ $si }}\n";
			#}
			
			if (exists($thisMap{$clvl})){
				for (my $kl=0;$kl<scalar(@{$specIprofiles{ $si }});$kl++){
					${$thisMap{$clvl}}[$kl] += ${$specIprofiles{ $si }}[$kl];
				}
			} else {
				$thisMap{$clvl} = [@{$specIprofiles{ $si }}];
			}
		}
		open Ot,">$oFx.$taxLs[$t]" or die "Can't openn $oFx.$taxLs[$t]\n" ;
		print Ot "$taxLs[$t]\t".join ("\t",@{matrixHeader()})."\n";
		print Ot "?\t".join ("\t",@bkgrnd1)."\n";
		foreach my $kk (sort keys %thisMap){
			chomp $kk;
			print Ot "$kk\t".join("\t",@{$thisMap{$kk}}) . "\n";
		}
		close Ot;
	}
	#print "@{$specIprofiles{specI_v2_Cluster34}}\n";
	
	
	
	#some stats on created abundances.. 
	%specIcnts = ();
	foreach my $k (keys %gene2specI){$specIcnts{ (keys %{$gene2specI{$k}})[0] }++;}
	my %histo;
	for my $k (sort {$specIcnts{$a} <=> $specIcnts{$b}} keys %specIcnts) {
		$histo{int $specIcnts{$k}/10}++;
	   # print "$k $specIcnts{$k} $NTax{$specItax{$k}}\n" ;#if ($specIcnts{$k}>=40);   # bbb c aaaa
	}
	foreach (sort{$a <=> $b} (keys %histo)){
		print "$_\t$histo{$_}\n";
	}


}

sub writeSpecItax{
	my ($hr1,$file) = @_;#\$specIid,"$SpecID/specI.tax2");
	print "Writing newly created  SpeciI tax to $file\n";
	my %sNTID = %{$hr1};
	open O,">$file" or die "can't open $file\n";
	foreach my $k (keys %sNTID){
		print O $k."\t".join("\t",@{$sNTID{$k}})."\n";
	}
	close O;
	die;
}

sub readSpecIids($){
	my ($inSImap ) = @_;
	#specItax is essentially only used to register what key is present..
	my %specIid;my %specItax;
	open I,"<$inSImap" or die "Can't open in specI map $inSImap\n";
	#specI v3 (progenomes2)
	#while (<I>){next if (m/^#/);chomp; my @xx = split /\t/;$xx[1] =~s/,//g; $specIid{$xx[1]} = $xx[0];$xx[1]=~m/^(\d+)\./; $specItax{$xx[0]}=$1;}
	#specI v4 (progenomes3)
	while (<I>){next if (m/^#/);chomp; 
		my @xx = split /\t/;
		foreach my $yy (split /;/,$xx[1]){
		$yy =~s/,//g; $specIid{$yy} = $xx[0];
		$xx[1]=~m/^(\d+)\./; $specItax{$xx[0]}=$1;
		}
	}
	close I;
	return (\%specIid, \%specItax);
}




#check which mean abundance profile single MG have, and how multi MGs correlate to this 
#then selects best correlating MG to be "the one" that just fits
sub getCorrs{
	my $dblAssi=0;
	my $dblA=0;my $singlA=0;my $singlMultA=0;my $skippedSIs=0; my $newAssigns=0;
	my $belowGeneIncl = 0;
	print "Correlations of MGs..\n";
	foreach my $k (keys %SpecIgenes){ #this is specI
		print "$k ; " if ($k =~ m/,/);
		my @tarGenes; #matrix vector summed over single MGs
		my $gcnt =0; my %selV;
		
		#commented, as gene needs to be stored in $SpecIgenes2
		#next if (exists($MGSlist{$k})); #this is an MGS and doesn't have Q2S etc and should be static in any case..
		
		foreach my $cog(keys %{$SpecIgenes{$k}}){ #this is COG
			$selV{$cog} =0;
			if (@{$SpecIgenes{$k}{$cog}}>1){next;}#multi copy, dont use this gene
			my $gid = ${$SpecIgenes{$k}{$cog}}[0];
			if (($Q2S{$gid} // 0) > 1){$dblA++;next;}
			$singlA++; 
			#print " $gid ";
			addMatrixRow(\@tarGenes,$gid);
			$gcnt++;				
			$dblAssi++ if (add2geneList($k,$cog,$gid));
			$selV{$cog} =1;
		}
		#print "T2: $gcnt   " . scalar(keys(%selV))." D" if ($k eq "TEC2");
		
		
		if ($gcnt ==0){$skippedSIs++;next;}	
		if (nonZero(\@tarGenes) < 3){ next;}#print "XX" ;
		#if ($gcnt < $minGenes){$belowGeneIncl++;next;}
		
		#doesn't need norm, since we do spearman correlation
		#but now corr and see which genes just fit best of the multi choices..
		#print "@tarGenes\n";
		foreach my $c(keys %{$SpecIgenes{$k}}){
			if (!$selV{$c}){ #$SpecIgenes{$k}{$c} =~ m/,/){#only look at multi assigned genes
				#die "$SpecIgenes{$k}{$c}\n";
				my @spl = @{$SpecIgenes{$k}{$c}};
				my @subCors; my $max=0;
				foreach my $sg (@spl){
					die "$sg doesn't exist in gene list \n" unless matrixRowExists($sg);
					my $corr = matrixCorrelation(\@tarGenes,$sg);
					if ($corr > $max){$max = $corr;}
					push (@subCors,  $corr  );
				}
				if ($max < $globalCorrThreshold){next;}
				$newAssigns++;
				for (my $i=0;$i<@spl;$i++){
					unless ($subCors[$i]>$max-0.03){next;}
					my $sg = $spl[$i];
					#this assignment is what I need, now I know this gene is blocked for assignment to other SpecI's
					$dblAssi++ if (add2geneList($k,$c,$sg));
					addMatrixRow(\@tarGenes,$sg);
					$gcnt++;
				}
			}
		}

		#norm vector
		for (my $j=0;$j<@tarGenes;$j++){$tarGenes[$j] /= $gcnt;}
		#and save the final specI profile..
		$specIprofiles{$k} = \@tarGenes;
		#print @tarGenes ."XX\n";
	}
	print "double assignment $dblAssi; assigned: $newAssigns   Stats in Run: mult. spec. $dblA $singlA $singlMultA $skippedSIs $belowGeneIncl\n";
	
	
}

sub resolveAssignments{
	my $auditFile = "$outD/assignment_resolution.tsv";
	open AR,">$auditFile" or die "Can't open $auditFile\n";
	print AR "gene\tCOG\tcandidates\tselected_owner\tselected_score\treason\n";

	# Rebuild both assignment indexes from one canonical set.  The old
	# SpecIgenes2 values can disagree with gene2specI after multi-hit handling.
	undef %SpecIgenes2;
	my %bestForSlot;
	my ($multiOwner, $missingProfile, $discarded, $duplicateCOG) = (0, 0, 0, 0);
	foreach my $gene (keys %gene2specI){
		my $owners = $gene2specI{$gene};
		delete $gene2specI{$gene}; # release the old nested map as the new one is built
		die "Can't find gene $gene in FMG matrix!\n" unless matrixRowExists($gene);

		my (@candidates, @auditCandidates);
		foreach my $owner (sort keys %{$owners}){
			unless (exists($specIprofiles{$owner})){
				push(@auditCandidates, "$owner:NA(profile_missing)");
				$missingProfile++;
				next;
			}
			my $score = matrixCorrelation($specIprofiles{$owner},$gene);
			push(@candidates, [$owner, $score]);
			push(@auditCandidates, "$owner:".sprintf("%.4f",$score));
		}
		if (!@candidates){
			my $cog = $gene2COG{$gene} // "";
			print AR "$gene\t$cog\t".join(";",@auditCandidates)."\t\t\tno_profile_candidate\n";
			$discarded++;
			next;
		}

		my ($bestOwner, $bestScore);
		foreach my $candidate (@candidates){
			my ($owner, $score) = @{$candidate};
			if (!defined($bestOwner) || $score > $bestScore ||
				($score == $bestScore &&
					((exists($MGSlist{$owner}) > exists($MGSlist{$bestOwner})) ||
					 (exists($MGSlist{$owner}) == exists($MGSlist{$bestOwner}) && ($owner cmp $bestOwner) < 0)))){
				($bestOwner, $bestScore) = ($owner, $score);
			}
		}
		my $reason = @candidates > 1 ? "multiple_candidates" : "single_candidate";
		$reason .= "+profileless_ignored" if @auditCandidates > @candidates;
		$multiOwner++ if @candidates > 1;
		if ($reason ne "single_candidate"){
			my $cog = $gene2COG{$gene} // "";
			print AR "$gene\t$cog\t".join(";",@auditCandidates)."\t$bestOwner\t".sprintf("%.4f",$bestScore)."\t$reason\n";
		}

		my $cog = $gene2COG{$gene} // "";
		my $slot = join("\x1e",$bestOwner,$cog);
		my $current = $bestForSlot{$slot};
		if (!defined($current) || $bestScore > $current->[2] ||
			($bestScore == $current->[2] && ($gene cmp $current->[0]) < 0)){
			if (defined($current)){
				print AR "$current->[0]\t$cog\t$current->[1]:".sprintf("%.4f",$current->[2])."\t\t\tduplicate_cog_replaced_by_$gene\n";
				$duplicateCOG++;
			}
			$bestForSlot{$slot} = [$gene,$bestOwner,$bestScore];
		} else {
			print AR "$gene\t$cog\t$bestOwner:".sprintf("%.4f",$bestScore)."\t\t\tduplicate_cog_discarded_by_$current->[0]\n";
			$duplicateCOG++;
		}
	}
	foreach my $slot (keys %bestForSlot){
		my ($owner,$cog) = split(/\x1e/,$slot,2);
		my $gene = $bestForSlot{$slot}[0];
		$gene2specI{$gene} = {$owner => 1};
		$SpecIgenes2{$owner}{$cog} = $gene;
	}
	close AR;
	print "Resolved $multiOwner multi-owner genes; ignored $missingProfile profileless candidates; discarded $discarded genes without a profiled owner; collapsed $duplicateCOG duplicate owner/COG assignments: $auditFile\n";
}


sub tree4FMGs{
	my $btout = "$GCd/${MGterm}/specIphylo/";
	system "mkdir -p $btout" unless (-d $btout);
	my $hr = readFasta("$GCd/${MGterm}/*.faa"); my %FAA = %{$hr};
	$hr = readFasta("$GCd/${MGterm}/*.fna"); my %FNA = %{$hr};
	my $SaSe = "|";
	open ON,">$btout/all.fna"; open OA,">$btout/all.faa"; 
	my %catT;
	foreach my $SI (keys %SpecIgenes){
		my @CGs = keys %{$SpecIgenes{$SI}};
		foreach my $cg (@CGs){
			my $gID = ${$SpecIgenes{$SI}{$cg}}[0];
			next unless (exists($gene2specI{$gID})); #just make sure this specI is also in the latest FMG counting...
			$gID =~ m/^([^,]+)/; $gID = $1;
			die "@{$SpecIgenes{$SI}{$cg}}\n" unless (exists($FNA{$gID}));
			print ON ">$SI$SaSe$cg\n$FNA{$gID}\n";
			print OA ">$SI$SaSe$cg\n$FAA{$gID}\n";
			push(@{$catT{$cg}},"$SI$SaSe$cg");
		}
	}
	close ON; close OA;
	open OC,">$btout/all.cats";
	foreach my $cg (keys %catT){
		print OC join("\t",@{$catT{$cg}})."\n";
	}
	close OC;
	print "Creating phylogeny for found specI's//\n";
	my $cmd= "$bts  -aa  $btout/all.faa -smplSep '\\$SaSe' -cats $btout/all.cats -outD $btout -runIQtree 1 -runFastTree 0 -runRaxMLng 1 -cores $BlastCores  -AAtree 1 -bootstrap 000 -NTfiltCount 300 -NTfilt 0.1 -NTfiltPerGene 0.5 -minOverlapMSA 2 -MSAprogram 2 -AutoModel 0 -iqFast 1 \n";
	systemW "$cmd\n";
	my $QSBoptHR = emptyQsubOpt(1,"");
	#my ($dep,$qcmd) = qsubSystem($btout."treeCmd.sh",$cmd,$BlastCores,"1G","FMGtree","","",1,[],$QSBoptHR);
}


sub add2geneList($ $ $){ #assign a gene ($sg) to a speci($k), and its COG ($c)
	my ($k,$c,$sg) = @_;
	my $ret=0;
	if (exists($gene2specI{$sg})){
		 $ret=1;#print "should not happen\n";#$gene2specI{$spl[$i]}   $k\n";
	} 
	$gene2specI{$sg}{$k} = 1;
	$SpecIgenes2{$k}{$c}=$sg;#set mark to block this MG in this specI...
	return $ret;
}


