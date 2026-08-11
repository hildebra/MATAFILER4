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

use Mods::GenoMetaAss qw(median mean);
use Mods::TamocFunc qw(readTable);
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Getopt::Long qw(GetOptions);
use File::Path qw(make_path);

#.16 delegates constrained sparse profile assignment to cc.bin --speci-assign
my $version = 0.16;

sub MGSassign; sub transferSI2MGS;
sub readMGS;
sub shellQuote;
sub writeSpeciAssignHandoffs;
sub runSpeciAssign;
sub importSpeciAssignments;
sub writeTaxonomyMatrices;

my $globalCorrThreshold = 0.6;
my $nearBestWindow = 0.03;
my $minGuideNonzero = 3;
my $matrixStorage;

if (@ARGV == 0){
	die "Not enough input args: use ./annotateMGwSpecIs3.pl -GCd <path> -cores <n> [options]\n";
}

my $GCd ="";#$ARGV[0]."/";
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
) or die "Invalid annotateMGwSpecIs3.pl options\n";
die "Unexpected positional arguments: @ARGV\n" if @ARGV;
die "Needs option -GCd $GCd\n" if ($GCd eq "");
die "Gene-catalog directory not found: $GCd\n" unless -d $GCd;
die "-cores and -minGenes must be positive\n" unless $BlastCores > 0 && $minGenes > 0;
die "-matrixStorage must be auto, dense, or sparse\n" unless $matrixStorage =~ /^(auto|dense|sparse)$/;


my $speciesGTDB = $useGTDBmg eq "GTDB" ? "GTDB_GTDB" : "specI_GTDB";
my $MGterm = $useGTDBmg eq "GTDB" ? "GTDBmg" : "FMG";
die "-MGset option has to be \"GTDB\" or \"FMG\"\n" unless ($useGTDBmg eq "GTDB" || $useGTDBmg eq "FMG");
my $MGdir = "$GCd/$MGterm/";
if ($outD eq "") {
	$outD = $MGSfile ne "" ? "$GCd/Anno/Tax/${MGterm}_MGS/" : "$GCd/Anno/Tax/$MGterm/";
}
my $GTDBspecI = getProgPaths($speciesGTDB);


print "-------------------------------------------------------------------------\nannotateMG script v $version\n-------------------------------------------------------------------------\n";
print "Using gene cat in $GCd\nMGS: $MGSfile\nMGStax: $MGStax\ncores: $BlastCores MGset: $useGTDBmg minGenes: $minGenes matrixStorage: $matrixStorage\n";


make_path($outD) unless -d $outD;
die "Can't create output directory $outD\n" unless -d $outD;

# Taxonomy is retained in Perl because the worker intentionally receives only
# explicit assignment records.
my %specIfullTax = %{readTable($GTDBspecI,"\t")};
make_path($MGdir) unless -d $MGdir;
die "Can't create marker-gene output directory $MGdir\n" unless -d $MGdir;
#assign each COG separately

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
print "Reading gene assignments (LCA files) from $MGdir\n";
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
if ($allOK==0){die "One or more Blast files were not ok.. restart procedure\n";}

# Preserve v2's MGS/specI transfer semantics before producing explicit worker
# inputs.  The C++ worker deliberately has no knowledge of MATAFILER paths or
# taxonomy databases.
my $SI2MGS = MGSassign();
transferSI2MGS($SI2MGS);
undef %speci2MGS;

my $handoffPrefix = "$outD/speci_assign";
my ($seedFile,$candidateFile,$geneCogFile,$seedCount,$candidateCount) = writeSpeciAssignHandoffs(
	"$handoffPrefix.seed_members.tsv", "$handoffPrefix.candidate_edges.tsv", "$handoffPrefix.gene_cogs.tsv");
die "No high-confidence specI/MGS seed members were available for cc.bin --speci-assign\n"
	unless $seedCount;
die "No specI/MGS candidate edges were available for cc.bin --speci-assign\n"
	unless $candidateCount;

# The C++ worker always uses compact CSR storage.  Keep accepting the v2
# option so existing callers retain the same CLI, but do not retain a second
# Perl matrix representation.
warn "-matrixStorage=$matrixStorage is ignored by annotateMGwSpecIs3.pl; cc.bin --speci-assign always uses CSR storage.\n"
	if $matrixStorage ne "auto";
undef $SI2MGS;
undef %SpecIgenes;
undef %Q2S;
undef %Gene2MGS;
undef %gene2COG;

my $matrixFile = "$GCd/Matrix.$MGterm.mat";
die "Marker-gene matrix not found: $matrixFile\n" unless -e $matrixFile;
my $speciAssignBin = getProgPaths("canopy");
runSpeciAssign($speciAssignBin,$matrixFile,$seedFile,$candidateFile,$geneCogFile,$handoffPrefix);

my $assignmentsFile = "$handoffPrefix.assignments.tsv";
my $workerAudit = "$handoffPrefix.assignment_resolution.tsv";
my $workerMatrix = "$handoffPrefix.speci.mat";
my $assignmentCount = importSpeciAssignments($assignmentsFile,"$MGdir/gene2specI.txt");
rename($workerAudit,"$outD/assignment_resolution.tsv")
	or die "Can't publish assignment audit $workerAudit: $!\n";
rename($workerMatrix,"$outD/specI.mat")
	or die "Can't publish specI matrix $workerMatrix: $!\n";
writeTaxonomyMatrices("$outD/specI.mat",\%specIfullTax);

print "Finished SpecI annotations & matrix: $outD ($assignmentCount canonical assignments)\n";
exit(0);
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
		# The selected specI recorded in MGS2speci.txt is authoritative over
		# a conflicting input MGS taxonomy.
		$specIfullTax{$MGS} = [ @{$specIfullTax{$valSI}} ]
			if exists($specIfullTax{$valSI});
		foreach my $COG (keys %{$SpecIgenes{$valSI}}){
			foreach my $gen (@{$SpecIgenes{$valSI}{$COG}}){
				$Q2S{$gen}{$MGS} = $Q2S{$gen}{$valSI};
				delete $Q2S{$gen}{$valSI};
			}
			my @adds = @{$SpecIgenes{$valSI}{$COG}};
			if (exists($SpecIgenes{$MGS}{$COG})){ 
				push(@adds,@{$SpecIgenes{$MGS}{$COG}});
				@{$SpecIgenes{$MGS}{$COG}} = do { my %seen; grep { !$seen{$_}++ } @adds };
			} else {
				@{$SpecIgenes{$MGS}{$COG}} = @adds;
			}
			delete $SpecIgenes{$valSI}{$COG};
		}
		delete $SpecIgenes{$valSI};
	}
}

sub writeSpeciAssignHandoffs{
	my ($seedFile,$candidateFile,$geneCogFile) = @_;
	open my $seedFH, ">", $seedFile or die "Can't write seed handoff $seedFile: $!\n";
	open my $candidateFH, ">", $candidateFile or die "Can't write candidate handoff $candidateFile: $!\n";
	open my $geneCogFH, ">", $geneCogFile or die "Can't write gene-to-COG handoff $geneCogFile: $!\n";
	print {$seedFH} "specI\tgene\tCOG\tsource\n";
	print {$candidateFH} "gene\tspecI\tCOG\tsource\n";
	print {$geneCogFH} "gene\tCOG\n";

	my $geneCogCount = 0;
	foreach my $gene (sort keys %gene2COG){
		my $cog = $gene2COG{$gene};
		next unless defined($cog) && length($cog);
		print {$geneCogFH} "$gene\t$cog\n";
		$geneCogCount++;
	}
	close $geneCogFH or die "Can't close gene-to-COG handoff $geneCogFile: $!\n";

	my (%candidateSeen,%seedSeen);
	my ($seedCount,$candidateCount) = (0,0);
	foreach my $owner (sort keys %SpecIgenes){
		foreach my $cog (sort keys %{$SpecIgenes{$owner}}){
			my %seenGenes;
			my @genes = sort grep { length($_) && !$seenGenes{$_}++ } @{$SpecIgenes{$owner}{$cog}};
			next unless @genes;
			foreach my $gene (@genes){
				my $key = join("\x1e",$gene,$owner,$cog);
				next if $candidateSeen{$key}++;
				my $source = exists($Gene2MGS{$gene}) && exists($MGSlist{$owner}) && $Gene2MGS{$gene} eq $owner ? "MGS" : "specI";
				print {$candidateFH} "$gene\t$owner\t$cog\t$source\n";
				$candidateCount++;
			}

			# Match v2's high-confidence seed rule: exactly one copy in the
			# owner/COG slot and no competing LCA owner for the marker gene.
			next unless @genes == 1;
			my $gene = $genes[0];
			my $ownerCount = exists($Q2S{$gene}) ? scalar(keys %{$Q2S{$gene}}) : 0;
			next if $ownerCount > 1;
			my $seedKey = join("\x1e",$owner,$gene,$cog);
			next if $seedSeen{$seedKey}++;
			my $source = exists($Gene2MGS{$gene}) && exists($MGSlist{$owner}) && $Gene2MGS{$gene} eq $owner ? "MGS" : "specI";
			print {$seedFH} "$owner\t$gene\t$cog\t$source\n";
			$seedCount++;
		}
	}
	close $seedFH or die "Can't close seed handoff $seedFile: $!\n";
	close $candidateFH or die "Can't close candidate handoff $candidateFile: $!\n";
	print "Wrote $seedCount seed members, $candidateCount constrained candidate edges, and $geneCogCount gene-to-COG background mappings\n";
	return ($seedFile,$candidateFile,$geneCogFile,$seedCount,$candidateCount);
}

sub shellQuote{
	my ($value) = @_;
	$value =~ s/'/'"'"'/g;
	return "'$value'";
}

sub runSpeciAssign{
	my ($binary,$matrixFile,$seedFile,$candidateFile,$geneCogFile,$outputPrefix) = @_;
	my @arguments = (
		"--speci-assign",
		"--matrix", $matrixFile,
		"--seed-members", $seedFile,
		"--candidate-edges", $candidateFile,
		"--gene-cogs", $geneCogFile,
		"--output-prefix", $outputPrefix,
		"--num-threads", $BlastCores,
		"--min-genes", $minGenes,
		"--min-guide-nonzero", $minGuideNonzero,
		"--correlation-threshold", $globalCorrThreshold,
		"--near-best-window", $nearBestWindow,
	);

	my $status;
	if ($binary =~ /\n/) {
		# getProgPaths may prepend a configured environment activation snippet.
		# Preserve that site configuration while shell-quoting every data argument.
		my @lines = grep { /\S/ } split /\n/, $binary;
		my $launch = $lines[-1] // "";
		$launch =~ s/^\s+|\s+$//g;
		my ($executable) = split /\s+/, $launch, 2;
		die "Configured canopy executable does not exist: $launch\n" unless length($executable) && -e $executable;
		die "Configured canopy executable is not executable: $executable\n" unless -x $executable;
		my $command = $binary." ".join(" ", map { shellQuote($_) } @arguments);
		print "Running constrained sparse specI assignment: $command\n";
		$status = system("bash", "-c", $command);
	} else {
		die "Configured canopy executable does not exist: $binary\n" unless -e $binary;
		die "Configured canopy executable is not executable: $binary\n" unless -x $binary;
		my @command = ($binary,@arguments);
		print "Running constrained sparse specI assignment: ".join(" ",@command)."\n";
		$status = system(@command);
	}
	die "Couldn't start cc.bin --speci-assign: $!\n" if $status == -1;
	die "cc.bin --speci-assign failed with exit status ".($status >> 8)."\n" if $status != 0;
	foreach my $output ("$outputPrefix.assignments.tsv", "$outputPrefix.assignment_resolution.tsv", "$outputPrefix.speci.mat"){
		die "cc.bin --speci-assign did not create expected output $output\n" unless -e $output;
	}
}

sub importSpeciAssignments{
	my ($inputFile,$legacyFile) = @_;
	open my $input, "<", $inputFile or die "Can't read speci_assign assignments $inputFile: $!\n";
	my $header = <$input>;
	die "speci_assign assignments are empty: $inputFile\n" unless defined($header);
	chomp $header;
	my @header = split /\t/, $header, -1;
	my %column;
	@column{@header} = (0..$#header);
	foreach my $required (qw(gene specI COG)){
		die "speci_assign assignments lack required $required column: $inputFile\n"
			unless exists($column{$required});
	}
	open my $legacy, ">", $legacyFile or die "Can't write legacy gene2specI output $legacyFile: $!\n";
	my $count = 0;
	my %seen;
	while (my $line = <$input>){
		chomp $line;
		next if $line eq "";
		my @fields = split /\t/, $line, -1;
		my ($gene,$owner,$cog) = @fields[$column{gene},$column{specI},$column{COG}];
		die "Malformed speci_assign assignment row in $inputFile: $line\n"
			unless defined($gene) && defined($owner) && defined($cog) && length($gene) && length($owner) && length($cog);
		die "Duplicate final assignment for marker gene $gene in $inputFile\n" if $seen{$gene}++;
		print {$legacy} "$gene\t$owner\t$cog\n";
		$count++;
	}
	close $input or die "Can't close speci_assign assignments $inputFile: $!\n";
	close $legacy or die "Can't close legacy gene2specI output $legacyFile: $!\n";
	return $count;
}

sub writeTaxonomyMatrices{
	my ($matrixFile,$taxonomy) = @_;
	open my $input, "<", $matrixFile or die "Can't read specI matrix $matrixFile: $!\n";
	my $header = <$input>;
	die "specI matrix is empty: $matrixFile\n" unless defined($header);
	chomp $header;
	my @header = split /\t/, $header, -1;
	die "Malformed specI matrix header in $matrixFile\n" unless @header >= 3 && $header[0] eq "SpecI";
	my @samples = @header[1..$#header];
	my $background = <$input>;
	die "specI matrix lacks the '?' background row: $matrixFile\n" unless defined($background);
	chomp $background;
	my @background = split /\t/, $background, -1;
	die "Malformed specI background row in $matrixFile\n"
		unless @background == @header && $background[0] eq "?";
	my @backgroundValues = @background[1..$#background];

	my @taxLevels = qw(superkingdom phylum class order family genus species);
	my @levelMaps = map { {} } @taxLevels;
	my $taxFile = $matrixFile;
	$taxFile =~ s/\.[^\.]*$/.tax/;
	open my $taxOut, ">", $taxFile or die "Can't write specI taxonomy $taxFile: $!\n";
	while (my $line = <$input>){
		chomp $line;
		next if $line eq "";
		my @fields = split /\t/, $line, -1;
		die "Malformed specI matrix row in $matrixFile: $line\n" unless @fields == @header;
		my $owner = shift @fields;
		my $tax = $taxonomy->{$owner};
		die "No taxonomy available for final specI/MGS owner $owner\n" unless defined($tax);
		die "Incomplete taxonomy for final specI/MGS owner $owner\n" unless @{$tax} >= @taxLevels;
		print {$taxOut} "$owner\t".join("\t",@{$tax})."\n";
		for (my $level = 0; $level < @taxLevels; $level++){
			my $label = join(";",@{$tax}[0..$level]);
			my $target = $levelMaps[$level]{$label};
			if (!defined($target)){
				$levelMaps[$level]{$label} = [map { $_ + 0 } @fields];
				next;
			}
			for (my $column = 0; $column < @fields; $column++){
				$target->[$column] += $fields[$column];
			}
		}
	}
	close $input or die "Can't close specI matrix $matrixFile: $!\n";
	close $taxOut or die "Can't close specI taxonomy $taxFile: $!\n";

	my $outputPrefix = $matrixFile;
	$outputPrefix =~ s/\.[^\.]*$//;
	for (my $level = 0; $level < @taxLevels; $level++){
		my $output = "$outputPrefix.$taxLevels[$level]";
		open my $out, ">", $output or die "Can't write taxonomy matrix $output: $!\n";
		print {$out} "$taxLevels[$level]\t".join("\t",@samples)."\n";
		print {$out} "?\t".join("\t",@backgroundValues)."\n";
		foreach my $label (sort keys %{$levelMaps[$level]}){
			print {$out} "$label\t".join("\t",@{$levelMaps[$level]{$label}})."\n";
		}
		close $out or die "Can't close taxonomy matrix $output: $!\n";
	}
}
