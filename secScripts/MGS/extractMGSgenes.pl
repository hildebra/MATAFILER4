#!/usr/bin/perl
#extracts FAA/FNA for each MGS
#perl extractMGSgenes.pl /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/Canopy4_AC/clusters.txt /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/Canopy4_AC//extr/ /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/ 200
#perl extractMGSgenes.pl /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/Binning/MetaBat/MB2.clusters.ext.can /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/Binning/MetaBat/extr/ /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/ 200 /scratch/bork/hildebra/MGStest/



use strict;
use warnings;
use Mods::GenoMetaAss qw(systemW readFasta);
use Mods::Binning qw(runCheckM);
use File::Path qw(make_path remove_tree);



sub readCluster{
	my ($cF) = @_;
	my %canos;
	open I,"<$cF" or die "can;t open $cF\n";
	while (my $line = <I>){
		chomp $line;
		$line =~ s/ //g;
		my @spl = split /\t/,$line;
		push(@{$canos{$spl[0]}}, $spl[1]);
	}
	close I;
	return \%canos;
}
my $ncore = 20;


die "Usage: $0 CLUSTERS OUT-DIR GC-DIR MIN-GENES [TMP-DIR] [CLUSTER-ID]\n"
	unless @ARGV >= 4 && @ARGV <= 6;
my $cluF = $ARGV[0];
my $oDir = $ARGV[1];
my $GCd = $ARGV[2];
my $tmpD = "$oDir/tmp/";
my $minGenes = $ARGV[3];
$tmpD = $ARGV[4] if @ARGV >= 5;
my $clusterID = @ARGV >= 6 ? $ARGV[5] : 95;
die "MIN-GENES must be a non-negative integer\n" unless $minGenes =~ /^\d+$/;
die "CLUSTER-ID must be between 1 and 100\n"
	unless $clusterID =~ /^\d+$/ && $clusterID >= 1 && $clusterID <= 100;
remove_tree($oDir) if -d $oDir;
make_path($oDir, $tmpD);

my $hr = readCluster($cluF);
my %clust = %{$hr};
my %wantedCatalogueGenes;
for my $cluster (keys %clust) {
	next if @{$clust{$cluster}} < $minGenes;
	$wantedCatalogueGenes{$_} = 1 for @{$clust{$cluster}};
}

print "Reading ref FNA..\n";
$hr = readFasta("$GCd/compl.incompl.$clusterID.fna", 1, "\\s", \%wantedCatalogueGenes, { fai => 1 });
my %FNA = %{$hr};
#my @test = keys %FNA; print "$test[0] $test[1] $test[123]\n"; print "$FNA{13220655}\n";
foreach my $cl (sort keys %clust){
	my $oF = "$oDir/$cl.fna";
	my @refG = @{$clust{$cl}};
	next if (scalar(@refG) < $minGenes);
	open O,">$oF" or die $!;
	foreach my $rg (@refG){
		chomp $rg;
		die "Can't find gene $rg in gene cat\n" unless(exists($FNA{$rg}));
		my $rn = $rg; $rn =~ s/://;
		print O ">$rn\n$FNA{$rg}\n";
	}
	close O;
}


print "Reading ref FAA..\n";
$hr = readFasta("$GCd/compl.incompl.$clusterID.prot.faa", 1, "\\s", \%wantedCatalogueGenes, { fai => 1 });
my %FAA = %{$hr};
foreach my $cl (sort keys %clust){
	my $oF = "$oDir/$cl.faa";
	my @refG = @{$clust{$cl}};
	next if (scalar(@refG) < $minGenes);
	open O,">$oF" or die $!;
	foreach my $rg (@refG){
		die "Can't find gene $rg in gene cat\n" unless(exists($FAA{$rg}));
		my $rn = $rg; $rn =~ s/://;
		print O ">$rn\n$FAA{$rg}\n";
	}
	close O;
}

print "Starting checkm...\n";
my $outFile = $cluF.".cm";
runCheckM($oDir,$outFile,$tmpD,$ncore) unless (-e $outFile);



print "Done\n";


