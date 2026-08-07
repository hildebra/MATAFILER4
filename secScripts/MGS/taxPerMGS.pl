#!/usr/bin/perl
#uses MGS to cluster genes and kraken tax assignments
#./taxPerMGS.pl /g/scb/bork/hildebra/SNP/GCs/DramaGCv5//Binning/MetaBat//MB2.clusters.ext.can.Rhcl /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/
use warnings; use strict;

use Mods::Binning qw(readMGSrev );


die "Usage: $0 MGS-file gene-catalog-dir output-prefix\n" unless @ARGV == 3;
my $refMGf = $ARGV[0];
my $GCd = $ARGV[1];
my $outF = $ARGV[2];



my $LCAout="$outF.LCA";
my $taxout = "$outF.tax";
my $krakF = "$GCd/Anno/Tax/krak2.txt"; #krak_0.01.txt

if (-s $LCAout && -s $taxout){exit(0);}
unless (-s $krakF) {
	warn "Optional Kraken input is missing or empty; skipping MGS Kraken taxonomy:\n$krakF\n";
	exit(0);
}

my $hr = readMGSrev($refMGf);
my %MGs = %{$hr};
my %allMGS = map { $_ => 1 } values %MGs;
my %tCnt;
#my $krakF = "$GCd/Anno/Tax/krak_0.01.txt"; #
open my $krakenIn, '<', $krakF or do {
	warn "Optional Kraken input could not be read; skipping MGS Kraken taxonomy:\n$krakF: $!\n";
	exit(0);
};
while (<$krakenIn>){
	chomp;my @spl=split /\t/;
	next unless (exists($MGs{$spl[0]}));
	#my @s2 = split /;/,$spl[1];
	unless (defined $spl[1]){
		next;
	}
	my $i=0;
	foreach my $t (split /;/,$spl[1]){
		$tCnt{$MGs{$spl[0]}}{$i}{$t}++;
		$i++;
	}
}
close $krakenIn;
my %tStat;
open OL,">$LCAout" or die "Cannot write $LCAout: $!\n";
open OC,">$taxout" or die "Cannot write $taxout: $!\n";
foreach my $mgs (sort keys %allMGS){
	my $tax=""; my $maxD=0;
	my $lineage_broken = 0;
	print OC "$mgs\t";
	print OL "$mgs\t";
	for (my $i=0;$i<8;$i++){
		if ($lineage_broken || !exists($tCnt{$mgs}{$i})) {
			$tax .= "?;";
			$lineage_broken = 1;
			next;
		}
		my %curT = %{$tCnt{$mgs}{$i}}; my $tSum=0;
		print OC "\t$i";
		foreach my $t (sort keys %curT){
			$tSum+=$curT{$t};
			print OC "$t:$curT{$t};";
		}
		my ($maxT) = sort { $curT{$b} <=> $curT{$a} || $a cmp $b } keys %curT;
		my $max = defined($maxT) ? $curT{$maxT} : 0;
		if ($tSum < 100) {
			$tax .= "?;";
			$lineage_broken = 1;
			next;
		}
		if (($max/$tSum) > 0.8){
			$tax .= "${maxT};";
			$maxD=$i;
		} else {
			$tax .= "?;";
			$lineage_broken = 1;
		}
	}
	$tStat{$maxD}++;
	print OC "\n";
	print OL "$tax\n";
}
close OL or die "Cannot close $LCAout: $!\n";
close OC or die "Cannot close $taxout: $!\n";
#print stats how many LCA levels were hit by MGS:
foreach my $d (sort(keys%tStat)){
	print "$d:$tStat{$d}\t";
}
print "\n";
