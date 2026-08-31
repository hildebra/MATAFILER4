#!/usr/bin/env perl
#takes as input motu cluster name and extracts the MG from Shini's DB

use warnings;
use strict;
use Mods::IO_Tamoc_progs qw(getProgPaths);
#use Mods::Subm qw(qsubSystem emptyQsubOpt);
#Cluster1088 = TEC6, test case

my $test="Cluster1088";
my $tarMotu = $test;
my $oDir = getProgPaths("legacyMotuMGOutput");
system "mkdir -p $oDir" unless (-d $oDir);
my $outF = "$oDir/$tarMotu.MG.fna";
unlink $outF if (-e $outF);
my $linkF = getProgPaths("legacyMotuLinkMap");
my $mapF = getProgPaths("legacyMotuMap");
my $fastaRef = getProgPaths("legacyMotuFasta");



my $gns = `grep '$tarMotu' $linkF  | cut -f1`;
chomp $gns;
my @genes1 = split /\n/,$gns;

for (my $i=0;$i<@genes1;$i++){
	$genes1[$i] =~ m/(COG\d+)\.(.*)/;
	print "$genes1[$i]\n";
	my $geneTarX = `grep -P '$1\\t$2\$' $mapF  | cut -f1`;
	my @tmp = split /\n/,$geneTarX;
	my $geneTar =  $tmp[0];
	
	print "$geneTar\n\n";
	my $smCmd = "samtools faidx  $fastaRef '$geneTar' >> $outF";
	if (system($smCmd)){die "Failed $smCmd\n";}
}
print "Finished grabbing ". @genes1 . " MGs\n$outF\n";
exit 0;