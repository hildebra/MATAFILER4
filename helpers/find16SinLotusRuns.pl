#!/usr/bin/perl
#takes predicted 16S and scans real 16S experiments for given 16S subset

use warnings;
use strict;
use Mods::IO_Tamoc_progs qw(getProgPaths);

my $mkBldbBin = getProgPaths("makeblastdb");
my $blastBin = getProgPaths("blastn");

my $lotusD = getProgPaths("legacyLotusRunsDir");
my $otuTar = $lotusD."otus.fa";
my $tar16sDB = getProgPaths("legacyLotus16S_DB");
my $taxblastf = getProgPaths("legacyLotusTaxBlast");
#my $TECdir = "/g/scb/bork/hildebra/SNP/GNMass2_singl/alien-11-374-0/Binning/";
#my $tar16sDB = $TECdir."TEC16.fa";
#my $taxblastf = $TECdir."HMPHits.blast";

my $cmd = "$mkBldbBin -in $tar16sDB -dbtype 'nucl'\n";
unless (-f $tar16sDB.".nhr"){	system($cmd);}
my $strand = "both";
#-perc_identity 75
$cmd = "$blastBin -query $otuTar -db $tar16sDB -out $taxblastf -perc_identity 95 -outfmt 6 -max_target_seqs 50 -evalue 0.001 -num_threads 60 -strand $strand \n"; #-strand plus both minus
system $cmd."\n";