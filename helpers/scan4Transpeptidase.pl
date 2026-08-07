#!/usr/bin/env perl

use warnings; use strict;
use Mods::IO_Tamoc_progs qw(getProgPaths);

#http://pfam.xfam.org/family/Transpeptidase#tabview=tab6
my $hmmBin3 = getProgPaths("hmmsearch");

my $hmmModel = getProgPaths("transpeptidaseHMM_DB");
my $faa = getProgPaths("legacyTranspeptidaseFAA");
print "$hmmBin3 -Z 11927849 -E 1000 --cpu 4 $hmmModel $faa >".getProgPaths("legacyTranspeptidaseOutput")."\n";

print "samtools faidx $faa 'fig|6666666.162524.peg.1641' 'fig|6666666.162524.peg.2390' >".getProgPaths("legacyTranspeptidaseOutput").".faa\n";