#!/usr/bin/env perl
use warnings;
use strict;
use File::Path qw(make_path);
use Mods::GenoMetaAss qw(systemW);
use Mods::IO_Tamoc_progs qw(getProgPaths);

die "Usage: $0 <contigs.fasta> <temporary-dir> <genes.gff> <threads>\n" unless @ARGV == 4;
my ($ctgFile, $tmpPath, $gffFile, $threads) = @ARGV;
die "Contig FASTA is missing or empty: $ctgFile\n" unless -s $ctgFile;
die "GFF is missing or empty: $gffFile\n" unless -s $gffFile;
die "Thread count must be a positive integer\n" unless $threads =~ /^\d+$/ && $threads > 0;

my $whok = getProgPaths('whokaryote');
my $outdir = "$tmpPath/whoK";
make_path($outdir) unless -d $outdir;
systemW("$whok --contigs $ctgFile --outdir $outdir --gff $gffFile --minsize 4000 --threads $threads\n");
die "Whokaryote did not create its output directory: $outdir\n" unless -d $outdir;
