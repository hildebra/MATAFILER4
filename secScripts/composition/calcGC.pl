#!/usr/bin/perl
#./calcGC.pl /g/bork1/hildebra/SNP/GNMassSimu2/AssmblGrp_1/metag/scaffolds.fasta.filt  /g/bork1/hildebra/SNP/GNMassSimu2/AssmblGrp_1/metag/ContigStats/sgc.2
use warnings;
use strict;

use Mods::GenoMetaAss qw(gzipopen);


my $inF = $ARGV[0];
my $outF = $ARGV[1];
my $isGenes = 0;
$isGenes = 1 if (@ARGV > 2);


#open I,"<$inF" or die "Can't open $inF";
my ($FAS ,$status) = gzipopen($inF,"calcGC infile",1);

open O,">$outF" or die "Can't open $outF";
print O "contig\tGC\n";

if ($isGenes){
	open O3,">${outF}3" or die "Can't open ${outF}3";
	print O3 "contig\tGC\n";
	my $ctgF = "${outF}3";
	$ctgF =~ s/\.pergene//;
	open OC3,">${ctgF}" or die "Can't open ${ctgF}";
	print OC3 "contig\tGC\n";
}

my ($curTag, $curCtg, $sequence) = ("", "", "");
my ($contigGC3, $contigAT3) = (0, 0);
while (my $line = <$FAS>) {
	if ($line =~ /^>(\S+)/) {
		my $tag = $1;
		(my $ctg = $tag) =~ s/_\d+$//;
		emit_sequence($ctg ne $curCtg) if $curTag ne "";
		($curTag, $curCtg, $sequence) = ($tag, $ctg, "");
	} else {
		chomp $line;
		$sequence .= uc($line);
	}
}
emit_sequence(1) if $curTag ne "";
close O; close $FAS;
if ($isGenes) { close O3; close OC3; }

sub gc_row {
	my ($handle, $tag, $gc, $at) = @_;
	my $value = $gc + $at ? sprintf('%.3f', 100 * $gc / ($gc + $at)) : -1;
	print {$handle} "$tag\t$value\n";
}

sub emit_sequence {
	my ($endContig) = @_;
	gc_row(\*O, $curTag, scalar($sequence =~ tr/GC//), scalar($sequence =~ tr/AT//));
	return unless $isGenes;
	die "Gene length is not a multiple of three: $curTag\n" if length($sequence) % 3;
	my $third = "";
	for (my $i = 2; $i < length($sequence); $i += 3) { $third .= substr($sequence, $i, 1); }
	my $gc = $third =~ tr/GC//;
	my $at = $third =~ tr/AT//;
	gc_row(\*O3, $curTag, $gc, $at);
	$contigGC3 += $gc; $contigAT3 += $at;
	if ($endContig) {
		gc_row(\*OC3, $curCtg, $contigGC3, $contigAT3);
		($contigGC3, $contigAT3) = (0, 0);
	}
}
