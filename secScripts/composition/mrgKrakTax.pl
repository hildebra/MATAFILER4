#!/usr/bin/env perl
# Merge per-sample Kraken read counts ("<lineage><TAB><reads>", from krak2_count_tax.pl) into one
# matrix: lineages in rows, samples in columns, 0 where a sample has no reads of a lineage.
#   mrgKrakTax.pl <suffix> <output matrix> <sample><suffix> ...
# The sample name is the file name without directory and <suffix> (e.g. ".0.01.krak.txt").
use strict;
use warnings;
use File::Basename qw(basename);

die "Usage: $0 <file suffix> <output matrix> <count files...>\n" unless @ARGV >= 3;
my ($suffix, $outF, @files) = @ARGV;
my (%counts, @samples);
for my $f (@files) {
	my $sample = basename($f);
	$sample =~ s/\Q$suffix\E$// or die "$f does not end in $suffix\n";
	push @samples, $sample;
	open my $in, '<', $f or die "Cannot open $f: $!\n";
	while (my $line = <$in>) {
		$line =~ s/[\r\n]+$//;
		next if $line eq '';
		my ($lineage, $n) = split /\t/, $line;
		die "Malformed count line in $f: $line\n" unless defined($n) && $n =~ /^\d+$/;
		$counts{$lineage}{$sample} += $n;
	}
	close $in;
}
open my $out, '>', "$outF.tmp" or die "Cannot write $outF.tmp: $!\n";
print {$out} join("\t", 'Taxon', @samples), "\n";
for my $lineage (sort keys %counts) {
	print {$out} join("\t", $lineage, map { $counts{$lineage}{$_} // 0 } @samples), "\n";
}
close $out or die "Cannot write $outF.tmp: $!\n";
rename("$outF.tmp", $outF) or die "Cannot rename $outF.tmp: $!\n";
