#!/usr/bin/env perl
# Best hit per sequence from a HMMER --domtblout table (replaces bmn-HMMerBestHit.py, which needed python2
# and kept the wrong hit for every gene but the first).
# Per target sequence (column 1) the line with the highest score in the score column is printed
# (default 14 = domain score; on ties the later line wins); hits below the minimum score are ignored.
# The input does not need to be sorted; output is sorted by sequence name.
use strict;
use warnings;
use Getopt::Long;

my $usage = "Usage: $0 [-c score_column (default 14)] [-s min_score (default 25)] <domtblout>\n";
my ($col, $minScore) = (14, 25);
GetOptions("c|column=i" => \$col, "s|minscore=f" => \$minScore) or die $usage;
die $usage unless @ARGV == 1;
die "-c must be at least 1\n" if $col < 1;
my $in = $ARGV[0];

open my $fh, '<', $in or die "HMM output file not found or unreadable: $in\n";
my (%best, %bestScore);
while (my $line = <$fh>) {
	next if $line =~ /^#/ || $line !~ /\S/;
	$line .= "\n" unless $line =~ /\n\z/;
	my @f = split ' ', $line;
	die "Line $. of $in has " . scalar(@f) . " columns, score column $col is missing\n" if @f < $col;
	my $score = $f[$col - 1];
	die "Non-numeric score '$score' in line $. of $in\n"
		unless $score =~ /^[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?$/;
	next if $score < $minScore;
	my $seq = $f[0];
	if (!exists $bestScore{$seq} || $score >= $bestScore{$seq}) {
		$best{$seq} = $line;
		$bestScore{$seq} = $score;
	}
}
close $fh or die "Cannot close $in: $!\n";
print $best{$_} for sort keys %best;
