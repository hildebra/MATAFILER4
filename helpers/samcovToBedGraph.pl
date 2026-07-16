#!/usr/bin/env perl
use strict;
use warnings;

# Convert samtools depth output to zero-based, half-open bedGraph intervals.
# A coordinate gap must start a new interval even if the depth is unchanged.

my ($chromosome, $start, $end, $depth);
while (my $line = <STDIN>) {
	$line =~ s/[\r\n]+$//;
	next if ($line eq '');
	my ($next_chromosome, $position, $next_depth) = split /\s+/, $line;
	die "Malformed samtools depth line: $line\n"
		unless (defined($next_depth) && $position =~ /^\d+$/ && $position > 0);
	my $zero_based = $position - 1;
	if (!defined $chromosome) {
		($chromosome, $start, $end, $depth) =
			($next_chromosome, $zero_based, $position, $next_depth);
		next;
	}
	if ($next_chromosome eq $chromosome
		&& $zero_based == $end && $next_depth == $depth) {
		$end = $position;
		next;
	}
	print join("\t", $chromosome, $start, $end, $depth), "\n";
	($chromosome, $start, $end, $depth) =
		($next_chromosome, $zero_based, $position, $next_depth);
}
print join("\t", $chromosome, $start, $end, $depth), "\n"
	if (defined $chromosome);
