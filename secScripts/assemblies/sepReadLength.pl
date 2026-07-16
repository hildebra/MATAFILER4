#!/usr/bin/env perl
use strict;
use warnings;

die "Usage: $0 <minimum-long-length> <fasta>\n" unless @ARGV == 2;
my ($minimum_length, $input) = @ARGV;
die "Length threshold must be a non-negative integer\n"
	unless $minimum_length =~ /\A\d+\z/;
die "Input FASTA is missing or empty: $input\n" unless -s $input;

my $long_path = "$input.long";
my $short_path = "$input.short.$$";
open my $in, '<', $input or die "Cannot open $input: $!\n";
open my $long, '>', $long_path or die "Cannot open $long_path: $!\n";
open my $short, '>', $short_path or die "Cannot open $short_path: $!\n";

my ($header, $sequence) = ('', '');
my ($long_count, $short_count) = (0, 0);
sub emit_record {
	return unless length $header;
	my $target;
	if (length($sequence) > $minimum_length) {
		$target = $long;
		$long_count++;
	} else {
		$target = $short;
		$short_count++;
	}
	print {$target} "$header\n$sequence\n" or die "Cannot write split FASTA: $!\n";
}

while (my $line = <$in>) {
	$line =~ s/[\r\n]+\z//;
	if ($line =~ /^>/) {
		emit_record();
		($header, $sequence) = ($line, '');
	} else {
		die "Sequence encountered before first FASTA header in $input\n"
			unless length $header;
		$line =~ s/\s+//g;
		$sequence .= $line;
	}
}
emit_record();
die "No FASTA records found in $input\n" unless $long_count + $short_count;
close $in or die "Cannot close $input: $!\n";
close $long or die "Cannot close $long_path: $!\n";
close $short or die "Cannot close $short_path: $!\n";
rename $short_path, $input or die "Cannot replace $input with $short_path: $!\n";
print "Long records: $long_count; retained short records: $short_count\n";
