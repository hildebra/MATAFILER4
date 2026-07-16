#!/usr/bin/env perl
use strict;
use warnings;

die "Usage: $0 <fasta> <maximum-fragment-length>\n" unless @ARGV == 2;
my ($input, $fragment_length) = @ARGV;
die "Fragment length must be a positive integer\n"
	unless $fragment_length =~ /\A\d+\z/ && $fragment_length > 0;
die "Input FASTA is missing or empty: $input\n" unless -s $input;

my $temporary = "$input.tmp.$$";
open my $in, '<', $input or die "Cannot open $input: $!\n";
open my $out, '>', $temporary or die "Cannot open $temporary: $!\n";

my ($header, $sequence) = ('', '');
my ($records, $fragments) = (0, 0);

sub emit_record {
	return unless length $header;
	die "FASTA record '$header' has an empty sequence\n" unless length $sequence;
	$records++;
	my $part = 0;
	for (my $start = 0; $start < length($sequence); $start += $fragment_length) {
		$part++;
		my $fragment = substr($sequence, $start, $fragment_length);
		my $end = $start + length($fragment);
		my $fragment_header = $header;
		if (length($sequence) > $fragment_length) {
			my $suffix = sprintf('_part%d_%d-%d', $part, $start + 1, $end);
			$fragment_header =~ s/^>(\S+)/>$1$suffix/;
		}
		print {$out} "$fragment_header\n$fragment\n"
			or die "Cannot write $temporary: $!\n";
		$fragments++;
	}
}

while (my $line = <$in>) {
	$line =~ s/[\r\n]+\z//;
	if ($line =~ /^>/) {
		emit_record();
		$header = $line;
		$sequence = '';
		next;
	}
	die "Input appears to be FASTQ, not FASTA\n" if $line =~ /^\@/ && !length($header);
	die "Sequence encountered before first FASTA header in $input\n"
		unless length $header;
	$line =~ s/\s+//g;
	die "Invalid empty sequence line in $input\n" unless length $line;
	$sequence .= $line;
}
emit_record();
die "No FASTA records found in $input\n" unless $records;

close $in or die "Cannot close $input: $!\n";
close $out or die "Cannot close $temporary: $!\n";
rename $temporary, $input or die "Cannot replace $input with $temporary: $!\n";
print "Split $records sequence(s) into $fragments fragment(s) of at most $fragment_length bp\n";
