#!/usr/bin/env perl
use strict;
use warnings;
use Mods::IO_Tamoc_progs qw(getProgPaths);

die "Usage: $0 <bam> <reference-name-marker> <rewrite-body:0|1>\n" unless @ARGV == 3;
my ($input_bam, $reference_marker, $rewrite_body) = @ARGV;
die "Input BAM is missing or empty: $input_bam\n" unless -s $input_bam;
die "Reference marker must not be empty\n" unless length $reference_marker;
die "rewrite-body must be 0 or 1\n" unless $rewrite_body =~ /\A[01]\z/;
my $samtools = getProgPaths('samtools');

sub keep_header_line {
	my ($line) = @_;
	return 1 unless $line =~ /^\@SQ\t/;
	my ($name) = $line =~ /(?:^|\t)SN:([^\t\r\n]+)/;
	return defined($name) && index($name, $reference_marker) >= 0;
}

open my $header, '-|', $samtools, 'view', '-H', $input_bam
	or die "Cannot read BAM header from $input_bam: $!\n";
my @filtered_header = grep { keep_header_line($_) } <$header>;
close $header or die "samtools failed while reading the header from $input_bam\n";

if (!$rewrite_body) {
	print @filtered_header or die "Cannot write filtered BAM header: $!\n";
	exit 0;
}

my $temporary_bam = "$input_bam.reheader.$$.bam";
open my $converter, '|-', $samtools, 'view', '-b', '-o', $temporary_bam, '-'
	or die "Cannot start samtools BAM writer: $!\n";
print {$converter} @filtered_header or die "Cannot write filtered header to samtools: $!\n";
open my $body, '-|', $samtools, 'view', $input_bam
	or die "Cannot read BAM records from $input_bam: $!\n";
while (my $line = <$body>) {
	print {$converter} $line or die "Cannot write BAM record to samtools: $!\n";
}
close $body or die "samtools failed while reading records from $input_bam\n";
close $converter or die "samtools failed while writing $temporary_bam\n";
system($samtools, 'index', $temporary_bam) == 0
	or die "samtools could not index $temporary_bam\n";
rename $temporary_bam, $input_bam
	or die "Cannot replace $input_bam with $temporary_bam: $!\n";
my $temporary_index = "$temporary_bam.bai";
rename $temporary_index, "$input_bam.bai"
	or die "Cannot replace $input_bam.bai with $temporary_index: $!\n";
open my $marker, '>', "$input_bam.hd" or die "Cannot create $input_bam.hd: $!\n";
close $marker or die "Cannot close $input_bam.hd: $!\n";
print "$input_bam\n";
