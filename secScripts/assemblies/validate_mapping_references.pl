#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);

my ($assembly, $samtools, $sample_dirs) = ('', 'samtools', '');
GetOptions(
	'assembly=s'    => \$assembly,
	'samtools=s'    => \$samtools,
	'sample-dirs=s' => \$sample_dirs,
) or die "Invalid options\n";
die "Unexpected positional arguments: @ARGV\n" if @ARGV;
die "--assembly and --sample-dirs are required\n"
	unless length($assembly) && length($sample_dirs);
die "Assembly is missing or empty: $assembly\n" unless -s $assembly;

my $fai = "$assembly.fai";
if (!-s $fai || (stat($fai))[9] < (stat($assembly))[9]) {
	unlink $fai if -e $fai;
	system($samtools, 'faidx', $assembly) == 0
		or die "Could not build assembly index $fai with $samtools faidx\n";
}

my (%reference, @reference_order);
open my $fai_fh, '<', $fai or die "Could not read $fai: $!\n";
while (my $line = <$fai_fh>) {
	chomp $line;
	my ($name, $length) = split /\t/, $line, 3;
	die "Malformed FASTA index line in $fai: $line\n"
		unless defined($name) && length($name)
		&& defined($length) && $length =~ /^\d+$/;
	die "Duplicate reference name '$name' in $fai\n" if exists $reference{$name};
	$reference{$name} = 0 + $length;
	push @reference_order, $name;
}
close $fai_fh;
die "Assembly index contains no references: $fai\n" unless @reference_order;

my (@alignments, %seen_alignment);
for my $sample_dir (split /,/, $sample_dirs) {
	$sample_dir =~ s{/+$}{};
	die "Empty sample directory in --sample-dirs\n" unless length $sample_dir;
	my $marker = "$sample_dir/mapping/done.sto";
	open my $marker_fh, '<', $marker
		or die "Missing mapping completion marker $marker: $!\n";
	my $mapped_name = <$marker_fh>;
	close $marker_fh;
	chomp $mapped_name if defined $mapped_name;
	die "Invalid mapping filename in $marker\n"
		unless defined($mapped_name) && $mapped_name =~ /^[^\\\/\r\n]+\.(?:bam|cram)$/i;

	my $primary = "$sample_dir/mapping/$mapped_name";
	my $supplemental = $primary;
	$supplemental =~ s/-smd\./.sup-smd./;
	my $primary_found = 0;
	for my $candidate ($primary, $supplemental) {
		if (!-s $candidate && $candidate =~ /\.bam$/i) {
			(my $alternate = $candidate) =~ s/\.bam$/.cram/i;
			$candidate = $alternate if -s $alternate;
		} elsif (!-s $candidate && $candidate =~ /\.cram$/i) {
			(my $alternate = $candidate) =~ s/\.cram$/.bam/i;
			$candidate = $alternate if -s $alternate;
		}
		next unless -s $candidate;
		$primary_found = 1 if $candidate !~ /\.sup-smd\./;
		push @alignments, $candidate unless $seen_alignment{$candidate}++;
	}
	die "No non-empty primary BAM or CRAM named by $marker\n" unless $primary_found;
}

die "No alignment files were found\n" unless @alignments;
my @failures;
for my $alignment (@alignments) {
	my (%alignment_reference, @alignment_order);
	open my $header_fh, '-|', $samtools, 'view', '-H', $alignment
		or die "Could not run $samtools view -H for $alignment: $!\n";
	while (my $line = <$header_fh>) {
		next unless $line =~ /^\@SQ\t/;
		chomp $line;
		my ($name) = $line =~ /(?:^|\t)SN:([^\t]+)/;
		my ($length) = $line =~ /(?:^|\t)LN:(\d+)/;
		next unless defined($name) && defined($length);
		$alignment_reference{$name} = 0 + $length;
		push @alignment_order, $name;
	}
	my $header_ok = close $header_fh;
	if (!$header_ok || !@alignment_order) {
		push @failures, "$alignment: could not read a non-empty sequence dictionary";
		next;
	}

	my @differences;
	for my $name (@alignment_order) {
		if (!exists $reference{$name}) {
			push @differences, "missing '$name'";
		} elsif ($reference{$name} != $alignment_reference{$name}) {
			push @differences, "length mismatch '$name' ($alignment_reference{$name} != $reference{$name})";
		}
		last if @differences >= 5;
	}
	if (!@differences && @alignment_order != @reference_order) {
		push @differences, 'different sequence counts (' . scalar(@alignment_order)
			. ' != ' . scalar(@reference_order) . ')';
	}
	push @failures, "$alignment: " . join('; ', @differences) if @differences;
}

if (@failures) {
	die "Binning mapping/reference mismatch. These mappings must be regenerated "
		."against the current assembly ($assembly):\n  "
		.join("\n  ", @failures)."\n";
}

print "Validated ".scalar(@alignments)." mapping file(s) against $assembly\n";

