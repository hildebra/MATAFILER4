#!/usr/bin/env perl
# Interleave the two mate files of a paired FASTQ library on stdout, so that both
# mates go through one search run (DIAMOND reads its queries from stdin when -q is
# omitted, and writes hits in query order). Each read name is cut at the first
# whitespace and ends in /1 or /2, so the hit parsers find the mates of a pair next
# to each other and recognise them as mates.
use strict;
use warnings;
use Mods::GenoMetaAss qw(gzipopen);

die "Usage: $0 <R1.fq[.gz]> <R2.fq[.gz]>\n" unless @ARGV == 2;
my ($r1, $r2) = @ARGV;
my ($fh1) = gzipopen($r1, 'mate 1 reads', 1, 0);
my ($fh2) = gzipopen($r2, 'mate 2 reads', 1, 0);

#one FASTQ record (4 lines), or () at the end of the file
sub next_record {
	my ($fh, $file) = @_;
	my $head = <$fh>;
	return () unless defined $head;
	my @rec = ($head, scalar(<$fh>), scalar(<$fh>), scalar(<$fh>));
	die "Truncated FASTQ record in $file\n" if grep { !defined } @rec;
	die "Not a FASTQ record in $file: $head" unless $head =~ /^@/ && $rec[2] =~ /^\+/;
	return @rec;
}

#read name without comment and mate suffix
sub base_name {
	my ($head) = @_;
	my ($name) = $head =~ /^@(\S+)/;
	die "FASTQ record without a read name: $head" unless defined $name;
	$name =~ s/\/[12]$//;
	return $name;
}

binmode STDOUT;
my $pairs = 0;
while (1) {
	my @m1 = next_record($fh1, $r1);
	my @m2 = next_record($fh2, $r2);
	last if !@m1 && !@m2;
	die "$r1 and $r2 hold different numbers of reads (one ends after $pairs pairs)\n" unless @m1 && @m2;
	my ($n1, $n2) = (base_name($m1[0]), base_name($m2[0]));
	if ($n1 ne $n2) {
		#names that differ only in a final mate digit (e.g. SRR1.5.1 / SRR1.5.2)
		my ($b1, $b2) = ($n1, $n2);
		die "Mates out of order at pair ".($pairs + 1).": $n1 vs $n2\n"
			unless $b1 =~ s/[._:]?1$// && $b2 =~ s/[._:]?2$// && $b1 eq $b2;
		$n1 = $b1;
	}
	print "\@$n1/1\n", @m1[1 .. 3], "\@$n1/2\n", @m2[1 .. 3];
	$pairs++;
}
close($fh1) or die "Reading $r1 failed\n";
close($fh2) or die "Reading $r2 failed\n";
close(STDOUT) or die "Writing interleaved reads failed: $!\n";
