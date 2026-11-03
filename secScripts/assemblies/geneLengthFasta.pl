#!/usr/bin/env perl
use strict;
use warnings;
use Mods::GenoMetaAss qw(gzipopen);

die "Usage: $0 <genes.fasta[.gz]> <output-lengths.tsv>\n" unless @ARGV == 2;
my ($input, $output) = @ARGV;
my ($in, $ok) = gzipopen($input, 'gene FASTA', 1);
die "Cannot open $input\n" unless $ok && $in;
#written under a temporary name: callers only test that the table exists, so a killed job must not leave a partial one
my $tmpOut = "$output.tmp.$$";
END { unlink $tmpOut if (defined($tmpOut) && -e $tmpOut); } #no leftovers after a failure
open my $out, '>', $tmpOut or die "Cannot open $tmpOut: $!\n";

my ($id, $length, $records) = ('', 0, 0);
while (my $line = <$in>) {
	$line =~ s/[\r\n]+\z//;
	if ($line =~ /^>(\S+)/) {
		print {$out} "$id\t$length\n" or die "Cannot write $output: $!\n" if length $id;
		($id, $length) = ($1, 0);
		$records++;
		next;
	}
	die "Sequence encountered before first FASTA header in $input\n" unless length $id;
	$line =~ s/\s+//g;
	$length += length($line);
}
die "No FASTA records found in $input\n" unless $records;
print {$out} "$id\t$length\n" or die "Cannot write $output: $!\n";
close $out or die "Cannot close $output: $!\n";
close $in or die "Cannot close $input: $!\n";
rename($tmpOut, $output) or die "Cannot rename $tmpOut to $output: $!\n";
print "Done\n";
