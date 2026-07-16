#!/usr/bin/env perl
use warnings;
use strict;

die "Usage: $0 <contigs.fasta> <sample-tag>\n" unless @ARGV == 2;
my ($inF, $tag) = @ARGV;
die "Input FASTA is missing or empty: $inF\n" unless -s $inF;
die "Sample tag must not be empty\n" unless length $tag;
die "Sample tag contains FASTA-delimiting whitespace or '>'\n" if $tag =~ /[\s>]/;

my $tmpOut   = "$inF.tmp";
my $transOut = "$inF.lnk";
open my $in,  '<', $inF      or die "Cannot open $inF: $!\n";
open my $out, '>', $tmpOut   or die "Cannot open $tmpOut: $!\n";
open my $lnk, '>', $transOut or die "Cannot open $transOut: $!\n";

my ($old_header, $seq) = ('', '');
my @circular;
my $count = 0;

sub write_record {
    return unless length $old_header;
    ++$count;
    my $new_header = '>' . $tag . '__C' . $count . '_L=' . length($seq) . '=';
    my $wrapped = $seq;
    $wrapped =~ s/(.{1,80})/$1\n/gs;
    print {$out} "$new_header\n$wrapped" or die "Cannot write $tmpOut: $!\n";
    print {$lnk} "$new_header\t$old_header\n" or die "Cannot write $transOut: $!\n";
    push @circular, $new_header
        if $old_header =~ /^>ctg[\d_x]+c$/
        || $old_header =~ /^>ctg[\d_x]+ .* circular=yes$/;
}

while (my $line = <$in>) {
    chomp $line;
    if ($line =~ /^>/) {
        write_record();
        $old_header = $line;
        $seq = '';
    } else {
        die "Sequence encountered before first FASTA header in $inF\n"
            unless length $old_header;
        $seq .= $line;
    }
}
write_record();

close $in  or die "Cannot close $inF: $!\n";
close $out or die "Cannot close $tmpOut: $!\n";
close $lnk or die "Cannot close $transOut: $!\n";

if (@circular) {
    open my $circ, '>', "$inF.circ" or die "Cannot open $inF.circ: $!\n";
    print {$circ} join("\n", @circular), "\n" or die "Cannot write $inF.circ: $!\n";
    close $circ or die "Cannot close $inF.circ: $!\n";
}

rename $tmpOut, $inF or die "Cannot replace $inF with $tmpOut: $!\n";
print "Done renaming $count contigs\n";
