#!/usr/bin/env perl
use warnings;
use strict;

die "Usage: $0 <fasta[,fasta...]> [min-size [secondary-min-size [output]]]\n" unless @ARGV;
my ($inputs, $min_size, $secondary_min, $output) = @ARGV;
$min_size     = 500 unless defined $min_size;
$secondary_min = 200 unless defined $secondary_min;
my $custom_output = defined $output;
$output       = "$inputs.filt" unless $custom_output;
my $secondary_output = $custom_output ? "${output}2" : "$inputs.filt2";

open my $primary,   '>', $output           or die "Cannot open $output: $!\n";
open my $secondary, '>', $secondary_output or die "Cannot open $secondary_output: $!\n";

sub emit_record {
    my ($header, $sequence) = @_;
    return unless length $header;
    my $length = length $sequence;
    my $target = $length >= $min_size ? $primary
               : $secondary_min > 0 && $length >= $secondary_min ? $secondary
               : undef;
    return unless $target;
    $sequence =~ s/(.{1,80})/$1\n/gs;
    print {$target} "$header\n$sequence" or die "Cannot write filtered sequence: $!\n";
}

for my $file (split /,/, $inputs) {
    open my $in, '<', $file or die "Cannot open $file: $!\n";
    my $first = <$in>;
    next unless defined $first;
    seek $in, 0, 0 or die "Cannot rewind $file: $!\n";

    if ($first =~ /^@/) {
        while (my $header = <$in>) {
            my $sequence = <$in>;
            my $plus     = <$in>;
            my $quality  = <$in>;
            die "Truncated FASTQ record in $file\n"
                unless defined $sequence && defined $plus && defined $quality;
            die "Invalid FASTQ header in $file: $header" unless $header =~ /^@/;
            chomp($header, $sequence);
            $header =~ s/^@/>/;
            emit_record($header, $sequence);
        }
    } else {
        my ($header, $sequence) = ('', '');
        while (my $line = <$in>) {
            chomp $line;
            if ($line =~ /^>/) {
                emit_record($header, $sequence);
                ($header, $sequence) = ($line, '');
            } else {
                die "Sequence encountered before first FASTA header in $file\n" unless length $header;
                $sequence .= $line;
            }
        }
        emit_record($header, $sequence);
    }
    close $in or die "Cannot close $file: $!\n";
}

close $primary   or die "Cannot close $output: $!\n";
close $secondary or die "Cannot close $secondary_output: $!\n";
