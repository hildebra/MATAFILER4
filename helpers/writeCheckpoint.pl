#!/usr/bin/env perl

use strict;
use warnings;
use FindBin qw($RealBin);
use lib "$RealBin/..";
use Getopt::Long qw(GetOptions);
use Mods::Checkpoint qw(write_checkpoint);

my $stone = '';
my @outputs;
my @parameter_pairs;
my $help = 0;

GetOptions(
	'stone=s'  => \$stone,
	'output=s@' => \@outputs,
	'param=s@'  => \@parameter_pairs,
	'help|h'    => \$help,
) or die "Try --help\n";

if ($help) {
	print "Usage: $0 --stone FILE [--output FILE ...] [--param KEY=VALUE ...]\n";
	exit 0;
}
die "--stone is required\n" unless length $stone;
die "Unexpected positional arguments: @ARGV\n" if @ARGV;

my %parameters;
for my $pair (@parameter_pairs) {
	die "Invalid --param '$pair' (expected KEY=VALUE)\n" unless $pair =~ /^([^=]+)=(.*)$/;
	$parameters{$1} = $2;
}

write_checkpoint($stone, parameters => \%parameters, outputs => \@outputs);
