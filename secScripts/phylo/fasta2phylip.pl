#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std qw(getopts);
use File::Temp qw(tempfile);

# Relaxed sequential PHYLIP; retain the existing 50-character name limit.
my %options;
my $usage = "Usage: $0 [-h] [-v] [-c 0..50] [aligned.fasta]\n"
	."Reads STDIN when no file is given. Names longer than 50 characters use\n"
	."the first -c characters (default 7) and the remaining characters from the end.\n";
getopts('hvc:', \%options) or die $usage;
if ($options{h}) { print $usage; exit 0; }
my $front = $options{c} // 7;
die "-c must be an integer between 0 and 50\n"
	unless $front =~ /^\d+$/ && $front <= 50;
my ($spool) = tempfile(UNLINK => 1);
my ($name, $sequence, $length, $taxa) = ('', '', undef, 0);
my %names;
my $flush = sub {
	return unless length($name);
	$length //= length($sequence);
	die "Unequal alignment lengths for '$name': ".length($sequence).", expected $length\n"
		unless length($sequence) == $length;
	die "Empty alignment sequence for '$name'\n" unless $length;
	print {$spool} "$name $sequence\n" or die "Cannot spool PHYLIP: $!\n";
	$taxa++;
};
while (<>) {
	s/[\r\n]+\z//;
	s/^\s+//;
	die "Empty FASTA identifier\n" if /^>\s*$/;
	next unless /\S/;
	if (/^>\s*(\S+)/) {
		$flush->();
		my $original = $1;
		$name = length($original) <= 50 ? $original
			: substr($original, 0, $front).($front == 50 ? '' : substr($original, -(50 - $front)));
		die "Duplicate PHYLIP name '$name' (input '$original'); shorten input identifiers uniquely\n"
			if $names{$name}++;
		print STDERR "$original => $name\n" if $options{v} && $name ne $original;
		$sequence = '';
	} else {
		die "Sequence data before a FASTA header\n" unless length($name);
		s/\s+//g;
		$sequence .= $_;
	}
}
$flush->();
die "No aligned sequences found\n" unless $taxa;
seek($spool, 0, 0) or die "Cannot rewind PHYLIP: $!\n";
print "$taxa   $length\n" or die "Cannot write PHYLIP header: $!\n";
while (<$spool>) { print or die "Cannot write PHYLIP: $!\n"; }
close $spool or die "Cannot close PHYLIP spool: $!\n";
