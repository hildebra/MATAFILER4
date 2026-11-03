#!/usr/bin/perl
# Average each gene's k-mer features over +/- radius neighbouring genes.
# Contig boundaries truncate the window; every input gene is emitted once.
use strict;
use warnings;
use Mods::GenoMetaAss qw(gzipopen gzipwrite);
use Mods::math qw(avgArray roundAr);

my ($inK, $radius) = @ARGV;
die "Usage: $0 input.4kmer.gz radius (nonnegative integer)\n"
	unless defined($inK) && defined($radius) && $radius =~ /^\d+$/;
(my $outF = $inK) =~ s/4kmer\.gz$/4kmer.pm$radius.gz/
	or die "Input filename must end in 4kmer.gz\n";
my ($input, $ok) = gzipopen($inK, "K-mer per gene", 1);
my ($output) = gzipwrite($outF, "averaged gene k-mers");
my (@rows, @genes);
my $next = 0;
my $contig = "";
while (my $line = <$input>) {
	chomp $line;
	next if $line eq "";
	my ($gene, @features) = split /\t/, $line;
	if ($gene eq "Contig") {
		print {$output} "$line\n";
		next;
	}
	my ($current) = $gene =~ /^(.*)_\d+$/;
	die "Cannot find contig info for $gene\n" unless defined $current;
	if ($current ne $contig) {
		emit_ready(1);
		@rows = (); @genes = (); $next = 0;
		$contig = $current;
	}
	push @rows, \@features;
	push @genes, $gene;
	emit_ready(0);
}
emit_ready(1);
close $input or die "Cannot close $inK: $!\n";
close $output or die "Cannot close $outF: $!\n";

sub emit_ready {
	my ($flush) = @_;
	while ($next < @rows && ($flush || $next + $radius < @rows)) {
		my $first = $next > $radius ? $next - $radius : 0;
		my $last = $next + $radius;
		$last = $#rows if $last > $#rows;
		my $mean = roundAr(avgArray([@rows[$first .. $last]]), 2);
		print {$output} join("\t", $genes[$next], @$mean), "\n";
		$next++;
		# Retain only the left context needed by the next output gene.
		my $discard = $next - $radius;
		if ($discard > 0) {
			splice @rows, 0, $discard;
			splice @genes, 0, $discard;
			$next -= $discard;
		}
	}
}
