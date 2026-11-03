#!/usr/bin/env perl
# Read counts per taxon from Kraken2 per-read output, for a series of confidence thresholds.
#
#   kraken2 --db DB --confidence 0 --output - ... | krak2_count_tax.pl -db DB -thresholds 0.01,0.1 -out prefix
#
# Kraken2 has no kraken-filter: the reads are classified once (confidence 0), and the call at each
# threshold t is recomputed as kraken2 does it (classify.cc ResolveTree): starting from the
# confidence-0 call, the taxon's own k-mer hits, then the hits of its clade, must reach
# ceil(t * k-mers of the read); otherwise the call moves to the parent.
# The taxonomy is the database's own (taxo.k2d), so it always matches the reported taxon IDs.
# Writes <prefix>.<t>.cnt.tax: "d__<domain>;<phylum>;<class>;<order>;<family>;<genus>;<species><TAB>reads",
# missing ranks as "?"; reads without a domain-level call are not counted.
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use POSIX qw(ceil);

my ($db, $thresholdList, $outPrefix) = ('', '', '');
GetOptions('db=s' => \$db, 'thresholds=s' => \$thresholdList, 'out=s' => \$outPrefix)
	or die "Error in command line arguments\n";
die "Usage: $0 -db <kraken2 DB dir|taxo.k2d> -thresholds t1,t2,... -out <prefix> [kraken2 output ...]\n"
	unless $db ne '' && $thresholdList ne '' && $outPrefix ne '';
my @thresholds = split /,/, $thresholdList;
for my $t (@thresholds) { die "Invalid confidence threshold $t\n" unless $t =~ /^(?:0|1|0?\.\d+)$/; }

# Kraken2 taxo.k2d (src/taxonomy.cc, Taxonomy::WriteToDisk; also read by secScripts/phylo/taxid2ranks.pl):
#   "K2TAXDAT", size_t node_count, size_t name_data_len, size_t rank_data_len,
#   node_count x TaxonomyNode {parent_id, first_child, child_count, name_offset,
#   rank_offset, external_id, godparent_id} (7 x uint64), name data, rank data.
# Node 0 is a zeroed dummy; the root is internal ID 1. Strings are NUL-terminated.
my $taxoF = -d $db ? "$db/taxo.k2d" : $db;
open my $fh, '<:raw', $taxoF or die "Cannot read Kraken2 taxonomy $taxoF: $!\n";
my $data = do { local $/; <$fh> };
close $fh;
die "$taxoF is not a Kraken2 taxonomy\n" unless defined($data) && length($data) >= 32 && substr($data, 0, 8) eq 'K2TAXDAT';
my ($nodeCount, $nameLen, $rankLen) = unpack('Q<3', substr($data, 8, 24));
my $nameStart = 32 + $nodeCount * 56;
my $rankStart = $nameStart + $nameLen;
die "$taxoF is truncated\n" if length($data) < $rankStart + $rankLen;
my (%internal, @external);
for my $i (1 .. $nodeCount - 1) {
	$external[$i] = unpack('Q<', substr($data, 32 + $i * 56 + 40, 8));
	$internal{$external[$i]} = $i;
}
sub k2d_string {
	my ($start, $limit, $offset) = @_;
	return '' if $offset >= $limit;
	my $end = index($data, "\0", $start + $offset);
	$end = $start + $limit if $end < 0 || $end > $start + $limit;
	return substr($data, $start + $offset, $end - $start - $offset);
}
my (%parentCache, %rankCache, %nameCache);
sub parent_of { #external taxon ID of the parent; 0 above the root
	my ($t) = @_;
	return $parentCache{$t} //= do {
		my $i = $internal{$t} // 0;
		my $p = $i > 1 ? unpack('Q<', substr($data, 32 + $i * 56, 8)) : 0;
		$p ? $external[$p] : 0;
	};
}
sub node_info {
	my ($t) = @_;
	my $i = $internal{$t};
	my (undef, undef, undef, $nameOff, $rankOff) = unpack('Q<5', substr($data, 32 + $i * 56, 40));
	$rankCache{$t} = k2d_string($rankStart, $rankLen, $rankOff);
	($nameCache{$t} = k2d_string($nameStart, $nameLen, $nameOff)) =~ s/\s+/_/g;
}

my %clade; #taxon -> {ancestor or self => 1}
sub ancestors_of {
	my ($t) = @_;
	return $clade{$t} //= do {
		my %set; my $x = $t; my $guard = 0;
		while ($x) { $set{$x} = 1; $x = parent_of($x); die "Taxonomy loop at $t\n" if ++$guard > 1000; }
		\%set;
	};
}
my @columns = qw(phylum class order family genus species);
my %lineageCache;
sub lineage_string {
	my ($t) = @_;
	return $lineageCache{$t} //= do {
		my (%byRank, $top); my $x = $t;
		while ($x && $x != 1) {
			node_info($x) unless exists $rankCache{$x};
			my $rank = lc $rankCache{$x};
			$rank = 'domain' if $rank eq 'superkingdom';
			$byRank{$rank} //= $nameCache{$x};
			$top = $nameCache{$x} if $rank eq 'acellular root'; #Viruses in current NCBI taxonomies
			$x = parent_of($x);
		}
		my $domain = $byRank{domain} // $top;
		defined($domain) ? join(';', "d__$domain", map { $byRank{$_} // '?' } @columns) : '';
	};
}

my %counts; #threshold -> lineage -> reads
$counts{$_} = {} for @thresholds;
my ($reads, $classified) = (0, 0);
while (my $line = <>) {
	$line =~ s/[\r\n]+$//;
	my @f = split /\t/, $line;
	die "Unexpected Kraken2 output line (5 columns expected): $line\n" unless @f >= 5;
	$reads++;
	next unless $f[0] eq 'C';
	my $call0 = $f[2];
	$call0 = $1 if $call0 =~ /\(taxid (\d+)\)\s*$/; #--use-names output
	die "Taxon $call0 is not in $taxoF\n" unless exists $internal{$call0};
	$classified++;
	#k-mer string: "taxid:count" runs, A = ambiguous, 0 = no hit, |:| = mate border (not a k-mer)
	my %hits; my $total = 0;
	for my $token (split / /, $f[4]) {
		next if $token eq '|:|' || $token eq '';
		my ($t, $n) = split /:/, $token;
		$total += $n;
		$hits{$t} += $n if $t ne 'A' && $t ne '0';
	}
	for my $thr (@thresholds) {
		my $required = ceil($thr * $total);
		my $call = $call0;
		my $score = $hits{$call} // 0;
		while ($call && $score < $required) {
			$score = 0;
			for my $t (keys %hits) { $score += $hits{$t} if ancestors_of($t)->{$call}; }
			last if $score >= $required;
			$call = parent_of($call);
		}
		next unless $call;
		my $lin = lineage_string($call);
		$counts{$thr}{$lin}++ if $lin ne '';
	}
}

for my $thr (@thresholds) {
	my $outF = "$outPrefix.$thr.cnt.tax";
	open my $out, '>', "$outF.tmp" or die "Cannot write $outF.tmp: $!\n";
	print {$out} "$_\t$counts{$thr}{$_}\n" for sort keys %{$counts{$thr}};
	close $out or die "Cannot write $outF.tmp: $!\n";
	rename("$outF.tmp", $outF) or die "Cannot rename $outF.tmp: $!\n";
}
print "Kraken2 reads: $reads, classified at confidence 0: $classified\n";
