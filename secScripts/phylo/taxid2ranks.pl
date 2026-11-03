#!/usr/bin/env perl
# Expand NCBI taxonomy IDs to the seven ranks used by the gene catalogue
# (geneCat.pl krakenTax): domain, phylum, class, order, family, genus, species.
#
# Replaces secScripts/phylo/get_ranks.py (ete3 NCBITaxa). The taxonomy is read
# from the Kraken2 database whose classifications are being expanded, so the
# names always match the IDs Kraken2 reported and no ete3/NCBI download is needed:
#   1. <db>/taxo.k2d                       (Kraken2 binary taxonomy, always present)
#   2. <db>/taxonomy/nodes.dmp + names.dmp (NCBI taxdump kept by kraken2-build)
#
# Usage: taxid2ranks.pl -db <kraken2 DB dir | taxo.k2d | taxdump dir> [-in ids.txt] [taxid ...]
# Output: one line per requested ID, in input order:
#   taxid<TAB>domain<TAB>phylum<TAB>class<TAB>order<TAB>family<TAB>genus<TAB>species
# Missing ranks are empty. Unknown IDs get an all-empty lineage (and a warning),
# instead of silently disappearing as with get_ranks.py.
#
# The first column accepts NCBI's "domain" (2025 rename) and the older
# "superkingdom"; viruses, whose top node is an "acellular root" in current
# NCBI taxonomies, keep that name ("Viruses") as before.

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $db = '';
my $inFile = '';
GetOptions('db=s' => \$db, 'in=s' => \$inFile)
	or die "Usage: $0 -db <kraken2 DB dir|taxo.k2d|taxdump dir> [-in ids.txt] [taxid ...]\n";
die "$0: -db is required\n" if $db eq '';

my @ids = @ARGV;
if ($inFile ne '') {
	open my $in, '<', $inFile or die "Cannot read $inFile: $!\n";
	while (my $line = <$in>) {
		$line =~ s/[\r\n]+$//;
		push @ids, grep { $_ ne '' } split /\s+/, $line;
	}
	close $in;
}
for my $id (@ids) {
	die "$0: taxonomy ID '$id' is not a non-negative integer\n" unless $id =~ /^\d+$/;
}
exit 0 unless @ids;

my @columns = ('domain', qw(phylum class order family genus species));
my %column = map { $columns[$_] => $_ } 0..$#columns;
$column{superkingdom} = 0;
my %topFallback = ('acellular root' => 1); # Viruses in current NCBI taxonomies

my %want = map { $_ => 1 } @ids;
my $lineages = -d $db && !-e "$db/taxo.k2d" && -e "$db/nodes.dmp" ? read_taxdump($db, \%want)
	: -d $db && !-e "$db/taxo.k2d" && -e "$db/taxonomy/nodes.dmp" ? read_taxdump("$db/taxonomy", \%want)
	: read_k2d(-d $db ? "$db/taxo.k2d" : $db, \%want);

my $unknown = 0;
for my $id (@ids) {
	my @selected = ('') x @columns;
	if (my $path = $lineages->{$id}) {
		my $top = '';
		for my $node (@{$path}) {
			my ($rank, $name) = @{$node};
			$rank = lc $rank;
			if (exists $column{$rank}) {
				$selected[$column{$rank}] = $name if $selected[$column{$rank}] eq '';
			} elsif ($topFallback{$rank}) {
				$top = $name;
			}
		}
		$selected[0] = $top if $selected[0] eq '';
	} else {
		$unknown++;
	}
	print join("\t", $id, @selected), "\n";
}
warn "$0: $unknown taxonomy ID(s) not found in $db; reported with empty ranks\n" if $unknown;
exit 0;

# Kraken2 taxo.k2d (src/taxonomy.cc, Taxonomy::WriteToDisk):
#   "K2TAXDAT", size_t node_count, size_t name_data_len, size_t rank_data_len,
#   node_count x TaxonomyNode {parent_id, first_child, child_count, name_offset,
#   rank_offset, external_id, godparent_id} (7 x uint64), name data, rank data.
# Node 0 is a zeroed dummy; the root is internal ID 1. Strings are NUL-terminated.
sub read_k2d {
	my ($file, $want) = @_;
	open my $fh, '<:raw', $file or die "Cannot read Kraken2 taxonomy $file: $!\n";
	local $/;
	my $data = <$fh>;
	close $fh;
	my $magic = 'K2TAXDAT';
	die "$file is not a Kraken2 taxonomy (missing $magic header)\n"
		unless defined($data) && length($data) >= 32 && substr($data, 0, 8) eq $magic;
	my ($nodeCount, $nameLen, $rankLen) = unpack('Q<3', substr($data, 8, 24));
	my $nodeBytes = 56;
	my $nameStart = 32 + $nodeCount * $nodeBytes;
	my $rankStart = $nameStart + $nameLen;
	die "$file is truncated (expected ".($rankStart + $rankLen)." bytes, found ".length($data).")\n"
		if length($data) < $rankStart + $rankLen;

	my $string = sub {
		my ($start, $limit, $offset) = @_;
		return '' if $offset >= $limit;
		my $end = index($data, "\0", $start + $offset);
		$end = $start + $limit if $end < 0 || $end > $start + $limit;
		return substr($data, $start + $offset, $end - $start - $offset);
	};
	my %internal; # external ID -> internal ID, only for requested IDs
	for my $i (1 .. $nodeCount - 1) {
		my $external = unpack('Q<', substr($data, 32 + $i * $nodeBytes + 40, 8));
		$internal{$external} = $i if $want->{$external};
	}
	my %lineage;
	for my $external (keys %internal) {
		my @path;
		my %seen;
		my $node = $internal{$external};
		while ($node && !$seen{$node}++) {
			my ($parent, undef, undef, $nameOff, $rankOff) =
				unpack('Q<5', substr($data, 32 + $node * $nodeBytes, 40));
			push @path, [$string->($rankStart, $rankLen, $rankOff), $string->($nameStart, $nameLen, $nameOff)]
				unless $node == 1; # root
			last if $node == 1;
			$node = $parent;
		}
		$lineage{$external} = \@path;
	}
	return \%lineage;
}

# NCBI taxdump: nodes.dmp "taxid\t|\tparent\t|\trank\t|..." and names.dmp
# "taxid\t|\tname\t|\tunique name\t|\tname class\t|" (scientific names only).
sub read_taxdump {
	my ($dir, $want) = @_;
	my (%parent, %rank);
	open my $nodes, '<', "$dir/nodes.dmp" or die "Cannot read $dir/nodes.dmp: $!\n";
	while (my $line = <$nodes>) {
		my ($id, $par, $rk) = split /\t\|\t/, $line, 4;
		next unless defined $rk;
		$parent{$id} = $par; $rank{$id} = $rk;
	}
	close $nodes;
	my %need;
	for my $id (keys %{$want}) {
		my $node = $id; my %seen;
		while (exists $parent{$node} && !$seen{$node}++) {
			$need{$node} = 1;
			last if $node eq '1';
			$node = $parent{$node};
		}
	}
	my %name;
	open my $names, '<', "$dir/names.dmp" or die "Cannot read $dir/names.dmp: $!\n";
	while (my $line = <$names>) {
		next unless $line =~ /\tscientific name\t\|/;
		my ($id, $nm) = split /\t\|\t/, $line, 3;
		$name{$id} = $nm if $need{$id};
	}
	close $names;
	my %lineage;
	for my $id (keys %{$want}) {
		next unless exists $parent{$id};
		my @path; my $node = $id; my %seen;
		while (exists $parent{$node} && !$seen{$node}++ && $node ne '1') {
			push @path, [$rank{$node}, $name{$node} // ''];
			$node = $parent{$node};
		}
		$lineage{$id} = \@path;
	}
	return \%lineage;
}
