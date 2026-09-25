use strict;
use warnings;
use Test::More;
use File::Temp qw(tempdir);
use File::Path qw(make_path);
use File::Spec;
use FindBin qw($Bin);

# secScripts/phylo/taxid2ranks.pl replaces the ete3-based get_ranks.py used by
# geneCat.pl krakenTax. It must read Kraken2's taxo.k2d and an NCBI taxdump.

my $root = File::Spec->rel2abs("$Bin/..");
my $script = "$root/secScripts/phylo/taxid2ranks.pl";
my $tmp = tempdir(CLEANUP => 1);

# [external id, parent external id, rank, name]; 2025 NCBI style ("domain",
# "acellular root") plus an old-style "superkingdom" branch.
my @taxa = (
	[1, 0, 'no rank', 'root'],
	[131567, 1, 'cellular root', 'cellular organisms'],
	[2, 131567, 'domain', 'Bacteria'],
	[1224, 2, 'phylum', 'Pseudomonadota'],
	[1236, 1224, 'class', 'Gammaproteobacteria'],
	[91347, 1236, 'order', 'Enterobacterales'],
	[543, 91347, 'family', 'Enterobacteriaceae'],
	[561, 543, 'genus', 'Escherichia'],
	[562, 561, 'species', 'Escherichia coli'],
	[83333, 562, 'strain', 'Escherichia coli K-12'],
	[2157, 131567, 'superkingdom', 'Archaea'],
	[10239, 1, 'acellular root', 'Viruses'],
	[2731341, 10239, 'realm', 'Duplodnaviria'],
);

sub write_k2d {
	my ($file) = @_;
	my %internal = map { $taxa[$_][0] => $_ + 1 } 0..$#taxa; # node 0 is the dummy
	my ($names, $ranks) = ('', '');
	my @nodes = (pack('Q<7', (0) x 7));
	for my $t (@taxa) {
		my ($ext, $par, $rank, $name) = @{$t};
		my ($nameOff, $rankOff) = (length($names), length($ranks));
		$names .= "$name\0"; $ranks .= "$rank\0";
		push @nodes, pack('Q<7', $par ? $internal{$par} : 0, 0, 0, $nameOff, $rankOff, $ext, 0);
	}
	open my $fh, '>:raw', $file or die "$file: $!";
	print {$fh} 'K2TAXDAT', pack('Q<3', scalar(@nodes), length($names), length($ranks)),
		@nodes, $names, $ranks;
	close $fh or die "$file: $!";
}
sub write_taxdump {
	my ($dir) = @_;
	make_path($dir);
	open my $n, '>', "$dir/nodes.dmp" or die;
	open my $m, '>', "$dir/names.dmp" or die;
	for my $t (@taxa) {
		my ($ext, $par, $rank, $name) = @{$t};
		print {$n} join("\t|\t", $ext, $par || 1, $rank, 'XX', 0), "\t|\n";
		print {$m} join("\t|\t", $ext, $name, '', 'scientific name'), "\t|\n";
		print {$m} join("\t|\t", $ext, "$name synonym", '', 'synonym'), "\t|\n";
	}
	close $n; close $m;
}
sub run {
	my (@args) = @_;
	my $out = qx{"$^X" "$script" @args 2>&1};
	return ($? >> 8, $out);
}

my $expected = join('', map { "$_\n" }
	"562\tBacteria\tPseudomonadota\tGammaproteobacteria\tEnterobacterales\tEnterobacteriaceae\tEscherichia\tEscherichia coli",
	"83333\tBacteria\tPseudomonadota\tGammaproteobacteria\tEnterobacterales\tEnterobacteriaceae\tEscherichia\tEscherichia coli",
	"2157\tArchaea\t\t\t\t\t\t",
	"2731341\tViruses\t\t\t\t\t\t",
	"999999\t\t\t\t\t\t\t",
);

my $k2 = "$tmp/k2db";
make_path($k2);
write_k2d("$k2/taxo.k2d");
my ($status, $out) = run('-db', $k2, qw(562 83333 2157 2731341 999999));
is($status, 0, 'taxo.k2d lookup succeeds');
(my $stdout = $out) =~ s/^.*not found.*\n//m;
is($stdout, $expected, 'taxo.k2d: seven ranks, domain or superkingdom first, unknown IDs kept with empty ranks');
like($out, qr/1 taxonomy ID\(s\) not found/, 'unknown IDs are reported');

open my $ids, '>', "$tmp/ids.txt" or die; print {$ids} "562\n83333\n2157 2731341\n999999\n"; close $ids;
($status, $out) = run('-db', "$k2/taxo.k2d", '-in', "$tmp/ids.txt");
($stdout = $out) =~ s/^.*not found.*\n//m;
is($stdout, $expected, 'IDs can be read from a file and the DB given as the taxo.k2d path');

my $dump = "$tmp/dumpdb";
write_taxdump("$dump/taxonomy");
($status, $out) = run('-db', $dump, qw(562 83333 2157 2731341 999999));
($stdout = $out) =~ s/^.*not found.*\n//m;
is($stdout, $expected, 'an NCBI taxdump in <db>/taxonomy/ gives the same lineages');

open my $bad, '>:raw', "$tmp/bad.k2d" or die; print {$bad} 'NOTKRAKEN' . ("\0" x 40); close $bad;
($status, $out) = run('-db', "$tmp/bad.k2d", 562);
isnt($status, 0, 'a file without the K2TAXDAT header is rejected');
($status, $out) = run('-db', $k2, 'abc');
isnt($status, 0, 'non-numeric IDs are rejected');

done_testing;
