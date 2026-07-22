use strict;
use warnings;

use File::Spec;
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::PhyloAlignment qw(filter_alignment_by_overlap);

my %nt = (
	a => 'AC-NR',
	b => 'AT-NN',
	c => 'A--NY',
);
my ($nt_filtered, $nt_length, $nt_removed) =
	filter_alignment_by_overlap(\%nt, 0, 2);
is_deeply(
	$nt_filtered,
	{ a => 'ACR', b => 'ATN', c => 'A-Y' },
	'nucleotide overlap counts IUPAC calls but not gaps or N',
);
is($nt_length, 3, 'nucleotide overlap reports retained columns');
is($nt_removed, 2, 'nucleotide overlap reports removed columns');

my %aa = (
	a => 'M-XQ',
	b => 'M.AQ',
	c => 'M?XQ',
);
my ($aa_filtered, $aa_length) =
	filter_alignment_by_overlap(\%aa, 1, 2);
is_deeply(
	$aa_filtered,
	{ a => 'MQ', b => 'MQ', c => 'MQ' },
	'amino-acid overlap treats X, gap, dot, and question mark as missing',
);
is($aa_length, 2, 'amino-acid overlap reports retained columns');

my ($unchanged, $unchanged_length, $unchanged_removed) =
	filter_alignment_by_overlap(\%nt, 0, 0);
is_deeply($unchanged, \%nt, 'zero overlap threshold leaves alignment unchanged');
is($unchanged_length, 5, 'unchanged alignment length is reported');
is($unchanged_removed, 0, 'unchanged alignment removes no columns');

my ($empty, $empty_length) =
	filter_alignment_by_overlap({ a => 'A', b => '-' }, 0, 3);
is_deeply($empty, { a => '', b => '' },
	'a threshold above taxon count yields an empty locus');
is($empty_length, 0, 'empty filtered locus reports zero length');

eval { filter_alignment_by_overlap({ a => 'AA', b => 'A' }, 0, 1) };
like($@, qr/unequal lengths/, 'unequal input alignment lengths fail explicitly');

done_testing();
