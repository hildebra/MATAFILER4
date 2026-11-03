use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use lib File::Spec->catdir($Bin, '..');
use lib File::Spec->catdir($Bin, 'lib');
use MFTestConfig;
use Test::More;

use Mods::phyloTools qw(fixHDs4Phylo);
use Mods::MGSLocus qw(protein_kmer_similarity);
use Mods::MosaicLoci qw(read_paf_stream);
use Mods::GenoMetaAss qw(readFasta readFastHD);

sub write_file {
	my ($path, $contents) = @_;
	open my $fh, '>', $path or die "Cannot write $path: $!";
	print {$fh} $contents;
	close $fh or die "Cannot close $path: $!";
}

sub quietly {
	my ($code) = @_;
	my $sink = '';
	open my $capture, '>', \$sink or die $!;
	my $old = select $capture;
	my $result = eval { $code->() };
	my $error = $@;
	select $old;
	die $error if $error;
	return $result;
}

my $tmp = tempdir(CLEANUP => 1);

# --- fixHDs4Phylo -----------------------------------------------------------
my $clean = File::Spec->catfile($tmp, 'clean.faa');
write_file($clean, ">s1|g1\nMKV\n>s2|g1\nMKI\n");
is(quietly(sub { fixHDs4Phylo($clean) }), $clean, 'headers without reserved characters are used unchanged');
ok(!-e "$clean.fix", 'no rewritten copy is created when nothing needs fixing');
is(fixHDs4Phylo(''), '', 'an empty input name is passed through');

my $dirty = File::Spec->catfile($tmp, 'dirty.faa');
write_file($dirty, ">b:sample(1)\nMKV\n>a,sample;[2]\nMKI\n");
my $fixed = quietly(sub { fixHDs4Phylo($dirty) });
is($fixed, "$dirty.fix", 'headers with reserved characters are rewritten to a .fix copy');
is_deeply(readFastHD($fixed), ['a|sample||2|1', 'b|sample|1|1'],
	'reserved characters become | and names are written in sorted order with a counter');
is_deeply([sort values %{ readFasta($fixed) }], ['MKI', 'MKV'], 'sequences are preserved');

my $long = File::Spec->catfile($tmp, 'long.faa');
my $prefix = 'x' x 40;
write_file($long, ">${prefix}:one\nAA\n>${prefix}:two\nCC\n");
my $long_fixed = quietly(sub { fixHDs4Phylo($long) });
is_deeply(readFastHD($long_fixed), ["${prefix}1", "${prefix}2"],
	'names truncated to 40 characters stay unique through the counter');
my $long_nt = File::Spec->catfile($tmp, 'long.fna');
write_file($long_nt, ">${prefix}:two\nCCCCCC\n>${prefix}:one\nAAAAAA\n");
my $nt_fixed = quietly(sub { fixHDs4Phylo($long_nt) });
is(readFasta($long_fixed)->{"${prefix}1"}, 'AA', 'protein record one keeps its renamed identifier');
is(readFasta($nt_fixed)->{"${prefix}1"}, 'AAAAAA',
	'paired nucleotide input is renamed identically regardless of record order');

# --- protein_kmer_similarity -----------------------------------------------
is(protein_kmer_similarity('MKVLAAG', 'MKVLAAG'), 1, 'identical proteins have similarity 1');
is(protein_kmer_similarity('MKVLAAG', 'WWWWWWW'), 0, 'unrelated proteins have similarity 0');
is(protein_kmer_similarity('mkvl-aag*', 'MKVLAAG'), 1, 'case, gaps and stops are ignored');
is(protein_kmer_similarity('MKV', 'MKV'), 0, 'sequences shorter than k have similarity 0');
is(protein_kmer_similarity(undef, 'MKVL'), 0, 'undefined input has similarity 0');
my $partial = protein_kmer_similarity('ABCDEF', 'ABCDXY');
is($partial, 2 * 1 / (3 + 3), 'similarity is the Dice coefficient of distinct k-mers');
is(protein_kmer_similarity('ABCDEF', 'ABCDXY', 2), 2 * 3 / (5 + 5), 'k can be chosen explicitly');

# --- read_paf_stream --------------------------------------------------------
my $paf = join('', map { join("\t", @{$_}) . "\n" } (
	# query qlen qs qe strand target tlen ts te matches alnlen mapq
	[qw(q1 100 0 100 + t1 200 0 100 90 100 60)],
	[qw(q1 100 0 100 + t1 200 50 150 95 100 60)],   # better duplicate for q1/t1
	[qw(q1 100 0 50 + t2 100 0 50 49 50 10)],       # low query coverage
	[qw(q2 100 0 100 - q2 100 0 100 100 100 60)],   # self alignment
	[qw(q3 100 0 100 + t3 100 0 100 70 100 60)],    # low identity
));
open my $fh, '<', \$paf or die $!;
my %statistics;
my $hits = read_paf_stream($fh, 'fixture', {
	exclude_self => 1, minimum_identity => 0.8, minimum_query_coverage => 0.8,
	statistics => \%statistics,
});
is_deeply([sort keys %{$hits}], ['q1'], 'only queries with passing alignments are returned');
is(scalar(@{ $hits->{q1} }), 1, 'duplicate query-target alignments collapse to one');
is($hits->{q1}[0]{matches}, 95, 'the alignment with most matches is retained');
is($hits->{q1}[0]{identity}, 0.95, 'identity is matches over alignment length');
is($hits->{q1}[0]{target_coverage}, 0.5, 'target coverage uses the target length');
is_deeply(\%statistics, {
	raw_alignments => 5, self_alignments_filtered => 1, threshold_alignments_filtered => 2,
	duplicate_alignments_replaced => 1, queries_with_retained_alignments => 1,
	retained_alignments => 1,
}, 'filter statistics account for every input alignment');

my $multi = join('', map { join("\t", @{$_}) . "\n" } (
	[qw(q1 100 0 100 + tB 100 0 100 90 100 60)],
	[qw(q1 100 0 100 + tA 100 0 100 90 100 60)],
));
open my $multi_fh, '<', \$multi or die $!;
my $multi_hits = read_paf_stream($multi_fh, 'fixture');
is_deeply([map { $_->{target} } @{ $multi_hits->{q1} }], ['tA', 'tB'],
	'several targets per query are returned in a reproducible (sorted) order');

my $bad = "q1\t100\t0\n";
open my $bad_fh, '<', \$bad or die $!;
ok(!eval { read_paf_stream($bad_fh, 'broken.paf'); 1 }, 'truncated PAF rows are rejected');
like($@, qr/Malformed PAF row 1 in broken\.paf/, 'the error names the row and source');

done_testing();
