use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use lib File::Spec->catdir($Bin, '..');
use lib File::Spec->catdir($Bin, 'lib');
use MFTestConfig;
use Test::More;

use Mods::GenoMetaAss qw(reverse_complement reverse_complement_IUPAC convertNT2AA
	lcp is_integer parse_duration quantile median readFasta writeFasta readFastHD
	prefixFAhd readTabByKey filsizeMB fileGZe fileGZs gzipopen);

sub write_file {
	my ($path, $contents) = @_;
	open my $fh, '>', $path or die "Cannot write $path: $!";
	print {$fh} $contents;
	close $fh or die "Cannot close $path: $!";
}

sub slurp {
	my ($path) = @_;
	open my $fh, '<', $path or die "Cannot read $path: $!";
	local $/;
	return <$fh>;
}

my $tmp = tempdir(CLEANUP => 1);

# --- sequence transforms ---------------------------------------------------
is(reverse_complement('AACGTt'), 'aACGTT', 'reverse complement keeps case');
is(reverse_complement('ACGTN'), 'NACGT', 'reverse complement leaves N untouched');
is(reverse_complement(reverse_complement('GATTACA')), 'GATTACA',
	'reverse complement is its own inverse');
is(reverse_complement_IUPAC('ACGTRYKMBDHVNSW'), 'WSNBDHVKMRYACGT',
	'IUPAC reverse complement maps every ambiguity code to its complement');
is(reverse_complement_IUPAC('kmry'), 'rykm', 'IUPAC reverse complement handles lower case');
my $iupac = 'ACGTRYKMBDHVNSW';
is(reverse_complement_IUPAC(reverse_complement_IUPAC($iupac)), $iupac,
	'IUPAC reverse complement is its own inverse');

{
	my @warnings;
	local $SIG{__WARN__} = sub { push @warnings, @_ };
	my $stdout = '';
	open my $capture, '>', \$stdout or die $!;
	my $old = select $capture;
	is(convertNT2AA('ATGAAATAG'), 'MK*', 'standard codons are translated with stop as *');
	is(convertNT2AA('atgaaa'), 'MK', 'lower-case codons are translated');
	is(convertNT2AA('ATGNNNTGG'), 'MXW', 'codons with ambiguous bases become X');
	is(convertNT2AA('ATGAANNN'), 'M', 'trailing N padding is removed and a partial codon is dropped');
	select $old;
	is_deeply(\@warnings, [], 'translation of ambiguous codons emits no warnings');
}

# --- small string helpers --------------------------------------------------
is(lcp('sample_R1.fq', 'sample_R2.fq'), 'sample_R', 'lcp returns the longest common prefix');
is(lcp('abc', 'xyz'), '', 'lcp of unrelated strings is empty');
is(lcp('only'), 'only', 'lcp of one string is the string');
ok(is_integer('42') && is_integer('-3') && is_integer('+7'), 'signed integers are recognised');
ok(!is_integer('4.2') && !is_integer('4e3') && !is_integer('') && !is_integer(undef),
	'decimals, exponents, empty and undefined values are not integers');
is(parse_duration(59), '00:00:59', 'durations below a minute are zero-padded');
is(parse_duration(3725), '01:02:05', 'durations are formatted as hh:mm:ss');
is(parse_duration(90061), '25:01:01', 'durations above one day keep counting hours');

# --- legacy GenoMetaAss statistics ----------------------------------------
is(quantile(0.5, 1 .. 5), 3, 'GenoMetaAss::quantile rounds the rank position');
is(quantile(0.9, 7), 7, 'GenoMetaAss::quantile of one value returns it');
ok(!defined(quantile(0.5)), 'GenoMetaAss::quantile of no values is undefined');
ok(!eval { quantile(2, 1, 2); 1 }, 'GenoMetaAss::quantile rejects fractions above 1');
is(median(1, 3), 2, 'GenoMetaAss::median of an even list averages the middle values');
is(median(), 0, 'GenoMetaAss::median keeps its legacy empty result of zero');

# --- FASTA IO --------------------------------------------------------------
my $fasta = File::Spec->catfile($tmp, 'records.fna');
write_file($fasta, ">g2 second gene\nCCCC\nGG\n>g1\nAAAA\n>g3 third\nTTTT\n");
is_deeply(readFasta($fasta, 1), { g1 => 'AAAA', g2 => 'CCCCGG', g3 => 'TTTT' },
	'readFasta joins wrapped sequence lines and shortens headers on request');
is_deeply(readFastHD($fasta), [qw(g2 g1 g3)], 'readFastHD returns short headers in file order');
is_deeply(readFastHD($fasta, 1), ['g2 second gene', 'g1', 'g3 third'],
	'readFastHD can return complete header lines');

my $written = File::Spec->catfile($tmp, 'written.fna');
writeFasta({ z => 'GG', a => 'CC', m => 'TT' }, $written);
is(slurp($written), ">a\nCC\n>m\nTT\n>z\nGG\n",
	'writeFasta output is sorted by identifier and therefore reproducible');
writeFasta({ '>x' => 'AC', y => 'GT' }, $written);
is(slurp($written), ">x\nAC\n>y\nGT\n", 'writeFasta accepts identifiers with or without >');
writeFasta({ a => 'A', b => 'C', c => 'G' }, $written, 2);
is(slurp($written), ">a\nA\n>b\nC\n", 'writeFasta writes at most the requested number of records');
is_deeply(readFasta($written), { a => 'A', b => 'C' }, 'written FASTA reads back unchanged');

my $prefixed = prefixFAhd({ '>g1' => 'AA', g2 => 'CC' }, 'S1');
is_deeply($prefixed, { 'S1.g1' => 'AA', 'S1.g2' => 'CC' },
	'prefixFAhd prefixes identifiers with and without a leading >');

# --- compressed-file helpers ----------------------------------------------
my $plain = File::Spec->catfile($tmp, 'table.tsv');
write_file($plain, "k1\tv1\nk2\tv2\n");
is(fileGZe($plain), 1, 'fileGZe finds a nonempty plain file');
is(fileGZe("$plain.gz"), 1, 'fileGZe falls back from a .gz name to the plain file');
is(fileGZe(File::Spec->catfile($tmp, 'missing.tsv')), 0, 'fileGZe reports a missing file');
my $empty = File::Spec->catfile($tmp, 'empty.tsv');
write_file($empty, '');
is(fileGZe($empty), 0, 'fileGZe treats an empty file as absent');
is(fileGZs($plain), -s $plain, 'fileGZs returns the plain file size');

my $gz = File::Spec->catfile($tmp, 'gz_table.tsv');
write_file($gz, "a\t1\nb\t2\nc\t3\n");
system('gzip', '-f', $gz) == 0 or die "gzip failed";
is(fileGZs($gz), 5 * (-s "$gz.gz"), 'fileGZs estimates uncompressed size of a gzip-only file');
my %table = readTabByKey("$gz.gz");
is_deeply(\%table, { a => 1, b => 2, c => 3 }, 'readTabByKey reads gzip-compressed two-column tables');
my ($handle, $ok) = gzipopen("$gz.gz", 'test table');
ok($ok, 'gzipopen opens gzip input');
my @lines = <$handle>;
close $handle;
is(scalar(@lines), 3, 'gzipopen streams all decompressed lines');

my $mb_file = File::Spec->catfile($tmp, 'one_mb.bin');
write_file($mb_file, 'x' x (1024 * 1024));
is(filsizeMB($mb_file), 1, 'filsizeMB reports sizes in MiB');
is(filsizeMB($tmp, 'one_mb.bin', 'one_mb.bin'), 2, 'filsizeMB resolves relative names and sums them');
is(filsizeMB(), 0, 'filsizeMB of nothing is zero');

done_testing();
