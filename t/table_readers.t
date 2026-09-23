use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use lib File::Spec->catdir($Bin, '..');
use lib File::Spec->catdir($Bin, 'lib');
use MFTestConfig;
use Test::More;

use Mods::TamocFunc qw(uniq readTabbed readTable);
use Mods::geneCat qw(read_matrix readGeneIdx);
use Mods::Binning qw(readMGS readMGSrev);
use Mods::IO_Tamoc_progs qw(convert2Gb greaterComputeSpace);

sub write_file {
	my ($path, $contents) = @_;
	open my $fh, '>', $path or die "Cannot write $path: $!";
	print {$fh} $contents;
	close $fh or die "Cannot close $path: $!";
}

# Several readers print progress; keep TAP output clean.
sub quietly {
	my ($code) = @_;
	my $sink = '';
	open my $capture, '>', \$sink or die $!;
	my $old = select $capture;
	my @result = eval { $code->() };
	my $error = $@;
	select $old;
	die $error if $error;
	return wantarray ? @result : $result[0];
}

my $tmp = tempdir(CLEANUP => 1);

# --- TamocFunc -------------------------------------------------------------
is_deeply([uniq(qw(b a b c a))], [qw(b a c)], 'uniq keeps the first occurrence in order');

my $two_col = File::Spec->catfile($tmp, 'two.tsv');
write_file($two_col, "g1\tK01\ng2\tK02\ng1\tK03\n");
is_deeply(readTabbed($two_col), { g1 => 'K03', g2 => 'K02' },
	'readTabbed maps column 1 to column 2, later rows replacing earlier ones');

my $table = File::Spec->catfile($tmp, 'table.tsv');
write_file($table, "#comment\nid1\ta\t\tc\nid2\tx\n");
is_deeply(readTable($table, "\t"), { id1 => ['a', '', 'c'], id2 => ['x'] },
	'readTable skips comments and keeps empty trailing fields');
is_deeply(readTable($table, "\t", ';'), { id1 => 'a;;c', id2 => 'x' },
	'readTable can merge the remaining columns into one string');

# --- geneCat matrix and index readers --------------------------------------
my $matrix = File::Spec->catfile($tmp, 'Matrix.mat');
write_file($matrix, "Gene\tS1\tS2\n1\t5\t0\n2\t0\t3\n3\t1\t1\n");
my $all = quietly(sub { read_matrix($matrix) });
is_deeply($all->{header}, ['S1', 'S2'], 'read_matrix returns the sample header');
is_deeply($all->{2}, [0, 3], 'read_matrix returns each row by gene ID');
my $subset = quietly(sub { read_matrix($matrix, "\t", { 1 => 1, 3 => 1 }) });
is_deeply([sort grep { $_ ne 'header' } keys %{$subset}], [1, 3],
	'read_matrix subset keeps only requested genes');
is_deeply($subset->{header}, ['S1', 'S2'],
	'read_matrix keeps the header even when a subset is requested');
my $short = File::Spec->catfile($tmp, 'short.mat');
write_file($short, "Gene\tS1\n1\t5\n");
ok(!eval { quietly(sub { read_matrix($short) }); 1 },
	'read_matrix refuses a matrix without enough rows');

my $index = File::Spec->catfile($tmp, 'genes2rows.txt');
write_file($index, "geneA\tx\t1\ngeneB\tx\t2\ngeneC\tx\t1\n");
my ($idx, $count) = quietly(sub { readGeneIdx($index) });
is($count, 3, 'readGeneIdx counts every index entry');
is_deeply($idx, { 1 => ['geneA', 'geneC'], 2 => ['geneB'] },
	'readGeneIdx groups gene names by catalogue row');

# --- MGS membership readers ------------------------------------------------
my $mgs = File::Spec->catfile($tmp, 'clusters.txt');
write_file($mgs, "MGS.1\tg1\nMGS.1\tg2\nMGS.2\tg3\n");
is_deeply(readMGS($mgs), { 'MGS.1' => ['g1', 'g2'], 'MGS.2' => ['g3'] },
	'readMGS groups genes by MGS in file order');
is_deeply(quietly(sub { readMGSrev($mgs) }), { g1 => 'MGS.1', g2 => 'MGS.1', g3 => 'MGS.2' },
	'readMGSrev maps each gene to its MGS');

# --- scheduler resource parsing -------------------------------------------
is(convert2Gb('0'), 0, 'zero scratch space stays zero');
is(convert2Gb('20G'), 20, 'gigabyte requests are kept');
is(convert2Gb('1.4G'), 1, 'fractional gigabytes are rounded to the nearest integer');
is(convert2Gb('2T'), 2048, 'terabyte requests are converted to gigabytes');
is(convert2Gb('512M'), 0.5, 'megabyte requests are converted to gigabytes');
is(convert2Gb('30'), 30, 'plain numbers are interpreted as gigabytes');
is(greaterComputeSpace('10G', '512M', '1T', '3'), 1024,
	'greaterComputeSpace returns the largest request in gigabytes');
is(greaterComputeSpace('512M'), 0.5, 'greaterComputeSpace converts megabytes');
ok(!eval { greaterComputeSpace('10GB'); 1 }, 'greaterComputeSpace rejects unknown units');

done_testing();
