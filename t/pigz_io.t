use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::GenoMetaAss qw(gzipopen gzipwrite);

my $root = tempdir(CLEANUP => 1);
my $fake_pigz = File::Spec->catfile($Bin, 'bin', 'pigz');
my $log = File::Spec->catfile($root, 'pigz.log');
my $fake_bgzip = File::Spec->catfile($Bin, 'bin', 'bgzip');
my $bgzip_log = File::Spec->catfile($root, 'bgzip.log');

local $ENV{MF4_TEST_PIGZ_LOG} = $log;
local $ENV{MF4_TEST_BGZIP_LOG} = $bgzip_log;
local $Mods::IO_Tamoc_progs::CONFIG_LOADED = 1;
local %Mods::IO_Tamoc_progs::CONFIG_HASH = (pigz => $fake_pigz, bgzip => $fake_bgzip);

my $first = File::Spec->catfile($root, 'first.txt.gz');
my $first_writer = gzipwrite(
	$first, 'first test member', { level => 3, threads => 2 },
);
print {$first_writer} "header\tS1\tS2\nrow1\t1\t2\n"
	or die "Cannot write first test member: $!\n";
ok(close($first_writer), 'bgzip-backed writer closes successfully');
ok(-s $first, 'bgzip-backed writer creates a nonempty gzip file');

my $second = File::Spec->catfile($root, 'second.txt.gz');
my $second_writer = gzipwrite($second, 'second test member');
print {$second_writer} "row2\t3\t4\n"
	or die "Cannot write second test member: $!\n";
ok(close($second_writer), 'a second bgzip-backed writer closes successfully');

open my $second_input, '<', $second or die "Cannot open $second: $!\n";
open my $first_append, '>>', $first or die "Cannot append to $first: $!\n";
binmode $second_input;
binmode $first_append;
my $buffer;
while (read($second_input, $buffer, 64 * 1024)) {
	print {$first_append} $buffer or die "Cannot append gzip member: $!\n";
}
close $second_input or die "Cannot close $second: $!\n";
close $first_append or die "Cannot close $first: $!\n";

my ($reader, $ok) = gzipopen($first, 'concatenated test input', 1, 0);
ok($ok && defined($reader), 'gzipopen launches configured pigz');
my $contents = do { local $/; <$reader> };
ok(close($reader), 'pigz-backed reader validates the complete stream on close');
is(
	$contents,
	"header\tS1\tS2\nrow1\t1\t2\nrow2\t3\t4\n",
	'gzipopen reads every concatenated gzip member',
);

open my $log_fh, '<', $log or die "Cannot open $log: $!\n";
my $invocations = do { local $/; <$log_fh> };
close $log_fh or die "Cannot close $log: $!\n";
open my $bgzip_log_fh, '<', $bgzip_log or die "Cannot open $bgzip_log: $!\n";
my $bgzip_invocations = do { local $/; <$bgzip_log_fh> };
close $bgzip_log_fh or die "Cannot close $bgzip_log: $!\n";
like($bgzip_invocations, qr/^-\@\t2\t-l\t3\t-c$/m,
	'gzipwrite passes threads and level to the configured bgzip');
open my $second_raw, '<:raw', $second or die "Cannot open $second: $!\n";
read($second_raw, my $header, 18) or die "Cannot read $second: $!\n";
close $second_raw;
is(substr($header, 0, 4) . substr($header, 12, 2), "\x1f\x8b\x08\x04BC",
	'gzipwrite writes BGZF (gzip members carrying the BC extra field)');
like($invocations, qr/-dc\t--\t\Q$first\E/,
	'gzipopen passes the compressed path to configured pigz without a shell');

my $matrix = File::Spec->catfile($root, 'Matrix.mat.gz');
my $matrix_writer = gzipwrite($matrix, 'large matrix fixture');
print {$matrix_writer} "gene\tS1\tS2\n"
	or die "Cannot write matrix header: $!\n";
for my $row (1 .. 20_000) {
	print {$matrix_writer} "gene$row\t$row\t", $row + 1, "\n"
		or die "Cannot write matrix row: $!\n";
}
ok(close($matrix_writer), 'large matrix fixture is compressed successfully');

my $mgs_path = File::Spec->catfile($Bin, '..', 'secScripts', 'MGS.pl');
open my $mgs_fh, '<', $mgs_path or die "Cannot open $mgs_path: $!\n";
my $mgs_source = do { local $/; <$mgs_fh> };
close $mgs_fh or die "Cannot close $mgs_path: $!\n";
my ($matrix_helper) =
	$mgs_source =~ /(sub _matrix_sample_count \{.*?^\})/ms;
ok(defined($matrix_helper), 'MGS matrix-header helper can be isolated');
my $helper_loaded = eval "$matrix_helper\n1;";
ok($helper_loaded, 'MGS matrix-header helper compiles independently') or diag($@);
is(main::_matrix_sample_count($matrix), 2,
	'MGS counts samples without failing when a large pigz stream is closed after its header');

done_testing();
