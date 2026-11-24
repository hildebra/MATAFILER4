use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use IO::Uncompress::Gunzip qw($GunzipError);
use Test::More;

use lib "$Bin/..", "$Bin/lib";
use MFTestConfig;
use Mods::IO_Tamoc_progs qw(bgzfArgs bgzfCmd bgzfFilesCmd bgzfTreeCmd);
use Mods::GenoMetaAss qw(gzipopen gzipwrite);
use Mods::WorkflowControl ();

# Pipeline output is compressed as BGZF through the configured bgzip (t/bin/bgzip
# in tests, which writes real BGZF). These checks run the generated job-script
# commands and read the results back the ways the pipeline does.

my $root = tempdir(CLEANUP => 1);

sub write_text {
	my ($path, $text) = @_;
	open my $fh, '>', $path or die "Cannot write $path: $!\n";
	print {$fh} $text;
	close $fh or die "Cannot close $path: $!\n";
}

sub slurp_raw {
	my ($path) = @_;
	open my $fh, '<:raw', $path or die "Cannot read $path: $!\n";
	local $/;
	my $data = <$fh>;
	close $fh;
	return $data;
}

# Number of BGZF blocks, or an error: every member must carry the BC extra field
# with its exact size, and the file must end with the BGZF end-of-file block.
sub bgzf_blocks {
	my ($path) = @_;
	my $data = slurp_raw($path);
	my ($pos, $blocks) = (0, 0);
	while ($pos < length $data) {
		return "no BGZF header at byte $pos"
			unless substr($data, $pos, 4) eq "\x1f\x8b\x08\x04" && substr($data, $pos + 12, 2) eq 'BC';
		$pos += unpack('v', substr($data, $pos + 16, 2)) + 1;
		$blocks++;
	}
	return "block sizes overrun the file" if $pos != length $data;
	return "missing BGZF end-of-file block"
		unless substr($data, -28) eq pack('H*', '1f8b08040000000000ff0600424302001b0003000000000000000000');
	return $blocks;
}

sub gunzip_text {
	my ($path) = @_;
	my ($fh, $ok) = gzipopen($path, 'test input', 1, 0);
	my $text = do { local $/; <$fh> };
	close $fh or die "Reading $path failed\n";
	return $text;
}

sub run_shell {
	my ($script) = @_;
	my $file = File::Spec->catfile($root, 'job.sh');
	write_text($file, $script);
	system('bash', $file);
	return $? >> 8;
}

# Options: threads and level as bgzip expects them, from the configured bgzip.
my @args = bgzfArgs(4, level => 6);
like($args[0], qr{t/bin/bgzip$}, 'bgzip comes from the bgzip config key');
is_deeply([@args[1 .. $#args]], ['-@', 4, '-l', 6], 'threads and compression level are bgzip options');
is_deeply([(bgzfArgs(undef))[1, 2]], ['-@', 1], 'a missing thread count means one thread');
ok(!eval { bgzfArgs(2, level => 'high'); 1 }, 'an invalid compression level is refused');

{
	local $Mods::IO_Tamoc_progs::CONFIG_HASH{bgzip} = 'bgzip';
	like((bgzfArgs(1))[0], qr{^/.*/bgzip$},
		'a bare bgzip name is resolved on PATH to an absolute path for job scripts');
	local $Mods::IO_Tamoc_progs::CONFIG_HASH{bgzip} = 'bgzip-not-installed-anywhere';
	is((bgzfArgs(1))[0], 'bgzip-not-installed-anywhere', 'an unresolvable bgzip name is kept as configured');
}

# Stream compression, larger than one BGZF block.
my $big = join('', map { "\@read$_\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n" } 1 .. 5000);
my $plain = File::Spec->catfile($root, 'reads.fq');
write_text($plain, $big);
my $stream = File::Spec->catfile($root, 'stream.fq.gz');
is(run_shell("cat $plain | " . bgzfCmd(2) . " -c > $stream\n"), 0, 'stream compression succeeds');
my $blocks = bgzf_blocks($stream);
ok($blocks =~ /^\d+$/ && $blocks > 2, "stream output is BGZF with several blocks ($blocks)");
is(gunzip_text($stream), $big, 'gzipopen reads every BGZF block back');
my $gunzip = IO::Uncompress::Gunzip->new($stream, MultiStream => 1) or die $GunzipError;
is(do { local $/; <$gunzip> }, $big, 'IO::Uncompress::Gunzip with MultiStream reads every block');
is(Mods::WorkflowControl::_fastq_file_has_records($stream), 1, 'FASTQ record check reads BGZF input');

# Empty output, as written for steps without hits.
my $empty = File::Spec->catfile($root, 'empty.gz');
is(run_shell(bgzfCmd(1) . " -c </dev/null > $empty\n"), 0, 'empty BGZF output is written');
is(bgzf_blocks($empty), 1, 'empty BGZF output holds only the end-of-file block');
is(gunzip_text($empty), '', 'empty BGZF output reads back empty');

# Perl-side writer.
my $table = File::Spec->catfile($root, 'table.tsv.gz');
my $writer = gzipwrite($table, 'test table', { threads => 2 });
print {$writer} $big;
ok(close($writer), 'gzipwrite closes successfully');
ok(bgzf_blocks($table) =~ /^\d+$/, 'gzipwrite output is BGZF');
is(gunzip_text($table), $big, 'gzipwrite output reads back unchanged');

# In place: explicit files and wildcards; existing .gz files and directories are left alone.
my $dir = File::Spec->catdir($root, 'inplace');
make_path("$dir/sub");
write_text("$dir/a.txt", "a\n");
write_text("$dir/b.txt", "b\n");
write_text("$dir/sub/c.txt", "c\n");
write_text("$dir/single.fna", ">x\nACGT\n");
my $kept = File::Spec->catfile($dir, 'old.gz');
is(run_shell("printf 'old\\n' | " . bgzfCmd(1) . " -c > $kept\n"), 0, 'fixture .gz file written');
my $kept_before = slurp_raw($kept);
is(run_shell(bgzfFilesCmd(2, "$dir/single.fna")), 0, 'in-place compression of one file succeeds');
ok(!-e "$dir/single.fna" && -s "$dir/single.fna.gz", 'the file is replaced by its .gz');
is(gunzip_text("$dir/single.fna.gz"), ">x\nACGT\n", 'in-place output reads back unchanged');
is(run_shell(bgzfFilesCmd(2, "$dir/*")), 0, 'in-place compression of a wildcard succeeds');
ok(-s "$dir/a.txt.gz" && -s "$dir/b.txt.gz" && !-e "$dir/a.txt" && !-e "$dir/b.txt",
	'every matching regular file is compressed');
is(slurp_raw($kept), $kept_before, 'a file already ending in .gz is left unchanged');
ok(!-e "$kept.gz", 'a .gz file is not compressed a second time');
ok(-d "$dir/sub" && -e "$dir/sub/c.txt", 'directories are left alone by the wildcard');
is(bgzf_blocks("$dir/a.txt.gz"), 2, 'wildcard output is BGZF');

# Recursive, replacing pigz -r.
is(run_shell(bgzfTreeCmd(2, "$dir/sub", $kept)), 0, 'recursive compression succeeds');
ok(-s "$dir/sub/c.txt.gz" && !-e "$dir/sub/c.txt", 'files below a directory are compressed');
is(slurp_raw($kept), $kept_before, 'recursive compression leaves .gz files alone');
is(bgzfTreeCmd(2), '', 'no paths, no command');

done_testing();
