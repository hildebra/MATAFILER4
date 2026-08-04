use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use IPC::Open3;
use Symbol qw(gensym);
use Test::More;

my $root = File::Spec->catdir($Bin, '..');
my $script = File::Spec->catfile(
	$root, 'secScripts', 'miTag', 'catchLSUSSU.pl',
);
my $tmp = tempdir(CLEANUP => 1);

sub write_file {
	my ($path, $contents) = @_;
	open my $fh, '>', $path or die "Cannot write $path: $!";
	print {$fh} $contents;
	close $fh or die "Cannot close $path: $!";
}

sub run_with_ssu_reference {
	my ($tag, $ssu_reference) = @_;
	my $case_root = File::Spec->catdir($tmp, $tag);
	make_path($case_root);
	my $reads = File::Spec->catfile($case_root, 'reads.fq');
	my $lsu = File::Spec->catfile($case_root, 'lsu.fasta');
	my $config = File::Spec->catfile($case_root, 'MATAFILERcfg.txt');
	write_file($reads, "read fixture\n");
	write_file($lsu, ">lsu\nACGT\n");
	write_file($config, join("\n",
		"MFLRDir\t$root",
		"BINDir\t$case_root",
		"Rpath\tR",
		"DBDir\t$case_root",
		"CONDcmd\tconda",
		"CONDA\tshell hook",
		"Rscript\tRscript",
		"sortmerna\t/bin/true",
		"SSUdbFAsrt\t$ssu_reference",
		"LSUdbFAsrt\t$lsu",
		"SSUidx\t$case_root",
		"LSUidx\t$case_root",
	)."\n");

	my $stderr = gensym;
	my $pid = open3(
		undef, my $stdout, $stderr,
		$^X, '-I'.$root, $script,
		'-R1', $reads, '-R2', $reads,
		'-alignDir', File::Spec->catdir($case_root, 'align'),
		'-tmpDir', File::Spec->catdir($case_root, 'scratch'),
		'-smplID', 'sample', '-cores', 1,
		'-config', $config,
	);
	my $output = do { local $/; <$stdout> } // '';
	my $errors = do { local $/; <$stderr> } // '';
	waitpid($pid, 0);
	return ($? >> 8, $output, $errors);
}

is(system($^X, '-I'.$root, '-c', $script), 0,
	'catchLSUSSU compiles');

my $missing = File::Spec->catfile($tmp, 'missing-ssu.fasta');
my ($status, undef, $errors) = run_with_ssu_reference('missing', $missing);
isnt($status, 0, 'a missing SortMeRNA reference is rejected');
like($errors, qr/does not exist: \Q$missing\E/,
	'the missing-reference diagnostic is distinct from access failure');

my $directory = File::Spec->catdir($tmp, 'ssu-directory');
make_path($directory);
($status, undef, $errors) = run_with_ssu_reference('directory', $directory);
isnt($status, 0, 'a directory cannot be used as a SortMeRNA reference');
like($errors, qr/is not a regular file: \Q$directory\E/,
	'the non-file diagnostic identifies the configured path');

my $empty = File::Spec->catfile($tmp, 'empty-ssu.fasta');
write_file($empty, '');
($status, undef, $errors) = run_with_ssu_reference('empty', $empty);
isnt($status, 0, 'an empty SortMeRNA reference is rejected');
like($errors, qr/is empty: \Q$empty\E/,
	'the empty-reference diagnostic is explicit');

my $unreadable = File::Spec->catfile($tmp, 'unreadable-ssu.fasta');
write_file($unreadable, ">ssu\nACGT\n");
chmod 0000, $unreadable or die "Cannot make $unreadable unreadable: $!";
SKIP: {
	skip 'filesystem or test identity does not enforce mode-000 read denial', 3
		if -r $unreadable;
	($status, undef, $errors) = run_with_ssu_reference('unreadable', $unreadable);
	isnt($status, 0, 'an existing inaccessible SortMeRNA reference is rejected');
	like($errors, qr/exists but cannot be opened for reading: \Q$unreadable\E: Permission denied/,
		'the access diagnostic preserves the operating-system error');
	like($errors, qr/effective uid=.*?effective groups=/,
		'the access diagnostic reports the Slurm-process credential context');
}
chmod 0600, $unreadable or die "Cannot restore $unreadable permissions: $!";

done_testing;
