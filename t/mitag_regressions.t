use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IPC::Open3;
use Symbol qw(gensym);
use Test::More;

my $root = File::Spec->catdir($Bin, '..');
my $scripts = File::Spec->catdir($root, 'secScripts', 'miTag');
my $tmp = tempdir(CLEANUP => 1);

sub write_file {
	my ($path, $contents) = @_;
	open my $fh, '>', $path or die "Cannot write $path: $!";
	print {$fh} $contents;
	close $fh;
}

sub read_gzip {
	my ($path) = @_;
	my $contents = '';
	gunzip($path => \$contents) or die "Cannot read $path: $GunzipError";
	return $contents;
}

sub run_script {
	my ($script, @arguments) = @_;
	my $error = gensym;
	my $pid = open3(undef, my $stdout, $error, $^X, '-I'.$root,
		File::Spec->catfile($scripts, $script), @arguments);
	my $output = do { local $/; <$stdout> } // '';
	my $errors = do { local $/; <$error> } // '';
	waitpid($pid, 0);
	return ($? >> 8, $output, $errors);
}

my $hierarchyDir = File::Spec->catdir($tmp, 'hierarchies');
mkdir $hierarchyDir or die $!;
my $hierarchy = join("\n",
	join("\t", qw(read domain phylum class order family genus)),
	join("\t", qw(r1 Bacteria Firmicutes Bacilli Bacillales Bacillaceae Bacillus)),
	join("\t", qw(r2 Bacteria)),
)."\n";
my $hierarchyPath = File::Spec->catfile($hierarchyDir, 'sample.hiera.txt.gz');
gzip(\$hierarchy => $hierarchyPath) or die "Cannot create fixture: $GzipError";

my $outputPrefix = File::Spec->catfile($tmp, 'merged');
my $unrelated = "$outputPrefix.keep";
write_file($unrelated, "preserve\n");
my ($status, $output, $errors) = run_script('miTagTaxTable.pl', 'family', $outputPrefix, $hierarchyDir);
is($status, 0, 'tax table supports a non-prefix subset of requested ranks');
my $familyOutput = read_gzip("$outputPrefix.family.txt.gz");
like($familyOutput, qr/^family\tsample$/m, 'gzip suffix is removed from the sample column name');
like($familyOutput, qr/^Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae\t1$/m,
	'family lineage is built from the header-derived family column');
like($familyOutput, qr/^Bacteria;\?;\?;\?;\?\t1$/m,
	'missing ranks are padded as separate semicolon-delimited fields');
ok(-e $unrelated, 'tax table cleanup preserves unrelated files sharing the output prefix');

my $readsDir = File::Spec->catdir($tmp, 'reads');
mkdir $readsDir or die $!;
my $r1 = File::Spec->catfile($readsDir, 'r1.fq');
my $r2a = File::Spec->catfile($readsDir, 'r2a.fq');
my $r2b = File::Spec->catfile($readsDir, 'r2b.fq');
write_file($_, "") for ($r1, $r2a, $r2b);
($status, $output, $errors) = run_script('catchLSUSSU.pl',
	'-R1', $r1, '-R2', "$r2a,$r2b", '-alignDir', File::Spec->catdir($tmp, 'align'),
	'-tmpDir', File::Spec->catdir($tmp, 'scratch'), '-smplID', 'sample', '-cores', 1);
isnt($status, 0, 'mismatched paired-read lists are rejected before tool execution');
like($errors, qr/same number of files/, 'paired-read mismatch has a clear diagnostic');

($status, $output, $errors) = run_script('lotus_LCA_blast3.pl',
	'-dir', $readsDir, '-DBdir', $readsDir, '-smplID', 'sample', '-lengthTolerance', 1.2);
isnt($status, 0, 'out-of-range LCA length tolerance is rejected');
like($errors, qr/lengthTolerance must be between 0 and 1/, 'length-tolerance validation is explicit');

done_testing();
