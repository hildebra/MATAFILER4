use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::GenoMetaAss qw(fileGZe fileGZs resolveExistingFile);
use Mods::IO_Tamoc_progs qw(
	buildMapperIdx checkMapsDoneSH inputFmtMegahit inputFmtMegahitRuntimeLibraries truePath
);
use Mods::ReadLibrary qw(newReadLibrary);
use Mods::Subm qw(qsubSystem2);

my $root = tempdir(CLEANUP => 1);
$root =~ s{\\}{/}g;

my @p1 = ('a.1.fq', 'b.1.fq', 'c.1.fq');
my @p2 = ('a.2.fq', 'b.2.fq', 'c.2.fq');
my @singletons = ('single.a.fq', 'single.b.fq');
is(
	inputFmtMegahit(\@p1, \@p2, \@singletons, $root),
	'-1 a.1.fq,b.1.fq,c.1.fq -2 a.2.fq,b.2.fq,c.2.fq -r single.a.fq,single.b.fq',
	'MEGAHIT accepts independent paired and singleton library counts',
);
my $error = '';
eval { inputFmtMegahit(['a.1.fq'], [], [], $root) };
$error = $@;
like($error, qr/Unequal paired read array lengths for MEGAHIT/,
	'MEGAHIT rejects genuinely unpaired mate arrays');

my $runtime_libraries = [
	newReadLibrary(
		id => 'lib0', scope => 'primary', phase => 'clean',
		files => {r1 => "$root/lib 0.1.fq.gz", r2 => "$root/lib 0.2.fq.gz",
			single => "$root/lib 0.singl.fq.gz"},
	),
	newReadLibrary(
		id => 'lib1', scope => 'primary', phase => 'clean',
		files => {r1 => "$root/lib1.1.fq.gz", r2 => "$root/lib1.2.fq.gz",
			single => "$root/lib1.singl.fq.gz"},
	),
];
my ($runtime_setup, $runtime_args) = inputFmtMegahitRuntimeLibraries(
	$runtime_libraries, 'mh_args',
);
is($runtime_args, '"${mh_args[@]}"',
	'runtime MEGAHIT arguments use a Bash array without word splitting');
like($runtime_setup, qr/Skipping missing or empty optional MEGAHIT singleton/,
	'runtime MEGAHIT setup treats projected singleton files as optional');
like($runtime_setup, qr/Missing or empty paired input for MEGAHIT/,
	'runtime MEGAHIT setup still requires both paired-read files');
my $runtime_script = "$root/megahit-inputs.sh";
for my $read ("$root/lib 0.1.fq.gz", "$root/lib 0.2.fq.gz",
	"$root/lib1.1.fq.gz", "$root/lib1.2.fq.gz", "$root/lib1.singl.fq.gz") {
	open my $read_fh, '>', $read or die $!;
	print {$read_fh} "non-empty\n";
	close $read_fh;
}
open my $runtime_fh, '>', $runtime_script or die $!;
print {$runtime_fh} "#!/bin/bash\n$runtime_setup\nprintf '<%s>\\n' $runtime_args\n";
close $runtime_fh;
is(system('bash', '-n', $runtime_script), 0,
	'generated runtime MEGAHIT input setup is valid Bash');
my $runtime_output = `bash "$runtime_script" 2>/dev/null`;
is($?, 0, 'runtime MEGAHIT input setup executes successfully');
like($runtime_output, qr/<-r>\n<\Q${root}\/lib1.singl.fq.gz\E>/,
	'the generated MEGAHIT -r list retains an existing singleton');
unlike($runtime_output, qr/lib 0\.singl/,
	'the generated MEGAHIT -r list omits a missing projected singleton');

{
	no warnings 'redefine';
	local *Mods::IO_Tamoc_progs::getProgPaths = sub { return 'bwa' };
	my (undef, undef, $check_path) = buildMapperIdx("$root/reference.fa", 2, 0, 2);
	is($check_path, "$root/reference.fa.pac", 'BWA index check uses the actual .pac path');
}

{
	my ($command, $index_path, $check_path) = buildMapperIdx("$root/reference.fa", 2, 0, -2);
	is($command, '', 'automatic strobealign selection emits no reference-index command');
	is($index_path, "$root/reference.fa", 'strobealign maps directly against the reference FASTA');
	is($check_path, "$root/reference.fa", 'strobealign reference check uses the FASTA itself');
}

my $empty = "$root/empty.txt";
open my $empty_fh, '>', $empty or die $!;
close $empty_fh;
ok(!fileGZe($empty), 'zero-byte artifacts are not complete');
my $plain = "$root/plain.txt";
open my $plain_fh, '>', $plain or die $!;
print {$plain_fh} '12345';
close $plain_fh;
ok(fileGZe("$plain.gz"), 'gzip-aware existence check resolves a plain alternative');
is(fileGZs("$plain.gz"), 5, 'gzip-aware size check resolves a plain alternative');
my $gzip = "$root/compressed.txt.gz";
open my $gzip_fh, '>', $gzip or die $!;
print {$gzip_fh} '1234';
close $gzip_fh;
is(fileGZs("$root/compressed.txt"), 20,
	'fileGZs retains the legacy five-times compressed-size estimate');
my ($resolved_plain, $resolved_stat) = resolveExistingFile("$plain.gz");
is($resolved_plain, $plain, 'plain/gzip resolution returns the existing alternative');
is($resolved_stat->[7], 5, 'plain/gzip resolution reuses the selected file stat');

my $mapping_check = checkMapsDoneSH(["$root/sample/"]);
like($mapping_check, qr/find .*?-smd\.bam.*?-smd\.cram.*?-size \+0c/,
	'directory mapping check requires a non-empty BAM or CRAM as well as its marker');
my $file_check = checkMapsDoneSH(["$root/map.cram"]);
like($file_check, qr/\[ ! -s \Q$root\/map.cram\E \]/,
	'direct mapping input must be non-empty');
my $mapping_script = "$root/check-mapping.sh";
open my $mapping_fh, '>', $mapping_script or die $!;
print {$mapping_fh} "#!/bin/bash\nset -eo pipefail\n$mapping_check$file_check";
close $mapping_fh;
is(system('bash', '-n', $mapping_script), 0,
	'generated mapping prerequisite checks are valid Bash');

{
	local $ENV{MF4_TEST_ROOT} = $root;
	is(truePath('$MF4_TEST_ROOT/results'), "$root/results",
		'truePath expands environment variables with path suffixes');
}
{
	local $ENV{MF4_TEST_MISSING};
	delete $ENV{MF4_TEST_MISSING};
	my $path_error = '';
	eval { truePath('$MF4_TEST_MISSING/results') };
	$path_error = $@;
	like($path_error, qr/Environment variable \$MF4_TEST_MISSING .* is not set/,
		'truePath reports unset environment variables');
}

my $script = "$root/direct.sh";
my $ran = "$root/direct.ran";
open my $script_fh, '>', $script or die $!;
print {$script_fh} "#!/bin/bash\nset -e\ntouch $ran\n";
close $script_fh;
qsubSystem2($script, { qmode => 'bash' });
ok(-e $ran, 'qsubSystem2 no longer aborts at its old debug die');

done_testing;
