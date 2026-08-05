use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use IO::Compress::Gzip qw(gzip $GzipError);
use lib File::Spec->catdir($Bin, '..');
use Mods::GenoMetaAss qw(fileGZe fileGZs resolveExistingFile);
use Mods::IO_Tamoc_progs qw(
	buildMapperIdx checkMapsDoneSH inputFmtMegahit inputFmtMegahitRuntimeLibraries truePath
);
use Mods::ReadLibrary qw(newReadLibrary);
use Mods::Subm qw(qsubSystem2);

my $root = tempdir(CLEANUP => 1);
$root =~ s{\\}{/}g;

my $internal_config_path = File::Spec->catfile($Bin, '..', 'Mods', 'config_internal.txt');
open my $internal_config_fh, '<', $internal_config_path or die $!;
my $internal_config = do { local $/; <$internal_config_fh> };
close $internal_config_fh;
my $iqtree_selector = 'iqtree' . "\t" .
	'$(command -v iqtree3 || command -v iqtree2)' . "\t" . 'env:MF4phylo';
like(
	$internal_config,
	qr/^\Q$iqtree_selector\E$/m,
	'IQ-TREE configuration prefers iqtree3 and falls back to iqtree2',
);
my $r_environment_path = File::Spec->catfile($Bin, '..', 'helpers', 'install', 'MGTK_R.yml');
open my $r_environment_fh, '<', $r_environment_path or die $!;
my $r_environment = do { local $/; <$r_environment_fh> };
close $r_environment_fh;
like(
	$r_environment,
	qr/bioconda::bioconductor-ggtree=3\.14/,
	'the active MF4 R environment installs the R 4.4-compatible ggtree release',
);
my $installer_path = File::Spec->catfile($Bin, '..', 'helpers', 'install', 'installer.sh');
open my $installer_fh, '<', $installer_path or die $!;
my $installer = do { local $/; <$installer_fh> };
close $installer_fh;
like(
	$installer,
	qr/run -n MF4_R Rscript --vanilla -e.*?library\(ggtree\)/s,
	'the installer verifies that ggtree loads in MF4_R',
);
my $phylo_tools_path = File::Spec->catfile($Bin, '..', 'Mods', 'phyloTools.pm');
open my $phylo_tools_fh, '<', $phylo_tools_path or die $!;
my $phylo_tools = do { local $/; <$phylo_tools_fh> };
close $phylo_tools_fh;
like(
	$phylo_tools,
	qr/my \$cmd = "\$iqTree -s \$inMSA \$threadOpts -pre \$treeOut -seed 678 -quiet "/,
	'IQ-TREE invocations use the selected bounded thread policy and quiet output',
);
unlike(
	$phylo_tools,
	qr/(?:^|\s)-keep-ident(?:\s|$)/m,
	'IQ-TREE uses its standard identical-sequence handling without -keep-ident',
);
like(
	$phylo_tools,
	qr/my \$threadOpts = "-T \$ncore"/,
	'IQ-TREE uses every core allocated to its cluster job without AUTO benchmarking',
);
unlike(
	$phylo_tools,
	qr/-T AUTO|threads-max|ntmax/,
	'IQ-TREE commands no longer repeat automatic thread-efficiency tests',
);
like(
	$phylo_tools,
	qr/my \$usePartitionModel = \$partiF ne "" && !\$iqPathogen.*?if \(!\$iqLegacy && \$iqMemMB > 0\).*?if \(\$usePartitionModel\).*?WARNING: IQ-TREE -mem disabled because partition models do not support.*?else \{\s*\$cmd \.= "-mem \$\{iqMemMB\}M "/s,
	'modern IQ-TREE uses the RAM cap only when it is compatible with the active model',
);
like(
	$phylo_tools,
	qr/\$cmd \.= " -p \$partiF " if \$usePartitionModel.*?-m MFP\+MERGE/s,
	'partitioned IQ-TREE ModelFinder receives the generated loci and optimizes their merge scheme',
);
like(
	$phylo_tools,
	qr/\$cmd \.= "--pathogen " if \$iqPathogen && !\$iqLegacy/,
	'IQ-TREE 3 pathogen mode is available to low-divergence callers',
);
like(
	$phylo_tools,
	qr/my \$cmapleLengthLimit = 32767.*?_fastaAlignmentLength\(\$inMSA\).*?WARNING: IQ-TREE --pathogen disabled.*?falling back to standard IQ-TREE mode/s,
	'overlong alignments warn and fall back from CMAPLE pathogen mode before IQ-TREE runs',
);
like(
	$phylo_tools,
	qr/\$cmd \.= \$iqLegacy \? "-m GTR\+F\+I\+G4 " : "-m GTR\+F\+G2 "/,
	'modern nucleotide trees use GTR+F+G2 while legacy mode retains the previous model',
);
like(
	$phylo_tools,
	qr/\$taxonCount >= 750.*?\$runSafe = 1.*?_iqtreeLogRequestsSafeKernel.*?restarting once with -safe/s,
	'large strain trees start with the safe kernel and explicit underflow triggers one safe retry',
);
like(
	$phylo_tools,
	qr/iqtreeOutputComplete\(\$treeOut, \$inMSA.*?cleanupIQTreeTransients\(\$treeOut\)/s,
	'IQ-TREE completion is taxon-validated before temporary artifacts are removed',
);

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
like($runtime_setup, qr/mf4_fastq_has_records/,
	'runtime MEGAHIT setup tests for actual FASTQ records rather than compressed byte size');
like($runtime_setup, qr/Skipping zero-record cleaned paired library for MEGAHIT/,
	'runtime MEGAHIT setup skips a cleaner output that contains no paired records');
like($runtime_setup, qr/No FASTQ records remain for MEGAHIT after cleaning/,
	'runtime MEGAHIT setup stops an all-empty assembly group before invoking MEGAHIT');
my $runtime_script = "$root/megahit-inputs.sh";
sub write_gzip {
	my ($path, $contents) = @_;
	gzip \$contents => $path or die "Cannot gzip $path: $GzipError";
}
write_gzip("$root/lib 0.1.fq.gz", '');
write_gzip("$root/lib 0.2.fq.gz", '');
my $fastq = "\@read\nACGT\n+\n!!!!\n";
for my $read ("$root/lib1.1.fq.gz", "$root/lib1.2.fq.gz", "$root/lib1.singl.fq.gz") {
	write_gzip($read, $fastq);
}
open my $runtime_fh, '>', $runtime_script or die $!;
print {$runtime_fh} "#!/bin/bash\nset -eo pipefail\n$runtime_setup\nprintf '<%s>\\n' $runtime_args\n";
close $runtime_fh;
is(system('bash', '-n', $runtime_script), 0,
	'generated runtime MEGAHIT input setup is valid Bash');
my $runtime_output = `bash "$runtime_script" 2>/dev/null`;
is($?, 0, 'runtime MEGAHIT input setup executes successfully with another nonempty library');
like($runtime_output, qr/<-1>\n<\Q${root}\/lib1.1.fq.gz\E>\n<-2>\n<\Q${root}\/lib1.2.fq.gz\E>/,
	'the generated MEGAHIT paired inputs exclude the header-only gzip library');
like($runtime_output, qr/<-r>\n<\Q${root}\/lib1.singl.fq.gz\E>/,
	'the generated MEGAHIT -r list retains an existing singleton with FASTQ records');
unlike($runtime_output, qr/lib 0\.(?:1|2)\.fq\.gz/,
	'the generated MEGAHIT command omits zero-record paired files');
unlink "$root/lib1.1.fq.gz" or die $!;
unlink "$root/lib1.2.fq.gz" or die $!;
unlink "$root/lib1.singl.fq.gz" or die $!;
write_gzip("$root/lib1.1.fq.gz", '');
write_gzip("$root/lib1.2.fq.gz", '');
write_gzip("$root/lib1.singl.fq.gz", '');
is(system('bash', $runtime_script) >> 8, 42,
	'an all-zero-record assembly group exits before MEGAHIT can create an empty FASTA');

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
