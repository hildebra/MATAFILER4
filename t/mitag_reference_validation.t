use strict;
use warnings;

use File::Basename qw(dirname);
use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IPC::Open3;
use Symbol qw(gensym);
use Test::More;

my $root = File::Spec->rel2abs(File::Spec->catdir($Bin, '..'));
my $scripts = File::Spec->catdir($root, 'secScripts', 'miTag');
my $tmp = tempdir(CLEANUP => 1);

sub write_file {
	my ($path, $contents) = @_;
	make_path(dirname($path));
	open my $handle, '>', $path or die "Cannot write $path: $!";
	print {$handle} $contents;
	close $handle or die "Cannot close $path: $!";
}

sub write_gzip {
	my ($path, $contents) = @_;
	make_path(dirname($path));
	gzip(\$contents => $path)
		or die "Cannot create $path: $GzipError";
}

sub read_gzip {
	my ($path) = @_;
	my $contents = '';
	gunzip($path => \$contents, MultiStream => 1)
		or die "Cannot read $path: $GunzipError";
	return $contents;
}

sub run_script {
	my ($script, @arguments) = @_;
	my $error = gensym;
	my $pid = open3(
		undef, my $stdout, $error,
		$^X, '-I'.$root, File::Spec->catfile($scripts, $script),
		@arguments,
	);
	my $output = do { local $/; <$stdout> } // '';
	my $errors = do { local $/; <$error> } // '';
	waitpid($pid, 0);
	return ($? >> 8, $output, $errors);
}

sub log_lines {
	my ($path) = @_;
	return () unless -e $path;
	open my $handle, '<', $path or die "Cannot read $path: $!";
	my @lines = <$handle>;
	close $handle or die "Cannot close $path: $!";
	chomp @lines;
	return @lines;
}

my $fakeSortmerna = File::Spec->catfile($tmp, 'fake-sortmerna');
write_file($fakeSortmerna, q{#!/usr/bin/env perl
use strict;
use warnings;
use File::Path qw(make_path);
use IO::Compress::Gzip qw(gzip $GzipError);

if (grep { $_ eq '--version' } @ARGV) {
	print "SortMeRNA test double 4.3.7\n";
	exit 0;
}
my ($aligned, $workdir, $zipOut) = ('', '', '');
my (@reads, @flags);
my $paired = 0;
for (my $index = 0; $index < @ARGV; $index++) {
	my $argument = $ARGV[$index];
	push @flags, $argument;
	if ($argument eq '--reads') {
		push @reads, $ARGV[++$index];
	} elsif ($argument eq '--aligned') {
		$aligned = $ARGV[++$index];
	} elsif ($argument eq '--workdir') {
		$workdir = $ARGV[++$index];
	} elsif ($argument eq '--zip-out') {
		$zipOut = $ARGV[++$index];
	} elsif ($argument eq '--out2') {
		$paired = 1;
	}
}
die "fake SortMeRNA needs --aligned\n" if $aligned eq '';
die "catchLSUSSU did not force gzip output\n" unless $zipOut eq '1';
make_path($workdir) if $workdir ne '' && !-d $workdir;
if (defined($ENV{FAKE_SORTMERNA_LOG}) && $ENV{FAKE_SORTMERNA_LOG} ne '') {
	open my $log, '>>', $ENV{FAKE_SORTMERNA_LOG} or die $!;
	print {$log} $aligned, "\t", join(' ', @flags), "\n";
	close $log or die $!;
}
if (defined($ENV{FAKE_SORTMERNA_MISSING_TAG})
		&& $aligned =~ /$ENV{FAKE_SORTMERNA_MISSING_TAG}/) {
	exit 0;
}
my $payload = $ENV{FAKE_SORTMERNA_ZERO}
	? ''
	: "\@fake\nACGT\n+\nIIII\n";
my @outputs = $paired
	? ($aligned.'_fwd.fq.gz', $aligned.'_rev.fq.gz')
	: ($aligned.'.fq.gz');
for my $output (@outputs) {
	gzip(\$payload => $output)
		or die "Cannot create $output: $GzipError\n";
}
exit 0;
});
chmod 0755, $fakeSortmerna
	or die "Cannot make $fakeSortmerna executable: $!";

my $fakeVsearch = File::Spec->catfile($tmp, 'fake-vsearch');
write_file($fakeVsearch, q{#!/usr/bin/env perl
use strict;
use warnings;
my ($indexOutput, $userOutput) = ('', '');
for (my $index = 0; $index < @ARGV; $index++) {
	if ($ARGV[$index] eq '--output') {
		$indexOutput = $ARGV[++$index];
	} elsif ($ARGV[$index] eq '--userout') {
		$userOutput = $ARGV[++$index];
	}
}
my $output = $indexOutput ne '' ? $indexOutput : $userOutput;
die "fake VSEARCH has no output argument\n" if $output eq '';
open my $handle, '>', $output or die $!;
print {$handle} "fake index\n" if $indexOutput ne '';
print {$handle} "fake\tref\t99\t4\t0\t0\t1\t4\t1\t4\t4\n"
	if $userOutput ne '';
close $handle or die $!;
if (defined($ENV{FAKE_VSEARCH_LOG}) && $ENV{FAKE_VSEARCH_LOG} ne '') {
	open my $log, '>>', $ENV{FAKE_VSEARCH_LOG} or die $!;
	print {$log} join(' ', @ARGV), "\n";
	close $log or die $!;
}
exit 0;
});
chmod 0755, $fakeVsearch
	or die "Cannot make $fakeVsearch executable: $!";

my $fakeLca = File::Spec->catfile($tmp, 'fake-lca');
write_file($fakeLca, q{#!/usr/bin/env perl
use strict;
use warnings;
my ($output, $inputs) = ('', '');
for (my $index = 0; $index < @ARGV; $index++) {
	if ($ARGV[$index] eq '-o') {
		$output = $ARGV[++$index];
	} elsif ($ARGV[$index] eq '-i') {
		$inputs = $ARGV[++$index];
	}
}
die "fake LCA has no output\n" if $output eq '';
open my $handle, '>', $output or die $!;
print {$handle} "read\tdomain\tphylum\nfake\tBacteria\tFirmicutes\n";
close $handle or die $!;
if (defined($ENV{FAKE_LCA_LOG}) && $ENV{FAKE_LCA_LOG} ne '') {
	open my $log, '>>', $ENV{FAKE_LCA_LOG} or die $!;
	print {$log} $inputs, "\n";
	close $log or die $!;
}
exit 0;
});
chmod 0755, $fakeLca or die "Cannot make $fakeLca executable: $!";

my $fakeFlash = File::Spec->catfile($tmp, 'fake-flash');
write_file($fakeFlash, q{#!/usr/bin/env perl
use strict;
use warnings;
use File::Path qw(make_path);
my ($tag, $directory) = ('', '');
for (my $index = 0; $index < @ARGV; $index++) {
	if ($ARGV[$index] eq '-o') {
		$tag = $ARGV[++$index];
	} elsif ($ARGV[$index] eq '-d') {
		$directory = $ARGV[++$index];
	}
}
die "fake FLASH needs -o and -d\n" if $tag eq '' || $directory eq '';
make_path($directory) unless -d $directory;
my %outputs = (
	'extendedFrags.fastq' => "\@merged\nACGTACGT\n+\nIIIIIIII\n",
	'notCombined_1.fastq' => "\@unmerged/1\nACGT\n+\nIIII\n",
	'notCombined_2.fastq' => "\@unmerged/2\nTGCA\n+\nIIII\n",
);
for my $suffix (keys %outputs) {
	my $path = $directory.'/'.$tag.'.'.$suffix;
	open my $handle, '>', $path or die $!;
	print {$handle} $outputs{$suffix};
	close $handle or die $!;
}
exit 0;
});
chmod 0755, $fakeFlash or die "Cannot make $fakeFlash executable: $!";

my $databaseDir = File::Spec->catdir($tmp, 'databases');
my $indexDir = File::Spec->catdir($tmp, 'sortmerna-index');
make_path($databaseDir, $indexDir);
my %database = (
	ssu_sort => File::Spec->catfile($databaseDir, 'sort-ssu.fasta'),
	lsu_sort => File::Spec->catfile($databaseDir, 'sort-lsu.fasta'),
	ssu => File::Spec->catfile($databaseDir, 'ssu.fasta'),
	ssu_tax => File::Spec->catfile($databaseDir, 'ssu.tax'),
	lsu => File::Spec->catfile($databaseDir, 'lsu.fasta'),
	lsu_tax => File::Spec->catfile($databaseDir, 'lsu.tax'),
);
write_file($_, '') for values %database;

my $config = File::Spec->catfile($tmp, 'matafiler.cfg');
write_file($config, join("\n",
	"MFLRDir\t$root",
	"BINDir\t$tmp",
	"Rpath\tR",
	"DBDir\t$databaseDir",
	"CONDcmd\tconda",
	"CONDA\tshell hook",
	"Rscript\tRscript",
	"sortmerna\t$fakeSortmerna",
	"SSUdbFAsrt\t$database{ssu_sort}",
	"LSUdbFAsrt\t$database{lsu_sort}",
	"SSUidx\t$indexDir",
	"LSUidx\t$indexDir",
	"vsearch\t$fakeVsearch",
	"flash\t$fakeFlash",
	"LCA\t$fakeLca",
	"SSUdbFA\t$database{ssu}",
	"SSUtax\t$database{ssu_tax}",
	"LSUdbFA\t$database{lsu}",
	"LSUtax\t$database{lsu_tax}",
)."\n");

my $readsDir = File::Spec->catdir($tmp, 'reads');
make_path($readsDir);
my $fastq = "\@read\nACGT\n+\nIIII\n";
my $r1 = File::Spec->catfile($readsDir, 'sample_R1.fastq');
my $r2 = File::Spec->catfile($readsDir, 'sample_R2.fastq');
my $single = File::Spec->catfile($readsDir, 'sample_single.fastq');
write_file($_, $fastq) for ($r1, $r2, $single);

my $align = File::Spec->catdir($tmp, 'aligned');
my $scratch = File::Spec->catdir($tmp, 'scratch');
my $sortLog = File::Spec->catfile($tmp, 'sortmerna.log');
{
	local $ENV{FAKE_SORTMERNA_LOG} = $sortLog;
	my ($status, $output, $errors) = run_script(
		'catchLSUSSU.pl',
		'-R1', $r1, '-R2', $r2, '-RS', $single,
		'-alignDir', $align, '-tmpDir', $scratch,
		'-smplID', 'sample', '-cores', 2, '-config', $config,
	);
	is($status, 0, 'catchLSUSSU completes for uncompressed paired and singleton FASTQ')
		or diag($output, $errors);
}
for my $tag (qw(SSU LSU)) {
	ok(-e File::Spec->catfile($align, $tag.'_pull.sto'),
		"$tag extraction checkpoint is created");
	for my $suffix (qw(r1.fq.gz r2.fq.gz fq.gz)) {
		my $path = File::Spec->catfile(
			$align, 'reads_'.$tag.'.'.$suffix,
		);
		ok(-s $path, "$tag $suffix is a nonempty gzip container");
		like(read_gzip($path), qr/(?:^$|^\@fake)/,
			"$tag $suffix is readable gzip output");
	}
}
my @sortRuns = log_lines($sortLog);
is(scalar(@sortRuns), 4,
	'one SortMeRNA process is used per marker and input library');
ok(!grep(/--kvdb|--readb/, @sortRuns),
	'SortMeRNA runs use isolated workdirs instead of manual internal stores');
ok(-e $database{ssu_sort} && !-s $database{ssu_sort},
	'existence-only reference validation accepts an empty test reference');

my $lcaDir = File::Spec->catdir($align, 'ltsLCA');
for my $tag (qw(SSU LSU)) {
	write_file(File::Spec->catfile($lcaDir, $tag.'_ass.sto'), "old\n");
	write_file(File::Spec->catfile(
		$lcaDir, $tag.'riboRun_bl.hiera.txt',
	), "old\n");
}
write_file(File::Spec->catfile($lcaDir, 'Assigned.sto'), "old\n");
unlink File::Spec->catfile($align, 'SSU_pull.sto')
	or die "Cannot remove SSU stone: $!";
unlink File::Spec->catfile($align, 'reads_SSU.fq.gz')
	or die "Cannot remove SSU single output: $!";
unlink $database{lsu_sort} or die "Cannot remove LSU test reference: $!";
{
	local $ENV{FAKE_SORTMERNA_LOG} = $sortLog;
	my ($status, $output, $errors) = run_script(
		'catchLSUSSU.pl',
		'-R1', $r1, '-R2', $r2, '-RS', $single,
		'-alignDir', $align, '-tmpDir', $scratch,
		'-smplID', 'sample', '-cores', 2, '-config', $config,
	);
	is($status, 0,
		'partial extraction recovery does not validate the completed LSU database')
		or diag($output, $errors);
}
@sortRuns = log_lines($sortLog);
is(scalar(@sortRuns), 6, 'partial recovery reruns only the missing SSU inputs');
ok(!-e File::Spec->catfile($lcaDir, 'SSU_ass.sto')
	&& !-e File::Spec->catfile($lcaDir, 'Assigned.sto'),
	'SSU extraction recovery invalidates stale marker and aggregate taxonomy state');
ok(-e File::Spec->catfile($lcaDir, 'LSU_ass.sto'),
	'SSU extraction recovery leaves completed LSU taxonomy state untouched');

unlink $database{ssu_sort} or die "Cannot remove SSU test reference: $!";
{
	local $ENV{FAKE_SORTMERNA_LOG} = $sortLog;
	my ($status, $output, $errors) = run_script(
		'catchLSUSSU.pl',
		'-R1', $r1, '-R2', $r2, '-RS', $single,
		'-alignDir', $align, '-tmpDir', $scratch,
		'-smplID', 'sample', '-cores', 2, '-config', $config,
	);
	is($status, 0,
		'a fully complete extraction exits before database or scheduler-side tool setup')
		or diag($output, $errors);
}
is(scalar(log_lines($sortLog)), 6,
	'a complete extraction performs no additional SortMeRNA calls');

write_file($database{ssu_sort}, '');
write_file($database{lsu_sort}, '');
my $missingAlign = File::Spec->catdir($tmp, 'missing-output');
{
	local $ENV{FAKE_SORTMERNA_MISSING_TAG} = 'SSU';
	my ($status, $output, $errors) = run_script(
		'catchLSUSSU.pl',
		'-RS', $single,
		'-alignDir', $missingAlign, '-tmpDir', $scratch,
		'-smplID', 'missing', '-cores', 1, '-config', $config,
	);
	isnt($status, 0, 'missing SortMeRNA output is a hard failure');
	like($errors, qr/produced no expected single output/,
		'missing SortMeRNA output has a precise diagnostic');
	ok(!-e File::Spec->catfile($missingAlign, 'SSU_pull.sto'),
		'missing SortMeRNA output cannot create a success checkpoint');
}

my $zeroAlign = File::Spec->catdir($tmp, 'zero-hit');
{
	local $ENV{FAKE_SORTMERNA_ZERO} = 1;
	my ($status, $output, $errors) = run_script(
		'catchLSUSSU.pl',
		'-RS', $single,
		'-alignDir', $zeroAlign, '-tmpDir', $scratch,
		'-smplID', 'zero', '-cores', 1, '-config', $config,
	);
	is($status, 0, 'a genuine zero-hit SortMeRNA run completes')
		or diag($output, $errors);
}
for my $tag (qw(SSU LSU)) {
	ok(-e File::Spec->catfile($zeroAlign, $tag.'_pull.sto'),
		"$tag zero-hit checkpoint is present");
	for my $suffix (qw(r1.fq.gz r2.fq.gz fq.gz)) {
		my $path = File::Spec->catfile(
			$zeroAlign, 'reads_'.$tag.'.'.$suffix,
		);
		ok(-s $path, "$tag zero-hit $suffix is a valid gzip container");
		is(read_gzip($path), '', "$tag zero-hit $suffix expands to no reads");
	}
}

my $vsearchLog = File::Spec->catfile($tmp, 'vsearch.log');
my $lcaLog = File::Spec->catfile($tmp, 'lca.log');
{
	local $ENV{FAKE_VSEARCH_LOG} = $vsearchLog;
	local $ENV{FAKE_LCA_LOG} = $lcaLog;
	my ($status, $output, $errors) = run_script(
		'lotus_LCA_blast3.pl',
		'-dir', $align, '-DBdir', $databaseDir,
		'-smplID', 'sample', '-pairedRds', 2,
		'-cores', 2, '-simMode', 4, '-config', $config,
	);
	is($status, 0, 'lotus LCA completes in paired-plus-single VSEARCH mode')
		or diag($output, $errors);
}
ok(-e File::Spec->catfile($lcaDir, 'Assigned.sto'),
	'LCA creates its aggregate checkpoint only after both markers finish');
for my $tag (qw(SSU LSU)) {
	ok(-e File::Spec->catfile($lcaDir, $tag.'_ass.sto'),
		"$tag assignment checkpoint is present");
	ok(-s File::Spec->catfile(
		$lcaDir, $tag.'riboRun_bl.hiera.txt',
	), "$tag hierarchy is nonempty");
}
ok(!-e $database{ssu}.'.vudb' && !-e $database{lsu}.'.vudb',
	'VSEARCH fallback indexes are built in job-local scratch, not beside shared databases');
my @lcaRuns = log_lines($lcaLog);
ok(@lcaRuns && !grep(/\.m8\.gz(?:,|$)/, @lcaRuns),
	'LCA receives plain .m8 files rather than falsely named uncompressed gzip files');

unlink File::Spec->catfile($lcaDir, 'SSU_ass.sto')
	or die "Cannot remove SSU assignment stone: $!";
unlink File::Spec->catfile($lcaDir, 'SSUriboRun_bl.hiera.txt')
	or die "Cannot remove SSU hierarchy: $!";
unlink $database{lsu} or die "Cannot remove LSU database: $!";
{
	local $ENV{FAKE_VSEARCH_LOG} = $vsearchLog;
	local $ENV{FAKE_LCA_LOG} = $lcaLog;
	my ($status, $output, $errors) = run_script(
		'lotus_LCA_blast3.pl',
		'-dir', $align, '-DBdir', $databaseDir,
		'-smplID', 'sample', '-pairedRds', 2,
		'-cores', 2, '-simMode', 4, '-config', $config,
	);
	is($status, 0,
		'partial LCA recovery does not validate the completed LSU database')
		or diag($output, $errors);
}
ok(-e File::Spec->catfile($lcaDir, 'Assigned.sto'),
	'partial LCA recovery restores the aggregate checkpoint');

my $badRibo = File::Spec->catdir($tmp, 'bad-ribo');
my $badFastq = "\@bad\nACGT\n+\nIII\n";
for my $tag (qw(SSU LSU)) {
	write_gzip(
		File::Spec->catfile($badRibo, 'reads_'.$tag.'.fq.gz'),
		$tag eq 'SSU' ? $badFastq : $fastq,
	);
}
{
	my ($status, $output, $errors) = run_script(
		'lotus_LCA_blast3.pl',
		'-dir', $badRibo, '-DBdir', $databaseDir,
		'-smplID', 'bad', '-pairedRds', 0,
		'-cores', 1, '-simMode', 4, '-config', $config,
	);
	isnt($status, 0, 'malformed FASTQ is rejected before similarity search');
	like($errors, qr/Sequence\/quality length mismatch/,
		'malformed FASTQ has a record-specific diagnostic');
	ok(!-e File::Spec->catfile(
		$badRibo, 'ltsLCA', 'SSU_ass.sto',
	), 'malformed FASTQ cannot create an assignment checkpoint');
}

is(system($^X, '-I'.$root, '-c',
	File::Spec->catfile($scripts, 'catchLSUSSU.pl')), 0,
	'catchLSUSSU compiles');
is(system($^X, '-I'.$root, '-c',
	File::Spec->catfile($scripts, 'lotus_LCA_blast3.pl')), 0,
	'lotus_LCA_blast3 compiles');

done_testing;
