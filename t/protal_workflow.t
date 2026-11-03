use strict;
use warnings;

use Test::More;
use File::Temp qw(tempdir);
use File::Path qw(make_path);
use File::Spec;
use FindBin qw($Bin);
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Compress::Bzip2 qw(bzip2 $Bzip2Error);
use lib "$Bin/..", "$Bin/lib";
use MFTestConfig;
use Mods::SampleCompletion qw(completion_request_signature);
use Mods::ReadLibrary qw(newReadLibrary readLibrariesByScope protalReadInput);
use Mods::WorkflowControl qw(normalise_job_dependencies sample_is_ignored);
use Mods::WorkflowResilience qw(atomic_write_text retry_unlink);

# The Protal subs of MATAF4.pl run here against a stand-in protal: per-sample
# jobs, local batch maps, the cohort merge and the strain MSAs from kept SAMs.

sub read_file { my ($p) = @_; open my $fh, '<', $p or die "$p: $!"; local $/; my $c = <$fh>; close $fh; return $c; }
sub write_file { my ($p, $c) = @_; open my $fh, '>', $p or die "$p: $!"; binmode $fh; print {$fh} $c; close $fh or die $!; }

my $root = File::Spec->catdir($Bin, '..');
my $source = read_file(File::Spec->catfile($root, 'MATAF4.pl'));

our (%MFopt, %MFconfig, %MFglobal, %map, %protalSampleJobs, %protalBatchQueue);
our ($QSBoptHR, $logDir, $JNUM, $selectedFrom, $selectedTo, $controllerBaseOut,
	$pigzBin, $protalBatchCount, $ignoredSamplesHR);
our @samples;
for my $name (qw(_shell_quote _shell_command pathIsStrictChild
		protalSamplePaths protalSampleComplete protalSkipMatches protalCohortDir
		protalVersionPreflight protalFileCompression protalSampleInput protalReadStream
		protalBaseCommand protalDatabaseCheck protalClearSampleCommand
		protalPublishSampleCommand _protalMapField protalSampleJob queueProtalBatch
		flushProtalBatches submitProtalBatch mergeProtalProfiles submitProtalStrains)) {
	my ($body) = $source =~ /(^sub \Q$name\E\b[^\n]*\{.*?^\})/ms;
	die "Missing sub $name in MATAF4.pl\n" unless defined $body;
	eval $body;
	die "$name: $@" if $@;
}

my $tmp = tempdir(CLEANUP => 1);
my $fake = "$tmp/protal";
write_file($fake, <<'FAKE');
#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
# Stand-in for protal: with -1/-2 or a map it "aligns" a sample by counting its
# reads into the SAM, writes a profile and its logs, profiles an existing SAM
# instead of aligning, and writes strain output unless --no_strains.
open my $log, '>>', $ENV{FAKE_PROTAL_LOG} or die "log: $!";
print {$log} join("\t", @ARGV), "\n";
close $log;
my (%opt, %flag);
my @args = @ARGV;
while (@args) {
	my $arg = shift @args;
	if ($arg eq '--version') { print "protal v$ENV{FAKE_PROTAL_VERSION} (commit 0)\n"; exit 0; }
	if ($arg =~ /^--(no_strains|force)$/) { $flag{$1} = 1; next; }
	$opt{$arg} = shift @args;
}
sub reads_in {
	my ($path) = @_;
	open my $fh, '<:raw', $path or die "fake protal cannot read $path: $!\n";
	local $/;
	my $data = <$fh>;
	close $fh;
	if (substr($data, 0, 2) eq "\x1f\x8b") {
		my $plain;
		gunzip(\$data => \$plain, MultiStream => 1) or die "gunzip: $GunzipError\n";
		$data = $plain;
	}
	return scalar(() = $data =~ /\n/g) / 4;
}
sub put { my ($p, $c) = @_; make_path(dirname($p)); open my $fh, '>', $p or die "$p: $!"; print {$fh} $c; close $fh; }
sub profile { my ($p) = @_; put($p, "1\td__Bacteria;s__Fake\t100\n"); put("$p$_", "log\n") for qw(.log .gene.log .genes.log); }
sub align {
	my ($type, $first, $second, $sam, $profile) = @_;
	my $reads = reads_in($first) + (defined($second) ? reads_in($second) : 0);
	put($sam, "\@CO\ttype=$type\treads=$reads\n");
	put("$sam.err", '');
	profile($profile);
}
if (defined $opt{'--map'}) {
	open my $fh, '<', $opt{'--map'} or die "map: $!";
	my (%var, @columns, @names, $failed);
	while (my $line = <$fh>) {
		chomp $line;
		my @cells = split /\t/, $line;
		if ($cells[0] eq '#SAMPLEID') { @columns = @cells; next; }
		if ($cells[0] =~ /^#/) { $var{$cells[0]} = $cells[1]; next; }
		my %row;
		@row{@columns} = @cells;
		my $out = $var{'#OUTPUT_DIR'};
		my $sam = $row{SAM} =~ m{^/} ? $row{SAM} : "$out/alignments/$row{SAM}";
		my $profile = $row{PROFILE} =~ m{^/} ? $row{PROFILE} : "$out/profiles/$row{PROFILE}";
		if (-e $sam) {
			profile($profile);
			put("$sam.err", '');
		} elsif ($row{FIRST} =~ /FAIL/) {
			$failed++;
			next;
		} else {
			my $second = defined($row{SECOND}) && $row{SECOND} ne '-' ? $row{SECOND} : undef;
			align($row{READ_TYPE} || 'pe', $row{FIRST}, $second, $sam, $profile);
		}
		push @names, $row{'#SAMPLEID'};
	}
	put(($var{'#STRAIN_OUTPUT_DIR'} || "$var{'#OUTPUT_DIR'}/strains").'/species.tsv', join("\n", @names)."\n")
		unless $flag{no_strains};
	exit($failed ? 1 : 0);
}
my $prefix = $opt{'--prefix'};
align($opt{'--read_type'}, $opt{'-1'}, $opt{'-2'}, "$prefix.sam.zst",
	($opt{'--profile_dir'} // dirname($prefix)).'/'.basename($prefix).'.profile');
FAKE
chmod 0755, $fake;
$ENV{FAKE_PROTAL_LOG} = "$tmp/protal.log";
$ENV{FAKE_PROTAL_VERSION} = '0.7.8';

my %progs = (protal => $fake, protal_db => '/db/protal', protalProfileUtils => 'protal_profile_utils',
	bzip2 => 'bzip2');
sub getProgPaths { my ($name) = @_; return exists $progs{$name} ? $progs{$name} : "prog_$name"; }
my @jobs;
sub qsubSystem {
	push @jobs, {script => $_[0], cmd => $_[1], cores => $_[2], memory => $_[3], name => $_[4],
		deps => $_[5], tmp_space => $QSBoptHR->{tmpSpace}};
	return ('J'.scalar(@jobs), '');
}
my %readsets;
sub sampleReadSet { return $readsets{$_[0]}{$_[1]}; }
my %downloaded;
sub sampleAccessionSource { return $downloaded{$_[0]} ? ('ena', 'ERR1') : ('', ''); }
my %locks;
sub recordSampleLockJobs { my ($lock, $jobsAR) = @_; push @{$locks{$lock}}, @{$jobsAR}; return 1; }
sub sampleCompletionRequestSignature { return "request-$_[0]"; }

# Mode 2 (strain MSAs) keeps the SAMs that the cohort subtest builds the MSAs from;
# the mode-1 subtest checks that they are discarded otherwise.
%MFopt = (DoProtal => 2, ProtalCores => 8, ProtalMem => 60, protalIgnoreErrors => 1, protalBatchSize => 0);
%MFconfig = (silent => 1, inspectState => 0);
%MFglobal = (runTmpDirGlobal => "$tmp/scratch/");
$QSBoptHR = {tmpSpace => '99G', LOCKfile => 'controller.lock'};
$logDir = "$tmp/log/";
$JNUM = 3;
$controllerBaseOut = "$tmp/out";
$pigzBin = 'gzip';
$protalBatchCount = 0;
$ignoredSamplesHR = {};
make_path($logDir, $MFglobal{runTmpDirGlobal});

# Read files and samples. Each FASTQ record is 4 lines.
my $reads = sub { my ($n) = @_; return join('', map { "\@r$_\nACGT\n+\nIIII\n" } 1 .. $n); };
my $gz = sub { my ($path, $n) = @_; my $data = $reads->($n); gzip(\$data => $path) or die $GzipError; return $path; };
make_path("$tmp/in");
sub add_sample {
	my ($key, %args) = @_;
	my $sampleRoot = "$tmp/out/$key/";
	make_path($sampleRoot);
	$map{$key} = {SmplID => $key, wrdir => $sampleRoot};
	my @libraries;
	my $i = 0;
	for my $library (@{$args{libraries}}) {
		push @libraries, newReadLibrary(
			id => "$key:primary:".$i++, sample => $key, scope => 'primary', phase => 'staged',
			technology => $library->{technology} || 'ill',
			files => {map { $_ => "$tmp/staged/".($library->{$_} =~ s{.*/}{}r) } grep { $_ ne 'technology' } keys %{$library}},
			source_files => {map { $_ => $library->{$_} } grep { $_ ne 'technology' } keys %{$library}},
		);
	}
	$readsets{$key}{raw} = {libraries => \@libraries};
	$downloaded{$key} = $args{downloaded} ? 1 : 0;
	push @samples, $key;
	return $sampleRoot;
}
sub run_job {
	my ($job) = @_;
	my $script = "$tmp/run.".($job->{name}).'.sh';
	write_file($script, "#!/bin/bash\nset -eo pipefail\n".$job->{cmd});
	return system('bash', '-c', 'bash "$1" >"$2" 2>&1', 'run', $script, "$script.log") >> 8;
}

subtest 'read streams: files as they are, joined runs, decompression' => sub {
	is(protalReadStream('/in/a.fq.gz'), q('/in/a.fq.gz'), 'a single gzip file is passed as it is');
	is(protalReadStream('/in/a.fq.gz', '/in/b.fq.gz'), q(<('cat' '/in/a.fq.gz' '/in/b.fq.gz')),
		'gzip runs are concatenated into one stream');
	is(protalReadStream('/in/a,b.fq'), q(<('cat' '/in/a,b.fq')),
		'a path with a comma goes through a pipe (protal splits lists at commas)');
	is(protalReadStream('/in/a.fq.gz', '/in/b.fq.bz2'), q(<('gzip' '-dc' '/in/a.fq.gz'; 'bzip2' '-dc' '/in/b.fq.bz2')),
		'mixed or unreadable compressions are decompressed');
	is(protalFileCompression('x.fastq.zst'), 'zstd', 'zstd is read by protal itself');
	is(protalFileCompression('x.fq.xz'), 'xz', 'xz needs a pipe');
};

my $localRoot = add_sample('LOC', libraries => [{r1 => $gz->("$tmp/in/LOC_1.fq.gz", 3), r2 => $gz->("$tmp/in/LOC_2.fq.gz", 3)}]);
my $twoRunRoot = add_sample('RUNS', downloaded => 1, libraries => [
	{r1 => $gz->("$tmp/in/R1_1.fq.gz", 2), r2 => $gz->("$tmp/in/R1_2.fq.gz", 2)},
	{r1 => $gz->("$tmp/in/R2_1.fq.gz", 5), r2 => $gz->("$tmp/in/R2_2.fq.gz", 5)},
]);
my $bz2Data = $reads->(4);
bzip2(\$bz2Data => "$tmp/in/BZ.fq.bz2") or die $Bzip2Error;
my $bz2Root = add_sample('BZ', libraries => [{single => "$tmp/in/BZ.fq.bz2"}]);
my $longRoot = add_sample('ONT', libraries => [{single => $gz->("$tmp/in/ONT.fq.gz", 2), technology => 'ONT'}]);
my $bamRoot = add_sample('BAM', libraries => [{single => "$tmp/in/cache.fq.gz", bam => "$tmp/in/x.bam"}]);

subtest 'input plan: local samples batch, downloaded samples run alone' => sub {
	my $local = protalSampleInput('LOC', 'LOC', $localRoot, 'request-LOC');
	ok($local->{batch}, 'a local sample with one gzip pair joins the batch');
	is($local->{read_type}, 'pe', 'paired-end reads');
	is_deeply($local->{r1}, ["$tmp/in/LOC_1.fq.gz"], 'the source file is read, not the staged copy');
	my $runs = protalSampleInput('RUNS', 'RUNS', $twoRunRoot, 'request-RUNS');
	ok(!$runs->{batch}, 'a downloaded sample gets a job of its own');
	is(scalar(@{$runs->{r1}}), 2, 'both runs of the downloaded sample are profiled');
	ok(!protalSampleInput('BZ', 'BZ', $bz2Root, 'request-BZ')->{batch},
		'a local bzip2 sample needs decompressing, so it runs alone');
	is(protalSampleInput('ONT', 'ONT', $longRoot, 'request-ONT')->{read_type}, 'ont',
		'long reads keep their kind');
	my $bam;
	my @warnings;
	{
		local $SIG{__WARN__} = sub { push @warnings, @_ };
		$bam = protalSampleInput('BAM', 'BAM', $bamRoot, 'request-BAM');
	}
	ok($bam->{skipped}, 'an alignment-only sample is skipped in tolerant mode');
	my $skip = protalSamplePaths($bamRoot, 'BAM')->{skip};
	ok(protalSkipMatches($skip, 'request-BAM'), 'the skip marker is scoped to the request');
	ok(!protalSkipMatches($skip, 'request-other'), 'another request does not accept the marker');
	like(join('', @warnings), qr/skipping Protal for this sample/, 'the skip is reported');
	make_path(protalSamplePaths($localRoot, 'LOC')->{dir});
	write_file(protalSamplePaths($localRoot, 'LOC')->{skip}, "old\nreason\n");
	protalSampleInput('LOC', 'LOC', $localRoot, 'request-LOC');
	ok(!-e protalSamplePaths($localRoot, 'LOC')->{skip}, 'a stale skip marker is removed');
	local $MFopt{protalIgnoreErrors} = 0;
	eval { protalSampleInput('BAM', 'BAM', $bamRoot, 'request-BAM') };
	like($@, qr/no read library protal can profile/, 'strict mode stops on unusable input');
};

subtest 'per-sample job (mode 2): SAM and profile kept in the sample output' => sub {
	@jobs = ();
	my $input = protalSampleInput('RUNS', 'RUNS', $twoRunRoot, 'request-RUNS');
	my $scratch = "$tmp/scratch/RUNS/Protal/";
	my $job = protalSampleJob('RUNS', $twoRunRoot, $input, $scratch);
	is($job, 'J1', 'the job is submitted');
	my ($submitted) = @jobs;
	is($submitted->{name}, 'PT3', 'named after the sample index');
	is($submitted->{deps}, '', 'it waits for no staging job');
	is($submitted->{tmp_space}, 0, 'it reserves no node-local scratch');
	is($QSBoptHR->{tmpSpace}, '99G', 'the scratch request is restored');
	my $paths = protalSamplePaths($twoRunRoot, 'RUNS');
	like($submitted->{cmd}, qr/'--read_type' 'pe' -1 <\('cat' '\Q$tmp\E\/in\/R1_1\.fq\.gz' '\Q$tmp\E\/in\/R2_1\.fq\.gz'\) -2 <\('cat'/,
		'the runs are joined mate by mate');
	like($submitted->{cmd}, qr/'--prefix' '\Q$paths->{dir}\E\/RUNS' '--outdir' '\Q$scratch\E' '--profile_dir' '\Q$paths->{profile_dir}\E' '--threads' '8' '--no_strains'/,
		'SAM and profile go to the sample output, the rest to scratch');
	unlike($submitted->{cmd}, qr/--force/, 'no --force: the outputs are removed instead');
	write_file($paths->{sam}, 'old');
	write_file($paths->{stone}, '');
	is(run_job($submitted), 0, 'the job runs');
	ok(-e $paths->{stone}, 'the stone marks the sample complete');
	like(read_file($paths->{sam}), qr/type=pe\treads=14/, 'the SAM is new and holds the reads of both runs');
	ok(-e $paths->{profile}, 'the profile is kept');
	ok(!-e $paths->{diagnostics}[0], 'the profile logs are removed');
	ok(!-e $paths->{sam_errors}, 'an empty read-error list is removed');
	ok(!-d $scratch, 'the scratch directory is removed');
	ok(protalSampleComplete($twoRunRoot, 'RUNS'), 'the sample counts as complete');

	@jobs = ();
	protalSampleJob('BZ', $bz2Root, protalSampleInput('BZ', 'BZ', $bz2Root, 'request-BZ'), "$tmp/scratch/BZ/Protal/");
	like($jobs[0]{cmd}, qr/'--read_type' 'se' -1 <\('bzip2' '-dc' /, 'bzip2 reads are decompressed through a pipe');
	unlike($jobs[0]{cmd}, qr/ -2 /, 'single-end reads have no second file');
	is(run_job($jobs[0]), 0, 'the single-end job runs');
	like(read_file(protalSamplePaths($bz2Root, 'BZ')->{sam}), qr/reads=4/, 'the decompressed reads were aligned');
	eval { protalSampleJob('BZ', $bz2Root, {}, "$tmp/elsewhere/") };
	like($@, qr/Unsafe Protal scratch directory/, 'the scratch directory must be a Protal directory');
};

subtest 'batch map of local samples (mode 2)' => sub {
	@jobs = (); %locks = ();
	my $failRoot = add_sample('FAILS', libraries => [{r1 => $gz->("$tmp/in/FAIL_1.fq.gz", 1), r2 => $gz->("$tmp/in/FAIL_2.fq.gz", 1)}]);
	my $seRoot = add_sample('SE', libraries => [{single => $gz->("$tmp/in/SE.fq.gz", 6)}]);
	for my $sample (['LOC', $localRoot], ['FAILS', $failRoot], ['SE', $seRoot]) {
		my ($key, $sampleRoot) = @{$sample};
		my $input = protalSampleInput($key, $key, $sampleRoot, "request-$key");
		ok($input->{batch}, "$key joins the batch");
		queueProtalBatch($key, $key, $sampleRoot, $input, "$sampleRoot/lock");
	}
	my @batchJobs = flushProtalBatches();
	is(scalar(@batchJobs), 1, 'one job for all local samples');
	is_deeply(\%protalBatchQueue, {}, 'the queue is empty afterwards');
	is($locks{"$localRoot/lock"}[0], $batchJobs[0], 'each member lock records the batch job');
	is($protalSampleJobs{SE}, $batchJobs[0], 'the merge waits for the batch job');
	my $job = $jobs[0];
	is($QSBoptHR->{LOCKfile}, 'controller.lock', 'the lock setting is restored');
	my ($mapFile) = $job->{cmd} =~ /'--map' '([^']+)'/;
	ok(defined($mapFile) && -s $mapFile, 'the map is written');
	my @map = split /\n/, read_file($mapFile);
	is($map[3], join("\t", '#SAMPLEID', qw(PREFIX FIRST SECOND SAM PROFILE READ_TYPE)), 'map columns');
	my %rows = map { my @c = split /\t/; ($c[0] => \@c) } grep { !/^#/ } @map;
	is($rows{LOC}[4], protalSamplePaths($localRoot, 'LOC')->{sam}, 'the SAM path is absolute, in the sample output');
	is($rows{LOC}[5], protalSamplePaths($localRoot, 'LOC')->{profile}, 'so is the profile');
	is_deeply([@{$rows{SE}}[3, 6]], ['-', 'se'], 'a single-end sample has no second file and its read type');
	like($job->{cmd}, qr/'--no_strains'/, 'no strain MSAs in the batch');
	unlike($job->{cmd}, qr/'--profile_dir'/, 'profiles are not moved to one directory');
	is(run_job($job), 1, 'the batch job fails when one of its samples fails');
	ok(protalSampleComplete($localRoot, 'LOC'), 'the samples that worked are published nevertheless');
	ok(protalSampleComplete($seRoot, 'SE'), 'also the single-end one');
	like(read_file(protalSamplePaths($seRoot, 'SE')->{sam}), qr/type=se\treads=6/, 'with its read type');
	ok(!-e protalSamplePaths($failRoot, 'FAILS')->{stone}, 'the failed sample is not');

	@jobs = ();
	local $MFopt{protalBatchSize} = 2;
	queueProtalBatch($_, $_, $map{$_}{wrdir}, protalSampleInput($_, $_, $map{$_}{wrdir}, "request-$_"), '')
		for qw(LOC FAILS SE);
	is(scalar(flushProtalBatches()), 2, '-protalBatchSize splits the batch');
	is(scalar(grep { /^[^#]/ } split /\n/, read_file(($jobs[1]{cmd} =~ /'--map' '([^']+)'/)[0])), 1,
		'the last batch takes the rest');
};

subtest 'without strain MSAs (mode 1) the SAM is not kept' => sub {
	local $MFopt{DoProtal} = 1;
	local @samples = @samples;
	local ($selectedFrom, $selectedTo) = ($selectedFrom, $selectedTo);
	@jobs = ();
	my $downloadRoot = add_sample('DL1', downloaded => 1, libraries => [
		{r1 => $gz->("$tmp/in/DL1_1.fq.gz", 2), r2 => $gz->("$tmp/in/DL1_2.fq.gz", 2)}]);
	my $paths = protalSamplePaths($downloadRoot, 'DL1');
	make_path($paths->{dir});
	write_file($paths->{sam}, 'kept by an earlier mode-2 run');
	my $scratch = "$tmp/scratch/DL1/Protal/";
	protalSampleJob('DL1', $downloadRoot,
		protalSampleInput('DL1', 'DL1', $downloadRoot, 'request-DL1'), $scratch);
	like($jobs[0]{cmd}, qr/'--prefix' '\Q$scratch\EDL1' /, 'the per-sample job writes the SAM to scratch');
	is(run_job($jobs[0]), 0, 'the per-sample job runs');
	ok(-e $paths->{profile} && -e $paths->{stone}, 'the profile is kept and the sample published');
	ok(!-e $paths->{sam}, 'no SAM is kept, and an earlier one is removed');
	ok(!-d $scratch, 'the SAM is deleted with the scratch');
	ok(protalSampleComplete($downloadRoot, 'DL1'), 'the sample is complete without a SAM');
	{
		local $MFopt{DoProtal} = 2;
		ok(!protalSampleComplete($downloadRoot, 'DL1'), 'mode 2 would align it again for its SAM');
	}

	@jobs = ();
	my $batchRoot = add_sample('LB1', libraries => [
		{r1 => $gz->("$tmp/in/LB1_1.fq.gz", 3), r2 => $gz->("$tmp/in/LB1_2.fq.gz", 3)}]);
	queueProtalBatch('LB1', 'LB1', $batchRoot,
		protalSampleInput('LB1', 'LB1', $batchRoot, 'request-LB1'), '');
	flushProtalBatches();
	my ($mapFile) = $jobs[0]{cmd} =~ /'--map' '([^']+)'/;
	my $map = read_file($mapFile);
	my ($workDir) = $map =~ /^#OUTPUT_DIR\t(.*)$/m;
	my ($row) = grep { /^LB1\t/ } split /\n/, $map;
	is((split /\t/, $row)[4], File::Spec->catfile($workDir, 'alignments', 'LB1.sam.zst'),
		'the batch writes the SAM to its scratch');
	is(run_job($jobs[0]), 0, 'the batch runs');
	ok(protalSampleComplete($batchRoot, 'LB1'), 'the batch sample is published');
	ok(!-e protalSamplePaths($batchRoot, 'LB1')->{sam}, 'without a SAM in its output');
	ok(!-d $workDir, 'the batch scratch, with the SAM, is removed');

	@jobs = ();
	@samples = qw(DL1 LB1);
	($selectedFrom, $selectedTo) = (0, 2);
	mergeProtalProfiles();
	is(scalar(@jobs), 1, 'mode 1 merges the profiles and builds no strain MSAs');
};

subtest 'cohort: merge and strain MSAs from the kept SAMs' => sub {
	@jobs = (); %protalSampleJobs = ();
	@samples = qw(LOC RUNS BZ);
	($selectedFrom, $selectedTo) = (0, 3);
	my $out = "$tmp/out/pseudoGC/protal";
	my $merge = mergeProtalProfiles();
	is(scalar(@jobs), 2, 'merge and strain jobs');
	like($jobs[0]{cmd}, qr/'protal_profile_utils' 'merge' '--input' .*LOC\.profile.*RUNS\.profile.*BZ\.profile/,
		'the profiles are merged');
	like($jobs[0]{cmd}, qr/'\Q$out\E\/Protal\.abundance\.tsv'/, 'into pseudoGC/protal');
	my $strains = $jobs[1];
	like($strains->{cmd}, qr/'--map' '[^']+'/, 'the strain job runs a map');
	unlike($strains->{cmd}, qr/--no_strains|--force/, 'with strains, without aligning again');
	my ($mapFile) = $strains->{cmd} =~ /'--map' '([^']+)'/;
	my @rows = grep { !/^#/ } split /\n/, read_file($mapFile);
	is(scalar(@rows), 3, 'one row per sample');
	my @loc = split /\t/, $rows[0];
	is($loc[4], protalSamplePaths($localRoot, 'LOC')->{sam}, 'the kept SAM');
	is($loc[3], '-', 'no reads are needed');
	is(run_job($strains), 0, 'the strain job runs');
	ok(-l "$out/strains/current" || -e "$out/strains/current", 'strains/current points at the cohort');
	my $signature = read_file("$out/Protal.strains.current");
	chomp $signature;
	is(read_file("$out/strains/$signature/species.tsv"), "LOC\nRUNS\nBZ\n", 'the MSAs cover the cohort');
	ok(-e "$out/Protal.strains.$signature.sto", 'the cohort stone is written');
	like(read_file(protalSamplePaths($localRoot, 'LOC')->{sam}), qr/reads=6\n/, 'the SAMs are left as they were');

	@jobs = ();
	mergeProtalProfiles();
	is(scalar(@jobs), 1, 'complete strain MSAs of the same cohort are not built again');

	@jobs = ();
	($selectedFrom, $selectedTo) = (0, 2);
	mergeProtalProfiles();
	is(scalar(@jobs), 1, 'no strain MSAs for part of the cohort');

	@jobs = ();
	($selectedFrom, $selectedTo) = (0, 3);
	add_sample('NEW', libraries => [{r1 => "$tmp/in/LOC_1.fq.gz", r2 => "$tmp/in/LOC_2.fq.gz"}]);
	@samples = qw(LOC RUNS NEW);
	my @warnings;
	{
		local $SIG{__WARN__} = sub { push @warnings, @_ };
		mergeProtalProfiles();
	}
	is(scalar(@jobs), 0, 'a sample without profile or job defers the cohort outputs');
	like(join('', @warnings), qr/Deferring Protal cohort outputs: 1 .*first: NEW/, 'and says which');
	$protalSampleJobs{NEW} = 'J77';
	mergeProtalProfiles();
	is($jobs[0]{deps}, 'J77', 'the merge waits for the pending job');
	@samples = qw(LOC RUNS BZ BAM);
	($selectedFrom, $selectedTo) = (0, 4);
	@jobs = ();
	mergeProtalProfiles();
	unlike($jobs[0]{cmd}, qr/BAM\.profile/, 'a sample skipped for its input is left out');
};

subtest 'protal version' => sub {
	local $ENV{FAKE_PROTAL_VERSION} = '0.6.0';
	eval { protalVersionPreflight() };
	like($@, qr/needs protal 0\.7\.8 or later, but .* is 0\.6\.0/, 'protal 0.6 is refused');
	$ENV{FAKE_PROTAL_VERSION} = '0.7.8';
	ok(eval { protalVersionPreflight(); 1 }, '0.7.8 is accepted');
	$ENV{FAKE_PROTAL_VERSION} = '0.10.0';
	ok(eval { protalVersionPreflight(); 1 }, 'later versions are accepted');
};

subtest 'wiring in the controller' => sub {
	like($source, qr/"protalBatchSize=i" => \\\$MFopt\{protalBatchSize\}/, '-protalBatchSize is parsed');
	like($source, qr/DoMOTU2 DoMetaPhlan DoProtal protalIgnoreErrors/,
		'the mode is part of the sample signature: mode 2 needs the SAMs mode 1 discards');
	like($source, qr/\$configuration\{protal_contract\} = 2;/, 'the signature records the Protal contract');
	unlike($source, qr/\$requested\{DoProtal\} =/, 'the mode is not normalised away');
	like($source, qr/if \(\$calcProtal\) \{\s*my \$protalInput = protalSampleInput\(.*?queueProtalBatch\(.*?protalSampleJob\(.*?add2SampleDeps\(\\\@sampleDeps, \[\$protalJob\]\).*?deferLoopProducerWave\(\s*'input staging'/s,
		'Protal is planned before the staging deferral, and its own jobs hold the sample scratch');
	unlike($source, qr/\$calcUnzip=1 if \([^;]*\$calcProtal/, 'a pending profile does not stage the reads');
	like($source, qr/if \(\$JNUM == \(\$to-1\)\)\{.*?push \@grandDeps, flushProtalBatches\(\);/s,
		'a loop pass submits and waits for its batch');
	like($source, qr/sub postprocess\{.*?mergeProtalProfiles\(\)/s, 'cohort outputs in postprocessing');
	like($source, qr/\{id => 'alignments', kind => 'nonempty', path => \$protal->\{sam\}\}\s*if \$MFopt\{DoProtal\} == 2/,
		'sample completion requires the SAM in mode 2 only');
};

done_testing();
