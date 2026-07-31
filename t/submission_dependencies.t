use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::Subm qw(
	qsubSystem qsubSystemJobAlive qsubSystemWaitMaxJobs numLiveUserJobs numActiveUserJobs
	reconcileSlurmDependencies MFnext recordSampleLockJobs sampleLockActiveJobs
	primeSampleLockJobSnapshot
	slurmJobFailureSummary
	submitSlurmWithDependencyRecovery
	deferredSubmissionDependency submissionDependencyDeferred
	handleSubmissionFailure
);

sub slurm_options {
	return {
		rTag => 'run', qmode => 'slurm', doSubmit => 0, doSync => 0,
		tmpSpaceTag => '', tmpSpace => 0, submissionConfig => '', constraint => [],
		LOCKfile => '', excludeNodes => '', medQueue => 'normal', medTime => '',
		longQueue => 'normal', highMemQueue => 'normal', gpuQueue => 'normal',
		netQueue => 'normal', shortQueue => 'normal', longTime => '',
		useLongQueue => 0, useHiMemQueue => 0, useGPUQueue => 0,
		useNetQueue => 0, useShortQueue => 0, gpuCount => 0,
		wcKeysForJob => '', xtraNodeCmds => '', afterAny => 0,
	};
}

sub bash_options {
	my $options = slurm_options();
	$options->{qmode} = 'bash';
	$options->{doSubmit} = 1;
	$options->{submittedJobs} = 0;
	return $options;
}

my $root = tempdir(CLEANUP => 1);

my $sample_lock = File::Spec->catfile($root, 'sample.lock');
my $lock_options = slurm_options();
$lock_options->{rTag} = 'run';
recordSampleLockJobs($sample_lock, ['run701', '702', 'not-a-job'], $lock_options);
open my $sample_lock_fh, '<', $sample_lock or die "Cannot read $sample_lock: $!";
my $sample_lock_text = do { local $/; <$sample_lock_fh> };
close $sample_lock_fh;
is($sample_lock_text, "701\n702\n",
	'sample locks persist only concrete scheduler job IDs');
my $lock_queries = 0;
$lock_options->{sampleLockCheckInterval} = 60;
$lock_options->{schedulerClock} = sub { 1000 };
$lock_options->{sampleLockJobRunner} = sub {
	$lock_queries++;
	return ("702\n999\n", 0);
};
is(sampleLockActiveJobs($sample_lock, $lock_options), 1,
	'a sample lock remains active while one recorded job is in the scheduler');
is(sampleLockActiveJobs($sample_lock, $lock_options), 1,
	'a scheduler snapshot is shared across repeated sample-lock checks');
is($lock_queries, 1, 'sample-lock checks use one cached scheduler query');
recordSampleLockJobs($sample_lock, ['703'], $lock_options);
is(sampleLockActiveJobs($sample_lock, $lock_options), 2,
	'a locally submitted job is merged safely into the cached active set');

my $finished_lock = File::Spec->catfile($root, 'finished.lock');
recordSampleLockJobs($finished_lock, ['704'], $lock_options);
is(sampleLockActiveJobs($finished_lock, $lock_options), 1,
	'a newly recorded job remains active even with an older scheduler snapshot');
delete $lock_options->{sampleLockActiveState};
$lock_options->{sampleLockJobRunner} = sub { return ("999\n", 0); };
is(sampleLockActiveJobs($finished_lock, $lock_options), 0,
	'a lock is releasable after all recorded jobs leave the scheduler');

my $legacy_lock = File::Spec->catfile($root, 'legacy.lock');
open my $legacy_fh, '>', $legacy_lock or die "Cannot create $legacy_lock: $!";
close $legacy_fh;
is(sampleLockActiveJobs($legacy_lock, $lock_options), undef,
	'an empty legacy lock is retained when its ownership cannot be verified');

my $batch_lock_a = File::Spec->catfile($root, 'batch-a.lock');
my $batch_lock_b = File::Spec->catfile($root, 'batch-b.lock');
recordSampleLockJobs($batch_lock_a, [901, 902], $lock_options);
recordSampleLockJobs($batch_lock_b, [903], $lock_options);
my ($batch_queue_calls, $batch_accounting_calls) = (0, 0);
my $batch_options = slurm_options();
$batch_options->{sampleLockJobRunner} = sub {
	$batch_queue_calls++;
	return ("901|RUNNING\n999|PENDING\n", 0);
};
$batch_options->{sampleLockAccountingRunner} = sub {
	my ($command, $ids) = @_;
	$batch_accounting_calls++;
	is_deeply($ids, [902, 903],
		'one accounting batch receives every tracked job absent from squeue');
	return ("902|COMPLETED|0:0\n903|FAILED|1:0\n", 0);
};
my $batch_snapshot = primeSampleLockJobSnapshot(
	[$batch_lock_a, $batch_lock_b], $batch_options,
);
is_deeply($batch_snapshot, { tracked => 3, queued => 1, accounted => 2 },
	'the pass snapshot combines scheduler and accounting state');
is(sampleLockActiveJobs($batch_lock_a, $batch_options), 1,
	'all samples reuse the primed active-job snapshot');
is(sampleLockActiveJobs($batch_lock_b, $batch_options), 0,
	'a terminal accounted job no longer keeps its sample lock active');
is($batch_queue_calls, 1, 'the pass performs one squeue query');
is($batch_accounting_calls, 1, 'the pass performs one batched sacct query');
my ($batch_reconciled, $batch_error) = reconcileSlurmDependencies(
	'902;903', 0, $batch_options,
);
is($batch_reconciled, '', 'cached successful and failed terminal dependencies are removed');
like($batch_error, qr/903 \(FAILED, exit 1:0\)/,
	'the cached failed prerequisite is still reported');
is($batch_accounting_calls, 1,
	'dependency reconciliation reuses the pass accounting cache');

my $mfnext_lock = File::Spec->catfile($root, 'mfnext.lock');
my $submitted_before_mfnext = $lock_options->{submittedJobs} || 0;
my @mfnext_dependencies = ('run705', '706');
MFnext($mfnext_lock, \@mfnext_dependencies, 42, $lock_options);
open my $mfnext_fh, '<', $mfnext_lock or die "Cannot read $mfnext_lock: $!";
my $mfnext_text = do { local $/; <$mfnext_fh> };
close $mfnext_fh;
is($mfnext_text, "705\n706\n",
	'MFnext records all sample dependencies in the lock ledger');
is($lock_options->{submittedJobs} || 0, $submitted_before_mfnext,
	'MFnext does not submit a scheduler job to release the lock');

my $accounting_options = slurm_options();
$accounting_options->{rTag} = 'abc';
$accounting_options->{jobAccountingRunner} = sub {
	return (join("\n",
		"801|abc_mA1492|OUT_OF_MEMORY|0:125|OutOfMemory",
		"802|abc_mA1493|TIMEOUT|0:0|TimeLimit",
		"803|abc_GP1492|FAILED|2:0|NonZeroExitCode",
		"804|abc_GP1493|COMPLETED|0:0|None",
		"805|abc_MAP1492|PENDING|0:0|Dependency",
	)."\n", 0);
};
my $failure_summary = slurmJobFailureSummary({
	801 => { requested_name => 'mA1492' },
	802 => { requested_name => 'mA1493' },
	803 => { requested_name => '_GP1492' },
	804 => { requested_name => '_GP1493' },
	805 => { requested_name => '_MAP1492' },
}, $accounting_options);
is($failure_summary->{failed}, 3,
	'Slurm accounting counts only terminal unsuccessful jobs');
is_deeply($failure_summary->{categories}{mA}{failures},
	{ OOM => 1, TIMEOUT => 1 },
	'OOM and timeout occurrences are grouped under the assembly category');
is_deeply($failure_summary->{categories}{_GP}{failures},
	{ FAILED => 1 },
	'non-zero exits are grouped under their MATAFILER job category');
ok(!exists $failure_summary->{categories}{_MAP},
	'pending dependency jobs are not reported as failures');

my $large_accounting_calls = 0;
my $large_accounting_options = slurm_options();
$large_accounting_options->{jobAccountingRunner} = sub {
	$large_accounting_calls++;
	return ("", 0);
};
my %moderately_large_job_set = map {
	(10_000 + $_) => { requested_name => "_GP$_" }
} 1 .. 1_200;
slurmJobFailureSummary(\%moderately_large_job_set, $large_accounting_options);
is($large_accounting_calls, 1,
	'a moderately large run is summarized with one sacct call');

$large_accounting_calls = 0;
my %very_large_job_set = map {
	(20_000 + $_) => { requested_name => "_MAP$_" }
} 1 .. 4_100;
slurmJobFailureSummary(\%very_large_job_set, $large_accounting_options);
is($large_accounting_calls, 2,
	'exceptionally large accounting requests are split into bounded sacct calls');

my $recent_options = {
	slurmDependencyMinAge => 300,
	slurmDependencySubmittedAt => { 201 => 900 },
};
ok(!Mods::Subm::_slurm_dependencies_need_reconciliation([201], $recent_options, 1000),
	'a recently submitted dependency does not query Slurm accounting');
ok(Mods::Subm::_slurm_dependencies_need_reconciliation([201], $recent_options, 1200),
	'a dependency is reconciled once it reaches the Slurm retention age');
ok(Mods::Subm::_slurm_dependencies_need_reconciliation([202], $recent_options, 1000),
	'an externally supplied dependency with unknown age is reconciled');

my $dependency_checks = {
	slurmDependencyAccountingLookup => sub {
		return {
			101 => { state => 'COMPLETED', exit_code => '0:0' },
			102 => { state => 'RUNNING', exit_code => '0:0' },
		};
	},
	slurmDependencyKnownCheck => sub { return $_[0] == 103; },
};
my ($reconciled, $dependency_error) = reconcileSlurmDependencies(
	'101;102;103', 0, $dependency_checks,
);
is($reconciled, '102;103',
	'successful completed Slurm jobs are omitted while active and newly known jobs remain');
is($dependency_error, '', 'valid reconciled dependencies have no error');

($reconciled, $dependency_error) = reconcileSlurmDependencies(
	'104;105', 0, {
		slurmDependencyAccountingLookup => sub {
			return { 104 => { state => 'FAILED', exit_code => '1:0' } };
		},
		slurmDependencyKnownCheck => sub { return 0; },
	},
);
is($reconciled, '', 'failed and genuinely unknown dependencies are not emitted');
like($dependency_error, qr/104 \(FAILED, exit 1:0\).*105 \(unknown to both sacct and slurmctld\)/,
	'failed and aged-out unknown dependencies produce a specific diagnostic');

($reconciled, $dependency_error) = reconcileSlurmDependencies(
	'106', 1, {
		slurmDependencyAccountingLookup => sub {
			return { 106 => { state => 'CANCELLED', exit_code => '0:15' } };
		},
	},
);
is($reconciled, '', 'a terminal unsuccessful job already fulfils afterany');
is($dependency_error, '', 'afterany accepts an unsuccessful terminal dependency');

($reconciled, $dependency_error) = reconcileSlurmDependencies(
	'301;302', 0, {
		slurmDependencyRequireController => 1,
		slurmDependencyAccountingLookup => sub {
			return {
				301 => { state => 'COMPLETED', exit_code => '0:0' },
				302 => { state => 'RUNNING', exit_code => '0:0' },
			};
		},
		slurmDependencyKnownCheck => sub { return $_[0] == 302; },
	},
);
is($reconciled, '302',
	'strict recovery removes a fulfilled aged-out dependency and retains a controller-known job');
is($dependency_error, '', 'strict recovery accepts the verified dependency set');

my $retry_script = File::Spec->catfile($root, 'retry-dependency.sh');
open my $retry_fh, '>', $retry_script or die $!;
print {$retry_fh} "#!/bin/bash\n#SBATCH --dependency=afterok:401:402\necho recovered\n";
close $retry_fh;
my $submission_attempts = 0;
my $retry_options = {
	slurmSubmissionRunner => sub {
		$submission_attempts++;
		return $submission_attempts == 1
			? ("sbatch: error: Batch job submission failed: Job dependency problem\n", 256)
			: ("Submitted batch job 999\n", 0);
	},
	slurmDependencyAccountingLookup => sub {
		return {
			401 => { state => 'COMPLETED', exit_code => '0:0' },
			402 => { state => 'RUNNING', exit_code => '0:0' },
		};
	},
	slurmDependencyKnownCheck => sub { return $_[0] == 402; },
};
my ($retry_output, $retry_status, $did_retry);
{
	local $SIG{__WARN__} = sub { };
	($retry_output, $retry_status, $did_retry) =
		submitSlurmWithDependencyRecovery('sbatch retry-dependency.sh', $retry_script, $retry_options);
}
is($retry_status, 0, 'dependency-problem recovery returns the successful retry status');
is($retry_output, "Submitted batch job 999\n", 'dependency-problem recovery returns retry output');
is($did_retry, 1, 'a Slurm dependency problem triggers exactly one controlled retry');
is($submission_attempts, 2, 'the recovered Slurm script is submitted twice in total');
open $retry_fh, '<', $retry_script or die $!;
my $retry_contents = do { local $/; <$retry_fh> };
close $retry_fh;
like($retry_contents, qr/^#SBATCH --dependency=afterok:402$/m,
	'the retry script retains only controller-present dependencies');

my $fulfilled_script = File::Spec->catfile($root, 'fulfilled-dependencies.sh');
open my $fulfilled_fh, '>', $fulfilled_script or die $!;
print {$fulfilled_fh} "#!/bin/bash\n#SBATCH --dependency=afterok:451:452\necho independent\n";
close $fulfilled_fh;
my $fulfilled_attempts = 0;
{
	local $SIG{__WARN__} = sub { };
	my (undef, $fulfilled_status, $fulfilled_retry) =
		submitSlurmWithDependencyRecovery(
			'sbatch fulfilled-dependencies.sh', $fulfilled_script, {
				slurmSubmissionRunner => sub {
					$fulfilled_attempts++;
					return $fulfilled_attempts == 1
						? ("sbatch: error: Job dependency problem\n", 256)
						: ("Submitted batch job 1000\n", 0);
				},
				slurmDependencyAccountingLookup => sub {
					return {
						451 => { state => 'COMPLETED', exit_code => '0:0' },
						452 => { state => 'COMPLETED', exit_code => '0:0' },
					};
				},
				slurmDependencyKnownCheck => sub { return 0; },
			},
		);
	is($fulfilled_status, 0, 'a job with only fulfilled aged-out dependencies is resubmitted');
	is($fulfilled_retry, 1, 'fulfilled aged-out dependencies trigger one retry');
}
open $fulfilled_fh, '<', $fulfilled_script or die $!;
my $fulfilled_contents = do { local $/; <$fulfilled_fh> };
close $fulfilled_fh;
unlike($fulfilled_contents, qr/^#SBATCH --dependency=/m,
	'the dependency directive is removed when every prerequisite is already fulfilled');

my $failed_retry_script = File::Spec->catfile($root, 'failed-retry-dependency.sh');
open my $failed_retry_fh, '>', $failed_retry_script or die $!;
print {$failed_retry_fh} "#!/bin/bash\n#SBATCH --dependency=afterok:501\necho blocked\n";
close $failed_retry_fh;
my $failed_attempts = 0;
my ($failed_output, $failed_status, $failed_did_retry) =
	submitSlurmWithDependencyRecovery(
		'sbatch failed-retry-dependency.sh', $failed_retry_script, {
			slurmSubmissionRunner => sub {
				$failed_attempts++;
				return ("sbatch: error: Job dependency problem\n", 256);
			},
			slurmDependencyAccountingLookup => sub {
				return { 501 => { state => 'FAILED', exit_code => '1:0' } };
			},
			slurmDependencyKnownCheck => sub { return 0; },
		},
	);
is($failed_status, 256, 'an unsuccessful prerequisite preserves the original submission failure');
is($failed_did_retry, 0, 'an unsuccessful prerequisite is never removed for a retry');
is($failed_attempts, 1, 'unsafe dependency recovery does not resubmit');
like($failed_output, qr/recovery aborted:.*501 \(FAILED, exit 1:0\)/s,
	'the failed prerequisite is identified in the recovery diagnostic');

my $script = File::Spec->catfile($root, 'short-dependency.sh');
my $options = slurm_options();
my ($job, $command) = qsubSystem(
	$script, 'echo ready', 1, '1G', 'consumer', '7;7;;', '', 1, [], $options,
);
open my $fh, '<', $script or die "Cannot read $script: $!";
my $contents = do { local $/; <$fh> };
close $fh;
like($contents, qr/^#SBATCH --dependency=afterok:7$/m,
	'a valid short Slurm dependency is emitted and duplicate entries are removed');
is($job, 'runconsumer', 'dry-run submission keeps its predictable tagged job name');
like($command, qr/^sbatch\s+\Q$script\E/m, 'dry-run returns the scheduler command');

my $empty_script = File::Spec->catfile($root, 'no-dependency.sh');
$options = slurm_options();
qsubSystem($empty_script, 'echo ready', 1, '1G', 'independent', ';;;', '', 1, [], $options);
open $fh, '<', $empty_script or die "Cannot read $empty_script: $!";
$contents = do { local $/; <$fh> };
close $fh;
unlike($contents, qr/^#SBATCH --dependency=/m,
	'an empty dependency collection does not create a scheduler dependency directive');

my $continued_script = File::Spec->catfile($root, 'continued-after-failure.sh');
$options = slurm_options();
$options->{doSubmit} = 1;
$options->{continueOnSubmitError} = 1;
$options->{submissionErrors} = [];
my ($continued_job, $continued_command) = qsubSystem(
	$continued_script, 'echo unsafe', 1, '1G', 'blocked',
	'__MF4_SUBMISSION_FAILED__', '', 1, [], $options,
);
is($continued_job, '__MF4_SUBMISSION_FAILED__',
	'a failed dependency is propagated without submitting its consumer');
is(scalar @{$options->{submissionErrors}}, 1,
	'the skipped dependent submission is recorded for the final summary');

$options = slurm_options();
is(qsubSystemJobAlive([], $options), undef,
	'waiting on an empty dependency set returns without querying the scheduler');

$options = slurm_options();
$options->{doSubmit} = 1;
$options->{submittedJobs} = 0;
$options->{pendingJobCheckInterval} = 60;
$options->{schedulerClock} = sub { return 1000; };
my $pending_queries = 0;
$options->{pendingJobRunner} = sub {
	$pending_queries++;
	return ($pending_queries == 1
		? "101\n102\n103\n"
		: "101\n102\n103\n104\n", 0);
};
qsubSystemWaitMaxJobs(10, 0, $options);
qsubSystemWaitMaxJobs(10, 0, $options);
is($pending_queries, 1,
	'repeated per-sample throttling reuses a recent scheduler count');
$options->{submittedJobs} = 6;
qsubSystemWaitMaxJobs(10, 0, $options);
is($pending_queries, 1,
	'locally submitted jobs are conservatively counted without a scheduler round trip');
$options->{submittedJobs} = 8;
qsubSystemWaitMaxJobs(10, 0, $options);
is($pending_queries, 2,
	'the scheduler is queried immediately when the conservative estimate exceeds the limit');

my $capacity_options = slurm_options();
$capacity_options->{doSubmit} = 1;
$capacity_options->{submittedJobs} = 0;
$capacity_options->{pendingJobCheckInterval} = 0;
$capacity_options->{nonblockingMaxConcurrentJobs} = 1;
my $capacity_queries = 0;
$capacity_options->{liveJobRunner} = sub {
	$capacity_queries++;
	return ("201\n202\n", 0);
};
my ($capacity_admitted, $capacity_output);
{
	local *STDOUT;
	open STDOUT, '>', \$capacity_output
		or die "Cannot capture capacity output: $!";
	$capacity_admitted = qsubSystemWaitMaxJobs(2, 0, $capacity_options);
}
is($capacity_admitted, 0,
	'loop-mode throttling yields immediately instead of waiting at the live-job cap');
ok($capacity_options->{capacityDeferred},
	'the non-blocking throttle records a pass-level capacity deferral');
like($capacity_output, qr/deferring further submissions until the next loop pass/,
	'the capacity deferral is reported explicitly');
is(qsubSystemWaitMaxJobs(2, 0, $capacity_options), 0,
	'later submissions in the same pass remain deferred');
is($capacity_queries, 1,
	'a capacity-deferred pass does not repeat scheduler queries for every sample');

my $capacity_script = File::Spec->catfile($root, 'capacity-deferred.sh');
my $capacity_submission_attempts = 0;
$capacity_options->{maxConcurrentJobs} = 2;
$capacity_options->{slurmSubmissionRunner} = sub {
	$capacity_submission_attempts++;
	return ("Submitted batch job 9999\n", 0);
};
my ($capacity_job) = qsubSystem(
	$capacity_script, 'echo later', 1, '1G', 'capacity', '', '', 1, [],
	$capacity_options,
);
is($capacity_job, deferredSubmissionDependency(),
	'central submission propagates the capacity-deferral dependency marker');
is($capacity_submission_attempts, 0,
	'capacity-deferred work is not passed to sbatch');
ok(submissionDependencyDeferred("17;".deferredSubmissionDependency()),
	'the capacity marker remains detectable when combined with real dependencies');

my $continued_failure = { continueOnSubmitError => 1, submissionErrors => [] };

$options = slurm_options();
my $live_command = '';
$options->{liveJobRunner} = sub {
	$live_command = $_[0];
	return ("11\n12\n13\n14\n", 0);
};
is(numLiveUserJobs($options), 4, 'live-job accounting counts running and pending queue entries together');
is(numLiveUserJobs($options, 0, [qw(11 13 99)]), 2,
	'live-job accounting can be restricted to IDs submitted by this invocation');
unlike($live_command, qr/-t\s+PENDING/, 'the maxConcurrentJobs query is not restricted to pending jobs');

$options = slurm_options();
$options->{jobStatusRunner} = sub {
	return ("101|RUNNING\n102|PENDING\n103|PENDING\n", 0);
};
my $threshold_output = '';
{
	local *STDOUT;
	open STDOUT, '>', \$threshold_output
		or die "Cannot capture threshold output: $!";
	qsubSystemJobAlive([qw(101 102 103)], $options, 0, 3);
}
like($threshold_output,
	qr/1 active job\(s\) remain among 3 queued dependencies; loop threshold 3 reached/,
	'dependency-pending jobs do not inflate the executing-job threshold');

$options = slurm_options();
$options->{activeJobRunner} = sub {
	return ("101\n777\n", 0);
};
is(numActiveUserJobs($options, 0, [qw(101 102)]), 1,
	'active-job counting ignores unrelated user jobs not submitted by this loop');

my $bash_script = File::Spec->catfile($root, 'counted-submission.sh');
$options = bash_options();
qsubSystem($bash_script, 'echo counted', 1, '1G', 'counted', '', '', 1, [], $options);
print "\n"; # qsubSystem's progress prefix intentionally has no trailing newline
is($options->{submittedJobs}, 1,
	'a successful immediate execution increments the submitted-job counter');

my $deferred_script = File::Spec->catfile($root, 'deferred-submission.sh');
qsubSystem($deferred_script, 'echo deferred', 1, '1G', 'deferred', '', '', 0, [], $options);
is($options->{submittedJobs}, 1,
	'a deferred command does not increment the submitted-job counter');

my $release_lock = File::Spec->catfile($root, 'release-sample.lock');
open my $release_fh, '>', $release_lock or die $!;
close $release_fh or die $!;
$options = bash_options();
@{$options}{qw(afterAny useShortQueue tmpSpace)} = (1, 1, '9G');
my @release_dependencies = ('upstream');
MFnext($release_lock, \@release_dependencies, 7, $options);
is_deeply(
	[@{$options}{qw(afterAny useShortQueue tmpSpace)}], [1, 1, '9G'],
	'lock release restores the caller scheduler options after submission');

done_testing;
