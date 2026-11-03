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

my ($accounting_retry_now, $accounting_retry_calls, $accounting_retry_sleeps) = (0, 0, 0);
my $accounting_retry_options = slurm_options();
$accounting_retry_options->{schedulerClock} = sub { $accounting_retry_now };
$accounting_retry_options->{schedulerSleeper} = sub {
	$accounting_retry_sleeps++;
	$accounting_retry_now += $_[0];
};
$accounting_retry_options->{slurmQueryRetrySeconds} = 300;
$accounting_retry_options->{slurmQueryMaxErrorSeconds} = 1_200;
$accounting_retry_options->{jobAccountingRunner} = sub {
	$accounting_retry_calls++;
	return $accounting_retry_calls <= 2
		? ("slurm_load_jobs error: Unexpected message received\n", 1)
		: ("806|abc_FT1|COMPLETED|0:0|None\n", 0);
};
my $accounting_retry_warning = '';
my $retried_failure_summary;
{
	local $SIG{__WARN__} = sub { $accounting_retry_warning .= shift };
	$retried_failure_summary = slurmJobFailureSummary({ 806 => 'FT1' }, $accounting_retry_options);
}
is($accounting_retry_calls, 3, 'Slurm failure accounting retries transient sacct errors');
is($accounting_retry_sleeps, 2, 'Slurm failure accounting waits between transient errors');
is($retried_failure_summary->{failed}, 0, 'successful accounting recovery preserves job results');

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

my $capture_script = File::Spec->catfile($root, 'slurm-stderr.pl');
open my $capture_fh, '>', $capture_script or die $!;
print {$capture_fh} "print STDERR qq{captured scheduler error\\n}; exit 7;\\n";
close $capture_fh or die $!;
my ($captured_error, $captured_status) = Mods::Subm::_run_slurm_submission(
	"$^X $capture_script\n", {});
is($captured_error, "captured scheduler error\n",
	'a newline-terminated scheduler command still captures stderr');
is($captured_status >> 8, 7,
	'scheduler stderr capture preserves the real non-zero exit status');

my ($transient_attempts, @transient_sleeps) = (0);
my @transient_warnings;
my ($transient_output, $transient_status, $transient_retried);
{
	local $SIG{__WARN__} = sub { push @transient_warnings, @_ };
	($transient_output, $transient_status, $transient_retried) =
		submitSlurmWithDependencyRecovery(
			"sbatch transient.sh\n", $capture_script, {
				slurmSubmitMaxAttempts => 3,
				slurmSubmitRetrySeconds => 2,
				slurmSubmitRetryMaxSeconds => 3,
				schedulerSleeper => sub { push @transient_sleeps, $_[0] },
				slurmSubmissionRunner => sub {
					$transient_attempts++;
					return $transient_attempts == 1
						? ("sbatch: error: Batch job submission failed: Unexpected message received\n", 256)
						: ("Submitted batch job 8111\n", 0);
				},
			},
		);
}
is($transient_status, 0, 'a transient Slurm transport response is retried successfully');
is($transient_output, "Submitted batch job 8111\n",
	'transient recovery returns the accepted submission response');
is($transient_retried, 1, 'transient submission recovery reports that it retried');
is($transient_attempts, 2, 'the exact Unexpected-message failure causes one resubmission');
is_deeply(\@transient_sleeps, [2], 'transient submission retry uses bounded backoff');
like(join('', @transient_warnings), qr/Unexpected message received.*retrying in 2s/s,
	'transient retry remains visible without terminating the controller');

my $permanent_attempts = 0;
my ($permanent_output, $permanent_status, $permanent_retried) =
	submitSlurmWithDependencyRecovery(
		'sbatch permanent.sh', $capture_script, {
			slurmSubmitMaxAttempts => 8,
			schedulerSleeper => sub { die 'permanent errors must not sleep' },
			slurmSubmissionRunner => sub {
				$permanent_attempts++;
				return ("sbatch: error: Invalid account or account/partition combination specified\n", 256);
			},
		},
	);
is($permanent_status, 256, 'a permanent Slurm configuration error remains a failure');
is($permanent_retried, 0, 'a permanent Slurm error is not retried');
is($permanent_attempts, 1, 'permanent submission errors make exactly one attempt');
like($permanent_output, qr/Invalid account/, 'permanent failure keeps its diagnostic');

my $warning_script = File::Spec->catfile($root, 'warning-response.sh');
my $warning_options = slurm_options();
$warning_options->{doSubmit} = 1;
$warning_options->{submittedJobs} = 0;
$warning_options->{slurmSubmissionRunner} = sub {
	return ("sbatch: warning: using default cluster\nSubmitted batch job 8123\n", 0);
};
my ($warning_job) = qsubSystem(
	$warning_script, 'echo accepted', 1, '1G', 'warningParse',
	'', '', 1, [], $warning_options,
);
is($warning_job, 'run8123',
	'a valid Slurm job-ID line is parsed despite surrounding warning text');

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
$options->{capacityCheckEverySubmissions} = 10;
my $pending_queries = 0;
$options->{pendingJobRunner} = sub {
	$pending_queries++;
	return ($pending_queries == 1
		? "101\n102\n103\n"
		: "101\n102\n103\n104\n", 0);
};
qsubSystemWaitMaxJobs(100, 0, $options);
qsubSystemWaitMaxJobs(100, 0, $options);
is($pending_queries, 1,
	'repeated per-submission throttling reuses the cached scheduler count');
$options->{submittedJobs} = 9;
qsubSystemWaitMaxJobs(100, 0, $options);
is($pending_queries, 1,
	'locally submitted jobs are conservatively counted without a scheduler round trip');
$options->{submittedJobs} = 10;
qsubSystemWaitMaxJobs(100, 0, $options);
is($pending_queries, 2,
	'the scheduler count is refreshed after the configured submission batch');

my $near_cap_options = slurm_options();
$near_cap_options->{doSubmit} = 1;
$near_cap_options->{submittedJobs} = 0;
$near_cap_options->{capacityCheckEverySubmissions} = 10;
my $near_cap_queries = 0;
$near_cap_options->{pendingJobRunner} = sub {
	$near_cap_queries++;
	return ("101\n102\n103\n", 0);
};
qsubSystemWaitMaxJobs(10, 0, $near_cap_options);
$near_cap_options->{submittedJobs} = 6;
qsubSystemWaitMaxJobs(10, 0, $near_cap_options);
is($near_cap_queries, 1,
	'the conservative count avoids a query while capacity remains below the cap');
$near_cap_options->{submittedJobs} = 7;
qsubSystemWaitMaxJobs(10, 0, $near_cap_options);
is($near_cap_queries, 2,
	'the scheduler is queried before the conservative count can exceed the cap');

my $capacity_options = slurm_options();
$capacity_options->{doSubmit} = 1;
$capacity_options->{submittedJobs} = 0;
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

my ($slurm_retry_now, $slurm_retry_calls, $slurm_retry_sleeps) = (0, 0, 0);
$options = slurm_options();
$options->{schedulerClock} = sub { $slurm_retry_now };
$options->{schedulerSleeper} = sub {
	$slurm_retry_sleeps++;
	$slurm_retry_now += $_[0];
};
$options->{slurmQueryRetrySeconds} = 300;
$options->{slurmQueryMaxErrorSeconds} = 1_200;
$options->{liveJobRunner} = sub {
	$slurm_retry_calls++;
	return $slurm_retry_calls <= 2
		? ("slurm_load_jobs error: Unexpected message received\n", 1)
		: ("501|RUNNING\n", 0);
};
my $slurm_retry_warning = '';
{
	local $SIG{__WARN__} = sub { $slurm_retry_warning .= shift };
	is(numLiveUserJobs($options), 1,
		'transient Slurm controller failures are retried without treating jobs as lost');
}
is($slurm_retry_sleeps, 2, 'transient Slurm failures wait one full retry interval each time');
like($slurm_retry_warning, qr/Transient Slurm scheduler query failure.*retrying in 300s/s,
	'transient Slurm diagnostics identify the retry delay');

my ($slurm_timeout_now, $slurm_timeout_sleeps) = (0, 0);
$options = slurm_options();
$options->{schedulerClock} = sub { $slurm_timeout_now };
$options->{schedulerSleeper} = sub {
	$slurm_timeout_sleeps++;
	$slurm_timeout_now += $_[0];
};
$options->{slurmQueryRetrySeconds} = 300;
$options->{slurmQueryMaxErrorSeconds} = 1_200;
$options->{liveJobRunner} = sub {
	return ("slurm_load_jobs error: Unexpected message received\n", 1);
};
my $slurm_timeout_error = '';
{
	local $SIG{__WARN__} = sub { };
	eval { numLiveUserJobs($options); 1 } or $slurm_timeout_error = $@;
}
like($slurm_timeout_error, qr/failed continuously for 1200s.*Unexpected message received/s,
	'continuous Slurm controller failures stop only after the 20-minute error budget');
is($slurm_timeout_sleeps, 4, 'the 20-minute error budget has four five-minute waits');
my ($shared_queue_calls, $wait_queue_calls) = (0, 0);
$options = slurm_options();
$options->{schedulerClock} = sub { 1000 };
$options->{liveJobRunner} = sub {
	$shared_queue_calls++;
	return ("21|RUNNING\n22|PENDING\n23|PENDING\n", 0);
};
$options->{jobStatusRunner} = sub {
	$wait_queue_calls++;
	return ("21|RUNNING\n22|PENDING\n23|PENDING\n", 0);
};
is(numLiveUserJobs($options), 3,
	'a live-job count records a reusable full queue snapshot');
{
	local *STDOUT;
	open STDOUT, '>', \my $shared_wait_output
		or die "Cannot capture shared-snapshot output: $!";
	qsubSystemJobAlive([qw(21 22 23)], $options, 0, 1);
}
is($shared_queue_calls, 1,
	'the live-job count performs one full scheduler query');
is($wait_queue_calls, 0,
	'the immediately following dependency wait reuses that scheduler snapshot');


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

# A bounded wait lets the caller escalate accounting-confirmed OOM outcomes while
# the rest of a submission wave is still queued, instead of blocking on its tail.
$options = slurm_options();
my $bounded_polls = 0;
$options->{jobPollSeconds} = 1;
$options->{jobStatusRunner} = sub {
	$bounded_polls++;
	return ("101|RUNNING\n102|PENDING\n", 0);
};
my ($bounded_remaining, $bounded_output);
{
	local *STDOUT;
	open STDOUT, '>', \$bounded_output
		or die "Cannot capture bounded-wait output: $!";
	$bounded_remaining = qsubSystemJobAlive([qw(101 102 103)], $options, 0, -1, 0);
}
is_deeply([sort @{$bounded_remaining}], [qw(101 102)],
	'a spent wait budget returns the dependencies that are still queued');
is($bounded_polls, 1,
	'the budget is checked after the poll, so one scheduler query is enough');
like($bounded_output, qr/2\/3 job\(s\) still queued after \d+s; handing control back/,
	'a bounded wait reports what it hands back to the caller');

$options = slurm_options();
$options->{jobPollSeconds} = 1;
$options->{jobStatusRunner} = sub { return ("", 0); };
is_deeply(qsubSystemJobAlive([qw(101 102)], $options, 0, -1, 3600), [],
	'a wave that has fully drained returns an empty set well inside its budget');

$options = slurm_options();
$options->{activeJobRunner} = sub {
	return ("101\n777\n", 0);
};
is(numActiveUserJobs($options, 0, [qw(101 102)]), 1,
	'active-job counting ignores unrelated user jobs not submitted by this loop');

# Slurm orders pending jobs by priority, so a recovery job submitted after a
# large wave is simply the youngest and starts last. Handicapping the wave with
# --nice, and submitting recovery at nice 0, is the only user-level lever that
# reorders jobs the scheduler has already accepted.
my $niced_script = File::Spec->catfile($root, 'niced-submission.sh');
$options = slurm_options();
$options->{jobNice} = 5000;
qsubSystem($niced_script, 'echo niced', 1, '1G', 'niced', '', '', 1, [], $options);
print "\n"; # qsubSystem's progress prefix intentionally has no trailing newline
open my $niced_fh, '<', $niced_script or die "Cannot read $niced_script: $!";
my $niced_text = do { local $/; <$niced_fh> };
close $niced_fh;
like($niced_text, qr/^#SBATCH --nice=5000$/m,
	'a bulk submission carries its configured priority handicap');

my $recovery_script = File::Spec->catfile($root, 'recovery-submission.sh');
$options->{jobNice} = 0;
qsubSystem($recovery_script, 'echo recovery', 1, '1G', 'recovery', '', '', 1, [], $options);
print "\n";
open my $recovery_fh, '<', $recovery_script or die "Cannot read $recovery_script: $!";
my $recovery_text = do { local $/; <$recovery_fh> };
close $recovery_fh;
unlike($recovery_text, qr/--nice=/,
	'a recovery submission asks for no handicap at all, outranking the niced wave');

my $negative_script = File::Spec->catfile($root, 'negative-nice-submission.sh');
$options->{jobNice} = -100;
qsubSystem($negative_script, 'echo negative', 1, '1G', 'negative', '', '', 1, [], $options);
print "\n";
open my $negative_fh, '<', $negative_script or die "Cannot read $negative_script: $!";
my $negative_text = do { local $/; <$negative_fh> };
close $negative_fh;
unlike($negative_text, qr/--nice=/,
	'a negative nice is dropped rather than submitted, since raising priority needs an operator');

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
