use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::Subm qw(
	qsubSystem qsubSystemJobAlive numActiveUserJobs reconcileSlurmDependencies
	submitSlurmWithDependencyRecovery
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

done_testing;
