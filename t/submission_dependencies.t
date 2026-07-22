use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::Subm qw(qsubSystem qsubSystemJobAlive reconcileSlurmDependencies);

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
