use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::Subm qw(qsubSystem qsubSystemJobAlive);

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

my $root = tempdir(CLEANUP => 1);
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

$options = slurm_options();
is(qsubSystemJobAlive([], $options), undef,
	'waiting on an empty dependency set returns without querying the scheduler');

done_testing;
