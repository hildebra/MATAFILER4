use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::Subm qw(qsubSystem);

my $tmpdir = tempdir(CLEANUP => 1);
my $script = File::Spec->catfile($tmpdir, 'download.sh');
$script =~ s{\\}{/}g;
my $options = {
	rTag => 'TST',
	doSubmit => 0,
	doSync => 0,
	medQueue => 'compute',
	medTime => '24:00:00',
	longQueue => 'long',
	longTime => '168:00:00',
	gpuQueue => 'gpu',
	netQueue => 'network',
	highMemQueue => 'highmem',
	shortQueue => 'short',
	useLongQueue => 0,
	useGPUQueue => 0,
	useNetQueue => 1,
	useShortQueue => 0,
	useHiMemQueue => 0,
	gpuCount => 0,
	tmpSpace => 0,
	tmpSpaceTag => '',
	excludeNodes => '',
	submissionConfig => '',
	constraint => [],
	LOCKfile => '',
	xtraNodeCmds => '',
	wcKeysForJob => '',
	qmode => 'slurm',
};

qsubSystem($script, 'echo download', 1, '1G', 'download', '', '', 0, [], $options);

open my $fh, '<', $script or die "Cannot read $script: $!";
my $contents = do { local $/; <$fh> };
close $fh;

like($contents, qr/^#SBATCH -p "network"$/m, 'uses the network queue');
like($contents, qr/^#SBATCH --time=168:00:00$/m, 'uses the long wall time');
like($contents,
	qr/^echo "SLURM job ID: \$SLURM_JOB_ID"\necho \$HOSTNAME;$/m,
	'prints the allocated Slurm job ID before executing the job payload');
is($options->{useNetQueue}, 0, 'network queue selection is one-shot');

my $scratch_script = File::Spec->catfile($tmpdir, 'scratch.sh');
$scratch_script =~ s{\\}{/}g;
$options->{tmpSpace} = '12G';
{
	no warnings 'redefine';
	local *Mods::Subm::getProgPaths = sub {
		my ($key) = @_;
		return '/node/local/tmp' if $key eq 'nodeTmpDir';
		return '';
	};
	qsubSystem($scratch_script, 'pwd', 1, '1G', 'scratch', '', '', 0, [], $options);
}
open my $scratch_fh, '<', $scratch_script
	or die "Cannot read $scratch_script: $!";
my $scratch_contents = do { local $/; <$scratch_fh> };
close $scratch_fh;
like($scratch_contents, qr/^node_tmp_root="\/node\/local\/tmp"$/m,
	'requested node scratch uses the configured node-local root');
like($scratch_contents, qr/^export TMPDIR="\$node_tmp_workdir"\ncd "\$node_tmp_workdir"$/m,
	'node-scratch jobs execute from their job-specific local directory');

done_testing;
