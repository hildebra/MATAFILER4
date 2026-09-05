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
	downloadQueue => 'nbi-download',
	highMemQueue => 'highmem',
	shortQueue => 'short',
	useLongQueue => 0,
	useGPUQueue => 0,
	useNetQueue => 1,
	useDownloadQueue => 0,
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
like($contents, qr/^#SBATCH --export=ALL,TMPDIR=\/tmp$/m,
	'Slurm jobs reset inherited TMPDIR before the shell starts');
like($contents, qr/^#SBATCH --chdir=\Q$tmpdir\E$/m,
	'Slurm jobs start beside their persistent generated script, not the submitter cwd');
like($contents,
	qr/^echo "SLURM job ID: \$SLURM_JOB_ID"\necho \$HOSTNAME;$/m,
	'prints the allocated Slurm job ID before executing the job payload');
like($contents, qr/^#SBATCH --nice=2500$/m,
	'a direct Slurm submission without a priority option uses the moderate default');
is($options->{useNetQueue}, 0, 'network queue selection is one-shot');

my $archive_script = File::Spec->catfile($tmpdir, 'archive-download.sh');
$archive_script =~ s{\\}{/}g;
$options->{useDownloadQueue} = 1;
qsubSystem(
	$archive_script, 'echo archive', 2, '16G', 'archive',
	'', '', 0, [], $options,
);
open my $archive_fh, '<', $archive_script
	or die "Cannot read $archive_script: $!";
my $archive_contents = do { local $/; <$archive_fh> };
close $archive_fh;
like($archive_contents, qr/^#SBATCH -p "nbi-download"$/m,
	'archive acquisition uses the dedicated download queue');
like($archive_contents, qr/^#SBATCH --time=24:00:00$/m,
	'download jobs retain the default queue wall time');

{
	no warnings 'redefine';
	my %config = (
		MFLRDir => $tmpdir,
		mediumQueue => 'compute',
		downloadQueue => '',
		qsubPEenv => 'smp',
	);
	local *Mods::Subm::getProgPaths = sub {
		my ($key) = @_;
		return exists($config{$key}) ? $config{$key} : '';
	};
	my $fallback = Mods::Subm::emptyQsubOpt(0, '', 'slurm');
	is($fallback->{downloadQueue}, 'compute',
		'an empty downloadQueue falls back to mediumQueue');
	is($fallback->{jobNice}, 2500,
		'every new scheduler option set uses the moderate priority default');
}
is($options->{useDownloadQueue}, 0,
	'download queue selection is one-shot');

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
like($scratch_contents, qr/^#SBATCH --chdir=\Q$tmpdir\E$/m,
	'node-scratch jobs still start in the persistent script directory before changing locally');
like($scratch_contents, qr/^export TMPDIR="\$node_tmp_workdir"\ncd "\$node_tmp_workdir"$/m,
	'node-scratch jobs execute from their job-specific local directory');

my $legacy_local = '/nbi/local/ssd/old-job/MF4/matafiler4.old-job';
my $legacy_script = File::Spec->catfile($tmpdir, 'legacy-cwd.sh');
$options->{tmpSpace} = 0;
qsubSystem($legacy_script, 'pwd', 1, '1G', 'legacy', '', $legacy_local, 0, [], $options);
open my $legacy_fh, '<', $legacy_script or die "Cannot read $legacy_script: $!";
my $legacy_contents = do { local $/; <$legacy_fh> };
close $legacy_fh;
like($legacy_contents, qr/^#SBATCH --chdir=\Q$tmpdir\E$/m,
	'a stale caller-provided node-local cwd cannot become Slurm\'s startup directory');
unlike($legacy_contents, qr/\Q$legacy_local\E/,
	'the stale node-local cwd is absent from the generated Slurm script');

done_testing;
