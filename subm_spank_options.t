use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::Subm qw(qsubSystem2);

# A job submitted from inside a job that has its own SPANK options (here: --localscratch=ssd:15) must not pass them
# on: some sites apply the inherited value instead of the new job's #SBATCH request.
local $ENV{SLURM_SPANK__SLURM_SPANK_OPTION_tmpvol_localscratch} = 'ssd:15';
local $ENV{_SLURM_SPANK_OPTION_tmpvol_localscratch} = 'ssd:15';
local $ENV{MF4_SPANK_TEST_KEEP} = 'kept';

my ($output, $status) = Mods::Subm::_run_slurm_submission('env', {});
is($status, 0, 'the submission command ran');
unlike($output, qr/SPANK_OPTION_/, 'inherited SPANK plugin options are not passed to sbatch');
like($output, qr/^MF4_SPANK_TEST_KEEP=kept$/m, 'all other environment variables are passed unchanged');
is($ENV{SLURM_SPANK__SLURM_SPANK_OPTION_tmpvol_localscratch}, 'ssd:15', 'the submitting process keeps its own environment');

# qsubSystem2 (resubmits saved scripts) drops them the same way; local mode runs the script directly
my $tmpdir = tempdir(CLEANUP => 1);
my $script = File::Spec->catfile($tmpdir, 'job.sh');
my $envF = File::Spec->catfile($tmpdir, 'env.txt');
open my $fh, '>', $script or die "Cannot write $script: $!";
print {$fh} "#!/bin/bash\nenv > $envF\n";
close $fh;
qsubSystem2($script, { qmode => 'bash', jobNice => 0 });
open my $in, '<', $envF or die "Cannot read $envF: $!";
my $jobEnv = do { local $/; <$in> };
close $in;
unlike($jobEnv, qr/SPANK_OPTION_/, 'qsubSystem2 does not pass inherited SPANK plugin options either');
like($jobEnv, qr/^MF4_SPANK_TEST_KEEP=kept$/m, 'qsubSystem2 passes all other environment variables');

done_testing;
