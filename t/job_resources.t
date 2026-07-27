use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::JobResources qw(initialize_job_resource_log timed_job_command);

my $tmp = tempdir(CLEANUP => 1);
my $resource_log = File::Spec->catfile($tmp, 'job_resources.tsv');
my $result = File::Spec->catfile($tmp, 'result.txt');

initialize_job_resource_log($resource_log);
initialize_job_resource_log($resource_log);
is(
	do {
		open my $input, '<', $resource_log or die $!;
		local $/;
		<$input>;
	},
	"stage\tjob_id\twall_seconds\tmax_rss_kb\n",
	'resource-log initialization is idempotent',
);

my $success_command = timed_job_command(
	'unit_stage', "printf 'completed\\n' > '$result'", $resource_log,
);
is(system('bash', '-c', $success_command), 0,
	'timed job wrapper preserves successful execution');
is(
	do {
		open my $input, '<', $result or die $!;
		local $/;
		<$input>;
	},
	"completed\n",
	'timed job wrapper executes the complete payload',
);

open my $resource_input, '<', $resource_log or die $!;
my @resource_rows = <$resource_input>;
close $resource_input;
is(scalar(@resource_rows), 2, 'one resource row is appended for one completed command');
like(
	$resource_rows[1],
	qr/^unit_stage\tlocal\t\d+(?:\.\d+)?\t(?:\d+|NA)\n$/,
	'resource row records stage, job ID, wall seconds, and peak RSS',
);

my $failure_command = timed_job_command('failing_stage', 'exit 7', $resource_log);
my $failure_status = system('bash', '-c', $failure_command);
is($failure_status >> 8, 7, 'timed job wrapper preserves the payload exit status');

eval { timed_job_command('bad stage', 'true', $resource_log) };
like($@, qr/Invalid resource-accounting stage/,
	'unsafe resource stage labels are rejected');

done_testing();
