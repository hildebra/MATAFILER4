package Mods::JobResources;

use strict;
use warnings;
use Exporter qw(import);

our @EXPORT_OK = qw(initialize_job_resource_log timed_job_command);

sub _shell_quote {
	my ($value) = @_;
	$value = '' unless defined $value;
	$value =~ s/'/'"'"'/g;
	return "'$value'";
}

sub initialize_job_resource_log {
	my ($path) = @_;
	die "Job resource log path is required\n" unless defined($path) && length($path);
	return if -s $path;
	open my $output, '>', $path or die "Cannot create job resource log $path: $!\n";
	print {$output} "stage\tjob_id\twall_seconds\tmax_rss_kb\n"
		or die "Cannot write job resource log $path: $!\n";
	close $output or die "Cannot close job resource log $path: $!\n";
}

sub timed_job_command {
	my ($stage, $command, $resource_log) = @_;
	die "Invalid resource-accounting stage '$stage'\n"
		unless defined($stage) && $stage =~ /^[A-Za-z0-9_.-]+$/;
	die "Job resource log path is required\n"
		unless defined($resource_log) && length($resource_log);

	my $payload = "set -eo pipefail\n" . ($command // '');
	my $quoted_payload = _shell_quote($payload);
	my $quoted_log = _shell_quote($resource_log);
	my $job_id = '${SLURM_JOB_ID:-${JOB_ID:-${LSB_JOBID:-local}}}';
	my $format = "$stage\t$job_id\t%e\t%M";
	return join("\n",
		'if [ -x /usr/bin/time ]; then',
		qq{  /usr/bin/time -f "$format" -a -o $quoted_log bash -c $quoted_payload},
		'else',
		'  resource_start_seconds=$(date +%s)',
		'  set +e',
		qq{  bash -c $quoted_payload},
		'  resource_status=$?',
		'  set -e',
		'  resource_wall_seconds=$(($(date +%s) - resource_start_seconds))',
		qq{  printf '%s\\t%s\\t%s\\t%s\\n' '$stage' "$job_id" "\$resource_wall_seconds" NA >> $quoted_log},
		'  exit $resource_status',
		'fi',
		'',
	);
}

1;
