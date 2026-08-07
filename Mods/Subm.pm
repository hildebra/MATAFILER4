package Mods::Subm;

use warnings;
use strict;
use Fcntl qw(:flock);
#use List::MoreUtils 'first_index'; 
use Mods::IO_Tamoc_progs qw(getProgPaths convert2Gb);
use Mods::WorkflowControl qw(normalise_job_dependencies);


use Exporter qw(import);
our @EXPORT_OK = qw( findQsubSys emptyQsubOpt qsubSystem qsubSystem2 qsubSystemJobAlive
		qsubSystemWaitMaxJobs MFnext add2SampleDeps numUserJobs numLiveUserJobs numActiveUserJobs reconcileSlurmDependencies
		recordSampleLockJobs sampleLockActiveJobs primeSampleLockJobSnapshot slurmJobFailureSummary
		submitSlurmWithDependencyRecovery deferredSubmissionDependency
		submissionDependencyDeferred handleSubmissionFailure);

my $FAILED_SUBMISSION_DEPENDENCY = '__MF4_SUBMISSION_FAILED__';
my $DEFERRED_SUBMISSION_DEPENDENCY = '__MF4_SUBMISSION_DEFERRED__';

sub deferredSubmissionDependency { return $DEFERRED_SUBMISSION_DEPENDENCY; }
sub submissionDependencyDeferred {
	return scalar grep { $_ eq $DEFERRED_SUBMISSION_DEPENDENCY }
		split /;/, normalise_job_dependencies(@_);
}


sub recordSampleLockJobs {
	my ($lock_file, $jobs, $optHR) = @_;
	return 0 unless defined($lock_file) && $lock_file ne '';
	my @job_ids = split /;/, normalise_job_dependencies($jobs);
	my $run_tag = $optHR->{rTag} || '';
	for (@job_ids) {
		s/^\Q$run_tag\E//;
	}
	@job_ids = grep { /^\d+$/ } @job_ids;
	return 0 unless @job_ids;

	open my $lock_fh, '>>', $lock_file
		or die "Cannot update sample job lock $lock_file: $!\n";
	flock($lock_fh, LOCK_EX)
		or die "Cannot lock sample job ledger $lock_file: $!\n";
	print {$lock_fh} "$_\n" for @job_ids;
	close $lock_fh
		or die "Cannot close sample job ledger $lock_file: $!\n";

	# A cached scheduler snapshot predating these submissions must still treat
	# every locally accepted job as active.
	my $state = $optHR->{sampleLockActiveState};
	if ($state) {
		$state->{jobs}{$_} = 1 for @job_ids;
	}
	my $queue_state = $optHR->{schedulerQueueState};
	if ($queue_state) {
		$queue_state->{jobs}{$_} = 1 for @job_ids;
		$queue_state->{states}{$_} ||= 'PENDING' for @job_ids;
	}
	return scalar @job_ids;
}

sub _sample_lock_job_ids {
	my ($lock_file) = @_;
	return undef unless defined($lock_file) && -e $lock_file;
	open my $lock_fh, '<', $lock_file
		or do {
			warn "Cannot inspect sample job lock $lock_file: $!\n";
			return undef;
		};
	flock($lock_fh, LOCK_SH)
		or do {
			warn "Cannot lock sample job ledger $lock_file for reading: $!\n";
			close $lock_fh;
			return undef;
		};
	my (%wanted, $malformed);
	while (my $line = <$lock_fh>) {
		$line =~ s/^\s+|\s+$//g;
		next if $line eq '';
		if ($line =~ /^\d+$/) {
			$wanted{$line} = 1;
		} else {
			$malformed = 1;
		}
	}
	close $lock_fh;
	return undef if $malformed || !%wanted;
	return \%wanted;
}

sub _run_scheduler_query {
	my ($runner, $command, @extra) = @_;
	if ($runner) {
		my ($output, $status) = $runner->($command, @extra);
		$status ||= 0;
		return ($output, $status);
	}
	my $output = `$command`;
	return ($output, $?);
}

sub _scheduler_queue_command {
	my ($qmode) = @_;
	return "squeue -h -u \$USER -o '%i|%T'" if $qmode eq 'slurm';
	return "qstat -u \$USER | awk 'NR > 2 {print \$1 \"|\" \$5}'"
		if $qmode eq 'sge';
	return "bjobs -noheader -o 'jobid stat'" if $qmode eq 'lsf';
	return '' if $qmode eq 'bash';
	die "Cannot inspect jobs for scheduler mode '$qmode'\n";
}

sub _scheduler_state_is_executing {
	my ($qmode, $state) = @_;
	$state = uc($state || '');
	return $state =~ /^(?:RUNNING|COMPLETING|CONFIGURING|RESIZING|SIGNALING|STAGE_OUT)$/
		if $qmode eq 'slurm';
	return $state =~ /^(?:R|T|RR|RT)$/ if $qmode eq 'sge';
	return $state eq 'RUN' if $qmode eq 'lsf';
	return 0;
}

sub _scheduler_queue_snapshot {
	my ($optHR, %args) = @_;
	$optHR ||= {};
	my $qmode = $optHR->{qmode} || 'slurm';
	return {
		checkedAt => time, jobs => {}, states => {}, qmode => $qmode,
	} if $qmode eq 'bash';

	my $clock = $optHR->{schedulerClock};
	my $now = $clock ? $clock->() : time;
	my $cached = $optHR->{schedulerQueueState};
	my $reuse_seconds = 0 + ($args{reuse_seconds} || 0);
	if ($cached && $reuse_seconds > 0
			&& $now - ($cached->{checkedAt} || 0) < $reuse_seconds) {
		return $cached;
	}

	my $command = _scheduler_queue_command($qmode);
	my ($output, $status) = _run_scheduler_query($args{runner}, $command);
	die($args{error} || "Failed to query scheduler jobs with: $command\n")
		if $status != 0;
	my (%jobs, %states);
	for my $line (split /\n/, $output) {
		$line =~ s/^\s+|\s+$//g;
		next if $line eq '';
		my ($job_id, $state);
		if ($line =~ /^(\d+)\|(\S+)$/) {
			($job_id, $state) = ($1, uc($2));
		} elsif ($line =~ /^(\d+)\s+(\S+)$/) {
			($job_id, $state) = ($1, uc($2));
		} elsif ($line =~ /^(\d+)$/) {
			($job_id, $state) = ($1, $args{bare_state} || '');
		} else {
			next;
		}
		$jobs{$job_id} = 1;
		$states{$job_id} = $state;
	}
	my $snapshot = {
		checkedAt => $clock ? $clock->() : time,
		jobs => \%jobs,
		states => \%states,
		qmode => $qmode,
	};
	$optHR->{schedulerQueueState} = $snapshot;
	return $snapshot;
}

sub _slurm_job_batches {
	my ($job_ids) = @_;
	my @remaining = grep { /^\d+$/ } @{$job_ids || []};
	my @batches;
	while (@remaining) {
		my (@batch, $batch_chars);
		while (@remaining && @batch < 4000) {
			my $next_chars = length($remaining[0]) + (@batch ? 1 : 0);
			last if @batch && $batch_chars + $next_chars > 48_000;
			my $job_id = shift @remaining;
			push @batch, $job_id;
			$batch_chars += $next_chars;
		}
		push @batches, \@batch;
	}
	return @batches;
}

sub _batched_slurm_accounting_lookup {
	my ($job_ids, $optHR) = @_;
	my %records;
	my $runner = $optHR->{sampleLockAccountingRunner}
		|| $optHR->{slurmDependencyAccountingRunner};
	for my $batch (_slurm_job_batches($job_ids)) {
		my $command = "sacct -X -n -P -j ".join(',', @{$batch})
			." --format=JobIDRaw,State,ExitCode";
		my ($output, $status) = _run_scheduler_query($runner, $command, $batch);
		return undef if $status != 0;
		for my $line (split /\n/, $output) {
			my ($job_id, $state, $exit_code) = split /\|/, $line, 3;
			next unless defined($job_id) && $job_id =~ /^\d+$/;
			for ($state, $exit_code) {
				$_ = '' unless defined $_;
				s/^\s+|\s+$//g;
			}
			$state =~ s/\+$//;
			$state =~ s/\s.*$//;
			$records{$job_id} = { state => uc($state), exit_code => $exit_code };
		}
		$records{$_} ||= { state => '', exit_code => '' }
			for @{$batch};
	}
	return \%records;
}

sub primeSampleLockJobSnapshot {
	my ($lock_files, $optHR) = @_;
	$optHR ||= {};
	return { tracked => 0, queued => 0, accounted => 0 }
		unless (($optHR->{qmode} || 'slurm') eq 'slurm');
	my %tracked;
	for my $lock_file (@{$lock_files || []}) {
		my $wanted = _sample_lock_job_ids($lock_file);
		$tracked{$_} = 1 for keys %{$wanted || {}};
	}
	$optHR->{accountingJobIds}{$_} = 1 for keys %tracked;
	my $state = eval {
		_scheduler_queue_snapshot(
			$optHR,
			runner => $optHR->{sampleLockJobRunner},
			error => "Failed to query scheduler while priming sample locks\n",
		);
	};
	if (!$state) {
		warn "Failed to query scheduler while priming sample locks; retaining existing lock state\n";
		return undef;
	}
	$state->{passSnapshot} = 1;
	$optHR->{sampleLockActiveState} = $state;
	my $throttle_live_jobs = scalar keys %{$state->{jobs}};
	$throttle_live_jobs--
		if (($ENV{SLURM_JOBID} || '') ne '' && $state->{jobs}{$ENV{SLURM_JOBID}});
	$optHR->{liveJobThrottleState} = {
		checkedAt => $state->{checkedAt}, liveJobs => $throttle_live_jobs,
		submittedJobs => $optHR->{submittedJobs} || 0,
	};
	$optHR->{slurmDependencyAccountingCache} ||= {};
	$optHR->{slurmDependencyAccountingCache}{$_} =
		{ state => $state->{states}{$_}, exit_code => '' }
		for grep { $tracked{$_} } keys %{$state->{jobs}};
	my @missing = grep { !$state->{jobs}{$_} } sort { $a <=> $b } keys %tracked;
	my $records = @missing ? _batched_slurm_accounting_lookup(\@missing, $optHR) : {};
	if (defined $records) {
		$optHR->{slurmDependencyAccountingCache} ||= {};
		@{$optHR->{slurmDependencyAccountingCache}}{keys %{$records}} = values %{$records};
	}
	return {
		tracked => scalar(keys %tracked),
		queued => scalar(grep { $tracked{$_} } keys %{$state->{jobs}}),
		accounted => defined($records) ? scalar(keys %{$records}) : 0,
	};
}

sub sampleLockActiveJobs {
	my ($lock_file, $optHR) = @_;
	my $wanted = _sample_lock_job_ids($lock_file);
	# An empty or legacy lock has no ownership evidence. Retain it unless the
	# caller explicitly requested lock removal.
	return undef unless $wanted;
	$optHR->{accountingJobIds}{$_} = 1 for keys %{$wanted};

	my $clock = $optHR->{schedulerClock};
	my $now = $clock ? $clock->() : time;
	my $cache_seconds = defined($optHR->{sampleLockCheckInterval})
		? 0 + $optHR->{sampleLockCheckInterval} : 30;
	$cache_seconds = 0 if $cache_seconds < 0;
	my $state = $optHR->{sampleLockActiveState};
	if (!$state || (!$state->{passSnapshot}
			&& $now - $state->{checkedAt} >= $cache_seconds)) {
		if (($optHR->{qmode} || 'slurm') eq 'bash') {
			return 0;
		}
		$state = eval {
			_scheduler_queue_snapshot(
				$optHR,
				runner => $optHR->{sampleLockJobRunner},
				error => "Failed to query scheduler while inspecting $lock_file\n",
			);
		};
		if (!$state) {
			warn "Failed to query scheduler while inspecting $lock_file; retaining lock\n";
			return undef;
		}
		$optHR->{sampleLockActiveState} = $state;
	}
	return scalar grep { $state->{jobs}{$_} } keys %{$wanted};
}

sub _job_category {
	my ($job_name, $requested_name, $run_tag) = @_;
	my $category = defined($requested_name) && $requested_name ne ''
		? $requested_name : ($job_name || 'unknown');
	$category =~ s/^\Q$run_tag\E//
		if defined($run_tag) && $run_tag ne '';
	$category =~ s/\s.*$//;
	$category =~ s/(?:[._-]?\d+)+(?:[._-].*)?$//;
	$category =~ s/[._-]+$//;
	return $category ne '' ? $category : 'unknown';
}

sub _failure_class {
	my ($state, $exit_code, $reason) = @_;
	$state ||= '';
	$exit_code ||= '';
	$reason ||= '';
	my %terminal = map { $_ => 1 } qw(
		BOOT_FAIL CANCELLED COMPLETED DEADLINE FAILED NODE_FAIL
		OUT_OF_MEMORY PREEMPTED REVOKED TIMEOUT
	);
	return '' unless $terminal{$state};
	return 'OOM' if $state eq 'OUT_OF_MEMORY' || $reason =~ /out.?of.?memory/i;
	return 'TIMEOUT' if $state eq 'TIMEOUT' || $reason =~ /time.?limit/i;
	return 'DEPENDENCY' if $reason =~ /dependenc/i;
	return 'NODE_FAIL' if $state =~ /^(?:BOOT_FAIL|NODE_FAIL)$/;
	return 'PREEMPTED' if $state =~ /^(?:PREEMPTED|REVOKED)$/;
	return 'CANCELLED' if $state eq 'CANCELLED';
	return 'FAILED'
		if $state eq 'FAILED' || ($exit_code ne '' && $exit_code !~ /^0:0$/);
	return '';
}

sub slurmJobFailureSummary {
	my ($jobs, $optHR) = @_;
	$jobs ||= {};
	$optHR ||= {};
	return { queried => 0, failed => 0, unresolved => 0, categories => {} }
		unless (($optHR->{qmode} || '') eq 'slurm');

	my %requested;
	if (ref($jobs) eq 'HASH') {
		for my $job_id (keys %{$jobs}) {
			next unless $job_id =~ /^\d+$/;
			my $record = $jobs->{$job_id};
			$requested{$job_id} = ref($record) eq 'HASH'
				? ($record->{requested_name} || '') : ($record || '');
		}
	} elsif (ref($jobs) eq 'ARRAY') {
		$requested{$_} = '' for grep { /^\d+$/ } @{$jobs};
	}
	my @job_ids = sort { $a <=> $b } keys %requested;
	return { queried => 0, failed => 0, unresolved => 0, categories => {} }
		unless @job_ids;

	my (%seen, %categories);
	my $failed = 0;
	for my $batch (_slurm_job_batches(\@job_ids)) {
		my $command = "sacct -X -n -P -j ".join(',', @{$batch})
			." --format=JobIDRaw,JobName,State,ExitCode,Reason";
		my ($output, $status);
		if (my $runner = $optHR->{jobAccountingRunner}) {
			($output, $status) = $runner->($command, $batch);
			$status ||= 0;
		} else {
			$output = `$command 2>/dev/null`;
			$status = $?;
		}
		return undef if $status != 0;
		for my $line (split /\n/, $output) {
			my ($job_id, $job_name, $state, $exit_code, $reason) =
				split /\|/, $line, 5;
			next unless defined($job_id) && exists $requested{$job_id};
			$seen{$job_id} = 1;
			for ($job_name, $state, $exit_code, $reason) {
				$_ = '' unless defined $_;
				s/^\s+|\s+$//g;
			}
			$state =~ s/\+$//;
			$state =~ s/\s.*$//;
			$state = uc($state);
			my $failure = _failure_class($state, $exit_code, $reason);
			next if $failure eq '';
			my $category = _job_category(
				$job_name, $requested{$job_id}, $optHR->{rTag} || '',
			);
			$categories{$category}{total}++;
			$categories{$category}{failures}{$failure}++;
			$failed++;
		}
	}
	return {
		queried => scalar(keys %requested),
		failed => $failed,
		unresolved => scalar(grep { !$seen{$_} } keys %requested),
		categories => \%categories,
	};
}

sub _slurm_dependency_accounting_lookup {
	my ($dependencies, $optHR) = @_;
	$optHR ||= {};
	my @job_ids = grep { /^\d+$/ } @{$dependencies};
	return {} unless (@job_ids);
	my $cache = $optHR->{slurmDependencyAccountingCache} ||= {};
	my @missing = grep { !exists $cache->{$_} } @job_ids;
	if (@missing) {
		my $records = _batched_slurm_accounting_lookup(\@missing, $optHR);
		return undef unless defined $records;
		@{$cache}{keys %{$records}} = values %{$records};
	}
	return {
		map { exists($cache->{$_}) ? ($_ => $cache->{$_}) : () } @job_ids
	};
}

sub _slurm_dependency_is_known {
	my ($job_id) = @_;
	return 1 unless ($job_id =~ /^\d+$/);
	system("scontrol show job $job_id >/dev/null 2>&1");
	return $? == 0;
}

sub _slurm_dependency_min_age {
	my ($optHR) = @_;
	return $optHR->{slurmDependencyMinAge}
		if defined $optHR->{slurmDependencyMinAge};

	my $output = `scontrol show config 2>/dev/null`;
	my $min_age = 300; # Slurm's default MinJobAge, in seconds.
	$min_age = $1 if ($? == 0 && $output =~ /^\s*MinJobAge\s*=\s*(\d+)/m);
	$optHR->{slurmDependencyMinAge} = $min_age;
	return $min_age;
}

sub _slurm_dependencies_need_reconciliation {
	my ($dependencies, $optHR, $now) = @_;
	$now = time unless defined $now;
	my $submitted_at = $optHR->{slurmDependencySubmittedAt} ||= {};
	my $min_age = _slurm_dependency_min_age($optHR);
	for my $dependency (@{$dependencies}) {
		return 1 unless exists $submitted_at->{$dependency};
		return 1 if $now - $submitted_at->{$dependency} >= $min_age;
	}
	return 0;
}

sub reconcileSlurmDependencies {
	my ($dependencies, $after_any, $optHR) = @_;
	$optHR ||= {};
	my @dependencies = split /;/, normalise_job_dependencies($dependencies);
	return ('', '') unless (@dependencies);
	my $require_controller = $optHR->{slurmDependencyRequireController} || 0;

	my $lookup = $optHR->{slurmDependencyAccountingLookup};
	my $records = $lookup
		? $lookup->(\@dependencies)
		: _slurm_dependency_accounting_lookup(\@dependencies, $optHR);
	# If accounting itself is unavailable, retain the original scheduler
	# dependencies.  This preserves the historical behaviour rather than
	# guessing that jobs have completed.
	return (join(';', @dependencies), '')
		unless defined($records) || $require_controller;
	$records ||= {};

	my $known_check = $optHR->{slurmDependencyKnownCheck};
	my (@remaining, @invalid);
	my %terminal = map { $_ => 1 } qw(
		BOOT_FAIL CANCELLED COMPLETED DEADLINE FAILED NODE_FAIL
		OUT_OF_MEMORY PREEMPTED REVOKED TIMEOUT
	);
	for my $dependency (@dependencies) {
		if ($require_controller) {
			my $known = $known_check
				? $known_check->($dependency)
				: _slurm_dependency_is_known($dependency);
			if ($known) {
				push @remaining, $dependency;
				next;
			}
		}

		my $record = $records->{$dependency};
		if ($record) {
			my $state = $record->{state} || '';
			my $exit_code = $record->{exit_code} || '';
			if ($state eq 'COMPLETED' && $exit_code =~ /^0:0$/) {
				# A successful prerequisite has already fulfilled both afterok and
				# afterany.  Omitting it also avoids Slurm's MinJobAge window.
				next;
			}
			if ($terminal{$state}) {
				# Any terminal job fulfils afterany, but unsuccessful terminal jobs
				# can never fulfil afterok.
				next if $after_any;
				push @invalid, "$dependency ($state, exit $exit_code)";
				next;
			}
			if ($require_controller) {
				push @invalid, "$dependency ($state in sacct, absent from slurmctld)";
			} else {
				push @remaining, $dependency;
			}
			next;
		}

		my $known = $known_check
			? $known_check->($dependency)
			: _slurm_dependency_is_known($dependency);
		if ($known) {
			# Accounting can lag briefly behind slurmctld for a new job.
			push @remaining, $dependency;
		} else {
			push @invalid, "$dependency (unknown to both sacct and slurmctld)";
		}
	}

	my $error = @invalid
		? 'Cannot submit a Slurm dependency chain: '.join(', ', @invalid)
		: '';
	return (join(';', @remaining), $error);
}

sub _run_slurm_submission {
	my ($command, $optHR) = @_;
	if (my $runner = $optHR->{slurmSubmissionRunner}) {
		return $runner->($command);
	}
	my $output = `$command 2>&1`;
	return ($output, $?);
}

sub _slurm_script_dependency {
	my ($script_path) = @_;
	open my $input, '<', $script_path
		or die "Cannot inspect Slurm job script $script_path: $!\n";
	while (my $line = <$input>) {
		if ($line =~ /^#SBATCH\s+--dependency=(afterok|afterany):([0-9:]+)\s*$/) {
			close $input;
			return ($1, split /:/, $2);
		}
	}
	close $input;
	return ('');
}

sub _rewrite_slurm_script_dependency {
	my ($script_path, $dependency_type, $dependencies) = @_;
	open my $input, '<', $script_path
		or die "Cannot read Slurm job script $script_path for dependency recovery: $!\n";
	my @lines = <$input>;
	close $input or die "Cannot close Slurm job script $script_path: $!\n";

	my $replacement = @{$dependencies}
		? "#SBATCH --dependency=$dependency_type:".join(':', @{$dependencies})."\n"
		: '';
	my $replaced = 0;
	for my $line (@lines) {
		next unless $line =~ /^#SBATCH\s+--dependency=/;
		$line = $replacement;
		$replaced++;
	}
	die "Cannot recover Slurm dependencies: expected one dependency directive in $script_path, found $replaced\n"
		unless $replaced == 1;

	my @script_stat = stat($script_path);
	my $temporary = "$script_path.dependency-rewrite.$$";
	open my $output, '>', $temporary
		or die "Cannot create dependency-recovery script $temporary: $!\n";
	print {$output} @lines
		or die "Cannot write dependency-recovery script $temporary: $!\n";
	close $output
		or die "Cannot close dependency-recovery script $temporary: $!\n";
	chmod($script_stat[2] & 07777, $temporary) if @script_stat;
	rename $temporary, $script_path
		or die "Cannot install recovered Slurm script $script_path: $!\n";
}

sub submitSlurmWithDependencyRecovery {
	my ($command, $script_path, $optHR) = @_;
	$optHR ||= {};
	my ($output, $status) = _run_slurm_submission($command, $optHR);
	return ($output, $status, 0)
		unless $status != 0 && $output =~ /Job dependency problem/i;

	my ($dependency_type, @dependencies) = _slurm_script_dependency($script_path);
	return ($output, $status, 0) unless $dependency_type ne '' && @dependencies;

	my %recovery_options = (%{$optHR}, slurmDependencyRequireController => 1);
	my ($reconciled, $dependency_error) = reconcileSlurmDependencies(
		join(';', @dependencies), $dependency_type eq 'afterany', \%recovery_options,
	);
	if ($dependency_error ne '') {
		$output .= "\nSlurm dependency recovery aborted: $dependency_error\n";
		return ($output, $status, 0);
	}

	my @reconciled = grep { length } split /;/, $reconciled;
	my %retained = map { $_ => 1 } @reconciled;
	my @removed = grep { !$retained{$_} } @dependencies;
	_rewrite_slurm_script_dependency($script_path, $dependency_type, \@reconciled);
	warn "Slurm rejected dependencies for $script_path; verified all prerequisites, "
		. "removed fulfilled job(s) ".(@removed ? join(',', @removed) : '<none>')
		. ", and retrying once with "
		. (@reconciled ? join(',', @reconciled) : 'no dependency directive')."\n";

	my ($retry_output, $retry_status) = _run_slurm_submission($command, $optHR);
	return ($retry_output, $retry_status, 1);
}

sub _continue_after_submission_failure {
	my ($optHR, $message) = @_;
	return 0 unless ($optHR->{continueOnSubmitError});
	$optHR->{submissionErrors} = []
		unless (ref($optHR->{submissionErrors}) eq 'ARRAY');
	push @{$optHR->{submissionErrors}}, $message;
	warn "$message\nMATAFILER will skip dependent jobs and continue with later work.\n";
	return 1;
}

sub handleSubmissionFailure {
	my ($optHR, $message) = @_;
	return $FAILED_SUBMISSION_DEPENDENCY
		if _continue_after_submission_failure($optHR, $message);
	die "$message\n";
}






sub randStr($){ #will be prefixed to jobname, to make jobs unique to each MF run
	my ($len) = @_;
	my @letters=('A'..'Z','a'..'z',1..9);
	my @letters2=('A'..'Z','a'..'z');
	my $total=scalar(@letters);
	my $newletter ="";
	$newletter = $letters2[int(rand scalar(@letters2))];
	for (my $i=1;$i<$len;$i++){
		$newletter .= $letters[rand $total];
	}
	return $newletter;
}



sub qsubSystem($ $ $ $ $ $ $ $ $ $){
	#args: 1[file to save bash & error & output] 2[actual bash cmd] 3[cores reserved for job]
	# 4["1G": Ram usage per core in GB] 5[0/1: synchronous execution] 6[name of job] 
	# 7[name of job dependencies, separated by ";"]
	# 8[0/1: excute in cwd?] 9[0/1: return qsub cmd or submit job to cluster]
	# Falk Hildebrand, may 2015
	my ($tmpsh,$cmd,$ncores,$memory,$jname,$waitJID,$cwd,$immSubm, $restrHostsAR, $optHR) = @_;
	my $requestedJobName = $jname;
	#$doSync, 5th arg
	#14,12G
	#die $tmpsh."\n";
	#my $jname = $tmpsh;
	#$jname =~ s/.*\///g;$jname =~ s/\.sh$//g;
	#\n#\$ -N $tmpsh
	return("") if ($cmd eq "");
	my $LSF = 0;
	my $qbin = "qsub";
	my $xtra = "";
	my $rTag = $optHR->{rTag};
	my $qmode = $optHR->{qmode};
	
	my $tmpScratchTag = $optHR->{tmpSpaceTag};
	#my $xtraNodeCmds = $optHR->{xtraNodeCmds};
	my $submissionConfig = $optHR->{submissionConfig};
	my @constrains = @{$optHR->{constraint}};# #SBATCH --constraint=
	#die "@constrains";
	my $lockFile = $optHR->{LOCKfile};
	my $nthreads= $ncores;
	if ($ncores =~ m/,/){my @spl = split /,/,$ncores;$ncores = $spl[1]; $nthreads=$spl[0];}
	#different format for bsub and slurm
	if ($memory =~ m/^[\.\d]+$/){$memory  = int($memory+0.5);}
	if ($memory =~ m/^0G$/){$memory  = "1G";} #most likely a rounding error from too many cores..
	if ($memory =~ s/G$//){$memory = int( ($memory* 1024 ) +0.5);};
	my $tmpSpace = convert2Gb( $optHR->{tmpSpace} );
	#die " $optHR->{tmpSpace}   $tmpSpace\n";
	#my $tmpSpace2 = $optHR->{tmpMinG};
	#my $wcKeysForJob = $optHR->{wcKeysForJob};
	my $exclNodes = $optHR->{excludeNodes};
	
	
	#die ($memory."\n");
	#my $queues = "\"".$optHR->{shortQueue}."\"";#"\"medium_priority\"";
	my $queues = "\"".$optHR->{medQueue}."\"";#"\"medium_priority\"";
	my $time = $optHR->{medTime};#"24:00:00";
	if ($optHR->{useHiMemQueue} == 1){
		$queues = "\"".$optHR->{highMemQueue}."\"";$optHR->{useHiMemQueue}=0;
	} elsif ($optHR->{useLongQueue} ==1){
		$queues = "\"".$optHR->{longQueue}."\"";#"\"medium_priority\"";
		#$time = "335:00:00";
		$optHR->{useLongQueue}=0;
	} elsif ($optHR->{useGPUQueue} ==1){
		$queues = "\"".$optHR->{gpuQueue}."\"";#"\"medium_priority\"";
		#$time = "23:00:00";
		$optHR->{useGPUQueue}=0;
	} elsif (defined $optHR->{useNetQueue} && $optHR->{useNetQueue} ==1){
		$queues = "\"".$optHR->{netQueue}."\"";
		$time = $optHR->{longTime};
		$optHR->{useNetQueue}=0;
	} elsif ($optHR->{useShortQueue} ==1){
		$queues = "\"".$optHR->{shortQueue}."\"";#"\"medium_priority\"";
		#$time = "00:45:00";
		$optHR->{useShortQueue}=0;
	}
	$waitJID = normalise_job_dependencies($waitJID);
	my @jspl = split /;/, $waitJID;
	my $has_failed_dependency = grep { $_ eq $FAILED_SUBMISSION_DEPENDENCY } @jspl;
	my $has_deferred_dependency = grep { $_ eq $DEFERRED_SUBMISSION_DEPENDENCY } @jspl;
	@jspl = grep { $_ ne $FAILED_SUBMISSION_DEPENDENCY
		&& $_ ne $DEFERRED_SUBMISSION_DEPENDENCY } @jspl;
	$waitJID = join(';', @jspl);
	my $dependency_error = '';
	if ($qmode eq 'slurm' && $optHR->{doSubmit} != 0 && $immSubm && @jspl) {
		for (@jspl) { s/^\Q$rTag\E//; }
		if (_slurm_dependencies_need_reconciliation(\@jspl, $optHR)) {
			my $reconciled;
			($reconciled, $dependency_error) = reconcileSlurmDependencies(
				join(';', @jspl), $optHR->{afterAny}, $optHR,
			);
			@jspl = split /;/, $reconciled;
			$waitJID = $reconciled;
			$has_failed_dependency = 1 if $dependency_error ne '';
			# A dependency confirmed as live cannot age out until at least
			# MinJobAge after it subsequently ends, so defer its next lookup.
			my $now = time;
			my $submitted_at = $optHR->{slurmDependencySubmittedAt} ||= {};
			for my $dependency (@jspl) {
				$submitted_at->{$dependency} = $now;
			}
		}
	}

	if ($cwd ne "" && !-d $cwd){system "mkdir -p $cwd";}
	#if ($memory > 250001){$queues = "\"scb\"";}
	$tmpsh =~ m/^(.*\/)[^\/]+$/;
	system "mkdir -p $1" unless (-d $1);
	open O,">",$tmpsh or die "Can't open qsub bash script $tmpsh\n";
	#die "$cmd\n";
	#print "$memory   $queues\n";
	#if (`hostname` !~ m/submaster/){
	if ($qmode eq "slurm"){$LSF = 2;$qbin="sbatch";
		#if ($memory > 250001){$queues = "\"bigmem\"";}
		##SBATCH --cpus-per-task=$ncores\n
		print O "#!/bin/bash\n#SBATCH -N 1\n#SBATCH --cpus-per-task=$ncores\n#SBATCH -o $tmpsh.otxt\n"; #\n#SBATCH -n  $ncores
		
		if ($nthreads != $ncores ){print O "#SBATCH --threads-per-core=1\n#SBATCH --hint=compute_bound\n";} #  specifically for iqtree/raxml
		print O "#SBATCH -e $tmpsh.etxt\n#SBATCH --mem=$memory\n#SBATCH --export=ALL\n";
		#print O "#SBATCH --kill-on-invalid-dep=yes\n";
		#print O "#SBATCH --tmp=$tmpSpace\n" if ($tmpSpace>0);#SBATCH --gres=ssd\n
		foreach my $subTerm ( split /;/, $submissionConfig){
			print O "#SBATCH $subTerm\n" if ($submissionConfig ne "");
		}
		if ($tmpSpace>0 && $tmpScratchTag ne ""){
			print O "#SBATCH $tmpScratchTag". int($tmpSpace+0.5) ."\n" ;
		}
		#"#SBATCH --gres=ssd"
		print O "#SBATCH -p $queues\n";
		print O "#SBATCH --gres=gpu:".$optHR->{gpuCount}."\n" if ($optHR->{gpuCount} > 0);
		#print O "#SBATCH --gres=tmp:${tmpSpace2}G\n" if ($tmpSpace2>0); #50g
		print O "#SBATCH --time=$time\n" unless ($time eq "");
		print O "#SBATCH --exclude=$exclNodes\n" unless ($exclNodes eq "");
		#print O "#SBATCH --localscratch=ssd:50\n"; #for EI cluster
		print O "#SBATCH --chdir=$cwd\n" if ($cwd ne "");
		print O "#SBATCH -J $rTag$jname\n" if ($jname ne "");
		print O "#SBATCH --wc=". $optHR->{wcKeysForJob} . "\n" if ($optHR->{wcKeysForJob} ne "");
		if (@constrains){
			print O "#SBATCH --constraint=". join(",",@constrains) ."\n" if (@constrains);
		}
		#foreach (@constrains){
	#		print O "#SBATCH --constraint=$_\n" if ($_ ne "");
		#}
		if (@jspl > 0) {
			for (@jspl) {s/^\Q$rTag\E//;}
			
			#$xtra .= "--dependency=afterok:".join(":",@jspl)." " if (@jspl > 0);
			if ($optHR->{afterAny}){
				print O "#SBATCH --dependency=afterany:".join(":",@jspl)."\n" ;
			} else {
				print O "#SBATCH --dependency=afterok:".join(":",@jspl)."\n" ;
			}
			#use this one for now, as slurm currently faults without a reason..
			#$xtra .= "--dependency=afterany:".join(":",@jspl)." " if (@jspl > 0);
		}

		#print O "#\$ -S /bin/bash\n#\$ -v LD_LIBRARY_PATH=".$optHR->{cpplib}."\n";##\$ -v TMPDIR=/dev/shm\n";
		#print O "#\$ -v PERL5LIB=".$optHR->{perl5lib}."\n"; #causes problems..
	} elsif ($qmode eq "bash"){
		$qbin="bash";$LSF=3;
		print O "#!/bin/bash\n";
	} elsif ($qmode eq "sge"){
		print O "#!/bin/bash\n#\$ -S /bin/bash\n#\$ -cwd\n#\$ -pe ".$optHR->{qsubPEenv}." $nthreads\n#\$ -o $tmpsh.otxt\n#\$ -e $tmpsh.etxt\n#\$ -l h_rss=$memory\n";#h_vmem=$mem\n";
		print O "#\$ -v LD_LIBRARY_PATH=".$optHR->{cpplib}."\n";#\$ -v TMPDIR=/dev/shm\n";
#		print O "#\$ -v PERL5LIB=".$optHR->{perl5lib}."\n";
		print O "#\$ -V\n";
	} else {
		$LSF = 1;$qbin="bsub";
		print O "#!/bin/bash\n";
		my $lsfLibraryPath = getProgPaths("lsfLDLibraryPath", 0);
		print O "export LD_LIBRARY_PATH=$lsfLibraryPath:\${LD_LIBRARY_PATH}\n\n"
			if $lsfLibraryPath ne "";
		#print O "export LD_LIBRARY_PATH=/g/bork3/home/hildebra/env/zlib-1.2.8:/g/bork3/x86_64/lib64:/lib:/lib64:/usr/lib64:\${LD_LIBRARY_PATH}\n\n";#:/g/software/linux/pack/python-2.7/lib/\nexport PATH=/g/bork3/home/zeller/py-virtualenvs/py2.7_bio1/bin/:\${PATH}\n\n";
		##BSUB -n $ncores\n#BSUB -o $tmpsh.otxt\n#BSUB -e $tmpsh.etxt\n#BSUB -M $mem\n#\$ -v LD_LIBRARY_PATH=/g/bork3/x86_64/lib64:/lib:/lib64:/usr/lib64\n#\$ -v TMPDIR=/dev/shm\n#BSUB -q medium_priority\n";
		my @restrHosts = @{$restrHostsAR};
		if ( @restrHosts > 0){
			$xtra .= " -m \"".join(" ",@restrHosts)."\" ";
			$queues = "\"medium_priority scb\"";
		}
		$xtra .= "-n $nthreads -oo $tmpsh.otxt -eo $tmpsh.etxt -q $queues -M $memory -R \"select[(mem>=$memory)] ";
		$xtra .= "rusage[tmp=$tmpSpace] " if ($tmpSpace>0);
		$xtra .= "span[hosts=1]\" -R \"rusage[mem=$memory]\" "; #
	}
	#set abortion on program fails
	print O "echo \"SLURM job ID: \$SLURM_JOB_ID\"\n" if ($qmode eq "slurm");
	print O "echo \$HOSTNAME;\n";
	print O "set -eo pipefail\n";
	print O "ulimit -c 0;\n";
	if ($tmpSpace > 0) {
		my $nodeTmpDir = getProgPaths("nodeTmpDir", 0);
		if (defined($nodeTmpDir) && $nodeTmpDir ne "") {
			die "Unsafe nodeTmpDir setting: $nodeTmpDir\n"
				if $nodeTmpDir =~ /[\r\n`;&|<>]/ || $nodeTmpDir =~ /\$\(/;
			$nodeTmpDir =~ s/"/\\"/g;
			print O "node_tmp_root=\"$nodeTmpDir\"\n";
			print O "node_tmp_job_id=\"\${SLURM_JOB_ID:-\${JOB_ID:-\${LSB_JOBID:-\$\$}}}\"\n";
			print O "node_tmp_workdir=\"\${node_tmp_root%/}/matafiler4.\${node_tmp_job_id}\"\n";
			print O "mkdir -p \"\$node_tmp_workdir\"\n";
			print O "export TMPDIR=\"\$node_tmp_workdir\"\n";
			print O "cd \"\$node_tmp_workdir\"\n";
		}
	}
	#any xtra commands (like module load perl?)
	print O "$optHR->{xtraNodeCmds}\n";
	#prevent core dump files
	#file location check availability
	#print O $optHR->{LocationCheckStrg};

	print O $cmd."\n";
	close O;
	#sleep (1);
	my $depSet=0;
	if ($LSF==2){#slurm
		if ($optHR->{doSync} == 1){$qbin = "srun";}
		
	} elsif ($LSF==3){ #bash
		$xtra = "";
	} elsif ($LSF==1){ #bsub #-M memLimit; -q queueName;  -m "host_name[@cluster_name]; -n minProcessors; 
		if ($optHR->{doSync} == 1){$xtra.="-K ";}
		if ($jname ne ""){$xtra.="-J $rTag$jname ";}
		if (@jspl > 0) {
			my @jspl = split(";",$waitJID);
			#remove empty elements
			@jspl = grep /\S/, @jspl;
			for (@jspl) { s/^\Q$rTag\E//; }
			if (@jspl > 0 ){
				$waitJID = join(") && done(",@jspl);
				$xtra.="-w \"done($waitJID)\" ";
			}
		}
		$tmpsh = " < ".$tmpsh;
	} else{ #qsub
		if ($optHR->{doSync} == 1){$xtra.="-sync y ";}
		if ($jname ne ""){$xtra.="-N $rTag$jname ";}
		if (@jspl > 0) {
			for (@jspl) { s/^\Q$rTag\E//; }
			if (@jspl > 0 ){$xtra.="-hold_jid ".join(",",@jspl) ." ";}
		}
			#$waitJID =~ s/;/,/g;$xtra.="-hold_jid $waitJID ";}
	}
	if ($cwd ne ""){if ($LSF==1) {$xtra.="-cwd $cwd"; }  elsif ($LSF == 0) {$xtra.="-wd $cwd";} }
	my $qcm = "$qbin $xtra $tmpsh \n";
	my $LOGhandle = "";
	if (exists $optHR->{LOG}){ $LOGhandle = $optHR->{LOG};}
	#if (@restrHosts > 0){die $qcm;}
	if ($optHR->{doSubmit} != 0 && $immSubm){
		return ($DEFERRED_SUBMISSION_DEPENDENCY, $qcm)
			if $has_deferred_dependency;

		if ($has_failed_dependency) {
			my $message = "Skipping submission for $tmpsh because an upstream submission failed";
			$message .= ": $dependency_error" if $dependency_error ne '';
			return ($FAILED_SUBMISSION_DEPENDENCY, $qcm)
				if (_continue_after_submission_failure($optHR, $message));
			die "$message\n";
		}
		my $capacityAvailable = qsubSystemWaitMaxJobs(
			$optHR->{maxConcurrentJobs} || 0,
			$optHR->{killDependencyNever} || 0,
			$optHR,
		);
		return ($DEFERRED_SUBMISSION_DEPENDENCY, $qcm)
			unless $capacityAvailable;
		system "rm -f $tmpsh.otxt $tmpsh.etxt";
		print $LOGhandle $qcm."\n" unless ($LOGhandle eq "" || !defined($LOGhandle) );
		#print("$qcm\n\n");
		#actual job excecution!
		my ($ret, $submit_status) = $LSF == 2
			? submitSlurmWithDependencyRecovery($qcm, $tmpsh, $optHR)
			: do {
				my $output = `$qcm`;
				($output, $?);
			};
		if ($submit_status != 0) {
			delete $optHR->{liveJobThrottleState};
			my $exit_code = $submit_status == -1 ? -1 : ($submit_status >> 8);
			my $message = "Job submission failed (exit $exit_code): $qcm$ret";
			return ($FAILED_SUBMISSION_DEPENDENCY, $qcm)
				if (_continue_after_submission_failure($optHR, $message));
			die $message;
		}
		if ($LSF == 2){#slurm get jobid
			chomp $ret;
			unless ($ret =~ /^Submitted batch job (\d+)\s*$/) {
				delete $optHR->{liveJobThrottleState};
				die "Could not parse Slurm job id from submission output: $ret\n";
			}
			$jname=$1;
			$optHR->{slurmDependencySubmittedAt} ||= {};
			$optHR->{slurmDependencySubmittedAt}{$jname} = time;
		} elsif ($LSF == 0) {
			die "Could not parse SGE job id from submission output: $ret\n"
				unless ($ret =~ /\bYour job(?:-array)?\s+(\d+)\b/);
			$jname=$1;
		} elsif ($LSF == 1) {
			die "Could not parse LSF job id from submission output: $ret\n"
				unless ($ret =~ /\bJob <(\d+)>/);
			$jname=$1;
		}
		$optHR->{submittedJobs} = 0 unless (defined $optHR->{submittedJobs});
		$optHR->{submittedJobs}++;
		if ($LSF == 3) {
			print "Completed local job $requestedJobName\n";
		} else {
			print "Submitted $requestedJobName as job $jname\n";
			$optHR->{submittedJobRecords}{$jname} = {
				requested_name => $requestedJobName,
			};
		}
		# Only record a lock owner after the scheduler has accepted the job.
		recordSampleLockJobs($lockFile, [$jname], $optHR) if $LSF != 3;
	}
	
	#die "$qcm\n";
	my $retJName = "$rTag$jname"; $retJName = "" if (!$immSubm); #return empty (for slurm), since no fwd job predictions..

	return ($retJName,$qcm);
}




sub _wanted_scheduler_jobs {
	my ($job_ids, $run_tag) = @_;
	my %wanted;
	return \%wanted unless defined $job_ids;
	for my $job_id (split /;/, normalise_job_dependencies($job_ids)) {
		$job_id =~ s/^\Q$run_tag\E// if defined($run_tag) && $run_tag ne '';
		$wanted{$job_id} = 1 if $job_id =~ /^\d+$/;
	}
	return \%wanted;
}

sub numUserJobs{
	return numLiveUserJobs(@_);
}

sub numLiveUserJobs {
	my ($optHR) = @_;
	my $rmSelf = @_ > 1 ? $_[1] : 0;
	my $jobIds = @_ > 2 ? $_[2] : undef;
	my $qmode = defined($optHR->{qmode}) ? $optHR->{qmode} : "slurm";
	return 0 if $qmode eq 'bash';
	my $wanted = _wanted_scheduler_jobs($jobIds, $optHR->{rTag} || '');
	return 0 if defined($jobIds) && !%{$wanted};
	my $snapshot = _scheduler_queue_snapshot(
		$optHR,
		runner => $optHR->{liveJobRunner} || $optHR->{pendingJobRunner},
		error => "Failed to count live user jobs\n",
	);
	my %live = %{$snapshot->{jobs}};
	delete $live{$ENV{SLURM_JOBID}}
		if ($rmSelf && $qmode eq "slurm" && ($ENV{SLURM_JOBID} || '') ne '');
	return scalar grep { !%{$wanted} || $wanted->{$_} } keys %live;
}

sub numActiveUserJobs{
	my ($optHR) = @_;
	my $rmSelf = @_ > 1 ? $_[1] : 0;
	my $jobIds = @_ > 2 ? $_[2] : undef;
	my $qmode = defined($optHR->{qmode}) ? $optHR->{qmode} : "slurm";
	return 0 if $qmode eq 'bash';
	my $wanted = _wanted_scheduler_jobs($jobIds, $optHR->{rTag} || '');
	return 0 if defined($jobIds) && !%{$wanted};
	my $snapshot = _scheduler_queue_snapshot(
		$optHR,
		runner => $optHR->{activeJobRunner},
		bare_state => $qmode eq 'slurm' ? 'RUNNING'
			: $qmode eq 'sge' ? 'R' : 'RUN',
		error => "Failed to count active user jobs\n",
	);
	my %active = map { $_ => 1 } grep {
		_scheduler_state_is_executing($qmode, $snapshot->{states}{$_})
	} keys %{$snapshot->{jobs}};
	delete $active{$ENV{SLURM_JOBID}}
		if ($rmSelf && $qmode eq "slurm" && ($ENV{SLURM_JOBID} || '') ne '');
	return scalar grep { !%{$wanted} || $wanted->{$_} } keys %active;
}


sub findQsubSys($){
	my $iniVal = "";
	$iniVal = $_[0] if (@_ > 0);
	#my $iniVal = "lsf";
	if ($iniVal ne ""){
		$iniVal = lc $iniVal; 
		$iniVal = "lsf" if ($iniVal eq "bsub");
		$iniVal = "sge" if ($iniVal eq "qsub");
		$iniVal = "slurm" if ($iniVal eq "sbatch");
	} else {
		$iniVal = "lsf";
		my $bpath = `which bsub  2>/dev/null`;chomp $bpath;my $bpresent=0; 
		$bpresent=1 if ($bpath !~ m/\n/ && -e $bpath);
		my $spath = `which sbatch  2>/dev/null`;chomp $spath;my $spresent=0; 
		$spresent=1 if ($spath !~ m/\n/ && -e $spath);
		my $qpath = `which qsub  2>/dev/null`; chomp $qpath;
		my $qpresent=0; $qpresent=1 if ($qpath !~ m/\n/ && -e $qpath);
		#print "$qpath\n";
		if ($spresent ){#slurm gets preference
			$iniVal="slurm";
		}elsif (!$bpresent && $qpresent){
			$iniVal = "sge";
		}elsif (!$qpresent && !$bpresent && !$spresent){
			die "No queueing system found (sbatch, qsub, or bsub). Use -qsubSystem bash for local execution.\n";
		}
	}
	#die;
	return $iniVal;
}
sub emptyQsubOpt{
	my ($doSubm) = $_[0];
	my $locChkStr = $_[1];
	my $qmode = "";
	$qmode = $_[2] if (@_ > 2);
	
	if (@_ > 2){$qmode = $_[2];}
	$qmode = findQsubSys($qmode);
	die "qsub system mode has to be \'lsf\', \'bash\', \'slurm\' or \'sge\'!\n" if ($qmode ne "lsf" &&$qmode ne "slurm" && $qmode ne "sge"&& $qmode ne "bash");
	my $MFdir = getProgPaths("MFLRDir");
	my $longQ = getProgPaths("longQueue",0); my $shortQ =  getProgPaths("shortQueue",0); my $medQ = getProgPaths("mediumQueue",1);
	#die "$shortQ\n";
	my $gpuQ = getProgPaths("gpuQueue",0);
	my $netQ = getProgPaths("netQueue",0);
	my $himemQ = getProgPaths("highMemQueue",0);
	if ($longQ eq ""){$longQ =  $medQ;}
	if ($medQ eq "" ){die "FATAL: no medium queue defined!\n";};
	if ($gpuQ eq "" ){$gpuQ = $medQ;};
	if ($netQ eq "" ){$netQ = $medQ;};
	if ($himemQ eq "" ){$himemQ = $medQ;};
	if ($shortQ eq "" ){$shortQ = $medQ;};
	my $xtraNodeCmds = getProgPaths("subXtraCmd",0);
	$xtraNodeCmds = "" unless (defined $xtraNodeCmds);
	my $medTime = getProgPaths("medTime",0);	my $shortTime = getProgPaths("shortTime",0);
	my $longTime = getProgPaths("longTime",0);
	my $subConfig = getProgPaths("submissionConfig",0);
	my @constr = ();
	if ($subConfig =~ s/--constraint=(\S+)//){
		#print "!!! $1\n";
		push(@constr, $1);
	}
	chomp($subConfig);
	@constr = grep(/\S/, @constr);
	#die "@constr\n$subConfig\nYW\n";
	
	#if ($qmode eq "slurm"){$shortQ = "htc"; $longQ="htc";}#$shortQ = "1day"; $longQ="1month";}
	my %ret = (
		rTag => randStr(3),
		doSubmit => $doSubm,
		LocationCheckStrg => $locChkStr,
		doSync => 0,
		longQueue => $longQ,
		gpuQueue => $gpuQ,
		netQueue => $netQ,
		highMemQueue => $himemQ,
		longTime => $longTime,#7days
		medQueue => $medQ,
		medTime => $medTime,#"24:00:00",
		shortQueue => $shortQ,
		shortTime => $shortTime, #2hrs
		useLongQueue => 0,
		useGPUQueue => 0,
		useNetQueue => 0,
		gpuCount => 0,
		useShortQueue => 0,
		useHiMemQueue => 0,
		submissionConfig => $subConfig,
		constraint => \@constr,
		qsubPEenv => getProgPaths("qsubPEenv"),
		perl5lib => "$MFdir:\$PERL5LIB",
		cpplib => "",
		tmpSpace => 15, #default was 15G; unit is G
		tmpSpaceTag => getProgPaths("nodeTmpDirTAG",0),
		LOCKfile => "",
		# Number of commands successfully handed to the configured execution
		# backend. Callers can snapshot this value to detect no-op passes.
		submittedJobs => 0,
		capacityCheckEverySubmissions => 10,
		jobPollSeconds => 20,
		#tmpMinG => 10,
		afterAny => 0,
		excludeNodes => "",
		xtraNodeCmds => $xtraNodeCmds,
		qmode => $qmode,
		wcKeysForJob => "",
		#LOG => undef,
	);
	#die "$MFdir\n";
	return \%ret;
}

sub qsubSystemJobAlive{
	my ($jAr,$optHR) = @_;
	my $killFailedJobs=0;
	$killFailedJobs = $_[2] if (@_ > 2);
	# A non-negative threshold makes loopTillComplete wait only until this many
	# submitted dependencies are actually executing.  Other callers retain the
	# historical behaviour of waiting for every queued dependency to disappear.
	my $activeThreshold=-1;
	$activeThreshold = $_[3] if (@_ > 3);
	my @jobs = split /;/, normalise_job_dependencies($jAr);
	@jobs = grep { $_ ne $FAILED_SUBMISSION_DEPENDENCY
		&& $_ ne $DEFERRED_SUBMISSION_DEPENDENCY } @jobs;
	return unless (@jobs);
	
	
	my $qmode = $optHR->{qmode};
	my $rTag = $optHR->{rTag};

	for (@jobs) {s/^\Q$rTag\E//;}
	return if $qmode eq "bash";
	my %wanted = map { $_ => 1 } @jobs;
	my $announced = 0;
	my $first_poll = 1;
	while (1) {
		my $snapshot = _scheduler_queue_snapshot(
			$optHR,
			runner => $optHR->{jobStatusRunner},
			reuse_seconds => $first_poll ? 2 : 0,
			error => "Failed to query active jobs\n",
		);
		$first_poll = 0;
		my @remaining = grep { $snapshot->{jobs}{$_} } keys %wanted;
		my @active = grep {
			$snapshot->{jobs}{$_}
				&& _scheduler_state_is_executing($qmode, $snapshot->{states}{$_})
		} keys %wanted;
		if ($activeThreshold >= 0 && @active <= $activeThreshold) {
			print scalar(@active)." active job(s) remain among ".scalar(@remaining).
				" queued dependencies; loop threshold $activeThreshold reached.\n";
			last;
		}
		last unless (@remaining);
		my $waitDescription = $activeThreshold >= 0
			? scalar(@active)." active job(s) among ".scalar(@remaining)."/".scalar(@jobs)." queued dependencies"
			: scalar(@remaining)."/".scalar(@jobs)." jobs";
		print "Waiting for $waitDescription to finish\n"
			unless ($announced++);
		if ($killFailedJobs){
			my $killed = qsubDepNeverKill();
			print " Killed $killed jobs with Dependency never completed\n" if ($killed > 0);
		}
		my $pollSeconds = defined($optHR->{jobPollSeconds})
			? 0 + $optHR->{jobPollSeconds} : 20;
		$pollSeconds = 1 if ($pollSeconds < 1);
		sleep($pollSeconds);
	}
	delete $optHR->{sampleLockActiveState};
	delete $optHR->{schedulerQueueState};
	return;
}

sub qsubDepNeverKill{
	my $srchCmd = "squeue -u \$USER -t PENDING -o \"\%8i \%.15R \%17E\"  | grep 'ependencyNe' | cut -f1 -d' ' | xargs  -t -i scancel {} | wc -l ";
	my $num = 0;
	$num = `$srchCmd`; chomp $num;
	return $num;
	
}


sub qsubSystemWaitMaxJobs{
	my ($checkMaxNumJobs) = @_;
	my $killPend = $_[1] if (@_ > 1);
	my $optHR = {}; $optHR = $_[2] if (@_ > 2);
	
	return 1 if ($checkMaxNumJobs <= 0);
	return 1 if exists($optHR->{doSubmit}) && !$optHR->{doSubmit};
	my $clock = $optHR->{schedulerClock};
	return 0 if $optHR->{nonblockingMaxConcurrentJobs}
		&& $optHR->{capacityDeferred};
	my $submittedNow = $optHR->{submittedJobs} || 0;
	my $refreshEvery = defined($optHR->{capacityCheckEverySubmissions})
		? int($optHR->{capacityCheckEverySubmissions}) : 10;
	$refreshEvery = 1 if $refreshEvery < 1;
	my $cached = $optHR->{liveJobThrottleState};
	if ($cached) {
		my $submittedSince = $submittedNow - ($cached->{submittedJobs} || 0);
		$submittedSince = 0 if $submittedSince < 0;
		# Treat every locally accepted job as live until the next scheduler
		# query. Refresh after a bounded submission batch, or sooner when this
		# conservative upper bound reaches the configured cap.
		my $estimatedLive = $cached->{liveJobs} + $submittedSince;
		return 1 if $estimatedLive < $checkMaxNumJobs
			&& $submittedSince < $refreshEvery;
	}
	my $num = numLiveUserJobs($optHR, 1);
	$optHR->{liveJobThrottleState} = {
		checkedAt => $clock ? $clock->() : time,
		liveJobs => $num,
		submittedJobs => $submittedNow,
	};
	if ($num >= $checkMaxNumJobs && $optHR->{nonblockingMaxConcurrentJobs}) {
		if ($killPend) {
			my $killed = qsubDepNeverKill();
			print " Killed $killed jobs with Dependency never completed\n"
				if $killed > 0;
		}
		$optHR->{capacityDeferred} = 1;
		unless ($optHR->{capacityDeferralAnnounced}) {
			print "running + pending jobs are at the $checkMaxNumJobs limit; "
				."deferring further submissions until the next loop pass.\n";
			$optHR->{capacityDeferralAnnounced} = 1;
		}
		return 0;
	}
	my $waitCnt = 0;
	while ($num >= $checkMaxNumJobs){
		if ($killPend){
			my $killed = qsubDepNeverKill();
			print " Killed $killed jobs with Dependency never completed\n" if ($killed > 0);
			
			#die;
		}
		print "waiting for running + pending jobs to fall below $checkMaxNumJobs...\n" if ($waitCnt == 0);
		sleep(40);
		$num = numLiveUserJobs($optHR, 1);
		$optHR->{liveJobThrottleState} = {
			checkedAt => $clock ? $clock->() : time,
			liveJobs => $num,
			submittedJobs => $optHR->{submittedJobs} || 0,
		};
		$waitCnt++;
		#print " $num ";
	}
	return 1;
}

sub qsubSystem2{
	my ($tmpsh,$optHR) = @_;
	my $hxr = {};
	$hxr = $_[2] if (@_ > 2 && defined $_[2]);
	my %xtras = %{$hxr}; 
	my $ncores = 0; 
	if (exists($xtras{cores})){$ncores = $xtras{cores};}
	my $nthreads= $ncores;
	if ($ncores =~ m/,/){my @spl = split /,/,$ncores;$ncores = $spl[1]; $nthreads=$spl[0];}
	if ($ncores != 0){#read in file, change it..
		open my $in,"<",$tmpsh or die "qsubSystem2: cant open $tmpsh\n";chomp(my @lines = <$in>); close $in;
		for (my $i=0;$i<@lines;$i++){
			if ($lines[$i] =~ m/--cpus-per-task/ || $lines[$i] =~ m/--mincpus/){
				$lines[$i] = "#SBATCH --cpus-per-task=$ncores";
				my $next_line = $i+1 < @lines ? $lines[$i+1] : "";
				$lines[$i] .= "\n#SBATCH --threads-per-core=1\n#SBATCH --hint=compute_bound" unless ($next_line =~ m/threads-per-core/);
			}
		}
		open my $out,">",$tmpsh or die "qsubSystem2: cant update $tmpsh\n";
		print {$out} join("\n",@lines), "\n";
		close $out or die "qsubSystem2: cant close updated $tmpsh\n";
	}
	my $xtra = "";
	my $qbin = "qsub";
	my $qmode = $optHR->{qmode};
	if ($qmode eq "slurm"){$qbin="sbatch";
	} elsif ($qmode eq "sge"){
	} elsif ($qmode eq "bash"){$qbin="bash";
	} else {$qbin="bsub";$xtra="<";
	}
	my $qcm = "$qbin $xtra $tmpsh \n";
	system($qcm) == 0 or die "qsubSystem2 failed: $qcm";
	return $qcm;
}



# Persist every job that currently owns this sample. MATAF4 removes the ledger
# itself on a later pass once none of these jobs remains in the scheduler.
sub MFnext($ $ $ $){
	my ($lckFile,$aR,$Jnum,$QSBoptHR) = @_;
	return if (! @{$aR});
	recordSampleLockJobs($lckFile, $aR, $QSBoptHR);
}


sub add2SampleDeps($ $){
	my ($ar1, $ar2) = @_;
	@{$ar1} = split /;/, normalise_job_dependencies($ar1, $ar2);
}
