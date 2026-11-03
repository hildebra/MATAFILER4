package Mods::JobGraph;
#Supervised job graphs. While a graph is active ($optHR->{jobGraph}), qsubSystem records each job instead of
#submitting it (jobGraphRecord); dependencies on recorded jobs become edges of the graph. runJobGraph then submits
#every job once all its prerequisites completed (no scheduler dependencies, so a failed job cannot leave the rest
#of the stage pending forever), waits for it, and decides its outcome from a completion marker the job writes as
#its last command. A job killed for lack of memory is resubmitted with more memory; a job that fails for good
#blocks only the jobs that depend on it.
#No compile-time dependency on Mods::Subm (which calls jobGraphRecord): its functions are called fully qualified.

use strict;
use warnings;
use File::Spec;
use Exporter qw(import);

our @EXPORT_OK = qw(newJobGraph jobGraphRecord runJobGraph jobGraphReport jobGraphOutcome);

#qsubSystem options that apply to the next job only; qsubSystem resets them after use
my @ONE_SHOT = qw(useHiMemQueue useLongQueue useGPUQueue useDownloadQueue useNetQueue useShortQueue);

sub newJobGraph{
	return { jobs => [], byId => {}, n => 0 };
}

#called by qsubSystem (same arguments) while a graph is active; returns the job's graph id like a job id
sub jobGraphRecord{
	my ($graph, $tmpsh, $cmd, $ncores, $memory, $jname, $waitJID, $cwd, $immSubm, $restrHostsAR, $optHR) = @_;
	return ("", "") if (!defined($cmd) || $cmd eq "");
	my %opts = %{$optHR};
	delete $opts{jobGraph};
	$opts{constraint} = [ @{ $optHR->{constraint} || [] } ];
	$optHR->{$_} = 0 foreach (grep { exists $optHR->{$_} } @ONE_SHOT);
	my $id = "jobgraph" . (++$graph->{n});
	my (@deps, @external);
	foreach my $d (grep { $_ ne '' } split /[;,\s]+/, (defined($waitJID) ? $waitJID : '')){
		if (exists $graph->{byId}{$d}){ push @deps, $d; } else { push @external, $d; }
	}
	my $job = { id => $id, script => File::Spec->rel2abs($tmpsh), cmd => $cmd, cores => $ncores,
		memory => $memory, name => (defined($jname) ? $jname : ''), deps => \@deps,
		external => join(';', @external), cwd => (defined($cwd) ? $cwd : ''),
		hosts => ($restrHostsAR || []), opts => \%opts };
	push @{$graph->{jobs}}, $job;
	$graph->{byId}{$id} = $job;
	return ($id, "");
}

sub _memGB{
	my ($m) = @_;
	return (defined($m) && $m =~ /^(\d+(?:\.\d+)?)G$/i) ? $1 : undef;
}

#scheduler accounting of a finished job (Slurm only); {} if unavailable
sub _accounting{
	my ($jobId, $optHR) = @_;
	return {} unless (($optHR->{qmode} // '') eq 'slurm' && defined($jobId) && $jobId =~ /^\d+$/);
	my $command = "sacct -n -P -j $jobId --format=JobIDRaw,State,ExitCode";
	my ($output, $status);
	if ($optHR->{jobAccountingRunner}){
		($output, $status) = $optHR->{jobAccountingRunner}->($command);
		$status //= 0;
	} else {
		$output = `$command 2>/dev/null`;
		$status = $?;
	}
	return {} if ($status || !defined($output));
	my %acc = (oom => 0, signal9 => 0);
	foreach my $line (split /\n/, $output){
		my ($raw, $state, $exit) = map { defined($_) ? $_ : '' } split /\|/, $line, 3;
		next unless ($raw =~ /^\Q$jobId\E(?:[._].*)?$/);
		$state = uc($state); $state =~ s/^\s+//; $state =~ s/[\s+].*$//; #"CANCELLED by 123", "OUT_OF_MEMORY+"
		$acc{oom} = 1 if ($state =~ /^OUT_OF_ME/);
		$acc{signal9} = 1 if ($exit =~ /:9$/ || $exit =~ /^137:/);
		if ($raw eq $jobId){ $acc{state} = $state; $acc{exit} = $exit; }
	}
	return {} unless (defined($acc{state}));
	$acc{summary} = "scheduler state $acc{state}, exit code $acc{exit}";
	return \%acc;
}

sub _tail{
	my ($path, $count) = @_;
	return () unless (defined($path) && -s $path);
	open my $in, '<', $path or return ();
	my @ring;
	while (my $l = <$in>){
		$l =~ s/[\r\n]+\z//;
		next unless ($l =~ /\S/);
		push @ring, (length($l) > 300 ? substr($l, 0, 300) . ' [...]' : $l);
		shift @ring if (@ring > $count);
	}
	close $in;
	return @ring;
}

#outcome of a job that left the queue: {ok} or {ok => 0, retryable, oom, reason, scheduler, output}
#%o: transient (regex): a failure whose error output ends with a matching line is retried with the same memory
sub jobGraphOutcome{
	my ($job, $jobId, $optHR, %o) = @_;
	my $marker = "$job->{script}.done";
	my $tries = defined($o{settle_tries}) ? $o{settle_tries} : 4;
	for my $try (0 .. $tries){
		return { ok => 1 } if (-e $marker);
		#shared filesystems can show a finished job's files a few seconds late
		($o{sleeper} || sub { sleep $_[0] })->(defined($o{settle_seconds}) ? $o{settle_seconds} : 5) if ($try < $tries);
	}
	my %r = (ok => 0, retryable => 0, oom => 0);
	my $acc = _accounting($jobId, $optHR);
	if (!defined($acc->{state})){
		#no accounting (other schedulers, sacct unavailable): the job died before its last command; as for an
		#out-of-memory kill, try again with more memory
		@r{qw(retryable oom)} = (1, 1);
		$r{reason} = "job left the queue without completing (no scheduler accounting)";
	} elsif ($acc->{oom} || ($acc->{state} eq 'FAILED' && $acc->{signal9})){
		@r{qw(retryable oom)} = (1, 1);
		$r{reason} = "killed for exceeding its memory";
	} elsif ($acc->{state} =~ /^(?:NODE_FAIL|PREEMPTED|BOOT_FAIL|REQUEUED)/){
		$r{retryable} = 1;
		$r{reason} = "lost to a scheduler/node problem";
	} else {
		$r{reason} = "failed";
		if (defined($o{transient}) && grep { $_ =~ $o{transient} } _tail("$job->{script}.etxt", 40)){
			$r{retryable} = 1;
			$r{reason} = "failed with a transient error";
		}
	}
	$r{scheduler} = $acc->{summary} if ($acc->{summary});
	my @out = _tail("$job->{script}.etxt", 8);
	@out = _tail("$job->{script}.otxt", 5) unless (@out);
	$r{output} = \@out;
	$r{output_file} = "$job->{script}.etxt";
	return \%r;
}

#submits the recorded jobs as their prerequisites complete and supervises them until the graph is finished.
#%o: max_attempts (3), memory_factor (1.5), label, transient (see jobGraphOutcome), and for tests: submit, wait, outcome,
#pause, settle_tries
sub runJobGraph{
	my ($graph, $optHR, %o) = @_;
	my $maxAttempts = $o{max_attempts} || 3;
	my $factor = $o{memory_factor} || 1.5;
	my $label = defined($o{label}) ? $o{label} : 'job graph';
	my $rtag = defined($optHR->{rTag}) ? $optHR->{rTag} : '';
	require Mods::Subm;
	my $submit = $o{submit} || sub { return &Mods::Subm::qsubSystem(@_); };
	my $wait = $o{wait} || sub { return Mods::Subm::qsubSystemJobAlive($_[0], $optHR, 0, -1, 300) || []; };
	my @outcomeOpts = map { defined($o{$_}) ? ($_ => $o{$_}) : () } qw(transient settle_tries);
	my $outcome = $o{outcome} || sub { return jobGraphOutcome($_[0], $_[1], $optHR, @outcomeOpts); };
	my $pause = $o{pause} || sub { sleep $_[0] };
	my @jobs = @{$graph->{jobs}};
	my %state = map { $_->{id} => 'pending' } @jobs;
	foreach my $job (@jobs){ $job->{attempts} = 0; $job->{failures} = []; $job->{jobIds} = []; $job->{mem} = $job->{memory}; }
	my %dependents;
	foreach my $job (@jobs){ push @{$dependents{$_}}, $job->{id} foreach (@{$job->{deps}}); }
	my $block; $block = sub {
		foreach my $d (@{$dependents{$_[0]} || []}){
			next unless ($state{$d} eq 'pending');
			$state{$d} = 'blocked'; $block->($d);
		}
	};
	my (%running, @retried);
	my $announced = "";
	print "$label: supervising ".scalar(@jobs)." job(s)\n";
	while (1){
		my @ready = grep { $state{$_->{id}} eq 'pending' && !grep { $state{$_} ne 'done' } @{$_->{deps}} } @jobs;
		my $deferred = 0;
		foreach my $job (@ready){
			if ($job->{attempts} > 0){ #keep the logs of the failed attempt; qsubSystem removes them
				foreach my $ext (qw(etxt otxt)){
					my $log = "$job->{script}.$ext";
					rename $log, "$log.attempt$job->{attempts}" if (-e $log);
				}
			}
			unlink "$job->{script}.done";
			my %opts = %{$job->{opts}};
			$opts{constraint} = [ @{ $job->{opts}{constraint} || [] } ];
			my ($id) = $submit->($job->{script}, "$job->{cmd}\ntouch $job->{script}.done\n", $job->{cores}, $job->{mem},
				$job->{name}, $job->{external}, $job->{cwd}, 1, $job->{hosts}, \%opts);
			if (defined($id) && $id ne '' && Mods::Subm::submissionDependencyDeferred($id)){
				$deferred++; next; #postponed by the concurrent-job limit: not an attempt
			}
			$job->{attempts}++;
			if (!defined($id) || $id eq '' || Mods::Subm::submissionDependencyFailed($id)){
				push @{$job->{failures}}, { ok => 0, retryable => 0, reason => 'the scheduler rejected the submission' };
				$state{$job->{id}} = 'failed'; $block->($job->{id});
				next;
			}
			(my $bare = $id) =~ s/^\Q$rtag\E//;
			push @{$job->{jobIds}}, $bare;
			$state{$job->{id}} = 'running';
			#keyed by graph id: local (bash) runs return the job name, which several jobs may share
			$running{$job->{id}} = [$job, $bare];
		}
		last if (!%running && !$deferred);
		if (%running){
			my %alive = map { $_ => 1 } @{ $wait->([ do { my %s; grep { !$s{$_}++ } map { $_->[1] } values %running } ]) || [] };
			foreach my $gid (sort keys %running){
				my ($job, $bare) = @{$running{$gid}};
				next if ($alive{$bare});
				delete $running{$gid};
				my $r = $outcome->($job, $bare);
				if ($r->{ok}){ $state{$job->{id}} = 'done'; next; }
				push @{$job->{failures}}, $r;
				warn "$label: $job->{name} (job $bare) $r->{reason} on attempt $job->{attempts}/$maxAttempts"
					. ($r->{scheduler} ? " ($r->{scheduler})" : "") . "\n";
				if ($r->{retryable} && $job->{attempts} < $maxAttempts){
					if ($r->{oom} && defined(_memGB($job->{mem}))){
						my $more = int(_memGB($job->{mem}) * $factor + 0.5) . "G";
						warn "  raising the memory request of $job->{name} from $job->{mem} to $more\n";
						$job->{mem} = $more;
					}
					warn "  resubmitting $job->{name} (attempt " . ($job->{attempts} + 1) . "/$maxAttempts)\n";
					push @retried, $job->{id};
					$state{$job->{id}} = 'pending';
				} else {
					$state{$job->{id}} = 'failed'; $block->($job->{id});
				}
			}
		} else {
			$pause->(30); #only submissions postponed by the job limit are left
		}
		my %n; $n{$state{$_}}++ foreach (keys %state);
		my $progress = join(", ", map { "$n{$_} $_" } grep { $n{$_} } qw(done running pending failed blocked));
		if ($progress ne $announced){ print "$label: $progress\n"; $announced = $progress; }
	}
	my @failed = grep { $state{$_->{id}} eq 'failed' } @jobs;
	my @blocked = grep { $state{$_->{id}} ne 'failed' && $state{$_->{id}} ne 'done' } @jobs;
	return { ok => (!@failed && !@blocked) ? 1 : 0, failed => \@failed, blocked => \@blocked,
		retried => [ do { my %s; grep { !$s{$_}++ } @retried } ], state => \%state, label => $label };
}

sub jobGraphReport{
	my ($res) = @_;
	my @lines = ("$res->{label}: " . scalar(@{$res->{failed}}) . " job(s) failed, " . scalar(@{$res->{blocked}})
		. " job(s) not run because a prerequisite failed:");
	foreach my $job (@{$res->{failed}}){
		my $last = $job->{failures}[-1] || {};
		push @lines, "  $job->{name} failed after $job->{attempts} attempt(s): " . ($last->{reason} // 'unknown');
		push @lines, "    script: $job->{script}";
		push @lines, "    job(s): " . join(', ', @{$job->{jobIds}}) if (@{$job->{jobIds}});
		push @lines, "    scheduler: $last->{scheduler}" if ($last->{scheduler});
		push @lines, "    not retried: this is not a memory or node problem" unless ($last->{retryable});
		if (@{$last->{output} || []}){
			push @lines, "    last lines of $last->{output_file}:";
			push @lines, map { "      | $_" } @{$last->{output}};
		}
	}
	push @lines, "  not run: " . join(', ', map { $_->{name} } @{$res->{blocked}}) if (@{$res->{blocked}});
	return join("\n", @lines) . "\n";
}

1;
