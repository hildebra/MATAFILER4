use strict;
use warnings;
use Test::More;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use lib "$Bin/..", "$Bin/lib";
use MFTestConfig;
use Mods::Subm qw(qsubSystem deferredSubmissionDependency);
use Mods::JobGraph qw(newJobGraph runJobGraph jobGraphReport jobGraphOutcome);

#supervised job graphs (geneCat functional-annotation stages): recording through qsubSystem, submission once the
#prerequisites completed, resubmission after an OOM kill with more memory, blocking of a failed job's dependents
my $tmp = tempdir(CLEANUP => 1);
sub write_file { open my $fh, '>', $_[0] or die $!; print {$fh} $_[1]; close $fh or die $!; }
sub read_file { open my $fh, '<', $_[0] or die $!; local $/; return <$fh>; }
sub queue_options {
	return {rTag=>'run',qmode=>'bash',doSubmit=>1,doSync=>0,
		tmpSpaceTag=>'',tmpSpace=>0,submissionConfig=>'',constraint=>[],LOCKfile=>'',
		excludeNodes=>'',medQueue=>'normal',medTime=>'',useHiMemQueue=>0,useLongQueue=>0,
		useGPUQueue=>0,useNetQueue=>0,useShortQueue=>0,wcKeysForJob=>'',xtraNodeCmds=>'',
		continueOnSubmitError=>1,submittedJobs=>0};
}

subtest 'qsubSystem records jobs while a graph is active' => sub {
	my $opt = queue_options();
	my $graph = newJobGraph();
	$opt->{jobGraph} = $graph;
	$opt->{useLongQueue} = 1; $opt->{constraint} = ['avx2'];
	my ($a) = qsubSystem("$tmp/rec/a.sh", "echo a", 4, "10G", "A", "", "", 1, [], $opt);
	is($a, 'jobgraph1', 'a graph id is returned');
	ok(!-e "$tmp/rec/a.sh", 'nothing is submitted (no script written)');
	is($opt->{useLongQueue}, 0, 'one-shot queue flags are reset as qsubSystem does');
	my ($b) = qsubSystem("$tmp/rec/b.sh", "echo b", 1, "1G", "B", "$a;run4711", "", 1, [], $opt);
	my ($ja, $jb) = @{$graph->{jobs}};
	is($ja->{opts}{useLongQueue}, 1, 'the job keeps the queue flag it was recorded with');
	is_deeply($ja->{opts}{constraint}, ['avx2'], 'and its constraints');
	ok(!exists $ja->{opts}{jobGraph}, 'recorded options do not record again');
	is_deeply($jb->{deps}, ['jobgraph1'], 'dependency on a recorded job is a graph edge');
	is($jb->{external}, 'run4711', 'other dependencies are passed on to the scheduler');
	is($ja->{memory}, '10G', 'memory request recorded');
	my ($none) = qsubSystem("$tmp/rec/c.sh", "", 1, "1G", "C", "", "", 1, [], $opt);
	is($none, '', 'an empty command records nothing, as qsubSystem submits nothing');
	is(scalar(@{$graph->{jobs}}), 2, 'two jobs recorded');
};

#a graph of A,B (chunks), C (needs A and B), D (chunk), E (needs D); outcomes scripted per job and attempt
sub scripted_graph {
	my $graph = newJobGraph();
	my $opt = { rTag => 'run', constraint => [] };
	$opt->{jobGraph} = $graph;
	my %id;
	$id{$_} = (Mods::JobGraph::jobGraphRecord($graph, "$tmp/g/$_.sh", "echo $_", 1, "10G", $_, "", "", 1, [], $opt))[0] for qw(A B D);
	$id{C} = (Mods::JobGraph::jobGraphRecord($graph, "$tmp/g/C.sh", "echo C", 1, "10G", 'C', "$id{A};$id{B}", "", 1, [], $opt))[0];
	$id{E} = (Mods::JobGraph::jobGraphRecord($graph, "$tmp/g/E.sh", "echo E", 1, "10G", 'E', $id{D}, "", 1, [], $opt))[0];
	return $graph;
}

subtest 'supervision: OOM retry with more memory, prerequisites first, failures block their dependents' => sub {
	my $graph = scripted_graph();
	my %plan = (A => [qw(oom ok)], B => [qw(ok)], C => [qw(ok)], D => [qw(fatal)], E => [qw(ok)]);
	my (@log, %attempt, %jobName);
	my $n = 100;
	my $res = runJobGraph($graph, { rTag => 'run' }, label => 'TestStage', max_attempts => 3,
		submit => sub {
			my ($script, $cmd, $cores, $mem, $name) = @_;
			$attempt{$name}++;
			push @log, "$name:$mem";
			like($cmd, qr/\ntouch \Q$script\E\.done\n$/, "$name writes its completion marker last") if ($attempt{$name} == 1);
			my $id = ++$n; $jobName{$id} = $name;
			return ("run$id", '');
		},
		wait => sub { return [] },
		pause => sub {},
		outcome => sub {
			my ($job, $id) = @_;
			my $step = $plan{$job->{name}}[$job->{attempts} - 1];
			return { ok => 1 } if ($step eq 'ok');
			return { ok => 0, retryable => 1, oom => 1, reason => 'killed for exceeding its memory' } if ($step eq 'oom');
			return { ok => 0, retryable => 0, oom => 0, reason => 'failed', output => ['Error: broken input'], output_file => 'D.etxt' };
		},
	);
	local $SIG{__WARN__} = sub {};
	ok(!$res->{ok}, 'stage not ok while a job failed');
	is_deeply([map { $_->{name} } @{$res->{failed}}], ['D'], 'the failing job is reported');
	is_deeply([map { $_->{name} } @{$res->{blocked}}], ['E'], 'its dependent is not run');
	ok(!grep({ /^E:/ } @log), 'E was never submitted');
	is_deeply([grep { /^A:/ } @log], ['A:10G', 'A:15G'], 'A resubmitted after the OOM with 1.5x memory');
	my ($cPos) = grep { $log[$_] =~ /^C:/ } 0 .. $#log;
	my ($a2Pos) = grep { $log[$_] eq 'A:15G' } 0 .. $#log;
	ok(defined($cPos) && $cPos > $a2Pos, 'C submitted only after its prerequisites completed');
	is($res->{state}{$graph->{jobs}[3]{id}}, 'done', 'C finished');
	my $report = jobGraphReport($res);
	like($report, qr/TestStage: 1 job\(s\) failed, 1 job\(s\) not run/, 'report summary');
	like($report, qr/D failed after 1 attempt\(s\): failed/, 'report names the job and reason');
	like($report, qr/\| Error: broken input/, 'and shows its error output');
	like($report, qr/not run: E/, 'and the jobs it blocked');
};

subtest 'retries are bounded; postponed submissions are not attempts' => sub {
	my $graph = newJobGraph();
	Mods::JobGraph::jobGraphRecord($graph, "$tmp/g/X.sh", "echo X", 1, "100G", 'X', "", "", 1, [], { constraint => [] });
	my ($defer, @mem) = (2);
	my $res = runJobGraph($graph, { rTag => '' }, max_attempts => 3,
		submit => sub { return (deferredSubmissionDependency(), '') if ($defer-- > 0); push @mem, $_[3]; return (scalar(@mem), ''); },
		wait => sub { [] }, pause => sub {},
		outcome => sub { return { ok => 0, retryable => 1, oom => 1, reason => 'killed for exceeding its memory' } },
	);
	local $SIG{__WARN__} = sub {};
	ok(!$res->{ok}, 'gives up after max_attempts');
	is_deeply(\@mem, ['100G', '150G', '225G'], 'three attempts with growing memory; deferrals did not count');
};

subtest 'outcome: completion marker and Slurm accounting' => sub {
	my $job = { script => "$tmp/o/job.sh", name => 'job' };
	mkdir "$tmp/o";
	my %fast = (settle_tries => 0);
	my $acc = sub { my $out = shift; return { qmode => 'slurm', jobAccountingRunner => sub { ($out, 0) } } };
	write_file("$tmp/o/job.sh.done", '');
	ok(jobGraphOutcome($job, '42', $acc->("42|OUT_OF_MEMORY|0:125\n"), %fast)->{ok},
		'completion marker present: success, even if Slurm reports an OOM event');
	unlink "$tmp/o/job.sh.done";
	my $r = jobGraphOutcome($job, '42', $acc->("42|FAILED|1:0\n42.batch|OUT_OF_MEMORY|0:125\n"), %fast);
	ok(!$r->{ok} && $r->{retryable} && $r->{oom}, 'OOM in a job step: retried with more memory');
	$r = jobGraphOutcome($job, '42', $acc->("42|FAILED|0:9\n"), %fast);
	ok($r->{oom}, 'killed by signal 9: treated as an OOM kill');
	$r = jobGraphOutcome($job, '42', $acc->("42|NODE_FAIL|0:0\n"), %fast);
	ok($r->{retryable} && !$r->{oom}, 'node failure: retried with the same memory');
	write_file("$tmp/o/job.sh.etxt", "Traceback\nValueError: bad\n");
	$r = jobGraphOutcome($job, '42', $acc->("42|FAILED|1:0\n42.batch|FAILED|1:0\n"), %fast);
	ok(!$r->{retryable}, 'a real error is not retried');
	is_deeply($r->{output}, ['Traceback', 'ValueError: bad'], 'its error output is kept for the report');
	like($r->{scheduler}, qr/FAILED, exit code 1:0/, 'with the scheduler state');
	$r = jobGraphOutcome($job, 'x', { qmode => 'sge' }, %fast);
	ok($r->{retryable} && $r->{oom}, 'no accounting: retried with more memory');
	write_file("$tmp/o/job.sh.etxt", "Traceback\nRuntimeError: Annotation worker timed out after 250s on a 125-seed sub-batch\n"
		. "WARNING: annotation workers did not exit within 120 s; force-terminating.\n");
	my $failed = $acc->("42|FAILED|1:0\n");
	$r = jobGraphOutcome($job, '42', $failed, %fast);
	ok(!$r->{retryable}, 'an error not declared transient is not retried');
	$r = jobGraphOutcome($job, '42', $failed, %fast, transient => qr/Annotation worker timed out/);
	ok($r->{retryable} && !$r->{oom}, 'a declared transient error is retried with the same memory');
	is($r->{reason}, 'failed with a transient error', 'and reported as such');
};

subtest 'runJobGraph passes the transient pattern to the outcome check' => sub {
	my $graph = newJobGraph();
	Mods::JobGraph::jobGraphRecord($graph, "$tmp/t/T.sh", "echo T", 1, "10G", 'T', "", "", 1, [], { constraint => [] });
	mkdir "$tmp/t";
	write_file("$tmp/t/T.sh.etxt", "RuntimeError: Annotation worker timed out after 250s\n");
	my $opt = { rTag => '', qmode => 'slurm', jobAccountingRunner => sub { ("1|FAILED|1:0\n", 0) } };
	my @mem;
	local $SIG{__WARN__} = sub {};
	my $res = runJobGraph($graph, $opt, max_attempts => 2, transient => qr/Annotation worker timed out/, settle_tries => 0,
		submit => sub { push @mem, $_[3]; return ('1', ''); }, wait => sub { [] }, pause => sub {});
	ok(!$res->{ok}, 'still failing after the last attempt');
	is_deeply(\@mem, ['10G', '10G'], 'resubmitted once, with the same memory');
};

subtest 'end to end through qsubSystem (local bash jobs)' => sub {
	my $opt = queue_options();
	my $graph = newJobGraph();
	$opt->{jobGraph} = $graph;
	mkdir "$tmp/e2e";
	my ($first) = qsubSystem("$tmp/e2e/first.sh", "echo one > $tmp/e2e/out.txt", 1, "1G", "same", "", "", 1, [], $opt);
	my ($second) = qsubSystem("$tmp/e2e/second.sh", "echo two >> $tmp/e2e/out.txt", 1, "1G", "same", $first, "", 1, [], $opt);
	delete $opt->{jobGraph};
	ok(!-e "$tmp/e2e/out.txt", 'nothing ran while recording');
	my $res;
	{ local *STDOUT; open STDOUT, '>', \my $sink; $res = runJobGraph($graph, $opt, label => 'E2E'); }
	ok($res->{ok}, 'both jobs completed (shared job names do not collide)');
	is(read_file("$tmp/e2e/out.txt"), "one\ntwo\n", 'in dependency order');
	ok(-e "$tmp/e2e/first.sh.done" && -e "$tmp/e2e/second.sh.done", 'completion markers written');
};

subtest 'geneCat runs both functional stages supervised' => sub {
	my $gc = read_file("$Bin/../secScripts/geneCat.pl");
	like($gc, qr/_submitStageOnce\('FuncAssign', sub \{ _superviseFuncStage\('FuncAssign', sub \{/, 'FuncAssign is supervised');
	like($gc, qr/_submitStageOnce\('FuncEMAP', sub \{ _superviseFuncStage\('FuncEMAP', sub \{/, 'FuncEMAP is supervised');
	like($gc, qr/useLongQueue\} = 1;[^\n]*\n\s*my \(\$dep,\$qcmd\) = qsubSystem\(\$qsubDir\."func_GC\.sh"/, 'FuncAssign controller on the long queue');
	like($gc, qr/useLongQueue\} = 1;[^\n]*\n\s*my \(\$dep,\$qcmd\) = qsubSystem\(\$qsubDir\."emap_GC\.sh"/, 'FuncEMAP controller on the long queue');
	like($gc, qr/local \$QSBoptHR->\{jobGraph\} = \$graph;/, 'jobs are recorded only while the stage is built');
	like($gc, qr/_inflightRecordJob\(\$stage, \$ENV\{SLURM_JOB_ID\}\)/, 'the in-flight marker names the live controller job');
	like($gc, qr/runJobGraph\([^;]*transient => qr\/Annotation worker timed out\//, 'eggNOG-mapper annotation timeouts are retried');
};

done_testing();
