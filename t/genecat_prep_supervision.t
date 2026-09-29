use strict;
use warnings;
no warnings 'once';
use FindBin qw($Bin);
use File::Spec;
use File::Path qw(make_path);
use File::Temp qw(tempdir);
use Test::More;
use lib "$Bin/lib";
use MFTestConfig;

my $root = File::Spec->rel2abs("$Bin/..");
sub read_file { open my $in, '<', $_[0] or die "$_[0]: $!"; local $/; return <$in>; }
sub write_file { open my $out, '>', $_[0] or die "$_[0]: $!"; print {$out} $_[1] or die $!; close $out or die $!; }

my $source = read_file("$root/secScripts/geneCat.pl");
my $code = <<'IMPORTS';
package GCPrep;
use strict;
use warnings;
use Fcntl qw(O_CREAT O_EXCL O_WRONLY);
use Errno qw(EEXIST ENOENT);
use IO::Handle;
use English;
use Mods::WorkflowResilience qw(retry_unlink write_workflow_record read_workflow_record);
use Mods::Subm qw(submissionDependencyDeferred submissionDependencyFailed);
our ($qsubDir);
IMPORTS
for my $name (qw(_sync_file _lock_stale_seconds _lock_heartbeat_seconds _lock_holder_label
        _read_append_lock _append_complete _break_stale_append_lock _append_file_locked
        _prep_record_paths _prep_batch_range _tail_lines _prep_job_accounting
        _prep_batch_outcome _prep_failure_report _supervise_prep_batches)) {
    my ($helper) = $source =~ /^(sub \Q$name\E\b[^\n]*\{.*?^\})/ms;
    die "Missing $name" unless defined $helper;
    $code .= "$helper\n";
}
eval $code;
die $@ if $@;
# fsync of a read-only handle is refused on some platforms (Cygwin); durability is not under test here.
{ no warnings "redefine"; *GCPrep::_sync_file = sub { 1 }; }

my $tmp = tempdir(CLEANUP => 1);
$GCPrep::qsubDir = "$tmp/LOGandSUB";
make_path("$tmp/LOGandSUB/preprocess");

subtest 'batch ranges cover every sample exactly once' => sub {
    my ($covered, $previous_end) = (0, 0);
    for my $batch (0 .. 10) {
        my ($from, $to) = GCPrep::_prep_batch_range($batch, 2008, 11);
        is($from, $previous_end, "batch $batch starts where the previous ended");
        $covered += $to - $from;
        $previous_end = $to;
    }
    is($covered, 2008, 'last sample is not dropped by rounding');
};

subtest 'outcome is decided by the heartbeat/failure records, not the scheduler' => sub {
    my $script = "$tmp/LOGandSUB/preprocess/Preprocess.4.sh";
    my ($heartbeat, $failure) = GCPrep::_prep_record_paths(4);
    my %fast = (settle_tries => 0, settle_seconds => 0, opts => {qmode => 'bash'});

    GCPrep::write_workflow_record($heartbeat, status => 'completed', stage => 'mode-subprepSmpls');
    ok(GCPrep::_prep_batch_outcome(4, '1', $script, %fast)->{ok}, 'completed heartbeat is success');

    unlink $heartbeat;
    write_file("$script.etxt", "noise\n" . join('', map { "line $_\n" } 1 .. 12));
    GCPrep::write_workflow_record($failure, status => 'failed', stage => 'mode-subprepSmpls',
        reason => "Timed out waiting for lock /x/compl.fna.gz.lock\n");
    my $lock = GCPrep::_prep_batch_outcome(4, '1', $script, %fast);
    ok(!$lock->{ok} && $lock->{retryable}, 'die message is transient and retried');
    like($lock->{reason}, qr/Timed out waiting for lock/, 'the recorded error is reported');
    is_deeply($lock->{output}, [map { "line $_" } 5 .. 12], 'last eight error lines are attached');

    GCPrep::write_workflow_record($failure, status => 'failed', stage => 'mode-subprepSmpls',
        reason => 'controller exit status 33');
    ok(!GCPrep::_prep_batch_outcome(4, '1', $script, %fast)->{retryable}, 'data errors (exit 33) are not retried');

    unlink $failure;
    GCPrep::write_workflow_record($heartbeat, status => 'running', stage => 'mode-subprepSmpls');
    my $killed = GCPrep::_prep_batch_outcome(4, '7', $script, %fast, opts => {
        qmode => 'slurm', jobAccountingRunner => sub { ("OUT_OF_MEMORY+|0:125|01:02:03\n", 0) },
    });
    ok($killed->{oom} && $killed->{retryable}, 'silent kill with OOM accounting requests more memory');
    like($killed->{reason}, qr/no completion record/, 'silent kill is explained');
    like($killed->{scheduler}, qr/OUT_OF_MEMORY.*01:02:03/, 'scheduler accounting is reported');
    unlink $heartbeat;
};

subtest 'supervisor resubmits failed batches and reports permanent failures' => sub {
    my (%launches, %memory_requested);
    my %script = (
        0 => [qw(ok)],
        1 => [qw(oom ok)],            # killed once, then fine
        2 => [qw(lock lock lock)],    # exhausts three attempts
        3 => [qw(fatal)],             # not retried
    );
    my $jobs = 100;
    my %job_batch;
    my $result = GCPrep::_supervise_prep_batches(
        count => 4, memory => 30, max_memory => 400, max_attempts => 3, rtag => 'ab',
        submit => sub {
            my ($batch, $attempt, $memory) = @_;
            $launches{$batch}++;
            $memory_requested{$batch}[$attempt - 1] = $memory;
            my $job = ++$jobs;
            $job_batch{$job} = [$batch, $attempt];
            return "ab$job";
        },
        cleanup => sub {},
        wait => sub { return [] },
        pause => sub {},
        outcome => sub {
            my ($batch, $job) = @_;
            my ($b, $attempt) = @{$job_batch{$job}};
            my $step = $script{$b}[$attempt - 1];
            return {ok => 1} if $step eq 'ok';
            return {ok => 0, retryable => 1, oom => 1, reason => 'killed'} if $step eq 'oom';
            return {ok => 0, retryable => 1, oom => 0, reason => 'Timed out waiting for lock',
                output => ['Timed out waiting for lock'], output_file => 'x.etxt'} if $step eq 'lock';
            return {ok => 0, retryable => 0, oom => 0, reason => 'exit 33'};
        },
    );
    local $SIG{__WARN__} = sub {};
    is_deeply($result->{failed}, [2, 3], 'only the exhausted and the non-retryable batch fail');
    is_deeply($result->{retried}, [1, 2], 'resubmitted batches are listed');
    is_deeply(\%launches, {0 => 1, 1 => 2, 2 => 3, 3 => 1}, 'attempt counts per batch');
    is_deeply($memory_requested{1}, [30, 45], 'an OOM raises the memory request for the retry');
    is_deeply($memory_requested{2}, [30, 30, 30], 'other failures keep the memory request');

    my $report = GCPrep::_prep_failure_report(2, $result->{state}{2}, 'samples 1-9');
    like($report, qr/batch 2 failed after 3 attempt\(s\) \(samples 1-9\)/, 'report names batch and attempts');
    like($report, qr/job\(s\): 103, 104, 107|job\(s\): \d+, \d+, \d+/, 'every job id is listed');
    like($report, qr/\| Timed out waiting for lock/, 'the job error output is shown to the user');
    like(GCPrep::_prep_failure_report(3, $result->{state}{3}, ''), qr/not a transient error/,
        'permanent failures say why they were not retried');
};

subtest 'submissions postponed by the job limit are not attempts' => sub {
    my (@submitted, $defer) = ();
    $defer = 2;
    my $result = GCPrep::_supervise_prep_batches(
        count => 1, memory => 30, max_attempts => 1, rtag => '',
        submit => sub {
            return Mods::Subm::deferredSubmissionDependency() if $defer-- > 0;
            push @submitted, $_[1];
            return '55';
        },
        cleanup => sub {}, wait => sub { [] }, pause => sub {}, outcome => sub { {ok => 1} },
    );
    ok($result->{ok}, 'batch finally runs');
    is_deeply(\@submitted, [1], 'deferrals did not consume the attempt budget');
};

subtest 'stale append lock is broken and its partial append rolled back' => sub {
    local $ENV{GENECAT_LOCK_STALE_SECONDS} = 60;
    my $dest = "$tmp/dest.gz";
    my $part = "$tmp/part.gz.0";
    my $lock = "$tmp/dest.lock";
    write_file($dest, "GOOD");
    write_file($part, "NEXT");
    # A killed job had appended garbage after recording size=4.
    write_file($dest, "GOODhalf-written");
    write_file($lock, "host=n1\npid=9\njob=5\nsize=4\ndest=$dest\nsource=$part\nmarker=$part.appended\n");
    my $old = time - 600;
    utime($old, $old, $lock) or die $!;
    local $SIG{__WARN__} = sub {};
    GCPrep::_append_file_locked($part, $dest, $lock, "$part.appended");
    is(read_file($dest), "GOODNEXT", 'partial bytes were removed before the new batch was appended');
    ok(!-e $lock && !-e $part, 'lock released and source consumed');
    ok(-e "$part.appended", 'completion marker written');

    # A lock whose append completed (marker present) must not lose data.
    write_file($dest, "GOODNEXT");
    write_file($part, "AGAIN");
    write_file("$tmp/other.gz.1.appended", "appended\n");
    write_file($lock, "host=n1\npid=9\nsize=4\ndest=$dest\nsource=$tmp/other.gz.1\nmarker=$tmp/other.gz.1.appended\n");
    utime($old, $old, $lock);
    GCPrep::_append_file_locked($part, $dest, $lock, "$part.appended2");
    is(read_file($dest), "GOODNEXTAGAIN", 'finished append of a dead holder is kept');
};

subtest 'live lock is respected and a waiting append proceeds after release' => sub {
    local $ENV{GENECAT_LOCK_STALE_SECONDS} = 3600;
    my $dest = "$tmp/dest2.gz";
    my $part = "$tmp/part2.gz.0";
    my $lock = "$tmp/dest2.lock";
    write_file($dest, "A");
    write_file($part, "B");
    write_file($lock, "host=n1\npid=1\nsize=1\ndest=$dest\nsource=x\nmarker=\n");
    my $pid = fork() // die $!;
    if (!$pid) { sleep 2; unlink $lock; exit 0; }
    local $SIG{__WARN__} = sub {};
    open my $quiet, '>', File::Spec->devnull or die $!;
    my $old_fh = select $quiet;
    GCPrep::_append_file_locked($part, $dest, $lock);
    select $old_fh;
    waitpid($pid, 0);
    is(read_file($dest), "AB", 'append waited for the fresh lock instead of failing after a fixed timeout');
};

done_testing();
