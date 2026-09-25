use strict;
use warnings;
use Test::More;
use JSON::PP qw(decode_json);
use File::Temp qw(tempdir);
use File::Path qw(make_path);
use FindBin qw($Bin);
use lib "$Bin/..", "$Bin/lib";
use MFTestConfig;
use Mods::Subm qw(qsubSystem qsubSystemWaitMaxJobs deferredSubmissionDependency
    submissionDependencyDeferred submissionDependencyFailed handleSubmissionFailure recordSampleLockJobs
    submitSlurmWithDependencyRecovery add2SampleDeps);
use Mods::WorkflowControl qw(normalise_job_dependencies deferred_command_dependencies
    augment_deferred_submission hybrid_group_ready hybrid_package_complete);
use Mods::GenoMetaAss qw(fileGZe);

my $tmp = tempdir(CLEANUP => 1);
sub read_file { open my $fh, '<', $_[0] or die $!; local $/; return <$fh>; }
sub write_file { open my $fh, '>', $_[0] or die $!; print {$fh} $_[1]; close $fh or die $!; }
my $source = read_file("$Bin/../MATAF4.pl");
my %runOptions = (submit => 1);
my %MFconfig = (checkMaxNumJobs => 0, killDepNever => 0, silent => 1);
my %MFopt = (DoAssembly => 5);
my %checkpointNames = (preAssemblyDone => 'pre.done', assemblyDone => 'final.done');
my (%AsGrps, %map);
my $curSmpl = 'terminal';
my $QSBoptHR;
my $baseOut = "$tmp/";
my $logDir = "$tmp/";
my $JNUM = 1;
my $warnings = '';
local $SIG{__WARN__} = sub { $warnings .= $_[0]; };
for my $name (qw(postSubmQsub prepPreAssmbl nopareil)) {
    my ($body) = $source =~ /(^sub \Q$name\E[^\n]*\{.*?^\})/ms;
    die "Missing $name" unless defined $body;
    $body =~ s/^(sub \Q$name\E)\(\)/$1/; # main forward declarations have no prototype
    eval $body; die $@ if $@;
}
sub movePreAssmData { die 'Incomplete inputs must not be packaged'; }
sub queue_options {
    return {rTag=>'run',qmode=>'bash',doSubmit=>1,doSync=>0,
        tmpSpaceTag=>'',tmpSpace=>0,submissionConfig=>'',constraint=>[],LOCKfile=>'',
        excludeNodes=>'',medQueue=>'normal',medTime=>'',useHiMemQueue=>0,useLongQueue=>0,
        useGPUQueue=>0,useNetQueue=>0,useShortQueue=>0,wcKeysForJob=>'',xtraNodeCmds=>'',
        continueOnSubmitError=>1,submittedJobs=>0};
}

for my $marker ('__MF4_SUBMISSION_FAILED__', deferredSubmissionDependency()) {
    $QSBoptHR = queue_options();
    my ($job, $command) = qsubSystem("$tmp/blocked.sh", "touch $tmp/unsafe", 1, '1G',
        'blocked', $marker, '', 0, [], $QSBoptHR);
    is($job, $marker, "deferred construction propagates $marker");
    is($command, '', 'blocked construction exposes no runnable deferred command');
    is($QSBoptHR->{submittedJobs}, 0, 'blocked construction does not submit');

    write_file("$tmp/consumer.sh", "#!/bin/bash\ntouch $tmp/unsafe\n");
    my $result = postSubmQsub("$tmp/audit.sh", "bash $tmp/consumer.sh\n", $marker);
    is($result, $marker, "deferred release propagates $marker");
    ok(!-e "$tmp/unsafe", 'blocked deferred consumer never executes');
    is($QSBoptHR->{submittedJobs}, 0, 'blocked release does not increment accepted jobs');
    unlink "$tmp/unsafe";
}
$QSBoptHR = queue_options();
$QSBoptHR->{continueOnSubmitError} = 0;
my $ok = eval { postSubmQsub("$tmp/audit.sh", "bash $tmp/consumer.sh\n", '__MF4_SUBMISSION_FAILED__'); 1 };
ok(!$ok, 'deferred release retains fail-fast policy when continuation is disabled');
like($@, qr/upstream submission failed/, 'fail-fast explains the failed prerequisite');
unlink "$tmp/unsafe";
$QSBoptHR = queue_options();
is(postSubmQsub("$tmp/audit.sh", "bash $tmp/consumer.sh\n", ''), '',
    'successful local deferred work returns no asynchronous job ID');
ok(-e "$tmp/unsafe", 'healthy local deferred work still executes');
is($QSBoptHR->{submittedJobs}, 1, 'healthy release is counted once');

for my $case (
    ['empty last primary', 1, 1, 1, 1],
    ['no-primary last', 0, 0, 1, 1],
    ['empty before missing package', 1, 1, 0, 0],
    ['nonempty missing preassembly', 1, 0, 1, 0],
) {
    my ($name, $has_primary, $empty, $packages, $ready) = @$case;
    %AsGrps = (group => {CntAimAss=>2,CntPreAss=>$packages,SupportReads=>'ONT:support'});
    %map = (terminal => {hasPrimaryRds=>$has_primary,inputFilesEmpty=>$empty,SupportReads=>'ONT:support'});
    my @state = prepPreAssmbl("$tmp/pre", "$tmp/package", "$tmp/map", "$tmp/scratch",
        "$tmp/stats", 'group', "$tmp/final/scaffolds.fasta.filt", "$tmp/final");
    is($state[2], $ready, "$name: final hybrid release reflects package readiness");
    is($state[1], ($empty || !$has_primary) ? 0 : 1,
        "$name: empty members do not request their own preassembly");
    is($AsGrps{group}{CntPreAssNoPrim}, ($empty || !$has_primary) ? 1 : 0,
        "$name: ineligible member is counted exactly once");
}
make_path("$tmp/final");
write_file("$tmp/final/scaffolds.fasta.filt", ">ctg\nACGT\n");
write_file("$tmp/final/final.done", '');
my @finished = prepPreAssmbl("$tmp/pre", "$tmp/package", "$tmp/map", "$tmp/scratch",
    "$tmp/stats", 'group', "$tmp/final/scaffolds.fasta.filt", "$tmp/final");
is_deeply(\@finished, [0,0,0,0], 'published final assembly bypasses preassembly gating');


# Exercise the release helper with a fake Slurm endpoint, including the
# original per-sample prerequisites embedded before the assembly job existed.
$QSBoptHR = queue_options();
$QSBoptHR->{qmode} = 'slurm';
my @accepted;
$QSBoptHR->{slurmSubmissionRunner} = sub {
    my ($command) = @_;
    push @accepted, $command;
    return ('Submitted batch job '.(800 + @accepted)."\n", 0);
};
my @commands;
for my $i (1, 2) {
    my (undef, $command) = qsubSystem("$tmp/map$i.sh", 'echo mapped', 1, '1G',
        "map$i", 'run'.(700+$i), '', 0, [], $QSBoptHR);
    push @commands, $command;
}
is(postSubmQsub("$tmp/group.sh", join('', @commands), 'run900;run900'), 'run801;run802',
    'group release returns the actual mapping job IDs');
like(read_file("$tmp/map1.sh"), qr/^#SBATCH --dependency=afterok:701:900$/m,
    'first mapping waits for its own staging and the shared assembly');
like(read_file("$tmp/map2.sh"), qr/^#SBATCH --dependency=afterok:702:900$/m,
    'second mapping waits for its own staging and the shared assembly');
is(scalar(@accepted), 2, 'independent group mappings each submit once');

# Evaluate the actual early-exit bookkeeping, rather than reconstructing it.
my @sampleDeps;
my @EBIjobs;
my $smplTmpDir = "$tmp/scratch/";
my $jdep = 'run500';
sub uploadRawFilePrep { return $_[3] ? 'run502' : 'run501'; }
my ($upload_block) = $source =~ /(\tmy \$uplJob = uploadRawFilePrep.*?\n\tadd2SampleDeps\([^\n]+\);)/s;
ok(defined($upload_block), 'locate upload bookkeeping in main controller');
eval $upload_block; die $@ if $@;
is_deeply(\@sampleDeps, ['run501','run502'], 'both upload scratch owners enter sample cleanup dependencies');
my $primaryCleanLibraries = [];
my $nonParDir = "$tmp/nonpareil/";
my $SmplName = 'sample';
my $primaryDep = 'run500;run503';
my $mergJbN = 'run504';
my $cAssGrp = 'group';
$AsGrps{group}{readDeps} = 'run505';
sub libraryFiles { return ['fixture.fastq']; }
sub getProgPaths { return 'unused-nonpareil'; }
my ($np_block) = $source =~ /(\t\tmy \$globalNPD = \$baseOut\."NonPareil\/";.*?)\t\tMFnext/s;
ok(defined($np_block), 'locate Nonpareil bookkeeping before its early exit');
my @np_submissions;
{
    no warnings 'redefine';
    local *qsubSystem = sub ($$$$$$$$$$) { push @np_submissions, [@_]; return ('run506', ''); };
    eval $np_block; die $@ if $@;
}
is_deeply(\@sampleDeps, [map { "run$_" } (501,502,500,503,504,505,506)],
    'Nonpareil early exit retains staging, filtering, merging, profiling, uploads and its own job');
is($np_submissions[0][5], $primaryDep, 'Nonpareil waits for its read producers');
make_path($nonParDir);
write_file("$nonParDir/sample.npo", "complete\n");
is(nopareil(['fixture.fastq'], $nonParDir, "$tmp/NonPareil/", 'sample', ''), '',
    'completed Nonpareil returns no invented scheduler dependency');
is(nopareil(['fixture.fastq'], $nonParDir, "$tmp/NonPareil/", 'sample', 'run503'), 'run503',
    'completed Nonpareil preserves a real caller dependency');


open my $probe, '-|', $^X, "$Bin/../docs/audits/2026-09-11/submission-release-probe.pl"
    or die "Cannot run group release fixture: $!";
my $release = decode_json(do { local $/; <$probe> });
ok(close($probe), 'actual terminal-group release fixture executes');
my ($finalizer) = grep { $_->{job} eq 'empty-finalizer' } @{$release->{events}};
my ($binner) = grep { $_->{job} eq 'binning' } @{$release->{events}};
ok(!defined($binner), 'terminal-empty release never starts binning without the normal statistics chain');
is($finalizer->{input_dependency}, join(';', @{$release->{sample_dependencies}}),
    'terminal cleanup receives every shared job dependency instead of running immediately');
like($finalizer->{input_dependency}, qr/(?:^|;)run202(?:;|$)/,
    'terminal cleanup still waits for released mapping jobs');
is($release->{completed_empty_member_case}{empty_members_accounted}, 1,
    'completed empty member is counted before its early return');
is($release->{completed_empty_member_case}{final_hybrid_ready}, 1,
    'a later complete package releases the final hybrid assembly');


# Count all validated empty terminal representations, but not successful or
# failed nonempty samples. Execute the actual early-return branch each time.
my ($closed_body) = $source =~ /(^\tif \(\$closedSample\) \{.*?^\t\})/ms;
my ($closedSample, %progStats);
my $smplLockF = "$tmp/sample.lock";
my %loopSampleCompleted;
my $curOutDir = "$tmp/closed-member/";
sub MFnext {}
sub loop2C_check {}
for my $case (
    [5,'skipped_empty_input',1], [5,'skipped_too_small',1],
    [5,'skipped_cleaned_empty',1], [5,'completed',0],
    [5,'skipped_sdm_warning',0], [2,'skipped_empty_input',0],
) {
    my ($mode,$status,$count) = @$case;
    $MFopt{DoAssembly} = $mode;
    $AsGrps{group}{CntPreAssNoPrim} = 2;
    $AsGrps{group}{CntAimAss} = 2;
    $AsGrps{group}{ClosedCompleted} = [];
    $closedSample = {outcome=>{status=>$status},components=>{}};
    eval "for (1) { $closed_body }"; die $@ if $@;
    is($AsGrps{group}{CntPreAssNoPrim}, 2+$count,
        "assembly mode $mode, $status: completed-member accounting is scoped correctly");
    # only completed members of a shared group registered no reads they should have
    is_deeply($AsGrps{group}{ClosedCompleted}, $status eq 'completed' ? [$curOutDir] : [],
        "assembly mode $mode, $status: fast-path member is recorded for the subset guard only when completed");
}

done_testing();
