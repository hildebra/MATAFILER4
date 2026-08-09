use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use IO::Compress::Gzip qw(gzip $GzipError);
use lib File::Spec->catdir($Bin, '..');
use Mods::GenoMetaAss qw(resetAsGrps contig_stats_coverage_complete);
use Mods::IO_Tamoc_progs qw(inputFmtSpades);
use Mods::WorkflowControl qw(
	advance_loop_window overlap_loop_window rolling_completed_frontier rolling_loop_transition priority_outputs_complete parse_loop_spec should_rerun_locked_window assembly_cores_for_input assembly_group_output_dirs balanced_parallel_batches hybrid_group_ready
	hybrid_package_complete hybrid_package_sample_id hybrid_local_scratch_gb missing_input_files source_input_files parse_ignored_samples
	sample_base_output_dir sample_is_ignored workflow_members_match
	normalise_job_dependencies append_job_dependencies deferred_command_dependencies augment_deferred_submission
	commands_are_lightweight_filesystem input_terminal_status cleaned_primary_libraries_empty sample_empty_marker_reason reconcile_sample_empty_marker cleanup_stage_barrier
);

is(assembly_cores_for_input(input_mb => 0, configured_cores => 0), 8,
	'automatic assembly cores use the minimum for empty or tiny groups');
is(assembly_cores_for_input(input_mb => 500, configured_cores => 0), 8,
	'automatic assembly cores remain at the minimum through 500 MiB');
is(assembly_cores_for_input(input_mb => 5370, configured_cores => 0), 28,
	'automatic assembly cores scale linearly through the input range');
is(assembly_cores_for_input(input_mb => 10 * 1024, configured_cores => 0), 48,
	'automatic assembly cores reach the maximum at 10 GiB');
is(assembly_cores_for_input(input_mb => 50 * 1024, configured_cores => 0), 48,
	'automatic assembly cores are capped for very large groups');
is(assembly_cores_for_input(input_mb => 50 * 1024, configured_cores => 24), 24,
	'an explicit assembly-core setting overrides automatic scaling');

is(normalise_job_dependencies('run12;;run7', ['run7', '', 'run3;run12']),
	'run12;run7;run3', 'job dependencies are flattened, deduplicated, and stable');
my $job_dependencies = 'run12;';
is(append_job_dependencies(\$job_dependencies, 'run7', 'run12'), 'run12;run7',
	'job dependencies are appended without empty or duplicate scheduler entries');
is(deferred_command_dependencies(
	dependencies => 'run12;run7', submitted => ['run20', 'run21'],
), 'run12;run7',
	'independent deferred commands share only their producer barrier by default');
is(deferred_command_dependencies(
	dependencies => 'run12;run7', submitted => ['run20', 'run21'],
	chain_previous => 1,
), 'run12;run7;run21',
	'genuine multi-stage deferred commands explicitly depend on their predecessor');
ok(commands_are_lightweight_filesystem(
	"rm -rf /tmp/sample/rawRds\nmkdir -p /tmp/sample/rawRds;\n"
	."sleep 1\nln -s /reads/a.fq.gz /tmp/sample/rawRds/a.fq.gz\n"
	."touch /tmp/sample/rawRds/done.sto\n"
), 'lightweight UZ filesystem setup is safe to execute locally');
ok(!commands_are_lightweight_filesystem(
	"mkdir -p /tmp/sample/rawRds\npigz -d -c /reads/a.fq.gz > /tmp/sample/rawRds/a.fq\n"
), 'decompression and redirection keep UZ work on the scheduler');
ok(!commands_are_lightweight_filesystem(
	"ln -s /reads/a.fq.gz /tmp/sample/rawRds/a.fq.gz && touch /tmp/sample/rawRds/done.sto\n"
), 'compound shell expressions are not classified as lightweight');

is(input_terminal_status(input_size_mb => 0, threshold_mb => 1),
	'skipped_empty_input', 'zero input has the explicit empty terminal status');
is(input_terminal_status(input_size_mb => 0.5, threshold_mb => 1),
	'skipped_too_small', 'nonzero input below the threshold is explicitly too small');
is(input_terminal_status(input_size_mb => 8272.7, threshold_mb => 1), '',
	'multi-gigabyte input cannot retain an empty or too-small outcome');
my $empty_state_root = tempdir(CLEANUP => 1);
my $empty_marker = File::Spec->catfile($empty_state_root, 'SMPL.empty');
open(my $empty_marker_fh, '>', $empty_marker) or die "Cannot create $empty_marker: $!";
close($empty_marker_fh) or die "Cannot close $empty_marker: $!";
is(reconcile_sample_empty_marker(
	sample_root => $empty_state_root, input_size_mb => 8272.7, threshold_mb => 1,
), 0, 'authoritative nonempty input clears the cached empty state');
ok(!-e $empty_marker, 'stale SMPL.empty is physically removed');
open($empty_marker_fh, '>', $empty_marker) or die "Cannot recreate $empty_marker: $!";
close($empty_marker_fh) or die "Cannot close $empty_marker: $!";
is(reconcile_sample_empty_marker(
	sample_root => $empty_state_root, input_size_mb => 0.5, threshold_mb => 1,
), 1, 'a genuinely too-small sample retains its empty marker');
ok(-e $empty_marker, 'valid SMPL.empty remains present for a too-small sample');
open($empty_marker_fh, '>', $empty_marker) or die "Cannot rewrite $empty_marker: $!";
print {$empty_marker_fh} "cleaned_primary_reads_empty\n";
close($empty_marker_fh) or die "Cannot close $empty_marker: $!";
is(sample_empty_marker_reason($empty_state_root), 'cleaned_primary_reads_empty',
	'a reasoned cleaner marker records that no primary reads survived');
is(reconcile_sample_empty_marker(
	sample_root => $empty_state_root, input_size_mb => 8272.7, threshold_mb => 1,
), 1, 'a cleaned-empty marker remains authoritative despite large raw input');
my $cleaned_r1 = File::Spec->catfile($empty_state_root, 'filtered.1.fq.gz');
my $cleaned_r2 = File::Spec->catfile($empty_state_root, 'filtered.2.fq.gz');
my $cleaned_single = File::Spec->catfile($empty_state_root, 'filtered.singl.fq.gz');
my $empty_fastq = '';
gzip \$empty_fastq => $cleaned_r1 or die "Cannot write $cleaned_r1: $GzipError";
gzip \$empty_fastq => $cleaned_r2 or die "Cannot write $cleaned_r2: $GzipError";
gzip \$empty_fastq => $cleaned_single or die "Cannot write $cleaned_single: $GzipError";
my $cleaned_libraries = [{
	files => {r1 => $cleaned_r1, r2 => $cleaned_r2, single => $cleaned_single},
}];
ok(cleaned_primary_libraries_empty($cleaned_libraries),
	'header-only gzip cleaner outputs are recognized as empty primary FASTQ libraries');
unlink $cleaned_single or die "Cannot replace $cleaned_single: $!";
my $one_fastq_record = "\@read\nACGT\n+\n!!!!\n";
gzip \$one_fastq_record => $cleaned_single
	or die "Cannot write $cleaned_single: $GzipError";
ok(!cleaned_primary_libraries_empty($cleaned_libraries),
	'a single surviving FASTQ record keeps a completed cleaner output eligible for assembly');
ok(-e $empty_marker, 'cleaned-empty marker is retained for the next submission pass');
my $cleanup_barrier = cleanup_stage_barrier(
	{name => 'contig stats', required => 1, complete => 1},
	{name => 'binning', required => 1, complete => 0, dependencies => 'run9;run9'},
	{name => 'consSNP', required => 1, complete => 0, dependencies => ''},
	{name => 'disabled stage', required => 0, complete => 0},
);
ok(!$cleanup_barrier->{ready}, 'cleanup waits when a required terminal stage is neither complete nor scheduled');
is_deeply($cleanup_barrier->{blocked}, ['consSNP'], 'cleanup barrier reports the unsatisfied terminal stage');
is($cleanup_barrier->{dependencies}, 'run9', 'cleanup barrier retains scheduled terminal dependencies');
$cleanup_barrier = cleanup_stage_barrier(
	{name => 'contig stats', required => 1, complete => 0, dependencies => 'run4'},
	{name => 'binning', required => 1, complete => 1},
	{name => 'consSNP', required => 1, complete => 0, dependencies => ['run8', 'run4']},
);
ok($cleanup_barrier->{ready}, 'cleanup may follow every incomplete stage when each has a terminal job');
is($cleanup_barrier->{dependencies}, 'run4;run8', 'terminal cleanup dependencies are normalized and deduplicated');
is(hybrid_local_scratch_gb(
	assembler_gb => 150, preassembly_bytes => 512 * 1024**2, max_synthetic_depth => 20,
), 160, 'hybrid scratch adds estimated synthetic reads to assembler workspace');
is_deeply(balanced_parallel_batches(21, 20), [11, 10],
	'synthetic simulator batches avoid a nearly empty final wave');
is_deeply(balanced_parallel_batches(45, 20), [15, 15, 15],
	'synthetic simulator batches are balanced and capped below the concurrency limit');
is_deeply(balanced_parallel_batches(40, 20), [20, 20],
	'exact concurrency multiples retain full balanced waves');

my $slurm = augment_deferred_submission(
	qmode => 'slurm', run_tag => 'run', dependencies => 'run12;run7',
	command => 'sbatch /tmp/map.sh',
	script => "#!/bin/bash\n#SBATCH --dependency=afterok:3:12\necho map\n",
);
like($slurm->{script}, qr/^#SBATCH --dependency=afterok:3:12:7$/m,
	'deferred Slurm submission preserves existing dependencies and adds final group jobs');
my $sge = augment_deferred_submission(
	qmode => 'sge', run_tag => 'run', dependencies => 'run12;run7',
	command => 'qsub -hold_jid 3,12 /tmp/map.sh', script => '',
);
like($sge->{command}, qr/-hold_jid 3,12,7\b/,
	'deferred SGE submission merges final dependencies into hold_jid');
my $lsf = augment_deferred_submission(
	qmode => 'lsf', run_tag => 'run', dependencies => 'run12;run7',
	command => 'bsub -w "done(3) && done(12)" < /tmp/map.sh', script => '',
);
like($lsf->{command}, qr/-w "done\(3\) && done\(12\) && done\(7\)"/,
	'deferred LSF submission merges final dependencies into its wait expression');

sub write_file {
	my ($path, $contents) = @_;
	(my $dir = $path) =~ s{/[^/]+$}{};
	make_path($dir);
	open my $fh, '>', $path or die "Cannot write $path: $!";
	print {$fh} $contents;
	close $fh;
}

is(sample_base_output_dir('/runs/S.1/archive/S.1/', 'S.1'), '/runs/S.1/archive/',
	'base output removes only the final literal sample directory');
eval { sample_base_output_dir('/runs/S10/', 'S1') };
like($@, qr/does not end in sample directory/, 'base output rejects a non-matching suffix');

my $ignored = parse_ignored_samples('S1, S.2');
ok(sample_is_ignored($ignored, 'S1'), 'exact ignored sample is selected');
ok(sample_is_ignored($ignored, 'S.2'), 'regex characters are treated literally');
ok(!sample_is_ignored($ignored, 'S10'), 'ignore matching does not select prefixes');

ok(workflow_members_match(['/run/A/', '/run/B'], ['/run/B/', '/run/A']),
	'assembly members compare exactly independent of ordering and trailing slash');
ok(!workflow_members_match(['/run/A', '/run/B'], ['/run/A', '/run/C']),
	'same-sized group with a replaced member is rejected');
ok(!workflow_members_match(['/run/A'], ['/run/A', '/run/B']),
	'group with an extra member is rejected');

my %map = (
	opt => { smpl_order => ['A', 'B', 'C'] },
	A => { AssGroup => 'gut', wrdir => '/run/A/' },
	B => { AssGroup => 'gut', wrdir => '/run/B/' },
	C => { AssGroup => '-1', wrdir => '/run/C/' },
);
is_deeply(assembly_group_output_dirs(\%map, 'gut'), ['/run/A/', '/run/B/'],
	'expected group membership comes from the complete map rather than loop progress');
is_deeply(assembly_group_output_dirs(\%map, 'C'), ['/run/C/'],
	'single-sample implicit assembly group is supported');

my $window = advance_loop_window(
	from => 0, to => 250, upper => 600, window_size => 250, initial_loops => 6,
);
is_deeply($window, {
	from => 250, to => 500, loop_count => 6, reset_index => 249, has_window => 1,
}, 'first completed window advances directly to the next disjoint range');
$window = advance_loop_window(
	from => $window->{from}, to => $window->{to}, upper => 600,
	window_size => 250, initial_loops => 6,
);
is_deeply($window, {
	from => 500, to => 600, loop_count => 6, reset_index => 499, has_window => 1,
}, 'final partial window starts at the previous upper boundary');
$window = advance_loop_window(
	from => $window->{from}, to => $window->{to}, upper => 600,
	window_size => 250, initial_loops => 6,
);
is_deeply($window, {
	from => 600, to => 600, loop_count => 0, reset_index => 599, has_window => 0,
}, 'window loop terminates without returning to earlier samples');

my $transition = rolling_loop_transition(
	from => 500, to => 600, upper => 600, loop_count => 3, window_size => 250,
);
is_deeply($transition, {action => 'repeat', to => 600},
	'a nonempty final partial window retains its remaining loop passes');
$transition = rolling_loop_transition(
	from => 600, to => 600, upper => 600, loop_count => 3, window_size => 250,
);
is_deeply($transition, {action => 'verify', to => 600},
	'an exhausted frontier cannot repeat an empty range even with loop passes remaining');
$transition = rolling_loop_transition(
	from => 250, to => 500, upper => 600, loop_count => 0, window_size => 250,
);
is_deeply($transition, {action => 'expand', to => 600},
	'a completed rolling range admits the capped final partial window');
$transition = rolling_loop_transition(
	from => 500, to => 600, upper => 600, loop_count => 0, window_size => 250,
);
is_deeply($transition, {action => 'verify', to => 600},
	'an exhausted nonempty final window proceeds to full-range verification');

my $overlap = overlap_loop_window(
	to => 250, upper => 600, window_size => 250, submitted_jobs => 12,
);
is_deeply($overlap, {to => 500, extended => 1, job_limit => 62},
	'a lightly loaded pass admits the next block before waiting');
$overlap = overlap_loop_window(
	to => 250, upper => 600, window_size => 250, submitted_jobs => 63,
);
is_deeply($overlap, {to => 250, extended => 0, job_limit => 62},
	'a normally loaded pass retains the current block boundary');
$overlap = overlap_loop_window(
	to => 250, upper => 600, window_size => 250, submitted_jobs => 1000,
	last_pass => 1,
);
is_deeply($overlap, {to => 500, extended => 1, job_limit => 62},
	'the final pass admits the next block even when it is not lightly loaded');
$overlap = overlap_loop_window(
	to => 250, upper => 600, window_size => 250, submitted_jobs => 0,
	last_pass => 1,
);
is_deeply($overlap, {to => 500, extended => 1, job_limit => 62},
	'an idle final pass merges the next block instead of closing the window first');
$overlap = overlap_loop_window(
	to => 500, upper => 600, window_size => 250, submitted_jobs => 1,
);
is_deeply($overlap, {to => 600, extended => 1, job_limit => 62},
	'an overlap extension is capped at the selected sample upper bound');
$overlap = overlap_loop_window(
	to => 250, upper => 600, window_size => 250, submitted_jobs => 0,
);
ok(!$overlap->{extended},
	'an idle non-final pass retains the normal completed-block transition');
$overlap = overlap_loop_window(
	to => 250, upper => 600, window_size => 250, submitted_jobs => 1,
	already_extended => 1, last_pass => 1,
);
ok(!$overlap->{extended},
	'only one new block is admitted per submission pass, including a final pass');
$overlap = overlap_loop_window(
	to => 600, upper => 600, window_size => 250, submitted_jobs => 1000,
	last_pass => 1,
);
ok(!$overlap->{extended}, 'the final pass cannot extend beyond the selected range');

my $rolling = rolling_completed_frontier(
	from => 0, upper => 10, window_size => 4,
	completed => { 0 => 1, 1 => 1, 3 => 1 },
);
is_deeply($rolling, { from => 2, to => 6, advanced => 2, finished => 0 },
	'the frontier advances only across the continuously completed prefix');
$rolling = rolling_completed_frontier(
	from => 2, upper => 10, window_size => 4,
	completed => { 2 => 1, 3 => 1, 4 => 1, 7 => 1 },
);
is_deeply($rolling, { from => 5, to => 9, advanced => 3, finished => 0 },
	'completed slots replenish the rolling range without skipping a gap');
$rolling = rolling_completed_frontier(
	from => 8, upper => 10, window_size => 4,
	completed => { 8 => 1, 9 => 1 },
);
ok($rolling->{finished} && $rolling->{from} == 10,
	'a fully completed prefix reaches the selected upper bound');

my @priority_visits;
my %priority_present = map { $_ => 1 } qw(mapping depth-2 assembly);
my $priority = priority_outputs_complete([
	{ name => 'mapping', required => 1, all => ['mapping'] },
	{ name => 'depth', required => 1, any => ['depth-1', 'depth-2'] },
	{ name => 'assembly', required => 1, all => ['assembly'] },
	{ name => 'binning', required => 0, all => ['binning'] },
], sub {
	my ($path) = @_;
	push @priority_visits, $path;
	return $priority_present{$path};
});
ok($priority->{complete}, 'ordered priority probing accepts every requested stage');
is_deeply(\@priority_visits, [qw(mapping depth-1 depth-2 assembly)],
	'priority probing checks stages in order and stops an alternative group once found');
@priority_visits = ();
delete $priority_present{mapping};
$priority = priority_outputs_complete([
	{ name => 'mapping', required => 1, all => ['mapping'] },
	{ name => 'assembly', required => 1, all => ['assembly'] },
], sub { push @priority_visits, $_[0]; return $priority_present{$_[0]}; });
is($priority->{missing_stage}, 'mapping', 'the first missing priority stage is reported');
is_deeply(\@priority_visits, ['mapping'],
	'a missing high-priority output prevents lower-priority filesystem checks');

is_deeply(parse_loop_spec('6:250'), {loop_count => 6, window_size => 250},
	'a loop window specification is parsed exactly');
is_deeply(parse_loop_spec('1'), {loop_count => 6, window_size => 0},
	'the historical boolean-style loop switch keeps its six-pass behavior');
eval { parse_loop_spec('0:250') };
like($@, qr/requires positive X:Y values/,
	'a zero loop count cannot silently truncate processing to the first window');
eval { parse_loop_spec('6:250garbage') };
like($@, qr/Invalid -loopTillComplete value/,
	'a partially matching loop specification is rejected');

ok(should_rerun_locked_window(
	active_jobs => 9, sample_count => 1000, remove_locks => 0,
), 'fewer than one percent active jobs reruns a retained-lock window');
ok(!should_rerun_locked_window(
	active_jobs => 10, sample_count => 1000, remove_locks => 0,
), 'exactly one percent active jobs does not satisfy the strict threshold');
ok(should_rerun_locked_window(
	active_jobs => 2, sample_count => 100, remove_locks => 0,
), 'fewer than the default active-job threshold reruns even when one percent is smaller');
ok(should_rerun_locked_window(
	active_jobs => 3, sample_count => 100, remove_locks => 0,
), 'the default threshold includes three active jobs');
ok(!should_rerun_locked_window(
	active_jobs => 3, active_job_threshold => 2,
	sample_count => 100, remove_locks => 0,
), 'a configured active-job threshold overrides the default');
ok(!should_rerun_locked_window(
	active_jobs => 0, sample_count => 1000, remove_locks => 0,
), 'zero active jobs cannot cause a retained-lock retry loop');
ok(!should_rerun_locked_window(
	active_jobs => 1, sample_count => 1000, remove_locks => 1,
), 'lock-removal mode does not use the retained-lock retry policy');

my $root = tempdir(CLEANUP => 1);
$root =~ s{\\}{/}g;
my $coverage_dir = "$root/ContigStats";
for my $suffix (qw(percontig median.percontig pergene count_pergene)) {
	write_file("$coverage_dir/Coverage.$suffix.gz", "$suffix\n");
}
write_file("$coverage_dir/Coverage.stone", "done\n");
ok(contig_stats_coverage_complete($coverage_dir, 'Coverage'),
	'planner and worker share a complete ContigStats coverage contract');
unlink "$coverage_dir/Coverage.percontig.gz";
ok(!contig_stats_coverage_complete($coverage_dir, 'Coverage'),
	'a completion stone cannot hide a missing coverage derivative');

my $package = "$root/preAssmblGrp_gut";
write_file("$package/scaffolds.fasta.filt", ">contig\nACGT\n");
write_file("$package/Coverage.percontig.gz", "coverage\n");
write_file("$package/Coverage.median.percontig", "median\n");
write_file("$package/mapping.coverage.gz", "window coverage\n");
write_file("$package/breakpoints.tsv.gz", "compressed breakpoint fixture\n");
write_file("$package/package.manifest.tsv", "key\tvalue\nschema_version\t2\nsample_id\tSample.A\n");
write_file("$package/moved.sto", "done\n");
ok(hybrid_package_complete($package),
	'versioned hybrid package requires all package-local outputs');
is(hybrid_package_sample_id($package), 'Sample.A',
	'synthetic-read naming uses the sample recorded in its package');
unlink "$package/Coverage.percontig.gz";
ok(!hybrid_package_complete($package), 'hybrid package with a missing output is incomplete');

ok(hybrid_group_ready(2, 1, 3),
	'no-primary member counts toward final hybrid group readiness');
ok(!hybrid_group_ready(1, 1, 3), 'hybrid group waits for remaining packages');

my $valid_input = "$root/input.1.fastq";
my $empty_input = "$root/input.2.fastq";
write_file($valid_input, "reads\n");
write_file($empty_input, '');
is_deeply(missing_input_files($valid_input), [], 'nonempty input passes validation');
is_deeply(missing_input_files($valid_input, $empty_input, "$root/missing.fastq"),
	[$empty_input, "$root/missing.fastq"],
	'all empty and missing libraries are reported, not only the first');
my $source_dir = "$root/source";
make_path($source_dir);
my $source_read = "$source_dir/read.1.fq.gz";
write_file($source_read, "reads\n");
my $resolved_sources = source_input_files($source_dir, 'read.1.fq.gz', '/external/read.2.fq.gz');
is_deeply($resolved_sources, [$source_read, '/external/read.2.fq.gz'],
	'retry validation resolves original relative inputs without rewriting absolute inputs');
is_deeply(missing_input_files($resolved_sources->[0]), [],
	'an existing source remains valid when its generated rawRds destination is absent');

my $spades_inputs = inputFmtSpades(
	['A.1.fq.gz', 'B.1.fq.gz'], ['A.2.fq.gz', 'B.2.fq.gz'],
	['A.s.fq.gz', 'B.s.fq.gz'], $root, ['ill', 'ill'],
);
is(scalar(() = $spades_inputs =~ /--pe1-1/g), 2,
	'metaSPAdes coassembly keeps all left-read files in one paired library');
is(scalar(() = $spades_inputs =~ /--pe1-2/g), 2,
	'metaSPAdes coassembly keeps all right-read files in one paired library');
unlike($spades_inputs, qr/--pe2-/, 'coassembly does not create an unsupported second metaSPAdes library');

done_testing;
