use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::GenoMetaAss qw(resetAsGrps contig_stats_coverage_complete);
use Mods::IO_Tamoc_progs qw(inputFmtSpades);
use Mods::WorkflowControl qw(
	advance_loop_window overlap_loop_window rolling_completed_frontier rolling_loop_transition priority_outputs_complete parse_loop_spec should_rerun_locked_window assembly_cores_for_input assembly_group_output_dirs balanced_parallel_batches hybrid_group_ready
	hybrid_package_complete hybrid_package_sample_id hybrid_local_scratch_gb missing_input_files source_input_files parse_ignored_samples
	sample_base_output_dir sample_is_ignored workflow_members_match
	normalise_job_dependencies append_job_dependencies deferred_command_dependencies augment_deferred_submission
	commands_are_lightweight_filesystem cleanup_stage_barrier
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

open my $source, '<', File::Spec->catfile($Bin, '..', 'MATAF4.pl')
	or die "Cannot inspect MATAF4.pl: $!";
my $mataf4 = do { local $/; <$source> };
close $source;
like($mataf4, qr/sub spadesAssembly.*?\$cmd \.= " --meta ".*?qsubSystem\([^\n]+\$nCores[^\n]+\$defTotMem/s,
	'SPAdes always uses metagenome mode and requests the resources used by its command');
unlike($mataf4, qr/sub spadesAssembly.*?\$cmd \.= " --sc ".*?sub movePreAssmData/s,
	'SPAdes no longer substitutes single-cell mode for metagenomic coassembly');
like($mataf4, qr/--nano-raw "\.join\(" ", \@inRds\)/,
	'Flye receives every long-read file collected for the assembly group');
like($mataf4, qr/\$primaryInputSize \*= \$totalSmpls \/ \$knownSmpls/,
	'missing group-member sizes are extrapolated in the correct direction');
like($mataf4, qr/\$baseOut\s*=\s*sample_base_output_dir\(\$curOutDir,\s*\$curSmpl\)/,
	'normal execution uses safe sample base-output derivation');
like($mataf4, qr/sample_is_ignored\(\$ignoredSamplesHR,\s*\$SmplName\)/,
	'normal execution uses exact ignored-sample lookup');
like($mataf4, qr/my \@missingInputs\s*=\s*\@\{missing_input_files\(\@sourceInputs\)\}/,
	'normal execution validates authoritative inputs instead of generated rawRds destinations');
like($mataf4, qr/if \(hybrid_package_complete\(\$mvD\)\)/,
	'hybrid execution uses package-local completeness');
like($mataf4, qr/sub movePreAssmData\s*\{.*?\$mvD\s*=~\s*s\{\/\+\$\}\{\};.*?my \$stage\s*=\s*"\$mvD\.stage\.\$packageTag"/s,
	'preassembly packaging strips trailing slashes before creating sibling staging directories');
like($mataf4, qr/SupportReads\}\s*=~\s*m\/\(\?:PB\|ONT\):\//,
	'hybrid preassembly recognizes both supported long-read technologies');
like($mataf4, qr/long_reads_detected.*?\(\?:PB\|ONT\)/s,
	'final hybrid assembly applies the same PB/ONT technology contract');
like($mataf4, qr/mixes ONT and PacBio reads/,
	'mixed ONT/PacBio hybrid groups are rejected explicitly');
like($mataf4, qr/sim_failed=0.*?wait \\\$sim_pid_/s,
	'each background synthetic-read simulator is waited on and validated');
like($mataf4, qr/balanced_parallel_batches\(scalar\(\@simulatorCommands\), 20\)/,
	'synthetic-read simulator command waves are capped at twenty processes');
is(scalar(() = $mataf4 =~ /--length-template/g), 1,
	'hybrid command construction emits one length-template argument');
like($mataf4, qr/\$packageSample\.synthetic\.fastq\.gz/,
	'synthetic read filenames identify their source package sample');
like($mataf4,
	qr/hybrid_local_scratch_gb\(.*?assembler_gb.*?preassembly_bytes.*?max_synthetic_depth/s,
	'metaMDBG local scratch combines assembler space with synthetic-read space');
like($mataf4,
	qr/sub spadesAssembly.*?\$QSBoptHR->\{tmpSpace\} = \$locDiskSpace;.*?if \(\$hostFilter \|\| \$MFopt\{SpadesAlwaysHDDnode\}\).*?qsubSystem.*?else \{.*?qsubSystem.*?\$QSBoptHR->\{tmpSpace\} = \$tmpSHDD/s,
	'both SPAdes host-selection branches request and then restore assembly scratch');
like($mataf4,
	qr/sub longRdAssembly.*?my \$assemblerScratchGB = \$HDDspace\{assembler\}.*?if \(\$nameProg eq "mMDBG"\).*?hybrid_local_scratch_gb.*?\$QSBoptHR->\{tmpSpace\} = \$assemblerScratchGB\."G".*?qsubSystem/s,
	'Flye and metaMDBG submissions explicitly request their node-local assembly space');
like($mataf4,
	qr/sub detectRibo.*?\$curSHFF = \$predefSHDD if \(\$curSHFF < \$predefSHDD\).*?\$QSBoptHR->\{tmpSpace\}= \$curSHFF \. "G".*?RiboFinder\.sh/s,
	'ribosomal extraction scratch keeps its configured floor while scaling with input');
like($mataf4,
	qr/sub metphlanMapping.*?inputFileSizeMB\} \* 6.*?qsubSystem\(\$qsubFile.*?\$QSBoptHR->\{tmpSpace\} = \$previousTmpSpace/s,
	'MetaPhlAn requests input-scaled scratch for its uncompressed SAM');
like($mataf4,
	qr/sub mOTU2Mapping.*?inputFileSizeMB\} \* 4.*?qsubSystem\(\$logDir\."mOTU2_prof\.sh".*?\$QSBoptHR->\{tmpSpace\} = \$previousTmpSpace/s,
	'mOTUs requests input-scaled scratch and restores the submission default');
like($mataf4,
	qr/sub krakenTaxEst.*?my \$requestedTmpSpace = \$HDDspace\{kraken\}.*?inputFileSizeMB\} \* 6.*?qsubSystem\(\$logDir\."KrkTax\.sh".*?\$QSBoptHR->\{tmpSpace\} = \$previousTmpSpace/s,
	'Kraken taxonomy requests space for raw and translated local intermediates');
unlike($mataf4, qr/zcat \$nodeTmp\/contigs\.fasta\.gz/,
	'metaMDBG recovery no longer uses an unchecked zcat/remove chain');
like($mataf4,
	qr/test -s \$nodeTmp\/contigs\.fasta\.gz.*?\$pigzBin -dc.*?test -s \$nodeTmp\/scaffolds\.fasta/s,
	'metaMDBG compressed output is validated before and after normalization');
like($mataf4, qr/HybridAssemblyComparison\.tsv/,
	'hybrid finalization requires a comparative assembly report');
like($mataf4,
	qr/my \$comparisonPreassembly = \$hybridPreassemblies\[0\];.*?--preassembly \$comparisonPreassembly --output/s,
	'hybrid comparison receives exactly one deterministic preassembly');
unlike($mataf4,
	qr/join\(" ", map \{ "--preassembly \$_" \} \@hybridPreassemblies\)/,
	'hybrid comparison does not repeat the preassembly option');
like($mataf4, qr{mapping/\$SmplN-smd\.bam\.breakpoints\.tsv\.gz},
	'metagStats reads the sample breakpoint report from its mapping directory');
like($mataf4, qr/BreakpointCount.*BreakpointContigs.*BreakpointBases.*BreakpointMeanLength.*BreakpointMaxLength/,
	'metagStats exposes the important per-sample breakpoint statistics');
unlike($mataf4, qr/die "Deleting previous results/,
	'unfinished-result repair no longer stops at a debug die');
unlike($mataf4, qr/die "now recalc/,
	'redo-failed repair no longer stops at a debug die');
like($mataf4,
	qr/CntAimAss\}\s*>\s*1\s*&&\s*!\$MFconfig\{OKtoRWassGrps\}.*?Refusing to rebuild shared assembly group/s,
	'shared assembly rebuild retains the established destructive authorization gate');
like($mataf4, qr/Submitting deferred assembly-group mapping jobs.*?postSubmQsub/s,
	'assembly-group mappings deferred before the final member are submitted after assembly scheduling');
like($mataf4,
	qr/deferred_command_dependencies\(.*?chain_previous\s*=>\s*\$options->\{chain_previous\}/s,
	'deferred command serialization requires an explicit caller option');
like($mataf4,
	qr/sub postSubmQsub.*?\{submittedJobs\}\+\+.*?\{slurmDependencySubmittedAt\}\{\$scheduler_job_id\}\s*=\s*time/s,
	'accepted deferred submissions update loop and Slurm dependency bookkeeping');
like($mataf4,
	qr/sub postSubmQsub.*?recordSampleLockJobs\(\s*\$QSBoptHR->\{LOCKfile\}, \[\$scheduler_job_id\]/s,
	'accepted deferred submissions record their scheduler ID in the sample lock ledger');
like($mataf4,
	qr/sampleLockActiveJobs\(\$smplLockF, \$QSBoptHR\).*?unlink \$smplLockF/s,
	'MATAF4 releases a sample lock only after its recorded jobs leave the scheduler');
like($mataf4,
	qr/my \$lightweightLocal = commands_are_lightweight_filesystem\(\$unzipcmd\).*?systemW \$unzipcmd.*?qsubSystem\(\$logDir\."UNZP\.sh"/s,
	'lightweight UZ setup runs locally while data-processing commands remain scheduled');
unlike($mataf4,
	qr/\$unzipcmd \.= .*?mkdir -p \$finDest\/rawRds\/;\\nsleep 1;/,
	'lightweight UZ setup no longer adds an unconditional one-second delay');
like($mataf4,
	qr/if \(\$MFconfig\{maxUnzpJobs\} > 0 && \@unzipjobs >= \$MFconfig\{maxUnzpJobs\}\).*?\$unzipjobs\[-\(\$MFconfig\{maxUnzpJobs\}\)\]/s,
	'heavy UZ submissions use the configured rolling concurrency cap');
like($mataf4,
	qr/if \(\$MFconfig\{precheckInputDirs\}\) \{.*?populateInputSizesFast\(\$_\) for \@samples;.*?if \(\$MFconfig\{inspectState\}\).*?for \(\$JNUM=\$from/s,
	'full input-directory sizing is opt-in before inspection or sample processing');
like($mataf4,
	qr/\$MFconfig\{precheckInputDirs\}\s*=\s*0;.*?"precheckInputDirs=i"\s*=>\s*\\\$MFconfig\{precheckInputDirs\}/s,
	'the expensive all-sample input precheck is disabled by default and exposed explicitly');
like($mataf4,
	qr/my %runReport = \(.*?samples => \{\}.*?empty_samples => \{\}.*?present_assemblies => 0.*?my %d2Inputs = \(.*?filtered_read1 => \[\].*?dependencies => ""/s,
	'run reporting and cross-sample distance state are grouped by responsibility');
like($mataf4,
	qr/my \$assemblyOutputsRequired = \$MFopt\{DoAssembly\} \? 1 : 0;.*?my \$assemblyWorkflowComplete = !\$assemblyOutputsRequired \|\| \$boolAssemblyOK/s,
	'assembly-independent completion does not require an assembly checkpoint');
like($mataf4,
	qr/my \$terminalOutputsComplete = !\$assemblyOutputsRequired \|\| \(.*?\$efinAssLoc.*?\$cleanupContigStatsComplete/s,
	'terminal assembly outputs are required only when assembly is enabled');
like($mataf4,
	qr/name => 'final assembly publication', required => \$assemblyOutputsRequired.*?name => 'contig stats', required => \$assemblyOutputsRequired/s,
	'assembly-independent cleanup skips assembly and ContigStats barriers');
like($mataf4,
	qr/sub postprocess.*?\$runReport\{present_assemblies\} = 0.*?\$closedSample->\{present_assembly\}/s,
	"only sentinel-confirmed assemblies increment the final report counter");
like($mataf4,
	qr/if \(\$MFopt\{DoAssembly\} && \$runReport\{present_assemblies\} > 0/s,
	'gene-catalog suggestions are disabled for assembly-independent runs');
like($mataf4,
	qr/binning requires -assembleMG.*?!\$MFopt\{DoAssembly\} && \$MFopt\{DoMetaBat2\}.*?assembly consensus SNP calling requires -assembleMG.*?!\$MFopt\{DoAssembly\}.*?structural-variant calling requires -assembleMG/s,
	'assembly-only binning and variant requests fail early without assembly');
unlike($mataf4, qr/^my (?:%jmp|%MFstats|\$baseDir|\$mvCmd)\b/m,
	'confirmed dead global variables are removed');
unlike($mataf4, qr/^sub (?:check_sdm_loc|setupInput|fastCovCalc)\b/m,
	'uncalled empty or explicitly defunct routines are removed');
like($mataf4,
	qr/foreach my \$sampleName \(keys %\{\$d2Inputs\{samples\}\}\).*?print O "\$sampleName\\t\$inputRawFQs\{\$sampleName\}/s,
	'raw-input reporting retains the sample-name key instead of substituting its boolean value');
like($mataf4,
	qr/sub discoverSampleInputs.*?my \$cached = \$map\{\$sample\}\{inputDiscovery\};.*?return \$cached.*?doDateFileCheck.*?SupportReads may specify one directory or a list of files.*?Unsupported SupportReads format/s,
	'one cached discovery covers date-filtered primary reads and both support-input forms');
like($mataf4,
	qr/sub populateInputSizesFast.*?discoverSampleInputs\(\$sample\).*?primary_bytes.*?support_bytes.*?sub spaceInAssGrp/s,
	'fast input sizing uses the shared discovery result');
like($mataf4,
	qr/Suppl: %.1f Mb.*?\{inputXFileSizeMB\}.*?if \(\@paX1 \|\| \@paXs \|\| \@paBamX\)/s,
	'input summary prints supplementary size only for discovered support reads');
unlike($mataf4,
	qr/Suppl: %.1f Mb", \$map\{\$curSmpl\}\{inputFileSizeMB\}/,
	'input summary cannot substitute the primary size for supplementary input');
like($mataf4,
	qr/Primary and supplementary inputs resolve to the same file\(s\).*?\@overlap/s,
	'input discovery rejects physical files assigned to both read scopes');
like($mataf4,
	qr/sub seedUnzip2tmp.*?discoverSampleInputs\(\$curSmpl, \$fastp\).*?\@pa1 = \@\{\$inputDiscovery->\{primary\}\{read1\}\}.*?\@paX1 = \@\{\$inputDiscovery->\{support\}\{read1\}\}/s,
	'input staging reuses the cached primary and support file selections');
my ($seed_unzip_source) = $mataf4 =~ /(sub seedUnzip2tmp\{.*?)(?=\nsub \w)/s;
ok(defined($seed_unzip_source), 'seedUnzip2tmp source can be isolated');
unlike($seed_unzip_source || "", qr/\b(?:discoverReadFiles|parseSupportReads)\s*\(/,
	'input staging contains no duplicate file-discovery implementation');
like($mataf4, qr/#4\.14:.*?cache one validated input discovery.*?#4\.15:.*?#4\.16:.*?#4\.17:.*?#4\.18:.*?#4\.19:.*?#4\.20:.*?#4\.21:.*?#4\.22:.*?#4\.23:.*?#4\.24:.*?#4\.25:.*?#4\.26:.*?#4\.27:.*?#4\.28:.*?#4\.29:.*?#4\.30:.*?#4\.31:.*?my \$MATFILER_ver = 4\.31;/s,
	"MATAFILER history retains shared input discovery through version 4.31");
like($mataf4,
	qr/return unless \$summary->\{failed\};.*?my \@failureColumns.*?Job_category/s,
	'the end-of-run Slurm failure report is an occurrence matrix shown only when failures exist');
like($mataf4,
	qr/my %runOptions = \(.*?operationMode.*?sharedTmpDir.*?nodeTmpDir.*?from.*?to.*?submit.*?loopCount.*?loopInitialCount.*?loopWindowSize/s,
	'runtime and command-line controls are consolidated in one named hash');
unlike($mataf4,
	qr/my \$(?:ARGV0|sharedTmpDirP|nodeTmpDirBase|FROM1|TO1|doSubmit|loop2completion|loop2c_winsize|loop2completion_ini)\b/,
	'retired standalone runtime option globals are not reintroduced');
like($mataf4,
	qr/my %checkpointNames = \(.*?preAssemblyDone.*?assemblyDone.*?my %sampleCheckpoints = \(.*?primaryMapping.*?supportMapping.*?mappingComplete.*?primaryConsensus.*?supportConsensus/s,
	'fixed and sample-resolved checkpoints are consolidated into named hashes');
like($mataf4,
	qr/sub metagAssemblyRun.*?spaceInAssGrp\(\$curSmpl, 1\).*?assembly_cores_for_input\(.*?local \$MFopt\{AssemblyCores\} = \$assemblyCores/s,
	'assembly submissions use group-wide primary and support input for automatic core scaling');
like($mataf4,
	qr/\$MFopt\{AssemblyCores\} = 0;.*?\$MFconfig\{schedulerPollSeconds\} = 20;.*?\$MFconfig\{schedulerCapacityCheckJobs\} = 10;.*?schedulerPollSeconds=i.*?schedulerCapacityCheckJobs=i/s,
	'automatic assembly scaling and batched scheduler capacity checks are enabled by default');
like($mataf4,
	qr/qsubSystemWaitMaxJobs\(\s*\$MFconfig\{checkMaxNumJobs\}, \$MFconfig\{killDepNever\}, \$QSBoptHR/s,
	'the main sample loop shares scheduler throttle state instead of querying once per sample');
my $finished_shortcut_position = index($mataf4, 'Sample already complete; no jobs submitted');
my $sample_scratch_creation_position = index(
	$mataf4, 'make_path($smplTmpDir) unless -d $smplTmpDir',
);
ok($finished_shortcut_position >= 0
		&& $sample_scratch_creation_position > $finished_shortcut_position,
	'completed samples bypass sample scratch creation');
my $sample_log_creation_position = index(
	$mataf4, 'make_path($logDir) unless -d $logDir',
);
ok($sample_log_creation_position > $finished_shortcut_position,
	'completed samples bypass creation of a missing per-sample LOGandSUB directory');
unlike($mataf4, qr/present: \$curOutDir/,
	'completed samples no longer print a repetitive output-directory message');
like($mataf4,
	qr/Submitting deferred assembly-group mapping jobs.*?my \$publicationDeps = normalise_job_dependencies\(.*?MultiContigStats\.sh/s,
	'assembly groups release mapping, producer publication, and contig statistics in stages');
like($mataf4,
	qr/params\{mappingCommand\}.*?my \$mappingCommand = delete\(\$mapparhr->\{mappingCommand\}\).*?qsubSystem\(\$combinedScript/s,
	'MAP and sort/depth are submitted as one combined scheduler job');
unlike($mataf4, qr/SBAM/,
	'the obsolete separate SBAM scheduler stage is no longer named or submitted');
unlike($mataf4, qr/sub (?:clean_tmp|manageFiles)\b|MapCopies|AssCopies|MapSupCopies|MapCopiesNoDel|PsAssCopies/,
	'result-moving cleanup routines and copy-state queues have been removed');
like($mataf4,
	qr/my \$publishStage = "\$finalD\/\.\$baseN\.mapping-stage".*?quickcheck.*?touch \$cramSTO/s,
	'mapping publication validates staged output and writes its completion marker last');
like($mataf4,
	qr/if \(\$outstat && \$outstat2 && !\$breakpointDone && \$mappingCommand eq ""\).*?my \$breakpointCoverage.*?return \(\$breakpointJob/s,
	'a missing breakpoint report is repaired from canonical coverage without rebuilding mapping');
like($mataf4,
	qr/sub check_map_done.*?Only canonical outputs count as complete.*?sub check_depth_done/s,
	'mapping completeness is determined only from canonical final outputs');
unlike($mataf4, qr/Migrating legacy scratch|legacyMetagDirs|moveMappings/,
	'outputs and partial files from older MATAFILER layouts are ignored');
like($mataf4,
	qr/my \$mappingArtifactsPresent\s*=.*?if \(!\$efinAssLoc && !\$ePreAssmbly && \$mappingArtifactsPresent\)/s,
	'a missing clean-run assembly does not masquerade as a mapping redo');
like($mataf4,
	qr/my \$variantCommonInputsPublished.*?my \$variantCommonInputsPending.*?my \$primaryVariantMappingPending.*?my \$supportVariantMappingPending.*?my \$variantInputsMayBePending.*?allowPendingInputs => \(\$variantInputsMayBePending/s,
	"variant calls accept each published input or its own same-pass producer");
like($mataf4,
	qr/my \$runConsensus =.*?callConsSNP.*?callConsSNPSupp.*?SNPconsensus_vcf\(\\%SNPinfo\) if \$runConsensus/s,
	'an SV-only request does not execute the SNP-consensus workflow');
like($mataf4,
	qr/my \$supportMappingPublished\s*=.*?\$eFinSupMapCovGZ.*?\$supportMappingPublished/s,
	'hybrid binning waits for support mapping to be published');
like($mataf4,
	qr/my \$binningArtifactsPresent\s*=.*?redoing binning due to support mapping not included/s,
	'support-map repair is reported only when binning artifacts actually exist');
like($mataf4,
	qr/my \$supportCoverageRequired\s*=\s*\$locMapSup2Assembly\s*&&\s*\$efinAssLoc.*?\$calcCoverage.*?\$supportCoverageRequired/s,
	'support coverage becomes required only after final hybrid assembly publication');
like($mataf4,
	qr/sub runContigStats.*?\$requireSupportCoverage.*?contig_stats_coverage_complete\(\$ContigStatsDir, "Cov\.sup"\)/s,
	'ContigStats completeness follows the active hybrid phase instead of map metadata alone');
like($mataf4,
	qr/indication that hybrid assembly is already done.*?return \(\$ePreAssmbly,\$doPreAssmFlag,0,0\)/s,
	'completed final assemblies are no longer gated by retained preassembly packages');
like($mataf4,
	qr/sub spadesAssembly.*?\$finalOut\.assembly-stage.*?mv \$stageOut \$finalOut.*?sub longRdAssembly.*?\$finalOut\.hybrid-stage.*?mv \$stageOut \$finalOut.*?sub megahitAssembly.*?\$finalOut\.assembly-stage.*?mv \$stageOut \$finalOut/s,
	'all assembly producers validate and rotate staged output into the canonical directory');
like($mataf4,
	qr/sub createPsAssLongReads.*?my \$psStage = "\$psFinal\.stage".*?mv -f \$psStage \$psFinal.*?touch \$pseudoAssFileFlag/s,
	'pseudoassembly is atomically published before its completion marker');
like($mataf4,
	qr/sub calcCoverage2nd.*?return \(\$jobName,\$tmpCmd\)/s,
	'completed secondary coverage does not generate a perpetual combined mapping job');
like($mataf4,
	qr/deferMappingCleanup => 1.*?my \$secondMapCmd = \$bigMap\."\\n"\.\$bigCov.*?sharedMapWork/s,
	'multi-reference mapping keeps shared node-local alignments until every reference is published');
like($mataf4, qr/if \(\$MFopt\{DoMetaBat2\} && !\$doPreAssmFlag.*?\$AssemblyGo/s,
	'binning can be scheduled on the first pass once the final group assembly is scheduled');
like($mataf4, qr/append_job_dependencies\(\\\$AsGrps\{\$cAssGrp\}\{BinDeps\}, \$contRun\)/,
	'all binners wait for the contig-stat job that creates their coverage inputs');
like($mataf4, qr/submitGenomeBinner\(\$binnerTmp,\$finAssLoc/,
	'first-pass binning consumes the producer-published final assembly path');
like($mataf4,
	qr/Submitting deferred assembly-group mapping jobs.*?MultiContigStats\.sh.*?Submitting deferred assembly-group Cons jobs.*?MultiConsensus\.sh/s,
	'assembly groups release mapping, ContigStats, and deferred Cons jobs in producer order');
like($mataf4,
	qr/immediateSubm => \(\$variantSubmissionDeferred \? 0 : 1\).*?PostConsCmd.*?variantSubmissionCommands/s,
	'first-pass Cons jobs are submitted immediately or retained until group producer ids are available');
like($mataf4,
	qr/sub deferLoopProducerWave.*?return 0 unless \$runOptions\{loopCount\}.*?MFnext\(.*?loop2C_check\(/s,
	'producer-wave deferral is loop-only and preserves lock and loop bookkeeping');
like($mataf4,
	qr/seedUnzip2tmp.*?append_job_dependencies\(\\\$AsGrps\{\$cAssGrp\}\{SeqClnDeps\}, \$jdep\).*?deferLoopProducerWave\(\s*'input staging'.*?sdmClean.*?append_job_dependencies\(\\\$AsGrps\{\$cAssGrp\}\{SeqClnDeps\}, \$sdmjN\).*?deferLoopProducerWave\(\s*'quality filtering'.*?removeHostSeqs.*?append_job_dependencies\(\\\$AsGrps\{\$cAssGrp\}\{SeqClnDeps\}, \$sdmjN\).*?deferLoopProducerWave\(\s*'host filtering'.*?mergeReads.*?deferLoopProducerWave\(\s*'read merging'/s,
	'input waves also hold back a multi-sample assembly until every member is ready');
like($mataf4,
	qr/SeqClnDeps\}, \$sdmjN.*?deferLoopProducerWave\(\s*'input preparation'.*?deferLoopProducerWave\(\s*'assembly-group input preparation'.*?metagAssemblyRun/s,
	'host cleaning and other input producers must finish before group assembly is submitted');
like($mataf4,
	qr/metagAssemblyRun.*?deferLoopProducerWave\(\s*'assembly'.*?genePredictions/s,
	'loop mode stops after assembly instead of creating an assembly-to-annotation chain');
like($mataf4,
	qr/my \$currentMappingDeps.*?append_job_dependencies\(\\\$currentMappingDeps, \$deferredDeps\).*?append_job_dependencies\(\\\$currentMappingDeps, \$map2Ctgs_2\).*?append_job_dependencies\(\\\$currentMappingDeps, \$mapSup2Ctgs_2\).*?my \$mappingWaveDeps = normalise_job_dependencies\(\s*\$currentMappingDeps, \$AsGrps\{\$cAssGrp\}\{MapDeps\}.*?deferLoopProducerWave\(\s*'assembly mapping', \$mappingWaveDeps.*?#---------------- producer barriers/s,
	'deferred, primary, and support mappings form an assembly-group wave before downstream work');
like($mataf4,
	qr/my \$currentContigStatsDeps.*?append_job_dependencies\(\\\$currentContigStatsDeps, \$jdep\).*?my \$contigStatsWaveDeps = normalise_job_dependencies\(\s*\$currentContigStatsDeps, \$AsGrps\{\$cAssGrp\}\{BinDeps\}.*?deferLoopProducerWave\(\s*'contig statistics', \$contigStatsWaveDeps.*?submitGenomeBinner.*?createConsSNPandSVs/s,
	'contig statistics complete in an earlier pass than binning and consensus');
like($mataf4,
	qr/for my \$SNPinfo \(\@pendingSecondMapSNP\).*?next unless -s \$secondReference && -s \$secondMapping;.*?createConsSNPandSVs\(\$SNPinfo\)/s,
	'secondary-reference ConsSNP also waits for its published reference and mapping');
like($mataf4,
	qr/my \$consensusFastasComplete.*?fileGZe\(\$contigsSNP\).*?fileGZe\(\$genePredSNP\).*?fileGZe\(\$genePredAASNP\).*?my \$primaryConsensusComplete.*?my \$supportConsensusComplete/s,
	"ConsSNP completion covers every requested consensus FASTA, VCF, and stone");
like($mataf4,
	qr/my \$genePredGff.*?my \$boolGenePredOK = fileGZe\(\$genePredProtein\) && fileGZe\(\$genePredGff\).*?sub genePredictions.*?fileGZs\("\$expectedD\/genes\$bacmark\.gff"\)/s,
	"gene prediction is complete only when its GFF is published with the protein output");
like($mataf4,
	qr/my \$consensusNeedsContigStats.*?!\$cleanupContigStatsComplete.*?\|\| \$consensusNeedsContigStats.*?my \$variantPrerequisiteDeps = normalise_job_dependencies\(\s*\$publicationDeps, \$currentContigStatsDeps, \$AsGrps\{\$cAssGrp\}\{BinDeps\}.*?jdeps => \$variantPrerequisiteDeps/s,
	"ConsSNP explicitly repairs and depends on the requested ContigStats products");
unlike($mataf4,
	qr/\{(?:AssemblJobName|MapDeps|BinDeps|SeqClnDeps|SeqUnZDeps|UnzpDeps|readDeps|DiamDeps|scndMapping|prodRun)\}\s*\.=/,
	'central workflow dependencies are not assembled with fragile string concatenation');

open my $subm_source, '<', File::Spec->catfile($Bin, '..', 'Mods', 'Subm.pm')
	or die "Cannot inspect Mods/Subm.pm: $!";
my $subm = do { local $/; <$subm_source> };
close $subm_source;
like($subm, qr/\$waitJID\s*=\s*normalise_job_dependencies\(\$waitJID\)/,
	'all scheduler submissions use the shared dependency normalizer');

unlike($mataf4, qr/CSfinJobName/,
	'per-sample ContigStats jobs are not serialized behind another sample ContigStats job');
unlike($subm, qr/length\(\$waitJID\)\s*>\s*3/,
	'short valid scheduler job ids are not silently discarded');
like($subm,
	qr/sub MFnext.*?recordSampleLockJobs\(\$lckFile, \$aR, \$QSBoptHR\)/s,
	'sample completion persists dependency IDs without submitting an RMLOCK job');
unlike($subm, qr/RMLCK|rmLock\.sh/,
	'the submission layer no longer creates scheduler jobs solely to release sample locks');
like($subm,
	qr/my \$pollSeconds = defined\(\$optHR->\{jobPollSeconds\}\).*?sleep\(\$pollSeconds\)/s,
	'loop waiting uses the configurable scheduler polling interval');

open my $group_source, '<', File::Spec->catfile($Bin, '..', 'Mods', 'GenoMetaAss.pm')
	or die "Cannot inspect Mods/GenoMetaAss.pm: $!";
my $groups = do { local $/; <$group_source> };
close $group_source;
like($groups, qr/sub resetAsGrps.*?\{UnzpDeps\}\s*=\s*""/s,
	'loop reset removes completed unzip job ids before the next submission pass');

my %loop_groups = (group => {
	UnzpDeps => 'old-job',
	FilterSeq1 => [], FilterSeq2 => [], FilterSeqS => [], ReadTec => [],
	CntPreAss => 2, CntPreAssMiss => 3,
	CntPreAssNoPrim => 4, preAsmblDir => ['old-package'],
	AssemblSmplDirs => "/old/sample\n", PostConsCmd => 'old-consensus',
});
resetAsGrps(\%loop_groups);
is($loop_groups{group}{UnzpDeps}, '',
	'loop reset behavior removes stale unzip dependencies');
is_deeply(
	[@{$loop_groups{group}}{qw(CntPreAss CntPreAssMiss CntPreAssNoPrim)}],
	[0, 0, 0],
	'hybrid readiness counters do not accumulate across loopTillComplete passes');
is_deeply($loop_groups{group}{preAsmblDir}, [],
	'hybrid package paths are rebuilt once per loop pass without duplicates');
is($loop_groups{group}{AssemblSmplDirs}, '',
	'assembly membership output is rebuilt once per loop pass without duplicates');
is_deeply(
	[@{$loop_groups{group}}{qw(PostConsCmd)}],
	[''],
	'deferred downstream submission state is cleared between loop passes');

like($mataf4,
	qr/sub longRdAssembly\s*\{.*?\$finalOut\s*=~\s*s\{\/\+\$\}\{\};.*?my \$stageOut\s*=\s*"\$finalOut\.hybrid-stage"/s,
	'hybrid publication strips trailing slashes before creating sibling staging directories');

like($mataf4,
	qr/\$mapAssFlag\s*=\s*1\s+if\s*\(\$map\{\$curSmpl\}\{hasPrimaryRds\}\s*&&\s*\$MFopt\{map2Assembly\}\s*&&\s*!\$eFinMapCovGZ/s,
	'primary mapping cannot trigger staging for support-only samples');
like($mataf4,
	qr/if\s*\(!\$presence\s*&&\s*\$runThis\)\{\s*print\s+"sdm'ing support reads\.\.\\n"\s+if\s*\(\$useXtras\)/s,
	'support SDM message is emitted only when an SDM job is submitted');

like($mataf4,
	qr/if \(\$is2ndMap && \$MFopt\{MapRewrite2nd\}\)/,
	'secondary-map rewrite cannot delete primary assembly mappings');
unlike($mataf4, qr/gunzip \$REF/,
	'compressed mapping references are never decompressed in place');

like($mataf4,
	qr/print "AssmblGrp: ".*?CntAimAss\} > 1 && !\$MFconfig\{silent\}/s,
	'assembly-group progress is shown only for multi-sample assembly groups');
like($mataf4,
	qr/print "MapGroup: ".*?CntAimMap\} > 1 && !\$MFconfig\{silent\}/s,
	'mapping-group progress is shown only for multi-sample mapping groups');
like($mataf4,
	qr/my \$hybridAssemblyRequested = \$MFopt\{DoAssembly\} == 5.*?print "precnt: .*?if \(\$hybridAssemblyRequested && !\$MFconfig\{silent\}\)/s,
	'preassembly count is reported only for requested hybrid assemblies');
unlike($mataf4, qr/Running Contig Stats on assembly \(\$immSubm\)/,
	'normal immediate ContigStats submission no longer emits a routine status line');
like($mataf4,
	qr/\$pigzBin -dc \$compressedReference > \$stagedReference.*?mappingReference/s,
	'compressed mapping references are staged in mapper-local storage');
like($mataf4,
	qr/canonical BAM\/CRAM.*?if \(\$outstat && !\$outstat2 && \$mappingCommand eq ""\).*?\$smtBin depth/s,
	'missing coverage is repaired directly from the canonical alignment');
like($mataf4,
	qr/\$mappingInputSizeMB = \$supportRds\s*\? \(\$map\{\$curSmpl\}\{inputXFileSizeMB\}/s,
	'support mapping resources are based on support-read input size');
unlike($mataf4, qr/if \(1\)\{\s*#always active #\(\$numLib > 1\)/,
	'single BAM segments no longer pass through an unconditional samtools cat');
like($mataf4,
	qr/if \(\@bamParts > 1\).*?\$smtBin cat.*?elsif \(\@bamParts == 1\).*?mv \$bamParts\[0\]/s,
	'samtools cat is reserved for multiple BAM segments');
like($mataf4,
	qr/my \$serialiseDuplicateSorts = \$locDoRmDup.*?my \$sortProcessCount = \(\$locDoRmDup && !\$serialiseDuplicateSorts\) \? 2 : 1.*?my \$sortMemoryMB = .*?\(\$numCore \* \$sortProcessCount\)/s,
	'samtools per-thread memory shares the total budget across concurrent sorts');
like($mataf4,
	qr/if \(\$serialiseDuplicateSorts\).*?sort -n .*?-o \$nameSortedBam.*?fixmate .*?\$nameSortedBam - \| .*?sort/s,
	'large mappings materialize name-sorted data instead of keeping two sorts resident');
like($mataf4,
	qr/my \$sortMemoryCapMB = \$serialiseDuplicateSorts \? 768 : 2048/s,
	'large mappings cap each sort thread at a conservative memory arena');
like($mataf4,
	qr/my \$controlledSortMB = \$sortMemoryMB \* \$numCore \* \$sortProcessCount.*?\$controlledSortMB \* 1\.25.*?4096/s,
	'the mapping scheduler request includes sort-process and non-arena headroom');
like($mataf4,
	qr/my %requiredMappers.*?next if \(\$mapper == 3 \|\| \$mapper == 5\).*?next if \(\$cmdDB eq ""\)/s,
	'assembly index jobs are deduplicated and omitted for direct-FASTA mappers');
like($mataf4,
	qr/my \$recordReadTechnology = libraryTechnology\(\$libraries,.*?my \$readTec = \$recordReadTechnology \|\| \$declaredReadTechnology/s,
	'automatic mapper selection uses the validated technology attached to the library records');

done_testing;
