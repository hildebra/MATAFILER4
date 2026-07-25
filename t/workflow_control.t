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
	advance_loop_window overlap_loop_window parse_loop_spec should_rerun_locked_window assembly_group_output_dirs balanced_parallel_batches hybrid_group_ready
	hybrid_package_complete hybrid_package_sample_id hybrid_local_scratch_gb missing_input_files source_input_files parse_ignored_samples
	sample_base_output_dir sample_is_ignored workflow_members_match
	normalise_job_dependencies append_job_dependencies deferred_command_dependencies augment_deferred_submission
	commands_are_lightweight_filesystem cleanup_stage_barrier
);

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
	qr/sub postSubmQsub.*?my \$lock_file = \$QSBoptHR->\{LOCKfile\}.*?open my \$lock_fh/s,
	'accepted deferred submissions create the active sample lock');
like($mataf4,
	qr/my \$lightweightLocal = commands_are_lightweight_filesystem\(\$unzipcmd\).*?systemW \$unzipcmd.*?qsubSystem\(\$logDir\."UNZP\.sh"/s,
	'lightweight UZ setup runs locally while data-processing commands remain scheduled');
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
	qr/my \$assemblyDownstreamScheduled.*?my \$assemblyDownstreamDeferred.*?my \$variantCommonInputsPublished.*?my \$primaryVariantInputsPublished.*?my \$supportVariantInputsPublished.*?my \$variantInputsMayBePending.*?\$variantCommonInputsReady.*?\$primaryVariantInputsReady.*?\$doMapping.*?\$supportVariantInputsReady.*?\$mapSuppAssFlag.*?allowPendingInputs => \(\$variantInputsMayBePending/s,
	'variant calls accept published inputs or same-pass producer jobs with scheduler dependencies');
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
	qr/for my \$SNPinfo \(\@pendingSecondMapSNP\).*?next unless -s \$secondReference && -s \$secondMapping;.*?createConsSNPandSVs\(\$SNPinfo\)/s,
	'secondary-reference ConsSNP also waits for its published reference and mapping');
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
like($subm, qr/push\(\@\{\$aR\},\s*\$jN\)/,
	'lock-release submission tracks the returned job id rather than its shell command');

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
