use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::GenoMetaAss qw(resetAsGrps);
use Mods::WorkflowControl qw(
	advance_loop_window assembly_group_output_dirs hybrid_group_ready
	hybrid_package_complete missing_input_files source_input_files parse_ignored_samples
	sample_base_output_dir sample_is_ignored workflow_members_match
	normalise_job_dependencies append_job_dependencies augment_deferred_submission
);

is(normalise_job_dependencies('run12;;run7', ['run7', '', 'run3;run12']),
	'run12;run7;run3', 'job dependencies are flattened, deduplicated, and stable');
my $job_dependencies = 'run12;';
is(append_job_dependencies(\$job_dependencies, 'run7', 'run12'), 'run12;run7',
	'job dependencies are appended without empty or duplicate scheduler entries');

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

my $root = tempdir(CLEANUP => 1);
$root =~ s{\\}{/}g;
my $package = "$root/preAssmblGrp_gut";
write_file("$package/scaffolds.fasta.filt", ">contig\nACGT\n");
write_file("$package/Coverage.percontig.gz", "coverage\n");
write_file("$package/Coverage.median.percontig", "median\n");
write_file("$package/mapping.coverage.gz", "window coverage\n");
write_file("$package/breakpoints.tsv.gz", "compressed breakpoint fixture\n");
write_file("$package/package.manifest.tsv", "key\tvalue\nschema_version\t2\n");
write_file("$package/moved.sto", "done\n");
ok(hybrid_package_complete($package),
	'versioned hybrid package requires all package-local outputs');
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

open my $source, '<', File::Spec->catfile($Bin, '..', 'MATAF4.pl')
	or die "Cannot inspect MATAF4.pl: $!";
my $mataf4 = do { local $/; <$source> };
close $source;
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
like($mataf4, qr/HybridAssemblyComparison\.tsv/,
	'hybrid finalization requires a comparative assembly report');
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
	qr/my \$command_dependencies\s*=\s*normalise_job_dependencies\(.*?\$submitted\[-1\]/s,
	'deferred combined mapping commands retain their submission order');
like($mataf4,
	qr/Submitting deferred assembly-group mapping jobs.*?my \$publicationDeps = normalise_job_dependencies\(.*?MultiContigStats\.sh.*?MultiConsensus\.sh/s,
	'assembly groups release mapping, producer publication, contig statistics, and consensus in stages');
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
	qr/my \$supportMappingPublished\s*=.*?\$eFinSupMapCovGZ.*?\$supportMappingPublished/s,
	'hybrid binning waits for support mapping to be published');
like($mataf4,
	qr/my \$binningArtifactsPresent\s*=.*?redoing binning due to support mapping not included/s,
	'support-map repair is reported only when binning artifacts actually exist');
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
like($mataf4, qr/jdeps => \$AsGrps\{\$cAssGrp\}\{BinDeps\}.*?deferRegionPlanning/s,
	'first-pass consensus receives real scheduler dependencies and defers input inspection');
unlike($mataf4,
	qr/\{(?:AssemblJobName|MapDeps|BinDeps|SeqClnDeps|SeqUnZDeps|UnzpDeps|readDeps|DiamDeps|scndMapping|prodRun)\}\s*\.=/,
	'central workflow dependencies are not assembled with fragile string concatenation');

open my $subm_source, '<', File::Spec->catfile($Bin, '..', 'Mods', 'Subm.pm')
	or die "Cannot inspect Mods/Subm.pm: $!";
my $subm = do { local $/; <$subm_source> };
close $subm_source;
like($subm, qr/\$waitJID\s*=\s*normalise_job_dependencies\(\$waitJID\)/,
	'all scheduler submissions use the shared dependency normalizer');
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

done_testing;
