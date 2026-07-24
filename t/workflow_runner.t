use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::WorkflowRunner qw(run_workflow_preflight);

sub write_file {
	my ($path, $contents) = @_;
	(my $dir = $path) =~ s{/[^/]+$}{};
	make_path($dir);
	open my $fh, '>', $path or die "Cannot write $path: $!";
	print {$fh} $contents;
	close $fh;
}

sub actions_by_id {
	my ($plan) = @_;
	return map { $_->{id} => $_ } @{$plan->{actions}};
}

my $root = tempdir(CLEANUP => 1);
$root =~ s{\\}{/}g;
my $sample_a = "$root/run/A/";
my $sample_b = "$root/run/B/";
my $run_tmp = "$root/tmp/hybrid-run";
my $group_dir = "$root/run/AssmblGrp_hybrid/metag";
my %map = (
	opt => { smpl_order => ['A', 'B'] },
	A => {
		SmplID => 'A', wrdir => $sample_a, AssGroup => 'hybrid',
		SupportReads => 'PB:/reads/a.fastq.gz', hasPrimaryRds => 1,
	},
	B => {
		SmplID => 'B', wrdir => $sample_b, AssGroup => 'hybrid',
		SupportReads => '', hasPrimaryRds => 1,
	},
);
my %groups = (hybrid => { CntAimAss => 2, SupportReads => ',PB:/reads/a.fastq.gz' });
my %options = (
	assembly_mode => 5,
	map_to_assembly => 1,
	map_support_to_assembly => 1,
	run_tmp_dir => $run_tmp,
);

# A represents an interrupted hybrid packaging job. B is already complete and
# must be preserved across subsequent loop iterations.
my $package_a = "$run_tmp/A/preAssmblGrp_hybrid";
my $package_b = "$run_tmp/B/preAssmblGrp_hybrid";
write_file("$package_a/scaffolds.fasta.filt", ">partial\nACGT\n");
write_file("$package_a/moved.sto", "done\n");
write_file("$package_b/scaffolds.fasta.filt", ">complete\nACGT\n");
write_file("$package_b/Coverage.percontig.gz", "coverage\n");
write_file("$package_b/Coverage.median.percontig.gz", "median\n");
write_file("$package_b/mapping.coverage.gz", "mapping coverage\n");
write_file("$package_b/breakpoints.tsv.gz", "breakpoints\n");
write_file("$package_b/package.manifest.tsv", "key\tvalue\nschema_version\t2\n");
write_file("$package_b/moved.sto", "done\n");
write_file("$sample_a/mapping/A-smd.cram.sto", "done\n");
write_file("$sample_a/mapping/A-smd.bam.coverage.gz", "partial coverage\n");

my $dry_run = run_workflow_preflight(
	map => \%map, groups => \%groups, options => \%options,
	apply_repairs => 0, allow_group_rewrite => 0, iteration => 0,
);
ok(-e "$package_a/moved.sto", '-submit 0 style preflight does not repair files');
cmp_ok($dry_run->{repairs}{would_remove_targets}, '>=', 4,
	'dry-run reports safe partial targets it would remove');

my $first = run_workflow_preflight(
	map => \%map, groups => \%groups, options => \%options,
	apply_repairs => 1, allow_group_rewrite => 0, iteration => 0,
);
ok(!-e "$package_a/scaffolds.fasta.filt" && !-e "$package_a/moved.sto",
	'interrupted hybrid package is invalidated automatically');
ok(-e "$package_b/scaffolds.fasta.filt" && -e "$package_b/moved.sto",
	'completed package for another assembly-group member is retained');
ok(!-e "$sample_a/mapping/A-smd.cram.sto",
	'interrupted downstream mapping marker is invalidated automatically');
my %first_actions = actions_by_id($first->{plan});
ok($first_actions{'sample:A:submit:preassembly'},
	'invalidated package is planned for resubmission');
ok(!$first_actions{'sample:B:submit:preassembly'},
	'completed package is not planned for resubmission');
ok(grep($_ eq 'sample:A:submit:preassembly',
	@{$first_actions{'group:hybrid:submit:assembly'}{depends_on}}),
	'combined hybrid group assembly waits for the remaining package');

# Simulate the scheduler completing this pass before loopTillComplete reaches
# its next preflight boundary.
write_file("$package_a/scaffolds.fasta.filt", ">complete\nACGT\n");
write_file("$package_a/Coverage.percontig.gz", "coverage\n");
write_file("$package_a/Coverage.median.percontig.gz", "median\n");
write_file("$package_a/mapping.coverage.gz", "mapping coverage\n");
write_file("$package_a/breakpoints.tsv.gz", "breakpoints\n");
write_file("$package_a/package.manifest.tsv", "key\tvalue\nschema_version\t2\n");
write_file("$package_a/moved.sto", "done\n");
write_file("$group_dir/scaffolds.fasta.filt", ">group\nACGT\n");
write_file("$group_dir/ass.done.sto", "done\n");
write_file("$group_dir/smpls_used.txt", "$sample_a\n$sample_b\n");
write_file("$group_dir/genePred/proteins.shrtHD.faa", ">gene\nM\n");

my $second = run_workflow_preflight(
	map => \%map, groups => \%groups, options => \%options,
	apply_repairs => 1, allow_group_rewrite => 0, iteration => 1,
);
my %second_actions = actions_by_id($second->{plan});
is($second->{iteration}, 1, 'loop iteration is recorded in the audit result');
ok(!$second_actions{'sample:A:submit:preassembly'}
	&& !$second_actions{'sample:B:submit:preassembly'},
	'next loop pass recognizes every completed preassembly package');
ok(!$second_actions{'group:hybrid:submit:assembly'},
	'next loop pass recognizes the completed final group assembly');
ok($second_actions{'sample:A:submit:mapping'}
	&& $second_actions{'sample:A:submit:support_mapping'},
	'next loop pass continues with primary and support mapping');
is($second_actions{'sample:A:submit:mapping'}{scope}{id}, 'A',
	'downstream work remains scoped to the correct group member');

# A changed group remains protected unless the established explicit gate is
# supplied, even though ordinary partial-file repairs are automatic.
write_file("$group_dir/smpls_used.txt", "$sample_a\n$root/run/DIFFERENT/\n");
my $protected = run_workflow_preflight(
	map => \%map, groups => \%groups, options => \%options,
	apply_repairs => 1, allow_group_rewrite => 0, iteration => 2,
);
ok(-e "$group_dir/scaffolds.fasta.filt",
	'assembly-group mismatch is not deleted without authorization');
cmp_ok($protected->{repairs}{blocked_repairs}, '>=', 1,
	'protected assembly-group repair is reported as blocked');

my $authorised = run_workflow_preflight(
	map => \%map, groups => \%groups, options => \%options,
	apply_repairs => 1, allow_group_rewrite => 1, iteration => 2,
);
ok(!-d $group_dir, 'authorized group-wide invalidation removes the stale group assembly');
cmp_ok($authorised->{repairs}{removed_targets}, '>=', 1,
	'authorized deletion is captured in the repair audit');

open my $source, '<', File::Spec->catfile($Bin, '..', 'MATAF4.pl')
	or die "Cannot inspect MATAF4.pl: $!";
my $mataf4 = do { local $/; <$source> };
close $source;
like($mataf4,
	qr/qsubSystemJobAlive\s*\([^;]+;\s*.*?runAutomaticWorkflowPreflight\(\$workflowIteration\).*?resetAsGrps/s,
	'loopTillComplete runs the repeated preflight after waiting and before resetting group state');
like($mataf4,
	qr/No jobs were submitted in the current iteration.*?elsif \(\$loop2completion\).*?\$JNUM = \$from - 1/s,
	'loopTillComplete only resets the sample index when another iteration is needed');
like($mataf4,
	qr/\$lastWindowPass = \$loop2completion == 1.*?overlap_loop_window\(.*?submitted_jobs => \$submittedThisIteration.*?last_pass => \$lastWindowPass.*?if \(\$overlapWindow->\{extended\}\).*?\$loop2completion = \$loop2completion_ini.*?return;/s,
	'light or final loop passes extend through one more sample block before waiting');
like($mataf4,
	qr/if \(\$overlapWindow->\{extended\}\).*?return;\s*}\s*\$loopIterationSubmissionStart =/s,
	'an extended pass retains its submission snapshot until both blocks reach the wait boundary');
like($mataf4,
	qr/\$doSubmit && !\$MFconfig\{rmSmplLocks\} && \$submittedThisIteration == 0.*?numActiveUserJobs\(\s*\$QSBoptHR,\s*1,\s*\[keys %loopSubmittedJobIds\].*?should_rerun_locked_window\(.*?active_job_threshold => \$MFconfig\{loopTillCompleteActiveJobs\}.*?if \(\$rerunLockedWindow.*?\$JNUM = \$from - 1;.*?return;/s,
	'a no-op retained-lock pass reruns its current window when few jobs remain active');
like($mataf4,
	qr/normalise_job_dependencies\(\\\@grandDeps\).*?\$loopJobId =~ s\/\^.*?\$loopSubmittedJobIds\{\$loopJobId\} = 1.*?qsubSystemJobAlive/s,
	'loopTillComplete retains submitted job ids for scoped active-job checks');
like($mataf4,
	qr/qsubSystemJobAlive\(\s*\\\@grandDeps,\s*\$QSBoptHR,\s*1,\s*\$MFconfig\{loopTillCompleteActiveJobs\}/s,
	'loopTillComplete starts its next pass at the configured active-job threshold');
like($mataf4,
	qr/\$loop2completion > 1 \|\| !\$loopFinalLockRetryUsed.*?\$loopFinalLockRetryUsed = 1/s,
	'the retained-lock policy permits only one extra retry beyond the configured pass budget');
like($mataf4,
	qr/getCmdLineOptions;.*?rewrite options cannot be combined with -loopTillComplete.*?setupHPC\(\)/s,
	'unsafe rewrite combinations fail before scheduler setup or job submission');
like($mataf4, qr/#4\.11:.*?loopTillComplete.*?#4\.12:.*?#4\.13:.*?my \$MATFILER_ver = 4\.13/s,
	'MATAFILER history records active-job loop triggering in version 4.13');

done_testing;
