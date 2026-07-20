use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

my $root = tempdir(CLEANUP => 1);
my $output = File::Spec->catdir($root, 'output');
my $scratch = File::Spec->catdir($root, 'scratch');
my $group_dir = File::Spec->catdir($output, 'AssmblGrp_group', 'metag');
my $state_dir = File::Spec->catdir($group_dir, '.cleanup-indexes');
my $assembly = File::Spec->catfile($group_dir, 'scaffolds.fasta.filt');
my $cleaner = File::Spec->catfile($Bin, '..', 'secScripts', 'cleanup_finished_sample.pl');
make_path($group_dir, $scratch);

sub write_file {
	my ($path, $content) = @_;
	my (undef, $dir) = File::Spec->splitpath($path);
	make_path($dir) unless -d $dir;
	open my $fh, '>', $path or die $!;
	print {$fh} defined($content) ? $content : "x\n";
	close $fh;
}

sub run_cleaner {
	my (%args) = @_;
	my @command = (
		$^X, $cleaner,
		'--sample', $args{sample},
		'--state-dir', $state_dir,
		'--allowed-root', $output,
	);
	push @command, map { ('--member', $_) } @{$args{members}};
	push @command, map { ('--member-lock', $_) } @{$args{member_locks} || []};
	push @command,
		'--mapping-dir', $args{mapping_dir},
		'--sample-temp', $args{sample_temp},
		'--scratch-root', $scratch,
		'--assembly', ($args{assembly} // $assembly),
		'--snp-log-dir', $args{snp_log_dir};
	return system(@command);
}

write_file($assembly, ">ctg\nACGT\n");
for my $index ("$assembly.fai", "$assembly.mmi", "$assembly.bw2.1.bt2", "$assembly.kma.seq.b", "$assembly.pac") {
	write_file($index);
}

my %sample_paths;
for my $sample (qw(S1 S2)) {
	my $sample_dir = File::Spec->catdir($output, $sample);
	my $mapping_dir = File::Spec->catdir($sample_dir, 'mapping');
	my $snp_log_dir = File::Spec->catdir($sample_dir, 'LOGandSUB', 'SNP');
	my $sample_temp = File::Spec->catdir($scratch, $sample);
	make_path($mapping_dir, $snp_log_dir, $sample_temp);
	write_file(File::Spec->catfile($mapping_dir, "$sample-smd.bam"));
	write_file(File::Spec->catfile($mapping_dir, "$sample-smd.bam.bai"));
	write_file(File::Spec->catfile($mapping_dir, "$sample.sup-smd.cram.crai"));
	write_file(File::Spec->catfile($snp_log_dir, "$sample.0.bed"));
	write_file(File::Spec->catfile($sample_temp, 'temporary.fq'));
	$sample_paths{$sample} = {
		mapping_dir => $mapping_dir,
		snp_log_dir => $snp_log_dir,
		sample_temp => $sample_temp,
	};
}

is(run_cleaner(
	sample => 'S1', members => [qw(S1 S2)],
	%{$sample_paths{S1}},
), 0, 'first completed sample cleanup succeeds');
ok(!-e File::Spec->catfile($sample_paths{S1}{mapping_dir}, 'S1-smd.bam.bai'),
	'sample BAM index is removed');
ok(!-e File::Spec->catfile($sample_paths{S1}{mapping_dir}, 'S1.sup-smd.cram.crai'),
	'sample supplementary CRAM index is removed');
ok(-s File::Spec->catfile($sample_paths{S1}{mapping_dir}, 'S1-smd.bam'),
	'canonical alignment is retained');
ok(!-d $sample_paths{S1}{sample_temp}, 'sample temporary directory is removed');
ok(-s "$assembly.mmi", 'shared mapper index remains until every group member completes');

is(run_cleaner(
	sample => 'S2', members => [qw(S1 S2)],
	%{$sample_paths{S2}},
), 0, 'last completed group member cleanup succeeds');
ok(-s $assembly, 'canonical assembly is retained');
ok(!-e "$assembly.fai" && !-e "$assembly.mmi" && !-e "$assembly.bw2.1.bt2"
	&& !-e "$assembly.kma.seq.b" && !-e "$assembly.pac",
	'group-owned FASTA and mapper indexes are removed after the group barrier');

write_file("$assembly.mmi");
my $active_lock = File::Spec->catfile($output, 'S1', 'LOGandSUB', 'MGTK.locked');
write_file($active_lock);
is(run_cleaner(
	sample => 'S2', members => [qw(S1 S2)], member_locks => [$active_lock],
	%{$sample_paths{S2}},
), 0, 'a completed member may recheck the barrier while another job is active');
ok(-s "$assembly.mmi", 'active member lock prevents premature shared-index cleanup');
unlink $active_lock or die $!;
is(run_cleaner(
	sample => 'S1', members => [qw(S1 S2)], member_locks => [$active_lock],
	%{$sample_paths{S1}},
), 0, 'group cleanup can run after the active job lock is released');
ok(!-e "$assembly.mmi", 'shared index is removed only after the active job finishes');

write_file($assembly, ">ctg\nACGTACGT\n");
write_file("$assembly.mmi");
is(run_cleaner(
	sample => 'S1', members => [qw(S1 S2)],
	%{$sample_paths{S1}},
), 0, 'first sample records completion against a rebuilt assembly');
ok(-s "$assembly.mmi", 'markers from an older assembly cannot satisfy the group barrier');
is(run_cleaner(
	sample => 'S2', members => [qw(S1 S2)],
	%{$sample_paths{S2}},
), 0, 'last sample records completion against the rebuilt assembly');
ok(!-e "$assembly.mmi", 'rebuilt assembly index is removed after fresh group completion');

my $external_root = File::Spec->catdir($root, 'external');
my $external_assembly = File::Spec->catfile($external_root, 'reference.fa');
make_path($external_root);
write_file($external_assembly);
write_file("$external_assembly.mmi");
my $external_state = $state_dir;
is(run_cleaner(
	sample => 'external', members => ['external'],
	assembly => $external_assembly,
	mapping_dir => $sample_paths{S1}{mapping_dir},
	snp_log_dir => $sample_paths{S1}{snp_log_dir},
	sample_temp => $sample_paths{S1}{sample_temp},
), 0, 'cleanup safely accepts an external assembly');
ok(-s "$external_assembly.mmi", 'indexes adjacent to an external reference are retained');

my $scratch_alias = File::Spec->catdir($root, 'scratch-alias');
symlink($scratch, $scratch_alias) or die "Cannot create scratch alias: $!";
my $alias_temp_real = File::Spec->catdir($scratch, 'alias-sample');
my $alias_temp_path = File::Spec->catdir($scratch_alias, 'alias-sample');
make_path($alias_temp_real);
write_file(File::Spec->catfile($alias_temp_real, 'temporary.fq'));
is(run_cleaner(
	sample => 'alias-sample', members => ['alias-sample'],
	mapping_dir => $sample_paths{S1}{mapping_dir},
	snp_log_dir => $sample_paths{S1}{snp_log_dir},
	sample_temp => $alias_temp_path,
), 0, 'cleanup accepts a sample temporary path through a scratch symlink alias');
ok(!-d $alias_temp_real, 'scratch alias cleanup removes the canonical sample temporary directory');

open my $mata_fh, '<', File::Spec->catfile($Bin, '..', 'MATAF4.pl') or die $!;
my $mata_source = do { local $/; <$mata_fh> };
close $mata_fh;
unlike($mata_source, qr/'--mode', 'invalidate'/,
	'cleanup script is not called before sample completion');
like($mata_source, qr/runFinishedCleanup\(finishedCleanupArguments\(/s,
	'fully completed samples invoke centralized cleanup');
like($mata_source,
	qr/MFnext\(\$smplLockF,\\\@sampleDeps.*?submitFinishedCleanup\(.*?\\\@sampleDeps.*?loop2C_check/s,
	'normal submissions enqueue cleanup after work completion and lock release');
like($mata_source,
	qr/sub submitFinishedCleanup.*?afterAny\} = 0;.*?qsubSystem\(.*?\$dependencyString/s,
	'cleanup uses successful scheduler dependencies rather than after-any execution');
like($mata_source, qr/sub submitFinishedCleanup.*?--kill-on-invalid-dep=yes/s,
	'failed Slurm dependency chains cancel cleanup instead of leaving it pending forever');
unlike($mata_source, qr/system "rm -f \$finalCommAssDir\/scaffolds\.fasta\.filt/s,
	'legacy inline mapper-index deletion is removed');

done_testing;
