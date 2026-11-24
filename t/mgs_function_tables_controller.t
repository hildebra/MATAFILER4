use strict;
use warnings;

# MGS.pl's function table stage (_submit_function_tables): which schemes it
# submits, how it waits for the geneCat annotation checkpoints, and when its
# own checkpoint lets a rerun skip the job. The controller cannot run without a
# cluster, so its subroutines run here with the scheduler replaced; the
# checkpoint is written by the job's real checkpoint command.

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;
use lib "$Bin/..", "$Bin/lib";
use MFTestConfig;

my $root = File::Spec->rel2abs("$Bin/..");
my $tmp = tempdir(CLEANUP => 1);

sub write_file {
	my ($path, $text) = @_;
	make_path((File::Spec->splitpath($path))[1]);
	open my $fh, '>', $path or die "$path: $!";
	print {$fh} $text;
	close $fh or die "$path: $!";
}
sub slurp {
	open my $fh, '<', $_[0] or die "$_[0]: $!";
	local $/;
	return <$fh>;
}

{
	package MGSFuncTables;
	use subs qw(sleep time); # the wait loop runs on a fake clock
	use Mods::Checkpoint qw(read_checkpoint checkpoint_valid write_checkpoint);
	use Mods::WorkflowResilience qw(retry_unlink);
	use Mods::FuncTables qw(
		function_table_scheme function_table_stage_stones plan_function_tables
		parse_function_table_list function_table_scheme_names
	);
	our ($clock, @submissions, @log);
	$clock = 1000;
	sub time() { return $clock; }          # prototype: "time - $start" must parse as (time) - $start
	sub sleep(;$) { $clock += $_[0]; return $_[0]; }
	sub printL { push @log, @_; }
	sub _mgs_workflow_stage {}
	sub getProgPaths { return "perl $root/secScripts/MGS/funcTablesMGS.pl" if $_[0] eq 'funcTablesMGS_scr'; die "unexpected tool $_[0]\n"; }
	sub qsubSystem { push @submissions, [@_]; return ('4711', $_[1]); }
}

package MGSFuncTables;
# the controller globals the subroutines below use
our ($GCd, $outD, $logDir, $checkpointWriter, $funcTablesSto, $funcTablesDir, $doSubmit,
	$wait4stone, $wait4funcStone, $wait4stoneTimeout, $MGSpipelineVersion, $binD, $binDctg, $binDctgFam,
	$doBinCtgsPerFam, $finalClustersFilt, $clock);
our (@funcTableSchemes, @funcTableModuleSets, @submissions, @log);
our (%QSBopt, %checkpointParameters, %stageCheckpointContract, %compatibleCheckpointReleases,
	%stageCheckpointCompatibilityLogged, %funcTablesCheckpointParameters);

my $source = main::slurp("$root/secScripts/MGS.pl");
my %stageContract = (eval { my %s = $source =~ /my %stageCheckpointContract = \((.*?)\);/s ? eval "($1)" : (); %s });
Test::More::is($stageContract{'mgs-function-tables'}, 1, 'the function table stage has its own checkpoint contract');
my @subs;
for my $name (qw(_checkpoint_parameters_for_stage _checkpoint_valid_for_stage _shell_quote
		_checkpoint_command _wait_for_stone _submit_function_tables)) {
	my ($sub) = $source =~ /(^sub \Q$name\E \{.*?^\})/ms;
	Test::More::ok(defined $sub, "MGS.pl defines $name") or Test::More::BAIL_OUT("cannot isolate $name");
	push @subs, $sub;
}
eval 'package MGSFuncTables; use strict; use warnings; ' . join("\n\n", @subs) . "; 1;"
	or Test::More::BAIL_OUT("MGS.pl function table subroutines do not compile: $@");

%stageCheckpointContract = %stageContract;
%compatibleCheckpointReleases = map { $_ => 1 } qw(0.54 0.55);
$MGSpipelineVersion = 0.56;
%checkpointParameters = (catalog_identity => 'catalog-A', cluster_id => 95, binner => 2,
	pipeline_version => $MGSpipelineVersion);
$GCd = "$tmp/GC";
$outD = "$GCd/Bin_SB";
$logDir = "$outD/LOGandSUB";
$funcTablesSto = "$logDir/checkpoints/funcTables.stone";
$funcTablesDir = "$outD/Annotation/FunctionTables/";
$checkpointWriter = "$root/helpers/writeCheckpoint.pl";
$wait4stoneTimeout = 30;
@funcTableSchemes = qw(KEGG CZy emNOG);
@funcTableModuleSets = qw(BSB GBM);
File::Path::make_path("$logDir/checkpoints");
main::write_file("$outD/MAGvsGC.txt.gz", "report");
main::write_file("$GCd/Anno/Func/DIAass_KGM.srt.gzgeneAss.gz", "1\tK00001\n");
main::write_file("$GCd/Anno/Func/emapper/eggNOGmapper_NOG.geneAss.gz", "1\tCOG0001\n");
my $stones = function_table_stage_stones($GCd);

sub reset_run { @submissions = (); @log = (); %funcTablesCheckpointParameters = (); %stageCheckpointCompatibilityLogged = (); }
sub logged { return join('', @log); }

# Dry run before any annotation stage finished: nothing to do, no waiting.
reset_run(); $doSubmit = 0; ($wait4stone, $wait4funcStone) = ('', '');
_submit_function_tables();
Test::More::is(scalar(@submissions), 0, 'no job before the annotation stages finished');
Test::More::like(logged(), qr/no annotation scheme is ready/, 'the controller says why');
Test::More::is($clock, 1000, 'a dry run does not wait');

# geneCat passes both stones; FuncAssign is done, eggNOG-mapper never finishes.
main::write_file($stones->{func}, "{\"stone\":1}\n");
reset_run(); $doSubmit = 1; ($wait4stone, $wait4funcStone) = ($stones->{emap}, $stones->{func});
_submit_function_tables();
Test::More::is($clock, 1030, 'a missing eggNOG-mapper stone is waited for until -wait4stoneTimeout');
Test::More::like(logged(), qr/stage 'emap' did not finish within 30s.*schemes are skipped/s,
	'its schemes are skipped instead of stopping the controller');
Test::More::like(logged(), qr/CZy skipped, no annotation file/, 'a database that was not annotated is skipped');
Test::More::is(scalar(@submissions), 1, 'one job is submitted');
my ($script, $cmd, $cores, $memory, $jobName, $dependency) = @{$submissions[0]};
Test::More::is($script, "$logDir/funcTablesMGS.sh", 'job script in the MGS log directory');
Test::More::is($dependency, '', 'the job has no scheduler dependency: its inputs are complete');
Test::More::like($cmd, qr{funcTablesMGS\.pl '-GCd' '\Q$GCd\E' '-binD' '\Q$outD\E' .*'-schemes' 'KEGG' '-modules' 'BSB,GBM'},
	'only the ready scheme is passed to the driver');
Test::More::like($cmd, qr{\Arm -f '[^']*function_tables\.manifest\.tsv'\n}, 'a stale manifest cannot pass the output checks');
Test::More::like($cmd, qr{test -s '[^']*FunctionTables//?KEGG\.tsv'}, 'the scheme table is checked before the checkpoint');

# The job succeeds: run its checkpoint command on the outputs.
main::write_file("$funcTablesDir/function_tables.manifest.tsv", "scheme\n");
main::write_file("$funcTablesDir/KEGG.tsv", "\tMGS.1\n");
my ($checkpointCmd) = $cmd =~ /^('perl' .*writeCheckpoint\.pl.*)$/m;
Test::More::ok(defined $checkpointCmd, 'the job writes the stage checkpoint');
Test::More::is(system('bash', '-c', $checkpointCmd), 0, 'the checkpoint command runs');
my $manifest = read_checkpoint($funcTablesSto);
Test::More::is($manifest->{parameters}{schemes}, 'KEGG', 'the checkpoint records the schemes');
Test::More::ok(!exists $manifest->{parameters}{pipeline_version}, 'and not the MGS release, so version bumps keep it');

reset_run(); $wait4stoneTimeout = 0; $wait4stone = ''; # eggNOG-mapper disabled: stop waiting for it
@funcTableSchemes = qw(KEGG);
_submit_function_tables();
Test::More::is(scalar(@submissions), 0, 'an unchanged rerun submits nothing');
Test::More::like(logged(), qr/\(KEGG\) are up to date/, 'and reports the tables as current');

# Re-annotation rewrites the FuncAssign checkpoint: the tables are made again.
$clock += 100;
main::write_file($stones->{func}, "{\"stone\":2, \"rerun\":true}\n");
reset_run();
_submit_function_tables();
Test::More::is(scalar(@submissions), 1, 'a re-annotated catalogue reruns the tables');
Test::More::ok(!-e $funcTablesSto, 'the stale checkpoint is removed before resubmitting');

# eggNOG-mapper finishes later: its scheme joins and the tables are rebuilt.
main::write_file("$funcTablesDir/KEGG.tsv", "\tMGS.1\n");
($checkpointCmd) = $submissions[0][1] =~ /^('perl' .*writeCheckpoint\.pl.*)$/m;
system('bash', '-c', $checkpointCmd) == 0 or die "checkpoint command failed";
main::write_file($stones->{emap}, "{\"stone\":1}\n");
reset_run(); @funcTableSchemes = qw(KEGG emNOG); $wait4stone = $stones->{emap};
_submit_function_tables();
Test::More::is(scalar(@submissions), 1, 'a newly finished annotation stage reruns the tables');
Test::More::like($submissions[0][1], qr/'-schemes' 'KEGG,emNOG'/, 'with both schemes');

# A new MAG report (re-clustering) invalidates the tables.
main::write_file("$funcTablesDir/emNOG.tsv", "\tMGS.1\n");
($checkpointCmd) = $submissions[0][1] =~ /^('perl' .*writeCheckpoint\.pl.*)$/m;
system('bash', '-c', $checkpointCmd) == 0 or die "checkpoint command failed";
reset_run(); _submit_function_tables();
Test::More::is(scalar(@submissions), 0, 'resumable once both schemes were made');
main::write_file("$outD/MAGvsGC.txt.gz", "new report with other MGS");
reset_run(); _submit_function_tables();
Test::More::is(scalar(@submissions), 1, 'a rewritten MAG report reruns the tables');

reset_run(); @funcTableSchemes = ();
_submit_function_tables();
Test::More::is(scalar(@submissions) + scalar(@log), 0, '-funcTables 0 does nothing');

Test::More::done_testing();
