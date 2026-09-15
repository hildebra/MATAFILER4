use strict;
use warnings;
use Test::More;
use File::Temp qw(tempdir);
use File::Path qw(make_path);
use File::Spec;
use FindBin qw($Bin);
use lib "$Bin/..", "$Bin/lib";
use MFTestConfig;

sub write_file {
	my ($path, $text) = @_;
	open my $fh, '>', $path or die "$path: $!";
	print {$fh} $text or die "$path: $!";
	close $fh or die "$path: $!";
}
sub slurp {
	open my $fh, '<', $_[0] or die "$_[0]: $!";
	local $/;
	my $text = <$fh>;
	close $fh or die $!;
	return $text;
}

# Execute the controller's actual option setup, Canopy helper, and abundance
# submission block. Only site configuration and the external scheduler/QC tools
# are substituted; filtering, protein extraction, and checkpoints run normally.
{
	package MGSFixes;
	use Getopt::Long qw(GetOptions);
	use Mods::Binning qw(createBinFAA);
	use Mods::Checkpoint qw(read_checkpoint checkpoint_valid write_checkpoint);
	use Mods::WorkflowResilience qw(retry_unlink retry_rename);
	our ($GCd, $clusterID, $cmSuffix, $rewrClusterMAGs, $useCheckM2,
		$canCore, $numCore, $nodeTmpD, $logDir, $failQC);
	our (%QSBopt, @submissions, @qualityArgs, %paths);
	sub getProgPaths { return $paths{$_[0]} // $_[0]; }
	sub checkMF {}
	sub printL {}
	sub _checkpoint_valid { return 0; }
	sub _checkpoint_command { return "CHECKPOINT\n"; }
	sub runCheckM { @qualityArgs = @_; return 'checkm1'; }
	sub runCheckM2 { @qualityArgs = @_; return 'checkm2'; }
	sub qsubSystem {
		my @args = @_;
		push @submissions, \@args;
		if ($args[1] =~ /^checkm[12]$/ && !$failQC) {
			main::write_file($qualityArgs[1], "Name\tCompleteness\tContamination\nCanoA\t95\t1\n");
		}
		return ('123', $args[1]);
	}
	sub qsubSystemJobAlive {}
}
my $harness = 'package MGSFixes; our ($GCd, $clusterID, $cmSuffix, $rewrClusterMAGs, $useCheckM2, $canCore, $numCore, $nodeTmpD, $logDir); our %QSBopt;';
my $source = slurp("$Bin/../secScripts/MGS.pl");
my ($option_setup) = $source =~ /(my \$inD = .*?)(?=\n%checkpointParameters =)/s;
die 'Cannot isolate MGS option setup' unless defined $option_setup;
my $parse = eval 'package MGSFixes; sub { local @ARGV = @_; my $clusterID = 95; '
	. $option_setup . '; return [0 + $rewrClusterMAGs, 0 + $rewrTAX, $strainRedo]; }';
die $@ unless $parse;
for my $case (
	['none', 0, 0, 'none'], ['cluster', 1, 0, 'none'],
	['tax', 0, 1, 'none'], ['tree', 0, 0, 'tree'],
	['input', 0, 0, 'input'], ['all', 1, 1, 'all'],
) {
	my ($mode, @expected) = @$case;
	is_deeply($parse->('-redo', $mode), \@expected, "-redo $mode selects the intended rebuild stages");
}
is_deeply($parse->(), [0, 0, 'none'], 'default mode preserves existing work');
for my $old (qw(redoCluster redoTax)) {
	my $warning = '';
	local $SIG{__WARN__} = sub { $warning .= join '', @_; };
	eval { $parse->("-$old", 1); };
	like($@, qr/Invalid MGS.pl options/, "-$old is rejected");
	like($warning, qr/Unknown option: \L$old\E/i, "-$old reports an unknown option");
}
eval { $parse->('-redo', 'typo'); };
like($@, qr/-redo must be one of: none, cluster, tax, tree, input, all/, 'invalid modes list the consolidated choices');

my ($canopy_helper) = $source =~ /(sub CanopyPrep\{.*?)(?=\nsub getGoodMBstats)/s;
my ($count_helpers) = $source =~ /(sub _mgs_ids \{.*?)(?=\nsub _write_single_mgs_observations)/s;
die 'Cannot isolate Canopy helper' unless defined $canopy_helper && defined $count_helpers;
eval $harness . $count_helpers . $canopy_helper;
die $@ if $@;

my $tmp = tempdir(CLEANUP => 1);
$MGSFixes::GCd = $tmp;
$MGSFixes::clusterID = 95;
$MGSFixes::cmSuffix = '.cm2';
$MGSFixes::useCheckM2 = 1;
$MGSFixes::rewrClusterMAGs = 0;
$MGSFixes::canCore = 7;
$MGSFixes::numCore = 3;
$MGSFixes::nodeTmpD = $tmp;
$MGSFixes::logDir = $tmp;
my $input = "$tmp/canopies.txt";
my $bins = "$tmp/Bins";
my $protein = "$tmp/compl.incompl.95.prot.faa";
my $stone = "$input.filt.checkpoint";
make_path($bins);
write_file($protein, join '', map { ">$_\nMACD\n" } 1..700);
write_file($input, join '', map { "CanoA\t$_\n" } 1..700);
write_file("$input.filt", "OLD\t1\n");
write_file("$input.filt.cm2", "obsolete QC\n");
write_file("$bins/OLD.faa", ">old\nM\n");
write_file("$bins/notes.txt", "keep\n");
is(MGSFixes::CanopyPrep($input, $bins), 1, 'legacy Canopy cache is rebuilt');
like(slurp("$input.filt"), qr/^CanoA\t700$/m, 'rebuilt filter contains current assignments');
ok(-s "$bins/CanoA.faa" && !-e "$bins/OLD.faa", 'QC input bins contain the current generation');
is(slurp("$bins/notes.txt"), "keep\n", 'Canopy cleanup preserves unrelated files');
ok(-s $stone, 'successful preparation publishes a manifest');
is(scalar @MGSFixes::submissions, 1, 'initial preparation runs one QC job');
is_deeply([@{$MGSFixes::submissions[-1]}[1..3]], ['checkm2', 7, '50G'], 'CheckM2 retains its cores and total memory');
is(MGSFixes::CanopyPrep($input, $bins), 1, 'matching Canopy cache is reusable');
is(scalar @MGSFixes::submissions, 1, 'matching inputs skip QC');

# A same-size rewrite with a changed timestamp must invalidate the manifest.
write_file($input, join '', map { "CanoB\t$_\n" } 1..700);
my $source_mtime = (stat $input)[9];
utime($source_mtime + 5, $source_mtime + 5, $input) or die $!;
is(MGSFixes::CanopyPrep($input, $bins), 1, 'changed source assignments rebuild the cache');
ok(-s "$bins/CanoB.faa" && !-e "$bins/CanoA.faa", 'removed Canopy is absent from the new QC input');
is(scalar @MGSFixes::submissions, 2, 'changed source reruns QC');
write_file($protein, join '', map { ">$_\nMCHANGED\n" } 1..700);
MGSFixes::CanopyPrep($input, $bins);
is(scalar @MGSFixes::submissions, 3, 'changed protein catalog reruns QC');
like(slurp("$bins/CanoB.faa"), qr/MCHANGED/, 'protein bins use the changed catalog');
$MGSFixes::rewrClusterMAGs = 1;
MGSFixes::CanopyPrep($input, $bins);
is(scalar @MGSFixes::submissions, 4, 'explicit cluster redo also rebuilds Canopy QC');
$MGSFixes::rewrClusterMAGs = 0;
$MGSFixes::useCheckM2 = 0;
$MGSFixes::cmSuffix = '.cm';
MGSFixes::CanopyPrep($input, $bins);
is_deeply([@{$MGSFixes::submissions[-1]}[1..3]], ['checkm1', 3, '200G'], 'checker change runs CheckM1 with its configured resources');
$MGSFixes::paths{checkm} = 'updated-checkm';
MGSFixes::CanopyPrep($input, $bins);
is(scalar @MGSFixes::submissions, 6, 'changed configured checker command invalidates QC');
unlink "$input.filt.cm" or die $!;
MGSFixes::CanopyPrep($input, $bins);
is(scalar @MGSFixes::submissions, 7, 'missing quality output invalidates QC');
write_file("$input.filt", '');
MGSFixes::CanopyPrep($input, $bins);
is(scalar @MGSFixes::submissions, 8, 'damaged filtered assignments invalidate the cache');

$MGSFixes::rewrClusterMAGs = 1;
$MGSFixes::failQC = 1;
eval { MGSFixes::CanopyPrep($input, $bins); };
like($@, qr/quality checking completed without producing/, 'failed QC is reported');
ok(!-e $stone, 'failed QC leaves no reusable manifest');
$MGSFixes::failQC = 0;
$MGSFixes::rewrClusterMAGs = 0;
MGSFixes::CanopyPrep($input, $bins);
is(scalar @MGSFixes::submissions, 10, 'resume retries failed QC');
write_file($input, "tiny\t1\n");
is(MGSFixes::CanopyPrep($input, $bins), 0, 'empty filtering result completes without QC');
ok(-e "$input.filt" && !-s "$input.filt" && -s $stone, 'zero-Canopy result is checkpointed');
ok(!-e "$input.filt.cm" && !-e "$bins/CanoB.faa", 'zero-Canopy result removes stale QC and protein bins');
is(MGSFixes::CanopyPrep($input, $bins), 0, 'unchanged zero-Canopy result resumes');
is(scalar @MGSFixes::submissions, 10, 'zero-Canopy cache does not run QC');

my ($abundance) = $source =~ /(my \$specIoutDir = .*?)(?=\nqsubSystemJobAlive\( \\\@annotation_jobs)/s;
die 'Cannot isolate abundance submission' unless defined $abundance;
my $submit_abundance = eval $harness . 'sub {
	my ($legacyV, $annoDir, $COGdir, $ABmgsSton, $finalClustersFilt, $GTDBtaxF,
		$useGTDBmg, $checkpointWriter) = (0, "anno", "GTDBmg", "ab.stone", "core", "tax", "GTDB", "writer");'
	. $abundance . '}';
die $@ unless $submit_abundance;
$submit_abundance->();
my $job = $MGSFixes::submissions[-1];
like($job->[1], qr/-cores 7\b/, 'abundance worker receives requested bottleneck cores');
is($job->[2], 7, 'abundance scheduler allocation matches worker cores');
is($job->[3], '64G', 'abundance memory remains a total 64G request');

done_testing();
