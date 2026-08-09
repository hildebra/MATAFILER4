use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::Binning ();
use Mods::Binning qw(
	MB2N50 binningOutputsComplete emptyBinnerAssignmentCommand
	runCheckM runCheckM2 runSCGBinner runSemiBin
);
use Mods::IO_Tamoc_progs qw(checkMapsDoneSH jgi_depth_cmd);

sub write_file {
	my ($path, $contents) = @_;
	open my $fh, '>', $path or die "Cannot write $path: $!";
	print {$fh} $contents;
	close $fh or die "Cannot close $path: $!";
}

sub slurp {
	my ($path) = @_;
	open my $fh, '<', $path or die "Cannot read $path: $!";
	local $/;
	return <$fh>;
}

my $root = tempdir(CLEANUP => 1);
$root =~ s{\\}{/}g;

my $base = "$root/sample";
write_file($base, '');
write_file("$base.assStat", "header\n");
write_file("$base.cm", '');
write_file("$base.cm2", "Name\tCompleteness\tContamination\n");
ok(binningOutputsComplete($base, 1, 1),
	'binning completion requires and accepts both requested quality reports');
unlink "$base.cm2" or die $!;
ok(!binningOutputsComplete($base, 1, 1),
	'one quality report cannot complete a run when both were requested');
ok(binningOutputsComplete($base, 1, 0),
	'unrequested quality formats do not block completion');

my $stats = MB2N50({bin1 => [
	'ctg1_L=40=', 'ctg2_L=30=', 'ctg3_L=20=', 'ctg4_L=10=',
]});
is($stats->{bin1}{N20}, 40, 'N20 is calculated from longest to shortest contig');
is($stats->{bin1}{N80}, 20, 'N80 is calculated from longest to shortest contig');
eval { MB2N50({bad => ['contig-without-length']}) };
like($@, qr/Cannot determine contig length/, 'malformed contig identifiers cannot reuse a stale regex capture');

my $sample_dir = "$root/S1";
make_path("$sample_dir/mapping");
write_file("$sample_dir/mapping/done.sto", "S1-smd.bam\n");
write_file("$sample_dir/mapping/S1-smd.bam", 'p');
write_file("$sample_dir/mapping/S1.sup-smd.bam", 's');
my ($conversion, $alignments) = Mods::Binning::createBams(
	[$sample_dir], "$root/tmp", "$root/out", 'bins', "$root/ref.fa", 2, 0, 0, 'bam',
);
is($conversion, '', 'small completed BAM files do not require conversion');
is_deeply($alignments, [
	"$sample_dir/mapping/S1-smd.bam",
	"$sample_dir/mapping/S1.sup-smd.bam",
], 'primary and supplemental mappings are both supplied to BAM-based binners');

my $support_only_dir = "$root/S2";
make_path("$support_only_dir/mapping");
write_file("$support_only_dir/mapping/done.sto", "S2.sup-smd.bam\n");
write_file("$support_only_dir/mapping/S2.sup-smd.bam", 's');
my ($support_conversion, $support_alignments) = Mods::Binning::createBams(
	[$support_only_dir], "$root/support-tmp", "$root/support-out", 'support-bins',
	"$root/ref.fa", 2, 0, 0, 'bam',
);
is($support_conversion, '', 'support-only BAM does not require conversion');
is_deeply($support_alignments, ["$support_only_dir/mapping/S2.sup-smd.bam"],
	'support-only marker supplies its mapping without constructing a sup.sup filename');

my ($empty_command, $empty_alignments) = Mods::Binning::createBams(
	[], "$root/tmp", "$root/empty", 'empty', "$root/ref.fa", 2, 1, 0, 'bam',
);
like($empty_command, qr/: > \Q$root\/empty\/empty\E/,
	'no-mapping handling returns an explicit empty-assignment command');
is_deeply($empty_alignments, [], 'no-mapping handling returns no fabricated alignments');
ok(!-e "$root/empty/empty", 'command construction does not publish outputs prematurely');
is($empty_command, emptyBinnerAssignmentCommand("$root/empty", 'empty'),
	'no-mapping and undersized-assembly handling share the empty assignment command');

my $depth_command;
{
	no warnings 'redefine';
	local *Mods::IO_Tamoc_progs::getProgPaths = sub {
		return $_[0] eq 'samtools' ? 'samtools' : 'jgi_summarize_bam_contig_depths';
	};
	$depth_command = jgi_depth_cmd([$sample_dir], "$root/depth", 95, 2, "$root/ref.fa");
}
like($depth_command, qr/\Q$sample_dir\/mapping\/S1-smd.bam\E/,
	'JGI depth includes the primary alignment');
like($depth_command, qr/\Q$sample_dir\/mapping\/S1.sup-smd.bam\E/,
	'JGI depth includes the supplemental alignment');
like($depth_command, qr/^set -e\n/,
	'JGI depth generation cannot hide a failed CRAM conversion or depth program');
my $support_depth_command;
{
	no warnings 'redefine';
	local *Mods::IO_Tamoc_progs::getProgPaths = sub {
		return $_[0] eq 'samtools' ? 'samtools' : 'jgi_summarize_bam_contig_depths';
	};
	$support_depth_command = jgi_depth_cmd(
		[$support_only_dir], "$root/support-depth", 95, 2, "$root/ref.fa");
}
like($support_depth_command, qr/\Q$support_only_dir\/mapping\/S2.sup-smd.bam\E/,
	'JGI depth accepts a support-only assembly mapping');
unlike($support_depth_command, qr/\.sup\.sup-smd\./,
	'support-only depth discovery does not invent a duplicate support suffix');
like(checkMapsDoneSH([$sample_dir]), qr/\Q$sample_dir\/mapping\/done.sto\E/,
	'mapping checks recognize existing sample directories without requiring a trailing slash');

my $direct_bam = "$root/direct.bam";
write_file($direct_bam, 'bam');
my $scgbinner_command;
{
	no warnings 'redefine';
	local *Mods::Binning::getProgPaths = sub { return 'scgbinner' };
	$scgbinner_command = runSCGBinner(
		'', "$root/scg-out", "$root/scg-tmp", 'sample',
		"$root/ref.fa", 4, [$direct_bam], 17, 42,
	);
}
like($scgbinner_command,
	qr/scgbinner .*?-b "\Q$direct_bam\E" -t 4 -p 17/,
	'SCGBinner command uses the preflight-derived batch size and quotes its BAM list');
is(system('bash', '-n', '-c', $scgbinner_command), 0,
	'SCGBinner failure-diagnostic command is valid Bash');

my $large_bam = "$root/large.bam";
open my $large_fh, '>', $large_bam or die "Cannot write $large_bam: $!";
seek($large_fh, 15 * 1024 * 1024, 0) or die "Cannot seek in $large_bam: $!";
print {$large_fh} 'x';
close $large_fh or die "Cannot close $large_bam: $!";
my $large_bam2 = "$root/large2.bam";
open my $large_fh2, '>', $large_bam2 or die "Cannot write $large_bam2: $!";
seek($large_fh2, 15 * 1024 * 1024, 0) or die "Cannot seek in $large_bam2: $!";
print {$large_fh2} 'x';
close $large_fh2 or die "Cannot close $large_bam2: $!";
my ($semibin_small, $semibin_default, $semibin_environment, $semibin_multi);
{
	no warnings 'redefine';
	local *Mods::Binning::getProgPaths = sub { return 'SemiBin2' };
	$semibin_small = runSemiBin('', "$root/sb-small", "$root/sbtmp-small", 'sample',
		"$root/ref.fa", 4, [$direct_bam], 'hiSeq', '');
	$semibin_default = runSemiBin('', "$root/sb", "$root/sbtmp", 'sample',
		"$root/ref.fa", 4, [$large_bam], 'hiSeq', '');
	$semibin_environment = runSemiBin('', "$root/sb2", "$root/sbtmp2", 'sample',
		"$root/ref.fa", 4, [$large_bam], 'ONT', 'ocean');
	$semibin_multi = runSemiBin('', "$root/sb-multi", "$root/sbtmp-multi", 'sample',
		"$root/ref.fa", 4, [$large_bam, $large_bam2], 'hiSeq', 'ocean');
}
like($semibin_small, qr/: > \Q$root\/sb-small\/sample\E/,
	'SemiBin excludes BAM files of 15 MiB or less and publishes an empty assignment');
unlike($semibin_small, qr/SemiBin2 single_easy_bin/,
	'SemiBin is not invoked with a crash-prone small BAM');
like($semibin_default, qr/--environment human_gut/,
	'an empty SB_env defaults SemiBin to the curated human-gut environment');
like($semibin_environment, qr/--environment ocean/,
	'an explicitly selected SemiBin environment is passed through');
like($semibin_environment, qr/--sequencing-type=long_read/,
	'long-read libraries select SemiBin long-read mode');
unlike($semibin_multi, qr/--environment/,
	'multiple usable BAMs omit pretrained environments for SemiBin multi-sample training');
eval {
	runSemiBin('', "$root/sb-invalid", "$root/sbtmp-invalid", 'sample',
		"$root/ref.fa", 4, [$large_bam], 'hiSeq', 'human_gut; touch BAD');
};
like($@, qr/Invalid SemiBin2 environment/,
	'SemiBin environment overrides cannot inject shell commands');

my $assignment = "$root/check-bins";
write_file($assignment, "contig\tbin1\n");
my ($checkm_command, $checkm2_command);
{
	no warnings 'redefine';
	local *Mods::Binning::getProgPaths = sub { return $_[0] };
	$checkm_command = runCheckM("$root/bins", "$assignment.cm", "$root/checkm", 3, 0, 'fna');
	$checkm2_command = runCheckM2("$root/bins", "$assignment.cm2", "$root/checkm2", 3, 0, 'fna');
}
like($checkm_command, qr/^set -e\n/, 'CheckM command propagation stops at the first failed program');
like($checkm_command, qr/checkm lineage_wf/, 'a single headerless assigned bin is quality checked');
like($checkm2_command, qr/^set -e\n/, 'CheckM2 command propagation stops at the first failed program');
like($checkm2_command, qr/test -s \Q$assignment.cm2\E/,
	'CheckM2 must publish a non-empty quality report');

done_testing();
