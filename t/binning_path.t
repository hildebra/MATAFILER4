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
	MB2N50 binningOutputsComplete runCheckM runCheckM2 runSemiBin
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

my ($empty_command, $empty_alignments) = Mods::Binning::createBams(
	[], "$root/tmp", "$root/empty", 'empty', "$root/ref.fa", 2, 1, 0, 'bam',
);
like($empty_command, qr/: > \Q$root\/empty\/empty\E/,
	'no-mapping handling returns an explicit empty-assignment command');
is_deeply($empty_alignments, [], 'no-mapping handling returns no fabricated alignments');
ok(!-e "$root/empty/empty", 'command construction does not publish outputs prematurely');

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
like(checkMapsDoneSH([$sample_dir]), qr/\Q$sample_dir\/mapping\/done.sto\E/,
	'mapping checks recognize existing sample directories without requiring a trailing slash');

my $direct_bam = "$root/direct.bam";
write_file($direct_bam, 'bam');
my ($semibin_default, $semibin_environment);
{
	no warnings 'redefine';
	local *Mods::Binning::getProgPaths = sub { return 'SemiBin2' };
	$semibin_default = runSemiBin('', "$root/sb", "$root/sbtmp", 'sample',
		"$root/ref.fa", 4, [$direct_bam], 'hiSeq', '');
	$semibin_environment = runSemiBin('', "$root/sb2", "$root/sbtmp2", 'sample',
		"$root/ref.fa", 4, [$direct_bam], 'ONT', 'ocean');
}
unlike($semibin_default, qr/--environment/,
	'an empty SB_env uses SemiBin self-training instead of silently forcing human_gut');
like($semibin_environment, qr/--environment ocean/,
	'an explicitly selected SemiBin environment is passed through');
like($semibin_environment, qr/--sequencing-type=long_read/,
	'long-read libraries select SemiBin long-read mode');

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

my $mataf4 = slurp(File::Spec->catfile($Bin, '..', 'MATAF4.pl'));
like($mataf4, qr/my \$binningComplete = binningOutputsComplete.*?&& !\$binningComplete/s,
	'MATAF4 uses the shared complete-output predicate when scheduling binning');
like($mataf4, qr/\$MBcmd = "" if \(-e \$MetaBat2out\)/,
	'quality-only repair reuses a valid empty bin assignment');
like($mataf4, qr/my \@binLibraries = .*?getRawLibrariesAssmGrp.*?0.*?getRawLibrariesAssmGrp.*?1/s,
	'binner sequencing mode is derived from primary and support library records');
unlike($mataf4, qr/\$postCmd \.= " -read[12S] /,
	'MATAF4 no longer passes ignored raw-read options to bin quality checking');

my $runner = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'assemblies', 'runBinners.pl'));
like($runner, qr/my \$resultCheck = "test -e '\$BinDir\/\$smplIDs1'".*?printf .*?> '\$stone'/s,
	'the binner assignment is verified before its completion stone is published');
like($runner, qr/if \(\$DoMetaBat2 == 2\).*?output_recluster_bins.*?output_bins/s,
	'SemiBin native output is accepted before standardized assignments are generated');

done_testing();
