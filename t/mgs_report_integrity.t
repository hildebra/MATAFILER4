use strict;
use warnings;
use Test::More;
use File::Temp qw(tempdir);
use File::Path qw(make_path);
use FindBin qw($Bin);
use IO::Compress::Gzip qw(gzip $GzipError);
use IPC::Open3 qw(open3);
use Symbol qw(gensym);
use lib "$Bin/..", "$Bin/lib";
use MFTestConfig;
use Mods::Binning ();

sub write_file {
	my ($path, $text) = @_;
	open my $fh, '>', $path or die "$path: $!";
	binmode $fh;
	print {$fh} $text or die "$path: $!";
	close $fh or die "$path: $!";
}
my $tmp = tempdir(CLEANUP => 1);
my @header = qw(MAG MGS Representative4MGS Match2MGS Uniqueness AssociatedMGS
	Completeness Contamination LCAcompleteness N50 N_Genes GC CodingDensity
	CentreScore CompoundScore Domain Phylum Class Order Family Genus Species COG1 other_genes);
my $report = join("\t", @header) . "\n"
	. join("\t", 'S1__binA', 'MGS1', '*', 1, 1, '', 95, 1, 90, 1000,
		700, 50, 90, 1, 1, ('?') x 7, 1, '') . "\n";
my $compressed;
gzip(\$report => \$compressed) or die $GzipError;
my $valid = "$tmp/report.gz";
my $broken = "$tmp/truncated.gz";
write_file($valid, $compressed);
# All table rows are still readable, but the CRC/size trailer is absent. Only
# checking the decompressor's close status detects this corruption.
write_file($broken, substr($compressed, 0, -8));
my $map = { opt => { smpl_order => ['S1'] }, S1 => { FamGroup => 'Fam', AssGroup => '' } };
for my $case (
	['representatives', sub { Mods::Binning::getRepresentBins($_[0]) }, { MGS1 => 'S1__binA' }],
	['family representatives', sub { Mods::Binning::getRepresentBinsPerFamily($_[0], $map) }, { 'Fam.MGS1' => 'S1__binA' }],
) {
	my ($name, $read, $expected) = @$case;
	is_deeply($read->($valid), $expected, "$name accept a complete gzip report");
	my $result = eval { $read->($broken) };
	like($@, qr/Cannot finish reading MAG report \Q$broken\E/, "$name reject a truncated gzip report");
	ok(!defined $result, "$name do not return partial data after decompression failure");
}
my $gc = "$tmp/GC";
my $bin_dir = "$gc/Bin_SB";
make_path("$bin_dir/Annotation");
write_file("$bin_dir/SB.clusters.core", "MGS1\t1\n");
write_file("$bin_dir/Annotation/marker2MGS.txt", "previous complete result\n");
my $err = gensym;
my $pid = open3(undef, my $out, $err, $^X, "-I$Bin/..", "-I$Bin/lib", '-MMFTestConfig',
	"$Bin/../secScripts/MGS/markersPerMGS.pl", '-GCd', $gc, '-MAGlogFile', $broken);
my $stdout = do { local $/; <$out> // '' };
my $stderr = do { local $/; <$err> // '' };
waitpid($pid, 0);
ok(($? >> 8) != 0, 'marker worker fails on a truncated gzip report');
like($stderr, qr/Cannot finish reading MAG report \Q$broken\E/, 'marker worker reports decompression failure at the input boundary');
open my $previous, '<', "$bin_dir/Annotation/marker2MGS.txt" or die $!;
is(do { local $/; <$previous> }, "previous complete result\n", 'marker worker preserves the previous complete output');
close $previous or die $!;
ok(!-e "$bin_dir/Annotation/marker2MGS.LCA.txt", 'marker worker does not publish new derivatives from a corrupt report');
done_testing();
