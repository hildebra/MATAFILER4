use strict;
use warnings;

# -keepRawReads and the run-level accession manifest.
#
# Archive downloads normally live in <globalTmpDir>/<RunID>/<SmplID>/accession_download/
# and are removed with the sample temp directory. With -keepRawReads they must
# live outside every sample temp directory, so no cleanup path can reach them,
# and cleanup must never be handed them as --download-temp.

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use JSON::PP qw(decode_json encode_json);
use Test::More;

my $root = File::Spec->catdir($Bin, '..');
my $tmp = tempdir(CLEANUP => 1);

sub read_file { my ($p) = @_; open my $fh, '<', $p or die "$p: $!"; local $/; my $c = <$fh>; close $fh; return $c; }
sub write_file { my ($p, $c) = @_; open my $fh, '>', $p or die "$p: $!"; print {$fh} $c; close $fh or die $!; }
sub source_sub {
	my ($src, $name) = @_;
	my ($s) = $src =~ /^(sub \Q$name\E\b[^\n]*\{.*?^\})/ms;
	die "cannot find sub $name\n" unless defined $s;
	return $s;
}

my $mfSrc = read_file("$root/MATAF4.pl");

package RetentionHarness;
use File::Spec;
use JSON::PP qw(decode_json);
our (%map, %MFconfig, %MFglobal, $baseOut, $keptAccessionReadsDir, @written);
sub atomic_write_text { my ($path, $text) = @_; push @written, [$path, $text];
	open my $fh, '>', $path or die $!; print {$fh} $text; close $fh; }
package main;

for my $name (qw(sampleAccessionSource accessionDownloadDirFor writeAccessionManifest)) {
	my $code = source_sub($mfSrc, $name);
	eval "package RetentionHarness; $code; 1" or die "compile $name: $@";
}

$RetentionHarness::keptAccessionReadsDir = 'accession_reads';
$RetentionHarness::MFglobal{runTmpDirGlobal} = "$tmp/scratch/RUN1/";
$RetentionHarness::MFconfig{silent} = 1;

subtest 'download location' => sub {
	local $RetentionHarness::MFconfig{keepRawReads} = 0;
	my $default = RetentionHarness::accessionDownloadDirFor('S1');
	is(File::Spec->canonpath($default), File::Spec->canonpath("$tmp/scratch/RUN1/S1/accession_download"),
		'default location is unchanged: inside the sample temp directory');

	$RetentionHarness::MFconfig{keepRawReads} = 1;
	my $kept = RetentionHarness::accessionDownloadDirFor('S1');
	is(File::Spec->canonpath($kept), File::Spec->canonpath("$tmp/scratch/RUN1/accession_reads/S1"),
		'-keepRawReads derives a run-level location without user input');
	my $sampleTemp = File::Spec->canonpath("$tmp/scratch/RUN1/S1");
	unlike(File::Spec->canonpath($kept), qr/^\Q$sampleTemp\E(?:\/|$)/,
		'kept reads are outside the sample temp directory');
	ok(!eval { RetentionHarness::accessionDownloadDirFor('accession_reads'); 1 }
		&& $@ =~ /collides/, 'a SmplID equal to the kept-reads directory is refused');
};

subtest 'cleanup is never handed kept reads' => sub {
	my $args = source_sub($mfSrc, 'finishedCleanupArguments');
	like($args, qr/'--download-temp'.*\n\s*if \$accessionProvider ne '' && !\$MFconfig\{keepRawReads\}/,
		'--download-temp is only passed when reads are not kept');
	my @uses = grep { !/^\s*#/ && !/accession_download\.sh/ } grep { /accession_download/ } split /\n/, $mfSrc;
	is(scalar(@uses), 2, 'accession_download is named only by the default location and the cleanup argument')
		or diag explain \@uses;
	like($mfSrc, qr/\$accessionOnlyCleanup = \$accessionProvider ne ''\s*&& !\$MFconfig\{keepRawReads\}/,
		'no accession-only cleanup job is scheduled for kept reads');
	like($mfSrc, qr/"keepRawReads=i" => \\\$MFconfig\{keepRawReads\}/, '-keepRawReads is a command-line option');
	like($mfSrc, qr/\$MFconfig\{keepRawReads\} = 0;/, '-keepRawReads defaults to off');
};

subtest 'accession manifest' => sub {
	local $RetentionHarness::MFconfig{keepRawReads} = 1;
	my $out = "$tmp/out";
	make_path("$out/ENA1", "$out/LOCAL1");
	local $RetentionHarness::baseOut = $out;
	local %RetentionHarness::map = (
		ENA1 => { SmplID => 'ENA1', wrdir => "$out/ENA1", ENA_download => 'ERR1' },
		LOCAL1 => { SmplID => 'LOCAL1', wrdir => "$out/LOCAL1", ENA_download => '', SRA_download => '' },
	);
	write_file("$out/ENA1/input_accession.json", encode_json({
		schema => 1, provider => 'ena', accession => 'ERR1',
		downloaded_at => '2026-10-05T10:00:00Z',
		files => [
			{ run => 'ERR1', role => 'r1', filename => 'ERR1_1.fastq.gz', bytes => 10, md5 => 'a' x 32, records => 5 },
			{ run => 'ERR1', role => 'r2', filename => 'ERR1_2.fastq.gz', bytes => 11, md5 => 'b' x 32, records => 5 },
		],
	}));
	my $keptDir = RetentionHarness::accessionDownloadDirFor('ENA1');
	make_path($keptDir);
	write_file(File::Spec->catfile($keptDir, 'ERR1_1.fastq.gz'), 'x' x 10);

	@RetentionHarness::written = ();
	RetentionHarness::writeAccessionManifest([qw(ENA1 LOCAL1)]);
	is(scalar(@RetentionHarness::written), 1, 'one manifest is written');
	is($RetentionHarness::written[0][0], "$out/accession_manifest.tsv", 'manifest sits in the run output directory');
	my @lines = split /\n/, $RetentionHarness::written[0][1];
	is($lines[0], join("\t", qw(SmplID provider accession run role filename bytes md5 records downloaded_at path present retained)),
		'manifest header');
	is(scalar(@lines), 3, 'one row per archive file; local samples are not listed');
	my @r1 = split /\t/, $lines[1];
	is_deeply([@r1[0..5, 7, 8]], ['ENA1', 'ena', 'ERR1', 'ERR1', 'r1', 'ERR1_1.fastq.gz', 'a' x 32, 5],
		'row carries archive provenance and checksum');
	is($r1[11], 'yes', 'present file is reported as present');
	is((split /\t/, $lines[2])[11], 'no', 'missing file is reported as absent');
	is($r1[12], 'yes', 'retention policy is recorded');

	@RetentionHarness::written = ();
	RetentionHarness::writeAccessionManifest([qw(LOCAL1)]);
	is(scalar(@RetentionHarness::written), 0, 'runs without accession samples write no manifest');
};

subtest 'finished-sample cleanup script is unchanged by this feature' => sub {
	my $cleanup = read_file("$root/secScripts/cleanup_finished_sample.pl");
	unlike($cleanup, qr/keepRawReads|accession_reads/, 'cleanup_finished_sample.pl has no retention logic of its own');
};

done_testing();
