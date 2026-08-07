use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

my $root = File::Spec->catdir($Bin, '..');
my $binary = File::Spec->catfile($root, 'bin', 'MSAfix');
my $tmp = tempdir(CLEANUP => 1);

sub write_alignment {
	my ($name, $sequences) = @_;
	my $path = File::Spec->catfile($tmp, "$name.fna");
	open my $fh, '>', $path or die "Cannot create $path: $!";
	my $index = 0;
	for my $sequence (@{$sequences}) {
		$index++;
		print {$fh} ">sample${index}|$name\n$sequence\n";
	}
	close $fh or die "Cannot close $path: $!";
	return $path;
}

sub sequence_with_c_block {
	my ($length, $offset, $count) = @_;
	my $sequence = 'A' x $length;
	substr($sequence, $offset, $count) = 'C' x $count;
	return $sequence;
}

my @alignments;
for my $index (1 .. 8) {
	push @alignments, write_alignment(
		"stable$index",
		[
			'A' x 120,
			'A' x 119 . 'C',
			'A' x 120,
			'A' x 120,
		],
	);
}
my $permissive = write_alignment(
	'permissive_strain_signal',
	[
		'A' x 120,
		'C' x 6 . 'A' x 114,
		'A' x 6 . 'C' x 6 . 'A' x 108,
		'-' x 60 . 'A' x 60,
	],
);
my $wrong_orthologue = write_alignment(
	'wrong_orthologue',
	['A' x 120, 'A' x 120, 'A' x 120, 'C' x 120],
);
my $low_occupancy = write_alignment(
	'low_occupancy',
	[('A' x 20 . '-' x 100) x 4],
);
push @alignments, $permissive, $wrong_orthologue, $low_occupancy;

my $manifest = File::Spec->catfile($tmp, 'manifest.txt');
my $report = File::Spec->catfile($tmp, 'report.tsv');
my $keep = File::Spec->catfile($tmp, 'keep.txt');
open my $manifest_fh, '>', $manifest or die "Cannot create $manifest: $!";
print {$manifest_fh} "$_\n" for @alignments;
close $manifest_fh or die "Cannot close $manifest: $!";

my $status = system(
	$binary,
	'-manifest', $manifest,
	'-report', $report,
	'-keep', $keep,
	'-sequenceType', 'nt',
);
is($status, 0, 'native MSAfix post-alignment locus QC completes');
ok(-s $report, 'QC writes an auditable per-locus report');
ok(-s $keep, 'QC writes a retained-locus manifest');
my @report_temporary = glob($report.'.tmp.*');
my @keep_temporary = glob($keep.'.tmp.*');
is(scalar(@report_temporary), 0, 'native QC leaves no report temporary output');
is(scalar(@keep_temporary), 0, 'native QC leaves no keep-list temporary output');

open my $keep_fh, '<', $keep or die "Cannot read $keep: $!";
chomp(my @kept = <$keep_fh>);
close $keep_fh;
is(scalar(@kept), 9, 'stable loci and a modest strain-divergence locus are retained');
ok(grep($_ eq $permissive, @kept),
	'a partially missing locus with modest strain divergence passes permissive defaults');
ok(!grep($_ eq $wrong_orthologue, @kept),
	'a grossly non-comparable orthologue is rejected');
ok(!grep($_ eq $low_occupancy, @kept),
	'a locus dominated by missing alignment cells is rejected');

open my $report_fh, '<', $report or die "Cannot read $report: $!";
my $report_text = do { local $/; <$report_fh> };
close $report_fh;
like($report_text,
	qr/\Q$wrong_orthologue\E\tREJECT\t[^\n]*high_p90_consensus_divergence/,
	'report explains the extreme sequence-divergence rejection');
like($report_text, qr/\Q$low_occupancy\E\tREJECT\t[^\n]*low_occupancy/,
	'report explains the low-occupancy rejection');
like($report_text, qr/\Q$permissive\E\tPASS\t\./,
	'report records retained strain-variable loci');

my $broad_manifest = File::Spec->catfile($tmp, 'broad-aa-manifest.txt');
my $broad_report = File::Spec->catfile($tmp, 'broad-aa-report.tsv');
my $broad_keep = File::Spec->catfile($tmp, 'broad-aa-keep.txt');
open my $broad_manifest_fh, '>', $broad_manifest
	or die "Cannot create $broad_manifest: $!";
print {$broad_manifest_fh} "$wrong_orthologue\n";
close $broad_manifest_fh or die "Cannot close $broad_manifest: $!";
my $broad_status = system(
	$binary,
	'-manifest', $broad_manifest,
	'-report', $broad_report,
	'-keep', $broad_keep,
	'-sequenceType', 'aa',
	'-maxMedianDivergence', 1,
	'-maxP90Divergence', 1,
	'-relativeModifiedZ', 1_000_001,
);
is($broad_status, 0, 'broad-AA locus QC profile completes');
open my $broad_keep_fh, '<', $broad_keep
	or die "Cannot read $broad_keep: $!";
chomp(my @broad_kept = <$broad_keep_fh>);
close $broad_keep_fh;
is_deeply(\@broad_kept, [$wrong_orthologue],
	'broad-AA profile retains a structurally valid locus despite deep between-species divergence');

my @relative_loci;
for my $index (1 .. 8) {
	my $third_distance = (70, 70, 90, 90, 110, 110, 130, 130)[$index - 1];
	my $fourth_distance = (156, 156, 161, 161, 167, 167, 173, 173)[$index - 1];
	push @relative_loci, write_alignment(
		"relative_stable$index",
		[
			'A' x 1000,
			sequence_with_c_block(1000, 0, 10),
			sequence_with_c_block(1000, 200, $third_distance),
			sequence_with_c_block(1000, 500, $fourth_distance),
		],
	);
}
my $relative_outlier = write_alignment(
	'relative_rate_outlier',
	[
		'A' x 1000,
		sequence_with_c_block(1000, 0, 10),
		sequence_with_c_block(1000, 200, 200),
		sequence_with_c_block(1000, 500, 271),
	],
);
my $relative_manifest = File::Spec->catfile($tmp, 'relative-manifest.txt');
open my $relative_manifest_fh, '>', $relative_manifest
	or die "Cannot create $relative_manifest: $!";
print {$relative_manifest_fh} "$_\n" for @relative_loci, $relative_outlier;
close $relative_manifest_fh or die "Cannot close $relative_manifest: $!";

my $permissive_relative_report = File::Spec->catfile($tmp, 'relative-z8-report.tsv');
my $permissive_relative_keep = File::Spec->catfile($tmp, 'relative-z8-keep.txt');
my $permissive_relative_status = system(
	$binary,
	'-manifest', $relative_manifest,
	'-report', $permissive_relative_report,
	'-keep', $permissive_relative_keep,
	'-sequenceType', 'nt',
	'-relativeModifiedZ', 8,
);
is($permissive_relative_status, 0, 'permissive relative-divergence QC completes');
open my $permissive_relative_keep_fh, '<', $permissive_relative_keep
	or die "Cannot read $permissive_relative_keep: $!";
chomp(my @permissive_relative_kept = <$permissive_relative_keep_fh>);
close $permissive_relative_keep_fh;
ok(grep($_ eq $relative_outlier, @permissive_relative_kept),
	'an intermediate divergence locus passes the former Z=8 cutoff');

my $strict_relative_report = File::Spec->catfile($tmp, 'relative-z5-report.tsv');
my $strict_relative_keep = File::Spec->catfile($tmp, 'relative-z5-keep.txt');
my $strict_relative_status = system(
	$binary,
	'-manifest', $relative_manifest,
	'-report', $strict_relative_report,
	'-keep', $strict_relative_keep,
	'-sequenceType', 'nt',
	'-relativeModifiedZ', 5,
);
is($strict_relative_status, 0, 'strict relative-divergence QC completes');
open my $strict_relative_keep_fh, '<', $strict_relative_keep
	or die "Cannot read $strict_relative_keep: $!";
chomp(my @strict_relative_kept = <$strict_relative_keep_fh>);
close $strict_relative_keep_fh;
is(scalar(@strict_relative_kept), scalar(@relative_loci),
	'the stricter cutoff removes exactly the divergence-outlier locus');
ok(!grep($_ eq $relative_outlier, @strict_relative_kept),
	'the stricter cutoff excludes the anomalously fast locus rather than a sample');
open my $strict_relative_report_fh, '<', $strict_relative_report
	or die "Cannot read $strict_relative_report: $!";
my $strict_relative_text = do { local $/; <$strict_relative_report_fh> };
close $strict_relative_report_fh;
like($strict_relative_text,
	qr/\Q$relative_outlier\E\tREJECT\t[^\n]*relative_p90_divergence_outlier/,
	'the strict report records the locus-level P90 divergence reason');

my $help = qx{"$binary" -help 2>&1};
like($help, qr/minOccupancy FLOAT.*0\.35/s,
	'help documents the permissive metagenomic occupancy default');
like($help, qr/maxP90Divergence FLOAT.*NT 0\.3/s,
	'help documents the permissive nucleotide divergence guard');

done_testing();
