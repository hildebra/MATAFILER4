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

my $help = qx{"$binary" -help 2>&1};
like($help, qr/minOccupancy FLOAT.*0\.35/s,
	'help documents the permissive metagenomic occupancy default');
like($help, qr/maxP90Divergence FLOAT.*NT 0\.3/s,
	'help documents the permissive nucleotide divergence guard');

done_testing();
