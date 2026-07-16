use strict;
use warnings;

use File::Find;
use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::WorkflowState qw(inspect_workflow_state encode_state_report);

sub write_file {
	my ($path, $contents) = @_;
	(my $dir = $path) =~ s{/[^/]+$}{};
	make_path($dir);
	open my $fh, '>', $path or die "Cannot write $path: $!";
	print {$fh} $contents;
	close $fh;
}

sub filesystem_snapshot {
	my ($root) = @_;
	my @entries;
	find(sub {
		return if ($File::Find::name eq $root);
		my $relative = File::Spec->abs2rel($File::Find::name, $root);
		my $size = -f $File::Find::name ? (-s $File::Find::name || 0) : -1;
		push @entries, "$relative:$size";
	}, $root);
	return [sort @entries];
}

my $root = tempdir(CLEANUP => 1);
$root =~ s{\\}{/}g;
my $sample_a = "$root/run/A/";
my $sample_b = "$root/run/B/";
my %map = (
	opt => { smpl_order => ['A', 'B'] },
	A => { SmplID => 'A', wrdir => $sample_a, AssGroup => 'gut', SupportReads => 'PB:/reads/a.fq.gz' },
	B => { SmplID => 'B', wrdir => $sample_b, AssGroup => 'gut', SupportReads => '' },
);
my %groups = (gut => { CntAimAss => 2 });
my $group_dir = "$root/run/AssmblGrp_gut/metag";

write_file("$group_dir/scaffolds.fasta.filt", ">contig\nACGT\n");
write_file("$group_dir/ass.done.sto", "");
write_file("$group_dir/smpls_used.txt", "$sample_a\n$sample_b\n");
write_file("$sample_a/mapping/A-smd.cram.sto", "");
write_file("$sample_a/mapping/A-smd.bam.coverage.gz", "coverage\n");
write_file("$sample_a/mapping/A-smd.bam.breakpoints.tsv.gz", "breakpoints\n");

my $before = filesystem_snapshot($root);
my $report = inspect_workflow_state(map => \%map, groups => \%groups);
my $after = filesystem_snapshot($root);

is_deeply($after, $before, 'inspection does not modify the filesystem');
is($report->{read_only}, 1, 'report identifies the read-only mode');
is($report->{summary}{samples}, 2, 'reports every sample');
is($report->{summary}{assembly_groups}, 1, 'reports the assembly group');
is($report->{assembly_groups}[0]{membership}{matches}, 1, 'compares exact group membership');
is($report->{samples}[0]{stages}{assembly}{status}, 'COMPLETE', 'valid assembly is complete');
is($report->{samples}[0]{stages}{mapping}{status}, 'PARTIAL', 'mapping stone without CRAM is partial');
ok(grep($_ eq 'MARKER_WITHOUT_VALID_OUTPUT', @{$report->{samples}[0]{stages}{mapping}{issues}}),
	'reports a marker without its output');
like(encode_state_report($report), qr/"schema_version"\s*:\s*1/, 'encodes a versioned JSON report');

write_file("$group_dir/smpls_used.txt", "$sample_a\n$root/run/C/\n");
my $mismatch = inspect_workflow_state(map => \%map, groups => \%groups);
is($mismatch->{assembly_groups}[0]{membership}{matches}, 0,
	'detects changed membership even when the member count is unchanged');
ok(grep($_ eq 'GROUP_MEMBERSHIP_MISMATCH', @{$mismatch->{assembly_groups}[0]{issues}}),
	'reports group membership mismatch');

my %ont_map = (
	opt => { smpl_order => ['ONT_SAMPLE'] },
	ONT_SAMPLE => {
		SmplID => 'ONT_SAMPLE', wrdir => "$root/run/ONT_SAMPLE/", AssGroup => 'ont_hybrid',
		SupportReads => 'ONT:/reads/ont.fastq.gz', hasPrimaryRds => 1,
	},
);
my %ont_groups = (ont_hybrid => { CntAimAss => 1 });
my $ont_report = inspect_workflow_state(
	map => \%ont_map,
	groups => \%ont_groups,
	options => { assembly_mode => 5, run_tmp_dir => "$root/tmp" },
);
is($ont_report->{assembly_groups}[0]{hybrid}, 1,
	'ONT support reads use the same hybrid state logic as PB support reads');

done_testing;
