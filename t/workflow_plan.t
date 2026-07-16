use strict;
use warnings;

use File::Find;
use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::WorkflowPlan qw(build_workflow_plan encode_workflow_plan validate_workflow_plan);
use Mods::WorkflowState qw(inspect_workflow_state);

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

sub actions_by_id {
	my ($plan) = @_;
	return map { $_->{id} => $_ } @{$plan->{actions}};
}

my $root = tempdir(CLEANUP => 1);
$root =~ s{\\}{/}g;
my $sample_a = "$root/ordinary/A/";
my %ordinary_map = (
	opt => { smpl_order => ['A'] },
	A => {
		SmplID => 'A', wrdir => $sample_a, AssGroup => 'ordinary',
		SupportReads => '', hasPrimaryRds => 1,
	},
);
my %ordinary_groups = (ordinary => { CntAimAss => 1 });
my $assembly_dir = "$sample_a/assemblies/metag";
write_file("$assembly_dir/scaffolds.fasta.filt", ">contig\nACGT\n");
write_file("$assembly_dir/ass.done.sto", "done\n");
write_file("$assembly_dir/smpls_used.txt", "$sample_a\n");
write_file("$sample_a/mapping/A-smd.cram.sto", "done\n");

my $before = filesystem_snapshot($root);
my $ordinary_state = inspect_workflow_state(
	map => \%ordinary_map,
	groups => \%ordinary_groups,
	options => {
		assembly_mode => 1,
		map_to_assembly => 1,
		map_support_to_assembly => 0,
		run_tmp_dir => "$root/tmp/run",
	},
);
my $ordinary_plan = build_workflow_plan($ordinary_state);
my $after = filesystem_snapshot($root);
is_deeply($after, $before, 'planning does not modify the filesystem');
is($ordinary_plan->{execution_supported}, 0, 'plan cannot execute actions');
is_deeply([validate_workflow_plan($ordinary_plan)], [], 'ordinary plan dependency graph is valid');
my %ordinary_actions = actions_by_id($ordinary_plan);
ok($ordinary_actions{'sample:A:repair:mapping'}, 'partial mapping gets an explicit repair action');
ok($ordinary_actions{'sample:A:submit:mapping'}, 'partial mapping gets a separate submit action');
is_deeply(
	$ordinary_actions{'sample:A:submit:mapping'}{depends_on},
	['sample:A:repair:mapping'],
	'mapping submission waits for its repair action',
);
ok($ordinary_actions{'group:ordinary:submit:gene_prediction'},
	'missing group gene prediction is planned separately');
ok(grep($_ eq 'group:ordinary:submit:gene_prediction',
	@{$ordinary_actions{'sample:A:submit:contig_stats'}{depends_on}}),
	'contig statistics wait for gene prediction');
like(encode_workflow_plan($ordinary_plan), qr/"mode"\s*:\s*"repair_submission_plan"/,
	'encodes a versioned plan document');

write_file("$assembly_dir/smpls_used.txt", "$root/ordinary/DIFFERENT/\n");
my $changed_state = inspect_workflow_state(
	map => \%ordinary_map,
	groups => \%ordinary_groups,
	options => { assembly_mode => 1, map_to_assembly => 1, map_support_to_assembly => 0 },
);
my $changed_plan = build_workflow_plan($changed_state);
my %changed_actions = actions_by_id($changed_plan);
is($changed_actions{'group:ordinary:repair:membership'}{authorization}, 'OKtoRWassGrps',
	'assembly-group invalidation requires the existing destructive authorization');
ok($changed_actions{'group:ordinary:repair:membership'}{requires_confirmation},
	'group-wide repair requires confirmation');
ok(grep($_ eq 'group:ordinary:repair:membership',
	@{$changed_actions{'group:ordinary:submit:assembly'}{depends_on}}),
	'changed group is repaired before its assembly is resubmitted');
ok($changed_actions{'sample:A:submit:mapping'},
	'group invalidation forces completed downstream mapping to be resubmitted');

my $hybrid_a = "$root/hybrid/A/";
my $hybrid_b = "$root/hybrid/B/";
my %hybrid_map = (
	opt => { smpl_order => ['A', 'B'] },
	A => {
		SmplID => 'A', wrdir => $hybrid_a, AssGroup => 'hybrid',
		SupportReads => 'PB:/reads/a.fastq.gz', hasPrimaryRds => 1,
	},
	B => {
		SmplID => 'B', wrdir => $hybrid_b, AssGroup => 'hybrid',
		SupportReads => '', hasPrimaryRds => 1,
	},
);
my %hybrid_groups = (hybrid => { CntAimAss => 2, SupportReads => ',PB:/reads/a.fastq.gz' });
my $hybrid_state = inspect_workflow_state(
	map => \%hybrid_map,
	groups => \%hybrid_groups,
	options => {
		assembly_mode => 5,
		map_to_assembly => 1,
		map_support_to_assembly => 1,
		run_tmp_dir => "$root/tmp/hybrid-run",
	},
);
my $hybrid_plan = build_workflow_plan($hybrid_state);
my %hybrid_actions = actions_by_id($hybrid_plan);
is($hybrid_state->{assembly_groups}[0]{hybrid}, 1, 'inspection identifies hybrid group logic');
ok($hybrid_actions{'sample:A:submit:preassembly'}, 'first hybrid preassembly is explicit');
ok($hybrid_actions{'sample:B:submit:preassembly'}, 'second hybrid preassembly is explicit');
ok(grep($_ eq 'sample:A:submit:preassembly',
	@{$hybrid_actions{'sample:B:submit:preassembly'}{depends_on}}),
	'hybrid preassemblies are serialized to protect the shared group staging area');
is($hybrid_actions{'group:hybrid:submit:assembly'}{operation}, 'submit_hybrid_group_assembly',
	'final hybrid assembly is a distinct group action');
ok(grep($_ eq 'sample:B:submit:preassembly',
	@{$hybrid_actions{'group:hybrid:submit:assembly'}{depends_on}}),
	'final hybrid assembly waits for all preassembly packages');
ok(grep($_ eq 'group:hybrid:submit:assembly',
	@{$hybrid_actions{'sample:A:submit:support_mapping'}{depends_on}}),
	'support mapping waits for final hybrid assembly');
is_deeply([validate_workflow_plan($hybrid_plan)], [], 'hybrid dependency graph is valid');

my $package_a = "$root/tmp/hybrid-run/A/preAssmblGrp_hybrid";
write_file("$package_a/scaffolds.fasta.filt", ">preassembly\nACGT\n");
write_file("$package_a/Coverage.percontig.gz", "coverage\n");
write_file("$package_a/Coverage.median.percontig.gz", "median\n");
write_file("$package_a/moved.sto", "done\n");
my $resumed_state = inspect_workflow_state(
	map => \%hybrid_map,
	groups => \%hybrid_groups,
	options => {
		assembly_mode => 5,
		map_to_assembly => 1,
		map_support_to_assembly => 1,
		run_tmp_dir => "$root/tmp/hybrid-run",
	},
);
my $resumed_plan = build_workflow_plan($resumed_state);
my %resumed_actions = actions_by_id($resumed_plan);
is($resumed_state->{samples}[0]{stages}{preassembly_package}{status}, 'COMPLETE',
	'completed hybrid preassembly package is detected');
ok(!$resumed_actions{'sample:A:submit:preassembly'},
	'completed hybrid preassembly is not resubmitted');
ok($resumed_actions{'sample:B:submit:preassembly'},
	'remaining hybrid preassembly is still planned');

done_testing;
