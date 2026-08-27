use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::geneCat qw(readGene2tax);
use Mods::MGSLocus qw(build_locus_groups choose_locus_candidate member_context_map robust_depth_mask);

sub write_file {
	my ($path, $contents) = @_;
	open my $fh, '>', $path or die "Cannot write $path: $!";
	print {$fh} $contents;
	close $fh or die "Cannot close $path: $!";
}

my $tmp = tempdir(CLEANUP => 1);
my $mapping = File::Spec->catfile($tmp, 'genes.gene2MGS');
write_file($mapping, "10\tMGS.1\tCOG1\n11\tMGS.1\tCOG1\n12\tMGS.1\tCOG2\n");
my ($si, $gene_to_cog, $gene_to_mgs, $priorities) = readGene2tax($mapping, 10, []);
is_deeply(
	$priorities->{'MGS.1'},
	['MGS.1|COG1|10', 'MGS.1|COG1|11', 'MGS.1|COG2|12'],
	'same-COG catalogue genes are retained in ranked order as distinct seed loci',
);
is($si->{'MGS.1'}{'MGS.1|COG1|11'}, 11, 'seed locus maps to its catalogue cluster');

my $filtered_mapping = File::Spec->catfile($tmp, 'filtered-genes.gene2MGS');
write_file($filtered_mapping, join('',
	"30\tMGS.3\tCOG_OTHER\n",
	"31\tMGS.3\tCOG_WANTED\n",
	"32\tMGS.3\tCOG_WANTED\n",
));
my (undef, undef, undef, $filtered_priorities) = readGene2tax(
	$filtered_mapping, 1, ['MGS.3'], undef,
	{ allowed_cogs_by_mgs => { 'MGS.3' => { COG_WANTED => 1 } } },
);
is_deeply($filtered_priorities->{'MGS.3'}, ['MGS.3|COG_WANTED|31'],
	'a COG demand is applied before the per-MGS gene cap');

my $duplicates = File::Spec->catfile($tmp, 'duplicate-genes.gene2MGS');
write_file($duplicates, join('',
	"20\tMGS.2\tCOG1\n",
	map { "20\tMGS.2\tCOG$_\n" } 2 .. 10,
));
my $duplicate_warnings = '';
{
	local $SIG{__WARN__} = sub { $duplicate_warnings .= join('', @_); };
	readGene2tax($duplicates, 20, []);
}
my @duplicate_examples = $duplicate_warnings =~ /Ignoring duplicate catalogue gene 20/g;
is(scalar(@duplicate_examples), 5,
	'duplicate catalogue gene diagnostics are limited to five detailed examples');
like($duplicate_warnings, qr/Suppressed 4 additional duplicate catalogue gene warnings.*9 total/,
	'duplicate catalogue gene diagnostics retain an aggregate suppressed count');

my @records = (
	{ mgs => 'MGS.1', cog => 'COG1', gene => '10', rank => 0 },
	{ mgs => 'MGS.1', cog => 'COG1', gene => '11', rank => 1 },
	{ mgs => 'MGS.1', cog => 'COG1', gene => '12', rank => 2 },
);
my %members = (
	10 => 'S1__ctg_1,S2__ctg_1',
	11 => 'S3__ctg_1,S4__ctg_1',
	12 => 'S1__ctg_5,S3__ctg_5',
);
my $protein = 'M' . ('ACDEFGHIKLMNPQRSTVWY' x 8);
my %proteins = (10 => $protein, 11 => $protein, 12 => $protein);
my $model = build_locus_groups(\@records, \%members, \%proteins);
is(scalar(@{$model->{groups}}), 2, 'mutually exclusive similar seeds merge but a co-occurring paralog stays separate');
is($model->{gene_to_locus}{10}, 'MGS.1|COG1|10', 'first ranked seed names the merged locus');
is($model->{gene_to_locus}{11}, 'MGS.1|COG1|10', 'alternative seed maps to the primary locus');
is($model->{gene_to_locus}{12}, 'MGS.1|COG1|12', 'co-occurring seed retains its own locus');
is($model->{member_to_seed}{'S1__ctg_1'}, '10',
	'default locus-model results retain the member-to-seed compatibility index');

my $lean_model = build_locus_groups(
	\@records, \%members, \%proteins,
	{ include_member_to_seed => 0, include_gene_to_locus => 0 },
);
is_deeply($lean_model->{member_to_seed}, {},
	'callers can omit the memory-heavy member-to-seed result');
is_deeply($lean_model->{gene_to_locus}, {},
	'callers can omit an unused gene-to-locus result');
is_deeply(
	[map { $_->{locus_id} } @{$lean_model->{groups}}],
	[map { $_->{locus_id} } @{$model->{groups}}],
	'omitting optional indexes does not change locus construction',
);

my %allow_only_10_11 = ("10\t11" => 1);
my $catalogue_gated = build_locus_groups(
	\@records, \%members, \%proteins,
	{allowed_merge_pairs => \%allow_only_10_11, require_complete_linkage => 1},
);
is(scalar(@{$catalogue_gated->{groups}}), 2,
	'catalogue-confirmed pair can merge while unconfirmed same-COG seeds stay separate');
my $mosaic_disabled = build_locus_groups(
	\@records, \%members, \%proteins,
	{allowed_merge_pairs => {}, require_complete_linkage => 1},
);
is(scalar(@{$mosaic_disabled->{groups}}), 3,
	'an empty mosaic allowlist safely retains every same-COG seed as a separate locus');
is($mosaic_disabled->{merged_seeds}, 0,
	'disabling mosaic checks performs no seed-cluster merging');

my @mixed_strain_records = @records[0, 1];
my %mixed_strain_members = (
	10 => 'S1__ctg_1,S2__ctg_1,S3__ctg_1',
	11 => 'S3__ctg_2,S4__ctg_1,S5__ctg_1',
);
my $mixed_strain_model = build_locus_groups(
	\@mixed_strain_records, \%mixed_strain_members, \%proteins,
	{
		allowed_merge_pairs => \%allow_only_10_11,
		require_complete_linkage => 1,
		allow_confirmed_cooccurrence => 1,
	},
);
is(scalar(@{$mixed_strain_model->{groups}}), 1,
	'a catalogue-confirmed pair can survive one mixed-strain sample co-occurrence');

my @chain_records = (
	{mgs => 'MGS.3', cog => 'COG1', gene => '20', rank => 0},
	{mgs => 'MGS.3', cog => 'COG1', gene => '21', rank => 1},
	{mgs => 'MGS.3', cog => 'COG1', gene => '22', rank => 2},
);
my %chain_members = (
	20 => 'S1__ctg_1',
	21 => 'S2__ctg_1',
	22 => 'S3__ctg_1',
);
my %chain_proteins = map { $_ => $protein } qw(20 21 22);
my %chain_edges = ("20\t21" => 1, "21\t22" => 1);
my $complete_linkage = build_locus_groups(
	\@chain_records, \%chain_members, \%chain_proteins,
	{allowed_merge_pairs => \%chain_edges, require_complete_linkage => 1},
);
is(scalar(@{$complete_linkage->{groups}}), 2,
	'complete linkage prevents transitive mosaic chaining without a 20-22 confirmation');
is($complete_linkage->{incomplete_linkage_rejections}, 1,
	'complete-linkage protection reports the rejected transitive merge');
my $transitive_mosaic = build_locus_groups(
	\@chain_records, \%chain_members, \%chain_proteins,
	{allowed_merge_pairs => \%chain_edges, require_complete_linkage => 0},
);
is(scalar(@{$transitive_mosaic->{groups}}), 1,
	'transitive Mosaic edges can represent three alternatives of one homologue');
is_deeply($transitive_mosaic->{groups}[0]{genes}, [qw(20 21 22)],
	'all three genes are retained in the merged Mosaic locus');

my $dominant = choose_locus_candidate([
	{ id => 'copyA', protein => $protein, depth => 12, seed => '10', context => {} },
	{ id => 'copyB', protein => $protein, depth => 3,  seed => '10', context => {} },
], { 10 => $protein }, {});
is($dominant->{status}, 'selected', 'a dominant duplicate can be resolved');
is($dominant->{candidate}{id}, 'copyA', 'dominant-depth candidate is selected within one seed');

my $ambiguous = choose_locus_candidate([
	{ id => 'copyA', protein => $protein, depth => 10, seed => '10', context => {} },
	{ id => 'copyB', protein => $protein, depth => 10, seed => '12', context => {} },
], { 10 => $protein, 12 => $protein }, {});
is($ambiguous->{status}, 'ambiguous', 'equivalent candidates from different seeds remain unresolved');

is_deeply(
	robust_depth_mask([10, 10, 11, 9, 10, 11, 9, 100]),
	[1, 1, 1, 1, 1, 1, 1, 0],
	'robust depth filtering removes a supported extreme outlier',
);
is_deeply(
	robust_depth_mask([10, 10, 100]),
	[1, 1, 1],
	'small locus sets are not quantile-filtered',
);


# --- Phase-I split-worker locus-model determinism -----------------------------
# Same-COG seeds merge using catalogue-wide co-occurrence and synteny context.
# A split worker only holds the members of its own samples, so a model built
# from that shard can disagree with the catalogue-wide one about locus identity.
# strain_within.pl therefore builds the model once in the parent and publishes
# it; these tests pin both halves of that contract.
my $slice_protein = 'MKAILVVLLYTFATANADTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCKL'
	. 'RGVAPLHLGKCNIAGWILGNPECESLSTASSWSYIVETSSSDNGTCYPGDFIDYEELREQL'
	. 'SSVSSFERFEIFPKESSWPNHNTNGVTAACSHEGKSSFYRNLLWLTEKEGSYPKLKNSYVN';
my $slice_variant = $slice_protein;
substr($slice_variant, $_ * 10, 1) = substr("ACDEFGHIKLMNPQRSTVWY", $_ % 20, 1)
	for 1 .. 18;

my @slice_records = (
	{ mgs => 'MGS.9', cog => 'COG1', gene => 'g1', rank => 0 },
	{ mgs => 'MGS.9', cog => 'COG1', gene => 'g2', rank => 1 },
	{ mgs => 'MGS.9', cog => 'COG2', gene => 'n1', rank => 2 },
	{ mgs => 'MGS.9', cog => 'COG3', gene => 'n2', rank => 3 },
);
my %slice_proteins = (
	g1 => $slice_protein, g2 => $slice_variant,
	n1 => $slice_protein, n2 => $slice_variant,
);
# smplA carries g1, smplB carries g2; both keep the same two neighbours.
my %catalogue_members = (
	g1 => 'smplA__ctgA_10',
	g2 => 'smplB__ctgB_10',
	n1 => 'smplA__ctgA_11,smplB__ctgB_11',
	n2 => 'smplA__ctgA_12,smplB__ctgB_12',
);
my %worker_members = (
	g1 => 'smplA__ctgA_10',
	g2 => 'smplB__ctgB_10',
	n1 => 'smplA__ctgA_11',
	n2 => 'smplA__ctgA_12',
);
my %slice_options = (
	allowed_merge_pairs => { "g1\tg2" => 1 },
	require_complete_linkage => 1,
	allow_confirmed_cooccurrence => 1,
	include_member_to_seed => 0,
	include_gene_to_locus => 0,
);

my $catalogue_model = build_locus_groups(
	\@slice_records, \%catalogue_members, \%slice_proteins, {%slice_options});
my $worker_model = build_locus_groups(
	\@slice_records, \%worker_members, \%slice_proteins, {%slice_options});
my @catalogue_loci = sort map { $_->{locus_id} }
	grep { $_->{cog} eq 'COG1' } @{$catalogue_model->{groups}};
my @worker_loci = sort map { $_->{locus_id} }
	grep { $_->{cog} eq 'COG1' } @{$worker_model->{groups}};
is_deeply(\@catalogue_loci, ['MGS.9|COG1|g1'],
	'catalogue-wide synteny context merges a confirmed same-COG pair into one locus');
is_deeply(\@worker_loci, ['MGS.9|COG1|g1', 'MGS.9|COG1|g2'],
	'a sample-restricted shard splits one locus in two, so the parent must publish the model');

my $context_only = build_locus_groups(
	\@slice_records, \%catalogue_members, \%slice_proteins,
	{%slice_options, include_member_context => 0});
is_deeply(
	[sort map { $_->{locus_id} } @{$context_only->{groups}}],
	[sort map { $_->{locus_id} } @{$catalogue_model->{groups}}],
	'omitting member contexts does not change the published locus model',
);
is_deeply($context_only->{member_context}, {},
	'callers can skip the catalogue-wide member contexts they never consume');

# Member contexts describe one member's own contig neighbourhood, so a worker
# can rebuild exactly its own share around the published model.
my $worker_contexts = member_context_map(\@slice_records, \%worker_members);
my %expected_worker_contexts = map {
	$_ => $catalogue_model->{member_context}{$_}
} grep { /^smplA__/ } keys %{$catalogue_model->{member_context}};
is_deeply($worker_contexts, \%expected_worker_contexts,
	'member contexts derived from one worker shard match the catalogue-wide contexts');
ok(scalar(keys %{$worker_contexts}),
	'the worker shard actually produced member contexts to compare');

done_testing();
