use strict;
use warnings;

use File::Copy qw(copy);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

my $root = File::Spec->catdir($Bin, q{..});
my $helper = File::Spec->catfile(
	$root, qw(secScripts MGS finalize_strain_tree_inputs.pl));
ok(-s $helper, 'standalone strain shard finalizer exists');

my $staging = tempdir(CLEANUP => 1);
my $published = tempdir(CLEANUP => 1);
my $fake_pigz = File::Spec->catfile($staging, 'fake-pigz');
write_file($fake_pigz, "#!/bin/sh\ncat\n");
chmod 0755, $fake_pigz or die "Cannot mark $fake_pigz executable: $!";

my %part = (
	'allFNAs.fna.0' => ">sA|COG2|g2\nATG2\n>sA|COG1|g1\nATG1\n",
	'allFAAs.faa.0' => ">sA|COG2|g2\nP2\n>sA|COG1|g1\nP1\n",
	'link2GC.txt.0' => "sA|COG2|g2\tgene2\nsA|COG1|g1\tgene1\n",
	'all.cat.tmp.0' => "MGS.2\tMGS.2|COG2|g2\tsA\tsA|COG2|g2\n"
		. "MGS.2\tMGS.2|COG1|g1\tsA\tsA|COG1|g1\n",
	'sampleQC.tsv.tmp.0' => "MGS.2\tsA\tbackbone\t0.1\t0.2\t2\n",
	'allFNAs.fna.1' => ">sB|COG1|g1\nBTG1\n",
	'allFAAs.faa.1' => ">sB|COG1|g1\nBP1\n",
	'link2GC.txt.1' => "sB|COG1|g1\tgene1\n",
	'all.cat.tmp.1' => "MGS.2\tMGS.2|COG1|g1\tsB\tsB|COG1|g1\n",
	'sampleQC.tsv.tmp.1' => "MGS.2\tsB\tplacement\t0.3\t0.1\t1\n",
);
write_file(File::Spec->catfile($staging, $_), $part{$_}) for keys %part;
write_file(File::Spec->catfile($staging, '.strain_tree_input.outgroup.fna'),
	">MGS.9|COG1|g1\nOTG1\n");
write_file(File::Spec->catfile($staging, '.strain_tree_input.outgroup.faa'),
	">MGS.9|COG1|g1\nOP1\n");
write_file(File::Spec->catfile($staging, '.strain_tree_input.outgroup.cat.tsv'),
	"MGS.2|COG1|g1\tMGS.9\tMGS.9|COG1|g1\n");
write_file(File::Spec->catfile($staging, '.strain_tree_input.plan.tsv'),
	"strain-staged-input-v1\noutgroup\tMGS.9\nmgs\tMGS.2\n");

my $manifest = File::Spec->catfile($staging, '.strain_tree_input.shards.tsv');
my @manifest = (
	'strain-shard-input-v1',
	"value\tmgs\tMGS.2", "value\toutgroup\tMGS.9",
	"value\tgeneration\ttest.1", "value\tseparator\t|",
	"value\texpected_records\t3", "value\texpected_ingroup_samples\t2",
	"value\texpected_loci\t2",
	"output\tfna\tallFNAs.fna", "output\tfaa\tallFAAs.faa",
	"output\tlink\tlink2GC.txt", "output\tcategory\tall.cat",
	"output\tqc\tsampleQC.tsv", "output\tdata_log\tdata.log",
);
for my $worker (0, 1) {
	push @manifest, join("\t", 'worker', $worker, 1, $worker ? 1 : 2);
	for my $entry (
		['fna', 'allFNAs.fna'], ['faa', 'allFAAs.faa'],
		['link', 'link2GC.txt'], ['category', 'all.cat.tmp'],
		['qc', 'sampleQC.tsv.tmp'],
	) {
		my $name = "$entry->[1].$worker";
		push @manifest, join("\t", 'part', $worker, $entry->[0],
			$name, -s File::Spec->catfile($staging, $name));
	}
}
write_file($manifest, join("\n", @manifest)."\n");

my @prepare = ($^X, $helper, '-staging', $staging, '-manifest', $manifest,
	'-mode', 'prepare', '-pigz', $fake_pigz, '-cores', 2);
is(system(@prepare), 0, 'standalone helper fuses a valid two-worker shard handoff');
ok(-s File::Spec->catfile($staging, '.strain_tree_input.prepared.tsv'),
	'helper publishes its prepared marker last');
ok(-s File::Spec->catfile($staging, 'allFNAs.fna.0')
		&& -s File::Spec->catfile($staging, 'allFNAs.fna.1'),
	'worker shards remain until publication is committed');

is(read_file(File::Spec->catfile($staging, 'allFNAs.fna.gz')),
	">MGS.9|COG1|g1\nOTG1\n>sA|COG1|g1\nATG1\n"
	. ">sB|COG1|g1\nBTG1\n>sA|COG2|g2\nATG2\n",
	'FNA shards and outgroup overlay are fused directly in locus/sample order');
is(read_file(File::Spec->catfile($staging, 'all.cat.gz')),
	"MGS.9|COG1|g1\tsA|COG1|g1\tsB|COG1|g1\n"
	. "sA|COG2|g2\n",
	'category shards are grouped directly with the outgroup overlay');
is(read_file(File::Spec->catfile($staging, 'sampleQC.tsv.gz')),
	"MGS\tsample\tstatus\tambiguous_fraction\tcsp_fraction\tvalidated_loci\n"
	. "MGS.2\tsA\tbackbone\t0.1\t0.2\t2\n"
	. "MGS.2\tsB\tplacement\t0.3\t0.1\t1\n",
	'QC shards are finalized deterministically');

my $changed = File::Spec->catfile($staging, 'allFNAs.fna.0');
write_file($changed, $part{'allFNAs.fna.0'}."\n");
my $failure = qx{@prepare 2>&1};
ok($? != 0, 'helper rejects a shard changed after manifest publication');
like($failure, qr/Shard part size changed/, 'changed-shard failure is explicit');
write_file($changed, $part{'allFNAs.fna.0'});

my @cleanup = ($^X, $helper, '-staging', $staging, '-manifest', $manifest,
	'-mode', 'cleanup', '-publishedDir', $published);
my $cleanup_failure = qx{@cleanup 2>&1};
my $cleanup_status = $?;
ok($cleanup_status != 0,
	'cleanup refuses to discard shards before every output is published');
like($cleanup_failure, qr/Cannot clean shard handoff before every finalized output is published/,
	'cleanup refusal identifies the incomplete publication contract');
ok(-s $manifest && -s File::Spec->catfile($staging, 'allFNAs.fna.0'),
	'refused cleanup preserves the manifest and worker shards');

for my $name (qw(allFNAs.fna.gz allFAAs.faa.gz all.cat.gz sampleQC.tsv.gz
	link2GC.txt.gz data.log.gz)) {
	copy(File::Spec->catfile($staging, $name), File::Spec->catfile($published, $name))
		or die "Cannot publish test fixture $name: $!";
}
write_file(
	File::Spec->catfile($staging, '.strain_tree_input.cleanup.tsv'),
	"strain-shard-cleanup-v1\ttest.1\tMGS.2\tMGS.9\n");
my $already_cleaned = File::Spec->catfile($staging, 'allFAAs.faa.0');
unlink $already_cleaned
	or die "Cannot simulate interrupted cleanup for $already_cleaned: $!";
my @resume_prepare = (@prepare, '-publishedDir', $published);
is(system(@resume_prepare), 0,
	'prepare mode finishes interrupted cleanup without rebuilding published artifacts');
ok(!-e $manifest && !-e File::Spec->catfile($staging, 'allFNAs.fna.0')
		&& !-e File::Spec->catfile($staging, 'sampleQC.tsv.tmp.1'),
	'committed cleanup removes the manifest and its exact worker shards');

sub write_file {
	my ($path, $contents) = @_;
	open my $output, '>', $path or die "Cannot write $path: $!";
	print {$output} $contents or die "Cannot write contents to $path: $!";
	close $output or die "Cannot close $path: $!";
}

sub read_file {
	my ($path) = @_;
	open my $input, '<', $path or die "Cannot read $path: $!";
	local $/;
	my $contents = <$input>;
	close $input or die "Cannot close $path: $!";
	return $contents;
}

done_testing();
