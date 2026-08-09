use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use IPC::Open3 qw(open3);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Symbol qw(gensym);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use lib File::Spec->catdir($Bin, 'lib');
use MFTestConfig;
use Mods::Binning qw(
	createBin2 createBinCtgs filterMGS_CM MB2assigns MB2assignedBinIds readCMquals
);
use Mods::geneCat qw(createGene2MGS);

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

my $tmp = tempdir(CLEANUP => 1);

my $sorted_mgs = File::Spec->catfile($tmp, 'ranked.srt');
write_file($sorted_mgs, "MGS1\t2,1\n");
{
	no warnings 'redefine';
	local *Mods::geneCat::readGene2Func = sub { return { 1 => 'COG1', 2 => 'COG2' } };
	my $mapping = createGene2MGS($sorted_mgs, $tmp);
	is(slurp($mapping), "2\tMGS1\tCOG2\n1\tMGS1\tCOG1\n",
		'sorted comma-separated MGS genes retain their priority order');
}

my $mgs = File::Spec->catfile($tmp, 'clusters.txt');
my $fasta = File::Spec->catfile($tmp, 'genes.fna');
my $bins = File::Spec->catdir($tmp, 'bins');
write_file($mgs, "MGS1\t1\nMGS2\t2\n");
write_file($fasta, ">1\nAAAA\n>2\nCCCC"); # deliberately no trailing newline
createBin2($bins, $mgs, $fasta);
like(slurp(File::Spec->catfile($bins, 'MGS1.fna')), qr/>1\nAAAA\n/, 'first FASTA record is retained');
like(slurp(File::Spec->catfile($bins, 'MGS2.fna')), qr/>2\nCCCC\n/, 'final FASTA record is flushed at EOF');

my $many_mgs = File::Spec->catfile($tmp, 'many.clusters.txt');
my $many_fasta = File::Spec->catfile($tmp, 'many.genes.fna');
my $many_bins = File::Spec->catdir($tmp, 'many-bins');
write_file($many_mgs,
	join('', map { "MGS$_\t$_\n" } 1 .. 70) . "MGS1\t71\n");
write_file($many_fasta,
	join('', map { ">$_\nSEQ$_\n" } 1 .. 71));
createBin2($many_bins, $many_mgs, $many_fasta);
is(slurp(File::Spec->catfile($many_bins, 'MGS1.fna')),
	">1\nSEQ1\n>71\nSEQ71\n",
	'bounded FASTA output handles reopen in append mode without losing records');
is_deeply(
	[glob(File::Spec->catfile($many_bins, '*.tmp.*'))], [],
	'streamed FASTA extraction publishes outputs without leftover temporary files',
);

my $cm2 = File::Spec->catfile($tmp, 'bins.cm2');
write_file($cm2, "Name\tCompleteness\tContamination\n1\t95.5\t1.2\n");
my $quality = readCMquals($cm2);
is($quality->{1}{compl}, '95.5', 'CheckM2 completeness is parsed');
is($quality->{1}{line}, "95.5\t1.2", 'quality row is retained for MGS replacement output');

my $boundary_cm2 = File::Spec->catfile($tmp, 'boundary.cm2');
write_file($boundary_cm2, "Name\tCompleteness\tContamination\nexact\t50\t5\nbelow\t49.9\t1\n");
my $passing = filterMGS_CM($boundary_cm2, 50, 5, 1);
ok(exists($passing->{exact}), 'quality values exactly on accepted thresholds are retained');
ok(!exists($passing->{below}), 'quality values below the completeness threshold are rejected');

my $assignments = File::Spec->catfile($tmp, 'assignments.tsv');
write_file($assignments, "contig1\t2\n");
eval { MB2assigns($assignments, $cm2) };
like($@, qr/No quality record for assigned bin '2'/, 'assigned bins without quality records fail explicitly');

my $empty_assignments = File::Spec->catfile($tmp, 'empty-assignments.tsv');
my $empty_quality = File::Spec->catfile($tmp, 'empty-quality.cm2');
write_file($empty_assignments, "Sequence ID\tBin\n");
write_file($empty_quality, "Name\tCompleteness\tContamination\n");
my ($empty_bins, $empty_bin_quality) = MB2assigns($empty_assignments, $empty_quality);
is_deeply($empty_bins, {}, 'a header-only bin assignment is a valid empty biological result');
is_deeply($empty_bin_quality, {}, 'an empty bin assignment does not require fabricated quality rows');

my $id_assignments = File::Spec->catfile($tmp, 'id-assignments.tsv');
write_file($id_assignments,
	"Sequence ID\tBin\ncontig1\t1\ncontig2\t1\nunassigned\t0\n");
my ($assigned_ids, $assigned_quality) = MB2assignedBinIds($id_assignments, $cm2);
is_deeply($assigned_ids, { 1 => 1 },
	'ID-only bin parsing retains unique assigned bin IDs without contig arrays');
is($assigned_quality->{1}{compl}, '95.5',
	'ID-only bin parsing retains validated quality records');

my $sample_dir = File::Spec->catdir($tmp, 'sample');
my $assembly_dir = File::Spec->catdir($tmp, 'assembly');
my $assembly_pointer_dir = File::Spec->catdir($sample_dir, 'assemblies', 'metag');
my $assembly_bin_dir = File::Spec->catdir($assembly_dir, 'Binning', 'SB');
my $contig_output_dir = File::Spec->catdir($tmp, 'representative-contigs');
make_path($assembly_pointer_dir, $assembly_bin_dir, $contig_output_dir);
write_file(File::Spec->catfile($assembly_pointer_dir, 'assembly.txt'), "$assembly_dir\n");
write_file(File::Spec->catfile($assembly_bin_dir, 'S1'),
	"needed1\t1.fa.gz\nunused\t2\nneeded2\t1.fa.gz\n");
write_file(File::Spec->catfile($assembly_dir, 'scaffolds.fasta.filt'),
	">needed1\nAAAA\n>unused\nNNNN\n>needed2\nCCCC\n");
my $representative_guide = File::Spec->catfile($tmp, 'MAGvsGC.txt');
write_file($representative_guide,
	"MAG\tMGS\tRepresentative4MGS\tCompleteness\tContamination\tLCAcompleteness\tN50\n"
	. "S1__1.fa.gz\tMGS.1\t*\t95\t1\t1\t1000\n");
my %representative_map = (
	opt => { smpl_order => ['S1'] },
	altNms => {},
	S1 => { wrdir => $sample_dir },
);
createBinCtgs(
	$contig_output_dir, \%representative_map, $representative_guide, 0, 'SB',
);
my $compressed_contigs =
	File::Spec->catfile($contig_output_dir, 'MGS.1.ctgs.S1__1.fa.gz');
my $contig_contents = '';
ok(gunzip($compressed_contigs => \$contig_contents),
	'representative contigs are published as a valid gzip stream')
	or diag($GunzipError);
is(
	$contig_contents,
	">needed1\nAAAA\n>needed2\nCCCC\n",
	'representative-contig extraction preserves the MAG FASTA suffix and selected contigs',
);
ok(!-e "$compressed_contigs.fna",
	'a pre-compressed MAG identifier does not receive a trailing .fna suffix');

my $gc = File::Spec->catdir($tmp, 'GC');
make_path(File::Spec->catdir($gc, 'Anno', 'Tax'));
my $tax_mgs = File::Spec->catfile($tmp, 'tax.clusters');
my $tax_rows = join('', map { "MGS1\t$_\n" } 1 .. 100) . "MGS2\t200\n";
write_file($tax_mgs, $tax_rows);
my $lineage = join(';', qw(Domain Phylum Class Order Family Genus Species Strain));
my $kraken = join('', map { "$_\t$lineage\n" } 1 .. 100);
write_file(File::Spec->catfile($gc, 'Anno', 'Tax', 'krak2.txt'), $kraken);
my $tax_prefix = File::Spec->catfile($tmp, 'taxonomy');
my $tax_script = File::Spec->catfile($Bin, '..', 'secScripts', 'MGS', 'taxPerMGS.pl');
is(system($^X, '-I'.File::Spec->catdir($Bin, '..'), $tax_script, $tax_mgs, $gc, $tax_prefix), 0,
	'taxPerMGS completes on classified and unclassified MGS');
my $lca = slurp("$tax_prefix.LCA");
like($lca, qr/^MGS1\tDomain;Phylum;Class;Order;Family;Genus;Species;Strain;$/m,
	'well-supported taxonomy contains all eight ranks');
like($lca, qr/^MGS2\t\?;\?;\?;\?;\?;\?;\?;\?;$/m,
	'MGS without Kraken hits is retained with eight unknown ranks');

my $missing_kraken_gc = File::Spec->catdir($tmp, 'missing-kraken-GC');
make_path(File::Spec->catdir($missing_kraken_gc, 'Anno', 'Tax'));
my $missing_kraken_prefix = File::Spec->catfile($tmp, 'missing-kraken-taxonomy');
my $missing_kraken_stdout = gensym;
my $missing_kraken_stderr = gensym;
my $missing_kraken_pid = open3(undef, $missing_kraken_stdout, $missing_kraken_stderr,
	$^X, '-I'.File::Spec->catdir($Bin, '..'),
	$tax_script, $tax_mgs, $missing_kraken_gc, $missing_kraken_prefix);
my $missing_kraken_output = do { local $/; <$missing_kraken_stdout> // '' };
my $missing_kraken_error = do { local $/; <$missing_kraken_stderr> // '' };
waitpid($missing_kraken_pid, 0);
my $missing_kraken_status = $? >> 8;
is($missing_kraken_status, 0, 'taxPerMGS treats missing optional Kraken input as a successful skip');
is($missing_kraken_output, '', 'missing optional Kraken input produces no normal output');
like($missing_kraken_error, qr/Optional Kraken input is missing or empty; skipping MGS Kraken taxonomy/,
	'taxPerMGS explains why optional Kraken taxonomy was skipped');
ok(!-e "$missing_kraken_prefix.LCA" && !-e "$missing_kraken_prefix.tax",
	'taxPerMGS does not manufacture taxonomy outputs without Kraken input');

done_testing();
