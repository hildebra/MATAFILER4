use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::Binning qw(createBin2 filterMGS_CM MB2assigns readCMquals);

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
my $mgs = File::Spec->catfile($tmp, 'clusters.txt');
my $fasta = File::Spec->catfile($tmp, 'genes.fna');
my $bins = File::Spec->catdir($tmp, 'bins');
write_file($mgs, "MGS1\t1\nMGS2\t2\n");
write_file($fasta, ">1\nAAAA\n>2\nCCCC"); # deliberately no trailing newline
createBin2($bins, $mgs, $fasta);
like(slurp(File::Spec->catfile($bins, 'MGS1.fna')), qr/>1\nAAAA\n/, 'first FASTA record is retained');
like(slurp(File::Spec->catfile($bins, 'MGS2.fna')), qr/>2\nCCCC\n/, 'final FASTA record is flushed at EOF');

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

my $between_source = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'MGS', 'phylo_MGS_between.pl'));
like($between_source, qr/open I,"<\$GCd\/FMG\.subset\.cats"/,
	'between-MGS phylogeny intentionally remains tied to the FMG marker set');
my $strain_source = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'MGS', 'strain_within.pl'));
like($strain_source, qr/my \$tree_sample_separator = quotemeta\(\$SaSe\)/,
	'within-MGS tree command escapes the pipe sample separator as a regular expression');
my $resort_source = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'MGS', 'resortMGSgenes4importance.pl'));
like($resort_source, qr/print O evalCurMGS\(""\) if \$curMGS ne "";/,
	'gene-priority resorting flushes its final MGS at EOF');

done_testing();
