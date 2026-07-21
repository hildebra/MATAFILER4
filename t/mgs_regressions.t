use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use IPC::Open3 qw(open3);
use Symbol qw(gensym);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::Binning qw(createBin2 filterMGS_CM MB2assigns readCMquals);
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

my $between_source = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'MGS', 'phylo_MGS_between.pl'));
like($between_source, qr/open I,"<\$GCd\/FMG\.subset\.cats"/,
	'between-MGS phylogeny intentionally remains tied to the FMG marker set');
my $strain_source = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'MGS', 'strain_within.pl'));
like($strain_source, qr/my \$tree_sample_separator = quotemeta\(\$SaSe\)/,
	'within-MGS tree command escapes the pipe sample separator as a regular expression');
like($strain_source, qr/\$nxtCmd \.= "-Hcores \$maxCores " if \$maxCores > 0;/,
	'within-MGS analysis only forwards a configured positive heavy-core limit');
my $strain2_source = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'MGS', 'strain_within_2.2.pl'));
like($strain2_source, qr/while \(!-s \$treePath && \$x < \@defTreeFiles\)/,
	'strain postprocessing selects only nonempty fallback trees');
like($strain2_source, qr/my \$jobCores = \$nCore;/,
	'strain postprocessing honours its cores option');
like($strain2_source, qr/qsubSystemWaitMaxJobs\(\$checkMaxNumJobs,0,\$QSBoptHR\) if \$doSubmit;/,
	'queue throttling uses the selected backend and is disabled for dry runs');
like($strain2_source, qr/"test -s "\.shellQuote\(\$analysisReport\).*?"touch "\.shellQuote\(\$analysisStone\)/s,
	'new R analyses validate their report before writing a success stone');
like($strain2_source, qr/my \$MGSd = dirname\(\$FMGpD\);/,
	'treeWAS receives a trailing-slash-independent parent directory');
like($strain2_source, qr/split \/\\Q\$SaSe\\E\/, \$gn, 2/,
	'strain postprocessing splits compound tree identifiers at the first separator');
my $resort_source = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'MGS', 'resortMGSgenes4importance.pl'));
like($resort_source, qr/print O evalCurMGS\(""\) if \$curMGS ne "";/,
	'gene-priority resorting flushes its final MGS at EOF');
like($resort_source, qr/compl\.incompl\.\$clusterID\.fna\.clstr\.idx/,
	'gene-priority resorting uses the propagated catalog identity');
my $mgs_source = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'MGS.pl'));
like($mgs_source, qr/"clusterID=i" => \\\$clusterID/,
	'MGS accepts a gene-catalog cluster identity');
like($mgs_source, qr/-MGset \$useGTDBmg -clusterID \$clusterID -cores \$numCore/,
	'MGS passes cluster identity to MAG clustering');
like($mgs_source, qr/-MGset \$useGTDBmg -clusterID \$clusterID -maxCores \$canCore/,
	'MGS passes cluster identity to strain analysis');
unlike($mgs_source, qr/die "Kraken MGS taxonomy stage incomplete/,
	'missing optional Kraken taxonomy no longer aborts MGS');
like($mgs_source, qr/warn "Optional Kraken MGS taxonomy stage incomplete; continuing without Kraken-derived MGS taxonomy/,
	'missing optional Kraken taxonomy is still reported');
like($mgs_source, qr/if \(!-s \$krakenInput\).*?skipping MGS Kraken taxonomy:.*?else \{.*?qsubSystem\(\$logDir\."\/krak2MGS\.sh"/s,
	'MGS does not submit the optional Kraken taxonomy job without its input');
unlike($mgs_source, qr/compl\.incompl\.95\.(?:fna|prot)/,
	'MGS has no active catalog path pinned to identity 95');
like($strain_source, qr/compl\.incompl\.\$clusterID\.fna\.clstr\.idx/,
	'within-MGS analysis reads the selected catalog index');
like($strain_source,
	qr/\$MGSfile = \$sortedMGS;\s+\$gene2taxF = createGene2MGS\(\$MGSfile,\$GCd\)/s,
	'within-MGS analysis builds its gene mapping from the sorted guide');
like($strain_source, qr/'-forceSNPcalls', \$forceVCF2FNA/,
	'parallel extraction workers inherit forced consensus regeneration');
like($strain_source, qr/'-SNPadaptiveQual', \$useAdaptiveQual/,
	'parallel extraction workers inherit adaptive SNP filtering');
like($strain_source, qr/rename \$mergeFile, \$outfile or die/,
	'part-file merging publishes completed output atomically');
like($strain_source,
	qr/!\$reSubmit && !\$repairCAT && !\$redoSubmissionData.*?&& -e \$treeStone/s,
	'explicit repair and resubmission modes bypass completed-tree skipping');
like($strain_source, qr/tooFewSamples\.sto/,
	'undersampled MGS are checkpointed separately from missing inputs');
like($strain_source, qr/falling back to on-the-fly generation/,
	'failed consensus precomputation has a local fallback');
like($strain_source, qr/my \$ng = "\$sd3\$SaSe" \. externalLocusName\(\$curLocus\{\$gX\}, \$MGS\)/,
	'within-MGS sequence identifiers carry sample, COG, and primary catalogue gene');
like($strain_source, qr/internalLocusName\(\$external_locus, \$MGS\)/,
	'external category loci are restored to exact MGS-qualified internal keys');
like($strain_source, qr/\@identifier_parts != 3/,
	'stale two- or four-part sequence identifiers trigger input regeneration');
like($strain_source, qr/!exists\(\$legacyLocusMGS\{\$MGS\}\).*?-e \$treeStone/s,
	'stale identifier formats cannot be hidden by an existing tree checkpoint');
like($strain_source, qr/stale sequence identifiers but no regenerated temporary input/,
	'stale final inputs are not silently reused when regeneration produced no files');
like($strain_source, qr/No prior conspecific-sample log found.*?continuing without historical exclusions/,
	'missing optional exclusion logs do not abort sparse-MGS resume runs');
like($strain_source, qr/robust_depth_mask\(\\\@abunGs\)/,
	'within-MGS abundance filtering uses a robust depth mask');

my $build_tree_source = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'phylo', 'buildTree5.pl'));
like($build_tree_source, qr/\(\?<sample>\.\*\?\).*?\(\?<gene>\.\+\)/,
	'tree sequence identifiers split at the first separator and retain compound locus names');
like($build_tree_source, qr/sub geneFileStem.*?sprintf\("_%02X", ord\(\$1\)\)/s,
	'compound locus names are encoded safely and deterministically for downstream filenames');

done_testing();
