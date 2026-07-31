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

my $between_source = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'MGS', 'phylo_MGS_between.pl'));
like($between_source,
	qr/"MGset=s" => \\\$markerSet.*?my \$markerTag = \$markerSet eq "GTDB" \? "GTDBmg" : "FMG".*?open I,"<\$markerSubset"/s,
	'between-MGS phylogeny selects predefined GTDB or FMG catalog markers');
like($between_source,
	qr/my \$first = delete \$MGSFMG.*?\$ambiguousMGSFMG.*?Retained \$usable_fmg unambiguous \$markerSet marker genes/s,
	'between-MGS phylogeny excludes ambiguous paralogous marker cells instead of choosing the first copy');
like($between_source, qr/foreach my \$mg \(sort keys %MGSFMG\).*?foreach my \$cg \(sort keys %catT\)/s,
	'between-MGS FASTA and category records are emitted deterministically');
like($between_source, qr/if \(\$mgs_with_fmg < 3\).*?SKIPPED=too_few_marker_bearing_MGS/s,
	'between-MGS phylogeny reports a successful cardinality skip before invoking a tree builder');
like($between_source, qr/if \(!-s \$treeFile\)/,
	'between-MGS phylogeny rebuilds an empty tree instead of treating it as complete');
like($between_source, qr/test -s \$treeFile.*?if \(\$visualize\).*?test -s \$treePdf/s,
	'standalone between-MGS jobs validate the tree and conditionally validate visualization');
like($between_source,
	qr/-bootstrap 1000.*?-NTfilt 0\.5.*?-NTfiltPerGene 0\.7.*?-GenesPerSpecies 0\.5.*?-MSAprogram \$MSAprog.*?-AutoModel 1/s,
	'between-MGS trees use broad-phylogeny filtering and partition-aware model selection');
unlike($between_source, qr/minimumColumnOccupancy|-minOverlapMSA/,
	'between-MGS trees do not impose a fixed taxon-count overlap filter on alignment columns');
like($between_source, qr/\$externalDep =~ s\/\^\\Q\$localTag\\E\/\/.*?WAITID=\$externalDep/s,
	'between-MGS launcher exports an untagged numeric dependency across processes');
my $strain_source = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'MGS', 'strain_within.pl'));
like($strain_source, qr/my \$tree_sample_separator = quotemeta\(\$SaSe\)/,
	'within-MGS tree command escapes the pipe sample separator as a regular expression');
like($strain_source, qr/\$nxtCmd \.= "-Hcores \$maxCores " if \$maxCores > 0;/,
	'within-MGS analysis only forwards a configured positive heavy-core limit');
my $strain2_source = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'MGS', 'strain_within_2.2.pl'));
like($strain2_source, qr/my \$version = 0\.34;/,
	'strain postprocessing output cleanup increments its workflow version');
like($strain2_source, qr/-withinSpecies 1 -AAtree 0/,
	'legacy within-MGS tree jobs explicitly enable within-species filtering');
like($strain2_source,
	qr/Strain postprocessing v\$version.*?Mode:.*?Inputs:.*?Paths:.*?Resources:.*?Metadata:.*?Association tests:/s,
	'strain postprocessing starts with a structured runtime configuration header');
like($strain2_source, qr/R-analysis plan:.*?existing result\(s\) reused/s,
	'strain postprocessing summarizes per-MGS submission planning');
unlike($strain2_source, qr/print "\$d: ";/,
	'strain postprocessing no longer emits an unlabelled fragment for every MGS');
unlike($strain2_source, qr/Running network of shared strains\.\.\\n\$cmd|print \$cmdPic/,
	'strain postprocessing does not echo full routine commands');
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
like($strain2_source,
	qr/if \(\$rewriteRanalysis\).*?if \(\$doSubmit\) \{.*?remove_tree\(\$within\).*?Dry run: existing strain2 results and checkpoints will be preserved/s,
	'a strain postprocessing dry run does not perform global rewrite cleanup');
like($strain2_source, qr/if \(\$doSubmit && -d \$destD\) \{\s+unlink/s,
	'a strain postprocessing dry run preserves incomplete per-MGS output');
like($strain2_source, qr/if \(!\$rewriteRanalysis && -s \$analysisReport/,
	'a dry rewrite still prepares jobs for previously completed analyses');
like($strain2_source,
	qr/rm -f "\.join\(" ", map \{ shellQuote\(\$_\) \} \(\$analysisLog, \$analysisReport, \$analysisStone\)\)/,
	'an attempted R-analysis job clears only its own stale completion files at execution');
like($strain2_source,
	qr/my \$wrHead=1;.*?\$analysisAttempted\{\$d\}.*?\$wrHead=0/s,
	'the first newly attempted R analysis owns the summary header');
like($strain2_source,
	qr/existingSummaryHeader.*?headerOwner.*?next if \$summaryHeader ne '' && \$line eq \$summaryHeader/s,
	'incremental summary assembly retains one canonical header and removes duplicates');
unlike($strain2_source, qr/\$wrHead=1 if \( \$cnt == 0\)/,
	'R-analysis header output is no longer tied to the first directory index');
like($strain2_source,
	qr/test -s "\.shellQuote\(\$treewasOutfile\).*?test -s "\.shellQuote\(\$summaryOutfile\).*?touch "\.shellQuote\(\$treeWasStone\)/s,
	'treeWAS validates both result files before publishing its checkpoint');
like($strain2_source,
	qr/rm -f "\.join\(" ", map \{ shellQuote\(\$_\) \}.*?\$treewasOutfile, \$summaryOutfile, \$treeWasStone/s,
	'treeWAS removes stale result files before validating a fresh attempt');
like($strain2_source,
	qr/my \(\$networkDep, \$networkStone\) = strainNetwork\(\);.*?my \(\$treeWasDep, \$treeWasStone\) = treeWas\(\);.*?qsubSystemJobAlive\(\\\@postAnalysisJobs, \$QSBoptHR\).*?missingPostAnalysisStones/s,
	'strain postprocessing waits for downstream jobs and validates their checkpoints');
like($strain2_source, qr/return \(\$dep, \$networkStone\).*?return \(\$dep, \$treeWasStone\)/s,
	'network and treeWAS helpers return dependencies and checkpoint paths');
my $resort_source = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'MGS', 'resortMGSgenes4importance.pl'));
like($resort_source, qr/print O evalCurMGS\(""\) if \$curMGS ne "";/,
	'gene-priority resorting flushes its final MGS at EOF');
like($resort_source, qr/compl\.incompl\.\$clusterID\.fna\.clstr\.idx/,
	'gene-priority resorting uses the propagated catalog identity');
my $mgs_source = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'MGS.pl'));
my ($mgs_main, $mgs_subroutines) = split /# Subroutines\n/, $mgs_source, 2;
ok(defined($mgs_subroutines), 'MGS has a distinct subroutine section after its main routing');
unlike($mgs_main, qr/^sub \w+\s*\{/m,
	'MGS keeps subroutine definitions below its main routing');
like($mgs_source, qr/my \$MGSpipelineVersion = 0\.51;/,
	'MGS version includes marker-set-aware between-MGS phylogeny');
like($mgs_source,
	qr/Starting MGS pipeline v\$MGSpipelineVersion.*?GetOptions\(.*?open LOG,.*?Configuration accepted; loading mapping and catalogue metadata.*?my \@checkpointInputs.*?getDirsPerAssmblGrp/s,
	'MGS displays startup configuration before loading input metadata');
like($mgs_source,
	qr/sub _exclude_empty_samples.*?SMPL\.empty.*?delete \$groups->\{\$group\}.*?delete \$map->\{\$_\}/s,
	'MGS removes empty samples and all-empty assembly groups before downstream path handling');
like($mgs_source,
	qr/_exclude_empty_samples\(\\%DOs, \\%map\).*?No non-empty samples were found/s,
	'MGS applies empty-sample exclusion before requiring eligible input');
like($mgs_source, qr/\$checkpointParameters\{empty_samples\} = join\(',', \@emptySamples\)/,
	'MGS checkpoint provenance records the excluded empty samples');
like($mgs_source,
	qr/_count_lines_up_to\(20, \@marker_lca_files\).*?last if \$count >= \$limit/s,
	'MGS marker-LCA preflight stops as soon as its minimum evidence is established');
like($mgs_source,
	qr/MB2assignedBinIds\(\$MBout.*?sub getGoodMBstats.*?MB2assignedBinIds\(\$MBf/s,
	'MGS validation and MAG statistics use ID-only bin assignment scans');
like($mgs_source,
	qr/sub CanopyPrep.*?my %canCnts.*?close \$canopy_input.*?open \$canopy_input.*?kept_canopies/s,
	'MGS Canopy filtering counts first and streams qualifying records in a second pass');
unlike($mgs_source, qr/sub CanopyPrep.*?%can2gene/s,
	'MGS Canopy filtering no longer stores every gene membership in memory');
unlike($mgs_source,
	qr/job_resources\.tsv|timed_job_command|\/usr\/bin\/time|resource_start_seconds/,
	'MGS does not add inline resource-accounting wrappers to submitted commands');
my $cluster_mags_source = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'MGS', 'clusterMAGs.pl'));
my $markers_per_mgs_source = slurp(File::Spec->catfile(
	$Bin, '..', 'secScripts', 'MGS', 'markersPerMGS.pl',
));
like($mgs_source,
	qr/"perlClusterMAGs!" => \\\$perlClusterMAGs.*?if \(\$perlClusterMAGs\).*?getProgPaths\("clusterMGS_scr"\).*?else \{.*?getProgPaths\("clusterMAGs"\)/s,
	'MGS exposes the Perl clustering implementation only through an explicit flag and defaults to the binary');
like($mgs_source,
	qr/perl_cluster_mags\s+=> \$perlClusterMAGs \? 1 : 0/,
	'MGS checkpoints record the selected MAG clustering engine');
like($cluster_mags_source,
	qr/"perlClusterMAGs!" => \\\$perlClusterMAGs.*?if \(!\$perlClusterMAGs\).*?systemW \$cmd.*?exit;.*?Entering the explicitly requested Perl clusterMAGs compatibility algorithm/s,
	'the clusterMAGs wrapper enters its Perl algorithm only when explicitly requested');
unlike($mgs_source . $cluster_mags_source . $markers_per_mgs_source,
	qr/LOGandSUB\/MAGvsGC\.txt\.gz/,
	'all MGS producers and consumers use the Bin directory MAG report path');
my @compressed_report_checks =
	($mgs_source . $cluster_mags_source) =~ /test -s \$outD\/MAGvsGC\.txt\.gz/g;
is(scalar(@compressed_report_checks), 2,
	'both binary clustering routes validate the compressed MAG report in the Bin directory');
like($cluster_mags_source,
	qr/open OX,">\$outD\/MAGvsGC\.txt".*?unless -s "\$outD\/MAGvsGC\.txt\.gz"/s,
	'the Perl compatibility path also publishes and validates the canonical report');
like($markers_per_mgs_source,
	qr/\$MAGrep = "\$GC_bin_dir\/MAGvsGC\.txt\.gz"/,
	'marker aggregation infers the canonical Bin directory report');
like($cluster_mags_source,
	qr/sub mapsWithoutEmptySamples.*?SMPL\.empty.*?next if exists \$empty\{\$sample\}/s,
	'MAG clustering removes marked empty samples from its derived input maps');
like($cluster_mags_source, qr/-map \$clusteringMapF/,
	'the clustering binary receives only the non-empty sample maps');
like($mgs_source,
	qr/MGS pipeline v\$MGSpipelineVersion.*?Mode:.*?Inputs:.*?Paths:.*?Clustering:.*?Quality:.*?Resources:.*?Optional analyses:/s,
	'MGS starts with a structured runtime configuration header');
unlike($mgs_source, qr/print \$cmd\."\\n"|printL \$ph2Cmd|print \$ph1OUT/,
	'MGS does not echo full routine or scheduler commands');
like($mgs_source,
	qr/if \(\$useGTDBmg eq "FMG"\).*?include custom reference genomes.*?\$baseTreeCmd -outD \$outD\/customRefs\/ -refGenos \[refs\].*?Custom reference genomes are unavailable for predefined GTDB-marker trees/s,
	'MGS offers custom references only for the FMG extraction path');
like($mgs_source, qr/MAG availability summary:.*?Missing\/empty MAG examples:.*?complete list is retained in pipeline\.log/s,
	'MGS summarizes missing MAGs while retaining their complete list in the log');
unlike($mgs_source, qr/printL "no MAG file: \$MBout/,
	'MGS no longer prints one missing-MAG path per assembly group');
like($mgs_source, qr/my \$useGTDBmg = "GTDB";/,
	'MGS uses its documented GTDB marker default');
like($mgs_source,
	qr/%checkpointParameters = \(.*?marker_set\s+=> \$useGTDBmg.*?binner\s+=> \$binSpeciesMG.*?quality_checker/s,
	'MGS checkpoints fingerprint marker-set, binner, and quality-checker choices');
like($mgs_source, qr/return 0 unless defined\(\$file\) && -s \$file;.*?checkpoint_valid/s,
	'MGS does not accept provenance-free empty checkpoint stones');
like($mgs_source, qr/checkpointInputs.*?catalog_fna.*?catalog_faa.*?map_/s,
	'MGS checkpoints fingerprint primary catalogue, matrix, and map inputs');
like($mgs_source,
	qr/stage1ProvenanceInvalid\s+=.*?!_checkpoint_valid\(\$st1ston\).*?if \(\$rewrClusterMAGs \|\| \$stage1ProvenanceInvalid\).*?my \$ph1flag =.*?_checkpoint_valid\(\$st1ston\)/s,
	'a clustering-input mismatch invalidates cluster products instead of blessing stale Stage-I output');
like($mgs_source, qr/'gtdb-taxonomy'.*?\$finalClustersFilt.*?'mgs-abundance'.*?\$finalClustersFilt.*?'marker-mgs-abundance'.*?\$finalClustersFilt/s,
	'downstream taxonomy and abundance checkpoints track the current MGS membership');
like($mgs_source, qr/binExtractionValid.*?remove_tree.*?geneBinFiles.*?contigBinFiles.*?_touch_checkpoint/s,
	'bin extraction rebuilds cleanly and records its concrete genome outputs');
like($mgs_source,
	qr/binExtractionValid.*?_representative_contig_outputs_valid.*?sub _representative_contig_outputs_valid.*?fa\|fna\|fasta.*?\\.gz/s,
	'bin extraction rejects stale uncompressed or malformed representative-contig outputs');
like($mgs_source, qr/if \(\$rewrClusterMAGs \|\| \$stage1ProvenanceInvalid\).*?invalidate downstream checkpoint.*?between_phylo.*?within_phylo/s,
	'reclustering invalidates dependent checkpoints and phylogenies');
like($mgs_source, qr/"clusterID=i" => \\\$clusterID/,
	'MGS accepts a gene-catalog cluster identity');
like($mgs_source,
	qr/-clusterID \$clusterID.*?-geneCatIdx \$GCd\/compl\.incompl\.\$clusterID\.fna\.clstr\.idx/s,
	'MGS passes cluster identity to both Perl and binary MAG clustering routes');
like($mgs_source, qr/-MGset \$useGTDBmg -clusterID \$clusterID -maxCores \$canCore/,
	'MGS passes cluster identity to strain analysis');
like($mgs_source,
	qr/\$cmdSI2 = "\$MMLscr -GCd \$GCd -cores \$numCore -MGset \$useGTDBmg -Binner \$BinnerShrt -binD \$outD/,
	'MGS forwards its selected marker set and configured small-core count to marker abundance');
like($mgs_source,
	qr/qsubSystem\(\$logDir\."\/abundMGS_core\.sh".*?\$cmdSI2,.*?\$numCore,"100G"/s,
	'marker abundance retains a fixed total-memory request as core count changes');
like($mgs_source,
	qr/qsubSystem\(\$logDir\."\/abundMGS\.sh".*?\$cmdSI,.*?1,/s,
	'MGS consistently runs specI abundance through the configured submission backend');
unlike($mgs_source, qr/my \@files = glob \("\$GCd\/FMG\/tax\/\*tmp\.m8"\)/,
	'MGS abundance scheduling is not gated by hard-coded FMG temporary alignments');
unlike($mgs_source, qr/systemW\s+\$cmdSI\b/,
	'MGS does not execute the expensive specI abundance command directly on the launcher');
like($mgs_source,
	qr/Stage I clustering done.*?_touch_checkpoint\(\$st1ston.*?Preparing between-MGS phylogeny immediately after MGS creation.*?#get checkM quality/s,
	'MGS starts between-MGS inference immediately after publishing the Stage-I MGS set');
like($mgs_source, qr/\$phyloBetween .*?-MGset \$useGTDBmg .*?-MSAprogram 4/,
	'MGS forwards its selected predefined marker set to the between-MGS phylogeny');
like($mgs_source,
	qr/-visualize 0.*?Visualization needs the abundance matrix.*?interMGSphyloViz\.sh.*?MGSphyloViz", \$treedep/s,
	'MGS defers only abundance-dependent tree visualization and preserves the tree dependency');
unlike($mgs_source, qr/die "Kraken MGS taxonomy stage incomplete/,
	'missing optional Kraken taxonomy no longer aborts MGS');
like($mgs_source, qr/warn "Optional Kraken MGS taxonomy stage incomplete; continuing without Kraken-derived MGS taxonomy/,
	'missing optional Kraken taxonomy is still reported');
like($mgs_source, qr/if \(!-s \$krakenInput\).*?skipping MGS Kraken taxonomy:.*?else \{.*?qsubSystem\(\$logDir\."\/krak2MGS\.sh"/s,
	'MGS does not submit the optional Kraken taxonomy job without its input');
like($mgs_source, qr/my \$profileSamples = _matrix_sample_count\("\$GCd\/Matrix\.mat\.gz"\)/,
	'MGS bases Canopy eligibility on the actual abundance-matrix columns');
like($mgs_source, qr/Requested Canopy assignments are missing or empty; continuing without Canopy MGS/,
	'MGS treats an absent optional Canopy result as a MAG-only run');
like($mgs_source, qr/sub _finish_without_mgs.*?no-usable-mgs.*?exit 0/s,
	'MGS has an explicit successful no-usable-MGS terminal state');
like($mgs_source, qr/no assigned MAG passed the minimum 60% completeness\/10% contamination screen/,
	'MGS stops before MAG clustering when sparse input contains no minimally usable bins');
like($mgs_source, qr/_write_single_mgs_observations\(\$finalClusters2, \$observation_file\)/,
	'a single MGS receives the observation table omitted by the clustering implementation');
like($mgs_source, qr/Weighted MGS assignments were not produced; retaining the valid unweighted assignments/,
	'a valid sparse unweighted result no longer requires a weighted alternative');
like($mgs_source, qr/my \$coreMGSCount = _mgs_count\(\$finalClustersFilt\).*?no MGS retained any core genes/s,
	'an empty post-filter result becomes a reported no-MGS outcome');
unlike($mgs_source, qr/scalar\s*\(?_mgs_ids/,
	'MGS counts list-returning assignment IDs without imposing scalar context on sort');
my ($mgs_count_helpers) = $mgs_source =~ /(sub _mgs_ids \{.*?^\}\s+sub _mgs_count \{.*?^\})/ms;
ok(defined $mgs_count_helpers, 'MGS assignment-count helpers can be isolated for testing');
my $mgs_helpers_loaded = eval "$mgs_count_helpers\n1;";
ok($mgs_helpers_loaded, 'MGS assignment-count helpers compile independently') or diag($@);
my $countable_mgs = File::Spec->catfile($tmp, 'countable.clusters');
write_file($countable_mgs, "Bin\tGene\nMGS1\tg1\nMGS1\tg2\n");
is(main::_mgs_count($countable_mgs), 1,
	'MGS assignment counting returns the number of distinct bins, not scalar-sort output');
like($mgs_source, qr/Activating the only available weighted MGS assignments.*?\$weightedMGSCount = 0;.*?\$activatedOnlyWeighted = 1;/s,
	'an already-activated weighted-only result is not moved a second time');
like($mgs_source, qr/if \(\$coreMGSCount < 3\).*?Skipping between-MGS phylogeny/s,
	'MGS does not launch a phylogeny with fewer than three taxa');
like($mgs_source, qr/\$ph2Cmd \.= "-MGSphylo \$iniTree " if -s \$iniTree \|\| \$treedep ne ""/,
	'strain analysis receives an existing or pending between-MGS tree but tolerates a sparse skip');
like($mgs_source, qr/make_path\(\$binD, \$binDctg, \$binDctgFam\)/,
	'MGS creates genome output directories without shell interpolation');

my $marker_source = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'MGS', 'markersPerMGS.pl'));
like($marker_source, qr/use File::Path qw\(make_path\)/,
	'marker abundance uses the filesystem API to create its output directory');
like($marker_source, qr/\@cntPerCat = sort \{ \$a <=> \$b \} \@cntPerCat/,
	'marker occurrence statistics are sorted numerically');
like($marker_source, qr/\$gene2MGS\{\$_\}\{\$MGS\} = 1/,
	'marker-to-MGS links are deduplicated before ambiguity filtering');
like($marker_source, qr/if \(\@MGSs > 1\) \{\s+\$ambGenes\+\+;\s+next;/s,
	'ambiguous markers are excluded instead of counted for every candidate MGS');
like($marker_source, qr/my \@MGSs = sort keys %\{\$gene2MGS\{\$mark\}\};\s+# Shared markers.*?next if \@MGSs > 1;/s,
	'ambiguous markers are also excluded from taxonomy correction');

my $cluster_wrapper_source = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'MGS', 'clusterMAGs.pl'));
like($cluster_wrapper_source, qr/my \$mapF = resolve_catalog_maps\(\$GCd\)/,
	'clusterMAGs resolves every map through the catalog map manifest');
unlike($cluster_wrapper_source, qr/else\{\s+die "can't find indir \$inD/s,
	'clusterMAGs no longer aborts before its catalog map resolution');
like($cluster_wrapper_source,
	qr/if \(\$logDir eq ""\)\{\$logDir = \$outD\."LOGandSUB\/";\}.*?make_path\(\$outD, \$logDir\)/s,
	'clusterMAGs defaults and creates its output and log directories');
like($cluster_wrapper_source,
	qr/if \(\$redo && -e .*?\) \{\s+unlink .*?or die/s,
	'clusterMAGs redo removes its prior cluster with checked filesystem operations');
like($mgs_source, qr/my \$mapF = resolve_catalog_maps\(\$GCd\)/,
	'MGS resolves single and multiple maps through inmap.txt');
like($mgs_source, qr/my \$specIoutDir = "\$annoDir\/\$\{COGdir\}_MGS"/,
	'MGS-specific annotation output belongs to the selected Bin directory');
like($mgs_source, qr/-outD \$specIoutDir/,
	'MGS passes its binner-local output directory to the annotation worker');
like($mgs_source, qr/\$tmpD \.= "\/" unless \$tmpD =~ m\{\/\$\}/,
	'MGS normalizes a custom temporary directory before adding assembly groups');
unlike($mgs_source, qr/Anno\/Tax\/\$\{COGdir\}_MGS/,
	'MGS no longer owns active MGS annotation output in a shared catalog directory');
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
like($strain_source, qr/rename \$mergeFileByName\{\$name\}, \$outfile or die/,
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
like($strain_source, qr/abundance_pattern_mask\(\\\@abunGs,/,
	'within-MGS abundance filtering uses the configurable robust pattern mask');

my $build_tree_source = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'phylo', 'buildTree5.pl'));
like($build_tree_source, qr/\(\?<sample>\.\*\?\).*?\(\?<gene>\.\+\)/,
	'tree sequence identifiers split at the first separator and retain compound locus names');
like($build_tree_source, qr/sub geneFileStem.*?sprintf\("_%02X", ord\(\$1\)\)/s,
	'compound locus names are encoded safely and deterministically for downstream filenames');

my $gene_cat_source = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'geneCat.pl'));
like($gene_cat_source, qr/my \$version = 0\.55;/,
	'geneCat output cleanup increments its workflow version');
like($gene_cat_source,
	qr/GeneCat pipeline v\$version.*?Mode:.*?Inputs:.*?Paths:.*?Clustering:.*?Resources:.*?Batches:.*?Downstream:/s,
	'geneCat starts with a structured runtime configuration header');
unlike($gene_cat_source, qr/print \$cmd\."\\n\\n"|print "\$qcmd\\n"/,
	'geneCat does not echo full routine or scheduler commands');
unlike($gene_cat_source, qr/getProgPaths\("specIphylo_scr"\)/,
	'geneCat no longer owns or resolves the obsolete interspecies-tree launcher');
like($gene_cat_source,
	qr/\$missingGffCount <= 5.*?No more missing-GFF examples.*?has \$missingGffCount gene\(s\) without GFF entries/s,
	'geneCat caps missing-GFF examples and reports their aggregate count');
like($gene_cat_source, qr/Per-sample completeness output is now shown every 25 assemblies.*?GeneCompleteness/s,
	'geneCat throttles console completeness progress while identifying the detailed report');
like($gene_cat_source, qr/No more incomplete-input examples will be shown.*?Found " \. scalar\(\@probSample\)/s,
	'geneCat bounds repetitive preflight diagnostics and retains an affected-sample summary');
like($gene_cat_source, qr/Protein extraction: \$proteinSamplesDone\/\$proteinSampleCount samples/s,
	'geneCat reports periodic aggregate protein-extraction progress');
like($gene_cat_source, qr/\$doubleSmplWarnCount <= 5.*?No more duplicate-sample-name examples/s,
	'geneCat bounds duplicate-sample-name conflicts');
like($gene_cat_source, qr/\$wrongSmplWarnCount <= 5.*?No more unexpected-sample-name examples/s,
	'geneCat bounds unexpected-sample-name conflicts');
unlike($gene_cat_source, qr/print STDERR \$errStr/,
	'geneCat does not print duplicate-sample details both immediately and again in its summary');
my $submission_source = slurp(File::Spec->catfile($Bin, '..', 'Mods', 'Subm.pm'));
like($submission_source, qr/Submitted \$requestedJobName as job \$jname/,
	'scheduler submissions emit one concise stage-and-job-ID line');
unlike($submission_source, qr/print "SUB:\$jname\\t"/,
	'scheduler submissions no longer emit incomplete inline fragments');
like($gene_cat_source, qr/\$matrixSampleCount = _matrix_sample_count\(\$matrixFile\)/,
	'geneCat uses produced matrix cardinality rather than raw map cardinality for Canopy');
like($gene_cat_source, qr/matrix_sample_count=\$\(gzip -cd .*?declutter-skipped-low-sample-count/s,
	'fire-and-forget decluttering defers cardinality checks until the matrix exists');
like($gene_cat_source, qr/canopy-skipped-low-sample-count.*?SKIPPED\.txt/s,
	'low-sample Canopy skips are checkpointed with a durable explanation');
like($gene_cat_source, qr/Canopy clustering completed but found no clusters.*?SKIPPED\.txt/s,
	'a successful Canopy run with no biological clusters is not treated as a tool crash');
like($gene_cat_source, qr/Kraken classified no catalog genes; leaving taxonomy outputs empty/s,
	'an empty Kraken classification does not invoke aggregation with an empty reference');

my $filter_source = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'R_scripts', 'filterMB2.R'));
like($filter_source, qr/if \(nrow\(M\) == 0\).*?quit\(save="no", status=0\)/s,
	'MGS core filtering publishes an empty biological result successfully');
like($filter_source, qr/Mext = M\[extCriteria,,drop=FALSE\]/,
	'MGS core filtering preserves table shape for a single retained row');
like($filter_source, qr/extCriteria\[is\.na\(extCriteria\)\] = FALSE/,
	'MGS extended-core filtering rejects incomplete logical rows explicitly');

done_testing();
