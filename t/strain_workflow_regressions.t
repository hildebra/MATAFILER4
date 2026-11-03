use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::GenoMetaAss qw(
	readClstrRev readFasta writeClstrRevBinaryShards readClstrRevBinaryShard
	writeSequenceBinaryCache readSequenceBinaryCache
);

sub write_file {
	my ($path, $contents) = @_;
	open my $fh, '>', $path or die "Cannot write $path: $!";
	print {$fh} $contents or die "Cannot write $path: $!";
	close $fh or die "Cannot close $path: $!";
}

sub slurp {
	my ($path) = @_;
	open my $fh, '<', $path or die "Cannot read $path: $!";
	local $/;
	return <$fh>;
}

my $tmp = tempdir(CLEANUP => 1);
my $fasta = File::Spec->catfile($tmp, 'records.fa');
write_file($fasta, <<'FASTA');
>keep1 D=8 CSP=0.01
AAAA
>drop D=2 CSP=0.50
CCCC
>keep3 D=5 CSP=0.02
GGGG
FASTA

my %wanted = (keep1 => 1, keep3 => 1);
is_deeply(
	readFasta($fasta, 1, '\\s', \%wanted),
	{ keep1 => 'AAAA', keep3 => 'GGGG' },
	'FASTA subset selection applies independently to intermediate and final records',
);
is_deeply(
	readFasta($fasta, 0, '\\s', \%wanted),
	{
		'keep1 D=8 CSP=0.01' => 'AAAA',
		'keep3 D=5 CSP=0.02' => 'GGGG',
	},
	'FASTA subset lookup can use short IDs while retaining full headers',
);

my $glob_dir = File::Spec->catdir($tmp, 'glob');
mkdir $glob_dir or die "Cannot create $glob_dir: $!";
write_file(File::Spec->catfile($glob_dir, 'a_empty.fa'), '');
write_file(File::Spec->catfile($glob_dir, 'b_records.fa'), ">later\nACGT\n");
is_deeply(
	readFasta(File::Spec->catfile($glob_dir, '*.fa'), 1, '\\s'),
	{ later => 'ACGT' },
	'an empty member of a FASTA glob does not suppress later files',
);

my $cluster_index = File::Spec->catfile($tmp, 'cluster.idx');
write_file($cluster_index, "seed1\tsample1__gene1,sample2__gene2\n");
my (undef, $empty_cluster_subset) = readClstrRev($cluster_index, 0, {}, {});
is_deeply(
	$empty_cluster_subset,
	{},
	'an explicitly empty cluster-member subset does not fall back to the complete catalogue',
);

write_file($cluster_index, join("\n",
	"seed1\tsample1__gene1,sample2__gene2,alias1__gene3",
	"seed2\tsample2__gene4,sample3__gene5",
	"drop\tsample1__ignored",
).'\n');
my @binary_shards = map { File::Spec->catfile($tmp, "cluster.worker.$_.bin") }
	0 .. 2;
my $shard_fingerprint = 'a' x 64;
my $shard_metadata = writeClstrRevBinaryShards(
	$cluster_index,
	{ seed1 => 1, seed2 => 1 },
	{ sample1 => 0, alias1 => 0, sample2 => 1, sample3 => 1 },
	\@binary_shards,
	$shard_fingerprint,
);
is_deeply(
	[map { $_->{records} } @{$shard_metadata}],
	[1, 2, 0],
	'binary cluster-index publication records only selected clusters in each worker partition',
);
is_deeply(
	readClstrRevBinaryShard($binary_shards[0], $shard_fingerprint, 0, 3),
	{ seed1 => 'sample1__gene1,alias1__gene3' },
	'binary cluster-index shard preserves member order and catalogue aliases',
);
is_deeply(
	readClstrRevBinaryShard($binary_shards[1], $shard_fingerprint, 1, 3),
	{ seed1 => 'sample2__gene2', seed2 => 'sample2__gene4,sample3__gene5' },
	'binary cluster-index shard contains only the assigned worker members',
);
is_deeply(
	readClstrRevBinaryShard($binary_shards[2], $shard_fingerprint, 2, 3),
	{},
	'an empty worker receives a valid binary shard rather than falling back to the full index',
);
my $wrong_shard_provenance = eval {
	readClstrRevBinaryShard($binary_shards[0], 'b' x 64, 0, 3);
	1;
};
ok(!$wrong_shard_provenance && $@ =~ /Invalid binary cluster-index shard header/,
	'binary cluster-index reader rejects a shard from another provenance generation');
my $corrupt_shard = File::Spec->catfile($tmp, 'cluster.worker.corrupt.bin');
my $corrupt_contents = slurp($binary_shards[0]);
substr($corrupt_contents, -1, 1) = chr(ord(substr($corrupt_contents, -1, 1)) ^ 1);
write_file($corrupt_shard, $corrupt_contents);
my $corrupt_shard_accepted = eval {
	readClstrRevBinaryShard($corrupt_shard, $shard_fingerprint, 0, 3);
	1;
};
ok(!$corrupt_shard_accepted && $@ =~ /payload digest mismatch/,
	'binary cluster-index reader rejects payload or trailer corruption');

my $protein_cache = File::Spec->catfile($tmp, 'catalog.proteins.bin');
my $protein_fingerprint = 'c' x 64;
my $protein_metadata = writeSequenceBinaryCache(
	{ seed1 => 'MPEPTIDE', seed2 => 'MSECOND', empty => '' },
	$protein_cache,
	$protein_fingerprint,
);
is_deeply(
	$protein_metadata,
	{ records => 2, bytes => -s $protein_cache },
	'binary sequence cache records only nonempty reference proteins and reports its durable size',
);
is_deeply(
	readSequenceBinaryCache($protein_cache, $protein_fingerprint),
	{ seed1 => 'MPEPTIDE', seed2 => 'MSECOND' },
	'binary sequence cache round-trips the common catalogue-protein subset exactly',
);
my $wrong_protein_provenance = eval {
	readSequenceBinaryCache($protein_cache, 'd' x 64);
	1;
};
ok(!$wrong_protein_provenance && $@ =~ /Invalid binary sequence-cache header/,
	'binary sequence cache rejects a different selected-catalogue generation');
my $corrupt_protein_cache = File::Spec->catfile($tmp, 'catalog.proteins.corrupt.bin');
my $corrupt_protein_contents = slurp($protein_cache);
substr($corrupt_protein_contents, -1, 1) =
	chr(ord(substr($corrupt_protein_contents, -1, 1)) ^ 1);
write_file($corrupt_protein_cache, $corrupt_protein_contents);
my $corrupt_protein_accepted = eval {
	readSequenceBinaryCache($corrupt_protein_cache, $protein_fingerprint);
	1;
};
ok(!$corrupt_protein_accepted && $@ =~ /payload digest mismatch/,
	'binary sequence cache rejects same-size protein-cache corruption');

my $strain = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'MGS', 'strain_within.pl'));
my $strain2 = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'MGS', 'strain_within_2.2.pl'));

# Exercise the Phase-I guide identity independently of the full controller. A
# split worker validates the contract before prepRun() assigns $outD, whereas
# the parent records it afterwards. Both moments must select the run-local
# sorted guide even when an older catalogue-side .srt also exists.
my ($phase1_fingerprint_helpers) = $strain =~
	/(sub phase1PathStatComponent \{.*?^\}\n\nsub phase1GuideStatFingerprint \{.*?^\})\n\nsub phase1CatalogStatFingerprint/ms;
BAIL_OUT('Cannot extract Phase-I guide fingerprint helpers')
	unless defined $phase1_fingerprint_helpers;
my $phase1_helpers_loaded = eval <<"PERL";
package TestPhase1GuideFingerprint;
use strict;
use warnings;
use Cwd qw(abs_path);
use Digest::SHA qw(sha256_hex);
use File::Basename qw(basename);
use File::Spec;
our \$phase1InputContractVersion = 3;
our \$outD = '';
sub resolveExistingFile {
	my (\$path) = \@_;
	return -f \$path ? \$path : undef;
}
$phase1_fingerprint_helpers
1;
PERL
ok($phase1_helpers_loaded, 'Phase-I guide fingerprint helpers load independently')
	or diag($@);
my $contract_catalogue = File::Spec->catdir($tmp, 'contract-catalogue');
my $contract_output = File::Spec->catdir($tmp, 'contract-output');
mkdir $contract_catalogue or die "Cannot create $contract_catalogue: $!";
mkdir $contract_output or die "Cannot create $contract_output: $!";
my $contract_guide = File::Spec->catfile(
	$contract_catalogue, 'SB.clusters.core');
write_file($contract_guide, "MGS.1\tgene1\n");
write_file(File::Spec->catfile($contract_catalogue, 'SB.clusters.obs'),
	"MGS.1\t1\n");
write_file("$contract_guide.srt", "MGS.catalogue\tgene-old\n");
write_file("$contract_guide.srt.gene2MGS", "gene-old\tMGS.catalogue\n");
my $staged_sorted = File::Spec->catfile(
	$contract_output, 'SB.clusters.core.srt');
write_file($staged_sorted, "MGS.1\tgene1\n");
write_file("$staged_sorted.gene2MGS", "gene1\tMGS.1\n");
$TestPhase1GuideFingerprint::outD = $contract_output;
my $parent_guide_fingerprint =
	TestPhase1GuideFingerprint::phase1GuideStatFingerprint($contract_guide);
$TestPhase1GuideFingerprint::outD = '';
my $worker_guide_fingerprint =
	TestPhase1GuideFingerprint::phase1GuideStatFingerprint(
		$contract_guide, undef, $contract_output);
is($worker_guide_fingerprint, $parent_guide_fingerprint,
	'pre-initialization worker and initialized parent fingerprint the same run-local guide');
isnt(
	TestPhase1GuideFingerprint::phase1GuideStatFingerprint($contract_guide),
	$parent_guide_fingerprint,
	'a catalogue-side sorted guide has a distinct identity and cannot be selected accidentally');

my ($durable_output_helper_source) = $strain =~
	/(sub strainOutputHasDurablePhaseIState \{.*?\n\})\n\nsub phase1PathStatComponent/s;
ok(defined($durable_output_helper_source),
	'existing-output evidence helper is available for isolated testing');
$durable_output_helper_source =~
	s/\Asub strainOutputHasDurablePhaseIState/sub/;
my $durable_output_helper = eval $durable_output_helper_source;
die "Cannot compile existing-output evidence helper: $@" if $@;

my $fresh_strain_output = File::Spec->catdir($tmp, 'fresh_strain_output');
mkdir $fresh_strain_output or die "Cannot create $fresh_strain_output: $!";
for my $operational (qw(LOGandSUB stones strainsScr1 .scratch)) {
	my $directory = File::Spec->catdir($fresh_strain_output, $operational);
	mkdir $directory or die "Cannot create $directory: $!";
}
ok(!$durable_output_helper->(
		$fresh_strain_output,
		File::Spec->catfile($fresh_strain_output, 'LOGandSUB', 'missing.summary'),
	),
	'a fresh output containing only operational directories permits a subset build');
my $existing_mgs_directory = File::Spec->catdir($fresh_strain_output, 'MGS.1');
mkdir $existing_mgs_directory
	or die "Cannot create $existing_mgs_directory: $!";
ok($durable_output_helper->($fresh_strain_output),
	'an existing per-MGS directory blocks a destructive subset rebuild');
my $summary_evidence = File::Spec->catfile($tmp, 'strainSampleStats.summary.tsv');
write_file($summary_evidence, "durable\n");
ok($durable_output_helper->(File::Spec->catdir($tmp, 'absent_output'), $summary_evidence),
	'a durable Phase-I summary blocks a destructive subset rebuild without a directory scan');

my @analysisPhaseMarkers = (
	q{$waitForAnalysis->('strainStats');},
	q{combineResults(0);},
	q{my ($networkDep, $networkStone) = strainNetwork();},
	q{my ($treeWasDep, $treeWasStone) = treeWas();},
	q{my ($phyloFigureDep, $phyloFigureStone) = visualizeSignPhylos();},
	q{$waitForAnalysis->('popGenStats');},
	q{combineResults(1);},
);
my $previousPhaseMarker = -1;
my $analysisPhasesInOrder = 1;
for my $phaseMarker (@analysisPhaseMarkers) {
	my $phaseMarkerPosition = index($strain2, $phaseMarker);
	$analysisPhasesInOrder = 0
		if $phaseMarkerPosition < 0 || $phaseMarkerPosition <= $previousPhaseMarker;
	$previousPhaseMarker = $phaseMarkerPosition;
}
ok($analysisPhasesInOrder,
	'strain summaries and submitted network/treeWAS/phylogeny work start before the independent population phase is awaited');

my %isolatedAnalysisCommand = (
	strainStats => 'strainStatsR', popGenStats => 'popGenCommand',
);
for my $analysis (sort keys %isolatedAnalysisCommand) {
	my $variable = $isolatedAnalysisCommand{$analysis};
	}
for my $bigTree ('$outD', '$scratchD', '$preConDir', '$outD2', '$scratch_mgs') {
	}

my ($phase1WorkerCommandSource) = $strain =~
	/(sub phase1WorkerCommand \{.*?^\})/ms;
ok(defined($phase1WorkerCommandSource),
	'the Phase-I worker command builder is available for option-contract auditing');
my %phase1WorkerFlag = map { $_ => 1 }
	$phase1WorkerCommandSource =~ /'(-[A-Za-z][A-Za-z0-9]*)'/g;
my @requiredPhase1WorkerFlags = qw(
	-GCd -outD -MGS -clusterID -submit -onlySubmit -maxSubJob
	-MGSminGenesPSmpl -multiGeneSmplMax -conspGeneSmplMax
	-minBadLociPSmpl -presortGenes -maxGenes -treeLocusBudget
	-taxonAwareLocusSelection -disableQC -breakpointGeneFlank
	-abundanceMinLoci -abundanceMinFold -abundanceMaxFold
	-abundanceMaxModifiedZ -prepareMosaicLoci -flushEvery -MGset
	-SNPcaller -minSNPDepth -minSNPCallQual -forceSNPcalls
	-preCompConsSNP -skipIndels -SNPadaptiveQual -SNPdepthFilterScale
	-SNPindelRangeFilt -tmpD -mosaicLoci -MGSabundance -MGSsubset
);
is_deeply(
	[grep { !$phase1WorkerFlag{$_} } @requiredPhase1WorkerFlags],
	[],
	'all mandatory and conditional Phase-I execution flags are propagated to workers',
);

# Stage-I worker memory. The output buffers hold one string per MGS the worker
# has touched, so a sample count cannot bound them; and once the locus model is
# built a worker consults neither the ranked seed index nor most catalogue
# proteins again, because it exits before Phase II.

# Two strain_within runs over one gene catalogue must not share, overwrite or
# delete each other's derived MGS guide. The sorter names its output after the
# guide it is handed, so the guide is staged inside the output directory and
# every product follows it there.

done_testing();
