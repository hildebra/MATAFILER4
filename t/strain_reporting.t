use strict;
use warnings;
no warnings 'once';
use FindBin qw($Bin);
use lib "$Bin/..", "$Bin/lib";
use MFTestConfig;
use File::Path qw(make_path);
use File::Temp qw(tempdir);
use IO::Compress::Gzip qw(gzip $GzipError);
use Test::More;
use Mods::StrainSampleStats qw(count_msa_samples);

sub write_file {
	my ($path, $text) = @_;
	make_path(File::Basename::dirname($path));
	open my $out, '>', $path or die "$path: $!";
	print {$out} $text;
	close $out or die "$path: $!";
}
sub slurp {
	open my $in, '<', $_[0] or die "$_[0]: $!";
	local $/;
	return <$in>;
}
sub read_table {
	my @lines = split /\n/, slurp($_[0]);
	my @columns = split /\t/, shift @lines;
	return [map {
		my @values = split /\t/, $_, -1;
		my %row; @row{@columns} = @values; \%row;
	} @lines];
}
sub compressed {
	my ($path, $contents) = @_;
	make_path(File::Basename::dirname($path));
	gzip \$contents => $path or die $GzipError;
	return $path;
}

my $tmp = tempdir(CLEANUP => 1);
my $locus1 = compressed("$tmp/locus1.fna.gz",
	">sampleA|cog|g1 annotation\nA--N\n>sampleA|cog|duplicate\nTTTT\n>masked|cog|g\nNN--\n>outgroup|cog|g\nAAAA\n");
my $locus2 = compressed("$tmp/locus2.fna.gz",
	">sampleA|cog|g2\nCCCC\n>sampleB|cog|g2\n--\naC\n>outgroup|cog|g2\nAAAA\n");
my $aa = compressed("$tmp/locus1.faa.gz", ">protein_only\nAAAA\n");
my $syn = compressed("$tmp/locus1.syn.fna.gz", ">subset_only\nAAAA\n");
my $merged = compressed("$tmp/MSAli.fna.gz", ">merged_only\nAAAA\n");
is_deeply(count_msa_samples([$locus1, $locus2, $aa, $syn, $merged], 'outgroup'),
	{ msa_samples => 2, msa_outgroup_samples => 1 },
	'MSA sample union excludes outgroup, masked records, AA duplicates and derived alignments');
is_deeply(count_msa_samples([$aa], ''),
	{ msa_samples => 1, msa_outgroup_samples => 0 }, 'AA-only alignment counts are supported');
is_deeply(count_msa_samples([$locus1], 'sampleA'),
	{ msa_samples => 1, msa_outgroup_samples => 1 }, 'outgroup comparison uses the parsed sample ID');

my $source = slurp("$Bin/../secScripts/MGS/strain_within.pl");
my @helpers;
for my $name (qw(mgsTreePath mgsOutputComplete preparedOutgroupLog msaOnlyArtifactsReady writeMGSSampleHistograms
	writeSelectionAttritionSummary writeGeneLengthSampleSummary writeStrainSummary printSampleStatsSummary)) {
	my ($sub) = $source =~ /(sub \Q$name\E \{.*?^\})/ms;
	BAIL_OUT("Cannot find $name") unless $sub;
	push @helpers, $sub;
}
my $loaded = eval <<'SETUP' . join("\n", @helpers) . "\n1;";
package StrainReportFixture;
use File::Spec;
use Digest::SHA qw(sha256_hex);
use Mods::GenoMetaAss qw(gzipopen fileGZe fileGZs mean);
use Mods::math qw(medianArray);
use Mods::StrainSampleStats qw(count_msa_samples);
use Mods::WorkflowResilience qw(retry_open retry_close retry_rename atomic_write_text);
our ($onlyMSA, $strictBackbone, $phyloProg, $LOGDIR, $outD, $version,
	$recoveryLogName, $summaryLogName, $sampleStatsLogName, $sampleStatsSummaryLogName,
	$scratchD, $doSubmit, $phase1SampleSummary);
our (@specis, @samples, %SIdirs, %ConspecificMGS, %PreferredOutgroup);
sub limitedNotice { }
sub stepProgress { }
SETUP
ok($loaded, 'reporting helpers load independently') or BAIL_OUT($@);

$StrainReportFixture::onlyMSA = 1;
$StrainReportFixture::strictBackbone = 0;
$StrainReportFixture::phyloProg = 1;
$StrainReportFixture::doSubmit = 1;
$StrainReportFixture::LOGDIR = "$tmp/LOGandSUB";
$StrainReportFixture::outD = $tmp;
$StrainReportFixture::scratchD = "$tmp/scratch";
$StrainReportFixture::version = 'test';
$StrainReportFixture::recoveryLogName = 'strainRecovery.tsv';
$StrainReportFixture::summaryLogName = 'strain_within.summary.log';
$StrainReportFixture::sampleStatsLogName = 'strainSampleStats.tsv';
$StrainReportFixture::sampleStatsSummaryLogName = 'strainSampleStats.summary.tsv';
@StrainReportFixture::specis = qw(MGS.cached MGS.legacy MGS.empty MGS.missing);
@StrainReportFixture::samples = qw(sampleA sampleB sampleC);
make_path($StrainReportFixture::LOGDIR);
for my $mgs (@StrainReportFixture::specis) {
	$StrainReportFixture::SIdirs{$mgs} = "$tmp/$mgs";
	make_path("$tmp/$mgs");
}
compressed("$tmp/MGS.cached/MSA/g.fna.gz", ">sampleA|cog|g\nACGT\n");
write_file("$tmp/MGS.cached/msaOnly.complete.tsv",
	"status\tmsa_complete\nmsa_samples\t4\nmsa_outgroup_samples\t1\n");
compressed("$tmp/MGS.legacy/MSA/g.fna.gz",
	">sampleA|cog|g\nACGT\n>sampleB|cog|g\nTTTT\n>outgroup|cog|g\nAAAA\n");
write_file("$tmp/MGS.legacy/msaOnly.complete.tsv", "status\tmsa_complete\n");
write_file("$tmp/MGS.legacy/data.log", "OG:outgroup\n");
write_file("$tmp/MGS.empty/tooFewSamples.sto", "too few samples\n");
# Stale tree metrics must not contaminate an MSA-only run.
write_file("$tmp/MGS.cached/phylo/selection_attrition.tsv",
	"metric\tvalue\nbackbone_samples\t999\nplacement_samples\t55\n");
my $report = StrainReportFixture::writeMGSSampleHistograms();
my $rows = read_table($report->{details});
is(scalar(@$rows), 4, 'every selected MGS has a sample-count row');
my %row = map { $_->{MGS} => $_ } @$rows;
is($row{'MGS.cached'}{msa_samples}, 4, 'new completion markers supply cached MSA counts');
is($row{'MGS.legacy'}{msa_samples}, 'NA', 'legacy outputs stay complete without an expensive alignment recount');
is($row{'MGS.cached'}{backbone_samples}, 'NA', 'MSA mode ignores stale backbone counts');
is($row{'MGS.cached'}{placement_samples}, 'NA', 'placement is inapplicable in MSA mode');
is($row{'MGS.cached'}{tree_status}, 'not_requested', 'MSA output is not called a missing tree');
is($row{'MGS.cached'}{output_status}, 'msa_complete', 'completion and sample measurement are separate');
is($row{'MGS.empty'}{msa_samples}, 0, 'a terminal MGS contributes a measured zero');
is($row{'MGS.missing'}{msa_samples}, 'NA', 'missing output is not counted as zero');
is($row{'MGS.missing'}{output_status}, 'output_missing', 'missing outputs remain visible');
is_deeply($report->{roles}, ['msa'], 'MSA-only histogram has one applicable role');
is_deeply($report->{statistics}{msa}, {
	count => 2, missing => 2, minimum => 0, maximum => 4, median => 2, mean => '2.00',
}, 'histogram denominator includes measured zeros and excludes missing measurements');
my $legacyCached = { %{$row{'MGS.legacy'}}, msa_samples => 2 };
my $hist = read_table($report->{histogram});
is($hist->[0]{MGS_count}, 1, 'zero bin contains the terminal MGS');
is($hist->[0]{fraction}, '0.500000', 'histogram fraction uses the number of measured MGS');
is($hist->[0]{missing_MGS}, 2, 'histogram records its missing denominator');

@StrainReportFixture::specis = ('MGS.missing');
$report = StrainReportFixture::writeMGSSampleHistograms();
is($report->{statistics}{msa}{mean}, 'NA', 'no observations produce NA rather than a zero mean');
is(read_table($report->{histogram})->[0]{fraction}, 'NA', 'empty histogram fractions are unavailable');

@StrainReportFixture::specis = ('MGS.cached');
$StrainReportFixture::onlyMSA = 0;
$StrainReportFixture::strictBackbone = 1;
$report = StrainReportFixture::writeMGSSampleHistograms();
is_deeply($report->{roles}, [qw(backbone placement)], 'backbone mode retains both distributions');
is($report->{statistics}{backbone}{mean}, '999.00', 'tree mode uses selection attrition');
is(read_table($report->{details})->[0]{excluded_samples}, 'NA', 'absent exclusion metrics remain unknown');
$StrainReportFixture::strictBackbone = 0;
$StrainReportFixture::phyloProg = 2;
write_file("$tmp/MGS.cached/treeDone.sto", "done\n");
write_file("$tmp/MGS.cached/phylo/VERYFASTTREE_allsites.nwk", "(a,b,c);\n");
$report = StrainReportFixture::writeMGSSampleHistograms();
is_deeply($report->{roles}, ['tree'], 'ordinary tree mode does not print a redundant placement histogram');
is(read_table($report->{details})->[0]{output_status}, 'tree_complete', 'configured non-IQ-TREE output is recognized');

$StrainReportFixture::onlyMSA = 1;
write_file("$tmp/LOGandSUB/strainRecovery.tsv", join("\n",
	"MGS\tsample\toutcome\treason\tretained_genes\tqc_status\tambiguous_failure\tconspecific_failure\tmosaic_loci",
	"MGS.cached\tsampleA\trecovered\tpassed_qc\t2000\tsingle_strain\t0\t0\t0",
	"MGS.cached\tsampleB\trecovered\tpassed_qc\t2001\tmixed_strain\t1\t0\t2",
	"MGS.cached\tsampleC\tfiltered\ttoo_few_after_abundance\t2\t\t0\t0\t0",
	"MGS.not_selected\tsampleA\trecovered\tpassed_qc\t9000\tsingle_strain\t0\t0\t0", '')."\n");
my $stdout = '';
{
	open my $capture, '>', \$stdout or die $!;
	local *STDOUT = $capture;
	StrainReportFixture::printSampleStatsSummary({samples => 3, processed_samples => 2});
	StrainReportFixture::writeStrainSummary({'eligible MSA-only job' => 1}, {});
}
like($stdout, qr/2 recovered \+ 1 filtered = 3 evaluated pairs/, 'summary scopes reused recovery ledger to selected MGS');
like($stdout, qr/>2000: 1\./, 'cumulative thresholds are strictly greater than N');
like($stdout, qr/2 processed; 1 skipped\/unavailable/, 'resume summary uses saved sample-processing counts');
unlike($stdout, qr/recovered_MAGs|backbone_samples_per_MGS|Placement samples per MGS/, 'summary omits obsolete MAG keys and inactive roles');
like(slurp("$tmp/LOGandSUB/strainSelectionAttrition.tsv"),
	qr/recovery\trecovered_pairs\.loci_gt_2000\t1/, 'exact threshold remains in machine-readable filter totals');
like(slurp("$tmp/LOGandSUB/strainSelectionAttrition.tsv"),
	qr/submission\teligible MSA-only job\t1/, 'submission accounting remains available without duplicate STDOUT');

# Reusing the existing run-wide table avoids legacy alignment rescans, while
# changing the marker invalidates the cached measurement.
{
	no warnings 'redefine';
	local *StrainReportFixture::count_msa_samples = sub { die "unexpected rescan" };
	my %reused;
	ok(StrainReportFixture::msaOnlyArtifactsReady("$tmp/MGS.legacy", \%reused, $legacyCached),
		'unchanged legacy completion marker reuses previously measured counts');
	is($reused{msa_samples}, 2, 'cached legacy ingroup count is preserved');
	write_file("$tmp/MGS.legacy/msaOnly.complete.tsv", "status\tmsa_complete\nupdated\t1\n");
	ok(StrainReportFixture::msaOnlyArtifactsReady("$tmp/MGS.legacy", \%reused, $legacyCached),
		'changed completion marker remains reusable without recounting alignments');
	is($reused{msa_samples}, 'NA', 'changed completion marker invalidates cached counts');
}

# Older outputs may lack the outgroup identity: do not call an unknown total
# an ingroup sample count.
unlink "$tmp/MGS.legacy/data.log" or die $!;
my %unknown;
ok(StrainReportFixture::msaOnlyArtifactsReady("$tmp/MGS.legacy", \%unknown),
	'legacy output remains complete without optional outgroup metadata');
is($unknown{msa_samples}, 'NA', 'absent historical counts are unavailable rather than inferred from extraction');

done_testing();
