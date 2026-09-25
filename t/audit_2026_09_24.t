use strict;
use warnings;
use Test::More;
use File::Temp qw(tempdir);
use File::Path qw(make_path);
use File::Spec;
use FindBin qw($Bin);
use IO::Compress::Gzip qw(gzip $GzipError);
use lib "$Bin/..", "$Bin/lib";
use MFTestConfig;
use Mods::GenoMetaAss qw(fileGZs gzipopen);
use Mods::StatsLogReader qw(read_stats_log_excerpt);

# Regressions for the 2026-09-24 audit fixes. Caller bodies are isolated from
# the real scripts, as in local_algorithm_regressions.t.

my $root = File::Spec->rel2abs("$Bin/..");
my $tmp = tempdir(CLEANUP => 1);
sub write_file {
	my ($path, $text) = @_;
	open my $fh, '>', $path or die "$path: $!";
	print {$fh} $text;
	close $fh or die "$path: $!";
}
sub read_file {
	open my $fh, '<', $_[0] or die "$_[0]: $!";
	local $/; return <$fh> // '';
}
sub source_sub {
	my ($source, $name) = @_;
	my ($code) = $source =~ /(^sub \Q$name\E\b[^\{;]*\{.*?^\})/ms;
	die "Cannot isolate $name" unless defined $code;
	return $code;
}

my $main = read_file("$root/MATAF4.pl");
our (%locStats, %MFconfig, %MFopt);
eval join("\n", map { source_sub($main, $_) } qw(
	getContamination getSNPStats getRgStr getGeneStats optiDups remComma
	_sdm_version _sdm_version_at_least _parse_sdm_stats_text
	_sdm_histogram_max_length sdmStatsMany
));
die $@ if $@;

# hostile reports a rounded fraction; the column holds percentages as for kraken
write_file("$tmp/KrakHS.sh.etxt", "hostile run\n");
write_file("$tmp/KrakHS.sh.otxt",
	qq{{\n "reads_in": 1000,\n "reads_out": 975,\n "reads_removed": 25,\n "reads_removed_proportion": 0.025,\n}\n});
my $conta = getContamination("$tmp/KrakHS.sh.etxt", "$tmp/KrakHS.sh.otxt", '');
is($conta->{FilteredContaRdsPerc}, '2.500', 'hostile removal fraction is reported as a percentage');
is($conta->{FilteredContaRds}, 25, 'hostile removed reads are parsed');

# vcf2fna >= 0.44 wording
write_file("$tmp/cons.sh.otxt", "Total bp that can be determined: 100 in 2 entries.\n"
	."  - Found 4 SNPs and 1 INDELS.\n"
	."  - Passing Filters: 4, 0; 1 entries (major, minor SNPs; INDELS). Conflicts resolved: 2\n");
my $snp = getSNPStats("$tmp/cons.sh");
is($snp->{SNP_Passed}, 4, 'passing SNPs are parsed from current vcf2fna output');
is($snp->{INDEL_Passed}, 1, 'passing INDELs are parsed, also with resolved conflicts');

my $rg = getRgStr('S1', 'lib0', 'lib0', 0, 3, 'PB');
like($rg, qr/ID:S1\\tSM:S1/, 'minimap2 read group ID carries the sample name');
unlike($rg, qr/\$smpl/, 'read group ID is not the literal variable name');

# separateContigs.pl publishes GeneStats.txt.gz
make_path("$tmp/cs");
gzip(\"GeneNumber\tAvgGeneLength\n100\t900\t950\t90000\t10000\t80\t5\t5\t10\n"
	=> "$tmp/cs/GeneStats.txt.gz") or die $GzipError;
is(getGeneStats("$tmp/cs/GeneStats.txt")->{GeneNumber}, 100, 'gzipped gene statistics are read');

# markdup statistics now come from the combined map.sh job
make_path("$tmp/log");
write_file("$tmp/log/map.sh.etxt", "samtools markdup stats\nWRITTEN: 900\nDUPLICATE PAIR: 20\n"
	."DUPLICATE SINGLE: 3\nDUPLICATE PAIR OPTICAL: 2\nDUPLICATE SINGLE OPTICAL: 1\nESTIMATED_LIBRARY_SIZE: 5000\n");
write_file("$tmp/log/map2.sh.etxt", "old mapping log without duplicate statistics\n");
my $dups = optiDups("$tmp/log");
is($dups->{PCRduplicates}, 23, 'duplicate statistics are read from map.sh.etxt');
is($dups->{EstLibSize}, 5000, 'library size is read from map.sh.etxt');

# each sdm log has its own length histogram
make_path("$tmp/sdm/LOGandSUB/sdm");
my $sdmDir = "$tmp/sdm/LOGandSUB/sdm";
write_file("$sdmDir/filter.S.log", "sdm 3.53 beta\n  Reads processed:  10\n  Rejected:  0 (0.0%)\n"
	."  Accepted (high quality):  10 (100.0%)\n5%/50%/95% quantiles\n     - sequence Length : 630/708/789\n");
write_file("$sdmDir/filter.S_lenHist.txt", "630\t1\t1\n890\t1\t1\n");
write_file("$sdmDir/filter_lenHist.txt", "150\t1\t1\n2000\t1\t1\n");
my $sdm = sdmStatsMany(["$sdmDir/filter.S.log"], "$tmp/sdm", '', 2000);
is($sdm->{MaxSeqLength}, 890, 'singleton sdm log uses its own length histogram');

# ENA: an extra unpaired file must not become read 1
my $ena = read_file("$root/secScripts/fileManage/ENASRAdl.pl");
{
	no warnings 'redefine';
	eval "use File::Basename qw(basename);\n".source_sub($ena, 'ena_role');
	die $@ if $@;
}
my @uris = ('x/SRR1.fastq.gz', 'x/SRR1_1.fastq.gz', 'x/SRR1_2.fastq.gz');
is_deeply([map { ena_role('PAIRED', $uris[$_], $_, 3) } 0 .. 2], [qw(single r1 r2)],
	'three-file paired ENA runs keep the unpaired file single');
is_deeply([map { ena_role('PAIRED', "x/a$_.fq.gz", $_, 2) } 0, 1], [qw(r1 r2)],
	'two unnamed paired files still fall back to their order');

# geneCat sample batches cover every sample
my $geneCat = read_file("$root/secScripts/geneCat.pl");
my ($toExpr) = $geneCat =~ /my \$locTo = (int\([^;]+\));/;
my ($fromExpr) = $geneCat =~ /my \$locFrom = (int\([^;]+\));/;
ok(defined($toExpr) && defined($fromExpr), 'geneCat batch bounds were found');
my $gaps = 0;
for my $case ([2008, 11], [2019, 11], [1000, 19], [1000, 99], [7, 3]) {
	our ($maxSmpls, $batchNum) = @{$case};
	my $next = 0;
	for our $batch (0 .. $batchNum - 1) {
		my $from = eval $fromExpr; my $to = eval $toExpr;
		$gaps++ if $from != $next;
		$next = $to;
	}
	$gaps++ if $next != $maxSmpls;
}
is($gaps, 0, 'batches are contiguous and end at the last sample');

# eggNOG split: '-' is missing, partial EC numbers are kept
SKIP: {
	skip 'bash/awk unavailable', 3 unless system('bash -c "command -v awk" >/dev/null 2>&1') == 0;
	make_path("$tmp/egg");
	my @header = ('#query', map { "c$_" } 2 .. 21);
	my %annotated = (11 => '2.7.7.-', 19 => 'GT2', 21 => 'ABC_tran');
	write_file("$tmp/egg/ann.tsv", join("\t", @header)."\n"
		.join("\t", 'g1', map { $annotated{$_} // 'x' } 2 .. 21)."\n"
		.join("\t", 'g2', map { '-' } 2 .. 21)."\n");
	system('bash', "$root/secScripts/GC/eggNOG_split.sh", "$tmp/egg/ann.tsv") == 0
		or die "eggNOG_split.sh failed\n";
	is(read_file("$tmp/egg/eggNOGmapper_EC.geneAss"), "g1\t2.7.7.-\n", 'partial EC numbers are kept');
	is(read_file("$tmp/egg/eggNOGmapper_CAZy.geneAss"), "g1\tGT2\n", 'CAZy rows are kept, missing ones dropped');
	is(read_file("$tmp/egg/eggNOGmapper_PFAM.geneAss"), "g1\tABC_tran\n", 'PFAM rows are kept, header dropped');
}

done_testing();
