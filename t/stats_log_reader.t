use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use IO::Compress::Gzip qw(gzip $GzipError);
use Test::More;

use lib File::Spec->catdir(File::Spec->rel2abs('.'));
use Mods::StatsLogReader qw(
	read_stats_log_excerpt
	reset_stats_log_sampling
	stats_log_sampling_summary
);

my $tmp = tempdir(CLEANUP => 1);

sub write_file {
	my ($path, $contents) = @_;
	open my $fh, '>:raw', $path or die "Cannot write $path: $!";
	print {$fh} $contents or die "Cannot populate $path: $!";
	close $fh or die "Cannot close $path: $!";
}

my $small = File::Spec->catfile($tmp, 'small.log');
write_file($small, "header\nsummary\n");
is(read_stats_log_excerpt($small), "header\nsummary\n",
	'small statistics logs are read completely');

my $large = File::Spec->catfile($tmp, 'large.log');
write_file($large,
	"HEAD_MARKER\n" . ('A' x (3 * 1024 * 1024)) . "\n"
	. "MIDDLE_MARKER\n" . ('B' x (3 * 1024 * 1024)) . "\nTAIL_MARKER\n");
reset_stats_log_sampling();
my $excerpt = read_stats_log_excerpt($large);
like($excerpt, qr/^HEAD_MARKER$/m, 'oversized plain-log excerpt retains its header');
like($excerpt, qr/^TAIL_MARKER$/m, 'oversized plain-log excerpt retains its tail');
unlike($excerpt, qr/MIDDLE_MARKER/, 'oversized plain-log excerpt skips its middle');
cmp_ok(length($excerpt), '<=', 1024 * 1024 + 2,
	'oversized plain-log retention is capped at the two configured edges');
my $summary = stats_log_sampling_summary();
is($summary->{large_files}, 1, 'oversized plain log is reported once');
cmp_ok($summary->{retained_bytes}, '<', $summary->{source_bytes},
	'sampling summary records the avoided materialization');
my $single_line = File::Spec->catfile($tmp, 'single-line.log');
write_file($single_line, 'SINGLE_HEAD' . ('X' x (20 * 1024)) . 'SINGLE_TAIL');
my $single_excerpt = read_stats_log_excerpt(
	$single_line, max_file_bytes => 8 * 1024,
);
like($single_excerpt, qr/^SINGLE_HEAD/,
	'oversized single-line excerpts retain their first edge');
like($single_excerpt, qr/SINGLE_TAIL$/,
	'oversized single-line excerpts retain their last edge with clamped defaults');

my $tail = read_stats_log_excerpt(
	$large, mode => 'tail', tail_lines => 1, edge_bytes => 4096,
);
is($tail, 'TAIL_MARKER', 'tail mode reads the requested final complete line');
unlike($tail, qr/HEAD_MARKER|MIDDLE_MARKER/,
	'tail mode does not materialize unrelated log regions');

my $compressed = File::Spec->catfile($tmp, 'large.log.gz');
my $compressed_text =
	"GZIP_HEAD\n" . ('C' x (32 * 1024)) . "\n"
	. "GZIP_MIDDLE\n" . ('D' x (32 * 1024)) . "\nGZIP_TAIL\n";
gzip(\$compressed_text => $compressed)
	or die "Cannot create $compressed: $GzipError";
reset_stats_log_sampling();
my $compressed_excerpt = read_stats_log_excerpt(
	$compressed, max_file_bytes => 8 * 1024, edge_bytes => 1024,
);
like($compressed_excerpt, qr/^GZIP_HEAD$/m,
	'oversized compressed-log excerpt retains its header');
like($compressed_excerpt, qr/^GZIP_TAIL$/m,
	'oversized compressed-log excerpt retains its tail');
unlike($compressed_excerpt, qr/GZIP_MIDDLE/,
	'oversized compressed-log excerpt retains bounded edge data only');
$summary = stats_log_sampling_summary();
is($summary->{compressed_files}, 1,
	'compressed oversized logs are reported as bounded-memory scans');

unlink $large or die "Cannot remove plain fallback fixture $large: $!";
my $resolved_tail = read_stats_log_excerpt(
	File::Spec->catfile($tmp, 'large.log'),
	mode => 'tail', tail_lines => 1, max_file_bytes => 8 * 1024, edge_bytes => 1024,
);
is($resolved_tail, 'GZIP_TAIL',
	'statistics reader resolves a compressed fallback when the requested plain log is absent');

open my $mataf4_fh, '<', 'MATAF4.pl' or die "Cannot read MATAF4.pl: $!";
local $/;
my $mataf4_stats = <$mataf4_fh>;
close $mataf4_fh or die "Cannot close MATAF4.pl: $!";
my ($postprocess_code) = $mataf4_stats =~ /(sub postprocess.*?)(?=\nsub spaceInAssGrp)/s;
ok(defined($postprocess_code), 'postprocessing can be isolated for performance checks');
like($postprocess_code,
	qr/read_sample_completion\(.*?\$closedSample->\{metagstats\}.*?my \$statsCollectionSeconds.*?my \$statsWriteStarted.*?_metag_stats_text.*?rename \$temporary, \$MGSfile.*?Created sample summary table/s,
	"summary timing covers sentinel reads plus table serialization and publication");
my ($completion_stats_code) = $mataf4_stats =~ /(sub createSampleCompletionSentinel.*?)(?=\nsub cleanupCompletionRequirements)/s;
ok(defined($completion_stats_code), "sample closure statistics can be isolated");
like($completion_stats_code,
	qr/reset_stats_log_sampling\(\).*?values\s*=>\s*smplStats\s*\(.*?write_sample_completion\(.*?stats_log_sampling_summary\(\).*?Statistics log safeguard/s,
	"sample closure caches statistics and reports oversized log sampling");
my ($sample_stats_code) = $mataf4_stats =~ /(sub smplStats .*?)(?=\n# smplStats is implemented)/s;
ok(defined($sample_stats_code), 'per-sample statistics implementation can be isolated');
my @sdm_log_globs = $sample_stats_code =~ /glob\("\$inD\/LOGandSUB\/sdm\/filter\*\.log"\)/g;
is(scalar(@sdm_log_globs), 1,
	'per-sample SDM statistics scan the log directory once');
unlike($sample_stats_code, qr/glob\("\$inD\/LOGandSUB\/sdm\/filterSuppl\*\.log"\)/,
	'support-read SDM logs are selected from the shared directory scan');
like($sample_stats_code,
	qr/my \$sdmHistogramMax = .*?_sdm_histogram_max_length.*?sdmStatsMany\(\\\@primary_logs, \$inD, '', \$sdmHistogramMax\).*?sdmStatsMany\(\\\@support_logs, \$inD, '_Sup', \$sdmHistogramMax\)/s,
	'primary and support SDM summaries reuse one histogram tail read');
like($mataf4_stats,
	qr/my \$alignStats = read_stats_log_excerpt\(\$inFi\);\s*if \(\$alignStats ne ''\)/,
	'mapping statistics accept a compressed fallback without rechecking the plain path');
like($mataf4_stats,
	qr/sub getGeneStats.*?while \(my \$line = <\$geneFH>\).*?close \$geneFH/s,
	'gene statistics are streamed instead of materialized and split');
like($mataf4_stats,
	qr/sub getASsemblyStats.*?while \(my \$header = <\$circFH>\).*?next unless \$header =~ \/\^>\//s,
	'circular-contig statistics stream FASTA headers and skip sequence bodies');

done_testing();
