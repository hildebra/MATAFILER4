use strict;
use warnings;

use File::Spec;
use FindBin qw($Bin);
use Test::More;

open my $source, '<', File::Spec->catfile($Bin, '..', 'MATAF4.pl')
	or die "Cannot inspect MATAF4.pl: $!";
my $mataf4 = do { local $/; <$source> };
close $source;

my ($rem_comma) = $mataf4 =~ /(sub remComma\(\$\)\s*\{.*?^\})/ms;
my ($parser) = $mataf4 =~ /(sub _sdm_version\s*\{.*?)(?=^sub sdmStats)/ms;
ok(defined($rem_comma) && defined($parser), 'SDM parser can be isolated for testing');
eval "package TestSDMStats; $rem_comma\n$parser";
is($@, '', 'isolated SDM parser compiles');

my $sdm340 = <<'SDM';
This is sdm (simple demultiplexer) 3.40 beta.

sdm run in No Map Mode.
sdm 3.40 beta
Processing summary:
  Reads processed:                                             1,613,103
  Rejected:                                                47,438 (2.9%)
  Accepted (high quality):                             1,565,319 (97.0%)
  Accepted (mid quality):                                     346 (0.0%)

Post-filter read characteristics:
Min/Avg/Max stats Pair 1
     - sequence Length : 250/6393.24/27154
     - Quality :   38/86.4912/93
     - Median sequence Length : 6170, Quality : 92
     - Accum. Error 3.49868
SDM

my $new = TestSDMStats::_parse_sdm_stats_text($sdm340, 0, '');
is($new->{SDMVersion}, '3.40', 'SDM 3.40 version is detected');
is($new->{totRds}, 1613103, '3.40 processed reads are parsed without separators');
is($new->{Rejected1}, 47438, '3.40 rejected reads ignore the percentage annotation');
is($new->{Accepted1}, 1565665, 'high- and mid-quality accepted reads are combined');
is($new->{Singl1}, 1565665, 'single-end accepted total is reported consistently');
is($new->{AvgSeqLen}, '6393.24', '3.40 average sequence length is parsed');
is($new->{MaxSeqLength}, '27154', '3.40 maximum sequence length is parsed');
is($new->{AvgSeqQual}, '86.4912', '3.40 average quality is parsed');
is($new->{accErr}, '3.49868', '3.40 accumulated error is parsed');

my $sdm340_paired = <<'SDM';
sdm run in No Map Mode. Using paired end sequencing files.
sdm 3.40 beta
Processing summary:                                               Read 1             Read 2
  Reads processed:                                            67,414,539          67,414,539
  Rejected:                                             4,714,512 (7.0%)    6,000,180 (8.9%)
  Accepted (high quality):                            62,700,027 (93.0%)  61,414,359 (91.1%)
  Accepted (mid quality):                                       0 (0.0%)            0 (0.0%)
  Reads trimmed at the end:                               761,191 (1.1%)    1,569,210 (2.3%)
  Recovered singleton reads:                                           0                   0

Post-filter read characteristics:                                 Read 1             Read 2
Min/Avg/Max stats Pair 1
     - sequence Length : 105/149.981/150
     - Quality :   31/36.121/37
     - Accum. Error 0.0795506
SDM
my $paired = TestSDMStats::_parse_sdm_stats_text($sdm340_paired, 0, '');
is($paired->{totRds}, 134829078, '3.40 paired processed totals include both reads');
is_deeply([$paired->{Rejected1}, $paired->{Rejected2}], [4714512, 6000180],
	'3.40 paired rejected counts are split by read');
is_deeply([$paired->{Accepted1}, $paired->{Accepted2}], [62700027, 61414359],
	'3.40 paired accepted counts are split by read');
is_deeply([$paired->{Singl1}, $paired->{Singl2}], [0, 0],
	'3.40 recovered singleton counts are split by read');
is($paired->{AvgSeqLen}, '149.981', '3.40 paired read characteristics are parsed');

my $legacy = <<'SDM';
Reads processed: 1,000
Rejected: 100
Accepted (High qual): 900
 - sequence Length : 50/125.5/250
 - Quality : 20/35.5/40
 - Accum. Error 0.25
SDM
my $old = TestSDMStats::_parse_sdm_stats_text($legacy, 300, '_Sup');
is($old->{totRds_Sup}, 1000, 'legacy single-end totals remain supported');
is($old->{Accepted1_Sup}, 900, 'legacy accepted reads remain supported');
is($old->{MaxSeqLength_Sup}, 300, 'histogram maximum still overrides the log maximum');
is($old->{SDMVersion_Sup}, '', 'unversioned legacy logs do not invent a version');

done_testing();
