use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

open my $sourceFH, '<', File::Spec->catfile($Bin, '..', 'MATAF4.pl')
	or die "Cannot inspect MATAF4.pl: $!";
my $mataf4 = do { local $/; <$sourceFH> };
close $sourceFH;

my ($optionHelpers) = $mataf4 =~ /(sub _shell_quote\s*\{.*?)(?=^sub sdmOptSet)/ms;
ok(defined $optionHelpers, 'SDM option helpers can be isolated for testing');
eval "package TestSDMCleaner; our \%MFopt; $optionHelpers";
is($@, '', 'isolated SDM option helpers compile');

is(TestSDMCleaner::_shell_quote("reads/a'b.fq"), q{'reads/a'"'"'b.fq'},
	'shell quoting safely preserves an apostrophe');
eval { TestSDMCleaner::_shell_quote("bad\nargument") };
like($@, qr/NUL or newline/, 'shell arguments containing command separators are rejected');

my $tmpdir = tempdir(CLEANUP => 1);
my $baseOptions = File::Spec->catfile($tmpdir, 'base.txt');
open my $baseFH, '>', $baseOptions or die "Cannot write $baseOptions: $!";
print {$baseFH} <<'OPTIONS';
minSeqLength	50
maxAccumulatedError	1.0
maxAmbiguousNT	2
*maxAmbiguousNT	99
TrimWindowWidth	15
TrimWindowThreshhold	20
BinErrorModelAlpha	0.01
OPTIONS
close $baseFH;

$TestSDMCleaner::MFopt{sdmProbabilisticFilter} = 0;
$TestSDMCleaner::MFopt{sdm_opt} = {minSeqLength => 42};
my $adapted = TestSDMCleaner::adaptSDMopt($baseOptions, $tmpdir, 80, 'hiSeq');
open my $adaptedFH, '<', $adapted or die "Cannot read $adapted: $!";
my $adaptedText = do { local $/; <$adaptedFH> };
close $adaptedFH;
like($adaptedText, qr/^maxAccumulatedError\t1\.2$/m,
	'short-read accumulated-error threshold is adapted');
like($adaptedText, qr/^maxAmbiguousNT\t1$/m,
	'reads shorter than 90 bases retain the stricter ambiguous-base limit');
unlike($adaptedText, qr/^maxAmbiguousNT\t2$/m,
	'the short-read ambiguous-base limit is not overwritten');
like($adaptedText, qr/^\*maxAmbiguousNT\t99$/m,
	'commented option examples are not rewritten');
like($adaptedText, qr/^BinErrorModelAlpha\t-1$/m,
	'probabilistic filtering can be disabled deterministically');
like($adaptedText, qr/^minSeqLength\t42$/m,
	'explicit SDM overrides are applied');

my $longAdapted = TestSDMCleaner::adaptSDMopt($baseOptions, $tmpdir, 0, 'PB');
open my $longFH, '<', $longAdapted or die "Cannot read $longAdapted: $!";
my $longText = do { local $/; <$longFH> };
close $longFH;
like($longText, qr/^minSeqLength\t42$/m,
	'explicit overrides also apply when read length is unknown or reads are long');
like($longText, qr/^maxAccumulatedError\t1\.0$/m,
	'long-read options do not receive short-read length heuristics');

$TestSDMCleaner::MFopt{sdm_opt} = {missingOption => 1};
eval { TestSDMCleaner::adaptSDMopt($baseOptions, $tmpdir, 81, 'hiSeq') };
like($@, qr/Expected one active 'missingOption'/,
	'unknown SDM overrides fail instead of being silently ignored');
eval { TestSDMCleaner::_set_sdm_option("name\told\n", 'name', "bad\nvalue") };
like($@, qr/Invalid value/, 'SDM override values cannot inject additional option lines');

my ($optionSelector) = $mataf4 =~ /(sub sdmOptSet\s*\{.*?)(?=^sub sdmClean)/ms;
like($optionSelector, qr/technology eq '454'.*?adaptSDMopt\([^\n]+?'pair'\).*?adaptSDMopt\([^\n]+?'single'\)/s,
	'paired and singleton 454 option files use distinct output names');

my ($cleaner) = $mataf4 =~ /(sub sdmClean\(\)\{.*?)(?=^sub mocat_reorder)/ms;
like($cleaner, qr/return '' unless \@\{\$libraries\}/,
	'an empty scope does not repeat an unrelated upstream dependency');
like($cleaner, qr/filterSupplDone\.stone.*?filterDone\.stone|filterDone\.stone.*?filterSupplDone\.stone/s,
	'primary and support cleaning use independent checkpoints');
like($cleaner, qr/_shell_command\(\s*\$sdmBin,/s,
	'SDM commands quote executable, input, output, and option arguments');
like($cleaner, qr/local \$QSBoptHR->\{tmpSpace\} = 0/,
	'scheduler scratch settings are restored even if submission fails');
like($cleaner, qr/grep \{ !-e \$_ \} \@requiredOutputs/,
	'a successful checkpoint accepts valid empty filtered outputs');
unlike($cleaner, qr/requiredNonEmpty/,
	'redundant size-based completion state has been removed');

done_testing();
