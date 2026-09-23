use strict;
use warnings;

use File::Spec;
use FindBin qw($Bin);
use lib File::Spec->catdir($Bin, '..');
use Test::More;

use Mods::math qw(avgArray medianArray quantileArray quantileArrayR nonZero
	meanArray roundAr roundF);

# Numeric helpers shared by abundance, strain and QC code. These tests pin the
# statistical conventions, because the three quantile flavours in the code base
# intentionally differ (see docs/audits/2026-09-11/report.md).

is(nonZero([0, 1, -2, 3.5, 0]), 2, 'nonZero counts strictly positive values only');
is(nonZero([]), 0, 'nonZero of an empty array is zero');

is_deeply(avgArray([[1, 2, 3]]), [1, 2, 3], 'avgArray of one array returns that array');
is_deeply(avgArray([[1, 2, 3], [3, 4, 5]]), [2, 3, 4], 'avgArray averages element-wise');
my $original = [2, 4];
avgArray([$original, [4, 8]]);
is_deeply($original, [2, 4], 'avgArray does not modify its first input array');
ok(!eval { avgArray([[1, 2], [1]]); 1 }, 'avgArray rejects arrays of unequal length');
ok(!eval { avgArray([]); 1 }, 'avgArray rejects an empty list of arrays');

is(meanArray([1, 2, 3, 4]), 2.5, 'meanArray computes the arithmetic mean');
ok(!defined(meanArray([])), 'meanArray of an empty array is undefined');

is(medianArray(3, 1, 2), 2, 'odd-length median is the middle value after sorting');
is(medianArray(4, 1, 3, 2), 2.5, 'even-length median averages the two middle values');
is(medianArray(7), 7, 'singleton median is the value itself');
is(medianArray(10, 9, 100), 10, 'median sorts numerically, not lexically');

is(quantileArray(0.5, 1 .. 10), 6, 'quantileArray uses floor(n*q) as zero-based index');
is(quantileArray(0, 5, 1, 3), 1, 'quantileArray 0 is the minimum');
is(quantileArray(1, 5, 1, 3), 5, 'quantileArray 1 is clamped to the maximum');
ok(!defined(quantileArray(0.5)), 'quantileArray of no values is undefined');
ok(!eval { quantileArray(1.5, 1, 2); 1 }, 'quantileArray rejects fractions above 1');

is(quantileArrayR([1, 2, 3, 4, 5], 0.5), 3, 'interpolated median of 1..5 is 3');
is(quantileArrayR([1, 2, 3, 4], 0.5), 2.5, 'interpolated median interpolates between neighbours');
is_deeply([quantileArrayR([10, 20, 30, 40, 50], 0.1, 0.9)], [14, 46],
	'interpolated quantiles follow R type 7 (linear) positions');
is(quantileArrayR([1, 2, 3], 1), 3, 'interpolated quantile 1 is the maximum');
ok(!eval { quantileArrayR([], 0.5); 1 }, 'interpolated quantile rejects an empty array');
ok(!eval { quantileArrayR([1], -0.1); 1 }, 'interpolated quantile rejects negative thresholds');

is_deeply(roundAr([1.234, -1.235, 2.5], 1), [1.2, -1.2, 2.5], 'roundAr rounds to the requested decimals');
is_deeply(roundAr([0.5, -0.5], 0), [1, -1], 'roundAr rounds halves away from zero');
is(roundF(2.345, 100), 2.35, 'roundF rounds with an explicit scaling factor');
is(roundF(-2.5, 1), -3, 'roundF rounds negative halves away from zero');

done_testing();
