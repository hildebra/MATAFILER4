use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::StrainQC qw(breakpoint_gene_mask abundance_pattern_mask);

sub write_file {
	my ($path, $text) = @_;
	open my $fh, '>', $path or die $!;
	print {$fh} $text;
	close $fh or die $!;
}

my $tmp = tempdir(CLEANUP => 1);
my $breakpoints = File::Spec->catfile($tmp, 'sample.breakpoints.tsv');
write_file($breakpoints, join('',
	"contig\tstart\tend\tlength\tmean_depth\tleft_depth\tright_depth\n",
	"ctgA\t1000\t1100\t100\t0\t12\t11\n",
));
my $gff = File::Spec->catfile($tmp, 'genes.gff');
write_file($gff, join('',
	"##gff-version 3\n",
	"ctgA\tcaller\tCDS\t850\t950\t.\t+\t0\tID=7_41;\n",
	"ctgA\tcaller\tCDS\t990\t1020\t.\t+\t0\tID=7_42;\n",
	"ctgA\tcaller\tCDS\t1200\t1300\t.\t+\t0\tID=7_43;\n",
));
my %wanted = map { $_ => 1 } qw(S1__ctgA_41 S1__ctgA_42 S1__ctgA_43);
my $masked = breakpoint_gene_mask($gff, $breakpoints, \%wanted, 25);
ok(!$masked->{'S1__ctgA_41'}, 'gene beyond the breakpoint flank remains usable');
ok($masked->{'S1__ctgA_42'}, 'gene overlapping a mapping breakpoint is masked');
ok(!$masked->{'S1__ctgA_43'}, 'downstream non-overlapping gene remains usable');

is_deeply(
	abundance_pattern_mask([10, 9, 11, 10, 10, 9, 11, 80]),
	[1, 1, 1, 1, 1, 1, 1, 0],
	'a locus with a different abundance pattern is masked when enough loci support the baseline',
);
is_deeply(
	abundance_pattern_mask([10, 10, 80], {minimum_count => 8}),
	[1, 1, 1],
	'sparse strain observations are retained when abundance outlier inference is underpowered',
);

done_testing();
