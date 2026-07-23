use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::geneCat qw(readGene2tax);
use Mods::MGSLocus qw(build_locus_groups choose_locus_candidate robust_depth_mask);

sub write_file {
	my ($path, $contents) = @_;
	open my $fh, '>', $path or die "Cannot write $path: $!";
	print {$fh} $contents;
	close $fh or die "Cannot close $path: $!";
}

my $tmp = tempdir(CLEANUP => 1);
my $mapping = File::Spec->catfile($tmp, 'genes.gene2MGS');
write_file($mapping, "10\tMGS.1\tCOG1\n11\tMGS.1\tCOG1\n12\tMGS.1\tCOG2\n");
my ($si, $gene_to_cog, $gene_to_mgs, $priorities) = readGene2tax($mapping, 10, []);
is_deeply(
	$priorities->{'MGS.1'},
	['MGS.1|COG1|10', 'MGS.1|COG1|11', 'MGS.1|COG2|12'],
	'same-COG catalogue genes are retained in ranked order as distinct seed loci',
);
is($si->{'MGS.1'}{'MGS.1|COG1|11'}, 11, 'seed locus maps to its catalogue cluster');

my $duplicates = File::Spec->catfile($tmp, 'duplicate-genes.gene2MGS');
write_file($duplicates, join('',
	"20\tMGS.2\tCOG1\n",
	map { "20\tMGS.2\tCOG$_\n" } 2 .. 10,
));
my $duplicate_warnings = '';
{
	local $SIG{__WARN__} = sub { $duplicate_warnings .= join('', @_); };
	readGene2tax($duplicates, 20, []);
}
my @duplicate_examples = $duplicate_warnings =~ /Ignoring duplicate catalogue gene 20/g;
is(scalar(@duplicate_examples), 5,
	'duplicate catalogue gene diagnostics are limited to five detailed examples');
like($duplicate_warnings, qr/Suppressed 4 additional duplicate catalogue gene warnings.*9 total/,
	'duplicate catalogue gene diagnostics retain an aggregate suppressed count');

my @records = (
	{ mgs => 'MGS.1', cog => 'COG1', gene => '10', rank => 0 },
	{ mgs => 'MGS.1', cog => 'COG1', gene => '11', rank => 1 },
	{ mgs => 'MGS.1', cog => 'COG1', gene => '12', rank => 2 },
);
my %members = (
	10 => 'S1__ctg_1,S2__ctg_1',
	11 => 'S3__ctg_1,S4__ctg_1',
	12 => 'S1__ctg_5,S3__ctg_5',
);
my $protein = 'M' . ('ACDEFGHIKLMNPQRSTVWY' x 8);
my %proteins = (10 => $protein, 11 => $protein, 12 => $protein);
my $model = build_locus_groups(\@records, \%members, \%proteins);
is(scalar(@{$model->{groups}}), 2, 'mutually exclusive similar seeds merge but a co-occurring paralog stays separate');
is($model->{gene_to_locus}{10}, 'MGS.1|COG1|10', 'first ranked seed names the merged locus');
is($model->{gene_to_locus}{11}, 'MGS.1|COG1|10', 'alternative seed maps to the primary locus');
is($model->{gene_to_locus}{12}, 'MGS.1|COG1|12', 'co-occurring seed retains its own locus');

my $dominant = choose_locus_candidate([
	{ id => 'copyA', protein => $protein, depth => 12, seed => '10', context => {} },
	{ id => 'copyB', protein => $protein, depth => 3,  seed => '10', context => {} },
], { 10 => $protein }, {});
is($dominant->{status}, 'selected', 'a dominant duplicate can be resolved');
is($dominant->{candidate}{id}, 'copyA', 'dominant-depth candidate is selected within one seed');

my $ambiguous = choose_locus_candidate([
	{ id => 'copyA', protein => $protein, depth => 10, seed => '10', context => {} },
	{ id => 'copyB', protein => $protein, depth => 10, seed => '12', context => {} },
], { 10 => $protein, 12 => $protein }, {});
is($ambiguous->{status}, 'ambiguous', 'equivalent candidates from different seeds remain unresolved');

is_deeply(
	robust_depth_mask([10, 10, 11, 9, 10, 11, 9, 100]),
	[1, 1, 1, 1, 1, 1, 1, 0],
	'robust depth filtering removes a supported extreme outlier',
);
is_deeply(
	robust_depth_mask([10, 10, 100]),
	[1, 1, 1],
	'small locus sets are not quantile-filtered',
);

done_testing();
