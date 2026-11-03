use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::GTDBTaxonomy qw(merge_gtdb_summaries read_gtdb_taxonomy);

sub write_file {
	my ($path, $contents) = @_;
	open my $fh, '>', $path or die "Cannot write $path: $!\n";
	print {$fh} $contents;
	close $fh or die "Cannot close $path: $!\n";
}

sub slurp {
	my ($path) = @_;
	open my $fh, '<', $path or die "Cannot read $path: $!\n";
	local $/;
	return <$fh>;
}

my $tmp = tempdir(CLEANUP => 1);
my $header = "user_genome\tclassification\tnote\n";
my $bac = File::Spec->catfile($tmp, 'gtdbtk.bac120.summary.tsv');
my $arc = File::Spec->catfile($tmp, 'gtdbtk.ar53.summary.tsv');
write_file($bac, $header . "MGS.1\td__Bacteria;p__Firmicutes;s__Example\tbacterial\n");
write_file($arc, $header . "MGS.2\td__Archaea;p__Thermoproteota;s__Example\tarchaeal\n");

my $merged = File::Spec->catfile($tmp, 'gtdbtk.summary.tsv');
my $taxonomy = File::Spec->catfile($tmp, 'GTDBTK.tax');
is(merge_gtdb_summaries([$bac, $arc], $merged, $taxonomy), 2,
	'merging domain summaries reports the number of taxonomy rows');
is(slurp($merged), $header
	. "MGS.1\td__Bacteria;p__Firmicutes;s__Example\tbacterial\n"
	. "MGS.2\td__Archaea;p__Thermoproteota;s__Example\tarchaeal\n",
	'bacterial and archaeal summaries share one canonical header');
is(slurp($taxonomy), "user_genome\tclassification\n"
	. "MGS.1\tBacteria;Firmicutes;Example\n"
	. "MGS.2\tArchaea;Thermoproteota;Example\n",
	'two-column taxonomy strips only GTDB rank prefixes');

my $parsed = read_gtdb_taxonomy($taxonomy);
is_deeply([sort keys %{$parsed}], [qw(MGS.1 MGS.2)],
	'taxonomy reader excludes the canonical header');
is($parsed->{'MGS.1'}, 'Bacteria;Firmicutes;Example;?;?;?;?',
	'taxonomy reader pads incomplete lineages to seven ranks');

my $legacy = File::Spec->catfile($tmp, 'legacy.GTDBTK.tax');
write_file($legacy, "user_genome\tclassification\nMGS.1\tBacteria;Firmicutes\n"
	. "user_genome\tclassification\nMGS.2\tArchaea;Thermoproteota;\n");
my $legacy_parsed = read_gtdb_taxonomy($legacy);
is_deeply([sort keys %{$legacy_parsed}], [qw(MGS.1 MGS.2)],
	'historical duplicate domain headers are ignored rather than becoming an MGS');
is($legacy_parsed->{'MGS.2'}, 'Archaea;Thermoproteota;?;?;?;?;?',
	'empty GTDB rank fields are represented as missing values');

my $empty_domain = File::Spec->catfile($tmp, 'empty-domain.tsv');
write_file($empty_domain, $header);
my $merged_with_empty = File::Spec->catfile($tmp, 'merged-with-empty.tsv');
my $tax_with_empty = File::Spec->catfile($tmp, 'tax-with-empty.tsv');
is(merge_gtdb_summaries([$empty_domain, $bac], $merged_with_empty, $tax_with_empty), 1,
	'a header-only unused domain is accepted when another domain has classifications');

my $all_empty_summary = File::Spec->catfile($tmp, 'all-empty-summary.tsv');
my $all_empty_tax = File::Spec->catfile($tmp, 'all-empty-tax.tsv');
eval { merge_gtdb_summaries([$empty_domain], $all_empty_summary, $all_empty_tax) };
like($@, qr/contain no taxonomy data rows/,
	'header-only GTDB output cannot satisfy taxonomy completion');
ok(!-e $all_empty_summary && !-e $all_empty_tax,
	'a failed empty merge does not publish partial taxonomy outputs');

my $bad_header = File::Spec->catfile($tmp, 'bad-header.tsv');
write_file($bad_header, "genome\ttaxonomy\nMGS.3\tBacteria\n");
eval {
	merge_gtdb_summaries([$bad_header],
		File::Spec->catfile($tmp, 'bad-merged.tsv'),
		File::Spec->catfile($tmp, 'bad-tax.tsv'));
};
like($@, qr/Invalid GTDB-Tk summary header/,
	'malformed GTDB summary headers fail explicitly');

my $duplicate = File::Spec->catfile($tmp, 'duplicate.tsv');
write_file($duplicate, $header . "MGS.1\td__Archaea\tduplicate\n");
eval {
	merge_gtdb_summaries([$bac, $duplicate],
		File::Spec->catfile($tmp, 'duplicate-merged.tsv'),
		File::Spec->catfile($tmp, 'duplicate-tax.tsv'));
};
like($@, qr/Duplicate genome 'MGS\.1'/,
	'domain summary overlap cannot silently overwrite an assignment');

my $header_only_tax = File::Spec->catfile($tmp, 'header-only.tax');
write_file($header_only_tax, "user_genome\tclassification\n");
eval { read_gtdb_taxonomy($header_only_tax) };
like($@, qr/contains no data rows/,
	'header-only taxonomy input fails before marker abundance construction');

done_testing();
