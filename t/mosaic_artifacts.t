use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

sub write_file {
	my ($path, $contents) = @_;
	open my $fh, '>', $path or die "Cannot write $path: $!";
	print {$fh} $contents or die "Cannot write $path: $!";
	close $fh or die "Cannot close $path: $!";
}

my $root = File::Spec->catdir($Bin, '..');
my $tmp = tempdir(CLEANUP => 1);
my $catalogue = File::Spec->catdir($tmp, 'GC');
my $annotation_dir = File::Spec->catdir(
	$catalogue, qw(Anno Func emapper),
);
make_path($annotation_dir);

write_file(
	File::Spec->catfile($catalogue, 'compl.incompl.95.fna'),
	">1\nACGTACGTACGT\n",
);
write_file(
	File::Spec->catfile($annotation_dir, 'eggNOGmapper_NOG.geneAss'),
	"1\tNOG1\n",
);
write_file(
	File::Spec->catfile($catalogue, 'Matrix.mat.gz'),
	"gene\tS1\n1\t1\n",
);

my $raw = File::Spec->catfile($tmp, 'SB.clusters');
my $core = "$raw.core";
write_file($raw,
	"Bin\tgene\tobservations\tmulticopy\tassignment\tmarker\n"
	."MGS.1\t1\t10\t0\t1\t1\n");
write_file($core, "MGS.1\t1\t10\t0\t1\t1\n");

my $supplied_paf = File::Spec->catfile($tmp, 'supplied.paf');
write_file($supplied_paf, '');
my $fake_rtk = File::Spec->catfile($tmp, 'fake-rtk2');
write_file($fake_rtk, <<'FAKE_RTK');
#!/usr/bin/env perl
use strict;
use warnings;
my $prefix = '';
for (my $i = 0; $i < @ARGV; $i++) {
	$prefix = $ARGV[$i + 1] if $ARGV[$i] eq '-o';
}
die "missing -o\n" unless length $prefix;
open my $report, '>', "$prefix.mosaic.tsv" or die $!;
print {$report} join("\t", qw(
	mgs group left right status reason identity query_coverage target_coverage
)), "\n";
close $report or die $!;
open my $summary, '>', "$prefix.mosaic.summary.tsv" or die $!;
print {$summary} "metric\tvalue\ncandidate_pairs\t0\nconfirmed_pairs\t0\n";
close $summary or die $!;
open my $concat, '>', "$prefix.concat.list" or die $!;
close $concat or die $!;
FAKE_RTK
chmod 0755, $fake_rtk or die "Cannot make $fake_rtk executable: $!";

my $output = File::Spec->catfile($tmp, 'mosaic.confirmed.tsv');
my @obsolete_suffixes = qw(
	.minimap2.paf .rtk.mosaic.tsv .rtk.mosaic.summary.tsv .rtk.concat.list
	.rejected.tsv .outgroups.tsv
);
write_file("$output$_", "legacy\n") for @obsolete_suffixes;
my $script = File::Spec->catfile(
	$root, 'secScripts', 'MGS', 'prepare_mosaic_loci.pl',
);
my $status = system(
	$^X, "-I$root", $script,
	'-GCd', $catalogue, '-MGS', $raw, '-coreMGS', $core,
	'-output', $output, '-paf', $supplied_paf, '-rtk', $fake_rtk,
	'-tmpD', $tmp,
);
is($status, 0, 'minimal Mosaic preparation succeeds with a supplied PAF');

for my $path ($output, "$output.candidates.tsv", "$output.summary.tsv") {
	ok(-s $path, "persistent Mosaic output is nonempty: $path");
}
ok(-e $supplied_paf, 'an explicitly supplied PAF is retained');
for my $suffix (@obsolete_suffixes) {
	ok(!-e "$output$suffix", "redundant Mosaic artifact is absent: $suffix");
}
opendir my $tmp_fh, $tmp or die "Cannot inspect $tmp: $!";
my @workspaces = grep { /^mosaic-loci-/ } readdir $tmp_fh;
closedir $tmp_fh;
is_deeply(\@workspaces, [], 'generated Mosaic workspace is removed after completion');

done_testing();
