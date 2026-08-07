use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::SNP qw(estimateConsensusCores);

is(estimateConsensusCores(0, 10), 1, 'unknown or empty consensus input uses one core');
is(estimateConsensusCores(299, 10), 1, 'alignment input below 300 MiB remains single-core');
is(estimateConsensusCores(300, 10), 2, 'alignment input at 300 MiB receives two cores');
is(estimateConsensusCores(600, 10), 3, 'alignment input at 600 MiB receives three cores');
is(estimateConsensusCores(1024, 10), 4, 'alignment input at 1 GiB receives four cores');
is(estimateConsensusCores(2048, 10), 5, 'alignment input at 2 GiB receives five cores');
is(estimateConsensusCores(4096, 10), 6, 'alignment input at 4 GiB receives six cores');
is(estimateConsensusCores(6144, 10), 7, 'alignment input at 6 GiB receives seven cores');
is(estimateConsensusCores(8192, 10), 8, 'alignment input at 8 GiB receives eight cores');
is(estimateConsensusCores(9216, 10), 9, 'alignment input at 9 GiB receives nine cores');
is(estimateConsensusCores(10240, 10), 10, 'alignment input at 10 GiB reaches ten cores');
is(estimateConsensusCores(10241, 8), 8, 'large consensus input is capped by SNPcores');

my $root = tempdir(CLEANUP => 1);
my $planner = File::Spec->catfile($Bin, '..', 'secScripts', 'SNP', 'plan_consensus_regions.pl');
my $fai = File::Spec->catfile($root, 'reference.fa.fai');
my $mapping = File::Spec->catfile($root, 'sample.bam');
my $depth = File::Spec->catfile($root, 'depth.percontig');
my $prefix = File::Spec->catfile($root, 'sample.');

open my $fh, '>', $fai or die $!;
print {$fh} "ctgA\t100\t0\t0\t0\nctgB\t100\t0\t0\t0\n";
close $fh;
open $fh, '>', $mapping or die $!;
print {$fh} "mapping placeholder\n";
close $fh;
open $fh, '>', $depth or die $!;
print {$fh} "ctgA\t20\nctgB\t0\n";
close $fh;

my $status = system(
	$^X, $planner,
	'--fai', $fai,
	'--mapping', $mapping,
	'--depth', $depth,
	'--jobs', 3,
	'--output-prefix', $prefix,
	'--samtools', 'unused-in-depth-mode',
);
is($status, 0, 'runtime SNP region planner succeeds with per-contig depth');

my %next_start = (ctgA => 0, ctgB => 0);
my @region_bases;
for my $index (0 .. 2) {
	my $bed = "$prefix$index.bed";
	ok(-s $bed, "runtime planner writes nonempty BED $index");
	open $fh, '<', $bed or die $!;
	my $bases = 0;
	while (my $line = <$fh>) {
		chomp $line;
		my ($contig, $start, $stop) = split /\t/, $line;
		is($start, $next_start{$contig}, "BED intervals remain contiguous for $contig");
		ok($stop > $start, "BED interval is nonempty for $contig");
		$next_start{$contig} = $stop;
		$bases += $stop - $start;
	}
	close $fh;
	push @region_bases, $bases;
}
is_deeply(\%next_start, {ctgA => 100, ctgB => 100},
	'runtime planning covers every reference base exactly once');
ok($region_bases[0] < $region_bases[2],
	'high-depth sequence receives fewer bases than low-depth sequence');

my $fake_samtools = File::Spec->catfile($root, 'samtools');
open $fh, '>', $fake_samtools or die $!;
print {$fh} <<'FAKE_SAMTOOLS';
#!/usr/bin/env perl
use strict;
use warnings;
die "expected idxstats\n" unless @ARGV == 2 && $ARGV[0] eq 'idxstats';
print "ctgA\t100\t1000\t0\nctgB\t100\t0\t0\n*\t0\t0\t0\n";
FAKE_SAMTOOLS
close $fh;
chmod 0755, $fake_samtools or die $!;
my $idx_prefix = File::Spec->catfile($root, 'support.');
$status = system(
	$^X, $planner,
	'--fai', $fai,
	'--mapping', $mapping,
	'--jobs', 2,
	'--output-prefix', $idx_prefix,
	'--samtools', $fake_samtools,
);
is($status, 0, 'runtime planner falls back to mapping index statistics');
ok(-s "${idx_prefix}0.bed" && -s "${idx_prefix}1.bed",
	'idxstats fallback writes every requested supplementary BED');

open my $snp_fh, '<', File::Spec->catfile($Bin, '..', 'Mods', 'SNP.pm') or die $!;
my $snp_source = do { local $/; <$snp_fh> };
close $snp_fh;
like($snp_source,
	qr/my \@snpScopes = \(.*?name => 'primary'.*?mapping => \$tar\[0\].*?name => 'supplementary'.*?mapping => \$tarS\[0\]/s,
	'primary and supplementary consensus products share one explicit scope model');
like($snp_source,
	qr/for my \$scope \(grep \{ \$_->\{run\} \} \@snpScopes\).*?\$regionPlanner.*?--mapping \$scope->\{mapping\}.*?--output-prefix \$outputPrefix/s,
	'all requested mappings use the same runtime region-planning path');
like($snp_source,
	qr/for my \$scope \(grep \{ \$_->\{run\} \} \@snpScopes\).*?pileupcall\(.*?\[\$scope->\{mapping\}\].*?\$scope->\{chunks\} = \$chunks/s,
	'all requested mappings use the same pileup and chunk-collection path');
like($snp_source,
	qr/for my \$scope \(grep \{ \$_->\{run\} \} \@snpScopes\).*?my \@chunks = \@\{\$scope->\{chunks\}\}.*?\$bcftBin concat.*?\$scope->\{vcf\}/s,
	'only scopes with planned chunks enter the shared VCF concatenation path');
like($snp_source,
	qr/sub _coverage_file_for_mapping.*?return "\$coverage\.gz" if \$allowPendingInputs/s,
	'pending consensus jobs use the canonical future compressed-coverage path');
like($snp_source,
	qr/pending SNP inputs require scheduler dependencies.*?test -s \$refFA.*?test -s \$_->\{mapping\}.*?test -s \$contigDepthF/s,
	'pending consensus inputs require an afterok chain and validate every requested scope inside the allocation');
like($snp_source,
	qr/SNP GFF is required for consensus statistics.*?unless \$allowPendingInputs \|\| -s \$gffF.*?test -s \$gffF.*?my \$vcf2fnaIns = "-ref \$refFA -gff \$gffF "/s,
	"vcf2fna always receives its GFF, including stats-only and pending-producer runs");
like($snp_source,
	qr/SNP ContigStats depth is missing or empty.*?test -s \$contigDepthF.*?--depth \$contigDepthF/s,
	"primary consensus validates ContigStats depth before runtime region planning");
like($snp_source,
	qr/my \@normalizeScopes = grep.*?for my \$scope \(\@normalizeScopes\).*?rm -f \$vcf \$vcf\.csi.*?touch \$normStone/s,
	"all requested scopes use the same normalization checkpoint path");
like($snp_source,
	qr/my \@consensusScopes = grep.*?-inVCF .*?map \{ \$_->\{vcf\} \} \@consensusScopes.*?-depthF .*?\$_->\{depth_file\}/s,
	"vcf2fna inputs and depth files are derived from the same requested scopes");
like($snp_source,
	qr/estimateConsensusCores\(\$consensusInputMB, \$maxSNPcores\).*?qsubSystem\([^;]+\$cmdFTag\.oSNPc\.sh[^;]+\$actualCores,/s,
	'oSNPc submission requests the automatically estimated runtime region count');
unlike($snp_source, qr/int\(\$actualCores\s*\*\s*1\.1\)/,
	'oSNPc no longer inflates or rounds its estimated core request');

done_testing;
