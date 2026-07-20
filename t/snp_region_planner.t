use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

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
	qr/if \(\$runSupport\).*?\$cmdAll \.= \$xtra2\.\$pilecmd if \(!\$onlyNormalize\)/s,
	'existing supplementary VCFs cannot trigger concat without planned chunks');
like($snp_source,
	qr/consVCF_region_planner.*?--mapping \$tar\[0\].*?--output-prefix \$bedPrefix/s,
	'primary region planning is emitted into the SNP allocation');
like($snp_source,
	qr/consVCF_region_planner.*?--mapping \$tarS\[0\].*?--output-prefix \$\{bedPrefix\}sup-/s,
	'supplementary mappings receive independent runtime region planning');

done_testing;
