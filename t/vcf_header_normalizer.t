use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Test::More;

sub write_file {
	my ($path, $contents) = @_;
	open my $fh, '>', $path or die "Cannot write $path: $!";
	print {$fh} $contents or die "Cannot write $path: $!";
	close $fh or die "Cannot close $path: $!";
}

sub slurp {
	my ($path) = @_;
	my $contents = '';
	gunzip $path => \$contents
		or die "Cannot decompress $path: $GunzipError";
	return $contents;
}

my $tmp = tempdir(CLEANUP => 1);
my $script = File::Spec->catfile(
	$Bin, '..', 'secScripts', 'SNP', 'normalizeVCFHeaders.pl',
);
my $header = "#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO\n";
my $input = File::Spec->catfile($tmp, 'duplicated.vcf');
my $output = File::Spec->catfile($tmp, 'normalized.vcf.gz');
write_file($input,
	"##fileformat=VCFv4.2\n".$header
	."ctg1\t1\t.\tA\tC\t30\tPASS\t.\n"
	."##fileformat=VCFv4.2\n".$header
	."ctg2\t2\t.\tG\tT\t30\tPASS\t.\n"
);

is(system($^X, $script, '-input', $input, '-output', $output), 0,
	'identical concatenated VCF headers are normalized successfully');
my $normalized = slurp($output);
is(scalar(() = $normalized =~ /^#CHROM/mg), 1,
	'normalized VCF contains one column header');
like($normalized, qr/^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO$/m,
	'normalized VCF publishes a canonical tab-delimited column header');
like($normalized, qr/^ctg1\t1/m, 'records before the repeated header are retained');
like($normalized, qr/^ctg2\t2/m, 'records after the repeated header are retained');

my $incompatible = File::Spec->catfile($tmp, 'incompatible.vcf');
my $rejected = File::Spec->catfile($tmp, 'rejected.vcf.gz');
write_file($incompatible,
	$header."ctg1\t1\t.\tA\tC\t30\tPASS\t.\n"
	."#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tSAMPLE\n"
);
ok(system($^X, $script, '-input', $incompatible, '-output', $rejected) != 0,
	'incompatible repeated VCF headers are rejected');
ok(!-e $rejected, 'an incompatible VCF does not publish partial output');

done_testing();
