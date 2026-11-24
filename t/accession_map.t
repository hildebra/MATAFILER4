use strict;
use warnings;

# makeAccessionMap.sh: one row per BioSample, amplicon / PCR / empty runs excluded.

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

my $script = File::Spec->catfile($Bin, '..', 'secScripts', 'fileManage', 'makeAccessionMap.sh');
my $tmp = tempdir(CLEANUP => 1);
my $report = "$tmp/report.tsv";
open my $fh, '>', $report or die $!;
print {$fh} join("\t", qw(run_accession sample_accession library_strategy library_selection read_count)), "\n";
print {$fh} "$_\n" for (
	"SRR1\tSAMN1\tAMPLICON\tPCR\t100",
	"SRR2\tSAMN1\tWGS\tRANDOM\t500",
	"SRR3\tSAMN2\tWGS\tRANDOM\t400",
	"SRR4\tSAMN2\tOTHER\tRANDOM\t300",
	"SRR5\tSAMN3\tWGS\tPCR\t50",
	"SRR6\tSAMN4\tWGS\tRANDOM\t0",
	"SRR7\tSAMN5\tWGA\tRANDOM\t10",
);
close $fh;

my $out = `bash '$script' '$report' RUN1 /out/ 2>'$tmp/err'`;
is($? >> 8, 0, 'script succeeds');
my @lines = split /\n/, $out;
is_deeply([@lines[0..2]], ["#SmplID\tENAdownload\tAssmblGrps", "#RunID\tRUN1", "#OutPath\t/out/"], 'map header');
is_deeply([@lines[3..$#lines]], ["SAMN1\tSRR2\tSAMN1", "SAMN2\tSRR3,SRR4\tSAMN2", "SAMN5\tSRR7\tSAMN5"],
	'one row per BioSample; runs of one sample merged; OTHER and WGA shotgun kept');
my $err = do { local (@ARGV, $/) = "$tmp/err"; <> };
like($err, qr/PCR.*SRR5/, 'mislabelled amplicon (WGS + PCR) excluded and reported');
like($err, qr/0 reads.*SRR6/, 'empty run excluded and reported');
like($err, qr/kept 4 run\(s\) in 3 sample\(s\); excluded AMPLICON 1, PCR 1, empty 1/, 'summary');

open $fh, '>', "$tmp/bad.tsv" or die $!; print {$fh} "run_accession\tlibrary_strategy\nSRR1\tWGS\n"; close $fh;
`bash '$script' '$tmp/bad.tsv' R /o/ 2>/dev/null`;
isnt($? >> 8, 0, 'missing columns are refused');

done_testing();
