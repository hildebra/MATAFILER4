use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::Binning ();
use Mods::geneCat qw(calculate_spearman_correlation correlation);
use Mods::GenoMetaAss qw(checkSeqTech readGFF renameFastaCnts renameFastqCnts systemW);
use Mods::IO_Tamoc_progs qw(mapperDBbuilt);
use Mods::math qw(meanArray quantileArray round);
use Mods::SNP ();
use Mods::Subm qw(qsubSystem);
use Mods::TamocFunc qw(bam2cram);

my $root = tempdir(CLEANUP => 1);
$root =~ s{\\}{/}g;

my $aviti_error = '';
eval { checkSeqTech('AVITI') };
$aviti_error = $@;
is($aviti_error, '', 'AVITI is accepted as a sequencing technology');

is(round(-1.5, 0), -2, 'round handles negative halves symmetrically');
is(quantileArray(1, 1, 2, 3), 3, 'legacy quantile clamps the upper endpoint');
ok(!defined meanArray([]), 'empty means are undefined instead of dividing by zero');
is(correlation([1, 2, 3], [2, 4, 6]), '1.0000', 'Pearson correlation uses every value and the correct means');
is(calculate_spearman_correlation([1, 1, 2], [1, 2, 3]), 0.866025403784439,
	'Spearman correlation correctly handles tied ranks');

my $fastq = "$root/reads.fq";
open my $fq, '>', $fastq or die $!;
print {$fq} "\@old1\nAC\n+\n!!\n\@old2\nGT\n+\n##\n";
close $fq;
renameFastqCnts($fastq, 'sample');
open $fq, '<', $fastq or die $!;
my $fastq_text = do { local $/; <$fq> };
close $fq;
like($fastq_text, qr/^\@sample_0\nAC/m, 'the first FASTQ record is renamed');
like($fastq_text, qr/^\@sample_1\nGT/m, 'subsequent FASTQ records are renamed');

my $fasta = "$root/seqs.fa";
my $rename_log = "$root/seqs.rename";
open my $fa, '>', $fasta or die $!;
print {$fa} ">one\nAA\nAA\n>two\nCC\n";
close $fa;
renameFastaCnts($fasta, 'contig', $rename_log);
open $fa, '<', $fasta or die $!;
my $fasta_text = do { local $/; <$fa> };
close $fa;
is($fasta_text, ">contig_0\nAA\nAA\n>contig_1\nCC\n", 'FASTA IDs count records rather than wrapped lines');

my $gff = "$root/features.gff";
open my $gf, '>', $gff or die $!;
print {$gf} "ctg\tsrc\tgene\t1\t2\t.\t+\t.\tID=1_7\nmalformed\nctg\tsrc\tgene\t3\t4\t.\t+\t.\tName=no_id\n";
close $gf;
my $gff_records;
{
	local $SIG{__WARN__} = sub { };
	$gff_records = readGFF($gff);
}
is_deeply([sort keys %$gff_records], ['ctg_7'], 'malformed GFF records cannot reuse a stale regex capture');

my $index = "$root/reference.fa.bw2";
for my $suffix (qw(1 2 3 4 rev.1 rev.2)) {
	open my $idx, '>', "$index.$suffix.bt2" or die $!;
	print {$idx} 'index';
	close $idx;
}
ok(mapperDBbuilt("$root/reference.fa", 1), 'a complete Bowtie2 index is accepted');
ok(mapperDBbuilt("$root/reference.fa", -1),
	'automatic short-read mapper completeness resolves to Bowtie2 only');
ok(mapperDBbuilt("$root/reference.fa", -2),
	'automatic strobealign completeness does not require an index');
unlink "$index.3.bt2" or die $!;
ok(!mapperDBbuilt("$root/reference.fa", 1), 'a partial Bowtie2 index is rejected');

my $bam = "$root/direct.bam";
open my $bf, '>', $bam or die $!;
print {$bf} 'bam-data';
close $bf;
my ($conversion, $bams) = Mods::Binning::createBams([$bam], $root, $root, 'sample', "$root/ref.fa", 1, 0, 0, 'bam');
is($conversion, '', 'a direct BAM needs no conversion');
is_deeply($bams, [$bam], 'a direct BAM is retained as binner input');

my $failure = '';
eval { systemW('exit 7') };
$failure = $@;
like($failure, qr/exit code 7/, 'systemW propagates the decoded child exit code');

my $script = "$root/failing.sh";
my $lock = "$root/failing.lock";
my $options = {
	rTag => 'TST', doSubmit => 1, doSync => 0,
	medQueue => 'compute', medTime => '', longQueue => 'compute', longTime => '',
	gpuQueue => 'compute', netQueue => 'compute', highMemQueue => 'compute', shortQueue => 'compute',
	useLongQueue => 0, useGPUQueue => 0, useNetQueue => 0, useShortQueue => 0, useHiMemQueue => 0,
	gpuCount => 0, tmpSpace => 0, tmpSpaceTag => '', excludeNodes => '', submissionConfig => '',
	constraint => [], LOCKfile => $lock, xtraNodeCmds => '', wcKeysForJob => '', qmode => 'bash',
};
my $submit_error = '';
{
	my $submission_output = '';
	open my $captured_stdout, '>', \$submission_output or die $!;
	local *STDOUT = $captured_stdout;
	eval { qsubSystem($script, 'exit 9', 1, '1G', 'fail', '', '', 1, [], $options) };
}
$submit_error = $@;
like($submit_error, qr/Job submission failed \(exit 9\)/, 'failed direct submissions are reported');
ok(!-e $lock, 'a failed submission cannot create its lock');

my $fai = "$root/reference.fa.fai";
open my $fai_fh, '>', $fai or die $!;
print {$fai_fh} "ctg\t100\t0\t0\t0\n";
close $fai_fh;
is_deeply([Mods::SNP::regionsFromFAI($fai)], ['ctg:1-100'], 'samtools regions are one-based');

my $pending_bam = '$SLURM_LOCAL_SCRATCH/MF4/sample/sample-smd.bam';
my ($cram_command, $pending_cram);
my $bam2cram_error = '';
{
	no warnings 'redefine';
	local *Mods::TamocFunc::getProgPaths = sub { return 'samtools' };
	eval {
		($cram_command, $pending_cram) = bam2cram(
			$pending_bam, "$root/reference.fa", 1, 1, "$root/map.cram.sto", 4
		);
	};
}
$bam2cram_error = $@;
is($bam2cram_error, '', 'bam2cram accepts a BAM produced later in scheduler-local scratch');
like($cram_command, qr/if \[ ! -s "\$SLURM_LOCAL_SCRATCH\/MF4\/sample\/sample-smd\.bam" \]/,
	'bam2cram defers its input check to the generated scheduler command');
like($cram_command, qr/samtools .*?-o \$SLURM_LOCAL_SCRATCH\/MF4\/sample\/sample-smd\.cram/s,
	'bam2cram retains scheduler variable expansion in its conversion command');
is($pending_cram, '$SLURM_LOCAL_SCRATCH/MF4/sample/sample-smd.cram',
	'bam2cram returns the pending CRAM path');

done_testing;
