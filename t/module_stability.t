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

my $split_reference = "$root/split-reference.fa";
open my $split_fa, '>', $split_reference or die $!;
print {$split_fa} ">ctgA\nAAAAAAAAAA\n>ctgB\nCCCCC\n";
close $split_fa;
open my $split_fai, '>', "$split_reference.fai" or die $!;
print {$split_fai} "ctgA\t10\t0\t0\t0\nctgB\t5\t0\t0\t0\n";
close $split_fai;
my ($base_regions) = Mods::SNP::getRegionsBam(6, $split_reference, $root);
is_deeply($base_regions, [
	"ctgA\t0\t6\n",
	"ctgA\t6\t10\nctgB\t0\t2\n",
	"ctgB\t2\t5\n",
], 'FAI region splitting emits complete zero-based half-open BED coverage');

my $depth_file = "$root/depth.percontig";
open my $depth_fh, '>', $depth_file or die $!;
print {$depth_fh} "S__C1_L=1000=\t10000\nS__C2_L=1000=\t10000\n";
close $depth_fh;
my ($depth_regions) = Mods::SNP::getRegionsBamDepth($depth_file, 4, 4);
my %next_start = ("S__C1_L=1000=" => 0, "S__C2_L=1000=" => 0);
for my $line (map { split /\n/ } @{$depth_regions}) {
	next unless length $line;
	my ($contig, $start, $stop) = split /\t/, $line;
	is($start, $next_start{$contig}, "depth-balanced BED is contiguous for $contig");
	ok($stop > $start, "depth-balanced BED interval is nonempty for $contig");
	$next_start{$contig} = $stop;
}
is_deeply(\%next_start, {"S__C1_L=1000=" => 1000, "S__C2_L=1000=" => 1000},
	'depth-balanced splitting neither drops nor duplicates contig bases');

my $pileup_bed_dir = "$root/pileup/";
my ($chunk_files, $pileup_command);
{
	no warnings 'redefine';
	local *Mods::SNP::getProgPaths = sub { return $_[0] };
	my (undef, $chunks, $command) = Mods::SNP::pileupcall(
		["$root/input.cram"], '', {
			assembly => $split_reference, nodeTmpD => $root, smpl => 'sample',
			qsubDir => $pileup_bed_dir, runLocal => 1, JNUM => 1,
			SNPcaller => 'FB', overwrite => 0, deferRegionPlanning => 0,
			SeqTech => 'ILL', run2ctg => 1, rdep => '', normIndels => 1,
		}, {}, $root, "$root/chunk", 1, ["ctgA\t0\t10\n"], 1,
	);
	$chunk_files = $chunks;
	$pileup_command = $command;
}
is_deeply($chunk_files, ["$root/chunk.0.vcf.gz"], 'pileup declares its exact BGZF chunk output');
my $quoted_chunk_file = quotemeta("$root/chunk.0.vcf.gz");
like($pileup_command, qr/freebayes.*?\| bcftools view -Oz -o $quoted_chunk_file -/s,
	'FreeBayes output is converted to BGZF instead of receiving a misleading .gz suffix');
like($pileup_command, qr/test -s $quoted_chunk_file && bcftools index -f $quoted_chunk_file && test -s $quoted_chunk_file\.csi/s,
	'each FreeBayes VCF chunk is indexed before concat');
unlike($pileup_command, qr/input\.cram\.(?:crai|bai)/,
	'pileup does not delete canonical mapping indexes');

my $restart_chunk = "$root/restart-chunk.0.vcf.gz";
open my $restart_fh, '>', $restart_chunk or die $!;
print {$restart_fh} "existing BGZF placeholder\n";
close $restart_fh;
my $restart_command;
{
	no warnings 'redefine';
	local *Mods::SNP::getProgPaths = sub { return $_[0] };
	my (undef, undef, $command) = Mods::SNP::pileupcall(
		["$root/input.cram"], '', {
			assembly => $split_reference, nodeTmpD => $root, smpl => 'sample',
			qsubDir => $pileup_bed_dir, runLocal => 1, JNUM => 1,
			SNPcaller => 'FB', overwrite => 0, deferRegionPlanning => 0,
			SeqTech => 'ILL', run2ctg => 1, rdep => '', normIndels => 1,
		}, {}, $root, "$root/restart-chunk", 1, ["ctgA\t0\t10\n"], 1,
	);
	$restart_command = $command;
}
my $quoted_restart_chunk = quotemeta($restart_chunk);
like($restart_command, qr/bcftools index -f $quoted_restart_chunk && test -s $quoted_restart_chunk\.csi/s,
	'a restart repairs an existing unindexed VCF chunk');
unlike($restart_command, qr/freebayes -f/,
	'a restart does not repeat an already completed variant call solely to create its index');

my $mpileup_command;
{
	no warnings 'redefine';
	local *Mods::SNP::getProgPaths = sub { return $_[0] };
	my (undef, undef, $command) = Mods::SNP::pileupcall(
		["$root/input.cram"], '', {
			assembly => $split_reference, nodeTmpD => $root, smpl => 'sample',
			qsubDir => $pileup_bed_dir, runLocal => 1, JNUM => 1,
			SNPcaller => 'MPI', overwrite => 0, deferRegionPlanning => 0,
			SeqTech => 'ILL', run2ctg => 1, rdep => '', normIndels => 1,
		}, {}, $root, "$root/mpileup-chunk", 1, ["ctgA\t0\t10\n"], 1,
	);
	$mpileup_command = $command;
}
like($mpileup_command, qr/bcftools mpileup .*? -Ou .*? -a FORMAT\/DP,FORMAT\/AD,FORMAT\/ADF,FORMAT\/ADR,FORMAT\/SP/s,
	'bcftools calling uses a binary pipe and supported optional annotations');
my $quoted_mpileup_chunk = quotemeta("$root/mpileup-chunk.0.vcf.gz");
like($mpileup_command, qr/test -s $quoted_mpileup_chunk && bcftools index -f $quoted_mpileup_chunk && test -s $quoted_mpileup_chunk\.csi/s,
	'each mpileup VCF chunk is indexed before concat');
unlike($mpileup_command, qr/INFO\/(?:PV4|FS|IDV|MQ0F|BQBZ|SCBZ|RPBZ|MQBZ)/,
	'bcftools command does not pass obsolete or automatic INFO fields as requested annotations');

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
