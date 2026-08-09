use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use IPC::Open3;
use Symbol qw(gensym);
use Test::More;

my $root = File::Spec->catdir($Bin, '..');
my $scripts = File::Spec->catdir($root, 'secScripts', 'assemblies');
my $tmp = tempdir(CLEANUP => 1);

sub write_file {
	my ($path, $contents) = @_;
	open my $fh, '>', $path or die "Cannot write $path: $!";
	print {$fh} $contents;
	close $fh;
}

sub read_file {
	my ($path) = @_;
	open my $fh, '<', $path or die "Cannot read $path: $!";
	local $/;
	return <$fh>;
}

sub run_script {
	my ($script, @arguments) = @_;
	my $error = gensym;
	my $pid = open3(undef, my $stdout, $error, $^X, '-I'.$root,
		File::Spec->catfile($scripts, $script), @arguments);
	my $output = do { local $/; <$stdout> } // '';
	my $errors = do { local $/; <$error> } // '';
	waitpid($pid, 0);
	return ($? >> 8, $output, $errors);
}

my $scaffold = File::Spec->catfile($tmp, 'scaffold.fa');
write_file($scaffold, ">scaf\n".('A' x 10).('N' x 5).('C' x 10)."\n");
my ($status, $output, $errors) = run_script('assemblathon_stats.pl', '-n', 5, $scaffold);
is($status, 0, 'assemblathon statistics accept a custom N threshold');
like($output, qr/Percentage of assembly in scaffolded contigs\s+80\.0%/,
	'scaffolded percentage uses scaffold totals and excludes gap bases');
like($output, qr/Number of contigs\s+2\b/, 'custom N threshold splits the scaffold into two contigs');
like($output, qr/Total size of contigs\s+20\b/, 'contig totals exclude the configured gap');
like($output, qr/Average length of break \(>=5 Ns\).*\s5\b/,
	'gap statistics use the configured threshold');

my $even = File::Spec->catfile($tmp, 'even.fa');
write_file($even, ">a\n".('A' x 10)."\n>b\n".('C' x 20)."\n");
($status, $output, $errors) = run_script('assemblathon_stats.pl', $even);
is($status, 0, 'assemblathon statistics accept two contigs');
like($output, qr/Median scaffold size\s+15\b/, 'even-sized median averages both middle values');

my $split = File::Spec->catfile($tmp, 'split.fa');
write_file($split, ">long description\nABCDEFGHIJ\n");
($status, $output, $errors) = run_script('splitFNAbyLength.pl', $split, 4);
is($status, 0, 'long FASTA record is split successfully');
my $split_text = read_file($split);
my @split_ids = $split_text =~ /^>(\S+)/mg;
is_deeply(\@split_ids,
	[qw(long_part1_1-4 long_part2_5-8 long_part3_9-10)],
	'split fragments receive unique coordinate-bearing identifiers');
is($split_text, ">long_part1_1-4 description\nABCD\n>long_part2_5-8 description\nEFGH\n>long_part3_9-10 description\nIJ\n",
	'split fragments preserve all sequence bases');

my $lengths = File::Spec->catfile($tmp, 'lengths.fa');
write_file($lengths, ">equal\nAAAAA\n>long\nCCCCCC\n");
($status, $output, $errors) = run_script('sepReadLength.pl', 5, $lengths);
is($status, 0, 'length separation succeeds');
is(read_file($lengths), ">equal\nAAAAA\n", 'record equal to threshold remains in short input');
is(read_file("$lengths.long"), ">long\nCCCCCC\n", 'only truly longer record enters long output');

my $bad_fastq = File::Spec->catfile($tmp, 'bad.fastq');
write_file($bad_fastq, "\@r\nAAAA\nnot-plus\nIIII\n");
($status, $output, $errors) = run_script('sizeFilterFas.pl', $bad_fastq, 1, 0);
isnt($status, 0, 'invalid FASTQ separator is rejected');
like($errors, qr/Invalid FASTQ separator/, 'FASTQ validation explains the failure');

my $assembly = File::Spec->catfile($tmp, 'assembly.fa');
write_file($assembly, ">a\nAAAA\n");
my $fake_samtools = File::Spec->catfile($tmp, 'fake-samtools');
write_file($fake_samtools, <<'FAKE_SAMTOOLS');
#!/usr/bin/env perl
use strict;
use warnings;
if ($ARGV[0] eq 'faidx') {
	open my $fasta, '<', $ARGV[1] or die $!;
	open my $fai, '>', "$ARGV[1].fai" or die $!;
	my ($name, $length) = ('', 0);
	while (my $line = <$fasta>) {
		if ($line =~ /^>(\S+)/) {
			print {$fai} "$name\t$length\t0\t0\t0\n" if length $name;
			($name, $length) = ($1, 0);
		} else {
			$line =~ s/\s+//g;
			$length += length $line;
		}
	}
	print {$fai} "$name\t$length\t0\t0\t0\n" if length $name;
	exit 0;
}
if ($ARGV[0] eq 'view' && $ARGV[1] eq '-H') {
	open my $alignment, '<', $ARGV[2] or die $!;
	print while <$alignment>;
	exit 0;
}
die "unsupported fake samtools invocation: @ARGV\n";
FAKE_SAMTOOLS
chmod 0755, $fake_samtools or die $!;
my $mapping_dir = File::Spec->catdir($tmp, 'mapping');
mkdir $mapping_dir or die $!;
write_file(File::Spec->catfile($mapping_dir, 'done.sto'), "sample-smd.cram\n");
my $alignment = File::Spec->catfile($mapping_dir, 'sample-smd.cram');
write_file($alignment, "\@HD\tVN:1.6\n\@SQ\tSN:a\tLN:4\n");
($status, $output, $errors) = run_script('validate_mapping_references.pl',
	'--assembly', $assembly, '--samtools', $fake_samtools,
	'--sample-dirs', $tmp);
is($status, 0, 'binning reference validator accepts a matching CRAM dictionary');
like($output, qr/Validated 1 mapping file/, 'successful validation reports its alignment count');
write_file($alignment, "\@HD\tVN:1.6\n\@SQ\tSN:preassembly_contig\tLN:4\n");
($status, $output, $errors) = run_script('validate_mapping_references.pl',
	'--assembly', $assembly, '--samtools', $fake_samtools,
	'--sample-dirs', $tmp);
isnt($status, 0, 'binning reference validator rejects a stale preassembly CRAM');
like($errors, qr/mapping\/reference mismatch.*missing 'preassembly_contig'/is,
	'reference mismatch failure identifies the stale contig and required remap');

unlink $alignment or die $!;
write_file(File::Spec->catfile($mapping_dir, 'done.sto'), "sample.sup-smd.cram\n");
my $support_alignment = File::Spec->catfile($mapping_dir, 'sample.sup-smd.cram');
write_file($support_alignment, "\@HD\tVN:1.6\n\@SQ\tSN:a\tLN:4\n");
($status, $output, $errors) = run_script('validate_mapping_references.pl',
	'--assembly', $assembly, '--samtools', $fake_samtools,
	'--sample-dirs', $tmp);
is($status, 0, 'binning reference validator accepts a support-only mapping');
like($output, qr/Validated 1 mapping file/,
	'support-only mapping is validated without inventing a primary mapping');

my $bin_dir = File::Spec->catdir($tmp, 'bins');
mkdir $bin_dir or die $!;
my $sentinel = File::Spec->catfile($bin_dir, 'keep.txt');
write_file($sentinel, "keep\n");
($status, $output, $errors) = run_script('runBinners.pl',
	'-binner', 0, '-binD', $bin_dir, '-tmpD', File::Spec->catdir($tmp, 'scratch'),
	'-smplID', 'sample', '-assmbl', $assembly, '-assmblGrp', 1,
	'-cores', 1, '-smplDirs', $tmp);
isnt($status, 0, 'unsupported binner fails before execution');
ok(-e $sentinel, 'invalid binner does not erase an existing output directory');

my $small_bin_dir = File::Spec->catdir($tmp, 'small-bins');
my $small_bin_tmp = File::Spec->catdir($tmp, 'small-bin-tmp');
my $sample_completion = File::Spec->catfile($tmp, "MF4.sentinel.small.json");
write_file($sample_completion, "closed\n");
($status, $output, $errors) = run_script('runBinners.pl',
	'-binner', 1, '-binD', $small_bin_dir, '-tmpD', $small_bin_tmp,
	'-smplID', 'small', '-assmbl', $assembly, '-assmblGrp', 1,
	'-cores', 1, '-smplDirs', $tmp, '-minAssemblySizeMB', 0.000005);
is($status, 0, 'undersized assembly completes without invoking a binner');
ok(!-e $sample_completion,
	"binner execution invalidates the sample completion sentinel before publishing output");
like($output, qr/Assembly has 4 bp.*publishing an empty bin assignment/s,
	'undersized assembly reports the measured sequence size and cutoff action');
my $small_assignment = File::Spec->catfile($small_bin_dir, 'small');
ok(-e $small_assignment && !-s $small_assignment,
	'undersized assembly publishes the standard empty bin assignment');
ok(-s File::Spec->catfile($small_bin_dir, 'Binning.stone'),
	'undersized assembly publishes the normal binner completion stone');

my $scg_small_bin_dir = File::Spec->catdir($tmp, 'scg-small-bins');
my $scg_small_tmp = File::Spec->catdir($tmp, 'scg-small-tmp');
($status, $output, $errors) = run_script('runBinners.pl',
	'-binner', 5, '-binD', $scg_small_bin_dir, '-tmpD', $scg_small_tmp,
	'-smplID', 'scg-small', '-assmbl', $assembly, '-assmblGrp', 1,
	'-cores', 1, '-smplDirs', $tmp);
is($status, 0,
	'SCGBinner input with fewer than two eligible contigs completes as an empty result');
like($output,
	qr/SCGBinner preflight: 1 total contigs; 0 contigs >=1000 bp.*?requires at least 2 contigs >=1000 bp.*?publishing an empty bin assignment/s,
	'SCGBinner reports the exact eligible-contig reason before skipping training');
ok(-e File::Spec->catfile($scg_small_bin_dir, 'scg-small')
		&& !-s File::Spec->catfile($scg_small_bin_dir, 'scg-small'),
	'SCGBinner technical-minimum guard publishes the standard empty assignment');
ok(-s File::Spec->catfile($scg_small_bin_dir, 'Binning.stone'),
	'SCGBinner technical-minimum guard publishes the normal completion stone');

($status, $output, $errors) = run_script('checkBinQual.pl',
	'-asm', $assembly, '-binF', $small_assignment,
	'-tmpD', File::Spec->catdir($tmp, 'small-quality-tmp'),
	'-ncore', 1, '-checkM2', 1, '-checkM1', 0, '-binner', 1);
is($status, 0, 'standard quality postprocessing accepts the undersized pseudo output');
ok(-s "$small_assignment.cm2" && -s "$small_assignment.assStat",
	'undersized assembly receives complete empty quality and assembly-stat outputs');

done_testing();
