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

my $separate_contigs = read_file(File::Spec->catfile($scripts, 'separateContigs.pl'));
like($separate_contigs,
	qr/my \$anyCoverageAvailable\s*=\s*.*?fileGZe\(\$primaryCoverage\).*?fileGZe\(\$supportCoverage\).*?unless \$anyCoverageAvailable/s,
	'contig statistics accepts either primary or supplementary mapped-read coverage');
unlike($separate_contigs,
	qr/die "Could not find required coverage file \$inF/,
	'missing primary coverage is not unconditionally fatal');
like($separate_contigs,
	qr/sub geneAbundance.*?contig_stats_coverage_complete\(\$outDab, \$oPrefix\).*?if \(!fileGZe\(\$inF\)\)/s,
	'completed coverage derivatives are recognized before requiring their source coverage file');
ok(index($separate_contigs, '$inD =~ s{[\\\\/]+$}{};') >= 0
		&& index($separate_contigs, '$inD .= "/";') >= 0,
	'input directories are normalized instead of requiring a caller-supplied trailing slash');
unlike($separate_contigs,
	qr/if \(-s \$outFfin.*?Gene abundance was already calculated/s,
	'incomplete uncompressed output subsets cannot bypass the shared completion contract');
unlike($separate_contigs,
	qr/\$readLength > 0 && \$readLengthSup > 0/,
	'a missing primary stream does not require an otherwise unused primary read length');
like($separate_contigs,
	qr/my \$kind = \$isSupport.*?read length must be a positive integer.*?unless \$readL > 0/s,
	'each available coverage stream validates only its own read length');

done_testing();
