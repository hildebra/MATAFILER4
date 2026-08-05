#!/usr/bin/env perl

use strict;
use warnings;

use File::Copy qw(copy);
use File::Path qw(make_path remove_tree);
use File::Spec;
use Getopt::Long qw(GetOptions);
use IO::Compress::Gzip qw(gzip $GzipError);
use Mods::GenoMetaAss qw(systemW);
use Mods::IO_Tamoc_progs qw(getProgPaths setConfigFile);

my $VERSION = '0.6';
my ($read1, $read2, $readS) = ('', '', '');
my ($alignPath, $tmpRoot, $sample, $configFile) = ('', '', '', '');
my $threads = 1;
my $doRiboAssembly = 0;

GetOptions(
	'R1=s'          => \$read1,
	'R2=s'          => \$read2,
	'RS=s'          => \$readS,
	'alignDir=s'    => \$alignPath,
	'smplID=s'      => \$sample,
	'tmpDir=s'      => \$tmpRoot,
	'cores=i'       => \$threads,
	'config=s'      => \$configFile,
	'assmblRibos=i' => \$doRiboAssembly,
) or die "Error in command line arguments\n";
die "Unexpected positional arguments: @ARGV\n" if @ARGV;
die "-alignDir is required\n" if $alignPath eq '';
die "-tmpDir is required\n" if $tmpRoot eq '';
die "-smplID is required\n" if $sample eq '';
die "-cores must be a positive integer\n" if $threads < 1;
die "Ribosomal assembly is no longer supported\n" if $doRiboAssembly;

my @r1 = parseReadList($read1);
my @r2 = parseReadList($read2);
my @single = parseReadList($readS);
die "At least one -R1 or -RS input is required\n" unless @r1 || @single;
die "-R2 was supplied without -R1\n" if @r2 && !@r1;
die "-R1 and -R2 must contain the same number of files\n"
	if @r2 && @r1 != @r2;
for my $readFile (@r1, @r2, @single) {
	die "Read input does not exist: $readFile\n" unless -e $readFile;
}

my (@pairedR1, @pairedR2);
if (@r2) {
	@pairedR1 = @r1;
	@pairedR2 = @r2;
} else {
	push @single, @r1;
}
my $hasPairs = @pairedR1 ? 1 : 0;
my $hasSingles = @single ? 1 : 0;

$alignPath = File::Spec->canonpath(File::Spec->rel2abs($alignPath));
$tmpRoot = File::Spec->canonpath(File::Spec->rel2abs($tmpRoot));
make_path($alignPath) unless -d $alignPath;
make_path($tmpRoot) unless -d $tmpRoot;
die "-alignDir is not a directory: $alignPath\n" unless -d $alignPath;
die "-tmpDir is not a directory: $tmpRoot\n" unless -d $tmpRoot;

my $safeSample = $sample;
$safeSample =~ s/[^A-Za-z0-9_.-]+/_/g;
my $tmpPath = File::Spec->catdir(
	$tmpRoot, "catchLSUSSU_${safeSample}_$$",
);
make_path($tmpPath);

for my $tag (qw(SSU LSU)) {
	repairInvalidCheckpoint(
		$alignPath, $tag, $hasPairs, $hasSingles,
	);
}
if (allMarkersComplete($alignPath, $hasPairs, $hasSingles)) {
	print "All RiboFind SortMeRNA targets are complete\n";
	remove_tree($tmpPath);
	exit 0;
}

setConfigFile($configFile);
my $sortmerna = getProgPaths('sortmerna');
announceSortmerna($sortmerna);

print "Skipping ITS\n";
for my $tag (qw(SSU LSU)) {
	next if markerComplete($alignPath, $tag, $hasPairs, $hasSingles);

	my $referenceKey = $tag eq 'SSU' ? 'SSUdbFAsrt' : 'LSUdbFAsrt';
	my $indexKey = $tag eq 'SSU' ? 'SSUidx' : 'LSUidx';
	my $reference = getProgPaths($referenceKey);
	my $index = getProgPaths($indexKey, 0);
	validateSortmernaReference($referenceKey, $reference);
	validateSortmernaIndex($indexKey, $index);

	invalidateMarkerState($alignPath, $tag);
	my %outputs = runMarker(
		sortmerna => $sortmerna,
		tag => $tag,
		reference => $reference,
		index => $index,
		paired_r1 => \@pairedR1,
		paired_r2 => \@pairedR2,
		single => \@single,
		work_dir => $tmpPath,
		threads => $threads,
	);
	publishMarkerOutputs($tmpPath, $alignPath, $tag, \%outputs);
	touchFile(File::Spec->catfile($alignPath, "${tag}_pull.sto"));
	die "$tag checkpoint could not be validated after publication\n"
		unless markerComplete($alignPath, $tag, $hasPairs, $hasSingles);
}

remove_tree($tmpPath) if -d $tmpPath;
print "Finished pulling out ribosomal reads\n";
exit 0;


sub parseReadList {
	my ($value) = @_;
	return grep { defined($_) && $_ ne '' && $_ ne '-1' }
		split(/,/, $value // '');
}


sub allMarkersComplete {
	my ($directory, $hasPaired, $hasSingle) = @_;
	for my $tag (qw(SSU LSU)) {
		return 0 unless markerComplete(
			$directory, $tag, $hasPaired, $hasSingle,
		);
	}
	return 1;
}


sub markerComplete {
	my ($directory, $tag, $hasPaired, $hasSingle) = @_;
	return 0 unless -e File::Spec->catfile(
		$directory, "${tag}_pull.sto",
	);
	if ($hasPaired) {
		return 0 unless -s File::Spec->catfile(
			$directory, "reads_${tag}.r1.fq.gz",
		);
		return 0 unless -s File::Spec->catfile(
			$directory, "reads_${tag}.r2.fq.gz",
		);
	}
	if ($hasSingle) {
		return 0 unless -s File::Spec->catfile(
			$directory, "reads_${tag}.fq.gz",
		);
	}
	return 1;
}


sub repairInvalidCheckpoint {
	my ($directory, $tag, $hasPaired, $hasSingle) = @_;
	my $stone = File::Spec->catfile($directory, "${tag}_pull.sto");
	return unless -e $stone;
	return if markerComplete($directory, $tag, $hasPaired, $hasSingle);
	warn "$tag profile checkpoint has missing or empty outputs; rerunning marker\n";
	invalidateMarkerState($directory, $tag);
}


sub invalidateMarkerState {
	my ($directory, $tag) = @_;
	my $lcaDir = File::Spec->catdir($directory, 'ltsLCA');
	for my $path (
		File::Spec->catfile($directory, "${tag}_pull.sto"),
		File::Spec->catfile($lcaDir, "${tag}_ass.sto"),
		File::Spec->catfile($lcaDir, 'Assigned.sto'),
		File::Spec->catfile($lcaDir, "${tag}riboRun_bl.hiera.txt"),
		File::Spec->catfile($lcaDir, "${tag}riboRun_bl.hiera.txt.gz"),
	) {
		unlinkChecked($path);
	}
}


sub runMarker {
	my (%args) = @_;
	my %outputs = (r1 => [], r2 => [], single => []);
	my $pairR1 = $args{paired_r1};
	my $pairR2 = $args{paired_r2};
	for my $index (0 .. $#{$pairR1}) {
		my $prefix = File::Spec->catfile(
			$args{work_dir}, "reads_$args{tag}.pair$index",
		);
		runSortmerna(
			%args,
			prefix => $prefix,
			read1 => $pairR1->[$index],
			read2 => $pairR2->[$index],
		);
		push @{$outputs{r1}}, findSortmernaOutput($prefix, 'fwd');
		push @{$outputs{r2}}, findSortmernaOutput($prefix, 'rev');
	}
	for my $index (0 .. $#{$args{single}}) {
		my $prefix = File::Spec->catfile(
			$args{work_dir}, "reads_$args{tag}.single$index",
		);
		runSortmerna(
			%args,
			prefix => $prefix,
			read1 => $args{single}->[$index],
			read2 => '',
		);
		push @{$outputs{single}}, findSortmernaOutput($prefix, 'single');
	}
	return %outputs;
}


sub runSortmerna {
	my (%args) = @_;
	my $referenceArgs = join ' ', map {
		'--ref '.shellQuote($_)
	} split(/:/, $args{reference});
	my $indexArgs = $args{index} ne ''
		? '--idx-dir '.shellQuote($args{index}).' --index 0'
		: '';
	my $workDir = "$args{prefix}.work";
	my $command = join ' ',
		$args{sortmerna},
		$referenceArgs,
		'--reads', shellQuote($args{read1}),
		($args{read2} ne ''
			? ('--reads', shellQuote($args{read2}))
			: ()),
		$indexArgs,
		'--workdir', shellQuote($workDir),
		'--aligned', shellQuote($args{prefix}),
		'--fastx',
		($args{read2} ne '' ? ('--paired_in', '--out2') : ()),
		'--no-best',
		'--zip-out', 1,
		'--threads', $args{threads},
		'-e', '1e-12',
		'--num_alignments', 1;
	systemW($command);
	remove_tree($workDir) if -d $workDir;
}


sub findSortmernaOutput {
	my ($prefix, $kind) = @_;
	my $stem = $kind eq 'fwd' ? "${prefix}_fwd"
		: $kind eq 'rev' ? "${prefix}_rev"
		: $prefix;
	my @candidates = map { "$stem.$_" } qw(fq.gz fa.gz);
	my @found = grep { -e $_ } @candidates;
	die "SortMeRNA produced no expected $kind output for $prefix\n"
		unless @found;
	die "SortMeRNA produced ambiguous $kind outputs for $prefix: @found\n"
		if @found > 1;
	assertGzipOutput($found[0]);
	return $found[0];
}


sub assertGzipOutput {
	my ($path) = @_;
	die "SortMeRNA output is empty or incomplete: $path\n" unless -s $path;
	open my $handle, '<', $path
		or die "Cannot inspect SortMeRNA output $path: $!\n";
	binmode $handle;
	my $magic = '';
	my $bytes = read($handle, $magic, 2);
	close $handle or die "Cannot close SortMeRNA output $path: $!\n";
	die "SortMeRNA output is not gzip-compressed despite --zip-out: $path\n"
		unless defined($bytes) && $bytes == 2 && $magic eq "\x1f\x8b";
}


sub publishMarkerOutputs {
	my ($workDir, $destinationDir, $tag, $outputs) = @_;
	my %suffix = (
		r1 => 'r1.fq.gz',
		r2 => 'r2.fq.gz',
		single => 'fq.gz',
	);
	for my $kind (qw(r1 r2 single)) {
		my $staged = File::Spec->catfile(
			$workDir, "reads_$tag.$suffix{$kind}.publish",
		);
		if (@{$outputs->{$kind}}) {
			concatenateFiles($outputs->{$kind}, $staged);
		} else {
			my $empty = '';
			gzip(\$empty => $staged)
				or die "Cannot create empty gzip output $staged: $GzipError\n";
		}
		assertGzipOutput($staged);
		my $destination = File::Spec->catfile(
			$destinationDir, "reads_$tag.$suffix{$kind}",
		);
		installAtomically($staged, $destination);
	}
}


sub concatenateFiles {
	my ($sources, $destination) = @_;
	open my $output, '>', $destination
		or die "Cannot create merged SortMeRNA output $destination: $!\n";
	binmode $output;
	for my $source (@{$sources}) {
		open my $input, '<', $source
			or die "Cannot read SortMeRNA output $source: $!\n";
		binmode $input;
		my $buffer;
		while (read($input, $buffer, 1024 * 1024)) {
			print {$output} $buffer
				or die "Cannot write merged SortMeRNA output $destination: $!\n";
		}
		die "Cannot finish reading SortMeRNA output $source: $!\n" if $!;
		close $input or die "Cannot close SortMeRNA output $source: $!\n";
	}
	close $output
		or die "Cannot close merged SortMeRNA output $destination: $!\n";
}


sub installAtomically {
	my ($source, $destination) = @_;
	my $temporary = "$destination.tmp.$$";
	unlinkChecked($temporary);
	copy($source, $temporary)
		or die "Cannot stage $source as $temporary: $!\n";
	rename($temporary, $destination)
		or die "Cannot install $destination atomically: $!\n";
}


sub touchFile {
	my ($path) = @_;
	open my $handle, '>', $path
		or die "Cannot create checkpoint $path: $!\n";
	close $handle or die "Cannot close checkpoint $path: $!\n";
}


sub unlinkChecked {
	my ($path) = @_;
	return unless -e $path || -l $path;
	unlink $path or die "Cannot remove stale RiboFind result $path: $!\n";
}


sub validateSortmernaReference {
	my ($configKey, $configuredReferences) = @_;
	die "SortMeRNA configuration '$configKey' has no reference path\n"
		if !defined($configuredReferences) || $configuredReferences eq '';
	for my $reference (split(/:/, $configuredReferences, -1)) {
		die "SortMeRNA configuration '$configKey' contains an empty reference path\n"
			if $reference eq '';
		die "Configured SortMeRNA reference '$configKey' does not exist: $reference\n"
			unless -e $reference;
	}
}


sub validateSortmernaIndex {
	my ($configKey, $indexDirectory) = @_;
	return if !defined($indexDirectory) || $indexDirectory eq '';
	die "Configured SortMeRNA index '$configKey' does not exist: $indexDirectory\n"
		unless -e $indexDirectory;
}


sub announceSortmerna {
	my ($sortmerna) = @_;
	print "catchLSUSSU v$VERSION\n";
	my $status = systemW("$sortmerna --version", 0);
	warn "Could not query the configured SortMeRNA version\n" if $status;
}


sub shellQuote {
	my ($value) = @_;
	$value = '' unless defined $value;
	$value =~ s/'/'"'"'/g;
	return "'$value'";
}
