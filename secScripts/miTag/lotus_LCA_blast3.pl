#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename qw(basename);
use File::Copy qw(move);
use File::Path qw(make_path remove_tree);
use File::Spec;
use Getopt::Long qw(GetOptions);
use Mods::GenoMetaAss qw(
	readFasta reverse_complement_IUPAC systemW
);
use Mods::IO_Tamoc_progs qw(getProgPaths setConfigFile);

my $lengthTolerance = 0.85;
my $maxReads = 50_000;
my $searchMode = 4;
my ($inputDir, $sample, $tmpRoot, $configFile, $databaseDir)
	= ('', '', '', '', '');
my $pairedMode = 0;
my $cores = 1;
my $keepReads = 0;

GetOptions(
	'dir=s'             => \$inputDir,
	'smplID=s'          => \$sample,
	'tmpD=s'            => \$tmpRoot,
	'pairedRds=i'       => \$pairedMode,
	'cores=i'           => \$cores,
	'simMode=i'         => \$searchMode,
	'config=s'          => \$configFile,
	'keepReads=i'       => \$keepReads,
	'maxReadNum=i'      => \$maxReads,
	'lengthTolerance=f' => \$lengthTolerance,
	'DBdir=s'           => \$databaseDir,
) or die "Error in command line arguments\n";
die "Unexpected positional arguments: @ARGV\n" if @ARGV;
die "Needs -dir, -smplID and -DBdir arguments\n"
	if $databaseDir eq '' || $inputDir eq '' || $sample eq '';
die "Input directory does not exist: $inputDir\n" unless -d $inputDir;
die "Database directory does not exist: $databaseDir\n" unless -d $databaseDir;
die "-cores must be a positive integer\n" if $cores < 1;
die "-pairedRds must be 0, 1, or 2\n"
	unless $pairedMode == 0 || $pairedMode == 1 || $pairedMode == 2;
die "-simMode must be 2 (LAMBDA), 3 (SortMeRNA), or 4 (VSEARCH)\n"
	unless $searchMode == 2 || $searchMode == 3 || $searchMode == 4;
die "-maxReadNum must be zero or greater\n" if $maxReads < 0;
die "-lengthTolerance must be between 0 and 1\n"
	if $lengthTolerance < 0 || $lengthTolerance > 1;

$inputDir = File::Spec->canonpath(File::Spec->rel2abs($inputDir));
$databaseDir = File::Spec->canonpath(File::Spec->rel2abs($databaseDir));
my $outputDir = File::Spec->catdir($inputDir, 'ltsLCA');
$tmpRoot = $tmpRoot eq ''
	? File::Spec->catdir($outputDir, 'tmp')
	: File::Spec->canonpath(File::Spec->rel2abs($tmpRoot));
make_path($outputDir) unless -d $outputDir;
make_path($tmpRoot) unless -d $tmpRoot;

my $safeSample = $sample;
$safeSample =~ s/[^A-Za-z0-9_.-]+/_/g;
my $workDir = File::Spec->catdir(
	$tmpRoot, 'lotusLCA_'.$safeSample.'_'.$$,
);
my $mergeDir = File::Spec->catdir($workDir, 'flashMerge');
make_path($mergeDir);

setConfigFile($configFile);
my $flash = $pairedMode ? getProgPaths('flash') : '';
my $lambda = $searchMode == 2 ? getProgPaths('lambda') : '';
my $sortmerna = $searchMode == 3 ? getProgPaths('sortmerna') : '';
my $vsearch = $searchMode == 4 ? getProgPaths('vsearch') : '';
my $lca = getProgPaths('LCA');

my @databaseKeys = qw(
	LSUdbFA LSUtax SSUdbFA SSUtax PR2dbFA PR2tax
);
my @configured = @{getProgPaths(\@databaseKeys, 0)};
for my $index (0 .. $#configured) {
	next if !defined($configured[$index]) || $configured[$index] eq '';
	$configured[$index] = File::Spec->catfile(
		$databaseDir, basename($configured[$index]),
	);
}
my ($lsuDb, $lsuTax, $ssuDb, $ssuTax, $pr2Db, $pr2Tax)
	= map { defined($_) ? $_ : '' } @configured;

my @tags = qw(SSU LSU);
my $allComplete = -e File::Spec->catfile($outputDir, 'Assigned.sto');
for my $tag (@tags) {
	$allComplete &&= markerAssignmentComplete($outputDir, $tag);
}
if ($allComplete) {
	print "All ribosomal assignments are already complete\n";
	remove_tree($workDir);
	exit 0;
}
unlinkChecked(File::Spec->catfile($outputDir, 'Assigned.sto'));

my $assignmentOK = 1;
for my $tag (@tags) {
	next if markerAssignmentComplete($outputDir, $tag);
	invalidateAssignmentOutput($outputDir, $tag);

	my (@databases, @taxonomies);
	if ($tag eq 'LSU') {
		@databases = ($lsuDb);
		@taxonomies = ($lsuTax);
	} else {
		@databases = ($pr2Db, $ssuDb);
		@taxonomies = ($pr2Tax, $ssuTax);
	}
	my (@usableDatabases, @usableTaxonomies);
	for my $index (0 .. $#databases) {
		next if !defined($databases[$index]) || $databases[$index] eq '';
		die "Missing $tag reference database: $databases[$index]\n"
			unless -e $databases[$index];
		die "Missing $tag taxonomy database: $taxonomies[$index]\n"
			unless defined($taxonomies[$index])
				&& $taxonomies[$index] ne ''
				&& -e $taxonomies[$index];
		push @usableDatabases, $databases[$index];
		push @usableTaxonomies, $taxonomies[$index];
	}
	die "No usable databases configured for $tag\n" unless @usableDatabases;

	my $singleInput = materializeReadInput(
		File::Spec->catfile($inputDir, 'reads_'.$tag.'.fq'),
		$workDir,
	);
	my $r1Input = materializeReadInput(
		File::Spec->catfile($inputDir, 'reads_'.$tag.'.r1.fq'),
		$workDir,
	);
	my $r2Input = materializeReadInput(
		File::Spec->catfile($inputDir, 'reads_'.$tag.'.r2.fq'),
		$workDir,
	);
	my $tagMode = $pairedMode;
	my $hasInput = 1;
	if ($tagMode > 0) {
		my $flashStatus = runFlash(
			$r1Input, $r2Input, $mergeDir, $tag, $inputDir,
		);
		if ($flashStatus == 1) {
			runLca(
				File::Spec->catfile($mergeDir, $tag),
				$singleInput,
				\@usableDatabases,
				\@usableTaxonomies,
				$tag.'riboRun_bl',
				$tag,
				$tagMode,
			);
		} elsif ($flashStatus == 3) {
			$tagMode = 0;
		} else {
			warn "Problem with $tag paired input; the profile must be repeated\n";
			$assignmentOK = 0;
			next;
		}
	}
	if ($tagMode == 0) {
		$hasInput = $singleInput ne '' && -s $singleInput ? 1 : 0;
		if ($hasInput) {
			runLca(
				$singleInput,
				'',
				\@usableDatabases,
				\@usableTaxonomies,
				$tag.'riboRun_bl',
				$tag,
				0,
			);
		} else {
			touchFile(File::Spec->catfile(
				$outputDir, $tag.'riboRun_bl.hiera.txt',
			));
		}
	}

	my $hierarchy = File::Spec->catfile(
		$outputDir, $tag.'riboRun_bl.hiera.txt',
	);
	unless (resultExists($hierarchy, !$hasInput)) {
		warn "Expected hierarchy output was not created for $tag: $hierarchy\n";
		$assignmentOK = 0;
		next;
	}
	touchFile(File::Spec->catfile($outputDir, $tag.'_ass.sto'));
}

remove_tree($workDir) if -d $workDir;
if ($assignmentOK) {
	for my $tag (@tags) {
		die "Cannot mark assignment complete; $tag output is incomplete\n"
			unless markerAssignmentComplete($outputDir, $tag);
	}
	touchFile(File::Spec->catfile($outputDir, 'Assigned.sto'));
	print File::Spec->catfile($outputDir, 'Assigned.sto'), "\n";
	exit 0;
}
die "Ribosomal LCA assignment was incomplete\n";


sub markerAssignmentComplete {
	my ($directory, $tag) = @_;
	return -e File::Spec->catfile($directory, $tag.'_ass.sto')
		&& resultExists(
			File::Spec->catfile(
				$directory, $tag.'riboRun_bl.hiera.txt',
			),
			1,
		);
}


sub invalidateAssignmentOutput {
	my ($directory, $tag) = @_;
	for my $path (
		File::Spec->catfile($directory, $tag.'_ass.sto'),
		File::Spec->catfile($directory, 'Assigned.sto'),
		File::Spec->catfile($directory, $tag.'riboRun_bl.hiera.txt'),
		File::Spec->catfile($directory, $tag.'riboRun_bl.hiera.txt.gz'),
	) {
		unlinkChecked($path);
	}
}


sub runFlash {
	my ($r1, $r2, $destination, $tag, $profileDir) = @_;
	return 3 if ($r1 eq '' || !-s $r1) && ($r2 eq '' || !-s $r2);
	if ($r1 eq '' || $r2 eq '' || !-s $r1 || !-s $r2) {
		warn "Incomplete paired input for $tag: '$r1' / '$r2'\n";
		unlinkChecked(File::Spec->catfile(
			$profileDir, $tag.'_pull.sto',
		));
		return 2;
	}
	my $r1Lines = countLines($r1);
	my $r2Lines = countLines($r2);
	if ($r1Lines != $r2Lines || $r1Lines % 4 != 0) {
		warn "Paired FASTQ files have incompatible record counts: "
			."$r1Lines / $r2Lines lines\n";
		unlinkChecked(File::Spec->catfile(
			$profileDir, $tag.'_pull.sto',
		));
		return 2;
	}
	print "Running FLASH for $tag\n";
	my $command = join ' ',
		$flash,
		'-M', 200,
		'-o', shellQuote($tag),
		'-d', shellQuote($destination),
		'-t', $cores,
		shellQuote($r1),
		shellQuote($r2);
	if (systemW($command, 0)) {
		warn "FLASH failed for $tag\n";
		unlinkChecked(File::Spec->catfile(
			$profileDir, $tag.'_pull.sto',
		));
		return 2;
	}
	for my $suffix (
		'extendedFrags.fastq',
		'notCombined_1.fastq',
		'notCombined_2.fastq',
	) {
		my $path = File::Spec->catfile(
			$destination, $tag.'.'.$suffix,
		);
		unless (-e $path) {
			warn "FLASH did not create expected output $path\n";
			unlinkChecked(File::Spec->catfile(
				$profileDir, $tag.'_pull.sto',
			));
			return 2;
		}
	}
	return 1;
}


sub fastqToFasta {
	my ($input, $deleteInput) = @_;
	return $input if !-s $input;
	open my $inputHandle, '<', $input
		or die "Input sequence file is not available: $input: $!\n";
	my $firstLine = <$inputHandle>;
	die "Input sequence file is empty: $input\n" unless defined $firstLine;
	if ($firstLine =~ /^>/) {
		close $inputHandle or die "Cannot close $input: $!\n";
		return $input;
	}
	seek($inputHandle, 0, 0) or die "Cannot rewind $input: $!\n";
	my $output = $input;
	$output =~ s/\.f(?:ast)?q$/.fa/i;
	die "Cannot derive a FASTA path from $input\n" if $output eq $input;
	open my $outputHandle, '>', $output
		or die "Cannot write FASTA output $output: $!\n";
	my $record = 0;
	while (my $header = <$inputHandle>) {
		my $sequence = <$inputHandle>;
		my $plus = <$inputHandle>;
		my $quality = <$inputHandle>;
		die "Truncated FASTQ record ".($record + 1)." in $input\n"
			unless defined($sequence) && defined($plus) && defined($quality);
		die "Invalid FASTQ header in record ".($record + 1)." of $input\n"
			unless $header =~ /^@/;
		die "Invalid FASTQ separator in record ".($record + 1)." of $input\n"
			unless $plus =~ /^\+/;
		chomp($header, $sequence, $quality);
		die "Sequence/quality length mismatch in record ".($record + 1)
			." of $input\n" unless length($sequence) == length($quality);
		$header =~ s/^@/>/;
		print {$outputHandle} "$header\n$sequence\n"
			or die "Cannot write FASTA output $output: $!\n";
		$record++;
	}
	close $inputHandle or die "Cannot close $input: $!\n";
	close $outputHandle or die "Cannot close $output: $!\n";
	unlinkChecked($input) if $deleteInput;
	return $output;
}


sub runLca {
	my (
		$queryBase, $extraSingle, $databases, $taxonomies,
		$outputName, $marker, $mode,
	) = @_;
	my $hierarchy = File::Spec->catfile(
		$outputDir, $outputName.'.hiera.txt',
	);
	my $query;
	if (!$mode) {
		die "Cannot find single-read input $queryBase\n" unless -e $queryBase;
		$query = fastqToFasta($queryBase, 0);
	} else {
		my $mergedFastq = $queryBase.'.extendedFrags.fastq';
		die "Missing FLASH merged output $mergedFastq\n"
			unless -e $mergedFastq;
		$query = fastqToFasta($mergedFastq, 1);

		my $unmergedR1 = $queryBase.'.notCombined_1.fastq';
		my $unmergedR2 = $queryBase.'.notCombined_2.fastq';
		my $interleaved = File::Spec->catfile(
			$workDir, 'inter'.$outputName.'.fa',
		);
		mergePairs($unmergedR1, $unmergedR2, $interleaved);
		if ($mode == 2 && $extraSingle ne '' && -s $extraSingle) {
			my $singleFasta = fastqToFasta($extraSingle, 0);
			appendFile($singleFasta, $interleaved);
		}
		appendFile($interleaved, $query) if -s $interleaved;
		unlinkChecked($interleaved);
	}
	unless (-s $query) {
		print "$outputName has empty input files\n";
		touchFile($hierarchy);
		return;
	}

	if ($maxReads > 0) {
		my $total = countFastaRecords($query);
		print "Found $total candidates in $query\n";
		if ($total > $maxReads) {
			print "Keeping the first $maxReads candidates\n";
			limitFastaRecords($query, $maxReads);
		}
	}

	my (@similarityOutputs, @usedTaxonomies);
	for my $index (0 .. $#{$databases}) {
		my $database = $databases->[$index];
		my $taxonomy = $taxonomies->[$index];
		my $similarity = File::Spec->catfile(
			$workDir, $outputName.'.'.$index.'.m8',
		);
		runSimilaritySearch(
			$query, $database, $similarity, $marker,
			scalar(@{$databases}), $index,
		);
		die "Similarity search produced no output file: $similarity\n"
			unless -e $similarity;
		push @similarityOutputs, $similarity;
		push @usedTaxonomies, $taxonomy;
	}

	my $flags = '-LCAfrac 0.8 -showHitRead -cover '
		.$lengthTolerance.' -minAlignLen 50';
	$flags .= ' -reportID' if $keepReads;
	my $command = join ' ',
		$lca,
		'-i', shellQuote(join(',', @similarityOutputs)),
		'-r', shellQuote(join(',', @usedTaxonomies)),
		'-o', shellQuote($hierarchy),
		$flags;
	print "Running LCA for $marker\n";
	systemW($command);
	die "LCA completed without producing $hierarchy\n" unless -e $hierarchy;

	extractReads($query, $hierarchy, File::Spec->catfile(
		$outputDir, $outputName.'.extr.fa',
	));
	unlinkChecked($_) for @similarityOutputs;
}


sub runSimilaritySearch {
	my (
		$query, $database, $output, $marker,
		$databaseCount, $databaseIndex,
	) = @_;
	if ($searchMode == 2) {
		my $index = $database.'.lba.gz';
		die "Cannot find required LAMBDA index $index\n" unless -e $index;
		my $columns = 'qseqid sseqid pident length mismatch gapopen '
			.'qstart qend sstart send qlen';
		my $command = join ' ',
			$lambda,
			'searchn',
			'-t', $cores,
			'--percent-identity', 75,
			'--num-matches', 200,
			'--e-value', '1e-8',
			'--output-columns', shellQuote($columns),
			'-q', shellQuote($query),
			'-i', shellQuote($index),
			'-o', shellQuote($output);
		systemW($command);
	} elsif ($searchMode == 3) {
		my $configuredIndex = $databaseCount == 1
			? getProgPaths($marker.'idx', 0)
			: '';
		if ($configuredIndex ne '') {
			die "Configured SortMeRNA index does not exist: $configuredIndex\n"
				unless -e $configuredIndex;
		}
		my $prefix = File::Spec->catfile(
			$workDir, 'sortmerna_'.$marker.'_'.$databaseIndex,
		);
		my $sortWork = $prefix.'.work';
		my $indexArgs = $configuredIndex ne ''
			? '--idx-dir '.shellQuote($configuredIndex).' --index 0'
			: '';
		my $command = join ' ',
			$sortmerna,
			'--ref', shellQuote($database),
			'--reads', shellQuote($query),
			$indexArgs,
			'--workdir', shellQuote($sortWork),
			'--blast', 1,
			'--aligned', shellQuote($prefix),
			'--zip-out', 0,
			'--threads', $cores,
			'-e', 0.1,
			'--num_alignments', 50;
		systemW($command);
		my $blastOutput = $prefix.'.blast';
		die "SortMeRNA produced no BLAST-format output $blastOutput\n"
			unless -e $blastOutput;
		move($blastOutput, $output)
			or die "Cannot move $blastOutput to $output: $!\n";
		remove_tree($sortWork) if -d $sortWork;
	} else {
		my $sharedIndex = $database.'.vudb';
		my $vsearchIndex = $sharedIndex;
		unless (-e $vsearchIndex) {
			$vsearchIndex = File::Spec->catfile(
				$workDir, basename($database).'.vudb',
			);
			unless (-e $vsearchIndex) {
				systemW(join ' ',
					$vsearch,
					'--makeudb_usearch', shellQuote($database),
					'--output', shellQuote($vsearchIndex),
				);
				die "VSEARCH did not create index $vsearchIndex\n"
					unless -e $vsearchIndex;
			}
		}
		my $fields = 'query+target+id+alnlen+mism+opens+'
			.'qlo+qhi+tlo+thi+ql';
		systemW(join ' ',
			$vsearch,
			'--usearch_global', shellQuote($query),
			'--db', shellQuote($vsearchIndex),
			'--id', 0.75,
			'--query_cov', $lengthTolerance,
			'--userfields', $fields,
			'--userout', shellQuote($output),
			'--maxaccepts', 100,
			'--maxrejects', 100,
			'--strand', 'both',
			'--threads', $cores,
		);
	}
}


sub mergePairs {
	my ($r1, $r2, $output) = @_;
	if ((!-e $r1 || !-s $r1) && (!-e $r2 || !-s $r2)) {
		touchFile($output);
		return;
	}
	die "FLASH produced only one nonempty unmerged mate: $r1 / $r2\n"
		unless -s $r1 && -s $r2;
	my $r1Fasta = fastqToFasta($r1, 1);
	my $r2Fasta = fastqToFasta($r2, 1);
	my %forward = %{readFasta($r1Fasta, 1)};
	my %reverse = %{readFasta($r2Fasta, 1)};
	my %reverseByPair;
	for my $id (keys %reverse) {
		my $pairId = $id;
		$pairId =~ s/\/2$//;
		die "Duplicate reverse-read identifier after pair normalization: $pairId\n"
			if exists $reverseByPair{$pairId};
		$reverseByPair{$pairId} = $reverse{$id};
	}
	open my $outputHandle, '>', $output
		or die "Cannot write combined reads $output: $!\n";
	for my $id (sort keys %forward) {
		my $pairId = $id;
		$pairId =~ s/\/1$//;
		die "Cannot find reverse mate for $id\n"
			unless exists $reverseByPair{$pairId};
		print {$outputHandle} '>'.$pairId."\n"
			.$forward{$id}
			.reverse_complement_IUPAC($reverseByPair{$pairId})
			."\n"
			or die "Cannot write combined reads $output: $!\n";
		delete $reverseByPair{$pairId};
	}
	die "Reverse-read file contains unmatched identifiers: "
		.join(',', sort keys %reverseByPair)."\n"
		if keys %reverseByPair;
	close $outputHandle or die "Cannot close combined reads $output: $!\n";
	unlinkChecked($r1Fasta);
	unlinkChecked($r2Fasta);
}


sub appendFile {
	my ($source, $destination) = @_;
	return unless -s $source;
	open my $input, '<', $source
		or die "Cannot read $source: $!\n";
	open my $output, '>>', $destination
		or die "Cannot append to $destination: $!\n";
	my $buffer;
	while (read($input, $buffer, 1024 * 1024)) {
		print {$output} $buffer
			or die "Cannot append to $destination: $!\n";
	}
	die "Cannot finish reading $source: $!\n" if $!;
	close $input or die "Cannot close $source: $!\n";
	close $output or die "Cannot close $destination: $!\n";
}


sub countFastaRecords {
	my ($path) = @_;
	open my $handle, '<', $path or die "Cannot read FASTA $path: $!\n";
	my $count = 0;
	while (my $line = <$handle>) {
		$count++ if $line =~ /^>/;
	}
	close $handle or die "Cannot close FASTA $path: $!\n";
	return $count;
}


sub limitFastaRecords {
	my ($path, $limit) = @_;
	my $temporary = $path.'.limit.'.$$;
	open my $input, '<', $path or die "Cannot read FASTA $path: $!\n";
	open my $output, '>', $temporary
		or die "Cannot write limited FASTA $temporary: $!\n";
	my $records = 0;
	while (my $line = <$input>) {
		if ($line =~ /^>/) {
			last if $records >= $limit;
			$records++;
		}
		print {$output} $line
			or die "Cannot write limited FASTA $temporary: $!\n";
	}
	close $input or die "Cannot close FASTA $path: $!\n";
	close $output or die "Cannot close limited FASTA $temporary: $!\n";
	rename($temporary, $path)
		or die "Cannot install limited FASTA $path: $!\n";
}


sub extractReads {
	my ($query, $hierarchy, $output) = @_;
	return unless $keepReads;
	my %sequences = %{readFasta($query, 1)};
	open my $input, '<', $hierarchy
		or die "Cannot read hierarchy file $hierarchy: $!\n";
	open my $outputHandle, '>', $output
		or die "Cannot write extracted reads $output: $!\n";
	while (my $line = <$input>) {
		chomp $line;
		next if $line eq '';
		my @fields = split(/\t/, $line, -1);
		next unless @fields && exists $sequences{$fields[0]};
		print {$outputHandle} '>'.$fields[0].' '.$fields[-1]."\n"
			.$sequences{$fields[0]}."\n"
			or die "Cannot write extracted reads $output: $!\n";
	}
	close $input or die "Cannot close hierarchy file $hierarchy: $!\n";
	close $outputHandle or die "Cannot close extracted reads $output: $!\n";
	systemW('gzip -f '.shellQuote($output));
}


sub materializeReadInput {
	my ($plainPath, $destinationRoot) = @_;
	return $plainPath if -e $plainPath;
	my $gzipPath = $plainPath.'.gz';
	return '' unless -e $gzipPath;
	my $inputDir = File::Spec->catdir($destinationRoot, 'inputs');
	make_path($inputDir) unless -d $inputDir;
	my $temporary = File::Spec->catfile(
		$inputDir, basename($plainPath),
	);
	systemW('gzip -dc '.shellQuote($gzipPath)
		.' > '.shellQuote($temporary));
	die "Could not materialize compressed read input $gzipPath\n"
		unless -e $temporary;
	return $temporary;
}


sub resultExists {
	my ($path, $allowEmpty) = @_;
	return 1 if -e $path && ($allowEmpty || -s $path);
	return 1 if -e $path.'.gz' && ($allowEmpty || -s $path.'.gz');
	return 0;
}


sub touchFile {
	my ($path) = @_;
	open my $handle, '>', $path
		or die "Cannot create checkpoint/output $path: $!\n";
	close $handle or die "Cannot close checkpoint/output $path: $!\n";
}


sub unlinkChecked {
	my ($path) = @_;
	return unless -e $path || -l $path;
	unlink $path or die "Cannot remove stale file $path: $!\n";
}


sub countLines {
	my ($path) = @_;
	open my $handle, '<', $path
		or die "Cannot count records in $path: $!\n";
	my $count = 0;
	$count++ while <$handle>;
	close $handle or die "Cannot close $path: $!\n";
	return $count;
}


sub shellQuote {
	my ($value) = @_;
	$value = '' unless defined $value;
	$value =~ s/'/'"'"'/g;
	return "'$value'";
}
