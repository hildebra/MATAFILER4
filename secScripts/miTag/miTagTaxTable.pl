#!/usr/bin/env perl

use strict;
use warnings;

use Mods::GenoMetaAss qw(gzipopen systemW);

sub shellQuote;

# The input is a directory of <sample>.hiera.txt[.gz] files, or a sample list
# with one "<column>\t<hierarchy file>" line per sample (MATAF4 writes one per
# marker, for exactly the samples of its map).
die "Usage: $0 tax_level[,tax_level...] output_prefix input_directory|sample_list\n" unless @ARGV == 3;
my $taxLevelArg = lc shift @ARGV;
my $outPrefix = shift @ARGV;
my $input = shift @ARGV;
my @levels = grep { $_ ne "" } split(/,/, $taxLevelArg);
die "At least one taxonomic level is required\n" unless @levels;
die "Output prefix must not be empty\n" if $outPrefix eq "";
die "Input directory or sample list does not exist: $input\n" unless -e $input;

my %seenLevel;
die "Duplicate taxonomic levels are not supported\n" if grep { $seenLevel{$_}++ } @levels;

for my $level (@levels){
	unlink "$outPrefix.$level.txt" if -e "$outPrefix.$level.txt";
	unlink "$outPrefix.$level.txt.gz" if -e "$outPrefix.$level.txt.gz";
}

my @inputs; # [sample column, hierarchy file]
if (-d $input) {
	opendir(my $dirHandle, $input) or die "Cannot open directory $input: $!\n";
	my @entries = grep { /\.hiera\.txt(?:\.gz)?$/ } readdir($dirHandle);
	closedir($dirHandle) or die "Cannot close directory $input: $!\n";
	# A link whose target is gone drops that sample from every table, so name it.
	my @dangling = sort grep { -l "$input/$_" && !-e "$input/$_" } @entries;
	warn "Skipping ".scalar(@dangling)." hierarchy link(s) whose sample output no longer exists: "
		.join(', ', @dangling)."\n" if @dangling;
	for my $file (sort grep { -f "$input/$_" } @entries) {
		(my $tag = $file) =~ s/\.hiera\.txt(?:\.gz)?$//;
		push @inputs, [$tag, "$input/$file"];
	}
} else {
	open my $listHandle, '<', $input or die "Cannot read sample list $input: $!\n";
	while (my $line = <$listHandle>) {
		$line =~ s/\r?\n$//;
		next if $line eq "";
		my ($tag, $path) = split /\t/, $line, 2;
		die "Sample list line $. needs a column name and a hierarchy file: $line\n"
			unless defined($path) && $tag ne "" && $path ne "";
		# every listed sample belongs in the tables; never drop one silently
		die "Hierarchy of $tag does not exist: $path\n" unless -f $path;
		push @inputs, [$tag, $path];
	}
	close $listHandle or die "Cannot close sample list $input: $!\n";
}

print "Detected ".scalar(@inputs)." input files in $input\n";
exit(0) unless @inputs;

my %sites;
my %taxa;
my %seenTag;

for my $entry (@inputs) {
	my ($tag, $file) = @{$entry};
	my ($inputHandle,$readOk) = gzipopen($file, "tax infile");
	die "Duplicate sample tag '$tag' derived from hierarchy inputs\n" if $seenTag{$tag}++;

	my %column;
	my $header = <$inputHandle>;
	# A zero-byte hierarchy is the producer's completed-empty representation.
	# Keep its sample tag even though there are no taxa or header to parse.
	unless (defined $header) {
		close $inputHandle or die "Cannot close taxonomy input $file: $!\n";
		next;
	}
	chomp $header;
	my @headerFields = split /\t/, $header, -1;
	for my $level (@levels){
		for (my $i=0; $i<@headerFields; $i++) {
			if (lc($headerFields[$i]) eq $level) {
				$column{$level} = $i - 1; # data rows discard the leading read identifier
				last;
			}
		}
		die "Could not find taxonomic level '$level' in $file\n"
			unless exists($column{$level}) && $column{$level} >= 0;
	}

	while (my $row=<$inputHandle>) {
		chomp $row;
		next if $row eq "";
		my @fields = split /\t/, $row, -1;
		shift @fields; # discard read identifier, matching the header offset above

		# PR2 contains an extra Opisthokonta supergroup. Preserve the historical
		# normalization while guarding short/malformed records.
		if (@fields > 1 && $fields[1] eq "Opisthokonta"){
			splice @fields, 1, 1;
			splice @fields, 4, 0, "?";
		}

		for my $level (@levels){
			my $lastColumn = $column{$level};
			my @lineage;
			for my $index (0..$lastColumn){
				my $value = $index < @fields && defined($fields[$index]) && $fields[$index] ne ""
					? $fields[$index] : "?";
				push @lineage, $value;
			}
			my $lineage = join(';', @lineage);
			$sites{$level}{$tag}{$lineage}++;
			$taxa{$level}{$lineage}++;
		}
	}
	close $inputHandle or die "Cannot close taxonomy input $file: $!\n";
}

print "Read input files.\n";
for my $level (@levels){
	my @taxaKeys = sort {
		$taxa{$level}{$b} <=> $taxa{$level}{$a} || $a cmp $b
	} keys %{$taxa{$level} // {}};
	my @siteKeys = sort keys %seenTag;
	my $output = "$outPrefix.$level.txt";
	open my $outputHandle, '>', $output or die "Cannot write $output: $!\n";
	print {$outputHandle} $level, map { "\t$_" } @siteKeys;
	print {$outputHandle} "\n";
	for my $lineage (@taxaKeys) {
		print {$outputHandle} $lineage;
		for my $site (@siteKeys) {
			print {$outputHandle} "\t", ($sites{$level}{$site}{$lineage} // 0);
		}
		print {$outputHandle} "\n";
	}
	close $outputHandle or die "Cannot close $output: $!\n";
	systemW("gzip -f ".shellQuote($output));
}


sub shellQuote{
	my ($value) = @_;
	$value = "" unless defined $value;
	$value =~ s/'/'"'"'/g;
	return "'$value'";
}
