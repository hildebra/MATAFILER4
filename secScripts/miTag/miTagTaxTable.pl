#!/usr/bin/env perl

use strict;
use warnings;

use Mods::GenoMetaAss qw(gzipopen systemW);

sub shellQuote;

die "Usage: $0 tax_level[,tax_level...] output_prefix input_directory\n" unless @ARGV == 3;
my $taxLevelArg = lc shift @ARGV;
my $outPrefix = shift @ARGV;
my $inputDir = shift @ARGV;
my @levels = grep { $_ ne "" } split(/,/, $taxLevelArg);
die "At least one taxonomic level is required\n" unless @levels;
die "Output prefix must not be empty\n" if $outPrefix eq "";
die "Input directory does not exist: $inputDir\n" unless -d $inputDir;

my %seenLevel;
die "Duplicate taxonomic levels are not supported\n" if grep { $seenLevel{$_}++ } @levels;

for my $level (@levels){
	unlink "$outPrefix.$level.txt" if -e "$outPrefix.$level.txt";
	unlink "$outPrefix.$level.txt.gz" if -e "$outPrefix.$level.txt.gz";
}

opendir(my $dirHandle, $inputDir) or die "Cannot open directory $inputDir: $!\n";
my @files = sort grep { /\.hiera\.txt(?:\.gz)?$/ && -f "$inputDir/$_" } readdir($dirHandle);
closedir($dirHandle) or die "Cannot close directory $inputDir: $!\n";

print "Detected ".scalar(@files)." input files in dir $inputDir\n";
exit(0) unless @files;

my %sites;
my %taxa;
my %seenTag;

for my $file (@files) {
	my ($inputHandle,$readOk) = gzipopen("$inputDir/$file", "tax infile");
	my $tag = $file;
	$tag =~ s/\.hiera\.txt(?:\.gz)?$//;
	die "Duplicate sample tag '$tag' derived from hierarchy inputs\n" if $seenTag{$tag}++;

	my %column;
	my $header = <$inputHandle>;
	die "Empty taxonomy hierarchy input: $inputDir/$file\n" unless defined $header;
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
	close $inputHandle or die "Cannot close taxonomy input $inputDir/$file: $!\n";
}

print "Read input files.\n";
for my $level (@levels){
	my @taxaKeys = sort {
		$taxa{$level}{$b} <=> $taxa{$level}{$a} || $a cmp $b
	} keys %{$taxa{$level}};
	my @siteKeys = sort keys %{$sites{$level}};
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
