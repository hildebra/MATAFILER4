#!/usr/bin/env perl

use strict;
use warnings;

use FindBin qw($Bin);
use Getopt::Long qw(GetOptions);
use File::Basename qw(dirname);
use File::Path qw(make_path);

use lib "$Bin/../..";
use Mods::GenoMetaAss qw(gzipopen gzipwrite);

my ($input, $output) = ('', '');
GetOptions(
	'input=s'  => \$input,
	'output=s' => \$output,
) or die "Usage: $0 -input input.vcf[.gz] -output output.vcf\n";
die "-input is required\n" unless length($input) && -e $input;
die "-output is required\n" unless length($output);
$output .= '.gz' unless $output =~ /\.gz\z/;

my $output_dir = dirname($output);
make_path($output_dir) unless -d $output_dir;
my $temporary = "$output.tmp.$$.gz";
END {
	unlink $temporary if defined($temporary) && -e $temporary;
}

my ($in, $ok) = gzipopen($input, "VCF header normalization", 0);
die "Cannot read VCF $input\n" unless $ok && defined($in);
my $out = gzipwrite($temporary, "normalized VCF");

my $chrom_header = '';
my $duplicate_headers = 0;
my $late_metadata = 0;
while (my $line = <$in>) {
	print "$line XX\n";
	if ($line =~ /^#CHROM(?:[ \t]+|$)/) {
		my $comparison = $line;
		$comparison =~ s/\r?\n\z//;
		die "$line\n";
		# VCF requires tab-delimited columns. Accept whitespace-delimited input
		# from older producers, but publish one canonical header representation.
		$comparison = join("\t", split /[ \t]+/, $comparison);
		if (!length($chrom_header)) {
			$chrom_header = $comparison;
			print {$out} "$chrom_header\n" or die "Cannot write $temporary: $!\n";
		} else {
			die "VCF $input contains incompatible #CHROM headers; refusing to merge sample columns\n"
				if $comparison ne $chrom_header;
			$duplicate_headers++;
		}
		next;
	}
	if (length($chrom_header) && $line =~ /^#/) {
		# Concatenated VCF chunks commonly repeat their metadata block before
		# the repeated #CHROM line. Metadata is valid only before #CHROM.
		$late_metadata++;
		next;
	}
	print {$out} $line or die "Cannot write $temporary: $!\n";
}
close $in or die "Cannot close $input: $!\n";
close $out or die "Cannot close $temporary: $!\n";
die "VCF $input has no #CHROM header\n" unless length($chrom_header);
rename $temporary, $output or die "Cannot replace $output: $!\n";



warn "Normalized $input: removed $duplicate_headers repeated #CHROM header(s)"
	.($late_metadata ? " and $late_metadata late metadata line(s)" : "")."\n"
	if $duplicate_headers || $late_metadata;
