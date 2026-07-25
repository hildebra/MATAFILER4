#!/usr/bin/env perl

use strict;
use warnings;

use FindBin qw($Bin);
use Getopt::Long qw(GetOptions);
use File::Basename qw(dirname);
use File::Path qw(make_path);

use lib "$Bin/../..";
use Mods::GenoMetaAss qw(gzipopen);

my ($input, $output) = ('', '');
GetOptions(
	'input=s'  => \$input,
	'output=s' => \$output,
) or die "Usage: $0 -input input.vcf[.gz] -output output.vcf\n";
die "-input is required\n" unless length($input) && -e $input;
die "-output is required\n" unless length($output);

my $output_dir = dirname($output);
make_path($output_dir) unless -d $output_dir;
my $temporary = "$output.tmp.$$";
END {
	unlink $temporary if defined($temporary) && -e $temporary;
}

my ($in, $ok) = gzipopen($input, "VCF header normalization", 0);
die "Cannot read VCF $input\n" unless $ok && defined($in);
open my $out, '>', $temporary or die "Cannot create $temporary: $!\n";

my $chrom_header = '';
my $duplicate_headers = 0;
my $late_metadata = 0;
while (my $line = <$in>) {
	if ($line =~ /^#CHROM(?:\t|$)/) {
		my $comparison = $line;
		$comparison =~ s/\r?\n\z//;
		if (!length($chrom_header)) {
			$chrom_header = $comparison;
			print {$out} $line or die "Cannot write $temporary: $!\n";
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
