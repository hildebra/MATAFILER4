#!/usr/bin/env perl
###############################################################################
#
#    calc.kmerfreq.pl
#
#	 Calculates kmer frequency in fasta files. 
#    Length restriction, variable kmer, and rc kmer to remove strand bias.
#    
#    Copyright (C) 2012 Mads Albertsen
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

#pragmas
use strict;
use warnings;

#core Perl modules
use Getopt::Long;

#locally-written modules
BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# get input params
my $global_options = checkParams();

my $fastafile;;
my $outfile;
my $minlength;
my $kmerlength;
my $split;

$fastafile = &overrideDefault("inputfile.fasta",'fastafile');
$outfile = &overrideDefault("outfile.tab",'outfile');
$minlength = &overrideDefault("500",'minlength');
$kmerlength = &overrideDefault("4",'kmer');
$split = &overrideDefault("0",'split');

my $kmerlengthprobe="";
my $header;
my $seq;
my $prevheader;
my $subseq;
my $sequence = "";
my $output = "Contig";
my @probes;
my @kmerheader;
my %kmer;

######################################################################
# CODE HERE
######################################################################

if ($split != 0){
	my $counttemp = 0;
	open(INtemp, $fastafile) or die("Cannot open $fastafile\n");
	open(OUTtemp, ">$fastafile.sub.$split.fa") or die("Cannot open $fastafile\n");
	while ( my $line = <INtemp>) {
		chomp $line;
		if ($line =~ m/>/) {		
			$prevheader = $header;
			$header = $line;	
			my $count2 = 0;
			if (($minlength <= length($seq)) and ($counttemp > 0)) {
				for ($count2 = 0; ($count2+$split*2) < length($seq); $count2=$count2+$split)  {
					my $enddist = $count2+$split;
					print OUTtemp "$prevheader.$count2.$enddist\n";
					$subseq = substr($seq,$count2,$split);
					print OUTtemp "$subseq\n";			
				}
			my $endlength = length($seq);
			$subseq = substr($seq,$count2,$endlength);			
			print OUTtemp "$prevheader.$count2.$endlength\n";
			print OUTtemp "$subseq\n";	
			}
			$counttemp++;
			$seq = "";
		}
		else{			
			$seq = $seq.$line;
		}
	}
	my $count2 = 0;
	if ($minlength <= length($seq)) {
		my $endlength2 = length($seq);
		print "$endlength2\n";
		print "here\n";
		for ($count2 = 0; ($count2+$split*2) < length($seq); $count2=$count2+$split)  {
			my $enddist = $count2+$split;
			print OUTtemp "$header.$count2.$enddist\n";
			$subseq = substr($seq,$count2,$split);
			print OUTtemp "$subseq\n";			
		}
		my $endlength = length($seq);
		$subseq = substr($seq,$count2,$endlength);			
		print OUTtemp "$header.$count2.$endlength\n";
		print OUTtemp "$subseq\n";	
	}
	close INtemp;
	close OUTtemp;
	$fastafile = "$fastafile.sub.$split.fa";	
	$minlength = $split-1;
}




open(IN, $fastafile) or die("Cannot open $fastafile\n");
open(OUT, ">$outfile") or die("Cannot create $outfile\n");

for (my $count2 = 0; $count2 < $kmerlength; $count2++)  {                                          #Create all possible kmers
	$kmerlengthprobe = $kmerlengthprobe."N";
}

push (@probes,$kmerlengthprobe);

foreach my $probe (@probes){
	if ($probe =~ m/N/) { 								
		my $temp1 = $probe;
		$temp1 =~ s/N/A/;											
		push (@probes, "$temp1");						
		$temp1 = $probe;
		$temp1 =~ s/N/T/;											
		push (@probes, "$temp1");						
		$temp1 = $probe;
		$temp1 =~ s/N/C/;											
		push (@probes, "$temp1");						
		$temp1 = $probe;
		$temp1 =~ s/N/G/;											
		push (@probes, "$temp1");								
	}				
}

foreach my $probe (@probes){
	if ($probe =~ m/N/) { 	
	}
	else{
		$output = $output."\t$probe";
		push (@kmerheader,$probe);
		$kmer{$probe} = 0;
	}
}
print OUT "$output\n";	

# Use exactly the same emission rule at a header boundary and at EOF.
$sequence = "";
$header = undef;
while (my $line = <IN>) {
	chomp $line;
	if ($line =~ /^>(.*)/) {
		emit_kmers($header, $sequence);
		$header = $1;
		$sequence = "";
	} else {
		$sequence .= uc($line);
	}
}
emit_kmers($header, $sequence);

sub emit_kmers {
	my ($name, $dna) = @_;
	return unless defined($name) && length($dna) >= $minlength;
	$kmer{$_} = 0 for @kmerheader;
	my $total = 0;
	for (my $i = 0; $i + $kmerlength <= length($dna); $i++) {
		my $probe = substr($dna, $i, $kmerlength);
		next unless exists $kmer{$probe};
		$kmer{$probe}++;
		my $rc = reverse($probe);
		$rc =~ tr/ACGT/TGCA/;
		$kmer{$rc}++;
		$total += 2;
	}
	# Preserve the existing non-final-record policy: omit uninformative records.
	return unless $total;
	print OUT join("\t", $name,
		map { sprintf("%.5f", $kmer{$_} / $total) * 100 } @kmerheader), "\n";
}


close IN;
close OUT;

exit;

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "fastafile|i:s", "outfile|o:s", "minlength|m:s", "kmer|k:s", "split|s:s");
    my %options;

    # Add any other command line options, and the code to handle them
    # 
    GetOptions( \%options, @standard_options );
    
	#if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # Compulsosy items
    #if(!exists $options{'infile'} ) { print "**ERROR: $0 : \n"; exec("pod2usage $0"); }

    return \%options;
}

sub overrideDefault
{
    #-----
    # Set and override default values for parameters
    #
    my ($default_value, $option_name) = @_;
    if(exists $global_options->{$option_name}) 
    {
        return $global_options->{$option_name};
    }
    return $default_value;
}

__DATA__

=head1 NAME

    calc.kmer.freq.pl

=head1 COPYRIGHT

   copyright (C) 2012 Mads Albertsen

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 DESCRIPTION

	Calculates kmer frequency in fasta files. 
	Length restriction, variable kmer, and rc kmer to remove strand bias.

=head1 SYNOPSIS

script.pl  -i [-h -o -m -k]

 [-help -h]           Displays this basic usage information
 [-fastafile -i]      Input fastafile file. 
 [-outputfile -o]     Outputfile. Tab separated tetranucleotide frequency.
 [-minlength -m]      Minimum contig length to use for calculations (default: 500).
 [-kmer -k]           kmer length (default: 4).
 [-split -s]          Split sequences in subsequences of length -s (default: no split).
 
=cut
