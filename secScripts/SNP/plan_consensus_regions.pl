#!/usr/bin/env perl
use strict;
use warnings;

use File::Basename qw(dirname);
use File::Path qw(make_path);
use Getopt::Long qw(GetOptions);

my %opt = (jobs => 1);
GetOptions(
	'fai=s'           => \$opt{fai},
	'mapping=s'       => \$opt{mapping},
	'depth=s'         => \$opt{depth},
	'jobs=i'          => \$opt{jobs},
	'output-prefix=s' => \$opt{output_prefix},
	'pigz=s'          => \$opt{pigz},
	'samtools=s'      => \$opt{samtools},
) or die "invalid region-planner arguments\n";

for my $required (qw(fai mapping output_prefix samtools)) {
	(my $flag = $required) =~ s/_/-/g;
	die "--$flag is required\n" unless defined($opt{$required}) && length($opt{$required});
}
die "--jobs must be a positive integer\n" unless $opt{jobs} && $opt{jobs} > 0;
die "FAI is missing or empty: $opt{fai}\n" unless -s $opt{fai};
die "mapping is missing or empty: $opt{mapping}\n" unless -s $opt{mapping};

my (@contigs, %length_for);
open my $fai_fh, '<', $opt{fai} or die "can't open $opt{fai}: $!\n";
while (my $line = <$fai_fh>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	my ($contig, $length) = split /\t/, $line;
	die "invalid FAI record in $opt{fai}: $line\n"
		unless defined($contig) && length($contig) && defined($length)
			&& $length =~ /^\d+$/ && $length > 0;
	die "duplicate contig '$contig' in $opt{fai}\n" if exists $length_for{$contig};
	push @contigs, $contig;
	$length_for{$contig} = 0 + $length;
}
close $fai_fh or die "can't close $opt{fai}: $!\n";
die "FAI contains no nonempty contigs: $opt{fai}\n" unless @contigs;

my %depth_for;
if (defined($opt{depth}) && length($opt{depth}) && -s $opt{depth}) {
	my $depth_fh;
	if ($opt{depth} =~ /\.gz$/) {
		die "--pigz is required for compressed depth input\n"
			unless defined($opt{pigz}) && length($opt{pigz});
		open $depth_fh, '-|', $opt{pigz}, '-dc', '--', $opt{depth}
			or die "can't start $opt{pigz} for $opt{depth}: $!\n";
	} else {
		open $depth_fh, '<', $opt{depth} or die "can't open $opt{depth}: $!\n";
	}
	while (my $line = <$depth_fh>) {
		chomp $line;
		next if $line =~ /^\s*$/;
		my ($contig, $depth) = split /\t/, $line;
		next unless exists $length_for{$contig};
		die "invalid depth for '$contig' in $opt{depth}: $depth\n"
			unless defined($depth) && $depth =~ /^\d+(?:\.\d+)?(?:[eE][+-]?\d+)?$/;
		$depth_for{$contig} = 0 + $depth;
	}
	close $depth_fh or die "can't close $opt{depth}: $!\n";
}

# Prefer measured per-contig depth. Supplementary calls commonly lack that
# file, so obtain a runtime workload estimate from their finished alignment.
if (!%depth_for) {
	open my $idx_fh, '-|', $opt{samtools}, 'idxstats', $opt{mapping}
		or die "can't run $opt{samtools} idxstats: $!\n";
	while (my $line = <$idx_fh>) {
		chomp $line;
		my ($contig, $length, $mapped) = split /\t/, $line;
		next unless defined($contig) && exists $length_for{$contig};
		die "invalid idxstats record: $line\n"
			unless defined($length) && $length =~ /^\d+$/
				&& defined($mapped) && $mapped =~ /^\d+$/;
		$depth_for{$contig} = $length > 0 ? $mapped / $length : 0;
	}
	close $idx_fh or die "$opt{samtools} idxstats failed for $opt{mapping}\n";
}

my $total_bases = 0;
$total_bases += $length_for{$_} for @contigs;
die "requested $opt{jobs} regions for only $total_bases reference bases\n"
	if $opt{jobs} > $total_bases;

# A constant reference cost retains zero-coverage contigs; measured depth
# apportions smaller regions to highly covered contigs.
my %weight_for = map { $_ => 1 + ($depth_for{$_} // 0) } @contigs;
my $remaining_cost = 0;
$remaining_cost += $length_for{$_} * $weight_for{$_} for @contigs;
my $remaining_bases = $total_bases;
my ($contig_index, $contig_start) = (0, 0);
my @regions;

for my $region_index (0 .. $opt{jobs} - 1) {
	my $remaining_jobs = $opt{jobs} - $region_index;
	my $target_cost = $remaining_cost / $remaining_jobs;
	my $region_cost = 0;
	while ($contig_index < @contigs) {
		my $contig = $contigs[$contig_index];
		my $available = $length_for{$contig} - $contig_start;
		if ($available <= 0) {
			$contig_index++;
			$contig_start = 0;
			next;
		}
		my $max_take = $remaining_bases - ($remaining_jobs - 1);
		my $take;
		if ($remaining_jobs == 1) {
			$take = $available;
		} else {
			my $needed = $target_cost - $region_cost;
			$take = int($needed / $weight_for{$contig});
			$take++ if $take * $weight_for{$contig} < $needed;
			$take = 1 if $take < 1;
			$take = $available if $take > $available;
			$take = $max_take if $take > $max_take;
		}
		last if $take < 1;
		my $stop = $contig_start + $take;
		$regions[$region_index] .= join("\t", $contig, $contig_start, $stop)."\n";
		my $cost = $take * $weight_for{$contig};
		$region_cost += $cost;
		$remaining_cost -= $cost;
		$remaining_bases -= $take;
		$contig_start = $stop;
		if ($contig_start == $length_for{$contig}) {
			$contig_index++;
			$contig_start = 0;
		}
		last if $remaining_jobs > 1 && $region_cost >= $target_cost;
	}
	die "failed to create nonempty region $region_index\n"
		unless defined($regions[$region_index]) && length($regions[$region_index]);
}
die "region planning omitted $remaining_bases reference bases\n" if $remaining_bases != 0;

my $output_dir = dirname($opt{output_prefix});
make_path($output_dir) unless -d $output_dir;
for my $region_index (0 .. $#regions) {
	my $bed = "$opt{output_prefix}$region_index.bed";
	open my $bed_fh, '>', $bed or die "can't write $bed: $!\n";
	print {$bed_fh} $regions[$region_index];
	close $bed_fh or die "can't close $bed: $!\n";
}

print "Planned ".scalar(@regions)." coverage-balanced SNP regions across $total_bases bases\n";
