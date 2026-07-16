#!/usr/bin/env perl
use strict;
use warnings;
use threads;
use File::Basename qw(dirname);
use File::Path qw(make_path);
use Mods::IO_Tamoc_progs qw(getProgPaths buildMapperIdx);

die "Usage: $0 <references[,references...]> <sample-dir> <output-db> <cores> <base-names> <final-dirs>\n"
	unless @ARGV == 6;
my ($reference_list, $sample_dir, $output_db, $cores, $base_names, $final_dirs) = @ARGV;
die "Core count must be a positive integer\n" unless $cores =~ /^\d+$/ && $cores > 0;
my @references = split /,/, $reference_list;
my @bases = split /,/, $base_names;
my @destinations = split /,/, $final_dirs;
die "Reference, base-name and destination counts differ\n"
	unless @references && @references == @bases && @references == @destinations;
for my $i (0 .. $#references) {
	die "Reference is missing or empty: $references[$i]\n" unless -s $references[$i];
	die "Destination directory does not exist: $destinations[$i]\n" unless -d $destinations[$i];
}
make_path(dirname($output_db)) unless -d dirname($output_db);

my $assembly_pointer = "$sample_dir/assemblies/metag/assembly.txt";
open my $pointer, '<', $assembly_pointer or die "Cannot open $assembly_pointer: $!\n";
my $assembly_dir = <$pointer>;
close $pointer;
die "Empty assembly pointer: $assembly_pointer\n" unless defined $assembly_dir;
$assembly_dir =~ s/[\r\n]+\z//;
my $assembly = "$assembly_dir/scaffolds.fasta.filt";
die "Cannot find sample assembly at $assembly\n" unless -s $assembly;

my (%contig_length, @contig_order);
open my $assembly_fh, '<', $assembly or die "Cannot open $assembly: $!\n";
my $current = '';
while (my $line = <$assembly_fh>) {
	$line =~ s/[\r\n]+\z//;
	if ($line =~ /^>(\S+)/) {
		$current = $1;
		die "Duplicate FASTA identifier '$current' in $assembly\n" if exists $contig_length{$current};
		$contig_length{$current} = 0;
		push @contig_order, $current;
		next;
	}
	die "Sequence before first FASTA header in $assembly\n" unless length $current;
	$line =~ s/\s+//g;
	$contig_length{$current} += length($line);
}
close $assembly_fh;
die "No FASTA records in $assembly\n" unless @contig_order;

my $blat = getProgPaths('blat');
my $pigz = getProgPaths('pigz');
my (@threads, @blast_files);
for my $i (0 .. $#references) {
	my $index = $i;
	my $blast_file = "$output_db.$i.b8";
	push @blast_files, $blast_file;
	$threads[$i] = threads->create(sub {
		my $status = system($blat, '-t=dna', '-q=dna', '-minIdentity=95', '-minScore=100',
			'-out=blast8', $references[$index], $assembly, $blast_file);
		return $status if $status != 0;
		my $compressed = "$destinations[$index]/$bases[$index].b8.gz";
		open my $gzip, '-|', $pigz, '-p', $cores, '-c', $blast_file or return 255;
		open my $out, '>', $compressed or return 255;
		while (my $chunk = <$gzip>) { print {$out} $chunk or return 255; }
		close $gzip or return 255;
		close $out or return 255;
		return 0;
	});
}
for my $i (0 .. $#threads) {
	my $status = $threads[$i]->join();
	die "Reference alignment failed for $references[$i] (status $status)\n" if $status != 0;
}

my $combined_blast = "$output_db.all.b8";
open my $combined, '>', $combined_blast or die "Cannot write $combined_blast: $!\n";
open my $sorted, '-|', 'sort', @blast_files or die "Cannot sort BLAT outputs: $!\n";
while (my $line = <$sorted>) { print {$combined} $line or die "Cannot write $combined_blast: $!\n"; }
close $sorted or die "Sorting BLAT outputs failed\n";
close $combined or die "Cannot close $combined_blast: $!\n";
my $removed = unlink @blast_files;
die "Cannot remove one or more temporary BLAT files: $!\n" unless $removed == @blast_files;

my %intervals;
open my $blast, '<', $combined_blast or die "Cannot open $combined_blast: $!\n";
my $line_number = 0;
while (my $line = <$blast>) {
	$line_number++;
	chomp $line;
	my @field = split /\t/, $line;
	die "Malformed BLAT output at line $line_number\n" unless @field >= 12;
	my ($query, $identity, $query_start, $query_end) = @field[0, 2, 6, 7];
	next unless $identity > 95;
	die "Unknown BLAT query '$query' at line $line_number\n" unless exists $contig_length{$query};
	my ($start, $end) = $query_start <= $query_end
		? ($query_start - 1, $query_end) : ($query_end - 1, $query_start);
	$start = 0 if $start < 0;
	$end = $contig_length{$query} if $end > $contig_length{$query};
	push @{$intervals{$query}}, [$start, $end] if $end > $start;
}
close $blast;

my %exclude;
for my $query (keys %intervals) {
	my @sorted_intervals = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$intervals{$query}};
	my ($covered, $start, $end) = (0, @{$sorted_intervals[0]});
	for my $interval (@sorted_intervals[1 .. $#sorted_intervals]) {
		if ($interval->[0] <= $end) {
			$end = $interval->[1] if $interval->[1] > $end;
		} else {
			$covered += $end - $start;
			($start, $end) = @$interval;
		}
	}
	$covered += $end - $start;
	$exclude{$query} = 1 if $covered >= $contig_length{$query} * 0.8;
}

open $assembly_fh, '<', $assembly or die "Cannot reopen $assembly: $!\n";
open my $database, '>', $output_db or die "Cannot write $output_db: $!\n";
my ($skip, $skipped_records, $skipped_bases) = (0, 0, 0);
while (my $line = <$assembly_fh>) {
	if ($line =~ /^>(\S+)/) {
		$skip = exists $exclude{$1};
		$skipped_records++ if $skip;
	}
	if ($skip) {
		if ($line !~ /^>/) { my $sequence = $line; $sequence =~ s/\s+//g; $skipped_bases += length($sequence); }
		next;
	}
	print {$database} $line or die "Cannot write $output_db: $!\n";
}
close $assembly_fh;
for my $reference (@references) {
	open my $reference_fh, '<', $reference or die "Cannot open $reference: $!\n";
	while (my $line = <$reference_fh>) { print {$database} $line or die "Cannot append $reference: $!\n"; }
	close $reference_fh;
}
close $database or die "Cannot close $output_db: $!\n";
print "Skipped $skipped_records FASTA entries ($skipped_bases bp); added ".scalar(@references)." reference FASTA(s).\n";

my ($build_command) = buildMapperIdx($output_db, $cores, 0, 0);
system($build_command) == 0 or die "Mapper index construction failed\n";
exit 0;
