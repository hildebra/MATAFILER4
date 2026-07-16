#!/usr/bin/perl
use warnings;
use strict;

sub usage {
	return "Usage: bamFilter.pl [max_edit_rate [min_query_coverage [min_mapq [min_end_clip]]]]\n";
}

sub fraction_arg {
	my ($value, $name) = @_;
	die usage() . "$name must be a number between 0 and 1 (received '$value')\n"
		unless defined($value)
		&& $value =~ /\A(?:\d+(?:\.\d*)?|\.\d+)\z/
		&& $value >= 0 && $value <= 1;
	return 0 + $value;
}

sub integer_arg {
	my ($value, $name, $maximum) = @_;
	die usage() . "$name must be a non-negative integer (received '$value')\n"
		unless defined($value) && $value =~ /\A\d+\z/
		&& (!defined($maximum) || $value <= $maximum);
	return 0 + $value;
}

sub parse_cigar {
	my ($cigar, $line_number) = @_;
	die "bamFilter: SAM line $line_number has an invalid mapped CIGAR '$cigar'\n"
		unless defined($cigar) && $cigar =~ /\A(?:[1-9]\d*[MIDNSHP=X])+\z/;

	my @operations;
	while ($cigar =~ /(\d+)([MIDNSHP=X])/g) {
		push @operations, [0 + $1, $2];
	}
	return @operations;
}

die usage() if @ARGV > 4;

# Positional arguments are retained for compatibility with existing MATAFILER
# command construction.  "max_edit_rate" is the maximum NM edit fraction,
# not percent identity (0.05 means at most 5% edits).
my $max_edit_rate = @ARGV >= 1 ? fraction_arg($ARGV[0], 'max_edit_rate') : 0.05;
my $min_query_coverage = @ARGV >= 2 ? fraction_arg($ARGV[1], 'min_query_coverage') : 0.8;
my $min_mapq = @ARGV >= 3 ? integer_arg($ARGV[2], 'min_mapq', 255) : 20;
my $min_end_clip = @ARGV >= 4 ? integer_arg($ARGV[3], 'min_end_clip') : 0;

my %counts = (
	records          => 0,
	retained         => 0,
	filtered         => 0,
	already_unmapped => 0,
	mapq             => 0,
	coverage         => 0,
	edit_rate        => 0,
	end_clipped      => 0,
);

my $input_line_number = 0;
while (my $line = <STDIN>) {
	$input_line_number++;
	if ($line =~ /^\@/) {
		print $line or die "bamFilter: cannot write SAM header: $!\n";
		next;
	}

	chomp $line;
	$line =~ s/\r\z//;
	die "bamFilter: upstream process reported an error at input line $input_line_number: $line\n"
		if $line =~ /^ERR:/;

	my @sam = split /\t/, $line, -1;
	die "bamFilter: malformed SAM at input line $input_line_number: expected at least 11 fields, found "
		. scalar(@sam) . "\n"
		if @sam < 11;

	$counts{records}++;
	die "bamFilter: SAM line $input_line_number has invalid FLAG '$sam[1]'\n"
		unless $sam[1] =~ /\A\d+\z/ && $sam[1] <= 65535;
	my $flag = 0 + $sam[1];

	# Records that were already unmapped do not contain enough alignment data to
	# evaluate and must pass through unchanged.
	if ($flag & 0x4) {
		print "$line\n" or die "bamFilter: cannot write SAM record: $!\n";
		$counts{already_unmapped}++;
		next;
	}

	die "bamFilter: mapped SAM line $input_line_number has no sequence\n"
		if $sam[9] eq '*';
	my $sequence_length = length($sam[9]);
	die "bamFilter: mapped SAM line $input_line_number has an empty sequence\n"
		if $sequence_length == 0;
	if ($sam[10] ne '*' && length($sam[10]) != $sequence_length) {
		die "bamFilter: SAM line $input_line_number has SEQ length $sequence_length but QUAL length "
			. length($sam[10]) . "\n";
	}
	die "bamFilter: SAM line $input_line_number has invalid MAPQ '$sam[4]'\n"
		unless $sam[4] =~ /\A\d+\z/ && $sam[4] <= 255;

	my @cigar_operations = parse_cigar($sam[5], $input_line_number);
	my ($cigar_query_length, $aligned_query_bases, $hard_clipped_bases, $deleted_bases) = (0, 0, 0, 0);
	for my $operation (@cigar_operations) {
		my ($length, $code) = @{$operation};
		$cigar_query_length += $length if $code =~ /[MIS=X]/;
		$aligned_query_bases += $length if $code =~ /[MI=X]/;
		$hard_clipped_bases += $length if $code eq 'H';
		$deleted_bases += $length if $code eq 'D';
	}
	die "bamFilter: SAM line $input_line_number has SEQ length $sequence_length but CIGAR consumes "
		. "$cigar_query_length query bases\n"
		if $cigar_query_length != $sequence_length;

	my $nm;
	for my $field (@sam[11 .. $#sam]) {
		if ($field =~ /\ANM:i:(\d+)\z/) {
			$nm = 0 + $1;
			last;
		}
	}
	die "bamFilter: mapped SAM line $input_line_number is missing a valid NM:i tag\n"
		unless defined($nm);

	# Hard-clipped bases are absent from SEQ but belong to the original query.
	# Insertions consume query sequence and therefore count as covered.  Deletions
	# add alignment columns for the NM edit-rate denominator; reference skips do not.
	my $original_query_length = $sequence_length + $hard_clipped_bases;
	my $query_coverage = $aligned_query_bases / $original_query_length;
	my $alignment_columns = $aligned_query_bases + $deleted_bases;
	die "bamFilter: mapped SAM line $input_line_number has no aligned query bases\n"
		if $alignment_columns == 0;
	my $edit_rate = $nm / $alignment_columns;
	my ($left_clip, $right_clip) = (0, 0);
	for my $operation (@cigar_operations) {
		last unless $operation->[1] =~ /[HS]/;
		$left_clip += $operation->[0];
	}
	for my $operation (reverse @cigar_operations) {
		last unless $operation->[1] =~ /[HS]/;
		$right_clip += $operation->[0];
	}

	my $failed = 0;
	if ($min_end_clip > 0
		&& $left_clip >= $min_end_clip
		&& $right_clip >= $min_end_clip) {
		$failed = 1;
		$counts{end_clipped}++;
	}
	if ($query_coverage < $min_query_coverage) {
		$failed = 1;
		$counts{coverage}++;
	}
	if ($edit_rate > $max_edit_rate) {
		$failed = 1;
		$counts{edit_rate}++;
	}
	if ($sam[4] < $min_mapq) {
		$failed = 1;
		$counts{mapq}++;
	}

	if ($failed) {
		$sam[1] = $flag | 0x4;
		$counts{filtered}++;
	} else {
		$counts{retained}++;
	}
	print join("\t", @sam), "\n" or die "bamFilter: cannot write SAM record: $!\n";
}

print STDERR "BamFilter\n"
	. "Input records: $counts{records}\n"
	. "Retained mapped records: $counts{retained}\n"
	. "Newly filtered records: $counts{filtered}\n"
	. "Already unmapped records: $counts{already_unmapped}\n"
	. "Failure reasons (non-exclusive):\n"
	. "  Mapping quality (<$min_mapq): $counts{mapq}\n"
	. "  Query coverage (<$min_query_coverage): $counts{coverage}\n"
	. "  Edit rate (>$max_edit_rate): $counts{edit_rate}\n"
	. "  Both ends clipped (>=$min_end_clip bases; 0 disables): $counts{end_clipped}\n"
	. "Done\n";

exit 0;
