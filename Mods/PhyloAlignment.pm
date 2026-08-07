package Mods::PhyloAlignment;

use strict;
use warnings;
use Exporter qw(import);

our @EXPORT_OK = qw(filter_alignment_by_overlap);

sub filter_alignment_by_overlap {
	my ($sequences, $is_amino_acid, $minimum_overlap) = @_;
	die "alignment sequences must be a hash reference\n"
		unless ref($sequences) eq 'HASH';
	die "minimum alignment overlap must be a nonnegative integer\n"
		unless defined($minimum_overlap)
			&& $minimum_overlap =~ /^\d+$/;

	my @ids = sort keys %{$sequences};
	return ({}, 0, 0) unless @ids;
	my $original_length = length($sequences->{$ids[0]} // '');
	for my $id (@ids) {
		die "alignment sequence '$id' is undefined\n"
			unless defined $sequences->{$id};
		die "alignment sequences have unequal lengths\n"
			unless length($sequences->{$id}) == $original_length;
	}

	my %filtered = map { $_ => $sequences->{$_} } @ids;
	return (\%filtered, $original_length, 0)
		if $minimum_overlap == 0 || $original_length == 0;

	my @called = (0) x $original_length;
	my $missing = $is_amino_acid ? qr/^[-Xx?.]$/ : qr/^[-Nn?.]$/;
	for my $id (@ids) {
		my @characters = split //, $sequences->{$id};
		for my $position (0 .. $#characters) {
			$called[$position]++ unless $characters[$position] =~ $missing;
		}
	}
	my @keep = grep { $called[$_] >= $minimum_overlap }
		0 .. $#called;
	for my $id (@ids) {
		my @characters = split //, $sequences->{$id};
		$filtered{$id} = join('', @characters[@keep]);
	}
	my $retained_length = scalar @keep;
	return (\%filtered, $retained_length,
		$original_length - $retained_length);
}

1;
