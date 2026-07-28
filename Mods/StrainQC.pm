package Mods::StrainQC;

use strict;
use warnings;

use Exporter qw(import);
use Mods::GenoMetaAss qw(gzipopen);
use Mods::MGSLocus qw(robust_depth_mask);

our @EXPORT_OK = qw(
	breakpoint_gene_mask
	abundance_pattern_mask
);

sub abundance_pattern_mask {
	my ($depths, $options) = @_;
	return robust_depth_mask($depths, $options);
}

sub breakpoint_gene_mask {
	my ($gff, $breakpoint_file, $wanted, $flank) = @_;
	$wanted ||= {};
	$flank //= 50;
	return {} unless defined($gff) && length($gff) && -e $gff;
	return {} unless defined($breakpoint_file) && length($breakpoint_file) && -e $breakpoint_file;
	die "Breakpoint flank must be non-negative\n" if $flank < 0;

	my %breakpoints;
	my ($breakpoint_fh) = gzipopen($breakpoint_file, "strain breakpoint report", 1);
	my $line_number = 0;
	while (my $line = <$breakpoint_fh>) {
		$line_number++;
		$line =~ s/[\r\n]+$//;
		next if $line eq '' || $line =~ /^#/;
		my @fields = split /\t/, $line;
		next if $line_number == 1 && $fields[0] eq 'contig';
		die "Malformed breakpoint row $line_number in $breakpoint_file\n"
			unless @fields >= 3 && $fields[1] =~ /^\d+$/ && $fields[2] =~ /^\d+$/
				&& $fields[2] > $fields[1];
		push @{$breakpoints{$fields[0]}}, [0 + $fields[1], 0 + $fields[2]];
	}
	close $breakpoint_fh or die "Cannot close breakpoint report $breakpoint_file: $!\n";
	return {} unless keys %breakpoints;

	my %wanted_alias;
	for my $wanted_id (keys %{$wanted}) {
		$wanted_alias{$wanted_id} = $wanted_id;
		my $assembly_gene = $wanted_id;
		$assembly_gene =~ s/^.*__//;
		$wanted_alias{$assembly_gene} = $wanted_id;
	}

	my %masked;
	my ($gff_fh) = gzipopen($gff, "strain gene coordinates", 1);
	$line_number = 0;
	while (my $line = <$gff_fh>) {
		$line_number++;
		next if $line =~ /^#/;
		$line =~ s/[\r\n]+$//;
		next if $line eq '';
		my @fields = split /\t/, $line, -1;
		die "Malformed GFF row $line_number in $gff\n" unless @fields >= 9;
		my ($contig, $start, $end, $attributes) = @fields[0, 3, 4, 8];
		next unless exists $breakpoints{$contig};
		next unless $start =~ /^\d+$/ && $end =~ /^\d+$/ && $end >= $start;
		my ($id) = $attributes =~ /(?:^|;)ID=([^;]+)/;
		($id) = $attributes =~ /(?:^|;)(?:gene_id|locus_tag)=?["']?([^;"']+)/ unless defined $id;
		next unless defined($id) && length($id);
		my @ids = ($id, "${contig}_${id}");
		push @ids, "${contig}_$1" if $id =~ /_(\d+)$/;
		my $wanted_id = '';
		for my $candidate (@ids) {
			if (!keys(%{$wanted})) {
				$wanted_id = $candidate;
				last;
			}
			if (exists($wanted_alias{$candidate})) {
				$wanted_id = $wanted_alias{$candidate};
				last;
			}
		}
		next unless length $wanted_id;

		# GFF is one-based inclusive; breakpoint reports are zero-based half-open.
		my $gene_start = $start - 1;
		my $gene_end = $end;
		for my $break (@{$breakpoints{$contig}}) {
			my $break_start = $break->[0] - $flank;
			$break_start = 0 if $break_start < 0;
			my $break_end = $break->[1] + $flank;
			if ($gene_start < $break_end && $gene_end > $break_start) {
				$masked{$wanted_id} = {
					contig => $contig,
					break_start => $break->[0],
					break_end => $break->[1],
				};
				last;
			}
		}
	}
	close $gff_fh or die "Cannot close GFF $gff: $!\n";
	return \%masked;
}

1;
