package Mods::StrainPlacement;

use strict;
use warnings;
use Exporter qw(import);
use Mods::GenoMetaAss qw(fileGZe gzipopen);

our @EXPORT_OK = qw(
	read_sample_qc
	split_strict_backbone
	nearest_backbone_placements
	write_placed_tree
);

sub _read_fasta {
	my ($file) = @_;
	open my $fh, '<', $file or die "Cannot read FASTA $file: $!\n";
	my (%seq, $id);
	while (my $line = <$fh>) {
		$line =~ s/[\r\n]+$//;
		next if $line eq '';
		if ($line =~ /^>(\S+)/) {
			$id = $1;
			die "Duplicate FASTA identifier '$id' in $file\n" if exists $seq{$id};
			$seq{$id} = '';
			next;
		}
		die "Sequence before a FASTA header in $file\n" unless defined $id;
		$seq{$id} .= $line;
	}
	close $fh or die "Cannot close FASTA $file: $!\n";
	return \%seq;
}

sub _write_fasta {
	my ($file, $seq, $ids) = @_;
	open my $fh, '>', $file or die "Cannot write FASTA $file: $!\n";
	for my $id (@{$ids}) {
		print {$fh} ">$id\n$seq->{$id}\n";
	}
	close $fh or die "Cannot close FASTA $file: $!\n";
}

sub _informative_count {
	my ($sequence, $is_aa) = @_;
	my $copy = $sequence // '';
	return $is_aa
		? ($copy =~ tr/A-WYZa-wyz//)
		: ($copy =~ tr/ACGTUacgtu//);
}

sub _quantile {
	my ($fraction, @values) = @_;
	return 0 unless @values;
	@values = sort { $a <=> $b } @values;
	my $position = $fraction * $#values;
	my $lower = int($position);
	my $upper = $lower == $position ? $lower : $lower + 1;
	return $values[$lower] if $lower == $upper;
	my $weight = $position - $lower;
	return $values[$lower] * (1 - $weight) + $values[$upper] * $weight;
}

sub read_sample_qc {
	my ($file) = @_;
	return {} unless defined($file) && length($file) && fileGZe($file);
	my ($fh) = gzipopen($file, "sample QC file", 1);
	my %status;
	while (my $line = <$fh>) {
		$line =~ s/[\r\n]+$//;
		next if $line eq '' || $line =~ /^#/ || $line =~ /^MGS\t/;
		my @field = split /\t/, $line, -1;
		die "Malformed sample QC row in $file: $line\n" unless @field >= 6;
		my ($sample, $sample_status) = @field[1, 2];
		die "Invalid sample QC status '$sample_status' for $sample\n"
			unless $sample_status eq 'backbone' || $sample_status eq 'placement';
		$status{$sample} = $sample_status
			if !exists($status{$sample}) || $sample_status eq 'placement';
	}
	close $fh or die "Cannot close sample QC file $file: $!\n";
	return \%status;
}

sub split_strict_backbone {
	my ($full_fasta, $backbone_fasta, $placement_fasta, $status, $options) = @_;
	$options ||= {};
	my $is_aa = $options->{is_aa} // 0;
	my $minimum_backbone = $options->{minimum_backbone} // 3;
	my $coverage_fraction = $options->{coverage_fraction} // 0.35;
	my $outgroup = $options->{outgroup} // '';
	my $seq = _read_fasta($full_fasta);
	my @ids = sort keys %{$seq};
	die "Strict-backbone input $full_fasta contains no sequences\n" unless @ids;
	my %informative = map { $_ => _informative_count($seq->{$_}, $is_aa) } @ids;
	my $q90 = _quantile(0.90, values %informative);
	my (%classification_reason, %requested_reason, @backbone, @placement);
	for my $id (@ids) {
		my $sampleLocusQC = ($status->{$id} // '') eq 'placement';
		my $lowCoverage = $q90 > 0
			&& $informative{$id} < $coverage_fraction * $q90;
		$lowCoverage = 0 if length($outgroup) && $id eq $outgroup;
		if ($lowCoverage) {
			my @reason = ('low_validated_coverage');
			push @reason, 'sample_locus_qc' if $sampleLocusQC;
			push @placement, $id;
			$classification_reason{$id} = join(',', @reason);
			$requested_reason{$id} = $classification_reason{$id};
		} else {
			push @backbone, $id;
			$classification_reason{$id} = 'retained_after_locus_qc_masking'
				if $sampleLocusQC;
		}
	}
	my $fallback = 0;
	if (@backbone < $minimum_backbone) {
		$fallback = 1;
		@backbone = @ids;
		@placement = ();
	}
	_write_fasta($backbone_fasta, $seq, \@backbone);
	_write_fasta($placement_fasta, $seq, \@placement);
	return {
		backbone => \@backbone,
		placement => \@placement,
		reason => \%classification_reason,
		informative => \%informative,
		q90_informative => $q90,
		fallback => $fallback,
		requested_reason => \%requested_reason,
	};
}

sub _pair_distance {
	my ($left, $right, $is_aa) = @_;
	die "Cannot compare unequal alignment lengths\n" unless length($left) == length($right);
	my ($overlap, $different) = (0, 0);
	for (my $i = 0; $i < length($left); $i++) {
		my $a = uc substr($left, $i, 1);
		my $b = uc substr($right, $i, 1);
		my $valid = $is_aa
			? ($a =~ /^[ACDEFGHIKLMNPQRSTVWY]$/ && $b =~ /^[ACDEFGHIKLMNPQRSTVWY]$/)
			: ($a =~ /^[ACGTU]$/ && $b =~ /^[ACGTU]$/);
		next unless $valid;
		$overlap++;
		$different++ if $a ne $b;
	}
	return ($overlap ? $different / $overlap : undef, $overlap);
}

sub nearest_backbone_placements {
	my ($backbone_fasta, $placement_fasta, $minimum_overlap, $is_aa) = @_;
	$minimum_overlap //= 400;
	my $backbone = _read_fasta($backbone_fasta);
	my $queries = _read_fasta($placement_fasta);
	my %placement;
	for my $query (sort keys %{$queries}) {
		my ($best_anchor, $best_distance, $best_overlap);
		for my $anchor (sort keys %{$backbone}) {
			my ($distance, $overlap) = _pair_distance(
				$queries->{$query}, $backbone->{$anchor}, $is_aa,
			);
			next unless defined($distance) && $overlap >= $minimum_overlap;
			if (!defined($best_distance)
				|| $distance < $best_distance
				|| ($distance == $best_distance && $overlap > $best_overlap)
				|| ($distance == $best_distance && $overlap == $best_overlap
					&& $anchor lt $best_anchor)) {
				($best_anchor, $best_distance, $best_overlap) = ($anchor, $distance, $overlap);
			}
		}
		$placement{$query} = {
			anchor => $best_anchor // '',
			distance => $best_distance,
			overlap => $best_overlap // 0,
			status => defined($best_distance) ? 'placed' : 'insufficient_overlap',
		};
	}
	return \%placement;
}

sub _newick_label_pattern {
	my ($label) = @_;
	my $quoted = quotemeta($label);
	return qr/(?<![A-Za-z0-9_.-])($quoted)(?=[:),;])/;
}

sub write_placed_tree {
	my ($backbone_tree, $output_tree, $placements) = @_;
	open my $in, '<', $backbone_tree or die "Cannot read backbone tree $backbone_tree: $!\n";
	local $/;
	my $tree = <$in>;
	close $in or die "Cannot close backbone tree $backbone_tree: $!\n";
	for my $query (sort keys %{$placements}) {
		my $entry = $placements->{$query};
		next unless ($entry->{status} // '') eq 'placed';
		my $anchor = $entry->{anchor};
		my $pattern = _newick_label_pattern($anchor);
		my $distance = sprintf('%.8g', $entry->{distance});
		my $replacement = "($anchor:0,$query:$distance)";
		my $count = ($tree =~ s/$pattern/$replacement/);
		die "Could not find unique placement anchor '$anchor' in $backbone_tree\n"
			unless $count == 1;
	}
	my $temporary_output = "$output_tree.tmp.$$";
	unlink $temporary_output
		or die "Cannot remove stale placed-tree temporary $temporary_output: $!\n"
		if -e $temporary_output;
	open my $out, '>', $temporary_output
		or die "Cannot write placed tree $temporary_output: $!\n";
	print {$out} $tree;
	close $out or die "Cannot close placed tree $temporary_output: $!\n";
	rename $temporary_output, $output_tree
		or die "Cannot publish placed tree $temporary_output as $output_tree: $!\n";
}

1;
