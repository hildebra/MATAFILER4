package Mods::MosaicLoci;

use strict;
use warnings;

use Exporter qw(import);
use Mods::GenoMetaAss qw(gzipopen);

our @EXPORT_OK = qw(
	pair_key
	discover_mosaic_candidates
	read_paf_hits
	confirm_mosaic_candidates
	select_outgroup_panel
	read_mosaic_catalogue
);

sub pair_key {
	my ($left, $right) = @_;
	return $left le $right ? "$left\t$right" : "$right\t$left";
}

sub _members {
	my ($value) = @_;
	return () unless defined $value;
	my @members = ref($value) eq 'ARRAY' ? @{$value} : split /,/, $value;
	for (@members) {
		s/^>//;
		s/^\s+|\s+$//g;
	}
	return grep { length } @members;
}

sub discover_mosaic_candidates {
	my ($records, $cluster_members, $sequences, $options) = @_;
	$options ||= {};
	my $minimum_length_ratio = $options->{minimum_length_ratio} // 0.80;
	my $maximum_sample_overlap_fraction =
		$options->{maximum_sample_overlap_fraction} // 0.15;
	my $maximum_overlap_samples = $options->{maximum_overlap_samples} // 1;
	my %by_mgs_cog;
	for my $record (@{$records || []}) {
		next unless defined($record->{gene}) && defined($record->{mgs}) && defined($record->{cog});
		next if $record->{cog} eq '' || $record->{cog} eq '-';
		push @{$by_mgs_cog{$record->{mgs}}{$record->{cog}}}, $record;
	}

	my @candidates;
	for my $mgs (sort keys %by_mgs_cog) {
		for my $cog (sort keys %{$by_mgs_cog{$mgs}}) {
			my @genes = sort {
				($a->{rank} // 0) <=> ($b->{rank} // 0) || $a->{gene} cmp $b->{gene}
			} @{$by_mgs_cog{$mgs}{$cog}};
			for my $i (0 .. $#genes - 1) {
				for my $j ($i + 1 .. $#genes) {
					my ($left, $right) = ($genes[$i]{gene}, $genes[$j]{gene});
					next unless defined($sequences->{$left}) && defined($sequences->{$right});
					my $left_length = length($sequences->{$left});
					my $right_length = length($sequences->{$right});
					next unless $left_length && $right_length;
					my $length_ratio = $left_length < $right_length
						? $left_length / $right_length : $right_length / $left_length;
					next if $length_ratio < $minimum_length_ratio;
					my %left_samples;
					for my $member (_members($cluster_members->{$left})) {
						my ($sample) = split /__/, $member, 2;
						$left_samples{$sample} = 1 if defined($sample) && length($sample);
					}
					my %right_samples;
					for my $member (_members($cluster_members->{$right})) {
						my ($sample) = split /__/, $member, 2;
						$right_samples{$sample} = 1 if defined($sample) && length($sample);
					}
					my $overlap = grep { exists $right_samples{$_} } keys %left_samples;
					next unless keys(%left_samples) && keys(%right_samples);
					my $smaller_sample_set = scalar(keys %left_samples) < scalar(keys %right_samples)
						? scalar(keys %left_samples) : scalar(keys %right_samples);
					my $overlap_fraction = $overlap / $smaller_sample_set;
					next if $overlap > $maximum_overlap_samples
						|| $overlap_fraction > $maximum_sample_overlap_fraction;
					push @candidates, {
						mgs => $mgs, cog => $cog, left => $left, right => $right,
						length_ratio => $length_ratio,
						left_samples => scalar(keys %left_samples),
						right_samples => scalar(keys %right_samples),
						overlap_samples => $overlap,
						overlap_fraction => $overlap_fraction,
					};
				}
			}
		}
	}
	return \@candidates;
}

sub read_paf_hits {
	my ($path) = @_;
	my %hits;
	my ($fh) = gzipopen($path, "mosaic whole-catalogue alignments", 1);
	my $line_number = 0;
	while (my $line = <$fh>) {
		$line_number++;
		$line =~ s/[\r\n]+$//;
		next if $line eq '' || $line =~ /^#/;
		my @fields = split /\t/, $line;
		die "Malformed PAF row $line_number in $path\n" unless @fields >= 12;
		my ($query, $query_length, $query_start, $query_end, $strand,
			$target, $target_length, $target_start, $target_end,
			$matches, $alignment_length, $mapq) = @fields[0 .. 11];
		die "Invalid numeric PAF row $line_number in $path\n"
			unless grep(!/^\d+$/, ($query_length, $query_start, $query_end,
				$target_length, $target_start, $target_end, $matches,
				$alignment_length, $mapq)) == 0
				&& $query_length > 0 && $target_length > 0 && $alignment_length > 0;
		push @{$hits{$query}}, {
			query => $query, target => $target, strand => $strand,
			query_length => 0 + $query_length, target_length => 0 + $target_length,
			query_coverage => ($query_end - $query_start) / $query_length,
			target_coverage => ($target_end - $target_start) / $target_length,
			identity => $matches / $alignment_length,
			matches => 0 + $matches, alignment_length => 0 + $alignment_length,
			mapq => 0 + $mapq,
		};
	}
	close $fh or die "Cannot close PAF $path: $!\n";
	return \%hits;
}

sub _directional_partner {
	my ($query, $partner, $hits, $options) = @_;
	my $minimum_identity = $options->{minimum_identity} // 0.90;
	my $minimum_query_coverage = $options->{minimum_query_coverage} // 0.80;
	my $minimum_target_coverage = $options->{minimum_target_coverage} // 0.80;
	my $minimum_score_margin = $options->{minimum_score_margin} // 0.02;
	my @raw_usable = grep {
		$_->{target} ne $query
			&& $_->{identity} >= $minimum_identity
			&& $_->{query_coverage} >= $minimum_query_coverage
			&& $_->{target_coverage} >= $minimum_target_coverage
	} @{$hits->{$query} || []};
	my %best_by_target;
	for my $hit (@raw_usable) {
		my $current = $best_by_target{$hit->{target}};
		if (!$current
			|| $hit->{matches} > $current->{matches}
			|| ($hit->{matches} == $current->{matches}
				&& $hit->{identity} > $current->{identity})
			|| ($hit->{matches} == $current->{matches}
				&& $hit->{identity} == $current->{identity}
				&& $hit->{query_coverage} > $current->{query_coverage})) {
			$best_by_target{$hit->{target}} = $hit;
		}
	}
	my @usable = values %best_by_target;
	@usable = sort {
		$b->{matches} <=> $a->{matches}
			|| $b->{identity} <=> $a->{identity}
			|| $b->{query_coverage} <=> $a->{query_coverage}
			|| $a->{target} cmp $b->{target}
	} @usable;
	my ($partner_hit) = grep { $_->{target} eq $partner } @usable;
	return (undef, 'partner_alignment_missing') unless $partner_hit;
	my ($best) = @usable;
	return (undef, 'partner_not_best') unless $best->{target} eq $partner;
	if (@usable > 1) {
		my $runner = $usable[1];
		my $required = $partner_hit->{matches} * (1 - $minimum_score_margin);
		return (undef, 'nonunique_catalogue_hit') if $runner->{matches} >= $required;
	}
	return ($partner_hit, 'confirmed');
}

sub confirm_mosaic_candidates {
	my ($candidates, $hits, $options) = @_;
	$options ||= {};
	my (@confirmed, @rejected);
	for my $candidate (@{$candidates || []}) {
		my ($left_hit, $left_reason) =
			_directional_partner($candidate->{left}, $candidate->{right}, $hits, $options);
		my ($right_hit, $right_reason) =
			_directional_partner($candidate->{right}, $candidate->{left}, $hits, $options);
		unless ($left_hit && $right_hit) {
			push @rejected, {
				%{$candidate},
				reason => $left_hit ? $right_reason : $left_reason,
			};
			next;
		}
		push @confirmed, {
			%{$candidate},
			identity => ($left_hit->{identity} + $right_hit->{identity}) / 2,
			query_coverage => ($left_hit->{query_coverage} + $right_hit->{query_coverage}) / 2,
			target_coverage => ($left_hit->{target_coverage} + $right_hit->{target_coverage}) / 2,
		};
	}
	return (\@confirmed, \@rejected);
}

sub _median {
	my @values = sort { $a <=> $b } @_;
	return 0 unless @values;
	my $middle = int(@values / 2);
	return @values % 2 ? $values[$middle]
		: ($values[$middle - 1] + $values[$middle]) / 2;
}

sub select_outgroup_panel {
	my ($records, $hits, $options) = @_;
	$options ||= {};
	my $minimum_identity = $options->{minimum_identity} // 0.80;
	my $maximum_identity = $options->{maximum_identity} // 0.95;
	my $minimum_coverage = $options->{minimum_coverage} // 0.75;
	my $minimum_loci = $options->{minimum_loci} // 10;
	my $target_identity = $options->{target_identity} // 0.88;
	my %gene_record = map { $_->{gene} => $_ } @{$records || []};
	my (%best_per_query_target_mgs, %aggregate);

	for my $query (keys %gene_record) {
		my $source = $gene_record{$query};
		for my $hit (@{$hits->{$query} || []}) {
			my $target = $gene_record{$hit->{target}} or next;
			next if $target->{mgs} eq $source->{mgs};
			next unless $target->{cog} eq $source->{cog};
			next if $hit->{identity} < $minimum_identity || $hit->{identity} > $maximum_identity;
			next if $hit->{query_coverage} < $minimum_coverage
				|| $hit->{target_coverage} < $minimum_coverage;
			my $slot = \$best_per_query_target_mgs{$query}{$target->{mgs}};
			if (!defined($$slot)
				|| $hit->{matches} > $$slot->{matches}
				|| ($hit->{matches} == $$slot->{matches}
					&& $hit->{target} lt $$slot->{target})) {
				$$slot = {%{$hit}, source_mgs => $source->{mgs},
					target_mgs => $target->{mgs}, cog => $source->{cog}};
			}
		}
	}

	for my $query (keys %best_per_query_target_mgs) {
		for my $target_mgs (keys %{$best_per_query_target_mgs{$query}}) {
			my $hit = $best_per_query_target_mgs{$query}{$target_mgs};
			my $key = "$hit->{source_mgs}\t$target_mgs";
			push @{$aggregate{$key}{hits}}, $hit;
			$aggregate{$key}{source_mgs} = $hit->{source_mgs};
			$aggregate{$key}{target_mgs} = $target_mgs;
		}
	}

	my (%preferred, %gene_map);
	for my $source_mgs (sort { $a cmp $b } map { $_->{mgs} } @{$records || []}) {
		next if exists $preferred{$source_mgs};
		my @panels;
		for my $entry (values %aggregate) {
			next unless $entry->{source_mgs} eq $source_mgs;
			my %source_genes = map { $_->{query} => 1 } @{$entry->{hits}};
			my $loci = scalar keys %source_genes;
			next if $loci < $minimum_loci;
			my $median_identity = _median(map { $_->{identity} } @{$entry->{hits}});
			push @panels, {
				%{$entry}, loci => $loci, median_identity => $median_identity,
				distance_from_target => abs($median_identity - $target_identity),
			};
		}
		@panels = sort {
			$b->{loci} <=> $a->{loci}
				|| $a->{distance_from_target} <=> $b->{distance_from_target}
				|| $a->{target_mgs} cmp $b->{target_mgs}
		} @panels;
		next unless @panels;
		my $best = $panels[0];
		$preferred{$source_mgs} = $best;
		for my $hit (@{$best->{hits}}) {
			$gene_map{$source_mgs}{$hit->{query}} = $hit;
		}
	}
	return (\%preferred, \%gene_map);
}

sub read_mosaic_catalogue {
	my ($path) = @_;
	my (%pairs, %outgroups, %outgroup_genes);
	return (\%pairs, \%outgroups, \%outgroup_genes)
		unless defined($path) && length($path) && -s $path;
	my ($fh) = gzipopen($path, "confirmed mosaic catalogue", 1);
	my $line_number = 0;
	while (my $line = <$fh>) {
		$line_number++;
		$line =~ s/[\r\n]+$//;
		next if $line eq '' || $line =~ /^#/;
		my @fields = split /\t/, $line, -1;
		if ($fields[0] eq 'MOSAIC') {
			die "Malformed MOSAIC row $line_number in $path\n" unless @fields >= 8;
			$pairs{pair_key($fields[3], $fields[4])} = 1;
		} elsif ($fields[0] eq 'OUTGROUP') {
			die "Malformed OUTGROUP row $line_number in $path\n" unless @fields >= 5;
			$outgroups{$fields[1]} = $fields[2];
		} elsif ($fields[0] eq 'OUTGROUP_GENE') {
			die "Malformed OUTGROUP_GENE row $line_number in $path\n" unless @fields >= 7;
			$outgroup_genes{$fields[1]}{$fields[3]} = $fields[4];
		} else {
			die "Unknown mosaic catalogue row type '$fields[0]' at $path line $line_number\n";
		}
	}
	close $fh or die "Cannot close mosaic catalogue $path: $!\n";
	return (\%pairs, \%outgroups, \%outgroup_genes);
}

1;
