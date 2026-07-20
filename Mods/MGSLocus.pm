package Mods::MGSLocus;

use strict;
use warnings;

use Exporter qw(import);

our @EXPORT_OK = qw(
	build_locus_groups
	choose_locus_candidate
	protein_kmer_similarity
	robust_depth_mask
);

sub _member_parts {
	my ($member) = @_;
	return unless defined $member;
	$member =~ s/^>//;
	my ($sample, $gene) = split /__/, $member, 2;
	return unless defined($sample) && length($sample) && defined($gene) && length($gene);
	my ($contig, $position) = $gene =~ /^(.*)_([0-9]+)$/;
	return ($sample, $contig, defined($position) ? 0 + $position : undef, $member);
}

sub _members {
	my ($value) = @_;
	return () unless defined $value;
	my @members = ref($value) eq 'ARRAY' ? @{$value} : split(/,/, $value);
	for (@members) {
		s/^>//;
		s/^\s+|\s+$//g;
	}
	return grep { length } @members;
}

sub _set_similarity {
	my ($left, $right) = @_;
	return 0 unless $left && $right && keys(%{$left}) && keys(%{$right});
	my $intersection = 0;
	$intersection++ for grep { exists $right->{$_} } keys %{$left};
	my %union = (%{$left}, %{$right});
	return $intersection / scalar(keys %union);
}

sub protein_kmer_similarity {
	my ($left, $right, $k) = @_;
	$k ||= 4;
	return 0 unless defined($left) && defined($right);
	$left = uc($left);
	$right = uc($right);
	$left =~ s/[^A-Z]//g;
	$right =~ s/[^A-Z]//g;
	return 0 if length($left) < $k || length($right) < $k;

	my (%left_kmers, %right_kmers);
	$left_kmers{substr($left, $_, $k)} = 1 for 0 .. length($left) - $k;
	$right_kmers{substr($right, $_, $k)} = 1 for 0 .. length($right) - $k;
	my $intersection = 0;
	$intersection++ for grep { exists $right_kmers{$_} } keys %left_kmers;
	my $denominator = scalar(keys %left_kmers) + scalar(keys %right_kmers);
	return $denominator ? (2 * $intersection / $denominator) : 0;
}

sub _find {
	my ($parent, $node) = @_;
	$parent->{$node} = _find($parent, $parent->{$node}) if $parent->{$node} ne $node;
	return $parent->{$node};
}

sub build_locus_groups {
	my ($records, $cluster_members, $proteins, $options) = @_;
	$options ||= {};
	my $context_distance = $options->{context_distance} // 5;
	my $min_length_ratio = $options->{min_length_ratio} // 0.75;
	my $min_sequence_with_context = $options->{min_sequence_with_context} // 0.55;
	my $min_sequence_without_context = $options->{min_sequence_without_context} // 0.72;
	my $min_context_similarity = $options->{min_context_similarity} // 0.25;

	my (%record_by_gene, %sample_set, %member_seed, %positions);
	for my $record (@{$records || []}) {
		my $gene = $record->{gene};
		next unless defined($gene) && length($gene);
		$record_by_gene{$gene} = $record;
		for my $member (_members($cluster_members->{$gene})) {
			my ($sample, $contig, $position, $clean_member) = _member_parts($member);
			next unless defined $sample;
			$sample_set{$gene}{$sample} = 1;
			$member_seed{$clean_member} = $gene;
			next unless defined($contig) && defined($position);
			push @{$positions{$sample}{$contig}}, {
				position => $position,
				member   => $clean_member,
				gene     => $gene,
				mgs      => $record->{mgs},
				cog      => $record->{cog},
			};
		}
	}

	my (%gene_context, %member_context);
	for my $sample (keys %positions) {
		for my $contig (keys %{$positions{$sample}}) {
			my @entries = sort {
				$a->{position} <=> $b->{position} || $a->{gene} cmp $b->{gene}
			} @{$positions{$sample}{$contig}};
			for my $i (0 .. $#entries) {
				for my $j (0 .. $#entries) {
					next if $i == $j;
					next if abs($entries[$i]{position} - $entries[$j]{position}) > $context_distance;
					my $token = join('|', $entries[$j]{mgs}, $entries[$j]{cog});
					$gene_context{$entries[$i]{gene}}{$token}++;
					$member_context{$entries[$i]{member}}{$token}++;
				}
			}
		}
	}

	my %by_mgs_cog;
	for my $record (@{$records || []}) {
		push @{$by_mgs_cog{$record->{mgs}}{$record->{cog}}}, $record;
	}

	my (@groups, %gene_to_locus, %locus_context, %locus_by_id);
	my $merged_seeds = 0;
	for my $mgs (sort keys %by_mgs_cog) {
		for my $cog (sort keys %{$by_mgs_cog{$mgs}}) {
			my @seeds = sort {
				$a->{rank} <=> $b->{rank} || $a->{gene} cmp $b->{gene}
			} @{$by_mgs_cog{$mgs}{$cog}};
			my (%parent, %component_samples);
			for my $seed (@seeds) {
				$parent{$seed->{gene}} = $seed->{gene};
				$component_samples{$seed->{gene}} = { %{$sample_set{$seed->{gene}} || {}} };
			}

			my @edges;
			for my $i (0 .. $#seeds - 1) {
				for my $j ($i + 1 .. $#seeds) {
					my ($left, $right) = ($seeds[$i]{gene}, $seeds[$j]{gene});
					my $cooccurs = grep { exists $sample_set{$right}{$_} } keys %{$sample_set{$left} || {}};
					next if $cooccurs;
					my ($left_seq, $right_seq) = ($proteins->{$left}, $proteins->{$right});
					next unless defined($left_seq) && defined($right_seq) && length($left_seq) && length($right_seq);
					my $length_ratio = length($left_seq) < length($right_seq)
						? length($left_seq) / length($right_seq)
						: length($right_seq) / length($left_seq);
					next if $length_ratio < $min_length_ratio;
					my $sequence_score = protein_kmer_similarity($left_seq, $right_seq);
					my $context_score = _set_similarity($gene_context{$left}, $gene_context{$right});
					my $has_context = scalar(keys %{$gene_context{$left} || {}}) >= 2
						&& scalar(keys %{$gene_context{$right} || {}}) >= 2;
					next if $has_context
						? ($sequence_score < $min_sequence_with_context || $context_score < $min_context_similarity)
						: ($sequence_score < $min_sequence_without_context);
					push @edges, {
						left => $left, right => $right,
						score => $sequence_score + 0.15 * $context_score,
					};
				}
			}

			for my $edge (sort {
				$b->{score} <=> $a->{score} || $a->{left} cmp $b->{left} || $a->{right} cmp $b->{right}
			} @edges) {
				my $left_root = _find(\%parent, $edge->{left});
				my $right_root = _find(\%parent, $edge->{right});
				next if $left_root eq $right_root;
				my $overlap = grep { exists $component_samples{$right_root}{$_} } keys %{$component_samples{$left_root}};
				next if $overlap;
				$parent{$right_root} = $left_root;
				$component_samples{$left_root}{$_} = 1 for keys %{$component_samples{$right_root}};
				delete $component_samples{$right_root};
				$merged_seeds++;
			}

			my %components;
			push @{$components{_find(\%parent, $_->{gene})}}, $_ for @seeds;
			for my $component (values %components) {
				my @ordered = sort {
					$a->{rank} <=> $b->{rank} || $a->{gene} cmp $b->{gene}
				} @{$component};
				my $primary = $ordered[0]{gene};
				my $locus_id = join('|', $mgs, $cog, $primary);
				my %context;
				for my $seed (@ordered) {
					$gene_to_locus{$seed->{gene}} = $locus_id;
					$context{$_} += $gene_context{$seed->{gene}}{$_} for keys %{$gene_context{$seed->{gene}} || {}};
				}
				$locus_context{$locus_id} = \%context;
				my $group = {
					mgs => $mgs, cog => $cog, primary_gene => $primary,
					genes => [map { $_->{gene} } @ordered],
					rank => $ordered[0]{rank}, locus_id => $locus_id,
				};
				push @groups, $group;
				$locus_by_id{$locus_id} = $group;
			}
		}
	}

	@groups = sort {
		$a->{mgs} cmp $b->{mgs} || $a->{rank} <=> $b->{rank} || $a->{locus_id} cmp $b->{locus_id}
	} @groups;
	return {
		groups => \@groups,
		gene_to_locus => \%gene_to_locus,
		locus_by_id => \%locus_by_id,
		member_to_seed => \%member_seed,
		member_context => \%member_context,
		locus_context => \%locus_context,
		merged_seeds => $merged_seeds,
	};
}

sub choose_locus_candidate {
	my ($candidates, $seed_proteins, $locus_context, $options) = @_;
	$options ||= {};
	my $sequence_margin = $options->{sequence_margin} // 0.08;
	my $context_margin = $options->{context_margin} // 0.25;
	my $depth_ratio = $options->{depth_ratio} // 2;
	return { status => 'missing' } unless $candidates && @{$candidates};

	my @scored;
	for my $candidate (@{$candidates}) {
		my $sequence_score = 0;
		for my $seed_sequence (values %{$seed_proteins || {}}) {
			my $score = protein_kmer_similarity($candidate->{protein}, $seed_sequence);
			$sequence_score = $score if $score > $sequence_score;
		}
		push @scored, {
			%{$candidate},
			sequence_score => $sequence_score,
			context_score => _set_similarity($candidate->{context} || {}, $locus_context || {}),
		};
	}
	@scored = sort {
		$b->{sequence_score} <=> $a->{sequence_score}
			|| $b->{context_score} <=> $a->{context_score}
			|| ($b->{depth} // 0) <=> ($a->{depth} // 0)
			|| $a->{id} cmp $b->{id}
	} @scored;
	return { status => 'selected', candidate => $scored[0], reason => 'unique' } if @scored == 1;

	my ($best, $runner_up) = @scored[0, 1];
	if ($best->{sequence_score} - $runner_up->{sequence_score} >= $sequence_margin) {
		return { status => 'selected', candidate => $best, reason => 'sequence' };
	}
	if ($best->{context_score} >= 0.5
		&& $best->{context_score} - $runner_up->{context_score} >= $context_margin) {
		return { status => 'selected', candidate => $best, reason => 'context' };
	}
	my %seeds = map { ($_->{seed} // '') => 1 } @scored;
	my $runner_depth = $runner_up->{depth} // 0;
	if (keys(%seeds) == 1 && $runner_depth > 0
		&& ($best->{depth} // 0) / $runner_depth >= $depth_ratio) {
		return { status => 'selected', candidate => $best, reason => 'depth' };
	}
	return { status => 'ambiguous', candidates => \@scored };
}

sub _median {
	my (@values) = sort { $a <=> $b } @_;
	return 0 unless @values;
	my $middle = int(@values / 2);
	return @values % 2 ? $values[$middle] : ($values[$middle - 1] + $values[$middle]) / 2;
}

sub robust_depth_mask {
	my ($depths, $options) = @_;
	$options ||= {};
	my $minimum_count = $options->{minimum_count} // 8;
	my $maximum_modified_z = $options->{maximum_modified_z} // 3.5;
	my $minimum_fold = $options->{minimum_fold} // (1 / 3);
	my $maximum_fold = $options->{maximum_fold} // 3;
	my @mask = (1) x scalar(@{$depths || []});
	return \@mask unless $depths && @{$depths} >= $minimum_count;
	my $median = _median(@{$depths});
	return \@mask unless $median > 0;
	my @deviations = map { abs($_ - $median) } @{$depths};
	my $mad = _median(@deviations);
	for my $i (0 .. $#{$depths}) {
		my $fold = $depths->[$i] / $median;
		my $outside_fold = $fold < $minimum_fold || $fold > $maximum_fold;
		next unless $outside_fold;
		my $modified_z = $mad > 0 ? 0.6745 * abs($depths->[$i] - $median) / $mad : 1e9;
		$mask[$i] = 0 if $modified_z > $maximum_modified_z;
	}
	return \@mask;
}

1;
