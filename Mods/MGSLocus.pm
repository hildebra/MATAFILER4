package Mods::MGSLocus;

use strict;
use warnings;

use Exporter qw(import);

our @EXPORT_OK = qw(
	build_locus_groups
	choose_locus_candidate
	member_context_map
	accumulate_locus_context
	merge_candidate_seeds
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

sub _pair_key {
	my ($left, $right) = @_;
	return $left le $right ? "$left\t$right" : "$right\t$left";
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

#Scan the catalogue members of every ranked seed once and derive the three
#summaries the grouping needs.  Seed-level summaries (sample_set, gene_context)
#aggregate across every sample that carries the seed, so they are only correct
#when $cluster_members covers the whole catalogue.  Member-level summaries
#(member_seed, member_context) describe a member's own contig neighbourhood and
#are therefore identical whether they are derived from the complete catalogue or
#from a shard holding every member of the requested samples.
sub _scan_members {
	my ($records, $cluster_members, $options, $accumulator) = @_;
	$options ||= {};
	# Seed-level summaries are only ever read for seeds that could merge. When the
	# caller names that set, everything else is scanned for its neighbour tokens
	# but never stored, which is what keeps a catalogue-wide scan affordable.
	my $context_seeds = $options->{context_seeds};
	# Sample-level summaries are additive over disjoint sample sets, because a
	# contig belongs to exactly one sample and no context ever crosses samples.
	# An accumulator therefore lets a caller stream the catalogue in slices.
	$accumulator ||= {};
	my $sample_set = $accumulator->{sample_set} ||= {};
	my $gene_context = $accumulator->{gene_context} ||= {};
	my $context_distance = $options->{context_distance} // 5;
	my $include_member_to_seed = $options->{include_member_to_seed} ? 1 : 0;
	my $include_member_context = $options->{include_member_context} ? 1 : 0;
	my $include_gene_context = $options->{include_gene_context} ? 1 : 0;
	my $want_positions = $include_member_context || $include_gene_context;
	# Callers that hold a throwaway catalogue-wide membership map can release each
	# comma-joined member string as it is consumed instead of keeping the whole
	# map resident beside the summaries built from it.  Ranked records name each
	# seed at most once, so a consumed entry is never needed again.
	my $consume = $options->{consume_cluster_members} ? 1 : 0;

	my (%member_seed, %positions);
	for my $record (@{$records || []}) {
		my $gene = $record->{gene};
		next unless defined($gene) && length($gene);
		my $wanted = !$context_seeds || $context_seeds->{$gene} ? 1 : 0;
		for my $member (_members($consume
				? delete($cluster_members->{$gene}) : $cluster_members->{$gene})) {
			my ($sample, $contig, $position, $clean_member) = _member_parts($member);
			next unless defined $sample;
			$sample_set->{$gene}{$sample} = 1 if $wanted;
			$member_seed{$clean_member} = $gene if $include_member_to_seed;
			next unless $want_positions;
			next unless defined($contig) && defined($position);
			# Entries are arrays rather than hashes: at catalogue scale this is the
			# single largest transient structure, and a five-key hash per member
			# costs several times what these four fields do.
			push @{$positions{$sample}{$contig}}, [
				$position, $gene, $record->{mgs}, $record->{cog},
				$include_member_context ? $clean_member : undef,
			];
		}
	}

	my %member_context;
	for my $sample (keys %positions) {
		for my $contig (keys %{$positions{$sample}}) {
			my @entries = sort {
				$a->[0] <=> $b->[0] || $a->[1] cmp $b->[1]
			} @{$positions{$sample}{$contig}};
			for my $i (0 .. $#entries) {
				my $recordGeneContext = $include_gene_context
					&& (!$context_seeds || $context_seeds->{$entries[$i][1]});
				next unless $recordGeneContext || $include_member_context;
				for my $j (0 .. $#entries) {
					next if $i == $j;
					next if abs($entries[$i][0] - $entries[$j][0]) > $context_distance;
					my $token = join('|', $entries[$j][2], $entries[$j][3]);
					$gene_context->{$entries[$i][1]}{$token}++ if $recordGeneContext;
					$member_context{$entries[$i][4]}{$token}++ if $include_member_context;
				}
			}
		}
	}
	# Position entries can dominate peak memory for large catalogues and are no
	# longer needed after their compact context summaries have been built.
	%positions = ();
	return ($sample_set, \%member_seed, $gene_context, \%member_context);
}

#Merge into $accumulator the seed-level summaries for one slice of the catalogue.
#Callers stream disjoint sample slices through this so the transient per-member
#structures stay bounded by one slice instead of the whole catalogue.
sub accumulate_locus_context {
	my ($accumulator, $records, $cluster_members, $options) = @_;
	die "Locus-context accumulator must be a hash reference\n"
		unless ref($accumulator) eq 'HASH';
	$options ||= {};
	_scan_members($records, $cluster_members, {
		context_distance => $options->{context_distance} // 5,
		include_member_to_seed => 0,
		include_member_context => 0,
		include_gene_context => 1,
		consume_cluster_members => $options->{consume_cluster_members} ? 1 : 0,
		context_seeds => $options->{context_seeds},
	}, $accumulator);
	return $accumulator;
}

#Seeds whose catalogue-wide summaries can change a locus boundary: a seed alone
#in its MGS/COG has nothing to merge with, and when a confirmed-pair allowlist is
#given only seeds named in such a pair are ever compared. Everything else needs
#no sample set and no gene context at all.
sub merge_candidate_seeds {
	my ($records, $allowed_merge_pairs) = @_;
	my %byGroup;
	for my $record (@{$records || []}) {
		my $gene = $record->{gene};
		next unless defined($gene) && length($gene);
		push @{$byGroup{$record->{mgs}}{$record->{cog}}}, $gene;
	}
	my %candidates;
	for my $mgs (keys %byGroup) {
		for my $cog (keys %{$byGroup{$mgs}}) {
			my @seeds = @{$byGroup{$mgs}{$cog}};
			next if @seeds < 2;
			unless (defined $allowed_merge_pairs) {
				$candidates{$_} = 1 for @seeds;
				next;
			}
			for my $i (0 .. $#seeds - 1) {
				for my $j ($i + 1 .. $#seeds) {
					next unless $allowed_merge_pairs->{_pair_key($seeds[$i], $seeds[$j])};
					$candidates{$seeds[$i]} = 1;
					$candidates{$seeds[$j]} = 1;
				}
			}
		}
	}
	return \%candidates;
}

#Member contexts for one sample slice, without the catalogue-wide seed summaries
#that grouping needs.  Split extraction workers use this to rebuild the context
#of their own members around a locus model published by the parent.
sub member_context_map {
	my ($records, $cluster_members, $options) = @_;
	$options ||= {};
	my (undef, undef, undef, $member_context) =
		_scan_members($records, $cluster_members, {
			context_distance => $options->{context_distance} // 5,
			include_member_to_seed => 0,
			include_member_context => 1,
			include_gene_context => 0,
		});
	return $member_context;
}

sub build_locus_groups {
	my ($records, $cluster_members, $proteins, $options) = @_;
	$options ||= {};
	my $include_member_context = exists($options->{include_member_context})
		? $options->{include_member_context} : 1;
	my $include_member_to_seed = exists($options->{include_member_to_seed})
		? $options->{include_member_to_seed} : 1;
	my $include_gene_to_locus = exists($options->{include_gene_to_locus})
		? $options->{include_gene_to_locus} : 1;
	my $context_distance = $options->{context_distance} // 5;
	my $min_length_ratio = $options->{min_length_ratio} // 0.75;
	my $min_sequence_with_context = $options->{min_sequence_with_context} // 0.55;
	my $min_sequence_without_context = $options->{min_sequence_without_context} // 0.72;
	my $min_context_similarity = $options->{min_context_similarity} // 0.25;
	my $allowed_merge_pairs = $options->{allowed_merge_pairs};
	my $require_complete_linkage = $options->{require_complete_linkage} // 0;
	my $allow_confirmed_cooccurrence =
		$options->{allow_confirmed_cooccurrence} // 0;
	my $consume_cluster_members = $options->{consume_cluster_members} ? 1 : 0;

	my ($sample_set, $member_seed, $gene_context, $member_context);
	if (my $precomputed = $options->{precomputed_context}) {
		# A caller that streamed the catalogue in slices has already derived the
		# seed-level summaries; $cluster_members is then not consulted at all.
		# Grouping releases each seed's entry as it publishes a locus, so take a
		# shallow copy: the per-seed hashes are shared, but the caller's own index
		# survives. Restricted to merge candidates this costs almost nothing.
		$sample_set = { %{$precomputed->{sample_set} || {}} };
		$gene_context = { %{$precomputed->{gene_context} || {}} };
		($member_seed, $member_context) = ({}, {});
	} else {
		($sample_set, $member_seed, $gene_context, $member_context) =
			_scan_members($records, $cluster_members, {
				context_distance => $context_distance,
				include_member_to_seed => $include_member_to_seed,
				include_member_context => $include_member_context,
				include_gene_context => 1,
				consume_cluster_members => $consume_cluster_members,
				context_seeds => $options->{context_seeds},
			});
	}

	my %by_mgs_cog;
	for my $record (@{$records || []}) {
		push @{$by_mgs_cog{$record->{mgs}}{$record->{cog}}}, $record;
	}

	my (@groups, %gene_to_locus, %locus_context, %locus_by_id);
	my ($merged_seeds, $incomplete_linkage_rejections) = (0, 0);
	for my $mgs (sort keys %by_mgs_cog) {
		for my $cog (sort keys %{$by_mgs_cog{$mgs}}) {
			my @seeds = sort {
				$a->{rank} <=> $b->{rank} || $a->{gene} cmp $b->{gene}
			} @{$by_mgs_cog{$mgs}{$cog}};
			my (%parent, %component_samples);
			for my $seed (@seeds) {
				$parent{$seed->{gene}} = $seed->{gene};
				$component_samples{$seed->{gene}} = { %{$sample_set->{$seed->{gene}} || {}} };
			}

			my @edges;
			my %edge_ok;
			for my $i (0 .. $#seeds - 1) {
				for my $j ($i + 1 .. $#seeds) {
					my ($left, $right) = ($seeds[$i]{gene}, $seeds[$j]{gene});
					my $pair_key = _pair_key($left, $right);
					my $pair_is_confirmed = defined($allowed_merge_pairs)
						&& $allowed_merge_pairs->{$pair_key};
					next if defined($allowed_merge_pairs)
						&& !$pair_is_confirmed;
					my $cooccurs = grep { exists $sample_set->{$right}{$_} } keys %{$sample_set->{$left} || {}};
					next if $cooccurs
						&& !($allow_confirmed_cooccurrence && $pair_is_confirmed);
					my ($left_seq, $right_seq) = ($proteins->{$left}, $proteins->{$right});
					next unless defined($left_seq) && defined($right_seq) && length($left_seq) && length($right_seq);
					my $length_ratio = length($left_seq) < length($right_seq)
						? length($left_seq) / length($right_seq)
						: length($right_seq) / length($left_seq);
					next if $length_ratio < $min_length_ratio;
					my $sequence_score = protein_kmer_similarity($left_seq, $right_seq);
					my $context_score = _set_similarity($gene_context->{$left}, $gene_context->{$right});
					my $has_context = scalar(keys %{$gene_context->{$left} || {}}) >= 2
						&& scalar(keys %{$gene_context->{$right} || {}}) >= 2;
					next if $has_context
						? ($sequence_score < $min_sequence_with_context || $context_score < $min_context_similarity)
						: ($sequence_score < $min_sequence_without_context);
					push @edges, {
						left => $left, right => $right,
						score => $sequence_score + 0.15 * $context_score,
					};
					$edge_ok{$pair_key} = 1;
				}
			}

			my %component_genes = map { $_->{gene} => { $_->{gene} => 1 } } @seeds;
			for my $edge (sort {
				$b->{score} <=> $a->{score} || $a->{left} cmp $b->{left} || $a->{right} cmp $b->{right}
			} @edges) {
				my $left_root = _find(\%parent, $edge->{left});
				my $right_root = _find(\%parent, $edge->{right});
				next if $left_root eq $right_root;
				my $overlap = grep { exists $component_samples{$right_root}{$_} } keys %{$component_samples{$left_root}};
				next if $overlap && !$allow_confirmed_cooccurrence;
				if ($require_complete_linkage) {
					my $all_compatible = 1;
					OUTER:
					for my $left_gene (keys %{$component_genes{$left_root}}) {
						for my $right_gene (keys %{$component_genes{$right_root}}) {
							unless ($edge_ok{_pair_key($left_gene, $right_gene)}) {
								$all_compatible = 0;
								last OUTER;
							}
						}
					}
					unless ($all_compatible) {
						$incomplete_linkage_rejections++;
						next;
					}
				}
				$parent{$right_root} = $left_root;
				$component_samples{$left_root}{$_} = 1 for keys %{$component_samples{$right_root}};
				delete $component_samples{$right_root};
				$component_genes{$left_root}{$_} = 1 for keys %{$component_genes{$right_root}};
				delete $component_genes{$right_root};
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
					$gene_to_locus{$seed->{gene}} = $locus_id if $include_gene_to_locus;
					$context{$_} += $gene_context->{$seed->{gene}}{$_} for keys %{$gene_context->{$seed->{gene}} || {}};
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
			# All comparisons for this MGS/COG are complete.  Release the
			# expanded per-seed inputs as the compact groups are published.
			for my $seed (@seeds) {
				delete $sample_set->{$seed->{gene}};
				delete $gene_context->{$seed->{gene}};
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
		member_to_seed => $member_seed,
		member_context => $member_context,
		locus_context => \%locus_context,
		merged_seeds => $merged_seeds,
		incomplete_linkage_rejections => $incomplete_linkage_rejections,
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
