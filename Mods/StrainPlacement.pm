package Mods::StrainPlacement;

use strict;
use warnings;
use Exporter qw(import);
use JSON::PP qw(decode_json);
use Mods::GenoMetaAss qw(fileGZe gzipopen);

our @EXPORT_OK = qw(
	read_sample_qc
	split_strict_backbone
	read_epa_jplace
	filter_epa_placement_outliers
	write_epa_placed_tree
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
sub _valid_mask {
	my ($sequence, $is_aa) = @_;
	my $mask = uc($sequence // '');
	if ($is_aa) {
		$mask =~ tr/ACDEFGHIKLMNPQRSTVWY/\x01/;
	} else {
		$mask =~ tr/ACGTU/\x01/;
	}
	$mask =~ tr/\x01/\x00/c;
	return $mask;
}

sub _dna_state_mask {
	my ($sequence) = @_;
	my $mask = uc($sequence // '');
	$mask =~ tr/ACGTU/\x01\x02\x04\x08\x08/;
	$mask =~ tr/\x01\x02\x04\x08/\x00/c;
	return $mask;
}

sub _mask_count {
	return $_[0] =~ tr/\x01/\x01/;
}

sub _alignment_locus_ranges {
	my ($partition_file, $alignment_length) = @_;
	return [[0, $alignment_length]]
		unless defined($partition_file) && length($partition_file);
	die "Strict-backbone overlap requires the alignment partition file $partition_file\n"
		unless -s $partition_file;
	open my $fh, '<', $partition_file
		or die "Cannot read alignment partition file $partition_file: $!\n";
	my @ranges;
	while (my $line = <$fh>) {
		while ($line =~ /(\d+)\s*-\s*(\d+)/g) {
			my ($start, $end) = ($1 - 1, $2);
			die "Invalid alignment range $1-$2 in $partition_file\n"
				if $start < 0 || $end <= $start || $end > $alignment_length;
			push @ranges, [$start, $end];
		}
	}
	close $fh or die "Cannot close alignment partition file $partition_file: $!\n";
	die "Alignment partition file $partition_file contains no locus ranges\n"
		unless @ranges;
	@ranges = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @ranges;
	my $next = 0;
	for my $range (@ranges) {
		die "Alignment partition ranges in $partition_file do not cover the "
			."concatenated alignment exactly at position ".($next + 1)."\n"
			unless $range->[0] == $next;
		$next = $range->[1];
	}
	die "Alignment partition ranges in $partition_file end at $next, expected "
		."$alignment_length\n" unless $next == $alignment_length;
	return \@ranges;
}

sub _backbone_overlap_metrics {
	my ($seq, $backbone, $queries, $ranges, $is_aa, $alignment_length) = @_;
	my $seen_once = "\x00" x $alignment_length;
	my $supported = "\x00" x $alignment_length;
	my $state_union = "\x00" x $alignment_length;
	for my $id (@{$backbone}) {
		my $valid = _valid_mask($seq->{$id}, $is_aa);
		$supported |= ($seen_once & $valid);
		$seen_once |= $valid;
		$state_union |= _dna_state_mask($seq->{$id}) unless $is_aa;
	}
	my $nt_factor = $is_aa ? 3 : 1;
	my %metrics;
	for my $id (@{$queries}) {
		my $valid = _valid_mask($seq->{$id}, $is_aa);
		my $overlap_mask = $valid & $supported;
		my $overlap_sites = _mask_count($overlap_mask);
		my $overlap_loci = 0;
		for my $range (@{$ranges}) {
			my $locus_mask = substr(
				$overlap_mask, $range->[0], $range->[1] - $range->[0]);
			$overlap_loci++ if _mask_count($locus_mask);
		}
		my $state_divergence;
		if (!$is_aa && $overlap_sites) {
			my $matching = _dna_state_mask($seq->{$id}) & $state_union;
			$matching =~ tr/\x01\x02\x04\x08/\x01/;
			my $matching_supported = $matching & $supported;
			my $matching_sites = _mask_count($matching_supported);
			$state_divergence = ($overlap_sites - $matching_sites) / $overlap_sites;
		}
		$metrics{$id} = {
			backbone_overlap_nt => $overlap_sites * $nt_factor,
			backbone_overlap_loci => $overlap_loci,
			backbone_state_divergence => $state_divergence,
		};
	}
	return \%metrics;
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
	my $backbone_eligible = $options->{backbone_eligible} || {};
	my $backbone_ineligible_reason = $options->{backbone_ineligible_reason} || {};
	my $placement_eligible = $options->{placement_eligible} || {};
	my $placement_ineligible_reason = $options->{placement_ineligible_reason} || {};
	my $minimum_backbone_overlap_nt =
		$options->{minimum_backbone_overlap_nt} // 0;
	my $minimum_backbone_overlap_loci =
		$options->{minimum_backbone_overlap_loci} // 1;
	die "Minimum backbone overlap NT and loci must be non-negative integers\n"
		unless $minimum_backbone_overlap_nt =~ /^\d+$/
			&& $minimum_backbone_overlap_loci =~ /^\d+$/;
	my $seq = _read_fasta($full_fasta);
	my @ids = sort keys %{$seq};
	die "Strict-backbone input $full_fasta contains no sequences\n" unless @ids;
	my $alignment_length = length($seq->{$ids[0]});
	my @unequal = grep { length($seq->{$_}) != $alignment_length } @ids;
	die "Strict-backbone input $full_fasta contains unequal sequence lengths\n"
		if @unequal;
	my %informative = map { $_ => _informative_count($seq->{$_}, $is_aa) } @ids;
	my $q90 = _quantile(0.90, values %informative);
	my (%classification_reason, %requested_reason, @backbone, @placement, @excluded);
	for my $id (@ids) {
		my $sampleLocusQC = ($status->{$id} // '') eq 'placement';
		if (($informative{$id} // 0) == 0) {
			push @excluded, $id;
			$classification_reason{$id} = 'no_informative_alignment_sites';
			$requested_reason{$id} = $classification_reason{$id};
			next;
		}

		my $lowCoverage = $q90 > 0
			&& $informative{$id} < $coverage_fraction * $q90;
		$lowCoverage = 0 if length($outgroup) && $id eq $outgroup;
		my $backboneRejected = exists($backbone_eligible->{$id})
			&& !$backbone_eligible->{$id};
		$backboneRejected = 0 if length($outgroup) && $id eq $outgroup;
		if ($lowCoverage || $backboneRejected) {
			my @reason;
			push @reason, 'low_validated_coverage' if $lowCoverage;
			push @reason, ($backbone_ineligible_reason->{$id}
				// 'backbone_coverage_not_met') if $backboneRejected;
			push @reason, 'sample_locus_qc' if $sampleLocusQC;
			if (exists($placement_eligible->{$id}) && !$placement_eligible->{$id}) {
				push @reason, ($placement_ineligible_reason->{$id}
					// 'placement_coverage_not_met');
				push @excluded, $id;
			} else {
				push @placement, $id;
			}
			$classification_reason{$id} = join(',', @reason);
			$requested_reason{$id} = $classification_reason{$id};
		} else {
			push @backbone, $id;
			$classification_reason{$id} = 'retained_after_locus_qc_masking'
				if $sampleLocusQC;
		}
	}
	my $fallback = 0;
	my %backbone_overlap;
	if (@backbone < $minimum_backbone) {
		$fallback = 1;
		# Placement-eligible samples can still support an underpowered inference,
		# but samples that failed the restored coverage gate stay excluded.
		@backbone = sort (@backbone, @placement);
		@placement = ();
	} elsif (@placement) {
		my $ranges = _alignment_locus_ranges(
			$options->{partition_file}, $alignment_length);
		%backbone_overlap = %{_backbone_overlap_metrics(
			$seq, \@backbone, \@placement, $ranges, $is_aa, $alignment_length)};
		my @retained_placement;
		for my $id (@placement) {
			my $metric = $backbone_overlap{$id};
			my $exclusion_reason;
			$exclusion_reason = 'below_backbone_overlap_nt'
				if $metric->{backbone_overlap_nt} < $minimum_backbone_overlap_nt;
			$exclusion_reason = 'below_backbone_overlap_loci'
				if !defined($exclusion_reason)
					&& $metric->{backbone_overlap_loci} < $minimum_backbone_overlap_loci;
			if (defined $exclusion_reason) {
				$classification_reason{$id} .= ','
					if length($classification_reason{$id} // '');
				$classification_reason{$id} .= $exclusion_reason;
				$requested_reason{$id} = $classification_reason{$id};
				push @excluded, $id;
			} else {
				push @retained_placement, $id;
			}
		}
		@placement = @retained_placement;
	}
	_write_fasta($backbone_fasta, $seq, \@backbone);
	_write_fasta($placement_fasta, $seq, \@placement);
	return {
		backbone => \@backbone,
		placement => \@placement,
		excluded => \@excluded,
		reason => \%classification_reason,
		informative => \%informative,
		q90_informative => $q90,
		fallback => $fallback,
		requested_reason => \%requested_reason,
		backbone_overlap => \%backbone_overlap,
	};
}

sub read_epa_jplace {
	my ($jplace_file, $expected_queries) = @_;
	open my $fh, '<', $jplace_file
		or die "Cannot read EPA-ng placement file $jplace_file: $!\n";
	local $/;
	my $text = <$fh>;
	close $fh or die "Cannot close EPA-ng placement file $jplace_file: $!\n";
	my $jplace = eval { decode_json($text) };
	die "Cannot parse EPA-ng placement file $jplace_file: $@\n"
		if $@ || ref($jplace) ne 'HASH';
	die "EPA-ng placement file $jplace_file has no edge-labelled reference tree\n"
		unless defined($jplace->{tree}) && length($jplace->{tree});
	my @fields = ref($jplace->{fields}) eq 'ARRAY'
		? @{$jplace->{fields}}
		: qw(edge_num likelihood like_weight_ratio distal_length pendant_length);
	my %field_index = map { $fields[$_] => $_ } 0 .. $#fields;
	my $weight_field = exists($field_index{like_weight_ratio})
		? 'like_weight_ratio' : 'likelihood_weight_ratio';
	for my $required (qw(edge_num likelihood distal_length pendant_length), $weight_field) {
		die "EPA-ng placement file $jplace_file lacks required field '$required'\n"
			unless exists $field_index{$required};
	}
	my $root = _parse_epa_tree($jplace->{tree});
	my %edges;
	_index_epa_edges($root, undef, \%edges);
	my $tree_length = _epa_total_tree_length($root);
	my %placements = map { $_ => {status => 'not_reported'} }
		@{$expected_queries // []};
	for my $record (@{$jplace->{placements} // []}) {
		next unless ref($record) eq 'HASH'
			&& ref($record->{p}) eq 'ARRAY' && @{$record->{p}};
		my @names;
		push @names, @{$record->{n}} if ref($record->{n}) eq 'ARRAY';
		push @names, map { ref($_) eq 'ARRAY' ? $_->[0] : () }
			@{$record->{nm} // []};
		next unless @names;
		my @ranked = sort {
			($b->{likelihood_weight_ratio} // -9e99)
				<=> ($a->{likelihood_weight_ratio} // -9e99)
			|| ($b->{likelihood} // -9e99) <=> ($a->{likelihood} // -9e99)
		} map {
			{
				edge => $_->[$field_index{edge_num}],
				likelihood => $_->[$field_index{likelihood}],
				likelihood_weight_ratio => $_->[$field_index{$weight_field}],
				distal_length => $_->[$field_index{distal_length}],
				pendant_length => $_->[$field_index{pendant_length}],
			}
		} grep {
			ref($_) eq 'ARRAY' && defined($_->[$field_index{edge_num}])
		} @{$record->{p}};
		next unless @ranked;
		my $best = $ranked[0];
		my $edpl = _epa_edpl(\@ranked, \%edges, $tree_length);
		for my $name (@names) {
			next unless defined($name) && length($name);
			$placements{$name} = {
				status => 'placed', %{$best}, edpl => $edpl,
				candidate_placements => scalar(@ranked),
			};
		}
	}
	return {tree => $jplace->{tree}, placements => \%placements};
}

sub _epa_terminal_lengths {
	my ($node, $outgroup, $lengths) = @_;
	if (!@{$node->{children}}) {
		push @{$lengths}, $node->{length}
			if defined($node->{length}) && $node->{length} >= 0
				&& (!defined($outgroup) || !length($outgroup)
					|| ($node->{name} // '') ne $outgroup);
		return;
	}
	_epa_terminal_lengths($_, $outgroup, $lengths) for @{$node->{children}};
}

sub filter_epa_placement_outliers {
	my ($epa_tree, $placements, $options) = @_;
	$options ||= {};
	my $factor = $options->{pendant_outlier_factor} // 5;
	my $minimum_threshold = $options->{pendant_minimum_threshold} // 0.02;
	die "EPA pendant outlier factor and minimum threshold must be non-negative\n"
		if $factor < 0 || $minimum_threshold < 0;
	return {
		enabled => 0, excluded => [], retained => [],
		backbone_q95 => undef, threshold => undef,
	} if $factor == 0;
	die "EPA placement outlier filtering requires a placement hash\n"
		unless ref($placements) eq 'HASH';
	my $root = _parse_epa_tree($epa_tree);
	my @terminal_lengths;
	_epa_terminal_lengths($root, $options->{outgroup}, \@terminal_lengths);
	die "EPA placement outlier filtering found no usable backbone terminal branches\n"
		unless @terminal_lengths;
	my $backbone_q95 = _quantile(0.95, @terminal_lengths);
	my $threshold = $factor * $backbone_q95;
	$threshold = $minimum_threshold if $threshold < $minimum_threshold;
	my (@excluded, @retained);
	for my $sample (sort keys %{$placements}) {
		my $placement = $placements->{$sample};
		next unless ($placement->{status} // '') eq 'placed';
		$placement->{pendant_outlier_limit} = $threshold;
		if (defined($placement->{pendant_length})
				&& $placement->{pendant_length} > $threshold) {
			$placement->{status} = 'excluded_outlier';
			$placement->{placement_filter_reason} = 'pendant_length_outlier';
			push @excluded, $sample;
		} else {
			push @retained, $sample;
		}
	}
	return {
		enabled => 1, excluded => \@excluded, retained => \@retained,
		backbone_q95 => $backbone_q95, threshold => $threshold,
	};
}

sub _skip_newick_space {
	my ($text, $position_ref) = @_;
	$$position_ref++ while $$position_ref < length($text)
		&& substr($text, $$position_ref, 1) =~ /\s/;
}

sub _parse_newick_label {
	my ($text, $position_ref) = @_;
	_skip_newick_space($text, $position_ref);
	return '' if $$position_ref >= length($text)
		|| substr($text, $$position_ref, 1) =~ /[:;,(){}]/;
	if (substr($text, $$position_ref, 1) eq "'") {
		$$position_ref++;
		my $label = '';
		while ($$position_ref < length($text)) {
			my $character = substr($text, $$position_ref, 1);
			$$position_ref++;
			if ($character eq "'") {
				if ($$position_ref < length($text)
					&& substr($text, $$position_ref, 1) eq "'") {
					$label .= "'";
					$$position_ref++;
					next;
				}
				return $label;
			}
			$label .= $character;
		}
		die "Unterminated quoted Newick label\n";
	}
	my $start = $$position_ref;
	$$position_ref++ while $$position_ref < length($text)
		&& substr($text, $$position_ref, 1) !~ /[\s:;,(){}]/;
	return substr($text, $start, $$position_ref - $start);
}

sub _parse_newick_branch {
	my ($text, $position_ref) = @_;
	_skip_newick_space($text, $position_ref);
	my $node;
	if (substr($text, $$position_ref, 1) eq '(') {
		$$position_ref++;
		my @children;
		while (1) {
			push @children, _parse_newick_branch($text, $position_ref);
			_skip_newick_space($text, $position_ref);
			my $separator = substr($text, $$position_ref, 1);
			if ($separator eq ',') {
				$$position_ref++;
				next;
			}
			die "Malformed Newick tree: expected ')'\n" unless $separator eq ')';
			$$position_ref++;
			$node = {children => \@children, name => _parse_newick_label($text, $position_ref)};
			last;
		}
	} else {
		my $name = _parse_newick_label($text, $position_ref);
		die "Malformed Newick tree: missing leaf label\n" unless length($name);
		$node = {children => [], name => $name};
	}
	_skip_newick_space($text, $position_ref);
	if (substr($text, $$position_ref, 1) eq ':') {
		$$position_ref++;
		_skip_newick_space($text, $position_ref);
		my $start = $$position_ref;
		$$position_ref++ while $$position_ref < length($text)
			&& substr($text, $$position_ref, 1) !~ /[\s,;(){}]/;
		my $length = substr($text, $start, $$position_ref - $start);
		die "Malformed Newick branch length '$length'\n"
			unless $length =~ /^[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?$/;
		$node->{length} = 0 + $length;
	}
	_skip_newick_space($text, $position_ref);
	if (substr($text, $$position_ref, 1) eq '{') {
		$$position_ref++;
		my $start = $$position_ref;
		$$position_ref++ while $$position_ref < length($text)
			&& substr($text, $$position_ref, 1) =~ /\d/;
		my $edge = substr($text, $start, $$position_ref - $start);
		die "Malformed EPA-ng edge label\n" unless length($edge)
			&& substr($text, $$position_ref, 1) eq '}';
		$$position_ref++;
		$node->{edge} = 0 + $edge;
	}
	return $node;
}

sub _parse_epa_tree {
	my ($text) = @_;
	my $position = 0;
	my $root = _parse_newick_branch($text, \$position);
	_skip_newick_space($text, \$position);
	die "Malformed EPA-ng reference tree: expected terminating ';'\n"
		unless substr($text, $position, 1) eq ';';
	$position++;
	_skip_newick_space($text, \$position);
	die "Malformed EPA-ng reference tree: trailing characters\n" if $position < length($text);
	return $root;
}

sub _index_epa_edges {
	my ($node, $parent, $edges) = @_;
	$node->{parent} = $parent if defined $parent;
	die "Duplicate EPA-ng edge label $node->{edge}\n"
		if defined($node->{edge}) && exists($edges->{$node->{edge}});
	$edges->{$node->{edge}} = $node if defined $node->{edge};
	_index_epa_edges($_, $node, $edges) for @{$node->{children}};
}
sub _epa_total_tree_length {
	my ($node) = @_;
	my $total = defined($node->{parent}) ? ($node->{length} // 0) : 0;
	$total += _epa_total_tree_length($_) for @{$node->{children}};
	return $total;
}

sub _epa_node_distance {
	my ($left, $right) = @_;
	my (%left_distance, $distance);
	$distance = 0;
	for (my $node = $left; defined($node); $node = $node->{parent}) {
		$left_distance{$node} = $distance;
		$distance += $node->{length} // 0 if defined $node->{parent};
	}
	$distance = 0;
	for (my $node = $right; defined($node); $node = $node->{parent}) {
		return $distance + $left_distance{$node}
			if exists $left_distance{$node};
		$distance += $node->{length} // 0 if defined $node->{parent};
	}
	die "EPA-ng placement tree contains disconnected nodes\n";
}

sub _epa_placement_distance {
	my ($left, $right, $edges) = @_;
	return undef unless defined($left->{edge}) && defined($right->{edge});
	my $left_child = $edges->{$left->{edge}};
	my $right_child = $edges->{$right->{edge}};
	return undef unless $left_child && $right_child;
	my $left_parent = $left_child->{parent};
	my $right_parent = $right_child->{parent};
	return undef unless $left_parent && $right_parent;
	my $left_length = $left_child->{length};
	my $right_length = $right_child->{length};
	return undef unless defined($left_length) && defined($right_length);
	my $left_distal = $left->{distal_length} // 0;
	my $right_distal = $right->{distal_length} // 0;
	$left_distal = 0 if $left_distal < 0;
	$right_distal = 0 if $right_distal < 0;
	$left_distal = $left_length if $left_distal > $left_length;
	$right_distal = $right_length if $right_distal > $right_length;
	return abs($left_distal - $right_distal)
		if $left->{edge} == $right->{edge};
	my @left_endpoint = (
		[$left_child, $left_distal],
		[$left_parent, $left_length - $left_distal],
	);
	my @right_endpoint = (
		[$right_child, $right_distal],
		[$right_parent, $right_length - $right_distal],
	);
	my $minimum;
	for my $left_end (@left_endpoint) {
		for my $right_end (@right_endpoint) {
			my $distance = $left_end->[1] + $right_end->[1]
				+ _epa_node_distance($left_end->[0], $right_end->[0]);
			$minimum = $distance
				if !defined($minimum) || $distance < $minimum;
		}
	}
	return $minimum;
}

sub _epa_edpl {
	my ($placements, $edges, $tree_length) = @_;
	return 0 unless $tree_length > 0 && @{$placements} > 1;
	my @usable = grep {
		defined($_->{likelihood_weight_ratio})
			&& $_->{likelihood_weight_ratio} > 0
			&& defined($_->{edge}) && exists($edges->{$_->{edge}})
	} @{$placements};
	return 0 unless @usable > 1;
	my $weight_sum = 0;
	$weight_sum += $_->{likelihood_weight_ratio} for @usable;
	return 0 unless $weight_sum > 0;
	my $edpl = 0;
	for my $left_index (0 .. $#usable - 1) {
		for my $right_index ($left_index + 1 .. $#usable) {
			my $distance = _epa_placement_distance(
				$usable[$left_index], $usable[$right_index], $edges);
			next unless defined $distance;
			my $left_weight =
				$usable[$left_index]{likelihood_weight_ratio} / $weight_sum;
			my $right_weight =
				$usable[$right_index]{likelihood_weight_ratio} / $weight_sum;
			$edpl += 2 * $left_weight * $right_weight
				* $distance / $tree_length;
		}
	}
	return $edpl;
}


sub _replace_child {
	my ($parent, $old, $new) = @_;
	for my $index (0 .. $#{$parent->{children}}) {
		next unless $parent->{children}[$index] == $old;
		$parent->{children}[$index] = $new;
		return;
	}
	die "Could not replace EPA-ng reference-tree branch\n";
}

sub _newick_label {
	my ($label) = @_;
	return $label if $label =~ /^[A-Za-z0-9_.|+-]+$/;
	$label =~ s/'/''/g;
	return "'$label'";
}

sub _render_newick {
	my ($node) = @_;
	my $text = @{$node->{children}}
		? '('.join(',', map { _render_newick($_) } @{$node->{children}}).')'
		: '';
	$text .= _newick_label($node->{name}) if defined($node->{name}) && length($node->{name});
	$text .= ':'.sprintf('%.12g', $node->{length}) if defined $node->{length};
	return $text;
}

sub write_epa_placed_tree {
	my ($epa_tree, $output_tree, $placements) = @_;
	my $root = _parse_epa_tree($epa_tree);
	my %edge;
	_index_epa_edges($root, undef, \%edge);
	my %on_edge;
	for my $query (sort keys %{$placements}) {
		my $placement = $placements->{$query};
		next unless ($placement->{status} // '') eq 'placed';
		for my $field (qw(edge distal_length pendant_length)) {
			die "EPA-ng placement for $query has no $field\n"
				unless defined($placement->{$field}) && $placement->{$field} =~ /^[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?$/;
		}
		push @{$on_edge{$placement->{edge}}}, [$query, $placement];
	}
	for my $edge_number (sort { $a <=> $b } keys %on_edge) {
		my $child = $edge{$edge_number}
			or die "EPA-ng placement refers to absent edge $edge_number\n";
		my $parent = $child->{parent}
			or die "EPA-ng placement refers to the root edge $edge_number\n";
		my $branch_length = $child->{length};
		die "EPA-ng reference edge $edge_number has no branch length\n"
			unless defined($branch_length) && $branch_length >= 0;
		my @points = sort {
			$a->[1]{distal_length} <=> $b->[1]{distal_length}
			|| $a->[0] cmp $b->[0]
		} @{$on_edge{$edge_number}};
		my $previous_distance = -1;
		my $attachment = $child;
		while (@points) {
			my $distance = $points[0][1]{distal_length};
			die "EPA-ng placement distance outside reference edge $edge_number\n"
				if $distance < -1e-8 || $distance > $branch_length + 1e-8;
			$distance = 0 if $distance < 0;
			$distance = $branch_length if $distance > $branch_length;
			my @at_point;
			while (@points && abs($points[0][1]{distal_length} - $distance) < 1e-10) {
				push @at_point, shift @points;
			}
			my $node;
			if ($previous_distance >= 0 && abs($distance - $previous_distance) < 1e-10) {
				$node = $attachment;
			} else {
				$node = {children => [$attachment]};
				$attachment->{parent} = $node;
				$attachment->{length} = $distance - ($previous_distance < 0 ? 0 : $previous_distance);
				$attachment = $node;
			}
			for my $entry (@at_point) {
				my ($query, $placement) = @{$entry};
				my $tip = {children => [], name => $query,
					length => 0 + $placement->{pendant_length}, parent => $node};
				push @{$node->{children}}, $tip;
			}
			$previous_distance = $distance;
		}
		$attachment->{length} = $branch_length - $previous_distance;
		_replace_child($parent, $child, $attachment);
		$attachment->{parent} = $parent;
	}
	my $temporary_output = "$output_tree.tmp.$$";
	unlink $temporary_output
		or die "Cannot remove stale placed-tree temporary $temporary_output: $!\n"
		if -e $temporary_output;
	open my $out, '>', $temporary_output
		or die "Cannot write placed tree $temporary_output: $!\n";
	print {$out} _render_newick($root), ";\n";
	close $out or die "Cannot close placed tree $temporary_output: $!\n";
	rename $temporary_output, $output_tree
		or die "Cannot publish placed tree $temporary_output as $output_tree: $!\n";
}

1;
