package Mods::StrainPlacement;

use strict;
use warnings;
use Exporter qw(import);
use JSON::PP qw(decode_json);
use Mods::GenoMetaAss qw(fileGZe gzipopen);

our @EXPORT_OK = qw(
	read_sample_qc
	canonical_sample_qc_status
	split_strict_backbone
	read_epa_jplace
	filter_epa_placement_outliers
	map_epa_placements_to_backbone
	write_epa_placed_tree
);

#The sample QC verdict records what extraction found about a sample, not what a
#later stage does with it: whether it is deleted, kept out of the strict
#backbone, or retained depends on -excludeFlaggedSamples and -placeOnBackbone.
#The pre-1.08 spellings named one such disposition ('placement') and so
#described an action the default configuration does not take; they are still
#accepted so existing sampleQC tables keep loading.
our %SAMPLE_QC_STATUS_ALIAS = (
	single_strain => 'single_strain',
	mixed_strain  => 'mixed_strain',
	backbone      => 'single_strain',
	placement     => 'mixed_strain',
);

sub canonical_sample_qc_status {
	my ($raw) = @_;
	return '' unless defined($raw) && length($raw);
	return $SAMPLE_QC_STATUS_ALIAS{$raw} // '';
}

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
		my ($sample, $raw_status) = @field[1, 2];
		my $sample_status = canonical_sample_qc_status($raw_status);
		die "Invalid sample QC status '$raw_status' for $sample\n"
			unless length $sample_status;
		# A sample seen under several MGS keeps its worst verdict.
		$status{$sample} = $sample_status
			if !exists($status{$sample}) || $sample_status eq 'mixed_strain';
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
	my @coverage_baseline = grep { !length($outgroup) || $_ ne $outgroup } @ids;
	@coverage_baseline = @ids unless @coverage_baseline;
	my $q90 = _quantile(0.90,
		map { $informative{$_} } @coverage_baseline);
	my (%classification_reason, %requested_reason, @backbone, @placement, @excluded);
	for my $id (@ids) {
		# read_sample_qc() returns canonical verdicts, but this exported helper is
		# also used directly by callers that may still carry the documented legacy
		# backbone/placement spellings. Apply the same compatibility mapping here so
		# the route does not depend on which public entry point supplied the hash.
		my $sampleLocusQC =
			canonical_sample_qc_status($status->{$id}) eq 'mixed_strain';
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
		backbone_q95 => undef, backbone_terminal_count => 0,
		factor => $factor, minimum_threshold => $minimum_threshold,
		scaled_threshold => undef, threshold => undef,
		threshold_source => 'disabled', placed_query_count => 0,
		query_pendant_count => 0, query_pendant_missing_count => 0,
	} if $factor == 0;
	die "EPA placement outlier filtering requires a placement hash\n"
		unless ref($placements) eq 'HASH';
	my $root = _parse_epa_tree($epa_tree);
	my @terminal_lengths;
	_epa_terminal_lengths($root, $options->{outgroup}, \@terminal_lengths);
	die "EPA placement outlier filtering found no usable backbone terminal branches\n"
		unless @terminal_lengths;
	my $backbone_q95 = _quantile(0.95, @terminal_lengths);
	my $scaled_threshold = $factor * $backbone_q95;
	my $threshold = $scaled_threshold;
	my $threshold_source = 'scaled_backbone_q95';
	if ($threshold < $minimum_threshold) {
		$threshold = $minimum_threshold;
		$threshold_source = 'minimum_floor';
	}
	my (@excluded, @retained, @query_pendant_lengths);
	my ($placed_query_count, $query_pendant_missing_count) = (0, 0);
	for my $sample (sort keys %{$placements}) {
		my $placement = $placements->{$sample};
		next unless ($placement->{status} // '') eq 'placed';
		$placed_query_count++;
		$placement->{pendant_outlier_limit} = $threshold;
		if (defined($placement->{pendant_length})) {
			push @query_pendant_lengths, $placement->{pendant_length};
		} else {
			$query_pendant_missing_count++;
		}
		if (defined($placement->{pendant_length})
				&& $placement->{pendant_length} > $threshold) {
			$placement->{status} = 'excluded_outlier';
			$placement->{placement_filter_reason} = 'pendant_length_outlier';
			push @excluded, $sample;
		} else {
			push @retained, $sample;
		}
	}
	my @sorted_pendant_lengths = sort { $a <=> $b } @query_pendant_lengths;
	return {
		enabled => 1, excluded => \@excluded, retained => \@retained,
		backbone_q95 => $backbone_q95,
		backbone_terminal_count => scalar(@terminal_lengths),
		factor => $factor, minimum_threshold => $minimum_threshold,
		scaled_threshold => $scaled_threshold, threshold => $threshold,
		threshold_source => $threshold_source,
		placed_query_count => $placed_query_count,
		query_pendant_count => scalar(@sorted_pendant_lengths),
		query_pendant_missing_count => $query_pendant_missing_count,
		query_pendant_min => @sorted_pendant_lengths ? $sorted_pendant_lengths[0] : undef,
		query_pendant_median => @sorted_pendant_lengths ? _quantile(0.5, @sorted_pendant_lengths) : undef,
		query_pendant_q95 => @sorted_pendant_lengths ? _quantile(0.95, @sorted_pendant_lengths) : undef,
		query_pendant_max => @sorted_pendant_lengths ? $sorted_pendant_lengths[-1] : undef,
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

sub _reference_name_key {
	return join('', map { length($_).':'.$_.';' } @_);
}

sub _collect_reference_nodes {
	my ($node, $parent, $nodes, $terminals, $label) = @_;
	$node->{parent} = $parent if defined $parent;
	my @descendants;
	if (@{$node->{children}}) {
		push @descendants,
			@{_collect_reference_nodes($_, $node, $nodes, $terminals, $label)}
			for @{$node->{children}};
		@descendants = sort @descendants;
	} else {
		my $name = $node->{name} // '';
		die "$label contains an unnamed terminal\n" unless length($name);
		die "$label contains duplicate terminal '$name'\n"
			if $terminals->{$name}++;
		@descendants = ($name);
	}
	push @{$nodes}, {
		node => $node,
		rooted_key => _reference_name_key(@descendants),
		rooted_taxa => \@descendants,
	} if defined $parent;
	return \@descendants;
}

sub _canonical_reference_split {
	my ($descendants, $all_taxa) = @_;
	my %distal = map { $_ => 1 } @{$descendants};
	my @proximal = grep { !$distal{$_} } @{$all_taxa};
	my @distal = @{$descendants};
	my $distal_key = _reference_name_key(@distal);
	my $proximal_key = _reference_name_key(@proximal);
	return @proximal < @distal
		|| (@proximal == @distal && $proximal_key lt $distal_key)
		? ($proximal_key, \@proximal)
		: ($distal_key, \@distal);
}

sub _index_reference_splits {
	my ($root, $label, $require_edge_labels) = @_;
	my (@nodes, %terminals);
	my $all_taxa = _collect_reference_nodes(
		$root, undef, \@nodes, \%terminals, $label);
	my (%by_split, %by_edge);
	for my $entry (@nodes) {
		my ($split_key, $split_taxa) = _canonical_reference_split(
			$entry->{rooted_taxa}, $all_taxa);
		$entry->{split_key} = $split_key;
		$entry->{split_taxa} = $split_taxa;
		push @{$by_split{$split_key}}, $entry;
		next unless $require_edge_labels;
		my $edge = $entry->{node}{edge};
		die "$label branch has no EPA edge number\n" unless defined $edge;
		die "$label contains duplicate EPA edge number $edge\n"
			if exists $by_edge{$edge};
		$entry->{edge} = $edge;
		$by_edge{$edge} = $entry;
	}
	return {
		all_taxa => $all_taxa,
		terminals => \%terminals,
		by_split => \%by_split,
		by_edge => \%by_edge,
	};
}

sub _match_epa_edges_to_backbone {
	my ($epa_index, $backbone_index) = @_;
	die "EPA-ng jplace and authoritative backbone terminal sets differ\n"
		unless _reference_name_key(@{$epa_index->{all_taxa}}) eq
			_reference_name_key(@{$backbone_index->{all_taxa}});
	my %edge_map;
	for my $split_key (keys %{$epa_index->{by_split}}) {
		my $epa_entries = $epa_index->{by_split}{$split_key};
		my $backbone_entries = $backbone_index->{by_split}{$split_key}
			or die "EPA-ng jplace contains a branch split absent from the "
				."authoritative backbone\n";
		die "EPA-ng jplace and authoritative backbone represent split "
			."$split_key a different number of times\n"
			unless @{$epa_entries} == @{$backbone_entries};
		for my $epa_entry (@{$epa_entries}) {
			my @matches = @{$backbone_entries} == 1
				? @{$backbone_entries}
				: grep {
					$_->{rooted_key} eq $epa_entry->{rooted_key}
				} @{$backbone_entries};
			die "Cannot disambiguate binary-root halves for EPA edge "
				."$epa_entry->{edge} on the authoritative backbone\n"
				unless @matches == 1;
			$edge_map{$epa_entry->{edge}} = {
				epa => $epa_entry,
				backbone => $matches[0],
			};
		}
	}
	for my $split_key (keys %{$backbone_index->{by_split}}) {
		die "Authoritative backbone contains a branch split absent from the "
			."EPA-ng jplace tree\n"
			unless exists $epa_index->{by_split}{$split_key};
	}
	return \%edge_map;
}

sub map_epa_placements_to_backbone {
	my ($epa_tree, $backbone_tree, $placements) = @_;
	die "EPA backbone graft mapping requires a placement hash\n"
		unless ref($placements) eq 'HASH';
	my $epa_root = _parse_epa_tree($epa_tree);
	my $backbone_root = _parse_epa_tree($backbone_tree);
	my $epa_index = _index_reference_splits(
		$epa_root, 'EPA-ng jplace tree', 1);
	my $backbone_index = _index_reference_splits(
		$backbone_root, 'authoritative backbone tree', 0);
	my $edge_map = _match_epa_edges_to_backbone(
		$epa_index, $backbone_index);

	my %placement_count;
	for my $sample (keys %{$placements}) {
		my $edge = $placements->{$sample}{edge};
		$placement_count{$edge}++ if defined $edge;
	}
	my (@rows, %row_by_edge);
	my ($different_lengths, $jplace_lengths_missing, $zero_backbone_edges,
		$max_difference) = (0, 0, 0, 0);
	for my $edge (sort { $a <=> $b } keys %{$edge_map}) {
		my $mapping = $edge_map->{$edge};
		my $epa_node = $mapping->{epa}{node};
		my $backbone_node = $mapping->{backbone}{node};
		die "Authoritative backbone branch mapped from EPA edge $edge "
			."has no non-negative length\n"
			unless defined($backbone_node->{length})
				&& $backbone_node->{length} >= 0;
		my $jplace_length = defined($epa_node->{length})
			? 0 + $epa_node->{length} : undef;
		my $backbone_length = 0 + $backbone_node->{length};
		my $difference = defined($jplace_length)
			? $jplace_length - $backbone_length : undef;
		my $is_different = defined($difference)
			&& $difference != 0 ? 1 : 0;
		$different_lengths += $is_different;
		$jplace_lengths_missing++ unless defined $jplace_length;
		$zero_backbone_edges++ if $backbone_length == 0;
		$max_difference = abs($difference)
			if defined($difference) && abs($difference) > $max_difference;
		my $split_taxa = $mapping->{backbone}{split_taxa};
		my $row = {
			edge => $edge,
			edge_type => @{$split_taxa} == 1 ? 'terminal' : 'internal',
			terminal => @{$split_taxa} == 1 ? $split_taxa->[0] : '',
			split_size => scalar(@{$split_taxa}),
			jplace_length => $jplace_length,
			backbone_length => $backbone_length,
			difference => $difference,
			changed => $is_different,
			placement_count => $placement_count{$edge} // 0,
			clamped_placement_count => 0,
		};
		push @rows, $row;
		$row_by_edge{$edge} = $row;
	}

	my ($mapped_placements, $clamped_placements) = (0, 0);
	for my $sample (sort keys %{$placements}) {
		my $placement = $placements->{$sample};
		next unless defined $placement->{edge};
		my $mapping = $edge_map->{$placement->{edge}}
			or die "EPA-ng placement for $sample refers to the root or an "
				."unmapped edge $placement->{edge}\n";
		die "EPA-ng placement for $sample has no numeric distal_length\n"
			unless defined($placement->{distal_length})
				&& $placement->{distal_length}
					=~ /^[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?$/;
		my $backbone_length = 0 + $mapping->{backbone}{node}{length};
		my $applied_distance = 0 + $placement->{distal_length};
		my $clamped = 0;
		if ($applied_distance < 0) {
			$applied_distance = 0;
			$clamped = 1;
		} elsif ($applied_distance > $backbone_length) {
			$applied_distance = $backbone_length;
			$clamped = 1;
		}
		$placement->{backbone_split_key} =
			$mapping->{backbone}{split_key};
		$placement->{backbone_rooted_key} =
			$mapping->{backbone}{rooted_key};
		$placement->{backbone_distal_length} = $applied_distance;
		$mapped_placements++;
		if ($clamped) {
			$clamped_placements++;
			$row_by_edge{$placement->{edge}}{clamped_placement_count}++;
		}
	}
	return {
		rows => \@rows,
		compared_edge_count => scalar(@rows),
		different_length_count => $different_lengths,
		jplace_length_missing_count => $jplace_lengths_missing,
		zero_backbone_edge_count => $zero_backbone_edges,
		max_absolute_difference => $max_difference,
		mapped_placement_count => $mapped_placements,
		clamped_placement_count => $clamped_placements,
	};
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

sub _backbone_entry_for_placement {
	my ($backbone_index, $placement, $query) = @_;
	for my $field (qw(backbone_split_key backbone_rooted_key)) {
		die "EPA-ng placement for $query has no mapped $field\n"
			unless defined($placement->{$field})
				&& length($placement->{$field});
	}
	my $entries =
		$backbone_index->{by_split}{$placement->{backbone_split_key}}
		or die "Mapped backbone split for EPA-ng placement $query is absent\n";
	my @matches = @{$entries} == 1
		? @{$entries}
		: grep {
			$_->{rooted_key} eq $placement->{backbone_rooted_key}
		} @{$entries};
	die "Mapped backbone branch for EPA-ng placement $query is ambiguous\n"
		unless @matches == 1;
	return $matches[0];
}

sub write_epa_placed_tree {
	my ($backbone_tree, $output_tree, $placements) = @_;
	my $root = _parse_epa_tree($backbone_tree);
	my $backbone_index = _index_reference_splits(
		$root, 'authoritative backbone tree', 0);
	my %on_branch;
	for my $query (sort keys %{$placements}) {
		my $placement = $placements->{$query};
		next unless ($placement->{status} // '') eq 'placed';
		for my $field (qw(backbone_distal_length pendant_length)) {
			die "EPA-ng placement for $query has no $field\n"
				unless defined($placement->{$field})
					&& $placement->{$field}
						=~ /^[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?$/;
		}
		my $entry = _backbone_entry_for_placement(
			$backbone_index, $placement, $query);
		my $branch_key = _reference_name_key(
			$entry->{split_key}, $entry->{rooted_key});
		$on_branch{$branch_key}{entry} = $entry;
		push @{$on_branch{$branch_key}{placements}},
			[$query, $placement, 0 + $placement->{backbone_distal_length}];
	}
	for my $branch_key (sort keys %on_branch) {
		my $group = $on_branch{$branch_key};
		my $child = $group->{entry}{node};
		my $parent = $child->{parent}
			or die "EPA-ng placement maps to the backbone root rather than a branch\n";
		my $branch_length = $child->{length};
		die "Authoritative backbone branch has no non-negative length\n"
			unless defined($branch_length) && $branch_length >= 0;
		my @points = sort {
			$a->[2] <=> $b->[2] || $a->[0] cmp $b->[0]
		} @{$group->{placements}};
		my $previous_distance = -1;
		my $attachment = $child;
		while (@points) {
			my $distance = $points[0][2];
			die "Mapped EPA-ng placement distance outside authoritative "
				."backbone branch\n"
				if $distance < -1e-8
					|| $distance > $branch_length + 1e-8;
			$distance = 0 if $distance < 0;
			$distance = $branch_length if $distance > $branch_length;
			my @at_point;
			while (@points && abs($points[0][2] - $distance) < 1e-10) {
				push @at_point, shift @points;
			}
			my $node;
			if ($previous_distance >= 0
					&& abs($distance - $previous_distance) < 1e-10) {
				$node = $attachment;
			} else {
				$node = {children => [$attachment]};
				$attachment->{parent} = $node;
				$attachment->{length} = $distance
					- ($previous_distance < 0 ? 0 : $previous_distance);
				$attachment = $node;
			}
			for my $entry (@at_point) {
				my ($query, $placement) = @{$entry};
				my $tip = {
					children => [],
					name => $query,
					length => 0 + $placement->{pendant_length},
					parent => $node,
				};
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
