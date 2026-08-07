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
	my $seq = _read_fasta($full_fasta);
	my @ids = sort keys %{$seq};
	die "Strict-backbone input $full_fasta contains no sequences\n" unless @ids;
	my %informative = map { $_ => _informative_count($seq->{$_}, $is_aa) } @ids;
	my $q90 = _quantile(0.90, values %informative);
	my (%classification_reason, %requested_reason, @backbone, @placement, @excluded);
	for my $id (@ids) {
		my $sampleLocusQC = ($status->{$id} // '') eq 'placement';
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
	if (@backbone < $minimum_backbone) {
		$fallback = 1;
		# Placement-eligible samples can still support an underpowered inference,
		# but samples that failed the restored coverage gate stay excluded.
		@backbone = sort (@backbone, @placement);
		@placement = ();
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
	die "Cannot parse EPA-ng placement file $jplace_file: $@\n" if $@ || ref($jplace) ne 'HASH';
	die "EPA-ng placement file $jplace_file has no edge-labelled reference tree\n"
		unless defined($jplace->{tree}) && length($jplace->{tree});
	my %placements = map { $_ => {status => 'not_reported'} } @{$expected_queries // []};
	for my $record (@{$jplace->{placements} // []}) {
		next unless ref($record) eq 'HASH' && ref($record->{p}) eq 'ARRAY' && @{$record->{p}};
		my @names;
		push @names, @{$record->{n}} if ref($record->{n}) eq 'ARRAY';
		push @names, map { ref($_) eq 'ARRAY' ? $_->[0] : () } @{$record->{nm} // []};
		next unless @names;
		my @ranked = sort {
			($b->[2] // -9e99) <=> ($a->[2] // -9e99)
			|| ($b->[1] // -9e99) <=> ($a->[1] // -9e99)
		} grep { ref($_) eq 'ARRAY' && defined($_->[0]) } @{$record->{p}};
		next unless @ranked;
		my $best = $ranked[0];
		for my $name (@names) {
			next unless defined($name) && length($name);
			$placements{$name} = {
				status => 'placed', edge => $best->[0], likelihood => $best->[1],
				likelihood_weight_ratio => $best->[2], distal_length => $best->[3],
				pendant_length => $best->[4],
			};
		}
	}
	return {tree => $jplace->{tree}, placements => \%placements};
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
