package Mods::WorkflowState;

use warnings;
use strict;

use Exporter qw(import);
use JSON::PP;

our @EXPORT_OK = qw(inspect_workflow_state encode_state_report);

sub _join_path {
	my @parts = @_;
	my $path = shift @parts;
	$path =~ s{/+$}{};
	for my $part (@parts) {
		$part =~ s{^/+}{};
		$part =~ s{/+$}{};
		$path .= "/$part";
	}
	return $path;
}

sub _base_output_dir {
	my ($output_dir, $sample_key) = @_;
	my $base = $output_dir;
	$base =~ s{/+$}{};
	return "$base/" if ($base =~ s{\Q$sample_key\E$}{});
	$base =~ s{/[^/]+$}{};
	return "$base/";
}

sub _normalise_path {
	my ($path) = @_;
	$path = "" unless (defined $path);
	$path =~ s{\\}{/}g;
	$path =~ s{/+}{/}g;
	$path =~ s{/$}{};
	return $path;
}

sub _file_state {
	my ($requested_path, $allow_gzip) = @_;
	my @candidates = ($requested_path);
	if ($allow_gzip) {
		if ($requested_path =~ m{\.gz$}) {
			(my $plain = $requested_path) =~ s{\.gz$}{};
			push @candidates, $plain;
		} else {
			push @candidates, "$requested_path.gz";
		}
	}

	my %seen;
	for my $path (grep { !$seen{$_}++ } @candidates) {
		next unless (-e $path);
		my $size = -s $path;
		$size = 0 unless (defined $size);
		return {
			requested_path => $requested_path,
			path => $path,
			exists => 1,
			nonempty => $size > 0 ? 1 : 0,
			size_bytes => 0 + $size,
		};
	}

	return {
		requested_path => $requested_path,
		path => undef,
		exists => 0,
		nonempty => 0,
		size_bytes => 0,
	};
}

sub _stage_state {
	my (%args) = @_;
	my @artifacts = map { _file_state($_, 0) } @{$args{artifacts} || []};
	push @artifacts, map { _file_state($_, 1) } @{$args{gzip_artifacts} || []};
	my @markers = map { _file_state($_, 0) } @{$args{markers} || []};

	my $artifacts_complete = @artifacts ? 1 : 0;
	$artifacts_complete = 0 if (grep { !$_->{nonempty} } @artifacts);
	my $markers_complete = 1;
	$markers_complete = 0 if (grep { !$_->{exists} } @markers);
	my $complete = $artifacts_complete && $markers_complete;
	my $anything_exists = grep { $_->{exists} } (@artifacts, @markers);

	my @issues;
	if ((grep { $_->{exists} } @markers) && grep { !$_->{nonempty} } @artifacts) {
		push @issues, 'MARKER_WITHOUT_VALID_OUTPUT';
	}
	if ((grep { $_->{nonempty} } @artifacts) && grep { !$_->{exists} } @markers) {
		push @issues, 'OUTPUT_WITHOUT_MARKER';
	}
	if (grep { $_->{exists} && !$_->{nonempty} } @artifacts) {
		push @issues, 'EMPTY_OUTPUT';
	}

	return {
		status => $complete ? 'COMPLETE' : ($anything_exists ? 'PARTIAL' : 'MISSING'),
		complete => $complete ? 1 : 0,
		artifacts => \@artifacts,
		markers => \@markers,
		issues => \@issues,
	};
}

sub _alternative_output_state {
	my (@paths) = @_;
	my @artifacts = map { _file_state($_, 1) } @paths;
	my $complete = grep { $_->{nonempty} } @artifacts;
	my $anything_exists = grep { $_->{exists} } @artifacts;
	return {
		status => $complete ? 'COMPLETE' : ($anything_exists ? 'PARTIAL' : 'MISSING'),
		complete => $complete ? 1 : 0,
		artifacts => \@artifacts,
		markers => [],
		issues => $anything_exists && !$complete ? ['EMPTY_OUTPUT'] : [],
	};
}

sub _not_applicable_state {
	return {
		status => 'NOT_APPLICABLE',
		complete => 1,
		artifacts => [],
		markers => [],
		issues => [],
	};
}

sub _unknown_state {
	my ($reason) = @_;
	return {
		status => 'UNKNOWN',
		complete => 0,
		artifacts => [],
		markers => [],
		issues => [$reason],
	};
}

sub _support_technology {
	my ($support_reads) = @_;
	return '' unless (defined $support_reads && $support_reads ne '');
	return uc($1) if ($support_reads =~ m{(?:^|[,;])(PB|ONT|mate):}i);
	return 'OTHER';
}

sub _read_membership_file {
	my ($path) = @_;
	return (0, []) unless (-e $path);
	open my $fh, '<', $path or return (1, []);
	my @members;
	while (my $line = <$fh>) {
		chomp $line;
		next if ($line =~ m{^\s*$});
		push @members, _normalise_path($line);
	}
	close $fh;
	return (1, \@members);
}

sub _same_members {
	my ($left, $right) = @_;
	return 0 unless (@{$left} == @{$right});
	my @left_sorted = sort @{$left};
	my @right_sorted = sort @{$right};
	for (my $i = 0; $i < @left_sorted; $i++) {
		return 0 if ($left_sorted[$i] ne $right_sorted[$i]);
	}
	return 1;
}

sub inspect_workflow_state {
	my (%args) = @_;
	my $map = $args{map} || die 'inspect_workflow_state requires map';
	my $groups = $args{groups} || die 'inspect_workflow_state requires groups';
	my $options = $args{options} || {};
	my $assembly_mode = exists($options->{assembly_mode}) ? 0 + $options->{assembly_mode} : 1;
	my $assembly_requested = $assembly_mode ? 1 : 0;
	my $map_to_assembly = exists($options->{map_to_assembly})
		? ($options->{map_to_assembly} ? 1 : 0) : 1;
	my $map_support_to_assembly = exists($options->{map_support_to_assembly})
		? ($options->{map_support_to_assembly} ? 1 : 0) : 1;
	my $run_tmp_dir = _normalise_path($options->{run_tmp_dir} || '');
	my @sample_order = @{$map->{opt}{smpl_order} || []};

	my %group_members;
	for my $sample_key (@sample_order) {
		my $group_id = $map->{$sample_key}{AssGroup};
		$group_id = $sample_key if (!defined $group_id || $group_id eq '-1');
		push @{$group_members{$group_id}}, $sample_key;
	}

	my @group_states;
	my %group_state_by_id;
	for my $group_id (sort keys %group_members) {
		my @member_keys = @{$group_members{$group_id}};
		my @expected_dirs = map { _normalise_path($map->{$_}{wrdir}) } @member_keys;
		my $first_sample = $member_keys[0];
		my $first_output = $map->{$first_sample}{wrdir};
		my $aim = $groups->{$group_id}{CntAimAss};
		$aim = scalar(@member_keys) unless (defined $aim && $aim > 0);
		my $assembly_dir = $aim > 1
			? _join_path(_base_output_dir($first_output, $first_sample), "AssmblGrp_$group_id", 'metag')
			: _join_path($first_output, 'assemblies', 'metag');
		my $membership_path = _join_path($assembly_dir, 'smpls_used.txt');
		my ($membership_exists, $actual_dirs) = _read_membership_file($membership_path);
		my $membership_matches = $membership_exists ? _same_members(\@expected_dirs, $actual_dirs) : 0;
		my $hybrid = 0;
		if ($assembly_mode == 5) {
			$hybrid = grep { _support_technology($map->{$_}{SupportReads}) =~ m{^(?:PB|ONT)$} } @member_keys;
			$hybrid = $hybrid ? 1 : 0;
		}
		my $assembly = _stage_state(
			artifacts => [_join_path($assembly_dir, 'scaffolds.fasta.filt')],
			markers => [_join_path($assembly_dir, 'ass.done.sto')],
		);
		my $gene_prediction = _alternative_output_state(
			_join_path($assembly_dir, 'genePred', 'proteins.shrtHD.faa'),
			_join_path($assembly_dir, 'genePred', 'proteins.bac.shrtHD.faa'),
		);
		my @issues;
		push @issues, 'GROUP_MEMBERSHIP_MISMATCH' if ($membership_exists && !$membership_matches);
		push @issues, 'GROUP_MEMBERSHIP_MISSING' if ($assembly->{complete} && !$membership_exists);
		my $state = {
			group_id => "$group_id",
			member_keys => \@member_keys,
			expected_output_dirs => \@expected_dirs,
			assembly_dir => $assembly_dir,
			hybrid => $hybrid,
			membership => {
				path => $membership_path,
				exists => $membership_exists ? 1 : 0,
				matches => $membership_matches ? 1 : 0,
				actual_output_dirs => $actual_dirs,
			},
			stages => {
				assembly => $assembly,
				gene_prediction => $gene_prediction,
			},
			issues => \@issues,
		};
		push @group_states, $state;
		$group_state_by_id{$group_id} = $state;
	}

	my @sample_states;
	my $issue_count = 0;
	for my $sample_key (@sample_order) {
		my $sample = $map->{$sample_key};
		my $sample_id = $sample->{SmplID};
		my $output_dir = $sample->{wrdir};
		my $group_id = $sample->{AssGroup};
		$group_id = $sample_key if (!defined $group_id || $group_id eq '-1');
		my $assembly_dir = $group_state_by_id{$group_id}{assembly_dir};
		my $mapping_dir = _join_path($output_dir, 'mapping');
		my $contig_stats_dir = _join_path($output_dir, 'assemblies', 'metag', 'ContigStats');

		my $assembly = $group_state_by_id{$group_id}{stages}{assembly};
		my $hybrid = $group_state_by_id{$group_id}{hybrid};
		my $has_primary_reads = exists($sample->{hasPrimaryRds})
			? ($sample->{hasPrimaryRds} ? 1 : 0) : 1;
		my $support_reads_declared = ($sample->{SupportReads} || '') ne '' ? 1 : 0;
		my $preassembly = $hybrid ? _stage_state(
			artifacts => [_join_path($assembly_dir, 'scaffolds.fasta.filt')],
			markers => [_join_path($assembly_dir, 'preassmblDone.sto')],
		) : _not_applicable_state();
		my $preassembly_package = _not_applicable_state();
		if ($hybrid && $has_primary_reads) {
			if ($run_tmp_dir ne '') {
				my $package_dir = _join_path($run_tmp_dir, $sample_id, "preAssmblGrp_$group_id");
				$preassembly_package = _stage_state(
					artifacts => [_join_path($package_dir, 'scaffolds.fasta.filt')],
					gzip_artifacts => [
						_join_path($package_dir, 'Coverage.percontig'),
						_join_path($package_dir, 'Coverage.median.percontig'),
					],
					markers => [_join_path($package_dir, 'moved.sto')],
				);
			} else {
				$preassembly_package = _unknown_state('INSPECTION_PATH_UNAVAILABLE');
			}
		}
		my $mapping = $map_to_assembly && $has_primary_reads ? _stage_state(
			artifacts => [
				_join_path($mapping_dir, "$sample_id-smd.cram"),
				_join_path($mapping_dir, "$sample_id-smd.bam.coverage.gz"),
			],
			markers => [_join_path($mapping_dir, "$sample_id-smd.cram.sto")],
		) : _not_applicable_state();
		my $support_mapping = $map_support_to_assembly && $support_reads_declared ? _stage_state(
			artifacts => [
				_join_path($mapping_dir, "$sample_id.sup-smd.cram"),
				_join_path($mapping_dir, "$sample_id.sup-smd.bam.coverage.gz"),
			],
			markers => [_join_path($mapping_dir, "$sample_id.sup-smd.cram.sto")],
		) : _not_applicable_state();
		my $coverage = $map_to_assembly && $has_primary_reads ? _stage_state(
			gzip_artifacts => [_join_path($contig_stats_dir, 'Coverage.percontig')],
		) : _not_applicable_state();
		my $support_coverage = $map_support_to_assembly && $support_reads_declared ? _stage_state(
			gzip_artifacts => [_join_path($contig_stats_dir, 'Cov.sup.percontig')],
		) : _not_applicable_state();
		my $gene_prediction = $group_state_by_id{$group_id}{stages}{gene_prediction};
		my $empty_marker = _file_state(_join_path($output_dir, 'SMPL.empty'), 0);
		my @issues = (
			@{$assembly->{issues}},
			@{$preassembly_package->{issues}},
			@{$mapping->{issues}},
			@{$support_mapping->{issues}},
			@{$gene_prediction->{issues}},
		);
		$issue_count += scalar(@issues);

		push @sample_states, {
			sample_key => $sample_key,
			sample_id => $sample_id,
			assembly_group => "$group_id",
			output_dir => $output_dir,
			sample_empty => $empty_marker->{exists} ? 1 : 0,
			has_primary_reads => $has_primary_reads,
			support_reads_declared => $support_reads_declared,
			support_read_technology => _support_technology($sample->{SupportReads}),
			stages => {
				assembly => $assembly,
				preassembly => $preassembly,
				preassembly_package => $preassembly_package,
				mapping => $mapping,
				support_mapping => $support_mapping,
				coverage => $coverage,
				support_coverage => $support_coverage,
				gene_prediction => $gene_prediction,
			},
			issues => \@issues,
		};
	}

	$issue_count += scalar(@{$_->{issues}}) for @group_states;
	return {
		schema_version => 1,
		mode => 'inspection',
		read_only => 1,
		workflow => {
			assembly_mode => $assembly_mode,
			assembly_requested => $assembly_requested,
			map_to_assembly => $map_to_assembly,
			map_support_to_assembly => $map_support_to_assembly,
			run_tmp_dir => $run_tmp_dir,
		},
		summary => {
			samples => scalar(@sample_states),
			assembly_groups => scalar(@group_states),
			issues => $issue_count,
		},
		samples => \@sample_states,
		assembly_groups => \@group_states,
	};
}

sub encode_state_report {
	my ($report) = @_;
	return JSON::PP->new->canonical(1)->pretty(1)->encode($report);
}

1;
