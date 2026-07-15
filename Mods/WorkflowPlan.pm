package Mods::WorkflowPlan;

use warnings;
use strict;

use Exporter qw(import);
use JSON::PP;

our @EXPORT_OK = qw(build_workflow_plan encode_workflow_plan validate_workflow_plan);

sub _existing_paths {
	my ($stage) = @_;
	my @paths;
	for my $entry (@{$stage->{artifacts} || []}, @{$stage->{markers} || []}) {
		push @paths, $entry->{path} if ($entry->{exists} && defined $entry->{path});
	}
	return @paths;
}

sub _expected_paths {
	my ($stage) = @_;
	my @paths;
	for my $entry (@{$stage->{artifacts} || []}, @{$stage->{markers} || []}) {
		push @paths, $entry->{requested_path} if (defined $entry->{requested_path});
	}
	return @paths;
}

sub _push_unique {
	my ($array, @values) = @_;
	my %seen = map { $_ => 1 } @{$array};
	push @{$array}, grep { defined $_ && $_ ne '' && !$seen{$_}++ } @values;
}

sub _add_action {
	my ($actions, $by_id, $action) = @_;
	die "Duplicate workflow plan action id: $action->{id}" if (exists $by_id->{$action->{id}});
	$action->{depends_on} ||= [];
	$action->{reason_codes} ||= [];
	$action->{targets} ||= [];
	$action->{automatic_targets} ||= [];
	$action->{auto_apply} = $action->{auto_apply} ? 1 : 0;
	$action->{automatic_policy} ||= 'none';
	$action->{expected_outputs} ||= [];
	$action->{requires_confirmation} = $action->{requires_confirmation} ? 1 : 0;
	push @{$actions}, $action;
	$by_id->{$action->{id}} = $action;
	return $action->{id};
}

sub _partial_stage_repair {
	my (%args) = @_;
	my $stage = $args{stage};
	return '' unless ($stage->{status} eq 'PARTIAL' && @{$stage->{issues} || []});
	my @targets = _existing_paths($stage);
	my %automatic_stage = map { $_ => 1 } qw(
		mapping support_mapping coverage support_coverage preassembly_package
	);
	my @automatic_targets;
	if ($automatic_stage{$args{stage_name}}) {
		@automatic_targets = @targets;
	} else {
		my %issues = map { $_ => 1 } @{$stage->{issues}};
		if ($issues{MARKER_WITHOUT_VALID_OUTPUT}) {
			push @automatic_targets, map { $_->{path} }
				grep { $_->{exists} && defined $_->{path} } @{$stage->{markers} || []};
		}
		push @automatic_targets, map { $_->{path} }
			grep { $_->{exists} && !$_->{nonempty} && defined $_->{path} }
			@{$stage->{artifacts} || []};
	}
	my $automatic = @automatic_targets ? 1 : 0;
	my $fully_automatic = $automatic && @automatic_targets == @targets ? 1 : 0;
	return _add_action($args{actions}, $args{by_id}, {
		id => $args{id},
		kind => 'repair',
		operation => 'invalidate_partial_stage',
		scope => $args{scope},
		stage => $args{stage_name},
		depends_on => $args{depends_on} || [],
		reason_codes => [@{$stage->{issues}}],
		targets => \@targets,
		automatic_targets => \@automatic_targets,
		auto_apply => $automatic,
		automatic_policy => $automatic ? 'invalidate_incomplete_stage' : 'none',
		expected_outputs => [],
		requires_confirmation => $fully_automatic ? 0 : 1,
		authorization => $fully_automatic ? 'automatic_safe_repair' : 'explicit_apply',
		risk => 'destructive',
	});
}

sub _stage_needs_work {
	my ($stage) = @_;
	return 0 if (!defined $stage || $stage->{status} eq 'NOT_APPLICABLE');
	return $stage->{complete} ? 0 : 1;
}

sub _stage_applicable {
	my ($stage) = @_;
	return defined($stage) && $stage->{status} ne 'NOT_APPLICABLE';
}

sub build_workflow_plan {
	my ($state) = @_;
	die 'build_workflow_plan requires an inspection report' unless (ref($state) eq 'HASH');
	die 'workflow plan requires a read-only inspection report'
		unless ($state->{mode} eq 'inspection' && $state->{read_only});

	my @actions;
	my %by_id;
	my %group_submit;
	my %group_repair;
	my %sample_by_key = map { $_->{sample_key} => $_ } @{$state->{samples} || []};
	my $workflow = $state->{workflow} || {};

	for my $group (@{$state->{assembly_groups} || []}) {
		my $group_id = $group->{group_id};
		my $scope = { type => 'assembly_group', id => $group_id };
		my @repair_dependencies;

		my @membership_issues = grep {
			$_ eq 'GROUP_MEMBERSHIP_MISMATCH' || $_ eq 'GROUP_MEMBERSHIP_MISSING'
		} @{$group->{issues} || []};
		if (@membership_issues) {
			my @targets = ($group->{assembly_dir});
			for my $sample_key (@{$group->{member_keys}}) {
				my $sample = $sample_by_key{$sample_key};
				if ($sample) {
					push @targets, "$sample->{output_dir}/mapping";
					push @targets, "$sample->{output_dir}/assemblies/metag/ContigStats";
				}
			}
			my $id = _add_action(\@actions, \%by_id, {
				id => "group:$group_id:repair:membership",
				kind => 'repair',
				operation => 'invalidate_changed_assembly_group',
				scope => $scope,
				stage => 'assembly',
				reason_codes => \@membership_issues,
				targets => \@targets,
				automatic_targets => [],
				auto_apply => 0,
				automatic_policy => 'none',
				requires_confirmation => 1,
				authorization => 'OKtoRWassGrps',
				risk => 'group_wide_destructive',
			});
			push @repair_dependencies, $id;
		}

		my $assembly_repair = _partial_stage_repair(
			actions => \@actions,
			by_id => \%by_id,
			id => "group:$group_id:repair:assembly",
			scope => $scope,
			stage_name => 'assembly',
			stage => $group->{stages}{assembly},
			depends_on => [@repair_dependencies],
		);
		push @repair_dependencies, $assembly_repair if ($assembly_repair ne '');
		$group_repair{$group_id} = [@repair_dependencies];
		my $assembly_needs_work = _stage_needs_work($group->{stages}{assembly})
			|| @repair_dependencies;

		unless ($workflow->{assembly_requested}) {
			if ($assembly_needs_work) {
				$group_submit{$group_id} = _add_action(\@actions, \%by_id, {
					id => "group:$group_id:manual:assembly_required",
					kind => 'manual',
					operation => 'provide_existing_assembly_or_enable_assembly',
					scope => $scope,
					stage => 'assembly',
					depends_on => [@repair_dependencies],
					reason_codes => ['ASSEMBLY_REQUIRED_BUT_NOT_REQUESTED'],
					expected_outputs => [_expected_paths($group->{stages}{assembly})],
					requires_confirmation => 0,
					authorization => 'user_decision',
					risk => 'blocking',
				});
			}
			next;
		}
		if ($assembly_needs_work) {
			my @assembly_dependencies = @repair_dependencies;
			if ($group->{hybrid}) {
				my $previous_preassembly = '';
				for my $sample_key (@{$group->{member_keys}}) {
					my $sample = $sample_by_key{$sample_key};
					next unless ($sample && $sample->{has_primary_reads} && !$sample->{sample_empty});
					my $package = $sample->{stages}{preassembly_package};
					my @pre_dependencies = @repair_dependencies;
					my $repair = _partial_stage_repair(
						actions => \@actions,
						by_id => \%by_id,
						id => "sample:$sample_key:repair:preassembly_package",
						scope => { type => 'sample', id => $sample_key },
						stage_name => 'preassembly_package',
						stage => $package,
						depends_on => [@pre_dependencies],
					);
					push @pre_dependencies, $repair if ($repair ne '');
					push @pre_dependencies, $previous_preassembly if ($previous_preassembly ne '');
					if (_stage_needs_work($package)) {
						$previous_preassembly = _add_action(\@actions, \%by_id, {
							id => "sample:$sample_key:submit:preassembly",
							kind => 'submit',
							operation => 'submit_preassembly_and_package',
							scope => { type => 'sample', id => $sample_key },
							stage => 'preassembly_package',
							depends_on => \@pre_dependencies,
							reason_codes => [$package->{status} eq 'UNKNOWN' ? 'STATE_UNINSPECTABLE' : 'STAGE_INCOMPLETE'],
							expected_outputs => [_expected_paths($package)],
							requires_confirmation => 0,
							authorization => 'scheduler_submit',
							risk => 'job_submission',
						});
					}
					push @assembly_dependencies, $previous_preassembly if ($previous_preassembly ne '');
				}
			}

			$group_submit{$group_id} = _add_action(\@actions, \%by_id, {
				id => "group:$group_id:submit:assembly",
				kind => 'submit',
				operation => $group->{hybrid} ? 'submit_hybrid_group_assembly' : 'submit_group_assembly',
				scope => $scope,
				stage => 'assembly',
				depends_on => \@assembly_dependencies,
				reason_codes => ['STAGE_INCOMPLETE'],
				expected_outputs => [_expected_paths($group->{stages}{assembly})],
				requires_confirmation => 0,
				authorization => 'scheduler_submit',
				risk => 'job_submission',
			});
		}

		if (_stage_needs_work($group->{stages}{gene_prediction}) || @repair_dependencies) {
			my @dependencies = @repair_dependencies;
			my $gene_repair = _partial_stage_repair(
				actions => \@actions,
				by_id => \%by_id,
				id => "group:$group_id:repair:gene_prediction",
				scope => $scope,
				stage_name => 'gene_prediction',
				stage => $group->{stages}{gene_prediction},
				depends_on => [@repair_dependencies],
			);
			push @dependencies, $gene_repair if ($gene_repair ne '');
			push @dependencies, $group_submit{$group_id} if ($group_submit{$group_id});
			$group_submit{"$group_id:gene_prediction"} = _add_action(\@actions, \%by_id, {
				id => "group:$group_id:submit:gene_prediction",
				kind => 'submit',
				operation => 'submit_gene_prediction',
				scope => $scope,
				stage => 'gene_prediction',
				depends_on => \@dependencies,
				reason_codes => ['STAGE_INCOMPLETE'],
				expected_outputs => [_expected_paths($group->{stages}{gene_prediction})],
				requires_confirmation => 0,
				authorization => 'scheduler_submit',
				risk => 'job_submission',
			});
		}
	}

	for my $sample (@{$state->{samples} || []}) {
		my $sample_key = $sample->{sample_key};
		my $group_id = $sample->{assembly_group};
		my $scope = { type => 'sample', id => $sample_key };
		my @group_repair_dependencies = @{$group_repair{$group_id} || []};
		my $group_invalidated = @group_repair_dependencies ? 1 : 0;
		my @assembly_dependencies = @group_repair_dependencies;
		push @assembly_dependencies, $group_submit{$group_id} if ($group_submit{$group_id});
		my %submit;
		if ($sample->{sample_empty}) {
			_add_action(\@actions, \%by_id, {
				id => "sample:$sample_key:manual:empty",
				kind => 'manual',
				operation => 'confirm_empty_sample_skip',
				scope => $scope,
				stage => 'sample_input',
				depends_on => [],
				reason_codes => ['SAMPLE_MARKED_EMPTY'],
				requires_confirmation => 0,
				authorization => 'user_review',
				risk => 'informational',
			});
			next;
		}

		for my $stage_name (qw(mapping support_mapping)) {
			my $stage = $sample->{stages}{$stage_name};
			next unless (_stage_applicable($stage));
			next unless (_stage_needs_work($stage) || $group_invalidated);
			my @repair_dependencies = @group_repair_dependencies;
			my $repair = _partial_stage_repair(
				actions => \@actions,
				by_id => \%by_id,
				id => "sample:$sample_key:repair:$stage_name",
				scope => $scope,
				stage_name => $stage_name,
				stage => $stage,
				depends_on => [@repair_dependencies],
			);
			my @dependencies = @assembly_dependencies;
			push @dependencies, $repair if ($repair ne '');
			$submit{$stage_name} = _add_action(\@actions, \%by_id, {
				id => "sample:$sample_key:submit:$stage_name",
				kind => 'submit',
				operation => $stage_name eq 'mapping' ? 'submit_primary_mapping' : 'submit_support_mapping',
				scope => $scope,
				stage => $stage_name,
				depends_on => \@dependencies,
				reason_codes => ['STAGE_INCOMPLETE'],
				expected_outputs => [_expected_paths($stage)],
				requires_confirmation => 0,
				authorization => 'scheduler_submit',
				risk => 'job_submission',
			});
		}

		my @coverage_stages = grep {
			_stage_applicable($sample->{stages}{$_})
				&& (_stage_needs_work($sample->{stages}{$_}) || $group_invalidated)
		} qw(coverage support_coverage);
		if (@coverage_stages) {
			my @coverage_repairs;
			for my $stage_name (@coverage_stages) {
				my $repair = _partial_stage_repair(
					actions => \@actions,
					by_id => \%by_id,
					id => "sample:$sample_key:repair:$stage_name",
					scope => $scope,
					stage_name => $stage_name,
					stage => $sample->{stages}{$stage_name},
					depends_on => [@group_repair_dependencies],
				);
				push @coverage_repairs, $repair if ($repair ne '');
			}
			my @dependencies = (@assembly_dependencies, @coverage_repairs);
			push @dependencies, $group_submit{"$group_id:gene_prediction"}
				if ($group_submit{"$group_id:gene_prediction"});
			push @dependencies, grep { defined $_ && $_ ne '' } @submit{qw(mapping support_mapping)};
			my @expected;
			_push_unique(\@expected, _expected_paths($sample->{stages}{$_})) for @coverage_stages;
			_add_action(\@actions, \%by_id, {
				id => "sample:$sample_key:submit:contig_stats",
				kind => 'submit',
				operation => 'submit_contig_stats',
				scope => $scope,
				stage => 'contig_stats',
				depends_on => \@dependencies,
				reason_codes => ['COVERAGE_INCOMPLETE'],
				expected_outputs => \@expected,
				requires_confirmation => 0,
				authorization => 'scheduler_submit',
				risk => 'job_submission',
			});
		}
	}

	my $repair_count = grep { $_->{kind} eq 'repair' } @actions;
	my $submit_count = grep { $_->{kind} eq 'submit' } @actions;
	my $manual_count = grep { $_->{kind} eq 'manual' } @actions;
	my $automatic_repair_count = grep { $_->{kind} eq 'repair' && $_->{auto_apply} } @actions;
	my $plan = {
		schema_version => 1,
		mode => 'repair_submission_plan',
		read_only => 1,
		execution_supported => 0,
		automatic_repairs_supported => 1,
		summary => {
			actions => scalar(@actions),
			repairs => 0 + $repair_count,
			submissions => 0 + $submit_count,
			manual_steps => 0 + $manual_count,
			automatic_repairs => 0 + $automatic_repair_count,
			confirmations_required => 0 + grep { $_->{requires_confirmation} } @actions,
		},
		actions => \@actions,
		source_state => $state,
	};
	my @validation_errors = validate_workflow_plan($plan);
	die 'Invalid workflow plan: '.join('; ', @validation_errors) if (@validation_errors);
	return $plan;
}

sub validate_workflow_plan {
	my ($plan) = @_;
	my @errors;
	my %actions = map { $_->{id} => $_ } @{$plan->{actions} || []};
	for my $action (@{$plan->{actions} || []}) {
		push @errors, "action without id" unless (defined $action->{id} && $action->{id} ne '');
		for my $dependency (@{$action->{depends_on} || []}) {
			push @errors, "$action->{id} depends on missing action $dependency"
				unless (exists $actions{$dependency});
		}
	}

	my (%visiting, %visited);
	my $visit;
	$visit = sub {
		my ($id) = @_;
		return if ($visited{$id});
		return unless (exists $actions{$id});
		if ($visiting{$id}) {
			push @errors, "dependency cycle at $id";
			return;
		}
		$visiting{$id} = 1;
		$visit->($_) for @{$actions{$id}{depends_on} || []};
		delete $visiting{$id};
		$visited{$id} = 1;
	};
	$visit->($_) for sort keys %actions;
	return @errors;
}

sub encode_workflow_plan {
	my ($plan) = @_;
	return JSON::PP->new->canonical(1)->pretty(1)->encode($plan);
}

1;
