package Mods::WorkflowRunner;

use warnings;
use strict;

use Exporter qw(import);

use Mods::WorkflowApply qw(apply_workflow_plan);
use Mods::WorkflowPlan qw(build_workflow_plan);
use Mods::WorkflowState qw(inspect_workflow_state);

our @EXPORT_OK = qw(run_workflow_preflight);

sub _allowed_roots {
	my ($state) = @_;
	my @roots = map { $_->{output_dir} } @{$state->{samples} || []};
	push @roots, map { $_->{assembly_dir} } @{$state->{assembly_groups} || []};
	push @roots, $state->{workflow}{run_tmp_dir}
		if (($state->{workflow}{run_tmp_dir} || '') ne '');
	my %seen;
	return [grep { defined $_ && $_ ne '' && !$seen{$_}++ } @roots];
}

sub run_workflow_preflight {
	my (%args) = @_;
	my $state_before = inspect_workflow_state(
		map => $args{map},
		groups => $args{groups},
		options => $args{options} || {},
	);
	my $plan_before = build_workflow_plan($state_before);
	my $repairs = apply_workflow_plan(
		plan => $plan_before,
		allowed_roots => _allowed_roots($state_before),
		dry_run => $args{apply_repairs} ? 0 : 1,
		allow_group_rewrite => $args{allow_group_rewrite} ? 1 : 0,
	);
	my $state_after = inspect_workflow_state(
		map => $args{map},
		groups => $args{groups},
		options => $args{options} || {},
	);
	my $plan_after = build_workflow_plan($state_after);
	$plan_after->{automatic_preflight} = $repairs;
	$plan_after->{iteration} = 0 + ($args{iteration} || 0);
	return {
		iteration => 0 + ($args{iteration} || 0),
		state_before => $state_before,
		plan_before => $plan_before,
		repairs => $repairs,
		state => $state_after,
		plan => $plan_after,
	};
}

1;
