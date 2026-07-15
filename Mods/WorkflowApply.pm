package Mods::WorkflowApply;

use warnings;
use strict;

use Exporter qw(import);
use File::Path qw(remove_tree);

our @EXPORT_OK = qw(apply_workflow_plan);

sub _normalise_path {
	my ($path) = @_;
	$path = '' unless (defined $path);
	$path =~ s{\\}{/}g;
	$path =~ s{/+}{/}g;
	$path =~ s{/$}{};
	return $path;
}

sub _safe_path {
	my ($path, $roots) = @_;
	$path = _normalise_path($path);
	return 0 if ($path eq '' || $path eq '/' || $path =~ m{(?:^|/)\.\.(?:/|$)});
	for my $root (@{$roots}) {
		$root = _normalise_path($root);
		next if ($root eq '' || $root eq '/');
		return 1 if ($path eq $root || index($path, "$root/") == 0);
	}
	return 0;
}

sub _remove_target {
	my (%args) = @_;
	my $path = $args{path};
	die "Unsafe workflow repair target: $path"
		unless (_safe_path($path, $args{allowed_roots}));
	return 'missing' unless (-e $path || -l $path);
	return 'would_remove' if ($args{dry_run});
	if (-d $path && !-l $path) {
		die "Directory repair requires group rewrite authorization: $path"
			unless ($args{allow_directory});
		my $errors = [];
		remove_tree($path, { error => \$errors });
		die "Failed to remove workflow repair directory $path"
			if (@{$errors});
	} else {
		unlink($path) or die "Failed to remove workflow repair target $path: $!";
	}
	return 'removed';
}

sub apply_workflow_plan {
	my (%args) = @_;
	my $plan = $args{plan} || die 'apply_workflow_plan requires plan';
	my @allowed_roots = map { _normalise_path($_) } @{$args{allowed_roots} || []};
	die 'apply_workflow_plan requires allowed_roots' unless (@allowed_roots);
	my $dry_run = $args{dry_run} ? 1 : 0;
	my $allow_group_rewrite = $args{allow_group_rewrite} ? 1 : 0;
	my %result_by_id;
	my @results;
	my ($removed, $would_remove, $blocked) = (0, 0, 0);

	for my $action (@{$plan->{actions} || []}) {
		next unless ($action->{kind} eq 'repair');
		my @repair_dependencies = grep { exists $result_by_id{$_} } @{$action->{depends_on} || []};
		my @failed_dependencies = grep {
			$result_by_id{$_}{status} eq 'blocked' || $result_by_id{$_}{status} eq 'skipped'
		} @repair_dependencies;
		my $result = {
			id => $action->{id},
			operation => $action->{operation},
			status => 'skipped',
			targets => [],
		};
		if (@failed_dependencies) {
			$result->{reason} = 'repair_dependency_not_applied';
			$blocked++;
		} else {
			my $group_rewrite = $action->{operation} eq 'invalidate_changed_assembly_group';
			my $automatic = $action->{auto_apply} ? 1 : 0;
			my $authorised_group = $group_rewrite && $allow_group_rewrite;
			if ($automatic || $authorised_group) {
				my @targets = $authorised_group
					? @{$action->{targets} || []}
					: @{$action->{automatic_targets} || []};
				$result->{status} = $dry_run ? 'dry_run' : 'applied';
				for my $target (@targets) {
					my $status = _remove_target(
						path => $target,
						allowed_roots => \@allowed_roots,
						dry_run => $dry_run,
						allow_directory => $authorised_group,
					);
					push @{$result->{targets}}, { path => $target, status => $status };
					$removed++ if ($status eq 'removed');
					$would_remove++ if ($status eq 'would_remove');
				}
			} else {
				$result->{status} = 'blocked';
				$result->{reason} = $group_rewrite
					? 'OKtoRWassGrps_required' : 'manual_repair_required';
				$blocked++;
			}
		}
		push @results, $result;
		$result_by_id{$action->{id}} = $result;
	}

	return {
		mode => $dry_run ? 'repair_dry_run' : 'automatic_repair',
		removed_targets => $removed,
		would_remove_targets => $would_remove,
		blocked_repairs => $blocked,
		actions => \@results,
	};
}

1;
