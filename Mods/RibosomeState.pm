package Mods::RibosomeState;

use strict;
use warnings;

use Exporter qw(import);
use File::Glob qw(bsd_glob GLOB_NOSORT);
use File::Path qw(remove_tree);
use File::Spec;

use Mods::SampleCompletion qw(invalidate_sample_completion);

our @EXPORT_OK = qw(
	normalise_ribosome_request
	prepare_ribosome_rerun
	ribosome_completion_evidence
);

sub _result_present {
	my ($path) = @_;
	return (-e $path || -e "$path.gz") ? 1 : 0;
}

sub _profile_outputs_present {
	my ($root, $tag) = @_;
	for my $suffix ('r1.fq.gz', 'r2.fq.gz', 'fq.gz') {
		return 0 unless -s File::Spec->catfile(
			$root, "reads_$tag.$suffix",
		);
	}
	return 1;
}

sub ribosome_completion_evidence {
	my (%args) = @_;
	my $requested = $args{requested} ? 1 : 0;
	return {requested => 0, complete => 1} unless $requested;

	my $ribo_root = $args{ribo_root};
	if (!defined($ribo_root) || $ribo_root eq '') {
		my $sample_root = $args{sample_root};
		die "ribosome_completion_evidence requires sample_root or ribo_root\n"
			unless defined($sample_root) && length($sample_root);
		$ribo_root = File::Spec->catdir($sample_root, 'ribos');
	}
	my $lca_root = File::Spec->catdir($ribo_root, 'ltsLCA');
	my $assembly_requested = $args{assembly_requested} ? 1 : 0;
	my %evidence = (
		requested => 1,
		assembly_requested => $assembly_requested,
		ssu_profile_complete => -e File::Spec->catfile($ribo_root, 'SSU_pull.sto') ? 1 : 0,
		lsu_profile_complete => -e File::Spec->catfile($ribo_root, 'LSU_pull.sto') ? 1 : 0,
		ssu_profile_outputs_complete => _profile_outputs_present($ribo_root, 'SSU'),
		lsu_profile_outputs_complete => _profile_outputs_present($ribo_root, 'LSU'),
		assembly_complete => (!$assembly_requested
			|| -e File::Spec->catfile($ribo_root, 'Ass', 'allAss.sto')) ? 1 : 0,
		assignment_complete_stone => -e File::Spec->catfile($lca_root, 'Assigned.sto') ? 1 : 0,
		ssu_assignment_complete => -e File::Spec->catfile($lca_root, 'SSU_ass.sto') ? 1 : 0,
		lsu_assignment_complete => -e File::Spec->catfile($lca_root, 'LSU_ass.sto') ? 1 : 0,
		ssu_hierarchy_complete => _result_present(
			File::Spec->catfile($lca_root, 'SSUriboRun_bl.hiera.txt'),
		),
		lsu_hierarchy_complete => _result_present(
			File::Spec->catfile($lca_root, 'LSUriboRun_bl.hiera.txt'),
		),
	);
	$evidence{profile_complete} = $evidence{ssu_profile_complete}
		&& $evidence{lsu_profile_complete}
		&& $evidence{ssu_profile_outputs_complete}
		&& $evidence{lsu_profile_outputs_complete}
		&& $evidence{assembly_complete} ? 1 : 0;
	$evidence{taxonomy_complete} = $evidence{assignment_complete_stone}
		&& $evidence{ssu_assignment_complete} && $evidence{lsu_assignment_complete}
		&& $evidence{ssu_hierarchy_complete} && $evidence{lsu_hierarchy_complete} ? 1 : 0;
	$evidence{complete} = $evidence{profile_complete}
		&& $evidence{taxonomy_complete} ? 1 : 0;
	return \%evidence;
}

sub normalise_ribosome_request {
	my ($options) = @_;
	die "normalise_ribosome_request requires an option hash\n"
		unless ref($options) eq 'HASH';
	$options->{DoRibofind} = 1
		if $options->{RedoRiboFind} || $options->{RedoRiboAssign};
	return $options->{DoRibofind} ? 1 : 0;
}

sub _remove_tree_checked {
	my ($path) = @_;
	return 0 unless -e $path || -l $path;
	my $errors;
	remove_tree($path, {error => \$errors});
	if ($errors && @{$errors}) {
		my @messages;
		for my $record (@{$errors}) {
			for my $failed_path (keys %{$record}) {
				my $message = $record->{$failed_path};
				push @messages, ($failed_path || $path).": $message";
			}
		}
		die "Cannot remove RiboFind results: ".join('; ', @messages)."\n";
	}
	return 1;
}

sub _unlink_if_present {
	my ($path) = @_;
	return 0 unless -e $path || -l $path;
	unlink $path or die "Cannot remove stale RiboFind result $path: $!\n";
	return 1;
}

sub prepare_ribosome_rerun {
	my (%args) = @_;
	my $redo_profile = $args{redo_profile} ? 1 : 0;
	my $redo_assignment = $args{redo_assignment} ? 1 : 0;
	return {removed => 0, profile => 0, assignment => 0}
		unless $redo_profile || $redo_assignment;

	my $sample_root = $args{sample_root};
	my $central_root = $args{central_root};
	my $sample = $args{sample};
	die "prepare_ribosome_rerun requires sample_root\n"
		unless defined($sample_root) && length($sample_root);
	die "prepare_ribosome_rerun requires central_root\n"
		unless defined($central_root) && length($central_root);
	die "prepare_ribosome_rerun requires sample\n"
		unless defined($sample) && length($sample);

	# Invalidate before MATAF4 assesses completion. Otherwise it can remember
	# the old stones as complete, remove them, and submit no replacement work.
	invalidate_sample_completion($sample_root);
	my $ribo_root = File::Spec->catdir($sample_root, 'ribos');
	my $removed = $redo_profile
		? _remove_tree_checked($ribo_root)
		: _remove_tree_checked(File::Spec->catdir($ribo_root, 'ltsLCA'));

	# Central links and merged tables derive from ltsLCA output. Invalidate them
	# so equal file sizes or an unchanged sample count cannot preserve old data.
	for my $tag (qw(SSU LSU)) {
		my $sample_result = File::Spec->catfile(
			$central_root, $tag, "$sample.$tag.hiera.txt",
		);
		$removed += _unlink_if_present($sample_result);
		$removed += _unlink_if_present("$sample_result.gz");

		my @aggregate_results = bsd_glob(
			File::Spec->catfile($central_root, "$tag.miTag*"),
			GLOB_NOSORT,
		);
		push @aggregate_results,
			File::Spec->catfile($central_root, "$tag.cnt.stone");
		for my $path (@aggregate_results) {
			next unless -f $path || -l $path;
			$removed += _unlink_if_present($path);
		}
	}

	return {
		removed => $removed,
		profile => $redo_profile,
		assignment => $redo_profile || $redo_assignment,
	};
}

1;
