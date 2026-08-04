package Mods::RibosomeState;

use strict;
use warnings;

use Exporter qw(import);
use File::Glob qw(bsd_glob GLOB_NOSORT);
use File::Path qw(remove_tree);
use File::Spec;

use Mods::SampleCompletion qw(invalidate_sample_completion);

our @EXPORT_OK = qw(prepare_ribosome_rerun);

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
