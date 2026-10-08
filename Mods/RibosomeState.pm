package Mods::RibosomeState;

use strict;
use warnings;

use Exporter qw(import);
use File::Basename qw(basename);
use File::Glob qw(bsd_glob GLOB_NOSORT);
use File::Path qw(remove_tree);
use File::Spec;

use Mods::SampleCompletion qw(
	completion_component_evidence completion_request_signature
	invalidate_sample_completion
);

our @EXPORT_OK = qw(
	compress_ribosome_hierarchies
	lca_reference_copy_current
	normalise_ribosome_request
	prepare_ribosome_rerun
	ribosome_completion_evidence
	ribosome_merge_cohort
	ribosome_merge_manifest
);
sub ribosome_completion_evidence {
	my (%args) = @_;
	my $requested = $args{requested} ? 1 : 0;
	return completion_component_evidence(requested => 0) unless $requested;

	my $ribo_root = $args{ribo_root};
	if (!defined($ribo_root) || $ribo_root eq '') {
		my $sample_root = $args{sample_root};
		die "ribosome_completion_evidence requires sample_root or ribo_root\n"
			unless defined($sample_root) && length($sample_root);
		$ribo_root = File::Spec->catdir($sample_root, 'ribos');
	}
	my $lca_root = File::Spec->catdir($ribo_root, 'ltsLCA');
	my $assembly_requested = $args{assembly_requested} ? 1 : 0;
	my @checks = (
		{id => 'ssu_profile_stone', kind => 'exists',
			path => File::Spec->catfile($ribo_root, 'SSU_pull.sto')},
		{id => 'lsu_profile_stone', kind => 'exists',
			path => File::Spec->catfile($ribo_root, 'LSU_pull.sto')},
		(map {
			+{id => "ssu_profile_$_", kind => 'nonempty',
				path => File::Spec->catfile($ribo_root, "reads_SSU.$_")}
		} ('r1.fq.gz', 'r2.fq.gz', 'fq.gz')),
		(map {
			+{id => "lsu_profile_$_", kind => 'nonempty',
				path => File::Spec->catfile($ribo_root, "reads_LSU.$_")}
		} ('r1.fq.gz', 'r2.fq.gz', 'fq.gz')),
		{id => 'assignment_stone', kind => 'exists',
			path => File::Spec->catfile($lca_root, 'Assigned.sto')},
		{id => 'ssu_assignment_stone', kind => 'exists',
			path => File::Spec->catfile($lca_root, 'SSU_ass.sto')},
		{id => 'lsu_assignment_stone', kind => 'exists',
			path => File::Spec->catfile($lca_root, 'LSU_ass.sto')},
		{id => 'ssu_hierarchy', kind => 'exists_any', paths => [
			File::Spec->catfile($lca_root, 'SSUriboRun_bl.hiera.txt'),
			File::Spec->catfile($lca_root, 'SSUriboRun_bl.hiera.txt.gz'),
		]},
		{id => 'lsu_hierarchy', kind => 'exists_any', paths => [
			File::Spec->catfile($lca_root, 'LSUriboRun_bl.hiera.txt'),
			File::Spec->catfile($lca_root, 'LSUriboRun_bl.hiera.txt.gz'),
		]},
	);
	push @checks, {
		id => 'ribosomal_assembly_stone', kind => 'exists',
		path => File::Spec->catfile($ribo_root, 'Ass', 'allAss.sto'),
	} if $assembly_requested;

	my $evidence = completion_component_evidence(
		requested => 1,
		checks => \@checks,
	);
	my $ok = sub {
		my ($id) = @_;
		return $evidence->{checks}{$id}{ok} ? 1 : 0;
	};
	$evidence->{assembly_requested} = $assembly_requested;
	$evidence->{ssu_profile_complete} = $ok->('ssu_profile_stone');
	$evidence->{lsu_profile_complete} = $ok->('lsu_profile_stone');
	$evidence->{ssu_profile_outputs_complete} =
		!grep { !$ok->("ssu_profile_$_") } ('r1.fq.gz', 'r2.fq.gz', 'fq.gz');
	$evidence->{lsu_profile_outputs_complete} =
		!grep { !$ok->("lsu_profile_$_") } ('r1.fq.gz', 'r2.fq.gz', 'fq.gz');
	$evidence->{assembly_complete} = !$assembly_requested
		|| $ok->('ribosomal_assembly_stone') ? 1 : 0;
	$evidence->{assignment_complete_stone} = $ok->('assignment_stone');
	$evidence->{ssu_assignment_complete} = $ok->('ssu_assignment_stone');
	$evidence->{lsu_assignment_complete} = $ok->('lsu_assignment_stone');
	$evidence->{ssu_hierarchy_complete} = $ok->('ssu_hierarchy');
	$evidence->{lsu_hierarchy_complete} = $ok->('lsu_hierarchy');
	$evidence->{profile_complete} = $evidence->{ssu_profile_complete}
		&& $evidence->{lsu_profile_complete}
		&& $evidence->{ssu_profile_outputs_complete}
		&& $evidence->{lsu_profile_outputs_complete}
		&& $evidence->{assembly_complete} ? 1 : 0;
	$evidence->{taxonomy_complete} = $evidence->{assignment_complete_stone}
		&& $evidence->{ssu_assignment_complete} && $evidence->{lsu_assignment_complete}
		&& $evidence->{ssu_hierarchy_complete} && $evidence->{lsu_hierarchy_complete} ? 1 : 0;
	$evidence->{complete} = $evidence->{profile_complete}
		&& $evidence->{taxonomy_complete} ? 1 : 0;
	$evidence->{status} = $evidence->{complete} ? 'complete' : 'incomplete';
	return $evidence;
}

sub normalise_ribosome_request {
	my ($options) = @_;
	die "normalise_ribosome_request requires an option hash\n"
		unless ref($options) eq 'HASH';
	# catchLSUSSU rejects assembly, and nothing writes the Ass/allAss.sto that
	# completion would then require: every sample would fail on every pass.
	die "-riobsomalAssembly is no longer supported; remove it or set it to 0\n"
		if $options->{doRiboAssembl};
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

	# Merged tables derive from ltsLCA output: remove them with their merge
	# signature, and the per-sample links that versions before 2026-10-08 kept.
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

# Compresses each marker's hierarchy in the sample directory, where the cohort
# merge reads it. Returns the markers without a hierarchy.
sub compress_ribosome_hierarchies {
	my (%args) = @_;
	die "compress_ribosome_hierarchies requires sample_root\n"
		unless defined($args{sample_root}) && length($args{sample_root});
	my @missing;
	for my $tag (qw(SSU LSU)) {
		my $source = File::Spec->catfile(
			$args{sample_root}, 'ribos', 'ltsLCA', "${tag}riboRun_bl.hiera.txt",
		);
		# -f: a .gz beside the plain file is left by an interrupted gzip
		if (-e $source) {
			system('gzip', '-f', $source) == 0
				or die "Cannot compress RiboFind hierarchy $source\n";
		}
		push @missing, $tag unless -e "$source.gz";
	}
	return \@missing;
}

# Sorts the samples of a map (samples => [{name, root}, ...]) for the SSU/LSU
# merge: "cohort" holds those with a complete profile and assignment, each with
# its hierarchy paths; "incomplete" names the others. Samples that excluded->()
# accepts (no reads, skipped, ignored) are in neither list.
sub ribosome_merge_cohort {
	my (%args) = @_;
	die "ribosome_merge_cohort requires a sample list\n"
		unless ref($args{samples}) eq 'ARRAY';
	my (@cohort, @incomplete);
	for my $sample (@{$args{samples}}) {
		my $evidence = ribosome_completion_evidence(
			sample_root => $sample->{root}, requested => 1,
		);
		if ($evidence->{complete}) {
			push @cohort, {
				%{$sample},
				hierarchy => {
					SSU => $evidence->{checks}{ssu_hierarchy}{matched_path},
					LSU => $evidence->{checks}{lsu_hierarchy}{matched_path},
				},
			};
		} elsif (!$args{excluded} || !$args{excluded}->($sample)) {
			push @incomplete, $sample->{name};
		}
	}
	return {cohort => \@cohort, incomplete => \@incomplete};
}

# The input of one marker's merge: a "<column>\t<hierarchy>" line per cohort
# sample, and a signature over the columns, paths, sizes and modification
# times. A merge is current only while its stored signature matches, so a new,
# removed or re-assigned sample leads to a new merge.
sub ribosome_merge_manifest {
	my (%args) = @_;
	die "ribosome_merge_manifest requires a cohort\n"
		unless ref($args{cohort}) eq 'ARRAY';
	die "ribosome_merge_manifest requires tag SSU or LSU\n"
		unless defined($args{tag}) && $args{tag} =~ /^(?:SSU|LSU)$/;
	my (@lines, @inputs);
	for my $sample (sort { $a->{name} cmp $b->{name} } @{$args{cohort}}) {
		my $path = $sample->{hierarchy}{$args{tag}};
		my @stat = defined($path) ? stat($path) : ();
		die "RiboFind sample $sample->{name} has no $args{tag} hierarchy\n"
			unless @stat;
		my $column = "$sample->{name}.$args{tag}";
		die "RiboFind merge input contains a tab or newline: $column $path\n"
			if "$column$path" =~ /[\t\r\n]/;
		push @lines, "$column\t$path\n";
		push @inputs, [$column, $path, $stat[7], $stat[9]];
	}
	return {
		text => join('', @lines),
		signature => completion_request_signature({
			ribofind_merge => 1, tag => $args{tag},
			levels => $args{levels} || [], inputs => \@inputs,
		}),
	};
}

# The LCA jobs read each reference (FASTA, its LAMBDA index and taxonomy) from
# a shared copy. That copy is current only when every file matches its source
# in size: existence alone accepts the truncated files an interrupted copy job
# leaves behind, and every later LCA job then fails on them.
sub lca_reference_copy_current {
	my (%args) = @_;
	for my $key (qw(fasta taxonomy directory)) {
		die "lca_reference_copy_current requires $key\n"
			unless defined($args{$key}) && length($args{$key});
	}
	for my $source ($args{fasta}, "$args{fasta}.lba.gz", $args{taxonomy}) {
		my @source = stat($source);
		return 0 unless @source;
		my @copy = stat(File::Spec->catfile($args{directory}, basename($source)));
		return 0 unless @copy && $copy[7] == $source[7];
	}
	return 1;
}

1;
