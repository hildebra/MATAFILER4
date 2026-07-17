package Mods::WorkflowControl;

use warnings;
use strict;

use Exporter qw(import);

our @EXPORT_OK = qw(
	advance_loop_window
	assembly_group_output_dirs
	hybrid_group_ready
	hybrid_package_complete
	hybrid_package_sample_id
	hybrid_local_scratch_gb
	missing_input_files
	normalise_job_dependencies
	append_job_dependencies
	augment_deferred_submission
	source_input_files
	parse_ignored_samples
	sample_base_output_dir
	sample_is_ignored
	workflow_members_match
);

sub hybrid_package_sample_id {
	my ($package_dir) = @_;
	return '' unless defined($package_dir) && length($package_dir);
	$package_dir =~ s{/$}{};
	open my $manifest, '<', "$package_dir/package.manifest.tsv" or return '';
	while (my $line = <$manifest>) {
		$line =~ s/[\r\n]+$//;
		if ($line =~ /^sample_id\t(.+)$/) {
			my $sample_id = $1;
			close $manifest;
			return $sample_id;
		}
	}
	close $manifest;
	return '';
}

sub hybrid_local_scratch_gb {
	my (%args) = @_;
	my $assembler_gb = 0 + ($args{assembler_gb} || 0);
	my $preassembly_bytes = 0 + ($args{preassembly_bytes} || 0);
	my $max_depth = 0 + ($args{max_synthetic_depth} || 0);
	die 'hybrid_local_scratch_gb requires non-negative inputs'
		if ($assembler_gb < 0 || $preassembly_bytes < 0 || $max_depth < 0);
	# Synthetic FASTQ is written compressed. Budget one byte per simulated base
	# at the configured depth cap, rounded up, in addition to assembler scratch.
	my $synthetic_gb = $preassembly_bytes * $max_depth / (1024 ** 3);
	$synthetic_gb = int($synthetic_gb) == $synthetic_gb
		? int($synthetic_gb) : int($synthetic_gb) + 1;
	return int($assembler_gb + $synthetic_gb);
}

sub normalise_job_dependencies {
	my (@values) = @_;
	my @pending = @values;
	my @dependencies;
	my %seen;
	while (@pending) {
		my $value = shift @pending;
		if (ref($value) eq 'ARRAY') {
			unshift @pending, @{$value};
			next;
		}
		next unless (defined $value);
		for my $dependency (split /[;,\s]+/, $value) {
			next if ($dependency eq '' || $seen{$dependency}++);
			push @dependencies, $dependency;
		}
	}
	return join(';', @dependencies);
}

sub append_job_dependencies {
	my ($target, @values) = @_;
	die 'append_job_dependencies requires a scalar reference'
		unless (ref($target) eq 'SCALAR');
	${$target} = normalise_job_dependencies(${$target}, @values);
	return ${$target};
}

sub augment_deferred_submission {
	my (%args) = @_;
	my $mode = $args{qmode} || '';
	my $command = $args{command} || '';
	my $script = defined($args{script}) ? $args{script} : '';
	my $rtag = $args{run_tag} || '';
	my $dependencies = normalise_job_dependencies($args{dependencies});
	return { command => $command, script => $script }
		if ($dependencies eq '' || $mode eq 'bash');
	my @dependencies = split /;/, $dependencies;
	for (@dependencies) { s/^\Q$rtag\E// if ($rtag ne ''); }

	if ($mode eq 'slurm') {
		my $new_dependencies = join(':', @dependencies);
		if ($script =~ /^#SBATCH --dependency=(afterok|afterany):([^\r\n]*)/m) {
			my ($policy, $existing) = ($1, $2);
			my $merged = normalise_job_dependencies(split(/:/, $existing), @dependencies);
			$merged =~ s/;/:/g;
			$script =~ s/^#SBATCH --dependency=(?:afterok|afterany):[^\r\n]*/#SBATCH --dependency=$policy:$merged/m;
		} else {
			$script =~ s/^(#![^\r\n]*\r?\n)/$1#SBATCH --dependency=afterok:$new_dependencies\n/;
		}
	} elsif ($mode eq 'sge') {
		my $new_dependencies = join(',', @dependencies);
		if ($command =~ /-hold_jid\s+(\S+)/) {
			my $existing = $1;
			my $merged = normalise_job_dependencies(split(/,/, $existing), @dependencies);
			$merged =~ s/;/,/g;
			$command =~ s/-hold_jid\s+\S+/-hold_jid $merged/;
		} else {
			$command =~ s/^(qsub\s+)/$1-hold_jid $new_dependencies /;
		}
	} elsif ($mode eq 'lsf') {
		my @existing;
		if ($command =~ /-w\s+"([^"]*)"/) {
			@existing = ($1 =~ /done\(([^)]+)\)/g);
		}
		my $merged = normalise_job_dependencies(@existing, @dependencies);
		my $expression = join(' && ', map { "done($_)" } split(/;/, $merged));
		if ($command =~ /-w\s+"[^"]*"/) {
			$command =~ s/-w\s+"[^"]*"/-w "$expression"/;
		} else {
			$command =~ s/^(bsub\s+)/$1-w "$expression" /;
		}
	} else {
		die "Cannot augment deferred submission for unknown queue mode '$mode'";
	}
	return { command => $command, script => $script };
}

sub _normalise_path {
	my ($path) = @_;
	$path = '' unless (defined $path);
	$path =~ s{\\}{/}g;
	$path =~ s{/+}{/}g;
	$path =~ s{/$}{};
	return $path;
}

sub sample_base_output_dir {
	my ($sample_output, $sample_key) = @_;
	die 'sample_base_output_dir requires a sample output path'
		unless (defined $sample_output && $sample_output ne '');
	die 'sample_base_output_dir requires a sample key'
		unless (defined $sample_key && $sample_key ne '');
	my $path = $sample_output;
	$path =~ s{[\\/]+$}{};
	die "Sample output '$sample_output' does not end in sample directory '$sample_key'"
		unless ($path =~ s{[\\/]\Q$sample_key\E$}{});
	return "$path/";
}

sub parse_ignored_samples {
	my ($value) = @_;
	$value = '' unless (defined $value);
	my %ignored = map { $_ => 1 }
		grep { $_ ne '' }
		map { my $name = $_; $name =~ s/^\s+|\s+$//g; $name }
		split /,/, $value;
	return \%ignored;
}

sub sample_is_ignored {
	my ($ignored, $sample_id) = @_;
	return 0 unless (ref($ignored) eq 'HASH' && defined $sample_id);
	return exists($ignored->{$sample_id}) ? 1 : 0;
}

sub workflow_members_match {
	my ($expected, $actual) = @_;
	return 0 unless (ref($expected) eq 'ARRAY' && ref($actual) eq 'ARRAY');
	my @expected = sort map { _normalise_path($_) } @{$expected};
	my @actual = sort map { _normalise_path($_) } @{$actual};
	return 0 unless (@expected == @actual);
	for (my $i = 0; $i < @expected; $i++) {
		return 0 if ($expected[$i] ne $actual[$i]);
	}
	return 1;
}

sub assembly_group_output_dirs {
	my ($map, $group_id) = @_;
	die 'assembly_group_output_dirs requires a map' unless (ref($map) eq 'HASH');
	my @dirs;
	for my $sample_key (@{$map->{opt}{smpl_order} || []}) {
		my $sample_group = $map->{$sample_key}{AssGroup};
		$sample_group = $sample_key
			if (!defined($sample_group) || $sample_group eq '-1');
		push @dirs, $map->{$sample_key}{wrdir} if ("$sample_group" eq "$group_id");
	}
	return [@dirs];
}

sub _nonempty_with_optional_gzip {
	my ($path) = @_;
	return 1 if (-s $path);
	return 1 if ($path =~ m{\.gz$} && -s substr($path, 0, -3));
	return 1 if ($path !~ m{\.gz$} && -s "$path.gz");
	return 0;
}

sub hybrid_package_complete {
	my ($package_dir) = @_;
	return 0 unless (defined $package_dir && $package_dir ne '');
	$package_dir =~ s{/$}{};
	return 0 unless (-s "$package_dir/scaffolds.fasta.filt");
	return 0 unless (_nonempty_with_optional_gzip("$package_dir/Coverage.percontig"));
	return 0 unless (_nonempty_with_optional_gzip("$package_dir/Coverage.median.percontig"));
	return 0 unless (-s "$package_dir/mapping.coverage.gz");
	return 0 unless (-s "$package_dir/breakpoints.tsv.gz");
	return 0 unless (-s "$package_dir/package.manifest.tsv");
	open my $manifest, '<', "$package_dir/package.manifest.tsv" or return 0;
	local $/;
	my $manifest_text = <$manifest>;
	close $manifest;
	return 0 unless defined($manifest_text) && $manifest_text =~ /^schema_version\t2$/m;
	return -e "$package_dir/moved.sto" ? 1 : 0;
}

sub hybrid_group_ready {
	my ($complete_packages, $members_without_primary, $target_members) = @_;
	return (0 + $complete_packages) + (0 + $members_without_primary)
		>= (0 + $target_members) ? 1 : 0;
}

sub missing_input_files {
	my (@files) = @_;
	return [grep { !defined($_) || $_ eq '' || !-e $_ || !-s $_ } @files];
}

sub source_input_files {
	my ($base_dir, @files) = @_;
	$base_dir = '' unless (defined $base_dir);
	$base_dir =~ s{/$}{};
	return [map {
		my $path = defined($_) ? $_ : '';
		($path eq '' || $path =~ m{^/} || $base_dir eq '') ? $path : "$base_dir/$path";
	} @files];
}

sub advance_loop_window {
	my (%args) = @_;
	my $from = 0 + $args{from};
	my $to = 0 + $args{to};
	my $upper = 0 + $args{upper};
	my $window = 0 + $args{window_size};
	my $initial_loops = 0 + $args{initial_loops};
	die 'advance_loop_window requires a positive window size' unless ($window > 0);

	my $next_from = $to;
	$next_from = $upper if ($next_from > $upper);
	my $next_to = $next_from + $window;
	$next_to = $upper if ($next_to > $upper);
	my $has_window = $next_from < $upper ? 1 : 0;
	return {
		from => $next_from,
		to => $next_to,
		loop_count => $has_window ? $initial_loops : 0,
		reset_index => $next_from - 1,
		has_window => $has_window,
	};
}

1;
