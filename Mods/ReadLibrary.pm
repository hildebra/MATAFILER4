package Mods::ReadLibrary;

use strict;
use warnings;
use Exporter qw(import);

our @EXPORT_OK = qw(
	newReadLibrary cloneReadLibraries validateReadLibraries
	readLibrariesFromArrays
	ensureSeqSetLibraries ensureCleanSeqSetLibraries
	syncSeqSetLegacy syncCleanSeqSetLegacy
	readLibraries readLibrariesByScope libraryFiles libraryPairs
	libraryTechnology legacyLibraryArrays replaceScopeLibraries
);

my @FILE_ROLES = qw(r1 r2 single bam);

sub _copy_hash {
	my ($hash) = @_;
	return { map { $_ => $hash->{$_} } keys %{$hash || {}} };
}

sub newReadLibrary {
	my (%args) = @_;
	my $scope = $args{scope} || 'primary';
	die "Read library scope must be 'primary' or 'support', not '$scope'\n"
		unless ($scope eq 'primary' || $scope eq 'support');
	my $files = _copy_hash($args{files});
	foreach my $role (@FILE_ROLES) {
		$files->{$role} = '' unless defined($files->{$role});
	}
	my $library = {
		id => defined($args{id}) ? $args{id} : '',
		sample => defined($args{sample}) ? $args{sample} : '',
		scope => $scope,
		technology => defined($args{technology}) ? $args{technology} : '',
		is_long => ($args{is_long} || (($args{technology} || '') =~ /^(?:PB|ONT)$/)) ? 1 : 0,
		label => defined($args{label}) ? $args{label} : '',
		phase => defined($args{phase}) ? $args{phase} : 'raw',
		files => $files,
	};
	$library->{source_files} = _copy_hash($args{source_files}) if ($args{source_files});
	$library->{metadata} = _copy_hash($args{metadata}) if ($args{metadata});
	validateReadLibraries([$library]);
	return $library;
}

sub cloneReadLibraries {
	my ($libraries) = @_;
	return [map {
		my %copy = %{$_};
		$copy{files} = _copy_hash($_->{files});
		$copy{source_files} = _copy_hash($_->{source_files}) if ($_->{source_files});
		$copy{metadata} = _copy_hash($_->{metadata}) if ($_->{metadata});
		\%copy;
	} @{$libraries || []}];
}

sub validateReadLibraries {
	my ($libraries) = @_;
	die "Read libraries must be an array reference\n" unless (ref($libraries) eq 'ARRAY');
	my %ids;
	foreach my $library (@{$libraries}) {
		die "Each read library must be a hash reference\n" unless (ref($library) eq 'HASH');
		die "Read library has invalid scope\n"
			unless (($library->{scope} || '') eq 'primary' || ($library->{scope} || '') eq 'support');
		die "Read library files must be a hash reference\n" unless (ref($library->{files}) eq 'HASH');
		my $r1 = $library->{files}{r1} || '';
		my $r2 = $library->{files}{r2} || '';
		die "Read library '$library->{id}' has only one mate ($r1 / $r2)\n" if (($r1 eq '') != ($r2 eq ''));
		die "Read library '$library->{id}' contains no reads\n"
			unless grep { defined($library->{files}{$_}) && $library->{files}{$_} ne '' } @FILE_ROLES;
		if (defined($library->{id}) && $library->{id} ne '') {
			die "Duplicate read library id '$library->{id}'\n" if ($ids{$library->{id}}++);
		}
	}
	return $libraries;
}

sub _value_at {
	my ($array, $index) = @_;
	return '' unless (ref($array) eq 'ARRAY' && defined($array->[$index]));
	return $array->[$index];
}

sub _legacy_scope_libraries {
	my ($object, $scope, $clean, $sample) = @_;
	my ($r1_key, $r2_key, $single_key, $label_key, $technology_key, $long_key);
	if ($clean) {
		($r1_key, $r2_key, $single_key, $label_key, $technology_key, $long_key) = $scope eq 'support'
			? qw(arpX1 arpX2 singArX matArX readTecX is3rdGenX)
			: qw(arp1 arp2 singAr matAr readTec is3rdGen);
	} else {
		($r1_key, $r2_key, $single_key, $label_key, $technology_key, $long_key) = $scope eq 'support'
			? qw(paX1 paX2 paXs libInfoX seqTechX is3rdGenX)
			: qw(pa1 pa2 pas libInfo seqTech is3rdGen);
	}
	my $r1 = ref($object->{$r1_key}) eq 'ARRAY' ? $object->{$r1_key} : [];
	my $r2 = ref($object->{$r2_key}) eq 'ARRAY' ? $object->{$r2_key} : [];
	my $single = ref($object->{$single_key}) eq 'ARRAY' ? $object->{$single_key} : [];
	my $labels = ref($object->{$label_key}) eq 'ARRAY' ? $object->{$label_key} : [];
	my $slots = @$r1;
	$slots = @$r2 if (@$r2 > $slots);
	$slots = @$single if (@$single > $slots);
	$slots = @$labels if (@$labels > $slots);
	my @libraries;
	for (my $i = 0; $i < $slots; $i++) {
		my %files = (
			r1 => _value_at($r1, $i), r2 => _value_at($r2, $i),
			single => _value_at($single, $i), bam => '',
		);
		next unless grep { $files{$_} ne '' } @FILE_ROLES;
		push @libraries, newReadLibrary(
			id => join(':', grep { $_ ne '' } ($sample || '', $scope, $i)),
			sample => $sample || '', scope => $scope,
			technology => $object->{$technology_key} || '',
			is_long => $object->{$long_key} || 0,
			label => _value_at($labels, $i) || "lib$i",
			phase => $clean ? 'clean' : 'raw', files => \%files,
		);
	}
	return \@libraries;
}

sub readLibrariesFromArrays {
	my (%args) = @_;
	my $scope = $args{scope} || 'primary';
	my $clean = $args{clean} ? 1 : 0;
	if ($args{separate_roles}) {
		my $r1 = $args{r1} || [];
		my $r2 = $args{r2} || [];
		my $single = $args{single} || [];
		my $labels = $args{labels} || [];
		my $pairTechnology = defined($args{pair_technology}) ? $args{pair_technology} : ($args{technology} || '');
		my $singleTechnology = defined($args{single_technology}) ? $args{single_technology} : ($args{technology} || '');
		my @libraries;
		my $recordIndex = 0;
		my $pairSlots = @$r1 > @$r2 ? @$r1 : @$r2;
		for (my $i = 0; $i < $pairSlots; $i++) {
			push @libraries, newReadLibrary(
				id => join(':', grep { $_ ne '' } ($args{sample} || '', $scope, $recordIndex++)),
				sample => $args{sample} || '', scope => $scope,
				technology => $pairTechnology,
				is_long => is_long_technology($pairTechnology, $args{is_long}),
				label => _value_at($labels, $i) || "lib$i",
				phase => $args{phase} || ($clean ? 'clean' : 'raw'),
				files => {r1 => _value_at($r1, $i), r2 => _value_at($r2, $i), single => '', bam => ''},
			);
		}
		for (my $i = 0; $i < @$single; $i++) {
			my $label = _value_at($labels, $i) || "lib$i";
			$label .= '.single' if $pairSlots;
			push @libraries, newReadLibrary(
				id => join(':', grep { $_ ne '' } ($args{sample} || '', $scope, $recordIndex++)),
				sample => $args{sample} || '', scope => $scope,
				technology => $singleTechnology,
				is_long => is_long_technology($singleTechnology, $args{is_long}),
				label => $label,
				phase => $args{phase} || ($clean ? 'clean' : 'raw'),
				files => {r1 => '', r2 => '', single => _value_at($single, $i), bam => ''},
			);
		}
		return \@libraries;
	}
	my %legacy;
	if ($clean) {
		if ($scope eq 'support') {
			@legacy{qw(arpX1 arpX2 singArX matArX readTecX is3rdGenX)} =
				($args{r1} || [], $args{r2} || [], $args{single} || [], $args{labels} || [], $args{technology} || '', $args{is_long} || 0);
		} else {
			@legacy{qw(arp1 arp2 singAr matAr readTec is3rdGen)} =
				($args{r1} || [], $args{r2} || [], $args{single} || [], $args{labels} || [], $args{technology} || '', $args{is_long} || 0);
		}
	} else {
		if ($scope eq 'support') {
			@legacy{qw(paX1 paX2 paXs libInfoX seqTechX is3rdGenX)} =
				($args{r1} || [], $args{r2} || [], $args{single} || [], $args{labels} || [], $args{technology} || '', $args{is_long} || 0);
		} else {
			@legacy{qw(pa1 pa2 pas libInfo seqTech is3rdGen)} =
				($args{r1} || [], $args{r2} || [], $args{single} || [], $args{labels} || [], $args{technology} || '', $args{is_long} || 0);
		}
	}
	my $libraries = _legacy_scope_libraries(\%legacy, $scope, $clean, $args{sample});
	foreach my $library (@{$libraries}) {
		$library->{phase} = $args{phase} if defined($args{phase});
	}
	return $libraries;
}

sub is_long_technology {
	my ($technology, $fallback) = @_;
	return 1 if defined($technology) && ($technology eq 'PB' || $technology eq 'ONT');
	return 0 if defined($technology) && $technology ne '';
	return $fallback ? 1 : 0;
}

sub ensureSeqSetLibraries {
	my ($seqset, $sample) = @_;
	die "seqSet must be a hash reference\n" unless (ref($seqset) eq 'HASH');
	if (ref($seqset->{libraries}) ne 'ARRAY') {
		$seqset->{libraries} = [
			@{_legacy_scope_libraries($seqset, 'primary', 0, $sample)},
			@{_legacy_scope_libraries($seqset, 'support', 0, $sample)},
		];
	}
	validateReadLibraries($seqset->{libraries}) if (@{$seqset->{libraries}});
	foreach my $library (@{$seqset->{libraries}}) {
		die "Raw seqSet contains non-raw library '$library->{id}' with phase '$library->{phase}'\n"
			unless (($library->{phase} || '') eq 'raw' || ($library->{phase} || '') eq 'staged');
		die "Raw seqSet for sample '$sample' contains library '$library->{id}' owned by '$library->{sample}'\n"
			if (defined($sample) && $sample ne '' && ($library->{sample} || '') ne ''
				&& $library->{sample} ne $sample);
	}
	return $seqset->{libraries};
}

sub ensureCleanSeqSetLibraries {
	my ($clean, $sample) = @_;
	die "cleanSeqSet must be a hash reference\n" unless (ref($clean) eq 'HASH');
	if (ref($clean->{libraries}) ne 'ARRAY') {
		$clean->{libraries} = [
			@{_legacy_scope_libraries($clean, 'primary', 1, $sample)},
			@{_legacy_scope_libraries($clean, 'support', 1, $sample)},
		];
	}
	if (ref($clean->{merged_library}) ne 'HASH' && ref($clean->{mrgHshHR}) eq 'HASH'
		&& ($clean->{mrgHshHR}{mrg} || '') ne '') {
		$clean->{merged_library} = newReadLibrary(
			id => join(':', grep { $_ ne '' } ($sample || '', 'primary', 'merged')),
			sample => $sample || '', scope => 'primary', technology => '', label => 'merged', phase => 'merged',
			files => {r1 => $clean->{mrgHshHR}{pair1} || '', r2 => $clean->{mrgHshHR}{pair2} || '',
				single => $clean->{mrgHshHR}{mrg} || '', bam => ''},
		);
	}
	validateReadLibraries($clean->{libraries}) if (@{$clean->{libraries}});
	foreach my $library (@{$clean->{libraries}}) {
		die "cleanSeqSet contains non-clean library '$library->{id}' with phase '$library->{phase}'\n"
			unless (($library->{phase} || '') eq 'clean');
		die "cleanSeqSet for sample '$sample' contains library '$library->{id}' owned by '$library->{sample}'\n"
			if (defined($sample) && $sample ne '' && ($library->{sample} || '') ne ''
				&& $library->{sample} ne $sample);
	}
	return $clean->{libraries};
}

sub readLibraries {
	my ($object, $clean, $sample) = @_;
	return $clean ? ensureCleanSeqSetLibraries($object, $sample) : ensureSeqSetLibraries($object, $sample);
}

sub readLibrariesByScope {
	my ($object, $scope, $clean, $sample) = @_;
	return [grep { ($_->{scope} || '') eq $scope } @{readLibraries($object, $clean, $sample)}];
}

sub libraryFiles {
	my ($libraries, $role) = @_;
	die "Unknown read-library role '$role'\n" unless grep { $_ eq $role } @FILE_ROLES;
	return [map { $_->{files}{$role} } grep { defined($_->{files}{$role}) && $_->{files}{$role} ne '' } @{$libraries || []}];
}

sub libraryPairs {
	my ($libraries) = @_;
	return [grep { ($_->{files}{r1} || '') ne '' && ($_->{files}{r2} || '') ne '' } @{$libraries || []}];
}

sub libraryTechnology {
	my ($libraries, $context, $requireAll) = @_;
	if ($requireAll) {
		my @missing = grep { !defined($_->{technology}) || $_->{technology} eq '' } @{$libraries || []};
		die "Read libraries for $context have no technology: ".join(', ', map { $_->{id} || '<unnamed>' } @missing)."\n"
			if @missing;
	}
	my %technologies = map { ($_->{technology} || '') => 1 } @{$libraries || []};
	delete $technologies{''};
	if (keys(%technologies) > 1) {
		my $where = defined($context) && $context ne '' ? " for $context" : '';
		die "Read libraries$where mix technologies: ".join(', ', sort keys %technologies)."\n";
	}
	my ($technology) = keys %technologies;
	return $technology || '';
}

sub legacyLibraryArrays {
	my ($libraries, $aligned) = @_;
	my (@r1, @r2, @single, @labels, @technologies);
	foreach my $library (@{$libraries || []}) {
		if ($aligned) {
			push @r1, $library->{files}{r1} || '';
			push @r2, $library->{files}{r2} || '';
			push @single, $library->{files}{single} || '';
		} else {
			push @r1, $library->{files}{r1} if ($library->{files}{r1});
			push @r2, $library->{files}{r2} if ($library->{files}{r2});
			push @single, $library->{files}{single} if ($library->{files}{single});
		}
		push @labels, $library->{label} || '';
		push @technologies, $library->{technology} || '';
	}
	return (\@r1, \@r2, \@single, \@labels, \@technologies);
}

sub _scope_metadata {
	my ($libraries) = @_;
	my ($technology) = map { $_->{technology} }
		grep { defined($_->{technology}) && $_->{technology} ne '' } @{$libraries || []};
	my $is_long = grep { $_->{is_long} } @{$libraries || []};
	return ($technology || '', $is_long ? 1 : 0);
}

sub syncSeqSetLegacy {
	my ($seqset) = @_;
	my $libraries = ensureSeqSetLibraries($seqset);
	foreach my $scope (qw(primary support)) {
		my $selected = [grep { $_->{scope} eq $scope } @{$libraries}];
		my ($r1, $r2, $single, $labels) = legacyLibraryArrays($selected, 0);
		my ($technology, $is_long) = _scope_metadata($selected);
		if ($scope eq 'primary') {
			@$seqset{qw(pa1 pa2 pas libInfo seqTech is3rdGen)} = ($r1, $r2, $single, $labels, $technology, $is_long);
		} else {
			@$seqset{qw(paX1 paX2 paXs libInfoX seqTechX is3rdGenX)} = ($r1, $r2, $single, $labels, $technology, $is_long);
		}
	}
	return $seqset;
}

sub syncCleanSeqSetLegacy {
	my ($clean) = @_;
	my $libraries = ensureCleanSeqSetLibraries($clean);
	foreach my $scope (qw(primary support)) {
		my $selected = [grep { $_->{scope} eq $scope } @{$libraries}];
		my ($r1, $r2, $single, $labels) = legacyLibraryArrays($selected, 0);
		my ($technology, $is_long) = _scope_metadata($selected);
		if ($scope eq 'primary') {
			@$clean{qw(arp1 arp2 singAr matAr readTec is3rdGen)} = ($r1, $r2, $single, $labels, $technology, $is_long);
		} else {
			@$clean{qw(arpX1 arpX2 singArX matArX readTecX is3rdGenX)} = ($r1, $r2, $single, $labels, $technology, $is_long);
		}
	}
	$clean->{mrgHshHR} ||= {};
	if (ref($clean->{merged_library}) eq 'HASH') {
		$clean->{mrgHshHR} = {
			mrg => $clean->{merged_library}{files}{single} || '',
			pair1 => $clean->{merged_library}{files}{r1} || '',
			pair2 => $clean->{merged_library}{files}{r2} || '',
		};
	}
	return $clean;
}

sub replaceScopeLibraries {
	my ($object, $scope, $replacement, $clean, $sample) = @_;
	validateReadLibraries($replacement) if (@{$replacement || []});
	my $libraries = readLibraries($object, $clean, $sample);
	$object->{libraries} = [(grep { $_->{scope} ne $scope } @{$libraries}), @{$replacement || []}];
	return $clean ? syncCleanSeqSetLegacy($object) : syncSeqSetLegacy($object);
}

1;
