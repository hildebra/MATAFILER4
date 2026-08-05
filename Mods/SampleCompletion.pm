package Mods::SampleCompletion;
use strict;
use warnings;

use Digest::SHA qw(sha256_hex);
use Exporter qw(import);
use File::Basename qw(dirname);
use File::Path qw(make_path);
use File::Spec;
use JSON::PP ();

our @EXPORT_OK = qw(
	sample_completion_path completion_request_signature
	completion_component_evidence
	read_sample_completion write_sample_completion
	invalidate_sample_completion
);

my $SENTINEL_NAME = 'MATAFILER.sample.complete.json';
my $SCHEMA_NAME = 'MATAFILER.sample-completion';
my $SCHEMA_VERSION = 3;
my $COMPONENT_SCHEMA_NAME = 'MATAFILER.sample-components';
my $COMPONENT_SCHEMA_VERSION = 1;
my $STATS_SCHEMA_NAME = 'MATAFILER.metagstats';
my $STATS_SCHEMA_VERSION = 1;

sub _sample_root {
	my ($root) = @_;
	die "sample completion root is required\n"
		unless defined($root) && length($root);
	my $canonical = File::Spec->canonpath(File::Spec->rel2abs($root));
	my $volume_root = File::Spec->rootdir();
	die "refusing unsafe sample completion root '$root'\n"
		if $canonical eq $volume_root;
	return $canonical;
}

sub sample_completion_path {
	my ($root) = @_;
	return File::Spec->catfile(_sample_root($root), $SENTINEL_NAME);
}

sub completion_request_signature {
	my ($request) = @_;
	die "completion request must be a hash reference\n"
		unless ref($request) eq 'HASH';
	my $json = JSON::PP->new->canonical(1)->utf8(1)->encode($request);
	return sha256_hex($json);
}

sub _json_codec {
	return JSON::PP->new->canonical(1)->pretty(1)->utf8(1);
}

sub _check_ok {
	my ($check) = @_;
	my $kind = $check->{kind} || '';
	if ($kind eq 'exists') {
		return -e $check->{path} ? 1 : 0;
	}
	if ($kind eq 'nonempty') {
		return -s $check->{path} ? 1 : 0;
	}
	if ($kind eq 'exists_any') {
		return scalar(grep { -e $_ } @{$check->{paths}}) ? 1 : 0;
	}
	if ($kind eq 'nonempty_any') {
		return scalar(grep { -s $_ } @{$check->{paths}}) ? 1 : 0;
	}
	if ($kind eq 'directory') {
		return -d $check->{path} ? 1 : 0;
	}
	die "unknown completion check kind '$kind'\n";
}

sub completion_component_evidence {
	my (%args) = @_;
	my $requested = $args{requested} ? 1 : 0;
	my $applicable = exists($args{applicable}) ? ($args{applicable} ? 1 : 0) : 1;
	my $checks = $args{checks} || [];
	die "completion component checks must be an array reference\n"
		unless ref($checks) eq 'ARRAY';

	my %evidence_checks;
	if ($requested && $applicable) {
		for my $check (@{$checks}) {
			die "completion component check must be a hash reference\n"
				unless ref($check) eq 'HASH';
			my $id = $check->{id} || '';
			die "completion component check id is required\n" if $id eq '';
			die "duplicate completion component check '$id'\n"
				if exists $evidence_checks{$id};
			my $kind = $check->{kind} || '';
			die "completion component check '$id' has no path\n"
				if $kind =~ /^(?:exists|nonempty|directory)$/
					&& (!defined($check->{path}) || $check->{path} eq '');
			die "completion component check '$id' has no candidate paths\n"
				if $kind =~ /^(?:exists_any|nonempty_any)$/
					&& (ref($check->{paths}) ne 'ARRAY' || !@{$check->{paths}});
			my %record = (
				kind => $kind,
				required => exists($check->{required}) ? ($check->{required} ? 1 : 0) : 1,
			);
			$record{path} = $check->{path} if defined($check->{path});
			$record{paths} = [@{$check->{paths}}] if ref($check->{paths}) eq 'ARRAY';
			$record{ok} = _check_ok({%{$check}, kind => $kind});
			$evidence_checks{$id} = \%record;
		}
	}

	my $complete = 1;
	if ($requested && $applicable) {
		$complete = !grep {
			$evidence_checks{$_}{required} && !$evidence_checks{$_}{ok}
		} keys %evidence_checks;
		$complete = $complete ? 1 : 0;
	}
	my $status = !$requested ? 'not_requested'
		: !$applicable ? 'not_applicable'
		: $complete ? 'complete' : 'incomplete';
	return {
		requested => $requested,
		applicable => $applicable,
		complete => $complete,
		status => $status,
		reason => $args{reason} || '',
		checks => \%evidence_checks,
	};
}

sub _contracts_error {
	my ($record) = @_;
	return 'sentinel contract registry is missing'
		unless ref($record->{contracts}) eq 'HASH';
	my $components = $record->{contracts}{components};
	return 'sentinel component contract is unsupported'
		unless ref($components) eq 'HASH'
			&& ($components->{schema} || '') eq $COMPONENT_SCHEMA_NAME
			&& ($components->{schema_version} || 0) == $COMPONENT_SCHEMA_VERSION;
	my $stats = $record->{contracts}{metagstats};
	return 'sentinel metagStats contract is unsupported'
		unless ref($stats) eq 'HASH'
			&& ($stats->{schema} || '') eq $STATS_SCHEMA_NAME
			&& ($stats->{schema_version} || 0) == $STATS_SCHEMA_VERSION;
	my @inventories = (
		['component', $components->{inventory},
			ref($record->{components}) eq 'HASH' ? [keys %{$record->{components}}] : []],
		['metagStats field', $stats->{fields},
			ref($record->{metagstats}{field_availability}) eq 'HASH'
				? [keys %{$record->{metagstats}{field_availability}}] : []],
		['metagStats family', $stats->{families},
			ref($record->{metagstats}{families}) eq 'HASH'
				? [keys %{$record->{metagstats}{families}}] : []],
	);
	for my $inventory (@inventories) {
		my ($label, $stored, $actual) = @{$inventory};
		return "sentinel $label contract inventory is missing"
			unless ref($stored) eq 'ARRAY';
		my $stored_json = JSON::PP->new->canonical(1)->encode([sort @{$stored}]);
		my $actual_json = JSON::PP->new->canonical(1)->encode([sort @{$actual}]);
		return "sentinel $label contract inventory does not match its evidence"
			if $stored_json ne $actual_json;
	}
	return '';
}

sub _component_revalidation_projection {
	my ($components) = @_;
	my %projection;
	for my $name (sort keys %{$components || {}}) {
		my $component = $components->{$name};
		next unless ref($component) eq 'HASH';
		my %checks;
		for my $check_name (sort keys %{$component->{checks} || {}}) {
			my $check = $component->{checks}{$check_name};
			next unless ref($check) eq 'HASH' && $check->{required};
			$checks{$check_name} = {%{$check}};
		}
		$projection{$name} = {
			requested => $component->{requested} ? 1 : 0,
			applicable => $component->{applicable} ? 1 : 0,
			complete => $component->{complete} ? 1 : 0,
			status => $component->{status} || '',
			checks => \%checks,
		};
	}
	return \%projection;
}

sub _metagstats_error {
	my ($stats) = @_;
	return 'sentinel metagStats record is incomplete'
		unless ref($stats) eq 'HASH' && defined($stats->{DIR})
			&& ref($stats->{values}) eq 'HASH'
			&& ref($stats->{field_availability}) eq 'HASH'
			&& ref($stats->{families}) eq 'HASH';
	for my $field (sort keys %{$stats->{field_availability}}) {
		my $available = $stats->{field_availability}{$field};
		return "sentinel metagStats availability for '$field' is malformed"
			unless defined($available) && $available =~ /^(?:0|1)$/;
		my $value = $stats->{values}{$field};
		my $derived = exists($stats->{values}{$field})
			&& defined($value) && $value ne '' && $value ne '-1' ? 1 : 0;
		return "sentinel metagStats availability for '$field' is inconsistent"
			if $available != $derived;
	}
	my %allowed_status = map { $_ => 1 } qw(
		not_applicable missing partial complete
	);
	for my $name (sort keys %{$stats->{families}}) {
		my $family = $stats->{families}{$name};
		return "sentinel metagStats family '$name' is malformed"
			unless ref($family) eq 'HASH'
				&& defined($family->{requested}) && $family->{requested} =~ /^(?:0|1)$/
				&& defined($family->{applicable}) && $family->{applicable} =~ /^(?:0|1)$/
				&& defined($family->{ok}) && $family->{ok} =~ /^(?:0|1)$/
				&& defined($family->{status})
				&& $allowed_status{$family->{status}}
				&& ref($family->{fields}) eq 'HASH';
		return "sentinel metagStats family '$name' has inconsistent applicability"
			if ($family->{requested} != $family->{applicable});
		my $expected_ok = $family->{status} =~ /^(?:not_applicable|complete)$/ ? 1 : 0;
		return "sentinel metagStats family '$name' has inconsistent status"
			if ($family->{ok} != $expected_ok);
		for my $field (sort keys %{$family->{fields}}) {
			my $available = $family->{fields}{$field};
			return "sentinel metagStats family '$name' field '$field' is malformed"
				unless defined($available) && $available =~ /^(?:0|1)$/;
			return "sentinel metagStats family '$name' field '$field' is not declared"
				unless exists($stats->{field_availability}{$field});
			return "sentinel metagStats family '$name' field '$field' is inconsistent"
				if $available != $stats->{field_availability}{$field};
		}
	}
	return '';
}

sub _completion_record_error {
	my ($record, $expected_components) = @_;
	return 'sentinel outcome record is missing'
		unless ref($record->{outcome}) eq 'HASH'
			&& defined($record->{outcome}{status})
			&& $record->{outcome}{status} ne '';
	my %allowed_status = map { $_ => 1 } qw(
		completed skipped_too_small skipped_empty_input skipped_sdm_warning
	);
	return "unknown sentinel outcome '$record->{outcome}{status}'"
		unless $allowed_status{$record->{outcome}{status}};
	return 'sentinel component evidence is missing'
		unless ref($record->{components}) eq 'HASH';
	for my $name (sort keys %{$record->{components}}) {
		my $component = $record->{components}{$name};
		return "sentinel component '$name' is malformed"
			unless ref($component) eq 'HASH'
				&& defined($component->{requested}) && $component->{requested} =~ /^(?:0|1)$/
				&& defined($component->{applicable}) && $component->{applicable} =~ /^(?:0|1)$/
				&& defined($component->{complete}) && $component->{complete} =~ /^(?:0|1)$/
				&& defined($component->{status})
				&& ref($component->{checks}) eq 'HASH';
		return "sentinel component '$name' has no completion checks"
			if $component->{requested} && $component->{applicable}
				&& !keys %{$component->{checks}};
		my $expected_status = !$component->{requested} ? 'not_requested'
			: !$component->{applicable} ? 'not_applicable'
			: $component->{complete} ? 'complete' : 'incomplete';
		return "sentinel component '$name' has inconsistent status"
			unless $component->{status} eq $expected_status;
		my $failed_required = 0;
		for my $check_name (sort keys %{$component->{checks}}) {
			my $check = $component->{checks}{$check_name};
			return "sentinel component '$name' check '$check_name' is malformed"
				unless ref($check) eq 'HASH'
					&& defined($check->{kind})
					&& $check->{kind} =~ /^(?:exists|nonempty|exists_any|nonempty_any|directory)$/
					&& defined($check->{required}) && $check->{required} =~ /^(?:0|1)$/
					&& defined($check->{ok}) && $check->{ok} =~ /^(?:0|1)$/;
			$failed_required = 1 if $check->{required} && !$check->{ok};
		}
		my $derived_complete = $component->{requested} && $component->{applicable}
			? !$failed_required : 1;
		return "sentinel component '$name' has inconsistent completeness"
			if ($component->{complete} != $derived_complete);
		if ($record->{outcome}{status} eq 'completed'
				&& $component->{requested} && $component->{applicable}
				&& !$component->{complete}) {
			return "sentinel says requested component '$name' is incomplete";
		}
	}
	return '' unless defined($expected_components);
	return 'expected component evidence is malformed'
		unless ref($expected_components) eq 'HASH';
	return '' if $record->{outcome}{status} ne 'completed';
	my $stored = JSON::PP->new->canonical(1)->encode(
		_component_revalidation_projection($record->{components}),
	);
	my $current = JSON::PP->new->canonical(1)->encode(
		_component_revalidation_projection($expected_components),
	);
	return 'required workflow component outputs have changed or are missing'
		if $stored ne $current;
	return '';
}

sub read_sample_completion {
	my (%args) = @_;
	my $path = defined($args{path}) && length($args{path})
		? $args{path} : sample_completion_path($args{root});
	return wantarray ? (undef, 'missing') : undef unless -e $path;
	return wantarray ? (undef, 'empty sentinel') : undef unless -s $path;

	open my $fh, '<:raw', $path
		or return wantarray ? (undef, "cannot read sentinel: $!") : undef;
	local $/;
	my $json = <$fh>;
	close $fh
		or return wantarray ? (undef, "cannot close sentinel: $!") : undef;
	my $record = eval { _json_codec()->decode($json) };
	return wantarray ? (undef, "invalid JSON: $@") : undef
		unless ref($record) eq 'HASH';

	my $error = '';
	if (($record->{schema} || '') ne $SCHEMA_NAME) {
		$error = 'unknown sentinel schema';
	} elsif (($record->{schema_version} || 0) != $SCHEMA_VERSION) {
		$error = 'unsupported sentinel schema version';
	} else {
		$error = _contracts_error($record);
	}
	if ($error eq '' && (!defined($record->{sample}) || $record->{sample} eq '')) {
		$error = 'sentinel sample is missing';
	} elsif ($error eq '' && defined($args{sample}) && $record->{sample} ne $args{sample}) {
		$error = "sentinel belongs to sample '$record->{sample}'";
	} elsif ($error eq '' && (!defined($record->{request_signature})
			|| $record->{request_signature} !~ /^[0-9a-f]{64}$/)) {
		$error = 'sentinel request signature is invalid';
	} elsif ($error eq '' && defined($args{request_signature})
			&& $record->{request_signature} ne $args{request_signature}) {
		$error = 'requested workflow differs from the closed workflow';
	} elsif ($error eq '') {
		$error = _metagstats_error($record->{metagstats});
		$error = _completion_record_error($record, $args{expected_components})
			if $error eq '';
	}
	return wantarray ? (undef, $error) : undef if $error ne '';
	return wantarray ? ($record, '') : $record;
}

sub write_sample_completion {
	my (%args) = @_;
	for my $required (qw(root sample request_signature metagstats components outcome)) {
		die "sample completion $required is required\n"
			unless exists($args{$required}) && defined($args{$required});
	}
	die "sample completion request signature is invalid\n"
		unless $args{request_signature} =~ /^[0-9a-f]{64}$/;
	my $stats_validation_error = _metagstats_error($args{metagstats});
	die "invalid sample completion metagStats evidence: $stats_validation_error\n"
		if $stats_validation_error ne '';
	die "sample completion components must be a hash reference\n"
		unless ref($args{components}) eq 'HASH';
	die "sample completion outcome must contain a status\n"
		unless ref($args{outcome}) eq 'HASH'
			&& defined($args{outcome}{status}) && $args{outcome}{status} ne '';
	my $validation_error = _completion_record_error({
		components => $args{components}, outcome => $args{outcome},
	});
	die "invalid sample completion evidence: $validation_error\n"
		if $validation_error ne '';

	my $path = sample_completion_path($args{root});
	my $directory = dirname($path);
	make_path($directory) unless -d $directory;
	my $record = {
		schema => $SCHEMA_NAME,
		schema_version => $SCHEMA_VERSION,
		contracts => {
			components => {
				schema => $COMPONENT_SCHEMA_NAME,
				schema_version => $COMPONENT_SCHEMA_VERSION,
				inventory => [sort keys %{$args{components}}],
			},
			metagstats => {
				schema => $STATS_SCHEMA_NAME,
				schema_version => $STATS_SCHEMA_VERSION,
				fields => [sort keys %{$args{metagstats}{field_availability}}],
				families => [sort keys %{$args{metagstats}{families}}],
			},
		},
		sample => $args{sample},
		request_signature => $args{request_signature},
		created_epoch => time,
		present_assembly => $args{present_assembly} ? 1 : 0,
		empty_sample => $args{empty_sample} ? 1 : 0,
		empty_input_size_mb => 0 + ($args{empty_input_size_mb} || 0),
		components => $args{components},
		outcome => $args{outcome},
		metagstats => $args{metagstats},
	};
	my $temporary = "$path.tmp.$$";
	my $json = _json_codec()->encode($record);
	open my $fh, '>:raw', $temporary
		or die "cannot write sample completion sentinel $temporary: $!\n";
	print {$fh} $json
		or die "cannot populate sample completion sentinel $temporary: $!\n";
	close $fh
		or die "cannot close sample completion sentinel $temporary: $!\n";
	rename $temporary, $path
		or die "cannot publish sample completion sentinel $path: $!\n";
	return $path;
}

sub invalidate_sample_completion {
	my ($root) = @_;
	my $path = sample_completion_path($root);
	return 0 unless -e $path;
	unlink $path or die "cannot invalidate sample completion sentinel $path: $!\n";
	return 1;
}

1;
