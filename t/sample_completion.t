use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use JSON::PP ();
use Test::More;

use lib File::Spec->catdir(File::Spec->rel2abs('.'));
use Mods::SampleCompletion qw(
	sample_completion_path completion_request_signature
	completion_component_evidence completion_record_needs_refresh
	read_sample_completion write_sample_completion invalidate_sample_completion
);

sub write_file {
	my ($path, $content) = @_;
	my (undef, $directory) = File::Spec->splitpath($path);
	make_path($directory) if $directory ne '' && !-d $directory;
	open my $fh, '>', $path or die "Cannot write $path: $!";
	print {$fh} $content;
	close $fh or die "Cannot close $path: $!";
}

sub read_json {
	my ($path) = @_;
	open my $fh, '<:raw', $path or die "Cannot read $path: $!";
	local $/;
	my $record = JSON::PP->new->decode(<$fh>);
	close $fh or die "Cannot close $path: $!";
	return $record;
}

sub stats_record {
	return {
		DIR => '/reads/S1',
		values => {RawInputSize => '1.250G', ScaffN50 => 12345},
		families => {input => 'complete', assembly => 'complete'},
	};
}

my $root = tempdir(CLEANUP => 1);
my $sample_root = File::Spec->catdir($root, 'S1');
my $request = {
	contract => 3,
	sample => {id => 'S1', assembly_group => ['S1']},
	requested => {assembly => 2, mapping => 1, binner => 0},
};
my $signature = completion_request_signature($request);
like($signature, qr/^[0-9a-f]{64}$/,
	'workflow request has a stable SHA-256 signature');
is(completion_request_signature({%{$request}}), $signature,
	'canonical request signatures do not depend on hash iteration order');
isnt(completion_request_signature({
	%{$request}, requested => {%{$request->{requested}}, binner => 2},
}), $signature, 'changing requested outputs changes the completion signature');

my $ribo_result = File::Spec->catfile($sample_root, 'ribos', 'Assigned.sto');
my $ribo_optional = File::Spec->catfile($sample_root, 'ribos', 'diagnostic.txt');
my $missing_alternative = File::Spec->catfile($sample_root, 'ribos', 'missing.txt');
write_file($ribo_optional, "diagnostic\n");
write_file($ribo_result, '');
my $components = {
	ribofind => completion_component_evidence(
		requested => 1,
		checks => [
			{id => 'assignment', kind => 'exists', path => $ribo_result},
			{id => 'diagnostic', kind => 'nonempty',
				path => $ribo_optional, required => 0},
			{id => 'alternative', kind => 'exists_any',
				paths => [$missing_alternative, $ribo_optional], required => 0},
		],
	),
	kraken => completion_component_evidence(requested => 0),
	read_cleaning => completion_component_evidence(
		requested => 1, applicable => 0,
		reason => 'transient output',
	),
};
# RiboFind exposes these derived values while scheduling; the generic checks are
# the authoritative persisted evidence, so the serializer must not duplicate them.
$components->{ribofind}{profile_complete} = 1;
$components->{ribofind}{taxonomy_complete} = 1;
my $outcome = {
	status => 'completed',
	input_size_mb => {primary => 10, supplementary => 2, total => 12},
};
my $path = write_sample_completion(
	root => $sample_root,
	sample => 'S1',
	request_signature => $signature,
	present_assembly => 1,
	components => $components,
	outcome => $outcome,
	metagstats => stats_record(),
);
is($path, sample_completion_path($sample_root, 'S1'),
	'sentinel filename includes the sample id');
like($path, qr/MF4\.sentinel\.S1\.json$/,
	'sentinel uses the requested MF4 filename');
ok(-s $path, 'completion sentinel is published nonempty');
ok(!-e "$path.tmp.$$", 'atomic publication leaves no temporary sentinel');

my ($record, $error) = read_sample_completion(
	root => $sample_root, sample => 'S1', request_signature => $signature,
	expected_components => $components,
);
is($error, '', 'matching sentinel is valid');
is($record->{schema_version}, 3, 'sentinel uses schema version 3');
is($record->{contracts}{components}{schema_version}, 1,
	'component evidence has an independent version');
is($record->{contracts}{metagstats}{schema_version}, 1,
	'statistics evidence has an independent version');
ok(!exists($record->{contracts}{components}{inventory}),
	'component keys are not repeated in a contract inventory');
ok(!exists($record->{contracts}{metagstats}{fields})
	&& !exists($record->{contracts}{metagstats}{families}),
	'statistics keys are not repeated in contract inventories');
is($record->{sample}, 'S1', 'sentinel records its sample identity');
is($record->{present_assembly}, 1, 'sentinel records assembly availability');
ok(!exists($record->{empty_sample}) && !exists($record->{empty_input_size_mb}),
	'empty-sample facts are not duplicated outside the outcome');
ok(!exists($record->{components}{ribofind}{status}),
	'component status is derived rather than persisted');
ok(!exists($record->{components}{ribofind}{profile_complete})
	&& !exists($record->{components}{ribofind}{taxonomy_complete}),
	'RiboFind derived booleans are not duplicated in the sentinel');
is($record->{components}{ribofind}{checks}{assignment}{size_bytes}, 0,
	'zero-byte completion stones record their file size');
is($record->{components}{ribofind}{checks}{diagnostic}{size_bytes},
	length("diagnostic\n"), 'nonempty outputs record their file size');
is($record->{components}{ribofind}{checks}{alternative}{matched_path},
	$ribo_optional, 'multi-path checks record which file matched');
is($record->{components}{ribofind}{checks}{alternative}{size_bytes},
	length("diagnostic\n"), 'the matched alternative records its file size');
ok(!exists($record->{components}{ribofind}{checks}{assignment}{required}),
	'the default required=true value is omitted');
is($record->{components}{ribofind}{checks}{diagnostic}{required}, 0,
	'optional checks remain explicit');
ok(!exists($record->{metagstats}{field_availability}),
	'metagStats values are not duplicated by an availability map');
is($record->{metagstats}{families}{input}, 'complete',
	'statistics families persist only their status');
ok(!completion_record_needs_refresh($record, $components),
	'newly written sentinel has complete file-size evidence');

my $without_size = read_json($path);
delete $without_size->{components}{ribofind}{checks}{assignment}{size_bytes};
write_file($path, JSON::PP->new->canonical(1)->encode($without_size));
($record, $error) = read_sample_completion(
	root => $sample_root, sample => 'S1', request_signature => $signature,
	expected_components => $components,
);
is($error, '', 'a valid completed sample remains closed while sizes are backfilled');
ok(completion_record_needs_refresh($record, $components),
	'a completed sentinel missing a file size requests a lightweight refresh');
write_sample_completion(
	root => $sample_root, sample => 'S1', request_signature => $signature,
	present_assembly => $record->{present_assembly}, components => $components,
	outcome => $record->{outcome}, metagstats => $record->{metagstats},
);
($record, $error) = read_sample_completion(
	root => $sample_root, sample => 'S1', request_signature => $signature,
	expected_components => $components,
);
is($record->{components}{ribofind}{checks}{assignment}{size_bytes}, 0,
	'lightweight refresh restores the missing file size');
ok(!completion_record_needs_refresh($record, $components),
	'refreshed sentinel does not request another rewrite');

unlink $ribo_optional or die $!;
my $optional_missing_components = {
	%{$components},
	ribofind => completion_component_evidence(
		requested => 1,
		checks => [
			{id => 'assignment', kind => 'exists', path => $ribo_result},
			{id => 'diagnostic', kind => 'nonempty',
				path => $ribo_optional, required => 0},
			{id => 'alternative', kind => 'exists_any',
				paths => [$missing_alternative, $ribo_optional], required => 0},
		],
	),
};
($record, $error) = read_sample_completion(
	root => $sample_root, sample => 'S1', request_signature => $signature,
	expected_components => $optional_missing_components,
);
is($error, '', 'removing optional diagnostics does not reopen a sample');
ok(completion_record_needs_refresh($record, $optional_missing_components),
	'removing an optional file clears its stale stored size on the next visit');

unlink $ribo_result or die $!;
my $missing_components = {
	%{$components},
	ribofind => completion_component_evidence(
		requested => 1,
		checks => [{id => 'assignment', kind => 'exists', path => $ribo_result}],
	),
};
($record, $error) = read_sample_completion(
	root => $sample_root, sample => 'S1', request_signature => $signature,
	expected_components => $missing_components,
);
ok(!defined($record), 'missing live component output reopens a completed sample');
like($error, qr/workflow component outputs have changed or are missing/,
	'component-evidence mismatch explains the reopen');
write_file($ribo_result, '');
write_file($ribo_optional, "diagnostic\n");

my $different_signature = completion_request_signature({contract => 4});
($record, $error) = read_sample_completion(
	root => $sample_root, sample => 'S1', request_signature => $different_signature,
);
ok(!defined($record), 'changed workflow does not accept an old sentinel');
like($error, qr/requested workflow differs/,
	'signature mismatch explains why the sample must reopen');
($record, $error) = read_sample_completion(
	root => $sample_root, sample => 'another', request_signature => $signature,
	path => $path,
);
ok(!defined($record), 'sentinel cannot close a different sample');
like($error, qr/belongs to sample 'S1'/,
	'sample mismatch is reported');

eval { sample_completion_path($sample_root, '../unsafe') };
like($@, qr/unsafe sample id/, 'unsafe sample ids cannot escape the sample directory');

my $missing_ribo = {
	ribofind => completion_component_evidence(
		requested => 1,
		checks => [{id => 'assignment', kind => 'exists',
			path => File::Spec->catfile($root, 'missing')}],
	),
};
my $invalid_root = File::Spec->catdir($root, 'invalid');
eval {
	write_sample_completion(
		root => $invalid_root, sample => 'S1', request_signature => $signature,
		components => $missing_ribo, outcome => $outcome,
		metagstats => stats_record(),
	);
};
like($@, qr/requested component 'ribofind' is incomplete/,
	'a normal completion sentinel cannot publish incomplete requested work');

ok(invalidate_sample_completion($sample_root, 'S1'),
	'exact invalidation removes an existing sentinel');
ok(!-e $path, 'invalidated sentinel is absent');
ok(!invalidate_sample_completion($sample_root, 'S1'),
	'exact invalidation is idempotent');

write_sample_completion(
	root => $sample_root, sample => 'S1', request_signature => $signature,
	components => $missing_ribo,
	outcome => {
		status => 'skipped_too_small',
		input_size_mb => {primary => 0.5, supplementary => 0.25, total => 0.75},
		small_sample_threshold_mb => 1,
	},
	metagstats => stats_record(),
);
($record, $error) = read_sample_completion(
	root => $sample_root, sample => 'S1', request_signature => $signature,
);
is($record->{outcome}{status}, 'skipped_too_small',
	'small samples retain their terminal reason');
is($record->{outcome}{input_size_mb}{total}, 0.75,
	'terminal input size adds primary and supplementary input');
is($record->{outcome}{small_sample_threshold_mb}, 1,
	'too-small outcomes retain the threshold that classified them');
ok(!exists($record->{outcome}{small_sample}),
	'too-small state is not duplicated as a boolean');
ok(invalidate_sample_completion($sample_root),
	'root-level invalidation finds the sample-specific sentinel');

write_sample_completion(
	root => $sample_root, sample => 'S1', request_signature => $signature,
	components => $missing_ribo,
	outcome => {
		status => 'skipped_sdm_warning',
		input_size_mb => {primary => 8, supplementary => 1, total => 9},
		sdm_warning => {type => 'invalid_paired_read', log => '/tmp/sdm.log'},
	},
	metagstats => stats_record(),
);
($record, $error) = read_sample_completion(
	root => $sample_root, sample => 'S1', request_signature => $signature,
);
is_deeply($record->{outcome}{sdm_warning},
	{type => 'invalid_paired_read', log => '/tmp/sdm.log'},
	'SDM terminal outcomes retain only the useful warning details');

write_sample_completion(
	root => $sample_root, sample => 'S1', request_signature => $signature,
	components => $missing_ribo,
	outcome => {
		status => 'skipped_cleaned_empty',
		input_size_mb => {primary => 8272.7, supplementary => 0, total => 8272.7},
		cleaned_input_scope => 'primary',
	},
	metagstats => stats_record(),
);
($record, $error) = read_sample_completion(
	root => $sample_root, sample => 'S1', request_signature => $signature,
);
is($record->{outcome}{status}, 'skipped_cleaned_empty',
	'a zero-record cleaned input retains its explicit terminal outcome');
is($record->{outcome}{cleaned_input_scope}, 'primary',
	'cleaned-empty outcome records the read scope that became empty');

my $decoded = read_json($path);
$decoded->{schema_version} = 2;
write_file($path, JSON::PP->new->canonical(1)->encode($decoded));
($record, $error) = read_sample_completion(root => $sample_root, sample => 'S1');
ok(!defined($record), 'non-v3 sentinel is rejected');
like($error, qr/unsupported sentinel schema version/,
	'schema rejection is diagnostic');
$decoded->{schema_version} = 3;
$decoded->{contracts}{components}{schema_version} = 99;
write_file($path, JSON::PP->new->canonical(1)->encode($decoded));
($record, $error) = read_sample_completion(root => $sample_root, sample => 'S1');
ok(!defined($record), 'unknown component-evidence version is rejected');
like($error, qr/component contract is unsupported/,
	'component-contract incompatibility is diagnostic');

write_file($path, "not JSON\n");
($record, $error) = read_sample_completion(root => $sample_root, sample => 'S1');
ok(!defined($record), 'malformed sentinel is rejected');
like($error, qr/invalid JSON/, 'malformed sentinel reports its parse failure');

my $snp_path = File::Spec->catfile(File::Spec->rel2abs('.'), 'Mods', 'SNP.pm');
open my $snp_fh, '<', $snp_path or die "Cannot read $snp_path: $!";
local $/;
my $snp_source = <$snp_fh>;
close $snp_fh or die "Cannot close $snp_path: $!";
like($snp_source,
	qr/sub SNPconsensus_vcf.*?invalidate_sample_completion\(\$SNPIHR->\{sampleRoot\}\)/s,
	'SNP consensus invalidates a detected sample sentinel');
like($snp_source,
	qr/sub SVcall_vcf.*?if \(\$mode ==0 \).*?invalidate_sample_completion\(\$SNPIHR->\{sampleRoot\}\)/s,
	'structural-variant work invalidates a detected sample sentinel');

done_testing;
