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
	completion_component_evidence
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

sub stats_record {
	return {
		DIR => '/reads/S1',
		values => {RawInputSize => '1.250G', ScaffN50 => 12345},
		field_availability => {RawInputSize => 1, ScaffN50 => 1},
		families => {
			input => {
				requested => 1, applicable => 1, ok => 1,
				status => 'complete',
				fields => {RawInputSize => 1},
			},
		},
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
write_file($ribo_optional, "diagnostic\n");
write_file($ribo_result, '');
my $components = {
	ribofind => completion_component_evidence(
		requested => 1,
		checks => [
			{id => 'assignment', kind => 'exists', path => $ribo_result},
			{id => 'diagnostic', kind => 'nonempty',
				path => $ribo_optional, required => 0},
		],
	),
	kraken => completion_component_evidence(requested => 0),
	read_cleaning => completion_component_evidence(
		requested => 1, applicable => 0,
		reason => 'transient output',
	),
};
my $outcome = {
	status => 'completed',
	input_size_mb => {primary => 10, supplementary => 2, total => 12},
	small_sample => 0,
	sdm_warning => {detected => 0, type => '', log => ''},
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
is($path, sample_completion_path($sample_root),
	'sentinel uses the fixed sample-local filename');
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
is_deeply($record->{contracts}{components}{inventory},
	[sort keys %{$components}],
	'component contract records its complete inventory');
is_deeply($record->{contracts}{metagstats}{fields},
	[qw(RawInputSize ScaffN50)],
	'statistics contract records its field inventory');
is_deeply($record->{contracts}{metagstats}{families}, ['input'],
	'statistics contract records its family inventory');
is($record->{sample}, 'S1', 'sentinel records its sample identity');
is($record->{present_assembly}, 1, 'sentinel records assembly availability');
is($record->{empty_sample}, 0, 'ordinary completion is not marked empty');
is($record->{components}{ribofind}{status}, 'complete',
	'requested persistent component records complete status');
is($record->{components}{ribofind}{checks}{assignment}{ok}, 1,
	'individual persistent file check is recorded');
is($record->{components}{kraken}{status}, 'not_requested',
	'unrequested workflow is explicit in the sentinel');
is($record->{components}{read_cleaning}{status}, 'not_applicable',
	'transient workflow is explicit without blocking closure');
is_deeply($record->{metagstats}{field_availability},
	{RawInputSize => 1, ScaffN50 => 1},
	'sentinel stores per-field statistics availability');
is($record->{metagstats}{families}{input}{status}, 'complete',
	'sentinel stores statistics-family status');
unlink $ribo_optional or die $!;
my $optional_missing_components = {
	%{$components},
	ribofind => completion_component_evidence(
		requested => 1,
		checks => [
			{id => 'assignment', kind => 'exists', path => $ribo_result},
			{id => 'diagnostic', kind => 'nonempty',
				path => $ribo_optional, required => 0},
		],
	),
};
($record, $error) = read_sample_completion(
	root => $sample_root, sample => 'S1', request_signature => $signature,
	expected_components => $optional_missing_components,
);
is($error, '', 'removing an optional diagnostic does not reopen a sample');
ok(defined($record), 'required component evidence remains authoritative');


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

my $different_signature = completion_request_signature({contract => 4});
($record, $error) = read_sample_completion(
	root => $sample_root, sample => 'S1', request_signature => $different_signature,
);
ok(!defined($record), 'changed workflow does not accept an old sentinel');
like($error, qr/requested workflow differs/,
	'signature mismatch explains why the sample must reopen');

($record, $error) = read_sample_completion(
	root => $sample_root, sample => 'another', request_signature => $signature,
);
ok(!defined($record), 'sentinel cannot close a different sample');
like($error, qr/belongs to sample 'S1'/,
	'sample mismatch is reported');

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
		root => $invalid_root, sample => 'S1',
		request_signature => $signature,
		components => $missing_ribo,
		outcome => {status => 'completed'},
		metagstats => stats_record(),
	);
};
like($@, qr/requested component 'ribofind' is incomplete/,
	'a normal completion sentinel cannot publish incomplete requested work');

ok(invalidate_sample_completion($sample_root),
	'invalidation removes an existing sentinel');
ok(!-e $path, 'invalidated sentinel is absent');
ok(!invalidate_sample_completion($sample_root),
	'invalidation is idempotent when no sentinel exists');

write_sample_completion(
	root => $sample_root,
	sample => 'S1',
	request_signature => $signature,
	empty_sample => 1,
	empty_input_size_mb => 0.75,
	components => $missing_ribo,
	outcome => {
		status => 'skipped_too_small',
		input_size_mb => {primary => 0.5, supplementary => 0.25, total => 0.75},
		small_sample => 1,
		sdm_warning => {detected => 0, type => '', log => ''},
	},
	metagstats => stats_record(),
);
($record, $error) = read_sample_completion(
	root => $sample_root, sample => 'S1', request_signature => $signature,
);
is($record->{empty_sample}, 1, 'terminal empty samples are represented in the sentinel');
is($record->{empty_input_size_mb}, 0.75, 'empty-sample input size survives serialization');
is($record->{outcome}{status}, 'skipped_too_small',
	'small samples retain their terminal reason');
is($record->{outcome}{input_size_mb}{total}, 0.75,
	'terminal input size adds primary and supplementary input');
ok(invalidate_sample_completion($sample_root), 'empty-sample sentinel is invalidated normally');

write_sample_completion(
	root => $sample_root, sample => 'S1', request_signature => $signature,
	components => $components, outcome => $outcome, metagstats => stats_record(),
);
open my $json_fh, '<', $path or die $!;
local $/;
my $json = <$json_fh>;
close $json_fh or die $!;
my $decoded = JSON::PP->new->decode($json);
$decoded->{schema_version} = 2;
write_file($path, JSON::PP->new->canonical(1)->encode($decoded));
($record, $error) = read_sample_completion(root => $sample_root, sample => 'S1');
ok(!defined($record), 'legacy schema sentinel is rejected for one-time regeneration');
like($error, qr/unsupported sentinel schema version/,
	'legacy schema rejection is diagnostic');

$decoded->{schema_version} = 3;
$decoded->{contracts}{components}{schema_version} = 99;
write_file($path, JSON::PP->new->canonical(1)->encode($decoded));
($record, $error) = read_sample_completion(root => $sample_root, sample => 'S1');
ok(!defined($record), 'unknown component-evidence version is rejected');
like($error, qr/component contract is unsupported/,
	'component-contract incompatibility is diagnostic');
$decoded->{contracts}{components}{schema_version} = 1;
pop @{$decoded->{contracts}{components}{inventory}};
write_file($path, JSON::PP->new->canonical(1)->encode($decoded));
($record, $error) = read_sample_completion(root => $sample_root, sample => 'S1');
ok(!defined($record), 'a partial component inventory is rejected');
like($error, qr/component contract inventory does not match/,
	'component-inventory corruption is diagnostic');


write_file($path, "not JSON\n");
($record, $error) = read_sample_completion(root => $sample_root, sample => 'S1');
ok(!defined($record), 'malformed sentinel is rejected');
like($error, qr/invalid JSON/, 'malformed sentinel reports its parse failure');

my $snp_path = File::Spec->catfile(File::Spec->rel2abs('.'), 'Mods', 'SNP.pm');
open my $snp_fh, '<', $snp_path or die "Cannot read $snp_path: $!";
my $snp_source = <$snp_fh>;
close $snp_fh or die "Cannot close $snp_path: $!";
like($snp_source,
	qr/sub SNPconsensus_vcf.*?invalidate_sample_completion\(\$SNPIHR->\{sampleRoot\}\)/s,
	'SNP consensus invalidates a detected sample sentinel');
like($snp_source,
	qr/sub SVcall_vcf.*?if \(\$mode ==0 \).*?invalidate_sample_completion\(\$SNPIHR->\{sampleRoot\}\)/s,
	'structural-variant work invalidates a detected sample sentinel');

done_testing;
