use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use Test::More;

use lib File::Spec->catdir(File::Spec->rel2abs('.'));
use Mods::SampleCompletion qw(
	sample_completion_path completion_request_signature
	read_sample_completion write_sample_completion invalidate_sample_completion
);

my $root = tempdir(CLEANUP => 1);
my $sample_root = File::Spec->catdir($root, 'S1');
my $request = {
	contract => 1,
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

my $path = write_sample_completion(
	root => $sample_root,
	sample => 'S1',
	request_signature => $signature,
	present_assembly => 1,
	metagstats => {
		DIR => '/reads/S1',
		values => {RawInputSize => '1.250G', ScaffN50 => 12345},
	},
);
is($path, sample_completion_path($sample_root),
	'sentinel uses the fixed sample-local filename');
ok(-s $path, 'completion sentinel is published nonempty');
ok(!-e "$path.tmp.$$", 'atomic publication leaves no temporary sentinel');

my ($record, $error) = read_sample_completion(
	root => $sample_root, sample => 'S1', request_signature => $signature,
);
is($error, '', 'matching sentinel is valid');
is($record->{sample}, 'S1', 'sentinel records its sample identity');
is($record->{present_assembly}, 1, 'sentinel records assembly availability');
is($record->{empty_sample}, 0, "ordinary completion is not marked empty");
is_deeply($record->{metagstats}, {
	DIR => '/reads/S1',
	values => {RawInputSize => '1.250G', ScaffN50 => 12345},
}, 'sentinel stores the complete metagStats sample object');

my $different_signature = completion_request_signature({contract => 2});
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

ok(invalidate_sample_completion($sample_root),
	'invalidation removes an existing sentinel');
ok(!-e $path, 'invalidated sentinel is absent');
ok(!invalidate_sample_completion($sample_root),
	'invalidation is idempotent when no sentinel exists');

write_sample_completion(
	root => $sample_root,
	sample => "S1",
	request_signature => $signature,
	empty_sample => 1,
	empty_input_size_mb => 0.75,
	metagstats => {DIR => "/reads/S1", values => {RawInputSize => "0.001G"}},
);
($record, $error) = read_sample_completion(
	root => $sample_root, sample => "S1", request_signature => $signature,
);
is($record->{empty_sample}, 1, "terminal empty samples are represented in the sentinel");
is($record->{empty_input_size_mb}, 0.75, "empty-sample input size survives serialization");
ok(invalidate_sample_completion($sample_root), "empty-sample sentinel is invalidated normally");

open my $bad, '>', $path or die $!;
print {$bad} "not JSON\n";
close $bad or die $!;
($record, $error) = read_sample_completion(root => $sample_root, sample => 'S1');
ok(!defined($record), 'malformed sentinel is rejected');
like($error, qr/invalid JSON/, 'malformed sentinel reports its parse failure');

my $snp_path = File::Spec->catfile(File::Spec->rel2abs("."), "Mods", "SNP.pm");
open my $snp_fh, "<", $snp_path or die "Cannot read $snp_path: $!";
local $/;
my $snp_source = <$snp_fh>;
close $snp_fh or die "Cannot close $snp_path: $!";
like($snp_source,
	qr/sub SNPconsensus_vcf.*?invalidate_sample_completion\(\$SNPIHR->\{sampleRoot\}\)/s,
	"SNP consensus invalidates a detected sample sentinel");
like($snp_source,
	qr/sub SVcall_vcf.*?if \(\$mode ==0 \).*?invalidate_sample_completion\(\$SNPIHR->\{sampleRoot\}\)/s,
	"structural-variant work invalidates a detected sample sentinel");

done_testing;
