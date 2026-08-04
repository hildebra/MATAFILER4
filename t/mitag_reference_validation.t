use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

my $root = File::Spec->catdir($Bin, '..');
my $script = File::Spec->catfile(
	$root, 'secScripts', 'miTag', 'catchLSUSSU.pl',
);

is(system($^X, '-I'.$root, '-c', $script), 0,
	'catchLSUSSU compiles');

open my $source_handle, '<', $script
	or die "Cannot read $script: $!";
my $source = do { local $/; <$source_handle> };
close $source_handle or die "Cannot close $script: $!";

my ($reference_validator) = $source =~
	/(sub validateSortmernaReference\{.*?^\})/ms;
my ($index_validator) = $source =~
	/(sub validateSortmernaIndex\{.*?^\})/ms;
ok(defined $reference_validator, 'reference validator is present');
ok(defined $index_validator, 'index validator is present');

my $loaded = eval join("\n",
	'package MitagExistenceValidation;',
	'use strict;',
	'use warnings;',
	$reference_validator // '',
	$index_validator // '',
	'1;',
);
ok($loaded, 'existence validators can be exercised independently')
	or diag($@);

sub trapped_error {
	my ($callback) = @_;
	local $@;
	eval { $callback->(); 1 };
	return $@;
}

my $tmp = tempdir(CLEANUP => 1);
my $empty_file = File::Spec->catfile($tmp, 'empty.fasta');
open my $empty_handle, '>', $empty_file
	or die "Cannot create $empty_file: $!";
close $empty_handle or die "Cannot close $empty_file: $!";
my $existing_directory = File::Spec->catdir($tmp, 'existing-index');
mkdir $existing_directory or die "Cannot create $existing_directory: $!";
my $missing = File::Spec->catfile($tmp, 'missing');

is(trapped_error(sub {
	MitagExistenceValidation::validateSortmernaReference(
		'SSUdbFAsrt', "$empty_file:$existing_directory",
	);
}), '', 'existing reference paths are accepted without type or access probes');
like(trapped_error(sub {
	MitagExistenceValidation::validateSortmernaReference(
		'SSUdbFAsrt', "$empty_file:$missing",
	);
}), qr/does not exist: \Q$missing\E/,
	'a missing reference path is rejected');

is(trapped_error(sub {
	MitagExistenceValidation::validateSortmernaIndex(
		'SSUidx', $empty_file,
	);
}), '', 'an existing index path is accepted without directory or access probes');
like(trapped_error(sub {
	MitagExistenceValidation::validateSortmernaIndex(
		'SSUidx', $missing,
	);
}), qr/does not exist: \Q$missing\E/,
	'a missing index path is rejected');

unlike($reference_validator, qr/(?:-[frsx]\s+\$reference|open\s+my\s+\$referenceHandle|read\s*\()/,
	'reference preflight contains no type, permission, open, or read probe');
unlike($index_validator, qr/-[drx]\s+\$indexDirectory/,
	'index preflight contains no directory, read, or search-permission probe');
like($reference_validator, qr/unless -e \$reference/,
	'reference preflight uses an existence check');
like($index_validator, qr/unless -e \$indexDirectory/,
	'index preflight uses an existence check');

my $lca_script = File::Spec->catfile(
	$root, 'secScripts', 'miTag', 'lotus_LCA_blast3.pl',
);
is(system($^X, '-I'.$root, '-c', $lca_script), 0,
	'lotus_LCA_blast3 compiles');
open my $lca_source_handle, '<', $lca_script
	or die "Cannot read $lca_script: $!";
my $lca_source = do { local $/; <$lca_source_handle> };
close $lca_source_handle or die "Cannot close $lca_script: $!";
ok(index($lca_source, 'unless -e $dbfa[$index]') >= 0,
	'downstream reference database preflight uses an existence check');
ok(index($lca_source, '-e $dbtax[$index]') >= 0,
	'downstream taxonomy database preflight uses an existence check');
ok(index($lca_source, 'unless -s $dbfa[$index]') < 0
	&& index($lca_source, '-s $dbtax[$index]') < 0,
	'downstream database preflight contains no size probe');

done_testing;
