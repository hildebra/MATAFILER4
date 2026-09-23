use strict;
use warnings;
use FindBin qw($Bin);
use File::Spec;
use File::Temp qw(tempdir);
use Test::More;

my $root = File::Spec->catdir($Bin, '..');
my $path = File::Spec->catfile($root, 'secScripts', 'geneCat.pl');
open my $fh, '<', $path or die "Cannot read $path: $!";
local $/;
my $source = <$fh>;
close $fh or die "Cannot close $path: $!";

my ($path_helper_source) = $source =~ /^(sub _catalog_path_for_compare \{.*?^\})/ms;
ok(defined($path_helper_source), 'catalogue comparison helper is present');
my $path_helper = defined($path_helper_source)
	? eval('package GeneCatPathCompareTest; use Cwd (); use File::Spec; '
		. $path_helper_source . '; \\&_catalog_path_for_compare')
	: undef;
is($@, '', 'catalogue comparison helper compiles in isolation');

SKIP: {
	skip 'catalogue comparison helper did not compile', 2 unless $path_helper;
	my $tmp = tempdir(CLEANUP => 1);
	my $missing = File::Spec->catdir($tmp, 'not-created');
	is($path_helper->("$missing/"), $path_helper->($missing),
		'non-existing catalogue paths use a stable lexical fallback');

	my $catalogue = File::Spec->catdir($tmp, 'catalogue');
	mkdir $catalogue or die "Cannot create $catalogue: $!";
	my $alias = File::Spec->catfile($tmp, 'catalogue-alias');
	symlink($catalogue, $alias)
		or skip "symlinks unavailable for catalogue identity regression: $!", 1;
	is($path_helper->($alias), $path_helper->($catalogue),
		'existing symlink spellings resolve to the same catalogue identity');
}

done_testing();
