use strict;
use warnings;

use File::Spec;
use FindBin qw($Bin);
use Test::More;

my $root = File::Spec->catdir($Bin, q{..});
my $script = File::Spec->catfile($root, q{secScripts}, q{phylo}, q{buildTree5.pl});

open my $script_handle, '<', $script or die "Cannot read $script: $!";
my $script_text = do { local $/; <$script_handle> };
close $script_handle or die "Cannot close $script: $!";
my $coordinate_bounds_checks = () = $script_text =~ /next if \$position >= length\(\$sequence\);/g;
cmp_ok($coordinate_bounds_checks, '>=', 2,
	'taxon-aware raw and alignment coordinate scorers safely skip uneven sequence tails');
is(system($^X, q{-I}.$root, q{-c}, $script), 0, q{buildTree5.pl compiles});

done_testing;
