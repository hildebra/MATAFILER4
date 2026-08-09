use strict;
use warnings;

use File::Spec;
use FindBin qw($Bin);
use Test::More;

my $root = File::Spec->catdir($Bin, q{..});
my $script = File::Spec->catfile($root, q{secScripts}, q{phylo}, q{buildTree5.pl});

is(system($^X, q{-I}.$root, q{-c}, $script), 0, q{buildTree5.pl compiles});

done_testing;
