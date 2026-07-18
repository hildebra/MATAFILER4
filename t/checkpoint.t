use strict;
use warnings;
use Test::More;
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use lib File::Spec->catdir($Bin, '..');
use Mods::Checkpoint qw(write_checkpoint checkpoint_valid read_checkpoint);

my $dir = tempdir(CLEANUP => 1);
my $output = File::Spec->catfile($dir, 'result.txt');
my $stone = File::Spec->catfile($dir, 'stage.stone');
open my $out, '>', $output or die $!;
print {$out} "complete\n";
close $out or die $!;

write_checkpoint($stone,
	parameters => { cluster_id => 97, stage => 'test' },
	outputs => [$output],
);
ok(checkpoint_valid($stone, parameters => { cluster_id => 97 }), 'matching manifest is valid');
ok(!checkpoint_valid($stone, parameters => { cluster_id => 95 }), 'parameter mismatch invalidates checkpoint');
is(read_checkpoint($stone)->{outputs}[0]{size}, -s $output, 'manifest records output size');

open $out, '>>', $output or die $!;
print {$out} "changed\n";
close $out or die $!;
ok(!checkpoint_valid($stone, parameters => { cluster_id => 97 }), 'changed output invalidates checkpoint');

my $legacy = File::Spec->catfile($dir, 'legacy.stone');
open my $legacy_fh, '>', $legacy or die $!;
close $legacy_fh or die $!;
ok(checkpoint_valid($legacy, parameters => { cluster_id => 97 }), 'empty legacy checkpoint remains compatible');

done_testing();
