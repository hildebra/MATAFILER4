use strict;
use warnings;

use File::Path qw(make_path remove_tree);
use File::Glob qw(bsd_glob);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

# The workflow deletes multi-GB scratch and output trees during initialization.
# Extract the real helper and exercise it against actual directories.
my $source = File::Spec->catfile($Bin, '..', 'secScripts', 'MGS', 'strain_within.pl');
open my $handle, '<', $source or die "Cannot read $source: $!";
my $text = do { local $/; <$handle> };
close $handle or die "Cannot close $source: $!";
my ($helper) = $text =~ /^(sub fastRemoveTree \{.*?^\})/ms;
BAIL_OUT('Cannot extract fastRemoveTree') unless defined $helper;
my ($probeCache) = $text =~ /^(my \$systemRemoveAvailable;)/ms;
BAIL_OUT('Cannot extract the rm capability cache') unless defined $probeCache;
my $loaded = eval "$probeCache\n$helper\n1;";
BAIL_OUT("Cannot load fastRemoveTree: $@") unless $loaded;

my $temporary = tempdir('fast-remove-XXXXXX', TMPDIR => 1, CLEANUP => 1);

sub build_tree {
	my ($root, $breadth, $depth) = @_;
	make_path($root);
	for my $index (1 .. $breadth) {
		my $child = File::Spec->catdir($root, "d$index");
		make_path($child);
		open my $file, '>', File::Spec->catfile($child, "f$index.txt")
			or die "Cannot populate $child: $!";
		print {$file} 'x' x 64;
		close $file or die "Cannot close file in $child: $!";
		build_tree($child, $breadth, $depth - 1) if $depth > 1;
	}
	return $root;
}

is(fastRemoveTree(undef), 0, 'an undefined target removes nothing');
is(fastRemoveTree(File::Spec->catdir($temporary, 'absent')), 0,
	'an absent target removes nothing and does not die');
my $plainFile = File::Spec->catfile($temporary, 'a-file');
open my $file, '>', $plainFile or die $!;
print {$file} 'not a directory';
close $file or die $!;
is(fastRemoveTree($plainFile), 0, 'a plain file is not treated as a tree');
ok(-e $plainFile, 'and is left untouched');

my $tree = build_tree(File::Spec->catdir($temporary, 'scratch'), 4, 3);
ok(-d $tree, 'a nested tree exists before removal');
is(fastRemoveTree($tree, wait => 1), 1, q{the tree is removed});
ok(!-e $tree, q{nothing remains at the original path});
is_deeply([bsd_glob(File::Spec->catfile($temporary, q{*.deleting.*}))], [],
	q{a synchronous removal leaves no parked directory behind});

# A directory parked by an interrupted earlier run must not accumulate.
my $again = build_tree(File::Spec->catdir($temporary, 'scratch'), 2, 2);
my $orphan = File::Spec->catdir($temporary, 'scratch.deleting.1.2.3');
build_tree($orphan, 2, 2);
ok(-d $orphan, 'an orphaned parked directory exists');
cmp_ok(fastRemoveTree($again, wait => 1), q{>=}, 2,
	'the current target and the orphan are both removed');
ok(!-e $again && !-e $orphan, 'neither the target nor the orphan remains');

# Removal must be atomic from the caller's point of view: the path disappears
# rather than being progressively emptied, so partial state is never visible.
my $atomic = build_tree(File::Spec->catdir($temporary, 'atomic'), 3, 2);
my $childBefore = scalar(bsd_glob(File::Spec->catfile($atomic, '*')));
ok($childBefore, 'the tree has children before removal');
fastRemoveTree($atomic);
ok(!-e $atomic,
	q{the target path is gone the moment the call returns, even though the unlink is backgrounded});

# Protected siblings must survive: only the named target is parked and removed.
my $keep = build_tree(File::Spec->catdir($temporary, 'keep'), 2, 1);
my $drop = build_tree(File::Spec->catdir($temporary, 'drop'), 2, 1);
fastRemoveTree($drop, wait => 1);
ok(!-e $drop, q{the requested tree is removed});
ok(-d $keep && scalar(bsd_glob(File::Spec->catfile($keep, '*'))),
	'a sibling directory with a shared parent is untouched');

# A backgrounded unlink must return promptly rather than waiting for the tree.
my $deep = build_tree(File::Spec->catdir($temporary, q{deep}), 6, 3);
my $started = time;
fastRemoveTree($deep);
cmp_ok(time - $started, q{<=}, 5,
	q{backgrounding returns without waiting for the recursive unlink});
ok(!-e $deep, q{and the caller's path is already free});
my $parked_wait = 0;
while ($parked_wait++ < 100
		&& scalar(bsd_glob(File::Spec->catfile($temporary, q{deep.deleting.*})))) {
	select(undef, undef, undef, 0.1);
}
is_deeply([bsd_glob(File::Spec->catfile($temporary, q{deep.deleting.*}))], [],
	q{the backgrounded unlink does complete on its own});

done_testing();
