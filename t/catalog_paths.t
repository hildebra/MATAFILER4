use strict;
use warnings;

use Cwd qw(abs_path);
use File::Copy qw(copy);
use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::CatalogPaths qw(
	catalog_identity catalog_map_manifest catalog_map_specs_match
	resolve_catalog_maps write_catalog_maps
);

my $tmp = tempdir(CLEANUP => 1);
my $catalog = File::Spec->catdir($tmp, 'catalog');
my $logs = File::Spec->catdir($catalog, 'LOGandSUB');
make_path($logs);

my @maps;
for my $index (0 .. 1) {
	my $map = File::Spec->catfile($logs, "map.$index.txt");
	open my $fh, '>', $map or die "Cannot write $map: $!";
	print {$fh} "#SmplID\tPath\nsample$index\t/tmp/sample$index\n";
	close $fh;
	push @maps, $map;
}

my $resolved = write_catalog_maps($catalog, \@maps);
is($resolved, join(',', map { abs_path($_) } @maps),
	'map manifest retains every catalog-local map in order');
is(resolve_catalog_maps($catalog), $resolved,
	'multi-map inmap.txt resolves to the complete map set');
ok(catalog_map_specs_match(join(',', @maps), $resolved),
	'comma-separated map lists match after canonicalizing every entry');
ok(catalog_map_specs_match("$logs/./map.0.txt,$logs/./map.1.txt", $resolved),
	'comma-separated map lists tolerate equivalent path spellings');
my $relocated_dir = File::Spec->catdir($tmp, 'relocated-maps');
make_path($relocated_dir);
my @relocated_maps;
for my $source (@maps) {
	my $target = File::Spec->catfile($relocated_dir, (File::Spec->splitpath($source))[2]);
	copy($source, $target) or die "Cannot copy $source to $target: $!";
	push @relocated_maps, $target;
}

ok(catalog_map_specs_match(join(',', @relocated_maps), $resolved),
	'map lists with identical ordered contents match after a map relocation');
open my $changed_fh, '>>', $relocated_maps[1]
	or die "Cannot modify relocated map: $!";
print {$changed_fh} "changed\n";
close $changed_fh;
ok(!catalog_map_specs_match(join(',', @relocated_maps), $resolved),
	'map lists with changed contents still fail the resume guard');
ok(catalog_map_specs_match(join(",\n", @maps), $resolved),
	'comma/newline map lists preserve every map entry');

ok(!catalog_map_specs_match($maps[0], $resolved),
	'map-list comparison rejects a missing map entry');
ok(!catalog_map_specs_match('', ''),
	'empty map specifications never match');

open my $manifest_fh, '<', catalog_map_manifest($catalog)
	or die "Cannot read map manifest: $!";
my @manifest_lines = <$manifest_fh>;
close $manifest_fh;
chomp @manifest_lines;
is_deeply(\@manifest_lines, ['map.0.txt', 'map.1.txt'],
	'inmap.txt stores relocatable catalog-local map paths');
my $legacy_log = File::Spec->catfile($logs, 'GCmaps.inf');
open my $legacy_fh, '<', $legacy_log
	or die "Cannot read compatibility map log: $!";
my $legacy_line = <$legacy_fh> // '';
close $legacy_fh;
chomp $legacy_line;
is($legacy_line, $resolved,
	'GCmaps.inf is still written as a compatibility log');

my $identity = catalog_identity($catalog);
like($identity, qr/\A[0-9a-f]{64}\z/, 'catalog identity is a SHA-256 value');
is(catalog_identity($catalog), $identity,
	'repeated identity retrieval returns the stored value');

my $moved = File::Spec->catdir($tmp, 'moved-catalog');
rename($catalog, $moved) or die "Cannot move test catalog: $!";
is(catalog_identity($moved), $identity,
	'catalog identity remains stable when the complete catalog is moved');
my $moved_maps = resolve_catalog_maps($moved);
is(scalar(split /,/, $moved_maps), 2,
	'catalog-local map manifests still resolve after the catalog is moved');

my $legacy = File::Spec->catdir($tmp, 'legacy-catalog');
my $legacy_logs = File::Spec->catdir($legacy, 'LOGandSUB');
make_path($legacy_logs);
for my $index (0 .. 1) {
	my $source = File::Spec->catfile($moved, 'LOGandSUB', "map.$index.txt");
	my $target = File::Spec->catfile($legacy_logs, "map.$index.txt");
	open my $in, '<', $source or die "Cannot read $source: $!";
	open my $out, '>', $target or die "Cannot write $target: $!";
	print {$out} $_ while <$in>;
	close $in;
	close $out;
}
my $legacy_resolved = resolve_catalog_maps($legacy);
ok(-s catalog_map_manifest($legacy),
	'catalog-local map copies are migrated automatically to inmap.txt');
is(scalar(split /,/, $legacy_resolved), 2,
	'automatic migration preserves every copied map');

done_testing;
