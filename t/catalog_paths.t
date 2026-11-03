use strict;
use warnings;

use Cwd qw(abs_path);
use File::Copy qw(copy);
use File::Glob qw(bsd_glob);
use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::CatalogPaths qw(
	catalog_identity catalog_map_manifest catalog_map_specs_match
	resolve_catalog_maps write_catalog_maps filter_catalog_maps
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

# MGS and the compatibility clustering wrapper share the same streaming map
# filter. Preserve map metadata and ordering, including maps emptied of samples.
sub read_map_text {
	my ($path) = @_;
	open my $fh, '<', $path or die "Cannot read map $path: $!";
	local $/;
	return <$fh>;
}
my $filter_source = File::Spec->catdir($tmp, 'filter-source');
my $filter_target = File::Spec->catdir($tmp, 'filtered maps');
make_path($filter_source);
my @filter_inputs = map { File::Spec->catfile($filter_source, "input.$_.txt") } 0 .. 1;
my @input_text = (
	"#SmplID\tPath\r\n#OutPath\t/data/\r\n\nkeep\t/data/keep\r\nempty\t/data/empty\r\n",
	"#SmplID\tPath\nempty2\t/data/empty2\n",
);
for my $index (0 .. $#filter_inputs) {
	open my $fh, '>', $filter_inputs[$index] or die $!;
	print {$fh} $input_text[$index];
	close $fh or die $!;
}
my $filter_spec = join(',', @filter_inputs);
is(filter_catalog_maps($filter_spec, [], $filter_target), $filter_spec,
	'no exclusions reuse the original map specification');
ok(!-e $filter_target, 'no exclusions create no filtered-map directory');
my @filtered = split /,/, filter_catalog_maps($filter_spec, [qw(empty empty2)], $filter_target);
is_deeply(\@filtered,
	[map { File::Spec->catfile($filter_target, "map.$_.txt") } 0 .. 1],
	'filtered maps preserve input order and deterministic names');
is(read_map_text($filtered[0]),
	"#SmplID\tPath\r\n#OutPath\t/data/\r\n\nkeep\t/data/keep\r\n",
	'filtering preserves metadata, blank lines, retained samples and line endings');
is(read_map_text($filtered[1]), "#SmplID\tPath\n",
	'a map with no eligible samples retains its header');
for my $index (0 .. $#filter_inputs) {
	is(read_map_text($filter_inputs[$index]), $input_text[$index],
		"filtering leaves source map $index unchanged");
}
filter_catalog_maps($filter_spec, ['keep'], $filter_target);
is(read_map_text($filtered[0]),
	"#SmplID\tPath\r\n#OutPath\t/data/\r\n\nempty\t/data/empty\r\n",
	'a subsequent filter replaces the previous selection');
is_deeply([bsd_glob(File::Spec->catfile($filter_target, '*.tmp.*'))], [],
	'successful publication leaves no temporary maps');
my $failed_target = File::Spec->catdir($tmp, 'failed-publication');
my $blocking_directory = File::Spec->catdir($failed_target, 'map.0.txt');
make_path($blocking_directory);
eval { filter_catalog_maps($filter_spec, ['empty'], $failed_target) };
like($@, qr/Cannot publish filtered map/, 'publication failures are reported');
ok(-d $blocking_directory, 'failed publication does not remove an existing destination');
is_deeply([bsd_glob(File::Spec->catfile($failed_target, '*.tmp.*'))], [],
	'failed publication cleans its temporary output');

done_testing;
