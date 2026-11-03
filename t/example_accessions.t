use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::GenoMetaAss qw(readMap);

my $repository = File::Spec->catdir($Bin, '..');
my $map_file = File::Spec->catfile(
	$repository, 'examples', 'maps', 'testAccessions.map',
);
my $run_script = File::Spec->catfile(
	$repository, 'examples', '3.runMATAFILER_accessions.mfc',
);

ok(-s $map_file, 'Example 3 mapping file is shipped');
ok(-s $run_script, 'Example 3 runner is shipped');

my $temporary_root = tempdir(CLEANUP => 1);
make_path(
	File::Spec->catdir($temporary_root, 'examples', 'data'),
	File::Spec->catdir($temporary_root, 'examples', 'output'),
);

my ($map, $groups);
{
	local $ENV{MF4DIR} = $temporary_root;
	($map, $groups) = readMap($map_file, -1, {}, {}, 0);
}

is_deeply(
	$map->{opt}{smpl_order},
	[qw(LocalEx1 ENApublic SRApublic)],
	'Example 3 preserves the local, ENA, and SRA sample order',
);
is($map->{LocalEx1}{prefix}, 'SRR8797712',
	'the local example row reuses an existing Example 1 file prefix');
is($map->{ENApublic}{ENA_download}, 'ERR10009595',
	'the ENA example uses the verified small paired WGS run');
is($map->{SRApublic}{SRA_download}, 'SRR10090860',
	'the SRA example uses the verified small paired shotgun WGS run');

for my $sample (qw(ENApublic SRApublic)) {
	is($map->{$sample}{dir}, '', "$sample has no Path designation");
	is($map->{$sample}{prefix}, '', "$sample has no SmplPrefix designation");
	is($map->{$sample}{SeqTechDeclared}, '',
		"$sample leaves SeqTech available for archive inference");
	ok($map->{$sample}{hasPrimaryRds},
		"$sample registers its accession as primary reads");
}

for my $group (qw(ENAreads LocalReads SRAreads)) {
	ok(exists($groups->{$group}) && $groups->{$group}{CntAimAss} == 1,
		"$group is a separate one-sample assembly group");
}

is(system('bash', '-n', $run_script), 0,
	'Example 3 runner passes Bash syntax validation');

done_testing();
