use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::IO_Tamoc_progs qw(getProgPaths setConfigFile truePath);

my $root = tempdir(CLEANUP => 1);
$root =~ s{\\}{/}g;
my $config_one = "$root/config-one.txt";
my $config_two = "$root/config-two.txt";

sub write_config {
	my ($path, $samtools) = @_;
	open my $fh, '>', $path or die "Cannot write $path: $!\n";
	print {$fh} join("\n",
		"customBefore\t[DBDir]/one/[DBDir]",
		"motus2_DB\t$root/user-motus",
		"emptyRequired\t",
		"samtools\t$samtools",
		"MFLRDir\t$root/mf",
		"BINDir\t$root/bin",
		"Rpath\t[MFLRDir]/R",
		"MGSTKDir\t$root/mgs",
		"CONDcmd\tmicromamba",
		"CONDA\t" . 'eval "$(micromamba shell hook -s bash)"',
		"CONDAbaseEnv\tMFF",
		"Rscript\tRscript",
		"DBDir\t$root/db",
	), "\n";
	close $fh or die "Cannot close $path: $!\n";
}

write_config($config_one, 'samtools-one');
write_config($config_two, 'samtools-two');

my ($samtools, $custom_before, $motus_db, $captured);
{
	open my $stdout, '>', \$captured or die $!;
	local *STDOUT = $stdout;
	setConfigFile($config_one);
	$samtools = getProgPaths('samtools');
	$custom_before = getProgPaths('customBefore');
	$motus_db = getProgPaths('motus2_DB');
}
is($samtools, 'samtools-one', 'selected config overrides the bundled program default');
is($custom_before, "$root/db/one/$root/db",
	'foundational placeholders resolve globally even when declared later');
is($motus_db, "$root/user-motus",
	'a blank lower-priority bundled entry cannot overwrite a user value');

write_config($config_one, 'samtools-one-modified');
is(getProgPaths('samtools'), 'samtools-one',
	'repeated lookups use the parsed in-memory configuration cache');

{
	open my $stdout, '>', \$captured or die $!;
	local *STDOUT = $stdout;
	setConfigFile($config_one);
	$samtools = getProgPaths('samtools');
}
is($samtools, 'samtools-one-modified',
	'explicitly reselecting a config invalidates and reloads the cache');

{
	open my $stdout, '>', \$captured or die $!;
	local *STDOUT = $stdout;
	setConfigFile($config_two);
	$samtools = getProgPaths('samtools');
}
is($samtools, 'samtools-two', 'switching config files cannot reuse stale values');

my $error = '';
eval { getProgPaths('emptyRequired') };
$error = $@;
like($error, qr/emptyRequired/, 'a required empty value is rejected');
is(getProgPaths('emptyRequired', 0), '', 'an optional empty value remains supported');

$error = '';
eval { getProgPaths(['samtools', 'missingRequired']) };
$error = $@;
like($error, qr/missingRequired/, 'required array lookups report missing keys');
my $optional = getProgPaths(['samtools', 'missingOptional'], 0);
is_deeply($optional, ['samtools-two', ''], 'optional array lookups preserve positions');

$error = '';
eval { getProgPaths('missingScalar') };
$error = $@;
like($error, qr/\Q$config_two\E/,
	'lookup diagnostics retain the selected user-config filename');

{
	local $ENV{MF4_CONFIG_ROOT} = "$root/env-root";
	local $ENV{MF4_CONFIG_LEAF} = 'leaf';
	is(truePath('$MF4_CONFIG_ROOT/${MF4_CONFIG_LEAF}/$MF4_CONFIG_LEAF'),
		"$root/env-root/leaf/leaf",
		'truePath expands repeated and braced environment variables');
}

$error = '';
eval { truePath('$MF4_CONFIG_MISSING/path') };
$error = $@;
like($error, qr/MF4_CONFIG_MISSING/, 'truePath reports an unset environment variable');

my $mataf_file = File::Spec->catfile($Bin, '..', 'MATAF4.pl');
open my $mataf_fh, '<', $mataf_file or die "Cannot read $mataf_file: $!\n";
my $mataf_source = do { local $/; <$mataf_fh> };
close $mataf_fh or die "Cannot close $mataf_file: $!\n";
my ($defaults_body) = $mataf_source =~ /sub setDefaultMFconfig\s*\{(.*?)\n\}\s*\nsub help\s*\{/s;
ok(defined($defaults_body), 'located setDefaultMFconfig for config-order regression');
unlike($defaults_body // '', qr/getProgPaths\s*\(/,
	'default initialization performs no config lookup before -config is parsed');
like($mataf_source,
	qr/setConfigFile\(\$MFconfig\{configFile\}\);.*?\$MFopt\{baseSDMopt\}\s*=\s*getProgPaths/s,
	'config-backed SDM defaults resolve after selecting the user config');

my $internal_config = File::Spec->catfile($Bin, '..', 'Mods', 'config_internal.txt');
open my $internal_fh, '<', $internal_config or die "Cannot read $internal_config: $!\n";
my $internal_source = do { local $/; <$internal_fh> };
close $internal_fh or die "Cannot close $internal_config: $!\n";
my %bin_tools = (
	rare => 'rtk2', sdm => 'sdm', LCA => 'LCA', readCov => 'rdCover',
	MSAfix => 'MSAfix', canopy => 'cc.bin', clusterMAGs => 'clusterMAGs',
	vcf2fna => 'vcf2fna', trimomatic => 'trimmomatic-0.36.jar',
	treeDistScr => 'distv9.pl', fna2faa => 'fna2faa',
);
for my $key (sort keys %bin_tools){
	like($internal_source, qr/^\Q$key\E\t[^\n]*\[MFLRDir\]\/bin\/[^\n]*\Q$bin_tools{$key}\E/m,
		"$key uses the MF4-owned bin artifact");
}
my $database_config = File::Spec->catfile($Bin, '..', 'Mods', 'config_DBs.txt');
open my $database_fh, '<', $database_config or die "Cannot read $database_config: $!\n";
my $database_source = do { local $/; <$database_fh> };
close $database_fh or die "Cannot close $database_config: $!\n";
like($database_source, qr/^FMGdir\t\[MFLRDir\]\/bin\/fetchMG\//m,
	'fetchMG uses the MF4-owned bin directory');

done_testing;
