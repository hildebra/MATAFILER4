use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use IO::Compress::Gzip qw(gzip $GzipError);
use Test::More;

my $root = File::Spec->catdir($Bin, '..');

sub slurp {
	my ($path) = @_;
	open my $fh, '<', $path or die "Cannot read $path: $!";
	my $text = do { local $/; <$fh> };
	close $fh or die "Cannot close $path: $!";
	return $text;
}

my $source = slurp(File::Spec->catfile($root, 'MATAF4.pl'));
my ($motus_mapping) = $source =~ /sub mOTU2Mapping\{(.*?)\n\}\n\n\nsub prepKraken/s;
ok(defined($motus_mapping), 'found the mOTUs mapping implementation');
like($motus_mapping, qr/getProgPaths\("motus2_DB",0\).*?getProgPaths\("motus2"\)/s,
	'mOTUs executable and database are resolved through getProgPaths');
like($motus_mapping, qr/join\(" ",\@car1\).*?join\(" ",\@car2\).*?join\(" ",\@sar\)/s,
	'mOTUs 4 receives multiple inputs as separate arguments');
like($motus_mapping, qr/\$m2Bin profile .*?-n \$smp .*?-t \$Ncore -o \$motuOut/s,
	'mOTUs 4 profile command uses its required output option');
like($motus_mapping, qr/my \$motuOut = "\$tmpD\/\$smp\.motus".*?gzip -c \$motuOut > \$finOutD/s,
	'mOTUs intermediates stay in scratch and only the compressed profile is published');
unlike($motus_mapping, qr/-q -c|-t \$Ncore -q/,
	'removed mOTUs 3 options are not emitted');

my ($metaphlan_prep) = $source =~ /sub prepMetaphlan\{(.*?)\n\}\n\n\nsub metphlanMapping/s;
my ($metaphlan_mapping) = $source =~ /sub metphlanMapping\{(.*?)(?=\nsub TaxaTarget\{)/s;
ok(defined($metaphlan_prep), 'found the MetaPhlAn version preparation');
ok(defined($metaphlan_mapping), 'found the MetaPhlAn mapping implementation');
like($metaphlan_prep, qr/MetaPhlanModernCLI.*?MPverL == 4.*?MPverParts\[1\].*?>= 2/s,
	'MetaPhlAn 4.2 and later select the current interface');
like($metaphlan_mapping, qr/--db_dir \$mpDB1 --index \$mpDB2/,
	'MetaPhlAn 4.2 receives the renamed database options');
like($metaphlan_mapping, qr/if \(\$MFglobal\{MetaPhlanModernCLI\}\).*?--nreads .*? -o.*?else.*?--unclassified_estimation --add_viruses/s,
	'removed MetaPhlAn options are confined to the legacy 4.0/4.1 branch');
like($metaphlan_mapping, qr/--input_type sam \$taxinfo \$vXparams/,
	'the current supported SAM interface remains in use');

my $checkm2_env = slurp(File::Spec->catfile($root, 'helpers', 'install', 'checkm2.yml'));
like($checkm2_env, qr/^\s*-\s+bioconda::metaphlan=4\.2\.6\s*$/m,
	'MetaPhlAn is pinned to 4.2.6');
my $motus_env = slurp(File::Spec->catfile($root, 'helpers', 'install', 'motus.yml'));
like($motus_env, qr/^name:\s+MF4motus\s*$/m, 'mOTUs has a separate environment');
like($motus_env, qr/^\s*-\s+bioconda::motus=4\.1\.0\s*$/m, 'mOTUs is pinned to 4.1.0');
like($motus_env, qr/^\s*-\s+bioconda::bwa>=0\.7\.19\s*$/m, 'mOTUs BWA requirement is explicit');
like($motus_env, qr/^\s*-\s+bioconda::vsearch>=2\.30\.4\s*$/m, 'mOTUs VSEARCH requirement is explicit');

my $gtdb_env = slurp(File::Spec->catfile($root, 'helpers', 'install', 'GTDBTK.yml'));
unlike($gtdb_env, qr/^\s*-\s+motus=/m, 'the stale mOTUs 3 pin is absent from GTDB-Tk');
my $internal = slurp(File::Spec->catfile($root, 'Mods', 'config_internal.txt'));
like($internal, qr/^motus2\tmotus\tenv:MF4motus$/m,
	'getProgPaths runs mOTUs from the dedicated environment');
my $databases = slurp(File::Spec->catfile($root, 'Mods', 'config_DBs.txt'));
like($databases, qr/^metPhl2_db\t\[DBDir\]\/MP4\/mpa_vJan25_CHOCOPhlAnSGB_202503$/m,
	'default MetaPhlAn prefix matches the installed database');
like($databases, qr/^motus2_DB\t\[DBDir\]\/mOTUs$/m,
	'default mOTUs path is the parent of db_mOTU');

my $installer = slurp(File::Spec->catfile($root, 'helpers', 'install', 'installer.sh'));
like($installer, qr/ensure_environment MF4motus .*?verify_environment_tools MF4motus/s,
	'the installer creates and verifies the mOTUs environment');
like($installer, qr/metaphlan --install --db_dir "\$MP4DB" --index "\$METAPHLAN_DB_INDEX"/,
	'the installer uses the MetaPhlAn 4.2 database interface');
like($installer, qr/motus downloadMGDB -db "\$MOTUSDB"/,
	'the installer downloads the mOTUs marker-gene database');

my $tutorial = slurp(File::Spec->catfile($root, 'docs', 'profiling_tutorial.md'));
for my $flag (qw(profileRibosome profileFunct profileMetaphlan profileMOTU2 profileProtal)) {
	like($tutorial, qr/-\Q$flag\E\b/, "tutorial documents -$flag");
}
for my $key (qw(LSUdbFA LSUtax SSUdbFA SSUtax metPhl2_db motus2_DB protal_db)) {
	like($tutorial, qr/\b\Q$key\E\b/, "tutorial documents database key $key");
}
like($tutorial, qr/eggNOG40_path_DB.*?KEGG_path_DB.*?URE_path_DB/s,
	'tutorial documents the functional database key range');

# mOTUs 4.1 writes a non-commented second header and d__ GTDB taxonomy.
# Exercise the actual MATAFILER merger against that schema.
my $tmp = tempdir(CLEANUP => 1);
my $input_dir = File::Spec->catdir($tmp, 'mOTU2');
make_path($input_dir);
my $profile = <<"PROFILE";
#tool_version=4.1.0\tdatabase_version=4.1\tmin_alignment_length=75\tmin_mgcs=3\tcount_mode=INSERT_SCALED\tvalue_type=counts
mOTU\tTaxonomy\tSampleA
mOTUv4.0_000001\td__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__Lactobacillus testii\t12
mOTUv4.0_000002\td__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli\t3
PROFILE
my $profile_path = File::Spec->catfile($input_dir, 'SampleA.motu2.tab.gz');
gzip(\$profile => $profile_path) or die "Cannot write $profile_path: $GzipError";
my $merger = File::Spec->catfile($root, 'secScripts', 'composition', 'mrgMotu2.pl');
my $test_config = File::Spec->catfile($root, 't', 'MATAFILERcfg.txt');
my $runner = 'use Mods::IO_Tamoc_progs qw(setConfigFile); '
	. 'my $config = shift @ARGV; my $script = shift @ARGV; '
	. 'setConfigFile($config); my $ok = do $script; '
	. 'die $@ if $@; die $! unless defined $ok;';
local $ENV{MF4_TEST_ROOT} = File::Spec->rel2abs($root);
is(system($^X, "-I$root", '-e', $runner, $test_config, $merger, $input_dir, 1), 0,
	'mOTUs merger accepts a 4.1 profile');
my $kingdom = slurp(File::Spec->catfile($tmp, 'm2.kingdom.txt'));
like($kingdom, qr/^Bacteria\t15$/m,
	'd__ is stripped and counts are summarized at the first rank');
my $species = slurp(File::Spec->catfile($tmp, 'm2.species.txt'));
like($species, qr/^Bacteria;Bacillota;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;Lactobacillus testii\t12$/m,
	'all seven mOTUs 4 GTDB ranks reach the species matrix');

done_testing();
