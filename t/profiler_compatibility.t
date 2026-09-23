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

for my $flag (qw(profileRibosome profileFunct profileMetaphlan profileMOTU2 profileProtal)) {
	}
for my $key (qw(LSUdbFA LSUtax SSUdbFA SSUtax metPhl2_db motus2_DB protal_db)) {
	}

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
