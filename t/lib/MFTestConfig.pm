package MFTestConfig;

use strict;
use warnings;

use File::Basename qw(dirname);
use File::Spec;
use Mods::IO_Tamoc_progs qw(getProgPaths setConfigFile);

my $test_dir = dirname(dirname(__FILE__));
my $repo_root = dirname($test_dir);
$ENV{MF4_TEST_ROOT} = $repo_root;
my $test_bin = File::Spec->catdir($test_dir, 'bin');
$ENV{PATH} = $test_bin . ':' . ($ENV{PATH} // '');
my $config_output = '';
{
	open my $config_stdout, '>', \$config_output
		or die "Cannot capture test config output: $!\n";
	local *STDOUT = $config_stdout;
	setConfigFile(File::Spec->catfile($test_dir, 'MATAFILERcfg.txt'));
	getProgPaths('pigz');
}

1;
