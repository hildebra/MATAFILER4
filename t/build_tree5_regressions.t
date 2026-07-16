use strict;
use warnings;

use File::Spec;
use FindBin qw($Bin);
use Test::More;

my $root = File::Spec->catdir($Bin, '..');
my $script = File::Spec->catfile($root, 'secScripts', 'phylo', 'buildTree5.pl');

open my $fh, '<', $script or die "Cannot read $script: $!";
my $source = do { local $/; <$fh> };
close $fh;

my $compile_status = system($^X, '-I'.$root, '-c', $script);
is($compile_status, 0, 'buildTree5.pl compiles');

like($source, qr/-outD is required/, 'an output directory is explicitly required');
like($source, qr/Refusing to use filesystem root/, 'filesystem roots are rejected as output directories');
like($source, qr/buildTree5_\$\{tmpTag\}_\$\$/, 'work is isolated in a process-owned temporary directory');
like($source, qr/safeRemoveTree\(\$tmpD, \$tmpBase\)/, 'cleanup is limited to the owned temporary directory');

unlike($source, qr/touch \$IQtreef/, 'an empty IQ-TREE checkpoint is not manufactured');
unlike($source, qr/\$calcSyn\s*=\s*0\s*;\s*\$calcNonSyn\s*=\s*0\s*;\s*\n\s*if \(\$cogCats/, 
	'synonymous and nonsynonymous tree options are retained');
like($source, qr/my \$lengthInNt = .*\? \$totalNTs\{\$sp\} : \$totalNTs\{\$sp\} \* 3;/,
	'NTfiltCount compares a nucleotide-equivalent length');
like($source, qr/systemW\(\$cmd1\."\\n"\.\$cmd2\."\\n"\)/,
	'alignment and post-filter commands use checked execution');
like($source, qr/MSA command completed without producing/, 'alignment output is verified');

like($source, qr/\$pigzBin -d .*\$partiF\.gz/, 'compressed partition restoration names the gzip file');
like($source, qr/for my \$disM \(\@subfls\)/, 'all discovered distance matrices are merged');
like($source, qr/\$ffd\{\$k\} = 4/, 'fourfold degeneracy is classified by codon family');

for my $environment_name (qw(
	MF4_GUBBINS_BIN MF4_CLONALFRAMEML_BIN MF4_FASTGEAR_BIN
	MF4_FASTGEAR_MATLAB_BIN MF4_FASTGEAR_PARAM_FILE
)) {
	like($source, qr/\Q$environment_name\E/, "$environment_name documents dormant-tool reactivation");
}

done_testing();
