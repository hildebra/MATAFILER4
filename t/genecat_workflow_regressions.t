use strict;
use warnings;
use FindBin qw($Bin);
use File::Spec;
use File::Temp qw(tempdir);
use Test::More;

my $root = File::Spec->catdir($Bin, '..');
my $path = File::Spec->catfile($root, 'secScripts', 'geneCat.pl');
open my $fh, '<', $path or die "Cannot read $path: $!";
local $/;
my $source = <$fh>;
close $fh or die "Cannot close $path: $!";

ok(index($source, q{"SNPcaller=s" => \$SNPcaller}) >= 0, 'geneCat accepts the strain SNP caller');
ok(index($source, q{-strains $doStrains -SNPcaller $SNPcaller}) >= 0, 'geneCat forwards the SNP caller to MGS');
ok(index($source, q{$cmd .= "$magPi -mem 150 -GCd " . _shell_quote($OutD)}) >= 0, 'the configured MGS command remains raw while its catalogue path is quoted');
ok(index($source, q{-canopies " . _shell_quote("$canopyExpectedDir/clusters.txt")}) >= 0, 'the Canopy handoff path is shell quoted');
ok(index($source, q{if ($mode eq 'geneCat' && $doMags)}) >= 0, 'MGS-only checker and binner constraints are conditional on doMags');
ok(index($source, q{if $mode eq 'geneCat' && $doStrains && !$doMags}) >= 0, 'strain analysis cannot be silently requested while MGS is disabled');
ok(index($source, q{clusterSingleStep($complStone,$incomplStone,$clnLnStone,$cogStone,$bdir,$OutD,"",$COGdep)}) >= 0, 'single-step fire mode does not duplicate the accumulated command buffer');
ok(index($source, q{clusterMultiStep($complStone,$incomplStone,$clnLnStone,$cogStone,$bdir,$OutD,"",$COGdep)}) >= 0, 'multi-step fire mode does not duplicate the accumulated command buffer');
ok(index($source, q{die "clustering failed\n" if $submitLocal}) >= 0, 'deferred clustering is not checked before its generated commands run');
ok(index($source, q{if (!_stone_valid($FMGstone, $cdhID))}) >= 0, 'fire mode includes marker extraction before marker LCA');
ok(index($source, q{FuncAssign.done.sh}) >= 0, 'functional assignment has a convergence checkpoint job');
ok(index($source, q{join(";", @funcDeps)}) >= 0, 'functional convergence waits for all database matrices');
like($source,
	qr/-mode FuncAssign .*?-fastaSplit \$fastaSplits .*?-FuncMinBitSc \$minBitSc .*?-FuncMinAlLeng \$minAlLeng .*?-FuncMinPercSbjCov \$minPercSbjCov .*?-FuncMinPerID \$minPerID .*?-FuncMinEVal \$minEVal/s,
	'functional child preserves the requested split size and filtering thresholds');
ok(index($source, q{retry_unlink($emapStone, label => 'invalidate stale eggNOG checkpoint')}) >= 0, 'MGS cannot observe an invalid stale eggNOG checkpoint while its replacement runs');
ok(index($source, q{-canopyAutoCorr $canopyAutoCorr -stone $canopyStone}) >= 0, 'deferred Canopy preserves output selection and completion ownership');
ok(index($source, q{my $outF = "$f.emapper.annotations"}) >= 0, 'eggNOG resume checks the output actually consumed by the combiner');
ok(index($source, q{my ($jobDep,$jobCmd) = qsubSystem}) >= 0, 'FOAM workers retain scheduler job IDs rather than command text');
like($source, qr/'publish-catalog',\s*"\$OutD\/\$primaryClusterFNA", "\$OutD\/\$primaryClusterCLS"/, 'publication checkpoint covers both core catalogue outputs');
ok(index($source, q{defined($oldBinner) && $oldBinner =~ /^\d+$/ && $oldBinner == $binSpeciesMG}) >= 0, 'prior MGS output is reused only for the current binner');
ok(index($source, q{$oldGCdIdentity eq $currentGCdIdentity}) >= 0, 'prior MGS output is reused only for the current catalogue identity');
ok(index($source, q{defined($oldClusterID) && $oldClusterID =~ /^\d+$/ && $oldClusterID == $cdhID}) >= 0, 'prior MGS output is reused only for the current cluster identity');
ok(index($source, q{my $oldOutD = _command_option_value($tmpS, '-outD')}) >= 0, 'quoted prior MGS output paths remain resumable');
ok(index($source, q{$MGSoutD = $oldOutD}) >= 0, 'the recovered MGS output spelling is retained rather than canonicalized');

my ($path_helper_source) = $source =~ /^(sub _catalog_path_for_compare \{.*?^\})/ms;
ok(defined($path_helper_source), 'catalogue comparison helper is present');
my $path_helper = defined($path_helper_source)
	? eval('package GeneCatPathCompareTest; use Cwd (); use File::Spec; '
		. $path_helper_source . '; \\&_catalog_path_for_compare')
	: undef;
is($@, '', 'catalogue comparison helper compiles in isolation');

SKIP: {
	skip 'catalogue comparison helper did not compile', 2 unless $path_helper;
	my $tmp = tempdir(CLEANUP => 1);
	my $missing = File::Spec->catdir($tmp, 'not-created');
	is($path_helper->("$missing/"), $path_helper->($missing),
		'non-existing catalogue paths use a stable lexical fallback');

	my $catalogue = File::Spec->catdir($tmp, 'catalogue');
	mkdir $catalogue or die "Cannot create $catalogue: $!";
	my $alias = File::Spec->catfile($tmp, 'catalogue-alias');
	symlink($catalogue, $alias)
		or skip "symlinks unavailable for catalogue identity regression: $!", 1;
	is($path_helper->($alias), $path_helper->($catalogue),
		'existing symlink spellings resolve to the same catalogue identity');
}

done_testing();
