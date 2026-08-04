use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::phyloTools qw(
	runQItree iqtreeOutputComplete cleanupIQTreeTransients
);

sub write_file {
	my ($path, $text) = @_;
	open my $fh, '>', $path or die "Cannot write $path: $!";
	print {$fh} $text;
	close $fh or die "Cannot close $path: $!";
}

sub slurp {
	my ($path) = @_;
	open my $fh, '<', $path or die "Cannot read $path: $!";
	local $/;
	my $text = <$fh> // '';
	close $fh or die "Cannot close $path: $!";
	return $text;
}

my $tmp = tempdir(CLEANUP => 1);
my $alignment = File::Spec->catfile($tmp, 'alignment.fna');
my $prefix = File::Spec->catfile($tmp, 'IQtree_allsites');
write_file($alignment, ">A\nAAAA\n>B\nAAAA\n>C\nAAAT\n");
write_file("$prefix.treefile", "(A:0,B:0);\n");
write_file("$prefix.log", join('',
	"NOTE: 1 identical sequences will be ignored for subsequent analysis\n",
	"ERROR: Numerical underflow (lh-derivative). Run again with -safe option\n",
));

my $reason = '';
ok(!iqtreeOutputComplete($prefix, $alignment, \$reason),
	'an underflowing partial IQ-TREE run is not complete');
like($reason, qr/no successful completion signature/,
	'the failed completion reason is diagnostic');

write_file("$prefix.log", "Analysis results written to:\n");
ok(!iqtreeOutputComplete($prefix, $alignment, \$reason),
	'a successful log cannot hide omitted identical samples');
like($reason, qr/1 missing \(C\)/,
	'taxon validation identifies the missing alignment sample');

write_file("$prefix.treefile", "((A:0,B:0)99:0.1,C:0.2);\n");
ok(iqtreeOutputComplete($prefix, $alignment, \$reason),
	'a completed tree with every input taxon is accepted');
is($reason, '', 'successful validation clears the diagnostic reason');

write_file("$prefix.log", join('',
	"Analysis results written to:\n",
	"ERROR: late failure\n",
));
ok(!iqtreeOutputComplete($prefix, $alignment, \$reason),
	'an error after a prior completion signature invalidates the output');
like($reason, qr/ends in an error/,
	'a late failure has a distinct diagnostic');

write_file("$prefix.log", "Analysis results written to:\n");
for my $suffix (qw(.bionj .ckp.gz .mldist .uniqueseq.phy .varsites)) {
	write_file($prefix.$suffix, "temporary\n");
}
cleanupIQTreeTransients($prefix);
for my $suffix (qw(.bionj .ckp.gz .mldist .uniqueseq.phy .varsites)) {
	ok(!-e $prefix.$suffix, "$suffix is removed after successful completion");
}
ok(-s "$prefix.treefile", 'the validated tree is retained');
ok(-s "$prefix.log", 'the completed IQ-TREE log is retained');

my $fakeIqtree = File::Spec->catfile($tmp, 'fake-iqtree');
write_file($fakeIqtree, <<'FAKE_IQTREE');
#!/bin/sh
if [ "$1" = "--version" ]; then
	printf 'IQ-TREE multicore version 3.0.1\n'
	exit 0
fi
arguments="$*"
prefix=''
safe=0
while [ "$#" -gt 0 ]; do
	case "$1" in
		-pre)
			shift
			prefix="$1"
			;;
		-safe)
			safe=1
			;;
	esac
	shift
done
printf '%s\n' "$arguments" >> "${prefix}.calls"
printf 'temporary\n' > "${prefix}.uniqueseq.phy"
if [ "$safe" -eq 0 ]; then
	printf '(A:0,B:0);\n' > "${prefix}.treefile"
	printf 'ERROR: Numerical underflow (lh-derivative). Run again with -safe option\n' > "${prefix}.log"
	exit 2
fi
printf '((A:0,B:0):0,C:0);\n' > "${prefix}.treefile"
printf 'Analysis results written to:\n' > "${prefix}.log"
exit 0
FAKE_IQTREE
chmod 0755, $fakeIqtree or die "Cannot make $fakeIqtree executable: $!";

my $retryPrefix = File::Spec->catfile($tmp, 'IQtree_retry');
{
	no warnings 'redefine';
	local *Mods::phyloTools::getProgPaths = sub { return $fakeIqtree };
	runQItree({
		inMSA => $alignment,
		IQtreeout => $retryPrefix,
		ncore => 1,
		outgr => '',
		bootStrap => 0,
		useAA => 0,
		iqtreeFast => 0,
		autoModel => 0,
		partition => '',
		runSafe => 0,
		iqMemMB => 0,
		iqPathogen => 0,
		iqLegacy => 0,
		constraintTree => '',
	});
}
my @calls = grep { length } split /\n/, slurp("$retryPrefix.calls");
is(scalar(@calls), 2, 'numerical underflow triggers exactly one retry');
unlike($calls[0], qr/(?:^|\s)-keep-ident(?:\s|$)/,
	'the IQ-TREE command does not use -keep-ident');
unlike($calls[0], qr/(?:^|\s)-safe(?:\s|$)/,
	'the small initial test attempt uses the normal likelihood kernel');
like($calls[1], qr/(?:^|\s)-safe(?:\s|$)/,
	'the numerical-underflow retry uses the safe likelihood kernel');
ok(-s "$retryPrefix.unsafe.log", 'the failed unsafe diagnostic is preserved');
ok(iqtreeOutputComplete($retryPrefix, $alignment, \$reason),
	'the safe retry completes with every input taxon');
ok(!-e "$retryPrefix.uniqueseq.phy",
	'the successful retry removes its unique-sequence temporary file');

done_testing();
