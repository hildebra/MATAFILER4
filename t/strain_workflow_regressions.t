use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::GenoMetaAss qw(readFasta);

sub write_file {
	my ($path, $contents) = @_;
	open my $fh, '>', $path or die "Cannot write $path: $!";
	print {$fh} $contents or die "Cannot write $path: $!";
	close $fh or die "Cannot close $path: $!";
}

sub slurp {
	my ($path) = @_;
	open my $fh, '<', $path or die "Cannot read $path: $!";
	local $/;
	return <$fh>;
}

my $tmp = tempdir(CLEANUP => 1);
my $fasta = File::Spec->catfile($tmp, 'records.fa');
write_file($fasta, <<'FASTA');
>keep1 D=8 CSP=0.01
AAAA
>drop D=2 CSP=0.50
CCCC
>keep3 D=5 CSP=0.02
GGGG
FASTA

my %wanted = (keep1 => 1, keep3 => 1);
is_deeply(
	readFasta($fasta, 1, '\\s', \%wanted),
	{ keep1 => 'AAAA', keep3 => 'GGGG' },
	'FASTA subset selection applies independently to intermediate and final records',
);
is_deeply(
	readFasta($fasta, 0, '\\s', \%wanted),
	{
		'keep1 D=8 CSP=0.01' => 'AAAA',
		'keep3 D=5 CSP=0.02' => 'GGGG',
	},
	'FASTA subset lookup can use short IDs while retaining full headers',
);

my $glob_dir = File::Spec->catdir($tmp, 'glob');
mkdir $glob_dir or die "Cannot create $glob_dir: $!";
write_file(File::Spec->catfile($glob_dir, 'a_empty.fa'), '');
write_file(File::Spec->catfile($glob_dir, 'b_records.fa'), ">later\nACGT\n");
is_deeply(
	readFasta(File::Spec->catfile($glob_dir, '*.fa'), 1, '\\s'),
	{ later => 'ACGT' },
	'an empty member of a FASTA glob does not suppress later files',
);

my $strain = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'MGS', 'strain_within.pl'));
like($strain, qr/sub consensusInputState .*?\$nt_ready && \$aa_ready.*?return 'regenerate' if \$vcf_ready/s,
	'consensus resume requires the paired NT and AA outputs and repairs from VCF');
like($strain, qr/readFasta\(\$fastaf,1,"\\\\s",\\%subG\).*?readFasta\(\$fastafAA,0,"\\\\s",\\%subG\)/s,
	'within-strain extraction reads only candidate consensus genes');
like($strain, qr/test -s "?\.shellQuote\(\$IQtreef\).*?touch "?\.shellQuote\(\$treeStone\)/s,
	'a tree completion stone is conditional on a nonempty tree');
like($strain, qr/unlink \$treeStone.*?if \$doSubmit && -e \$treeStone/s,
	'a submitted tree retry cannot pass through a stale completion stone');
like($strain, qr/unlink \$IQtreef.*?stale tree output/s,
	'a submitted tree retry must publish a fresh nonempty tree');
like($strain, qr/unlink \$checkF.*?stale split-worker checkpoint.*?if \$doSubmit && -e \$checkF/s,
	'a submitted split-worker retry clears its stale completion stone');
like($strain, qr/ConspecificMGS\.\$subJob\.log.*?sub mergeConspecificLogs/s,
	'split workers write isolated conspecific logs that are explicitly merged');
like($strain, qr/qsubSystemJobAlive\([^\n]+QSBoptHR[^\n]+if [^\n]+doSubmit/,
	'dry runs do not poll scheduler jobs that were never submitted');
like($strain, qr/\$nxtCmd \.= "-submit \$doSubmit ";.*?-qsubSystem/s,
	'postprocessing inherits submission state and the selected queue backend');
like($strain, qr/sub assertSafeWorkflowRemoval .*?Refusing to remove filesystem root.*?matafiler-strain-workdir/s,
	'custom recursive output removal requires a workflow-owned directory');

my $build_tree = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'phylo', 'buildTree5.pl'));
like($build_tree, qr/if \(\$numSeq < 3\)/,
	'three-sample MGS accepted by the wrapper are retained for a minimal tree');

done_testing();
