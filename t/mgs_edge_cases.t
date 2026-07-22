use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use IPC::Open3 qw(open3);
use Symbol qw(gensym);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::Binning ();

sub write_file {
	my ($path, $contents) = @_;
	open my $fh, '>', $path or die "Cannot write $path: $!";
	print {$fh} $contents;
	close $fh or die "Cannot close $path: $!";
}

sub slurp {
	my ($path) = @_;
	open my $fh, '<', $path or die "Cannot read $path: $!";
	local $/;
	return <$fh>;
}

my $tmp = tempdir(CLEANUP => 1);

# A marker in the final comma/tab field used to retain its line ending, so it
# was silently ranked as an ordinary gene. Run the complete sorter to cover
# both the marker parsing and its final-MGS EOF flush.
my $gc = File::Spec->catdir($tmp, 'GC');
make_path($gc);
write_file(File::Spec->catfile($gc, 'FMG.subset.cats'), "COG1\tunused\t markerZ \r\n");
write_file(
	File::Spec->catfile($gc, 'compl.incompl.95.fna.clstr.idx'),
	"geneA\t>S1__geneA,>S2__geneA\nmarkerZ\t>S1__markerZ,>S2__markerZ\n",
);
my $mgs = File::Spec->catfile($tmp, 'groups.core');
write_file($mgs, "MGS1\tgeneA\t10\t0\t1\tunused\nMGS1\tmarkerZ\t10\t0\t1\tunused\n");
write_file(File::Spec->catfile($tmp, 'groups.obs'), "MGS1\t10\n");

my $sorter = File::Spec->catfile(
	$Bin, '..', 'secScripts', 'MGS', 'resortMGSgenes4importance.pl',
);
my $sorter_err = gensym;
my $sorter_pid = open3(
	undef, my $sorter_out, $sorter_err,
	$^X, '-I'.File::Spec->catdir($Bin, '..'),
	$sorter, $gc, $mgs, 'FMG', 'test', 95,
);
my $sorter_stdout = do { local $/; <$sorter_out> // '' };
my $sorter_stderr = do { local $/; <$sorter_err> // '' };
waitpid($sorter_pid, 0);
is($? >> 8, 0, 'gene-priority sorter accepts a CRLF-terminated final marker field')
	or diag($sorter_stdout, $sorter_stderr);
is(
	slurp("$mgs.srt"),
	"MGS1\tmarkerZ,geneA\n",
	'trimmed marker ID is prioritised and the final MGS is flushed at EOF',
);

# The per-family representative scan formerly flushed only at MGS transitions
# and manufactured an undefined key/value when every candidate failed QC.
my $guide = File::Spec->catfile($tmp, 'MAGvsGC.txt');
write_file($guide, join('',
	"Bin\tMGS\tCompleteness\tContamination\tRepresentative4MGS\tLCAcompleteness\tN50\n",
	"S2__binB\tMGS1\t90\t1\t\t90\t1000\n",
	"S1__binA\tMGS1\t90\t1\t*\t90\t1000\n",
	"S3__last\tMGS2\t80\t2\t*\t80\t900\n",
	"S4__bad\tMGS2\t60\t1\t\t60\t800\n",
));
my $mapping = {
	opt => { smpl_order => [qw(S1 S2 S3 S4)] },
	S1 => { FamGroup => 'FamA', AssGroup => '' },
	S2 => { FamGroup => 'FamA', AssGroup => '' },
	S3 => { FamGroup => 'FamB', AssGroup => '' },
	S4 => { FamGroup => 'FamBad', AssGroup => '' },
};
my ($representatives, $warnings);
{
	local $SIG{__WARN__} = sub { $warnings .= join('', @_); };
	$representatives = Mods::Binning::getRepresentBinsPerFamily($guide, $mapping);
}
is_deeply(
	$representatives,
	{
		'FamA.MGS1' => 'S1__binA',
		'FamB.MGS2' => 'S3__last',
	},
	'per-family representatives include the final MGS and deterministic tie breaking',
);
like(
	$warnings,
	qr/No eligible representative bin for family 'FamBad' in MGS 'MGS2'; skipping/,
	'a family with no QC-eligible candidate is reported and skipped',
);
unlike($warnings, qr/uninitialized/i, 'ineligible families do not trigger undefined-value warnings');

done_testing;
