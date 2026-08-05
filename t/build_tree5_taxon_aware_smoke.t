use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

my $root = File::Spec->rel2abs(File::Spec->catdir($Bin, '..'));
my $temporary = tempdir('buildtree-taxon-aware-XXXXXX', TMPDIR => 1, CLEANUP => 1);
my $output = File::Spec->catdir($temporary, 'output');
make_path($output);

sub write_file {
	my ($path, $contents) = @_;
	open my $handle, '>', $path or die "Cannot write $path: $!";
	print {$handle} $contents or die "Cannot populate $path: $!";
	close $handle or die "Cannot close $path: $!";
}

my $mafft = File::Spec->catfile($temporary, 'mafft-pass-through');
write_file($mafft, <<'SH');
#!/bin/sh
for argument do
	input="$argument"
done
exec /bin/cat "$input"
SH
chmod 0755, $mafft or die "Cannot make $mafft executable: $!";

my $trimal = File::Spec->catfile($temporary, 'trimal-backtranslate');
write_file($trimal, <<'SH');
#!/bin/sh
while [ "$#" -gt 0 ]; do
	case "$1" in
		-out) output="$2"; shift 2 ;;
		-backtrans) input="$2"; shift 2 ;;
		*) shift ;;
	esac
done
exec /bin/cp "$input" "$output"
SH
chmod 0755, $trimal or die "Cannot make $trimal executable: $!";

my $config = File::Spec->catfile($temporary, 'MATAFILERcfg.txt');
write_file($config, join("\n",
	"MFLRDir\t$root",
	"BINDir\t$root/bin",
	"DBDir\t$temporary",
	"MGSTKDir\t$temporary",
	"Rpath\tR",
	"SINGcmd\ttrue",
	"CONDcmd\ttrue",
	"CONDA\tshell hook",
	"CONDAbaseEnv\tbase",
	"PY3cmd\tpython3",
	"Rscript\tRscript",
	"pigz\t$root/t/bin/pigz",
	"mafft\t$mafft",
	"trimal\t$trimal",
)."\n");

my (%aa, %nt);
for my $gene (qw(g1 g2 g4)) {
	for my $sample (qw(s1 s2 s3 s4)) {
		$aa{"$sample|$gene"} = 'MKTAAAVVVQ';
		$nt{"$sample|$gene"} = 'ATGAAAACTGCTGCTGCTGTTGTTGTTCAG';
	}
}
for my $sample (qw(s1 s2 s5)) {
	$aa{"$sample|g3"} = 'MKTAAAVVVQ';
	$nt{"$sample|g3"} = 'ATGAAAACTGCTGCTGCTGTTGTTGTTCAG';
}
my $faa = File::Spec->catfile($temporary, 'input.faa');
my $fna = File::Spec->catfile($temporary, 'input.fna');
write_file($faa, join('', map { ">$_\n$aa{$_}\n" } sort keys %aa));
write_file($fna, join('', map { ">$_\n$nt{$_}\n" } sort keys %nt));
my $categories = File::Spec->catfile($temporary, 'input.cat');
write_file($categories, join('',
	map {
		my $gene = $_;
		join("\t", sort grep { /\|\Q$gene\E\z/ } keys %aa)."\n"
	} qw(g1 g2 g3 g4)
));

my $wrapper = File::Spec->catfile($temporary, 'run-buildtree.pl');
write_file($wrapper, <<'PERL');
use strict;
use warnings;
use Mods::IO_Tamoc_progs qw(setConfigFile);
my $config = shift @ARGV;
my $script = shift @ARGV;
setConfigFile($config);
my $result = do $script;
die $@ if $@;
die "Cannot execute $script: $!\n" unless defined $result;
PERL

my $script = File::Spec->catfile($root, 'secScripts', 'phylo', 'buildTree5.pl');
my @command = (
	$^X, '-I'.$root, $wrapper, $config, $script,
	'-fna', $fna, '-aa', $faa, '-cats', $categories,
	'-outD', $output, '-smplSep', '\\|', '-AAtree', 0,
	'-MSAprogram', 2, '-runLengthCheck', 0, '-postAlignmentLocusQC', 0,
	'-taxonAwareMaxLoci', 3,
	'-taxonAwareCoreLoci', 2, '-taxonAwareCandidateExtra', 1,
	'-taxonAwareMinSequenceNT', 9, '-taxonAwareTargetLoci', 2,
	'-taxonAwareTargetNT', 30, '-placementMinOverlap', 6,
);
is(system(@command), 0, 'taxon-aware buildTree smoke workflow completes');

my $candidateAudit = File::Spec->catfile(
	$output, 'phylo', 'taxon_aware_locus_candidates.tsv');
my $finalAudit = File::Spec->catfile(
	$output, 'phylo', 'taxon_aware_locus_selection.tsv');
my $sampleAudit = File::Spec->catfile(
	$output, 'phylo', 'taxon_aware_sample_selection.tsv');
ok(-s $candidateAudit, 'candidate locus audit is written');
ok(-s $finalAudit, 'final locus audit is written');
ok(-s $sampleAudit, 'final sample audit is written');

open my $candidateHandle, '<', $candidateAudit or die $!;
my $candidateText = do { local $/; <$candidateHandle> };
close $candidateHandle;
like($candidateText, qr/^candidate\tg4\t1\t4\tqc_backfill\t/m,
	'extra candidate is explicitly marked as QC backfill');

open my $finalHandle, '<', $finalAudit or die $!;
my $finalText = do { local $/; <$finalHandle> };
close $finalHandle;
like($finalText, qr/^final\tg3\t1\t3\ttaxon_rescue\t/m,
	'rare-sample locus is retained by the final taxon-rescue pass');
like($finalText, qr/^final\tg4\t0\t/m,
	'less complementary robust locus is left out of the bounded final set');

open my $sampleHandle, '<', $sampleAudit or die $!;
my $sampleText = do { local $/; <$sampleHandle> };
close $sampleHandle;
like($sampleText, qr/^s5\t1\t30\t1\t30\tplacement_candidate\tusable_sparse_anchor$/m,
	'rare but anchored sample is retained for placement');

my $mergedAlignment = File::Spec->catfile($output, 'MSA', 'MSAli.fna');
ok(-s $mergedAlignment, 'bounded final alignment is produced');
open my $alignmentHandle, '<', $mergedAlignment or die $!;
my $alignmentText = do { local $/; <$alignmentHandle> };
close $alignmentHandle;
like($alignmentText, qr/^>s5$/m, 'rescued sparse sample remains in the merged alignment');

done_testing();
