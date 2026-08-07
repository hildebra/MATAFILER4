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
		$nt{"$sample|$gene"} = $gene eq 'g1' ? ('ATG' x 10)
			: $gene eq 'g2' ? ('GCG' x 10) : ('AAA' x 10);
	}
}
for my $sample (qw(s1 s2 s5)) {
	$aa{"$sample|g3"} = 'MKTAAAVVVQ';
	$nt{"$sample|g3"} = $sample eq 's5' ? ('CTG'.('ATG' x 9)) : ('ATG' x 10);
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
	'-MSAprogram', 2, '-runLengthCheck', 0, '-postAlignmentLocusQC', 1,
	'-taxonAwareMaxLoci', 3,
	'-taxonAwareCoreLoci', 2, '-taxonAwareCandidateExtra', 1,
	'-taxonAwareMinSequenceNT', 9, '-taxonAwareTargetLoci', 2,
	'-taxonAwareTargetNT', 30, '-placementMinOverlap', 6,
	'-rateMergePartitions', 1, '-rateMergeMaxBins', 4, '-rateMergeTargetSites', 1,
	'-rateMergeMinLoci', 1, '-rateMergeMinSites', 1,
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

my $partitionFile = $mergedAlignment.'.partition.RAXML';
my $rateAudit = File::Spec->catfile($output, 'phylo', 'rate_merged_partitions.tsv');
ok(-s $partitionFile, 'deterministically grouped partition file is produced');
ok(-s $rateAudit, 'rate/GC partition assignments are audited');
open my $partitionHandle, '<', $partitionFile or die $!;
my @partitionLines = grep { /\S/ } <$partitionHandle>;
close $partitionHandle;
cmp_ok(scalar(@partitionLines), '>=', 2,
	'divergence and GC differences produce more than one partition bin');
cmp_ok(scalar(@partitionLines), '<=', 4,
	'scaled bin count respects the configured maximum');
open my $rateAuditHandle, '<', $rateAudit or die $!;
my $rateAuditText = do { local $/; <$rateAuditHandle> };
close $rateAuditHandle;
like($rateAuditText, qr/^locus\talignment\tselection_phase\tstart\tend\talignment_sites\teffective_called_sites\trate_proxy/m,
	'rate-partition audit has stable coordinate and metric columns');
like($rateAuditText, qr/\tp90_consensus_divergence\t/,
	'MSAfix post-alignment divergence supplies the rate proxy');
is(scalar(() = $rateAuditText =~ /^g[123]\t/gm), 3,
	'every selected locus has exactly one audited partition assignment');
like($rateAuditText, qr/^g3\t[^\t]+\ttaxon_rescue\t[^\n]+\ttaxon_rescue_to_/m,
	'taxon-rescue locus is attached to a robust bin and cannot seed its own partition');

my $collapsedOutput = File::Spec->catdir($temporary, 'collapsed-output');
make_path($collapsedOutput);
my @collapsedCommand = @command;
for my $index (0 .. $#collapsedCommand - 1) {
	$collapsedCommand[$index + 1] = $collapsedOutput
		if $collapsedCommand[$index] eq '-outD';
	$collapsedCommand[$index + 1] = 20
		if $collapsedCommand[$index] eq '-rateMergeMinLoci';
	$collapsedCommand[$index + 1] = 20_000
		if $collapsedCommand[$index] eq '-rateMergeMinSites';
	$collapsedCommand[$index + 1] = 1
		if $collapsedCommand[$index] eq '-postAlignmentLocusQC';
}
is(system(@collapsedCommand), 0,
	'rate-merging smoke workflow completes with production small-bin thresholds');
my $collapsedPartition = File::Spec->catfile(
	$collapsedOutput, 'MSA', 'MSAli.fna.partition.RAXML');
open my $collapsedHandle, '<', $collapsedPartition or die $!;
my @collapsedLines = grep { /\S/ } <$collapsedHandle>;
close $collapsedHandle;
is(scalar(@collapsedLines), 1,
	'undersized rate/GC bins collapse into a single pooled partition');
like($collapsedLines[0], qr/= \d+-\d+(?:, \d+-\d+){2}$/,
	'collapsed IQ-TREE partition retains all three non-contiguous locus ranges');
my $collapsedAudit = File::Spec->catfile(
	$collapsedOutput, 'phylo', 'rate_merged_partitions.tsv');
open my $collapsedAuditHandle, '<', $collapsedAudit or die $!;
my $collapsedAuditText = do { local $/; <$collapsedAuditHandle> };
close $collapsedAuditHandle;
like($collapsedAuditText, qr/\tp90_consensus_divergence\t/,
	'native post-alignment P90 consensus divergence is the preferred rate proxy');

done_testing();
