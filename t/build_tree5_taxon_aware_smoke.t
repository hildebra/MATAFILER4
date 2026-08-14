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

sub slurp {
	my ($path) = @_;
	open my $handle, '<', $path or die "Cannot read $path: $!";
	local $/;
	return <$handle>;
}

my $mafft = File::Spec->catfile($temporary, 'mafft-pass-through');
write_file($mafft, <<'SH');
#!/bin/sh
if [ -n "$MATAFILER_TEST_MAFFT_COUNT" ]; then
	printf "%s\n" "$*" >> "$MATAFILER_TEST_MAFFT_COUNT"
fi
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
my $msaFixShim = File::Spec->catfile($temporary, 'MSAfix-v2.15-shim');
write_file($msaFixShim, <<'PERL');
#!/usr/bin/env perl
use strict;
use warnings;

my (@forward, $report, $threads, $singleAlignmentMode, %recovery);
while (@ARGV) {
    my $argument = shift @ARGV;
    $singleAlignmentMode = 1 if $argument eq '-i';
    if ($argument eq '-minOverlapMSA') {
        shift @ARGV;
        next;
    }
    if ($argument eq '-threads') {
        $threads = shift @ARGV;
        next;
    }
    if ($argument eq '-recoverTechnicalOffsets') {
        $recovery{enabled} = 1;
        next;
    }
    if ($argument =~ /^-(codingFrame|geneticCode|recoveryBand)$/) {
        $recovery{$1} = shift @ARGV;
        next;
    }
    $report = $ARGV[0] if $argument eq '-report' && @ARGV;
    push @forward, $argument;
}
die "MSAfix shim requires a positive -threads value\n"
    unless defined($threads) && $threads =~ /^\d+$/ && $threads > 0;
if ($singleAlignmentMode) {
    die "MSAfix shim requires coding-NT recovery\n" unless $recovery{enabled};
    die "MSAfix shim requires frame 1, genetic code 11, and recovery band 3\n"
        unless ($recovery{codingFrame} // '') eq '1'
            && ($recovery{geneticCode} // '') eq '11'
            && ($recovery{recoveryBand} // '') eq '3';
}

my $real = $ENV{MATAFILER_TEST_MSAFIX_REAL}
	or die "MATAFILER_TEST_MSAFIX_REAL is unset\n";
system($real, @forward);
exit($? >> 8) if $?;
exit 0 unless defined $report;
die "MSAfix shim report is missing or empty\n" unless -s $report;
open my $input, '<', $report or die "Cannot read $report: $!\n";
my $header = <$input> // die "Empty MSAfix report\n";
chomp $header;
exit 0 if $header =~ /(?:^|\t)effective_sites(?:\t|$)/;
my @columns = split /\t/, $header, -1;
my %index = map { $columns[$_] => $_ } 0 .. $#columns;
die "Legacy MSAfix report has no alignment column\n" unless exists $index{alignment};
my @rows = <$input>;
close $input or die "Cannot close $report: $!\n";
open my $output, '>', $report or die "Cannot update $report: $!\n";
print {$output} $header,
	"\tcalled_cells\tgc_cells\tgc_fraction\teffective_sites\n";
for my $row (@rows) {
	chomp $row;
	my @field = split /\t/, $row, -1;
	my $alignment = $field[$index{alignment}];
	open my $alignmentFH, '<', $alignment
		or die "Cannot read alignment $alignment: $!\n";
	my ($sequences, $called, $gc) = (0, 0, 0);
	while (my $line = <$alignmentFH>) {
		if ($line =~ /^>/) {
			$sequences++;
			next;
		}
		$line = uc($line);
		$called += ($line =~ tr/ACGT//);
		$gc += ($line =~ tr/GC//);
	}
	close $alignmentFH or die "Cannot close $alignment: $!\n";
	my $gcFraction = $called ? $gc / $called : 0;
	my $effective = $sequences ? $called / $sequences : 0;
	print {$output} join("\t", @field, $called, $gc,
		sprintf('%.8g', $gcFraction), sprintf('%.8g', $effective)), "\n";
}
close $output or die "Cannot close updated $report: $!\n";
PERL
chmod 0755, $msaFixShim or die "Cannot make $msaFixShim executable: $!";
local $ENV{MATAFILER_TEST_MSAFIX_REAL} = File::Spec->catfile($root, 'bin', 'MSAfix');

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
	"MSAfix\t$msaFixShim",
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
my $mafftCount = File::Spec->catfile($temporary, 'mafft.calls');
local $ENV{MATAFILER_TEST_MAFFT_COUNT} = $mafftCount;
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

my $terminalOutput = File::Spec->catdir($temporary, 'terminal-output');
my $terminalCategories = File::Spec->catfile($temporary, 'terminal.cat');
write_file($terminalCategories, "s1|g1\ts2|g1\n");
my @terminalCommand = @command;
for my $index (0 .. $#terminalCommand - 1) {
	$terminalCommand[$index + 1] = $terminalOutput
		if $terminalCommand[$index] eq '-outD';
	$terminalCommand[$index + 1] = $terminalCategories
		if $terminalCommand[$index] eq '-cats';
}
is(system(@terminalCommand), 0,
	'a category set without any three-sample locus exits as a successful terminal outcome');
my $terminalMarker = File::Spec->catfile($terminalOutput, 'noTree.sto');
ok(-s $terminalMarker, 'the taxon-aware terminal outcome is persisted');
open my $terminalHandle, '<', $terminalMarker or die $!;
my $terminalText = do { local $/; <$terminalHandle> };
close $terminalHandle;
like($terminalText,
	qr/^status\tvalid_no_tree\nreason\ttaxon_aware_no_category_with_three_usable_samples\n/m,
	'the terminal marker records the exact stable selection reason');
my $terminalState = File::Spec->catfile($terminalOutput, 'buildTree.state.tsv');
ok(-s $terminalState,
	'the valid no-tree outcome retains the consolidated workflow state');
unlike(slurp($terminalState), qr/^status\tfailed$/m,
	'the valid no-tree outcome is not recorded as a workflow failure');

my $candidateAudit = File::Spec->catfile(
	$output, 'phylo', 'taxon_aware_locus_candidates.tsv');
my $finalAudit = File::Spec->catfile(
	$output, 'phylo', 'taxon_aware_locus_selection.tsv');
my $sampleAudit = File::Spec->catfile(
	$output, 'phylo', 'taxon_aware_sample_selection.tsv');
my $attritionAudit = File::Spec->catfile(
	$output, 'phylo', 'selection_attrition.tsv');
ok(-s $candidateAudit, 'candidate locus audit is written');
ok(-s $finalAudit, 'final locus audit is written');
ok(-s $sampleAudit, 'final sample audit is written');
ok(-s $attritionAudit, 'compact end-to-end selection attrition audit is written');

open my $candidateHandle, '<', $candidateAudit or die $!;
my $candidateText = do { local $/; <$candidateHandle> };
close $candidateHandle;
like($candidateText, qr/^candidate\tg4\t1\t4\tqc_backfill\t/m,
	'extra candidate is explicitly marked as QC backfill');

open my $finalHandle, '<', $finalAudit or die $!;
my $finalText = do { local $/; <$finalHandle> };
close $finalHandle;
like($finalText,
	qr/^stage\tgene\tselected\trank\tphase\tquality_score\trobust_score\toccupancy\tprevalence\tinformation_score\tinformation_density\tvariable_density\texcess_variation_penalty/m,
	'final locus audit exposes prevalence, bounded information, and excess-variation scoring');
like($finalText, qr/^final\tg3\t1\t3\ttaxon_rescue\t/m,
	'rare-sample locus is retained by the final taxon-rescue pass');
like($finalText, qr/^final\tg4\t0\t/m,
	'less complementary robust locus is left out of the bounded final set');

open my $sampleHandle, '<', $sampleAudit or die $!;
my $sampleText = do { local $/; <$sampleHandle> };
close $sampleHandle;
like($sampleText, qr/^s5\t1\t30\t1\t30\tplacement_candidate\tusable_sparse_anchor$/m,
	'rare but anchored sample is retained for placement');

open my $attritionHandle, '<', $attritionAudit or die $!;
my $attritionText = do { local $/; <$attritionHandle> };
close $attritionHandle;
like($attritionText, qr/^input_loci\t4$/m, 'attrition audit records all input loci');
like($attritionText, qr/^candidate_loci\t4$/m, 'attrition audit records the bounded alignment candidates');
like($attritionText, qr/^final_loci\t3$/m, 'attrition audit records the bounded final locus set');
like($attritionText, qr/^backbone_samples\t5$/m, 'attrition audit records final tree samples');

my $mergedAlignment = File::Spec->catfile($output, 'MSA', 'MSAli.fna');
my $compressedMergedAlignment = "$mergedAlignment.gz";
ok(-s $compressedMergedAlignment, 'bounded final alignment is retained in compressed form');
open my $alignmentHandle, '-|', 'gzip', '-cd', $compressedMergedAlignment or die $!;
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


my $initialMafftRuns = (() = slurp($mafftCount) =~ /^/gm);
cmp_ok($initialMafftRuns, q{>}, 0,
	q{the initial workflow ran per-locus MSA jobs});
my @downstreamResume = (@command, q{-continue}, 1, q{-iqLegacy}, 1);
is(system(@downstreamResume), 0,
	q{a downstream-only option change resumes from the retained selected MSA});
is((() = slurp($mafftCount) =~ /^/gm), $initialMafftRuns,
	q{a downstream-only option change does not run the per-locus MSAs again});
my $validCheckpointMafftRuns = (() = slurp($mafftCount) =~ /^/gm);
write_file($mergedAlignment, ">corrupt1\nACGT!\n>corrupt2\nACGT!\n");
my @invalidCheckpointResume = (@command, q{-continue}, 1, q{-iqLegacy}, 1);
is(system(@invalidCheckpointResume), 0,
	q{an invalid retained alignment checkpoint is rebuilt before tree inference});
cmp_ok((() = slurp($mafftCount) =~ /^/gm), q{>}, $validCheckpointMafftRuns,
	q{an invalid retained alignment checkpoint reruns the per-locus MSAs});
ok(!-e $mergedAlignment && -s $compressedMergedAlignment,
	q{the rebuilt retained alignment replaces the corrupt plain checkpoint});
my $fullAlignment = File::Spec->catfile($output, q{MSA}, q{MSAli.full.fna});
my $compressedFullAlignment = "$fullAlignment.gz";
write_file($fullAlignment, ">corrupt1\nACGT!\n>corrupt2\nACGT!\n");
my $strictCheckpointMafftRuns = (() = slurp($mafftCount) =~ /^/gm);
my @invalidFullCheckpointResume = (
	@command, q{-continue}, 1, q{-iqLegacy}, 1, q{-strictBackbone}, 1);
is(system(@invalidFullCheckpointResume), 0,
	q{an invalid retained full alignment is replaced before strict-backbone inference});
is((() = slurp($mafftCount) =~ /^/gm), $strictCheckpointMafftRuns,
	q{replacing an invalid full alignment reuses the valid primary checkpoint});
ok(!-e $fullAlignment && -s $compressedFullAlignment,
	q{the corrected full alignment is retained in compressed form});
open my $fullAlignmentHandle, q{-|}, q{gzip}, q{-cd}, $compressedFullAlignment or die $!;
my $fullAlignmentText = do { local $/; <$fullAlignmentHandle> };
close $fullAlignmentHandle;
unlike($fullAlignmentText, qr/corrupt1/,
	q{strict-backbone does not propagate the corrupt full alignment to inference});
my $workflowState = File::Spec->catfile($output, q{buildTree.state.tsv});
ok(-s $workflowState,
	q{one consolidated workflow state records checkpoint and lifecycle data});
my $workflowStateText = slurp($workflowState);
like($workflowStateText, qr/^msa_selection_policy\t.*minimum_overlap=/m,
	q{the state records the MSA-selection policy});
like($workflowStateText, qr/^tree_stage_policy\t.*iqtree_legacy=1/m,
	q{the state records the downstream tree-stage policy});
ok(!-e File::Spec->catfile($output, q{MSA}, q{alignment_work.policy.tsv})
	&& !-e File::Spec->catfile($output, q{phylo}, q{post_alignment.policy.tsv})
	&& !-e File::Spec->catfile($output, q{phylo}, q{post_alignment_locus_qc.policy.tsv}),
	q{the consolidated state replaces the three legacy policy files});
my @selectionResume = (@command, q{-continue}, 1, q{-minOverlapMSA}, 0.1);
is(system(@selectionResume), 0,
	q{an MSA-selection option change rebuilds the selected alignment});
cmp_ok((() = slurp($mafftCount) =~ /^/gm), q{>}, $initialMafftRuns,
	q{an MSA-selection option change runs the per-locus MSAs again});


done_testing();
