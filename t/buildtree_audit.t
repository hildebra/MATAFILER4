use strict;
use warnings;
no warnings 'once';
use FindBin qw($Bin);
use lib "$Bin/..", "$Bin/lib";
use MFTestConfig;
use File::Path qw(make_path);
use File::Basename qw(dirname);
use File::Temp qw(tempdir);
use File::Spec;
use Text::ParseWords qw(shellwords);
use Test::More;
use Mods::phyloTools qw(MSA runFasttree runVeryFasttree getTreeLeafs);
use Mods::StrainPlacement qw(read_epa_jplace write_epa_placed_tree);
my $root = File::Spec->rel2abs("$Bin/..");
my $tmp = tempdir('buildtree-audit-XXXXXX', TMPDIR => 1, CLEANUP => 1);
sub put {
	my ($path, $text) = @_;
	make_path(dirname($path));
	open my $out, '>', $path or die "$path: $!";
	print {$out} $text or die $!;
	close $out or die $!;
}
sub slurp { open my $in, '<', $_[0] or die "$! $_[0]"; local $/; return <$in>; }
sub run {
	my ($log, @cmd) = @_;
	my $pid = fork(); die $! unless defined $pid;
	if (!$pid) {
		open STDOUT, '>', $log or die $!;
		open STDERR, '>', "$log.err" or die $!;
		exec @cmd; die "Cannot execute $cmd[0]: $!";
	}
	waitpid($pid, 0); return $?;
}
sub sq { my $s = shift; $s =~ s/'/'"'"'/g; return "'$s'"; }
my $fake = "$tmp/fasttree";
put($fake, <<'FAKE');
use strict; use warnings;
if (($ENV{TREE_TEST_MODE} // '') eq 'fail') { print '(partial:'; exit 7; }
exit 0 if ($ENV{TREE_TEST_MODE} // '') eq 'empty';
my $input = $ARGV[-1];
open my $in, '<', $input or die "$input: $!";
my @names; while (<$in>) { push @names, $1 if /^>(\S+)/; }
if ($ENV{TREE_TEST_COUNT}) {
	open my $count, '>>', $ENV{TREE_TEST_COUNT} or die $!;
	print {$count} "run\n"; close $count;
}
print '(', join(',', map { "$_:0.1" } @names), ");\n";
FAKE
my $fasta = "$tmp/input with 'quotes.fna";
put($fasta, ">A\nACGT\n>B\nACGT\n>C\nACGT\n");
subtest 'Shared FastTree wrapper publishes only successful nonempty output' => sub {
	no warnings 'redefine';
	local *Mods::phyloTools::getProgPaths = sub { return "env TREE_WRAPPER=1 $^X ".sq($fake); };
	for my $runner (\&runFasttree, \&runVeryFasttree) {
		my $tree = "$tmp/tree with 'quotes.nwk";
		put($tree, "(old:0.1,tree:0.1);\n");
		for my $mode ('fail', 'empty') {
			local $ENV{TREE_TEST_MODE} = $mode;
			my $ok = eval { $runner->($fasta, $tree, 0, 2); 1 };
			ok(!$ok, "$mode engine output is rejected");
			is(slurp($tree), "(old:0.1,tree:0.1);\n", 'prior valid tree is preserved');
		}
		$runner->($fasta, $tree, 0, 2);
		is(slurp($tree), "(A:0.1,B:0.1,C:0.1);\n", 'paths and environment wrapper survive atomic publication');
	}
	is(scalar(glob("$tmp/.fasttree-*")), undef, 'failed runs leave no temporary checkpoint');
	for my $mode (0, 1, 2, 4, 5) {
		my $command = MSA($fasta, "$tmp/aligned ' file.fna", 2, $mode, 400);
		$command =~ s/;\z//;
		my @args = shellwords($command);
		is(scalar(grep { $_ eq $fasta } @args), 1, "MSA mode $mode preserves the input as one argument");
		ok(grep({ $_ eq "$tmp/aligned ' file.fna" } @args), "MSA mode $mode preserves the output argument");
	}
};
subtest 'Tree leaf lookup reuses the tip parser without internal support labels' => sub {
	my $tree = "$tmp/labeled tree.nwk";
	put($tree, "(('sample A':0.1,'sample''B':0.2)99:0.3,C:0.4)root;\n");
	is_deeply([sort keys %{getTreeLeafs($tree)}], ["C", 'sample A', "sample'B"],
		'quoted tips are retained; root names and support values are not taxa');
};
subtest 'PHYLIP conversion preserves aligned sequences and validates before output' => sub {
	my $converter = "$root/secScripts/phylo/fasta2phylip.pl";
	my $input = "$tmp/converter.fna";
	my $log = "$tmp/converter.out";
	put($input, ">a description\r\nAC\r\nGT\r\n>long_sample_name\r\nTGCA\r\n");
	is(run($log, $^X, $converter, '-c', 50, $input), 0, 'wrapped CRLF FASTA converts');
	is(slurp($log), "2   4\na ACGT\nlong_sample_name TGCA\n", 'no extra padding and header matches all rows');
	put($input, ">only\nACGTAA\n");
	is(run($log, $^X, $converter, $input), 0, 'one-sequence input converts');
	like(slurp($log), qr/^1   6\n/, 'the final sequence contributes to the length');
	put($input, ">a\nACGT\n>b\nACGTA\n");
	ok(run($log, $^X, $converter, $input), 'unequal alignment lengths are rejected');
	is(-s $log, 0, 'invalid alignment does not produce a partial PHYLIP file');
	put($input, '>'.('a' x 50)."1\nACGT\n>".('a' x 50)."2\nACGT\n");
	ok(run($log, $^X, $converter, '-c', 50, $input), 'identifier truncation collisions are rejected');
};

my $source = slurp("$root/secScripts/phylo/buildTree5.pl");
my @helpers;
for my $name (qw(coreHyPhy selecAnalysis geneFileStem shellQuote finalizeStagedStrainCategory
	writeEpaPlacementReport publishEpaPlacement readEpaFilterBackboneTree
	writeEpaBackboneGraftAudit printEpaBackboneGraftSummary epaFilterMetricValue
	writeEpaPlacementFilterSummary printEpaPlacementFilterSummary)) {
	my ($sub) = $source =~ /(sub \Q$name\E\s*\{.*?^\})/ms;
	BAIL_OUT("Cannot extract $name") unless $sub;
	push @helpers, $sub;
}
my $setup = <<'SETUP';
package BuildTreeAuditFixture;
use File::Spec;
use File::Path qw(make_path remove_tree);
use Mods::GenoMetaAss qw(readFasta gzipopen fileGZs);
use Mods::phyloTools qw(getTreeLeafs);
use Mods::StrainPlacement qw(map_epa_placements_to_backbone filter_epa_placement_outliers write_epa_placed_tree);
use Mods::WorkflowResilience qw(retry_open retry_close retry_rename retry_unlink atomic_write_text);
our ($outgroup, $MsaWorkD, $MSAsubsD, $subsetPopgenStats, $reparseHyphyJson,
	$epaPendantOutlierFactor, $epaPendantMinThreshold, $parser);
our (@runs, @prepared, @pruned);
sub getProgPaths { return $parser; }
sub parseSeqId { return ($_[0]); }
sub limitedWarn { }
sub pruneTree { @pruned = @{$_[1]}; main::put($_[2], '(A,B,C);'); }
sub hyphy {
	push @runs, $_[2];
	@prepared = (main::slurp($_[0]));
	main::put($_[4], 'done'); main::put("$_[4].json.gz", 'fixture');
}
sub safeRemoveTree { remove_tree($_[0]); }
SETUP
my $loaded = eval $setup . join("\n", @helpers) . "\n1;";
BAIL_OUT($@) unless $loaded;
subtest 'Selection analysis retains ingroup samples and aligned codon lengths' => sub {
	my $msa = "$tmp/selection";
	my $tree = "$tmp/selection.tree";
	my $work = "$tmp/selection-work"; make_path($work);
	put($tree, '(A:0.1,B:0.1,C:0.1,OG:0.1,OGplus:0.1);');
	put("$msa/g.0.fna", ">A\nATGTAA\n>B\nATGTGG\n>C\nATGTGA\n>OG\nATGCCC\n>OGplus\nATGTAG\n");
	$BuildTreeAuditFixture::outgroup = '';
	BuildTreeAuditFixture::coreHyPhy($msa, 'g', '', $tree, $work, "$work/emptyOG.log");
	is(scalar(@BuildTreeAuditFixture::pruned), 5, 'an empty outgroup excludes no samples');
	like($BuildTreeAuditFixture::prepared[0], qr/>B\nATGTGG\n/, 'terminal tryptophan is not a stop codon');
	like($BuildTreeAuditFixture::prepared[0], qr/>A\nATGNNN\n/, 'terminal stop is masked without shortening the alignment');
	$BuildTreeAuditFixture::outgroup = 'OG';
	BuildTreeAuditFixture::coreHyPhy($msa, 'g', '', $tree, $work, "$work/exactOG.log");
	ok(!grep({ $_ eq 'OG' } @BuildTreeAuditFixture::pruned), 'the exact outgroup is excluded');
	ok(grep({ $_ eq 'OGplus' } @BuildTreeAuditFixture::pruned), 'similar sample names are retained');
	put("$msa/small.0.fna", ">A\nATGCCC\n>B\nATGCCC\n");
	my $afterReturn = 0;
	for (1) {
		BuildTreeAuditFixture::coreHyPhy($msa, 'small', '', $tree, $work, "$work/small.log");
		$afterReturn++;
	}
	is($afterReturn, 1, 'an underpowered locus returns normally rather than jumping out of the caller loop');
};
subtest 'Shared selection summaries keep column alignment and preserve MSA workspace' => sub {
	my $work = "$tmp/summary-work";
	put("$work/MSA/g.0.fna", ">A\nATGCCC\n>B\nATGCCC\n>C\nATGCCC\n");
	$BuildTreeAuditFixture::MsaWorkD = "$work/MSA";
	$BuildTreeAuditFixture::MSAsubsD = "$work/MSA/clnd";
	$BuildTreeAuditFixture::subsetPopgenStats = '';
	$BuildTreeAuditFixture::reparseHyphyJson = 0;
	my $parser = "$tmp/fubar-parser.pl";
	put($parser, 'print "\t3\t2\t0\t0\t1\t1\t1\t1\n";');
	$BuildTreeAuditFixture::parser = "$^X ".sq($parser);
	make_path("$tmp/selection-summary");
	BuildTreeAuditFixture::selecAnalysis(['g'], "$tmp/selection.tree", "$tmp/selection-summary", $work);
	like(slurp("$tmp/selection-summary/hyphy.fubar.txt"), qr/^g\t3\t2\t0\t0\t1\t1\t1\t1$/m,
		'legacy parser leading tab does not shift summary columns');
	ok(-s "$work/MSA/g.0.fna", 'selection cleanup retains the caller alignment workspace');
	ok(!-d "$work/hyphy", 'only the analysis temporary directory is removed');
};
subtest 'Legacy staged categories reject conflicts like the shard finalizer' => sub {
	my $raw = "$tmp/category.raw"; my $overlay = "$tmp/category.overlay";
	put($raw, "MGS\tL1\tA\tA|g1\nMGS\tL1\tA\tA|other\n");
	my $ok = eval { BuildTreeAuditFixture::finalizeStagedStrainCategory($raw, $overlay, "$tmp/category.final", 'MGS'); 1 };
	ok(!$ok && $@ =~ /Conflicting category/, 'contradictory identifiers cannot be silently overwritten');
	put($raw, "MGS\tL1\tA\tA|g1\n"); put($overlay, "absent\tOG\tOG|g1\n");
	$ok = eval { BuildTreeAuditFixture::finalizeStagedStrainCategory($raw, $overlay, "$tmp/category.final", 'MGS'); 1 };
	ok(!$ok && $@ =~ /absent locus/, 'overlay cannot introduce an ingroup-free locus');
};
subtest 'EPA publication preserves reports and rejects foreign or duplicate queries' => sub {
	my $jplace = "$tmp/test.jplace";
	my $json = '{"tree":"(A:0.1{0},B:0.1{1},C:0.1{2});","placements":[{"p":[[0,-10,1,0.05,0.01]],"n":["Q"]}]}';
	put($jplace, $json);
	my $ok = eval { read_epa_jplace($jplace, ['different']); 1 };
	ok(!$ok && $@ =~ /Unexpected EPA query/, 'retained placement for a different query set is rejected');
	(my $duplicate = $json) =~ s/"Q"/"Q","Q"/;
	put($jplace, $duplicate);
	$ok = eval { read_epa_jplace($jplace, ['Q']); 1 };
	ok(!$ok && $@ =~ /Duplicate EPA query/, 'duplicate query entries cannot overwrite each other');
	put($jplace, $json);
	my $epa = read_epa_jplace($jplace, ['Q']);
	my $backbone = "$tmp/publish/backbone.treefile"; my $tree = "$tmp/publish/primary.treefile";
	put($backbone, '(A:0.1,B:0.1,C:0.1);');
	$BuildTreeAuditFixture::epaPendantOutlierFactor = 10;
	$BuildTreeAuditFixture::epaPendantMinThreshold = 0.1;
	my $split = { backbone_overlap => { Q => {backbone_overlap_nt => 120} }, reason => {Q => 'sparse'} };
	BuildTreeAuditFixture::publishEpaPlacement($epa, $backbone, $tree, $split, "$tmp/publish");
	ok(getTreeLeafs($tree)->{Q}, 'shared publication grafts the query onto the saved backbone');
	my @rows = split /\n/, slurp("$tmp/publish/strict_backbone.epa_placements.tsv");
	is(scalar(split /\t/, $rows[0], -1), 16, 'placement report retains its 16 columns');
	like($rows[1], qr/^Q\tplaced\t120\tNA\tNA\t.*\tsparse$/, 'original overlap and reason metadata are preserved');
	for my $distance (-1e-9, 0.1 + 1e-9) {
		local $epa->{placements}{Q}{backbone_distal_length} = $distance;
		local $SIG{ALRM} = sub { die "grafting stalled\n"; };
		alarm 2;
		$ok = eval { write_epa_placed_tree(slurp($backbone), $tree, $epa->{placements}); 1 };
		alarm 0;
		ok($ok, "tolerated endpoint overshoot $distance finishes") or diag($@);
		ok(getTreeLeafs($tree)->{Q}, 'query is retained at the clamped endpoint');
	}
	my %collision = (A => { %{$epa->{placements}{Q}} });
	$ok = eval { write_epa_placed_tree(slurp($backbone), $tree, \%collision); 1 };
	ok(!$ok && $@ =~ /already a backbone tip/, 'query cannot duplicate a reference tip');
	{
		local $epa->{placements}{Q}{pendant_length} = -0.1;
		$ok = eval { write_epa_placed_tree(slurp($backbone), $tree, $epa->{placements}); 1 };
		ok(!$ok && $@ =~ /negative pendant/, 'negative pendant lengths are rejected');
	}
	my $old = slurp($tree);
	{
		no warnings 'redefine';
		local *BuildTreeAuditFixture::write_epa_placed_tree = sub { die "publication failure\n"; };
		$ok = eval { BuildTreeAuditFixture::publishEpaPlacement($epa, $backbone, $tree, $split, "$tmp/publish"); 1 };
	}
	ok(!$ok, 'publication failure propagates');
	is(slurp($tree), $old, 'failed filter redo preserves the prior published tree');
};

subtest 'Single-locus continuation respects input and downstream policy changes' => sub {
	my $config = "$tmp/config.txt";
	put($config, slurp("$Bin/MATAFILERcfg.txt")."\nfasttree\t$^X ".sq($fake)."\n");
	my $wrapper = "$tmp/run-buildtree.pl";
	put($wrapper, 'use Mods::IO_Tamoc_progs qw(setConfigFile); setConfigFile(shift @ARGV); my $s = shift @ARGV; do $s; die $@ if $@;');
	my $input = "$tmp/single.fna"; my $output = "$tmp/single output";
	put($input, join('', map { ">$_\n".('ACGT' x 30)."\n" } qw(A B C)));
	my $originalInput = slurp($input);
	local $ENV{TREE_TEST_COUNT} = "$tmp/single.calls";
	my @cmd = ($^X, "-I$root", $wrapper, $config, "$root/secScripts/phylo/buildTree5.pl",
		'-fna', $input, '-outD', $output, '-isAligned', 1, '-runFastTree', 1,
		'-runLengthCheck', 0, '-postAlignmentLocusQC', 0, '-postAlignmentSequenceOutlierMask', 0,
		'-taxonAwareLocusSelection', 0, '-rateMergePartitions', 0,
		'-completionMarker', "$output/treeDone.sto");
	my $log = "$tmp/single.log";
	my $status = run($log, @cmd);
	is($status, 0, 'initial single-locus tree succeeds') or diag(slurp("$log.err"));
	return unless $status == 0;
	is(slurp($input), $originalInput, 'input alignment is not rewritten through its staging symlink');
	unless (-s $ENV{TREE_TEST_COUNT}) { diag(slurp($log), slurp("$log.err")); fail('inference ran'); return; }
	is(run($log, @cmd, '-continue', 1), 0, 'unchanged single-locus output remains reusable') or diag(slurp("$log.err"));
	is(scalar(() = slurp($ENV{TREE_TEST_COUNT}) =~ /run/g), 1, 'unchanged resume does not rerun inference');
	like(slurp($log), qr/durable completion marker and current policy match/, 'unchanged resume takes the cheap completion path');
	put($input, join('', map { ">$_\n".('ACGT' x 30)."\n" } qw(A B C D)));
	is(run($log, @cmd, '-continue', 1), 0, 'changed single-locus input rebuilds');
	ok(getTreeLeafs("$output/phylo/FASTTREE_allsites.nwk")->{D}, 'new sample is present in the rebuilt tree');
	is(scalar(() = slurp($ENV{TREE_TEST_COUNT}) =~ /run/g), 2, 'input change reruns inference exactly once');
	is(run($log, @cmd, '-continue', 1, '-iqLegacy', 1), 0, 'downstream policy change rebuilds tree stages');
	is(scalar(() = slurp($ENV{TREE_TEST_COUNT}) =~ /run/g), 3, 'downstream setting is not hidden by single-locus completion');
};

done_testing();
