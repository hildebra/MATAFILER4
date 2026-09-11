use strict;
use warnings;
no warnings 'once';
use FindBin qw($Bin);
use lib "$Bin/..", "$Bin/lib";
use MFTestConfig;
use File::Path qw(make_path);
use File::Temp qw(tempdir);
use File::Basename qw(dirname);
use Text::ParseWords qw(shellwords);
use Test::More;

sub write_file {
	my ($path, $content) = @_;
	make_path(dirname($path));
	open my $out, '>', $path or die "$path: $!";
	print {$out} $content;
	close $out or die "$path: $!";
}
sub slurp {
	open my $in, '<', $_[0] or die "$_[0]: $!";
	local $/;
	return <$in>;
}
sub read_table {
	my @lines = split /\n/, slurp($_[0]);
	my @columns = split /\t/, shift @lines;
	return [map {
		my @values = split /\t/, $_, -1;
		my %row; @row{@columns} = @values; \%row;
	} @lines];
}
my $source = slurp("$Bin/../secScripts/MGS/strain_within.pl");
my @helpers;
for my $name (qw(retainedMSAInputMB mgsTreePath mgsOutputComplete
	msaOnlyArtifactsReady preparedOutgroupLog lifecycleMarkerReason
	epaOnlyRetryReady prepareEpaOnlyRetryState persistentMGSInputState
	scratchMGSInputState validateTreeInputResolution writeTooFewMarker
	writeNoRecoverableLociMarker writeStrainWorkflowState
	cleanupLegacyStrainWorkflowStateFiles writeStrainWorkflowHeartbeat
	writeStrainWorkflowFailure createConsFastas reduceSeqTech shellQuote)) {
	my ($sub) = $source =~ /(sub \Q$name\E\s*\{.*?^\})/ms;
	BAIL_OUT("Missing controller helper $name") unless $sub;
	push @helpers, $sub;
}
my ($preparedReuse) = $source =~ /(sub addOutgroup2MGS\{.*?)(?=\n\tmy \$rawCategory =)/s;
BAIL_OUT('Missing prepared-input reuse block') unless $preparedReuse;
push @helpers, $preparedReuse . "\n\treturn ('raw');\n}\n";
my ($lean) = $source =~ /(\tmy %resumeEntry;.*?)(?=\tmy \$publishedInputsReady)/s;
BAIL_OUT('Missing active lean resume block') unless $lean;
my $setup = <<'SETUP';
package StrainControllerFixture;
use File::Spec;
use Digest::SHA qw(sha256_hex);
use File::Path qw(make_path);
use File::Glob qw(bsd_glob);
use Mods::GenoMetaAss qw(fileGZe fileGZs gzipopen);
use Mods::WorkflowResilience qw(atomic_write_text write_workflow_record
	retry_open retry_close retry_rename retry_unlink);
our ($onlyMSA, $strictBackbone, $phyloProg, $onlySubmit, $recalcTrees, $reSubmit,
	$repairCAT, $deepRepair, $redoSubmissionData, $outD, $scratchD, $LOGDIR,
	$FNAstdof, $FAAstdof, $CATstdof, $LINKstdof, $QCstdof,
	$workflowStatePath, $workflowStage, $workflowStatus, $workflowReason,
	$SaSe, $leanOnlySubmitResume,
	$legacyWorkflowHeartbeatPath, $legacyWorkflowFailurePath,
	$lSNPdir, $lConsVCF, $lConsVCFsup, $lMAPdir, $bamDepthFsuffix,
	$bamDepthFsuffixSup, $noIndels, $minSNPDepth, $minSNPCallQual,
	$useAdaptiveQual, $depthFilterScale, $indelRange, $SNPconsLOGs);
our (@specis, %SIdirs, %ConspecificMGS, %persistentMGSInputStateCache,
	%scratchMGSInputStateCache, %map, %staged, %legacyLocusMGS, %stagedShardHandoff);
our ($stagedProbes, $sequenceScans, $lastCommand);
sub stagedMGSInputsReady { $stagedProbes++; return $staged{$_[0]} || 0; }
sub count_msa_samples { $sequenceScans++; die 'Unexpected alignment sequence scan'; }
sub limitedNotice { }
sub getProgPaths { return 'vcf2fna'; }
sub getAssemblContigs { return "$_[0]/assembly/contigs.fna"; }
sub getAssemblGFF { return "$_[0]/assembly/genes.gff"; }
sub systemW { $lastCommand = $_[0]; }
SETUP
my $wrapper = <<'LEAN';
sub lean_decision {
	my ($outD2, $ensureLocusMSAs) = @_;
	my $leanOnlySubmitResume = 1;
	my (%treeDisposition, %MGSepaOnlyRetry);
	my ($epaOnlyRetry, $epaRecovery, $epaOnlyRetryCount) = (0, 0, 0);
	MGS_SUBMISSION: for my $MGS ('test') {
LEAN
$wrapper .= $lean . <<'END_LEAN';
	}
	return { disposition => \%treeDisposition, epa => \%MGSepaOnlyRetry,
		epaOnlyRetry => $epaOnlyRetry, epaRecovery => $epaRecovery,
		epaOnlyRetryCount => $epaOnlyRetryCount };
}
1;
END_LEAN
my $loaded = eval $setup . join("\n", @helpers) . $wrapper;
BAIL_OUT($@) unless $loaded;
my $tmp = tempdir(CLEANUP => 1);
$StrainControllerFixture::LOGDIR = "$tmp/logs";
$StrainControllerFixture::outD = "$tmp/out";
$StrainControllerFixture::scratchD = "$tmp/scratch";
$StrainControllerFixture::onlyMSA = 0;
$StrainControllerFixture::phyloProg = 1;
$StrainControllerFixture::strictBackbone = 1;
$StrainControllerFixture::onlySubmit = 1;
($StrainControllerFixture::recalcTrees, $StrainControllerFixture::reSubmit,
	$StrainControllerFixture::repairCAT, $StrainControllerFixture::deepRepair,
	$StrainControllerFixture::redoSubmissionData) = (0) x 5;
($StrainControllerFixture::FNAstdof, $StrainControllerFixture::FAAstdof,
	$StrainControllerFixture::CATstdof, $StrainControllerFixture::LINKstdof,
	$StrainControllerFixture::QCstdof) = qw(all.fna all.faa all.cat all.link all.qc);

subtest 'Retained-alignment sizing uses MiB and only file metadata' => sub {
	my $dir = "$tmp/sizing";
	for my $suffix ('', '.gz') {
		my $path = "$dir/MSA/MSAli.fna$suffix";
		write_file($path, '');
		open my $out, '+<', $path or die $!;
		truncate($out, 8 * 1024 * 1024) or die $!;
		close $out;
		is(StrainControllerFixture::retainedMSAInputMB($dir), 8,
			"8 MiB $suffix file stays 8 MiB; compressed data need not be read");
		unlink $path;
	}
	is(StrainControllerFixture::retainedMSAInputMB($dir), 1,
		'missing input retains the one-MiB scheduling floor');
};

subtest 'Current completion files and legacy nonempty tree markers are reusable' => sub {
	my $dir = "$tmp/completed";
	for my $program (1 .. 3) {
		$StrainControllerFixture::phyloProg = $program;
		my $tree = StrainControllerFixture::mgsTreePath($dir);
		write_file($tree, '(a,b,c);');
		write_file("$dir/treeDone.sto", '');
		ok(!StrainControllerFixture::mgsOutputComplete($dir), 'empty marker is not completion');
		my $decision = StrainControllerFixture::lean_decision($dir, 0);
		is_deeply($decision->{disposition}, {}, 'lean dispatch does not skip an empty completion marker');
		write_file("$dir/treeDone.sto", "done\n");
		ok(StrainControllerFixture::mgsOutputComplete($dir), "phyloProg $program accepts the existing plain marker");
		$decision = StrainControllerFixture::lean_decision($dir, 0);
		is($decision->{disposition}{'valid tree already present'}, 1, 'lean dispatch reuses the completed tree');
		rename $tree, "$tree.gz" or die $!;
		ok(StrainControllerFixture::mgsOutputComplete($dir), 'compressed tree remains usable without reading it');
		unlink "$tree.gz";
		ok(!StrainControllerFixture::mgsOutputComplete($dir), 'a marker without its tree stays incomplete');
	}
	$StrainControllerFixture::phyloProg = 1;
	$StrainControllerFixture::onlyMSA = 1;
	# The readiness check must not decompress these deliberately opaque bytes.
	write_file("$dir/MSA/locus.fna.gz", 'retained alignment');
	for my $metadata ("status\tmsa_complete\n",
		"status\tmsa_complete\nmsa_samples\t3\nmsa_outgroup_samples\t1\n") {
		write_file("$dir/msaOnly.complete.tsv", $metadata);
		ok(StrainControllerFixture::mgsOutputComplete($dir), 'current and pre-count MSA markers are accepted');
		my %counts;
		ok(StrainControllerFixture::msaOnlyArtifactsReady($dir, \%counts),
			'sample reporting also avoids reading retained alignment contents');
		is($counts{msa_samples}, $metadata =~ /msa_samples/ ? 3 : 'NA',
			'only previously measured counts are reported');
		is(StrainControllerFixture::lean_decision($dir, 0)->{disposition}{'valid MSA already present'},
			1, 'lean dispatch reuses the retained MSA without scanning sequences');
	}
	$StrainControllerFixture::onlyMSA = 0;
};

subtest 'Lean resume ignores interrupted empty terminal markers' => sub {
	my $dir = "$tmp/terminal";
	for my $marker (qw(tooFewSamples.sto noRecoverableLoci.sto noTree.sto)) {
		write_file("$dir/$marker", '');
		is_deeply(StrainControllerFixture::lean_decision($dir, 0)->{disposition}, {},
			"empty $marker does not suppress a retry");
		write_file("$dir/$marker", "reason\tinsufficient_input\n");
		is(StrainControllerFixture::lean_decision($dir, 0)->{disposition}{'valid no-tree: insufficient_input'},
			1, "current $marker contents retain their meaning");
		unlink "$dir/$marker";
	}
	StrainControllerFixture::writeTooFewMarker($dir, 2, 12);
	is(slurp("$dir/tooFewSamples.sto"), "reason\ttoo_few_samples\nsamples\t2\ngenes\t12\n",
		'atomic publication retains the too-few marker format');
	StrainControllerFixture::writeNoRecoverableLociMarker($dir);
	is(slurp("$dir/noRecoverableLoci.sto"), "reason\tempty_extraction\n",
		'atomic publication retains the no-locus marker format');
};

my $epa = "$tmp/epa";
for my $file (qw(MSA/MSAli.fna.gz MSA/MSAli.placement.fna.gz
	phylo/IQtree_allsites.backbone.treefile phylo/IQtree_allsites.backbone.log
	phylo/strict_backbone.samples.tsv)) {
	write_file("$epa/$file", 'retained');
}
subtest 'Active lean dispatch recognizes legacy EPA recovery without a pending marker' => sub {
	my $decision = StrainControllerFixture::lean_decision($epa, 0);
	is($decision->{epa}{test}, 'legacy_missing_final', 'legacy retained inputs select EPA-only recovery');
	is_deeply([@{$decision}{qw(epaOnlyRetry epaRecovery epaOnlyRetryCount)}], [1, 1, 1],
		'dispatch actually enables the isolated retry, not merely the readiness helper');
	write_file("$epa/treeDone.sto", "done\n");
	ok(StrainControllerFixture::prepareEpaOnlyRetryState($epa, $decision->{epa}{test}),
		'legacy recovery is prepared');
	ok(!-e "$epa/treeDone.sto", 'stale completion is cleared when its final tree is absent');
	is(StrainControllerFixture::lean_decision($epa, 0)->{epa}{test}, 'explicit_pending',
		'current placement marker is also accepted');
	$StrainControllerFixture::onlyMSA = 1;
	is_deeply(StrainControllerFixture::lean_decision($epa, 0)->{epa}, {},
		'MSA-only dispatch never selects tree placement recovery');
	$StrainControllerFixture::onlyMSA = 0;
	write_file("$epa/phylo/epa-ng/epa_result.jplace", 'retained placement');
	is_deeply(StrainControllerFixture::lean_decision($epa, 0)->{epa}, {},
		'existing placement output prevents unnecessary EPA inference');
	unlink "$epa/phylo/epa-ng/epa_result.jplace";
};

subtest 'Final resolution trusts completed outputs and refreshes only unfinished inputs' => sub {
	my $dir = "$tmp/resolution";
	@StrainControllerFixture::specis = qw(done missing published removed terminal placement);
	%StrainControllerFixture::SIdirs = map { $_ => "$dir/$_" } @StrainControllerFixture::specis;
	$StrainControllerFixture::SIdirs{placement} = $epa;
	write_file("$dir/done/treeDone.sto", "done\n");
	write_file("$dir/done/phylo/IQtree_allsites.treefile", '(a,b,c);');
	write_file("$dir/terminal/noTree.sto", "reason\tmasked_alignment\n");
	# Model states from before jobs published and then cleaned their inputs.
	%StrainControllerFixture::persistentMGSInputStateCache = (published => 'missing', removed => 'complete');
	%StrainControllerFixture::scratchMGSInputStateCache = (published => 'missing', removed => 'complete');
	write_file("$dir/published/$_", 'records') for qw(all.fna all.faa all.cat);
	$StrainControllerFixture::stagedProbes = 0;
	is(StrainControllerFixture::validateTreeInputResolution(), 2, 'only missing and removed inputs require repair');
	my $rows = read_table("$StrainControllerFixture::LOGDIR/tree_input_resolution.tsv");
	my %row = map { $_->{MGS} => $_ } @$rows;
	is($row{done}{resolution}, 'tree_output_complete', 'cleaned-up inputs do not reopen a completed tree');
	is($row{done}{persistent_state}, 'not_required', 'completed output needs no input-sidecar lookup');
	is($row{published}{resolution}, 'tree_input_ready', 'new publication replaces a cached missing state');
	is($row{removed}{resolution}, 'repair_required', 'removed files invalidate an old complete state');
	is($row{placement}{resolution}, 'placement_retry_ready', 'retained EPA inputs need no extraction repair');
	is($StrainControllerFixture::stagedProbes, 3, 'scratch is checked only for unfinished, nonterminal MGS');
	my $queue = read_table("$StrainControllerFixture::LOGDIR/tree_input_repair.queue.tsv");
	is_deeply([map { $_->{MGS} } @$queue], [qw(missing removed)], 'repair queue contains only unresolved MGS');
	@StrainControllerFixture::specis = ('done');
	$StrainControllerFixture::onlyMSA = 1;
	write_file("$dir/done/msaOnly.complete.tsv", "status\tmsa_complete\n");
	write_file("$dir/done/MSA/locus.fna.gz", 'opaque retained MSA');
	$StrainControllerFixture::stagedProbes = 0;
	is(StrainControllerFixture::validateTreeInputResolution(), 0, 'legacy MSA completion is also resolved without inputs');
	is($StrainControllerFixture::stagedProbes, 0, 'a completed MSA resume never probes scratch');
	is($StrainControllerFixture::sequenceScans || 0, 0, 'all readiness checks avoid sequence scans');
	ok(!-e "$StrainControllerFixture::LOGDIR/tree_input_repair.queue.tsv", 'resolved runs clear stale repair queues');
	$StrainControllerFixture::onlyMSA = 0;
};

subtest 'Published and legacy prepared scratch inputs share one reuse path' => sub {
	my $published = "$tmp/prepared/published";
	my $scratch = "$tmp/prepared/scratch";
	$StrainControllerFixture::SIdirs{prepared} = $published;
	$StrainControllerFixture::SaSe = '|';
	$StrainControllerFixture::leanOnlySubmitResume = 0;
	for my $dir ($published, $scratch) {
		write_file("$dir/$_", 'records') for qw(all.fna all.faa all.link all.qc merge.complete.tsv);
		write_file("$dir/data.log", "OG:outgroup\r\n");
		write_file("$dir/all.cat", "s1|c|g1\ts2|c|g1\toutgroup|c|g1\r\n\r\ns1|c|g2\ts2|c|g2\r\n");
	}
	is_deeply([StrainControllerFixture::addOutgroup2MGS('prepared', 'old_outgroup', $scratch)],
		[3, 2, 'outgroup', 0, 1, 2], 'published input wins and counts exclude its saved outgroup');
	unlink "$published/all.faa";
	is_deeply([StrainControllerFixture::addOutgroup2MGS('prepared', '', $scratch)],
		[3, 2, 'outgroup', 1, 1, 2], 'complete old scratch input is reused with a copy request');
	write_file("$scratch/merge.complete.tsv", '');
	is_deeply([StrainControllerFixture::addOutgroup2MGS('prepared', '', $scratch)], ['raw'],
		'uncommitted prepared scratch input is not trusted');
	write_file("$scratch/merge.complete.tsv", "complete\n");
	unlink "$scratch/all.faa";
	is_deeply([StrainControllerFixture::addOutgroup2MGS('prepared', '', $scratch)], ['raw'],
		'full auditing requires the complete scratch input set');
	$StrainControllerFixture::leanOnlySubmitResume = 1;
	is_deeply([StrainControllerFixture::addOutgroup2MGS('prepared', '', $scratch)],
		[3, 2, 'outgroup', 1, 1, 2], 'lean reuse retains its existing checkpoint and category fast path');
	$StrainControllerFixture::stagedShardHandoff{prepared} = {};
	is_deeply([StrainControllerFixture::addOutgroup2MGS('prepared', '', $scratch)], ['raw'],
		'new worker shards take precedence over old scratch publication');
	delete $StrainControllerFixture::stagedShardHandoff{prepared};
	$StrainControllerFixture::repairCAT = 1;
	is_deeply([StrainControllerFixture::addOutgroup2MGS('prepared', '', $scratch)], ['raw'],
		'explicit category repair bypasses prepared inputs');
	$StrainControllerFixture::repairCAT = 0;
	$StrainControllerFixture::leanOnlySubmitResume = 0;
};

subtest 'Workflow state separates planned, partial, completed, and failed outcomes' => sub {
	$StrainControllerFixture::workflowStatePath = "$tmp/workflow.state.tsv";
	$StrainControllerFixture::legacyWorkflowHeartbeatPath = '';
	$StrainControllerFixture::legacyWorkflowFailurePath = '';
	StrainControllerFixture::writeStrainWorkflowHeartbeat('phase1');
	for my $status (qw(planned partial)) {
		StrainControllerFixture::writeStrainWorkflowHeartbeat(undef, $status, 'worker repair or submission needed');
		my $row = read_table($StrainControllerFixture::workflowStatePath)->[0];
		is($row->{status}, $status, "$status is persisted explicitly");
		is($row->{stage}, 'phase1', 'the unfinished stage remains visible');
		is($row->{reason}, 'worker repair or submission needed', 'the reason is retained');
	}
	StrainControllerFixture::writeStrainWorkflowHeartbeat('complete');
	is(read_table($StrainControllerFixture::workflowStatePath)->[0]{status}, 'completed', 'ordinary completion is unchanged');
	StrainControllerFixture::writeStrainWorkflowFailure('test failure');
	is(read_table($StrainControllerFixture::workflowStatePath)->[0]{status}, 'failed', 'failure reporting is unchanged');
	# Inspect each actual controlled early exit so none falls through END as success.
	my @repairExits = $source =~ /(\$completionMessage = "Phase I (?:requires|generation).*?exit\(0\);)/sg;
	is(scalar(@repairExits), 3, 'all three Phase-I repair exits are covered');
	like($_, qr/writeStrainWorkflowHeartbeat\(undef, 'partial', \$completionMessage\)/,
		'Phase-I repair records partial status before exiting') for @repairExits;
};

subtest 'Consensus command simplification preserves primary and hybrid options' => sub {
	($StrainControllerFixture::lSNPdir, $StrainControllerFixture::lConsVCF,
		$StrainControllerFixture::lConsVCFsup, $StrainControllerFixture::lMAPdir,
		$StrainControllerFixture::bamDepthFsuffix, $StrainControllerFixture::bamDepthFsuffixSup) =
		('SNP', 'primary.vcf', 'support.vcf', '/mapping', '.depth', '.support.depth');
	($StrainControllerFixture::minSNPDepth, $StrainControllerFixture::minSNPCallQual,
		$StrainControllerFixture::useAdaptiveQual, $StrainControllerFixture::depthFilterScale,
		$StrainControllerFixture::indelRange) = (1, 3, 0, 0.1, 5);
	$StrainControllerFixture::SNPconsLOGs = "$tmp/consensus log.txt";
	for my $case (['', '', 'ill'], ['miSeq', 'PB:reads', 'ill,PB'],
		['ill', 'ONT:reads', 'ill,ONT'], ['ONT', '', 'ONT'],
		['PB', '', 'PB'], ['unknown', '', 'unspecified'], ['0', '', 'unspecified']) {
		my ($primary, $support, $platforms) = @$case;
		$StrainControllerFixture::map{sample} = { SeqTech => $primary, SupportReads => $support };
		for my $skip (0, 1) {
			$StrainControllerFixture::noIndels = $skip;
			my $dir = "$tmp/input with spaces'and quote";
			my $cmd = StrainControllerFixture::createConsFastas($dir, 'sample', 'genes.fna', 'genes.faa', 1, 1);
			my @args = shellwords($cmd);
			is(shift @args, 'vcf2fna', 'configured executable is preserved');
			my (%options, $hasSkip);
			while (@args) {
				my $key = shift @args;
				if ($key eq '-skipINDELs') { $hasSkip = 1; next; }
				$options{$key} = shift @args;
			}
			is($options{'-seqPlatform'}, $platforms, 'platform ordering and normalization are unchanged');
			is($hasSkip || 0, $skip, 'indel flag follows the requested policy');
			is_deeply([@options{qw(-minCallDepth -minCallQual -minCallQualAdaptive -depthFilterScale -indelRange)}],
				[1, 3, 0, 0.1, 5], 'both command paths retain identical SNP thresholds');
			is($options{'-inVCF'}, "$dir/SNP/primary.vcf".($support ? ",$dir/SNP/support.vcf" : ''),
				'VCF paths remain a correctly quoted paired argument');
			is($options{'-depthF'}, "$dir/mapping/sample.depth".($support ? ",$dir/mapping/sample.support.depth" : ''),
				'depth paths remain paired in the same order');
			is($options{'-oCtg'} // '', $support ? '/dev/null' : '', 'hybrid contig output behavior is unchanged');
			is($options{'>>'}, $StrainControllerFixture::SNPconsLOGs, 'log destination survives shell quoting');
		}
	}
	my $expected = StrainControllerFixture::createConsFastas('/input', 'sample', 'nt', 'aa', 0, 1);
	StrainControllerFixture::createConsFastas('/input', 'sample', 'nt', 'aa', 0, 0);
	is($StrainControllerFixture::lastCommand, $expected, 'execution and generated scripts use the same command');
};

done_testing();
