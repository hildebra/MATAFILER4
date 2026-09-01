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
use Mods::Checkpoint qw(write_checkpoint checkpoint_valid);

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
my $sorter_source = slurp($sorter);
like($sorter_source, qr/Can't open temporary outfile \$tmpout: \$!\\n/,
	'sorter open failures retain the operating-system reason');
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
like($sorter_stdout, qr/^MGS1 \(2\)::/m,
	'sorter retains a per-MGS detail preview');
like($sorter_stdout,
	qr/MGS detail preview: 1\/1 shown.*?Sorting summary: MGS=1, input_rows=2.*?elapsed=\d+s/s,
	'sorter follows its preview with aggregate totals and elapsed time');

# Catalogue-scale runs used to print one long diagnostic row for every MGS.
# Exercise more groups than the preview limit and pin both suppression and the
# aggregate replacement so later workflow stages remain visible in stdout.
my $preview_mgs = File::Spec->catfile($tmp, 'preview.core');
write_file($preview_mgs, join('', map {
	"MGS$_\tpreviewGene$_\t10\t0\t1\tunused\n"
} 1 .. 7));
write_file(File::Spec->catfile($tmp, 'preview.obs'), join('', map {
	"MGS$_\t10\n"
} 1 .. 7));
my $preview_err = gensym;
my $preview_pid = open3(
	undef, my $preview_out, $preview_err,
	$^X, '-I'.File::Spec->catdir($Bin, '..'),
	$sorter, $gc, $preview_mgs, 'FMG', 'test', 95,
);
my $preview_stdout = do { local $/; <$preview_out> // '' };
my $preview_stderr = do { local $/; <$preview_err> // '' };
waitpid($preview_pid, 0);
is($? >> 8, 0, 'gene-priority sorter completes a multi-MGS preview run')
	or diag($preview_stdout, $preview_stderr);
my @preview_details = $preview_stdout =~ /^(MGS\d+ \([^\n]+)$/mg;
is(scalar(@preview_details), 5,
	'sorter reports only the first five per-MGS detail rows');
unlike($preview_stdout, qr/^MGS6 \(/m,
	'sorter suppresses per-MGS detail after the preview');
like($preview_stdout,
	qr/MGS detail preview: 5\/7 shown; 2 omitted.*?Sorting summary: MGS=7, input_rows=7, unique_genes=7, ranked_genes=7, emitted_genes=7/s,
	'sorter summarizes every MGS after the bounded preview');
is(scalar(split /\n/, slurp("$preview_mgs.srt")), 7,
	'output suppression does not remove any sorted guide records');

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

# The between-MGS worker must consume the predefined GTDB catalog files when
# requested, without falling back to an FMG subset. Two marker-bearing MGS are
# enough to exercise selection and then reach the intentional sparse-tree skip.
my $gtdb_gc = File::Spec->catdir($tmp, 'GTDB-GC');
my $gtdb_out = File::Spec->catdir($tmp, 'GTDB-tree');
make_path($gtdb_gc);
write_file(
	File::Spec->catfile($gtdb_gc, 'GTDBmg.subset.cats'),
	"bac120_marker\t2\tgeneG1,geneG2\n",
);
my $gtdb_mgs = File::Spec->catfile($tmp, 'gtdb.core');
write_file($gtdb_mgs, "MGS1\tgeneG1\nMGS2\tgeneG2\n");
my $between = File::Spec->catfile(
	$Bin, '..', 'secScripts', 'MGS', 'phylo_MGS_between.pl',
);
my $between_err = gensym;
my $between_pid = open3(
	undef, my $between_out, $between_err,
	$^X, '-I'.File::Spec->catdir($Bin, '..'),
	$between, '-GCd', $gtdb_gc, '-MGS', $gtdb_mgs,
	'-MGset', 'GTDB', '-outD', $gtdb_out,
);
my $between_stdout = do { local $/; <$between_out> // '' };
my $between_stderr = do { local $/; <$between_err> // '' };
waitpid($between_pid, 0);
is($? >> 8, 0, 'between-MGS worker accepts predefined GTDB marker genes')
	or diag($between_stdout, $between_stderr);
like(
	$between_stdout,
	qr/Found 2 GTDB marker genes.*?SKIPPED=too_few_marker_bearing_MGS:2/s,
	'GTDB selection reaches the expected sparse-tree skip using GTDB marker-bearing MGS',
);

# MGS 0.55 changed only the strain launcher. A global release fingerprint used
# to invalidate otherwise unchanged GTDB outputs written by 0.54.
my $resume_core = File::Spec->catfile($tmp, 'resume.core');
my $resume_tax = File::Spec->catfile($tmp, 'GTDBTK.tax');
my $resume_summary = File::Spec->catfile($tmp, 'gtdbtk.summary.tsv');
write_file($resume_core, "MGS1\tgene1\n");
write_file($resume_tax, "MGS1\td__Bacteria\n");
write_file($resume_summary, "user_genome\tclassification\nMGS1\td__Bacteria\n");
my $legacy_gtdb_stone = File::Spec->catfile($tmp, 'GTDBTK.stone');
my %legacy_gtdb_parameters = (
	pipeline_version => '0.54',
	stage => 'gtdb-taxonomy',
	catalog_identity => 'catalog-A',
	cluster_id => 95,
	binner => 'SC',
	catalog_fna_size => 123,
	catalog_fna_mtime => 456,
);
write_checkpoint(
	$legacy_gtdb_stone,
	parameters => \%legacy_gtdb_parameters,
	outputs => [$resume_tax, $resume_summary, $resume_core],
);
my %new_global_parameters = (%legacy_gtdb_parameters, pipeline_version => '0.55');
ok(
	!checkpoint_valid($legacy_gtdb_stone, parameters => \%new_global_parameters),
	'a controller-only version bump invalidates the old whole-pipeline predicate',
);
my %gtdb_stage_parameters = map {
	$_ => $legacy_gtdb_parameters{$_}
} qw(stage catalog_identity cluster_id binner catalog_fna_size catalog_fna_mtime);
ok(
	checkpoint_valid($legacy_gtdb_stone, parameters => \%gtdb_stage_parameters),
	'unchanged GTDB direct inputs and outputs remain resumable across that release bump',
);

# The same compatibility rule applies to downstream abundance checkpoints:
# every current parameter except the release number must match, and the
# checkpoint output records remain authoritative.
my $legacy_ab_output = File::Spec->catfile($tmp, "MGS.matL7.txt");
write_file($legacy_ab_output, "MGS1\t1\n");
my $legacy_ab_stone = File::Spec->catfile($tmp, "abund.mgs_core.stone");
my %legacy_ab_parameters = (
	pipeline_version => "0.54",
	stage => "marker-mgs-abundance",
	catalog_identity => "catalog-A",
	cluster_id => 95,
	marker_set => "GTDB",
	binner => "SC",
	quality_checker => "checkm2",
);
write_checkpoint(
	$legacy_ab_stone,
	parameters => \%legacy_ab_parameters,
	outputs => [$legacy_ab_output],
);
my %current_ab_parameters = (%legacy_ab_parameters, pipeline_version => "0.55");
ok(
	!checkpoint_valid($legacy_ab_stone, parameters => \%current_ab_parameters),
	"the exact predicate rejects the otherwise compatible 0.54 abundance checkpoint",
);
delete $current_ab_parameters{pipeline_version};
ok(
	checkpoint_valid($legacy_ab_stone, parameters => \%current_ab_parameters),
	"removing only the proven-compatible release field retains all abundance parameters",
);
write_file($legacy_ab_output, "MGS1\tchanged\n");
ok(
	!checkpoint_valid($legacy_ab_stone, parameters => \%current_ab_parameters),
	"compatible-release reuse still rejects a changed recorded abundance output",
);

my $mgs_entrypoint = slurp(File::Spec->catfile($Bin, '..', 'secScripts', 'MGS.pl'));
like(
	$mgs_entrypoint,
	qr/my %stageCheckpointContract = \(.*?'extract-bin-contigs'\s*=>\s*1.*?'gtdb-taxonomy'\s*=>\s*1.*?sub _checkpoint_parameters_for_stage/s,
	'MGS declares explicit contracts for bin extraction and GTDB taxonomy',
);
like(
	$mgs_entrypoint,
	qr/sub _checkpoint_parameters_for_stage .*?catalog_identity cluster_id binner genomes_per_family.*?catalog_fna_size catalog_fna_mtime.*?catalog_faa_size catalog_faa_mtime.*?mag_report_size.*?stage_contract/s,
	'bin-extraction resume tracks its direct catalogue, map, MAG-report, and contract inputs',
);
like(
	$mgs_entrypoint,
	qr/sub _checkpoint_valid_for_stage .*?read_checkpoint\(\$file\).*?pipeline_version.*?0\\\.\(\?:54\|55\).*?created_epoch.*?checkpoint_valid\(\$file, parameters => \\%legacyParameters\)/s,
	'legacy checkpoint adoption is limited to compatible JSON releases and stable recorded outputs',
);
like(
	$mgs_entrypoint,
	qr/my \$binExtractionValid = _checkpoint_valid_for_stage\(\$BinExtrSto, 'extract-bin-contigs'\).*?my \$gtdbTaxonomyValid = .*?_checkpoint_valid_for_stage\(\$GTDBtaxSto, 'gtdb-taxonomy'\)/s,
	'bin extraction and GTDB resume use their stage-specific validators',
);

like(
	$mgs_entrypoint,
	qr/sub _checkpoint_valid_for_compatible_release .*?MGSpipelineVersion.*?0\.55.*?pipeline_version.*?0\.54.*?delete \$compatibleParameters\{pipeline_version\}.*?checkpoint_valid\(\$file, parameters => \\%compatibleParameters\).*?sub _checkpoint_valid .*?_checkpoint_valid_for_compatible_release.*?sub _checkpoint_valid_for_resume .*?_checkpoint_valid_for_compatible_release/s,
	"global checkpoint reuse is restricted to the proven 0.54-to-0.55 transition",
);
like(
	$mgs_entrypoint,
	qr/my \@stage1AssignmentFiles = \(\$finalClusters2, "\$\{finalClusters2\}UW", \$finalClustersW\);.*?my \$stage1AssignmentsPresent = grep \{ _mgs_count\(\$_, 1\) \} \@stage1AssignmentFiles;.*?my \$stage1ResumeValid = !\$rewrClusterMAGs && \$stage1AssignmentsPresent/s,
	"Stage I resumes from primary or recoverable weighted assignments without requiring core",
);

my ($mgs_assignment_helpers) = $mgs_entrypoint =~
	/(sub _mgs_ids \{.*?\n\}\n\nsub _mgs_count \{.*?\n\})/s;
ok(defined($mgs_assignment_helpers), 'MGS assignment helpers can be exercised independently');
my $assignment_helpers_loaded = defined($mgs_assignment_helpers)
	&& eval "package MGSResumeAssignmentProbe; use strict; use warnings; $mgs_assignment_helpers; 1;";
ok($assignment_helpers_loaded, 'MGS assignment helpers compile for the resume regression')
	or diag($@);
if ($assignment_helpers_loaded) {
	# Header emitted by bin/clusterMAGs (the default engine) for .clusters and
	# .Wclusters, and by secScripts/MGS/clusterMAGs.pl for the -perlClusterMAGs
	# compatibility path.  Both must be recognised: counting either as an MGS
	# inflates every assignment count by one.
	my $binary_header = "MGS\tGene\tOcc\tMultiCopy\tMultiBin\tisMarkerGene\n";
	my $weighted_header = "MGS\tGene\tOcc\tW_MultiCopy\tW_MultiBin\tisMarkerGene\n";
	my $compat_header = "Bin\tGene\tOcc\tMultiCopy\tMultiBin\tisMarkerGene\n";
	my $assignment_row = sub {
		my ($mgs, $gene) = @_;
		return "$mgs\t$gene\t10\t0\t1\t0\n";
	};

	my %header_only = (
		'clusterMAGs binary'      => $binary_header,
		'weighted clusterMAGs'    => $weighted_header,
		'Perl compatibility path' => $compat_header,
	);
	for my $engine (sort keys %header_only) {
		my $header_only_assignments = File::Spec->catfile($tmp, "header-only.$engine.MGS");
		write_file($header_only_assignments, $header_only{$engine});
		is(
			MGSResumeAssignmentProbe::_mgs_count($header_only_assignments, 1),
			0,
			"a header-only $engine assignment file cannot satisfy the Stage I resume gate",
		);
	}

	# A genuine singleton must report exactly one MGS: MGS.pl synthesizes the
	# missing observation table only for $activeMGSCount == 1, and
	# _write_single_mgs_observations dies unless _mgs_ids returns one id.
	my $singleton_assignments = File::Spec->catfile($tmp, 'singleton.MGS');
	write_file($singleton_assignments, $binary_header.$assignment_row->('MGS.1', 'gene1'));
	is(
		MGSResumeAssignmentProbe::_mgs_count($singleton_assignments),
		1,
		'a single clustered MGS is counted once, not alongside its header',
	);

	my $multi_assignments = File::Spec->catfile($tmp, 'multi.MGS');
	write_file($multi_assignments, $binary_header
		.$assignment_row->('MGS.1', 'gene1')
		.$assignment_row->('MGS.1', 'gene2')
		.$assignment_row->('MGS.2', 'gene3'));
	is(
		MGSResumeAssignmentProbe::_mgs_count($multi_assignments),
		2,
		'distinct MGS are counted without the header contributing an extra id',
	);

	# The post-filtered .core table carries no header, so its leading row is a
	# real assignment and must never be skipped.
	my $core_assignments = File::Spec->catfile($tmp, 'core.MGS');
	write_file($core_assignments, $assignment_row->('MGS.1', 'gene1')
		.$assignment_row->('MGS.2', 'gene2'));
	is(
		MGSResumeAssignmentProbe::_mgs_count($core_assignments),
		2,
		'a headerless core table keeps its first assignment',
	);

	my $early_stop_assignments = File::Spec->catfile($tmp, 'early-stop.MGS');
	write_file($early_stop_assignments,
		$binary_header.$assignment_row->('MGS.1', 'gene1')."malformed\n");
	is(
		MGSResumeAssignmentProbe::_mgs_count($early_stop_assignments, 1),
		1,
		'the resume probe stops after the first assignment instead of rescanning the file',
	);
}
like(
	$mgs_entrypoint,
	qr/my \$recoverMissingActiveMGS = sub .*?\$recoverMissingActiveMGS->\(\) unless \$ph1flag;/s,
	"an interrupted weighted handoff is recovered before clustering can restart",
);
like(
	$mgs_entrypoint,
	qr/systemW \$postCmd if !-s \$finalClustersFilt;.*?_touch_checkpoint\(\$st1ston, .stage-1., \$finalClusters2, \$finalClustersFilt\)\s+if \$ph1flag;/s,
	"the cheap core is regenerated but only newly clustered primary/core outputs receive provenance",
);
unlike(
	$mgs_entrypoint,
	qr/_touch_checkpoint\(\$st1ston, .stage-1., .*?unless _checkpoint_valid\(\$st1ston\)/s,
	"a mismatched recovered Stage I product is not laundered into a current checkpoint",
);
like(
	$mgs_entrypoint,
	qr/my \@geneBinFiles = grep \{ -f \$_ && \/\\\.\(\?:fna\|faa\)\\z\/i \} glob\("\$binD\/\*"\);/,
	"BinExtr records only its FASTA products, excluding downstream GTDB archives",
);
like(
	$mgs_entrypoint,
	qr/my \@strainArguments = \(.*?\x27-SNPcaller\x27, \$SNPcaller.*?map \{ _shell_quote\(\$_\) \} \@strainArguments/s,
	"MGS forwards the selected SNP caller with shell-quoted strain arguments",
);

like(
	$mgs_entrypoint,
	qr/my \$ph1flag = \$stage1ResumeValid \? 0 : 1;.*?if \(\$ph1flag\) \{.*?glob\("\$outD\/\$BinnerShrt\.clusters\*"\).*?\$invalidateMGSDerivatives->\(\);.*?\$recoverMissingActiveMGS->\(\) unless \$ph1flag;/s,
	"every Stage I rebuild removes stale core/downstream derivatives before clustering",
);
like(
	$mgs_entrypoint,
	qr/my \$invalidateMGSDerivatives = sub \{.*?\$strainLockPath = "\$lockBase\.strain_within\.lock";.*?acquire_workflow_lock\(.*?for my \$derived .*?for my \$phyloDir/s,
	"MGS invalidation acquires the strain sibling lock before deleting any derivative",
);
like(
	$mgs_entrypoint,
	qr/sub _legacy_bin_extraction_checkpoint_valid .*?GTDBtk\.tar\.gz.*?next if \$canonicalPath eq \$archivePath;.*?else \{\s*return 0;.*?my \@stat = stat\(\$path\);.*?return \$coreSeen && \$geneFnaSeen && \$geneFaaSeen && \$contigSeen.*?sub _checkpoint_valid_for_stage .*?_legacy_bin_extraction_checkpoint_valid\(\$manifest, \\%legacyParameters\)/s,
	"legacy BinExtr reuse ignores only the downstream GTDB archive and validates all primary outputs",
);

done_testing;
