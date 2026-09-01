use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use File::Path qw(remove_tree);
use FindBin qw($Bin);
use lib File::Spec->catdir($Bin, q{..});
use Test::More;

my $root = File::Spec->catdir($Bin, q{..});
my $script = File::Spec->catfile($root, q{secScripts}, q{phylo}, q{buildTree5.pl});

open my $script_handle, '<', $script or die "Cannot read $script: $!";
my $script_text = do { local $/; <$script_handle> };
close $script_handle or die "Cannot close $script: $!";

my ($epa_resource_helper) = $script_text =~
	/(sub epaResourcePlan \{.*?return \(\$threads, \$memoryMB\);\n\})/s;
my ($iqtree_explicit_helper, $iqtree_model_helper, $iqtree_partition_helper) = $script_text =~
	/(sub iqtreeExplicitEpaModel \{.*?\n\})\n\n(sub iqtreePlacementModel \{.*?\n\})\n\n(sub iqtreeGtrPartitionCount \{.*?\n\})\n\nsub epaModelArtifact/s;
BAIL_OUT('Cannot extract EPA-ng helper functions')
	unless defined($epa_resource_helper) && defined($iqtree_explicit_helper)
		&& defined($iqtree_model_helper) && defined($iqtree_partition_helper);
my $epa_helpers = "$epa_resource_helper\n$iqtree_explicit_helper\n$iqtree_model_helper\n$iqtree_partition_helper";
my $helpers_loaded = eval "package TestBuildTreeEpaHelpers; $epa_helpers; 1;";
ok($helpers_loaded, 'EPA-ng model and resource helpers load independently')
	or diag($@);
my ($classification_helper) = $script_text =~
	/(sub readStrictBackboneClassification \{.*?\n\})\n\nsub runEpaOnlyPlacement/s;
BAIL_OUT('Cannot extract strict-backbone classification reader')
	unless defined $classification_helper;
my $classification_helper_loaded = eval
	"package TestBuildTreeEpaClassification; $classification_helper; 1;";
ok($classification_helper_loaded,
	'EPA-only strict-backbone classification reader loads independently')
	or diag($@);
my ($msa_lifecycle_helpers) = $script_text =~
	/(sub copyStreamToPathAtomically \{.*?\n\})\n\nsub requireConfiguredTool/s;
BAIL_OUT('Cannot extract MSA finalization helpers')
	unless defined($msa_lifecycle_helpers);
my $msa_helpers = <<'PERL';
package TestBuildTreeMSAFinalizer;
use Cwd qw(abs_path);
use File::Basename qw(dirname);
use File::Copy qw(copy);
use File::Path qw(make_path remove_tree);
use File::Spec;
our ($pigzBin, $ncore, $removeMSA, $partiExt);
sub retry_unlink {
	my ($path, %options) = @_;
	return 1 unless -e $path || -l $path;
	unlink $path or die "Cannot unlink $path: $!\n";
	return 1;
}
sub retry_operation {
	my (%options) = @_;
	return $options{code}->();
}
sub retry_rename {
	my ($from, $to, %options) = @_;
	rename $from, $to or die "Cannot rename $from to $to: $!\n";
	return 1;
}
sub retry_open {
	my ($mode, $path, %options) = @_;
	open my $handle, $mode, $path or die "Cannot open $path: $!\n";
	return $handle;
}
sub retry_close {
	my ($handle, $label, %options) = @_;
	close $handle or die "Cannot $label: $!\n";
	return 1;
}
sub gzipopen {
	my ($path) = @_;
	my $resolved = -s $path ? $path : -s "$path.gz" ? "$path.gz" : '';
	return (undef, 0, undef) unless length $resolved;
	open my $handle, '<', $resolved or die "Cannot read $resolved: $!\n";
	return ($handle, 1, $resolved);
}
sub safeRemoveTree {
	my ($path, $parent) = @_;
	remove_tree($path);
	return 1;
}
PERL
$msa_helpers .= "$msa_lifecycle_helpers\n1;";
my $msa_helpers_loaded = eval $msa_helpers;
ok($msa_helpers_loaded, 'MSA finalization helpers load independently')
	or diag($@);
my ($staged_category_helper, $staged_qc_helper) = $script_text =~
	/(sub finalizeStagedStrainCategory \{.*?
})

(sub finalizeStagedSampleQC \{.*?
})

sub stagedOverlayRecordsPresent/s;
BAIL_OUT('Cannot extract staged strain finalization helpers')
	unless defined($staged_category_helper) && defined($staged_qc_helper);
my $staged_helpers = <<'PERL';
package TestBuildTreeStagedInputs;
use Mods::GenoMetaAss qw(fileGZs gzipopen);
use Mods::StrainPlacement qw(canonical_sample_qc_status);
sub retry_rename {
	my ($from, $to, %options) = @_;
	rename $from, $to or die "Cannot rename $from to $to: $!\n";
	return 1;
}
PERL
$staged_helpers .= "$staged_category_helper\n$staged_qc_helper\n1;";
my $staged_helpers_loaded = eval $staged_helpers;
ok($staged_helpers_loaded, 'staged strain category/QC finalizers load independently')
	or diag($@);
my ($staged_files_helper) = $script_text =~
	/(sub stagedTreeInputFiles \{.*?
})

sub publishStagedTreeInputs/s;
BAIL_OUT('Cannot extract staged tree-file listing helper') unless defined $staged_files_helper;
my $staged_files_helper_loaded = eval
	"package TestBuildTreeStagedFiles; $staged_files_helper; 1;";
ok($staged_files_helper_loaded, 'staged tree-file listing helper loads independently')
	or diag($@);
ok(index($script_text, q{use Mods::StrainParts qw(append_fasta_records_atomic);}) >= 0
	&& index($script_text, 'sub prepareStagedStrainInputs {') >= 0
	&& index($script_text, 'sub finalizeStagedStrainCategory {') >= 0
	&& index($script_text, 'sub finalizeStagedSampleQC {') >= 0
	&& index($script_text, q{my $plannedMGS = '';}) >= 0
	&& index($script_text, 'my $mgsLine = <$planIn>;') >= 0
	&& index($script_text, q{does not match MGS $plannedMGS}) >= 0,
	'buildTree5 owns staged strain finalization using only lightweight input helpers');
ok(index($script_text, q{my $stagedPlan =}) >= 0
	&& index($script_text, q{my $shardManifest =}) >= 0
	&& index($script_text, q{my $stagedPrimaryInput =}) >= 0
	&& index($script_text, q{my $stagedResidualInput = -s $stagedPlan && (@stagedFiles || $hasShardHandoff);}) >= 0
	&& index($script_text, q{runStagedStrainShardHelper(}) >= 0
	&& index($script_text, q{'-publishedDir', $output}) >= 0
	&& index($script_text, q{$hasShardHandoff = -s $shardManifest ? 1 : 0;}) >= 0
	&& index($script_text, q{runStagedStrainShardHelper('cleanup'}) >= 0
	&& index($script_text, q{unless (@missing || $stagedPrimaryInput || $stagedResidualInput)}) >= 0,
	q{a staged plan resumes aggregate or worker-shard publication through the external helper});
ok(index($script_text, q{die "Staged sample QC input is missing: $rawSampleQC\n" unless fileGZs($rawSampleQC);}) >= 0,
	'tree-side preparation will not mark a missing sample-QC sidecar complete');
unlike($script_text, qr/use Mods::geneCat|readGene2tax|catalogProteins/,
	'buildTree5 staged finalization does not load gene-catalogue indexes');


my $temporary = tempdir(CLEANUP => 1);
my $iqtree_prefix = File::Spec->catfile($temporary, 'IQtree_allsites.backbone');
sub write_test_file {
	my ($path, $contents) = @_;
	open my $handle, '>', $path or die "Cannot write $path: $!";
	print {$handle} $contents;
	close $handle or die "Cannot close $path: $!";
}
my $fake_pigz = File::Spec->catfile($temporary, 'pigz');
write_test_file($fake_pigz, <<'PIGZ');
#!/bin/sh
input=''
for argument in "$@"; do
	input="$argument"
done
cat "$input"
PIGZ
chmod 0755, $fake_pigz or die "Cannot mark $fake_pigz executable: $!";
my $msa_directory = File::Spec->catdir($temporary, 'MSA');
my $msa_cleaned = File::Spec->catdir($msa_directory, 'clnd');
my $msa_work_directory = File::Spec->catdir($temporary, 'MSA-work');
mkdir $msa_directory or die "Cannot create $msa_directory: $!";
mkdir $msa_cleaned or die "Cannot create $msa_cleaned: $!";
mkdir $msa_work_directory or die "Cannot create $msa_work_directory: $!";
write_test_file(File::Spec->catfile($msa_directory, 'COG0001.0.fna'), "gene locus\n");
write_test_file(File::Spec->catfile($msa_directory, 'COG0001.0.faa.gz'), "gene locus\n");
write_test_file(File::Spec->catfile($msa_cleaned, 'COG0001.fna'), "cleaned locus\n");
my $placement_work = File::Spec->catfile($msa_work_directory, 'MSAli.placement.fna');
write_test_file($placement_work, "placement source\n");
my $placement_link = File::Spec->catfile($msa_directory, 'MSAli.placement.fna');
my $retained_alignment = File::Spec->catfile($msa_directory, 'MSAli.fna');
my $retained_work = File::Spec->catfile($msa_work_directory, 'MSAli.fna');
write_test_file($retained_work, "concatenated alignment\n");
write_test_file("$retained_work.nxs", "NEXUS alignment\n");
write_test_file("$retained_work.partition.RAXML", "DNA, part1 = 1-10\n");
my $already_compressed = File::Spec->catfile($msa_directory, 'MSAli.syn.fna.gz');
write_test_file($already_compressed, "retained compressed alignment\n");
{
	no warnings 'once';
	$TestBuildTreeMSAFinalizer::pigzBin = $fake_pigz;
	$TestBuildTreeMSAFinalizer::ncore = 1;
	$TestBuildTreeMSAFinalizer::removeMSA = 1;
	$TestBuildTreeMSAFinalizer::partiExt = '.partition.RAXML';
}
is(TestBuildTreeMSAFinalizer::finalizeMSAArtifacts(
	$msa_directory, $msa_work_directory), 2,
	'finalization publishes both scratch MSAli alignments as compressed checkpoints');
ok(!-e File::Spec->catfile($msa_directory, 'COG0001.0.fna')
		&& !-e File::Spec->catfile($msa_directory, 'COG0001.0.faa.gz'),
	'finalization removes root-level single-locus nucleotide and amino-acid alignments');
ok(!-d $msa_cleaned, 'finalization removes legacy persistent cleaned-locus artifacts');
ok(!-e $retained_alignment && -s "$retained_alignment.gz"
		&& !-e $placement_link && -s "$placement_link.gz"
		&& -s $already_compressed,
	'finalization leaves only durable compressed MSAli FASTA artifacts in persistent storage');
ok(-s "$retained_alignment.nxs.gz"
		&& !-e "$retained_alignment.partition.RAXML"
		&& -s "$retained_alignment.partition.RAXML.gz",
	'NEXUS data and the partition checkpoint are stored only in compressed form');
ok(-s $retained_work && -s $placement_work,
	'publishing and finalizing do not consume the active scratch alignments');
my $restored_directory = File::Spec->catdir($temporary, 'MSA-restored');
mkdir $restored_directory or die "Cannot create $restored_directory: $!";
my $restored_alignment = File::Spec->catfile($restored_directory, 'MSAli.fna');
is(TestBuildTreeMSAFinalizer::restoreCompressedMSAArtifact(
	$retained_alignment, $restored_alignment), $restored_alignment,
	'continuation restores the retained concatenated alignment into fresh scratch');
ok(-s $restored_alignment && -s "$retained_alignment.gz"
		&& !-e $retained_alignment && !-l $retained_alignment,
	'restoration preserves the gzip checkpoint without exposing a persistent plain path');

my $retained_locus_directory = File::Spec->catdir($temporary, 'MSA-retained-loci');
mkdir $retained_locus_directory or die "Cannot create $retained_locus_directory: $!";
my $retained_locus = File::Spec->catfile($retained_locus_directory, 'COG0002.0.fna');
my $discarded_protein = File::Spec->catfile($retained_locus_directory, 'COG0002.0.faa');
my $retained_locus_alignment = File::Spec->catfile($retained_locus_directory, 'MSAli.fna');
write_test_file($retained_locus, "nucleotide locus\n");
write_test_file($discarded_protein, "protein locus\n");
write_test_file($retained_locus_alignment, "concatenated alignment\n");
$TestBuildTreeMSAFinalizer::removeMSA = 0;
is(TestBuildTreeMSAFinalizer::finalizeMSAArtifacts($retained_locus_directory), 1,
	'MSA finalization still compresses the retained concatenated alignment in retention mode');
ok(!-e $retained_locus && -s "$retained_locus.gz" && !-e $discarded_protein,
	'-rmMSA 0 retains nucleotide sub-alignments as gzip checkpoints for later '
		.'population-genetics work while removing unneeded protein alignments');
like($script_text,
	qr/my \$finOutMSA = File::Spec->catfile\(\$MsaWorkD.*?if \(\$endFileExists\).*?restoreCompressedMSAArtifact\(\$publishedOutMSA, \$finOutMSA\).*?if \(!\$endFileExists\).*?runMSAFix\(\$tmpOutMSA, \$maxGapPerCol\).*?publishCompressedMSAArtifact\(\$finOutMSA, \$publishedOutMSA\)/s,
	'per-locus continuation restores gzip checkpoints to scratch, while fresh '
		.'nucleotide MSAs run MSAfix before their first compressed publication');
like($script_text,
	qr/my \$multAli = File::Spec->catfile\(\$MsaWorkD, 'MSAli\.fna'\).*?publishMSAArtifactSet\(\$multAli, \$multAliArtifact\)/s,
	'the concatenated plain MSA is created on scratch and published as a durable artifact');
unlike($script_text, qr/\$pigzBin[^;\n]*\s-d(?:\s|["'])/,
	'buildTree5 no longer decompresses retained MSA checkpoints in place');
unlike($script_text, qr/msaArtifactByWork|exposeMSAWorkingArtifact|copyMSASidecarAtomically/,
	'the scratch lifecycle does not retain path maps, persistent aliases, or a separate sidecar copier');
my $publication_staging = File::Spec->catdir($temporary, 'publication-resume');
like($script_text, qr/sub createTreeOpt\{.*?my \$partiF=\$multF\.\$partiExt;.*?inMSA => \$multF/s,
	'tree programs receive the restored scratch alignment and scratch partition paths directly');

mkdir $publication_staging or die "Cannot create $publication_staging: $!";
write_test_file(File::Spec->catfile($publication_staging, '.strain_tree_input.plan.tsv'),
	"strain-staged-input-v1\noutgroup\tMGS.2643\n");
write_test_file(File::Spec->catfile($publication_staging, 'data.log'), "OG:MGS.2643\n");
write_test_file(File::Spec->catfile($publication_staging, 'sampleQC.tsv'), "MGS\tsample\n");
write_test_file(File::Spec->catfile($publication_staging, 'sampleQC.tsv.tmp'), "incomplete\n");
write_test_file(File::Spec->catfile($publication_staging, 'allFNAs.fna.rewrite.123'), "partial\n");
write_test_file(File::Spec->catfile($publication_staging, 'allFAAs.faa.sort.124'), "partial\n");
write_test_file(File::Spec->catfile($publication_staging, 'allGeneCats.cat.write.125'), "partial\n");
write_test_file(File::Spec->catfile($publication_staging, 'allFNAs.fna.0'), "worker shard\n");
write_test_file(File::Spec->catfile($publication_staging, 'all.cat.tmp.0'), "worker shard\n");
is_deeply(
	[map { (File::Spec->splitpath($_))[2] }
		TestBuildTreeStagedFiles::stagedTreeInputFiles($publication_staging)],
	['data.log', 'sampleQC.tsv'],
	'publication resume detects auxiliary staged files while ignoring unfinished sidecars and worker shards');
my $staged_directory = File::Spec->catdir($temporary, 'staged-strain-inputs');
mkdir $staged_directory or die "Cannot create $staged_directory: $!";
my $raw_category = File::Spec->catfile($staged_directory, 'all.cat.tmp');
my $overlay_category = File::Spec->catfile($staged_directory, '.strain_tree_input.outgroup.cat.tsv');
my $final_category = File::Spec->catfile($staged_directory, 'all.cat');
write_test_file($raw_category,
	"MGS.2\tMGS.2|COG2|gene2\tsampleB\tsampleB|COG2|gene2\n"
	. "MGS.2\tMGS.2|COG1|gene1\tsampleC\tsampleC|COG1|gene1\n"
	. "MGS.2\tMGS.2|COG1|gene1\tsampleA\tsampleA|COG1|gene1\n");
write_test_file($overlay_category,
	"MGS.2|COG1|gene1\tMGS.2643\tMGS.2643|COG1|gene1\n");
is_deeply(
	[TestBuildTreeStagedInputs::finalizeStagedStrainCategory(
		$raw_category, $overlay_category, $final_category, 'MGS.2')],
	[2, 4],
	'tree-side category finalization groups raw Stage-I rows and applies the compact outgroup overlay');
open my $final_category_fh, '<', $final_category or die "Cannot read $final_category: $!";
my $final_category_text = do { local $/; <$final_category_fh> };
close $final_category_fh;
is($final_category_text,
	"MGS.2643|COG1|gene1\tsampleA|COG1|gene1\tsampleC|COG1|gene1\n"
	. "sampleB|COG2|gene2\n",
	'tree-side category finalization writes the deterministic locus/sample format consumed by buildTree5');
my $mixed_category = File::Spec->catfile($staged_directory, 'mixed.cat.tmp');
write_test_file($mixed_category,
	"MGS.2\tMGS.2|COG1|gene1\tsampleA\tsampleA|COG1|gene1\n"
	. "MGS.3\tMGS.3|COG1|gene1\tsampleB\tsampleB|COG1|gene1\n");
my $mixed_category_ok = eval {
	TestBuildTreeStagedInputs::finalizeStagedStrainCategory(
		$mixed_category, $overlay_category, $final_category, 'MGS.2');
	1;
};
ok(!$mixed_category_ok, 'tree-side category finalization rejects mixed-MGS staged input');
like($@, qr/mixes MGS MGS\.2 and MGS\.3/,
	'mixed-MGS category failure identifies both conflicting MGS values');
my $wrong_category = File::Spec->catfile($staged_directory, 'wrong.cat.tmp');
write_test_file($wrong_category,
	"MGS.3\tMGS.3|COG1|gene1\tsampleB\tsampleB|COG1|gene1\n");
my $wrong_category_ok = eval {
	TestBuildTreeStagedInputs::finalizeStagedStrainCategory(
		$wrong_category, $overlay_category, $final_category, 'MGS.2');
	1;
};
ok(!$wrong_category_ok, 'tree-side category finalization rejects an MGS that differs from the plan');
like($@, qr/belongs to MGS\.3, expected MGS\.2/,
	'category MGS mismatch identifies the planned and observed MGS values');
my $raw_qc = File::Spec->catfile($staged_directory, 'sampleQC.tsv.tmp');
my $final_qc = File::Spec->catfile($staged_directory, 'sampleQC.tsv');
write_test_file($raw_qc,
	"MGS.2\tsampleA\tbackbone\t0.10\t0.20\t7\n"
	. "MGS.2\tsampleA\tplacement\t0.30\t0.10\t9\n"
	. "MGS.2\tsampleB\tbackbone\t0.01\t0.02\t8\n");
is(TestBuildTreeStagedInputs::finalizeStagedSampleQC($raw_qc, $final_qc, 'MGS.2'), 2,
	'tree-side QC finalization deduplicates the raw Stage-I sample rows');
open my $final_qc_fh, '<', $final_qc or die "Cannot read $final_qc: $!";
my $final_qc_text = do { local $/; <$final_qc_fh> };
close $final_qc_fh;
is($final_qc_text,
	"MGS\tsample\tstatus\tambiguous_fraction\tcsp_fraction\tvalidated_loci\n"
	. "MGS.2\tsampleA\tplacement\t0.30\t0.2\t9\n"
	. "MGS.2\tsampleB\tbackbone\t0.01\t0.02\t8\n",
	'tree-side QC finalization preserves the controller aggregation policy without catalogue data');
my $mixed_qc = File::Spec->catfile($staged_directory, 'mixed.sampleQC.tsv.tmp');
write_test_file($mixed_qc,
	"MGS.2\tsampleA\tbackbone\t0.1\t0.2\t7\n"
	. "MGS.3\tsampleB\tbackbone\t0.1\t0.2\t7\n");
my $mixed_qc_ok = eval {
	TestBuildTreeStagedInputs::finalizeStagedSampleQC($mixed_qc, $final_qc, 'MGS.2');
	1;
};
ok(!$mixed_qc_ok, 'tree-side QC finalization rejects mixed-MGS staged input');
like($@, qr/mixes MGS MGS\.2 and MGS\.3/,
	'mixed-MGS QC failure identifies both conflicting MGS values');
my $wrong_qc = File::Spec->catfile($staged_directory, 'wrong.sampleQC.tsv.tmp');
write_test_file($wrong_qc,
	"MGS.3\tsampleB\tbackbone\t0.1\t0.2\t7\n");
my $wrong_qc_ok = eval {
	TestBuildTreeStagedInputs::finalizeStagedSampleQC($wrong_qc, $final_qc, 'MGS.2');
	1;
};
ok(!$wrong_qc_ok, 'tree-side QC finalization rejects an MGS that differs from the plan');
like($@, qr/belongs to MGS\.3, expected MGS\.2/,
	'QC MGS mismatch identifies the planned and observed MGS values');
write_test_file("$iqtree_prefix.iqtree", <<'IQTREE');
Best-fit model according to BIC: GTR+F+G2
Rate parameter R:
  A-C: 0.9
  A-G: 3.1
  A-T: 1.2
  C-G: 0.8
  C-T: 4.2
  G-T: 1.0
State frequencies: (empirical counts from alignment)
  pi(A) = 0.29
  pi(C) = 0.21
  pi(G) = 0.23
  pi(T) = 0.27
Gamma shape alpha: 0.73
IQTREE
is(TestBuildTreeEpaHelpers::iqtreePlacementModel($iqtree_prefix),
	'GTR{0.9/3.1/1.2/0.8/4.2/1}+FU{0.29/0.21/0.23/0.27}+G2{0.73}',
	'IQ-TREE fitted GTR rates, base frequencies, and gamma shape are passed explicitly');

my $invariantReport = <<'IQTREE';
Rate parameter R:
A-C: 1
A-G: 2
A-T: 3
C-G: 4
C-T: 5
G-T: 1
pi(A) = 0.1
pi(C) = 0.2
pi(G) = 0.3
pi(T) = 0.4
Proportion of invariable sites: 0.15
Gamma shape alpha: 0.6
IQTREE
is(TestBuildTreeEpaHelpers::iqtreeExplicitEpaModel(
	'GTR+F+I+G4', $invariantReport),
	'GTR{1/2/3/4/5/1}+FU{0.1/0.2/0.3/0.4}+I{0.15}+G4{0.6}',
	'invariant-site proportion and gamma categories are retained in the explicit descriptor');

write_test_file("$iqtree_prefix.iqtree",
	"Model of substitution: GTR+F+G2\n");
write_test_file("$iqtree_prefix.log", <<'IQTREE_LOG');
Rate parameter R:
A-C: 0.7
A-G: 2.8
A-T: 1.1
C-G: 0.6
C-T: 3.9
G-T: 1
pi(A) = 0.22
pi(C) = 0.28
pi(G) = 0.31
pi(T) = 0.19
Gamma shape alpha: 0.44
IQTREE_LOG
is(TestBuildTreeEpaHelpers::iqtreePlacementModel($iqtree_prefix),
	'GTR{0.7/2.8/1.1/0.6/3.9/1}+FU{0.22/0.28/0.31/0.19}+G2{0.44}',
	'IQ-TREE model labels and fitted parameters can be combined across report artifacts');

my $compactReport = <<'IQTREE_COMPACT';
Substitution rates (ML): 1 2 3 4 5 1
Base frequencies (empirical): 0.1 0.2 0.3 0.4
Gamma shape alpha: 0.5
IQTREE_COMPACT
is(TestBuildTreeEpaHelpers::iqtreeExplicitEpaModel(
	'GTR+F+G2', $compactReport),
	'GTR{1/2/3/4/5/1}+FU{0.1/0.2/0.3/0.4}+G2{0.5}',
	'compact IQ-TREE rate and frequency vectors are parsed in documented state order');
my $iqtree3SinglePartition = <<'IQTREE3';
SUBSTITUTION PROCESS
--------------------
  ID  Model           Speed  Parameters
   1  GTR+F+G2       1.0000  GTR{1.65448,9.12038,1.58968,0.848991,17.5795}+F{0.259129,0.215778,0.267947,0.257146}+G2{0.0201622}
IQTREE3
write_test_file("$iqtree_prefix.iqtree", $iqtree3SinglePartition);
write_test_file("$iqtree_prefix.log", '');
is(TestBuildTreeEpaHelpers::iqtreePlacementModel($iqtree_prefix),
	'GTR{1.65448/9.12038/1.58968/0.848991/17.5795/1}+FU{0.259129/0.215778/0.267947/0.257146}+G2{0.0201622}',
	'IQ-TREE 3 compact single-partition GTR table is converted to an EPA-ng descriptor');
is(TestBuildTreeEpaHelpers::iqtreeGtrPartitionCount($iqtree_prefix), 1,
	'IQ-TREE 3 compact single-partition report is identified as one fitted GTR model');

my @iqtree3PartitionRows = map {
	sprintf(' %2d  GTR+F+G2       1.0000  GTR{1.65448,9.12038,1.58968,0.848991,17.5795}+F{0.259129,0.215778,0.267947,0.257146}+G2{0.0201622}', $_)
} 1 .. 8;
write_test_file("$iqtree_prefix.iqtree", join("\n",
	'SUBSTITUTION PROCESS',
	'  ID  Model           Speed  Parameters',
	@iqtree3PartitionRows,
	'',
));
write_test_file("$iqtree_prefix.log",
	'Command: iqtree3 -s alignment.fna -m GTR+F+G2 -T 12');
my $partitioned_model_warning = '';
{
	local $SIG{__WARN__} = sub { $partitioned_model_warning .= $_[0] };
	is(TestBuildTreeEpaHelpers::iqtreePlacementModel($iqtree_prefix), 'GTR+F+G2',
		'partitioned IQ-TREE 3 report retains the generic EPA-ng model rather than one parameter row');
}
like($partitioned_model_warning,
	qr/reported 8 fitted GTR parameter sets for separate partitions/,
	'partitioned IQ-TREE 3 fallback reports the actual number of fitted model rows');
is(TestBuildTreeEpaHelpers::iqtreeGtrPartitionCount($iqtree_prefix), 8,
	'IQ-TREE 3 compact partitioned report is identified for an unpartitioned EPA-ng refit');
like($script_text,
	qr/sub epaModelArtifact \{.*?return epaRefitIqtreeModel\(.*?\n\t\t\tif iqtreeGtrPartitionCount\(.*?\) > 1;.*?sub epaRefitIqtreeModel \{.*?\$refitOpts\{partition\} = '';.*?\$refitOpts\{fixedTree\} = \$backboneTree;.*?my \$model = iqtreePlacementModel\(\$refitPrefix\);/s,
	'partitioned IQ-TREE GTR backbones are refit unpartitioned on their fixed topology before EPA-ng');
like($partitioned_model_warning, qr/will not silently use one partition's rates/,
	'partitioned IQ-TREE 3 fallback makes the EPA-ng single-model limitation explicit');

write_test_file("$iqtree_prefix.iqtree",
	"Model of substitution: GTR+F+G2\n");
write_test_file("$iqtree_prefix.log", "IQ-TREE log without fitted parameters\n");
my $incomplete_model_warning = '';
{
	local $SIG{__WARN__} = sub { $incomplete_model_warning .= $_[0] };
	is(TestBuildTreeEpaHelpers::iqtreePlacementModel($iqtree_prefix),
		'GTR+F+G2',
		'an incomplete fitted GTR report falls back to the generic descriptor');
}
like($incomplete_model_warning,
	qr/using the generic symbolic descriptor instead\. EPA-ng does not refit missing parameters/,
	'GTR parsing failure is reported as a visible warning');

write_test_file("$iqtree_prefix.iqtree", "IQ-TREE report without a model label\n");
write_test_file("$iqtree_prefix.log",
	"Command: iqtree3 -s alignment.fna -m HKY+F+G4 -T 12\n");
my $generic_model_warning = '';
{
	local $SIG{__WARN__} = sub { $generic_model_warning .= $_[0] };
	is(TestBuildTreeEpaHelpers::iqtreePlacementModel($iqtree_prefix), 'HKY+F+G4',
		'IQ-TREE command-line model is the final parser fallback');
}
like($generic_model_warning, qr/could not serialize fitted IQ-TREE parameters/,
	'non-GTR generic model fallback is also reported as a warning');
is_deeply(
	[TestBuildTreeEpaHelpers::epaResourcePlan(12, 12, -1, 4500)],
	[2, 2700],
	'EPA threads are capped to one per GB of the derived planning budget');
is_deeply(
	[TestBuildTreeEpaHelpers::epaResourcePlan(4, 12, 0, 4500)],
	[4, 0],
	'zero planning memory keeps the requested core-capped thread count');

my $classification = File::Spec->catfile($temporary, 'strict_backbone.samples.tsv');
write_test_file($classification, join("\n",
	join("\t", qw(sample tree_role reason informative_positions q90_informative
		backbone_overlap_nt backbone_overlap_loci backbone_state_divergence)),
	join("\t", qw(backbone1 backbone validated 1000 900 1000 8 0)),
	join("\t", qw(query1 placement sparse 450 900 425 4 0.03)),
	join("\t", qw(outlier1 excluded divergent 200 900 100 1 0.8)),
)."\n");
my $classification_state =
	TestBuildTreeEpaClassification::readStrictBackboneClassification($classification);
is_deeply($classification_state->{placement}, ['query1'],
	'EPA-only recovery reuses exactly the samples classified for placement');
is($classification_state->{backbone_overlap}{query1}{backbone_overlap_nt}, 425,
	'EPA-only placement reporting retains the original backbone-overlap metric');
like($script_text,
	qr/print STDERR "EPA-ng command: \$command";.*?systemW\(\$command\)/s,
	'the exact EPA-ng command is written to STDERR before execution');

like($script_text,
	qr/"epaOnly=i" => \\\$epaOnly.*?if \(\$epaOnly\).*?runEpaOnlyPlacement\(.*?exit\(0\)/s,
	'EPA-only mode exits through its dedicated placement path before ordinary MSA and inference work');
like($script_text,
	qr/sub runEpaOnlyPlacement.*?requires a validated IQ-TREE backbone.*?map_epa_placements_to_backbone\(.*?write_epa_placed_tree\(\$backboneTreeText, \$primaryTree.*?backbone retained=\$backboneTree/s,
	'EPA-only mode maps jplace edges and grafts placements onto the retained backbone');
like($script_text,
	qr/if \(\$subsetSmpls >0\).*?if \(\$redoEPAfilter\).*?runRedoEpaFilter\(.*?exit\(0\).*?warn "MSAprobs.*?prepGenoDirs/s,
	'forced EPA filtering exits before sequence inputs, alignment, and inference startup');
like($script_text, qr/my \$redoEPAfilter =\s*\(\$ENV\{MATAFILER_REDO_EPA_FILTER\}/,
	'BuildTree inherits forced filtering when an older saved command is resubmitted');
my ($redo_epa_body) = $script_text =~
	/(sub runRedoEpaFilter .*?)(?=sub readEpaFilterBackboneTree)/s;
ok(defined($redo_epa_body), 'focused forced-EPA publication helper is available');
like($redo_epa_body,
	qr/readStrictBackboneClassification.*?read_epa_jplace.*?map_epa_placements_to_backbone.*?filter_epa_placement_outliers.*?write_epa_placed_tree.*?writeCompletionMarker/s,
	'forced EPA filtering reads only retained publication artifacts and republishes lifecycle state');
unlike($redo_epa_body, qr/runEpaNgPlacement|prepGenoDirs|mergeMSAs|treeAtHeart/,
	'forced EPA filtering cannot start EPA-ng, alignment, or tree inference');
my $redo_output = File::Spec->catdir($temporary, 'redo_output');
my $redo_phylo = File::Spec->catdir($redo_output, 'phylo');
my $redo_epa = File::Spec->catdir($redo_phylo, 'epa-ng');
mkdir $redo_output or die "Cannot create $redo_output: $!";
mkdir $redo_phylo or die "Cannot create $redo_phylo: $!";
mkdir $redo_epa or die "Cannot create $redo_epa: $!";
write_test_file(
	File::Spec->catfile($redo_phylo, 'IQtree_allsites.backbone.treefile'),
	"(A:0.1,B:0.1,C:0.1);\n",
);
write_test_file(
	File::Spec->catfile($redo_phylo, 'strict_backbone.samples.tsv'),
	join("\n",
		join("\t", qw(sample tree_role reason informative_positions
			q90_informative backbone_overlap_nt backbone_overlap_loci
			backbone_state_divergence)),
		(map { join("\t", $_, 'backbone', 'validated', 1000, 900, 1000, 8, 0) }
			qw(A B C)),
		join("\t", qw(query1 placement sparse 450 900 425 4 0.03)),
	)."\n",
);
write_test_file(
	File::Spec->catfile($redo_epa, 'epa_result.jplace'),
	'{"tree":"(A:0.1{0},B:0.1{1},C:0.1{2});",'
	.'"placements":[{"p":[[0,-10,1,0.05,0.01]],"n":["query1"]}],'
	.'"metadata":{},"version":3,'
	.'"fields":["edge_num","likelihood","like_weight_ratio",'
	.'"distal_length","pendant_length"]}',
);
my $redo_completion = File::Spec->catfile($redo_output, 'treeDone.sto');
my $redo_wrapper = File::Spec->catfile($temporary, 'run-redo-buildtree.pl');
write_test_file($redo_wrapper, <<'PERL');
use strict;
use warnings;
use Mods::IO_Tamoc_progs qw(setConfigFile);
my $root = shift @ARGV;
my $config = shift @ARGV;
my $script = shift @ARGV;
$ENV{MF4_TEST_ROOT} = $root;
setConfigFile($config);
my $result = do $script;
die $@ if $@;
die "Cannot execute $script: $!\n" unless defined $result;
PERL
is(system(
	$^X, "-I$root", $redo_wrapper, $root,
	File::Spec->catfile($Bin, 'MATAFILERcfg.txt'), $script,
	'-fna', File::Spec->catfile($redo_output, 'deliberately_missing.fna'),
	'-aa', File::Spec->catfile($redo_output, 'deliberately_missing.faa'),
	'-cats', File::Spec->catfile($redo_output, 'deliberately_missing.cat'),
	'-outD', $redo_output, '-runIQtree', 1, '-cores', 1,
	'-withinSpecies', 1, '-strictBackbone', 1,
	'-rateMergePartitions', 0, '-continue', 1, '-redoEPAfilter', 1,
	'-completionMarker', $redo_completion,
), 0, 'forced EPA filtering succeeds without opening missing sequence inputs');
my $redo_primary =
	File::Spec->catfile($redo_phylo, 'IQtree_allsites.treefile');
ok(-s $redo_primary && -s $redo_completion,
	'forced EPA filtering republishes the primary tree and completion marker');
open my $redo_tree_handle, '<', $redo_primary
	or die "Cannot read $redo_primary: $!";
my $redo_tree_text = do { local $/; <$redo_tree_handle> };
close $redo_tree_handle or die "Cannot close $redo_primary: $!";
like($redo_tree_text, qr/query1/,
	'forced EPA filtering grafts the retained query before exiting');


my $coordinate_bounds_checks = () = $script_text =~ /next if \$position >= length\(\$sequence\);/g;
cmp_ok($coordinate_bounds_checks, '>=', 2,
	'taxon-aware raw and alignment coordinate scorers safely skip uneven sequence tails');
like($script_text,
	qr/if \(defined\(\$candidateSelection->\{terminal_reason\}\).*?writeOutcomeMarker\(\$terminalMarker, 'valid_no_tree', \$reason.*?exit\(0\)/s,
	'the candidate-selection terminal result is published as a persistent valid no-tree outcome');
like($script_text,
	qr/unless \(keys %metrics\).*?terminal_reason => 'taxon_aware_no_category_with_three_usable_samples'/s,
	'a taxon-aware candidate set with fewer than three usable samples returns a stable terminal reason');
my $placement_outlier_calls = () = $script_text =~ /filter_epa_placement_outliers\(/g;
cmp_ok($placement_outlier_calls, '>=', 3,
	'fresh, EPA-only, and forced-redo publication apply pendant-branch outlier QC');
my $backbone_mapping_calls =
	() = $script_text =~ /map_epa_placements_to_backbone\(/g;
cmp_ok($backbone_mapping_calls, '>=', 3,
	'normal, EPA-only, and forced-redo publication map jplace edges onto the backbone');
like($script_text, qr/strict_backbone\.epa_backbone_grafts\.tsv/,
	'EPA backbone grafting publishes a per-edge mapping report');
unlike($script_text, qr/write_epa_placed_tree\(\$epaResult->\{tree\}/,
	'the jplace Newick tree is never used as the publication template');
like($script_text,
	qr/pendant_outlier_limit placement_filter_reason reason/s,
	'EPA placement reports expose the applied cutoff and exclusion reason');
is(system($^X, q{-I}.$root, q{-c}, $script), 0, q{buildTree5.pl compiles});

done_testing;
