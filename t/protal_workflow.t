use strict;
use warnings;

use File::Spec;
use FindBin qw($Bin);
use Test::More;

my $root = File::Spec->catdir($Bin, '..');
sub slurp {
	my ($path) = @_;
	open my $fh, '<', $path or die "Cannot read $path: $!";
	my $text = do { local $/; <$fh> };
	close $fh or die "Cannot close $path: $!";
	return $text;
}

my $source = slurp(File::Spec->catfile($root, 'MATAF4.pl'));
like($source, qr/"profileProtal=i"\s*=>\s*\\\$MFopt\{DoProtal\}/,
	'profileProtal is a parsed MATAF4 option');
like($source, qr/"protalIgnoreErrors=i"\s*=>\s*\\\$MFopt\{protalIgnoreErrors\}/,
	'protalIgnoreErrors is a parsed MATAF4 option');
like($source, qr/\$MFopt\{DoProtal\}=0.*?\$MFopt\{ProtalCores\}=4;\s*\$MFopt\{ProtalMem\}=100/s,
	'Protal defaults include explicit CPU and memory resources');
like($source, qr/\$MFopt\{protalIgnoreErrors\}=1/,
	'tolerant Protal input selection is enabled by default');
like($source, qr/DoProtal protalIgnoreErrors DoTaxaTarget/,
	'the tolerant-input setting participates in sample completion signatures');
like($source, qr/-protalIgnoreErrors must be 0 or 1/,
	'the tolerant-input setting is constrained to a boolean value');
like($source, qr/redoFails.*?\(\$MFopt\{DoProtal\} == 1 && \$calcProtal && !\$MFopt\{protalIgnoreErrors\}\)/s,
	'only strict singular Protal failure participates in per-sample redo cleanup');
like($source, qr/\$configuration\{protal_database\}\s*=\s*\$protalDB/,
	'the Protal database identity participates in sample completion signatures');
like($source, qr/\$requireRawReadsFlag\s*=\s*1\s+if\s*\([^;]*\$calcProtal/s,
	'raw input cleanup is held while a Protal profile is pending');
like($source, qr/# Protal deliberately consumes.*?protalMapping\(.*?\$jdep.*?sdmClean\(/s,
	'Protal is submitted from staged raw input before SDM processing');
like($source, qr/\$protalProfileJobs\{\$profile\}\s*=\s*\$protalJob.*?add2SampleDeps\(\\\@sampleDeps,\s*\[\$protalJob\]\)/s,
	'per-sample Protal jobs are retained for both merge and sample cleanup dependencies');

my ($mapping) = $source =~ /sub protalMapping \{(.*?)\n\}\n\nsub mergeProtalProfiles/s;
ok(defined($mapping), 'found the complete protalMapping implementation');
my ($selector) = $source =~ /sub protalRawPairSelection \{(.*?)\n\}\n\nsub registerCombinedProtalSample/s;
ok(defined($selector), 'found the shared Protal raw-pair selector');
like($selector, qr/sampleReadSet\(\$sampleKey, 'raw'\).*?readLibrariesByScope\(\s*\$rawReadSet, 'primary', 0, \$sampleKey\).*?singleShortReadPair/s,
	'both Protal modes consume the structured primary raw-read library');
like($selector, qr/singleShortReadPair\(.*?ignore_incompatible\s*=>\s*\$MFopt\{protalIgnoreErrors\}/s,
	'both modes delegate tolerant first-pair selection to the existing read helper');
like($mapping, qr/protalRawPairSelection\(\s*\$curSmpl/s,
	'singular Protal reuses the common raw-pair selector');
like($mapping, qr/\$selection->\{skipped\}.*?atomic_write_text\(\$skip, "\$requestSignature\\n\$reason\\n".*?return ''/s,
	'an incompatible sample receives a durable request-scoped skip marker');
like($mapping, qr/\$selection->\{warning\}.*?retry_unlink\(\$skip.*?make_path\(\$profileDir\)/s,
	'fallback selection warns and clears stale skip evidence before normal profiling');
unlike($mapping, qr/sampleReadSet\([^\n]*['"]clean['"]/,
	'Protal does not consume the cleaned-read set');
like($mapping, qr/my \$protal = getProgPaths\('protal'\)/,
	'the Protal executable is resolved through shared configuration');
like($mapping, qr/'-1', \$read1, '-2', \$read2.*?'--outdir', \$tmpD.*?'--profile_dir', \$profileDir.*?'--no_strains'/s,
	'Protal receives one raw pair, scratch output, durable profile output, and no strain request');
like($mapping, qr/getProgPaths\('protal_db', 0\).*?PROTAL_DB_PATH/s,
	'Protal supports configured database paths and the documented environment fallback');
like($mapping, qr/qw\(\.log \.gene\.log \.genes\.log \.truth_annotated\).*?_shell_command\('rm', '-f', '--', \@diagnostics\)/s,
	'profile diagnostic files are removed after successful profiling');
like($mapping, qr/_shell_command\('rm', '-rf', '--', \$tmpD\).*?_shell_command\('touch', \$stone\)/s,
	'node-local alignment output is removed before the durable completion stone is written');
like($mapping, qr/_shell_command\('test', '-e', \$profile\)/,
	'a successful zero-detection profile is accepted when the file exists');

my ($merge) = $source =~ /sub mergeProtalProfiles \{(.*?)\n\}\n\nsub mergeMP2Table/s;
ok(defined($merge), 'found the complete mergeProtalProfiles implementation');
like($merge, qr/\$selectedFrom < \$selectedTo.*?sample_is_ignored.*?SMPL\.empty/s,
	'the merge cohort follows the selected range and excludes intentionally skipped samples');
like($merge, qr/sampleCompletionRequestSignature.*?protalSkipMatches\(/s,
	'the merge excludes only incompatible samples skipped for the current request');
like($merge, qr/neither a complete profile nor a .*?submitted Protal job.*?return;/s,
	'the merge is deferred instead of silently emitting a partial cohort');
like($merge, qr/\$controllerBaseOut, 'pseudoGC', 'protal_singular'/,
	'the singular merged table uses its requested stable controller output root');
like($merge, qr/getProgPaths\('protalProfileUtils'\).*?_shell_command\(\$profileUtils, 'merge', '--input', \@profiles\)/s,
	'the final job invokes protal_profile_utils merge with the selected profiles');
like($merge, qr/normalise_job_dependencies\(\\\@dependencies\)/,
	'the final merge is scheduler-dependent on the per-sample Protal jobs');
like($merge, qr/Protal\.abundance\.tsv.*?_shell_command\('mv', '-f', '--', \$temporary, \$merged\)/s,
	'the merged abundance table is atomically published from a temporary file');
like($source, qr/my \$useProtalSkip = \$MFopt\{protalIgnoreErrors\}.*?protalSkipMatches.*?incompatible_input_skip/s,
	'sample completion accepts only a signature-valid tolerant skip marker');
like($source, qr/!\(-e \$protalProfile && -e \$protalStone\).*?\@protalChecks/s,
	'a complete Protal profile takes precedence over skip evidence');
like($source, qr/sub postprocess\{.*?mergeProtalProfiles\(\)/s,
	'Protal profile merging is part of normal postprocessing');
my ($combined) = $source =~ /sub submitCombinedProtal \{(.*?)\n\}\n\nsub protalMapping/s;
ok(defined($combined), 'found the complete combined-Protal implementation');
like($source, qr/-profileProtal must be 0 .*? 1 .*? or 2 .*?\$MFopt\{DoProtal\} >= 0.*?<= 2/s,
	'profileProtal accepts only off, singular, and combined modes');
like($source, qr/profileProtal 2 requires the complete mapped cohort.*?\$selectedFrom == 0 && \$selectedTo == \@samples/s,
	'combined mode requires the full mapped cohort');
like($source, qr/sub protalCombinedOutputDir \{.*?\$controllerBaseOut, 'pseudoGC', 'protal'/s,
	'combined results use pseudoGC/protal');
like($combined, qr/my \$profileDir = File::Spec->catdir\(\$outDir, 'profiles'\).*?my \$strainDir = \$protalCombinedStrains/s,
	'profiles and strain MSAs are routed to durable combined-result directories');
like($combined, qr/#SAMPLEID.*?FIRST.*?SECOND.*?SAM.*?PREFIX.*?PROFILE/s,
	'the generated combined map declares Protal input and output columns');
like($combined, qr/getProgPaths\('protal'\).*?getProgPaths\('protalProfileUtils'\).*?getProgPaths\('protal_db', 0\)/s,
	'combined mode resolves all Protal commands and its database through getProgPaths');
like($combined, qr/'--map', \$mapFile.*?'--profile_dir', \$profileDir.*?'--threads'.*?'--force'/s,
	'combined mode invokes Protal once with the generated map');
unlike($combined, qr/--no_strains/,
	'combined mode leaves Protal strain analysis enabled so MSAs are retained');
like($combined, qr/my \$workDir.*?nodeTmpDir.*?#SAM_OUTPUT_DIR.*?alignments.*?#MISC_OUTPUT_DIR.*?misc/s,
	'large alignments and miscellaneous files remain in temporary work space');
like($combined, qr/if \(\$protalDB eq '' && \@profiles\)/,
	'an all-incompatible tolerant cohort can finish without a database');
like($source, qr/my \$combinedProtalHoldsScratch =.*?DoProtal\} == 2 && !\$protalCombinedComplete/s,
	'normal sample cleanup is held until combined Protal completes');
like($combined, qr/pathIsStrictChild\(\$target, \$MFglobal\{runTmpDirGlobal\}\).*?ProtalScratchCleanup\.sh/s,
	'the final cleanup script contains only validated run-scratch children');
like($combined, qr/protalProfileUtils.*?merge.*?_shell_command\('bash', \$cleanupScript\).*?\$protalCombinedCurrent.*?\$protalCombinedStone/s,
	'profiles are merged and scratch is cleaned before combined completion is published');
like($combined, qr/keys %\{\$QSBoptHR->\{submittedJobRecords\}.*?normalise_job_dependencies/s,
	'the final job waits for all jobs submitted by this controller pass');
like($combined, qr/recordSampleLockJobs\(\$lockFile, \[\$job\], \$QSBoptHR\)/,
	'each participating sample lock records the combined job');
like($source, qr/sub protalCombinedRunComplete \{.*?protalCombinedCurrent.*?protalCombinedStrains/s,
	'combined completion checks the active cohort marker and retained MSA directory');

my $internal = slurp(File::Spec->catfile($root, 'Mods', 'config_internal.txt'));
like($internal, qr/^protal\tprotal$/m, 'shared config exposes the Protal executable');
like($internal, qr/^protalProfileUtils\tprotal_profile_utils$/m,
	'shared config exposes the profile merge utility');
my $databases = slurp(File::Spec->catfile($root, 'Mods', 'config_DBs.txt'));
like($databases, qr/^protal_db\s*$/m,
	'database config provides an optional Protal DB override');
my $environment = slurp(File::Spec->catfile($root, 'helpers', 'install', 'MF4.yml'));
like($environment, qr/^\s*-\s+bioconda::protal=0\.6\.0a\s*$/m,
	'the base MF4 environment installs the Protal package and profile utility');
my $installer = slurp(File::Spec->catfile($root, 'helpers', 'install', 'installer.sh'));
like($installer, qr/\$tool" --version.*?&& !.*?\$tool" --help/s,
	'the installer can verify commands which expose help but no version option');
like($installer, qr/verify_environment_tools MF4.*?Protal profiling.*?protal protal_profile_utils/s,
	'the installer verifies both Protal commands in the base environment');
my $installCheck = slurp(File::Spec->catfile($root, 'Mods', 'TamocFunc.pm'));
like($installCheck, qr/checkProg\("Protal profiler","protal",0\).*?checkProg\("Protal profile merger","protalProfileUtils",0\)/s,
	'MATAFILER checkInstall validates both configured Protal commands');

done_testing();
