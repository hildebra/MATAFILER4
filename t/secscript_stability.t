use strict;
use warnings;
use Test::More;
use File::Spec;
use File::Temp qw(tempdir);
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IPC::Open3;
use Symbol qw(gensym);

my $root = File::Spec->rel2abs('.');
my $tmp = tempdir(CLEANUP => 1);

sub write_file {
    my ($path, $contents) = @_;
    open my $fh, '>', $path or die "Cannot open $path: $!";
    print {$fh} $contents or die "Cannot write $path: $!";
    close $fh or die "Cannot close $path: $!";
}

sub read_file {
    my ($path) = @_;
    open my $fh, '<', $path or die "Cannot open $path: $!";
    local $/;
    my $contents = <$fh>;
    close $fh or die "Cannot close $path: $!";
    return $contents;
}

my $rename_input = File::Spec->catfile($tmp, 'rename.fasta');
write_file($rename_input, ">first\nACGT\n>second\nAACCGG\n");
is(system($^X, File::Spec->catfile($root, 'secScripts', 'assemblies', 'renameCtgs.pl'), $rename_input, 'sample'), 0,
   'renameCtgs completes');
my $renamed = read_file($rename_input);
like($renamed, qr/^>sample__C1_L=4=\nACGT\n>sample__C2_L=6=\nAACCGG\n$/,
     'renameCtgs records the correct length for the final contig');

my $filter_input = File::Spec->catfile($tmp, 'filter.fasta');
write_file($filter_input, ">long\nAAAAAAAAAA\n>last-secondary\nCCCCCC\n");
is(system($^X, File::Spec->catfile($root, 'secScripts', 'assemblies', 'sizeFilterFas.pl'), $filter_input, 8, 5), 0,
   'sizeFilterFas completes');
like(read_file("$filter_input.filt"), qr/>long\nAAAAAAAAAA\n/, 'primary output contains long sequence');
like(read_file("$filter_input.filt2"), qr/>last-secondary\nCCCCCC\n/,
     'secondary output contains the final short FASTA record');

my $synthetic_fasta = File::Spec->catfile($tmp, 'synthetic.fasta');
write_file($synthetic_fasta,
		   ">ctgA descriptive header\n" . ('A' x 30_000)
		   . "\n>splitShort\n" . ('G' x 4_000) . "\n>last\n" . ('C' x 1_000) . "\n");
my $synthetic_coverage = File::Spec->catfile($tmp, 'mapping.coverage.gz');
my $coverage_text = "ctgA\t0\t10000\t4\nctgA\t10000\t12000\t0\nctgA\t12000\t30000\t1\n"
	. "splitShort\t0\t1800\t2\nsplitShort\t1800\t2200\t0\nsplitShort\t2200\t4000\t2\n"
	. "last\t0\t1000\t1\n";
gzip(\$coverage_text => $synthetic_coverage)
    or die "Cannot create $synthetic_coverage: $GzipError";
my $synthetic_fastq = File::Spec->catfile($tmp, 'synthetic.fastq.gz');
my $simulator = File::Spec->catfile($root, 'secScripts', 'assemblies', 'split_fasta4metaMDBG.pl');
my $simulator_source = read_file($simulator);
unlike($simulator_source, qr/my %intervals = map/,
	'synthetic-read simulation does not preallocate per-contig interval arrays');
like($simulator_source, qr/push \@\{\$intervals\{\$id\}\}, 0, \$lengths/,
	'coverage intervals use flat storage instead of per-interval array objects');
like($simulator_source, qr/delete\(\$coverage_runs\{\$header\}\)/,
	'coverage state is released as each contig is simulated');
my $breakpoint_detector = File::Spec->catfile($root, 'secScripts', 'assemblies', 'breakpoints.pl');
my $breakpoint_tsv = File::Spec->catfile($tmp, 'breakpoints.tsv.gz');
is(system($^X, '-I' . $root, $breakpoint_detector,
          '--assembly', $synthetic_fasta, '--coverage', $synthetic_coverage,
          '--output', $breakpoint_tsv, '--breakpoint-depth', 0.1,
          '--min-breakpoint-length', 100, '--smooth-gap', 100,
          '--flank-length', 500, '--min-flank-depth', 1,
          '--max-flank-fraction', 0.1), 0,
   'standalone breakpoint detector completes');
my $breakpoint_text = '';
gunzip($breakpoint_tsv => \$breakpoint_text)
    or die "Cannot read $breakpoint_tsv: $GunzipError";
like($breakpoint_text, qr/^ctgA\t10000\t12000\t2000\t/m,
     'breakpoint TSV contains the supported internal low-coverage interval');
my $smooth_fasta = File::Spec->catfile($tmp, 'smooth-break.fasta');
write_file($smooth_fasta, ">smooth\n" . ('A' x 5_000) . "\n");
my $smooth_coverage = File::Spec->catfile($tmp, 'smooth-break.coverage.gz');
my $smooth_coverage_text = join '',
    "smooth\t0\t2000\t5\n", "smooth\t2000\t2500\t0\n",
    "smooth\t2500\t2550\t1\n", "smooth\t2550\t3000\t0\n",
    "smooth\t3000\t4500\t5\n", "smooth\t4500\t5000\t0\n";
gzip(\$smooth_coverage_text => $smooth_coverage)
    or die "Cannot create $smooth_coverage: $GzipError";
my $smooth_tsv = File::Spec->catfile($tmp, 'smooth-breakpoints.tsv.gz');
is(system($^X, '-I' . $root, $breakpoint_detector,
          '--assembly', $smooth_fasta, '--coverage', $smooth_coverage,
          '--output', $smooth_tsv, '--smooth-gap', 100), 0,
   'breakpoint smoothing run completes');
my $smooth_result = '';
gunzip($smooth_tsv => \$smooth_result)
    or die "Cannot read $smooth_tsv: $GunzipError";
like($smooth_result, qr/^smooth\t2000\t3000\t1000\t/m,
     'a short noisy mapping island is smoothed into one supported breakpoint');
unlike($smooth_result, qr/^smooth\t4500\t5000\t/m,
       'a terminal low-coverage run is rejected without support on both sides');

# Breakpoint support is relative to the weaker flank.  A genuinely supported
# zero-depth interval must therefore be treated identically on a low-depth
# contig and a high-depth contig, provided both clear the 1x evidence floor.
my $depth_scale_fasta = File::Spec->catfile($tmp, 'depth-scale.fasta');
write_file($depth_scale_fasta,
    ">low2x\n" . ('A' x 5_000) . "\n>high20x\n" . ('C' x 5_000) . "\n");
my $depth_scale_coverage = File::Spec->catfile($tmp, 'depth-scale.coverage.gz');
my $depth_scale_text = join '',
    "low2x\t0\t2000\t2\n", "low2x\t2000\t3000\t0\n", "low2x\t3000\t5000\t2\n",
    "high20x\t0\t2000\t20\n", "high20x\t2000\t3000\t0\n", "high20x\t3000\t5000\t20\n";
gzip(\$depth_scale_text => $depth_scale_coverage)
    or die "Cannot create $depth_scale_coverage: $GzipError";
my $depth_scale_tsv = File::Spec->catfile($tmp, 'depth-scale.breakpoints.tsv.gz');
is(system($^X, '-I' . $root, $breakpoint_detector,
          '--assembly', $depth_scale_fasta, '--coverage', $depth_scale_coverage,
          '--output', $depth_scale_tsv), 0,
   'breakpoint detector accepts matched low- and high-depth contigs');
my $depth_scale_result = '';
gunzip($depth_scale_tsv => \$depth_scale_result)
    or die "Cannot read $depth_scale_tsv: $GunzipError";
like($depth_scale_result, qr/^low2x\t2000\t3000\t1000\t0\.0000\t2\.0000\t2\.0000$/m,
     'a breakpoint on a 2x contig is retained');
like($depth_scale_result, qr/^high20x\t2000\t3000\t1000\t0\.0000\t20\.0000\t20\.0000$/m,
     'the equivalent breakpoint on a 20x contig is retained');

my $simulator_err = gensym;
my $simulator_pid = open3(undef, my $simulator_out, $simulator_err,
    $^X, '-I' . $root, $simulator,
    '--assembly', $synthetic_fasta, '--coverage', $synthetic_coverage,
    '--breakpoints', $breakpoint_tsv,
    '--output', $synthetic_fastq, '--mean-read-length', 3_000,
    '--read-length-sd', 300, '--seed', 7);
my $simulator_stdout = do { local $/; <$simulator_out> };
my $simulator_stderr = do { local $/; <$simulator_err> };
waitpid($simulator_pid, 0);
is($? >> 8, 0, 'metaMDBG preparation accepts its flag-based interface');
is($simulator_stderr, '', 'metaMDBG simulation emits no warnings for valid input');
like($simulator_stdout, qr/Synthetic read simulation summary/,
     'simulation reports a readable end-of-run summary');
like($simulator_stdout, qr/Breakpoints identified:\s+2 across 2 contig/,
     'simulation summary reports identified breakpoints');
like($simulator_stdout, qr/Simulated reads:\s+19\b/,
     'simulation summary reports its output read count');
my $synthetic_text = '';
gunzip($synthetic_fastq => \$synthetic_text)
    or die "Cannot read $synthetic_fastq: $GunzipError";
my @synthetic_headers = ($synthetic_text =~ /^\@([^\n]+)/mg);
is(scalar @synthetic_headers, 19,
   'coverage integrals determine the number of randomly placed reads');
my @ctga_coordinates = map {
    /^ctgA_SIM_\d+_START_(\d+)_END_(\d+)_ANCHOR_(\d+)$/ ? [$1, $2, $3] : ()
} @synthetic_headers;
is(scalar @ctga_coordinates, 19, 'the two covered ctgA blocks produce their expected reads');
ok(!grep({ !($_->[1] <= 10_000 || $_->[0] >= 12_000) } @ctga_coordinates),
   'no simulated read overlaps or crosses the zero-coverage breakpoint');
is(scalar(grep { $_->[2] < 10_000 } @ctga_coordinates), 13,
   'the four-fold high-coverage block receives proportionally more read anchors');
is(scalar(grep { $_->[2] >= 12_000 } @ctga_coordinates), 6,
   'the lower-coverage block receives proportionally fewer read anchors');
is(scalar(grep { /^(?:last|splitShort)_SIM_/ } @synthetic_headers), 0,
	'short contigs and breakpoint-shortened fragments emit no synthetic reads');
my @synthetic_sequences = ($synthetic_text =~ /^\@[^\n]+\n([^\n]+)\n\+\n/mg);
my @ctga_lengths = map { $_->[1] - $_->[0] } @ctga_coordinates;
ok(scalar(keys %{ { map { $_ => 1 } @ctga_lengths } }) > 3,
   'simulated read lengths vary around the requested mean');
my $total_simulated_length = 0;
$total_simulated_length += $_ for @ctga_lengths;
my $mean_simulated_length = $total_simulated_length / @ctga_lengths;
cmp_ok($mean_simulated_length, '>', 2_600, 'simulated mean length remains near the request');
cmp_ok($mean_simulated_length, '<', 3_400, 'simulated mean length remains near the request');

my $legacy_coverage = File::Spec->catfile($tmp, 'Coverage.median.percontig.gz');
my $legacy_coverage_text = "ctgA\t1\nlast\t1\n";
gzip(\$legacy_coverage_text => $legacy_coverage)
    or die "Cannot create $legacy_coverage: $GzipError";
my $legacy_fastq = File::Spec->catfile($tmp, 'legacy.fastq.gz');
my $empty_breakpoints = File::Spec->catfile($tmp, 'no-breakpoints.tsv');
write_file($empty_breakpoints, "contig\tstart\tend\tlength\tmean_depth\tleft_depth\tright_depth\n");
is(system($^X, '-I' . $root,
          $simulator,
          '--assembly', $synthetic_fasta, '--coverage', $legacy_coverage,
          '--breakpoints', $empty_breakpoints,
          '--output', $legacy_fastq, '--mean-read-length', 3_000,
          '--read-length-sd', 300, '--seed', 7), 0,
   'older hybrid packages retain contig-wide coverage compatibility');
my $legacy_text = '';
gunzip($legacy_fastq => \$legacy_text)
    or die "Cannot read $legacy_fastq: $GunzipError";
is(scalar(() = $legacy_text =~ /^\@/mg), 10,
	'contig-wide fallback skips contigs shorter than the synthetic-read minimum');
my $comparison = File::Spec->catfile($tmp, 'HybridAssemblyComparison.tsv');
is(system($^X, File::Spec->catfile($root, 'secScripts', 'assemblies', 'compare_hybrid_assemblies.pl'),
          '--preassembly', $synthetic_fasta, '--final', $smooth_fasta,
          '--output', $comparison), 0,
   'comparative hybrid assembly report completes');
like(read_file($comparison), qr/^N50\t30000\t5000\t-25000\t/m,
     'comparative report contains preassembly and final N50 values');

my $depth_converter = File::Spec->catfile($root, 'helpers', 'samcovToBedGraph.pl');
my $converter_err = gensym;
my $converter_pid = open3(my $converter_in, my $converter_out, $converter_err, $^X, $depth_converter);
print {$converter_in} "ctg\t1\t2\nctg\t2\t2\nctg\t4\t2\nctg\t5\t3\n";
close $converter_in;
my $converted = do { local $/; <$converter_out> };
my $converter_errors = do { local $/; <$converter_err> };
waitpid($converter_pid, 0);
is($? >> 8, 0, 'samtools-depth converter completes');
is($converter_errors, '', 'samtools-depth converter emits no warnings');
is($converted, "ctg\t0\t2\t2\nctg\t3\t4\t2\nctg\t4\t5\t3\n",
   'bedGraph conversion preserves coordinate gaps and emits the final interval');

my $gene_cat = read_file(File::Spec->catfile($root, 'secScripts', 'geneCat.pl'));
unlike($gene_cat, qr/rm -rf \$GCdir\/\* \$tmpDir\*/, 'geneCat has no wildcard clean-start deletion');
unlike($gene_cat, qr/system "rm -r \$metaGD\/\$path2GPdir/, 'geneCat does not delete predictions while inspecting them');
unlike($gene_cat, qr/system "rm -rf \$metaGD\/\$path2CS/, 'geneCat does not delete contig stats while inspecting them');
like($gene_cat, qr/genemat\.done\.sh/, 'matrix completion uses a convergence job');
like($gene_cat, qr/No usable assembly.*if \$requireAllAssemblies/s,
     'missing assemblies fail only when requireAllAssemblies is enabled');
like($gene_cat, qr/\$map\{\$smpl\}\{assFinSmpl\} eq \$smpl/,
     'assembly-group precheck recognizes the explicitly final assembly sample');
unlike($gene_cat, qr/!\s*fileGZe\("\$metaGD\/scaffolds\.fasta\.filt"\) \|\| !-e "\$metaGD\/longReads/,
       'hybrid assembly precheck does not require both short- and long-read assemblies');
unlike($gene_cat, qr/my \$cmd \.= "\$kaijBin/, 'Kaiju command is initialized before concatenation');

my $parse = read_file(File::Spec->catfile($root, 'secScripts', 'functions', 'parseBlastFunct2.pl'));
like($parse, qr/CNT_\$\{minBLE\}_\$\{minPID\}/, 'functional result checks use threshold and percent identity');
like($parse, qr/\.\$normMethod\.gene\.cnts\.gz/, 'functional result checks include normalization in output names');

my $abr_db = File::Spec->catdir($tmp, 'abr-db');
mkdir $abr_db or die "Cannot create $abr_db: $!";
write_file(File::Spec->catfile($abr_db, 'ardb.tabs.parsed'),
           "SUB\tx\tSYM\tCAT\tx\tx\tx\t80\n");
write_file(File::Spec->catfile($abr_db, 'ardb_and_reforg_mapping'),
           "x\tSYM\tx\tDRUG\n");
write_file(File::Spec->catfile($abr_db, 'ardb_vs_reforg9f.overlap90shortest_famthres_or_symbol.sorted.besthit'),
           "ALT\tSUB\n");
my $abr_blast = File::Spec->catfile($tmp, 'abr.srt.gz');
my $blast_text = join '',
    "paired/1\tSUB\t90\t50\t0\t0\t1\t50\t1\t50\t1e-20\t100\n",
    "paired/2\tSUB\t92\t50\t0\t0\t1\t50\t51\t100\t1e-20\t105\n",
    "second-only/2\tSUB\t95\t80\t0\t0\t1\t80\t1\t80\t1e-20\t120\n";
gzip(\$blast_text => $abr_blast) or die "Cannot create $abr_blast: $GzipError";
my $abr_genes = File::Spec->catfile($tmp, 'abr.genes.txt');
my $abr_cats = File::Spec->catfile($tmp, 'abr.cats.txt');
is(system($^X, '-I' . $root,
          File::Spec->catfile($root, 'secScripts', 'functions', 'ABRblastFilter2.pl'),
          $abr_blast, $abr_genes, $abr_cats, $abr_db), 0,
   'ABR filter handles paired and read-2-only hits');
my $abr_output = read_file($abr_cats);
like($abr_output, qr/^paired\/1\tSUB\t/m, 'ABR filter combines a paired hit without losing its subject key');
like($abr_output, qr/^second-only\/2\tSUB\t/m, 'ABR filter retains read-2-only hits');
ok(-e "$abr_blast.stone", 'ABR completion marker is written after successful output');

my $mgs = read_file(File::Spec->catfile($root, 'secScripts', 'MGS.pl'));
like($mgs, qr/Select exactly one quality checker/, 'MGS rejects ambiguous CheckM/CheckM2 configuration');
like($mgs, qr/runCheckM\(\$binCanDir,\$ChkMevalF/, 'MGS supports CheckM1 for canopy quality checks');
like($mgs, qr/my \$sco = \$spl\[12\]-\(\$spl\[13\]\*2\)/,
     'MAG replacement score uses completeness and contamination columns');
unlike($mgs, qr/my \$testKey = ">\$\{cc\}_\$cnt";\s*my \$curGene = ""; my \$cnt=0/s,
       'contig gene lookup does not use the outer counter before initializing its own counter');
like($mgs, qr/if \(\$LOGstr =~ m\/:::Correct:\/\)/, 'Rhcl success output is recognized without a stray quote');
unlike($mgs, qr/m\/\\":::Correct:/, 'Rhcl parser no longer requires an impossible leading quote');
like($mgs, qr/test -s \$GTDBtaxF.*touch \$GTDBtaxSto/s,
     'GTDB checkpoint follows validation of final taxonomy outputs');
like($mgs, qr/test -s \$annoDir\/specI\.tax\\n";\s*\$cmdSI \.= "touch \$ABmgsSton/s,
     'MGS abundance checkpoint follows final output validation');
like($mgs, qr/_touch_checkpoint\(\$iniMB2sto\) unless -e \$iniMB2sto \|\| \@missedMAGs/,
     'missing MAG groups prevent the global MAG checkpoint from becoming sticky');
unlike($mgs, qr/foreach my \$Doo \(\@DoosD\)\{\s*last if \(-e "\$iniMB2sto"\)/s,
       'MGS validates MAG outputs even when a previous global checkpoint exists');

my $mataf4_stats = read_file(File::Spec->catfile($root, 'MATAF4.pl'));
like($mataf4_stats, qr/sub _smpl_stats_columns.*?sub _metag_stats_text/s,
     'sample statistics use one central ordered schema and final serializer');
like($mataf4_stats, qr/return \{ SNP_TotalResolvedBp=>/,
     'statistics helpers return named values instead of tab-delimited fragments');
like($mataf4_stats, qr/ref\(\$seq_set->\{pa1\}\) eq 'ARRAY'/,
     'sample statistics validate optional read arrays before dereferencing');
like($mataf4_stats, qr/\$map\{\$SmplN\}\{inputFileSizeMB\}/,
     'sample statistics use their sample argument for input size');
like($mataf4_stats, qr/\$values\{RawInputSizeSub\}.*?inputXFileSizeMB/s,
     'sample statistics report supplementary raw input size separately');
unlike($mataf4_stats, qr/my \@sdm = qw\(SDMVersion/,
       'metagStats schema does not expose the internal SDM version');
unlike($mataf4_stats, qr/system "rm -rf \$inD\/assemblies\/metag\/corrected"/,
       'sample statistics do not delete assembly data');
unlike($mataf4_stats, qr/sub smplStats\(\)/,
       'sample statistics no longer declare a misleading zero-argument prototype');
like($mataf4_stats, qr/getContamination\([^;]+prepEBI[^;]+\);/s,
     'EBI contamination fields are emitted unconditionally for a stable schema');
like($mataf4_stats, qr/\$value =~ s\/\[\\t\\r\\n\]\+\/ \/g/,
     'central serialization prevents embedded delimiters from corrupting metagStats');
like($mataf4_stats, qr/my %sampleStats;.*?my \@sampleStatsOrder/s,
     'statistics are retained in a central per-sample object');
like($mataf4_stats, qr/grep \{ \$observed\{\$_\} \} \@preferred/,
     'metagStats emits only columns containing an observed value');
unlike($mataf4_stats, qr/my \$statStr\b/,
       'sample-wise tab-string accumulation has been removed');
like($mataf4_stats, qr/sub getHybridAssemblyStats.*?HybridAssemblyComparison\.tsv/s,
     'hybrid comparative assembly metrics are merged into sample statistics');

my ($central_stats_code) = $mataf4_stats =~ /(sub _smpl_stats_columns.*?)(?=\nsub sdmStats)/s;
ok(defined($central_stats_code), 'central statistics implementation can be isolated for testing');
my $central_eval = "package TestCentralSampleStats; sub getBinSubdirName { return 'B'.\$_[0]; }\n"
    . $central_stats_code;
eval $central_eval;
is($@, '', 'central statistics implementation compiles independently');
my @preferred_columns = TestCentralSampleStats::_smpl_stats_columns();
my %preferred_position;
@preferred_position{@preferred_columns} = (0 .. $#preferred_columns);
my @pipeline_markers = qw(
	RawInputSize FilteredContaRdsPerc_EBI totRds SDMAcceptedPercent
	FilteredContaRdsPerc Merged AvgGenomeSizeEst ContigN50 HybridPreassemblyCount
	ReadsPaired OpticalDuplicates BreakpointCount GeneNumber HQ_bins_B1 SNP_TotalResolvedBp
);
ok(!grep({ !exists $preferred_position{$_} } @pipeline_markers),
	'pipeline-order statistics markers all occur in the preferred schema');
is_deeply(
	[sort { $preferred_position{$a} <=> $preferred_position{$b} } @pipeline_markers],
	\@pipeline_markers,
	'metagStats program blocks follow workflow submission order');
like($mataf4_stats, qr/AssemblyBreakpointPercent.*?GeneCodingPercent.*?SNPsPerMbp.*?INDELsPerMbp/s,
	'useful assembly, gene, and variant-density statistics are reported');
like($mataf4_stats, qr/"\$\{SCdir\}_total_bins".*?\$totBins/s,
	'binner statistics expose total bins without requiring downstream column arithmetic');
my %central_fixture = (
	A => { DIR => '/sample/A', values => { RawInputSize => '1.000G', RawInputSizeSub => '0.250G', BreakpointCount => '' } },
	B => { DIR => '/sample/B', values => { RawInputSize => '2.000G', HybridFinalN50 => 5000 } },
);
my $central_text = TestCentralSampleStats::_metag_stats_text(\%central_fixture, [qw(A B)]);
like($central_text, qr/^SMPLID\tDIR\tRawInputSize\tRawInputSizeSub\tHybridFinalN50$/m,
     'final statistics header contains only populated columns in preferred order');
unlike($central_text, qr/BreakpointCount/,
       'globally empty statistics columns are omitted');
my @central_lines = split /\n/, $central_text;
my @central_widths = map { scalar(split /\t/, $_, -1) } @central_lines;
is_deeply(\@central_widths, [5, 5, 5],
          'central serialization keeps every sample aligned to the selected columns');

my ($hybrid_stats_code) = $mataf4_stats =~ /(sub getHybridAssemblyStats.*?)(?=\nsub smplStats)/s;
ok(defined($hybrid_stats_code), 'hybrid statistics parser can be isolated for testing');
eval "package TestHybridSampleStats; $hybrid_stats_code";
is($@, '', 'hybrid statistics parser compiles independently');
my $hybrid_values = TestHybridSampleStats::getHybridAssemblyStats($comparison);
is($hybrid_values->{HybridPreassemblyCount}, 1,
   'hybrid statistics report its preassembly count');
is($hybrid_values->{HybridPreassemblyN50}, 30_000,
   'hybrid statistics include the source assembly N50');
is($hybrid_values->{HybridFinalN50}, 5_000,
   'hybrid statistics include the final assembly N50');

done_testing();
