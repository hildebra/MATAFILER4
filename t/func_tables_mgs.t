use strict;
use warnings;

# Per-MAG function tables (yoke) for every MGS: input preparation in
# Mods::FuncTables and the secScripts/MGS/funcTablesMGS.pl driver, run against
# the yoke stand-in t/bin/yoke, which fails when yoke's input assumptions are
# broken.

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IPC::Open3 qw(open3);
use Symbol qw(gensym);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use lib File::Spec->catdir($Bin, 'lib');
use MFTestConfig;
use Mods::FuncTables qw(
	parse_function_table_list function_table_scheme_names default_function_table_schemes
	function_table_stage_stones plan_function_tables natural_mgs_sort
	normalise_magvsgc normalise_gene_annotations write_function_descriptions
);

my $root = File::Spec->rel2abs(File::Spec->catdir($Bin, '..'));
chmod 0755, "$root/t/bin/pigz", "$root/t/bin/yoke";
my $tmp = tempdir(CLEANUP => 1);

sub write_file {
	my ($path, $text) = @_;
	make_path((File::Spec->splitpath($path))[1]);
	if ($path =~ /\.gz$/) { gzip(\$text => $path) or die "$path: $GzipError"; return; }
	open my $fh, '>', $path or die "Cannot write $path: $!";
	print {$fh} $text;
	close $fh or die "Cannot close $path: $!";
}

sub read_file {
	my ($path) = @_;
	my $text;
	if ($path =~ /\.gz$/) { gunzip($path => \$text) or die "$path: $GunzipError"; return $text; }
	open my $fh, '<', $path or die "Cannot read $path: $!";
	local $/;
	return <$fh>;
}

sub run_script {
	my (@args) = @_;
	my $err = gensym;
	my $pid = open3(undef, my $out, $err, $^X, "-I$root", "-I$root/t/lib", '-MMFTestConfig', @args);
	my $stdout = do { local $/; <$out> // '' };
	my $stderr = do { local $/; <$err> // '' };
	waitpid($pid, 0);
	return ($? >> 8, "$stdout$stderr");
}

# --- scheme lists ---------------------------------------------------------
is_deeply([parse_function_table_list('0', 'scheme', function_table_scheme_names())], [],
	'"0" selects no scheme');
is_deeply([parse_function_table_list('none', 'scheme', function_table_scheme_names())], [],
	'"none" selects no scheme');
is_deeply([parse_function_table_list('KEGG, emKO,KEGG', 'scheme', function_table_scheme_names())],
	['KEGG', 'emKO'], 'scheme lists keep order and drop duplicates');
is_deeply([parse_function_table_list(default_function_table_schemes(), 'scheme', function_table_scheme_names())],
	[function_table_scheme_names()], 'the default selects every scheme');
ok(!eval { parse_function_table_list('KEGG,KGG', 'scheme', function_table_scheme_names()); 1 },
	'a misspelt scheme is an error, not a silent skip');
like($@, qr/Unknown scheme entry 'KGG'/, 'the error names the unknown scheme');
is_deeply([natural_mgs_sort(qw(MGS.10 MGS.2 MGS.1 other))], [qw(MGS.1 MGS.2 MGS.10 other)],
	'MGS IDs sort by their number');

# --- MAGvsGC normalisation ------------------------------------------------
# Layout of the Perl clusterMAGs path: lower-case other_genes, no taxonomy,
# rows in sample order (MGS.10 interleaved), an unassigned MAG and MAGs that
# yoke's quality filter would drop.
my $perlMAG = "$tmp/perl/MAGvsGC.txt.gz";
write_file($perlMAG, join('',
	"MAG\tMGS\tRepresentative4MGS\tCompleteness\tContamination\tLCAcompleteness\tCOG1\tother_genes\n",
	"S1__1\tMGS.10\t*\t90\t1\t1\t1\t1,2,,3\n",
	"S1__2\tMGS.2\t*\t80\t2\t1\t\t4,5\n",
	"S1__3\t?\t\t70\t1\t1\t\t6\n",
	"S2__1\tMGS.10\t\t85\t1\t1\t1\t1,3\n",
	"S2__2\tMGS.2\t\t0\t1\t1\t\t7\n",
	"S2__3\tMGS.2\t\t75\t100\t1\t\t8\n",
	"Cano__4\tMGS.4\t*\t0\t0\t1\t\t9\n",
));
my $perlNorm = "$tmp/perl/norm.txt.gz";
my $stats = normalise_magvsgc($perlMAG, $perlNorm, collect_genes => 1);
my @rows = split /\n/, read_file($perlNorm);
my $header = shift @rows;
like($header, qr/\tother_Genes\z/, 'the gene list column is renamed to other_Genes');
is_deeply([map { (split /\t/)[0] } @rows], [qw(S1__2 S1__1 S2__1)],
	'rows are grouped by MGS in MGS order, input order within an MGS');
is($stats->{mgs}, 2, 'MGS left after filtering');
is($stats->{mags}, 3, 'MAGs kept');
is($stats->{unassigned}, 1, 'MAGs without MGS are counted and left out');
is($stats->{filtered}, 3, 'MAGs yoke would filter out are left out, so no MGS is left empty');
is($stats->{taxonomy}, 0, 'the Perl layout has no taxonomy');
is_deeply([sort { $a <=> $b } keys %{$stats->{genes}}], [1 .. 5], 'genes of kept MAGs, without contig separators');

# Layout of the clusterMAGs binary: 22 fixed columns including taxonomy.
my @nativeCols = qw(MAG MGS Representative4MGS Match2MGS Uniqueness AssociatedMGS Completeness
	Contamination LCAcompleteness N50 N_Genes GC CodingDensity CentreScore CompoundScore
	Domain Phylum Class Order Family Genus Species);
sub native_row {
	my ($mag, $mgs, $compl, $genes, $species) = @_;
	my %v = (MAG => $mag, MGS => $mgs, Representative4MGS => '', Completeness => $compl,
		Contamination => 1, Domain => 'd__Bacteria', Phylum => 'p__Firm', Class => 'c__C',
		Order => 'o__O', Family => 'f__F', Genus => 'g__G', Species => $species);
	return join("\t", (map { $v{$_} // 0 } @nativeCols), 'g1', $genes) . "\n";
}
my $nativeHeader = join("\t", @nativeCols, 'COG0001', 'other_Genes') . "\n";
my $nativeMAG = "$tmp/native/MAGvsGC.txt.gz";
write_file($nativeMAG, $nativeHeader
	. native_row('S1__1', 'MGS.1', 95, '1,2,,3', 's__A')
	. native_row('S1__2', 'MGS.2', 90, '4,5', 's__B')
	. native_row('S2__1', 'MGS.1', 92, '1,3,6', 's__A'));
my $nativeStats = normalise_magvsgc($nativeMAG, "$tmp/native/norm.txt.gz");
is($nativeStats->{taxonomy}, 1, 'the binary layout carries taxonomy');
ok(!exists $nativeStats->{genes}, 'genes are only collected on request');

write_file("$tmp/bad/MAGvsGC.txt.gz", "MAG\tCompleteness\n");
ok(!eval { normalise_magvsgc("$tmp/bad/MAGvsGC.txt.gz", "$tmp/bad/norm.txt.gz"); 1 },
	'a report without the needed columns is rejected');
like($@, qr/lacks column\(s\): MGS, Representative4MGS, Contamination, other_Genes/,
	'the error names every missing column');

# --- gene annotations -----------------------------------------------------
my $anno = "$tmp/anno/DIAass_KGM.srt.gzgeneAss.gz";
write_file($anno, join('',
	"# comment\n",
	"1\tK00001,K00002\n",
	"2\tK00003|K00004\n",
	"3\tK00005;M00001,M00002\n",
	"4\t-\n",
	"5\tK00006\tcategory\tdefinition\n",
	"9\tK00007\n",
	"6\tK00001,K00001\n",
));
my $annoStats = normalise_gene_annotations($anno, "$tmp/anno/norm.gz", { map { $_ => undef } 1 .. 6 });
is(read_file("$tmp/anno/norm.gz"),
	"1\tK00001\n1\tK00002\n2\tK00003\n2\tK00004\n3\tK00005\n5\tK00006\n6\tK00001\n",
	'one row per function: "," and "|" split, first ";" level kept, extra columns, "-" and other genes dropped');
is($annoStats->{genes}, 5, 'annotated genes counted');
is_deeply([sort keys %{$annoStats->{functions}}], [map { sprintf 'K%05d', $_ } 1 .. 6], 'function IDs collected');

# --- descriptions ---------------------------------------------------------
write_file("$tmp/descr/KEGG.descr", "K00001\tadh\talcohol \"dehydrogenase\"\nK00002\t\n#K00003\tcomment\n");
my $described = write_function_descriptions(["$tmp/descr/KEGG.descr", "$tmp/descr/missing"],
	{ K00001 => undef, K00002 => undef, K00003 => undef }, "$tmp/descr/out.descr");
is($described, 1, 'only functions with a non-empty description count as described');
is(read_file("$tmp/descr/out.descr"),
	"K00001\talcohol 'dehydrogenase'\nK00002\tK00002\nK00003\tK00003\n",
	'two columns; last field as description, the ID otherwise; quotes made safe for pandas');

# --- planning against geneCat checkpoints ---------------------------------
my $GCd = "$tmp/GC";
my $stones = function_table_stage_stones($GCd);
write_file("$GCd/Anno/Func/DIAass_KGM.srt.gzgeneAss.gz", "1\tK00001,K00002\n2\tK00003\n3\tK00001\n4\tK00009\n6\tK00003\n");
write_file("$GCd/Anno/Func/KEGG.descr", "K00001\tfirst\nK00002\tsecond\n");
write_file("$GCd/Anno/Func/emapper/eggNOGmapper_NOG.geneAss.gz", "1\tCOG0001\n4\tCOG0002\n");
my ($ready, $skipped) = plan_function_tables($GCd, ['KEGG', 'CZy', 'emNOG'], $stones);
is_deeply($ready, [], 'nothing runs before the annotation stages finished');
write_file($stones->{func}, "{}\n");
($ready, $skipped) = plan_function_tables($GCd, ['KEGG', 'CZy', 'emNOG'], $stones);
is_deeply($ready, ['KEGG'], 'diamond schemes run once FuncAssign finished');
is_deeply([map { $_->[0] } @$skipped], ['CZy', 'emNOG'], 'the rest are skipped');
like($skipped->[0][1], qr/no annotation file/, 'a database that was not annotated is skipped with its reason');
like($skipped->[1][1], qr/stage 'emap' has not finished/, 'eggNOG-mapper schemes wait for their stage');
write_file($stones->{emap}, "{}\n");
($ready) = plan_function_tables($GCd, ['KEGG', 'emNOG'], $stones);
is_deeply($ready, ['KEGG', 'emNOG'], 'both stages finished: both schemes run');

# --- driver, end to end with the yoke stand-in ----------------------------
my $binD = "$tmp/GC/Bin_SB";
write_file("$binD/MAGvsGC.txt.gz", $nativeHeader
	. native_row('S1__1', 'MGS.1', 95, '1,2,,3', 's__A')
	. native_row('S1__2', 'MGS.2', 90, '4,6', 's__B')
	. native_row('S2__1', 'MGS.1', 92, '1,3', 's__A')
	. native_row('S3__1', '?', 92, '1', ''));
my $driver = "$root/secScripts/MGS/funcTablesMGS.pl";
my $outD = "$binD/Annotation/FunctionTables";
local $ENV{MF4_TEST_YOKE_LOG} = "$tmp/yoke.log";
my ($exit, $output) = run_script($driver, '-GCd', $GCd, '-binD', $binD, '-schemes', 'KEGG,emNOG');
is($exit, 0, 'the driver finishes') or diag($output);
ok(-s "$outD/function_tables.manifest.tsv", 'the manifest is written');
# Gene 1 has K00001 and K00002; yoke keeps one function per gene (the last row).
is(read_file("$outD/KEGG/MGS.1.tsv.gz"), "function\tS1__1\tS2__1\nK00001\t1\t1\nK00002\t1\t1\nK00003\t1\t0\n",
	'MGS.1 counts cover both MAGs of the interleaved MGS');
is(read_file("$outD/KEGG/MGS.2.tsv.gz"), "function\tS1__2\nK00003\t1\nK00009\t1\n", 'MGS.2 counts');
ok(-s "$outD/KEGG.tsv" && -s "$outD/emNOG.tsv", 'per-MGS mean tables for every scheme');
for my $set (qw(BSB GBM)) {
	ok(-s "$outD/modules/KEGG/MGS.1.$set.MODscore" && -s "$outD/modules/KEGG/MGS.2.$set.MODscore",
		"$set module scores for every KEGG MGS table");
}
ok(!-e "$outD/modules/emNOG", 'no module scores for a non-KEGG scheme');
ok(-s "$outD/summarised/KEGG/Species.tsv" && -s "$outD/summarised/emNOG/Domain.tsv",
	'rank summaries are made when MAGvsGC has taxonomy');
ok(-x "$outD/input/bin/rtk", 'yoke finds rtk through the launcher on PATH');
my @calls = split /\n/, read_file("$tmp/yoke.log");
is(scalar(grep { /^function_counts\t/ } @calls), 2, 'one function_counts call per scheme');
ok(!grep({ /^module_score\t.*presence\./ } @calls), 'presence tables are not module-scored');

# A rerun after re-clustering removes tables of MGS that no longer exist.
write_file("$binD/MAGvsGC.txt.gz", $nativeHeader . native_row('S1__1', 'MGS.1', 95, '1,2', 's__A'));
($exit, $output) = run_script($driver, '-GCd', $GCd, '-binD', $binD, '-schemes', 'KEGG', '-modules', '0');
is($exit, 0, 'the driver reruns') or diag($output);
ok(!-e "$outD/KEGG/MGS.2.tsv.gz", 'tables of a removed MGS are deleted');
ok(!-e "$outD/modules/KEGG", '-modules 0 scores no modules');

# Failures stop the job before the manifest is written.
unlink "$outD/function_tables.manifest.tsv";
{
	local $ENV{MF4_TEST_YOKE_FAIL} = 'module_score';
	($exit, $output) = run_script($driver, '-GCd', $GCd, '-binD', $binD, '-schemes', 'KEGG');
}
isnt($exit, 0, 'a failed yoke call fails the driver');
like($output, qr/KEGG BSB module scores failed/, 'the failing step is named');
ok(!-e "$outD/function_tables.manifest.tsv", 'no manifest after a failure');
($exit, $output) = run_script($driver, '-GCd', $GCd, '-binD', $binD, '-schemes', 'CZy');
isnt($exit, 0, 'a scheme without annotation fails the driver');
like($output, qr{Annotation for scheme CZy not found; looked for: .*DIAass_CZy}, 'the missing file is named');
($exit, $output) = run_script($driver, '-GCd', $GCd, '-binD', $binD, '-schemes', '0');
isnt($exit, 0, 'the driver needs at least one scheme');

done_testing();
