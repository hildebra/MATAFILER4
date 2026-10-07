use strict;
use warnings;
no warnings qw(once redefine prototype); #subroutine stubs for the extracted pipeline code
use Test::More;
use File::Temp qw(tempdir);
use File::Path qw(make_path remove_tree);
use File::Spec;
use FindBin qw($Bin);
use IO::Compress::Gzip qw(gzip $GzipError);
use lib "$Bin/..", "$Bin/lib";
use MFTestConfig;
use Mods::Checkpoint ();
use Mods::FuncTools ();
use Mods::IO_Tamoc_progs ();

#regression tests for the functional-annotation audit of 1 October 2026 (docs/audits/2026-10-01/report.md)
our $root = File::Spec->rel2abs("$Bin/..");
my $tmp = tempdir(CLEANUP => 1);
$ENV{PERL5OPT} = "-I$root -I$root/t/lib -MMFTestConfig";

sub write_file { my ($p, $c) = @_; open my $fh, '>', $p or die "$p: $!"; print {$fh} $c; close $fh or die $!; }
sub read_file { my ($p) = @_; open my $fh, '<', $p or die "$p: $!"; local $/; my $c = <$fh>; close $fh; return $c; }
sub touchf { for my $f (@_) { make_path($1) if $f =~ m{^(.*)/[^/]+$}; write_file($f, "x\n"); } }
sub source_sub {
	my ($src, $name) = @_;
	my ($s) = $src =~ /^(sub \Q$name\E\b[^\n]*\{.*?^\})/ms;
	die "cannot find sub $name\n" unless defined $s;
	return $s;
}
my $gcSrc = read_file("$root/secScripts/geneCat.pl");
my $mfSrc = read_file("$root/MATAF4.pl");

# ---------------- HMMER best hit (geneCat -mode FOAM / ABR) ----------------
subtest 'HMMER best hit helper' => sub {
	sub dom { my ($g, $h, $s) = @_;
		return join(' ', $g, '-', 100, $h, '-', 200, '1e-10', $s, 0.1, 1, 1, '1e-10', '1e-10', $s, 0.1, 1, 100, 1, 100, 1, 100, 0.9, 'some description') . "\n"; }
	my $in = "$tmp/dom.txt";
	write_file($in, "# target name accession tlen query name\n" . join('',
		dom('g2', 'K00010', 300), dom('g1', 'K1', 300), dom('g1', 'K2', 30), dom('g2', 'K00020', 30),
		dom('g3', 'K00100', 200), dom('g3', 'K00200', 26), dom('g3', 'K00300', 90), dom('g4', 'K5', 10)));
	my $out = `$^X $root/secScripts/functions/hmmBestHit.pl $in`;
	is($?, 0, 'helper runs');
	my %best = map { my @f = split ' '; ($f[0] => $f[3]) } split /\n/, $out;
	is_deeply(\%best, {g1 => 'K1', g2 => 'K00010', g3 => 'K00100'},
		'highest-scoring hit per gene, also for genes after the first; hits below 25 ignored');
	write_file("$tmp/dom.empty", "# no hits\n");
	$out = `$^X $root/secScripts/functions/hmmBestHit.pl $tmp/dom.empty`;
	ok($? == 0 && $out eq '', 'a chunk without hits gives empty output, not a crash');
	like(Mods::IO_Tamoc_progs::getProgPaths('hmmBestHit_scr'), qr/^perl \S*hmmBestHit\.pl$/, 'configured helper is the Perl port');
	my $foam = source_sub($gcSrc, 'FOAMassign');
	like($foam, qr/\$\{DB\}split_/, 'FOAM and ABR use separate catalog splits');
	unlike($foam, qr/\$rmBin -f -r \$tmpOut\* \$subFls/, 'chunk jobs no longer delete their (shared) input chunk');
};

# ---------------- VFDB table ----------------
subtest 'VFDB table build and staleness' => sub {
	my $vd = "$tmp/vfdb"; make_path($vd);
	my $a = ">VFG000001(gb|WP_1) (fdxA) ferredoxin [2Fe-2S] family protein [Iron uptake (VF0123) - Nutritional/Metabolic factor (VFC0272)] [Escherichia coli]\nMKK\n"
		. ">VFG000002(gb|WP_2) (iucA) aerobactin synthase [Iron uptake (VF0123) - Nutritional/Metabolic factor (VFC0272)] [Escherichia coli]\nMKK\n";
	write_file("$vd/VFDB_setA_pro.fas", $a);
	write_file("$vd/VFDB_setB_pro.fas", $a);
	ok(Mods::FuncTools::vfTabStale($vd), 'missing VF.tab must be built');
	is(system("$^X $root/secScripts/functions/prepVFDB.pl $vd >/dev/null"), 0, 'prepVFDB.pl runs');
	my %L1 = map { my @c = split /\t/; ($c[0] => $c[3]) } split /\n/, read_file("$vd/VF.tab");
	is($L1{'VFG000001(gb|WP_1)'}, $L1{'VFG000002(gb|WP_2)'}, 'a bracket in the description does not change the VF name');
	unlike($L1{'VFG000001(gb|WP_1)'}, qr/2Fe/, 'VF name taken from the VF bracket');
	my $old = time - 1000;
	utime($old, $old, "$vd/VFDB_setA_pro.fas", "$vd/VFDB_setB_pro.fas");
	ok(!Mods::FuncTools::vfTabStale($vd), 'current VF.tab is reused');
	utime(time + 10, time + 10, "$vd/VFDB_setB_pro.fas");
	ok(Mods::FuncTools::vfTabStale($vd), 'VF.tab older than a VFDB FASTA is rebuilt');
	utime($old, $old, "$vd/VFDB_setB_pro.fas");
	write_file("$vd/VF.tab", "VFG000001(gb|WP_1)\tfdxA\tdesc\tVF0123_x\n");
	ok(Mods::FuncTools::vfTabStale($vd), 'VF.tab in the old 4/5-column layout is rebuilt');
};

# ---------------- .length table ----------------
subtest 'length table is written atomically' => sub {
	my $fa = "$tmp/len.faa";
	write_file($fa, ">a\nMKK\n>b desc\nMK\nKKK\n");
	is(system("$^X $root/secScripts/assemblies/geneLengthFasta.pl $fa $tmp/len.length >/dev/null"), 0, 'length script runs');
	is(read_file("$tmp/len.length"), "a\t3\nb\t5\n", 'lengths');
	is(scalar(() = glob("$tmp/len.length.tmp*")), 0, 'no temporary file left');
	write_file("$tmp/bad.faa", "MKK\n");
	isnt(system("$^X $root/secScripts/assemblies/geneLengthFasta.pl $tmp/bad.faa $tmp/bad.length 2>/dev/null"), 0, 'invalid FASTA fails');
	ok(!-e "$tmp/bad.length" && !glob("$tmp/bad.length.tmp*"), 'a failed run leaves no table behind');
};

# ---------------- parseBlastFunct2.pl ----------------
subtest 'hit parser: description separators, CAZy families, truncated input' => sub {
	my $pbSrc = read_file("$root/secScripts/functions/parseBlastFunct2.pl");
	my $hs = source_sub($pbSrc, 'hieraSafe');
	my $hieraSafe = eval "sub { $hs; return hieraSafe(\@_); }" or die $@;
	is($hieraSafe->('Porters (uniporters, symporters; antiporters)|x'), 'Porters (uniporters_ symporters_ antiporters)_x',
		'rtk separators in descriptions are replaced');
	like($pbSrc, qr/join\(",",\@curCOGs\)\."\\t"\.join\(",",uniq\(\@subs2\)\)/, 'CAZy families joined with the AND separator');
	unlike($pbSrc, qr/join\(";",\@curCOGs\)/, 'no CAZy family joined with the hierarchy separator');
	#truncated gzip (deflate data complete, trailer missing): must fail, no stone, no per-gene file
	my $d = "$tmp/trunc"; make_path("$d/tmp");
	write_file("$d/db.length", "x:A\t1000\n");
	my $txt = join('', map { "q$_\tx:A\t60\t300\t0\t0\t1\t300\t1\t300\t1e-50\t300\n" } 1 .. 200);
	my $gz; gzip(\$txt => \$gz) or die $GzipError;
	my $in = "$d/DIAass_ACL.srt.gz";
	write_file($in, substr($gz, 0, length($gz) - 8));
	my $rc = system("$^X $root/secScripts/functions/parseBlastFunct2.pl -i $in -DB ACL -mode 2 -singleSpecies 1 -calcGeneLengthNorm 0 -percID 25 -minAlignLen 30 -minBitScore 45 -eval 1e-8 -LF $d/db.length -tmp $d/tmp/ -queryType genes >$d/log 2>&1");
	isnt($rc, 0, 'truncated input fails the parse');
	ok(!-e "$in.stone" && !-e "${in}geneAss.gz", 'no stone and no per-gene file for truncated input');
};

# ---------------- combine_DIA.pl ----------------
subtest 'cohort merge of read-based tables' => sub {
	my $P = "$tmp/comb";
	my $DBD = "$P/DiamDB"; make_path($DBD, "$P/in/S0", "$P/in/S1", "$P/tmp");
	write_file("$DBD/tcdb.faa.length", "gnl|TC-DB|P1|1.A.1.1.1\t100\n");
	write_file("$DBD/TCDBhir.txt", "");
	write_file("$P/map.txt", "#SmplID\tPath\n#OutPath\t$P/out/\n#RunID\tR\n#DirPath\t$P/in/\nS0\tS0\nS1\tS1\n");
	my $dia = "$P/out/R/S1/diamond"; make_path($dia);
	my $hits = "r1/1\tgnl|TC-DB|P1|1.A.1.1.1\t90\t60\t0\t0\t1\t180\t1\t60\t1e-30\t200\n"
		. "r2/1\tgnl|TC-DB|P1|1.A.1.1.1\t90\t60\t0\t0\t1\t180\t1\t60\t1e-30\t200\n";
	gzip(\$hits => "$dia/dia.TCDB.blast.srt.gz") or die $GzipError;
	is(system("$^X $root/secScripts/functions/parseBlastFunct2.pl -i $dia/dia.TCDB.blast.srt.gz -DB TCDB -eval 1e-7 -percID 40 -minAlignLen 20 -minPercSbjCov 0.1 -mode 0 -queryType reads -LF $DBD/tcdb.faa.length -reportDomains 0 -DButil $DBD/ -tmp $P/tmp >$P/parse.log 2>&1"),
		0, 'TCDB read-mode parse');
	is(system("cd $P && $^X $root/secScripts/functions/combine_DIA.pl $P/out TCDB $P/map.txt >$P/comb.log 2>&1"), 0,
		'merge works when the first map sample has no diamond/ directory');
	my ($cat) = glob("$P/out/pseudoGC/FUNCT/TCDB/TCDB.CAT.mat.ALL.cnt.*.txt");
	ok(defined($cat) && read_file($cat) =~ /\n1\.A\.1\t/, 'TCDB category matrix has the class rows (.CATcnts is read)');
};

# ---------------- MATAF4.pl read-based search caller ----------------
our (%MFopt, %progStats, @jobs, %LIBS, $curSmpl, $pigzBin, $logDir, $avx2Constr, $JNUM, %HDDspace, $QSBoptHR);
$curSmpl = 'S'; $pigzBin = 'pigz'; $logDir = "$tmp/log/"; $avx2Constr = ''; $JNUM = 1; %HDDspace = (diamond => 1);
$QSBoptHR = {constraint => [], tmpSpace => 0, General_Hosts => []};
{
	no warnings 'redefine';
	*main::sampleReadSet = sub { return {merged_library => {}} };
	*main::readLibrariesByScope = sub { return [] };
	*main::getRdLibraries = sub { return %LIBS };
	*main::prepDiamondDB = sub { my ($db) = @_; return ("$db.ref", $db, '') };
	*main::getSpecificDBpaths = sub { return ("$tmp/abrdb/", 'ardb.fa', $_[0]) };
	*main::getProgPaths = sub { return "$^X $root/secScripts/functions/parseBlastFunct2.pl" if $_[0] eq 'secCogBin_scr'; return "prog_$_[0]" };
	*main::qsubSystem = sub { push @jobs, [@_]; return ($_[4], $_[1]) };
}
for my $name (qw(runDiamond prepareDiamondRerun IsDiaRunFinished)) {
	my $s = source_sub($mfSrc, $name);
	eval $s; die $@ if $@;
}
subtest 'read-based search: finished databases, temp dirs, ABR tables' => sub {
	%MFopt = (diaCores => 1, diaRunSensitive => 0, diaFrameshift => 0, diaEVal => '1e-7', DiaPercID => 40, DiaMinAlignLen => 20,
		DiaMinFracQueryCov => 0.1, DiaRmRawHits => 1, globalDiamondDependence => {NOG => 'NOG-1', KGM => 'KGM-1', ABR => 'ABR-1'},
		diamondMem => 7, PABtaxChk => 0);
	my $o = "$tmp/rd/diamond/"; make_path($o);
	touchf("$o/dia.NOG.blast.srt.gz.stone", "$o/dia.NOG.blast.srt.gz.read-counts-v1.stone");
	%LIBS = (1 => ['r2.fq.gz'], 2 => ['r1.fq.gz']); @jobs = (); %progStats = ();
	main::runDiamond($o, "$tmp/db/", "$tmp/scr", 'dep', 'NOG,KGM');
	is_deeply([map { $_->[4] } @jobs], ['_DKGM1', '_DPKGM1'], '-rmRawDiamondHits: a parsed database is not realigned');
	my ($search) = grep { $_->[4] eq '_DKGM1' } @jobs;
	like($search->[1], qr{-t \Q$tmp\E/scr/KGM/}, 'per-database temp dir');
	like($search->[1], qr/LC_ALL=C sort/, 'hits are sorted in byte order');
	@jobs = ();
	my $ab = "$tmp/rdabr/diamond/"; make_path($ab);
	main::runDiamond($ab, "$tmp/db/", "$tmp/scr", 'dep', 'ABR');
	my ($abrParse) = grep { $_->[4] eq '_DPABR1' } @jobs;
	like($abrParse->[1], qr{ABR\.cats\.txt \Q$tmp\E/abrdb/\n}, 'ABR filter reads the ARDB tables from the database dir');
};
subtest 'read-based reparse keeps other databases' => sub {
	my $d = "$tmp/rr";
	make_path("$d/diamond/CNT_1e-7_40");
	for my $db (qw(NOG KGM CZy ABRc TCDB PTV MOH)) {
		touchf("$d/diamond/dia.$db.blast.srt.gz", "$d/diamond/dia.$db.blast.srt.gz.stone", "$d/diamond/CNT_1e-7_40/${db}parse.$db.ALL.cnt.cat.cnts.gz");
	}
	%MFopt = (DoDiamond => 1, maxReqDiaDB => 6, rewriteDiamond => 0, redoDiamondParse => 1, diaEVal => '1e-7', DiaPercID => 40,
		reqDiaDB => 'NOG,KGM,CZy,ABRc,TCDB,PTV');
	main::prepareDiamondRerun("$d/");
	ok(-e "$d/diamond/dia.MOH.blast.srt.gz.stone" && -e "$d/diamond/CNT_1e-7_40/MOHparse.MOH.ALL.cnt.cat.cnts.gz",
		'six requested databases: an unrequested database keeps tables and marker');
	ok(!-e "$d/diamond/dia.NOG.blast.srt.gz.stone" && !glob("$d/diamond/CNT_1e-7_40/NOGparse*"),
		'requested database: marker and tables of CNT_<eval>_<DiaPercID> removed');
	ok(-e "$d/diamond/dia.NOG.blast.srt.gz", 'raw hits kept for reparsing');
};
subtest 'old merged-read searches request the search stage' => sub {
	my $d = "$tmp/prov";
	touchf("$d/diamond/dia.NOG.blast.srt.gz");
	%MFopt = (DoDiamond => 1, reqDiaDB => 'NOG', rewriteDiamond => 0, rewriteAllIfAnyDiamond => 0, doReadMerge => 1);
	is_deeply([main::IsDiaRunFinished($d)], [1, 1], 'reparse of hits without read-count provenance schedules the search stage');
	touchf("$d/diamond/dia.NOG.blast.srt.gz.read-counts-v1.stone");
	is_deeply([main::IsDiaRunFinished($d)], [0, 1], 'hits with provenance: parse only');
	unlink "$d/diamond/dia.NOG.blast.srt.gz.read-counts-v1.stone";
	$MFopt{doReadMerge} = 0;
	is_deeply([main::IsDiaRunFinished($d)], [0, 1], 'without read merging: parse only');
};

# ---------------- geneCat.pl FuncAssign: per-database parameters ----------------
my $pkg = <<'PKG';
package GCProbe;
use strict; use warnings;
use File::Path qw(make_path remove_tree);
use Digest::MD5 ();
use Cwd ();
use Mods::FuncTools qw(assignFuncPerGene calc_modules bigFuncDB);
our ($cdhID, $funcAligner, $qsubDir, $minEVal, $minPerID, $minPercSbjCov, $minPercQueryCov,
	$fastaSplits, $fastaSplitsBig, $redoFunc, $minAlLeng, $minBitSc, %funcDBcutoffs, $rmBin, $pigzBin, $sedBin,
	$rareBin, $countMatrixF, $rtkFunDelims, $touchBin, $GLBtmp, $tmpDir, $GCdir, $curDB_o, @Q);
sub qsubSystem { push @Q, [@_]; return ("J" . scalar(@Q), ""); }
PKG
$pkg .= source_sub($gcSrc, $_) . "\n" for qw(geneCatFunc _gcTmpTag _funcSplitSize _funcSplitDir _funcTmpDir _funcStoneParams _qsbCopy _readFuncParams _writeFuncParams);
eval "$pkg\n1;" or die $@;
my @FQ;
my $dbdir = "$tmp/kegg/"; make_path($dbdir, "$tmp/moddb");
{
	no warnings 'redefine';
	*Mods::FuncTools::qsubSystem = sub { push @FQ, [@_]; return ("F" . scalar(@FQ), ""); };
	*Mods::FuncTools::getProgPaths = sub { my $n = shift; return "$tmp/node" if $n eq 'nodeTmpDir';
		return "$tmp/moddb/" if $n eq 'Module_path_DB'; return "BIN_$n"; };
	*Mods::FuncTools::getSpecificDBpaths = sub { return ($dbdir, 'euk_pro.pep', $_[0]); };
}
touchf(map { "$dbdir$_" } qw(euk_pro.pep euk_pro.pep.db.dmnd euk_pro.pep.length));
{
	no strict 'refs';
	${"GCProbe::$_->[0]"} = $_->[1] for (
		[cdhID => 95], [funcAligner => 'diamond'], [minEVal => 1e-8], [minPerID => 25], [minPercSbjCov => 0.5],
		[minPercQueryCov => 0.8], [fastaSplits => 2], [fastaSplitsBig => 4], [redoFunc => 0], [minAlLeng => 30], [minBitSc => 45],
		[rmBin => 'rm'], [pigzBin => 'pigz'], [sedBin => 'sed'], [rareBin => 'rtk'], [countMatrixF => 'Matrix.mat'],
		[rtkFunDelims => ''], [touchBin => 'touch'], [curDB_o => 'KGM']);
}
my $gc = "$tmp/GC/"; make_path("$gc/Anno/Func", "$gc/LOGandSUB");
$GCProbe::GCdir = $gc; $GCProbe::qsubDir = "$gc/LOGandSUB/";
$GCProbe::GLBtmp = $GCProbe::tmpDir = "$tmp/glob/GC/ID/";
write_file("$gc/compl.incompl.95.prot.faa", join('', map { ">$_\nMAAAAAAAAA\n" } 1 .. 400));
my $outD = "$gc/Anno/Func/";
my $qs = {constraint => [], tmpSpace => 15, qsubDir => "$gc/LOGandSUB/"};
sub runGCF { @GCProbe::Q = (); @FQ = (); return GCProbe::geneCatFunc($gc, GCProbe::_funcTmpDir(), 'KGM', 4, $qs); }
sub chunkJobs { return grep { $_->[0] =~ m{/DKGM\.\d+\.sh$} } @FQ; }
my @products = ("$outD/DIAass_KGM.srt.gz", "$outD/DIAass_KGM.srt.gzgeneAss.gz", "$outD/.KGM.matrix.done");

sub record { return GCProbe::_readFuncParams("$outD/.KGM.params") || {}; }
sub colJob { my ($c) = grep { $_->[0] =~ /colDIAKGM\.sh$/ } @FQ; return $c; }
subtest 'FuncAssign: changed cutoffs or aligner recompute only what they affect' => sub {
	touchf(@products);
	runGCF();
	is(scalar(@GCProbe::Q) + scalar(@FQ), 0, 'finished database without a parameter record (older run): kept, nothing submitted');
	is_deeply([@{record()}{qw(aligner alnEval eval percID minBitScore minAlignLen minPercSbjCov minPercQueryCov)}],
		['diamond', 1e-8, 1e-8, 25, 45, 30, 0.5, 0.8], 'record written for the current settings');
	$GCProbe::minPerID = 40;
	runGCF();
	ok(!-e "$outD/DIAass_KGM.srt.gzgeneAss.gz" && -e "$outD/DIAass_KGM.srt.gz", 'changed identity cutoff: assignments removed, alignments kept');
	is(scalar(chunkJobs()), 0, 'no realignment for a parse cutoff');
	ok(defined(colJob()), 're-interpretation job submitted');
	is(scalar(grep { $_->[0] =~ /KGM_matrix\.sh$/ } @GCProbe::Q), 1, 'matrix rebuilt');
	is(record()->{percID}, 40, 'record follows the new cutoffs');
	touchf(@products);
	$GCProbe::minEVal = 1e-10;
	runGCF();
	is(scalar(chunkJobs()), 0, 'stricter e-value: no realignment');
	ok(-e "$outD/DIAass_KGM.srt.gz" && !-e "$outD/DIAass_KGM.srt.gzgeneAss.gz", 'stricter e-value: alignments kept, assignments redone');
	like(colJob()->[1], qr/-eval 1e-10\n/, 'parser filters with the stricter e-value');
	is_deeply([@{record()}{qw(alnEval eval)}], [1e-8, 1e-10], 'record keeps the alignment e-value');
	#alignment still incomplete (no DIAass): remaining chunks use the alignment e-value of the existing ones
	unlink glob("$outD/DIAass_KGM*"), "$outD/.KGM.matrix.done";
	runGCF();
	my ($chunk) = chunkJobs();
	like($chunk->[1], qr/ -e 1e-08 /, 'chunks still to align use the recorded (looser) alignment e-value');
	like(colJob()->[1], qr/-eval 1e-10\n/, 'and the parser the current one');
	like($chunk->[1], qr/-o \S+\.tmp\.\$\$\.gz /, 'chunk output written under a job-unique temporary name');
	touchf(@products);
	$GCProbe::minEVal = 1e-6;
	runGCF();
	is(scalar(chunkJobs()), 4, 'less strict e-value: realigned');
	ok(!-e "$outD/DIAass_KGM.srt.gz", 'old alignments removed');
	like((chunkJobs())[0][1], qr/ -e 1e-06 /, 'with the new e-value');
	is_deeply([@{record()}{qw(alnEval eval)}], [1e-6, 1e-6], 'record follows the new alignment e-value');
	touchf(@products);
	$GCProbe::funcAligner = 'foldseek';
	touchf("$dbdir/euk_pro.pep.DB3di");
	runGCF();
	is(scalar(chunkJobs()), 4, 'changed aligner: realigned');
	$GCProbe::funcAligner = 'diamond';
	runGCF(); #back to diamond: realigned again
	touchf(@products);
	runGCF();
	is(scalar(@GCProbe::Q) + scalar(@FQ), 0, 'unchanged settings: nothing submitted');
	unlink "$outD/DIAass_KGM.srt.gzgeneAss.gz";
	runGCF();
	ok(defined(colJob()) && !chunkJobs(), 'marker without per-gene assignments: re-interpreted instead of dead-ending the stone job');
};
subtest 'FuncAssign: chunk outputs from an earlier split, collection dependencies' => sub {
	unlink glob("$outD/DIAass_KGM*"), "$outD/.KGM.matrix.done";
	remove_tree(GCProbe::_funcTmpDir());
	my $old = GCProbe::_funcTmpDir() . "KGM/R0/DiaAs.sub.0.KGM.gz";
	touchf($old);
	utime(time - 100000, time - 100000, $old);
	runGCF();
	ok(!-e $old, 'chunk output older than its catalog chunk is discarded');
	ok(scalar(grep { $_->[0] =~ m{/DKGM\.0\.sh$} } @FQ), 'and chunk 0 is realigned');
	#all alignments present, per-gene file and length table missing: the collection must wait for the length job
	touchf("$outD/DIAass_KGM.srt.gz");
	unlink "$dbdir/euk_pro.pep.length";
	runGCF();
	my ($len) = grep { $_->[0] =~ /DBlength_KGM\.sh$/ } @FQ;
	my ($col) = grep { $_->[0] =~ /colDIAKGM\.sh$/ } @FQ;
	ok(defined($len) && defined($col), 'length and collection jobs submitted');
	is($col->[5], 'F1', 'collection depends on the DB length job (first submitted job, id F1)');
	touchf("$dbdir/euk_pro.pep.length");
};
#KEGG/eggNOG get smaller catalog chunks (DIAMOND temp files grow with the chunk), in a split of their own
subtest 'FuncAssign: KEGG/eggNOG use -fastaSplitBigDB' => sub {
	is(GCProbe::_funcSplitSize($_), 4, "$_ uses -fastaSplitBigDB") for qw(KGM KGE KGB NOG);
	is(GCProbe::_funcSplitSize($_), 2, "$_ uses -fastaSplit") for qw(CZy TCDB ABRc VFA);
	my $base = GCProbe::_funcSplitDir(2);
	unlike($base, qr/_2\/$/, 'the -fastaSplit split keeps its directory');
	(my $big = $base) =~ s{/$}{_4/};
	is(GCProbe::_funcSplitDir(4), $big, 'the -fastaSplitBigDB split has its own directory');
	unlink glob("$outD/DIAass_KGM*"), "$outD/.KGM.matrix.done";
	runGCF();
	my @c = chunkJobs();
	is(scalar(@c), 4, 'KGM: one chunk job per -fastaSplitBigDB chunk');
	like($c[0][1], qr/ -q \Q$big\E\S+ /, 'KGM chunks are read from the -fastaSplitBigDB split');
	like(source_sub($gcSrc, 'geneCatFunc'), qr/fastaSplits => \$splitSize/, 'per-database chunk size reaches assignFuncPerGene');
	like($gcSrc, qr/-fastaSplit \$fastaSplits -fastaSplitBigDB \$fastaSplitsBig /, 'main run forwards both chunk sizes to FuncAssign');
	like($gcSrc, qr/my \$doneCmd = "\$rmBin -rf "\.join\(" ", \@splitDirs\)/, 'the final FuncAssign job removes every split');
};
subtest 'module definitions missing: modules skipped, matrix job not failed' => sub {
	my $cmd = Mods::FuncTools::calc_modules("$tmp/x.L0.txt", "$tmp/mods/", 0.5, 0.5, 0);
	is($cmd, '', 'no module command without definitions');
	touchf("$tmp/moddb/module_new.list", "$tmp/moddb/mod.descr", "$tmp/moddb/mod_hiera.txt");
	$cmd = Mods::FuncTools::calc_modules("$tmp/x.L0.txt", "$tmp/mods/", 0.5, 0.5, 0);
	like($cmd, qr/module_new\.list/, 'installed module set is used');
	unlike($cmd, qr/modg\.list|module_s\.list|modGBM\.list/, 'missing module sets skipped');
};

# ---------------- geneCat.pl: in-progress markers of the functional stages ----------------
my $ifl = <<'PKG';
package GCInflight;
use strict; use warnings;
use File::Basename qw(dirname);
use File::Path qw(make_path);
use Fcntl qw(O_CREAT O_EXCL O_WRONLY);
use Sys::Hostname qw(hostname);
our ($QSBoptHR, $GCdir, $live);
sub numLiveUserJobs { return $live; }
PKG
$ifl .= source_sub($gcSrc, $_) . "\n" for qw(_inflightMarker _readInflight _stageJobState _stageInFlight _inflightAcquire _inflightRecordJob _submitStageOnce);
eval "$ifl\n1;" or die $@;
subtest 'functional stages are not submitted twice while in flight' => sub {
	my $bin = "$tmp/fakebin"; make_path($bin);
	write_file("$bin/squeue", "#!/usr/bin/env perl\nprint \$ENV{FAKE_SQUEUE_OUT} // ''; exit(\$ENV{FAKE_SQUEUE_RC} // 0);\n");
	chmod 0755, "$bin/squeue";
	local $ENV{PATH} = "$bin:$ENV{PATH}";
	$GCInflight::QSBoptHR = {qmode => 'slurm', rTag => 'MF_'};
	$GCInflight::GCdir = "$tmp/inflightGC";
	$GCInflight::live = 0;
	my $marker = GCInflight::_inflightMarker('FuncAssign');
	my $calls = 0;
	my $submit = sub { $calls++; return 'MF_4242' };
	GCInflight::_submitStageOnce('FuncAssign', $submit);
	is($calls, 1, 'first submission runs');
	like(read_file($marker), qr/^job\tMF_4242$/m, 'final job recorded in the marker');
	local $ENV{FAKE_SQUEUE_OUT} = "PENDING|Dependency\n";
	GCInflight::_submitStageOnce('FuncAssign', $submit);
	is($calls, 1, 'not submitted again while the final job waits on running jobs');
	$ENV{FAKE_SQUEUE_OUT} = "RUNNING|None\n";
	GCInflight::_submitStageOnce('FuncAssign', $submit);
	is($calls, 1, 'not submitted again while the final job runs');
	$ENV{FAKE_SQUEUE_OUT} = "PENDING|DependencyNeverSatisfied\n";
	GCInflight::_submitStageOnce('FuncAssign', $submit);
	is($calls, 2, 'final job stuck on a failed dependency: stage resubmitted');
	$ENV{FAKE_SQUEUE_OUT} = '';
	GCInflight::_submitStageOnce('FuncAssign', $submit);
	is($calls, 3, 'final job gone: stage resubmitted');
	{
		local $ENV{FAKE_SQUEUE_OUT} = "slurm_load_jobs error: Invalid job id specified\n";
		local $ENV{FAKE_SQUEUE_RC} = 1;
		is(GCInflight::_stageInFlight('FuncAssign'), '', 'purged job id counts as finished');
	}
	{
		local $ENV{FAKE_SQUEUE_OUT} = "slurm_load_jobs error: Socket timed out\n";
		local $ENV{FAKE_SQUEUE_RC} = 1;
		$GCInflight::live = 1;
		GCInflight::_inflightAcquire('FuncAssign'); GCInflight::_inflightRecordJob('FuncAssign', 'MF_4243');
		isnt(GCInflight::_stageInFlight('FuncAssign'), '', 'transient squeue error: falls back to the live-job count');
		$GCInflight::live = 0;
		unlink $marker;
	}
	#submission in progress in another process on this host / a dead submitter
	write_file($marker, "host\t" . Sys::Hostname::hostname() . "\npid\t$$\ntime\t" . time . "\n");
	GCInflight::_submitStageOnce('FuncAssign', $submit);
	is($calls, 3, 'not submitted while another live process is submitting');
	my $pid = fork() // die; if (!$pid) { exit 0 } waitpid($pid, 0);
	write_file($marker, "host\t" . Sys::Hostname::hostname() . "\npid\t$pid\ntime\t" . time . "\n");
	GCInflight::_submitStageOnce('FuncAssign', $submit);
	is($calls, 4, 'marker of a submitter that died is stale');
	unlink $marker;
	eval { GCInflight::_submitStageOnce('FuncAssign', sub { die "submission failed\n" }) };
	like($@, qr/submission failed/, 'submission errors propagate');
	ok(!-e $marker, 'and leave no marker behind');
	GCInflight::_submitStageOnce('FuncAssign', sub { unlink $marker; return 'local' });
	ok(!-e $marker, 'final job already ran (local execution): no marker recreated');
	ok(GCInflight::_inflightAcquire('FuncEMAP') && !GCInflight::_inflightAcquire('FuncEMAP'), 'marker is created exclusively');
};

# ---------------- geneCat.pl stage wiring (static) ----------------
subtest 'geneCat functional stage wiring' => sub {
	my ($stage) = $gcSrc =~ /^(\s*\$stageCmd \.= "\$GCscr -mode FuncAssign[^\n]*)$/m;
	like($stage, qr/-redoFunc \$redoFunc\$funcFwd/, 'FuncAssign stage forwards -redoFunc, tmp dirs and scheduler');
	like($gcSrc, qr/-mode FuncEMAP [^\n]*\$funcFwd/, 'FuncEMAP stage forwards tmp dirs and scheduler');
	like($gcSrc, qr/\$redoFunc \|\| !\(-s \$funcStone && checkpoint_valid/, 'legacy empty functional stone no longer valid');
	like($gcSrc, qr/-s \$emapStone && _stone_valid/, 'legacy empty eggNOG stone no longer valid');
	like($gcSrc, qr/! -name 'funcSplit_\*' ! -name 'GCanno_\*' ! -name 'eggNOGmapper_\*'/, 'final tmp cleanup keeps queued functional inputs');
	unlike($gcSrc, qr/-extHiera -hieraSrtDown/, 'eggNOG-mapper matrices without -extHiera');
	like($gcSrc, qr/\|\| \[ \\\$\? -eq 1 \]; \} > \$tarAnno3\.tmp/, 'eggNOG chunk merge fails on a missing chunk');
	like($gcSrc, qr/push \@funcOuts, "\$GCdir\/Anno\/Func\/\.\$\{curDB\}\.matrix\.done"/, 'functional stone requires the matrix markers');
	like($gcSrc, qr/\$fastaSplits =~ \/\^\\d\+\[MG\]\?\$\//, '-fastaSplit accepts only what splitFastas understands');
	like($gcSrc, qr/_submitStageOnce\('FuncAssign', sub \{/, 'FuncAssign mode submits under its in-progress marker');
	like($gcSrc, qr/_submitStageOnce\('FuncEMAP', sub \{/, 'FuncEMAP mode submits under its in-progress marker');
	like($gcSrc, qr/\$doneCmd \.= "\$rmBin -f "\._inflightMarker\('FuncAssign'\)/, 'FuncAssign final job removes its marker');
	like(source_sub($gcSrc, 'geneCatFunc_emapper'), qr/_inflightMarker\('FuncEMAP'\).*?return \$jobName;/s,
		'eggNOG final job removes its marker and is returned for recording');
	like($gcSrc, qr/if \(\$funcBusy eq "" && \(\$redoFunc/, 'main flow skips a FuncAssign stage in flight');
	like($gcSrc, qr/if \(\$emapBusy eq "" && /, 'main flow skips a FuncEMAP stage in flight');
};

done_testing();
