use strict;
use warnings;
use Test::More;
use File::Spec;
use File::Temp qw(tempdir);
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

#regression tests for gene-catalogue functional annotation:
#eggNOG_split.sh category extraction and parseBlastFunct2.pl hit filtering
my $root = File::Spec->rel2abs('.');
my $tmp = tempdir(CLEANUP => 1);
my $test_lib = File::Spec->catdir($root, 't', 'lib');
local $ENV{PERL5OPT} = join ' ', grep { defined($_) && length($_) }
	"-I$root", "-I$test_lib", '-MMFTestConfig', $ENV{PERL5OPT};

sub write_file { my ($p,$c)=@_; open my $fh,'>',$p or die "$p: $!"; print {$fh} $c; close $fh; }
sub read_file { my ($p)=@_; open my $fh,'<',$p or die "$p: $!"; local $/; my $c=<$fh>; close $fh; return $c; }
sub read_gz { my ($p)=@_; my $c; gunzip($p => \$c) or die "gunzip $p: $GunzipError"; return $c; }
sub two_cols { my ($txt)=@_; return { map { my @s = split /\t/, $_, 2; ($s[0] => $s[1]) } grep { length } split /\n/, $txt }; }

# ---------------- eggNOG_split.sh ----------------
my @hdr = qw(#query seed_ortholog evalue score eggNOG_OGs max_annot_lvl COG_category Description Preferred_name GOs EC KEGG_ko KEGG_Pathway KEGG_Module KEGG_Reaction KEGG_rclass BRITE KEGG_TC CAZy BiGG_Reaction PFAMs);
sub emrow { my %v = @_; my @r = ('-') x 21; $r[0]=$v{q}; $r[4]=$v{og}//'-'; $r[9]=$v{go}//'-'; $r[10]=$v{ec}//'-';
	$r[11]=$v{ko}//'-'; $r[12]=$v{pw}//'-'; $r[13]=$v{mod}//'-'; $r[18]=$v{cazy}//'-'; $r[19]=$v{bigg}//'-'; $r[20]=$v{pfam}//'-'; return join("\t",@r); }
my $emDir = File::Spec->catdir($tmp,'emap'); mkdir $emDir;
my $anno = File::Spec->catfile($emDir,'MF.emapper.annotations');
write_file($anno, join("\n", join("\t",@hdr),
	emrow(q=>'1', og=>'COG0366@1|root,COG0366@2|Bacteria', go=>'GO:0005975', ec=>'3.2.1.1', ko=>'ko:K01176', pw=>'ko00500,map00500,ko01100,map01100', cazy=>'GH13', pfam=>'Alpha-amylase'),
	emrow(q=>'2', og=>'COG0438@1|root', ec=>'2.4.1.-', ko=>'ko:K00754', cazy=>'GT2', pfam=>'Glycos_transf_1,Glycos_transf_2'),
	emrow(q=>'MM2__C1_L=500;_2', og=>'COG0001@1|root', ko=>'ko:K01845,ko:K00002', mod=>'M00121', bigg=>'GSA2')) . "\n");
my $split = File::Spec->catfile($root,'secScripts','GC','eggNOG_split.sh');
is(system('bash', $split, $anno), 0, 'eggNOG_split.sh completes');
my %em = map { $_ => two_cols(read_file(File::Spec->catfile($emDir,"eggNOGmapper_$_.geneAss"))) } qw(CAZy EC GO BIGG PFAM KO KGM KGP NOG);
is_deeply($em{CAZy}, {1=>'GH13', 2=>'GT2'}, 'CAZy families without the digit 1 are kept (old awk filter bug)');
is_deeply($em{GO}, {1=>'GO:0005975'}, 'GO terms kept, "-" dropped');
is_deeply($em{EC}, {1=>'3.2.1.1', 2=>'2.4.1.-'}, 'EC numbers kept');
is_deeply($em{PFAM}, {1=>'Alpha-amylase', 2=>'Glycos_transf_1,Glycos_transf_2'}, 'PFAM names kept');
is_deeply($em{BIGG}, {'MM2__C1_L=500;_2'=>'GSA2'}, 'BiGG kept, non-numeric gene IDs supported');
is_deeply($em{KO}, {1=>'K01176', 2=>'K00754', 'MM2__C1_L=500;_2'=>'K01845,K00002'}, 'KO table without ko: prefix');
is($em{KGP}{1}, 'K01176;ko00500,ko01100', 'KEGG pathways stay comma separated after removing map ids');
is($em{KGM}{'MM2__C1_L=500;_2'}, 'K01845,K00002;M00121', 'KO;module hierarchy');
is($em{NOG}{'MM2__C1_L=500;_2'}, 'COG0001', 'root OG extracted for non-numeric gene IDs');
ok(!exists $em{NOG}{'#query'} && !exists $em{CAZy}{'#query'}, 'header line skipped');
my $gzAnno = "$anno.gz"; gzip($anno => $gzAnno) or die $GzipError; unlink $anno;
is(system('bash', $split, $gzAnno), 0, 'eggNOG_split.sh reads gzipped annotations');
my $bad = File::Spec->catfile($tmp,'bad','MF.emapper.annotations'); mkdir File::Spec->catdir($tmp,'bad');
write_file($bad, "1\tonly_two_columns\n");
isnt(system("bash $split $bad 2>/dev/null"), 0, 'truncated annotation rows are rejected');
ok(!-e File::Spec->catfile($tmp,'bad','eggNOGmapper_CAZy.geneAss'), 'no outputs written for rejected input');

# ---------------- parseBlastFunct2.pl ----------------
my $pb = File::Spec->catfile($root,'secScripts','functions','parseBlastFunct2.pl');
my $lenF = File::Spec->catfile($tmp,'db.length');
write_file($lenF, join("", map {"$_\n"} "x:LONG\t1000","x:SHORT\t320","x:SHORT2\t200","x:HUGE\t5000","x:A\t1000","x:B\t1000","x:C\t320"));
my $n = 0;
sub parse_hits {
	my ($rows, @extra) = @_;
	my $d = File::Spec->catdir($tmp, "pb".(++$n)); mkdir $d; mkdir "$d/tmp";
	my $in = "$d/DIAass_ACL.srt.gz";
	my $txt = join("", map { join("\t",@$_)."\n" } @$rows);
	gzip(\$txt => $in) or die $GzipError;
	my $rc = system("$^X $pb -i $in -DB ACL -mode 1 -singleSpecies 1 -calcGeneLengthNorm 0 -percID 25 -minAlignLen 30 -minBitScore 45 -eval 1e-8 -LF $lenF -tmp $d/tmp/ -summaryTbls 0 @extra >$d/log 2>&1");
	return ($rc, -e "${in}geneAss.gz" ? two_cols(read_gz("${in}geneAss.gz")) : undef);
}
#12 columns (read-based layout): subject coverage only, length of each hit's own subject
my ($rc, $res) = parse_hits([
	[qw(g1 x:LONG 60 300 0 0 1 300 1 300 1e-50 300)],
	[qw(g1 x:SHORT 55 300 0 0 1 300 1 300 1e-45 250)],
	[qw(g2 x:SHORT2 60 150 0 0 1 150 1 150 1e-40 200)],
	[qw(g2 x:HUGE 60 160 0 0 1 160 1 160 1e-45 260)]], '-minPercSbjCov 0.5');
is($rc, 0, 'parser completes on 12-column input');
is_deeply($res, {g1=>'SHORT', g2=>'SHORT2'}, 'subject coverage uses each hit\'s own subject length');
#14 columns (qlen slen): subject OR query coverage
($rc, $res) = parse_hits([
	[qw(p1 x:A 60 90 0 0 1 90 1 90 1e-50 200 100 1000)],
	[qw(p2 x:B 60 100 0 0 1 100 1 100 1e-50 200 400 1000)],
	[qw(p3 x:C 60 300 0 0 1 300 1 300 1e-50 300 1000 320)]], '-minPercSbjCov 0.5 -minPercQueryCov 0.8');
is($rc, 0, 'parser completes on 14-column input');
is_deeply($res, {p1=>'A', p3=>'C'}, 'hit passes on query coverage OR subject coverage');
($rc, $res) = parse_hits([
	[qw(p1 x:A 60 90 0 0 1 90 1 90 1e-50 200 100 1000)],
	[qw(p3 x:C 60 300 0 0 1 300 1 300 1e-50 300 1000 320)]], '-minPercSbjCov 0.5');
is_deeply($res, {p3=>'C'}, 'without -minPercQueryCov only subject coverage applies');

# ---------------- geneCat.pl structure ----------------
my $gc = read_file(File::Spec->catfile($root,'secScripts','geneCat.pl'));
like($gc, qr/-mode FuncAssign [^\n]*-stone \$funcStone/, 'functional stone is written by the FuncAssign job graph');
unlike($gc, qr/_checkpoint_command\(\$checkpointWriter, \$funcStone/, 'no functional stone at submission time');
like($gc, qr/\.emapper\.annotations"; #the file eggNOG-mapper writes/, 'eggNOG-mapper resume checks the real output file');
like($gc, qr/eggNOGmapper_KO/, 'eggNOG-mapper KO table is summarised');

# ---------------- VFDB (VFA / VFB) ----------------
my $vfDir = File::Spec->catdir($tmp,'VIRDB'); mkdir $vfDir;
my $hA1 = '>VFG000076(gb|NP_460360) (ssaQ) type III secretion system protein SsaQ [TTSS (SPI-2 encode) (VF0036) - Effector delivery system (VFC0086)] [Salmonella enterica]';
my $hA2 = '>VFG037176(gb|WP_001081735) (plc1) phospholipase C [Phospholipase C (VF0470) - Exotoxin (VFC0235)] [Acinetobacter baumannii ACICU]';
my $hB1 = '>VFG001234(gb|WP_000000001) (fimA) type 1 fimbrial major subunit FimA [Type 1 fimbriae (VF0221) - Adherence (VFC0001)] [Escherichia coli CFT073]';
my $hB2 = '>VFG009999(gi:12345) (xyz) old style header [Some VF (VF9999)] [Old organism]';
write_file("$vfDir/VFDB_setA_pro.fas", "$hA1\nMSTAAAA\n$hA2\nMKKKKKK\n");
write_file("$vfDir/VFDB_setB_pro.fas", "$hA1\nMSTAAAA\n$hA2\nMKKKKKK\n$hB1\nMFFFFFF\n$hB2\nMGGGGGG\n");
is(system("$^X ".File::Spec->catfile($root,'secScripts','functions','prepVFDB.pl')." $vfDir >/dev/null"), 0, 'prepVFDB.pl builds VF.tab');
my %vf = map { my @c = split /\t/; ($c[0] => \@c) } grep { length } split /\n/, read_file("$vfDir/VF.tab");
is_deeply([@{$vf{'VFG000076(gb|NP_460360)'}}[1,3,4,6]], ['ssaQ_VF0036','VF0036_TTSS_(SPI-2_encode)','VFC0086_Effector_delivery_system','A'], 'VF.tab: gene, VF and category levels, set A member');
is($vf{'VFG001234(gb|WP_000000001)'}[6], 'B', 'VF.tab: set B only entry flagged');
is($vf{'VFG009999(gi:12345)'}[4], 'VFC_unclassified', 'VF.tab: header without category handled');
my $vfLen = "$vfDir/VFDB_setB_pro.fas.length";
write_file($vfLen, join("", map {"$_\t400\n"} keys %vf));
{
	my $d = File::Spec->catdir($tmp,'vfa'); mkdir $d; mkdir "$d/tmp";
	my $in = "$d/DIAass_VFA.srt.gz";
	my $txt = join("", map { join("\t",@$_)."\n" }
		[qw(v1), 'VFG000076(gb|NP_460360)', qw(85 380 0 0 1 380 1 380 1e-60 400 390 400)],
		[qw(v2), 'VFG037176(gb|WP_001081735)', qw(45 380 0 0 1 380 1 380 1e-30 200 390 400)]);
	gzip(\$txt => $in) or die $GzipError;
	my $rc = system("$^X $pb -i $in -DB VFA -mode 2 -singleSpecies 1 -calcGeneLengthNorm 0 -percID 60 -minAlignLen 50 -minBitScore 60 -eval 1e-10 -minPercSbjCov 0.7 -minPercQueryCov 0.8 -LF $vfLen -DButil $vfDir/ -tmp $d/tmp/ -summaryTbls 0 >$d/log 2>&1");
	is($rc, 0, 'parser runs in VFA mode');
	my $vfa = two_cols(read_gz("${in}geneAss.gz"));
	is_deeply($vfa, {v1 => "ssaQ_VF0036\tVF0036_TTSS_(SPI-2_encode)\tVFC0086_Effector_delivery_system"}, 'VFA: 3-level hierarchy, 45% identity hit rejected');
}
like($gc, qr/VFA => \{percID => 60, minPercSbjCov => 0\.7, minPercQueryCov => 0\.8/, 'VFDB-specific cutoffs defined');

done_testing();
