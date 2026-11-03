use strict;
use warnings;
use Test::More;
use File::Temp qw(tempdir);
use File::Path qw(make_path);
use File::Spec;
use FindBin qw($Bin);
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IPC::Open3;
use Symbol qw(gensym);
use lib "$Bin/..", "$Bin/lib";
use MFTestConfig;
use Mods::FuncTools qw(mergeBlastPair);
use Mods::GenoMetaAss qw(coverage_derivative_paths coverage_derivatives_complete);
use Mods::StatsLogReader qw(parse_bam_filter_counters);
use Mods::IO_Tamoc_progs qw(getProgPaths);

my $root = File::Spec->rel2abs("$Bin/..");
my $tmp = tempdir(CLEANUP => 1);
local $ENV{PERL5OPT} = "-I$root -I$root/t/lib -MMFTestConfig";
sub write_file {
    my ($path, $text) = @_;
    open my $fh, '>', $path or die "$path: $!";
    print {$fh} $text;
    close $fh or die "$path: $!";
}
sub read_file {
    open my $fh, '<', $_[0] or die "$_[0]: $!";
    local $/; return <$fh> // '';
}
sub write_gzip { gzip(\$_[1] => $_[0]) or die $GzipError; }
sub read_gzip {
    my $text = '';
    gunzip($_[0] => \$text, MultiStream => 1) or die $GunzipError;
    return $text;
}
sub run_script {
    my ($script, @args) = @_;
    my $err = gensym;
    my $pid = open3(undef, my $out, $err, $^X, "$root/$script", @args);
    my $stdout = do { local $/; <$out> } // '';
    my $stderr = do { local $/; <$err> } // '';
    waitpid($pid, 0);
    return ($? >> 8, $stdout, $stderr);
}
sub hit {
    my ($query, $subject, $length, $score, $start, $end) = @_;
    $length //= 50; $score //= 100; $start //= 1; $end //= $start + $length - 1;
    return [$query, $subject, 90, $length, 0, 0, 1, $length * 3, $start, $end, '1e-20', $score];
}
sub hit_text { return join("\t", @{hit(@_)})."\n"; }
sub functional {
    my ($name, $text, $options) = @_;
    $options //= {};
    my $dir = "$tmp/$name";
    make_path($dir);
    write_gzip("$dir/hits.srt.gz", $text);
    write_file("$dir/lengths", $options->{lengths} // "A\t1000\nB\t100\nC\t10000\n");
    my $db = $options->{db} // 'TEST';
    if ($db eq 'ABRc') {
        write_file("$dir/card.parsed.f11.tab.map", "A\tx\tGeneA\tDrugA\t100\nB\tx\tGeneB\tDrugB\t300\n");
    }
    my ($status, $stdout, $stderr) = run_script('secScripts/functions/parseBlastFunct2.pl',
        '-i', "$dir/hits.srt.gz", '-tmp', "$dir/scratch", '-LF', "$dir/lengths",
        '-DB', $db, '-DButil', "$dir/", '-eval', '1e-7', @{$options->{args} // []});
    if ($options->{fails}) {
        isnt($status, 0, "$name rejects invalid input");
        like($stderr, $options->{fails}, "$name identifies the invalid input");
        ok(!-e "$dir/hits.srt.gz.stone", "$name does not publish success");
        return;
    }
    is($status, 0, "$name completes") or diag("$stdout\n$stderr");
    return {} if $status;
    my %tables;
    my $kind = $db eq 'ABRc' ? 'cat' : 'gene';
    for my $norm (qw(cnt GLN)) {
        $tables{$norm} = read_gzip("$dir/CNT_1e-7_20/${db}parse.$db.ALL.$norm.$kind.cnts.gz");
    }
    return \%tables;
}

my $r = functional('candidate-length', hit_text('r/1','A',50,50).hit_text('r/1','B').hit_text('r/1','C',50,40));
is($r->{GLN}, "B\t0.5\n", 'normalization uses the selected subject, neither first nor last candidate');
$r = functional('candidate-coverage', hit_text('r/1','A',50,50).hit_text('r/1','B'), {args=>['-minPercSbjCov',0.4]});
is($r->{cnt}, "B\t1\n", 'candidate coverage uses its own subject length');
$r = functional('legacy-coverage-alias', hit_text('r/1','A',50,50).hit_text('r/1','B'), {args=>['-minFractQueryCov',0.4]});
is($r->{GLN}, "B\t0.5\n", 'legacy coverage alias retains its subject-coverage meaning');
$r = functional('card-candidate-length', hit_text('r/1','A',50,50).hit_text('r/1','B'), {db=>'ABRc'});
is($r->{cnt}, "GeneA\t1\n", 'CARD rejects B below its own length-scaled bit-score threshold');
for my $mate (1,2) {
    $r = functional("mate-$mate", hit_text("r/$mate",'B'));
    is_deeply($r, {cnt=>"B\t1\n",GLN=>"B\t0.5\n"}, "mate $mate alone contributes one hit");
}
$r = functional('disjoint-mates', hit_text('r/1','A',50,50).hit_text('r/2','B'));
is($r->{cnt}, "B\t1\n", 'different subjects compete independently as singleton candidates');
$r = functional('unsuffixed-single', hit_text('r','B'));
is($r->{cnt}, "B\t1\n", 'unsuffixed single read has multiplicity one');
$r = functional('gene-identifiers', hit_text('gene/1','B').hit_text('gene/2','B'), {args=>['-queryType','genes']});
is($r->{cnt}, "B\t2\n", 'gene IDs that resemble mates are counted as independent genes');
is($r->{GLN}, "B\t1\n", 'independent genes are not overlap-collapsed');
$r = functional('explicit-merged', hit_text('r','B'), {args=>['-queryType','merged']});
is($r->{cnt}, "B\t2\n", 'explicit merged query contributes two reads');
my $mixed = hit_text('a','B').hit_text('b/1','B').hit_text('b/2','B');
(my $merged = hit_text('c','B')) =~ s/\n$/\tMF4:read_count=2\n/;
$r = functional('mixed-provenance', $mixed.$merged);
is($r->{cnt}, "B\t5\n", 'single, combined pair, and merged library retain 1+2+2 counts');
is($r->{GLN}, "B\t1.5\n", 'GLN uses aligned length without multiplying pre-merged overlap');
functional('missing-candidate-length', hit_text('r/1','B').hit_text('r/1','missing'), {fails=>qr/subject length.*missing/});
functional('nonnumeric-candidate-length', hit_text('r/1','B'), {lengths=>"B\tinvalid\n",fails=>qr/subject length.*B/});
functional('zero-candidate-length', hit_text('r/1','B'), {lengths=>"B\t0\n",fails=>qr/subject length.*B/});

for my $case (
    ['contained',1,200,50,100,200], ['identical',1,50,1,50,50],
    ['one-base-overlap',1,50,50,99,99], ['adjacent',1,50,51,100,100],
    ['disjoint',1,50,71,120,100], ['reversed',200,1,100,50,200],
) {
    my ($name,$s1,$e1,$s2,$e2,$expected) = @$case;
    my $a = hit('r/1','B',abs($e1-$s1)+1,abs($e1-$s1)+1,$s1,$e1);
    my $b = hit('r/2','B',abs($e2-$s2)+1,abs($e2-$s2)+1,$s2,$e2);
    my $before = [@$a];
    my $merged = mergeBlastPair($a,$b);
    is($merged->[3], $expected, "$name union length");
    is($merged->[11], $expected, "$name scaled bit score");
    is(mergeBlastPair($b,$a)->[3],$expected,"$name symmetric length");
    is_deeply($a,$before,"$name leaves original evidence intact");
}
my $gapped = mergeBlastPair(hit('r/1','B',210,210,1,200),hit('r/2','B',60,60,50,100));
cmp_ok($gapped->[3],'>=',210,'gapped estimate does not shorten the larger hit');
cmp_ok($gapped->[3],'<=',270,'gapped estimate does not exceed the summed lengths');
$r = functional('paired-containment', hit_text('r/1','B',200,200,1,200).hit_text('r/2','B',51,51,50,100), {lengths=>"B\t200\n"});
is_deeply($r,{cnt=>"B\t2\n",GLN=>"B\t1\n"},'functional caller interprets corrected union and pair multiplicity');

# Exercise the other shared-overlap consumer through its real script.
my $abr = "$tmp/abr"; make_path($abr);
write_file("$abr/ardb.tabs.parsed", "B\tx\tSYM\tCAT\tx\tx\tx\t80\n");
write_file("$abr/ardb_and_reforg_mapping", "x\tSYM\tx\tDRUG\n");
write_file("$abr/ardb_vs_reforg9f.overlap90shortest_famthres_or_symbol.sorted.besthit", "ALT\tB\n");
write_gzip("$abr/hits.srt.gz", hit_text('r/1','B',200,200,1,200).hit_text('r/2','B',51,51,50,100));
my ($status, $stdout, $stderr) = run_script('secScripts/functions/ABRblastFilter2.pl', "$abr/hits.srt.gz", "$abr/genes", "$abr/cats", $abr);
is($status,0,'ABR overlap consumer completes') or diag($stderr);
my @abrFields = split /\t/, read_file("$abr/genes");
is($abrFields[3],200,'ABR reports corrected contained alignment length');
is(0+$abrFields[11],200,'ABR reports corrected contained bit score');

# mOTUs must add an already-known taxon to every sample at all ranks.
my $motus = "$tmp/mOTU2"; make_path($motus);
my $tax = 'd__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__test';
write_gzip("$motus/$_.motu2.tab.gz", "mOTU\tTaxonomy\t$_\nm1\t$tax\t".($_ eq 'A' ? 12 : 3)."\n") for qw(A B);
($status,$stdout,$stderr) = run_script('secScripts/composition/mrgMotu2.pl', $motus, 2);
is($status,0,'two-sample mOTUs merge completes');
for my $rank (qw(kingdom phylum class order family genus species)) {
    like(read_file("$tmp/m2.$rank.txt"),qr/\t12\t3\n\z/,"$rank retains both samples' abundance");
}

# All combinations of producer-empty, header-only and populated hierarchies.
my $hiera = "$tmp/hiera"; make_path($hiera);
write_file("$hiera/empty.hiera.txt",'');
write_gzip("$hiera/header.hiera.txt.gz","read\tdomain\tphylum\n");
write_file("$hiera/full.hiera.txt","read\tdomain\tphylum\nr\tBacteria\tBacillota\n");
($status,$stdout,$stderr) = run_script('secScripts/miTag/miTagTaxTable.pl','domain,phylum',"$tmp/ribo",$hiera);
is($status,0,'empty hierarchies merge successfully') or diag($stderr);
is(read_gzip("$tmp/ribo.domain.txt.gz"),"domain\tempty\tfull\theader\nBacteria\t0\t1\t0\n",'empty samples survive as explicit zero columns');
unlink "$hiera/full.hiera.txt";
($status,$stdout,$stderr) = run_script('secScripts/miTag/miTagTaxTable.pl','domain,phylum',"$tmp/ribo",$hiera);
is($status,0,'all-empty cohort merges successfully') or diag($stderr);
is(read_gzip("$tmp/ribo.phylum.txt.gz"),"phylum\tempty\theader\n",'all-empty cohort retains roster and zero taxa');
write_file("$hiera/bad.hiera.txt","read\twrong-rank\n");
($status,$stdout,$stderr) = run_script('secScripts/miTag/miTagTaxTable.pl','domain',"$tmp/ribo",$hiera);
isnt($status,0,'malformed nonempty hierarchy still fails');

# Window oracle: centered +/- radius, clipped independently at each contig edge.
my $features = "Contig\tAAAA\tCCCC\n";
for my $ctg (['long',24], ['short',2], ['one',1]) {
    $features .= join('', map { "$ctg->[0]_$_\t$_\t".(2*$_)."\n" } 1..$ctg->[1]);
}
write_gzip("$tmp/genes.4kmer.gz",$features);
for my $radius (0,1,5) {
    ($status,$stdout,$stderr) = run_script('secScripts/composition/kmer_Ngenes.pl',"$tmp/genes.4kmer.gz",$radius);
    is($status,0,"gene window radius $radius completes") or diag($stderr);
    my $expected = "Contig\tAAAA\tCCCC\n";
    for my $ctg (['long',24], ['short',2], ['one',1]) {
        for my $i (1..$ctg->[1]) {
            my $lo = $i-$radius < 1 ? 1 : $i-$radius;
            my $hi = $i+$radius > $ctg->[1] ? $ctg->[1] : $i+$radius;
            my $mean = ($lo+$hi)/2;
            $expected .= "$ctg->[0]_$i\t$mean\t".(2*$mean)."\n";
        }
    }
    is(read_gzip("$tmp/genes.4kmer.pm$radius.gz"),$expected,"radius $radius emits every gene once with the centered mean");
}
for my $order (0,1) {
    my @seqs = (">valid\n".('acgt'x30)."\n",">ambiguous\n".('N'x120)."\n");
    @seqs = reverse @seqs if $order;
    write_file("$tmp/kmer.fa",join('',@seqs));
    ($status,$stdout,$stderr) = run_script('secScripts/composition/calc.kmerfreq.pl','-i',"$tmp/kmer.fa",'-o',"$tmp/kmer.$order",'-m',100);
    is($status,0,"uninformative k-mer record in position $order does not crash") or diag($stderr);
    unlike(read_file("$tmp/kmer.$order"),qr/^ambiguous\t/m,'uninformative k-mer record is omitted');
}
is(read_file("$tmp/kmer.0"),read_file("$tmp/kmer.1"),'k-mer emission is independent of final-record position');
write_file("$tmp/gc.fa",">ambiguous\nNNNN\n>lowercase\ngcgc\n>at\natat\n");
($status,$stdout,$stderr) = run_script('secScripts/composition/calcGC.pl',"$tmp/gc.fa","$tmp/gc");
is($status,0,'GC computation completes');
is(read_file("$tmp/gc"),"contig\tGC\nambiguous\t-1\nlowercase\t100.000\nat\t0.000\n",'GC handles missing values, lowercase and newlines');
write_file("$tmp/genes.fa",">ctg_1\nNNN\n>ctg_2\ngcg\n>ctg2_1\nNNN\n");
($status,$stdout,$stderr) = run_script('secScripts/composition/calcGC.pl',"$tmp/genes.fa","$tmp/gc.pergene",'genes');
is($status,0,'GC3 computation completes') or diag($stderr);
is(read_file("$tmp/gc.pergene3"),"contig\tGC\nctg_1\t-1\nctg_2\t100.000\nctg2_1\t-1\n",'GC3 has separate rows and normalized case');
is(read_file("$tmp/gc3"),"contig\tGC\nctg\t100.000\nctg2\t-1\n",'contig GC3 uses exact contig IDs, including prefix-related names');

# Invoke the actual main-script caller bodies with local logs and scheduler state.
my $main = read_file("$root/MATAF4.pl");
sub source_sub {
    my ($name) = @_;
    my ($code) = $main =~ /(^sub \Q$name\E\b[^\{;]*\{.*?^\})/ms;
    die "Cannot isolate $name" unless defined $code;
    return $code;
}
my (%locStats, $logText);
my $QSBoptHR = {}; my $logDir = "$tmp/"; my $JNUM = 0;
sub read_stats_log_excerpt { return $logText; }
eval source_sub('bwtLogRd')."\n".source_sub('getMapStats')."\n".source_sub('calcCoverage2nd');
die $@ if $@;
my $current = "BamFilter\nInput records: 4\nRetained mapped records: 1\nNewly filtered records: 1\nAlready unmapped records: 1\nMalformed SAM records skipped: 1\n";
my $legacy = "Inentries: 4\nTotalRetained: 1\nTotalRm: 3\n";
for my $mapper ('This is strobealign', '[M::worker_pipeline::0.1] mapped 4 sequences') {
    for my $log ($current,$legacy) {
        $logText = "$mapper\n$log";
        my $stats = getMapStats($tmp);
        is($stats->{AlignedReads},1,"$mapper parses retained records");
        is($stats->{OverallAlignment},25,"$mapper uses retained/input percentage");
        is($stats->{UniqueAlgned},'',"$mapper does not invent unique alignment statistics");
    }
}
is_deeply(parse_bam_filter_counters($current.$legacy),{records=>8,retained=>2,filtered=>4,unmapped=>1,malformed=>1},'multiple filter processes sum once each');
ok(!defined(parse_bam_filter_counters("Input records: 4\n")),'partial filter block is unavailable');
ok(!defined(parse_bam_filter_counters('')),'missing filter counters are unavailable');
$logText = "2 reads; of these:\n  2 (100.00%) were unpaired; of these:\n  0 (0.00%) aligned 0 times\n  2 (100.00%) aligned exactly 1 time\n  0 (0.00%) aligned >1 times\n100.00% overall alignment rate\n";
my $stats = getMapStats($tmp);
is($stats->{ReadsPaired},2,'Bowtie total survives without filter counters');
cmp_ok($stats->{OverallAlignment},'==',100,'Bowtie rate remains available without filter counters');
$logText .= "Input records: 2\nRetained mapped records: 1\nNewly filtered records: 1\nAlready unmapped records: 0\nMalformed SAM records skipped: 0\n";
$stats = getMapStats($tmp);
is($stats->{ReadsPaired},2,'Bowtie total survives current filter counters');
is($stats->{OverallAlignment},50,'Bowtie overall rate reflects filtering');

my $cov = "$tmp/sample-smd.bam.coverage.gz";
write_gzip($cov,"sample__C1_L=1000=\t0\t1000\t10\n");
write_file("$tmp/genes.gff","sample__C1_L=1000=\tProdigal\tCDS\t1\t300\t.\t+\t0\tID=1_1;partial=00\n");
ok(!coverage_derivatives_complete($cov),'missing derivatives require work');
my (undef,$command) = calcCoverage2nd($cov,"$tmp/genes.gff",100,'fixture','dependency',{submit=>0});
like($command,qr/\Q$cov\E/,'secondary coverage command reads the compressed input');
SKIP: {
    skip 'bundled calculator requires Linux', 3 unless $^O eq 'linux' && -x "$root/bin/rdCover";
    is(system("$root/bin/rdCover",$cov,"$tmp/genes.gff",100),0,'bundled calculator creates fixture derivatives');
    ok(coverage_derivatives_complete($cov),'completion helper accepts real calculator output');
    my ($dep,$repeat) = calcCoverage2nd($cov,"$tmp/genes.gff",100,'fixture','dependency',{submit=>0});
    is_deeply([$dep,$repeat],['dependency',''],'successful secondary coverage is not rescheduled');
}
for my $suffix (qw(pergene percontig median.percontig)) {
    unlink $_ for @{coverage_derivative_paths($cov,$suffix)};
    write_file("$cov.$suffix",'legacy');
}
ok(coverage_derivatives_complete($cov),'historical derivative spelling remains accepted');
unlink "$cov.pergene";
write_gzip("$cov.pergene.gz",'compressed legacy');
ok(coverage_derivatives_complete($cov),'compressed derivatives remain accepted');
unlink "$cov.percontig";
ok(!coverage_derivatives_complete($cov),'partial legacy products remain incomplete');

# Caller evidence must use the same set of derivative candidates as scheduling.

done_testing();
