use strict;
use warnings;
use Test::More;
use File::Temp qw(tempdir);
use File::Path qw(make_path);
use File::Spec;
use FindBin qw($Bin);
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use lib "$Bin/..", "$Bin/lib";
use MFTestConfig;

my $root = File::Spec->rel2abs("$Bin/..");
my $tmp = tempdir(CLEANUP => 1);
local $ENV{PERL5OPT} = "-I$root -I$root/t/lib -MMFTestConfig";
sub write_file {
    open my $fh, '>', $_[0] or die $!;
    print {$fh} $_[1]; close $fh or die $!;
}
sub read_file {
    open my $fh, '<', $_[0] or die "$! $_[0]";
    local $/; return <$fh>;
}
sub read_gzip {
    my $out = '';
    gunzip($_[0] => \$out, MultiStream=>1) or die $GunzipError;
    return $out;
}
my $fakeDiamond = "$tmp/diamond.pl";
write_file($fakeDiamond, <<'FAKE');
use strict; use warnings; use File::Copy qw(copy);
my ($query,$out);
while (@ARGV) {
    my $arg = shift @ARGV;
    $query = shift @ARGV if $arg eq '-q';
    $out = shift @ARGV if $arg eq '-o';
}
copy($query, "$out.gz") or die "Cannot copy fixture hits: $!";
FAKE
my %libraries;
for my $library ([0,'single'], [1,'pair/2'], [2,'pair/1'], [3,'merged']) {
    my ($key,$query) = @$library;
    my $text = "$query\tB\t90\t50\t0\t0\t1\t150\t1\t50\t1e-20\t100\n";
    my $path = "$tmp/input$key.gz";
    gzip(\$text => $path) or die $GzipError;
    $libraries{$key} = [$path];
}
my $curSmpl = 'sample';
my $pigzBin = Mods::IO_Tamoc_progs::getProgPaths('pigz');
my $logDir = "$tmp/log/";
my $avx2Constr = ''; my $JNUM = 1;
my %HDDspace = (diamond=>1);
my %progStats;
my $QSBoptHR = {constraint=>[],tmpSpace=>0,General_Hosts=>[]};
my %MFopt = (diaCores=>1,diaRunSensitive=>0,diaFrameshift=>0,diaEVal=>'1e-7',
    DiaPercID=>20,DiaMinAlignLen=>30,DiaMinFracQueryCov=>0,DiaRmRawHits=>0,
    globalDiamondDependence=>{TEST=>'TEST-1'},diamondMem=>1);
my @jobs;
sub sampleReadSet { return {merged_library=>{}}; }
sub readLibrariesByScope { return []; }
sub getRdLibraries { return %libraries; }
sub prepDiamondDB { return ('ref','TEST',''); }
sub getProgPaths {
    my ($name) = @_;
    return "$^X $fakeDiamond" if $name eq 'diamond';
    return "$^X $root/secScripts/functions/parseBlastFunct2.pl" if $name eq 'secCogBin_scr';
    return 'unused-mmseqs' if $name eq 'mmseqs2';
    die "Unexpected program lookup: $name";
}
sub qsubSystem { push @jobs, [@_]; return ($_[4], $_[1]); }
my $source = read_file("$root/MATAF4.pl");
my ($runDiamond) = $source =~ /(^sub runDiamond\(\)\{.*?^\})/ms;
ok(defined($runDiamond),'locate actual raw-read search caller');
eval $runDiamond; die $@ if $@;
make_path("$tmp/db","$tmp/diamond");
write_file("$tmp/db/ref.length", "B\t100\n");
runDiamond("$tmp/diamond/","$tmp/db/","$tmp/scratch",'','TEST');
is(scalar(@jobs),2,'caller schedules search and interpretation');
for my $job (@jobs) {
    my $script = "$tmp/$job->[4].sh";
    write_file($script,$job->[1]);
    is(system('bash','-e','-o','pipefail',$script),0,"generated $job->[4] command executes with fixture aligner");
}
my $raw = read_gzip("$tmp/diamond/dia.TEST.blast.srt.gz");
unlike($raw,qr/^single\t.*\tMF4:read_count=/m,'single library uses parser default without rewriting hits');
like($raw,qr/^merged\t.*\tMF4:read_count=2$/m,'merged library records explicit count two');
unlike($raw,qr/^pair\/2\t.*\tMF4:read_count=/m,'unmerged mates use parser default without rewriting hits');
is(read_gzip("$tmp/diamond/CNT_1e-7_20/TESTparse.TEST.ALL.cnt.gene.cnts.gz"),"B\t5\n",
    'caller and parser preserve mixed-library count 1 + 2 + 2');
is(read_gzip("$tmp/diamond/CNT_1e-7_20/TESTparse.TEST.ALL.GLN.gene.cnts.gz"),"B\t1.5\n",
    'caller and parser preserve alignment-length normalization');
ok(-e "$tmp/diamond/dia.TEST.blast.srt.gz.read-counts-v1.stone",'new search records provenance format');
@jobs = ();
runDiamond("$tmp/diamond/","$tmp/db/","$tmp/scratch",'','TEST');
is(scalar(@jobs),0,'completed search with provenance is reused');
my ($prepareRerun) = $source =~ /(^sub prepareDiamondRerun\([^\n]*\{.*?^\})/ms;
eval $prepareRerun; die $@ if $@;
@MFopt{qw(DoDiamond reqDiaDB maxReqDiaDB redoDiamondParse rewriteDiamond)} = (1,'TEST',1,1,0);
prepareDiamondRerun($tmp);
ok(-e "$tmp/diamond/dia.TEST.blast.srt.gz.read-counts-v1.stone",'reparse cleanup preserves search provenance');
@jobs = ();
runDiamond("$tmp/diamond/","$tmp/db/","$tmp/scratch",'','TEST');
is(scalar(@jobs),1,'explicit reparse reuses the marked search without redundant alignment');
unlink "$tmp/diamond/dia.TEST.blast.srt.gz.read-counts-v1.stone";
@jobs = ();
runDiamond("$tmp/diamond/","$tmp/db/","$tmp/scratch",'','TEST');
is(scalar(@jobs),2,'old cached merged-library hits without provenance are regenerated and reinterpreted');
done_testing();
