#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../../..";
use Mods::WorkflowControl qw(append_job_dependencies normalise_job_dependencies hybrid_group_ready hybrid_package_complete);
use Mods::GenoMetaAss qw(fileGZe);
use File::Temp qw(tempdir);
use File::Path qw(make_path);
use Mods::Subm qw(add2SampleDeps);
use JSON::PP;

# Execute the actual terminal-empty assembly release branch with recording
# stubs. This inspects orchestration only: no biological tool or scheduler runs.
my $path = shift @ARGV || "$Bin/../../../MATAF4.pl";
open my $fh, '<', $path or die $!;
my $source = do { local $/; <$fh> };
my ($body) = $source =~ /(\t\tif \(\$releaseSharedAssembly\) \{.*?)\n\t\t\$finalizeEmptySample->\(/s;
die 'Cannot locate terminal group release branch' unless $body;
my ($releaseSharedAssembly, $jdep, $SmplName, $cAssGrp, $curOutDir, $curSmpl) =
    (1, '', 'empty', 'g', '/fixture/empty', 'empty');
my ($nodeSpTmpD, $metagAssDir, $geneDir, $metaGscaffDir, $finalCommAssDir) =
    qw(/scratch/empty /assembly/pre /assembly/genes /assembly/scaff /assembly/metag);
my ($assemblyFlag, $AssemblyGo, $ePreAssmbly, $doPreAssmFlag, $postPreAssmblGo) = (1,1,0,0,1);
my ($smplLockF, $JNUM, $QSBoptHR, $logDir) = ('/fixture/lock', 2, {}, '/fixture/log');
my ($binningComplete, $ePreAssmblPck, $smplTmpDir, $finAssLoc, $BinningOut, $cleanedEmpty) =
    (0,0,'/scratch/empty','/assembly/metag/scaffolds.fasta.filt','/assembly/bins',0);
my %MFconfig = (silent=>1);
my (%map, %progStats);
my $closedSample;
my %checkpointNames = (preAssemblyDone=>'pre.done', assemblyDone=>'final.done');
my %MFopt = (DoAssembly=>5,DoMetaBat2=>4,useBinnerScratch=>1);
my %AsGrps = (g => {CntAss=>2,AssemblSmplDirs=>"/fixture/first\n/fixture/empty\n",
    PostAssemblCmd=>'mapping command',PostClnCmd=>'statistics command',PostConsCmd=>'variant command',
    MapDeps=>'',BinDeps=>'',SeqClnDeps=>'run100'});
my (%loopSampleCompleted, @sampleDeps);
my @smplIDs = ('first', 'empty');
my @events;
sub metagAssemblyRun { $AsGrps{g}{AssemblJobName}='run200'; push @events, {job=>'assembly',deps=>$AsGrps{g}{SeqClnDeps}}; return ''; }
sub genePredictions { push @events, {job=>'genes',deps=>$_[2]}; return 'run201'; }
sub postSubmQsub { push @events, {job=>'mapping',deps=>$_[2]}; return 'run202'; }
sub submitGenomeBinner { push @events, {job=>'binning',deps=>$AsGrps{g}{BinDeps},scratch=>$_[0]}; return 'run203'; }
sub MFnext {}
sub loop2C_check {}
my $finalizeEmptySample = sub { my %args=@_; push @events, {job=>'empty-finalizer',%args}; };
my $messages = '';
{
    open my $capture, '>', \$messages or die $!;
    local *STDOUT = $capture;
    eval "for (1) { $body }"; die $@ if $@;
}
my $unreleased_statistics = $AsGrps{g}{PostClnCmd};
my $unreleased_variants = $AsGrps{g}{PostConsCmd};
# Run the actual completed-terminal early return, then visit a ready real
# member. The completed member must remain in hybrid exclusion accounting.
my $tmp = tempdir(CLEANUP=>1);
make_path("$tmp/package");
for my $name (qw(scaffolds.fasta.filt Coverage.percontig Coverage.median.percontig mapping.coverage.gz breakpoints.tsv.gz moved.sto)) {
    open my $out, '>', "$tmp/package/$name" or die $!; print {$out} "ready\n"; close $out;
}
open my $manifest, '>', "$tmp/package/package.manifest.tsv" or die $!;
print {$manifest} "schema_version\t2\n"; close $manifest;
%AsGrps = (g=>{CntAss=>1,CntAimAss=>2,SupportReads=>'ONT:support'});
$closedSample = {outcome=>{status=>'skipped_empty_input'},components=>{}};
my ($closed_body) = $source =~ /(^\tif \(\$closedSample\) \{.*?^\t\})/ms;
die 'Missing completed-sample branch' unless $closed_body;
{
    open my $capture, '>>', \$messages or die $!;
    local *STDOUT = $capture;
    eval "for (1) { $closed_body }"; die $@ if $@;
}
$AsGrps{g}{CntAss}++;
%map = (real=>{hasPrimaryRds=>1,inputFilesEmpty=>0,SupportReads=>'ONT:support'});
$curSmpl = 'real';
my ($prep) = $source =~ /(^sub prepPreAssmbl[^\n]*\{.*?^\})/ms;
die 'Missing prepPreAssmbl' unless $prep;
eval $prep; die $@ if $@;
my @hybrid_state;
{
    open my $capture, '>>', \$messages or die $!;
    local *STDOUT = $capture;
    @hybrid_state = prepPreAssmbl("$tmp/pre", "$tmp/package", "$tmp/map", "$tmp/scratch", "$tmp/stats", 'g', "$tmp/final/fasta", "$tmp/final");
}
print JSON::PP->new->canonical->pretty->encode({events=>\@events,
    completed_empty_member_case=>{visited_members=>2,target_members=>2,
        packages=>$AsGrps{g}{CntPreAss},empty_members_accounted=>$AsGrps{g}{CntPreAssNoPrim},
        final_hybrid_ready=>$hybrid_state[2]},
    sample_dependencies=>\@sampleDeps,
    unreleased_statistics=>$unreleased_statistics,unreleased_variants=>$unreleased_variants,
    trace=>$messages});
