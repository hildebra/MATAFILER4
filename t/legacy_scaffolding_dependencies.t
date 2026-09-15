use strict;
use warnings;
use Test::More;
use FindBin qw($Bin);
use File::Temp qw(tempdir);
use File::Path qw(make_path);
use lib "$Bin/..";
use Mods::WorkflowControl qw(normalise_job_dependencies append_job_dependencies assembly_cores_for_input);
use Mods::ReadLibrary qw(readLibrariesFromArrays libraryPairs libraryTechnology);
use Mods::IO_Tamoc_progs qw(createGapFillopt);
use Mods::Subm qw(add2SampleDeps);

my $tmp = tempdir(CLEANUP=>1);
sub read_file { open my $fh, '<', $_[0] or die $!; local $/; return <$fh>; }
sub write_file { open my $fh, '>', $_[0] or die $!; print {$fh} $_[1]; close $fh or die $!; }
my $source = read_file("$Bin/../MATAF4.pl");
my %MFopt = (DoAssembly=>0,AssemblyCores=>4,MapperCores=>2,MapperProg=>1,scaffoldMinSize=>500);
my %MFglobal = (shortAssembly=>'');
my %MFconfig = (mateInsertLength=>10000);
my %MFcontstants = (bwt2IdxFileSuffix=>'.bw2');
my %AsGrps;
my $curSmpl = 'sample';
my $logDir = "$tmp/log/";
my $JNUM = 1;
my $QSBoptHR = {tmpSpace=>8};
my $smtBin = 'fixture-samtools';
my $scaffTarExternal = "$tmp/external.fna";
my $scaffTarExternalName = 'ext';
my @scaffTarExternalOLib1 = ("$tmp/gap.R1.fq");
my @scaffTarExternalOLib2 = ("$tmp/gap.R2.fq");
my $libraries = readLibrariesFromArrays(sample=>'sample',scope=>'primary',phase=>'raw',
    technology=>'hiSeq',r1=>["$tmp/mate.R1.fq"],r2=>["$tmp/mate.R2.fq"],labels=>['mate']);
my @jobs;
sub getRawLibrariesAssmGrp { return $libraries; }
sub getCleanLibrariesAssmGrp { return $_[2] ? [] : $libraries; }
sub spaceInAssGrp { return 100; }
sub getProgPaths { return 'fixture-'.$_[0]; }
sub buildMapperIdx { return ("fixture-index $_[0]\n", "$_[0].bw2", "$_[0].bw2.1.bt2"); }
sub qsubSystem { push @jobs, [@_]; return ('run'.(600+@jobs), ''); }
sub megahitAssembly { return 'run200'; }
for my $name (qw(metagAssemblyRun scaffoldCtgs GapFillCtgs)) {
    my ($body) = $source =~ /(^sub \Q$name\E[^\n]*\{.*?^\})/ms;
    die "Missing $name" unless $body;
    eval $body; die $@ if $@;
}
make_path($logDir);
write_file($scaffTarExternal, ">ref\nACGT\n");
my ($cAssGrp, $nodeSpTmpD, $metagAssDir, $geneDir, $SmplNameX, $scaffoldFlag, $metaGscaffDir) =
    ('g', "$tmp/node", "$tmp/metag/", "$tmp/genes/", 'sample', 0, "$tmp/scaff/");
my ($assemblyFlag, $AssemblyGo, $ePreAssmbly, $doPreAssmFlag, $postPreAssmblGo, $finalCommAssDir) =
    (0,1,0,0,0,"$tmp/final/");
sub reset_group {
    @jobs = ();
    %AsGrps = (g=>{CntAimAss=>2,AssemblJobName=>$_[0],UnzpDeps=>'run101;run102;run101',SupportReads=>''});
}
sub run_assembly {
    my $messages = '';
    open my $capture, '>', \$messages or die $!;
    local *STDOUT = $capture;
    return metagAssemblyRun($cAssGrp,"$nodeSpTmpD/ass",$metagAssDir,$geneDir,$SmplNameX,
        $scaffoldFlag,$metaGscaffDir,$assemblyFlag,$AssemblyGo,$ePreAssmbly,$doPreAssmFlag,
        $postPreAssmblGo,$finalCommAssDir);
}

reset_group('run200');
my $external;
my $ok = eval { $external = run_assembly(); 1 };
ok($ok, 'legacy scaffolding and gap filling prepare from a fresh output directory') or diag $@;
is(scalar(@jobs), 2, 'legacy path still submits scaffolding and gap filling');
is($jobs[0][5], 'run200;run101;run102', 'scaffolding waits for assembly and deduplicated raw staging');
is($jobs[1][5], 'run601', 'gap filling waits for the scaffold producer');
is($external, 'run601;run602', 'both legacy consumer IDs return to the controller');
is($AsGrps{g}{AssemblJobName}, 'run200', 'external scaffolding does not serialize unrelated assembly mapping');
like($jobs[0][1], qr{fixture-bwt2.*\Q$tmp/mate.R1.fq\E.*\Q$tmp/mate.R2.fq\E}s,
    'legacy scaffolding command still consumes the raw mate libraries');
like($jobs[1][1] || '', qr/fixture-gapfiller.*Scaffolds_pass2\.fa/s,
    'legacy gap filling command still consumes the scaffold output');
is($QSBoptHR->{tmpSpace}, 8, 'gap filling restores the caller scratch reservation');

@scaffTarExternalOLib1 = (); @scaffTarExternalOLib2 = ();
reset_group('');
is(run_assembly(), 'run601', 'an existing assembly still returns its external scaffold job');
is($jobs[0][5], 'run101;run102', 'published assembly does not erase pending raw staging dependencies');

$scaffTarExternal = '';
reset_group('run200');
is(run_assembly(), '', 'disabled external scaffolding returns no consumer dependencies');
is(scalar(@jobs), 0, 'disabled legacy path submits no extra jobs');
is($AsGrps{g}{AssemblJobName}, 'run200', 'ordinary assembly dependencies are unchanged');

# Ordinary assembly scaffolding shares the staging fix, and remains a producer
# of the assembly that downstream mappings must wait for.
$scaffoldFlag = 1;
reset_group('run200');
is(run_assembly(), '', 'ordinary scaffolding is not classified as external work');
is($jobs[0][5], 'run200;run101;run102', 'ordinary scaffolding also waits for staged raw libraries');
is($AsGrps{g}{AssemblJobName}, 'run200;run601', 'ordinary scaffolding stays in assembly publication dependencies');
$scaffoldFlag = 0;

$scaffTarExternal = "$tmp/external.fna";
@scaffTarExternalOLib1 = ("$tmp/gap.R1.fq"); @scaffTarExternalOLib2 = ("$tmp/gap.R2.fq");
my $original_libraries = $libraries;
$libraries = [];
reset_group('run200');
my $nothing = run_assembly();
is(scalar(@jobs), 0, 'missing eligible group libraries do not schedule gap filling without a scaffold');
is($nothing, 'run200;run101;run102', 'no-work return retains real upstream prerequisites');
$libraries = $original_libraries;

# Verify the main caller adds returned consumer IDs before the next producer
# wave can leave this sample; a return value alone cannot protect cleanup.
my @sampleDeps;
my ($dispatch) = $source =~ /(\t\tmy \$externalScaffoldingDeps = metagAssemblyRun\( \$cAssGrp,.*?\n\t\tadd2SampleDeps[^\n]*;)/s;
ok(defined($dispatch), 'locate the actual main-loop assembly dispatch and dependency capture');
reset_group('run200');
{
    my $messages = '';
    open my $capture, '>', \$messages or die $!;
    local *STDOUT = $capture;
    eval $dispatch; die $@ if $@;
}
is_deeply(\@sampleDeps, ['run601','run602'], 'main loop retains both external consumers for sample waits and cleanup');
is($AsGrps{g}{AssemblJobName}, 'run200', 'capturing cleanup consumers still leaves assembly mapping independent');

done_testing();
