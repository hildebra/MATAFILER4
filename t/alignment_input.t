use strict;
use warnings;
use Test::More;
use File::Temp qw(tempdir);
use File::Path qw(make_path);
use File::Basename qw(basename dirname);
use File::Spec;
use File::Find ();
use Cwd qw(abs_path);
use FindBin qw($Bin);
use List::Util qw(max sum);
use lib "$Bin/..", "$Bin/lib";
use MFTestConfig;
use Mods::SampleCompletion qw(completion_request_signature);
use Mods::GenoMetaAss qw(iniCleanSeqSetHR checkSeqTech is3rdGenSeqTech parseSupportReads discoverReadFiles);
use Mods::ReadLibrary qw(newReadLibrary readLibrariesFromArrays syncSeqSetLegacy readLibrariesByScope
    replaceScopeLibraries libraryFiles libraryPairs libraryTechnology);
use Mods::IO_Tamoc_progs qw(decideMapper);
use Mods::WorkflowControl qw(source_input_files missing_input_files commands_are_lightweight_filesystem normalise_job_dependencies);

sub read_file { open my $fh, '<', $_[0] or die $!; local $/; return <$fh>; }
sub write_file { open my $fh, '>', $_[0] or die $!; print {$fh} $_[1]; close $fh or die $!; }
my $sourceText = read_file("$Bin/../MATAF4.pl");
our (%MFconfig, %MFopt, %map, %AsGrps, %MFcontstants, %HDDspace, %make2ndMapDecoy, %map2ndTogRefDB);
our ($curSmpl, $QSBoptHR, $logDir, $JNUM, $pigzBin, $smtBin);
for my $name (qw(_shell_quote _shell_command _staged_read_files_present _validate_sdm_integer_setting
    outfiles_Bam alignmentCacheRequired alignmentFileStamp alignmentCacheIdentity alignmentCacheComplete
    alignmentFastqCommand alignmentCacheCommand alignmentMappingCommand complexGunzCpMv
    discoverSampleInputs seedUnzip2tmp sdmClean getAlgnCmdBase mapReadsToRef)) {
    my ($body) = $sourceText =~ /(^sub \Q$name\E[^\n]*\{.*?^\})/ms;
    ($body) = $sourceText =~ /(^sub sdmClean\(\)\{.*?)(?=^sub mocat_reorder)/ms if $name eq 'sdmClean';
    die "Missing helper $name" unless defined $body;
    $body =~ s/^(sub \Q$name\E)\(\)/$1/;
    eval $body; die "$name: $@" if $@;
}
my $tmp = tempdir(CLEANUP => 1);
my $fastq = "\@read\nACGT\n+\nIIII\n";
my $alignment = "$tmp/input 'quoted'.sam";
write_file($alignment, $fastq);
my $fakeSamtools = "$tmp/samtools";
write_file($fakeSamtools, <<'PY');
#!/usr/bin/env python3
import sys, os, gzip
args=sys.argv[1:]
if args[0]=='cat':
    for f in args[1:]: sys.stdout.buffer.write(open(f,'rb').read())
    sys.exit(0)
assert args[0]=='fastq'
assert args[args.index('-1')+1]=='/dev/null'
assert args[args.index('-2')+1]=='/dev/null'
with open(os.environ['EXTRACTION_LOG'],'a') as f: f.write('extract\n')
data=open(args[-1],'rb').read()
out=args[args.index('-0')+1]
if out=='-': sys.stdout.buffer.write(data)
else:
    with gzip.open(out,'wb') as f: f.write(data)
if os.environ.get('EXTRACTION_FAIL'): sys.exit(17)
PY
chmod 0755, $fakeSamtools;
$ENV{EXTRACTION_LOG} = "$tmp/extractions";
write_file($ENV{EXTRACTION_LOG}, '');
my $fakeMapper = "$tmp/mapper";
write_file($fakeMapper, <<'PY');
#!/usr/bin/env python3
import sys, gzip, os
args=sys.argv[1:]
if '-U' in args: query=args[args.index('-U')+1]
elif '-i' in args: query=args[args.index('-i')+1]
else: query=args[-1]
if query in ['-','--']: data=sys.stdin.buffer.read()
else:
    with gzip.open(query,'rb') if query.endswith('.gz') else open(query,'rb') as f: data=f.read()
with open(os.environ['MAPPING_LOG'],'a') as f: f.write(query+'\n')
sys.stdout.buffer.write(data)
if os.environ.get('MAPPING_FAIL'): sys.exit(18)
PY
chmod 0755, $fakeMapper;
$ENV{MAPPING_LOG} = "$tmp/mappings";
write_file($ENV{MAPPING_LOG}, '');
my %progs = (samtools=>$fakeSamtools, sdm=>abs_path("$Bin/../bin/sdm"));
sub getProgPaths { return $progs{$_[0]} || $fakeMapper; }
my @submissions;
sub qsubSystem { push @submissions, [@_]; return ('job'.scalar(@submissions), ''); }
sub systemW { die 'local staging failed' if system('bash', '-eo', 'pipefail', '-c', $_[0]); }
my %readsets;
sub sampleReadSet { my ($s,$phase,$value)=@_; $readsets{$s}{$phase}=$value if @_>2; return $readsets{$s}{$phase}; }
sub sdmOptSet { return ("$Bin/../data/sdm_PacBio.txt", "$Bin/../data/sdm_PacBio.txt"); }
sub cleaned_primary_libraries_empty { return 0; }
sub check_map_done { return 0; }
sub getRgStr { return 'RG'; }
sub getMapProgNm { return 'testmapper'; }
sub alignPostTreat {
    my ($p,$i,$k)=@_;
    my $file="$p->{nodeTmp}/part.$i.$k.bam";
    my @parts=@{$p->{subBamsAR}};
    $parts[0] = ($parts[0] || '')." $file";
    return (" | cat > $file\n", \@parts);
}
my $runIndex=0;
sub run_command {
    my ($cmd)=@_;
    my $script="$tmp/run".(++$runIndex).'.sh';
    write_file($script, "#!/bin/bash\nset -eo pipefail\n".$cmd);
    return system('bash','-c','bash "$1" >"$2" 2>&1','test',$script,"$script.log") >> 8;
}
sub extraction_count { return scalar(() = read_file($ENV{EXTRACTION_LOG}) =~ /extract/g); }

ok(!alignmentCacheRequired({primary_passes=>1},'primary'), 'one raw mapping can stream');
ok(alignmentCacheRequired({primary_passes=>2},'primary'), 'repeated raw mappings retain cache');
ok(!alignmentCacheRequired({primary_passes=>2,support_passes=>1},'support'), 'support policy is independent');
ok(alignmentCacheRequired({upload=>1},'support'), 'upload needs a file');
ok(alignmentCacheRequired({unfiltered_files=>1},'primary'), 'SDM bypass needs files for clean consumers');
is(outfiles_Bam('/cache','a.BAM'),'/cache/a.unbam.fq.gz','BAM cache naming remains compatible');
is(outfiles_Bam('/cache','a.sam'),'/cache/a.unsam.fq.gz','SAM has a distinct cache name');
is(outfiles_Bam('/cache','a.cram'),'/cache/a.uncram.fq.gz','CRAM has a distinct cache name');
my $cache="$tmp/cache dir/reads.fq.gz";
my ($cmd,$complete)=alignmentCacheCommand($alignment,'',$fakeSamtools,0,$cache);
ok(!$complete,'new cache is incomplete');
is(run_command($cmd),0,'quoted source and destination extract successfully');
my ($signature)=alignmentCacheIdentity($alignment,'',$fakeSamtools);
ok(alignmentCacheComplete($cache,$signature),'published cache passes local validation');
my $count=extraction_count();
is(run_command($cmd),0,'cache reuse succeeds');
is(extraction_count(),$count,'cache reuse skips extraction');
unlink $cache;
is(run_command($cmd),0,'missing cache is recreated despite its marker');
is(extraction_count(),++$count,'missing cache triggers extraction');
write_file($cache,'truncated');
is(run_command($cmd),0,'changed cache is regenerated');
is(extraction_count(),++$count,'changed cache triggers extraction');
write_file($alignment,$fastq.$fastq);
isnt(run_command($cmd),0,'queued job rejects a changed source');
($cmd,$complete)=alignmentCacheCommand($alignment,'',$fakeSamtools,0,$cache);
ok(!$complete,'source change invalidates cache locally');
is(run_command($cmd),0,'newly prepared job refreshes changed source');
unlink "$cache.source.stone";
{
    local $ENV{EXTRACTION_FAIL}=1;
    isnt(run_command($cmd),0,'failed extraction fails cache command');
    ok(!-e "$cache.source.stone",'failed extraction publishes no marker');
    ok(!-e "$cache.partial.gz",'failed extraction removes partial output');
}
write_file($alignment,$fastq);
my $decodeRef="$tmp/decode.fa";
write_file($decodeRef,">decode\nACGT\n"); write_file("$decodeRef.fai","decode\t4\t8\t4\t5\n");
like(alignmentFastqCommand($alignment,$decodeRef,$fakeSamtools,0,'-'),
    qr/'--reference' '\Q$decodeRef\E'/,'samtools receives the decode reference');
my ($before)=alignmentCacheIdentity($alignment,$decodeRef,$fakeSamtools);
write_file($decodeRef,">decode\nACGTACGT\n");
my ($after)=alignmentCacheIdentity($alignment,$decodeRef,$fakeSamtools);
isnt($after,$before,'decode reference participates in cache identity');

# Use the real discovery and staging functions, mocking only job submission.
make_path("$tmp/primary","$tmp/support","$tmp/log","$tmp/staged");
write_file("$tmp/primary/a.sam",$fastq);
write_file("$tmp/primary/a.fq",$fastq);
write_file("$tmp/support/s.cram",$fastq);
$curSmpl='sample'; $JNUM=1; $logDir="$tmp/log/"; $pigzBin='gzip'; $smtBin=$fakeSamtools;
$QSBoptHR={tmpSpace=>0,General_Hosts=>[]};
%MFconfig=(readsRpairs=>0,rawFileSrchStrSingl=>'\\.fq$',rawFileBamSrchSing=>'\\.(?:bam|sam|cram)$',
    rawFileSrchStrXtra1=>'',rawFileSrchStrXtra2=>'',prefSinglFQgreps=>0,doDateFileCheck=>0,
    filterFromSource=>0,splitFastaInput=>0,abortOnEmptyInput=>1,defaultReadLength=>400,defaultReadLengthX=>400,XfirstReads=>-1);
%MFopt=(unzipCores=>1,useUnmapped=>0,sdmCores=>1,sdmMem=>'1G',gzipSDMOut=>1,trimAdapters=>0,SDMlogQualvsLen=>0);
%map=(sample=>{hasPrimaryRds=>1,rddir=>"$tmp/primary",prefix=>'',SeqTech=>'PB',SupportReads=>"PB:$tmp/support",clip=>'',
    cut5pR1=>0,cut5pR2=>0,firstXrdsRd=>0,firstXrdsWr=>0});
my $found=discoverSampleInputs($curSmpl,"$tmp/primary");
is_deeply($found->{primary}{bam},['a.sam'],'primary SAM discovered by inputBAMregex');
is_deeply($found->{support}{bam},["$tmp/support/s.cram"],'support directory discovers CRAM');
$map{sample}{SupportReads}="PB:$tmp/support/s.cram";
$found=discoverSampleInputs($curSmpl,"$tmp/primary");
is_deeply($found->{support}{bam},["$tmp/support/s.cram"],'explicit support CRAM classified as alignment');
sub stage {
    my ($policy)=@_;
    @submissions=();
    seedUnzip2tmp("$tmp/primary",'sample','',"$tmp/node","$tmp/staged",1,'',0,"$tmp/inputs.txt",$policy);
}
stage({primary_passes=>1,support_passes=>1});
my $raw=sampleReadSet('sample','raw');
is(scalar(@{$raw->{libraries}}),3,'mixed FASTQ, primary SAM, support CRAM keep separate libraries');
is($raw->{libraries}[1]{files}{bam},"$tmp/primary/a.sam",'original alignment retained in raw library');
is($raw->{libraries}[2]{scope},'support','support library scope retained');
ok(!-e $raw->{libraries}[1]{files}{single},'direct mode creates no raw FASTQ cache');
ok(!grep($_->[1] =~ /'fastq'/,@submissions),'direct staging submits no extraction');
stage({primary_passes=>2,support_passes=>1});
is(scalar(@submissions),1,'switch to repeated mapping invalidates staging marker');
like($submissions[0][1],qr/'fastq'/,'repeated mapping schedules extraction');
is(run_command($submissions[0][1]),0,'generated mixed staging command executes');
$raw=sampleReadSet('sample','raw');
ok(-s $raw->{libraries}[1]{files}{single},'primary repeated-mapping cache materialized');
ok(!-e $raw->{libraries}[2]{files}{single},'single-pass support still avoids extraction');
$count=extraction_count();
stage({primary_passes=>2,support_passes=>1});
is(scalar(@submissions),0,'completed cache and files reuse staging marker');
unlink "$tmp/staged/rawRds/done.sto";
stage({primary_passes=>2,support_passes=>1});
is(run_command($submissions[0][1]),0,'restaging preserves reusable alignment cache');
is(extraction_count(),$count,'restaging other inputs does not re-extract a valid cache');
stage({primary_passes=>1,support_passes=>1});
ok(sampleReadSet('sample','raw')->{libraries}[1]{metadata}{alignment_cache_required},
    'a valid cache remains selected when only one mapping pass is left');
is(scalar(@submissions),0,'remaining mapping pass reuses existing staging');
$MFconfig{inputCramReference}=$decodeRef;
stage({upload=>1});
is(run_command($submissions[0][1]),0,'upload staging materializes both scopes');
$raw=sampleReadSet('sample','raw');
ok(-s $raw->{libraries}[2]{source_files}{single},'upload receives an existing support FASTQ');

# SDM must prefer original alignment even when a raw cache is also present.
is($raw->{libraries}[2]{metadata}{cram_reference},$decodeRef,'support decode reference inherits the primary setting');
@submissions=();
make_path("$tmp/result");
sdmClean("$tmp/result","$tmp/clean",'upstream',1,1);
like($submissions[0][1],qr{'-i' '\Q$tmp/support/s.cram\E'},'SDM consumes original CRAM');
like($submissions[0][1],qr/'-cramRef' '\Q$decodeRef\E'/,'SDM receives the decode reference');
unlike($submissions[0][1],qr/\.uncram\.fq/,'SDM ignores the raw cache');
like($submissions[0][1],qr/filtered\.suppl\.s\.fq\.gz/,'support clean output name stays compatible');

@submissions=();
sdmClean("$tmp/result","$tmp/clean",'upstream',1,0);
like($submissions[0][1],qr{'-i' '\Q$tmp/primary/a.sam\E'},'mixed primary SDM also selects original alignment');
like($submissions[0][1],qr/filtered\.lib1\.s\.fq\.gz/,'mixed-library filtered output indexing stays compatible');

# Exercise complete mapping command construction and its runtime branches.
my $ref="$tmp/reference.fa"; write_file($ref,">ref\nACGT\n");
write_file("$ref.idx.1.bt2",'index');
%MFcontstants=(bwt2IdxFileSuffix=>'.idx',kmaIdxFileSuffix=>'.kma');
my $lib=newReadLibrary(id=>'alignment',sample=>'sample',scope=>'primary',phase=>'staged',technology=>'PB',
    label=>'same',files=>{single=>"$tmp/nonexistent.fq.gz",bam=>$alignment},metadata=>{alignment_cache_required=>0});
sub mapping {
    my ($mapper,$refs,$libs)=@_;
    $MFopt{MapperProg}=$mapper; $MFopt{MapperMemory}=1; $MFopt{mapModeTogether}=0;
    $MFopt{DoMapModeDecoy}=0; $MFopt{largeMapperDB}=0;
    my @refs=split /,/,$refs;
    my (undef,undef,$params)=mapReadsToRef({smplName=>join(',',map {"sample$_"} 0..$#refs),assGrp=>'group',
        is2ndMap=>0,cramAlig=>0,submNow=>0,unalDir=>'',mapCores=>1,mapSupport=>0,sbj=>$refs,libraries=>$libs,
        glbMapDir=>join(',',map {"$tmp/map$_"} 0..$#refs),nodeTmp=>"$tmp/mapnode",readTec=>'PB',
        glbTmp=>"$tmp/mapwork",outDir=>join(',',map {"$tmp/mapfinal$_"} 0..$#refs)},'stage');
    return $params->{mappingCommand};
}
for my $mapper (1,3,4,5) {
    my $mapping=mapping($mapper,$ref,[$lib]);
    is(run_command($mapping),0,"mapper $mapper complete streaming command executes");
    is(read_file("$tmp/mapwork/sample0.iniAlignment.bam"),$fastq,"mapper $mapper receives raw reads through stdin");
    ok(!-e "$tmp/mapnode_map/alignment.0.fq.gz","mapper $mapper avoids cache on eligible pass");
    unlink "$tmp/mapwork/sample0.iniAlignment.bam";
    {
        local $ENV{EXTRACTION_FAIL}=1;
        isnt(run_command($mapping),0,"mapper $mapper propagates extraction failure");
        ok(!-e "$tmp/mapwork/sample0.iniAlignment.bam","mapper $mapper cannot publish failed extraction");
    }
    {
        local $ENV{MAPPING_FAIL}=1;
        isnt(run_command($mapping),0,"mapper $mapper propagates mapper failure through post-treatment");
        ok(!-e "$tmp/mapwork/sample0.iniAlignment.bam","mapper $mapper cannot publish failed mapping");
    }
}
# A large sparse FASTA exercises the runtime size limit without allocating 1 GB.
open my $large, '>', $ref or die $!;
print {$large} ">large\n";
truncate($large,1000000001) or die $!;
close $large;
my $largeMapping=mapping(3,$ref,[$lib]);
is(run_command($largeMapping),0,'large FASTA uses minimap2 file fallback');
ok(-s "$tmp/mapnode_map/alignment.0.fq.gz",'large FASTA gets replayable reads');
unlink "$tmp/mapwork/sample0.iniAlignment.bam";
write_file($ref,"MMI\0binary-index\n");
my $mapping=mapping(3,$ref,[$lib]);
is(run_command($mapping),0,'minimap2 unknown/binary reference uses file fallback');
ok(-s "$tmp/mapnode_map/alignment.0.fq.gz",'minimap2 fallback materializes replayable FASTQ');
is(read_file("$tmp/mapwork/sample0.iniAlignment.bam"),$fastq,'file fallback receives same raw reads');
unlink "$tmp/mapwork/sample0.iniAlignment.bam";
write_file($ref,">ref\nACGT\n");
my $ref2="$tmp/ref2.fa"; write_file($ref2,">ref2\nACGT\n"); write_file("$ref2.idx.1.bt2",'index');
$count=extraction_count();
$mapping=mapping(1,"$ref,$ref2",[$lib]);
is(run_command($mapping),0,'multiple references map from one extracted cache');
is(extraction_count(),$count+1,'multiple references extract only once');
is(read_file("$tmp/mapwork/sample1.iniAlignment.bam"),$fastq,'second reference reuses raw reads');
unlink "$tmp/mapwork/sample0.iniAlignment.bam", "$tmp/mapwork/sample1.iniAlignment.bam";

my $plain=newReadLibrary(id=>'plain',sample=>'sample',scope=>'primary',phase=>'staged',technology=>'PB',
    label=>'same',files=>{single=>"$tmp/plain.fq"});
write_file("$tmp/plain.fq",$fastq);
$count=extraction_count();
$mapping=mapping(1,$ref,[$plain,$lib]);
is(run_command($mapping),0,'Bowtie2 handles same-labelled FASTQ and alignment libraries separately');
is(read_file("$tmp/mapwork/sample0.iniAlignment.bam"),$fastq.$fastq,'mixed-library mapping retains both inputs');
is(extraction_count(),$count+1,'mixed-library mapping extracts only the alignment library');
unlink "$tmp/mapwork/sample0.iniAlignment.bam";

# Real bundled SDM: verify the actual generated cleaner command, with no cache.
SKIP: {
    skip 'bundled SDM executable unavailable', 4 unless -x $progs{sdm};
    my $seq='ACGTTGCAAGTC' x 34;
    my $qual='I' x length($seq);
    my $sam="$tmp/real.sam";
    write_file($sam,"\@HD\tVN:1.6\tSO:unknown\n\@SQ\tSN:ref\tLN:1000\nreal\t4\t*\t0\t0\t*\t*\t0\t0\t$seq\t$qual\n");
    my $real=newReadLibrary(id=>'real',sample=>'sample',scope=>'primary',phase=>'staged',technology=>'PB',label=>'real',
        files=>{single=>"$tmp/never-created.fq.gz",bam=>$sam});
    sampleReadSet('sample','raw',{libraries=>[$real],samplReadLength=>length($seq)});
    sampleReadSet('sample','clean',iniCleanSeqSetHR(sampleReadSet('sample','raw')));
    @submissions=();
    sdmClean("$tmp/result","$tmp/real-clean",'',1,0);
    is(run_command($submissions[0][1]),0,'real SDM executes generated direct-SAM cleaner command');
    ok(-s "$tmp/real-clean/filtered.s.fq.gz",'real SDM produces expected filtered output');
    ok(-e "$tmp/real-clean/filterDone.stone",'real SDM publishes normal filtering checkpoint');
    ok(!-e "$tmp/never-created.fq.gz",'real SDM needs no extraction cache');
}
done_testing();
