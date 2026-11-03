use strict;
use warnings;
no warnings qw(once redefine prototype); #subroutine stubs for the extracted pipeline code
use Test::More;
use File::Temp qw(tempdir);
use File::Path qw(make_path);
use File::Spec;
use FindBin qw($Bin);
use IO::Compress::Gzip qw(gzip $GzipError);
use lib "$Bin/..", "$Bin/lib";
use MFTestConfig;
use Mods::IO_Tamoc_progs qw(decideMapper buildMapperIdx mapperDBbuilt);
use Mods::ReadLibrary qw(newReadLibrary libraryPairs libraryFiles);

#regression tests for the mapper/search-tool options audit of 8 October 2026
#(docs/audits/2026-10-08/mappers.md); mapping commands are tested in t/alignment_input.t
our $root = File::Spec->rel2abs("$Bin/..");
my $tmp = tempdir(CLEANUP => 1);
$ENV{PERL5OPT} = "-I$root -I$root/t/lib -MMFTestConfig";

sub write_file { my ($p, $c) = @_; open my $fh, '>', $p or die "$p: $!"; print {$fh} $c; close $fh or die $!; }
sub read_file { my ($p) = @_; open my $fh, '<', $p or die "$p: $!"; local $/; my $c = <$fh>; close $fh; return $c; }
sub source_sub {
	my ($src, $name) = @_;
	my ($s) = $src =~ /^(sub \Q$name\E\b[^\n]*\{.*?^\})/ms;
	die "cannot find sub $name\n" unless defined $s;
	return $s;
}
sub executable { my ($p, $c) = @_; write_file($p, $c); chmod 0755, $p; }
my $mfSrc = read_file("$root/MATAF4.pl");

# ---------------- kma removed, mapper indexes ----------------
subtest 'mapper selection and indexes' => sub {
	ok(!eval { decideMapper(4, ''); 1 } && $@ =~ /kma\) is no longer supported/, '-mapper 4 (kma) is refused');
	is(decideMapper($_, ''), $_, "mapper $_ still accepted") for (1, 2, 3, 5);
	like($mfSrc, qr/-mapper 4 \(kma\) was removed/, 'MATAF4 refuses -mapper 4 at startup');
	unlike(read_file("$root/Mods/config_internal.txt"), qr/^kma\t/m, 'no kma program path');

	my $ref = "$tmp/idx/ref.fa"; make_path("$tmp/idx"); write_file($ref, ">r\nACGT\n");
	my ($cmd, $target, $chk) = buildMapperIdx($ref, 2, 0, 3);
	is($cmd, '', 'minimap2: no prebuilt .mmi (it would override the -x preset)');
	is($target, $ref, 'minimap2 maps against the FASTA');
	ok(mapperDBbuilt($ref, 3), 'minimap2 needs no index');

	local *Mods::IO_Tamoc_progs::getProgPaths = sub { return $_[0] eq 'bwt2' ? "$tmp/fakebt2" : 'bwa' };
	write_file("$ref.pac", 'x');
	($cmd, undef, $chk) = buildMapperIdx($ref, 2, 0, 2);
	isnt($cmd, '', 'bwa: an index with only .pac (interrupted build) is rebuilt');
	ok(!mapperDBbuilt($ref, 2), 'bwa: .pac alone is not a complete index');
	is($chk, "$ref.sa", 'bwa completeness is judged by .sa');
	write_file("$ref.sa", 'x');
	($cmd) = buildMapperIdx($ref, 2, 0, 2);
	is($cmd, '', 'bwa: complete index is reused');

	executable("$tmp/fakebt2-build", "#!/bin/sh\ntouch $tmp/bt2.built\n");
	($cmd) = buildMapperIdx($ref, 2, 0, 1);
	write_file("$ref.bw2.$_.bt2l", 'x') for qw(1 2 3 4 rev.1 rev.2);
	is(system('bash', '-c', $cmd), 0, 'bowtie2 index command runs');
	ok(!-e "$tmp/bt2.built", 'an index that bowtie2-build wrote as .bt2l is not rebuilt');
	unlink "$ref.bw2.rev.2.bt2l";
	system('bash', '-c', $cmd);
	ok(-e "$tmp/bt2.built", 'an incomplete index is rebuilt');

	my $out = `$^X $root/secScripts/assemblies/deployMapDB.pl a b c 1 d e 4 2>&1`;
	ok($? != 0 && $out =~ /Mapper must be 1, 2, 3 or 5/, 'deployMapDB.pl takes the mapper of the sample, not kma');
	like(read_file("$root/secScripts/assemblies/deployMapDB.pl"), qr/buildMapperIdx\(\$output_db, \$cores, 0, \$mapper\)/,
		'deployMapDB.pl indexes the decoy database for that mapper');
	like(source_sub($mfSrc, 'scaffoldCtgs'), qr/buildMapperIdx\("\$refCtgs",\$Ncore,0,1\);\s*\}/,
		'scaffolding maps mate libraries with bowtie2, so it builds a bowtie2 index');
};

# ---------------- mate interleaving for DIAMOND ----------------
subtest 'interleaveMates.pl' => sub {
	my $il = "$^X $root/secScripts/functions/interleaveMates.pl";
	my $fq = sub { join '', map { "\@$_\nACGT\n+\nIIII\n" } @_ };
	write_file("$tmp/a1.fq", $fq->('SRR1.1 length=4', 'SRR1.2 length=4'));
	write_file("$tmp/a2.fq", $fq->('SRR1.1 length=4', 'SRR1.2 length=4'));
	is(`$il $tmp/a1.fq $tmp/a2.fq`, $fq->('SRR1.1/1', 'SRR1.1/2', 'SRR1.2/1', 'SRR1.2/2'),
		'mates without a suffix (fasterq-dump) are interleaved and named /1, /2');
	my $gz1 = $fq->('X 1:N:0:ACGT', 'Y/1', 'Z.5.1'); my $gz2 = $fq->('X 2:N:0:ACGT', 'Y/2', 'Z.5.2');
	gzip(\$gz1 => "$tmp/b1.fq.gz") or die $GzipError; gzip(\$gz2 => "$tmp/b2.fq.gz") or die $GzipError;
	is(`$il $tmp/b1.fq.gz $tmp/b2.fq.gz`, $fq->('X/1', 'X/2', 'Y/1', 'Y/2', 'Z.5/1', 'Z.5/2'),
		'gzip input; Illumina comments, /1 /2 and .1 .2 mate names are recognised');
	write_file("$tmp/c2.fq", $fq->('SRR1.2', 'SRR1.1'));
	my $out = `$il $tmp/a1.fq $tmp/c2.fq 2>&1`;
	ok($? != 0 && $out =~ /Mates out of order/, 'mates out of order stop the search');
	write_file("$tmp/d2.fq", $fq->('SRR1.1'));
	$out = `$il $tmp/a1.fq $tmp/d2.fq 2>&1`;
	ok($? != 0 && $out =~ /different numbers of reads/, 'mate files of different length stop the search');
};

# ---------------- MATAF4.pl read-based DIAMOND search ----------------
our (%MFopt, %MFglobal, %map, %progStats, @jobs, %LIBS, $curSmpl, $pigzBin, $logDir, $avx2Constr, $JNUM, %HDDspace, $QSBoptHR);
$curSmpl = 'S'; $pigzBin = 'gzip'; $logDir = "$tmp/log/"; $avx2Constr = ''; $JNUM = 1; %HDDspace = (diamond => 1);
$QSBoptHR = {constraint => [], tmpSpace => 0, General_Hosts => []};
*main::sampleReadSet = sub { return {merged_library => {}} };
*main::getRdLibraries = sub { return %LIBS };
*main::prepDiamondDB = sub { my ($db) = @_; return ("$db.ref", $db, '') };
*main::qsubSystem = sub { push @jobs, [@_]; return ($_[4], $_[1]) };
my %progs;
*main::getProgPaths = sub { return $progs{$_[0]} // "prog_$_[0]" };
eval source_sub($mfSrc, 'runDiamond'); die $@ if $@;
subtest 'read-based DIAMOND options' => sub {
	*main::readLibrariesByScope = sub { return [] };
	my $base = {diaCores => 4, diaRunSensitive => 0, diaFrameshift => 0, diaEVal => '1e-7', DiaPercID => 40,
		DiaMinAlignLen => 20, DiaMinFracQueryCov => 0.1, DiaRmRawHits => 0, globalDiamondDependence => {KGM => 'KGM-1'},
		diamondMem => 7, PABtaxChk => 0};
	my $search = sub {
		my ($opt) = @_;
		%MFopt = (%$base, %$opt); @jobs = (); %progStats = ();
		my $o = "$tmp/rd" . (++$JNUM) . "/diamond/"; make_path($o);
		main::runDiamond($o, "$tmp/db/", "$tmp/scr", 'dep', 'KGM');
		my ($job) = grep { $_->[4] =~ /^_DKGM/ } @jobs;
		return $job->[1];
	};
	%LIBS = (0 => ['s.fq.gz'], 1 => ['r2.fq.gz'], 2 => ['r1.fq.gz']);
	my $cmd = $search->({});
	like($cmd, qr/^prog_interleaveMates_scr r1\.fq\.gz r2\.fq\.gz \| prog_diamond blastx (?!.* -q ).*$/m,
		'both mates of a pair go through one DIAMOND run, read from stdin');
	like($cmd, qr/^prog_diamond blastx .* -q s\.fq\.gz$/m, 'single reads are searched from their file');
	unlike($cmd, qr/\bsort\b/, 'no re-sorting: DIAMOND writes the hits of a pair together');
	like($cmd, qr/ -e 1e-4 /, 'default search e-value 1e-4');
	like($cmd, qr/ --min-orf 25 /, 'short-read ORF filter kept without frameshift alignment');
	$cmd = $search->({diaEVal => '1e-7,1e-3'});
	like($cmd, qr/ -e 1e-3 /, 'the search keeps hits up to the most permissive -DiaParseEvals value');
	$cmd = $search->({diaFrameshift => 15});
	unlike($cmd, qr/--min-orf/, 'frameshift mode keeps DIAMOND\'s own ORF handling (no --min-orf)');
	like($cmd, qr/ -F 15 /, 'frameshift penalty passed');
	%LIBS = (1 => ['r2.fq.gz']);
	ok(!eval { $search->({}); 1 } && $@ =~ /unequal numbers of mate files/, 'a mate file without its partner is refused');
};

# ---------------- MetaPhlAn ----------------
subtest 'MetaPhlAn mapping' => sub {
	my $db = "$tmp/mpdb/mpa_vTest";
	make_path("$tmp/mpdb");
	for my $part (['1.bt2l', 2], ['rev.1.bt2l', 1]) {
		open my $fh, '>', "$db.$part->[0]" or die $!; truncate($fh, $part->[1] * 1024**3) or die $!; close $fh;
	}
	my $fq = sub { join '', map { "\@r$_\nACGT\n+\nIIII\n" } 1 .. $_[0] };
	write_file("$tmp/m1.fq", $fq->(3)); write_file("$tmp/m2.fq", $fq->(3)); write_file("$tmp/ms.fq", $fq->(1));
	my $lib = newReadLibrary(id => 'l', sample => 'S', scope => 'primary', phase => 'clean', technology => 'hiSeq', label => 'l',
		files => {r1 => "$tmp/m1.fq", r2 => "$tmp/m2.fq", single => "$tmp/ms.fq"});
	*main::readLibrariesByScope = sub { return [$lib] };
	executable("$tmp/fakebowtie2", <<'FAKE');
#!/usr/bin/env perl
my ($in, $sam); while (@ARGV) { my $a = shift; $in = shift if $a eq '-U'; $sam = shift if $a eq '-S'; }
my $n = 0; for my $f (split /,/, $in) { open my $fh, '<', $f or die; $n++ while <$fh>; } $n /= 4;
open my $o, '>', $sam or die; print {$o} "r1\t0\tm\t1\t42\t4M\t*\t0\t0\tACGT\tIIII\n"; close $o;
print STDERR "Warning: test\n$n reads; of these:\n  $n (100.00%) were unpaired; of these:\n";
FAKE
	executable("$tmp/fakemetaphlan", "#!/bin/sh\necho \"\$@\" > $tmp/metaphlan.args\nwhile [ \$# -gt 1 ]; do shift; done\necho x > \"\$1\"\n");
	%progs = (bwt2 => "$tmp/fakebowtie2", metPhl2 => "$tmp/fakemetaphlan", metPhl2_db => $db, metPhl2Merge => 'merge');
	%MFopt = (DoMetaPhlan => 4); %MFglobal = (MetaPhlanModernCLI => 1); %map = (S => {inputFileSizeMB => 10});
	@jobs = ();
	eval source_sub($mfSrc, 'metphlanMapping'); die $@ if $@;
	main::metphlanMapping("$tmp/mpwork", "$tmp/mpout/", 'S', 2, '');
	my ($job) = @jobs;
	is($job->[3], '9G', 'memory request covers the bowtie2 index (3 GiB) plus 6 GB');
	like($job->[1], qr/ -U \Q$tmp\E\/m1\.fq,\Q$tmp\E\/m2\.fq,\Q$tmp\E\/ms\.fq 2> /, 'all mates mapped unpaired, as MetaPhlAn does itself');
	unlike($job->[1], qr/ -1 | -2 |\.etxt/, 'no paired mapping, no read count from the scheduler log');
	write_file("$tmp/mp.sh", "set -eo pipefail\n$job->[1]");
	is(system('bash', "$tmp/mp.sh"), 0, 'generated MetaPhlAn job runs (bash submission mode, no .etxt)');
	like(read_file("$tmp/metaphlan.args"), qr/--nreads 7 /, '--nreads counts every mate (3 pairs + 1 single = 7 reads)');

	executable("$tmp/fakemp", "#!/bin/sh\necho MetaPhlAn version 4.2.6\n");
	%progs = (metPhl2 => "[[ -n x ]] && $tmp/fakemp");
	%MFopt = (DoMetaPhlan => 1); %MFglobal = ();
	eval source_sub($mfSrc, 'prepMetaphlan'); die $@ if $@;
	main::prepMetaphlan();
	is($MFopt{DoMetaPhlan}, 4, 'version probe runs in bash ([[ ]] of the env: prologue)');
	ok($MFglobal{MetaPhlanModernCLI}, 'MetaPhlAn 4.2 command line detected');
};

# ---------------- other search-tool calls ----------------
subtest 'gene catalogue, LAMBDA, ABR' => sub {
	my $gc = read_file("$root/secScripts/geneCat.pl");
	like($gc, qr/clusterFNA\(\$FMGfileList\{\$cog\},[^\n]*,\$cogThreads,1,[^\n]*,\$cogMem\)/,
		'marker-gene mmseqs clustering gets the threads and memory of its own job');
	like($gc, qr/\(\$cogThreads, \$cogMem\) = \$submitLocal \? \(\$cogCores, \$memCOG\/2\)/, '... the cogCluster.sh request');
	like(read_file("$root/Mods/FuncTools.pm"), qr/\$idxThreads = \$BlastCores < 2 \? 2 : \$BlastCores;\s*\n\s*\$cmdIdx \.=\s*"\$lambdaBin mkindexn -t \$idxThreads/,
		'lambda3 mkindexn gets at least 2 threads (3.1 rejects -t 1)');
	like(read_file("$root/secScripts/functions/ABRblastFilter2.pl"), qr/\$words\[0\] =~ m\/\\\/2\$\//,
		'ABR: only IDs ending in /2 are mate 2');
};

done_testing();
