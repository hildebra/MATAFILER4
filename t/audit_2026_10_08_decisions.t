use strict;
use warnings;
no warnings qw(once redefine prototype); #subroutine stubs for the extracted pipeline code
use Test::More;
use File::Temp qw(tempdir);
use File::Path qw(make_path);
use File::Spec;
use FindBin qw($Bin);
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use lib "$Bin/..", "$Bin/lib";
use MFTestConfig;
use Mods::ReadLibrary qw(newReadLibrary libraryPairs libraryFiles);

#regression tests for the decisions taken on the mapper audit of 8 October 2026
#(docs/audits/2026-10-08/mappers.md, "Decided and implemented")
our $root = File::Spec->rel2abs("$Bin/..");
my $tmp = tempdir(CLEANUP => 1);
$ENV{PERL5OPT} = "-I$root -I$root/t/lib -MMFTestConfig";

sub write_file { my ($p, $c) = @_; open my $fh, '>', $p or die "$p: $!"; print {$fh} $c; close $fh or die $!; }
sub read_file { my ($p) = @_; open my $fh, '<', $p or die "$p: $!"; local $/; my $c = <$fh>; close $fh; return $c; }
sub read_gzip { my $out = ''; gunzip($_[0] => \$out, MultiStream => 1) or die $GunzipError; return $out; }
sub source_sub {
	my ($src, $name) = @_;
	my ($s) = $src =~ /^(sub \Q$name\E\b[^\n]*\{.*?^\})/ms;
	die "cannot find sub $name\n" unless defined $s;
	return $s;
}
sub executable { my ($p, $c) = @_; write_file($p, $c); chmod 0755, $p; }
sub counts { my %c = map { split /\t/ } split /\n/, read_file($_[0]); return \%c; }
my $mfSrc = read_file("$root/MATAF4.pl");

# ---------------- Kraken2 ----------------
#taxo.k2d as kraken2 writes it: header, 56-byte nodes {parent, first_child, child_count, name_offset,
#rank_offset, external_id, godparent}, NUL-terminated names and ranks
sub write_k2d {
	my ($path, @nodes) = @_; #[external, parent external, rank, name]
	my %internal = map { $nodes[$_][0] => $_ + 1 } 0 .. $#nodes;
	my ($names, $ranks, $nodeData) = ('', '', pack('Q<7', (0) x 7));
	for my $n (@nodes) {
		my ($ext, $parent, $rank, $name) = @{$n};
		$nodeData .= pack('Q<7', $parent ? $internal{$parent} : 0, 0, 0, length($names), length($ranks), $ext, 0);
		$names .= "$name\0"; $ranks .= "$rank\0";
	}
	open my $fh, '>:raw', $path or die $!;
	print {$fh} 'K2TAXDAT', pack('Q<3', scalar(@nodes) + 1, length($names), length($ranks)), $nodeData, $names, $ranks;
	close $fh;
}
my $kdb = "$tmp/krakenDB/testDB"; make_path($kdb);
write_k2d("$kdb/taxo.k2d", [1, 0, 'no rank', 'root'], [2, 1, 'superkingdom', 'Bacteria'], [1239, 2, 'phylum', 'Bacillota'],
	[91061, 1239, 'class', 'Bacilli'], [1385, 91061, 'order', 'Bacillales'], [186817, 1385, 'family', 'Bacillaceae'],
	[1386, 186817, 'genus', 'Bacillus'], [1423, 1386, 'species', 'Bacillus subtilis'], [1396, 1386, 'species', 'Bacillus cereus'],
	[224308, 1423, 'strain', 'Bacillus subtilis 168']);
my $kraken2Out = join('', map { join("\t", @{$_}) . "\n" }
	['C', 'r1', 1423, 150, '1423:8 1386:1 0:1'],           #species at every threshold
	['C', 'r2', 1423, 150, '1423:2 1396:2 2:2 0:4'],      #species up to 0.2, genus at 0.3
	['C', 'r3', 224308, '150|150', '224308:3 |:| A:7'],   #pair called at a strain: counted at its species
	['U', 'r4', 0, 150, '0:10'],
	['C', 'r5', 2, 150, '2:1 0:9']);                      #domain only; unclassified from 0.2 on
my $sp = 'd__Bacteria;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus_subtilis';
my $genus = 'd__Bacteria;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;?';
my $domain = 'd__Bacteria;?;?;?;?;?;?';

subtest 'Kraken2 counts per confidence threshold' => sub {
	write_file("$tmp/k2.out", $kraken2Out);
	make_path("$tmp/kc");
	my $out = `$^X $root/secScripts/composition/krak2_count_tax.pl -db $kdb -thresholds 0.01,0.2,0.3 -out $tmp/kc/krak $tmp/k2.out 2>&1`;
	is($?, 0, 'counting script runs') or diag $out;
	is_deeply(counts("$tmp/kc/krak.0.01.cnt.tax"), {$sp => 3, $domain => 1}, '0.01: three species reads, one domain-only read');
	is_deeply(counts("$tmp/kc/krak.0.2.cnt.tax"), {$sp => 3}, '0.2: the domain-only read loses its call (as kraken2 --confidence 0.2)');
	is_deeply(counts("$tmp/kc/krak.0.3.cnt.tax"), {$sp => 2, $genus => 1}, '0.3: a read with 2 of 10 k-mers on its species moves to the genus');

	write_file("$tmp/S1.0.3.krak.txt", "$sp\t2\n$genus\t1\n");
	write_file("$tmp/S2.0.3.krak.txt", "$sp\t5\n");
	is(system("$^X $root/secScripts/composition/mrgKrakTax.pl .0.3.krak.txt $tmp/Krak.0.3.mat $tmp/S1.0.3.krak.txt $tmp/S2.0.3.krak.txt"), 0,
		'cohort merge runs');
	is(read_file("$tmp/Krak.0.3.mat"), "Taxon\tS1\tS2\n$genus\t1\t0\n$sp\t2\t5\n", 'merged matrix: lineages x samples');
};

our (%MFopt, %MFglobal, %map, %HDDspace, $curSmpl, $logDir, $JNUM, $krakDeps, $QSBoptHR, @jobs, %progs);
$curSmpl = 'S'; $logDir = "$tmp/log/"; $JNUM = 1; $krakDeps = ''; $QSBoptHR = {constraint => [], tmpSpace => 0, General_Hosts => []};
*main::qsubSystem = sub { push @jobs, [@_]; return ("job" . scalar(@jobs), $_[1]) };
*main::getProgPaths = sub { return $progs{$_[0]} // "prog_$_[0]" };
*main::sampleReadSet = sub { return {merged_library => {}} };

subtest 'Kraken2 profiling job' => sub {
	open my $h, '>', "$kdb/hash.k2d" or die $!; truncate($h, 2 * 1024**3) or die $!; close $h;
	write_file("$tmp/fake.out", $kraken2Out);
	executable("$tmp/kraken2", "#!/bin/sh\necho \"\$@\" >> $tmp/kraken2.args\ncat $tmp/fake.out\n");
	%progs = (kraken2 => "$tmp/kraken2", krakCnts_scr => "$^X $root/secScripts/composition/krak2_count_tax.pl");
	%MFopt = (krakenCores => 3, globalKraTaxkDB => 'testDB'); %MFglobal = (krakenDBDirGlobal => "$tmp/krakenDB/");
	my @libs = (newReadLibrary(id => 'p', sample => 'S', scope => 'primary', phase => 'clean', technology => 'hiSeq', label => 'p',
		files => {r1 => "$tmp/r1.fq.gz", r2 => "$tmp/r2.fq.gz", single => "$tmp/s.fq"}));
	*main::readLibrariesByScope = sub { return \@libs };
	eval source_sub($mfSrc, 'krakenTaxEst'); die $@ if $@;
	@jobs = ();
	main::krakenTaxEst("$tmp/kout", "$tmp/ktmp/", 'S', '');
	my ($job) = @jobs;
	is($job->[3], '6G', 'memory request: hash.k2d (2 GiB) + 4 GB');
	like($job->[1], qr{kraken2 --db \Q$tmp\E/krakenDB//testDB --threads 3 --confidence 0 --paired --gzip-compressed \Q$tmp\E/r1\.fq\.gz \Q$tmp\E/r2\.fq\.gz\n},
		'paired library: kraken2 at confidence 0, gzip input');
	like($job->[1], qr{--confidence 0 \Q$tmp\E/s\.fq\n}, 'singletons: kraken2 without --paired');
	unlike($job->[1], qr/kraken-filter|kraken-translate|--preload|--fastq-input/, 'no Kraken 1 options or tools');
	write_file("$tmp/kraken.sh", "set -eo pipefail\n$job->[1]");
	is(system('bash', "$tmp/kraken.sh"), 0, 'generated Kraken2 job runs');
	is_deeply(counts("$tmp/kout/krak.0.3.cnt.tax"), {$sp => 4, $genus => 2}, 'counts summed over both libraries');
	ok(-e "$tmp/kout/krakDone.sto" && -e "$tmp/kout/krak.0.01.cnt.tax", 'all thresholds written, sample closed');
	like($mfSrc, qr/getProgPaths\("mrgKrak_scr"\)/, 'cohort tables merged by mrgKrakTax.pl, not merge_metaphlan_tables.py');
};

# ---------------- long reads in read-based DIAMOND ----------------
subtest 'long reads: DIAMOND range culling and per-range assignment' => sub {
	%progs = ();
	our %LIBS = (0 => ['long.fq.gz', 'short.fq.gz']);
	my @libs = ({files => {single => 'long.fq.gz'}, is_long => 1}, {files => {single => 'short.fq.gz'}, is_long => 0});
	*main::readLibrariesByScope = sub { return \@libs };
	*main::getRdLibraries = sub { return %LIBS };
	*main::prepDiamondDB = sub { my ($db) = @_; return ("$db.ref", $db, '') };
	our (%progStats, $pigzBin, $avx2Constr);
	$pigzBin = 'gzip'; $avx2Constr = ''; %HDDspace = (diamond => 1);
	%MFopt = (diaCores => 4, diaRunSensitive => 0, diaFrameshift => 0, diaEVal => '1e-7', DiaPercID => 40, DiaMinAlignLen => 20,
		DiaMinFracQueryCov => 0.1, DiaRmRawHits => 0, globalDiamondDependence => {KGM => 'KGM-1'}, diamondMem => 16);
	eval source_sub($mfSrc, 'runDiamond'); die $@ if $@;
	@jobs = (); make_path("$tmp/rd/diamond");
	main::runDiamond("$tmp/rd/diamond/", "$tmp/db/", "$tmp/scr", 'dep', 'KGM');
	my ($search) = grep { $_->[4] =~ /^_DKGM/ } @jobs;
	like($search->[1], qr/^prog_diamond blastx (?!.* -k )(?!.*--min-orf).* --long-reads -d .* -q long\.fq\.gz$/m,
		'long reads: --long-reads (range culling, --top 10, -F 15), no -k 5, no --min-orf');
	like($search->[1], qr/^prog_diamond blastx .* --min-orf 25 -k 5 -d .* -q short\.fq\.gz$/m, 'short single reads unchanged');
	like($search->[1], qr/MF4:ranges=1/, 'long-read hits are tagged for the parser');
	is($search->[3], '16G', 'DIAMOND read jobs request 16 GB by default');
	like($mfSrc, qr/\$MFopt\{diamondMem\} = 16;/, '-DiaMem default 16');

	#parser: one read with two genes (B and C); B2 overlaps B and scores lower
	my $dbd = "$tmp/pdb"; make_path($dbd);
	write_file("$dbd/ref.length", "B\t100\nB2\t100\nC\t100\n");
	my @hits = ("L1\tB\t90\t100\t0\t0\t1\t300\t1\t100\t1e-30\t200", "L1\tB2\t85\t95\t0\t0\t10\t290\t1\t95\t1e-25\t150",
		"L1\tC\t80\t100\t0\t0\t900\t600\t1\t100\t1e-20\t120");
	for my $case (['ranges', "\tMF4:ranges=1", {B => 1, C => 1}], ['plain', '', {B => 1}]) {
		my ($name, $tag, $expected) = @{$case};
		my $d = "$tmp/parse_$name"; make_path($d);
		my $text = join('', map { "$_$tag\n" } @hits);
		gzip(\$text => "$d/dia.TEST.blast.srt.gz") or die $GzipError;
		my $log = `$^X $root/secScripts/functions/parseBlastFunct2.pl -i $d/dia.TEST.blast.srt.gz -DB TEST -eval 1e-7 -percID 20 -minAlignLen 30 -minPercSbjCov 0 -mode 0 -queryType reads -LF $dbd/ref.length -reportDomains 0 -DButil $dbd/ -tmp $d/tmp 2>&1`;
		is($?, 0, "parser runs ($name)") or diag $log;
		is_deeply({map { split /\t/ } split /\n/, read_gzip("$d/CNT_1e-7_20/TESTparse.TEST.ALL.cnt.gene.cnts.gz")}, $expected,
			$name eq 'ranges' ? 'a tagged long read counts each gene it holds' : 'untagged hits: one assignment per read');
	}
};

# ---------------- TaxaTarget, mOTUs ----------------
subtest 'TaxaTarget' => sub {
	my $tt = "$tmp/taxaTarget"; make_path("$tt/run_pipeline_scripts", "$tt/data");
	write_file("$tt/$_", "x\n") for qw(kaijux diamond data/marker_geneDB.fasta.kaiju.fmi data/phylogroup_total_mgLen.txt);
	write_file("$tt/run_pipeline_scripts/environment.txt", "# paths\nkaiju='$tt/kaijux'\ndiamond='$tt/diamond'\ntaxatarget='$tt/'\n");
	#stand-in: records the first read name it gets; TT_CASE selects success, "no reads mapped" or a silent failure
	write_file("$tt/run_pipeline_scripts/run_protist_pipeline_fda.py", <<'FAKE');
use IO::Uncompress::Gunzip qw($GunzipError);
my %a; while (@ARGV) { my $k = shift; $a{$k} = $k eq '--tmp' ? 1 : shift; }
my $z = IO::Uncompress::Gunzip->new($a{'-r'}) or die $GunzipError; my $h = <$z>;
open my $o, '>>', "$a{'-o'}/../headers.txt" or die; print {$o} $h; close $o;
if ($ENV{TT_CASE} eq 'noreads') { print STDERR "No reads mapped to the marker genes with Kaiju. Analysis ended!\n"; exit 1; }
exit 0 if $ENV{TT_CASE} eq 'silent';
open my $r, '>', "$a{'-o'}/Taxonomic_report.txt" or die; print {$r} "Lineage\tTaxa\tRank\tRead_count\tAbundance\tBusco_count\n"; close $r;
FAKE
	%progs = (TaxaTarget => "$^X $tt/run_pipeline_scripts/run_protist_pipeline_fda.py");
	%MFopt = (DoTaxaTarget => 1); %map = (S => {inputFileSizeMB => 100});
	our $pigzBin = 'gzip';
	eval source_sub($mfSrc, $_) for qw(taxaTargetDir prepTaxaTarget TaxaTarget); die $@ if $@;
	ok(eval { main::prepTaxaTarget(); 1 }, 'complete install passes the startup check') or diag $@;
	unlink "$tt/data/phylogroup_total_mgLen.txt";
	ok(!eval { main::prepTaxaTarget(); 1 } && $@ =~ /phylogroup_total_mgLen\.txt/, 'missing data file (no profile, exit 0 in TaxaTarget) stops MF4');
	write_file("$tt/data/phylogroup_total_mgLen.txt", "x\n");
	%progs = (TaxaTarget => 'python run_protist_pipeline_fda.py');
	ok(!eval { main::prepTaxaTarget(); 1 } && $@ =~ /full path|path of/, 'command without the install path is refused');
	%progs = (TaxaTarget => "$^X $tt/run_pipeline_scripts/run_protist_pipeline_fda.py");

	my $fq = "\@r1/1\nACGT\n+\nIIII\n";
	for my $f (qw(t1 t2 ts)) { my $z = $fq; gzip(\$z => "$tmp/$f.fq.gz") or die; }
	my $lib = newReadLibrary(id => 'l', sample => 'S', scope => 'primary', phase => 'clean', technology => 'hiSeq', label => 'l',
		files => {r1 => "$tmp/t1.fq.gz", r2 => "$tmp/t2.fq.gz", single => "$tmp/ts.fq.gz"});
	*main::readLibrariesByScope = sub { return [$lib] };
	for my $case (['ok', 0], ['noreads', 0], ['silent', 5]) {
		my ($name, $exit) = @{$case};
		my $out = "$tmp/tt_$name/"; make_path($out);
		@jobs = ();
		main::TaxaTarget("$tmp/tt_tmp_$name/", $out, 'S', 2, '');
		my ($job) = @jobs;
		is($job->[3], '6G', "memory request 6G ($name)") if $name eq 'ok';
		write_file("$tmp/tt_$name.sh", "set -eo pipefail\n$job->[1]");
		local $ENV{TT_CASE} = $name;
		is(system("bash $tmp/tt_$name.sh 2>/dev/null") >> 8, $exit, "$name: job exit $exit");
		if ($name eq 'ok') {
			is(read_file("$out/S/headers.txt"), "\@r1\n\@r1\n", 'paired and single-end runs; read names without /1 (kaiju vs extractor)');
			ok(-e "$out/S.TaxTar.sto", 'profile complete');
		} elsif ($name eq 'noreads') {
			ok(-e "$out/S.TaxTar.sto" && -e "$out/S/lib0/no_reads_mapped.txt", 'no protist reads: valid empty result, sample closes');
		} else {
			ok(!-e "$out/S.TaxTar.sto", 'exit 0 without Taxonomic_report.txt is a failure');
		}
	}
};

subtest 'mOTUs memory' => sub {
	my $db = "$tmp/motusdb"; make_path("$db/db_mOTU");
	open my $h, '>', "$db/db_mOTU/mOTUsv4.1.db.fna.gz.bwt" or die $!; truncate($h, 12 * 1024**3) or die $!; close $h;
	%progs = (motus2_DB => $db, motus2 => 'motus');
	%MFopt = (DoMOTU2 => 1); %map = (S => {inputFileSizeMB => 100});
	my $lib = newReadLibrary(id => 'l', sample => 'S', scope => 'primary', phase => 'clean', technology => 'hiSeq', label => 'l',
		files => {r1 => "$tmp/m1.fq", r2 => "$tmp/m2.fq"});
	*main::readLibrariesByScope = sub { return [$lib] };
	eval source_sub($mfSrc, 'mOTU2Mapping'); die $@ if $@;
	@jobs = ();
	main::mOTU2Mapping("$tmp/motutmp", "$tmp/motuout/", 'S', 4, '');
	is($jobs[0][3], '18G', 'mOTUs: bwa index (12 GiB) + 6 GB');
	unlink "$db/db_mOTU/mOTUsv4.1.db.fna.gz.bwt";
	@jobs = ();
	main::mOTU2Mapping("$tmp/motutmp", "$tmp/motuout/", 'S', 4, '');
	is($jobs[0][3], '16G', 'mOTUs: at least 16 GB');
};

# ---------------- other decisions ----------------
subtest 'mosaic loci, geneCat DIAMOND resources' => sub {
	like(read_file("$root/secScripts/MGS/prepare_mosaic_loci.pl"), qr/'-x', \$DEFAULT\{minimap_preset\}, '-s', 40,/,
		'mosaic self-alignment: -s 40 after the asm preset');
	like(read_file("$root/Mods/FuncTools.pm"), qr/--compress 1 --quiet \$diaSensFlag-t /, 'geneCat DIAMOND chunks get the recorded sensitivity');
	like(read_file("$root/secScripts/geneCat.pl"), qr/my \$funcDiaSensitivity = "mid-sensitive";/, 'geneCat uses --mid-sensitive');
};

done_testing();
