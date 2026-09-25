use strict;
use warnings;

# End-to-end tests of helpers/install/get_gtdb.pl. The real script is run in
# separate processes on the stub data in helpers/install/get_gtdb/ (--test
# mode), on hand-made archives (extract mode) and with a fake wget on PATH
# (download mode), always against temporary MATAFILER directories so the
# repository configuration is never modified.

use Cwd qw(getcwd);
use File::Basename ();
use File::Copy qw(copy);
use File::Find ();
use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use IO::Compress::Gzip qw(gzip $GzipError);
use JSON::PP ();
use POSIX ();
use Test::More;

my $root = File::Spec->rel2abs(File::Spec->catdir($Bin, '..'));
my $script = "$root/helpers/install/get_gtdb.pl";
my $stubs = "$root/helpers/install/get_gtdb";
plan skip_all => 'needs GNU tar and gzip' unless have('tar') && have('gzip');

my $tmp = tempdir('get_gtdb_t_XXXXXX', TMPDIR => 1, CLEANUP => 1);

sub have {
	my ($p) = @_;
	return scalar grep { -x "$_/$p" } File::Spec->path;
}

sub slurp {
	my ($f) = @_;
	open my $fh, '<:raw', $f or die "Cannot read $f: $!";
	local $/;
	my $t = <$fh>;
	close $fh;
	return $t;
}

sub spew {
	my ($f, $t) = @_;
	open my $fh, '>:raw', $f or die "Cannot write $f: $!";
	print {$fh} $t;
	close $fh or die "Cannot write $f: $!";
}

sub gz_spew {
	my ($f, $t) = @_;
	gzip(\$t => $f) or die "gzip $f: $GzipError";
}

# Run the script in its own process: stdin from /dev/null (non-interactive),
# MF4DIR/MGTKDIR removed unless given, other %env entries set.
sub run_gtdb {
	my ($env, @args) = @_;
	my $out = "$tmp/.stdout";
	my $err = "$tmp/.stderr";
	my $pid = fork();
	die "fork: $!" unless defined $pid;
	if (!$pid) {
		delete @ENV{qw(MF4DIR MGTKDIR MAMBA_EXE)};
		for my $k (keys %$env) {
			if (defined $env->{$k}) { $ENV{$k} = $env->{$k} } else { delete $ENV{$k} }
		}
		chdir($env->{_cwd} // $tmp) or POSIX::_exit(126);
		delete $ENV{_cwd};
		open STDIN, '<', File::Spec->devnull or POSIX::_exit(127);
		open STDOUT, '>', $out or POSIX::_exit(127);
		open STDERR, '>', $err or POSIX::_exit(127);
		my $prog = $env->{_script} // $script;
		exec($^X, $prog, @args) or POSIX::_exit(127);
	}
	waitpid($pid, 0);
	return { exit => $? >> 8, out => slurp($out), err => slurp($err) };
}

sub tree {
	my ($dir) = @_;
	my %t;
	File::Find::find({ no_chdir => 1, wanted => sub {
		return unless -f $_;
		(my $rel = $_) =~ s{^\Q$dir\E/}{};
		$t{$rel} = slurp($_);
	} }, $dir);
	return \%t;
}

sub stub_members {
	my ($archive) = @_;
	my $d = tempdir(DIR => $tmp);
	system('tar', '-xzf', $archive, '--strip-components', '2', '-C', $d) == 0 or die "tar $archive";
	return tree($d);
}

my $LIN = join('',
	"s__Escherichia coli\t0\td__Bacteria\tp__Pseudomonadota\tc__Gammaproteobacteria\to__Enterobacterales\tf__Enterobacteriaceae\tg__Escherichia\ts__Escherichia coli\n",
	"s__Methanocatella smithii\t0\td__Archaea\tp__Methanobacteriota\tc__Methanobacteria\to__Methanobacteriales\tf__Methanobacteriaceae\tg__Methanocatella\ts__Methanocatella smithii\n");
my $CLUST226 = "s__Escherichia coli\tRS_GCF_009898805.1\ns__Methanocatella smithii\tRS_GCF_945873965.1\n";
my $MGTAX = join('',
	"RS_GCF_009898805.1\td__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli\n",
	"RS_GCF_945873965.1\td__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanocatella;s__Methanocatella smithii\n");
my $BASE = 'https://data.ace.uq.edu.au/public/gtdb/data/releases/';

# A temporary MATAFILER directory with copies of the shipped configuration.
sub fake_mf4 {
	my ($name, %o) = @_;
	my $d = "$tmp/$name";
	make_path("$d/Mods");
	copy("$root/Mods/$_", "$d/Mods/$_") or die "copy $_: $!" for qw(config_DBs.txt config_internal.txt);
	spew("$d/config.txt", $o{config_txt}) if defined $o{config_txt};
	return $d;
}

sub cfg_value {   # the value regex of IO_Tamoc_progs::loadConfigs
	my ($text, $key) = @_;
	for my $l (split /\n/, $text) {
		next if $l =~ /^#/;
		return $1 if $l =~ m/^\Q$key\E\t([^#^\t]+)/;
	}
	return undef;
}

###########################################################################
# 1. r226 (current default) "all" in test mode, split GTDB-Tk package
###########################################################################
my $dbroot = "$tmp/DBs";
my $dl226 = "$tmp/dl226";
my $out226 = "$dbroot/MarkerG/GTDB_r226_MGTK";
my $mf4 = fake_mf4('mf4', config_txt => "MFLRDir\t/opt/mf4/\nDBDir\t\$GTDB_TEST_DBROOT/\t#site DBs\n");
my %env226 = (MF4DIR => "$mf4/", GTDB_TEST_DBROOT => $dbroot);

my $r = run_gtdb(\%env226, 'all', '-v', '226', '-t', $dl226, '-d', $out226, '--test');
is($r->{exit}, 0, 'r226 all --test succeeds') or diag($r->{err});
like($r->{out}, qr/get_gtdb/, 'greeting printed');
like($r->{err}, qr/Run in non-interactive terminal, skipping config update/, 'non-interactive run skips config update');

for my $f (qw(bac120_marker_genes_all_r226.tar.gz ar53_marker_genes_all_r226.tar.gz
	bac120_taxonomy_r226.tsv.gz ar53_taxonomy_r226.tsv.gz bac120_metadata_r226.tsv.gz ar53_metadata_r226.tsv.gz)) {
	ok(-f "$dl226/$f", "downloaded $f");
}
is(slurp("$dl226/bac120_marker_genes_all_r226.tar.gz"), slurp("$stubs/bac_markers.tar.gz"),
	'bacterial marker stub used for bac_markers (label based, not the ".tar" substring)');
is(slurp("$dl226/tk_database_parts/gtdbtk_dummy.tar.gz.part_$_"), slurp("$stubs/gtdbtk_dummy.tar.gz.part_$_"),
	"split part $_ downloaded") for qw(aa ab ac);
ok(!-e "$dl226/gtdbtk_r226_data.tar.gz", 'split mode leaves no concatenated archive');
like(slurp("$dl226/.download.finished"), qr/^\d{4}-\d\d-\d\d \d\d:\d\d:\d\d(?:\.\d{6})?$/, 'download marker holds a Python style timestamp');

my $meta = slurp("$dl226/meta.json");
my $u = "${BASE}/release226/226.0";
like($meta, qr/\A\{\n    "version": 226\.0,\n    "urls": \{\n        "bac_markers": "\Q$u\E\/genomic_files_all\/bac120_marker_genes_all_r226\.tar\.gz",\n/,
	'download meta.json layout (json.dump indent=4, key order, float version)');
like($meta, qr/"tk_database_parts": "\Q$u\E\/auxillary_files\/gtdbtk_package\/split_package\/"\n    \},\n    "date": "[^"]+",\n    "script_version": 0\.1,\n    "summary": "GTDB formatted for MATAFILER using get_gtdb\.pl",\n    "dir": "\Q$dl226\E",\n    "tk": "split"\n\}\z/,
	'download meta.json tail with dir and tk');
my $m = JSON::PP->new->decode($meta);
is_deeply([sort keys %{ $m->{urls} }], [sort qw(bac_markers bac_taxonomy bac_metadata arc_markers arc_taxonomy arc_metadata tk_database tk_database_parts)], 'all URL labels recorded');
is($m->{urls}{tk_database}, "$u/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz", 'full package URL');

my $mg = "$out226/markerGenes";
my %want = (%{ stub_members("$stubs/bac_markers.tar.gz") }, %{ stub_members("$stubs/arc_markers.tar.gz") });
is_deeply([sort grep { !/^\./ } map { s{^.*/}{}r } glob("$mg/*")], [sort(keys(%want), 'GTDBmg.tax')], 'markerGenes holds all bacterial and archaeal markers');
is(slurp("$mg/$_"), $want{$_}, "marker $_ content") for sort keys %want;
ok(-e "$mg/.marker.finished", 'marker completion flag');
ok(!-e "$mg/arc", 'temporary arc directory removed');
is(slurp("$out226/gtdb_r226_lineageGTDB.tab"), $LIN, 'r226 lineage table');
is(slurp("$out226/gtdb_r226_clustering.tab"), $CLUST226, 'r226 clustering table lists every accession');
is(slurp("$mg/GTDBmg.tax"), $MGTAX, 'r226 GTDBmg.tax');
is(slurp("$out226/gtdb/gtdbtk_dummy/aa/aa.txt"), slurp_member('aa/aa.txt'), 'GTDB-Tk split parts extracted (aa)');
ok(-f "$out226/gtdb/gtdbtk_dummy/$_/$_.txt", "GTDB-Tk file $_ extracted") for qw(ab ac);
ok(-e "$out226/gtdb/.tk.finished", 'GTDB-Tk completion flag');
my $dm = slurp("$out226/meta.json");
like($dm, qr/\A\{\n    "version": 226\.0,\n.*\n    "summary": "GTDB formatted for MATAFILER using get_gtdb\.pl"\n\}\z/s, 'extraction meta.json has no dir/tk');

my $readme = slurp("$out226/README.txt");
like($readme, qr/\AWe recommend moving the output directory \(\Q$out226\E\) to your MATAFILER DBDir/, 'README starts with the recommendation');
like($readme, qr/This is in \Q$mf4\E\/Mods\/config_DBs\.txt\.\n/, 'README names MF4DIR/Mods/config_DBs.txt');
like($readme, qr/^GTDBPath\t\[DBDir\]\/MarkerG\/GTDB_r226_MGTK\/markerGenes\/$/m, 'README GTDBPath line');
like($readme, qr/^GTDB_GTDB\t\[DBDir\]\/MarkerG\/GTDB_r226_MGTK\/gtdb_r226_lineageGTDB\.tab$/m, 'README GTDB_GTDB line');
like($readme, qr/^GTDB_lnks\t\[DBDir\]\/MarkerG\/GTDB_r226_MGTK\/gtdb_r226_clustering\.tab$/m, 'README GTDB_lnks line');
like($readme, qr/^get_gtdb\.pl configure -d \Q$out226\E$/m, 'README names get_gtdb.pl configure');
like($readme, qr/^GTDBtk_DB\t\[DBDir\]\/MarkerG\/GTDB_r226_MGTK\/gtdb\/gtdbtk_dummy$/m, 'README GTDBtk_DB points to the release directory');
unlike($readme, qr/GTDBtk_mash/, 'no mash line for GTDB-Tk >= 2.5');
like($readme, qr/between 2\.4\.1 and 2\.6\.1 .*'MF4gtdbtk'.*\n\nmicromamba install --name MF4gtdbtk gtdbtk==2\.6\.1\n\z/s, 'README GTDB-Tk install hint');
like($r->{out}, qr/\Q$readme\E\n\nThese instructions were also written to:\n\Q$out226\E\/README\.txt\n\z/, 'instructions printed');

sub slurp_member {
	my ($m) = @_;
	my $d = tempdir(DIR => $tmp);
	open my $cat, '>:raw', "$d/x.tar.gz" or die;
	print {$cat} slurp("$stubs/gtdbtk_dummy.tar.gz.part_$_") for qw(aa ab ac);
	close $cat;
	system('tar', '-xzf', "$d/x.tar.gz", '-C', $d) == 0 or die 'tar';
	return slurp("$d/gtdbtk_dummy/$m");
}

# Idempotent rerun: same outputs, completed steps are skipped
my $before = tree($out226);
$r = run_gtdb(\%env226, 'all', '-v', '226', '-t', $dl226, '-d', $out226, '--test');
is($r->{exit}, 0, 'rerun succeeds') or diag($r->{err});
like($r->{err}, qr/Download completion marker found/, 'rerun skips downloads');
like($r->{err}, qr/Marker extraction already complete/, 'rerun skips marker extraction');
like($r->{err}, qr/GTDBtk extraction already completed/, 'rerun skips GTDB-Tk extraction');
my $after = tree($out226);
delete @{$_}{qw(meta.json)} for $before, $after;
is_deeply($after, $before, 'rerun leaves identical outputs');

###########################################################################
# 2. configure: config_DBs.txt update in a temporary MATAFILER directory
###########################################################################
my $orig_cfg = slurp("$mf4/Mods/config_DBs.txt");
$r = run_gtdb(\%env226, 'configure', '-d', $out226);
is($r->{exit}, 0, 'configure succeeds') or diag($r->{err});
is(slurp("$mf4/Mods/config_DBs.bup1"), $orig_cfg, 'backup config_DBs.bup1 holds the previous config');
my $new = slurp("$mf4/Mods/config_DBs.txt");
my $tag = '#Updated by get_gtdb.pl';
like($new, qr/^GTDBPath\t\[DBDir\]\/MarkerG\/GTDB_r226_MGTK\/markerGenes\t\Q$tag\E$/m, 'GTDBPath updated with [DBDir]');
like($new, qr/^GTDB_GTDB\t\[DBDir\]\/MarkerG\/GTDB_r226_MGTK\/gtdb_r226_lineageGTDB\.tab\t\Q$tag\E$/m, 'GTDB_GTDB updated');
like($new, qr/^GTDB_lnks\t\[DBDir\]\/MarkerG\/GTDB_r226_MGTK\/gtdb_r226_clustering\.tab\t\Q$tag\E$/m, 'GTDB_lnks updated');
like($new, qr/^GTDBtk_DB\t\[DBDir\]\/MarkerG\/GTDB_r226_MGTK\/gtdb\/gtdbtk_dummy\t\Q$tag\E$/m, 'GTDBtk_DB updated to the release directory');
unlike($new, qr/^GTDBtk_mash/m, 'no active GTDBtk_mash line for GTDB-Tk >= 2.5');
my @orig_l = split /\n/, $orig_cfg, -1;
my @new_l = split /\n/, $new, -1;
is(scalar @new_l, scalar @orig_l, 'line count preserved');
my @changed = grep { $orig_l[$_] ne $new_l[$_] } 0 .. $#orig_l;
my %keys_changed = map { ($new_l[$_] =~ /^#?([^\t]+)/)[0] => 1 } @changed;
ok(!grep({ !/^(?:GTDBPath|GTDB_GTDB|GTDB_lnks|GTDBtk_DB|GTDBtk_mash)$/ } keys %keys_changed), 'only GTDB keys changed')
	or diag(join("\n", map { $new_l[$_] } @changed));
for my $k (qw(GTDBPath GTDB_GTDB GTDB_lnks GTDBtk_DB)) {
	my $v = cfg_value($new, $k);
	ok(defined $v, "$k parseable by the pipeline value regex");
	(my $exp = $v) =~ s/\[DBDir\]/$dbroot\//;
	ok(-e $exp, "$k resolves to an existing path") or diag($exp);
}
$r = run_gtdb(\%env226, 'configure', '-d', $out226);
is($r->{exit}, 0, 'second configure succeeds');
ok(-e "$mf4/Mods/config_DBs.bup2", 'second backup is config_DBs.bup2');
is(slurp("$mf4/Mods/config_DBs.txt"), $new, 'configure is idempotent');

# MGTKDIR fallback, absolute paths, mash restore for GTDB-Tk < 2.5 (r214),
# user config override warning, CRLF input line kept.
my $mf4b = "$tmp/mf4b";
make_path("$mf4b/Mods");
spew("$mf4b/Mods/config_DBs.txt", "#DBs\n#GTDBtk_mash\t/old/mashD\t#Updated by get_gtdb.py\n"
	. "GTDBtk_DB\t/old/gtdbtk\nGTDBPath\t/old/markerGenes/\nOther\tkeep\r\n\nGTDB_GTDB\t/old/lin\t#x\nGTDB_lnks\t/old/cl\n");
spew("$mf4b/config.txt", "DBDir\t/nonexistent/db/\nGTDBPath\t/user/override\n");

###########################################################################
# 3. Older release r214 with the full GTDB-Tk package
###########################################################################
my $dl214 = "$tmp/dl214";
my $out214 = "$tmp/out214";
$r = run_gtdb({ MGTKDIR => $mf4b }, 'all', '-v', '214', '-t', $dl214, '-d', $out214, '--test', '--tk', 'full');
is($r->{exit}, 0, 'r214 all --test --tk full succeeds') or diag($r->{err});
like($r->{err}, qr/untested for 214\.0/, 'untested version warning');
ok(-f "$dl214/gtdbtk_r214_data.tar.gz", 'full GTDB-Tk package downloaded');
ok(!-e "$dl214/tk_database_parts", 'no split parts in full mode');
like(slurp("$dl214/meta.json"), qr/"version": 214\.0,.*"tk": "full"\n\}\z/s, 'r214 meta.json');
like(slurp("$dl214/meta.json"), qr/"bac_markers": "\Q${BASE}\E\/release214\/214\.0\/genomic_files_all\/bac120_marker_genes_all_r214\.tar\.gz"/, 'r214 URLs');
is(slurp("$out214/gtdb_r214_lineageGTDB.tab"), $LIN, 'r214 lineage table');
is(slurp("$out214/gtdb_r214_clustering.tab"), '', 'r214 clustering only lists representatives (none in stubs)');
ok(!-e "$out214/markerGenes/GTDBmg.tax", 'no GTDBmg.tax before r226');
ok(-f "$out214/gtdb/gtdbtk_dummy/ab/ab.txt", 'full GTDB-Tk package extracted');
my $readme214 = slurp("$out214/README.txt");
like($readme214, qr/^GTDBtk_mash\t\[DBDir\]\/MarkerG\/GTDB_r214_MGTK\/gtdb\/gtdbtk_dummy\/mashD$/m, 'mash line for GTDB-Tk < 2.5');
like($readme214, qr/between 2\.1\.0 and 2\.3\.2 .*gtdbtk==2\.3\.2\n\z/s, 'r214 GTDB-Tk version');
like($readme214, qr/This is in \Q$mf4b\E\/Mods\/config_DBs\.txt/, 'README uses MGTKDIR when MF4DIR is unset');

$r = run_gtdb({ MGTKDIR => $mf4b }, 'configure', '-d', $out214);
is($r->{exit}, 0, 'configure via MGTKDIR succeeds') or diag($r->{err});
like($r->{err}, qr/will use absolute paths/, 'database outside DBDir uses absolute paths');
like($r->{err}, qr/config\.txt defines GTDBPath, which override/, 'user config override reported');
my $abs = Cwd::realpath($out214);
is(slurp("$mf4b/Mods/config_DBs.txt"), "#DBs\n"
	. "GTDBtk_mash\t$abs/gtdb/gtdbtk_dummy/mashD\t$tag\n"
	. "GTDBtk_DB\t$abs/gtdb/gtdbtk_dummy\t$tag\n"
	. "GTDBPath\t$abs/markerGenes\t$tag\n"
	. "Other\tkeep\r\n\n"
	. "GTDB_GTDB\t$abs/gtdb_r214_lineageGTDB.tab\t$tag\n"
	. "GTDB_lnks\t$abs/gtdb_r214_clustering.tab\t$tag\n", 'r214 config: absolute paths, mash restored');

# Script-location fallback (no MF4DIR/MGTKDIR): a copy of the script in a
# temporary repository layout updates that layout's config.
my $fake_repo = "$tmp/repo";
make_path("$fake_repo/helpers/install", "$fake_repo/Mods");
copy($script, "$fake_repo/helpers/install/get_gtdb.pl") or die;
spew("$fake_repo/Mods/config_DBs.txt", "GTDBPath\t/old\nGTDBtk_mash\t/old/mashD\n");
$r = run_gtdb({ _script => "$fake_repo/helpers/install/get_gtdb.pl" }, 'configure', '-d', $out226);
is($r->{exit}, 0, 'configure via script location succeeds') or diag($r->{err});
my $abs226 = Cwd::realpath($out226);
is(slurp("$fake_repo/Mods/config_DBs.txt"),
	"GTDBPath\t$abs226/markerGenes\t$tag\n#GTDBtk_mash\t/old/mashD $tag\n",
	'script-location config updated, mash commented once');

###########################################################################
# 4. r220 with --tk skip, then extract of hand-made archives
###########################################################################
my $dl220 = "$tmp/dl220";
$r = run_gtdb({ MF4DIR => $mf4 }, 'all', '-v', '220', '-t', $dl220, '-d', "$tmp/out220", '--test', '--tk', 'skip');
is($r->{exit}, 0, 'r220 --tk skip succeeds (Python crashed here)') or diag($r->{err});
unlike($r->{err}, qr/untested/, 'no untested warning for r220');
ok(!-e "$tmp/out220/gtdb", 'no GTDB-Tk directory with --tk skip');
ok(!grep({ /gtdbtk/ } glob("$dl220/*")), 'nothing GTDB-Tk downloaded');
like(slurp("$tmp/out220/README.txt"), qr/These were not downloaded by this script/, 'README for skipped GTDB-Tk');
is(slurp("$tmp/out220/gtdb_r220_lineageGTDB.tab"), $LIN, 'r220 lineage table');

# Hand-made r220 download directory without completion flag (offline check).
sub make_marker_tar {
	my ($file, $top, %members) = @_;
	my $d = tempdir(DIR => $tmp);
	for my $m (keys %members) {
		my $p = "$d/$top/$m";
		make_path(File::Basename::dirname($p));
		spew($p, $members{$m});
	}
	system('tar', '-czf', $file, '-C', $d, $top) == 0 or die "tar -czf $file";
}
my $dlx = "$tmp/dlx";
make_path($dlx);
make_marker_tar("$dlx/bac120_marker_genes_all_r220.tar.gz", 'bac120_marker_genes_all_r220',
	'faa/M1.faa' => ">B1\nAAA\n", 'fna/M1.fna' => ">B1\nGGG\n", 'faa/B2.faa' => ">B2\nMMM\n");
make_marker_tar("$dlx/ar53_marker_genes_all_r220.tar.gz", 'ar53_marker_genes_all_r220',
	'faa/M1.faa' => ">A1\nCCC\n", 'faa/A2.faa' => ">A2\nKKK\n");
gz_spew("$dlx/bac120_taxonomy_r220.tsv.gz", "G1\td__B;p__P;s__S one\nG2\td__B;p__P;s__S one\nG3\td__B;p__Q;s__S two");
gz_spew("$dlx/ar53_taxonomy_r220.tsv.gz", "\nA1\td__A;p__R;s__S three\n");
my $hdr = "accession\tgtdb_representative\tgtdb_taxonomy\textra\n";
gz_spew("$dlx/bac120_metadata_r220.tsv.gz", $hdr . "G1\tt\td__B;p__P;s__S one\tx\nG2\tf\td__B;p__P;s__S one\tx\nG3\tt\td__B;p__Q;s__S two\tx\n");
gz_spew("$dlx/ar53_metadata_r220.tsv.gz", $hdr . "A1\tt\td__A;p__R;s__S three\tx\n");
spew("$dlx/meta.json", qq({\n    "version": 220.0,\n    "tk": "skip"\n}));

my $outx = "$tmp/outx";
$r = run_gtdb({}, 'extract', '-t', $dlx, '-d', $outx);
is($r->{exit}, 0, 'extract of hand-made r220 archives succeeds') or diag($r->{err});
is(slurp("$outx/markerGenes/M1.faa"), ">B1\nAAA\n>A1\nCCC\n", 'shared marker: archaeal sequences appended to bacterial file');
is(slurp("$outx/markerGenes/M1.fna"), ">B1\nGGG\n", 'bacteria-only marker kept');
is(slurp("$outx/markerGenes/A2.faa"), ">A2\nKKK\n", 'single archaea-only marker moved (Python needed > 1)');
is(slurp("$outx/markerGenes/B2.faa"), ">B2\nMMM\n", 'bacterial marker present');
ok(!-e "$outx/markerGenes/arc", 'arc directory removed');
is(slurp("$outx/gtdb_r220_lineageGTDB.tab"),
	"s__S one\t0\td__B\tp__P\ts__S one\ns__S two\t0\td__B\tp__Q\ts__S two\ns__S three\t0\td__A\tp__R\ts__S three\n",
	'lineage: one row per species, blank lines skipped, missing final newline handled');
is(slurp("$outx/gtdb_r220_clustering.tab"), "s__S one\tG1\ns__S two\tG3\ns__S three\tA1\n",
	'clustering: representatives only, archaeal header dropped');

# Failures: corrupted and missing archives
my $good = slurp("$dlx/bac120_marker_genes_all_r220.tar.gz");
spew("$dlx/bac120_marker_genes_all_r220.tar.gz", "this is not a gzip tarball\n");
$r = run_gtdb({}, 'extract', '-t', $dlx, '-d', "$tmp/outbad");
isnt($r->{exit}, 0, 'corrupted marker archive fails');
like($r->{err}, qr/Extraction of .*bac120_marker_genes_all_r220\.tar\.gz failed/, 'corruption reported');
ok(!-e "$tmp/outbad/markerGenes/.marker.finished", 'no completion flag after failure');
spew("$dlx/bac120_marker_genes_all_r220.tar.gz", $good);
$r = run_gtdb({}, 'extract', '-t', $dlx, '-d', "$tmp/outbad");
is($r->{exit}, 0, 'extract recovers after the archive is replaced') or diag($r->{err});
is(slurp("$tmp/outbad/markerGenes/M1.faa"), ">B1\nAAA\n>A1\nCCC\n", 'recovered output correct');
unlink "$dlx/ar53_taxonomy_r220.tsv.gz";
$r = run_gtdb({}, 'extract', '-t', $dlx, '-d', "$tmp/outmiss");
isnt($r->{exit}, 0, 'missing download fails');
like($r->{err}, qr/Extract mode, but .*ar53_taxonomy_r220\.tsv\.gz not found/, 'missing file reported');
$r = run_gtdb({}, 'extract', '-t', "$tmp/nowhere", '-d', "$tmp/outmiss");
isnt($r->{exit}, 0, 'extract without meta.json fails');
like($r->{err}, qr/No meta\.json found/, 'missing meta.json reported');

# Split package with a corrupted part
my $dlc = "$tmp/dlc";
$r = run_gtdb({}, 'download', '-v', '226', '-t', $dlc, '--test');
is($r->{exit}, 0, 'download subcommand in test mode');
spew("$dlc/tk_database_parts/gtdbtk_dummy.tar.gz.part_ab", "garbage");
$r = run_gtdb({}, 'extract', '-t', $dlc, '-d', "$tmp/outc");
isnt($r->{exit}, 0, 'corrupted split part fails');
like($r->{err}, qr/GTDB-Tk parts failed|tar stopped reading/, 'split extraction failure reported');
ok(!-e "$tmp/outc/gtdb/.tk.finished", 'no GTDB-Tk completion flag');

###########################################################################
# 5. Real download code path with a fake wget
###########################################################################
my $fakebin = "$tmp/fakebin";
make_path($fakebin);
spew("$fakebin/wget", <<'EOF');
#!/usr/bin/env perl
use strict; use warnings; use File::Copy;
my $src = $ENV{FAKE_SRC};
open my $log, '>>', "$src/wget.log"; print {$log} "@ARGV\n"; close $log;
my ($o, $p, $url);
for (my $i = 0; $i < @ARGV; $i++) {
	if ($ARGV[$i] eq '-O') { $o = $ARGV[++$i] } elsif ($ARGV[$i] eq '-P') { $p = $ARGV[++$i] }
	elsif ($ARGV[$i] =~ m{^https?://}) { $url = $ARGV[$i] }
}
if (defined $o) {
	my $name = (split m{/}, $url)[-1];
	if ($ENV{FAKE_FULLY} && -e $o) { print STDERR "The file is already fully retrieved; nothing to do.\n"; exit 1 }
	if (!-e "$src/$name") { print STDERR "ERROR 404: Not Found.\n"; exit 8 }
	copy("$src/$name", $o) or exit 3; exit 0;
}
mkdir $p;
my @parts = map { "gtdbtk_r226_data.tar.gz.part_$_" } qw(aa ab ac);
my $t = time - 100;
for my $i (0 .. $#parts) {
	my $to = "$p/$parts[$i]";
	next if -e $to;   # -nc
	copy("$src/parts/$parts[$i]", $to) or exit 3;
	if ($ENV{FAKE_FAIL_LAST} && $i == $#parts) { truncate($to, 40); utime($t + 50, $t + 50, $to); print STDERR "Connection reset\n"; exit 4 }
	utime($t + $i, $t + $i, $to);
}
exit 0;
EOF
chmod 0755, "$fakebin/wget";
my $src = "$tmp/src";
make_path("$src/parts");
my %stub_name = (bac120_marker_genes_all_r226 => 'bac_markers.tar.gz', ar53_marker_genes_all_r226 => 'arc_markers.tar.gz',
	bac120_taxonomy_r226 => 'bac_taxonomy.tsv.gz', ar53_taxonomy_r226 => 'arc_taxonomy.tsv.gz',
	bac120_metadata_r226 => 'bac_metadata.tsv.gz', ar53_metadata_r226 => 'arc_metadata.tsv.gz');
for my $k (keys %stub_name) {
	my $ext = $stub_name{$k} =~ /markers/ ? '.tar.gz' : '.tsv.gz';
	copy("$stubs/$stub_name{$k}", "$src/$k$ext") or die;
}
copy("$stubs/gtdbtk_dummy.tar.gz.part_$_", "$src/parts/gtdbtk_r226_data.tar.gz.part_$_") or die for qw(aa ab ac);
my %wenv = (PATH => "$fakebin:$ENV{PATH}", FAKE_SRC => $src);

$r = run_gtdb({ %wenv, FAKE_FAIL_LAST => 1 }, 'download', '-v', '226', '-t', "$tmp/dlw");
isnt($r->{exit}, 0, 'failed split-package download fails');
like($r->{err}, qr/Download of .*split_package\/ failed/, 'download failure reported');
ok(!-e "$tmp/dlw/tk_database_parts/gtdbtk_r226_data.tar.gz.part_ac", 'incomplete last part removed');
ok(-e "$tmp/dlw/tk_database_parts/gtdbtk_r226_data.tar.gz.part_$_", "complete part $_ kept") for qw(aa ab);
ok(!-e "$tmp/dlw/.download.finished", 'no completion flag after failed download');
like(slurp("$src/wget.log"), qr/^-c --progress=dot:giga \Q$BASE\E\/release226\/226\.0\/genomic_files_all\/bac120_marker_genes_all_r226\.tar\.gz -O \Q$tmp\E\/dlw\/bac120_marker_genes_all_r226\.tar\.gz$/m,
	'single file wget command (resume with -c)');
like(slurp("$src/wget.log"), qr/^-r -np -nH -nd -R index\.html\* --level=1 -nc -c --progress=dot:giga \S+split_package\/ -P \Q$tmp\E\/dlw\/tk_database_parts$/m,
	'directory wget command (-nc -c)');

$r = run_gtdb({ %wenv, FAKE_FULLY => 1 }, 'download', '-v', '226', '-t', "$tmp/dlw");
is($r->{exit}, 0, 'resumed download succeeds ("fully retrieved" accepted)') or diag($r->{err});
is(slurp("$tmp/dlw/tk_database_parts/gtdbtk_r226_data.tar.gz.part_ac"), slurp("$stubs/gtdbtk_dummy.tar.gz.part_ac"), 'missing part fetched on rerun');
ok(-e "$tmp/dlw/.download.finished", 'download completion flag');
$r = run_gtdb({}, 'extract', '-t', "$tmp/dlw", '-d', "$tmp/outw");
is($r->{exit}, 0, 'extract of the wget download succeeds') or diag($r->{err});
is(slurp("$tmp/outw/gtdb_r226_clustering.tab"), $CLUST226, 'clustering from downloaded files');
ok(-f "$tmp/outw/gtdb/gtdbtk_dummy/ac/ac.txt", 'GTDB-Tk data from downloaded parts');

# pigz code path (the MF4 environments ship pigz): a shim delegating to gzip
my $pigzbin = "$tmp/pigzbin";
make_path($pigzbin);
spew("$pigzbin/pigz", "#!/usr/bin/env perl\nopen my \$l, '>>', \$ENV{PIGZ_LOG}; print {\$l} \"\@ARGV\\n\"; close \$l;\nexec 'gzip', \@ARGV;\n");
chmod 0755, "$pigzbin/pigz";
$r = run_gtdb({ PATH => "$pigzbin:$ENV{PATH}", PIGZ_LOG => "$tmp/pigz.log" }, 'extract', '-t', "$tmp/dlw", '-d', "$tmp/outp");
is($r->{exit}, 0, 'extract with pigz succeeds') or diag($r->{err});
like(slurp("$tmp/pigz.log"), qr/^-d$/m, 'tar used pigz for decompression');
like(slurp("$tmp/pigz.log"), qr/^-dc \S+bac120_metadata_r226\.tsv\.gz$/m, 'tables streamed through pigz');
my ($tw, $tp) = (tree("$tmp/outw"), tree("$tmp/outp"));
delete @{$_}{qw(meta.json markerGenes/.marker.finished gtdb/.tk.finished)} for $tw, $tp;
is_deeply($tp, $tw, 'pigz and gzip paths give identical outputs');

unlink "$src/ar53_metadata_r226.tsv.gz";
$r = run_gtdb(\%wenv, 'download', '-v', '226', '-t', "$tmp/dl404", '--tk', 'skip');
isnt($r->{exit}, 0, 'HTTP error fails the download');
like($r->{err}, qr/ERROR 404/, 'wget output dumped');
unlike($r->{err}, qr/Removing incomplete file/, 'no directory cleanup for single files');

###########################################################################
# 6. Command line handling
###########################################################################
$r = run_gtdb({}, 'download', '-t', "$tmp/x");
is($r->{exit}, 2, 'missing --version is a usage error');
like($r->{err}, qr/the following arguments are required: -v\/--version/, 'usage error message');
$r = run_gtdb({}, 'all', '-v', '226', '--tk', 'half');
is($r->{exit}, 2, 'invalid --tk choice');
$r = run_gtdb({}, 'download', '-v', 'r226');
is($r->{exit}, 2, 'invalid version');
$r = run_gtdb({}, 'bogus');
is($r->{exit}, 2, 'invalid subcommand');
$r = run_gtdb({});
is($r->{exit}, 0, 'no subcommand prints help');
like($r->{out}, qr/usage: get_gtdb\.pl \[-h\] \[--debug\] \{download,extract,configure,all\}/, 'main usage');
$r = run_gtdb({}, 'configure', '-h');
like($r->{out}, qr/usage: get_gtdb\.pl configure \[-h\] \[-d DEST\] \[--install-tk\]/, 'configure usage');
$r = run_gtdb({}, '--debug', 'extract', '-t', "$tmp/dlw", '-d', "$tmp/outdbg");
like($r->{err}, qr/^DEBUG \[/m, '--debug before the subcommand enables debug output');

done_testing();
