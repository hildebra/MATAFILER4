use strict;
use warnings;
no warnings 'once';
use FindBin qw($Bin);
use File::Spec;
use File::Temp qw(tempdir);
use File::Path qw(make_path remove_tree);
use File::Copy qw(copy);
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Test::More;
use lib "$Bin/lib";
use MFTestConfig;

my $root = File::Spec->rel2abs("$Bin/..");
sub read_file {
    open my $in, '<', $_[0] or die "$_[0]: $!";
    local $/;
    return <$in>;
}
sub write_file {
    open my $out, '>', $_[0] or die "$_[0]: $!";
    print {$out} $_[1] or die $!;
    close $out or die $!;
}
sub unzip {
    my $text = '';
    gunzip($_[0] => \$text, MultiStream => 1) or die $GunzipError;
    return $text;
}
my $source = read_file("$root/secScripts/geneCat.pl");
my $imports = <<'IMPORTS';
package GCRecovery;
use strict;
use warnings;
use File::Path qw(make_path remove_tree);
use File::Spec;
use Fcntl qw(O_CREAT O_EXCL O_WRONLY);
use Errno qw(EEXIST);
use IO::Handle;
use English;
use Mods::GenoMetaAss qw(gzipopen systemW readFasta);
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::Checkpoint qw(write_checkpoint checkpoint_valid read_checkpoint);
use Mods::WorkflowResilience qw(retry_unlink retry_rename);
use Mods::geneCat qw(readGeneIdxSpl attachProteins3);
our ($primaryClusterFNA, $primaryClusterCLS, $pigzBin, $catBin, $cpBin, $rmBin,
    $mkdirBin, $mvBin, $tmpDir, $cdhID, $numCor0, $totMem, $submitLocal, $COGdir,
    $GCdir, $GCscr, $checkpointWriter, $QSBoptHR, $clustMMseq, $avx2Constr,
    $useGTDBmg, $mapF, $qsubDir, $rareBin, $countMatrixF, $rtkFunDelims,
    $mmseqs2Bin, $cdhitBin, $GLBtmp, $headBin, $tailBin, $clusterCov);
IMPORTS
my @helpers = qw(_shell_quote _checkpoint_command _stone_valid _sync_file
    _for_each_fasta_record _safe_reset_dir _append_file_locked _reset_collation_outputs
    gzifelscat _catalog_backup_command _merged_catalog_backup_valid clusterFNA
    clusterSingleStep rewriteClusNumbers addCOGgenes mergeClsSam rewriteFastaHdIdx
    combineClstr krakenTax geneCatFunc_emapper);
for my $name (@helpers) {
    my ($helper) = $source =~ /^(sub \Q$name\E\b[^\n]*\{.*?^\})/ms;
    die "Missing $name" unless defined $helper;
    $imports .= "$helper\n";
}
eval $imports;
die $@ if $@;
$GCRecovery::pigzBin = GCRecovery::getProgPaths('pigz');
$GCRecovery::catBin = 'env GENECAT_CAT_WRAPPER=1 cat';
$GCRecovery::cpBin = 'cp';
$GCRecovery::rmBin = 'rm';
$GCRecovery::mkdirBin = 'mkdir';
$GCRecovery::mvBin = 'mv';
$GCRecovery::cdhID = 97;
$GCRecovery::clusterCov = 0.9;
$GCRecovery::primaryClusterFNA = 'compl.incompl.97.fna';
$GCRecovery::primaryClusterCLS = 'compl.incompl.97.fna.clstr';
$GCRecovery::COGdir = 'COG';

# A real child process runs the production marker merge from generated commands.
if (@ARGV && $ARGV[0] eq '--merge') {
    shift @ARGV;
    my %options = @ARGV;
    die 'configured geneCat wrapper was lost' unless $ENV{GENECAT_CHILD_WRAPPER};
    $GCRecovery::GCdir = $options{'-o'};
    GCRecovery::mergeClsSam($options{'-tmp'}, $options{'-clusterID'}, $options{'-o'});
    exit 0;
}

my $tmp = tempdir(CLEANUP => 1);
subtest 'interrupted collation resets owned aggregates before replay' => sub {
    my ($bucket, $scratch, $logs) = map { "$tmp/$_" } qw(B0 scratch logs);
    make_path($bucket, "$scratch/COG", $logs);
    my $part = "$scratch/compl.fna.gz.0";
    gzip(\">s__gene\nACGT\n" => $part) or die $GzipError;
    GCRecovery::_append_file_locked($part, "$bucket/compl.fna.gz", "$scratch/compl.fna.gz.lock");
    write_file("$scratch/compl.fna.gz.lock", "interrupted worker\n");
    write_file("$scratch/COG/preclus.marker.fna.0", ">s__marker\nATG\n");
    write_file("$logs/GeneCompleteness.txt.0", "partial\n");
    write_file("$logs/Missed_samples.txt", "stale\n");
    write_file("$logs/inmap.txt", "keep mapping\n");
    GCRecovery::_reset_collation_outputs($bucket, $scratch, $logs, 'COG');
    ok(!-e "$scratch/compl.fna.gz.lock", 'stale append lock removed');
    ok(!-e "$scratch/COG/preclus.marker.fna.0", 'partial marker batch removed');
    ok(!-e "$logs/GeneCompleteness.txt.0" && !-e "$logs/Missed_samples.txt", 'stale batch reports removed');
    is(read_file("$logs/inmap.txt"), "keep mapping\n", 'mapping outside batch outputs survives');
    for my $batch (0, 1) {
        my $data = $batch ? ">s__other\nTTTT\n" : ">s__gene\nACGT\n";
        gzip(\$data => $part) or die $GzipError;
        GCRecovery::_append_file_locked($part, "$bucket/compl.fna.gz", "$scratch/compl.fna.gz.lock");
    }
    is(unzip("$bucket/compl.fna.gz"), ">s__gene\nACGT\n>s__other\nTTTT\n",
        'replaying all batches contains each gene exactly once');
};

subtest 'external genes and repeat cluster conversion' => sub {
    my $fasta = "$tmp/external.fna";
    my $rows = "$tmp/genes2rows.txt";
    write_file($fasta, ">s__gene description\nATG\nAAA\n>reference_gene external description\nCCC\nGGG\n");
    write_file($rows, "#Gene\tCluster\tRepresentative\n1\t>Cluster 0\ts__gene\n2\t>Cluster 100\treference_gene\n");
    write_file("$fasta.clstr", ">Cluster 0\n0 6nt, >s__gene... *\n>Cluster 100\n0 6nt, >reference_gene... *\n");
    my ($index) = GCRecovery::readGeneIdxSpl($rows, '__');
    GCRecovery::rewriteFastaHdIdx($fasta, $index, '__');
    is(read_file($fasta), ">1\nATGAAA\n>2\nCCCGGG\n", 'both source and external representatives receive numeric IDs');
    GCRecovery::rewriteFastaHdIdx($fasta, $index, '__');
    is(read_file($fasta), ">1\nATGAAA\n>2\nCCCGGG\n", 'numbered nucleotide catalogue survives retry unchanged');
    my $external = "$tmp/external.faa";
    my $protein = "$tmp/selected.faa";
    write_file($external, ">unused\nZZZ\n>reference_gene external description\nPG\nG\n");
    GCRecovery::attachProteins3('', $protein, $external, $index->{xtraSmpls}, '', {require_all => 1});
    is(read_file($protein), ">2\nPGG\n", 'existing streaming helper extracts and renames external proteins');
    eval { GCRecovery::attachProteins3('', "$tmp/missing.faa", $external, {absent => [3]}, '', {require_all => 1}) };
    like($@, qr/Cannot find required protein absent/, 'missing selected external protein fails explicitly');
    GCRecovery::combineClstr($fasta, $rows);
    my $clusters = read_file("$fasta.clstr");
    my $members = read_file("$fasta.clstr.idx");
    GCRecovery::combineClstr($fasta, $rows);
    is(read_file("$fasta.clstr"), $clusters, 'second conversion preserves numeric cluster headers');
    is(read_file("$fasta.clstr.idx"), $members, 'second conversion preserves cluster membership');
    unlink "$fasta.clstr.idx" or die $!;
    GCRecovery::combineClstr($fasta, $rows);
    is(read_file("$fasta.clstr.idx"), $members, 'missing member index is regenerated after interrupted conversion');
    write_file($rows, "#Gene\tCluster\tRepresentative\n9\t>Cluster 99\ts__gene\n");
    eval { GCRecovery::combineClstr($fasta, $rows) };
    like($@, qr/Cannot match cluster/, 'unrelated index is rejected');
    is(read_file("$fasta.clstr"), $clusters, 'failed conversion does not replace the original cluster file');
};

for my $local_mode (0, 1) {
subtest 'single-step backups and retries, local submission=' . $local_mode => sub {
    local $GCRecovery::tmpDir = "$tmp/cluster-work-$local_mode"; # intentionally no trailing slash
    local $GCRecovery::GCdir = "$tmp/catalogue-$local_mode";
    local $GCRecovery::QSBoptHR = {tmpSpace => '10G', constraint => []};
    local $GCRecovery::numCor0 = 1;
    local $GCRecovery::totMem = 1;
    local $GCRecovery::submitLocal = $local_mode;
    local $GCRecovery::qsubDir = "$tmp/";
    local $GCRecovery::useGTDBmg = 0;
    local $GCRecovery::clustMMseq = 1;
    local $GCRecovery::mapF = "$tmp/map.txt";
    local $GCRecovery::GCscr = 'env GENECAT_CHILD_WRAPPER=1 ' . GCRecovery::_shell_quote($^X)
        . ' ' . GCRecovery::_shell_quote(__FILE__) . ' --merge';
    local $GCRecovery::checkpointWriter = 'env GENECAT_CHECKPOINT_WRAPPER=1 '
        . GCRecovery::_shell_quote($^X) . ' ' . GCRecovery::_shell_quote("$root/helpers/writeCheckpoint.pl");
    local $ENV{PERL5LIB} = join(':', $root, $ENV{PERL5LIB} // ());
    my $scratch = $GCRecovery::tmpDir;
    my $catalogue = $GCRecovery::GCdir;
    my $bucket = "$catalogue/B0";
    make_path("$scratch/COG", "$catalogue/LOGandSUB", $bucket);
    my $regular = "$tmp/regular.fna";
    write_file($regular, ">s__gene\nATG\n");
    write_file("$regular.clstr", ">Cluster 0\n0 3nt, >s__gene... *\n");
    for my $marker (1..40) {
        write_file("$scratch/COG/COG$marker.97.fna", ">s__marker$marker\nCCC\n");
        write_file("$scratch/COG/COG$marker.97.fna.clstr", ">Cluster 0\n0 3nt, >s__marker$marker... *\n");
    }
    for my $kind (qw(compl incompl 5Pcompl 3Pcompl)) {
        my $data = $kind eq 'compl' ? read_file($regular) : '';
        gzip(\$data => "$bucket/$kind.fna.gz") or die $GzipError;
    }
    my ($complete, $merged, $clean, $cog) = map { "$catalogue/$_.stone" } qw(complete merged clean cog);
    GCRecovery::write_checkpoint($cog, parameters => {cluster_id => 97});
    no warnings 'redefine';
    local *GCRecovery::qsubSystem = sub {
        GCRecovery::systemW($_[1]);
        return ('finished-fixture-job', '');
    };
    local *GCRecovery::qsubSystemJobAlive = sub {};
    local *GCRecovery::clusterFNA = sub {
        return 'cp ' . GCRecovery::_shell_quote($regular) . ' ' . GCRecovery::_shell_quote($_[1]) . "\n"
            . 'cp ' . GCRecovery::_shell_quote("$regular.clstr") . ' ' . GCRecovery::_shell_quote("$_[1].clstr") . "\n";
    };
    my $build = sub { GCRecovery::clusterSingleStep($complete, $merged, $clean, $cog, $bucket, $catalogue, '', '') };
    is(GCRecovery::systemW($build->(), 0) // 0, 0, 'deferred clustering, actual marker merge and backups run successfully');
    my $name = $GCRecovery::primaryClusterFNA;
    my $merged_fasta = unzip("$bucket/$name.gz");
    is(scalar(() = $merged_fasta =~ /^>/mg), 41, 'published backup contains regular representative and all 40 markers');
    is(unzip("$bucket/unmerged.$name.gz"), read_file($regular), 'pre-merge backup is stored separately');
    ok(GCRecovery::_merged_catalog_backup_valid($merged, 97), 'merged checkpoint validates marker-inclusive backups');
    is($build->(), '', 'valid merged backup avoids repeating clustering or marker merge');

    unlink "$scratch/$name", "$scratch/$name.clstr";
    # Exercise the restoration branch from geneCatFlow, not a test-only restore implementation.
    my ($restore) = $source =~ /(#restore files for publication\/post-processing\n.*?)(?=\n\t\})/s;
    die 'Missing restoration branch' unless defined $restore;
    my $cmd = '';
    {
        my $pigzBin = $GCRecovery::pigzBin;
        my $bdir = $bucket;
        my $tmpDir = $scratch;
        my $primaryClusterFNA = $name;
        my $primaryClusterCLS = "$name.clstr";
        my $submitLocal = 0;
        eval $restore;
        die $@ if $@;
    }
    is(GCRecovery::systemW($cmd, 0), 0, 'restoration commands run after temporary core files are lost');
    is(read_file("$scratch/$name"), $merged_fasta, 'restored catalogue retains every marker');

    # An old completion stone has no evidence that backups include marker genes.
    GCRecovery::write_checkpoint($merged, parameters => {cluster_id => 97, stage => 'incomplete-clustering'});
    ok(!GCRecovery::_merged_catalog_backup_valid($merged, 97), 'old pre-merge completion is not trusted as a merged backup');
    is(GCRecovery::systemW($build->(), 0) // 0, 0, 'retry reconstructs from the separate pre-merge backup');
    is(unzip("$bucket/$name.gz"), $merged_fasta, 'retry neither loses nor duplicates markers');

    copy("$scratch/$name", "$catalogue/$name") or die $!;
    copy("$scratch/$name.clstr", "$catalogue/$name.clstr") or die $!;
    my $moveStone = "$catalogue/move.stone";
    GCRecovery::write_checkpoint($moveStone, parameters => {cluster_id => 97},
        outputs => ["$catalogue/$name", "$catalogue/$name.clstr"]);
    remove_tree($bucket);
    my ($resume) = $source =~ /(my \$published_catalog = .*?;)/s;
    my ($OutD, $primaryClusterFNA, $primaryClusterCLS, $protStone, $cdhID)
        = ($catalogue, $name, "$name.clstr", "$catalogue/prot.stone", 97);
    my $published = eval("package GCRecovery; $resume\n\$published_catalog");
    die $@ if $@;
    ok($published, 'published core outputs remain resumable after temporary backups are cleaned');
};
}

subtest 'Kraken accepts no classifications but propagates program failure' => sub {
    my $fake = "$tmp/kraken.pl";
    write_file($fake, 'print "U\tgene\t0\t100\t0:100\n"; exit($ENV{FAIL_KRAKEN} ? 7 : 0);');
    no warnings 'redefine';
    my $real_get = \&GCRecovery::getProgPaths;
    local *GCRecovery::getProgPaths = sub {
        return 'env GENECAT_KRAKEN_WRAPPER=1 ' . GCRecovery::_shell_quote($^X) . ' ' . GCRecovery::_shell_quote($fake) if $_[0] eq 'kraken2';
        return "$tmp/db" if $_[0] eq 'Kraken2_path_DB';
        return 'mini' if $_[0] eq 'Kraken2_mini';
        return 'false' if $_[0] eq 'taxid2tax_scr'; # must never be called with no tax IDs
        return $real_get->(@_);
    };
    local $GCRecovery::rareBin = 'false'; # no aggregation for zero classified genes
    local $GCRecovery::countMatrixF = 'Matrix.mat';
    local $GCRecovery::rtkFunDelims = '';
    local $ENV{FAIL_KRAKEN} = 0;
    my @warnings;
    local $SIG{__WARN__} = sub {push @warnings, @_};
    eval { GCRecovery::krakenTax("$tmp/kraken-empty", "$tmp/kraken-work", 1) };
    is($@, '', 'successful unclassified run does not fail under pipefail');
    ok(-e "$tmp/kraken-empty/Anno/Tax/krak2.out" && -z "$tmp/kraken-empty/Anno/Tax/krak2.out", 'empty successful output is published');
    like(join('', @warnings), qr/classified no catalog genes/, 'empty result is reported clearly');
    $ENV{FAIL_KRAKEN} = 1;
    eval { GCRecovery::krakenTax("$tmp/kraken-failed", "$tmp/kraken-work", 1) };
    ok(length($@), 'failed Kraken process still fails the pipeline');
    ok(!-e "$tmp/kraken-failed/Anno/Tax/krak2.out", 'partial failed output is not published as a completed run');
};

subtest 'program registrations and configured command wrappers' => sub {
    my %site_settings = map { $_ => 1 } qw(avx2_constraint globalTmpDir nodeTmpDir);
    my $active = join("\n", grep { !/^\s*#/ } split /\n/, $source);
    my %keys = map { $_ => 1 } $active =~ /getProgPaths\(["']([^"']+)["']/g;
    for my $key (sort keys %keys) {
        next if $site_settings{$key};
        my $value = eval { GCRecovery::getProgPaths($key, $key eq 'Kaiju_path_DB' ? 0 : 1) };
        is($@, '', "configuration registers $key");
        ok(length($value), "$key has a default command or data path") unless $key eq 'Kaiju_path_DB';
    }
    for my $key (qw(geneCat_scr writeCheckpoint_scr extre100_scr calcGC_scr genelength_scr
            hmmBestHit_scr decluterGC_scr kmerPerGene_scr)) {
        my $command = GCRecovery::getProgPaths($key);
        like($command, qr/^(?:perl|python2) /, "$key explicitly selects its interpreter");
        my ($file) = $command =~ /\s(\S+\.(?:pl|py))$/;
        ok(defined($file) && -f $file, "$key points to an existing script");
    }
    local $GCRecovery::cdhitBin = 'env GENECAT_CDHIT_WRAPPER=1 cd-hit-est';
    my $cluster_command = GCRecovery::clusterFNA('input.fna', 'output.fna', 0, 0, 97, 2, 0, "$tmp/cdhit", 0, 1);
    like($cluster_command, qr/^env GENECAT_CDHIT_WRAPPER=1 cd-hit-est -i /,
        'CD-HIT nucleotide command uses the configured executable without suffix manipulation');

    my $catalogue = "$tmp/emapper";
    make_path($catalogue);
    write_file("$catalogue/compl.incompl.97.prot.faa", ">1\nMKK\n");
    my $fake = "$tmp/configured-emapper.pl";
    write_file($fake, 'die "wrapper lost" unless $ENV{GENECAT_EMAPPER_WRAPPER}; '
        . 'open my $out, ">", $ENV{GENECAT_EMAPPER_LOG} or die $!; print {$out} join("\n", @ARGV); close $out;');
    local $ENV{GENECAT_EMAPPER_LOG} = "$tmp/emapper-args.txt";
    my $configured = 'env GENECAT_EMAPPER_WRAPPER=1 ' . GCRecovery::_shell_quote($^X)
        . ' ' . GCRecovery::_shell_quote($fake);
    no warnings 'redefine';
    my $real_get = \&GCRecovery::getProgPaths;
    local *GCRecovery::getProgPaths = sub { return $configured if $_[0] eq 'emapper'; return $real_get->(@_) };
    local *GCRecovery::splitFastas = sub { return [$_[0]] };
    my @jobs;
    local *GCRecovery::qsubSystem = sub { push @jobs, [@_]; return ('job-' . scalar(@jobs), '') };
    local $GCRecovery::GLBtmp = "$tmp/emapper-split";
    local $GCRecovery::qsubDir = "$tmp/emapper-log";
    local $GCRecovery::QSBoptHR = {tmpSpace => '0', constraint => []};
    local $GCRecovery::avx2Constr = '';
    local $GCRecovery::headBin = 'head';
    local $GCRecovery::tailBin = 'tail';
    local $GCRecovery::rareBin = 'rare';
    local $GCRecovery::rtkFunDelims = '';
    local $GCRecovery::countMatrixF = 'Matrix.mat';
    local $GCRecovery::checkpointWriter = GCRecovery::getProgPaths('writeCheckpoint_scr');
    GCRecovery::geneCatFunc_emapper($catalogue, "$tmp/emapper-tmp", 3, 0, 1, "$catalogue/emap.stone");
    my ($worker) = grep { $_->[4] eq 'eMAP0' } @jobs;
    ok(defined($worker), 'eggNOG worker is submitted');
    is(GCRecovery::systemW($worker->[1], 0), 0, 'eggNOG worker executes the configured wrapper');
    like(read_file($ENV{GENECAT_EMAPPER_LOG}), qr/--cpu\n3\n-i\n\Q$catalogue\E\/compl\.incompl\.97\.prot\.faa/,
        'configured eggNOG command receives requested cores and catalogue identity');
};

subtest 'marker extraction reads the selected catalogue identity' => sub {
    my $marker_source = read_file("$root/secScripts/GC/extrAllE100GC.pl");
    my ($extract) = $marker_source =~ /^(sub getGeneSeqsSubGenes \{.*?^\})/ms;
    my $GCd = "$tmp/markers-97";
    my $clusterID = 97;
    make_path($GCd);
    write_file("$GCd/FMG.subset.cats", "COG1\t1\t1\n");
    for my $file (['fna', 'ATG'], ['prot.faa', 'MKK']) {
        my $path = "$GCd/compl.incompl.97.$file->[0]";
        write_file($path, ">1\n$file->[1]\n");
        # Valid prebuilt .fai keeps this small fixture independent of samtools.
        write_file("$path.fai", "1\t3\t3\t3\t4\n");
    }
    my $extractor = eval("package GCRecovery; $extract\n\\&getGeneSeqsSubGenes");
    die $@ if $@;
    $extractor->('FMG');
    is(read_file("$GCd/FMG/COG1.fna"), ">1\nATG\n", 'marker nucleotides come from identity 97');
    is(read_file("$GCd/FMG/COG1.faa"), ">1\nMKK\n", 'marker proteins come from identity 97');
    };

done_testing();
