use strict;
use warnings;
use Test::More;
use File::Spec;
use File::Temp qw(tempdir);
use IO::Compress::Gzip qw(gzip $GzipError);

my $root = File::Spec->rel2abs('.');
my $tmp = tempdir(CLEANUP => 1);

sub write_file {
    my ($path, $contents) = @_;
    open my $fh, '>', $path or die "Cannot open $path: $!";
    print {$fh} $contents or die "Cannot write $path: $!";
    close $fh or die "Cannot close $path: $!";
}

sub read_file {
    my ($path) = @_;
    open my $fh, '<', $path or die "Cannot open $path: $!";
    local $/;
    my $contents = <$fh>;
    close $fh or die "Cannot close $path: $!";
    return $contents;
}

my $rename_input = File::Spec->catfile($tmp, 'rename.fasta');
write_file($rename_input, ">first\nACGT\n>second\nAACCGG\n");
is(system($^X, File::Spec->catfile($root, 'secScripts', 'assemblies', 'renameCtgs.pl'), $rename_input, 'sample'), 0,
   'renameCtgs completes');
my $renamed = read_file($rename_input);
like($renamed, qr/^>sample__C1_L=4=\nACGT\n>sample__C2_L=6=\nAACCGG\n$/,
     'renameCtgs records the correct length for the final contig');

my $filter_input = File::Spec->catfile($tmp, 'filter.fasta');
write_file($filter_input, ">long\nAAAAAAAAAA\n>last-secondary\nCCCCCC\n");
is(system($^X, File::Spec->catfile($root, 'secScripts', 'assemblies', 'sizeFilterFas.pl'), $filter_input, 8, 5), 0,
   'sizeFilterFas completes');
like(read_file("$filter_input.filt"), qr/>long\nAAAAAAAAAA\n/, 'primary output contains long sequence');
like(read_file("$filter_input.filt2"), qr/>last-secondary\nCCCCCC\n/,
     'secondary output contains the final short FASTA record');

my $gene_cat = read_file(File::Spec->catfile($root, 'secScripts', 'geneCat.pl'));
unlike($gene_cat, qr/rm -rf \$GCdir\/\* \$tmpDir\*/, 'geneCat has no wildcard clean-start deletion');
unlike($gene_cat, qr/system "rm -r \$metaGD\/\$path2GPdir/, 'geneCat does not delete predictions while inspecting them');
unlike($gene_cat, qr/system "rm -rf \$metaGD\/\$path2CS/, 'geneCat does not delete contig stats while inspecting them');
like($gene_cat, qr/genemat\.done\.sh/, 'matrix completion uses a convergence job');
like($gene_cat, qr/No usable assembly.*if \$requireAllAssemblies/s,
     'missing assemblies fail only when requireAllAssemblies is enabled');
like($gene_cat, qr/\$map\{\$smpl\}\{assFinSmpl\} eq \$smpl/,
     'assembly-group precheck recognizes the explicitly final assembly sample');
unlike($gene_cat, qr/!\s*fileGZe\("\$metaGD\/scaffolds\.fasta\.filt"\) \|\| !-e "\$metaGD\/longReads/,
       'hybrid assembly precheck does not require both short- and long-read assemblies');
unlike($gene_cat, qr/my \$cmd \.= "\$kaijBin/, 'Kaiju command is initialized before concatenation');

my $parse = read_file(File::Spec->catfile($root, 'secScripts', 'functions', 'parseBlastFunct2.pl'));
like($parse, qr/CNT_\$\{minBLE\}_\$\{minPID\}/, 'functional result checks use threshold and percent identity');
like($parse, qr/\.\$normMethod\.gene\.cnts\.gz/, 'functional result checks include normalization in output names');

my $abr_db = File::Spec->catdir($tmp, 'abr-db');
mkdir $abr_db or die "Cannot create $abr_db: $!";
write_file(File::Spec->catfile($abr_db, 'ardb.tabs.parsed'),
           "SUB\tx\tSYM\tCAT\tx\tx\tx\t80\n");
write_file(File::Spec->catfile($abr_db, 'ardb_and_reforg_mapping'),
           "x\tSYM\tx\tDRUG\n");
write_file(File::Spec->catfile($abr_db, 'ardb_vs_reforg9f.overlap90shortest_famthres_or_symbol.sorted.besthit'),
           "ALT\tSUB\n");
my $abr_blast = File::Spec->catfile($tmp, 'abr.srt.gz');
my $blast_text = join '',
    "paired/1\tSUB\t90\t50\t0\t0\t1\t50\t1\t50\t1e-20\t100\n",
    "paired/2\tSUB\t92\t50\t0\t0\t1\t50\t51\t100\t1e-20\t105\n",
    "second-only/2\tSUB\t95\t80\t0\t0\t1\t80\t1\t80\t1e-20\t120\n";
gzip(\$blast_text => $abr_blast) or die "Cannot create $abr_blast: $GzipError";
my $abr_genes = File::Spec->catfile($tmp, 'abr.genes.txt');
my $abr_cats = File::Spec->catfile($tmp, 'abr.cats.txt');
is(system($^X, '-I' . $root,
          File::Spec->catfile($root, 'secScripts', 'functions', 'ABRblastFilter2.pl'),
          $abr_blast, $abr_genes, $abr_cats, $abr_db), 0,
   'ABR filter handles paired and read-2-only hits');
my $abr_output = read_file($abr_cats);
like($abr_output, qr/^paired\/1\tSUB\t/m, 'ABR filter combines a paired hit without losing its subject key');
like($abr_output, qr/^second-only\/2\tSUB\t/m, 'ABR filter retains read-2-only hits');
ok(-e "$abr_blast.stone", 'ABR completion marker is written after successful output');

my $mgs = read_file(File::Spec->catfile($root, 'secScripts', 'MGS.pl'));
like($mgs, qr/Select exactly one quality checker/, 'MGS rejects ambiguous CheckM/CheckM2 configuration');
like($mgs, qr/runCheckM\(\$binCanDir,\$ChkMevalF/, 'MGS supports CheckM1 for canopy quality checks');
like($mgs, qr/my \$sco = \$spl\[12\]-\(\$spl\[13\]\*2\)/,
     'MAG replacement score uses completeness and contamination columns');
unlike($mgs, qr/my \$testKey = ">\$\{cc\}_\$cnt";\s*my \$curGene = ""; my \$cnt=0/s,
       'contig gene lookup does not use the outer counter before initializing its own counter');
like($mgs, qr/if \(\$LOGstr =~ m\/:::Correct:\/\)/, 'Rhcl success output is recognized without a stray quote');
unlike($mgs, qr/m\/\\":::Correct:/, 'Rhcl parser no longer requires an impossible leading quote');
like($mgs, qr/test -s \$GTDBtaxF.*touch \$GTDBtaxSto/s,
     'GTDB checkpoint follows validation of final taxonomy outputs');
like($mgs, qr/test -s \$annoDir\/specI\.tax\\n";\s*\$cmdSI \.= "touch \$ABmgsSton/s,
     'MGS abundance checkpoint follows final output validation');
like($mgs, qr/_touch_checkpoint\(\$iniMB2sto\) unless -e \$iniMB2sto \|\| \@missedMAGs/,
     'missing MAG groups prevent the global MAG checkpoint from becoming sticky');
unlike($mgs, qr/foreach my \$Doo \(\@DoosD\)\{\s*last if \(-e "\$iniMB2sto"\)/s,
       'MGS validates MAG outputs even when a previous global checkpoint exists');

done_testing();
