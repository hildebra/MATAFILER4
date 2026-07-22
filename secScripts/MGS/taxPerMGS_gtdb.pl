#!/usr/bin/perl
#uses MGS to cluster genes and kraken tax assignments
#./taxPerMGS.pl /g/scb/bork/hildebra/SNP/GCs/DramaGCv5//Binning/MetaBat//MB2.clusters.ext.can.Rhcl /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/
use warnings; use strict;

use Mods::Binning qw(readMGSrev );
use Mods::IO_Tamoc_progs qw(getProgPaths );
use Mods::GenoMetaAss qw( systemW gzipopen);
use Mods::GTDBTaxonomy qw(merge_gtdb_summaries);
use File::Path qw(make_path remove_tree);
use File::Temp qw(tempdir);


#my $COND = getProgPaths("CONDA"); #source conda..
#my $py3activate = getProgPaths("py3activate"); #source conda.. 
my $GTDBtkBin = getProgPaths("GTDBtk");
my $GTDBtkDB = getProgPaths("GTDBtk_DB");
my $GTDBtkMash = getProgPaths("GTDBtk_mash",0);

#die "$GTDBtkDB\n$GTDBtkMash\n";
die "Usage: $0 genome-dir cores temporary-root output-dir\n" unless @ARGV == 4;
my $refMGd = $ARGV[0];
my $ncore = $ARGV[1];
my $tmpD = $ARGV[2];
my $Bdir = $ARGV[3];
die "Genome directory does not exist: $refMGd\n" unless -d $refMGd;
die "Core count must be positive\n" unless defined($ncore) && $ncore =~ /^\d+$/ && $ncore > 0;
die "Temporary root and output directory are required\n" unless length($tmpD) && length($Bdir);

my $pplacer_cores = $ncore;
$pplacer_cores = 2 if ($ncore > 2);

$tmpD.="/GTDB/";
make_path($tmpD, $Bdir);
# Each attempt receives a private output directory.  Reusing GTDBTK/ allowed a
# failed rerun to discover summaries left by an earlier attempt (and made two
# concurrent jobs share intermediate state).
my $runD = tempdir("run-XXXXXX", DIR => $tmpD, CLEANUP => 0);
my $oDir = "$runD/GTDBTK/";
make_path($oDir);

my $summaryOut = "$Bdir/gtdbtk.summary.tsv";
my $taxonomyOut = "$Bdir/GTDBTK.tax";
# Do not let the surrounding shell's output tests accept results from a failed
# earlier attempt if this invocation exits before publishing new files.
for my $oldOutput ($summaryOut, $taxonomyOut) {
	next unless -e $oldOutput;
	unlink $oldOutput or die "Cannot remove stale GTDB taxonomy output $oldOutput: $!\n";
}



if ($GTDBtkBin =~ m/ activate /){
	$GTDBtkBin =~ s/activate (\S+)/activate $1\nexport GTDBTK_DATA_PATH=$GTDBtkDB\n/;
}


#get GTDBtk version
my $verSt = `$GTDBtkBin --version`;
die "Could not run GTDB-Tk version command\n" if $? != 0;
my ($version_text) = $verSt =~ /(\d+(?:\.\d+)+)/;
die "Could not parse GTDB-Tk version from: $verSt\n" unless defined $version_text;
my ($major, $minor) = split /\./, $version_text;
my $GTDBver = "$major.$minor" + 0;

#print "$verSt\n\n$1\n";

print "Using GTDBtk ver $GTDBver\n";
#print "$GTDBtkBin --version\n";
#die;

my $cmd = "";
#$cmd .= "$COND\n$py3activate\n";
#--scratch_dir $tmpD 
#$cmd .= "export GTDBTK_DATA_PATH=$GTDBtkDB\n";

my $mashArg="" ;my $hook = "";
if ($GTDBtkMash ne "" && $GTDBver >= 2.1){ #for newer GTDBtk versions not supported
	$mashArg = "--mash_db $GTDBtkMash/\$MVERSION/";
	$hook = "MVERSION=\$(mash --version | grep -oE '[0-9]+(\\.[0-9]+)+' | head -n1)\n";
	$hook .= "test -n \"\$MVERSION\"\n";
}

$cmd .= $hook;
$cmd .= "$GTDBtkBin classify_wf -x fna $mashArg --cpus $ncore --pplacer_cpus $pplacer_cores --genome_dir $refMGd --out_dir $oDir"; #--scratch_dir $tmpD/GTtmp/ --genes

print "\n\n".$cmd."\n\n";
#die;
systemW $cmd;

my @summaries = sort glob("$oDir/gtdbtk.*.summary.tsv");
die "GTDB-Tk produced no summary files in $oDir\n" unless @summaries;
my $taxonomyRows = merge_gtdb_summaries(\@summaries, $summaryOut, $taxonomyOut);
print "Published $taxonomyRows GTDB taxonomy assignments\n";
systemW "tar -zcvf $refMGd/GTDBtk.tar.gz $oDir";
remove_tree($runD) if -d $runD;

#transfer files
