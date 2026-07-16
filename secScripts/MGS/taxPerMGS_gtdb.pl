#!/usr/bin/perl
#uses MGS to cluster genes and kraken tax assignments
#./taxPerMGS.pl /g/scb/bork/hildebra/SNP/GCs/DramaGCv5//Binning/MetaBat//MB2.clusters.ext.can.Rhcl /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/
use warnings; use strict;

use Mods::Binning qw(readMGSrev );
use Mods::IO_Tamoc_progs qw(getProgPaths );
use Mods::GenoMetaAss qw( systemW gzipopen);
use File::Path qw(make_path remove_tree);


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
my $oDir = "$tmpD/GTDBTK/";
make_path($oDir);



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
open my $summary_out, '>', "$Bdir/gtdbtk.summary.tsv" or die "Cannot write GTDB summary: $!\n";
for my $summary (@summaries) {
	open my $summary_in, '<', $summary or die "Cannot read $summary: $!\n";
	while (my $line = <$summary_in>) { print {$summary_out} $line; }
	close $summary_in or die "Cannot close $summary: $!\n";
}
close $summary_out or die "Cannot close $Bdir/gtdbtk.summary.tsv: $!\n";
systemW "cut -f1,2 $Bdir/gtdbtk.summary.tsv | sed 's/.__//g' > $Bdir/GTDBTK.tax";
systemW "tar -zcvf $refMGd/GTDBtk.tar.gz $oDir";
remove_tree($tmpD) if -d $tmpD;

#transfer files
