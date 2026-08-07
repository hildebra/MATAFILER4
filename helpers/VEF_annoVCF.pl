#!/usr/bin/perl
use strict;
use Mods::IO_Tamoc_progs qw(getProgPaths);
use warnings;

my $gff = getProgPaths("legacyVEF_gff");
my $fna = getProgPaths("legacyVEF_fna");
my $gbk = getProgPaths("legacyVEF_gbk");

#my $inVCF = getProgPaths("legacyVEF_inputVCF");
my $inVCF = getProgPaths("legacyVEF_inputVCF");
my $outVCF = getProgPaths("legacyVEF_outputVCF");

my $outDir = getProgPaths("legacyVEF_outputDir");
system "mkdir -p $outDir";
my $gffgz = $gff.".gz";

#VEF
my $VEFdir = getProgPaths("variantEffectPredictorDir");
my $VEFbin = "$VEFdir./variant_effect_predictor.pl";
my $CACHbin = "$VEFdir./gtf2vep.pl";
my $CAHCdir = "$VEFdir/custCache/";
#system "mkdir -p $CAHCdir";
#die "sort -k1,1 -k2,2n -k3,3n $gff | /g/bork3/home/hildebra/bin/samtools-1.2/tabix-0.2.6/./bgzip > $gffgz";
#die "/g/bork3/home/hildebra/bin/samtools-1.2/tabix-0.2.6/tabix -p gff $gffgz";
#build custom cache
#die "$CACHbin -i $gff -f $fna -d 81 -s TEC2 --dir $CAHCdir";
#die "$VEFbin --custom $gff,myFeatures,gff,overlap,0 -i $inVCF -o $outVCF --force_overwrite --offline --species TEC2 --everything --format vcf --vcf --dir $CAHCdir";


#coovar: problems with creating indices
# die "/g/bork3/home/hildebra/bin/coovar-0.07/./coovar.pl -e $gff -r $fna -v $inVCF -o $outDir\n";
 
 #snpeff
 my $SEFFbin = getProgPaths("snpEffJar");
 my $SEFFcfg = getProgPaths("snpEffConfig");
 my $SEFFdata = getProgPaths("snpEffData");
 my $geno = "TEC2.1"; my $Gname = "TEC2";
 my $custDir = "$SEFFdata/$geno/";
 if (0){ #build DB anew..
	my $cfgAdd = "# $Gname, version 1, test
	$geno.genome : $Gname ";
	if (`cat $SEFFcfg` =~ /$geno\.genome : $Gname /){
		print "Entry for such a genome already exists in DB"; exit(3);
	} else {
		#system "echo \"$cfgAdd\" >> $SEFFcfg";
	}
	system "mkdir -p $custDir; cp $fna $custDir/sequences.fa; cp $gff $custDir/genes.gff; cp $gbk $custDir/genes.gbk";
	my $cmd = "java -jar $SEFFbin build -genbank -v $geno\n"; 
 }
 my $cmd = "java -Xmx4g -jar $SEFFbin $geno $inVCF > $outVCF";
 die $cmd."\n";
 
 
 
 
 