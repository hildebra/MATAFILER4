#!/usr/bin/bash
#download example data for MG-TK
#(c) Falk Hildebrand

if [ -z ${MF4DIR+x} ]; then echo "MF4DIR is unset"; exit; else echo "MF4DIR is set to '$MF4DIR', downloading data to $MF4DIR/examples/data/"; fi
rm -rf $MF4DIR/examples/data
mkdir $MF4DIR/examples/data








echo ""
echo ""
echo "downloading PacBio example metagenomes"
echo ""

#curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR154/009/SRR15489009/SRR15489009_subreads.fastq.gz --output $MF4DIR/examples/data/SRR15489009_sub.fastq.gz
#zcat  $MF4DIR/examples/data/SRR15489009_sub.fastq.gz | head -n 100000 >  $MF4DIR/examples/data/SRR15489009_sub.fq.gz


curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR154/009/SRR15489009/SRR15489009_subreads.fastq.gz | zcat | head -n 300000 | gzip > $MF4DIR/examples/data/SRR15489009_sub.fq.gz

#curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR154/013/SRR15489013/SRR15489013_subreads.fastq.gz --output $MF4DIR/examples/data/SRR15489013_sub.fastq.gz
#zcat  $MF4DIR/examples/data/SRR15489013_sub.fastq.gz | head -n 100000 >  $MF4DIR/examples/data/SRR15489013_sub.fq.gz

curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR154/013/SRR15489013/SRR15489013_subreads.fastq.gz | zcat | head -n 300000 | gzip > $MF4DIR/examples/data/SRR15489013_sub.fq.gz


#rm -f $MF4DIR/examples/data/SRR15489013_sub.fastq.gz $MF4DIR/examples/data/SRR15489009_sub.fastq.gz




echo ""
echo "downloading illumina example metagenomes (PRJNA529586)"
echo ""
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR879/002/SRR8797712/SRR8797712_1.fastq.gz --output $MF4DIR/examples/data/SRR8797712_1.fq.gz
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR879/002/SRR8797712/SRR8797712_2.fastq.gz --output $MF4DIR/examples/data/SRR8797712_2.fq.gz
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR879/003/SRR8797713/SRR8797713_1.fastq.gz --output $MF4DIR/examples/data/SRR8797713_1.fq.gz
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR879/003/SRR8797713/SRR8797713_2.fastq.gz --output $MF4DIR/examples/data/SRR8797713_2.fq.gz
 
 
