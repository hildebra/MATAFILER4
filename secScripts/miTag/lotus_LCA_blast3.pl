#!/usr/bin/env perl
#assigns LCA tax to LSU / SSU / ITS seq fragments (extracted from metag)
#ex ./lotus_LCA_blast.pl /g/scb/bork/hildebra/Tamoc/FinSoil/Sample_G2677_115/ribos FS31 20
#  ./lotus_LCA_blast.pl /g/scb/bork/hildebra/SNP/SimuG/simulated_FungiSA_5/ribos S5 20 /tmp/hildebra/simu/S5/
use strict;
use warnings;
#use threads;
use Mods::GenoMetaAss qw(reverse_complement_IUPAC readFasta systemW writeFasta);
use Mods::IO_Tamoc_progs qw(getProgPaths setConfigFile);
use Getopt::Long qw( GetOptions );
use File::Basename qw(basename);
use File::Copy qw(move);
use File::Path qw(make_path remove_tree);
use File::Spec;



sub getTaxForOTUfromRefBlast;
sub splitBlastTax;
sub writeBlastHiera; sub runBlastLCA;
sub merge;
sub fastq2fna; sub flashed;
sub rework_tmpLines;
sub extractReads;
sub materializeReadInput;
sub resultExists;
sub touchFile;
sub shellQuote;
sub countLines;



#some more pars for LCA	
my @idThr = (97,95,93,91,88,78);
my $lengthTolerance = 0.85;
my $maxRds = 50000;
my $maxHitOnly = 0;
my $DoPar = 0;
my $blMode = 4; #2=lambda,1=Blast,3=sortmeRNA, 4=vsearch


if (@ARGV ==0 ){die"no Input arguments given\n";}
my $inD = ""; #$ARGV[0]; #expects TAMOC output files in this dir
my $SmplID = "";#$ARGV[1];#Name of the Sample
my $BlastCores = 1;#$ARGV[2]; #parrallel execution
my $tmpD = ""; #tmp dir for faster I/O
my $readsRpairs=0;
my $cfgFile="";
my $extractDNA = 0;
my $DBdir = "";

GetOptions(
	"dir=s" => \$inD,
	"smplID=s"      => \$SmplID,
	"tmpD=s"      => \$tmpD,
	"pairedRds=i"      => \$readsRpairs, #has to be 1 or 2
	"cores=i" => \$BlastCores,
	"simMode=i" => \$blMode,
	"config=s" => \$cfgFile, #main MATAFILER config file
	"keepReads=i"      => \$extractDNA,
	"maxReadNum=i"      => \$maxRds, #how many reads to actually use in LCA
	"lengthTolerance=f"	=> \$lengthTolerance, #what fraction needs to align
	"DBdir=s" => \$DBdir,
) or die("Error in command line arguments\n");
die "Unexpected positional arguments: @ARGV\n" if @ARGV;

die "Needs at least -dir, -smplID and -DBdir arguments\n" if ($DBdir eq "" || $inD eq "" || $SmplID eq "");
die "Input directory does not exist: $inD\n" unless -d $inD;
die "Database directory does not exist: $DBdir\n" unless -d $DBdir;
die "-cores must be a positive integer\n" if $BlastCores < 1;
die "-pairedRds must be 0, 1, or 2\n" unless $readsRpairs == 0 || $readsRpairs == 1 || $readsRpairs == 2;
die "-simMode must be 2 (LAMBDA), 3 (SortMeRNA), or 4 (VSEARCH)\n"
	unless $blMode == 2 || $blMode == 3 || $blMode == 4;
die "-maxReadNum must be zero or greater\n" if $maxRds < 0;
die "-lengthTolerance must be between 0 and 1\n" if $lengthTolerance < 0 || $lengthTolerance > 1;

$inD .= "/" unless ($inD =~ m/\/$/);
#dir to store LCA
my $outdir = $inD."ltsLCA/";
my $tmpRoot = $tmpD eq "" ? "$outdir/tmp/" : $tmpD;
$tmpRoot = File::Spec->canonpath(File::Spec->rel2abs($tmpRoot));
make_path($outdir) unless -d $outdir;
make_path($tmpRoot) unless -d $tmpRoot;
my $safeSample = $SmplID;
$safeSample =~ s/[^A-Za-z0-9_.-]+/_/g;
$tmpD = File::Spec->catdir($tmpRoot, "lotusLCA_${safeSample}_$$");
make_path($tmpD);

#dir to store merged reads in
my $mergeD = $tmpD."flashMerge/";
make_path($mergeD);


#extractReads("/tmp/MATAFILER//FS2ITS//LCA/flashMerge/SSU.extendedFrags.fna","/g/scb/bork/hildebra/Tamoc/SoiSciReb/FS2/ribos/ltsLCA//SSUriboRun_bl.hiera.txt");
#die "USA\n";

setConfigFile($cfgFile);
#binaries
my $flashBin = $readsRpairs ? getProgPaths("flash") : "";
my $lambdaBin = $blMode == 2 ? getProgPaths("lambda") : "";
my $vsearchBin= $blMode == 4 ? getProgPaths("vsearch") : "";
#my $lambdaIdxBin = $lambdaBin."_indexer";#getProgPaths("");#"/g/bork3/home/hildebra/dev/lotus//bin//lambda/lambda_indexer";
my $LCAbin = getProgPaths("LCA");#"/g/bork3/home/hildebra/dev/C++/LCA/./LCA";
my $smrnaBin = $blMode == 3 ? getProgPaths("sortmerna") : "";


my @DBn = ("LSUdbFA","LSUtax","SSUdbFA","SSUtax","ITSdbFA","ITStax","PR2dbFA","PR2tax");
my $LCAar = getProgPaths(\@DBn,0);
my @LCAdbs = @{$LCAar}; 
#reduce to files
for (my $i=0;$i<@LCAdbs;$i++){
	next if ($LCAdbs[$i] eq "");
	$LCAdbs[$i] = basename($LCAdbs[$i]);
}

#datbases
#datbases
my ($LSUdbFA,$LSUtax) = ("","");
unless ($LCAdbs[0] eq ""){$LSUdbFA = "$DBdir/$LCAdbs[0]";$LSUtax = "$DBdir/$LCAdbs[1]";}
my ($SSUdbFA,$SSUtax) = ("","");
unless ($LCAdbs[2]  eq ""){$SSUdbFA = "$DBdir/$LCAdbs[2]";$SSUtax = "$DBdir/$LCAdbs[3]";}
my ($ITSdbFA,$ITStax) = ("","");
unless ($LCAdbs[4] eq ""){$ITSdbFA = "$DBdir/$LCAdbs[4]";$ITStax = "$DBdir/$LCAdbs[5]";}
my ($PR2dbFA,$PR2tax) = ("","");
unless ($LCAdbs[6] eq ""){$PR2dbFA = "$DBdir/$LCAdbs[6]";$PR2tax = "$DBdir/$LCAdbs[7]";}
#my $ITSdbFA = getProgPaths("ITSdbFA"); $ITSdbFA =~ m/([^\/]+)$/; $ITSdbFA = $1; $ITSdbFA = "$DBdir/$ITSdbFA";
#my $ITStax = getProgPaths("ITStax");$ITStax =~ m/([^\/]+)$/; $ITStax = $1;$ITStax = "$DBdir/$ITStax";
#my $PR2dbFA = "$DBdir/gb203_pr2_all_10_28_99p.fasta";
#my $PR2tax = "$DBdir/PR2_taxonomy.txt";
#die "$SSUdbFA\n";
#some cleanups (includes prev runs)
unlink $_ for (glob("$inD/*.blast"), glob("$outdir/*riboRun_bl"));

my @tags=("SSU","LSU");
my $allAlreadyComplete = -e "$outdir/Assigned.sto";
for my $tag (@tags){
	$allAlreadyComplete &&= -e "$outdir/${tag}_ass.sto"
		&& resultExists("$outdir/${tag}riboRun_bl.hiera.txt", 1);
}
if ($allAlreadyComplete){
	print "All assigned already\n";
	remove_tree($tmpD);
	exit(0);
}

my $inputOK=1;
for my $tag (@tags){
	my $checkpoint = "$outdir/${tag}_ass.sto";
	my $hierarchy = "$outdir/${tag}riboRun_bl.hiera.txt";
	next if -e $checkpoint && resultExists($hierarchy, 1);
	unlink $checkpoint if -e $checkpoint;

	my (@dbfa,@dbtax);
	if ($tag eq "LSU") {@dbfa = ($LSUdbFA); @dbtax = ($LSUtax);
	} elsif ($tag eq "SSU") {@dbfa = ($PR2dbFA,$SSUdbFA); @dbtax = ($PR2tax,$SSUtax);
	} else {die "Unrecognized tag $tag\n";}
	my (@usableDb,@usableTax);
	for my $index (0..$#dbfa){
		next if !defined($dbfa[$index]) || $dbfa[$index] eq "";
		die "Missing $tag reference database: $dbfa[$index]\n" unless -s $dbfa[$index];
		die "Missing $tag taxonomy database: $dbtax[$index]\n" unless defined($dbtax[$index]) && -s $dbtax[$index];
		push @usableDb, $dbfa[$index];
		push @usableTax, $dbtax[$index];
	}
	die "No usable databases configured for $tag\n" unless @usableDb;

	my $singleInput = materializeReadInput("$inD/reads_${tag}.fq", $tmpD);
	my $r1Input = materializeReadInput("$inD/reads_${tag}.r1.fq", $tmpD);
	my $r2Input = materializeReadInput("$inD/reads_${tag}.r2.fq", $tmpD);
	my $markerMode = $readsRpairs;
	my $go = 1;
	if ($markerMode > 0){
		$go = flashed($r1Input,$r2Input,$mergeD,$tag,$inD);
		if ($go == 1){
			runBlastLCA($mergeD.$tag,$singleInput,\@usableDb,\@usableTax,"${tag}riboRun_bl",$blMode,$SmplID,$tag,$go,$markerMode);
		} elsif ($go == 3){
			$markerMode = 0; # marker-local fallback; do not change later markers
		} else {
			warn "Problem with $tag paired input; run needs to be repeated\n";
			$inputOK=0;
			next;
		}
	}
	if ($markerMode == 0){
		$go = 0 if $singleInput eq "" || !-s $singleInput;
		runBlastLCA($singleInput,"",\@usableDb,\@usableTax,"${tag}riboRun_bl",$blMode,$SmplID,$tag,$go,$markerMode);
	}

	my $allowEmpty = $go == 0;
	unless (resultExists($hierarchy, $allowEmpty)){
		warn "Expected hierarchy output was not created for $tag: $hierarchy\n";
		$inputOK=0;
		next;
	}
	touchFile($checkpoint);
}

remove_tree($tmpD) if -d $tmpD;
if ($inputOK){
	for my $tag (@tags){
		die "Cannot mark assignment complete; $tag output is incomplete\n"
			unless -e "$outdir/${tag}_ass.sto" && resultExists("$outdir/${tag}riboRun_bl.hiera.txt", 1);
	}
	touchFile("$outdir/Assigned.sto");
	print "$outdir/Assigned.sto\n";
	exit(0);
}
die "Ribosomal LCA assignment was incomplete\n";


#flash merge & some file checks
sub flashed($ $ $ $ $){
	my ($r1,$r2,$outD,$outT,$primD) = @_;
	return 3 if ($r1 eq "" || !-s $r1) && ($r2 eq "" || !-s $r2);
	if ($r1 eq "" || $r2 eq "" || !-s $r1 || !-s $r2){
		warn "Incomplete paired input for $outT: '$r1' / '$r2'\n";
		unlink "$primD/${outT}_pull.sto" if -e "$primD/${outT}_pull.sto";
		return 2;
	}
	my $r1Lines = countLines($r1);
	my $r2Lines = countLines($r2);
	if ($r1Lines != $r2Lines || $r1Lines % 4 != 0){
		warn "Paired FASTQ files have incompatible record counts: $r1Lines / $r2Lines lines\n";
		unlink "$primD/${outT}_pull.sto" if -e "$primD/${outT}_pull.sto";
		return 2;
	}
	print "running flash..\n";
	my $mergCmd = "$flashBin -M 200 -o ".shellQuote($outT)." -d ".shellQuote($outD)." -t $BlastCores ".shellQuote($r1)." ".shellQuote($r2);
	if (systemW($mergCmd, 0)){
		print STDERR "\n$mergCmd\nfailed\n";
		unlink "$primD/${outT}_pull.sto" if -e "$primD/${outT}_pull.sto";
		return 2;
	}
	return 1;
}

sub fastq2fna($ $){
	my ($in,$doDel) = @_;
	#deactivate, since lambda can just read fq...
	#reactivate for merging of reads
	#return $in;
	#print $in."\n";
	return $in if (-z $in);
	my $out = $in;
	$out =~ s/\.f[^\.]*q$/\.fa/g;
	die "Couldn't convert $in to .fa ending\n" if ($in eq $out);
	#die $out;
	#$out =~ s/\.fq$/\.fna/g;
	open my $inputHandle,"<",$in or die "Input fastq file $in not available: $!\n";
	my $firstLine = <$inputHandle>;
	die "Input sequence file is empty: $in\n" unless defined $firstLine;
	if ($firstLine =~ /^>/){close $inputHandle; return $in;}
	seek($inputHandle, 0, 0) or die "Cannot rewind $in: $!\n";
	open my $outputHandle, ">", $out or die "Cannot write FASTA output $out: $!\n";
	my $record = 0;
	while (my $header = <$inputHandle>){
		my $sequence = <$inputHandle>;
		my $plus = <$inputHandle>;
		my $quality = <$inputHandle>;
		die "Truncated FASTQ record ".($record + 1)." in $in\n"
			unless defined($sequence) && defined($plus) && defined($quality);
		die "Invalid FASTQ header in record ".($record + 1)." of $in\n" unless $header =~ /^@/;
		die "Invalid FASTQ separator in record ".($record + 1)." of $in\n" unless $plus =~ /^\+/;
		chomp($header, $sequence);
		$header =~ s/^@/>/;
		print {$outputHandle} "$header\n$sequence\n";
		$record++;
	}
	close $inputHandle or die "Cannot close $in: $!\n";
	close $outputHandle or die "Cannot close $out: $!\n";
	unlink $in or die "Cannot remove temporary FASTQ $in: $!\n" if $doDel && $in ne $out;
	return $out;
}

#find reads that could not match any LCA and redo with different DB
sub findUnassigned($ $ $ ){
	my ($BTr,$Fr,$outF) = @_;
	my %BT = %{$BTr}; my %Fas = %{$Fr};
	#my @t =keys %Fas; print "$t[0]\n";
	my $cnt =0; my $dcn =0 ;
	my @kk = keys %Fas;
	my $totFas = @kk;
	
	if ($outF eq ""){
		foreach my $k (keys %BT){
			my @curT = @{$BT{$k}};
			if ( @curT==0 || $curT[0] eq "?" ){$dcn++;} #|| $curT[1] eq "?"
			$cnt++;
		}
		if (@kk == 0){
			print "Total of ". ($cnt-$dcn)." / $cnt reads have LCA assignments\n";
		} else {
			print "$dcn / $cnt reads failed LCA assignments, checked $totFas reads.\n";
		}
		return;
	}
	foreach my $k (keys %BT){
		my @curT = @{$BT{$k}};
		#print $k."\t${$BT{$k}}[0]   ${$BT{$k}}[2]\n" 
		if ( @curT==0 || $curT[0] eq "?" ){#|| $curT[1] eq "?"){
			delete $BT{$k};
			$dcn++;
			
			#print ">".$k."\n".$Fas{$k}."\n";
		} #else {print $k."\t${$BT{$k}}[0]   ${$BT{$k}}[2]\n" ;}
		else {
			die "Can't find fasta entry for $k\n" unless (exists $Fas{$k});
			delete $Fas{$k};
		}
		$cnt ++;
		#die if ($cnt ==100);
	}
	@kk = keys %Fas;
	print "$dcn / $cnt reads failed LCA assignments\nWriting ".@kk." of previous $totFas reads for next iteration.\n";
	open O,">$outF" or die "can;t open unassigned fasta file $outF\n";
	foreach my $k(@kk){
		print O ">".$k."\n".$Fas{$k}."\n";
	}
	close O;
	return ($outF,\%BT,$dcn,\%Fas);
}

#main routine that does sim search & starts the LCA
sub runBlastLCA{
	my ($queryO,$queryXtrSingl,$DBar,$DBtaxar,$id,$doblast,$SmplID,$MKname,$go,$pairedMode) =@_;
	my $taxblastf_base = $outdir."$id";
	if ($tmpD ne ""){
		$taxblastf_base = $tmpD."$id";
	}
	#die "$queryXtrSingl\n";
	my $doInter=1; my $doQuery=1; #fine control of what subparts to do..
	my $r1 = $queryO; my $r2 = $queryO;
	my $interLeaveO = "";
	my $hof1 = "$outdir/$id.hiera.txt"; my $hof2 = "$outdir/$id"."_inter.hiera.txt";
	#check for empty input
	if ($go == 2){ return ;}
	if (!$go){
		print "$id has empty input files\n";
		touchFile($hof1);
		return;
	}

	if (!$pairedMode){
		die "can't find single read input file: $queryO\n" unless (-e $queryO);
		#jsut check by default if this is fna
		$queryO = fastq2fna($queryO,0);
			#die "$queryO\n";
		#$queryO = fastq2fna($queryO);
		#die "$queryO\n";
	}elsif ($queryO !~ m/\.fa$/){#merged fastq that need to be processed
		my $queryOx = "$queryO.extendedFrags.fastq";
		$queryO = fastq2fna($queryOx,1);
		
		$r1 .= ".notCombined_1.fastq";
		$r2 .= ".notCombined_2.fastq";
		$interLeaveO = $tmpD."inter$id.fa";
		#print $r1." V\n";
		#$r1 = fastq2fna($r1); $r2 = fastq2fna($r2); $queryO = fastq2fna($queryO);
		
		merge($r1,$r2,$interLeaveO) if ($doInter);
		#die "$queryXtrSingl\n";
		if ($pairedMode==2 && -e $queryXtrSingl){
			print "attaching single reads to interleave";
			$queryXtrSingl = fastq2fna($queryXtrSingl,0);
			systemW "cat $queryXtrSingl >> $interLeaveO \nrm $queryXtrSingl";
			}
	} else {
		print "no interLeaveO\n";
		$doInter=0;
	}
	#die "YY\n";
	if (-z $interLeaveO){print "No interleaved files $id\n";$doInter=0;}
	if (-z $queryO){print "No merged files $id\n";$doQuery=0;}
	#my $BlastCores = 20;
#die "$queryO\n";
	#my $fasrA = readFasta( $queryO,1); 
	#my %fas=%{$fasrA}; my @t = keys %fas; die "$queryO\n@t\n$t[0]\n";
	#my $fasrI = readFasta($interLeaveO,1); 
	#print $DBar."\n";
	my @DBa = @{$DBar}; my @DBtaxa = @{$DBtaxar};
	my $query= $queryO; my $interLeave = $interLeaveO;
	my $BlastTaxRi = {}; my $BlastTaxR = {};
	my $fullBlastTaxRi = {}; my $fullBlastTaxR = {};
	my $simName = ""; my $leftover = 111;
	
	#takes too long, first check how many reads (and if this can be reduced)
	if ($interLeave ne "" && -e $interLeave){$doInter=1;} else {$doInter=0;}
	my $totRds1 = 0; #merged
	my $totRds2 = 0; #interleaved
	if ($maxRds > 0){
		my $hr = readFasta($query,0);
		$totRds1 = scalar(keys(%{$hr}));
		print "Found $totRds1 merged candidates in $query\n";
		print "Using max $maxRds of $totRds1 (total).\n" if ($totRds1 > $maxRds);
		writeFasta($hr,$query,$maxRds - 1) if $totRds1 > $maxRds;
	}

	if ($doInter && ($maxRds == 0 || $maxRds > $totRds1)){
		if ($maxRds > 0){ #prefer merged reads and cap the combined total
			my $hr = readFasta($interLeave,0);
			$totRds2 = scalar(keys(%{$hr}));
			my $remaining = $maxRds-$totRds1;
			print "..and $totRds2 interleave candidates in $query\n";
			print "Using $remaining of $totRds2 (+ $totRds1 total).\n" if $totRds2 > $remaining;
			writeFasta($hr,$interLeave,$remaining - 1) if $totRds2 > $remaining;
		}
		systemW("cat ".shellQuote($interLeave)." >> ".shellQuote($query));
		unlink $interLeave or die "Cannot remove appended interleaved reads $interLeave: $!\n";
		#$cmd .= "$lambdaBin $defLopt -q $interLeave -i $DB.lambda -o $taxblastf2\n";
		#$cmd .= "\nmv $tmptaxblastf $taxblastf2\n";
		#unless ( -e $taxblastf2){	
		print "Also comparing interleaved files\n";
	}
	my @taxouts=(); my @taxouts_inter = (); my @DBtaxa_used = ();
	for (my $DBi=0;$DBi<@DBa; $DBi++){
		my $DB = $DBa[$DBi]; my $DBtax = $DBtaxa[$DBi];
		next if ($DB eq "");
		#print "Running sim search $doblast..\n";
		my $taxblastf = $taxblastf_base.".$DBi.m8.gz";
		my $taxblastf2=$taxblastf_base.".$DBi.i.m8.gz";
		push @taxouts, $taxblastf;
		push @taxouts_inter, $taxblastf2;
		push @DBtaxa_used, $DBtaxa[$DBi];
		
		if (!-e $query){$doQuery=0;} else {$doQuery=1;}
		#die "$doQuery $doInter\n$interLeave\n";
		if ($doblast == 1){
			die "using blast is deprecated and no longer supported!\n";
			$taxblastf.=".blast"; $simName="blast";
			print "Running Blast\n";
			my $mkBldbBin = getProgPaths("makeblastdb");
			my $blastBin = getProgPaths("blastn");
			my $cmd = "$mkBldbBin -in $DB -dbtype 'nucl'\n";
			unless (-f $DB.".nhr"){	systemW($cmd);}
			my $strand = "both";
			#-perc_identity 75
			$cmd = "";
			if ($doQuery ){
				$cmd .= "$blastBin -query $query -db $DB -out $taxblastf -outfmt 6 -max_target_seqs 50 -evalue 0.1 -num_threads $BlastCores -strand $strand \n"; #-strand plus both minus
			}
			unless ( -e $taxblastf){	print "Running blast on combined files\n";
				systemW($cmd);
			} else {	print "Blast output $taxblastf does exist\n";}
			
		} elsif ($doblast==2){
			$simName = "lambda";
			if (0 && !-f "$DB.lba.gz"){ #lambda3 # !-f $DB.".dna5.fm.sa.val"  ) { #don't do this at all from nodes
				print "Building LAMBDA index anew (may take some time)..\n";
#				my $cmdIdx = "$lambdaIdxBin -p blastn -t $BlastCores -d $DB";
				my $cmdIdx = "$lambdaBin mkindexn -t $BlastCores -d $DB";
				$cmdIdx .= "gzip $DB.lba;\n";
				if (systemW($cmdIdx)){die ("Lamdba ref DB build failed\n$cmdIdx\n");}
			} elsif (!-f "$DB.lba.gz"){#!-d "$DB.lambda" || !-f "$DB.lambda/index.lf.drp"){
				die "Can not find required lambda index dir at $DB.lambda\n";
			}
			print "Starting LAMBDA similarity search..\n";
			my $tmptaxblastf = "$outdir/tax.m8";
			$tmptaxblastf = "$tmpD/tax.m8" unless ($tmpD eq "");
			my $cmd = "";
#			my $outcols = "'qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen'";
			my $outcols = "'qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen'";
			my $defLopt = " searchn -t $BlastCores --percent-identity 75 --num-matches 200 --e-value 1e-8 --output-columns $outcols ";
			#my $defLopt = "-t $BlastCores -id 75 -nm 200 -p blastn -e 1e-8 -so 7 -sl 16 -sd 1 -b 5 -pd on -oc $outcols ";
			#my $defLopt = "searchn -t $BlastCores --percent-identity 75 -n 200 -e 1e-12 --seed-offset 7 --seed-length 14 --seed-delta 1 -b -3 --filter-putative-duplicates on ";
			#$lambdaBin -t $BlastCores -id 75 -nm 200 -p blastn -e 1e-8 -so 7 -sl 16 -sd 1 -b 5 -pd on -q $query -d $DB -o $tmptaxblastf -oc $outcols;
			if ($doQuery || $doInter){
				$cmd .= "$lambdaBin $defLopt -q $query -i $DB.lba.gz  -o $taxblastf\n";
			}
			
			
			#$cmd .= "\nmv $tmptaxblastf $taxblastf\n";
			#die "$doQuery \n$cmd\n";
			print "\n\n$cmd\n\n";
			systemW($cmd);
			#} else {	print "Blast output $taxblastf2 does exist\n";}}
			#systemW $cmd ;#or die "\n$cmd\n failed\n";
		} elsif ($doblast==3) {
			$simName = "smRNA";

			my $idxDir = getProgPaths("${MKname}idx", 0);
			my $idxStr = ($idxDir ne "") ? "--idx-dir '$idxDir' --index 0" : "";
			my $cmd = "";
			if ($doInter && -e $interLeave && !-z $interLeave) {
				print "Running sortmerna on interleaved\n";
				my $pfx = "$tmpD/smrna_i$DBi";
				$cmd .= "$smrnaBin --ref '$DB' --reads '$interLeave' $idxStr --blast 1 --threads $BlastCores -e 0.1 --num_alignments 50 --paired_in --aligned '$pfx' --kvdb '${pfx}.kvdb' --readb '${pfx}.readb'\n";
				$cmd .= "rm -rf '${pfx}.kvdb' '${pfx}.readb'\n";
				$cmd .= "mv '${pfx}.blast' '$taxblastf2'\n";
			}
			if ($doQuery) {
				print "Running sortmerna on merged reads\n";
				my $pfx = "$tmpD/smrna_q$DBi";
				$cmd .= "$smrnaBin --ref '$DB' --reads '$query' $idxStr --blast 1 --threads $BlastCores -e 0.1 --num_alignments 50 --aligned '$pfx' --kvdb '${pfx}.kvdb' --readb '${pfx}.readb'\n";
				$cmd .= "rm -rf '${pfx}.kvdb' '${pfx}.readb'\n";
				$cmd .= "mv '${pfx}.blast' '$taxblastf'\n";
			}
			systemW $cmd;
			print "done sortmerna assignment\n";
		} elsif ($doblast==4) {
			my $udbDB = $DB . ".vudb";
			unless (-e $udbDB){
				my $temporaryUdb = "$udbDB.$$";
				systemW("$vsearchBin --makeudb_usearch ".shellQuote($DB)." --output ".shellQuote($temporaryUdb));
				if (-e $udbDB){
					unlink $temporaryUdb or die "Cannot remove redundant VSEARCH index $temporaryUdb: $!\n";
				} else {
					move($temporaryUdb, $udbDB) or die "Cannot install VSEARCH index $udbDB: $!\n";
				}
			}
			my $cmd = "$vsearchBin ";
			my $outCols="query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+ql";
			$cmd .= "--usearch_global ".shellQuote($query)." --db ".shellQuote($udbDB)." --id 0.75 --query_cov $lengthTolerance --userfields $outCols --userout ".shellQuote($taxblastf)." --maxaccepts 100 --maxrejects 100 --strand both --threads $BlastCores;";
			systemW $cmd;
		}

	}
		#die;
	print "Running LCA..\n";
	my $LCcmd = "";
	#$hof1 = writeBlastHiera($fullBlastTaxR,$id,$simName);  undef $BlastTaxR ;
	my $LCAflags = "-LCAfrac 0.8 -showHitRead -cover $lengthTolerance -minAlignLen 50 ";
	$LCAflags .= "-reportID " if ($extractDNA);
	$LCcmd .= "$LCAbin -i ". join(",",@taxouts) . " -r ".join(",",@DBtaxa_used)." -o $hof1 $LCAflags\n";
	
	
	#merging of interleave results
	if (0&&$doInter){
		$LCcmd .= "$LCAbin -i ". join(",",@taxouts_inter) . " -r ".join(",",@DBtaxa_used)." -o $hof2 $LCAflags\n";
		$LCcmd.= "tail -n+2 $hof2 >> $hof1; rm $hof2;"; #remove header from LCA
	}
	
	print $LCcmd."\n";
	systemW $LCcmd;
	die "LCA command completed without producing $hof1\n" unless -e $hof1;

	#extract also reads?
	extractReads($query,$hof1,"$outdir/$id.extr.fa");
	unlink $_ for grep { -e $_ } (@taxouts, @taxouts_inter);
	unlink $_ for glob($queryO."__U*");
	if ($interLeaveO ne ""){
		unlink $_ for glob($interLeaveO."*");
	}
}

#file operations on paired end files
sub merge($ $ $){
	my ($r1,$r2,$interleave) = @_;
	$r2 = fastq2fna($r2,1);
	$r1 = fastq2fna($r1,1); 
	print "Combining unmerged read pairs..\n";
	my %forward = %{readFasta($r1,1)};
	my %reverse = %{readFasta($r2,1)};
	my %reverseByPair;
	for my $id (keys %reverse){
		my $pairId = $id;
		$pairId =~ s/\/2$//;
		die "Duplicate reverse-read identifier after pair normalization: $pairId\n"
			if exists $reverseByPair{$pairId};
		$reverseByPair{$pairId} = $reverse{$id};
	}
	open my $combinedHandle, ">", $interleave or die "Cannot write combined reads $interleave: $!\n";
	for my $id (sort keys %forward){
		my $pairId = $id;
		$pairId =~ s/\/1$//;
		die "Cannot find reverse mate for $id\n" unless exists $reverseByPair{$pairId};
		print {$combinedHandle} ">$pairId\n$forward{$id}", reverse_complement_IUPAC($reverseByPair{$pairId}), "\n";
		delete $reverseByPair{$pairId};
	}
	die "Reverse-read file contains unmatched identifiers: ".join(",", sort keys %reverseByPair)."\n"
		if keys %reverseByPair;
	close $combinedHandle or die "Cannot close combined reads $interleave: $!\n";
}

sub extractReads($ $ $){
	my ($query,$hof1,$ofile) = @_;
	#die "$query $hof1\n";
	return if (!$extractDNA);
	my %FNA ;
	if ($extractDNA){
		my $hr = readFasta($query);%FNA = %{$hr};
	}
	open I,"<$hof1" or die "Can't open hier file $hof1\n";
	open O,">$ofile" or die "Can not open stream to $ofile\n";
	while (my $l = <I>){
		chomp $l;
		my @spl = split /\t/,$l;
		unless (exists($FNA{$spl[0]})){print "no entry for $spl[0]\n" ; next;}
		print O ">$spl[0] $spl[$#spl]\n$FNA{$spl[0]}\n";
	}
	close I;
	close O;
	systemW("gzip -f ".shellQuote($ofile));
	#die "$ofile\n";
}


sub materializeReadInput{
	my ($plainPath, $workDir) = @_;
	return $plainPath if -e $plainPath;
	my $gzipPath = "$plainPath.gz";
	return "" unless -e $gzipPath;
	my $inputDir = File::Spec->catdir($workDir, "inputs");
	make_path($inputDir) unless -d $inputDir;
	my $temporary = File::Spec->catfile($inputDir, basename($plainPath));
	systemW("gzip -dc ".shellQuote($gzipPath)." > ".shellQuote($temporary));
	die "Could not materialize compressed read input $gzipPath\n" unless -e $temporary;
	return $temporary;
}


sub resultExists{
	my ($path, $allowEmpty) = @_;
	return 1 if -e $path && ($allowEmpty || -s $path);
	return 1 if -e "$path.gz" && ($allowEmpty || -s "$path.gz");
	return 0;
}


sub touchFile{
	my ($path) = @_;
	open my $touchHandle, ">", $path or die "Cannot create checkpoint $path: $!\n";
	close $touchHandle or die "Cannot close checkpoint $path: $!\n";
}


sub countLines{
	my ($path) = @_;
	open my $lineHandle, "<", $path or die "Cannot count records in $path: $!\n";
	my $count = 0;
	$count++ while <$lineHandle>;
	close $lineHandle or die "Cannot close $path: $!\n";
	return $count;
}


sub shellQuote{
	my ($value) = @_;
	$value = "" unless defined $value;
	$value =~ s/'/'"'"'/g;
	return "'$value'";
}

#function that goes line by line through Blast m8 table
sub getTaxForOTUfromRefBlast($ $ $ $){
	my ($blastout,$GGref,$interLMode,$BlastCores) = @_;
	#sp,ge,fa,or,cl,ph
	my %GG = %{$GGref};
	my @ggk = keys(%GG);
	my $maxGGdep = scalar(@{$GG{$ggk[0]}});
	@ggk = ();
	open B,"<",$blastout or die "Could not read $blastout\n";
	my $sotu = "";my $sID=0 ; my $sLength=0; 
	#my @sTax=(); my @sMaxTaxNum = ();
	my %retRef ;# my $retDRef ={};
	my $cnt=0;
	my @tmpLines = (); #stores Blast lines
	my @spl; #temp line delim
	my $minBit = 120; my $minEval = 1e-14;
	my %prevQueries = ();
	my $line = "";
	my @thrs; my $thrsCnt = 0;
	#for (my $i=0;$i<$BlastCores;$i++){$thrs[$i] = threads->create(\&empty_proc);}
	while ($line = <B>){
		$cnt++; chomp $line; #$line2 = $line;
		my @spl = split("\t",$line);
		my $totu = $spl[0]; #line otu
		$totu =~ s/^>//;
		if ($cnt == 1) {$sotu = $totu;}
		#print $line." XX $spl[11] $spl[10]\n"; die () if ($cnt > 50);
		
		#check if this is a 2 read-hit type of match (interleaved mode) & merge subsequently
		if ($interLMode){
			if (@tmpLines>0 && exists($prevQueries{ $spl[1]}) ){
				my @prevHit = @{$prevQueries{$spl[1]}};
				die "something went wrong with the inter matching: $prevHit[1] - $spl[1]\n" unless ( $prevHit[1] eq $spl[1] );
				$prevHit[11] += $spl[11];#bit score
				$prevHit[3] += $spl[3];$prevHit[5] += $spl[5];$prevHit[4] += $spl[4];#alignment length,mistmatches,gap openings
				$prevHit[2] = ($prevHit[2] + $spl[2]) / 2;
				#$tmpLines[-1] = \@prevHit;
				@spl = @prevHit;
				#next;
			} else {$prevQueries{ $spl[1] } = \@spl;}
		}
		if (($spl[11] < $minBit) || ($spl[10] > $minEval) ){ #just filter out..
			#print "ss\n"; 
			#next;
			$spl[2] =0; #simply deactivate this way...
		}
		
		if ($sotu eq $totu){
			push(@tmpLines,\@spl);
			if ($spl[2] > $sID && $spl[3] > $sLength){$sID = $spl[2]; $sLength = $spl[3];}
			if ($spl[3] > ($sLength*1.4) && $spl[2] > ($sID*0.9)) {$sID = $spl[2]; $sLength = $spl[3];} #longer alignment is worth it..
			#print $sID."\n";
		} else {
			#print "DD";
			#print "Maybe\n";
			
#			if (0 && $sotu ne ""){				if ($thrsCnt>= $BlastCores){$thrsCnt=0;}				#get result of old job
#				my $ret = $thrs[$thrsCnt]->join();				foreach (keys %{$ret}){$retRef{$_} = ${$ret}{$_};}				$thrs[$thrsCnt] = threads->create(\&rework_tmpLines,\@tmpLines,$sotu,$sID,$sLength,\%GG,$maxGGdep,$sMaxTax) ;
#				$thrsCnt++; print $thrsCnt." ";			}
			if ($sotu ne ""){
				my ($AR) = rework_tmpLines(\@tmpLines,$sID,$sLength,\%GG,$maxGGdep,); 
				$retRef{$sotu} = $AR;
			}
			$sotu = $totu; undef @tmpLines ; undef %prevQueries ;
			push(@tmpLines,\@spl);$prevQueries{ $spl[1] } = \@spl;
			$sID =  $spl[2]; $sLength = $spl[3];
		}
	}
	#last OTU in extra request
	if ($sotu ne ""){
		#my @spl = split("\t",$line);
		
		my ($AR) = rework_tmpLines(\@tmpLines,$sID,$sLength,\%GG,$maxGGdep);
		$retRef{$sotu} = $AR;
#		my($sTaxX,$sMaxTaxX) = LCA(\@sTax,\@sMaxTaxNum,$maxGGdep);
#		$ret{$sotu} = $sTaxX;
#		$retD{$sotu} = $sMaxTaxX;
	}
	close B;
	#debug 
	#my %ret = %{$retRef};	my @tmp = @{$ret{$sotu}};print "\n@tmp  $sotu\n";
	undef  @tmpLines;undef %prevQueries ;
	#print "Assigned $refDBname Taxonomy to OTU's\n",0;
	
	return (\%retRef);
}
