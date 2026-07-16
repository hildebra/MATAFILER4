package Mods::SNP;
use warnings;
use strict;

use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::GenoMetaAss qw( gzipopen systemW filsizeMB readFasta readGFF writeFasta reverse_complement_IUPAC fileGZs fileGZe );
use Mods::Subm qw(qsubSystem emptyQsubOpt qsubSystem2);
use Mods::TamocFunc qw (cram2bsam);
use List::Util qw/shuffle/;


use Exporter qw(import);
our @EXPORT_OK = qw(SNPconsensus_vcf SVcall_vcf);

sub regionsFromFAI($){
	my ($inF ) =@_;
	my @ret;
	open I,"<$inF" or die "can't open fai $inF\n";
	while ( my $line = <I>){
		chomp $line;
		my @fields = split /\t/,$line;
		push(@ret,$fields[0] . ":1-" . $fields[1]);
	}
	close I;
	return (@ret);
}

sub getRegionsBamDepth{
	my ($depthPC,$totalSpl,$maxSNPcores) = @_;
	my %LperC; my %DperC; my $tDep=0; my %contigNum; my $tLen=0;
	my $cnt=0;
	#print "$depthPC\n";
	#open I ,"<$depthPC" or die "can;t oopen depth file $depthPC\n";
	my ($IN ,$status) = gzipopen($depthPC,"contig depth file",0);
	return ([],[]) unless ($status && defined $IN);
	while (<$IN>){
		chomp; my @spl = split /\t/;
		
		if ($spl[0] =~ m/__C(\d+)_L=(\d+)=/){
			$contigNum{$1} = $spl[0];
			$LperC{$spl[0]} = $2;
			$tLen += $2;
			$tDep += $2 * $spl[1] ;
		} else {
			$contigNum{$cnt} = $spl[0];
			$LperC{$spl[0]} = 1000;
			$tLen += 1000;
			$tDep += 1000 * $spl[1] ;
		}
		$DperC{$spl[0]} = $spl[1];
		$cnt++;
	}
	close $IN;
	#some hardcoded rules to make small samples into smaller job numbers..
	if ($tDep ==0 || $tLen ==0){
		return [],[];
	}
	
	if ($tDep <1e6){
		$totalSpl = 1;
	} elsif ($tDep <5e6){
		$totalSpl = int($maxSNPcores/3);
	} elsif ($tDep <20e6){
		$totalSpl = int($maxSNPcores/2);
	} else {$totalSpl = int($maxSNPcores/2);}
	$totalSpl = 1 if ($totalSpl < 1);
	#print "totalSpl $totalSpl $maxSNPcores $tDep\n";
	#expected depth per bin
	my $exD = $tDep/$totalSpl;
	#print "$tDep/$tLen\n";
	my $avgD = $tDep/$tLen;
	#now count up per bin to get to this number..
	my $curD = 0; my @regions; $cnt =0;
	my @regOrd;
	my @srtK = sort {$a <=> $b} keys %contigNum;
	my $startCtg =0;
	my $idx=0;
	foreach my $id (@srtK){
		push (@regOrd,$contigNum{$id});
	}
	while ( $idx<@srtK ){
		#print "$idx ";
		my $k1 = $srtK[$idx];
		my $k = $contigNum{$k1};
		my $thisD = $DperC{$k} * ($LperC{$k} - $startCtg);
		#print "$idx ".@srtK." $curD+$thisD < $exD $startCtg\n";
		
		if ($DperC{$k} > $avgD*0.85){ #deep sample, uses more processing power..
			$thisD *=  1.2;
		} 
		if ($DperC{$k} < $avgD*0.1){ #deep sample, uses more processing power..
			$thisD *=  0.8;
		}
		#if ($k =~ m/MM2__C152_L=33826=/){die "$curD + $thisD) < $exD \n";}
		if ( ($curD + $thisD) > $exD ){
			#die "$curD + $thisD) < $exD \n";
			if ($LperC{$k} < 1000){
				$regions[$cnt] .= "$k\t0\t$LperC{$k}\n";
				$curD = 0; $cnt ++; $startCtg = 0; $idx++;
			} else {
				#how much into the contig?
				my $stopCtg = int($LperC{$k} * (($exD - $curD) / $thisD) + $startCtg);
				if ($stopCtg < 150){#just don't take this, reset counter..
				#print "A";
					#$curD = $thisD; $cnt ++;$regions[$cnt] .= "$k\t0\t$LperC{$k}\n";
					$startCtg = 0;$curD =0; $idx ++; 
				} elsif ( ($LperC{$k} - $stopCtg) < 150){ #just take whole...
				#print "B";
					$regions[$cnt] .= "$k\t$startCtg\t$LperC{$k}\n"; $curD = 0; $cnt ++; $startCtg = 0; $idx++;
				} else { #take part of contig, reset counter..
				#print "C";
					#die "$stopCtg : $LperC{$k} * ($exD - $curD) / $thisD;\n$k\t0\t$stopCtg\n$k\t$stopCtg\t$LperC{$k}\n";
					$regions[$cnt] .= "$k\t$startCtg\t$stopCtg\n";$startCtg = $stopCtg;
					$curD = 0; $cnt ++; 
				}
			}
			#reset
		} else {
			$regions[$cnt] .= "$k\t$startCtg\t$LperC{$k}\n"; $idx++;
			$curD += $thisD; $startCtg = 0;
		}
	}
	#die @regions." regions\n";
	return (\@regions,\@regOrd);
}

sub getRegionsBam{
	my ($splitFAsizeL,$refFA,$tmpD) = @_;
	
	my $smtBin = getProgPaths("samtools");
	my $py_genRegs = getProgPaths("genRegions_scr"); #python $frDir/fasta_generate_regions.py
	my $regionFile = "$refFA.reg";
	systemW "$smtBin faidx $refFA" unless (-e "$refFA.fai");
	my @curReg = regionsFromFAI("$refFA.fai");
	#die "@curReg\n";
	#if (!-e $regionFile ){		systemW "python $py_genRegs $refFA.fai $splitFAsizeL > $regionFile\n"	}
	#open my $handle, '<', $regionFile;	chomp(@curReg = <$handle>);	close $handle;
	
	#now comes the real work: translate to bed format
	#@curReg = shuffle @curReg;
	my @regions; my @regOrd;
	my $curCnt = 0; my $cnt=0;
	foreach my $reg (@curReg){
		my @spl = split(/:/,$reg);
		push(@regOrd,$spl[0]);
		my @spl2 = split(/-/,$spl[1]);
		my $sum =  $spl2[1]-$spl2[0];
		if ($curCnt + $sum > $splitFAsizeL){
			#how much can still be added?
			my $max = $curCnt + $sum - $splitFAsizeL;
			my $start = $spl2[0];
			if ( (($spl2[0]+$max-1)) + 100 > $spl2[1]){
				#print $spl2[0]+$max-1 .":".$spl2[1]."\n";
				$regions[$cnt] .= "$spl[0]\t$spl2[0]\t$spl2[1]\n";
				$cnt++;
				$curCnt=0; 
			} else {
				$start = ($spl2[0]+$max);
				$regions[$cnt] .= "$spl[0]\t$spl2[0]\t". ($start) ."\n";
				$cnt++;
				$regions[$cnt] .= "$spl[0]\t".($start)."\t".$spl2[1]."\n";
				$curCnt = $spl2[1] - $start;
			}
			#print "$regions{$cnt}\n";
			#die "new\n$regions{$cnt}\n";
		} else {
			$regions[$cnt] .= "$spl[0]\t$spl2[0]\t$spl2[1]\n";
			$curCnt += $sum;
			
		}
	}
#	die "@regions\n";
	return (\@regions,\@regOrd);
}


sub pileupcall{
	
	my ($tarR,$tag,$SNPIHR,$QSBoptHR,$scrDir,$tmpOut,$myParL,$crAR, $varsOnly) = @_;
	
	my @curReg = @{$crAR};
	
	#die "@curReg\n";
	my $bcftBin = getProgPaths("bcftools");
	my $refFA = $SNPIHR->{assembly};
	my $tmpdir = $SNPIHR->{nodeTmpD};
	my $smplNm = $SNPIHR->{smpl};
	my $qsubDirE = $SNPIHR->{qsubDir};
	my $runLocalTmp = $SNPIHR->{runLocal};
	my $x = $SNPIHR->{JNUM};
	my $threads = 0;#$SNPIHR->{threads};
	my $run2ctg = $SNPIHR->{run2ctg};
	my $rdep = $SNPIHR->{rdep} ;
	my $normalizeIndel = 1; $normalizeIndel = $SNPIHR->{normIndels} if (exists($SNPIHR->{normIndels})) ;


	my $useFB = 1;	$useFB = 0 if (uc($SNPIHR->{SNPcaller}) eq "MPI");
	my $overwrite = $SNPIHR->{overwrite};
	my $deferRegionPlanning = $SNPIHR->{deferRegionPlanning} || 0;

	#basic caller options..
	my $minBQ=30; my $minMQ=30.0;
	#freebayes std options
	my $frAllOpts= "-u -i -m $minMQ -q $minBQ -C 1 -F 0.1 -k -X --pooled-continuous --report-monomorphic  --min-repeat-entropy 1 --use-best-n-alleles 2 -G 1 ";
	#bcftools options #-q = map qual -Q = base qual
	my $bcfAllOpts = " --min-BQ $minBQ -d 12000 --threads $threads  --min-MQ $minMQ -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP,INFO/FS,INFO/IDV,INFO/MQ0F,INFO/BQBZ,INFO/SCBZ,INFO/RPBZ,INFO/MQBZ"; 
	#--skip-indels #--count-orphans -a DP,AD,ADF,ADR,SP
	my $bcftCallOpts = " --ploidy 1 -c -M --threads 0 -a INFO/PV4 "; #--multiallelic-caller  -> replaced with -c (consensus caller)
	$bcftCallOpts .= "--variants-only " if ($varsOnly);
	#$bcftCallOpts .= "-O z " ;# if ($normalizeIndel);
	my $bcftViewOpts = "-Oz  ";
	if ($tag eq "sup-"){
		if ($SNPIHR->{SeqTechSuppl} eq "ONT"){$bcfAllOpts.=" -X ont ";
		} elsif ($SNPIHR->{SeqTechSuppl} eq "PB"){$bcfAllOpts.=" -X pacbio-ccs ";
		} else{$bcfAllOpts.=" -X illumina ";}
	} else {
		if ($SNPIHR->{SeqTech} eq "ONT"){$bcfAllOpts.=" -X ont ";
		} elsif ($SNPIHR->{SeqTech} eq "PB"){$bcfAllOpts.=" -X pacbio-ccs ";
		} else{$bcfAllOpts.=" -X illumina ";}
	}
	
	my $cmd = "";  
	my $locXtrCmd = ""; $locXtrCmd = " &\npids+=(\$!)" if ($runLocalTmp);
	#my $tag = "primary";
	#
	#$cmd .= "mkdir -p $tmpdir\n";
	my $cmdAll2 = "";
	$cmdAll2 .= "echo \"Processing bams - mpileup $tag\"\n";
	$cmdAll2 .= "pids=()\n" if ($runLocalTmp);
	if ($useFB){
		my $frbBin = getProgPaths("freebayes");
		$cmd = "ulimit -s unlimited\n$frbBin -f $refFA  $frAllOpts ";
	} else {
		$cmd = "$bcftBin mpileup --fasta-ref $refFA $bcfAllOpts ";
	}
	my $normPreTag = ""; #$normPreTag = ".pre" if ($normalizeIndel);
	my @allDeps2; my @checkF=();
	#implement in parallel as too slow in single core mode :/
	system "rm -f $tmpOut.$tag*" if ($overwrite);
	system "mkdir -p $qsubDirE" unless (-d $qsubDirE);
	my $bedJobs = 0;
	for (my $i=0;$i<@curReg;$i++){ #go over regions in bed file, submit a job for each "region"
		#$tar[0] = bam file;  $bedF = bedfile with regions
		next if (!$run2ctg);
		$bcftViewOpts .= "--no-header " if (scalar(@checkF)==1);
		system "rm -f $qsubDirE/$smplNm.$tag*.bed" if ($i==0 && !$deferRegionPlanning) ;
		my $cmd2 = $cmd ;
		my $bedF = $qsubDirE."$smplNm.${tag}$i.bed";
		push @checkF, $bedF;
		next if (!$deferRegionPlanning && !-e $bedF && -e "$tmpOut.$tag$i" && !$overwrite);
		if ($myParL){
			unless ($deferRegionPlanning) {
				open O,">",$bedF or die $!;print O $curReg[$i];close O;
			}
			if ($useFB){
				$cmd2 .= " -t $bedF $tarR->[0] > $tmpOut.$tag$i && rm $bedF $locXtrCmd\n"; #--region '$curReg[$i]'
			}else{
				$cmd2 .= " -R $bedF $tarR->[0] | $bcftBin call $bcftCallOpts | $bcftBin view $bcftViewOpts  > $tmpOut.$tag$i$normPreTag.gz  && rm $bedF $locXtrCmd\n"; 
			}
		} else {
			die "incomplete control structure SNP.pm\n";
		}
		$cmdAll2 .= $cmd2."\n" if ($run2ctg);
		if (!$runLocalTmp){
			my ($dep,$qcmd) = qsubSystem($qsubDirE. $SNPIHR->{cmdFileTag} . ".ac$tag.$smplNm.$i.sh",$cmd2,1,"15G","FBC$x.$i",$rdep,"",1,[],$QSBoptHR);
			push (@allDeps2,$dep);
		}
		$bedJobs++;
		#last if ($i == 1000);
	}
	#$bedJobs =1 if ($bedJobs<1);
	$cmdAll2 .= "status=0\nfor pid in \"\${pids[\@]}\"; do if ! wait \"\$pid\"; then status=1; fi; done\ntest \"\$status\" -eq 0\n"
		if ($run2ctg && $runLocalTmp);
	$cmdAll2 .= "rm -f $tarR->[0].crai $tarR->[0].bai;\n";
	
	
	
	$cmdAll2 .= "\necho \"Finished mpileup $tag\"\n\n";	
	#$cmdAll .= "rm -f $tarR->[0];\n" if ($run2ctg && $tarR->[0] ne "");
	if ($bedJobs ==0){$cmdAll2="";}
	return (\@allDeps2, \@checkF, $cmdAll2);
}


sub SNPconsensus_vcf{
	my ($SNPIHR)  = @_;
	#my $kpathBin = getProgPaths("kpath");
	my $smtBin = getProgPaths("samtools");
	my $vcfcnsScr = getProgPaths("vcfCons_scr");
	#my $tabixBin = getProgPaths("tabix");
	#my $vcfLD = getProgPaths("vcfLib_dir");#"/g/bork3/home/hildebra/bin/vcflib/bin/";
	my $ctg2fas = getProgPaths("contig2fast_scr");
	my $pigzBin  = getProgPaths("pigz");
	#my $py3 = getProgPaths("py3activate",0);my $py3d = getProgPaths("pydeacti",0);
	my $bcftBin = getProgPaths("bcftools");
	my $vcf2fnaBin = getProgPaths("vcf2fna");
	
	my $memPJob =0;
	#get parameteres
	my $samcores = 12;
	#my %SNPinfo = %{$SNPIHR};
	my $QSBoptHR = $SNPIHR->{QSHR};
	my $immediateSubm = exists($SNPIHR->{immediateSubm}) ? $SNPIHR->{immediateSubm} : 1;
	my $submissionCommands = "";
	my $x = $SNPIHR->{JNUM};
	my $jdep = ""; $jdep = $SNPIHR->{jdeps} if (exists($SNPIHR->{jdeps}));
	my $tmpdir = $SNPIHR->{nodeTmpD};
	my $smplNm = $SNPIHR->{smpl};
	my $refFA = $SNPIHR->{assembly};
	my $qsubDirE = $SNPIHR->{qsubDir};
	my $scrDir = $SNPIHR->{scratch};
	my $bamcram = $SNPIHR->{bamcram};
	my $splitFAsize = $SNPIHR->{bpSplit};
	my $overwrite = $SNPIHR->{overwrite};
	my $runLocalTmp = $SNPIHR->{runLocal};
	my $maxSNPcores= $SNPIHR->{maxCores};
	my $hasPrimaryRds = $SNPIHR->{hasPrimaryRds};
	my $normalizeIndel = 1; $normalizeIndel = $SNPIHR->{normIndels} if (exists($SNPIHR->{normIndels})) ;
	my $deferRegionPlanning = $SNPIHR->{deferRegionPlanning} || 0;
	my $onlyNormalize = 0;

	$memPJob = $SNPIHR->{memPJob};
	my $SNPstone = $SNPIHR->{STOconSNP}; my $SNPsuppStone = $SNPIHR->{STOconSNPsupp};
	my $minDepth = 0;
	$minDepth =$SNPIHR->{minDepth}  if (exists($SNPIHR->{minDepth} ));
	my $minCallQual = $SNPIHR->{minCallQual};#20;
	#my $SNPstone = $ofasConsDir."SNP.cons.stone";
	#my $memReq = "20G";
	#my $memReq = $SNPIHR->{memReq};
	my $vcfFile = ""; $vcfFile = $SNPIHR->{vcfFile} if (exists ($SNPIHR->{vcfFile}));
	my $vcfFileS = ""; $vcfFileS = $SNPIHR->{vcfFileSupp} if (exists ($SNPIHR->{vcfFileSupp}));
	my $cmdFTag = $SNPIHR->{cmdFileTag};
	my $firstInSample = 0;$firstInSample = $SNPIHR->{firstInSample} if (exists($SNPIHR->{firstInSample}));
	my $useFB = 1;	$useFB = 0 if (uc($SNPIHR->{SNPcaller}) eq "MPI");
	$vcfcnsScr = getProgPaths("vcfCons_FB_scr") if ($useFB);
	my $actualCores  = $maxSNPcores;
	if ($runLocalTmp){
		$scrDir = $tmpdir;
	}
	

	my $saveVCF=1;
	if ($vcfFile eq ""){
		$saveVCF=0;
		$vcfFile = "$scrDir/$smplNm.fin.vcf.gz"; #switch back to vcf for vcf2fna
		$vcfFileS = "$scrDir/$smplNm.fin-sup.vcf.gz";
	}
	if ($normalizeIndel){
		my $SNPstone2 = $SNPstone; $SNPstone2 =~ s/\.norm\.stone/\.stone/;
		if (! -e $SNPstone && -e $SNPstone2){#can just do the normalization..
			$onlyNormalize = 1 ; system "rm -f $SNPstone2"; $maxSNPcores=1;
		}
		#die "$onlyNormalize $SNPstone2 $SNPstone\n";
	}	
	


	
	#key change for C++ program vcf2fasta
	my $reportVCFonly=1; #was 0 before


	my $ofasCons = $SNPIHR->{ofas};
	$ofasCons =~ m/(^.*)\/[^\/]+$/;
	my $ofasConsDir = $1."/";
	#$ofasCons .= ".gz" unless ($ofasCons =~ m/\.gz$/); #don't change, needed without..
	my $run2ctg=1; #flag to determine if I run the cram to bam, mpileup, consensus contig steps..
	system "rm -f $ofasConsDir/*" if ($overwrite);
	#die "$ofasCons\n";
	if (fileGZe ($vcfFile) || ( -e "$ofasCons.gz" && -s "$ofasCons.gz" > 200)){
		$run2ctg =0 ;
	} else {
		system "rm -f $ofasConsDir/*"; #better safe than sorry..
	}
	#first all important regions on finalDir
	my @curReg = ("1"); my @regOrd;
	my $myParL=0;
	if ($splitFAsize>0){$myParL=1;}
	if ($myParL && $run2ctg && $deferRegionPlanning) {
		my $runtimeJobs = $SNPIHR->{split_jobs} || $maxSNPcores || 1;
		$runtimeJobs = $maxSNPcores if ($runtimeJobs > $maxSNPcores);
		$runtimeJobs = 1 if ($runtimeJobs < 1);
		@curReg = (('runtime') x $runtimeJobs);
	} elsif ($myParL && $run2ctg){ #no, don't redo freebayes part
		my ($refAR,$refAR2);
		if ($hasPrimaryRds && -e $SNPIHR->{depthF} ne ""){
			 ($refAR,$refAR2) = getRegionsBamDepth($SNPIHR->{depthF},$SNPIHR->{split_jobs},$maxSNPcores);
		} else {
			 ($refAR,$refAR2) = getRegionsBam($splitFAsize,$refFA,$tmpdir);
		}
		@curReg = @{$refAR}; 
		#@regOrd = @{$refAR2};  #use .fai instead for this..
		if (@curReg == 0){return("");}
		
		#open O,">$refFA.reg" or die "can't open region file $refFA.reg\n";		print O join("\n",@regOrd);		close O;
		#die "$refFA.reg\n@curReg\n";
	}
	if ($runLocalTmp){
		$actualCores = scalar(@curReg);
		$samcores = $actualCores;#$SNPIHR->{split_jobs};
	}
	
	my $rdep=$jdep;
	#prepare files..
	my $cleanCmd = ""; 
	my $xtra = "";
	$xtra .= "echo \"Preparing data\"\n";
	$xtra .= "mkdir -p $scrDir;\n";
	$xtra .= "$smtBin faidx $refFA;\n" unless (-e "$refFA.fai");
	if ($deferRegionPlanning && $myParL && $run2ctg) {
		my $runtimeJobs = scalar(@curReg);
		my $bedPrefix = "$qsubDirE/$smplNm.";
		$xtra .= "rm -f ${bedPrefix}*.bed\n";
		$xtra .= "awk -v OFS='\\t' -v n=$runtimeJobs -v chunk=$splitFAsize -v prefix='$bedPrefix' 'BEGIN { i=0 } { for (s=0; s<\$2; s+=chunk) { e=s+chunk; if (e>\$2) e=\$2; f=prefix (i%n) \".bed\"; print \$1, s, e >> f; i++ } }' $refFA.fai\n";
		$xtra .= "for i in \$(seq 0 ".($runtimeJobs - 1)."); do test -s ${bedPrefix}\${i}.bed || exit 42; done\n";
		if ($SNPsuppStone ne "") {
			$xtra .= "for i in \$(seq 0 ".($runtimeJobs - 1)."); do cp ${bedPrefix}\${i}.bed ${bedPrefix}sup-\${i}.bed; done\n";
		}
	}
	#$xtra .= "exit\n"; #DEBUG
	#$xtra .= "cp $refFA $refFA.fai $scrDir;\n";$refFA =~ m/\/([^\/]+$)/;$refFA = "$scrDir/$1";
	#my $preTar = 
	
	$xtra .= "echo \"Creating c/bams indexes primary reads\"\n";
	#my @tar = ("");$tar[0] = ${$SNPIHR->{MAR}}[0]; #$preTar;
	
	my @tar= @{$SNPIHR->{MAR}}[0] if ($hasPrimaryRds);my $cmdAll = "";my @allDeps2; my @checkF;
	my $tmpOut = "$scrDir/$smplNm.cons.vcf";my $depthFile ="";

	#supplementary mappings?
	my @tarS = ();@tarS = @{$SNPIHR->{MARsupp}} if ($SNPsuppStone ne "" && exists ($SNPIHR->{MARsupp} ) );
	my $tmpOut2 = "$scrDir/$smplNm.X.cons.vcf";my $depthFileS  = "";
	
	my $memReqGB = 20; #memory requested overall
	my $bcramSiz = 0; 
	if ($hasPrimaryRds){$bcramSiz = filsizeMB($tar[0]);} 
	if (@tarS){my $tmpSI = filsizeMB($tarS[0]); $bcramSiz = $tmpSI if ($tmpSI > $bcramSiz);}
	my $refSize = filsizeMB($refFA);
	my @limits = (1500,3500,5500,7500,10000,12000); my @memRperLimit = (15,20,30,40,60,120);
	if ($bamcram eq "bam"){ for (my $i=0;$i<@limits;$i++){$limits[$i] *= 1.7;}} #increase limits for bams..
	for (my $i=0;$i<@limits;$i++){ last if (($bcramSiz+$refSize)<$limits[$i]); $memReqGB = $memRperLimit[$i];}
	$memReqGB = $memPJob if ($memPJob > 0);
	
	
	#die "hasPrimaryRds: $hasPrimaryRds\n";

	if ($hasPrimaryRds){ #primary reads SNP call
		
		die "Can't find input file $tar[0] (SNP.pm)\n"
			unless ($deferRegionPlanning || -e $tar[0]);
		if ($bamcram eq "cram"){ #create index for bam/cram
			$xtra .= "if [ ! -e $tar[0].crai ] || [ ! -s $tar[0].crai ]; then rm -f $tar[0].crai; $smtBin index -@ $samcores  $tar[0]; fi\n";
		} else {
			$xtra .= "if [ ! -e $tar[0].bai ] || [ ! -s $tar[0].bai ]; then rm -f $tar[0].bai; $smtBin index -@ $samcores  $tar[0]; fi\n";
		}
		
		#find depthfil for input bam (primary)
		$depthFile = $tar[0];$depthFile =~ s/\.cram$|\.bam$/\.bam\.coverage\.gz/;
		if (!-e $depthFile && !$deferRegionPlanning){$depthFile =~ s/\.gz//; die "no depth file found (SNP.pm): $depthFile\n$tar[0]\n" if (!-e $depthFile);}
		if (!$runLocalTmp && $run2ctg && (!-e $tar[0] || !-e $refFA) ){
			my ($dep,$qcmd) = qsubSystem($qsubDirE."$cmdFTag.CramToBam$x.sh",$xtra,2,"17G","CtB$x",$jdep,"",$immediateSubm ? $samcores : 0,[],$QSBoptHR);
			$submissionCommands .= $qcmd unless ($immediateSubm);
			$cleanCmd .= "rm -r $scrDir\n";$rdep = $dep;$xtra = "";
		}
		$SNPIHR->{run2ctg} = $run2ctg;$SNPIHR->{rdep} = $rdep;

		#$SNPIHR->{assembly} = $refFA;
		$cmdAll .= $xtra if ($run2ctg && !$onlyNormalize);
		my ($dAR,$cAR,$pilecmd) =  pileupcall(\@tar,"",$SNPIHR,$QSBoptHR,$scrDir,$tmpOut,$myParL,\@curReg,$reportVCFonly);
		@allDeps2 = @{$dAR}; @checkF = @{$cAR};
		$cmdAll .= $pilecmd if (!$onlyNormalize);
	}
	
	if (@tarS){ #supplementary reads SNP call
		my $xtra2 = "echo \"Creating c/bams indexes supplemental reads\"\n";
		
		$depthFileS = $tarS[0];$depthFileS =~ s/\.cram$|\.bam$/\.bam\.coverage\.gz/;
		if (!-e $depthFileS && !$deferRegionPlanning){$depthFileS =~ s/\.gz//; die "no suppl depth file found (SNP.pm): $depthFileS\n$tarS[0]\n" if (!-e $depthFileS);}
		#die "$depthFileS\n";
		if ($bamcram eq "cram"){ #create index for bam/cram
			$xtra2 .= "if [ ! -e $tarS[0].crai ] || [ ! -s $tarS[0].crai ]; then rm -f $tarS[0].crai; $smtBin index -@ $samcores  $tarS[0]; fi\n";
		} else {
			$xtra2 .= "if [ ! -e $tarS[0].bai ] || [ ! -s $tarS[0].bai ]; then rm -f $tarS[0].bai; $smtBin index -@ $samcores  $tarS[0]; fi\n";
		}
		my ($dAR,$cAR,$pilecmd) =  pileupcall(\@tarS,"sup-",$SNPIHR,$QSBoptHR,$scrDir,$tmpOut2,$myParL,\@curReg,$reportVCFonly);
		$cmdAll .= $xtra2.$pilecmd if (!$onlyNormalize);
		push(@allDeps2, @{$dAR}); push(@checkF, @{$cAR});
	}
	
	#all bed files should have been removing inside the two pileupcall() subs
	$cmdAll .= "if ls $qsubDirE/$smplNm.*.bed 1> /dev/null 2>&1 ;then echo \"Bed files still present, probably incorrect run\"; exit 33; else echo \"bed files deleted, looks good\"; fi\n\n" if (@curReg);
		

	
	#from here on: merge XX vcf's into one
	my $sortCmd = "";
	my $vcfSuff=""; #indicates some extra step in merging .. not needed any longer
	$vcfFile .= ".gz" unless ($vcfFile =~ m/\.gz$/);
	if ($myParL && $cmdAll ne ""){
		#this string simply sorts all output files in correct numerical order.. doesn't touch file contents!
		#can create quite a bit of data..
		my $sortedFileList = " | awk -F '.' '{print \$(NF-1),\$0}'  | sort -n -k1 | cut -f2 -d' '";
		$sortCmd .= "mkdir -p $ofasConsDir;\n";
		if ($hasPrimaryRds){
			$sortCmd .= "cat `ls $tmpOut.* $sortedFileList` >$vcfFile ;\nrm -f $tmpOut.*;\n";
		}
		if ($SNPsuppStone ne "" ){
			$vcfFileS .= ".gz" unless ($vcfFileS =~ m/\.gz$/);
			$sortCmd .= "cat `ls $tmpOut2.* $sortedFileList` >$vcfFileS ;\nrm -f $tmpOut2.*;\n";
			
		}
		$cmdAll .= $sortCmd;
	}
	
	my $bcfNormOpts = " -O z --force  -f $refFA "; #-a -m '+both'
	my $cmd3 = "";
	if($normalizeIndel){#bcftools norm -f ref.fa in.vcf
		$cmd3.="\necho \"left-normalizing indels\"\n";
		if ($hasPrimaryRds){
			$cmd3 .= "$bcftBin index $vcfFile\n$bcftBin norm $bcfNormOpts $vcfFile > $vcfFile.norm\n"; # && rm $tmpOut.$tag$i$normPreTag.gz* & 
			$cmd3 .= "rm $vcfFile $vcfFile.csi; mv $vcfFile.norm $vcfFile;\n";
		}
		
		if ($SNPsuppStone ne "" ){
			$cmd3 .= "$bcftBin index $vcfFileS\n$bcftBin norm $bcfNormOpts $vcfFileS > $vcfFileS.norm\n"; 
			$cmd3 .= "rm $vcfFileS $vcfFileS.csi; mv $vcfFileS.norm $vcfFileS;\n";
		}
		
		$cmd3 .= "\necho \"Done normalizing\"\n"; #wait \$(jobs -p);\n
		#die "$cmd3";
		$cmdAll .= $cmd3;
	}
	
	my $postcmd = "";
	if (!$onlyNormalize && ($hasPrimaryRds && fileGZs($vcfFile)) && ($SNPsuppStone eq "" && fileGZs($vcfFileS)) ){$cmdAll="";}
	#die "$cmdAll\n!$onlyNormalize && ($hasPrimaryRds || fileGZs($vcfFile)) && ($SNPsuppStone eq  || fileGZs($vcfFileS)) \n";
	my $vcf2fnaOpt = "";
	my $vcf2fnaOuts = "-oCtg $ofasCons.gz";
	if (!-e $SNPIHR->{genefna}){
		$vcf2fnaOuts .= " -oGeneNT $SNPIHR->{genefna} -oGeneAA $SNPIHR->{genefaa} ";
	}
	die "SNP.pm::SNPconsensus_vcf: gff file not defined" unless (exists($SNPIHR->{gffFile}));
	my $gffF = $SNPIHR->{gffFile}; $gffF .= ".gz" if (-e $gffF . ".gz");
	my $vcf2fnaIns = "-ref $refFA -gff $gffF ";
	
	if (!$hasPrimaryRds){ #only support available..
		my $tmpST = $SNPIHR->{SeqTechSuppl}; if ($tmpST eq ""){$tmpST = "ill";}
		$vcf2fnaOpt = "-seqPlatform $tmpST -t 1 -minCallDepth $minDepth -minCallQual $minCallQual ";
		$vcf2fnaIns .= "-inVCF $vcfFileS -depthF $depthFileS ";
	} elsif ($SNPsuppStone eq "" ){#only primary reads available..
		my $tmpST = $SNPIHR->{SeqTech}; if ($tmpST eq ""){$tmpST = "ill";}
		$vcf2fnaOpt = "-seqPlatform $tmpST -t 1 -minCallDepth $minDepth -minCallQual $minCallQual ";
		$vcf2fnaIns .= "-inVCF $vcfFile -depthF $depthFile ";
	} else {#and for two vcfs..
		$vcf2fnaOpt = "-seqPlatform $SNPIHR->{SeqTech},$SNPIHR->{SeqTechSuppl} -t 1 -minCallDepth $minDepth,$minDepth -minCallQual $minCallQual ";
		$vcf2fnaIns .= "-inVCF $vcfFile,$vcfFileS -depthF $depthFile,$depthFileS ";
	}
	if (!$SNPIHR->{createFastas}){
		$postcmd.="\n##In case you want to create consensus fastas, use (uncomment):\n##$vcf2fnaBin $vcf2fnaOpt $vcf2fnaIns $vcf2fnaOuts\n";
		$postcmd.="#create stats only of hypothetical consensus generations:\n$vcf2fnaBin $vcf2fnaOpt $vcf2fnaIns; \n\n"
	} else {
		$postcmd.="\n# Create consensus fastas\n";
		$postcmd.="$vcf2fnaBin $vcf2fnaOpt $vcf2fnaIns $vcf2fnaOuts;\n\n";
	}

	
	
	$postcmd .= "rm -f $vcfFile$vcfSuff $vcfFileS.csi $vcfFile.csi;\n" if ($vcfSuff ne ""); #functionality no longer used..
	$postcmd .= "rm -f $vcfFileS $vcfFile;\n" if (!$saveVCF && $SNPIHR->{createFastas});
	#} else {
	#	if ($SNPsuppStone ne "" ){die "support reads activated. combined SNP calling only works current with use of the \"-SNPsaveVCF 1\" MG-TK option. Aborting\n";}
		#$postcmd .= "#DEBUG\ncp $tmpOut.lz4 $ofasConsDir\n\n";
	#	$postcmd .= "zcat `ls $tmpOut.*.gz $sortedFileList`  |   $vcfcnsScr $ofasCons.depStat $minDepth  $minCallQual | $pigzBin -p $samcores -c >$ofasCons.gz ;\n\n"; #$refFA.fai 
	#}
	
	$postcmd .= "\necho \"Finished depthStat\"\n\n";
	$postcmd .= "rm -f $tmpOut*;\n";
	#$postcmd .= "$pigzBin -p $samcores $ofasCons;\n";

	$cmdAll .= "\n$postcmd\n" ;#always needs to run, to create seq stats in log file#if ($run2ctg != 0);


	if ($cmdAll ne ""){
		$cmdAll .= "touch $SNPstone\n" if ($hasPrimaryRds);# unless (-e $SNPstone);
		$cmdAll .= "touch $SNPsuppStone\n" if ($SNPsuppStone ne "");
	}
	
	#die "$run2ctg\n$cmdAll\n";
	if ($myParL && !$runLocalTmp && $cmdAll ne ""){
		if ( ($overwrite || !-e "$ofasCons")){
			my ($dep,$qcmd) = qsubSystem($qsubDirE."$cmdFTag.cacSNP.sh",$postcmd,1,$memReqGB."G","Cons$x",join(";",@allDeps2),"",$immediateSubm,[],$QSBoptHR);
			$submissionCommands .= $qcmd unless ($immediateSubm);
			$rdep =$dep;
		}
	}
	#die "$cmdAll\n $SNPIHR->{genefna}\n";
	
	if ($runLocalTmp && $cmdAll ne ""){#qsub all together now
		#this is the new way of doing this
		my $tmpS = $QSBoptHR->{tmpSpace};
		$QSBoptHR->{tmpSpace} = 3*$memReqGB; #in GB
		my ($dep,$qcmd) = qsubSystem($qsubDirE."$cmdFTag.oSNPc.sh",$cmdAll,int($actualCores*1.1),$memReqGB."G","Cons$x",join(";",$jdep,@allDeps2),"",$immediateSubm,[],$QSBoptHR);
		$submissionCommands .= $qcmd unless ($immediateSubm);
		$rdep =$dep;
		$QSBoptHR->{tmpSpace} = $tmpS;
	}
	#$SNPIHR->{intermedVCF} = $oVcfCons;
	$SNPIHR->{cleanCmd} = $cleanCmd;
	#die;
	return wantarray ? ($rdep,$submissionCommands) : $rdep;
}






sub SVcall_vcf{
	my ($SNPIHR) = @_;
	my $mode = $SNPIHR->{callSVs};
	if ($mode ==0 ){return ("");}
	
	my $SVcallerFlag = "";
	if ($mode == 1){	$SVcallerFlag = "DL"; #delly
	}elsif ($mode == 2){	$SVcallerFlag = "GY"; # gridss
	}else {die"Invalid callSVs option: $mode\n";}
	my $bcftBin = getProgPaths("bcftools");

	
	#my $QSBoptHR = $SNPIHR->{QSHR};
	#my $x = $SNPIHR->{JNUM};
	my $jdep = ""; $jdep = $SNPIHR->{jdeps} if (exists($SNPIHR->{jdeps}));
	my $tmpdir = $SNPIHR->{nodeTmpD};
	#my $smplNm = $SNPIHR->{smpl};
	my $refFA = $SNPIHR->{assembly};
	#my $qsubDirE = $SNPIHR->{qsubDir};
	#my $scrDir = $SNPIHR->{scratch};
	#my $bamcram = $SNPIHR->{bamcram};
	#my $splitFAsize = $SNPIHR->{bpSplit};
	#my $overwrite = $SNPIHR->{overwrite};
	#my $runLocalTmp = $SNPIHR->{runLocal};
	#my $SVout= $SNPIHR->{vcfSVfile};
	#my $cmdFTag = $SNPIHR->{cmdFileTag};
	my $maxSNPcores= $SNPIHR->{maxCores};
	$maxSNPcores=1 if ($mode ==1); #currently limiting SVcalls to 1
	my $actualCores = $maxSNPcores;
	my $samCores = $maxSNPcores;
	#infer outdir
	my $outD = $SNPIHR->{vcfSVfile};$outD =~ s/\/[^\/]+$/\//;

#path to tmp files
	my $bamTmp = "$tmpdir/$SNPIHR->{smpl}-smd.bam";
	my $tmpVCF = "$tmpdir/SV.$SNPIHR->{smpl}.bcf";
	my $bamTmpS = "";
	my $tmpVCFS = ""; 
	my @tar = @{$SNPIHR->{MAR}}[0]; 
#	$tar[0] = ${$SNPIHR->{MAR}}[0]; #$preTar;
	my @tarS = ();
	if ($SNPIHR->{STOconSNPsupp} ne "" && exists ($SNPIHR->{MARsupp} )){
		$bamTmpS = "$tmpdir/$SNPIHR->{smpl}-smd.suppl.bam";
		$tmpVCFS = "$tmpdir/SV.sup-$SNPIHR->{smpl}.bcf";
		@tarS = @{$SNPIHR->{MARsupp}}[0];
	}
	

	my $smtBin = getProgPaths("samtools");


	my $xtra = "";$xtra .= "echo \"Preparing data\"\n";
	$xtra .= "mkdir -p $tmpdir;\n";#$SNPIHR->{scratch};\n";
	$xtra .= "$smtBin faidx $refFA;\n" unless (-e "$refFA.fai");
	$xtra .= "echo \"Creating c/bams indexes primary reads\"\n";
	
	if ($SNPIHR->{bamcram} eq "cram"){ #create index for bam/cram
		#$xtra .= "if [ ! -e $tar[0].crai ] || [ ! -s $tar[0].crai ]; then rm -f $tar[0].crai; $smtBin index -@ $samcores  $tar[0]; fi\n";
		$xtra .= "#uncramming already stored results..\n" . cram2bsam("$tar[0]",$refFA,$bamTmp,1,$samCores) ."\n" ;
		#create index..
		$xtra .= "$smtBin index -@ $samCores  $bamTmp;\n";
		if (@tarS ){
			##also consider creating bams for suppl mappings
			$xtra .= "#uncramming supplemental mappings ..\n" . cram2bsam("$tarS[0]",$refFA,$bamTmpS,1,$samCores) ."\n" ;
			$xtra .= "$smtBin index -@ $samCores  $bamTmpS;\n";
		}

	}  else {
		$xtra .= "ln -s $tar[0] $bamTmp;\n";
		$xtra .= "$smtBin index -@ $samCores  $bamTmp;\n";
		if (@tarS ){
			$xtra .= "ln -s $tarS[0] $bamTmpS;\n";
			$xtra .= "$smtBin index -@ $samCores  $bamTmpS;\n";
		}
	}
	
	
		

	my $cmd = "";
	if ($mode == 1){ #delly..
		my $dellyBin = getProgPaths("delly");
		#delly call -g hg38.fa input.bam > delly.vcf
		my $dmode = "call"; $dmode = "lr" if ($SNPIHR->{SeqTech} eq "PB" || $SNPIHR->{SeqTech} eq "ONT");
		$cmd .= "echo \"main delly2 call\";\n$dellyBin $dmode -g $refFA $bamTmp | $bcftBin view -O b -o $tmpVCF -\n";#> $tmpVCF\n ";
		
		if (@tarS){#suppl mappings..
			$dmode = "call"; $dmode = "lr" if ($SNPIHR->{SeqTechSuppl} eq "PB" || $SNPIHR->{SeqTechSuppl} eq "ONT");
			$cmd .= "echo \"supplemental delly call\";\n$dellyBin $dmode -g $refFA $bamTmpS | $bcftBin view -O b -o $tmpVCFS -\n ";
		}

	} else{ #gridss..
		my $gridssBin = getProgPaths("gridss");
		die "Gridss not implemented yet";
	}


	#cleanup..
	$cmd .= "if [ ! -d $outD ] ; then mkdir -p $outD; fi\n"; #create final outdir..
	$cmd .= "\n#Finalize by moving to final location\nmv $tmpVCF $SNPIHR->{vcfSVfile}\n";
	if (@tarS){
		$cmd .= "mv $tmpVCFS $SNPIHR->{vcfSVfileS}\n";
	}
	
	$cmd .= "#cleanup..\nrm -rf $tmpdir\n";
	
	my $cmdAll = $xtra ."\n\n".$cmd;
	#print $cmdAll;

	my $immediateSubm = exists($SNPIHR->{immediateSubm}) ? $SNPIHR->{immediateSubm} : 1;
	my ($dep,$qcmd) = qsubSystem($SNPIHR->{qsubDir} . "$SNPIHR->{cmdFileTag}.SV.sh",$cmdAll,int($actualCores),"20G","SV$SNPIHR->{JNUM}",$jdep,"",$immediateSubm,[],$SNPIHR->{QSHR});
	
	return wantarray ? ($dep, $immediateSubm ? "" : $qcmd) : $dep;
}



