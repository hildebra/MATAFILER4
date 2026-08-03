package Mods::SNP;
use warnings;
use strict;

use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::GenoMetaAss qw( gzipopen systemW filsizeMB readFasta readGFF writeFasta reverse_complement_IUPAC fileGZs fileGZe );
use Mods::Subm qw(qsubSystem emptyQsubOpt qsubSystem2);
use Mods::TamocFunc qw (cram2bsam);
use Mods::SampleCompletion qw(invalidate_sample_completion);
use File::Glob qw(bsd_glob);
use File::Path qw(make_path);


use Exporter qw(import);
our @EXPORT_OK = qw(SNPconsensus_vcf SVcall_vcf estimateConsensusCores);

sub estimateConsensusCores {
	my ($inputSizeMB, $maxCores) = @_;
	$inputSizeMB = 0 unless defined($inputSizeMB) && $inputSizeMB > 0;
	$maxCores = int($maxCores || 1);
	$maxCores = 1 if $maxCores < 1;

	# Tiered scaling gives modest mappings useful parallelism early, then adds
	# workers more gradually for large metagenomes. The final value remains
	# bounded by the user-configured SNP core ceiling.
	my @sizeThresholdsMB = (300, 600, 1024, 2048, 4096, 6144, 8192, 9216, 10240);
	my $cores = 1;
	$cores++ for grep { $inputSizeMB >= $_ } @sizeThresholdsMB;
	$cores = $maxCores if $cores > $maxCores;
	return $cores;
}

sub regionsFromFAI($){
	my ($inF) = @_;
	my @regions;
	open my $faiFH, '<', $inF or die "can't open fai $inF: $!\n";
	while (my $line = <$faiFH>){
		chomp $line;
		my @fields = split /\t/, $line;
		die "invalid FAI record in $inF: $line\n"
			unless @fields >= 2 && length($fields[0])
				&& $fields[1] =~ /^\d+$/ && $fields[1] > 0;
		push @regions, "$fields[0]:1-$fields[1]";
	}
	close $faiFH or die "can't close fai $inF: $!\n";
	return @regions;
}

sub getRegionsBamDepth{
	my ($depthPC, $totalSpl, $maxSNPcores) = @_;
	my (%lengthFor, %depthFor);
	my @contigs;
	my ($totalDepthBases, $totalLength, $allLengthsKnown) = (0, 0, 1);
	my ($depthFH, $status) = gzipopen($depthPC, "contig depth file", 0);
	return ([], []) unless $status && defined $depthFH;
	while (my $line = <$depthFH>){
		chomp $line;
		next if $line =~ /^\s*$/;
		my @fields = split /\t/, $line;
		die "invalid contig-depth record in $depthPC: $line\n"
			unless @fields >= 2 && length($fields[0])
				&& $fields[1] =~ /^\d+(?:\.\d+)?(?:[eE][+-]?\d+)?$/;
		my ($contig, $depth) = @fields[0, 1];
		die "duplicate contig '$contig' in $depthPC\n" if exists $lengthFor{$contig};
		my $length;
		if ($contig =~ /__C\d+_L=(\d+)=/) {
			$length = $1;
		} else {
			$length = 1;
			$allLengthsKnown = 0;
		}
		push @contigs, $contig;
		$lengthFor{$contig} = $length;
		$depthFor{$contig} = 0 + $depth;
		$totalLength += $length;
		$totalDepthBases += $length * $depth;
	}
	close $depthFH or die "can't close contig depth file $depthPC: $!\n";
	# A guessed length can omit or duplicate bases. Let the caller fall back to
	# the authoritative FASTA index when contig names do not encode lengths.
	return ([], []) unless @contigs && $totalLength > 0 && $allLengthsKnown;

	$totalSpl = int($totalSpl || 1);
	$maxSNPcores = int($maxSNPcores || 1);
	$totalSpl = $maxSNPcores if $totalSpl > $maxSNPcores;
	$totalSpl = 1 if $totalSpl < 1 || $totalDepthBases < 1e6;
	if ($totalDepthBases < 5e6) {
		my $cap = int($maxSNPcores / 3) || 1;
		$totalSpl = $cap if $totalSpl > $cap;
	} elsif ($totalDepthBases < 20e6) {
		my $cap = int($maxSNPcores / 2) || 1;
		$totalSpl = $cap if $totalSpl > $cap;
	}

	my $averageDepth = $totalDepthBases / $totalLength;
	my $targetCost = $totalDepthBases > 0 ? $totalDepthBases / $totalSpl : 0;
	my @regions;
	my ($regionIndex, $regionCost) = (0, 0);
	for my $contig (@contigs) {
		my $length = $lengthFor{$contig};
		my $weight = $depthFor{$contig};
		$weight *= 1.2 if $averageDepth > 0 && $depthFor{$contig} > $averageDepth * 0.85;
		$weight *= 0.8 if $averageDepth > 0 && $depthFor{$contig} < $averageDepth * 0.1;
		my $start = 0;
		while ($start < $length) {
			my $take = $length - $start;
			if ($weight > 0 && $regionIndex < $totalSpl - 1) {
				my $remainingCost = $targetCost - $regionCost;
				if ($remainingCost <= 0) {
					$regionIndex++;
					$regionCost = 0;
					next;
				}
				my $capacity = int($remainingCost / $weight);
				$capacity = 1 if $capacity < 1;
				$take = $capacity if $capacity < $take;
			}
			my $stop = $start + $take;
			$regions[$regionIndex] .= "$contig\t$start\t$stop\n";
			$regionCost += $take * $weight;
			$start = $stop;
			if ($regionIndex < $totalSpl - 1 && $targetCost > 0 && $regionCost >= $targetCost) {
				$regionIndex++;
				$regionCost = 0;
			}
		}
	}
	@regions = grep { defined($_) && length($_) } @regions;
	return (\@regions, \@contigs);
}

sub getRegionsBam{
	my ($splitFAsizeL, $refFA, $tmpD) = @_;
	die "reference FASTA is missing or empty: $refFA\n" unless -s $refFA;
	die "SNP region size must be positive\n"
		unless defined($splitFAsizeL) && $splitFAsizeL > 0;
	unless (-s "$refFA.fai") {
		my $smtBin = getProgPaths("samtools");
		systemW "$smtBin faidx $refFA";
	}
	my @indexedRegions = regionsFromFAI("$refFA.fai");
	my (@regions, @contigOrder);
	my ($regionBases, $regionIndex) = (0, 0);
	for my $indexedRegion (@indexedRegions) {
		my ($contig, $length) = $indexedRegion =~ /^(.+):1-(\d+)$/
			or die "invalid indexed region '$indexedRegion' for $refFA\n";
		push @contigOrder, $contig;
		my $start = 0;
		while ($start < $length) {
			my $capacity = $splitFAsizeL - $regionBases;
			my $remaining = $length - $start;
			my $take = $remaining < $capacity ? $remaining : $capacity;
			my $stop = $start + $take;
			$regions[$regionIndex] .= "$contig\t$start\t$stop\n";
			$regionBases += $take;
			$start = $stop;
			if ($regionBases == $splitFAsizeL) {
				$regionIndex++;
				$regionBases = 0;
			}
		}
	}
	return (\@regions, \@contigOrder);
}


sub pileupcall{
	
	my ($tarR,$tag,$SNPIHR,$QSBoptHR,$scrDir,$tmpOut,$myParL,$crAR, $varsOnly) = @_;
	die "SNP pileup mapping is not defined\n"
		unless ref($tarR) eq 'ARRAY' && @{$tarR} && defined($tarR->[0]) && length($tarR->[0]);
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

	#basic caller options..
	my $minBQ=30; my $minMQ=30.0;
	#freebayes std options
	my $frAllOpts= "-u -i -m $minMQ -q $minBQ -C 1 -F 0.1 -k -X --pooled-continuous --report-monomorphic  --min-repeat-entropy 1 --use-best-n-alleles 2 -G 1 ";
	#bcftools options #-q = map qual -Q = base qual
	# Keep the pipe in uncompressed BCF and request only annotations accepted by
	# current bcftools mpileup. Bias INFO fields are emitted by the caller when
	# applicable and are not valid values for mpileup's --annotate option.
	my $bcfAllOpts = " -Ou --min-BQ $minBQ -d 12000 --threads $threads --min-MQ $minMQ -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP";
	#--skip-indels #--count-orphans -a DP,AD,ADF,ADR,SP
	my $bcftCallOpts = " --ploidy 1 -c -M --threads 0 -Ou "; #--multiallelic-caller  -> replaced with -c (consensus caller)
	$bcftCallOpts .= "--variants-only " if ($varsOnly);
	#$bcftCallOpts .= "-O z " ;# if ($normalizeIndel);
	if ($tag eq "sup-"){
		my $technology = $SNPIHR->{SeqTechSuppl} // "";
		if ($technology eq "ONT"){$bcfAllOpts.=" -X ont ";
		} elsif ($technology eq "PB"){$bcfAllOpts.=" -X pacbio-ccs ";
		} else{$bcfAllOpts.=" -X illumina ";}
	} else {
		my $technology = $SNPIHR->{SeqTech} // "";
		if ($technology eq "ONT"){$bcfAllOpts.=" -X ont ";
		} elsif ($technology eq "PB"){$bcfAllOpts.=" -X pacbio-ccs ";
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
	my @allDeps2; my @chunkFiles;
	#implement in parallel as too slow in single core mode :/
	my @oldChunks = (
		bsd_glob("$tmpOut.$tag*.vcf.gz"),
		bsd_glob("$tmpOut.$tag*.vcf.gz.csi"),
		bsd_glob("$tmpOut.$tag*.vcf.gz.tbi"),
	);
	unlink @oldChunks if $overwrite && @oldChunks;
	make_path($qsubDirE) unless -d $qsubDirE;
	my $bedJobs = 0;
	for (my $i=0;$i<@curReg;$i++){ #go over regions in bed file, submit a job for each "region"
		#$tar[0] = bam file;  $bedF = bedfile with regions
		next if (!$run2ctg);
		my $cmd2 = $cmd ;
		my $bedF = $qsubDirE."$smplNm.${tag}$i.bed";
		my $chunkFile = "$tmpOut.$tag$i.vcf.gz";
		push @chunkFiles, $chunkFile;
		next if (-s $chunkFile && (-s "$chunkFile.csi" || -s "$chunkFile.tbi") && !$overwrite);
		if ($myParL){
			if (-s $chunkFile && !$overwrite) {
				# A prior run may have completed the expensive call but failed before
				# concat. Repair that restart state by indexing the existing BGZF chunk.
				$cmd2 = "$bcftBin index -f $chunkFile && test -s $chunkFile.csi && rm -f $bedF $locXtrCmd\n";
			} elsif ($useFB){
				$cmd2 .= " -t $bedF $tarR->[0] | $bcftBin view -Oz -o $chunkFile - && test -s $chunkFile && $bcftBin index -f $chunkFile && test -s $chunkFile.csi && rm -f $bedF $locXtrCmd\n";
			}else{
				$cmd2 .= " -R $bedF $tarR->[0] | $bcftBin call $bcftCallOpts | $bcftBin view -Oz -o $chunkFile - && test -s $chunkFile && $bcftBin index -f $chunkFile && test -s $chunkFile.csi && rm -f $bedF $locXtrCmd\n";
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
	$cmdAll2 .= "\necho \"Finished mpileup $tag\"\n\n";	
	#$cmdAll .= "rm -f $tarR->[0];\n" if ($run2ctg && $tarR->[0] ne "");
	if ($bedJobs ==0){$cmdAll2="";}
	return (\@allDeps2, \@chunkFiles, $cmdAll2);
}

sub _coverage_file_for_mapping {
	my ($mapping, $label, $allowPendingInputs) = @_;
	my $coverage = $mapping;
	die "$label mapping has an unsupported suffix: $mapping\n"
		unless $coverage =~ s/\.(?:cram|bam)$/\.bam.coverage/;
	return "$coverage.gz" if -s "$coverage.gz";
	return $coverage if -s $coverage;
	# Canonical assembly mappings publish compressed coverage. When the mapping
	# job is an afterok dependency, use that known future path rather than
	# inspecting a file which cannot exist yet.
	return "$coverage.gz" if $allowPendingInputs;
	die "$label coverage is missing for $mapping (expected $coverage or $coverage.gz)\n";
}


sub SNPconsensus_vcf{
	my ($SNPIHR)  = @_;
	die "SNP consensus options must be a hash reference\n" unless ref($SNPIHR) eq 'HASH';
	for my $required (qw(assembly nodeTmpD smpl qsubDir scratch bamcram bpSplit maxCores SNPcaller ofas QSHR)) {
		die "SNP consensus option '$required' is missing\n"
			unless exists($SNPIHR->{$required}) && defined($SNPIHR->{$required}) && length($SNPIHR->{$required});
	}
	invalidate_sample_completion($SNPIHR->{sampleRoot}) if $SNPIHR->{sampleRoot};
	my $smtBin = getProgPaths("samtools");
	my $bcftBin = getProgPaths("bcftools");
	my $vcf2fnaBin = getProgPaths("vcf2fna");
	my $regionPlanner = getProgPaths("consVCF_region_planner");
	my $pigzBin = getProgPaths("pigz");
	
	my $memPJob =0;
	#get parameteres
	my $samcores = 12;
	#my %SNPinfo = %{$SNPIHR};
	my $QSBoptHR = $SNPIHR->{QSHR};
	my $immediateSubm = exists($SNPIHR->{immediateSubm}) ? $SNPIHR->{immediateSubm} : 1;
	my $submissionCommands = "";
	my $x = $SNPIHR->{JNUM};
	my $jdep = ""; $jdep = $SNPIHR->{jdeps} if (exists($SNPIHR->{jdeps}));
	my $allowPendingInputs = $SNPIHR->{allowPendingInputs} ? 1 : 0;
	my $tmpdir = $SNPIHR->{nodeTmpD};
	my $smplNm = $SNPIHR->{smpl};
	my $refFA = $SNPIHR->{assembly};
	my $qsubDirE = $SNPIHR->{qsubDir};
	my $scrDir = $SNPIHR->{scratch};
	my $bamcram = $SNPIHR->{bamcram};
	my $splitFAsize = $SNPIHR->{bpSplit};
	my $overwrite = $SNPIHR->{overwrite} ? 1 : 0;
	my $runLocalTmp = $SNPIHR->{runLocal} ? 1 : 0;
	die "SNP consensus requires runLocal mode so preparation, chunks, and finalization share one checked allocation\n"
		unless $runLocalTmp;
	die "pending SNP inputs require scheduler dependencies for immediate submission\n"
		if ($allowPendingInputs && $immediateSubm && $jdep !~ /\S/);
	my $maxSNPcores= int($SNPIHR->{maxCores});
	die "maxCores must be positive\n" unless $maxSNPcores > 0;
	# A sample may have primary reads while this invocation requests only the
	# supplementary consensus. Keep physical read availability separate from
	# the requested SNP product so an already-complete primary call is not rerun.
	my $hasPrimaryRds = exists($SNPIHR->{callConsSNP})
		? ($SNPIHR->{callConsSNP} ? 1 : 0)
		: (exists($SNPIHR->{hasPrimaryRds}) ? ($SNPIHR->{hasPrimaryRds} ? 1 : 0) : 1);
	my $normalizeIndel = 1; $normalizeIndel = $SNPIHR->{normIndels} if (exists($SNPIHR->{normIndels})) ;
	my $onlyNormalize = 0;

	$memPJob = $SNPIHR->{memPJob} || 0;
	my $SNPstone = $SNPIHR->{STOconSNP} // "";
	my $SNPsuppStone = $SNPIHR->{STOconSNPsupp} // "";
	my $minDepth = 0;
	$minDepth =$SNPIHR->{minDepth}  if (exists($SNPIHR->{minDepth} ));
	my $minCallQual = $SNPIHR->{minCallQual} // 20;
	#my $SNPstone = $ofasConsDir."SNP.cons.stone";
	#my $memReq = "20G";
	#my $memReq = $SNPIHR->{memReq};
	my $vcfFile = ""; $vcfFile = $SNPIHR->{vcfFile} if (exists ($SNPIHR->{vcfFile}));
	my $vcfFileS = ""; $vcfFileS = $SNPIHR->{vcfFileSupp} if (exists ($SNPIHR->{vcfFileSupp}));
	my $cmdFTag = $SNPIHR->{cmdFileTag} // "ConsSNP";
	my $firstInSample = 0;$firstInSample = $SNPIHR->{firstInSample} if (exists($SNPIHR->{firstInSample}));
	my $useFB = uc($SNPIHR->{SNPcaller}) eq "MPI" ? 0 : 1;
	my $actualCores  = $maxSNPcores;
	my $consensusInputMB = 0;
	for my $mappingSet (
		[$hasPrimaryRds, $SNPIHR->{MAR}],
		[length($SNPsuppStone), $SNPIHR->{MARsupp}],
	) {
		next unless $mappingSet->[0] && ref($mappingSet->[1]) eq 'ARRAY';
		for my $mapping (@{$mappingSet->[1]}) {
			next unless defined($mapping) && length($mapping) && -s $mapping;
			my $sizeMB = filsizeMB($mapping);
			$consensusInputMB = $sizeMB if $sizeMB > $consensusInputMB;
		}
	}
	if ($consensusInputMB <= 0 && ($SNPIHR->{inputSizeMB} || 0) > 0) {
		$consensusInputMB = $SNPIHR->{inputSizeMB};
	}
	

	my $saveVCF = exists($SNPIHR->{saveVCF}) ? ($SNPIHR->{saveVCF} ? 1 : 0) : 1;
	if ($vcfFile eq ""){
		$saveVCF=0;
		$vcfFile = "$scrDir/$smplNm.fin.vcf.gz"; #switch back to vcf for vcf2fna
		$vcfFileS = "$scrDir/$smplNm.fin-sup.vcf.gz";
	}
	$vcfFile .= ".gz" unless $vcfFile =~ /\.gz$/;
	$vcfFileS .= ".gz" if length($vcfFileS) && $vcfFileS !~ /\.gz$/;
	if ($normalizeIndel && length($SNPstone)){
		my $SNPstone2 = $SNPstone; $SNPstone2 =~ s/\.norm\.stone/\.stone/;
		if (! -e $SNPstone && -e $SNPstone2){#can just do the normalization..
			$onlyNormalize = 1;
			unlink $SNPstone2 or die "can't remove obsolete SNP marker $SNPstone2: $!\n";
			$maxSNPcores=1;
		}
		#die "$onlyNormalize $SNPstone2 $SNPstone\n";
	}	
	


	
	#key change for C++ program vcf2fasta
	my $reportVCFonly=1; #was 0 before


	my $ofasCons = $SNPIHR->{ofas};
	my ($ofasConsDir) = $ofasCons =~ m{^(.*)/[^/]+$};
	die "SNP consensus output must include a directory: $ofasCons\n" unless defined($ofasConsDir) && length($ofasConsDir);
	$ofasConsDir .= "/";
	#$ofasCons .= ".gz" unless ($ofasCons =~ m/\.gz$/); #don't change, needed without..
	my $supportRequested = length($SNPsuppStone) ? 1 : 0;
	my $primaryVcfReady = !$hasPrimaryRds || fileGZe($vcfFile);
	my $supportVcfReady = !$supportRequested || fileGZe($vcfFileS);
	my $runPrimary = $hasPrimaryRds && ($overwrite || !$primaryVcfReady);
	my $runSupport = $supportRequested && ($overwrite || !$supportVcfReady);
	my $run2ctg = $runPrimary || $runSupport;
	if ($overwrite) {
		my @oldScratchVcfs = grep { -f $_ } bsd_glob("$scrDir/$smplNm*.vcf.gz");
		unlink @oldScratchVcfs if @oldScratchVcfs;
	}
	#first all important regions on finalDir
	my @curReg = ("1");
	my $myParL=0;
	if ($splitFAsize>0){$myParL=1;}
	if ($myParL && $run2ctg) {
		my $configuredJobs = int($SNPIHR->{split_jobs} || 0);
		my $runtimeJobs = $configuredJobs > 0
			? $configuredJobs
			: estimateConsensusCores($consensusInputMB, $maxSNPcores);
		$runtimeJobs = $maxSNPcores if ($runtimeJobs > $maxSNPcores);
		$runtimeJobs = 1 if ($runtimeJobs < 1);
		@curReg = (('runtime') x $runtimeJobs);
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
	$xtra .= "echo \"Consensus SNP allocation: $actualCores cores (input estimate: "
		.int($consensusInputMB + 0.5)." MB)\"\n";
	$xtra .= "mkdir -p $scrDir;\n";
	#$xtra .= "exit\n"; #DEBUG
	#$xtra .= "cp $refFA $refFA.fai $scrDir;\n";$refFA =~ m/\/([^\/]+$)/;$refFA = "$scrDir/$1";
	#my $preTar = 
	
	#my @tar = ("");$tar[0] = ${$SNPIHR->{MAR}}[0]; #$preTar;
	
	my @tar = $hasPrimaryRds && ref($SNPIHR->{MAR}) eq 'ARRAY' ? ($SNPIHR->{MAR}->[0]) : ();
	my $cmdAll = "";my @allDeps2; my (@primaryChunks, @supportChunks);
	my $tmpOut = "$scrDir/$smplNm.cons.vcf";my $depthFile ="";

	#supplementary mappings?
	my @tarS = $supportRequested && ref($SNPIHR->{MARsupp}) eq 'ARRAY' ? ($SNPIHR->{MARsupp}->[0]) : ();
	die "primary SNP mapping is missing\n"
		if $hasPrimaryRds && (!@tar || !defined($tar[0])
			|| (!$allowPendingInputs && !-s $tar[0]));
	die "supplementary SNP mapping is missing\n"
		if $supportRequested && (!@tarS || !defined($tarS[0])
			|| (!$allowPendingInputs && !-s $tarS[0]));
	die "SNP reference is missing or empty: $refFA\n"
		unless $allowPendingInputs || -s $refFA;
	if ($allowPendingInputs) {
		$xtra .= "test -s $refFA\n";
		$xtra .= "test -s $tar[0]\n" if $hasPrimaryRds;
		$xtra .= "test -s $tarS[0]\n" if $supportRequested;
	}
	$xtra .= "$smtBin faidx $refFA;\n" unless (-s "$refFA.fai");
	my $tmpOut2 = "$scrDir/$smplNm.X.cons.vcf";my $depthFileS  = "";
	my ($primaryRegionCmd, $supportRegionCmd) = ("", "");
	if ($myParL && $run2ctg) {
		my $runtimeJobs = scalar(@curReg);
		my $bedPrefix = "$qsubDirE/$smplNm.";
		$xtra .= "rm -f ${bedPrefix}*.bed\n";
		if ($runPrimary) {
			my $depthArg = defined($SNPIHR->{depthF}) && length($SNPIHR->{depthF})
				? " --depth $SNPIHR->{depthF}" : "";
			$primaryRegionCmd = "$regionPlanner --fai $refFA.fai --mapping $tar[0]$depthArg --jobs $runtimeJobs --output-prefix $bedPrefix --samtools $smtBin --pigz $pigzBin\n";
		}
		if ($runSupport) {
			$supportRegionCmd = "$regionPlanner --fai $refFA.fai --mapping $tarS[0] --jobs $runtimeJobs --output-prefix ${bedPrefix}sup- --samtools $smtBin --pigz $pigzBin\n";
		}
	}
	
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
		$depthFile = _coverage_file_for_mapping($tar[0], "primary SNP", $allowPendingInputs);
		$xtra .= "test -s $depthFile\n" if $allowPendingInputs;
	}
	if ($runPrimary){
		$xtra .= "echo \"Creating c/bams indexes primary reads\"\n";
		if ($bamcram eq "cram"){ #create index for bam/cram
			$xtra .= "if [ ! -e $tar[0].crai ] || [ ! -s $tar[0].crai ]; then rm -f $tar[0].crai; $smtBin index -@ $samcores  $tar[0]; fi\n";
		} else {
			$xtra .= "if [ ! -e $tar[0].bai ] || [ ! -s $tar[0].bai ]; then rm -f $tar[0].bai; $smtBin index -@ $samcores  $tar[0]; fi\n";
		}
		$xtra .= $primaryRegionCmd;
		
		$SNPIHR->{run2ctg} = 1;$SNPIHR->{rdep} = $rdep;

		#$SNPIHR->{assembly} = $refFA;
		$cmdAll .= $xtra if (!$onlyNormalize);
		my ($dAR,$cAR,$pilecmd) =  pileupcall(\@tar,"",$SNPIHR,$QSBoptHR,$scrDir,$tmpOut,$myParL,\@curReg,$reportVCFonly);
		@allDeps2 = @{$dAR}; @primaryChunks = @{$cAR};
		$cmdAll .= $pilecmd if (!$onlyNormalize);
	}
	
	if (@tarS){
		$depthFileS = _coverage_file_for_mapping($tarS[0], "supplementary SNP", $allowPendingInputs);
	}
	if ($runSupport){ #supplementary reads SNP call
		my $xtra2 = "echo \"Creating c/bams indexes supplemental reads\"\n";
		$xtra2 .= "test -s $depthFileS\n" if $allowPendingInputs;
		$cmdAll .= $xtra if !$runPrimary && !$onlyNormalize;
		$SNPIHR->{run2ctg} = 1;
		#die "$depthFileS\n";
		if ($bamcram eq "cram"){ #create index for bam/cram
			$xtra2 .= "if [ ! -e $tarS[0].crai ] || [ ! -s $tarS[0].crai ]; then rm -f $tarS[0].crai; $smtBin index -@ $samcores  $tarS[0]; fi\n";
		} else {
			$xtra2 .= "if [ ! -e $tarS[0].bai ] || [ ! -s $tarS[0].bai ]; then rm -f $tarS[0].bai; $smtBin index -@ $samcores  $tarS[0]; fi\n";
		}
		$xtra2 .= $supportRegionCmd;
		my ($dAR,$cAR,$pilecmd) =  pileupcall(\@tarS,"sup-",$SNPIHR,$QSBoptHR,$scrDir,$tmpOut2,$myParL,\@curReg,$reportVCFonly);
		$cmdAll .= $xtra2.$pilecmd if (!$onlyNormalize);
		push(@allDeps2, @{$dAR}); @supportChunks = @{$cAR};
	}
	
	# Every successful chunk removes its own BED file.
	$cmdAll .= "if ls $qsubDirE/$smplNm.*.bed 1> /dev/null 2>&1 ;then echo \"Bed files still present, probably incorrect run\"; exit 33; else echo \"bed files deleted, looks good\"; fi\n\n" if (@curReg);
		

	
	#from here on: merge XX vcf's into one
	my $sortCmd = "";
	if ($myParL && $cmdAll ne ""){
		$sortCmd .= "mkdir -p $ofasConsDir;\n";
		if ($runPrimary){
			die "no primary VCF chunks were planned\n" unless @primaryChunks;
			$sortCmd .= "$bcftBin concat -a -Oz -o $vcfFile ".join(" ", @primaryChunks)."\n";
			my @primaryChunkArtifacts = map { ($_, "$_.csi", "$_.tbi") } @primaryChunks;
			$sortCmd .= "test -s $vcfFile\nrm -f ".join(" ", @primaryChunkArtifacts)."\n";
		}
		if ($runSupport){
			die "no supplementary VCF chunks were planned\n" unless @supportChunks;
			$sortCmd .= "$bcftBin concat -a -Oz -o $vcfFileS ".join(" ", @supportChunks)."\n";
			my @supportChunkArtifacts = map { ($_, "$_.csi", "$_.tbi") } @supportChunks;
			$sortCmd .= "test -s $vcfFileS\nrm -f ".join(" ", @supportChunkArtifacts)."\n";
		}
		$cmdAll .= $sortCmd;
	}
	
	my $bcfNormOpts = " -O z -f $refFA "; #-a -m '+both'
	my $cmd3 = "";
	if($normalizeIndel){#bcftools norm -f ref.fa in.vcf
		$cmd3.="\necho \"left-normalizing indels\"\n";
		if ($hasPrimaryRds){
			$cmd3 .= "$bcftBin index -f $vcfFile\n$bcftBin norm $bcfNormOpts -o $vcfFile.norm $vcfFile\n";
			$cmd3 .= "test -s $vcfFile.norm\nrm -f $vcfFile $vcfFile.csi; mv $vcfFile.norm $vcfFile;\n";
		}
		
		if ($supportRequested){
			$cmd3 .= "$bcftBin index -f $vcfFileS\n$bcftBin norm $bcfNormOpts -o $vcfFileS.norm $vcfFileS\n";
			$cmd3 .= "test -s $vcfFileS.norm\nrm -f $vcfFileS $vcfFileS.csi; mv $vcfFileS.norm $vcfFileS;\n";
		}
		
		$cmd3 .= "\necho \"Done normalizing\"\n"; #wait \$(jobs -p);\n
		#die "$cmd3";
		$cmdAll .= $cmd3;
	}
	
	my $postcmd = "";
	my $vcf2fnaOpt = "";
	my $createFastas = $SNPIHR->{createFastas} ? 1 : 0;
	my $vcf2fnaOuts = $createFastas ? "-oCtg $ofasCons.gz" : "";
	my $createGeneFastas = $createFastas && defined($SNPIHR->{genefna}) && length($SNPIHR->{genefna})
		&& defined($SNPIHR->{genefaa}) && length($SNPIHR->{genefaa});
	if ($createGeneFastas && !fileGZe($SNPIHR->{genefna})){
		$vcf2fnaOuts .= " -oGeneNT $SNPIHR->{genefna} -oGeneAA $SNPIHR->{genefaa} ";
	}
	my $vcf2fnaIns = "-ref $refFA ";
	if (defined($SNPIHR->{gffFile}) && length($SNPIHR->{gffFile})) {
		my $gffF = $SNPIHR->{gffFile};
		if (!-s $gffF && $gffF =~ /\.gz$/) {
			(my $plainGff = $gffF) =~ s/\.gz$//;
			$gffF = $plainGff if -s $plainGff;
		} elsif (!-s $gffF && -s "$gffF.gz") {
			$gffF .= ".gz";
		}
		my $gffAvailable = -s $gffF || ($allowPendingInputs && $createGeneFastas);
		$vcf2fnaIns .= "-gff $gffF " if $gffAvailable;
		die "SNP GFF is missing or empty: $gffF\n"
			if $createGeneFastas && !$gffAvailable;
		$postcmd .= "test -s $gffF\n" if $allowPendingInputs && $createGeneFastas;
	} elsif ($createGeneFastas) {
		die "SNP GFF is required when gene consensus FASTAs are requested\n";
	}
	
	if (!$hasPrimaryRds){ #only support available..
		my $tmpST = $SNPIHR->{SeqTechSuppl} // ""; if ($tmpST eq ""){$tmpST = "ill";}
		$vcf2fnaOpt = "-seqPlatform $tmpST -t 1 -minCallDepth $minDepth -minCallQual $minCallQual ";
		$vcf2fnaIns .= "-inVCF $vcfFileS -depthF $depthFileS ";
	} elsif ($SNPsuppStone eq "" ){#only primary reads available..
		my $tmpST = $SNPIHR->{SeqTech} // ""; if ($tmpST eq ""){$tmpST = "ill";}
		$vcf2fnaOpt = "-seqPlatform $tmpST -t 1 -minCallDepth $minDepth -minCallQual $minCallQual ";
		$vcf2fnaIns .= "-inVCF $vcfFile -depthF $depthFile ";
	} else {#and for two vcfs..
		$vcf2fnaOpt = "-seqPlatform $SNPIHR->{SeqTech},$SNPIHR->{SeqTechSuppl} -t 1 -minCallDepth $minDepth,$minDepth -minCallQual $minCallQual ";
		$vcf2fnaIns .= "-inVCF $vcfFile,$vcfFileS -depthF $depthFile,$depthFileS ";
	}
	if (!$createFastas){
		$postcmd.="\n##In case you want to create consensus fastas, use (uncomment):\n##$vcf2fnaBin $vcf2fnaOpt $vcf2fnaIns $vcf2fnaOuts\n";
		$postcmd.="#create stats only of hypothetical consensus generations:\n$vcf2fnaBin $vcf2fnaOpt $vcf2fnaIns; \n\n"
	} else {
		$postcmd.="\n# Create consensus fastas\n";
		$postcmd.="$vcf2fnaBin $vcf2fnaOpt $vcf2fnaIns $vcf2fnaOuts;\n\n";
	}

	
	
	$postcmd .= "test -s $ofasCons.gz\n" if $createFastas;
	$postcmd .= "test -s $SNPIHR->{genefna}\ntest -s $SNPIHR->{genefaa}\n" if $createGeneFastas;
	$postcmd .= "rm -f $vcfFile.csi $vcfFileS.csi\n";
	$postcmd .= "rm -f $vcfFile $vcfFileS\n" if !$saveVCF;
	#} else {
	#	if ($SNPsuppStone ne "" ){die "support reads activated. combined SNP calling only works current with use of the \"-SNPsaveVCF 1\" MATAFILER option. Aborting\n";}
		#$postcmd .= "#DEBUG\ncp $tmpOut.lz4 $ofasConsDir\n\n";
	#	$postcmd .= "zcat `ls $tmpOut.*.gz $sortedFileList`  |   $vcfcnsScr $ofasCons.depStat $minDepth  $minCallQual | $pigzBin -p $samcores -c >$ofasCons.gz ;\n\n"; #$refFA.fai 
	#}
	
	$postcmd .= "\necho \"Finished depthStat\"\n\n";
	$postcmd .= "rm -f $tmpOut.*.vcf.gz $tmpOut2.*.vcf.gz\n";
	#$postcmd .= "$pigzBin -p $samcores $ofasCons;\n";

	$cmdAll .= "\n$postcmd\n" ;#always needs to run, to create seq stats in log file#if ($run2ctg != 0);


	if ($cmdAll ne ""){
		$cmdAll .= "test -s $vcfFile\n" if $hasPrimaryRds && $saveVCF;
		$cmdAll .= "test -s $vcfFileS\n" if $supportRequested && $saveVCF;
		$cmdAll .= "touch $SNPstone\n" if $hasPrimaryRds && length($SNPstone);
		$cmdAll .= "touch $SNPsuppStone\n" if $supportRequested && length($SNPsuppStone);
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
		my ($dep,$qcmd) = qsubSystem($qsubDirE."$cmdFTag.oSNPc.sh",$cmdAll,$actualCores,$memReqGB."G","Cons$x",join(";",$jdep,@allDeps2),"",$immediateSubm,[],$QSBoptHR);
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
	my $allowPendingInputs = $SNPIHR->{allowPendingInputs} ? 1 : 0;
	my $hasPrimary = exists($SNPIHR->{hasPrimaryRds}) ? ($SNPIHR->{hasPrimaryRds} ? 1 : 0) : 1;
	my $callPrimary = $hasPrimary && ($SNPIHR->{callSVs} || 0);
	my $callSupport = $SNPIHR->{callSVsSupp} || 0;
	my $mode = $callPrimary || $callSupport;
	if ($mode ==0 ){return ("");}
	die "primary and supplementary SV callers must use the same mode\n"
		if $callPrimary && $callSupport && $callPrimary != $callSupport;
	
	my $SVcallerFlag = "";
	if ($mode == 1){	$SVcallerFlag = "DL"; #delly
	}elsif ($mode == 2){	$SVcallerFlag = "GY"; # gridss
	}else {die"Invalid callSVs option: $mode\n";}
	invalidate_sample_completion($SNPIHR->{sampleRoot}) if $SNPIHR->{sampleRoot};
	my $bcftBin = getProgPaths("bcftools");

	
	#my $QSBoptHR = $SNPIHR->{QSHR};
	#my $x = $SNPIHR->{JNUM};
	my $jdep = ""; $jdep = $SNPIHR->{jdeps} if (exists($SNPIHR->{jdeps}));
	my $tmpdir = $SNPIHR->{nodeTmpD}."/SV.$SNPIHR->{smpl}";
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
	my $finalSV = $callPrimary ? $SNPIHR->{vcfSVfile} : $SNPIHR->{vcfSVfileS};
	die "SV output file is not defined\n" unless defined($finalSV) && length($finalSV);
	my $outD = $finalSV;$outD =~ s/\/[^\/]+$/\//;

#path to tmp files
	my $bamTmp = "$tmpdir/$SNPIHR->{smpl}-smd.bam";
	my $tmpVCF = "$tmpdir/SV.$SNPIHR->{smpl}.bcf";
	my $bamTmpS = "";
	my $tmpVCFS = ""; 
	my @tar = $callPrimary && ref($SNPIHR->{MAR}) eq 'ARRAY' ? ($SNPIHR->{MAR}->[0]) : ();
#	$tar[0] = ${$SNPIHR->{MAR}}[0]; #$preTar;
	my @tarS = ();
	if ($callSupport && exists ($SNPIHR->{MARsupp} )){
		$bamTmpS = "$tmpdir/$SNPIHR->{smpl}-smd.suppl.bam";
		$tmpVCFS = "$tmpdir/SV.sup-$SNPIHR->{smpl}.bcf";
		@tarS = @{$SNPIHR->{MARsupp}}[0];
	}
	die "primary SV mapping is missing\n"
		if $callPrimary && (!@tar || (!$allowPendingInputs && !-s $tar[0]));
	die "supplementary SV mapping is missing\n"
		if $callSupport && (!@tarS || (!$allowPendingInputs && !-s $tarS[0]));
	die "SV reference is missing or empty: $refFA\n"
		unless defined($refFA) && ($allowPendingInputs || -s $refFA);
	

	my $smtBin = getProgPaths("samtools");


	my $xtra = "";$xtra .= "echo \"Preparing data\"\n";
	$xtra .= "rm -rf $tmpdir\nmkdir -p $tmpdir;\n";
	if ($allowPendingInputs) {
		$xtra .= "test -s $refFA\n";
		$xtra .= "test -s $tar[0]\n" if $callPrimary;
		$xtra .= "test -s $tarS[0]\n" if $callSupport;
	}
	$xtra .= "$smtBin faidx $refFA;\n" unless (-s "$refFA.fai");
	$xtra .= "echo \"Creating c/bams indexes primary reads\"\n";
	
	if ($SNPIHR->{bamcram} eq "cram"){ #create index for bam/cram
		#$xtra .= "if [ ! -e $tar[0].crai ] || [ ! -s $tar[0].crai ]; then rm -f $tar[0].crai; $smtBin index -@ $samcores  $tar[0]; fi\n";
		if ($callPrimary) {
			$xtra .= "#uncramming already stored results..\n" . cram2bsam("$tar[0]",$refFA,$bamTmp,1,$samCores) ."\n" ;
			$xtra .= "$smtBin index -@ $samCores  $bamTmp;\n";
		}
		if (@tarS ){
			##also consider creating bams for suppl mappings
			$xtra .= "#uncramming supplemental mappings ..\n" . cram2bsam("$tarS[0]",$refFA,$bamTmpS,1,$samCores) ."\n" ;
			$xtra .= "$smtBin index -@ $samCores  $bamTmpS;\n";
		}

	}  else {
		if ($callPrimary) {
			$xtra .= "ln -sf $tar[0] $bamTmp;\n";
			$xtra .= "$smtBin index -@ $samCores  $bamTmp;\n";
		}
		if (@tarS ){
			$xtra .= "ln -s $tarS[0] $bamTmpS;\n";
			$xtra .= "$smtBin index -@ $samCores  $bamTmpS;\n";
		}
	}
	
	
		

	my $cmd = "";
	if ($mode == 1){ #delly..
		my $dellyBin = getProgPaths("delly");
		#delly call -g hg38.fa input.bam > delly.vcf
		my $primaryTechnology = $SNPIHR->{SeqTech} // "";
		my $dmode = "call"; $dmode = "lr" if ($primaryTechnology eq "PB" || $primaryTechnology eq "ONT");
		$cmd .= "echo \"main delly2 call\";\n$dellyBin $dmode -g $refFA $bamTmp | $bcftBin view -O b -o $tmpVCF -\ntest -s $tmpVCF\n" if $callPrimary;
		
		if (@tarS){#suppl mappings..
			my $supportTechnology = $SNPIHR->{SeqTechSuppl} // "";
			$dmode = "call"; $dmode = "lr" if ($supportTechnology eq "PB" || $supportTechnology eq "ONT");
			$cmd .= "echo \"supplemental delly call\";\n$dellyBin $dmode -g $refFA $bamTmpS | $bcftBin view -O b -o $tmpVCFS -\ntest -s $tmpVCFS\n";
		}

	} else{ #gridss..
		die "Gridss not implemented yet";
	}


	#cleanup..
	$cmd .= "if [ ! -d $outD ] ; then mkdir -p $outD; fi\n"; #create final outdir..
	$cmd .= "\n#Finalize by moving to final location\nmv $tmpVCF $SNPIHR->{vcfSVfile}\ntest -s $SNPIHR->{vcfSVfile}\n" if $callPrimary;
	if (@tarS){
		$cmd .= "mv $tmpVCFS $SNPIHR->{vcfSVfileS}\ntest -s $SNPIHR->{vcfSVfileS}\n";
	}
	
	$cmd .= "#cleanup..\nrm -rf $tmpdir\n";
	
	my $cmdAll = $xtra ."\n\n".$cmd;
	#print $cmdAll;

	my $immediateSubm = exists($SNPIHR->{immediateSubm}) ? $SNPIHR->{immediateSubm} : 1;
	die "pending SV inputs require scheduler dependencies for immediate submission\n"
		if ($allowPendingInputs && $immediateSubm && $jdep !~ /\S/);
	my ($dep,$qcmd) = qsubSystem($SNPIHR->{qsubDir} . "$SNPIHR->{cmdFileTag}.SV.sh",$cmdAll,int($actualCores),"20G","SV$SNPIHR->{JNUM}",$jdep,"",$immediateSubm,[],$SNPIHR->{QSHR});
	
	return wantarray ? ($dep, $immediateSubm ? "" : $qcmd) : $dep;
}
