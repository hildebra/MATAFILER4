package Mods::Binning;
# 2026-07 sparse-MGS hardening: accept header-only bin assignments as an empty result.
use Exporter qw(import);
our @EXPORT_OK = qw(
				runMetaBat runSemiBin  runMetaDecoder  runGenomeFace runSCGBinner
				runCheckM runCheckM2 MB2N50
				getBinSubdirName binningOutputsComplete emptyBinnerAssignmentCommand
				createBin2 createBinFAA createBinCtgs
				readMGS readMGSrev deNovo16S readMGSrevRed minQualFilter 
				filterMGS_CM MB2assigns MB2assignedBinIds calcLCAcompl readCMquals);

use warnings;
use strict;
use File::Basename;
use File::Path qw(make_path);
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::GenoMetaAss qw (systemW readFasta gzipopen getAssemblPath reverse_complement_IUPAC);
use Mods::TamocFunc qw (cram2bsam);



sub getBinSubdirName{
	my ($BinSel) = @_;
	my $BinnerName = "None";
	if ($BinSel == 1){
		$BinnerName = "MB2";
	} elsif ($BinSel == 2){#SemiBin
		$BinnerName = "SB";
	}elsif ($BinSel == 3){#MetaDecoder
		$BinnerName = "MD";
	}elsif ($BinSel == 4){#GenomeFace
		$BinnerName = "GF";
	}elsif ($BinSel == 5){#SCGBinner
		$BinnerName = "SC";
	}
	return $BinnerName;
}

sub binningOutputsComplete {
	my ($base, $useCheckM1, $useCheckM2) = @_;
	return 0 unless defined($base);
	my @baseStat = stat($base);
	my @assemblyStat = stat("$base.assStat");
	return 0 unless @baseStat && @assemblyStat && $assemblyStat[7] > 0;
	if ($useCheckM1) {
		my @checkM1Stat = stat("$base.cm");
		return 0 unless @checkM1Stat;
	}
	if ($useCheckM2) {
		my @checkM2Stat = stat("$base.cm2");
		return 0 unless @checkM2Stat && $checkM2Stat[7] > 0;
	}
	return 1;
}

sub emptyBinnerAssignmentCommand {
	my ($outDir, $name) = @_;
	die "emptyBinnerAssignmentCommand requires an output directory and name\n"
		unless defined($outDir) && length($outDir) && defined($name) && length($name);
	return "mkdir -p $outDir\n: > $outDir/$name\n";
}


sub readCMquals{
	my ($IQ) = @_;
	
	my %rQ;
	my $CM2mode=0; $CM2mode = 1 if ($IQ =~m/\.cm2/);
	open I,"<$IQ" or die "Can't open maxbin2 quality $IQ\n";
	while (<I>){
		chomp; my @spl  = split /\t/, $_, -1;
		next if ($spl[0] eq "Bin Id" || $spl[0] eq "Name");
		#die "can't find Bin \"$spl[0]\"\n" unless (exists($ret{$spl[0]}));
		my $Bin = shift @spl;
		if ($CM2mode){
			die "Malformed CheckM2 row for $Bin in $IQ\n" unless @spl >= 2 && $spl[0] =~ /^\d+(?:\.\d+)?$/ && $spl[1] =~ /^\d+(?:\.\d+)?$/;
			$rQ{$Bin}{compl} = $spl[0];
			$rQ{$Bin}{conta} = $spl[1];
			$rQ{$Bin}{hetero} = 0;
			$rQ{$Bin}{line} = join("\t",@spl);
		} else {
			die "Malformed CheckM row for $Bin in $IQ\n" unless @spl >= 13 && $spl[10] =~ /^\d+(?:\.\d+)?$/ && $spl[11] =~ /^\d+(?:\.\d+)?$/;
			$rQ{$Bin}{compl} = $spl[10];
			$rQ{$Bin}{conta} = $spl[11];
			$rQ{$Bin}{hetero} = $spl[12];
			$rQ{$Bin}{line} = join("\t",@spl);
		}
	}
	close I;
	return (\%rQ);
}

sub minQualFilter($ $ $ $ $){
	my ($hr1,$hr2,$Compl, $Conta, $LCAcompl) = @_;
	#%MB = %{$hr1}; %MBQ = %{$hr2};
	#foreach my $bin (keys %MB){
	my $prevS = scalar(keys %{$hr1});
	#print "$prevS , " . scalar(keys %{$hr2}) . "size\n";
	foreach my $bin (keys %{$hr1}){
		my $fails_lca = $LCAcompl > 0
			&& exists($hr2->{$bin}{LCAcompl})
			&& $hr2->{$bin}{LCAcompl} < $LCAcompl;
		if ($hr2->{$bin}{compl}< $Compl || $hr2->{$bin}{conta}> $Conta || $fails_lca || scalar(@{$hr1->{$bin}}) == 0){
			delete ($hr1->{$bin});
			delete ($hr2->{$bin});
		}
	}
	my $delEntries = $prevS - scalar(keys %{$hr1});
	print "Deleted $delEntries MAGs due to not meeting min qual criteria (Compl:$Compl, Conta:$Conta, LCA:$LCAcompl where available)\n" if ($delEntries);
	return ($hr1, $hr2);
}

sub MB2assigns($ $){
	my ($inF,$IQ) = @_;
	my %ret;
	#print "$inF\n";
	open I,"<$inF" or die "Can't open Binner output $inF\n";
	while (<I>){
		chomp; next if /^\s*$/;
		my @spl  = split /\t/, $_, -1;
		die "Malformed binner assignment in $inF at line $.\n" unless @spl >= 2 && length($spl[0]) && length($spl[1]);
		next if $spl[0] eq 'Sequence ID';
		next if ($spl[1] eq "0");
		push(@{$ret{$spl[1]}}, $spl[0]);
	}
	close I;
	
	my $rQHR = readCMquals($IQ);
	foreach my $bin (keys %ret) {
		die "No quality record for assigned bin '$bin' in $IQ\n" unless exists $rQHR->{$bin};
	}
	#print "$inF, $IQ ". scalar(keys %{$rQHR}) ."\n";

	return (\%ret,$rQHR);
}

sub MB2assignedBinIds {
	my ($inF, $IQ) = @_;
	my %assigned;
	open my $input, '<', $inF or die "Can't open Binner output $inF\n";
	while (my $line = <$input>) {
		chomp $line;
		next if $line =~ /^\s*$/;
		my @fields = split /\t/, $line, -1;
		die "Malformed binner assignment in $inF at line $.\n"
			unless @fields >= 2 && length($fields[0]) && length($fields[1]);
		next if $fields[0] eq 'Sequence ID';
		next if $fields[1] eq '0';
		$assigned{$fields[1]} = 1;
	}
	close $input or die "Cannot close Binner output $inF: $!\n";

	my $quality = readCMquals($IQ);
	for my $bin (keys %assigned) {
		die "No quality record for assigned bin '$bin' in $IQ\n"
			unless exists $quality->{$bin};
	}
	return (\%assigned, $quality);
}


#return how congruent LCA assignments of marker genes in MAG are
sub calcLCAcompl{
	my ($genesAR,$LCAHR)  = @_;
	my @genes = @{$genesAR}; my %LCA  = %{$LCAHR};
	#print "@genes\n";
	my %LCAcnts; my $hits =0 ;
	foreach my $gen (@genes){
		next unless (exists($LCA{$gen}));
		$hits++;
		my @locLCA = @{$LCA{$gen}};
		for (my $i=0;$i<@locLCA;$i++){
			$LCAcnts{$i}{$locLCA[$i]}++;
		}
	}
	#calc completeness per LCA level
	my @lvls = sort { $a <=> $b } keys %LCAcnts;
	my @complPerLvl = (); my @HitsPerLvl = ();
	foreach my $LV (@lvls){
		my %entr = %{$LCAcnts{$LV}};
		#my $ignCnt = $entr{"?"} if (exists($entr{"?"}));
		my $maxCnt = 0; my $sumCnt=0;
		foreach my $k (keys %entr){
			#print "$k $entr{$k} ; ";
			if ($k eq "?"){ next;}
			$sumCnt += $entr{$k};
			$maxCnt = $entr{$k} if ($maxCnt < $entr{$k});
		}
		if ($sumCnt > 0){
		$complPerLvl[$LV] = $maxCnt/$sumCnt;
		} else {$complPerLvl[$LV] =  0;}
		$HitsPerLvl[$LV] = $sumCnt;
		#print "  .. $maxCnt $sumCnt $complPerLvl[$LV] \n";
	}
	#scale error on higher tax levels for compound scoring
	my $maxL = @lvls;
	if ($maxL ==0) {return 0;}
	my $maxScore = 0;my $LCAcompl = 0;
	foreach my $LV (@lvls){
		$LCAcompl += ($maxL-$LV) * $HitsPerLvl[$LV]/$hits * $complPerLvl[$LV];
		$maxScore+=($maxL-$LV) * $HitsPerLvl[$LV]/$hits;
		#print "L$LCAcompl $maxScore :: ";
	}
	if ($maxScore == 0) {return 0;}
	#print "	$LCAcompl /= $maxScore;";
	$LCAcompl /= $maxScore;
	#print "Found $hits genes with LCA\n";
	#print $LCAcompl."\n";
	#die;
	return $LCAcompl;
}



#get list of MGS, based on checkM filtering
sub filterMGS_CM{
	my $CMfile=$_[0];
	my $complT=50; my $contT=5; my $retBetter=1;
	$complT=$_[1] if (@_ > 1); 
	$contT=$_[2] if (@_ > 2); 
	$retBetter=$_[3] if (@_ > 3); 
	my %ret; my $cnt =0;
	my $CMtag = "CM";
	my $contIdx = 12; my $complIdx = 11;
	if ($CMfile =~ m/\.cm2/){
		$CMtag = "CM2";
		$contIdx = 2; $complIdx = 1;
	}
	open I,"<$CMfile" or die "Binning.pm::filterMGS_CM: can't open $CMtag file $CMfile\n";
	while (<I>){
		$cnt++; next if ($cnt == 1);
		chomp;my @spl=split/\t/;
		#next if ($spl[0] eq "Bin Id");
		if (($spl[$complIdx] < $complT || $spl[$contIdx] > $contT) ){
			$ret{$spl[0]} = [$spl[$complIdx],$spl[$contIdx]] if (!$retBetter);
			#die "$spl[11],$spl[12]\n" if ($spl[0] eq "MB2bin12");
		} elsif ($retBetter) {
			$ret{$spl[0]} = [$spl[$complIdx],$spl[$contIdx]];
		}
	}
	close I;
	return(\%ret);
}

 sub deNovo16S{
	#returns 1)hash with rDNAs and 2) name of best seq
	my ($fasFile,$outfile) = @_;
	die "Fasta $fasFile doesn't exist\n" unless (-e $fasFile);
	my %ret; my $refSize = 1580;
	print "Detecting de novo 16S\n";
	
	my $newGff = "$outfile.ribo.gff";
	system "rm -f $newGff" if (-e $newGff);
	my $rnaBin = getProgPaths("rnammer");
	my $cmd = "$rnaBin -S bac -m ssu -gff $newGff < $fasFile";
	systemW $cmd."\n";
	#FP929038        RNAmmer-1.2     rRNA    507757  508732  1202.9  -       .       16s_rRNA
	my $genoR = readFasta($fasFile,1); my %geno = %{$genoR};
	#my @k = keys(%geno); die @k.join(" ",@k)."\n";
	open II,"<",$newGff or die "Can't open RNAmmer output $newGff: $!\n"; my $gffcnt=0;
	while(<II>){next if (/^#/ || length($_) < 1);my @spl = split(/\s+/);
		#print $_;
		next unless (@spl >= 7 && exists $geno{$spl[0]});
		my $newS = substr($geno{$spl[0]},$spl[3]-1,$spl[4]-$spl[3]+1);
		if ($spl[6] eq "-"){$newS = reverse_complement_IUPAC($newS);}
		#die $newS;
		$ret{$spl[0]."_rrn_".$gffcnt} = $newS;
		$gffcnt++;
	}
	close II;
	#check for new better 16S
	my $bestV=100000;
	my @head = keys(%ret);	my $best = 0; 	my $cnt =0;
	if (@head == 0){
		#die "Empty array\n$newGff\n\n$cmd\n";
		print "Empty 16S array\n$newGff\n";
		return (\%ret,"");
	} else {
		foreach my $hd (@head){
			my $tmp = abs($refSize - length($ret{$hd}));
			if ($tmp < $bestV){$bestV = $tmp; $best = $cnt;}
			$cnt++;
			#print $tmp."\n";
		}
	}
	#print "16S ".$bestV."\n";
	my %ret3 = ($head[$best],$ret{$head[$best]});
	return (\%ret,$head[$best]);
}

sub readBinSB($){
	my ($inF) = @_;
	my %can2gene;
	#print "reading SemiBin file: $inF\n";
	open I,"<$inF" or die "can't open canopy file $inF\n";
	while (<I>){
		chomp; my @spl = split /\t/;
		#die "$spl[0] $spl[1]\n";
		push(@{$can2gene{$spl[1]}},$spl[0]);
	}
	close I;
	return (\%can2gene);

}

sub readMGS($){
	my ($inF) = @_;
	my %can2gene;
	#print "reading MGS file: $inF\n";
	open I,"<$inF" or die "can't open canopy file $inF\n";
	while (<I>){
		chomp; my @spl = split /\t/;
		#die "$spl[0] $spl[1]\n";
		push(@{$can2gene{$spl[0]}},$spl[1]);
	}
	close I;
	return (\%can2gene);
}
sub readMGSrev($){
	my ($inF) = @_;
	my %g2c;
	my $dbl=0;
	#print "reading MGS reverse: $inF\n";
	open I,"<$inF" or die "can't open canopy file $inF\n";
	while (<I>){
		chomp; my @spl = split /\t/;
		#die "$spl[0] $spl[1]\n";
		#next unless (defined $spl[0] && defined $spl[1]);
		$dbl++ if (exists($g2c{$spl[1]}));
		$g2c{$spl[1]} = $spl[0];
	}
	close I;
	print "$dbl double gene assignments\n";
	return (\%g2c);
}

sub readMGSrevRed($){
	my ($inF) = @_;
	my %g2c;
	my $dbl=0;
	#print "reading MGS reverse: $inF\n";
	open I,"<$inF" or die "can't open canopy file $inF\n";
	while (<I>){
		chomp; my @spl = split /\t/;
		#die "$spl[0] $spl[1]\n";
		#next unless (defined $spl[0] && defined $spl[1]);
		$dbl++ if (exists($g2c{$spl[1]}));
		$g2c{$spl[1]}{ $spl[0]} = 1;
	}
	close I;
	print "$dbl double gene assignments\n";
	return (\%g2c);
}


sub getRepresentBins{
	my ($guide) = @_;
	print "Reading guidMGS file $guide\n";
	my ($I,$OK) = gzipopen($guide,"MGSvsGC",1);
	my @hd = (); my $lastMGS = "";
	my $repIdx = 0; my $ComplIdx=0; my $ContaIdx=0; my $LCAidx=0;my $N50idx=0;
	my %ret;
	while (my $line = <$I>){
		chomp $line; my @spl =  split /\t/,$line;
		if (!@hd){
			@hd = @spl;
			while ($hd[$repIdx] ne "Representative4MGS" && ($repIdx+2) < @hd){$repIdx++;}
			die "Couldn't find \"Representative4MGS\" string in @hd\n!\n" if ($repIdx >= @hd) ;
			while ($hd[$ComplIdx] ne "Completeness" && ($ComplIdx+2) < @hd){$ComplIdx++;}
			while ($hd[$ContaIdx] ne "Contamination" && ($ContaIdx+2) < @hd){$ContaIdx++;}
			while ($hd[$LCAidx] ne "LCAcompleteness" && ($LCAidx+2) < scalar(@hd)){$LCAidx++;}
			while ($hd[$N50idx] ne "N50" && ($N50idx+2) < @hd){$N50idx++;}
			#die "$repIdx = 0; my $ComplIdx=0; my $ContaIdx=0; my $LCAidx=0;my $N50idx=0;\n@hd \n";
			next;
		}
		
		my $curMGS = $spl[1];
		if ($curMGS ne $lastMGS){ #new cycle
			$lastMGS = $curMGS;
			$ret{$curMGS} = $spl[0] if ($spl[0] !~ m/^Cano__/); #set default to last one..
		}
		if ($spl[0] !~ m/^Cano__/ && ($spl[$repIdx] eq "*" || !defined($ret{$curMGS}) ) ){
			$ret{$curMGS} = $spl[0] ;
			#print "MGS $curMGS repr: $spl[0]\n";
		}
	}
	close $I;
	
	print "Found representative for ". scalar(keys %ret). " MGS\n";
	return \%ret;
}




sub getRepresentBinsPerFamily{ #needs some work
	my ($guide,$hrMap) = @_;
	print "Reading guidMGS file $guide\n";
	my ($I,$OK) = gzipopen($guide,"MGSvsGC",1);
	my @hd = (); my $lastMGS = "";
	my $repIdx = 0; my $ComplIdx=0; my $ContaIdx=0; my $LCAidx=0;my $N50idx=0;
	my $rejQuali = 0; my $famFnd=0; my $assmGrpFnd=0; my $noGrpFnd=0;
	my %ret;
	my %map = %{$hrMap};
	
	
	my %conta; my %compl; my %fam; my %score;
	
	
	#read in family info for each sample
	my %famSmpl;
	my @smpls = @{$map{opt}{smpl_order}};
	foreach my $smpl (@smpls){
		#$famSmpl{$smpl} = $map{$smpl}{AssGroup};
		if ($map{$smpl}{FamGroup} ne ""){
			$famSmpl{$smpl} = $map{$smpl}{FamGroup};$famFnd++;
		} elsif ($map{$smpl}{AssGroup} ne ""){#fallback assembly group
			$famSmpl{$smpl} = $map{$smpl}{AssGroup}; $assmGrpFnd++;
		} else {
			$famSmpl{$smpl} = $smpl; $noGrpFnd++;
		}
	}

	my $flushFamilyRepresentatives = sub {
		return unless length($lastMGS);
		foreach my $ff (sort keys %fam){
			my @eligible;
			foreach my $lBin (sort @{$fam{$ff}}){
				if ($compl{$lBin} > 60 && $conta{$lBin} < 20){
					push @eligible, $lBin;
				} else {
					$rejQuali++;
				}
			}
			unless (@eligible){
				warn "No eligible representative bin for family '$ff' in MGS '$lastMGS'; skipping\n";
				next;
			}
			my ($bestBin) = sort {
				$score{$b} <=> $score{$a} || $a cmp $b
			} @eligible;
			my $MGSfam = $ff . "." . $lastMGS;
			$ret{$MGSfam} = $bestBin;
		}
	};

	while (my $line = <$I>){
		chomp $line; my @spl =  split /\t/,$line;
		if (!@hd){
			@hd = @spl;
			while ($hd[$repIdx] ne "Representative4MGS" && ($repIdx+2) < @hd){$repIdx++;}
			die "Couldn't find \"Representative4MGS\" string in @hd\n!\n" if ($repIdx >= @hd) ;
			while ($hd[$ComplIdx] ne "Completeness" && ($ComplIdx+2) < @hd){$ComplIdx++;}
			while ($hd[$ContaIdx] ne "Contamination" && ($ContaIdx+2) < @hd){$ContaIdx++;}
			while ($hd[$LCAidx] ne "LCAcompleteness" && ($LCAidx+2) < scalar(@hd)){$LCAidx++;}
			while ($hd[$N50idx] ne "N50" && ($N50idx+2) < @hd){$N50idx++;}
			#die "$repIdx = 0; my $ComplIdx=0; my $ContaIdx=0; my $LCAidx=0;my $N50idx=0;\n@hd \n";
			next;
		}
		
		
		
		my $curMGS = $spl[1]; my $curBin = $spl[0];
		if ($curMGS ne $lastMGS){
			# Evaluate the preceding MGS before resetting its candidate state.
			$flushFamilyRepresentatives->();
			
			#reset params
			$lastMGS = $curMGS;
			%conta = (); %compl = (); %fam=(); %score = ();
		}
		
		if ($spl[0] !~ m/^Cano__/){
			#$ret{$curMGS} = $spl[0] ;
			$conta{$curBin} = $spl[$ContaIdx];$compl{$curBin} = $spl[$ComplIdx];
			$score{$curBin} = $spl[$ComplIdx] - (2. * $spl[$ContaIdx]);
			unless ($spl[0] =~ m/^(.+)__(.+)$/ && defined($famSmpl{$1})){
				warn "Cannot assign representative candidate '$curBin' in MGS '$curMGS' to a known sample family; skipping\n";
				next;
			}
			my $cFam = $famSmpl{$1};
			push(@{$fam{$cFam}}, $curBin);
		}	
	}
	# The loop transition flushes every group except the final MGS in the file.
	$flushFamilyRepresentatives->();
	close $I;
	
	print "Found representative for ". scalar(keys %ret). " MGS, $rejQuali rejected on Quali\n";
	print "Identified $famFnd families, $assmGrpFnd assembly groups, $noGrpFnd fallbacks\n";
	return \%ret;
}

#extract 1 reference genome per MGS, choosing the default rep and extracting its contigs from original assembly file
sub createBinCtgs{
	#$binDctg,$hrM,"$binDir/MAGvsGC.txt.gz
	my ($outD,$hrMap,$guideF,$perFam,$BinShrt) = @_;
	die "Output directory is required\n" unless defined($outD) && length($outD);
	die "Mapping data are required\n" unless ref($hrMap) eq 'HASH';
	die "Representative-bin guide is missing: $guideF\n" unless defined($guideF) && -s $guideF;
	make_path($outD);
	
	my $hr;
	if ($perFam){
		print "Getting per family ref genomes for MAGs\n";
		$hr = getRepresentBinsPerFamily($guideF,$hrMap);
	} else {
		print "Getting per AssmblGrp ref genomes for MAGs\n";
		$hr = getRepresentBins($guideF);
	}
	my %repBins = %{$hr};
	my %map = %{$hrMap};
	my %representatives_by_sample;

	for my $MGS (sort keys %repBins) {
		my $MAG = $repBins{$MGS};
		if ($MAG =~ m/^Cano__/) {
			print STDERR "Could not retrieve MAG for $MGS : $MAG , because is Canopy\n";
			next;
		}
		unless ($MAG =~ m/^(.+)__(.+)$/) {
			print STDERR "Could not match \"$MAG\" to sample and contig!\n";
			next;
		}
		my ($smpl, $bin) = ($1, $2);
		$smpl = $map{altNms}{$smpl} if defined($map{altNms}{$smpl});
		die "No assembly mapping for representative sample '$smpl'\n"
			unless exists($map{$smpl}) && defined($map{$smpl}{wrdir});
		push @{$representatives_by_sample{$smpl}{$bin}}, {
			mgs => $MGS,
			mag => $MAG,
		};
	}

	for my $smpl (sort keys %representatives_by_sample) {
		my $dirIn = $map{$smpl}{wrdir};
		my $assDir = getAssemblPath($dirIn);
		my $BinFile = "$assDir/Binning/$BinShrt/$smpl";
		my %wanted_bins = map { $_ => 1 } keys %{$representatives_by_sample{$smpl}};
		my (%contigs_by_bin, %wanted_contigs);

		open my $bin_input, '<', $BinFile or die "can't open bin file $BinFile\n";
		while (my $line = <$bin_input>) {
			chomp $line;
			next if $line =~ /^\s*$/;
			my @fields = split /\t/, $line, -1;
			die "Malformed bin assignment in $BinFile at line $.\n"
				unless @fields >= 2 && length($fields[0]) && length($fields[1]);
			next if $fields[0] eq 'Sequence ID';
			next unless $wanted_bins{$fields[1]};
			push @{$contigs_by_bin{$fields[1]}}, $fields[0];
			$wanted_contigs{$fields[0]} = 1;
		}
		close $bin_input or die "Cannot close bin file $BinFile: $!\n";

		for my $bin (keys %wanted_bins) {
			die "Representative bin '$bin' was not found for sample '$smpl'\n"
				unless exists($contigs_by_bin{$bin}) && @{$contigs_by_bin{$bin}};
		}
		my $sequences = readFasta(
			"$assDir/scaffolds.fasta.filt", 1, "\\s", \%wanted_contigs,
		);

		for my $bin (sort keys %{$representatives_by_sample{$smpl}}) {
			for my $representative (@{$representatives_by_sample{$smpl}{$bin}}) {
				my $outF = "$outD/$representative->{mgs}.ctgs.$representative->{mag}.fna";
				my $temporary = "$outF.tmp.$$";
				open my $output, '>', $temporary or die "Couldn't open $temporary\n";
				for my $ctg (@{$contigs_by_bin{$bin}}) {
					die "Contig '$ctg' from bin '$bin' is absent from the assembly for '$smpl'\n"
						unless exists($sequences->{$ctg});
					print {$output} ">$ctg\n$sequences->{$ctg}\n"
						or die "Cannot write $temporary: $!\n";
				}
				close $output or die "Cannot close $temporary: $!\n";
				unlink $outF or die "Cannot replace existing output $outF: $!\n" if -e $outF;
				rename $temporary, $outF
					or die "Cannot publish representative contigs $outF: $!\n";
			}
		}
	}

	print "----------------------\nDone\nWrote representative MAGs (contigs) to $outD\n----------------------\n";
}

sub createBin2{
	my ($binD,$cnopyF,$refFA) = @_;
	my $hr = readMGSrevRed($cnopyF);
	my ($I,$OK) = gzipopen($refFA,"reference gene cat",1);
	my $seq=""; my $hd="";
	my $geneCnt=0;
	my $fileEnd = ".fna"; $fileEnd = ".faa" if ($refFA =~ m/\.faa$/);
	my $max_open_outputs = 64;
	my (%open_outputs, %last_used, %temporary_outputs, %mgs_written);
	my $access_counter = 0;

	make_path($binD);
	my $output_handle = sub {
		my ($MGS) = @_;
		$access_counter++;
		if (exists $open_outputs{$MGS}) {
			$last_used{$MGS} = $access_counter;
			return $open_outputs{$MGS};
		}
		if (keys(%open_outputs) >= $max_open_outputs) {
			my ($oldest) = sort {
				$last_used{$a} <=> $last_used{$b} || $a cmp $b
			} keys %open_outputs;
			close $open_outputs{$oldest}
				or die "Cannot close temporary MGS output $temporary_outputs{$oldest}: $!\n";
			delete $open_outputs{$oldest};
			delete $last_used{$oldest};
		}
		my $temporary = $temporary_outputs{$MGS} ||= "$binD/$MGS$fileEnd.tmp.$$";
		my $mode = $mgs_written{$MGS} ? '>>' : '>';
		open my $output, $mode, $temporary
			or die "Cannot open temporary MGS output $temporary: $!\n";
		$open_outputs{$MGS} = $output;
		$last_used{$MGS} = $access_counter;
		return $output;
	};

	my $store_record = sub {
		return unless $hd =~ m/^>(\d+)/;
		my $gene_id = $1;
		return unless exists($hr->{$gene_id});
		$geneCnt++;
		for my $MGS (sort keys %{$hr->{$gene_id}}) {
			my $output = $output_handle->($MGS);
			print {$output} "$hd\n$seq\n"
				or die "Cannot write temporary MGS output $temporary_outputs{$MGS}: $!\n";
			$mgs_written{$MGS} = 1;
		}
	};
	while (my $line = <$I>){
		chomp $line;
		if ($line =~ m/^>/){
			$store_record->();
			$hd = $line; $seq = "";  
			next;
		}
		$seq .= $line;
	}
	$store_record->();
	close $I or die "Cannot close reference gene catalogue $refFA: $!\n";
	for my $MGS (keys %open_outputs) {
		close $open_outputs{$MGS}
			or die "Cannot close temporary MGS output $temporary_outputs{$MGS}: $!\n";
	}
	my $mgs_count = scalar keys %mgs_written;
	die "No genes from $cnopyF were found in $refFA\n" unless $mgs_count;
	print "Found $geneCnt genes in $mgs_count MGS (avg " . int($geneCnt/$mgs_count*100)/100  . " per MGS). Writing to $binD\n";
	for my $MGS (sort keys %mgs_written) {
		my $output = "$binD/$MGS$fileEnd";
		unlink $output or die "Cannot replace existing MGS output $output: $!\n" if -e $output;
		rename $temporary_outputs{$MGS}, $output
			or die "Cannot publish MGS output $output: $!\n";
	}
	print "----------------------\nDone\nWrote representative MAGs (genes) to $binD\n----------------------\n";
	
}

sub createBinFAA{
	my ($binD,$cnopyF,$refFA) = @_;
	my $suffix = "faa";
	$suffix = $_[3] if (@_ > 3);
	print "Reading reference MGS $cnopyF\n";
	my $clustHR = readMGS($cnopyF);#my %clust = %{$hr};
	print "Reading ref FAA $refFA\n";
	my $faaHR = readFasta($refFA,1);
	#my %FAA = %{$faaHR};
	
	system "mkdir -p $binD" unless (-d $binD);
	foreach my $cl (sort keys %{$clustHR}){
		my $oF = "$binD/$cl.$suffix";
		my @refG = @{$clustHR->{$cl}};
		my %alreadySeen;
		my $ostr = ""; my $gcnt=0;
		foreach my $rg (@refG){
			unless(exists($faaHR->{$rg})){
				print "Can't find gene $rg in gene cat\n" ;
				next;
			}
			my $rn = $rg; $rn =~ s/://;
			next if (exists($alreadySeen{$rn}));
			$ostr .= ">$rn\n$faaHR->{$rg}\n";
			$alreadySeen{$rn} = 1;
			$gcnt ++;
		}
		if ($ostr ne ""){
			open O,">$oF" or die $!;
			print O $ostr;
			close O;
		}
	}
}




sub MB2N50($){
	my ($hr) = @_;
	my %ret;
	my %M = %{$hr};
	my %sizes_to_shorthand = (1000     => '1K',
							  10000    => '10K',
							  100000   => '100K',
							  1000000  => '1M',
							  10000000 => '10M');
	foreach my $k(keys(%M)){
		my @mem = @{$M{$k}};
		my $totL=0; my $ctgs=0; my @lengs;
		#die @mem;
		foreach my $x (@mem){
			die "Cannot determine contig length from bin member '$x'\n"
				unless $x =~ m/_L=(\d+)=/;
			my $length = $1;
			push(@lengs,$length);
			$totL += $length; $ctgs++;
		}
		my $meanL = $totL/$ctgs;
		
		$ret{$k}{tL} = $totL; $ret{$k}{meanL} = $meanL;$ret{$k}{cN} = $ctgs;
		# find number of sequences above certain sizes
		foreach my $size (1000,10000,100000,1000000){
			$ret{$k}{$sizes_to_shorthand{$size}}=0;
			foreach my $l (@lengs){
				$ret{$k}{$sizes_to_shorthand{$size}} ++ if ($l>=$size);
			}
		}
		#and find N50
		@lengs = sort { $b <=> $a } @lengs;
		my $N20 = int ($totL *0.2); my $N50 = int ($totL *0.5);my $N80 = int ($totL *0.8);
		my $cumL=0;
		foreach my $l (@lengs){
			$cumL += $l;
			if (!exists($ret{$k}{N20}) && $cumL >= $N20){$ret{$k}{N20} = $l;}
			if (!exists($ret{$k}{N50}) && $cumL >= $N50){$ret{$k}{N50} = $l;}
			if (!exists($ret{$k}{N80}) && $cumL >= $N80){$ret{$k}{N80} = $l;}
		}

	}
	return \%ret;
}



sub runCheckM{#runs checkM on *.faa files (each file one Bin)
	my ($binD,$outFile,$tmpD,$ncore) = ($_[0],$_[1],$_[2],$_[3]);
	my $runNow = 1; $runNow = $_[4] if (@_ > 4);
	my $ext = "faa"; $ext = $_[5] if (@_ > 5);
	my $gtag = "--genes"; $gtag = "" if ($ext eq "fna");
	#system "rm -rf $tmpD/CM/";
	#system "mkdir -p $tmpD/tmp/" unless(-d "$tmpD/tmp/");
	my $cmC = "set -e\n";
	$cmC .= "rm -rf $tmpD/CM/;mkdir -p $tmpD/tmp/\n";
	#my $p2a = getProgPaths("py2activate");
	#my $pd = getProgPaths("pydeacti");
	my $checkMBin = getProgPaths("checkm");
	#$cmC .= "$p2a\n";
	$cmC .= "$checkMBin lineage_wf $gtag -x $ext -t $ncore --tab_table -f $outFile -q --pplacer_threads 3 --tmpdir $tmpD/tmp/ $binD $tmpD/CM/\n";
	$cmC .= "test -s $outFile\n";
	#$cmC .= "$pd\n";
	$cmC .= "rm -rf $tmpD/CM/ $tmpD/tmp/\n";
	
	#first check there's something in the bin file..
	my $outFile2 = $outFile; $outFile2 =~ s/\.cm$//;
	if (-e $outFile2){
		open I,"<$outFile2" or die "No input for runCheckM function:$outFile\n";
		my %binsFnd;
		while (my $l =<I>){
			chomp $l;
			next if $l =~ /^\s*$/;
			my @spl = split /\t/,$l, -1;
			next if @spl >= 2 && $spl[0] eq 'Sequence ID';
			die "Malformed bin assignment in $outFile2: $l\n"
				unless @spl >= 2 && length($spl[1]);
			next if $spl[1] eq '0';
			$binsFnd{$spl[1]} = 1;
		}
		close I;
		my @bins = keys %binsFnd;
		print "Found ".scalar(@bins)." metag Bins\n";
		if (!@bins){$cmC = "set -e\ntouch $outFile\n";}
	}
	print "$cmC\n";
	
	systemW "$cmC" if ($runNow > 0);
	return $cmC;
}

sub runCheckM2{#runs checkM2 on *.faa files (each file one Bin)
	my ($binD,$outFile,$tmpD,$ncore) = ($_[0],$_[1],$_[2],$_[3]);
	my $runNow = 1; $runNow = $_[4] if (@_ > 4);
	my $ext = "faa"; $ext = $_[5] if (@_ > 5);
	$tmpD .= "/CM2";
	my $gtag = "--genes"; $gtag = "" if ($ext eq "fna");
	#system "rm -rf $tmpD/CM/";
	#system "mkdir -p $tmpD/tmp/" unless(-d "$tmpD/tmp/");
	my $cmC = "set -e\n";
	$cmC .= "rm -rf $tmpD/;  mkdir -p $tmpD/\n";
	my $outD = $outFile; $outD =~ s/\/[^\/]+$/\//; $outD .= "/CHM2/";
	#my $pd = getProgPaths("pydeacti");
	#my $condaA = getProgPaths("CONDA");
	my $checkM2Bin = getProgPaths("checkm2");
	my $chm2DB = getProgPaths("checkm2DB");
	
	#$cmC .= "$condaA\nconda activate checkm2\n";
	if ($checkM2Bin =~ m/ activate /){
		$checkM2Bin =~ s/activate (\S+)/activate $1\nexport CHECKM2DB=$chm2DB\n/;
	}
	
	
	$cmC .= "$checkM2Bin predict --input $binD $gtag -x $ext --force -t $ncore  --output-directory $tmpD \n";#--tmpdir $tmpD/tmp/ $binD $tmpD/CM/\n";
	#$cmC .= "$pd\n";
	#debugging only	
#	$cmC .= "mkdir -p $outD;cp -r $tmpD $outD\n"; #replace with more targeted function later
	$cmC .= "cp $tmpD/quality_report.tsv $outFile\n"; #replace with more targeted function later
	$cmC .= "test -s $outFile\n";
	$cmC .= "rm -rf $tmpD\n";
	
	if ($runNow > 0){
		print "running checkM2 local..\n";
		systemW "$cmC" ;
	}
	return $cmC;
}



#actual binners start here



sub createBams{
	my ($dirsAR,$tmpDir,$outDir,$nm,$fna,$cores,$fakeEmpty,$minBamSiz,$fmt) = @_;
	my @dirSS = @{$dirsAR};
	my $numBams = 0;
	my @BAMS;
	my %seen_inputs;
	my $uncramCmd = "";
	if ($fmt ne "bam" && $fmt ne "sam"){die"createBams:: fmt has to be either bam or sam\nAborting..\n";}
	foreach my $DDI (@dirSS){
		$numBams++;
		my @inputs;
		if (-f $DDI && $DDI =~ /\.(?:bam|cram)$/i) {
			push @inputs, $DDI if (-s $DDI > $minBamSiz);
		} else {
			$DDI =~ s{/$}{};
			my $marker = "$DDI/mapping/done.sto";
			unless (-s $marker){
				warn "createBams: missing non-empty $marker; skipping $DDI for $nm\n";
				next;
			}
			open my $marker_fh, '<', $marker or die "Cannot read $marker: $!\n";
			my $mapped_name = <$marker_fh>;
			close $marker_fh;
			chomp $mapped_name;
			my $named_mapping = "$DDI/mapping/$mapped_name";
			my @candidates = ($named_mapping);
			if ($mapped_name !~ /\.sup-smd\./i) {
				(my $supplemental = $named_mapping) =~ s/-smd\./.sup-smd./i;
				push @candidates, $supplemental if $supplemental ne $named_mapping;
			}
			for my $candidate (@candidates) {
				if (!-e $candidate && $candidate =~ /\.bam$/) {
					(my $cram = $candidate) =~ s/\.bam$/.cram/;
					$candidate = $cram if (-e $cram);
				} elsif (!-e $candidate && $candidate =~ /\.cram$/) {
					(my $bam = $candidate) =~ s/\.cram$/.bam/;
					$candidate = $bam if (-e $bam);
				}
				push @inputs, $candidate
					if (-s $candidate && -s $candidate > $minBamSiz && !$seen_inputs{$candidate}++);
			}
		}
		for (my $input_no = 0; $input_no < @inputs; $input_no++) {
			my $input = $inputs[$input_no];
			if ($fmt eq 'bam' && $input =~ /\.bam$/i) {
				push @BAMS, $input;
				next;
			}
			my $suffix = $input_no == 0 ? '' : '.sup';
			my $output = "$tmpDir/$nm.$numBams$suffix.$fmt";
			$uncramCmd .= cram2bsam($input,$fna,$output,$fmt eq 'bam' ? 1 : 2,$cores)
				unless (-s $output);
			push @BAMS, $output;
		}
	}
	#die "@dirSS\n@BAMS\n";
	if (@BAMS == 0 && $fakeEmpty){
		warn "No non-empty mapping files found for $nm; publishing an empty bin assignment\n";
		return (emptyBinnerAssignmentCommand($outDir, $nm), []);
	}
	return ($uncramCmd,\@BAMS);
}



sub runSemiBin{
	my ($jgO,$outDir, $tmpDir, $nm, $fna, $cores, $dirsAR, $seqTec, $giveSBenv ) = @_;
	#human_gut/dog_gut/ocean/soil/cat_gut/human_oral/mouse_gut/pig_gut/built_environment/wastewater/global
	#my $semibinGTDB = getProgPaths("semibinGTDB");
	#get list of bams/crams..

#die;
	# SemiBin2 can crash on very small alignment files.  Treat mappings of
	# 15 MiB or less as unusable for this binner; the other binners do not share
	# this restriction and keep their zero-byte-only cutoff below.
	my $fakeEmpty=1;my $minBamSiz = 15*1024*1024;
	my ($uncramCmd,$BAMSar) = createBams($dirsAR,$tmpDir,$outDir,$nm,$fna,$cores,$fakeEmpty,$minBamSiz,"bam");
	my @BAMS = @{$BAMSar};
	return $uncramCmd unless (@BAMS);
	my $numBams = @BAMS;
	
	
	# --environment human_gut, dog_gut, ocean, soil, cat_gut, human_oral, mouse_gut, pig_gut, built_environment, wastewater, chicken_caecum, global
	# A curated environment is valid only for single-sample mode.  Multiple
	# usable BAMs make this a multi-sample/coassembly run, for which SemiBin2
	# must train from the supplied samples instead of loading a pretrained
	# environment.  For one BAM, default to human_gut unless explicitly set.
	my $selectedEnvironment = $giveSBenv ne "" ? $giveSBenv : "human_gut";
	die "Invalid SemiBin2 environment '$selectedEnvironment'\n"
		unless $selectedEnvironment =~ /^[A-Za-z0-9_-]+$/;
	my $senv = $numBams == 1 ? "--environment $selectedEnvironment" : "";
	my $SBbin = getProgPaths("SemiBin2");
	my $smode = "single_easy_bin ";
	my $dflags = " --random-seed 555 --tmpdir $tmpDir -p $cores";
	my $seqType = "--sequencing-type=short_read ";
	$seqType = "--sequencing-type=long_read " if ($seqTec eq "PB" || $seqTec eq "ONT" || $seqTec eq "hybrid");#PAcBIo/ONT
	my $cmd1 = "\necho \"preparing BAMs..\"\n$uncramCmd\n\n";
	$cmd1 .= "echo \"CRAM->BAM finished\"\n";
#	my $cmd = "";
	#my $output = "$outD/$nm.semibin";
	my $cmd = "\n\n###Running SemiBin...\n";
	if (@BAMS == 1 && $jgO ne ""){
		$cmd .= "$SBbin $smode --depth-metabat2 $jgO -i $fna $senv $seqType --output $outDir $dflags\n";
#	} elsif (@BAMS == 1){ #by default run with def env, if only 1 bam.. too slow otherwise
#		$cmd .= "$SBbin $smode -i $fna -b  ". join(" ",@BAMS). " $senv $seqType --output $outDir $dflags\n"; 

	} else {
		#--reference-db-data-dir $semibinGTDB --training-type self 
		$cmd .= "$SBbin $smode -i $fna -b  ". join(" ",@BAMS). " $senv $seqType --output $outDir $dflags\n";
	}
	#die $cmd;
	return $cmd1.$cmd;

}


sub runMetaDecoder{
	my ($jgO,$outDir, $tmpDir, $nm, $fna, $cores, $dirsAR ) = @_;
	#human_gut/dog_gut/ocean/soil/cat_gut/human_oral/mouse_gut/pig_gut/built_environment/wastewater/global
	my $MDbin = getProgPaths("MetaDecoder");
	my $baseN = "$tmpDir/$nm";
	#get list of bams/crams..
	my $fakeEmpty=1;my $minBamSiz = 0;
	my ($uncramCmd,$BAMSar) = createBams($dirsAR,$tmpDir,$outDir,$nm,$fna,$cores,$fakeEmpty,$minBamSiz,"sam");
	my @SAMS = @{$BAMSar};
	return $uncramCmd unless (@SAMS);
	

	my $cmd = "###preparing SAMs..\n$uncramCmd\n\n";
	$cmd .= "$MDbin coverage -s ". join(" ",@SAMS). " -o $baseN.coverage --threads $cores --mapq 10\n";
	$cmd .= "$MDbin seed  --threads $cores -f $fna -o $baseN.SEED\n";
	$cmd .= "$MDbin cluster -f $fna -c $baseN.coverage -s $baseN.SEED -o $baseN --random_number 511 --no_clusters\n"; 
	$cmd .= "mkdir -p $outDir;\n \ncp $baseN.cluster $outDir/$nm;\n";
	
	#still need to transfer files from here..
	#die $cmd;
	
	return $cmd;

}


sub runSCGBinner{
	my ($jgO,$outDir, $tmpDir, $nm, $fna, $cores, $dirsAR,
		$batchSize, $eligibleContigs) = @_;
	$batchSize //= 1024;
	$eligibleContigs //= "unknown";
	die "SCGBinner batch size must be a positive integer\n"
		unless $batchSize =~ /^\d+$/ && $batchSize > 0;
	my $fakeEmpty=1;my $minBamSiz = 0;
	my ($uncramCmd,$BAMSar) = createBams($dirsAR,$tmpDir,$outDir,$nm,$fna,$cores,$fakeEmpty,$minBamSiz,"bam");
	my @BAMS = @{$BAMSar};
	return $uncramCmd unless (@BAMS);
	#die "runSCGBinner::@BAMS\n";
	my $SCGbin = getProgPaths("SCGBinner");
	my $cmd = "###preparing BAMs..\n$uncramCmd\n\n";
	$cmd .= "###Running SCGBinner...\n";
	$cmd .= "mkdir -p $outDir\n";
	$cmd .= "set +e\n";
	$cmd .= "$SCGbin -a $fna -o $tmpDir -b \"" . join(" ",@BAMS)
		. "\" -t $cores -p $batchSize\n";
	$cmd .= "scgbinner_status=\$?\nset -e\n";
	$cmd .= "if [[ \$scgbinner_status -ne 0 ]]; then\n";
	$cmd .= "  echo \"ERROR: SCGBinner failed after preflight found $eligibleContigs "
		. "contigs >=1000 bp and selected batch size $batchSize. "
		. "Inspect the preceding SCGBinner marker and feature-generation messages.\" >&2\n";
	$cmd .= "  exit \$scgbinner_status\nfi\n";
	$cmd .= "if [[ ! -e $tmpDir/scgbinner_res/SCGBINNER_result.tsv ]]; then\n";
	$cmd .= "  echo \"ERROR: SCGBinner exited successfully but did not create "
		. "scgbinner_res/SCGBINNER_result.tsv\" >&2\n  exit 1\nfi\n";
	# move result TSV, then clean intermediate dirs
	$cmd .= "mv $tmpDir/scgbinner_res/SCGBINNER_result.tsv $outDir/$nm\n";
	$cmd .= "rm -rf $tmpDir\n";##$outDir/scgbinner_res $outDir/data_augmentation\n";
	return $cmd;
}

sub runGenomeFace{
	my ($jgO,$outD,$nm, $fna, $ncore, $tmpDir ) = @_;

	##genomeface -i $FASTA -a jgi_depths.tsv -g gfmarkers.txt -m 1000 -o GF
	##sed -E 's|_[0-9]+ |\t|' FMGids.txt | awk -F'=' '$2>1000' > gfmarkers.txt
	##awk -F'\t' '{print $2"\t"$1}' GF/bins.tsv > GF/outfile && rm GF/bins.tsv

	if (-e "$outD/GeFa.sto"){ return ""; }

	my $GFbin = getProgPaths("GenomeFace");
	my $outFile = "$outD/$nm";
	if (!-e $fna){die "Can't find required scaffold file in GenomeFace routine: $fna\n"; }
	my $cmd = "";
	my $assD = dirname($fna);
	my $FMGdir = "$assD/ContigStats/FMG";
	my $gfRawD = "$outD/gfOut";
	$cmd .= "mkdir -p $tmpDir $outD\n";
	$cmd .= "rm -rf $gfRawD\n"; # clear any partial previous run; genomeface creates this dir itself
	$cmd .= "sed -E 's|_[0-9]+ |\\t|' $FMGdir/FMGids.txt | awk -F'=' '\$2>1000' > $tmpDir/gfmarkers.txt\n";
	$cmd .= "[[ ! -s $tmpDir/gfmarkers.txt ]] && echo 'No markers, exiting early.' && exit 1\n";
	$cmd .= "$GFbin -i $fna -a $jgO -g $tmpDir/gfmarkers.txt -m 1000 -o $gfRawD\n";
	$cmd .= "awk -F'\\t' '{print \$2 \"\\t\" \$1}' $gfRawD/bins.tsv > $outFile && rm -rf $gfRawD\n";
	$cmd .= "echo \"$nm\" > $outD/GeFa.sto\n";
	return $cmd;

}

sub runMetaBat{
	my ($jgO,$outD,$nm, $fna ) = @_;
	#if (-e "$outD/MeBa.sto"){		return "";	}
	my $numC = 0;
	$numC = $_[4] if (@_ > 4);
	#print "Running MetaBat..\n";
	my $metab2Bin = getProgPaths("metabat2");
	#my $mbBin = "/g/bork3/home/hildebra/bin/metabat/./metabat";
	my $outCtgNms = "$outD/$nm.ctgs.txt";
	my $outFna = "$outD/$nm.fasta.fna";
	my $outD2 = "$outD/$nm";
	my $outMat = "$outD/$nm.mat.txt";
	#my $jgO = "$tmp/depth.jgi";
	if (!-e $fna){die "Can't find requried scaffold file in metabat routine: $fna\n"; }
	#start metabat
	my $cmd = "";#$before."\n";
	#$cmd .= "$mbBin -i $fna -a $jgO.depth.txt -o $outFna -p $jgO.pairs.sparse -l $outCtgNms --minCVSum  10 -t $numCore --saveCls $outMat -v\n";
	$cmd .= "mkdir -p $outD\n";
	$cmd .= "$metab2Bin -i $fna -a $jgO -o $outD2 -l --saveCls -t $numC --noBinOut -m 2500 --minCVSum 1\n"; #$outMat
	#$cmd .= "echo \"$nm\" > $outD/MeBa.sto\n";
	return $cmd;
	#die $cmd."\n";
	#my $jobName = "MBat";
	#my ($jobNameX, $tmpCmd) = qsubSystem($logDir."metaBat.sh",$cmd,$numCore,"60G",$jobName,"","",1,[],$QSBoptHR);

}
