package Mods::IO_Tamoc_progs;
use warnings;
use Cwd 'abs_path';
use strict;

use vars qw($CONFIG_FILE @CONFIG_TEXT %CONFIG_HASH);
$CONFIG_FILE="";
@CONFIG_TEXT = ();
%CONFIG_HASH = ();
sub setConfigFile;

#TAMOC programs related to IO to other programs, program paths .. not real subroutines that do anything

use Exporter qw(import);
our @EXPORT_OK = qw(getProgPaths truePath
					inputFmtSpades inputFmtMegahit jgi_depth_cmd createGapFillopt setConfigFile 
					buildMapperIdx mapperDBbuilt decideMapper checkMapsDoneSH greaterComputeSpace convert2Gb);


sub decideMapper($ $){
	my ($MapperProg,$readTec) = @_;
	if ($MapperProg < 0 ) {
		if ($readTec eq "PB" || $readTec eq "ONT"){  #default to minimap2 for 3rd gen seq
			$MapperProg = 3;
		} else {#otherwise use bowtie2
			if ($MapperProg == -1){
				$MapperProg = 1 ;
			} elsif ($MapperProg == -2){
				$MapperProg = 5 ; #.. or strobealign
			}
		}
	}
	if ($MapperProg<1 || $MapperProg>5){
		die "IO_Tamoc_progs.pm::decideMapper:: unknown MapperProg provided: $MapperProg\n! Aborting..\n";
	}
	return $MapperProg;
}

sub checkMapsDoneSH{
	my ($inAR) = @_;
	my @dirSS = @{$inAR};
	my $ctrlStr = "";
	foreach my $DDI (@dirSS){
		if ( $DDI =~ m/\/$/  ){
			$ctrlStr .= "if [ ! -e $DDI/mapping/done.sto ] || ! find $DDI/mapping -maxdepth 1 -type f \\( -name '*-smd.bam' -o -name '*-smd.cram' \\) -size +0c -print -quit | grep -q .; then echo \"Can't find a completed mapping in $DDI/mapping !! Aborting .. \"; exit 1; fi \n";
		} else {
			$ctrlStr .= "if [ ! -s $DDI ]; then echo \"Can't find non-empty mapping file $DDI !! Aborting ..\"; exit 1; fi \n";
		}
	}
	return $ctrlStr;
}

#computes which string in the style of "120G" or "120" is greatest and returns this (in "G" format)
sub greaterComputeSpace{
	my @inA  = @_;
	my $ret = 0; #unit is G (gigabyte)
	foreach my $in (@inA){
		if ($in =~ m/\D/){
			#print "NON\n";
			if ($in =~ /^\d+G$/){
				$in =~ s/G//;$in = $in + 0;
			}elsif ($in =~ /^\d+M$/){
				$in =~ s/M//; $in = 0+$in; $in /= 1024;
			}elsif ($in =~ /^\d+T$/){
				$in =~ s/T//; $in = 0+$in; $in *= 1024;
			} else {
				die "Unrecognized compute space: $in\n";
			}
		} else{
			$in = 0+$in; #convert to num
		}
		$ret = $in if ($in > $ret);
	}
	#$ret .= "G";
	#die "$ret\n";
	return $ret;
}

#converts string of style "120G" or "120" to GB. If no "G/M/T" given, assumes already in Gb ("G")
sub convert2Gb($){
	my ($tmpSpace) = @_;
	if ($tmpSpace !~ m/\D/ && ($tmpSpace eq "0" || ($tmpSpace+0) == 0) ){$tmpSpace = 0 ;
	}elsif ($tmpSpace =~ s/G$//){$tmpSpace = int($tmpSpace+0.5);
	}elsif ($tmpSpace =~ s/T$//){$tmpSpace *= 1024 ;
	}elsif ($tmpSpace =~ s/M$//){$tmpSpace /= 1024 ;
	} else {$tmpSpace = int($tmpSpace+0.5);}
	return $tmpSpace;
}

sub setConfigFile{
	my @var = @_;
	my $customCfg = 0;
	if (@var == 1 && $var[0] eq "internal"){
		my $modDir = $INC{"Mods/IO_Tamoc_progs.pm"};
		$modDir =~ s/IO_Tamoc_progs.pm//;
		$CONFIG_FILE = "$modDir/../Mods/config_internal.txt";
	} elsif (@var == 1 && $var[0] eq "DBconfig"){
		my $modDir = $INC{"Mods/IO_Tamoc_progs.pm"};
		$modDir =~ s/IO_Tamoc_progs.pm//;
		$CONFIG_FILE = "$modDir/../Mods/config_DBs.txt";
	} elsif (@var == 1 && $var[0] ne ""){
		$CONFIG_FILE = $var[0];
		$customCfg = 1;
	} else {#default value
		my $modDir = $INC{"Mods/IO_Tamoc_progs.pm"};
		$modDir =~ s/IO_Tamoc_progs.pm//;
		$CONFIG_FILE = "$modDir/MATAFILERcfg.txt";
	}
	die "Can't find MATAFILER config file: $CONFIG_FILE\nConsider changing path to config file via \"-config\" argument.\n Aborting..\n" unless (-e $CONFIG_FILE);
	print "Using config file : $CONFIG_FILE\n" if ($customCfg);
}

sub truePath{
	my ($TMCpath) = $_[0];
	my $enforce=0; $enforce = $_[1] if (@_ > 1);
	
	if ($enforce){
		if ($TMCpath =~ m/\$([^\$^\/^\\]+)/){
			my $envName = $1;
			die "Environment variable \$$envName used in path '$TMCpath' is not set\n"
				unless (exists($ENV{$envName}) && defined($ENV{$envName}) && $ENV{$envName} ne '');
			my $envVar = $ENV{$envName};
			$TMCpath =~ s/\$([^\$^\/^\\]+)/$envVar/;
		}
		#die "$TMCpath\n";
	}
	
	if ($TMCpath =~ m/^\$/){
		$TMCpath =~ s/^\$//; 
		my ($envName, $suffix) = $TMCpath =~ m{^([^/\\]+)(.*)$};
		die "Environment variable \$$envName used in path '\$$TMCpath' is not set\n"
			unless (exists($ENV{$envName}) && defined($ENV{$envName}) && $ENV{$envName} ne '');
		$TMCpath = $ENV{$envName} . $suffix;
	}
	return $TMCpath;

}

sub loadConfigs{
	#loads once in every program run the entire config file(s) into hash %CONFIG_HASH
	if (scalar @CONFIG_TEXT == 0){
		setConfigFile() if ($CONFIG_FILE eq "");
		print "READING config files \"$CONFIG_FILE\" .. ";
		if ($CONFIG_FILE eq "" ){die "IO_Tamoc_progs.pm::loadConfigs: CONFIG_FILE not set!\n";}
		open I,"<$CONFIG_FILE" or die "Can't open $CONFIG_FILE\n";
		chomp(@CONFIG_TEXT = <I>);
		close I;
		setConfigFile("internal") ;
		if ($CONFIG_FILE eq "" ){die "IO_Tamoc_progs.pm::loadConfigs: CONFIG_FILE internal not set!\n";}
		open I,"<$CONFIG_FILE" or die "Can't open internal $CONFIG_FILE\n";
		my @INTtmp;
		chomp(@INTtmp = <I>);
		close I;
		push(@CONFIG_TEXT,@INTtmp);
		#DB config read..
		setConfigFile("DBconfig") ;
		if ($CONFIG_FILE eq "" ){die "IO_Tamoc_progs.pm::loadConfigs: CONFIG_FILE DB not set!\n";}
		open I,"<$CONFIG_FILE" or die "Can't open DBconfig $CONFIG_FILE\n";
		@INTtmp=();
		chomp(@INTtmp = <I>);
		close I;
		push(@CONFIG_TEXT,@INTtmp);
	}
	#my $condaA = getProgPaths("CONDA");
	#die "@CONFIG_TEXT\n";
	print "converting config files.. ";
	my $TMCpath = "";my $Tset=0; my $BINpath = "";my $Bset=0; my $Rpath=""; my $RpathSet=0;
	my $DBpath = "";my $DBset=0; my $SINGcmd = ""; my $SINGset=0; 
	my $CONDset = 0; my $CONDcmd="";my $CONDset2 = 0; my $CONDA="";my $CONDset3 = 0; my $CONDAbaseEnv="MFF";
	my $PY3cmd = ""; my $PY3set=0; my $Rscriptcmd = ""; my $Rscriptset=0; my $MGSTKDir = ""; my $MGSTKDirset=0;
	foreach my $l (@CONFIG_TEXT){
		next if ($l =~ m/^#/ || length($l) == 0);
		if (!$Tset && $l =~ m/^MFLRDir\t([^#]+)/){
			$Tset=1;$TMCpath = truePath($1);
#			next;
		} elsif (!$Bset && $l =~ m/^BINDir\t([^#]+)/){
			$BINpath = truePath($1); $Bset=1;
		} elsif (!$RpathSet && $l =~ m/^Rpath\t([^#]+)/){
			my $prePath = $1;
			if (!$Tset){die"Problem in configs: MFLRDir needs to be set before Rpath\n";}
			$prePath =~s/\[MFLRDir\]/$TMCpath/;
			$Rpath= truePath($prePath);
			$RpathSet=1;
			#die "\n\n$prePath\n$Rpath\n";
		} elsif (!$DBset && $l =~ m/^DBDir\t([^#]+)/){
			$DBpath = truePath($1); $DBset=1;
		} elsif (!$MGSTKDirset && $l =~ m/^MGSTKDir\t([^#]+)/){
			$MGSTKDir = truePath($1); $MGSTKDirset=1;
		} elsif (!$SINGset && $l =~ m/^SINGcmd\t([^#]+)/){
			$SINGcmd = truePath($1); $SINGset=1;
			#die "SINGcmd no longer supported\nPlease remove from Config\n";
		} elsif (!$CONDset && $l =~ m/^CONDcmd\t([^#]+)/){
			$CONDcmd = truePath($1); $CONDset=1;
		} elsif (!$CONDset2 && $l =~ m/^CONDA\t([^#]+)/){
			$CONDA = $1;
			if ($CONDA =~ m/\[CONDcmd\]/){
				die "Ensure \"CONDcmd\" is set in config before \"CONDA\"\n" if (!$CONDset);
				$CONDA =~ s/\[CONDcmd\]/$CONDcmd/;
			}
			if ($CONDA !~ m/shell hook/){
				$CONDA = truePath($CONDA,1); 
				my $Ctmp = $CONDA; $Ctmp =~ s/^[\.\s]+//g;
				if (!-s $Ctmp){die "\n\nWARNING:\nCould not find conda config at $CONDA !\n please ensure \"micromamba.sh\" or \"mambda.sh\" exist at this location\n\n";}
				if ($CONDA !~ m/^\./ || $CONDA !~ m/mamba.sh/){
					die "Your \"CONDA\" seems to be wrongly setup. Ensure this is configured in \"[MG-TK-dir]/config.txt\" and has a form similar to:\"\nCONDA\t. \$MAMBA_ROOT_PREFIX/etc/profile.d/mamba.sh\n\"\n";
				}
			}
			$CONDset2=1;
		} elsif (!$CONDset3 && $l =~ m/^CONDAbaseEnv\t([^#]+)/){
			$CONDAbaseEnv = truePath($1); $CONDset3=1;
			
		} elsif (!$PY3set && $l =~ m/^PY3cmd\t([^#]+)/){
			$PY3cmd = $1; $PY3set=1;
			if ($SINGset){$PY3cmd =~ s/\[SINGcmd\]/$SINGcmd/;}
		} elsif (!$Rscriptset && $l =~ m/^Rscript\t([^#]+)/){
			$Rscriptcmd = $1; $Rscriptset=1;
			#print "Rscript set!! : $Rscriptcmd  \n\n\n";
		
		
		} else { #bit strange way of doing this.. but compliant with old style
			my $XVar = "";
			#print "$l\n";
			my @spl = split (/\t/,$l);
			$XVar = $spl[0];
			if (@spl == 1) {$CONFIG_HASH{$XVar} = "";next;}
			
			if ($l !~ m/^$XVar\t([^#^\t]+)/){next;}
			my $reV = $1;
			
			#die "$reV  $XVar  $l\n";
			die "$reV\n" if (!defined($reV));
			$reV =~ s/\[MFLRDir\]/$TMCpath/ if ($Tset);
			if ($MGSTKDirset){
				$reV =~ s/\[MGSTKDir\]/$MGSTKDir/ ;
			}
			$reV =~ s/\[BINDir\]/$BINpath/ if ($Bset);
			$reV =~ s/\[DBDir\]/$DBpath/ if ($DBset);
			$reV =~ s/\[SINGcmd\]/$SINGcmd/ if ($SINGset);
			$reV =~ s/\[PY3\]/$PY3cmd/ if ($PY3set);
			$reV =~ s/\[Rscript\]/$Rscriptcmd/ if ($Rscriptset);
			$reV =~ s/\[Rpath\]/$Rpath/ if ($RpathSet);
			if ($l =~ m/env:([^#^\t]+)/){
				my $tarEnv = $1;
				#$reV = "$CONDA;$CONDcmd activate $1\n$reV";
				$reV = "$CONDA;if [[ \$CONDA_DEFAULT_ENV != $tarEnv ]]; then $CONDcmd activate $tarEnv; fi\n$reV";
			}
			
			#return $reV;
			$CONFIG_HASH{$XVar} = $reV;
		}
	}
	#some check about basic params being set..
	if (!$Rscriptset){die "Could not find \"Rscript\" correctly configured in config file, please check your local config! Aborting..\n";}
	if (!$CONDset){die "Could not find \"CONDcmd\" correctly configured in config file, please check your local config! Aborting..\n";}
	if (!$DBset){die "Could not find \"DBDir\" correctly configured in config file, please check your local config! Aborting..\n";}
	if (!$Tset){die "Could not find \"MFLRDir\" correctly configured in config file, please check your local config! Aborting..\n";}
	if (!$CONDset2){die "Could not find \"CONDA\" correctly configured in config file, please check your local config! Aborting..\n";}
	$CONFIG_HASH{"activateBase"} = "$CONDA;$CONDcmd activate $CONDAbaseEnv\n";
	$CONFIG_HASH{"CONDAbaseEnv"} = $CONDAbaseEnv;
	$CONFIG_HASH{"MFLRDir"} = $TMCpath;
	$CONFIG_HASH{"BINDir"} = $BINpath;
	$CONFIG_HASH{"DBDir"} = $DBpath;
	$CONFIG_HASH{"Rscript"} = $Rscriptcmd;
	$CONFIG_HASH{"Rpath"} = $Rpath;
	$CONFIG_HASH{"MGSTKDir"} = $MGSTKDir;
	print "  Done. ";
}


sub getProgPaths{
	my @var = @_;
	my $srchVar = $var[0] ;
	my $required=1;
	if (@var > 1){$required = $var[1];}
	#die "$required\n";

	my @multVars = ();
	if (ref $srchVar eq 'ARRAY') {
		#print "ARRAY\n";
		@multVars = @{$srchVar};
	}
	if (scalar(keys %CONFIG_HASH) == 0){
		#read in config hash _once_
		loadConfigs();
	}
	if (scalar(keys %CONFIG_HASH) == 0){
		die "Something went wrong loading MATAFILER configs.. aborting\n";
	}
	
	
	if (@multVars > 0){
		my @retA;
		for (my$i=0;$i<scalar(@multVars);$i++){if (exists($CONFIG_HASH{$multVars[$i]})) { $retA[$i] = $CONFIG_HASH{$multVars[$i]};}} 
		return \@retA;
	}
	if (exists($CONFIG_HASH{$srchVar})){
		return $CONFIG_HASH{$srchVar};
	} else {
		die "Can't find configuration for $srchVar in MATAFILER config ($CONFIG_FILE)\n" if ($required != 0);
	}
	
	return "";
}

sub activateBase{
	die "Mods/IO_Tamoc_progs.pm::activateBase:: not used";
	my $ret = "";
	#my $condaA = getProgPaths("CONDA");
	#$ret = "$condaA\n$CONDcmd activate $1\n$reV";
	return $ret;
}


sub mapperDBbuilt( $ $){
	my ($DBbtRef, $MapperProg2) = @_;
	my $bwt2IdxFileSuffix = ".bw2";my $mini2IdxFileSuffix = ".mmi";
	my $kmaIdxFileSuffix = ".kma";
	if ($MapperProg2 == 5){return 1;} #strobealign doesn't need index..
	#print "($MapperProg2 == 1 || $MapperProg2 == -1) && !-s $DBbtRef$bwt2IdxFileSuffix.rev.2.bt2\n";
	my @bt2_small = map { "$DBbtRef$bwt2IdxFileSuffix.$_.bt2" } qw(1 2 3 4 rev.1 rev.2);
	my @bt2_large = map { "$DBbtRef$bwt2IdxFileSuffix.$_.bt2l" } qw(1 2 3 4 rev.1 rev.2);
	my $bowtie_complete = !(grep { !-s $_ } @bt2_small) || !(grep { !-s $_ } @bt2_large);
	if ( 
		($MapperProg2 ==0 && !-e "$DBbtRef$bwt2IdxFileSuffix.0.sa") 
		|| ( ($MapperProg2 == 1 || $MapperProg2 == -1) && !$bowtie_complete ) #bowtie2
		||( $MapperProg2 == 2 && !-s "$DBbtRef.pac" ) #bwa
		||( ($MapperProg2 == 3 || $MapperProg2 == -1 ) && !-s "$DBbtRef$mini2IdxFileSuffix" ) #minimap2
		||( ($MapperProg2 == 4 ) && !-s "$DBbtRef$kmaIdxFileSuffix.seq.b" )#kma
	) {
		return 0;
	}
	return 1;
}

sub buildMapperIdx($ $ $ $){
	my ($REF,$ncore,$lrgDB,$MapperProg) = @_;
	#1=bowtie2, 2=bwa, 3=minimap2
	$MapperProg = decideMapper($MapperProg,"");
	if ($MapperProg == 5){return ("",$REF,$REF);} #strobealign doesn't need index..
	my $bwt2IdxFileSuffix = ".bw2";my $mini2IdxFileSuffix = ".mmi";
	my $kmaIdxFileSuffix = ".kma";
	my $bwtIdx = $REF.$bwt2IdxFileSuffix;
	my $chkFi = $bwtIdx;
	my @required_index_files;
	if ($MapperProg==1){
		my $extension = $lrgDB ? 'bt2l' : 'bt2';
		@required_index_files = map { "$bwtIdx.$_.$extension" } qw(1 2 3 4 rev.1 rev.2);
		$chkFi = $required_index_files[-1];
	}elsif ($MapperProg==2){$chkFi = $REF.".pac";
	} elsif ($MapperProg == 3){$chkFi = $REF.$mini2IdxFileSuffix;
	} elsif ($MapperProg == 4){$chkFi = $REF.$kmaIdxFileSuffix.".seq.b";
	}
	my $dbCmd ="";
	my $missing_test = @required_index_files
		? join(' || ', map { "[ ! -s $_ ]" } @required_index_files)
		: "[ ! -s $chkFi ]";
	$dbCmd .= "if $missing_test; then \n";
	$dbCmd .= "echo \"Building index for mapper $MapperProg\"\n";
	if ($MapperProg==1){
		my $bwt2Bin = getProgPaths("bwt2");
		$dbCmd .= $bwt2Bin."-build ";
		$dbCmd .= " --large-index "if ($lrgDB);
		$dbCmd .= "-q $REF --threads $ncore $bwtIdx\n";
		$dbCmd = "" if (mapperDBbuilt($REF,1));
	} elsif($MapperProg==2) { 
		my $bwaBin  = getProgPaths("bwa");
		$dbCmd .= $bwaBin." index $REF\n";
		if (-s $REF.".pac"){$dbCmd = "";} 
		#die $dbCmd."\n";
	} elsif ($MapperProg == 3){
		$bwtIdx = $REF.$mini2IdxFileSuffix;
		my $mini2  = getProgPaths("minimap2");
		$dbCmd .= "$mini2 -t $ncore -H -d $bwtIdx $REF\n";
		$dbCmd = "" if (-s $bwtIdx);
	} elsif ($MapperProg==4){ 			
		$bwtIdx = $REF.$kmaIdxFileSuffix;
		my $kmaBin = getProgPaths("kma");
		$dbCmd .= "$kmaBin index -i $REF -o $bwtIdx 2>/dev/null \n"; #-t $ncore 
	}

	$dbCmd .= "fi\n" unless ($dbCmd eq "");
	if ($MapperProg == 1 && $dbCmd eq "") {
		$chkFi = -s "$bwtIdx.rev.2.bt2l" ? "$bwtIdx.rev.2.bt2l" : "$bwtIdx.rev.2.bt2";
	}
	#die "$dbCmd\n";
	#my $jobN = "_ASDB$JNUM"; my $tmpCmd;
	#($jobN, $tmpCmd) = qsubSystem($logDir."BAM2CRAMxtra.sh",$dbCmd,1,"10G",$jobN,"","",1,[],\%QSBopt);
	return ($dbCmd,$bwtIdx,$chkFi);
}


sub inputFmtSpades($ $ $ $ $){
	my ($p1ar,$p2ar,$singlAr,$logDir,$cReadTecAr) = @_;
	my @p1 = @{$p1ar}; my @p2 = @{$p2ar};
	my @singl = @{$singlAr};
	my @readTec = @{$cReadTecAr};
	my $doYAML = 0;
	if (@p1 > 9){$doYAML=1;}
	if (@p1 != @p2){print "Unequal paired read array lengths arrays for Spades\n"; exit(2);}
	my $sprds = "";
	if ($doYAML==0){
		for (my $i =0; $i<@p1;$i++){
			next if ($p1[$i] eq "");
			my $peTerm = "--pe";$peTerm = "--gemcode" if ($readTec[$i] =~ m/SLR/);
			$sprds .= " ${peTerm}".($i+1) ."-1 $p1[$i] ${peTerm}".($i+1) ."-2 $p2[$i]";
		}
		for (my $i=0;$i<@singl;$i++){
			next if ($singl[$i] eq "");
			my $peTerm = "--pe";$peTerm = "--gemcode" if ($readTec[$i] =~ m/SLR/);
			$sprds .= " ${peTerm}".($i+1) ."-s $singl[$i]";
		}
	} else {
		open O,">$logDir/spadesInput.yaml" or die "Can't write $logDir/spadesInput.yaml\n";
		print O "  [\n      {\n        orientation: \"fr\",\n        type: \"paired-end\",\n        left reads: [\n";
		my @valid_pairs = grep { $p1[$_] ne "" } 0..$#p1;
		print O join(",\n", map { "          \"$p1[$_]\"" } @valid_pairs), "\n";
		print O "     ],\n        right reads: [\n";
		print O join(",\n", map { "          \"$p2[$_]\"" } @valid_pairs), "\n";
		print O "       ]\n      }";
		my @valid_singletons = grep { $_ ne "" } @singl;
		if (@valid_singletons > 0){
			print O ",\n      {\n        type: \"single\",\n        single reads: [\n";
			print O join(",\n", map { "          \"$_\"" } @valid_singletons), "\n";
			print O "       ]\n      }";
		}
		#end of yaml file
		print O "\n     ]\n";
		close O;
		$sprds = " --dataset $logDir/spadesInput.yaml";
	}
	#die "$sprds\n\ninputFmtSpades\n";
	return $sprds;
 }

 sub inputFmtMegahit($ $ $ $){
	my ($p1ar,$p2ar,$singlAr,$logDir) = @_;
	my @p1 = @{$p1ar}; my @p2 = @{$p2ar};
	my @singl = @{$singlAr};
	@p1 = grep !/^$/, @p1;
	@p2 = grep !/^$/, @p2;
	@singl = grep !/^$/, @singl;

	die "Unequal paired read array lengths for MEGAHIT\np1:@p1\np2:@p2\n"
		if (@p1 != @p2);
	die "No read inputs supplied to MEGAHIT\n" if (!@p1 && !@singl);
	my $sprds = "";
	
	if (@p1 > 0){ 
		$sprds .= "-1 ".join(",",@p1) . " -2 ".join(",",@p2)." ";
	}
	if (@singl > 0){
		$sprds .= "-r ".join(",",@singl);
	}
	return $sprds;
 }


sub jgi_depth_cmd{
	my $dirsAR = $_[0];
	my $out = $_[1];
	my $perID = $_[2];#) = @_;
	my $numCores = 1;
	$numCores = $_[3] if (@_ > 3);
	my $refFA = "";
	$refFA = $_[4] if (@_ > 4);
	#die $dirs."\n";
	my $smtBin = getProgPaths("samtools");
	my $jgiScr = getProgPaths("jgiDepth");

	my @dirSS = @{$dirsAR};#split(',',$dirs);
	#go through each dir and find sample name
	my $comBAM = "";
	my $isCram=0;
	foreach my $DDI (@dirSS){
		if (-f $DDI && $DDI =~ /\.(?:bam|cram)$/i) {
			$isCram=1 if ($DDI =~ /\.cram$/i);
			$comBAM .= "$DDI ";
		} else {
			$DDI =~ s{/$}{};
			my $marker = "$DDI/mapping/done.sto";
			die "jgi_depth_cmd:::Missing mapping marker $marker\n" unless (-s $marker);
			open my $marker_fh, '<', $marker or die "Cannot read $marker: $!\n";
			my $SmplNm = <$marker_fh>;
			close $marker_fh;
			chomp $SmplNm;
			my $tbam = "$DDI/mapping/$SmplNm";
			if (!-s $tbam && $tbam =~ /\.bam$/){(my $cram = $tbam) =~ s/\.bam$/.cram/; $tbam = $cram if (-s $cram);}
			die "jgi_depth_cmd:::Can't find a non-empty BAM or CRAM at $DDI\n" unless (-s $tbam);
			$isCram=1 if ($tbam =~ /\.cram$/i);
			$comBAM .= "$tbam ";
		}
	}
	#my $comBAM = join("/mapping/Align_ment-smd.bam ",@dirSS);
	# Split conda activation prefix from the binary name so the activation can be
	# emitted before any pipe, keeping just the bare binary after the pipe.
	my ($jgiActivate, $jgiBin) = ("", $jgiScr);
	if ($jgiScr =~ /\n/) {
		($jgiActivate, $jgiBin) = ($jgiScr =~ /^(.*\n)(.+)$/s);
	}

	my $covCmd = "";
	my @temporary_bams;
	$covCmd .= "rm -f $out.jgi.*\n";
	if ($isCram){
		die "jgi_depth_cmd:::No reference Fasta given for @dirSS\n" if ($refFA eq "");
		# Convert all CRAMs to temp BAMs before activation (samtools is in the base env, not MF4binners)
		my @splSS = split /\s/,$comBAM;
		$comBAM="";
		for (my $i=0;$i<@splSS;$i++){
			next if ($splSS[$i] eq "");
			my $tmpBam = "$out.jgi.tmp.$i.bam";
			$covCmd .= "$smtBin view -T $refFA -@ $numCores -b $splSS[$i] > $tmpBam\n";
			$comBAM .= "$tmpBam ";
			push @temporary_bams, $tmpBam;
		}
	}
	$covCmd .= $jgiActivate; # conda activation after samtools, before jgi
	$covCmd .= $jgiBin;
	#--pairedContigs $out.jgi.pairs.sparse
	$covCmd .= " --outputDepth $out.jgi.depth.txt  --percentIdentity $perID  $comBAM\n";
	$covCmd .= "test -s $out.jgi.depth.txt\n";
	$covCmd .= "rm -f ".join(" ", @temporary_bams)."\n" if (@temporary_bams);
	#$covCmd .= "gzip $out.jgi*\n";
	if (-s "$out.jgi.depth.txt"){$covCmd="";}

	#$covCmd .= "gzip $nxtBAM.jgi*\n";
	return $covCmd;
}

#2nd: arrray of files, paired sep by ","
sub createGapFillopt($ $ $){
 my ($ofile, $arFiles, $insertSizAr) = @_;
 my @Files = @{$arFiles};
 my @insertSiz = @{$insertSizAr};
 my $opt = "";
 for (my $cnt=0;$cnt<@Files; $cnt++){
	my $line = "LIB$cnt bwa ";
	my $lineMode = 2; 
	if ($lineMode==2){ #PE; only available option at the moment
		my @curFils = split(",",$Files[$cnt]);
		die ("Only paired reads accepatble to GapFiler.\n") if (@curFils != 2);
		$line .= $curFils[0]." ".$curFils[1];
	}
	$line .= " $insertSiz[$cnt] 0.3 FR\n";
	$opt .= $line;
 }
 #print $opt."\n";
	 open O,">",$ofile or die "Cannot write GapFiller options $ofile: $!\n";
	 print O $opt;
	 close O or die "Cannot close GapFiller options $ofile: $!\n";
}
