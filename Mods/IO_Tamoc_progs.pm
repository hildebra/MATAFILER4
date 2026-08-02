package Mods::IO_Tamoc_progs;
use warnings;
use Cwd 'abs_path';
use strict;
use Mods::ReadLibrary qw(legacyLibraryArrays);

use vars qw($CONFIG_FILE @CONFIG_TEXT %CONFIG_HASH $CONFIG_LOADED);
$CONFIG_FILE="";
@CONFIG_TEXT = ();
%CONFIG_HASH = ();
$CONFIG_LOADED = 0;
sub setConfigFile;

#TAMOC programs related to IO to other programs, program paths .. not real subroutines that do anything

use Exporter qw(import);
our @EXPORT_OK = qw(getProgPaths truePath
					inputFmtSpades inputFmtSpadesLibraries inputFmtMegahit inputFmtMegahitLibraries inputFmtMegahitRuntimeLibraries jgi_depth_cmd createGapFillopt setConfigFile
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
		if (-d $DDI || $DDI =~ m/\/$/){
			$DDI =~ s{/$}{};
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

sub _bundledConfigFile{
	my ($fileName) = @_;
	my $modDir = $INC{"Mods/IO_Tamoc_progs.pm"} || __FILE__;
	$modDir =~ s{IO_Tamoc_progs\.pm$}{};
	return $modDir.$fileName;
}

sub setConfigFile{
	my @var = @_;
	my $customCfg = 0;
	my $newConfig;
	if (@var == 1 && $var[0] eq "internal"){
		$newConfig = _bundledConfigFile("config_internal.txt");
	} elsif (@var == 1 && $var[0] eq "DBconfig"){
		$newConfig = _bundledConfigFile("config_DBs.txt");
	} elsif (@var == 1 && $var[0] ne ""){
		$newConfig = $var[0];
		$customCfg = 1;
	} else {#default value
		$newConfig = _bundledConfigFile("MATAFILERcfg.txt");
	}
	$CONFIG_FILE = $newConfig;
	die "Can't find MATAFILER config file: $CONFIG_FILE\nConsider changing path to config file via \"-config\" argument.\n Aborting..\n" unless (-e $CONFIG_FILE);
	# An explicit selection starts a fresh configuration generation. This matters
	# to test harnesses and long-lived callers that select another site config
	# after an earlier lookup.
	@CONFIG_TEXT = ();
	%CONFIG_HASH = ();
	$CONFIG_LOADED = 0;
	print "Using config file : $CONFIG_FILE\n" if ($customCfg);
}

sub truePath{
	my ($TMCpath, $enforce) = @_;
	$enforce = 0 unless defined $enforce;
	return $TMCpath unless ($enforce || $TMCpath =~ /^\$/);

	# Expand every $NAME or ${NAME} occurrence. Historically only the first
	# occurrence was expanded and braced names were interpreted incorrectly.
	$TMCpath =~ s{
		\$(?:\{([A-Za-z_][A-Za-z0-9_]*)\}|([A-Za-z_][A-Za-z0-9_]*))
	}{
		my $envName = defined($1) ? $1 : $2;
		die "Environment variable \$$envName used in path '$TMCpath' is not set\n"
			unless (exists($ENV{$envName}) && defined($ENV{$envName}) && $ENV{$envName} ne '');
		$ENV{$envName};
	}gex;
	return $TMCpath;
}

sub loadConfigs{
	#loads once in every program run the entire config file(s) into hash %CONFIG_HASH
	return if $CONFIG_LOADED;
	setConfigFile() if ($CONFIG_FILE eq "");
	my $selectedConfig = $CONFIG_FILE;
	print "READING config files \"$selectedConfig\" .. ";
	@CONFIG_TEXT = ();
	%CONFIG_HASH = ();
	my @configFiles = (
		[$selectedConfig, "selected"],
		[_bundledConfigFile("config_internal.txt"), "internal"],
		[_bundledConfigFile("config_DBs.txt"), "database"],
	);
	my %seenConfig;
	for my $configSpec (@configFiles){
		my ($configPath, $configKind) = @{$configSpec};
		next if $seenConfig{$configPath}++;
		open my $configFH, "<", $configPath
			or die "Can\x27t open $configKind config $configPath: $!\n";
		my @configLines = <$configFH>;
		close $configFH or die "Can\x27t close $configKind config $configPath: $!\n";
		chomp @configLines;
		s/\r$// for @configLines;
		push @CONFIG_TEXT, @configLines;
	}
	# Resolve foundational keys before ordinary entries. Source order is retained
	# within each key, so selected-user values still precede bundled defaults,
	# while placeholder expansion no longer depends on line order.
	my @foundationOrder = qw(MFLRDir BINDir DBDir MGSTKDir SINGcmd CONDcmd CONDA CONDAbaseEnv PY3cmd Rscript Rpath);
	my %foundationRank;
	@foundationRank{@foundationOrder} = (0 .. $#foundationOrder);
	my @foundationLines = map { [] } @foundationOrder;
	my @ordinaryLines;
	for my $line (@CONFIG_TEXT){
		my ($key) = $line =~ /^([^\t]+)/;
		if (defined($key) && exists($foundationRank{$key})){
			push @{$foundationLines[$foundationRank{$key}]}, $line;
		} else {
			push @ordinaryLines, $line;
		}
	}
	# Stable buckets keep this initialization linear in the number of config
	# lines. It is paid once; subsequent getProgPaths calls are hash lookups.
	@CONFIG_TEXT = ((map { @{$_} } @foundationLines), @ordinaryLines);
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
			$prePath =~s/\[MFLRDir\]/$TMCpath/g;
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
					die "Your \"CONDA\" seems to be wrongly setup. Ensure this is configured in \"[MATAFILER-dir]/config.txt\" and has a form similar to:\"\nCONDA\t. \$MAMBA_ROOT_PREFIX/etc/profile.d/mamba.sh\n\"\n";
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
			if (@spl == 1) {
				$CONFIG_HASH{$XVar} = "" unless exists($CONFIG_HASH{$XVar});
				next;
			}
			
			if ($l !~ m/^$XVar\t([^#^\t]+)/){next;}
			my $reV = $1;
			
			#die "$reV  $XVar  $l\n";
			die "$reV\n" if (!defined($reV));
			$reV =~ s/\[MFLRDir\]/$TMCpath/g if ($Tset);
			if ($MGSTKDirset){
				$reV =~ s/\[MGSTKDir\]/$MGSTKDir/g;
			}
			$reV =~ s/\[BINDir\]/$BINpath/g if ($Bset);
			$reV =~ s/\[DBDir\]/$DBpath/g if ($DBset);
			$reV =~ s/\[SINGcmd\]/$SINGcmd/g if ($SINGset);
			$reV =~ s/\[PY3\]/$PY3cmd/g if ($PY3set);
			$reV =~ s/\[Rscript\]/$Rscriptcmd/g if ($Rscriptset);
			$reV =~ s/\[Rpath\]/$Rpath/g if ($RpathSet);
			if ($l =~ m/env:([^#^\t]+)/){
				my $tarEnv = $1;
				#$reV = "$CONDA;$CONDcmd activate $1\n$reV";
				$reV = "$CONDA;if [[ \$CONDA_DEFAULT_ENV != $tarEnv ]]; then $CONDcmd activate $tarEnv; fi\n$reV";
			}
			
			#return $reV;
			# Configuration sources are loaded from most to least specific:
			# the selected user config, internal defaults, then database defaults.
			# Preserve the first definition so a user can override shipped values.
			$CONFIG_HASH{$XVar} = $reV unless exists($CONFIG_HASH{$XVar});
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
	$CONFIG_LOADED = 1;
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
	# Parsing and placeholder expansion are paid once per selected config. The
	# many pipeline lookups after this point are ordinary hash reads.
	loadConfigs() unless $CONFIG_LOADED;
	die "Something went wrong loading MATAFILER configs.. aborting\n" unless $CONFIG_LOADED;
	
	
	if (@multVars > 0){
		my @missing = grep { !exists($CONFIG_HASH{$_}) || ($required != 0 && $CONFIG_HASH{$_} eq "") } @multVars;
		die "Can't find configuration for ".join(", ", @missing)." in MATAFILER config ($CONFIG_FILE)\n"
			if $required != 0 && @missing;
		my @retA = map { exists($CONFIG_HASH{$_}) ? $CONFIG_HASH{$_} : "" } @multVars;
		return \@retA;
	}
	if (exists($CONFIG_HASH{$srchVar}) && ($required == 0 || $CONFIG_HASH{$srchVar} ne "")){
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
	$MapperProg2 = decideMapper($MapperProg2, "");
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
			# metaSPAdes accepts one paired-end short-read library. Repeated
			# --pe1 arguments are multiple files belonging to that same library;
			# do not turn coassembly samples into separate libraries.
			my $library = $peTerm eq "--pe" ? 1 : $i + 1;
			$sprds .= " ${peTerm}${library}-1 $p1[$i] ${peTerm}${library}-2 $p2[$i]";
		}
		for (my $i=0;$i<@singl;$i++){
			next if ($singl[$i] eq "");
			my $peTerm = "--pe";$peTerm = "--gemcode" if ($readTec[$i] =~ m/SLR/);
			my $library = $peTerm eq "--pe" ? 1 : $i + 1;
			$sprds .= " ${peTerm}${library}-s $singl[$i]";
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
	my @mapping_files;
	my %seen_mapping;
	my $isCram=0;
	foreach my $DDI (@dirSS){
		if (-f $DDI && $DDI =~ /\.(?:bam|cram)$/i) {
			die "jgi_depth_cmd:::Empty mapping file $DDI\n" unless -s $DDI;
			push @mapping_files, $DDI unless $seen_mapping{$DDI}++;
		} else {
			$DDI =~ s{/$}{};
			my $marker = "$DDI/mapping/done.sto";
			die "jgi_depth_cmd:::Missing mapping marker $marker\n" unless (-s $marker);
			open my $marker_fh, '<', $marker or die "Cannot read $marker: $!\n";
			my $SmplNm = <$marker_fh>;
			close $marker_fh;
			chomp $SmplNm;
			my $named_mapping = "$DDI/mapping/$SmplNm";
			my @candidates = ($named_mapping);
			if ($SmplNm !~ /\.sup-smd\./i) {
				(my $supplemental = $named_mapping) =~ s/-smd\./.sup-smd./i;
				push @candidates, $supplemental if $supplemental ne $named_mapping;
			}
			my $mapping_found = 0;
			for my $candidate (@candidates) {
				if (!-s $candidate && $candidate =~ /\.bam$/i) {
					(my $cram = $candidate) =~ s/\.bam$/.cram/i;
					$candidate = $cram if -s $cram;
				} elsif (!-s $candidate && $candidate =~ /\.cram$/i) {
					(my $bam = $candidate) =~ s/\.cram$/.bam/i;
					$candidate = $bam if -s $bam;
				}
				next unless -s $candidate;
				$mapping_found = 1;
				push @mapping_files, $candidate unless $seen_mapping{$candidate}++;
			}
			die "jgi_depth_cmd:::Can't find a non-empty BAM or CRAM named by $marker\n"
				unless $mapping_found;
		}
	}
	$isCram = scalar grep { /\.cram$/i } @mapping_files;
	#my $comBAM = join("/mapping/Align_ment-smd.bam ",@dirSS);
	# Split conda activation prefix from the binary name so the activation can be
	# emitted before any pipe, keeping just the bare binary after the pipe.
	my ($jgiActivate, $jgiBin) = ("", $jgiScr);
	if ($jgiScr =~ /\n/) {
		($jgiActivate, $jgiBin) = ($jgiScr =~ /^(.*\n)(.+)$/s);
	}

	my $covCmd = "set -e\n";
	my @temporary_bams;
	$covCmd .= "rm -f $out.jgi.*\n";
	if ($isCram){
		die "jgi_depth_cmd:::No reference Fasta given for @dirSS\n" if ($refFA eq "");
		# Convert all CRAMs to temp BAMs before activation (samtools is in the base env, not MF4binners)
		for (my $i=0;$i<@mapping_files;$i++){
			next unless $mapping_files[$i] =~ /\.cram$/i;
			my $tmpBam = "$out.jgi.tmp.$i.bam";
			$covCmd .= "$smtBin view -T $refFA -@ $numCores -b $mapping_files[$i] > $tmpBam\n";
			$mapping_files[$i] = $tmpBam;
			push @temporary_bams, $tmpBam;
		}
	}
	$covCmd .= $jgiActivate; # conda activation after samtools, before jgi
	$covCmd .= $jgiBin;
	#--pairedContigs $out.jgi.pairs.sparse
	$covCmd .= " --outputDepth $out.jgi.depth.txt  --percentIdentity $perID  ".join(' ', @mapping_files)."\n";
	$covCmd .= "test -s $out.jgi.depth.txt\n";
	$covCmd .= "rm -f ".join(" ", @temporary_bams)."\n" if (@temporary_bams);
	#$covCmd .= "gzip $out.jgi*\n";
	if (-s "$out.jgi.depth.txt"){$covCmd="";}

	#$covCmd .= "gzip $nxtBAM.jgi*\n";
	return $covCmd;
}

sub inputFmtSpadesLibraries {
	my ($libraries, $logDir) = @_;
	my ($r1, $r2, $single, $labels, $technologies) = legacyLibraryArrays($libraries, 1);
	return inputFmtSpades($r1, $r2, $single, $logDir, $technologies);
}

sub inputFmtMegahitLibraries {
	my ($libraries, $logDir) = @_;
	my ($r1, $r2, $single) = legacyLibraryArrays($libraries, 1);
	return inputFmtMegahit($r1, $r2, $single, $logDir);
}

sub _shell_quote_megahit_input {
	my ($value) = @_;
	die "Cannot quote an undefined MEGAHIT input\n" unless defined $value;
	die "MEGAHIT input contains a NUL or newline\n" if $value =~ /[\0\r\n]/;
	$value =~ s/'/'"'"'/g;
	return "'$value'";
}

# Clean-read records describe outputs before their cleaning jobs run.  In
# particular, SDM's paired-read orphan file is optional and may be absent in
# older/stale clean directories.  Resolve those optional inputs in the
# assembly job, after its cleaning dependencies have completed, instead of
# baking every projected path into MEGAHIT's -r argument.
sub inputFmtMegahitRuntimeLibraries {
	my ($libraries, $arrayName) = @_;
	$arrayName ||= 'megahit_inputs';
	die "Invalid MEGAHIT shell-array name '$arrayName'\n"
		unless $arrayName =~ /^[A-Za-z_][A-Za-z0-9_]*$/;
	my ($r1, $r2, $single) = legacyLibraryArrays($libraries, 1);
	my (@left, @right, @singletons);
	for (my $i = 0; $i < @{$r1}; $i++) {
		my $left = $r1->[$i] || '';
		my $right = $r2->[$i] || '';
		die "Read library index $i has only one mate for MEGAHIT\n"
			if (($left eq '') != ($right eq ''));
		if ($left ne '') {
			die "MEGAHIT input paths cannot contain commas: $left / $right\n"
				if $left =~ /,/ || $right =~ /,/;
			push @left, $left;
			push @right, $right;
		}
		my $orphan = $single->[$i] || '';
		if ($orphan ne '') {
			die "MEGAHIT input paths cannot contain commas: $orphan\n" if $orphan =~ /,/;
			push @singletons, $orphan;
		}
	}
	die "No read inputs supplied to MEGAHIT\n" unless @left || @singletons;

	my $command = "$arrayName=()\n";
	if (@left) {
		for (my $i = 0; $i < @left; $i++) {
			my $leftQ = _shell_quote_megahit_input($left[$i]);
			my $rightQ = _shell_quote_megahit_input($right[$i]);
			my $messageQ = _shell_quote_megahit_input(
				"Missing or empty paired input for MEGAHIT: $left[$i] / $right[$i]"
			);
			$command .= "if [[ ! -s $leftQ || ! -s $rightQ ]]; then "
				."printf '%s\\n' $messageQ >&2; exit 41; fi\n";
		}
		my $leftCSV = _shell_quote_megahit_input(join(',', @left));
		my $rightCSV = _shell_quote_megahit_input(join(',', @right));
		$command .= "$arrayName+=( -1 $leftCSV -2 $rightCSV )\n";
	}

	if (@singletons) {
		my $singleArray = "${arrayName}_singletons";
		my $singleCSV = "${arrayName}_singleton_csv";
		$command .= "$singleArray=()\n";
		for my $orphan (@singletons) {
			my $orphanQ = _shell_quote_megahit_input($orphan);
			my $messageQ = _shell_quote_megahit_input(
				"Skipping missing or empty optional MEGAHIT singleton: $orphan"
			);
			$command .= "if [[ -s $orphanQ ]]; then $singleArray+=( $orphanQ ); "
				."else printf '%s\\n' $messageQ >&2; fi\n";
		}
		my $singleLength = '${#'.$singleArray.'[@]}';
		my $singleValues = '${'.$singleArray.'[@]}';
		$command .= "if (( $singleLength )); then "
			."printf -v $singleCSV '%s,' \"$singleValues\"; "
			."$singleCSV=\${$singleCSV%,}; $arrayName+=( -r \"\$$singleCSV\" ); fi\n";
	}
	my $argumentCount = '${#'.$arrayName.'[@]}';
	$command .= "if (( $argumentCount == 0 )); then printf '%s\\n' "
		."'No non-empty read inputs remain for MEGAHIT' >&2; exit 42; fi\n";
	return ($command, '"${'.$arrayName.'[@]}"');
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
