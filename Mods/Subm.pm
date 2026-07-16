package Mods::Subm;

use warnings;
use strict;
#use List::MoreUtils 'first_index'; 
use Mods::IO_Tamoc_progs qw(getProgPaths convert2Gb);
use Mods::WorkflowControl qw(normalise_job_dependencies);


use Exporter qw(import);
our @EXPORT_OK = qw( findQsubSys emptyQsubOpt qsubSystem qsubSystem2 qsubSystemJobAlive
		qsubSystemWaitMaxJobs MFnext add2SampleDeps numUserJobs);

my $FAILED_SUBMISSION_DEPENDENCY = '__MF4_SUBMISSION_FAILED__';

sub _continue_after_submission_failure {
	my ($optHR, $message) = @_;
	return 0 unless ($optHR->{continueOnSubmitError});
	$optHR->{submissionErrors} = []
		unless (ref($optHR->{submissionErrors}) eq 'ARRAY');
	push @{$optHR->{submissionErrors}}, $message;
	warn "$message\nMATAFILER will skip dependent jobs and continue with later work.\n";
	return 1;
}





sub randStr($){ #will be prefixed to jobname, to make jobs unique to each MF run
	my ($len) = @_;
	my @letters=('A'..'Z','a'..'z',1..9);
	my @letters2=('A'..'Z','a'..'z');
	my $total=scalar(@letters);
	my $newletter ="";
	$newletter = $letters2[int(rand scalar(@letters2))];
	for (my $i=1;$i<$len;$i++){
		$newletter .= $letters[rand $total];
	}
	return $newletter;
}



sub qsubSystem($ $ $ $ $ $ $ $ $ $){
	#args: 1[file to save bash & error & output] 2[actual bash cmd] 3[cores reserved for job]
	# 4["1G": Ram usage per core in GB] 5[0/1: synchronous execution] 6[name of job] 
	# 7[name of job dependencies, separated by ";"]
	# 8[0/1: excute in cwd?] 9[0/1: return qsub cmd or submit job to cluster]
	# Falk Hildebrand, may 2015
	my ($tmpsh,$cmd,$ncores,$memory,$jname,$waitJID,$cwd,$immSubm, $restrHostsAR, $optHR) = @_;
	#$doSync, 5th arg
	#14,12G
	#die $tmpsh."\n";
	#my $jname = $tmpsh;
	#$jname =~ s/.*\///g;$jname =~ s/\.sh$//g;
	#\n#\$ -N $tmpsh
	return("") if ($cmd eq "");
	my $LSF = 0;
	my $qbin = "qsub";
	my $xtra = "";
	my $rTag = $optHR->{rTag};
	my $qmode = $optHR->{qmode};
	
	my $tmpScratchTag = $optHR->{tmpSpaceTag};
	#my $xtraNodeCmds = $optHR->{xtraNodeCmds};
	my $submissionConfig = $optHR->{submissionConfig};
	my @constrains = @{$optHR->{constraint}};# #SBATCH --constraint=
	#die "@constrains";
	my $lockFile = $optHR->{LOCKfile};
	my $nthreads= $ncores;
	if ($ncores =~ m/,/){my @spl = split /,/,$ncores;$ncores = $spl[1]; $nthreads=$spl[0];}
	#different format for bsub and slurm
	if ($memory =~ m/^[\.\d]+$/){$memory  = int($memory+0.5);}
	if ($memory =~ m/^0G$/){$memory  = "1G";} #most likely a rounding error from too many cores..
	if ($memory =~ s/G$//){$memory = int( ($memory* 1024 ) +0.5);};
	my $tmpSpace = convert2Gb( $optHR->{tmpSpace} );
	#die " $optHR->{tmpSpace}   $tmpSpace\n";
	#my $tmpSpace2 = $optHR->{tmpMinG};
	#my $wcKeysForJob = $optHR->{wcKeysForJob};
	my $exclNodes = $optHR->{excludeNodes};
	
	
	#die ($memory."\n");
	#my $queues = "\"".$optHR->{shortQueue}."\"";#"\"medium_priority\"";
	my $queues = "\"".$optHR->{medQueue}."\"";#"\"medium_priority\"";
	my $time = $optHR->{medTime};#"24:00:00";
	if ($optHR->{useHiMemQueue} == 1){
		$queues = "\"".$optHR->{highMemQueue}."\"";$optHR->{useHiMemQueue}=0;
	} elsif ($optHR->{useLongQueue} ==1){
		$queues = "\"".$optHR->{longQueue}."\"";#"\"medium_priority\"";
		#$time = "335:00:00";
		$optHR->{useLongQueue}=0;
	} elsif ($optHR->{useGPUQueue} ==1){
		$queues = "\"".$optHR->{gpuQueue}."\"";#"\"medium_priority\"";
		#$time = "23:00:00";
		$optHR->{useGPUQueue}=0;
	} elsif (defined $optHR->{useNetQueue} && $optHR->{useNetQueue} ==1){
		$queues = "\"".$optHR->{netQueue}."\"";
		$time = $optHR->{longTime};
		$optHR->{useNetQueue}=0;
	} elsif ($optHR->{useShortQueue} ==1){
		$queues = "\"".$optHR->{shortQueue}."\"";#"\"medium_priority\"";
		#$time = "00:45:00";
		$optHR->{useShortQueue}=0;
	}
	$waitJID = normalise_job_dependencies($waitJID);
	my @jspl = split /;/, $waitJID;
	my $has_failed_dependency = grep { $_ eq $FAILED_SUBMISSION_DEPENDENCY } @jspl;
	@jspl = grep { $_ ne $FAILED_SUBMISSION_DEPENDENCY } @jspl;
	$waitJID = join(';', @jspl);

	if ($cwd ne "" && !-d $cwd){system "mkdir -p $cwd";}
	#if ($memory > 250001){$queues = "\"scb\"";}
	$tmpsh =~ m/^(.*\/)[^\/]+$/;
	system "mkdir -p $1" unless (-d $1);
	open O,">",$tmpsh or die "Can't open qsub bash script $tmpsh\n";
	#die "$cmd\n";
	#print "$memory   $queues\n";
	#if (`hostname` !~ m/submaster/){
	if ($qmode eq "slurm"){$LSF = 2;$qbin="sbatch";
		#if ($memory > 250001){$queues = "\"bigmem\"";}
		##SBATCH --cpus-per-task=$ncores\n
		print O "#!/bin/bash\n#SBATCH -N 1\n#SBATCH --cpus-per-task=$ncores\n#SBATCH -o $tmpsh.otxt\n"; #\n#SBATCH -n  $ncores
		
		if ($nthreads != $ncores ){print O "#SBATCH --threads-per-core=1\n#SBATCH --hint=compute_bound\n";} #  specifically for iqtree/raxml
		print O "#SBATCH -e $tmpsh.etxt\n#SBATCH --mem=$memory\n#SBATCH --export=ALL\n";
		#print O "#SBATCH --kill-on-invalid-dep=yes\n";
		#print O "#SBATCH --tmp=$tmpSpace\n" if ($tmpSpace>0);#SBATCH --gres=ssd\n
		foreach my $subTerm ( split /;/, $submissionConfig){
			print O "#SBATCH $subTerm\n" if ($submissionConfig ne "");
		}
		if ($tmpSpace>0 && $tmpScratchTag ne ""){
			print O "#SBATCH $tmpScratchTag". int($tmpSpace+0.5) ."\n" ;
		}
		#"#SBATCH --gres=ssd"
		print O "#SBATCH -p $queues\n";
		print O "#SBATCH --gres=gpu:".$optHR->{gpuCount}."\n" if ($optHR->{gpuCount} > 0);
		#print O "#SBATCH --gres=tmp:${tmpSpace2}G\n" if ($tmpSpace2>0); #50g
		print O "#SBATCH --time=$time\n" unless ($time eq "");
		print O "#SBATCH --exclude=$exclNodes\n" unless ($exclNodes eq "");
		#print O "#SBATCH --localscratch=ssd:50\n"; #for EI cluster
		print O "#SBATCH --chdir=$cwd\n" if ($cwd ne "");
		print O "#SBATCH -J $rTag$jname\n" if ($jname ne "");
		print O "#SBATCH --wc=". $optHR->{wcKeysForJob} . "\n" if ($optHR->{wcKeysForJob} ne "");
		if (@constrains){
			print O "#SBATCH --constraint=". join(",",@constrains) ."\n" if (@constrains);
		}
		#foreach (@constrains){
	#		print O "#SBATCH --constraint=$_\n" if ($_ ne "");
		#}
		if (@jspl > 0) {
			for (@jspl) {s/^\Q$rTag\E//;}
			
			#$xtra .= "--dependency=afterok:".join(":",@jspl)." " if (@jspl > 0);
			if ($optHR->{afterAny}){
				print O "#SBATCH --dependency=afterany:".join(":",@jspl)."\n" ;
			} else {
				print O "#SBATCH --dependency=afterok:".join(":",@jspl)."\n" ;
			}
			#use this one for now, as slurm currently faults without a reason..
			#$xtra .= "--dependency=afterany:".join(":",@jspl)." " if (@jspl > 0);
		}

		#print O "#\$ -S /bin/bash\n#\$ -v LD_LIBRARY_PATH=".$optHR->{cpplib}."\n";##\$ -v TMPDIR=/dev/shm\n";
		#print O "#\$ -v PERL5LIB=".$optHR->{perl5lib}."\n"; #causes problems..
	} elsif ($qmode eq "bash"){
		$qbin="bash";$LSF=3;
		print O "#!/bin/bash\n";
	} elsif ($qmode eq "sge"){
		print O "#!/bin/bash\n#\$ -S /bin/bash\n#\$ -cwd\n#\$ -pe ".$optHR->{qsubPEenv}." $nthreads\n#\$ -o $tmpsh.otxt\n#\$ -e $tmpsh.etxt\n#\$ -l h_rss=$memory\n";#h_vmem=$mem\n";
		print O "#\$ -v LD_LIBRARY_PATH=".$optHR->{cpplib}."\n";#\$ -v TMPDIR=/dev/shm\n";
#		print O "#\$ -v PERL5LIB=".$optHR->{perl5lib}."\n";
		print O "#\$ -V\n";
	} else {
		$LSF = 1;$qbin="bsub";
		print O "#!/bin/bash\n";
		print O "export LD_LIBRARY_PATH=/g/bork3/home/hildebra/env/env/miniconda/lib/:/g/bork3/home/hildebra/env/zlib-1.2.8/:/g/bork8/costea/boost_1_53_0/:/shared/ibm/platform_lsf/9.1/linux2.6-glibc2.3-x86_64/lib:/g/bork3/x86_64/lib64:/g/bork3/x86_64/lib:\${LD_LIBRARY_PATH}\n\n";
		#print O "export LD_LIBRARY_PATH=/g/bork3/home/hildebra/env/zlib-1.2.8:/g/bork3/x86_64/lib64:/lib:/lib64:/usr/lib64:\${LD_LIBRARY_PATH}\n\n";#:/g/software/linux/pack/python-2.7/lib/\nexport PATH=/g/bork3/home/zeller/py-virtualenvs/py2.7_bio1/bin/:\${PATH}\n\n";
		##BSUB -n $ncores\n#BSUB -o $tmpsh.otxt\n#BSUB -e $tmpsh.etxt\n#BSUB -M $mem\n#\$ -v LD_LIBRARY_PATH=/g/bork3/x86_64/lib64:/lib:/lib64:/usr/lib64\n#\$ -v TMPDIR=/dev/shm\n#BSUB -q medium_priority\n";
		my @restrHosts = @{$restrHostsAR};
		if ( @restrHosts > 0){
			$xtra .= " -m \"".join(" ",@restrHosts)."\" ";
			$queues = "\"medium_priority scb\"";
		}
		$xtra .= "-n $nthreads -oo $tmpsh.otxt -eo $tmpsh.etxt -q $queues -M $memory -R \"select[(mem>=$memory)] ";
		$xtra .= "rusage[tmp=$tmpSpace] " if ($tmpSpace>0);
		$xtra .= "span[hosts=1]\" -R \"rusage[mem=$memory]\" "; #
	}
	#set abortion on program fails
	print O "echo \$HOSTNAME;\n";
	print O "set -eo pipefail\n";
	print O "ulimit -c 0;\n";
	#any xtra commands (like module load perl?)
	print O "$optHR->{xtraNodeCmds}\n";
	#prevent core dump files
	#file location check availability
	#print O $optHR->{LocationCheckStrg};

	print O $cmd."\n";
	close O;
	#sleep (1);
	my $depSet=0;
	if ($LSF==2){#slurm
		if ($optHR->{doSync} == 1){$qbin = "srun";}
		
	} elsif ($LSF==3){ #bash
		$xtra = "";
	} elsif ($LSF==1){ #bsub #-M memLimit; -q queueName;  -m "host_name[@cluster_name]; -n minProcessors; 
		if ($optHR->{doSync} == 1){$xtra.="-K ";}
		if ($jname ne ""){$xtra.="-J $rTag$jname ";}
		if (@jspl > 0) {
			my @jspl = split(";",$waitJID);
			#remove empty elements
			@jspl = grep /\S/, @jspl;
			for (@jspl) { s/^\Q$rTag\E//; }
			if (@jspl > 0 ){
				$waitJID = join(") && done(",@jspl);
				$xtra.="-w \"done($waitJID)\" ";
			}
		}
		$tmpsh = " < ".$tmpsh;
	} else{ #qsub
		if ($optHR->{doSync} == 1){$xtra.="-sync y ";}
		if ($jname ne ""){$xtra.="-N $rTag$jname ";}
		if (@jspl > 0) {
			for (@jspl) { s/^\Q$rTag\E//; }
			if (@jspl > 0 ){$xtra.="-hold_jid ".join(",",@jspl) ." ";}
		}
			#$waitJID =~ s/;/,/g;$xtra.="-hold_jid $waitJID ";}
	}
	if ($cwd ne ""){if ($LSF==1) {$xtra.="-cwd $cwd"; }  elsif ($LSF == 0) {$xtra.="-wd $cwd";} }
	my $qcm = "$qbin $xtra $tmpsh \n";
	my $LOGhandle = "";
	if (exists $optHR->{LOG}){ $LOGhandle = $optHR->{LOG};}
	#if (@restrHosts > 0){die $qcm;}
	if ($optHR->{doSubmit} != 0 && $immSubm){
		if ($has_failed_dependency) {
			my $message = "Skipping submission for $tmpsh because an upstream submission failed";
			return ($FAILED_SUBMISSION_DEPENDENCY, $qcm)
				if (_continue_after_submission_failure($optHR, $message));
			die "$message\n";
		}
		system "rm -f $tmpsh.otxt $tmpsh.etxt";
		print $LOGhandle $qcm."\n" unless ($LOGhandle eq "" || !defined($LOGhandle) );
		#print("$qcm\n\n");
		print "SUB:$jname\t";
		#actual job excecution!
		my $ret = `$qcm`;
		my $submit_status = $?;
		if ($submit_status != 0) {
			my $exit_code = $submit_status == -1 ? -1 : ($submit_status >> 8);
			my $message = "Job submission failed (exit $exit_code): $qcm$ret";
			return ($FAILED_SUBMISSION_DEPENDENCY, $qcm)
				if (_continue_after_submission_failure($optHR, $message));
			die $message;
		}
		if ($LSF == 2){#slurm get jobid
			chomp $ret;
			die "Could not parse Slurm job id from submission output: $ret\n"
				unless ($ret =~ /^Submitted batch job (\d+)\s*$/);
			$jname=$1;
		} elsif ($LSF == 0) {
			die "Could not parse SGE job id from submission output: $ret\n"
				unless ($ret =~ /\bYour job(?:-array)?\s+(\d+)\b/);
			$jname=$1;
		} elsif ($LSF == 1) {
			die "Could not parse LSF job id from submission output: $ret\n"
				unless ($ret =~ /\bJob <(\d+)>/);
			$jname=$1;
		}
		# Only record a lock after the scheduler has accepted the job.
		if ($lockFile ne "" && !-e $lockFile){
			open my $lock, ">", $lockFile or die "Cannot create lock $lockFile: $!\n";
			close $lock or die "Cannot close lock $lockFile: $!\n";
		}
	}
	
	#die "$qcm\n";
	my $retJName = "$rTag$jname"; $retJName = "" if (!$immSubm); #return empty (for slurm), since no fwd job predictions..

	return ($retJName,$qcm);
}



sub numPendingJobs($){
	my ($optHR) = @_;
	my $qmode = "slurm"; $qmode = $optHR->{qmode} if (defined($optHR->{qmode}));
	my $srchCmd="" ;#= "squeue -u \$USER  -t PENDING | wc -l";
	my $num = 0;
	if ($qmode eq "slurm"){
		$srchCmd = "squeue -h -u \$USER -t PENDING | wc -l";
	} elsif ($qmode eq "sge"){
		$srchCmd = "qstat | grep \$USER  | wc -l";
		die "Subm.pm::numPendingJobs() not implemented for sge!\n";
	} elsif ($qmode eq "bash"){
		return 0;
	} else {$srchCmd="bjobs -p -noheader | wc -l";
		die "Subm.pm::numPendingJobs() not implemented for bsub!\n";
	}
	$num = `$srchCmd`; chomp $num;
	die "Failed to count pending jobs with: $srchCmd\n" if ($? != 0);
	return $num;
}
sub numUserJobs{
	my ($optHR) = $_[0];
	my $rmSelf=0; $rmSelf = $_[1] if (@_>1);
	my $qmode = "slurm"; $qmode = $optHR->{qmode} if (defined($optHR->{qmode}));
	my $srchCmd ="";#= "squeue -u \$USER   | wc -l";
	if ($qmode eq "slurm"){
		$srchCmd = "squeue -h -u \$USER | wc -l";
	} elsif ($qmode eq "sge"){
		$srchCmd = "qstat | grep \$USER  | wc -l";
	} elsif ($qmode eq "bash"){
		return 0;
	} else {$srchCmd="bjobs -noheader | wc -l";
	}
	my $num = 0;
	$num = `$srchCmd`; chomp $num;
	die "Failed to count user jobs with: $srchCmd\n" if ($? != 0);


	if ($rmSelf && $qmode eq "slurm"){
		my $SjobID = `echo \$SLURM_JOBID`; chomp $SjobID;
		#print "\"$SjobID\"\n";
		if ($SjobID ne ""){$num --;}
	}

	return $num;
}


sub findQsubSys($){
	my $iniVal = "";
	$iniVal = $_[0] if (@_ > 0);
	#my $iniVal = "lsf";
	if ($iniVal ne ""){
		$iniVal = lc $iniVal; 
		$iniVal = "lsf" if ($iniVal eq "bsub");
		$iniVal = "sge" if ($iniVal eq "qsub");
		$iniVal = "slurm" if ($iniVal eq "sbatch");
	} else {
		$iniVal = "lsf";
		my $bpath = `which bsub  2>/dev/null`;chomp $bpath;my $bpresent=0; 
		$bpresent=1 if ($bpath !~ m/\n/ && -e $bpath);
		my $spath = `which sbatch  2>/dev/null`;chomp $spath;my $spresent=0; 
		$spresent=1 if ($spath !~ m/\n/ && -e $spath);
		my $qpath = `which qsub  2>/dev/null`; chomp $qpath;
		my $qpresent=0; $qpresent=1 if ($qpath !~ m/\n/ && -e $qpath);
		#print "$qpath\n";
		if ($spresent ){#slurm gets preference
			$iniVal="slurm";
		}elsif (!$bpresent && $qpresent){
			$iniVal = "sge";
		}elsif (!$qpresent && !$bpresent && !$spresent){
			die "No queueing system found (sbatch, qsub, or bsub). Use -qsubSystem bash for local execution.\n";
		}
	print "Using qsubsystem: $iniVal\n";
	}
	#die;
	return $iniVal;
}
sub emptyQsubOpt{
	my ($doSubm) = $_[0];
	my $locChkStr = $_[1];
	my $qmode = "";
	$qmode = $_[2] if (@_ > 2);
	
	if (@_ > 2){$qmode = $_[2];}
	$qmode = findQsubSys($qmode);
	die "qsub system mode has to be \'lsf\', \'bash\', \'slurm\' or \'sge\'!\n" if ($qmode ne "lsf" &&$qmode ne "slurm" && $qmode ne "sge"&& $qmode ne "bash");
	my $MFdir = getProgPaths("MFLRDir");
	my $longQ = getProgPaths("longQueue",0); my $shortQ =  getProgPaths("shortQueue",0); my $medQ = getProgPaths("mediumQueue",1);
	#die "$shortQ\n";
	my $gpuQ = getProgPaths("gpuQueue",0);
	my $netQ = getProgPaths("netQueue",0);
	my $himemQ = getProgPaths("highMemQueue",0);
	if ($longQ eq ""){$longQ =  $medQ;}
	if ($medQ eq "" ){die "FATAL: no medium queue defined!\n";};
	if ($gpuQ eq "" ){$gpuQ = $medQ;};
	if ($netQ eq "" ){$netQ = $medQ;};
	if ($himemQ eq "" ){$himemQ = $medQ;};
	if ($shortQ eq "" ){$shortQ = $medQ;};
	my $xtraNodeCmds = getProgPaths("subXtraCmd",0);
	$xtraNodeCmds = "" unless (defined $xtraNodeCmds);
	my $medTime = getProgPaths("medTime",0);	my $shortTime = getProgPaths("shortTime",0);
	my $longTime = getProgPaths("longTime",0);
	my $subConfig = getProgPaths("submissionConfig",0);
	my @constr = ();
	if ($subConfig =~ s/--constraint=(\S+)//){
		#print "!!! $1\n";
		push(@constr, $1);
	}
	chomp($subConfig);
	@constr = grep(/\S/, @constr);
	#die "@constr\n$subConfig\nYW\n";
	
	#if ($qmode eq "slurm"){$shortQ = "htc"; $longQ="htc";}#$shortQ = "1day"; $longQ="1month";}
	my %ret = (
		rTag => randStr(3),
		doSubmit => $doSubm,
		LocationCheckStrg => $locChkStr,
		doSync => 0,
		longQueue => $longQ,
		gpuQueue => $gpuQ,
		netQueue => $netQ,
		highMemQueue => $himemQ,
		longTime => $longTime,#7days
		medQueue => $medQ,
		medTime => $medTime,#"24:00:00",
		shortQueue => $shortQ,
		shortTime => $shortTime, #2hrs
		useLongQueue => 0,
		useGPUQueue => 0,
		useNetQueue => 0,
		gpuCount => 0,
		useShortQueue => 0,
		useHiMemQueue => 0,
		submissionConfig => $subConfig,
		constraint => \@constr,
		qsubPEenv => getProgPaths("qsubPEenv"),
		perl5lib => "$MFdir:\$PERL5LIB",
		cpplib => "",
		tmpSpace => 15, #default was 15G; unit is G
		tmpSpaceTag => getProgPaths("nodeTmpDirTAG",0),
		LOCKfile => "",
		#tmpMinG => 10,
		afterAny => 0,
		excludeNodes => "",
		xtraNodeCmds => $xtraNodeCmds,
		qmode => $qmode,
		wcKeysForJob => "",
		#LOG => undef,
	);
	#die "$MFdir\n";
	return \%ret;
}

sub qsubSystemJobAlive{
	my ($jAr,$optHR) = @_;
	my $killFailedJobs=0;
	$killFailedJobs = $_[2] if (@_ > 2);
	my @jobs = split /;/, normalise_job_dependencies($jAr);
	@jobs = grep { $_ ne $FAILED_SUBMISSION_DEPENDENCY } @jobs;
	return unless (@jobs);
	
	
	my $qmode = $optHR->{qmode};
	my $cmd1="";
	my $rTag = $optHR->{rTag};

	for (@jobs) {s/^\Q$rTag\E//;}
	if ($qmode eq "slurm"){
		$cmd1 = "squeue -h -u \$USER -o '%i'";
	} elsif ($qmode eq "sge"){
		$cmd1 = "qstat -u \$USER | awk 'NR > 2 {print \$1}'"
	} elsif ($qmode eq "bash"){
		return;
	} else {$cmd1="bjobs -noheader -o jobid";
	}
	my %wanted = map { $_ => 1 } @jobs;
	my $announced = 0;
	while (1) {
		my $output = `$cmd1`;
		die "Failed to query active jobs with: $cmd1\n" if ($? != 0);
		my %active = map { $_ => 1 }
			grep { /^\d+$/ }
			map { my $id = $_; $id =~ s/^\s+|\s+$//g; $id }
			split /\n/, $output;
		my @remaining = grep { $active{$_} } keys %wanted;
		last unless (@remaining);
		print "Waiting for ".scalar(@remaining)."/".scalar(@jobs)." jobs to finish\n"
			unless ($announced++);
		if ($killFailedJobs){
			my $killed = qsubDepNeverKill();
			print " Killed $killed jobs with Dependency never completed\n" if ($killed > 0);
		}
		sleep (60);
	}
	#print "returning\n";
	return;
}

sub qsubDepNeverKill{
	my $srchCmd = "squeue -u \$USER -t PENDING -o \"\%8i \%.15R \%17E\"  | grep 'ependencyNe' | cut -f1 -d' ' | xargs  -t -i scancel {} | wc -l ";
	my $num = 0;
	$num = `$srchCmd`; chomp $num;
	return $num;
	
}


sub qsubSystemWaitMaxJobs{
	my ($checkMaxNumJobs) = @_;
	my $killPend = $_[1] if (@_ > 1);
	my $optHR = {}; $optHR = $_[2] if (@_ > 2);
	
	return if ($checkMaxNumJobs <= 0);
	#my $srchCmd = "squeue |grep \$USER | grep PD |wc -l";
	my $num = numPendingJobs($optHR);
	my $waitCnt = 0;
	while ($num > $checkMaxNumJobs){
		if ($killPend){
			my $killed = qsubDepNeverKill();
			print " Killed $killed jobs with Dependency never completed\n" if ($killed > 0);
			
			#die;
		}
		print "waiting for jobs to finish (>$checkMaxNumJobs, qsubSystemWaitMaxJobs)...\n" if ($waitCnt == 0);
		sleep(40);
		$num =  numPendingJobs($optHR);;#`$srchCmd`; chomp $num;
		$waitCnt++;
		#print " $num ";
	}
	return;
}

sub qsubSystem2{
	my ($tmpsh,$optHR) = @_;
	my $hxr = {};
	$hxr = $_[2] if (@_ > 2 && defined $_[2]);
	my %xtras = %{$hxr}; 
	my $ncores = 0; 
	if (exists($xtras{cores})){$ncores = $xtras{cores};}
	my $nthreads= $ncores;
	if ($ncores =~ m/,/){my @spl = split /,/,$ncores;$ncores = $spl[1]; $nthreads=$spl[0];}
	if ($ncores != 0){#read in file, change it..
		open my $in,"<",$tmpsh or die "qsubSystem2: cant open $tmpsh\n";chomp(my @lines = <$in>); close $in;
		for (my $i=0;$i<@lines;$i++){
			if ($lines[$i] =~ m/--cpus-per-task/ || $lines[$i] =~ m/--mincpus/){
				$lines[$i] = "#SBATCH --cpus-per-task=$ncores";
				my $next_line = $i+1 < @lines ? $lines[$i+1] : "";
				$lines[$i] .= "\n#SBATCH --threads-per-core=1\n#SBATCH --hint=compute_bound" unless ($next_line =~ m/threads-per-core/);
			}
		}
		open my $out,">",$tmpsh or die "qsubSystem2: cant update $tmpsh\n";
		print {$out} join("\n",@lines), "\n";
		close $out or die "qsubSystem2: cant close updated $tmpsh\n";
	}
	my $xtra = "";
	my $qbin = "qsub";
	my $qmode = $optHR->{qmode};
	if ($qmode eq "slurm"){$qbin="sbatch";
	} elsif ($qmode eq "sge"){
	} elsif ($qmode eq "bash"){$qbin="bash";
	} else {$qbin="bsub";$xtra="<";
	}
	my $qcm = "$qbin $xtra $tmpsh \n";
	system($qcm) == 0 or die "qsubSystem2 failed: $qcm";
	return $qcm;
}



#handles deleting of lock file, if all jobs have finished for current sample
sub MFnext($ $ $ $){
	my ($lckFile,$aR,$Jnum,$QSBoptHR) = @_;
	return if (! @{$aR});
	my $logF = $lckFile; $logF =~ s/\/[^\/]+$/rmLock.sh/;
	my $cmd = "echo \"all smpl associated jobs seem to have quit. Releasing lock..\"\nrm -f $lckFile\n";
	#my @jobs = @{$aR};
	my $jDepe = normalise_job_dependencies($aR);
	my $jobN = "RMLCK$Jnum";
	#print "$logF\n$jDepe\n\n"; 
	$QSBoptHR->{afterAny}=1;
	my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
	$QSBoptHR->{useShortQueue} =1;
	my ($jN,$submitCommand) = qsubSystem($logF,$cmd,1,"1G",$jobN,$jDepe,"",1,[],$QSBoptHR);
	$QSBoptHR->{afterAny}=0;$QSBoptHR->{useShortQueue}=0;
	$QSBoptHR->{tmpSpace} =$tmpSHDD;
	push(@{$aR}, $jN) if ($jN ne "");
}


sub add2SampleDeps($ $){
	my ($ar1, $ar2) = @_;
	@{$ar1} = split /;/, normalise_job_dependencies($ar1, $ar2);
}
