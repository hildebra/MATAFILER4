package Mods::GenoMetaAss;
use warnings;
#use Cwd 'abs_path';
use strict;
use Fcntl qw(S_ISDIR S_ISREG);
use IO::Compress::Gzip ();
use IO::Uncompress::Gunzip ();
use Time::HiRes ();
#use List::MoreUtils 'first_index'; 
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::ReadLibrary qw(
	cloneReadLibraries ensureSeqSetLibraries ensureCleanSeqSetLibraries
	syncCleanSeqSetLegacy readLibrariesByScope legacyLibraryArrays
);

use Exporter qw(import);
our @EXPORT_OK = qw(
		gzipwrite gzipopen lcp prefix_find
		fileGZe fileGZs filsizeMB resolveExistingFile contig_stats_coverage_complete
		
		readMap 
		systemW
		
		
		readMapS getDirsPerAssmblGrp checkSeqTech is3rdGenSeqTech 
		parseSupportReads normaliseSupportReads discoverReadFiles
		resetAsGrps checkAssmblGrp
		getRawSeqsAssmGrp getCleanSeqsAssmGrp 
		getRawLibrariesAssmGrp getCleanLibrariesAssmGrp
		addFileLocs2AssmGrp iniCleanSeqSetHR hasSuppRds
		
		getAssemblPath getAssemblGFF getAssemblContigs
		
		renameFastHD  prefixFAhd parse_duration resolve_path
		clenSplitFastas 
		
		readClstrRev  readClstrRevGenes readClstrRevContigSubset readClstrRevSmplCtgGenSubset
		unzipFileARezip  is_integer 
		readGFF reverse_complement reverse_complement_IUPAC
		
		readFasta
		writeFasta readFastHD splitFastas  renameFastaCnts renameFastqCnts
		
		readTabByKey convertNT2AA runDiamond median mean quantile
		 );#Binning Related

my %readDirectoryCache;



sub parse_duration {
	my $seconds = shift;
	my $hours = int( $seconds / (60*60) );
	my $mins = ( $seconds / 60 ) % 60;
	my $secs = $seconds % 60;
	return sprintf("00:00:%02d", $seconds) if $seconds < 60;
	return sprintf("%02d:%02d:%02d", $hours,$mins,$secs);
}

#sums up filsizes for 1 or several files; first arg can be path to files, all subsequent args need to be filenames
sub filsizeMB{
	my $totalMapSize=0;
	if (!@_){return  $totalMapSize;}
	my $path="";
	if ($_[0] ne "" && -d $_[0]){$path=shift @_;}
	foreach my $fh (@_){
		next unless (defined($fh) && $fh ne "");
		my $file = $fh;
		$file = "$path/$fh" if ($path ne "" && !_path_is_absolute($fh));
		my @fileStat = stat($file);
		$totalMapSize += $fileStat[7] / (1024 * 1024)
			if (@fileStat && S_ISREG($fileStat[2]));
	}
	return $totalMapSize;
}

sub _path_is_absolute {
	my ($path) = @_;
	return defined($path) && $path =~ m{^(?:/|[A-Za-z]:[\\/]|\\\\)};
}

sub _join_path {
	my ($base, $path) = @_;
	return $path if (!defined($base) || $base eq "" || _path_is_absolute($path));
	$base =~ s{[\\/]+$}{};
	$path =~ s{^[\\/]+}{};
	return "$base/$path";
}

# SupportReads historically accepted both "TECH:file1,file2" and
# "TECH:file1;TECH:file2". Parse the technology once, preserve colons in paths
# (notably Windows drive letters), and anchor relative paths at #DirPath.
sub parseSupportReads {
	my ($spec, $base_dir) = @_;
	return ("", []) unless (defined($spec) && $spec =~ /\S/);
	$base_dir = "" unless defined($base_dir);

	my $technology = "";
	my @paths;
	foreach my $case (split /;/, $spec) {
		$case =~ s/^\s+|\s+$//g;
		next if ($case eq "");

		my $payload = $case;
		if ($case =~ /^([^:,\s]+):(.*)$/s) {
			my ($case_technology, $case_payload) = ($1, $2);
			if ($technology ne "" && $case_technology ne $technology) {
				die "SupportReads mixes technologies '$technology' and '$case_technology' in '$spec'. Use one support technology per sample.\n";
			}
			$technology = $case_technology;
			$payload = $case_payload;
		} elsif ($technology eq "") {
			die "SupportReads entry '$case' has no sequencing-technology prefix (for example PB:path/to/reads.fq.gz).\n";
		}

		foreach my $path (split /,/, $payload) {
			$path =~ s/^\s+|\s+$//g;
			next if ($path eq "");
			# Accept the malformed-but-previously-emitted TECH:a,TECH:b form.
			$path =~ s/^\Q$technology\E:// if ($technology ne "");
			$path = resolve_path($path);
			$path = _join_path($base_dir, $path);
			push @paths, $path;
		}
	}

	die "SupportReads '$spec' does not contain any input paths.\n" if (!@paths);
	return ($technology, \@paths);
}

sub normaliseSupportReads {
	my ($spec, $base_dir) = @_;
	return "" unless (defined($spec) && $spec =~ /\S/);
	my ($technology, $paths) = parseSupportReads($spec, $base_dir);
	return $technology.":".join(",", @{$paths});
}

sub _natural_cmp {
	my ($left, $right) = @_;
	my @left_parts = split /(\d+)/, lc($left);
	my @right_parts = split /(\d+)/, lc($right);
	while (@left_parts && @right_parts) {
		my ($a, $b) = (shift(@left_parts), shift(@right_parts));
		my $cmp = ($a =~ /^\d+$/ && $b =~ /^\d+$/) ? ($a <=> $b) : ($a cmp $b);
		return $cmp if $cmp;
	}
	return @left_parts <=> @right_parts || $left cmp $right;
}

sub _mate_key {
	my ($file) = @_;
	my $key = lc($file);
	my $changed = ($key =~ s/([._-])r[12]([._-])/$1r#$2/i);
	$changed ||= ($key =~ s/([._-])[12](?=[._-](?:f(?:ast)?q|sequence))/$1#/i);
	$changed ||= ($key =~ s/([._-])[12](?=\.f(?:ast)?q(?:\.(?:gz|bz2))?$)/$1#/i);
	return ($key, $changed ? 1 : 0);
}

sub _compile_input_regex {
	my ($label, $pattern) = @_;
	return undef unless (defined($pattern) && $pattern ne "");
	my $regex = eval { qr/$pattern/ };
	die "Invalid $label input-file regular expression '$pattern': $@" if (!$regex);
	return $regex;
}

# Read a directory once, classify each file once, then sort paired reads by a
# shared mate key. This avoids independent greps silently pairing different
# lanes when lexical ordering differs.
sub discoverReadFiles {
	my ($dir, $prefix, $patterns) = @_;
	die "discoverReadFiles expects a pattern hash\n" unless (ref($patterns) eq "HASH");
	my @dirStat = Time::HiRes::stat($dir);
	die "Input directory does not exist: $dir\n"
		unless (@dirStat && S_ISDIR($dirStat[2]));
	$prefix = "" unless defined($prefix);
	my $prefix_regex = qr/^\Q$prefix\E/;
	my %regex = (
		read1  => _compile_input_regex("read-1", $patterns->{read1}),
		read2  => _compile_input_regex("read-2", $patterns->{read2}),
		single => _compile_input_regex("single-read", $patterns->{single}),
		bam    => _compile_input_regex("BAM", $patterns->{bam}),
	);

	my $dirSignature = join(':', @dirStat[0, 1, 7, 9, 10]);
	my $cached = $readDirectoryCache{$dir};
	my @entries;
	if ($cached && $cached->{signature} eq $dirSignature) {
		@entries = @{$cached->{entries}};
	} else {
		opendir(my $dh, $dir) or die "Could not open input directory $dir: $!\n";
		@entries = readdir($dh);
		closedir($dh);
		$readDirectoryCache{$dir} = {
			signature => $dirSignature,
			entries => [@entries],
		};
	}

	# Apply the cheap sample-prefix filter before statting a directory entry.
	# Shared raw-read directories can contain files for thousands of samples.
	my (%fileSizes, @files);
	for my $file (grep { /$prefix_regex/ } @entries) {
		my @fileStat = stat("$dir/$file");
		next unless (@fileStat && S_ISREG($fileStat[2]));
		$fileSizes{$file} = 0 + $fileStat[7];
		push @files, $file;
	}

	my %found = (read1 => [], read2 => [], single => [], bam => [], file_sizes => \%fileSizes);
	foreach my $file (@files) {
		push @{$found{read1}}, $file if ($regex{read1} && $file =~ /$regex{read1}/);
		push @{$found{read2}}, $file if ($regex{read2} && $file =~ /$regex{read2}/);
		push @{$found{single}}, $file if ($regex{single} && $file =~ /$regex{single}/);
		push @{$found{bam}}, $file if ($regex{bam} && $file =~ /$regex{bam}/);
	}

	if ($regex{single}) {
		my %exclude;
		if ($patterns->{prefer_single}) {
			@exclude{@{$found{single}}} = (1) x @{$found{single}};
			@{$found{read1}} = grep { !$exclude{$_} } @{$found{read1}};
			@{$found{read2}} = grep { !$exclude{$_} } @{$found{read2}};
		} else {
			@exclude{@{$found{read1}}, @{$found{read2}}} = (1) x (@{$found{read1}} + @{$found{read2}});
			@{$found{single}} = grep { !$exclude{$_} } @{$found{single}};
		}
	}
	my %read1_names = map { $_ => 1 } @{$found{read1}};
	my @mate_overlap = grep { $read1_names{$_} } @{$found{read2}};
	die "Files match both read-1 and read-2 patterns in $dir: ".join(", ", @mate_overlap)."\n"
		if (@mate_overlap);

	foreach my $type (qw(single bam)) {
		@{$found{$type}} = sort { _natural_cmp($a, $b) } @{$found{$type}};
	}
	if (@{$found{read1}} != @{$found{read2}}) {
		die "Unequal paired-read counts in $dir: ".scalar(@{$found{read1}})." read-1 and ".scalar(@{$found{read2}})." read-2 files.\n";
	}
	my @r1 = map { [$_, _mate_key($_)] } @{$found{read1}};
	my @r2 = map { [$_, _mate_key($_)] } @{$found{read2}};
	@r1 = sort { _natural_cmp($a->[1], $b->[1]) || _natural_cmp($a->[0], $b->[0]) } @r1;
	@r2 = sort { _natural_cmp($a->[1], $b->[1]) || _natural_cmp($a->[0], $b->[0]) } @r2;
	for (my $i = 0; $i < @r1; $i++) {
		if ($r1[$i][2] && $r2[$i][2] && $r1[$i][1] ne $r2[$i][1]) {
			die "Paired-read names do not match in $dir: '$r1[$i][0]' and '$r2[$i][0]'.\n";
		}
	}
	$found{read1} = [map { $_->[0] } @r1];
	$found{read2} = [map { $_->[0] } @r2];
	return \%found;
}

#check if file or file.gz exists
sub resolveExistingFile {
	my ($file) = @_;
	return unless (defined($file) && $file ne '');
	my @candidates = ($file);
	if ($file =~ m/\.gz$/) {
		(my $plain = $file) =~ s/\.gz$//;
		push @candidates, $plain;
	} else {
		push @candidates, "$file.gz";
	}
	my %seen;
	for my $candidate (grep { !$seen{$_}++ } @candidates) {
		my @fileStat = stat($candidate);
		return wantarray ? ($candidate, \@fileStat) : $candidate if @fileStat;
	}
	return;
}

sub fileGZe{
	my ($fil) = @_;
	return 0 unless (defined($fil) && $fil ne '');
	my @candidates = ($fil);
	if ($fil =~ m/\.gz$/) {
		push @candidates, substr($fil, 0, -3);
	} else {
		push @candidates, "$fil.gz";
	}
	my %seen;
	for my $candidate (grep { !$seen{$_}++ } @candidates) {
		my @fileStat = stat($candidate);
		return 1 if (@fileStat && $fileStat[7] > 0);
	}
	return 0;
}


#report file size, check if file or file.gz exists
sub fileGZs{
	my ($fil) = @_;
	return 0 unless (defined($fil) && $fil ne '');
	my @exactStat = stat($fil);
	return $exactStat[7] if @exactStat;
	if ($fil =~ m/\.gz$/) {
		my @plainStat = stat(substr($fil, 0, -3));
		return $plainStat[7] if @plainStat;
	} else {
		my @gzipStat = stat("$fil.gz");
		return $gzipStat[7] * 5 if @gzipStat;
	}
	return 0;
}

# Keep the workflow planner and separateContigs worker on one completion
# contract. The stone alone is insufficient because older/interrupted jobs may
# have published only a subset of the coverage derivatives.
sub contig_stats_coverage_complete{
	my ($dir, $prefix) = @_;
	return 0 unless defined($dir) && defined($prefix) && length($prefix);
	$dir =~ s{/+$}{};
	return 0 unless -e "$dir/$prefix.stone";
	foreach my $suffix (qw(percontig median.percontig pergene count_pergene)){
		return 0 unless fileGZe("$dir/$prefix.$suffix");
	}
	return 1;
}

sub prefixFAhd{
	my ($hr,$nm)= @_;
	my %FNA=%{$hr};
	my %rFNA;
	foreach my $k (keys %FNA){
		my $n;#=$k;
		if ($k =~ m/^>(.*)/){
			$n = $nm.".".$1;
		} else {
			$n = $nm.".".$k;
		}
		$rFNA{$n} = $FNA{$k};
	}
	return(\%rFNA);
}

		
sub gzipwrite{
	my ($outF,$descr) = @_;
	$outF .= ".gz" if ( $outF !~ m/\.gz$/);
	my $O = IO::Compress::Gzip->new($outF)
		or die "error opening gzip output $outF: $IO::Compress::Gzip::GzipError\n";
	#open my $O, ':>gzip', $outF or die "error starting gzip pipe $outF\n$!\n\n";
	#my $pigzBin = getProgPaths("piz");
	#open (my $O, "| $pigzBin -c > $outF") or die "error starting gzip pipe $outF\n$!";
	return $O;
}
#cached path to a working pigz binary; resolved at most once per process.
my $PIGZ_BIN; my $PIGZ_TRIED = 0;
sub _pigzBinCached {
	unless ($PIGZ_TRIED) {
		$PIGZ_TRIED = 1;
		my $p = eval { getProgPaths("pigz") };
		$p = "" if $@ || !defined($p);
		$PIGZ_BIN = (length($p) && -x $p) ? $p : "";
	}
	return $PIGZ_BIN;
}

sub gzipopen{
	my ($inF,$descr) = @_;
	my $dodie = 1;
	if (@_ > 2){$dodie = $_[2];}
	my $verbose=1;
	if (@_ > 3){$verbose = $_[3];}
	my $alreadyResolved = @_ > 4 ? $_[4] : 0;
	my $knownStat = @_ > 5 ? $_[5] : undef;
	my ($resolved, $fileStat);
	if ($alreadyResolved) {
		my @resolvedStat = ref($knownStat) eq 'ARRAY' ? @{$knownStat} : stat($inF);
		($resolved, $fileStat) = ($inF, \@resolvedStat) if @resolvedStat;
	} else {
		($resolved, $fileStat) = resolveExistingFile($inF);
	}
	$inF = $resolved if defined($resolved);
	
	my $ISTR; my $OK = 1;
	my $msg = "Can't open $descr file $inF\n";
	#print "$dodie  $verbose  $inF\n";
	#if (!-e $inF){{if ($dodie){die $msg;} else { $OK=0;print $msg if ($verbose);}}}
	#my $pigzBin = getProgPaths("pigz");

	if (!defined($resolved)) {
		$OK=0;
		if ($dodie) { die $msg; }
		print $msg if ($verbose);
	} elsif($inF =~ m/\.gz$/ ){
		$msg = "Can't open a pipe to $descr file $inF\n";
		my $pigz = _pigzBinCached();
		my $usedPigz = 0;
		if (length($pigz)) {
			#NOTE: unlike IO::Uncompress::Gunzip->new, a successfully-opened pipe doesn't
			#guarantee the gzip stream itself is valid -- corruption would only surface
			#once the caller reads/closes the handle (pigz's own error goes to its stderr,
			#not to $ISTR). We already know $inF exists at this point (resolveExistingFile
			#succeeded above), so the main remaining risk is a truncated/corrupt archive
			#silently reading as if it ended early rather than dying loudly. This is a
			#deliberate trade-off for a large, consistent speedup on the multi-GB gzip
			#files this pipeline reads; ask if you'd like close()-time exit-status checks
			#added on top for stricter corruption detection.
			if (open($ISTR, "-|", $pigz, "-dc", "--", $inF)) {
				$usedPigz = 1;
			}
		}
		if (!$usedPigz) {
			$ISTR = IO::Uncompress::Gunzip->new($inF);
			if (!$ISTR) {if ($dodie){die "$msg$IO::Uncompress::Gunzip::GunzipError\n";} else {$OK=0; print $msg if ($verbose);}}
		}
	} else{
		if (!open($ISTR, "<", "$inF") ) {if ($dodie){die $msg;} else {$OK=0; print $msg if ($verbose);}}
	}
	#print "$OK $inF\n";
	return ($ISTR,$OK,$resolved,$fileStat);
}




sub first_index (&@) {
    my $f = shift;
    for my $i (0 .. $#_) {
	local *_ = \$_[$i];	
	return $i if $f->();
    }
    return -1;
}

#makes sure really all previous split fasta files are removed
sub clenSplitFastas($ $){
	my ($inF , $path) = @_;
	$inF =~ m/\/([^\/]+)$/;
	my $inF2 = $1;
	system "rm -f $path/$inF2.*";
}
sub splitFastas($ $ $){
	my ($inF,$num , $path) = @_;
	system "mkdir -p $path" unless (-d $path);
	my $fCnt = 0; my $curCnt=0;
	my $preNum = $num; 
	#check if a file size is given instead of splits..
	if ($preNum =~ m/[GM]$/){
		my $fileSM =  (-s $inF) / (1024 * 1024);
		$num = $preNum; $num =~ s/[MG]//; 
		$num *= 1024 if ($preNum =~ m/G$/);
		$num = int(1+$fileSM / $num);
	}
	$inF =~ m/\/([^\/]+)$/;
	my $inF2 = $1;
	my @nFiles = ("$path/$inF2.$fCnt.$num");
	#print "$nFiles[-1]\n";
	if ($num < 2){
		print "No split required!\n";
		unlink $nFiles[-1] if (-e $nFiles[-1] || -l $nFiles[-1]);
		symlink $inF, $nFiles[-1] or die "Can't link $nFiles[-1] to $inF: $!\n";
		return \@nFiles;
	}
	if (-e $nFiles[-1] && -e "$path/$inF2.".($num-1).".$num" && !-e "$path/$inF2.$num.$num"){
		print "seems to exist already\n";
		for (my $i=1;$i<$num;$i++){
			push(@nFiles,"$path/$inF2.$i.$num");
		}
		return \@nFiles;
	}
	unlink $nFiles[-1] or die "Can't remove old split $nFiles[-1]: $!\n" if (-e $nFiles[-1]);
	open my $count_fh, '<', $inF or die "Can't open FASTA $inF: $!\n";
	my $protN = 0;
	while (my $line = <$count_fh>) { $protN++ if ($line =~ /^>/); }
	close $count_fh;
	my $pPerFile = int($protN/$num)+10;
	open I,"<$inF" or die "Can't open FASTA $inF: $!\n";
	open my $out,">".$nFiles[-1] or die "Can't write FASTA split $nFiles[-1]: $!\n";
	while (my $l = <I>){
		if ($l =~ m/^>/){
			$curCnt++;
			if ($curCnt > $pPerFile){
				$fCnt++; close $out; 
				push(@nFiles,"$path/$inF2.$fCnt.$num");
				open $out,">$nFiles[-1]" or die "Can't write FASTA split $nFiles[-1]: $!\n";
				$curCnt=0;
			}
		}
		print $out $l;
	}
	close I; close $out;
	return \@nFiles;
	
}


sub readFasta{
	my $fils = $_[0];
	my $cutHd=0;
	my $sepChr= "\\s";
	my $subs; my $doSubs=0;
	$cutHd = $_[1] if (@_ > 1);
	$sepChr = $_[2] if (@_ > 2);
	if (@_ > 3){
		$subs = $_[3];
		die "readFasta subset must be a hash reference\n" unless ref($subs) eq 'HASH';
		$doSubs = 1;
	}
	my $Hseq = {};

	my @files = glob $fils;
	#if (-z $fil){ return \%Hseq;}
	foreach my $fil (@files){
		#next unless (-e $fil);
		#my $FAS; my $status=0; 
		#my $doAdd = 1;
		#open($FAS,"<","$fil") || die ("Couldn't open FASTA file $fil\n");
		my ($FAS ,$status) = gzipopen($fil,"fasta file to readFasta",0);
		if (@files == 1 && $status == 0){die "Can't open fasta file $fil\n";}
		next if ($status==0);
		my $first_line = <$FAS>;
		if (!defined $first_line){#could be empty file
			print "Empty:: $fil $status\n";
			close $FAS;
			next;
		}
		die "Malformed FASTA file $fil: first record has no header\n"
			unless $first_line =~ /^>/;

		my $sepRe = qr/$sepChr/; #compile once, reused for every header in this file
		my $prepare_header = sub {
			my ($line) = @_;
			chomp $line;
			my $full_header = substr($line, 1);
			my $short_header = $full_header;
			$short_header =~ s/$sepRe.*//;
			my $stored_header = $cutHd ? $short_header : $full_header;
			return ($stored_header, $short_header);
		};
		my ($trHe, $srcHe) = $prepare_header->($first_line);
		my $temp = "";
		my $store_record = sub {
			return if $doSubs && !exists($subs->{$trHe}) && !exists($subs->{$srcHe});
			$Hseq->{$trHe} = $temp;
		};

		while(my $line = <$FAS>){
			chomp($line);
			if ($line =~ m/^>/){
				$store_record->();
				($trHe, $srcHe) = $prepare_header->($line);
				$temp = "";
				next;
			}
			$temp .= ($line);
		}
		$store_record->();
		close ($FAS) or die "Cannot close FASTA file $fil: $!\n";
	}
	return $Hseq;
}


sub median
{
    my @vals = sort {$a <=> $b} @_;
	return 0 if (@vals == 0);
    my $len = (scalar @vals)-1;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}

sub quantile #format: quantile(0.25,@values);
{
	my $cut = shift @_;
	if ($cut<0 || $cut > 1){die "GenoMetaAss::quantile: requested quantile <0 or >1: $cut\nAborting..";}
	if (@_ == 0){return;}
	if (@_ == 1){return ($_[0]);}
	#print "X@_\n\n" if (!defined($_[0]));
    my @vals = sort {$a <=> $b} @_;
    my $len = (scalar @vals)-1;
    return $vals[int(($len*$cut) + 0.5)];
}
sub mean{
    my @vals = @_;
    my $len = scalar @vals;
	return 0 if (@vals == 0);
	my $sum =0; foreach my $a (@vals){$sum+=$a;}
	return ($sum/$len);
}

sub convertNT2AA($){
	my ($text) = @_;
	#print $text."\n";
	 # translate a DNA 3-character codon to an amino acid. We take three letter groups from
	# the incoming string of C A T and G and translate them via a Hash.  It's listed vertically
	# so we can label each of the amino acids as we set it up.

	#$text = "aaatgaccgatcagctacgatcagctataaaaaccccggagctacgatcatcg";

	$text =~ s/[N]*$//i;
	if (length($text) % 3 != 0){
		print "Input not correct length!\n";
		my $ltr = length($text) % 3;
		my $lt = length($text);
		$text = substr $text,0,($lt -$ltr);
	}



	my %convertor = (
		'TCA' => 'S',    # Serine
		'TCC' => 'S',    # Serine
		'TCG' => 'S',    # Serine
		'TCT' => 'S',    # Serine
		'TTC' => 'F',    # Phenylalanine
		'TTT' => 'F',    # Phenylalanine
		'TTA' => 'L',    # Leucine
		'TTG' => 'L',    # Leucine
		'TAC' => 'Y',    # Tyrosine
		'TAT' => 'Y',    # Tyrosine
		'TAA' => '*',    # Stop
		'TAG' => '*',    # Stop
		'TGC' => 'C',    # Cysteine
		'TGT' => 'C',    # Cysteine
		'TGA' => '*',    # Stop
		'TGG' => 'W',    # Tryptophan
		'CTA' => 'L',    # Leucine
		'CTC' => 'L',    # Leucine
		'CTG' => 'L',    # Leucine
		'CTT' => 'L',    # Leucine
		'CCA' => 'P',    # Proline
		'CCC' => 'P',    # Proline
		'CCG' => 'P',    # Proline
		'CCT' => 'P',    # Proline
		'CAC' => 'H',    # Histidine
		'CAT' => 'H',    # Histidine
		'CAA' => 'Q',    # Glutamine
		'CAG' => 'Q',    # Glutamine
		'CGA' => 'R',    # Arginine
		'CGC' => 'R',    # Arginine
		'CGG' => 'R',    # Arginine
		'CGT' => 'R',    # Arginine
		'ATA' => 'I',    # Isoleucine
		'ATC' => 'I',    # Isoleucine
		'ATT' => 'I',    # Isoleucine
		'ATG' => 'M',    # Methionine
		'ACA' => 'T',    # Threonine
		'ACC' => 'T',    # Threonine
		'ACG' => 'T',    # Threonine
		'ACT' => 'T',    # Threonine
		'AAC' => 'N',    # Asparagine
		'AAT' => 'N',    # Asparagine
		'AAA' => 'K',    # Lysine
		'AAG' => 'K',    # Lysine
		'AGC' => 'S',    # Serine
		'AGT' => 'S',    # Serine
		'AGA' => 'R',    # Arginine
		'AGG' => 'R',    # Arginine
		'GTA' => 'V',    # Valine
		'GTC' => 'V',    # Valine
		'GTG' => 'V',    # Valine
		'GTT' => 'V',    # Valine
		'GCA' => 'A',    # Alanine
		'GCC' => 'A',    # Alanine
		'GCG' => 'A',    # Alanine
		'GCT' => 'A',    # Alanine
		'GAC' => 'D',    # Aspartic Acid
		'GAT' => 'D',    # Aspartic Acid
		'GAA' => 'E',    # Glutamic Acid
		'GAG' => 'E',    # Glutamic Acid
		'GGA' => 'G',    # Glycine
		'GGC' => 'G',    # Glycine
		'GGG' => 'G',    # Glycine
		'GGT' => 'G',    # Glycine
		);

	# We don't actually know where the groups of 3 will start in our sample piece of DNA, so we've got
	# three ways of doing the coding ... here's a loop to work out each of the possibilities in turn,
	# leaving the odd extra letters on the beginning or end.
	#for ($s=0; $s<3; $s++) {
	#last; #bs, just need on code
	#        $scrap = substr($text,0,$s);
	#        $main = substr($text,$s);
	#        $main =~ s/(...)/"$convertor{uc $1}" || "?"/eg;
	#        print "$scrap$main\n";
	#        }
	$text =~ s/(...)/"$convertor{uc $1}" || "?"/eg;
	$text =~ s/[^ACDEFGHIKLMNPQRSTVWY*]/X/g;
	#die $text;
	return $text;
}



sub lcp {#
	(join("\0", @_) =~ /^ ([^\0]*) [^\0]* (?:\0 \1 [^\0]*)* $/sx)[0];
}
sub prefix_find($){
	my ($ar) = @_;
	my @rds = @{$ar}; my @newRds = @rds;
	my $first = $rds[0]; 
	#die "$first FFF\n";
	for (my $i=0;$i<@rds; $i++){

		my $second = $rds[$i];
		my @matches;
		next unless (defined $first);
		for (my $start = 0; $start < length($first); $start++) {
			for (my $len = $start+1; $len< length($first); $len++) {
				my $substr = substr($first, $start, $len);
				push @matches, $second =~ m[($substr)]g;
			}
		}
		#print "@matches\n";
		my ($len, $longest) = 0;
		length > $len and ($longest, $len) = ($_, length) for @matches;
		$first = $longest;
		#print "$longest $i\n";
		#"$first\0$second" =~ m/^(.*)\0\1/s;
		#$first = $1;
	}
	
	#print "$first common prefix\n";
	#for (my $i=0; $i<@newRds; $i++){
	#	$newRds[$i] =~ s/^$first//;
	#	$newRds[$i] .= "/" if ($newRds[$i] !~ m/\/$/ && length($newRds[$i] ) != 0 );
	#}
	return ($first);
}

 
sub is_integer {
   defined $_[0] && $_[0] =~ /^[+-]?\d+$/;
}




sub readFastHD{#only reads headers
	my $inF = $_[0];
	my $complH=0;
	if (@_ > 1){$complH = $_[1];}
	#die "$complH\n";
	my @ret;
	open I,"<$inF" or die "can t open $inF\n";
	while (my $l = <I>){
		chomp $l;
		if ($l =~ m/^>(\S+)/){
			if (!$complH){
				push @ret,$1;
			} else {
				$l =~ m/^>(.*)/;
				push @ret,$1;
			}
		}
	}
	close I;
	return \@ret;
}

sub renameFastHD($ $ $){ #set a new name for headers in fasta files
	my ($inF,$hr,$extr) = @_;
	my %COGid = %{$hr};
	my $oS = "";
	my %CidCnt;
	open I,"<$inF" or die "can t open renameFastHD  $inF\n";
	while (my $l = <I>){
		if ($l =~ m/^>/){
			chomp $l;
			unless (exists $COGid{$l}){die "Can't find $l in COGs\n";}
			if (exists $CidCnt{$COGid{$l}}){$CidCnt{$COGid{$l}}++;
				$oS .= ">".$extr."_".$COGid{$l}."..".$CidCnt{$COGid{$l}}."\n";
			} else {$CidCnt{$COGid{$l}} = 0;
				$oS .= ">".$extr."_".$COGid{$l}."\n";
			}
			

		} else {
			$oS .= $l;
		}
	}
	close I;
	open O,">$inF" or die "Can't open rename Fasta HD out file $inF\n";
	print O $oS;
	close O;
}
sub renameFastqCnts($ $){ #set a new name for headers in fastq files, using a simple scheme of prefix and just counts afterwards
	my ($inF,$prefix) = @_;
	open I,"<$inF" or die "can t open $inF\n";
	my $cnt  = 0; my $line_in_record = 0;
	open O,">$inF.tmp" or die "can't open $inF.tmp\n";
	while (my $l = <I>){
		if ($line_in_record == 0){
			die "Malformed FASTQ record in $inF: expected header, got $l"
				unless ($l =~ /^@/);
			print O "@".$prefix."_".$cnt."\n";
			$cnt++;
		} else {
			print O $l;
		}
		$line_in_record = ($line_in_record + 1) % 4;
	}
	die "Malformed FASTQ file $inF: incomplete final record\n" if ($line_in_record != 0);
	close I; close O;
	rename "$inF.tmp", $inF or die "Can't replace $inF: $!\n";
}

sub renameFastaCnts($ $ $){ #set a new name for headers in fasta files, using a simple scheme of prefix and just counts afterwards
	my ($inF,$prefix,$logF) = @_;
	open I,"<$inF" or die "can t open $inF\n";
	my $cnt  = 0;
	open O,">$inF.tmp" or die "can't open $inF.tmp\n";
	open L,">$logF" or die "can't open $logF\n";
	while (my $l = <I>){
		if ($l =~ m/^>/){
			print O ">".$prefix."_".$cnt."\n";
			print L ">".$prefix."_".$cnt."\t$l";
			$cnt ++;
		} else {
			print O $l;
		}
	}
	close I; close O; close L;
	rename "$inF.tmp", $inF or die "Can't replace $inF: $!\n";
}

sub readClstrRevSmplCtgGenSubset{
	my $inF = $_[0];
	my $subsHR = {};
	$subsHR = $_[1] if (@_ > 1);
	#my %subs = %{$subsHR};
	my $retR={}; my %retF;
	
	#my @key = keys %{$subsHR}; print "$key[0] $key[1]\n";
	print "Reading clstr $inF  .. \n";
	my ($I,$OK) = gzipopen($inF,"ClstrFile");
	#open I,"<$inF" or die "Can't find rev clustering file $inF\n"; 
	my $curCl=""; my $cnts=0;
	while (my $lin = <$I>){
		chomp $lin; 
		next if ($lin =~ m/^#/ || length($lin) < 5);
		my @arr = split /\t/,$lin,-1;
		#my $pos = index($_, "\t");
		$curCl = ($arr[0]);#substr($_,0,$pos);
		
		#extract reverse: contig_gene to gene cat 
		my @tmpArr = split /,/,$arr[1];
		foreach my $gene (@tmpArr) {
			my $ctg = substr($gene,1); $ctg =~ s/_(\d+)$//; my $geneN=$1;
			#print "$ctg = $gene = $geneN\n";
			 if (exists($$subsHR{$ctg})){
				 my @spl2 = split(/__/,$ctg);
				#${$retR}{$ctg}{$geneN} = $curCl;
				${$retR}{$spl2[0]}{$spl2[1]}{$geneN} = $curCl;
				$cnts++;
			 }
		}
	}
	close $I;
	print "Found $cnts genes.\n";
	#die;
	return ($retR);
}

sub readClstrRevContigSubset{ #gets the exact assembled genes clustered in GC genes, subset by contigs (used for clusterMAGs.pl)
	my $inF = $_[0];
	my $subsHR = {};
	$subsHR = $_[1] if (@_ > 1);
	#my %subs = %{$subsHR};
	my $retR={}; my %retF;
	
	#my @key = keys %{$subsHR}; print "$key[0] $key[1]\n";
	print "Reading clstr $inF  .. ";
	my ($I,$OK) = gzipopen($inF,"ClstrFile");
	#open I,"<$inF" or die "Can't find rev clustering file $inF\n"; 
	my $curCl=""; my $cnts=0;
	while (my $lin = <$I>){
		chomp $lin; 
		next if ($lin =~ m/^#/ || length($lin) < 5);
		my @arr = split /\t/,$lin,-1;
		#my $pos = index($_, "\t");
		$curCl = ($arr[0]);#substr($_,0,$pos);
		
		#extract reverse: contig_gene to gene cat 
		my @tmpArr = split /,/,$arr[1];
		foreach my $gene (@tmpArr) {
			my $ctg = substr($gene,1); $ctg =~ s/_\d+$//; 
			#print "$ctg = $gene\n";
			 if (exists($$subsHR{$ctg})){
				${$retR}{$gene} = $curCl;
				$cnts++;
			 }
		}
	}
	close $I;
	print "Found $cnts genes.\n";
	#die;
	return ($retR);
}


sub readClstrRev{ #gets the exact assembled genes clustered in GC genes #version for my shortened index file
	my $inF = $_[0];
	my $createR = 1;
	$createR = $_[1] if (@_ > 1);
	my $subsHR = {};
	$subsHR = $_[2] if (@_ > 2);
	#optional 4th arg: hashref of sample IDs (the part of a member name before "__").
	#when given, member lists stored in %retF are pre-filtered to only these samples,
	#so callers that only need a subset of samples (e.g. one split-worker out of many)
	#never have to hold the other samples' members in memory in the first place.
	my $memberSubsHR = undef;
	#An explicitly supplied empty hash is a valid "keep none" partition.  Do
	#not turn it into undef, which means "do not filter" and would make an
	#otherwise idle split worker load the complete catalogue.
	$memberSubsHR = $_[3] if (@_ > 3 && ref($_[3]) eq 'HASH');
	#my %subs = %{$subsHR};
	my $doSubset=0;
	if (scalar(keys(%{$subsHR})) > 0 ){
		$doSubset = 1 ;
		if ($createR == 2){
			$doSubset = 2 ;
		}
	}
	my $retR={}; my %retF;
	print "Reading clstr $inF  .. \n";
	print "Restricting cluster members to a subset of " . scalar(keys %{$memberSubsHR}) . " samples\n"
		if $memberSubsHR;
	my ($I,$OK) = gzipopen($inF,"ClstrFile");
	#open I,"<$inF" or die "Can't find rev clustering file $inF\n"; 
	my $curCl=0;
	while (my $lin = <$I>){
		chomp $lin; 
		next if ($lin =~ m/^#/ || length($lin) < 5);
		my @arr = split /\t/,$lin,-1;
		#my $pos = index($_, "\t");
		$curCl = ($arr[0]);#substr($_,0,$pos);
		
		if ($doSubset == 1){
			next unless (exists($$subsHR{$curCl}));
		}
		
		#my $rem = $arr[1]; #substr($_,$pos+1);
		#print $curCl."XX$rem\n";
		#foreach (split /,/,$rem){$retR{$_} = $curCl;}
		my @tmpArr; my $tmpArrSplit = 0;
		if ($createR > 0){
			@tmpArr = split /,/,$arr[1]; $tmpArrSplit = 1;
			if ($doSubset == 1){
				my $gogo=0;
				foreach (@tmpArr) {
					 if (exists($$subsHR{$_})){
						 $gogo=1;
						 last;
					 }
				}
				next unless ($gogo);
			}
			@{$retR}{@tmpArr} = ($curCl) x @tmpArr;
		}
		if ($createR != 2){
			if ($memberSubsHR){
				@tmpArr = split /,/,$arr[1] unless $tmpArrSplit;
				my @keep;
				foreach my $mem (@tmpArr) {
					my $m2 = $mem; $m2 =~ s/^>//;
					my ($smp) = split /__/, $m2, 2;
					push @keep, $mem if defined($smp) && exists($memberSubsHR->{$smp});
				}
				$retF{$curCl} = join(",", @keep) if @keep;
			} else {
				$retF{$curCl} = $arr[1];
			}
		}
		@arr = ();

	}
	close $I;
	print "done reading ClStr\n";
	return ($retR,\%retF);
}


sub readClstrRevGenes{ #gets the exact assembled genes clustered in GC genes #version for my shortened index file
	my $inF = $_[0];
	my $subsHR = {};
	$subsHR = $_[1] if (@_ > 1);
	#my %subs = %{$subsHR};
	my $doSubset=0;
	$doSubset = 1 if (scalar(keys(%{$subsHR})) > 0 );
	my %retR; 
	print "Reading clstr2 $inF  .. \n";
	my ($I,$OK) = gzipopen($inF,"ClstrFile");
	#open I,"<$inF" or die "Can't find rev clustering file $inF\n"; 
	my $curCl=0;
	while (my $lin = <$I>){
		chomp $lin; 
		next if ($lin =~ m/^#/ || length($lin) < 5);
		my @arr = split /\t/,$lin,-1;
		#my $pos = index($_, "\t");
		$curCl = int($arr[0]);#substr($_,0,$pos);
		
		if ($doSubset){
			next unless (exists($subsHR->{$curCl}));
		}
		my @tmpArr = split /,/,$arr[1];
		foreach (@tmpArr) {
			next unless (exists($subsHR->{$_}));
			$retR{$_} = $curCl;
		}
#			@retR{@tmpArr} = ($curCl) x @tmpArr;
	}
	close $I;
	print "done reading ClStr gene 2\n";
	return (\%retR);
}



sub unzipFileARezip($){
	my ($inFar) = @_;
	my @inFs = @{$inFar}; my $totCnt=0;
	my $bef = "gunzip "; my $aft="gzip ";
	#print "$inF.gz\n";
	foreach my $inF (@inFs){
		if (-e "$inF.gz"){
			$bef .= " $inF.gz"; $aft .= " $inF"; $totCnt++;
		}
	}
	if ($totCnt==0) {
		return ("","");
	} else {
		return $bef."\n",$aft."\n";
	}
}
sub reverse_complement_IUPAC ($) {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
        return $revcomp;
}


sub readGFF($){
	my ($inF) =@_;
	my %ret;my $sbcnt=1;
	my $lcnt=0;
	#if (!-e $inF && -e "$inF.gz"){$inF .= ".gz";}
	my $entries=0;
	my ($GFF ,$status) = gzipopen($inF,"gff file to readGFF",1);

	#open I,"<$inF";
	while (<$GFF>){
		if (m/^#/){$sbcnt=1;next;}
		chomp;
		$lcnt++;
		my @spl = split(/\t/, $_, -1);
		if (@spl < 9) {
			warn "readGFF::too short: line $lcnt, $inF\n";
			next;
		}
		my ($feature_id) = $spl[8] =~ /(?:^|;)ID=\d+_(\d+)(?:;|$)/;
		unless (defined $feature_id) {
			warn "readGFF::missing supported ID attribute: line $lcnt, $inF\n";
			next;
		}
		#my $k = ">".$spl[0]."_$1";
		my $k = $spl[0]."_$feature_id";
		#print $k;
		$ret{$k}=$_;
		$entries++;
	}
	
	close $GFF;
	#die "Found $entries gff entries\n";
	return \%ret;
}



sub getAssemblPath{
	my $cD = $_[0];
	my $newpath = "";
	$newpath = $_[1] if (@_ > 1);
	$newpath =~ s/\/\//\//g;$newpath =~ s/\/\//\//g;
	my $dieOnFail=1;
	$dieOnFail = $_[2] if (@_ > 2);
	my $firstF = "$cD/assemblies/metag/assembly.txt";
	if (!-e $firstF){
		print "GenoMetaAss::getAssemblPath::Assembly path missing: $cD\n";
		return ("") if (!$dieOnFail);
		
	}
	open my $I,"<$firstF" or die "GenoMetaAss::getAssemblPath::Assembly path missing: $cD\n";
	my $metaGD = <$I>;#= `cat $cD/assemblies/metag/assembly.txt`; 
	chomp $metaGD;
	my $metaGD2 = $metaGD;$metaGD2 =~ s/\/\//\//g;$metaGD2 =~ s/\/\//\//g;
	
	close $I;
	if ($newpath ne ""){
		if ($newpath ne $metaGD2 || !-d $metaGD  ){#needs replacement?
			print "replacing $metaGD with $newpath\n" ;
			open my $I,">$firstF"; print $I $newpath; close $I;
			$metaGD = $newpath;
		}
	}
	
	return($metaGD);
}



sub getAssemblGFF{
	my $cD = $_[0];
	my $dieOnFail = 1; $dieOnFail = $_[1] if (@_ > 1);
	my $assmD = getAssemblPath($cD,"",$dieOnFail);
	return "" if ($assmD eq "");

	my $gffF = "$assmD//genePred/genes.gff";
	$gffF .= ".gz" if (-e $gffF . ".gz");
	if (!-s $gffF){print STDERR "getAssemblGFF::Could not find gff in dir $cD!:\n$gffF\n";die if ($dieOnFail);}
	return $gffF;
}

sub getAssemblContigs{
	my $cD = $_[0];
	my $dieOnFail = 1; $dieOnFail = $_[1] if (@_ > 1);
	my $assmD = getAssemblPath($cD,"",$dieOnFail);
	return "" if ($assmD eq "");
	my $ctgF = "$assmD/scaffolds.fasta.filt";
	$ctgF .= ".gz" if (!-e $ctgF && -e "$ctgF.gz");
	if (!-s $ctgF){print STDERR "getAssemblContigs::Could not find assembly in dir $cD!:\n$ctgF\n";die if ($dieOnFail);}
	return $ctgF;
}


sub iniCleanSeqSetHR{
	my ($seqSetHR) = @_;
	my $libraries = cloneReadLibraries(ensureSeqSetLibraries($seqSetHR));
	foreach my $library (@{$libraries}) {$library->{phase} = 'clean';}
	my $HR = {libraries => $libraries, mrgHshHR => {}};
	syncCleanSeqSetLegacy($HR);
	return $HR;
}

sub addFileLocs2AssmGrp{
	my ($AsGrpsHR, $cAssGrp,$SmplName, $cleanSeqSetHR, $seqSetHR) = @_;
	ensureCleanSeqSetLibraries($cleanSeqSetHR, $SmplName);
	ensureSeqSetLibraries($seqSetHR, $SmplName);
	$AsGrpsHR->{$cAssGrp}{InputOrder} ||= [];
	push @{$AsGrpsHR->{$cAssGrp}{InputOrder}}, $SmplName
		unless grep { $_ eq $SmplName } @{$AsGrpsHR->{$cAssGrp}{InputOrder}};
	${${$AsGrpsHR}{$cAssGrp}{CleanSeqs}}{$SmplName} = $cleanSeqSetHR;
	${${$AsGrpsHR}{$cAssGrp}{RawSeqs}}{$SmplName} = $seqSetHR;
	return $AsGrpsHR;
}

sub _ordered_group_samples {
	my ($asG, $grp, $store, $specific) = @_;
	my %by_sample = %{$asG->{$grp}{$store} || {}};
	return ($specific) if (defined($specific) && $specific ne '');
	my @samples = grep { exists($by_sample{$_}) } @{$asG->{$grp}{InputOrder} || []};
	@samples = sort keys %by_sample if (!@samples);
	return @samples;
}

sub getCleanLibrariesAssmGrp {
	my ($asG, $grp, $support, $specific) = @_;
	my @libraries;
	foreach my $sample (_ordered_group_samples($asG, $grp, 'CleanSeqs', $specific)) {
		my $clean = $asG->{$grp}{CleanSeqs}{$sample}
			or die "Clean read libraries for sample $sample are unavailable in group $grp\n";
		push @libraries, @{readLibrariesByScope($clean, $support ? 'support' : 'primary', 1, $sample)};
	}
	return \@libraries;
}

sub getRawLibrariesAssmGrp {
	my ($asG, $grp, $support, $specific) = @_;
	my @libraries;
	foreach my $sample (_ordered_group_samples($asG, $grp, 'RawSeqs', $specific)) {
		my $raw = $asG->{$grp}{RawSeqs}{$sample}
			or die "Raw read libraries for sample $sample are unavailable in group $grp\n";
		push @libraries, @{readLibrariesByScope($raw, $support ? 'support' : 'primary', 0, $sample)};
	}
	return \@libraries;
}


sub getCleanSeqsAssmGrp{
	my ($asG, $grp, $support) = @_;
	my $specSmpl = "";$specSmpl = $_[3] if (@_ >= 4); #specific sample only??
	my @smpls = ();@smpls = ($specSmpl) if ($specSmpl ne "");
	
	my $libraries = getCleanLibrariesAssmGrp($asG, $grp, $support, $specSmpl);
	my ($pa1, $pa2, $pas, $labels, $readTec) = legacyLibraryArrays($libraries, 1);
	#die "getCleanSeqsAssmGrp::\n@pa1 - $grp - $support\n@pa2\n@pas\n";
	return ($pa1, $pa2, $pas, $readTec);
}

sub hasSuppRds{
	my ($asG, $grp,$smpl) = @_;
	return scalar(@{getRawLibrariesAssmGrp($asG, $grp, 1, $smpl)}) ? 1 : 0;
}

sub getRawSeqsAssmGrp{
	my ($asG, $grp, $support) = @_;
	
	my $specSmpl = "";$specSmpl = $_[3] if (@_ >= 4); #specific sample only??
	my @smpls = ();@smpls = ($specSmpl) if ($specSmpl ne "");
	
	my $libraries = getRawLibrariesAssmGrp($asG, $grp, $support, $specSmpl);
	my ($pa1, $pa2, $pas, $libs, $readTec) = legacyLibraryArrays($libraries, 0);
	#die "@pa1\n@pa2\n@pas\n";
	return ($pa1, $pa2, $pas, $libs, $readTec);
}


sub emptyAssGrpsObj($){
	my ($hr) = @_;
	my %AsGrps = %{$hr};
	foreach my $k (keys %AsGrps){
		$AsGrps{$k}{CntAss} = 0;
		$AsGrps{$k}{CntMap} = 0;
		#dependencies
		$AsGrps{$k}{ITSDeps} = "";
		$AsGrps{$k}{MapDeps} = "";
		$AsGrps{$k}{DiamDeps} = "";
		$AsGrps{$k}{SeqClnDeps} = "";
		#print $AsGrps{$k}{CntAimAss} ."\n";
		#filteredSequenceFiles for assembly
		$AsGrps{$k}{PostAssemblCmd} = "";
		$AsGrps{$k}{PostClnCmd} = "";
		$AsGrps{$k}{PostConsCmd} = "";
		$AsGrps{$k}{AssemblSmplDirs} = "";
		$AsGrps{$k}{scndMapping} = "";
		#complex hashes replaces FilterSeq1 ..
		$AsGrps{$k}{CleanSeqs} = {};
		$AsGrps{$k}{RawSeqs} = {};
		$AsGrps{$k}{InputOrder} = [];
		#@{$AsGrps{$k}{FilterSeq1}} = (); @{$AsGrps{$k}{FilterSeq2}} = (); @{$AsGrps{$k}{FilterSeqS}} = ();
		#@{$AsGrps{$k}{RawSeq1}} = (); @{$AsGrps{$k}{RawSeq2}} = (); @{$AsGrps{$k}{Libs}} = ();
		
		#print $k."\n";;
	}
	$AsGrps{global}{DiamCln} = "";
	$AsGrps{global}{DiamDeps} = "";
	#die;
	return \%AsGrps;
}



sub resetAsGrps{
	my ($AsGrps) = @_;
	foreach my $cAssGrp (keys %{$AsGrps}){
		$AsGrps->{$cAssGrp}{CntMap} = 0;
		$AsGrps->{$cAssGrp}{CntAss} = 0;
		$AsGrps->{$cAssGrp}{CntPreAss} = 0;
		$AsGrps->{$cAssGrp}{CntPreAssMiss} = 0;
		$AsGrps->{$cAssGrp}{CntPreAssNoPrim} = 0;
		$AsGrps->{$cAssGrp}{preAsmblDir} = [];
		$AsGrps->{$cAssGrp}{AssemblSmplDirs} = "";
		$AsGrps->{$cAssGrp}{BinDeps} = "";
		$AsGrps->{$cAssGrp}{readDeps} = "";
		$AsGrps->{$cAssGrp}{scndMapping} = "";
		$AsGrps->{$cAssGrp}{SeqUnZDeps} = "";
		# UnzpDeps is copied into readDeps at the start of each pass. Keeping
		# completed job IDs here makes the next loop depend on stale Slurm jobs.
		$AsGrps->{$cAssGrp}{UnzpDeps} = "";
		@{$AsGrps->{$cAssGrp}{FilterSeq1}} = ();
		@{$AsGrps->{$cAssGrp}{FilterSeq2}} = ();
		@{$AsGrps->{$cAssGrp}{FilterSeqS}} = ();
		@{$AsGrps->{$cAssGrp}{ReadTec}} = ();
		$AsGrps->{$cAssGrp}{CleanSeqs} = {};
		$AsGrps->{$cAssGrp}{RawSeqs} = {};
		$AsGrps->{$cAssGrp}{InputOrder} = [];
		$AsGrps->{$cAssGrp}{SeqClnDeps} = "";
		$AsGrps->{$cAssGrp}{prodRun} = "";
		$AsGrps->{$cAssGrp}{AssemblJobName} = "";
		$AsGrps->{$cAssGrp}{PostAssemblCmd} = "";
		$AsGrps->{$cAssGrp}{MapDeps} = "";
		$AsGrps->{$cAssGrp}{BinDeps} = "";
		$AsGrps->{$cAssGrp}{PostClnCmd} = "";
		$AsGrps->{$cAssGrp}{PostConsCmd} = "";
	}
}



sub readMapS{
	my ($inF,$folderStrClassical) = ($_[0],$_[1]);
	my $xtraCols = ""; $xtraCols = $_[2] if (@_ > 2);
	my @spl = split /,/,$inF;
	#my %ret; my %agbp;
	my @outDirs = ();
	#my $hr1 = \%ret;my $hr2 = \%agbp;  my $cnt = -1;
	my $hr1 = {};my $hr2 = {};  my $cnt = -1;
	foreach my $map (@spl){
		($hr1,$hr2) = readMap($map,$cnt,$hr1,$hr2,$folderStrClassical,$xtraCols);
		#%ret = %{$hr1};
		$cnt = $hr1->{opt}{totSmpls};
		push(@outDirs,$hr1->{opt}{outDir});
		#print "\n".$cnt."\n";
	}
	#%ret = %{$hr1};
	#print keys %ret;
	$hr1->{opt}{outDir} = join(",",@outDirs);
	#die;
	return ($hr1,$hr2);
}

#infer Assembly dirs & corrsponding bams with several Samples (compound assemblies)
sub getDirsPerAssmblGrp{
	my %AsGrps; my %map;
	my $chkDirs=0;
	if (@_ == 1){
		my ($mapF) = @_;
		my ($hrm,$asGrpObj) = readMapS($mapF);
		%AsGrps = %{$asGrpObj};		%map = %{$hrm};
		$chkDirs=1;
	} elsif (@_ == 2) {
		my ($hrm,$asGrpObj) = @_;
		%AsGrps = %{$asGrpObj};		%map = %{$hrm};
	} else {
		die "Too many/few args given in GenoMetaAss::getDirsPerAssmblGrp\n";
	}
	my @smpls = @{$map{opt}{smpl_order}};
	my $cnt=0;
	my %DOs;
	foreach my $smpl (@smpls){
		my $cntAim = $AsGrps{ $map{$smpl}{MapGroup} }{CntAimMap};
		my $dir2rd = $map{$smpl}{wrdir};
		my $cAssGrp = $map{$smpl}{AssGroup};
		my $cMapGrp = $map{$smpl}{MapGroup};
		$AsGrps{$cMapGrp}{CntMap} ++;
		#next if ($AsGrps{$cMapGrp}{CntMap}  < $AsGrps{$cMapGrp}{CntAimMap} );

		my $tar = "AssmblGrp_$cAssGrp";
		if ($chkDirs && !-d "$dir2rd"){
			print "Can't read $dir2rd\nSkipping..\n";
			next;
		}
		#(-s "$curOutDir/mapping/Align_ment-smd.bam" || -s "$curOutDir/mapping/Align_ment-smd.cram")
		#my $bam = $inD."mapping/Align_ment-smd.bam";
		#print "$cAssGrp\n";
		push (@{$DOs{$cAssGrp}{wrdir}}, $dir2rd);
		push(@{$DOs{$cAssGrp}{SmplID}},$map{$smpl}{SmplID});
		$cnt++;
		#last if ($cnt > 50);
	}
	return(\%DOs,\%map);
}


#should only be used for dirs!
sub resolve_path($){
	my ($inP) = @_;
	return "" if ($inP eq "");
	my $originalPath = $inP;
	$inP =~ s{
		\$(?:\{([A-Za-z_][A-Za-z0-9_]*)\}|([A-Za-z_][A-Za-z0-9_]*))
	}{
		my $envName = defined($1) ? $1 : $2;
		die "Environment variable \$$envName used in path '$originalPath' is not set\n"
			unless (exists($ENV{$envName}) && defined($ENV{$envName}) && $ENV{$envName} ne '');
		$ENV{$envName};
	}gex;
	#doesn't work: dir might not exist yet!
	#$inP = `realpath $inP` ; chomp $inP; $inP .= "/";
	#if (-d $inP && $inP !~ m/\/$/){$inP .= "/";}
	$inP =~ s/\/\//\//g;$inP =~ s/\/\//\//g;
	return $inP;
}

sub is3rdGenSeqTech{
	my $curReadTec = $_[0];
	if ($curReadTec eq ""){return -1;}
	my $is3rdGen = 0; #important flag for long reads (Oxford Nanopore / PacBio)
	$is3rdGen = 1 if ($curReadTec eq "ONT" || $curReadTec eq "PB");
	return $is3rdGen;
}

sub checkSeqTech{
	#ONT,PB,proto,miSeq,AVITI,GAII
	my $inT = $_[0];
	my $msg = "Mapping file";
	$msg = $_[1] if (@_ > 1);
	if ($inT ne "" && $inT ne "ONT"&& $inT ne "ill"&& $inT ne "hiSeq" && $inT ne "454" && $inT ne "SLR" && $inT ne "PB" && $inT ne "proto" && $inT ne "miSeq" && $inT ne "AVITI" && $inT ne "GAII" && $inT ne "GAII_solexa"){
		die "$msg: Can't recognize SeqTech: \"$inT\"\nHas to be one of \"ONT\",\"PB\", \"proto\",\"SLR\", \"ill\", \"miSeq\", \"AVITI\", \"hiSeq\", \"454\", \"GAII\",\"GAII_solexa\" or \"\"\n\n";
	}
}


sub checkAssmblGrp{
	my ($mapHR) = @_;
	my %map = %{$mapHR};
	
	my %memberAGs = %{$map{opt}{asGrpMemsHr}};
	foreach my $k (keys %memberAGs){
		my %curSeqTech; my %curSeqTechGen;
		my %curSeqTechSupp; my %curSeqTechGenSupp;
		foreach (@{$memberAGs{$k}}){  
			$curSeqTech{ $map{$_}{SeqTech} } ++;
			$curSeqTechGen{ is3rdGenSeqTech($map{$_}{SeqTech})  } ++;
			if ($map{$_}{SupportReads} ne ""){
				$map{$_}{SupportReads} =~ m/^([^:]+):/;
				$curSeqTechSupp{$1}++;
				$curSeqTechGenSupp{ is3rdGenSeqTech($1)  } ++;
			}
		}
		#DEBUG
		if ( 0  #just added to allow for || statements easily removable
			#|| scalar(keys(%curSeqTechGen)) > 1 ||  scalar(keys(%curSeqTechGenSupp)) > 1 
			|| ( exists($curSeqTechGenSupp{0}) && exists($curSeqTechGen{0}))  
			|| ( exists($curSeqTechGenSupp{1}) && exists($curSeqTechGen{1})) 
			|| ( exists($curSeqTechGenSupp{0}) && exists($curSeqTechGen{1}))
			){
			my %genNames = ("0"=>"2nd", "1"=>"3rd");	
			print STDERR "Fatal error in map:\nIt seems that assembly group $k contains a mix of second and third gen seq techs:\nPrimary input: ";
			foreach my $k(keys %curSeqTech){print "$k(N=$curSeqTech{$k}) " if ($k ne "");}print "\n  Generation : ";
			foreach my $k(keys %curSeqTechGen){print "$genNames{$k}(N=$curSeqTechGen{$k}) " if ($k != -1);}print "\nSuppl input: ";
			foreach my $k(keys %curSeqTechSupp){print "$k(N=$curSeqTechSupp{$k}) " if ($k ne "");}print "\n  Generation: ";
			foreach my $k(keys %curSeqTechGenSupp){print "$genNames{$k}(N=$curSeqTechGenSupp{$k}) " if ($k != -1);}print "\n";
			print "Assuming this is for hybrid assemblies, MATAFILER only supports 2nd gen for primary reads (e.g. illumina) and 3rd gen for suppl reads (e.g. PacBio). Will exit until fixed\n";
			die ;
		}
		
		
	}
}

sub readMap{
	my $inF = $_[0];
	my $Scnt = defined $_[1] ? $_[1] : 0;
	my %ret = defined $_[2] ? %{$_[2]} : (); 
	my %agBP = defined $_[3] ? %{$_[3]} : ();
	my $folderStrClassical = defined $_[4] ? $_[4] : -1;
	my $xtraColStr = defined $_[5] ? $_[5] : "";
#die "$folderStrClassical\n";
	my @order = exists $ret{opt}{smpl_order} ? @{$ret{opt}{smpl_order}} : ();
	my %oldAssmGrps = ();
	if (@order > 0){#something already exists
		my @asgps= keys (%agBP);
		$oldAssmGrps{$_}++ for (@asgps);
	}
	my $dirCol = -1; my $SmplPrefixCol = -1; #to find location of primary reads either has to be defined.
	my $smplCol = 0; my $rLenCol = -1; 
	my $SeqTech = -1; #mate = mate pairs; SLR = synthetic long reads, e.g. 10X, TELseq..
	my $SeqTechS = -1; 
	my $AssGroupCol = -1; my $EstCovCol = -1; my $MapGroupCol = -1; my $SupRdsCol = -1;my $ExcludeAssemble = -1;
	my $cut5pR1 = -1;my $cut5pR2 = -1; my $FamGroupCol = -1; my $firstXrdsRd = -1;my $firstXrdsWr = -1;
	my $ENA_DLcol=-1; my $SRA_DLcol=-1;
	#some global params
	my $dir2dirs = ""; #dir on file system, where all dirs specified in map can be found (enables different indirs with different mapping files)
	my $dir2out = "";
	my $baseID = ""; my $mocatFiltPath = "";
	my $inDirSet = 0;
	my $cnt = -1;
	my $illuminaClip ="";
	my $GlbTmpD = "";
	my $NodeTmpD = "";
	my $infFoldClass = -1;
	my $xtraCol = -1;
	my $relaxSmplID = 0;
	my @dir2dirsA;
	my %trackMGs; my %trackAGs; #hashes to track the last (final) sample in each mapgroup.. important to know this to check if assembly / mapping is done 
	my %memberMGs ; my %memberAGs ;
	my %trackDirs;
	my %trackPrefixs;
	
	my $DOWARN = 1;
	my $warnDeactivateMsg = "In case you want to continue, insert \"#WARNING OFF\" underneath the header of your map file.\n";
	
	#print $inF."\n";
	open I,"<$inF" or die "Can't open map: $inF\n";
	#AssGrps
	while (my $line = <I>){
		#use chomp and s/\R//g; to catch also windows formated lines..
		$cnt++;$line =~ s/\R//g;chomp $line;
		next if (length($line) <= 1);
		if ($line =~ m/^#/ && $cnt > 0 ){#check for ssome global parameters
			if ($line =~ m/^#DirPath\s(\S+).*$/){$dir2dirs = resolve_path($1); $dir2dirs.="/" unless ($dir2dirs=~m/\/$/); push(@dir2dirsA,$dir2dirs);}
			if ($line =~ m/^#OutPath\s+(\S+)/){
				$dir2out = resolve_path($1); $inDirSet=0;$dir2out .= "/" if ($dir2out !~ m/\/$/);
			}
			if ($line =~ m/^#RunID\s+(\S+)/){$baseID = $1; $inDirSet=0;}
			if ($line =~ m/^#mocatFiltPath\s+(\S+)/){$mocatFiltPath = resolve_path($1);}
			if ($line =~ m/^#illuminaClip\s+(\S+)/){$illuminaClip = resolve_path($1);}
			if ($line =~ m/^#NodeTmpDir\s+(\S+)/){$NodeTmpD = $1;}
			if ($line =~ m/^#GlobalTmpDir\s+(\S+)/){$GlbTmpD = $1;}
			if ($line =~ m/^#WARNING\sOFF/){$DOWARN=0;print "Warning: Deactivated Warnings while reading the map file! ..\n";}
			if ($line =~ m/^#RelaxSMPLID\sTRUE/){$relaxSmplID = 1;print "Relaxed SMPLIDs\n";}
			if (!$inDirSet && $dir2out ne "" && $baseID ne ""){
				$dir2out.=$baseID unless ($dir2out =~ m/$baseID[\/]*$/); $inDirSet =1;
				$dir2out .= "/" if ($dir2out !~ m/\/$/);
			}
			next;
		}
		my @spl = split(/\t/,$line,-1);
		if ($cnt == 0){
			#read column headers -> check for MATAFILER specific instructions
			#die "@spl\n";
			$smplCol = first_index { /^#SmplID$/ } @spl;
			$dirCol = first_index { /^Path$/ } @spl;
			$SmplPrefixCol = first_index { /^SmplPrefix$/ } @spl;
			$SeqTech = first_index { /^SeqTech$/ } @spl;
			$SeqTechS = first_index { /^SeqTechSingl$/ } @spl;
			$rLenCol = first_index { /^ReadLength$/ } @spl;
			$AssGroupCol = first_index { /^AssmblGrps$/ } @spl;
			$FamGroupCol = first_index { /^FamilyGrps$/ } @spl;
			$EstCovCol = first_index { /^EstCoverage$/ } @spl;
			$ENA_DLcol = first_index { /^ENAdownload$/ } @spl;
			$SRA_DLcol = first_index { /^SRAdownload$/ } @spl;
			$MapGroupCol = first_index { /^MapGrps$/ } @spl;
			$SupRdsCol = first_index { /^SupportReads$/ } @spl;
			$ExcludeAssemble = first_index { /^ExcludeAssembly$/ } @spl;
			$cut5pR2 = first_index { /^cut5PR2$/ } @spl;
			$cut5pR1 = first_index { /^cut5PR1$/ } @spl;
			$firstXrdsRd = first_index {/^firstXreadsRd/} @spl;
			$firstXrdsWr = first_index {/^firstXreadsWr/} @spl;
			if ($xtraColStr ne ""){
				$xtraCol = first_index { /$xtraColStr/ } @spl;
			}
			#die "MAP: $AssGroupCol\n";
			#die "Only \"Path\" or \"SmplPrefix\" can be defined in mapping file. Both is not supported.\n" if ($dirCol != -1 && $SmplPrefixCol != -1);
			die "Could not find \"#SmplID\" in input map\n" unless($smplCol> -1);
			die "Expected to find at least \"SmplPrefix\" or \"Path\" as column headers in .map\n" if ($dirCol == -1 && $SmplPrefixCol == -1);
			die "Either \"SmplPrefix\" or \"Path\" has to be second column in .map\n" unless ($dirCol == 1 || $SmplPrefixCol == 1);
			next;
		} #maybe later check for col labels etc
		my $smplPrefixUsed=0; my $samplePathUsed=0;
		$Scnt++;
		#die "1 $GlbTmpD 2 $NodeTmpD map\n";
		#die $spl[0]." ".$spl[1]."\n";
		die "inPath has to be set in mapping file!\n" if ($dir2dirs eq "");
		#die "$dir2out\n";
		die "Provide tag \"#OutPath\" in map!\n" if ($dir2out eq "");
		die "Provide tag \"#RunID\" in map!\n" if ($baseID eq "");
		die "Provide tag \"#DirPath\" in map!\n" if ($dir2dirs eq "");
		my $curSmp = $spl[$smplCol];
		
		die"Error in .map, found empty sample id on line $cnt,: \"$line\"\n" if ($curSmp eq "");
		my $altCurSmp = "";
		#print $curSmp." ";
		die "\"opt\" is reserved keyword and can't be used as sample name (.map line $cnt)\n" if ($curSmp eq "opt");
		die "\"altNms\" is reserved keyword and can't be used as samples\n" if ($curSmp eq "altNms");
		#die "\"smpl_order\" is reserved keyword and can't be used as samples\n" if ($curSmp eq "smpl_order");
		#die "\"totSmpls\" is reserved keyword and can't be used as samples\n" if ($curSmp eq "totSmpls");
		die "Double sample ID \"$curSmp\"\n$line\n" if (exists $ret{$curSmp});
		die "Can't use character \"$\" in sampleID: $curSmp\n" if ($curSmp =~ m/\$/);
		die "Can't use character \"_\" in sampleID: $curSmp\n" if ($curSmp =~ m/_/);
		die "Can't use character \",\" in sampleID: $curSmp\n" if ($curSmp =~ m/,/);
		die "Can't use character \"-\" in sampleID: $curSmp\n" if (!$relaxSmplID && $curSmp =~ m/-/);
		die "Recommended not to use numeric as first char \"0-9\" in sampleID: $curSmp\n" if (!$relaxSmplID && $curSmp =~ m/^[0-9]/);
		die "Use alphanumeric characters (a-zA-Z0-9.) for sampleID: $curSmp\n" if (!$relaxSmplID && $curSmp =~ m/[^a-zA-Z0-9\.]/);
		
		my $cdir = ""; 
		#read in the path to sample
		if ($dirCol >= 0 && @spl > $dirCol && $spl[$dirCol] ne ""){
			$cdir = resolve_path($spl[$dirCol]);
			$samplePathUsed=1;
			my $input_dir = _join_path($dir2dirs, $cdir);
			if ($DOWARN && exists($trackDirs{$input_dir}) ){ die "Warning: Found the sample path \"$cdir\" more than once. This would lead to using reads twice, aborting.\n $warnDeactivateMsg";}
			$trackDirs{$input_dir} = 1;
		}
		my $cdir2= $cdir;
		$ret{$curSmp}{dir} = $cdir;#this one should stay without a tag
		$ret{$curSmp}{rddir} = _join_path($dir2dirs, $cdir);
		$ret{$curSmp}{clip} = $illuminaClip;
		$ret{$curSmp}{rddir} .="/" unless ($ret{$curSmp}{rddir} =~ m/\/$/);
		#die "$ret{$curSmp}{rddir} $dirCol $cdir $curSmp\n $smplCol $dirCol\n";
		if ($SmplPrefixCol>=0 && @spl > $SmplPrefixCol && $spl[$SmplPrefixCol] ne ""){
			$cdir2 = $spl[$SmplPrefixCol];
			$ret{$curSmp}{prefix} = $cdir2;
			if ($DOWARN && exists($trackPrefixs{"$dir2dirs/$cdir2"}) ){ die "Warning: Found the sample path \"$cdir2\" more than once. This would lead to using reads twice, aborting.\n $warnDeactivateMsg";}
			$trackPrefixs{"$dir2dirs/$cdir2"} = 1;
			$smplPrefixUsed =1;
		} else {
			$ret{$curSmp}{prefix} = "";
		}
		$cdir2.="/" unless ($cdir2 =~ m/\/$/);
		
		#basic check that no twice usage..
		if ($DOWARN && $smplPrefixUsed && $samplePathUsed){die"Warning: in mapping file both \"SmplPrefix\" and \"Path\" are set for sample $curSmp!\nThis is not supported\n$warnDeactivateMsg";}

		
		#old MATAFILER versions, deprecated in MATAFILER
		if ($folderStrClassical== -1){
			if (-d  $dir2out.$cdir2 && !-d $dir2out.$curSmp){
				if ($infFoldClass==0){die "readMap: Inferring old/new folder structure failed, as both folders seem to be valid\n";}
				$ret{$curSmp}{wrdir} = $dir2out.$cdir2;
				$infFoldClass = 1;
			} else {
				if ($infFoldClass==1){die "readMap: Inferring old/new folder structure failed, as both folders seem to be valid\n";}
				$ret{$curSmp}{wrdir} = $dir2out.$curSmp."/";
				$infFoldClass = 0;
			}
		}elsif ($folderStrClassical == 1){
			$ret{$curSmp}{wrdir} = $dir2out.$cdir2;
		} else {
			$ret{$curSmp}{wrdir} = $dir2out.$curSmp."/";
		}
		#die "$ret{$curSmp}{wrdir}\n";
		
		#base info for sample on current map line
		$ret{$curSmp}{SmplID} = $curSmp;
		$ret{$curSmp}{mapFinSmpl} = $curSmp;
		$ret{$curSmp}{assFinSmpl} = $curSmp;
		#ONT,PB,proto,miSeq,GAII etc
		if ($SeqTech >= 0) { 
			my $RT=$spl[$SeqTech];checkSeqTech($RT);$ret{$curSmp}{SeqTech} = $RT; 
			if ($RT ne "" && $ret{$curSmp}{prefix} eq "" && $ret{$curSmp}{dir} eq ""){
				die "For sample $curSmp, found read tech, but no primary input file location given. This can lead to undescribed behaviour, please remove readTech entry (set to \"\"), before proceeding.\n";
			}
		} else {$ret{$curSmp}{SeqTech} = "ill";} #set by default to illumina
		if ($SeqTechS >= 0) { my $RT=$spl[$SeqTechS];checkSeqTech($RT);$ret{$curSmp}{SeqTechSingl} = $RT; } else {$ret{$curSmp}{SeqTechSingl} = "";}
		
		$ret{$curSmp}{hasPrimaryRds}= 1;$ret{$curSmp}{hasPrimaryRds} = 0 if ($ret{$curSmp}{prefix} eq "" && $ret{$curSmp}{dir} eq "");
		
		if ($rLenCol >= 0){$ret{$curSmp}{readLength} = $spl[$rLenCol];} else {$ret{$curSmp}{readLength} = 0;}
		if ($ExcludeAssemble >= 0){$ret{$curSmp}{ExcludeAssem} = $spl[$ExcludeAssemble];} else {$ret{$curSmp}{ExcludeAssem} = 0;}
		#cut first basepairs from each read.. (read1 and read2)
		
		#instructions to remove X bases at 5' read 2
		$ret{$curSmp}{cut5pR2} = 0;
		if ($cut5pR2 >= 0){
			if (defined ($spl[$cut5pR2]) && length($spl[$cut5pR2]) > 0){
				$ret{$curSmp}{cut5pR2} = $spl[$cut5pR2];
			}
		} 
		#instructions to remove X bases at 5' read 1
		$ret{$curSmp}{cut5pR1} = 0;
		if ($cut5pR1 >= 0){
			if (defined ($spl[$cut5pR1]) && length($spl[$cut5pR1]) > 0){
				$ret{$curSmp}{cut5pR1} = $spl[$cut5pR1];
			}
		} 
		
		#read only the first few reads in each sample..
		$ret{$curSmp}{firstXrdsRd} = 0;
		if ($firstXrdsRd >= 0){
			if (defined ($spl[$firstXrdsRd]) && length($spl[$firstXrdsRd]) > 0){
				$ret{$curSmp}{firstXrdsRd} = $spl[$firstXrdsRd];
			}
		} 
		#write only the first few reads in each sample..
		$ret{$curSmp}{firstXrdsWr} = 0;
		if ($firstXrdsWr >= 0){
			if (defined ($spl[$firstXrdsWr]) && length($spl[$firstXrdsWr]) > 0){
				$ret{$curSmp}{firstXrdsWr} = $spl[$firstXrdsWr];
			}
		} 
		
		if ($xtraCol != -1){$ret{$curSmp}{$xtraColStr} = $spl[$xtraCol];}
		if ($EstCovCol >= 0){$ret{$curSmp}{DoEstCoverage} = $spl[$EstCovCol];} else {$ret{$curSmp}{DoEstCoverage} = 0;}
		
		#download samples? this is the ENA/SRA ID
		if ($ENA_DLcol >= 0){$ret{$curSmp}{ENA_download} = $spl[$ENA_DLcol];} else {$ret{$curSmp}{ENA_download} = "";}
		if ($SRA_DLcol >= 0){$ret{$curSmp}{SRA_download} = $spl[$SRA_DLcol];} else {$ret{$curSmp}{SRA_download} = "";}

		#create artifical assmblgrp based on counts..
		my $curAG = $Scnt;
		#print "$curAG ";
		if ($AssGroupCol >= 0 && $spl[$AssGroupCol] ne ""){
			$curAG = $spl[$AssGroupCol];
			if ($DOWARN && $curAG !~ /\D/){die "Error in .map, sample $curSmp: AssmblGrps are not allowed to be purely numeric characters!\ncurrent AssmblGrps: $curAG\n\n";}
			
			$ret{$curSmp}{AssGroup} = $curAG ;
			if (!exists($agBP{$curAG}{CntAimAss})){$agBP{$curAG}{CntAimAss}=0;}
			$agBP{$curAG}{CntAimAss}++;
			$altCurSmp = $curSmp."M".$agBP{$curAG}{CntAimAss};
			$memberAGs{$curAG} = [] unless (exists $memberAGs{$curAG});
			$trackAGs{$curAG} = $curSmp;  push(@{$memberAGs{$curAG}},$curSmp);
			$agBP{$curAG}{prodRun} = "";
			if (exists($oldAssmGrps{$curAG})){
				die "Warning: assembly group \"$curAG\" seems to exist in one of the previous maps. This is currently not supported for MATAFILER, please make sure all assembly groups are uniquely named across maps!\n";
			}
			#print $agBP{$spl[$AssGroupCol]}{CntAimAss}. " :$spl[$AssGroupCol]\n" ;
		} else {
			if (exists ($trackAGs{$curAG})){ #should not exist!
				die "Numerical AssmblGrps defined in .map ($curAG)? Please amend as this is clashing with autogenerated IDs\n";
			}
			$ret{$curSmp}{AssGroup} = $curAG; $agBP{$curAG}{CntAimAss}=0;$agBP{$curAG}{prodRun} = "";
		}
		
		#die "$SupRdsCol\t@spl\n$spl[$SupRdsCol]\n";
		$agBP{$curAG}{SupportReads} = "" unless (exists($agBP{$curAG}{SupportReads}));
		if ($SupRdsCol >= 0 && $SupRdsCol < @spl) { 
			$ret{$curSmp}{SupportReads} = normaliseSupportReads($spl[$SupRdsCol], $dir2dirs);
			if ($ret{$curSmp}{SupportReads} ne "") {
				my ($support_technology) = parseSupportReads($ret{$curSmp}{SupportReads});
				checkSeqTech($support_technology, "Mapping file SupportReads for sample $curSmp");
				$agBP{$curAG}{SupportReads} .= "," if ($agBP{$curAG}{SupportReads} ne "");
				$agBP{$curAG}{SupportReads} .= $ret{$curSmp}{SupportReads};
			}
			#print "\n\n$ret{$curSmp}{SupportReads} \n\n";
			
		} else {
			$ret{$curSmp}{SupportReads} = "";
		}

		
		#print $ret{$curSmp}{AssGroup}."\n\n" if ($curSmp eq "NEGE");
		if ($MapGroupCol >= 0 && $spl[$MapGroupCol] ne ""){
			my $curM = "M_".$spl[$MapGroupCol];
			$ret{$curSmp}{MapGroup} = $curM;
			if (!exists($agBP{$curM}{CntAimMap})){$agBP{$curM}{CntAimMap}=0;}
			$agBP{$curM}{CntAimMap}++;
			$memberMGs{$curM} = [] unless (exists $memberMGs{$curM});
			$trackMGs{$curM} = $curSmp; push(@{$memberMGs{$curM}},$curSmp);
			#print $agBP{$spl[$MapGroupCol]}{CntAimMap}. " :$spl[$MapGroupCol]\n" ;
		} else {$ret{$curSmp}{MapGroup} = $Scnt; $agBP{$Scnt}{CntAimMap}=0;}
		
		if ($FamGroupCol >= 0 && scalar(@spl) > $FamGroupCol && $spl[$FamGroupCol] ne "" ){
			my $curF = $spl[$FamGroupCol];
			$ret{$curSmp}{FamGroup} = $curF;
			#if (!exists($agBP{$curF}{CntAimFam})){$agBP{$curF}{CntAimFam}=0;}
			$agBP{$curF}{CntAimFam}++;
			$memberMGs{$curF} = [] unless (exists $memberMGs{$curF});
			$trackMGs{$curF} = $curSmp; push(@{$memberMGs{$curF}},$curSmp);
		} else {$ret{$curSmp}{FamGroup} = ""; $agBP{$Scnt}{CntAimFam}=0;}
		
		if ($altCurSmp ne ""){$ret{altNms}{$altCurSmp} = $curSmp;} 
		push(@order,$curSmp);
		#print $spl[0]."\n";
	}
	close I; #done reading mapping file
	
	#insert final sample destination for all AssGroups and MapGroups
	foreach my $k (keys %memberMGs){
		foreach (@{$memberMGs{$k}}){  $ret{$_}{mapFinSmpl} = $trackMGs{$k};  }
	}
	foreach my $k (keys %memberAGs){
		foreach (@{$memberAGs{$k}}){  
			$ret{$_}{assFinSmpl} = $trackAGs{$k};  
			#insert members of each AG back into object..
			$ret{$_}{AG_members} = $memberAGs{$k};
		}
	}
	
	#my @forbiddenSmplIDs = qw(opt totSmpls smpl_order inDir outDir baseID mocatFiltPath);
	
	#die();
	#if (exists($ret{opt})) {die "Sample ID \"opt\" seems to have been used.. this is a reserved keyword, please rename sample\n";}
	$ret{opt}{folderStruct} = $infFoldClass;
	$ret{opt}{inDir} = join(",",@dir2dirsA) ;#if ($dir2dirs ne "");
	$ret{opt}{outDir} = $dir2out ;#if ($dir2out ne "");
	$ret{opt}{baseID} = $baseID ;#if ($baseID ne "");
	$ret{opt}{mocatFiltPath} = $mocatFiltPath ;#if ($mocatFiltPath ne "");
	$ret{opt}{GlbTmpD} =$GlbTmpD;
	$ret{opt}{NodeTmpD} =	$NodeTmpD ;
	$ret{opt}{totSmpls} = $Scnt;
	$ret{opt}{smpl_order} = \@order;
	#make this redundant at later point...
	#$ret{inDir} = join(",",@dir2dirsA) if ($dir2dirs ne "");
	#$ret{outDir} = $dir2out if ($dir2out ne "");
	#$ret{baseID} = $baseID if ($baseID ne "");
	#$ret{mocatFiltPath} = $mocatFiltPath if ($mocatFiltPath ne "");
	#@order = keys %agBP;die "@order\n";
	my $asGrpHr = emptyAssGrpsObj(\%agBP);
	$ret{opt}{asGrpHr} = $asGrpHr;
	$ret{opt}{asGrpMemsHr} = \%memberAGs;
	

	return (\%ret,$asGrpHr);
}





sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}

sub systemW{
	my ($cmddd) = $_[0];
	if ($cmddd eq ""){return;}
	my $killOnDead = 1;
	$killOnDead = $_[1] if (@_ > 1);
	$cmddd = "set -e\nulimit -c 0;\n$cmddd"; #make sure only single line excecuted.. (and not a huge core dump)
	my $stat= system('bash', '-o', 'pipefail', '-c', $cmddd);
	#can use $? also instead for status
	if ($stat){
		my $detail = $stat == -1
			? "could not execute: $!"
			: ($stat & 127)
				? "terminated by signal ".($stat & 127)
				: "exit code ".($stat >> 8);
		my $message = "system call \n$cmddd\nfailed ($detail).\n";
		die $message if ($killOnDead);
		warn $message;
	}
	#chomp $stat;
	return $stat;
}

sub readTabByKey{
	my ($inF) = @_;
	my %ret;
	my ($I,$OK) = gzipopen($inF,"tab file");
	my $maxTabs = 0;
	while (my $l = <$I>){
		chomp $l; my @tmp = split /\t/,$l;
		$ret {$tmp[0]} = $tmp[1];
		if (@tmp > $maxTabs){$maxTabs = @tmp;}
	}
	close $I;
	if ($maxTabs > 2){print "Warning in Tab reader: more than 2 columns were present ($maxTabs)\n";}
	return %ret;
}

sub writeFasta{
	my ($hr,$of) = ($_[0],$_[1]);
	my $maxFs = -1;
	$maxFs = $_[2] if (@_ > 2);
	my %FA = %{$hr};
	
	if ($maxFs <0){	$maxFs =  scalar(keys %FA); } #$maxFs+=1000;}
	#die "$maxFs\n";
	my $cnt=0;
	open O,">$of" or die "can't open out fasta $of\n";
	foreach my $k (keys %FA){
		$cnt++; 
		if ($k =~ m/^>/){
			print O $k."\n".$FA{$k}."\n";
		} else {
			print O ">".$k."\n".$FA{$k}."\n";
		}
		last if ( $cnt > $maxFs);
	}
	close O;
	#die $of;
}

1;
