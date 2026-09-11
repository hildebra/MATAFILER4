#use strict;
#use warnings;
use Mods::TamocFunc qw(sortgzblast uniq);
use Mods::GenoMetaAss qw(gzipopen);
use Mods::FuncTools qw(mergeBlastPair);
use Mods::IO_Tamoc_progs qw(getProgPaths);
use File::Basename qw(dirname);
use File::Path qw(make_path);
#use List::MoreUtils qw(uniq);
use strict; 
die "Usage: $0 <blast.gz> <gene-output> <category-output> [ABR-database-dir]\n" unless @ARGV >= 3;
my ($inputfile, $outputfile, $outputfilecats, $dbdir) = @ARGV;
$dbdir //= getProgPaths('ABRfors_path_DB');
$dbdir =~ s{/+$}{};
my $ardbfile = "$dbdir/ardb.tabs.parsed";
my $mapfile = "$dbdir/ardb_and_reforg_mapping";
my $besthitfile = "$dbdir/ardb_vs_reforg9f.overlap90shortest_famthres_or_symbol.sorted.besthit";
for my $required ($inputfile, $ardbfile, $mapfile, $besthitfile) {
	die "Required ABR input is missing or empty: $required\n" unless -s $required;
}



$inputfile = sortgzblast($inputfile);
 
# die $inputfile."\n";
 
 
my $outD = dirname($outputfilecats);
make_path($outD) unless -d $outD;
 
#reads ardb tabs
my %ssym = (); #cat 1
my %scat = (); #cat 2
my %sthres = (); #cat 3

#read CAT DB
open my $ardb_fh, '<', $ardbfile or die "Cannot open $ardbfile: $!\n";
while (<$ardb_fh>) {

    my $aLine = $_;
    chomp ($aLine);

    my @words = split (/\t/, $aLine);

	if (uc ($words [3]) eq "BACA") { $words [7] = "80"; } # assume/fix typo in ardb file
    $ssym {$words [0]}{uc ($words [2])} = "";
    $scat {$words [0]}{uc ($words [3])} = "";
    $sthres {$words [0]}{uc ($words [7])} = "";
}
close $ardb_fh or die "Cannot close $ardbfile: $!\n";

open my $map_fh, '<', $mapfile or die "Cannot open $mapfile: $!\n";
my %sym2drug = (); #specific drug resistance
while (<$map_fh>) {
    my $aLine = $_;
    chomp ($aLine);
    my @words = split (/\t/, $aLine);
    $sym2drug {$words [1]} = $words [3];
}
close $map_fh or die "Cannot close $mapfile: $!\n";


#read in ID cutoffs
open my $besthit_fh, '<', $besthitfile or die "Cannot open $besthitfile: $!\n";
while (<$besthit_fh>) {
	my $aLine = $_;
	chomp ($aLine);
	my @words = split (/\t/, $aLine);
	foreach my $asym (keys % {$ssym {$words [1]}}) {
		$ssym {$words [0]}{$asym} = "";
	}
	foreach my $acat (keys % {$scat {$words [1]}}) {
		$scat {$words [0]}{$acat} = "";
	}
	foreach my $athres (keys % {$sthres {$words [1]}}) {
		$sthres {$words [0]}{$athres} = "";
	}
}

close $besthit_fh or die "Cannot close $besthitfile: $!\n";

open my $gene_out, '>', $outputfile or die "Cannot open $outputfile: $!\n";
open my $cat_out, '>', $outputfilecats or die "Cannot open $outputfilecats: $!\n";
my ($blast_fh, $blast_ok) = gzipopen($inputfile, 'ABR blast input', 1);
die "Cannot open ABR blast input $inputfile\n" unless $blast_ok;
my $quOld = "";
my %wordv1; my %wordv2;my ($okhit,$retstr,$jnLine) ;
while (<$blast_fh>) {
	my $aLine = $_;
	chomp ($aLine);

	my @words = split (/\t/, $aLine);
	my $query = $words [0];
	$query =~ s/\/\d$//;
	
	if ($quOld eq "" ){$quOld = $query;
	} elsif ($quOld ne $query){  $quOld = $query;
		($okhit,$retstr,$jnLine) = workwords(\%wordv1,\%wordv2);
		
		if ($okhit){
			print {$gene_out} "$jnLine\n";
			print {$cat_out} "$retstr";
		}
		
		undef %wordv2; undef %wordv1;

	}
	
	#fill in array
	if (@words > 10 && $words[0] =~ m/2$/){
		$wordv2{$words [1]} = \@words;
	} elsif (@words > 10) {$wordv1{ $words [1]} = \@words;}
}

($okhit,$retstr,$jnLine) = workwords(\%wordv1,\%wordv2);

if ($okhit){
	print {$gene_out} "$jnLine\n";
	print {$cat_out} "$retstr";
}
close $blast_fh or die "Cannot finish reading $inputfile: $!\n";
close $gene_out or die "Cannot close $outputfile: $!\n";
close $cat_out or die "Cannot close $outputfilecats: $!\n";

open my $stone, '>', "$inputfile.stone" or die "Cannot create $inputfile.stone: $!\n";
close $stone or die "Cannot close $inputfile.stone: $!\n";




sub combineBlasts($ $){
	my ($wh1,$wh2) = @_;
	my %bl1 = %{$wh1}; my %bl2 = %{$wh2};
	my %ret;
	
	my @allKs = uniq ( keys %bl1, keys %bl2); #
	#die "@allKs\n";
	foreach my $k (@allKs){
		my $ex1 = exists ($bl1{$k});
		unless ($ex1 && exists ($bl2{$k}) ){
			if ($ex1){$ret{$k} = $bl1{$k};
			} else { $ret{$k} = $bl2{$k};}
			next;
		}
		$ret{$k} = mergeBlastPair($bl1{$k}, $bl2{$k});
	}
	
	return \%ret;
}

sub bestBlHit($){
	my ($hr) = @_;
	#my $emode=0; my $scomode = 1;
	my %blasts = %{$hr}; my $bestBit=0; my $bestk=""; my $bestScore=0;my $bestID=0; my $bestLen=0; my $bestIDever=0;
	my $k = "";
	foreach $k (keys %blasts){
	#print $k."\n";
		#my ($Query,$Subject,$id,$AlLen,$mistmatches,$gapOpe,$qstart,$qend,$sstart,$send,$eval,$bitSc) = ${$blRes{$k}}[11];
		#print $eval."\n";
		#sort by eval #changed from bestE -> bestScore
		#if ( ( ($emode && $bestE > $eval) || ($scomode && $bestScore< ($id * $AlLen)) ) 
		#			#&& ($eval <= $minBLE || $bitSc >= $minScore)
		#			&& ($quCovFrac == 0 || $AlLen > $DBlen{$Subject}*$quCovFrac) 
		#			&& (($fndCat || $noHardCatCheck) || exists $c2CAT{$Subject}) ) {
		#			#print "Y";
		
		#my $curScore = ${$blasts{$k}}[2] * ${$blasts{$k}}[3];
		#if ($bestScore < $curScore){$curScore = $bestScore;$bestk = $k;	}
		
		#	$bestAlLen=$AlLen;$bestE = $eval;#$bestQuery = $Query;
		
		#just sort by bitscore
		#print "@{$blasts{$k}}\n";
		#if ($bestBit < ${$blasts{$k}}[11]){
		#	$bestBit = ${$blasts{$k}}[11];  $bestk = $k;
		#}
		my $cID = ${$blasts{$k}}[2]; my $cLe = ${$blasts{$k}}[3]; my $cSc=${$blasts{$k}}[11];
		if (($bestID -5)< $cID){
			if ( ( $cID >= ($bestID *0.97) && $cID >= $bestIDever*0.9   && ( $cLe >= $bestLen * 1.15 ) )  ||  #length is just better (15%+)
					(   ($cID >= $bestID *1.03) && ( $cLe >= $bestLen * 0.9) ) ||#id is just better, while length is not too much off
					$cSc > $bestBit*0.8){   #convincing score
				$bestID = $cID;  $bestk = $k; $bestLen = $cLe; $bestBit = $cSc;
				if ($bestID > $bestIDever){$bestIDever=$bestID;}
			}
		}
	}
	if ($bestk eq ""){die "something went wrong with ABR blast:\n$k:$blasts{$k}\n";}
	return $blasts{$bestk};
}

sub workwords(){
	my ($wh1,$wh2) = @_;
	if (keys %{$wh2} != 0 && keys %{$wh1}==0){my $tmp = $wh1; $wh1 = $wh2; $wh2 = $tmp;}
	if (keys %{$wh1} == 0){return (0,"");}
	#1st combine scores
	my $whX = $wh1;
	if (keys %{$wh2} != 0){$whX = combineBlasts($wh1,$wh2);}
	#2nd: sort these, find best hit
	#print %{$whX}."\n";
	my $arBhit = bestBlHit($whX);
	#print "found\n";
	my @words = @{$arBhit};
	#contains already combined blast scores
	my $query = $words [0];
	my $subject = $words [1];
#	$qlen = $words [3];
#	$slen = $words [4];
	my $id = $words [2];
#	$bitscore = $words [5];
#	$all = $words [6];
	my $evalue = $words [10];

##	$overlap = $all / $llen;
	my @thress = sort {$a <=> $b} keys % {$sthres {$subject}};
	my $thres = $thress [0];
	my $okhit=0;my $retstr2 = ""; my $jnLine = "";
	if ($evalue < 0.00001 && $id >= $thres) {
		my %ldrugs = ();
		foreach my $sym (keys % {$ssym {$subject}}) {
			foreach my $drug (split (/\,/, $sym2drug {$sym})) {
				$ldrugs {$drug} = "";
			}
		}
		$okhit = 1;
		$jnLine = join ("\t",@words);
	#print "@words\n";
		$retstr2 = "$query\t$subject\t$id\t".join (",", keys % {$ssym {$subject}})."\t".join (",", keys % {$scat {$subject}})."\t".join (",", keys % ldrugs)."\n";
		#print FH2 "$aLine\n";
		#print FH3 "$query\t$subject\t$id\t".join (",", keys % {$ssym {$subject}})."\t".join (",", keys % {$scat {$subject}})."\t".join (",", keys % ldrugs)."\n";
	}
	return($okhit,$retstr2,$jnLine);
}
