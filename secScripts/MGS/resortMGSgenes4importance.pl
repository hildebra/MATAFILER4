#!/usr/bin/perl
#script that takes a selection of MGS genes (canopy format) and sorts them based on a) marker genes b) copy unmber c) overall occurence
#./resortMGSgenes4importance.pl /g/bork3/home/hildebra/data/SNP/GCs/DramaGCv5/ /g/bork3/home/hildebra/data/SNP/GCs/DramaGCv5/Binning/MetaBat/MB2.clusters.ext.can.Rhcl.mgs
use warnings;
use strict;
use Mods::GenoMetaAss qw(gzipopen);
use Mods::TamocFunc qw(readTabbed);
use Mods::math qw(meanArray medianArray);

sub evalCurMGS;
sub collectRequestedGenes;
sub coreRowFields;
sub readInformativeNT;
sub elapsedText;


#v0.1: adopt .core MGS files to get additional info for sorting genes by importance
#v0.11: 9.2.24: retain more genes/MGS
#v0.12: 11.2.24: adopted to weighted multiBin scores; more subs to make script more modifiable
#v0.13: flush the final MGS and handle marker-free groups safely
#v0.14: count distinct samples, use multicopy evidence, and make ranking deterministic
#v0.15: trim marker identifiers and remove an unused external-tool config dependency
#v0.16: 28.8.26: publish the sorted guide atomically; reject blank, header and
#	repeated MGS rows; index only the catalogue genes this MGS set references;
#	rank on precomputed keys; use the MGS observation table as the prevalence
#	prior and compare gene prevalence on MGS-restricted counts; rank genes
#	without catalogue occurrence last instead of first; rank on expected
#	informative nucleotides rather than on locus count alone
#v0.17: preview only the first few per-MGS detail rows, then report aggregate
#	ranking totals and elapsed time instead of writing thousands of near-identical lines
my $version = 0.17;

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Gene rejection thresholds.
# The MGS core table has already passed filterMB2.R, which keeps only
# MultiBin < 3 and MultiCopy/Occ < 0.1 (< 0.2 for marker genes). Thresholds
# here have to sit inside those bounds to have any effect at all.
# MultiBin is an integer count of the bins a gene is confidently assigned to
# (clusterMAGs.pl), so a core gene carries either 1 or 2: a gene at 2 recruits
# reads from a second MGS and would give a chimeric consensus.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
my $maxMultiBin = 2;		#reject genes confidently assigned to >= this many bins
my $maxCopyFraction = 0.05;	#reject above this fraction of multi-copy observations
my $minObsCopyCheck = 5;	#..but only with this many observations to judge on
my $minMultiCopyHits = 2;	#..and this much repeated duplication evidence
my $prevLowFactor = 0.25;	#prevalence window around the expected MGS prevalence
my $prevHighFactor = 4;
my $minExpectedPrevalence = 4;	#below this the prevalence window is not informative

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Informative sequence length.
# Downstream both caps (-presortGenes here, -treeLocusBudget at tree time) count
# loci, but tree resolution accrues per informative nucleotide: a 300bp and a
# 3000bp locus otherwise cost the same slot. Genes are therefore ranked on the
# usable sequence they are expected to contribute, prevalence x informative NT,
# using the same definition of "informative" as buildTree5.pl's
# informativeSequenceLength (unambiguous residues only, AA counted x3).
# The prevalence window above still applies first, so this only ever reorders
# genes that already have acceptable prevalence - length can never buy a locus
# past the callability requirement.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
my $rankOnInformativeNT = 1;	#0 falls back to ranking on prevalence alone
my $minInformativeNT = 150;	#reject loci with less usable sequence (0 disables)

#Core marker genes stay phylogenetically informative down to the strain level and
#are emitted ahead of everything else, so they are held to the chimera and
#paralogy checks (a chimeric marker is actively harmful) but not to the
#prevalence window, and only to a much lower length floor. This mirrors how
#filterMB2.R already treats markers: more tolerant, but never exempt from the
#copy-number and multi-bin checks.
my $markerMinInformativeNT = 60;	#length floor applied to marker genes
my $markerExemptPrevalence = 1;		#markers ignore the prevalence window

#Cross-taxon sharing: geneOcc (samples the gene cluster was assembled in, across
#all taxa) divided by the MGS's own prevalence. ~1 means the gene is private to
#this species; well above 1 means it is also found where the species is absent,
#which is what horizontal transfer, shared mobile elements and conserved domains
#all look like. Such genes are demoted behind every clean gene rather than
#rejected, so a sparse MGS with nothing better can still use them.
my $maxSharingRatio = 2;	#demote above this ratio (0 disables the demotion)

#The candidate pool handed downstream is budgeted in informative nucleotides
#rather than in loci, so a set of short genes no longer yields a thin alignment
#just because it filled the locus count. The locus floor keeps the pool safely
#above every configured downstream cap: strain_within.pl reads this guide with
#-presortGenes (default 1200, twice) and -outgroupReferenceGeneCap (default
#2500), so 5000 is 2x the largest of them. Raise it alongside those options.
my $candidateBudgetNT = 2_500_000;	#target cumulative informative NT per MGS (0 disables)
my $minCandidateLoci = 5000;	#never emit fewer ranked genes than this
my $maxCandidateLoci = 0;	#hard ceiling on emitted genes (0 = none)
my $downstreamLocusCap = 1200;	#-presortGenes default, reported on for visibility


# set up some base variables
die "Usage: $0 GC-dir MGS-file GTDB|FMG mode [cluster-ID]\n" unless @ARGV == 4 || @ARGV == 5;
my $GCd = $ARGV[0];
my $MGSfile = $ARGV[1];
my $useGTDBmg = $ARGV[2];
my $mode = $ARGV[3];	#accepted for caller compatibility (strain_within.pl), unused here
my $clusterID = @ARGV == 5 ? $ARGV[4] : 95;
die "cluster-ID must be between 1 and 100\n"
	unless $clusterID =~ /^\d+$/ && $clusterID >= 1 && $clusterID <= 100;
my $obsFile = $MGSfile; #$ARGV[3];# if (@ARGV > 3);
$obsFile =~ s/\.core$//; $obsFile.=".obs";
die "ARG 2 option has to be \"GTDB\" or \"FMG\"\n" unless ($useGTDBmg eq "GTDB" || $useGTDBmg eq "FMG");

#main output file. Written to a temporary and renamed only once the last MGS is
#flushed: a partially written guide is indistinguishable from a complete one to
#the -onlySubmit resume path in strain_within.pl, which only tests for size.
my $finout = "$MGSfile.srt";
my $tmpout = "$finout.tmp.$$";
my $published = 0;
#do not leave a half-written temporary behind for a later run to trip over
END { unlink $tmpout if (!$published && defined($tmpout) && -e $tmpout); }
if (-e $finout && -s $finout){
	#print $finout."\n"; exit(0);
	print "Overwriting $finout\n";
}


my $sortStarted = time;
print "\n--------------------------------------------------\nResorting MGS genes for importance in strain phylo ver $version\n--------------------------------------------------\n";
#my @FMG40 = ("COG0012","COG0016","COG0018","COG0048","COG0049","COG0052","COG0080","COG0081","COG0085","COG0087","COG0088","COG0090","COG0091","COG0092","COG0093","COG0094","COG0096","COG0097","COG0098","COG0099","COG0100","COG0102","COG0103","COG0124","COG0172","COG0184","COG0185","COG0186","COG0197","COG0200","COG0201","COG0202","COG0215","COG0256","COG0495","COG0522","COG0525","COG0533","COG0541","COG0552");
#my %FMG40 = map { $_ => 1 } @FMG40;


die "Can't find main infile $MGSfile\n" unless (-s $MGSfile);

#read MGS occurrence to understand distribution. This is the primary prevalence
#prior: it is measured per MGS, while the gene-derived estimate used as fallback
#is only inferred from whichever marker genes this MGS happens to carry.
my $MGSobs = {};
if (-s $obsFile) {
	$MGSobs = readTabbed($obsFile);
} else {
	warn "MGS occurrence table $obsFile is unavailable; falling back to the marker-gene prevalence estimate\n";
}

#load GTDB/FMG genes directly..
my $inMGFile="$GCd/FMG.subset.cats";
if ($useGTDBmg eq "GTDB"){
	$inMGFile="$GCd/GTDBmg.subset.cats";
}

#load list of reference marker genes (to mark these later as important genes)
my %MGset=();
open I,"<$inMGFile" or die "resortMGSgenes4importance.pl: Couldn't open $inMGFile\n";
my $totMGSgenes=0;
while (<I>){
	chomp;
	s/\r$//;
	my @spl1 = split /\t/, $_, -1;
	next unless defined($spl1[2]) && length($spl1[2]);
	my @spl2 = split /,/,$spl1[2];
	foreach my $gene (@spl2){
		$gene =~ s/^\s+|\s+$//g;
		next unless length($gene);
		$MGset{$gene} = 1;
		$totMGSgenes++;
	}
}
close I;
print STDERR "Loaded $totMGSgenes $useGTDBmg marker genes from $inMGFile\n";



#alt: go with compl.incompl.95.fna.clstr.idx to calc gene occurrences..
#Only genes this MGS set references are ever looked up, so the cheap pre-pass
#below keeps the rest of the gene catalogue out of memory entirely.
my $requested = collectRequestedGenes($MGSfile);
print STDERR "MGS file references " . scalar(keys %{$requested}) . " catalogue genes\n";
my %geneOcc;
my ($I,$ST) = gzipopen("$GCd/compl.incompl.$clusterID.fna.clstr.idx","gene cat index file");
#my $maxOcc = 0;
while (my $line = <$I>){
	next if (substr($line,0,1) eq "#"); #index header: "#Gene\tmembers"
	my $tab = index($line,"\t");
	next if ($tab < 1); #die $1." no tab char\n";
	my $gene = substr($line,0,$tab);
	next unless exists($requested->{$gene});
	my $members = substr($line,$tab+1); chomp $members;
	my %samples;
	for my $member (split /,/, $members) {
		my $sample = substr($member,0,1) eq ">" ? substr($member,1) : $member;
		my $sep = index($sample,"__");
		$sample = substr($sample,0,$sep) if ($sep >= 0);
		$samples{$sample} = undef if length($sample);
	}
	$geneOcc{$gene} = scalar(keys %samples);
	#$maxOcc = $spl[1] if ($maxOcc < $spl[1]);
}
close $I;
print STDERR "Counted gene occurrence for " . scalar(keys(%geneOcc)) . " genes\n";

#usable sequence per gene, so the downstream locus budget can be spent on the
#loci that actually carry alignable sequence
my $infoNT = $rankOnInformativeNT ? readInformativeNT($GCd,$clusterID,$requested) : {};
$requested = {};

#my $gene2taxF = createGene2MGS($MGSfile,$GCd);
#my %gen2Bin ; my %gene2COG;
#open I,"<$gene2taxF" or die "Cant' open $gene2taxF\n";
#while (<I>){
#	chomp; my @spl = split /\t/;
#	push (@{$gen2Bin{$spl[1]}}, $spl[0]);
#	$gene2COG{$spl[0]} = $spl[2] if (defined ($spl[2]));
#}
#close I;


#my $sortedOutFile = "$MGSfile.srt";
open O,">$tmpout" or die "Can't open temporary outfile $tmpout: $!\n";
open I,"<$MGSfile" or die "cant open infil $MGSfile\n";
my $cn=0; my $rowCnt=0; my $MGScnt = 0; my $geneCnt=0;
my $dupGenes=0; my $dupWarnLimit=5;
my $detailPreviewLimit = 5;
my $detailLinesReported = 0;
my ($summaryRankedGenes, $summaryMarkerGenes, $summarySharedDemotions,
	$summaryInformativeNT, $summaryCapInformativeNT) = (0, 0, 0, 0, 0);
my (@summaryEmittedLoci, @summaryMedianSharing, @summaryMedianInformativeNT);
my $curMGS=""; my %doneMGS;
my %occ; my %multiCp; my %markers; my %multiBin;
#foreach my $mg (keys %gen2Bin){
while (my $line = <I>){
	my $spl = coreRowFields($line);
	next unless defined $spl;
	$rowCnt++;
	my $MGS = $spl->[0];
	$curMGS = $MGS if ($curMGS eq "");

	if ($MGS ne $curMGS){
		#Grouping by MGS is assumed throughout: a second block for an already
		#flushed MGS would lose its best genes to the -presortGenes cap downstream.
		die "MGS $MGS reappears after it was completed; the MGS file must be grouped by MGS\n"
			if exists($doneMGS{$MGS});
		my $retS = evalCurMGS($MGS);
		print O $retS;
	}
	my $gene = $spl->[1];
	#push @genes,$gene;
	if (exists($multiBin{$gene})) {
		$dupGenes++;
		warn "Ignoring duplicate gene $gene in MGS $MGS\n"
			if ($dupGenes <= $dupWarnLimit);
		next;
	}
	if (exists($MGset{$gene})){
		$markers{$gene}=$spl->[2];
	} else {
		$occ{$gene} = $spl->[2];
	}
	$multiCp{$gene} = $spl->[3]; $multiBin{$gene} = $spl->[4];
	$cn ++;


}
print O evalCurMGS("") if $curMGS ne "";
close O or die "Can't close temporary outfile $tmpout: $!\n";
close I;
rename $tmpout, $finout or die "Can't publish $finout: $!\n";
$published = 1;
warn "Suppressed ".($dupGenes - $dupWarnLimit)." further duplicate-gene warnings "
	."($dupGenes total)\n" if ($dupGenes > $dupWarnLimit);
#Report a compact preview contract. The detailed values remain visible for the
#first few MGS, while catalogue-scale runs end with totals rather than thousands
#of lines that obscure the following workflow stages.
my $detailLinesOmitted = $MGScnt - $detailLinesReported;
print "MGS detail preview: $detailLinesReported/$MGScnt shown";
print "; $detailLinesOmitted omitted" if $detailLinesOmitted > 0;
print "\n";
my $medianEmitted = @summaryEmittedLoci ? medianArray(@summaryEmittedLoci) : 0;
my $medianSharing = @summaryMedianSharing ? medianArray(@summaryMedianSharing) : 0;
my $medianInformativeNT = @summaryMedianInformativeNT
	? medianArray(@summaryMedianInformativeNT) : 0;
print "Sorting summary: MGS=$MGScnt, input_rows=$rowCnt, unique_genes=$cn, "
	."ranked_genes=$summaryRankedGenes, emitted_genes=$geneCnt, "
	."marker_genes=$summaryMarkerGenes, sharing_demotions=$summarySharedDemotions\n";
print "Sorting summary: median_emitted_loci=$medianEmitted, "
	."median_sharing_ratio=".sprintf('%.2f', $medianSharing)
	.", median_informative_NT=$medianInformativeNT, informative_NT_emitted="
	."$summaryInformativeNT, informative_NT_at_first_$downstreamLocusCap="
	."$summaryCapInformativeNT, elapsed=".elapsedText(time - $sortStarted)."\n";
print "Saved in $finout\n";


exit(0);


#validate and split one MGS core row; returns undef for rows that carry no gene
sub coreRowFields{
	my ($line) = @_;
	$line =~ s/[\r\n]+$//;
	return undef if ($line =~ m/^\s*$/); #blank lines are not an error
	return undef if ($line =~ m/^#/); #commented line
	my @spl = split /\t/,$line;
	#header of an unfiltered .clusters table: "Bin Gene Occ MultiCopy MultiBin isMarkerGene"
	return undef if (defined($spl[0]) && $spl[0] eq "Bin");
	die "Undefine MGS on line $line\n" unless (defined($spl[0]) && $spl[0] ne "");
	die "Malformed MGS row: $line\n" unless @spl >= 6 && defined $spl[1] && length $spl[1];
	return \@spl;
}

#informative sequence length per catalogue gene, matching buildTree5.pl's
#informativeSequenceLength: unambiguous residues only. The protein catalogue is
#preferred because it holds the same information at a third of the I/O; its
#counts are converted to nucleotides the same way the tree-time selection does.
sub readInformativeNT{
	my ($GCd,$clusterID,$wanted) = @_;
	my $protF = "$GCd/compl.incompl.$clusterID.prot.faa";
	my $nucF = "$GCd/compl.incompl.$clusterID.fna";
	my ($seqF,$isAA) = (-s $protF) ? ($protF,1) : ((-s $nucF) ? ($nucF,0) : (undef,0));
	unless (defined $seqF){
		warn "No gene catalogue sequences at $protF or $nucF; "
			."ranking on prevalence alone and skipping the informative-length filter\n";
		return {};
	}
	print STDERR "Measuring informative sequence length from $seqF\n";
	my ($IN,$OK) = gzipopen($seqF,"gene catalogue sequence file");
	my %infoNT; my $cur; my $cnt=0;
	while (my $line = <$IN>){
		if (substr($line,0,1) eq ">"){
			$infoNT{$cur} = $isAA ? $cnt*3 : $cnt if defined $cur;
			my $id = substr($line,1);
			$id =~ s/[\r\n]+$//; $id =~ s/\s.*$//;
			$cur = exists($wanted->{$id}) ? $id : undef;
			$cnt = 0;
		} elsif (defined $cur){
			#same character classes informativeSequenceLength() keeps after uc()
			$cnt += $isAA ? ($line =~ tr/ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy//)
				: ($line =~ tr/ACGTUacgtu//);
		}
	}
	$infoNT{$cur} = $isAA ? $cnt*3 : $cnt if defined $cur;
	close $IN;
	my $found = scalar(keys %infoNT); my $asked = scalar(keys %{$wanted});
	print STDERR "Informative length known for $found/$asked requested genes\n";
	warn "Only $found of $asked MGS genes were found in $seqF; the rest rank "
		."neutrally on length\n" if ($asked && $found < $asked * 0.9);
	return \%infoNT;
}

#which catalogue genes does this MGS file ask about? Reading the guide once up
#front is far cheaper than indexing the whole gene catalogue.
sub collectRequestedGenes{
	my ($inFile) = @_;
	my %wanted;
	open my $in,"<",$inFile or die "cant open infil $inFile\n";
	while (my $line = <$in>){
		my $spl = coreRowFields($line);
		next unless defined $spl;
		$wanted{$spl->[1]} = undef;
	}
	close $in;
	return \%wanted;
}


sub evalCurMGS{
	#this routine decides which MGS genes (already pre-filtered for core genes) will be handed on to strain phylo construction.. should be "certain" cutoffs for removing genes (intra phylo will do another round of filtering)
	my ($MGS) = @_;
	my $mrkCnt=0;
	my @all_genes = (keys(%markers), keys(%occ));
	my $obsOf = sub { my ($g) = @_; return exists($markers{$g}) ? $markers{$g} : $occ{$g}; };

	#Expected prevalence on the MGS-restricted scale of the .core Occ column.
	#The MAG observation count is measured per MGS and is the better prior; the
	#marker median only estimates it from whichever markers this MGS carries.
	my $recordedObs = $MGSobs->{$curMGS};
	my $expected_occ;
	if (defined($recordedObs) && $recordedObs =~ /^\d+(?:\.\d+)?$/ && $recordedObs > 0){
		$expected_occ = $recordedObs;
	} else {
		my @marker_obs = grep { $_ > 0 } map { $obsOf->($_) } keys %markers;
		if (@marker_obs >= 3){
			$expected_occ = medianArray(@marker_obs);
		} else {
			my @known_obs = grep { $_ > 0 } map { $obsOf->($_) } @all_genes;
			$expected_occ = @known_obs ? medianArray(@known_obs) : 1;
		}
	}

	#Cross-taxon sharing. geneOcc counts every sample the gene cluster was
	#assembled in, including samples without this MGS, while $expected_occ is how
	#many this species occupies. Their ratio is therefore ~1 for a gene private to
	#this species and rises the more it is found where the species is absent: the
	#signature of a conserved domain, a shared mobile element, or recent transfer.
	#Needs no annotation, so it also works where eggNOG is missing or wrong.
	my $sharingOf = sub {
		my ($g) = @_;
		return undef unless defined($geneOcc{$g}) && $expected_occ > 0;
		return $geneOcc{$g} / $expected_occ;
	};

	#Usable sequence per gene. A gene the catalogue has no sequence for ranks
	#neutrally on this axis rather than being punished for a catalogue gap.
	my @known_nt = grep { $_ > 0 } map { $infoNT->{$_} // 0 } @all_genes;
	my $median_nt = @known_nt ? medianArray(@known_nt) : 1;
	my $ntOf = sub {
		my ($g) = @_;
		return defined($infoNT->{$g}) ? $infoNT->{$g} : $median_nt;
	};

	#One pass per gene: the comparator used to recompute all of this on every
	#single comparison, which dominated the runtime for large MGS.
	my %copyFrac; my %rankKey; my $sharedCnt = 0; my @sharingSeen;
	for my $gene (@all_genes){
		my $observations = $obsOf->($gene);
		$copyFrac{$gene} = ($observations && $observations > 0)
			? ($multiCp{$gene} || 0) / $observations : 0;
		#a gene without catalogue occurrence carries no evidence at all, so it
		#ranks behind every measured gene instead of tying with a typical one
		my $unmeasured = defined($geneOcc{$gene}) ? 0 : 1;
		#Demote genes that look shared across taxa behind every clean gene. This
		#leads the ranking because a shared locus is wide for the wrong reason:
		#it recruits reads from another species and yields false strain variation,
		#so it must not ride high on the breadth term below. Demotion rather than
		#rejection keeps it available to a sparse MGS that has nothing better.
		#A small prevalence makes the ratio too noisy to act on, so the same
		#guard as the prevalence window applies.
		my $sharing = $sharingOf->($gene);
		push @sharingSeen, $sharing if defined $sharing;
		my $shared = ($maxSharingRatio > 0 && defined($sharing)
			&& $expected_occ >= $minExpectedPrevalence
			&& $sharing > $maxSharingRatio) ? 1 : 0;
		$sharedCnt += $shared;
		#Expected informative nucleotides this locus contributes to the alignment.
		#A gene that is both widely distributed across the MGS and long is the
		#most useful for tracking strains, so this leads the ranking. Copy
		#fraction and multi-bin are already enforced as hard rejections above;
		#using them again as dominant sort keys would let a negligible paralogy
		#difference outrank a far more informative locus, so they follow.
		my $yield = ($observations || 0) * $ntOf->($gene);
		$rankKey{$gene} = [ $shared, -$yield, $copyFrac{$gene},
			-($observations || 0), $multiBin{$gene} // 0, $unmeasured,
			$sharing // 0 ];	#unmeasured genes are already separated by key 5
	}
	my $medSharing = @sharingSeen ? medianArray(@sharingSeen) : 0;

	my $reject_gene = sub {
		my ($gene) = @_;
		my $observations = $obsOf->($gene);
		my $isMarker = exists($markers{$gene});
		# Chimera and paralogy checks apply to markers too: a marker that
		# recruits reads from a second MGS corrupts every tree it enters.
		return 1 if defined($multiBin{$gene}) && $multiBin{$gene} >= $maxMultiBin;
		# Repeated duplication across several MAGs is paralogy evidence.
		return 1 if $observations >= $minObsCopyCheck
			&& ($multiCp{$gene} || 0) >= $minMultiCopyHits
			&& $copyFrac{$gene} > $maxCopyFraction;
		# Only reject gross prevalence mismatches when enough samples exist.
		# Both sides are MGS-restricted counts: testing catalogue-wide breadth
		# against an MGS prevalence compares two different quantities.
		if (!($isMarker && $markerExemptPrevalence)
				&& $expected_occ >= $minExpectedPrevalence && defined($observations)) {
			return 1 if $observations < $prevLowFactor * $expected_occ
				|| $observations > $prevHighFactor * $expected_occ;
		}
		# Too little alignable sequence to contribute anything but noise. Only
		# genes the catalogue actually has a sequence for can fail this.
		my $lenFloor = $isMarker ? $markerMinInformativeNT : $minInformativeNT;
		if ($lenFloor > 0 && defined($infoNT->{$gene})) {
			return 1 if $infoNT->{$gene} < $lenFloor;
		}
		return 0;
	};
	my $rank_genes = sub {
		my ($left, $right) = @_;
		my ($l,$r) = ($rankKey{$left}, $rankKey{$right});
		return $l->[0] <=> $r->[0]
			|| $l->[1] <=> $r->[1]
			|| $l->[2] <=> $r->[2]
			|| $l->[3] <=> $r->[3]
			|| $l->[4] <=> $r->[4]
			|| $l->[5] <=> $r->[5]
			|| $l->[6] <=> $r->[6]
			|| $left cmp $right;
	};

	#markers first: they keep loci comparable across samples. The order within
	#each group is the order the -presortGenes cap downstream will honour.
	my @ranked;
	my @mrks = sort { $rank_genes->($a, $b) } keys %markers;
	for my $gn (@mrks) {
		next if $reject_gene->($gn);
		push @ranked, $gn;
		$mrkCnt++;
	}
	my @srtedGenes = sort { $rank_genes->($a, $b) } keys %occ;
	for my $gn (@srtedGenes) {
		next if $reject_gene->($gn);
		push @ranked, $gn;
	}

	#spend the candidate budget in informative nucleotides, not in locus slots
	my @finalGs; my $budgetNT = 0; my $capNT = 0;
	for my $gn (@ranked) {
		last if ($maxCandidateLoci && @finalGs >= $maxCandidateLoci);
		last if ($candidateBudgetNT && $budgetNT >= $candidateBudgetNT
			&& @finalGs >= $minCandidateLoci);
		push @finalGs, $gn;
		$budgetNT += $ntOf->($gn);
		$capNT = $budgetNT if (@finalGs <= $downstreamLocusCap);
	}

	my $medMBi = @all_genes ? medianArray(map { $multiBin{$_} } @all_genes) : 0;
	my $avgMBi = @all_genes ? meanArray([map { $multiBin{$_} } @all_genes]) : 0;
	my $detailLine = "${curMGS} (".scalar(keys %multiBin)."):: "
		.scalar(@finalGs) ."/". scalar(@ranked) ." genes, $mrkCnt markerGs used"
		. ", expected prevalence: " . int($expected_occ*100)/100
		. ", median sharing ratio: " . int($medSharing*100)/100
		. " ($sharedCnt demoted), median informative NT: $median_nt"
		. ", informative NT emitted/at cap $downstreamLocusCap: $budgetNT/$capNT"
		. ", median/avg multiBin $medMBi/" . int($avgMBi*100)/100;
	if ($detailLinesReported < $detailPreviewLimit) {
		print "$detailLine\n";
		$detailLinesReported++;
	}
	$summaryRankedGenes += scalar(@ranked);
	$summaryMarkerGenes += $mrkCnt;
	$summarySharedDemotions += $sharedCnt;
	$summaryInformativeNT += $budgetNT;
	$summaryCapInformativeNT += $capNT;
	push @summaryEmittedLoci, scalar(@finalGs);
	push @summaryMedianSharing, $medSharing;
	push @summaryMedianInformativeNT, $median_nt;

	$geneCnt+= scalar(@finalGs);

	#my @finalGs = @finalList{@keys};
	#die "\n@finalGs\n";

	my $retStr= $curMGS."\t".join(",", @finalGs)."\n";


	$doneMGS{$curMGS} = 1;
	$curMGS = $MGS;
	%occ = (); %multiCp= (); %markers= (); %multiBin= ();
	$MGScnt++;
	return $retStr;
}

sub elapsedText {
	my ($seconds) = @_;
	$seconds = int($seconds || 0);
	my $hours = int($seconds / 3600);
	my $minutes = int(($seconds % 3600) / 60);
	my $remaining = $seconds % 60;
	return $hours ? "${hours}h${minutes}m${remaining}s"
		: $minutes ? "${minutes}m${remaining}s" : "${remaining}s";
}
