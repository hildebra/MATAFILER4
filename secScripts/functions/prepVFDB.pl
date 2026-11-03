#!/usr/bin/env perl
#prepares the VFDB annotation table (VF.tab) used by parseBlastFunct2.pl for -functDB VFA / VFB (and legacy VDB)
#usage: prepVFDB.pl <VFDB dir> [out file, default <VFDB dir>/VF.tab]
#expects VFDB_setA_pro.fas and VFDB_setB_pro.fas (optionally .gz) from http://www.mgc.ac.cn/VFs/download.htm
#VFDB protein header format:
#>VFG000076(gb|NP_460360) (ssaQ) type III secretion system protein SsaQ [TTSS (SPI-2 encode) (VF0036) - Effector delivery system (VFC0086)] [Salmonella enterica ...]
#VF.tab columns (tab separated, no header):
# 0 DB sequence ID (first word of header, as reported by diamond)
# 1 L0: gene_VFID      (e.g. ssaQ_VF0036)
# 2 description        (protein description [organism])
# 3 L1: VFID_VFname    (e.g. VF0036_TTSS_(SPI-2_encode))
# 4 L2: VFCID_category (e.g. VFC0086_Effector_delivery_system)
# 5 organism
# 6 set                (A = core, experimentally verified; B = full dataset only)
use warnings;
use strict;
use Mods::GenoMetaAss qw(gzipopen);

die "usage: $0 <VFDB dir> [out VF.tab]\n" if (@ARGV < 1);
my $dir = $ARGV[0]; $dir =~ s/\/+$//;
my $outF = @ARGV > 1 ? $ARGV[1] : "$dir/VF.tab";

sub clean{ #names become rtk features: no hierarchy (;), AND (,) or OR (|) separators, no whitespace
	my ($s) = @_;
	$s = "" unless defined $s;
	$s =~ s/^\s+|\s+$//g;
	$s =~ s/[\s;,|\/]+/_/g;
	$s =~ s/[^A-Za-z0-9_.()\-]//g;
	return $s;
}

sub readSet{
	my ($f) = @_;
	my %ret;
	my ($I,$OK) = gzipopen($f,"VFDB fasta",1);
	while (my $l = <$I>){
		next unless ($l =~ m/^>/);
		$l =~ s/[\r\n]+$//;
		$l =~ m/^>(\S+)/; $ret{$1} = $l;
	}
	close $I;
	return \%ret;
}

my ($fA) = grep { -e $_ } ("$dir/VFDB_setA_pro.fas","$dir/VFDB_setA_pro.fas.gz");
my ($fB) = grep { -e $_ } ("$dir/VFDB_setB_pro.fas","$dir/VFDB_setB_pro.fas.gz");
die "Can't find VFDB_setB_pro.fas(.gz) in $dir\n" unless (defined $fB);
die "Can't find VFDB_setA_pro.fas(.gz) in $dir\n" unless (defined $fA);
my %A = %{readSet($fA)};
my %B = %{readSet($fB)};
foreach my $id (keys %A){ $B{$id} = $A{$id} unless (exists $B{$id}); } #set A should be a subset of B; be safe

my $unparsed = 0; my $cnt = 0; my %VFs; my %VFCs;
my $tmpF = "$outF.tmp.$$"; #VFA and VFB DB-prep jobs may run concurrently
open my $O, ">", $tmpF or die "Can't write $tmpF\n";
foreach my $id (sort keys %B){
	my $h = $B{$id};
	my ($gene,$desc,$vfN,$vfID,$vfcN,$vfcID,$org) = ("","","","","","","");
	#VF/VFC names cannot contain brackets: a "[2Fe-2S]" in the description must not start the VF bracket
	if ($h =~ m/^>\S+\s+\(([^)]*)\)\s+(.*?)\s*\[([^\[\]]+?)\s+\((VF\d+)\)(?:\s+-\s+([^\[\]]+?)\s+\((VFC\d+)\))?\]\s*\[([^\]]*)\]\s*$/){
		($gene,$desc,$vfN,$vfID,$vfcN,$vfcID,$org) = ($1,$2,$3,$4,$5 // "",$6 // "",$7);
	} else {
		$unparsed++;
		$h =~ m/^>\S+\s*(.*)$/; $desc = $1;
	}
	$gene = $id if ($gene eq "" || $gene eq "-");
	my $L0 = clean($gene) . ($vfID ne "" ? "_$vfID" : "");
	my $L1 = $vfID ne "" ? clean("${vfID}_$vfN") : "VF_unassigned";
	my $L2 = $vfcID ne "" ? clean("${vfcID}_$vfcN") : "VFC_unclassified";
	$desc =~ s/\t/ /g; $org =~ s/\t/ /g;
	print $O join("\t",$id,$L0,"$desc [$org]",$L1,$L2,$org,(exists $A{$id} ? "A" : "B"))."\n";
	$VFs{$L1}=1; $VFCs{$L2}=1; $cnt++;
}
close $O;
rename $tmpF, $outF or die "Can't move $tmpF to $outF\n";
print "prepVFDB: wrote $cnt entries (".scalar(keys %A)." in set A) to $outF; ".scalar(keys %VFs)." VFs, ".scalar(keys %VFCs)." VF categories\n";
print "prepVFDB: WARNING $unparsed headers did not match the expected VFDB format (assigned to VF_unassigned/VFC_unclassified)\n" if ($unparsed);
