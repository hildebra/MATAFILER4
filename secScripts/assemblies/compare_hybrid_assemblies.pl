#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my ($final, $output, $help, @preassemblies);
Getopt::Long::Configure(qw(no_auto_abbrev no_ignore_case));
GetOptions('final=s'=>\$final, 'preassembly=s@'=>\@preassemblies,
	'output=s'=>\$output, 'help|h'=>\$help) or die usage();
if ($help) { print usage(); exit 0; }
die usage('--final, at least one --preassembly, and --output are required')
	unless defined($final) && @preassemblies && defined($output);
die usage('unexpected positional arguments') if @ARGV;

my @pre_lengths;
my ($pre_gc, $pre_bases) = (0,0);
for my $path (@preassemblies) { my ($l,$gc,$bases)=fasta_stats($path); push @pre_lengths,@$l; $pre_gc+=$gc; $pre_bases+=$bases; }
my ($final_lengths,$final_gc,$final_bases)=fasta_stats($final);
my $pre = metrics(\@pre_lengths,$pre_gc,$pre_bases);
my $fin = metrics($final_lengths,$final_gc,$final_bases);
open my $out, '>', $output or die "Cannot write '$output': $!\n";
print {$out} "metric\tpreassemblies\thybrid_final\tdelta\tratio_final_to_pre\n";
for my $metric (qw(contigs total_bp N50 N90 longest GC_percent)) {
	my ($a,$b)=($pre->{$metric},$fin->{$metric}); my $delta=$b-$a; my $ratio=$a ? $b/$a : 0;
	print {$out} join("\t",$metric,fmt($a),fmt($b),fmt($delta),sprintf('%.4f',$ratio)),"\n";
}
print {$out} "source_preassembly_count\t",scalar(@preassemblies),"\t1\t\t\n";
close $out or die "Cannot close '$output': $!\n";
print "Comparative assembly report written to $output\n";

sub usage { my($e)=@_; return (defined($e)?"Error: $e\n\n":'')."Usage: compare_hybrid_assemblies.pl --final FILE --preassembly FILE [--preassembly FILE ...] --output FILE\n"; }
sub fasta_stats {
	my($path)=@_; open my $fh,'<',$path or die "Cannot open '$path': $!\n"; my(@lengths,$seq); $seq=''; my($gc,$bases)=(0,0);
	while(my $line=<$fh>){$line=~s/[\r\n]+$//; if($line=~/^>/){if(length $seq){push @lengths,length($seq);$gc+=($seq=~tr/GCgc//);$bases+=length($seq)}$seq='';}else{$seq.=$line}}
	if(length $seq){push @lengths,length($seq);$gc+=($seq=~tr/GCgc//);$bases+=length($seq)} close $fh; die "No sequences in '$path'\n" unless @lengths; return(\@lengths,$gc,$bases);
}
sub nx { my($l,$fraction)=@_; my$total=0;$total+=$_ for @$l;my$target=$total*$fraction;my$c=0;for my$x(sort{$b<=>$a}@$l){$c+=$x;return$x if$c >= $target}return 0; }
sub metrics { my($l,$gc,$bases)=@_; my$longest=0; for my $length (@$l) { $longest=$length if $length>$longest; } return {contigs=>scalar(@$l),total_bp=>$bases,N50=>nx($l,.5),N90=>nx($l,.9),longest=>$longest,GC_percent=>$bases?100*$gc/$bases:0}; }
sub fmt { my($v)=@_; return int($v)==$v ? int($v) : sprintf('%.4f',$v); }
