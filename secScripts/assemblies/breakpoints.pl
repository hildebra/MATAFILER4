#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use Mods::GenoMetaAss qw(gzipopen gzipwrite);

# Detect well-supported assembly breakpoints from bedGraph coverage. Short
# high-coverage islands inside a low-coverage run are bridged before testing
# the candidate, which makes the result less sensitive to individual noisy
# mappings. A breakpoint is reported only when both flanks have real support.
my ($assembly, $coverage, $output, $help);
my $break_depth = 0.10;
my $min_length = 100;
my $smooth_gap = 100;
my $flank_length = 500;
my $min_flank_depth = 1;
my $max_flank_fraction = 0.10;
Getopt::Long::Configure(qw(no_auto_abbrev no_ignore_case));
GetOptions(
	'assembly=s' => \$assembly,
	'coverage=s' => \$coverage,
	'output=s' => \$output,
	'breakpoint-depth=f' => \$break_depth,
	'min-breakpoint-length=i' => \$min_length,
	'smooth-gap=i' => \$smooth_gap,
	'flank-length=i' => \$flank_length,
	'min-flank-depth=f' => \$min_flank_depth,
	'max-flank-fraction=f' => \$max_flank_fraction,
	'help|h' => \$help,
) or die usage();
if ($help) { print usage(); exit 0; }
die usage('unexpected positional arguments: '.join(' ', @ARGV)) if @ARGV;
die usage('--assembly, --coverage and --output are required')
	unless defined($assembly) && defined($coverage) && defined($output);
$output .= '.gz' unless $output =~ /\.gz$/;
die "Invalid breakpoint parameters\n" unless $break_depth >= 0 && $min_length > 0
	&& $smooth_gap >= 0 && $flank_length > 0 && $min_flank_depth >= 0
	&& $max_flank_fraction >= 0;

my ($lengths, $order) = fasta_lengths($assembly);
my %raw = map { $_ => [] } @$order;
my ($fh) = gzipopen($coverage, 'mapping coverage', 1);
my ($line_no, $matched, $unknown) = (0, 0, 0);
while (my $line = <$fh>) {
	$line =~ s/[\r\n]+$//;
	next if $line eq '' || $line =~ /^(?:#|track\b|browser\b)/;
	$line_no++;
	my @f = split /\s+/, $line;
	die "Malformed coverage line $line_no: $line\n" unless @f >= 4;
	my ($id, $start, $end, $depth) = @f[0..3];
	$id = canonical_id($id);
	die "Invalid coverage line $line_no: $line\n"
		unless $start =~ /^\d+$/ && $end =~ /^\d+$/ && $end > $start
		&& $depth =~ /^(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/ && $depth >= 0;
	if (!exists $lengths->{$id}) { $unknown++; next; }
	next if $start >= $lengths->{$id};
	$end = $lengths->{$id} if $end > $lengths->{$id};
	push @{$raw{$id}}, [0+$start, 0+$end, 0+$depth];
	$matched++;
}
close $fh or die "Cannot close coverage '$coverage': $!\n";
die "Coverage '$coverage' has no intervals matching the assembly\n" unless $matched;
warn "Ignored $unknown coverage interval(s) for unknown contigs\n" if $unknown;

my $out = gzipwrite($output, 'breakpoint TSV');
print {$out} join("\t", qw(contig start end length mean_depth left_depth right_depth)), "\n";
my ($reported, $reported_bases, $candidate_count) = (0, 0, 0);
for my $id (@$order) {
	my $runs = complete_runs($id, $lengths->{$id}, $raw{$id});
	my @low = map { [$_->[0], $_->[1]] } grep { $_->[2] <= $break_depth } @$runs;
	my @smoothed;
	for my $low (@low) {
		if (@smoothed && $low->[0] - $smoothed[-1][1] <= $smooth_gap) {
			$smoothed[-1][1] = $low->[1];
		} else { push @smoothed, [@$low]; }
	}
	for my $candidate (@smoothed) {
		my ($start, $end) = @$candidate;
		next if $end - $start < $min_length;
		$candidate_count++;
		# Terminal low-coverage tails are not internal contig breakpoints because
		# they cannot be supported by sequence on both sides.
		next if $start == 0 || $end == $lengths->{$id};
		my $left_start = $start > $flank_length ? $start - $flank_length : 0;
		my $right_end = $end + $flank_length < $lengths->{$id}
			? $end + $flank_length : $lengths->{$id};
		my $candidate_depth = mean_depth($runs, $start, $end);
		my $left_depth = mean_depth($runs, $left_start, $start);
		my $right_depth = mean_depth($runs, $end, $right_end);
		next if $left_depth < $min_flank_depth || $right_depth < $min_flank_depth;
		my $supported_limit = $max_flank_fraction
			* ($left_depth < $right_depth ? $left_depth : $right_depth);
		next if $candidate_depth > $break_depth || $candidate_depth > $supported_limit;
		print {$out} join("\t", $id, $start, $end, $end-$start,
			map { sprintf('%.4f', $_) } ($candidate_depth, $left_depth, $right_depth)), "\n";
		$reported++; $reported_bases += $end - $start;
	}
}
close $out or die "Cannot close breakpoint TSV '$output': $!\n";
print "Breakpoint detection summary\n",
	"  Low-coverage candidates: $candidate_count\n",
	"  Supported breakpoints:   $reported\n",
	"  Breakpoint bases:         $reported_bases\n",
	"  Output TSV:               $output\n";

sub usage {
	my ($error) = @_;
	return (defined($error) ? "Error: $error\n\n" : '').<<'USAGE';
Usage: breakpoints.pl --assembly FILE --coverage FILE --output FILE [options]
  --output FILE                   Gzipped TSV output (.gz is added if omitted)
  --breakpoint-depth FLOAT       Low-depth cutoff [0.1]
  --min-breakpoint-length INT    Minimum low region [100]
  --smooth-gap INT               Bridge noisy high-depth islands up to bp [100]
  --flank-length INT             Support window on each side [500]
  --min-flank-depth FLOAT        Minimum mean depth on each side [1]
  --max-flank-fraction FLOAT     Maximum break/flank depth ratio [0.1]
USAGE
}

sub canonical_id { my ($v)=@_; $v =~ s/^>//; $v =~ s/^\s+|\s+$//g; $v =~ s/\s.*$//; return $v; }
sub fasta_lengths {
	my ($path)=@_; open my $f, '<', $path or die "Cannot open '$path': $!\n";
	my (%length, @order); my $id='';
	while (my $line=<$f>) { $line =~ s/[\r\n]+$//;
		if ($line =~ /^>(.*)$/) { $id=canonical_id($1); die "Duplicate/empty FASTA id '$id'\n" if $id eq '' || exists $length{$id}; $length{$id}=0; push @order,$id; next; }
		die "Sequence before FASTA header in '$path'\n" if $id eq '' && $line ne '';
		$length{$id} += length($line) if $id ne '';
	}
	close $f; die "No FASTA records in '$path'\n" unless @order; return (\%length, \@order);
}
sub complete_runs {
	my ($id,$length,$intervals)=@_; my @sorted=sort {$a->[0]<=>$b->[0] || $a->[1]<=>$b->[1]} @$intervals;
	my (@runs,$cursor); $cursor=0;
	for my $r (@sorted) { die "Overlapping coverage for '$id' at $r->[0]-$r->[1]\n" if $r->[0] < $cursor;
		push @runs, [$cursor,$r->[0],0] if $r->[0]>$cursor; push @runs, [@$r]; $cursor=$r->[1]; }
	push @runs, [$cursor,$length,0] if $cursor<$length; return \@runs;
}
sub mean_depth {
	my ($runs,$start,$end)=@_; return 0 if $end <= $start; my $sum=0;
	for my $r (@$runs) { last if $r->[0] >= $end; next if $r->[1] <= $start;
		my $s=$r->[0]>$start?$r->[0]:$start; my $e=$r->[1]<$end?$r->[1]:$end; $sum += ($e-$s)*$r->[2]; }
	return $sum/($end-$start);
}
