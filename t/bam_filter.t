use strict;
use warnings;

use File::Spec;
use FindBin qw($Bin);
use IPC::Open3;
use Symbol qw(gensym);
use Test::More;

my $script = File::Spec->catfile($Bin, '..', 'secScripts', 'assemblies', 'bamFilter.pl');

sub sam_record {
	my (%args) = @_;
	my @fields = (
		$args{name}  // 'read',
		$args{flag}  // 0,
		$args{rname} // 'ref',
		$args{pos}   // 1,
		$args{mapq}  // 60,
		$args{cigar} // '10M',
		'*', 0, 0,
		$args{seq}  // ('A' x 10),
		$args{qual} // ('I' x 10),
	);
	push @fields, "NM:i:$args{nm}" if exists $args{nm};
	return join("\t", @fields) . "\n";
}

sub run_filter {
	my ($input, @arguments) = @_;
	my $error = gensym;
	my $pid = open3(my $stdin, my $stdout, $error, $^X, $script, @arguments);
	print {$stdin} $input;
	close $stdin;
	my $output = do { local $/; <$stdout> };
	my $errors = do { local $/; <$error> };
	waitpid($pid, 0);
	return ($? >> 8, $output // '', $errors // '');
}

my $header = "\@SQ\tSN:ref\tLN:1000\n";
my $passing = sam_record(name => 'pass', nm => 0);
my $low_mapq = sam_record(name => 'low-mapq', mapq => 5, nm => 0);
my $unmapped = sam_record(
	name => 'unmapped', flag => 4, rname => '*', pos => 0, mapq => 0,
	cigar => '*', seq => '*', qual => '*',
);
my ($status, $output, $errors) = run_filter(
	$header . $passing . $low_mapq . $unmapped, 0.05, 0.8, 20, 0,
);
is($status, 0, 'valid SAM stream is filtered successfully');
my @output_lines = split /\n/, $output;
is($output_lines[0], "\@SQ\tSN:ref\tLN:1000", 'SAM header is preserved');
is((split /\t/, $output_lines[1])[1], 0, 'passing mapped record remains mapped');
is((split /\t/, $output_lines[2])[1], 4, 'failed record is actually marked unmapped');
is($output_lines[3], substr($unmapped, 0, -1), 'already-unmapped record passes through unchanged');
like($errors, qr/Newly filtered records: 1/, 'summary counts newly filtered records');
like($errors, qr/Failure reasons \(non-exclusive\)/, 'summary labels overlapping reason counts');

my $no_quality = sam_record(name => 'no-quality', qual => '*', nm => 0);
($status, $output, $errors) = run_filter($no_quality);
is($status, 0, 'QUAL=* is accepted for a mapped record');
is((split /\t/, $output)[1], 0, 'record without qualities can pass filtering');

my $hard_clipped = sam_record(
	name => 'hard-clipped', cigar => '10H90M', seq => ('A' x 90),
	qual => ('I' x 90), nm => 0,
);
($status, $output, $errors) = run_filter($hard_clipped, 0.05, 0.9, 20, 0);
is($status, 0, 'hard-clipped record is evaluated successfully');
is((split /\t/, $output)[1], 0,
	'hard clips extend original query length instead of being subtracted from SEQ');

my $insertion = sam_record(
	name => 'insertion', cigar => '10S85M5I', seq => ('A' x 100),
	qual => ('I' x 100), nm => 5,
);
($status, $output, $errors) = run_filter($insertion, 0.06, 0.9, 20, 0);
is($status, 0, 'insertion-containing record is evaluated successfully');
is((split /\t/, $output)[1], 0, 'insertions count as aligned query coverage');

my $exact_end_clip = sam_record(
	name => 'end-clipped', cigar => '3S94M3S', seq => ('A' x 100),
	qual => ('I' x 100), nm => 0,
);
($status, $output, $errors) = run_filter($exact_end_clip, 0.05, 0.8, 20, 3);
is($status, 0, 'end-clipped record is evaluated successfully');
is((split /\t/, $output)[1], 4, 'end clipping threshold is inclusive');

my $combined_end_clip = sam_record(
	name => 'combined-end-clip', cigar => '1H2S94M2S1H', seq => ('A' x 98),
	qual => ('I' x 98), nm => 0,
);
($status, $output, $errors) = run_filter($combined_end_clip, 0.05, 0.8, 20, 3);
is($status, 0, 'combined hard/soft-clipped record is evaluated successfully');
is((split /\t/, $output)[1], 4, 'adjacent hard and soft clipping is combined at each end');

my $deletion = sam_record(
	name => 'deletion', cigar => '90M10D', seq => ('A' x 90),
	qual => ('I' x 90), nm => 10,
);
($status, $output, $errors) = run_filter($deletion, 0.05, 0.8, 20, 0);
is($status, 0, 'deletion-containing record is evaluated successfully');
is((split /\t/, $output)[1], 4, 'deletions contribute to edit rate using alignment columns');

my $missing_nm = sam_record(name => 'missing-nm');
($status, $output, $errors) = run_filter($missing_nm);
isnt($status, 0, 'missing NM tag fails closed');
like($errors, qr/missing a valid NM:i tag/, 'missing NM failure is actionable');

($status, $output, $errors) = run_filter("broken\tSAM\n");
isnt($status, 0, 'malformed SAM returns a failure status');
like($errors, qr/expected at least 11 fields/, 'malformed SAM reports its field count problem');

my $bad_cigar_length = sam_record(name => 'bad-length', cigar => '9M', nm => 0);
($status, $output, $errors) = run_filter($bad_cigar_length);
isnt($status, 0, 'CIGAR and SEQ length mismatch returns a failure status');
like($errors, qr/CIGAR consumes 9 query bases/, 'CIGAR mismatch reports both data sources');

($status, $output, $errors) = run_filter('', 1.1);
isnt($status, 0, 'invalid threshold returns a failure status');
like($errors, qr/max_edit_rate must be a number between 0 and 1/, 'invalid argument has clear usage guidance');

done_testing();
