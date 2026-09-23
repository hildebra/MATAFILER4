use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

# The shipped sdm option files, after MATAF4's read-length adaptation, must not
# discard full-length reads of the read lengths each technology commonly
# produces. The option helpers are loaded from MATAF4.pl so the test exercises
# the real adaptation code.

my $root = File::Spec->catdir($Bin, '..');
open my $sourceFH, '<', File::Spec->catfile($root, 'MATAF4.pl')
	or die "Cannot inspect MATAF4.pl: $!";
my $mataf4 = do { local $/; <$sourceFH> };
close $sourceFH;
my ($optionHelpers) = $mataf4 =~ /(sub _shell_quote\s*\{.*?)(?=^sub sdmOptSet)/ms;
BAIL_OUT('Cannot isolate the SDM option helpers from MATAF4.pl') unless defined $optionHelpers;
eval "package TestSDMOptions; our \%MFopt; $optionHelpers; 1" or BAIL_OUT("SDM helpers: $@");
$TestSDMOptions::MFopt{sdmProbabilisticFilter} = 1;
$TestSDMOptions::MFopt{sdm_opt} = {};

sub options_of {
	my ($path) = @_;
	open my $fh, '<', $path or die "Cannot read $path: $!";
	my %options;
	while (my $line = <$fh>) {
		next if $line =~ /^\s*(?:#|\*|$)/;
		chomp $line;
		my ($key, $value) = split /\t/, $line;
		$options{$key} = $value if defined $value;
	}
	return \%options;
}

my $data = File::Spec->catdir($root, 'data');
opendir my $dh, $data or die "Cannot list $data: $!";
my @files = sort grep { /^sdm_.*\.txt$/ } readdir $dh;
closedir $dh;
ok(@files >= 8, 'shipped sdm option files were found');
for my $file (@files) {
	my $options = options_of(File::Spec->catfile($data, $file));
	ok(exists $options->{minSeqLength}, "$file defines minSeqLength");
	like($options->{minSeqLength} // '', qr/^\d+(?:\.\d+)?$/, "$file minSeqLength is numeric");
}

# technology => [option file, read lengths that are routinely sequenced]
my %technology = (
	hiSeq => ['sdm_opt_inifilter.txt', [100, 150]],
	ill   => ['sdm_opt_inifilter_relaxed.txt', [100, 150]],
	miSeq => ['sdm_opt_miSeq.txt', [150, 250, 300]],
	AVITI => ['sdm_AVITI.txt', [150, 300]],
);
my %known_issue = (
	'miSeq:150' => 1,
	'AVITI:150' => 1,
);

my $tmp = tempdir(CLEANUP => 1);
for my $tech (sort keys %technology) {
	my ($file, $lengths) = @{ $technology{$tech} };
	for my $length (@{$lengths}) {
		my $adapted = TestSDMOptions::adaptSDMopt(File::Spec->catfile($data, $file), $tmp, $length, $tech);
		my $options = options_of($adapted);
		my $minimum = $options->{minSeqLength};
		my $effective = $minimum <= 1 ? $minimum * $length : $minimum;
		{
			local $TODO = $known_issue{"$tech:$length"}
				? "absolute minSeqLength $minimum in $file is not rescaled for $length bp reads"
				: undef;
			ok($effective <= $length,
				"$tech $length bp: full-length reads pass minSeqLength ($minimum -> $effective bp)");
		}
		my $truncate = $options->{TruncateSequenceLength};
		if (defined $truncate && $truncate > 0) {
			ok($truncate >= $effective,
				"$tech $length bp: truncation length $truncate keeps reads above minSeqLength");
		}
	}
}

done_testing();
