use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

# The R-analysis batches are generated shell, so the properties that matter are
# properties of the emitted script: one failing MGS must not stop the batch, a
# batch that lost an MGS must still report failure, and a retry must skip the
# MGS that already produced a store. Extract the real generators and run them.
my $root = File::Spec->catdir($Bin, '..');
my $source = File::Spec->catfile($root, 'secScripts', 'MGS', 'strain_within_2.2.pl');
open my $handle, '<', $source or die "Cannot read $source: $!";
my $text = do { local $/; <$handle> };
close $handle or die "Cannot close $source: $!";

my ($prelude) = $text =~ /^(my \$cmdPrelude = .*?;\n)/ms;
my ($epilogue) = $text =~ /^(my \$cmdEpilogue = .*?;\n)/ms;
my ($isolator) = $text =~ /^(my \$isolatedAnalysisBlock = sub \{.*?^\};\n)/ms;
my ($quoter) = $text =~ /^(sub shellQuote \{.*?^\})/ms;
BAIL_OUT('Cannot extract the R-analysis batch generators')
	unless defined($prelude) && defined($epilogue)
		&& defined($isolator) && defined($quoter);

my ($cmdPrelude, $cmdEpilogue, $isolatedAnalysisBlock) = eval
	"$quoter\n$prelude\n$epilogue\n$isolator\n"
	. '($cmdPrelude, $cmdEpilogue, $isolatedAnalysisBlock);';
BAIL_OUT("Cannot load the R-analysis batch generators: $@") if $@;
ok(defined($cmdPrelude) && defined($cmdEpilogue)
	&& ref($isolatedAnalysisBlock) eq 'CODE',
	'batch generators load standalone');

unlike($cmdPrelude, qr/set -\w*u/,
	'the prelude does not enable nounset, which conda activation hooks do not survive');
like($cmdPrelude, qr/set -e/, 'the prelude still stops on unexpected errors');
like($cmdPrelude, qr/ulimit[^\n]*\|\| true/,
	'a stack limit the node refuses does not fail the job');

my $bash = -x '/bin/bash' ? '/bin/bash' : 'bash';
my $temporary = tempdir('strain-batch-XXXXXX', TMPDIR => 1, CLEANUP => 1);

# Two analyses in one batch: the first cannot produce its store, the second can.
sub build_batch {
	my ($firstSucceeds) = @_;
	my $cmd = $cmdPrelude;
	for my $case (['A', $firstSucceeds], ['B', 1]) {
		my ($name, $succeeds) = @{$case};
		my $store = File::Spec->catfile($temporary, "store$name");
		my $log = File::Spec->catfile($temporary, "log$name");
		my $body = $succeeds
			? "printf 'result\\n' > ".shellQuote($store)." 2> ".shellQuote($log)."\n"
			: "false > ".shellQuote($log)."\n";
		$cmd .= "echo ".shellQuote("At tree MGS.$name")."\n";
		$cmd .= "if ! test -s ".shellQuote($store)."; then\n";
		# The real generator announces the start inside the completion guard, so
		# a skipped analysis prints nothing between its "At tree" line and the next.
		$cmd .= "echo ".shellQuote("Starting strainStats for MGS.$name")."\n";
		$cmd .= $isolatedAnalysisBlock->($body, $store, "strainStats for MGS.$name", $log);
		$cmd .= "fi\n";
	}
	$cmd .= "echo ".shellQuote('BATCH END')."\n";
	return $cmd . $cmdEpilogue;
}

sub run_batch {
	my ($cmd) = @_;
	my $script = File::Spec->catfile($temporary, 'batch.sh');
	open my $out, '>', $script or die "Cannot write $script: $!";
	print {$out} $cmd or die "Cannot populate $script: $!";
	close $out or die "Cannot close $script: $!";
	my $output = `"$bash" "$script" 2>&1`;
	return ($? >> 8, $output);
}

unlink glob(File::Spec->catfile($temporary, 'store*'));
my ($status, $output) = run_batch(build_batch(0));
like($output, qr/FAILED strainStats for MGS\.A/,
	'a failing analysis is reported by name');
like($output, qr/Completed strainStats for MGS\.B/,
	'the next analysis in the batch still runs after a failure');
like($output, qr/BATCH END/, 'the batch runs to completion instead of aborting');
isnt($status, 0,
	'a batch that lost an analysis still exits non-zero for Slurm accounting');
ok(-s File::Spec->catfile($temporary, 'storeB'),
	'the surviving analysis produced its durable store');
ok(!-e File::Spec->catfile($temporary, 'storeA'),
	'the failing analysis produced no store, so the parent will retry its batch');

# A retry reruns the same command: the completed analysis must be skipped.
($status, $output) = run_batch(build_batch(0));
like($output, qr/Starting strainStats for MGS\.A/,
	'a retry reattempts the analysis that failed');
unlike($output, qr/Starting strainStats for MGS\.B/,
	'a retry skips the analysis whose store already exists, so retries make progress');
isnt($status, 0, 'a retry that still cannot finish keeps reporting failure');

unlink glob(File::Spec->catfile($temporary, 'store*'));
($status, $output) = run_batch(build_batch(1));
is($status, 0, 'a batch whose analyses all succeed exits zero');
like($output, qr/Completed strainStats for MGS\.A/,
	'the recovered analysis is reported as complete');
unlike($output, qr/FAILED/, 'a fully successful batch reports no failures');

done_testing();
