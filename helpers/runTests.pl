#!/usr/bin/env perl
# Run the MATAFILER4 unit-test suite (t/*.t) and summarise the results.
#
#   perl helpers/runTests.pl                 # run everything, 4 parallel jobs
#   perl helpers/runTests.pl -jobs 1 -v      # sequential, verbose TAP output
#   perl helpers/runTests.pl strain mgs      # only tests whose name matches
#   perl helpers/runTests.pl -list           # list available test files
#   perl helpers/runTests.pl -report out.tsv # also write a per-file summary
#
# The runner can be started from any directory. It puts the repository and
# t/lib on @INC, restores the executable bit of the test shims in t/bin
# (lost e.g. on Windows checkouts), and reports which compiled helpers from
# bin/ are unavailable, since tests that need them skip themselves.
# Exit status: 0 if all tests pass, 1 if any fail, 2 on usage errors.

use strict;
use warnings;

use Cwd qw(abs_path);
use File::Basename qw(dirname basename);
use File::Spec;
use FindBin qw($RealBin);
use Getopt::Long qw(GetOptions);
use TAP::Harness;
use Time::HiRes qw(time);

my $root = abs_path(File::Spec->catdir($RealBin, '..'));
my $test_dir = File::Spec->catdir($root, 't');

my %opt = (jobs => 4, verbose => 0, list => 0, report => '', timer => 0, help => 0);
GetOptions(
	'jobs|j=i'  => \$opt{jobs},
	'verbose|v' => \$opt{verbose},
	'list'      => \$opt{list},
	'report=s'  => \$opt{report},
	'timer'     => \$opt{timer},
	'help|h'    => \$opt{help},
) or usage(2);
usage(0) if $opt{help};
if ($opt{jobs} < 1) { print STDERR "-jobs must be at least 1\n"; exit 2; }

sub usage {
	my ($status) = @_;
	print <<"USAGE";
Usage: perl helpers/runTests.pl [options] [pattern ...]

Runs the MATAFILER4 unit tests in t/ and prints a summary.

Options:
  -jobs N      run N test files in parallel (default 4)
  -v           show the full TAP output of every test
  -timer       print the run time of each test file
  -list        list the test files (after pattern filtering) and exit
  -report F    write a tab-separated per-file summary to F
  -h           this help

Patterns are case-insensitive regular expressions matched against the test
file names, e.g. "strain" or "^mgs_".
USAGE
	exit $status;
}

opendir my $dh, $test_dir or die "Cannot list $test_dir: $!\n";
my @tests = sort grep { /\.t$/ && -f File::Spec->catfile($test_dir, $_) } readdir $dh;
closedir $dh;
if (@ARGV) {
	my @patterns = map { qr/$_/i } @ARGV;
	@tests = grep { my $name = $_; grep { $name =~ $_ } @patterns } @tests;
}
unless (@tests) { print STDERR "No test files match: @ARGV\n"; exit 2; }

if ($opt{list}) {
	print "$_\n" for @tests;
	exit 0;
}

# Tests locate fixtures relative to the repository; some use the working directory.
chdir $root or die "Cannot change to $root: $!\n";

# Test shims must be executable; a checkout without file modes loses this bit.
my $shim_dir = File::Spec->catdir($test_dir, 'bin');
if (opendir my $shims, $shim_dir) {
	for my $shim (grep { !/^\./ } readdir $shims) {
		my $path = File::Spec->catfile($shim_dir, $shim);
		next unless -f $path && !-x $path;
		my $mode = (stat $path)[2] & 07777;
		chmod($mode | 0100, $path) or warn "Cannot make test shim $path executable: $!\n";
	}
	closedir $shims;
}

# Compiled helpers some tests need; they skip themselves when these are absent.
my @compiled = qw(sdm rtk2 rdCover MSAfix vcf2fna LCA clusterMAGs cc.bin);
my @unavailable = grep { !-x File::Spec->catfile($root, 'bin', $_) } @compiled;

my @lib = ($root, File::Spec->catdir($test_dir, 'lib'));
$ENV{PERL5LIB} = join(':', @lib, grep { defined && length } $ENV{PERL5LIB});

printf "MATAFILER4 unit tests: %d file(s) from %s\n", scalar(@tests), $test_dir;
printf "Perl %s (%s)\n", sprintf('%vd', $^V), $^X;
print "Compiled helpers not available (dependent tests skip): @unavailable\n" if @unavailable;
print "\n";

my $started = time;
my $harness = TAP::Harness->new({
	lib       => \@lib,
	jobs      => $opt{jobs},
	verbosity => $opt{verbose} ? 1 : 0,
	timer     => $opt{timer},
	color     => (-t STDOUT ? 1 : 0),
	merge     => 0,
});
my $aggregate = $harness->runtests(map { File::Spec->catfile('t', $_) } @tests);
my $elapsed = time - $started;

# Per-file classification for the summary and optional report.
my (@failed, @skipped, @todo_passed, @rows);
for my $file (map { File::Spec->catfile('t', $_) } @tests) {
	my ($parser) = $aggregate->parsers($file);
	my $name = basename($file);
	unless ($parser) {
		push @failed, "$name (did not run)";
		push @rows, [$name, 'not-run', 0, 0, 0, ''];
		next;
	}
	my $status = 'pass';
	my $note = '';
	if ($parser->skip_all) {
		$status = 'skipped';
		$note = $parser->skip_all;
		push @skipped, "$name: $note";
	} elsif ($parser->has_problems) {
		$status = 'fail';
		my @bad = $parser->failed;
		$note = @bad ? 'failed test(s) ' . join(',', @bad) : 'exit status ' . ($parser->exit // '?');
		$note .= '; parse errors' if $parser->parse_errors;
		push @failed, "$name: $note";
	}
	if (my @unexpected = $parser->todo_passed) {
		push @todo_passed, "$name: TODO test(s) now pass: " . join(',', @unexpected);
	}
	push @rows, [$name, $status, scalar($parser->tests_run), scalar($parser->passed),
		scalar($parser->failed), $note];
}

print "\n", '=' x 72, "\n";
printf "Files: %d   Tests: %d   Passed: %d   Failed: %d   TODO: %d   Time: %.1fs\n",
	scalar(@tests), $aggregate->total, scalar($aggregate->passed),
	scalar($aggregate->failed), scalar($aggregate->todo), $elapsed;
if (@skipped) {
	print "\nSkipped test files:\n";
	print "  $_\n" for @skipped;
}
if (@todo_passed) {
	print "\nKnown issues that appear fixed (remove the TODO marker):\n";
	print "  $_\n" for @todo_passed;
}
if (@failed) {
	print "\nFAILED test files:\n";
	print "  $_\n" for @failed;
	print "\nRe-run one file verbosely with: perl helpers/runTests.pl -jobs 1 -v <name>\n";
}
print "\nResult: ", (@failed ? 'FAIL' : 'PASS'), "\n";

if ($opt{report}) {
	open my $out, '>', $opt{report} or die "Cannot write report $opt{report}: $!\n";
	print {$out} join("\t", qw(file status tests passed failed note)), "\n";
	print {$out} join("\t", @{$_}), "\n" for @rows;
	close $out or die "Cannot close report $opt{report}: $!\n";
	print "Per-file report written to $opt{report}\n";
}

exit(@failed ? 1 : 0);
