use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::WorkflowResilience qw(
	retry_operation retry_unlink retry_rename
	preflight_executable preflight_directory filesystem_capacity
);

my ($attempts, @sleeps);
my $result = retry_operation(
	label => 'injected transient operation',
	attempts => 3,
	delays => [0, 0],
	sleeper => sub { push @sleeps, $_[0] },
	code => sub { ++$attempts >= 3 ? 'ready' : 0 },
);
is($result, 'ready', 'retry_operation returns the eventual successful result');
is($attempts, 3, 'retry_operation retries only until success');
is_deeply(\@sleeps, [0, 0], 'retry_operation uses the configured delays');

my $warning = '';
my $nonfatal;
{
	local $SIG{__WARN__} = sub { $warning .= shift };
	$nonfatal = retry_operation(
		label => 'injected permanent operation',
		attempts => 2,
		delays => [0],
		sleeper => sub { },
		fatal => 0,
		code => sub { 0 },
	);
}
ok(!defined($nonfatal), 'nonfatal exhausted retries return undef');
like($warning, qr/failed after 2 attempt/, 'nonfatal exhausted retries remain visible');

my $directory = tempdir(CLEANUP => 1);
ok(preflight_directory($directory, 'test workflow directory'),
	'preflight_directory validates create/write/rename/cleanup behavior');
my $source = File::Spec->catfile($directory, 'source');
my $destination = File::Spec->catfile($directory, 'destination');
open my $fh, '>', $source or die "Cannot create $source: $!";
print {$fh} "payload\n";
close $fh or die "Cannot close $source: $!";
ok(retry_rename($source, $destination), 'retry_rename publishes a file');
ok(-s $destination, 'retry_rename leaves a nonempty destination');
ok(retry_unlink($destination), 'retry_unlink removes the published file');
ok(!-e $destination, 'retry_unlink verifies the file is gone');

is(preflight_executable($^X, 'Perl'), $^X,
	'preflight_executable accepts a configured absolute executable');
my $capacity = filesystem_capacity($directory);
ok(ref($capacity) eq 'HASH', 'filesystem_capacity returns structured capacity state');
ok(exists($capacity->{available_kb}) && exists($capacity->{available_inodes}),
	'filesystem capacity reports blocks and inodes when df is available');

done_testing();
