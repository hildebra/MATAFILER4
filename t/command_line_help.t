use strict;
use warnings;

use Cwd qw(getcwd);
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

my $root = "$Bin/..";
my $script = "$root/MATAF4.pl";
my $reference = "$root/docs/flag_reference.md";

open(my $script_fh, '<', $script) or die "Cannot read $script: $!";
my $source = do { local $/; <$script_fh> };
close($script_fh);

my ($get_options) = $source =~
	/sub getCmdLineOptions\s*\{(.*?)\n\s*\)(?:\s+or\s+die\s+[^;]+)?;/s;
ok(defined($get_options), 'found the MATAF4.pl GetOptions block');

my %accepted;
while (defined($get_options) && $get_options =~ /^\s*"([^"]+)"\s*=>/mg) {
	my $spec = $1;
	$spec =~ s/=[sif]$//;
	$accepted{$_} = 1 for split(/\|/, $spec);
}

open(my $reference_fh, '<', $reference) or die "Cannot read $reference: $!";
my %documented;
my $in_mataf4 = 0;
while (my $line = <$reference_fh>) {
	if ($line =~ /^##\s+MATAF4\.pl\s*$/i) {
		$in_mataf4 = 1;
		next;
	}
	last if ($in_mataf4
		&& $line =~ /^##\s+Flag comparison against previous manual\.md\s*$/i);
	next unless ($in_mataf4 && $line =~ /^\|\s*(`-[^|]+?)\s*\|/);
	my $aliases = $1;
	while ($aliases =~ /`-([A-Za-z0-9_?]+)`/g) {
		$documented{$1} = 1;
	}
}
close($reference_fh);

open($reference_fh, '<', $reference) or die "Cannot read $reference: $!";
my $reference_text = do { local $/; <$reference_fh> };
close($reference_fh);
unlike(
	$reference_text,
	qr/Accepted by [A-Za-z0-9_.]+; (?:inspect|see) source\/help/i,
	'flag descriptions omit repetitive acceptance boilerplate',
);
like(
	$reference_text,
	qr/`-reProfileRibosome`.*?Remove existing RiboFind extraction and assignment results/s,
	'flag reference describes RiboFind profile invalidation',
);

is_deeply(
	[sort keys %documented],
	[sort keys %accepted],
	'docs/flag_reference.md documents every accepted MATAF4.pl option and alias',
);

my $outside = tempdir(CLEANUP => 1);
my $original_dir = getcwd();
chdir($outside) or die "Cannot enter $outside: $!";
my ($help, $help_status);
{
	local $ENV{PERL5LIB};
	delete $ENV{PERL5LIB};
	$help = qx{"$^X" "$script" --help 2>&1};
	$help_status = $? >> 8;
}
chdir($original_dir) or die "Cannot return to $original_dir: $!";
is($help_status, 0, '--help loads checkout-local modules and exits successfully outside the repository');
like($help, qr/Usage:\s+MATAF4\.pl -map <mapping-file> \[options\]/,
	'help shows command syntax');
like($help, qr/read from docs\/flag_reference\.md/,
	'help identifies the flag reference as its source');

my @missing_from_help = grep {
	$help !~ /(?<![A-Za-z0-9_])-\Q$_\E(?![A-Za-z0-9_])/
} sort keys %accepted;
is_deeply(\@missing_from_help, [], '--help displays every accepted option and alias');

done_testing();
