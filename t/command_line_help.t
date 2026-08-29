use strict;
use warnings;

use Cwd qw(getcwd);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::FlagReference qw(readFlagReference scriptOptionNames helpRequested);

my $root = "$Bin/..";
my $reference = "$root/docs/flag_reference.md";

#Every entry point renders -help from docs/flag_reference.md via
#Mods::FlagReference, so the option tables exist exactly once. These checks keep
#the three views in step: the GetOptions block, the reference section, and the
#help text a user actually sees.
my @entryPoints = (
	{ name => 'MATAF4.pl',            path => "$root/MATAF4.pl",
		options_from => 'sub getCmdLineOptions' },
	{ name => 'geneCat.pl',           path => "$root/secScripts/geneCat.pl" },
	{ name => 'MGS.pl',               path => "$root/secScripts/MGS.pl" },
	{ name => 'strain_within.pl',     path => "$root/secScripts/MGS/strain_within.pl" },
	{ name => 'strain_within_2.2.pl', path => "$root/secScripts/MGS/strain_within_2.2.pl" },
	{ name => 'buildTree5.pl',        path => "$root/secScripts/phylo/buildTree5.pl" },
);

#Option names documented in this script's section of the reference.
sub documentedOptions {
	my ($script) = @_;
	my %documented;
	for my $section (readFlagReference($script, $reference)) {
		for my $entry (@{$section->{options}}, @{$section->{legacy}}) {
			while ($entry->{aliases} =~ /-([A-Za-z0-9_?]+)/g) {
				$documented{$1} = 1;
			}
		}
	}
	return %documented;
}

for my $entry (@entryPoints) {
	my $name = $entry->{name};
	my %accepted = scriptOptionNames($entry->{path});
	ok(scalar(keys %accepted), "found the $name GetOptions block");

	my %documented = documentedOptions($name);
	is_deeply([sort keys %documented], [sort keys %accepted],
		"docs/flag_reference.md documents every $name option and alias");
}

#The shared help-request test must accept both dash spellings and nothing else.
ok(helpRequested('-help'),  'helpRequested accepts -help');
ok(helpRequested('--help'), 'helpRequested accepts --help');
ok(helpRequested('-h'),     'helpRequested accepts -h');
ok(helpRequested('-?'),     'helpRequested accepts -?');
ok(!helpRequested('-GCd', '/tmp/x', '-onlyMSA', 1),
	'helpRequested ignores ordinary options');
ok(!helpRequested('-helper'), 'helpRequested does not match option prefixes');

#The reference must not fall back to placeholders instead of real descriptions.
open(my $reference_fh, '<', $reference) or die "Cannot read $reference: $!";
my $reference_text = do { local $/; <$reference_fh> };
close($reference_fh);
unlike($reference_text,
	qr/Accepted by [A-Za-z0-9_.]+; (?:inspect|see) source\/help/i,
	'flag descriptions omit repetitive acceptance boilerplate');
unlike($reference_text, qr/See source\/help for details/i,
	'flag descriptions are not placeholders');
like($reference_text,
	qr/`-reProfileRibosome`.*?delete RiboFind extraction, assignments and merged profiles/s,
	'flag reference describes RiboFind profile invalidation');

#Each script must actually render its own section, with every option in it.
my $outside = tempdir(CLEANUP => 1);
my $original_dir = getcwd();
for my $entry (@entryPoints) {
	my $name = $entry->{name};
	my %accepted = scriptOptionNames($entry->{path});

	chdir($outside) or die "Cannot enter $outside: $!";
	my ($help, $status);
	{
		#the scripts are always invoked with the checkout on the module path
		local $ENV{PERL5LIB} = File::Spec->rel2abs($root, $original_dir);
		$help = qx{"$^X" "$entry->{path}" -help 2>&1};
		$status = $? >> 8;
	}
	chdir($original_dir) or die "Cannot return to $original_dir: $!";

	is($status, 0, "$name -help exits successfully outside the repository");
	like($help, qr/\Q$name\E command-line help/, "$name -help names itself");
	like($help, qr/read from docs\/flag_reference\.md/,
		"$name -help identifies the flag reference as its source");
	like($help, qr/Usage:\s+\Q$name\E /, "$name -help shows command syntax");

	my @missing = grep {
		$help !~ /(?<![A-Za-z0-9_])-\Q$_\E(?![A-Za-z0-9_])/
	} sort keys %accepted;
	is_deeply(\@missing, [], "$name -help displays every accepted option and alias");
}

done_testing();
