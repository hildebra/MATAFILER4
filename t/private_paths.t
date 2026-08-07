use strict;
use warnings;

use File::Find;
use File::Spec;
use FindBin qw($Bin);
use Test::More;

my $root = File::Spec->catdir($Bin, '..');
my @scan_roots = map { File::Spec->catdir($root, $_) }
	qw(MATAF4.pl Mods helpers secScripts docs examples);
my $private_path = qr{
	(?:
		/hpc-home/ |
		/ei/(?:projects|workarea/users)/ |
		/qib/ |
		/nbi/ |
		/g/(?:bork\d*|scb|korbel)/ |
		/(?:local|scratch)/bork/ |
		/tmp/(?:hildebra|rel25vaz)(?:/|\b) |
		\\Users\\(?:hildebra|rel25vaz)\\
	)
}x;
my @violations;

find(
	{
		no_chdir => 1,
		wanted => sub {
			return unless -f $_;
			return if $_ =~ m{/\.[^/]*\.swp$};
			return if -B $_; # prebuilt/generated binaries are not maintained source
			open my $fh, '<:raw', $_ or die "Cannot read $_: $!\n";
			my $line_number = 0;
			while (my $line = <$fh>){
				$line_number++;
				# Historic examples may remain in comments. Only executable/config
				# content is guarded, as requested.
				$line =~ s/#.*$//;
				push @violations, "$_:$line_number"
					if $line =~ $private_path;
			}
			close $fh or die "Cannot close $_: $!\n";
		},
	},
	@scan_roots,
);

is_deeply(\@violations, [],
	'active source, config, and documentation examples contain no private paths')
	or diag join("\n", @violations);

done_testing;
