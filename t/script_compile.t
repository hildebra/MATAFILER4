use strict;
use warnings;

use File::Find ();
use File::Spec;
use FindBin qw($Bin);
use Test::More;

# Every Perl entry point and module must at least compile. Scripts that are
# only reached through generated job scripts are otherwise never loaded by the
# test suite, so a syntax error or a broken import would first surface on the
# cluster. Failures caused only by an optional CPAN module that is not
# installed in the current environment are reported as skips, not passes.

my $root = File::Spec->rel2abs(File::Spec->catdir($Bin, '..'));
my @files;
File::Find::find({ no_chdir => 1, wanted => sub {
	my $path = $File::Find::name;
	if (-d $path && $path =~ m{/(?:\.git|t)$}) { $File::Find::prune = 1; return; }
	return unless -f $path && $path =~ /\.(?:pl|pm)$/;
	push @files, $path;
} }, $root);
@files = sort @files;
ok(@files > 100, 'Perl sources were discovered (' . scalar(@files) . ')');

local $ENV{PERL5LIB} = join(':', $root, File::Spec->catdir($root, 't', 'lib'),
	grep { defined && length } $ENV{PERL5LIB});

for my $file (@files) {
	my $name = File::Spec->abs2rel($file, $root);
	my $output = qx{"$^X" -c "$file" 2>&1};
	my $status = $?;
	if ($status == 0) {
		pass("$name compiles");
		next;
	}
	my @missing = $output =~ /Can't locate (\S+\.pm) in \@INC/g;
	my @external = grep { !m{^Mods/} } @missing;
	SKIP: {
		skip("$name needs optional module(s) not installed here: @external", 1)
			if @external && @external == @missing;
		fail("$name compiles");
		diag($output);
	}
}

done_testing();
