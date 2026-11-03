use strict;
use warnings;

use File::Spec;
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::GenoMetaAss qw(resolve_path);

{
	local $ENV{MF4DIR} = '/opt/matafiler4';
	is(resolve_path('$MF4DIR/examples//output/'), '/opt/matafiler4/examples/output/',
		'expands a configured environment variable and normalizes separators');
	is(resolve_path('${MF4DIR}/examples/data'), '/opt/matafiler4/examples/data',
		'supports braced environment-variable paths');
}

{
	local $ENV{MATAFILER_TEST_UNSET};
	delete $ENV{MATAFILER_TEST_UNSET};
	my $error = '';
	eval { resolve_path('$MATAFILER_TEST_UNSET/examples/output') };
	$error = $@;
	like($error, qr/Environment variable \$MATAFILER_TEST_UNSET .* is not set/,
		'unset variables produce an actionable error instead of a root-level path');
}

done_testing;
