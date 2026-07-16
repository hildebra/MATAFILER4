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

open my $map_fh, '<', File::Spec->catfile($Bin, '..', 'examples', 'maps', 'testAG.map')
	or die "Cannot read testAG.map: $!";
my $map = do { local $/; <$map_fh> };
close $map_fh;
like($map, qr/^#OutPath\s+\$MF4DIR\/examples\/output\//m,
	'Illumina example output uses the documented MF4DIR variable');
like($map, qr/^#DirPath\s+\$MF4DIR\/examples\/\/data/m,
	'Illumina example input uses the documented MF4DIR variable');
unlike($map, qr/\$MGTKDIR/, 'Illumina example no longer uses the obsolete variable');

done_testing;
