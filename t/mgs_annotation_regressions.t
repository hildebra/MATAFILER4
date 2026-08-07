use strict;
use warnings;

use File::Spec;
use FindBin qw($Bin);
use Test::More;

my $script = File::Spec->catfile(
	$Bin, '..', 'secScripts', 'GC', 'annotateMGwSpecIs2.pl',
);
open my $fh, '<', $script or die "Cannot read $script: $!";
my $source = do { local $/; <$fh> };
close $fh or die "Cannot close $script: $!";

like($source, qr/my \$version = 0\.14;/,
	'MGS/specI annotation records the shared-marker behavior change');
like($source, qr/next if \$spl\[0\] eq "Bin";.*?Malformed MGS row/s,
	'MGS input headers are skipped and malformed rows fail explicitly');
like($source,
	qr/\$markerOwners\{\$gen\}\{\$curMGS\} = 1.*?my \@owners = sort keys.*?if \(\@owners != 1\).*?next;/s,
	'markers are assigned only after unique MGS ownership is established');
unlike($source, qr/\$Gene2MGS\{\$gen\} = \$curMGS;\s+\$MGSlist\{\$curMGS\}/,
	'MGS ownership is not overwritten in input order');
like($source, qr/foreach my \$MGS \(sort keys %MGSlist\).*?\$curS =0/s,
	'marker-count summaries retain MGS with zero unambiguous markers');
like($source,
	qr/if \(exists\(\$tmpTax\{\$curMGS\}\)\)\{\s+\$specIfullTax\{\$curMGS\} = \$tmpTax\{\$curMGS\}/s,
	'explicit MGS taxonomy takes precedence over any colliding reference identifier');

done_testing();
