use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, "..");

# geneCat.pl -binSpeciesMG -1 (the default) bases MGS.pl on the binner that
# MATAF4.pl ran in the assembly phase, found in <assembly>/Binning/<binner>/.

my $root = File::Spec->rel2abs(File::Spec->catdir($Bin, '..'));
my $script = File::Spec->catfile($root, 'secScripts', 'geneCat.pl');
open my $fh, '<', $script or die "Cannot read $script: $!";
my $source = do { local $/; <$fh> };
close $fh;
my ($detect) = $source =~ /^(sub _detect_assembly_binner\s*\{.*?^\})/ms;
BAIL_OUT('Cannot isolate _detect_assembly_binner from geneCat.pl') unless defined $detect;
like($source, qr/^my \$binSpeciesMG = -1;/m, '-binSpeciesMG defaults to auto-detection');
like($source, qr/"\$OutD\/Bin_\$\{binnerShrt\}\/"/, 'the MGS output folder follows the selected binner');

{
	package GeneCatBinnerTest;
	use Mods::Binning qw(getBinSubdirName);
	use Mods::GenoMetaAss qw(getAssemblPath);
	eval "$detect; 1" or die $@;
}

sub touch {
	my ($path) = @_;
	open my $out, '>', $path or die "Cannot write $path: $!";
	close $out or die "Cannot close $path: $!";
}

my $tmp = tempdir(CLEANUP => 1);
$tmp =~ s{\\}{/}g;
my %map; my @samples;
sub add_sample {
	my ($smpl, $binner, %opt) = @_;
	my $dir = "$tmp/$smpl";
	my $metaGD = "$dir/assemblies/metag";
	make_path("$metaGD/Binning/$binner");
	open my $out, '>', "$metaGD/assembly.txt" or die $!;
	print {$out} "$metaGD\n";
	close $out or die $!;
	touch("$metaGD/Binning/$binner/$smpl") unless $opt{noOutput};
	$map{$smpl} = { wrdir => $dir, prefix => '', ExcludeAssem => $opt{exclude} ? '1' : '0' };
	push @samples, $smpl;
}

sub detect {
	my @warnings;
	local $SIG{__WARN__} = sub { push @warnings, $_[0] };
	my $stdout = '';
	open my $capture, '>', \$stdout or die $!;
	my $old = select $capture;
	my $binner = GeneCatBinnerTest::_detect_assembly_binner(\%map, \@samples);
	select $old;
	return ($binner, join('', @warnings), $stdout);
}

add_sample('S1', 'SC');
add_sample('S2', 'SC');
my ($binner, $warn, $out) = detect();
is($binner, 5, 'SCGBinner outputs select -binSpeciesMG 5');
like($out, qr/Bin_SC\//, 'the run log names the Bin_SC output folder');
like($out, qr/2\/2 assembly groups/, 'the count is reported per assembly group');
is($warn, '', 'a consistent binner gives no warning');

add_sample('S3', 'SB');
($binner, $warn) = detect();
is($binner, 5, 'the binner found in most assembly groups wins');
like($warn, qr/several binners.*SC=2, SB=1/, 'mixed binner outputs are reported');
like($warn, qr/1 of 3 assembly groups have no SC/, 'assembly groups without the chosen binner are reported');

add_sample('S4', 'SB', exclude => 1);
add_sample('S5', 'SB', exclude => 1);
($binner) = detect();
is($binner, 5, 'excluded assemblies are not counted');

%map = (); @samples = ();
add_sample('E1', 'MB2', noOutput => 1);
($binner, $warn) = detect();
is($binner, 2, 'without binner outputs the SemiBin default is kept');
like($warn, qr/no binner output found/, 'missing binner outputs are reported');

done_testing();
