use strict;
use warnings;
use FindBin qw($Bin);
use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use JSON::PP qw(decode_json);
use Test::More;
use Mods::GenoMetaAss qw(systemW);
use Mods::Checkpoint qw(checkpoint_valid);

my $root = File::Spec->rel2abs(File::Spec->catdir($Bin, '..'));
my $source = read_file("$root/secScripts/geneCat.pl");
sub read_file {
	open my $fh, '<', $_[0] or die "Cannot read $_[0]: $!";
	local $/;
	return <$fh>;
}
sub write_file {
	open my $fh, '>', $_[0] or die "Cannot write $_[0]: $!";
	print {$fh} $_[1] or die "Cannot write $_[0]: $!";
	close $fh or die "Cannot close $_[0]: $!";
}
for my $name (qw(_shell_quote _checkpoint_command)) {
	my ($helper) = $source =~ /^(sub \Q$name\E \{.*?^\})/ms;
	die "Missing $name" unless defined $helper;
	eval $helper;
	die $@ if $@;
}
my ($matrix_helper) = $source =~ /^(sub _gene_matrix_commands \{.*?^\})/ms;
die 'Missing matrix command helper' unless defined $matrix_helper;
$matrix_helper =~ s/^sub _gene_matrix_commands/sub/;

my $tmp = tempdir(CLEANUP => 1);
my $fake_rare = "$tmp/fake rare.pl";
write_file($fake_rare, <<'FAKE');
use strict;
use warnings;
use JSON::PP qw(encode_json);
die 'configured command wrapper was lost' unless $ENV{GENECAT_TEST_WRAPPED};
my @args = @ARGV;
die 'wrong command' unless shift(@ARGV) eq 'geneMat';
my %options;
while (@ARGV) {
	my $option = shift @ARGV;
	$options{$option} = $option =~ /^-(?:i|t|map|refD|o)$/ ? shift(@ARGV) : 1;
}
my $output = $options{'-o'} // die 'missing output';
open my $matrix, '>', "$output.mat.gz" or die $!;
print {$matrix} encode_json(\@args) or die $!;
close $matrix or die $!;
open my $rows, '>', "$output.genes2rows.txt" or die $!;
print {$rows} "#Gene\n1\n" or die $!;
close $rows or die $!;
FAKE

for my $deferred (0, 1) {
	subtest($deferred ? 'combined deferred script' : 'separate matrix workers', sub {
		my $OutD = "$tmp/catalogue '$deferred with spaces";
		make_path($OutD);
		my $assDirs = "$tmp/assemblies with spaces";
		my $mapF = "$tmp/map 'one.txt,$tmp/map two.txt";
		my $primaryClusterCLS = 'compl.incompl.97.fna.clstr';
		my $countMatrixP = 'Matrix';
		my $rmBin = 'rm';
		my $countMatrixF = 'Matrix.mat';
		my $oldNameFolders = $deferred;
		my $CalcgGneMatSuppl = !$deferred;
		my $rareBin = 'env GENECAT_TEST_WRAPPED=1 ' . _shell_quote($^X) . ' ' . _shell_quote($fake_rare);
		my $builder = eval $matrix_helper;
		die $@ if $@;
		my @commands = $builder->($OutD, $assDirs, 4);
		is(scalar @commands, 3, 'builds all three abundance variants');

		my $cdhID = 97;
		my $checkpointWriter = _shell_quote($^X) . " " . _shell_quote("$root/helpers/writeCheckpoint.pl");
		my $matrixSton = "$OutD/matrix.stone";
		my ($done_source) = $source =~ /(my \$doneCommand = _checkpoint_command\(.*?\);)/s;
		die 'Missing matrix checkpoint command' unless defined $done_source;
		my $done = eval($done_source . "\n\$doneCommand");
		die $@ if $@;
		local $ENV{PERL5LIB} = join(':', $root, $ENV{PERL5LIB} // ());
		if ($deferred) {
			is(systemW(join("\n", @commands, $done), 0), 0,
				'all matrix commands and completion run successfully as one script');
		} else {
			is(systemW($_, 0), 0, 'matrix worker succeeds independently') for @commands;
			is(systemW($done, 0), 0, 'matrix convergence succeeds');
		}
		for my $variant (['Matrix', ''], ['Mat.cov', '-useCoverage'], ['Mat.med', '-useCovMedian']) {
			my ($prefix, $coverage_option) = @$variant;
			my $path = "$OutD/$prefix.mat.gz";
			ok(-s $path, "$prefix output was created");
			next unless -s $path;
			my $args = decode_json(read_file($path));
			my %seen = map { $_ => 1 } @$args;
			ok($seen{$mapF} && $seen{$assDirs} && $seen{"$OutD/$primaryClusterCLS"},
				"$prefix preserves literal paths and the selected catalogue identity");
			is(!!$seen{'-oldMapStyle'}, !!$oldNameFolders, "$prefix preserves the folder style");
			is(!!$seen{'-calcSupplCov'}, !!$CalcgGneMatSuppl, "$prefix preserves supplementary coverage selection");
			my @coverage = grep { /^-useCov/ } @$args;
			is_deeply(\@coverage, length($coverage_option) ? [$coverage_option] : [],
				"$prefix uses its intended abundance statistic");
			is(!!(-e "$OutD/$prefix.genes2rows.txt"), !length($coverage_option),
				"$prefix retains only the primary gene-to-row map");
		}
		ok(checkpoint_valid($matrixSton, parameters => {cluster_id => 97}),
			'completed outputs have a valid matrix checkpoint');
		unlink "$OutD/Mat.med.mat.gz" or die $!;
		ok(!checkpoint_valid($matrixSton, parameters => {cluster_id => 97}),
			'losing the median matrix invalidates the convergence checkpoint');
	});
}
done_testing();
