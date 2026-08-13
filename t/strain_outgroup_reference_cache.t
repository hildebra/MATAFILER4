use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::GenoMetaAss qw(ensureFastaIndex readFasta writeFasta);

sub write_file {
	my ($path, $contents) = @_;
	open my $output, '>', $path or die "Cannot create $path: $!";
	print {$output} $contents or die "Cannot write $path: $!";
	close $output or die "Cannot close $path: $!";
}

sub slurp {
	my ($path) = @_;
	open my $input, '<', $path or die "Cannot read $path: $!";
	local $/;
	my $contents = <$input> // '';
	close $input or die "Cannot close $path: $!";
	return $contents;
}

sub run_helper {
	my (@arguments) = @_;
	open my $process, '-|', $^X, @arguments
		or die "Cannot start outgroup-reference helper: $!";
	local $/;
	my $output = <$process> // '';
	my $closed = close $process;
	return ($closed && $? == 0, $output, $?);
}

my $temporary = tempdir(CLEANUP => 1);
my $fake_samtools = File::Spec->catfile($temporary, 'samtools');
write_file($fake_samtools, <<'FAKE');
#!/usr/bin/env perl
use strict;
use warnings;
my ($operation, @arguments) = @ARGV;
die "Only faidx is supported\n" unless ($operation // '') eq 'faidx';
if (@arguments == 1) {
	my ($fasta) = @arguments;
	open my $input, '<', $fasta or die $!;
	open my $index, '>', "$fasta.fai" or die $!;
	my ($identifier, $sequence) = ('', '');
	my $flush = sub {
		return unless length $identifier;
		print {$index} join("\t", $identifier, length($sequence), 0,
			length($sequence), length($sequence) + 1), "\n";
	};
	while (my $line = <$input>) {
		chomp $line;
		if ($line =~ /^>(\S+)/) {
			$flush->();
			($identifier, $sequence) = ($1, '');
		} else {
			$sequence .= $line;
		}
	}
	$flush->();
	close $input or die $!;
	close $index or die $!;
	exit 0;
}
if (@arguments == 3 && $arguments[0] eq '-r') {
	my (undef, $ids, $fasta) = @arguments;
	open my $input, '<', $fasta or die $!;
	my (%sequence, $identifier);
	while (my $line = <$input>) {
		chomp $line;
		if ($line =~ /^>(\S+)/) {
			$identifier = $1;
			$sequence{$identifier} = '';
		} elsif (defined $identifier) {
			$sequence{$identifier} .= $line;
		}
	}
	close $input or die $!;
	open my $wanted, '<', $ids or die $!;
	while (my $id = <$wanted>) {
		chomp $id;
		next unless exists $sequence{$id};
		print ">$id\n$sequence{$id}\n";
	}
	close $wanted or die $!;
	exit 0;
}
die "Unsupported fake samtools invocation: @ARGV\n";
FAKE
chmod 0755, $fake_samtools or die "Cannot make fake samtools executable: $!";

my $adopted = File::Spec->catfile($temporary, 'adopted.fna');
write_file($adopted, ">adopt1\nAAAA\n>adopt2\nCC\n");
my $adopted_index = join("",
	"adopt1\t4\t8\t4\t5\n",
	"adopt2\t2\t21\t2\t3\n",
);
write_file("$adopted.fai", $adopted_index);
is(
	ensureFastaIndex($adopted),
	"$adopted.fai",
	'a complete existing plain-FASTA index is adopted without rerunning samtools',
);
is(slurp("$adopted.fai"), $adopted_index,
	'adopting an existing complete index leaves its content unchanged');
ok(-s "$adopted.fai.complete",
	'adopted direct-samtools indexes receive a shared completion fingerprint');

my $nt = File::Spec->catfile($temporary, 'catalogue.fna');
my $aa = File::Spec->catfile($temporary, 'catalogue.faa');
write_file($nt, ">g1\nAAAA\n>g2\nCCCC\n>g3\nGGGG\n");
write_file($aa, ">g1\nKK\n>g2\nPP\n>g3\nGG\n");
my $nt_ids = File::Spec->catfile($temporary, 'nt.ids');
my $aa_ids = File::Spec->catfile($temporary, 'aa.ids');
write_file($nt_ids, "g2\nmissing\n");
write_file($aa_ids, "g1\ng2\nmissing\n");
my $cache = File::Spec->catdir($temporary, 'cache');
my $helper = File::Spec->catfile($Bin, '..', 'secScripts', 'MGS',
	'prepare_strain_outgroup_refs.pl');
my @arguments = ($helper, '-nt', $nt, '-aa', $aa, '-ntIDs', $nt_ids,
	'-aaIDs', $aa_ids, '-outD', $cache, '-samtools', $fake_samtools);
my ($ok, $output, $status) = run_helper(@arguments);
ok($ok, "indexed outgroup-reference helper succeeds (status=$status)");
like($output, qr/status=created requested_nt=2 retained_nt=1 requested_aa=3 retained_aa=2/,
	'creation summary reports requested, retained, and implicitly missing references');
is(slurp(File::Spec->catfile($cache, 'references.fna')), ">g2\nCCCC\n",
	'compact nucleotide cache contains only an available requested outgroup gene');
is(slurp(File::Spec->catfile($cache, 'references.faa')), ">g1\nKK\n>g2\nPP\n",
	'compact protein cache contains outgroup and target similarity proteins only');
for my $path (
	"$nt.fai", "$nt.fai.complete", "$aa.fai", "$aa.fai.complete",
	File::Spec->catfile($cache, 'references.fna.fai'),
	File::Spec->catfile($cache, 'references.faa.fai'),
	File::Spec->catfile($cache, 'manifest.tsv'),
	File::Spec->catfile($cache, 'complete.sto'),
) {
	ok(-s $path, "indexed cache artifact is nonempty: ".File::Spec->abs2rel($path, $temporary));
}

($ok, $output, $status) = run_helper(@arguments);
ok($ok, "matching outgroup-reference cache is reusable (status=$status)");
like($output, qr/status=reused requested_nt=2 retained_nt=1 requested_aa=3 retained_aa=2/,
	'an identical requirement set reuses the completed indexed cache');

# Simulate an interrupted source-index build: a nonempty partial index without
# its completion fingerprint must be rebuilt, never adopted.
write_file("$nt.fai", "g1\t4\t0\t4\t5\n");
unlink "$nt.fai.complete" or die "Cannot remove source-index completion marker: $!";

write_file($nt_ids, "g3\n");
write_file($aa_ids, "g3\n");
($ok, $output, $status) = run_helper(@arguments);
ok($ok, "changed outgroup requirements rebuild the cache (status=$status)");
like($output, qr/status=created requested_nt=1 retained_nt=1 requested_aa=1 retained_aa=1/,
	'a changed requirement digest invalidates the prior cache');
is(slurp(File::Spec->catfile($cache, 'references.fna')), ">g3\nGGGG\n",
	'rebuilt nucleotide cache does not retain obsolete requested genes');
is(slurp(File::Spec->catfile($cache, 'references.faa')), ">g3\nGG\n",
	'rebuilt protein cache does not retain obsolete requested genes');


my @indexed_progress;
my %indexed_wanted = (g1 => 1, g3 => 1, missing => 1);
is_deeply(
	readFasta($nt, 1, '\\s', \%indexed_wanted,
		sub { push @indexed_progress, { %{$_[0]} } },
		{ fai => 1, require_fai => 1, samtools => $fake_samtools }),
	{ g1 => 'AAAA', g3 => 'GGGG' },
	'readFasta uses the shared index for an explicitly indexed subset read',
);
is($indexed_progress[-1]{mode}, 'fai',
	'indexed read progress identifies the FAI path');
is($indexed_progress[-1]{records_retained}, 2,
	'indexed read progress reports available requested records');

my @large_subset_progress;
is_deeply(
	readFasta($nt, 1, '\\s', { g1 => 1, g2 => 1, g3 => 1 },
		sub { push @large_subset_progress, { %{$_[0]} } },
		{ fai => 1, fai_max_fraction => 0.25, samtools => $fake_samtools }),
	{ g1 => 'AAAA', g2 => 'CCCC', g3 => 'GGGG' },
	'readFasta keeps a large catalogue subset on the faster streaming path',
);
ok(scalar(grep { ($_->{mode} // '') eq 'fai_streaming_preferred' }
	@large_subset_progress),
	'indexed subset planning reports when sequential streaming is preferred');

my $described = File::Spec->catfile($temporary, 'described.fna');
write_file($described, ">g1 description retained\nACGT\n>g2 another description\nTGCA\n");
my %described_wanted = (g1 => 1);
is_deeply(
	readFasta($described, 0, '\\s', \%described_wanted,
		{ fai => 1, samtools => $fake_samtools }),
	{ 'g1 description retained' => 'ACGT' },
	'full-header reads retain streaming semantics when FAI output cannot preserve headers',
);

my $written = File::Spec->catfile($temporary, 'written.fna');
writeFasta({ w1 => 'AAAA', w2 => 'CCCC' }, $written,
	{ fai => 1, samtools => $fake_samtools, quiet => 1 });
ok(-s "$written.fai" && -s "$written.fai.complete",
	'writeFasta creates a reusable index when explicitly requested');
is_deeply(
	readFasta($written, 1, '\\s', { w2 => 1 },
		{ fai => 1, require_fai => 1, samtools => $fake_samtools }),
	{ w2 => 'CCCC' },
	'writeFasta output is immediately available through indexed readFasta',
);
writeFasta({ w1 => 'TTTT', w2 => 'GGGG' }, $written,
	{ fai => 1, samtools => $fake_samtools, quiet => 1 });
is_deeply(
	readFasta($written, 1, '\\s', { w1 => 1 },
		{ fai => 1, require_fai => 1, samtools => $fake_samtools }),
	{ w1 => 'TTTT' },
	'writeFasta forcibly refreshes an index even when size and timestamp granularity match',
);

my $unindexed_large = File::Spec->catfile($temporary, 'unindexed-large.fna');
write_file($unindexed_large, ">u1\nAAAA\n>u2\nCCCC\n");
is_deeply(
	readFasta($unindexed_large, 1, '\\s', { u1 => 1, u2 => 1 },
		{ fai => 1, fai_max_ids_without_index => 1,
			samtools => "$temporary/not-needed-samtools" }),
	{ u1 => 'AAAA', u2 => 'CCCC' },
	'a large first read streams instead of building an index and then rereading the FASTA',
);
ok(!-e "$unindexed_large.fai",
	'the first-read guard does not pay the index construction cost for a large subset');

my $fallback_source = File::Spec->catfile($temporary, 'fallback.fna');
write_file($fallback_source, ">fallback\nAACCGG\n");
my @fallback_warnings;
my $fallback_sequences;
{
	local $SIG{__WARN__} = sub { push @fallback_warnings, @_ };
	$fallback_sequences = readFasta(
		$fallback_source, 1, '\\s', { fallback => 1 },
		{ fai => 1, samtools => "$temporary/missing-samtools" },
	);
}
is_deeply($fallback_sequences, { fallback => 'AACCGG' },
	'an unavailable indexed reader falls back to the existing streaming parser');
like(join('', @fallback_warnings), qr/using streaming scan/,
	'indexed-reader fallback is visible without aborting an automatic workflow');
done_testing;
