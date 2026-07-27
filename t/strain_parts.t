use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Test::More;

use lib File::Spec->catdir($Bin, '..');
use Mods::StrainParts qw(
	exact_worker_parts write_split_generation write_worker_completion
	split_generation_complete clear_split_generation
	resolve_fasta_artifact append_fasta_records_atomic
	sort_fasta_by_locus
);

sub write_file {
	my ($path, $contents) = @_;
	open my $fh, '>', $path or die "Cannot write $path: $!";
	print {$fh} $contents;
	close $fh or die "Cannot close $path: $!";
}

my $tmp = tempdir(CLEANUP => 1);
my $prefix = File::Spec->catfile($tmp, 'allFNAs.fna');
write_file("$prefix.0", "worker zero\n");
write_file("$prefix.2", "worker two\n");
write_file("$prefix.3", "old out-of-range worker\n");
write_file("$prefix.merge.1234", "abandoned partial merge\n");
write_file("$prefix.not-a-worker", "noise\n");

is_deeply(
	[exact_worker_parts($prefix, 3)],
	["$prefix.0", "$prefix.2"],
	'only exact in-range numeric worker suffixes are merge inputs',
);

my $manifest = File::Spec->catfile($tmp, 'mainExtr.generation');
my $stone_prefix = File::Spec->catfile($tmp, 'mainExtr');
write_split_generation($manifest, 'gen.1', 3);
write_worker_completion("$stone_prefix.0.stone", 'gen.1');
write_worker_completion("$stone_prefix.1.stone", 'gen.1');
ok(!split_generation_complete($manifest, $stone_prefix, 3),
	'a generation is incomplete while a worker completion is missing');

write_worker_completion("$stone_prefix.2.stone", 'older.gen');
ok(!split_generation_complete($manifest, $stone_prefix, 3),
	'a completion from another generation cannot complete the set');

write_worker_completion("$stone_prefix.2.stone", 'gen.1');
ok(split_generation_complete($manifest, $stone_prefix, 3),
	'all workers from the same generation complete the set');
ok(!split_generation_complete($manifest, $stone_prefix, 2),
	'a manifest made for another worker count is rejected');

write_file("$stone_prefix.not-a-worker.stone", "keep\n");
clear_split_generation($manifest, $stone_prefix);
ok(!-e $manifest && !-e "$stone_prefix.0.stone" && !-e "$stone_prefix.2.stone",
	'clearing generation state removes its manifest and numeric worker stones');
ok(-e "$stone_prefix.not-a-worker.stone",
	'clearing generation state leaves unrelated files alone');

my $plain_fasta = File::Spec->catfile($tmp, 'plain.fna');
write_file($plain_fasta, ">sample\nACGT");
is(append_fasta_records_atomic($plain_fasta, ">outgroup\nTGCA\n"), $plain_fasta,
	'plain FASTA append keeps the canonical artifact path');
open my $plain_fh, '<', $plain_fasta or die "Cannot read $plain_fasta: $!";
my $plain_contents = do { local $/; <$plain_fh> };
close $plain_fh;
is($plain_contents, ">sample\nACGT\n>outgroup\nTGCA\n",
	'plain FASTA rewrite separates and preserves existing records');

my $gzip_nominal = File::Spec->catfile($tmp, 'compressed.faa');
gzip \">sample\nMPEP\n" => "$gzip_nominal.gz"
	or die "Cannot create gzip fixture: $GzipError";
is(resolve_fasta_artifact($gzip_nominal), "$gzip_nominal.gz",
	'FASTA resolution exposes the existing compressed artifact for duplicate checks');
is(append_fasta_records_atomic($gzip_nominal, ">outgroup\nMTEST\n"), "$gzip_nominal.gz",
	'compressed FASTA append resolves and retains the gzip artifact');
ok(!-e $gzip_nominal, 'compressed append does not create a shadowing plain sidecar');
my $gzip_contents = '';
gunzip "$gzip_nominal.gz" => \$gzip_contents
	or die "Cannot read gzip fixture: $GunzipError";
is($gzip_contents, ">sample\nMPEP\n>outgroup\nMTEST\n",
	'compressed FASTA rewrite preserves original and appended records');

my $ambiguous = File::Spec->catfile($tmp, 'ambiguous.fna');
write_file($ambiguous, ">plain\nAAAA\n");
gzip \">gzip\nCCCC\n" => "$ambiguous.gz"
	or die "Cannot create ambiguous gzip fixture: $GzipError";
eval { append_fasta_records_atomic($ambiguous, ">new\nGGGG\n") };
like($@, qr/Ambiguous FASTA sidecars/, 'ambiguous plain/gzip sidecars fail instead of corrupting either copy');

my $unsorted_fasta = File::Spec->catfile($tmp, 'unsorted.fna');
write_file($unsorted_fasta,
	">sampleB|COG2|gene9\nCCCC\n"
	. ">sampleC|COG1|gene2\nGG\nGG\n"
	. ">sampleA|COG1|gene2\nAAAA\n"
	. ">sampleA|COG1|gene1\nTTTT\n");
is(sort_fasta_by_locus($unsorted_fasta, '|'), 4,
	'locus sorting reports its FASTA record count');
open my $sorted_fh, '<', $unsorted_fasta or die "Cannot read $unsorted_fasta: $!";
my $sorted_contents = do { local $/; <$sorted_fh> };
close $sorted_fh;
is($sorted_contents,
	">sampleA|COG1|gene1\nTTTT\n"
	. ">sampleA|COG1|gene2\nAAAA\n"
	. ">sampleC|COG1|gene2\nGG\nGG\n"
	. ">sampleB|COG2|gene9\nCCCC\n",
	'FASTA records sort by eggNOG, gene-catalogue ID, then sample without changing multiline sequences');

done_testing();
