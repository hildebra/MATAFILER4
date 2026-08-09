use strict;
use warnings;
use Test::More;

use lib '.';
use Mods::ReadLibrary qw(
	newReadLibrary readLibrariesFromArrays validateReadLibraries
	ensureSeqSetLibraries ensureCleanSeqSetLibraries
	syncSeqSetLegacy syncCleanSeqSetLegacy replaceScopeLibraries
	readLibrariesByScope libraryFiles libraryPairs legacyLibraryArrays
	libraryTechnology
);

my $primary = readLibrariesFromArrays(
	sample => 'S1', scope => 'primary', phase => 'staged', technology => 'ill',
	r1 => [qw(S1.L1.R1.fq.gz S1.L2.R1.fq.gz)],
	r2 => [qw(S1.L1.R2.fq.gz S1.L2.R2.fq.gz)],
	single => ['S1.single.fq.gz'], labels => [qw(lib1 lib2)],
);
is(scalar(@{$primary}), 2, 'parallel legacy inputs become explicit library records');
is($primary->[0]{files}{r1}, 'S1.L1.R1.fq.gz', 'record owns its read-1 path');
is($primary->[0]{files}{r2}, 'S1.L1.R2.fq.gz', 'record owns the corresponding read-2 path');
is($primary->[0]{files}{single}, 'S1.single.fq.gz', 'a library can retain its associated singleton stream');
is($primary->[1]{technology}, 'ill', 'technology is attached to each library rather than a side array');

my $role_separated = readLibrariesFromArrays(
	sample => 'mixed', scope => 'primary', phase => 'staged', separate_roles => 1,
	pair_technology => 'hiSeq', single_technology => 'ONT',
	r1 => ['mixed.R1.fq.gz'], r2 => ['mixed.R2.fq.gz'],
	single => ['mixed.ont.fq.gz'], labels => ['source'],
);
is(scalar(@{$role_separated}), 2, 'paired and singleton sources become distinct records when roles are separate');
is_deeply($role_separated->[0]{files}, {
	r1 => 'mixed.R1.fq.gz', r2 => 'mixed.R2.fq.gz', single => '', bam => '',
}, 'paired record contains only its two linked mates');
is_deeply($role_separated->[1]{files}, {
	r1 => '', r2 => '', single => 'mixed.ont.fq.gz', bam => '',
}, 'singleton record does not inherit unrelated paired files');
is_deeply([map { $_->{technology} } @{$role_separated}], [qw(hiSeq ONT)],
	'each separated input retains its own sequencing technology');
like($role_separated->[1]{label}, qr/\.single$/, 'separated singleton receives an unambiguous library label');
my $mixed_seqset = {libraries => $role_separated};
syncSeqSetLegacy($mixed_seqset);
is($mixed_seqset->{seqTech}, 'hiSeq', 'legacy scalar technology remains the first-library compatibility value');
eval { libraryTechnology($role_separated, 'mapping mixed') };
like($@, qr/mapping mixed.*hiSeq.*ONT|mapping mixed.*ONT.*hiSeq/,
	'mixed-technology consumers fail explicitly instead of silently choosing a mapper');
my $unknown_technology = [newReadLibrary(
	id => 'unknown-tech', scope => 'primary', files => {single => 'unknown.fq.gz'},
)];
eval { libraryTechnology($unknown_technology, 'mapping unknown', 1) };
like($@, qr/mapping unknown.*unknown-tech/, 'mapping rejects records whose technology is missing');

my $support = [newReadLibrary(
	id => 'S1:support:0', sample => 'S1', scope => 'support', technology => 'ONT',
	is_long => 1, label => 'ont', phase => 'staged',
	files => {single => 'S1.ont.fq.gz'},
)];
my $derived_long = newReadLibrary(
	id => 'S1:support:derived', sample => 'S1', scope => 'support', technology => 'ONT',
	label => 'ont-derived', files => {single => 'derived.ont.fq.gz'},
);
is($derived_long->{is_long}, 1, 'ONT/PacBio technology automatically marks a record as long-read');
my $seqset = {
	libraries => [@{$primary}, @{$support}], rawReads => 'inventory',
	totalInputSizeMB => 10, inputXFileSizeMB => 20,
};
syncSeqSetLegacy($seqset);
is_deeply($seqset->{pa1}, [qw(S1.L1.R1.fq.gz S1.L2.R1.fq.gz)], 'legacy primary read-1 view is regenerated');
is_deeply($seqset->{paXs}, ['S1.ont.fq.gz'], 'legacy support singleton view is regenerated');
is($seqset->{seqTechX}, 'ONT', 'legacy support technology remains available');
is($seqset->{rawReads}, 'inventory', 'non-read map compatibility metadata is retained');

my $from_legacy = {
	pa1 => ['old.R1.fq.gz'], pa2 => ['old.R2.fq.gz'], pas => [],
	libInfo => ['old'], seqTech => 'miSeq', is3rdGen => 0,
	paX1 => [], paX2 => [], paXs => [], libInfoX => [], seqTechX => '', is3rdGenX => 0,
};
my $converted = ensureSeqSetLibraries($from_legacy, 'oldSample');
is(scalar(@{$converted}), 1, 'legacy seqSet is converted lazily for compatibility');
is($converted->[0]{sample}, 'oldSample', 'compatibility conversion records the sample');
is($converted->[0]{technology}, 'miSeq', 'compatibility conversion preserves technology');

my $clean = {
	arp1 => ['clean.R1.fq.gz'], arp2 => ['clean.R2.fq.gz'], singAr => [], matAr => ['clean'],
	readTec => 'ill', is3rdGen => 0,
	arpX1 => [], arpX2 => [], singArX => [], matArX => [], readTecX => '', is3rdGenX => 0,
};
my $clean_records = ensureCleanSeqSetLibraries($clean, 'S1');
is($clean_records->[0]{phase}, 'clean', 'legacy clean state becomes clean-phase records');
my $replacement = readLibrariesFromArrays(
	sample => 'S1', scope => 'support', clean => 1, phase => 'clean', technology => 'PB', is_long => 1,
	single => [qw(pb.1.fq.gz pb.2.fq.gz)], labels => [qw(pb1 pb2)],
);
replaceScopeLibraries($clean, 'support', $replacement, 1, 'S1');
is_deeply(libraryFiles(readLibrariesByScope($clean, 'support', 1, 'S1'), 'single'),
	[qw(pb.1.fq.gz pb.2.fq.gz)], 'a scope is replaced as records');
is_deeply($clean->{singArX}, [qw(pb.1.fq.gz pb.2.fq.gz)], 'clean compatibility view follows record replacement');
$clean->{merged_library} = newReadLibrary(
	id => 'S1:primary:merged', sample => 'S1', scope => 'primary', technology => 'ill',
	label => 'merged', phase => 'merged',
	files => {r1 => 'unmerged.R1.fq.gz', r2 => 'unmerged.R2.fq.gz', single => 'merged.fq.gz'},
);
syncCleanSeqSetLegacy($clean);
is_deeply($clean->{mrgHshHR}, {
	pair1 => 'unmerged.R1.fq.gz', pair2 => 'unmerged.R2.fq.gz', mrg => 'merged.fq.gz',
}, 'merged library has a regenerated legacy merge view');

is(scalar(@{libraryPairs($primary)}), 2, 'paired-library selection operates on records');
my ($r1, $r2, $single, $labels, $technologies) = legacyLibraryArrays($primary, 1);
is_deeply($technologies, [qw(ill ill)], 'aligned compatibility projection repeats per-record technology');

eval { newReadLibrary(id => 'broken', scope => 'primary', files => {r1 => 'only.R1.fq.gz'}) };
like($@, qr/only one mate/, 'half-paired records are rejected at construction');
eval { validateReadLibraries([$primary->[0], $primary->[0]]) };
like($@, qr/Duplicate read library id/, 'duplicate record identities are rejected');
my $wrong_phase = {libraries => [newReadLibrary(
	id => 'wrong-phase', sample => 'S1', scope => 'primary', phase => 'clean',
	technology => 'hiSeq', files => {single => 'wrong.fq.gz'},
)]};
eval { ensureSeqSetLibraries($wrong_phase, 'S1') };
like($@, qr/Raw seqSet contains non-raw library/, 'raw consumers reject cleaned library records');
my $wrong_sample = {libraries => [newReadLibrary(
	id => 'wrong-sample', sample => 'S2', scope => 'primary', phase => 'staged',
	technology => 'hiSeq', files => {single => 'wrong-sample.fq.gz'},
)]};
eval { ensureSeqSetLibraries($wrong_sample, 'S1') };
like($@, qr/owned by 'S2'/, 'sample-scoped consumers reject a library owned by another sample');

done_testing();
