use strict;
use warnings;
use Test::More;
use File::Path qw(make_path);
use File::Temp qw(tempdir);

use lib '.';
use Mods::GenoMetaAss qw(
	parseSupportReads normaliseSupportReads discoverReadFiles readMap
	addFileLocs2AssmGrp getRawSeqsAssmGrp getCleanSeqsAssmGrp hasSuppRds
	getRawLibrariesAssmGrp getCleanLibrariesAssmGrp
);
use Mods::WorkflowControl qw(source_input_files);

my $tmp = tempdir(CLEANUP => 1);
my $read_dir = "$tmp/reads";
my $support_dir = "$tmp/support";
my $out_dir = "$tmp/out";
make_path($read_dir, $support_dir, $out_dir);

sub touch_file {
	my ($path, $contents) = @_;
	open(my $fh, '>', $path) or die "Cannot create $path: $!";
	print {$fh} defined($contents) ? $contents : "x\n";
	close($fh);
}

for my $file (qw(
	sample+1_L10_R1_001.fastq.gz sample+1_L10_R2_001.fastq.gz
	sample+1_L2_R1_001.fastq.gz sample+1_L2_R2_001.fastq.gz
	sample+1_single.fastq.gz unrelated_R1.fastq.gz
)) {
	touch_file("$read_dir/$file");
}

my $found = discoverReadFiles($read_dir, 'sample+1', {
	read1 => '.*R1.*\.fastq\.gz$',
	read2 => '.*R2.*\.fastq\.gz$',
	single => '.*single\.fastq\.gz$',
	bam => '', prefer_single => 0,
});
is_deeply(
	$found->{read1},
	[qw(sample+1_L2_R1_001.fastq.gz sample+1_L10_R1_001.fastq.gz)],
	'input discovery treats sample prefixes literally and naturally sorts lanes',
);
is_deeply(
	$found->{read2},
	[qw(sample+1_L2_R2_001.fastq.gz sample+1_L10_R2_001.fastq.gz)],
	'paired mates use the same natural lane order',
);
is_deeply($found->{single}, ['sample+1_single.fastq.gz'], 'unrelated files are excluded');
is(
	$found->{file_sizes}{'sample+1_L2_R1_001.fastq.gz'},
	-s "$read_dir/sample+1_L2_R1_001.fastq.gz",
	'input discovery returns the metadata gathered while classifying files',
);
ok(!exists $found->{file_sizes}{'unrelated_R1.fastq.gz'},
	'input discovery applies the sample prefix before gathering file metadata');

touch_file("$read_dir/sample+1_L3_R1_001.fastq.gz");
my $pair_error = eval {
	discoverReadFiles($read_dir, 'sample+1', {
		read1 => '.*R1.*\.fastq\.gz$', read2 => '.*R2.*\.fastq\.gz$',
		single => '', bam => '', prefer_single => 0,
	});
	1;
};
ok(!$pair_error, 'unequal mate counts are rejected instead of silently dropping all pairs');
like($@, qr/Unequal paired-read counts/, 'mate-count error identifies the input problem');

my $overlap_dir = "$tmp/overlap";
make_path($overlap_dir);
touch_file("$overlap_dir/sample_R1.fastq.gz");
eval {
	discoverReadFiles($overlap_dir, 'sample', {
		read1 => '.*fastq\.gz$', read2 => '.*fastq\.gz$',
		single => '', bam => '', prefer_single => 0,
	});
};
like($@, qr/match both read-1 and read-2 patterns/, 'ambiguous mate regular expressions are rejected');

my ($support_tech, $support_paths) = parseSupportReads(
	'PB:long/a.fastq.gz;PB:long/b.fastq.gz', $tmp,
);
is($support_tech, 'PB', 'support technology is parsed once');
is_deeply(
	$support_paths,
	["$tmp/long/a.fastq.gz", "$tmp/long/b.fastq.gz"],
	'repeated support specifications are normalized and relative paths use the map input root',
);
is(
	normaliseSupportReads('ONT:C:\\reads\\ont.fastq.gz'),
	'ONT:C:\\reads\\ont.fastq.gz',
	'colons in absolute drive paths are preserved',
);
eval { parseSupportReads('PB:a.fastq.gz;ONT:b.fastq.gz', $tmp) };
like($@, qr/mixes technologies/, 'mixed support technologies in one sample are rejected clearly');

touch_file("$support_dir/a.fastq.gz");
touch_file("$support_dir/b.fastq.gz");
my $map_file = "$tmp/test.map";
open(my $map_fh, '>', $map_file) or die "Cannot create map: $!";
print {$map_fh} join("\t", '#SmplID', qw(Path SeqTech SupportReads AssmblGrps)), "\n";
print {$map_fh} "#DirPath $tmp\n#OutPath $out_dir\n#RunID run\n";
print {$map_fh} join("\t", 'S1', 'reads', 'ill',
	'PB:support/a.fastq.gz;PB:support/b.fastq.gz', 'groupA'), "\n";
close($map_fh);
my ($map, $groups) = readMap($map_file, -1, {}, {}, 0);
is(
	$map->{S1}{SupportReads},
	"PB:$tmp/support/a.fastq.gz,$tmp/support/b.fastq.gz",
	'mapping-file support inputs are stored in one canonical representation',
);
is(
	$groups->{groupA}{SupportReads},
	$map->{S1}{SupportReads},
	'assembly-group support summary has no synthetic leading comma',
);

my $transfer_parent = "$tmp/Ileal_samples";
my $transfer_dir = "$transfer_parent/Transfer";
make_path($transfer_dir);
touch_file("$transfer_dir/IL7_1.fq.gz");
touch_file("$transfer_dir/IL7_2.fq.gz");
my $prefix_map_file = "$tmp/prefix_input_root.map";
open(my $prefix_map_fh, '>', $prefix_map_file) or die "Cannot create map: $!";
print {$prefix_map_fh} join("\t", '#SmplID', qw(SmplPrefix SeqTech)), "\n";
print {$prefix_map_fh} "#DirPath\t$transfer_dir/\n#OutPath\t$out_dir\n#RunID\tprefix_run\n";
print {$prefix_map_fh} join("\t", 'IL7', 'IL7_', 'ill'), "\n";
close($prefix_map_fh);
my ($prefix_map) = readMap($prefix_map_file, -1, {}, {}, 0);
is($prefix_map->{IL7}{rddir}, "$transfer_dir/",
	'#DirPath remains the primary input directory for SmplPrefix-only map rows');
my $prefix_found = discoverReadFiles($prefix_map->{IL7}{rddir}, $prefix_map->{IL7}{prefix}, {
	read1 => '.*1\\.fq\\.gz$', read2 => '.*2\\.fq\\.gz$',
	single => '', bam => '', prefer_single => 0,
});
is_deeply(
	source_input_files($prefix_map->{IL7}{rddir}, @{$prefix_found->{read1}}),
	["$transfer_dir/IL7_1.fq.gz"],
	'unzip source resolution retains the complete map #DirPath including its final directory',
);

my $bam_dir = "$tmp/bam_reads";
make_path($bam_dir);
touch_file("$bam_dir/QIBS1.hifi_reads.bam");
touch_file("$bam_dir/QIBS10.hifi_reads.bam");
my $bam_map_file = "$tmp/bam_only.map";
open(my $bam_map_fh, '>', $bam_map_file) or die "Cannot create map: $!";
print {$bam_map_fh} join("\t", '#SmplID', qw(SmplPrefix SeqTech SupportReads AssmblGrps)), "\r\n";
print {$bam_map_fh} "#DirPath $bam_dir\r\n#OutPath $out_dir\r\n#RunID bam_run\r\n";
print {$bam_map_fh} join("\t", 'QIBS1.', 'QIBS1.hifi_reads.bam', 'PB'), "\r\n";
print {$bam_map_fh} join("\t", 'QIBS10', 'QIBS10.hifi_reads.bam', 'PB'), "\r\n";
close($bam_map_fh);
my $bam_map_warnings = '';
my ($bam_map, $bam_groups);
{
	local $SIG{__WARN__} = sub { $bam_map_warnings .= join('', @_); };
	($bam_map, $bam_groups) = readMap($bam_map_file, -1, {}, {}, 0);
}
is($bam_map_warnings, '', 'BAM-only map rows may omit unused trailing columns without warnings');
is($bam_map->{'QIBS1.'}{SupportReads}, '', 'omitted SupportReads is treated as empty');
is($bam_map->{'QIBS1.'}{AssGroup}, 0, 'omitted AssmblGrps gets an automatic group');
my $bam_found = discoverReadFiles($bam_dir, 'QIBS1.', {
	read1 => '', read2 => '', single => '', bam => '.*\.bam$', prefer_single => 0,
});
is_deeply(
	$bam_found->{bam},
	['QIBS1.hifi_reads.bam'],
	'BAM-only primary input is discovered from its sample prefix',
);

my %assembly_groups = (groupA => { CleanSeqs => {}, RawSeqs => {}, InputOrder => [] });
my $raw_b = {
	pa1 => ['B.1.fq.gz'], pa2 => ['B.2.fq.gz'], pas => [], libInfo => ['B'], seqTech => 'ill',
	paX1 => [], paX2 => [], paXs => [], libInfoX => [], seqTechX => '',
};
my $clean_b = {
	arp1 => ['B.1.clean.fq.gz'], arp2 => ['B.2.clean.fq.gz'], singAr => [], readTec => 'ill',
	arpX1 => [], arpX2 => [], singArX => [], readTecX => '',
};
my $raw_a = {
	pa1 => ['A.1.fq.gz'], pa2 => ['A.2.fq.gz'], pas => [], libInfo => ['A'], seqTech => 'ill',
	paX1 => [], paX2 => [], paXs => ['A.long.fq.gz'], libInfoX => ['ONT'], seqTechX => 'ONT',
};
my $clean_a = {
	arp1 => ['A.lane1.1.clean.fq.gz', 'A.lane2.1.clean.fq.gz'],
	arp2 => ['A.lane1.2.clean.fq.gz', 'A.lane2.2.clean.fq.gz'],
	singAr => [], readTec => 'ill', arpX1 => [], arpX2 => [],
	singArX => ['A.long.clean.fq.gz'], readTecX => 'ONT',
};
addFileLocs2AssmGrp(\%assembly_groups, 'groupA', 'B', $clean_b, $raw_b);
addFileLocs2AssmGrp(\%assembly_groups, 'groupA', 'A', $clean_a, $raw_a);
my $raw_libraries = getRawLibrariesAssmGrp(\%assembly_groups, 'groupA', 0);
is_deeply([map { $_->{sample} } @{$raw_libraries}], [qw(B A)],
	'assembly groups expose ordered explicit raw-library records');
my $clean_libraries = getCleanLibrariesAssmGrp(\%assembly_groups, 'groupA', 0);
is_deeply([map { $_->{files}{r1} } @{$clean_libraries}],
	['B.1.clean.fq.gz', 'A.lane1.1.clean.fq.gz', 'A.lane2.1.clean.fq.gz'],
	'assembly groups aggregate explicit clean-library records without truncation');
my ($raw_r1) = getRawSeqsAssmGrp(\%assembly_groups, 'groupA', 0);
is_deeply($raw_r1, ['B.1.fq.gz', 'A.1.fq.gz'], 'assembly-group raw inputs retain map/insertion order');
my ($clean_r1, $clean_r2, undef, $clean_tech) = getCleanSeqsAssmGrp(\%assembly_groups, 'groupA', 0);
is_deeply(
	$clean_r1,
	['B.1.clean.fq.gz', 'A.lane1.1.clean.fq.gz', 'A.lane2.1.clean.fq.gz'],
	'multiple cleaned libraries are retained rather than truncating to the first file',
);
is_deeply($clean_r2, ['B.2.clean.fq.gz', 'A.lane1.2.clean.fq.gz', 'A.lane2.2.clean.fq.gz'], 'clean mate arrays remain aligned');
is_deeply($clean_tech, [qw(ill ill ill)], 'read technology remains aligned with expanded clean libraries');
ok(!hasSuppRds(\%assembly_groups, 'groupA', 'B'), 'empty support arrays do not count as support reads');
ok(hasSuppRds(\%assembly_groups, 'groupA', 'A'), 'non-empty support input is detected');

done_testing();
