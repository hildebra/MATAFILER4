use strict;
use warnings;
use Test::More;
use File::Temp qw(tempdir);
use File::Path qw(make_path);
use File::Spec;
use FindBin qw($Bin);
use lib "$Bin/..", "$Bin/lib";
use MFTestConfig;
use Mods::GenoMetaAss qw(splitFastas median);
use Mods::WorkflowState qw(inspect_workflow_state);

# Regressions for the 2026-09-25 audit pass (second pass over the pipeline).
# Sub bodies are isolated from the real scripts, as in audit_2026_09_24.t.

my $root = File::Spec->rel2abs("$Bin/..");
my $tmp = tempdir(CLEANUP => 1);
$tmp =~ s{\\}{/}g;
sub write_file {
	my ($path, $text) = @_;
	(my $dir = $path) =~ s{/[^/]+$}{};
	make_path($dir) unless -d $dir;
	open my $fh, '>', $path or die "$path: $!";
	print {$fh} $text;
	close $fh or die "$path: $!";
}
sub read_file {
	open my $fh, '<', $_[0] or die "$_[0]: $!";
	local $/; return <$fh> // '';
}
sub source_sub {
	my ($source, $name) = @_;
	my ($code) = $source =~ /(^sub \Q$name\E\b[^\{;]*\{.*?^\})/ms;
	die "Cannot isolate $name" unless defined $code;
	return $code;
}
sub fasta_records { my $n = () = read_file($_[0]) =~ /^>/mg; return $n; }

# ---------------------------------------------------------------- splitFastas
# A split is reused only when its stone (written last) matches the input and
# every chunk is intact; truncated or stale chunk sets are rebuilt.
my $faa = "$tmp/q.faa";
write_file($faa, join('', map { ">g$_\nMKV\n" } 1..100));
my $splits = splitFastas($faa, 4, "$tmp/split");
my $total = 0; $total += fasta_records($_) for @{$splits};
is($total, 100, 'all records are distributed over the chunks');
ok(-s "$tmp/split/q.faa.4.split.sto", 'a finished split writes its stone');
is_deeply(splitFastas($faa, 4, "$tmp/split"), $splits, 'an intact split is reused');
write_file($splits->[-1], ">g100\nMKV\n"); # simulate a chunk truncated by an interrupted split
my $rebuilt = splitFastas($faa, 4, "$tmp/split");
$total = 0; $total += fasta_records($_) for @{$rebuilt};
is($total, 100, 'a truncated chunk set is split again');
write_file($faa, join('', map { ">g$_\nMKV\n" } 1..101));
utime(time - 100, time - 100, $faa);
$rebuilt = splitFastas($faa, 4, "$tmp/split");
$total = 0; $total += fasta_records($_) for @{$rebuilt};
is($total, 101, 'a changed input is split again');

# ------------------------------------------------------------ WorkflowState
my $sample_a = "$tmp/run/A/";
my $sample_b = "$tmp/run/B/";
my %map = (
	opt => { smpl_order => ['A', 'B'] },
	A => { SmplID => 'A', wrdir => $sample_a, AssGroup => 'gut', SupportReads => 'PB:/reads/a.fq.gz', hasPrimaryRds => 1 },
	B => { SmplID => 'B', wrdir => $sample_b, AssGroup => 'gut', SupportReads => '', hasPrimaryRds => 1 },
);
my %groups = (gut => { CntAimAss => 2 });
my $group_dir = "$tmp/run/AssmblGrp_gut/metag";
write_file("$group_dir/scaffolds.fasta.filt", ">contig\nACGT\n");
write_file("$group_dir/ass.done.sto", "");
write_file("$group_dir/smpls_used.txt", "$sample_a\n$sample_b\n");
write_file("$sample_a/mapping/A-smd.cram.sto", "");
write_file("$sample_a/mapping/A-smd.bam.coverage.gz", "coverage\n");
write_file("$sample_a/mapping/A-smd.bam.breakpoints.tsv.gz", "breakpoints\n");
my $kept = inspect_workflow_state(map => \%map, groups => \%groups);
is($kept->{samples}[0]{stages}{mapping}{status}, 'PARTIAL',
	'a missing CRAM is a partial mapping when CRAMs are kept');
my $removed = inspect_workflow_state(map => \%map, groups => \%groups,
	options => { mapping_cram_kept => 0 });
is($removed->{samples}[0]{stages}{mapping}{status}, 'COMPLETE',
	'a CRAM removed by finished-sample cleanup (-mapSaveCRAM 0 with a binner) is not a repair target');
is($removed->{workflow}{mapping_cram_kept}, 0, 'the CRAM policy is reported');

my $package = "${sample_a}assemblies/preAssmblGrp_gut";
write_file("$package/scaffolds.fasta.filt", ">preassembly\nACGT\n");
write_file("$package/Coverage.percontig.gz", "coverage\n");
write_file("$package/Coverage.median.percontig.gz", "median\n");
write_file("$package/mapping.coverage.gz", "mapping coverage\n");
write_file("$package/breakpoints.tsv.gz", "breakpoints\n");
write_file("$package/package.manifest.tsv", "key\tvalue\nschema_version\t2\n");
write_file("$package/moved.sto", "done\n");
my $hybrid = inspect_workflow_state(map => \%map, groups => \%groups,
	options => { assembly_mode => 5, run_tmp_dir => "$tmp/scratch" });
is($hybrid->{samples}[0]{stages}{preassembly_package}{status}, 'COMPLETE',
	'the hybrid handoff package is inspected under the sample output, not in run scratch');

# ------------------------------------------------------- MATAF4.pl subs
my $main = read_file("$root/MATAF4.pl");
our (%progStats, %MFopt);
eval source_sub($main, 'reduceProgStats');
die $@ if $@;
%MFopt = (DoMetaPhlan => 1, DoMOTU2 => 1, DoKraken => 1, DoProtal => 0);
%progStats = (metaPhl2FailCnts => 0, mOTU2FailCnts => 2, KrakTaxFailCnts => 1, protalFailCnts => 3);
reduceProgStats();
is_deeply(\%progStats, { metaPhl2FailCnts => 0, mOTU2FailCnts => 1, KrakTaxFailCnts => 0, protalFailCnts => 3 },
	'undoing an empty sample never pushes a failure counter below zero');

# command construction that cannot be executed here: guard the fixed fragments
like($main, qr/\$recoveredPattern = "\$prefix\.\[12\]\.singl\.\$fEnd"/,
	'sdm singleton recovery only globs the two mate files of this library');
like($main, qr/genePredictions\(\$refDB\[\$i\],\$gpDir,/,
	'secondary-reference gene prediction runs in a private subdirectory (genePredictions wipes its outDir)');
like($main, qr/smplName => join\(",",\@bamBaseNameS\)/,
	'secondary mapping names one intermediate BAM per reference');
like($main, qr/\} elsif \(\$porechopFlag && \$is3rdGen\)\{/,
	'porechop only runs for long-read samples whose staged path is recorded');
like($main, qr/getRgStr\(\$outNms\[0\]/, 'the read group uses a single output name');
like($main, qr/symlink\(\$refAbs, "\$bwt2outDl\/\$bwt2Name\[\$i\]\.fa"\)/,
	'secondary-mapping SNP calling finds the reference under <name>.fa');
like($main, qr/fileGZe\("\$dir_RibFind\/SSU\.miTag\.\$lvl\.txt"\)/,
	'the RiboFind merge skip-check accepts the gzipped tables the merger writes');
like($main, qr/\$MBcmd = "" if \(-e \$MetaBat2out && -s "\$BinDir\/Binning\.stone"\)/,
	'a binner assignment is only reused together with the binner stone');
my $strain = read_file("$root/secScripts/MGS/strain_within.pl");
like($strain, qr/if \(\$doSubmit\) \{\n\s*resetMGSTreeOutputs\(\$outD2, \$MGS\);/,
	'-submit 0 does not reset finished strain trees under -redo tree');

# --------------------------------------- annotateMGwSpecIs3.pl readMGS
{
	package MGSTaxProbe;
	use Mods::GenoMetaAss qw(median);
	our ($MGStax, %MGSlist, %gene2COG, %Gene2MGS, %SpecIgenes, %specIfullTax);
	my $src = main::read_file("$root/secScripts/GC/annotateMGwSpecIs3.pl");
	eval main::source_sub($src, 'readMGS');
	die $@ if $@;
	$MGStax = "$tmp/GTDBTK.tax";
	main::write_file($MGStax, "user_genome\tclassification\nMGS.1\tBacteria;Bacillota;Bacilli;Bacillales;;;\nMGS.2\tBacteria;;;;;;\n");
	main::write_file("$tmp/guide.MGS", "MGS.1\tg1\nMGS.2\tg2\n");
	%gene2COG = (g1 => 'COG0012', g2 => 'COG0016');
	{
		local *STDOUT; open STDOUT, '>', File::Spec->devnull;
		readMGS("$tmp/guide.MGS");
	}
	Test::More::is_deeply($specIfullTax{'MGS.1'}, [qw(Bacteria Bacillota Bacilli Bacillales ? ? ?)],
		'a family-level novel MGS gets "?" for every missing rank');
	Test::More::is_deeply($specIfullTax{'MGS.2'}, [qw(Bacteria ? ? ? ? ? ?)],
		'a phylum-level novel MGS gets "?" for every missing rank');
}

# --------------------------------------- phylo_MGS_between.pl inputs
my $between = read_file("$root/secScripts/MGS/phylo_MGS_between.pl");
like($between, qr/publishIfChanged->\("\$btout\/all\.faa\.tmp", "\$btout\/all\.faa"\)/,
	'the between-MGS launcher only replaces all.faa when its content changed');
unlike($between, qr/\(\$policy\{schema\} \/\/ ""\) eq "13"/,
	'the launcher no longer compares against one historical buildTree5 schema number');

done_testing;
