use strict;
use warnings;
use Test::More;
use File::Basename qw(dirname);
use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IPC::Open3;
use Symbol qw(gensym);
use lib "$Bin/..", "$Bin/lib";
use MFTestConfig;
use Mods::GenoMetaAss qw(fileGZe);
use Mods::RibosomeState qw(
	compress_ribosome_hierarchies lca_reference_copy_current normalise_ribosome_request
	ribosome_completion_evidence ribosome_merge_cohort ribosome_merge_manifest
);
use Mods::SampleCompletion qw(read_sample_completion sample_completion_path write_sample_completion);
use Mods::WorkflowControl qw(sample_is_ignored);
use Mods::WorkflowResilience qw(atomic_write_text);

# Regressions for the 2026-10-08 RiboFind audit: read extraction
# (catchLSUSSU.pl), LCA assignment (lotus_LCA_blast3.pl), the cohort merge
# (miTagTaxTable.pl) and their MATAF4 orchestration.

my $root = File::Spec->rel2abs("$Bin/..");
my $scripts = "$root/secScripts/miTag";
my $tmp = tempdir(CLEANUP => 1);
$tmp =~ s{\\}{/}g;
local $ENV{PERL5OPT} = join ' ', grep { defined($_) && length($_) }
	"-I$root", "-I$root/t/lib", '-MMFTestConfig', $ENV{PERL5OPT};

sub write_file {
	my ($path, $text) = @_;
	make_path(dirname($path)) unless -d dirname($path);
	open my $fh, '>', $path or die "$path: $!";
	print {$fh} $text;
	close $fh or die "$path: $!";
}
sub write_gzip {
	my ($path, $text) = @_;
	make_path(dirname($path)) unless -d dirname($path);
	gzip(\$text => $path) or die "$path: $GzipError";
}
sub read_file {
	open my $fh, '<', $_[0] or die "$_[0]: $!";
	local $/; return <$fh> // '';
}
sub read_gzip {
	my $text = '';
	gunzip($_[0] => \$text, MultiStream => 1) or die "$_[0]: $GunzipError";
	return $text;
}
sub lines { return -e $_[0] ? scalar(() = read_file($_[0]) =~ /\n/g) : 0 }
sub run_script {
	my ($script, @arguments) = @_;
	my $error = gensym;
	my $pid = open3(undef, my $stdout, $error,
		$^X, "-I$root", "$scripts/$script", @arguments);
	my $output = do { local $/; <$stdout> } // '';
	my $errors = do { local $/; <$error> } // '';
	waitpid($pid, 0);
	return ($? >> 8, $output, $errors);
}
sub tool {
	my ($name, $body) = @_;
	write_file("$tmp/bin/$name", $body);
	chmod 0755, "$tmp/bin/$name" or die "$name: $!";
	return "$tmp/bin/$name";
}

my $sortLog = "$tmp/sortmerna.log";
my $queryLog = "$tmp/queries.log";
my $sortmerna = tool('sortmerna', qq{#!/usr/bin/env perl
use IO::Compress::Gzip qw(gzip);
exit 0 if grep { \$_ eq '--version' } \@ARGV;
my (\$aligned, \$paired, \$blast, \$reads1) = ('', 0, 0, '');
for (my \$i = 0; \$i < \@ARGV; \$i++) {
	\$aligned = \$ARGV[++\$i] if \$ARGV[\$i] eq '--aligned';
	\$reads1 ||= \$ARGV[\$i + 1] if \$ARGV[\$i] eq '--reads';
	\$paired = 1 if \$ARGV[\$i] eq '--out2';
	\$blast = 1 if \$ARGV[\$i] eq '--blast';
}
if (\$blast) { # BLAST -m8 as SortMeRNA writes it: 12 columns, e-value and bit score last
	open my \$q, '<', \$reads1 or die; my (\$id) = <\$q> =~ /^>(\\S+)/;
	open my \$b, '>', "\$aligned.blast" or die;
	print {\$b} "\$id\\tref\\t99\\t4\\t0\\t0\\t1\\t4\\t1\\t4\\t1e-05\\t8\\n"; close \$b;
	exit 0;
}
open my \$log, '>>', '$sortLog' or die; print {\$log} "\$aligned\\n"; close \$log;
my \$reads = "\\\@new\\nACGT\\n+\\nIIII\\n";
gzip(\\\$reads => \$_) for \$paired ? ("\${aligned}_fwd.fq.gz", "\${aligned}_rev.fq.gz") : ("\$aligned.fq.gz");
});
my $vsearch = tool('vsearch', qq{#!/usr/bin/env perl
my (\$out, \$query) = ('', '');
for (my \$i = 0; \$i < \@ARGV; \$i++) {
	\$out = \$ARGV[\$i + 1] if \$ARGV[\$i] eq '--userout' || \$ARGV[\$i] eq '--output';
	\$query = \$ARGV[\$i + 1] if \$ARGV[\$i] eq '--usearch_global';
}
if (\$query ne '') {
	open my \$q, '<', \$query or die; my \$first = <\$q>;
	open my \$log, '>>', '$queryLog' or die; print {\$log} \$first; close \$log;
}
open my \$o, '>', \$out or die; close \$o;
});
my $lcaInputLog = "$tmp/lca_inputs.log";
my $lca = tool('LCA', qq{#!/usr/bin/env perl
my (\$out, \$in) = ('', '');
for (my \$i = 0; \$i < \@ARGV; \$i++) {
	\$out = \$ARGV[\$i + 1] if \$ARGV[\$i] eq '-o';
	\$in = \$ARGV[\$i + 1] if \$ARGV[\$i] eq '-i';
}
open my \$log, '>>', '$lcaInputLog' or die;
for my \$file (split /,/, \$in) { open my \$f, '<', \$file or die; print {\$log} <\$f>; }
close \$log;
open my \$o, '>', \$out or die; print {\$o} "read\\tdomain\\nr\\tBacteria\\n"; close \$o;
});

my $db = "$tmp/db";
write_file("$db/$_", '') for qw(sort-ssu.fasta sort-lsu.fasta ssu.fasta ssu.tax lsu.fasta lsu.tax);
my $config = "$tmp/matafiler.cfg";
write_file($config, join("\n",
	"MFLRDir\t$root", "BINDir\t$tmp/bin", "DBDir\t$db", "Rpath\tR",
	"CONDcmd\tconda", "CONDA\tshell hook", "Rscript\tRscript",
	"sortmerna\t$sortmerna", "SSUdbFAsrt\t$db/sort-ssu.fasta", "LSUdbFAsrt\t$db/sort-lsu.fasta",
	"SSUidx\t$db", "LSUidx\t$db", "vsearch\t$vsearch", "flash\t/bin/false", "LCA\t$lca",
	"SSUdbFA\t$db/ssu.fasta", "SSUtax\t$db/ssu.tax",
	"LSUdbFA\t$db/lsu.fasta", "LSUtax\t$db/lsu.tax",
)."\n");
my $fastq = "\@r\nACGT\n+\nIIII\n";
write_file("$tmp/in/R1.fq", $fastq);
write_file("$tmp/in/R2.fq", $fastq);
write_file("$tmp/in/S.fq", $fastq);

# ----------------------------------------------------- extraction contract
# Before catchLSUSSU 0.6 a pair-only marker published no singleton file. MATAF4
# requires all three files, catchLSUSSU only the roles of the current input:
# the job exited "complete" and MATAF4 resubmitted it on every pass.
{
	my $ribo = "$tmp/legacy/ribos";
	for my $tag (qw(SSU LSU)) {
		write_file("$ribo/${tag}_pull.sto", '');
		write_gzip("$ribo/reads_$tag.r$_.fq.gz", $fastq) for 1, 2;
	}
	write_file("$ribo/reads_LSU.fq.gz", ''); # zero-byte placeholder of older versions
	ok(!ribosome_completion_evidence(ribo_root => $ribo, requested => 1)->{profile_complete},
		'MATAF4 rejects a profile without the singleton file');
	my ($status, $output, $errors) = run_script('catchLSUSSU.pl',
		'-R1', "$tmp/in/R1.fq", '-R2', "$tmp/in/R2.fq", '-RS', '-1',
		'-alignDir', $ribo, '-tmpDir', "$tmp/scratch", '-smplID', 'legacy',
		'-cores', 1, '-config', $config);
	is($status, 0, 'catchLSUSSU completes the legacy pair-only profile') or diag($errors);
	is(lines($sortLog), 0, 'without repeating the SortMeRNA search');
	for my $tag (qw(SSU LSU)) {
		my $single = "$ribo/reads_$tag.fq.gz";
		ok(-s $single && read_gzip($single) eq '',
			"$tag singleton role is an empty gzip container");
	}
	ok(ribosome_completion_evidence(ribo_root => $ribo, requested => 1)->{profile_complete},
		'MATAF4 now accepts the profile, so the sample can close');
}
{
	# A role the input does have is never backfilled: the marker is rerun.
	my $ribo = "$tmp/partial/ribos";
	for my $tag (qw(SSU LSU)) {
		write_file("$ribo/${tag}_pull.sto", '');
		write_gzip("$ribo/reads_$tag.r$_.fq.gz", $fastq) for 1, 2;
		write_gzip("$ribo/reads_$tag.fq.gz", $fastq) if $tag eq 'LSU';
	}
	unlink $sortLog;
	my ($status, $output, $errors) = run_script('catchLSUSSU.pl',
		'-R1', "$tmp/in/R1.fq", '-R2', "$tmp/in/R2.fq", '-RS', "$tmp/in/S.fq",
		'-alignDir', $ribo, '-tmpDir', "$tmp/scratch", '-smplID', 'partial',
		'-cores', 1, '-config', $config);
	is($status, 0, 'catchLSUSSU reruns a marker missing an input role') or diag($errors);
	is(lines($sortLog), 2, 'only that marker is searched again (pairs and singletons)');
	like(read_gzip("$ribo/reads_SSU.fq.gz"), qr/^\@new/, 'its singleton reads are extracted');
}
{
	# Older lotus versions unpacked the published reads in place. A rerun
	# removes those plain files, so they cannot outlive the new extraction.
	my $ribo = "$tmp/plain/ribos";
	for my $tag (qw(SSU LSU)) {
		write_file("$ribo/${tag}_pull.sto", '');
		write_file("$ribo/reads_$tag.$_", "\@old\nACGT\n+\nIIII\n") for qw(r1.fq r2.fq fq);
	}
	my ($status, $output, $errors) = run_script('catchLSUSSU.pl',
		'-R1', "$tmp/in/R1.fq", '-R2', "$tmp/in/R2.fq", '-RS', "$tmp/in/S.fq",
		'-alignDir', $ribo, '-tmpDir', "$tmp/scratch", '-smplID', 'plain',
		'-cores', 1, '-config', $config);
	is($status, 0, 'catchLSUSSU re-extracts a profile published as plain FASTQ') or diag($errors);
	ok(!grep({ -e "$ribo/reads_SSU.$_" || -e "$ribo/reads_LSU.$_" } qw(r1.fq r2.fq fq)),
		'the legacy plain read files are removed');
}

# ----------------------------------------------------------- LCA assignment
{
	# A plain file next to the published .gz is stale: the .gz is assigned.
	my $ribo = "$tmp/stale/ribos";
	for my $tag (qw(SSU LSU)) {
		write_file("$ribo/${tag}_pull.sto", '');
		write_file("$ribo/reads_$tag.fq", "\@old\nACGT\n+\nIIII\n");
		write_gzip("$ribo/reads_$tag.fq.gz", "\@new\nACGT\n+\nIIII\n");
		write_gzip("$ribo/reads_$tag.r$_.fq.gz", '') for 1, 2;
	}
	my ($status, $output, $errors) = run_script('lotus_LCA_blast3.pl',
		'-dir', $ribo, '-DBdir', $db, '-smplID', 'stale', '-pairedRds', 0,
		'-cores', 1, '-simMode', 4, '-config', $config);
	is($status, 0, 'LCA assigns a sample with leftover plain reads') or diag($errors);
	is(read_file($queryLog), ">new\n>new\n", 'both markers use the published reads');
}
{
	# Without its extraction checkpoint, missing reads looked like a sample
	# without ribosomal reads: empty hierarchies were checkpointed as complete.
	my $ribo = "$tmp/unextracted/ribos";
	make_path($ribo);
	my ($status, $output, $errors) = run_script('lotus_LCA_blast3.pl',
		'-dir', $ribo, '-DBdir', $db, '-smplID', 'unextracted', '-pairedRds', 2,
		'-cores', 1, '-simMode', 4, '-config', $config);
	isnt($status, 0, 'LCA refuses a marker that was never extracted');
	like($errors, qr/SSU extraction checkpoint .* is missing/, 'and names the missing checkpoint');
	ok(!-e "$ribo/ltsLCA/Assigned.sto" && !-e "$ribo/ltsLCA/SSU_ass.sto",
		'no assignment checkpoint is written');
}

{
	# SortMeRNA as similarity search (-simMode 3): LCA reads column 11 as the
	# query length for -cover; SortMeRNA's table has the e-value there.
	my $ribo = "$tmp/smrsearch/ribos";
	for my $tag (qw(SSU LSU)) {
		write_file("$ribo/${tag}_pull.sto", '');
		write_gzip("$ribo/reads_$tag.fq.gz", "\@hit7 extra\nACGTA\n+\nIIIII\n");
		write_gzip("$ribo/reads_$tag.r$_.fq.gz", '') for 1, 2;
	}
	unlink $lcaInputLog;
	my ($status, $output, $errors) = run_script('lotus_LCA_blast3.pl',
		'-dir', $ribo, '-DBdir', $db, '-smplID', 'smr', '-pairedRds', 0,
		'-cores', 1, '-simMode', 3, '-config', $config);
	is($status, 0, 'LCA runs with SortMeRNA as similarity search') or diag($errors);
	is(read_file($lcaInputLog), "hit7\tref\t99\t4\t0\t0\t1\t4\t1\t4\t5\n" x 2,
		'and receives the query length, not the e-value, in column 11');
}

# ------------------------------------------------ LCA reference copy (detectRibo)
{
	my $source = "$tmp/ref";
	my $copy = "$tmp/LCADB";
	write_file("$source/lsu.fasta", ">a\nACGT\n");
	write_file("$source/lsu.fasta.lba.gz", "index-bytes");
	write_file("$source/lsu.tax", "a\tBacteria\n");
	my %reference = (fasta => "$source/lsu.fasta", taxonomy => "$source/lsu.tax", directory => $copy);
	ok(!lca_reference_copy_current(%reference), 'a missing copy is not current');
	write_file("$copy/$_", read_file("$source/$_")) for qw(lsu.fasta lsu.fasta.lba.gz lsu.tax);
	ok(lca_reference_copy_current(%reference), 'a complete copy is current');
	write_file("$copy/lsu.fasta.lba.gz", "index");
	ok(!lca_reference_copy_current(%reference), 'a truncated index copy is copied again');
	write_file("$copy/lsu.fasta.lba.gz", "index-bytes");
	unlink "$copy/lsu.tax";
	ok(!lca_reference_copy_current(%reference), 'a copy without its taxonomy is copied again');
	write_file("$copy/lsu.tax", read_file("$source/lsu.tax"));
	unlink "$source/lsu.fasta.lba.gz";
	ok(!lca_reference_copy_current(%reference), 'an unbuilt source index is built and copied');
}

# ---------------------------------------------------- hierarchies (RiboMeta)
{
	my $sample = "$tmp/samples/S1";
	write_file("$sample/ribos/ltsLCA/SSUriboRun_bl.hiera.txt", "read\tdomain\nr\tBacteria\n");
	write_file("$sample/ribos/ltsLCA/LSUriboRun_bl.hiera.txt", "read\tdomain\n");
	write_file("$sample/ribos/ltsLCA/LSUriboRun_bl.hiera.txt.gz", "partial"); # interrupted gzip
	is_deeply(compress_ribosome_hierarchies(sample_root => $sample), [],
		'both hierarchies are compressed in the sample directory');
	ok(!-e "$sample/ribos/ltsLCA/SSUriboRun_bl.hiera.txt", 'no plain copy is left');
	is(read_gzip("$sample/ribos/ltsLCA/LSUriboRun_bl.hiera.txt.gz"), "read\tdomain\n",
		'a leftover partial .gz is replaced');
	unlink "$sample/ribos/ltsLCA/LSUriboRun_bl.hiera.txt.gz";
	is_deeply(compress_ribosome_hierarchies(sample_root => $sample), ['LSU'],
		'a missing hierarchy is reported');
}

# ------------------------------------------------- merge from a sample list
{
	write_gzip("$tmp/list/one/A.gz", "read\tdomain\nr\tBacteria\n");
	write_file("$tmp/list/two/B.txt", "read\tdomain\nr\tArchaea\n");
	write_file("$tmp/list/samples.tsv", "A.SSU\t$tmp/list/one/A.gz\nB.SSU\t$tmp/list/two/B.txt\n");
	my ($status, $output, $errors) = run_script('miTagTaxTable.pl', 'domain',
		"$tmp/list/merged", "$tmp/list/samples.tsv");
	is($status, 0, 'the merge reads a sample list') or diag($errors);
	is(read_gzip("$tmp/list/merged.domain.txt.gz"), "domain\tA.SSU\tB.SSU\nArchaea\t0\t1\nBacteria\t1\t0\n",
		'with the listed column names, from any folder');
	write_file("$tmp/list/samples.tsv", "A.SSU\t$tmp/list/one/A.gz\nC.SSU\t$tmp/list/gone.txt\n");
	($status, $output, $errors) = run_script('miTagTaxTable.pl', 'domain',
		"$tmp/list/merged", "$tmp/list/samples.tsv");
	isnt($status, 0, 'a listed sample without hierarchy stops the merge');
	like($errors, qr/Hierarchy of C\.SSU does not exist/, 'and is named');
}
SKIP: {
	# Directory input still names hierarchy links it cannot read.
	my $probe = "$tmp/symlink-probe";
	write_file("$probe.target", '');
	my $linked = eval { symlink("$probe.target", $probe) } && -l $probe;
	unlink "$probe.target";
	skip 'POSIX symbolic links are not available', 3
		unless $linked && -l $probe && !-e $probe;
	write_file("$tmp/merge/A.hiera.txt", "read\tdomain\nr\tBacteria\n");
	write_gzip("$tmp/gone.hiera.txt.gz", "read\tdomain\n");
	symlink("$tmp/gone.hiera.txt.gz", "$tmp/merge/B.SSU.hiera.txt.gz") or die $!;
	unlink "$tmp/gone.hiera.txt.gz";
	my ($status, $output, $errors) = run_script('miTagTaxTable.pl', 'domain', "$tmp/merged", "$tmp/merge");
	is($status, 0, 'a directory merge still runs');
	like($errors, qr/Skipping 1 hierarchy link.*B\.SSU\.hiera\.txt\.gz/, 'and names the dangling link');
	like(read_gzip("$tmp/merged.domain.txt.gz"), qr/^domain\tA\n/, 'the readable sample is merged');
}

# ------------------------------------------- riboSummary over the whole map
# The real subs from MATAF4.pl, with qsubSystem capturing the merge jobs,
# which are then run here.
our (%MFopt, %MFconfig, %map, %preDIRs, $controllerBaseOut, $ignoredSamplesHR, $QSBoptHR);
our @samples;
my @mergeJobs;
{
	no warnings qw(redefine once);
	*main::qsubSystem = sub { push @mergeJobs, {script => $_[0], cmd => $_[1]}; return ('job'.@mergeJobs, '') };
	*main::getProgPaths = sub {
		die "unexpected getProgPaths($_[0])\n" unless $_[0] eq 'mrgMiTag_scr';
		return "$^X -I$root $scripts/miTagTaxTable.pl";
	};
}
my $main = read_file("$root/MATAF4.pl");
my %subBody;
for my $name (qw(_shell_quote _shell_command riboCohortDir riboSampleExcluded riboSummary)) {
	($subBody{$name}) = $main =~ /(^sub \Q$name\E\b[^\n]*\{.*?^\})/ms;
	die "Missing sub $name in MATAF4.pl\n" unless defined $subBody{$name};
	eval $subBody{$name};
	die "$name: $@" if $@;
}
sub complete_ribo {
	my ($sampleRoot, $taxon) = @_;
	my $ribo = "$sampleRoot/ribos";
	for my $tag (qw(SSU LSU)) {
		write_file("$ribo/${tag}_pull.sto", '');
		write_gzip("$ribo/reads_$tag.$_", $fastq) for qw(r1.fq.gz r2.fq.gz fq.gz);
		write_file("$ribo/ltsLCA/${tag}_ass.sto", '');
		write_gzip("$ribo/ltsLCA/${tag}riboRun_bl.hiera.txt.gz",
			join("\t", qw(read domain phylum class order family genus species hit2db))."\n"
			.join("\t", 'r', $taxon, ('?') x 7)."\n");
	}
	write_file("$ribo/ltsLCA/Assigned.sto", '');
}
sub summary_output {
	my $printed = '';
	open my $capture, '>', \$printed or die;
	my $previous = select $capture;
	riboSummary();
	select $previous;
	return $printed;
}
sub run_merge_jobs {
	for my $job (@mergeJobs) {
		write_file("$tmp/mergejob.sh", "set -eo pipefail\n$job->{cmd}");
		system('bash', "$tmp/mergejob.sh") == 0 or die "merge job failed: $job->{cmd}\n";
	}
	@mergeJobs = ();
}
{
	%MFopt = (DoRibofind => 1);
	%MFconfig = (silent => 0);
	%preDIRs = (dir2RiboF => 'pseudoGC/Phylo/RiboFind/');
	$controllerBaseOut = "$tmp/runA/";
	$ignoredSamplesHR = {D => 1};
	# Two #OutPath folders (or map files): A, C, D, E in runA, B in runB.
	%map = map { $_->[0] => {SmplID => $_->[0], wrdir => "$tmp/$_->[1]/$_->[0]/"} }
		([A => 'runA'], [B => 'runB'], [C => 'runA'], [D => 'runA'], [E => 'runA']);
	@samples = qw(A B C D E);
	complete_ribo("$tmp/runA/A", 'Bacteria');
	complete_ribo("$tmp/runB/B", 'Archaea');
	write_file("$tmp/runA/C/SMPL.empty", '');               # no reads
	write_sample_completion(                                # skipped, with its outcome kept
		root => "$tmp/runA/E/", sample => 'E', request_signature => 'e' x 64,
		present_assembly => 0, components => {},
		outcome => {status => 'skipped_too_small', input_size_mb => {primary => 0, supplementary => 0, total => 0}},
		metagstats => {DIR => 'E', values => {}, families => {}},
	);
	my $cohort = "$tmp/runA/pseudoGC/Phylo/RiboFind";

	summary_output();
	is(scalar(@mergeJobs), 2, 'the SSU and LSU tables are merged');
	is(read_file("$cohort/SSU.miTag.samples.tsv"),
		"A.SSU\t$tmp/runA/A/ribos/ltsLCA/SSUriboRun_bl.hiera.txt.gz\n"
		."B.SSU\t$tmp/runB/B/ribos/ltsLCA/SSUriboRun_bl.hiera.txt.gz\n",
		'from exactly the profiled map samples, across output folders');
	run_merge_jobs();
	is(read_gzip("$cohort/SSU.miTag.domain.txt.gz"), "domain\tA.SSU\tB.SSU\nArchaea\t0\t1\nBacteria\t1\t0\n",
		'into one table next to the first sample\'s output');
	ok(-s "$cohort/SSU.miTag.sto" && -s "$cohort/LSU.miTag.sto", 'each merge records its input');
	like(summary_output(), qr/SSU table is current for 2 samples/, 'an unchanged map is not merged again');
	is(scalar(@mergeJobs), 0, 'no merge job is submitted');

	# A run over another map (or a sample dropped from it) merges its own samples.
	@samples = qw(A C D E);
	summary_output();
	is(scalar(@mergeJobs), 2, 'a changed map is merged again');
	run_merge_jobs();
	is(read_gzip("$cohort/SSU.miTag.domain.txt.gz"), "domain\tA.SSU\nBacteria\t1\n",
		'and the table holds only its samples');

	# A re-assigned sample changes its hierarchy and the merge input with it.
	utime(time + 10, time + 10, "$tmp/runA/A/ribos/ltsLCA/SSUriboRun_bl.hiera.txt.gz");
	summary_output();
	is(scalar(@mergeJobs), 1, 'a changed hierarchy is merged again');
	@mergeJobs = ();

	# Every sample of the map must be profiled first, -from/-to or not.
	push @samples, 'F';
	$map{F} = {SmplID => 'F', wrdir => "$tmp/runB/F/"};
	like(summary_output(), qr/1 of 4 map samples have no complete RiboFind profile \(first: F\)/,
		'an unprofiled map sample holds the merge back (ignored samples do not count)');
	is(scalar(@mergeJobs), 0, 'no partial table is merged');
	ok($subBody{riboSummary} !~ /selected(?:From|To)|runOptions\{(?:from|to)\}/,
		'riboSummary does not depend on the -from/-to range');
}

# ----------------------------------------------------------- options
{
	my %options = (doRiboAssembl => 1);
	ok(!eval { normalise_ribosome_request(\%options); 1 }, '-riobsomalAssembly 1 stops at startup');
	like($@, qr/no longer supported/, 'with an explicit reason');
	%options = (doRiboAssembl => 0, RedoRiboAssign => 1);
	ok(normalise_ribosome_request(\%options), '-riobsomalAssembly 0 is accepted');
}

# ----------------------------------------------------------- MATAF4 wiring
{
	my ($detectRibo) = $main =~ /(^sub detectRibo\(\)\{.*?^\})/ms;
	ok($main =~ /\|\| \$calcKraken \|\| \$calcRibofind \|\| \$calcMOTU2 /
		&& $main =~ /my \$dowstreamAnalysisFlag = \(\$stagedReadsAnalysisFlag \|\| \$calcProtal \|\| \$riboAssignOnly\)/,
		'a pending LCA step alone neither stages nor cleans reads, but keeps the sample open');
	ok($main =~ /if \(\$riboAssignOnly\) \{\s*my \$riboLCAJob = detectRibo\(\$nodeSpTmpD\."ITS\/",\$curOutDir\."ribos\/","",.*?deferLoopProducerWave\(\s*'input staging'/s,
		'the LCA step is submitted before input staging, without a read dependency');
	ok($main =~ /if \(\$calcRibofind\)\{ #extraction, then assignment/,
		'the read-dependent submission is left to extraction');
	my $readSetCalls = () = $detectRibo =~ /sampleReadSet\(/g;
	ok($detectRibo =~ /if \(\$calcRiboFind\) \{\s*#[^\n]*\n\s*my \$cleanSeqSetHR = sampleReadSet\(/
		&& $readSetCalls == 1,
		'detectRibo reads the cleaned-read set only to extract');
	ok($detectRibo =~ /-pairedRds 2 / && $detectRibo =~ /normalise_job_dependencies\(\$jobName, \$DBdep\)/,
		'the LCA job takes the read layout from the extraction and waits for its database copy');
	ok($detectRibo =~ /if \(!exists\(\$MFopt\{globalRiboDependence\}->\{\$DBrna2\}\) \)/,
		'the LCA database is prepared once per output folder');
	ok($main =~ /my \@DBn = \("LSUdbFA","LSUtax","SSUdbFA","SSUtax","PR2dbFA","PR2tax"\)/,
		'detectRibo copies a configured PR2 reference for the LCA jobs');
	ok($main =~ /next if \(lca_reference_copy_current\(fasta => \$DB, taxonomy => \$taxDB, directory => \$DBrna2\)\);/,
		'detectRibo checks the LCA reference copy by size');
	ok($main !~ /riboFind(?:Fail|Compl)Cnts/, 'no per-pass RiboFind counters remain');
	ok($main !~ qr/my \$mrgCmd2? = [^;]*unless/, 'riboSummary has no "my ... unless" declarations');
}

done_testing;
