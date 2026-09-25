use strict;
use warnings;
use Test::More;
use File::Temp qw(tempdir);
use File::Path qw(make_path);
use File::Spec;
use FindBin qw($Bin);
use lib "$Bin/..", "$Bin/lib";
use MFTestConfig;
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::SampleCompletion qw(sample_completion_path invalidate_sample_completion);

# Regressions for the follow-up decisions of 25 September 2026 (see
# docs/audits/2026-09-25/report.md, "Decisions applied").

my $root = File::Spec->rel2abs("$Bin/..");
my $tmp = tempdir(CLEANUP => 1);
sub read_file { open my $fh, '<', $_[0] or die "$_[0]: $!"; local $/; return <$fh> // ''; }
sub write_file { my ($p, $t) = @_; open my $fh, '>', $p or die "$p: $!"; print {$fh} $t; close $fh; }
sub source_sub {
	my ($source, $name) = @_;
	my ($code) = $source =~ /(^sub \Q$name\E\b[^\{;]*\{.*?^\})/ms;
	die "Cannot isolate $name" unless defined $code;
	return $code;
}
my $main = read_file("$root/MATAF4.pl");

# ---- item 7: a shared assembly is never built from a subset of its members
our %AsGrps;
eval source_sub($main, 'deferGroupAssemblyForClosedMembers');
die $@ if $@;
my @members = map { my $d = "$tmp/run/$_/"; make_path($d); $d } qw(A B);
write_file(sample_completion_path($members[$_], (qw(A B))[$_]), "{}\n") for 0, 1;
%AsGrps = (grp => {ClosedCompleted => [@members]});
my $said = '';
{
	open my $out, '>', \$said or die; local *STDOUT = $out;
	is(deferGroupAssemblyForClosedMembers('grp'), 1, 'assembly start is deferred while fast-path members registered no reads');
}
ok(!-e sample_completion_path($members[0], 'A') && !-e sample_completion_path($members[1], 'B'),
	'the stale completion sentinels of those members are removed');
like($said, qr/2 member\(s\) were closed as complete/, 'the deferral is explained');
is(deferGroupAssemblyForClosedMembers('grp'), 0, 'once reopened, the next start proceeds');
%AsGrps = (grp => {});
is(deferGroupAssemblyForClosedMembers('grp'), 0, 'groups without fast-path members are not delayed');
like($main, qr/for grep \{ -d \$_ \} \@\{assembly_group_output_dirs\(\\%map, \$cAssGrp\)\};/,
	'a group rewrite reopens every member of the shared assembly');
like($main, qr/\$AsGrps->\{\$cAssGrp\}\{ClosedCompleted\}|ClosedCompleted/, 'fast-path members are recorded');
like(read_file("$root/Mods/GenoMetaAss.pm"), qr/\{ClosedCompleted\} = \[\];/, 'the record is reset every loop pass');

# ---- items 5 and 17: database locations follow the installer (DBDir)
my $dbDir = getProgPaths('DBDir');
like(getProgPaths('checkm2DB'), qr{^\Q$dbDir\E/*CM2/CheckM2_database/uniref100\.KO\.1\.dmnd$},
	'CheckM2 database is where installer.sh downloads it (DBDir/CM2)');
like(getProgPaths('hostileDB'), qr{^\Q$dbDir\E/*hostile/?$}, 'hostile index is under DBDir');
like(getProgPaths('PtostT5_Weights'), qr{^\Q$dbDir\E/*PtostT5_W$}, 'ProstT5 weights are under DBDir');
like(read_file("$root/helpers/install/installer.sh"), qr/CM2DB="\$DBdir\/CM2"/, 'installer still uses DBdir/CM2');

# ---- item 14: abundance checkpoints ignore the empty-sample list
my $mgs = read_file("$root/secScripts/MGS.pl");
unlike($mgs, qr/_checkpoint_valid\(\$ABmgsSton2?\)/, 'abundance stones are checked with the resume validator');
like($mgs, qr/delete \$parameters\{empty_samples\}\s*if \$stage eq 'mgs-abundance' \|\| \$stage eq 'marker-mgs-abundance';/,
	'and written without the empty-sample list');

# ---- item 15: marker matrices exist before geneCat writes the marker stone
my $e100 = read_file("$root/secScripts/GC/extrAllE100GC.pl");
unlike($e100, qr/fileGZe\("\$GCd\/Mattrix/, 'the misspelled skip test is gone');
like($e100, qr/qsubSystemJobAlive\(\\\@matrixJobs, \$QSBoptHR\)/, 'the script waits for its matrix-subset jobs');
like($e100, qr/Marker gene matrix jobs finished without publishing/, 'and verifies their outputs');

# ---- item 16: no per-gene k-mers -> clean skip instead of a crash
like(read_file("$root/secScripts/GC/kmerPerGene.pl"), qr/Skipping the catalogue k-mer table/, 'k-mer step skips without MATAF4 k-mers');
like(read_file("$root/secScripts/geneCat.pl"), qr/if \[ -e \$OutD\/\$primaryClusterFNA\.kmer \]; then \$pigzBin/,
	'geneCat compresses the k-mer table only when it was written');

# ---- item 18: -perlClusterMAGs uses the sample__bin convention
my $cluster = read_file("$root/secScripts/MGS/clusterMAGs.pl");
is(scalar(() = $cluster =~ /^\t+(?:my )?\$uniqMBid = "\$smplIDs\[-1\]__\$bin";/mg), 2, 'both MAG id sites use sample__bin');
unlike($cluster, qr/^\t+(?:my )?\$uniqMBid = "\$smplIDs\[-1\]\.\$bin";/m, 'no live sample.bin ids remain');

done_testing;
