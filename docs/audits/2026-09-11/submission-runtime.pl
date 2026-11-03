#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Spec;
use File::Temp qw(tempdir);
use File::Path qw(make_path);
use Time::HiRes qw(time);
use JSON::PP;
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use lib "$Bin/../../..", "$Bin/../../../t/lib";
use MFTestConfig;
use Mods::GenoMetaAss qw(coverage_derivative_paths);

# Compare the sources saved immediately before this dependency audit with the
# current checkout. No aligner, cluster, or production dataset is used.
my $before = shift @ARGV or die "usage: $0 BEFORE_SOURCE_DIRECTORY\n";
my $root = File::Spec->rel2abs("$Bin/../../..");
my $tmp = tempdir(CLEANUP => 1);
sub read_file { open my $f, '<', $_[0] or die $!; local $/; return <$f>; }
sub write_file { open my $f, '>', $_[0] or die $!; print {$f} $_[1]; close $f or die $!; }
sub median { my @s = sort {$a <=> $b} @_; return $s[int(@s/2)]; }
my ($curSmpl, $JNUM, $avx2Constr) = ('sample', 1, '');
my $pigzBin = Mods::IO_Tamoc_progs::getProgPaths('pigz');
my $logDir = "$tmp/log/";
my %HDDspace = (diamond => 1);
my %MFopt = (diaCores=>1,diaRunSensitive=>0,diaFrameshift=>0,diaEVal=>'1e-7',
    DiaPercID=>20,DiaMinAlignLen=>30,DiaMinFracQueryCov=>0,DiaRmRawHits=>0,
    globalDiamondDependence=>{TEST=>'TEST-1'},diamondMem=>1);
my %progStats;
my $QSBoptHR = {constraint=>[],tmpSpace=>0,General_Hosts=>[]};
my (%libraries, @jobs);
sub sampleReadSet { return {}; }
sub readLibrariesByScope { return []; }
sub getRdLibraries { return %libraries; }
sub prepDiamondDB { return ('ref','TEST',''); }
sub getProgPaths { return 'fixture-diamond'; }
sub qsubSystem { push @jobs, [@_]; return ($_[4], $_[1]); }

sub extract {
    my ($file, $name) = @_;
    my ($body) = read_file($file) =~ /(^sub \Q$name\E[^\n]*\{.*?^\})/ms;
    die "Missing $name in $file" unless $body;
    $body =~ s/^sub \Q$name\E(?:\(\))?/sub /;
    my $sub = eval $body; die $@ if $@; return $sub;
}
my %results = (description => 'Warm local filesystem and generated shell command microbenchmarks; not HPC wall time');
my @coverage = map { extract($_, 'coverage_derivatives_complete') }
    ("$before/GenoMetaAss.pm", "$root/Mods/GenoMetaAss.pm");
for my $suffix (qw(pergene percontig median.percontig)) {
    write_file("$tmp/coverage.$suffix", 'complete');
}
for my $i (0, 1) {
    my @timings;
    for (1..3) {
        my $start = time;
        $coverage[$i]->("$tmp/coverage.gz") for 1..30000;
        push @timings, time - $start;
    }
    $results{coverage}{$i ? 'after' : 'before'} = {
        checks => 30000, seconds => \@timings, median_seconds => median(@timings),
    };
}
{
    package AuditFileProbe;
    use overload '-X' => sub { ${$_[0][0]}++; return $_[0][1]; }, fallback => 1;
}
for my $i (0, 1) {
    my $probes = 0;
    no warnings 'redefine';
    local *coverage_derivative_paths = sub {
        return [map { bless [\$probes, $_ == 0 ? 1 : 0], 'AuditFileProbe' } 0..3];
    };
    die 'Unexpected coverage result' unless $coverage[$i]->('canonical');
    $results{coverage}{$i ? 'after' : 'before'}{file_tests_per_complete_check} = $probes;
}

my @search = map { extract($_, 'runDiamond') } ("$before/MATAF4.pl", "$root/MATAF4.pl");
my %input;
for my $key (0..3) {
    my $text = join '', map {
        my $query = $key == 1 ? "pair$_/2" : $key == 2 ? "pair$_/1" : "lib${key}_$_";
        "$query\tB\t90\t50\t0\t0\t1\t150\t1\t50\t1e-20\t100\n"
    } 1..50000;
    gzip(\$text => \$input{$key}) or die $GzipError;
}
for my $case ([unmerged => [0,1,2]], [mixed => [0,1,2,3]], [merged => [3]]) {
    my ($name, $keys) = @$case;
    %libraries = map { $_ => ["fixture$_.gz"] } @$keys;
    my @outputs;
    for my $i (0,1) {
        @jobs = ();
        my $out = "$tmp/out_${name}_$i/";
        my $scratch = "$tmp/scratch_${name}_$i";
        make_path($out);
        $search[$i]->($out, "$tmp/db/", $scratch, '', 'TEST');
        my $command = $jobs[0][1];
        # Keep all actual collection commands, excluding only the aligner.
        $command =~ s/^fixture-diamond[^\n]*\n//mg;
        write_file("$tmp/worker.sh", $command);
        my @timings;
        for (1..3) {
            make_path($scratch);
            write_file("$scratch/DiaAssignment.sub.TEST.$_.0.gz", $input{$_}) for @$keys;
            my $start = time;
            system('bash', '-e', '-o', 'pipefail', "$tmp/worker.sh") == 0 or die 'Collection failed';
            push @timings, time - $start;
        }
        my $hits = '';
        gunzip("${out}dia.TEST.blast.srt.gz" => \$hits, MultiStream=>1) or die $GunzipError;
        $hits =~ s/\tMF4:read_count=1(?=\n)//g;
        push @outputs, $hits;
        $results{diamond_collection}{$name}{$i ? 'after' : 'before'} = {
            seconds => \@timings, median_seconds => median(@timings),
            hit_rows => 50000 * @$keys,
            counted_intermediate_files => scalar(() = $command =~ / > [^\n]*\.counted\.gz/g),
        };
    }
    $results{diamond_collection}{$name}{same_hits_and_multiplicity} = $outputs[0] eq $outputs[1] ? JSON::PP::true : JSON::PP::false;
    die 'Changed hit content' unless $outputs[0] eq $outputs[1];
}
print JSON::PP->new->canonical->pretty->encode(\%results);
