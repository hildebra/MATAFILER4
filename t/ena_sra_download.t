use strict;
use warnings;

use Digest::MD5;
use File::Copy qw(copy);
use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use IO::Compress::Gzip qw(gzip $GzipError);
use JSON::PP;
use FindBin qw($Bin);
use Test::More;

my $root = tempdir(CLEANUP => 1);
my $fixtures = File::Spec->catdir($root, 'fixtures');
my $tools = File::Spec->catdir($root, 'tools');
make_path($fixtures, $tools);
my $downloader = File::Spec->catfile(
	$Bin, '..', 'secScripts', 'fileManage', 'ENASRAdl.pl',
);
my $mataf = File::Spec->catfile($Bin, '..', 'MATAF4.pl');
open my $mataf_fh, '<', $mataf or die "Cannot read $mataf: $!";
my $mataf_source = do { local $/; <$mataf_fh> };
close $mataf_fh;
like(
	$mataf_source,
	qr/my \$downloadConfig = \$MFconfig\{configFile\};.*?File::Spec->rel2abs\(\$downloadConfig\).*?\('--config', \$downloadConfig\).*?\@configArguments/s,
	'the download job receives the selected config as a scheduler-safe absolute path');

sub write_file {
	my ($path, $content) = @_;
	open my $fh, '>', $path or die "Cannot write $path: $!";
	print {$fh} $content;
	close $fh;
}

sub write_tool {
	my ($name, $content) = @_;
	my $path = File::Spec->catfile($tools, $name);
	write_file($path, $content);
	chmod 0755, $path or die "Cannot chmod $path: $!";
	return $path;
}

sub read_json {
	my ($path) = @_;
	open my $fh, '<', $path or die "Cannot read $path: $!";
	local $/;
	my $decoded = JSON::PP->new->decode(<$fh>);
	close $fh;
	return $decoded;
}

sub md5_file {
	my ($path) = @_;
	open my $fh, '<', $path or die $!;
	binmode $fh;
	my $md5 = Digest::MD5->new->addfile($fh)->hexdigest;
	close $fh;
	return $md5;
}

my $fake_wget = write_tool('fake-wget', <<'PERL');
#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename qw(basename);
use File::Copy qw(copy);
exit 73 if $ENV{FAKE_WGET_FAIL};
my ($output, $url);
for (my $i = 0; $i < @ARGV; $i++) {
	$output = $ARGV[++$i] if $ARGV[$i] eq '--output-document';
}
$url = $ARGV[-1];
my $source;
if ($url =~ /ena\/portal\/api\/filereport/) {
	$source = $ENV{FAKE_ENA_REPORT};
} elsif ($url =~ /trace\.ncbi\.nlm\.nih\.gov/) {
	$source = $ENV{FAKE_SRA_RUNINFO};
} else {
	$source = "$ENV{FAKE_ARCHIVE_DIR}/".basename($url);
}
die "No fixture for $url\n" unless defined($source) && -f $source;
copy($source, $output) or die "copy $source -> $output failed: $!\n";
PERL

my %ena_fastq = (
	'ERR900001_1.fastq.gz' => "\@r1/1\nACGT\n+\nIIII\n",
	'ERR900001_2.fastq.gz' => "\@r1/2\nTGCA\n+\nIIII\n",
);
for my $name (keys %ena_fastq) {
	my $path = File::Spec->catfile($fixtures, $name);
	gzip \$ena_fastq{$name} => $path
		or die "Cannot gzip fixture: $GzipError";
}
my @ena_names = sort keys %ena_fastq;
my @ena_md5 = map { md5_file(File::Spec->catfile($fixtures, $_)) } @ena_names;
my @ena_bytes = map { -s File::Spec->catfile($fixtures, $_) } @ena_names;
my $ena_report = File::Spec->catfile($fixtures, 'ena.tsv');
write_file(
	$ena_report,
	join("\t", qw(run_accession instrument_platform library_layout fastq_ftp fastq_md5 fastq_bytes))."\n"
	.join("\t",
		'ERR900001', 'ILLUMINA', 'PAIRED',
		join(';', map { "ftp.sra.ebi.ac.uk/vol1/fastq/$_" } @ena_names),
		join(';', @ena_md5), join(';', @ena_bytes),
	)."\n",
);

my $ena_output = File::Spec->catdir($root, 'ena-download');
my $ena_metadata = File::Spec->catfile($root, 'ena-metadata.json');
{
	local $ENV{FAKE_ENA_REPORT} = $ena_report;
	local $ENV{FAKE_ARCHIVE_DIR} = $fixtures;
	is(system(
		$^X, $downloader,
		'--provider', 'ena', '--accession', 'ERR900001',
		'--output-dir', $ena_output, '--metadata-out', $ena_metadata,
		'--wget', $fake_wget, '--threads', 2,
	), 0, 'ENA download and validation succeeds');
}
my $ena = read_json($ena_metadata);
is($ena->{seqtech}, 'ill', 'ENA ILLUMINA platform maps to MATAFILER ill');
is($ena->{provider}, 'ena', 'ENA provenance records the selected provider');
is_deeply([map { $_->{role} } @{$ena->{files}}], [qw(r1 r2)],
	'ENA paired files receive explicit manifest roles');
is(scalar(grep { !-s File::Spec->catfile($ena_output, $_->{filename}) } @{$ena->{files}}),
	0, 'all ENA manifest files are present and non-empty');

{
	local $ENV{FAKE_WGET_FAIL} = 1;
	is(system(
		$^X, $downloader,
		'--provider', 'ena', '--accession', 'ERR900001',
		'--output-dir', $ena_output, '--metadata-out', $ena_metadata,
		'--wget', $fake_wget,
	), 0, 'a fully validated ENA download is reused without network access');
}

my $corrupt = File::Spec->catfile($ena_output, $ena_names[0]);
write_file($corrupt, "corrupt\n");
{
	local $ENV{FAKE_ENA_REPORT} = $ena_report;
	local $ENV{FAKE_ARCHIVE_DIR} = $fixtures;
	is(system(
		$^X, $downloader,
		'--provider', 'ena', '--accession', 'ERR900001',
		'--output-dir', $ena_output, '--metadata-out', $ena_metadata,
		'--wget', $fake_wget,
	), 0, 'a corrupt ENA file is detected and downloaded again');
}
is(md5_file($corrupt), md5_file(File::Spec->catfile($fixtures, $ena_names[0])),
	'recovered ENA file matches the archive checksum');

my $sra_runinfo = File::Spec->catfile($fixtures, 'sra.csv');
write_file($sra_runinfo,
	"Run,Platform,LibraryLayout\nSRR900001,OXFORD_NANOPORE,SINGLE\n");
my $fake_prefetch = write_tool('prefetch', <<'PERL');
#!/usr/bin/env perl
use strict;
use warnings;
use File::Path qw(make_path);
my ($output, $run);
for (my $i = 0; $i < @ARGV; $i++) {
	$output = $ARGV[++$i] if $ARGV[$i] eq '--output-directory';
}
$run = $ARGV[-1];
make_path($output);
die "prefetch resume cache was discarded\n"
	unless -e "$output/.partial";
open my $fh, '>', "$output/$run.sra" or die $!;
print {$fh} "validated sra object\n";
close $fh;
PERL
my $fake_validator = write_tool('vdb-validate', <<'PERL');
#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename qw(basename);
my $directory = $ARGV[0];
my $run = basename($directory);
exit((-s "$directory/$run.sra") ? 0 : 1);
PERL
my $fake_fasterq = write_tool('fasterq-dump', <<'PERL');
#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename qw(basename);
my ($output, $temporary, $threads, $source);
for (my $i = 0; $i < @ARGV; $i++) {
	if ($ARGV[$i] eq '--outdir') {
		$output = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-t') {
		$temporary = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-e') {
		$threads = $ARGV[++$i];
	}
}
$source = $ARGV[-1];
die "fasterq-dump did not receive an accession directory\n"
	unless -d $source && basename($source) =~ /^[SED]RR\d+$/;
die "fasterq-dump thread or temporary-directory options are invalid\n"
	unless defined($threads) && $threads > 0
		&& defined($temporary) && $temporary eq $output;
my $run = basename($source);
open my $fh, '>', "$output/$run.fastq" or die $!;
print {$fh} "\@ont1\nACGTA\n+\nIIIII\n";
close $fh;
PERL
my $fake_gzip = write_tool('gzip', <<'PERL');
#!/usr/bin/env perl
use strict;
use warnings;
use IO::Compress::Gzip qw(gzip $GzipError);
my $input = $ARGV[-1];
gzip $input => "$input.gz" or die $GzipError;
unlink $input or die $!;
PERL
my $tool_config = File::Spec->catfile($root, 'downloader-config.txt');
write_file(
	$tool_config,
	join("\n",
		"MFLRDir\t$root",
		"BINDir\t$root/bin",
		"DBDir\t$root/db",
		"MGSTKDir\t$root/mgs",
		"CONDcmd\tmicromamba",
		'CONDA' . "\t" . 'eval "$(micromamba shell hook -s bash)"',
		"CONDAbaseEnv\tMF4",
		"Rscript\tRscript",
		"Rpath\t[MFLRDir]/R",
		"wget\t$fake_wget",
		"prefetch\t$fake_prefetch",
		"fasterq-dump\t$fake_fasterq",
		"vdb-validate\t$fake_validator",
		"pigz\t$fake_gzip",
	), "\n",
);

my $sra_output = File::Spec->catdir($root, 'sra-download');
my $sra_metadata = File::Spec->catfile($root, 'sra-metadata.json');
my $sra_resume_cache = File::Spec->catdir(
	$sra_output, '.ncbi-cache', 'SRR900001',
);
make_path($sra_resume_cache);
write_file(File::Spec->catfile($sra_resume_cache, '.partial'), 'resume');
{
	local $ENV{FAKE_SRA_RUNINFO} = $sra_runinfo;
	is(system(
		$^X, $downloader,
		'--provider', 'sra', '--accession', 'SRR900001',
		'--output-dir', $sra_output, '--metadata-out', $sra_metadata,
		'--config', $tool_config,
	), 0, 'SRA tools resolve through getProgPaths and the selected config');
}
my $sra = read_json($sra_metadata);
is($sra->{seqtech}, 'ONT', 'NCBI OXFORD_NANOPORE platform maps to ONT');
is($sra->{files}[0]{role}, 'single', 'split-3 singleton output is explicit');
ok(-s File::Spec->catfile($sra_output, 'SRR900001.fastq.gz'),
	'validated SRA FASTQ is published as gzip');
ok(!-d File::Spec->catdir($sra_output, '.ncbi-cache')
	&& !-d File::Spec->catdir($sra_output, '.fasterq'),
	'SRA object and conversion intermediates are removed after validation');

done_testing();
