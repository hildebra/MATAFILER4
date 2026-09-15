#!/usr/bin/env perl
use strict;
use warnings;

use Digest::MD5;
use File::Basename qw(basename dirname);
use File::Copy qw(move);
use File::Find ();
use File::Path qw(make_path remove_tree);
use File::Spec;
use FindBin qw($RealBin);
use Getopt::Long qw(GetOptions);
use IO::Uncompress::Gunzip qw($GunzipError);
use JSON::PP;
use POSIX qw(strftime);
use Text::ParseWords qw(parse_line);

use lib "$RealBin/../..";
use Mods::IO_Tamoc_progs qw(getProgPaths setConfigFile);

my %opt = (threads => 4, retries => 3);
GetOptions(
	'provider=s' => \$opt{provider}, 'accession=s' => \$opt{accession},
	'output-dir=s' => \$opt{output_dir}, 'metadata-out=s' => \$opt{metadata_out},
	'config=s' => \$opt{config},
	'threads=i' => \$opt{threads}, 'retries=i' => \$opt{retries},
	'wget=s' => \$opt{wget}, 'prefetch=s' => \$opt{prefetch},
	'fasterq-dump=s' => \$opt{fasterq_dump},
	'vdb-validate=s' => \$opt{vdb_validate},
	'compressor=s' => \$opt{compressor}, 'help' => \$opt{help},
) or die "Invalid arguments; use --help for usage\n";

setConfigFile($opt{config})
	if defined($opt{config}) && $opt{config} ne '';

if ($opt{help}) {
	print <<'USAGE';
Usage:
  ENASRAdl.pl --provider ena|sra --accession ID[,ID...]
               --output-dir DIR --metadata-out FILE [--threads N]
               [--config MATAFILER_CONFIG]

ENA downloads use the ENA Portal API and verify archive-supplied byte counts
and MD5 checksums. SRA downloads use NCBI prefetch, vdb-validate and
fasterq-dump. Every gzip FASTQ is decompressed and structurally validated
before an atomic completion manifest is published.
USAGE
	exit 0;
}
for my $required (qw(provider accession output_dir metadata_out)) {
	(my $flag = $required) =~ s/_/-/g;
	die "--$flag is required\n"
		unless defined($opt{$required}) && $opt{$required} ne '';
}
$opt{provider} = lc $opt{provider};
die "--provider must be 'ena' or 'sra'\n"
	unless $opt{provider} eq 'ena' || $opt{provider} eq 'sra';
die "--threads must be a positive integer\n"
	unless $opt{threads} =~ /^\d+$/ && $opt{threads} > 0;
die "--retries must be a positive integer\n"
	unless $opt{retries} =~ /^\d+$/ && $opt{retries} > 0;

sub normalise_accessions {
	my ($value) = @_;
	my @ids = grep { length } split /[\s,;]+/, $value;
	die "No archive accessions were supplied\n" unless @ids;
	my %seen;
	for my $id (@ids) {
		die "Unsafe or unsupported accession '$id'\n"
			unless $id =~ /^[A-Za-z0-9][A-Za-z0-9_.-]*$/;
	}
	return [grep { !$seen{$_}++ } @ids];
}

my $accessions = normalise_accessions($opt{accession});
my $accession_key = join(',', @{$accessions});
make_path($opt{output_dir}) unless -d $opt{output_dir};
die "Output directory is not writable: $opt{output_dir}\n"
	unless -d $opt{output_dir} && -w $opt{output_dir};

sub atomic_write {
	my ($path, $text) = @_;
	my $dir = dirname($path);
	make_path($dir) unless -d $dir;
	my $temporary = "$path.tmp.$$";
	open my $fh, '>', $temporary or die "Cannot write $temporary: $!\n";
	print {$fh} $text or die "Cannot write $temporary: $!\n";
	close $fh or die "Cannot close $temporary: $!\n";
	rename $temporary, $path or die "Cannot replace $path: $!\n";
}

sub read_json {
	my ($path) = @_;
	return unless -s $path;
	open my $fh, '<', $path or return;
	local $/;
	my $text = <$fh>;
	close $fh;
	my $decoded = eval { JSON::PP->new->decode($text) };
	return ref($decoded) eq 'HASH' ? $decoded : undef;
}

sub executable {
	my (@names) = grep { defined($_) && $_ ne '' } @_;
	for my $name (@names) {
		if (File::Spec->file_name_is_absolute($name) || $name =~ m{[\\/]}) {
			return $name if -f $name && -x $name;
			next;
		}
		for my $dir (File::Spec->path()) {
			my $candidate = File::Spec->catfile($dir, $name);
			return $candidate if -f $candidate && -x $candidate;
		}
	}
	return;
}

sub configured_executable {
	my ($override, $config_key, $description) = @_;
	my @configured;
	if (defined($override) && $override ne '') {
		@configured = ($override);
	} elsif (ref($config_key) eq 'ARRAY') {
		@configured = @{getProgPaths($config_key, 0)};
	} else {
		@configured = (getProgPaths($config_key));
	}
	my $resolved = executable(@configured);
	die "Configured $description executable is unavailable: "
		.join(', ', grep { defined($_) && $_ ne '' } @configured)."\n"
		unless defined $resolved;
	return $resolved;
}

sub run_checked {
	my (@command) = @_;
	my $status = system(@command);
	return if $status == 0;
	my $detail = $status == -1 ? "could not execute: $!"
		: ($status & 127) ? "signal ".($status & 127)
		: "exit code ".($status >> 8);
	die "Command failed ($detail): ".join(' ', @command)."\n";
}

sub md5_file {
	my ($path) = @_;
	open my $fh, '<', $path or die "Cannot read $path: $!\n";
	binmode $fh;
	my $digest = Digest::MD5->new->addfile($fh)->hexdigest;
	close $fh;
	return $digest;
}

sub validate_fastq {
	my ($path) = @_;
	my $fh;
	if ($path =~ /\.gz$/i) {
		$fh = IO::Uncompress::Gunzip->new($path, MultiStream => 1)
			or die "Cannot decompress $path: $GunzipError\n";
	} else {
		open $fh, '<', $path or die "Cannot read $path: $!\n";
	}
	my $records = 0;
	while (1) {
		my $header = <$fh>;
		last unless defined $header;
		my $sequence = <$fh>;
		my $plus = <$fh>;
		my $quality = <$fh>;
		die "Truncated FASTQ record in $path after $records complete record(s)\n"
			unless defined($sequence) && defined($plus) && defined($quality);
		die "Invalid FASTQ header in $path at record ".($records + 1)."\n"
			unless $header =~ /^\@/;
		die "Invalid FASTQ separator in $path at record ".($records + 1)."\n"
			unless $plus =~ /^\+/;
		chomp($sequence, $quality);
		$sequence =~ s/\r$//;
		$quality =~ s/\r$//;
		die "Sequence/quality length mismatch in $path at record ".($records + 1)."\n"
			unless length($sequence) == length($quality);
		$records++;
	}
	close $fh or die "Cannot finish reading $path: $!\n";
	die "FASTQ contains no records: $path\n" unless $records;
	return $records;
}

sub validate_record {
	my ($record, $directory) = @_;
	return 0 unless ref($record) eq 'HASH';
	my $filename = $record->{filename} // '';
	return 0 if $filename eq '' || basename($filename) ne $filename;
	my $path = File::Spec->catfile($directory, $filename);
	my @stat = stat($path);
	return 0 unless @stat && $stat[7] > 0;
	return 0 if defined($record->{bytes}) && $record->{bytes} ne ''
		&& $stat[7] != $record->{bytes};
	return 0 if ($record->{md5} || '') ne ''
		&& lc(md5_file($path)) ne lc($record->{md5});
	my $records = eval { validate_fastq($path) };
	return 0 unless defined($records) && $records > 0 && !$@;
	$record->{records} = $records;
	return 1;
}

sub metadata_reusable {
	my ($metadata) = @_;
	return 0 unless ref($metadata) eq 'HASH'
		&& ($metadata->{schema} || 0) == 1
		&& ($metadata->{provider} || '') eq $opt{provider}
		&& ($metadata->{accession} || '') eq $accession_key
		&& ref($metadata->{files}) eq 'ARRAY' && @{$metadata->{files}};
	for my $record (@{$metadata->{files}}) {
		return 0 unless validate_record($record, $opt{output_dir});
	}
	return 1;
}

my $internal_manifest = File::Spec->catfile($opt{output_dir}, 'download.json');
for my $candidate ($internal_manifest, $opt{metadata_out}) {
	my $metadata = read_json($candidate);
	next unless metadata_reusable($metadata);
	my $json = JSON::PP->new->canonical->pretty->encode($metadata);
	atomic_write($internal_manifest, $json);
	atomic_write($opt{metadata_out}, $json);
	print "Validated existing $opt{provider} download for $accession_key in $opt{output_dir}\n";
	exit 0;
}

my $wget;
sub wget_tool {
	$wget ||= configured_executable($opt{wget}, 'wget', 'wget');
	return $wget;
}

sub url_encode {
	my ($value) = @_;
	$value =~ s/([^A-Za-z0-9_.~-])/sprintf("%%%02X", ord($1))/ge;
	return $value;
}

sub fetch_url {
	my ($url, $destination) = @_;
	my $temporary = "$destination.part.$$";
	unlink $temporary if -e $temporary;
	run_checked(
		wget_tool(), '--quiet', '--https-only', '--tries', $opt{retries},
		'--timeout', '60', '--output-document', $temporary, $url,
	);
	die "Download produced an empty file: $url\n" unless -s $temporary;
	rename $temporary, $destination
		or die "Cannot publish download $destination: $!\n";
}

sub seqtech_for_platform {
	my ($platform) = @_;
	$platform = uc($platform // '');
	return 'PB'  if $platform =~ /PACBIO/;
	return 'ONT' if $platform =~ /OXFORD[_ ]NANOPORE/;
	return '454' if $platform =~ /(?:^|_)LS454$/ || $platform eq '454';
	return 'ill';
}

sub infer_seqtech {
	my (@platforms) = @_;
	my %technologies = map { seqtech_for_platform($_) => 1 }
		grep { defined($_) && $_ ne '' } @platforms;
	die "Archive metadata mixes incompatible sequencing technologies: "
		.join(', ', sort keys %technologies)."\n" if keys(%technologies) > 1;
	return (keys(%technologies))[0] || 'ill';
}

sub ena_role {
	my ($layout, $uri, $index, $count) = @_;
	return 'single' unless uc($layout || '') eq 'PAIRED';
	my $name = basename($uri);
	return 'r1' if $name =~ /(?:^|[_.])1(?:[_.]|\.f(?:ast)?q)/i;
	return 'r2' if $name =~ /(?:^|[_.])2(?:[_.]|\.f(?:ast)?q)/i;
	return 'r1' if $index == 0 && $count >= 2;
	return 'r2' if $index == 1 && $count >= 2;
	return 'single';
}

sub download_ena {
	my (@files, @platforms);
	my (%runs, %filenames);
	for my $requested (@{$accessions}) {
		my $report = File::Spec->catfile(
			$opt{output_dir}, ".ena-report-$requested-$$.tsv",
		);
		my $url = 'https://www.ebi.ac.uk/ena/portal/api/filereport?'
			.'accession='.url_encode($requested)
			.'&result=read_run&fields=run_accession,instrument_platform,library_layout,fastq_ftp,fastq_md5,fastq_bytes'
			.'&format=tsv&download=true';
		fetch_url($url, $report);
		open my $fh, '<', $report or die "Cannot read ENA report $report: $!\n";
		my $header = <$fh>;
		die "ENA returned no metadata for accession $requested\n" unless defined $header;
		chomp $header;
		$header =~ s/\r$//;
		my @headers = split /\t/, $header, -1;
		my %column;
		@column{@headers} = (0 .. $#headers);
		for my $required (qw(run_accession instrument_platform library_layout fastq_ftp fastq_md5 fastq_bytes)) {
			die "ENA report for $requested lacks '$required'\n"
				unless exists $column{$required};
		}
		my $rows = 0;
		while (my $line = <$fh>) {
			chomp $line;
			$line =~ s/\r$//;
			next if $line eq '';
			my @values = split /\t/, $line, -1;
			my $run = $values[$column{run_accession}] || '';
			next if $run eq '';
			$rows++;
			next if $runs{$run}++;
			my $platform = $values[$column{instrument_platform}] || '';
			my $layout = $values[$column{library_layout}] || '';
			my @uris = split /;/, ($values[$column{fastq_ftp}] || '');
			my @md5s = split /;/, ($values[$column{fastq_md5}] || '');
			my @bytes = split /;/, ($values[$column{fastq_bytes}] || '');
			die "ENA exposes no public FASTQ files for run $run (requested as $requested)\n"
				unless @uris && $uris[0] ne '';
			die "ENA metadata has inconsistent file/checksum/size counts for run $run\n"
				unless @uris == @md5s && @uris == @bytes;
			push @platforms, $platform if $platform ne '';
			for my $index (0 .. $#uris) {
				my $uri = $uris[$index];
				$uri =~ s{^ftp://}{};
				my $download_url = $uri =~ m{^https://} ? $uri : "https://$uri";
				die "ENA returned a non-EBI FASTQ location for $run: $download_url\n"
					unless $download_url =~ m{^https://(?:[^/]+\.)?(?:ebi\.ac\.uk|sra\.ebi\.ac\.uk)/}i;
				my $filename = basename($uri);
				die "Unsafe ENA FASTQ filename for $run: $filename\n"
					if $filename eq '' || $filename =~ m{[\\/]};
				die "ENA FASTQ filename collision: $filename\n" if $filenames{$filename}++;
				die "Invalid ENA byte count for $run/$filename\n"
					unless $bytes[$index] =~ /^\d+$/ && $bytes[$index] > 0;
				die "Invalid ENA MD5 for $run/$filename\n"
					unless $md5s[$index] =~ /^[0-9a-fA-F]{32}$/;
				my $record = {
					run => $run, filename => $filename,
					role => ena_role($layout, $uri, $index, scalar @uris),
					layout => uc($layout || 'UNKNOWN'),
					bytes => 0 + $bytes[$index], md5 => lc($md5s[$index]),
					url => $download_url,
				};
				my $destination = File::Spec->catfile($opt{output_dir}, $filename);
				if (!validate_record($record, $opt{output_dir})) {
					unlink $destination if -e $destination;
					my $last_error = '';
					for my $attempt (1 .. $opt{retries}) {
						my $ok = eval {
							fetch_url($download_url, $destination);
							die "ENA byte-count mismatch for $filename\n"
								unless -s $destination == $record->{bytes};
							die "ENA MD5 mismatch for $filename\n"
								unless lc(md5_file($destination)) eq $record->{md5};
							$record->{records} = validate_fastq($destination);
							1;
						};
						last if $ok;
						$last_error = $@ || 'unknown validation failure';
						unlink $destination if -e $destination;
					}
					die "Could not download and validate ENA file $filename: $last_error"
						unless -s $destination;
				}
				push @files, $record;
			}
		}
		close $fh;
		unlink $report;
		die "ENA returned no read runs for accession $requested\n" unless $rows;
	}
	return (\@files, \@platforms);
}

sub ncbi_run_info {
	my @records;
	for my $requested (@{$accessions}) {
		my $report = File::Spec->catfile(
			$opt{output_dir}, ".sra-runinfo-$requested-$$.csv",
		);
		my $url = 'https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?'
			.'save=efetch&db=sra&rettype=runinfo&term='.url_encode($requested);
		my $fetched = eval { fetch_url($url, $report); 1 };
		if (!$fetched) {
			unlink $report if -e $report;
			die "NCBI RunInfo lookup failed for non-run accession $requested: $@"
				unless $requested =~ /^[SED]RR\d+$/i;
			warn "NCBI RunInfo lookup failed for $requested; using the run accession and default sequencing technology\n";
			push @records, {run => uc($requested), platform => '', layout => 'UNKNOWN'};
			next;
		}
		open my $fh, '<', $report or die "Cannot read NCBI RunInfo $report: $!\n";
		my $header = <$fh>;
		die "NCBI returned no RunInfo for $requested\n" unless defined $header;
		chomp $header;
		$header =~ s/\r$//;
		my @headers = parse_line(',', 0, $header);
		my %column;
		@column{@headers} = (0 .. $#headers);
		die "NCBI RunInfo lacks a Run column for $requested\n" unless exists $column{Run};
		my $rows = 0;
		while (my $line = <$fh>) {
			chomp $line;
			$line =~ s/\r$//;
			next if $line eq '';
			my @values = parse_line(',', 0, $line);
			my $run = $values[$column{Run}] || '';
			next unless $run =~ /^[SED]RR\d+$/i;
			push @records, {
				run => uc($run),
				platform => exists($column{Platform}) ? ($values[$column{Platform}] || '') : '',
				layout => exists($column{LibraryLayout}) ? ($values[$column{LibraryLayout}] || 'UNKNOWN') : 'UNKNOWN',
			};
			$rows++;
		}
		close $fh;
		unlink $report;
		die "NCBI returned no run accessions for $requested\n" unless $rows;
	}
	my %seen;
	return [grep { !$seen{$_->{run}}++ } @records];
}

sub find_sra_path {
	my ($cache, $run) = @_;
	for my $candidate (
		File::Spec->catfile($cache, "$run.sra"),
		File::Spec->catfile($cache, $run),
		File::Spec->catfile($cache, $run, "$run.sra"),
	) {
		return $candidate if -f $candidate && -s $candidate;
	}
	my @matches;
	File::Find::find({
		no_chdir => 1,
		wanted => sub {
			return unless -f $_ && -s _;
			my $name = basename($_);
			push @matches, $File::Find::name if $name eq $run || $name eq "$run.sra";
		},
	}, $cache);
	die "prefetch completed but no SRA object was found for $run in $cache\n"
		unless @matches == 1;
	return $matches[0];
}

sub compress_fastq {
	my ($path, $compressor) = @_;
	my $name = basename($compressor);
	if ($name =~ /pigz/i) {
		run_checked($compressor, '-f', '-p', $opt{threads}, $path);
	} else {
		run_checked($compressor, '-f', $path);
	}
	my $compressed = "$path.gz";
	die "Compression did not produce $compressed\n" unless -s $compressed;
	return $compressed;
}

sub sra_role {
	my ($run, $filename) = @_;
	return 'r1' if $filename eq $run.'_1.fastq';
	return 'r2' if $filename eq $run.'_2.fastq';
	return 'single' if $filename eq $run.'.fastq';
	die "Unexpected fasterq-dump output for $run: $filename\n";
}

sub download_sra {
	my $prefetch = configured_executable(
		$opt{prefetch}, 'prefetch', 'NCBI SRA prefetch',
	);
	my $fasterq = configured_executable(
		$opt{fasterq_dump}, 'fasterq-dump', 'NCBI fasterq-dump',
	);
	my $validator = configured_executable(
		$opt{vdb_validate}, 'vdb-validate', 'NCBI vdb-validate',
	);
	my $compressor = configured_executable(
		$opt{compressor}, [qw(pigz gzip)], 'FASTQ compressor',
	);
	my $runs = ncbi_run_info();
	my $cache = File::Spec->catdir($opt{output_dir}, '.ncbi-cache');
	my $conversion_root = File::Spec->catdir($opt{output_dir}, '.fasterq');
	make_path($cache, $conversion_root);
	my (@files, @platforms);
	for my $record (@{$runs}) {
		my $run = $record->{run};
		push @platforms, $record->{platform} if $record->{platform} ne '';
		my $run_cache = File::Spec->catdir($cache, $run);
		run_checked($prefetch, '--max-size', 'u', '--output-directory', $run_cache, $run);
		my $sra_path = find_sra_path($run_cache, $run);
		my $validated = eval { run_checked($validator, $sra_path); 1 };
		if (!$validated) {
			my $validation_error = $@ || 'unknown validation failure';
			remove_tree($run_cache);
			die "NCBI SRA validation failed for $run; removed the invalid cache: "
				.$validation_error;
		}
		my $conversion = File::Spec->catdir($conversion_root, $run);
		remove_tree($conversion) if -d $conversion;
		make_path($conversion);
		run_checked(
			$fasterq, '--split-3', '-e', $opt{threads},
			'--outdir', $conversion, '-t', $conversion, $sra_path,
		);
		opendir my $dh, $conversion
			or die "Cannot inspect fasterq-dump output $conversion: $!\n";
		my @fastqs = sort grep {
			/^\Q$run\E(?:_[12])?\.fastq$/
				&& -f File::Spec->catfile($conversion, $_)
		} readdir $dh;
		closedir $dh;
		die "fasterq-dump produced no FASTQ files for $run\n" unless @fastqs;
		my %role_counts;
		for my $filename (@fastqs) {
			my $role = sra_role($run, $filename);
			my $plain = File::Spec->catfile($conversion, $filename);
			my $records = validate_fastq($plain);
			my $compressed = compress_fastq($plain, $compressor);
			my $destination_name = "$filename.gz";
			my $destination = File::Spec->catfile($opt{output_dir}, $destination_name);
			unlink $destination if -e $destination;
			move($compressed, $destination)
				or die "Cannot publish $destination: $!\n";
			my $file_record = {
				run => $run, filename => $destination_name, role => $role,
				layout => uc($record->{layout} || 'UNKNOWN'),
				bytes => 0 + (-s $destination),
				md5 => md5_file($destination), records => $records,
			};
			die "Published SRA FASTQ failed validation: $destination\n"
				unless validate_record($file_record, $opt{output_dir});
			push @files, $file_record;
			$role_counts{$role} = $records;
		}
		die "SRA paired FASTQs contain different record counts for $run\n"
			if exists($role_counts{r1}) != exists($role_counts{r2})
				|| (exists($role_counts{r1}) && $role_counts{r1} != $role_counts{r2});
	}
	remove_tree($cache, $conversion_root);
	return (\@files, \@platforms);
}

my ($files, $platforms) = $opt{provider} eq 'ena'
	? download_ena() : download_sra();
die "Archive download produced no validated FASTQ files\n"
	unless ref($files) eq 'ARRAY' && @{$files};

my %platform_seen;
my @unique_platforms = grep { !$platform_seen{$_}++ } @{$platforms};
my $total_bytes = 0;
$total_bytes += $_->{bytes} for @{$files};
my $metadata = {
	schema => 1, provider => $opt{provider}, accession => $accession_key,
	requested_accessions => $accessions,
	downloaded_at => strftime('%Y-%m-%dT%H:%M:%SZ', gmtime()),
	seqtech => infer_seqtech(@unique_platforms),
	platforms => \@unique_platforms, total_bytes => 0 + $total_bytes,
	files => $files,
};
my $json = JSON::PP->new->canonical->pretty->encode($metadata);
atomic_write($internal_manifest, $json);
atomic_write($opt{metadata_out}, $json);
printf "Downloaded and validated %d FASTQ file(s) for %s:%s (%.2f GiB, SeqTech=%s)\n",
	scalar(@{$files}), $opt{provider}, $accession_key,
	$total_bytes / (1024 ** 3), $metadata->{seqtech};
