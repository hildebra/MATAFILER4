#!/usr/bin/env perl
use strict;
use warnings;

use Digest::SHA;
use Fcntl qw(:flock);
use File::Path qw(make_path remove_tree);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Getopt::Long qw(GetOptions);

use lib File::Spec->catdir($Bin, '..', '..');
use Mods::GenoMetaAss qw(ensureFastaIndex);
use Mods::IO_Tamoc_progs qw(getProgPaths);

my ($nt_source, $aa_source, $nt_ids, $aa_ids, $out_dir, $mosaic, $samtools, $help);
Getopt::Long::Configure(qw(no_auto_abbrev no_ignore_case));
GetOptions(
	'nt=s' => \$nt_source,
	'aa=s' => \$aa_source,
	'ntIDs=s' => \$nt_ids,
	'aaIDs=s' => \$aa_ids,
	'outD=s' => \$out_dir,
	'mosaic=s' => \$mosaic,
	'samtools=s' => \$samtools,
	'help|h' => \$help,
) or die usage();
if ($help) { print usage(); exit 0; }
die usage('unexpected positional arguments: '.join(' ', @ARGV)) if @ARGV;
for my $required (
	['-nt', $nt_source], ['-aa', $aa_source], ['-ntIDs', $nt_ids],
	['-aaIDs', $aa_ids], ['-outD', $out_dir],
) {
	die usage("$required->[0] is required")
		unless defined($required->[1]) && length($required->[1]);
}
for my $source ([$nt_source, 'nucleotide FASTA'], [$aa_source, 'protein FASTA'],
		[$nt_ids, 'nucleotide ID list'], [$aa_ids, 'protein ID list']) {
	die "$source->[1] is missing or empty: $source->[0]\n" unless -s $source->[0];
}
$samtools ||= getProgPaths('samtools');
die "Configured samtools command must be one executable path, not shell syntax: $samtools\n"
	if $samtools =~ /[\r\n;]/;

$out_dir = File::Spec->rel2abs($out_dir);
make_path($out_dir) unless -d $out_dir;
open my $cache_lock, '>>', File::Spec->catfile($out_dir, '.cache.lock')
	or die "Cannot open outgroup-reference cache lock in $out_dir: $!\n";
flock($cache_lock, LOCK_EX)
	or die "Cannot lock outgroup-reference cache $out_dir: $!\n";

my ($nt_count, $nt_digest) = validate_id_file($nt_ids, 'nucleotide');
my ($aa_count, $aa_digest) = validate_id_file($aa_ids, 'protein');
my $manifest_text = manifest_text(
	nt_source => $nt_source, aa_source => $aa_source,
	nt_count => $nt_count, aa_count => $aa_count,
	nt_digest => $nt_digest, aa_digest => $aa_digest,
	mosaic => $mosaic,
);
my $manifest = File::Spec->catfile($out_dir, 'manifest.tsv');
my $complete = File::Spec->catfile($out_dir, 'complete.sto');
my $nt_output = File::Spec->catfile($out_dir, 'references.fna');
my $aa_output = File::Spec->catfile($out_dir, 'references.faa');
if (-s $complete && -s $manifest && -s $nt_output && -s "$nt_output.fai"
		&& -s $aa_output && -s "$aa_output.fai"
		&& slurp($manifest) eq $manifest_text) {
	my $retainedNT = count_index_records("$nt_output.fai");
	my $retainedAA = count_index_records("$aa_output.fai");
	print "OUTGROUP_REFERENCE_CACHE status=reused requested_nt=$nt_count retained_nt=$retainedNT "
		."requested_aa=$aa_count retained_aa=$retainedAA dir=$out_dir\n";
	exit 0;
}

ensureFastaIndex($nt_source, { samtools => $samtools });
ensureFastaIndex($aa_source, { samtools => $samtools });
my $build_dir = tempdir('.build-XXXXXX', DIR => $out_dir, CLEANUP => 0);
my $ok = eval {
	my $new_nt = File::Spec->catfile($build_dir, 'references.fna');
	my $new_aa = File::Spec->catfile($build_dir, 'references.faa');
	my $availableNT = File::Spec->catfile($build_dir, 'available.nt.ids');
	my $availableAA = File::Spec->catfile($build_dir, 'available.aa.ids');
	my ($retainedNT, $missingNT) = filter_available_ids(
		$nt_ids, "$nt_source.fai", $availableNT, 'nucleotide');
	my ($retainedAA, $missingAA) = filter_available_ids(
		$aa_ids, "$aa_source.fai", $availableAA, 'protein');
	die "None of the requested nucleotide outgroup references exist in $nt_source\n"
		unless $retainedNT;
	die "None of the requested protein outgroup references exist in $aa_source\n"
		unless $retainedAA;
	extract_indexed_fasta($samtools, $nt_source, $availableNT, $new_nt, 'nucleotide');
	extract_indexed_fasta($samtools, $aa_source, $availableAA, $new_aa, 'protein');
	run_command($samtools, 'faidx', $new_nt);
	run_command($samtools, 'faidx', $new_aa);
	validate_extracted_index($availableNT, "$new_nt.fai", $retainedNT, 'nucleotide');
	validate_extracted_index($availableAA, "$new_aa.fai", $retainedAA, 'protein');
	warn "Selective reference cache omitted unavailable IDs: nucleotide=$missingNT, protein=$missingAA\n"
		if $missingNT || $missingAA;
	write_file(File::Spec->catfile($build_dir, 'manifest.tsv'), $manifest_text);
	write_file(File::Spec->catfile($build_dir, 'complete.sto'),
		join("\t", 'outgroup-reference-cache-v1', time, $$)."\n");
	unlink $complete if -e $complete;
	for my $name (qw(references.fna references.fna.fai references.faa references.faa.fai manifest.tsv)) {
		my $source = File::Spec->catfile($build_dir, $name);
		my $destination = File::Spec->catfile($out_dir, $name);
		rename $source, $destination
			or die "Cannot publish $source as $destination: $!\n";
	}
	my $new_complete = File::Spec->catfile($build_dir, 'complete.sto');
	rename $new_complete, $complete
		or die "Cannot publish outgroup-reference completion marker $complete: $!\n";
	1;
};
my $error = $@;
remove_tree($build_dir) if -d $build_dir;
die $error unless $ok;
my $retainedNT = count_index_records("$nt_output.fai");
my $retainedAA = count_index_records("$aa_output.fai");
print "OUTGROUP_REFERENCE_CACHE status=created requested_nt=$nt_count retained_nt=$retainedNT "
	."requested_aa=$aa_count retained_aa=$retainedAA dir=$out_dir\n";

exit 0;


sub extract_indexed_fasta {
	my ($binary, $source, $ids, $output, $label) = @_;
	open my $result, '>', $output or die "Cannot create $label cache $output: $!\n";
	my $pid = open my $pipe, '-|', $binary, 'faidx', '-r', $ids, $source;
	die "Cannot start samtools faidx for $label references: $!\n" unless defined $pid;
	while (my $chunk = <$pipe>) {
		print {$result} $chunk or die "Cannot write $label cache $output: $!\n";
	}
	my $closed = close $pipe;
	my $status = $?;
	close $result or die "Cannot close $label cache $output: $!\n";
	die command_error("samtools faidx extraction for $label", $status)
		unless $closed && $status == 0;
	die "samtools created an empty $label reference cache: $output\n" unless -s $output;
}

sub validate_id_file {
	my ($path, $label) = @_;
	open my $input, '<', $path or die "Cannot read $label ID list $path: $!\n";
	my $digest = Digest::SHA->new(256);
	my ($count, $previous) = (0, undef);
	while (my $line = <$input>) {
		$digest->add($line);
		$line =~ s/[\r\n]+\z//;
		die "Blank identifier in $label ID list $path\n" unless length $line;
		die "Unsafe identifier '$line' in $label ID list $path\n"
			if $line =~ /[\s\x00-\x1f\x7f]/;
		die "Duplicate adjacent identifier '$line' in $label ID list $path\n"
			if defined($previous) && $line eq $previous;
		$previous = $line;
		$count++;
	}
	close $input or die "Cannot close $label ID list $path: $!\n";
	die "The $label ID list is empty: $path\n" unless $count;
	return ($count, $digest->hexdigest);
}

sub filter_available_ids {
	my ($requested, $source_index, $available, $label) = @_;
	open my $request_fh, '<', $requested
		or die "Cannot read requested $label IDs $requested: $!\n";
	my %wanted;
	while (my $line = <$request_fh>) {
		$line =~ s/[\r\n]+\z//;
		$wanted{$line} = 1 if length $line;
	}
	close $request_fh or die "Cannot close requested $label IDs $requested: $!\n";
	open my $index_fh, '<', $source_index
		or die "Cannot read source $label FASTA index $source_index: $!\n";
	open my $available_fh, '>', $available
		or die "Cannot create available $label ID list $available: $!\n";
	my $retained = 0;
	while (my $line = <$index_fh>) {
		my ($identifier) = split /\t/, $line, 2;
		next unless exists $wanted{$identifier};
		print {$available_fh} "$identifier\n"
			or die "Cannot write available $label ID list $available: $!\n";
		delete $wanted{$identifier};
		$retained++;
	}
	close $index_fh or die "Cannot close source $label FASTA index $source_index: $!\n";
	close $available_fh or die "Cannot close available $label ID list $available: $!\n";
	return ($retained, scalar(keys %wanted));
}

sub count_index_records {
	my ($index) = @_;
	open my $input, '<', $index or die "Cannot read FASTA index $index: $!\n";
	my $count = 0;
	$count++ while <$input>;
	close $input or die "Cannot close FASTA index $index: $!\n";
	return $count;
}

sub validate_extracted_index {
	my ($ids, $index, $expected_count, $label) = @_;
	open my $wanted, '<', $ids or die "Cannot reopen $label ID list $ids: $!\n";
	open my $observed, '<', $index or die "Cannot read extracted $label index $index: $!\n";
	my $count = 0;
	while (my $wanted_line = <$wanted>) {
		$wanted_line =~ s/[\r\n]+\z//;
		my $index_line = <$observed>;
		die "Extracted $label FASTA is missing requested identifier '$wanted_line'\n"
			unless defined $index_line;
		my ($observed_id) = split /\t/, $index_line, 2;
		die "Extracted $label FASTA order/content mismatch: wanted '$wanted_line', observed '$observed_id'\n"
			unless defined($observed_id) && $observed_id eq $wanted_line;
		$count++;
	}
	die "Extracted $label FASTA contains unexpected extra records\n" if defined(<$observed>);
	close $wanted or die "Cannot close $label ID list $ids: $!\n";
	close $observed or die "Cannot close extracted $label index $index: $!\n";
	die "Extracted $label record count changed: expected=$expected_count observed=$count\n"
		unless $count == $expected_count;
}

sub manifest_text {
	my %arg = @_;
	my @lines = (join("\t", qw(key value)));
	push @lines, join("\t", 'format', 'outgroup-reference-cache-v1');
	for my $kind (qw(nt aa)) {
		my $source = File::Spec->rel2abs($arg{"${kind}_source"});
		my @stat = stat($source);
		push @lines,
			join("\t", "${kind}_source", $source),
			join("\t", "${kind}_source_size", $stat[7]),
			join("\t", "${kind}_source_mtime", $stat[9]),
			join("\t", "${kind}_ids_sha256", $arg{"${kind}_digest"}),
			join("\t", "${kind}_count", $arg{"${kind}_count"});
	}
	if (defined($arg{mosaic}) && length($arg{mosaic}) && -s $arg{mosaic}) {
		my $path = File::Spec->rel2abs($arg{mosaic});
		my @stat = stat($path);
		push @lines, join("\t", 'mosaic', $path),
			join("\t", 'mosaic_size', $stat[7]),
			join("\t", 'mosaic_mtime', $stat[9]);
	} else {
		push @lines, join("\t", 'mosaic', 'none');
	}
	return join("\n", @lines)."\n";
}

sub run_command {
	my @command = @_;
	my $status = system { $command[0] } @command;
	die "Cannot execute $command[0]: $!\n" if $status == -1;
	die command_error(join(' ', @command), $status) if $status != 0;
}

sub command_error {
	my ($label, $status) = @_;
	return "$label failed from signal ".($status & 127)."\n" if $status & 127;
	return "$label failed with exit code ".($status >> 8)."\n";
}


sub write_file {
	my ($path, $contents) = @_;
	open my $output, '>', $path or die "Cannot create $path: $!\n";
	print {$output} $contents or die "Cannot write $path: $!\n";
	close $output or die "Cannot close $path: $!\n";
}

sub slurp {
	my ($path) = @_;
	open my $input, '<', $path or die "Cannot read $path: $!\n";
	local $/;
	my $contents = <$input> // '';
	close $input or die "Cannot close $path: $!\n";
	return $contents;
}

sub usage {
	my ($error) = @_;
	my $message = <<'USAGE';
Usage: prepare_strain_outgroup_refs.pl -nt FILE -aa FILE
       -ntIDs FILE -aaIDs FILE -outD DIR [-mosaic FILE] [-samtools PATH]

Builds an atomic, restartable NT/AA outgroup-reference cache by indexed FASTA
lookup. Source and compact cache FASTAs receive samtools .fai indexes.
USAGE
	return defined($error) ? "$error\n$message" : $message;
}
