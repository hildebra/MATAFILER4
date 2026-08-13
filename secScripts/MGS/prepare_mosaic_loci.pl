#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use File::Basename qw(dirname);
use File::Copy qw(copy);
use File::Path qw(make_path);
use File::Temp qw(tempdir);

use Mods::FuncTools qw(readGene2Func);
use Mods::GenoMetaAss qw(
	readFasta resolveExistingFile
);
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::MosaicLoci qw(
	read_paf_hits
	select_outgroup_panel
);

my $VERSION = '0.17';
my %DEFAULT = (
	cluster_id => 95,
	threads => 20,
	min_identity => 0.90,
	min_query_coverage => 0.80,
	min_target_coverage => 0.80,
	min_score_margin => 0.02,
	min_length_ratio => 0.80,
	max_sample_overlap_fraction => 0.15,
	max_abundance_overlap_fraction => 0.15,
	outgroup_min_identity => 0.80,
	outgroup_max_identity => 0.95,
	outgroup_target_identity => 0.88,
	outgroup_min_coverage => 0.75,
	outgroup_min_loci => 10,
	minimap_preset => 'asm20',
	max_secondary_hits => 0,
);

my ($GCd, $mgs_file, $core_mgs_file, $output, $candidate_output, $paf, $matrix,
	$rtk_override, $tmp_base, $help);
my $cluster_id = $DEFAULT{cluster_id};
my $threads = $DEFAULT{threads};
Getopt::Long::Configure(qw(no_auto_abbrev no_ignore_case));
GetOptions(
	'GCd=s' => \$GCd,
	'MGS=s' => \$mgs_file,
	'coreMGS=s' => \$core_mgs_file,
	'output=s' => \$output,
	'candidates=s' => \$candidate_output,
	'paf=s' => \$paf,
	'matrix=s' => \$matrix,
	'rtk=s' => \$rtk_override,
	'tmpD=s' => \$tmp_base,
	'clusterID=i' => \$cluster_id,
	'threads=i' => \$threads,
	'minIdentity=f' => \$DEFAULT{min_identity},
	'minQueryCoverage=f' => \$DEFAULT{min_query_coverage},
	'minTargetCoverage=f' => \$DEFAULT{min_target_coverage},
	'minScoreMargin=f' => \$DEFAULT{min_score_margin},
	'minLengthRatio=f' => \$DEFAULT{min_length_ratio},
	'maxSampleOverlapFraction=f' => \$DEFAULT{max_sample_overlap_fraction},
	'maxAbundanceOverlapFraction=f' => \$DEFAULT{max_abundance_overlap_fraction},
	'outgroupMinIdentity=f' => \$DEFAULT{outgroup_min_identity},
	'outgroupMaxIdentity=f' => \$DEFAULT{outgroup_max_identity},
	'outgroupTargetIdentity=f' => \$DEFAULT{outgroup_target_identity},
	'outgroupMinCoverage=f' => \$DEFAULT{outgroup_min_coverage},
	'outgroupMinLoci=i' => \$DEFAULT{outgroup_min_loci},
	'minimapPreset=s' => \$DEFAULT{minimap_preset},
	'maxTargets=i' => \$DEFAULT{max_secondary_hits},
	'help|h' => \$help,
) or die usage();
if ($help) { print usage(); exit 0; }
die usage('unexpected positional arguments: '.join(' ', @ARGV)) if @ARGV;
die usage('-GCd, -MGS, and -output are required')
	unless defined($GCd) && defined($mgs_file) && defined($output);
die "Gene catalogue directory is missing: $GCd\n" unless -d $GCd;
die "MGS file is missing or empty: $mgs_file\n" unless -s $mgs_file;
$core_mgs_file ||= $mgs_file;
die "Core MGS file is missing or empty: $core_mgs_file\n" unless -s $core_mgs_file;
die "-clusterID must be between 1 and 100\n"
	unless $cluster_id >= 1 && $cluster_id <= 100;
die "-threads must be positive\n" unless $threads > 0;
die "Mosaic identity/coverage/margin options must be between zero and one\n"
	if grep { $_ < 0 || $_ > 1 } @DEFAULT{qw(
		min_identity min_query_coverage min_target_coverage min_score_margin
		min_length_ratio max_sample_overlap_fraction max_abundance_overlap_fraction
		outgroup_min_identity outgroup_max_identity
		outgroup_target_identity outgroup_min_coverage
	)};
die "Outgroup identity interval is reversed\n"
	if $DEFAULT{outgroup_min_identity} >= $DEFAULT{outgroup_max_identity};
die "-outgroupTargetIdentity must lie within the permitted outgroup interval\n"
	if $DEFAULT{outgroup_target_identity} < $DEFAULT{outgroup_min_identity}
		|| $DEFAULT{outgroup_target_identity} > $DEFAULT{outgroup_max_identity};
die "-outgroupMinLoci must be positive\n" unless $DEFAULT{outgroup_min_loci} > 0;
die "-minimapPreset must be asm5, asm10, or asm20\n"
	unless $DEFAULT{minimap_preset} =~ /\Aasm(?:5|10|20)\z/;
die "-maxTargets must be non-negative\n" unless $DEFAULT{max_secondary_hits} >= 0;

$candidate_output ||= "$output.candidates.tsv";
$matrix ||= "$GCd/Matrix.mat.gz";
die "Gene abundance matrix is missing or empty: $matrix\n" unless -s $matrix;
my $rtk = defined($rtk_override) && length($rtk_override)
	? $rtk_override : getProgPaths('rare');
die "rtk executable was not resolved\n" unless defined($rtk) && length($rtk);
make_path(dirname($output), dirname($candidate_output));

my $records = read_mgs_records($mgs_file);
die "No MGS genes were found in raw assignment file $mgs_file\n" unless @{$records};
my $core_records = read_mgs_records($core_mgs_file, $GCd, 1);
my @outgroup_records = grep {
	!$_->{ambiguous} && defined($_->{cog}) && length($_->{cog})
} @{$core_records};
print "Loaded ".scalar(@{$records})." raw MGS genes for comprehensive mosaic "
	."self-alignment; ".scalar(@outgroup_records)
	." annotated core genes remain eligible for outgroup selection\n";
my %selected = map { $_->{gene} => 1 } @{$records};
my $catalogue_fasta = resolveExistingFile("$GCd/compl.incompl.$cluster_id.fna")
	or die "Catalogue nucleotide FASTA is missing for identity $cluster_id\n";
my $sequences = readFasta($catalogue_fasta, 1, "\\s", \%selected, { fai => 1 });
die "No raw MGS genes had catalogue sequences in $catalogue_fasta\n"
	unless keys %{$sequences};

my $work_dir = tempdir('mosaic-loci-XXXXXX',
	defined($tmp_base) && length($tmp_base) ? (DIR => $tmp_base) : (),
	CLEANUP => 1);
my $query_fasta = "$work_dir/raw_mgs_genes.fna";
write_query_fasta($query_fasta, $records, $sequences);
my %alignment_statistics;
my $minimum_stream_identity = $DEFAULT{min_identity} < $DEFAULT{outgroup_min_identity}
	? $DEFAULT{min_identity} : $DEFAULT{outgroup_min_identity};
my $minimum_stream_query_coverage =
	$DEFAULT{min_query_coverage} < $DEFAULT{outgroup_min_coverage}
	? $DEFAULT{min_query_coverage} : $DEFAULT{outgroup_min_coverage};
my $minimum_stream_target_coverage =
	$DEFAULT{min_target_coverage} < $DEFAULT{outgroup_min_coverage}
	? $DEFAULT{min_target_coverage} : $DEFAULT{outgroup_min_coverage};
my %paf_options = (
	minimum_identity => $minimum_stream_identity,
	minimum_query_coverage => $minimum_stream_query_coverage,
	minimum_target_coverage => $minimum_stream_target_coverage,
	exclude_self => 1,
	statistics => \%alignment_statistics,
);
my $paf_path;
if (defined($paf) && length($paf)) {
	die "Raw-MGS self-alignment is missing: $paf\n" unless -e $paf;
	$paf_path = $paf;
	print "Reusing raw-MGS self-alignment PAF $paf_path\n";
} else {
	$paf_path = "$work_dir/raw_mgs.minimap2.paf";
	unless (-e $paf_path) {
		my $minimap2 = getProgPaths('minimap2');
		my $maximum_targets = $DEFAULT{max_secondary_hits} > 0
			? $DEFAULT{max_secondary_hits}
			: scalar(keys %{$sequences}) - 1;
		$maximum_targets = 1 if $maximum_targets < 1;
		print "Self-aligning ".scalar(keys %{$sequences})." raw MGS genes before rtk2 "
			."using $threads threads and at most $maximum_targets secondary targets\n";
		my @command = (
			$minimap2, '-x', $DEFAULT{minimap_preset}, '-c', '-D',
			'--secondary=yes', '-p', '0', '-N', $maximum_targets,
			'-t', $threads, $query_fasta, $query_fasta,
		);
		my $temporary_paf = "$paf_path.tmp.$$";
		open my $paf_fh, '>', $temporary_paf
			or die "Cannot create temporary PAF $temporary_paf: $!\n";
		open my $minimap_fh, '-|', @command
			or die "Cannot start minimap2 raw-MGS self-alignment: $!\n";
		while (my $line = <$minimap_fh>) {
			print {$paf_fh} $line or die "Cannot write $temporary_paf: $!\n";
		}
		my $output_closed = close $paf_fh;
		my $minimap_closed = close $minimap_fh;
		my $minimap_status = $?;
		unless ($output_closed && $minimap_closed) {
			unlink $temporary_paf;
			die "minimap2 raw-MGS self-alignment failed with status $minimap_status\n";
		}
		rename $temporary_paf, $paf_path
			or die "Cannot install PAF $paf_path: $!\n";
	} else {
		print "Reusing existing raw-MGS minimap2 PAF $paf_path\n";
	}
}

my $rtk_prefix = "$work_dir/rtk";
my @rtk_command = (
	$rtk, 'mosaic', '-i', $matrix, '-reference', $mgs_file,
	'-paf', $paf_path, '-o', $rtk_prefix,
	'-t', $threads,
	'-mosaicMinIdentity', $DEFAULT{min_identity},
	'-mosaicMinQueryCoverage', $DEFAULT{min_query_coverage},
	'-mosaicMinTargetCoverage', $DEFAULT{min_target_coverage},
	'-mosaicMinScoreMargin', $DEFAULT{min_score_margin},
	'-mosaicMinLengthRatio', $DEFAULT{min_length_ratio},
	'-mosaicMaxSampleOverlapFraction', $DEFAULT{max_sample_overlap_fraction},
	'-mosaicMaxAbundanceOverlapFraction', $DEFAULT{max_abundance_overlap_fraction},
	'-mosaicRequireCigar',
);
print "Running rtk2 mosaic against $matrix after minimap2\n";
system(@rtk_command) == 0
	or die "rtk2 mosaic failed with status $?\n";
my $rtk_report = "$rtk_prefix.mosaic.tsv";
my $rtk_summary_path = "$rtk_prefix.mosaic.summary.tsv";
die "rtk2 mosaic did not create $rtk_report\n" unless -s $rtk_report;
copy_atomic($rtk_report, $candidate_output);
my ($confirmed, $rejected) = read_rtk_mosaic_results($rtk_report);
my $rtk_statistics = read_metric_table($rtk_summary_path);

my $hits = read_paf_hits($paf_path, \%paf_options);
warn "minimap2 found no non-diagonal raw-MGS alignments; continuing with empty biological results\n"
	unless ($alignment_statistics{raw_alignments} || 0) > 0;
my ($outgroups, $outgroup_genes) = select_outgroup_panel(
	\@outgroup_records, $hits,
	{
		minimum_identity => $DEFAULT{outgroup_min_identity},
		maximum_identity => $DEFAULT{outgroup_max_identity},
		target_identity => $DEFAULT{outgroup_target_identity},
		minimum_coverage => $DEFAULT{outgroup_min_coverage},
		minimum_loci => $DEFAULT{outgroup_min_loci},
	},
);
write_confirmed_catalogue($output, $confirmed, $outgroups, $outgroup_genes, \%DEFAULT);
my $alignment_count = 0;
$alignment_count += scalar(@{$_}) for values %{$hits};
my $outgroup_gene_count = 0;
$outgroup_gene_count += scalar(keys %{$_}) for values %{$outgroup_genes};
write_summary("$output.summary.tsv", {
	version => $VERSION,
	%alignment_statistics,
	raw_mgs_records => scalar(@{$records}),
	core_outgroup_records => scalar(@outgroup_records),
	selected_sequences => scalar(keys %{$sequences}),
	queries_with_alignments => scalar(keys %{$hits}),
	alignments => $alignment_count,
	confirmed_mosaics => scalar(@{$confirmed}),
	rejected_mosaics => scalar(@{$rejected}),
	rtk_candidate_pairs => $rtk_statistics->{candidate_pairs} || 0,
	rtk_confirmed_pairs => $rtk_statistics->{confirmed_pairs} || 0,
	outgroup_mgs => scalar(keys %{$outgroups}),
	outgroup_genes => $outgroup_gene_count,
	minimap_preset => $DEFAULT{minimap_preset},
	threads => $threads,
	max_secondary_hits => $DEFAULT{max_secondary_hits},
});
remove_legacy_intermediates($output,
	defined($paf) && length($paf) ? $paf_path : undef);

print "\nMosaic preprocessing summary (v$VERSION)\n";
print "  Raw MGS input genes:         ".scalar(@{$records})."\n";
print "  Core outgroup genes:         ".scalar(@outgroup_records)."\n";
print "  Minimap queries/alignments:  ".scalar(keys %{$hits})."/$alignment_count "
	."($threads threads, $DEFAULT{minimap_preset}, "
	.($DEFAULT{max_secondary_hits} ? "max targets $DEFAULT{max_secondary_hits}" : "all raw-MGS targets").")\n";
print "    Raw PAF records:           ".($alignment_statistics{raw_alignments} || 0)."\n";
print "    Self/threshold filtered:   ".($alignment_statistics{self_alignments_filtered} || 0)
	."/".($alignment_statistics{threshold_alignments_filtered} || 0)."\n";
print "    Duplicate target records:  "
	.(($alignment_statistics{duplicate_alignments_replaced} || 0)
		+ ($alignment_statistics{duplicate_alignments_filtered} || 0))."\n";
print "  Confirmed/rejected mosaics: ".scalar(@{$confirmed})."/".scalar(@{$rejected})."\n";
print "  rtk2 candidate/confirmed:   ".($rtk_statistics->{candidate_pairs} || 0)."/"
	.($rtk_statistics->{confirmed_pairs} || 0)."\n";
print "  Unique MGS-outgroup links:  ".scalar(keys %{$outgroups})."\n";
print "  Proposed outgroup gene links: $outgroup_gene_count\n";
print "  Candidate table:            $candidate_output\n";
print "  Confirmed catalogue:        $output\n";
print "  Diagnostic metrics:         $output.summary.tsv\n";
warn "No mosaic pairs or outgroups passed selection; inspect $output.summary.tsv for the filtering stage responsible\n"
	unless @{$confirmed} || keys %{$outgroups};
exit 0;

sub read_mgs_records {
	my ($path, $catalogue_dir, $with_annotations) = @_;
	my (@assignments, %record_for_gene, %rank, %wanted);
	open my $fh, '<', $path or die "Cannot open MGS file $path: $!\n";
	while (my $line = <$fh>) {
		$line =~ s/[\r\n]+$//;
		next if $line eq '' || $line =~ /^#/;
		my @fields = split /\t/, $line, -1;
		die "Malformed MGS row in $path: $line\n" unless @fields >= 2 && length($fields[0]);
		next if $fields[0] =~ /\A(?:bin|mgs)\z/i
			&& $fields[1] =~ /gene/i;
		my $mgs = $fields[0];
		my @genes = split /,/, $fields[1];
		for my $gene (@genes) {
			$gene =~ s/^\s+|\s+$//g;
			next unless length($gene);
			if (my $existing = $record_for_gene{$gene}) {
				$existing->{ambiguous} = 1 if $existing->{mgs} ne $mgs;
				next;
			}
			my $assignment = {
				mgs => $mgs, gene => $gene, rank => $rank{$mgs}++, ambiguous => 0,
			};
			push @assignments, $assignment;
			$record_for_gene{$gene} = $assignment;
			$wanted{$gene} = 1 if $with_annotations;
		}
	}
	close $fh or die "Cannot close MGS file $path: $!\n";
	my $gene_to_cog = $with_annotations
		? readGene2Func($catalogue_dir, 'NOG', \%wanted) : {};
	my @records;
	for my $assignment (@assignments) {
		my $cog = $gene_to_cog->{$assignment->{gene}};
		push @records, {%{$assignment}, cog => defined($cog) ? $cog : ''};
	}
	return \@records;
}

sub write_query_fasta {
	my ($path, $records, $sequences) = @_;
	open my $fh, '>', $path or die "Cannot create $path: $!\n";
	for my $record (@{$records}) {
		my $sequence = $sequences->{$record->{gene}};
		next unless defined($sequence) && length($sequence);
		print {$fh} ">$record->{gene}\n$sequence\n";
	}
	close $fh or die "Cannot close $path: $!\n";
	die "No selected catalogue sequences could be written to $path\n" unless -s $path;
}

sub copy_atomic {
	my ($source, $destination) = @_;
	my $temporary = "$destination.tmp.$$";
	copy($source, $temporary)
		or die "Cannot copy $source to $temporary: $!\n";
	rename $temporary, $destination
		or die "Cannot install $destination: $!\n";
}

sub remove_legacy_intermediates {
	my ($prefix, $protected_paf) = @_;
	my @legacy = map { $prefix.$_ } qw(
		.minimap2.paf
		.rtk.mosaic.tsv
		.rtk.mosaic.summary.tsv
		.rtk.concat.list
		.rejected.tsv
		.outgroups.tsv
	);
	my $removed = 0;
	for my $path (@legacy) {
		next if defined($protected_paf) && $path eq $protected_paf;
		next unless -e $path;
		if (unlink $path) {
			$removed++;
		} else {
			warn "Cannot remove obsolete Mosaic intermediate $path: $!\n";
		}
	}
	print "Removed $removed obsolete persistent Mosaic intermediate(s)\n" if $removed;
}

sub read_metric_table {
	my ($path) = @_;
	die "rtk2 summary is missing or empty: $path\n" unless -s $path;
	open my $fh, '<', $path or die "Cannot open rtk2 summary $path: $!\n";
	my %metrics;
	my $header = <$fh>;
	die "rtk2 summary has no header: $path\n"
		unless defined($header) && $header =~ /^metric\tvalue/;
	while (my $line = <$fh>) {
		$line =~ s/[\r\n]+$//;
		next unless length($line);
		my ($metric, $value) = split /\t/, $line, 2;
		$metrics{$metric} = $value if defined($metric) && defined($value);
	}
	close $fh or die "Cannot close rtk2 summary $path: $!\n";
	return \%metrics;
}

sub read_rtk_mosaic_results {
	my ($path) = @_;
	open my $fh, '<', $path or die "Cannot open rtk2 mosaic report $path: $!\n";
	my $header = <$fh>;
	die "rtk2 mosaic report has no header: $path\n" unless defined($header);
	$header =~ s/[\r\n]+$//;
	my @header = split /\t/, $header, -1;
	my %field_index;
	@field_index{@header} = (0 .. $#header);
	for my $required (qw(
		mgs group left right status reason identity query_coverage target_coverage
	)) {
		die "rtk2 mosaic report lacks column $required: $path\n"
			unless exists($field_index{$required});
	}
	my (@confirmed, @rejected);
	my $line_number = 1;
	while (my $line = <$fh>) {
		$line_number++;
		$line =~ s/[\r\n]+$//;
		next unless length($line);
		my @fields = split /\t/, $line, -1;
		die "Short rtk2 mosaic row $line_number in $path\n" if @fields < @header;
		my $entry = {
			mgs => $fields[$field_index{mgs}],
			cog => $fields[$field_index{group}],
			left => $fields[$field_index{left}],
			right => $fields[$field_index{right}],
			reason => $fields[$field_index{reason}],
			identity => 0 + $fields[$field_index{identity}],
			query_coverage => 0 + $fields[$field_index{query_coverage}],
			target_coverage => 0 + $fields[$field_index{target_coverage}],
		};
		if ($fields[$field_index{status}] eq 'confirmed') {
			push @confirmed, $entry;
		} else {
			push @rejected, $entry;
		}
	}
	close $fh or die "Cannot close rtk2 mosaic report $path: $!\n";
	return (\@confirmed, \@rejected);
}

sub write_confirmed_catalogue {
	my ($path, $confirmed, $outgroups, $gene_map, $defaults) = @_;
	my $temporary = "$path.tmp.$$";
	open my $fh, '>', $temporary or die "Cannot create $temporary: $!\n";
	print {$fh} "# MATAFILER confirmed mosaic-locus and outgroup catalogue v$VERSION\n";
	print {$fh} "# min_identity=$defaults->{min_identity}; min_query_coverage=$defaults->{min_query_coverage}; "
		."min_target_coverage=$defaults->{min_target_coverage}; min_score_margin=$defaults->{min_score_margin}\n";
	for my $row (@{$confirmed}) {
		print {$fh} join("\t", 'MOSAIC', $row->{mgs}, $row->{cog},
			$row->{left}, $row->{right}, map { sprintf('%.5f', $_) }
				@{$row}{qw(identity query_coverage target_coverage)}), "\n";
	}
	for my $source (sort keys %{$outgroups}) {
		my $entry = $outgroups->{$source};
		print {$fh} join("\t", 'OUTGROUP', $source, $entry->{target_mgs},
			$entry->{loci}, sprintf('%.5f', $entry->{median_identity})), "\n";
		for my $query (sort keys %{$gene_map->{$source} || {}}) {
			my $hit = $gene_map->{$source}{$query};
			print {$fh} join("\t", 'OUTGROUP_GENE', $source, $entry->{target_mgs},
				$query, $hit->{target}, sprintf('%.5f', $hit->{identity}),
				sprintf('%.5f', $hit->{query_coverage})), "\n";
		}
	}
	close $fh or die "Cannot close $temporary: $!\n";
	rename $temporary, $path or die "Cannot install $path: $!\n";
}

sub write_summary {
	my ($path, $statistics) = @_;
	my $temporary = "$path.tmp.$$";
	open my $fh, '>', $temporary or die "Cannot create $temporary: $!\n";
	print {$fh} "metric\tvalue\n";
	for my $metric (qw(
		version raw_mgs_records core_outgroup_records selected_sequences
		queries_with_alignments alignments raw_alignments self_alignments_filtered threshold_alignments_filtered
		duplicate_alignments_replaced duplicate_alignments_filtered
		queries_with_retained_alignments retained_alignments
		confirmed_mosaics rejected_mosaics rtk_candidate_pairs rtk_confirmed_pairs
		outgroup_mgs outgroup_genes
		minimap_preset threads max_secondary_hits
	)) {
		my $value = defined($statistics->{$metric}) ? $statistics->{$metric} : 0;
		print {$fh} "$metric\t$value\n";
	}
	close $fh or die "Cannot close $temporary: $!\n";
	rename $temporary, $path or die "Cannot install $path: $!\n";
}

sub usage {
	my ($error) = @_;
	my $prefix = defined($error) ? "Error: $error\n\n" : '';
	return $prefix.<<"USAGE";
Usage: prepare_mosaic_loci.pl -GCd DIR -MGS FILE -output FILE [options]

Extracts every gene assigned in raw SB.clusters, aligns that raw MGS gene set
against itself, and runs rtk2 mosaic without a NOG/group restriction. Thus every
same-MGS sequence relationship can be assessed using abundance evidence. Only
rtk-confirmed pairs are published to FILE. Perl reuses the same PAF after rtk2
to choose an outgroup MGS represented across at least
$DEFAULT{outgroup_min_loci} loci and with median identity between
$DEFAULT{outgroup_min_identity} and $DEFAULT{outgroup_max_identity}.

Persistent outputs are FILE (the strain rerun checkpoint), FILE.candidates.tsv
(the complete auditable rtk2 decision table), and FILE.summary.tsv (small run
metrics). Generated FASTA, minimap2 PAF, and native rtk2 report/summary/concat
files stay in the temporary workspace and are removed when the run finishes.
An explicitly supplied -paf file is read but never removed.

  -coreMGS FILE                  Optional core table used only for outgroup genes
  -clusterID INT                 Gene-catalogue clustering identity [$DEFAULT{cluster_id}]
  -threads INT                   minimap2 threads [$DEFAULT{threads}]
  -paf FILE                      Reuse a raw-MGS self-alignment PAF
  -matrix FILE                   Gene abundance matrix [GCd/Matrix.mat.gz]
  -rtk FILE                      Override the configured rtk2 executable
  -minIdentity FLOAT             Partner nucleotide identity [$DEFAULT{min_identity}]
  -minQueryCoverage FLOAT        Query alignment coverage [$DEFAULT{min_query_coverage}]
  -minTargetCoverage FLOAT       Target alignment coverage [$DEFAULT{min_target_coverage}]
  -minScoreMargin FLOAT          Required lead over another raw-MGS hit [$DEFAULT{min_score_margin}]
  -minLengthRatio FLOAT          rtk2 candidate length ratio [$DEFAULT{min_length_ratio}]
  -maxSampleOverlapFraction FLOAT  Maximum co-occurrence fraction [$DEFAULT{max_sample_overlap_fraction}]
  -maxAbundanceOverlapFraction FLOAT  Maximum abundance overlap [$DEFAULT{max_abundance_overlap_fraction}]
  -outgroupMinIdentity FLOAT     Closest permitted outgroup identity [$DEFAULT{outgroup_min_identity}]
  -outgroupMaxIdentity FLOAT     Most similar permitted outgroup identity [$DEFAULT{outgroup_max_identity}]
  -outgroupTargetIdentity FLOAT  Preferred median outgroup identity [$DEFAULT{outgroup_target_identity}]
  -outgroupMinCoverage FLOAT     Per-locus outgroup alignment coverage [$DEFAULT{outgroup_min_coverage}]
  -outgroupMinLoci INT           Required homologous outgroup loci [$DEFAULT{outgroup_min_loci}]
  -minimapPreset asm5|asm10|asm20  Raw-MGS self-alignment preset [$DEFAULT{minimap_preset}]
  -maxTargets INT                Secondary targets per query; 0 retains all
                                 raw-MGS targets [$DEFAULT{max_secondary_hits}]
USAGE
}
