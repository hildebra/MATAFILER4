#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use File::Basename qw(dirname);
use File::Path qw(make_path);
use File::Temp qw(tempdir);

use Mods::FuncTools qw(readGene2Func);
use Mods::GenoMetaAss qw(
	readClstrRev readFasta resolveExistingFile
);
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::MosaicLoci qw(
	select_interesting_records discover_mosaic_candidates
	read_paf_stream read_paf_hits confirm_mosaic_candidates
	select_outgroup_panel
);

my $VERSION = '0.14';
my %DEFAULT = (
	cluster_id => 95,
	threads => 20,
	min_identity => 0.90,
	min_query_coverage => 0.80,
	min_target_coverage => 0.80,
	min_score_margin => 0.02,
	min_length_ratio => 0.80,
	max_sample_overlap_fraction => 0.15,
	max_overlap_samples => 1,
	outgroup_min_identity => 0.80,
	outgroup_max_identity => 0.95,
	outgroup_target_identity => 0.88,
	outgroup_min_coverage => 0.75,
	outgroup_min_loci => 10,
	minimap_preset => 'asm20',
	max_secondary_hits => 50,
);

my ($GCd, $mgs_file, $output, $candidate_output, $paf, $tmp_base, $help);
my $cluster_id = $DEFAULT{cluster_id};
my $threads = $DEFAULT{threads};
Getopt::Long::Configure(qw(no_auto_abbrev no_ignore_case));
GetOptions(
	'GCd=s' => \$GCd,
	'MGS=s' => \$mgs_file,
	'output=s' => \$output,
	'candidates=s' => \$candidate_output,
	'paf=s' => \$paf,
	'tmpD=s' => \$tmp_base,
	'clusterID=i' => \$cluster_id,
	'threads=i' => \$threads,
	'minIdentity=f' => \$DEFAULT{min_identity},
	'minQueryCoverage=f' => \$DEFAULT{min_query_coverage},
	'minTargetCoverage=f' => \$DEFAULT{min_target_coverage},
	'minScoreMargin=f' => \$DEFAULT{min_score_margin},
	'minLengthRatio=f' => \$DEFAULT{min_length_ratio},
	'maxSampleOverlapFraction=f' => \$DEFAULT{max_sample_overlap_fraction},
	'maxOverlapSamples=i' => \$DEFAULT{max_overlap_samples},
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
die "-clusterID must be between 1 and 100\n"
	unless $cluster_id >= 1 && $cluster_id <= 100;
die "-threads must be positive\n" unless $threads > 0;
die "Mosaic identity/coverage/margin options must be between zero and one\n"
	if grep { $_ < 0 || $_ > 1 } @DEFAULT{qw(
		min_identity min_query_coverage min_target_coverage min_score_margin
		min_length_ratio max_sample_overlap_fraction
		outgroup_min_identity outgroup_max_identity
		outgroup_target_identity outgroup_min_coverage
	)};
die "-maxOverlapSamples must be non-negative\n"
	unless $DEFAULT{max_overlap_samples} >= 0;
die "Outgroup identity interval is reversed\n"
	if $DEFAULT{outgroup_min_identity} >= $DEFAULT{outgroup_max_identity};
die "-outgroupTargetIdentity must lie within the permitted outgroup interval\n"
	if $DEFAULT{outgroup_target_identity} < $DEFAULT{outgroup_min_identity}
		|| $DEFAULT{outgroup_target_identity} > $DEFAULT{outgroup_max_identity};
die "-outgroupMinLoci must be positive\n" unless $DEFAULT{outgroup_min_loci} > 0;
die "-minimapPreset must be asm5, asm10, or asm20\n"
	unless $DEFAULT{minimap_preset} =~ /\Aasm(?:5|10|20)\z/;
die "-maxTargets must be positive\n" unless $DEFAULT{max_secondary_hits} > 0;

$candidate_output ||= "$output.candidates.tsv";
make_path(dirname($output), dirname($candidate_output));

my ($all_records) = read_mgs_records($mgs_file, $GCd);
die "No annotated MGS genes were found in $mgs_file\n" unless @{$all_records};
my ($records, $interest_statistics) = select_interesting_records($all_records);
print "Selected ".scalar(@{$records})."/".scalar(@{$all_records})
	." annotated genes with mosaic or outgroup comparison potential\n";
unless (@{$records}) {
	write_candidate_table($candidate_output, []);
	write_confirmed_catalogue($output, [], {}, {}, \%DEFAULT);
	write_rejections("$output.rejected.tsv", []);
	write_summary("$output.summary.tsv", {
		version => $VERSION,
		%{$interest_statistics},
		annotated_records => scalar(@{$all_records}),
		minimap_preset => $DEFAULT{minimap_preset},
		threads => $threads,
		max_secondary_hits => $DEFAULT{max_secondary_hits},
	});
	write_outgroup_table("$output.outgroups.tsv", {}, {});
	print "\nMosaic preprocessing summary (v$VERSION)\n";
	print "  Annotated input genes:       ".scalar(@{$all_records})."\n";
	print "  Comparison-interest genes:  0\n";
	print "    Mosaic-interest genes:     ".($interest_statistics->{mosaic_interest_records} || 0)."\n";
	print "    Outgroup-interest genes:   ".($interest_statistics->{outgroup_interest_records} || 0)."\n";
	print "    Unshared genes excluded:   ".($interest_statistics->{unshared_records} || 0)."\n";
	print "  Minimap queries/alignments:  0/0 (not run)\n";
	print "  Confirmed/rejected mosaics: 0/0\n";
	print "  Unique MGS-outgroup links:  0\n";
	print "  Proposed outgroup gene links: 0\n";
	print "  Candidate table:            $candidate_output\n";
	print "  Confirmed catalogue:        $output\n";
	print "  Outgroup proposals:         $output.outgroups.tsv\n";
	print "  Diagnostic metrics:         $output.summary.tsv\n";
	warn "No MGS genes shared an annotation within or between MGS; no minimap2 comparison was required\n";
	exit 0;
}
$all_records = [];
my %selected = map { $_->{gene} => 1 } @{$records};
my $cluster_index = resolveExistingFile("$GCd/compl.incompl.$cluster_id.fna.clstr.idx")
	or die "Catalogue cluster index is missing for identity $cluster_id\n";
my (undef, $cluster_members) = readClstrRev($cluster_index, 0, \%selected, undef);
my $catalogue_fasta = resolveExistingFile("$GCd/compl.incompl.$cluster_id.fna")
	or die "Catalogue nucleotide FASTA is missing for identity $cluster_id\n";
my $sequences = readFasta($catalogue_fasta, 1, "\\s", \%selected);

my %candidate_statistics;
my $candidates = discover_mosaic_candidates(
	$records, $cluster_members, $sequences,
	{
		minimum_length_ratio => $DEFAULT{min_length_ratio},
		maximum_sample_overlap_fraction => $DEFAULT{max_sample_overlap_fraction},
		maximum_overlap_samples => $DEFAULT{max_overlap_samples},
		statistics => \%candidate_statistics,
	},
);
write_candidate_table($candidate_output, $candidates);

my $work_dir = tempdir('mosaic-loci-XXXXXX',
	defined($tmp_base) && length($tmp_base) ? (DIR => $tmp_base) : (),
	CLEANUP => 1);
my $query_fasta = "$work_dir/selected_genes.fna";
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
my $hits;
if (defined($paf) && length($paf)) {
	die "Whole-catalogue alignment is missing or empty: $paf\n" unless -s $paf;
	$hits = read_paf_hits($paf, \%paf_options);
} else {
	my $minimap2 = getProgPaths('minimap2');
	print "Aligning ".scalar(@{$records})." interesting genes together in one minimap2 run "
		."using $threads threads and at most $DEFAULT{max_secondary_hits} secondary targets\n";
	my @command = (
		$minimap2, '-x', $DEFAULT{minimap_preset}, '-c', '-D',
		'--secondary=yes', '-N', $DEFAULT{max_secondary_hits},
		'-t', $threads, $catalogue_fasta, $query_fasta,
	);
	open my $minimap_fh, '-|', @command
		or die "Cannot start minimap2 whole-catalogue alignment: $!\n";
	$hits = read_paf_stream($minimap_fh, 'minimap2 stdout', \%paf_options);
	close $minimap_fh
		or die "minimap2 whole-catalogue alignment failed with status $?\n";
}
warn "minimap2 found no non-diagonal catalogue alignments; continuing with empty biological results\n"
	unless ($alignment_statistics{raw_alignments} || 0) > 0;
my ($confirmed, $rejected) = confirm_mosaic_candidates(
	$candidates, $hits,
	{
		minimum_identity => $DEFAULT{min_identity},
		minimum_query_coverage => $DEFAULT{min_query_coverage},
		minimum_target_coverage => $DEFAULT{min_target_coverage},
		minimum_score_margin => $DEFAULT{min_score_margin},
	},
);
my ($outgroups, $outgroup_genes) = select_outgroup_panel(
	$records, $hits,
	{
		minimum_identity => $DEFAULT{outgroup_min_identity},
		maximum_identity => $DEFAULT{outgroup_max_identity},
		target_identity => $DEFAULT{outgroup_target_identity},
		minimum_coverage => $DEFAULT{outgroup_min_coverage},
		minimum_loci => $DEFAULT{outgroup_min_loci},
	},
);
write_confirmed_catalogue($output, $confirmed, $outgroups, $outgroup_genes, \%DEFAULT);
write_rejections("$output.rejected.tsv", $rejected);
write_outgroup_table("$output.outgroups.tsv", $outgroups, $outgroup_genes);
my $alignment_count = 0;
$alignment_count += scalar(@{$_}) for values %{$hits};
my $outgroup_gene_count = 0;
$outgroup_gene_count += scalar(keys %{$_}) for values %{$outgroup_genes};
write_summary("$output.summary.tsv", {
	version => $VERSION,
	%{$interest_statistics},
	%candidate_statistics,
	%alignment_statistics,
	annotated_records => $interest_statistics->{input_records},
	selected_sequences => scalar(keys %{$sequences}),
	queries_with_alignments => scalar(keys %{$hits}),
	alignments => $alignment_count,
	confirmed_mosaics => scalar(@{$confirmed}),
	rejected_mosaics => scalar(@{$rejected}),
	outgroup_mgs => scalar(keys %{$outgroups}),
	outgroup_genes => $outgroup_gene_count,
	minimap_preset => $DEFAULT{minimap_preset},
	threads => $threads,
	max_secondary_hits => $DEFAULT{max_secondary_hits},
});

print "\nMosaic preprocessing summary (v$VERSION)\n";
print "  Annotated input genes:       $interest_statistics->{input_records}\n";
print "  Comparison-interest genes:  ".scalar(@{$records})."\n";
print "    Mosaic-interest genes:     ".($interest_statistics->{mosaic_interest_records} || 0)."\n";
print "    Outgroup-interest genes:   ".($interest_statistics->{outgroup_interest_records} || 0)."\n";
print "    Unshared genes excluded:   ".($interest_statistics->{unshared_records} || 0)."\n";
print "  Same-NOG pairs considered:  ".($candidate_statistics{same_cog_pairs} || 0)."\n";
print "    Length filtered:           ".($candidate_statistics{length_filtered_pairs} || 0)."\n";
print "    Membership missing:        ".($candidate_statistics{missing_membership_pairs} || 0)."\n";
print "    Sample-overlap filtered:   ".($candidate_statistics{overlap_filtered_pairs} || 0)."\n";
print "    Mosaic candidates:         ".scalar(@{$candidates})."\n";
print "  Minimap queries/alignments:  ".scalar(keys %{$hits})."/$alignment_count "
	."($threads threads, $DEFAULT{minimap_preset}, max targets $DEFAULT{max_secondary_hits})\n";
print "    Raw PAF records:           ".($alignment_statistics{raw_alignments} || 0)."\n";
print "    Self/threshold filtered:   ".($alignment_statistics{self_alignments_filtered} || 0)
	."/".($alignment_statistics{threshold_alignments_filtered} || 0)."\n";
print "    Duplicate target records:  "
	.(($alignment_statistics{duplicate_alignments_replaced} || 0)
		+ ($alignment_statistics{duplicate_alignments_filtered} || 0))."\n";
print "  Confirmed/rejected mosaics: ".scalar(@{$confirmed})."/".scalar(@{$rejected})."\n";
print "  Unique MGS-outgroup links:  ".scalar(keys %{$outgroups})."\n";
print "  Proposed outgroup gene links: $outgroup_gene_count\n";
print "  Candidate table:            $candidate_output\n";
print "  Confirmed catalogue:        $output\n";
print "  Outgroup proposals:         $output.outgroups.tsv\n";
print "  Diagnostic metrics:         $output.summary.tsv\n";
warn "No mosaic pairs or outgroups passed selection; inspect $output.summary.tsv for the filtering stage responsible\n"
	unless @{$confirmed} || keys %{$outgroups};
exit 0;

sub read_mgs_records {
	my ($path, $catalogue_dir) = @_;
	my (@assignments, %seen, %rank, %wanted);
	open my $fh, '<', $path or die "Cannot open MGS file $path: $!\n";
	while (my $line = <$fh>) {
		$line =~ s/[\r\n]+$//;
		next if $line eq '' || $line =~ /^#/;
		my @fields = split /\t/, $line, -1;
		die "Malformed MGS row in $path: $line\n" unless @fields >= 2 && length($fields[0]);
		my ($mgs, @genes);
		if ($fields[1] =~ /,/) {
			$mgs = $fields[0];
			@genes = grep { /^\d+$/ } split /,/, $fields[1];
		} elsif ($fields[1] =~ /^\d+$/) {
			($mgs, @genes) = ($fields[0], $fields[1]);
		} else {
			next;
		}
		for my $gene (@genes) {
			next if $seen{$gene}++;
			push @assignments, {
				mgs => $mgs, gene => $gene, rank => $rank{$mgs}++,
			};
			$wanted{$gene} = 1;
		}
	}
	close $fh or die "Cannot close MGS file $path: $!\n";
	my $gene_to_cog = readGene2Func($catalogue_dir, 'NOG', \%wanted);
	my @records;
	for my $assignment (@assignments) {
		my $cog = $gene_to_cog->{$assignment->{gene}};
		next unless defined($cog) && length($cog);
		push @records, {%{$assignment}, cog => $cog};
	}
	return \@records;
}

sub write_candidate_table {
	my ($path, $rows) = @_;
	my $temporary = "$path.tmp.$$";
	open my $fh, '>', $temporary or die "Cannot create $temporary: $!\n";
	print {$fh} join("\t", qw(
		MGS COG gene1 gene2 length_ratio gene1_samples gene2_samples
		overlap_samples overlap_fraction status
	)), "\n";
	for my $row (@{$rows}) {
		print {$fh} join("\t",
			$row->{mgs}, $row->{cog}, $row->{left}, $row->{right},
			sprintf('%.5f', $row->{length_ratio}), $row->{left_samples},
			$row->{right_samples}, $row->{overlap_samples},
			sprintf('%.5f', $row->{overlap_fraction}), 'candidate'), "\n";
	}
	close $fh or die "Cannot close $temporary: $!\n";
	rename $temporary, $path or die "Cannot install $path: $!\n";
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
		version annotated_records interesting_records mosaic_interest_records
		outgroup_interest_records unshared_records selected_sequences input_records same_cog_pairs
		missing_sequence_pairs length_filtered_pairs missing_membership_pairs
		overlap_filtered_pairs candidates queries_with_alignments alignments
		raw_alignments self_alignments_filtered threshold_alignments_filtered
		duplicate_alignments_replaced duplicate_alignments_filtered
		queries_with_retained_alignments retained_alignments
		confirmed_mosaics rejected_mosaics outgroup_mgs outgroup_genes
		minimap_preset threads max_secondary_hits
	)) {
		my $value = defined($statistics->{$metric}) ? $statistics->{$metric} : 0;
		print {$fh} "$metric\t$value\n";
	}
	close $fh or die "Cannot close $temporary: $!\n";
	rename $temporary, $path or die "Cannot install $path: $!\n";
}

sub write_outgroup_table {
	my ($path, $outgroups, $gene_map) = @_;
	my $temporary = "$path.tmp.$$";
	open my $fh, '>', $temporary or die "Cannot create $temporary: $!\n";
	print {$fh} join("\t", qw(
		source_MGS target_MGS homologous_loci median_identity proposed_gene_pairs
	)), "\n";
	for my $source (sort keys %{$outgroups}) {
		my $entry = $outgroups->{$source};
		my @pairs = map {
			$_.'->'.$gene_map->{$source}{$_}{target}
		} sort keys %{$gene_map->{$source} || {}};
		print {$fh} join("\t",
			$source, $entry->{target_mgs}, $entry->{loci},
			sprintf('%.5f', $entry->{median_identity}), join(',', @pairs),
		), "\n";
	}
	close $fh or die "Cannot close $temporary: $!\n";
	rename $temporary, $path or die "Cannot install $path: $!\n";
}

sub write_rejections {
	my ($path, $rows) = @_;
	my $temporary = "$path.tmp.$$";
	open my $fh, '>', $temporary or die "Cannot create $temporary: $!\n";
	print {$fh} join("\t", qw(MGS COG gene1 gene2 reason)), "\n";
	for my $row (@{$rows}) {
		print {$fh} join("\t", @{$row}{qw(mgs cog left right reason)}), "\n";
	}
	close $fh or die "Cannot close $temporary: $!\n";
	rename $temporary, $path or die "Cannot install $path: $!\n";
}

sub usage {
	my ($error) = @_;
	my $prefix = defined($error) ? "Error: $error\n\n" : '';
	return $prefix.<<"USAGE";
Usage: prepare_mosaic_loci.pl -GCd DIR -MGS FILE -output FILE [options]

Creates FILE.candidates.tsv first, aligns every selected gene against the
complete nucleotide gene catalogue, and writes only reciprocal, unique,
well-covered >=$DEFAULT{min_identity}-identity mosaic pairs to FILE. The same alignments are used
to choose an outgroup MGS represented across at least $DEFAULT{outgroup_min_loci} loci and with median
identity between $DEFAULT{outgroup_min_identity} and $DEFAULT{outgroup_max_identity}.

  -clusterID INT                 Gene-catalogue clustering identity [$DEFAULT{cluster_id}]
  -threads INT                   minimap2 threads [$DEFAULT{threads}]
  -paf FILE                      Reuse a PAF alignment instead of running minimap2
  -minIdentity FLOAT             Partner nucleotide identity [$DEFAULT{min_identity}]
  -minQueryCoverage FLOAT        Query alignment coverage [$DEFAULT{min_query_coverage}]
  -minTargetCoverage FLOAT       Target alignment coverage [$DEFAULT{min_target_coverage}]
  -minScoreMargin FLOAT          Required lead over another catalogue hit [$DEFAULT{min_score_margin}]
  -minLengthRatio FLOAT          Candidate prefilter length ratio [$DEFAULT{min_length_ratio}]
  -maxSampleOverlapFraction FLOAT  Maximum co-occurrence fraction [$DEFAULT{max_sample_overlap_fraction}]
  -maxOverlapSamples INT         Maximum rare co-occurring samples [$DEFAULT{max_overlap_samples}]
  -outgroupMinIdentity FLOAT     Closest permitted outgroup identity [$DEFAULT{outgroup_min_identity}]
  -outgroupMaxIdentity FLOAT     Most similar permitted outgroup identity [$DEFAULT{outgroup_max_identity}]
  -outgroupTargetIdentity FLOAT  Preferred median outgroup identity [$DEFAULT{outgroup_target_identity}]
  -outgroupMinCoverage FLOAT     Per-locus outgroup alignment coverage [$DEFAULT{outgroup_min_coverage}]
  -outgroupMinLoci INT           Required homologous outgroup loci [$DEFAULT{outgroup_min_loci}]
  -minimapPreset asm5|asm10|asm20  Whole-catalogue alignment preset [$DEFAULT{minimap_preset}]
  -maxTargets INT                 Maximum secondary targets per query [$DEFAULT{max_secondary_hits}]
USAGE
}
