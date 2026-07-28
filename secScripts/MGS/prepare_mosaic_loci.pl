#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use File::Basename qw(dirname);
use File::Path qw(make_path);
use File::Temp qw(tempdir);

use Mods::FuncTools qw(readGene2Func);
use Mods::GenoMetaAss qw(
	readClstrRev readFasta resolveExistingFile systemW
);
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::MosaicLoci qw(
	discover_mosaic_candidates read_paf_hits confirm_mosaic_candidates
	select_outgroup_panel
);

my $VERSION = '0.10';
my %DEFAULT = (
	cluster_id => 95,
	threads => 4,
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

$candidate_output ||= "$output.candidates.tsv";
make_path(dirname($output), dirname($candidate_output));

my ($records, $gene_to_record) = read_mgs_records($mgs_file, $GCd);
die "No annotated MGS genes were found in $mgs_file\n" unless @{$records};
my %selected = map { $_->{gene} => 1 } @{$records};
my $cluster_index = resolveExistingFile("$GCd/compl.incompl.$cluster_id.fna.clstr.idx")
	or die "Catalogue cluster index is missing for identity $cluster_id\n";
my (undef, $cluster_members) = readClstrRev($cluster_index, 0, \%selected, undef);
my $catalogue_fasta = resolveExistingFile("$GCd/compl.incompl.$cluster_id.fna")
	or die "Catalogue nucleotide FASTA is missing for identity $cluster_id\n";
my $sequences = readFasta($catalogue_fasta, 1, "\\s", \%selected);

my $candidates = discover_mosaic_candidates(
	$records, $cluster_members, $sequences,
	{
		minimum_length_ratio => $DEFAULT{min_length_ratio},
		maximum_sample_overlap_fraction => $DEFAULT{max_sample_overlap_fraction},
		maximum_overlap_samples => $DEFAULT{max_overlap_samples},
	},
);
write_candidate_table($candidate_output, $candidates);

my $work_dir = tempdir('mosaic-loci-XXXXXX',
	defined($tmp_base) && length($tmp_base) ? (DIR => $tmp_base) : (),
	CLEANUP => 1);
my $query_fasta = "$work_dir/selected_genes.fna";
write_query_fasta($query_fasta, $records, $sequences);
my $paf_file = $paf;
unless (defined($paf_file) && length($paf_file)) {
	$paf_file = "$work_dir/selected_vs_catalogue.paf";
	my $minimap2 = getProgPaths('minimap2');
	my $command = join(' ',
		shell_quote($minimap2), '-x', 'asm10', '--secondary=yes', '-N', '100',
		'-t', $threads, shell_quote($catalogue_fasta), shell_quote($query_fasta),
		'>', shell_quote($paf_file));
	systemW("$command\n");
}
die "Whole-catalogue alignment is missing or empty: $paf_file\n" unless -s $paf_file;
my $hits = read_paf_hits($paf_file);
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

print "Mosaic preprocessing v$VERSION complete: ".scalar(@{$candidates})." candidate pair(s), "
	.scalar(@{$confirmed})." confirmed, ".scalar(keys %{$outgroups})
	." MGS with catalogue-derived outgroups\n";
print "Candidate table: $candidate_output\nConfirmed catalogue: $output\n";
exit 0;

sub read_mgs_records {
	my ($path, $catalogue_dir) = @_;
	my $gene_to_cog = readGene2Func($catalogue_dir, 'NOG');
	my (@records, %seen, %gene_record, %rank);
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
			my $cog = $gene_to_cog->{$gene};
			next unless defined($cog) && length($cog);
			my $record = {
				mgs => $mgs, gene => $gene, cog => $cog, rank => $rank{$mgs}++,
			};
			push @records, $record;
			$gene_record{$gene} = $record;
		}
	}
	close $fh or die "Cannot close MGS file $path: $!\n";
	return (\@records, \%gene_record);
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
	print {$fh} "# MATAFILER confirmed mosaic-locus and outgroup catalogue\n";
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

sub shell_quote {
	my ($value) = @_;
	$value = '' unless defined $value;
	$value =~ s/'/'"'"'/g;
	return "'$value'";
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
USAGE
}
