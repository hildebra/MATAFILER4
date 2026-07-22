package Mods::GTDBTaxonomy;

use strict;
use warnings;
use Exporter qw(import);

our @EXPORT_OK = qw(merge_gtdb_summaries read_gtdb_taxonomy);

sub _temporary_path {
	my ($path) = @_;
	return "$path.tmp.$$";
}

sub _remove_if_present {
	my ($path) = @_;
	return unless -e $path;
	unlink $path or die "Cannot remove temporary taxonomy file $path: $!\n";
}

sub merge_gtdb_summaries {
	my ($summaries, $summary_out, $taxonomy_out) = @_;
	die "GTDB summary paths must be an array reference\n"
		unless ref($summaries) eq 'ARRAY' && @{$summaries};
	die "GTDB summary and taxonomy output paths are required\n"
		unless defined($summary_out) && length($summary_out)
			&& defined($taxonomy_out) && length($taxonomy_out);

	my $summary_tmp = _temporary_path($summary_out);
	my $taxonomy_tmp = _temporary_path($taxonomy_out);
	_remove_if_present($summary_tmp);
	_remove_if_present($taxonomy_tmp);

	my $header;
	my $row_count = 0;
	my %seen_genome;
	my $ok = eval {
		open my $summary_fh, '>', $summary_tmp
			or die "Cannot write temporary GTDB summary $summary_tmp: $!\n";
		open my $taxonomy_fh, '>', $taxonomy_tmp
			or die "Cannot write temporary GTDB taxonomy $taxonomy_tmp: $!\n";

		for my $summary (@{$summaries}) {
			die "GTDB-Tk summary is missing or empty: $summary\n" unless -s $summary;
			open my $input_fh, '<', $summary
				or die "Cannot read GTDB-Tk summary $summary: $!\n";
			my $file_header = <$input_fh>;
			die "GTDB-Tk summary has no header: $summary\n" unless defined $file_header;
			$file_header =~ s/\r?\n\z//;
			my @header_fields = split /\t/, $file_header, -1;
			die "Invalid GTDB-Tk summary header in $summary\n"
				unless @header_fields >= 2
					&& $header_fields[0] eq 'user_genome'
					&& $header_fields[1] eq 'classification';

			if (!defined $header) {
				$header = $file_header;
				print {$summary_fh} "$header\n";
				print {$taxonomy_fh} "user_genome\tclassification\n";
			} elsif ($file_header ne $header) {
				die "GTDB-Tk summary headers differ between domain outputs ($summary)\n";
			}

			while (my $line = <$input_fh>) {
				$line =~ s/\r?\n\z//;
				next if $line eq '';
				my @fields = split /\t/, $line, -1;
				die "Malformed GTDB-Tk summary row in $summary: $line\n"
					unless @fields >= 2 && length($fields[0]) && length($fields[1]);
				my ($genome, $classification) = @fields[0, 1];
				die "Duplicate genome '$genome' across GTDB-Tk summaries\n"
					if $seen_genome{$genome}++;
				$classification =~ s/(^|;)[a-z]__/$1/gi;
				print {$summary_fh} "$line\n";
				print {$taxonomy_fh} "$genome\t$classification\n";
				$row_count++;
			}
			close $input_fh or die "Cannot close GTDB-Tk summary $summary: $!\n";
		}

		die "GTDB-Tk summaries contain no taxonomy data rows\n" unless $row_count;
		close $summary_fh or die "Cannot close temporary GTDB summary $summary_tmp: $!\n";
		close $taxonomy_fh or die "Cannot close temporary GTDB taxonomy $taxonomy_tmp: $!\n";
		1;
	};

	if (!$ok) {
		my $error = $@ || "Unknown GTDB taxonomy merge failure\n";
		_remove_if_present($summary_tmp);
		_remove_if_present($taxonomy_tmp);
		die $error;
	}

	rename $summary_tmp, $summary_out
		or die "Cannot publish GTDB summary $summary_out: $!\n";
	rename $taxonomy_tmp, $taxonomy_out
		or die "Cannot publish GTDB taxonomy $taxonomy_out: $!\n";
	return $row_count;
}

sub read_gtdb_taxonomy {
	my ($path) = @_;
	die "GTDB taxonomy path is required\n" unless defined($path) && length($path);
	open my $fh, '<', $path or die "Cannot read GTDB taxonomy $path: $!\n";

	my %taxonomy;
	my $row_count = 0;
	while (my $line = <$fh>) {
		$line =~ s/\r?\n\z//;
		next if $line eq '';
		my @fields = split /\t/, $line, -1;
		if ($fields[0] eq 'user_genome') {
			die "Invalid GTDB taxonomy header in $path\n"
				unless @fields >= 2 && $fields[1] eq 'classification';
			next; # tolerate duplicate headers from historical bac120/ar53 concatenation
		}
		die "Malformed GTDB taxonomy row in $path: $line\n"
			unless @fields >= 2 && length($fields[0]) && length($fields[1]);
		my ($id, $classification) = @fields[0, 1];
		die "Duplicate GTDB taxonomy assignment for '$id' in $path\n"
			if exists $taxonomy{$id};
		$classification = '?' if $classification =~ /Unclassified Bacteria/;
		my @ranks = split /;/, $classification, -1;
		push @ranks, '?' while @ranks < 7;
		$taxonomy{$id} = join(';', @ranks);
		$row_count++;
	}
	close $fh or die "Cannot close GTDB taxonomy $path: $!\n";
	die "GTDB taxonomy contains no data rows: $path\n" unless $row_count;
	return \%taxonomy;
}

1;
