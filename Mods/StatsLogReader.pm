package Mods::StatsLogReader;
use strict;
use warnings;

use Exporter qw(import);
use IO::Uncompress::Gunzip qw($GunzipError);
use Mods::GenoMetaAss qw(resolveExistingFile);

our @EXPORT_OK = qw(
	read_stats_log_excerpt
	reset_stats_log_sampling
	stats_log_sampling_summary
);

my %sampling_summary;

sub reset_stats_log_sampling {
	%sampling_summary = (
		large_files => 0,
		compressed_files => 0,
		source_bytes => 0,
		retained_bytes => 0,
	);
}

sub stats_log_sampling_summary {
	return {%sampling_summary};
}

sub _record_sampling {
	my ($source_bytes, $retained_bytes, $compressed) = @_;
	$sampling_summary{large_files}++;
	$sampling_summary{compressed_files}++ if $compressed;
	$sampling_summary{source_bytes} += $source_bytes;
	$sampling_summary{retained_bytes} += $retained_bytes;
}

sub _read_bytes {
	my ($fh, $maximum) = @_;
	my $text = '';
	while (length($text) < $maximum) {
		my $wanted = $maximum - length($text);
		my $read = read($fh, my $chunk, $wanted > 64 * 1024 ? 64 * 1024 : $wanted);
		die "Cannot read statistics log: $!\n" unless defined $read;
		last if $read == 0;
		$text .= $chunk;
	}
	return $text;
}

sub _append_bounded_tail {
	my ($tail_ref, $chunk, $maximum) = @_;
	${$tail_ref} .= $chunk;
	my $excess = length(${$tail_ref}) - $maximum;
	substr(${$tail_ref}, 0, $excess, '') if $excess > 0;
}

sub _complete_head {
	my ($text, $truncated) = @_;
	$text =~ s/[^\n]*\z// if $truncated && $text =~ /\n/ && $text !~ /\n\z/;
	return $text;
}

sub _complete_tail {
	my ($text, $truncated) = @_;
	$text =~ s/\A[^\n]*\n// if $truncated && $text =~ /\n/;
	return $text;
}

sub _last_lines {
	my ($text, $count) = @_;
	return $text unless defined($count) && $count > 0;
	my @lines = split /\n/, $text, -1;
	pop @lines if @lines && $lines[-1] eq '';
	@lines = @lines[-$count .. -1] if @lines > $count;
	return join("\n", @lines);
}

sub _read_compressed_excerpt {
	my ($resolved, $mode, $maximum_bytes, $edge_bytes) = @_;
	my $fh = IO::Uncompress::Gunzip->new($resolved, MultiStream => 1);
	unless ($fh) {
		warn "Cannot open compressed statistics log '$resolved': $GunzipError\n";
		return ('', 0, 0);
	}

	my ($all, $head, $tail) = ('', '', '');
	my ($total_bytes, $sampled) = (0, 0);
	while (1) {
		my $read = read($fh, my $chunk, 64 * 1024);
		die "Cannot read compressed statistics log '$resolved': $!\n"
			unless defined $read;
		last if $read == 0;
		$total_bytes += $read;
		if ($mode eq 'tail') {
			_append_bounded_tail(\$tail, $chunk, $edge_bytes);
			next;
		}
		if (!$sampled) {
			$all .= $chunk;
			if (length($all) > $maximum_bytes) {
				$sampled = 1;
				$head = substr($all, 0, $edge_bytes);
				$tail = substr($all, -$edge_bytes);
				$all = '';
			}
		} else {
			_append_bounded_tail(\$tail, $chunk, $edge_bytes);
		}
	}
	close $fh or warn "Cannot close compressed statistics log '$resolved': $!\n";

	$sampled = 1 if $mode eq 'tail' && $total_bytes > $maximum_bytes;
	return ($all, $total_bytes, 0) if $mode ne 'tail' && !$sampled;

	my $text;
	if ($mode eq 'tail') {
		$text = _complete_tail($tail, $total_bytes > length($tail));
	} else {
		$head = _complete_head($head, $total_bytes > length($head));
		$tail = _complete_tail($tail, $total_bytes > length($tail));
		$text = $head . ($head ne '' && $tail ne '' ? "\n" : '') . $tail;
	}
	return ($text, $total_bytes, $sampled);
}

sub read_stats_log_excerpt {
	my ($path, %options) = @_;
	my $mode = $options{mode} || 'head_tail';
	die "Unknown statistics-log read mode '$mode'\n"
		unless $mode eq 'head_tail' || $mode eq 'tail';
	my $maximum_bytes = $options{max_file_bytes} // 5 * 1024 * 1024;
	my $edge_bytes = $options{edge_bytes} // 512 * 1024;
	die "Statistics-log size limits must be positive integers\n"
		unless $maximum_bytes =~ /^\d+$/ && $maximum_bytes > 0
			&& $edge_bytes =~ /^\d+$/ && $edge_bytes > 0;
	$edge_bytes = int($maximum_bytes / 2) if $edge_bytes * 2 > $maximum_bytes;
	$edge_bytes = 1 if $edge_bytes < 1;

	my ($resolved, $file_stat) = resolveExistingFile($path);
	return '' unless defined $resolved;
	my $compressed = $resolved =~ /\.gz\z/ ? 1 : 0;
	my ($text, $source_bytes, $sampled);

	if ($compressed) {
		($text, $source_bytes, $sampled) = _read_compressed_excerpt(
			$resolved, $mode, $maximum_bytes, $edge_bytes,
		);
	} else {
		open my $fh, '<:raw', $resolved or do {
			warn "Cannot read statistics log '$resolved': $!\n";
			return '';
		};
		$source_bytes = $file_stat->[7];
		if ($mode eq 'tail' && $source_bytes > $edge_bytes) {
			my $start = $source_bytes - $edge_bytes;
			seek($fh, $start, 0) or die "Cannot seek in statistics log '$resolved': $!\n";
			$text = _complete_tail(_read_bytes($fh, $edge_bytes), 1);
			$sampled = $source_bytes > $maximum_bytes ? 1 : 0;
		} elsif ($mode eq 'head_tail' && $source_bytes > $maximum_bytes) {
			my $head = _complete_head(_read_bytes($fh, $edge_bytes), 1);
			my $start = $source_bytes - $edge_bytes;
			seek($fh, $start, 0) or die "Cannot seek in statistics log '$resolved': $!\n";
			my $tail = _complete_tail(_read_bytes($fh, $edge_bytes), 1);
			$text = $head . ($head ne '' && $tail ne '' ? "\n" : '') . $tail;
			$sampled = 1;
		} else {
			$text = _read_bytes($fh, $source_bytes);
			$sampled = 0;
		}
		close $fh or warn "Cannot close statistics log '$resolved': $!\n";
	}

	$text = _last_lines($text, $options{tail_lines}) if $mode eq 'tail';
	_record_sampling($source_bytes, length($text), $compressed) if $sampled;
	return $text;
}

reset_stats_log_sampling();

1;
