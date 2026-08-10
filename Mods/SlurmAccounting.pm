package Mods::SlurmAccounting;

use strict;
use warnings;
use Exporter qw(import);

our @EXPORT_OK = qw(
	slurm_tree_memory_summary format_slurm_tree_memory_summary
	next_oom_retry_memory_mb
);

sub next_oom_retry_memory_mb {
	my ($current_mb, $maximum_mb) = @_;
	die "Current and maximum OOM-retry memory must be positive numbers\n"
		unless defined($current_mb) && defined($maximum_mb)
			&& $current_mb =~ /\A(?:\d+(?:\.\d*)?|\.\d+)\z/
			&& $maximum_mb =~ /\A(?:\d+(?:\.\d*)?|\.\d+)\z/
			&& $current_mb > 0 && $maximum_mb > 0;
	return undef if $current_mb >= $maximum_mb;
	my $next_mb = int($current_mb * 2 + 0.5);
	$next_mb = int($maximum_mb + 0.5) if $next_mb > $maximum_mb;
	return $next_mb;
}


sub _run_sacct {
	my ($command, $options) = @_;
	my $retrySeconds = defined($options->{retry_seconds}) ? $options->{retry_seconds} : 300;
	my $maximumSeconds = defined($options->{maximum_error_seconds})
		? $options->{maximum_error_seconds} : 1_200;
	my $sleeper = $options->{sleeper} || sub { sleep($_[0]); };
	my $clock = $options->{clock} || sub { time };
	my $started;
	while (1) {
		my ($output, $status);
		if (my $runner = $options->{runner}) {
			($output, $status) = $runner->($command);
			$status //= 0;
		} else {
			$output = `$command 2>/dev/null`;
			$status = $?;
		}
		return ($output, 0) if $status == 0;
		$started //= $clock->();
		my $elapsed = $clock->() - $started;
		return ($output, $status) if $elapsed >= $maximumSeconds;
		my $diagnostic = $output // ''; $diagnostic =~ s/\s+\z//;
		warn "Transient Slurm accounting failure after ${elapsed}s; retrying in ${retrySeconds}s"
			.($diagnostic ne '' ? ": $diagnostic" : '')."\n";
		$sleeper->($retrySeconds);
	}
}
sub _memory_mb {
	my ($value) = @_;
	return unless defined $value;
	$value =~ s/^\s+|\s+$//g;
	return unless $value =~ /^([\d.]+)\s*([KMGTPE]?)(?:i?B)?$/i;
	my ($number, $unit) = (0 + $1, uc($2));
	my %factor = ('' => 1 / (1024 * 1024), K => 1 / 1024, M => 1,
		G => 1024, T => 1024**2, P => 1024**3, E => 1024**4);
	return $number * $factor{$unit};
}

sub _quantile {
	my ($fraction, @values) = @_;
	return unless @values;
	@values = sort { $a <=> $b } @values;
	my $index = int(@values * $fraction);
	$index = $#values if $index > $#values;
	return $values[$index];
}

sub slurm_tree_memory_summary {
	my ($records, $options) = @_;
	$options ||= {};
	my %requested = map { $_->{job_id} => $_ } grep {
		defined($_->{job_id}) && $_->{job_id} =~ /^\d+$/
	} @{$records || []};
	my @ids = sort { $a <=> $b } keys %requested;
	my (%max_rss, %oom, %seen);

	while (@ids) {
		my @batch = splice(@ids, 0, 1000);
		my $command = "sacct -n -P -j ".join(',', @batch)
			." --format=JobIDRaw,State,ExitCode,MaxRSS";
		my ($output, $status) = _run_sacct($command, $options);
		return { available => 0, error => "sacct failed with status $status" }
			if $status != 0;
		for my $line (split /\n/, $output // '') {
			my ($raw_id, $state, $exit_code, $rss) = split /\|/, $line, 5;
			next unless defined $raw_id;
			my ($job_id) = $raw_id =~ /^(\d+)(?:[._].*)?$/;
			next unless defined($job_id) && exists $requested{$job_id};
			$seen{$job_id} = 1;
			$state = uc($state // '');
			$state =~ s/\+$//;
			$oom{$job_id} = 1
				if $state eq 'OUT_OF_MEMORY' || $state eq 'OOM'
					|| $state =~ /OUT.?OF.?MEMORY/;
			my $rss_mb = _memory_mb($rss);
			$max_rss{$job_id} = $rss_mb
				if defined($rss_mb)
					&& (!defined($max_rss{$job_id}) || $rss_mb > $max_rss{$job_id});
		}
	}

	my (@headroom_mb, @headroom_percent, @details);
	for my $job_id (sort { $a <=> $b } keys %requested) {
		my $record = $requested{$job_id};
		my $used_mb = $max_rss{$job_id};
		if (defined $used_mb && $record->{requested_mb} > 0) {
			my $headroom = $record->{requested_mb} - $used_mb;
			push @headroom_mb, $headroom;
			push @headroom_percent, 100 * $headroom / $record->{requested_mb};
		}
		push @details, {
			%{$record}, max_rss_mb => $used_mb, oom => $oom{$job_id} ? 1 : 0,
			seen => $seen{$job_id} ? 1 : 0,
		};
	}
	my $mean_mb;
	my $mean_percent;
	if (@headroom_mb) {
		$mean_mb = 0; $mean_mb += $_ for @headroom_mb; $mean_mb /= @headroom_mb;
		$mean_percent = 0; $mean_percent += $_ for @headroom_percent;
		$mean_percent /= @headroom_percent;
	}
	return {
		available => 1, jobs => \@details,
		accounted => scalar(@headroom_mb),
		missing => scalar(keys(%requested)) - scalar(@headroom_mb),
		oom_jobs => [grep { $_->{oom} } @details],
		mean_headroom_mb => $mean_mb,
		mean_headroom_percent => $mean_percent,
		q95_headroom_mb => _quantile(0.95, @headroom_mb),
		q95_headroom_percent => _quantile(0.95, @headroom_percent),
	};
}

sub format_slurm_tree_memory_summary {
	my ($summary) = @_;
	return "Tree memory accounting unavailable: ".($summary->{error} || 'unknown error')."\n"
		unless $summary->{available};
	my $text = "\nSLURM tree memory accounting:\n";
	my @oom = @{$summary->{oom_jobs} || []};
	if (@oom) {
		$text .= "  OOM events: ".scalar(@oom)." ("
			.join(', ', map { "$_->{mgs} [job $_->{job_id}]" } @oom).")\n";
	} else {
		$text .= "  OOM events: 0\n";
	}
	if ($summary->{accounted}) {
		$text .= sprintf(
			"  MaxRSS available: %d job(s); unavailable: %d\n"
			."  Estimated mean memory oversubscription (requested - MaxRSS): %.1f MiB (%.1f%% of request)\n"
			."  95th quantile memory under-/over-usage: %.1f MiB (%.1f%%; positive means under-used, negative means over-used)\n",
			$summary->{accounted}, $summary->{missing},
			$summary->{mean_headroom_mb}, $summary->{mean_headroom_percent},
			$summary->{q95_headroom_mb}, $summary->{q95_headroom_percent},
		);
	} else {
		$text .= "  MaxRSS is not yet available for any tracked tree job; utilization statistics cannot be estimated.\n";
	}
	return $text;
}

1;
