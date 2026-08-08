package Mods::WorkflowResilience;

use strict;
use warnings;

use Exporter qw(import);
use File::Basename qw(dirname);
use File::Path qw(make_path);

our @EXPORT_OK = qw(
	retry_operation retry_unlink retry_rename retry_open retry_close
	preflight_executable preflight_directory filesystem_capacity
);

sub retry_operation {
	my (%options) = @_;
	my $code = $options{code};
	die "retry_operation requires a code reference\n" unless ref($code) eq 'CODE';
	my $label = $options{label} || 'operation';
	my $attempts = defined($options{attempts}) ? 0 + $options{attempts} : 5;
	$attempts = 1 if $attempts < 1;
	my $delays = $options{delays} || [1, 2, 5, 10];
	my $sleeper = $options{sleeper} || sub { sleep($_[0]); };
	my ($last_error, $result);
	for my $attempt (1 .. $attempts) {
		my $ok = eval {
			$result = $code->($attempt);
			die "$label returned an unsuccessful result\n" unless $result;
			1;
		};
		return $result if $ok;
		$last_error = $@ || "$label failed: $!";
		$last_error =~ s/\s+\z//;
		last if $attempt == $attempts;
		my $delay = $delays->[$attempt - 1];
		$delay = $delays->[-1] unless defined $delay;
		$delay = 1 unless defined($delay) && $delay >= 0;
		warn "$label failed on attempt $attempt/$attempts; retrying in ${delay}s: $last_error\n";
		$sleeper->($delay);
	}
	my $message = "$label failed after $attempts attempt(s): $last_error\n";
	if ($options{fatal} // 1) {
		die $message;
	}
	warn $message;
	return;
}

sub retry_unlink {
	my ($path, %options) = @_;
	return 1 unless defined($path) && (-e $path || -l $path);
	return retry_operation(
		%options,
		label => $options{label} || "remove $path",
		code => sub { return 1 unless -e $path || -l $path; return unlink($path) },
	);
}

sub retry_rename {
	my ($source, $destination, %options) = @_;
	return retry_operation(
		%options,
		label => $options{label} || "publish $destination",
		code => sub {
			return 1 if -e $destination && !-e $source;
			return 0 unless rename($source, $destination);
			return -e $destination ? 1 : 0;
		},
	);
}

sub retry_open {
	my ($mode, $path, %options) = @_;
	my $handle;
	retry_operation(
		%options,
		label => $options{label} || "open $path",
		code => sub { return open($handle, $mode, $path) },
	);
	return $handle;
}

sub retry_close {
	my ($handle, $label, %options) = @_;
	return retry_operation(
		%options,
		label => $label || 'close file',
		code => sub { return close($handle) },
	);
}

sub preflight_executable {
	my ($program, $label) = @_;
	die "Missing configured program for $label\n" unless defined($program) && length($program);
	if ($program =~ m{[\\/]}) {
		die "Configured $label executable does not exist: $program\n" unless -e $program;
		die "Configured $label executable is not executable: $program\n" unless -x $program;
		return $program;
	}
	for my $directory (split /:/, $ENV{PATH} // '') {
		next unless length $directory;
		my $candidate = "$directory/$program";
		return $candidate if -f $candidate && -x $candidate;
	}
	die "Configured $label executable is not available on PATH: $program\n";
}

sub preflight_directory {
	my ($directory, $label) = @_;
	$label ||= 'workflow directory';
	make_path($directory) unless -d $directory;
	die "$label is not a directory: $directory\n" unless -d $directory;
	my $probe = "$directory/.matafiler-preflight-$$";
	my $published = "$probe.published";
	my $handle = retry_open('>', $probe, label => "$label write probe");
	print {$handle} "preflight\n" or die "Cannot write $label probe $probe: $!\n";
	retry_close($handle, "$label probe close");
	retry_rename($probe, $published, label => "$label rename probe");
	retry_unlink($published, label => "$label probe cleanup", fatal => 0);
	return 1;
}

sub _df_value {
	my ($flag, $directory) = @_;
	open my $df, '-|', 'df', $flag, $directory or return;
	my @lines = <$df>;
	close $df or return;
	return unless @lines >= 2;
	my @field = split /\s+/, $lines[-1];
	@field = grep { length } @field;
	return unless @field >= 4 && $field[-3] =~ /^\d+$/;
	return 0 + $field[-3];
}

sub filesystem_capacity {
	my ($directory) = @_;
	my $available_kb = _df_value('-Pk', $directory);
	my $available_inodes = _df_value('-Pi', $directory);
	return {
		available_kb => $available_kb,
		available_inodes => $available_inodes,
	};
}

1;
