package Mods::Checkpoint;

use strict;
use warnings;
use Exporter qw(import);
use File::Basename qw(dirname);
use File::Path qw(make_path);
use IO::Handle;
use JSON::PP qw(decode_json);

our @EXPORT_OK = qw(write_checkpoint checkpoint_valid read_checkpoint);

my $FORMAT = 'matafiler-checkpoint-v1';

sub _stringify_parameters {
	my ($parameters) = @_;
	$parameters ||= {};
	die "checkpoint parameters must be a hash reference\n"
		unless ref($parameters) eq 'HASH';
	return { map { $_ => defined($parameters->{$_}) ? "$parameters->{$_}" : '' } keys %$parameters };
}

sub read_checkpoint {
	my ($file, $knownStat) = @_;
	return unless defined $file;
	my @fileStat = ref($knownStat) eq 'ARRAY' ? @{$knownStat} : stat($file);
	return unless @fileStat && $fileStat[7] > 0;
	open my $fh, '<', $file or return;
	local $/;
	my $json = <$fh>;
	close $fh or return;
	my $manifest = eval { decode_json($json) };
	return if $@ || ref($manifest) ne 'HASH' || ($manifest->{format} // '') ne $FORMAT;
	return $manifest;
}

sub checkpoint_valid {
	my ($file, %options) = @_;
	return 0 unless defined $file;
	my @fileStat = stat($file);
	return 0 unless @fileStat;
	return 1 unless $fileStat[7] > 0; # legacy empty stones remain resumable

	my $manifest = read_checkpoint($file, \@fileStat) or return 0;
	my $expected = _stringify_parameters($options{parameters});
	my $actual = $manifest->{parameters};
	return 0 unless ref($actual) eq 'HASH';
	for my $key (keys %$expected) {
		return 0 unless exists $actual->{$key} && "$actual->{$key}" eq $expected->{$key};
	}

	return 0 unless ref($manifest->{outputs}) eq 'ARRAY';
	for my $record (@{$manifest->{outputs}}) {
		return 0 unless ref($record) eq 'HASH' && defined $record->{path};
		my $path = $record->{path};
		my @stat = stat($path);
		return 0 unless @stat && $stat[7] == ($record->{size} // -1);
	}
	return 1;
}

sub write_checkpoint {
	my ($file, %options) = @_;
	die "checkpoint path is required\n" unless defined $file && length $file;
	my $parameters = _stringify_parameters($options{parameters});
	my $outputs = $options{outputs} || [];
	die "checkpoint outputs must be an array reference\n" unless ref($outputs) eq 'ARRAY';

	my @output_records;
	for my $path (@$outputs) {
		die "checkpoint output path is undefined\n" unless defined $path && length $path;
		my @stat = stat($path);
		die "Cannot checkpoint missing output $path\n" unless @stat;
		push @output_records, { path => "$path", size => 0 + $stat[7], mtime => 0 + $stat[9] };
	}

	my $parent = dirname($file);
	make_path($parent) unless -d $parent;
	my $manifest = {
		format => $FORMAT,
		created_epoch => 0 + time,
		pid => 0 + $$,
		parameters => $parameters,
		outputs => \@output_records,
	};
	my $json = JSON::PP->new->ascii->canonical->pretty->encode($manifest);
	my $partial = "$file.part.$$";
	open my $fh, '>', $partial or die "Cannot create checkpoint partial $partial: $!\n";
	print {$fh} $json or die "Cannot write checkpoint partial $partial: $!\n";
	$fh->flush() or die "Cannot flush checkpoint partial $partial: $!\n";
	$fh->sync() or die "Cannot synchronize checkpoint partial $partial: $!\n";
	close $fh or die "Cannot close checkpoint partial $partial: $!\n";
	unless (rename $partial, $file) {
		my $rename_error = $!;
		# POSIX replaces atomically; retain a Windows compatibility fallback.
		if (-e $file) {
			unlink $file or die "Cannot replace checkpoint $file: $!\n";
			rename $partial, $file or die "Cannot publish checkpoint $file: $!\n";
		} else {
			die "Cannot publish checkpoint $file: $rename_error\n";
		}
	}
	return $file;
}

1;
