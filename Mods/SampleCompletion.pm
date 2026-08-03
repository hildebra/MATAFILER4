package Mods::SampleCompletion;
use strict;
use warnings;

use Digest::SHA qw(sha256_hex);
use Exporter qw(import);
use File::Basename qw(dirname);
use File::Path qw(make_path);
use File::Spec;
use JSON::PP ();

our @EXPORT_OK = qw(
	sample_completion_path completion_request_signature
	read_sample_completion write_sample_completion
	invalidate_sample_completion
);

my $SENTINEL_NAME = 'MATAFILER.sample.complete.json';
my $SCHEMA_NAME = 'MATAFILER.sample-completion';
my $SCHEMA_VERSION = 1;

sub _sample_root {
	my ($root) = @_;
	die "sample completion root is required\n"
		unless defined($root) && length($root);
	my $canonical = File::Spec->canonpath(File::Spec->rel2abs($root));
	my $volume_root = File::Spec->rootdir();
	die "refusing unsafe sample completion root '$root'\n"
		if $canonical eq $volume_root;
	return $canonical;
}

sub sample_completion_path {
	my ($root) = @_;
	return File::Spec->catfile(_sample_root($root), $SENTINEL_NAME);
}

sub completion_request_signature {
	my ($request) = @_;
	die "completion request must be a hash reference\n"
		unless ref($request) eq 'HASH';
	my $json = JSON::PP->new->canonical(1)->utf8(1)->encode($request);
	return sha256_hex($json);
}

sub _json_codec {
	return JSON::PP->new->canonical(1)->pretty(1)->utf8(1);
}

sub read_sample_completion {
	my (%args) = @_;
	my $path = defined($args{path}) && length($args{path})
		? $args{path} : sample_completion_path($args{root});
	return wantarray ? (undef, 'missing') : undef unless -e $path;
	return wantarray ? (undef, 'empty sentinel') : undef unless -s $path;

	open my $fh, '<:raw', $path
		or return wantarray ? (undef, "cannot read sentinel: $!") : undef;
	local $/;
	my $json = <$fh>;
	close $fh
		or return wantarray ? (undef, "cannot close sentinel: $!") : undef;
	my $record = eval { _json_codec()->decode($json) };
	return wantarray ? (undef, "invalid JSON: $@") : undef
		unless ref($record) eq 'HASH';

	my $error = '';
	if (($record->{schema} || '') ne $SCHEMA_NAME) {
		$error = 'unknown sentinel schema';
	} elsif (($record->{schema_version} || 0) != $SCHEMA_VERSION) {
		$error = 'unsupported sentinel schema version';
	} elsif (!defined($record->{sample}) || $record->{sample} eq '') {
		$error = 'sentinel sample is missing';
	} elsif (defined($args{sample}) && $record->{sample} ne $args{sample}) {
		$error = "sentinel belongs to sample '$record->{sample}'";
	} elsif (!defined($record->{request_signature})
			|| $record->{request_signature} !~ /^[0-9a-f]{64}$/) {
		$error = 'sentinel request signature is invalid';
	} elsif (defined($args{request_signature})
			&& $record->{request_signature} ne $args{request_signature}) {
		$error = 'requested workflow differs from the closed workflow';
	} elsif (ref($record->{metagstats}) ne 'HASH'
			|| !defined($record->{metagstats}{DIR})
			|| ref($record->{metagstats}{values}) ne 'HASH') {
		$error = 'sentinel metagStats record is incomplete';
	}
	return wantarray ? (undef, $error) : undef if $error ne '';
	return wantarray ? ($record, '') : $record;
}

sub write_sample_completion {
	my (%args) = @_;
	for my $required (qw(root sample request_signature metagstats)) {
		die "sample completion $required is required\n"
			unless exists($args{$required}) && defined($args{$required});
	}
	die "sample completion request signature is invalid\n"
		unless $args{request_signature} =~ /^[0-9a-f]{64}$/;
	die "sample completion metagStats record must contain DIR and values\n"
		unless ref($args{metagstats}) eq 'HASH'
			&& defined($args{metagstats}{DIR})
			&& ref($args{metagstats}{values}) eq 'HASH';

	my $path = sample_completion_path($args{root});
	my $directory = dirname($path);
	make_path($directory) unless -d $directory;
	my $record = {
		schema => $SCHEMA_NAME,
		schema_version => $SCHEMA_VERSION,
		sample => $args{sample},
		request_signature => $args{request_signature},
		created_epoch => time,
		present_assembly => $args{present_assembly} ? 1 : 0,
		empty_sample => $args{empty_sample} ? 1 : 0,
		empty_input_size_mb => 0 + ($args{empty_input_size_mb} || 0),
		metagstats => $args{metagstats},
	};
	my $temporary = "$path.tmp.$$";
	my $json = _json_codec()->encode($record);
	open my $fh, '>:raw', $temporary
		or die "cannot write sample completion sentinel $temporary: $!\n";
	print {$fh} $json
		or die "cannot populate sample completion sentinel $temporary: $!\n";
	close $fh
		or die "cannot close sample completion sentinel $temporary: $!\n";
	rename $temporary, $path
		or die "cannot publish sample completion sentinel $path: $!\n";
	return $path;
}

sub invalidate_sample_completion {
	my ($root) = @_;
	my $path = sample_completion_path($root);
	return 0 unless -e $path;
	unlink $path or die "cannot invalidate sample completion sentinel $path: $!\n";
	return 1;
}

1;
