#!/usr/bin/env perl
use strict;
use warnings;

use Cwd qw(abs_path);
use Digest::SHA qw(sha1_hex);
use File::Glob qw(bsd_glob);
use File::Path qw(make_path remove_tree);
use File::Spec;
use Getopt::Long qw(GetOptions);

my %opt;
my @members;
my @member_locks;
GetOptions(
	'sample=s'        => \$opt{sample},
	'member=s@'       => \@members,
	'member-lock=s@'  => \@member_locks,
	'state-dir=s'     => \$opt{state_dir},
	'allowed-root=s'  => \$opt{allowed_root},
	'mapping-dir=s'   => \$opt{mapping_dir},
	'sample-temp=s'   => \$opt{sample_temp},
	'scratch-root=s'  => \$opt{scratch_root},
	'assembly=s'      => \$opt{assembly},
	'snp-log-dir=s'   => \$opt{snp_log_dir},
) or die "invalid cleanup arguments\n";

for my $required (qw(sample state_dir allowed_root mapping_dir sample_temp scratch_root assembly snp_log_dir)) {
	(my $flag = $required) =~ s/_/-/g;
	die "--$flag is required\n" unless defined($opt{$required}) && length($opt{$required});
}

sub existing_absolute {
	my ($path, $label) = @_;
	my $absolute = abs_path($path);
	die "$label does not exist: $path\n" unless defined($absolute);
	return File::Spec->canonpath($absolute);
}

sub prospective_absolute {
	my ($path) = @_;
	return File::Spec->canonpath(File::Spec->rel2abs($path));
}

sub is_below {
	my ($path, $root) = @_;
	my $relative = File::Spec->abs2rel($path, $root);
	return $relative ne File::Spec->updir()
		&& $relative !~ m{^\.\.(?:[\\/]|$)}
		&& !File::Spec->file_name_is_absolute($relative);
}

sub require_below {
	my ($path, $root, $label) = @_;
	die "$label escapes its allowed root: $path (root $root)\n"
		unless is_below($path, $root) && $path ne $root;
}

sub marker_for {
	my ($state_dir, $sample) = @_;
	return File::Spec->catfile($state_dir, sha1_hex($sample).'.complete');
}

sub marker_matches_assembly {
	my ($marker_path, $signature) = @_;
	return 0 unless -s $marker_path;
	open my $fh, '<', $marker_path or return 0;
	my $line = <$fh>;
	close $fh;
	return 0 unless defined($line);
	chomp $line;
	my (undef, $stored_signature) = split /\t/, $line, 2;
	return defined($stored_signature) && $stored_signature eq $signature;
}

sub remove_files {
	my (@files) = @_;
	my %seen;
	for my $file (@files) {
		next unless defined($file) && !$seen{$file}++ && -f $file;
		unlink $file or die "can't remove $file: $!\n";
		print "Removed temporary file $file\n";
	}
}

my $allowed_root = existing_absolute($opt{allowed_root}, 'allowed root');
my $state_dir = prospective_absolute($opt{state_dir});
require_below($state_dir, $allowed_root, 'cleanup state directory');
my $marker = marker_for($state_dir, $opt{sample});
@members = ($opt{sample}) unless @members;

my $mapping_dir = prospective_absolute($opt{mapping_dir});
require_below($mapping_dir, $allowed_root, 'mapping directory');
my @alignment_indexes;
for my $stem ("$opt{sample}-smd", "$opt{sample}.sup-smd") {
	push @alignment_indexes,
		File::Spec->catfile($mapping_dir, "$stem.bam.bai"),
		File::Spec->catfile($mapping_dir, "$stem.bam.csi"),
		File::Spec->catfile($mapping_dir, "$stem.cram.crai"),
		File::Spec->catfile($mapping_dir, "$stem.cram.csi");
}
remove_files(@alignment_indexes);

my $snp_log_dir = prospective_absolute($opt{snp_log_dir});
require_below($snp_log_dir, $allowed_root, 'SNP log directory');
remove_files(bsd_glob(File::Spec->catfile($snp_log_dir, '*.bed')));

my $scratch_root = existing_absolute($opt{scratch_root}, 'scratch root');
my $sample_temp = prospective_absolute($opt{sample_temp});
require_below($sample_temp, $scratch_root, 'sample temporary directory');
if (-d $sample_temp) {
	remove_tree($sample_temp, {error => \my $errors});
	if (@{$errors}) {
		my @messages;
		for my $entry (@{$errors}) {
			my ($path, $message) = %{$entry};
			push @messages, "$path: $message";
		}
		die "failed to remove sample temporary directory: ".join('; ', @messages)."\n";
	}
	print "Removed sample temporary directory $sample_temp\n";
}

my $assembly = existing_absolute($opt{assembly}, 'assembly');
my @assembly_stat = stat($assembly);
my $assembly_signature = join(':', @assembly_stat[0, 1, 7, 9, 10]);
make_path($state_dir) unless -d $state_dir;
open my $marker_fh, '>', $marker or die "can't write $marker: $!\n";
print {$marker_fh} "$opt{sample}\t$assembly_signature\n";
close $marker_fh or die "can't close $marker: $!\n";

my @missing_members = grep {
	!marker_matches_assembly(marker_for($state_dir, $_), $assembly_signature)
} @members;
my @active_locks;
for my $lock (@member_locks) {
	my $absolute_lock = prospective_absolute($lock);
	require_below($absolute_lock, $allowed_root, 'assembly-group sample lock');
	push @active_locks, $absolute_lock if -e $absolute_lock;
}
if (@missing_members || @active_locks) {
	print "Retaining assembly indexes; waiting for ".scalar(@missing_members)
		." group completion marker(s) and ".scalar(@active_locks)." active job lock(s)\n";
	exit 0;
}

if (!is_below($assembly, $allowed_root)) {
	print "Retaining mapper indexes for external assembly $assembly\n";
	exit 0;
}

my @assembly_indexes = (
	"$assembly.fai",
	"$assembly.gzi",
	"$assembly.mmi",
	(map { "$assembly.$_" } qw(amb ann bwt pac sa 0123 bwt.2bit.64)),
	bsd_glob("$assembly.bw2.*.bt2"),
	bsd_glob("$assembly.bw2.*.bt2l"),
	bsd_glob("$assembly.bw2.*.ebwt"),
	bsd_glob("$assembly.bw2.*.ebwtl"),
	bsd_glob("$assembly.kma*"),
);
remove_files(@assembly_indexes);
print "All assembly-group members complete; mapper and FASTA indexes are clean\n";
