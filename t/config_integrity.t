use strict;
use warnings;

use File::Find ();
use File::Spec;
use FindBin qw($Bin);
use Test::More;

# Behavioural checks on the shipped configuration: every lookup the code makes
# must resolve, and every repository file the configuration points to must
# exist. This catches dangling entries after scripts are renamed or retired.

my $root = File::Spec->rel2abs(File::Spec->catdir($Bin, '..'));
my @config_files = map { File::Spec->catfile($root, 'Mods', $_) }
	qw(config_internal.txt config_DBs.txt config.old);

# Large compiled/binary assets that can be absent from lightweight checkouts.
sub is_binary_asset {
	my ($relative) = @_;
	return $relative =~ m{^bin/(?!.*\.(?:pl|pm|sh)$)} || $relative =~ /\.(?:hmm|jar)$/;
}

my (%keys_by_file, %values);
for my $config (@config_files) {
	open my $fh, '<', $config or die "Cannot read $config: $!";
	my $file = (File::Spec->splitpath($config))[2];
	while (my $line = <$fh>) {
		next if $line =~ /^\s*(?:#|$)/;
		chomp $line;
		my ($key, $value) = split /\t/, $line, 3;
		next unless defined $key && length $key;
		$key =~ s/\s+$//;
		push @{ $keys_by_file{$file}{$key} }, $value // '';
		$values{$key} //= $value // '';
	}
	close $fh;
}

for my $file (sort keys %keys_by_file) {
	my @duplicates = grep { @{ $keys_by_file{$file}{$_} } > 1 } sort keys %{ $keys_by_file{$file} };
	is_deeply(\@duplicates, [], "$file defines every key once")
		or diag("duplicated keys: @duplicates");
}

# Every [MFLRDir]-relative path must exist in the repository.
my @missing_paths;
my $checked_paths = 0;
for my $key (sort keys %values) {
	my $value = $values{$key};
	while ($value =~ m{\[MFLRDir\]/+([^\s\t#]+)}g) {
		my $relative = $1;
		$relative =~ s{/+$}{};
		$relative =~ s{/\./}{/}g;
		next if $relative =~ m{^(?:data/DBs|gits)(?:/|$)}; # installed by helpers/install
		my $path = File::Spec->catfile($root, split m{/}, $relative);
		next if !-e $path && is_binary_asset($relative);
		$checked_paths++;
		push @missing_paths, "$key -> $relative" unless -e $path;
	}
}
ok($checked_paths > 50, "configuration references were found ($checked_paths checked)");
is_deeply(\@missing_paths, [], 'every repository path named in the configuration exists')
	or diag(join("\n", @missing_paths));

# Every literal required getProgPaths() lookup in pipeline code must resolve.
my @sources;
File::Find::find({ no_chdir => 1, wanted => sub {
	my $path = $File::Find::name;
	if (-d $path && $path =~ m{/(?:\.git|t|docs|doc|examples)$}) { $File::Find::prune = 1; return; }
	push @sources, $path if -f $path && $path =~ /\.(?:pl|pm)$/;
} }, $root);

my %unresolved;
for my $source (@sources) {
	open my $fh, '<', $source or die "Cannot read $source: $!";
	while (my $line = <$fh>) {
		next if $line =~ /^\s*#/;
		$line =~ s/\s#.*$//;
		while ($line =~ /getProgPaths\(\s*["']([^"'\$]+)["']\s*(,\s*0\s*)?\)/g) {
			my ($key, $optional) = ($1, $2);
			# Keys may be present with an empty value as a site-configurable placeholder.
			next if $optional || exists $values{$key};
			push @{ $unresolved{$key} }, File::Spec->abs2rel($source, $root) . ":$.";
		}
	}
	close $fh;
}
is_deeply([sort keys %unresolved], [], 'every required getProgPaths() key is configured')
	or diag(join("\n", map { "$_ <- @{ $unresolved{$_} }" } sort keys %unresolved));

done_testing();
