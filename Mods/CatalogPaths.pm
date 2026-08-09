package Mods::CatalogPaths;

use warnings;
use strict;

use Cwd qw(abs_path);
use Digest::SHA qw(sha256_hex);
use Exporter qw(import);
use Fcntl qw(:DEFAULT :flock);
use Errno qw(EEXIST);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use File::Spec;

our @EXPORT_OK = qw(
	catalog_identity
	catalog_identity_file
	catalog_map_manifest
	catalog_map_specs_match
	resolve_catalog_maps
	write_catalog_maps
);

sub _catalog_root {
	my ($catalog_dir) = @_;
	die "A gene-catalog directory is required\n"
		unless defined($catalog_dir) && $catalog_dir =~ /\S/;
	$catalog_dir =~ s{[\\/]+$}{};
	my $absolute = abs_path($catalog_dir);
	return defined($absolute) ? $absolute : File::Spec->rel2abs($catalog_dir);
}

sub catalog_identity_file {
	my ($catalog_dir) = @_;
	return File::Spec->catfile(_catalog_root($catalog_dir), 'LOGandSUB', 'catalog.sha256');
}

sub catalog_map_manifest {
	my ($catalog_dir) = @_;
	return File::Spec->catfile(_catalog_root($catalog_dir), 'LOGandSUB', 'inmap.txt');
}

sub catalog_identity {
	my ($catalog_dir) = @_;
	my $root = _catalog_root($catalog_dir);
	my $identity_file = catalog_identity_file($root);
	my $log_dir = dirname($identity_file);
	make_path($log_dir) unless -d $log_dir;

	if (-e $identity_file) {
		open my $identity_fh, '<', $identity_file
			or die "Cannot read catalog identity $identity_file: $!\n";
		flock($identity_fh, LOCK_SH)
			or die "Cannot lock catalog identity $identity_file: $!\n";
		my $identity = <$identity_fh> // '';
		close $identity_fh
			or die "Cannot close catalog identity $identity_file: $!\n";
		$identity =~ s/\s+//g;
		die "Invalid catalog identity in $identity_file\n"
			unless $identity =~ /\A[0-9a-f]{64}\z/;
		return $identity;
	}

	# Generate once and store locally. Copying the catalog together with
	# LOGandSUB preserves its identity, while rebuilding a removed catalog
	# at the same path receives a new identity.
	my $identity = sha256_hex(join("\0",
		'MATAFILER4 gene catalog', $root, time, $$, rand(), {}));
	my $identity_fh;
	if (!sysopen($identity_fh, $identity_file, O_WRONLY | O_CREAT | O_EXCL)) {
		return catalog_identity($root) if $! == EEXIST;
		die "Cannot create catalog identity $identity_file: $!\n";
	}
	flock($identity_fh, LOCK_EX)
		or die "Cannot lock catalog identity $identity_file: $!\n";
	print {$identity_fh} "$identity\n"
		or die "Cannot write catalog identity $identity_file: $!\n";
	close $identity_fh
		or die "Cannot close catalog identity $identity_file: $!\n";
	return $identity;
}

sub _normalise_map_paths {
	my ($catalog_dir, @entries) = @_;
	my $manifest_dir = dirname(catalog_map_manifest($catalog_dir));
	my @maps;
	for my $entry (@entries) {
		next unless defined $entry;
		for my $map (split /,/, $entry) {
			$map =~ s/^\s+|\s+$//g;
			next if $map eq '';
			$map = File::Spec->catfile($manifest_dir, $map)
				unless File::Spec->file_name_is_absolute($map);
			my $absolute = abs_path($map);
			$map = defined($absolute) ? $absolute : File::Spec->rel2abs($map);
			die "Catalog mapping file does not exist: $map\n" unless -f $map;
			push @maps, $map;
		}
	}
	die "Catalog map manifest contains no mapping files\n" unless @maps;
	return @maps;
}

sub catalog_map_specs_match {
	my ($left_spec, $right_spec) = @_;
	return 0 unless defined($left_spec) && defined($right_spec);

	my @left = _canonical_map_spec($left_spec);
	my @right = _canonical_map_spec($right_spec);
	return 0 unless @left && @right;
	return 0 unless @left == @right;
	for my $index (0 .. $#left) {
		return 0 unless _map_entries_match($left[$index], $right[$index]);
	}
	return 1;
}

sub _map_entries_match {
	my ($left, $right) = @_;
	return 1 if $left eq $right;
	return 0 unless -f $left && -f $right;
	return 0 unless -s $left == -s $right;
	return _map_file_digest($left) eq _map_file_digest($right);
}

sub _map_file_digest {
	my ($path) = @_;
	open my $fh, "<:raw", $path
		or die "Cannot read mapping file $path: $!\n";
	my $digest = Digest::SHA->new(256);
	my $buffer;
	while (read($fh, $buffer, 1024 * 1024)) {
		$digest->add($buffer);
	}
	close $fh or die "Cannot close mapping file $path: $!\n";
	return $digest->hexdigest;
}

sub _canonical_map_spec {
	my ($spec) = @_;
	my @maps;
	for my $map (split /,\s*|\r?\n/, $spec, -1) {
		$map =~ s/^\s+|\s+$//g;
		return () if $map eq '';
		my $absolute = abs_path($map);
		push @maps, defined($absolute) ? $absolute : File::Spec->rel2abs($map);
	}
	return @maps;
}

sub _write_compatibility_map_log {
	my ($manifest, @maps) = @_;
	my $legacy_log = File::Spec->catfile(dirname($manifest), 'GCmaps.inf');
	my $legacy_temporary = "$legacy_log.$$";
	open my $legacy_fh, '>', $legacy_temporary
		or die "Cannot write compatibility map log $legacy_temporary: $!\n";
	print {$legacy_fh} join(',', @maps), "\n"
		or die "Cannot write compatibility map log $legacy_temporary: $!\n";
	close $legacy_fh
		or die "Cannot close compatibility map log $legacy_temporary: $!\n";
	rename($legacy_temporary, $legacy_log)
		or die "Cannot publish compatibility map log $legacy_log: $!\n";
}

sub resolve_catalog_maps {
	my ($catalog_dir) = @_;
	my $manifest = catalog_map_manifest($catalog_dir);
	if (!-f $manifest) {
		my $log_dir = dirname($manifest);
		my @copied_maps = sort {
			my ($an) = $a =~ /map\.(\d+)\.txt\z/;
			my ($bn) = $b =~ /map\.(\d+)\.txt\z/;
			($an // 0) <=> ($bn // 0);
		} grep { -f $_ } glob(File::Spec->catfile($log_dir, 'map.*.txt'));
		return write_catalog_maps($catalog_dir, \@copied_maps) if @copied_maps;
		die "Catalog map file is missing: $manifest\n";
	}

	open my $map_fh, '<', $manifest
		or die "Cannot read catalog map file $manifest: $!\n";
	my @lines = <$map_fh>;
	close $map_fh or die "Cannot close catalog map file $manifest: $!\n";
	chomp @lines;
	my @nonempty = grep { /\S/ } @lines;
	die "Catalog map file is empty: $manifest\n" unless @nonempty;

	# New catalogs use inmap.txt as a manifest (one copied map per line).
	# Older single-map catalogs stored the mapping table itself at this path.
	my $looks_like_mapping_table = grep { /^\s*#/ || /\t/ } @nonempty;
	if ($looks_like_mapping_table) {
		_write_compatibility_map_log($manifest, $manifest);
		return $manifest;
	}

	my @maps = _normalise_map_paths($catalog_dir, @nonempty);
	_write_compatibility_map_log($manifest, @maps);
	return join(',', @maps);
}

sub write_catalog_maps {
	my ($catalog_dir, $maps) = @_;
	die "write_catalog_maps expects an array reference\n" unless ref($maps) eq 'ARRAY';
	my @normalised = _normalise_map_paths($catalog_dir, @$maps);
	my $manifest = catalog_map_manifest($catalog_dir);
	my $manifest_dir = dirname($manifest);
	make_path($manifest_dir) unless -d $manifest_dir;
	my $temporary = "$manifest.$$";
	open my $map_fh, '>', $temporary
		or die "Cannot write catalog map manifest $temporary: $!\n";
	for my $map (@normalised) {
		my $entry = dirname($map) eq $manifest_dir ? basename($map) : $map;
		print {$map_fh} "$entry\n";
	}
	close $map_fh
		or die "Cannot close catalog map manifest $temporary: $!\n";
	rename($temporary, $manifest)
		or die "Cannot publish catalog map manifest $manifest: $!\n";

	# Write-only compatibility/log artifact. Current code must resolve maps
	# through inmap.txt and never depend on this legacy single-line file.
	_write_compatibility_map_log($manifest, @normalised);
	return join(',', @normalised);
}

1;
