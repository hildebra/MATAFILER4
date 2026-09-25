#!/usr/bin/env perl
#
# Produce a version of GTDB compatible with MATAFILER4.
#
# Perl port of the former helpers/install/get_gtdb.py (Anthony Duncan, script
# version 0.1). It downloads a given GTDB release, formats the marker genes and
# taxonomy tables into the structure MATAFILER expects, optionally extracts the
# GTDB-Tk reference data, and can update Mods/config_DBs.txt to use the result.
#
# Subcommands: download, extract, configure, all (run with -h for details).
#
# Adding a GTDB release works as in the Python version: subclass GTDBVersion
# (a package below), override new() for different URLs, _extract_markers for a
# different marker archive layout or _format_taxonomy for different metadata,
# and add the dispatch to GTDBVersion::from_version().
#
# Large archives (tens of GB) are never read into memory: marker and GTDB-Tk
# archives are streamed through the system tar (with pigz when available),
# metadata tables through gzip/pigz, and every exit status is checked.
# Only core Perl modules are used (perl >= 5.32 in the MF4 environment).
#
# Intentional deviations from get_gtdb.py (each is also marked "DEVIATION" at
# the place it applies):
#  1. The MATAFILER directory is taken from MF4DIR (exported by installer.sh),
#     then the legacy MGTKDIR, then the location of this script (../..).
#     The Python only read MGTKDIR and crashed (str(None)) when it was unset.
#  2. DBDir is looked up like Mods/IO_Tamoc_progs.pm does (config.txt, then
#     Mods/config_internal.txt, then Mods/config_DBs.txt; first definition
#     wins; $VAR expansion as in truePath()). The Python took the last DBDir in
#     config.txt only and did not expand the default "$MF4DIR/data/DBs/", so it
#     always fell back to absolute paths.
#  3. [DBDir]-relative config values are written as "[DBDir]/rel/path" instead
#     of the Python's "[DBDir]//rel/path"; updated lines are tagged
#     "#Updated by get_gtdb.pl". Keys that the user config (config.txt)
#     overrides are reported, since the pipeline gives config.txt precedence.
#  4. GTDBtk_mash handling is idempotent: an already commented line is not
#     commented again and the "#Updated by" tag is not appended repeatedly (the
#     Python added one '#' and one tag per run and also commented unrelated
#     lines that merely mentioned GTDBtk_mash).
#  5. GTDBtk_DB points to the release directory inside gtdb/ (e.g.
#     gtdb/release226) when there is exactly one, as GTDBTK_DATA_PATH requires;
#     GTDBtk_DB/GTDBtk_mash are left untouched when no GTDB-Tk data was
#     extracted (--tk skip) instead of being pointed at a missing directory.
#  6. The GTDB-Tk environment is MF4gtdbtk (installer.sh, config_internal.txt),
#     not MGTKgtdbtk; micromamba is taken from MAMBA_EXE when set.
#  7. GTDB-Tk compatibility table updated from the GTDB-Tk documentation
#     (checked 2026-09-25): R226 2.4.1-2.6.1, R220 2.4.0-2.6.1, R232 2.7.0+.
#     The Python capped R220/R226 at 2.5.2, so --install-tk would have
#     downgraded the gtdbtk=2.6 environment that installer.sh creates.
#  8. Output tables use LF line endings. Python's csv.writer wrote CRLF, which
#     the Perl consumers (MG_LCA.pl reformatGTDBtax, TamocFunc::readTabbed3)
#     do not strip. Content is otherwise byte-identical. Rows are split on tabs
#     without CSV quote interpretation (GTDB tables are not quoted).
#  9. --tk skip no longer crashes extraction (the string "skip" was truthy in
#     Python, so _extract_tk ran and died on a missing key). Split/full GTDB-Tk
#     handling follows the recorded --tk mode: "full" used to truncate the
#     downloaded full archive (cat of non-existent parts into it). Split parts
#     are streamed into tar instead of first concatenating ~110GB to disk.
# 10. Exactly one archaea-only marker file is now moved too (Python required
#     more than one, then failed on rmdir); an archive that yields no marker
#     files is an error rather than silently producing an empty database.
# 11. Failed-directory cleanup only runs for directory (split package)
#     downloads; for a failed single-file download the Python crashed inside
#     the cleanup (statistics.mode of an empty list), hiding the real error.
# 12. wget/curl are only required for "download"/"all"; the Python required
#     one of them even for "extract" and "configure" (offline systems).
#     wget uses --progress=dot:giga and curl -f -L, so log volume stays small
#     and HTTP errors or redirects are not saved as data; resume (-c / -C -)
#     and no-clobber behaviour are unchanged.
# 13. --test copies stubs selected by download label, from the get_gtdb/
#     directory next to this script (any working directory). The Python matched
#     "ar" in the URL, which also matches ".tar.gz", so the bacterial marker
#     stub was replaced by the archaeal one; its full-package stub was an empty
#     file, now it is the concatenated stub parts so extraction is testable.
# 14. Help, README.txt and the meta.json summary name get_gtdb.pl;
#     README.txt names Mods/config_DBs.txt (not config_DB.txt) and "gtdbtk"
#     (typo "gtbtk"), recommends the GTDB-Tk version --install-tk would
#     install, lists GTDBtk_mash only for GTDB-Tk < 2.5, and takes the
#     GTDB-Tk release subdirectory from gtdb/ (the Python picked an arbitrary
#     subdirectory of the output directory). --debug is accepted before or
#     after the subcommand.

use strict;
use warnings;

use Cwd qw(getcwd realpath);
use File::Basename qw(dirname basename);
use File::Copy qw(copy);
use File::Path qw(make_path);
use File::Spec;
use File::Temp ();
use Getopt::Long qw(GetOptionsFromArray);
use JSON::PP ();
use POSIX qw(strftime);
use Time::HiRes ();

our $VERSION = '0.1';          # script_version recorded in meta.json
our $SCRIPT_VERSION = 0.1;
our $DEBUG = 0;
our $PROG = 'get_gtdb.pl';
our $SCRIPT_DIR = dirname(File::Spec->rel2abs(__FILE__));
our $UPDATED_TAG = '#Updated by get_gtdb.pl';
our $GTDBTK_ENV = 'MF4gtdbtk';

###########################################################################
# Logging (format of the Python logger: 'LEVEL [dd/mm/YYYY HH:MM:SS]: msg')
###########################################################################

sub _log {
	my ($level, $fmt, @args) = @_;
	my $msg = @args ? sprintf($fmt, @args) : $fmt;
	print STDERR "$level [" . strftime('%d/%m/%Y %I:%M:%S', localtime) . "]: $msg\n";
}
sub log_info  { _log('INFO', @_) }
sub log_warn  { _log('WARNING', @_) }
sub log_error { _log('ERROR', @_) }
sub log_debug { _log('DEBUG', @_) if $DEBUG }

###########################################################################
# Small helpers
###########################################################################

# str(datetime.now()) as Python prints it
sub py_now {
	my ($s, $us) = Time::HiRes::gettimeofday();
	my $str = strftime('%Y-%m-%d %H:%M:%S', localtime($s));
	$str .= sprintf('.%06d', $us) if $us;
	return $str;
}

# str(Path(p)) normalisation: collapse '//', drop '.' parts and trailing '/'
sub norm_path {
	my ($p) = @_;
	return '.' if !defined($p) || $p eq '';
	my $abs = $p =~ m{^/};
	my @parts = grep { $_ ne '' && $_ ne '.' } split m{/+}, $p;
	my $r = join('/', @parts);
	return $abs ? "/$r" : ($r eq '' ? '.' : $r);
}

sub join_path { return norm_path(join('/', @_)) }

# Path.resolve() (non-strict): realpath of the deepest existing ancestor
sub resolve_path {
	my ($p) = @_;
	$p = File::Spec->rel2abs($p);
	my @rest;
	my $cur = $p;
	while (!-e $cur) {
		my $parent = dirname($cur);
		unshift @rest, basename($cur);
		last if $parent eq $cur;
		$cur = $parent;
	}
	my $real = realpath($cur);
	$real = $cur unless defined $real;
	my @parts;
	for my $seg (split(m{/+}, $real), @rest) {
		next if $seg eq '' || $seg eq '.';
		if ($seg eq '..') { pop @parts; next; }
		push @parts, $seg;
	}
	return '/' . join('/', @parts);
}

sub write_text {
	my ($file, $text) = @_;
	open my $fh, '>', $file or die "Cannot write $file: $!\n";
	print {$fh} $text or die "Cannot write $file: $!\n";
	close $fh or die "Cannot close $file: $!\n";
}

sub shell_quote {
	return join ' ', map { m{^[\w/.,:=+\@%^-]+$} ? $_ : "'" . s/'/'\\''/gr . "'" } @_;
}

sub program_avail {
	my ($cmd) = @_;
	for my $d (File::Spec->path) {
		my $p = File::Spec->catfile($d, $cmd);
		return 1 if -f $p && -x _;
	}
	return 0;
}

my %PIGZ;
sub have_pigz { $PIGZ{x} //= program_avail('pigz'); return $PIGZ{x} }

# Run a command without a shell. stdout/stderr go to temporary files (wget/tar
# logs are never held in memory while the command runs).
sub run_capture {
	my (@cmd) = @_;
	log_debug('Command: %s', shell_quote(@cmd));
	my $out = File::Temp->new(TMPDIR => 1);
	my $err = File::Temp->new(TMPDIR => 1);
	my $pid = fork();
	die "Cannot fork: $!\n" unless defined $pid;
	if ($pid == 0) {
		open STDIN, '<', File::Spec->devnull or POSIX::_exit(127);
		open STDOUT, '>&', $out or POSIX::_exit(127);
		open STDERR, '>&', $err or POSIX::_exit(127);
		exec { $cmd[0] } @cmd or do { print STDERR "Cannot execute $cmd[0]: $!\n"; POSIX::_exit(127) };
	}
	waitpid($pid, 0);
	my $status = $?;
	return ($status, _slurp_tail("$out"), _slurp_tail("$err"));
}

sub _slurp_tail {
	my ($f) = @_;
	my $max = 1 << 20;
	open my $fh, '<', $f or return '';
	binmode $fh;
	my $size = -s $fh;
	seek($fh, $size - $max, 0) if $size > $max;
	local $/;
	my $t = <$fh>;
	close $fh;
	return defined $t ? $t : '';
}

sub run_checked {
	my ($what, @cmd) = @_;
	my ($st, $out, $err) = run_capture(@cmd);
	if ($st != 0) {
		log_error('%s failed (exit status %d): %s', $what, $st >> 8, shell_quote(@cmd));
		log_error('Stderr: %s', $err) if length $err;
		log_error('Stdout: %s', $out) if length $out;
		die "$what failed\n";
	}
	return ($out, $err);
}

# Line reader over a gzip file, streamed through pigz/gzip (fallback: core
# IO::Uncompress::Gunzip). Multi-member files are read completely like zcat.
sub open_gz {
	my ($file) = @_;
	die "File not found: $file\n" unless -e $file;
	my $fh;
	if (have_pigz() || program_avail('gzip')) {
		my @cmd = (have_pigz() ? 'pigz' : 'gzip', '-dc', $file);
		log_debug('Command: %s', shell_quote(@cmd));
		open($fh, '-|', @cmd) or die "Cannot run $cmd[0] on $file: $!\n";
	} else {
		require IO::Uncompress::Gunzip;
		$fh = IO::Uncompress::Gunzip->new($file, MultiStream => 1)
			or die "Cannot open $file: "
				. do { no warnings 'once'; $IO::Uncompress::Gunzip::GunzipError } . "\n";
	}
	return $fh;
}

sub close_gz {
	my ($fh, $file) = @_;
	if (ref($fh) ne q{GLOB}) { $fh->close; return; }
	close($fh) or die "Decompression of $file failed (exit status " . ($? >> 8) . ")\n";
}

# Iterate over the rows of the byte concatenation of several gzip files, as the
# Python did with "cat a b > comb" (and "zcat x | tail -n +2" for skip_first).
# A last line without newline is joined with the next stream's first line.
sub each_concat_row {
	my ($sources, $cb) = @_;
	my $carry;
	for my $src (@$sources) {
		my ($file, $skip_first) = @$src;
		my $fh = open_gz($file);
		my $first = 1;
		while (defined(my $line = <$fh>)) {
			if ($first && $skip_first) { $first = 0; next; }
			$first = 0;
			if (defined $carry) { $line = $carry . $line; undef $carry; }
			if ($line !~ /\n\z/) { $carry = $line; next; }
			_row_cb($line, $cb);
		}
		close_gz($fh, $file);
	}
	_row_cb($carry, $cb) if defined $carry && length $carry;
}

sub _row_cb {
	my ($line, $cb) = @_;
	$line =~ s/\r?\n\z//;
	$line =~ s/\r\z//;
	return if $line eq '';   # csv readers skip empty rows
	$cb->([split /\t/, $line, -1]);
}

sub py_split_semicolon {
	my ($s) = @_;
	return ('') if $s eq '';
	return split /;/, $s, -1;
}

sub py_strip { my ($s) = @_; $s =~ s/^\s+//; $s =~ s/\s+\z//; return $s }

# csv.writer(delimiter="\t", quoting=QUOTE_NONE) without escapechar
sub row_quote_none {
	my @f = @_;
	for (@f) {
		die "Field needs escaping but quoting is disabled: $_\n" if /[\t"\r\n]/;
	}
	return join("\t", @f) . "\n";
}

# csv.writer(delimiter="\t") default QUOTE_MINIMAL
sub row_quote_minimal {
	my @f = map { /[\t"\r\n]/ ? '"' . s/"/""/gr . '"' : $_ } @_;
	return join("\t", @f) . "\n";
}

# Python float repr for version numbers (226.0, 207.2)
sub py_float {
	my ($v) = @_;
	return sprintf('%.1f', $v) if $v == int($v) && abs($v) < 1e16;
	for my $p (1 .. 17) {
		my $s = sprintf("%.${p}g", $v);
		return $s if $s + 0 == $v;
	}
	return "$v";
}

# json.dump(..., indent=4) with Python's default ensure_ascii escaping
sub py_json {
	my ($val, $ind) = @_;
	$ind //= 0;
	my $pad = ' ' x (4 * ($ind + 1));
	my $end = ' ' x (4 * $ind);
	if (ref($val) eq 'ARRAY') {  # ordered pairs => JSON object
		my @kv = @$val;
		return '{}' unless @kv;
		my @items;
		while (@kv) {
			my ($k, $v) = splice(@kv, 0, 2);
			push @items, $pad . py_json_str($k) . ': ' . py_json($v, $ind + 1);
		}
		return "{\n" . join(",\n", @items) . "\n$end}";
	}
	if (ref($val) eq 'SCALAR') { return $$val }   # raw number literal
	return 'null' unless defined $val;
	return py_json_str($val);
}

sub py_json_str {
	my ($s) = @_;
	my $u = $s;
	utf8::decode($u) unless utf8::is_utf8($u);
	my %esc = ('"' => '\\"', '\\' => '\\\\', "\n" => '\\n', "\r" => '\\r',
		"\t" => '\\t', "\b" => '\\b', "\f" => '\\f');
	my $o = '';
	for my $c (split //, $u) {
		if (exists $esc{$c}) { $o .= $esc{$c}; next; }
		my $n = ord($c);
		if ($n >= 0x20 && $n <= 0x7e) { $o .= $c; next; }
		if ($n > 0xFFFF) {
			$n -= 0x10000;
			$o .= sprintf('\\u%04x\\u%04x', 0xD800 + ($n >> 10), 0xDC00 + ($n & 0x3FF));
		} else {
			$o .= sprintf('\\u%04x', $n);
		}
	}
	return qq{"$o"};
}

sub version_as_tuple { return map { int($_) } split /\./, $_[0] }

sub version_ge {
	my @a = version_as_tuple($_[0]);
	my @b = version_as_tuple($_[1]);
	while (@a || @b) {
		return 1 if !@b;           # (2,5,2) >= (2,5)
		return 0 if !@a;
		my ($x, $y) = (shift @a, shift @b);
		return $x > $y if $x != $y;
	}
	return 1;
}

sub is_url_dir { return substr($_[0], -1) eq '/' }

sub mkdir_parents {
	my ($dest) = @_;
	log_debug('Making directory %s', $dest);
	my $target = basename($dest) =~ /\./ ? dirname($dest) : $dest;
	make_path($target) unless -d $target;
}

sub prompt_input_set {
	my ($prompt, $valid, $case) = @_;
	my @valid = map { $case ? lc($_) : $_ } @$valid;
	local $| = 1;
	while (1) {
		print "$prompt (" . join('/', @valid) . "): ";
		my $x = <STDIN>;
		if (!defined $x) {
			print "\n";
			log_warn('No answer received (end of input), assuming "no".');
			return 'no';
		}
		chomp $x;
		$x = lc $x if $case;
		my @m = grep { substr($_, 0, length $x) eq $x } @valid;
		if (@m == 1) {
			log_debug('User entered %s, selected %s', $x, $m[0]);
			return $m[0];
		}
		log_warn('%s not valid, must be one of [%s]', $x, join(', ', map {"'$_'"} @valid));
	}
}

sub greet {
	print "\n\033[96mget_gtdb\033[0m\n\033[96m" . ("\xe2\x94\x80" x 8) . "\033[0m\n";
}

###########################################################################
# MATAFILER installation and configuration lookup
###########################################################################

# DEVIATION 1: MF4DIR, then MGTKDIR, then this script's repository (../..).
sub mf4_dir_candidates {
	my @c;
	for my $env (qw(MF4DIR MGTKDIR)) {
		push @c, [$ENV{$env}, "environment variable $env"]
			if defined $ENV{$env} && length $ENV{$env};
	}
	push @c, [File::Spec->catdir($SCRIPT_DIR, '..', '..'), 'location of get_gtdb.pl'];
	return @c;
}

sub get_mf4_dir {
	for my $c (mf4_dir_candidates()) {
		my ($d, $src) = @$c;
		my $dn = $d;
		$dn =~ s{(?<=.)/+\z}{};
		if (-f File::Spec->catfile($dn, 'Mods', 'config_DBs.txt')) {
			log_debug('MATAFILER directory %s (from %s)', $dn, $src);
			return resolve_path($dn);
		}
		log_warn('%s (from %s) does not contain Mods/config_DBs.txt, ignoring it.', $d, $src);
	}
	return undef;
}

# Directory named in the README instructions (target system may differ).
sub instructions_mf4_dir {
	for my $env (qw(MF4DIR MGTKDIR)) {
		if (defined $ENV{$env} && length $ENV{$env}) {
			(my $d = $ENV{$env}) =~ s{(?<=.)/+\z}{};
			return $d;
		}
	}
	return '[Your MATAFILER Dir]';
}

# Read key -> value (second tab column) of a MATAFILER config; first wins as
# in IO_Tamoc_progs::loadConfigs.
sub parse_mgtk_config {
	my ($cfg) = @_;
	my %config;
	my @order;
	open my $fh, '<', $cfg or die "Cannot read $cfg: $!\n";
	while (my $line = <$fh>) {
		$line =~ s/\r?\n\z//;
		next if $line =~ /^\s*$/;
		next if $line =~ /^\s*#/;
		my @parts = split /\t/, $line, -1;
		next if exists $config{$parts[0]};
		$config{$parts[0]} = @parts > 1 ? $parts[1] : undef;
		push @order, $parts[0];
	}
	close $fh;
	return \%config;
}

# IO_Tamoc_progs::truePath(): expand $NAME / ${NAME} if the value starts with $
sub true_path {
	my ($p) = @_;
	return $p unless $p =~ /^\$/;
	$p =~ s{\$(?:\{([A-Za-z_][A-Za-z0-9_]*)\}|([A-Za-z_][A-Za-z0-9_]*))}{
		my $n = defined $1 ? $1 : $2;
		die "Environment variable \$$n used in DBDir '$_[0]' is not set\n"
			unless defined $ENV{$n} && $ENV{$n} ne '';
		$ENV{$n};
	}gex;
	return $p;
}

# DEVIATION 2: DBDir as the pipeline resolves it.
sub find_dbdir {
	my ($mf4) = @_;
	my @files = (File::Spec->catfile($mf4, 'config.txt'));
	push @files, File::Spec->catfile($mf4, 'Mods', 'MATAFILERcfg.txt') unless -e $files[0];
	push @files, map { File::Spec->catfile($mf4, 'Mods', $_) } qw(config_internal.txt config_DBs.txt);
	for my $f (@files) {
		next unless -e $f;
		log_debug('Look for DBDir in %s', $f);
		my $cfg = parse_mgtk_config($f);
		next unless exists $cfg->{DBDir} && defined $cfg->{DBDir};
		my $v = $cfg->{DBDir};
		$v =~ s/#.*//;
		$v = py_strip($v);
		next if $v eq '';
		my $exp = eval { true_path($v) };
		if (!defined $exp) {
			log_warn('Cannot expand DBDir "%s" from %s: %s', $v, $f, $@ =~ s/\n\z//r);
			return (undef, $v, $f);
		}
		return ($exp, $v, $f);
	}
	return (undef, undef, undef);
}

###########################################################################
# Downloaders
###########################################################################

sub cmd_wget {
	my ($src, $dest) = @_;
	mkdir_parents($dest);
	if (is_url_dir($src)) {
		return ('wget', '-r', '-np', '-nH', '-nd', '-R', 'index.html*', '--level=1',
			'-nc', '-c', '--progress=dot:giga', $src, '-P', $dest);
	}
	return ('wget', '-c', '--progress=dot:giga', $src, '-O', $dest);
}

sub cmd_curl {
	my ($src, $dest) = @_;
	mkdir_parents($dest);
	if (is_url_dir($src)) {
		die "Directory download is not supported by curl. "
			. "Please install wget, or download full GTDBtk package (--tk full).\n";
	}
	return ('curl', '-f', '-L', '-o', $dest, '-C', '-', $src);
}

# DEVIATION 13: stub selection by label, stubs next to this script.
sub fetch_test {
	my ($label, $src, $dest) = @_;
	mkdir_parents($dest);
	my $root = File::Spec->catdir($SCRIPT_DIR, 'get_gtdb');
	my @parts = map { "gtdbtk_dummy.tar.gz.part_$_" } qw(aa ab ac);
	my %stub = map { $_ => "$_.tar.gz" } qw(bac_markers arc_markers);
	$stub{$_} = "$_.tsv.gz" for qw(bac_taxonomy arc_taxonomy bac_metadata arc_metadata);
	my $cp = sub {
		my ($from, $to) = @_;
		copy($from, $to) or die "Cannot copy $from to $to: $!\n";
	};
	if (exists $stub{$label}) {
		$cp->(File::Spec->catfile($root, $stub{$label}), $dest);
	} elsif ($label eq 'tk_database_parts') {
		$cp->(File::Spec->catfile($root, $_), File::Spec->catfile($dest, $_)) for @parts;
	} elsif ($label eq 'tk_database') {
		open my $o, '>', $dest or die "Cannot write $dest: $!\n";
		binmode $o;
		for my $p (@parts) {
			open my $i, '<', File::Spec->catfile($root, $p) or die "Cannot read stub $p: $!\n";
			binmode $i;
			local $/;
			print {$o} scalar(<$i>);
			close $i;
		}
		close $o or die "Cannot write $dest: $!\n";
	} else {
		log_warn('Defaulting to empty file for %s in %s', $src, $dest);
		write_text($dest, '');
	}
	return (0, '', '');
}

###########################################################################
package GTDBVersion;
###########################################################################

use File::Path qw(make_path);

BEGIN { *log_info = \&main::log_info; *log_warn = \&main::log_warn;
	*log_error = \&main::log_error; *log_debug = \&main::log_debug; }

# DEVIATION 7: compatibility table from the GTDB-Tk documentation
# (https://ecogenomics.github.io/GTDBTk/installing/index.html, 2026-09-25).
# The Python had 226/220 => 2.4.0-2.5.2 and no 232.
our @GTDBTK_VERSIONS = (
	[232.0, '2.7.0', 'Current'],
	[226.0, '2.4.1', '2.6.1'],
	[220.0, '2.4.0', '2.6.1'],
	[214.0, '2.1.0', '2.3.2'],
	[207.2, '2.1.0', '2.3.2'],
	[207.0, '2.0.0', '2.0.0'],
	[202.0, '1.5.0', '1.7.0'],
);

our $DEFAULT_BASE_URL = 'https://data.ace.uq.edu.au/public/gtdb/data/releases/';

sub new {
	my ($class, $version, %kw) = @_;
	my $self = bless {
		version  => $version + 0,
		base_url => $kw{base_url} // $DEFAULT_BASE_URL,
		test     => $kw{test} ? 1 : 0,
	}, $class;
	$self->{required_urls} = $kw{required_urls} // $self->default_urls;
	if ($self->{test}) {
		log_warn('Using test mode for downloads. This will create dummy files '
			. 'rather than performing downloads.');
	}
	return $self;
}

sub v0 { sprintf('%.0f', $_[0]{version}) }
sub v1 { sprintf('%.1f', $_[0]{version}) }

# Ordered label => url pairs (dict insertion order of the Python)
sub default_urls {
	my ($self) = @_;
	my ($v0, $v1) = ($self->v0, $self->v1);
	my $dir = "$self->{base_url}/release$v0/$v1";
	return [
		bac_markers  => "$dir/genomic_files_all/bac120_marker_genes_all_r$v0.tar.gz",
		bac_taxonomy => "$dir/bac120_taxonomy_r$v0.tsv.gz",
		bac_metadata => "$dir/bac120_metadata_r$v0.tsv.gz",
		arc_markers  => "$dir/genomic_files_all/ar53_marker_genes_all_r$v0.tar.gz",
		arc_taxonomy => "$dir/ar53_taxonomy_r$v0.tsv.gz",
		arc_metadata => "$dir/ar53_metadata_r$v0.tsv.gz",
		tk_database  => "$dir/auxillary_files/gtdbtk_package/full_package/gtdbtk_r${v0}_data.tar.gz",
		tk_database_parts => "$dir/auxillary_files/gtdbtk_package/split_package/",
	];
}

sub url_pairs {
	my ($self) = @_;
	my @kv = @{ $self->{required_urls} };
	my @pairs;
	push @pairs, [splice(@kv, 0, 2)] while @kv;
	return @pairs;
}

sub local_name {
	my ($label, $url) = @_;
	return main::is_url_dir($url) ? $label : (split m{/}, $url)[-1];
}

# DEVIATION 12: the download tool is only resolved when downloading.
sub _downloader {
	my ($self) = @_;
	return $self->{dl} if $self->{dl};
	if ($self->{test}) {
		$self->{dl} = 'test';
	} elsif (main::program_avail('wget')) {
		$self->{dl} = 'wget';
	} elsif (main::program_avail('curl')) {
		$self->{dl} = 'curl';
		log_warn('Using curl for downloads. This does not support directory '
			. 'downloads so will not be able to download GTDBtk database in parts.');
	} else {
		log_error('wget or curl must be available to perform downloads');
		die "Missing wget or curl\n";
	}
	return $self->{dl};
}

sub _fetch {
	my ($self, $label, $url, $local) = @_;
	my $dl = $self->_downloader;
	return main::fetch_test($label, $url, $local) if $dl eq 'test';
	my @cmd = $dl eq 'wget' ? main::cmd_wget($url, $local) : main::cmd_curl($url, $local);
	log_debug("Download %s\nSource: %s\nDestination: %s\nCommand: %s",
		$label, $url, $local, main::shell_quote(@cmd));
	return main::run_capture(@cmd);
}

sub download {
	my ($self, %a) = @_;
	my $dl_dir = $a{download_dir};
	my $offline = $a{offline} ? 1 : 0;
	my $tk = $a{download_tk};    # 'split', 'full', or false/'skip'
	my $dl_check = main::join_path($dl_dir, '.download.finished');
	if (-e $dl_check) {
		log_info('Download completion marker found (%s). If you want to redo '
			. 'downloads delete this file.', $dl_check);
		return { map { $_->[0] => main::join_path($dl_dir, local_name(@$_)) } $self->url_pairs };
	}

	make_path($dl_dir) if !$offline && !-d $dl_dir;
	log_info('Downloads required:');
	my %tk_remove = map { $_ => 1 } qw(tk_database tk_database_parts);
	if (defined $tk && $tk eq 'full') { %tk_remove = (tk_database_parts => 1) }
	elsif (defined $tk && $tk eq 'split') { %tk_remove = (tk_database => 1) }
	my @filt = grep { !$tk_remove{ $_->[0] } } $self->url_pairs;
	log_info("* %s\t%s", @$_) for @filt;

	my %downloaded;
	for my $p (@filt) {
		my ($label, $url) = @$p;
		my $local = main::join_path($dl_dir, local_name($label, $url));
		die "File exists: $local\n" if -e $local && $a{error_existing};
		if ($offline) {
			if (!-e $local) {
				log_error('Extract mode, but %s not found', $local);
				die "File not found: $local\n";
			}
			log_info('Extract mode, download skipped for %s', $label);
			$downloaded{$label} = $local;
			next;
		}
		log_info('Download %s', $label);
		my ($st, $out, $err) = $self->_fetch($label, $url, $local);
		if ($st != 0 || !-e $local) {
			# wget reports an already complete file with "fully retrieved"
			if ($err !~ /fully retrieved/) {
				log_error('Download of %s failed. Dumping output.', $url);
				log_error('Stderr');
				log_error('%s', $err);
				log_error('Stdout');
				log_error('%s', $out);
				# DEVIATION 11: only directory downloads cannot be resumed
				$self->_clean_failed_gtdbtk_dir($local) if main::is_url_dir($url) && -d $local;
				die "Download of $url failed (exit status " . ($st >> 8) . ")\n";
			}
		}
		$downloaded{$label} = $local;
	}
	if (!$offline) {
		$self->_processing_metadata($dl_dir,
			dir => main::norm_path($dl_dir), tk => main::py_str_tk($tk));
		main::write_text($dl_check, main::py_now());
	}
	return \%downloaded;
}

# Remove the last, incomplete part of a failed split-package download.
# wget -nc does not resume these, so only parts known to be complete are kept.
sub _clean_failed_gtdbtk_dir {
	my ($self, $dir) = @_;
	log_info('Directory download failed. Removing incomplete file '
		. 'as cannot resume directory download.');
	my @gzs = sort grep { -f $_ } glob(main::join_path($dir, '*.tar.gz*'));
	return unless @gzs;
	my (%count, @sizes);
	for my $g (@gzs) { my $s = -s $g; push @sizes, $s unless $count{$s}++; }
	my $mode_size = $sizes[0];
	for my $s (@sizes) { $mode_size = $s if $count{$s} > $count{$mode_size}; }
	my $last = $gzs[0];
	for my $g (@gzs) { $last = $g if (stat $g)[9] > (stat $last)[9]; }
	log_debug('Last file: %s', $last);
	log_debug('Mode size: %s', $mode_size);
	if ((-s $last) != $mode_size || @gzs < 3) {
		log_info('Deleting %s', $last);
		unlink $last or die "Cannot delete $last: $!\n";
	}
}

sub _tar_decompress_opt {
	my ($src) = @_;
	return ('--use-compress-program=pigz') if $src =~ /\.gz\z/i && main::have_pigz();
	return ();
}

# Extract all FAA/FNA of a marker tarball flat into dest (tar streams it).
sub _extract_seqs {
	my ($self, $src, $dest) = @_;
	log_debug('Source %s', $src);
	log_debug('Destination: %s', $dest);
	make_path($dest) unless -d $dest;
	die "Marker archive not found: $src\n" unless -f $src;
	# --strip-components 2 flattens <archive>/{faa,fna}/<marker>
	main::run_checked("Extraction of $src", 'tar', _tar_decompress_opt($src), '-xf', $src,
		'--strip-components', '2', '-C', $dest);
	log_debug('Extraction complete');
	my @f = ((sort glob(main::join_path($dest, '*.fna'))), (sort glob(main::join_path($dest, '*.faa'))));
	return \@f;
}

sub _extract_markers {
	my ($self, $dest, $files) = @_;
	my $finished = main::join_path($dest, '.marker.finished');
	log_debug('Check for finished marker file: %s', $finished);
	if (-e $finished) {
		log_info('Marker extraction already complete.');
		return;
	}
	log_info('Extracting archaeal markers');
	my $arc_dest = main::join_path($dest, 'arc');
	make_path($arc_dest) unless -d $arc_dest;
	my $arc = $self->_extract_seqs($files->{arc_markers}, $arc_dest);
	log_info('Extracting bacterial markers');
	my $bac_dest = $dest;
	my $bac = $self->_extract_seqs($files->{bac_markers}, $bac_dest);
	# DEVIATION 10: an archive without marker sequences is an error
	die "No marker sequences (*.faa/*.fna) extracted from $files->{arc_markers}\n" unless @$arc;
	die "No marker sequences (*.faa/*.fna) extracted from $files->{bac_markers}\n" unless @$bac;

	log_info('Combining archaeal and bacterial markers');
	my (%by_name, @names);
	for my $f (@$bac, @$arc) {
		my $n = main::basename($f);
		push @names, $n unless $by_name{$n};
		push @{ $by_name{$n} }, $f;
	}
	for my $n (grep { @{ $by_name{$_} } > 1 } @names) {
		my ($bac_f, $arc_f) = @{ $by_name{$n} }[0, 1];
		log_debug('Concatenate marker %s', $n);
		main::append_file($arc_f, $bac_f);
		unlink $arc_f or die "Cannot remove $arc_f: $!\n";
	}
	# DEVIATION 10: move remaining archaea-specific markers (Python: only if >1)
	my @rest = sort grep { -e $_ } glob(main::join_path($arc_dest, '*.f*'));
	if (@rest) {
		log_info('Moving remaining archaea specific markers');
		for my $f (@rest) {
			my $to = main::join_path($bac_dest, main::basename($f));
			rename($f, $to) or die "Cannot move $f to $to: $!\n";
		}
	}
	rmdir $arc_dest or die "Cannot remove $arc_dest (not empty?): $!\n";
	main::write_text($finished, main::py_now());
}

# DEVIATION 9: mode taken from --tk; split parts are streamed into tar.
sub _extract_tk {
	my ($self, $dest, $files, $mode) = @_;
	my $finished = main::join_path($dest, '.tk.finished');
	if (-e $finished) {
		log_info('GTDBtk extraction already completed.');
		return;
	}
	make_path($dest) unless -d $dest;
	if ($mode eq 'split') {
		my $part_src = $files->{tk_database_parts}
			// die "GTDB-Tk split package was not downloaded\n";
		my @parts = sort grep { -f $_ } glob(main::join_path($part_src, '*.tar.gz.part*'));
		die "No GTDB-Tk database parts (*.tar.gz.part*) found in $part_src\n" unless @parts;
		log_info('GTDBtk database parts (%d) will be streamed into tar for '
			. 'extraction. Parts are retained.', scalar @parts);
		log_info('Extracting GTDBtk database to %s', $dest);
		my @tar = ('tar', (main::have_pigz() ? '--use-compress-program=pigz' : '-z'),
			'-xf', '-', '-C', $dest);
		log_debug('Command: cat %s | %s', join(' ', @parts), main::shell_quote(@tar));
		local $SIG{PIPE} = 'IGNORE';
		open(my $tfh, '|-', @tar) or die "Cannot run tar: $!\n";
		binmode $tfh;
		for my $p (@parts) {
			open my $in, '<', $p or die "Cannot read $p: $!\n";
			binmode $in;
			my $buf;
			while (1) {
				my $n = sysread($in, $buf, 8 << 20);
				die "Read error on $p: $!\n" unless defined $n;
				last if $n == 0;
				print {$tfh} $buf or do {
					close $tfh;
					die "tar stopped reading while extracting $p: $!\n";
				};
			}
			close $in;
		}
		close($tfh) or die "Extraction of GTDB-Tk parts failed (tar exit status "
			. ($? >> 8) . ")\n";
	} else {
		my $src = $files->{tk_database} // die "GTDB-Tk full package was not downloaded\n";
		die "GTDB-Tk archive not found: $src\n" unless -f $src;
		log_info('Extracting GTDBtk database to %s', $dest);
		log_debug('Source: %s', $src);
		main::run_checked("Extraction of $src", 'tar', _tar_decompress_opt($src), '-xf', $src, '-C', $dest);
	}
	log_info('GTDBtk extraction complete');
	main::write_text($finished, main::py_now());
}

sub _taxonomy_sources {
	my ($self, $f) = @_;
	for my $k (qw(bac_taxonomy arc_taxonomy bac_metadata arc_metadata)) {
		die "Missing downloaded file for $k\n" unless defined $f->{$k};
		die "File not found: $f->{$k}\n" unless -e $f->{$k};
	}
	# "cat bac_tax arc_tax" and "cat bac_md <(zcat arc_md | tail -n +2)"
	return ([[$f->{bac_taxonomy}, 0], [$f->{arc_taxonomy}, 0]],
		[[$f->{bac_metadata}, 0], [$f->{arc_metadata}, 1]]);
}

sub _open_out {
	my ($file) = @_;
	open my $fh, '>', $file or die "Cannot write $file: $!\n";
	return $fh;
}

sub _close_out {
	my ($fh, $file) = @_;
	close $fh or die "Cannot write $file: $!\n";
}

sub lineage_file    { main::join_path($_[1], "gtdb_r" . $_[0]->v0 . "_lineageGTDB.tab") }
sub clustering_file { main::join_path($_[1], "gtdb_r" . $_[0]->v0 . "_clustering.tab") }

sub _write_lineage {
	my ($self, $dest, $tax_src, $strip) = @_;
	log_info('Making species lineage');
	my $out = $self->lineage_file($dest);
	my $o = _open_out($out);
	my %species;
	main::each_concat_row($tax_src, sub {
		my ($r) = @_;
		my ($acc, $lineage) = @$r;
		die "Taxonomy row without taxonomy column: " . join("\t", @$r) . "\n" unless defined $lineage;
		my @parts = main::py_split_semicolon($lineage);
		@parts = map { main::py_strip($_) } @parts if $strip;
		return if $species{ $parts[-1] }++;
		print {$o} main::row_quote_none($parts[-1], 0, @parts);
	});
	_close_out($o, $out);
}

sub _md_columns {
	my ($hdr, @need) = @_;
	my %idx;
	$idx{ $hdr->[$_] } = $_ for 0 .. $#$hdr;    # DictReader: last duplicate wins
	for my $n (@need) { die "Metadata column '$n' not found\n" unless exists $idx{$n} }
	return @idx{@need};
}

# Make lineage and clustering tables (r220 format)
sub _format_taxonomy {
	my ($self, $dest, $files) = @_;
	log_info('Producing taxonomy tables');
	log_debug('Destination: %s', $dest);
	log_info('Concatenate source archaeal and bacterial tables');
	my ($tax_src, $md_src) = $self->_taxonomy_sources($files);
	$self->_write_lineage($dest, $tax_src, 0);

	log_info('Making species to representative accession mapping');
	my $out = $self->clustering_file($dest);
	my $o = _open_out($out);
	my ($hdr, $i_acc, $i_rep, $i_tax);
	main::each_concat_row($md_src, sub {
		my ($r) = @_;
		if (!$hdr) {
			$hdr = $r;
			($i_acc, $i_rep, $i_tax) = _md_columns($hdr, qw(accession gtdb_representative gtdb_taxonomy));
			return;
		}
		my $rep = $r->[$i_rep];
		return unless defined $rep && $rep eq 't';
		my ($tax, $acc) = ($r->[$i_tax], $r->[$i_acc]);
		die "Metadata row for $acc lacks gtdb_taxonomy\n" unless defined $tax;
		print {$o} main::row_quote_none((main::py_split_semicolon($tax))[-1], $acc);
	});
	_close_out($o, $out);
	log_info('Taxonomy tables finished');
}

# Accessions of all extracted marker sequences (kept from the Python API).
sub _accessions_from_markers {
	my ($self, $marker_dir) = @_;
	my %accs;
	require File::Find;
	File::Find::find({ no_chdir => 1, wanted => sub {
		return unless -f $_ && /\.faa\z/;
		open my $fh, '<', $_ or die "Cannot read $_: $!\n";
		while (my $l = <$fh>) {
			next unless substr($l, 0, 1) eq '>';
			$accs{ main::py_strip(substr($l, 1)) } = 1;
		}
		close $fh;
	} }, $marker_dir);
	return \%accs;
}

sub required_gtdbtk {
	my ($version) = @_;
	for my $e (@GTDBTK_VERSIONS) {
		return [$e->[1], $e->[2]] if $e->[0] == $version;
	}
	return undef;
}

# GTDB-Tk version this release is configured for: the lowest compatible one
# if the maximum is "Current" (no guarantee this table is up to date),
# otherwise the highest compatible one.
sub gtdbtk_version_use {
	my ($self) = @_;
	my $vt = required_gtdbtk($self->{version});
	my $known = defined $vt;
	$vt //= ['2.5', '2.5'];
	my ($min, $max) = @$vt;
	return ($max eq 'Current' ? $min : $max, $min, $max, $known);
}

sub _processing_metadata {
	my ($self, $dest, %kw) = @_;
	log_info('Writing metadata to meta.json');
	my @meta = (
		version        => \main::py_float($self->{version}),
		urls           => [@{ $self->{required_urls} }],
		date           => main::py_now(),
		script_version => \main::py_float($main::SCRIPT_VERSION),
		# DEVIATION 14: names this script
		summary        => 'GTDB formatted for MATAFILER using get_gtdb.pl',
	);
	for my $k (qw(dir tk)) { push @meta, $k => $kw{$k} if exists $kw{$k} }
	main::write_text(main::join_path($dest, 'meta.json'), main::py_json(\@meta));
}

sub extract {
	my ($self, %a) = @_;
	my $dl_dir = $a{dl_dir} // main::getcwd();
	my $dest = $a{dest_dir} // main::join_path(main::getcwd(), 'output');
	my $tk = $a{tk} // 'skip';
	$tk = 'skip' if !$tk || $tk eq 'False';
	my $files = $self->download(download_dir => $dl_dir, offline => 1,
		download_tk => ($tk eq 'skip' ? 0 : $tk));
	make_path($dest) unless -d $dest;
	$self->_extract_markers(main::join_path($dest, 'markerGenes'), $files);
	$self->_format_taxonomy($dest, $files);
	# DEVIATION 9: "skip" means skip
	$self->_extract_tk(main::join_path($dest, 'gtdb'), $files, $tk) if $tk ne 'skip';
	$self->_processing_metadata($dest);
}

sub ask_config_update {
	my ($do_update) = @_;
	return $do_update if defined $do_update;
	if (!-t STDIN) {
		log_info('Run in non-interactive terminal, skipping config update');
		return 0;
	}
	my $mf4 = main::get_mf4_dir();
	if (!defined $mf4) {
		log_warn("MATAFILER installation not found. \n"
			. "MATAFILER configuration will not be updated, as neither MF4DIR (or "
			. "MGTKDIR) nor the location of this script point to a MATAFILER "
			. "directory. \n"
			. "Rerun get_gtdb.pl configure on the system with MATAFILER "
			. "installed to complete configuration. \n"
			. "Additionally, instructions on updating the configuration will "
			. "be written to README.txt in the output "
			. "directory and to the terminal at the end of the program.");
		return 0;
	}
	log_info('MATAFILER installation found in %s', $mf4);
	my $x = main::prompt_input_set(
		'Do you want this script to automatically update MATAFILER '
		. 'configuration to use the downloaded databases?',
		[qw(yes no cancel)], 1);
	exit 0 if $x eq 'cancel';
	return $x eq 'yes';
}

# Directory for GTDBtk_DB: gtdb/<release> when exactly one subdirectory exists.
sub _tk_data_dir {
	my ($gtdb) = @_;
	return undef unless -d $gtdb;
	opendir my $dh, $gtdb or die "Cannot read $gtdb: $!\n";
	my @sub = sort grep { !/^\./ && -d main::join_path($gtdb, $_) } readdir $dh;
	closedir $dh;
	return main::join_path($gtdb, $sub[0]) if @sub == 1;
	log_warn('%s contains %d subdirectories; GTDBtk_DB will point to %s itself.',
		$gtdb, scalar @sub, $gtdb) if @sub > 1;
	return $gtdb;
}

sub update_config {
	my ($self, %a) = @_;
	my $mf4 = main::get_mf4_dir();
	if (!defined $mf4) {
		log_error('MATAFILER installation not found. Please ensure environmental '
			. 'variable MF4DIR contains the path to the MATAFILER directory.');
		die "MATAFILER installation not found\n";
	}
	log_info('MATAFILER directory: %s', $mf4);
	my $out_res = main::resolve_path($a{out_dir});
	my ($db_exp, $db_raw, $db_file) = main::find_dbdir($mf4);
	my ($use_abs, $db_res) = (1, undef);
	if (defined $db_exp) {
		$db_res = main::resolve_path($db_exp);
		log_info('MATAFILER DBDir: %s (%s)', $db_raw, $db_file);
		log_debug('MATAFILER DBDir resolved: %s', $db_res);
		log_debug('GTDB Database resolved: %s', $out_res);
		my $prefix = $db_res eq '/' ? '/' : "$db_res/";
		if (index($out_res, $prefix) == 0 && $out_res ne $db_res) {
			$use_abs = 0;
			log_debug('Data in DBDir, using [DBDir] format.');
		} else {
			log_info('Downloaded database is not in DBDir, will use absolute paths.');
		}
	} else {
		log_info('DBDir not found in MATAFILER configuration, will use absolute '
			. 'paths to database directories.');
	}

	my ($ver_use_s, $vmin, $vmax, $known) = $self->gtdbtk_version_use;
	log_warn('Cannot determine suitable GTDB-Tk version. Paths will be updated '
		. 'assuming a GTDB-Tk version of 2.5') unless $known;
	log_debug('Compatible GTDB-TK - Min: %s, Max %s', $vmin, $vmax);
	log_debug('Assuming GTDB-TK version %s', $ver_use_s);
	# GTDB-Tk >= 2.5 does not use mash, so its config line is commented out
	my $no_mash = main::version_ge($ver_use_s, '2.5');
	log_debug('Do not use mash: %s', $no_mash ? 'True' : 'False');

	my %to_change = (
		GTDBPath  => main::join_path($out_res, 'markerGenes'),
		GTDB_GTDB => $self->lineage_file($out_res),
		GTDB_lnks => $self->clustering_file($out_res),
	);
	# DEVIATION 5
	my $tk_dir = _tk_data_dir(main::join_path($out_res, 'gtdb'));
	if (defined $tk_dir) {
		$to_change{GTDBtk_DB} = $tk_dir;
		$to_change{GTDBtk_mash} = main::join_path($tk_dir, 'mashD');
	} else {
		log_info('No GTDB-Tk data in %s/gtdb; GTDBtk_DB and GTDBtk_mash are left '
			. 'unchanged.', $out_res);
	}
	for my $k (sort keys %to_change) {
		if (!$use_abs) {
			(my $rel = $to_change{$k}) =~ s{^\Q$db_res\E/*}{};
			$to_change{$k} = "[DBDir]/$rel";      # DEVIATION 3
		}
		die "Path for $k contains characters MATAFILER config values cannot hold "
			. "(tab, '#', '^'): $to_change{$k}\n" if $to_change{$k} =~ /[\t#^\r\n]/;
	}

	my $db_config = main::join_path($mf4, 'Mods', 'config_DBs.txt');
	my $i = 1;
	my $bup;
	while (1) {
		$bup = main::join_path($mf4, 'Mods', "config_DBs.bup$i");
		last unless -e $bup;
		$i++;
	}
	log_info('Backing up current MATAFILER config to %s', $bup);
	File::Copy::copy($db_config, $bup) or die "Cannot back up $db_config to $bup: $!\n";
	chmod((stat $db_config)[2] & 07777, $bup);

	my @mod;
	open my $fh, '<', $db_config or die "Cannot read $db_config: $!\n";
	while (my $line = <$fh>) {
		if (length(main::py_strip($line)) < 1) { push @mod, $line; next; }
		# DEVIATION 4: idempotent mash handling, only on the key itself
		if (exists $to_change{GTDBtk_mash}) {
			if (!$no_mash && $line =~ /^#+GTDBtk_mash(?:\t|\s|$)/) {
				$line =~ s/^#+//;
				log_debug('Restored mash path in config');
			} elsif ($no_mash && $line =~ /^GTDBtk_mash(?:\t|\s|$)/) {
				log_debug('Removed mash path in config');
				$line = '#' . main::py_strip($line);
				$line .= " $UPDATED_TAG" unless $line =~ /#Updated by get_gtdb\.p[ly]/;
				$line .= "\n";
			}
		}
		if (substr(main::py_strip($line), 0, 1) eq '#') { push @mod, $line; next; }
		my @parts = split /\t/, $line, -1;
		(my $key = $parts[0]) =~ s/\r?\n\z//;
		if (exists $to_change{$key}) {
			push @mod, "$key\t$to_change{$key}\t$UPDATED_TAG\n";
			log_debug('Updated %s to %s', $key, $to_change{$key});
		} else {
			push @mod, $line;
		}
	}
	close $fh;
	log_info('Writing updated config');
	main::write_text($db_config, join('', @mod));

	# DEVIATION 3: the user config (config.txt) takes precedence in the pipeline
	for my $user (main::join_path($mf4, 'config.txt')) {
		next unless -e $user;
		my $cfg = main::parse_mgtk_config($user);
		my @over = grep { exists $cfg->{$_} } sort keys %to_change;
		log_warn('%s defines %s, which override Mods/config_DBs.txt. Update or '
			. 'remove these lines there to use the new database.', $user, join(', ', @over)) if @over;
	}

	if (!$a{install_tk}) {
		log_info("No changes made to GTDB-Tk installed in $GTDBTK_ENV "
			. 'environment. You can run this script with --install-tk to '
			. 'install a version of GTDB-Tk known to work with this release.');
		return;
	}
	my $mamba = (defined $ENV{MAMBA_EXE} && length $ENV{MAMBA_EXE}) ? $ENV{MAMBA_EXE} : 'micromamba';
	log_info('Installing gtdbtk ==%s into environment %s', $ver_use_s, $GTDBTK_ENV);
	my ($st, $out, $err) = main::run_capture($mamba, 'install', '--name', $GTDBTK_ENV,
		'-c', 'bioconda', '-c', 'conda-forge', '-y', '--channel-priority', 'flexible',
		"gtdbtk==$ver_use_s");
	if ($st != 0) {
		log_error('gtdbtk installation failed. Dumping output.');
		log_error('Stderr');
		log_error('%s', $err);
		log_error('Stdout');
		log_error('%s', $out);
		die "gtdbtk installation failed\n";
	}
}

sub instructions {
	my ($self, $dest, $tk) = @_;
	my $text = $self->_mgtk_instructions($dest) . "\n\n" . $self->_gtdbtk_instructions($dest, $tk);
	main::write_text(main::join_path($dest, 'README.txt'), $text);
	print "$text\n";
	print "\nThese instructions were also written to:\n";
	print main::join_path($dest, 'README.txt'), "\n";
	return $text;
}

sub _mgtk_instructions {
	my ($self, $dest) = @_;
	my $mf4 = main::instructions_mf4_dir();
	my $v = $self->v0;
	return "We recommend moving the output directory ($dest) to your MATAFILER "
		. "DBDir, and to a version specific subdirectory i.e. to\n\n"
		. "<DBDir>/MarkerG/GTDB_r${v}_MGTK\n\n"
		. "To use this version in MATAFILER, you must update config_DBs.txt. This "
		. "is in $mf4/Mods/config_DBs.txt.\n"
		. "Update the lines:\n\n"
		. "GTDBPath\t[DBDir]/MarkerG/GTDB_r${v}_MGTK/markerGenes/\n"
		. "GTDB_GTDB\t[DBDir]/MarkerG/GTDB_r${v}_MGTK/gtdb_r${v}_lineageGTDB.tab\n"
		. "GTDB_lnks\t[DBDir]/MarkerG/GTDB_r${v}_MGTK/gtdb_r${v}_clustering.tab\n\n"
		. "Alternatively, you can run \n"
		. "get_gtdb.pl configure -d $dest\n"
		. "on the system with MATAFILER installed.";
}

sub _gtdbtk_instructions {
	my ($self, $dest, $tk) = @_;
	my $mf4 = main::instructions_mf4_dir();
	my $v = $self->v0;
	my ($use, $vmin, $vmax, $known) = $self->gtdbtk_version_use;
	my $no_mash = main::version_ge($use, '2.5');
	my $tk_inst = $known
		? "GTDBtk version between $vmin and $vmax is required for "
			. "this release of GTDB. Please install this using micromamba into "
			. "environment '$GTDBTK_ENV' i.e.\n\n"
			. "micromamba install --name $GTDBTK_ENV gtdbtk==$use\n"
		: "We could not identify which version of GTDBtk is required "
			. "for this release. Please check "
			. "https://ecogenomics.github.io/GTDBTk/installing/"
			. "index.html#gtdb-tk-reference-data and install the correct "
			. "required version into the $GTDBTK_ENV environment i.e.\n\n"
			. "micromamba install --name $GTDBTK_ENV gtdbtk==ver";
	my $tk_db;
	if ($tk && $tk ne 'skip') {
		my $sub = _tk_data_dir(main::join_path($dest, 'gtdb'));
		my $suffix = (defined $sub && $sub ne main::join_path($dest, 'gtdb'))
			? '/' . main::basename($sub) : '';
		$tk_db = "GTDBtk requires release specific databases. These have been "
			. "download and extracted by this script. \n"
			. "To use this version, you must update config_DBs.txt. This is "
			. "in $mf4/Mods/config_DBs.txt.\n"
			. "Update the line" . ($no_mash ? '' : 's') . ":\n\n"
			. "GTDBtk_DB\t[DBDir]/MarkerG/GTDB_r${v}_MGTK/gtdb$suffix"
			. ($no_mash ? '' : "\nGTDBtk_mash\t[DBDir]/MarkerG/GTDB_r${v}_MGTK/gtdb$suffix/mashD");
	} else {
		$tk_db = "GTDBtk requires release specific databases. These were not "
			. "downloaded by this script. \n"
			. "Please see https://ecogenomics.github.io/GTDBTk/installing/"
			. "index.html#gtdb-tk-reference-data for how to download this "
			. "data. \n"
			. "Additionally, GTDBtk has a script which will download the "
			. "most current version, 'download-db.sh'. \n\n"
			. "Once you have downloaded the correct version, you must update "
			. "config_DBs.txt. This is "
			. "in $mf4/Mods/config_DBs.txt.\n"
			. "Update the line" . ($no_mash ? '' : 's') . ":\n\n"
			. "GTDBtk_DB\t[DBDir]/MarkerG/GTDB_r${v}_MGTK/gtdbtk"
			. ($no_mash ? '' : "\nGTDBtk_mash\t[DBDir]/MarkerG/GTDB_r${v}_MGTK/gtdbtk/mashD");
	}
	return "$tk_db\n\n$tk_inst";
}

# Factory: a version specific subclass where one exists.
sub from_version {
	my ($version, %kw) = @_;
	return GTDB226->new($version, %kw) if $version == 226;
	if ($version != 220) {
		log_warn("Using default download URLs and extraction methods\n"
			. "This method was tested to work with r220, but are "
			. "untested for " . main::py_float($version) . ".");
	}
	return GTDBVersion->new($version, %kw);
}

###########################################################################
package GTDB226;
###########################################################################
# v226: marker genes can come from non-representative genomes, so the
# clustering table lists every accession and GTDBmg.tax keeps duplicates.

our @ISA = ('GTDBVersion');
BEGIN { *log_info = \&main::log_info; *log_debug = \&main::log_debug; }

sub _make_gtdbmg_tax {
	my ($self, $dest, $tax, $lineage) = @_;
	my $out = main::join_path($dest, 'GTDBmg.tax');
	log_info('Creating GTDBmg.tax with duplicates for species');
	log_debug('Output to %s', $out);
	my %lineage_map;
	open my $lf, '<', $lineage or die "Cannot read $lineage: $!\n";
	while (my $l = <$lf>) {
		$l =~ s/\r?\n\z//;
		my @r = split /\t/, $l, -1;
		$lineage_map{ $r[0] } = join(';', @r[2 .. $#r]);
	}
	close $lf;
	open my $tf, '<', $tax or die "Cannot read $tax: $!\n";
	my $o = GTDBVersion::_open_out($out);
	log_debug('Write to %s', $out);
	while (my $l = <$tf>) {
		$l =~ s/\r?\n\z//;
		my @r = split /\t/, $l, -1;
		die "Malformed clustering row in $tax: '$l'\n" if @r < 2;
		my ($species, $acc) = @r;
		die "Species '$species' from $tax is missing in $lineage\n" unless exists $lineage_map{$species};
		print {$o} main::row_quote_minimal(main::py_strip($acc), main::py_strip($lineage_map{$species}));
	}
	close $tf;
	GTDBVersion::_close_out($o, $out);
}

sub _format_taxonomy {
	my ($self, $dest, $files) = @_;
	log_info('Producing taxonomy tables');
	log_debug('Destination: %s', $dest);
	log_info('Concatenate source archaeal and bacterial tables');
	my ($tax_src, $md_src) = $self->_taxonomy_sources($files);
	$self->_write_lineage($dest, $tax_src, 1);

	log_info('Making species to accession mapping');
	log_info('For v226, we output species for all accessions as some genes '
		. 'are not from representatives');
	my $out = $self->clustering_file($dest);
	my $o = GTDBVersion::_open_out($out);
	my ($hdr, $i_acc, $i_tax);
	main::each_concat_row($md_src, sub {
		my ($r) = @_;
		if (!$hdr) {
			$hdr = $r;
			($i_acc, $i_tax) = GTDBVersion::_md_columns($hdr, qw(accession gtdb_taxonomy));
			return;
		}
		my ($acc, $tax) = ($r->[$i_acc], $r->[$i_tax]);
		die "Metadata row for $acc lacks gtdb_taxonomy\n" unless defined $tax;
		print {$o} main::row_quote_none(main::py_strip((main::py_split_semicolon($tax))[-1]),
			main::py_strip($acc));
	});
	GTDBVersion::_close_out($o, $out);

	$self->_make_gtdbmg_tax(main::join_path($dest, 'markerGenes'),
		$self->clustering_file($dest), $self->lineage_file($dest));
	log_info('Taxonomy tables finished');
}

###########################################################################
package main;
###########################################################################

sub append_file {
	my ($from, $to) = @_;
	open my $in, '<', $from or die "Cannot read $from: $!\n";
	open my $out, '>>', $to or die "Cannot append to $to: $!\n";
	binmode $in; binmode $out;
	my $buf;
	while (1) {
		my $n = sysread($in, $buf, 1 << 20);
		die "Read error on $from: $!\n" unless defined $n;
		last if $n == 0;
		print {$out} $buf or die "Cannot append to $to: $!\n";
	}
	close $in;
	close $out or die "Cannot append to $to: $!\n";
}

sub py_str_tk {
	my ($tk) = @_;
	return 'False' if !defined $tk || $tk eq '' || $tk eq '0';
	return $tk;
}

sub meta_from_dir {
	my ($dir) = @_;
	my $json = join_path($dir, 'meta.json');
	if (!-e $json) {
		log_error('No meta.json found in %s, unable to determine version. Please '
			. 'provide a directory which was created by this tool.', $dir);
		die "No meta.json in $dir\n";
	}
	open my $fh, '<', $json or die "Cannot read $json: $!\n";
	local $/;
	my $txt = <$fh>;
	close $fh;
	my $meta = eval { JSON::PP->new->decode($txt) };
	die "Cannot parse $json: $@" unless $meta;
	die "$json has no version\n" unless defined $meta->{version};
	return $meta;
}

###########################################################################
# Command line interface
###########################################################################

my %HELP = (
	debug => "Produce debug messages. This will output any shell commands being run.",
	dest => "Directory for MATAFILER formatted database. When extracting, contents will be overwritten. (default: ./output)",
	version => "Version of GTDB to download. Can include minor versions (202.1)",
	tmp => "Temporary directory to download files from GTDB to. (default: current directory)",
	tk => "Method used to download the GTDBtk database. This is very large (~110GB for r220). By default this will be downloaded in parts which are then combined during extraction ('split'). You can instead download the full file ('full'). This can also be be skipped using 'skip' if you already have this data or will acquire it another way. (default: split)",
	test => "Do not download any files, instead create dummy files in the destination. Included to allow test of script logic without needing to download and extract large files.",
	install_tk => "Install the version of GTDBtk required for this GTDB release into the $GTDBTK_ENV environment. This is not done by default, as officially MATAFILER only supports one release of GTDB and installs the correct GTDBtk using installer shell scripts.",
);

my %SUB = (
	download => {
		help => 'Download GTDB and GTDBtk databases.',
		desc => 'Download GTDB and GTDBtk databases.',
		opts => [qw(tmp version tk test)],
	},
	extract => {
		help => 'Extract downloaded GTDB and GTDBtk databases to the format and structure required by MATAFILER.',
		desc => 'Extract downloaded GTDB and GTDBtk databases to the format and structure required by MATAFILER.',
		opts => [qw(tmp dest)],
	},
	configure => {
		help => 'Update MATAFILER configuration and install correct version of GTDBtk to use a downloaded and extracted database version.',
		desc => 'Update MATAFILER configuration and install correct version of GTDBtk to use a downloaded and extracted database version. Will require internet access for MATAFILER download. Should be run by the user with MATAFILER installed.',
		opts => [qw(dest install_tk)],
	},
	all => {
		help => 'Download and extract GTDB and GTDBtk, and configure MATAFILER to use download version. Equivalent to runining download, extract, then configure subcommands. Will required internet access for download, and GTDBtk install.',
		desc => 'Download and extract GTDB and GTDBtk, and configure MATAFILER to use download version. Equivalent to runining download, extract, then configure subcommands. Will required internet access for download, and GTDBtk install.',
		opts => [qw(tmp dest version tk install_tk test)],
	},
);
my @SUB_ORDER = qw(download extract configure all);

my %OPT = (
	tmp        => ['tmp|t=s',     '-t TMP, --tmp TMP',              '[-t TMP]'],
	dest       => ['dest|d=s',    '-d DEST, --dest DEST',           '[-d DEST]'],
	version    => ['version|v=s', '-v VERSION, --version VERSION',  '-v VERSION'],
	tk         => ['tk=s',        '--tk {split,full,skip}',         '[--tk {split,full,skip}]'],
	test       => ['test',        '--test',                         '[--test]'],
	install_tk => ['install-tk',  '--install-tk',                   '[--install-tk]'],
);

sub wrap {
	my ($text, $indent, $width) = @_;
	$width //= 78;
	my @lines;
	my $cur = '';
	for my $w (split ' ', $text) {
		if (length($cur) && length($indent) + length($cur) + 1 + length($w) > $width) {
			push @lines, $indent . $cur;
			$cur = $w;
		} else {
			$cur = length($cur) ? "$cur $w" : $w;
		}
	}
	push @lines, $indent . $cur if length $cur;
	return join("\n", @lines) . "\n";
}

sub help_entry {
	my ($name, $text) = @_;
	my $pad = ' ' x 24;
	my $head = "  $name";
	my $body = wrap($text, $pad);
	return length($head) <= 22 ? sprintf('%-24s', $head) . substr($body, 24) : "$head\n$body";
}

sub main_usage { "usage: $PROG [-h] [--debug] {" . join(',', @SUB_ORDER) . "} ...\n" }

sub main_help {
	my $h = main_usage() . "\n";
	$h .= wrap('Download and format versions of the GTDB database for use in the '
		. 'pipeline MATAFILER. This program expects GNU tar and either wget or curl to be available.', '');
	$h .= "\npositional arguments:\n  {" . join(',', @SUB_ORDER) . "}\n" . (' ' x 24) . "Subcommand help\n";
	$h .= help_entry("  $_", $SUB{$_}{help}) for @SUB_ORDER;
	$h .= "\noptions:\n" . help_entry('-h, --help', 'show this help message and exit')
		. help_entry('--debug', $HELP{debug});
	$h .= "\n" . wrap('This is separated into three steps, which can be run separately as '
		. "subcommands 'download', 'extract', and 'configure'. They can also be all run at once using subcommand 'all'.", '');
	$h .= "\n" . wrap('This program expects GNU tar and either wget or curl to be available. '
		. 'Downloads try to resume if they failed, and should be skipped for existing files '
		. '(of the same size) when the script is rerun.', '');
	$h .= "\n" . wrap('If the system you want to put the created database on lacks internet '
		. "access, use 'download' subcommand to download files, then transfer to target "
		. "system, then use 'extract' to process files.", '');
	$h .= "\n" . wrap("'configure' will update MATAFILER configuration to use an extracted "
		. 'database, and install an appropriate version of gtdbtk. If your group is sharing '
		. 'copies of the database, once it has been extracted, others should be able to '
		. "configure MATAFILER to use it using the 'configure' subcommand.", '');
	$h .= "\n" . wrap('The MATAFILER directory is taken from MF4DIR (or MGTKDIR), otherwise '
		. 'from the location of this script.', '');
	return $h;
}

sub sub_usage {
	my ($cmd) = @_;
	return "usage: $PROG $cmd [-h] " . join(' ', map { $OPT{$_}[2] } @{ $SUB{$cmd}{opts} }) . "\n";
}

sub sub_help {
	my ($cmd) = @_;
	my $h = sub_usage($cmd) . "\n" . wrap($SUB{$cmd}{desc}, '') . "\noptions:\n";
	$h .= help_entry('-h, --help', 'show this help message and exit');
	$h .= help_entry($OPT{$_}[1], $HELP{$_}) for @{ $SUB{$cmd}{opts} };
	$h .= help_entry('--debug', $HELP{debug});
	return $h;
}

sub usage_error {
	my ($usage, $msg) = @_;
	print STDERR $usage, "$PROG: error: $msg\n";
	exit 2;
}

sub parse_cli {
	my (@argv) = @_;
	my $debug = 0;
	while (@argv && $argv[0] =~ /^-/) {
		my $a = shift @argv;
		if ($a eq '--debug') { $debug = 1 }
		elsif ($a eq '-h' || $a eq '--help') { print main_help(); exit 0 }
		else { usage_error(main_usage(), "unrecognized arguments: $a") }
	}
	if (!@argv) { print main_help(); exit 0 }
	my $cmd = shift @argv;
	usage_error(main_usage(), "argument subcommand: invalid choice: '$cmd' (choose from "
		. join(', ', map {"'$_'"} @SUB_ORDER) . ')') unless $SUB{$cmd};
	my %o = (debug => $debug);
	my @spec = ('debug' => \$o{debug}, 'help|h' => \$o{help});
	for my $k (@{ $SUB{$cmd}{opts} }) { push @spec, $OPT{$k}[0] => \$o{$k} }
	Getopt::Long::Configure(qw(gnu_getopt no_ignore_case));
	my @warn;
	my $ok = do {
		local $SIG{__WARN__} = sub { push @warn, $_[0] };
		GetOptionsFromArray(\@argv, @spec);
	};
	usage_error(sub_usage($cmd), join('', @warn) =~ s/\n+\z//r) unless $ok;
	usage_error(sub_usage($cmd), "unrecognized arguments: @argv") if @argv;
	if ($o{help}) { print sub_help($cmd); exit 0 }
	my %has = map { $_ => 1 } @{ $SUB{$cmd}{opts} };
	if ($has{version}) {
		usage_error(sub_usage($cmd), 'the following arguments are required: -v/--version')
			unless defined $o{version};
		usage_error(sub_usage($cmd), "argument -v/--version: invalid float value: '$o{version}'")
			unless $o{version} =~ /^\s*[+]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?\s*$/;
		$o{version} += 0;
	}
	if ($has{tk}) {
		$o{tk} //= 'split';
		usage_error(sub_usage($cmd), "argument --tk: invalid choice: '$o{tk}' (choose from 'split', 'full', 'skip')")
			unless $o{tk} =~ /^(?:split|full|skip)$/;
	}
	$o{tmp} = norm_path($has{tmp} && defined $o{tmp} ? $o{tmp} : getcwd());
	$o{dest} = norm_path($has{dest} && defined $o{dest} ? $o{dest} : join_path(getcwd(), 'output'));
	return ($cmd, \%o);
}

sub cli_download {
	my ($a) = @_;
	my $d = GTDBVersion::from_version($a->{version}, test => $a->{test});
	$d->download(download_dir => $a->{tmp}, download_tk => $a->{tk});
}

sub cli_extract {
	my ($a) = @_;
	my $meta = meta_from_dir($a->{tmp});
	my $d = GTDBVersion::from_version($meta->{version});
	$d->extract(dl_dir => $a->{tmp}, dest_dir => $a->{dest}, tk => $meta->{tk});
}

sub cli_configure {
	my ($a) = @_;
	my $meta = meta_from_dir($a->{dest});
	my $d = GTDBVersion::from_version($meta->{version});
	$d->update_config(out_dir => $a->{dest}, install_tk => $a->{install_tk});
}

sub cli_all {
	my ($a) = @_;
	my $d = GTDBVersion::from_version($a->{version}, test => $a->{test});
	my $config = GTDBVersion::ask_config_update();
	$d->download(download_dir => $a->{tmp}, download_tk => $a->{tk});
	$d->extract(dl_dir => $a->{tmp}, dest_dir => $a->{dest}, tk => $a->{tk});
	if ($config) {
		$d->update_config(out_dir => $a->{dest}, install_tk => $a->{install_tk});
	} else {
		log_warn('Did not update MATAFILER configuration.');
		log_warn("To update MATAFILER to use new database version, either run "
			. "'get_gtdb.pl configure -d %s' on the system with MATAFILER installed "
			. 'or follow the instructions in %s.', $a->{dest}, join_path($a->{dest}, 'README.txt'));
		$d->instructions($a->{dest}, $a->{tk});
	}
}

sub main {
	my ($cmd, $a) = parse_cli(@ARGV);
	$DEBUG = $a->{debug} ? 1 : 0;
	greet();
	my %dispatch = (download => \&cli_download, extract => \&cli_extract,
		configure => \&cli_configure, all => \&cli_all);
	my $ok = eval { $dispatch{$cmd}->($a); 1 };
	if (!$ok) {
		my $err = $@ // 'unknown error';
		$err =~ s/\n+\z//;
		log_error('%s', $err);
		exit 1;
	}
	exit 0;
}

main() unless caller;
