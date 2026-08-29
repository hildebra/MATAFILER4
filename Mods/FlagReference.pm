package Mods::FlagReference;
#Renders command-line help for the MATAFILER4 entry points straight from
#docs/flag_reference.md, so the option tables exist exactly once. Every script
#that accepts flags calls printFlagHelp() instead of carrying its own copy of
#the option list, which is what let the copies drift apart in the past.
use warnings;
use strict;
use Exporter qw(import);
use File::Basename qw(dirname);
use Cwd qw(abs_path);
use Text::Wrap qw(wrap);

our $VERSION = 0.11;
our @EXPORT_OK = qw(printFlagHelp flagReferencePath readFlagReference
	helpRequested scriptOptionSpecs scriptOptionNames);

#docs/ sits next to Mods/, so the reference is found from the module itself and
#not from the calling script, which lives at three different depths.
sub flagReferencePath {
	my $modulePath = abs_path(__FILE__) || __FILE__;
	return dirname(dirname($modulePath)) . "/docs/flag_reference.md";
}

#One definition of what counts as a help request, for every entry point:
#-help/-h/-? with one or two leading dashes.
sub helpRequested {
	my (@argv) = @_;
	@argv = @ARGV unless (@argv);
	return scalar(grep { defined($_) && /^--?(?:help|h|\?)$/ } @argv);
}

#Parses a script's GetOptions(...) call and returns one record per option:
#{ spec, name, aliases => [...], type, comment }. Used by the documentation
#tests, and by anything else that needs to know what a script actually accepts
#without duplicating the scanning logic.
sub scriptOptionSpecs {
	my ($path) = @_;
	open(my $fh, '<', $path) or die "Cannot read $path: $!";
	my @lines = <$fh>;
	close($fh);

	my $start = -1;
	for my $i (0 .. $#lines) {
		if ($lines[$i] =~ /^\s*GetOptions\s*\(/) { $start = $i; last; }
	}
	return () if ($start < 0);

	my $quote = chr(34);
	my @specs;
	for my $i ($start .. $#lines) {
		my $line = $lines[$i];
		if ($line =~ /^\s*$quote([^$quote]+)$quote\s*=>/) {
			my $spec = $1;
			my ($names, $argument) = $spec =~ /^([^=!:]+)([=!:].*)?$/;
			my ($comment) = $line =~ /#\s*(.*?)\s*$/;
			my @aliases = split(/\|/, $names);
			push @specs, {
				spec    => $spec,
				name    => $aliases[0],
				aliases => [@aliases],
				type    => _argumentType($argument),
				comment => (defined($comment) ? $comment : ""),
			};
		}
		#the call closes on a line holding only ")" / ") or die ..." / ");"
		last if ($i > $start && $line =~ /^\s*\)\s*(?:or\b|;|$)/);
	}
	return @specs;
}

sub _argumentType {
	my ($argument) = @_;
	return "flag" unless (defined($argument) && length($argument));
	return "flag"             if ($argument eq '!');
	return "optional integer" if ($argument =~ /^:/);
	return "integer"          if ($argument =~ /=i/);
	return "float"            if ($argument =~ /=f/);
	return "string"           if ($argument =~ /=s/);
	return $argument;
}

#Flat set of every accepted name, aliases included.
sub scriptOptionNames {
	my ($path) = @_;
	my %names;
	for my $spec (scriptOptionSpecs($path)) {
		$names{$_} = 1 for @{$spec->{aliases}};
	}
	return %names;
}

sub _plainCell {
	my ($text) = @_;
	$text = "" unless (defined($text));
	$text =~ s/^\s+|\s+$//g;
	$text =~ s/`//g;
	$text =~ s/\*\*//g;
	$text =~ s/\[([^\]]+)\]\([^\)]+\)/$1/g;
	return $text;
}

#A "## <something>.pl" heading starts another script's section; "## Legacy flag
#names" ends the MATAF4.pl tables. Everything else at ## or ### is a subsection.
sub _isSectionBoundary {
	my ($line) = @_;
	return 1 if ($line =~ /^##\s+\S+\.pl\s*$/);
	return 1 if ($line =~ /^##\s+Legacy flag names\s*$/i);
	return 0;
}

#Returns the parsed sections for one script heading, or an empty list when the
#reference cannot be read.
sub readFlagReference {
	my ($script, $reference) = @_;
	$reference = flagReferencePath() unless (defined($reference) && length($reference));
	open(my $fh, '<', $reference) or return ();

	my @sections;
	my $current;
	my $inScript = 0;
	while (my $line = <$fh>) {
		if ($line =~ /^##\s+\Q$script\E\s*$/) { $inScript = 1; next; }
		next unless ($inScript);
		last if (_isSectionBoundary($line));

		if ($line =~ /^#{2,3}\s+(.+?)\s*$/) {
			$current = { title => $1, options => [], legacy => [], bullets => [] };
			push @sections, $current;
			next;
		}
		#markdown list items inside a section (e.g. the geneCat run-mode list)
		#are rendered too, so such enumerations also live only in the reference
		if (defined($current) && $line =~ /^[-*]\s+(\S.*?)\s*$/) {
			push @{$current->{bullets}}, _plainCell($1);
			next;
		}
		next unless ($line =~ /^\|\s*`-/);
		next unless (defined($current));

		my @cells = split(/\|/, $line, -1);
		if (@cells >= 7) {
			my ($aliases, $type, $default, $status, $description)
				= map { _plainCell($_) } @cells[1 .. 5];
			push @{$current->{options}}, {
				aliases => $aliases, type => $type, default => $default,
				status => $status, description => $description,
			};
		} elsif (@cells == 5) {
			#the "internal and deprecated" tables carry three columns
			my ($aliases, $status, $purpose) = map { _plainCell($_) } @cells[1 .. 3];
			push @{$current->{legacy}}, {
				aliases => $aliases, status => $status, purpose => $purpose,
			};
		}
	}
	close($fh);
	return @sections;
}

sub _printOption {
	my ($option) = @_;
	my $signature = $option->{aliases};
	$signature .= " <$option->{type}>"
		if ($option->{type} ne "" && lc($option->{type}) ne "flag");

	my $details = $option->{description};
	my @metadata;
	push @metadata, "default: $option->{default}" if ($option->{default} ne "");
	push @metadata, "status: $option->{status}"
		if ($option->{status} ne "" && lc($option->{status}) ne "stable");
	$details .= " " if ($details ne "" && @metadata);
	$details .= "(" . join("; ", @metadata) . ")" if (@metadata);

	local $Text::Wrap::columns = 100;
	local $Text::Wrap::huge = 'overflow';
	print "  $signature\n";
	print wrap("      ", "      ", $details) . "\n" if ($details ne "");
}

sub _printBullet {
	my ($text) = @_;
	local $Text::Wrap::columns = 100;
	local $Text::Wrap::huge = 'overflow';
	print wrap("  ", "      ", $text) . "\n";
}

sub _printLegacy {
	my ($entry) = @_;
	local $Text::Wrap::columns = 100;
	local $Text::Wrap::huge = 'overflow';
	print "  $entry->{aliases}\n";
	my $details = $entry->{purpose};
	$details .= " ($entry->{status})" if ($entry->{status} ne "");
	print wrap("      ", "      ", $details) . "\n" if ($details ne "");
}

#printFlagHelp(script => "geneCat.pl", version => "0.58",
#              usage => ["geneCat.pl -GCd DIR [options]"],
#              summary => "...", notes => ["..."], exit => 1)
sub printFlagHelp {
	my (%args) = @_;
	my $script = $args{script} or die "printFlagHelp requires a script name\n";
	my $reference = $args{reference} || flagReferencePath();
	my @sections = readFlagReference($script, $reference);

	print "\n$script command-line help";
	print " (version $args{version})" if (defined($args{version}) && length($args{version}));
	print "\n";
	if (ref($args{usage}) eq 'ARRAY') {
		my $label = "Usage: ";
		for my $line (@{$args{usage}}) {
			print "$label$line\n";
			$label = "       ";
		}
	}
	if (defined($args{summary}) && length($args{summary})) {
		local $Text::Wrap::columns = 100;
		local $Text::Wrap::huge = 'overflow';
		print "\n" . wrap("", "", $args{summary}) . "\n";
	}
	print "\nOptions may be written with one or two leading dashes.\n";

	if (@sections && grep { @{$_->{options}} || @{$_->{legacy}} } @sections) {
		print "The complete option list below is read from docs/flag_reference.md.\n";
		foreach my $section (@sections) {
			next unless (@{$section->{options}} || @{$section->{legacy}}
				|| @{$section->{bullets}});
			print "\n$section->{title}:\n";
			_printBullet($_) foreach (@{$section->{bullets}});
			_printOption($_) foreach (@{$section->{options}});
			_printLegacy($_) foreach (@{$section->{legacy}});
		}
	} else {
		print "\nUnable to read the $script option tables from:\n  $reference\n";
		print "See that file for the complete flag reference.\n";
	}

	if (ref($args{notes}) eq 'ARRAY' && @{$args{notes}}) {
		local $Text::Wrap::columns = 100;
		local $Text::Wrap::huge = 'overflow';
		print "\n";
		print wrap("", "", $_) . "\n\n" foreach (@{$args{notes}});
	}

	print "\nFull reference: $reference\n\n";
	exit(0) if ($args{exit});
	return 1;
}

1;
