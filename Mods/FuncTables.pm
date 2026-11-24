package Mods::FuncTables;

# Per-MAG function tables for every MGS, made with yoke (github.com/falklab/yoke).
#
# This module holds the scheme registry shared by MGS.pl (which decides what to
# submit) and secScripts/MGS/funcTablesMGS.pl (which runs yoke), and prepares
# yoke's inputs from MATAFILER files:
#   - normalise_magvsgc: yoke reads MAGvsGC as a stream and starts a new MGS
#     whenever the MGS column changes, so the rows of one MGS must be adjacent.
#     It also treats "?" (unassigned) as an MGS and needs the column other_Genes.
#   - normalise_gene_annotations: MATAFILER *.geneAss files hold all functions of
#     a gene in one field (rtk syntax: ";" levels, "," and "|" between
#     functions). yoke expects one function ID per row.
#   - write_function_descriptions: yoke needs a headerless two-column table.

use strict;
use warnings;
use Exporter qw(import);
use File::Basename qw(dirname);
use File::Path qw(make_path);
use Mods::GenoMetaAss qw(gzipopen gzipwrite);
use Mods::WorkflowResilience qw(retry_rename retry_unlink);

our @EXPORT_OK = qw(
	function_table_scheme_names default_function_table_schemes
	parse_function_table_list function_table_scheme
	function_table_stage_stones resolve_function_table_scheme plan_function_tables
	normalise_magvsgc normalise_gene_annotations write_function_descriptions
	natural_mgs_sort
);

# stage: the geneCat annotation stage whose checkpoint marks the inputs complete
#   (func = diamond FuncAssign, emap = eggNOG-mapper).
# annotation: candidate gene->function files, relative to the gene catalogue.
# descriptions: candidate function description files; without one the
#   function ID is used as its description.
# modules: 1 when the functions are KEGG orthologs, which yoke module_score needs.
my %SCHEMES = (
	KEGG => {
		stage => 'func', modules => 1,
		annotation => ['Anno/Func/DIAass_KGM.srt.gzgeneAss.gz'],
		descriptions => ['Anno/Func/KEGG.descr'],
	},
	CZy => {
		stage => 'func', modules => 0,
		annotation => ['Anno/Func/DIAass_CZy.srt.gzgeneAss.gz'],
		descriptions => [],
	},
	TCDB => {
		stage => 'func', modules => 0,
		annotation => ['Anno/Func/DIAass_TCDB.srt.gzgeneAss.gz'],
		descriptions => ['Anno/Func/TCDB.anno'],
	},
	ABRc => {
		stage => 'func', modules => 0,
		annotation => ['Anno/Func/DIAass_ABRc.srt.gzgeneAss.gz'],
		descriptions => [],
	},
	VFA => {
		stage => 'func', modules => 0,
		annotation => ['Anno/Func/DIAass_VFA.srt.gzgeneAss.gz'],
		descriptions => [],
	},
	VFB => {
		stage => 'func', modules => 0,
		annotation => ['Anno/Func/DIAass_VFB.srt.gzgeneAss.gz'],
		descriptions => [],
	},
	emKO => {
		stage => 'emap', modules => 1,
		annotation => ['Anno/Func/emapper/eggNOGmapper_KO.geneAss.gz',
			'Anno/Func/emapper/eggNOGmapper_KO.geneAss'],
		descriptions => ['Anno/Func/KEGG.descr'],
	},
	emNOG => {
		stage => 'emap', modules => 0,
		annotation => ['Anno/Func/emapper/eggNOGmapper_NOG.geneAss.gz',
			'Anno/Func/emapper/eggNOGmapper_NOG.geneAss'],
		descriptions => ['Anno/Func/NOG.descr'],
	},
);
my @SCHEME_ORDER = qw(KEGG CZy TCDB ABRc VFA VFB emKO emNOG);

# Rows yoke would drop (its default filter keeps Completeness > 0 and
# Contamination < 100). Dropping them here also keeps yoke from building an
# MGS with no MAGs, which it cannot handle.
my @REQUIRED_MAGVSGC_COLUMNS = qw(MAG MGS Representative4MGS Completeness Contamination);
my @TAXONOMY_COLUMNS = qw(Domain Phylum Class Order Family Genus Species);
my %UNASSIGNED_MGS = map { $_ => 1 } ('', '?', '-');

sub function_table_scheme_names { return @SCHEME_ORDER; }

sub default_function_table_schemes { return join(',', @SCHEME_ORDER); }

sub function_table_scheme {
	my ($name) = @_;
	return defined($name) && exists($SCHEMES{$name}) ? { %{$SCHEMES{$name}}, name => $name } : undef;
}

# "0", "none" or "" select nothing. Unknown names are an error, so a typo does
# not silently skip a scheme. Order is kept, duplicates removed.
sub parse_function_table_list {
	my ($spec, $label, @valid) = @_;
	$label //= 'list';
	$spec //= '';
	$spec =~ s/^\s+|\s+$//g;
	return () if $spec eq '' || $spec eq '0' || lc($spec) eq 'none';
	my %known = map { $_ => 1 } @valid;
	my (@selected, %seen);
	for my $item (split /\s*,\s*/, $spec) {
		next unless length $item;
		die "Unknown $label entry '$item'; choose from: " . join(', ', @valid) . "\n"
			if @valid && !$known{$item};
		die "Invalid $label entry '$item'\n" unless $item =~ /\A[A-Za-z0-9_.-]+\z/;
		push @selected, $item unless $seen{$item}++;
	}
	return @selected;
}

sub function_table_stage_stones {
	my ($GCd) = @_;
	return {
		func => "$GCd/checkpoints/10.func.stone",
		emap => "$GCd/checkpoints/10.emap.stone",
	};
}

sub resolve_function_table_scheme {
	my ($GCd, $name) = @_;
	my $scheme = function_table_scheme($name)
		or die "Unknown function table scheme '$name'\n";
	my ($annotation) = grep { -s $_ } map { "$GCd/$_" } @{$scheme->{annotation}};
	my ($descriptions) = grep { -s $_ } map { "$GCd/$_" } @{$scheme->{descriptions}};
	return {
		name => $name,
		stage => $scheme->{stage},
		modules => $scheme->{modules},
		annotation => $annotation,
		descriptions => $descriptions,
	};
}

# Splits the requested schemes into those that can run now and those skipped
# (with the reason): a scheme runs when the geneCat stage that annotates it has
# finished (its checkpoint is non-empty) and its annotation file exists.
sub plan_function_tables {
	my ($GCd, $schemes, $stageStones) = @_;
	my (@ready, @skipped);
	for my $name (@{$schemes}) {
		my $scheme = resolve_function_table_scheme($GCd, $name);
		my $stone = $stageStones->{$scheme->{stage}};
		if (!defined($stone) || !-s $stone) {
			push @skipped, [$name, "annotation stage '$scheme->{stage}' has not finished"
				. (defined $stone ? " (no $stone)" : "")];
		} elsif (!defined $scheme->{annotation}) {
			push @skipped, [$name, "no annotation file (" . join(', ',
				@{function_table_scheme($name)->{annotation}}) . " under $GCd)"];
		} else {
			push @ready, $name;
		}
	}
	return (\@ready, \@skipped);
}

# MGS IDs in numeric order (MGS.2 before MGS.10); other IDs sort as text.
sub natural_mgs_sort {
	return sort {
		my ($ap, $an) = $a =~ /\A(.*?)(\d+)\z/ ? ($1, $2) : ($a, undef);
		my ($bp, $bn) = $b =~ /\A(.*?)(\d+)\z/ ? ($1, $2) : ($b, undef);
		($ap cmp $bp)
			|| (defined($an) && defined($bn) ? ($an <=> $bn) : (defined($an) <=> defined($bn)))
			|| ($a cmp $b)
	} @_;
}

sub _numeric {
	my ($value) = @_;
	return defined($value) && $value =~ /\A\s*[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?\s*\z/;
}

sub _close_reader {
	my ($fh, $path) = @_;
	close $fh or die "Cannot finish reading $path (status=$?): $!\n";
}

sub _publish {
	my ($temporary, $path) = @_;
	retry_unlink($path, label => "replace $path") if -e $path;
	retry_rename($temporary, $path, label => "publish $path")
		or die "Cannot publish $path\n";
}

# Writes a copy of MAGvsGC that yoke can read: rows grouped by MGS (MGS in
# natural order, rows of an MGS in input order), unassigned MAGs and MAGs yoke
# would filter out removed, the gene list column named other_Genes, gzipped.
# Returns counts, whether taxonomy columns are present, and (with
# collect_genes) the set of gene IDs in the kept MAGs.
sub normalise_magvsgc {
	my ($in, $out, %options) = @_;
	my $collect = $options{collect_genes} ? 1 : 0;
	make_path(dirname($out));
	my $rows = "$out.rows.$$.tmp";
	my ($I) = gzipopen($in, "MAGvsGC report", 1, 0);
	open my $R, '+>', $rows or die "Cannot create temporary file $rows: $!\n";

	my $header = <$I>;
	die "MAGvsGC report $in is empty\n" unless defined $header;
	$header =~ s/\r?\n\z//;
	my @columns = split /\t/, $header, -1;
	my %index;
	for my $i (0 .. $#columns) { $index{$columns[$i]} //= $i; }
	my ($geneColumn) = grep { lc($columns[$_]) eq 'other_genes' } 0 .. $#columns;
	my @missing = grep { !exists $index{$_} } @REQUIRED_MAGVSGC_COLUMNS;
	push @missing, 'other_Genes' unless defined $geneColumn;
	die "MAGvsGC report $in lacks column(s): " . join(', ', @missing) . "\n" if @missing;
	$columns[$geneColumn] = 'other_Genes';
	my $hasTaxonomy = (grep { exists $index{$_} } @TAXONOMY_COLUMNS) == @TAXONOMY_COLUMNS ? 1 : 0;

	my ($mgsColumn, $complColumn, $contaColumn) = @index{qw(MGS Completeness Contamination)};
	my (%rowsByMGS, %genes);
	my %count = (mags => 0, unassigned => 0, filtered => 0);
	while (my $line = <$I>) {
		$line =~ s/\r?\n\z//;
		next unless length $line;
		my @fields = split /\t/, $line, -1;
		my $mgs = $fields[$mgsColumn] // '';
		$mgs =~ s/^\s+|\s+$//g;
		if ($UNASSIGNED_MGS{$mgs}) { $count{unassigned}++; next; }
		my ($completeness, $contamination) = @fields[$complColumn, $contaColumn];
		unless (_numeric($completeness) && _numeric($contamination)
				&& $completeness > 0 && $contamination < 100) {
			$count{filtered}++;
			next;
		}
		my $offset = tell($R);
		print {$R} "$line\n" or die "Cannot write temporary file $rows: $!\n";
		push @{$rowsByMGS{$mgs}}, $offset;
		$count{mags}++;
		if ($collect) {
			for my $gene (split /,/, $fields[$geneColumn] // '') {
				$genes{$gene} = undef if length($gene) && $gene ne '?';
			}
		}
	}
	_close_reader($I, $in);

	my $temporary = "$out.$$.tmp.gz";
	my $O = gzipwrite($temporary, "normalised MAGvsGC report");
	print {$O} join("\t", @columns), "\n";
	my @mgsIDs = natural_mgs_sort(keys %rowsByMGS);
	for my $mgs (@mgsIDs) {
		for my $offset (@{$rowsByMGS{$mgs}}) {
			seek($R, $offset, 0) or die "Cannot seek in temporary file $rows: $!\n";
			my $line = <$R>;
			print {$O} $line or die "Cannot write $temporary: $!\n";
		}
	}
	close $O or die "Cannot finish writing $temporary (status=$?): $!\n";
	close $R;
	unlink $rows;
	_publish($temporary, $out);

	return {
		%count,
		mgs => scalar(@mgsIDs),
		taxonomy => $hasTaxonomy,
		($collect ? (genes => \%genes) : ()),
	};
}

# Rewrites a MATAFILER gene->function file as one "gene<TAB>function" row per
# function, keeping only the first (most specific) hierarchy level, and only
# genes in $keepGenes when it is given. Returns counts and the function IDs.
sub normalise_gene_annotations {
	my ($in, $out, $keepGenes) = @_;
	make_path(dirname($out));
	my ($I) = gzipopen($in, "gene annotation", 1, 0);
	my $temporary = "$out.$$.tmp.gz";
	my $O = gzipwrite($temporary, "normalised gene annotation");
	my (%functions, $genes, $rows);
	$genes = $rows = 0;
	while (my $line = <$I>) {
		next if $line =~ /^\s*(?:#|$)/;
		$line =~ s/\r?\n\z//;
		my ($gene, $value) = split /\t/, $line, 3;
		next unless defined($value) && length($gene);
		next if defined($keepGenes) && !exists($keepGenes->{$gene});
		my ($level) = split /;/, $value, 2;
		next unless defined $level;
		my (%seen, @ids);
		for my $id (split /[,|]/, $level) {
			$id =~ s/^\s+|\s+$//g;
			next if $id eq '' || $id eq '-' || $seen{$id}++;
			push @ids, $id;
		}
		next unless @ids;
		$genes++;
		for my $id (@ids) {
			print {$O} "$gene\t$id\n" or die "Cannot write $temporary: $!\n";
			$functions{$id} = undef;
			$rows++;
		}
	}
	_close_reader($I, $in);
	close $O or die "Cannot finish writing $temporary (status=$?): $!\n";
	_publish($temporary, $out);
	return { genes => $genes, rows => $rows, functions => \%functions };
}

# Writes "function<TAB>description" for every function in $functions. The
# description is the last non-empty field of the function's row in the first
# source file that names it; functions without one are described by their ID.
# Returns the number of functions that found a description.
sub write_function_descriptions {
	my ($sources, $functions, $out) = @_;
	my %description;
	for my $source (grep { defined($_) && -s $_ } @{$sources || []}) {
		my ($I) = gzipopen($source, "function descriptions", 1, 0);
		while (my $line = <$I>) {
			next if $line =~ /^\s*(?:#|$)/;
			$line =~ s/\r?\n\z//;
			my ($id, @fields) = split /\t/, $line;
			next unless defined $id;
			$id =~ s/^\s+|\s+$//g;
			next unless exists($functions->{$id}) && !exists($description{$id});
			my ($text) = grep { defined($_) && /\S/ } reverse @fields;
			$description{$id} = $text if defined $text;
		}
		_close_reader($I, $source);
	}
	make_path(dirname($out));
	my $temporary = "$out.$$.tmp";
	open my $O, '>', $temporary or die "Cannot write $temporary: $!\n";
	for my $id (sort keys %{$functions}) {
		my $text = $description{$id} // $id;
		$text =~ s/[\t\r\n]+/ /g;
		$text =~ tr/"/'/; # an unbalanced quote would make pandas merge rows
		$text =~ s/^\s+|\s+$//g;
		$text = $id unless length $text;
		print {$O} "$id\t$text\n" or die "Cannot write $temporary: $!\n";
	}
	close $O or die "Cannot finish writing $temporary: $!\n";
	_publish($temporary, $out);
	return scalar(keys %description);
}

1;
