#!/usr/bin/env perl
# Per-MAG function tables for every MGS, made with yoke (github.com/falklab/yoke).
# Submitted by MGS.pl once the MGS and the gene catalogue annotations are final.
#
# For each scheme (e.g. KEGG) it writes to <outD>:
#   <scheme>/<MGS>.tsv.gz, <scheme>/presence.<MGS>.tsv.gz  yoke function_counts
#   <scheme>.tsv                                           mean count per MGS
#   modules/<scheme>/<MGS>.<set>.*                         yoke module_score (KEGG ortholog schemes)
#   summarised/<scheme>/<rank>/...                         yoke function_summarise (needs MAGvsGC taxonomy)
#   input/                                                 the yoke inputs prepared from MATAFILER files
#   function_tables.manifest.tsv                           what was made; written last
#
# Usage: funcTablesMGS.pl -GCd <gene catalogue> -binD <GCd>/Bin_<binner> -schemes KEGG,emKO [-modules BSB,GBM]
use warnings;
use strict;
use Getopt::Long qw( GetOptions );
use File::Path qw(make_path remove_tree);
use File::Spec;
use Mods::IO_Tamoc_progs qw(getProgPaths);
use Mods::FuncTables qw(
	function_table_scheme function_table_scheme_names parse_function_table_list resolve_function_table_scheme
	normalise_magvsgc normalise_gene_annotations write_function_descriptions
);
use Mods::WorkflowResilience qw(atomic_write_text);

my ($GCd, $binD, $outD, $schemeSpec, $moduleSpec) = ("", "", "", "", "BSB,GBM");
my $tablesPerCall = 200; #module_score is called on batches of tables: one shell argument is limited to 128 kB
GetOptions(
	"GCd=s"     => \$GCd,       #gene catalogue
	"binD=s"    => \$binD,      #MGS output dir (<GCd>/Bin_<binner>/), holds MAGvsGC.txt.gz
	"outD=s"    => \$outD,      #default <binD>/Annotation/FunctionTables/
	"schemes=s" => \$schemeSpec,#annotation schemes, see Mods/FuncTables.pm
	"modules=s" => \$moduleSpec,#yoke bundled module sets scored for KEGG ortholog schemes; 0 = none
) or die "Invalid funcTablesMGS.pl options\n";
die "Unexpected positional arguments: @ARGV\n" if @ARGV;
die "-GCd must be an existing directory\n" unless length($GCd) && -d $GCd;
die "-binD must be an existing directory\n" unless length($binD) && -d $binD;
my @schemes = parse_function_table_list($schemeSpec, 'scheme', function_table_scheme_names());
die "-schemes selects no annotation scheme\n" unless @schemes;
my @moduleSets = parse_function_table_list($moduleSpec, 'module set');
$outD = "$binD/Annotation/FunctionTables/" if $outD eq "";

my $magvsgc = "$binD/MAGvsGC.txt.gz";
die "MAGvsGC report not found: $magvsgc\n" unless -s $magvsgc;
my $yoke = getProgPaths("yoke");
my $rtk = getProgPaths("rare");
my @ranks = qw(Domain Phylum Class Order Family Genus Species);

#resolve every input before any work, so a missing file fails fast
my @resolved;
for my $name (@schemes) {
	my $scheme = resolve_function_table_scheme($GCd, $name);
	die "Annotation for scheme $name not found; looked for: "
		. join(', ', map { "$GCd/$_" } @{function_table_scheme($name)->{annotation}}) . "\n"
		unless defined $scheme->{annotation};
	push @resolved, $scheme;
}

my $inputD = "$outD/input/";
make_path($outD, $inputD);
my @manifest = ("scheme\tstep\tdetail\tpath");

#yoke module_score calls "rtk module" from PATH: point "rtk" at MATAFILER's rtk
my $shimD = "$inputD/bin/";
make_path($shimD);
atomic_write_text("$shimD/rtk", "#!/usr/bin/env bash\nexec $rtk \"\$@\"\n", label => 'rtk launcher for yoke');
chmod 0755, "$shimD/rtk" or die "Cannot make $shimD/rtk executable: $!\n";
$ENV{PATH} = "$shimD:$ENV{PATH}";

run("$yoke --version", "yoke version");

print "Preparing MAGvsGC for yoke\n";
my $normMAG = "$inputD/MAGvsGC.byMGS.txt.gz";
my $mag = normalise_magvsgc($magvsgc, $normMAG, collect_genes => 1);
my $nGenes = scalar keys %{$mag->{genes}};
print "MAGvsGC: $mag->{mags} MAGs in $mag->{mgs} MGS with $nGenes genes; "
	. "$mag->{unassigned} MAGs without MGS and $mag->{filtered} MAGs outside yoke's quality filter left out\n";
die "No MAG in $magvsgc belongs to an MGS\n" unless $mag->{mgs} > 0;
print "MAGvsGC has no taxonomy columns (Domain..Species): function_summarise is skipped\n"
	unless $mag->{taxonomy};

for my $scheme (@resolved) {
	my $name = $scheme->{name};
	print "\n== $name: $scheme->{annotation}\n";
	my $annoF = "$inputD/$name.geneAss.gz";
	my $descrF = "$inputD/$name.descr";
	my $anno = normalise_gene_annotations($scheme->{annotation}, $annoF, $mag->{genes});
	my $nFunctions = scalar keys %{$anno->{functions}};
	print "$anno->{genes} MAG genes with $anno->{rows} annotations to $nFunctions functions\n";
	my $described = write_function_descriptions([$scheme->{descriptions}], $anno->{functions}, $descrF);
	print "Descriptions: $described of $nFunctions functions described"
		. (defined $scheme->{descriptions} ? " from $scheme->{descriptions}" : "; no description file, IDs used") . "\n";

	#yoke overwrites files but leaves tables of MGS that no longer exist
	my $tableD = "$outD/$name/";
	remove_tree($tableD, "$outD/modules/$name", "$outD/summarised/$name");
	unlink "$outD/$name.tsv";
	run(join(' ', $yoke, 'function_counts', '-m', shq($normMAG), '-a', shq($annoF), '-d', shq($descrF),
		'-o', shq($outD), '-n', shq($name), '--compress'), "$name function counts");
	die "yoke function_counts did not write $outD/$name.tsv\n" unless -s "$outD/$name.tsv";
	my @tables = list_files($tableD, qr/\A(?!presence\.).+\.tsv\.gz\z/);
	die "yoke function_counts wrote no MGS table to $tableD\n" unless @tables;
	push @manifest, "$name\tfunction_counts\t" . scalar(@tables) . " MGS\t$outD/$name.tsv";

	if ($scheme->{modules}) {
		for my $set (@moduleSets) {
			my $modD = "$outD/modules/$name/";
			make_path($modD);
			for (my $i = 0; $i < @tables; $i += $tablesPerCall) {
				my $last = $i + $tablesPerCall - 1; $last = $#tables if $last > $#tables;
				run(join(' ', $yoke, 'module_score', '-i', (map { shq($_) } @tables[$i .. $last]),
					'-o', shq($modD), '-m', shq($set)), "$name $set module scores");
			}
			my @scores = list_files($modD, qr/\.\Q$set\E\.MODscore\z/);
			die "yoke module_score wrote " . scalar(@scores) . " of " . scalar(@tables)
				. " $set score files to $modD\n" unless @scores == @tables;
			push @manifest, "$name\tmodule_score\t$set\t$modD";
		}
	}

	if ($mag->{taxonomy}) {
		my $sumD = "$outD/summarised/$name/";
		run(join(' ', $yoke, 'function_summarise', '-m', shq($normMAG), '-c', shq($tableD),
			'-r', @ranks, '-o', shq($sumD)), "$name rank summaries");
		my @missing = grep { !-s "$sumD/$_.tsv" } @ranks;
		die "yoke function_summarise did not write rank table(s) @missing to $sumD\n" if @missing;
		push @manifest, "$name\tfunction_summarise\t" . join(',', @ranks) . "\t$sumD";
	}
}

#written last: its presence means every requested table was made
atomic_write_text("$outD/function_tables.manifest.tsv", join("\n", @manifest) . "\n",
	label => 'function table manifest');
print "\nFunction tables finished: $outD\n";
exit(0);


sub shq {
	my ($value) = @_;
	$value =~ s/'/'"'"'/g;
	return "'$value'";
}

sub run {
	my ($cmd, $label) = @_;
	print "Running $label:\n$cmd\n";
	system("bash", "-c", "set -o pipefail; $cmd") == 0
		or die "$label failed (exit " . ($? >> 8) . ")\n";
}

#files in $dir whose name matches $pattern, as full paths in name order
sub list_files {
	my ($dir, $pattern) = @_;
	opendir(my $dh, $dir) or die "Cannot read directory $dir: $!\n";
	my @files = map { File::Spec->catfile($dir, $_) } sort grep { $_ =~ $pattern && -f "$dir/$_" } readdir $dh;
	closedir $dh;
	return @files;
}
