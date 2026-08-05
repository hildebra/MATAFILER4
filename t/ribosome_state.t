use strict;
use warnings;

use File::Basename qw(dirname);
use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use Test::More;

use lib File::Spec->catdir(File::Spec->rel2abs('.'));
use Mods::RibosomeState qw(
	normalise_ribosome_request
	prepare_ribosome_rerun
	ribosome_completion_evidence
);
use Mods::SampleCompletion qw(sample_completion_path);

sub write_file {
	my ($path, $content) = @_;
	make_path(dirname($path));
	open my $fh, '>', $path or die "Cannot write $path: $!";
	print {$fh} defined($content) ? $content : "result\n";
	close $fh or die "Cannot close $path: $!";
}

sub populate_results {
	my ($sample_root, $central_root, $sample) = @_;
	my $ribo_root = File::Spec->catdir($sample_root, 'ribos');
	my $lca_root = File::Spec->catdir($ribo_root, 'ltsLCA');
	write_file(File::Spec->catfile($ribo_root, 'SSU_pull.sto'));
	write_file(File::Spec->catfile($ribo_root, 'LSU_pull.sto'));
	write_file(File::Spec->catfile($lca_root, 'Assigned.sto'));
	write_file(File::Spec->catfile($lca_root, 'SSU_ass.sto'));
	write_file(File::Spec->catfile($lca_root, 'LSU_ass.sto'));
	write_file(sample_completion_path($sample_root), "closed\n");

	for my $tag (qw(SSU LSU)) {
		for my $suffix (qw(r1.fq.gz r2.fq.gz fq.gz)) {
			write_file(File::Spec->catfile(
				$ribo_root, 'reads_'.$tag.'.'.$suffix,
			), "gzip container\n");
		}
		my $source = File::Spec->catfile(
			$lca_root, $tag."riboRun_bl.hiera.txt.gz",
		);
		write_file($source);
		my $tag_dir = File::Spec->catdir($central_root, $tag);
		make_path($tag_dir);
		my $sample_link = File::Spec->catfile(
			$tag_dir, "$sample.$tag.hiera.txt.gz",
		);
		symlink $source, $sample_link
			or die "Cannot create test symlink $sample_link: $!";
		write_file(File::Spec->catfile(
			$tag_dir, "other.$tag.hiera.txt.gz",
		));
		write_file(File::Spec->catfile(
			$central_root, "$tag.miTag.domain.txt",
		));
		write_file(File::Spec->catfile(
			$central_root, "$tag.cnt.stone",
		));
	}
}

my $root = tempdir(CLEANUP => 1);
my %profile_redo = (
	DoRibofind => 0,
	RedoRiboFind => 1,
	RedoRiboAssign => 0,
);
ok(normalise_ribosome_request(\%profile_redo),
	'profile redo enables ribosomal profiling');
is($profile_redo{DoRibofind}, 1,
	'profile redo cannot fall through the no-work completion path');

my %assignment_redo = (
	DoRibofind => 0,
	RedoRiboFind => 0,
	RedoRiboAssign => 1,
);
ok(normalise_ribosome_request(\%assignment_redo),
	'assignment redo enables ribosomal profiling');

my %ordinary_request = (
	DoRibofind => 0,
	RedoRiboFind => 0,
	RedoRiboAssign => 0,
);
ok(!normalise_ribosome_request(\%ordinary_request),
	'ordinary runs do not enable ribosomal profiling');

my $sample_root = File::Spec->catdir($root, 'S1');
my $central_root = File::Spec->catdir($root, 'central results', 'RiboFind');
my $sample = 'S1';
populate_results($sample_root, $central_root, $sample);

my $evidence = ribosome_completion_evidence(
	sample_root => $sample_root, requested => 1, assembly_requested => 0,
);
ok($evidence->{complete},
	'RiboFind completion requires both extraction and taxonomy artifacts');
ok($evidence->{ssu_hierarchy_complete} && $evidence->{lsu_hierarchy_complete},
	'compressed SSU and LSU hierarchy outputs are accepted');
my $ssu_profile_output = File::Spec->catfile(
	$sample_root, 'ribos', 'reads_SSU.fq.gz',
);
unlink $ssu_profile_output
	or die "Cannot remove $ssu_profile_output: $!";
$evidence = ribosome_completion_evidence(
	sample_root => $sample_root, requested => 1, assembly_requested => 0,
);
ok(!$evidence->{profile_complete} && !$evidence->{complete},
	'a profile stone cannot mask a missing extracted-read output');
write_file($ssu_profile_output, "gzip container\n");
my $lsu_hierarchy = File::Spec->catfile(
	$sample_root, 'ribos', 'ltsLCA', 'LSUriboRun_bl.hiera.txt.gz',
);
write_file($lsu_hierarchy, '');
$evidence = ribosome_completion_evidence(
	sample_root => $sample_root, requested => 1, assembly_requested => 0,
);
ok($evidence->{complete},
	'a zero-hit hierarchy remains complete when its producer stones are present');
unlink $lsu_hierarchy or die "Cannot remove $lsu_hierarchy: $!";
$evidence = ribosome_completion_evidence(
	sample_root => $sample_root, requested => 1, assembly_requested => 0,
);
ok(!$evidence->{taxonomy_complete} && !$evidence->{complete},
	'missing LSU hierarchy output makes RiboFind incomplete');
write_file($lsu_hierarchy);

$evidence = ribosome_completion_evidence(
	sample_root => $sample_root, requested => 1, assembly_requested => 1,
);
ok(!$evidence->{profile_complete},
	'requested ribosomal assembly requires its own completion stone');
write_file(File::Spec->catfile(
	$sample_root, 'ribos', 'Ass', 'allAss.sto',
));
$evidence = ribosome_completion_evidence(
	sample_root => $sample_root, requested => 1, assembly_requested => 1,
);
ok($evidence->{complete},
	'complete extraction, assembly, and taxonomy evidence closes RiboFind');

my $assignment_result = prepare_ribosome_rerun(
	sample_root => $sample_root,
	central_root => $central_root,
	sample => $sample,
	redo_profile => 0,
	redo_assignment => 1,
);
ok(!$assignment_result->{profile},
	'assignment rerun does not require a new extraction');
ok($assignment_result->{assignment},
	'assignment rerun explicitly requires replacement LCA work');
ok($assignment_result->{removed} > 0,
	'assignment rerun reports removed results');
ok(-d File::Spec->catdir($sample_root, 'ribos'),
	'assignment rerun retains extracted ribosomal reads');
ok(!-e File::Spec->catdir($sample_root, 'ribos', 'ltsLCA'),
	'assignment rerun removes the old LCA results');
ok(-e File::Spec->catfile($sample_root, 'ribos', 'SSU_pull.sto'),
	'assignment rerun retains the RiboFind extraction stones');
ok(!-e sample_completion_path($sample_root),
	'assignment rerun invalidates the sample completion sentinel');

for my $tag (qw(SSU LSU)) {
	my $sample_link = File::Spec->catfile(
		$central_root, $tag, "$sample.$tag.hiera.txt.gz",
	);
	ok(!-e $sample_link && !-l $sample_link,
		"$tag assignment rerun removes the per-sample central result");
	ok(!-e File::Spec->catfile($central_root, "$tag.miTag.domain.txt"),
		"$tag assignment rerun invalidates the merged table");
	ok(!-e File::Spec->catfile($central_root, "$tag.cnt.stone"),
		"$tag assignment rerun invalidates the merged-table count");
	ok(-e File::Spec->catfile(
		$central_root, $tag, "other.$tag.hiera.txt.gz",
	), "$tag assignment rerun retains other samples");
}

populate_results($sample_root, $central_root, $sample);
my $profile_result = prepare_ribosome_rerun(
	sample_root => $sample_root,
	central_root => $central_root,
	sample => $sample,
	redo_profile => 1,
	redo_assignment => 0,
);
ok($profile_result->{assignment},
	'profile rerun explicitly requires replacement LCA work');
ok($profile_result->{profile}, 'profile rerun records extraction invalidation');
ok(!-e File::Spec->catdir($sample_root, 'ribos'),
	'profile rerun removes the complete sample-local RiboFind directory');
ok(!-e sample_completion_path($sample_root),
	'profile rerun invalidates the completion sentinel');

write_file(File::Spec->catfile($sample_root, 'ribos', 'keep.txt'));
write_file(sample_completion_path($sample_root), "closed\n");
my $no_redo_result = prepare_ribosome_rerun(
	sample_root => $sample_root,
	central_root => $central_root,
	sample => $sample,
	redo_profile => 0,
	redo_assignment => 0,
);
is($no_redo_result->{removed}, 0, 'ordinary runs perform no redo cleanup');
ok(-e File::Spec->catfile($sample_root, 'ribos', 'keep.txt'),
	'ordinary runs retain RiboFind results');
ok(-e sample_completion_path($sample_root),
	'ordinary runs retain the completion sentinel');

done_testing;
