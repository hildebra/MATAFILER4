use strict;
use warnings;

use File::Basename qw(dirname);
use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use Test::More;

use lib File::Spec->catdir(File::Spec->rel2abs('.'));
use Mods::RibosomeState qw(prepare_ribosome_rerun);
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
my $sample_root = File::Spec->catdir($root, 'S1');
my $central_root = File::Spec->catdir($root, 'central results', 'RiboFind');
my $sample = 'S1';
populate_results($sample_root, $central_root, $sample);

my $assignment_result = prepare_ribosome_rerun(
	sample_root => $sample_root,
	central_root => $central_root,
	sample => $sample,
	redo_profile => 0,
	redo_assignment => 1,
);
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
