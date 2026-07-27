use strict;
use warnings;
use Test::More;
use FindBin qw($Bin);
use lib "$Bin/..";
use Mods::SlurmAccounting qw(slurm_tree_memory_summary format_slurm_tree_memory_summary);

my @records = (
	{ job_id => 101, mgs => 'MGS.1', requested_mb => 10000 },
	{ job_id => 102, mgs => 'MGS.2', requested_mb => 20000 },
	{ job_id => 103, mgs => 'MGS.3', requested_mb => 8000 },
);
my $summary = slurm_tree_memory_summary(\@records, {
	runner => sub {
		return (join("\n",
			'101|COMPLETED|0:0|',
			'101.batch|COMPLETED|0:0|6G',
			'102|OUT_OF_MEMORY|0:125|',
			'102.batch|OUT_OF_MEMORY|0:125|21000M',
			'103|COMPLETED|0:0|',
		)."\n", 0);
	},
});

ok($summary->{available}, 'SLURM accounting output is accepted');
is($summary->{accounted}, 2, 'only jobs with MaxRSS contribute to utilization');
is($summary->{missing}, 1, 'missing MaxRSS is reported');
is_deeply([map { $_->{mgs} } @{$summary->{oom_jobs}}], ['MGS.2'],
	'OOM states are associated with their MGS');
cmp_ok(abs($summary->{mean_headroom_mb} - 1428), '<', 0.01,
	'mean signed requested-minus-used memory is calculated');
is($summary->{q95_headroom_mb}, 3856,
	'95th quantile uses signed memory headroom');
like(format_slurm_tree_memory_summary($summary),
	qr/OOM events: 1 \(MGS\.2 \[job 102\]\).*positive means under-used, negative means over-used/s,
	'report explains OOM identity and signed utilization');

done_testing();
