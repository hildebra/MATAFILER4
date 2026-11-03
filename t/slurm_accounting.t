use strict;
use warnings;
use Test::More;
use FindBin qw($Bin);
use lib "$Bin/..";
use Mods::SlurmAccounting qw(
	slurm_tree_memory_summary format_slurm_tree_memory_summary
	next_oom_retry_memory_mb slurm_oom_retry_plan
	slurm_job_id_from_dependency
);

is(next_oom_retry_memory_mb(10_240, 1_536_000), 20_480,
	'OOM retry memory doubles the previous request');
is(next_oom_retry_memory_mb(900_000, 1_536_000), 1_536_000,
	'OOM retry memory is clamped to the configured ceiling');
ok(!defined(next_oom_retry_memory_mb(1_536_000, 1_536_000)),
	'an OOM job already at the ceiling is not retried');

is(slurm_job_id_from_dependency('run123', 'run'), 123,
	'Slurm job IDs are normalized from tagged submission dependencies');
ok(!defined(slurm_job_id_from_dependency('submission-failed', 'run')),
	'non-job submission dependencies are rejected');

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
my $retry_plan = slurm_oom_retry_plan(\@records, 30_000, {
	runner => sub {
		return (join("\n",
			'101|COMPLETED|0:0|',
			'102|OUT_OF_MEMORY|0:125|',
			'103|FAILED|1:0|',
		)."\n", 0);
	},
});
is_deeply([sort { $a <=> $b } keys %{$retry_plan->{by_job_id}}], [102],
	'the shared retry planner selects only accounting-confirmed OOM jobs');
is($retry_plan->{by_job_id}{102}{next_mb}, 30_000,
	'the shared retry planner doubles memory and clamps it to the ceiling');
ok(!$retry_plan->{by_job_id}{102}{ceiling_reached},
	'a clamped increase remains retryable');

my $ceiling_plan = slurm_oom_retry_plan(
	[{ job_id => 104, requested_mb => 30_000 }], 30_000, {
		runner => sub { return ("104|OUT_OF_MEMORY|0:125|\n", 0); },
	},
);
ok($ceiling_plan->{by_job_id}{104}{ceiling_reached},
	'the shared retry planner marks OOM jobs already at the ceiling');
ok(!defined($ceiling_plan->{by_job_id}{104}{next_mb}),
	'a ceiling-limited OOM plan has no unsafe retry request');

my ($settle_calls, $settle_now) = (0, 0);
my $settled_plan = slurm_oom_retry_plan(
	[{ job_id => 105, requested_mb => 10_000 }], 40_000, {
		runner => sub {
			$settle_calls++;
			return $settle_calls == 1
				? ("", 0) : ("105|OUT_OF_MEMORY|0:125|\n", 0);
		},
		clock => sub { $settle_now },
		sleeper => sub { $settle_now += $_[0] },
		settle_retry_seconds => 2,
		settle_maximum_seconds => 10,
	},
);
is($settle_calls, 2,
	'the shared retry planner waits for briefly delayed accounting rows');
is($settled_plan->{by_job_id}{105}{next_mb}, 20_000,
	'a delayed OOM accounting row still escalates the next retry');

cmp_ok(abs($summary->{mean_headroom_mb} - 1428), '<', 0.01,
	'mean signed requested-minus-used memory is calculated');
is($summary->{q95_headroom_mb}, 3856,
	'95th quantile uses signed memory headroom');
like(format_slurm_tree_memory_summary($summary),
	qr/OOM events: 1 \(MGS\.2 \[job 102\]\).*positive means under-used, negative means over-used/s,
	'report explains OOM identity and signed utilization');

done_testing();
