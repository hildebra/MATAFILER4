#!/usr/bin/env python3
"""Small local fixtures for the MATAF4 algorithm audit; no external aligners.

Run from any directory with Python 3, Perl, and gzip installed. All fixture
inputs/outputs live in a temporary directory. Results describe actual behavior;
this is an audit demonstrator, not a passing regression-test specification.
"""
import gzip
import hashlib
import json
import re
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
PERL = ['perl', f'-I{ROOT}', f'-I{ROOT / "t/lib"}', '-MMFTestConfig']


def invoke(script, *args, stdin=None):
    return subprocess.run(PERL + [str(ROOT / script), *map(str, args)],
                          input=stdin, text=True, capture_output=True, cwd=ROOT)


def gzwrite(path, text):
    with gzip.open(path, 'wt') as out:
        out.write(text)


def gzread(path):
    with gzip.open(path, 'rt') as src:
        return src.read()


def subroutine(source, name):
    match = re.search(r'^sub ' + re.escape(name) + r'\b[^\n]*\{.*?^\}',
                      source, re.S | re.M)
    if not match:
        raise RuntimeError(f'Cannot extract {name}')
    return match.group()


def hit(query, subject, length=50, score=100, start=1, end=None):
    return '\t'.join(map(str, [query, subject, 90, length, 0, 0, 1,
                               3 * length, start, end or start + length - 1,
                               '1e-20', score])) + '\n'


with tempfile.TemporaryDirectory(prefix='mataf4-algorithm-audit-') as tmp:
    temp = Path(tmp)
    results = {}

    def functional(name, text, fraction=0):
        work = temp / name
        work.mkdir()
        blast = work / 'hits.srt.gz'
        gzwrite(blast, text)
        lengths = work / 'lengths.tsv'
        lengths.write_text('A\t1000\nB\t100\n')
        proc = invoke('secScripts/functions/parseBlastFunct2.pl',
                      '-i', blast, '-DB', 'TEST', '-tmp', work / 'scratch',
                      '-LF', lengths, '-eval', '1e-7', '-minFractQueryCov', fraction)
        if proc.returncode:
            raise RuntimeError(proc.stderr + proc.stdout)
        tables = work / 'CNT_1e-7_20'
        return {norm: gzread(tables / f'TESTparse.TEST.ALL.{norm}.gene.cnts.gz').strip()
                for norm in ('cnt', 'GLN')}

    text = hit('read/1', 'A', score=50) + hit('read/1', 'B', score=100)
    results['wrong_subject_length'] = {
        'actual': functional('length', text),
        'expected_GLN': 'B\t0.5',
        'actual_with_subject_fraction_0_4': functional('length-filter', text, 0.4),
        'expected_with_subject_fraction_0_4': 'B retained (50 / 100 > 0.4)',
    }
    results['mate_2_only'] = {
        'mate_1_only': functional('mate1', hit('read/1', 'B')),
        'mate_2_only': functional('mate2', hit('read/2', 'B')),
        'expected': 'Both orientations retain B with cnt=1 and GLN=0.5',
    }
    results['unsuffixed_single_read'] = {
        'actual': functional('single', hit('read', 'B')),
        'expected_cnt': 'B\t1 for a single-end read (no merged mate)',
    }

    source = (ROOT / 'secScripts/functions/parseBlastFunct2.pl').read_text()
    merged_code = 'use JSON::PP; use Mods::TamocFunc qw(uniq); use Mods::FuncTools qw(mergeBlastPair);\n' + subroutine(source, 'combineBlasts')
    merged_code += r'''
my @a = ("read/1", "B", 90, 200, 0, 0, 1, 600, 1, 200, 1e-20, 200);
my @b = ("read/2", "B", 90, 51, 0, 0, 1, 153, 50, 100, 1e-20, 51);
my $out = combineBlasts([\@a], [\@b]);
print encode_json({length=>$out->[0][3], bitscore=>$out->[0][11]});
'''
    proc = subprocess.run(PERL + ['-e', merged_code], text=True, capture_output=True, cwd=ROOT, check=True)
    results['contained_pair_overlap'] = {'actual': json.loads(proc.stdout),
                                         'expected_length': 200,
                                         'correct_overlap': 51}

    sam = ('pass\t0\tref\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\tNM:i:0\n'
           'fail\t0\tref\t1\t0\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\tNM:i:0\n')
    filtering = invoke('secScripts/assemblies/bamFilter.pl', stdin=sam)
    if filtering.returncode:
        raise RuntimeError(filtering.stderr)
    log = temp / 'map.log'
    log.write_text('[M::worker_pipeline::0.1] mapped 2 sequences\n' + filtering.stderr)
    main = (ROOT / 'MATAF4.pl').read_text()
    stats_code = r'''
use strict; use warnings; use JSON::PP; use Mods::StatsLogReader qw(parse_bam_filter_counters);
my %locStats;
my $logpath = shift;
sub read_stats_log_excerpt { open my $f, '<', $logpath or die $!; local $/; return <$f>; }
'''
    stats_code += subroutine(main, 'bwtLogRd') + '\n' + subroutine(main, 'getMapStats')
    stats_code += '\nprint encode_json(getMapStats("unused"));\n'
    proc = subprocess.run(PERL + ['-e', stats_code, str(log)], text=True, capture_output=True, cwd=ROOT, check=True)
    results['bam_filter_log_contract'] = {
        'actual': json.loads(proc.stdout),
        'expected': {'AlignedReads': 1, 'OverallAlignment': 50},
    }

    hierarchies = temp / 'hierarchies'
    hierarchies.mkdir()
    header = 'read\tdomain\n'
    (hierarchies / 'populated.hiera.txt').write_text(header + 'r1\tBacteria\n')
    (hierarchies / 'empty.hiera.txt').write_text('')
    proc = invoke('secScripts/miTag/miTagTaxTable.pl', 'domain', temp / 'merged', hierarchies)
    results['empty_ribo_profile'] = {'exit_code': proc.returncode, 'stderr': proc.stderr.strip(),
                                      'expected': 'Successful table with populated=1 and empty=0'}
    (hierarchies / 'empty.hiera.txt').write_text(header)
    proc = invoke('secScripts/miTag/miTagTaxTable.pl', 'domain', temp / 'merged-header', hierarchies)
    results['header_only_ribo_profile'] = {
        'exit_code': proc.returncode,
        'actual': gzread(temp / 'merged-header.domain.txt.gz').strip(),
        'expected': 'Empty sample retained as an explicit zero column',
    }

    kmers = temp / 'input.4kmer.gz'
    gzwrite(kmers, 'Contig\tAAAA\n' + ''.join(f'ctg_{i}\t{i}\n' for i in range(1, 9)))
    proc = invoke('secScripts/composition/kmer_Ngenes.pl', kmers, 5)
    if proc.returncode:
        raise RuntimeError(proc.stderr)
    results['rolling_gene_windows'] = {
        'actual': gzread(temp / 'input.4kmer.pm5.gz').strip().splitlines(),
        'expected': 'Each of ctg_1 ... ctg_8 occurs once, in order',
    }

    fasta = temp / 'ambiguous.fa'
    fasta.write_text('>ambiguous\nNNNN\n>valid\nGCGC\n')
    proc = invoke('secScripts/composition/calcGC.pl', fasta, temp / 'gc.tsv')
    results['ambiguous_gc_row'] = {'exit_code': proc.returncode,
                                  'actual': (temp / 'gc.tsv').read_text().strip(),
                                  'expected': 'Separate rows for ambiguous and valid records'}
    fasta.write_text('>valid\n' + 'ACGT' * 30 + '\n>ambiguous\n' + 'N' * 120 + '\n')
    proc = invoke('secScripts/composition/calc.kmerfreq.pl', '-i', fasta, '-m', 100, '-o', temp / 'kmer.tsv')
    results['ambiguous_final_kmer'] = {'exit_code': proc.returncode,
                                      'stderr': proc.stderr.strip(),
                                      'expected': 'A zero-informative sequence does not cause division by zero'}

    motus_dir = temp / 'mOTU2'
    motus_dir.mkdir()
    taxonomy = ('d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;'
                'f__Lactobacillaceae;g__Lactobacillus;s__Lactobacillus testii')
    for sample, count in [('A', 12), ('B', 3)]:
        gzwrite(motus_dir / f'{sample}.motu2.tab.gz',
                f'#tool_version=4.1.0\nmOTU\tTaxonomy\t{sample}\n'
                f'mOTUv4.0_000001\t{taxonomy}\t{count}\n')
    proc = invoke('secScripts/composition/mrgMotu2.pl', motus_dir, 2)
    if proc.returncode:
        raise RuntimeError(proc.stderr)
    results['motus_repeated_taxon_across_samples'] = {
        'motu_matrix': (temp / 'm2.motu.txt').read_text().strip(),
        'kingdom_matrix': (temp / 'm2.kingdom.txt').read_text().strip(),
        'expected_kingdom': 'Bacteria\t12\t3',
    }

    coverage_dir = temp / 'coverage'
    coverage_dir.mkdir()
    coverage = coverage_dir / 'sample-smd.bam.coverage.gz'
    gzwrite(coverage, 'sample__C1_L=1000=\t0\t1000\t10\n')
    gff = coverage_dir / 'genes.gff'
    gff.write_text('sample__C1_L=1000=\tProdigal\tCDS\t1\t300\t.\t+\t0\tID=1_1;partial=00\n')
    calculator = ROOT / 'bin/rdCover'
    proc = subprocess.run([str(calculator), str(coverage), str(gff), '100'],
                          capture_output=True, text=True, cwd=ROOT, check=True)
    completion_code = ('use JSON::PP; use Mods::GenoMetaAss qw(coverage_derivative_paths); '
                       'my $cov = shift; my %accepted; '
                       'for my $suffix (qw(pergene percontig median.percontig)) { '
                       '$accepted{$suffix} = (grep { -s $_ } @{coverage_derivative_paths($cov,$suffix)}) '
                       '? JSON::PP::true : JSON::PP::false; } print encode_json(\\%accepted);')
    completion = subprocess.run(PERL + ['-e', completion_code, str(coverage)],
                                capture_output=True, text=True, cwd=ROOT, check=True)
    results['secondary_coverage_output_names'] = {
        'calculator_sha256': hashlib.sha256(calculator.read_bytes()).hexdigest(),
        'actual_outputs': sorted(p.name for p in coverage_dir.iterdir() if p not in (coverage, gff)),
        'accepted_by_current_caller': json.loads(completion.stdout),
        'mean_gene_depth': (coverage_dir / 'sample-smd.bam.coverage.pergene').read_text().strip(),
        'expected_mean_gene_depth': 10,
    }

    median_code = 'use JSON::PP; use Mods::GenoMetaAss qw(median); print encode_json([median(1,3), median(1,3,9)]);'
    proc = subprocess.run(PERL + ['-e', median_code], text=True, capture_output=True, cwd=ROOT, check=True)
    results['shared_median_after_refactor'] = {'actual': json.loads(proc.stdout), 'expected': [2, 3]}

    print(json.dumps(results, indent=2))
