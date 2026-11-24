#!/bin/bash
# makeAccessionMap.sh - write a MATAFILER4 mapping file for public reads.
#
# Usage: makeAccessionMap.sh <study accession | ENA filereport TSV> <RunID> <OutPath> > study.map
#
# A study accession is resolved through the ENA Portal API (needs internet).
# A TSV must carry at least run_accession, sample_accession, library_strategy,
# library_selection and read_count.
#
# One map row per BioSample; all its kept runs are listed comma separated in
# ENAdownload, so MATAF4 downloads them into one sample. Runs are excluded,
# and reported on stderr, when library_strategy is AMPLICON, when
# library_selection is PCR (amplicons mislabelled as WGS or OTHER), or when
# ENA reports 0 reads. Never select by library_strategy WGS: shotgun runs are
# often submitted as OTHER or WGA.
set -euo pipefail
[ $# -eq 3 ] || { sed -n '2,15p' "$0" >&2; exit 2; }
src=$1; runid=$2; outpath=$3
if [ -f "$src" ]; then
	report=$src
else
	report=$(mktemp); trap 'rm -f "$report"' EXIT
	curl -sSf "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$src&result=read_run&fields=run_accession,sample_accession,library_strategy,library_selection,read_count&format=tsv" > "$report"
fi
printf '#SmplID\tENAdownload\tAssmblGrps\n#RunID\t%s\n#OutPath\t%s\n' "$runid" "$outpath"
awk -F'\t' '
NR==1 {
	for (i = 1; i <= NF; i++) c[$i] = i
	split("run_accession sample_accession library_strategy library_selection read_count", req, " ")
	for (i in req) if (!(req[i] in c)) { print "missing column " req[i] > "/dev/stderr"; bad = 1; exit 2 }
	next
}
{
	run = $c["run_accession"]; smp = $c["sample_accession"]
	strat = $c["library_strategy"]; sel = $c["library_selection"]; reads = $c["read_count"] + 0
	if (strat == "AMPLICON") { amp++; next }
	if (sel == "PCR")        { pcr++; print "excluded (library_selection PCR, strategy " strat "): " run > "/dev/stderr"; next }
	if (reads == 0)          { empty++; print "excluded (0 reads): " run > "/dev/stderr"; next }
	if (smp in runs) runs[smp] = runs[smp] "," run
	else { order[++n] = smp; runs[smp] = run }
	kept++
}
END {
	if (bad) exit 2
	for (i = 1; i <= n; i++) print order[i] "\t" runs[order[i]] "\t" order[i]
	printf "kept %d run(s) in %d sample(s); excluded AMPLICON %d, PCR %d, empty %d\n", kept, n, amp, pcr, empty > "/dev/stderr"
	if (n == 0) exit 1
}' "$report"
