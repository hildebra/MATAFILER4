#!/usr/bin/env bash
# Format eggNOG-mapper output (v2.1.x) for rtk sumMat: one "gene<TAB>annotation" file per category,
# written next to the input as eggNOGmapper_<category>.geneAss
#
# Usage: eggNOG_split.sh <MF.emapper.annotations[.gz]>
#
# Column reference (eggNOG-mapper v2.1.5-2.1.12):
# https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.12#user-content-Output_format
#   1 query  5 eggNOG_OGs  10 GOs  11 EC  12 KEGG_ko  13 KEGG_Pathway  14 KEGG_Module
#  19 CAZy  20 BiGG_Reaction  21 PFAMs
#
# Categories (rtk separators: ";" hierarchy, "," summed):
#   CAZy EC GO BIGG PFAM  column as is
#   KO                    KEGG KOs, "ko:" removed           (same layout as the diamond KGM table; used for modules)
#   KGM                   KOs;modules
#   KGP                   KOs;pathways (ko* ids, map* duplicates removed)
#   NOG                   first (root level) OG, "@taxid|name" removed
# Genes without an annotation ("-" or empty) are left out of that category.
set -euo pipefail

IN="${1:?usage: $0 <emapper.annotations[.gz]>}"
[[ -s "$IN" ]] || { echo "eggNOG_split: input '$IN' missing or empty" >&2; exit 1; }
OUTDIR="$(dirname "$IN")"
P="$OUTDIR/eggNOGmapper"
CATS=(CAZy EC GO BIGG PFAM KO KGM KGP NOG)

# outputs are written to .tmp files and only renamed if the whole input was processed;
# on failure the .tmp files are removed and existing outputs are left untouched
trap 'rm -f "${P}"_*.geneAss.tmp' EXIT
for c in "${CATS[@]}"; do : > "${P}_${c}.geneAss.tmp"; done

gzip -cdf -- "$IN" | awk -F'\t' -v OFS='\t' -v P="$P" '
	# keep only elements of a comma separated list that do not match re
	function drop(list, re,    n, a, i, out) {
		n = split(list, a, ","); out = ""
		for (i = 1; i <= n; i++) if (a[i] !~ re) out = out (out == "" ? "" : ",") a[i]
		return out
	}
	function has(x) { return x != "" && x != "-" }

	/^#/ { next }   # "#query" header and "##" comment/stat lines
	NF < 21 {
		printf("eggNOG_split: line %d has %d columns (expected >= 21)\n", NR, NF) > "/dev/stderr"
		bad = 1; exit 1
	}
	{
		g = $1
		if (has($19)) print g, $19 > (P "_CAZy.geneAss.tmp")
		if (has($11)) print g, $11 > (P "_EC.geneAss.tmp")
		if (has($10)) print g, $10 > (P "_GO.geneAss.tmp")
		if (has($20)) print g, $20 > (P "_BIGG.geneAss.tmp")
		if (has($21)) print g, $21 > (P "_PFAM.geneAss.tmp")

		if (has($12)) {
			ko = $12; gsub(/ko:/, "", ko)
			mod = has($14) ? $14 : ""
			pw  = has($13) ? drop($13, "^map") : ""
			print g, ko            > (P "_KO.geneAss.tmp")
			print g, ko ";" mod    > (P "_KGM.geneAss.tmp")
			print g, ko ";" pw     > (P "_KGP.geneAss.tmp")
		}

		if (has($5)) {
			split($5, og, ","); sub(/@.*/, "", og[1])
			print g, og[1] > (P "_NOG.geneAss.tmp")
		}
		n++
	}
	END {
		if (bad) exit 1
		if (n == 0) { print "eggNOG_split: no annotation rows found" > "/dev/stderr"; exit 1 }
		printf("eggNOG_split: %d annotated genes processed\n", n) > "/dev/stderr"
	}
'

for c in "${CATS[@]}"; do
	mv -f "${P}_${c}.geneAss.tmp" "${P}_${c}.geneAss"
done
