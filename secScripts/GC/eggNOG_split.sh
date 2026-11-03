#!/usr/bin/env bash
# Format eggNOG-mapper output (v3.x with eggNOG 7, or v2.1.x with eggNOG 5) for rtk sumMat:
# one "gene<TAB>annotation" file per category, written next to the input as eggNOGmapper_<category>.geneAss
#
# Usage: eggNOG_split.sh <MF.emapper.annotations[.gz]>
#
# Columns are taken by name from the "#query" header line. Without a recognised header the layout is
# inferred: 22 columns ending in the v3 confidence code (h/m/l/-) = v3, otherwise v2 positions.
#   v3 (22 columns, https://github.com/eggnogdb/eggnog-mapper/blob/v3.0.0-beta6/USAGE.md#the-annotations-columns):
#      1 query  5 eggNOG_OGs  8 COG_category  10 GOs  11 EC  12 KEGG_ko  13 KEGG_Pathway  14 KEGG_Module
#     19 CAZy  20 BiGG_Reaction  21 PFAMs  22 annotation_confidence
#   v2.1.5-2.1.12 (21 columns): as v3 for columns 10-21, but 7 COG_category (letters)
# All rows must have as many columns as the header (or the first row): a merge of chunks annotated by
# different eggNOG-mapper versions is rejected.
#
# Categories (rtk separators: ";" hierarchy, "," summed), same content for v2 and v3:
#   CAZy GO BIGG          column as is
#   EC                    "ec:" prefix (v3) removed
#   PFAM                  domain names; v3 domain coordinates ("_<start>_<end>") and repeated domains removed
#   KO                    KEGG KOs, "ko:" (v2) removed     (same layout as the diamond KGM table; used for modules)
#   KGM                   KOs;modules
#   KGP                   KOs;pathways as ko* ids (v2 map* duplicates removed, v3 bare numbers prefixed "ko")
#   NOG                   v2: first (root level) OG, "@taxid|name" removed
#                         v3: the COG in COG_category; genes without a COG get their first (broadest) eggNOG 7 OG,
#                             as name@taxid.cluster ("|" is an rtk separator)
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
	function fail(msg) { printf("eggNOG_split: %s\n", msg) > "/dev/stderr"; bad = 1; exit 1 }
	# keep only elements of a comma separated list that do not match re
	function drop(list, re,    n, a, i, out) {
		n = split(list, a, ","); out = ""
		for (i = 1; i <= n; i++) if (a[i] !~ re) out = out (out == "" ? "" : ",") a[i]
		return out
	}
	function has(x) { return x != "" && x != "-" }
	# v2 positions; v3 differs only in COG_category and the trailing confidence column
	function layout(v3) {
		isV3 = v3
		C["query"] = 1; C["eggNOG_OGs"] = 5; C["COG_category"] = v3 ? 8 : 7
		C["GOs"] = 10; C["EC"] = 11; C["KEGG_ko"] = 12; C["KEGG_Pathway"] = 13; C["KEGG_Module"] = 14
		C["CAZy"] = 19; C["BiGG_Reaction"] = 20; C["PFAMs"] = 21
	}
	function fromHeader(    i, k, n, need) {
		for (i = 1; i <= NF; i++) { k = $i; sub(/^#/, "", k); H[k] = i }
		if (!("eggNOG_OGs" in H)) return 0   # not an eggNOG-mapper header: infer from the rows
		n = split("query eggNOG_OGs COG_category GOs EC KEGG_ko KEGG_Pathway KEGG_Module CAZy BiGG_Reaction PFAMs", need, " ")
		for (i = 1; i <= n; i++) {
			if (!(need[i] in H)) fail("header lacks column " need[i])
			C[need[i]] = H[need[i]]
		}
		isV3 = ("annotation_confidence" in H) || ("tax_ceiling" in H)
		return 1
	}
	# v3: drop "ec:" prefixes
	function ecs(x,    n, a, i, out) {
		n = split(x, a, ","); out = ""
		for (i = 1; i <= n; i++) { sub(/^ec:/, "", a[i]); out = out (out == "" ? "" : ",") a[i] }
		return out
	}
	# ko* pathway ids: v2 lists ko and map duplicates, v3 bare numbers
	function pathways(x,    n, a, i, out) {
		n = split(x, a, ","); out = ""
		for (i = 1; i <= n; i++) {
			if (a[i] ~ /^map/) continue
			if (a[i] ~ /^[0-9]+$/) a[i] = "ko" a[i]
			out = out (out == "" ? "" : ",") a[i]
		}
		return out
	}
	# v3 PFAMs are "name_start_end" per domain hit: names only, each once
	function pfams(x,    n, a, i, out, seen) {
		if (!isV3) return x
		n = split(x, a, ","); out = ""; split("", seen)
		for (i = 1; i <= n; i++) {
			sub(/_[0-9]+_[0-9]+$/, "", a[i])
			if (a[i] in seen) continue
			seen[a[i]] = 1; out = out (out == "" ? "" : ",") a[i]
		}
		return out
	}
	function nog(ogs, cog,    og, n, a) {
		if (isV3) {
			n = split(cog, a, ",")
			if (a[1] ~ /^[A-Za-z]*COG[0-9]+$/) return a[1]
			if (!has(ogs)) return ""
			split(ogs, og, ","); gsub(/!/, "", og[1]); gsub(/[|;]/, ".", og[1])
			return og[1]
		}
		if (!has(ogs)) return ""
		split(ogs, og, ","); sub(/@.*/, "", og[1])
		return og[1]
	}

	/^##/ { next }   # comment/stat lines
	/^#/ {           # "#query" header (once per file; merged chunks keep only the first)
		if (hdrLine == "") {
			hdrLine = $0; ncol = NF
			named = fromHeader()
		} else if ($0 != hdrLine) {
			fail("line " NR ": header differs from the first header (annotations from different eggNOG-mapper versions?)")
		}
		next
	}
	!done {          # first data row: column layout
		if (!named) {
			if (ncol == 0) ncol = NF
			layout(NF == 22 && $22 ~ /^[hml-]+$/)
		}
		done = 1
	}
	NF != ncol || NF < 21 {
		fail(sprintf("line %d has %d columns (expected %d; a merge of different eggNOG-mapper versions?)", NR, NF, (ncol >= 21 ? ncol : 21)))
	}
	{
		g = $(C["query"])
		if (has($(C["CAZy"]))) print g, $(C["CAZy"]) > (P "_CAZy.geneAss.tmp")
		if (has($(C["EC"]))) print g, ecs($(C["EC"])) > (P "_EC.geneAss.tmp")
		if (has($(C["GOs"]))) print g, $(C["GOs"]) > (P "_GO.geneAss.tmp")
		if (has($(C["BiGG_Reaction"]))) print g, $(C["BiGG_Reaction"]) > (P "_BIGG.geneAss.tmp")
		if (has($(C["PFAMs"]))) print g, pfams($(C["PFAMs"])) > (P "_PFAM.geneAss.tmp")

		if (has($(C["KEGG_ko"]))) {
			ko = $(C["KEGG_ko"]); gsub(/ko:/, "", ko)
			mod = has($(C["KEGG_Module"])) ? $(C["KEGG_Module"]) : ""
			pw  = has($(C["KEGG_Pathway"])) ? pathways($(C["KEGG_Pathway"])) : ""
			print g, ko            > (P "_KO.geneAss.tmp")
			print g, ko ";" mod    > (P "_KGM.geneAss.tmp")
			print g, ko ";" pw     > (P "_KGP.geneAss.tmp")
		}

		og = nog($(C["eggNOG_OGs"]), $(C["COG_category"]))
		if (og != "") print g, og > (P "_NOG.geneAss.tmp")
		n++
	}
	END {
		if (bad) exit 1
		if (n == 0) { print "eggNOG_split: no annotation rows found" > "/dev/stderr"; exit 1 }
		printf("eggNOG_split: %d annotated genes processed (eggNOG-mapper %s output)\n", n, isV3 ? "v3" : "v2") > "/dev/stderr"
	}
'

for c in "${CATS[@]}"; do
	mv -f "${P}_${c}.geneAss.tmp" "${P}_${c}.geneAss"
done
