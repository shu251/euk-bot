#!/usr/bin/env bash
# Summarize Kaiju per-read output to superkingdom in kaiju2table's 5-column format.
# kaiju2table >= 1.10 rejects `-r superkingdom`, so the rollup is done here: each classified read's
# taxid is walked up nodes.dmp to its superkingdom ancestor. Reads assigned above superkingdom
# (e.g. "cellular organisms") are reported as "cannot be assigned to a superkingdom".
#
# Usage: kaiju_superkingdom.sh nodes.dmp kaiju.out > summary.tsv
set -euo pipefail
nodes=$1
kaiju_out=$2

awk -v out="$(basename "$kaiju_out")" '
BEGIN { FS = "\t" }
FNR == NR { split($0, f, "\t[|]\t"); parent[f[1]] = f[2]; rank[f[1]] = f[3]; next }
{
    total++
    if ($1 != "C") { n["unclassified"]++; next }
    t = $3
    if (!(t in sk)) {
        x = t
        while (x != "" && x != "1" && rank[x] != "superkingdom") x = parent[x]
        sk[t] = (rank[x] == "superkingdom") ? x : "NA"
    }
    n[sk[t]]++
}
END {
    nm["2"] = "Bacteria"; nm["2157"] = "Archaea"; nm["2759"] = "Eukaryota"; nm["10239"] = "Viruses"
    nm["unclassified"] = "unclassified"; nm["NA"] = "cannot be assigned to a superkingdom"
    print "file\tpercent\treads\ttaxon_id\ttaxon_name"
    for (k in n) {
        id = (k == "unclassified" || k == "NA") ? "NA" : k
        name = (k in nm) ? nm[k] : k
        printf "%s\t%.6f\t%d\t%s\t%s\n", out, (total ? 100 * n[k] / total : 0), n[k], id, name
    }
}' "$nodes" "$kaiju_out"
