#!/usr/bin/env python3
"""Per-sample read profile for judging read quality and choosing co-assembly groups.

Combines fastp reports, RiboDetector rRNA fractions, Nonpareil coverage/diversity, strandedness and
Kaiju superkingdom composition with sourmash read similarity. Suggests groups by average-linkage
clustering of samples on sourmash similarity, cut at --group-similarity.

Outputs:
  read_profile.tsv            one row per sample; the my_group / my_notes columns are left empty to
                              fill in by hand, then apply_groups.py writes them into a samplesheet
  samplesheet_suggested.csv   the input samplesheet with group = suggested_group
"""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

KAIJU_TAXA = {"Eukaryota": "euk", "Bacteria": "bac", "Archaea": "arc", "Viruses": "vir", "unclassified": "unclassified"}


def read_tsvs(paths, key="sample"):
    rows = {}
    for path in paths:
        with open(path, newline="") as fh:
            for rec in csv.DictReader(fh, delimiter="\t"):
                rows[rec[key]] = rec
    return rows


def fastp_stats(paths):
    out = {}
    for path in paths:
        d = json.loads(Path(path).read_text())
        before, after = d["summary"]["before_filtering"], d["summary"]["after_filtering"]
        out[Path(path).name.removesuffix(".fastp.json")] = {
            "reads_raw": before["total_reads"],
            "reads_after_fastp": after["total_reads"],
            "pct_kept_fastp": round(100 * after["total_reads"] / before["total_reads"], 2) if before["total_reads"] else "",
            "duplication_pct": round(100 * d.get("duplication", {}).get("rate", 0), 2),
            "insert_size_peak": d.get("insert_size", {}).get("peak", ""),
            "q30_pct": round(100 * after["q30_rate"], 2),
            "gc_pct": round(100 * after["gc_content"], 2),
        }
    return out


def kaiju_stats(paths):
    out = defaultdict(dict)
    for path in paths:
        sample = Path(path).name.removesuffix(".kaiju_superkingdom.tsv")
        with open(path, newline="") as fh:
            for rec in csv.DictReader(fh, delimiter="\t"):
                short = KAIJU_TAXA.get(rec["taxon_name"])
                if short:
                    out[sample][f"kaiju_pct_{short}"] = rec["percent"]
    return out


def read_similarity(path):
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh))
    if not rows:
        return [], []
    names = [Path(n).name.removesuffix(".sig") for n in rows[0]]
    return names, [[float(x) for x in row] for row in rows[1:]]


def average_linkage_groups(names, sim, threshold):
    """Merge the most similar clusters (mean pairwise similarity) until none are >= threshold."""
    clusters = [[i] for i in range(len(names))]

    def link(a, b):
        return sum(sim[i][j] for i in a for j in b) / (len(a) * len(b))

    while len(clusters) > 1:
        best, pair = -1.0, None
        for x in range(len(clusters)):
            for y in range(x + 1, len(clusters)):
                s = link(clusters[x], clusters[y])
                if s > best:
                    best, pair = s, (x, y)
        if best < threshold:
            break
        x, y = pair
        clusters[x] += clusters.pop(y)

    clusters.sort(key=lambda c: (-len(c), min(names[i] for i in c)))
    return {names[i]: f"sim_group{n}" for n, c in enumerate(clusters, 1) for i in c}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--samplesheet", type=Path, required=True)
    parser.add_argument("--fastp", type=Path, nargs="+", required=True)
    parser.add_argument("--rrna", type=Path, nargs="*", default=[])
    parser.add_argument("--similarity", type=Path, required=True)
    parser.add_argument("--nonpareil", type=Path, nargs="*", default=[])
    parser.add_argument("--strandedness", type=Path, nargs="*", default=[])
    parser.add_argument("--kaiju", type=Path, nargs="*", default=[])
    parser.add_argument("--group-similarity", type=float, default=0.5)
    args = parser.parse_args()

    with open(args.samplesheet, newline="") as fh:
        sheet = list(csv.DictReader(fh))
    sheet_group = {row["sample"]: row.get("group", "") or "" for row in sheet}
    samples = sorted(sheet_group)

    fastp = fastp_stats(args.fastp)
    rrna = read_tsvs(args.rrna)
    nonpareil = read_tsvs(args.nonpareil)
    strand = read_tsvs(args.strandedness)
    kaiju = kaiju_stats(args.kaiju)
    names, sim = read_similarity(args.similarity)
    idx = {n: i for i, n in enumerate(names)}
    suggested = average_linkage_groups(names, sim, args.group_similarity) if names else {}

    def nearest(s):
        if s not in idx or len(names) < 2:
            return "", ""
        i = idx[s]
        j = max((j for j in range(len(names)) if j != i), key=lambda j: sim[i][j])
        return names[j], round(sim[i][j], 3)

    def mean_within_group(s):
        g = sheet_group.get(s)
        mates = [t for t in samples if t != s and g and sheet_group[t] == g and t in idx]
        if s not in idx or not mates:
            return ""
        return round(sum(sim[idx[s]][idx[t]] for t in mates) / len(mates), 3)

    kaiju_cols = [f"kaiju_pct_{v}" for v in KAIJU_TAXA.values()] if args.kaiju else []
    np_cols = ["np_coverage", "np_diversity", "np_effort_for_95pct_bp"] if args.nonpareil else []
    strand_cols = ["strandedness", "frac_r1_sense"] if args.strandedness else []
    header = [
        "sample", "samplesheet_group", "suggested_group", "nearest_sample", "nearest_similarity",
        "mean_similarity_within_samplesheet_group",
        "reads_raw", "reads_after_fastp", "pct_kept_fastp", "duplication_pct", "insert_size_peak", "q30_pct", "gc_pct",
        "pct_rrna", "non_rrna_pairs", *np_cols, *strand_cols, *kaiju_cols,
        "my_group", "my_notes",
    ]
    with open("read_profile.tsv", "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header)
        for s in samples:
            near, near_sim = nearest(s)
            row = {
                "sample": s, "samplesheet_group": sheet_group[s], "suggested_group": suggested.get(s, ""),
                "nearest_sample": near, "nearest_similarity": near_sim,
                "mean_similarity_within_samplesheet_group": mean_within_group(s),
                **fastp.get(s, {}),
                "pct_rrna": rrna.get(s, {}).get("pct_rrna", ""),
                "non_rrna_pairs": rrna.get(s, {}).get("non_rrna_pairs", ""),
                **{c: nonpareil.get(s, {}).get(c, "") for c in np_cols},
                **{c: strand.get(s, {}).get(c, "") for c in strand_cols},
                **kaiju.get(s, {}),
            }
            writer.writerow([row.get(c, "") for c in header])

    with open("samplesheet_suggested.csv", "w", newline="") as fh:
        fields = [f for f in sheet[0].keys() if f != "group"] + ["group"] if sheet else ["sample", "group"]
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in sheet:
            writer.writerow({**row, "group": suggested.get(row["sample"], "")})


if __name__ == "__main__":
    main()
