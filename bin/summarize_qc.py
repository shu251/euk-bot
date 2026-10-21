#!/usr/bin/env python3
"""Per-sample read accounting for the preliminary QC step (fastp -> BBDuk -> RiboDetector, then BBNorm).

Every raw read pair ends up in exactly one of: removed by fastp, removed by BBDuk (PhiX / artifacts), rRNA
(RiboDetector), or non-rRNA (the clean reads used by every later step). Counts are given as pairs and as a
percent of the sample's raw pairs, so pct_removed_fastp + pct_removed_bbduk + pct_rrna + pct_non_rrna = 100.
A step that was skipped has empty columns and removes nothing. BBNorm only thins the assembler input, so
normalized_pairs / pct_normalized_of_non_rrna sit outside that sum.

The fastp detail columns split pct_removed_fastp by reason (as a percent of raw reads). qc_flag marks samples
that lost more than --max-pct-removed-fastp of their pairs in fastp.

Output: one row per sample, sorted by sample.
"""
import argparse
import csv
import json
from pathlib import Path

COLUMNS = [
    "sample", "raw_pairs",
    "removed_fastp_pairs", "pct_removed_fastp", "removed_bbduk_pairs", "pct_removed_bbduk",
    "rrna_pairs", "pct_rrna", "non_rrna_pairs", "pct_non_rrna",
    "normalized_pairs", "pct_normalized_of_non_rrna",
    "pct_low_quality", "pct_too_many_n", "pct_too_short", "pct_adapter_trimmed", "pct_rrna_of_rrna_input",
    "mean_read_length_raw", "mean_read_length_trimmed", "q30_pct_raw", "q30_pct_trimmed", "gc_pct_trimmed",
    "duplication_pct", "insert_size_peak", "qc_flag",
]


def pct(part, whole):
    return round(100 * part / whole, 2) if whole else ""


def fastp_row(path):
    d = json.loads(Path(path).read_text())
    before, after = d["summary"]["before_filtering"], d["summary"]["after_filtering"]
    filt = d.get("filtering_result", {})
    raw_reads = before["total_reads"]
    return {
        "raw_pairs": raw_reads // 2,
        "passed_pairs": after["total_reads"] // 2,
        "pct_low_quality": pct(filt.get("low_quality_reads", 0), raw_reads),
        "pct_too_many_n": pct(filt.get("too_many_N_reads", 0), raw_reads),
        "pct_too_short": pct(filt.get("too_short_reads", 0), raw_reads),
        "pct_adapter_trimmed": pct(d.get("adapter_cutting", {}).get("adapter_trimmed_reads", 0), raw_reads),
        "mean_read_length_raw": before["read1_mean_length"],
        "mean_read_length_trimmed": after["read1_mean_length"],
        "q30_pct_raw": round(100 * before["q30_rate"], 2),
        "q30_pct_trimmed": round(100 * after["q30_rate"], 2),
        "gc_pct_trimmed": round(100 * after["gc_content"], 2),
        "duplication_pct": round(100 * d.get("duplication", {}).get("rate", 0), 2),
        "insert_size_peak": d.get("insert_size", {}).get("peak", ""),
    }


def read_tsvs(paths):
    rows = {}
    for path in paths:
        with open(path, newline="") as fh:
            for rec in csv.DictReader(fh, delimiter="\t"):
                rows[rec["sample"]] = rec
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fastp", type=Path, nargs="+", required=True, help="<sample>.fastp.json files")
    parser.add_argument("--bbduk", type=Path, nargs="*", default=[], help="<sample>.bbduk.tsv files")
    parser.add_argument("--rrna", type=Path, nargs="*", default=[], help="<sample>.ribodetector.tsv files")
    parser.add_argument("--bbnorm", type=Path, nargs="*", default=[], help="<sample>.bbnorm.tsv files")
    parser.add_argument("--max-pct-removed-fastp", type=float, default=30.0)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    bbduk, rrna, bbnorm = read_tsvs(args.bbduk), read_tsvs(args.rrna), read_tsvs(args.bbnorm)

    rows = []
    for path in args.fastp:
        sample = path.name.removesuffix(".fastp.json")
        f = fastp_row(path)
        raw, after_fastp = f.pop("raw_pairs"), f.pop("passed_pairs")
        after_bbduk = int(bbduk[sample]["kept_pairs"]) if sample in bbduk else after_fastp
        if sample in rrna:
            rrna_input, non_rrna = int(rrna[sample]["input_pairs"]), int(rrna[sample]["non_rrna_pairs"])
            n_rrna = rrna_input - non_rrna
        else:
            rrna_input, non_rrna, n_rrna = None, after_bbduk, None
        n_norm = int(bbnorm[sample]["normalized_pairs"]) if sample in bbnorm else None
        removed_fastp = raw - after_fastp
        pct_removed = pct(removed_fastp, raw)

        def blank_if_none(value, fn=lambda v: v):
            return "" if value is None else fn(value)

        rows.append({
            "sample": sample, "raw_pairs": raw,
            "removed_fastp_pairs": removed_fastp, "pct_removed_fastp": pct_removed,
            "removed_bbduk_pairs": after_fastp - after_bbduk if sample in bbduk else "",
            "pct_removed_bbduk": pct(after_fastp - after_bbduk, raw) if sample in bbduk else "",
            "rrna_pairs": blank_if_none(n_rrna), "pct_rrna": blank_if_none(n_rrna, lambda n: pct(n, raw)),
            "non_rrna_pairs": non_rrna, "pct_non_rrna": pct(non_rrna, raw),
            "normalized_pairs": blank_if_none(n_norm),
            "pct_normalized_of_non_rrna": blank_if_none(n_norm, lambda n: pct(n, non_rrna)),
            "pct_rrna_of_rrna_input": blank_if_none(n_rrna, lambda n: pct(n, rrna_input)),
            "qc_flag": (f"removed_fastp>{args.max_pct_removed_fastp:g}%"
                        if pct_removed != "" and pct_removed > args.max_pct_removed_fastp else ""),
            **f,
        })

    with open(args.out, "w", newline="") as fh:
        out = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        out.writeheader()
        out.writerows(sorted(rows, key=lambda r: r["sample"]))


if __name__ == "__main__":
    main()
