#!/usr/bin/env python3
"""Per-sample read accounting for the preliminary QC step (fastp + RiboDetector).

Every raw read pair ends up in exactly one of: removed by fastp, rRNA (RiboDetector), or non-rRNA (kept for
assembly and mapping). Their counts are given as pairs and as a percent of the sample's raw pairs, so
pct_removed_fastp + pct_rrna + pct_non_rrna = 100. The fastp detail columns split pct_removed_fastp by reason
(as a percent of raw reads). Without RiboDetector (--skip_rrna_removal) the rRNA columns are empty and every
fastp-passed pair counts as non-rRNA.

Output: one row per sample, sorted by sample.
"""
import argparse
import csv
import json
from pathlib import Path

COLUMNS = [
    "sample", "raw_pairs",
    "removed_fastp_pairs", "pct_removed_fastp", "rrna_pairs", "pct_rrna", "non_rrna_pairs", "pct_non_rrna",
    "pct_low_quality", "pct_too_many_n", "pct_too_short", "pct_adapter_trimmed", "pct_rrna_of_fastp_passed",
    "mean_read_length_raw", "mean_read_length_trimmed", "q30_pct_raw", "q30_pct_trimmed", "gc_pct_trimmed",
    "duplication_pct", "insert_size_peak",
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


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fastp", type=Path, nargs="+", required=True, help="<sample>.fastp.json files")
    parser.add_argument("--rrna", type=Path, nargs="*", default=[], help="<sample>.ribodetector.tsv files")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    rrna = {}
    for path in args.rrna:
        with open(path, newline="") as fh:
            for rec in csv.DictReader(fh, delimiter="\t"):
                rrna[rec["sample"]] = rec

    rows = []
    for path in args.fastp:
        sample = path.name.removesuffix(".fastp.json")
        f = fastp_row(path)
        raw, passed = f.pop("raw_pairs"), f.pop("passed_pairs")
        if sample in rrna:
            non_rrna = int(rrna[sample]["non_rrna_pairs"])
            n_rrna = int(rrna[sample]["input_pairs"]) - non_rrna
        else:
            non_rrna, n_rrna = passed, None
        rows.append({
            "sample": sample, "raw_pairs": raw,
            "removed_fastp_pairs": raw - passed, "pct_removed_fastp": pct(raw - passed, raw),
            "rrna_pairs": "" if n_rrna is None else n_rrna, "pct_rrna": "" if n_rrna is None else pct(n_rrna, raw),
            "non_rrna_pairs": non_rrna, "pct_non_rrna": pct(non_rrna, raw),
            "pct_rrna_of_fastp_passed": "" if n_rrna is None else pct(n_rrna, passed),
            **f,
        })

    with open(args.out, "w", newline="") as fh:
        out = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        out.writeheader()
        out.writerows(sorted(rows, key=lambda r: r["sample"]))


if __name__ == "__main__":
    main()
