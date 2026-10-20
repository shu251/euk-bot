#!/usr/bin/env python3
"""Combine assembly QC results into one table, one row per assembly (<sample|coasm_group>.<assembler>).

Inputs: seqkit stats table (file column = <assembly>.contigs.fa), <assembly>.maprate.tsv files from
ASSEMBLY_MAPRATE, and optionally <assembly>.<lineage>.busco.txt BUSCO short summaries.
"""
import argparse
import csv
import re
from pathlib import Path

STATS_COLS = ["num_seqs", "sum_len", "min_len", "avg_len", "max_len", "N50", "GC(%)"]
BUSCO_RE = re.compile(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)")


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stats", type=Path, required=True)
    parser.add_argument("--maprate", type=Path, nargs="+", required=True)
    parser.add_argument("--busco", type=Path, nargs="*", default=[])
    parser.add_argument("--out", default="assembly_qc.tsv")
    args = parser.parse_args()

    rows = {}
    for rec in read_tsv(args.stats):
        name = Path(rec["file"]).name.removesuffix(".contigs.fa")
        rows[name] = {col: rec.get(col, "") for col in STATS_COLS}
    for path in args.maprate:
        for rec in read_tsv(path):
            rows.setdefault(rec["assembly"], {}).update(
                reads_processed=rec["reads_processed"], pct_reads_mapped=rec["pct_reads_mapped"])

    busco_cols = []
    for path in sorted(args.busco):
        name, lineage = path.name.removesuffix(".busco.txt").rsplit(".", 1)
        match = BUSCO_RE.search(path.read_text())
        cols = [f"busco_{lineage}_{k}" for k in ("complete", "single", "duplicated", "fragmented", "missing", "n")]
        busco_cols += [c for c in cols if c not in busco_cols]
        rows.setdefault(name, {}).update(zip(cols, match.groups() if match else [""] * len(cols)))

    header = ["assembly", "source", "assembler", *STATS_COLS, "reads_processed", "pct_reads_mapped", *busco_cols]
    with open(args.out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header)
        for name in sorted(rows):
            source, _, assembler = name.rpartition(".")
            writer.writerow([name, source, assembler, *(rows[name].get(col, "") for col in header[3:])])


if __name__ == "__main__":
    main()
