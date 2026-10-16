#!/usr/bin/env python3
"""Merge per-sample salmon quant.sf files into wide NumReads and TPM tables.

Each input is a salmon output directory named after its sample.
Writes <prefix>_numreads.tsv and <prefix>_tpm.tsv (rows = contigs, columns = samples).
"""
import argparse
import csv
from pathlib import Path

COLUMNS = ("NumReads", "TPM")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("quant_dirs", nargs="+", type=Path)
    parser.add_argument("--prefix", default="salmon")
    args = parser.parse_args()

    samples = []
    lengths = {}
    tables = {col: {} for col in COLUMNS}
    for quant_dir in sorted(args.quant_dirs):
        sample = quant_dir.name
        samples.append(sample)
        with open(quant_dir / "quant.sf") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                name = row["Name"]
                lengths[name] = row["Length"]
                for col in COLUMNS:
                    tables[col].setdefault(name, {})[sample] = row[col]

    for col, table in tables.items():
        with open(f"{args.prefix}_{col.lower()}.tsv", "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["Name", "Length", *samples])
            for name in sorted(table):
                writer.writerow([name, lengths[name], *(table[name].get(s, "0") for s in samples)])


if __name__ == "__main__":
    main()
