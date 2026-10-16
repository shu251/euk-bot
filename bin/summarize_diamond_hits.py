#!/usr/bin/env python3
"""Count best DIAMOND blastx hits per reference protein per sample.

Inputs are <sample>.blastx.tsv files (outfmt 6: qseqid sseqid pident length evalue bitscore).
Mates of a pair share a read ID, so each pair is counted once, using its highest-bitscore hit.
Optionally joins a taxonomy/annotation TSV (with header; first column = reference protein ID).

Output: long table  sample, sseqid, n_reads[, taxonomy columns...]
"""
import argparse
import csv
import re
from collections import Counter
from pathlib import Path

MATE_SUFFIX = re.compile(r"/[12]$")


def best_hit_counts(path):
    best = {}
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            read_id = MATE_SUFFIX.sub("", fields[0])
            sseqid, bitscore = fields[1], float(fields[5])
            if read_id not in best or bitscore > best[read_id][1]:
                best[read_id] = (sseqid, bitscore)
    return Counter(sseqid for sseqid, _ in best.values())


def read_taxonomy(path):
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        return header[1:], {row[0]: row[1:] for row in reader}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("hit_tables", nargs="+", type=Path)
    parser.add_argument("--taxonomy", type=Path)
    parser.add_argument("--out", default="reference_hit_counts.tsv")
    args = parser.parse_args()

    tax_cols, taxonomy = read_taxonomy(args.taxonomy) if args.taxonomy else ([], {})
    missing = [""] * len(tax_cols)

    with open(args.out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["sample", "sseqid", "n_reads", *tax_cols])
        for path in sorted(args.hit_tables):
            sample = path.name.removesuffix(".blastx.tsv")
            for sseqid, n in sorted(best_hit_counts(path).items()):
                writer.writerow([sample, sseqid, n, *taxonomy.get(sseqid, missing)])


if __name__ == "__main__":
    main()
