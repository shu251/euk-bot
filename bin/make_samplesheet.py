#!/usr/bin/env python3
"""Build a euk-bot samplesheet from a directory of Illumina paired-end FASTQ.

Expects bcl2fastq-style names: <sample>_S<n>_L<lane>_R{1,2}_001.fastq.gz.
Writes one row per lane; the pipeline concatenates lanes that share a sample name.
"""
import argparse
import csv
import re
import sys
from pathlib import Path

PATTERN = re.compile(r"^(?P<sample>.+)_S\d+_L(?P<lane>\d{3})_R1_001\.f(ast)?q\.gz$")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fastq_dir", type=Path)
    parser.add_argument("-o", "--output", default="samplesheet.csv")
    parser.add_argument("-n", "--max-samples", type=int, default=None,
                        help="keep only the first N samples (sorted by name), for test runs")
    args = parser.parse_args()

    fastq_dir = args.fastq_dir.resolve()
    rows = []
    problems = []
    for r1 in sorted(fastq_dir.iterdir()):
        match = PATTERN.match(r1.name)
        if not match:
            continue
        r2 = r1.with_name(r1.name.replace("_R1_001.", "_R2_001."))
        if not r2.exists():
            problems.append(f"missing R2 for {r1.name}")
            continue
        rows.append((match["sample"], str(r1), str(r2)))

    if problems:
        sys.exit("\n".join(problems))
    if not rows:
        sys.exit(f"No *_S*_L*_R1_001.fastq.gz files found in {fastq_dir}")

    samples = sorted({row[0] for row in rows})
    if args.max_samples:
        keep = set(samples[: args.max_samples])
        rows = [row for row in rows if row[0] in keep]
        samples = sorted(keep)

    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample", "fastq_1", "fastq_2"])
        writer.writerows(rows)

    print(f"Wrote {args.output}: {len(samples)} samples, {len(rows)} lane pairs", file=sys.stderr)


if __name__ == "__main__":
    main()
