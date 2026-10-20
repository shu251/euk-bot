#!/usr/bin/env python3
"""Write co-assembly groups chosen in read_profile.tsv into a samplesheet.

Takes the my_group column you filled in (blank = no co-assembly for that sample). Use
--fallback suggested or --fallback samplesheet to fill blanks from those columns instead.

Usage:
  python3 bin/apply_groups.py read_profile.tsv samplesheet.csv -o samplesheet_grouped.csv
"""
import argparse
import csv


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("read_profile")
    parser.add_argument("samplesheet")
    parser.add_argument("-o", "--out", required=True)
    parser.add_argument("--fallback", choices=["none", "suggested", "samplesheet"], default="none",
                        help="where to take a group from when my_group is blank")
    args = parser.parse_args()

    fallback_col = {"suggested": "suggested_group", "samplesheet": "samplesheet_group"}.get(args.fallback)
    with open(args.read_profile, newline="") as fh:
        chosen = {}
        for rec in csv.DictReader(fh, delimiter="\t"):
            group = rec.get("my_group", "").strip() or (rec.get(fallback_col, "").strip() if fallback_col else "")
            chosen[rec["sample"]] = group

    with open(args.samplesheet, newline="") as fh:
        rows = list(csv.DictReader(fh))
    fields = [f for f in rows[0].keys() if f != "group"] + ["group"]
    missing = sorted({r["sample"] for r in rows} - chosen.keys())
    if missing:
        print(f"warning: not in read profile, left ungrouped: {', '.join(missing)}")

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({**row, "group": chosen.get(row["sample"], "")})

    groups = sorted({g for g in chosen.values() if g})
    print(f"wrote {args.out}: {len(groups)} groups ({', '.join(groups)})")


if __name__ == "__main__":
    main()
