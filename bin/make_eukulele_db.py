#!/usr/bin/env python3
"""EUKulele taxonomy files for a euk-bot reference DB version (tax-table.txt, prot-map.json).

EUKulele looks up each hit's organism in prot-map.json (protein ID -> Source_ID) and that organism's ranks in
tax-table.txt (Source_ID + Domain, Supergroup, Division, Class, Order, Family, Genus, Species: the ranks of its
default tax-cutoffs.yaml). euk-bot lineages map onto those ranks as:

  PR2 (MarFERReT; EukProt by PR2 genus): Domain;Supergroup;Division;Subdivision;Class;Order;Family;Genus[;Species]
      Subdivision is dropped (EUKulele has no such rank); a lineage without a species takes the organism name.
  PhyloDB (Bacteria, Archaea): Domain;Division;Class;Order;Family;Genus;Species
      Supergroup repeats Domain, as in EUKulele's own PhyloDB table.
  UniEuk (EukProt genera not in PR2): unranked, variable depth. Only Domain (Eukaryota) and, when the lineage
      has the organism's genus followed by a binomial ("Genus species"), Genus and Species are filled; the other
      ranks are left empty, which EUKulele reads as "No <rank>".
  No lineage (NA): Domain from the source (MarFERReT and EukProt are all eukaryotes), the rest empty.

Source_IDs are T000001.. (one per organism + lineage); tax-table.txt also lists organism, source and lineage_source.
EUKulele swaps "." for "N" in hit IDs before the prot-map lookup, so prot-map keys are written that way too.
"""
import argparse
import csv
import gzip
import json
import sys
from collections import Counter

RANKS = ["Domain", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species"]
EUKARYOTE_SOURCES = {"MFT", "EUK"}


def ranks_from_lineage(source, organism, lineage, lineage_source):
    """The 8 EUKulele ranks for one organism ('' = unknown), and which rule produced them."""
    fields = [f.strip() for f in lineage.split(";")] if lineage and lineage != "NA" else []
    empty = [""] * len(RANKS)

    if not fields:
        ranks = empty[:]
        ranks[0] = "Eukaryota" if source in EUKARYOTE_SOURCES else ""
        return ranks, "no_lineage"

    if lineage_source.startswith("PR2"):
        if len(fields) in (8, 9):
            kept = fields[:3] + fields[4:]                 # drop Subdivision
            if len(kept) == 7:
                kept.append(organism)
            return kept, "pr2"
    elif lineage_source == "PhyloDB" and len(fields) == 7:
        return [fields[0], fields[0]] + fields[1:], "phylodb"
    elif lineage_source == "UniEuk":
        ranks = empty[:]
        ranks[0] = fields[0]
        genus = organism.split("_")[0]
        if genus in fields and any(f.startswith(genus + " ") for f in fields):   # a binomial: genus is a genus
            ranks[6] = genus
            ranks[7] = organism
        return ranks, "unieuk"

    # Anything else: keep the domain only, so it is never placed below what is known
    ranks = empty[:]
    ranks[0] = fields[0]
    return ranks, f"other:{lineage_source}:{len(fields)}"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("taxonomy", help="euk-bot-ref.<version>.taxonomy.tsv.gz")
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args()

    source_ids = {}       # (organism, source, lineage, lineage_source) -> Source_ID
    tax_rows = []
    rules = Counter()
    n_proteins = 0

    with gzip.open(args.taxonomy, "rt", newline="") as fh, open(f"{args.outdir}/prot-map.json", "w") as out:
        reader = csv.DictReader(fh, delimiter="\t")
        out.write("{")
        for row in reader:
            key = (row["organism"], row["source"], row["lineage"], row["lineage_source"])
            sid = source_ids.get(key)
            if sid is None:
                sid = f"T{len(source_ids) + 1:06d}"
                source_ids[key] = sid
                ranks, rule = ranks_from_lineage(row["source"], row["organism"], row["lineage"], row["lineage_source"])
                rules[rule] += 1
                tax_rows.append([sid, *ranks, row["organism"], row["source"], row["lineage_source"]])
            out.write(("," if n_proteins else "") + json.dumps(row["protein_id"].replace(".", "N")) + ":" + json.dumps(sid))
            n_proteins += 1
        out.write("}\n")

    with open(f"{args.outdir}/tax-table.txt", "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
        writer.writerow(["Source_ID", *RANKS, "organism", "source", "lineage_source"])
        writer.writerows(tax_rows)

    print(f"proteins: {n_proteins}")
    print(f"organisms (Source_IDs): {len(source_ids)}")
    for rule, n in sorted(rules.items()):
        print(f"  {rule}: {n} organisms")
    if any(rule.startswith("other") for rule in rules):
        print("warning: some lineages matched no rule and were kept at Domain only", file=sys.stderr)


if __name__ == "__main__":
    main()
