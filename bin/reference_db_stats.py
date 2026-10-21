#!/usr/bin/env python3
"""Summary tables for a euk-bot reference DB (the .faa.gz + .taxonomy.tsv.gz from build_reference_db.py),
for plotting elsewhere (e.g. Quarto/R). Reads both files once, in step: build_reference_db.py writes them in the
same order, and the script stops if the IDs ever disagree.

Outputs (TSV, in --outdir):
  size_by_source.tsv        source (MFT/EUK/PDB/all): proteins, amino acids, length mean/median/N50/min/max,
                            organisms, proteins without a lineage
  length_histogram.tsv      source, length, n_proteins (exact protein length distribution)
  organisms.tsv             one row per organism: lineage split into ranks, proteins, amino acids, annotation
                            coverage. Roll up to any rank from here for taxonomic coverage.
  annotation_by_source.tsv  source: proteins with a Pfam, distinct Pfam families (+ eggNOG/KEGG coverage)
  pfam_counts.tsv           source, pfam_id, pfam_name, n_proteins (Pfam shipped with MarFERReT)
With --eggnog (merged annotations from annotate_reference_db.sbatch), also:
  kegg_ko_counts.tsv        source, ko, n_proteins
  kegg_pathway_counts.tsv   source, pathway (mapNNNNN), n_proteins
  eggnog_pfam_counts.tsv    source, pfam_name, n_proteins (Pfam from eggNOG orthologs, for every source)
and eggNOG hit / KEGG KO / any-Pfam (MarFERReT or eggNOG) columns in organisms.tsv and annotation_by_source.tsv.

Rank columns: PR2 lineages (MarFERReT, and EukProt matched by genus) have 9 ranks, domain..species. PhyloDB
lineages have 7 (domain;phylum;class;order;family;genus;species); the phylum goes in `division` and
supergroup/subdivision stay blank. EukProt's UniEuk fallback lineages have no fixed ranks, so only `domain` is
filled; the full string is always in `lineage`.
"""
import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

from build_reference_db import open_text, read_fasta

RANKS = ["domain", "supergroup", "division", "subdivision", "class", "order", "family", "genus", "species"]
PHYLODB_RANKS = ["domain", "division", "class", "order", "family", "genus", "species"]
HAS_KO, HAS_PFAM = 1, 2


def split_lineage(lineage, lineage_source):
    parts = [p.strip() for p in lineage.split(";")] if lineage else []
    if lineage_source in ("PR2", "PR2_genus") and len(parts) == len(RANKS):
        return dict(zip(RANKS, parts))
    if lineage_source == "PhyloDB" and len(parts) == len(PHYLODB_RANKS):
        return dict(zip(PHYLODB_RANKS, parts))
    return {"domain": parts[0]} if parts else {}


def read_eggnog(path):
    """protein_id -> HAS_KO/HAS_PFAM flags for every protein with an eggNOG hit, plus per-source
    (source, ko), (source, pathway) and (source, pfam_name) counts. Source is the protein ID prefix."""
    flags, ko, pathway, pfam = {}, Counter(), Counter(), Counter()
    col = None
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                col = {name.lstrip("#"): i for i, name in enumerate(fields)}
                continue
            pid = fields[col["query"]]
            source = pid.split("_", 1)[0]
            flag = 0
            if fields[col["KEGG_ko"]] != "-":
                flag |= HAS_KO
                for k in fields[col["KEGG_ko"]].split(","):
                    ko[(source, k.removeprefix("ko:"))] += 1
                for p in fields[col["KEGG_Pathway"]].split(","):
                    if p.startswith("map"):   # each pathway is listed as both koNNNNN and mapNNNNN
                        pathway[(source, p)] += 1
            if fields[col["PFAMs"]] != "-":
                flag |= HAS_PFAM
                for name in fields[col["PFAMs"]].split(","):
                    pfam[(source, name)] += 1
            flags[pid] = flag
    print(f"{len(flags):,} proteins with an eggNOG hit", file=sys.stderr)
    return flags, ko, pathway, pfam


def length_stats(hist):
    """Mean, median, N50, min, max from a {length: count} histogram."""
    n = sum(hist.values())
    total = sum(length * count for length, count in hist.items())
    lengths = sorted(hist)
    seen, median = 0, None
    for length in lengths:
        seen += hist[length]
        if seen >= n / 2:
            median = length
            break
    covered, n50 = 0, None
    for length in reversed(lengths):
        covered += length * hist[length]
        if covered >= total / 2:
            n50 = length
            break
    return {"n_aa": total, "mean_length": round(total / n, 1), "median_length": median,
            "n50_length": n50, "min_length": lengths[0], "max_length": lengths[-1]}


def pct(part, whole):
    return round(100 * part / whole, 2) if whole else 0.0


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        out = csv.writer(fh, delimiter="\t", lineterminator="\n")
        out.writerow(header)
        out.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--faa", type=Path, required=True)
    ap.add_argument("--taxonomy", type=Path, required=True)
    ap.add_argument("--eggnog", type=Path, help="merged eggNOG-mapper annotations (.tsv or .tsv.gz)")
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    eggnog = read_eggnog(args.eggnog) if args.eggnog else None

    hist = defaultdict(Counter)       # source -> length -> proteins
    no_lineage = Counter()            # source -> proteins without a lineage
    organisms = defaultdict(Counter)  # (source, organism, taxid, lineage, lineage_source) -> counts
    pfam = Counter()                  # (source, pfam_id, pfam_name) -> proteins

    with open_text(args.taxonomy) as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        n = 0
        for (pid, seq), row in zip(read_fasta(args.faa), rows):
            if pid != row["protein_id"]:
                sys.exit(f"FASTA and taxonomy out of step at record {n + 1}: {pid} vs {row['protein_id']}")
            source, length = row["source"], len(seq.rstrip("*"))
            hist[source][length] += 1
            if not row["lineage"]:
                no_lineage[source] += 1
            org = organisms[(source, row["organism"], row["ncbi_taxid"], row["lineage"], row["lineage_source"])]
            org["n_proteins"] += 1
            org["n_aa"] += length
            if row["pfam_id"]:
                org["n_pfam"] += 1
                pfam[(source, row["pfam_id"], row["pfam_name"])] += 1
            if eggnog:
                flag = eggnog[0].get(pid)
                if flag is not None:
                    org["n_eggnog"] += 1
                    org["n_kegg_ko"] += bool(flag & HAS_KO)
                org["n_any_pfam"] += bool(row["pfam_id"] or (flag or 0) & HAS_PFAM)
            n += 1
            if n % 10_000_000 == 0:
                print(f"{n:,} proteins read", file=sys.stderr)
        if next(rows, None) is not None:
            sys.exit("taxonomy table has more rows than the FASTA")
    print(f"{n:,} proteins read", file=sys.stderr)

    sources = sorted(hist)
    hist["all"] = sum(hist.values(), Counter())
    no_lineage["all"] = sum(no_lineage.values())
    n_orgs = Counter(key[0] for key in organisms)
    n_orgs["all"] = len(organisms)
    by_source = defaultdict(Counter)
    for key, counts in organisms.items():
        by_source[key[0]].update(counts)

    size_cols = ["n_aa", "mean_length", "median_length", "n50_length", "min_length", "max_length"]
    write_tsv(args.outdir / "size_by_source.tsv",
              ["source", "n_proteins", *size_cols, "n_organisms", "n_proteins_no_lineage"],
              [[s, sum(hist[s].values()), *length_stats(hist[s]).values(), n_orgs[s], no_lineage[s]]
               for s in [*sources, "all"]])

    write_tsv(args.outdir / "length_histogram.tsv", ["source", "length", "n_proteins"],
              [[s, length, count] for s in sources for length, count in sorted(hist[s].items())])

    # count columns -> each is followed by its percent of n_proteins
    annot_cols = ["n_pfam", *(["n_eggnog", "n_kegg_ko", "n_any_pfam"] if eggnog else [])]

    def annot_values(c):
        return [v for col in annot_cols for v in (c[col], pct(c[col], c["n_proteins"]))]

    annot_header = [h for col in annot_cols for h in (col, "pct" + col[1:])]

    org_rows = []
    for (source, organism, taxid, lineage, lineage_source), c in sorted(organisms.items()):
        ranks = split_lineage(lineage, lineage_source)
        org_rows.append([source, organism, taxid, lineage_source, *(ranks.get(r, "") for r in RANKS), lineage,
                         c["n_proteins"], c["n_aa"], round(c["n_aa"] / c["n_proteins"], 1), *annot_values(c)])
    write_tsv(args.outdir / "organisms.tsv",
              ["source", "organism", "ncbi_taxid", "lineage_source", *RANKS, "lineage",
               "n_proteins", "n_aa", "mean_length", *annot_header], org_rows)

    n_distinct = defaultdict(Counter)   # source -> distinct Pfam IDs / KOs
    for source, *_ in pfam:
        n_distinct[source]["pfam"] += 1
    if eggnog:
        for source, _ko in eggnog[1]:
            n_distinct[source]["kegg_ko"] += 1
    distinct_cols = ["pfam", *(["kegg_ko"] if eggnog else [])]
    write_tsv(args.outdir / "annotation_by_source.tsv",
              ["source", "n_proteins", *annot_header, *(f"n_distinct_{d}" for d in distinct_cols)],
              [[s, by_source[s]["n_proteins"], *annot_values(by_source[s]), *(n_distinct[s][d] for d in distinct_cols)]
               for s in sources])

    def by_count(counts):
        return sorted(counts.items(), key=lambda kv: (kv[0][0], -kv[1], kv[0][1:]))

    write_tsv(args.outdir / "pfam_counts.tsv", ["source", "pfam_id", "pfam_name", "n_proteins"],
              [[*key, count] for key, count in by_count(pfam)])
    if eggnog:
        _flags, ko, pathway, eggnog_pfam = eggnog
        for name, header, counts in [("kegg_ko_counts.tsv", "ko", ko),
                                     ("kegg_pathway_counts.tsv", "pathway", pathway),
                                     ("eggnog_pfam_counts.tsv", "pfam_name", eggnog_pfam)]:
            write_tsv(args.outdir / name, ["source", header, "n_proteins"],
                      [[*key, count] for key, count in by_count(counts)])


if __name__ == "__main__":
    main()
