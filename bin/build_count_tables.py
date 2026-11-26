#!/usr/bin/env python3
"""Count tables from the contig table: per-contig counts with taxonomy + KO, and counts summed by taxon and KO.

Taxonomy is the contig-level EUKulele consensus (contig_eukulele_full_classification), split into rank columns
(domain, supergroup, division, class, order, family, genus, species); ranks below the consensus depth are empty.
Function is the KEGG KO of the representative protein's own eggNOG-mapper annotation (eggnog_KEGG_ko). A contig
with several KOs has its reads split evenly between them; one without a KO counts as KO "unassigned".

Library sizes per sample (from each salmon quant dir):
  mapped   sum of salmon NumReads over all contigs (CPM columns then sum to 1e6 per sample)
  input    read pairs salmon was given (aux_info/meta_info.json num_processed): the clean reads, mapped or not
CPM = reads / library size * 1e6, written for each library size in --libsize (none: raw reads only).

Writes (with --prefix):
  library_sizes.tsv          sample, input_reads, mapped_reads, percent_mapped
  contig_counts.tsv          one row per contig: taxonomy, KO, eggNOG name/description/COG, reads per sample,
                             then <sample>_cpm_<libsize> per sample
  ko_counts.tsv              long: sample, KO, reads, cpm_<libsize>...
  taxon_counts.tsv           long: sample, all rank columns, reads, cpm...        (finest consensus lineage)
  taxon_ko_counts.tsv        long: sample, all rank columns, KO, reads, cpm...
  taxon_counts.<rank>.tsv    as above, lineages cut at --rank (contigs resolved above it keep their shorter lineage)
  taxon_ko_counts.<rank>.tsv
Long tables only have rows with reads > 0. Because counts are sums of raw reads, the all-rank tables can be summed
to any coarser group; CPM sums the same way.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

RANKS = ["domain", "supergroup", "division", "class", "order", "family", "genus", "species"]
LIBSIZES = ["mapped", "input"]
TAX_COL = "contig_eukulele_full_classification"
LEVEL_COL = "contig_eukulele_classification_level"
SUPPORT_COL = "contig_eukulele_support"
KO_COL = "eggnog_KEGG_ko"
FUNC_COLS = ["eggnog_Preferred_name", "eggnog_Description", "eggnog_COG_category"]
UNASSIGNED = "unassigned"

csv.field_size_limit(sys.maxsize)


def lineage(value):
    """Consensus lineage -> 8 rank values (empty below its depth)."""
    names = [n for n in value.split(";")] if value else []
    names = names[:len(RANKS)]
    return tuple(names + [""] * (len(RANKS) - len(names)))


def kos(value):
    """'ko:K00001,ko:K00002' -> ('K00001', 'K00002'); none -> ('unassigned',)."""
    out = []
    for item in value.split(","):
        item = item.strip()
        if item.startswith("ko:"):
            item = item[3:]
        if item and item != "-" and item not in out:
            out.append(item)
    return tuple(out) or (UNASSIGNED,)


def read_library_sizes(quant_dirs):
    """{sample: num_processed} from salmon aux_info/meta_info.json, sample = dir name."""
    sizes = {}
    for d in quant_dirs:
        meta = d / "aux_info" / "meta_info.json"
        if not meta.exists():
            sys.exit(f"build_count_tables.py: {meta} not found")
        with open(meta) as fh:
            sizes[d.name] = int(json.load(fh)["num_processed"])
    return sizes


def fmt(x):
    return f"{x:.6g}"


def add(table, key, values):
    acc = table.get(key)
    if acc is None:
        table[key] = list(values)
    else:
        for i, v in enumerate(values):
            acc[i] += v


def write_long(path, header_keys, table, samples, scale):
    """Long table: sample, *key columns, reads, cpm per library size; rows with reads > 0 only."""
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample", *header_keys, "reads", *(f"cpm_{lib}" for lib in scale)])
        for i, sample in enumerate(samples):
            for key in sorted(table):
                reads = table[key][i]
                if reads > 0:
                    w.writerow([sample, *key, fmt(reads), *(fmt(reads * f[i]) for f in scale.values())])


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--contig-table", required=True, type=Path, help="contig_table.tsv (build_contig_table.py)")
    parser.add_argument("--quant-dirs", required=True, nargs="+", type=Path, help="salmon quant dirs, named by sample")
    parser.add_argument("--rank", default="class", choices=RANKS, help="rank for the cut taxon tables (default: class)")
    parser.add_argument("--libsize", default="mapped,input",
                        help="comma-separated library sizes for CPM: mapped, input, or none (default: mapped,input)")
    parser.add_argument("--prefix", default="", help="output file name prefix")
    args = parser.parse_args()

    libsizes = [] if args.libsize.strip() == "none" else [x.strip() for x in args.libsize.split(",") if x.strip()]
    bad = set(libsizes) - set(LIBSIZES)
    if bad:
        sys.exit(f"build_count_tables.py: --libsize must be mapped, input or none (got {', '.join(sorted(bad))})")
    cut = RANKS.index(args.rank) + 1

    input_reads = read_library_sizes(args.quant_dirs)
    samples = sorted(input_reads)
    n = len(samples)

    def contig_rows():
        """Yield (contig, length, level, support, lineage, kos, func, reads) for each contig-table row."""
        with open(args.contig_table, newline="") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader)
            col = {c: i for i, c in enumerate(header)}
            missing = [s for s in samples if s not in col]
            if missing:
                sys.exit(f"build_count_tables.py: samples missing from {args.contig_table}: {', '.join(missing)}")
            sample_idx = [col[s] for s in samples]
            idx = {c: col.get(c) for c in (TAX_COL, LEVEL_COL, SUPPORT_COL, KO_COL, *FUNC_COLS)}
            get = lambda row, c: row[idx[c]] if idx[c] is not None else ""
            for row in reader:
                yield (row[col["contig"]], row[col["length"]], get(row, LEVEL_COL), get(row, SUPPORT_COL),
                       lineage(get(row, TAX_COL)), kos(get(row, KO_COL)), [get(row, c) for c in FUNC_COLS],
                       [float(row[i]) for i in sample_idx])

    with open(args.contig_table, newline="") as fh:
        header = set(next(csv.reader(fh, delimiter="\t")))
    for c in (TAX_COL, KO_COL):
        if c not in header:
            print(f"warning: no {c} column in {args.contig_table}; treated as unassigned", file=sys.stderr)

    # Pass 1: totals and sums (pass 2 below writes the per-contig table, once the CPM factors are known)
    n_contigs = 0
    mapped = [0.0] * n
    by_ko, by_tax, by_tax_ko, by_cut, by_cut_ko = {}, {}, {}, {}, {}
    for _name, _length, _level, _support, lin, ko, _func, reads in contig_rows():
        n_contigs += 1
        if not any(reads):
            continue
        for i, r in enumerate(reads):
            mapped[i] += r
        lin_cut = lin[:cut]
        add(by_tax, lin, reads)
        add(by_cut, lin_cut, reads)
        share = [r / len(ko) for r in reads]
        for k in ko:
            add(by_ko, (k,), share)
            add(by_tax_ko, (*lin, k), share)
            add(by_cut_ko, (*lin_cut, k), share)

    totals = {"mapped": mapped, "input": [float(input_reads[s]) for s in samples]}
    scale = {lib: [1e6 / t if t else 0.0 for t in totals[lib]] for lib in libsizes}

    with open(f"{args.prefix}library_sizes.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample", "input_reads", "mapped_reads", "percent_mapped"])
        for i, s in enumerate(samples):
            pct = 100 * mapped[i] / input_reads[s] if input_reads[s] else 0.0
            w.writerow([s, input_reads[s], fmt(mapped[i]), f"{pct:.2f}"])

    with open(f"{args.prefix}contig_counts.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["contig", "length", "taxonomy_level", "taxonomy_support", *RANKS, "KO", "n_KO",
                    *(c.replace("eggnog_", "") for c in FUNC_COLS), *samples,
                    *(f"{s}_cpm_{lib}" for lib in scale for s in samples)])
        for name, length, level, support, lin, ko, func, reads in contig_rows():
            n_ko = 0 if ko == (UNASSIGNED,) else len(ko)
            w.writerow([name, length, level or "none", support, *lin, ",".join(ko), n_ko, *func,
                        *(fmt(r) for r in reads),
                        *(fmt(r * f[i]) for f in scale.values() for i, r in enumerate(reads))])

    cut_ranks = RANKS[:cut]
    write_long(f"{args.prefix}ko_counts.tsv", ["KO"], by_ko, samples, scale)
    write_long(f"{args.prefix}taxon_counts.tsv", RANKS, by_tax, samples, scale)
    write_long(f"{args.prefix}taxon_ko_counts.tsv", [*RANKS, "KO"], by_tax_ko, samples, scale)
    write_long(f"{args.prefix}taxon_counts.{args.rank}.tsv", cut_ranks, by_cut, samples, scale)
    write_long(f"{args.prefix}taxon_ko_counts.{args.rank}.tsv", [*cut_ranks, "KO"], by_cut_ko, samples, scale)

    print(f"contigs: {n_contigs}; samples: {n}; KOs: {len(by_ko)}; lineages: {len(by_tax)}; "
          f"{args.rank}-level lineages: {len(by_cut)}")
    for i, s in enumerate(samples):
        print(f"{s}: input {input_reads[s]}, mapped {mapped[i]:.0f}")


if __name__ == "__main__":
    main()
