#!/usr/bin/env python3
"""Contig-level taxonomy + function + counts table for the assembly branch.

Joins, per predicted protein (TransDecoder and/or MetaEuk):
  best DIAMOND hit against the reference (highest bitscore)    ref_hit, ref_pident, ref_evalue, ref_bitscore
  the hit's reference taxonomy (--taxonomy)                    ref_<column>  (source, organism, lineage, ...)
  the hit's eggNOG annotation (--reference-annotations)        ref_eggnog_<column>
  EUKulele consensus taxonomy of the protein                   eukulele_<column>
  eggNOG-mapper annotation of the protein itself               eggnog_<column>
and writes
  <prefix>protein_table.tsv  one row per protein
  <prefix>contig_table.tsv   one row per clustered contig: its representative protein's columns, n_proteins,
                             then the salmon NumReads of every sample
The representative protein of a contig is the one with the best reference hit (bitscore), then one with an
eggNOG annotation, then the longest. Contigs without a predicted protein are kept, with empty annotation.

Inputs are recognised by file name:
  *.transdecoder.pep, *.metaeuk.fas      proteins (ORF caller transdecoder / metaeuk)
  <caller>.proteins_vs_reference.tsv     DIAMOND blastp, outfmt 6 (12 standard columns)
  <caller>.emapper.annotations           eggNOG-mapper
  <caller>.eukulele_taxonomy.tsv         EUKulele taxonomy estimate
Reference tables are read for the hit proteins only, so they can be much larger than memory.
"""
import argparse
import csv
import gzip
import re
import sys
from pathlib import Path

EGGNOG_COLS = ["Preferred_name", "Description", "COG_category", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
               "CAZy", "PFAMs", "max_annot_lvl"]
EUKULELE_COLS = ["classification_level", "classification", "full_classification", "ambiguous"]
HIT_COLS = ["ref_hit", "ref_pident", "ref_evalue", "ref_bitscore"]
TRANSDECODER_SUFFIX = re.compile(r"\.p\d+$")

csv.field_size_limit(sys.maxsize)


def open_text(path):
    return gzip.open(path, "rt", newline="") if str(path).endswith(".gz") else open(path, newline="")


def classify_inputs(paths):
    """{kind: {caller: path}} from file names."""
    kinds = {"proteins": {}, "hits": {}, "eggnog": {}, "eukulele": {}}
    rules = [
        ("proteins", re.compile(r".*\.(transdecoder)\.pep$")),
        ("proteins", re.compile(r".*\.(metaeuk)\.fas$")),
        ("hits", re.compile(r"(.+)\.proteins_vs_reference\.tsv$")),
        ("eggnog", re.compile(r"(.+)\.emapper\.annotations$")),
        ("eukulele", re.compile(r"(.+)\.eukulele_taxonomy\.tsv$")),
    ]
    for path in paths:
        for kind, rule in rules:
            match = rule.match(path.name)
            if match:
                kinds[kind][match.group(1)] = path
                break
        else:
            sys.exit(f"build_contig_table.py: don't know what {path} is (see --help for the file names)")
    return kinds


def contig_of(caller, protein_id):
    """TransDecoder: <contig>.p<n>; MetaEuk: <target>|<contig>|<strand>|..."""
    if caller == "metaeuk" or "|" in protein_id:
        return protein_id.split("|")[1]
    return TRANSDECODER_SUFFIX.sub("", protein_id)


def read_proteins(path, caller):
    """[(protein_id, contig, length)] in file order."""
    proteins, pid, length = [], None, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if pid is not None:
                    proteins.append((pid, contig_of(caller, pid), length))
                pid, length = line[1:].split(None, 1)[0], 0
            else:
                length += len(line.strip().rstrip("*"))
    if pid is not None:
        proteins.append((pid, contig_of(caller, pid), length))
    return proteins


def best_hits(path):
    """qseqid -> [sseqid, pident, evalue, bitscore] of its highest-bitscore hit (the first one on ties)."""
    best = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            q, bitscore = f[0], float(f[11])
            if q not in best or bitscore > best[q][1]:
                best[q] = ([f[1], f[2], f[10], f[11]], bitscore)
    return {q: hit for q, (hit, _) in best.items()}


def read_keyed_table(path, keep_ids, columns=None, comment="##"):
    """Rows keyed by first column, only for keep_ids (None = all). Header may start with '#'.
    Returns (column names, {id: [values]}); with `columns`, only those that exist, in that order."""
    with open_text(path) as fh:
        lines = (line for line in fh if not (comment and line.startswith(comment)))
        reader = csv.reader(lines, delimiter="\t")
        header = next(reader, None)
        if header is None:
            return (columns or []), {}
        header[0] = header[0].lstrip("#")
        idx = [header.index(c) for c in columns if c in header] if columns else list(range(1, len(header)))
        names = [header[i] for i in idx]
        rows = {}
        for row in reader:
            if row and (keep_ids is None or row[0] in keep_ids):
                rows[row[0]] = [row[i] if i < len(row) else "" for i in idx]
    return names, rows


def read_eukulele(path):
    """transcript_name -> EUKULELE_COLS values (EUKulele writes a pandas index as the first column)."""
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return {row["transcript_name"]: [row.get(c, "") for c in EUKULELE_COLS] for row in reader}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("inputs", nargs="+", type=Path, help="protein, hit, eggNOG and EUKulele files (see above)")
    parser.add_argument("--numreads", required=True, type=Path, help="salmon_numreads.tsv (Name, Length, samples...)")
    parser.add_argument("--taxonomy", type=Path, help="reference taxonomy TSV(.gz), first column = protein ID")
    parser.add_argument("--reference-annotations", type=Path,
                        help="eggNOG-mapper annotations of the reference proteins (.tsv.gz)")
    parser.add_argument("--prefix", default="", help="output file name prefix")
    args = parser.parse_args()

    inputs = classify_inputs(args.inputs)
    callers = sorted(inputs["proteins"])
    for kind in ("hits", "eggnog", "eukulele"):
        extra = set(inputs[kind]) - set(callers)
        if extra:
            sys.exit(f"build_contig_table.py: {kind} for ORF caller(s) without proteins: {', '.join(sorted(extra))}")

    # Contigs and counts
    with open(args.numreads, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        samples = header[2:]
        contigs = {row[0]: row for row in reader}

    # Per protein: caller, id, contig, length, best hit, then the per-caller annotations
    proteins = []
    hits, eggnog, eukulele = {}, {}, {}
    eggnog_names = []
    for caller in callers:
        for pid, contig, length in read_proteins(inputs["proteins"][caller], caller):
            proteins.append((caller, pid, contig, length))
        if caller in inputs["hits"]:
            hits[caller] = best_hits(inputs["hits"][caller])
        if caller in inputs["eggnog"]:
            names, eggnog[caller] = read_keyed_table(inputs["eggnog"][caller], None, EGGNOG_COLS)
            eggnog_names = eggnog_names or names
        if caller in inputs["eukulele"]:
            eukulele[caller] = read_eukulele(inputs["eukulele"][caller])
    if eggnog and not eggnog_names:
        eggnog_names = EGGNOG_COLS

    unknown = sum(1 for _, _, contig, _ in proteins if contig not in contigs)
    if unknown:
        print(f"warning: {unknown} proteins are on contigs missing from {args.numreads}", file=sys.stderr)

    hit_ids = {hit[0] for caller_hits in hits.values() for hit in caller_hits.values()}
    tax_names, taxonomy = read_keyed_table(args.taxonomy, hit_ids, comment=None) if args.taxonomy else ([], {})
    ref_names, ref_annot = (read_keyed_table(args.reference_annotations, hit_ids, EGGNOG_COLS)
                            if args.reference_annotations else ([], {}))

    columns = (["orf_caller", "protein_id", "contig", "protein_length"]
               + (HIT_COLS if hits else [])
               + ["ref_" + c for c in tax_names]
               + ["ref_eggnog_" + c for c in ref_names]
               + (["eukulele_" + c for c in EUKULELE_COLS] if eukulele else [])
               + ["eggnog_" + c for c in eggnog_names])
    empty = {"hit": [""] * len(HIT_COLS), "tax": [""] * len(tax_names), "ref": [""] * len(ref_names),
             "euk": [""] * len(EUKULELE_COLS), "egg": [""] * len(eggnog_names)}

    def protein_row(caller, pid, contig, length):
        hit = hits.get(caller, {}).get(pid)
        egg = eggnog.get(caller, {}).get(pid)
        row = [caller, pid, contig, length]
        if hits:
            row += hit or empty["hit"]
        row += taxonomy.get(hit[0], empty["tax"]) if hit else empty["tax"]
        row += ref_annot.get(hit[0], empty["ref"]) if hit else empty["ref"]
        if eukulele:
            row += eukulele.get(caller, {}).get(pid, empty["euk"])
        if eggnog_names:
            row += egg or empty["egg"]
        rank = (float(hit[3]) if hit else -1.0, egg is not None, length)
        return row, rank

    best = {}           # contig -> (rank, row)
    n_proteins = {}
    with open(f"{args.prefix}protein_table.tsv", "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for caller, pid, contig, length in proteins:
            row, rank = protein_row(caller, pid, contig, length)
            writer.writerow(row)
            n_proteins[contig] = n_proteins.get(contig, 0) + 1
            if contig not in best or rank > best[contig][0]:
                best[contig] = (rank, row)

    rep_columns = [c for c in columns if c != "contig"]
    no_protein = [""] * len(rep_columns)
    with open(f"{args.prefix}contig_table.tsv", "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig", "length", "n_proteins", *rep_columns, *samples])
        for name, row in contigs.items():
            rep = best.get(name)
            rep_values = [v for c, v in zip(columns, rep[1]) if c != "contig"] if rep else no_protein
            writer.writerow([name, row[1], n_proteins.get(name, 0), *rep_values, *row[2:]])

    n_hit = sum(1 for (r, _) in best.values() if r[0] >= 0)
    print(f"contigs: {len(contigs)}; with a protein: {len(best)}; with a reference hit: {n_hit}")
    print(f"proteins: {len(proteins)} ({', '.join(callers)})")


if __name__ == "__main__":
    main()
