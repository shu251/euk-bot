#!/usr/bin/env python3
"""Merge MarFERReT, EukProt and the prokaryotes of PhyloDB into one protein FASTA plus a taxonomy table,
the inputs for euk-bot's --reference_dmnd (after diamond makedb) and --reference_taxonomy.

Protein IDs get a source prefix (MFT_, EUK_, PDB_) and '|' becomes '_'. A sequence identical (ignoring a
trailing '*') to one from an earlier source is dropped, which removes the MMETSP data shared by MarFERReT and
EukProt; sources are read in priority order MarFERReT, EukProt, PhyloDB. Duplicates within a source are kept. Eukaryote lineages are PR2-style ranks: MarFERReT's own PR2 taxonomy, and for EukProt the PR2
lineage of the same genus (from the PR2 taxonomy table), falling back to EukProt's UniEuk taxonomy.

Outputs: <out>.faa.gz, <out>.taxonomy.tsv.gz
  taxonomy columns: protein_id, source, organism, ncbi_taxid, lineage, lineage_source, pfam_name, pfam_id
"""
import argparse
import csv
import gzip
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

TAX_COLS = ["protein_id", "source", "organism", "ncbi_taxid", "lineage", "lineage_source", "pfam_name", "pfam_id"]


def open_text(path):
    return gzip.open(path, "rt", newline="") if str(path).endswith(".gz") else open(path, newline="")


def read_fasta(path):
    header, seq = None, []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:].split()[0], []
            elif line:
                seq.append(line)
    if header is not None:
        yield header, "".join(seq)


class Writer:
    """Writes FASTA + taxonomy rows, skipping sequences already written by an earlier source."""

    def __init__(self, out):
        self.faa = gzip.open(f"{out}.faa.gz", "wt", compresslevel=3)
        self.tax_fh = gzip.open(f"{out}.taxonomy.tsv.gz", "wt", newline="", compresslevel=3)
        self.tax = csv.writer(self.tax_fh, delimiter="\t", lineterminator="\n")
        self.tax.writerow(TAX_COLS)
        self.seen = set()       # digests from earlier sources
        self.current = set()    # digests from the source being added
        self.stats = Counter()

    def add(self, source, protein_id, seq, row):
        digest = hashlib.md5(seq.rstrip("*").encode()).digest()
        if digest in self.seen:
            self.stats[f"{source} duplicate of earlier source"] += 1
            return
        self.current.add(digest)
        pid = f"{source}_{protein_id.replace('|', '_')}"
        self.faa.write(f">{pid}\n{seq}\n")
        self.tax.writerow([pid, source, *row])
        self.stats[f"{source} written"] += 1

    def end_source(self):
        self.seen |= self.current
        self.current = set()

    def close(self):
        self.faa.close()
        self.tax_fh.close()


def pr2_genus_lineages(path):
    """Genus -> most common PR2 lineage (ranks down to genus). Finds the 'Eukaryota;...' field on each line."""
    by_genus = defaultdict(Counter)
    with open_text(path) as fh:
        for line in fh:
            for field in line.rstrip("\n").split("\t"):
                if field.startswith("Eukaryota;"):
                    ranks = [r for r in field.strip(";").split(";") if r]
                    if len(ranks) >= 8:
                        by_genus[ranks[7]][";".join(ranks[:8])] += 1
                    break
    return {genus: lineages.most_common(1)[0][0] for genus, lineages in by_genus.items()}


def add_marferret(writer, mft_dir):
    meta = {}
    with open_text(mft_dir / "MarFERReT.v1.1.1.metadata.csv") as fh:
        for row in csv.DictReader(fh):
            meta[row["entry_id"]] = (row["marferret_name"], row["tax_id"], row["pr2_taxonomy"])
    entry_of = {}
    with open_text(mft_dir / "MarFERReT.v1.1.1.proteins_info.tab.gz") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            entry_of[row["aa_id"]] = row["entry_id"]
    pfam = {}
    with open_text(mft_dir / "MarFERReT.v1.1.1.best_pfam_annotations.csv.gz") as fh:
        for row in csv.DictReader(fh):
            pfam[row["aa_id"]] = (row["pfam_name"], row["pfam_id"])

    for aa_id, seq in read_fasta(mft_dir / "MarFERReT.v1.1.1.proteins.faa.gz"):
        name, taxid, lineage = meta.get(entry_of.get(aa_id), ("", "", ""))
        writer.add("MFT", aa_id, seq, [name, taxid, lineage, "PR2" if lineage else "", *pfam.get(aa_id, ("", ""))])


def add_eukprot(writer, eukprot_dir, pr2_lineage_of_genus):
    datasets = {}
    table = next(eukprot_dir.glob("EukProt_included_data_sets*.txt"))
    with open_text(table) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            genus, name = row["Genus_UniEuk"].strip('"'), row["Name_to_Use"]
            if genus in pr2_lineage_of_genus:
                lineage, lineage_source = f"{pr2_lineage_of_genus[genus]};{name}", "PR2_genus"
            else:
                lineage, lineage_source = row["Taxonomy_UniEuk"], "UniEuk"
            datasets[row["EukProt_ID"]] = (name, lineage, lineage_source)

    for fasta in sorted((eukprot_dir / "proteins").glob("*.fasta")):
        name, lineage, lineage_source = datasets.get(fasta.name[:7], (fasta.stem[8:], "", ""))
        for pid, seq in read_fasta(fasta):
            writer.add("EUK", pid, seq, [name, "", lineage, lineage_source, "", ""])


def add_phylodb_prokaryotes(writer, phylodb_dir, phylodb_faa):
    def norm(strain):
        return re.sub(r"[\"']", "", strain).strip()

    lineage_of_strain = {}
    with open_text(phylodb_dir / "taxonomy-table.txt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["taxonomy"].startswith(("Bacteria", "Archaea")):
                lineage_of_strain[norm(row["strain_name"])] = row["taxonomy"]
    print(f"PhyloDB prokaryote strains: {len(lineage_of_strain)}", file=sys.stderr)

    with open(phylodb_dir / "prot-map.json") as fh:
        strain_of = {pid: norm(s) for pid, s in json.load(fh).items() if norm(s) in lineage_of_strain}

    for pid, seq in read_fasta(phylodb_faa):
        strain = strain_of.get(pid)
        if strain is None:
            writer.stats["PDB not prokaryote / unmapped"] += 1
            continue
        writer.add("PDB", pid, seq, [strain, "", lineage_of_strain[strain], "PhyloDB", "", ""])


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--marferret", type=Path, required=True, help="MarFERReT data/ directory")
    parser.add_argument("--eukprot", type=Path, required=True, help="EukProt v3 directory (proteins/ + data set table)")
    parser.add_argument("--phylodb", type=Path, required=True, help="PhyloDB directory (prot-map.json, taxonomy-table.txt)")
    parser.add_argument("--phylodb-faa", type=Path, help="PhyloDB protein FASTA (default: <phylodb>/reference.pep.fa.gz)")
    parser.add_argument("--pr2-taxonomy", type=Path, required=True, help="PR2 taxonomy table, for EukProt genus lookup")
    parser.add_argument("--out", required=True, help="output prefix")
    args = parser.parse_args()

    pr2 = pr2_genus_lineages(args.pr2_taxonomy)
    print(f"PR2 genera: {len(pr2)}", file=sys.stderr)

    writer = Writer(args.out)
    for label, step in [
        ("MarFERReT", lambda: add_marferret(writer, args.marferret)),
        ("EukProt", lambda: add_eukprot(writer, args.eukprot, pr2)),
        ("PhyloDB", lambda: add_phylodb_prokaryotes(writer, args.phylodb, args.phylodb_faa or args.phylodb / "reference.pep.fa.gz")),
    ]:
        print(f"== adding {label}", file=sys.stderr)
        step()
        writer.end_source()
        print("   " + ", ".join(f"{k}: {v}" for k, v in sorted(writer.stats.items())), file=sys.stderr)
    writer.close()


if __name__ == "__main__":
    main()
