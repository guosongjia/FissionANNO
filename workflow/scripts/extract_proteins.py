#!/usr/bin/env python
"""Extract protein or CDS nucleotide sequences from merged GFF3 + genome FASTA.

Output: gzipped FASTA, one record per mRNA.
Header: >{mrna_id} gene={gene_id} [provenance tags]

Modes:
  protein (default) — translated amino-acid sequence (to first stop)
  cds               — spliced coding nucleotide sequence (ATG to stop, no intron)
"""
import argparse
import gzip
import sys
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq


def parse_attrs(attr_str):
    d = {}
    for field in attr_str.rstrip(";").split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def parse_gff(path):
    gene_attrs = {}
    mrna_parent = {}
    mrna_attrs = {}
    cds_list = defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, _, ftype, start, end, _, strand, phase, attrs = parts
            a = parse_attrs(attrs)
            if ftype == "gene":
                gene_attrs[a.get("ID", "")] = a
            elif ftype == "mRNA":
                mid = a.get("ID", "")
                mrna_parent[mid] = a.get("Parent", "")
                mrna_attrs[mid] = a
            elif ftype == "CDS":
                parent = a.get("Parent", "")
                cds_list[parent].append(
                    (seqid, int(start) - 1, int(end), strand, int(phase) if phase != "." else 0)
                )
    return gene_attrs, mrna_parent, mrna_attrs, cds_list


def extract_cds_seq(cds_records, genome):
    """Return spliced CDS nucleotide sequence (phase-trimmed, codon-aligned)."""
    strand = cds_records[0][3]
    cds_sorted = sorted(cds_records, key=lambda x: x[1])

    seq_parts = []
    for seqid, start, end, _, _ in cds_sorted:
        if seqid not in genome:
            return None
        seq_parts.append(genome[seqid][start:end])

    cds_seq = "".join(seq_parts)
    if strand == "-":
        cds_seq = str(Seq(cds_seq).reverse_complement())

    first_in_mrna = cds_sorted[0] if strand == "+" else cds_sorted[-1]
    first_phase = first_in_mrna[4]
    cds_seq = cds_seq[first_phase:]

    trim = len(cds_seq) % 3
    if trim:
        cds_seq = cds_seq[:-trim]

    if not cds_seq:
        return None

    return cds_seq


def extract_protein(cds_records, genome):
    cds_seq = extract_cds_seq(cds_records, genome)
    if not cds_seq:
        return None
    return str(Seq(cds_seq).translate(to_stop=True))


def provenance_str(mrna_a, gene_a):
    tags = []
    for key in ("source", "relation", "HGT_call", "SOG_id", "evidence"):
        val = mrna_a.get(key) or gene_a.get(key)
        if val:
            tags.append(f"{key}={val}")
    return " ".join(tags)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gff", required=True)
    p.add_argument("--genome", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--mode", choices=["protein", "cds"], default="protein",
                   help="protein = translated aa; cds = spliced nucleotide coding sequence")
    args = p.parse_args()

    genome = {r.id: str(r.seq) for r in SeqIO.parse(args.genome, "fasta")}
    gene_attrs, mrna_parent, mrna_attrs, cds_list = parse_gff(args.gff)

    extract_fn = extract_protein if args.mode == "protein" else extract_cds_seq
    label = "proteins" if args.mode == "protein" else "CDS sequences"

    written = skipped = filtered = 0
    with gzip.open(args.output, "wt") as out:
        for mrna_id, cds_records in cds_list.items():
            ma = mrna_attrs.get(mrna_id, {})
            if ma.get("valid_orf", "").lower() == "false":
                filtered += 1
                continue
            seq = extract_fn(cds_records, genome)
            if not seq:
                skipped += 1
                continue
            gene_id = mrna_parent.get(mrna_id, "")
            prov = provenance_str(ma, gene_attrs.get(gene_id, {}))
            header = f">{mrna_id} gene={gene_id}"
            if prov:
                header += f" {prov}"
            out.write(f"{header}\n{seq}\n")
            written += 1

    print(f"{args.sample}: {written} {label} written, {skipped} skipped, "
          f"{filtered} filtered (valid_orf=False)", file=sys.stderr)


if __name__ == "__main__":
    main()
