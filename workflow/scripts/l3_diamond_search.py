#!/usr/bin/env python
"""Run DIAMOND blastp on ANNEVO predictions against UniRef50.

Outputs a TSV that downstream rules (l2_rescue, l3_uniref50_filter) consume.
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict


def parse_braker_gff(gff_path):
    current_gene = None
    gene_lines = defaultdict(list)
    gene_order = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] == "gene":
                m = re.search(r"ID=([^;]+)", parts[8])
                if m:
                    current_gene = m.group(1)
                    gene_lines[current_gene].append(parts)
                    gene_order.append(current_gene)
            elif current_gene:
                gene_lines[current_gene].append(parts)
    return gene_lines, gene_order


def extract_cds_sequences(gene_lines, genome_path):
    from Bio import SeqIO
    from Bio.Seq import Seq

    genome = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))
    proteins = {}
    for gene_id, parts_list in gene_lines.items():
        cds_parts = [(p[0], int(p[3]), int(p[4]), p[6])
                     for p in parts_list if p[2] == "CDS"]
        if not cds_parts:
            continue
        seqid = cds_parts[0][0]
        strand = cds_parts[0][3]
        cds_parts.sort(key=lambda x: x[1])
        seq = "".join(str(genome[sid].seq[s-1:e]) for sid, s, e, _ in cds_parts)
        if strand == "-":
            seq = str(Seq(seq).reverse_complement())
        prot = str(Seq(seq).translate())
        if prot.endswith("*"):
            prot = prot[:-1]
        proteins[gene_id] = prot
    return proteins


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gff", required=True)
    p.add_argument("--genome", required=True)
    p.add_argument("--diamond-db", required=True)
    p.add_argument("--evalue", type=float, required=True)
    p.add_argument("--min-prot-len", type=int, default=200)
    p.add_argument("--threads", type=int, required=True)
    p.add_argument("--out-diamond", required=True)
    args = p.parse_args()

    gene_lines, gene_order = parse_braker_gff(args.gff)
    if not gene_order:
        open(args.out_diamond, "w").close()
        print("No genes in input GFF; wrote empty DIAMOND TSV", file=sys.stderr)
        return

    proteins = extract_cds_sequences(gene_lines, args.genome)
    proteins = {gid: seq for gid, seq in proteins.items() if len(seq) >= args.min_prot_len}
    print(f"Querying {len(proteins)} proteins (≥{args.min_prot_len} aa) vs UniRef50",
          file=sys.stderr)

    if not proteins:
        open(args.out_diamond, "w").close()
        return

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as f:
        query_fa = f.name
        for gid, seq in proteins.items():
            f.write(f">{gid}\n{seq}\n")

    cmd = [
        "diamond", "blastp",
        "--db", args.diamond_db,
        "--query", query_fa,
        "--out", args.out_diamond,
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length",
        "evalue", "bitscore", "stitle",
        "--evalue", str(args.evalue),
        "--threads", str(args.threads),
        "--max-target-seqs", "5",
    ]
    subprocess.run(cmd, check=True)
    os.unlink(query_fa)


if __name__ == "__main__":
    main()
