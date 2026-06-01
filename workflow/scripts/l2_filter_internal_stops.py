#!/usr/bin/env python
"""Post-filter existing l2_kept.gff3 files: remove genes with internal stop codons.

Reads each sample's l2_kept.gff3 + genome FASTA, translates CDS, removes genes
containing internal stops, appends removed entries to the sidecar TSV.
Overwrites l2_kept.gff3 in place.

Usage:
    python l2_filter_internal_stops.py --manifest config/manifest.tsv \
        --outdir /path/to/annotation_Spombe
"""
import argparse
import os
import sys
from collections import defaultdict
from typing import Dict, List, Tuple

_CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}
_COMP = str.maketrans('ACGTacgt', 'TGCAtgca')


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def load_genome(fa_path: str) -> Dict[str, str]:
    genome: Dict[str, str] = {}
    cur_id = None
    chunks: List[str] = []
    with open(fa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_id is not None:
                    genome[cur_id] = "".join(chunks)
                cur_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if cur_id is not None:
        genome[cur_id] = "".join(chunks)
    return genome


def translate_cds_from_lines(cds_lines: List[str], strand: str, genome: Dict[str, str]) -> str:
    parts = []
    for line in cds_lines:
        cols = line.split("\t")
        parts.append((int(cols[3]), int(cols[4]), cols[7]))
    if strand == "+":
        parts.sort(key=lambda x: x[0])
    else:
        parts.sort(key=lambda x: x[0], reverse=True)
    nuc = ""
    seqid = cds_lines[0].split("\t")[0]
    contig_seq = genome.get(seqid, "")
    for i, (s, e, phase) in enumerate(parts):
        seg = contig_seq[s - 1:e]
        if strand == "-":
            seg = revcomp(seg)
        if i == 0:
            ph = int(phase) if phase.isdigit() else 0
            seg = seg[ph:]
        nuc += seg
    protein = []
    for i in range(0, len(nuc) - 2, 3):
        codon = nuc[i:i+3].upper()
        protein.append(_CODON_TABLE.get(codon, 'X'))
    return "".join(protein)


def filter_gff(gff_path: str, sample: str, sidecar_path: str,
               genome: Dict[str, str], record_to_sidecar: bool) -> Tuple[int, int]:
    """Filter one GFF file in place. Returns (kept_count, removed_count)."""
    if not os.path.exists(gff_path):
        return 0, 0

    gene_blocks: List[Tuple[str, List[str]]] = []
    header_lines: List[str] = []
    current_gene_id = None
    current_lines: List[str] = []

    with open(gff_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                header_lines.append(line)
                continue
            if not line.strip():
                continue
            cols = line.split("\t")
            if cols[2] == "gene":
                if current_gene_id:
                    gene_blocks.append((current_gene_id, current_lines))
                attrs = dict(x.split("=", 1) for x in cols[8].split(";") if "=" in x)
                current_gene_id = attrs.get("ID", "")
                current_lines = [line]
            else:
                current_lines.append(line)
    if current_gene_id:
        gene_blocks.append((current_gene_id, current_lines))

    kept_blocks: List[Tuple[str, List[str]]] = []
    removed = []
    for gene_id, lines in gene_blocks:
        cds_lines = [l for l in lines if l.split("\t")[2] == "CDS"]
        if not cds_lines:
            kept_blocks.append((gene_id, lines))
            continue
        strand = cds_lines[0].split("\t")[6]
        protein = translate_cds_from_lines(cds_lines, strand, genome)
        if "*" in protein[:-1]:
            n_stops = protein[:-1].count("*")
            mrna_line = next((l for l in lines if l.split("\t")[2] == "mRNA"), lines[0])
            mcols = mrna_line.split("\t")
            attrs = dict(x.split("=", 1) for x in mcols[8].split(";") if "=" in x)
            coord = f"{mcols[0]}:{mcols[3]}-{mcols[4]}({mcols[6]})"
            target = attrs.get("Target", "")
            tparts = target.split()
            src_full = tparts[0] if tparts else ""
            sp, prot = (src_full.split("|", 1) if "|" in src_full else ("", src_full))
            identity = attrs.get("Identity", "NA")
            relation = attrs.get("relation", "")
            removed.append(f"{sample}\t{coord}\t{sp}\t{prot}\t{identity}\t"
                           f"internal_stop\trelation={relation};n_internal_stops={n_stops}\n")
        else:
            kept_blocks.append((gene_id, lines))

    with open(gff_path, "w") as f:
        for h in header_lines:
            f.write(h + "\n")
        for gene_id, lines in kept_blocks:
            for l in lines:
                f.write(l + "\n")

    if removed and record_to_sidecar:
        with open(sidecar_path, "a") as f:
            for row in removed:
                f.write(row)

    return len(kept_blocks), len(removed)


def process_sample(sample: str, genome_fa: str, outdir: str) -> Tuple[int, int]:
    sidecar_path = os.path.join(outdir, "results", sample, "sidecar",
                                f"{sample}.intra_genus_HGT_candidates.tsv")
    kept_gff = os.path.join(outdir, "results", sample, "L2", f"{sample}.l2_kept.gff3")
    rescued_gff = os.path.join(outdir, "results", sample, "L2",
                               f"{sample}.l2_kept_rescued.gff3")

    if not os.path.exists(kept_gff) and not os.path.exists(rescued_gff):
        return 0, 0

    genome = load_genome(genome_fa)

    # Filter both files; only record to sidecar from the first one (kept) to avoid dupes
    kept1, removed1 = filter_gff(kept_gff, sample, sidecar_path, genome,
                                  record_to_sidecar=True)
    kept2, removed2 = filter_gff(rescued_gff, sample, sidecar_path, genome,
                                  record_to_sidecar=False)

    # Return the rescued file's totals (it's the one merge consumes)
    return kept2 if os.path.exists(rescued_gff) else kept1, max(removed1, removed2)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--manifest", required=True)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    samples = []
    with open(args.manifest) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            samples.append((parts[0], parts[1]))

    total_kept, total_removed = 0, 0
    for sample, genome_fa in samples:
        kept, removed = process_sample(sample, genome_fa, args.outdir)
        if removed:
            print(f"  {sample}: kept={kept}, removed={removed}")
        total_kept += kept
        total_removed += removed

    print(f"\nTotal: kept={total_kept}, removed={total_removed} (across {len(samples)} strains)")


if __name__ == "__main__":
    main()
