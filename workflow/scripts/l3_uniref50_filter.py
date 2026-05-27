#!/usr/bin/env python
"""Filter BRAKER4 predictions by UniRef50 DIAMOND hits.

Rules (CLAUDE.md §3.4):
  - no UniRef50 hit -> drop
  - any hit -> keep
  - top hit not Schizosaccharomyces -> tag HGT_call=putative_<top_taxon>
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict


def parse_braker_gff(gff_path):
    """Parse BRAKER4 GFF3 into gene hierarchies.

    Returns dict: gene_id -> {line_group: [lines], mRNA_ids: [...]}
    """
    genes = {}
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
            ftype = parts[2]
            attrs = parts[8]
            if ftype == "gene":
                m = re.search(r"ID=([^;]+)", attrs)
                if m:
                    current_gene = m.group(1)
                    gene_lines[current_gene].append(parts)
                    gene_order.append(current_gene)
            elif current_gene:
                gene_lines[current_gene].append(parts)

    return gene_lines, gene_order


def extract_cds_sequences(gene_lines, genome_path):
    """Extract and translate CDS sequences per gene from genome FASTA."""
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
        seq = ""
        for sid, s, e, st in cds_parts:
            seq += str(genome[sid].seq[s-1:e])
        if strand == "-":
            seq = str(Seq(seq).reverse_complement())
        prot = str(Seq(seq).translate())
        if prot.endswith("*"):
            prot = prot[:-1]
        proteins[gene_id] = prot

    return proteins


def run_diamond(proteins, diamond_db, evalue, threads, out_tsv):
    """Run diamond blastp and return path to results."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as f:
        query_fa = f.name
        for gid, seq in proteins.items():
            f.write(f">{gid}\n{seq}\n")

    cmd = [
        "diamond", "blastp",
        "--db", diamond_db,
        "--query", query_fa,
        "--out", out_tsv,
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length",
        "evalue", "bitscore", "stitle",
        "--evalue", str(evalue),
        "--threads", str(threads),
        "--max-target-seqs", "5",
    ]
    subprocess.run(cmd, check=True)
    os.unlink(query_fa)


def parse_diamond_hits(tsv_path, min_aln_len=80):
    """Parse DIAMOND output, return top hit per gene (filtered by min alignment length)."""
    top_hits = {}
    with open(tsv_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            qseqid = parts[0]
            aln_len = int(parts[3])
            if aln_len < min_aln_len:
                continue
            if qseqid not in top_hits:
                top_hits[qseqid] = {
                    "sseqid": parts[1],
                    "pident": parts[2],
                    "length": parts[3],
                    "evalue": parts[4],
                    "bitscore": parts[5],
                    "stitle": parts[6] if len(parts) > 6 else "",
                }
    return top_hits


def extract_taxon(stitle):
    """Extract taxon from UniRef50 stitle: '... Tax=Genus species TaxID=...'"""
    m = re.search(r"Tax=(.+?)\s+TaxID=", stitle)
    if m:
        return m.group(1).strip().replace(" ", "_")
    return "unknown"


SCHIZO_PATTERN = re.compile(r"_SCH(PO|PM|CR|JP)\b|Schizosaccharomyces", re.IGNORECASE)
TE_PATTERN = re.compile(
    r"\b(transposon|transposase|retrotransposon|integrase|polyprotein"
    r"|reverse[ _-]?transcriptase|RNA-directed DNA polymerase|Tf2[-_]?\d*)\b",
    re.IGNORECASE,
)
ANIMAL_PLANT_TAXA = {
    "metazoa", "viridiplantae", "mammalia", "aves", "eutheria",
    "insecta", "actinopteri", "magnoliopsida", "vertebrata",
    "tetrapoda", "amniota", "boreoeutheria", "primates", "rodentia",
    "carnivora", "passeriformes", "diptera", "lepidoptera",
    "embryophyta", "spermatophyta", "tracheophyta",
}
# Common SwissProt species codes for animal/plant contamination.
# Format: 5-letter codes (HUMAN, MOUSE) or 9XXXX (higher taxa within animal/plant).
ANIMAL_PLANT_REPID_CODES = {
    # Mammals
    "HUMAN", "MOUSE", "RAT", "BOVIN", "PIG", "SHEEP", "CAPHI", "CANLF", "FELCA",
    "HORSE", "RABIT", "MACMU", "PANTR", "GORGO",
    # Birds
    "CHICK", "TAEGU",
    # Other vertebrates
    "DANRE", "XENLA", "XENTR",
    # Insects / invertebrates
    "DROME", "CAEEL", "ANOGA", "AEDAE", "BOMMO",
    # Plants
    "ARATH", "ORYSJ", "ORYSI", "ZEAMA", "MAIZE", "SOYBN", "MEDTR", "SOLTU", "SOLLC",
    # Higher taxa codes (9XXXX) for animals/plants
    "9MAMM", "9EUTH", "9PRIM", "9CARN", "9PASS", "9ROSI", "9MONO", "9MAGN",
    "9METZ", "9VIRI", "9TETR", "9CETA", "9ARTI", "9FALC", "9MARS", "9PINI",
}


def extract_repid_code(stitle):
    """Extract organism code from UniRef RepID (e.g. RepID=VIME_HUMAN -> HUMAN)."""
    m = re.search(r"RepID=\S+_([A-Z0-9]{4,5})\b", stitle)
    return m.group(1) if m else None


def is_schpo_hit(hit):
    return bool(SCHIZO_PATTERN.search(hit["stitle"]))


def is_te_hit(hit):
    return bool(TE_PATTERN.search(hit["stitle"]))


def is_animal_plant_hit(hit):
    taxon_low = extract_taxon(hit["stitle"]).lower()
    if any(t in taxon_low for t in ANIMAL_PLANT_TAXA):
        return True
    code = extract_repid_code(hit["stitle"])
    if code and code in ANIMAL_PLANT_REPID_CODES:
        return True
    return False


def load_l1l2_seqid_set(gff_paths):
    """Return set of seqids that contain any L1/L2 gene (incl. pseudogene/truncated)."""
    seqids = set()
    for path in gff_paths:
        with open(path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9 or parts[2] != "gene":
                    continue
                seqids.add(parts[0])
    return seqids


def compute_bam_median_coverage(bam_path, region=None, threads=1):
    """Compute median read depth via `samtools depth -a`.

    region: None for whole-genome; or "seqid:start-end" string.
    """
    cmd = ["samtools", "depth", "-a"]
    if region:
        cmd += ["-r", region]
    cmd.append(bam_path)
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    depths = []
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3:
            depths.append(int(parts[2]))
    if not depths:
        return 0.0
    depths.sort()
    n = len(depths)
    return float(depths[n // 2]) if n % 2 else (depths[n // 2 - 1] + depths[n // 2]) / 2.0


def load_l1l2_intervals(gff_paths):
    """Load gene intervals from L1+L2 GFFs for overlap checking."""
    from collections import defaultdict
    intervals = defaultdict(list)
    for path in gff_paths:
        with open(path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9 or parts[2] != "gene":
                    continue
                intervals[parts[0]].append((int(parts[3]), int(parts[4])))
    for seqid in intervals:
        intervals[seqid].sort()
    return intervals


def overlaps_l1l2(gene_parts, l1l2_intervals):
    """Check if a BRAKER4 gene overlaps any L1/L2 gene by >= 50%."""
    gene_row = gene_parts[0]
    seqid = gene_row[0]
    g_start, g_end = int(gene_row[3]), int(gene_row[4])
    g_len = g_end - g_start + 1
    for s, e in l1l2_intervals.get(seqid, []):
        if s > g_end:
            break
        if e < g_start:
            continue
        ov = min(g_end, e) - max(g_start, s) + 1
        if ov >= 0.5 * g_len:
            return True
    return False


def add_provenance(parts_list, hit, taxon, hgt_call):
    """Add provenance attributes to gene-level features."""
    sseqid = hit["sseqid"]
    for parts in parts_list:
        if parts[2] in ("gene", "mRNA"):
            attrs = parts[8].rstrip(";")
            attrs += f";source=annevo_L3;evidence=UniRef50_{sseqid}"
            if hgt_call:
                attrs += f";HGT_call={hgt_call}"
            parts[8] = attrs


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--braker-gff", required=True)
    p.add_argument("--genome", required=True)
    p.add_argument("--l1-gff", required=True)
    p.add_argument("--l2-gff", required=True)
    p.add_argument("--diamond-tsv", required=True,
                   help="Pre-computed DIAMOND TSV from l3_diamond_search.py")
    p.add_argument("--schizo-keyword", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--out-gff", required=True)
    p.add_argument("--min-prot-len", type=int, default=200)
    p.add_argument("--min-aln-len", type=int, default=200)
    p.add_argument("--bam", default=None)
    p.add_argument("--coverage-floor", type=float, default=0.2)
    args = p.parse_args()

    gene_lines, gene_order = parse_braker_gff(args.braker_gff)

    if not gene_order:
        with open(args.out_gff, "w") as f:
            f.write("##gff-version 3\n")
        with open(args.out_diamond, "w") as f:
            pass
        print("No genes in BRAKER4 output; wrote empty files", file=sys.stderr)
        return

    # Remove genes overlapping L1+L2 annotations
    l1l2_ivs = load_l1l2_intervals([args.l1_gff, args.l2_gff])
    overlap_removed = 0
    filtered_order = []
    for gene_id in gene_order:
        if overlaps_l1l2(gene_lines[gene_id], l1l2_ivs):
            overlap_removed += 1
        else:
            filtered_order.append(gene_id)
    gene_order = filtered_order
    print(f"Removed {overlap_removed} genes overlapping L1/L2", file=sys.stderr)

    proteins = extract_cds_sequences(gene_lines, args.genome)

    # Filter by minimum protein length
    short_removed = 0
    long_proteins = {}
    for gid, seq in proteins.items():
        if len(seq) >= args.min_prot_len:
            long_proteins[gid] = seq
        else:
            short_removed += 1
    gene_order = [g for g in gene_order if g in long_proteins]
    proteins = long_proteins
    print(f"Removed {short_removed} genes < {args.min_prot_len} aa", file=sys.stderr)

    top_hits = parse_diamond_hits(args.diamond_tsv, min_aln_len=args.min_aln_len)

    # Pre-compute supporting state for downstream filters
    schpo_seqids = load_l1l2_seqid_set([args.l1_gff, args.l2_gff])
    assembly_median_cov = None  # lazily computed when needed

    kept = []
    counts = {"no_hit": 0, "schpo": 0, "te": 0,
              "metazoa_no_schizo_neighbor": 0, "metazoa_low_coverage": 0,
              "kept_schizo_neighbor_hgt": 0, "kept_other_hgt": 0}

    for gene_id in gene_order:
        hit = top_hits.get(gene_id)
        if not hit:
            counts["no_hit"] += 1
            continue
        if is_schpo_hit(hit):
            counts["schpo"] += 1
            continue
        if is_te_hit(hit):
            counts["te"] += 1
            continue

        gene_row = gene_lines[gene_id][0]
        seqid = gene_row[0]
        g_start, g_end = int(gene_row[3]), int(gene_row[4])
        taxon = extract_taxon(hit["stitle"])

        if is_animal_plant_hit(hit):
            # Test 1: contig must contain any L1/L2 Schizo gene
            if seqid not in schpo_seqids:
                counts["metazoa_no_schizo_neighbor"] += 1
                print(f"  drop {gene_id}: animal/plant hit on contig with no L1/L2 gene "
                      f"({taxon})", file=sys.stderr)
                continue
            # Test 2: coverage check
            if args.bam:
                if assembly_median_cov is None:
                    print("Computing assembly median coverage...", file=sys.stderr)
                    assembly_median_cov = compute_bam_median_coverage(args.bam,
                                                                      threads=args.threads)
                    print(f"  assembly median coverage = {assembly_median_cov:.1f}",
                          file=sys.stderr)
                region = f"{seqid}:{g_start}-{g_end}"
                gene_cov = compute_bam_median_coverage(args.bam, region=region,
                                                      threads=args.threads)
                threshold = args.coverage_floor * assembly_median_cov
                if gene_cov < threshold:
                    counts["metazoa_low_coverage"] += 1
                    print(f"  drop {gene_id}: animal/plant hit with low coverage "
                          f"({gene_cov:.1f} < {threshold:.1f}) ({taxon})",
                          file=sys.stderr)
                    continue
            # Passed both tests: keep with HGT tag
            add_provenance(gene_lines[gene_id], hit, taxon, f"putative_{taxon}_contam_check_passed")
            counts["kept_schizo_neighbor_hgt"] += 1
            kept.append(gene_id)
        else:
            # Other non-Schizo, non-TE, non-animal/plant hit: keep as HGT candidate
            add_provenance(gene_lines[gene_id], hit, taxon, f"putative_{taxon}")
            counts["kept_other_hgt"] += 1
            kept.append(gene_id)

    with open(args.out_gff, "w") as f:
        f.write("##gff-version 3\n")
        for n, gene_id in enumerate(kept, 1):
            new_gid = f"L3G_{args.sample}_{n:05d}"
            new_tid = f"L3T_{args.sample}_{n:05d}"
            for parts in gene_lines[gene_id]:
                p2 = list(parts)
                ftype = p2[2]
                attrs = p2[8]
                if ftype == "gene":
                    attrs = re.sub(r"ID=[^;]+", f"ID={new_gid}", attrs)
                elif ftype == "mRNA":
                    attrs = re.sub(r"ID=[^;]+", f"ID={new_tid}", attrs)
                    attrs = re.sub(r"Parent=[^;]+", f"Parent={new_gid}", attrs)
                else:
                    attrs = re.sub(r"Parent=[^;]+", f"Parent={new_tid}", attrs)
                p2[8] = attrs
                f.write("\t".join(p2) + "\n")

    print(f"\nKept {len(kept)} genes", file=sys.stderr)
    print(f"Filtered breakdown:", file=sys.stderr)
    for k, v in counts.items():
        print(f"  {k}: {v}", file=sys.stderr)


if __name__ == "__main__":
    main()
