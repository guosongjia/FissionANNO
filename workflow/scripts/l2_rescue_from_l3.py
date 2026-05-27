#!/usr/bin/env python
"""Rescue L2 singleton_no_sog candidates validated by L3 ANNEVO predictions.

Reads:
  - L2 sidecar TSV (intra_genus_HGT_candidates)
  - L2 miniprot raw GFF (to retrieve full gene structure for rescued loci)
  - ANNEVO remapped GFF (all predictions, not just UniRef50-filtered)
  - DIAMOND TSV (from L3 UniRef50 search)

For each sidecar entry tagged 'singleton_no_sog':
  - find any ANNEVO gene overlapping ≥ 50% reciprocally
  - if that ANNEVO gene has a DIAMOND hit (aln_len ≥ min_aln_len, query ≥ min_prot_len)
    AND the best miniprot hit at this locus has identity ≥ min_identity
    → rescue: append the L2 miniprot gene structure to l2_kept GFF with
              relation=non_reference_gene;sog_status=no_sog;validated_by=annevo_L3
              and remove the row from the sidecar

Outputs:
  - updated l2_kept GFF (input + rescued genes)
  - updated sidecar TSV (rescued rows removed)
  - rescue log TSV (one row per rescued locus)
"""
import argparse
import re
import sys
from collections import defaultdict


def parse_locus_string(locus_str):
    """Parse 'seqid:start-end(strand)' → (seqid, start, end, strand)."""
    m = re.match(r"^(.+):(\d+)-(\d+)\(([+-])\)$", locus_str)
    if not m:
        return None
    return m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)


def reciprocal_overlap(a_start, a_end, b_start, b_end):
    if a_end < b_start or b_end < a_start:
        return 0.0
    inter = min(a_end, b_end) - max(a_start, b_start) + 1
    a_len = a_end - a_start + 1
    b_len = b_end - b_start + 1
    return min(inter / a_len, inter / b_len)


def parse_gff_genes(gff_path):
    """Return list of (seqid, start, end, strand, gene_id, all_lines)."""
    genes = []
    current = None
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] == "gene":
                if current:
                    genes.append(current)
                m = re.search(r"ID=([^;]+)", parts[8])
                gid = m.group(1) if m else None
                current = {
                    "seqid": parts[0], "start": int(parts[3]), "end": int(parts[4]),
                    "strand": parts[6], "gid": gid, "lines": [parts],
                }
            elif current:
                current["lines"].append(parts)
    if current:
        genes.append(current)
    return genes


def parse_miniprot_loci(gff_path):
    """Group miniprot mRNA lines by locus. Returns dict: locus_key → [mRNA dicts].
    locus_key = (seqid, strand). Each mRNA dict has start/end/identity/src_protein/lines.
    """
    by_seq = defaultdict(list)
    pending_lines = []
    current_mrna = None
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] == "mRNA":
                if current_mrna:
                    current_mrna["lines"] = pending_lines
                    by_seq[(current_mrna["seqid"], current_mrna["strand"])].append(current_mrna)
                pending_lines = [parts]
                identity = re.search(r"Identity=([\d.]+)", parts[8])
                target = re.search(r"Target=(\S+)", parts[8])
                src = target.group(1) if target else ""
                current_mrna = {
                    "seqid": parts[0], "start": int(parts[3]), "end": int(parts[4]),
                    "strand": parts[6],
                    "identity": float(identity.group(1)) if identity else 0.0,
                    "src_protein": src,
                }
            elif current_mrna:
                pending_lines.append(parts)
    if current_mrna:
        current_mrna["lines"] = pending_lines
        by_seq[(current_mrna["seqid"], current_mrna["strand"])].append(current_mrna)
    return by_seq


def parse_diamond_hits(tsv_path, min_aln_len):
    hits = set()
    with open(tsv_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            try:
                aln = int(parts[3])
            except ValueError:
                continue
            if aln >= min_aln_len:
                hits.add(parts[0])
    return hits


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sidecar", required=True)
    p.add_argument("--miniprot-gff", required=True)
    p.add_argument("--annevo-gff", required=True,
                   help="Remapped ANNEVO GFF (full prediction, not UniRef50-filtered)")
    p.add_argument("--diamond-tsv", required=True)
    p.add_argument("--l2-kept-in", required=True)
    p.add_argument("--l2-kept-out", required=True)
    p.add_argument("--sidecar-out", required=True)
    p.add_argument("--rescue-log", required=True)
    p.add_argument("--min-aln-len", type=int, default=200)
    p.add_argument("--min-overlap", type=float, default=0.5)
    p.add_argument("--min-miniprot-identity", type=float, default=0.9)
    args = p.parse_args()

    # Read inputs
    annevo_genes = parse_gff_genes(args.annevo_gff)
    diamond_hit_qids = parse_diamond_hits(args.diamond_tsv, args.min_aln_len)
    miniprot_loci = parse_miniprot_loci(args.miniprot_gff)

    # Index ANNEVO by seqid (strand-agnostic for overlap test)
    annevo_by_seq = defaultdict(list)
    for g in annevo_genes:
        annevo_by_seq[g["seqid"]].append(g)

    # Read sidecar
    with open(args.sidecar) as f:
        header = f.readline().rstrip("\n")
        sidecar_rows = [line.rstrip("\n") for line in f if line.strip()]

    rescued_loci = []
    kept_sidecar = []
    next_gid_num = 100000  # IDs for rescued genes

    for row in sidecar_rows:
        cols = row.split("\t")
        if len(cols) < 6 or cols[5] != "singleton_no_sog":
            kept_sidecar.append(row)
            continue
        sample = cols[0]
        locus_str = cols[1]
        src_species = cols[2]
        src_protein = cols[3]
        identity = float(cols[4])

        loc = parse_locus_string(locus_str)
        if loc is None:
            kept_sidecar.append(row)
            continue
        seqid, l_start, l_end, l_strand = loc

        # Find ANNEVO gene overlapping this locus with reciprocal ≥ min_overlap
        # AND that ANNEVO gene has a DIAMOND hit with aln_len ≥ min_aln_len
        rescued = False
        for ag in annevo_by_seq.get(seqid, []):
            ovr = reciprocal_overlap(l_start, l_end, ag["start"], ag["end"])
            if ovr < args.min_overlap:
                continue
            if ag["gid"] not in diamond_hit_qids:
                continue
            # Find best miniprot mRNA at this locus to construct gene record
            best_m = None
            for m in miniprot_loci.get((seqid, l_strand), []):
                if m["start"] == l_start and m["end"] == l_end:
                    if best_m is None or m["identity"] > best_m["identity"]:
                        best_m = m
            if best_m is None or best_m["identity"] < args.min_miniprot_identity:
                continue

            next_gid_num += 1
            new_gid = f"L2G_{sample}_RESCUE_{next_gid_num}"
            new_tid = f"L2T_{sample}_RESCUE_{next_gid_num}"
            attrs_extra = (
                f";relation=non_reference_gene"
                f";sog_status=no_sog"
                f";validated_by=annevo_L3"
                f";SOG_id=NA"
                f";ref_id=NA"
                f";best_other_id={best_m['identity']:.4f}"
                f";best_other_sp={src_species}"
                f";source=miniprot_L2"
            )

            rescued_lines = []
            for parts in best_m["lines"]:
                p2 = list(parts)
                if p2[2] == "mRNA":
                    p2[8] = f"ID={new_tid};Parent={new_gid}{attrs_extra};orf_status=complete"
                else:
                    # Rewrite Parent to new mRNA id
                    new_attrs = re.sub(r"Parent=[^;]+", f"Parent={new_tid}", p2[8])
                    p2[8] = new_attrs
                rescued_lines.append(p2)
            # Synthesize a gene line (miniprot output has no gene feature)
            mrna_parts = best_m["lines"][0]
            gene_parts = [
                mrna_parts[0], "miniprot", "gene",
                str(best_m["start"]), str(best_m["end"]),
                mrna_parts[5], best_m["strand"], ".",
                f"ID={new_gid}{attrs_extra};orf_status=complete",
            ]
            rescued_lines.insert(0, gene_parts)
            rescued_loci.append({
                "lines": rescued_lines, "locus_str": locus_str,
                "src_protein": src_protein, "identity": identity,
                "annevo_gid": ag["gid"], "new_gid": new_gid,
            })
            rescued = True
            break

        if not rescued:
            kept_sidecar.append(row)

    # Write outputs
    with open(args.l2_kept_in) as fin, open(args.l2_kept_out, "w") as fout:
        fout.write(fin.read())
        if not fout.tell() or fout.tell() and fin.read() == "":
            pass
        for entry in rescued_loci:
            for parts in entry["lines"]:
                fout.write("\t".join(parts) + "\n")

    with open(args.sidecar_out, "w") as f:
        f.write(header + "\n")
        for row in kept_sidecar:
            f.write(row + "\n")

    with open(args.rescue_log, "w") as f:
        f.write("locus\tsrc_protein\tminiprot_identity\tannevo_gene_id\tnew_gene_id\n")
        for r in rescued_loci:
            f.write(f"{r['locus_str']}\t{r['src_protein']}\t{r['identity']:.4f}\t"
                    f"{r['annevo_gid']}\t{r['new_gid']}\n")

    print(f"Rescued {len(rescued_loci)} singleton_no_sog loci as non_reference_gene",
          file=sys.stderr)


if __name__ == "__main__":
    main()
