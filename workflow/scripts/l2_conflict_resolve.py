#!/usr/bin/env python
"""L2 conflict resolution: classify miniprot hits not overlapping L1 genes.

Architecture (v6):
  1. Pre-filter mRNAs by alignment length and identity.
  2. Drop any mRNA that overlaps an L1 gene (≥ overlap_min fraction of L covered).
     L1 lifton already handles these loci; miniprot adds nothing reliable here.
  3. Remaining (non-overlapping) candidates are collapsed into loci and classified
     by multi-species SOG comparison into:
       non_reference_gene / missing_lift / intra_genus_HGT_from_X /
       diverged_paralog_or_misannot
     Any intra_genus_HGT call is gated on candidate identity exceeding the
     SOG's precomputed max Lcon × donor pairwise identity (sog_lcon_max_id).
     Otherwise demoted to non_reference_gene (candidate is not more similar
     to donor than Lcon's own ortholog already is).
  4. ORF completeness filter: mRNAs without a stop_codon, OR whose translated CDS
     is shorter than orf_min_coverage * query_protein_length, go to sidecar.
  5. GFF output has full gene/mRNA/exon/CDS structure (miniprot omits gene+exon).
"""
import argparse
import logging
import os
import pickle
import re
from collections import defaultdict, OrderedDict
from typing import Dict, List, Tuple, Optional, Set


ATTR_RE = re.compile(r"([^=;]+)=([^;]*)")

def parse_attrs(s: str) -> Dict[str, str]:
    return {m.group(1): m.group(2) for m in ATTR_RE.finditer(s)}


def read_l1_genes(path: str) -> Dict[str, Dict]:
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] not in ("gene", "pseudogene"):
                continue
            attrs = parse_attrs(cols[8])
            gid = attrs.get("ID")
            if not gid:
                continue
            out[gid] = {
                "seqid": cols[0], "start": int(cols[3]),
                "end": int(cols[4]), "strand": cols[6], "type": cols[2],
            }
    return out


def read_miniprot_mrnas(path: str, min_aln_aa: int, min_identity: float) -> Dict[str, Dict]:
    out: Dict[str, Dict] = {}
    children: Dict[str, List[str]] = defaultdict(list)
    n_filtered = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            attrs = parse_attrs(cols[8])
            if cols[2] == "mRNA":
                mid = attrs.get("ID")
                if not mid:
                    continue
                identity = float(attrs.get("Identity", "0") or 0.0)
                target = attrs.get("Target", "")
                parts = target.split()
                src_full = parts[0] if parts else ""
                aa_start = int(parts[1]) if len(parts) >= 3 else 1
                aln_aa = (int(parts[2]) - int(parts[1]) + 1) if len(parts) >= 3 else 0
                if aln_aa < min_aln_aa or identity < min_identity:
                    n_filtered += 1
                    continue
                if "|" in src_full:
                    src_species, src_protein = src_full.split("|", 1)
                else:
                    src_species, src_protein = "", src_full
                out[mid] = {
                    "seqid": cols[0], "start": int(cols[3]), "end": int(cols[4]),
                    "strand": cols[6], "score": cols[5],
                    "rank": int(attrs.get("Rank", "0") or 0),
                    "identity": identity,
                    "aa_start": aa_start,
                    "aln_aa": aln_aa,
                    "src_species": src_species, "src_protein": src_protein,
                    "raw_line": line, "child_lines": [],
                }
            elif cols[2] in ("CDS", "stop_codon"):
                pid = attrs.get("Parent")
                if pid:
                    children[pid].append(line)
    for mid, lines in children.items():
        if mid in out:
            out[mid]["child_lines"] = lines
    logging.info(f"  pre-filter removed {n_filtered} mRNAs (min_aln_aa={min_aln_aa}, min_id={min_identity})")
    return out


def l_coverage(l_start: int, l_end: int, m_start: int, m_end: int) -> float:
    """Fraction of L1 gene length covered by M (one-sided)."""
    if l_end < m_start or m_end < l_start:
        return 0.0
    inter = min(l_end, m_end) - max(l_start, m_start) + 1
    return inter / (l_end - l_start + 1)


def build_l1_index(l1_genes: Dict[str, Dict]) -> Dict[Tuple[str,str], List[Tuple[int,int,str]]]:
    idx: Dict[Tuple[str,str], List[Tuple[int,int,str]]] = defaultdict(list)
    for gid, g in l1_genes.items():
        idx[(g["seqid"], g["strand"])].append((g["start"], g["end"], gid))
    for k in idx:
        idx[k].sort()
    return idx


def overlaps_l1(seqid: str, strand: str, start: int, end: int,
                l1_index, overlap_min: float) -> bool:
    """Return True if M overlaps any L1 gene by ≥ overlap_min fraction of L."""
    for s, e, gid in l1_index.get((seqid, strand), []):
        if e < start:
            continue
        if s > end:
            break
        if l_coverage(s, e, start, end) >= overlap_min:
            return True
    return False


def find_adjacent_l(seqid: str, strand: str, start: int, end: int,
                    l1_index, max_bp: int) -> Optional[str]:
    closest, closest_dist = None, None
    for s, e, gid in l1_index.get((seqid, strand), []):
        if e < start:
            d = start - e
            if d <= max_bp and (closest_dist is None or d < closest_dist):
                closest, closest_dist = gid, d
        elif s > end:
            d = s - end
            if d <= max_bp and (closest_dist is None or d < closest_dist):
                closest, closest_dist = gid, d
            else:
                break
    return closest


def collapse_loci(mrnas: List[Dict]) -> List[List[Dict]]:
    """Collapse mRNAs into loci by spatial overlap (per seqid+strand)."""
    by_key: Dict[Tuple[str,str], List[Tuple[int,int,Dict]]] = defaultdict(list)
    for m in mrnas:
        by_key[(m["seqid"], m["strand"])].append((m["start"], m["end"], m))
    groups: List[List[Dict]] = []
    for key, items in by_key.items():
        items.sort(key=lambda x: (x[0], x[1]))
        cur: List[Dict] = []
        cur_end = -1
        for s, e, m in items:
            if cur and s <= cur_end:
                cur.append(m)
                cur_end = max(cur_end, e)
            else:
                if cur:
                    groups.append(cur)
                cur = [m]
                cur_end = e
        if cur:
            groups.append(cur)
    return groups


def sog_of_protein(protein: str, sog_idx: Dict) -> Optional[str]:
    return sog_idx["protein_to_sog"].get(protein)


def classify_locus(locus_hits: List[Dict], ref_present_in_og_m: bool,
                   id_adv_pp: float, diverged_max_id: float,
                   og_m: str, sog_idx: Dict, short_to_full: Dict[str, str],
                   ref_short: str) -> Tuple[str, Dict]:
    best_per_species: Dict[str, Dict] = {}
    for h in locus_hits:
        sp = h["src_species"]
        if not sp:
            continue
        if sp not in best_per_species or h["identity"] > best_per_species[sp]["identity"]:
            best_per_species[sp] = h

    ref_id = best_per_species[ref_short]["identity"] if ref_short in best_per_species else None
    others = {sp: h["identity"] for sp, h in best_per_species.items() if sp != ref_short}

    if not ref_present_in_og_m:
        relation = "non_reference_gene"
        ev = {"ref_id": ref_id,
              "best_other_id": max(others.values()) if others else None,
              "best_other_sp": max(others, key=others.get) if others else None}
    elif ref_id is None:
        if others:
            best_sp = max(others, key=others.get)
            relation = f"intra_genus_HGT_from_{best_sp}"
            ev = {"ref_id": None, "best_other_id": others[best_sp], "best_other_sp": best_sp}
        else:
            relation = "diverged_paralog_or_misannot"
            ev = {"ref_id": None, "best_other_id": None, "best_other_sp": None}
    elif not others:
        relation = "missing_lift"
        ev = {"ref_id": ref_id, "best_other_id": None, "best_other_sp": None}
    else:
        best_sp = max(others, key=others.get)
        best_other_id = others[best_sp]
        advantage_pp = (best_other_id - ref_id) * 100.0
        if max(ref_id, best_other_id) < diverged_max_id:
            relation = "diverged_paralog_or_misannot"
            ev = {"ref_id": ref_id, "best_other_id": best_other_id, "best_other_sp": best_sp}
        elif advantage_pp >= id_adv_pp:
            relation = f"intra_genus_HGT_from_{best_sp}"
            ev = {"ref_id": ref_id, "best_other_id": best_other_id, "best_other_sp": best_sp}
        else:
            relation = "missing_lift"
            ev = {"ref_id": ref_id, "best_other_id": best_other_id, "best_other_sp": best_sp}

    # Pairwise cutoff gate: HGT requires candidate identity > max(ref × donor) in SOG
    if relation.startswith("intra_genus_HGT"):
        best_sp = ev.get("best_other_sp")
        sog_pw = sog_idx.get("sog_ref_max_id", {}).get(og_m, {})
        cutoff = sog_pw.get(best_sp)
        candidate_id = ev.get("best_other_id") or 0.0
        if cutoff is None or candidate_id <= cutoff:
            relation = "non_reference_gene"

    return relation, ev


def emit_gene_structure(m: Dict, gene_id: str, mrna_id: str,
                        extra_attrs: OrderedDict) -> List[str]:
    """Emit gene / mRNA / exon / CDS [/ stop_codon] lines for one miniprot hit.

    miniprot GFF3 has mRNA + CDS + stop_codon but no gene or exon rows.
    We synthesise gene and exon so downstream tools (e.g. gffutils) see a valid
    gene model.  Exon coordinates mirror each CDS feature (no introns in
    miniprot output).
    """
    raw_mrna = m["raw_line"].rstrip("\n").split("\t")
    seqid, src, strand = raw_mrna[0], raw_mrna[1], raw_mrna[6]
    start, end = m["start"], m["end"]
    score = raw_mrna[5]

    # gene line
    gene_attrs = f"ID={gene_id};" + ";".join(f"{k}={v}" for k, v in extra_attrs.items())
    gene_line = "\t".join([seqid, src, "gene", str(start), str(end), score, strand, ".", gene_attrs])

    # mRNA line: keep original attributes, patch ID/Parent, append extra
    mrna_attrs_str = raw_mrna[8]
    mrna_attrs_str = re.sub(r"(?:^|;)ID=[^;]*", "", mrna_attrs_str).lstrip(";")
    mrna_attrs_str = re.sub(r"(?:^|;)Parent=[^;]*", "", mrna_attrs_str).lstrip(";")
    mrna_suffix = ";".join(f"{k}={v}" for k, v in extra_attrs.items())
    mrna_attrs = f"ID={mrna_id};Parent={gene_id};{mrna_attrs_str};{mrna_suffix}"
    raw_mrna[8] = mrna_attrs
    mrna_line = "\t".join(raw_mrna)

    lines = [gene_line, mrna_line]

    cds_lines = []
    stop_lines = []
    for cl in m["child_lines"]:
        cols = cl.rstrip("\n").split("\t")
        if len(cols) < 9:
            continue
        child_attrs = parse_attrs(cols[8])
        child_attrs["Parent"] = mrna_id
        new_attr_str = ";".join(f"{k}={v}" for k, v in child_attrs.items())
        cols[8] = new_attr_str
        if cols[2] == "CDS":
            cds_lines.append("\t".join(cols))
            # exon mirrors CDS (miniprot produces single-exon models)
            exon_cols = cols[:]
            exon_cols[2] = "exon"
            exon_cols[7] = "."  # exon has no phase
            exon_attrs = {k: v for k, v in child_attrs.items()
                          if k not in ("Identity", "StopCodon", "Rank")}
            exon_cols[8] = ";".join(f"{k}={v}" for k, v in exon_attrs.items())
            lines.append("\t".join(exon_cols))
        elif cols[2] == "stop_codon":
            stop_lines.append("\t".join(cols))

    lines.extend(cds_lines)
    lines.extend(stop_lines)
    return lines


def has_stop_codon(m: Dict) -> bool:
    return any(
        cl.split("\t")[2] == "stop_codon"
        for cl in m["child_lines"]
        if len(cl.split("\t")) > 2
    )


def cds_aa_length(m: Dict) -> int:
    bp = 0
    for cl in m["child_lines"]:
        cols = cl.split("\t")
        if len(cols) >= 5 and cols[2] == "CDS":
            bp += int(cols[4]) - int(cols[3]) + 1
    return bp // 3


def load_protein_lengths(fa_path: str) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    cur_id, cur_len = None, 0
    with open(fa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_id is not None:
                    lengths[cur_id] = cur_len
                cur_id = line[1:].split()[0]
                cur_len = 0
            else:
                cur_len += len(line.strip())
        if cur_id is not None:
            lengths[cur_id] = cur_len
    return lengths


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample", required=True)
    p.add_argument("--l1-gff", required=True)
    p.add_argument("--l2-gff", required=True)
    p.add_argument("--sog-index", required=True)
    p.add_argument("--overlap-min", type=float, required=True)
    p.add_argument("--id-adv-pp", type=float, required=True)
    p.add_argument("--diverged-max-id", type=float, required=True)
    p.add_argument("--adjacent-max-bp", type=int, required=True)
    p.add_argument("--min-aln-aa", type=int, default=50)
    p.add_argument("--min-identity", type=float, default=0.3)
    p.add_argument("--hgt-min-identity", type=float, default=0.9)
    p.add_argument("--lcon-species", default="Schizosaccharomyces_pombe",
                   help="full species name of the reference species in the SOG table")
    p.add_argument("--protein-fa", required=True,
                   help="combined 9-species protein fasta (for query lengths)")
    p.add_argument("--orf-min-coverage", type=float, default=0.95,
                   help="translated CDS aa / query protein aa lower bound")
    p.add_argument("--out-gff", required=True)
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--out-sidecar", required=True)
    p.add_argument("--log-level", default="INFO")
    args = p.parse_args()

    ref_full = args.lcon_species
    ref_short = ref_full.replace("Schizosaccharomyces_", "S_")

    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO),
                        format="%(levelname)s: %(message)s")

    logging.info(f"reading L1 {args.l1_gff}")
    l1_genes = read_l1_genes(args.l1_gff)
    l1_index = build_l1_index(l1_genes)
    logging.info(f"  {len(l1_genes)} L1 genes")

    logging.info(f"reading L2 {args.l2_gff}")
    mrnas = read_miniprot_mrnas(args.l2_gff, args.min_aln_aa, args.min_identity)
    logging.info(f"  {len(mrnas)} mRNAs after pre-filter")

    with open(args.sog_index, "rb") as fh:
        sog_idx = pickle.load(fh)
    logging.info(f"loaded SOG index: {sog_idx['n_sogs']} SOGs")

    sog_ref_members: Dict[str, Set[str]] = {}
    for sog_id, sp2pids in sog_idx["sog_to_proteins"].items():
        sog_ref_members[sog_id] = set(sp2pids.get(ref_full, []))

    short_to_full = {sp.replace("Schizosaccharomyces_", "S_"): sp
                     for sp in sog_idx["species"]}

    logging.info(f"reading protein lengths {args.protein_fa}")
    protein_lengths = load_protein_lengths(args.protein_fa)
    logging.info(f"  {len(protein_lengths)} protein lengths loaded")

    counters: Dict[str, int] = defaultdict(int)

    # --- Phase 1: drop mRNAs overlapping L1 genes ---
    non_overlapping: List[Dict] = []
    n_overlap_dropped = 0
    for mid, m in mrnas.items():
        if overlaps_l1(m["seqid"], m["strand"], m["start"], m["end"],
                       l1_index, args.overlap_min):
            n_overlap_dropped += 1
        else:
            non_overlapping.append(m)
    logging.info(f"  L1-overlap dropped: {n_overlap_dropped}; non-overlapping candidates: {len(non_overlapping)}")

    # --- Phase 2: collapse into loci, classify ---
    loci = collapse_loci(non_overlapping)
    logging.info(f"  {len(loci)} loci after collapse")

    kept_data: List[Tuple[Dict, str, OrderedDict]] = []  # (best_hit, relation, tags)
    tsv_rows: List[List[str]] = []
    sidecar_rows: List[List[str]] = []

    for locus in loci:
        seqid = locus[0]["seqid"]
        strand = locus[0]["strand"]
        locus_start = min(h["start"] for h in locus)
        locus_end = max(h["end"] for h in locus)
        locus_coord = f"{seqid}:{locus_start}-{locus_end}({strand})"

        # Best SOG for this locus (lowest rank, then highest identity)
        og_m, og_m_source = None, None
        for h in sorted(locus, key=lambda x: (x["rank"], -x["identity"])):
            cand = sog_of_protein(h["src_protein"], sog_idx)
            if cand:
                og_m, og_m_source = cand, h["src_protein"]
                break

        adj_l = find_adjacent_l(seqid, strand, locus_start, locus_end,
                                l1_index, args.adjacent_max_bp)

        best_hit = max(locus, key=lambda x: (-x["rank"], x["identity"]))

        if og_m is None:
            counters["singleton_no_sog_to_sidecar"] += 1
            sidecar_rows.append([
                args.sample, locus_coord,
                best_hit["src_species"], best_hit["src_protein"],
                f"{best_hit['identity']:.4f}", "singleton_no_sog",
                f"adjacent_L={adj_l or '.'}",
            ])
            continue

        ref_present = bool(sog_ref_members.get(og_m, set()))
        relation, ev = classify_locus(locus, ref_present, args.id_adv_pp, args.diverged_max_id,
                                      og_m, sog_idx, short_to_full, ref_short)
        counters[relation] += 1

        # HGT identity gate
        if relation.startswith("intra_genus_HGT") and best_hit["identity"] < args.hgt_min_identity:
            counters["hgt_below_threshold_to_sidecar"] += 1
            sidecar_rows.append([
                args.sample, locus_coord,
                best_hit["src_species"], best_hit["src_protein"],
                f"{best_hit['identity']:.4f}", f"low_conf_{relation}",
                f"adjacent_L={adj_l or '.'}",
            ])
            continue

        ref_id_str = f"{ev['ref_id']:.4f}" if ev.get("ref_id") is not None else "NA"
        best_other_id_str = f"{ev['best_other_id']:.4f}" if ev.get("best_other_id") is not None else "NA"

        tags = OrderedDict([
            ("source", "miniprot_L2"),
            ("relation", relation),
            ("SOG_id", og_m),
            ("ref_id", ref_id_str),
            ("best_other_id", best_other_id_str),
            ("best_other_sp", ev.get("best_other_sp") or "."),
            ("adjacent_L", adj_l or "."),
        ])
        kept_data.append((best_hit, relation, tags))

        tsv_rows.append([
            args.sample, locus_coord, relation,
            og_m, og_m_source or ".",
            f"{best_hit['identity']:.4f}", str(best_hit["aln_aa"]),
            ref_id_str, best_other_id_str, ev.get("best_other_sp") or ".",
            adj_l or ".",
        ])

    logging.info(f"  classified {len(kept_data) + len(sidecar_rows)} loci with SOG")

    # --- Phase 3: ORF completeness + coverage filter ---
    gff_entries: List[Tuple[Dict, str, OrderedDict]] = []
    for m, relation, tags in kept_data:
        locus_coord = f"{m['seqid']}:{m['start']}-{m['end']}({m['strand']})"
        if not has_stop_codon(m):
            orf_status = "partial_no_stop" if m["aa_start"] > 1 else "full_aln_no_stop"
            counters[f"orf_{orf_status}_to_sidecar"] += 1
            sidecar_rows.append([
                args.sample, locus_coord,
                m["src_species"], m["src_protein"],
                f"{m['identity']:.4f}", orf_status,
                f"relation={relation}",
            ])
            continue
        aa_len = cds_aa_length(m)
        q_id = f"{m['src_species']}|{m['src_protein']}"
        q_len = protein_lengths.get(q_id, 0)
        cov = aa_len / q_len if q_len else 0.0
        if q_len and cov < args.orf_min_coverage:
            counters["orf_too_short_to_sidecar"] += 1
            sidecar_rows.append([
                args.sample, locus_coord,
                m["src_species"], m["src_protein"],
                f"{m['identity']:.4f}", "orf_too_short",
                f"relation={relation};orf_aa={aa_len};query_aa={q_len};cov={cov:.3f}",
            ])
            continue
        tags_with_orf = OrderedDict(tags)
        tags_with_orf["orf_status"] = "complete"
        gff_entries.append((m, relation, tags_with_orf))
        counters["orf_complete_kept"] += 1

    logging.info(f"  ORF filter: {counters['orf_complete_kept']} kept, "
                 f"{counters.get('orf_partial_no_stop_to_sidecar', 0)} partial_no_stop, "
                 f"{counters.get('orf_full_aln_no_stop_to_sidecar', 0)} full_aln_no_stop, "
                 f"{counters.get('orf_too_short_to_sidecar', 0)} too_short")

    # --- Write outputs ---
    for path in (args.out_gff, args.out_tsv, args.out_sidecar):
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)

    with open(args.out_gff, "w") as out:
        out.write("##gff-version 3\n")
        for i, (m, relation, tags) in enumerate(gff_entries):
            gene_id = f"L2G_{args.sample}_{i+1:05d}"
            mrna_id = f"L2T_{args.sample}_{i+1:05d}"
            lines = emit_gene_structure(m, gene_id, mrna_id, tags)
            out.write("\n".join(lines) + "\n")

    with open(args.out_tsv, "w") as out:
        out.write("\t".join([
            "sample", "locus", "relation", "SOG_id", "best_protein",
            "identity", "aln_aa", "ref_id", "best_other_id", "best_other_sp",
            "adjacent_L",
        ]) + "\n")
        for row in tsv_rows:
            out.write("\t".join(row) + "\n")

    with open(args.out_sidecar, "w") as out:
        out.write("\t".join([
            "sample", "locus", "src_species", "src_protein",
            "identity", "category", "note",
        ]) + "\n")
        for row in sidecar_rows:
            while len(row) < 7:
                row.append("")
            out.write("\t".join(row) + "\n")

    logging.info("classification breakdown:")
    for k in sorted(counters):
        logging.info(f"  {k:55s} {counters[k]:>8d}")
    logging.info(f"GFF entries (genes):   {len(gff_entries)}")
    logging.info(f"TSV rows:              {len(tsv_rows)}")
    logging.info(f"sidecar rows:          {len(sidecar_rows)}")


if __name__ == "__main__":
    main()
