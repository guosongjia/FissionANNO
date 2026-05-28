#!/usr/bin/env python
"""Build protein_id -> SOG_id index from the curated SOG table.

Input
-----
A tab-separated file with header row. Columns:
  SOG_id  Number  OG_category  <species_1> ... <species_N>  of_og_num

Each species cell is empty or a comma-space-separated list of protein IDs.
The species columns are auto-detected as everything between the third
column (OG_category) and the last column (of_og_num).

Output
------
A pickle file containing a dict with keys:
  - "protein_to_sog": {protein_id: SOG_id}
  - "sog_to_proteins": {SOG_id: {species_full_name: [protein_id, ...]}}
  - "sog_to_category": {SOG_id: OG_category}
  - "sog_ref_max_id": {SOG_id: {other_sp_short: max_identity_float}}
  - "species": [species_full_name, ...]   # column order preserved
  - "source_path": str
  - "n_sogs": int
  - "n_proteins": int

When --protein-fa is provided, pairwise global alignments are computed for
each SOG between reference species (ref) members and each other species' members. The maximum
identity (matches / aligned_length incl. gaps) per (SOG, other_species) pair
is stored in "sog_ref_max_id".
"""
import argparse
import logging
import os
import pickle
import sys
from typing import Dict, List


def parse_sog_table(path: str, strict: bool = True):
    species: List[str] = []
    protein_to_sog: Dict[str, str] = {}
    sog_to_proteins: Dict[str, Dict[str, List[str]]] = {}
    sog_to_category: Dict[str, str] = {}
    duplicates: Dict[str, List[str]] = {}
    malformed: List[int] = []

    with open(path, "r") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if len(header) < 5 or header[0] != "SOG_id":
            raise ValueError(f"Unexpected header: {header[:5]}")
        # species columns: between OG_category (idx 2) and last column (of_og_num)
        species_idx_range = list(range(3, len(header) - 1))
        species = [header[i] for i in species_idx_range]

        for ln, line in enumerate(fh, start=2):
            cols = line.rstrip("\n").split("\t")
            if len(cols) != len(header):
                logging.warning(f"line {ln}: column count {len(cols)} != header {len(header)}; "
                                f"first cols={cols[:3]}")
                malformed.append(ln)
                continue
            sog_id = cols[0]
            category = cols[2]
            sog_to_category[sog_id] = category
            sog_to_proteins[sog_id] = {sp: [] for sp in species}
            for sp_name, idx in zip(species, species_idx_range):
                cell = cols[idx].strip()
                if not cell:
                    continue
                pids = [p.strip() for p in cell.split(",") if p.strip()]
                sog_to_proteins[sog_id][sp_name] = pids
                for pid in pids:
                    if pid in protein_to_sog and protein_to_sog[pid] != sog_id:
                        duplicates.setdefault(pid, [protein_to_sog[pid]]).append(sog_id)
                    protein_to_sog[pid] = sog_id

    if duplicates:
        logging.warning(f"{len(duplicates)} protein IDs appear in multiple SOGs (using last seen)")
        for pid, sogs in list(duplicates.items())[:5]:
            logging.warning(f"  {pid} -> {sogs}")

    if malformed and strict:
        raise SystemExit(
            f"build_sog_index: {len(malformed)} malformed rows in {path}; "
            f"first offending line numbers: {malformed[:5]}. "
            f"Fix the source TSV or pass --no-strict to ignore."
        )

    return {
        "protein_to_sog": protein_to_sog,
        "sog_to_proteins": sog_to_proteins,
        "sog_to_category": sog_to_category,
        "species": species,
        "source_path": os.path.abspath(path),
        "n_sogs": len(sog_to_proteins),
        "n_proteins": len(protein_to_sog),
        "malformed_lines": malformed,
    }


def load_protein_seqs(fa_path: str) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    cur_id, cur_parts = None, []
    with open(fa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(cur_parts)
                cur_id = line[1:].split()[0]
                cur_parts = []
            else:
                cur_parts.append(line.strip())
        if cur_id is not None:
            seqs[cur_id] = "".join(cur_parts)
    return seqs


def compute_ref_pairwise(sog_to_proteins: Dict, species: List[str],
                          fa_path: str, ref_full: str) -> Dict[str, Dict[str, float]]:
    """For each SOG, compute max identity (matches / aligned_length incl gaps)
    between every reference species protein and every other-species protein. Identity definition
    matches miniprot's Identity field for cross-comparability.

    Returns {SOG_id: {other_sp_short: max_identity}}.
    """
    from Bio.Align import PairwiseAligner, substitution_matrices

    full_to_short = {sp: sp.replace("Schizosaccharomyces_", "S_") for sp in species}
    ref_short = ref_full.replace("Schizosaccharomyces_", "S_")
    seqs = load_protein_seqs(fa_path)
    logging.info(f"  loaded {len(seqs)} protein sequences from {fa_path}")

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1

    out: Dict[str, Dict[str, float]] = {}
    n_sogs_with_pair, n_alignments = 0, 0
    for sog_id, sp2pids in sog_to_proteins.items():
        ref_pids = sp2pids.get(ref_full, [])
        if not ref_pids:
            continue
        per_other: Dict[str, float] = {}
        for sp_full, pids in sp2pids.items():
            if sp_full == ref_full or not pids:
                continue
            sp_short = full_to_short[sp_full]
            best = 0.0
            for a in ref_pids:
                sa = seqs.get(f"{ref_short}|{a}")
                if not sa:
                    continue
                for b in pids:
                    sb = seqs.get(f"{sp_short}|{b}")
                    if not sb:
                        continue
                    aln = aligner.align(sa, sb)[0]
                    aligned_a, aligned_b = str(aln[0]), str(aln[1])
                    aln_len = len(aligned_a)
                    if aln_len == 0:
                        continue
                    matches = sum(1 for x, y in zip(aligned_a, aligned_b)
                                  if x == y and x != "-")
                    ident = matches / aln_len
                    if ident > best:
                        best = ident
                    n_alignments += 1
            if best > 0:
                per_other[sp_short] = best
        if per_other:
            out[sog_id] = per_other
            n_sogs_with_pair += 1
    logging.info(f"  computed pairwise identities for {n_sogs_with_pair} SOGs "
                 f"({n_alignments} alignments)")
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sog", required=True, help="path to SOG table TSV")
    p.add_argument("--protein-fa", default=None,
                   help="combined 9-species protein fasta for pairwise identity computation")
    p.add_argument("--output", required=True, help="output pickle path")
    p.add_argument("--no-strict", action="store_true",
                   help="proceed even if rows have wrong column count (default: fail loud)")
    p.add_argument("--ref-species", default="Schizosaccharomyces_pombe",
                   help="full species name of the reference species in the SOG table")
    p.add_argument("--log-level", default="INFO")
    args = p.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO),
                        format="%(levelname)s: %(message)s")
    logging.info(f"reading {args.sog}")
    idx = parse_sog_table(args.sog, strict=not args.no_strict)
    logging.info(f"parsed {idx['n_sogs']} SOGs, {idx['n_proteins']} proteins, {len(idx['species'])} species")

    if args.protein_fa:
        logging.info(f"computing reference pairwise identities from {args.protein_fa}")
        idx["sog_ref_max_id"] = compute_ref_pairwise(
            idx["sog_to_proteins"], idx["species"], args.protein_fa, args.ref_species)
    else:
        idx["sog_ref_max_id"] = {}

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "wb") as fh:
        pickle.dump(idx, fh, protocol=pickle.HIGHEST_PROTOCOL)
    logging.info(f"wrote {args.output}")


if __name__ == "__main__":
    main()
