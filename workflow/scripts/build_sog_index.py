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
A pickle file containing a dict with three keys:
  - "protein_to_sog": {protein_id: SOG_id}
  - "sog_to_proteins": {SOG_id: {species_full_name: [protein_id, ...]}}
  - "sog_to_category": {SOG_id: OG_category}
  - "species": [species_full_name, ...]   # column order preserved
  - "source_path": str
  - "n_sogs": int
  - "n_proteins": int

A warning is emitted for any protein_id that appears in more than one SOG
(should not happen with the curated table, but worth catching).
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


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sog", required=True, help="path to SOG table TSV")
    p.add_argument("--output", required=True, help="output pickle path")
    p.add_argument("--no-strict", action="store_true",
                   help="proceed even if rows have wrong column count (default: fail loud)")
    p.add_argument("--log-level", default="INFO")
    args = p.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO),
                        format="%(levelname)s: %(message)s")
    logging.info(f"reading {args.sog}")
    idx = parse_sog_table(args.sog, strict=not args.no_strict)
    logging.info(f"parsed {idx['n_sogs']} SOGs, {idx['n_proteins']} proteins, {len(idx['species'])} species")

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "wb") as fh:
        pickle.dump(idx, fh, protocol=pickle.HIGHEST_PROTOCOL)
    logging.info(f"wrote {args.output}")


if __name__ == "__main__":
    main()
