#!/usr/bin/env python
"""Build per-strain unmapped TSV from lifton output and refine stat.

Output schema:
  sample_name<TAB>gene_id<TAB>reason<TAB>parent_seq<TAB>note

Reason values:
  - lifton_unmapped              : lifton produced no transfer for the gene
                                   (gene_id appears in stats/unmapped_features.txt)
  - refine_pseudogene            : refine flagged the gene as pseudogene
                                   (refine stat: pseudogene=Y)
  - refine_truncated_at_contig_end : refine could not produce a valid ORF;
                                     the gene sits at a contig boundary
                                     (refine stat: valid_orf=N AND pseudogene=N AND truncated_orf=Y)
  - refine_frame_disrupted       : refine could not produce a valid ORF and
                                   the locus is not at a contig end — the
                                   reading frame is too disrupted for refine
                                   to repair, likely caused by sequence
                                   divergence (frameshifting indels / scattered
                                   substitutions that destroy stop codons)
                                   (refine stat: valid_orf=N AND pseudogene=N AND truncated_orf=N)

The first reason draws from `<lifton_outdir>/stats/unmapped_features.txt`
(falling back to `liftoff/unmapped_features.txt`). Each unmapped gene has
no parent_seq.

The other three reasons draw from refine stat (must contain the
truncated_orf column). parent_seq is the contig where the (failed) gene
landed in the refined GFF; note may carry extra context (e.g. lifton flag).
"""
import argparse
import logging
import os
import re
import sys
from typing import Dict, Tuple, Optional


def read_unmapped_set(lifton_outdir: str) -> set:
    candidates = [
        os.path.join(lifton_outdir, "stats", "unmapped_features.txt"),
        os.path.join(lifton_outdir, "liftoff", "unmapped_features.txt"),
    ]
    for p in candidates:
        if os.path.exists(p):
            unmapped = set()
            with open(p) as fh:
                for line in fh:
                    cols = line.rstrip("\n").split("\t")
                    if cols and cols[0]:
                        unmapped.add(cols[0])
            logging.info(f"read {len(unmapped)} unmapped gene IDs from {p}")
            return unmapped
    logging.warning(f"no unmapped_features.txt under {lifton_outdir}; assuming none")
    return set()


def read_refine_stat(path: str) -> Dict[str, Dict[str, str]]:
    rows: Dict[str, Dict[str, str]] = {}
    expected = ["gene_id", "cds_exon_match", "up_trigger", "up_len",
                "down_trigger", "down_len", "valid_orf", "pseudogene", "truncated_orf"]
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if header != expected:
            raise SystemExit(
                f"refine stat header mismatch in {path}\n"
                f"  expected: {expected}\n"
                f"  got:      {header}\n"
                f"This script requires the truncated_orf column. "
                f"Re-run lifton_gff3_refine.py with the current version."
            )
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < len(expected):
                continue
            rows[cols[0]] = dict(zip(header, cols))
    logging.info(f"read {len(rows)} stat rows from {path}")
    return rows


def gene_to_contig(gff_path: str) -> Dict[str, str]:
    g2c: Dict[str, str] = {}
    top_types = {"gene", "pseudogene"}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] not in top_types:
                continue
            m = re.search(r"ID=([^;]+)", cols[8])
            if not m:
                continue
            g2c[m.group(1)] = cols[0]
    return g2c


def read_score_flags(score_path: str) -> Dict[str, str]:
    """gene_id -> joined flags (from any of its transcripts) for note column."""
    if not os.path.exists(score_path):
        return {}

    def base_gene(tid: str) -> str:
        parts = tid.rsplit(".", 1)
        if len(parts) == 2 and parts[1].isdigit():
            return parts[0]
        return tid

    flags: Dict[str, set] = {}
    with open(score_path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 7:
                continue
            g = base_gene(cols[0])
            flags.setdefault(g, set()).update(cols[6].split(";"))
    return {g: ";".join(sorted(f)) for g, f in flags.items()}


def classify(stat_row: Dict[str, str]) -> Optional[str]:
    if stat_row["pseudogene"] == "Y":
        return "refine_pseudogene"
    if stat_row["valid_orf"] == "N":
        if stat_row["truncated_orf"] == "Y":
            return "refine_truncated_at_contig_end"
        return "refine_frame_disrupted"
    return None  # valid_orf=Y AND pseudogene=N -> not unmapped


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample", required=True)
    p.add_argument("--lifton-outdir", required=True,
                   help="path to <sample>.lifton_output (containing stats/, liftoff/, score.txt)")
    p.add_argument("--refine-stat", required=True)
    p.add_argument("--refined-gff", required=True,
                   help="refine GFF; used to look up parent_seq for non-unmapped reasons")
    p.add_argument("--output", required=True)
    p.add_argument("--log-level", default="INFO")
    args = p.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO),
                        format="%(levelname)s: %(message)s")

    unmapped = read_unmapped_set(args.lifton_outdir)
    stats = read_refine_stat(args.refine_stat)
    g2c = gene_to_contig(args.refined_gff)
    score_flags = read_score_flags(os.path.join(args.lifton_outdir, "score.txt"))

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    n_per_reason: Dict[str, int] = {}
    with open(args.output, "w") as out:
        out.write("sample_name\tgene_id\treason\tparent_seq\tnote\n")
        for gid in sorted(unmapped):
            n_per_reason["lifton_unmapped"] = n_per_reason.get("lifton_unmapped", 0) + 1
            out.write(f"{args.sample}\t{gid}\tlifton_unmapped\t.\t.\n")
        for gid, row in stats.items():
            reason = classify(row)
            if reason is None:
                continue
            parent = g2c.get(gid, ".")
            note_bits = []
            if gid in score_flags:
                note_bits.append(f"lifton_flags={score_flags[gid]}")
            note = ";".join(note_bits) if note_bits else "."
            n_per_reason[reason] = n_per_reason.get(reason, 0) + 1
            out.write(f"{args.sample}\t{gid}\t{reason}\t{parent}\t{note}\n")

    for r, n in sorted(n_per_reason.items()):
        logging.info(f"{args.sample}\t{r}\t{n}")


if __name__ == "__main__":
    main()
