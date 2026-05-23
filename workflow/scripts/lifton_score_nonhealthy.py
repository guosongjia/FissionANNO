#!/usr/bin/env python
"""Filter lifton score.txt to non-healthy transcripts.

Input is the score.txt produced by lifton (one row per transcript):
  tid  identity  ?  ?  ?  tool  flags  region

Output (sorted by gene_id, transcript_id):
  gene_id  transcript_id  identity  flags  region  classification

`flags` is the lifton column 7 verbatim (semicolon-separated). Each row is
classified by the strongest flag present:
  damaging > suspect

Healthy rows (any of identical/synonymous/nonsynonymous) are skipped
*unless* they also carry a damaging or suspect token (the lifton output
sometimes mixes synonymous with downgraded variants — keep those).

This sidecar is purely diagnostic: it does not feed downstream layers,
just preserves lifton's per-transcript flag info for QC use.
"""
import argparse
import logging
import os
import re
import sys


DAMAGING = {"frameshift", "start_lost", "stop_missing", "stop_codon_gain"}
SUSPECT = {"inframe_insertion", "inframe_deletion"}
HEALTHY = {"identical", "synonymous", "nonsynonymous"}


def base_gene(tid: str) -> str:
    parts = tid.rsplit(".", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    return tid


def classify(flag_set: set) -> str:
    if flag_set & DAMAGING:
        return "damaging"
    if flag_set & SUSPECT:
        return "suspect"
    return ""  # healthy / unrecognized


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--score", required=True, help="path to lifton score.txt")
    p.add_argument("--output", required=True)
    p.add_argument("--log-level", default="INFO")
    args = p.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO),
                        format="%(levelname)s: %(message)s")

    if not os.path.exists(args.score):
        logging.warning(f"score.txt missing: {args.score}; emitting empty output")
        os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
        with open(args.output, "w") as out:
            out.write("gene_id\ttranscript_id\tidentity\tflags\tregion\tclassification\n")
        return

    rows = []
    n_total = 0
    n_damaging = 0
    n_suspect = 0
    with open(args.score) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            n_total += 1
            tid, identity, _, _, _, _, flags, region = cols[:8]
            flag_set = set(flags.split(";"))
            cls = classify(flag_set)
            if not cls:
                continue
            if cls == "damaging":
                n_damaging += 1
            else:
                n_suspect += 1
            rows.append((base_gene(tid), tid, identity, flags, region, cls))

    rows.sort(key=lambda r: (r[0], r[1]))

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w") as out:
        out.write("gene_id\ttranscript_id\tidentity\tflags\tregion\tclassification\n")
        for r in rows:
            out.write("\t".join(r) + "\n")

    logging.info(f"total transcripts in score.txt: {n_total}")
    logging.info(f"damaging: {n_damaging} ({100*n_damaging/n_total:.1f}%)")
    logging.info(f"suspect: {n_suspect} ({100*n_suspect/n_total:.1f}%)")
    logging.info(f"wrote {len(rows)} rows to {args.output}")


if __name__ == "__main__":
    main()
