#!/usr/bin/env python
"""L1 vs L2 conflict resolution per CLAUDE.md §4.

Output:
  --out-gff       L2 hits kept (with provenance tags)
  --out-conflict  per-conflict decision log
  --out-sidecar   intra-genus HGT candidates (singletons / non-main hits)

NOT YET IMPLEMENTED — placeholder created during scaffold pass.
"""
import argparse
import sys


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sample", required=True)
    p.add_argument("--l1-gff", required=True)
    p.add_argument("--l2-gff", required=True)
    p.add_argument("--sog-index", required=True)
    p.add_argument("--overlap-min", type=float, required=True)
    p.add_argument("--id-adv-pp", type=float, required=True)
    p.add_argument("--bit-adv", type=float, required=True)
    p.add_argument("--diverged-max-id", type=float, required=True)
    p.add_argument("--adjacent-max-bp", type=int, required=True)
    p.add_argument("--out-gff", required=True)
    p.add_argument("--out-conflict", required=True)
    p.add_argument("--out-sidecar", required=True)
    args = p.parse_args()
    sys.stderr.write("l2_conflict_resolve.py: not implemented yet\n")
    sys.exit(2)


if __name__ == "__main__":
    main()
