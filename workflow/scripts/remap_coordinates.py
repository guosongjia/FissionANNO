#!/usr/bin/env python
"""Remap ANNEVO prediction GFF coordinates from residual fragments back to
original genome coordinates using the mapping TSV from extract_residuals.py.
"""
import argparse
import sys


def load_mapping(map_path):
    """Load residual_id → (orig_seqid, orig_start) mapping."""
    mapping = {}
    with open(map_path) as f:
        next(f)  # skip header
        for line in f:
            rid, seqid, start, end = line.rstrip("\n").split("\t")
            mapping[rid] = (seqid, int(start))
    return mapping


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gff", required=True)
    p.add_argument("--map", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    mapping = load_mapping(args.map)
    n_features = 0

    with open(args.gff) as fin, open(args.output, "w") as fout:
        fout.write("##gff-version 3\n")
        for line in fin:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            rid = parts[0]
            if rid not in mapping:
                continue
            orig_seqid, offset = mapping[rid]
            parts[0] = orig_seqid
            parts[3] = str(int(parts[3]) + offset - 1)
            parts[4] = str(int(parts[4]) + offset - 1)
            fout.write("\t".join(parts) + "\n")
            n_features += 1

    print(f"Remapped {n_features} features", file=sys.stderr)


if __name__ == "__main__":
    main()
