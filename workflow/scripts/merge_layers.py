#!/usr/bin/env python
"""Merge L1 refined + L2 kept_rescued + L3 kept GFFs into one per-strain GFF.

All gene/mRNA features get source= tag per CLAUDE.md §7.4.
"""
import argparse
import sys


def read_gff_features(path, default_source=None):
    features = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if default_source and parts[2] in ("gene", "mRNA"):
                if "source=" not in parts[8]:
                    parts[8] = parts[8].rstrip(";") + f";source={default_source}"
            features.append(parts)
    return features


def group_by_gene(features):
    groups, current = [], []
    for parts in features:
        if parts[2] == "gene":
            if current:
                groups.append(current)
            current = [parts]
        else:
            current.append(parts)
    if current:
        groups.append(current)
    return groups


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sample", required=True)
    p.add_argument("--l1", required=True)
    p.add_argument("--l2", required=True)
    p.add_argument("--l3", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    g1 = group_by_gene(read_gff_features(args.l1, "lifton"))
    g2 = group_by_gene(read_gff_features(args.l2, "miniprot_L2"))
    g3 = group_by_gene(read_gff_features(args.l3, "annevo_L3"))

    all_groups = g1 + g2 + g3
    all_groups.sort(key=lambda g: (g[0][0], int(g[0][3])))

    with open(args.output, "w") as f:
        f.write("##gff-version 3\n")
        for group in all_groups:
            for parts in group:
                f.write("\t".join(parts) + "\n")

    print(
        f"{args.sample}: L1={len(g1)} L2={len(g2)} L3={len(g3)} total={len(all_groups)} genes",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
