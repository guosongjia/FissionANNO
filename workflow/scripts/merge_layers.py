#!/usr/bin/env python
"""Merge L1 refined + L2 kept + L3 kept GFFs into one per-strain GFF.

All output features carry source= tag per CLAUDE.md §7.4.
"""
import argparse
import sys


def read_gff_features(path, default_source=None):
    """Read GFF features, optionally adding source attribute."""
    features = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if default_source and parts[2] in ("gene", "mRNA"):
                if "source=" not in parts[8]:
                    parts[8] = parts[8].rstrip(";") + f";source={default_source}"
            features.append(parts)
    return features


def sort_key(parts):
    return (parts[0], int(parts[3]))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sample", required=True)
    p.add_argument("--l1", required=True)
    p.add_argument("--l2", required=True)
    p.add_argument("--l3", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    all_features = []
    all_features.extend(read_gff_features(args.l1, "lifton"))
    all_features.extend(read_gff_features(args.l2, "miniprot_L2"))
    all_features.extend(read_gff_features(args.l3, "annevo_L3"))

    gene_groups = []
    current_group = []
    for parts in all_features:
        if parts[2] == "gene":
            if current_group:
                gene_groups.append(current_group)
            current_group = [parts]
        else:
            current_group.append(parts)
    if current_group:
        gene_groups.append(current_group)

    gene_groups.sort(key=lambda g: sort_key(g[0]))

    with open(args.output, "w") as f:
        f.write("##gff-version 3\n")
        for group in gene_groups:
            for parts in group:
                f.write("\t".join(parts) + "\n")

    n_genes = sum(1 for g in gene_groups if g[0][2] == "gene")
    print(f"Merged {n_genes} genes for {args.sample}", file=sys.stderr)


if __name__ == "__main__":
    main()
