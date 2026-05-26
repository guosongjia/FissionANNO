#!/usr/bin/env python
"""Hard-mask regions covered by L1+L2 GFFs (replace with N).

Produces a genome where all annotated gene regions are N-masked,
leaving only residual (unannotated) regions as real sequence.
"""
import argparse
import sys
from collections import defaultdict


def parse_gene_intervals(gff_paths):
    intervals = defaultdict(list)
    for path in gff_paths:
        with open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                if parts[2] == "gene":
                    intervals[parts[0]].append((int(parts[3]) - 1, int(parts[4])))
    return intervals


def merge_intervals(ivs):
    if not ivs:
        return []
    ivs.sort()
    merged = [ivs[0]]
    for s, e in ivs[1:]:
        if s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return merged


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fasta", required=True)
    p.add_argument("--gffs", nargs="+", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    intervals = parse_gene_intervals(args.gffs)
    for seqid in intervals:
        intervals[seqid] = merge_intervals(intervals[seqid])

    with open(args.fasta) as fin, open(args.output, "w") as fout:
        seqid = None
        seq_parts = []

        def flush():
            if seqid is None:
                return
            seq = "".join(seq_parts).upper()
            seq_arr = list(seq)
            for start, end in intervals.get(seqid, []):
                for i in range(start, min(end, len(seq_arr))):
                    seq_arr[i] = "N"
            masked = "".join(seq_arr)
            fout.write(f">{seqid}\n")
            for i in range(0, len(masked), 80):
                fout.write(masked[i:i+80] + "\n")

        for line in fin:
            if line.startswith(">"):
                flush()
                seqid = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        flush()

    n_intervals = sum(len(v) for v in intervals.values())
    print(f"Hard-masked {n_intervals} gene regions", file=sys.stderr)


if __name__ == "__main__":
    main()
