#!/usr/bin/env python
"""Extract residual (non-N) intervals ≥ min_bp from a hard-masked genome.

Outputs:
  - residual FASTA with entries named {seqid}:{start}-{end} (1-based)
  - coordinate mapping TSV for remapping predictions back to original coords
"""
import argparse
import sys


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fasta", required=True)
    p.add_argument("--min-bp", type=int, default=1000)
    p.add_argument("--out-fasta", required=True)
    p.add_argument("--out-map", required=True)
    args = p.parse_args()

    n_intervals = 0
    total_bp = 0

    with open(args.fasta) as fin, \
         open(args.out_fasta, "w") as fout, \
         open(args.out_map, "w") as fmap:
        fmap.write("residual_id\torig_seqid\torig_start\torig_end\n")

        seqid = None
        seq_parts = []

        def flush():
            nonlocal n_intervals, total_bp
            if seqid is None:
                return
            seq = "".join(seq_parts)
            i = 0
            while i < len(seq):
                if seq[i] == "N":
                    i += 1
                    continue
                j = i
                while j < len(seq) and seq[j] != "N":
                    j += 1
                length = j - i
                if length >= args.min_bp:
                    start_1based = i + 1
                    end_1based = j
                    rid = f"{seqid}:{start_1based}-{end_1based}"
                    fout.write(f">{rid}\n")
                    segment = seq[i:j]
                    for k in range(0, len(segment), 80):
                        fout.write(segment[k:k+80] + "\n")
                    fmap.write(f"{rid}\t{seqid}\t{start_1based}\t{end_1based}\n")
                    n_intervals += 1
                    total_bp += length
                i = j

        for line in fin:
            if line.startswith(">"):
                flush()
                seqid = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        flush()

    print(f"Extracted {n_intervals} residual intervals ({total_bp:,} bp total)", file=sys.stderr)


if __name__ == "__main__":
    main()
