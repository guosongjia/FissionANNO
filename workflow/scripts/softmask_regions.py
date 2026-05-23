#!/usr/bin/env python
"""Soft-mask regions covered by L1+L2 GFFs (lowercase the masked bases).

Used as input to BRAKER4 so it has HMM context but skips already-annotated
loci.

NOT YET IMPLEMENTED — placeholder created during scaffold pass.
"""
import argparse
import sys


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fasta", required=True)
    p.add_argument("--gffs", nargs="+", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()
    sys.stderr.write("softmask_regions.py: not implemented yet\n")
    sys.exit(2)


if __name__ == "__main__":
    main()
