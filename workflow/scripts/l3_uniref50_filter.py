#!/usr/bin/env python
"""Filter BRAKER4 predictions by UniRef50 DIAMOND hits.

Rules (CLAUDE.md §3.4):
  - no UniRef50 hit -> drop
  - any hit -> keep
  - top hit not Schizosaccharomyces -> tag HGT_call=putative_<top_taxon>

NOT YET IMPLEMENTED — placeholder created during scaffold pass.
"""
import argparse
import sys


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--braker-gff", required=True)
    p.add_argument("--genome", required=True)
    p.add_argument("--diamond-db", required=True)
    p.add_argument("--evalue", type=float, required=True)
    p.add_argument("--schizo-keyword", required=True)
    p.add_argument("--threads", type=int, required=True)
    p.add_argument("--out-gff", required=True)
    p.add_argument("--out-diamond", required=True)
    args = p.parse_args()
    sys.stderr.write("l3_uniref50_filter.py: not implemented yet\n")
    sys.exit(2)


if __name__ == "__main__":
    main()
