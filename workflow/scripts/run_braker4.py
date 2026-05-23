#!/usr/bin/env python
"""Drive the external BRAKER4 Snakemake workflow for a single sample.

Strategy: write a per-sample CSV into the BRAKER4 workflow's input slot,
invoke `bash run_snakemake.sh` (or equivalent), then locate the output GFF
and copy it to --output-gff.

NOT YET IMPLEMENTED — placeholder created during scaffold pass.
"""
import argparse
import sys


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sample", required=True)
    p.add_argument("--masked-fasta", required=True)
    p.add_argument("--protein", required=True)
    p.add_argument("--busco-lineage", required=True)
    p.add_argument("--braker4-workflow", required=True)
    p.add_argument("--threads", type=int, required=True)
    p.add_argument("--output-gff", required=True)
    args = p.parse_args()
    sys.stderr.write("run_braker4.py: not implemented yet\n")
    sys.exit(2)


if __name__ == "__main__":
    main()
