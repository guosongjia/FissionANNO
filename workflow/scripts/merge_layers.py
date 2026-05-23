#!/usr/bin/env python
"""Merge L1 refined + L2 kept + L3 kept GFFs into one per-strain GFF.

All output features carry source= / SOG_id= / HGT_call= / relation= tags
per CLAUDE.md §7.4.

NOT YET IMPLEMENTED — placeholder created during scaffold pass.
"""
import argparse
import sys


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sample", required=True)
    p.add_argument("--l1", required=True)
    p.add_argument("--l2", required=True)
    p.add_argument("--l3", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()
    sys.stderr.write("merge_layers.py: not implemented yet\n")
    sys.exit(2)


if __name__ == "__main__":
    main()
