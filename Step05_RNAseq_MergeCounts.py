#!/usr/bin/env python3
"""
step5_merge_counts.py  —  Merge per-embryo UMI count files into one matrix.

Reads all *_counts.txt files for a sample and concatenates them into a
genes × embryos count matrix saved to:
    results/<sample>/counts/<sample>_counts_matrix.tsv

Optionally annotate gene IDs with FlyBase gene symbols.

Usage:
    python step5_merge_counts.py --sample sample_4embryos
    python step5_merge_counts.py --sample sample_4embryos \
                                 --annotate gene_id_to_symbol.tsv
"""

import argparse
import glob
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import config

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--sample",   required=True)
    p.add_argument("--annotate", default=None,
                   help="Optional: path to gene_id_to_symbol.tsv (produced by "
                        "convert_flybase_mapping.py) to add a GeneSymbol column")
    return p.parse_args()

def main():
    args       = parse_args()
    counts_dir = os.path.join(config.OUTPUT_ROOT, args.sample, "counts")

    count_files = sorted(glob.glob(os.path.join(counts_dir, f"{args.sample}_*_counts.txt")))
    if not count_files:
        sys.exit(f"[step5] No count files found in {counts_dir}")

    print(f"[step5] Sample  : {args.sample}")
    print(f"[step5] Files   : {len(count_files)}")

    dfs = []
    for f in count_files:
        base   = os.path.basename(f)          # sample_4embryos_embryo1_counts.txt
        embryo = (base
                  .replace(f"{args.sample}_", "")
                  .replace("_counts.txt", ""))
        df = pd.read_csv(f, sep="\t", header=None, names=["GeneID", embryo])
        df.set_index("GeneID", inplace=True)
        dfs.append(df)

    matrix = pd.concat(dfs, axis=1, join="outer").fillna(0).astype(int)
    matrix.index.name = "GeneID"

    # Optional: add gene symbols
    if args.annotate:
        sym = pd.read_csv(args.annotate, sep="\t", index_col="GeneID")
        matrix.insert(0, "GeneSymbol", matrix.index.map(sym["Symbol"]).fillna(""))

    out_path = os.path.join(counts_dir, f"{args.sample}_counts_matrix.tsv")
    matrix.to_csv(out_path, sep="\t")

    print(f"[step5] Matrix shape : {matrix.shape[0]:,} genes × {len(dfs)} embryos")
    print(f"[step5] Saved → {out_path}")

    # Quick QC summary
    embryo_cols = [c for c in matrix.columns if c != "GeneSymbol"]
    total_umis   = matrix[embryo_cols].sum()
    detected     = (matrix[embryo_cols] > 0).sum()
    print("\n[step5] Per-embryo QC:")
    print(f"  {'Embryo':<25} {'Total UMIs':>12} {'Detected genes':>15}")
    for col in embryo_cols:
        print(f"  {col:<25} {int(total_umis[col]):>12,} {int(detected[col]):>15,}")

if __name__ == "__main__":
    main()
