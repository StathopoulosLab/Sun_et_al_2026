#!/usr/bin/env python3
"""
step5b_merge_samples.py  —  Merge per-sample count matrices into one cross-sample matrix.

After step5 produces per-sample count matrices, this step combines them into
a single genes × embryos matrix across ALL samples.  This is essential for
cross-sample normalization (quantile norm, TPM) so that WT and mutant embryos
are on the same scale.

Output:
    results/combined/<prefix>_counts_matrix.tsv

Usage:
    python step5b_merge_samples.py \
        --samples WT_seRNAseq_4embryos_set1 Mutant_seRNAseq_19embryos_set2 \
        --prefix WT_vs_Mutant

    python step5b_merge_samples.py \
        --samples WT_seRNAseq_4embryos_set1 Mutant_seRNAseq_19embryos_set2 \
        --prefix WT_vs_Mutant \
        --annotate gene_id_to_symbol.tsv
"""

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import config

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--samples",  required=True, nargs="+",
                   help="Sample names (space-separated) to merge")
    p.add_argument("--prefix",   required=True,
                   help="Name for the combined output (e.g. WT_vs_Mutant)")
    p.add_argument("--annotate", default=None,
                   help="Optional: gene_id_to_symbol.tsv for gene symbol column")
    return p.parse_args()

def main():
    args = parse_args()

    combined_dir = os.path.join(config.OUTPUT_ROOT, "combined")
    os.makedirs(combined_dir, exist_ok=True)

    print(f"[step5b] Merging {len(args.samples)} samples into combined matrix")

    all_dfs = []
    embryo_sample_map = {}  # track which embryo belongs to which sample

    for sample in args.samples:
        matrix_path = os.path.join(
            config.OUTPUT_ROOT, sample, "counts", f"{sample}_counts_matrix.tsv")

        if not os.path.exists(matrix_path):
            sys.exit(f"[step5b] Count matrix not found: {matrix_path}\n"
                      f"         Run step5_merge_counts.py --sample {sample} first.")

        df = pd.read_csv(matrix_path, sep="\t", index_col=0)

        # Remove GeneSymbol if present (we'll re-add at the end)
        if "GeneSymbol" in df.columns:
            df = df.drop(columns="GeneSymbol")

        # Prefix embryo columns with sample name to avoid collisions
        # e.g. "embryo1" → "WT_set1__embryo1"
        # But also keep a clean version for display
        short_sample = sample  # could shorten further if needed
        new_cols = {}
        for col in df.columns:
            new_name = f"{short_sample}__{col}"
            new_cols[col] = new_name
            embryo_sample_map[new_name] = sample

        df = df.rename(columns=new_cols)
        all_dfs.append(df)

        print(f"[step5b]   {sample}: {df.shape[1]} embryos, {df.shape[0]:,} genes")

    # Merge on gene IDs (outer join keeps all genes)
    combined = pd.concat(all_dfs, axis=1, join="outer").fillna(0).astype(int)
    combined.index.name = "GeneID"

    # Optional: add gene symbols
    if args.annotate:
        sym = pd.read_csv(args.annotate, sep="\t", index_col="GeneID")
        combined.insert(0, "GeneSymbol",
                        combined.index.map(sym["Symbol"]).fillna(""))

    out_path = os.path.join(combined_dir, f"{args.prefix}_counts_matrix.tsv")
    combined.to_csv(out_path, sep="\t")

    # Save sample membership for downstream steps
    membership_path = os.path.join(combined_dir, f"{args.prefix}_sample_groups.tsv")
    with open(membership_path, "w") as fh:
        fh.write("Embryo\tSample\n")
        for embryo, sample in sorted(embryo_sample_map.items()):
            fh.write(f"{embryo}\t{sample}\n")

    embryo_cols = [c for c in combined.columns if c != "GeneSymbol"]
    n_embryos   = len(embryo_cols)

    print(f"\n[step5b] Combined matrix : {combined.shape[0]:,} genes × {n_embryos} embryos")
    print(f"[step5b] Saved → {out_path}")
    print(f"[step5b] Sample groups → {membership_path}")

    # Per-sample QC
    total_umis = combined[embryo_cols].sum()
    detected   = (combined[embryo_cols] > 0).sum()

    print(f"\n[step5b] Per-embryo QC:")
    print(f"  {'Embryo':<55} {'Sample':<35} {'UMIs':>10} {'Genes':>8}")
    for col in embryo_cols:
        print(f"  {col:<55} {embryo_sample_map[col]:<35} "
              f"{int(total_umis[col]):>10,} {int(detected[col]):>8,}")

    # Per-sample summary
    print(f"\n[step5b] Per-sample summary:")
    for sample in args.samples:
        sample_cols = [c for c in embryo_cols if embryo_sample_map[c] == sample]
        median_umis = total_umis[sample_cols].median()
        median_genes = detected[sample_cols].median()
        print(f"  {sample}: {len(sample_cols)} embryos, "
              f"median {int(median_umis):,} UMIs, "
              f"median {int(median_genes):,} genes detected")

if __name__ == "__main__":
    main()
