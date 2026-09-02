#!/usr/bin/env python3
"""
step6_tpm.py  —  Normalize UMI counts to TPM (+ optional quantile norm).

Changes vs previous version
----------------------------
* Added --audit-goi: after TPM calculation, prints per-embryo log2(TPM+1)
  for a configurable set of key diagnostic genes so you can immediately
  spot genes that are being lost by the gene-length intersection step.
* Added --min-gene-expr-frac: drop genes where fewer than this fraction of
  embryos have TPM > 0 after normalization. Default 0 (no filter) — set to
  e.g. 0.05 to remove extreme singletons that inflate quantile normalization.
* Quantile normalization now also writes a log-scale QN matrix
  (<label>_log2tpm_quantile_matrix.tsv) for tools that expect log-space input
  directly (avoids repeated log transformation in downstream scripts).
* Per-embryo QC table now shows group/stage from sample_groups when available.
* Warns explicitly when a key gene (run, tll, eve, sna, HLH54F) is absent
  from the gene-lengths file and will be dropped from the TPM matrix.

Usage
-----
  python step6_tpm.py --combined all_genotypes --qc --quantile
  python step6_tpm.py --combined all_genotypes --qc --quantile --audit-goi
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import config

try:
    import Config_SampleMetadata as SC
    _HAS_SC = True
except ImportError:
    _HAS_SC = False

# Genes that are biologically important and should never silently disappear
_SENTINEL_GENES = [
    "run", "eve", "H", "ftz", "tll", "hkb", "hb", "Kr", "kni",
    "sna", "twi", "HLH54F", "FoxL1", "bcd", "cic", "sog",
]

# GOI for --audit-goi (symbols — resolved from gene_lengths index or symbol col)
_AUDIT_GENES = [
    "run", "eve", "H", "tll", "hkb", "hb", "Kr", "kni", "gt",
    "sna", "twi", "HLH54F", "FoxL1", "bcd", "cic",
]


def parse_args():
    p = argparse.ArgumentParser()
    grp = p.add_mutually_exclusive_group(required=True)
    grp.add_argument("--sample",   help="Per-sample mode")
    grp.add_argument("--combined", help="Combined mode: prefix from step5b")
    p.add_argument("--gene-lengths", default=config.GENE_LENGTHS)
    p.add_argument("--symbol-map", default="gene_id_to_symbol.tsv",
                   help="Fallback Symbol/GeneID map, used when the input matrix "
                        "has no inline GeneSymbol column (e.g. RUVg-corrected "
                        "matrices, which drop it) — needed for the sentinel-gene "
                        "check and --audit-goi to resolve symbols correctly. "
                        "Default: gene_id_to_symbol.tsv")
    p.add_argument("--quantile",     action="store_true",
                   help="Also output quantile-normalised TPM (and log2-QN)")
    p.add_argument("--qc",           action="store_true",
                   help="Use QC-filtered matrix from step5c")
    p.add_argument("--min-gene-expr-frac", type=float, default=0.0,
                   help="Drop genes expressed (TPM>0) in fewer than this fraction "
                        "of embryos (default: 0 = keep all)")
    p.add_argument("--audit-goi",    action="store_true",
                   help="Print per-embryo log2(TPM+1) for key diagnostic genes")
    return p.parse_args()


def resolve_paths(args):
    if args.combined:
        base   = os.path.join(config.OUTPUT_ROOT, "combined")
        suffix = "_counts_matrix_qc.tsv" if args.qc else "_counts_matrix.tsv"
        label  = f"{args.combined}_qc"   if args.qc else args.combined
        return os.path.join(base, f"{args.combined}{suffix}"), base, label
    else:
        base = os.path.join(config.OUTPUT_ROOT, args.sample, "counts")
        return (os.path.join(base, f"{args.sample}_counts_matrix.tsv"),
                base, args.sample)


def counts_to_tpm(counts: pd.DataFrame, lengths_kb: pd.Series) -> pd.DataFrame:
    rpk   = counts.div(lengths_kb, axis=0)
    scale = rpk.sum(axis=0) / 1_000_000
    return rpk.div(scale, axis=1)


def quantile_normalize(df: pd.DataFrame) -> pd.DataFrame:
    """limma-style quantile normalization."""
    vals   = df.values.astype(float)
    sorted_vals = np.sort(vals, axis=0)
    ref    = sorted_vals.mean(axis=1)
    result = np.empty_like(vals)
    for j in range(vals.shape[1]):
        order = np.argsort(vals[:, j], kind="mergesort")
        result[order, j] = ref
    return pd.DataFrame(result, index=df.index, columns=df.columns)


def _load_groups(combined_dir: str, prefix: str) -> "pd.DataFrame | None":
    for fname in (f"{prefix}_sample_groups_qc.tsv",
                  f"{prefix}_sample_groups.tsv"):
        path = os.path.join(combined_dir, fname)
        if os.path.exists(path):
            df = pd.read_csv(path, sep="\t")
            # Normalise: ensure "Embryo" is the index
            if "Embryo" in df.columns:
                return df.set_index("Embryo")
    return None


def main():
    args = parse_args()
    counts_path, out_dir, label = resolve_paths(args)

    if not os.path.exists(counts_path):
        sys.exit(f"[step6] Not found: {counts_path}\n"
                  "        Run step5b / step5c first.")
    os.makedirs(out_dir, exist_ok=True)

    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    gene_symbols = counts.pop("GeneSymbol") if "GeneSymbol" in counts.columns else None

    # Fallback: matrices with no inline GeneSymbol column (e.g. RUVg-corrected
    # output — ruvg_nextgen.R strips GeneSymbol before processing) still need
    # symbol resolution for the sentinel-gene check and --audit-goi. Without
    # this, both would silently check literal symbol strings ("run", "eve",
    # ...) against an FBgn-indexed file/matrix, which can never match,
    # producing false "gene missing" warnings even when the gene is present.
    if gene_symbols is None and os.path.exists(args.symbol_map):
        sym_df = pd.read_csv(args.symbol_map, sep="\t").dropna(subset=["Symbol", "GeneID"])
        sym_df = sym_df[sym_df["GeneID"].isin(counts.index)]
        if len(sym_df):
            gene_symbols = pd.Series(sym_df["Symbol"].values, index=sym_df["GeneID"].values)
            # Reindex to cover every gene in the matrix (blank for genes with
            # no symbol entry), matching how the inline GeneSymbol column
            # always behaved — avoids a KeyError later at tpm.insert(...,
            # gene_symbols.loc[tpm.index]) for any gene the map doesn't cover.
            gene_symbols = gene_symbols.reindex(counts.index).fillna("")
            print(f"[step6] No inline GeneSymbol column — resolved "
                  f"{len(sym_df):,} symbols from {args.symbol_map}")
    elif gene_symbols is None:
        print(f"[step6] No inline GeneSymbol column and {args.symbol_map} not "
              f"found — sentinel-gene check and --audit-goi symbol resolution "
              f"will be unreliable.")

    gene_lengths = pd.read_csv(args.gene_lengths, sep="\t", index_col=0)

    # ── Sentinel gene check ───────────────────────────────────────────────────
    # Build a symbol → GeneID map from gene_symbols if available
    sym_to_id: dict[str, str] = {}
    if gene_symbols is not None:
        sym_to_id = dict(zip(gene_symbols.values, gene_symbols.index))

    for sym in _SENTINEL_GENES:
        gid = sym_to_id.get(sym, sym)
        if gid not in gene_lengths.index and sym not in gene_lengths.index:
            print(f"[step6] WARNING: sentinel gene '{sym}' (id: '{gid}') "
                  f"absent from gene-lengths — will be dropped from TPM output. "
                  f"Check {args.gene_lengths} is complete.")

    # ── Intersect with lengths ────────────────────────────────────────────────
    common  = counts.index.intersection(gene_lengths.index)
    dropped = len(counts) - len(common)
    if dropped:
        print(f"[step6] {dropped} genes without length entry dropped from TPM")

    counts     = counts.loc[common]
    lengths_kb = gene_lengths.loc[common, "Length"] / 1000

    # ── TPM ───────────────────────────────────────────────────────────────────
    tpm          = counts_to_tpm(counts, lengths_kb)
    embryo_cols  = list(tpm.columns)

    # ── Optional sparsity filter ──────────────────────────────────────────────
    if args.min_gene_expr_frac > 0:
        expr_frac = (tpm > 0).mean(axis=1)
        keep_genes = expr_frac >= args.min_gene_expr_frac
        n_dropped  = (~keep_genes).sum()
        if n_dropped:
            print(f"[step6] Sparsity filter (frac>{args.min_gene_expr_frac}): "
                  f"dropping {n_dropped} genes, keeping {keep_genes.sum():,}")
            tpm    = tpm.loc[keep_genes]
            counts = counts.loc[keep_genes]

    expressed = (tpm > 1).any(axis=1).sum()

    if gene_symbols is not None:
        tpm.insert(0, "GeneSymbol", gene_symbols.loc[tpm.index])

    out_path = os.path.join(out_dir, f"{label}_tpm_matrix.tsv")
    tpm.to_csv(out_path, sep="\t")

    print(f"[step6] Mode     : {'combined' if args.combined else 'per-sample'}")
    print(f"[step6] Label    : {label}")
    print(f"[step6] Embryos  : {len(embryo_cols)}")
    print(f"[step6] Genes    : {len(tpm):,}  (TPM>1 in ≥1 embryo: {expressed:,})")
    print(f"[step6] Saved → {out_path}")

    # ── Load group metadata for reporting ─────────────────────────────────────
    groups_df = None
    if args.combined:
        groups_df = _load_groups(out_dir, args.combined)

    # ── GOI audit ─────────────────────────────────────────────────────────────
    if args.audit_goi:
        tpm_numeric = tpm[embryo_cols]
        audit_ids: dict[str, str] = {}
        if gene_symbols is not None:
            sym_map = dict(zip(gene_symbols.values, gene_symbols.index))
            for sym in _AUDIT_GENES:
                gid = sym_map.get(sym, sym)
                if gid in tpm_numeric.index:
                    audit_ids[sym] = gid
                else:
                    print(f"[step6] audit: '{sym}' not found in TPM matrix")
        else:
            for sym in _AUDIT_GENES:
                if sym in tpm_numeric.index:
                    audit_ids[sym] = sym

        if audit_ids:
            print(f"\n[step6] GOI audit — log2(TPM+1):")
            header = f"  {'Gene':<12}" + "".join(
                f"{c.split('__')[-1][:12]:>14}" for c in embryo_cols[:20])
            print(header)
            for sym, gid in audit_ids.items():
                vals = np.log2(tpm_numeric.loc[gid].values[:20] + 1)
                row  = f"  {sym:<12}" + "".join(f"{v:>14.2f}" for v in vals)
                print(row)
            if len(embryo_cols) > 20:
                print(f"  (showing first 20 of {len(embryo_cols)} embryos)")

    # ── Quantile-normalized TPM ───────────────────────────────────────────────
    if args.quantile:
        tpm_numeric = tpm[embryo_cols].copy()

        # QN in log2 space for variance stabilization
        log_tpm  = np.log2(tpm_numeric + 1)
        log_qn   = quantile_normalize(log_tpm)
        tpm_qn   = (2 ** log_qn) - 1
        tpm_qn   = tpm_qn.clip(lower=0)

        # ── Raw QN TPM ────────────────────────────────────────────────────────
        if gene_symbols is not None:
            tpm_qn.insert(0, "GeneSymbol", gene_symbols.loc[tpm_qn.index])
        qn_path = os.path.join(out_dir, f"{label}_tpm_quantile_matrix.tsv")
        tpm_qn.to_csv(qn_path, sep="\t")

        # ── Log2-QN matrix (convenience for downstream scripts) ───────────────
        log_qn_out = log_qn.copy()
        if gene_symbols is not None:
            log_qn_out.insert(0, "GeneSymbol", gene_symbols.loc[log_qn_out.index])
        log_qn_path = os.path.join(out_dir, f"{label}_log2tpm_quantile_matrix.tsv")
        log_qn_out.to_csv(log_qn_path, sep="\t")

        # Write R-safe version without GeneSymbol (read.table chokes on quoted symbols)
        log_qn_nolabel = log_qn_out.drop(columns="GeneSymbol", errors="ignore")
        log_qn_nolabel.to_csv(log_qn_path.replace(".tsv", "_nolabel.tsv"), sep="\t")

        expressed_qn = (tpm_qn[embryo_cols] > 1).any(axis=1).sum()
        print(f"\n[step6] QN-TPM genes    : {len(tpm_qn):,}  (TPM>1: {expressed_qn:,})")
        print(f"[step6] Saved → {qn_path}")
        print(f"[step6] Saved → {log_qn_path}")

        # Per-embryo median comparison
        print(f"\n[step6] Per-embryo median TPM (detected genes, TPM>0):")
        has_grp = groups_df is not None and "group" in groups_df.columns
        hdr_cols = "  " + ("" if not has_grp else f"{'Group':<16} {'Stage':<8}  ")
        print(f"  {'Embryo':<52} {'Group':<16} {'Stage':<8} "
              f"{'Raw med':>10} {'QN med':>10}")
        for col in embryo_cols:
            bare = col.split("__")[-1]
            raw_med = float(tpm_numeric[col][tpm_numeric[col] > 0].median())
            qn_col  = tpm_qn[col] if isinstance(tpm_qn.columns[0], str) else tpm_qn.iloc[:, embryo_cols.index(col)]
            qn_med  = float(qn_col[qn_col > 0].median()) if col in tpm_qn.columns else 0.0
            grp  = str(groups_df.loc[col, "group"]) if has_grp and col in groups_df.index else ""
            stg  = str(groups_df.loc[col, "stage"]) if has_grp and col in groups_df.index and "stage" in groups_df.columns else ""
            print(f"  {bare:<52} {grp:<16} {stg:<8} {raw_med:>10.1f} {qn_med:>10.1f}")


if __name__ == "__main__":
    main()
