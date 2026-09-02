#!/usr/bin/env python3
"""
plot_heatmaps_v7_informative.py — nc14b heatmap generator, informative-gene subset

Identical to plot_heatmaps_v7.py except GENES_BY_CATEGORY has been trimmed to
remove genes that show little or no contrast between FoxL1_BOTv and HLH54F_BOTCv
(identified by visual inspection of ventralized_clustered and
ventralized_no_wt_clustered outputs):

  Removed:  en, wg, sog, D, Opa   (flat / uninformative across both genotypes)
  Optional: cic, ths               (kept — borderline; remove here if desired)

The "other" and "signaling" categories are dropped entirely as a result.
TF category reduced to: HLH54F, foxl1, Oc.
Maternal category reduced to: bcd, cad, cic, dl  (sog removed).

Output filenames use the same suffix scheme as v7 but are written to a
subdirectory:  results/combined_figures/validation_heatmaps/informative/
so they don't overwrite the full-gene outputs.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.cluster.hierarchy import linkage, dendrogram

import Config_SampleMetadata as SC
import Config_DataLoader as DL

# =============================================================================
# CLI
# =============================================================================

_p = argparse.ArgumentParser()
_p.add_argument("--prefix", default="all_genotypes_plus_set3",
                help="Combined-matrix prefix, matching whatever --combined you "
                     "used for step5b/step5c/step6. Default: all_genotypes_plus_set3")
args = _p.parse_args()

# =============================================================================
# CONFIG
# =============================================================================

OUTDIR = "results/combined_figures/validation_heatmaps/informative"
os.makedirs(OUTDIR, exist_ok=True)

SYMBOL_MAP = "gene_id_to_symbol.tsv"

# Gene list — uninformative genes (en, wg, sog, D, Opa) removed.
# Edit here to further trim or restore genes.
GENES_BY_CATEGORY = {
    "gap":       ["hb", "Kr", "kni", "gt", "tll", "hkb"],
    "pair-rule": ["run", "eve", "H", "ftz", "odd", "prd"],
    "maternal":  ["cad", "cic", "dl"],
    "mesoderm":  ["sna", "twi", "HLH54F", "NetA"],
    "TF":        ["HLH54F", "foxl1", "Oc"],
    "signaling": ["pyr"],
}

# Colour per category (left-side strip + italic label colour)
CATEGORY_COLOURS = {
    "gap":       "#6BAED6",
    "pair-rule": "#E07B54",
    "maternal":  "#9E9AC8",
    "mesoderm":  "#74C476",
    "TF":        "#F768A1",
    "signaling": "#E6AB02",
}

# Group colours — override BOTv/BOTCv from sample_config
GROUP_COLOURS = {
    **SC.GROUP_COLOURS,
    "FoxL1_BOTv":   "#238b45",   # green
    "HLH54F_BOTCv": "#e377c2",   # pink
}
GROUP_ORDER = SC.GROUP_ORDER

# ── Outlier exclusions ────────────────────────────────────────────────────────
# (group, barcode_substring) pairs — applied to ordered_cols before subsetting.
EXCLUDE_EMBRYOS: list[tuple[str, str]] = [
    ("FoxL1_BOTv", "bc27"),
]

def _is_excluded(sample_id: str) -> bool:
    grp = group_of.get(sample_id, "")
    return any(grp == eg and eb in sample_id for eg, eb in EXCLUDE_EMBRYOS)

# Row-clustering settings for the "_clustered" heatmap versions.
# euclidean is used (rather than correlation) because rows are already
# z-scored (mean 0, std 1), where euclidean distance is monotonically
# equivalent to correlation distance, and it avoids NaNs for any
# zero-variance gene rows (constant expression -> std=0 -> z=0 everywhere).
CLUSTER_METHOD = "average"     # linkage method: average, complete, ward, single
CLUSTER_METRIC = "euclidean"   # distance metric

# =============================================================================
# LOAD SYMBOL MAP
# =============================================================================

symbol_to_fbgn: dict = {}
fbgn_to_symbol: dict = {}

if os.path.exists(SYMBOL_MAP):
    sym_df = pd.read_csv(SYMBOL_MAP, sep="\t").dropna(subset=["Symbol", "GeneID"])
    sym_df = sym_df[sym_df["Symbol"].str.strip() != ""]
    symbol_to_fbgn = dict(zip(sym_df["Symbol"], sym_df["GeneID"]))
    fbgn_to_symbol = dict(zip(sym_df["GeneID"],  sym_df["Symbol"]))
    print(f"[heatmap] Loaded {len(symbol_to_fbgn):,} gene symbol mappings")
else:
    print("[heatmap] WARNING: SYMBOL_MAP not found — using raw index only")

# =============================================================================
# LOAD DATA
# Strategy:
#   tpm_for_z : RUVg QN matrix (3,358 genes) — used for z-scores where available
#   tpm_wide  : pre-RUVg QC TPM (5,000 genes) — fallback for genes dropped by RUVg
#               loaded directly from results/combined/all_genotypes_qc_tpm_matrix.tsv
# This recovers genes like HLH54F, foxl1, bcd, dl, slp1, en, pyr, NetA which are
# present in the QC-filtered TPM but dropped during RUVg count correction.
# =============================================================================

TPM_WIDE_PATH = f"results/combined/{args.prefix}_qc_tpm_matrix.tsv"

print(f"[heatmap] Loading data (prefix='{args.prefix}')...")
_data = DL.load(prefix=args.prefix, use_qc=True, use_quantile=True)

tpm_for_z = _data.tpm_qn if _data.tpm_qn is not None else _data.tpm
group_of  = _data.group_of

# Load the wider TPM matrix for gene resolution + fallback values
if os.path.exists(TPM_WIDE_PATH):
    tpm_wide = pd.read_csv(TPM_WIDE_PATH, sep="	", index_col=0)
    print(f"[heatmap] Loaded wide TPM: {tpm_wide.shape[0]:,} genes x {tpm_wide.shape[1]} embryos")
else:
    print(f"[heatmap] WARNING: {TPM_WIDE_PATH} not found — falling back to RUVg matrix for resolution")
    tpm_wide = tpm_for_z

# Use wide matrix for gene resolution (has full gene index before RUVg filtering)
tpm_all = tpm_wide

by_group     = DL.embryos_by_group(_data)
ordered_cols = DL.ordered_cols(_data)

# ── Apply outlier exclusions (v6) ─────────────────────────────────────────────
_before = len(ordered_cols)
ordered_cols = [c for c in ordered_cols if not _is_excluded(c)]
_n_excl = _before - len(ordered_cols)
if _n_excl:
    print(f"[heatmap] Excluded {_n_excl} outlier embryo(s) per EXCLUDE_EMBRYOS")

n_removed = 49 - len(ordered_cols)
qc_note   = (f"{n_removed} embryo(s) removed by QC" if n_removed > 0
             else "all embryos pass QC")

# =============================================================================
# GENE RESOLUTION
# Resolves against tpm_all.index so genes filtered from QN matrix are still found.
# Returns dict {symbol: fbgn_id}.
# =============================================================================

def resolve_genes(symbols, index):
    resolved = {}
    index_set = set(index)
    for sym in symbols:
        if sym in index_set:
            resolved[sym] = sym
            continue
        if sym in symbol_to_fbgn and symbol_to_fbgn[sym] in index_set:
            resolved[sym] = symbol_to_fbgn[sym]
            continue
        hits = [k for k in symbol_to_fbgn
                if isinstance(k, str) and k.lower() == sym.lower()
                and symbol_to_fbgn[k] in index_set]
        if hits:
            canonical = hits[0]
            resolved[sym] = symbol_to_fbgn[canonical]
            if canonical != sym:
                print(f"[heatmap]   '{sym}' matched via '{canonical}'")
            continue
        print(f"[heatmap] WARNING missing gene: {sym}")
    return resolved


# Build flat ordered gene list preserving category order, deduplicating
seen = set()
all_gene_symbols = []
for cat, genes in GENES_BY_CATEGORY.items():
    for g in genes:
        if g not in seen:
            all_gene_symbols.append(g)
            seen.add(g)

resolved = resolve_genes(all_gene_symbols, tpm_all.index)
found_symbols = [g for g in all_gene_symbols if g in resolved]

# category lookup for each found symbol
cat_of = {}
for cat, genes in GENES_BY_CATEGORY.items():
    for g in genes:
        if g not in cat_of:   # first category wins for duplicates
            cat_of[g] = cat

print(f"[heatmap] Resolved {len(found_symbols)}/{len(all_gene_symbols)} genes")
print(f"[heatmap] Missing: {[g for g in all_gene_symbols if g not in resolved]}")

if not found_symbols:
    sys.exit("[heatmap] ERROR: No genes resolved. Check SYMBOL_MAP and TPM matrix.")

# =============================================================================
# SAMPLE SUBSETS
# =============================================================================

def subset_all_no_runt(cols):
    return [c for c in cols if group_of[c] != "runt"]

def subset_ventralized(cols):
    return [c for c in cols if group_of[c] in ("WT", "FoxL1_BOTv", "HLH54F_BOTCv")]

def subset_non_ventralized(cols):
    return [c for c in cols if group_of[c] not in ("runt", "FoxL1_BOTv", "HLH54F_BOTCv")]


_BASE_SUBSETS = {
    "all_no_runt":      subset_all_no_runt(ordered_cols),
    "ventralized":      subset_ventralized(ordered_cols),
    "non_ventralized":  subset_non_ventralized(ordered_cols),
}

# For each base subset, also generate a "_no_wt" variant with WT embryos
# excluded, so e.g. ventralized_no_wt compares FoxL1_BOTv vs HLH54F_BOTCv only.
SUBSETS = {}
for _name, _cols in _BASE_SUBSETS.items():
    SUBSETS[_name] = _cols
    SUBSETS[f"{_name}_no_wt"] = [c for c in _cols if group_of[c] != "WT"]

# =============================================================================
# Z-SCORE HELPER
# =============================================================================

def zscore_df(df):
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        z = df.subtract(df.mean(axis=1), axis=0).divide(df.std(axis=1), axis=0)
    return z.fillna(0)


def cluster_cols_within_groups(z, ordered_samples, group_of):
    """
    Reorder columns by hierarchical clustering independently within each group.
    Group boundaries (and GROUP_ORDER) are preserved — only the order of
    replicates within each group changes.  Returns a new column order list.
    """
    new_order = []
    # Walk groups in the order they appear in ordered_samples (preserves GROUP_ORDER)
    seen_groups = []
    for c in ordered_samples:
        g = group_of[c]
        if g not in seen_groups:
            seen_groups.append(g)

    for grp in seen_groups:
        grp_cols = [c for c in ordered_samples if group_of[c] == grp]
        if len(grp_cols) <= 2:
            # Nothing to cluster with <3 replicates; keep as-is
            new_order.extend(grp_cols)
            continue
        sub = z[grp_cols]
        try:
            # Cluster on transposed matrix (samples as rows, genes as columns)
            Z_col = linkage(sub.T.values, method=CLUSTER_METHOD, metric=CLUSTER_METRIC)
            col_order = dendrogram(Z_col, no_plot=True)["leaves"]
            new_order.extend([grp_cols[i] for i in col_order])
        except Exception as e:
            print(f"[heatmap] WARNING: col clustering failed for group {grp} ({e}); "
                  f"keeping original order")
            new_order.extend(grp_cols)
    return new_order

# =============================================================================
# PLOT
# =============================================================================

def plot_heatmap(name, subset_cols, cluster_genes=False, cluster_cols=False):
    # Build output label / filename suffix
    if cluster_genes and cluster_cols:
        label = f"{name}_both_clustered"
    elif cluster_genes:
        label = f"{name}_clustered"
    elif cluster_cols:
        label = f"{name}_col_clustered"
    else:
        label = name
    print(f"[heatmap] Building {label}  (n={len(subset_cols)} embryos)...")

    # ── Build ordered row/col matrices ───────────────────────────────────────
    ordered_rows_fbgn = [resolved[g] for g in found_symbols]

    # Columns: respect GROUP_ORDER within the subset
    ordered_samples = []
    for grp in GROUP_ORDER:
        ordered_samples += [c for c in subset_cols if group_of[c] == grp]

    # Build matrix: prefer RUVg QN values (tpm_for_z); fall back to wide QC TPM
    # for genes dropped by RUVg (HLH54F, foxl1, bcd, dl, slp1, en, pyr, NetA etc.)
    # Columns must exist in both matrices; restrict to samples present in wide TPM
    available_samples = [c for c in ordered_samples if c in tpm_wide.columns]
    if len(available_samples) < len(ordered_samples):
        missing_cols = set(ordered_samples) - set(available_samples)
        print(f"[heatmap] WARNING: {len(missing_cols)} sample(s) not in wide TPM, skipping: {missing_cols}")
    ordered_samples = available_samples

    mat_rows = []
    for fbgn in ordered_rows_fbgn:
        if fbgn in tpm_for_z.index:
            # Use RUVg-corrected QN values for normalisation-quality z-scores
            row = tpm_for_z.loc[fbgn, [c for c in ordered_samples if c in tpm_for_z.columns]]
            # Fill any samples missing from tpm_for_z from wide TPM
            if len(row) < len(ordered_samples):
                full = tpm_wide.loc[fbgn, ordered_samples].copy()
                full[row.index] = row.values
                mat_rows.append(full)
            else:
                mat_rows.append(row)
        elif fbgn in tpm_wide.index:
            # Gene was dropped by RUVg — use pre-RUVg QC TPM directly
            mat_rows.append(tpm_wide.loc[fbgn, ordered_samples])
        else:
            print(f"[heatmap] WARNING: {fbgn} not in any matrix, filling zeros")
            mat_rows.append(pd.Series(0.0, index=ordered_samples))

    mat = pd.DataFrame(mat_rows, index=found_symbols, columns=ordered_samples)
    z   = zscore_df(mat)

    # ── Optional within-group column clustering (v7) ──────────────────────────
    if cluster_cols and len(ordered_samples) > 1:
        new_col_order = cluster_cols_within_groups(z, ordered_samples, group_of)
        ordered_samples = new_col_order
        z   = z[ordered_samples]
        mat = mat[ordered_samples]

    n_genes = len(found_symbols)
    n_emb   = len(ordered_samples)

    # ── Optional row clustering ────────────────────────────────────────────
    # Reorders genes by similarity of their z-scored expression profile across
    # the embryos in this subset, instead of by the fixed category order.
    Z_linkage = None
    if cluster_genes and n_genes > 2:
        try:
            Z_linkage = linkage(z.values, method=CLUSTER_METHOD, metric=CLUSTER_METRIC)
            row_order = dendrogram(Z_linkage, no_plot=True)["leaves"]
            z   = z.iloc[row_order]
            mat = mat.iloc[row_order]
        except Exception as e:
            print(f"[heatmap] WARNING: clustering failed for {label} ({e}); "
                  f"falling back to category order")
            cluster_genes = False
            Z_linkage = None
    elif cluster_genes:
        print(f"[heatmap] WARNING: too few genes to cluster for {label}; "
              f"falling back to category order")
        cluster_genes = False

    genes_for_plot = z.index.tolist()

    # ── Figure layout ─────────────────────────────────────────────────────────
    fig_w = max(14, n_emb   * 0.60) + (2.6 if cluster_genes else 0)
    fig_h = max(8,  n_genes * 0.55)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    # left margin is wider for clustered plots to make room for the dendrogram
    left_margin = 0.34 if cluster_genes else 0.28
    fig.subplots_adjust(left=left_margin, right=0.78, top=0.80, bottom=0.20)

    # ── Main heatmap ──────────────────────────────────────────────────────────
    im = ax.imshow(z.values, aspect="auto", cmap="RdBu_r",
                   vmin=-2.5, vmax=2.5, interpolation="nearest")

    # ── x-axis: hide individual embryo labels (group strip carries this info) ──
    ax.set_xticks([])

    # ── y-axis labels — bold italic, coloured by category ────────────────────
    ax.set_yticks(range(n_genes))
    ax.set_yticklabels(genes_for_plot, fontstyle="italic",
                       fontweight="bold", fontsize=10)
    for tick, gene in zip(ax.get_yticklabels(), genes_for_plot):
        tick.set_color(CATEGORY_COLOURS.get(cat_of.get(gene, "other"), "black"))
    # small pad — category strip is placed further left so no collision
    ax.tick_params(axis="y", pad=4)

    # ── Dotted horizontal lines between gene categories ───────────────────────
    # Only meaningful when genes are still in fixed category order — once
    # rows are hierarchically clustered, categories are no longer contiguous.
    if not cluster_genes:
        cumulative = 0
        for cat, genes in GENES_BY_CATEGORY.items():
            found_in_cat = [g for g in genes if g in genes_for_plot]
            cumulative += len(found_in_cat)
            if 0 < cumulative < n_genes:
                ax.axhline(cumulative - 0.5, color="black", lw=1.0, ls=":")

    # ── Vertical black lines between sample groups ────────────────────────────
    pos = 0
    grp_boundaries = []
    grp_label_positions = []
    for grp in GROUP_ORDER:
        gc = [c for c in ordered_samples if group_of[c] == grp]
        if not gc:
            continue
        if pos > 0:
            grp_boundaries.append(pos - 0.5)
        grp_label_positions.append((pos + len(gc) / 2, grp))
        pos += len(gc)
    for b in grp_boundaries:
        ax.axvline(b, color="black", lw=2.0)

    # ── Colour strip above heatmap (group colours) ────────────────────────────
    strip_ax = ax.inset_axes([0, 1.01, 1, 0.04], transform=ax.transAxes)
    strip_ax.set_xlim(0, n_emb)
    strip_ax.set_ylim(0, 1)
    strip_ax.axis("off")
    for i, col in enumerate(ordered_samples):
        strip_ax.add_patch(plt.Rectangle(
            (i, 0), 1, 1,
            color=GROUP_COLOURS.get(group_of[col], "#888888"),
            transform=strip_ax.transData))

    # ── Group name labels above the colour strip ──────────────────────────────
    label_ax = ax.inset_axes([0, 1.055, 1, 0.05], transform=ax.transAxes)
    label_ax.set_xlim(0, n_emb)
    label_ax.set_ylim(0, 1)
    label_ax.axis("off")
    for xpos, lbl in grp_label_positions:
        label_ax.text(xpos, 0.2, lbl, ha="center", va="bottom",
                      fontsize=8, fontweight="bold",
                      color=GROUP_COLOURS.get(lbl, "black"))

    # ── Category colour strip to the LEFT of y-axis labels ───────────────────
    # -0.18 in axes fraction pushes the strip well clear of the gene name text
    cat_ax = ax.inset_axes([-0.18, 0, 0.04, 1], transform=ax.transAxes)
    cat_ax.set_xlim(0, 1)
    cat_ax.set_ylim(0, n_genes)
    cat_ax.axis("off")
    for i, gene in enumerate(genes_for_plot):
        cat_ax.add_patch(plt.Rectangle(
            (0, n_genes - i - 1), 1, 1,
            color=CATEGORY_COLOURS.get(cat_of.get(gene, "other"), "grey"),
            transform=cat_ax.transData))

    # ── Dendrogram to the LEFT of the category strip (clustered plots only) ──
    if cluster_genes and Z_linkage is not None:
        R = dendrogram(Z_linkage, no_plot=True)
        icoord, dcoord = R["icoord"], R["dcoord"]
        max_dist = max((max(d) for d in dcoord), default=1.0) or 1.0

        def _leaf_to_y(p):
            # scipy leaf positions run 5, 15, 25, ... left-to-right; convert
            # to the row-center coordinate used by cat_ax/heatmap (row 0 = top)
            leaf_rank = (p - 5.0) / 10.0
            return n_genes - leaf_rank - 0.5

        dendro_ax = ax.inset_axes([-0.34, 0, 0.13, 1], transform=ax.transAxes)
        dendro_ax.set_xlim(max_dist * 1.05, 0)   # root on the left, leaves on the right
        dendro_ax.set_ylim(0, n_genes)
        dendro_ax.axis("off")
        for ic, dc in zip(icoord, dcoord):
            dendro_ax.plot(dc, [_leaf_to_y(p) for p in ic], color="black", lw=0.9)

    # ── Colorbar ──────────────────────────────────────────────────────────────
    cbar = fig.colorbar(im, ax=ax, fraction=0.018, pad=0.01)
    cbar.set_label("Z-score", fontsize=9)

    # ── Legend (category + group) ─────────────────────────────────────────────
    cat_patches = [mpatches.Patch(color=v, label=k)
                   for k, v in CATEGORY_COLOURS.items()
                   if k in set(cat_of.values())]
    grp_patches = [mpatches.Patch(color=GROUP_COLOURS.get(g, "#888"), label=g)
                   for g in GROUP_ORDER
                   if any(group_of[c] == g for c in ordered_samples)]
    ax.legend(handles=cat_patches + grp_patches,
              loc="lower left", bbox_to_anchor=(1.02, 0.0),
              frameon=False, fontsize=8, title="category / group")

    # ── Title ─────────────────────────────────────────────────────────────────
    if cluster_genes and cluster_cols:
        subtitle = (f"Z-scored TPM  |  rows & cols hierarchically clustered "
                    f"({CLUSTER_METHOD} linkage, {CLUSTER_METRIC})  |  cols within groups")
    elif cluster_genes:
        subtitle = (f"Z-scored TPM  |  rows hierarchically clustered "
                    f"({CLUSTER_METHOD} linkage, {CLUSTER_METRIC})  |  group strip top")
    elif cluster_cols:
        subtitle = (f"Z-scored TPM  |  category strip left  |  "
                    f"cols clustered within groups ({CLUSTER_METHOD}/{CLUSTER_METRIC})")
    else:
        subtitle = "Z-scored TPM  |  category strip left  |  group strip top"
    ax.set_title(
        f"nc14b heatmap — {label}  (n={n_emb})  [{qc_note}]\n{subtitle}",
        fontsize=10, pad=40)

    # ── Save ──────────────────────────────────────────────────────────────────
    out_png = os.path.join(OUTDIR, f"heatmap_{label}.png")
    out_pdf = os.path.join(OUTDIR, f"heatmap_{label}.pdf")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf,           bbox_inches="tight")
    plt.close(fig)
    print(f"[heatmap] Saved {label}  ->  {out_png}")


# =============================================================================
# RUN ALL SUBSETS (each in category-order + clustered versions)
# =============================================================================

for name, cols in SUBSETS.items():
    plot_heatmap(name, cols, cluster_genes=False, cluster_cols=False)
    plot_heatmap(name, cols, cluster_genes=True,  cluster_cols=False)
    plot_heatmap(name, cols, cluster_genes=False, cluster_cols=True)
    plot_heatmap(name, cols, cluster_genes=True,  cluster_cols=True)

print("[heatmap] Done")
