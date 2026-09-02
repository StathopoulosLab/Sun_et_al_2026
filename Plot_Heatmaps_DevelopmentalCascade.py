#!/usr/bin/env python3
"""
plot_heatmaps_v8_devcascade.py — developmental-cascade heatmap generator

Based on plot_heatmaps_v7_informative.py, with two changes:

  1. GENES_BY_CATEGORY replaced with the full AP/DV patterning cascade
     (maternal -> gap -> pair-rule -> segment polarity -> homeotic, plus
     DV/mesoderm genes), organized into biologically-standard categories.

  2. Columns now span ALL developmental timepoints (not nc14b only), ordered
     chronologically within each genotype group. Rows (genes) are still
     hierarchically clustered by expression-pattern similarity (same linkage
     as before), but the dendrogram is *oriented*: at every merge, whichever
     branch has the earlier average timepoint of peak expression is drawn on
     top. This keeps genes with similar temporal/expression patterns grouped
     together (true clustering is unchanged) while giving the plot an overall
     top-to-bottom "early -> late" read, per your request.

===============================================================================
THINGS TO VERIFY / EDIT BEFORE RUNNING (marked with >>> below)
===============================================================================
This script assumes DL.load() / DL.ordered_cols() / DL.embryos_by_group()
return embryos across ALL timepoints (not just nc14b). If your data_loader
filters to a single stage internally, you'll need the multi-timepoint
equivalent — check Config_DataLoader.py / Config_SampleMetadata.py.

It also assumes there's *some* way to know each sample's developmental
timepoint. get_timepoint() below tries a few common conventions (an
attribute on the loaded data object, an attribute on sample_config, then a
regex match against the sample ID string). If none of those match your
actual pipeline, edit get_timepoint() directly — that's the one function
that has to be correct for the "early on top" ordering to make sense.
===============================================================================
"""

import os
import sys
import re
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.cluster.hierarchy import linkage, to_tree

import Config_SampleMetadata as SC
import Config_DataLoader as DL

# =============================================================================
# CLI
# =============================================================================

_p = argparse.ArgumentParser()
_p.add_argument("--prefix", default="all_genotypes_plus_set3",
                help="Combined-matrix prefix, matching whatever --combined you "
                     "used for step5b/step5c/step6. Determines both the "
                     "RUVg-corrected matrix (via Config_DataLoader.py) and the wide "
                     "pre-RUVg fallback TPM matrix. Default: all_genotypes_plus_set3")
args = _p.parse_args()

# =============================================================================
# CONFIG
# =============================================================================

OUTDIR = "results/combined_figures/validation_heatmaps/dev_cascade"
os.makedirs(OUTDIR, exist_ok=True)

SYMBOL_MAP = "gene_id_to_symbol.tsv"

# >>> EDIT if your stage labels/spelling differ. Order = chronological.
TIMEPOINT_ORDER = ["nc14b", "nc14late", "gastrulation"]

# Full developmental-cascade gene panel, grouped by classical category.
# Category membership still drives colour-coding (strip + label colour), but
# row ORDER for the un-clustered "base" plots now comes from
# MANUAL_GENE_ORDER below, not from walking this dict.
GENES_BY_CATEGORY = {
    "maternal":          ["cad", "D"],
    "mesoderm_TF":       ["FoxL1", "HLH54F"],
    "head_gap":          ["btd", "oc"],
    "gap":               ["hb", "gt", "kni", "Kr", "tll", "hkb"],
    "DV_patterning":     ["brk", "sog", "hry", "byn"],
    "pair-rule":         ["eve", "run", "ftz", "opa", "prd", "odd"],
    "segment_polarity":  ["en", "wg", "hh"],
    "homeotic":          ["Antp", "Dfd", "Ubx", "abd-A", "Abd-B"],
}

# Manually curated base row order (top -> bottom) for the un-clustered plots.
# This is independent of GENES_BY_CATEGORY grouping — categories interleave
# here, which is why the base-plot category dividers are now drawn wherever
# the category actually changes between consecutive rows, rather than at
# fixed block boundaries (see plot_heatmap). The clustered/dendrogram plots
# are unaffected by this list — they keep computing their own order.
MANUAL_GENE_ORDER = [
    "FoxL1", "Kr", "gt", "ftz", "hry", "eve", "prd", "odd", "opa", "Ubx",
    "D", "sog", "brk", "abd-A", "Abd-B", "Antp", "Dfd", "cad", "run", "en",
    "hh", "wg", "kni", "oc", "hb", "btd", "byn", "tll", "hkb", "HLH54F",
]

# Colour per category (left-side strip + italic label colour)
CATEGORY_COLOURS = {
    "maternal":          "#9E9AC8",
    "mesoderm_TF":       "#74C476",
    "head_gap":          "#41B6C4",
    "gap":               "#6BAED6",
    "DV_patterning":     "#FDAE6B",
    "pair-rule":         "#E07B54",
    "segment_polarity":  "#F768A1",
    "homeotic":          "#B2182B",
}

# Group colours — override from sample_config as needed
GROUP_COLOURS = {
    **SC.GROUP_COLOURS,
    "FoxL1_BOTv":   "#238b45",
    "HLH54F_BOTCv": "#e377c2",
}
GROUP_ORDER = SC.GROUP_ORDER

# ── Outlier exclusions ────────────────────────────────────────────────────────
EXCLUDE_EMBRYOS: list[tuple[str, str]] = [
    ("FoxL1_BOTv", "bc27"),
]

def _is_excluded(sample_id: str) -> bool:
    grp = group_of.get(sample_id, "")
    return any(grp == eg and eb in sample_id for eg, eb in EXCLUDE_EMBRYOS)

CLUSTER_METHOD = "average"
CLUSTER_METRIC = "euclidean"

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
# =============================================================================

TPM_WIDE_PATH = f"results/combined/{args.prefix}_qc_tpm_matrix.tsv"

print(f"[heatmap] Loading data (prefix='{args.prefix}')...")
_data = DL.load(prefix=args.prefix, use_qc=True, use_quantile=True)

tpm_for_z = _data.tpm_qn if _data.tpm_qn is not None else _data.tpm
group_of  = _data.group_of

if os.path.exists(TPM_WIDE_PATH):
    tpm_wide = pd.read_csv(TPM_WIDE_PATH, sep="\t", index_col=0)
    print(f"[heatmap] Loaded wide TPM: {tpm_wide.shape[0]:,} genes x {tpm_wide.shape[1]} embryos")
else:
    print(f"[heatmap] WARNING: {TPM_WIDE_PATH} not found — falling back to RUVg matrix for resolution")
    tpm_wide = tpm_for_z

# Force numeric dtype (any stray non-numeric cell becomes NaN rather than
# silently making the whole row/column dtype='object', which breaks log2()
# and other numeric ops downstream), and collapse duplicate FBgn rows if any
# (keeps the first occurrence — same behaviour as a plain .loc[fbgn] lookup
# would have silently assumed anyway).
tpm_wide = tpm_wide.apply(pd.to_numeric, errors="coerce")
_dupe_mask = tpm_wide.index.duplicated()
if _dupe_mask.any():
    print(f"[heatmap] WARNING: {int(_dupe_mask.sum())} duplicate FBgn row(s) in "
          f"{TPM_WIDE_PATH} — keeping first occurrence of each")
    tpm_wide = tpm_wide[~_dupe_mask]

tpm_all = tpm_wide

by_group     = DL.embryos_by_group(_data)
ordered_cols = DL.ordered_cols(_data)

_before = len(ordered_cols)
ordered_cols = [c for c in ordered_cols if not _is_excluded(c)]
_n_excl = _before - len(ordered_cols)
if _n_excl:
    print(f"[heatmap] Excluded {_n_excl} outlier embryo(s) per EXCLUDE_EMBRYOS")

# =============================================================================
# TIMEPOINT RESOLUTION
# >>> EDIT get_timepoint() if none of these strategies match your pipeline.
# =============================================================================

TIMEPOINT_INDEX = {tp: i for i, tp in enumerate(TIMEPOINT_ORDER)}

def get_timepoint(sample_id: str):
    # Strategy 1: an attribute the data loader itself might expose
    tp_map = getattr(_data, "timepoint_of", None)
    if tp_map and sample_id in tp_map:
        return tp_map[sample_id]
    # Strategy 2: sample_config might expose it instead
    tp_map = getattr(SC, "timepoint_of", None)
    if tp_map and sample_id in tp_map:
        return tp_map[sample_id]
    # Strategy 3: regex against the sample ID string itself
    sid = sample_id.lower()
    if "gast" in sid:
        return "gastrulation"
    if re.search(r"14.?late|14d\b|late", sid):
        return "nc14late"
    if "14b" in sid:
        return "nc14b"
    return None

timepoint_of = {c: get_timepoint(c) for c in ordered_cols}
_unresolved = [c for c, tp in timepoint_of.items() if tp is None]
if _unresolved:
    print(f"[heatmap] WARNING: could not resolve timepoint for {len(_unresolved)} "
          f"sample(s), e.g. {_unresolved[:5]} — edit get_timepoint()")

n_removed = 49 - len(ordered_cols)
qc_note   = (f"{n_removed} embryo(s) removed by QC" if n_removed > 0
             else "all embryos pass QC")

# =============================================================================
# GENE RESOLUTION
# =============================================================================

def resolve_genes(symbols, index):
    resolved = {}
    index_set = set(index)
    for sym in symbols:
        if sym in index_set:
            resolved[sym] = sym
            continue
        if sym in symbol_to_fbgn:
            fbgn = symbol_to_fbgn[sym]
            if fbgn in index_set:
                resolved[sym] = fbgn
                continue
            else:
                print(f"[heatmap] WARNING missing gene: {sym}  "
                      f"(symbol map has {sym} -> {fbgn}, but {fbgn} is not in "
                      f"the TPM matrix — likely filtered out upstream, e.g. by QC)")
                continue
        hits = [k for k in symbol_to_fbgn
                if isinstance(k, str) and k.lower() == sym.lower()]
        if hits:
            canonical = hits[0]
            fbgn = symbol_to_fbgn[canonical]
            if fbgn in index_set:
                resolved[sym] = fbgn
                if canonical != sym:
                    print(f"[heatmap]   '{sym}' matched via '{canonical}'")
                continue
            else:
                print(f"[heatmap] WARNING missing gene: {sym}  "
                      f"(matched symbol map entry '{canonical}' -> {fbgn}, but "
                      f"{fbgn} is not in the TPM matrix — likely filtered out upstream)")
                continue
        print(f"[heatmap] WARNING missing gene: {sym}  "
              f"(no entry for '{sym}' in {SYMBOL_MAP}, case-insensitive — "
              f"check spelling/aliases, e.g. FlyBase synonyms)")
    return resolved


all_gene_symbols = list(MANUAL_GENE_ORDER)

# Sanity check: every gene in the manual order should have a known category
# (for colour-coding); every gene in GENES_BY_CATEGORY should also appear in
# the manual order (otherwise it'd silently vanish from the base plots).
_cat_genes = {g for genes in GENES_BY_CATEGORY.values() for g in genes}
_manual_set = set(MANUAL_GENE_ORDER)
_uncategorized = [g for g in MANUAL_GENE_ORDER if g not in _cat_genes]
_orphaned = [g for g in _cat_genes if g not in _manual_set]
if _uncategorized:
    print(f"[heatmap] WARNING: no category assigned for {_uncategorized} "
          f"(add to GENES_BY_CATEGORY for colour-coding)")
if _orphaned:
    print(f"[heatmap] WARNING: {_orphaned} in GENES_BY_CATEGORY but missing "
          f"from MANUAL_GENE_ORDER — they won't appear in any plot")

resolved = resolve_genes(all_gene_symbols, tpm_all.index)
found_symbols = [g for g in all_gene_symbols if g in resolved]

cat_of = {}
for cat, genes in GENES_BY_CATEGORY.items():
    for g in genes:
        if g not in cat_of:
            cat_of[g] = cat

print(f"[heatmap] Resolved {len(found_symbols)}/{len(all_gene_symbols)} genes")
print(f"[heatmap] Missing: {[g for g in all_gene_symbols if g not in resolved]}")

if not found_symbols:
    sys.exit("[heatmap] ERROR: No genes resolved. Check SYMBOL_MAP and TPM matrix.")

# =============================================================================
# GLOBAL REFERENCES FOR ALTERNATE DISPLAY MODES (log2FC, log2TPM)
#
# >>> ASSUMPTION TO VERIFY: log2FC is computed per gene, per sample, against
# the mean raw TPM across ALL WT embryos at the SAME timepoint (pooled from
# the whole dataset, not just whichever WT samples happen to be in a given
# subset) — pseudocount of 1 on both sides: log2((TPM+1)/(WT_mean_TPM+1)).
# If you wanted a different reference (e.g. vs nc14b baseline instead of WT,
# or vs a specific genotype), tell me and I'll change compute_wt_reference().
# =============================================================================

def compute_wt_reference():
    wt_cols = [c for c in ordered_cols if group_of.get(c) == "WT"]
    ref = {}
    n_missing = 0
    for sym in found_symbols:
        fbgn = resolved[sym]
        per_tp = {}
        row = tpm_wide.loc[fbgn] if fbgn in tpm_wide.index else None
        for tp in TIMEPOINT_ORDER:
            cols_tp = [c for c in wt_cols
                       if timepoint_of.get(c) == tp and row is not None and c in row.index]
            if cols_tp:
                per_tp[tp] = float(row[cols_tp].mean())
            else:
                per_tp[tp] = None
                n_missing += 1
        ref[sym] = per_tp
    if n_missing:
        print(f"[heatmap] NOTE: {n_missing} gene/timepoint WT-reference value(s) unavailable "
              f"(no WT embryos resolved at that timepoint) — those cells will be blank in log2FC plots")
    return ref


def compute_log2fc_matrix(mat, wt_reference):
    """log2((TPM+1)/(WT_mean_TPM_at_matching_timepoint+1)). NaN where no WT reference exists."""
    out = pd.DataFrame(index=mat.index, columns=mat.columns, dtype=float)
    for gene in mat.index:
        ref_by_tp = wt_reference.get(gene, {})
        for c in mat.columns:
            wt_val = ref_by_tp.get(timepoint_of.get(c))
            out.loc[gene, c] = (np.nan if wt_val is None
                                 else np.log2((mat.loc[gene, c] + 1) / (wt_val + 1)))
    return out


WT_REFERENCE = compute_wt_reference()
LOG2FC_RANGE = 3.0   # +/- range for the log2FC colour scale; edit if your data needs more headroom

# Fixed colour-scale ceiling for log2(TPM+1) plots (99th percentile across the
# whole panel x whole dataset), so absolute-expression plots are comparable
# across subsets rather than each auto-scaling to its own max.
_all_tpm_vals = []
for _sym in found_symbols:
    _fbgn = resolved[_sym]
    if _fbgn in tpm_wide.index:
        _cols_present = [c for c in ordered_cols if c in tpm_wide.columns]
        _vals = np.asarray(tpm_wide.loc[_fbgn, _cols_present].values, dtype=float)
        _all_tpm_vals.append(_vals)
_all_tpm_vals = np.concatenate(_all_tpm_vals) if _all_tpm_vals else np.array([0.0])
_all_tpm_vals = _all_tpm_vals[~np.isnan(_all_tpm_vals)]
if _all_tpm_vals.size == 0:
    _all_tpm_vals = np.array([0.0])
LOG2TPM_VMAX = float(np.ceil(np.nanpercentile(np.log2(_all_tpm_vals + 1), 99)))
print(f"[heatmap] log2(TPM+1) colour scale fixed to 0-{LOG2TPM_VMAX:g} (99th percentile across panel)")

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

# =============================================================================
# TIME-ORIENTED CLUSTERING
#
# order_rows_by_time():
#   1. Runs standard hierarchical clustering on the z-scored matrix (this
#      determines WHICH genes group together — unchanged from before).
#   2. Walks the resulting tree and, at every merge, puts whichever branch
#      has the earlier mean "expression-weighted timepoint" on top. This does
#      NOT change which genes cluster together, only which of the two
#      branches is drawn first at each fork — so the final order is
#      simultaneously pattern-clustered AND roughly early-to-late top-to-bottom.
# =============================================================================

def compute_time_centroids(mat, ordered_samples, timepoint_of, timepoint_index):
    """
    Per-gene weighted-average timepoint index, using (non-negative) raw
    expression as the weight. Lower value = expression skewed toward earlier
    timepoints. Samples with unresolved timepoints are excluded.
    """
    centroids = {}
    for gene in mat.index:
        weights, positions = [], []
        for c in ordered_samples:
            tp = timepoint_of.get(c)
            if tp not in timepoint_index:
                continue
            w = max(float(mat.loc[gene, c]), 0.0)
            weights.append(w)
            positions.append(timepoint_index[tp])
        if not positions:
            centroids[gene] = 0.0
        elif sum(weights) == 0:
            centroids[gene] = float(np.mean(positions))
        else:
            centroids[gene] = float(np.average(positions, weights=weights))
    return centroids


def orient_and_layout_dendrogram(Z, leaf_names, centroids):
    """
    Walk the linkage tree; at each merge, place the lower-centroid branch
    first (top). Returns:
      row_order   : leaf_names reordered top-to-bottom
      line_segs   : list of ((x0, x1), (y0, y1)) segments to draw the
                    dendrogram, in the SAME row-position coordinates used for
                    the heatmap (y=0 at bottom leaf, increasing upward),
                    x = linkage distance (root = max distance)
    """
    tree = to_tree(Z)
    n_leaves = len(leaf_names)
    row_order = []           # filled in leaf-visit order (top to bottom)
    leaf_y = {}               # leaf.id -> y position (assigned once row_order is known)
    line_segs = []

    # Cache of subtree mean centroid, keyed by node id, so we only compute
    # each subtree's mean once (used both to decide branch order and to
    # report the parent's combined mean up the tree).
    _mean_cache = {}
    def _peek_mean(node):
        if node.id in _mean_cache:
            return _mean_cache[node.id]
        if node.is_leaf():
            val = centroids[leaf_names[node.id]]
        else:
            n_l, n_r = node.left.get_count(), node.right.get_count()
            val = (_peek_mean(node.left) * n_l + _peek_mean(node.right) * n_r) / (n_l + n_r)
        _mean_cache[node.id] = val
        return val

    def _walk(node):
        if node.is_leaf():
            row_order.append(node.id)
            return centroids[leaf_names[node.id]]
        l_mean = _peek_mean(node.left)
        r_mean = _peek_mean(node.right)
        if l_mean <= r_mean:
            first, second = node.left, node.right
        else:
            first, second = node.right, node.left
        _walk(first)
        _walk(second)
        n_l, n_r = node.left.get_count(), node.right.get_count()
        l_mean_full = _peek_mean(node.left)
        r_mean_full = _peek_mean(node.right)
        return (l_mean_full * n_l + r_mean_full * n_r) / (n_l + n_r)

    _walk(tree)

    # row_order currently lists leaf ids top-to-bottom (visit order).
    # Convert to (row_index -> y position), with row 0 (top) = highest y.
    for i, leaf_id in enumerate(row_order):
        leaf_y[leaf_id] = n_leaves - i - 1  # bottom-up y, matches imshow row math later

    max_dist = tree.dist if tree.dist > 0 else 1.0

    def _draw(node):
        if node.is_leaf():
            return leaf_y[node.id], 0.0
        y_l, h_l = _draw(node.left)
        y_r, h_r = _draw(node.right)
        h = node.dist
        line_segs.append(((h_l, h), (y_l, y_l)))
        line_segs.append(((h_r, h), (y_r, y_r)))
        line_segs.append(((h, h), (y_l, y_r)))
        return (y_l + y_r) / 2.0, h

    _draw(tree)

    row_order_names = [leaf_names[i] for i in row_order]
    return row_order_names, line_segs, max_dist

# =============================================================================
# PLOT
# =============================================================================

def build_ordered_samples(subset_cols):
    """GROUP_ORDER first, then chronological by timepoint within each group."""
    ordered_samples = []
    for grp in GROUP_ORDER:
        grp_cols = [c for c in subset_cols if group_of[c] == grp]
        grp_cols.sort(key=lambda c: (TIMEPOINT_INDEX.get(timepoint_of.get(c), 999), c))
        ordered_samples += grp_cols

    available_samples = [c for c in ordered_samples if c in tpm_wide.columns]
    if len(available_samples) < len(ordered_samples):
        missing_cols = set(ordered_samples) - set(available_samples)
        print(f"[heatmap] WARNING: {len(missing_cols)} sample(s) not in wide TPM, skipping: {missing_cols}")
    return available_samples


def build_matrices(ordered_samples):
    """Build the raw-TPM matrix and its row z-scores for the given sample columns."""
    ordered_rows_fbgn = [resolved[g] for g in found_symbols]
    mat_rows = []
    for fbgn in ordered_rows_fbgn:
        if fbgn in tpm_for_z.index:
            row = tpm_for_z.loc[fbgn, [c for c in ordered_samples if c in tpm_for_z.columns]]
            if len(row) < len(ordered_samples):
                full = tpm_wide.loc[fbgn, ordered_samples].copy()
                full[row.index] = row.values
                mat_rows.append(full)
            else:
                mat_rows.append(row)
        elif fbgn in tpm_wide.index:
            mat_rows.append(tpm_wide.loc[fbgn, ordered_samples])
        else:
            print(f"[heatmap] WARNING: {fbgn} not in any matrix, filling zeros")
            mat_rows.append(pd.Series(0.0, index=ordered_samples))

    mat = pd.DataFrame(mat_rows, index=found_symbols, columns=ordered_samples)
    z   = zscore_df(mat)
    return mat, z


def compute_gene_order(mat, z, ordered_samples, label):
    """
    Run clustering + time-orientation (see orient_and_layout_dendrogram above)
    and return (row_order_names, line_segs, max_dist), or None if clustering
    isn't possible (e.g. too few genes, or it errors out).
    """
    n_genes = len(z.index)
    if n_genes <= 2:
        print(f"[heatmap] WARNING: too few genes to cluster for {label}; "
              f"falling back to category order")
        return None
    try:
        Z_linkage = linkage(z.values, method=CLUSTER_METHOD, metric=CLUSTER_METRIC)
        centroids = compute_time_centroids(mat, ordered_samples, timepoint_of, TIMEPOINT_INDEX)
        return orient_and_layout_dendrogram(Z_linkage, list(z.index), centroids)
    except Exception as e:
        print(f"[heatmap] WARNING: clustering failed for {label} ({e}); "
              f"falling back to category order")
        return None


def plot_heatmap(name, subset_cols, cluster_genes=True, shared_order=None,
                  mode="zscore", zscore_range=2.5):
    """
    shared_order, if given, is a (row_order_names, line_segs, max_dist) tuple
    computed elsewhere (typically from the WT-included "parent" subset) —
    passing it in makes this plot reuse that gene order/dendrogram instead of
    re-clustering on its own (possibly smaller) sample set, so paired plots
    like ventralized / ventralized_no_wt show genes in the same top-to-bottom
    order for direct, apples-to-apples visual comparison.

    mode: "zscore" (default; row-wise z-score, range set by zscore_range),
          "log2tpm" (log2(TPM+1), fixed 0-LOG2TPM_VMAX colour scale, shared
          across all plots for comparability), or "log2fc" (log2 fold-change
          vs the pooled WT mean at the matching timepoint, see WT_REFERENCE).
    Row order (clustering / MANUAL_GENE_ORDER) is always computed from the
    z-scored data regardless of display mode, so gene order stays consistent
    across modes for the same subset.
    """
    mode_suffix = {"zscore": "" if zscore_range == 2.5 else f"_z{zscore_range:g}",
                   "log2tpm": "_log2tpm",
                   "log2fc": "_log2fc"}[mode]
    label = f"{name}_devcascade" + ("_clustered" if cluster_genes else "")
    file_label = label + mode_suffix
    print(f"[heatmap] Building {file_label}  (n={len(subset_cols)} embryos)...")

    ordered_samples = build_ordered_samples(subset_cols)
    mat, z = build_matrices(ordered_samples)

    n_genes = len(found_symbols)
    n_emb   = len(ordered_samples)

    line_segs, max_dist = None, 1.0
    if cluster_genes:
        order_bundle = shared_order if shared_order is not None else \
            compute_gene_order(mat, z, ordered_samples, label)
        if order_bundle is not None:
            row_order_names, line_segs, max_dist = order_bundle
            # Keep this robust even if a shared order and this subset's gene
            # list don't line up exactly (shouldn't happen with a fixed panel).
            row_order_names = [g for g in row_order_names if g in z.index]
            row_order_names += [g for g in z.index if g not in row_order_names]
            z   = z.loc[row_order_names]
            mat = mat.loc[row_order_names]
        else:
            cluster_genes = False

    genes_for_plot = z.index.tolist()

    fig_w = max(14, n_emb   * 0.55) + (2.6 if cluster_genes else 0)
    fig_h = max(8,  n_genes * 0.45)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    left_margin = 0.34 if cluster_genes else 0.28
    fig.subplots_adjust(left=left_margin, right=0.78, top=0.80, bottom=0.20)

    if mode == "zscore":
        values = z.values
        vmin, vmax = -zscore_range, zscore_range
        cmap = "RdBu_r"
        cbar_label = "Z-score"
        value_desc = f"Z-scored TPM (scale \u00b1{zscore_range:g})"
    elif mode == "log2tpm":
        values = np.log2(mat.values + 1)
        vmin, vmax = 0, LOG2TPM_VMAX
        cmap = "viridis"
        cbar_label = "log2(TPM + 1)"
        value_desc = "log2(TPM+1), fixed scale across subsets"
    elif mode == "log2fc":
        values = compute_log2fc_matrix(mat, WT_REFERENCE).values
        vmin, vmax = -LOG2FC_RANGE, LOG2FC_RANGE
        cmap = "RdBu_r"
        cbar_label = "log2FC vs WT (same timepoint)"
        value_desc = "log2FC vs pooled WT mean at matching timepoint (pseudocount=1)"
    else:
        raise ValueError(f"Unknown mode: {mode}")

    im = ax.imshow(values, aspect="auto", cmap=cmap,
                   vmin=vmin, vmax=vmax, interpolation="nearest")

    ax.set_xticks([])

    ax.set_yticks(range(n_genes))
    ax.set_yticklabels(genes_for_plot, fontstyle="italic",
                       fontweight="bold", fontsize=10)
    for tick, gene in zip(ax.get_yticklabels(), genes_for_plot):
        tick.set_color(CATEGORY_COLOURS.get(cat_of.get(gene, "other"), "black"))
    ax.tick_params(axis="y", pad=4)

    if not cluster_genes:
        for i in range(n_genes - 1):
            cat_i   = cat_of.get(genes_for_plot[i], "other")
            cat_ip1 = cat_of.get(genes_for_plot[i + 1], "other")
            if cat_i != cat_ip1:
                ax.axhline(i + 0.5, color="black", lw=0.6, ls=":")

    pos = 0
    grp_boundaries = []
    grp_label_positions = []
    # timepoint boundaries drawn as thin dotted lines within each group
    tp_boundaries = []
    for grp in GROUP_ORDER:
        gc = [c for c in ordered_samples if group_of[c] == grp]
        if not gc:
            continue
        if pos > 0:
            grp_boundaries.append(pos - 0.5)
        grp_label_positions.append((pos + len(gc) / 2, grp))
        prev_tp = None
        for i, c in enumerate(gc):
            tp = timepoint_of.get(c)
            if prev_tp is not None and tp != prev_tp:
                tp_boundaries.append(pos + i - 0.5)
            prev_tp = tp
        pos += len(gc)
    for b in grp_boundaries:
        ax.axvline(b, color="black", lw=2.0)
    for b in tp_boundaries:
        ax.axvline(b, color="grey", lw=0.6, ls=":")

    strip_ax = ax.inset_axes([0, 1.01, 1, 0.04], transform=ax.transAxes)
    strip_ax.set_xlim(0, n_emb)
    strip_ax.set_ylim(0, 1)
    strip_ax.axis("off")
    for i, col in enumerate(ordered_samples):
        strip_ax.add_patch(plt.Rectangle(
            (i, 0), 1, 1,
            color=GROUP_COLOURS.get(group_of[col], "#888888"),
            transform=strip_ax.transData))

    label_ax = ax.inset_axes([0, 1.055, 1, 0.05], transform=ax.transAxes)
    label_ax.set_xlim(0, n_emb)
    label_ax.set_ylim(0, 1)
    label_ax.axis("off")
    for xpos, lbl in grp_label_positions:
        label_ax.text(xpos, 0.2, lbl, ha="center", va="bottom",
                      fontsize=8, fontweight="bold",
                      color=GROUP_COLOURS.get(lbl, "black"))

    cat_ax = ax.inset_axes([-0.18, 0, 0.04, 1], transform=ax.transAxes)
    cat_ax.set_xlim(0, 1)
    cat_ax.set_ylim(0, n_genes)
    cat_ax.axis("off")
    for i, gene in enumerate(genes_for_plot):
        cat_ax.add_patch(plt.Rectangle(
            (0, n_genes - i - 1), 1, 1,
            color=CATEGORY_COLOURS.get(cat_of.get(gene, "other"), "grey"),
            transform=cat_ax.transData))

    # ── Time-oriented dendrogram to the LEFT of the category strip ───────────
    if cluster_genes and line_segs is not None:
        dendro_ax = ax.inset_axes([-0.34, 0, 0.13, 1], transform=ax.transAxes)
        dendro_ax.set_xlim(max_dist * 1.05, 0)   # root on the left, leaves on the right
        dendro_ax.set_ylim(0, n_genes)
        dendro_ax.axis("off")
        for (x0, x1), (y0, y1) in line_segs:
            dendro_ax.plot([x0, x1], [y0 + 0.5, y1 + 0.5], color="black", lw=0.9)

    cbar = fig.colorbar(im, ax=ax, fraction=0.018, pad=0.01)
    cbar.set_label(cbar_label, fontsize=9)

    cat_patches = [mpatches.Patch(color=v, label=k)
                   for k, v in CATEGORY_COLOURS.items()
                   if k in set(cat_of.values())]
    grp_patches = [mpatches.Patch(color=GROUP_COLOURS.get(g, "#888"), label=g)
                   for g in GROUP_ORDER
                   if any(group_of[c] == g for c in ordered_samples)]
    ax.legend(handles=cat_patches + grp_patches,
              loc="lower left", bbox_to_anchor=(1.02, 0.0),
              frameon=False, fontsize=8, title="category / group")

    tp_chain = " \u2192 ".join(TIMEPOINT_ORDER)  # precomputed to avoid backslash inside f-string expr
    if cluster_genes:
        subtitle = (f"{value_desc}  |  rows clustered by expression pattern "
                    f"({CLUSTER_METHOD}/{CLUSTER_METRIC}), oriented early\u2192late top\u2192bottom  |  "
                    f"cols: {tp_chain} within group")
    else:
        subtitle = (f"{value_desc}  |  category strip left  |  "
                    f"cols: {tp_chain} within group")
    ax.set_title(
        f"Developmental cascade heatmap — {file_label}  (n={n_emb})  [{qc_note}]\n{subtitle}",
        fontsize=10, pad=40)

    out_png = os.path.join(OUTDIR, f"heatmap_{file_label}.png")
    out_pdf = os.path.join(OUTDIR, f"heatmap_{file_label}.pdf")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf,           bbox_inches="tight")
    plt.close(fig)
    print(f"[heatmap] Saved {file_label}  ->  {out_png}")


# =============================================================================
# RUN ALL SUBSETS
# =============================================================================

for base_name, base_cols in _BASE_SUBSETS.items():
    # Compute the reference gene order ONCE per family, from the WT-included
    # subset, and reuse it for the "_no_wt" variant too — so e.g. ventralized
    # and ventralized_no_wt show genes in the same top-to-bottom order and
    # are directly comparable side by side. This order (from z-scored data)
    # is reused for every display mode below, so gene order is consistent
    # across the zscore / log2tpm / log2fc versions of the same subset too.
    ref_samples = build_ordered_samples(SUBSETS[base_name])
    ref_mat, ref_z = build_matrices(ref_samples)
    shared_order = compute_gene_order(ref_mat, ref_z, ref_samples, base_name)

    for variant_name in (base_name, f"{base_name}_no_wt"):
        cols = SUBSETS[variant_name]

        # Existing default outputs — unchanged filenames, unchanged behaviour
        plot_heatmap(variant_name, cols, cluster_genes=False)
        plot_heatmap(variant_name, cols, cluster_genes=True, shared_order=shared_order)

        # NEW: wider z-score colour scale (+/-4) — same data/clustering, just
        # less saturated so within-genotype variability reads more clearly
        plot_heatmap(variant_name, cols, cluster_genes=False,
                     mode="zscore", zscore_range=4.0)
        plot_heatmap(variant_name, cols, cluster_genes=True, shared_order=shared_order,
                     mode="zscore", zscore_range=4.0)

        # NEW: log2(TPM+1) — absolute expression level, not row-normalized
        plot_heatmap(variant_name, cols, cluster_genes=False, mode="log2tpm")
        plot_heatmap(variant_name, cols, cluster_genes=True, shared_order=shared_order,
                     mode="log2tpm")

        # NEW: log2FC vs pooled WT mean at the matching timepoint
        plot_heatmap(variant_name, cols, cluster_genes=False, mode="log2fc")
        plot_heatmap(variant_name, cols, cluster_genes=True, shared_order=shared_order,
                     mode="log2fc")

print("[heatmap] Done")
