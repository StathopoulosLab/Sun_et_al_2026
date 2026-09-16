#!/usr/bin/env python3
"""
Analysis_GeneClusters_BOTvBOTCv.py — Exploratory expression clustering to find
candidate "anchor genes" distinguishing FoxL1_BOTv (BOTv) vs HLH54F_BOTCv (BOTCv)
across developmental time.

WHAT THIS DOES (high level)
----------------------------
Your existing scripts (Plot_Heatmaps_Informative.py, Plot_PositionalClassification.py)
visualize a hand-picked gene list. This script flips that around: it starts from
the broad set of expressed genes (median TPM > 1 in at least one genotype x stage
group) and asks the data which genes cluster together, with no prior gene list
required. The goal is to surface new candidate "anchor genes" — genes with a
consistent, time-resolved difference between BOTv and BOTCv — that you can then
manually inspect, validate, and potentially fold into the curated gene lists used
by your other scripts.

STRATEGY
--------
1. Gene filter
   Median TPM > 1 within at least one (genotype x stage) cell. Using the median
   (not mean, not "any single sample") makes the filter robust to one dropout
   or outlier embryo without needing a manual exclude list, consistent with your
   "still vetting replicates" note. Computed on RUVg-corrected, QC-filtered TPM
   with fallback to the wider pre-RUVg QC TPM for genes RUVg correction drops
   (matches the convention in Plot_Heatmaps_Informative.py).

2. Collapse to genotype x stage profiles
   For each gene, take the per-(genotype, stage) median log2(TPM+1) across
   replicates -> a genes x 6 matrix (2 genotypes x 3 stages: nc14b, nc14d, gastr).
   This is the "relaxed" part of the design: a single outlier replicate cannot
   flip a gene's profile, because the median absorbs it. Z-score each gene's row
   across these 6 cells so clustering reflects *pattern shape*, not just absolute
   expression level (a highly expressed flat gene and a lowly expressed flat gene
   should cluster together).

3. Combined clustering (Step A)
   Hierarchical clustering (average linkage, correlation distance — i.e. cluster
   by *shape* of the profile) of all filtered genes on their 6-cell z-score
   profile. Number of clusters (k) is chosen data-drivenly by sweeping k and
   picking the silhouette-maximizing value within a sane range, then the
   dendrogram is cut at that k. Output: one big clustered heatmap (genes x
   6 genotype/stage cells) with a colour-coded cluster strip, similar styling
   to Plot_Heatmaps_Informative.py.

4. Cluster characterization
   Each cluster is auto-labeled using simple, transparent rules applied to its
   mean z-score profile (e.g. "BOTCv-biased / early-high", "BOTv-biased /
   late-rising", "shared temporal / concordant", "divergent / discordant").
   Within each cluster, genes are ranked by a genotype-contrast score
   (mean |z_BOTCv - z_BOTv| across stages) to surface anchor-gene candidates —
   the genes most responsible for that cluster's genotype-specific signature.

5. Temporal / trajectory view (Step B)
   For each cluster, plot the mean +/- IQR trajectory of BOTv vs BOTCv across
   nc14b -> nc14d -> gastr, so you can see *when* the two genotypes diverge or
   converge, not just whether they differ on average.

OUTPUT
------
  results/combined_figures/explore_clusters_botv_botcv/
    gene_profile_matrix.tsv             genes x 6 cells, z-scored medians (full filtered set)
    cluster_assignments.tsv             gene, cluster, contrast_score, auto_label
    anchor_gene_candidates.tsv          top-N contrast genes per cluster, ranked
    silhouette_sweep.tsv                k vs silhouette score (for transparency/QC)
    heatmap_all_clusters.pdf/png        Step A: full clustered heatmap, all filtered genes
    heatmap_cluster_<N>_<label>.pdf/png per-cluster zoom-in heatmap (genes x replicates, not just medians)
    trajectory_cluster_<N>_<label>.pdf/png  Step B: BOTv vs BOTCv mean trajectory per cluster
    trajectory_overview.pdf/png         all cluster trajectories on small multiples, for quick scanning

USAGE
-----
  python Analysis_GeneClusters_BOTvBOTCv.py
  python Analysis_GeneClusters_BOTvBOTCv.py --min-tpm 2 --k-max 20
  python Analysis_GeneClusters_BOTvBOTCv.py --groups FoxL1_BOTv HLH54F_BOTCv WT
  python Analysis_GeneClusters_BOTvBOTCv.py --top-anchors 15

NOTES
-----
- Restricted by default to FoxL1_BOTv (BOTv) vs HLH54F_BOTCv (BOTCv) only, per
  your request to focus first on this ventralized-genotype contrast. WT can be
  added with --groups for a 3-genotype version (e.g. to see which BOTv/BOTCv
  differences are also different from WT vs. specific to the mutant pair) but
  is NOT included by default, since WT lacks the toll allele and is not staged
  the same way as the two ventralized genotypes.
- This is intentionally exploratory / hypothesis-generating. Cluster identities
  and anchor gene calls should be treated as candidates for follow-up (e.g.
  in Plot_Heatmaps_Informative.py / Plot_PositionalClassification.py),
  not as a final result.
"""

import os
import sys
import argparse
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import silhouette_score

try:
    import Config_SampleMetadata as SC
except ImportError:
    sys.exit("[explore] Config_SampleMetadata.py not found.")
try:
    import Config_DataLoader as DL
except ImportError:
    sys.exit("[explore] Config_DataLoader.py not found.")

# =============================================================================
# CONFIG
# =============================================================================

OUT_DIR = os.path.join("results", "combined_figures", "explore_clusters_botv_botcv")
SYMBOL_MAP = "gene_id_to_symbol.tsv"
DEFAULT_PREFIX = "all_genotypes_plus_set3"

STAGE_ORDER_FULL = ["nc14b", "nc14d", "gastr", "late"]  # superset; "late" only
                                                          # appears for merged groups

# Which genotype groups have nc14d+gastr merged into "late". Set at runtime
# by main() from --merge-late-stages-for. Module-level so all helper
# functions (which call get_stage()) see the choice without it being
# threaded through every call.
#
# Rationale: FoxL1_BOTv and HLH54F_BOTCv often have very different replicate
# counts per stage (e.g. FoxL1_BOTv nc14b=3/nc14d=2/gastr=1 vs HLH54F_BOTCv
# nc14b=3/nc14d=1/gastr=1). Merging is only needed for whichever genotype is
# actually replicate-starved at nc14d/gastr — forcing the merge on both
# genotypes throws away real, well-replicated temporal resolution in the
# genotype that didn't need it. This is per-group rather than a single
# on/off switch for that reason.
MERGE_LATE_FOR: set[str] = set()


def get_stage(col: str, group: str) -> str:
    """
    Wrapper around SC.get_stage() that optionally merges nc14d + gastr into a
    single "late" bin, per genotype group (controlled by MERGE_LATE_FOR /
    --merge-late-stages-for).

    nc14d and gastr are developmentally close (minutes apart), so for a
    genotype with too few replicates at one or both of those stages (e.g.
    n=1), merging trades temporal resolution for replicate count in just
    that genotype, while leaving any well-replicated genotype at full
    3-stage resolution.
    """
    raw = SC.get_stage(col)
    if group in MERGE_LATE_FOR and raw in ("nc14d", "gastr"):
        return "late"
    return raw


def group_stage_order(group: str) -> list[str]:
    """Ordered stage list for one genotype group, respecting MERGE_LATE_FOR."""
    if group in MERGE_LATE_FOR:
        return ["nc14b", "late"]
    return ["nc14b", "nc14d", "gastr"]

DEFAULT_GROUPS = ["FoxL1_BOTv", "HLH54F_BOTCv"]

GROUP_COLOURS = {
    **SC.GROUP_COLOURS,
    "FoxL1_BOTv":   "#238b45",
    "HLH54F_BOTCv": "#e377c2",
}

# Same outlier convention as Plot_Heatmaps_Informative.py / Plot_PCA_UMAP.py
EXCLUDE_EMBRYOS: list[tuple[str, str]] = [
    ("FoxL1_BOTv", "bc27"),
]

CLUSTER_METHOD = "average"     # linkage method
CLUSTER_METRIC = "correlation"  # cluster by SHAPE of profile, not absolute level

# Cluster-labeling thresholds (applied to mean cluster z-score profile)
LABEL_BIAS_THRESH = 0.5     # |mean genotype difference| above this -> genotype-biased label
LABEL_TREND_THRESH = 0.4    # |late - early| above this -> "rising"/"falling" qualifier


def _is_excluded(sample_id: str, group: str) -> bool:
    for ex_group, ex_bc in EXCLUDE_EMBRYOS:
        if group == ex_group and ex_bc in sample_id:
            return True
    return False


# =============================================================================
# CLI
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--groups", nargs="+", default=DEFAULT_GROUPS,
                   help=f"Genotype groups to include (default: {DEFAULT_GROUPS}).")
    p.add_argument("--min-tpm", type=float, default=1.0,
                   help="Median-TPM threshold for the expression filter (default: 1.0).")
    p.add_argument("--min-stage-n", type=int, default=2,
                   help="Minimum replicates required in a (group, stage) cell for it "
                        "to count toward the filter / median profile (default: 2). "
                        "Cells with fewer replicates are set to NaN in the profile "
                        "matrix rather than trusted on n=1.")
    p.add_argument("--merge-late-stages-for", nargs="+", default=None, metavar="GROUP",
                   help="Merge nc14d and gastr into a single 'late' stage bin, "
                        "only for the named genotype group(s) (e.g. "
                        "--merge-late-stages-for HLH54F_BOTCv). Groups not "
                        "listed keep full nc14b/nc14d/gastr resolution. Useful "
                        "when one genotype's nc14d/gastr replicate counts are "
                        "too low (e.g. n=1) for a reliable per-stage median "
                        "while another genotype has enough replicates to keep "
                        "the full 3-stage trajectory — nc14d and gastr are "
                        "developmentally close together, so merging trades "
                        "temporal resolution for replicate count only in the "
                        "genotype(s) that need it. Off by default. Must be a "
                        "subset of --groups.")
    p.add_argument("--k-min", type=int, default=4)
    p.add_argument("--k-max", type=int, default=24)
    p.add_argument("--k-override", type=int, default=None,
                   help="Skip the silhouette sweep and force a specific k.")
    p.add_argument("--top-anchors", type=int, default=10,
                   help="Number of top contrast genes to report per cluster (default: 10).")
    p.add_argument("--outdir", default=OUT_DIR)
    p.add_argument("--prefix", default=DEFAULT_PREFIX,
                   help=f"Combined-matrix prefix to load, matching whatever "
                        f"--combined you used for step5b/step5c/step6 "
                        f"(default: {DEFAULT_PREFIX}). Determines both the "
                        f"RUVg-corrected matrix (via data_loader.py) and the "
                        f"wide pre-RUVg fallback TPM matrix.")
    return p.parse_args()


# =============================================================================
# Symbol map (optional)
# =============================================================================

def load_symbol_maps():
    symbol_to_fbgn, fbgn_to_symbol = {}, {}
    if os.path.exists(SYMBOL_MAP):
        sym_df = pd.read_csv(SYMBOL_MAP, sep="\t").dropna(subset=["Symbol", "GeneID"])
        sym_df = sym_df[sym_df["Symbol"].str.strip() != ""]
        symbol_to_fbgn = dict(zip(sym_df["Symbol"], sym_df["GeneID"]))
        fbgn_to_symbol = dict(zip(sym_df["GeneID"], sym_df["Symbol"]))
        print(f"[explore] Loaded {len(symbol_to_fbgn):,} gene symbol mappings")
    else:
        print("[explore] NOTE: SYMBOL_MAP not found — gene IDs will be reported as-is")
    return symbol_to_fbgn, fbgn_to_symbol


def label_for(gene_id, fbgn_to_symbol):
    return fbgn_to_symbol.get(gene_id, gene_id)


# =============================================================================
# Data loading
# =============================================================================

def load_data(groups, prefix):
    print(f"[explore] Loading data (prefix='{prefix}')...")
    data = DL.load(prefix=prefix, use_qc=True, use_quantile=True)
    tpm_for_z = data.tpm_qn if data.tpm_qn is not None else data.tpm
    group_of = data.group_of

    tpm_wide_path = f"results/combined/{prefix}_qc_tpm_matrix.tsv"
    if os.path.exists(tpm_wide_path):
        tpm_wide = pd.read_csv(tpm_wide_path, sep="\t", index_col=0)
        if "GeneSymbol" in tpm_wide.columns:
            tpm_wide = tpm_wide.drop(columns="GeneSymbol")
        print(f"[explore] Loaded wide TPM: {tpm_wide.shape[0]:,} genes x {tpm_wide.shape[1]} embryos")
    else:
        print(f"[explore] WARNING: {tpm_wide_path} not found — using RUVg matrix only "
              f"(some genes dropped by RUVg correction will be missing)")
        tpm_wide = tpm_for_z

    ordered_cols = DL.ordered_cols(data)
    ordered_cols = [c for c in ordered_cols if group_of.get(c) in groups]
    ordered_cols = [c for c in ordered_cols if not _is_excluded(c, group_of.get(c, ""))]

    print(f"[explore] {len(ordered_cols)} embryos across groups: {groups}")
    for g in groups:
        n = sum(1 for c in ordered_cols if group_of[c] == g)
        print(f"[explore]   {g}: {n} embryos")

    return tpm_for_z, tpm_wide, group_of, ordered_cols


def build_combined_tpm(tpm_for_z, tpm_wide, cols):
    """
    Build a single log2(TPM+1) matrix over `cols`, genes = union of both
    matrices' indices, preferring RUVg-corrected values and falling back to
    wide QC TPM for genes RUVg correction dropped (same logic as
    Plot_Heatmaps_Informative.py, generalized to all genes rather than a
    curated list).
    """
    cols_z = [c for c in cols if c in tpm_for_z.columns]
    cols_w = [c for c in cols if c in tpm_wide.columns]
    missing = set(cols) - set(cols_z) - set(cols_w)
    if missing:
        print(f"[explore] WARNING: {len(missing)} sample(s) not found in either "
              f"matrix, dropping: {missing}")

    genes_z = set(tpm_for_z.index)
    genes_w = set(tpm_wide.index)
    all_genes = sorted(genes_z | genes_w)

    # Start from wide matrix (broader gene coverage), then overwrite with
    # RUVg-corrected values wherever both gene and sample are available there.
    base = tpm_wide.reindex(index=all_genes, columns=cols)
    if not tpm_for_z.empty:
        z_sub = tpm_for_z.reindex(index=all_genes, columns=cols)
        base = base.combine_first(z_sub)  # fills NaN slots in `base` from z_sub
        # Prefer RUVg values where both present:
        common_idx = z_sub.dropna(how="all").index
        common_idx = common_idx.intersection(base.index)
        for col in cols:
            if col in z_sub.columns:
                vals = z_sub[col]
                mask = vals.notna()
                base.loc[mask.index[mask], col] = vals[mask]

    log_tpm = np.log2(base.astype(float) + 1)
    return log_tpm


# =============================================================================
# Filtering: median TPM > threshold within at least one (group, stage) cell
# =============================================================================

def group_stage_cells(cols, group_of):
    """Return {(group, stage): [cols]} for the given column list."""
    cells = {}
    for c in cols:
        g = group_of.get(c, "unknown")
        s = get_stage(c, g)
        cells.setdefault((g, s), []).append(c)
    return cells


def filter_expressed_genes(log_tpm, cells, min_tpm, min_stage_n):
    """
    Keep genes where median TPM (linear scale) > min_tpm in at least one
    (group, stage) cell that has >= min_stage_n replicates.
    """
    linear = (2 ** log_tpm) - 1
    keep_mask = pd.Series(False, index=log_tpm.index)
    cell_medians = {}

    for (g, s), cell_cols in cells.items():
        if len(cell_cols) < min_stage_n:
            continue
        med = linear[cell_cols].median(axis=1)
        cell_medians[(g, s)] = med
        keep_mask = keep_mask | (med > min_tpm)

    kept = log_tpm.index[keep_mask]
    print(f"[explore] Expression filter: median TPM > {min_tpm} in >=1 "
          f"(group,stage) cell with n>={min_stage_n} reps")
    print(f"[explore]   {len(kept):,} / {len(log_tpm):,} genes pass")
    return kept.tolist(), cell_medians


# =============================================================================
# Build genes x (group,stage) profile matrix, z-scored per gene
# =============================================================================

def build_profile_matrix(log_tpm, kept_genes, groups, cells, min_stage_n):
    """
    genes x (group,stage) matrix of median log2(TPM+1), restricted to kept_genes.

    (group,stage) cells with fewer than min_stage_n replicates are dropped
    from the matrix entirely rather than kept as an all-NaN column: cluster_genes()
    zero-fills before computing correlation distance, so an all-NaN column would
    act as a phantom "flat" dimension shared by every gene, diluting real signal
    across every pairwise comparison. Use --merge-late-stages-for to fold a
    thin stage into an adjacent one instead of dropping it, when appropriate.
    Returns (raw_profile, z_profile, cell_order, col_labels).
    """
    all_cell_order = [(g, s) for g in groups for s in group_stage_order(g) if (g, s) in cells]

    cell_order, dropped = [], []
    for (g, s) in all_cell_order:
        n = len(cells[(g, s)])
        if n < min_stage_n:
            dropped.append((g, s, n))
            continue
        cell_order.append((g, s))

    if dropped:
        print(f"[explore] Dropping {len(dropped)} (group,stage) cell(s) below "
              f"min_stage_n={min_stage_n} from the profile matrix entirely "
              f"(avoids an all-NaN column diluting every gene's correlation "
              f"profile with a phantom flat dimension):")
        for g, s, n in dropped:
            print(f"[explore]   {g}|{s}  (n={n} rep{'s' if n != 1 else ''}) — "
                  f"dropped. If {s} should instead be folded into an adjacent "
                  f"stage rather than dropped, add {g} to "
                  f"--merge-late-stages-for.")

    col_labels = [f"{g}|{s}" for g, s in cell_order]

    raw = pd.DataFrame(index=kept_genes, columns=col_labels, dtype=float)
    for (g, s), lbl in zip(cell_order, col_labels):
        cell_cols = cells[(g, s)]
        raw[lbl] = log_tpm.loc[kept_genes, cell_cols].median(axis=1)

    # Drop genes that ended up all-NaN (shouldn't happen given the filter, but safe)
    raw = raw.dropna(how="all")

    # Z-score each gene's row across available cells (ignoring NaN), so
    # clustering reflects pattern shape rather than absolute expression level.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        row_mean = raw.mean(axis=1, skipna=True)
        row_std = raw.std(axis=1, skipna=True)
    z = raw.subtract(row_mean, axis=0).divide(row_std.replace(0, np.nan), axis=0)
    # Genes with zero variance (flat across cells) -> all-zero z row rather than NaN
    flat_genes = row_std[row_std == 0].index
    z.loc[flat_genes] = 0.0

    return raw, z, cell_order, col_labels


# =============================================================================
# Clustering: silhouette-based k selection
# =============================================================================

def cluster_genes(z, k_min, k_max, k_override=None):
    """
    Hierarchical clustering of genes (rows of z) by profile shape.
    Returns (linkage_matrix, cluster_labels (1-indexed), chosen_k, sweep_df).
    Rows with any NaN are imputed to 0 (= "no deviation from gene mean" at that
    cell) purely for the distance computation. Since build_profile_matrix()
    now drops any (group,stage) column with fewer than min_stage_n replicates
    entirely (rather than keeping it as an all-NaN column), this should only
    ever affect isolated single-gene dropouts, not whole columns.

    The pairwise correlation distance matrix is computed and validated
    explicitly here (rather than letting linkage(X, metric=...) compute it
    internally), because with few columns (e.g. 3, after dropping
    under-replicated cells) floating-point noise can produce a tiny negative
    distance or an occasional non-finite value that would otherwise crash
    deep inside scipy. Validating it here lets us clip/repair it explicitly
    and report how many pairs were affected.
    """
    z_filled = z.fillna(0.0)
    nan_frac = z.isna().values.mean()
    if nan_frac > 0.01:
        print(f"[explore] WARNING: {nan_frac:.1%} of the z-score matrix is "
              f"NaN going into clustering (being zero-filled). With "
              f"build_profile_matrix() dropping under-replicated columns "
              f"entirely, this should normally be near 0% — if it's not, "
              f"something is generating widespread per-gene NaNs that "
              f"deserves a look before trusting the clustering.")
    X = z_filled.values.copy()   # .copy() avoids a read-only view on some
                                  # pandas/numpy versions

    # Correlation distance is undefined for zero-variance (flat) rows.
    # X is already z-scored per-row, so a flat row should be exactly zero —
    # but use a small tolerance rather than an exact ==0 comparison, since
    # with few columns (e.g. 3), floating-point residue from the earlier
    # divide-by-std step can leave a "flat" row at ~1e-15 instead of
    # literally 0. Jitter both cases so pdist never sees a zero-norm vector.
    row_norms = np.linalg.norm(X, axis=1)
    flat_idx = np.where(row_norms < 1e-8)[0]
    if len(flat_idx) > 0:
        print(f"[explore] {len(flat_idx)} gene(s) have a flat (zero-variance) "
              f"profile — adding negligible jitter so correlation distance "
              f"is defined for them (they'll cluster together, near-randomly, "
              f"which is appropriate since they have no real shape signal).")
        rng = np.random.default_rng(0)
        X[flat_idx] += rng.normal(scale=1e-6, size=(len(flat_idx), X.shape[1]))

    print(f"[explore] Computing pairwise {CLUSTER_METRIC} distances "
          f"({X.shape[0]:,} genes)...")
    dist_condensed = pdist(X, metric=CLUSTER_METRIC)

    # Correlation distance is mathematically bounded to [0, 2]. Clip tiny
    # floating-point overshoot (e.g. -1e-16 from a correlation rounding to
    # just over 1.0), and repair any genuinely non-finite leftovers by
    # treating them as maximally dissimilar (2.0) rather than crashing —
    # an undefined correlation between two near-degenerate vectors carries
    # no real shape information anyway.
    n_negative = int((dist_condensed < 0).sum())
    n_nonfinite = int((~np.isfinite(dist_condensed)).sum())
    if n_negative:
        print(f"[explore] clipping {n_negative} distance(s) with tiny "
              f"floating-point overshoot below 0.")
    if n_nonfinite:
        print(f"[explore] WARNING: {n_nonfinite} pairwise distance(s) were "
              f"non-finite (undefined correlation) — setting to the maximum "
              f"distance (2.0, 'maximally dissimilar') rather than crashing. "
              f"Worth spot-checking anchor_gene_candidates.tsv near these "
              f"genes if this number is large.")
        dist_condensed = np.where(np.isfinite(dist_condensed),
                                   dist_condensed, 2.0)
    dist_condensed = np.clip(dist_condensed, 0.0, 2.0)

    Z_link = linkage(dist_condensed, method=CLUSTER_METHOD)

    n = X.shape[0]
    k_max_eff = min(k_max, n - 1)
    if k_max_eff < k_min:
        k_min, k_max_eff = 2, max(2, min(k_max, n - 1))

    if k_override is not None:
        labels = fcluster(Z_link, t=k_override, criterion="maxclust")
        return Z_link, labels, k_override, None

    print(f"[explore] Sweeping k = {k_min}..{k_max_eff} for silhouette score "
          f"(distance metric: {CLUSTER_METRIC})...")
    dist_sq = squareform(dist_condensed)

    sweep_rows = []
    best_k, best_score = None, -np.inf
    for k in range(k_min, k_max_eff + 1):
        labels_k = fcluster(Z_link, t=k, criterion="maxclust")
        if len(set(labels_k)) < 2:
            continue
        try:
            score = silhouette_score(dist_sq, labels_k, metric="precomputed")
        except Exception:
            continue
        sweep_rows.append({"k": k, "silhouette": score})
        if score > best_score:
            best_k, best_score = k, score

    sweep_df = pd.DataFrame(sweep_rows)
    if best_k is None:
        print("[explore] WARNING: silhouette sweep failed; defaulting to k=8")
        best_k = min(8, k_max_eff)

    print(f"[explore] Selected k={best_k} (silhouette={best_score:.3f})")
    labels = fcluster(Z_link, t=best_k, criterion="maxclust")
    return Z_link, labels, best_k, sweep_df


# =============================================================================
# Cluster characterization / auto-labeling
# =============================================================================

def _canonical_late_value(mean_profile, group, stages_present):
    """
    Value to use for group's 'late' timepoint on a canonical (nc14b, late)
    comparison axis, regardless of whether that group was actually merged.

    - If the group has a literal 'late' column (it was merged), use it directly.
    - Otherwise, average whichever of nc14d/gastr are present, so a
      non-merged group can still be compared against a merged group on the
      same early-vs-late axis without losing its own full-resolution data
      elsewhere (the profile matrix / heatmap / per-cluster plots still keep
      nc14d and gastr separate for that group).
    """
    key_late = f"{group}|late"
    if key_late in mean_profile.index:
        return mean_profile[key_late]
    vals = [mean_profile[f"{group}|{s}"] for s in ("nc14d", "gastr")
            if f"{group}|{s}" in mean_profile.index]
    return float(np.mean(vals)) if vals else np.nan


def characterize_clusters(z, labels, cell_order, groups):
    """
    For each cluster, compute mean z-profile and an auto-generated descriptive
    label based on:
      - genotype bias: mean(group_A) - mean(group_B) at nc14b and at a
        canonical "late" timepoint (literal 'late' if the group was merged,
        else the mean of whichever of nc14d/gastr it has — see
        _canonical_late_value). This lets two genotypes be compared on a
        common 2-point timeline even when one is merged and the other isn't.
      - temporal trend: canonical late value - nc14b value, per genotype.
    Returns dict: cluster_id -> {label, mean_profile (Series), n_genes}
    """
    info = {}
    g_a, g_b = groups[0], groups[1] if len(groups) > 1 else (groups[0], None)

    for cl in sorted(set(labels)):
        idx = z.index[labels == cl]
        mean_profile = z.loc[idx].mean(axis=0)

        parts = []
        if g_b is not None:
            stage_pairs = []  # (label, value_a, value_b)
            if f"{g_a}|nc14b" in mean_profile.index and f"{g_b}|nc14b" in mean_profile.index:
                stage_pairs.append(("nc14b", mean_profile[f"{g_a}|nc14b"], mean_profile[f"{g_b}|nc14b"]))
            late_a = _canonical_late_value(mean_profile, g_a, cell_order)
            late_b = _canonical_late_value(mean_profile, g_b, cell_order)
            if not np.isnan(late_a) and not np.isnan(late_b):
                stage_pairs.append(("late", late_a, late_b))

            if stage_pairs:
                diffs = [va - vb for _, va, vb in stage_pairs]
                mean_diff = float(np.mean(diffs))
                if mean_diff > LABEL_BIAS_THRESH:
                    parts.append(f"{g_a}-high")
                elif mean_diff < -LABEL_BIAS_THRESH:
                    parts.append(f"{g_b}-high")
                else:
                    parts.append("concordant")

                # Does the bias grow, shrink, or flip from nc14b -> late?
                if len(diffs) >= 2:
                    if (diffs[-1] - diffs[0]) > LABEL_TREND_THRESH and abs(diffs[0]) < abs(diffs[-1]):
                        parts.append("diverging")
                    elif (diffs[0] - diffs[-1]) > LABEL_TREND_THRESH and abs(diffs[-1]) < abs(diffs[0]):
                        parts.append("converging")
                    elif np.sign(diffs[0] + 1e-9) != np.sign(diffs[-1] + 1e-9) and \
                         min(abs(diffs[0]), abs(diffs[-1])) > LABEL_BIAS_THRESH * 0.5:
                        parts.append("crossover")

        # Temporal trend pooled across genotypes present (nc14b -> canonical late)
        early_vals, late_vals = [], []
        for g in groups:
            if f"{g}|nc14b" in mean_profile.index:
                early_vals.append(mean_profile[f"{g}|nc14b"])
            lv = _canonical_late_value(mean_profile, g, cell_order)
            if not np.isnan(lv):
                late_vals.append(lv)
        if early_vals and late_vals:
            trend = float(np.mean(late_vals)) - float(np.mean(early_vals))
            if trend > LABEL_TREND_THRESH:
                parts.append("rising")
            elif trend < -LABEL_TREND_THRESH:
                parts.append("falling")
            else:
                parts.append("stable")

        auto_label = "_".join(parts) if parts else "uncategorized"
        info[cl] = {
            "label": auto_label,
            "mean_profile": mean_profile,
            "n_genes": int((labels == cl).sum()),
        }
    return info


def compute_contrast_scores(z, cell_order, groups):
    """
    Per-gene genotype-contrast score: mean |z_A - z_B| at nc14b and at a
    canonical "late" timepoint (see _canonical_late_value), across both
    genotypes' available stages. Higher = bigger, more consistent difference
    between the two genotypes -> stronger anchor-gene candidate. Works the
    same whether both genotypes are at full 3-stage resolution, both merged,
    or only one merged.
    """
    if len(groups) < 2:
        return pd.Series(0.0, index=z.index)
    g_a, g_b = groups[0], groups[1]

    diff_cols = []
    if f"{g_a}|nc14b" in z.columns and f"{g_b}|nc14b" in z.columns:
        diff_cols.append((z[f"{g_a}|nc14b"] - z[f"{g_b}|nc14b"]).abs())

    def _late_series(group):
        key_late = f"{group}|late"
        if key_late in z.columns:
            return z[key_late]
        present = [s for s in ("nc14d", "gastr") if f"{group}|{s}" in z.columns]
        if not present:
            return None
        return pd.concat([z[f"{group}|{s}"] for s in present], axis=1).mean(axis=1)

    late_a, late_b = _late_series(g_a), _late_series(g_b)
    if late_a is not None and late_b is not None:
        diff_cols.append((late_a - late_b).abs())

    if not diff_cols:
        return pd.Series(0.0, index=z.index)
    return pd.concat(diff_cols, axis=1).mean(axis=1)


# =============================================================================
# Plotting helpers
# =============================================================================

def save_fig(fig, outdir, name):
    png = os.path.join(outdir, f"{name}.png")
    pdf = os.path.join(outdir, f"{name}.pdf")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"[explore] Saved {name}")


CLUSTER_CMAP_NAME = "tab20"


def cluster_colour_map(cluster_ids):
    cmap = plt.get_cmap(CLUSTER_CMAP_NAME)
    ids = sorted(cluster_ids)
    return {cl: cmap(i % 20) for i, cl in enumerate(ids)}


def plot_full_clustered_heatmap(z, labels, Z_link, cell_order, col_labels,
                                 cluster_info, fbgn_to_symbol, outdir):
    order = dendrogram(Z_link, no_plot=True)["leaves"]
    z_ord = z.iloc[order]
    labels_ord = labels[order]
    gene_names = [label_for(g, fbgn_to_symbol) for g in z_ord.index]

    n_genes, n_cells = z_ord.shape
    fig_h = max(8, min(40, n_genes * 0.045))
    fig_w = max(10, n_cells * 1.1 + 6)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.16, right=0.72, top=0.92, bottom=0.08)

    im = ax.imshow(z_ord.values, aspect="auto", cmap="RdBu_r",
                   vmin=-2.5, vmax=2.5, interpolation="nearest")

    ax.set_xticks(range(n_cells))
    ax.set_xticklabels(col_labels, rotation=45, ha="right", fontsize=9)

    # Suppress individual y labels if too many genes; show every Nth
    if n_genes <= 120:
        ax.set_yticks(range(n_genes))
        ax.set_yticklabels(gene_names, fontsize=5)
    else:
        ax.set_yticks([])
        ax.set_ylabel(f"{n_genes:,} genes (labels hidden — see per-cluster heatmaps)",
                      fontsize=9)

    # Cluster colour strip on the left
    cl_colours = cluster_colour_map(set(labels))
    strip_ax = ax.inset_axes([-0.05, 0, 0.02, 1], transform=ax.transAxes)
    strip_ax.set_xlim(0, 1)
    strip_ax.set_ylim(0, n_genes)
    strip_ax.axis("off")
    for i, cl in enumerate(labels_ord):
        strip_ax.add_patch(plt.Rectangle(
            (0, n_genes - i - 1), 1, 1,
            color=cl_colours[cl], transform=strip_ax.transData))

    # Cluster boundary lines + labels
    boundaries = []
    pos = 0
    label_positions = []
    for cl in pd.unique(labels_ord):
        n_in = int((labels_ord == cl).sum())
        label_positions.append((pos + n_in / 2, cl))
        pos += n_in
        if pos < n_genes:
            boundaries.append(pos - 0.5)
    for b in boundaries:
        ax.axhline(b, color="black", lw=0.6, alpha=0.5)

    label_ax = ax.inset_axes([1.02, 0, 0.45, 1], transform=ax.transAxes)
    label_ax.set_xlim(0, 1)
    label_ax.set_ylim(0, n_genes)
    label_ax.axis("off")
    for ypos, cl in label_positions:
        lbl = cluster_info[cl]["label"]
        n_g = cluster_info[cl]["n_genes"]
        label_ax.text(0, n_genes - ypos, f"C{cl}: {lbl} (n={n_g})",
                      fontsize=6.5, va="center", ha="left",
                      color=cl_colours[cl])

    cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
    cbar.set_label("Z-score (within-gene)", fontsize=9)

    ax.set_title(
        f"Exploratory gene clustering — {n_genes:,} genes x "
        f"{n_cells} (genotype, stage) cells\n"
        f"median log2(TPM+1), z-scored per gene  |  {CLUSTER_METHOD} linkage, "
        f"{CLUSTER_METRIC} distance  |  k={len(set(labels))}",
        fontsize=10)

    save_fig(fig, outdir, "heatmap_all_clusters")


def plot_per_cluster_heatmap(cl, gene_ids, log_tpm, cols, group_of, fbgn_to_symbol,
                             cluster_info, outdir, suffix=""):
    """Zoom-in heatmap for one cluster at full replicate resolution (not just medians)."""
    sub = log_tpm.loc[gene_ids, cols]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        zsub = sub.subtract(sub.mean(axis=1), axis=0).divide(sub.std(axis=1), axis=0)
    zsub = zsub.fillna(0)

    # order columns by group then stage then sample id (stable, readable)
    def sort_key(c):
        g = group_of.get(c, "")
        s = get_stage(c, g)
        return (g, STAGE_ORDER_FULL.index(s) if s in STAGE_ORDER_FULL else 99, c)
    ordered_cols_local = sorted(cols, key=sort_key)
    zsub = zsub[ordered_cols_local]

    gene_names = [label_for(g, fbgn_to_symbol) for g in gene_ids]
    n_genes, n_emb = zsub.shape

    # Guard against runaway figure size for large clusters: cap height and
    # drop per-gene text labels above a threshold (full gene list with
    # contrast scores is always available in cluster_assignments.tsv /
    # anchor_gene_candidates.tsv — this plot is for visual pattern-checking).
    MAX_GENES_LABELED = 80
    MAX_FIG_H = 30
    show_labels = n_genes <= MAX_GENES_LABELED

    fig_w = max(8, n_emb * 0.35)
    fig_h = min(MAX_FIG_H, max(4, n_genes * (0.22 if show_labels else 0.05)))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=0.22 if show_labels else 0.08, right=0.85, top=0.85, bottom=0.18)

    im = ax.imshow(zsub.values, aspect="auto", cmap="RdBu_r",
                   vmin=-2.5, vmax=2.5, interpolation="nearest")
    ax.set_xticks([])
    if show_labels:
        ax.set_yticks(range(n_genes))
        ax.set_yticklabels(gene_names, fontstyle="italic", fontsize=8)
    else:
        ax.set_yticks([])
        ax.set_ylabel(f"{n_genes:,} genes (labels hidden — see "
                      f"cluster_assignments.tsv for full list)", fontsize=8)

    # group/stage strip
    pos = 0
    boundaries = []
    label_positions = []
    for c in ordered_cols_local:
        pass
    seen_groups = []
    for c in ordered_cols_local:
        g = group_of.get(c, "")
        if g not in seen_groups:
            seen_groups.append(g)
    pos = 0
    for g in seen_groups:
        g_cols = [c for c in ordered_cols_local if group_of.get(c) == g]
        if pos > 0:
            boundaries.append(pos - 0.5)
        label_positions.append((pos + len(g_cols) / 2, g))
        pos += len(g_cols)
    for b in boundaries:
        ax.axvline(b, color="black", lw=1.5)

    strip_ax = ax.inset_axes([0, 1.01, 1, 0.04], transform=ax.transAxes)
    strip_ax.set_xlim(0, n_emb)
    strip_ax.set_ylim(0, 1)
    strip_ax.axis("off")
    for i, c in enumerate(ordered_cols_local):
        strip_ax.add_patch(plt.Rectangle(
            (i, 0), 1, 1, color=GROUP_COLOURS.get(group_of.get(c), "#888"),
            transform=strip_ax.transData))
    for xpos, g in label_positions:
        strip_ax.text(xpos, 1.3, g, ha="center", va="bottom", fontsize=7,
                      fontweight="bold", color=GROUP_COLOURS.get(g, "black"),
                      transform=strip_ax.transData)

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Z-score", fontsize=8)

    lbl = cluster_info[cl]["label"]
    title_extra = " (top anchor genes)" if suffix else ""
    ax.set_title(f"Cluster {cl}: {lbl}{title_extra}  (n={n_genes} genes shown, replicate-level)",
                fontsize=10, pad=22)

    save_fig(fig, outdir, f"heatmap_cluster_{cl}_{lbl}{suffix}")


# Canonical x-axis positions so genotypes with different stage resolutions
# (e.g. one merged to nc14b/late, one at full nc14b/nc14d/gastr) can be
# plotted on the same trajectory axis. "late" is placed at the midpoint of
# where nc14d/gastr would fall, so a merged genotype's single late point
# sits visually between a non-merged genotype's nc14d and gastr points.
STAGE_X_POSITION = {"nc14b": 0, "nc14d": 1, "gastr": 2, "late": 1.5}
STAGE_X_TICKS = [0, 1, 1.5, 2]
STAGE_X_TICKLABELS = ["nc14b", "nc14d", "late", "gastr"]


def plot_cluster_trajectory(cl, mean_profile, groups, cluster_info, outdir):
    lbl = cluster_info[cl]["label"]
    n_g = cluster_info[cl]["n_genes"]
    fig, ax = plt.subplots(figsize=(5, 4))
    any_late = False
    for g in groups:
        xs, ys = [], []
        for s in group_stage_order(g):
            key = f"{g}|{s}"
            if key in mean_profile.index and not np.isnan(mean_profile[key]):
                xs.append(STAGE_X_POSITION[s])
                ys.append(mean_profile[key])
                if s == "late":
                    any_late = True
        if xs:
            ax.plot(xs, ys, marker="o", lw=2.2, ms=7,
                   color=GROUP_COLOURS.get(g, "black"), label=g)
    if any_late:
        ax.set_xticks(STAGE_X_TICKS)
        ax.set_xticklabels(STAGE_X_TICKLABELS)
    else:
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["nc14b", "nc14d", "gastr"])
    ax.axhline(0, color="grey", lw=0.8, ls=":")
    ax.set_ylabel("Mean z-score (cluster)")
    ax.set_title(f"Cluster {cl}: {lbl}\n(n={n_g} genes)", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save_fig(fig, outdir, f"trajectory_cluster_{cl}_{lbl}")


def plot_trajectory_overview(cluster_info, groups, outdir):
    cl_ids = sorted(cluster_info.keys())
    ncols = min(4, len(cl_ids))
    nrows = (len(cl_ids) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 2.6 * nrows),
                              squeeze=False)
    axes_flat = axes.flatten()

    for ax, cl in zip(axes_flat, cl_ids):
        mp = cluster_info[cl]["mean_profile"]
        any_late = False
        for g in groups:
            xs, ys = [], []
            for s in group_stage_order(g):
                key = f"{g}|{s}"
                if key in mp.index and not np.isnan(mp[key]):
                    xs.append(STAGE_X_POSITION[s])
                    ys.append(mp[key])
                    if s == "late":
                        any_late = True
            if xs:
                ax.plot(xs, ys, marker="o", ms=4, lw=1.6,
                       color=GROUP_COLOURS.get(g, "black"))
        ax.axhline(0, color="grey", lw=0.6, ls=":")
        if any_late:
            ax.set_xticks(STAGE_X_TICKS)
            ax.set_xticklabels(STAGE_X_TICKLABELS, fontsize=6)
        else:
            ax.set_xticks([0, 1, 2])
            ax.set_xticklabels(["nc14b", "nc14d", "gastr"], fontsize=6)
        ax.set_title(f"C{cl}: {cluster_info[cl]['label']}\n(n={cluster_info[cl]['n_genes']})",
                    fontsize=7)
        ax.tick_params(labelsize=6)
        ax.spines[["top", "right"]].set_visible(False)

    for ax in axes_flat[len(cl_ids):]:
        ax.set_visible(False)

    handles = [mpatches.Patch(color=GROUP_COLOURS.get(g, "black"), label=g) for g in groups]
    fig.legend(handles=handles, loc="upper center", ncol=len(groups),
              frameon=False, fontsize=9, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("All cluster trajectories — quick scan", fontsize=11, y=1.06)
    fig.tight_layout()
    save_fig(fig, outdir, "trajectory_overview")


# =============================================================================
# Main
# =============================================================================

def main():
    global MERGE_LATE_FOR

    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    MERGE_LATE_FOR = set(args.merge_late_stages_for or [])
    bad_groups = MERGE_LATE_FOR - set(args.groups)
    if bad_groups:
        sys.exit(f"[explore] --merge-late-stages-for group(s) {bad_groups} "
                 f"not in --groups {args.groups}.")
    if MERGE_LATE_FOR:
        print(f"[explore] --merge-late-stages-for active for: {sorted(MERGE_LATE_FOR)} "
              f"(nc14d + gastr -> 'late' for these groups only; other groups "
              f"keep full nc14b/nc14d/gastr resolution)")

    symbol_to_fbgn, fbgn_to_symbol = load_symbol_maps()
    tpm_for_z, tpm_wide, group_of, cols = load_data(args.groups, args.prefix)

    if len(cols) == 0:
        sys.exit("[explore] No embryos found for the requested groups. Check --groups.")

    log_tpm = build_combined_tpm(tpm_for_z, tpm_wide, cols)
    cells = group_stage_cells(cols, group_of)

    print("[explore] (group, stage) cell sizes:")
    for (g, s), c in sorted(cells.items()):
        flag = "" if len(c) >= args.min_stage_n else "  <-- below min_stage_n, cell DROPPED from profile matrix"
        print(f"[explore]   {g:14s} {s:6s} n={len(c)}{flag}")

    kept_genes, _ = filter_expressed_genes(log_tpm, cells, args.min_tpm, args.min_stage_n)
    if len(kept_genes) < 10:
        sys.exit(f"[explore] Only {len(kept_genes)} genes passed the filter — "
                 f"too few to cluster meaningfully. Try lowering --min-tpm.")

    raw_profile, z_profile, cell_order, col_labels = build_profile_matrix(
        log_tpm, kept_genes, args.groups, cells, args.min_stage_n)

    print(f"[explore] Profile matrix: {z_profile.shape[0]:,} genes x "
          f"{z_profile.shape[1]} (group,stage) cells")

    raw_profile.to_csv(os.path.join(args.outdir, "gene_profile_matrix_raw_log2tpm.tsv"), sep="\t")
    z_profile.to_csv(os.path.join(args.outdir, "gene_profile_matrix_zscored.tsv"), sep="\t")

    Z_link, labels, chosen_k, sweep_df = cluster_genes(
        z_profile, args.k_min, args.k_max, k_override=args.k_override)
    if sweep_df is not None:
        sweep_df.to_csv(os.path.join(args.outdir, "silhouette_sweep.tsv"), sep="\t", index=False)

    cluster_info = characterize_clusters(z_profile, labels, cell_order, args.groups)
    contrast_scores = compute_contrast_scores(z_profile, cell_order, args.groups)

    # ── Write cluster assignment table ───────────────────────────────────────
    assign_df = pd.DataFrame({
        "gene_id": z_profile.index,
        "symbol": [label_for(g, fbgn_to_symbol) for g in z_profile.index],
        "cluster": labels,
        "contrast_score": contrast_scores.reindex(z_profile.index).values,
    })
    assign_df["cluster_label"] = assign_df["cluster"].map(lambda c: cluster_info[c]["label"])
    assign_df = assign_df.sort_values(["cluster", "contrast_score"], ascending=[True, False])
    assign_df.to_csv(os.path.join(args.outdir, "cluster_assignments.tsv"), sep="\t", index=False)

    # ── Anchor gene candidates: top-N by contrast score per cluster ────────────
    anchor_rows = []
    for cl in sorted(set(labels)):
        sub = assign_df[assign_df["cluster"] == cl].head(args.top_anchors)
        anchor_rows.append(sub)
    anchor_df = pd.concat(anchor_rows, ignore_index=True)
    anchor_df.to_csv(os.path.join(args.outdir, "anchor_gene_candidates.tsv"), sep="\t", index=False)

    print(f"\n[explore] === Cluster summary ===")
    for cl in sorted(set(labels)):
        info = cluster_info[cl]
        top_genes = assign_df[assign_df["cluster"] == cl]["symbol"].head(5).tolist()
        print(f"[explore]   C{cl:>2} ({info['n_genes']:>4} genes) {info['label']:30s} "
              f"top: {', '.join(top_genes)}")

    # ── Plots ─────────────────────────────────────────────────────────────────
    plot_full_clustered_heatmap(z_profile, labels, Z_link, cell_order, col_labels,
                                cluster_info, fbgn_to_symbol, args.outdir)

    for cl in sorted(set(labels)):
        gene_ids = z_profile.index[labels == cl].tolist()
        plot_per_cluster_heatmap(cl, gene_ids, log_tpm, cols, group_of,
                                 fbgn_to_symbol, cluster_info, args.outdir)
        # For large clusters the full heatmap hides gene labels (see
        # plot_per_cluster_heatmap) — also emit a small, labeled heatmap of
        # just the top anchor candidates for that cluster so there's always
        # at least one readable, gene-labeled view per cluster.
        if len(gene_ids) > 80:
            top_ids = (assign_df[assign_df["cluster"] == cl]
                      .head(args.top_anchors)["gene_id"].tolist())
            plot_per_cluster_heatmap(cl, top_ids, log_tpm, cols, group_of,
                                     fbgn_to_symbol, cluster_info, args.outdir,
                                     suffix="_top_anchors")
        plot_cluster_trajectory(cl, cluster_info[cl]["mean_profile"], args.groups,
                                cluster_info, args.outdir)

    plot_trajectory_overview(cluster_info, args.groups, args.outdir)

    print(f"\n[explore] Done. {len(set(labels))} clusters from {len(kept_genes):,} "
          f"filtered genes. Outputs in {args.outdir}/")
    print("[explore] Suggested next step: open anchor_gene_candidates.tsv and "
          "cross-check top genes against known biology / your curated gene "
          "lists, then spot-check candidates with heatmap_cluster_<N>_*.png "
          "before adding any to Plot_Heatmaps_Informative.py's GENES_BY_CATEGORY.")


if __name__ == "__main__":
    main()
