#!/usr/bin/env python3
"""
plot_pca_split_v5.py  —  PCA + UMAP plots for seRNAseq pipeline

Subsets (matching heatmap groupings):
  • all_no_runt      — WT + all mutants except runt
  • ventralized      — WT + FoxL1_BOTv + HLH54F_BOTCv
  • non_ventralized  — WT + everything except FoxL1_BOTv + HLH54F_BOTCv + runt

For each subset, produces both PCA and UMAP in two ellipse variants:
  • with_ellipses  — 95% confidence ellipses (shaded)
  • no_ellipses    — points only

Axes limits are ALWAYS set to the with-ellipses extent so plot size
is identical between variants (more zoomed-out wins).

Outputs: PDF (+ PNG at 300 dpi) for every plot.

v4 additions
------------
--label-qc   Generate an additional labeled QC pass for every subset
             (PCA and UMAP) with:
               • every embryo labeled with group + barcode
               • larger font, dark outline for legibility on any background
               • no ellipses (cleaner canvas for label reading)
               • slightly larger points (s=100) and reduced alpha on points
                 so overlapping labels remain distinct
             Output files are suffixed _qc_labeled (PDF + PNG).
             Compatible with --label-embryos (both sets are produced).

v5 additions
------------
EXCLUDE_EMBRYOS  Hard-coded set of (group, barcode) pairs to drop before any
                 dimensionality reduction or plotting.
                 Edit EXCLUDE_EMBRYOS below to adjust.
"""

import argparse
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.stats import chi2

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from Config_SampleMetadata import (
    GROUP_COLOURS,
    GROUP_ORDER,
    STAGE_ALPHA,
    build_group_map,
    get_stage,
)

# ── Override colours to match heatmap script ─────────────────────────────────
GROUP_COLOURS = {
    **GROUP_COLOURS,
    "FoxL1_BOTv":   "#238b45",
    "HLH54F_BOTCv": "#e377c2",
}

# ── Outlier exclusions ─────────────────────────────────────────────────────────
# Each entry is a (group, barcode_substring) pair.  A sample is excluded if its
# group matches AND its column ID contains the barcode string (e.g. "bc27").
# Edit this list to add or remove outliers — it is applied once at data load.
EXCLUDE_EMBRYOS: list[tuple[str, str]] = [
    ("FoxL1_BOTv", "bc27"),
]

def _is_excluded(sample_id: str, group: str) -> bool:
    for ex_group, ex_bc in EXCLUDE_EMBRYOS:
        if ex_bc not in sample_id:
            continue
        if group == ex_group:
            return True
        # Barcode matches but group name doesn't — likely a naming mismatch
        # between build_group_map and EXCLUDE_EMBRYOS. Warn and still exclude,
        # since barcodes are unique within the sample set.
        print(f"[plot] WARNING: {sample_id} matches barcode '{ex_bc}' but "
              f"has group='{group}' (expected '{ex_group}'). "
              f"Excluding anyway — check build_group_map naming.")
        return True
    return False

# ── Subset definitions (mirror heatmap code) ─────────────────────────────────
VENTRALIZED_GROUPS    = {"FoxL1_BOTv", "HLH54F_BOTCv"}

SUBSETS = {
    "all_no_runt": dict(
        label="WT + all mutants (excl. runt)",
        keep=lambda g: g != "runt",
    ),
    "all_no_runt_no_wt": dict(
        label="All mutants (excl. runt, excl. WT)",
        keep=lambda g: g != "runt" and g != "WT",
    ),
    "ventralized": dict(
        label="WT + ventralized (FoxL1_BOTv, HLH54F_BOTCv)",
        keep=lambda g: g in VENTRALIZED_GROUPS or g == "WT",
    ),
    "ventralized_no_wt": dict(
        label="Ventralized only (FoxL1_BOTv, HLH54F_BOTCv)",
        keep=lambda g: g in VENTRALIZED_GROUPS,
    ),
    "non_ventralized": dict(
        label="WT + non-ventralized (excl. runt, excl. vent.)",
        keep=lambda g: g not in VENTRALIZED_GROUPS and g != "runt",
    ),
    "non_ventralized_no_wt": dict(
        label="Non-ventralized only (excl. runt, excl. WT, excl. vent.)",
        keep=lambda g: g not in VENTRALIZED_GROUPS and g != "runt" and g != "WT",
    ),
}

STAGE_MARKERS = {
    "nc14b": "o",
    "nc14d": "s",
    "gastr": "^",
}

# =============================================================================
# Confidence ellipse
# =============================================================================

def confidence_ellipse(x, y, ax, n_std=2.0, facecolor="none", **kwargs):
    if len(x) < 3:
        return None
    cov = np.cov(x, y)
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    order = eigenvalues.argsort()[::-1]
    eigenvalues  = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    angle  = np.degrees(np.arctan2(*eigenvectors[:, 0][::-1]))
    scale  = np.sqrt(chi2.ppf(0.95, df=2))
    width, height = 2 * scale * np.sqrt(eigenvalues)
    ellipse = Ellipse(
        xy=(np.mean(x), np.mean(y)),
        width=width, height=height, angle=angle,
        facecolor=facecolor, **kwargs,
    )
    ax.add_patch(ellipse)
    return ellipse

# =============================================================================
# Data loading
# =============================================================================

def load_tpm(prefix: str, results_dir: str = "results/combined"):
    candidates = [
        os.path.join(results_dir, f"{prefix}_qc_tpm_quantile_matrix.tsv"),
        os.path.join(results_dir, f"{prefix}_qc_tpm_matrix.tsv"),
        os.path.join(results_dir, f"{prefix}_tpm_matrix.tsv"),
        f"{prefix}_qc_tpm_quantile_matrix.tsv",
        f"{prefix}_qc_tpm_matrix.tsv",
        f"{prefix}_tpm_matrix.tsv",
    ]
    for path in candidates:
        if os.path.exists(path):
            print(f"[plot] Loading TPM from: {path}")
            df = pd.read_csv(path, sep="\t", index_col=0)
            if "GeneSymbol" in df.columns:
                df = df.drop(columns="GeneSymbol")
            return np.log2(df + 1)
    raise FileNotFoundError(f"Cannot find TPM matrix for prefix '{prefix}'.")

# =============================================================================
# Dimensionality reduction
# =============================================================================

def run_pca(log_tpm, n_top=2000, n_components=5):
    var      = log_tpm.var(axis=1)
    top_genes = var.nlargest(min(n_top, len(var))).index
    X  = log_tpm.loc[top_genes].T
    Xs = StandardScaler().fit_transform(X)
    n_comp = min(n_components, Xs.shape[0], Xs.shape[1])
    pca    = PCA(n_components=n_comp)
    coords = pca.fit_transform(Xs)
    df_coords = pd.DataFrame(
        coords, index=X.index,
        columns=[f"PC{i+1}" for i in range(n_comp)],
    )
    return df_coords, pca


def run_umap(log_tpm, n_top=2000, n_neighbors=15, min_dist=0.3, random_state=42):
    try:
        import umap as umap_lib
    except ImportError:
        print("[plot] WARNING: umap-learn not installed — skipping UMAP. "
              "Install with: pip install umap-learn")
        return None

    var       = log_tpm.var(axis=1)
    top_genes = var.nlargest(min(n_top, len(var))).index
    X  = log_tpm.loc[top_genes].T
    Xs = StandardScaler().fit_transform(X)

    n_neighbors = min(n_neighbors, Xs.shape[0] - 1)
    reducer = umap_lib.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        n_components=2,
        random_state=random_state,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        coords = reducer.fit_transform(Xs)

    df_umap = pd.DataFrame(coords, index=X.index, columns=["UMAP1", "UMAP2"])
    return df_umap

# =============================================================================
# Group / stage helpers
# =============================================================================

def assign_groups(index):
    group_map = build_group_map(list(index))
    return pd.Series(group_map, index=index, name="group")

def assign_stages(index):
    return pd.Series({sid: get_stage(sid) for sid in index}, name="stage")

def make_label(sample_id, group):
    m  = re.search(r"bc(\d+)", sample_id)
    bc = f"bc{m.group(1)}" if m else sample_id.split("_")[-1]
    short = {
        "WT": "WT", "runt": "runt", "BOTR": "BOTR",
        "FoxL1_BOTv": "FoxL1", "B6_BOTC": "B6",
        "HLH54F_BOTCv": "HLH54F",
    }.get(group, group)
    return f"{short} {bc}"

# =============================================================================
# Generic scatter panel (shared by PCA and UMAP)
# =============================================================================

def _scatter_panel(
    df2d,          # DataFrame with two coordinate columns
    col_x, col_y, # column names for x / y
    groups,
    stages,
    title,
    xlabel, ylabel,
    draw_ellipses=True,
    label_embryos=False,
    use_stage_shapes=False,
    fixed_lims=None,   # (xlim, ylim) tuple to force axes — pass from with-ellipses run
    ax=None,
    standalone=True,
    qc_label_mode=False,  # v4: dedicated QC-label styling (larger font, dark outline)
):
    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 7))
    else:
        fig = ax.figure

    present = [g for g in GROUP_ORDER if g in groups.values]
    extras  = sorted(set(groups.values) - set(GROUP_ORDER))
    present = present + extras

    # ── Draw ellipses first (behind points) ──────────────────────────────────
    if draw_ellipses:
        for grp in present:
            mask   = (groups == grp)
            x, y   = df2d.loc[mask, col_x].values, df2d.loc[mask, col_y].values
            colour = GROUP_COLOURS.get(grp, "#888888")
            if len(x) >= 3:
                confidence_ellipse(
                    x, y, ax, n_std=2.0,
                    facecolor=colour, edgecolor=colour,
                    alpha=0.15, linewidth=1.2, linestyle="--", zorder=1,
                )

    # ── Scatter points ────────────────────────────────────────────────────────
    for grp in present:
        mask   = (groups == grp)
        colour = GROUP_COLOURS.get(grp, "#888888")
        pt_size = 100 if qc_label_mode else 75
        pt_alpha = 0.75 if qc_label_mode else 1.0
        for sid, xi, yi in zip(df2d.loc[mask].index,
                                df2d.loc[mask, col_x].values,
                                df2d.loc[mask, col_y].values):
            stage = stages.get(sid, "nc14b")
            marker = STAGE_MARKERS.get(stage, "o") if use_stage_shapes else "o"
            alpha  = pt_alpha if (qc_label_mode or not use_stage_shapes) else STAGE_ALPHA.get(stage, 1.0)
            ax.scatter(xi, yi, c=colour, marker=marker,
                       s=pt_size, zorder=3, edgecolors="white",
                       linewidths=0.6, alpha=alpha)
            if label_embryos:
                ax.annotate(make_label(sid, grp), (xi, yi),
                            fontsize=6.5, ha="left", va="bottom",
                            xytext=(3, 3), textcoords="offset points",
                            color=colour, zorder=4)
            if qc_label_mode:
                lbl = make_label(sid, grp)
                ax.annotate(
                    lbl, (xi, yi),
                    fontsize=8, ha="left", va="bottom",
                    xytext=(4, 4), textcoords="offset points",
                    color=colour, zorder=5, fontweight="bold",
                    bbox=dict(
                        boxstyle="round,pad=0.15",
                        facecolor="white", edgecolor=colour,
                        alpha=0.7, linewidth=0.6,
                    ),
                )

    # ── Reference lines ───────────────────────────────────────────────────────
    ax.axhline(0, color="lightgrey", linewidth=0.6, zorder=0)
    ax.axvline(0, color="lightgrey", linewidth=0.6, zorder=0)

    # ── Axis limits: apply fixed if provided, else auto then record ───────────
    # We must call draw to let matplotlib auto-scale before reading limits.
    fig.canvas.draw()
    auto_xlim = ax.get_xlim()
    auto_ylim = ax.get_ylim()

    if fixed_lims is not None:
        # Use the pre-computed "zoomed-out" limits (from the with-ellipses run)
        xlim, ylim = fixed_lims
        # Take the more zoomed-out of the two (union)
        xlim = (min(xlim[0], auto_xlim[0]), max(xlim[1], auto_xlim[1]))
        ylim = (min(ylim[0], auto_ylim[0]), max(ylim[1], auto_ylim[1]))
    else:
        xlim, ylim = auto_xlim, auto_ylim

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    # ── Labels / title ────────────────────────────────────────────────────────
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=10)

    # ── Group legend ──────────────────────────────────────────────────────────
    grp_handles = [
        mpatches.Patch(color=GROUP_COLOURS.get(g, "#888"), label=g)
        for g in present
    ]
    leg1 = ax.legend(handles=grp_handles, title="Group",
                     fontsize=9, title_fontsize=9,
                     loc="upper right", framealpha=0.85)
    ax.add_artist(leg1)

    # ── Stage legend (if multiple stages present) ─────────────────────────────
    stages_present = [s for s in ["nc14b", "nc14d", "gastr"] if s in stages.values]
    if len(stages_present) > 1:
        from matplotlib.lines import Line2D
        if use_stage_shapes:
            stage_handles = [
                Line2D([0], [0], marker=STAGE_MARKERS.get(s, "o"), linestyle="",
                       markerfacecolor="#444", markeredgecolor="white",
                       markersize=8, label=s)
                for s in stages_present
            ]
        else:
            stage_handles = [
                Line2D([0], [0], marker="o", linestyle="",
                       markerfacecolor="#444", markeredgecolor="white",
                       alpha=STAGE_ALPHA.get(s, 1.0), markersize=8, label=s)
                for s in stages_present
            ]
        ax.legend(handles=stage_handles, title="Stage",
                  fontsize=8, title_fontsize=8,
                  loc="lower right", framealpha=0.85)

    if standalone:
        fig.tight_layout()

    return fig, ax, (xlim, ylim)


# =============================================================================
# High-level save helper — produces both PDF and PNG
# =============================================================================

def save_fig(fig, outdir, stem):
    """Save as PDF (vector) and PNG (300 dpi raster)."""
    for ext in ("pdf", "png"):
        path = os.path.join(outdir, f"{stem}.{ext}")
        fig.savefig(path, bbox_inches="tight",
                    dpi=(300 if ext == "png" else None))
        print(f"[plot] Saved: {path}")
    plt.close(fig)


# =============================================================================
# Per-subset PCA outputs
# =============================================================================

def make_pca_subset(
    log_tpm, groups, stages,
    subset_name, subset_label,
    outdir, prefix, n_top, label_embryos=False, label_qc=False,
):
    """
    For one subset:
      0. Subset to this group's samples, then re-run PCA on JUST those
         samples — top-variable genes, scaling, and variance-explained
         are all recomputed from the subset, not inherited from the
         full-dataset PCA. (Previously this just masked rows out of a
         shared full-dataset embedding, so excluded groups like WT were
         still shaping which genes were selected and what the axes mean.)
      1. Run with_ellipses → capture axis limits
      2. Run no_ellipses   → force same limits
    Saves PDF + PNG for both.
    """
    mask       = groups.apply(SUBSETS[subset_name]["keep"])
    sample_ids = groups.index[mask]

    if len(sample_ids) < 2:
        print(f"[plot] Skipping PCA {subset_name}: not enough samples")
        return

    sub_log_tpm        = log_tpm[sample_ids]
    sub_coords, sub_pca = run_pca(sub_log_tpm, n_top=n_top)
    sub_groups          = groups.loc[sample_ids]
    sub_stages          = stages.loc[sample_ids]

    var_x = sub_pca.explained_variance_ratio_[0] * 100
    var_y = sub_pca.explained_variance_ratio_[1] * 100

    # ── with ellipses ─────────────────────────────────────────────────────────
    fig, ax, lims = _scatter_panel(
        sub_coords, "PC1", "PC2",
        sub_groups, sub_stages,
        title=f"PCA — {subset_label} (n={len(sub_coords)})",
        xlabel=f"PC1 ({var_x:.1f}% variance)",
        ylabel=f"PC2 ({var_y:.1f}% variance)",
        draw_ellipses=True,
        label_embryos=label_embryos,
        use_stage_shapes=True,
    )
    save_fig(fig, outdir, f"{prefix}_pca_{subset_name}_ellipses")

    # ── no ellipses (same limits) ─────────────────────────────────────────────
    fig2, ax2, _ = _scatter_panel(
        sub_coords, "PC1", "PC2",
        sub_groups, sub_stages,
        title=f"PCA — {subset_label} (n={len(sub_coords)})",
        xlabel=f"PC1 ({var_x:.1f}% variance)",
        ylabel=f"PC2 ({var_y:.1f}% variance)",
        draw_ellipses=False,
        label_embryos=label_embryos,
        fixed_lims=lims,
        use_stage_shapes=True,
    )
    save_fig(fig2, outdir, f"{prefix}_pca_{subset_name}_no_ellipses")

    # ── QC labeled (v4): no ellipses, every embryo labeled ───────────────────
    if label_qc:
        fig3, ax3, _ = _scatter_panel(
            sub_coords, "PC1", "PC2",
            sub_groups, sub_stages,
            title=f"PCA QC labels — {subset_label} (n={len(sub_coords)})",
            xlabel=f"PC1 ({var_x:.1f}% variance)",
            ylabel=f"PC2 ({var_y:.1f}% variance)",
            draw_ellipses=False,
            label_embryos=False,
            qc_label_mode=True,
            fixed_lims=lims,
            use_stage_shapes=True,
        )
        save_fig(fig3, outdir, f"{prefix}_pca_{subset_name}_qc_labeled")


# =============================================================================
# Per-subset UMAP outputs
# =============================================================================

def make_umap_subset(
    log_tpm, groups, stages,
    subset_name, subset_label,
    outdir, prefix, n_top, umap_neighbors, umap_min_dist, label_embryos=False, label_qc=False,
):
    mask       = groups.apply(SUBSETS[subset_name]["keep"])
    sample_ids = groups.index[mask]

    if len(sample_ids) < 2:
        print(f"[plot] Skipping UMAP {subset_name}: not enough samples")
        return

    sub_log_tpm = log_tpm[sample_ids]
    sub_umap = run_umap(
        sub_log_tpm, n_top=n_top,
        n_neighbors=umap_neighbors,
        min_dist=umap_min_dist,
    )
    if sub_umap is None:
        return

    sub_groups = groups.loc[sample_ids]
    sub_stages = stages.loc[sample_ids]

    # ── with ellipses ─────────────────────────────────────────────────────────
    fig, ax, lims = _scatter_panel(
        sub_umap, "UMAP1", "UMAP2",
        sub_groups, sub_stages,
        title=f"UMAP — {subset_label} (n={len(sub_umap)})",
        xlabel="UMAP1", ylabel="UMAP2",
        draw_ellipses=True,
        label_embryos=label_embryos,
        use_stage_shapes=True,
    )
    save_fig(fig, outdir, f"{prefix}_umap_{subset_name}_ellipses")

    # ── no ellipses (same limits) ─────────────────────────────────────────────
    fig2, ax2, _ = _scatter_panel(
        sub_umap, "UMAP1", "UMAP2",
        sub_groups, sub_stages,
        title=f"UMAP — {subset_label} (n={len(sub_umap)})",
        xlabel="UMAP1", ylabel="UMAP2",
        draw_ellipses=False,
        label_embryos=label_embryos,
        fixed_lims=lims,
        use_stage_shapes=True,
    )
    save_fig(fig2, outdir, f"{prefix}_umap_{subset_name}_no_ellipses")

    # ── QC labeled (v4): no ellipses, every embryo labeled ───────────────────
    if label_qc:
        fig3, ax3, _ = _scatter_panel(
            sub_umap, "UMAP1", "UMAP2",
            sub_groups, sub_stages,
            title=f"UMAP QC labels — {subset_label} (n={len(sub_umap)})",
            xlabel="UMAP1", ylabel="UMAP2",
            draw_ellipses=False,
            label_embryos=False,
            qc_label_mode=True,
            fixed_lims=lims,
            use_stage_shapes=True,
        )
        save_fig(fig3, outdir, f"{prefix}_umap_{subset_name}_qc_labeled")


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix",      default="all_genotypes_ruvg")
    parser.add_argument("--tpm-file",    default=None,
                        help="Direct path to a TPM or log2TPM matrix TSV "
                             "(bypasses prefix-based file search). "
                             "If the matrix max > 30 it is assumed to be "
                             "linear TPM and log2(x+1) is applied; "
                             "otherwise it is used as-is (already log2).")
    parser.add_argument("--results-dir", default="results/combined")
    parser.add_argument("--outdir",      default="results/combined_figures/pca_umap")
    parser.add_argument("--label-embryos", action="store_true")
    parser.add_argument("--label-qc", action="store_true",
                        help="Generate additional _qc_labeled outputs with every embryo "
                             "labeled (bold, outlined text) for outlier identification.")
    parser.add_argument("--n-top",       type=int, default=2000)
    parser.add_argument("--pcs",         nargs="+", type=int, default=[1, 2, 3])
    parser.add_argument("--umap-neighbors", type=int, default=15)
    parser.add_argument("--umap-min-dist",  type=float, default=0.3)
    parser.add_argument("--skip-umap",   action="store_true",
                        help="Skip UMAP (useful if umap-learn is not installed)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ── Load & reduce ─────────────────────────────────────────────────────────
    if args.tpm_file:
        print(f"[plot] Loading TPM from: {args.tpm_file}")
        log_tpm = pd.read_csv(args.tpm_file, sep="\t", index_col=0)
        if "GeneSymbol" in log_tpm.columns:
            log_tpm = log_tpm.drop(columns="GeneSymbol")
        # Detect whether file is already log2-transformed.
        # Post-RUVg log2TPM files (e.g. *_log2tpm_ruvg_corrected.tsv) are
        # already on the log2 scale (typical max ~15-20); raw TPM can be >1e5.
        if log_tpm.values.max() > 30:
            print("[plot] Matrix max > 30 — assuming linear TPM; applying log2(x+1).")
            log_tpm = np.log2(log_tpm + 1)
        else:
            print("[plot] Matrix max ≤ 30 — assuming already log2-transformed; using as-is.")
    else:
        log_tpm = load_tpm(args.prefix, args.results_dir)

    # Full-dataset PCA is still needed for the "all samples" multi-panel
    # figure at the end of main(); per-subset PCA/UMAP below are each
    # computed independently from log_tpm, not derived from this.
    df_coords, pca = run_pca(log_tpm, n_top=args.n_top)
    groups = assign_groups(log_tpm.columns)
    stages = assign_stages(log_tpm.columns)

    # ── Apply outlier exclusions (v5) ─────────────────────────────────────────
    keep_mask = pd.Series(
        [not _is_excluded(sid, groups[sid]) for sid in log_tpm.columns],
        index=log_tpm.columns,
    )
    n_excluded = (~keep_mask).sum()
    if n_excluded:
        excluded_ids = log_tpm.columns[~keep_mask].tolist()
        print(f"[plot] Excluding {n_excluded} outlier embryo(s): {excluded_ids}")
        log_tpm = log_tpm.loc[:, keep_mask]
        groups  = groups[keep_mask]
        stages  = stages[keep_mask]
        # Recompute full-dataset PCA on filtered set
        df_coords, pca = run_pca(log_tpm, n_top=args.n_top)

    # ── Per-subset PCA + UMAP (each subset re-runs PCA/UMAP on only its own
    #    samples, so excluded groups — e.g. WT in the *_no_wt subsets — don't
    #    influence gene selection, scaling, or variance explained) ──────────
    for subset_name, cfg in SUBSETS.items():
        print(f"\n[plot] === {subset_name} ===")
        make_pca_subset(
            log_tpm, groups, stages,
            subset_name, cfg["label"],
            args.outdir, args.prefix, args.n_top,
            label_embryos=args.label_embryos,
            label_qc=args.label_qc,
        )
        if not args.skip_umap:
            make_umap_subset(
                log_tpm, groups, stages,
                subset_name, cfg["label"],
                args.outdir, args.prefix, args.n_top,
                args.umap_neighbors, args.umap_min_dist,
                label_embryos=args.label_embryos,
                label_qc=args.label_qc,
            )

    # ── Multi-panel all-data PCA (retained from v2) ───────────────────────────
    pcs   = args.pcs
    pairs = [(pcs[i], pcs[j]) for i in range(len(pcs)) for j in range(i+1, len(pcs))]
    ncols = min(len(pairs), 2)
    nrows = (len(pairs) + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(9*ncols, 7*nrows))
    axes = np.array(axes).flatten() if len(pairs) > 1 else [axes]

    for idx, (pc_x, pc_y) in enumerate(pairs):
        var_x = pca.explained_variance_ratio_[pc_x-1] * 100
        var_y = pca.explained_variance_ratio_[pc_y-1] * 100
        _scatter_panel(
            df_coords, f"PC{pc_x}", f"PC{pc_y}",
            groups, stages,
            title=f"PCA — {args.prefix} PC{pc_x} vs PC{pc_y}",
            xlabel=f"PC{pc_x} ({var_x:.1f}%)", ylabel=f"PC{pc_y} ({var_y:.1f}%)",
            draw_ellipses=True,
            use_stage_shapes=True,
            ax=axes[idx], standalone=False,
        )
    for ax in axes[len(pairs):]:
        ax.set_visible(False)

    fig.suptitle(f"PCA — {args.prefix} (all samples)", fontsize=13, y=1.01)
    fig.tight_layout()
    save_fig(fig, args.outdir, f"{args.prefix}_pca_all_panels")

    print(f"\n[plot] All outputs written to: {args.outdir}")
    counts = groups.value_counts().to_dict()
    print(f"[plot] {len(df_coords)} embryos | groups: {counts}")


if __name__ == "__main__":
    main()
