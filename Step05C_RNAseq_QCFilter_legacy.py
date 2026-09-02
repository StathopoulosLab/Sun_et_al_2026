#!/usr/bin/env python3
"""
step5c_qc_filter.py  —  QC and filter embryo libraries before normalization.

Group-aware and stage-aware threshold system
--------------------------------------------
The most important change in this version: hard QC thresholds are applied
differently depending on biological group and developmental stage.

The motivation is biological: HLH54F_BOTCv and FoxL1_BOTv embryos — especially
at nc14b — have genuine transcriptomic differences from WT that make them appear
"lower quality" by naive metrics even when the library is technically fine:

  * FoxL1_BOTv nc14b: run is high (good), but many AP patterning genes that
    dominate WT libraries are absent/reduced due to the BOT maternal background.
    This lowers detected-gene counts relative to WT nc14b.
  * HLH54F_BOTCv nc14b: similar BOTCv background effect. tll/hkb are elevated
    but many mid-embryo genes are absent.
  * Gastrula-stage embryos (any genotype): naturally lower total UMIs because
    many maternally-loaded genes are degraded. A gastr embryo with 4,000 UMIs
    may be a perfectly good library.

Threshold tiers (applied in order, first match wins per embryo):
  1. PROTECTED groups: FoxL1_BOTv and HLH54F_BOTCv at nc14b are never failed
     on UMI/gene thresholds — only on mito% and the GOI zero-expression check.
     These are your rarest and most biologically critical libraries.
  2. EXPANDED thresholds (QC_EXPANDED_* from sample_config): applied to
     FoxL1_BOTv / HLH54F_BOTCv at nc14d and gastr, and to BOTR at nc14b.
  3. STANDARD thresholds (--min-umis / --min-genes): applied to everything else
     (WT, runt, B6_BOTC at any stage).

The GOI zero-expression check is applied to ALL groups including protected ones
— an embryo where every single diagnostic gene has 0 counts is a demux failure
regardless of genotype.

You can override the protection with --no-protect-critical if needed.

Usage
-----
  # Step 1: inspect
  python step5c_qc_filter.py --combined all_genotypes --report-only

  # Step 2: filter with group-aware thresholds
  python step5c_qc_filter.py --combined all_genotypes \\
      --min-umis 3000 --min-genes 2000 --mad 3 --goi-check

  # Step 3: manual exclusion on top of thresholds
  python step5c_qc_filter.py --combined all_genotypes \\
      --min-umis 3000 --min-genes 2000 \\
      --exclude runt_Dm_nc14b_bc16

  # Override protection (use only if you have a specific reason)
  python step5c_qc_filter.py --combined all_genotypes \\
      --min-umis 3000 --min-genes 2000 --no-protect-critical
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
import config

try:
    import Config_SampleMetadata as SC
    _HAS_SC = True
except ImportError:
    _HAS_SC = False

MITO_PATTERNS = ["mt:", "mito", "Mt:"]

# Key diagnostic genes for GOI expression check (symbols — resolved below)
GOI_CHECK_SYMBOLS = ["run", "eve", "tll", "sna", "twi", "hb"]

# ── Group-aware threshold tiers ───────────────────────────────────────────────
# PROTECTED: FoxL1_BOTv and HLH54F_BOTCv at nc14b.
# These carry the BOT/BOTCv maternal background which depletes many AP genes
# that dominate WT library diversity.  Failing them on UMI/gene counts would
# discard biologically critical and irreplaceable nc14b data.  Only mito%
# and the GOI all-zero check are applied.
PROTECTED_GROUPS = {"FoxL1_BOTv", "HLH54F_BOTCv"}
PROTECTED_STAGES = {"nc14b"}

# EXPANDED: lower thresholds for groups/stages where the BOT/BOTCv background
# or developmental stage naturally produces smaller libraries than WT nc14b.
EXPANDED_GROUPS_STAGES: set[tuple[str, str]] = {
    ("FoxL1_BOTv",    "nc14d"),
    ("FoxL1_BOTv",    "gastr"),
    ("HLH54F_BOTCv",  "nc14d"),
    ("HLH54F_BOTCv",  "gastr"),
    ("BOTR",          "nc14b"),
    ("BOTR",          "nc14d"),
    ("B6_BOTC",       "nc14d"),
    ("B6_BOTC",       "gastr"),
}


def _get_threshold_tier(group: str, stage: str,
                        protect_critical: bool = True) -> str:
    """Return 'protected', 'expanded', or 'standard' for an embryo."""
    if protect_critical:
        if group in PROTECTED_GROUPS and stage in PROTECTED_STAGES:
            return "protected"
    if (group, stage) in EXPANDED_GROUPS_STAGES:
        return "expanded"
    return "standard"


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--combined",      required=True)
    p.add_argument("--report-only",   action="store_true")

    g = p.add_argument_group("Thresholds (0 = disabled)")
    g.add_argument("--min-umis",           type=int,   default=0)
    g.add_argument("--min-genes",          type=int,   default=0)
    g.add_argument("--max-mito",           type=float, default=0)
    g.add_argument("--min-complexity",     type=float, default=0)
    g.add_argument("--min-median-log2tpm", type=float, default=0,
                   help="Minimum per-embryo median log2(TPM+1) over detected genes. "
                        "Use sample_config.QC_MIN_MEDIAN_LOG2TPM as a starting point.")
    g.add_argument("--mad",               type=float, default=0,
                   help="Flag embryos >N MADs below group/global median UMIs+genes "
                        "(try 3). Applied per-group when group n≥4, else globally.")
    g.add_argument("--goi-check",         action="store_true",
                   help="Flag embryos where all key diagnostic genes have 0 counts "
                        "(likely a demux/alignment failure)")

    p.add_argument("--exclude", nargs="*", default=[],
                   help="Bare embryo IDs to manually exclude (without sample prefix)")
    p.add_argument("--no-protect-critical", action="store_true",
                   help="Disable protection of FoxL1_BOTv and HLH54F_BOTCv nc14b "
                        "embryos from UMI/gene thresholds. Use only if you have a "
                        "specific reason — these libraries are genuinely lower-complexity "
                        "due to the BOT maternal background, not technical failure.")
    return p.parse_args()


def _load_sample_groups(combined_dir: str, prefix: str) -> pd.DataFrame | None:
    """Load enriched sample-groups table from step5b.

    The groups file index uses full prefixed column names
    (e.g. 'WT_seRNAseq_4embryos_set1__WT_embryo_bc2') but the count
    matrix columns are bare names ('WT_embryo_bc2').  This function
    normalises the index to bare names so the join in compute_qc_metrics
    succeeds.
    """
    # Prefer the full, authoritative sample_groups.tsv over any existing
    # _sample_groups_qc.tsv. The _qc.tsv is a byproduct written fresh by
    # every successful filtering run (see bottom of main()), containing
    # ONLY that run's passing embryos. If an earlier run was more
    # restrictive (e.g. a bug, or stricter thresholds) than the current
    # one, loading its stale _qc.tsv here would silently drop group/stage
    # metadata for every embryo that wasn't in that earlier survivor list
    # — even ones that legitimately pass now. _qc.tsv is kept only as a
    # last-resort fallback for datasets where the full file was never
    # generated.
    for fname in (f"{prefix}_sample_groups.tsv",
                  f"{prefix}_sample_groups_qc.tsv"):
        path = os.path.join(combined_dir, fname)
        if os.path.exists(path):
            df = pd.read_csv(path, sep="\t")
            if "Embryo" not in df.columns:
                return None
            # Normalise index: strip sample prefix (everything before __)
            df["Embryo"] = df["Embryo"].apply(
                lambda e: e.split("__")[-1] if "__" in str(e) else str(e)
            )
            # Drop duplicates that arise when two sample prefixes map to
            # the same bare name (shouldn't happen after step5b fixes,
            # but defensive)
            df = df.drop_duplicates(subset="Embryo", keep="first")
            return df.set_index("Embryo")
    return None


def _identify_mito_genes(counts: pd.DataFrame,
                          gene_symbols: "pd.Series | None") -> pd.Index:
    mask = pd.Series(False, index=counts.index)
    if gene_symbols is not None:
        for pat in MITO_PATTERNS:
            mask |= gene_symbols.str.contains(pat, na=False)
    for pat in MITO_PATTERNS:
        mask |= counts.index.str.contains(pat, case=False)
    return counts.index[mask]


def _resolve_goi(counts_index: pd.Index,
                 gene_symbols: "pd.Series | None") -> dict[str, str]:
    """Return {symbol: row_id} for GOI_CHECK_SYMBOLS found in the matrix."""
    resolved: dict[str, str] = {}
    sym_to_id: dict[str, str] = {}
    if gene_symbols is not None:
        sym_to_id = dict(zip(gene_symbols.values, gene_symbols.index))

    for sym in GOI_CHECK_SYMBOLS:
        if sym in counts_index:
            resolved[sym] = sym
        elif sym in sym_to_id and sym_to_id[sym] in counts_index:
            resolved[sym] = sym_to_id[sym]
    return resolved


def compute_qc_metrics(counts: pd.DataFrame,
                        gene_symbols: "pd.Series | None",
                        groups_df: "pd.DataFrame | None") -> pd.DataFrame:
    embryo_cols = counts.columns.tolist()
    m = pd.DataFrame(index=embryo_cols)
    m.index.name = "Embryo"

    # Derive group and stage directly from sample_config (authoritative source).
    # Never read these from the groups TSV — it may be an old version without
    # those columns, which is why the tier lookup was silently returning standard.
    try:
        import Config_SampleMetadata as _SC
        _gmap = _SC.build_group_map(embryo_cols)
        m["group"] = [_gmap.get(c, "unknown") for c in embryo_cols]
        m["stage"] = [_SC.get_stage(c) for c in embryo_cols]
    except Exception:
        pass
    if groups_df is not None:
        for col in ("Sample", "set"):
            if col in groups_df.columns and col not in m.columns:
                m[col] = groups_df[col].reindex(m.index)

    m["Total_UMIs"]     = counts.sum(axis=0)
    m["Detected_Genes"] = (counts > 0).sum(axis=0)
    m["Complexity"]     = m["Detected_Genes"] / m["Total_UMIs"].replace(0, np.nan)

    # Median log2(CPM+1): scale raw counts to counts-per-million using each
    # embryo's own total UMIs, THEN log2. Previously this took log2 of the
    # raw count directly (no library-size scaling at all), which meant the
    # metric only ever reflected "does the median detected gene have 1 or 2
    # raw UMIs" — true of nearly every single-embryo CEL-Seq2 library
    # regardless of quality, since most detected genes sit at very low raw
    # counts even in deep libraries. That made QC_MIN_MEDIAN_LOG2TPM (default
    # 4.0) fail almost everything once actually invoked, independent of any
    # real quality difference between embryos.
    cpm = counts.div(m["Total_UMIs"].replace(0, np.nan), axis=1) * 1e6
    log2cpm   = np.log2(cpm + 1)
    expressed = counts > 0
    m["Median_log2CPM"] = pd.Series(
        {col: float(log2cpm.loc[expressed[col], col].median())
         for col in embryo_cols},
        dtype=float,
    )

    mito_genes = _identify_mito_genes(counts, gene_symbols)
    if len(mito_genes):
        m["Mito_Pct"] = 100 * counts.loc[mito_genes].sum(axis=0) / m["Total_UMIs"]
        print(f"[step5c] Mitochondrial genes: {len(mito_genes)}")
    else:
        m["Mito_Pct"] = 0.0
        print("[step5c] WARNING: no mitochondrial genes detected — mito filter disabled")

    return m


def flag_embryos(metrics: pd.DataFrame, counts: pd.DataFrame,
                  gene_symbols: "pd.Series | None", args) -> pd.Series:
    passes  = pd.Series(True,  index=metrics.index)
    reasons = pd.Series("",    index=metrics.index)
    tiers   = pd.Series("standard", index=metrics.index)  # for reporting

    protect_critical = not getattr(args, "no_protect_critical", False)

    # ── Resolve per-embryo threshold tier ─────────────────────────────────────
    has_group = "group" in metrics.columns
    has_stage = "stage" in metrics.columns

    # Pull expanded thresholds from sample_config when available
    exp_min_umis  = SC.QC_EXPANDED_MIN_UMIS  if _HAS_SC else 0
    exp_min_genes = SC.QC_EXPANDED_MIN_GENES if _HAS_SC else 0

    for embryo in metrics.index:
        grp   = str(metrics.loc[embryo, "group"]) if has_group else "unknown"
        stage = str(metrics.loc[embryo, "stage"]) if has_stage else "unknown"
        tiers[embryo] = _get_threshold_tier(grp, stage, protect_critical)

    # Report tier breakdown
    tier_counts = tiers.value_counts().to_dict()
    print(f"[step5c] Threshold tiers: "
          + "  ".join(f"{t}: {n}" for t, n in sorted(tier_counts.items())))
    if "protected" in tier_counts:
        prot_embryos = metrics.index[tiers == "protected"].tolist()
        print(f"[step5c] Protected (UMI/gene filters skipped):")
        for e in prot_embryos:
            bare = e.split("__")[-1]
            grp  = metrics.loc[e, "group"] if has_group else ""
            stg  = metrics.loc[e, "stage"] if has_stage else ""
            print(f"           {bare}  ({grp} / {stg})")

    # ── Manual exclusions ─────────────────────────────────────────────────────
    for excl in args.exclude:
        matched = [e for e in passes.index
                   if e.split("__")[-1] == excl or e == excl]
        if matched:
            for m in matched:
                passes[m]  = False
                reasons[m] += "manual_exclude; "
        else:
            print(f"[step5c] WARNING: --exclude '{excl}' not matched. "
                  f"Available bare IDs (first 5): "
                  f"{[e.split('__')[-1] for e in metrics.index[:5]]} ...")

    # ── Hard thresholds — tiered ──────────────────────────────────────────────
    # For each embryo we apply a different threshold depending on its tier:
    #   protected → skip UMI + gene + complexity + median-log2 thresholds
    #   expanded  → use QC_EXPANDED_* values from sample_config
    #   standard  → use the --min-umis / --min-genes arguments

    std_min_umis       = args.min_umis
    std_min_genes      = args.min_genes
    std_min_complexity = args.min_complexity
    std_min_medlog2    = args.min_median_log2tpm

    for embryo in metrics.index:
        tier = tiers[embryo]
        if tier == "protected":
            # Only mito and GOI-zero apply — skip everything else here
            continue

        row = metrics.loc[embryo]

        if tier == "expanded":
            min_u = exp_min_umis  if exp_min_umis  > 0 else std_min_umis
            min_g = exp_min_genes if exp_min_genes > 0 else std_min_genes
            # Expanded: no complexity or median-log2 threshold (too noisy at small n)
            min_c = 0.0
            min_m = 0.0
        else:  # standard
            min_u = std_min_umis
            min_g = std_min_genes
            min_c = std_min_complexity
            min_m = std_min_medlog2

        if min_u > 0 and row["Total_UMIs"] < min_u:
            passes[embryo]  = False
            reasons[embryo] += f"UMIs<{min_u}({tier}); "

        if min_g > 0 and row["Detected_Genes"] < min_g:
            passes[embryo]  = False
            reasons[embryo] += f"genes<{min_g}({tier}); "

        if min_c > 0 and row["Complexity"] < min_c:
            passes[embryo]  = False
            reasons[embryo] += f"cmplx<{min_c}; "

        if min_m > 0 and row["Median_log2CPM"] < min_m:
            passes[embryo]  = False
            reasons[embryo] += f"medlog2<{min_m}; "

    # ── Mito threshold — applied to ALL tiers ─────────────────────────────────
    if args.max_mito > 0:
        fail = metrics["Mito_Pct"] > args.max_mito
        passes[fail]  = False
        reasons[fail] += f"mito>{args.max_mito}%; "

    # ── MAD-based outlier detection — NEVER applied to protected embryos ──────
    if args.mad > 0:
        non_protected = metrics.index[tiers != "protected"]
        group_col = metrics["group"] if has_group else None

        for metric_col in ("Total_UMIs", "Detected_Genes"):
            vals_log = np.log10(metrics[metric_col].clip(lower=1))

            # Per-group MAD (only within non-protected embryos of each group)
            if group_col is not None:
                for grp in group_col.unique():
                    gmask = (group_col == grp) & (tiers != "protected")
                    gvals = vals_log[gmask]
                    if gmask.sum() < 4:
                        continue
                    med = gvals.median()
                    mad = np.median(np.abs(gvals - med))
                    if mad > 0:
                        lower = med - args.mad * 1.4826 * mad
                        fail  = gmask & (vals_log < lower)
                        passes[fail]  = False
                        reasons[fail] += f"{metric_col}_MAD_{grp}; "

            # Global MAD across non-protected only
            np_vals = vals_log[non_protected]
            med = np_vals.median()
            mad = np.median(np.abs(np_vals - med))
            if mad > 0:
                lower = med - args.mad * 1.4826 * mad
                fail  = (vals_log < lower) & (tiers != "protected")
                new_fail = fail & passes
                passes[new_fail]  = False
                reasons[new_fail] += f"{metric_col}_globalMAD; "

    # ── GOI zero-expression check — ALL tiers including protected ─────────────
    if args.goi_check:
        goi_ids = _resolve_goi(counts.index, gene_symbols)
        if goi_ids:
            goi_rows = list(goi_ids.values())
            goi_sums = counts.loc[goi_rows].sum(axis=0)
            all_zero = goi_sums == 0
            passes[all_zero]  = False
            reasons[all_zero] += "all_GOI_zero(demux_failure?); "
            n_zero = all_zero.sum()
            if n_zero:
                print(f"[step5c] GOI check: {n_zero} embryo(s) with 0 counts "
                      f"across all of {list(goi_ids.keys())} — "
                      f"likely demux/alignment failure")
            else:
                print(f"[step5c] GOI check: all embryos express ≥1 diagnostic gene ✓")
        else:
            print("[step5c] GOI check: no diagnostic genes found in matrix — skipping")

    metrics["QC_Pass"]   = passes
    metrics["QC_Reason"] = reasons.str.rstrip("; ")
    metrics["QC_Tier"]   = tiers   # expose tier in output TSV
    return passes

    metrics["QC_Pass"]   = passes
    metrics["QC_Reason"] = reasons.str.rstrip("; ")
    return passes


def plot_qc(metrics: pd.DataFrame, out_path: str, thresholds: dict):
    """Four-panel QC plot coloured by biological group when available."""
    try:
        if _HAS_SC:
            grp_colours = SC.GROUP_COLOURS
        else:
            grp_colours = {}
    except Exception:
        grp_colours = {}

    has_group = "group" in metrics.columns
    if has_group:
        default_col = "#AAAAAA"
        point_colors = [grp_colours.get(g, default_col)
                        for g in metrics["group"]]
        bar_colors   = [grp_colours.get(g, default_col)
                        for g in metrics["group"]]
    else:
        pass_colors  = ["#2171B5" if p else "#CB181D" for p in metrics["QC_Pass"]]
        point_colors = pass_colors
        bar_colors   = pass_colors

    # Overlay fail with red edge
    edge_colors = ["#CB181D" if not p else "none"
                   for p in metrics["QC_Pass"]]
    edge_widths = [2.0 if not p else 0.0 for p in metrics["QC_Pass"]]

    xlabels = [e.split("__")[-1] if "__" in e else e for e in metrics.index]

    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    fig.suptitle("Per-embryo QC metrics\n(red border = fail)", fontsize=13)

    def _bar(ax, col, ylabel, title, thresh_val=None, thresh_label=None):
        vals = metrics[col].values
        bars = ax.bar(range(len(metrics)), vals,
                      color=bar_colors, edgecolor=edge_colors,
                      linewidth=edge_widths)
        if thresh_val:
            ax.axhline(thresh_val, color="red", ls="--", lw=1.2,
                       label=thresh_label or str(thresh_val))
            ax.legend(fontsize=8)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_xticks(range(len(metrics)))
        ax.set_xticklabels(xlabels, rotation=90, fontsize=6)

    _bar(axes[0, 0], "Total_UMIs",     "Total UMIs",   "Total UMI counts",
         thresholds.get("min_umis"),  f"min: {thresholds.get('min_umis'):,}" if thresholds.get("min_umis") else None)
    _bar(axes[0, 1], "Detected_Genes", "Detected genes", "Genes detected (count>0)",
         thresholds.get("min_genes"), f"min: {thresholds.get('min_genes'):,}" if thresholds.get("min_genes") else None)
    _bar(axes[1, 0], "Mito_Pct",       "Mitochondrial %", "Mitochondrial fraction",
         thresholds.get("max_mito"),  f"max: {thresholds.get('max_mito')}%" if thresholds.get("max_mito") else None)

    # UMIs vs Genes scatter
    ax = axes[1, 1]
    for i, (embryo, row) in enumerate(metrics.iterrows()):
        ax.scatter(row["Total_UMIs"], row["Detected_Genes"],
                   c=point_colors[i],
                   edgecolors=edge_colors[i] if edge_colors[i] != "none" else "#333",
                   linewidths=edge_widths[i] if edge_widths[i] > 0 else 0.5,
                   marker="o" if row["QC_Pass"] else "X", s=50)
    max_umi = metrics["Total_UMIs"].max()
    for ratio in (0.05, 0.10, 0.20):
        x = np.array([0, max_umi])
        ax.plot(x, x * ratio, "--", color="gray", alpha=0.35, lw=0.8)
        ax.annotate(f"{ratio:.0%}", (max_umi * 0.93, max_umi * ratio),
                    fontsize=7, color="gray")
    ax.set_xlabel("Total UMIs")
    ax.set_ylabel("Detected genes")
    ax.set_title("Library complexity (UMIs vs genes)")

    # Group colour legend
    if has_group and grp_colours:
        import matplotlib.patches as mp
        handles = [mp.Patch(color=grp_colours.get(g, "#AAA"), label=g)
                   for g in metrics["group"].unique()
                   if not pd.isna(g)]
        axes[1, 1].legend(handles=handles, fontsize=7, loc="lower right",
                          title="Group", title_fontsize=7)

    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    args      = parse_args()
    cdir      = os.path.join(config.OUTPUT_ROOT, "combined")
    counts_path = os.path.join(cdir, f"{args.combined}_counts_matrix.tsv")

    if not os.path.exists(counts_path):
        sys.exit(f"[step5c] Not found: {counts_path}\n"
                  "         Run step5b first.")

    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    gene_symbols = counts.pop("GeneSymbol") if "GeneSymbol" in counts.columns else None

    # Normalise column names to bare embryo IDs (strip sample__ prefix if present).
    # step5b writes full prefixed names (e.g. 'WT_seRNAseq_4embryos_set1__WT_embryo_bc2')
    # into the count matrix; step5c, sample_groups.tsv, and all display code
    # work in bare-name space ('WT_embryo_bc2').  Strip here once so the rest
    # of the script never has to think about it.
    if any("__" in c for c in counts.columns):
        bare_map = {c: c.split("__")[-1] for c in counts.columns}
        counts = counts.rename(columns=bare_map)
        print(f"[step5c] Stripped sample prefixes from {sum(1 for c in bare_map if '__' in c)} columns")

    groups_df = _load_sample_groups(cdir, args.combined)
    if groups_df is not None:
        print(f"[step5c] Loaded group metadata ({len(groups_df)} embryos)")
    else:
        print("[step5c] WARNING: sample_groups.tsv not found — no group metadata")

    print(f"[step5c] Matrix: {counts.shape[0]:,} genes × {counts.shape[1]} embryos")

    metrics = compute_qc_metrics(counts, gene_symbols, groups_df)
    passes  = flag_embryos(metrics, counts, gene_symbols, args)

    n_pass = passes.sum()
    n_fail = (~passes).sum()

    # ── Save metrics ──────────────────────────────────────────────────────────
    metrics_path = os.path.join(cdir, f"{args.combined}_qc_metrics.tsv")
    metrics.to_csv(metrics_path, sep="\t", float_format="%.4f")

    # ── Print summary table ───────────────────────────────────────────────────
    print(f"\n[step5c] QC Summary: {len(metrics)} total  |  "
          f"{n_pass} PASS  |  {n_fail} FAIL")

    has_grp = "group" in metrics.columns
    hdr = (f"{'Embryo':<48} {'Group':<16} {'Stage':<8} {'Tier':<10} "
           f"{'UMIs':>10} {'Genes':>7} {'Mito%':>6} {'MedLog':>7} {'Pass':>5}  Reason"
           if has_grp else
           f"{'Embryo':<48} {'Tier':<10} {'UMIs':>10} {'Genes':>7} "
           f"{'Mito%':>6} {'MedLog':>7} {'Pass':>5}  Reason")
    print(f"\n{hdr}")
    print("─" * (len(hdr) + 20))

    for embryo, row in metrics.iterrows():
        bare   = embryo.split("__")[-1]
        status = "✓" if row["QC_Pass"] else "✗"
        reason = row.get("QC_Reason", "")
        tier   = row.get("QC_Tier", "standard")
        # Highlight protected embryos clearly
        tier_display = f"[{tier}]" if tier != "standard" else tier
        if has_grp:
            print(f"{bare:<48} {str(row.get('group','')):<16} "
                  f"{str(row.get('stage','')):<8} {tier_display:<10} "
                  f"{int(row['Total_UMIs']):>10,} {int(row['Detected_Genes']):>7,} "
                  f"{row['Mito_Pct']:>5.1f}% {row['Median_log2CPM']:>7.2f} "
                  f"{status:>5}  {reason}")
        else:
            print(f"{bare:<48} {tier_display:<10} {int(row['Total_UMIs']):>10,} "
                  f"{int(row['Detected_Genes']):>7,} "
                  f"{row['Mito_Pct']:>5.1f}% {row['Median_log2CPM']:>7.2f} "
                  f"{status:>5}  {reason}")

    if has_grp:
        print(f"\n[step5c] Per-group breakdown:")
        for grp in metrics["group"].dropna().unique():
            gmask = metrics["group"] == grp
            kept  = (passes & gmask).sum()
            total = gmask.sum()
            med_u = int(metrics.loc[gmask, "Total_UMIs"].median())
            med_g = int(metrics.loc[gmask, "Detected_Genes"].median())
            print(f"  {grp:<18}: {kept:>2}/{total:<2} pass  "
                  f"(med UMIs {med_u:,}  med genes {med_g:,})")

    # ── Plot ──────────────────────────────────────────────────────────────────
    thresholds = {
        "min_umis":  args.min_umis  or None,
        "min_genes": args.min_genes or None,
        "max_mito":  args.max_mito  or None,
    }
    plot_path = os.path.join(cdir, f"{args.combined}_qc_summary.png")
    plot_qc(metrics, plot_path, thresholds)
    print(f"\n[step5c] QC plot → {plot_path}")
    print(f"[step5c] Metrics  → {metrics_path}")

    if args.report_only:
        print(f"\n[step5c] Report-only — no filtering applied.")
        print(f"\n[step5c] Threshold tier that would be applied per embryo:")
        protect = not getattr(args, "no_protect_critical", False)
        has_g = "group" in metrics.columns
        has_s = "stage" in metrics.columns
        tier_summary: dict[str, list[str]] = {
            "protected": [], "expanded": [], "standard": []}
        for e, row in metrics.iterrows():
            grp   = str(row["group"]) if has_g else "unknown"
            stage = str(row["stage"]) if has_s else "unknown"
            tier  = _get_threshold_tier(grp, stage, protect)
            tier_summary[tier].append(e.split("__")[-1])
        for tier, embryos in tier_summary.items():
            if not embryos:
                continue
            print(f"\n  {tier.upper()} ({len(embryos)} embryo(s)):")
            if tier == "protected":
                print(f"    → UMI/gene/complexity thresholds SKIPPED")
                print(f"    → Only mito% and GOI-zero check applied")
            elif tier == "expanded":
                exp_u = SC.QC_EXPANDED_MIN_UMIS  if _HAS_SC else "n/a"
                exp_g = SC.QC_EXPANDED_MIN_GENES if _HAS_SC else "n/a"
                print(f"    → min_umis={exp_u:,}  min_genes={exp_g:,}  "
                      f"(no complexity/median-log2 filter)")
            else:
                print(f"    → min_umis={args.min_umis or 'not set'}  "
                      f"min_genes={args.min_genes or 'not set'}  "
                      f"min_median_log2tpm={args.min_median_log2tpm or 'not set'}")
            for e in embryos:
                print(f"    {e}")

        print(f"\n[step5c] Sample_config defaults:")
        if _HAS_SC:
            print(f"  QC_MIN_TOTAL_UMIS     = {SC.QC_MIN_TOTAL_UMIS:,}")
            print(f"  QC_MIN_DETECTED_GENES = {SC.QC_MIN_DETECTED_GENES:,}")
            print(f"  QC_MIN_MEDIAN_LOG2TPM = {SC.QC_MIN_MEDIAN_LOG2TPM}")
            print(f"  QC_EXPANDED_MIN_UMIS  = {SC.QC_EXPANDED_MIN_UMIS:,}  ← used for expanded tier")
            print(f"  QC_EXPANDED_MIN_GENES = {SC.QC_EXPANDED_MIN_GENES:,}  ← used for expanded tier")
        print(f"\n[step5c] Suggested run command:")
        print(f"  python step5c_qc_filter.py --combined {args.combined} \\")
        if _HAS_SC:
            print(f"      --min-umis {SC.QC_MIN_TOTAL_UMIS} "
                  f"--min-genes {SC.QC_MIN_DETECTED_GENES} \\")
            print(f"      --min-median-log2tpm {SC.QC_MIN_MEDIAN_LOG2TPM} \\")
        print(f"      --mad 3 --goi-check")
        return

    # ── Filter ────────────────────────────────────────────────────────────────
    keep_cols = passes[passes].index.tolist()
    filtered  = counts[keep_cols]
    if gene_symbols is not None:
        filtered.insert(0, "GeneSymbol", gene_symbols)

    out_path = os.path.join(cdir, f"{args.combined}_counts_matrix_qc.tsv")
    filtered.to_csv(out_path, sep="\t")

    # Exclusion list (TSV with reason)
    excl_rows = [(e, metrics.loc[e, "QC_Reason"]) for e in passes[~passes].index]
    excl_path = os.path.join(cdir, f"{args.combined}_excluded_embryos.tsv")
    pd.DataFrame(excl_rows, columns=["Embryo", "Reason"]).to_csv(
        excl_path, sep="\t", index=False)

    # Update sample-groups for downstream scripts.
    # group/stage come from `metrics` (freshly computed via
    # sample_config.build_group_map for every currently-passing embryo,
    # already verified correct above) rather than from groups_df — this
    # keeps the written file correct even if groups_df is stale, missing
    # entries, or predates a Config_SampleMetadata.py fix. Sample/set are
    # best-effort passthrough from groups_df since metrics doesn't carry them.
    has_group_stage = "group" in metrics.columns
    if has_group_stage or groups_df is not None:
        cols = [c for c in ("group", "stage") if c in metrics.columns]
        out_rows = metrics.loc[keep_cols, cols].copy() if cols else pd.DataFrame(index=keep_cols)
        if groups_df is not None:
            for col in ("Sample", "set"):
                if col in groups_df.columns:
                    out_rows[col] = groups_df[col].reindex(out_rows.index)
        out_rows.index.name = "Embryo"
        out_rows.reset_index().to_csv(
            os.path.join(cdir, f"{args.combined}_sample_groups_qc.tsv"),
            sep="\t", index=False)

    print(f"\n[step5c] Filtered: {len(keep_cols)} embryos retained, "
          f"{n_fail} removed")
    if excl_rows:
        for e, r in excl_rows:
            print(f"  ✗ {e.split('__')[-1]}  ({r})")
    print(f"[step5c] Filtered matrix → {out_path}")
    print(f"\n[step5c] Next: python step6_tpm.py --combined {args.combined} "
          f"--quantile --qc")


if __name__ == "__main__":
    main()
