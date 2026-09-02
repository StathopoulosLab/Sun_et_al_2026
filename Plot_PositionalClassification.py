#!/usr/bin/env python3
"""
plot_positional_classification.py  —  Stage-aware positional / DV classification
of FoxL1_BOTv and HLH54F_BOTCv embryos vs WT.

Major changes vs previous version
----------------------------------
Pattern gene expansion
  HLH54F_BOTCv (posterior mesoderm / tail identity):
    HIGH: hb, gt, tll, hkb, fkh, byn, Abd-B  (hkb added — behaves like tll)
    LOW : Kr, kni, Antp, Ubx, abd-A

  FoxL1_BOTv (anterior/mid positional shift):
    HIGH: gt, Antp, Abd-B, cad, prd, ftz, Dichaete, opa, run (higher at nc14d+), FoxL1
    LOW : hb, tll, Ubx, abd-A

DV axis — differentiate fully ventralized → unventralized:
  Ventral markers : dl, sna, twi, NetA
  Lateral markers : sog, brk, pyr, elav
  Dorsal markers  : Doc1, anc, Abd-B (also in positional)
  Note: Abd-B present dorsally and posteriorly — scored in both contexts

Stage-aware z-score reference
  For FoxL1_BOTv and HLH54F_BOTCv, z-scores are computed within stage bins
  (nc14b / nc14d / gastr) when ≥ 2 embryos per bin, then merged.  This means
  an nc14b embryo is not penalised for expressing stage-progressive genes at
  lower levels than gastrula embryos.

Stage-adjusted weights
  Genes whose signal strengthens or weakens with time carry stage-specific
  weight multipliers (see STAGE_PROGRESSIVE_GENES).

Outputs
-------
  pattern_scores_HLH54F.pdf/png
  pattern_scores_foxl1.pdf/png
  heatmap_HLH54F_vs_WT.pdf/png
  heatmap_foxl1_vs_WT.pdf/png
  heatmap_all_groups.pdf/png
  heatmap_HLH54F_by_call.pdf/png
  heatmap_foxl1_by_call.pdf/png
  heatmap_match_partial_combined.pdf/png
  violin_HLH54F_vs_WT.pdf/png          ← now includes pair-rule context genes
  violin_foxl1_vs_WT.pdf/png
  violin_HLH54F_by_call.pdf/png
  violin_foxl1_by_call.pdf/png
  violin_all_groups.pdf/png
  violin_match_partial_combined.pdf/png
  violin_DV_axis.pdf/png               ← NEW: ventral/lateral/dorsal per group
  matched_embryos_HLH54F_BOTCv.tsv
  matched_embryos_FoxL1_BOTv.tsv
  positional_classification.tsv

Usage
-----
  python plot_positional_classification_v2.py
  python plot_positional_classification_v2.py --match-threshold 0.65
  python plot_positional_classification_v2.py --no-all-groups

Edit EXCLUDE_EMBRYOS near the top of the script to adjust which embryos
are dropped before scoring, z-matrix construction, and plotting.
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
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

try:
    import Config_SampleMetadata as SC
except ImportError:
    sys.exit("[posclass] Config_SampleMetadata.py not found.")
try:
    import Config_DataLoader as DL
except ImportError:
    sys.exit("[posclass] Config_DataLoader.py not found.")

# =============================================================================
#  CONFIGURATION
# =============================================================================

OUT_DIR    = os.path.join("results", "combined_figures", "positional_classification")
SYMBOL_MAP = "gene_id_to_symbol.tsv"

# ── Stage-aware groups ────────────────────────────────────────────────────────
STAGE_AWARE_GROUPS = {"FoxL1_BOTv", "HLH54F_BOTCv"}

# Genes whose signal strengthens (+) or weakens (-) with stage.
# Values are weight multipliers: base_weight × multiplier.
STAGE_PROGRESSIVE_GENES: dict[str, dict[str, dict[str, float]]] = {
    "HLH54F_BOTCv": {
        "tll":   {"nc14b": 0.6, "nc14d": 1.0, "gastr": 1.4},
        "hkb":   {"nc14b": 0.6, "nc14d": 1.0, "gastr": 1.4},
        "Abd-B": {"nc14b": 0.5, "nc14d": 1.0, "gastr": 1.3},
        "byn":   {"nc14b": 0.5, "nc14d": 1.0, "gastr": 1.4},
        "fkh":   {"nc14b": 0.7, "nc14d": 1.0, "gastr": 1.2},
        "HLH54F":{"nc14b": 0.8, "nc14d": 1.0, "gastr": 1.0},
    },
    "FoxL1_BOTv": {
        "tll":   {"nc14b": 1.4, "nc14d": 1.0, "gastr": 0.7},  # LOW in foxl1 — stronger signal early
        "run":   {"nc14b": 0.8, "nc14d": 1.0, "gastr": 1.2},  # run HIGH in later stages
        "FoxL1": {"nc14b": 0.6, "nc14d": 1.0, "gastr": 1.3},
        "Kr":    {"nc14b": 0.8, "nc14d": 1.0, "gastr": 1.1},
        "Abd-B": {"nc14b": 0.5, "nc14d": 1.0, "gastr": 1.2},
    },
}

def _stage_weight(gene: str, genotype_group: str, stage: str,
                  base_weight: float) -> float:
    mod = STAGE_PROGRESSIVE_GENES.get(genotype_group, {}).get(gene, {})
    return base_weight * mod.get(stage, 1.0)

# ── HLH54F_BOTCv pattern ──────────────────────────────────────────────────────
# direction: +1 = expected HIGH vs WT, -1 = expected LOW
HLH54F_PATTERN: list[tuple[str, int]] = [
    ("hb",    +1),  # hunchback — HIGH
    ("gt",    +1),  # giant — HIGH
    ("tll",   +1),  # tailless — HIGH (strong specific marker)
    ("hkb",   +1),  # huckebein — HIGH (similar to tll in posterior)
    ("fkh",   +1),  # forkhead — HIGH (posterior gut)
    ("byn",   +1),  # brachyenteron — HIGH (posterior specific, BOTCv background)
    ("Abd-B", +1),  # Abdominal-B — HIGH
    ("HLH54F",+1),  # genotype gene
    ("Kr",    -1),  # Krüppel — LOW
    ("kni",   -1),  # knirps — LOW
    ("Ubx",   -1),  # Ultrabithorax — LOW
    ("abd-A", -1),  # abdominal-A — LOW
]

HLH54F_WEIGHTS: dict[str, float] = {
    "tll":    3.0,
    "hkb":    2.5,  # similar to tll — strong posterior marker
    "HLH54F": 2.0,
    "fkh":    1.5,
    "byn":    1.5,
    "hb":     1.5,
    "gt":     1.5,
    "Abd-B":  1.5,
    "Kr":     1.0,
    "abd-A":  1.0,
    "kni":    0.5,
    "Ubx":    0.5,
}

# ── FoxL1_BOTv pattern ────────────────────────────────────────────────────────
FOXL1_PATTERN: list[tuple[str, int]] = [
    ("gt",        +1),  # giant — HIGH
    ("Antp",      +1),  # Antennapedia — HIGH (mid A-P shift)
    ("Abd-B",     +1),  # Abdominal-B — HIGH
    ("cad",       +1),  # caudal — HIGH (BOT background elevates cad)
    ("prd",       +1),  # paired — HIGH (pair-rule, run-dependent)
    ("ftz",       +1),  # fushi tarazu — HIGH
    ("Dichaete",  +1),  # Dichaete (Sox gene) — HIGH in foxl1 positional shift
    ("opa",       +1),  # odd-paired — HIGH
    ("run",       +1),  # runt — HIGH in well-ventralized foxl1 (esp. nc14d+)
    ("FoxL1",     +1),  # genotype gene
    ("hb",        -1),  # hunchback — LOW
    ("tll",       -1),  # tailless — LOW (strong contrast vs HLH54F)
    ("Ubx",       -1),  # Ultrabithorax — LOW
    ("abd-A",     -1),  # abdominal-A — LOW
]

FOXL1_WEIGHTS: dict[str, float] = {
    "Kr":       3.0,    # strong marker — HIGH in foxl1 (positional shift)
    "FoxL1":    2.0,
    "tll":      2.0,    # LOW in foxl1 — good contrast
    "run":      1.5,
    "gt":       1.5,
    "Abd-B":    1.5,
    "cad":      1.5,
    "prd":      1.2,
    "ftz":      1.2,
    "Dichaete": 1.0,
    "opa":      1.0,
    "Antp":     0.8,
    "hb":       1.0,
    "Ubx":      0.5,
    "abd-A":    0.5,
}

# ── DV axis marker sets ───────────────────────────────────────────────────────
DV_VENTRAL_MARKERS  = ["dl", "sna", "twi", "NetA"]
DV_LATERAL_MARKERS  = ["sog", "brk", "pyr", "elav"]
DV_DORSAL_MARKERS   = ["Doc1", "anc", "Abd-B"]  # Abd-B marks dorsal pouch

# ── Outlier exclusions ─────────────────────────────────────────────────────────
# (group, barcode_substring) pairs — applied after DL.load(), before any
# scoring, z-matrix construction, or plotting. This is a secondary,
# script-specific exclusion and is not the primary QC filter (see
# sample_config.MANUAL_EXCLUDE / step5c_qc_filter.py for that).
EXCLUDE_EMBRYOS: list[tuple[str, str]] = [
    ("FoxL1_BOTv", "bc27"),
]

# ── Pair-rule context genes (shown in GOI violins) ────────────────────────────
PAIRRULE_CONTEXT = ["run", "eve", "H", "ftz", "odd", "prd", "slp1"]

# ── Context genes (shown in heatmaps / violins but not scored) ────────────────
CONTEXT_GENES = [
    "run", "cic", "HLH54F", "FoxL1", "bcd", "cad",
    "eve", "H",   # pair-rule context
]

# ── Full heatmap gene set (mirrors plot_heatmaps_v7_informative.py) ──────────
# Used by _draw_violins_heatmap_genes() to show the same gene complement as
# the informative heatmaps, grouped by category.
HEATMAP_GENES_BY_CATEGORY: dict[str, list[str]] = {
    "gap":       ["hb", "Kr", "kni", "gt", "tll", "hkb"],
    "pair-rule": ["run", "eve", "H", "ftz", "odd", "prd"],
    "maternal":  ["cad", "cic", "dl"],
    "mesoderm":  ["sna", "twi", "HLH54F", "NetA"],
    "TF":        ["HLH54F", "foxl1", "Oc"],
    "signaling": ["pyr"],
}
HEATMAP_GENES: list[str] = list(dict.fromkeys(
    g for genes in HEATMAP_GENES_BY_CATEGORY.values() for g in genes
))
HEATMAP_CATEGORY_COLOURS: dict[str, str] = {
    "gap":       "#6BAED6",
    "pair-rule": "#E07B54",
    "maternal":  "#9E9AC8",
    "mesoderm":  "#74C476",
    "TF":        "#F768A1",
    "signaling": "#E6AB02",
}

# ── Match thresholds ──────────────────────────────────────────────────────────
MATCH_THRESHOLD   = 0.65
PARTIAL_THRESHOLD = 0.45

# ── Colours ───────────────────────────────────────────────────────────────────
MATCH_COLOURS = {
    "match":    "#2ecc71",
    "partial":  "#f39c12",
    "mismatch": "#e74c3c",
    "WT":       "#555555",
    "other":    "#AAAAAA",
}
DIRECTION_COLOURS = {+1: "#C0392B", -1: "#2980B9"}

# Stage-aware colour palette — hue = genotype, shade = stage
# Mirrors the ATAC-seq colour scheme so figures are cross-analysis comparable.
#
#   WT                     dark grey   (nc14b only)
#   FoxL1_BOTv (BOTv)      green       nc14b=light, nc14d/gastr=dark
#   HLH54F_BOTCv (BOTCv)   pink/mag    nc14b=light, nc14d/gastr=dark
#   B6_BOTC (BOT)          blue        nc14d=medium, gastr=dark  (no nc14b)
#   BOTR                   purple      nc14b=medium, gastr=dark purple
#   runt                   red-orange  (not in ATAC scheme — kept distinct)

STAGE_AWARE_COLOURS: dict[tuple[str, str], str] = {
    # WT — only nc14b present
    ("WT",           "nc14b"): "#555555",
    ("WT",           "nc14d"): "#555555",
    ("WT",           "gastr"): "#555555",
    # FoxL1_BOTv — BOTv = green family
    ("FoxL1_BOTv",   "nc14b"): "#2ca02c",   # BOTv
    ("FoxL1_BOTv",   "nc14d"): "#1a6b1a",   # BOTv_late
    ("FoxL1_BOTv",   "gastr"): "#145214",   # BOTv_latest (extrapolated)
    # HLH54F_BOTCv — BOTCv = pink/magenta family
    ("HLH54F_BOTCv", "nc14b"): "#e377c2",   # BOTCv
    ("HLH54F_BOTCv", "nc14d"): "#b5369a",   # BOTCv_late
    ("HLH54F_BOTCv", "gastr"): "#7a1f6b",   # BOTCv_latest (extrapolated)
    # B6_BOTC — BOT = blue family (no nc14b in dataset)
    ("B6_BOTC",      "nc14b"): "#1f77b4",   # BOT (if nc14b ever added)
    ("B6_BOTC",      "nc14d"): "#1f77b4",   # BOT
    ("B6_BOTC",      "gastr"): "#145380",   # BOT_late
    # BOTR — purple family
    ("BOTR",         "nc14b"): "#9467bd",   # BOTR
    ("BOTR",         "nc14d"): "#7a50a8",   # BOTR_mid (interpolated)
    ("BOTR",         "gastr"): "#5c3585",   # BOTR_gastr
    # runt — kept red-orange (not in ATAC scheme)
    ("runt",         "nc14b"): "#d62728",
    ("runt",         "nc14d"): "#a61520",
    ("runt",         "gastr"): "#7a0e16",
}

def embryo_colour(col: str) -> str:
    """Return stage-aware colour for a single embryo column."""
    grp   = group_of.get(col, "unknown") if "group_of" in dir() else "unknown"
    stage = SC.get_stage(col)
    return STAGE_AWARE_COLOURS.get((grp, stage),
           SC.GROUP_COLOURS.get(grp, "#AAAAAA"))

# Legacy STAGE_COLOURS kept for heatmap annotation strips
STAGE_COLOURS = {"nc14b": "#FF6B6B", "nc14d": "#FFD93D", "gastr": "#6BCB77"}

GENOTYPE_TIER_COLOURS = {
    "HLH54F_BOTCv_match":    "#b5369a",   # BOTCv_late
    "HLH54F_BOTCv_partial":  "#e377c2",   # BOTCv
    "HLH54F_BOTCv_mismatch": "#f0b8e2",   # BOTCv light
    "FoxL1_BOTv_match":      "#1a6b1a",   # BOTv_late
    "FoxL1_BOTv_partial":    "#2ca02c",   # BOTv
    "FoxL1_BOTv_mismatch":   "#9ed99e",   # BOTv light
    "WT":                    "#555555",
}

DV_AXIS_COLOURS = {
    "ventral": "#C0392B",
    "lateral": "#F39C12",
    "dorsal":  "#2980B9",
}

# =============================================================================

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--match-threshold",   type=float, default=MATCH_THRESHOLD)
    p.add_argument("--partial-threshold", type=float, default=PARTIAL_THRESHOLD)
    p.add_argument("--no-all-groups",     action="store_true")
    p.add_argument("--prefix", default="all_genotypes_plus_set3",
                    help="Combined-matrix prefix, matching whatever --combined "
                         "you used for step5b/step5c/step6. Default: "
                         "all_genotypes_plus_set3")
    return p.parse_args()


args = parse_args()
MATCH_THRESHOLD   = args.match_threshold
PARTIAL_THRESHOLD = args.partial_threshold

os.makedirs(OUT_DIR, exist_ok=True)
sns.set_theme(style="white", font_scale=1.05)


def save(fig, name: str):
    fig.savefig(os.path.join(OUT_DIR, f"{name}.png"), dpi=200, bbox_inches="tight")
    fig.savefig(os.path.join(OUT_DIR, f"{name}.pdf"),           bbox_inches="tight")
    plt.close(fig)
    print(f"[posclass]   Saved → {name}.png / .pdf")


def short_label(col: str) -> str:
    import re
    grp   = group_of.get(col, "")
    stage = SC.get_stage(col)
    m     = re.search(r"(bc\d+)$", col)
    bc    = m.group(1) if m else col.split("_")[-1]
    if grp in STAGE_AWARE_GROUPS:
        return f"{grp} {bc}\n({stage})"
    return f"{grp} {bc}"


# =============================================================================
#  Data loading + symbol map
# =============================================================================

symbol_to_fbgn: dict[str, str] = {}
fbgn_to_symbol: dict[str, str] = {}

if SYMBOL_MAP and os.path.exists(SYMBOL_MAP):
    sym_df = pd.read_csv(SYMBOL_MAP, sep="\t").dropna(subset=["Symbol", "GeneID"])
    sym_df = sym_df[sym_df["Symbol"].str.strip() != ""]
    symbol_to_fbgn = dict(zip(sym_df["Symbol"], sym_df["GeneID"]))
    fbgn_to_symbol = dict(zip(sym_df["GeneID"],  sym_df["Symbol"]))
    print(f"[posclass] Loaded {len(symbol_to_fbgn):,} symbol mappings")


def resolve(symbols: list[str]) -> dict[str, str]:
    """Return {symbol: row_id} for symbols present in tpm."""
    out: dict[str, str] = {}
    for sym in symbols:
        if sym in tpm.index:
            out[sym] = sym
        elif sym in symbol_to_fbgn and symbol_to_fbgn[sym] in tpm.index:
            out[sym] = symbol_to_fbgn[sym]
        else:
            hits = [k for k in symbol_to_fbgn
                    if k.lower() == sym.lower() and symbol_to_fbgn[k] in tpm.index]
            if hits:
                out[sym] = symbol_to_fbgn[hits[0]]
                if hits[0] != sym:
                    print(f"[posclass]   '{sym}' matched via '{hits[0]}'")
            else:
                print(f"[posclass]   WARNING: '{sym}' not found — skipping")
    return out


print("[posclass] Loading data ...")
_data    = DL.load(prefix=args.prefix, use_qc=True, use_quantile=True)
tpm      = _data.tpm
tpm_qn   = _data.tpm_qn if _data.tpm_qn is not None else _data.tpm
group_of = _data.group_of
by_group = DL.embryos_by_group(_data)
all_cols = DL.ordered_cols(_data)

# ── Apply outlier exclusions (v2) ─────────────────────────────────────────────
def _is_excluded(col: str) -> bool:
    grp = group_of.get(col, "")
    return any(grp == eg and eb in col for eg, eb in EXCLUDE_EMBRYOS)

_before  = len(all_cols)
all_cols = [c for c in all_cols if not _is_excluded(c)]
_n_excl  = _before - len(all_cols)
if _n_excl:
    _excl_ids = [c for c in DL.ordered_cols(_data) if _is_excluded(c)]
    print(f"[posclass] Excluded {_n_excl} outlier embryo(s): {_excl_ids}")

wt_cols    = [c for c in by_group.get("WT",            []) if c in all_cols]
hlh_cols   = [c for c in by_group.get("HLH54F_BOTCv",  []) if c in all_cols]
foxl1_cols = [c for c in by_group.get("FoxL1_BOTv",    []) if c in all_cols]

print(f"[posclass] WT:{len(wt_cols)}  HLH54F:{len(hlh_cols)}  FoxL1:{len(foxl1_cols)}")
print(f"[posclass] All groups: {len(all_cols)} embryos")

# Collect all unique genes we need to resolve
all_pattern_genes_hlh  = [g for g, _ in HLH54F_PATTERN]
all_pattern_genes_fox  = [g for g, _ in FOXL1_PATTERN]
all_dv_genes           = DV_VENTRAL_MARKERS + DV_LATERAL_MARKERS + DV_DORSAL_MARKERS
all_genes_needed       = list(dict.fromkeys(
    all_pattern_genes_hlh + all_pattern_genes_fox +
    all_dv_genes + CONTEXT_GENES + PAIRRULE_CONTEXT
))

gene_ids      = resolve(all_genes_needed)
found_hlh_pat = [g for g, _ in HLH54F_PATTERN if g in gene_ids]
found_fox_pat = [g for g, _ in FOXL1_PATTERN   if g in gene_ids]
found_ctx     = [g for g in CONTEXT_GENES       if g in gene_ids]
found_pr      = [g for g in PAIRRULE_CONTEXT    if g in gene_ids]
found_dv      = {
    "ventral": [g for g in DV_VENTRAL_MARKERS if g in gene_ids],
    "lateral": [g for g in DV_LATERAL_MARKERS if g in gene_ids],
    "dorsal":  [g for g in DV_DORSAL_MARKERS  if g in gene_ids],
}

if not found_hlh_pat and not found_fox_pat:
    sys.exit("[posclass] No pattern genes resolved — check SYMBOL_MAP.")

print(f"[posclass] HLH54F pattern: {len(found_hlh_pat)}/{len(all_pattern_genes_hlh)} genes")
print(f"[posclass] FoxL1 pattern : {len(found_fox_pat)}/{len(all_pattern_genes_fox)} genes")
print(f"[posclass] DV genes: {sum(len(v) for v in found_dv.values())} resolved")

# =============================================================================
#  Expression matrices
# =============================================================================

all_found = list(gene_ids.keys())
log_mat = np.log2(
    tpm.loc[[gene_ids[g] for g in all_found], all_cols] + 1
)
log_mat.index = all_found

# ── Stage-stratified z-score matrix ───────────────────────────────────────────
def build_stage_stratified_zmat() -> pd.DataFrame:
    """
    For STAGE_AWARE_GROUPS: compute z-scores within stage bins when ≥2 embryos.
    For all other groups: use global z-score (existing behaviour).
    """
    base_rows = [gene_ids[g] for g in all_found]
    tpm_sub   = tpm_qn.loc[base_rows, all_cols].copy()
    tpm_sub.index = all_found

    # Start with global z-score as default
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        global_z = (tpm_sub
                    .subtract(tpm_sub.mean(axis=1), axis=0)
                    .divide(tpm_sub.std(axis=1).replace(0, np.nan), axis=0)
                    .fillna(0))

    z_out = global_z.copy()

    for stage in ("nc14b", "nc14d", "gastr"):
        stage_cols = [c for c in all_cols
                      if SC.get_stage(c) == stage
                      and group_of.get(c) in STAGE_AWARE_GROUPS]
        if len(stage_cols) < 2:
            continue  # fall back to global
        sub = tpm_sub[stage_cols]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sub_z = (sub
                     .subtract(sub.mean(axis=1), axis=0)
                     .divide(sub.std(axis=1).replace(0, np.nan), axis=0)
                     .fillna(0))
        z_out[stage_cols] = sub_z

    return z_out


z_mat = build_stage_stratified_zmat()

# =============================================================================
#  Scoring
# =============================================================================

def score_pattern(col: str,
                  pattern: list[tuple[str, int]],
                  base_weights: dict[str, float],
                  genotype_group: str,
                  z_threshold: float = 0.25) -> tuple[float, dict]:
    """
    Weighted pattern score with stage-adjusted gene weights.

    Returns (score ∈ [0,1], details dict).
    """
    stage    = SC.get_stage(col)
    w_match  = 0.0
    w_total  = 0.0
    details: dict[str, tuple | None] = {}

    for gene, direction in pattern:
        if gene not in z_mat.index:
            details[gene] = None
            continue
        z = float(z_mat.loc[gene, col])
        w = _stage_weight(gene, genotype_group, stage,
                          base_weights.get(gene, 1.0))
        correct = (z > z_threshold) if direction == +1 else (z < -z_threshold)
        details[gene] = (z, correct, direction, w)
        w_match += w * int(correct)
        w_total += w

    score = w_match / w_total if w_total > 0 else float("nan")
    return score, details


def call_match(score: float) -> str:
    if np.isnan(score):      return "unclear"
    if score >= MATCH_THRESHOLD:   return "match"
    if score >= PARTIAL_THRESHOLD: return "partial"
    return "mismatch"


# Score all embryos
hlh_scores:   dict[str, float] = {}
foxl1_scores: dict[str, float] = {}
hlh_details:  dict[str, dict]  = {}
foxl1_details:dict[str, dict]  = {}

for col in all_cols:
    grp = group_of.get(col, "unknown")
    s, d = score_pattern(col, HLH54F_PATTERN, HLH54F_WEIGHTS, "HLH54F_BOTCv")
    hlh_scores[col]   = s
    hlh_details[col]  = d
    s, d = score_pattern(col, FOXL1_PATTERN, FOXL1_WEIGHTS, "FoxL1_BOTv")
    foxl1_scores[col]  = s
    foxl1_details[col] = d

hlh_calls   = {c: call_match(hlh_scores[c])   for c in all_cols}
foxl1_calls = {c: call_match(foxl1_scores[c]) for c in all_cols}

# Print score summary
print(f"\n[posclass] Pattern scores (mean ± std):")
for grp, cols_g in by_group.items():
    if not cols_g: continue
    hv = [hlh_scores[c]   for c in cols_g if not np.isnan(hlh_scores.get(c, np.nan))]
    fv = [foxl1_scores[c] for c in cols_g if not np.isnan(foxl1_scores.get(c, np.nan))]
    print(f"  {grp:<18}  HLH54F {np.mean(hv):.2f}±{np.std(hv):.2f}  "
          f"FoxL1 {np.mean(fv):.2f}±{np.std(fv):.2f}")

# =============================================================================
#  Expression validation (mandatory gene checks)
# =============================================================================

HLH54F_EXPR_MIN = 0.5
FOXL1_RUN_MIN   = 1.0
FOXL1_GENE_MIN  = 0.3

def _log2tpm(gene_sym: str, col: str) -> float:
    if gene_sym not in gene_ids: return float("nan")
    rid = gene_ids[gene_sym]
    if col not in log_mat.columns or rid not in log_mat.index:
        return float("nan")
    v = float(log_mat.loc[rid, col]) if rid in log_mat.index else float("nan")
    return v


for col in all_cols:
    if hlh_calls[col] in ("match", "partial"):
        expr = _log2tpm("HLH54F", col)
        if not np.isnan(expr) and expr < HLH54F_EXPR_MIN:
            hlh_calls[col] = "mismatch"

    if foxl1_calls[col] in ("match", "partial"):
        run_expr   = _log2tpm("run",   col)
        foxl1_expr = _log2tpm("FoxL1", col)
        fail = False
        if not np.isnan(run_expr)   and run_expr   < FOXL1_RUN_MIN:   fail = True
        if not np.isnan(foxl1_expr) and foxl1_expr < FOXL1_GENE_MIN:  fail = True
        if fail:
            foxl1_calls[col] = "mismatch"

# =============================================================================
#  Helpers: heatmap drawing
# =============================================================================

def _draw_heatmap(cols: list[str], title: str, filename: str,
                  pattern: list[tuple[str, int]],
                  scores: dict[str, float], calls: dict[str, str],
                  show_context: bool = True):
    if not cols:
        print(f"[posclass]   No embryos for {filename} — skipping")
        return

    high_genes = [g for g, d in pattern if d == +1 and g in z_mat.index]
    low_genes  = [g for g, d in pattern if d == -1 and g in z_mat.index]
    ctx_here   = [g for g in found_ctx if g in z_mat.index] if show_context else []
    row_order  = high_genes + low_genes + ctx_here
    n_pat      = len(high_genes) + len(low_genes)

    if not (high_genes + low_genes):
        print(f"[posclass]   No pattern genes for {filename}")
        return

    # Sort columns: WT first, then target group by score desc
    wt_here  = [c for c in cols if group_of.get(c) == "WT"]
    mut_here = [c for c in cols if group_of.get(c) != "WT"]
    mut_sorted: list[str] = []
    for _grp in SC.GROUP_ORDER:
        _gc = sorted([c for c in mut_here if group_of.get(c) == _grp],
                     key=lambda c: -scores.get(c, 0))
        mut_sorted.extend(_gc)
    col_order = wt_here + mut_sorted

    z_plot   = z_mat.loc[row_order, col_order]
    log_plot = log_mat.loc[row_order, col_order]
    col_lbls = [short_label(c) for c in col_order]
    n_genes  = len(row_order)
    n_emb    = len(col_order)

    fig_w = max(10, n_emb * 0.75 + 5)
    fig_h = max(5,  n_genes * 0.55 + 3)
    fig, (ax_z, ax_log) = plt.subplots(
        1, 2, figsize=(fig_w * 2, fig_h),
        gridspec_kw={"wspace": 0.35})

    log_vmax = float(np.percentile(log_plot.values[log_plot.values > 0], 98)) \
               if (log_plot.values > 0).any() else 8

    pat_gene_set = {g for g, _ in pattern}

    for ax, data, cmap, cbar_lbl, vmin, vmax, center in [
        (ax_z,   z_plot,   "RdBu_r", "Z-score (QN-TPM)", -2.5, 2.5,  0),
        (ax_log, log_plot, "YlOrRd", "log₂(TPM+1)",       0, log_vmax, None),
    ]:
        kw = dict(cmap=cmap, xticklabels=col_lbls,
                  yticklabels=row_order, linewidths=0.4, linecolor="white",
                  cbar_kws={"label": cbar_lbl, "shrink": 0.5})
        if center is not None:
            kw.update(vmin=vmin, vmax=vmax, center=center)
        else:
            kw.update(vmin=vmin, vmax=vmax)
        sns.heatmap(data, ax=ax, **kw)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)

        for tick, gene in zip(ax.get_yticklabels(), row_order):
            tick.set_fontstyle("italic"); tick.set_fontsize(9)
            if gene in pat_gene_set:
                d = next((dd for gg, dd in pattern if gg == gene), 0)
                tick.set_color(DIRECTION_COLOURS.get(d, "black"))
            else:
                tick.set_color("#777777")

        # Direction strip on Z-score panel
        if ax is ax_z:
            ds = ax.inset_axes([-0.06, 0, 0.04, 1], transform=ax.transAxes)
            ds.set_xlim(0, 1); ds.set_ylim(0, n_genes); ds.axis("off")
            for i, gene in enumerate(row_order):
                col_d = (DIRECTION_COLOURS.get(
                    next((dd for gg, dd in pattern if gg == gene), 0), "#AAA")
                    if gene in pat_gene_set else "#DDDDDD")
                ds.add_patch(plt.Rectangle((0, n_genes - i - 1), 1, 1,
                                           color=col_d, transform=ds.transData))

        # Match-call colour strip + stage strip above heatmap
        strip = ax.inset_axes([0, 1.01, 1, 0.035], transform=ax.transAxes)
        strip.set_xlim(0, n_emb); strip.set_ylim(0, 1); strip.axis("off")
        for j, col in enumerate(col_order):
            c = (MATCH_COLOURS["WT"] if group_of.get(col) == "WT"
                 else MATCH_COLOURS.get(calls.get(col, "other"), "#AAA"))
            strip.add_patch(plt.Rectangle((j, 0), 1, 1, color=c,
                                          transform=strip.transData))

        # Stage strip (thin, above match strip) — for stage-aware groups
        stg_strip = ax.inset_axes([0, 1.05, 1, 0.02], transform=ax.transAxes)
        stg_strip.set_xlim(0, n_emb); stg_strip.set_ylim(0, 1); stg_strip.axis("off")
        for j, col in enumerate(col_order):
            stg = SC.get_stage(col)
            stg_strip.add_patch(plt.Rectangle(
                (j, 0), 1, 1,
                color=STAGE_COLOURS.get(stg, "#CCCCCC"),
                transform=stg_strip.transData))

        # Score annotations
        sc_ax = ax.inset_axes([0, 1.08, 1, 0.07], transform=ax.transAxes)
        sc_ax.set_xlim(0, n_emb); sc_ax.set_ylim(0, 1); sc_ax.axis("off")
        for j, col in enumerate(col_order):
            if group_of.get(col) != "WT":
                sv = scores.get(col, float("nan"))
                if not np.isnan(sv):
                    sc_ax.text(j + 0.5, 0.1, f"{sv:.2f}", ha="center",
                               va="bottom", fontsize=6.5,
                               color=MATCH_COLOURS.get(calls.get(col, "other"), "#888"))

        # WT boundary + HIGH/LOW boundary + context boundary
        n_wt = len(wt_here)
        if 0 < n_wt < n_emb:
            ax.axvline(n_wt, color="black", lw=2.0)
        if 0 < len(high_genes) < n_pat:
            ax.axhline(len(high_genes), color="black", lw=1.5, ls="--")
        if ctx_here and n_pat > 0:
            ax.axhline(n_pat, color="#666", lw=1.5, ls=":")

    # Group header
    hdr = ax_z.inset_axes([0, 1.16, 1, 0.07], transform=ax_z.transAxes)
    hdr.set_xlim(0, n_emb); hdr.set_ylim(0, 1); hdr.axis("off")
    if wt_here:
        hdr.text(len(wt_here) / 2, 0.1, "WT", ha="center", va="bottom",
                 fontsize=9, fontweight="bold",
                 color=MATCH_COLOURS["WT"])
    _pos = len(wt_here)
    for _grp in SC.GROUP_ORDER:
        _gc = [c for c in mut_sorted if group_of.get(c) == _grp]
        if _gc:
            hdr.text(_pos + len(_gc) / 2, 0.1, _grp, ha="center", va="bottom",
                     fontsize=8, fontweight="bold",
                     color=STAGE_AWARE_COLOURS.get((_grp, "nc14b"), SC.GROUP_COLOURS.get(_grp, "#888")))
            _pos += len(_gc)

    # Legend
    leg_handles = [
        mpatches.Patch(color=MATCH_COLOURS["match"],   label=f"match (≥{MATCH_THRESHOLD:.0%})"),
        mpatches.Patch(color=MATCH_COLOURS["partial"], label=f"partial (≥{PARTIAL_THRESHOLD:.0%})"),
        mpatches.Patch(color=MATCH_COLOURS["mismatch"],label="mismatch"),
        mpatches.Patch(color=MATCH_COLOURS["WT"],      label="WT"),
        mpatches.Patch(color=DIRECTION_COLOURS[+1],    label="expected HIGH"),
        mpatches.Patch(color=DIRECTION_COLOURS[-1],    label="expected LOW"),
        mpatches.Patch(color="#DDDDDD",                label="context gene"),
        mpatches.Patch(color=STAGE_COLOURS["nc14b"],   label="nc14b"),
        mpatches.Patch(color=STAGE_COLOURS["nc14d"],   label="nc14d"),
        mpatches.Patch(color=STAGE_COLOURS["gastr"],   label="gastr"),
    ]
    ax_z.legend(handles=leg_handles, bbox_to_anchor=(1.02, 1.0),
                loc="upper left", frameon=False, fontsize=7.5, title="legend")

    ctx_note = (f"\nContext genes (grey italic): {', '.join(ctx_here)}"
                if ctx_here else "")
    fig.suptitle(
        f"{title}  [{len(col_order)} embryos]\n"
        f"Left: Z-scored QN-TPM  |  Right: log₂(TPM+1)\n"
        f"Top strip: match call | stage (red=nc14b, yellow=nc14d, green=gastr)"
        f"{ctx_note}",
        fontsize=9, y=1.02)
    save(fig, filename)


# =============================================================================
#  Helper: violin plots (pattern genes + pair-rule context)
# =============================================================================

def _draw_violins(cols: list[str], title: str, filename: str,
                  pattern: list[tuple[str, int]],
                  show_context: bool = True,
                  show_pairrule: bool = True,
                  exclude_wt: bool = False):
    if exclude_wt:
        cols = [c for c in cols if group_of.get(c) != "WT"]
    pat_genes_present = [g for g, _ in pattern if g in z_mat.index]
    ctx_present       = [g for g in found_ctx if g in z_mat.index] if show_context else []
    pr_present        = [g for g in found_pr  if g in z_mat.index and
                         g not in pat_genes_present] if show_pairrule else []
    all_plot          = pat_genes_present + ctx_present + pr_present
    if not pat_genes_present or not cols:
        return

    pat_gene_set = {g for g, _ in pattern}
    pr_gene_set  = set(PAIRRULE_CONTEXT)

    records = []
    for col in cols:
        grp   = group_of.get(col, "unknown")
        stage = SC.get_stage(col)
        for gene in all_plot:
            direction = next((d for g, d in pattern if g == gene), None)
            is_ctx    = gene in ctx_present and gene not in pat_gene_set
            is_pr     = gene in pr_gene_set  and gene not in pat_gene_set
            records.append({
                "embryo":     col,
                "label":      short_label(col),
                "group":      grp,
                "stage":      stage,
                "gene":       gene,
                "log2TPM":    float(log_mat.loc[gene, col]),
                "direction":  ("HIGH" if direction == +1 else
                               "LOW"  if direction == -1 else
                               "context" if is_ctx else "pair-rule"),
                "section":    ("pair-rule" if is_pr else
                               "context"   if is_ctx else
                               "pattern"),
            })
    df = pd.DataFrame(records)

    # Order: pattern HIGH → pattern LOW → context → pair-rule
    high_g  = [g for g, d in pattern if d == +1 and g in pat_genes_present]
    low_g   = [g for g, d in pattern if d == -1 and g in pat_genes_present]
    gene_order = high_g + low_g + ctx_present + pr_present

    groups_present = [g for g in SC.GROUP_ORDER if g in df["group"].unique()]
    palette        = {g: STAGE_AWARE_COLOURS.get((g, "nc14b"),
                           SC.GROUP_COLOURS.get(g, "#888"))
                       for g in groups_present}

    ncols_v = min(5, len(gene_order))
    nrows_v = int(np.ceil(len(gene_order) / ncols_v))
    fig, axes = plt.subplots(nrows_v, ncols_v,
                             figsize=(4.5 * ncols_v, 4.2 * nrows_v), sharey=False)
    af = np.array(axes).flatten() if len(gene_order) > 1 else [axes]

    SECTION_BG = {
        "HIGH":       "#FFF5F5",
        "LOW":        "#F0F5FF",
        "context":    "#F8F8F8",
        "pair-rule":  "#F5FFF5",
    }

    for ax, gene in zip(af[:len(gene_order)], gene_order):
        sub       = df[df["gene"] == gene]
        d_txt     = sub["direction"].iloc[0] if len(sub) else "context"
        dir_col   = (DIRECTION_COLOURS[+1] if d_txt == "HIGH" else
                     DIRECTION_COLOURS[-1] if d_txt == "LOW" else "#777777")
        bg_col    = SECTION_BG.get(d_txt, "#F8F8F8")

        sns.violinplot(data=sub, x="group", y="log2TPM",
                       hue="group", order=groups_present, palette=palette,
                       inner=None, cut=0, linewidth=0.8, legend=False, ax=ax)
        # Per-point scatter — stage-aware colour per embryo
        for _ec in [c for c in all_cols if group_of.get(c, "") in groups_present]:
            _xi  = groups_present.index(group_of.get(_ec, ""))
            _xj  = _xi + float(np.random.uniform(-0.15, 0.15))
            if gene not in log_mat.index: continue
            ax.scatter(_xj, float(log_mat.loc[gene, _ec]),
                       c=embryo_colour(_ec), s=45, alpha=0.9,
                       zorder=3, edgecolors="white", linewidths=0.4)

        # Stage markers: overlay with different shapes for nc14d/gastr
        for stage, marker, sz in [("nc14d", "s", 35), ("gastr", "^", 35)]:
            stg_sub = sub[sub["stage"] == stage]
            if not stg_sub.empty:
                x_positions = [groups_present.index(g)
                                if g in groups_present else 0
                                for g in stg_sub["group"]]
                x_jitter = np.array(x_positions, dtype=float) + \
                           np.random.uniform(-0.15, 0.15, len(stg_sub))
                ax.scatter(x_jitter, stg_sub["log2TPM"].values,
                           marker=marker, s=sz, c="black", alpha=0.5,
                           zorder=4, linewidths=0)

        # WT median reference
        wt_sub = sub[sub["group"] == "WT"]
        if not wt_sub.empty:
            wt_med = wt_sub["log2TPM"].median()
            ax.axhline(wt_med, color=MATCH_COLOURS["WT"],
                       lw=1.5, ls="--", alpha=0.7, label=f"WT med={wt_med:.1f}")

        ax.set_facecolor(bg_col)
        ax.set_title(f"{gene}  [{d_txt}]", fontsize=10,
                     fontstyle="italic", fontweight="bold", color=dir_col)
        ax.set_xlabel("")
        ax.set_ylabel("log₂(TPM+1)", fontsize=8)
        ax.set_xticks(range(len(groups_present)))
        ax.set_xticklabels(groups_present, rotation=35, ha="right", fontsize=7.5)
        ax.spines[["top", "right"]].set_visible(False)
        if not wt_sub.empty:
            ax.legend(frameon=False, fontsize=6.5, loc="upper right")

    for ax in af[len(gene_order):]:
        ax.set_visible(False)

    # Legend patches
    # Stage-aware legend: show nc14b and late-stage colour per group
    grp_leg = []
    for g in groups_present:
        c_early = STAGE_AWARE_COLOURS.get((g,"nc14b"), SC.GROUP_COLOURS.get(g,"#888"))
        c_late  = STAGE_AWARE_COLOURS.get((g,"nc14d"), c_early)
        grp_leg.append(mpatches.Patch(color=c_early, label=f"{g} nc14b"))
        if c_late != c_early:
            grp_leg.append(mpatches.Patch(color=c_late, label=f"{g} nc14d/gastr"))
    from matplotlib.lines import Line2D
    stage_leg = [
        Line2D([0],[0], marker="s", linestyle="", color="black",
               alpha=0.5, markersize=6, label="nc14d (■)"),
        Line2D([0],[0], marker="^", linestyle="", color="black",
               alpha=0.5, markersize=6, label="gastr (▲)"),
    ]
    fig.legend(handles=grp_leg + stage_leg, title="group / stage",
               loc="lower right", ncol=3, frameon=False, fontsize=8)

    fig.suptitle(
        f"{title}\n"
        f"Pattern: red=HIGH, blue=LOW | Context: grey | Pair-rule: green\n"
        f"Stage: circle=nc14b, ■=nc14d, ▲=gastr  |  Dashed = WT median",
        fontsize=9, y=1.02)
    plt.tight_layout()
    save(fig, filename)


# =============================================================================
#  NEW: DV axis violin plot — ventral / lateral / dorsal per group
# =============================================================================

def _draw_dv_violins():
    """Full DV axis violin grid — one panel per gene, coloured by group."""
    all_dv_symbols = (found_dv["ventral"] + found_dv["lateral"] +
                      found_dv["dorsal"])
    if not all_dv_symbols:
        print("[posclass]   No DV genes resolved — skipping DV violin")
        return

    records = []
    for col in all_cols:
        grp   = group_of.get(col, "unknown")
        stage = SC.get_stage(col)
        for gene in all_dv_symbols:
            axis_cat = ("ventral" if gene in found_dv["ventral"] else
                        "lateral" if gene in found_dv["lateral"] else "dorsal")
            records.append({
                "embryo":   col,
                "group":    grp,
                "stage":    stage,
                "gene":     gene,
                "DV_axis":  axis_cat,
                "log2TPM":  float(log_mat.loc[gene, col]),
            })
    df = pd.DataFrame(records)

    groups_present = [g for g in SC.GROUP_ORDER if g in df["group"].unique()]
    palette        = {g: STAGE_AWARE_COLOURS.get((g, "nc14b"),
                           SC.GROUP_COLOURS.get(g, "#888"))
                       for g in groups_present}

    ncols_dv = min(5, len(all_dv_symbols))
    nrows_dv = int(np.ceil(len(all_dv_symbols) / ncols_dv))
    fig, axes = plt.subplots(nrows_dv, ncols_dv,
                             figsize=(4.5 * ncols_dv, 4.2 * nrows_dv), sharey=False)
    af = np.array(axes).flatten() if len(all_dv_symbols) > 1 else [axes]

    for ax, gene in zip(af[:len(all_dv_symbols)], all_dv_symbols):
        sub      = df[df["gene"] == gene]
        axis_cat = sub["DV_axis"].iloc[0] if len(sub) else "ventral"
        dir_col  = DV_AXIS_COLOURS.get(axis_cat, "#888")
        bg_col   = {"ventral": "#FFF5F5", "lateral": "#FFFBF0",
                    "dorsal":  "#F0F5FF"}.get(axis_cat, "#F8F8F8")

        sns.violinplot(data=sub, x="group", y="log2TPM",
                       hue="group", order=groups_present, palette=palette,
                       inner=None, cut=0, linewidth=0.8, legend=False, ax=ax)
        # Per-point scatter — stage-aware colour per embryo
        for _ec in [c for c in all_cols if group_of.get(c, "") in groups_present]:
            _xi  = groups_present.index(group_of.get(_ec, ""))
            _xj  = _xi + float(np.random.uniform(-0.15, 0.15))
            if gene not in log_mat.index: continue
            ax.scatter(_xj, float(log_mat.loc[gene, _ec]),
                       c=embryo_colour(_ec), s=45, alpha=0.9,
                       zorder=3, edgecolors="white", linewidths=0.4)

        # Stage shapes
        for stage, marker, sz in [("nc14d", "s", 35), ("gastr", "^", 35)]:
            stg_sub = sub[sub["stage"] == stage]
            if not stg_sub.empty:
                xp = [groups_present.index(g) if g in groups_present else 0
                      for g in stg_sub["group"]]
                xj = np.array(xp, float) + np.random.uniform(-0.15, 0.15, len(stg_sub))
                ax.scatter(xj, stg_sub["log2TPM"].values,
                           marker=marker, s=sz, c="black", alpha=0.5,
                           zorder=4, linewidths=0)

        wt_sub = sub[sub["group"] == "WT"]
        if not wt_sub.empty:
            ax.axhline(wt_sub["log2TPM"].median(),
                       color=MATCH_COLOURS["WT"],
                       lw=1.5, ls="--", alpha=0.7)

        ax.set_facecolor(bg_col)
        ax.set_title(f"{gene}\n[{axis_cat}]", fontsize=10,
                     fontstyle="italic", fontweight="bold", color=dir_col)
        ax.set_xlabel("")
        ax.set_ylabel("log₂(TPM+1)", fontsize=8)
        ax.set_xticks(range(len(groups_present)))
        ax.set_xticklabels(groups_present, rotation=35, ha="right", fontsize=7.5)
        ax.spines[["top", "right"]].set_visible(False)

    for ax in af[len(all_dv_symbols):]:
        ax.set_visible(False)

    # Stage-aware legend: show nc14b and late-stage colour per group
    grp_leg = []
    for g in groups_present:
        c_early = STAGE_AWARE_COLOURS.get((g,"nc14b"), SC.GROUP_COLOURS.get(g,"#888"))
        c_late  = STAGE_AWARE_COLOURS.get((g,"nc14d"), c_early)
        grp_leg.append(mpatches.Patch(color=c_early, label=f"{g} nc14b"))
        if c_late != c_early:
            grp_leg.append(mpatches.Patch(color=c_late, label=f"{g} nc14d/gastr"))
    axis_leg = [mpatches.Patch(color=DV_AXIS_COLOURS[a], label=a)
                for a in ("ventral", "lateral", "dorsal")]
    fig.legend(handles=grp_leg + axis_leg, title="group / DV axis",
               loc="lower right", ncol=3, frameon=False, fontsize=8)

    fig.suptitle(
        "D/V axis marker genes — per-group expression\n"
        "Red=ventral (dl sna twi NetA) | Orange=lateral (sog brk pyr elav) | "
        "Blue=dorsal (Doc1 anc Abd-B)\n"
        "Stage: circle=nc14b, ■=nc14d, ▲=gastr  |  Dashed = WT median",
        fontsize=9, y=1.02)
    plt.tight_layout()
    save(fig, "violin_DV_axis")


# =============================================================================
#  NEW: Stage-split violin — x-axis = group×stage, one panel per gene
# =============================================================================

# Stage ordering for x-axis
STAGE_ORDER = ["nc14b", "nc14d", "gastr"]

def _stage_group_label(grp: str, stage: str) -> str:
    """Short x-axis label combining group abbreviation and stage."""
    abbrev = {
        "WT":           "WT",
        "runt":         "runt",
        "BOTR":         "BOTR",
        "FoxL1_BOTv":   "Fox",
        "HLH54F_BOTCv": "HLH",
        "B6_BOTC":      "B6",
    }.get(grp, grp[:4])
    stage_abbrev = {"nc14b": "b", "nc14d": "d", "gastr": "g"}.get(stage, stage)
    return f"{abbrev}\n{stage_abbrev}"


def _draw_violins_by_stage(cols: list[str], title: str, filename: str,
                            pattern: list[tuple[str, int]],
                            show_context: bool = True,
                            show_pairrule: bool = True,
                            exclude_wt: bool = False):
    if exclude_wt:
        cols = [c for c in cols if group_of.get(c) != "WT"]
    """
    Stage-split violin plot.

    X-axis: each group × stage combination that has ≥1 embryo.
    Each panel shows one gene. Violins are coloured by group; stage is
    encoded on the x-axis label.  WT median reference line shown.
    Pair-rule context genes included in a separate section.
    """
    pat_genes_present = [g for g, _ in pattern if g in z_mat.index]
    ctx_present = [g for g in found_ctx if g in z_mat.index] if show_context else []
    pr_present  = [g for g in found_pr  if g in z_mat.index and
                   g not in pat_genes_present] if show_pairrule else []
    all_plot = pat_genes_present + ctx_present + pr_present
    if not pat_genes_present or not cols:
        return

    pat_gene_set = {g for g, _ in pattern}
    pr_gene_set  = set(PAIRRULE_CONTEXT)

    # Build records
    records = []
    for col in cols:
        grp   = group_of.get(col, "unknown")
        stage = SC.get_stage(col)
        for gene in all_plot:
            direction = next((d for g, d in pattern if g == gene), None)
            is_ctx = gene in ctx_present and gene not in pat_gene_set
            is_pr  = gene in pr_gene_set  and gene not in pat_gene_set
            records.append({
                "embryo":    col,
                "group":     grp,
                "stage":     stage,
                "gs_label":  _stage_group_label(grp, stage),
                "gene":      gene,
                "log2TPM":   float(log_mat.loc[gene, col]),
                "direction": ("HIGH" if direction == +1 else
                              "LOW"  if direction == -1 else
                              "context" if is_ctx else "pair-rule"),
                "section":   ("pair-rule" if is_pr else
                              "context"   if is_ctx else "pattern"),
            })
    df = pd.DataFrame(records)

    # Build ordered x-axis: group order × stage order, keep only combos present
    groups_present = [g for g in SC.GROUP_ORDER if g in df["group"].unique()]
    xs_ordered: list[tuple[str, str]] = []  # (group, stage) pairs
    for grp in groups_present:
        for stage in STAGE_ORDER:
            if not df[(df["group"] == grp) & (df["stage"] == stage)].empty:
                xs_ordered.append((grp, stage))

    xs_labels = [_stage_group_label(g, s) for g, s in xs_ordered]

    # Colour map: use group colour for each xs position
    xs_palette = {_stage_group_label(g, s):
                  STAGE_AWARE_COLOURS.get((g, s), SC.GROUP_COLOURS.get(g, "#888"))
                  for g, s in xs_ordered}
    # Alpha per stage (lighter for later stages)
    xs_alpha = {_stage_group_label(g, s): SC.STAGE_ALPHA.get(s, 0.7)
                for g, s in xs_ordered}

    high_g = [g for g, d in pattern if d == +1 and g in pat_genes_present]
    low_g  = [g for g, d in pattern if d == -1 and g in pat_genes_present]
    gene_order = high_g + low_g + ctx_present + pr_present

    SECTION_BG = {
        "HIGH":      "#FFF5F5",
        "LOW":       "#F0F5FF",
        "context":   "#F8F8F8",
        "pair-rule": "#F5FFF5",
    }

    ncols_v = min(5, len(gene_order))
    nrows_v = int(np.ceil(len(gene_order) / ncols_v))
    fig, axes = plt.subplots(nrows_v, ncols_v,
                             figsize=(max(6, len(xs_ordered) * 0.9 + 1) * ncols_v,
                                      4.5 * nrows_v),
                             sharey=False)
    af = np.array(axes).flatten() if len(gene_order) > 1 else [axes]

    for ax, gene in zip(af[:len(gene_order)], gene_order):
        sub      = df[df["gene"] == gene]
        d_txt    = sub["direction"].iloc[0] if len(sub) else "context"
        dir_col  = (DIRECTION_COLOURS[+1] if d_txt == "HIGH" else
                    DIRECTION_COLOURS[-1] if d_txt == "LOW" else "#777777")
        bg_col   = SECTION_BG.get(d_txt, "#F8F8F8")

        # Only plot xs that have data for this gene
        xs_here = [lbl for lbl in xs_labels if lbl in sub["gs_label"].values]
        sub_here = sub[sub["gs_label"].isin(xs_here)]

        if sub_here.empty:
            ax.set_visible(False)
            continue

        # Violin — need ≥2 points per category for a violin shape
        xs_counts = sub_here.groupby("gs_label").size()
        xs_violin = [x for x in xs_here if xs_counts.get(x, 0) >= 2]
        xs_strip  = xs_here  # all get strip points

        # Draw violins using full xs_here as order so positions stay aligned
        # with the scatter points below. Groups with n<2 are skipped (no violin)
        # but their scatter points still appear at the correct x position.
        if xs_violin:
            sns.violinplot(data=sub_here[sub_here["gs_label"].isin(xs_violin)],
                           x="gs_label", y="log2TPM",
                           hue="gs_label", order=xs_here,
                           palette=xs_palette,
                           inner=None, cut=0, linewidth=0.7,
                           legend=False, ax=ax)

        # Strip — all points regardless of count
        for j, lbl in enumerate(xs_here):
            pts = sub_here[sub_here["gs_label"] == lbl]["log2TPM"].values
            if len(pts) == 0:
                continue
            jitter = np.random.uniform(-0.18, 0.18, len(pts))
            ax.scatter(np.full(len(pts), xs_here.index(lbl)) + jitter,
                       pts,
                       color=xs_palette.get(lbl, "#888"),
                       alpha=xs_alpha.get(lbl, 0.8),
                       s=45, zorder=4, edgecolors="white", linewidths=0.4)

        # WT median reference
        wt_pts = sub_here[sub_here["group"] == "WT"]["log2TPM"]
        if not wt_pts.empty:
            ax.axhline(wt_pts.median(),
                       color=MATCH_COLOURS["WT"],
                       lw=1.5, ls="--", alpha=0.6,
                       label=f"WT med={wt_pts.median():.1f}")
            ax.legend(frameon=False, fontsize=6.5, loc="upper right")

        # Group boundary lines (vertical) between groups
        prev_grp = None
        for j, (g, s) in enumerate(xs_ordered):
            if prev_grp is not None and g != prev_grp:
                ax.axvline(j - 0.5, color="#AAAAAA", lw=1.0, ls=":")
            prev_grp = g

        ax.set_facecolor(bg_col)
        ax.set_title(f"{gene}  [{d_txt}]", fontsize=9.5,
                     fontstyle="italic", fontweight="bold", color=dir_col)
        ax.set_xlabel("")
        ax.set_ylabel("log₂(TPM+1)", fontsize=8)
        ax.set_xticks(range(len(xs_here)))
        ax.set_xticklabels(xs_here, fontsize=7.5)
        ax.spines[["top", "right"]].set_visible(False)

    for ax in af[len(gene_order):]:
        ax.set_visible(False)

    # Legend
    # Stage-aware legend: show nc14b and late-stage colour per group
    grp_leg = []
    for g in groups_present:
        c_early = STAGE_AWARE_COLOURS.get((g,"nc14b"), SC.GROUP_COLOURS.get(g,"#888"))
        c_late  = STAGE_AWARE_COLOURS.get((g,"nc14d"), c_early)
        grp_leg.append(mpatches.Patch(color=c_early, label=f"{g} nc14b"))
        if c_late != c_early:
            grp_leg.append(mpatches.Patch(color=c_late, label=f"{g} nc14d/gastr"))
    from matplotlib.lines import Line2D
    stage_leg = [
        mpatches.Patch(color="#DDDDDD", alpha=1.0, label="nc14b (full)"),
        mpatches.Patch(color="#DDDDDD", alpha=0.7, label="nc14d (mid)"),
        mpatches.Patch(color="#DDDDDD", alpha=0.45,label="gastr (faint)"),
    ]
    fig.legend(handles=grp_leg + stage_leg,
               title="group / stage", loc="lower right",
               ncol=3, frameon=False, fontsize=8)

    fig.suptitle(
        f"{title}  [stage-split x-axis]\n"
        f"X = group × stage  |  red=HIGH, blue=LOW  |  grey=context  |  "
        f"green=pair-rule  |  dashed=WT median",
        fontsize=9, y=1.02)
    plt.tight_layout()
    save(fig, filename)




# =============================================================================
#  Heatmap-gene violin plots — full gene complement matching informative heatmap
# =============================================================================

def _draw_violins_heatmap_genes(cols: list[str], title: str, filename: str,
                                 exclude_wt: bool = False):
    """
    Violin plot showing every gene in HEATMAP_GENES (mirrors the informative
    heatmap gene set), coloured by category strip colour and grouped by
    HEATMAP_GENES_BY_CATEGORY.  One panel per gene, laid out in category order.

    exclude_wt: if True, WT embryos are dropped before plotting.
    """
    if exclude_wt:
        cols = [c for c in cols if group_of.get(c) != "WT"]
    if not cols:
        return

    # Resolve which heatmap genes are present in log_mat
    genes_present = [g for g in HEATMAP_GENES if g in log_mat.index]
    if not genes_present:
        print(f"[posclass]   No heatmap genes in log_mat — skipping {filename}")
        return

    # Build category lookup for colouring
    cat_of_gene: dict[str, str] = {}
    for cat, genes in HEATMAP_GENES_BY_CATEGORY.items():
        for g in genes:
            if g not in cat_of_gene:
                cat_of_gene[g] = cat

    # Build records
    records = []
    for col in cols:
        grp   = group_of.get(col, "unknown")
        stage = SC.get_stage(col)
        for gene in genes_present:
            records.append({
                "embryo":   col,
                "group":    grp,
                "stage":    stage,
                "gene":     gene,
                "category": cat_of_gene.get(gene, "other"),
                "log2TPM":  float(log_mat.loc[gene, col]),
            })
    df = pd.DataFrame(records)

    groups_present = [g for g in SC.GROUP_ORDER if g in df["group"].unique()]
    palette = {g: STAGE_AWARE_COLOURS.get((g, "nc14b"),
                   SC.GROUP_COLOURS.get(g, "#888"))
               for g in groups_present}

    ncols_v = min(6, len(genes_present))
    nrows_v = int(np.ceil(len(genes_present) / ncols_v))
    fig, axes = plt.subplots(nrows_v, ncols_v,
                             figsize=(4.2 * ncols_v, 4.0 * nrows_v),
                             sharey=False)
    af = np.array(axes).flatten() if len(genes_present) > 1 else [axes]

    wt_medians: dict[str, float] = {}
    for gene in genes_present:
        wt_vals = df[(df["gene"] == gene) & (df["group"] == "WT")]["log2TPM"]
        if not wt_vals.empty:
            wt_medians[gene] = float(wt_vals.median())

    for ax, gene in zip(af[:len(genes_present)], genes_present):
        sub     = df[df["gene"] == gene]
        cat     = cat_of_gene.get(gene, "other")
        cat_col = HEATMAP_CATEGORY_COLOURS.get(cat, "#888888")
        bg_col  = cat_col + "18"   # very light tint of category colour

        sns.violinplot(data=sub, x="group", y="log2TPM",
                       hue="group", order=groups_present, palette=palette,
                       inner=None, cut=0, linewidth=0.8, legend=False, ax=ax)

        # Per-embryo scatter with stage-aware colour
        for col in [c for c in cols if group_of.get(c, "") in groups_present]:
            xi  = groups_present.index(group_of.get(col, ""))
            xj  = xi + float(np.random.uniform(-0.15, 0.15))
            ax.scatter(xj, float(log_mat.loc[gene, col]),
                       c=embryo_colour(col), s=40, alpha=0.85,
                       zorder=3, edgecolors="white", linewidths=0.4)

        # Stage shape overlays
        for stage, marker, sz in [("nc14d", "s", 32), ("gastr", "^", 32)]:
            stg_sub = sub[sub["stage"] == stage]
            if not stg_sub.empty:
                xp = [groups_present.index(g) if g in groups_present else 0
                      for g in stg_sub["group"]]
                xj = np.array(xp, float) + np.random.uniform(-0.15, 0.15, len(stg_sub))
                ax.scatter(xj, stg_sub["log2TPM"].values,
                           marker=marker, s=sz, c="black", alpha=0.45,
                           zorder=4, linewidths=0)

        # WT median reference (only when WT is included)
        if gene in wt_medians and not exclude_wt:
            ax.axhline(wt_medians[gene], color=MATCH_COLOURS["WT"],
                       lw=1.4, ls="--", alpha=0.7,
                       label=f"WT med={wt_medians[gene]:.1f}")
            ax.legend(frameon=False, fontsize=6, loc="upper right")

        ax.set_facecolor(bg_col)
        ax.set_title(f"{gene}", fontsize=10, fontstyle="italic",
                     fontweight="bold", color=cat_col)
        ax.text(0.02, 0.98, cat, transform=ax.transAxes,
                fontsize=6.5, color=cat_col, va="top", alpha=0.8)
        ax.set_xlabel("")
        ax.set_ylabel("log₂(TPM+1)", fontsize=8)
        ax.set_xticks(range(len(groups_present)))
        ax.set_xticklabels(groups_present, rotation=35, ha="right", fontsize=7.5)
        ax.spines[["top", "right"]].set_visible(False)

    for ax in af[len(genes_present):]:
        ax.set_visible(False)

    # Legend
    cat_leg = [mpatches.Patch(color=HEATMAP_CATEGORY_COLOURS.get(c, "#888"), label=c)
               for c in HEATMAP_GENES_BY_CATEGORY]
    grp_leg = []
    for g in groups_present:
        c_early = STAGE_AWARE_COLOURS.get((g, "nc14b"), SC.GROUP_COLOURS.get(g, "#888"))
        c_late  = STAGE_AWARE_COLOURS.get((g, "nc14d"), c_early)
        grp_leg.append(mpatches.Patch(color=c_early, label=f"{g} nc14b"))
        if c_late != c_early:
            grp_leg.append(mpatches.Patch(color=c_late, label=f"{g} nc14d/gastr"))
    from matplotlib.lines import Line2D
    stage_leg = [
        Line2D([0],[0], marker="s", linestyle="", color="black",
               alpha=0.5, markersize=6, label="nc14d (■)"),
        Line2D([0],[0], marker="^", linestyle="", color="black",
               alpha=0.5, markersize=6, label="gastr (▲)"),
    ]
    fig.legend(handles=cat_leg + grp_leg + stage_leg,
               title="category / group / stage",
               loc="lower right", ncol=3, frameon=False, fontsize=8)

    no_wt_note = " [no WT]" if exclude_wt else ""
    fig.suptitle(
        f"{title}{no_wt_note}\n"
        f"Full heatmap gene set — coloured by category  |  "
        f"Stage: circle=nc14b, ■=nc14d, ▲=gastr"
        + ("  |  dashed = WT median" if not exclude_wt else ""),
        fontsize=9, y=1.02)
    plt.tight_layout()
    save(fig, filename)


def _draw_violins_heatmap_genes_by_stage(cols: list[str], title: str, filename: str,
                                          exclude_wt: bool = False):
    """
    Stage-split version of _draw_violins_heatmap_genes.
    X-axis = group × stage combinations present in data (same layout as
    _draw_violins_by_stage).  One panel per heatmap gene, category order.
    """
    if exclude_wt:
        cols = [c for c in cols if group_of.get(c) != "WT"]
    if not cols:
        return

    genes_present = [g for g in HEATMAP_GENES if g in log_mat.index]
    if not genes_present:
        print(f"[posclass]   No heatmap genes in log_mat — skipping {filename}")
        return

    cat_of_gene: dict[str, str] = {}
    for cat, genes in HEATMAP_GENES_BY_CATEGORY.items():
        for g in genes:
            if g not in cat_of_gene:
                cat_of_gene[g] = cat

    # Build records
    records = []
    for col in cols:
        grp   = group_of.get(col, "unknown")
        stage = SC.get_stage(col)
        for gene in genes_present:
            records.append({
                "embryo":    col,
                "group":     grp,
                "stage":     stage,
                "gs_label":  _stage_group_label(grp, stage),
                "gene":      gene,
                "category":  cat_of_gene.get(gene, "other"),
                "log2TPM":   float(log_mat.loc[gene, col]),
            })
    df = pd.DataFrame(records)

    # Build ordered x-axis: group × stage, only combos present
    groups_present = [g for g in SC.GROUP_ORDER if g in df["group"].unique()]
    xs_ordered: list[tuple[str, str]] = []
    for grp in groups_present:
        for stage in STAGE_ORDER:
            if not df[(df["group"] == grp) & (df["stage"] == stage)].empty:
                xs_ordered.append((grp, stage))
    xs_labels = [_stage_group_label(g, s) for g, s in xs_ordered]

    xs_palette = {_stage_group_label(g, s):
                  STAGE_AWARE_COLOURS.get((g, s), SC.GROUP_COLOURS.get(g, "#888"))
                  for g, s in xs_ordered}
    xs_alpha   = {_stage_group_label(g, s): SC.STAGE_ALPHA.get(s, 0.7)
                  for g, s in xs_ordered}

    # WT medians (computed before exclude_wt filtering, so reference is valid)
    wt_medians: dict[str, float] = {}
    if not exclude_wt:
        for gene in genes_present:
            wt_vals = df[(df["gene"] == gene) & (df["group"] == "WT")]["log2TPM"]
            if not wt_vals.empty:
                wt_medians[gene] = float(wt_vals.median())

    ncols_v = min(6, len(genes_present))
    nrows_v = int(np.ceil(len(genes_present) / ncols_v))
    fig, axes = plt.subplots(nrows_v, ncols_v,
                             figsize=(max(5, len(xs_ordered) * 0.85 + 1) * ncols_v,
                                      4.2 * nrows_v),
                             sharey=False)
    af = np.array(axes).flatten() if len(genes_present) > 1 else [axes]

    for ax, gene in zip(af[:len(genes_present)], genes_present):
        sub     = df[df["gene"] == gene]
        cat     = cat_of_gene.get(gene, "other")
        cat_col = HEATMAP_CATEGORY_COLOURS.get(cat, "#888888")
        bg_col  = cat_col + "18"

        xs_here = [lbl for lbl in xs_labels if lbl in sub["gs_label"].values]
        sub_here = sub[sub["gs_label"].isin(xs_here)]
        if sub_here.empty:
            ax.set_visible(False)
            continue

        xs_violin = [x for x in xs_here
                     if sub_here.groupby("gs_label").size().get(x, 0) >= 2]
        if xs_violin:
            sns.violinplot(data=sub_here[sub_here["gs_label"].isin(xs_violin)],
                           x="gs_label", y="log2TPM",
                           hue="gs_label", order=xs_here,
                           palette=xs_palette,
                           inner=None, cut=0, linewidth=0.7,
                           legend=False, ax=ax)

        # Strip scatter for all points
        for j, lbl in enumerate(xs_here):
            pts = sub_here[sub_here["gs_label"] == lbl]["log2TPM"].values
            if len(pts) == 0:
                continue
            jitter = np.random.uniform(-0.18, 0.18, len(pts))
            ax.scatter(np.full(len(pts), j) + jitter, pts,
                       color=xs_palette.get(lbl, "#888"),
                       alpha=xs_alpha.get(lbl, 0.8),
                       s=40, zorder=4, edgecolors="white", linewidths=0.4)

        # WT median reference
        if gene in wt_medians:
            ax.axhline(wt_medians[gene], color=MATCH_COLOURS["WT"],
                       lw=1.4, ls="--", alpha=0.7,
                       label=f"WT med={wt_medians[gene]:.1f}")
            ax.legend(frameon=False, fontsize=6, loc="upper right")

        # Group boundary lines
        prev_grp = None
        for j, (g, s) in enumerate(xs_ordered):
            if prev_grp is not None and g != prev_grp:
                if _stage_group_label(g, s) in xs_here:
                    ax.axvline(xs_here.index(_stage_group_label(g, s)) - 0.5,
                               color="#AAAAAA", lw=1.0, ls=":")
            prev_grp = g

        ax.set_facecolor(bg_col)
        ax.set_title(f"{gene}", fontsize=10, fontstyle="italic",
                     fontweight="bold", color=cat_col)
        ax.text(0.02, 0.98, cat, transform=ax.transAxes,
                fontsize=6.5, color=cat_col, va="top", alpha=0.8)
        ax.set_xlabel("")
        ax.set_ylabel("log₂(TPM+1)", fontsize=8)
        ax.set_xticks(range(len(xs_here)))
        ax.set_xticklabels(xs_here, rotation=40, ha="right", fontsize=7)
        ax.spines[["top", "right"]].set_visible(False)

    for ax in af[len(genes_present):]:
        ax.set_visible(False)

    # Legend
    cat_leg = [mpatches.Patch(color=HEATMAP_CATEGORY_COLOURS.get(c, "#888"), label=c)
               for c in HEATMAP_GENES_BY_CATEGORY]
    grp_leg = []
    for g in groups_present:
        for s in STAGE_ORDER:
            if (g, s) in xs_ordered:
                grp_leg.append(mpatches.Patch(
                    color=STAGE_AWARE_COLOURS.get((g, s), SC.GROUP_COLOURS.get(g, "#888")),
                    label=f"{g} {s}"))
    fig.legend(handles=cat_leg + grp_leg,
               title="category / group×stage",
               loc="lower right", ncol=3, frameon=False, fontsize=8)

    no_wt_note = " [no WT]" if exclude_wt else ""
    fig.suptitle(
        f"{title}{no_wt_note}  [stage-split x-axis]\n"
        f"Full heatmap gene set — coloured by category  |  X = group × stage"
        + ("  |  dashed = WT median" if not exclude_wt else ""),
        fontsize=9, y=1.02)
    plt.tight_layout()
    save(fig, filename)


# =============================================================================
#  NEW: Stage-faceted heatmap — columns sorted by stage within each group
# =============================================================================

def _draw_heatmap_by_stage(cols: list[str], title: str, filename: str,
                            pattern: list[tuple[str, int]],
                            scores: dict[str, float],
                            calls: dict[str, str],
                            show_context: bool = True):
    """
    Heatmap with columns sorted: WT | then each genotype × stage block.
    Within each genotype, stages run nc14b → nc14d → gastr.
    A two-row colour strip above the heatmap shows genotype (top) and stage (bottom).
    """
    if not cols:
        return

    high_genes = [g for g, d in pattern if d == +1 and g in z_mat.index]
    low_genes  = [g for g, d in pattern if d == -1 and g in z_mat.index]
    ctx_here   = [g for g in found_ctx if g in z_mat.index] if show_context else []
    row_order  = high_genes + low_genes + ctx_here
    n_pat      = len(high_genes) + len(low_genes)

    if not (high_genes + low_genes):
        return

    # Sort: WT first, then each group × stage in order
    wt_here = [c for c in cols if group_of.get(c) == "WT"]
    mut_here = [c for c in cols if group_of.get(c) != "WT"]

    col_order: list[str] = list(wt_here)
    stage_boundaries: list[int] = []  # x positions where stage changes
    group_boundaries: list[int] = [len(wt_here)]  # x positions where group changes

    for grp in SC.GROUP_ORDER:
        if grp == "WT":
            continue
        grp_cols = [c for c in mut_here if group_of.get(c) == grp]
        if not grp_cols:
            continue
        prev_stage = None
        for stage in STAGE_ORDER:
            stage_cols = sorted([c for c in grp_cols if SC.get_stage(c) == stage],
                                 key=lambda c: -scores.get(c, 0))
            if stage_cols:
                if prev_stage is not None:
                    stage_boundaries.append(len(col_order))
                col_order.extend(stage_cols)
                prev_stage = stage
        if len(col_order) > group_boundaries[-1]:
            group_boundaries.append(len(col_order))

    z_plot   = z_mat.loc[row_order, col_order]
    log_plot = log_mat.loc[row_order, col_order]
    col_lbls = [short_label(c) for c in col_order]
    n_genes  = len(row_order)
    n_emb    = len(col_order)
    pat_gene_set = {g for g, _ in pattern}

    fig_w = max(10, n_emb * 0.75 + 5)
    fig_h = max(5, n_genes * 0.55 + 3)
    fig, (ax_z, ax_log) = plt.subplots(
        1, 2, figsize=(fig_w * 2, fig_h),
        gridspec_kw={"wspace": 0.35})

    log_vmax = float(np.percentile(log_plot.values[log_plot.values > 0], 98)) \
               if (log_plot.values > 0).any() else 8

    for ax, data, cmap, cbar_lbl, vmin, vmax, center in [
        (ax_z,   z_plot,   "RdBu_r", "Z-score (QN-TPM)", -2.5, 2.5,  0),
        (ax_log, log_plot, "YlOrRd", "log₂(TPM+1)",       0, log_vmax, None),
    ]:
        kw = dict(cmap=cmap, xticklabels=col_lbls,
                  yticklabels=row_order, linewidths=0.4, linecolor="white",
                  cbar_kws={"label": cbar_lbl, "shrink": 0.5})
        if center is not None:
            kw.update(vmin=vmin, vmax=vmax, center=center)
        else:
            kw.update(vmin=vmin, vmax=vmax)
        sns.heatmap(data, ax=ax, **kw)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=7.5)

        for tick, gene in zip(ax.get_yticklabels(), row_order):
            tick.set_fontstyle("italic"); tick.set_fontsize(9)
            if gene in pat_gene_set:
                d = next((dd for gg, dd in pattern if gg == gene), 0)
                tick.set_color(DIRECTION_COLOURS.get(d, "black"))
            else:
                tick.set_color("#777777")

        # Direction strip (left of z-score panel)
        if ax is ax_z:
            ds = ax.inset_axes([-0.06, 0, 0.04, 1], transform=ax.transAxes)
            ds.set_xlim(0, 1); ds.set_ylim(0, n_genes); ds.axis("off")
            for i, gene in enumerate(row_order):
                col_d = (DIRECTION_COLOURS.get(
                    next((dd for gg, dd in pattern if gg == gene), 0), "#AAA")
                    if gene in pat_gene_set else "#DDDDDD")
                ds.add_patch(plt.Rectangle((0, n_genes - i - 1), 1, 1,
                                           color=col_d, transform=ds.transData))

        # ── Two-row annotation strip above heatmap ────────────────────────────
        # Row 1 (bottom): match call / WT colour
        strip1 = ax.inset_axes([0, 1.01, 1, 0.030], transform=ax.transAxes)
        strip1.set_xlim(0, n_emb); strip1.set_ylim(0, 1); strip1.axis("off")
        for j, col in enumerate(col_order):
            c = (MATCH_COLOURS["WT"] if group_of.get(col) == "WT"
                 else MATCH_COLOURS.get(calls.get(col, "other"), "#AAA"))
            strip1.add_patch(plt.Rectangle((j, 0), 1, 1, color=c,
                                           transform=strip1.transData))

        # Row 2: stage colour
        strip2 = ax.inset_axes([0, 1.04, 1, 0.020], transform=ax.transAxes)
        strip2.set_xlim(0, n_emb); strip2.set_ylim(0, 1); strip2.axis("off")
        for j, col in enumerate(col_order):
            stg = SC.get_stage(col)
            strip2.add_patch(plt.Rectangle((j, 0), 1, 1,
                                           color=STAGE_COLOURS.get(stg, "#CCC"),
                                           transform=strip2.transData))

        # Row 3: group colour
        strip3 = ax.inset_axes([0, 1.065, 1, 0.020], transform=ax.transAxes)
        strip3.set_xlim(0, n_emb); strip3.set_ylim(0, 1); strip3.axis("off")
        for j, col in enumerate(col_order):
            grp_c = embryo_colour(col)
            strip3.add_patch(plt.Rectangle((j, 0), 1, 1, color=grp_c,
                                           transform=strip3.transData))

        # Score annotations
        sc_ax = ax.inset_axes([0, 1.09, 1, 0.06], transform=ax.transAxes)
        sc_ax.set_xlim(0, n_emb); sc_ax.set_ylim(0, 1); sc_ax.axis("off")
        for j, col in enumerate(col_order):
            if group_of.get(col) != "WT":
                sv = scores.get(col, float("nan"))
                if not np.isnan(sv):
                    sc_ax.text(j + 0.5, 0.1, f"{sv:.2f}", ha="center",
                               va="bottom", fontsize=6,
                               color=MATCH_COLOURS.get(calls.get(col, "other"), "#888"))

        # Group & stage boundary lines
        for bnd in group_boundaries[:-1]:
            if 0 < bnd < n_emb:
                ax.axvline(bnd, color="black", lw=2.0)
        for bnd in stage_boundaries:
            if 0 < bnd < n_emb:
                ax.axvline(bnd, color="#888888", lw=1.0, ls="--")

        # HIGH/LOW boundary
        if 0 < len(high_genes) < n_pat:
            ax.axhline(len(high_genes), color="black", lw=1.5, ls="--")
        if ctx_here and n_pat > 0:
            ax.axhline(n_pat, color="#666", lw=1.5, ls=":")

    # Group × stage header text
    hdr = ax_z.inset_axes([0, 1.16, 1, 0.07], transform=ax_z.transAxes)
    hdr.set_xlim(0, n_emb); hdr.set_ylim(0, 1); hdr.axis("off")
    # WT label
    if wt_here:
        hdr.text(len(wt_here) / 2, 0.5, "WT", ha="center", va="center",
                 fontsize=9, fontweight="bold",
                 color=MATCH_COLOURS["WT"])
    # Genotype + stage labels
    pos = len(wt_here)
    for grp in SC.GROUP_ORDER:
        if grp == "WT":
            continue
        for stage in STAGE_ORDER:
            stage_cols_here = [c for c in col_order[pos:]
                                if group_of.get(c) == grp
                                and SC.get_stage(c) == stage]
            # find their actual positions
            actual = [i for i, c in enumerate(col_order)
                      if group_of.get(c) == grp and SC.get_stage(c) == stage]
            if actual:
                mid = np.mean(actual)
                hdr.text(mid + 0.5, 0.5,
                         f"{grp}\n{stage}",
                         ha="center", va="center", fontsize=6.5,
                         color=STAGE_AWARE_COLOURS.get((grp, "nc14b"), SC.GROUP_COLOURS.get(grp, "#888")))

    # Legend
    leg_handles = [
        mpatches.Patch(color=MATCH_COLOURS["match"],    label=f"match (≥{MATCH_THRESHOLD:.0%})"),
        mpatches.Patch(color=MATCH_COLOURS["partial"],  label=f"partial (≥{PARTIAL_THRESHOLD:.0%})"),
        mpatches.Patch(color=MATCH_COLOURS["mismatch"], label="mismatch"),
        mpatches.Patch(color=MATCH_COLOURS["WT"],       label="WT"),
        mpatches.Patch(color=DIRECTION_COLOURS[+1],     label="expected HIGH"),
        mpatches.Patch(color=DIRECTION_COLOURS[-1],     label="expected LOW"),
        mpatches.Patch(color=STAGE_COLOURS["nc14b"],    label="nc14b"),
        mpatches.Patch(color=STAGE_COLOURS["nc14d"],    label="nc14d"),
        mpatches.Patch(color=STAGE_COLOURS["gastr"],    label="gastr"),
    ]
    ax_z.legend(handles=leg_handles, bbox_to_anchor=(1.02, 1.0),
                loc="upper left", frameon=False, fontsize=7.5, title="legend")

    ctx_note = f"\nContext: {', '.join(ctx_here)}" if ctx_here else ""
    fig.suptitle(
        f"{title}  [stage-sorted columns]\n"
        f"Columns: WT | genotype×stage blocks (nc14b→nc14d→gastr)\n"
        f"Strips: match call | stage colour | group colour{ctx_note}",
        fontsize=9, y=1.02)
    save(fig, filename)


# =============================================================================
#  NEW: Stage-split DV violin
# =============================================================================

def _draw_dv_violins_by_stage():
    """DV axis violins with x-axis split by group × stage."""
    all_dv_symbols = (found_dv["ventral"] + found_dv["lateral"] +
                      found_dv["dorsal"])
    if not all_dv_symbols:
        return

    records = []
    for col in all_cols:
        grp   = group_of.get(col, "unknown")
        stage = SC.get_stage(col)
        for gene in all_dv_symbols:
            axis_cat = ("ventral" if gene in found_dv["ventral"] else
                        "lateral" if gene in found_dv["lateral"] else "dorsal")
            records.append({
                "embryo":   col,
                "group":    grp,
                "stage":    stage,
                "gs_label": _stage_group_label(grp, stage),
                "gene":     gene,
                "DV_axis":  axis_cat,
                "log2TPM":  float(log_mat.loc[gene, col]),
            })
    df = pd.DataFrame(records)

    groups_present = [g for g in SC.GROUP_ORDER if g in df["group"].unique()]

    xs_ordered: list[tuple[str, str]] = []
    for grp in groups_present:
        for stage in STAGE_ORDER:
            if not df[(df["group"] == grp) & (df["stage"] == stage)].empty:
                xs_ordered.append((grp, stage))

    xs_labels  = [_stage_group_label(g, s) for g, s in xs_ordered]
    xs_palette = {_stage_group_label(g, s):
                  STAGE_AWARE_COLOURS.get((g, s), SC.GROUP_COLOURS.get(g, "#888"))
                  for g, s in xs_ordered}
    xs_alpha   = {_stage_group_label(g, s): SC.STAGE_ALPHA.get(s, 0.7)
                  for g, s in xs_ordered}

    ncols_dv = min(5, len(all_dv_symbols))
    nrows_dv = int(np.ceil(len(all_dv_symbols) / ncols_dv))
    fig, axes = plt.subplots(nrows_dv, ncols_dv,
                             figsize=(max(6, len(xs_ordered) * 0.9 + 1) * ncols_dv,
                                      4.5 * nrows_dv),
                             sharey=False)
    af = np.array(axes).flatten() if len(all_dv_symbols) > 1 else [axes]

    for ax, gene in zip(af[:len(all_dv_symbols)], all_dv_symbols):
        sub      = df[df["gene"] == gene]
        axis_cat = sub["DV_axis"].iloc[0] if len(sub) else "ventral"
        dir_col  = DV_AXIS_COLOURS.get(axis_cat, "#888")
        bg_col   = {"ventral": "#FFF5F5", "lateral": "#FFFBF0",
                    "dorsal":  "#F0F5FF"}.get(axis_cat, "#F8F8F8")

        xs_here = [lbl for lbl in xs_labels if lbl in sub["gs_label"].values]
        sub_here = sub[sub["gs_label"].isin(xs_here)]
        xs_violin = [x for x in xs_here
                     if sub_here[sub_here["gs_label"] == x].shape[0] >= 2]

        # Draw violins using full xs_here as order so positions stay aligned
        # with the scatter points below. Groups with n<2 are skipped (no violin)
        # but their scatter points still appear at the correct x position.
        if xs_violin:
            sns.violinplot(data=sub_here[sub_here["gs_label"].isin(xs_violin)],
                           x="gs_label", y="log2TPM",
                           hue="gs_label", order=xs_here,
                           palette=xs_palette,
                           inner=None, cut=0, linewidth=0.7,
                           legend=False, ax=ax)

        for j, lbl in enumerate(xs_here):
            pts = sub_here[sub_here["gs_label"] == lbl]["log2TPM"].values
            if len(pts) == 0:
                continue
            jitter = np.random.uniform(-0.18, 0.18, len(pts))
            ax.scatter(np.full(len(pts), xs_here.index(lbl)) + jitter,
                       pts,
                       color=xs_palette.get(lbl, "#888"),
                       alpha=xs_alpha.get(lbl, 0.8),
                       s=45, zorder=4, edgecolors="white", linewidths=0.4)

        wt_pts = sub_here[sub_here["group"] == "WT"]["log2TPM"]
        if not wt_pts.empty:
            ax.axhline(wt_pts.median(),
                       color=MATCH_COLOURS["WT"],
                       lw=1.5, ls="--", alpha=0.6)

        # Group boundaries
        prev_grp = None
        for j, (g, s) in enumerate(xs_ordered):
            if prev_grp is not None and g != prev_grp:
                ax.axvline(j - 0.5, color="#AAAAAA", lw=1.0, ls=":")
            prev_grp = g

        ax.set_facecolor(bg_col)
        ax.set_title(f"{gene}\n[{axis_cat}]", fontsize=9.5,
                     fontstyle="italic", fontweight="bold", color=dir_col)
        ax.set_xlabel("")
        ax.set_ylabel("log₂(TPM+1)", fontsize=8)
        ax.set_xticks(range(len(xs_here)))
        ax.set_xticklabels(xs_here, fontsize=7.5)
        ax.spines[["top", "right"]].set_visible(False)

    for ax in af[len(all_dv_symbols):]:
        ax.set_visible(False)

    # Stage-aware legend: show nc14b and late-stage colour per group
    grp_leg = []
    for g in groups_present:
        c_early = STAGE_AWARE_COLOURS.get((g,"nc14b"), SC.GROUP_COLOURS.get(g,"#888"))
        c_late  = STAGE_AWARE_COLOURS.get((g,"nc14d"), c_early)
        grp_leg.append(mpatches.Patch(color=c_early, label=f"{g} nc14b"))
        if c_late != c_early:
            grp_leg.append(mpatches.Patch(color=c_late, label=f"{g} nc14d/gastr"))
    axis_leg = [mpatches.Patch(color=DV_AXIS_COLOURS[a], label=a)
                for a in ("ventral", "lateral", "dorsal")]
    fig.legend(handles=grp_leg + axis_leg, title="group / DV axis",
               loc="lower right", ncol=3, frameon=False, fontsize=8)

    fig.suptitle(
        "D/V axis — stage-split  |  X = group × stage\n"
        "Red=ventral (dl sna twi NetA)  |  Orange=lateral (sog brk pyr elav)  |  "
        "Blue=dorsal (Doc1 anc Abd-B)\n"
        "Dashed = WT median  |  lighter alpha = later stage",
        fontsize=9, y=1.02)
    plt.tight_layout()
    save(fig, "violin_DV_axis_by_stage")

def _draw_score_bars(target_cols, wt_local, title, filename, scores, calls, pat_name):
    cols_plot = (sorted(wt_local,     key=lambda c: -scores.get(c, 0)) +
                 sorted(target_cols,  key=lambda c: -scores.get(c, 0)))
    if not cols_plot: return

    fig, ax = plt.subplots(figsize=(max(8, len(cols_plot) * 0.6), 5.5))
    for i, col in enumerate(cols_plot):
        grp   = group_of.get(col, "")
        score = scores.get(col, float("nan"))
        call  = calls.get(col, "other")
        colour = embryo_colour(col) if grp != "WT" else MATCH_COLOURS["WT"]
        # Desaturate mismatches slightly for readability
        if call == "mismatch" and grp != "WT":
            import colorsys as _cs
            _h,_s,_v = _cs.rgb_to_hsv(*[int(colour[i:i+2],16)/255 for i in (1,3,5)])
            colour = "#{:02x}{:02x}{:02x}".format(
                *[int(c*255) for c in _cs.hsv_to_rgb(_h, _s*0.45, min(_v*1.3, 1.0))])

        ax.bar(i, score if not np.isnan(score) else 0,
               color=colour, edgecolor="white", width=0.75, zorder=3)
        if not np.isnan(score):
            ax.text(i, score + 0.01, f"{score:.2f}",
                    ha="center", va="bottom", fontsize=7)
        # Stage marker as x-tick colour hint
        stg = SC.get_stage(col)
        ax.get_xticklabels()  # ensure ticks exist

    ax.axhline(MATCH_THRESHOLD,   color="#2ecc71", lw=1.5, ls="--",
               label=f"match ({MATCH_THRESHOLD:.0%})")
    ax.axhline(PARTIAL_THRESHOLD, color="#f39c12", lw=1.2, ls=":",
               label=f"partial ({PARTIAL_THRESHOLD:.0%})")
    if wt_local:
        ax.axvline(len(wt_local) - 0.5, color="black", lw=2.0)

    ax.set_xticks(range(len(cols_plot)))
    xlbls = []
    for col in cols_plot:
        bare  = col.split("__")[-1]
        stage = SC.get_stage(col)
        stage_sym = {"nc14b": "", "nc14d": "▪", "gastr": "▲"}.get(stage, "")
        import re
        m = re.search(r"(bc\d+)$", bare)
        bc = m.group(1) if m else bare.split("_")[-1]
        xlbls.append(f"{group_of.get(col,'')}\n{bc}{stage_sym}")
    ax.set_xticklabels(xlbls, rotation=40, ha="right", fontsize=8)
    ax.set_ylabel(f"{pat_name} pattern score", fontsize=10)
    ax.set_ylim(0, 1.1)
    ax.set_title(title, fontsize=11)
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)

    n_match = sum(1 for c in target_cols if calls[c] == "match")
    n_part  = sum(1 for c in target_cols if calls[c] == "partial")
    n_mis   = sum(1 for c in target_cols if calls[c] == "mismatch")
    ax.set_xlabel(f"n={len(target_cols)}  match={n_match}  "
                  f"partial={n_part}  mismatch={n_mis}", fontsize=9)
    plt.tight_layout()
    save(fig, filename)


# =============================================================================
#  Generate all figures
# =============================================================================

print("\n[posclass] ── HLH54F analysis ──")
_draw_score_bars(hlh_cols, wt_cols,
    title=f"HLH54F_BOTCv pattern score  [{len(hlh_cols)} embryos, {len(wt_cols)} WT]",
    filename="pattern_scores_HLH54F",
    scores=hlh_scores, calls=hlh_calls, pat_name="HLH54F_BOTCv")

_draw_heatmap(wt_cols + hlh_cols,
    title="HLH54F_BOTCv vs WT — positional pattern genes",
    filename="heatmap_HLH54F_vs_WT",
    pattern=HLH54F_PATTERN, scores=hlh_scores, calls=hlh_calls)

_draw_heatmap_by_stage(wt_cols + hlh_cols,
    title="HLH54F_BOTCv vs WT — stage-sorted columns",
    filename="heatmap_HLH54F_vs_WT_by_stage",
    pattern=HLH54F_PATTERN, scores=hlh_scores, calls=hlh_calls)

_draw_violins(wt_cols + hlh_cols,
    title=f"HLH54F_BOTCv vs WT  [{len(hlh_cols)} embryos]",
    filename="violin_HLH54F_vs_WT",
    pattern=HLH54F_PATTERN, show_pairrule=True)

_draw_violins_by_stage(wt_cols + hlh_cols,
    title=f"HLH54F_BOTCv vs WT — stage-split",
    filename="violin_HLH54F_vs_WT_by_stage",
    pattern=HLH54F_PATTERN, show_pairrule=True)

_draw_violins_heatmap_genes(wt_cols + hlh_cols,
    title=f"HLH54F_BOTCv vs WT — heatmap gene set  [{len(hlh_cols)} embryos]",
    filename="violin_HLH54F_vs_WT_heatmap_genes")

_draw_violins_heatmap_genes_by_stage(wt_cols + hlh_cols,
    title=f"HLH54F_BOTCv vs WT — heatmap gene set  [{len(hlh_cols)} embryos]",
    filename="violin_HLH54F_vs_WT_heatmap_genes_by_stage")

_draw_violins_heatmap_genes(wt_cols + hlh_cols,
    title=f"HLH54F_BOTCv only — heatmap gene set  [{len(hlh_cols)} embryos]",
    filename="violin_HLH54F_no_WT_heatmap_genes", exclude_wt=True)

_draw_violins_heatmap_genes_by_stage(wt_cols + hlh_cols,
    title=f"HLH54F_BOTCv only — heatmap gene set  [{len(hlh_cols)} embryos]",
    filename="violin_HLH54F_no_WT_heatmap_genes_by_stage", exclude_wt=True)


print("\n[posclass] ── FoxL1 analysis ──")
_draw_score_bars(foxl1_cols, wt_cols,
    title=f"FoxL1_BOTv pattern score  [{len(foxl1_cols)} embryos, {len(wt_cols)} WT]",
    filename="pattern_scores_foxl1",
    scores=foxl1_scores, calls=foxl1_calls, pat_name="FoxL1_BOTv")

_draw_heatmap(wt_cols + foxl1_cols,
    title="FoxL1_BOTv vs WT — positional pattern genes",
    filename="heatmap_foxl1_vs_WT",
    pattern=FOXL1_PATTERN, scores=foxl1_scores, calls=foxl1_calls)

_draw_heatmap_by_stage(wt_cols + foxl1_cols,
    title="FoxL1_BOTv vs WT — stage-sorted columns",
    filename="heatmap_foxl1_vs_WT_by_stage",
    pattern=FOXL1_PATTERN, scores=foxl1_scores, calls=foxl1_calls)

_draw_violins(wt_cols + foxl1_cols,
    title=f"FoxL1_BOTv vs WT  [{len(foxl1_cols)} embryos]",
    filename="violin_foxl1_vs_WT",
    pattern=FOXL1_PATTERN, show_pairrule=True)

_draw_violins_by_stage(wt_cols + foxl1_cols,
    title=f"FoxL1_BOTv vs WT — stage-split",
    filename="violin_foxl1_vs_WT_by_stage",
    pattern=FOXL1_PATTERN, show_pairrule=True)

_draw_violins_heatmap_genes(wt_cols + foxl1_cols,
    title=f"FoxL1_BOTv vs WT — heatmap gene set  [{len(foxl1_cols)} embryos]",
    filename="violin_foxl1_vs_WT_heatmap_genes")

_draw_violins_heatmap_genes_by_stage(wt_cols + foxl1_cols,
    title=f"FoxL1_BOTv vs WT — heatmap gene set  [{len(foxl1_cols)} embryos]",
    filename="violin_foxl1_vs_WT_heatmap_genes_by_stage")

_draw_violins_heatmap_genes(wt_cols + foxl1_cols,
    title=f"FoxL1_BOTv only — heatmap gene set  [{len(foxl1_cols)} embryos]",
    filename="violin_foxl1_no_WT_heatmap_genes", exclude_wt=True)

_draw_violins_heatmap_genes_by_stage(wt_cols + foxl1_cols,
    title=f"FoxL1_BOTv only — heatmap gene set  [{len(foxl1_cols)} embryos]",
    filename="violin_foxl1_no_WT_heatmap_genes_by_stage", exclude_wt=True)


print("\n[posclass] ── DV axis violins ──")
_draw_dv_violins()
_draw_dv_violins_by_stage()


if not args.no_all_groups:
    print("\n[posclass] ── All groups combined ──")
    combined_pat = list(dict.fromkeys(HLH54F_PATTERN + FOXL1_PATTERN))
    _draw_heatmap(all_cols,
        title="All groups — positional pattern genes",
        filename="heatmap_all_groups",
        pattern=combined_pat, scores=hlh_scores, calls=hlh_calls)
    _draw_heatmap_by_stage(all_cols,
        title="All groups — stage-sorted columns",
        filename="heatmap_all_groups_by_stage",
        pattern=combined_pat, scores=hlh_scores, calls=hlh_calls)
    _draw_violins(all_cols,
        title="All groups — positional gene expression",
        filename="violin_all_groups",
        pattern=combined_pat, show_pairrule=True)
    _draw_violins(all_cols,
        title="All groups (no WT) — positional gene expression",
        filename="violin_all_groups_no_WT",
        pattern=combined_pat, show_pairrule=True, exclude_wt=True)
    _draw_violins_by_stage(all_cols,
        title="All groups — stage-split",
        filename="violin_all_groups_by_stage",
        pattern=combined_pat, show_pairrule=True)
    _draw_violins_by_stage(all_cols,
        title="All groups (no WT) — stage-split",
        filename="violin_all_groups_no_WT_by_stage",
        pattern=combined_pat, show_pairrule=True, exclude_wt=True)
    _draw_violins_heatmap_genes(all_cols,
        title="All groups — heatmap gene set",
        filename="violin_all_groups_heatmap_genes")
    _draw_violins_heatmap_genes_by_stage(all_cols,
        title="All groups — heatmap gene set",
        filename="violin_all_groups_heatmap_genes_by_stage")
    _draw_violins_heatmap_genes(all_cols,
        title="All groups (no WT) — heatmap gene set",
        filename="violin_all_groups_no_WT_heatmap_genes", exclude_wt=True)
    _draw_violins_heatmap_genes_by_stage(all_cols,
        title="All groups (no WT) — heatmap gene set",
        filename="violin_all_groups_no_WT_heatmap_genes_by_stage", exclude_wt=True)


# ── Ventralized (FoxL1_BOTv + HLH54F_BOTCv) vs WT — combined ────────────────
print("\n[posclass] ── Ventralized combined vs WT ──")

vent_cols   = foxl1_cols + hlh_cols   # both ventralized genotypes
vent_wt     = wt_cols + vent_cols
combined_vent_pat = list(dict.fromkeys(HLH54F_PATTERN + FOXL1_PATTERN))

# Use HLH54F scores for the combined figure (scores span all embryos)
# — embryos outside each pattern's target group will have low scores,
# which is informative as a contrast
_vent_scores = {**hlh_scores, **foxl1_scores}  # per-embryo best score
_vent_calls  = {}
for col in vent_wt:
    grp = group_of.get(col, "")
    if grp == "WT":
        _vent_calls[col] = "WT"
    elif grp == "HLH54F_BOTCv":
        _vent_calls[col] = hlh_calls[col]
    elif grp == "FoxL1_BOTv":
        _vent_calls[col] = foxl1_calls[col]
    else:
        _vent_calls[col] = hlh_calls.get(col, "other")

_draw_heatmap(vent_wt,
    title="Ventralized (FoxL1_BOTv + HLH54F_BOTCv) vs WT",
    filename="heatmap_ventralized_combined_vs_WT",
    pattern=combined_vent_pat, scores=_vent_scores, calls=_vent_calls)

_draw_heatmap_by_stage(vent_wt,
    title="Ventralized vs WT — stage-sorted columns",
    filename="heatmap_ventralized_combined_vs_WT_by_stage",
    pattern=combined_vent_pat, scores=_vent_scores, calls=_vent_calls)

_draw_violins(vent_wt,
    title=f"Ventralized combined vs WT  [Fox:{len(foxl1_cols)} HLH:{len(hlh_cols)} WT:{len(wt_cols)}]",
    filename="violin_ventralized_combined_vs_WT",
    pattern=combined_vent_pat, show_pairrule=True)

_draw_violins(vent_wt,
    title=f"Ventralized combined (no WT)  [Fox:{len(foxl1_cols)} HLH:{len(hlh_cols)}]",
    filename="violin_ventralized_combined_no_WT",
    pattern=combined_vent_pat, show_pairrule=True, exclude_wt=True)

_draw_violins_by_stage(vent_wt,
    title="Ventralized combined vs WT — stage-split",
    filename="violin_ventralized_combined_vs_WT_by_stage",
    pattern=combined_vent_pat, show_pairrule=True)

_draw_violins_by_stage(vent_wt,
    title="Ventralized combined (no WT) — stage-split",
    filename="violin_ventralized_combined_no_WT_by_stage",
    pattern=combined_vent_pat, show_pairrule=True, exclude_wt=True)

_draw_violins_heatmap_genes(vent_wt,
    title=f"Ventralized combined vs WT — heatmap gene set  [Fox:{len(foxl1_cols)} HLH:{len(hlh_cols)} WT:{len(wt_cols)}]",
    filename="violin_ventralized_combined_vs_WT_heatmap_genes")

_draw_violins_heatmap_genes_by_stage(vent_wt,
    title=f"Ventralized combined vs WT — heatmap gene set  [Fox:{len(foxl1_cols)} HLH:{len(hlh_cols)} WT:{len(wt_cols)}]",
    filename="violin_ventralized_combined_vs_WT_heatmap_genes_by_stage")

_draw_violins_heatmap_genes(vent_wt,
    title=f"Ventralized combined (no WT) — heatmap gene set  [Fox:{len(foxl1_cols)} HLH:{len(hlh_cols)}]",
    filename="violin_ventralized_combined_no_WT_heatmap_genes", exclude_wt=True)

_draw_violins_heatmap_genes_by_stage(vent_wt,
    title=f"Ventralized combined (no WT) — heatmap gene set  [Fox:{len(foxl1_cols)} HLH:{len(hlh_cols)}]",
    filename="violin_ventralized_combined_no_WT_heatmap_genes_by_stage", exclude_wt=True)

# Score bars — show both HLH54F and FoxL1 scores side by side per embryo
_fig_sb, (_ax_hlh, _ax_fox) = plt.subplots(2, 1,
    figsize=(max(8, len(vent_wt) * 0.65), 9), sharex=True)

for _ax, _scores, _calls, _pat_name, _target in [
    (_ax_hlh, hlh_scores,   hlh_calls,   "HLH54F_BOTCv", hlh_cols),
    (_ax_fox, foxl1_scores, foxl1_calls, "FoxL1_BOTv",   foxl1_cols),
]:
    _cols_plot = wt_cols + sorted(vent_cols, key=lambda c: -_scores.get(c, 0))
    for _i, _col in enumerate(_cols_plot):
        _grp   = group_of.get(_col, "")
        _score = _scores.get(_col, float("nan"))
        _call  = _calls.get(_col, "other")
        _colour = (MATCH_COLOURS["WT"] if _grp == "WT"
                   else embryo_colour(_col))
        if _call == "mismatch":
            import colorsys as _cs2
            _h2,_s2,_v2 = _cs2.rgb_to_hsv(
                *[int(_colour[j:j+2],16)/255 for j in (1,3,5)])
            _colour = "#{:02x}{:02x}{:02x}".format(
                *[int(c*255) for c in _cs2.hsv_to_rgb(_h2, _s2*0.4, min(_v2*1.3,1))])
        _ax.bar(_i, _score if not np.isnan(_score) else 0,
                color=_colour, edgecolor="white", width=0.75, zorder=3)
        if not np.isnan(_score):
            _ax.text(_i, _score + 0.01, f"{_score:.2f}",
                     ha="center", va="bottom", fontsize=6)
    _ax.axhline(MATCH_THRESHOLD,   color="#2ecc71", lw=1.5, ls="--",
                label=f"match ({MATCH_THRESHOLD:.0%})")
    _ax.axhline(PARTIAL_THRESHOLD, color="#f39c12", lw=1.2, ls=":",
                label=f"partial ({PARTIAL_THRESHOLD:.0%})")
    if wt_cols:
        _ax.axvline(len(wt_cols) - 0.5, color="black", lw=2.0)
    # boundary between foxl1 and hlh54f
    _n_fox = sum(1 for c in vent_cols if group_of.get(c) == "FoxL1_BOTv")
    if _n_fox > 0 and _n_fox < len(vent_cols):
        _ax.axvline(len(wt_cols) + _n_fox - 0.5, color="#888", lw=1.2, ls="--")
    _ax.set_ylabel(f"{_pat_name}\npattern score", fontsize=9)
    _ax.set_ylim(0, 1.15)
    _ax.legend(frameon=False, fontsize=8, loc="upper right")
    _ax.spines[["top","right"]].set_visible(False)
    _n_m = sum(1 for c in _target if _calls[c] == "match")
    _n_p = sum(1 for c in _target if _calls[c] == "partial")
    _n_x = sum(1 for c in _target if _calls[c] == "mismatch")
    _ax.set_title(f"{_pat_name} score  |  match={_n_m} partial={_n_p} mismatch={_n_x}",
                  fontsize=9)

import re as _re
_xlbls = []
for _col in (wt_cols + sorted(vent_cols, key=lambda c: -hlh_scores.get(c, 0))):
    _bare  = _col.split("__")[-1]
    _stage = SC.get_stage(_col)
    _ssym  = {"nc14b":"","nc14d":"▪","gastr":"▲"}.get(_stage,"")
    _m2    = _re.search(r"(bc\d+)$", _bare)
    _bc    = _m2.group(1) if _m2 else _bare.split("_")[-1]
    _xlbls.append(f"{group_of.get(_col,'')[:4]}\n{_bc}{_ssym}")

_ax_fox.set_xticks(range(len(_xlbls)))
_ax_fox.set_xticklabels(_xlbls, rotation=40, ha="right", fontsize=7.5)
_fig_sb.suptitle(
    "Ventralized combined vs WT — both pattern scores\n"
    "Top: HLH54F_BOTCv score  |  Bottom: FoxL1_BOTv score\n"
    "Stage: ▪=nc14d  ▲=gastr  |  Dashed=WT boundary  |  Grey=FoxL1/HLH boundary",
    fontsize=9)
plt.tight_layout()
save(_fig_sb, "pattern_scores_ventralized_combined")

# ── Per-call heatmaps and violins ─────────────────────────────────────────────
for pat_name, target, pattern, scores, calls in [
    ("HLH54F_BOTCv", hlh_cols,   HLH54F_PATTERN, hlh_scores,   hlh_calls),
    ("FoxL1_BOTv",   foxl1_cols, FOXL1_PATTERN,  foxl1_scores, foxl1_calls),
]:
    print(f"\n[posclass] ── {pat_name} by call ──")

    call_order = {
        "WT":                     wt_cols,
        f"{pat_name}_match":      sorted([c for c in target if calls[c] == "match"],
                                         key=lambda c: -scores.get(c, 0)),
        f"{pat_name}_partial":    sorted([c for c in target if calls[c] == "partial"],
                                         key=lambda c: -scores.get(c, 0)),
        f"{pat_name}_mismatch":   sorted([c for c in target if calls[c] == "mismatch"],
                                         key=lambda c: -scores.get(c, 0)),
    }
    ordered_call_cols: list[str] = []
    for tier_cols in call_order.values():
        ordered_call_cols.extend(tier_cols)

    n_m = len(call_order[f"{pat_name}_match"])
    n_p = len(call_order[f"{pat_name}_partial"])
    n_x = len(call_order[f"{pat_name}_mismatch"])
    print(f"[posclass]   {pat_name}: {n_m} match, {n_p} partial, {n_x} mismatch")

    if ordered_call_cols:
        _draw_heatmap(ordered_call_cols,
            title=f"{pat_name} — embryos grouped by pattern call",
            filename=f"heatmap_{pat_name}_by_call",
            pattern=pattern, scores=scores, calls=calls)
        _draw_heatmap_by_stage(ordered_call_cols,
            title=f"{pat_name} — by call, stage-sorted",
            filename=f"heatmap_{pat_name}_by_call_by_stage",
            pattern=pattern, scores=scores, calls=calls)
        _draw_violins(ordered_call_cols,
            title=f"{pat_name} — expression by call",
            filename=f"violin_{pat_name}_by_call",
            pattern=pattern, show_pairrule=True)
        _draw_violins_by_stage(ordered_call_cols,
            title=f"{pat_name} — by call, stage-split",
            filename=f"violin_{pat_name}_by_call_by_stage",
            pattern=pattern, show_pairrule=True)


# =============================================================================
#  Write output tables
# =============================================================================

def _write_matched(target_cols, scores, calls, pattern, filename, pat_name):
    rows = []
    for col in target_cols:
        rows.append({
            "embryo":  col,
            "label":   short_label(col),
            "group":   group_of.get(col, ""),
            "stage":   SC.get_stage(col),
            "score":   round(scores.get(col, float("nan")), 3),
            "call":    calls.get(col, ""),
            **{f"z_{g}": round(float(z_mat.loc[g, col]), 3)
               for g, _ in pattern if g in z_mat.index},
        })
    df = pd.DataFrame(rows).sort_values("score", ascending=False)
    df.to_csv(os.path.join(OUT_DIR, filename), sep="\t", index=False)
    matched = (df["call"] == "match").sum()
    partial = (df["call"] == "partial").sum()
    print(f"[posclass]   {pat_name}: {matched} match, {partial} partial, "
          f"{len(df)-matched-partial} mismatch")
    return df


_write_matched(hlh_cols,   hlh_scores,   hlh_calls,   HLH54F_PATTERN,
               "matched_embryos_HLH54F_BOTCv.tsv", "HLH54F_BOTCv")
_write_matched(foxl1_cols, foxl1_scores, foxl1_calls, FOXL1_PATTERN,
               "matched_embryos_FoxL1_BOTv.tsv",   "FoxL1_BOTv")

all_rows = []
for col in all_cols:
    all_rows.append({
        "embryo":       col,
        "label":        short_label(col),
        "group":        group_of.get(col, ""),
        "stage":        SC.get_stage(col),
        "HLH54F_score": round(hlh_scores.get(col, float("nan")), 3),
        "HLH54F_call":  hlh_calls.get(col, ""),
        "foxl1_score":  round(foxl1_scores.get(col, float("nan")), 3),
        "foxl1_call":   foxl1_calls.get(col, ""),
    })
pd.DataFrame(all_rows).to_csv(
    os.path.join(OUT_DIR, "positional_classification.tsv"), sep="\t", index=False)

print(f"\n[posclass] Done.  All outputs in {OUT_DIR}/")
print("  pattern_scores_HLH54F.pdf/png")
print("  pattern_scores_foxl1.pdf/png")
print("  heatmap_HLH54F_vs_WT.pdf/png")
print("  heatmap_HLH54F_vs_WT_by_stage.pdf/png       ← NEW stage-sorted columns")
print("  heatmap_foxl1_vs_WT.pdf/png")
print("  heatmap_foxl1_vs_WT_by_stage.pdf/png         ← NEW stage-sorted columns")
print("  heatmap_HLH54F_by_call.pdf/png")
print("  heatmap_HLH54F_by_call_by_stage.pdf/png      ← NEW")
print("  heatmap_foxl1_by_call.pdf/png")
print("  heatmap_foxl1_by_call_by_stage.pdf/png       ← NEW")
print("  heatmap_all_groups.pdf/png")
print("  heatmap_all_groups_by_stage.pdf/png           ← NEW")
print("  violin_HLH54F_vs_WT.pdf/png")
print("  violin_HLH54F_vs_WT_by_stage.pdf/png         ← NEW stage-split x-axis")
print("  violin_foxl1_vs_WT.pdf/png")
print("  violin_foxl1_vs_WT_by_stage.pdf/png           ← NEW stage-split x-axis")
print("  violin_HLH54F_by_call.pdf/png")
print("  violin_HLH54F_by_call_by_stage.pdf/png       ← NEW")
print("  violin_foxl1_by_call.pdf/png")
print("  violin_foxl1_by_call_by_stage.pdf/png         ← NEW")
print("  violin_all_groups.pdf/png")
print("  violin_all_groups_by_stage.pdf/png             ← NEW")
print("  violin_DV_axis.pdf/png")
print("  violin_DV_axis_by_stage.pdf/png               ← NEW stage-split DV")
print("  matched_embryos_HLH54F_BOTCv.tsv")
print("  matched_embryos_FoxL1_BOTv.tsv")
print("  positional_classification.tsv")
print(f"\n  Match threshold : {MATCH_THRESHOLD:.0%}   Partial : {PARTIAL_THRESHOLD:.0%}")
print(f"  Stage-aware z-scoring active for: {sorted(STAGE_AWARE_GROUPS)}")
