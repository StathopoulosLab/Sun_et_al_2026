#!/usr/bin/env python3
"""
plot_coverage_tracks.py
───────────────────────────────────────────────────────────────────────────────
Generate genome-browser-style ATAC-seq coverage track plots from CPM-
normalized BigWig files, with optional RNA-seq violin panels and gene-body
annotation.

Modes (--mode):
  all              All genotypes, nc14b/nc14late anchors only  [default]
  ventralized      BOTv nc14b, BOTv nc14late, BOTCv nc14b, BOTCv nc14late
  non_ventralized  BOT, BOT_hR, BOTR, BOTC, BOTC_oR  (nc14b only)
  all_timepoints            All genotypes across EVERY confirmed
                             timepoint (nc14b, nc14late, nc14d, gastr,
                             st10_11, late) -- see sample_manifest.sh.
  ventralized_timepoints    Ventralized genotypes, all timepoints
  non_ventralized_timepoints Non-ventralized genotypes, all timepoints
  CVM              BOTv + BOTCv across nc14b, nc14late, and st10_11
  CVM_solo         BOTCv only across nc14b, nc14late, and st10_11

Use --bg_sub with any CVM mode (or any other mode) to swap all BigWig
stems from raw CPM (_CPM / _mean_CPM) to background-subtracted CPM
(_bgSubCPM / _mean_bgSubCPM) at runtime. Without --bg_sub the same
raw-CPM files used by all other modes are loaded.

Each genotype group shows:
  • Individual replicate tracks (light fill, thin outline) stacked first
  • Merged (mean) track below (bold fill) with genotype label

DEPENDENCIES:
  pip install pyBigWig matplotlib numpy pandas

USAGE EXAMPLES:

  # All genotypes, merged + replicates, with gene body
  python plot_coverage_tracks.py \\
      --gene sna \\
      --chrom chr2L --start 15_450_000 --end 15_460_000 \\
      --bigwig_dir ./BigWig_Files_CPM \\
      --annotation_bed dmel_r6_genes.bed \\
      --out sna_all.pdf

  # Ventralized only
  python plot_coverage_tracks.py \\
      --gene sna --chrom chr2L --start 15_450_000 --end 15_460_000 \\
      --bigwig_dir ./BigWig_Files_CPM \\
      --mode ventralized --out sna_vent.pdf

  # Non-ventralized only
  python plot_coverage_tracks.py \\
      --gene sna --chrom chr2L --start 15_450_000 --end 15_460_000 \\
      --bigwig_dir ./BigWig_Files_CPM \\
      --mode non_ventralized --out sna_nonvent.pdf

  # Batch: generates all three mode PDFs per gene into subfolders
  python plot_coverage_tracks.py \\
      --gene_bed key_genes.bed \\
      --bigwig_dir ./BigWig_Files_CPM \\
      --annotation_bed dmel_r6_genes.bed \\
      --out_dir ./coverage_plots

See --help for all options.
───────────────────────────────────────────────────────────────────────────────
"""

import argparse
import sys
import os
import math
from pathlib import Path
from typing import Optional, List, Dict, Tuple

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import matplotlib.ticker as mticker

try:
    import pyBigWig
except ImportError:
    sys.exit("ERROR: pyBigWig not installed.  Run: pip install pyBigWig")


# ─────────────────────────────────────────────────────────────────────────────
# GENOTYPE CONFIG
#
# MERGED_CONFIG  — one entry per genotype×timepoint
#   key → (display_label, hex_color, merged_bw_stem, timepoint_group, ventralized)
#
# REPLICATE_CONFIG — per-sample BigWigs
#   key → list of (rep_label, bw_stem)   (same key as MERGED_CONFIG)
# ─────────────────────────────────────────────────────────────────────────────

MERGED_CONFIG: Dict[str, tuple] = {
    # key          label                    color      bw_stem                      tp_grp      vent
    "BOTv"      : ("BOTv nc14b",           "#2ca02c", "BOTv_nc14b_mean_CPM",       "nc14b",    True ),
    "BOTv_late" : ("BOTv nc14late",            "#1a6b1a", "BOTv_nc14late_mean_CPM",    "nc14late", True ),
    "BOTCv"     : ("BOTCv nc14b",          "#e377c2", "BOTCv_nc14b_mean_CPM",      "nc14b",    True ),
    "BOTCv_late": ("BOTCv nc14late",           "#b5369a", "BOTCv_nc14late_mean_CPM",   "nc14late", True ),
    # ^ label intentionally omits "(aka nc14d)" -- HLH54F_3/4_nc14d_IR are
    # confirmed to be the same source BAMs as this pair (not a distinct
    # timepoint), so they were removed rather than merged/aliased in. This
    # is the single, only late timepoint for BOTCv.
    "BOT"       : ("BOT nc14b",            "#1f77b4", "BOT_nc14b_mean_CPM",        "nc14b",    False),
    "BOT_hR"    : ("BOT_hR nc14b",         "#17becf", "BOT_hR_nc14b_mean_CPM",     "nc14b",    False),
    "BOTR"      : ("BOTR nc14b",           "#9467bd", "BOTR_nc14b_mean_CPM",       "nc14b",    False),
    "BOTC"      : ("BOTC nc14b",           "#ff7f0e", "BOTC_nc14b_mean_CPM",       "nc14b",    False),
    "BOTC_oR"   : ("BOTC_oR nc14b",        "#bcbd22", "BOTC_oR_nc14b_mean_CPM",    "nc14b",    False),
    "yw_WT"     : ("yw_WT nc14b",          "#555555", "yw_WT_nc14b_mean_CPM",      "nc14b",    False),
    # ── CVM (super-late stage 10/11) ────────────────────────────────────────
    # Stems are raw-CPM base; --bg_sub swaps _CPM → _bgSubCPM at runtime
    "BOTv_st10" : ("BOTv st10",            "#2ca02c", "BOTv_st10_mean_CPM",        "st10",     True ),
    # NOTE: renamed BOTCv_st10 -> BOTCv_st10_11 and fixed the bw_stem below.
    # The old stem "BOTCv_st10_mean_CPM" pointed at a bigwig that
    # generate_diagnostic_bigwigs_{vent_v3,all_v2}.sh never actually
    # produces -- their timepoint key for these samples is "st10_11", not
    # "st10", so the real merged file is "BOTCv_st10_11_mean_*CPM.bigWig".
    "BOTCv_st10_11": ("BOTCv st10/11",     "#7a1f6b", "BOTCv_st10_11_mean_CPM", "st10_11",  True ),

    # ── Additional confirmed timepoints (from sample_manifest.sh) ──────────
    "BOTv_gastr"     : ("BOTv gastr",             "#0d3d0d", "BOTv_gastr_mean_CPM",     "gastr",   True ),
    "BOT_nc14d"      : ("BOT nc14late",               "#145380", "BOT_nc14d_mean_CPM",      "nc14d",   False),
    "BOT_hR_gastr"   : ("BOT_hR gastr",           "#0e8a8a", "BOT_hR_gastr_mean_CPM",   "gastr",   False),
    "BOTR_gastr"     : ("BOTR gastr",             "#5c3585", "BOTR_gastr_CPM",          "gastr",   False),

    # ── Additional confirmed timepoints (from sample_manifest.sh /
    #    sampleinfo.rtf). Color for BOTCv_gastr is my own pick (not
    #    in your genotype_colors palette) -- feel free to adjust. BOTC_nc14late
    #    uses your confirmed genotype_colors value (#d45e00). BOTCv nc14d
    #    was removed entirely: HLH54F_3/4_nc14d_IR are confirmed to be the
    #    same source BAMs as BOTCv_nc14late above (an alternate naming of
    #    identical sequencing data), not a distinct timepoint -- including it
    #    would have double-counted those reads.
    "BOTCv_gastr"    : ("BOTCv gastr",       "#8f2f7a", "BOTCv_gastr_CPM",       "gastr", True ),
    "BOTC_nc14late"  : ("BOTC nc14late",     "#d45e00", "BOTC_nc14late_mean_CPM", "nc14late", False),
    "BOTC_oR_late"   : ("BOTC_oR nc14late",      "#8c8c00", "BOTC_oR_late_mean_CPM", "late", False),
}

# Per-replicate BigWig stems.  Keys match MERGED_CONFIG.
# Edit rep stems to match your actual filenames in BigWig_Files_CPM/.
REPLICATE_CONFIG: Dict[str, List[Tuple[str, str]]] = {
    "BOTv"      : [("rep1", "BOTv_nc14b_rep1_CPM"),
                   ("rep2", "BOTv_nc14b_rep2_CPM")],
    "BOTv_late" : [("rep1", "BOTv_nc14late_rep1_CPM"),
                   ("rep2", "BOTv_nc14late_rep2_CPM")],
    "BOTCv"     : [("rep1", "BOTCv_nc14b_rep1_CPM"),
                   ("rep2", "BOTCv_nc14b_rep2_CPM")],
    "BOTCv_late": [("rep1", "BOTCv_nc14late_rep1_CPM"),
                   ("rep2", "BOTCv_nc14late_rep2_CPM")],
    "BOT"       : [("rep1", "BOT_nc14b_rep1_CPM"),
                   ("rep2", "BOT_nc14b_rep2_CPM")],
    "BOT_hR"    : [("rep1", "BOT_hR_nc14b_rep1_CPM"),
                   ("rep2", "BOT_hR_nc14b_rep2_CPM")],
    "BOTR"      : [("rep1", "BOTR_nc14b_rep1_CPM"),
                   ("rep2", "BOTR_nc14b_rep2_CPM")],
    "BOTC"      : [("rep1", "BOTC_nc14b_rep1_CPM"),
                   ("rep2", "BOTC_nc14b_rep2_CPM")],
    "BOTC_oR"   : [("rep1", "BOTC_oR_nc14b_rep1_CPM"),
                   ("rep2", "BOTC_oR_nc14b_rep2_CPM")],
    "yw_WT"     : [("rep1", "yw_WT_nc14b_rep1_CPM"),
                   ("rep2", "yw_WT_nc14b_rep2_CPM"),
                   ("rep3", "yw_WT_nc14b_rep3_CPM")],
    # ── CVM (super-late stage 10/11) ────────────────────────────────────────
    # Stems are raw-CPM base; --bg_sub swaps _CPM → _bgSubCPM at runtime
    # BOTv st10: add rep stems here once files are confirmed
    "BOTv_st10" : [],
    "BOTCv_st10_11": [("rep1", "BOTCv_st10_11_rep1_CPM"),
                       ("rep2", "BOTCv_st10_11_rep2_CPM"),
                       ("rep3", "BOTCv_st10_11_rep3_CPM")],

    # ── Additional confirmed timepoints ─────────────────────────────────────
    "BOTv_gastr"     : [("rep1", "BOTv_gastr_rep1_CPM"),
                        ("rep2", "BOTv_gastr_rep2_CPM")],
    "BOT_nc14d"      : [("rep1", "BOT_nc14d_rep1_CPM"),
                        ("rep2", "BOT_nc14d_rep2_CPM")],
    "BOT_hR_gastr"   : [("rep1", "BOT_hR_gastr_rep1_CPM"),
                        ("rep2", "BOT_hR_gastr_rep2_CPM")],
    "BOTR_gastr"     : [("rep1", "BOTR_gastr_rep1_CPM")],

    # ── Additional confirmed timepoints ─────────────────────────────────────
    "BOTCv_gastr"    : [("rep1", "BOTCv_gastr_rep1_CPM")],
    "BOTC_nc14late"  : [("rep1", "BOTC_nc14late_rep1_CPM"),
                        ("rep2", "BOTC_nc14late_rep2_CPM")],
    "BOTC_oR_late"   : [("rep1", "BOTC_oR_late_rep1_CPM"),
                        ("rep2", "BOTC_oR_late_rep2_CPM")],
}

# ── Mode presets ───────────────────────────────────────────────────────────
MODE_GENOTYPES: Dict[str, List[str]] = {
    # "all": UNCHANGED default -- nc14b/nc14late confirmed anchors only, one
    # row per genotype, same as before. Use "all_timepoints" below for the
    # expanded view across every confirmed timepoint.
    "all"            : ["BOTv", "BOTv_late", "BOTCv", "BOTCv_late",
                        "BOT", "BOT_hR", "BOTR", "BOTC", "BOTC_oR"],
    "ventralized"    : ["BOTv", "BOTv_late", "BOTCv", "BOTCv_late"],
    "non_ventralized": ["BOT", "BOT_hR", "BOTR", "BOTC", "BOTC_oR"],

    # NEW: every genotype x every timepoint in sample_manifest.sh (all
    # confirmed).
    "all_timepoints" : [
        "yw_WT",
        "BOTv", "BOTv_late", "BOTv_gastr",
        "BOTCv", "BOTCv_late", "BOTCv_gastr",
        "BOT", "BOT_nc14d",
        "BOT_hR", "BOT_hR_gastr",
        "BOTR", "BOTR_gastr",
        "BOTC", "BOTC_nc14late",
        "BOTC_oR", "BOTC_oR_late",
    ],
    "ventralized_timepoints"     : ["BOTv", "BOTv_late", "BOTv_gastr",
                                    "BOTCv", "BOTCv_late",
                                    "BOTCv_gastr"],
    "non_ventralized_timepoints" : ["yw_WT", "BOT", "BOT_nc14d",
                                    "BOT_hR", "BOT_hR_gastr",
                                    "BOTR", "BOTR_gastr",
                                    "BOTC", "BOTC_nc14late",
                                    "BOTC_oR", "BOTC_oR_late"],

    # CVM: BOTv + BOTCv across nc14b / nc14late / st10_11
    # CVM_solo: BOTCv only across all three timepoints
    # Use --bg_sub to switch all stems to _bgSubCPM variants
    "CVM"        : ["BOTv", "BOTv_late", "BOTCv", "BOTCv_late",
                    "BOTCv_st10_11"],
    "CVM_solo"   : ["BOTCv", "BOTCv_late", "BOTCv_st10_11"],
}


def apply_stage_filter(genotypes: List[str], stage_filter: str) -> List[str]:
    """Restrict a resolved genotype-key list to one developmental stage,
    using each key's tp_group in MERGED_CONFIG. 'nc14d' matches 'nc14late',
    'nc14d', AND 'late' tp_groups since they're all the same stage under
    different naming conventions (see --stage_filter help). 'all' is a
    no-op. Unknown genotype keys are dropped silently here -- build_row_list
    already warns about those when it tries to resolve a BigWig path."""
    if stage_filter == "all":
        return genotypes
    keep = []
    for g in genotypes:
        cfg = MERGED_CONFIG.get(g)
        if cfg is None:
            continue
        tp_group = cfg[3]
        if stage_filter == "nc14b" and tp_group == "nc14b":
            keep.append(g)
        elif stage_filter == "nc14d" and tp_group in ("nc14late", "nc14d", "late"):
            keep.append(g)
    return keep


# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def _find_bw(bigwig_dir: str, stem: str) -> str:
    """Return the first existing .bigWig / .bw path, or the .bigWig path if absent."""
    for ext in (".bigWig", ".bw"):
        p = os.path.join(bigwig_dir, stem + ext)
        if os.path.isfile(p):
            return p
    return os.path.join(bigwig_dir, stem + ".bigWig")   # will trigger warning


def _apply_bgsub(stem: str) -> str:
    """
    Swap the CPM suffix for its background-subtracted equivalent.

    Mapping (order matters — most specific first):
      _mean_CPM  →  _mean_bgSubCPM
      _CPM       →  _bgSubCPM

    Stems that already contain 'bgSub' are returned unchanged, so the
    function is idempotent and safe to call unconditionally.
    """
    if "bgSub" in stem:
        return stem
    if stem.endswith("_mean_CPM"):
        return stem[:-len("_mean_CPM")] + "_mean_bgSubCPM"
    if "_mean_CPM" in stem:
        return stem.replace("_mean_CPM", "_mean_bgSubCPM")
    if stem.endswith("_CPM"):
        return stem[:-len("_CPM")] + "_bgSubCPM"
    return stem.replace("_CPM", "_bgSubCPM")


def fetch_coverage(bw_path: str, chrom: str, start: int, end: int,
                   n_bins: int = 1000) -> np.ndarray:
    """Bin-average CPM across [start, end).  Returns zeros on any failure."""
    if not os.path.isfile(bw_path):
        print(f"  WARNING: BigWig not found: {bw_path}", file=sys.stderr)
        return np.zeros(n_bins)
    bw = pyBigWig.open(bw_path)
    try:
        vals = bw.stats(chrom, start, end, type="mean", nBins=n_bins)
        arr = np.array([v if v is not None else 0.0 for v in vals], dtype=float)
    except Exception as e:
        print(f"  WARNING: {chrom}:{start}-{end} failed in {bw_path}: {e}", file=sys.stderr)
        arr = np.zeros(n_bins)
    finally:
        bw.close()
    return arr


def _smooth_coverage(arr: np.ndarray, sigma: float) -> np.ndarray:
    """
    Gaussian-smooth a coverage array in-place (returns a new array).

    sigma is in bins (not bp).  At 1000 bins over a typical ~20 kb locus,
    1 bin ≈ 20 bp, so:
        sigma=3   →  ~60 bp  (light denoise, preserves narrow peaks)
        sigma=5   →  ~100 bp (good default for bgSub spikiness)
        sigma=10  →  ~200 bp (browser-style smooth)
        sigma=25  →  ~500 bp (heavy smoothing, broad shapes only)

    Uses reflect padding to avoid edge artifacts.
    scipy is optional — falls back to a pure-numpy running average if absent.
    """
    if sigma <= 0:
        return arr
    try:
        from scipy.ndimage import gaussian_filter1d
        return gaussian_filter1d(arr, sigma=sigma, mode="reflect")
    except ImportError:
        # numpy fallback: simple box-car of width ~2*sigma
        w = max(1, int(round(sigma * 2)))
        kernel = np.ones(w) / w
        pad = w // 2
        padded = np.pad(arr, pad, mode="reflect")
        return np.convolve(padded, kernel, mode="valid")[:len(arr)]


def load_rnaseq_tsv(tsv_path: Optional[str]) -> Optional[Dict[str, list]]:
    """genotype<TAB>CPM_value rows → {genotype: [values]}."""
    if not tsv_path or not os.path.isfile(tsv_path):
        return None
    import pandas as pd
    try:
        df = pd.read_csv(tsv_path, sep="\t", comment="#",
                         names=["genotype", "value"])
        return {g: grp["value"].tolist() for g, grp in df.groupby("genotype")}
    except Exception as e:
        print(f"  WARNING: Could not load {tsv_path}: {e}", file=sys.stderr)
        return None


def load_annotation_bed(bed_path: Optional[str], chrom: str,
                        start: int, end: int) -> List[dict]:
    """
    BED6 or BED12 → list of feature dicts overlapping [start, end).
    BED12 block columns parsed into exon list.
    """
    if not bed_path or not os.path.isfile(bed_path):
        return []
    feats = []
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if cols[0] != chrom:
                continue
            fs, fe = int(cols[1]), int(cols[2])
            if fe <= start or fs >= end:
                continue
            name   = cols[3] if len(cols) > 3 else ""
            strand = cols[5] if len(cols) > 5 else "+"
            feat   = {"start": fs, "end": fe, "name": name, "strand": strand, "exons": []}
            if len(cols) >= 12:
                try:
                    n       = int(cols[9])
                    sizes   = [int(x) for x in cols[10].rstrip(",").split(",")]
                    offsets = [int(x) for x in cols[11].rstrip(",").split(",")]
                    feat["exons"] = [(fs + offsets[i], fs + offsets[i] + sizes[i])
                                     for i in range(n)]
                except Exception:
                    pass
            feats.append(feat)
    return feats
# GTF INSERTED


def load_annotation_gtf(gtf_path, chrom, start, end, gene_name=None, debug=False):
    """
    Parse a FlyBase GFF3 (or GTF) and return feature dicts for gene_name
    overlapping [start, end).

    Two-pass: pass 1 collects mRNA/transcript lines → transcript_id→gene map;
    pass 2 collects exon/CDS/UTR via Parent= (GFF3) or transcript_id= (GTF).
    Handles chr-prefix mismatch automatically.
    """
    if not gtf_path or not os.path.isfile(gtf_path):
        return []

    def _pa(attr):
        d = {}
        if "=" in attr and '"' not in attr:
            for f in attr.strip().rstrip(";").split(";"):
                f = f.strip()
                if "=" in f:
                    k, _, v = f.partition("="); d[k.strip()] = v.strip()
        else:
            for f in attr.strip().rstrip(";").split(";"):
                f = f.strip()
                if not f: continue
                p = f.split(None, 1)
                if len(p) == 2: d[p[0]] = p[1].strip('"\'')
        return d

    def _cm(lc, qc):
        if lc == qc: return True
        if qc.startswith("chr") and lc == qc[3:]: return True
        if not qc.startswith("chr") and lc == "chr" + qc: return True
        return False

    def _gs(attrs, ftype):
        raw = attrs.get("gene_name") or attrs.get("gene_symbol") or ""
        if not raw:
            n = attrs.get("Name", "")
            if ftype in ("mRNA","transcript","ncRNA","lncRNA") and "-R" in n:
                raw = n.rsplit("-R", 1)[0]
            elif ftype in ("mRNA","transcript","ncRNA","lncRNA"):
                raw = n
        if not raw:
            raw = attrs.get("gene_id","")
        if raw.startswith("FBgn") or raw.startswith("FBtr"): raw = ""
        return raw

    TTYPES = {"mRNA","transcript","ncRNA","lncRNA","pre_miRNA",
              "snoRNA","snRNA","tRNA","rRNA","pseudogene"}
    tid_info = {}
    try:
        with open(gtf_path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip(): continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9: continue
                lc, _, ft, s0, e0, _, strand, _, att = cols[:9]
                if ft not in TTYPES or not _cm(lc, chrom): continue
                fs, fe = int(s0), int(e0)
                if fe < start or fs > end: continue
                attrs = _pa(att)
                tid = attrs.get("transcript_id") or attrs.get("ID","")
                if not tid: continue
                gname = _gs(attrs, ft)
                if gene_name:
                    if gname.lower() != gene_name.lower() and gene_name.lower() not in gname.lower():
                        continue
                tid_info[tid] = {"name": gname or gene_name or "",
                                 "strand": strand, "start": fs, "end": fe}
    except Exception as e:
        print(f"  WARNING: GTF pass-1 error: {e}", file=sys.stderr); return []

    if debug:
        print(f"  [debug] Pass 1: {len(tid_info)} transcript(s) for {gene_name!r} on {chrom}", file=sys.stderr)
        for t,v in list(tid_info.items())[:6]: print(f"    {t}: {v}", file=sys.stderr)

    if not tid_info:
        print(f"  WARNING: no transcripts found for {gene_name!r} on {chrom} in {gtf_path}", file=sys.stderr)
        return []

    feats = {tid: {**info, "transcript_id": tid,
                   "exons":[], "utrs":[], "has_cds":False, "efb":[]}
             for tid, info in tid_info.items()}

    FTYPES = {"exon","CDS","UTR","five_prime_utr","three_prime_utr"}
    try:
        with open(gtf_path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip(): continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9: continue
                lc, _, ft, s0, e0, _, _, _, att = cols[:9]
                if ft not in FTYPES or not _cm(lc, chrom): continue
                fs, fe = int(s0), int(e0)
                if fe < start or fs > end: continue
                attrs = _pa(att)
                parent = attrs.get("transcript_id") or attrs.get("Parent","")
                for tid in [t.strip() for t in parent.split(",")]:
                    if tid not in feats: continue
                    f = feats[tid]
                    f["start"] = min(f["start"], fs); f["end"] = max(f["end"], fe)
                    if ft == "CDS":
                        f["exons"].append((fs,fe)); f["has_cds"] = True
                    elif ft in ("UTR","five_prime_utr","three_prime_utr"):
                        f["utrs"].append((fs,fe))
                    elif ft == "exon":
                        f["efb"].append((fs,fe))
    except Exception as e:
        print(f"  WARNING: GTF pass-2 error: {e}", file=sys.stderr); return []

    result = list(feats.values())
    for f in result:
        if not f["has_cds"] and f["efb"]: f["exons"] = f["efb"]
        f.pop("has_cds",None); f.pop("efb",None)
        f["exons"] = sorted(set(f["exons"]))
        f["utrs"]  = sorted(set(f["utrs"]))

    if debug:
        print(f"  [debug] Pass 2:", file=sys.stderr)
        for f in result:
            print(f"    {f['name']} ({f['transcript_id']}) "
                  f"exons={len(f['exons'])} utrs={len(f['utrs'])} strand={f['strand']}", file=sys.stderr)

    result.sort(key=lambda f: (len(f["exons"]), f["end"]-f["start"]), reverse=True)
    return result[:6]

# ─────────────────────────────────────────────────────────────────────────────
# MACS2 PEAKS  — load one or more narrowPeak / BED files and draw a combined row
# ─────────────────────────────────────────────────────────────────────────────

def load_macs2_peaks(peak_files, chrom: str, start: int, end: int) -> List[dict]:
    """
    Load peaks overlapping [start, end) from one or more MACS2 narrowPeak or
    BED files.  Returns a list of dicts: {start, end, score, name, source}.

    Accepted formats
    ----------------
    • MACS2 narrowPeak  (BED6+4): cols 0-2 mandatory; col 4 = score (0-1000);
      col 6 = signalValue; col 7 = -log10(pval); col 8 = -log10(qval)
    • Plain BED3 / BED6

    peak_files can be a single path string, or a list of path strings, or None.
    """
    if not peak_files:
        return []
    if isinstance(peak_files, str):
        peak_files = [peak_files]

    peaks = []
    for path in peak_files:
        if not path or not os.path.isfile(path):
            print(f"  WARNING: peaks file not found: {path}", file=sys.stderr)
            continue
        src = os.path.basename(path)
        with open(path) as fh:
            for line in fh:
                if line.startswith(("#", "track", "browser")) or not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 3:
                    continue
                if cols[0] != chrom:
                    # handle chr-prefix mismatch
                    c = cols[0]
                    if chrom.startswith("chr") and c == chrom[3:]:
                        pass
                    elif not chrom.startswith("chr") and c == "chr" + chrom:
                        pass
                    else:
                        continue
                ps, pe = int(cols[1]), int(cols[2])
                if pe <= start or ps >= end:
                    continue
                score = float(cols[4]) if len(cols) > 4 else 500.0
                name  = cols[3]         if len(cols) > 3 else ""
                # For narrowPeak, prefer signalValue (col 6) as display score
                signal = float(cols[6]) if len(cols) > 6 else score
                peaks.append({"start": ps, "end": pe,
                               "score": score, "signal": signal,
                               "name": name, "source": src})
    return peaks


def load_named_peaks(named_peak_entries, chrom: str, start: int, end: int
                     ) -> List[Tuple[str, str, List[dict]]]:
    """
    Load per-genotype / per-sample MACS2 peak files for multi-row display.

    Parameters
    ----------
    named_peak_entries : list of str
        Each entry is  "Label:color=path"  or  "Label=path"  (color falls back
        to a tab10 cycle).
        Examples:
          "BOTv=#2ca02c=/data/BOTv_peaks.narrowPeak"
          "BOTCv=#e377c2=/data/BOTCv_peaks.narrowPeak"
          "BOTv=/data/BOTv_peaks.narrowPeak"   ← color auto-assigned

    Returns
    -------
    List of (label, hex_color, peaks_list) tuples in input order.
    """
    if not named_peak_entries:
        return []
    if isinstance(named_peak_entries, str):
        named_peak_entries = [named_peak_entries]

    import itertools
    _color_cycle = itertools.cycle([
        "#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
        "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf"])

    result = []
    for entry in named_peak_entries:
        if not entry:
            continue
        # Parse "Label=#hex=/path"  or  "Label=/path"
        parts = entry.split("=", 2)
        if len(parts) == 3:
            label, color, path = parts
        elif len(parts) == 2:
            label, path = parts
            color = next(_color_cycle)
        else:
            # bare path — use filename stem as label
            path  = parts[0]
            label = os.path.splitext(os.path.basename(path))[0]
            color = next(_color_cycle)

        label = label.strip()
        color = color.strip()
        path  = path.strip()

        peaks = load_macs2_peaks(path, chrom, start, end)
        result.append((label, color, peaks))

    return result


def _draw_peaks_row(ax, peaks: List[dict], chrom: str, start: int, end: int,
                    label: str = "MACS2 peaks",
                    color: str = "#333333",
                    score_threshold: float = 0.0,
                    label_color: Optional[str] = None,
                    ticks: Optional[List[int]] = None,
                    uniform_height: bool = False):
    """
    Draw a compact MACS2-peaks row (like a browser track).
    Peaks are drawn as filled rectangles; taller = higher signal.
    score_threshold filters by the MACS2 'score' column (0-1000 scale).

    uniform_height : if True, every visible peak is drawn at the same fixed
        height instead of being min/max-scaled against the 'signal' values of
        whatever else happens to be visible in this window. Curated locus
        lists (Merged DARs, Rescued sets, etc.) don't carry a meaningful
        signal gradient, and per-window rescaling makes the same locus look
        a different size in different plots/tracks -- use this for those.
        True MACS2 peak calls (real signal strength) should keep it False.
    """
    ax.set_xlim(start, end)
    ax.set_ylim(0, 1)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(axis="both", left=False, bottom=False,
                   labelleft=False, labelbottom=False)
    if ticks:
        _draw_tick_gridlines(ax, ticks, zorder=0)
    ax.axhline(0.05, color="#cccccc", linewidth=0.5, zorder=0)

    lc = label_color or color
    ax.text(-0.008, 0.55, label,
            transform=ax.transAxes,
            fontsize=6.5, color=lc, fontweight="bold",
            ha="right", va="center", clip_on=False)

    visible = [p for p in peaks if p["score"] >= score_threshold]
    if not visible:
        return

    if uniform_height:
        # Fixed block height -- no per-window rescaling, so the same locus
        # renders identically everywhere this track is drawn.
        UNIFORM_H = 0.75
        for p in visible:
            ps2 = max(p["start"], start)
            pe2 = min(p["end"],   end)
            if pe2 <= ps2:
                continue
            ax.add_patch(mpatches.Rectangle(
                (ps2, 0.05), pe2 - ps2, UNIFORM_H,
                facecolor=color, edgecolor="none", alpha=0.80, zorder=3))
        return

    # Normalise signal to [0.15, 0.90] height range
    sigs = [p["signal"] for p in visible]
    sig_min, sig_max = min(sigs), max(sigs)
    sig_range = sig_max - sig_min if sig_max != sig_min else 1.0

    for p in visible:
        ps2 = max(p["start"], start)
        pe2 = min(p["end"],   end)
        if pe2 <= ps2:
            continue
        h = 0.15 + 0.75 * (p["signal"] - sig_min) / sig_range
        ax.add_patch(mpatches.Rectangle(
            (ps2, 0.05), pe2 - ps2, h,
            facecolor=color, edgecolor="none", alpha=0.80, zorder=3))


# ─────────────────────────────────────────────────────────────────────────────
# TF BINDING SITES  — load a multi-TF BED and draw a colour-coded tick row
# ─────────────────────────────────────────────────────────────────────────────

# Default TF colour palette (matches the reference figure legend).
# Users can override via --tf_colors  TF=hex  TF=hex ...
TF_COLORS_DEFAULT: Dict[str, str] = {
    "Doc"   : "#d62728",   # red
    "Twist" : "#ff7f0e",   # orange
    "Sna"   : "#1f77b4",   # blue
    "Twi"   : "#2ca02c",   # green
    "Zfh1"  : "#9467bd",   # purple
    # extend freely
}


def load_tf_sites(tf_files, chrom: str, start: int, end: int,
                  name_col: int = 3) -> Dict[str, List[dict]]:
    """
    Load TF binding sites from one or more BED files.

    Two input modes
    ---------------
    1. Single multi-TF BED where col `name_col` (0-based, default 3) encodes
       the TF name  →  {'TF1': [{start, end, score}, ...], ...}
    2. Per-TF files supplied as  "TFname=path"  strings
       e.g.  ["Doc=/data/Doc_peaks.bed", "Sna=/data/Sna_peaks.bed"]
       The TF name is taken from the prefix before '='.

    tf_files : None | str | List[str]
        Each entry is either a plain BED path (multi-TF) or "Name=path".
    """
    if not tf_files:
        return {}
    if isinstance(tf_files, str):
        tf_files = [tf_files]

    result: Dict[str, List[dict]] = {}

    for entry in tf_files:
        if not entry:
            continue
        # Detect "Name=path" syntax
        if "=" in entry and not os.path.isfile(entry):
            tf_name, _, path = entry.partition("=")
            tf_name = tf_name.strip()
        else:
            tf_name = None
            path = entry

        if not os.path.isfile(path):
            print(f"  WARNING: TF sites file not found: {path}", file=sys.stderr)
            continue

        with open(path) as fh:
            for line in fh:
                if line.startswith(("#", "track", "browser")) or not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 3:
                    continue
                c = cols[0]
                if c != chrom:
                    if chrom.startswith("chr") and c == chrom[3:]: pass
                    elif not chrom.startswith("chr") and c == "chr" + chrom: pass
                    else: continue
                ps, pe = int(cols[1]), int(cols[2])
                if pe <= start or ps >= end:
                    continue
                score = float(cols[4]) if len(cols) > 4 else 500.0
                # Determine TF name
                if tf_name:
                    tf = tf_name
                elif len(cols) > name_col and cols[name_col].strip():
                    tf = cols[name_col].strip()
                else:
                    tf = os.path.splitext(os.path.basename(path))[0]
                result.setdefault(tf, []).append(
                    {"start": ps, "end": pe, "score": score})

    return result


def _draw_tf_row(ax, tf_sites: Dict[str, List[dict]],
                 chrom: str, start: int, end: int,
                 tf_colors: Optional[Dict[str, str]] = None,
                 label: str = "TF binding sites",
                 ticks: Optional[List[int]] = None):
    """
    Draw a multi-TF tick row.

    Each TF is assigned a distinct row within the panel; ticks are coloured per
    TF using tf_colors (falls back to TF_COLORS_DEFAULT, then a tab10 cycle for
    unknown TFs).
    """
    if not tf_sites:
        return

    ax.set_xlim(start, end)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(axis="both", left=False, bottom=False,
                   labelleft=False, labelbottom=False)
    if ticks:
        _draw_tick_gridlines(ax, ticks, zorder=0)

    # label
    ax.text(-0.008, 0.50, label,
            transform=ax.transAxes,
            fontsize=6.5, color="#444444",
            ha="right", va="center", clip_on=False)

    colors = dict(TF_COLORS_DEFAULT)
    if tf_colors:
        colors.update(tf_colors)

    # Colour cycle for unknown TFs
    import itertools
    _cycle = itertools.cycle([
        "#e41a1c","#377eb8","#4daf4a","#984ea3",
        "#ff7f00","#a65628","#f781bf","#999999"])
    _assigned: Dict[str, str] = {}
    def _get_color(tf):
        if tf in colors: return colors[tf]
        if tf not in _assigned: _assigned[tf] = next(_cycle)
        return _assigned[tf]

    tfs_ordered = sorted(tf_sites.keys())
    n = len(tfs_ordered)
    if n == 0: return

    ax.set_ylim(0, n)

    for row_i, tf in enumerate(tfs_ordered):
        y_center = row_i + 0.5
        c = _get_color(tf)

        # faint horizontal rule
        ax.axhline(y_center, color="#eeeeee", linewidth=0.4, zorder=0)

        # TF name on right margin
        ax.text(1.002, y_center / n, tf,
                transform=ax.transAxes,
                fontsize=5.8, color=c,
                ha="left", va="center", clip_on=False,
                fontweight="bold")

        for site in tf_sites[tf]:
            ps2 = max(site["start"], start)
            pe2 = min(site["end"],   end)
            if pe2 <= ps2: continue
            w = pe2 - ps2
            # Minimum visible width = 0.3% of locus
            w = max(w, (end - start) * 0.003)
            ax.add_patch(mpatches.Rectangle(
                (ps2, y_center - 0.38), w, 0.76,
                facecolor=c, edgecolor="none", alpha=0.85, zorder=3))


def _draw_combined_annotation_and_tracks(
    fig, gs, row_start: int,
    gene: str, chrom: str, start: int, end: int,
    annotation_features: Optional[List[dict]],
    macs2_peaks: Optional[List[dict]],
    named_peaks: Optional[List[Tuple[str, str, List[dict]]]] = None,
    tf_sites: Optional[Dict[str, List[dict]]] = None,
    tf_colors: Optional[Dict[str, str]] = None,
    peaks_color: str = "#333333",
    peaks_score_threshold: float = 0.0,
    ticks: Optional[List[int]] = None,
    tick_labels: Optional[List[str]] = None,
):
    """
    Draw the bottom annotation block.

    Layout (rows allocated in order):
      row_start+0        : RefGene annotation
      +1…+N              : one row per named peaks set (--macs2_named_peaks)
      +N+1               : merged/union MACS2 peaks row  (--macs2_peaks, optional)
      last               : TF binding sites row           (--tf_sites, optional)

    `ticks`/`tick_labels` are the shared, adaptively-formatted x-axis tick
    positions computed once in plot_tracks(); passing them through here keeps
    every sub-row's gridlines aligned to the same positions as the coverage
    tracks above and the labeled ruler on the annotation row itself.
    """
    idx = row_start
    ax_ann = fig.add_subplot(gs[idx, 0])
    _draw_annotation_row(ax_ann, gene, chrom, start, end, annotation_features,
                        ticks=ticks, tick_labels=tick_labels)
    idx += 1

    # Per-genotype named peak rows (colored, one row each)
    if named_peaks:
        for label, color, peaks in named_peaks:
            ax_pk = fig.add_subplot(gs[idx, 0])
            _draw_peaks_row(ax_pk, peaks, chrom, start, end,
                            label=label, color=color,
                            label_color=color,
                            # Named tracks (Merged DARs / Rescued sets) are
                            # curated locus lists, not raw peak calls -- their
                            # 'score' column often encodes direction (signed
                            # fold-change), not confidence. Don't apply the
                            # generic peaks_score_threshold here, or every
                            # negative-scored ("closing") locus silently
                            # vanishes. The grey MACS2-peaks row below keeps
                            # the real threshold since that's an actual
                            # confidence cutoff.
                            score_threshold=float("-inf"),
                            # Same reasoning extends to height: these rows'
                            # 'score' isn't a comparable signal-strength
                            # value, and per-window min/max rescaling made
                            # the same locus look a different size depending
                            # on what else was in view. Equal-sized blocks
                            # instead -- presence/absence is the information
                            # these tracks carry, not magnitude.
                            uniform_height=True,
                            ticks=ticks)
            idx += 1

    # Single merged/union peaks row (grey by default)
    if macs2_peaks is not None:
        ax_pk = fig.add_subplot(gs[idx, 0])
        _draw_peaks_row(ax_pk, macs2_peaks, chrom, start, end,
                        label="MACS2 peaks", color=peaks_color,
                        label_color="#444444",
                        score_threshold=peaks_score_threshold,
                        ticks=ticks)
        idx += 1

    if tf_sites:
        ax_tf = fig.add_subplot(gs[idx, 0])
        _draw_tf_row(ax_tf, tf_sites, chrom, start, end,
                     tf_colors=tf_colors, ticks=ticks)
        idx += 1

    return idx


def _tf_legend(fig, tf_sites: Dict[str, List[dict]],
               tf_colors: Optional[Dict[str, str]] = None,
               x: float = 0.235, y: float = 0.0,
               fontsize: float = 5.5):
    """Draw a compact TF colour legend below the figure."""
    if not tf_sites:
        return
    colors = dict(TF_COLORS_DEFAULT)
    if tf_colors:
        colors.update(tf_colors)
    import itertools
    _cycle = itertools.cycle([
        "#e41a1c","#377eb8","#4daf4a","#984ea3",
        "#ff7f00","#a65628","#f781bf","#999999"])
    _assigned: Dict[str, str] = {}
    def _gc(tf):
        if tf in colors: return colors[tf]
        if tf not in _assigned: _assigned[tf] = next(_cycle)
        return _assigned[tf]

    handles = [mpatches.Patch(facecolor=_gc(tf), label=tf)
               for tf in sorted(tf_sites.keys())]
    fig.legend(handles=handles, loc="lower left",
               bbox_to_anchor=(x, y),
               ncol=len(handles),
               fontsize=fontsize,
               frameon=False,
               handlelength=1.2,
               handleheight=0.8,
               columnspacing=0.8)


def _clean_ymax(raw: float) -> float:
    """Round up to a visually clean value."""
    if raw <= 0:
        return 1.0
    mag = 10 ** int(np.floor(np.log10(raw)))
    for mult in [1, 2, 5, 10]:
        c = mult * mag
        if c >= raw:
            return float(c)
    return float(10 * mag)


def _draw_axis_break(ax, color: str = "#444444", at: str = "top") -> None:
    """Draw a small diagonal '//' break mark on the left spine — at the top
    of the axes (at="top", the outlier's own bottom segment sits below) or
    at the bottom (at="bottom", used on the compressed upper segment of a
    split/broken axis, where the break sits just above where the lower
    segment picks up). Positioned in axes-fraction coordinates so it looks
    identical regardless of that track's data range.
    """
    kwargs = dict(transform=ax.transAxes, color=color, linewidth=0.9,
                  clip_on=False, solid_capstyle="butt")
    centers = (0.90, 1.04) if at == "top" else (-0.04, 0.10)
    for y_center in centers:
        ax.plot([-0.018, 0.018], [y_center - 0.05, y_center + 0.05], **kwargs)


def _hex_lighten(hex_color: str, factor: float = 0.45) -> str:
    """Blend a hex colour toward white by `factor` (0=original, 1=white)."""
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    r2 = int(r + (255 - r) * factor)
    g2 = int(g + (255 - g) * factor)
    b2 = int(b + (255 - b) * factor)
    return f"#{r2:02x}{g2:02x}{b2:02x}"


# ─────────────────────────────────────────────────────────────────────────────
# TRACK-ROW BUILDER
# Expands a list of genotype keys into an ordered list of row descriptors,
# interleaving replicate rows (if show_replicates=True) before each merged row.
#
# Each row descriptor is a dict:
#   kind        "replicate" | "merged"
#   geno_key    key into MERGED_CONFIG
#   label       display string
#   color       hex fill color
#   bw_path     resolved BigWig path
#   rep_label   e.g. "rep1"  (replicates only)
#   is_last_in_group  True for the merged row at the bottom of each genotype block
# ─────────────────────────────────────────────────────────────────────────────

def build_row_list(genotypes: List[str], bigwig_dir: str,
                   show_replicates: bool = True,
                   bg_sub: bool = False) -> List[dict]:
    rows = []
    for geno in genotypes:
        cfg = MERGED_CONFIG.get(geno)
        if cfg is None:
            print(f"  WARNING: Unknown genotype key '{geno}' — skipping", file=sys.stderr)
            continue
        label, color, merged_stem, tp_grp, _ = cfg
        light_color = _hex_lighten(color, 0.50)

        if bg_sub:
            merged_stem = _apply_bgsub(merged_stem)

        if show_replicates:
            rep_list = REPLICATE_CONFIG.get(geno, [])
            for rep_label, rep_stem in rep_list:
                if bg_sub:
                    rep_stem = _apply_bgsub(rep_stem)
                rows.append({
                    "kind"            : "replicate",
                    "geno_key"        : geno,
                    "label"           : f"  {rep_label}",   # indented
                    "color"           : light_color,
                    "outline_color"   : color,
                    "bw_path"         : _find_bw(bigwig_dir, rep_stem),
                    "rep_label"       : rep_label,
                    "tp_group"        : tp_grp,
                    "is_last_in_group": False,
                })

        rows.append({
            "kind"            : "merged",
            "geno_key"        : geno,
            "label"           : label,
            "color"           : color,
            "outline_color"   : color,
            "bw_path"         : _find_bw(bigwig_dir, merged_stem),
            "rep_label"       : None,
            "tp_group"        : tp_grp,
            "is_last_in_group": True,
        })

    # Mark timepoint-group boundaries for separator lines
    for i in range(1, len(rows)):
        if rows[i]["tp_group"] != rows[i - 1]["tp_group"]:
            rows[i]["sep_above"] = True
        else:
            rows[i]["sep_above"] = False
    if rows:
        rows[0]["sep_above"] = False

    return rows


# ─────────────────────────────────────────────────────────────────────────────
# ADAPTIVE X-AXIS TICKS
#
# Picks "nice" round-number tick positions (1/2/5/10 x 10^n) scaled to the
# window size, and formats them in bp / kb / Mb depending on zoom level —
# instead of always showing e.g. "30.848 Mb" even when zoomed into a 2kb
# window. Ticks are computed once per plot and reused as shared vertical
# gridlines across every track (coverage, peaks, TF sites, gene model) so a
# peak's x-position can be traced straight down through the whole figure.
# ─────────────────────────────────────────────────────────────────────────────

def _nice_step(span: float, target_n: int = 6) -> float:
    """Return a 'nice' step size (1/2/5/10 x 10^n) giving ~target_n ticks."""
    if span <= 0:
        return 1.0
    raw_step = span / target_n
    mag = 10 ** math.floor(math.log10(raw_step))
    for m in (1, 2, 2.5, 5, 10):
        step = m * mag
        if step >= raw_step:
            return step
    return 10 * mag


def _compute_ticks(start: int, end: int, target_n: int = 6) -> List[int]:
    """Evenly spaced 'nice' tick positions covering [start, end]."""
    span = end - start
    step = _nice_step(span, target_n)
    first = math.ceil(start / step) * step
    ticks = []
    t = first
    while t <= end + 1e-6:
        ticks.append(int(round(t)))
        t += step
    if len(ticks) < 2:
        ticks = [start, end]
    return ticks


def _format_ticks(ticks: List[int], span: int) -> List[str]:
    """
    Format tick values adaptively, genome-browser style.

    IMPORTANT: these are absolute chromosome coordinates (tens of millions
    of bp for Drosophila), not offsets from the window start. Below ~200kb
    zoom we show the full comma-separated bp coordinate (e.g. "30,850,000")
    rather than truncating it into "kb" units, since dividing a ~30Mb
    absolute position by 1,000 produces an unreadable, non-genomic-looking
    number (e.g. "30850 kb"). Above ~200kb we switch to Mb, with decimal
    precision scaled to the tick spacing so adjacent ticks stay distinct.
    """
    if span < 200_000:
        return [f"{v:,}" for v in ticks]

    step = (ticks[1] - ticks[0]) if len(ticks) > 1 else span
    step_mb = step / 1_000_000
    decimals = max(0, min(6, -math.floor(math.log10(step_mb)))) if step_mb > 0 else 3
    return [f"{v / 1_000_000:.{decimals}f} Mb" for v in ticks]


def _draw_tick_gridlines(ax, ticks: List[int], zorder: float = 0):
    """Faint vertical gridline at each shared tick position."""
    for t in ticks:
        ax.axvline(t, color="#eeeeee", linewidth=0.5, zorder=zorder)


# ─────────────────────────────────────────────────────────────────────────────
# GENE ANNOTATION ROW DRAWING  (extracted so it can be called in both plot fns)
# ─────────────────────────────────────────────────────────────────────────────

def _draw_annotation_row(ax_ann, gene: str, chrom: str,
                         start: int, end: int,
                         annotation_features: Optional[List[dict]],
                         ticks: Optional[List[int]] = None,
                         tick_labels: Optional[List[str]] = None):
    """
    IGV-style gene annotation row.
    • Thin line + chevrons = intron / strand
    • Tall dark block = CDS exon
    • Short lighter block = UTR
    Multiple transcripts stacked in greedy non-overlapping lanes.
    """
    locus_len = end - start
    n_feats   = len(annotation_features) if annotation_features else 0
    lane_h    = 1.0
    y_top     = max(1.2, n_feats * lane_h) + 0.6
    y_bot     = -1.1

    ax_ann.set_xlim(start, end)
    ax_ann.set_ylim(y_bot, y_top)
    for sp in ax_ann.spines.values():
        sp.set_visible(False)
    ax_ann.tick_params(axis="y", left=False, right=False, labelleft=False)
    ax_ann.tick_params(axis="x", bottom=True, labelbottom=True,
                       labelsize=6.5, colors="#555555", length=2, width=0.5, pad=1)

    if ticks is None:
        ticks = _compute_ticks(start, end)
    if tick_labels is None:
        tick_labels = _format_ticks(ticks, locus_len)
    ax_ann.set_xticks(ticks)
    ax_ann.set_xticklabels(tick_labels)
    _draw_tick_gridlines(ax_ann, ticks, zorder=0)

    ax_ann.spines["bottom"].set_visible(True)
    ax_ann.spines["bottom"].set_color("#cccccc")
    ax_ann.spines["bottom"].set_linewidth(0.5)

    CDS_COLOR  = "#1a1a8c"
    UTR_COLOR  = "#8888cc"
    LINE_COLOR = "#333377"
    CHEV_COLOR = "#7777bb"
    CDS_HALF   = 0.38
    UTR_HALF   = 0.22

    def _chevrons(ax, x0, x1, y0, strand):
        seg = x1 - x0
        if seg <= 0: return
        n = max(1, int(seg / max(locus_len / 18, 500)))
        xs = np.linspace(x0 + seg * 0.15, x1 - seg * 0.15, n)
        dx = locus_len * 0.009 * (1 if strand == "+" else -1)
        dy = 0.12
        for cx in xs:
            px = [cx - dx, cx, cx - dx] if strand == "+" else [cx + dx, cx, cx + dx]
            py = [y0 - dy, y0, y0 + dy]
            ax.plot(px, py, color=CHEV_COLOR, linewidth=0.7,
                    solid_capstyle="round", solid_joinstyle="round", zorder=4)

    if annotation_features:
        assigned: list = []  # (fs, fe, lane)

        def _lane(fs_, fe_):
            for li in range(20):
                if not any(lfs < fe_ and lfe > fs_ and ll == li
                           for lfs, lfe, ll in assigned):
                    return li
            return 0

        for feat in annotation_features:
            fs     = max(feat["start"], start)
            fe     = min(feat["end"],   end)
            if fe <= fs: continue
            strand = feat.get("strand", "+")
            name   = feat.get("name", gene)
            exons  = feat.get("exons", [])
            utrs   = feat.get("utrs",  [])

            lane = _lane(fs, fe)
            assigned.append((fs, fe, lane))
            y0 = lane * lane_h

            # backbone
            ax_ann.plot([fs, fe], [y0, y0], color=LINE_COLOR,
                        linewidth=1.4, solid_capstyle="butt", zorder=2)

            # chevrons on intron gaps
            all_blks = sorted(set(exons + utrs)) if (exons or utrs) else [(fs, fe)]
            for i in range(len(all_blks) - 1):
                _chevrons(ax_ann, all_blks[i][1], all_blks[i+1][0], y0, strand)
            if all_blks:
                if fs < all_blks[0][0]:
                    _chevrons(ax_ann, fs, all_blks[0][0], y0, strand)
                if all_blks[-1][1] < fe:
                    _chevrons(ax_ann, all_blks[-1][1], fe, y0, strand)

            # UTR blocks (drawn first so CDS overlays)
            for us, ue in utrs:
                us2, ue2 = max(us, start), min(ue, end)
                if ue2 > us2:
                    ax_ann.add_patch(mpatches.Rectangle(
                        (us2, y0 - UTR_HALF), ue2 - us2, 2 * UTR_HALF,
                        facecolor=UTR_COLOR, edgecolor="none", linewidth=0, zorder=3))

            # CDS exon blocks
            if exons:
                for es, ee in exons:
                    es2, ee2 = max(es, start), min(ee, end)
                    if ee2 > es2:
                        ax_ann.add_patch(mpatches.Rectangle(
                            (es2, y0 - CDS_HALF), ee2 - es2, 2 * CDS_HALF,
                            facecolor=CDS_COLOR, edgecolor="none", linewidth=0, zorder=4))
            elif not utrs:
                ax_ann.add_patch(mpatches.Rectangle(
                    (fs, y0 - UTR_HALF), fe - fs, 2 * UTR_HALF,
                    facecolor=UTR_COLOR, alpha=0.7, edgecolor="none", linewidth=0, zorder=3))

            # strand terminal arrow
            ax = fe if strand == "+" else fs
            dx = locus_len * 0.013 * (1 if strand == "+" else -1)
            ax_ann.annotate("", xy=(ax + dx, y0), xytext=(ax, y0),
                            arrowprops=dict(arrowstyle="-|>", color=LINE_COLOR,
                                            lw=0.9, mutation_scale=7), zorder=5)

            # label
            lx = np.clip((fs + fe) / 2, start, end)
            ax_ann.text(lx, y0 + CDS_HALF + 0.06, name,
                        ha="center", va="bottom", fontsize=7.5,
                        style="italic", color="#1a1a1a", fontweight="medium",
                        zorder=6, clip_on=True)
    else:
        ax_ann.text((start + end) / 2, 0.15, gene,
                    ha="center", va="bottom", fontsize=8.5,
                    style="italic", color="#444444", fontweight="medium", clip_on=False)

    # scale bar
    for target in [10000, 5000, 2000, 1000, 500, 200, 100]:
        if locus_len / target >= 3:
            bar_len = target; break
    else:
        bar_len = locus_len // 4
    sb_x1 = end - bar_len - locus_len * 0.015
    sb_x2 = sb_x1 + bar_len
    sb_y  = y_bot + 0.12
    
    ax_ann.annotate("", xy=(sb_x2, sb_y), xytext=(sb_x1, sb_y),
                    arrowprops=dict(arrowstyle="-", color="#444444", lw=1.2,
                                    shrinkA=0, shrinkB=0))
    ax_ann.text((sb_x1+sb_x2)/2, sb_y + 0.06,
                f"{bar_len//1000} kb" if bar_len >= 1000 else f"{bar_len} bp",
                ha="center", va="bottom", fontsize=6.5, color="#444444")
    ax_ann.text(start, sb_y - 0.08, chrom,
                ha="left", va="top", fontsize=5.5, color="#bbbbbb", clip_on=False)


# ─────────────────────────────────────────────────────────────────────────────
# MAIN PLOTTING FUNCTION
# ─────────────────────────────────────────────────────────────────────────────

def plot_tracks(
    gene: str,
    chrom: str,
    start: int,
    end: int,
    bigwig_dir: str,
    genotypes: List[str],
    out_path: str,
    rnaseq_data: Optional[Dict[str, list]] = None,
    annotation_features: Optional[List[dict]] = None,
    annotation_gtf: Optional[str] = None,
    n_bins: int = 1000,
    shared_ymax: Optional[float] = None,
    per_track_ymax: bool = False,
    show_replicates: bool = True,
    track_height: float = 0.72,
    rep_track_height: float = 0.52,
    fig_width: float = 8.5,
    fill_alpha: float = 0.80,
    rep_fill_alpha: float = 0.55,
    title_override: Optional[str] = None,
    # ── NEW: MACS2 peak + TF annotations ───────────────────────────────────
    macs2_peak_files=None,          # merged/union: single grey row
    macs2_named_peaks=None,         # per-genotype: list of "Label=color=path" strings → one colored row each
    peaks_color: str = "#333333",
    peaks_score_threshold: float = 0.0,
    peaks_row_height: float = 0.30,
    tf_files=None,
    tf_colors: Optional[Dict[str, str]] = None,
    tf_name_col: int = 3,
    tf_row_height: float = 0.22,
    # ── background subtraction ──────────────────────────────────────────────
    bg_sub: bool = False,           # swap _CPM → _bgSubCPM in all stem lookups
    # ── smoothing ───────────────────────────────────────────────────────────
    smooth_sigma: float = 0.0,      # Gaussian sigma in bins; 0 = no smoothing
    # ── broken/clipped y-axis for one or more outlier genotypes ────────────
    # Genotype key(s) whose MERGED track should be clipped to the same
    # shared ceiling everyone else uses, rather than that outlier setting
    # the ceiling for the whole panel. The shared ceiling is computed from
    # every OTHER track only. The clipped track gets a small "//" break
    # mark on its left spine plus a numeric label showing its true,
    # unclipped peak value, so the real magnitude is never lost, just
    # kept out of the shared scale.
    break_genotypes: Optional[List[str]] = None,
    break_label_color: str = "#cc0000",
):
    """
    Stacked coverage tracks with optional per-replicate sub-tracks.

    Track order within each genotype block (when show_replicates=True):
      rep1  (light fill, outline only, smaller height)
      rep2
      ...
      merged  (full color, bolder, taller — carries the genotype label)

    A dashed separator spans the full width between nc14b and nc14late blocks.
    """

    rows = build_row_list(genotypes, bigwig_dir, show_replicates=show_replicates,
                          bg_sub=bg_sub)
    if not rows:
        print(f"  ERROR: no valid genotype rows for {gene}", file=sys.stderr)
        return

    # Resolve annotation: GTF takes priority over pre-loaded BED features
    if annotation_gtf and not annotation_features:
        annotation_features = load_annotation_gtf(
            annotation_gtf, chrom, start, end, gene_name=gene)

    # ── Load MACS2 peaks and TF sites ─────────────────────────────────────
    macs2_peaks = load_macs2_peaks(macs2_peak_files, chrom, start, end) \
                  if macs2_peak_files else None
    named_peaks = load_named_peaks(macs2_named_peaks, chrom, start, end) \
                  if macs2_named_peaks else None
    tf_sites    = load_tf_sites(tf_files, chrom, start, end,
                                name_col=tf_name_col) \
                  if tf_files else None

    # ── shared adaptive x-axis ticks (bp/kb/Mb depending on window size) ───
    # Computed once and reused as gridlines across every track so a peak's
    # position can be traced from the coverage tracks down through the
    # peaks/TF rows to the labeled ruler on the annotation row.
    shared_ticks = _compute_ticks(start, end)
    shared_tick_labels = _format_ticks(shared_ticks, end - start)

    # ── broken-axis setup ────────────────────────────────────────────────
    break_set = set(break_genotypes) if break_genotypes else set()
    if break_set and per_track_ymax:
        print("  WARNING: --break_genotypes requires a shared y-axis ceiling "
              "to mean anything (that's the whole point — everyone else's "
              "scale, with the outlier clipped to it). Ignoring "
              "--per_track_ymax for this plot.", file=sys.stderr)
        per_track_ymax = False

    # ── load all coverage ──────────────────────────────────────────────────
    x = np.linspace(start, end, n_bins)
    coverages: Dict[str, np.ndarray] = {}
    global_max = 0.0            # ceiling driver — excludes broken genotypes
    break_true_max: Dict[str, float] = {}   # geno_key -> real unclipped max
    for row in rows:
        path = row["bw_path"]
        key  = path   # use path as cache key (same file may appear via merged + rep)
        if key not in coverages:
            raw = fetch_coverage(path, chrom, start, end, n_bins)
            coverages[key] = _smooth_coverage(raw, smooth_sigma) if smooth_sigma > 0 else raw
        row_max = coverages[key].max()
        if row["geno_key"] in break_set:
            break_true_max[row["geno_key"]] = max(
                break_true_max.get(row["geno_key"], 0.0), row_max)
        else:
            global_max = max(global_max, row_max)

    ymax_shared: Optional[float]
    if shared_ymax is not None:
        ymax_shared = shared_ymax
    elif not per_track_ymax:
        ymax_shared = _clean_ymax(global_max * 1.05)
    else:
        ymax_shared = None

    # ── figure geometry ────────────────────────────────────────────────────
    has_rna   = rnaseq_data is not None
    rna_frac  = 0.17 if has_rna else 0.0
    ann_ratio = 1.2   # annotation row relative to a merged-track height

    # Compute total figure height — include optional annotation sub-rows
    total_h = sum(rep_track_height if r["kind"] == "replicate" else track_height
                  for r in rows)
    total_h += ann_ratio * track_height   # RefGene row
    if named_peaks:
        total_h += peaks_row_height * len(named_peaks)
    if macs2_peaks is not None:
        total_h += peaks_row_height
    if tf_sites:
        n_tfs = max(len(tf_sites), 1)
        total_h += tf_row_height * n_tfs
    total_h += 0.55  # title headroom

    fig = plt.figure(figsize=(fig_width, total_h))
    fig.patch.set_facecolor("white")

    LEFT  = 0.235
    RIGHT = 0.97
    TOP   = 1.0 - 0.40 / total_h
    BOT   = 0.0

    n_rows = len(rows) + 1   # ATAC rows + RefGene row
    height_ratios = [rep_track_height if r["kind"] == "replicate" else track_height
                     for r in rows] + [ann_ratio * track_height]

    # Append optional annotation sub-rows
    if named_peaks:
        n_rows += len(named_peaks)
        height_ratios.extend([peaks_row_height] * len(named_peaks))
    if macs2_peaks is not None:
        n_rows += 1
        height_ratios.append(peaks_row_height)
    if tf_sites:
        n_tfs = max(len(tf_sites), 1)
        n_rows += 1
        height_ratios.append(tf_row_height * n_tfs)

    n_cols       = 2 if has_rna else 1
    width_ratios = [1.0 - rna_frac, rna_frac] if has_rna else [1.0]

    gs = GridSpec(
        n_rows, n_cols,
        figure=fig,
        height_ratios=height_ratios,
        width_ratios=width_ratios,
        hspace=0.0,
        wspace=0.03,
        left=LEFT, right=RIGHT,
        top=TOP,   bottom=BOT,
    )

    # ── draw rows ─────────────────────────────────────────────────────────
    axes_atac = []
    for i, row in enumerate(rows):
        color   = row["color"]
        arr     = coverages[row["bw_path"]]
        is_rep  = row["kind"] == "replicate"

        ymax_i = ymax_shared if ymax_shared is not None else _clean_ymax(arr.max() * 1.05)
        alpha  = rep_fill_alpha if is_rep else fill_alpha
        lw     = 0.25 if is_rep else 0.40

        # A broken-axis outlier gets a real two-segment split: a compressed
        # upper segment (its own linear scale from ymax_i up to its true
        # max, so its shape above the shared ceiling stays visible) stacked
        # on the normal lower segment (0 to ymax_i, on the same scale as
        # every other track). Only applies to merged rows whose real max
        # actually exceeds the ceiling -- a genotype being listed in
        # break_genotypes doesn't guarantee it's the outlier at every
        # locus (see the tll case), so this only fires when there's
        # something real to show above the break.
        true_max = None
        is_broken_row = False
        if (not is_rep) and row["geno_key"] in break_set:
            true_max = break_true_max.get(row["geno_key"], arr.max())
            is_broken_row = true_max > ymax_i * 1.001

        if is_broken_row:
            upper_ceiling = _clean_ymax(true_max * 1.08)
            inner = GridSpecFromSubplotSpec(
                2, 1, subplot_spec=gs[i, 0],
                height_ratios=[0.32, 0.68], hspace=0.14)
            ax_upper = fig.add_subplot(inner[0, 0])
            ax = fig.add_subplot(inner[1, 0])   # "ax" stays the row's main/lower axes
        else:
            ax = fig.add_subplot(gs[i, 0])
            ax_upper = None

        axes_atac.append(ax)

        ax.fill_between(x, arr, color=color, alpha=alpha, linewidth=0)
        ax.plot(x, arr, color=row["outline_color"], linewidth=lw, alpha=0.9)

        ax.set_xlim(start, end)
        ax.set_ylim(0, ymax_i)

        # Shared vertical gridlines (same positions used by every track below,
        # down through peaks/TF rows and the labeled ruler on the annotation
        # row) so a peak here can be traced straight down to what's binding it.
        _draw_tick_gridlines(ax, shared_ticks, zorder=0)

        # ── Y-AXIS (ADD HERE) ──
        ax.spines["left"].set_visible(True)
        ax.spines["left"].set_color("#444444")
        ax.spines["left"].set_linewidth(0.8)

        ax.yaxis.set_ticks_position("left")
        ax.yaxis.set_tick_params(
            length=2.5,
            width=0.6,
            colors="#444444",
            labelsize=5
        )

        ax.set_yticks(np.linspace(0, ymax_i, 3))
        ax.set_yticklabels([])
        # ───────────────────────

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)

        ax.spines["left"].set_visible(True)
        ax.spines["left"].set_color("#444444")
        ax.spines["left"].set_linewidth(0.8)

        ax.axhline(0, color="#dddddd", linewidth=0.4, zorder=0)
        ax.tick_params(axis="x",
               bottom=False, top=False,
               labelbottom=False)

        ax.tick_params(axis="y",
                    left=True, right=False,
                    length=2.5,
                    width=0.6,
                    colors="#444444",
                    labelleft=False)

        # CPM max label — only on merged rows to avoid clutter
        if not is_rep:
            ax.text(0.004, 0.93, f"{ymax_i:.0f}",
                    transform=ax.transAxes,
                    fontsize=5.5, color="#999999", va="top", ha="left",
                    clip_on=False)

        # Broken-axis outlier: draw the compressed upper segment showing
        # this track's real shape above the shared ceiling, with break
        # marks at the junction and the true peak value labeled.
        if is_broken_row:
            ax_upper.fill_between(x, arr, color=color, alpha=alpha, linewidth=0)
            ax_upper.plot(x, arr, color=row["outline_color"], linewidth=lw, alpha=0.9)
            ax_upper.set_xlim(start, end)
            ax_upper.set_ylim(ymax_i, upper_ceiling)

            _draw_tick_gridlines(ax_upper, shared_ticks, zorder=0)

            ax_upper.spines["top"].set_visible(False)
            ax_upper.spines["right"].set_visible(False)
            ax_upper.spines["bottom"].set_visible(False)
            ax_upper.spines["left"].set_visible(True)
            ax_upper.spines["left"].set_color("#444444")
            ax_upper.spines["left"].set_linewidth(0.8)

            ax_upper.tick_params(axis="x", bottom=False, top=False, labelbottom=False)
            ax_upper.tick_params(axis="y", left=True, right=False,
                                 length=2.5, width=0.6, colors="#444444",
                                 labelleft=False)
            ax_upper.set_yticks([ymax_i, upper_ceiling])

            # Break marks at the junction: bottom of the compressed upper
            # segment and top of the normal lower segment, so the seam
            # reads as a deliberate break rather than a rendering glitch.
            _draw_axis_break(ax_upper, color=row["outline_color"], at="bottom")
            _draw_axis_break(ax, color=row["outline_color"], at="top")

            # Precise numeric readout -- the compressed scale makes the
            # upper segment's own tick spacing coarse, so this keeps the
            # exact value legible regardless.
            ax_upper.text(0.996, 0.90, f"peak = {true_max:.1f}",
                          transform=ax_upper.transAxes,
                          fontsize=5.5, color=break_label_color, fontweight="bold",
                          va="top", ha="right", clip_on=False)
            ax_upper.text(0.004, 0.90, f"{upper_ceiling:.0f}",
                          transform=ax_upper.transAxes,
                          fontsize=5.5, color="#999999", va="top", ha="left",
                          clip_on=False)

        # Genotype label: merged rows get full bold label; rep rows get "rep N" in grey
        if is_rep:
            ax.text(-0.008, 0.50, row["label"],
                    transform=ax.transAxes,
                    fontsize=6.0, color="#aaaaaa",
                    ha="right", va="center", clip_on=False)
        else:
            ax.text(-0.008, 0.50, row["label"],
                    transform=ax.transAxes,
                    fontsize=7.2, color=color, fontweight="bold",
                    ha="right", va="center", clip_on=False)

        # Timepoint separator — dashed line above the first row of a new group
        if row.get("sep_above", False):
            ax.axhline(ymax_i * 0.998, color="#aaaaaa", linewidth=0.9,
                       linestyle=(0, (5, 3)), alpha=0.85, zorder=6)

    # Mirror labeled ticks on top of the topmost track — the labeled ruler at
    # the bottom of the figure is easy to lose track of once TF/peaks rows
    # are stacked underneath it, so a light-weight copy up top gives an
    # at-a-glance reference right next to the tallest peaks.
    if axes_atac:
        ax_top = axes_atac[0]
        ax_top.spines["top"].set_visible(True)
        ax_top.spines["top"].set_color("#cccccc")
        ax_top.spines["top"].set_linewidth(0.5)
        ax_top.set_xticks(shared_ticks)
        ax_top.set_xticklabels(shared_tick_labels)
        ax_top.tick_params(axis="x", top=True, bottom=False,
                           labeltop=True, labelbottom=False,
                           labelsize=6, colors="#555555",
                           length=2, width=0.5, pad=1)

    # ── gene annotation + MACS2 peaks + TF binding rows ──────────────────
    ann_row_start = len(rows)
    _draw_combined_annotation_and_tracks(
        fig=fig, gs=gs, row_start=ann_row_start,
        gene=gene, chrom=chrom, start=start, end=end,
        annotation_features=annotation_features,
        macs2_peaks=macs2_peaks,
        named_peaks=named_peaks,
        tf_sites=tf_sites,
        tf_colors=tf_colors,
        peaks_color=peaks_color,
        peaks_score_threshold=peaks_score_threshold,
        ticks=shared_ticks,
        tick_labels=shared_tick_labels,
    )

    # ── RNA panel ─────────────────────────────────────────────────────────
    if has_rna:
        all_vals = [v for vals in rnaseq_data.values() for v in vals]
        rna_ymax = _clean_ymax(max(all_vals) * 1.05) if all_vals else 1.0

        # Map geno_key → RNA axes index (one RNA violin per merged row)
        rna_row_idx = [i for i, r in enumerate(rows) if r["kind"] == "merged"]

        for idx in rna_row_idx:
            row   = rows[idx]
            geno  = row["geno_key"]
            color = row["color"]
            vals  = rnaseq_data.get(geno, [0.0])

            ax_r = fig.add_subplot(gs[idx, 1])
            if len(vals) >= 4:
                parts = ax_r.violinplot([vals], positions=[0], widths=0.65,
                                        showmedians=True, showextrema=False)
                for pc in parts["bodies"]:
                    pc.set_facecolor(color); pc.set_alpha(0.72); pc.set_edgecolor("none")
                parts["cmedians"].set_color(color)
                parts["cmedians"].set_linewidth(1.4)
            else:
                ax_r.scatter([0] * len(vals), vals, color=color, s=8, alpha=0.85, zorder=3)
            ax_r.set_xlim(-0.9, 0.9)
            ax_r.set_ylim(0, rna_ymax)
            ax_r.axis("off")

        # Replicate rows in RNA column — blank
        for i, row in enumerate(rows):
            if row["kind"] == "replicate":
                ax_r = fig.add_subplot(gs[i, 1])
                ax_r.axis("off")

        # Blank RNA cells for annotation sub-rows
        for extra_i in range(ann_row_start, n_rows):
            fig.add_subplot(gs[extra_i, 1]).axis("off")

        # RNA title on very first RNA axes
        if rna_row_idx:
            axes_atac[rna_row_idx[0]].set_title(
                "GEX", fontsize=6.5, color="#666666", pad=2)

        # Annotation row — blank RNA cell (legacy; loop above already covers extra rows)

    # ── TF colour legend ───────────────────────────────────────────────────
    if tf_sites:
        _tf_legend(fig, tf_sites, tf_colors=tf_colors,
                   x=LEFT, y=max(BOT - 0.06, -0.04))

    # ── figure title ───────────────────────────────────────────────────────
    title = title_override or f"{gene}   {chrom}:{start:,}–{end:,}"
    atac_right = RIGHT - rna_frac * (RIGHT - LEFT)
    fig.text((LEFT + atac_right) / 2, TOP + 0.28 / total_h, title,
             ha="center", va="bottom",
             fontsize=9.5, color="#1a1a1a", fontweight="medium")

    plt.savefig(out_path, dpi=200, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
# BATCH MODE
# Generates three PDFs per gene: all / ventralized / non_ventralized
# ─────────────────────────────────────────────────────────────────────────────

def run_batch_from_bed(
    bed_path: str,
    bigwig_dir: str,
    out_dir: str,
    modes: List[str],
    annotation_bed: Optional[str] = None,
    annotation_gtf: Optional[str] = None,
    rnaseq_dir: Optional[str] = None,
    stage_filter: str = "all",
    **kwargs,
):
    """
    BED file (tab-separated): chrom  start  end  gene_name  [score  strand]
    For each gene × each requested mode, one PDF is written to:
      {out_dir}/{mode}/{gene}_coverage_{mode}.pdf
    """
    for mode in modes:
        Path(os.path.join(out_dir, mode)).mkdir(parents=True, exist_ok=True)

    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols  = line.strip().split("\t")
            chrom = cols[0]
            start = int(cols[1])
            end   = int(cols[2])
            gene  = cols[3] if len(cols) > 3 else f"{chrom}_{start}"

            rna_data  = None
            if rnaseq_dir:
                rna_data = load_rnaseq_tsv(
                    os.path.join(rnaseq_dir, f"{gene}_expression.tsv"))

            ann_feats = load_annotation_bed(annotation_bed, chrom, start, end)

            for mode in modes:
                genotypes = apply_stage_filter(MODE_GENOTYPES[mode], stage_filter)
                if not genotypes:
                    print(f"  [{mode}] {gene}: SKIPPED -- --stage_filter "
                          f"{stage_filter} left no genotypes for this mode",
                          file=sys.stderr)
                    continue
                out_path  = os.path.join(out_dir, mode,
                                          f"{gene}_coverage_{mode}.pdf")
                print(f"  [{mode}] {gene} {chrom}:{start}-{end}")
                plot_tracks(
                    gene=gene, chrom=chrom, start=start, end=end,
                    bigwig_dir=bigwig_dir,
                    genotypes=genotypes,
                    out_path=out_path,
                    rnaseq_data=rna_data,
                    annotation_features=ann_feats if ann_feats else None,
                    annotation_gtf=annotation_gtf,
                    **kwargs,
                )


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def build_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    sg = p.add_argument_group("Single locus")
    sg.add_argument("--gene",  help="Gene name (title + output filename)")
    sg.add_argument("--chrom", help="Chromosome  e.g. chr2L")
    sg.add_argument("--start", help="Start coord (underscores OK: 15_450_000)")
    sg.add_argument("--end",   help="End coord")

    bg = p.add_argument_group("Batch mode")
    bg.add_argument("--gene_bed", help="BED file of loci (chrom start end gene [score strand])")
    bg.add_argument("--out_dir",  default="./coverage_plots",
                    help="Output directory for batch mode  [./coverage_plots]")

    dg = p.add_argument_group("Data")
    dg.add_argument("--bigwig_dir", required=True,
                    help="Root BigWig directory (contains both merged/ and per-sample files)")
    _MODE_CHOICES = ["all", "ventralized", "non_ventralized",
                      "all_timepoints", "ventralized_timepoints",
                      "non_ventralized_timepoints", "CVM", "CVM_solo"]
    dg.add_argument("--mode", choices=_MODE_CHOICES,
                    default="all",
                    help="Genotype subset preset  [all]. Use 'all_timepoints' "
                         "(or the _timepoints variants) to include every "
                         "confirmed timepoint per genotype, not "
                         "just the nc14b/nc14late anchors.")
    dg.add_argument("--modes", nargs="+",
                    choices=_MODE_CHOICES,
                    default=["all", "ventralized", "non_ventralized"],
                    help="(Batch only) list of modes to generate  [all three]")
    dg.add_argument("--genotypes", nargs="+", default=None,
                    help="Manual genotype key list (overrides --mode)")
    dg.add_argument("--stage_filter", choices=["all", "nc14b", "nc14d"],
                    default="all",
                    help="Restrict the genotype/timepoint tracks in this plot to "
                         "one developmental stage, applied AFTER --mode/--genotypes "
                         "resolves the genotype list. 'nc14b' keeps only tracks "
                         "whose MERGED_CONFIG tp_group is 'nc14b'. 'nc14d' keeps "
                         "tracks whose tp_group is 'nc14late', 'nc14d', OR 'late' -- "
                         "all three are the SAME stage, just named differently "
                         "depending on which genotype/labmate convention set them "
                         "up (see MERGED_CONFIG), so all three are treated as one "
                         "bucket here. 'all' (default) applies no filtering.")
    dg.add_argument("--no_replicates", action="store_true",
                    help="Show merged tracks only (no per-replicate sub-tracks)")
    dg.add_argument("--bg_sub", action="store_true",
                    help="Use background-subtracted BigWigs (_bgSubCPM) instead of "
                         "raw CPM. Applies to all tracks by swapping stem suffix. "
                         "Particularly intended for CVM / CVM_solo modes.")
    dg.add_argument("--rnaseq_tsv",
                    help="TSV: genotype<TAB>CPM_value (single locus)")
    dg.add_argument("--rnaseq_dir",
                    help="(Batch) directory of per-gene TSVs: {gene}_expression.tsv")
    dg.add_argument("--annotation_bed",
                    help="BED6 or BED12 for gene-body annotation row")
    dg.add_argument("--annotation_gtf",
                    help="GTF or GFF3 file for gene-body annotation (FlyBase GFF3 recommended)")

    ag = p.add_argument_group("MACS2 peaks & TF binding sites")
    ag.add_argument("--macs2_peaks", nargs="+", default=None,
                    metavar="FILE",
                    help="One or more MACS2 narrowPeak / BED files combined into a single "
                         "grey 'MACS2 peaks' row below the gene model.")
    ag.add_argument("--macs2_named_peaks", nargs="+", default=None,
                    metavar="LABEL=COLOR=FILE",
                    help="Per-genotype MACS2 peaks — one colored row per entry. "
                         "Format:  'Label=#hexcolor=/path/to.narrowPeak'\n"
                         "Color can be omitted (auto-assigned): 'Label=/path/to.narrowPeak'\n"
                         "Example:\n"
                         "  'BOTv=#2ca02c=/data/BOTv_peaks.narrowPeak'\n"
                         "  'BOTCv=#e377c2=/data/BOTCv_peaks.narrowPeak'\n"
                         "Rows appear in input order between RefGene and any --macs2_peaks row.")
    ag.add_argument("--peaks_color", default="#333333",
                    help="Fill colour for MACS2 peak rectangles  [#333333]")
    ag.add_argument("--peaks_score_threshold", type=float, default=0.0,
                    help="Minimum MACS2 score (col 5, 0-1000) to display  [0]")
    ag.add_argument("--peaks_row_height", type=float, default=0.30,
                    help="Height in inches of the MACS2 peaks row  [0.30]")
    ag.add_argument("--tf_sites", nargs="+", default=None,
                    metavar="FILE_OR_NAME=FILE",
                    help="TF binding-site BED files. Two formats accepted:\n"
                         "  (a) Single multi-TF BED where col 3 (--tf_name_col) "
                             "encodes the TF name.\n"
                         "  (b) Per-TF files as  'TFname=path'  pairs, "
                             "e.g. Doc=/data/Doc.bed Sna=/data/Sna.bed\n"
                         "A colour-coded tick row is drawn below MACS2 peaks.")
    ag.add_argument("--tf_name_col", type=int, default=3,
                    help="0-based BED column index holding the TF name (multi-TF BED)  [3]")
    ag.add_argument("--tf_colors", nargs="+", default=None,
                    metavar="TF=HEX",
                    help="Override per-TF colours, e.g.  Doc=#d62728 Sna=#1f77b4")
    ag.add_argument("--tf_row_height", type=float, default=0.22,
                    help="Height in inches per TF row  [0.22]")

    cg = p.add_argument_group("Cosmetics")
    cg.add_argument("--ymax",           type=float, default=None,
                    help="Shared y-axis ceiling (default: auto)")
    cg.add_argument("--per_track_ymax", action="store_true",
                    help="Independent y-axis per track")
    cg.add_argument("--break_genotypes", nargs="+", default=None,
                    help="Genotype key(s) (e.g. BOTR) whose MERGED track "
                         "should be clipped to the shared ceiling computed "
                         "from every OTHER track, instead of that outlier "
                         "setting the ceiling for the whole panel. The "
                         "clipped track gets a '//' break mark and a "
                         "numeric label with its true peak value. Overrides "
                         "--per_track_ymax if both are given.")
    cg.add_argument("--break_label_color", default="#cc0000",
                    help="Text color for the true-peak-value label on "
                         "broken-axis tracks  [#cc0000]")
    cg.add_argument("--smooth",         type=float, default=0.0,
                    metavar="SIGMA",
                    help="Gaussian smoothing sigma in bins applied after fetching "
                         "coverage (0 = off).  Useful for bgSub tracks.  "
                         "Typical values: 3 (light), 5 (default for bgSub), "
                         "10 (browser-style), 25 (heavy).  "
                         "Requires scipy for best results; falls back to numpy box-car.")
    cg.add_argument("--n_bins",         type=int,   default=1000)
    cg.add_argument("--fig_width",      type=float, default=8.5)
    cg.add_argument("--track_height",   type=float, default=0.72,
                    help="Merged track height in inches  [0.72]")
    cg.add_argument("--rep_track_height", type=float, default=0.52,
                    help="Replicate track height in inches  [0.52]")
    cg.add_argument("--fill_alpha",     type=float, default=0.80)
    cg.add_argument("--rep_fill_alpha", type=float, default=0.55)

    og = p.add_argument_group("Output")
    og.add_argument("--out", default=None, help="Output file (single locus)")

    return p


def main():
    args = build_parser().parse_args()

    # Parse TF colour overrides: ["Doc=#d62728", "Sna=#1f77b4"] → dict
    tf_colors_dict: Optional[Dict[str, str]] = None
    if args.tf_colors:
        tf_colors_dict = {}
        for item in args.tf_colors:
            if "=" in item:
                k, _, v = item.partition("=")
                tf_colors_dict[k.strip()] = v.strip()

    common = dict(
        bigwig_dir             = args.bigwig_dir,
        shared_ymax            = args.ymax,
        per_track_ymax         = args.per_track_ymax,
        break_genotypes        = args.break_genotypes,
        break_label_color      = args.break_label_color,
        show_replicates        = not args.no_replicates,
        bg_sub                 = args.bg_sub,
        smooth_sigma           = args.smooth,
        n_bins                 = args.n_bins,
        fig_width              = args.fig_width,
        track_height           = args.track_height,
        rep_track_height       = args.rep_track_height,
        fill_alpha             = args.fill_alpha,
        rep_fill_alpha         = args.rep_fill_alpha,
        annotation_gtf         = args.annotation_gtf,
        # MACS2 / TF
        macs2_peak_files       = args.macs2_peaks,
        macs2_named_peaks      = args.macs2_named_peaks,
        peaks_color            = args.peaks_color,
        peaks_score_threshold  = args.peaks_score_threshold,
        peaks_row_height       = args.peaks_row_height,
        tf_files               = args.tf_sites,
        tf_colors              = tf_colors_dict,
        tf_name_col            = args.tf_name_col,
        tf_row_height          = args.tf_row_height,
    )

    # ── batch mode ─────────────────────────────────────────────────────────
    if args.gene_bed:
        run_batch_from_bed(
            bed_path       = args.gene_bed,
            out_dir        = args.out_dir,
            modes          = args.modes,
            annotation_bed = args.annotation_bed,
            rnaseq_dir     = args.rnaseq_dir,
            stage_filter   = args.stage_filter,
            **common,
        )
        return

    # ── single locus ───────────────────────────────────────────────────────
    if not all([args.gene, args.chrom, args.start, args.end]):
        print("ERROR: provide --gene --chrom --start --end  (or --gene_bed)",
              file=sys.stderr)
        sys.exit(1)

    chrom = args.chrom
    start = int(args.start.replace("_", ""))
    end   = int(args.end.replace("_", ""))

    if args.genotypes:
        genotypes = args.genotypes
    else:
        genotypes = MODE_GENOTYPES[args.mode]

    genotypes = apply_stage_filter(genotypes, args.stage_filter)
    if not genotypes:
        print(f"  WARNING: --stage_filter {args.stage_filter} left no genotypes "
              f"to plot for this --mode/--genotypes selection -- skipping.",
              file=sys.stderr)
        return

    rna_data  = load_rnaseq_tsv(args.rnaseq_tsv)
    ann_feats = load_annotation_bed(args.annotation_bed, chrom, start, end)

    out_path = args.out or f"{args.gene}_coverage_{args.mode}.pdf"

    plot_tracks(
        gene                = args.gene,
        chrom               = chrom,
        start               = start,
        end                 = end,
        genotypes           = genotypes,
        out_path            = out_path,
        rnaseq_data         = rna_data,
        annotation_features = ann_feats if ann_feats else None,
        **common,
    )


if __name__ == "__main__":
    main()
