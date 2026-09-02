"""
Config_DataLoader.py  —  Shared data loading for all downstream plotting scripts.

Reads from the combined matrices produced by the upstream pipeline
(step5b → step5c → step6, with optional RUVg batch correction) and
provides a consistent interface across all plotting scripts.

RUVg batch correction
---------------------
The pipeline has a confirmed set1 vs set2 batch effect (PC1 = 27.9% before
correction, reduced to 16.5% after RUVg k=1). All downstream analyses that
compare across the two batches SHOULD use RUVg-corrected data. This loader
defaults to use_ruvg=True and will warn loudly if the corrected files are
absent.

File naming convention
----------------------
Uncorrected (prefix = e.g. "all_genotypes" or "all_genotypes_plus_set3"):
  results/combined/<prefix>_counts_matrix_qc.tsv
  results/combined/<prefix>_qc_tpm_matrix.tsv
  results/combined/<prefix>_qc_tpm_quantile_matrix.tsv

RUVg-corrected — everything below derived from <prefix>_ruvg:
  results/combined/<prefix>_ruvg_counts_matrix_qc.tsv
  results/combined/<prefix>_ruvg_qc_tpm_matrix.tsv
  results/combined/<prefix>_ruvg_qc_tpm_quantile_matrix.tsv

To generate RUVg files (prefix substituted for whatever --combined you used):
  python make_ruvg_metadata.py --combined <prefix>
  Rscript ruvg_nextgen.R --counts results/combined/<prefix>_counts_matrix_qc.tsv \\
      --metadata results/combined/sample_metadata.tsv \\
      --outdir results/combined/ruvg_nextgen_<prefix>
  cp results/combined/ruvg_nextgen_<prefix>/corrected_counts.tsv \\
     results/combined/<prefix>_ruvg_counts_matrix_qc.tsv
  python step6_tpm.py --combined <prefix>_ruvg --qc --quantile

Note: ruvg_covariate_design_v3.R is NOT the correction step — it's a
downstream DESeq2/limma script that assumes RUVg has already run and reads
its W_factors.tsv as a covariate. Run ruvg_nextgen.R first.

Column naming
-------------
Pipeline columns are prefixed with the sample name:
  WT_seRNAseq_4embryos_set1__WT_embryo_bc2
  28232_28233_28271_28272_28273_S1__runt_Dm_nc14b_bc07

After loading, the sample prefix is stripped:
  WT_embryo_bc2
  runt_Dm_nc14b_bc07

Usage
-----
    import Config_DataLoader as DL

    # Default: RUVg-corrected, QC-filtered, quantile-normalised
    data = DL.load()

    # Explicitly request uncorrected data (e.g. for QC plots)
    data = DL.load(use_ruvg=False)

    # data.counts    — count matrix
    # data.tpm       — raw TPM (for absolute expression values)
    # data.tpm_qn    — quantile-normalised TPM (for Z-scores, PCA, clustering)
    # data.group_of  — {col: group_label}
    # data.ruvg      — True/False (was RUVg correction applied?)

    by_group = DL.embryos_by_group(data)   # {group: [col, ...]}
    cols     = DL.ordered_cols(data)       # all cols in SC.GROUP_ORDER order
"""

import os
import sys
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

try:
    import Config_SampleMetadata as SC
    _HAS_SC = True
except ImportError:
    _HAS_SC = False
    print("[data_loader] WARNING: Config_SampleMetadata.py not found — "
          "group assignment will be limited")

# ── Defaults ─────────────────────────────────────────────────────────────────

RESULTS_ROOT   = "results"
COMBINED_DIR   = os.path.join(RESULTS_ROOT, "combined")
DEFAULT_PREFIX = "all_genotypes"


@dataclass
class CombinedData:
    """Container for all loaded matrices and metadata."""
    counts        : pd.DataFrame
    tpm           : pd.DataFrame
    tpm_qn        : pd.DataFrame | None
    embryo_cols   : list[str]
    col_to_sample : dict[str, str]
    group_of      : dict[str, str]
    qc_metrics    : pd.DataFrame | None = None
    ruvg          : bool = False
    _raw_col_map  : dict[str, str] = field(default_factory=dict, repr=False)


# ── Column name helpers ───────────────────────────────────────────────────────

def _strip_sample_prefix(raw_col: str) -> str:
    """
    Strip the sample prefix from a pipeline column name.

    WT_seRNAseq_4embryos_set1__WT_embryo_bc2     →  WT_embryo_bc2
    Mutant_seRNAseq_19embryos_set2__embryo_bc1   →  Mutant_embryo_bc1
    28232_..._S1__runt_Dm_nc14b_bc07             →  runt_Dm_nc14b_bc07
    """
    if "__" not in raw_col:
        return raw_col

    sample_part, embryo_part = raw_col.split("__", 1)

    # Set1 WT
    if "WT_seRNAseq" in sample_part or sample_part.startswith("WT"):
        return embryo_part if embryo_part.startswith("WT_") else f"WT_{embryo_part}"

    # Set1 Mutant — embryo IDs were renamed to canonical form
    # (FoxL1_BOTv_embryo_bcN / runt_embryo_bcN) by rename_matrix_columns.py.
    # Pass through as-is. Legacy "Mutant_embryo_bcN" names (pre-rename) are
    # preserved unchanged for backward compatibility, but should no longer
    # appear in matrices after the canonical rename pass.
    if "Mutant_seRNAseq" in sample_part:
        return embryo_part

    # Set2 — embryo_part already encodes genotype + stage
    return embryo_part


def _sample_of(raw_col: str) -> str:
    if "__" in raw_col:
        return raw_col.split("__")[0]
    return "unknown"


def _assign_group(clean_col: str) -> str:
    """
    Assign a biological group label to a clean column name.

    Priority:
      1. VALIDATED_BY_GROUP explicit lookup  — handles Mutant_ ambiguity
      2. GENOTYPE_PREFIX_TO_GROUP prefix match — set2 embryos
      3. WT prefix check
      4. "unknown" — dropped downstream
    """
    if not _HAS_SC:
        if "WT_embryo" in clean_col and not clean_col.startswith("Mutant"):
            return "WT"
        return "unknown"

    # 1. Explicit validated list (set1 Mutant embryos and any set2 in lists)
    if hasattr(SC, "VALIDATED_BY_GROUP"):
        bare = clean_col.replace("Mutant_", "")
        for grp, validated in SC.VALIDATED_BY_GROUP.items():
            if clean_col in validated or bare in validated:
                return grp

    # 2. Prefix match for set2 embryos and set1 WT
    if hasattr(SC, "GENOTYPE_PREFIX_TO_GROUP"):
        # Sort by prefix length descending so longer prefixes match first
        for prefix, grp in sorted(SC.GENOTYPE_PREFIX_TO_GROUP.items(),
                                   key=lambda x: -len(x[0])):
            if clean_col.startswith(prefix):
                return grp

    # 3. WT fallback
    if clean_col.startswith("WT_embryo") or clean_col.startswith("WT_WT"):
        return "WT"

    # 4. Set1 Mutant not in any validated list
    if clean_col.startswith("Mutant_"):
        return "unknown"

    return "unknown"


def _read_matrix(path: str, label: str, required: bool = True) -> pd.DataFrame | None:
    if not os.path.exists(path):
        if required:
            print(f"[data_loader] {label} not found: {path}")
        return None
    df = pd.read_csv(path, sep="\t", index_col=0)
    if "GeneSymbol" in df.columns:
        df = df.drop(columns="GeneSymbol")
    return df


# ── Main load function ────────────────────────────────────────────────────────

def load(prefix       : str  = DEFAULT_PREFIX,
         combined_dir : str  = COMBINED_DIR,
         use_qc       : bool = True,
         use_quantile : bool = True,
         use_ruvg     : bool = True,
         verbose      : bool = True) -> CombinedData:
    """
    Load combined matrices from the pipeline output.

    Parameters
    ----------
    prefix : str
        Base combined matrix prefix (default: "all_genotypes").
    combined_dir : str
        Directory containing combined outputs.
    use_qc : bool
        Load QC-filtered matrix. Default True.
    use_quantile : bool
        Also load quantile-normalised TPM. Default True.
    use_ruvg : bool
        Prefer RUVg batch-corrected matrices when available. Default True.
        Falls back to uncorrected with a warning if RUVg files are absent.
    verbose : bool
        Print loading status.
    """
    # ── Decide which prefix to use ───────────────────────────────────────────
    ruvg_applied = False

    if use_ruvg:
        ruvg_prefix  = f"{prefix}_ruvg"
        _ruvg_counts = os.path.join(
            combined_dir, f"{ruvg_prefix}_counts_matrix_qc.tsv")
        if os.path.exists(_ruvg_counts):
            active_prefix = ruvg_prefix
            ruvg_applied  = True
        else:
            if verbose:
                print("[data_loader] ⚠  RUVg-corrected matrix not found — "
                      "falling back to uncorrected data.")
                print("[data_loader]    A batch effect may be present between "
                      "sequencing runs.")
                print("[data_loader]    Cross-batch comparisons may be confounded.")
                print("[data_loader]    To generate corrected data:")
                print(f"[data_loader]      python make_ruvg_metadata.py --combined {prefix}")
                print(f"[data_loader]      Rscript ruvg_nextgen.R "
                      f"--counts {combined_dir}/{prefix}_counts_matrix_qc.tsv \\")
                print(f"[data_loader]          --metadata {combined_dir}/sample_metadata.tsv "
                      f"--outdir {combined_dir}/ruvg_nextgen_{prefix}")
                print(f"[data_loader]      cp {combined_dir}/ruvg_nextgen_{prefix}/"
                      f"corrected_counts.tsv \\")
                print(f"[data_loader]         {_ruvg_counts}")
                print(f"[data_loader]      python step6_tpm.py "
                      f"--combined {ruvg_prefix} --qc --quantile")
            active_prefix = prefix
    else:
        active_prefix = prefix

    # ── Build file paths ─────────────────────────────────────────────────────
    qc_suffix    = "_qc" if use_qc else ""
    label        = f"{active_prefix}{qc_suffix}"
    counts_path  = os.path.join(combined_dir,
                                f"{active_prefix}_counts_matrix_qc.tsv" if use_qc
                                else f"{active_prefix}_counts_matrix.tsv")
    tpm_path     = os.path.join(combined_dir, f"{label}_tpm_matrix.tsv")
    qn_path      = os.path.join(combined_dir, f"{label}_tpm_quantile_matrix.tsv")
    metrics_path = os.path.join(combined_dir, f"{prefix}_qc_metrics.tsv")

    # Fallback: no _qc suffix version
    if use_qc and not os.path.exists(counts_path):
        alt = os.path.join(combined_dir, f"{active_prefix}_counts_matrix.tsv")
        if os.path.exists(alt):
            if verbose:
                print(f"[data_loader] QC matrix not found; using unfiltered.")
            counts_path = alt
            label    = active_prefix
            tpm_path = os.path.join(combined_dir, f"{label}_tpm_matrix.tsv")
            qn_path  = os.path.join(combined_dir, f"{label}_tpm_quantile_matrix.tsv")

    # ── Load matrices ────────────────────────────────────────────────────────
    counts = _read_matrix(counts_path, "Counts matrix")
    if counts is None:
        sys.exit(f"[data_loader] Cannot proceed without counts matrix.\n"
                 f"  Expected: {counts_path}")

    tpm = _read_matrix(tpm_path, "TPM matrix")
    if tpm is None:
        sys.exit(f"[data_loader] Cannot proceed without TPM matrix.\n"
                 f"  Expected: {tpm_path}\n"
                 f"  Run: python step6_tpm.py --combined {active_prefix} --qc first.")

    tpm_qn = None
    if use_quantile:
        tpm_qn = _read_matrix(qn_path, "Quantile-normalized TPM", required=False)
        if tpm_qn is None and verbose:
            print("[data_loader] Quantile-normalized TPM not available — "
                  "using raw TPM for all analyses.")

    qc_metrics = None
    if os.path.exists(metrics_path):
        qc_metrics = pd.read_csv(metrics_path, sep="\t", index_col=0)

    # ── Rename columns ───────────────────────────────────────────────────────
    raw_cols     = list(counts.columns)
    raw_to_clean = {rc: _strip_sample_prefix(rc) for rc in raw_cols}
    clean_to_raw = {v: k for k, v in raw_to_clean.items()}

    counts = counts.rename(columns=raw_to_clean)
    tpm    = tpm.rename(columns=raw_to_clean)
    if tpm_qn is not None:
        tpm_qn = tpm_qn.rename(columns=raw_to_clean)

    clean_cols    = list(counts.columns)
    col_to_sample = {_strip_sample_prefix(rc): _sample_of(rc) for rc in raw_cols}
    group_of      = {c: _assign_group(c) for c in clean_cols}

    # ── Manual exclusions ────────────────────────────────────────────────────
    if _HAS_SC and hasattr(SC, "apply_manual_exclude") and hasattr(SC, "MANUAL_EXCLUDE"):
        before     = len(clean_cols)
        clean_cols = SC.apply_manual_exclude(clean_cols, verbose=verbose)
        if len(clean_cols) < before:
            counts = counts[clean_cols]
            tpm    = tpm[clean_cols]
            if tpm_qn is not None:
                tpm_qn = tpm_qn[clean_cols]
            col_to_sample = {c: col_to_sample[c] for c in clean_cols
                             if c in col_to_sample}
            group_of = {c: group_of[c] for c in clean_cols}

    # ── Drop unvalidated ─────────────────────────────────────────────────────
    unknown = [c for c in clean_cols if group_of.get(c) == "unknown"]
    if unknown:
        if verbose:
            print(f"[data_loader] Dropping {len(unknown)} unvalidated "
                  f"embryo(s): {unknown}")
        clean_cols = [c for c in clean_cols if group_of.get(c) != "unknown"]
        counts = counts[clean_cols]
        tpm    = tpm[clean_cols]
        if tpm_qn is not None:
            tpm_qn = tpm_qn[clean_cols]
        col_to_sample = {c: col_to_sample[c] for c in clean_cols
                         if c in col_to_sample}
        group_of = {c: group_of[c] for c in clean_cols}

    # ── Verbose report ───────────────────────────────────────────────────────
    if verbose:
        src  = "RUVg-corrected" if ruvg_applied else "uncorrected"
        qcs  = "QC-filtered" if use_qc else "unfiltered"
        qns  = " + quantile-norm" if tpm_qn is not None else ""
        print(f"[data_loader] Loaded {src} {qcs} combined matrices{qns}")
        print(f"[data_loader]   {len(counts):,} genes × {len(clean_cols)} embryos")
        if _HAS_SC and hasattr(SC, "GROUP_ORDER"):
            parts = [f"{g}: {sum(1 for v in group_of.values() if v == g)}"
                     for g in SC.GROUP_ORDER
                     if any(v == g for v in group_of.values())]
        else:
            from collections import Counter
            parts = [f"{g}: {n}" for g, n in Counter(group_of.values()).items()]
        print(f"[data_loader]   " + "  |  ".join(parts))

    return CombinedData(
        counts=counts, tpm=tpm, tpm_qn=tpm_qn,
        embryo_cols=clean_cols, col_to_sample=col_to_sample,
        group_of=group_of, qc_metrics=qc_metrics,
        ruvg=ruvg_applied, _raw_col_map=clean_to_raw,
    )


# ── Convenience functions ────────────────────────────────────────────────────

def embryos_by_group(data: CombinedData) -> dict[str, list[str]]:
    """Return {group_label: [col, ...]} in GROUP_ORDER order."""
    groups: dict[str, list[str]] = {}
    if _HAS_SC and hasattr(SC, "GROUP_ORDER"):
        for grp in SC.GROUP_ORDER:
            groups[grp] = []
    for col in data.embryo_cols:
        grp = data.group_of.get(col, "unknown")
        groups.setdefault(grp, []).append(col)
    return {g: v for g, v in groups.items() if v}


def ordered_cols(data: CombinedData) -> list[str]:
    """Return embryo columns in SC.GROUP_ORDER order."""
    by_grp = embryos_by_group(data)
    if _HAS_SC and hasattr(SC, "GROUP_ORDER"):
        result = []
        seen   = set()
        for grp in SC.GROUP_ORDER:
            for col in by_grp.get(grp, []):
                if col not in seen:
                    result.append(col)
                    seen.add(col)
        # Any remaining not in GROUP_ORDER
        for col in data.embryo_cols:
            if col not in seen:
                result.append(col)
        return result
    return list(data.embryo_cols)
