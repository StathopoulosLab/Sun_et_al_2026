#!/usr/bin/env python3
"""
step5b_merge_samples_expanded.py  —  Merge set1 validated embryos + set2.

Changes vs previous version
----------------------------
* Exports a richer sample_groups.tsv with group, stage, and set columns —
  consumed by step5c, step6, step7, and all downstream plotting scripts so
  they never have to re-infer these from column names.
* Set1 Mutant columns are canonicalised to FoxL1_BOTv_embryo_bcN /
  runt_embryo_bcN at merge time so data_loader never sees ambiguous bare names.
* Applies SC.apply_manual_exclude() to Set2 at merge time (previously this
  only happened inside data_loader, making the raw combined matrix dirty).
* Added --set2-filter: restrict which Set2 genotype prefixes to include.
* Warns clearly when a validated embryo is absent from the actual count matrix.
* Saves counts as int32 (halves RAM for large matrices).

Output
------
  results/combined/<prefix>_counts_matrix.tsv
  results/combined/<prefix>_sample_groups.tsv   (group / stage / set columns added)
"""

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import config

try:
    import Config_SampleMetadata as SC
    _HAS_SC = True
except ImportError:
    _HAS_SC = False
    print("[step5b] WARNING: Config_SampleMetadata.py not found — "
          "no validated-list filtering will be applied to set1")

# ── Set1 canonical name lookup ──────────────────────────────────────────────
# bare column → canonical name that sample_config / data_loader expect.
_SET1_CANON: dict[str, str] = {
    "WT_embryo_bc2": "WT_embryo_bc2",
    "WT_embryo_bc3": "WT_embryo_bc3",
    "WT_embryo_bc4": "WT_embryo_bc4",
    "embryo_bc1":    "FoxL1_BOTv_embryo_bc1",
    "embryo_bc2":    "FoxL1_BOTv_embryo_bc2",
    "embryo_bc8":    "FoxL1_BOTv_embryo_bc8",
    "FoxL1_BOTv_embryo_bc1": "FoxL1_BOTv_embryo_bc1",
    "FoxL1_BOTv_embryo_bc2": "FoxL1_BOTv_embryo_bc2",
    "FoxL1_BOTv_embryo_bc8": "FoxL1_BOTv_embryo_bc8",
    "embryo_bc13":   "runt_embryo_bc13",
    "embryo_bc14":   "runt_embryo_bc14",
    "embryo_bc18":   "runt_embryo_bc18",
    "runt_embryo_bc13": "runt_embryo_bc13",
    "runt_embryo_bc14": "runt_embryo_bc14",
    "runt_embryo_bc18": "runt_embryo_bc18",
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--wt-sample",     default="WT_seRNAseq_4embryos_set1")
    p.add_argument("--mutant-sample", default="Mutant_seRNAseq_19embryos_set2")
    p.add_argument("--set2-sample",   default="28232_28233_28271_28272_28273_S1")
    p.add_argument("--prefix",        default="all_genotypes")
    p.add_argument("--annotate",      default=None,
                   help="Optional gene_id_to_symbol.tsv for GeneSymbol column")
    p.add_argument("--no-set1",       action="store_true")
    p.add_argument("--no-set2",       action="store_true")
    p.add_argument("--set2-filter",   default=None,
                   help="Comma-separated genotype prefixes to keep from set2, "
                        "e.g. 'runt_Dm,BOTR_Dm'")
    return p.parse_args()


def _load_matrix(sample: str) -> pd.DataFrame:
    path = os.path.join(
        config.OUTPUT_ROOT, sample, "counts", f"{sample}_counts_matrix.tsv")
    if not os.path.exists(path):
        sys.exit(f"[step5b] Count matrix not found: {path}\n"
                 f"         Run step5_merge_counts.py --sample {sample} first.")
    df = pd.read_csv(path, sep="\t", index_col=0)
    if "GeneSymbol" in df.columns:
        df = df.drop(columns="GeneSymbol")
    return df


def _apply_set1_validated_filter(df: pd.DataFrame, sample: str) -> pd.DataFrame:
    """Keep only validated set1 embryos and canonicalise column names.

    Also drops any column whose canonical name resolves to a non-WT group
    when called on the WT sample matrix — these are demux bleed-through embryos
    that belong in the Mutant sample file.
    """
    if not _HAS_SC:
        return df

    all_validated = set(SC.get_all_validated())
    rename_map: dict[str, str] = {}
    keep: list[str] = []
    dropped: list[str] = []
    bleedthrough: list[str] = []

    for col in df.columns:
        canon = _SET1_CANON.get(col)
        if canon is None:
            stripped = col.removeprefix("Mutant_").removeprefix("WT_")
            canon    = _SET1_CANON.get(stripped)

        if canon and canon in all_validated:
            # Extra check: if this is the WT sample file but the canonical name
            # resolves to a non-WT group, it's a demux bleed-through — drop it.
            resolved_group = SC.build_group_map([canon]).get(canon, "unknown")
            is_wt_file = "WT_seRNAseq" in sample or sample.startswith("WT")
            if is_wt_file and resolved_group not in ("WT", "unknown"):
                bleedthrough.append(col)
                print(f"[step5b]   {sample}: dropping bleed-through embryo "
                      f"'{col}' → '{canon}' (group: {resolved_group}) — "
                      f"this embryo belongs in the Mutant sample file")
                continue
            rename_map[col] = canon
            keep.append(col)
        else:
            dropped.append(col)

    if dropped:
        print(f"[step5b]   {sample}: dropping {len(dropped)} non-validated "
              f"embryo(s): {dropped}")

    # Warn about validated embryos absent from THIS sample's matrix
    # (only check embryos that belong to this sample to avoid noise)
    is_wt_file   = "WT_seRNAseq" in sample or sample.startswith("WT")
    is_mut_file  = "Mutant_seRNAseq" in sample
    present_canon = set(rename_map.values())

    check_ids = (SC.VALIDATED_WT if is_wt_file
                 else SC.VALIDATED_RUNT + SC.VALIDATED_FOXL1 if is_mut_file
                 else [])
    for vid in check_ids:
        if vid not in present_canon:
            print(f"[step5b]   WARNING: validated '{vid}' absent from "
                  f"{sample} matrix — check demux / file path")

    if not keep:
        print(f"[step5b]   WARNING: no validated embryos for {sample}. "
              f"Matrix cols: {list(df.columns)[:8]} ...")
        return df

    return df[keep].rename(columns=rename_map)


def _get_metadata(bare: str, sample_name: str, set_id: str) -> dict:
    """Return group/stage/set metadata for one embryo."""
    group = "unknown"
    stage = "nc14b"
    if _HAS_SC:
        gmap  = SC.build_group_map([bare])
        group = gmap.get(bare, "unknown")
        stage = SC.get_stage(bare)
    else:
        for prefix, g in [
            ("WT_embryo", "WT"), ("runt_embryo", "runt"),
            ("runt_Dm", "runt"), ("BOTR_Dm", "BOTR"),
            ("FoxL1_BOTv", "FoxL1_BOTv"),
            ("HLH54F_BOTCv", "HLH54F_BOTCv"), ("B6_BOTC", "B6_BOTC"),
        ]:
            if prefix in bare:
                group = g
                break
        for s in ("gastr", "nc14d", "nc14b"):
            if s in bare:
                stage = s
                break
    return {"sample": sample_name, "group": group, "stage": stage, "set": set_id}


def main():
    args = parse_args()
    combined_dir = os.path.join(config.OUTPUT_ROOT, "combined")
    os.makedirs(combined_dir, exist_ok=True)

    all_dfs:   list[pd.DataFrame] = []
    meta_rows: list[dict]         = []

    set2_prefixes = (
        [p.strip() for p in args.set2_filter.split(",")]
        if args.set2_filter else []
    )

    # ── Set1: WT ──────────────────────────────────────────────────────────────
    if not args.no_set1:
        print(f"[step5b] Loading set1 WT: {args.wt_sample}")
        wt_df = _load_matrix(args.wt_sample)
        wt_df = _apply_set1_validated_filter(wt_df, args.wt_sample)
        renamed: dict[str, str] = {}
        for col in wt_df.columns:
            # Only add WT_ prefix to columns that are genuinely WT embryos.
            # The WT sample matrix occasionally contains a stray FoxL1/runt
            # column (demux bleed-through) — those have already been mapped
            # to their canonical name by _apply_set1_validated_filter() and
            # must NOT get a WT_ prefix prepended.
            is_wt = col.startswith("WT_") or col.startswith("WT_embryo")
            # Also check: if the canonical name resolves to a non-WT group,
            # it should not get the WT_ prefix regardless of column name.
            resolved_group = SC.build_group_map([col]).get(col, "unknown") if _HAS_SC else "unknown"
            if is_wt or resolved_group == "WT":
                clean = col if col.startswith("WT_") else f"WT_{col}"
            else:
                # Non-WT embryo in the WT matrix — use canonical name as-is
                clean = col
            new = f"{args.wt_sample}__{clean}"
            renamed[col] = new
            meta_rows.append({"embryo": new,
                               **_get_metadata(clean, args.wt_sample, "set1")})
        wt_df = wt_df.rename(columns=renamed)
        all_dfs.append(wt_df)
        print(f"[step5b]   WT: {wt_df.shape[1]} embryo(s)")

    # ── Set1: Mutant ──────────────────────────────────────────────────────────
    if not args.no_set1:
        print(f"[step5b] Loading set1 Mutant: {args.mutant_sample}")
        mut_df = _load_matrix(args.mutant_sample)
        mut_df = _apply_set1_validated_filter(mut_df, args.mutant_sample)
        renamed = {}
        for col in mut_df.columns:
            new = f"{args.mutant_sample}__{col}"
            renamed[col] = new
            meta_rows.append({"embryo": new,
                               **_get_metadata(col, args.mutant_sample, "set1")})
        mut_df = mut_df.rename(columns=renamed)
        all_dfs.append(mut_df)
        print(f"[step5b]   Mutant: {mut_df.shape[1]} embryo(s)")

    # ── Set2 ──────────────────────────────────────────────────────────────────
    if not args.no_set2:
        print(f"[step5b] Loading set2: {args.set2_sample}")
        s2_df = _load_matrix(args.set2_sample)

        if set2_prefixes:
            keep  = [c for c in s2_df.columns
                     if any(c.startswith(p) for p in set2_prefixes)]
            print(f"[step5b]   set2 filter: {len(keep)} / {s2_df.shape[1]} kept")
            s2_df = s2_df[keep]

        # Apply manual exclusions at the source
        if _HAS_SC:
            s2_df = s2_df[SC.apply_manual_exclude(list(s2_df.columns), verbose=True)]

        renamed = {}
        for col in s2_df.columns:
            new = f"{args.set2_sample}__{col}"
            renamed[col] = new
            meta_rows.append({"embryo": new,
                               **_get_metadata(col, args.set2_sample, "set2")})
        s2_df = s2_df.rename(columns=renamed)
        all_dfs.append(s2_df)
        print(f"[step5b]   Set2: {s2_df.shape[1]} embryo(s)")

    if not all_dfs:
        sys.exit("[step5b] No data loaded")

    # ── Merge ─────────────────────────────────────────────────────────────────
    combined = (pd.concat(all_dfs, axis=1, join="outer")
                  .fillna(0)
                  .astype("int32"))
    combined.index.name = "GeneID"

    if args.annotate and os.path.exists(args.annotate):
        sym = pd.read_csv(args.annotate, sep="\t", index_col="GeneID")
        combined.insert(0, "GeneSymbol",
                        combined.index.map(sym.get("Symbol", sym.iloc[:, 0]))
                                      .fillna(""))

    embryo_cols = [c for c in combined.columns if c != "GeneSymbol"]

    # ── Save count matrix ─────────────────────────────────────────────────────
    out_path = os.path.join(combined_dir, f"{args.prefix}_counts_matrix.tsv")
    combined.to_csv(out_path, sep="\t")

    # ── Save enriched sample-groups table ────────────────────────────────────
    total_umis = combined[embryo_cols].sum()
    detected   = (combined[embryo_cols] > 0).sum()

    meta_df = pd.DataFrame(meta_rows).set_index("embryo")
    meta_df["total_umis"]     = total_umis
    meta_df["detected_genes"] = detected

    grp_path = os.path.join(combined_dir, f"{args.prefix}_sample_groups.tsv")
    # Keep legacy "Embryo / Sample" column names so step5c / step7 still work
    meta_df.reset_index().rename(
        columns={"embryo": "Embryo", "sample": "Sample"}
    ).to_csv(grp_path, sep="\t", index=False)

    # ── Console summary ───────────────────────────────────────────────────────
    print(f"\n[step5b] Combined: {combined.shape[0]:,} genes × {len(embryo_cols)} embryos")
    print(f"[step5b] Saved → {out_path}")
    print(f"[step5b] Groups → {grp_path}")

    print(f"\n[step5b] Group breakdown:")
    for grp, cnt in meta_df["group"].value_counts().items():
        stages = meta_df[meta_df["group"] == grp]["stage"].value_counts().to_dict()
        stage_str = "  ".join(f"{s}:{n}" for s, n in stages.items())
        print(f"  {grp:<18}: n={cnt:>2}   {stage_str}")

    print(f"\n[step5b] Per-embryo QC summary:")
    print(f"  {'Embryo':<60} {'Group':<16} {'Stage':<8} "
          f"{'UMIs':>10} {'Genes':>8}")
    for row in meta_rows:
        e = row["embryo"]
        bare = e.split("__")[-1]
        print(f"  {bare:<60} {row['group']:<16} {row['stage']:<8} "
              f"{int(total_umis[e]):>10,} {int(detected[e]):>8,}")

    print(f"\n[step5b] Next:")
    print(f"  python step5c_qc_filter.py --combined {args.prefix} --report-only")
    print(f"  python step6_tpm.py --combined {args.prefix} --qc --quantile")


if __name__ == "__main__":
    main()
