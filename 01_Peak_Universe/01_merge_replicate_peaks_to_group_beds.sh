#!/usr/bin/env bash
set -euo pipefail

################################################################################
# GENERATE MERGED PEAK BEDs — ALL GENOTYPE GROUPS
#
# Produces two merged BED files per genotype group:
#
#   <group>_merged_union.bed          — any replicate has the peak (permissive)
#   <group>_merged_reproducible.bed   — peak present in ≥2 replicates (stringent)
#
# Groups and their confirmed replicates:
#
#   Yw               — yellow white WT anchor (3 reps)
#   BOT              — BOT_D7 nc14b (2 reps)
#   BOT_hR           — Run-D7 confirmed hets, nc14b (2 reps)
#   BOT_hR_gastr     — BOTR_rD7_2 nc14gastr, confirmed het (1 rep; union only)
#   BOTR             — Run-D7 confirmed nulls, nc14b (2 reps)
#   BOTR_gastr       — BOTR_rD7_1 nc14gastr, confirmed null (1 rep; union only)
#   BOTC             — HLH_B6 Nc14b, Cic deleted non-vent (2 reps)
#   BOTC_nc14late    — HLH_B6 Nc14late confirmed rep1 + Mat_Run_B6_4 nc14b interim (2 reps)
#                      Mat_Run_B6_4 BOTC-like by UMAP; nc14b used as acceptable proxy.
#   BOTC_oR          — Mat_Run_B6_3 + Mat_Run_B6_5 nc14b confirmed pair (2 reps)
#                      Mat_Run_B6_4 excluded from BOTC_oR (clusters with BOTC; incomplete rescue)
#   BOTv_nc14b       — FoxL1-High Nc14b confirmed pair (2 reps)
#   BOTv_nc14late    — FoxL1-High Nc14late confirmed pair (2 reps)
#   BOTCv_nc14b      — HLH54F-High Nc14b rep4+5; rep1/rep2 excluded (2 reps)
#   BOTCv_nc14late   — HLH54F-High Nc14late confirmed pair (2 reps)
#   BOTv_nc14late    — FoxL1-High_3_nc14d_IR (rep1 NEW — replaces _4_nc14d) + Nc14late_rep2 (2 reps)
#   BOTv_gastr       — FoxL1-High_2_gastr_IR + FoxL1-High_3_gastr_IR (2 reps) [NEW GROUP]
#   BOT_nc14d        — D7_1_nc14d_IR (rep1 NEW 28366) + BOT_D7_2_nc14d (rep2) — UPGRADED to pair
#   BOT_hR_gastr     — BOTR_rD7_2 (rep1) + RunD7_3_gastr_IR (rep2 NEW 28368 confirmed het)
#
# Usage:
#   chmod +x 01_merge_replicate_peaks_to_group_beds.sh
#   bash 01_merge_replicate_peaks_to_group_beds.sh
#
# Run from: Enhanced_pseudopink_RuntD7integration/
# Requires: bedtools
################################################################################

# ── PATHS ─────────────────────────────────────────────────────────────────────
NP="${NARROWPEAK_DIR:-./ATAC_NarrowPeaks}"   # override with env var if needed
OUT_DIR="./Merged_Peak_BEDs"
LOG_DIR="./Merged_Peak_BEDs/logs"

mkdir -p "$OUT_DIR" "$LOG_DIR"

# ── PARAMETERS ────────────────────────────────────────────────────────────────
REPRO_MIN=2    # peaks in ≥ this many replicates → reproducible
MERGE_DIST=0   # bp between peaks to merge; 0 = touch-only

# ── VALIDATION ────────────────────────────────────────────────────────────────
echo "======================================================="
echo "GENERATE MERGED PEAK BEDs"
echo "======================================================="
echo ""

if ! command -v bedtools &>/dev/null; then
  echo "ERROR: bedtools not found in PATH"; exit 1
fi
if [ ! -d "$NP" ]; then
  echo "ERROR: NarrowPeak directory not found: $NP"; exit 1
fi

echo "NarrowPeak dir : $NP"
echo "Output dir     : $OUT_DIR"
echo "Reproducible   : peak in ≥${REPRO_MIN} replicates"
echo "Merge distance : ${MERGE_DIST} bp"
echo ""

# ── HELPER ────────────────────────────────────────────────────────────────────
# merge_group LABEL FILE [FILE ...]
merge_group() {
  local label="$1"; shift
  local files=("$@")
  local log="${LOG_DIR}/${label}.log"

  {
    echo "--- ${label} (${#files[@]} file(s) specified) ---"

    local found=()
    for f in "${files[@]}"; do
      if [ -f "$f" ]; then
        found+=("$f")
        echo "  ✓ $(basename "$f")"
      else
        echo "  ✗ MISSING: $f"
      fi
    done

    if [ ${#found[@]} -eq 0 ]; then
      echo "  SKIP — no files found for ${label}"
      echo ""; return 0
    fi

    local union_bed="${OUT_DIR}/${label}_merged_union.bed"
    local repro_bed="${OUT_DIR}/${label}_merged_reproducible.bed"

    # UNION
    cat "${found[@]}" \
      | awk 'BEGIN{OFS="\t"} !/^#/ && NF>=3 {print $1,$2,$3}' \
      | sort -k1,1 -k2,2n \
      | bedtools merge -d "$MERGE_DIST" \
      | awk -v lbl="$label" 'BEGIN{OFS="\t"} {print $1,$2,$3, lbl"_union_"NR, ".", "."}' \
      > "$union_bed"
    echo "  Union:         $(wc -l < "$union_bed") peaks  → $(basename "$union_bed")"

    # REPRODUCIBLE
    if [ ${#found[@]} -lt "$REPRO_MIN" ]; then
      echo "  Reproducible:  SKIP (only ${#found[@]} rep(s), need ≥${REPRO_MIN})"
      > "$repro_bed"
    else
      cat "${found[@]}" \
        | awk 'BEGIN{OFS="\t"} !/^#/ && NF>=3 {print $1,$2,$3}' \
        | sort -k1,1 -k2,2n \
        | bedtools merge -d "$MERGE_DIST" -c 1 -o count \
        | awk -v min="$REPRO_MIN" '$4 >= min {print $1,$2,$3}' \
        | awk -v lbl="$label" 'BEGIN{OFS="\t"} {print $1,$2,$3, lbl"_repro_"NR, ".", "."}' \
        > "$repro_bed"
      echo "  Reproducible:  $(wc -l < "$repro_bed") peaks  → $(basename "$repro_bed")"
    fi
    echo ""
  } 2>&1 | tee -a "$log"
}

# ── GROUPS — explicit file lists ──────────────────────────────────────────────

# Yw — yellow white WT anchor
merge_group "Yw" \
  "$NP/Yw_04_AB_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak" \
  "$NP/Yw_07_AB_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak" \
  "$NP/Yw_08_AB_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"

# BOT — BOT_D7 nc14b
merge_group "BOT" \
  "$NP/BOT_D7_1_nc14b_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak" \
  "$NP/BOT_D7_2_nc14b_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"

# BOT_hR — Run-D7 confirmed hets, nc14b only
merge_group "BOT_hR" \
  "$NP/Run-D7_Dm_ATAC_Nc14b_01_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak" \
  "$NP/Run-D7_Dm_ATAC_Nc14b_02_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"

# BOT_hR_gastr — confirmed het pair (rD7_2 rep1 + RunD7_3_gastr_IR rep2 NEW 28368)
merge_group "BOT_hR_gastr" \
  "$NP/BOTR_rD7_2_nc14gastr_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak" \
  "$NP/RunD7_3_gastr_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"

# BOTR — Run-D7 confirmed nulls, nc14b only
merge_group "BOTR" \
  "$NP/Run-D7_Dm_ATAC_Nc14b_03_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak" \
  "$NP/Run-D7_Dm_ATAC_Nc14b_07_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"

# BOTR_gastr — BOTR_rD7_1 nc14gastr, confirmed null (single rep; union only)
merge_group "BOTR_gastr" \
  "$NP/BOTR_rD7_1_nc14gastr_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"

# BOTC — HLH B6 Nc14b, Cic deleted non-ventralized
merge_group "BOTC" \
  "$NP/HLH_B6_Nc14b_01_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak" \
  "$NP/HLH_B6_Nc14b_02_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"

# BOTC_nc14late — confirmed rep1 (HLH_B6_Nc14late_01) + Mat_Run_B6_4 nc14b interim proxy
# Mat_Run_B6_4 clusters BOTC-like by UMAP (excluded from BOTC_oR as incomplete rescue),
# which makes it an acceptable interim pair for BOTC nc14late.
# nc14b timepoint used as proxy; staging difference is negligible (~few minutes).
# Replace rep2 with a confirmed nc14late BOTC sample when available.
merge_group "BOTC_nc14late" \
  "$NP/HLH_B6_Nc14late_01_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak" \
  "$NP/Mat_Run_B6_4_nc14b_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"

# BOTC_oR — Mat_Run_B6_3 + Mat_Run_B6_5 nc14b confirmed pair
# Mat_Run_B6_4 excluded from BOTC_oR (clusters with BOTC; incomplete Runt rescue).
# n=2 → reproducible BED is generated; union BED unchanged.
merge_group "BOTC_oR" \
  "$NP/Mat_Run_B6_3_nc14b_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak" \
  "$NP/Mat_Run_B6_5_nc14b_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"

# BOTv_nc14b — FoxL1-High Nc14b confirmed pair (r=0.9691)
merge_group "BOTv_nc14b" \
  "$NP/FoxL1-High_Dm_ATAC_Nc14b_rep3_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak" \
  "$NP/FoxL1-High_Dm_ATAC_Nc14b_07_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"

# BOTv_nc14late — rep1 REPLACED: FoxL1-High_3_nc14d_IR (28361) confirmed over FoxL1-High_4_nc14d
# FoxL1-High_4_nc14d had lower intra-group r in QC UMAP; 28361 showed better concordance.
merge_group "BOTv_nc14late" \
  "$NP/FoxL1-High_3_nc14d_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak" \
  "$NP/FoxL1-High_Dm_ATAC_Nc14late_rep2_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"

# BOTv_gastr — NEW confirmed pair (28362 + 28363)
merge_group "BOTv_gastr" \
  "$NP/FoxL1-High_2_gastr_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak" \
  "$NP/FoxL1-High_3_gastr_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"

# BOTCv_nc14b — HLH54F-High Nc14b confirmed pair (r=0.9683)
# rep1 excluded (non-ventralized), rep2 excluded (UMAP outlier)
merge_group "BOTCv_nc14b" \
  "$NP/HLH54F-High_Dm_ATAC_Nc14b_rep4_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak" \
  "$NP/HLH54F-High_Dm_ATAC_Nc14b_rep5_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"

# BOTCv_nc14late — HLH54F-High Nc14late confirmed pair (r=0.9371)
merge_group "BOTCv_nc14late" \
  "$NP/HLH54F-High_Dm_ATAC_Nc14late_rep1_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak" \
  "$NP/HLH54F-High_Dm_ATAC_Nc14late_rep2_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"

# BOT_nc14d — UPGRADED to confirmed pair (28366 + existing rep2)
# rep1 = D7_1_nc14d_IR (28366 confirmed), rep2 = BOT_D7_2_nc14d (was singleton)
# n=2 → reproducible BED now generated
merge_group "BOT_nc14d" \
  "$NP/D7_1_nc14d_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak" \
  "$NP/BOT_D7_2_nc14d_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"

# ── SUMMARY ───────────────────────────────────────────────────────────────────
echo "======================================================="
echo "DONE"
echo "======================================================="
echo ""
echo "Output files:"
ls -lh "$OUT_DIR"/*.bed 2>/dev/null \
  | awk '{printf "  %-55s  %s\n", $NF, $5}' \
  || echo "  (none produced)"
echo ""
echo "Next: pass these BEDs to computeMatrix -R alongside your DAR BEDs"
echo "      to show called peaks as a reference track in your heatmaps."
