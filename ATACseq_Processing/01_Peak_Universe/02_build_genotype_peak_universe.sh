#!/usr/bin/env bash
################################################################################
# 02_build_genotype_peak_universe.sh
# macOS bash 3.2 compatible (no associative arrays)
#
# Run from: Enhanced_pseudopink_RuntD7integration/
#
# REPLICATE UPDATE:
#   BOTC_oR nc14b: Mat_Run_B6_5 added as confirmed rep2.
#                  GroupPeaks_BOTC_oR.bed and FullUniverse now include both
#                  rep1 (Mat_Run_B6_3) and rep2 (Mat_Run_B6_5).
#   Mat_Run_B6_4:  EXCLUDED from BOTC_oR nc14b (clusters with BOTC; incomplete rescue).
#                  Re-used as BOTC nc14late interim rep2 (see below) — its BOTC-like
#                  clustering is appropriate there. nc14b timepoint used as proxy;
#                  staging difference from nc14late is negligible (~few minutes).
#   BOTv nc14late:  rep1 REPLACED — FoxL1-High_4_nc14d → FoxL1-High_3_nc14d_IR (28361 confirmed).
#                   FoxL1-High_4_nc14d had lower intra-group correlation in QC.
#   BOTv gastr:     NEW confirmed pair added.
#                   rep1 = FoxL1-High_2_gastr_IR (28362), rep2 = FoxL1-High_3_gastr_IR (28363).
#   BOT nc14d:      UPGRADED from singleton to confirmed pair.
#                   rep1 = D7_1_nc14d_IR (28366 confirmed), rep2 = BOT_D7_2_nc14d (existing).
#   BOT_hR_gastr:   UPGRADED from singleton to confirmed het pair.
#                   rep1 = BOTR_rD7_2 (existing), rep2 = RunD7_3_gastr_IR (28368 confirmed het).
#                  No reproducible BED generated for this group.
#   BOTC nc14late: NEW group added.
#                  rep1 = HLH_B6_Nc14late_01 (confirmed).
#                  rep2 = Mat_Run_B6_4 nc14b (interim; BOTC-like by UMAP, acceptable
#                         proxy until a second confirmed nc14late BOTC rep is available).
#
# NOTE: BOTCv st10_11 entries are intentionally commented out pending
#       completion of the expanded universe QC (Build_FullUniverse_withSt10_11.sh).
#
# OUTPUTS (all in ./ATAC_NarrowPeaks/):
#
#   FullUniverse_nc14b_AllGeno_union_peaks.bed      — all files merged
#   FullUniverse_nc14b_NonVent_union_peaks.bed       — BOT/BOT_hR/BOTR/BOTC/BOTC_oR/BOTC_nc14late + yw_WT
#   FullUniverse_nc14b_Vent_union_peaks.bed          — BOTv/BOTCv + yw_WT
#
#   Per-genotype group BEDs:
#     GroupPeaks_yw_WT.bed
#     GroupPeaks_BOTv.bed
#     GroupPeaks_BOTCv.bed
#     GroupPeaks_BOT.bed
#     GroupPeaks_BOT_hR.bed
#     GroupPeaks_BOTR.bed
#     GroupPeaks_BOTC.bed
#     GroupPeaks_BOTC_nc14late.bed   ← NEW: HLH_B6_Nc14late_01 + Mat_Run_B6_4 nc14b (interim)
#     GroupPeaks_BOTC_oR.bed         ← includes rep2 (Mat_Run_B6_5)
#
# SAMPLE MANIFEST (32 files = 28 prev + 2 BOTv_gastr + 1 BOT_nc14d_rep1 + 1 BOT_hR_gastr_rep2):
#   yw_WT:          Yw_07, Yw_08, Yw_04
#   BOTv nc14b:     Nc14b_rep3, Nc14b_07
#   BOTv nc14late:  FoxL1-High_3_nc14d (rep1 NEW — replaces FoxL1-High_4_nc14d), Nc14late_rep2 (rep2)
#   BOTv gastr:     FoxL1-High_2_gastr_IR (rep1 NEW), FoxL1-High_3_gastr_IR (rep2 NEW) [NEW GROUP]
#   BOTCv:          Nc14b_rep4, Nc14b_rep5, Nc14late_rep1, Nc14late_rep2
#   BOT nc14b:      D7_1_nc14b, D7_2_nc14b
#   BOT nc14d:      D7_1_nc14d_IR (rep1 NEW — 28366), D7_2_nc14d (rep2, was singleton)
#   BOT_hR:         Nc14b_01, Nc14b_02
#   BOT_hR_gastr:   BOTR_rD7_2 (confirmed het rep1), RunD7_3_gastr_IR (confirmed het rep2 NEW)
#   BOTR:           Nc14b_03, Nc14b_07, gastr_rD7_1 (confirmed null)
#   BOTC:           Nc14b_01, Nc14b_02
#   BOTC_nc14late:  HLH_B6_Nc14late_01 (confirmed) + Mat_Run_B6_4_nc14b (interim proxy)
#   BOTC_oR:        Mat_Run_B6_3 (rep1) + Mat_Run_B6_5 (rep2)
#   BOTC_oR nc14d:  Mat_Run_B6_1_nc14d, Mat_Run_B6_2_nc14d
#
# NOTE on gastr identity swap:
#   rD7_2 filename prefix says BOTR but is the confirmed Runt NULL → GROUP=BOTR
#   rD7_1 filename prefix says BOTR but is the confirmed Runt HET  → GROUP=BOT_hR
################################################################################

set -euo pipefail

PEAK_DIR="./ATAC_NarrowPeaks"
LOG_FILE="${PEAK_DIR}/FullUniverse_build.log"

command -v bedtools >/dev/null 2>&1 || { echo "ERROR: bedtools not found on PATH"; exit 1; }

echo "=======================================================" | tee    "$LOG_FILE"
echo "PEAK UNIVERSE BUILD -- $(date)"                          | tee -a "$LOG_FILE"
echo "=======================================================" | tee -a "$LOG_FILE"
echo ""                                                        | tee -a "$LOG_FILE"

################################################################################
# MANIFEST — LABELS, FILES, GROUP (must stay in exact positional sync)
################################################################################

LABELS=(
  # ── yw_WT ─────────────────────────────────────────────────────────────────
  "yw_WT_07"
  "yw_WT_08"
  "yw_WT_04"
  # ── BOTv Nc14b ────────────────────────────────────────────────────────────
  "BOTv_Nc14b_rep3"
  "BOTv_Nc14b_07"
  # ── BOTv Nc14late ─────────────────────────────────────────────────────────
  "BOTv_Nc14late_rep1"    # FoxL1-High_3_nc14d_IR (NEW — replaces FoxL1-High_4_nc14d)
  "BOTv_Nc14late_rep2"    # FoxL1-High_Dm_ATAC_Nc14late_rep2 (unchanged)
  # ── BOTv gastr — NEW confirmed pair ───────────────────────────────────────
  "BOTv_gastr_rep1"       # FoxL1-High_2_gastr_IR (28362 confirmed)
  "BOTv_gastr_rep2"       # FoxL1-High_3_gastr_IR (28363 confirmed)
  # ── BOTCv Nc14b ───────────────────────────────────────────────────────────
  "BOTCv_Nc14b_rep4"
  "BOTCv_Nc14b_rep5"
  # ── BOTCv Nc14late ────────────────────────────────────────────────────────
  "BOTCv_Nc14late_rep1"
  "BOTCv_Nc14late_rep2"
  # ── BOT Nc14b ─────────────────────────────────────────────────────────────
  "BOT_nc14b_rep1"
  "BOT_nc14b_rep2"
  # ── BOT Nc14d — UPGRADED to confirmed pair (28366 + existing rep2) ────────
  "BOT_nc14d_rep1"        # D7_1_nc14d_IR (28366 NEW confirmed)
  "BOT_nc14d_rep2"        # BOT_D7_2_nc14d_IR (was singleton, now rep2)
  # ── BOT_hR Nc14b ──────────────────────────────────────────────────────────
  "BOT_hR_nc14b_rep1"
  "BOT_hR_nc14b_rep2"
  # ── BOTR Nc14b ────────────────────────────────────────────────────────────
  "BOTR_nc14b_rep1"
  "BOTR_nc14b_rep2"
  # ── BOTR gastrulation — confirmed Runt null (rD7_2) ──────────────────────
  "BOTR_gastr_rep1"
  # ── BOT_hR gastrulation — confirmed Runt het (rD7_1) ─────────────────────
  # ── BOT_hR gastr — confirmed het pair (rD7_2 + RunD7_3_gastr confirmed het) ─
  "BOT_hR_gastr_rep1"     # BOTR_rD7_2 (confirmed het)
  "BOT_hR_gastr_rep2"     # RunD7_3_gastr_IR (28368 confirmed het)
  # ── BOTC Nc14b ────────────────────────────────────────────────────────────
  "BOTC_nc14b_rep1"
  "BOTC_nc14b_rep2"
  # ── BOTC Nc14late — rep1 confirmed; rep2 = Mat_Run_B6_4 nc14b interim proxy ─
  "BOTC_nc14late_rep1"   # HLH_B6_Nc14late_01 — confirmed nc14late
  "BOTC_nc14late_rep2"   # Mat_Run_B6_4 nc14b — interim; BOTC-like UMAP; acceptable proxy
  # ── BOTC_oR Nc14b — rep1 + rep2 ──────────────────────────────────────────
  "BOTC_oR_nc14b_rep1"
  "BOTC_oR_nc14b_rep2"   # Mat_Run_B6_5 — newly confirmed
  # ── BOTC_oR Nc14d ─────────────────────────────────────────────────────────
  "BOTC_oR_nc14d_rep1"
  "BOTC_oR_nc14d_rep2"
)

FILES=(
  # ── yw_WT ─────────────────────────────────────────────────────────────────
  "${PEAK_DIR}/Yw_04_AB_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  "${PEAK_DIR}/Yw_07_AB_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  "${PEAK_DIR}/Yw_08_AB_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  # ── BOTv Nc14b ────────────────────────────────────────────────────────────
  "${PEAK_DIR}/FoxL1-High_Dm_ATAC_Nc14b_rep3_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  "${PEAK_DIR}/FoxL1-High_Dm_ATAC_Nc14b_07_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  # ── BOTv Nc14late ─────────────────────────────────────────────────────────
  "${PEAK_DIR}/FoxL1-High_3_nc14d_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  "${PEAK_DIR}/FoxL1-High_Dm_ATAC_Nc14late_rep2_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  # ── BOTv gastr — NEW confirmed pair ───────────────────────────────────────
  "${PEAK_DIR}/FoxL1-High_2_gastr_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  "${PEAK_DIR}/FoxL1-High_3_gastr_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  # ── BOTCv Nc14b ───────────────────────────────────────────────────────────
  "${PEAK_DIR}/HLH54F-High_Dm_ATAC_Nc14b_rep4_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  "${PEAK_DIR}/HLH54F-High_Dm_ATAC_Nc14b_rep5_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  # ── BOTCv Nc14late ────────────────────────────────────────────────────────
  "${PEAK_DIR}/HLH54F-High_Dm_ATAC_Nc14late_rep1_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  "${PEAK_DIR}/HLH54F-High_Dm_ATAC_Nc14late_rep2_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  # ── BOT Nc14b ─────────────────────────────────────────────────────────────
  "${PEAK_DIR}/BOT_D7_1_nc14b_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  "${PEAK_DIR}/BOT_D7_2_nc14b_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  # ── BOT Nc14d — UPGRADED to confirmed pair ────────────────────────────────
  "${PEAK_DIR}/D7_1_nc14d_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  "${PEAK_DIR}/BOT_D7_2_nc14d_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  # ── BOT_hR Nc14b ──────────────────────────────────────────────────────────
  "${PEAK_DIR}/Run-D7_Dm_ATAC_Nc14b_01_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  "${PEAK_DIR}/Run-D7_Dm_ATAC_Nc14b_02_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  # ── BOTR Nc14b ────────────────────────────────────────────────────────────
  "${PEAK_DIR}/Run-D7_Dm_ATAC_Nc14b_03_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  "${PEAK_DIR}/Run-D7_Dm_ATAC_Nc14b_07_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  # ── BOTR gastrulation — confirmed Runt null (rD7_2) ──────────────────────
  "${PEAK_DIR}/BOTR_rD7_2_nc14gastr_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  "${PEAK_DIR}/RunD7_3_gastr_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  # ── BOT_hR gastrulation — confirmed Runt het (rD7_1) ─────────────────────
  "${PEAK_DIR}/BOTR_rD7_1_nc14gastr_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  # ── BOTC Nc14b ────────────────────────────────────────────────────────────
  "${PEAK_DIR}/HLH_B6_Nc14b_01_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  "${PEAK_DIR}/HLH_B6_Nc14b_02_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  # ── BOTC Nc14late — confirmed rep1 + Mat_Run_B6_4 nc14b interim proxy ────
  "${PEAK_DIR}/HLH_B6_Nc14late_01_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  "${PEAK_DIR}/Mat_Run_B6_4_nc14b_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  # ── BOTC_oR Nc14b — rep1 + rep2 (Mat_Run_B6_4 excluded: clusters with BOTC)
  "${PEAK_DIR}/Mat_Run_B6_3_nc14b_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  "${PEAK_DIR}/Mat_Run_B6_5_nc14b_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  # ── BOTC_oR Nc14d ─────────────────────────────────────────────────────────
  "${PEAK_DIR}/Mat_Run_B6_1_nc14d_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
  "${PEAK_DIR}/Mat_Run_B6_2_nc14d_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak"
)

# GROUP — one value per file, same order as FILES (26 entries total)
GROUP=(
  yw_WT   yw_WT   yw_WT      # yw_WT (3)
  BOTv    BOTv               # BOTv nc14b (2)
  BOTv    BOTv               # BOTv nc14late (2; rep1=FoxL1-High_3_nc14d NEW)
  BOTv_gastr  BOTv_gastr     # BOTv gastr (2 NEW; 28362 + 28363)
  BOTCv   BOTCv              # BOTCv nc14b (2)
  BOTCv   BOTCv              # BOTCv nc14late (2)
  BOT     BOT                # BOT nc14b (2)
  BOT     BOT                # BOT nc14d — UPGRADED confirmed pair (28366 + rep2)
  BOT_hR  BOT_hR             # BOT_hR nc14b (2)
  BOTR    BOTR               # BOTR nc14b (2)
  BOTR                       # BOTR gastr rD7_1 — confirmed null (1)
  BOT_hR_gastr  BOT_hR_gastr # BOT_hR gastr — confirmed het pair (rD7_2 + RunD7_3 NEW)
  BOTC    BOTC               # BOTC nc14b (2)
  BOTC_nc14late BOTC_nc14late  # BOTC nc14late: confirmed rep1 + Mat_Run_B6_4 interim (2)
  BOTC_oR BOTC_oR            # BOTC_oR nc14b rep1 + rep2 (2)
  BOTC_oR BOTC_oR            # BOTC_oR nc14d (2)
)

################################################################################
# SANITY CHECK — array lengths must match
################################################################################

N_LABELS=${#LABELS[@]}
N_FILES=${#FILES[@]}
N_GROUP=${#GROUP[@]}

if [[ $N_LABELS -ne $N_FILES || $N_FILES -ne $N_GROUP ]]; then
  echo "ERROR: array length mismatch" | tee -a "$LOG_FILE"
  echo "  LABELS=${N_LABELS}  FILES=${N_FILES}  GROUP=${N_GROUP}" | tee -a "$LOG_FILE"
  exit 1
fi
echo "Array check: LABELS=${N_LABELS}  FILES=${N_FILES}  GROUP=${N_GROUP}  ✓" \
  | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

################################################################################
# VALIDATE INPUT FILES
################################################################################

echo "Validating input files..." | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

FOUND=0; MISSING=0
for i in "${!LABELS[@]}"; do
  label="${LABELS[$i]}"; f="${FILES[$i]}"
  if [[ -f "$f" ]]; then
    printf "  OK      [%-22s]  %s\n" "$label" "$(basename "$f")" | tee -a "$LOG_FILE"
    FOUND=$((FOUND+1))
  else
    printf "  MISSING [%-22s]  %s\n" "$label" "$f" | tee -a "$LOG_FILE"
    MISSING=$((MISSING+1))
  fi
done

echo "" | tee -a "$LOG_FILE"
echo "Found: ${FOUND}  |  Missing: ${MISSING}" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"
[[ $MISSING -gt 0 ]] && \
  echo "WARNING: ${MISSING} file(s) missing — they will be skipped." | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

################################################################################
# HELPER: merge_peaks <output_bed> <file1> [file2 ...]
# Strips chrM, sorts, bedtools merge, writes chr:start-end as name column.
################################################################################

merge_peaks() {
  local out_bed="$1"; shift
  local tmp_cat; tmp_cat=$(mktemp /tmp/peak_cat.XXXXXX)
  local tmp_srt; tmp_srt=$(mktemp /tmp/peak_srt.XXXXXX)
  local tmp_mrg; tmp_mrg=$(mktemp /tmp/peak_mrg.XXXXXX)

  for f in "$@"; do
    [[ -f "$f" ]] || continue
    grep -v "^chrM" "$f" | cut -f1-3 >> "$tmp_cat"
  done

  local raw_n; raw_n=$(wc -l < "$tmp_cat" | tr -d ' ')
  if [[ $raw_n -eq 0 ]]; then
    printf "  WARNING: no peaks for %s — skipping\n" \
      "$(basename "$out_bed")" | tee -a "$LOG_FILE"
    rm -f "$tmp_cat" "$tmp_srt" "$tmp_mrg"
    return 0
  fi

  sort -k1,1 -k2,2n "$tmp_cat" > "$tmp_srt"
  bedtools merge -i "$tmp_srt"   > "$tmp_mrg"

  local merged_n; merged_n=$(wc -l < "$tmp_mrg" | tr -d ' ')
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$1":"$2"-"$3}' "$tmp_mrg" > "$out_bed"
  rm -f "$tmp_cat" "$tmp_srt" "$tmp_mrg"

  printf "  raw=%-8s merged=%-8s → %s\n" \
    "$raw_n" "$merged_n" "$(basename "$out_bed")" | tee -a "$LOG_FILE"
}

################################################################################
# STEP 1: FULL UNIVERSE
################################################################################

echo "STEP 1: Full universe (all 28 files)..." | tee -a "$LOG_FILE"
ALL_FILES=()
for i in "${!FILES[@]}"; do
  [[ -f "${FILES[$i]}" ]] && ALL_FILES+=("${FILES[$i]}")
done
OUT_ALL="${PEAK_DIR}/FullUniverse_nc14b_AllGeno_union_peaks.bed"
merge_peaks "$OUT_ALL" "${ALL_FILES[@]}"
echo "" | tee -a "$LOG_FILE"

################################################################################
# STEP 2: NON-VENT SUB-UNIVERSE (BOT / BOT_hR / BOTR / BOTC / BOTC_oR / yw_WT)
################################################################################

echo "STEP 2: Non-vent sub-universe..." | tee -a "$LOG_FILE"
NV_FILES=()
for i in "${!FILES[@]}"; do
  [[ -f "${FILES[$i]}" ]] || continue
  case "${GROUP[$i]}" in BOT|BOT_hR|BOTR|BOTC|BOTC_nc14late|BOTC_oR|yw_WT)
    NV_FILES+=("${FILES[$i]}") ;; esac
done
OUT_NV="${PEAK_DIR}/FullUniverse_nc14b_NonVent_union_peaks.bed"
merge_peaks "$OUT_NV" "${NV_FILES[@]}"
echo "" | tee -a "$LOG_FILE"

################################################################################
# STEP 3: VENT SUB-UNIVERSE (BOTv / BOTCv / yw_WT)
################################################################################

echo "STEP 3: Vent sub-universe..." | tee -a "$LOG_FILE"
V_FILES=()
for i in "${!FILES[@]}"; do
  [[ -f "${FILES[$i]}" ]] || continue
  case "${GROUP[$i]}" in BOTv|BOTv_gastr|BOTCv|yw_WT)
    V_FILES+=("${FILES[$i]}") ;; esac
done
OUT_V="${PEAK_DIR}/FullUniverse_nc14b_Vent_union_peaks.bed"
merge_peaks "$OUT_V" "${V_FILES[@]}"
echo "" | tee -a "$LOG_FILE"

################################################################################
# STEP 4: PER-GENOTYPE GROUP BEDs
################################################################################

echo "STEP 4: Per-genotype group BEDs..." | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

for grp in yw_WT BOTv BOTv_gastr BOTCv BOT BOT_hR BOTR BOT_hR_gastr BOTC BOTC_nc14late BOTC_oR; do
  grp_files=()
  for i in "${!FILES[@]}"; do
    [[ "${GROUP[$i]}" == "$grp" && -f "${FILES[$i]}" ]] && \
      grp_files+=("${FILES[$i]}")
  done
  if [[ ${#grp_files[@]} -gt 0 ]]; then
    out_grp="${PEAK_DIR}/GroupPeaks_${grp}.bed"
    printf "  %-12s (%d files):  " "$grp" "${#grp_files[@]}" | tee -a "$LOG_FILE"
    merge_peaks "$out_grp" "${grp_files[@]}"
  else
    echo "  ${grp}: no files found — skipping" | tee -a "$LOG_FILE"
  fi
done

echo "" | tee -a "$LOG_FILE"

################################################################################
# SUMMARY
################################################################################

echo "=======================================================" | tee -a "$LOG_FILE"
echo "COMPLETE -- $(date)"                                     | tee -a "$LOG_FILE"
echo "=======================================================" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

for bed in \
    "$OUT_ALL" "$OUT_NV" "$OUT_V" \
    "${PEAK_DIR}/GroupPeaks_yw_WT.bed"   \
    "${PEAK_DIR}/GroupPeaks_BOTv.bed"    \
    "${PEAK_DIR}/GroupPeaks_BOTCv.bed"   \
    "${PEAK_DIR}/GroupPeaks_BOT.bed"     \
    "${PEAK_DIR}/GroupPeaks_BOT_hR.bed"  \
    "${PEAK_DIR}/GroupPeaks_BOTR.bed"    \
    "${PEAK_DIR}/GroupPeaks_BOTC.bed"    \
    "${PEAK_DIR}/GroupPeaks_BOTC_nc14late.bed" \
    "${PEAK_DIR}/GroupPeaks_BOTC_oR.bed" ; do
  if [[ -f "$bed" ]]; then
    n=$(wc -l < "$bed" | tr -d ' ')
    printf "  %-58s  %s peaks\n" "$(basename "$bed")" "$n" | tee -a "$LOG_FILE"
  else
    printf "  %-58s  MISSING\n" "$(basename "$bed")" | tee -a "$LOG_FILE"
  fi
done

echo "" | tee -a "$LOG_FILE"
echo "Peaks per chromosome (full universe):" | tee -a "$LOG_FILE"
cut -f1 "$OUT_ALL" | sort | uniq -c | sort -k2,2V | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"
echo "Log: $LOG_FILE"
