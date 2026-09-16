#!/bin/bash
set -euo pipefail

################################################################################
# 01_frip_qc.sh
#
# PURPOSE:
#   Local-machine FRiP (fraction of reads in peaks) calculator, using your
#   already-called MACS2 narrow/broad peaks and the BAMs sitting in
#   Control_Bams/ (same layout as 01_dar_calling_seven_genotype_nc14b.r
#   and 02_dar_calling_ventralized_temporal.r).
#
#   Does NOT re-call peaks or touch the cluster -- just flagstat + bedtools
#   intersect against files you already have locally.
#
# BAM LOOKUP:
#   Checks, in order:
#     ./Control_Bams/Fully_flattened_potential/<BASENAME><bamsuffix>
#     ./Control_Bams/<BASENAME><bamsuffix>
#   (default <bamsuffix> = _noq_rmdup.noChrM.bam, matching the DAR scripts)
#
# PEAK LOOKUP:
#   Globs <peakdir> for any file containing BASENAME and ending in
#   .narrowPeak or .broadPeak (whichever --peaktype you pick) -- this
#   doesn't assume your local files follow the cluster's exact naming
#   convention from Step04_ATAC_MACS2_Narrow-Broad_Peakcalling.sub.
#
# USAGE:
#   bash 01_frip_qc.sh --peakdir ./ATAC_NarrowPeaks BASENAME1 [BASENAME2 ...]
#   bash 01_frip_qc.sh --peakdir ./ATAC_NarrowPeaks --list samples.txt
#
#   Optional flags:
#     --peaktype  narrow|broad     (default: narrow)
#     --bamsuffix <suffix>         (default: _noq_rmdup.noChrM.bam)
#
# Run from the directory that contains Control_Bams/ (same working dir as
# the DAR scripts), or the BAM lookup paths above won't resolve.
################################################################################

PEAKDIR=""
PEAKTYPE="narrow"
BAMSUFFIX="_noq_rmdup.noChrM.bam"
BASENAMES=()
LIST_FILE=""

BAMDIR1="./Control_Bams/Fully_flattened_potential"
BAMDIR2="./Control_Bams"
BAMDIR3="../Enhanced_Run_wHets/INDV_Bams"   # yw control BAMs live here, not in Control_Bams

while [[ $# -gt 0 ]]; do
  case "$1" in
    --peakdir)
      PEAKDIR="$2"; shift 2 ;;
    --peaktype)
      PEAKTYPE="$2"; shift 2 ;;
    --bamsuffix)
      BAMSUFFIX="$2"; shift 2 ;;
    --list)
      LIST_FILE="$2"; shift 2 ;;
    -h|--help)
      grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *)
      BASENAMES+=("$1"); shift ;;
  esac
done

if [[ -z "$PEAKDIR" ]]; then
  echo "ERROR: --peakdir is required (folder containing your narrow/broad MACS2 peak files)"
  exit 1
fi
if [[ ! -d "$PEAKDIR" ]]; then
  echo "ERROR: peak directory not found: $PEAKDIR"
  exit 1
fi

if [[ -n "$LIST_FILE" ]]; then
  if [[ ! -f "$LIST_FILE" ]]; then
    echo "ERROR: --list file not found: $LIST_FILE"; exit 1
  fi
  while IFS= read -r line; do
    [[ -z "$line" ]] && continue
    BASENAMES+=("$line")
  done < "$LIST_FILE"
fi

if [[ ${#BASENAMES[@]} -eq 0 ]]; then
  echo "ERROR: no BASENAMEs given."
  echo "Usage: bash 01_frip_qc.sh --peakdir <dir> BASENAME1 [BASENAME2 ...]   or   --list samples.txt"
  exit 1
fi

if [[ "$PEAKTYPE" != "narrow" && "$PEAKTYPE" != "broad" ]]; then
  echo "ERROR: --peaktype must be 'narrow' or 'broad'"
  exit 1
fi
PEAKEXT="narrowPeak"
[[ "$PEAKTYPE" == "broad" ]] && PEAKEXT="broadPeak"

for TOOL in samtools bedtools; do
  command -v "$TOOL" &>/dev/null || { echo "ERROR: $TOOL not on PATH -- activate your atac conda env first"; exit 1; }
done

OUT_TSV="./qc_frip_local_${PEAKTYPE}.tsv"
printf "sample\tpeaktype\tbam_path\tpeak_path\ttotal_mapped\treads_in_peaks\tFRiP\tn_peaks\n" > "$OUT_TSV"

echo "Computing FRiP (${PEAKTYPE} peaks) against local BAMs..."
echo ""

for BASENAME in "${BASENAMES[@]}"; do
  echo "-- ${BASENAME}"

  # Find BAM in either known Control_Bams location
  BAM=""
  for CAND in "${BAMDIR1}/${BASENAME}${BAMSUFFIX}" "${BAMDIR2}/${BASENAME}${BAMSUFFIX}" "${BAMDIR3}/${BASENAME}${BAMSUFFIX}"; do
    if [[ -f "$CAND" ]]; then BAM="$CAND"; break; fi
  done
  if [[ -z "$BAM" ]]; then
    echo "   SKIP: BAM not found (looked for ${BASENAME}${BAMSUFFIX} in Control_Bams locations and ${BAMDIR3})"
    continue
  fi

  # Find peak file: any file in PEAKDIR containing BASENAME, right extension
  PEAKFILE=$(find "$PEAKDIR" -maxdepth 1 -type f -iname "*${BASENAME}*.${PEAKEXT}" | head -n 1)
  if [[ -z "$PEAKFILE" ]]; then
    echo "   SKIP: no *${BASENAME}*.${PEAKEXT} file found in $PEAKDIR"
    continue
  fi

  TOTAL=$(samtools flagstat "$BAM" | awk '/ mapped \(/{print $1; exit}')
  IN_PEAKS=$(bedtools intersect -u -a "$BAM" -b "$PEAKFILE" | samtools view -c -)
  NPEAKS=$(wc -l < "$PEAKFILE")

  FRIP=$(awk -v a="$IN_PEAKS" -v b="$TOTAL" \
    'BEGIN{ if (b>0) printf "%.4f", a/b; else print "NA" }')

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$BASENAME" "$PEAKTYPE" "$BAM" "$PEAKFILE" "$TOTAL" "$IN_PEAKS" "$FRIP" "$NPEAKS" \
    >> "$OUT_TSV"
done

echo ""
echo "Done -> $OUT_TSV"
echo ""
column -t "$OUT_TSV"
