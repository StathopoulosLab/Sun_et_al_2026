################################################################################
# SEVEN-GENOTYPE nc14b DAR ANALYSIS  —  v4
#
# Previous enhancements (v2):
#   Tier 1 named outputs | baseMean/AveExpr | per-contrast extension_bp |
#   priority-scored deduplication
#
# Changes in v3 (5 structural improvements — preserved in v4):
#   [L1]-[L5]+[L7] as documented below.
#
# New in v4 (performance — requires Precount_AllCounts_Matrices.r run first):
#   - All featureCounts calls replaced by cache loads from counts_cache/
#   - FC_THREADS read from SLURM_CPUS_PER_TASK (used in fallback path only)
#   - grp_ext_cache pre-loaded; on-the-fly extended Tier 1 counting removed
#     from run_contrast() — eliminates ~19 redundant featureCounts calls
#   - use_existing_counts always TRUE; narrow counts loaded from disk
#
# New in v3 (5 structural improvements):
#
#  [L1] GROUP-SPECIFIC PEAKS FOR TIER 1
#       Tier 1 condition-specific detection now counts reads in peaks called
#       from each group's own BAMs (GroupPeaks_*.bed from the shell script)
#       rather than the full union. This matches the Comprehensive DAR
#       approach and recovers peaks that exist in one condition but whose
#       boundaries in the union are dominated by the other condition.
#
#  [L2] CPM-SCALED TIER 1 THRESHOLDS
#       min_reads and max_reads thresholds are now expressed as CPM and
#       scaled per-contrast to replicate count and median library size,
#       so n=1 singleton groups are treated appropriately and large-library
#       samples don't inflate apparent signal.
#
#  [L3] CONTRAST-SPECIFIC RUV CONTROL SELECTION
#       RUVg control peaks are selected from only the two groups in each
#       contrast (pairwise edgeR), not globally across all 7 genotypes.
#       This avoids over-constraining the normalisation and recovers signal
#       that is stable within a contrast but variable across the full design.
#
#  [L4] BOUNDARY-BASED PEAK EXTENSION
#       Extended peaks grow from the peak start/end boundaries rather than
#       the center. This better captures regulatory elements where the
#       ATAC summit is offset toward the NFR, and avoids including equal
#       amounts of flanking nucleosomal DNA on both sides when the accessible
#       region is asymmetric. Peaks that become adjacent after extension are
#       merged before counting.
#
#  [L5] CONTINUOUS CONFIDENCE SCORE
#       Each DAR carries a conf_score = method_agreement × tier_score ×
#       −log10(p). method_agreement is 1.0 (single method), 1.5 (DESeq2
#       and limma agree directionally), or 2.0 (all four calls agree: narrow
#       DESeq2 + narrow limma + extended DESeq2 + extended limma). This score
#       propagates into background subtraction via scored_consistent_overlaps(),
#       which weights overlaps by the minimum conf_score of the pair rather
#       than treating all overlaps as equally valid.
#
#  [L7] TIER 3 DIAGNOSTIC
#       For each contrast, Tier 3-only DARs (found in extended but not narrow
#       analysis) are classified by pile distribution pattern into:
#         tight_boundary  — narrow peak exists in group-specific peaks;
#                           extension merely captured its full extent
#         adjacent_merged — a second narrow peak within ext_bp bp contributed
#                           to the extended region signal
#         diffuse         — genuinely broad accessibility; no narrow peak nearby
#       Results written to Diagnostics/ subdirectory.
#
# Genotype naming:
#   BOTv    FoxL1-High           triple + tolrm9, Cic intact, Runt GoF     [vent]
#   BOTCv   HLH54F-High vent.    triple + tolrm9, Cic deleted, Runt absent [vent]
#   BOT     D7 WT control        triple, Cic intact, Runt intact (2 copies)[non-vent]
#   BOT_hR  Run-D7 het           triple, Cic intact, Runt het              [non-vent]
#   BOTR    Run-D7 homo          triple, Cic intact, Runt null             [non-vent]
#   BOTC    HLH B6 non-vent      triple, Cic deleted, Runt transient       [non-vent]
#   BOTC_oR mat-tub>Run in B6    triple, Cic deleted, Runt overexpressed   [non-vent]
#
# REPLICATE UPDATE vs previous version:
#   BOTC_oR nc14b: Mat_Run_B6_5 confirmed as rep2 (QC r vs rep1=0.9405;
#                  nearest confirmed neighbour = BOTC_oR_nc14b_rep1).
#                  BOTC_oR is now n=2 → DESeq2 + limma both active for
#                  BOTC_oR_vs_BOTC and BOTC_oR_vs_BOTR contrasts.
#                  (Previously n=1 caused DESeq2 to return NA p-values;
#                  results were limma-only in practice.)
#   Mat_Run_B6_4:  EXCLUDED — clusters with BOTC (incomplete Runt rescue);
#                  excluded in all previous versions, remains excluded.
################################################################################

suppressPackageStartupMessages({
  library(Rsubread)
  library(DESeq2)
  library(limma)
  library(edgeR)
  library(RUVSeq)
  library(GenomicRanges)
  library(rtracklayer)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

################################################################################
# CONFIGURATION
################################################################################

# Thread count — reads SLURM allocation automatically; falls back to 4.
# Only used in the fallback featureCounts paths (first-run or missing cache).
FC_THREADS <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))

# Cache directory written by Precount_AllCounts_Matrices.r.
# Must be run once before this script if counts_cache/ does not yet exist.
CACHE_DIR <- "./counts_cache/SevenGeno_nc14b"

base_dir   <- "../Generate_fresh_counts"
output_dir <- file.path(base_dir, "Output/SevenGeno_nc14b_v3")
tier1_dir  <- file.path(output_dir, "Tier1_ConditionSpecific")
diag_dir   <- file.path(output_dir, "Diagnostics")
for (d in c(output_dir, tier1_dir, diag_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# BAM FILES
# ---------------------------------------------------------------------------
bamfiles <- c(
  "./Control_Bams/Fully_flattened_potential/FoxL1-High_Dm_ATAC_Nc14b_rep3_noq_rmdup.noChrM.bam",
  "./Control_Bams/Fully_flattened_potential/FoxL1-High_Dm_ATAC_Nc14b_07_noq_rmdup.noChrM.bam",
  "./Control_Bams/Fully_flattened_potential/HLH54F-High_Dm_ATAC_Nc14b_rep4_noq_rmdup.noChrM.bam",
  "./Control_Bams/Fully_flattened_potential/HLH54F-High_Dm_ATAC_Nc14b_rep5_noq_rmdup.noChrM.bam",
  "./Control_Bams/BOT_D7_1_nc14b_IR_noq_rmdup.noChrM.bam",
  "./Control_Bams/BOT_D7_2_nc14b_IR_noq_rmdup.noChrM.bam",
  "./Control_Bams/Run-D7_Dm_ATAC_Nc14b_01_noq_rmdup.noChrM.bam",
  "./Control_Bams/Run-D7_Dm_ATAC_Nc14b_02_noq_rmdup.noChrM.bam",
  "./Control_Bams/Run-D7_Dm_ATAC_Nc14b_03_noq_rmdup.noChrM.bam",
  "./Control_Bams/Run-D7_Dm_ATAC_Nc14b_07_noq_rmdup.noChrM.bam",
  "./Control_Bams/HLH_B6_Nc14b_01_noq_rmdup.noChrM.bam",
  "./Control_Bams/HLH_B6_Nc14b_02_noq_rmdup.noChrM.bam",
  "./Control_Bams/Mat_Run_B6_3_nc14b_IR_noq_rmdup.noChrM.bam",  # BOTC_oR rep1 (confirmed)
  # Mat_Run_B6_4 EXCLUDED — clusters with BOTC (incomplete Runt rescue)
  "./Control_Bams/Mat_Run_B6_5_nc14b_IR_noq_rmdup.noChrM.bam"   # BOTC_oR rep2 (confirmed; r vs rep1=0.9405)
)
sample_names <- basename(bamfiles)

# BOTC_oR is now n=2 — DESeq2 dispersion estimation is valid for both
# BOTC_oR_vs_BOTC and BOTC_oR_vs_BOTR contrasts. Previously n=1 caused
# DESeq2 to return NA p-values; limma was the only active method.
coldata <- data.frame(
  sample   = sample_names,
  genotype = factor(
    c("BOTv","BOTv","BOTCv","BOTCv","BOT","BOT",
      "BOT_hR","BOT_hR","BOTR","BOTR","BOTC","BOTC","BOTC_oR","BOTC_oR"),
    levels = c("BOTv","BOTCv","BOT","BOT_hR","BOTR","BOTC","BOTC_oR")),
  row.names = sample_names, stringsAsFactors = FALSE)
coldata$group <- factor(coldata$genotype,
  levels = c("BOTv","BOTCv","BOT","BOT_hR","BOTR","BOTC","BOTC_oR"))

# ---------------------------------------------------------------------------
# PEAK UNIVERSES
# Full union  — used for Tier 2/3 differential analysis (shared coordinate space)
# Group BEDs  — used for Tier 1 condition-specific detection (L1)
# ---------------------------------------------------------------------------
peak_dir              <- "./ATAC_NarrowPeaks"
full_universe_bed     <- file.path(peak_dir, "FullUniverse_nc14b_AllGeno_union_peaks.bed")
group_peak_beds       <- list(
  BOTv    = file.path(peak_dir, "GroupPeaks_BOTv.bed"),
  BOTCv   = file.path(peak_dir, "GroupPeaks_BOTCv.bed"),
  BOT     = file.path(peak_dir, "GroupPeaks_BOT.bed"),
  BOT_hR  = file.path(peak_dir, "GroupPeaks_BOT_hR.bed"),
  BOTR    = file.path(peak_dir, "GroupPeaks_BOTR.bed"),
  BOTC    = file.path(peak_dir, "GroupPeaks_BOTC.bed"),
  BOTC_oR = file.path(peak_dir, "GroupPeaks_BOTC_oR.bed")  # now contains rep1 + rep2 (n=2)
)
existing_counts_file  <- file.path(output_dir, "SevenGeno_nc14b_counts_matrix.txt")
# NOTE: use_existing_counts is no longer used in v4 — cache loading is always
# attempted from CACHE_DIR first, with featureCounts as fallback.
# existing_counts_file is retained only so Precount_AllCounts_Matrices.r can
# write there for compatibility with other downstream scripts.

# ---------------------------------------------------------------------------
# GLOBAL THRESHOLDS
# ---------------------------------------------------------------------------
log2fc_threshold     <- 0.5
pvalue_threshold     <- 0.05
# Tier 1 CPM thresholds (L2) — expressed as CPM, scaled per-contrast at runtime
min_cpm_present      <- 1.0   # minimum CPM in the "present" condition
max_cpm_absent       <- 0.25  # maximum CPM per replicate in the "absent" condition
extension_bp_global  <- 200
use_ruv              <- TRUE
ruv_k                <- 1

################################################################################
# HELPER FUNCTIONS — CORE
################################################################################

import_bed <- function(bed_file) {
  if (!file.exists(bed_file)) stop(paste("File not found:", bed_file))
  d  <- read.table(bed_file, sep="\t", header=FALSE,
                   stringsAsFactors=FALSE, comment.char="")
  gr <- GRanges(seqnames=d[,1],
                ranges=IRanges(start=d[,2]+1, end=d[,3]),
                strand="*")
  if (ncol(d) >= 4) names(gr) <- d[,4]
  else names(gr) <- paste0("peak_", seq_along(gr))
  gr
}

create_saf <- function(peaks_gr, name_prefix="peak") {
  if (is.null(names(peaks_gr)))
    names(peaks_gr) <- paste0(name_prefix, "_", seq_along(peaks_gr))
  peaks_gr <- trim(peaks_gr)
  saf <- data.frame(GeneID=names(peaks_gr),
                    Chr=as.character(seqnames(peaks_gr)),
                    Start=start(peaks_gr), End=end(peaks_gr),
                    Strand=".", stringsAsFactors=FALSE)
  bad <- apply(saf, 1, function(x) any(is.na(x)))
  if (any(bad)) saf <- saf[!bad,]
  saf$GeneID <- make.unique(as.character(saf$GeneID))
  saf$Start  <- as.integer(saf$Start); saf$End <- as.integer(saf$End)
  inv <- saf$Start < 1 | saf$End <= saf$Start | is.na(saf$Start) | is.na(saf$End)
  if (any(inv)) saf <- saf[!inv,]
  bc  <- grepl("\\s", saf$Chr) | saf$Chr == "" | is.na(saf$Chr)
  if (any(bc)) saf <- saf[!bc,]
  saf[, c("GeneID","Chr","Start","End","Strand")]
}

count_reads_in_peaks <- function(peaks_gr, bamfiles, name_prefix="peak") {
  cat("  Counting reads in", length(peaks_gr), "peaks across",
      length(bamfiles), "BAMs...\n")
  saf <- create_saf(peaks_gr, name_prefix)
  fc  <- featureCounts(files=bamfiles, annot.ext=saf, isPairedEnd=TRUE,
                       countMultiMappingReads=FALSE, primaryOnly=TRUE, nthreads=FC_THREADS)
  counts <- fc$counts
  colnames(counts) <- basename(bamfiles)
  list(counts=counts, valid_peak_ids=saf$GeneID)
}

# [L4] Boundary-based extension — extends from start/end, then reduces overlaps
# This captures asymmetric NFR positioning better than center-based extension.
build_extended_peaks_boundary <- function(peaks_gr, ext_bp) {
  # Shift start left and end right by ext_bp
  extended <- GRanges(
    seqnames = seqnames(peaks_gr),
    ranges   = IRanges(
      start = pmax(1L, start(peaks_gr) - ext_bp),
      end   = end(peaks_gr) + ext_bp),
    strand = "*")
  names(extended) <- names(peaks_gr)
  # Trim to chromosome limits, then merge overlapping windows
  extended <- trim(extended)
  # Preserve original names by taking the first name in each merged group
  merged <- reduce(extended, with.revmap=TRUE)
  # Name each merged region after the first constituent peak
  orig_names  <- names(peaks_gr)
  merged_names <- sapply(mcols(merged)$revmap, function(idx) orig_names[idx[1]])
  names(merged) <- paste0("bext", ext_bp, "_", merged_names)
  merged
}

################################################################################
# CACHE COLUMN ALIGNMENT
# The precount script stores counts with short superset keys (e.g. "BOTv_rep1").
# This script uses BAM basenames as sample_names. CACHE_KEY_MAP translates
# short keys → BAM basenames so every cache load gets correctly named columns.
#
# If you add or rename BAMs, update CACHE_KEY_MAP to match.
################################################################################

CACHE_KEY_MAP <- c(
  BOTv_rep1          = "FoxL1-High_Dm_ATAC_Nc14b_rep3_noq_rmdup.noChrM.bam",
  BOTv_rep2          = "FoxL1-High_Dm_ATAC_Nc14b_07_noq_rmdup.noChrM.bam",
  BOTCv_rep1         = "HLH54F-High_Dm_ATAC_Nc14b_rep4_noq_rmdup.noChrM.bam",
  BOTCv_rep2         = "HLH54F-High_Dm_ATAC_Nc14b_rep5_noq_rmdup.noChrM.bam",
  BOT_nc14b_rep1     = "BOT_D7_1_nc14b_IR_noq_rmdup.noChrM.bam",
  BOT_nc14b_rep2     = "BOT_D7_2_nc14b_IR_noq_rmdup.noChrM.bam",
  BOT_hR_rep1        = "Run-D7_Dm_ATAC_Nc14b_01_noq_rmdup.noChrM.bam",
  BOT_hR_rep2        = "Run-D7_Dm_ATAC_Nc14b_02_noq_rmdup.noChrM.bam",
  BOTR_rep1          = "Run-D7_Dm_ATAC_Nc14b_03_noq_rmdup.noChrM.bam",
  BOTR_rep2          = "Run-D7_Dm_ATAC_Nc14b_07_noq_rmdup.noChrM.bam",
  BOTC_rep1          = "HLH_B6_Nc14b_01_noq_rmdup.noChrM.bam",
  BOTC_rep2          = "HLH_B6_Nc14b_02_noq_rmdup.noChrM.bam",
  BOTC_oR_rep1       = "Mat_Run_B6_3_nc14b_IR_noq_rmdup.noChrM.bam",
  BOTC_oR_rep2       = "Mat_Run_B6_5_nc14b_IR_noq_rmdup.noChrM.bam"
)

realign_cache_cols <- function(counts_mat, target_names,
                                key_map = CACHE_KEY_MAP) {
  if (is.null(counts_mat)) return(NULL)
  cache_cols <- colnames(counts_mat)
  # Already aligned — nothing to do
  if (identical(cache_cols, target_names)) return(counts_mat)
  # Direct name match (e.g. cache was written with BAM basenames directly)
  matched <- match(target_names, cache_cols)
  if (!anyNA(matched)) return(counts_mat[, matched, drop=FALSE])
  # Translate via key_map: rename cache columns to BAM basenames, then match
  translated <- key_map[cache_cols]          # NA for any key not in map
  known      <- !is.na(translated)
  colnames(counts_mat)[known] <- translated[known]
  matched <- match(target_names, colnames(counts_mat))
  if (anyNA(matched)) {
    warning("realign_cache_cols: still unmatched after key_map translation:\n  ",
            paste(target_names[is.na(matched)], collapse="\n  "),
            "\nCache columns (after translation): ",
            paste(colnames(counts_mat), collapse=", "))
    return(NULL)
  }
  counts_mat[, matched, drop=FALSE]
}

# Apply to all count matrices in a grp_cache or grp_ext_cache list
realign_grp_cache <- function(cache_obj, target_names,
                               key_map = CACHE_KEY_MAP) {
  if (is.null(cache_obj)) return(NULL)
  cache_obj$counts <- lapply(cache_obj$counts, realign_cache_cols,
                              target_names=target_names, key_map=key_map)
  cache_obj
}

################################################################################
# [L1] GROUP-SPECIFIC PEAK LOADING AND COUNTING
# Reads are counted in each group's own called peaks, not the union.
# This cache is used exclusively for Tier 1 condition-specific detection.
################################################################################

load_group_peak_counts <- function(group_peak_beds, bamfiles, sample_names) {
  cat("Loading and counting group-specific peak sets (L1)...\n")
  grp_counts <- list()
  grp_peaks  <- list()
  for (grp in names(group_peak_beds)) {
    bed <- group_peak_beds[[grp]]
    if (!file.exists(bed)) {
      cat("  WARNING:", grp, "GroupPeaks BED not found —", bed, "\n")
      cat("    Run Build_FullUniverse_AllGeno.sh first\n")
      grp_counts[[grp]] <- NULL
      grp_peaks[[grp]]  <- NULL
      next
    }
    peaks_gr <- import_bed(bed)
    cat(" ", grp, ":", length(peaks_gr), "peaks\n")
    cr <- count_reads_in_peaks(peaks_gr, bamfiles,
                               paste0("grp_", grp))
    cnts <- cr$counts
    # Align to valid peak ids
    if (nrow(cnts) < length(peaks_gr))
      peaks_gr <- peaks_gr[names(peaks_gr) %in% cr$valid_peak_ids]
    rownames(cnts)    <- names(peaks_gr)
    colnames(cnts)    <- sample_names
    grp_counts[[grp]] <- cnts
    grp_peaks[[grp]]  <- peaks_gr
  }
  list(counts=grp_counts, peaks=grp_peaks)
}

################################################################################
# [L2] CPM-SCALED TIER 1 THRESHOLDS
# Returns per-contrast thresholds in raw read space after scaling to library size.
################################################################################

compute_tier1_thresholds <- function(counts_all_samples,
                                     grpA_idx, grpB_idx,
                                     min_cpm = min_cpm_present,
                                     max_cpm_per_rep = max_cpm_absent) {
  # Use median library size of the two groups in this contrast
  contrast_cols <- c(grpA_idx, grpB_idx)
  lib_sizes     <- colSums(counts_all_samples[, contrast_cols, drop=FALSE])
  median_lib    <- median(lib_sizes)

  # min_reads: CPM × median_lib / 1e6, but must be ≥ 1
  min_reads <- max(1, round(min_cpm * median_lib / 1e6))

  # max_reads_other: per-replicate CPM ceiling summed across replicates
  # (so a group with n=2 has 2x the raw read headroom compared to n=1)
  n_A <- length(grpA_idx); n_B <- length(grpB_idx)
  # Absent condition = the one with fewer expected reads; take the smaller group
  n_absent <- min(n_A, n_B)
  max_reads <- max(1, round(max_cpm_per_rep * median_lib / 1e6 * n_absent))

  cat(sprintf("    Tier1 thresholds: min_reads=%d  max_reads_other=%d",
              min_reads, max_reads),
      sprintf("  (median_lib=%.0f, min_cpm=%.2f, max_cpm/rep=%.2f)\n",
              median_lib, min_cpm, max_cpm_per_rep))
  list(min_reads=min_reads, max_reads_other=max_reads)
}

################################################################################
# [L1 + L2] identify_condition_specific_pair
# Now operates on group-specific peak counts (grp_counts_A/B) rather than
# union counts, and uses CPM-scaled thresholds.
################################################################################

identify_condition_specific_pair <- function(
    # Group-specific peaks and counts (L1)
    grp_peaks_A,   grp_counts_A,   # peaks and counts for grpA's own peaks
    grp_peaks_B,   grp_counts_B,   # peaks and counts for grpB's own peaks
    # Union counts for CPM threshold scaling (L2)
    union_counts,
    grpA_idx, grpB_idx,
    grpA_label, grpB_label,
    min_reads, max_reads_other,
    mode_label, contrast_name) {

  # grpA_idx / grpB_idx are integer positions in the full coldata / union_counts.
  # grp_counts_A/B only contain columns for that genotype's BAMs, so we must
  # match by sample name rather than position to avoid out-of-bounds indexing.
  sA <- colnames(union_counts)[grpA_idx]
  sB <- colnames(union_counts)[grpB_idx]

  sA_in_grpA <- sA[sA %in% colnames(grp_counts_A)]
  sB_in_grpA <- sB[sB %in% colnames(grp_counts_A)]
  sB_in_grpB <- sB[sB %in% colnames(grp_counts_B)]
  sA_in_grpB <- sA[sA %in% colnames(grp_counts_B)]

  if (length(sA_in_grpA) == 0 && length(sB_in_grpA) == 0) {
    empty_df <- data.frame(peak_id=character(), chr=character(),
      start=integer(), end=integer(), log2FC=numeric(),
      pvalue=numeric(), padj=numeric(), baseMean=numeric(),
      AveExpr=numeric(), method=character(), contrast=character(),
      mode=character(), stringsAsFactors=FALSE)
    return(list(A_specific=empty_df, B_specific=empty_df))
  }

  A_in_A <- if (length(sA_in_grpA)>0)
    rowSums(grp_counts_A[, sA_in_grpA, drop=FALSE]) else rep(0, nrow(grp_counts_A))
  A_in_B <- if (length(sB_in_grpA)>0)
    rowSums(grp_counts_A[, sB_in_grpA, drop=FALSE]) else rep(0, nrow(grp_counts_A))
  B_in_B <- if (length(sB_in_grpB)>0)
    rowSums(grp_counts_B[, sB_in_grpB, drop=FALSE]) else rep(0, nrow(grp_counts_B))
  B_in_A <- if (length(sA_in_grpB)>0)
    rowSums(grp_counts_B[, sA_in_grpB, drop=FALSE]) else rep(0, nrow(grp_counts_B))

  contrast_bam_idx <- c(grpA_idx, grpB_idx)
  lib_size_est <- median(colSums(union_counts[, contrast_bam_idx, drop=FALSE])) + 1

  mk <- function(peaks_gr, present_counts, absent_counts,
                 n_present, lfc, label) {
    idx <- present_counts >= min_reads & absent_counts <= max_reads_other
    if (sum(idx) == 0) return(data.frame(
      peak_id=character(), chr=character(), start=integer(), end=integer(),
      log2FC=numeric(), pvalue=numeric(), padj=numeric(),
      baseMean=numeric(), AveExpr=numeric(),
      method=character(), contrast=character(), mode=character(),
      stringsAsFactors=FALSE))
    present_mean <- present_counts[idx] / max(n_present, 1)
    ave_expr_est <- log2((present_mean / lib_size_est * 1e6) + 0.5)
    data.frame(
      peak_id  = names(peaks_gr)[idx],
      chr      = as.character(seqnames(peaks_gr[idx])),
      start    = start(peaks_gr[idx]) - 1,
      end      = end(peaks_gr[idx]),
      log2FC   = lfc,
      pvalue   = 1e-10,
      padj     = 1e-10,
      baseMean = present_mean,
      AveExpr  = ave_expr_est,
      method   = paste0(mode_label, "_", label, "_specific_", contrast_name),
      contrast = contrast_name,
      mode     = mode_label,
      stringsAsFactors = FALSE)
  }

  list(
    A_specific = mk(grp_peaks_A, A_in_A, A_in_B,
                    length(grpA_idx), -5, paste0(grpA_label, "_open_lost")),
    B_specific = mk(grp_peaks_B, B_in_B, B_in_A,
                    length(grpB_idx),  5, paste0(grpB_label, "_open_gained"))
  )
}

################################################################################
# [L3] CONTRAST-SPECIFIC RUV CONTROL SELECTION
# Finds peaks stable between the two specific groups in this contrast,
# not across all 7 genotypes. Uses pairwise edgeR quasi-likelihood test.
################################################################################

select_ruv_controls_pairwise <- function(counts, grpA_idx, grpB_idx,
                                          grpA_label, grpB_label) {
  cat("    [L3] Selecting RUV controls for",
      grpA_label, "vs", grpB_label, "...\n")

  # Subset to only the two groups in this contrast
  pair_idx    <- c(grpA_idx, grpB_idx)
  counts_pair <- counts[, pair_idx, drop=FALSE]
  group_pair  <- factor(c(rep(grpA_label, length(grpA_idx)),
                          rep(grpB_label, length(grpB_idx))))

  peak_means  <- rowMeans(counts_pair)
  peak_cv     <- apply(counts_pair, 1, sd) / (peak_means + 0.5)
  exp_cut     <- quantile(peak_means, 0.4)
  cv_cut      <- quantile(peak_cv[peak_means > exp_cut], 0.3)

  # Quick pairwise edgeR to find non-DE peaks between these two groups only
  dge  <- DGEList(counts=counts_pair, group=group_pair)
  dge  <- calcNormFactors(dge)
  des  <- model.matrix(~ group_pair)
  dge  <- estimateDisp(dge, des)
  fit  <- glmQLFit(dge, des)
  test <- glmQLFTest(fit, coef=2)
  res  <- topTags(test, n=Inf)$table

  non_de_idx <- which(res$PValue > 0.5 & res$logCPM > 1)
  stable_idx <- which(peak_means > exp_cut & peak_cv < cv_cut)
  int_idx    <- intersect(stable_idx, non_de_idx)

  if (length(int_idx) >= 50)        { ctrl <- int_idx;    m <- "intersection" }
  else if (length(stable_idx) >= 50){ ctrl <- stable_idx; m <- "empirical_stable" }
  else if (length(non_de_idx) >= 50){ ctrl <- non_de_idx; m <- "non_DE" }
  else {
    n    <- min(200, max(100, round(nrow(counts_pair) * 0.1)))
    ctrl <- order(peak_cv[peak_means > quantile(peak_means, 0.2)])[seq_len(n)]
    m    <- "top_stable_CV"
  }
  cat("      Pairwise RUV controls:", length(ctrl), "(", m, ")\n")
  rownames(counts)[ctrl]
}

################################################################################
# run_one_contrast — uses pairwise RUV controls (L3)
# Accepts ruv_control_ids pre-computed per contrast
################################################################################

run_one_contrast <- function(counts, coldata, peaks_gr, mode_label,
                             contrast_name, contrast_vec, contrast_deseq2,
                             use_ruv=TRUE, ruv_k=1, ruv_control_ids=NULL) {
  cat("  Running", contrast_name, "—", mode_label, "...\n")
  counts <- round(counts)
  if (use_ruv) {
    valid_ctrl <- intersect(ruv_control_ids, rownames(counts))
    cat("    RUVg:", length(valid_ctrl), "control peaks\n")
    set     <- newSeqExpressionSet(
      counts     = as.matrix(counts),
      phenoData  = data.frame(condition=coldata$group, row.names=colnames(counts)))
    set_ruv <- RUVg(set, valid_ctrl, k=ruv_k)
    ruv_f   <- pData(set_ruv)[, grep("^W_", colnames(pData(set_ruv))), drop=FALSE]
    coldata_ruv <- cbind(coldata, ruv_f)
    counts_use  <- counts(set_ruv)
  } else {
    coldata_ruv <- coldata
    counts_use  <- counts
  }
  df <- if (use_ruv && ruv_k > 0)
    as.formula(paste("~ 0 + group +",
                     paste(paste0("W_", seq_len(ruv_k)), collapse=" + ")))
  else ~ 0 + group

  # ── DESeq2 ──────────────────────────────────────────────────────────────────
  cat("    DESeq2...\n")
  dds    <- DESeqDataSetFromMatrix(countData=counts_use,
                                   colData=coldata_ruv, design=df)
  dds    <- DESeq(dds, quiet=TRUE)
  res_d2 <- results(dds, contrast=contrast_deseq2)
  d2_df  <- data.frame(
    peak_id  = rownames(res_d2),
    chr      = as.character(seqnames(peaks_gr)),
    start    = start(peaks_gr) - 1, end = end(peaks_gr),
    log2FC   = res_d2$log2FoldChange,
    pvalue   = res_d2$pvalue,
    padj     = res_d2$padj,
    baseMean = res_d2$baseMean,
    method   = paste0(mode_label, "_DESeq2_", contrast_name),
    contrast = contrast_name, mode = mode_label,
    stringsAsFactors=FALSE)

  # ── limma-voom ──────────────────────────────────────────────────────────────
  cat("    limma-voom...\n")
  dge  <- DGEList(counts=counts_use, group=coldata$group)
  dge  <- calcNormFactors(dge, method="TMM")
  des  <- model.matrix(df, data=coldata_ruv)
  fcv  <- setNames(rep(0, ncol(des)), colnames(des))
  for (nm in names(contrast_vec)) {
    cn <- paste0("group", nm)
    if (cn %in% names(fcv)) fcv[cn] <- contrast_vec[nm]
    else if (nm %in% names(fcv)) fcv[nm] <- contrast_vec[nm]
    else warning("Contrast term '", nm, "' not found")
  }
  v    <- voom(dge, des, plot=FALSE)
  fit  <- lmFit(v, des)
  fit2 <- contrasts.fit(fit, contrasts=fcv)
  fit2 <- eBayes(fit2)
  res_lm <- topTable(fit2, coef=1, number=Inf, sort.by="none")
  lm_df  <- data.frame(
    peak_id  = rownames(res_lm),
    chr      = as.character(seqnames(peaks_gr)),
    start    = start(peaks_gr) - 1, end = end(peaks_gr),
    log2FC   = res_lm$logFC,
    pvalue   = res_lm$P.Value,
    padj     = res_lm$adj.P.Val,
    AveExpr  = res_lm$AveExpr,
    method   = paste0(mode_label, "_limma_", contrast_name),
    contrast = contrast_name, mode = mode_label,
    stringsAsFactors=FALSE)

  list(deseq2=d2_df, limma=lm_df)
}

filter_sig <- function(df, lfc, pv)
  df %>% filter(!is.na(pvalue), pvalue <= pv, abs(log2FC) >= lfc)

################################################################################
# [L5] CONTINUOUS CONFIDENCE SCORE
# conf_score = method_agreement × tier_score × −log10(p + eps)
# method_agreement encoding:
#   1.0 — single method (Tier 1 condition-specific, or only one stat method)
#   1.5 — two methods agree directionally (DESeq2 + limma, same resolution)
#   2.0 — all four calls agree: narrow DESeq2 + narrow limma + ext DESeq2 + ext limma
################################################################################

# Build a method agreement map for one contrast before calling df_to_gr
# Returns a named vector: peak_id → method_agreement score
compute_method_agreement <- function(nd2_df, nlm_df, bd2_df, blm_df) {
  # Collect all significant peak IDs and track which methods called them
  sig_peaks <- unique(c(nd2_df$peak_id, nlm_df$peak_id,
                        bd2_df$peak_id, blm_df$peak_id))
  if (length(sig_peaks) == 0) return(setNames(numeric(0), character(0)))

  in_nd2 <- sig_peaks %in% nd2_df$peak_id
  in_nlm <- sig_peaks %in% nlm_df$peak_id
  in_bd2 <- sig_peaks %in% bd2_df$peak_id
  in_blm <- sig_peaks %in% blm_df$peak_id

  n_methods <- in_nd2 + in_nlm + in_bd2 + in_blm

  agreement <- ifelse(n_methods == 4, 2.0,
               ifelse(n_methods >= 2, 1.5, 1.0))
  setNames(agreement, sig_peaks)
}

df_to_gr <- function(df, tier="differential",
                     method_agreement_map=NULL) {
  if (nrow(df) == 0) return(GRanges())
  eps            <- 1e-300
  tier_score_val <- switch(tier,
    "condition_specific" = 3,
    "differential"       = ifelse(grepl("^NARROW", df$mode[1]), 2, 1),
    1)

  # Method agreement: 1.0 for Tier 1 (only one detection method), else look up
  m_agree <- if (tier == "condition_specific") {
    rep(1.0, nrow(df))
  } else if (!is.null(method_agreement_map) && length(method_agreement_map) > 0) {
    ma <- method_agreement_map[df$peak_id]
    ma[is.na(ma)] <- 1.0
    as.numeric(ma)
  } else { rep(1.0, nrow(df)) }

  pval_use   <- ifelse(is.na(df$pvalue), 1, df$pvalue)
  conf_score <- m_agree * tier_score_val * (-log10(pval_use + eps))

  GRanges(
    seqnames    = df$chr,
    ranges      = IRanges(start=df$start+1, end=df$end),
    log2FC      = df$log2FC,
    pvalue      = pval_use,
    padj        = if ("padj"     %in% names(df)) df$padj     else NA_real_,
    baseMean    = if ("baseMean" %in% names(df)) df$baseMean else NA_real_,
    AveExpr     = if ("AveExpr"  %in% names(df)) df$AveExpr  else NA_real_,
    method      = df$method,
    contrast    = df$contrast,
    tier        = tier,
    tier_score  = tier_score_val,
    m_agree     = m_agree,
    conf_score  = conf_score,
    mode        = df$mode)
}

# Priority deduplication — now uses conf_score instead of raw tier_score
priority_deduplicate <- function(gr) {
  if (length(gr) <= 1) return(gr)
  cat("    Deduplicating", length(gr), "DAR entries...\n")

  mcols(gr)$sign_dir <- sign(mcols(gr)$log2FC)

  dedup_one <- function(g) {
    if (length(g) <= 1) return(g)
    ovl    <- findOverlaps(g, g, minoverlap=1L)
    ovl_df <- as.data.frame(ovl)
    ovl_df <- ovl_df[ovl_df$queryHits != ovl_df$subjectHits, ]
    if (nrow(ovl_df) == 0) return(g)
    parent <- seq_len(length(g))
    find_root <- function(x) {
      while (parent[x] != x) { parent[x] <<- parent[parent[x]]; x <- parent[x] }
      x
    }
    for (i in seq_len(nrow(ovl_df))) {
      a <- find_root(ovl_df$queryHits[i])
      b <- find_root(ovl_df$subjectHits[i])
      if (a != b) parent[a] <- b
    }
    clusters <- sapply(seq_len(length(g)), find_root)
    scores   <- mcols(g)$conf_score
    keep <- sapply(unique(clusters), function(cl) {
      m <- which(clusters == cl); m[which.max(scores[m])]
    })
    g[sort(keep)]
  }

  result <- c(dedup_one(gr[mcols(gr)$sign_dir > 0]),
              dedup_one(gr[mcols(gr)$sign_dir < 0]))
  cat("    After deduplication:", length(result), "unique DARs\n")
  result
}

################################################################################
# [L5] SCORED OVERLAP FUNCTIONS
# Replaces binary find_consistent_overlaps for background subtraction steps.
# Returns overlapping DARs weighted by min(conf_score) of each pair.
################################################################################

# Binary version retained for simple set operations
find_consistent_overlaps <- function(gr1, gr2) {
  if (length(gr1)==0||length(gr2)==0) return(GRanges())
  c(subsetByOverlaps(gr1[gr1$log2FC>0], gr2[gr2$log2FC>0], minoverlap=1L),
    subsetByOverlaps(gr1[gr1$log2FC<0], gr2[gr2$log2FC<0], minoverlap=1L))
}
find_opposite_overlaps <- function(gr1, gr2) {
  if (length(gr1)==0||length(gr2)==0) return(GRanges())
  c(subsetByOverlaps(gr1[gr1$log2FC>0], gr2[gr2$log2FC<0], minoverlap=1L),
    subsetByOverlaps(gr1[gr1$log2FC<0], gr2[gr2$log2FC>0], minoverlap=1L))
}

# Scored version — used for biologically meaningful intersections
# Returns gr1 entries that overlap gr2 in the same direction,
# with conf_score updated to min(conf_score_gr1, conf_score_gr2)
# and a cross_conf_score column added for downstream filtering.
scored_consistent_overlaps <- function(gr1, gr2,
                                        min_cross_conf=0) {
  if (length(gr1)==0 || length(gr2)==0) return(GRanges())

  # Ensure conf_score exists
  if (!"conf_score" %in% names(mcols(gr1)))
    mcols(gr1)$conf_score <- mcols(gr1)$tier_score * (-log10(mcols(gr1)$pvalue + 1e-300))
  if (!"conf_score" %in% names(mcols(gr2)))
    mcols(gr2)$conf_score <- mcols(gr2)$tier_score * (-log10(mcols(gr2)$pvalue + 1e-300))

  process_dir <- function(g1_sub, g2_sub) {
    if (length(g1_sub)==0||length(g2_sub)==0) return(GRanges())
    ovl <- findOverlaps(g1_sub, g2_sub, minoverlap=1L)
    if (length(ovl)==0) return(GRanges())
    q_idx <- queryHits(ovl); s_idx <- subjectHits(ovl)
    # For each gr1 hit, take the max conf_score partner in gr2
    cross_scores <- tapply(
      mcols(g2_sub)$conf_score[s_idx],
      q_idx,
      max)
    hit_idx   <- as.integer(names(cross_scores))
    result_gr <- g1_sub[hit_idx]
    mcols(result_gr)$cross_conf_score <- as.numeric(cross_scores)
    # Cross confidence = geometric mean of the two scores
    mcols(result_gr)$conf_score <-
      sqrt(mcols(result_gr)$conf_score * mcols(result_gr)$cross_conf_score)
    result_gr[mcols(result_gr)$cross_conf_score >= min_cross_conf]
  }

  c(process_dir(gr1[gr1$log2FC>0], gr2[gr2$log2FC>0]),
    process_dir(gr1[gr1$log2FC<0], gr2[gr2$log2FC<0]))
}

################################################################################
# [L7] TIER 3 DIAGNOSTIC
# Classifies Tier 3-only DARs by what drove the extended-peak signal.
################################################################################

diagnose_tier3_dars <- function(contrast_name,
                                 all_gr,          # full deduplicated GRanges
                                 counts_narrow,   # narrow counts matrix
                                 counts_broad,    # extended counts matrix
                                 peaks_narrow,    # narrow GRanges
                                 peaks_broad,     # extended GRanges
                                 grp_peaks_A,     # group-specific peaks grpA
                                 grp_peaks_B,     # group-specific peaks grpB
                                 ext_bp,
                                 out_dir) {

  tier3_gr <- all_gr[mcols(all_gr)$tier == "differential" &
                     grepl("EXTENDED", mcols(all_gr)$mode)]

  cat("  [L7] Tier 3 diagnostic:", length(tier3_gr), "extended-only DARs\n")
  if (length(tier3_gr) == 0) return(invisible(NULL))

  # (a) Fraction of extended region covered by any narrow peak
  ovl_narrow <- findOverlaps(tier3_gr, peaks_narrow, minoverlap=1L)
  covered <- tabulate(queryHits(ovl_narrow), nbins=length(tier3_gr)) > 0

  # (b) Narrow peak in grpA-specific peaks?
  in_A <- if (!is.null(grp_peaks_A))
    countOverlaps(tier3_gr, grp_peaks_A, minoverlap=1L) > 0 else rep(FALSE, length(tier3_gr))
  in_B <- if (!is.null(grp_peaks_B))
    countOverlaps(tier3_gr, grp_peaks_B, minoverlap=1L) > 0 else rep(FALSE, length(tier3_gr))

  # (c) Adjacent narrow peak within ext_bp (would merge after extension)?
  flanked <- GRanges(
    seqnames = seqnames(tier3_gr),
    ranges   = IRanges(
      start = pmax(1L, start(tier3_gr) - ext_bp),
      end   = end(tier3_gr) + ext_bp),
    strand = "*")
  n_adjacent <- countOverlaps(flanked, peaks_narrow, minoverlap=1L)

  # Classify
  classification <- dplyr::case_when(
    covered & (in_A | in_B) ~ "tight_boundary",   # narrow peak exists; extension helped
    n_adjacent >= 2          ~ "adjacent_merged",  # two nearby peaks merged
    covered                  ~ "tight_boundary",
    TRUE                     ~ "diffuse"            # genuinely broad accessibility
  )

  diag_df <- data.frame(
    chr            = as.character(seqnames(tier3_gr)),
    start          = start(tier3_gr) - 1,
    end            = end(tier3_gr),
    log2FC         = mcols(tier3_gr)$log2FC,
    pvalue         = mcols(tier3_gr)$pvalue,
    conf_score     = mcols(tier3_gr)$conf_score,
    direction      = ifelse(mcols(tier3_gr)$log2FC > 0, "gained-open", "gained-close"),
    has_narrow_peak      = covered,
    in_grpA_peaks        = in_A,
    in_grpB_peaks        = in_B,
    n_adjacent_narrow    = n_adjacent,
    classification       = classification,
    stringsAsFactors = FALSE)

  diag_out <- file.path(out_dir,
    paste0(contrast_name, "_Tier3_diagnostic.txt"))
  write.table(diag_df, diag_out, sep="\t", quote=FALSE,
              row.names=FALSE, col.names=TRUE)

  tbl <- table(classification)
  cat("    Tier3 classification:", paste(names(tbl), tbl, sep="=", collapse=" | "), "\n")
  cat("    →", diag_out, "\n")
  invisible(diag_df)
}

################################################################################
# OUTPUT WRITERS
################################################################################

write_dar_outputs <- function(gr, tag, out_dir) {
  if (length(gr) == 0) { cat("  Skipping", tag, "(0 DARs)\n"); return(invisible(NULL)) }
  mc <- mcols(gr)
  df <- data.frame(
    chr        = as.character(seqnames(gr)),
    start      = start(gr) - 1, end = end(gr),
    name       = paste0(tag, "_DAR_", seq_along(gr)),
    score      = round(-log10(mc$pvalue + 1e-300) * 10),
    strand     = ".",
    log2FC     = mc$log2FC,
    pvalue     = mc$pvalue,
    padj       = if ("padj"       %in% names(mc)) mc$padj       else NA_real_,
    baseMean   = if ("baseMean"   %in% names(mc)) mc$baseMean   else NA_real_,
    AveExpr    = if ("AveExpr"    %in% names(mc)) mc$AveExpr    else NA_real_,
    conf_score = if ("conf_score" %in% names(mc)) mc$conf_score else NA_real_,
    m_agree    = if ("m_agree"    %in% names(mc)) mc$m_agree    else NA_real_,
    direction  = ifelse(mc$log2FC > 0, "gained-open", "gained-close"),
    tier       = if ("tier"       %in% names(mc)) mc$tier       else NA_character_,
    tier_score = if ("tier_score" %in% names(mc)) mc$tier_score else NA_real_,
    mode       = if ("mode"       %in% names(mc)) mc$mode       else NA_character_,
    method     = if ("method"     %in% names(mc)) mc$method     else NA_character_,
    contrast   = if ("contrast"   %in% names(mc)) mc$contrast   else NA_character_,
    stringsAsFactors = FALSE)
  write.table(df[, c("chr","start","end","name","score","strand")],
    file.path(out_dir, paste0(tag, "_DARs.bed")),
    sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(df,
    file.path(out_dir, paste0(tag, "_DARs_annotated.txt")),
    sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  cat("  Saved", nrow(df), "DARs →", paste0(tag, "_DARs.bed\n"))
  invisible(df)
}

write_tier1_outputs <- function(cs_narrow, cs_ext, contrast_name, out_dir) {
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  write_one <- function(df, tag) {
    if (nrow(df) == 0) { cat("  T1 skip:", tag, "(0)\n"); return(invisible(NULL)) }
    df$direction <- ifelse(df$log2FC > 0, "gained-open", "gained-close")
    df$name <- paste0(tag, "_", seq_len(nrow(df))); df$strand <- "."
    write.table(df[, c("chr","start","end","name","pvalue","strand")],
      file.path(out_dir, paste0(tag, "_DARs.bed")),
      sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(df,
      file.path(out_dir, paste0(tag, "_DARs_annotated.txt")),
      sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
    cat("  T1 saved", nrow(df), "→", paste0(tag, "_DARs.bed\n"))
  }
  A_lost <- rbind(
    cs_narrow$A_specific[cs_narrow$A_specific$log2FC < 0, , drop=FALSE],
    cs_ext$A_specific   [cs_ext$A_specific$log2FC    < 0, , drop=FALSE])
  A_lost <- A_lost[!duplicated(paste(A_lost$chr, A_lost$start, A_lost$end)), ]
  B_gained <- rbind(
    cs_narrow$B_specific[cs_narrow$B_specific$log2FC > 0, , drop=FALSE],
    cs_ext$B_specific   [cs_ext$B_specific$log2FC    > 0, , drop=FALSE])
  B_gained <- B_gained[!duplicated(paste(B_gained$chr, B_gained$start, B_gained$end)), ]
  write_one(A_lost,   paste0("TIER1_", contrast_name, "_A_open_lost"))
  write_one(B_gained, paste0("TIER1_", contrast_name, "_B_open_gained"))
  cat(sprintf("  TIER1 %s: A_lost=%d  B_gained=%d\n",
              contrast_name, nrow(A_lost), nrow(B_gained)))
}

write_method_outputs <- function(narrow_diff, broad_diff, contrast_name, out_dir) {
  mdir <- file.path(out_dir, "Method_Split")
  dir.create(mdir, recursive=TRUE, showWarnings=FALSE)
  fmt <- function(n, b) {
    cm <- rbind(n, b); if (nrow(cm) == 0) return(NULL)
    out <- data.frame(
      seqnames  = cm$chr, start=cm$start, end=cm$end,
      log2FC    = cm$log2FC, pvalue=cm$pvalue, padj=cm$padj,
      baseMean  = if ("baseMean" %in% names(cm)) cm$baseMean else NA_real_,
      AveExpr   = if ("AveExpr"  %in% names(cm)) cm$AveExpr  else NA_real_,
      direction = ifelse(cm$log2FC > 0, "gained-open", "gained-close"),
      mode=cm$mode, method=cm$method, contrast=cm$contrast,
      stringsAsFactors=FALSE)
    out <- out[order(out$pvalue, na.last=TRUE), ]
    out[!duplicated(paste(out$seqnames, out$start, out$end)), ]
  }
  for (meth in c("deseq2", "limma")) {
    ml  <- ifelse(meth == "deseq2", "DESeq2", "limma")
    all <- fmt(narrow_diff[[meth]], broad_diff[[meth]])
    if (!is.null(all) && nrow(all) > 0) {
      write.table(all,
        file.path(mdir, paste0(ml, "_RUV_", contrast_name, "_DARs.txt")),
        sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
      sig <- all[!is.na(all$pvalue) & all$pvalue <= pvalue_threshold &
                 abs(all$log2FC) >= log2fc_threshold, ]
      write.table(sig,
        file.path(mdir, paste0(ml, "_RUV_", contrast_name, "_sig_DARs.txt")),
        sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
      cat(" ", ml, ":", contrast_name, "— all:", nrow(all), "| sig:", nrow(sig), "\n")
    }
  }
}

################################################################################
# MAIN ANALYSIS
################################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("SEVEN-GENOTYPE nc14b DAR ANALYSIS — v4\n")
cat("L1:group-specific Tier1 | L2:CPM thresholds | L3:pairwise RUV |\n")
cat("L4:boundary extension  | L5:confidence score | L7:Tier3 diagnostic\n")
cat("v4: all counts loaded from cache (Precount_AllCounts_Matrices.r)\n")
cat(strrep("=", 80), "\n\n")

cat("Sample layout:\n"); print(coldata[, c("genotype","group")]); cat("\n")

# ---------------------------------------------------------------------------
# LOAD UNION PEAK UNIVERSE
# ---------------------------------------------------------------------------
cat("Loading union peak universe...\n")
existing_peaks <- import_bed(full_universe_bed)
if (is.null(names(existing_peaks)))
  names(existing_peaks) <- paste0("peak_", seq_along(existing_peaks))
cat("  Total peaks:", length(existing_peaks), "\n\n")

# ---------------------------------------------------------------------------
# NARROW COUNTS — load from cache or count on first run
# ---------------------------------------------------------------------------
narrow_cache_file <- file.path(CACHE_DIR, "narrow_counts.txt")
# Also accept the path the Precount script writes to directly
direct_narrow_file <- file.path(output_dir, "SevenGeno_nc14b_counts_matrix.txt")

if (file.exists(narrow_cache_file)) {
  cat("Loading narrow counts from cache...\n")
  counts_narrow <- as.matrix(read.table(narrow_cache_file,
                                         header=TRUE, row.names=1, sep="\t"))
} else if (file.exists(direct_narrow_file)) {
  cat("Loading narrow counts (output dir copy)...\n")
  counts_narrow <- as.matrix(read.table(direct_narrow_file,
                                         header=TRUE, row.names=1, sep="\t"))
} else {
  cat("Cache not found — running featureCounts for narrow peaks (first run).\n")
  cat("Run Precount_AllCounts_Matrices.r first to avoid this.\n")
  cr            <- count_reads_in_peaks(existing_peaks, bamfiles, "peak")
  counts_narrow <- cr$counts
  if (nrow(counts_narrow) < length(existing_peaks))
    existing_peaks <- existing_peaks[names(existing_peaks) %in% cr$valid_peak_ids]
  rownames(counts_narrow) <- names(existing_peaks)
  dir.create(CACHE_DIR, recursive=TRUE, showWarnings=FALSE)
  write.table(counts_narrow, narrow_cache_file,
              sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
  write.table(counts_narrow, direct_narrow_file,
              sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
  cat("Counts saved for future runs.\n\n")
}
counts_narrow <- realign_cache_cols(counts_narrow, sample_names)
if (is.null(counts_narrow))
  stop("Could not align narrow counts cache columns to sample_names. ",
       "Check that the cache was built with the same BAMs as this script.")

# Sync existing_peaks to counts rows
if (!all(rownames(counts_narrow) %in% names(existing_peaks))) {
  existing_peaks <- existing_peaks[names(existing_peaks) %in% rownames(counts_narrow)]
}
existing_peaks <- existing_peaks[rownames(counts_narrow)]

keep_n             <- rowSums(counts_narrow >= 5) >= 2
counts_narrow_filt <- counts_narrow[keep_n, ]
peaks_narrow_filt  <- existing_peaks[keep_n]
cat("Narrow peaks after filtering:", length(peaks_narrow_filt), "\n\n")

# ---------------------------------------------------------------------------
# [L1] GROUP-SPECIFIC PEAK COUNTS — load from cache or count on first run
# ---------------------------------------------------------------------------
cat("=== L1: Loading group-specific peak counts ===\n")
grp_cache_file <- file.path(CACHE_DIR, "grp_counts_cache.rds")
if (file.exists(grp_cache_file)) {
  cat("  Loading from cache:", grp_cache_file, "\n")
  grp_cache <- readRDS(grp_cache_file)
  grp_cache <- realign_grp_cache(grp_cache, sample_names)
} else {
  cat("  Cache not found — counting (first run). Run Precount script first.\n")
  grp_cache <- load_group_peak_counts(group_peak_beds, bamfiles, sample_names)
  dir.create(CACHE_DIR, recursive=TRUE, showWarnings=FALSE)
  saveRDS(grp_cache, grp_cache_file)
  cat("  Saved:", grp_cache_file, "\n")
}
cat("\n")

# ---------------------------------------------------------------------------
# [L4] EXTENDED UNION PEAK COUNTS — load from cache or count on first run
# ---------------------------------------------------------------------------
cat("=== L4: Loading boundary-extended peak counts ===\n")
ext_cache <- list()
for (ebp in sort(unique(c(extension_bp_global, 300L)))) {
  key        <- as.character(ebp)
  cache_file <- file.path(CACHE_DIR, paste0("ext", ebp, "_counts.rds"))
  if (file.exists(cache_file)) {
    cat("  Loading ext", ebp, "bp from cache\n")
    cached <- readRDS(cache_file)
    cached$counts <- realign_cache_cols(cached$counts, sample_names)
    ext_cache[[key]] <- cached
  } else {
    cat("  Cache missing for ext", ebp, "bp — counting (first run).\n")
    p_ext <- build_extended_peaks_boundary(existing_peaks, ebp)
    cr    <- count_reads_in_peaks(p_ext, bamfiles, paste0("bext", ebp))
    cnts  <- cr$counts
    if (nrow(cnts) < length(p_ext))
      p_ext <- p_ext[names(p_ext) %in% cr$valid_peak_ids]
    rownames(cnts) <- names(p_ext)
    cnts  <- realign_cache_cols(cnts, sample_names)
    keep  <- rowSums(cnts >= 5) >= 2
    ext_cache[[key]] <- list(counts=cnts[keep, ], peaks=p_ext[keep])
    dir.create(CACHE_DIR, recursive=TRUE, showWarnings=FALSE)
    saveRDS(ext_cache[[key]], cache_file)
    cat("  Saved:", cache_file, "\n")
  }
}
cat("\n")

# ---------------------------------------------------------------------------
# [L1+L4] EXTENDED GROUP-SPECIFIC PEAK COUNTS — load from cache
# Eliminates ~19 redundant featureCounts calls from within run_contrast().
# ---------------------------------------------------------------------------
cat("=== L1+L4: Loading extended group-specific peak counts ===\n")
grp_ext_cache <- list()
for (ebp in sort(unique(c(extension_bp_global, 300L)))) {
  key        <- as.character(ebp)
  cache_file <- file.path(CACHE_DIR, paste0("grp_ext", ebp, "_counts.rds"))
  if (file.exists(cache_file)) {
    cat("  Loading grp_ext", ebp, "bp from cache\n")
    grp_ext_cache[[key]] <- readRDS(cache_file)
    grp_ext_cache[[key]] <- realign_grp_cache(grp_ext_cache[[key]], sample_names)
  } else {
    cat("  Cache missing for grp_ext", ebp, "bp — counting (first run).\n")
    grp_counts_ext <- list(); grp_peaks_ext <- list()
    for (grp in names(group_peak_beds)) {
      bed <- group_peak_beds[[grp]]
      if (!file.exists(bed)) next
      gp     <- import_bed(bed)
      gp_ext <- build_extended_peaks_boundary(gp, ebp)
      cr     <- count_reads_in_peaks(gp_ext, bamfiles, paste0("grpext", ebp, "_", grp))
      cnts   <- cr$counts; colnames(cnts) <- sample_names
      gp_ext <- gp_ext[names(gp_ext) %in% cr$valid_peak_ids]
      rownames(cnts)          <- names(gp_ext)
      grp_counts_ext[[grp]]   <- cnts
      grp_peaks_ext[[grp]]    <- gp_ext
    }
    grp_ext_cache[[key]] <- list(counts=grp_counts_ext, peaks=grp_peaks_ext)
    dir.create(CACHE_DIR, recursive=TRUE, showWarnings=FALSE)
    saveRDS(grp_ext_cache[[key]], cache_file)
    cat("  Saved:", cache_file, "\n")
  }
}
cat("\n")

# ---------------------------------------------------------------------------
# DEFINE ALL 13 PAIRWISE CONTRASTS
# ---------------------------------------------------------------------------

idx_BOTv    <- which(coldata$group == "BOTv")
idx_BOTCv   <- which(coldata$group == "BOTCv")
idx_BOT     <- which(coldata$group == "BOT")
idx_BOT_hR  <- which(coldata$group == "BOT_hR")
idx_BOTR    <- which(coldata$group == "BOTR")
idx_BOTC    <- which(coldata$group == "BOTC")
idx_BOTC_oR <- which(coldata$group == "BOTC_oR")

comparisons <- list(
  list(name="BOTv_vs_BOTR",
       label="FoxL1 (Cic+,Runt GoF,vent) vs Runt-null",
       contrast_vec=c(BOTv=1,BOTR=-1), deseq2=c("group","BOTv","BOTR"),
       grpA_label="BOTv", grpB_label="BOTR",
       grpA_idx=idx_BOTv, grpB_idx=idx_BOTR),

  list(name="BOTCv_vs_BOTR",
       label="HLH54F-vent (Cic-,vent) vs Runt-null",
       contrast_vec=c(BOTCv=1,BOTR=-1), deseq2=c("group","BOTCv","BOTR"),
       grpA_label="BOTCv", grpB_label="BOTR",
       grpA_idx=idx_BOTCv, grpB_idx=idx_BOTR),

  list(name="BOTv_vs_BOTCv",
       label="FoxL1 vs HLH54F-vent (within tolrm9 background)",
       contrast_vec=c(BOTv=1,BOTCv=-1), deseq2=c("group","BOTv","BOTCv"),
       grpA_label="BOTv", grpB_label="BOTCv",
       grpA_idx=idx_BOTv, grpB_idx=idx_BOTCv),

  list(name="BOTCv_vs_BOTC",
       label="Cic-del vent vs non-vent — EMPIRICAL tolrm9",
       contrast_vec=c(BOTCv=1,BOTC=-1), deseq2=c("group","BOTCv","BOTC"),
       grpA_label="BOTCv", grpB_label="BOTC",
       grpA_idx=idx_BOTCv, grpB_idx=idx_BOTC),

  list(name="BOTC_vs_BOTR",
       label="Cic-deleted vs Runt-null, both non-vent",
       contrast_vec=c(BOTC=1,BOTR=-1), deseq2=c("group","BOTC","BOTR"),
       grpA_label="BOTC", grpB_label="BOTR",
       grpA_idx=idx_BOTC, grpB_idx=idx_BOTR),

  list(name="BOT_hR_vs_BOTR",
       label="Runt het vs null (1 copy vs 0; Cic intact, non-vent)",
       contrast_vec=c(BOT_hR=1,BOTR=-1), deseq2=c("group","BOT_hR","BOTR"),
       grpA_label="BOT_hR", grpB_label="BOTR",
       grpA_idx=idx_BOT_hR, grpB_idx=idx_BOTR),

  list(name="BOT_hR_vs_BOTCv",
       label="Runt het (non-vent) vs HLH54F-vent",
       contrast_vec=c(BOT_hR=1,BOTCv=-1), deseq2=c("group","BOT_hR","BOTCv"),
       grpA_label="BOT_hR", grpB_label="BOTCv",
       grpA_idx=idx_BOT_hR, grpB_idx=idx_BOTCv),

  list(name="BOTv_vs_BOT_hR",
       label="Runt GoF (vent) vs Runt het (non-vent)",
       contrast_vec=c(BOTv=1,BOT_hR=-1), deseq2=c("group","BOTv","BOT_hR"),
       grpA_label="BOTv", grpB_label="BOT_hR",
       grpA_idx=idx_BOTv, grpB_idx=idx_BOT_hR),

  list(name="BOT_vs_BOTR",
       label="Runt WT vs Runt null — CLEANEST; same Cic+ non-vent background",
       contrast_vec=c(BOT=1,BOTR=-1), deseq2=c("group","BOT","BOTR"),
       grpA_label="BOT", grpB_label="BOTR",
       grpA_idx=idx_BOT, grpB_idx=idx_BOTR,
       extension_bp=300),

  list(name="BOT_vs_BOT_hR",
       label="Runt WT vs Runt het (dosage step; Cic intact, non-vent)",
       contrast_vec=c(BOT=1,BOT_hR=-1), deseq2=c("group","BOT","BOT_hR"),
       grpA_label="BOT", grpB_label="BOT_hR",
       grpA_idx=idx_BOT, grpB_idx=idx_BOT_hR),

  list(name="BOT_vs_BOTC",
       label="Runt+ Cic+ vs Cic-deleted (non-vent; Cic effect with Runt intact)",
       contrast_vec=c(BOT=1,BOTC=-1), deseq2=c("group","BOT","BOTC"),
       grpA_label="BOT", grpB_label="BOTC",
       grpA_idx=idx_BOT, grpB_idx=idx_BOTC),

  list(name="BOTC_oR_vs_BOTC",
       label="Runt overexpression rescue in Cic-deleted context",
       contrast_vec=c(BOTC_oR=1,BOTC=-1), deseq2=c("group","BOTC_oR","BOTC"),
       grpA_label="BOTC_oR", grpB_label="BOTC",
       grpA_idx=idx_BOTC_oR, grpB_idx=idx_BOTC),

  list(name="BOTC_oR_vs_BOTR",
       label="oRunt (Cic-del) vs Runt-null (Cic+): cross-background validation",
       contrast_vec=c(BOTC_oR=1,BOTR=-1), deseq2=c("group","BOTC_oR","BOTR"),
       grpA_label="BOTC_oR", grpB_label="BOTR",
       grpA_idx=idx_BOTC_oR, grpB_idx=idx_BOTR)
)

################################################################################
# RUN ALL 13 CONTRASTS
################################################################################

cat(strrep("=", 80), "\n")
cat("RUNNING 13 CONTRASTS — L1-L5+L7 active\n")
cat(strrep("=", 80), "\n\n")

dar_list <- list()

run_contrast <- function(cmp) {
  cat(strrep("-", 60), "\n")
  cat("CONTRAST:", cmp$name, "\n  (", cmp$label, ")\n\n")

  ebp    <- if (!is.null(cmp$extension_bp)) cmp$extension_bp else extension_bp_global
  cached <- ext_cache[[as.character(ebp)]]
  cat("  extension_bp =", ebp, "(boundary-based, L4)\n")

  # ── [L2] Compute CPM-scaled Tier 1 thresholds ─────────────────────────────
  t1_thr <- compute_tier1_thresholds(
    counts_narrow_filt, cmp$grpA_idx, cmp$grpB_idx)
  min_r  <- t1_thr$min_reads
  max_r  <- t1_thr$max_reads_other

  # ── [L1] Tier 1 with group-specific peaks ─────────────────────────────────
  grp_A_peaks  <- grp_cache$peaks[[cmp$grpA_label]]
  grp_A_counts <- grp_cache$counts[[cmp$grpA_label]]
  grp_B_peaks  <- grp_cache$peaks[[cmp$grpB_label]]
  grp_B_counts <- grp_cache$counts[[cmp$grpB_label]]

  if (is.null(grp_A_counts) || is.null(grp_B_counts)) {
    cat("  WARNING: group-specific peaks missing — falling back to union peaks for Tier 1\n")
    # Fallback: use union counts subset (original behaviour)
    nc <- list(
      A_specific = data.frame(peak_id=character(),chr=character(),
        start=integer(),end=integer(),log2FC=numeric(),pvalue=numeric(),
        padj=numeric(),baseMean=numeric(),AveExpr=numeric(),
        method=character(),contrast=character(),mode=character(),
        stringsAsFactors=FALSE),
      B_specific = data.frame(peak_id=character(),chr=character(),
        start=integer(),end=integer(),log2FC=numeric(),pvalue=numeric(),
        padj=numeric(),baseMean=numeric(),AveExpr=numeric(),
        method=character(),contrast=character(),mode=character(),
        stringsAsFactors=FALSE))
    bc <- nc
  } else {
    nc <- identify_condition_specific_pair(
      grp_A_peaks, grp_A_counts, grp_B_peaks, grp_B_counts,
      counts_narrow_filt,
      cmp$grpA_idx, cmp$grpB_idx,
      cmp$grpA_label, cmp$grpB_label,
      min_r, max_r, "NARROW", cmp$name)

    # Extended version of Tier 1: extended group-specific peaks
    # Build on-the-fly from group peaks using boundary extension
    # ── [L1+L4] Extended group-specific peaks — loaded from grp_ext_cache ─────
    # grp_ext_cache is pre-populated before the contrast loop, eliminating
    # 2 featureCounts calls per contrast (~26 calls saved across 13 contrasts).
    ebp_key    <- as.character(ebp)
    grp_A_ext  <- grp_ext_cache[[ebp_key]]$peaks[[cmp$grpA_label]]
    cnts_A_ext <- grp_ext_cache[[ebp_key]]$counts[[cmp$grpA_label]]
    grp_B_ext  <- grp_ext_cache[[ebp_key]]$peaks[[cmp$grpB_label]]
    cnts_B_ext <- grp_ext_cache[[ebp_key]]$counts[[cmp$grpB_label]]

    if (is.null(cnts_A_ext) || is.null(cnts_B_ext)) {
      cat("  WARNING: grp_ext cache missing for", cmp$grpA_label, "/",
          cmp$grpB_label, "— falling back to featureCounts\n")
      grp_A_ext <- build_extended_peaks_boundary(grp_A_peaks, ebp)
      grp_B_ext <- build_extended_peaks_boundary(grp_B_peaks, ebp)
      cr_A_ext  <- count_reads_in_peaks(grp_A_ext, bamfiles,
                                         paste0("grpext_", cmp$grpA_label))
      cr_B_ext  <- count_reads_in_peaks(grp_B_ext, bamfiles,
                                         paste0("grpext_", cmp$grpB_label))
      cnts_A_ext <- cr_A_ext$counts
      cnts_B_ext <- cr_B_ext$counts
      colnames(cnts_A_ext) <- sample_names
      colnames(cnts_B_ext) <- sample_names
      grp_A_ext <- grp_A_ext[names(grp_A_ext) %in% cr_A_ext$valid_peak_ids]
      grp_B_ext <- grp_B_ext[names(grp_B_ext) %in% cr_B_ext$valid_peak_ids]
    }

    bc <- identify_condition_specific_pair(
      grp_A_ext, cnts_A_ext, grp_B_ext, cnts_B_ext,
      counts_narrow_filt,
      cmp$grpA_idx, cmp$grpB_idx,
      cmp$grpA_label, cmp$grpB_label,
      min_r, max_r,
      paste0("EXTENDED", ebp), cmp$name)
  }

  cat("  Tier1 (narrow group peaks):", nrow(nc$A_specific), "+",
      nrow(nc$B_specific), "\n")
  cat("  Tier1 (ext", ebp, "bp group peaks):", nrow(bc$A_specific), "+",
      nrow(bc$B_specific), "\n")

  write_tier1_outputs(nc, bc, cmp$name, tier1_dir)

  # ── [L3] Contrast-specific RUV controls ───────────────────────────────────
  ruv_ctrl_narrow_pair <- select_ruv_controls_pairwise(
    counts_narrow_filt, cmp$grpA_idx, cmp$grpB_idx,
    cmp$grpA_label, cmp$grpB_label)

  ruv_ctrl_broad_pair  <- select_ruv_controls_pairwise(
    cached$counts, cmp$grpA_idx, cmp$grpB_idx,
    cmp$grpA_label, cmp$grpB_label)

  # ── Tier 2/3: differential analysis ────────────────────────────────────────
  nd <- run_one_contrast(counts_narrow_filt, coldata, peaks_narrow_filt,
          "NARROW", cmp$name, cmp$contrast_vec, cmp$deseq2,
          use_ruv, ruv_k, ruv_ctrl_narrow_pair)
  bd <- run_one_contrast(cached$counts, coldata, cached$peaks,
          paste0("EXTENDED", ebp), cmp$name,
          cmp$contrast_vec, cmp$deseq2,
          use_ruv, ruv_k, ruv_ctrl_broad_pair)

  nd2 <- filter_sig(nd$deseq2, log2fc_threshold, pvalue_threshold)
  nlm <- filter_sig(nd$limma,  log2fc_threshold, pvalue_threshold)
  bd2 <- filter_sig(bd$deseq2, log2fc_threshold, pvalue_threshold)
  blm <- filter_sig(bd$limma,  log2fc_threshold, pvalue_threshold)

  cat("  Narrow  DESeq2:", nrow(nd2), " limma:", nrow(nlm), "\n")
  cat("  Extended DESeq2:", nrow(bd2), " limma:", nrow(blm), "\n")

  write_method_outputs(nd, bd, cmp$name, output_dir)

  # ── [L5] Method agreement map ─────────────────────────────────────────────
  ma_map <- compute_method_agreement(nd2, nlm, bd2, blm)

  # ── Assemble with confidence scores ───────────────────────────────────────
  all_gr <- c(
    df_to_gr(nc$A_specific, "condition_specific"),
    df_to_gr(nc$B_specific, "condition_specific"),
    df_to_gr(bc$A_specific, "condition_specific"),
    df_to_gr(bc$B_specific, "condition_specific"),
    df_to_gr(nd2, "differential", ma_map),
    df_to_gr(nlm, "differential", ma_map),
    df_to_gr(bd2, "differential", ma_map),
    df_to_gr(blm, "differential", ma_map))

  if (length(all_gr) > 1) all_gr <- priority_deduplicate(all_gr)
  mcols(all_gr)$direction <- ifelse(mcols(all_gr)$log2FC > 0,
                                    "gained-open", "gained-close")

  # ── [L7] Tier 3 diagnostic ────────────────────────────────────────────────
  diagnose_tier3_dars(
    contrast_name = cmp$name,
    all_gr        = all_gr,
    counts_narrow = counts_narrow_filt,
    counts_broad  = cached$counts,
    peaks_narrow  = peaks_narrow_filt,
    peaks_broad   = cached$peaks,
    grp_peaks_A   = grp_cache$peaks[[cmp$grpA_label]],
    grp_peaks_B   = grp_cache$peaks[[cmp$grpB_label]],
    ext_bp        = ebp,
    out_dir       = diag_dir)

  cat("  Total DARs (after dedup):", length(all_gr), "\n\n")
  all_gr
}

for (cmp in comparisons) {
  dar_list[[cmp$name]] <- run_contrast(cmp)
}

################################################################################
# BACKGROUND AND SPECIFICITY ANALYSIS
# Uses scored_consistent_overlaps (L5) for all biologically meaningful
# intersections; binary find_consistent_overlaps retained for simple set ops
################################################################################

cat(strrep("=", 80), "\n")
cat("BACKGROUND AND SPECIFICITY ANALYSIS — scored overlaps (L5)\n")
cat(strrep("=", 80), "\n\n")

gr_BOTv_BOTR    <- dar_list[["BOTv_vs_BOTR"]]
gr_BOTCv_BOTR   <- dar_list[["BOTCv_vs_BOTR"]]
gr_BOTv_BOTCv   <- dar_list[["BOTv_vs_BOTCv"]]
gr_BOTCv_BOTC   <- dar_list[["BOTCv_vs_BOTC"]]
gr_BOTC_BOTR    <- dar_list[["BOTC_vs_BOTR"]]
gr_BOT_hR_BOTR  <- dar_list[["BOT_hR_vs_BOTR"]]
gr_BOT_hR_BOTCv <- dar_list[["BOT_hR_vs_BOTCv"]]
gr_BOTv_BOT_hR  <- dar_list[["BOTv_vs_BOT_hR"]]
gr_BOT_BOTR     <- dar_list[["BOT_vs_BOTR"]]
gr_BOT_BOT_hR   <- dar_list[["BOT_vs_BOT_hR"]]
gr_BOT_BOTC     <- dar_list[["BOT_vs_BOTC"]]
gr_BOTCoR_BOTC  <- dar_list[["BOTC_oR_vs_BOTC"]]
gr_BOTCoR_BOTR  <- dar_list[["BOTC_oR_vs_BOTR"]]

# Background characterisation — use scored overlaps
tolrm9_inferred  <- scored_consistent_overlaps(gr_BOTv_BOTR, gr_BOTCv_BOTR)
tolrm9_empirical <- gr_BOTCv_BOTC
tolrm9_validated <- scored_consistent_overlaps(tolrm9_inferred, tolrm9_empirical)
BOTv_BOTCv_div   <- find_opposite_overlaps(gr_BOTv_BOTR, gr_BOTCv_BOTR)
cat("Inferred tolrm9 background :", length(tolrm9_inferred), "\n")
cat("Empirical tolrm9 background:", length(tolrm9_empirical), "\n")
cat("Validated (scored)         :", length(tolrm9_validated),
    sprintf("(%.0f%% of inferred)\n\n",
            100 * length(tolrm9_validated) / max(length(tolrm9_inferred), 1)))

BOTv_above_inferred_bg  <- if (length(tolrm9_inferred) > 0)
  gr_BOTv_BOTR[countOverlaps(gr_BOTv_BOTR, tolrm9_inferred) == 0] else gr_BOTv_BOTR
BOTCv_above_inferred_bg <- if (length(tolrm9_inferred) > 0)
  gr_BOTCv_BOTR[countOverlaps(gr_BOTCv_BOTR, tolrm9_inferred) == 0] else gr_BOTCv_BOTR
BOTCv_above_empirical_bg <- if (length(tolrm9_empirical) > 0)
  gr_BOTCv_BOTR[countOverlaps(gr_BOTCv_BOTR, tolrm9_empirical) == 0] else gr_BOTCv_BOTR
BOTv_BOTCv_filtered     <- if (length(tolrm9_inferred) > 0)
  gr_BOTv_BOTCv[countOverlaps(gr_BOTv_BOTCv, tolrm9_inferred) == 0] else gr_BOTv_BOTCv

cat("BOTv above inferred bg :", length(BOTv_above_inferred_bg), "\n")
cat("BOTCv above inferred bg:", length(BOTCv_above_inferred_bg), "\n")
cat("BOTCv above empirical bg:", length(BOTCv_above_empirical_bg), "\n\n")

BOTv_only <- if (length(BOTv_above_inferred_bg) > 0 && length(gr_BOTv_BOTCv) > 0) {
  scored_consistent_overlaps(BOTv_above_inferred_bg, gr_BOTv_BOTCv)
} else GRanges()

BOTCv_only <- if (length(BOTCv_above_inferred_bg) > 0 && length(gr_BOTv_BOTCv) > 0) {
  c(subsetByOverlaps(BOTCv_above_inferred_bg[BOTCv_above_inferred_bg$log2FC>0],
                     gr_BOTv_BOTCv[gr_BOTv_BOTCv$log2FC<0], minoverlap=1L),
    subsetByOverlaps(BOTCv_above_inferred_bg[BOTCv_above_inferred_bg$log2FC<0],
                     gr_BOTv_BOTCv[gr_BOTv_BOTCv$log2FC>0], minoverlap=1L))
} else GRanges()

cat("BOTv_only cross-validated :", length(BOTv_only), "\n")
cat("BOTCv_only cross-validated:", length(BOTCv_only), "\n\n")

# Scored overlaps for all biological category assignments
Cic_cross_bg           <- scored_consistent_overlaps(gr_BOTC_BOTR, gr_BOTCv_BOTR)
Runt_het_vs_null       <- gr_BOT_hR_BOTR
Runt_WT_vs_null        <- gr_BOT_BOTR
Runt_WT_vs_het         <- gr_BOT_BOT_hR
Runt_GoF_vs_het        <- gr_BOTv_BOT_hR
Runt_het_vs_HLH        <- gr_BOT_hR_BOTCv
Runt_WT_gradient       <- scored_consistent_overlaps(gr_BOT_BOT_hR, gr_BOT_hR_BOTR)
Runt_full_gradient     <- scored_consistent_overlaps(gr_BOTv_BOT_hR, gr_BOT_hR_BOTR)
Cic_effect_Runt_intact <- gr_BOT_BOTC
BOTC_oR_rescue         <- gr_BOTCoR_BOTC
Runt_rescue_cross_bg   <- scored_consistent_overlaps(gr_BOTCoR_BOTC, gr_BOT_BOTR)

cat("Cic cross-background :", length(Cic_cross_bg), "\n")
cat("Runt WT gradient     :", length(Runt_WT_gradient), "\n")
cat("Runt full gradient   :", length(Runt_full_gradient), "\n")
cat("Cic w/ Runt intact   :", length(Cic_effect_Runt_intact), "\n")
cat("BOTC_oR rescue       :", length(BOTC_oR_rescue), "\n")
cat("Rescue cross-context :", length(Runt_rescue_cross_bg), "\n\n")

################################################################################
# SAVE ALL OUTPUTS
################################################################################

cat(strrep("=", 80), "\n"); cat("SAVING RESULTS\n"); cat(strrep("=", 80), "\n\n")

for (nm in names(dar_list)) write_dar_outputs(dar_list[[nm]], nm, output_dir)

write_dar_outputs(tolrm9_inferred,          "tolrm9_inferred_bg",              output_dir)
write_dar_outputs(tolrm9_empirical,         "tolrm9_empirical_bg",             output_dir)
write_dar_outputs(tolrm9_validated,         "tolrm9_validated_both_methods",   output_dir)
write_dar_outputs(BOTv_BOTCv_div,           "BOTv_BOTCv_divergent_vs_BOTR",    output_dir)
write_dar_outputs(BOTv_above_inferred_bg,   "BOTv_above_inferred_bg",          output_dir)
write_dar_outputs(BOTCv_above_inferred_bg,  "BOTCv_above_inferred_bg",         output_dir)
write_dar_outputs(BOTCv_above_empirical_bg, "BOTCv_above_empirical_bg",        output_dir)
write_dar_outputs(BOTv_BOTCv_filtered,      "BOTv_BOTCv_bg_subtracted",        output_dir)
write_dar_outputs(BOTv_only,                "BOTv_only_crossvalidated",         output_dir)
write_dar_outputs(BOTCv_only,               "BOTCv_only_crossvalidated",        output_dir)
write_dar_outputs(Cic_cross_bg,             "Cic_effect_cross_background",     output_dir)
write_dar_outputs(Cic_effect_Runt_intact,   "Cic_effect_Runt_intact_BOT_BOTC", output_dir)
write_dar_outputs(Runt_het_vs_null,         "Runt_dosage_het_vs_null",         output_dir)
write_dar_outputs(Runt_WT_vs_null,          "Runt_dosage_WT_vs_null",          output_dir)
write_dar_outputs(Runt_WT_vs_het,           "Runt_dosage_WT_vs_het",           output_dir)
write_dar_outputs(Runt_GoF_vs_het,          "Runt_dosage_GoF_vs_het",          output_dir)
write_dar_outputs(Runt_het_vs_HLH,          "Runt_het_vs_HLH54F_vent",         output_dir)
write_dar_outputs(Runt_WT_gradient,         "Runt_WT_dosage_gradient_non_vent",output_dir)
write_dar_outputs(Runt_full_gradient,       "Runt_full_gradient_GoF_to_null",  output_dir)
write_dar_outputs(BOTC_oR_rescue,           "BOTC_oR_Runt_rescue",             output_dir)
write_dar_outputs(Runt_rescue_cross_bg,     "Runt_rescue_cross_Cic_context",   output_dir)

################################################################################
# SUMMARY
################################################################################

cat(strrep("=", 80), "\n"); cat("SUMMARY\n"); cat(strrep("=", 80), "\n\n")

contrast_summary <- data.frame(
  Comparison   = names(dar_list),
  Total_DARs   = sapply(dar_list, length),
  Gained_open  = sapply(dar_list, function(gr) sum(mcols(gr)$direction == "gained-open")),
  Gained_close = sapply(dar_list, function(gr) sum(mcols(gr)$direction == "gained-close")),
  Tier1_n      = sapply(dar_list, function(gr) sum(mcols(gr)$tier == "condition_specific")),
  Tier2_n      = sapply(dar_list, function(gr)
                   sum(mcols(gr)$tier == "differential" & grepl("^NARROW", mcols(gr)$mode))),
  Tier3_n      = sapply(dar_list, function(gr)
                   sum(mcols(gr)$tier == "differential" & grepl("^EXTENDED", mcols(gr)$mode))),
  Consensus_4method = sapply(dar_list, function(gr)
                   sum(!is.na(mcols(gr)$m_agree) & mcols(gr)$m_agree >= 2.0)),
  Consensus_2method = sapply(dar_list, function(gr)
                   sum(!is.na(mcols(gr)$m_agree) & mcols(gr)$m_agree >= 1.5))
)
print(contrast_summary); cat("\n")

write.table(contrast_summary,
  file.path(output_dir, "Summary_per_contrast.txt"),
  sep="\t", quote=FALSE, row.names=FALSE)

# Tier 3 diagnostic summary across contrasts
diag_files <- list.files(diag_dir, pattern="_Tier3_diagnostic.txt", full.names=TRUE)
if (length(diag_files) > 0) {
  cat("TIER 3 DIAGNOSTIC SUMMARY:\n")
  diag_summary <- do.call(rbind, lapply(diag_files, function(f) {
    d <- read.table(f, header=TRUE, sep="\t")
    tbl <- table(d$classification)
    data.frame(contrast=gsub("_Tier3_diagnostic.txt", "", basename(f)),
               total_tier3=nrow(d),
               tight_boundary=sum(tbl["tight_boundary"], na.rm=TRUE),
               adjacent_merged=sum(tbl["adjacent_merged"], na.rm=TRUE),
               diffuse=sum(tbl["diffuse"], na.rm=TRUE))
  }))
  print(diag_summary)
  write.table(diag_summary,
    file.path(diag_dir, "Tier3_summary_all_contrasts.txt"),
    sep="\t", quote=FALSE, row.names=FALSE)
  cat("\n")
}

################################################################################
# PLOTS
################################################################################

plot_data <- contrast_summary %>%
  select(Comparison, Gained_open, Gained_close) %>%
  pivot_longer(cols=c(Gained_open, Gained_close),
               names_to="direction", values_to="count") %>%
  mutate(direction=factor(direction,
    levels=c("Gained_open","Gained_close"),
    labels=c("Gained accessibility","Lost accessibility")),
    Comparison=factor(Comparison, levels=unique(Comparison)))

p1 <- ggplot(plot_data, aes(x=Comparison, y=count, fill=direction)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values=c("Gained accessibility"="deepskyblue",
                              "Lost accessibility"  ="firebrick3")) +
  theme_minimal(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top") +
  labs(title="DAR Counts Per Contrast — v3 (L1-L5+L7)",
       x=NULL, y="Number of DARs", fill=NULL)
ggsave(file.path(output_dir, "Plot1_DARs_per_contrast.pdf"), p1, width=14, height=5)

tier_data <- contrast_summary %>%
  select(Comparison, Tier1_n, Tier2_n, Tier3_n) %>%
  pivot_longer(cols=c(Tier1_n, Tier2_n, Tier3_n),
               names_to="tier", values_to="n") %>%
  mutate(tier=factor(tier,
    levels=c("Tier1_n","Tier2_n","Tier3_n"),
    labels=c("Tier 1: cond-specific","Tier 2: narrow diff","Tier 3: ext diff")))

p2 <- ggplot(tier_data, aes(x=Comparison, y=n, fill=tier)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("Tier 1: cond-specific"="#d62728",
                              "Tier 2: narrow diff"  ="#1f77b4",
                              "Tier 3: ext diff"     ="#aec7e8")) +
  theme_minimal(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top") +
  labs(title="DAR Tier Breakdown Per Contrast",
       x=NULL, y="Number of DARs", fill=NULL)
ggsave(file.path(output_dir, "Plot2_tier_breakdown.pdf"), p2, width=14, height=5)

conf_data <- contrast_summary %>%
  select(Comparison, Consensus_4method, Consensus_2method, Tier1_n) %>%
  pivot_longer(cols=c(Consensus_4method, Consensus_2method, Tier1_n),
               names_to="category", values_to="n") %>%
  mutate(category=factor(category,
    levels=c("Tier1_n","Consensus_4method","Consensus_2method"),
    labels=c("Tier 1 (cond-specific)","4-method consensus","2-method consensus")))

p3 <- ggplot(conf_data, aes(x=Comparison, y=n, fill=category)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values=c("Tier 1 (cond-specific)"="#d62728",
                              "4-method consensus"    ="#2ca02c",
                              "2-method consensus"    ="#98df8a")) +
  theme_minimal(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top") +
  labs(title="High-Confidence DARs by Method Agreement (L5)",
       x=NULL, y="Number of DARs", fill=NULL)
ggsave(file.path(output_dir, "Plot3_confidence_breakdown.pdf"), p3, width=14, height=5)

cat("Plots saved.\n\n")
cat(strrep("=", 80), "\n"); cat("ANALYSIS COMPLETE — v3\n"); cat(strrep("=", 80), "\n\n")
cat("Output directory:", output_dir, "\n")
cat("Tier 1 outputs  :", tier1_dir, "\n")
cat("Diagnostics     :", diag_dir, "\n\n")
cat("New columns in annotated outputs:\n")
cat("  conf_score — continuous confidence (method_agreement × tier_score × -log10p)\n")
cat("  m_agree    — method agreement: 1.0=single, 1.5=2-method, 2.0=4-method\n")
cat("  cross_conf_score — (scored overlap outputs only) conf score of overlap partner\n\n")
