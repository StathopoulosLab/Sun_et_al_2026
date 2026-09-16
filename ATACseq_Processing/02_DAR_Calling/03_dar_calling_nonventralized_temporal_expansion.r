################################################################################
# STEP 1 — NON-VENT TEMPORAL EXPANSION  —  v4
#
# Changes in v3 (applied inside run_temporal_dar, Section 8) — preserved in v4:
#   [L1]-[L5]+[L7] as documented below.
#
# New in v4 (performance — requires Precount_AllCounts_Matrices.r run first):
#   - FC_THREADS read from SLURM_CPUS_PER_TASK (used in fallback/UMAP only)
#   - CACHE_DIR configured for NonVent_temporal cache block
#   - run_temporal_dar() loads narrow counts, grp_cache, and ext_cache from
#     disk instead of calling featureCounts; grp_ext_cache pre-loaded to
#     eliminate repeated Tier 1 extended counting per contrast
#   - All featureCounts fallback paths retained for first-run safety
#
# Previous enhancements (v2):
#   Tier 1 named outputs | baseMean/AveExpr | per-contrast extension_bp |
#   priority-scored deduplication
#
# New in v3 — applied inside run_temporal_dar() (Section 8):
#
#  [L1] GROUP-SPECIFIC PEAKS FOR TIER 1
#       GroupPeaks_*.bed files (from Build_FullUniverse_AllGeno.sh) are loaded
#       and used as Tier 1 counting substrates instead of the union peaks.
#
#  [L2] CPM-SCALED TIER 1 THRESHOLDS
#       min_reads / max_reads computed from CPM × median library size,
#       scaled to replicate count. Singleton groups (n=1) handled correctly.
#
#  [L3] CONTRAST-SPECIFIC RUV CONTROL SELECTION
#       Pairwise edgeR between the two groups in each contrast selects
#       control peaks stable for that specific comparison.
#
#  [L4] BOUNDARY-BASED PEAK EXTENSION
#       Extended peaks grow from start/end boundaries, not peak center.
#       Adjacent peaks within ext_bp merge before counting.
#
#  [L5] CONTINUOUS CONFIDENCE SCORE
#       conf_score = method_agreement × tier_score × −log10(p). Propagates
#       into all output files and write_dar_temporal().
#
#  [L7] TIER 3 DIAGNOSTIC
#       Extended-only DARs classified as tight_boundary / adjacent_merged /
#       diffuse per contrast. Written to Diagnostics/ subdirectory.
#
# Sections 1–6 (BAM registry, replicate decisions, UMAP QC, staging test,
# coldata builder, peak universe extension) are UNCHANGED from v2.
# Only Section 8 (run_temporal_dar) is upgraded.
#
# Adds later timepoints to the non-ventralized genotype series and resolves
# two ambiguous/problematic replicates before downstream DAR analysis.
#
# NEW SAMPLES IN THIS SCRIPT:
#
#   BOTR_gastr_rep1  — Run-D7 true homozygous null, gastrulation
#   BOTR_gastr_rep2  — Run-D7 het (BOT_hR background), gastrulation
#                      NOTE: rep2 is genotypically BOT_hR not BOTR — assigned
#                      to BOT_hR_gastr group, not BOTR_gastr.
#   BOT_nc14d_rep2   — D7 WT control, nc14d (singleton; rep1 removed: failed UMAP)
#   BOTC_nc14late_rep1/2 — HLH B6 nc14late (rep1 confirmed; rep2 = Mat_Run_B6_4 nc14b interim)
#   BOTC_oR_late_rep1/2 — mat-tub>Run in B6, late/gastr (confirmed pair)
#
# REPLICATE DECISIONS FORMALIZED HERE:
#
#   BOTC_oR nc14b:
#     rep1 (Mat_Run_B6_3) — CONFIRMED: clusters with BOT (expected — Runt rescue)
#     rep2 (Mat_Run_B6_4) — EXCLUDED: still clusters with BOTC (incomplete rescue)
#     Action: use rep1 as singleton; note in all downstream outputs
#
#   BOT_late:
#     rep1 — clusters with BOT (consistent with nc14b identity → possibly younger)
#     rep2 — clusters with BOTC (consistent with nc14late → Runt expression fading)
#     Action: treat as INDIVIDUALS (BOT_late_rep1 / BOT_late_rep2) in QC pass;
#             run staging test (Section 4) to evaluate whether rep1 is nc14b misassigned.
#
# OUTPUT TRACKS ADDED TO BigWig / BED PIPELINE:
#   BOTR_gastr, BOT_hR_gastr, BOT_nc14d, BOTC_oR_late (as singletons or pairs)
#
# This script:
#   Section 1  — BAM registry: all new + revised candidates
#   Section 2  — BOTC_oR replicate decision and singleton DAR fallback
#   Section 3  — Temporal series QC (UMAP placement + pairwise correlations)
#   Section 4  — BOT_late staging test (nuclear cycle marker genes via ATAC proxy)
#   Section 5  — Updated coldata for downstream DAR scripts
#   Section 6  — Peak universe extension (add new samples to union BED)
#   Section 7  — Recommended contrasts for each new timepoint group
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
  library(uwot)       # for UMAP
  library(ggrepel)
})

# Thread count — reads SLURM allocation automatically; falls back to 4.
# Used only in UMAP featureCounts call and fallback paths.
FC_THREADS <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))

# Cache directory written by Precount_AllCounts_Matrices.r
CACHE_DIR <- "./counts_cache/NonVent_temporal"

# ---------------------------------------------------------------------------
# RUN_QC_UMAP: set FALSE once all candidate replicates are confirmed.
# The UMAP/correlation section runs a full featureCounts call on all 23 BAMs
# and exists only to evaluate whether new candidate samples cluster correctly.
# With all replicates confirmed, set this to FALSE to skip it entirely and
# go straight to the DAR analysis.
# ---------------------------------------------------------------------------
RUN_QC_UMAP <- FALSE

################################################################################
# SECTION 1 — BAM REGISTRY
################################################################################

# ---------------------------------------------------------------------------
# DIRECTORY ROOTS
# ---------------------------------------------------------------------------
BAM_CTRL     <- "./Control_Bams"
BAM_CTRL_POT <- "./Control_Bams/Fully_flattened_potential"
BAM_MAIN     <- "../../ATAC_ChIP_Integration_08.04.25/INDV_Bams/Fully_flattened"

output_dir  <- "./Temporal_Expansion_NonVent"
figures_dir <- file.path(output_dir, "figures")
tier1_dir   <- file.path(output_dir, "Tier1_ConditionSpecific")
diag_dir    <- file.path(output_dir, "Diagnostics")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tier1_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(diag_dir,    recursive = TRUE, showWarnings = FALSE)

# Global default extension for broad peak analysis; override per-contrast below
extension_bp_global <- 200
# [L2] CPM-based Tier 1 thresholds (used in Section 8)
min_cpm_present    <- 1.0    # CPM floor for "present" condition
max_cpm_absent     <- 0.25   # CPM ceiling per replicate for "absent" condition
# Tier 1 thresholds (legacy — used in old identify_condition_specific_pair)
min_reads_per_peak <- 8
max_reads_other    <- 3

# Group-specific peak BEDs (produced by Build_FullUniverse_AllGeno.sh)
peak_dir         <- "./ATAC_NarrowPeaks"
group_peak_beds  <- list(
  BOT     = file.path(peak_dir, "GroupPeaks_BOT.bed"),
  BOT_hR  = file.path(peak_dir, "GroupPeaks_BOT_hR.bed"),
  BOTR    = file.path(peak_dir, "GroupPeaks_BOTR.bed"),
  BOTC    = file.path(peak_dir, "GroupPeaks_BOTC.bed"),
  BOTC_oR = file.path(peak_dir, "GroupPeaks_BOTC_oR.bed"),
  BOTv    = file.path(peak_dir, "GroupPeaks_BOTv.bed"),
  BOTCv   = file.path(peak_dir, "GroupPeaks_BOTCv.bed")
)

# ---------------------------------------------------------------------------
# CONFIRMED nc14b PAIRS (from existing Step1_SevenGeno_nc14b — reference anchors)
# These are included in ALL temporal QC plots so new samples can be placed
# relative to the known coordinate system.
# ---------------------------------------------------------------------------
confirmed_nc14b_bams <- list(
  BOTv_nc14b_rep1   = file.path(BAM_CTRL_POT, "FoxL1-High_Dm_ATAC_Nc14b_rep3_noq_rmdup.noChrM.bam"),
  BOTv_nc14b_rep2   = file.path(BAM_CTRL_POT, "FoxL1-High_Dm_ATAC_Nc14b_07_noq_rmdup.noChrM.bam"),
  BOTCv_nc14b_rep1  = file.path(BAM_CTRL_POT, "HLH54F-High_Dm_ATAC_Nc14b_rep4_noq_rmdup.noChrM.bam"),
  BOTCv_nc14b_rep2  = file.path(BAM_CTRL_POT, "HLH54F-High_Dm_ATAC_Nc14b_rep5_noq_rmdup.noChrM.bam"),
  BOT_nc14b_rep1    = file.path(BAM_CTRL, "BOT_D7_1_nc14b_IR_noq_rmdup.noChrM.bam"),
  BOT_nc14b_rep2    = file.path(BAM_CTRL, "BOT_D7_2_nc14b_IR_noq_rmdup.noChrM.bam"),
  BOT_hR_nc14b_rep1 = file.path(BAM_CTRL, "Run-D7_Dm_ATAC_Nc14b_01_noq_rmdup.noChrM.bam"),
  BOT_hR_nc14b_rep2 = file.path(BAM_CTRL, "Run-D7_Dm_ATAC_Nc14b_02_noq_rmdup.noChrM.bam"),
  BOTR_nc14b_rep1   = file.path(BAM_CTRL, "Run-D7_Dm_ATAC_Nc14b_03_noq_rmdup.noChrM.bam"),
  BOTR_nc14b_rep2   = file.path(BAM_CTRL, "Run-D7_Dm_ATAC_Nc14b_07_noq_rmdup.noChrM.bam"),
  BOTC_nc14b_rep1   = file.path(BAM_CTRL, "HLH_B6_Nc14b_01_noq_rmdup.noChrM.bam"),
  BOTC_nc14b_rep2   = file.path(BAM_CTRL, "HLH_B6_Nc14b_02_noq_rmdup.noChrM.bam"),
  # BOTC_oR nc14b: confirmed pair (rep1 + rep2; r=0.9405)
  # Mat_Run_B6_4 remains excluded — clusters with BOTC (incomplete rescue)
  BOTC_oR_nc14b_rep1 = file.path(BAM_CTRL, "Mat_Run_B6_3_nc14b_IR_noq_rmdup.noChrM.bam"),
  BOTC_oR_nc14b_rep2 = file.path(BAM_CTRL, "Mat_Run_B6_5_nc14b_IR_noq_rmdup.noChrM.bam")
)

# ---------------------------------------------------------------------------
# NEW CANDIDATE BAMS — UPDATE PATHS BELOW WITH YOUR ACTUAL FILE NAMES
# These paths follow the naming convention from generate_CPM_bigwigs.sh.
# Adjust filenames to match your actual BAM names.
# ---------------------------------------------------------------------------
new_candidate_bams <- list(

  # BOTR gastrulation — confirmed true null
  BOTR_gastr_rep1   = file.path(BAM_CTRL, "BOTR_rD7_1_nc14gastr_IR_noq_rmdup.noChrM.bam"),

  # BOT_hR gastrulation — confirmed het pair
  # rep1 = BOTR_rD7_2 (original confirmed het)
  # rep2 = RunD7_3_gastr_IR (28368, confirmed het by IGV BigWig inspection)
  BOT_hR_gastr_rep1 = file.path(BAM_CTRL, "BOTR_rD7_2_nc14gastr_IR_noq_rmdup.noChrM.bam"),
  BOT_hR_gastr_rep2 = file.path(BAM_CTRL, "RunD7_3_gastr_IR_noq_rmdup.noChrM.bam"),

  # BOT nc14d — UPGRADED to confirmed pair (28366 rep1 + existing rep2)
  BOT_nc14d_rep1    = file.path(BAM_CTRL, "D7_1_nc14d_IR_noq_rmdup.noChrM.bam"),
  BOT_nc14d_rep2    = file.path(BAM_CTRL, "BOT_D7_2_nc14d_IR_noq_rmdup.noChrM.bam"),

  # BOTC nc14late — confirmed rep1 + Mat_Run_B6_4 nc14b interim proxy
  # Mat_Run_B6_4 clusters BOTC-like by UMAP (excluded from BOTC_oR); acceptable
  # interim pair. Replace rep2 with confirmed nc14late BOTC rep when available.
  BOTC_nc14late_rep1 = file.path(BAM_CTRL, "HLH_B6_Nc14late_01_noq_rmdup.noChrM.bam"),
  BOTC_nc14late_rep2 = file.path(BAM_CTRL, "Mat_Run_B6_4_nc14b_IR_noq_rmdup.noChrM.bam"),

  # BOTC_oR late — confirmed pair (your description: 2 potential late reps confirmed)
  BOTC_oR_late_rep1 = file.path(BAM_CTRL, "Mat_Run_B6_1_nc14d_IR_noq_rmdup.noChrM.bam"),
  BOTC_oR_late_rep2 = file.path(BAM_CTRL, "Mat_Run_B6_2_nc14d_IR_noq_rmdup.noChrM.bam")

)

# Combined: all anchors + candidates for temporal QC
all_bams <- c(confirmed_nc14b_bams, new_candidate_bams)

# ---------------------------------------------------------------------------
# COLDATA — full sample metadata for QC UMAP
# Add rows here as you confirm each new sample.
# ---------------------------------------------------------------------------
temporal_coldata <- data.frame(
  sample   = names(all_bams),
  bam_path = unlist(all_bams),
  genotype = c(
    # confirmed nc14b
    "BOTv","BOTv","BOTCv","BOTCv",
    "BOT","BOT","BOT_hR","BOT_hR","BOTR","BOTR","BOTC","BOTC",
    "BOTC_oR","BOTC_oR",          # rep1 + rep2 (both confirmed)
    # new candidates
    "BOTR","BOT_hR","BOT_hR",         # gastr: null rep1, het rep1+rep2
    "BOT","BOT",                      # nc14d confirmed pair (28366 + rep2)
    "BOTC","BOTC",                    # nc14late: confirmed + interim proxy
    "BOTC_oR","BOTC_oR"               # late pair
  ),
  timepoint = c(
    # confirmed nc14b
    rep("nc14b", 14),
    # new candidates
    "gastr","gastr","gastr",          # BOTR_gastr + BOT_hR_gastr rep1 + rep2
    "nc14d","nc14d",                  # BOT nc14d confirmed pair
    "nc14late","nc14late",            # BOTC nc14late
    "late","late"
  ),
  status = c(
    rep("confirmed", 14),
    rep("candidate", 9)               # 3 gastr + 2 nc14d + 2 nc14late + 2 late
  ),
  stringsAsFactors = FALSE
)

# Derived group label for plotting
temporal_coldata$group <- paste(temporal_coldata$genotype, temporal_coldata$timepoint, sep = "_")

# ---------------------------------------------------------------------------
# PROJECT-WIDE COLOR PALETTE
# ---------------------------------------------------------------------------
genotype_colors <- c(
  "BOT"      = "#1f77b4",   # blue
  "BOT_hR"   = "#17becf",   # teal
  "BOTR"     = "#9467bd",   # purple
  "BOTC"     = "#ff7f0e",   # orange
  "BOTC_oR"  = "#bcbd22",   # yellow-green
  "BOTv"     = "#2ca02c",   # green
  "BOTCv"    = "#e377c2"    # pink
)

timepoint_shapes <- c(
  "nc14b"        = 16,   # filled circle
  "nc14d"        = 17,   # filled triangle up
  "nc14late"     = 15,   # filled square
  "late"         = 18,   # filled diamond
  "gastr"        = 8,    # asterisk / starburst
  "unknown_late" = 1     # open circle — unconfirmed
)

################################################################################
# SECTION 2 — BOTC_oR REPLICATE DECISION
################################################################################

# DECISION: BOTC_oR nc14b is now a confirmed PAIR (n=2)
# ---------------------------------------------------------------------------
# BOTC_oR nc14b:
#   rep1 (Mat_Run_B6_3): confirmed — r vs BOT_nc14b high; clusters correctly
#   rep2 (Mat_Run_B6_5): confirmed — r vs rep1 = 0.9405; QC passes
#   rep2 (Mat_Run_B6_4): EXCLUDED — clusters with BOTC (incomplete rescue)
#
# Consequence for downstream DAR analysis:
#   - BOTC_oR_vs_BOTC and BOTC_oR_vs_BOTR contrasts now use n=2 for BOTC_oR
#   - DESeq2 + limma-voom BOTH active (no longer limma-only)
#   - Method agreement scoring (L5) can now reach 1.5 or 2.0 for BOTC_oR DARs
#   - Confidence scores for BOTC_oR DARs meaningfully upgraded vs v1 results
# ---------------------------------------------------------------------------

botcor_nc14b_decision <- list(
  confirmed_rep1   = "Mat_Run_B6_3_nc14b_IR_noq_rmdup.noChrM.bam",
  confirmed_rep2   = "Mat_Run_B6_5_nc14b_IR_noq_rmdup.noChrM.bam",
  excluded_rep     = "Mat_Run_B6_4_nc14b_IR_noq_rmdup.noChrM.bam",
  exclusion_reason = "Clusters with BOTC in UMAP; Runt rescue signature absent",
  analysis_mode    = "DESeq2_plus_limma",
  use_limma_only   = FALSE,
  note             = "n=2 — both DESeq2 and limma active for BOTC_oR contrasts"
)

cat("BOTC_oR nc14b decision:\n")
cat("  Confirmed rep1:  ", botcor_nc14b_decision$confirmed_rep1, "\n")
cat("  Confirmed rep2:  ", botcor_nc14b_decision$confirmed_rep2, "\n")
cat("  Excluded:        ", botcor_nc14b_decision$excluded_rep, "\n")
cat("  Reason:          ", botcor_nc14b_decision$exclusion_reason, "\n")
cat("  Analysis mode:   ", botcor_nc14b_decision$analysis_mode, "\n\n")

################################################################################
# SECTION 3 — TEMPORAL QC: UMAP + PAIRWISE CORRELATIONS
################################################################################

# Peak universe BED — use the full 7-genotype union from Step1_SevenGeno
existing_peaks_bed <- "./ATAC_NarrowPeaks/FullUniverse_nc14b_AllGeno_union_peaks.bed"

import_bed <- function(f) {
  d <- read.table(f, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                  comment.char = "")
  GRanges(
    seqnames = d[,1],
    ranges   = IRanges(start = d[,2] + 1, end = d[,3]),
    strand   = "*"
  )
}

# ---------------------------------------------------------------------------
# identify_condition_specific_pair — with baseMean/AveExpr proxies
# ---------------------------------------------------------------------------
identify_condition_specific_pair <- function(counts, peaks_gr,
    grpA_idx, grpB_idx, grpA_label, grpB_label,
    min_reads, max_reads_other, mode_label, contrast_name) {
  A_tot <- rowSums(counts[,grpA_idx,drop=FALSE])
  B_tot <- rowSums(counts[,grpB_idx,drop=FALSE])
  n_A   <- length(grpA_idx); n_B <- length(grpB_idx)
  lib_size_est <- median(colSums(counts)) + 1
  mk <- function(idx, lfc, label) {
    if (sum(idx)==0) return(data.frame(
      peak_id=character(),chr=character(),start=integer(),end=integer(),
      log2FC=numeric(),pvalue=numeric(),padj=numeric(),
      AveExpr=numeric(),baseMean=numeric(),
      method=character(),contrast=character(),mode=character(),
      stringsAsFactors=FALSE))
    present_mean <- ifelse(lfc<0, A_tot[idx]/max(n_A,1), B_tot[idx]/max(n_B,1))
    ave_expr_est <- log2((present_mean/lib_size_est*1e6) + 0.5)
    data.frame(
      peak_id  = names(peaks_gr)[idx],
      chr      = as.character(seqnames(peaks_gr[idx])),
      start    = start(peaks_gr[idx])-1,
      end      = end(peaks_gr[idx]),
      log2FC   = lfc,
      pvalue   = 1e-10,
      padj     = 1e-10,
      AveExpr  = ave_expr_est,
      baseMean = present_mean,
      method   = paste0(mode_label,"_",label,"_specific_",contrast_name),
      contrast = contrast_name,
      mode     = mode_label,
      stringsAsFactors=FALSE)
  }
  list(
    A_specific = mk(A_tot>=min_reads & B_tot<=max_reads_other, -5,
                    paste0(grpA_label,"_open_lost")),
    B_specific = mk(B_tot>=min_reads & A_tot<=max_reads_other,  5,
                    paste0(grpB_label,"_open_gained"))
  )
}

# ---------------------------------------------------------------------------
# df_to_gr — with tier_score for deduplication
# tier mapping: condition_specific=3, NARROW differential=2, EXTENDED=1
# ---------------------------------------------------------------------------
df_to_gr_temporal <- function(df, tier = "differential") {
  if (nrow(df)==0) return(GRanges())
  tier_score_val <- switch(tier,
    "condition_specific" = 3,
    "differential"       = ifelse(grepl("^NARROW",df$mode[1]), 2, 1),
    1)
  GRanges(seqnames = df$chr,
          ranges   = IRanges(start=df$start+1, end=df$end),
          log2FC   = df$log2FC,
          pvalue   = df$pvalue,
          padj     = if ("padj"    %in% names(df)) df$padj    else NA_real_,
          AveExpr  = if ("AveExpr" %in% names(df)) df$AveExpr else NA_real_,
          baseMean = if ("baseMean"%in% names(df)) df$baseMean else NA_real_,
          method   = if ("method"  %in% names(df)) df$method   else tier,
          contrast = if ("contrast"%in% names(df)) df$contrast else NA_character_,
          tier     = tier,
          tier_score = tier_score_val,
          mode     = if ("mode"    %in% names(df)) df$mode    else tier)
}

# ---------------------------------------------------------------------------
# priority_deduplicate — direction-aware, score-based cluster deduplication
# ---------------------------------------------------------------------------
priority_deduplicate <- function(gr) {
  if (length(gr) <= 1) return(gr)
  eps <- 1e-300
  mcols(gr)$dedup_score <- mcols(gr)$tier_score - log10(mcols(gr)$pvalue + eps)
  mcols(gr)$sign_dir    <- sign(mcols(gr)$log2FC)
  dedup_one <- function(g) {
    if (length(g) <= 1) return(g)
    ovl    <- findOverlaps(g, g, minoverlap=1L)
    ovl_df <- as.data.frame(ovl)
    ovl_df <- ovl_df[ovl_df$queryHits != ovl_df$subjectHits, ]
    if (nrow(ovl_df)==0) return(g)
    parent <- seq_len(length(g))
    find_root <- function(x) {
      while (parent[x]!=x) { parent[x] <<- parent[parent[x]]; x <- parent[x] }
      x
    }
    for (i in seq_len(nrow(ovl_df))) {
      a <- find_root(ovl_df$queryHits[i]); b <- find_root(ovl_df$subjectHits[i])
      if (a!=b) parent[a] <- b
    }
    clusters <- sapply(seq_len(length(g)), find_root)
    scores   <- mcols(g)$dedup_score
    keep <- sapply(unique(clusters), function(cl) {
      m <- which(clusters==cl); m[which.max(scores[m])]
    })
    g[sort(keep)]
  }
  result <- c(dedup_one(gr[mcols(gr)$sign_dir > 0]),
              dedup_one(gr[mcols(gr)$sign_dir < 0]))
  result
}

# ---------------------------------------------------------------------------
# write_tier1_outputs — saves named BED + annotated TXT per contrast per direction
# ---------------------------------------------------------------------------
write_tier1_outputs <- function(cs_narrow, cs_ext, contrast_name, t1dir) {
  dir.create(t1dir, recursive=TRUE, showWarnings=FALSE)
  write_one <- function(df, tag) {
    if (nrow(df)==0) return(invisible(NULL))
    df$direction <- ifelse(df$log2FC>0,"gained-open","gained-close")
    df$name <- paste0(tag,"_",seq_len(nrow(df))); df$strand <- "."
    write.table(df[,c("chr","start","end","name","pvalue","strand")],
      file.path(t1dir,paste0(tag,"_DARs.bed")),
      sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
    write.table(df,
      file.path(t1dir,paste0(tag,"_DARs_annotated.txt")),
      sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
    cat("  T1 saved",nrow(df),"→",paste0(tag,"_DARs.bed\n"))
  }
  A_lost <- rbind(
    cs_narrow$A_specific[cs_narrow$A_specific$log2FC<0,,drop=FALSE],
    cs_ext$A_specific   [cs_ext$A_specific$log2FC<0,,drop=FALSE])
  A_lost <- A_lost[!duplicated(paste(A_lost$chr,A_lost$start,A_lost$end)),]
  B_gained <- rbind(
    cs_narrow$B_specific[cs_narrow$B_specific$log2FC>0,,drop=FALSE],
    cs_ext$B_specific   [cs_ext$B_specific$log2FC>0,,drop=FALSE])
  B_gained <- B_gained[!duplicated(paste(B_gained$chr,B_gained$start,B_gained$end)),]
  write_one(A_lost,   paste0("TIER1_",contrast_name,"_A_open_lost"))
  write_one(B_gained, paste0("TIER1_",contrast_name,"_B_open_gained"))
  cat(sprintf("  TIER1 %s: A_lost=%d  B_gained=%d\n",
              contrast_name, nrow(A_lost), nrow(B_gained)))
}

# ---------------------------------------------------------------------------
# write_dar_temporal — includes baseMean, AveExpr, padj, tier_score
# ---------------------------------------------------------------------------
write_dar_temporal <- function(gr, tag, out_dir) {
  if (length(gr)==0) { cat("  Skip",tag,"(0)\n"); return(invisible(NULL)) }
  mc <- mcols(gr)
  df <- data.frame(
    chr       = as.character(seqnames(gr)),
    start     = start(gr)-1,
    end       = end(gr),
    name      = paste0(tag,"_",seq_along(gr)),
    score     = round(-log10(mc$pvalue+1e-300)*10),
    strand    = ".",
    log2FC    = mc$log2FC,
    pvalue    = mc$pvalue,
    padj      = if ("padj"      %in% names(mc)) mc$padj      else NA_real_,
    AveExpr   = if ("AveExpr"   %in% names(mc)) mc$AveExpr   else NA_real_,
    baseMean  = if ("baseMean"  %in% names(mc)) mc$baseMean  else NA_real_,
    direction = ifelse(mc$log2FC>0,"gained-open","gained-close"),
    tier      = if ("tier"      %in% names(mc)) mc$tier      else NA_character_,
    tier_score= if ("tier_score"%in% names(mc)) mc$tier_score else NA_real_,
    mode      = if ("mode"      %in% names(mc)) mc$mode      else NA_character_,
    method    = if ("method"    %in% names(mc)) mc$method    else NA_character_,
    contrast  = if ("contrast"  %in% names(mc)) mc$contrast  else NA_character_,
    stringsAsFactors=FALSE)
  write.table(df[,1:6],
    file.path(out_dir,paste0(tag,"_DARs.bed")),
    sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
  write.table(df,
    file.path(out_dir,paste0(tag,"_DARs_annotated.txt")),
    sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
  cat("  Saved",nrow(df),"→",paste0(tag,"_DARs.bed\n"))
}

run_temporal_umap <- function(bam_list, coldata_sub, peaks_bed_path, out_dir) {

  peaks_gr <- import_bed(peaks_bed_path)

  # ── UMAP counts cache ───────────────────────────────────────────────────────
  # featureCounts is expensive — cache the raw counts so re-sourcing the script
  # skips counting when the BAM list hasn't changed.
  umap_cache_file <- file.path(out_dir, "umap_counts_cache.txt")

  if (file.exists(umap_cache_file)) {
    cat("Loading UMAP counts from cache...\n")
    counts_cached <- as.matrix(read.table(umap_cache_file,
                                           header=TRUE, row.names=1, sep="\t"))
    if (all(coldata_sub$sample %in% colnames(counts_cached))) {
      counts <- counts_cached[, coldata_sub$sample, drop=FALSE]
      cat("  Loaded", ncol(counts), "samples from cache.\n")
    } else {
      missing_s <- coldata_sub$sample[!coldata_sub$sample %in% colnames(counts_cached)]
      cat("  Cache missing:", paste(missing_s, collapse=", "), "— recounting.\n")
      counts <- NULL
    }
  } else {
    counts <- NULL
  }

  if (is.null(counts)) {
    cat("Counting reads in peaks for", nrow(coldata_sub), "samples...\n")
    bed_df <- read.table(peaks_bed_path, header = FALSE, sep = "\t",
                         stringsAsFactors = FALSE)
    saf_df <- data.frame(
      GeneID = if (ncol(bed_df) >= 4) bed_df[,4] else paste0("peak_", seq_len(nrow(bed_df))),
      Chr    = bed_df[,1],
      Start  = bed_df[,2],
      End    = bed_df[,3],
      Strand = if (ncol(bed_df) >= 6) bed_df[,6] else "*",
      stringsAsFactors = FALSE
    )
    counts <- featureCounts(
      files           = bam_list,
      annot.ext       = saf_df,
      isGTFAnnotationFile = FALSE,
      useMetaFeatures = FALSE,
      isPairedEnd     = TRUE,
      nthreads        = FC_THREADS,
      allowMultiOverlap = TRUE,
      countMultiMappingReads = FALSE
    )$counts
    colnames(counts) <- coldata_sub$sample
    dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
    write.table(counts, umap_cache_file,
                sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
    cat("  UMAP counts cached →", umap_cache_file, "\n")
  }

  # CPM normalization + log2
  lib_sizes <- colSums(counts)
  cpm_mat   <- sweep(counts, 2, lib_sizes / 1e6, "/")
  log_cpm   <- log2(cpm_mat + 0.5)

  # Filter low-count peaks
  keep      <- rowSums(cpm_mat >= 0.5) >= 2
  log_cpm   <- log_cpm[keep, ]
  cat("Peaks after CPM filter:", sum(keep), "\n")

  # ── Pairwise Pearson correlations ─────────────────────────────────────────
  cor_mat <- cor(log_cpm, method = "pearson")
  cor_out <- file.path(out_dir, "temporal_sample_correlations.txt")
  write.table(round(cor_mat, 4), cor_out, sep = "\t", quote = FALSE)
  cat("Saved correlations →", cor_out, "\n")

  # ── UMAP ──────────────────────────────────────────────────────────────────
  set.seed(42)
  umap_res <- umap(
    t(log_cpm),
    n_neighbors  = min(15, floor(ncol(log_cpm) / 2)),
    min_dist     = 0.1,
    metric       = "cosine"
  )

  umap_df <- data.frame(
    UMAP1     = umap_res[, 1],
    UMAP2     = umap_res[, 2],
    sample    = coldata_sub$sample,
    genotype  = coldata_sub$genotype,
    timepoint = coldata_sub$timepoint,
    status    = coldata_sub$status
  )

  # ── Plot ──────────────────────────────────────────────────────────────────
  p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2,
                                 color = genotype,
                                 shape = timepoint,
                                 alpha = status)) +
    geom_point(size = 4, stroke = 1.2) +
    geom_text_repel(
      data    = subset(umap_df, status == "candidate"),
      mapping = aes(label = sample),
      size    = 3, max.overlaps = 20
    ) +
    scale_color_manual(values = genotype_colors) +
    scale_shape_manual(values = timepoint_shapes) +
    scale_alpha_manual(values = c("confirmed" = 0.75, "candidate" = 1.0),
                       guide = "none") +
    # Ellipses for confirmed nc14b groups
    stat_ellipse(
      data  = subset(umap_df, status == "confirmed"),
      aes(group = genotype),
      level = 0.80, linetype = "dashed", linewidth = 0.4, alpha = 0.4
    ) +
    labs(
      title    = "Temporal expansion QC — UMAP all non-vent samples",
      subtitle = "Candidates labeled; dashed ellipses = confirmed nc14b groups (80% CI)",
      x        = "UMAP1", y = "UMAP2",
      color    = "Genotype", shape = "Timepoint"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")

  ggsave(file.path(out_dir, "temporal_UMAP_allSamples.pdf"),
         p_umap, width = 10, height = 7)
  cat("Saved UMAP →", file.path(out_dir, "temporal_UMAP_allSamples.pdf"), "\n")

  # ── Nearest-neighbor table for candidates ─────────────────────────────────
  cand_samples <- coldata_sub$sample[coldata_sub$status == "candidate"]
  anch_samples <- coldata_sub$sample[coldata_sub$status == "confirmed"]

  nn_rows <- lapply(cand_samples, function(cand) {
    cors_to_anchors <- cor_mat[cand, anch_samples]
    top3            <- sort(cors_to_anchors, decreasing = TRUE)[1:min(3, length(cors_to_anchors))]
    data.frame(
      candidate     = cand,
      timepoint     = coldata_sub$timepoint[coldata_sub$sample == cand],
      nearest_1     = names(top3)[1], r1 = round(top3[1], 4),
      nearest_2     = names(top3)[2], r2 = round(top3[2], 4),
      nearest_3     = names(top3)[3], r3 = round(top3[3], 4),
      stringsAsFactors = FALSE
    )
  })

  nn_df <- do.call(rbind, nn_rows)
  nn_out <- file.path(out_dir, "temporal_candidate_nearest_neighbors.txt")
  write.table(nn_df, nn_out, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Nearest-neighbor table →", nn_out, "\n")
  print(nn_df)

  list(umap = umap_df, cor_mat = cor_mat, nn_table = nn_df, log_cpm = log_cpm)
}

# Build available_coldata (filter to BAMs that exist on disk)
existing_mask     <- file.exists(temporal_coldata$bam_path)
if (any(!existing_mask)) {
  cat("WARNING: Missing BAMs (will be skipped in QC):\n")
  cat(paste(" ", temporal_coldata$bam_path[!existing_mask]), sep = "\n")
}
available_coldata <- temporal_coldata[existing_mask, ]
available_bams    <- available_coldata$bam_path

# Run QC only if BAMs are available AND RUN_QC_UMAP is TRUE
temporal_qc <- NULL
if (RUN_QC_UMAP && nrow(available_coldata) >= 4) {
  temporal_qc <- run_temporal_umap(
    bam_list    = available_bams,
    coldata_sub = available_coldata,
    peaks_bed_path = existing_peaks_bed,
    out_dir     = figures_dir
  )
} else if (!RUN_QC_UMAP) {
  cat("Skipping temporal UMAP (RUN_QC_UMAP = FALSE — all replicates confirmed).\n")
} else {
  cat("Skipping temporal UMAP — fewer than 4 BAMs available.\n",
      "Update new_candidate_bams paths in Section 1 and re-run.\n")
}

################################################################################
# SECTION 4 — BOT_LATE STAGING TEST
################################################################################
#
# Problem: BOT_late_indvA clusters with BOT (nc14b-like); BOT_late_indvB clusters
# with BOTC (nc14late-like). Two hypotheses:
#
#   H1: indvA is actually nc14b (younger embryo, misassigned at collection).
#       Prediction: high r with BOT_nc14b pair, UMAP overlaps confirmed BOT nc14b.
#
#   H2: Both are nc14late but indvA has unusually high Runt signal (e.g. from a
#       slightly earlier collection within nc14late window).
#       Prediction: r with BOT_nc14b is lower than confirmed BOT pairs; some nc14late
#       marker accessibility differs between indvA and indvB.
#
# Test strategy:
#   (a) Pearson r vs confirmed BOT_nc14b pair — threshold: r >= 0.95 → H1 likely
#   (b) PC1 projection: project indvA/indvB onto PC1 axis defined by the
#       confirmed nc14b samples; if indvA PC1 score falls within the confirmed
#       nc14b distribution → H1.
#   (c) ATAC signal at Runt-target loci: if indvA has nc14b-level Runt-target
#       accessibility (using BOTR-open DARs as proxy for Runt-maintained sites),
#       that further supports H1.
#   (d) If H1 supported: reclassify indvA as BOT_nc14b and fold into nc14b analysis.
#       If neither hypothesis clearly wins: keep as individuals, flag as
#       "staging_ambiguous" and exclude from temporal regression but include
#       in exploratory figures.
# ---------------------------------------------------------------------------

run_staging_test <- function(qc_result, out_dir) {

  if (is.null(qc_result)) {
    cat("Staging test skipped — temporal QC not yet run.\n")
    return(invisible(NULL))
  }

  log_cpm <- qc_result$log_cpm
  cor_mat <- qc_result$cor_mat

  # (a) Pearson r vs confirmed BOT nc14b pair
  bot_nc14b_anchors <- c("BOT_nc14b_rep1", "BOT_nc14b_rep2")
  bot_nc14b_anchors <- intersect(bot_nc14b_anchors, colnames(cor_mat))
  candidates        <- c("BOT_late_indvA", "BOT_late_indvB")
  candidates        <- intersect(candidates, colnames(cor_mat))

  if (length(candidates) == 0) {
    cat("BOT_late candidates not yet available; staging test deferred.\n")
    return(invisible(NULL))
  }

  cat("\n=== BOT_LATE STAGING TEST ===\n")
  cat("Correlation with confirmed BOT nc14b anchors:\n")
  for (cand in candidates) {
    for (anch in bot_nc14b_anchors) {
      cat(sprintf("  %-20s vs %-20s  r = %.4f\n",
                  cand, anch, cor_mat[cand, anch]))
    }
  }

  # ── (b) PCA projection ───────────────────────────────────────────────────
  # Use only confirmed nc14b samples + the two BOT_late individuals
  staging_samples <- c(bot_nc14b_anchors,
                       intersect(c("BOTR_nc14b_rep1","BOTR_nc14b_rep2",
                                   "BOTC_nc14b_rep1","BOTC_nc14b_rep2"),
                                 colnames(log_cpm)),
                       candidates)

  log_sub  <- log_cpm[, staging_samples, drop = FALSE]
  pca_res  <- prcomp(t(log_sub), center = TRUE, scale. = FALSE)
  pca_df   <- data.frame(
    sample    = rownames(pca_res$x),
    PC1       = pca_res$x[, 1],
    PC2       = pca_res$x[, 2],
    is_candidate = rownames(pca_res$x) %in% candidates
  )

  # Confidence interval for confirmed BOT nc14b
  bot_nc14b_pc1 <- pca_df$PC1[pca_df$sample %in% bot_nc14b_anchors]
  bot_mean_pc1  <- mean(bot_nc14b_pc1)
  bot_sd_pc1    <- sd(bot_nc14b_pc1)

  cat("\nPC1 scores (staging axis):\n")
  for (s in staging_samples) {
    row   <- pca_df[pca_df$sample == s, ]
    flag  <- ""
    if (s %in% candidates) {
      z <- (row$PC1 - bot_mean_pc1) / (bot_sd_pc1 + 1e-6)
      flag <- sprintf("  [z vs BOT_nc14b = %.2f]", z)
    }
    cat(sprintf("  %-25s  PC1 = %7.3f%s\n", s, row$PC1, flag))
  }

  # Plot PCA staging
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2,
                               color = is_candidate,
                               label = sample)) +
    geom_point(size = 4) +
    geom_text_repel(size = 3) +
    # Confidence band for BOT nc14b
    annotate("rect",
             xmin = bot_mean_pc1 - 2*bot_sd_pc1,
             xmax = bot_mean_pc1 + 2*bot_sd_pc1,
             ymin = -Inf, ymax = Inf,
             alpha = 0.1, fill = "#1f77b4") +
    scale_color_manual(values = c("FALSE" = "#888888", "TRUE" = "#d62728"),
                       labels = c("FALSE" = "Confirmed nc14b", "TRUE" = "BOT_late candidate"),
                       name   = "") +
    labs(
      title    = "BOT_late staging test — PCA projection",
      subtitle = sprintf("Blue band = BOT nc14b 95%% CI (mean ± 2SD, PC1)\n  indvA within band → likely nc14b | indvA outside → likely nc14late"),
      x = sprintf("PC1 (%.1f%% var)", summary(pca_res)$importance[2,1]*100),
      y = sprintf("PC2 (%.1f%% var)", summary(pca_res)$importance[2,2]*100)
    ) +
    theme_bw(base_size = 12)

  ggsave(file.path(out_dir, "BOT_late_staging_test_PCA.pdf"),
         p_pca, width = 8, height = 6)
  cat("\nStaging test PCA →", file.path(out_dir, "BOT_late_staging_test_PCA.pdf"), "\n")

  # Return a staging call per candidate
  staging_calls <- lapply(candidates, function(cand) {
    pc1_val  <- pca_df$PC1[pca_df$sample == cand]
    z_score  <- (pc1_val - bot_mean_pc1) / (bot_sd_pc1 + 1e-6)
    r_mean   <- mean(cor_mat[cand, bot_nc14b_anchors])

    if (abs(z_score) < 2.0 && r_mean >= 0.94) {
      call <- "H1_likely_nc14b"
    } else if (abs(z_score) >= 2.0 && r_mean < 0.94) {
      call <- "H2_nc14late"
    } else {
      call <- "ambiguous"
    }
    data.frame(sample = cand, pc1 = pc1_val, z_vs_nc14b = round(z_score, 3),
               r_vs_nc14b_anchors = round(r_mean, 4), staging_call = call,
               stringsAsFactors = FALSE)
  })

  staging_df <- do.call(rbind, staging_calls)
  staging_out <- file.path(out_dir, "BOT_late_staging_calls.txt")
  write.table(staging_df, staging_out, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\nStaging calls:\n")
  print(staging_df)
  cat("\n  If staging_call == 'H1_likely_nc14b': reclassify indvA as BOT_nc14b_rep3\n",
      "  and include in the main Step1_SevenGeno_nc14b_DAR_Analysis.r BAM list.\n",
      "  If 'H2_nc14late' or 'ambiguous': treat as individual BOT_late exploratory\n",
      "  replicates; exclude from temporal regression but include in visual QC.\n\n")

  staging_df
}

staging_result <- run_staging_test(temporal_qc, figures_dir)

################################################################################
# SECTION 5 — UPDATED COLDATA FOR DOWNSTREAM DAR SCRIPTS
################################################################################
#
# After running Section 3+4:
#   1. Update 'staging_reclassify_indvA' below based on staging_result
#   2. This block generates the final temporal coldata that feeds
#      Step1_NonVent_Temporal_DAR_Analysis.r (written in Section 7)
# ---------------------------------------------------------------------------

# Set this TRUE if staging test supports H1 (indvA is nc14b)
staging_reclassify_indvA <- FALSE  # <-- update after reviewing staging_result

build_final_temporal_coldata <- function(reclassify_bot_late_A = FALSE) {

  # Confirmed nc14b (always included)
  rows <- list(
    # BOTv (ventralized anchors — included for cross-context comparison)
    data.frame(sample="BOTv_nc14b_rep1",  genotype="BOTv",   timepoint="nc14b",  bam=confirmed_nc14b_bams[["BOTv_nc14b_rep1"]],   status="confirmed"),
    data.frame(sample="BOTv_nc14b_rep2",  genotype="BOTv",   timepoint="nc14b",  bam=confirmed_nc14b_bams[["BOTv_nc14b_rep2"]],   status="confirmed"),
    # Non-vent nc14b
    data.frame(sample="BOT_nc14b_rep1",   genotype="BOT",    timepoint="nc14b",  bam=confirmed_nc14b_bams[["BOT_nc14b_rep1"]],    status="confirmed"),
    data.frame(sample="BOT_nc14b_rep2",   genotype="BOT",    timepoint="nc14b",  bam=confirmed_nc14b_bams[["BOT_nc14b_rep2"]],    status="confirmed"),
    data.frame(sample="BOT_hR_nc14b_rep1",genotype="BOT_hR", timepoint="nc14b",  bam=confirmed_nc14b_bams[["BOT_hR_nc14b_rep1"]], status="confirmed"),
    data.frame(sample="BOT_hR_nc14b_rep2",genotype="BOT_hR", timepoint="nc14b",  bam=confirmed_nc14b_bams[["BOT_hR_nc14b_rep2"]], status="confirmed"),
    data.frame(sample="BOTR_nc14b_rep1",  genotype="BOTR",   timepoint="nc14b",  bam=confirmed_nc14b_bams[["BOTR_nc14b_rep1"]],   status="confirmed"),
    data.frame(sample="BOTR_nc14b_rep2",  genotype="BOTR",   timepoint="nc14b",  bam=confirmed_nc14b_bams[["BOTR_nc14b_rep2"]],   status="confirmed"),
    data.frame(sample="BOTC_nc14b_rep1",  genotype="BOTC",   timepoint="nc14b",  bam=confirmed_nc14b_bams[["BOTC_nc14b_rep1"]],   status="confirmed"),
    data.frame(sample="BOTC_nc14b_rep2",  genotype="BOTC",   timepoint="nc14b",  bam=confirmed_nc14b_bams[["BOTC_nc14b_rep2"]],   status="confirmed"),
    # BOTC_oR — confirmed pair (n=2; Mat_Run_B6_3 + Mat_Run_B6_5)
    data.frame(sample="BOTC_oR_nc14b_rep1",genotype="BOTC_oR",timepoint="nc14b",bam=confirmed_nc14b_bams[["BOTC_oR_nc14b_rep1"]],status="confirmed"),
    data.frame(sample="BOTC_oR_nc14b_rep2",genotype="BOTC_oR",timepoint="nc14b",bam=confirmed_nc14b_bams[["BOTC_oR_nc14b_rep2"]],status="confirmed"),
    # New timepoints
    data.frame(sample="BOTR_gastr_rep1",  genotype="BOTR",   timepoint="gastr",  bam=new_candidate_bams[["BOTR_gastr_rep1"]],     status="candidate"),
    # BOT_hR_gastr UPGRADED to confirmed het pair
    data.frame(sample="BOT_hR_gastr_rep1",genotype="BOT_hR", timepoint="gastr",  bam=new_candidate_bams[["BOT_hR_gastr_rep1"]],  status="candidate"),
    data.frame(sample="BOT_hR_gastr_rep2",genotype="BOT_hR", timepoint="gastr",  bam=new_candidate_bams[["BOT_hR_gastr_rep2"]],  status="candidate"),
    # BOT nc14d UPGRADED to confirmed pair
    data.frame(sample="BOT_nc14d_rep1",   genotype="BOT",    timepoint="nc14d",  bam=new_candidate_bams[["BOT_nc14d_rep1"]],      status="candidate"),
    data.frame(sample="BOT_nc14d_rep2",   genotype="BOT",    timepoint="nc14d",  bam=new_candidate_bams[["BOT_nc14d_rep2"]],      status="candidate"),
    # BOTC nc14late — confirmed rep1 + Mat_Run_B6_4 nc14b interim proxy
    data.frame(sample="BOTC_nc14late_rep1",genotype="BOTC",  timepoint="nc14late",bam=new_candidate_bams[["BOTC_nc14late_rep1"]], status="candidate"),
    data.frame(sample="BOTC_nc14late_rep2",genotype="BOTC",  timepoint="nc14late",bam=new_candidate_bams[["BOTC_nc14late_rep2"]], status="interim_proxy"),
    data.frame(sample="BOTC_oR_late_rep1",genotype="BOTC_oR",timepoint="late",   bam=new_candidate_bams[["BOTC_oR_late_rep1"]],  status="candidate"),
    data.frame(sample="BOTC_oR_late_rep2",genotype="BOTC_oR",timepoint="late",   bam=new_candidate_bams[["BOTC_oR_late_rep2"]],  status="candidate")
  )

  # BOT_late: reclassify indvA if staging test supports H1
  if (reclassify_bot_late_A) {
    rows <- c(rows, list(
      data.frame(sample="BOT_nc14b_rep3",   genotype="BOT",    timepoint="nc14b",
                 bam=new_candidate_bams[["BOT_late_indvA"]],
                 status="reclassified_from_late", stringsAsFactors=FALSE)
    ))
    cat("NOTE: BOT_late_indvA reclassified as BOT_nc14b_rep3 (staging test H1 supported)\n")
    cat("      Add to Step1_SevenGeno_nc14b_DAR_Analysis.r BAM list — BOT group now has n=3\n\n")
  } else if ("BOT_late_indvA" %in% names(new_candidate_bams) &&
             "BOT_late_indvB" %in% names(new_candidate_bams)) {
    # BOT_late BAMs defined but staging test not yet resolved — hold as ambiguous
    rows <- c(rows, list(
      data.frame(sample="BOT_late_indvA",   genotype="BOT",    timepoint="unknown_late",
                 bam=new_candidate_bams[["BOT_late_indvA"]], status="staging_ambiguous",
                 stringsAsFactors=FALSE),
      data.frame(sample="BOT_late_indvB",   genotype="BOT",    timepoint="unknown_late",
                 bam=new_candidate_bams[["BOT_late_indvB"]], status="staging_ambiguous",
                 stringsAsFactors=FALSE)
    ))
  } else {
    # BOT_late BAMs not yet defined — skip (add paths to new_candidate_bams when ready)
    cat("NOTE: BOT_late_indvA/B not in new_candidate_bams — skipping staging rows.\n")
  }

  do.call(rbind, lapply(rows, function(r) {
    r$stringsAsFactors <- NULL; as.data.frame(r, stringsAsFactors = FALSE)
  }))
}

final_temporal_coldata <- build_final_temporal_coldata(staging_reclassify_indvA)

# Save coldata for downstream scripts
write.table(final_temporal_coldata,
            file.path(output_dir, "temporal_coldata_final.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Final temporal coldata →", file.path(output_dir, "temporal_coldata_final.txt"), "\n")
print(final_temporal_coldata[, c("sample","genotype","timepoint","status")])

################################################################################
# SECTION 6 — PEAK UNIVERSE EXTENSION
################################################################################
#
# The existing peak universe (FullUniverse_nc14b_AllGeno_union_peaks.bed) was
# built from nc14b samples only. Temporal samples at nc14d / late / gastr may
# access chromatin regions not open at nc14b.
#
# Strategy: extend the universe by unioning in the new timepoint peaks, but
# keep the original nc14b universe intact as a subset label so you can always
# restrict analyses to nc14b-accessible sites for apples-to-apples comparisons.
# ---------------------------------------------------------------------------

extend_peak_universe <- function(nc14b_bed_path, new_narrowpeak_paths, out_dir) {

  cat("\n=== EXTENDING PEAK UNIVERSE ===\n")

  new_nps <- new_narrowpeak_paths[file.exists(new_narrowpeak_paths)]
  if (length(new_nps) == 0) {
    cat("No new NarrowPeak files found — update new_narrowpeak_paths in Section 6.\n")
    return(invisible(NULL))
  }

  extended_bed <- file.path(out_dir, "FullUniverse_temporal_AllGeno_union_peaks.bed")
  cat("Merging", length(new_nps), "new peak files with nc14b universe...\n")

  # Shell call to bedtools (most efficient for large BEDs)
  np_str   <- paste(new_nps, collapse = " ")
  merge_cmd <- sprintf(
    "cat %s %s | awk 'BEGIN{OFS=\"\\t\"} !/^#/ && NF>=3 {print $1,$2,$3}' | sort -k1,1 -k2,2n | bedtools merge -d 200 | awk 'BEGIN{OFS=\"\\t\"} {print $1,$2,$3,\"temporal_union_\"NR,\".\",\".\"}' > %s",
    nc14b_bed_path, np_str, extended_bed
  )
  system(merge_cmd)

  orig_n    <- nrow(read.table(nc14b_bed_path, sep = "\t", header = FALSE,
                                comment.char = ""))
  extended_n <- nrow(read.table(extended_bed, sep = "\t", header = FALSE,
                                comment.char = ""))
  cat(sprintf("  nc14b universe:  %d peaks\n  Extended:        %d peaks (+%d new)\n",
              orig_n, extended_n, extended_n - orig_n))
  cat("  Extended universe →", extended_bed, "\n\n")
  extended_bed
}

# Paths to new NarrowPeak files — update to match your actual filenames
new_narrowpeak_paths <- c(
  # BOT nc14d — UPGRADED to confirmed pair
  "./ATAC_NarrowPeaks/D7_1_nc14d_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak",
  "./ATAC_NarrowPeaks/BOT_D7_2_nc14d_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak",
  "./ATAC_NarrowPeaks/Run-D7_Dm_ATAC_gastr_rep1_noq_rmdup_Macs2_peaks.narrowPeak",   # BOTR_gastr
  "./ATAC_NarrowPeaks/Run-D7_Dm_ATAC_gastr_rep2_noq_rmdup_Macs2_peaks.narrowPeak",   # BOT_hR_gastr rep1
  # BOT_hR_gastr rep2 — RunD7_3_gastr confirmed het (28368)
  "./ATAC_NarrowPeaks/RunD7_3_gastr_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak",
  "./ATAC_NarrowPeaks/HLH_B6_Nc14late_01_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak",
  "./ATAC_NarrowPeaks/Mat_Run_B6_4_nc14b_IR_q0.05_noq_rmdup_noChrM_Macs2_narrow_peaks.narrowPeak",
  "./ATAC_NarrowPeaks/Mat_Run_B6_3_late_noq_rmdup_Macs2_peaks.narrowPeak",
  "./ATAC_NarrowPeaks/Mat_Run_B6_4_late_noq_rmdup_Macs2_peaks.narrowPeak"
)

extended_universe_bed <- extend_peak_universe(
  nc14b_bed_path       = existing_peaks_bed,
  new_narrowpeak_paths = new_narrowpeak_paths,
  out_dir              = output_dir
)

################################################################################
# SECTION 7 — RECOMMENDED CONTRASTS FOR NEW TIMEPOINT GROUPS
# Now includes an extension_bp column: NA = use extension_bp_global (200)
################################################################################

cat("\n=== RECOMMENDED NEW CONTRASTS ===\n\n")

cat("NON-VENT TEMPORAL SERIES (once all candidates confirmed):\n")
new_contrasts <- data.frame(
  contrast_id  = 14:24,
  label        = c(
    "BOT_nc14d_vs_BOT_nc14b",
    "BOTR_gastr_vs_BOTR_nc14b",
    "BOT_hR_gastr_vs_BOT_hR_nc14b",
    "BOTR_gastr_vs_BOT_hR_gastr",
    "BOT_nc14d_vs_BOTR_nc14b",
    "BOTC_vs_BOTC_nc14late",
    "BOTC_oR_late_vs_BOTC_oR_nc14b",
    "BOTC_oR_late_vs_BOTC_nc14b",
    "BOT_late_indvA_vs_BOT_nc14b",
    "BOT_late_indvB_vs_BOTC_nc14b",
    "BOT_hR_gastr_rep2_vs_BOT_hR_nc14b"
  ),
  grpA = c("BOT_nc14d","BOTR_gastr","BOT_hR_gastr","BOTR_gastr",
           "BOT_nc14d","BOTC_nc14late","BOTC_oR_late","BOTC_oR_late",
           "BOT_late_indvA","BOT_late_indvB","BOT_hR_gastr_rep2"),
  grpB = c("BOT_nc14b","BOTR_nc14b","BOT_hR_nc14b","BOT_hR_gastr",
           "BOTR_nc14b","BOTC_nc14b","BOTC_oR_nc14b","BOTC_nc14b",
           "BOT_nc14b","BOTC_nc14b","BOT_hR_nc14b"),
  analysis_note = c(
    "temporal progression in WT non-vent",
    "Runt-null chromatin at gastrulation",
    "Runt-het chromatin at gastrulation (n=2 pair)",
    "Runt dosage at gastrulation (null vs het)",
    "Runt effect at nc14d: cross-timepoint",
    "BOTC temporal: nc14b → nc14late (rep2=interim proxy)",
    "Runt rescue progression to late stage",
    "late rescue vs non-rescue BOTC",
    "EXPLORATORY — pending staging call",
    "EXPLORATORY — pending staging call",
    "BOT_hR_gastr rep2 concordance check (28368)"
  ),
  method = c(
    "DESeq2+limma",   # BOT nc14d: UPGRADED to n=2
    "limma_only",     # BOTR_gastr: n=1
    "DESeq2+limma",   # BOT_hR_gastr: UPGRADED to n=2
    "limma_only",     # BOTR_gastr vs BOT_hR_gastr: BOTR n=1
    "limma_only",
    "DESeq2+limma",
    "DESeq2+limma","limma_only",
    "limma_only","limma_only",
    "limma_only"
  ),
  extension_bp = c(200, NA, NA, NA, 200, 200, 200, NA, NA, NA, NA),
  stringsAsFactors = FALSE
)

write.table(new_contrasts,
            file.path(output_dir, "recommended_temporal_contrasts.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(new_contrasts[, c("contrast_id","label","method","extension_bp","analysis_note")])

cat("\n",
    "NOTE on methods:\n",
    "  DESeq2+limma: both methods available (n>=2 per group in contrast)\n",
    "  limma_only:   use limma-voom only (n=1 in at least one group)\n",
    "  EXPLORATORY:  do not include in primary manuscript figures\n",
    "  BOT nc14d: UPGRADED to n=2 (28366 confirmed) → DESeq2+limma\n",
    "  BOT_hR_gastr: UPGRADED to n=2 (28368 confirmed het) → DESeq2+limma\n",
    "  BOTC nc14late: n=2 (rep2=Mat_Run_B6_4 nc14b interim proxy) → DESeq2+limma\n",
    "    Flag rep2 as interim; replace when confirmed nc14late BOTC rep available\n",
    "  BOTC_oR nc14b: n=2 → DESeq2+limma active\n\n")

################################################################################
# SECTION 8 — TEMPORAL DAR ANALYSIS — v3
# All five structural improvements active: L1 L2 L3 L4 L5 L7
################################################################################

# ── v3 helper functions (scoped inside Section 8) ───────────────────────────

.build_extended_boundary <- function(peaks_gr, ext_bp) {
  ext <- GRanges(seqnames=seqnames(peaks_gr),
    ranges=IRanges(start=pmax(1L,start(peaks_gr)-ext_bp), end=end(peaks_gr)+ext_bp),
    strand="*")
  names(ext) <- names(peaks_gr); ext <- trim(ext)
  merged <- reduce(ext, with.revmap=TRUE)
  orig_n <- names(peaks_gr)
  names(merged) <- paste0("bext",ext_bp,"_",
    sapply(mcols(merged)$revmap, function(i) orig_n[i[1]]))
  merged
}

.compute_tier1_thresholds <- function(counts_all, grpA_idx, grpB_idx) {
  lib  <- colSums(counts_all[,c(grpA_idx,grpB_idx),drop=FALSE])
  mlib <- median(lib)
  list(
    min_reads  = max(1L, round(min_cpm_present * mlib / 1e6)),
    max_reads  = max(1L, round(max_cpm_absent  * mlib / 1e6 *
                                min(length(grpA_idx),length(grpB_idx))))
  )
}

################################################################################
# CACHE COLUMN ALIGNMENT
# Maps precount superset keys (NV_COLS names in Precount_AllCounts_Matrices_v2.r)
# → avail$sample values (build_final_temporal_coldata sample column).
# These differ because the precount script uses shorter nc14b key names
# (e.g. "BOTv_rep1") while this script uses timepoint-qualified names
# (e.g. "BOTv_nc14b_rep1").
################################################################################

CACHE_KEY_MAP <- c(
  # nc14b: precount uses short keys, Step1 uses timepoint-qualified keys
  BOTv_rep1          = "BOTv_nc14b_rep1",
  BOTv_rep2          = "BOTv_nc14b_rep2",
  BOTCv_rep1         = "BOTCv_nc14b_rep1",     # not in final_temporal_coldata but harmless
  BOTCv_rep2         = "BOTCv_nc14b_rep2",
  BOT_nc14b_rep1     = "BOT_nc14b_rep1",
  BOT_nc14b_rep2     = "BOT_nc14b_rep2",
  BOT_hR_rep1        = "BOT_hR_nc14b_rep1",
  BOT_hR_rep2        = "BOT_hR_nc14b_rep2",
  BOTR_rep1          = "BOTR_nc14b_rep1",
  BOTR_rep2          = "BOTR_nc14b_rep2",
  BOTC_rep1          = "BOTC_nc14b_rep1",
  BOTC_rep2          = "BOTC_nc14b_rep2",
  BOTC_oR_rep1       = "BOTC_oR_nc14b_rep1",
  BOTC_oR_rep2       = "BOTC_oR_nc14b_rep2",
  # New timepoints: precount and Step1 keys already match
  BOTR_gastr_rep1    = "BOTR_gastr_rep1",
  BOT_hR_gastr_rep1  = "BOT_hR_gastr_rep1",
  BOT_hR_gastr_rep2  = "BOT_hR_gastr_rep2",
  BOT_nc14d_rep1     = "BOT_nc14d_rep1",
  BOT_nc14d_rep2     = "BOT_nc14d_rep2",
  BOTC_nc14late_rep1 = "BOTC_nc14late_rep1",
  BOTC_nc14late_rep2 = "BOTC_nc14late_rep2",
  BOTC_oR_late_rep1  = "BOTC_oR_late_rep1",
  BOTC_oR_late_rep2  = "BOTC_oR_late_rep2"
)

realign_cache_cols <- function(counts_mat, target_names,
                                key_map = CACHE_KEY_MAP) {
  if (is.null(counts_mat)) return(NULL)
  cache_cols <- colnames(counts_mat)
  if (identical(cache_cols, target_names)) return(counts_mat)
  matched <- match(target_names, cache_cols)
  if (!anyNA(matched)) return(counts_mat[, matched, drop=FALSE])
  translated <- key_map[cache_cols]
  known      <- !is.na(translated)
  colnames(counts_mat)[known] <- translated[known]
  matched <- match(target_names, colnames(counts_mat))
  if (anyNA(matched)) {
    warning("realign_cache_cols: unmatched after key_map translation:\n  ",
            paste(target_names[is.na(matched)], collapse="\n  "),
            "\nCache columns (translated): ",
            paste(colnames(counts_mat), collapse=", "))
    return(NULL)
  }
  counts_mat[, matched, drop=FALSE]
}

realign_grp_cache <- function(cache_obj, target_names,
                               key_map = CACHE_KEY_MAP) {
  if (is.null(cache_obj)) return(NULL)
  cache_obj$counts <- lapply(cache_obj$counts, realign_cache_cols,
                              target_names=target_names, key_map=key_map)
  cache_obj
}

.load_group_peaks <- function(bamlist, sample_col) {
  cat("  [L1] Loading group-specific peak counts...\n")
  gc <- list(); gp <- list()
  for (grp in names(group_peak_beds)) {
    bed <- group_peak_beds[[grp]]
    if (!file.exists(bed)) { gc[[grp]]<-NULL; gp[[grp]]<-NULL; next }
    d  <- read.table(bed,sep="\t",header=FALSE,stringsAsFactors=FALSE,comment.char="")
    gr <- GRanges(seqnames=d[,1],ranges=IRanges(start=d[,2]+1,end=d[,3]),strand="*")
    if (ncol(d)>=4) names(gr)<-d[,4] else names(gr)<-paste0("peak_",seq_along(gr))
    saf <- data.frame(GeneID=names(gr),Chr=as.character(seqnames(gr)),
      Start=start(gr),End=end(gr),Strand=".",stringsAsFactors=FALSE)
    saf <- saf[!apply(saf,1,function(x)any(is.na(x))),]
    saf$GeneID <- make.unique(saf$GeneID)
    saf$Start  <- as.integer(saf$Start); saf$End <- as.integer(saf$End)
    saf <- saf[saf$Start>=1&saf$End>saf$Start&!is.na(saf$Start)&!is.na(saf$End),]
    fc <- featureCounts(files=bamlist,annot.ext=saf,isPairedEnd=TRUE,
                        countMultiMappingReads=FALSE,primaryOnly=TRUE,nthreads=FC_THREADS)
    cnts <- fc$counts; colnames(cnts) <- sample_col
    gr   <- gr[names(gr)%in%saf$GeneID]; rownames(cnts) <- names(gr)
    gc[[grp]] <- cnts; gp[[grp]] <- gr
    cat("   ",grp,":",length(gr),"peaks\n")
  }
  list(counts=gc, peaks=gp)
}

.ruv_pairwise <- function(counts, grpA_idx, grpB_idx, grpA_label, grpB_label) {
  cat("    [L3] Pairwise RUV for",grpA_label,"vs",grpB_label,"...\n")
  pair   <- c(grpA_idx, grpB_idx)
  cnts_p <- counts[,pair,drop=FALSE]
  grp_p  <- factor(c(rep(grpA_label,length(grpA_idx)),rep(grpB_label,length(grpB_idx))))
  pm     <- rowMeans(cnts_p); cv <- apply(cnts_p,1,sd)/(pm+0.5)
  dge    <- DGEList(counts=cnts_p,group=grp_p); dge <- calcNormFactors(dge)
  des    <- model.matrix(~grp_p)
  dge    <- estimateDisp(dge,des); fit <- glmQLFit(dge,des)
  res    <- topTags(glmQLFTest(fit,coef=2),n=Inf)$table
  non_de <- which(res$PValue>0.5&res$logCPM>1)
  stable <- which(pm>quantile(pm,0.4)&cv<quantile(cv[pm>quantile(pm,0.4)],0.3))
  int    <- intersect(stable,non_de)
  ctrl   <- if(length(int)>=50) int else if(length(stable)>=50) stable
            else if(length(non_de)>=50) non_de
            else order(cv[pm>quantile(pm,0.2)])[seq_len(min(200,max(100,round(nrow(cnts_p)*0.1))))]
  cat("      Controls:",length(ctrl),"\n")
  rownames(counts)[ctrl]
}

.tier1_grp <- function(gp_A, gc_A, gp_B, gc_B, union_counts,
                        iA, iB, lA, lB, min_r, max_r, mode_lbl, cname) {
  # iA / iB are integer positions in union_counts (full sample set).
  # gc_A / gc_B only contain columns for that genotype's BAMs — subset by name.
  sA <- colnames(union_counts)[iA]
  sB <- colnames(union_counts)[iB]
  sA_in_A <- sA[sA %in% colnames(gc_A)]; sB_in_A <- sB[sB %in% colnames(gc_A)]
  sB_in_B <- sB[sB %in% colnames(gc_B)]; sA_in_B <- sA[sA %in% colnames(gc_B)]
  A_in_A <- if (length(sA_in_A)>0) rowSums(gc_A[,sA_in_A,drop=FALSE]) else rep(0,nrow(gc_A))
  A_in_B <- if (length(sB_in_A)>0) rowSums(gc_A[,sB_in_A,drop=FALSE]) else rep(0,nrow(gc_A))
  B_in_B <- if (length(sB_in_B)>0) rowSums(gc_B[,sB_in_B,drop=FALSE]) else rep(0,nrow(gc_B))
  B_in_A <- if (length(sA_in_B)>0) rowSums(gc_B[,sA_in_B,drop=FALSE]) else rep(0,nrow(gc_B))
  lib_e  <- median(colSums(union_counts[,c(iA,iB),drop=FALSE]))+1
  mk <- function(peaks_gr,present,absent,n_p,lfc,lbl) {
    idx <- present>=min_r & absent<=max_r
    if (sum(idx)==0) return(data.frame(peak_id=character(),chr=character(),
      start=integer(),end=integer(),log2FC=numeric(),pvalue=numeric(),padj=numeric(),
      baseMean=numeric(),AveExpr=numeric(),method=character(),contrast=character(),
      mode=character(),stringsAsFactors=FALSE))
    pm <- present[idx]/max(n_p,1)
    data.frame(peak_id=names(peaks_gr)[idx],chr=as.character(seqnames(peaks_gr[idx])),
      start=start(peaks_gr[idx])-1,end=end(peaks_gr[idx]),log2FC=lfc,pvalue=1e-10,
      padj=1e-10,baseMean=pm,AveExpr=log2((pm/lib_e*1e6)+0.5),
      method=paste0(mode_lbl,"_",lbl,"_specific_",cname),
      contrast=cname,mode=mode_lbl,stringsAsFactors=FALSE)
  }
  list(A_specific=mk(gp_A,A_in_A,A_in_B,length(iA),-5,paste0(lA,"_open_lost")),
       B_specific=mk(gp_B,B_in_B,B_in_A,length(iB), 5,paste0(lB,"_open_gained")))
}

.df_to_gr_v3 <- function(df, tier="differential", ma_map=NULL) {
  if (nrow(df)==0) return(GRanges())
  eps <- 1e-300
  ts  <- switch(tier,"condition_specific"=3,"differential"=ifelse(grepl("^NARROW",df$mode[1]),2,1),1)
  ma  <- if (tier=="condition_specific") rep(1.0,nrow(df))
         else if (!is.null(ma_map)&&length(ma_map)>0) {
           v<-ma_map[df$peak_id]; v[is.na(v)]<-1.0; as.numeric(v)
         } else rep(1.0,nrow(df))
  pv  <- ifelse(is.na(df$pvalue),1,df$pvalue)
  GRanges(seqnames=df$chr,ranges=IRanges(start=df$start+1,end=df$end),
    log2FC=df$log2FC,pvalue=pv,
    padj    =if("padj"    %in%names(df))df$padj    else NA_real_,
    baseMean=if("baseMean"%in%names(df))df$baseMean else NA_real_,
    AveExpr =if("AveExpr" %in%names(df))df$AveExpr  else NA_real_,
    method  =if("method"  %in%names(df))df$method   else tier,
    contrast=if("contrast"%in%names(df))df$contrast else NA_character_,
    tier=tier,tier_score=ts,m_agree=ma,
    conf_score=ma*ts*(-log10(pv+eps)),
    mode=if("mode"%in%names(df))df$mode else tier)
}

.dedup_v3 <- function(gr) {
  if (length(gr)<=1) return(gr)
  mcols(gr)$sign_dir <- sign(mcols(gr)$log2FC)
  one <- function(g) {
    if (length(g)<=1) return(g)
    ov <- as.data.frame(findOverlaps(g,g,minoverlap=1L))
    ov <- ov[ov$queryHits!=ov$subjectHits,]
    if (nrow(ov)==0) return(g)
    p <- seq_len(length(g))
    fr <- function(x){while(p[x]!=x){p[x]<<-p[p[x]];x<-p[x]};x}
    for (i in seq_len(nrow(ov))){a<-fr(ov$queryHits[i]);b<-fr(ov$subjectHits[i]);if(a!=b)p[a]<-b}
    cl<-sapply(seq_len(length(g)),fr); sc<-mcols(g)$conf_score
    keep<-sapply(unique(cl),function(c){m<-which(cl==c);m[which.max(sc[m])]})
    g[sort(keep)]
  }
  c(one(gr[mcols(gr)$sign_dir>0]),one(gr[mcols(gr)$sign_dir<0]))
}

.diag_tier3 <- function(cname, all_gr, peaks_n, grp_pA, grp_pB, ext_bp, out_dir) {
  t3 <- all_gr[mcols(all_gr)$tier=="differential"&grepl("EXTENDED",mcols(all_gr)$mode)]
  cat("  [L7] Tier3:",length(t3),"extended-only DARs\n")
  if (length(t3)==0) return(invisible(NULL))
  covered <- countOverlaps(t3,peaks_n,minoverlap=1L)>0
  in_A    <- if(!is.null(grp_pA)) countOverlaps(t3,grp_pA,minoverlap=1L)>0 else rep(FALSE,length(t3))
  in_B    <- if(!is.null(grp_pB)) countOverlaps(t3,grp_pB,minoverlap=1L)>0 else rep(FALSE,length(t3))
  fl      <- GRanges(seqnames=seqnames(t3),
    ranges=IRanges(start=pmax(1L,start(t3)-ext_bp),end=end(t3)+ext_bp),strand="*")
  n_adj   <- countOverlaps(fl,peaks_n,minoverlap=1L)
  cls     <- dplyr::case_when(covered&(in_A|in_B)~"tight_boundary",
               n_adj>=2~"adjacent_merged",covered~"tight_boundary",TRUE~"diffuse")
  df <- data.frame(chr=as.character(seqnames(t3)),start=start(t3)-1,end=end(t3),
    log2FC=mcols(t3)$log2FC,pvalue=mcols(t3)$pvalue,conf_score=mcols(t3)$conf_score,
    direction=ifelse(mcols(t3)$log2FC>0,"gained-open","gained-close"),
    has_narrow_peak=covered,in_grpA=in_A,in_grpB=in_B,n_adjacent=n_adj,
    classification=cls,stringsAsFactors=FALSE)
  write.table(df,file.path(out_dir,paste0(cname,"_Tier3_diagnostic.txt")),
              sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
  cat("    Classification:",paste(names(table(cls)),table(cls),sep="=",collapse=" | "),"\n")
  invisible(df)
}

.write_dar_v3 <- function(gr, tag, out_dir) {
  if (length(gr)==0) { cat("  Skip",tag,"(0)\n"); return(invisible(NULL)) }
  mc <- mcols(gr)
  df <- data.frame(
    chr=as.character(seqnames(gr)),start=start(gr)-1,end=end(gr),
    name=paste0(tag,"_",seq_along(gr)),
    score=round(-log10(mc$pvalue+1e-300)*10),strand=".",
    log2FC=mc$log2FC,pvalue=mc$pvalue,
    padj      =if("padj"      %in%names(mc))mc$padj       else NA_real_,
    baseMean  =if("baseMean"  %in%names(mc))mc$baseMean   else NA_real_,
    AveExpr   =if("AveExpr"   %in%names(mc))mc$AveExpr    else NA_real_,
    conf_score=if("conf_score"%in%names(mc))mc$conf_score else NA_real_,
    m_agree   =if("m_agree"   %in%names(mc))mc$m_agree    else NA_real_,
    direction =ifelse(mc$log2FC>0,"gained-open","gained-close"),
    tier      =if("tier"      %in%names(mc))mc$tier       else NA_character_,
    tier_score=if("tier_score"%in%names(mc))mc$tier_score else NA_real_,
    mode      =if("mode"      %in%names(mc))mc$mode       else NA_character_,
    method    =if("method"    %in%names(mc))mc$method     else NA_character_,
    contrast  =if("contrast"  %in%names(mc))mc$contrast   else NA_character_,
    stringsAsFactors=FALSE)
  write.table(df[,1:6],file.path(out_dir,paste0(tag,"_DARs.bed")),
              sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
  write.table(df,file.path(out_dir,paste0(tag,"_DARs_annotated.txt")),
              sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
  cat("  Saved",nrow(df),"→",paste0(tag,"_DARs.bed\n"))
}

run_temporal_dar <- function(final_coldata, peaks_bed, contrasts_df,
                              ext_bp_global = extension_bp_global,
                              out_dir       = output_dir,
                              t1_dir        = tier1_dir,
                              d_dir         = diag_dir) {

  cat("\n=== TEMPORAL DAR ANALYSIS v4 (L1-L5+L7, cache-accelerated) ===\n\n")
  cat("  Working directory:", getwd(), "\n")
  cat("  Cache directory:  ", CACHE_DIR, "\n")
  cat("  Cache files present:\n")
  cache_files <- list.files(CACHE_DIR, full.names=FALSE)
  if (length(cache_files) == 0) cat("    (none)\n") else
    cat(paste0("    ", cache_files, "\n"), sep="")

  peaks_all <- local({
    d  <- read.table(peaks_bed,sep="\t",header=FALSE,stringsAsFactors=FALSE,comment.char="")
    gr <- GRanges(seqnames=d[,1],ranges=IRanges(start=d[,2]+1,end=d[,3]),strand="*")
    if(ncol(d)>=4) names(gr)<-d[,4] else names(gr)<-paste0("peak_",seq_along(gr))
    gr
  })

  avail   <- final_coldata[file.exists(final_coldata$bam),]
  bamlist <- avail$bam
  if (nrow(avail)<4) {
    cat("Fewer than 4 BAMs available — deferred.\n"); return(invisible(NULL))
  }

  # ---------------------------------------------------------------------------
  # Narrow counts — load from cache, count only new peaks if needed, merge
  # ---------------------------------------------------------------------------
  narrow_cache_file <- file.path(CACHE_DIR, "narrow_counts.txt")
  new_peaks_cache   <- file.path(CACHE_DIR, "narrow_counts_newpeaks.txt")

  if (file.exists(narrow_cache_file)) {
    cat("  Loading narrow counts from cache...\n")
    counts_cached <- as.matrix(read.table(narrow_cache_file,
                                           header=TRUE, row.names=1, sep="\t"))
    counts_cached <- realign_cache_cols(counts_cached, avail$sample)
    if (is.null(counts_cached))
      stop("Could not align narrow counts cache columns to avail$sample.")

    # Find peaks in peaks_all that are not yet in the cache
    new_peak_names <- setdiff(names(peaks_all), rownames(counts_cached))

    if (length(new_peak_names) == 0) {
      cat("  All", length(peaks_all), "peaks found in cache.\n")
      counts_n <- counts_cached[names(peaks_all), , drop=FALSE]
    } else {
      cat("  Cache has", nrow(counts_cached), "peaks;",
          length(new_peak_names), "new peaks need counting.\n")

      # Load or compute counts for new peaks only
      if (file.exists(new_peaks_cache)) {
        cat("  Loading new-peak counts from cache...\n")
        counts_new <- as.matrix(read.table(new_peaks_cache,
                                            header=TRUE, row.names=1, sep="\t"))
        counts_new <- realign_cache_cols(counts_new, avail$sample)
      } else {
        cat("  Counting", length(new_peak_names), "new peaks only",
            "(much faster than full recount)...\n")
        new_peaks_gr <- peaks_all[new_peak_names]
        saf_new <- data.frame(
          GeneID = names(new_peaks_gr),
          Chr    = as.character(seqnames(new_peaks_gr)),
          Start  = start(new_peaks_gr),
          End    = end(new_peaks_gr),
          Strand = ".", stringsAsFactors = FALSE)
        saf_new <- saf_new[!apply(saf_new, 1, anyNA), ]
        saf_new$Start <- as.integer(saf_new$Start)
        saf_new$End   <- as.integer(saf_new$End)
        saf_new <- saf_new[saf_new$Start >= 1 & saf_new$End > saf_new$Start, ]
        fc_new <- featureCounts(files=bamlist, annot.ext=saf_new,
                                isPairedEnd=TRUE, countMultiMappingReads=FALSE,
                                primaryOnly=TRUE, nthreads=FC_THREADS)
        counts_new <- fc_new$counts
        colnames(counts_new) <- avail$sample
        rownames(counts_new) <- saf_new$GeneID
        write.table(counts_new, new_peaks_cache,
                    sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
        cat("  New-peak counts saved →", new_peaks_cache, "\n")
      }

      # Merge: cached peaks + new peaks, ordered to match peaks_all
      counts_n <- rbind(counts_cached, counts_new)
      counts_n <- counts_n[names(peaks_all)[names(peaks_all) %in% rownames(counts_n)],
                           , drop=FALSE]
      cat("  Merged:", nrow(counts_n), "total peaks.\n")
      # Promote merged matrix to main cache so future runs skip split-file logic
      write.table(counts_n, narrow_cache_file,
                  sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
      cat("  Promoted merged counts → ", narrow_cache_file, "\n")
    }
    peaks_n_all <- peaks_all[rownames(counts_n)]

  } else {
    cat("  Cache not found — running featureCounts for all peaks (first run).\n")
    saf <- data.frame(GeneID=names(peaks_all), Chr=as.character(seqnames(peaks_all)),
      Start=start(peaks_all), End=end(peaks_all), Strand=".", stringsAsFactors=FALSE)
    fc  <- featureCounts(files=bamlist, annot.ext=saf, isPairedEnd=TRUE,
                         countMultiMappingReads=FALSE, primaryOnly=TRUE, nthreads=FC_THREADS)
    counts_n <- fc$counts; colnames(counts_n) <- avail$sample
    rownames(counts_n) <- names(peaks_all)
    dir.create(CACHE_DIR, recursive=TRUE, showWarnings=FALSE)
    write.table(counts_n, narrow_cache_file,
                sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
    cat("  Saved:", narrow_cache_file, "\n")
    peaks_n_all <- peaks_all
  }
  kp_n     <- rowSums(counts_n >= 5) >= 2
  counts_n <- counts_n[kp_n, ]; peaks_n <- peaks_n_all[kp_n]

  # ---------------------------------------------------------------------------
  # [L1] Group-specific peak counts — load from cache or count on first run
  # ---------------------------------------------------------------------------
  grp_cache_file <- file.path(CACHE_DIR, "grp_counts_cache.rds")
  if (file.exists(grp_cache_file)) {
    cat("  [L1] Loading group-specific peak counts from cache...\n")
    grp_cache <- readRDS(grp_cache_file)
    grp_cache <- realign_grp_cache(grp_cache, avail$sample)
  } else {
    cat("  [L1] Cache not found — counting group peaks (first run).\n")
    grp_cache <- .load_group_peaks(bamlist, avail$sample)
    dir.create(CACHE_DIR, recursive=TRUE, showWarnings=FALSE)
    saveRDS(grp_cache, grp_cache_file)
    cat("  Saved:", grp_cache_file, "\n")
  }

  # ---------------------------------------------------------------------------
  # [L4] Extended union peak counts — load from cache or count on first run
  # ---------------------------------------------------------------------------
  cat("  [L4] Loading boundary-extended peak counts...\n")
  ext_bps <- sort(unique(c(ext_bp_global,
    contrasts_df$extension_bp[!is.na(contrasts_df$extension_bp)])))
  ext_cache <- list()
  for (ebp in ext_bps) {
    key        <- as.character(ebp)
    cache_file <- file.path(CACHE_DIR, paste0("ext", ebp, "_counts.rds"))
    if (file.exists(cache_file)) {
      cat("   Loading ext", ebp, "bp from cache\n")
      cached_ext <- readRDS(cache_file)
      cached_ext$counts <- realign_cache_cols(cached_ext$counts, avail$sample)
      ext_cache[[key]] <- cached_ext
    } else {
      cat("   Cache missing for ext", ebp, "bp — counting (first run).\n")
      p_ext <- .build_extended_boundary(peaks_all, ebp)
      saf_e <- data.frame(GeneID=names(p_ext),Chr=as.character(seqnames(p_ext)),
        Start=start(p_ext),End=end(p_ext),Strand=".",stringsAsFactors=FALSE)
      saf_e <- saf_e[!apply(saf_e,1,function(x)any(is.na(x))),]
      saf_e$GeneID <- make.unique(saf_e$GeneID)
      saf_e$Start <- as.integer(saf_e$Start); saf_e$End <- as.integer(saf_e$End)
      saf_e <- saf_e[saf_e$Start>=1&saf_e$End>saf_e$Start,]
      fc_e  <- featureCounts(files=bamlist,annot.ext=saf_e,isPairedEnd=TRUE,
                             countMultiMappingReads=FALSE,primaryOnly=TRUE,nthreads=FC_THREADS)
      cnts_e <- fc_e$counts; colnames(cnts_e) <- avail$sample
      p_ext  <- p_ext[names(p_ext)%in%saf_e$GeneID]; rownames(cnts_e) <- names(p_ext)
      kp_e   <- rowSums(cnts_e>=5)>=2
      ext_cache[[key]] <- list(counts=cnts_e[kp_e,],peaks=p_ext[kp_e])
      dir.create(CACHE_DIR, recursive=TRUE, showWarnings=FALSE)
      saveRDS(ext_cache[[key]], cache_file)
      cat("   Saved:", cache_file, "\n")
    }
    cat("   ", ebp, "bp:", nrow(ext_cache[[key]]$counts), "peaks\n")
  }

  # ---------------------------------------------------------------------------
  # [L1+L4] Extended group-specific peak counts — load from cache
  # Eliminates repeated featureCounts calls per contrast for Tier 1 extended.
  # ---------------------------------------------------------------------------
  cat("  [L1+L4] Loading extended group-specific peak counts...\n")
  grp_ext_cache <- list()
  for (ebp in ext_bps) {
    key        <- as.character(ebp)
    cache_file <- file.path(CACHE_DIR, paste0("grp_ext", ebp, "_counts.rds"))
    if (file.exists(cache_file)) {
      cat("   Loading grp_ext", ebp, "bp from cache\n")
      gec <- readRDS(cache_file)
      gec <- realign_grp_cache(gec, avail$sample)
      grp_ext_cache[[key]] <- gec
    } else {
      cat("   Cache missing for grp_ext", ebp, "bp — counting (first run).\n")
      grp_counts_ext <- list(); grp_peaks_ext <- list()
      for (grp in names(group_peak_beds)) {
        bed <- group_peak_beds[[grp]]
        if (!file.exists(bed)) next
        d   <- read.table(bed,sep="\t",header=FALSE,stringsAsFactors=FALSE,comment.char="")
        gp  <- GRanges(seqnames=d[,1],ranges=IRanges(start=d[,2]+1,end=d[,3]),strand="*")
        if (ncol(d)>=4) names(gp)<-d[,4] else names(gp)<-paste0("peak_",seq_along(gp))
        gp_ext <- .build_extended_boundary(gp, ebp)
        saf_g  <- data.frame(GeneID=names(gp_ext),Chr=as.character(seqnames(gp_ext)),
          Start=start(gp_ext),End=end(gp_ext),Strand=".",stringsAsFactors=FALSE)
        saf_g  <- saf_g[!apply(saf_g,1,function(x)any(is.na(x))),]
        saf_g$GeneID<-make.unique(saf_g$GeneID)
        saf_g$Start<-as.integer(saf_g$Start); saf_g$End<-as.integer(saf_g$End)
        saf_g  <- saf_g[saf_g$Start>=1&saf_g$End>saf_g$Start,]
        fc_g   <- featureCounts(files=bamlist,annot.ext=saf_g,isPairedEnd=TRUE,
                                countMultiMappingReads=FALSE,primaryOnly=TRUE,nthreads=FC_THREADS)
        cnts_g <- fc_g$counts; colnames(cnts_g) <- avail$sample
        gp_ext <- gp_ext[names(gp_ext)%in%saf_g$GeneID]; rownames(cnts_g) <- names(gp_ext)
        grp_counts_ext[[grp]] <- cnts_g; grp_peaks_ext[[grp]] <- gp_ext
      }
      grp_ext_cache[[key]] <- list(counts=grp_counts_ext, peaks=grp_peaks_ext)
      dir.create(CACHE_DIR, recursive=TRUE, showWarnings=FALSE)
      saveRDS(grp_ext_cache[[key]], cache_file)
      cat("   Saved:", cache_file, "\n")
    }
  }
  cat("\n")

  dar_list <- list()

  for (i in seq_len(nrow(contrasts_df))) {
    cmp <- contrasts_df[i,]
    cat(strrep("-",50),"\n"); cat("CONTRAST:",cmp$label,"\n")

    iA <- which(grepl(cmp$grpA,avail$sample,fixed=TRUE))
    iB <- which(grepl(cmp$grpB,avail$sample,fixed=TRUE))
    if (length(iA)==0||length(iB)==0) {
      cat("  Skipping — samples not yet available\n\n"); next
    }

    ebp    <- if (!is.na(cmp$extension_bp)) cmp$extension_bp else ext_bp_global
    cached <- ext_cache[[as.character(ebp)]]
    iA_b   <- which(grepl(cmp$grpA,colnames(cached$counts),fixed=TRUE))
    iB_b   <- which(grepl(cmp$grpB,colnames(cached$counts),fixed=TRUE))
    cat("  extension_bp =",ebp,"(boundary, L4)\n")

    # [L2] CPM thresholds
    thr <- .compute_tier1_thresholds(counts_n, iA, iB)
    cat(sprintf("  Tier1 thresholds: min_reads=%d  max_reads=%d\n",
                thr$min_reads, thr$max_reads))

    # [L1] Group-specific Tier 1
    # Parse genotype from group label (e.g. "BOT_nc14d" → "BOT")
    parse_geno <- function(grp_label) {
      # Match against known genotype names longest-first
      genos <- c("BOTC_oR","BOT_hR","BOTCv","BOTv","BOTC","BOTR","BOT")
      for (g in genos) if (startsWith(grp_label,g)) return(g)
      return(grp_label)
    }
    gA <- parse_geno(cmp$grpA); gB <- parse_geno(cmp$grpB)
    gpA <- grp_cache$peaks[[gA]];  gcA <- grp_cache$counts[[gA]]
    gpB <- grp_cache$peaks[[gB]];  gcB <- grp_cache$counts[[gB]]

    if (is.null(gcA)||is.null(gcB)) {
      cat("  WARNING: group peaks missing — falling back for Tier 1\n")
      nc <- list(A_specific=data.frame(peak_id=character(),chr=character(),
        start=integer(),end=integer(),log2FC=numeric(),pvalue=numeric(),padj=numeric(),
        baseMean=numeric(),AveExpr=numeric(),method=character(),contrast=character(),
        mode=character(),stringsAsFactors=FALSE),B_specific=data.frame())
      nc$B_specific <- nc$A_specific; bc <- nc
    } else {
      nc <- .tier1_grp(gpA,gcA,gpB,gcB,counts_n,iA,iB,
                       cmp$grpA,cmp$grpB,thr$min_reads,thr$max_reads,"NARROW",cmp$label)

      # [L1+L4] Extended Tier 1 — load from grp_ext_cache (v4)
      # Replaces per-contrast featureCounts calls for gpA_ext / gpB_ext.
      ebp_key <- as.character(ebp)
      gpA_ext  <- grp_ext_cache[[ebp_key]]$peaks[[gA]]
      cnts_Ae  <- grp_ext_cache[[ebp_key]]$counts[[gA]]
      gpB_ext  <- grp_ext_cache[[ebp_key]]$peaks[[gB]]
      cnts_Be  <- grp_ext_cache[[ebp_key]]$counts[[gB]]

      if (is.null(cnts_Ae) || is.null(cnts_Be)) {
        cat("  WARNING: grp_ext cache missing for", gA, "/", gB,
            "— falling back to featureCounts\n")
        gpA_ext <- .build_extended_boundary(gpA, ebp)
        gpB_ext <- .build_extended_boundary(gpB, ebp)
        saf_Ae <- data.frame(GeneID=names(gpA_ext),Chr=as.character(seqnames(gpA_ext)),
          Start=start(gpA_ext),End=end(gpA_ext),Strand=".",stringsAsFactors=FALSE)
        saf_Be <- data.frame(GeneID=names(gpB_ext),Chr=as.character(seqnames(gpB_ext)),
          Start=start(gpB_ext),End=end(gpB_ext),Strand=".",stringsAsFactors=FALSE)
        for (saf_tmp in list(saf_Ae, saf_Be)) {
          saf_tmp <- saf_tmp[!apply(saf_tmp,1,function(x)any(is.na(x))),]
          saf_tmp$GeneID<-make.unique(saf_tmp$GeneID)
          saf_tmp$Start<-as.integer(saf_tmp$Start); saf_tmp$End<-as.integer(saf_tmp$End)
        }
        fc_Ae <- featureCounts(files=bamlist,annot.ext=saf_Ae,isPairedEnd=TRUE,
                                countMultiMappingReads=FALSE,primaryOnly=TRUE,nthreads=FC_THREADS)
        fc_Be <- featureCounts(files=bamlist,annot.ext=saf_Be,isPairedEnd=TRUE,
                                countMultiMappingReads=FALSE,primaryOnly=TRUE,nthreads=FC_THREADS)
        cnts_Ae <- fc_Ae$counts; colnames(cnts_Ae) <- avail$sample
        cnts_Be <- fc_Be$counts; colnames(cnts_Be) <- avail$sample
        gpA_ext <- gpA_ext[names(gpA_ext)%in%saf_Ae$GeneID]
        gpB_ext <- gpB_ext[names(gpB_ext)%in%saf_Be$GeneID]
      }

      bc <- .tier1_grp(gpA_ext,cnts_Ae,gpB_ext,cnts_Be,counts_n,iA,iB,
                       cmp$grpA,cmp$grpB,thr$min_reads,thr$max_reads,
                       paste0("EXTENDED",ebp),cmp$label)
    }
    cat("  Tier1 narrow:",nrow(nc$A_specific),"+",nrow(nc$B_specific),"\n")
    cat("  Tier1 ext",ebp,":",nrow(bc$A_specific),"+",nrow(bc$B_specific),"\n")
    write_tier1_outputs(nc,bc,cmp$label,t1_dir)

    # [L3] Pairwise RUV controls
    ruv_n <- .ruv_pairwise(counts_n,    iA,   iB,   cmp$grpA, cmp$grpB)
    ruv_b <- .ruv_pairwise(cached$counts,iA_b,iB_b, cmp$grpA, cmp$grpB)

    # Tier 2/3 — limma-voom (+ DESeq2 where n>=2 in both groups)
    run_s8 <- function(cnts, p_gr, mode_lbl, iA_i, iB_i, ruv_ctrl) {
      cdat  <- avail[c(iA_i,iB_i),]
      cnt_s <- round(cnts[,c(iA_i,iB_i)])
      valid_ctrl <- intersect(ruv_ctrl, rownames(cnt_s))
      use_ruv    <- length(valid_ctrl) >= 10
      # For singleton groups (n=1 per group) RUV W_1 still saturates the model
      # (2 group params + 1 W = 3 params, 2 samples = -1 df). Skip RUV if either
      # group has only 1 sample.
      n_per_grp <- min(sum(cdat$sample %in% avail$sample[iA_i]),
                       sum(cdat$sample %in% avail$sample[iB_i]))
      if (use_ruv && n_per_grp < 2) use_ruv <- FALSE
      # Temporal contrasts compare same genotype across timepoints — genotype is
      # constant so use timepoint as the grouping variable instead.
      n_geno <- length(unique(cdat$genotype))
      grp_col <- if (n_geno < 2) "timepoint" else "genotype"
      cdat$grp_fac <- as.factor(cdat[[grp_col]])
      if (use_ruv) {
        set     <- newSeqExpressionSet(counts=as.matrix(cnt_s),
                     phenoData=data.frame(condition=cdat$sample,row.names=colnames(cnt_s)))
        set_ruv <- RUVg(set,valid_ctrl,k=1)
        ruv_f   <- pData(set_ruv)[,grep("^W_",colnames(pData(set_ruv))),drop=FALSE]
        cdat    <- cbind(cdat,ruv_f); cnt_s <- counts(set_ruv)
        df_f    <- as.formula(paste("~ 0 + grp_fac +",
                     paste(paste0("W_",seq_len(ncol(ruv_f))),collapse="+")))
      } else {
        df_f <- ~ 0 + grp_fac
      }
      dge <- DGEList(counts=cnt_s, group=cdat$grp_fac)
      dge <- calcNormFactors(dge, method="TMM")
      des <- model.matrix(df_f, data=cdat)
      # Strip "grp_fac" prefix from colnames so contrast vector assignment works
      colnames(des) <- sub("^grp_fac", grp_col, colnames(des))
      cv_vec <- setNames(rep(0,ncol(des)),colnames(des))
      a_vals <- unique(cdat[[grp_col]][cdat$sample %in% avail$sample[iA_i]])
      b_vals <- unique(cdat[[grp_col]][cdat$sample %in% avail$sample[iB_i]])
      a_nms  <- paste0(grp_col, a_vals)
      b_nms  <- paste0(grp_col, b_vals)
      cv_vec[intersect(a_nms,names(cv_vec))] <-  1
      cv_vec[intersect(b_nms,names(cv_vec))] <- -1
      v    <- voom(dge,des,plot=FALSE); fit <- lmFit(v,des)
      fit2 <- contrasts.fit(fit,contrasts=cv_vec); fit2 <- eBayes(fit2)
      res  <- topTable(fit2,coef=1,number=Inf,sort.by="none")
      sig  <- !is.na(res$P.Value)&res$P.Value<=0.05&abs(res$logFC)>=0.5
      data.frame(peak_id=rownames(res),chr=as.character(seqnames(p_gr)),
        start=start(p_gr)-1,end=end(p_gr),
        log2FC=res$logFC,pvalue=res$P.Value,padj=res$adj.P.Val,
        AveExpr=res$AveExpr,baseMean=rowMeans(cnt_s),
        is_sig=sig,contrast=cmp$label,mode=mode_lbl,stringsAsFactors=FALSE)
    }

    res_n <- run_s8(counts_n,    peaks_n,        "NARROW",            iA,   iB,   ruv_n)
    res_b <- run_s8(cached$counts,cached$peaks,  paste0("EXTENDED",ebp), iA_b, iB_b, ruv_b)
    cat("  Narrow sig:",sum(res_n$is_sig)," Extended sig:",sum(res_b$is_sig),"\n")

    # [L5] Method agreement (limma narrow + limma extended for singleton contrasts)
    all_sig_ids <- unique(c(res_n$peak_id[res_n$is_sig], res_b$peak_id[res_b$is_sig]))
    in_n <- all_sig_ids %in% res_n$peak_id[res_n$is_sig]
    in_b <- all_sig_ids %in% res_b$peak_id[res_b$is_sig]
    ma_map <- setNames(ifelse(in_n&in_b,1.5,1.0), all_sig_ids)

    all_gr <- c(
      .df_to_gr_v3(nc$A_specific,"condition_specific"),
      .df_to_gr_v3(nc$B_specific,"condition_specific"),
      .df_to_gr_v3(bc$A_specific,"condition_specific"),
      .df_to_gr_v3(bc$B_specific,"condition_specific"),
      .df_to_gr_v3(res_n[res_n$is_sig,],"differential",ma_map),
      .df_to_gr_v3(res_b[res_b$is_sig,],"differential",ma_map))

    if (length(all_gr)>1) all_gr <- .dedup_v3(all_gr)
    mcols(all_gr)$direction <- ifelse(mcols(all_gr)$log2FC>0,"gained-open","gained-close")
    dar_list[[cmp$label]] <- all_gr

    # [L7] Tier 3 diagnostic
    .diag_tier3(gsub("[^A-Za-z0-9_]","_",cmp$label),
                all_gr, peaks_n, gpA, gpB, ebp, d_dir)

    safe_tag <- gsub("[^A-Za-z0-9_]","_",cmp$label)
    .write_dar_v3(all_gr, safe_tag, out_dir)
    cat("  Total DARs (after dedup):",length(all_gr),"\n\n")
  }

  if (length(dar_list)>0) {
    sumdf <- data.frame(
      Contrast    = names(dar_list),
      Total_DARs  = sapply(dar_list,length),
      Tier1_n     = sapply(dar_list,function(g) sum(mcols(g)$tier=="condition_specific")),
      Tier2_n     = sapply(dar_list,function(g)
                      sum(mcols(g)$tier=="differential"&grepl("^NARROW",mcols(g)$mode))),
      Tier3_n     = sapply(dar_list,function(g)
                      sum(mcols(g)$tier=="differential"&grepl("^EXTENDED",mcols(g)$mode))),
      Consensus_2 = sapply(dar_list,function(g)
                      sum(!is.na(mcols(g)$m_agree)&mcols(g)$m_agree>=1.5)),
      Gained      = sapply(dar_list,function(g) sum(mcols(g)$direction=="gained-open")),
      Lost        = sapply(dar_list,function(g) sum(mcols(g)$direction=="gained-close"))
    )
    cat("\n=== TEMPORAL DAR SUMMARY v3 ===\n"); print(sumdf); cat("\n")
    write.table(sumdf,file.path(out_dir,"temporal_DAR_v3_summary.txt"),
                sep="\t",quote=FALSE,row.names=FALSE)

    # Tier 3 diagnostic summary
    dfs <- list.files(d_dir,pattern="_Tier3_diagnostic.txt",full.names=TRUE)
    if (length(dfs)>0) {
      ds <- do.call(rbind,lapply(dfs,function(f){
        d<-read.table(f,header=TRUE,sep="\t"); t<-table(d$classification)
        data.frame(contrast=gsub("_Tier3_diagnostic.txt","",basename(f)),
          total=nrow(d),
          tight_boundary  =sum(t["tight_boundary"],  na.rm=TRUE),
          adjacent_merged =sum(t["adjacent_merged"], na.rm=TRUE),
          diffuse         =sum(t["diffuse"],         na.rm=TRUE))
      }))
      cat("Tier3 diagnostic summary:\n"); print(ds)
      write.table(ds,file.path(d_dir,"Tier3_summary.txt"),sep="\t",quote=FALSE,row.names=FALSE)
    }
    cat("Tier 1 outputs →",t1_dir,"\n")
    cat("Diagnostics    →",d_dir, "\n\n")
  }
  invisible(dar_list)
}

# Use extended universe if built, otherwise nc14b universe.
# The narrow counts loading block handles counting any new peaks automatically.
peaks_for_dar <- if (exists("extended_universe_bed") &&
                      !is.null(extended_universe_bed) &&
                      file.exists(extended_universe_bed)) {
  cat("Using extended temporal universe for DAR analysis.\n")
  extended_universe_bed
} else {
  existing_peaks_bed
}

temporal_dar_results <- run_temporal_dar(
  final_coldata = final_temporal_coldata,
  peaks_bed     = peaks_for_dar,
  contrasts_df  = new_contrasts
)

cat("=== NEXT STEPS ===\n")
cat("1. Confirm BAM paths in Section 1 new_candidate_bams\n")
cat("2. Run Section 3 temporal QC UMAP\n")
cat("3. Review nearest-neighbor table and confirm candidates\n")
cat("4. Run Section 4 staging test for BOT_late_indvA/B\n")
cat("5. Set staging_reclassify_indvA <- TRUE/FALSE in Section 5\n")
cat("6. Re-run build_final_temporal_coldata() to get final coldata\n")
cat("7. Re-run script — Section 8 DAR analysis executes automatically\n")
cat("8. Update 01_merge_replicate_peaks_to_group_beds.sh and generate_CPM_bigwigs.sh\n",
    "   to add nc14d, late, gastr groups\n")
