################################################################################
# VENT TEMPORAL DAR ANALYSIS — v4
#
# Changes in v3 (5 structural improvements — preserved in v4):
#   [L1]-[L5]+[L7] as documented below.
#
# New in v4 (performance — requires Precount_AllCounts_Matrices.r run first):
#   - All featureCounts calls replaced by cache loads from counts_cache/
#   - FC_THREADS read from SLURM_CPUS_PER_TASK (used in fallback path only)
#   - grp_ext_cache pre-loaded; on-the-fly extended Tier 1 counting removed
#     from contrast loop — eliminates 12 redundant featureCounts calls
#   - Narrow union counts always loaded from disk
#
# Previous enhancements (v2):
#   Tier 1 named outputs | baseMean/AveExpr | per-contrast extension_bp |
#   priority-scored deduplication
#
# New in v3 (5 structural improvements):
#
#  [L1] GROUP-SPECIFIC PEAKS FOR TIER 1
#       Tier 1 uses GroupPeaks_BOTv.bed and GroupPeaks_BOTCv.bed (produced by
#       Build_FullUniverse_AllGeno.sh) as counting substrates, not the full
#       union. For temporal contrasts within a genotype (BOTv_temporal,
#       BOTCv_temporal) the same genotype's group peaks are used for both
#       timepoints — the timepoint difference is captured in which BAM
#       columns are assigned to grpA vs grpB.
#
#  [L2] CPM-SCALED TIER 1 THRESHOLDS
#       min_reads and max_reads thresholds are computed from CPM × median
#       library size, scaled to the replicate count in the absent condition.
#       All groups here are n=2, so thresholds are symmetric, but the CPM
#       anchoring means they track actual sequencing depth.
#
#  [L3] CONTRAST-SPECIFIC RUV CONTROL SELECTION
#       RUVg control peaks selected by pairwise edgeR between the two
#       specific groups in each contrast, not globally. For the temporal
#       contrasts this is particularly important: nc14b vs nc14late peaks
#       stable within BOTv may differ from those stable within BOTCv.
#
#  [L4] BOUNDARY-BASED PEAK EXTENSION
#       Extended peaks grow from start/end boundaries, not the peak center.
#       Nearby peaks within ext_bp of each other are merged before counting.
#
#  [L5] CONTINUOUS CONFIDENCE SCORE
#       conf_score = method_agreement × tier_score × −log10(p).
#       method_agreement: 1.0 (one method), 1.5 (limma narrow + limma ext
#       agree), 2.0 (all four calls agree). Propagates into derived temporal
#       category assignments via scored_consistent_overlaps().
#
#  [L7] TIER 3 DIAGNOSTIC
#       Extended-only DARs classified as tight_boundary / adjacent_merged /
#       diffuse per contrast. Written to Diagnostics/ subdirectory.
#
# Sample layout:
#   BOTv  nc14b:   Nc14b_rep3 + Nc14b_rep1  (r=0.9691)
#   BOTv  nc14late: FoxL1-High_3_nc14d_IR (28361, rep1 NEW) + Nc14late_rep2 (rep2)
#                   NOTE: FoxL1-High_4_nc14d retired — lower intra-group r in QC
#   BOTv  gastr:    FoxL1-High_2_gastr_IR (28362) + FoxL1-High_3_gastr_IR (28363) [NEW GROUP]
#   BOTCv nc14b:   Nc14b_rep4 + Nc14b_rep5  (r=0.9683)
#   BOTCv nc14late: Nc14late_rep1 + Nc14late_rep2 (r=0.9371)
#   n=2 per group — DESeq2 + limma-voom both run
################################################################################

suppressPackageStartupMessages({
  library(DESeq2); library(limma); library(edgeR)
  library(RUVSeq); library(GenomicRanges); library(rtracklayer)
  library(ggplot2); library(dplyr); library(tidyr)
})

# Thread count — reads SLURM allocation automatically; falls back to 4.
FC_THREADS <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))

# Cache directory written by Precount_AllCounts_Matrices.r
CACHE_DIR <- "./counts_cache/Vent_temporal"

base_dir   <- "../Generate_fresh_counts"
output_dir <- file.path(base_dir, "Output/Split2_Vent_temporal_v3")
tier1_dir  <- file.path(output_dir, "Tier1_ConditionSpecific")
diag_dir   <- file.path(output_dir, "Diagnostics")
for (d in c(output_dir, tier1_dir, diag_dir))
  dir.create(d, recursive=TRUE, showWarnings=FALSE)

peak_dir           <- "./ATAC_NarrowPeaks"
full_universe_bed  <- file.path(peak_dir, "FullUniverse_nc14b_AllGeno_union_peaks.bed")
group_peak_beds    <- list(
  BOTv  = file.path(peak_dir, "GroupPeaks_BOTv.bed"),
  BOTCv = file.path(peak_dir, "GroupPeaks_BOTCv.bed")
)

bamfiles <- c(
  # BOTv nc14b — confirmed pair (r=0.9691)
  "./Control_Bams/Fully_flattened_potential/FoxL1-High_Dm_ATAC_Nc14b_rep3_noq_rmdup.noChrM.bam",
  "./Control_Bams/Fully_flattened_potential/FoxL1-High_Dm_ATAC_Nc14b_07_noq_rmdup.noChrM.bam",
  # BOTv nc14late — rep1 REPLACED: FoxL1-High_3_nc14d_IR (28361) over FoxL1-High_4_nc14d
  "./Control_Bams/Fully_flattened_potential/FoxL1-High_3_nc14d_IR_noq_rmdup.noChrM.bam",
  "./Control_Bams/Fully_flattened_potential/FoxL1-High_Dm_ATAC_Nc14late_rep2_noq_rmdup.noChrM.bam",
  # BOTv gastr — NEW confirmed pair (28362 + 28363)
  "./Control_Bams/Fully_flattened_potential/FoxL1-High_2_gastr_IR_noq_rmdup.noChrM.bam",
  "./Control_Bams/Fully_flattened_potential/FoxL1-High_3_gastr_IR_noq_rmdup.noChrM.bam",
  # BOTCv nc14b — confirmed pair (r=0.9683)
  # original rep1 EXCLUDED (non-ventralized) | original rep2 EXCLUDED (UMAP outlier)
  "./Control_Bams/Fully_flattened_potential/HLH54F-High_Dm_ATAC_Nc14b_rep4_noq_rmdup.noChrM.bam",
  "./Control_Bams/Fully_flattened_potential/HLH54F-High_Dm_ATAC_Nc14b_rep5_noq_rmdup.noChrM.bam",
  # BOTCv nc14late — confirmed pair (r=0.9371)
  "./Control_Bams/Fully_flattened_potential/HLH54F-High_Dm_ATAC_Nc14late_rep1_noq_rmdup.noChrM.bam",
  "./Control_Bams/Fully_flattened_potential/HLH54F-High_Dm_ATAC_Nc14late_rep2_noq_rmdup.noChrM.bam"
)
sample_names <- basename(bamfiles)

coldata <- data.frame(
  sample    = sample_names,
  genotype  = factor(c("BOTv","BOTv","BOTv","BOTv","BOTv","BOTv",
                        "BOTCv","BOTCv","BOTCv","BOTCv"),
                     levels=c("BOTv","BOTCv")),
  timepoint = factor(c("nc14b","nc14b","nc14late","nc14late","gastr","gastr",
                       "nc14b","nc14b","nc14late","nc14late"),
                     levels=c("nc14b","nc14late","gastr")),
  row.names = sample_names, stringsAsFactors=FALSE)
coldata$group <- factor(paste(coldata$genotype, coldata$timepoint, sep="_"),
  levels=c("BOTv_nc14b","BOTv_nc14late","BOTv_gastr",
           "BOTCv_nc14b","BOTCv_nc14late"))

log2fc_threshold    <- 0.5
pvalue_threshold    <- 0.05
min_cpm_present     <- 1.0    # [L2] CPM floor for "present" condition
max_cpm_absent      <- 0.25   # [L2] CPM ceiling per replicate for "absent"
extension_bp_global <- 200

################################################################################
# SHARED HELPER FUNCTIONS
################################################################################

import_bed <- function(f) {
  d  <- read.table(f, sep="\t", header=FALSE, stringsAsFactors=FALSE, comment.char="")
  gr <- GRanges(seqnames=d[,1], ranges=IRanges(start=d[,2]+1, end=d[,3]), strand="*")
  if (ncol(d)>=4) names(gr)<-d[,4] else names(gr)<-paste0("peak_",seq_along(gr))
  gr
}

create_saf <- function(peaks_gr, name_prefix="peak") {
  if (is.null(names(peaks_gr)))
    names(peaks_gr) <- paste0(name_prefix,"_",seq_along(peaks_gr))
  peaks_gr <- trim(peaks_gr)
  saf <- data.frame(GeneID=names(peaks_gr), Chr=as.character(seqnames(peaks_gr)),
                    Start=start(peaks_gr), End=end(peaks_gr), Strand=".",
                    stringsAsFactors=FALSE)
  bad <- apply(saf,1,function(x) any(is.na(x))); if (any(bad)) saf <- saf[!bad,]
  saf$GeneID <- make.unique(saf$GeneID)
  saf$Start  <- as.integer(saf$Start); saf$End <- as.integer(saf$End)
  inv <- saf$Start<1|saf$End<=saf$Start|is.na(saf$Start)|is.na(saf$End)
  if (any(inv)) saf <- saf[!inv,]
  bc <- grepl("\\s",saf$Chr)|saf$Chr==""|is.na(saf$Chr)
  if (any(bc)) saf <- saf[!bc,]
  saf[,c("GeneID","Chr","Start","End","Strand")]
}

count_reads_in_peaks <- function(peaks_gr, bamfiles, name_prefix="peak") {
  cat("  Counting reads in",length(peaks_gr),"peaks...\n")
  saf <- create_saf(peaks_gr, name_prefix)
  fc  <- featureCounts(files=bamfiles, annot.ext=saf, isPairedEnd=TRUE,
                       countMultiMappingReads=FALSE, primaryOnly=TRUE, nthreads=FC_THREADS)
  counts <- fc$counts; colnames(counts) <- basename(bamfiles)
  list(counts=counts, valid_peak_ids=saf$GeneID)
}

# [L4] Boundary-based extension — extends from start/end, merges overlaps
build_extended_peaks_boundary <- function(peaks_gr, ext_bp) {
  extended <- GRanges(
    seqnames = seqnames(peaks_gr),
    ranges   = IRanges(start=pmax(1L, start(peaks_gr)-ext_bp),
                       end  =end(peaks_gr)+ext_bp),
    strand   = "*")
  names(extended) <- names(peaks_gr)
  extended <- trim(extended)
  merged   <- reduce(extended, with.revmap=TRUE)
  orig_names   <- names(peaks_gr)
  merged_names <- sapply(mcols(merged)$revmap, function(idx) orig_names[idx[1]])
  names(merged) <- paste0("bext",ext_bp,"_",merged_names)
  merged
}

# [L2] CPM-scaled Tier 1 thresholds
compute_tier1_thresholds <- function(counts_all, grpA_idx, grpB_idx,
                                      min_cpm=min_cpm_present,
                                      max_cpm_per_rep=max_cpm_absent) {
  lib_sizes  <- colSums(counts_all[,c(grpA_idx,grpB_idx),drop=FALSE])
  median_lib <- median(lib_sizes)
  min_reads  <- max(1, round(min_cpm * median_lib / 1e6))
  n_absent   <- min(length(grpA_idx), length(grpB_idx))
  max_reads  <- max(1, round(max_cpm_per_rep * median_lib / 1e6 * n_absent))
  cat(sprintf("    Tier1 thresholds: min_reads=%d  max_reads_other=%d",
              min_reads, max_reads),
      sprintf("  (median_lib=%.0f)\n", median_lib))
  list(min_reads=min_reads, max_reads_other=max_reads)
}

################################################################################
# CACHE COLUMN ALIGNMENT
# Explicit map: precount superset key → BAM basename used by this script.
################################################################################

CACHE_KEY_MAP <- c(
  BOTv_rep1           = "FoxL1-High_Dm_ATAC_Nc14b_rep3_noq_rmdup.noChrM.bam",
  BOTv_rep2           = "FoxL1-High_Dm_ATAC_Nc14b_07_noq_rmdup.noChrM.bam",
  BOTv_nc14late_rep1  = "FoxL1-High_3_nc14d_IR_noq_rmdup.noChrM.bam",
  BOTv_nc14late_rep2  = "FoxL1-High_Dm_ATAC_Nc14late_rep2_noq_rmdup.noChrM.bam",
  BOTv_gastr_rep1     = "FoxL1-High_2_gastr_IR_noq_rmdup.noChrM.bam",
  BOTv_gastr_rep2     = "FoxL1-High_3_gastr_IR_noq_rmdup.noChrM.bam",
  BOTCv_rep1          = "HLH54F-High_Dm_ATAC_Nc14b_rep4_noq_rmdup.noChrM.bam",
  BOTCv_rep2          = "HLH54F-High_Dm_ATAC_Nc14b_rep5_noq_rmdup.noChrM.bam",
  BOTCv_nc14late_rep1 = "HLH54F-High_Dm_ATAC_Nc14late_rep1_noq_rmdup.noChrM.bam",
  BOTCv_nc14late_rep2 = "HLH54F-High_Dm_ATAC_Nc14late_rep2_noq_rmdup.noChrM.bam"
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

# [L1] Load group-specific peak counts
load_group_peak_counts <- function(group_peak_beds, bamfiles, sample_names) {
  cat("Loading group-specific peak counts [L1]...\n")
  grp_counts <- list(); grp_peaks <- list()
  for (grp in names(group_peak_beds)) {
    bed <- group_peak_beds[[grp]]
    if (!file.exists(bed)) {
      cat("  WARNING:", grp, "GroupPeaks BED not found —", bed, "\n")
      grp_counts[[grp]] <- NULL; grp_peaks[[grp]] <- NULL; next
    }
    gr <- import_bed(bed)
    cat(" ", grp, ":", length(gr), "peaks\n")
    cr <- count_reads_in_peaks(gr, bamfiles, paste0("grp_",grp))
    cnts <- cr$counts
    if (nrow(cnts)<length(gr)) gr <- gr[names(gr)%in%cr$valid_peak_ids]
    rownames(cnts) <- names(gr); colnames(cnts) <- sample_names
    grp_counts[[grp]] <- cnts; grp_peaks[[grp]] <- gr
  }
  list(counts=grp_counts, peaks=grp_peaks)
}

# [L1+L2] Condition-specific detection on group-specific peaks
identify_condition_specific_pair <- function(
    grp_peaks_A, grp_counts_A, grp_peaks_B, grp_counts_B,
    union_counts, grpA_idx, grpB_idx,
    grpA_label, grpB_label, min_reads, max_reads_other,
    mode_label, contrast_name) {

  # grpA_idx / grpB_idx are integer positions in the full coldata / counts_filt.
  # grp_counts_A/B only contain columns for that genotype's BAMs, so we must
  # match by sample name rather than position.
  sA <- colnames(union_counts)[grpA_idx]   # sample names for group A
  sB <- colnames(union_counts)[grpB_idx]   # sample names for group B

  # DIAGNOSTIC — remove after confirming Tier 1 works
  cat("    [T1 debug] grpA samples:", paste(sA, collapse=", "), "\n")
  cat("    [T1 debug] grpB samples:", paste(sB, collapse=", "), "\n")
  cat("    [T1 debug] grp_counts_A cols:", paste(colnames(grp_counts_A), collapse=", "), "\n")
  cat("    [T1 debug] grp_counts_B cols:", paste(colnames(grp_counts_B), collapse=", "), "\n")

  sA_in_grpA <- sA[sA %in% colnames(grp_counts_A)]
  sB_in_grpA <- sB[sB %in% colnames(grp_counts_A)]
  sB_in_grpB <- sB[sB %in% colnames(grp_counts_B)]
  sA_in_grpB <- sA[sA %in% colnames(grp_counts_B)]

  cat("    [T1 debug] sA_in_grpA:", paste(sA_in_grpA, collapse=", "), "\n")
  cat("    [T1 debug] sB_in_grpA:", paste(sB_in_grpA, collapse=", "), "\n")
  cat("    [T1 debug] sB_in_grpB:", paste(sB_in_grpB, collapse=", "), "\n")
  cat("    [T1 debug] sA_in_grpB:", paste(sA_in_grpB, collapse=", "), "\n")
  sB_in_grpB <- sB[sB %in% colnames(grp_counts_B)]
  sA_in_grpB <- sA[sA %in% colnames(grp_counts_B)]

  # If none of the contrast samples appear in the group counts, return empty
  if (length(sA_in_grpA) == 0 && length(sB_in_grpA) == 0) {
    empty_df <- data.frame(peak_id=character(),chr=character(),
      start=integer(),end=integer(),log2FC=numeric(),
      pvalue=numeric(),padj=numeric(),baseMean=numeric(),
      AveExpr=numeric(),method=character(),contrast=character(),
      mode=character(),stringsAsFactors=FALSE)
    return(list(A_specific=empty_df, B_specific=empty_df))
  }

  A_in_A <- if (length(sA_in_grpA)>0)
    rowSums(grp_counts_A[, sA_in_grpA, drop=FALSE]) else
    rep(0, nrow(grp_counts_A))
  A_in_B <- if (length(sB_in_grpA)>0)
    rowSums(grp_counts_A[, sB_in_grpA, drop=FALSE]) else
    rep(0, nrow(grp_counts_A))
  B_in_B <- if (length(sB_in_grpB)>0)
    rowSums(grp_counts_B[, sB_in_grpB, drop=FALSE]) else
    rep(0, nrow(grp_counts_B))
  B_in_A <- if (length(sA_in_grpB)>0)
    rowSums(grp_counts_B[, sA_in_grpB, drop=FALSE]) else
    rep(0, nrow(grp_counts_B))

  lib_size_est <- median(colSums(union_counts[,c(grpA_idx,grpB_idx),drop=FALSE])) + 1

  mk <- function(peaks_gr, present_counts, absent_counts, n_present, lfc, label) {
    idx <- present_counts >= min_reads & absent_counts <= max_reads_other
    # DIAGNOSTIC
    cat("      [T1 debug mk]", label,
        "— peaks:", length(present_counts),
        "| present>=", min_reads, ":", sum(present_counts >= min_reads),
        "| absent<=", max_reads_other, ":", sum(absent_counts <= max_reads_other),
        "| both:", sum(idx),
        "| present range:", paste(range(present_counts), collapse="-"),
        "| absent range:", paste(range(absent_counts), collapse="-"), "\n")
    if (sum(idx)==0) return(data.frame(
      peak_id=character(),chr=character(),start=integer(),end=integer(),
      log2FC=numeric(),pvalue=numeric(),padj=numeric(),
      baseMean=numeric(),AveExpr=numeric(),
      method=character(),contrast=character(),mode=character(),
      stringsAsFactors=FALSE))
    pm <- present_counts[idx] / max(n_present, 1)
    data.frame(
      peak_id  = names(peaks_gr)[idx],
      chr      = as.character(seqnames(peaks_gr[idx])),
      start    = start(peaks_gr[idx])-1, end=end(peaks_gr[idx]),
      log2FC   = lfc, pvalue=1e-10, padj=1e-10,
      baseMean = pm,
      AveExpr  = log2((pm/lib_size_est*1e6)+0.5),
      method   = paste0(mode_label,"_",label,"_specific_",contrast_name),
      contrast = contrast_name, mode=mode_label,
      stringsAsFactors=FALSE)
  }
  list(
    A_specific = mk(grp_peaks_A, A_in_A, A_in_B, length(grpA_idx), -5,
                    paste0(grpA_label,"_open_lost")),
    B_specific = mk(grp_peaks_B, B_in_B, B_in_A, length(grpB_idx),  5,
                    paste0(grpB_label,"_open_gained"))
  )
}

# [L3] Contrast-specific RUV control selection
select_ruv_controls_pairwise <- function(counts, grpA_idx, grpB_idx,
                                          grpA_label, grpB_label) {
  cat("    [L3] Pairwise RUV controls for",grpA_label,"vs",grpB_label,"...\n")
  pair_idx    <- c(grpA_idx, grpB_idx)
  counts_pair <- counts[,pair_idx,drop=FALSE]
  group_pair  <- factor(c(rep(grpA_label,length(grpA_idx)),
                          rep(grpB_label,length(grpB_idx))))
  pm   <- rowMeans(counts_pair); cv <- apply(counts_pair,1,sd)/(pm+0.5)
  exp_cut <- quantile(pm,0.4); cv_cut <- quantile(cv[pm>exp_cut],0.3)
  dge  <- DGEList(counts=counts_pair,group=group_pair)
  dge  <- calcNormFactors(dge)
  des  <- model.matrix(~group_pair)
  dge  <- estimateDisp(dge,des); fit <- glmQLFit(dge,des)
  res  <- topTags(glmQLFTest(fit,coef=2),n=Inf)$table
  non_de <- which(res$PValue>0.5&res$logCPM>1)
  stable <- which(pm>exp_cut&cv<cv_cut)
  int    <- intersect(stable,non_de)
  if (length(int)>=50)       { ctrl<-int;    m<-"intersection" }
  else if(length(stable)>=50){ ctrl<-stable; m<-"empirical_stable" }
  else if(length(non_de)>=50){ ctrl<-non_de; m<-"non_DE" }
  else { n<-min(200,max(100,round(nrow(counts_pair)*0.1)))
         ctrl<-order(cv[pm>quantile(pm,0.2)])[seq_len(n)]; m<-"top_stable_CV" }
  cat("      Controls:",length(ctrl),"(",m,")\n")
  rownames(counts)[ctrl]
}

# [L5] Method agreement map
compute_method_agreement <- function(res_n_sig, res_b_sig) {
  # For vent script: two limma calls (narrow + extended); agreement = 1.5 if both sig
  all_peaks <- unique(c(res_n_sig$peak_id, res_b_sig$peak_id))
  if (length(all_peaks)==0) return(setNames(numeric(0),character(0)))
  in_n <- all_peaks %in% res_n_sig$peak_id
  in_b <- all_peaks %in% res_b_sig$peak_id
  agreement <- ifelse(in_n & in_b, 1.5, 1.0)
  setNames(agreement, all_peaks)
}

# [L5] df_to_gr with conf_score
df_to_gr <- function(df, tier="differential", method_agreement_map=NULL) {
  if (nrow(df)==0) return(GRanges())
  eps <- 1e-300
  tier_score_val <- switch(tier,
    "condition_specific"=3,
    "differential"=ifelse(grepl("^NARROW",df$mode[1]),2,1), 1)
  m_agree <- if (tier=="condition_specific") {
    rep(1.0,nrow(df))
  } else if (!is.null(method_agreement_map)&&length(method_agreement_map)>0) {
    ma <- method_agreement_map[df$peak_id]; ma[is.na(ma)]<-1.0; as.numeric(ma)
  } else rep(1.0,nrow(df))
  pval_use   <- ifelse(is.na(df$pvalue),1,df$pvalue)
  conf_score <- m_agree * tier_score_val * (-log10(pval_use+eps))
  GRanges(seqnames=df$chr, ranges=IRanges(start=df$start+1,end=df$end),
          log2FC=df$log2FC, pvalue=pval_use,
          padj    =if("padj"    %in%names(df)) df$padj    else NA_real_,
          baseMean=if("baseMean"%in%names(df)) df$baseMean else NA_real_,
          AveExpr =if("AveExpr" %in%names(df)) df$AveExpr  else NA_real_,
          method  =if("method"  %in%names(df)) df$method   else tier,
          contrast=if("contrast"%in%names(df)) df$contrast else NA_character_,
          tier=tier, tier_score=tier_score_val, m_agree=m_agree,
          conf_score=conf_score,
          mode=if("mode"%in%names(df)) df$mode else tier)
}

# Priority deduplication using conf_score
priority_deduplicate <- function(gr) {
  if (length(gr)<=1) return(gr)
  cat("    Deduplicating",length(gr),"entries...\n")
  mcols(gr)$sign_dir <- sign(mcols(gr)$log2FC)
  dedup_one <- function(g) {
    if (length(g)<=1) return(g)
    ovl    <- findOverlaps(g,g,minoverlap=1L)
    ovl_df <- as.data.frame(ovl)
    ovl_df <- ovl_df[ovl_df$queryHits!=ovl_df$subjectHits,]
    if (nrow(ovl_df)==0) return(g)
    parent <- seq_len(length(g))
    find_root <- function(x) {
      while(parent[x]!=x){parent[x]<<-parent[parent[x]];x<-parent[x]}; x
    }
    for (i in seq_len(nrow(ovl_df))) {
      a<-find_root(ovl_df$queryHits[i]); b<-find_root(ovl_df$subjectHits[i])
      if(a!=b) parent[a]<-b
    }
    clusters <- sapply(seq_len(length(g)),find_root)
    scores   <- mcols(g)$conf_score
    keep <- sapply(unique(clusters),function(cl){m<-which(clusters==cl); m[which.max(scores[m])]})
    g[sort(keep)]
  }
  result <- c(dedup_one(gr[mcols(gr)$sign_dir>0]), dedup_one(gr[mcols(gr)$sign_dir<0]))
  cat("    After deduplication:",length(result),"\n")
  result
}

# [L5] Scored consistent overlaps for derived temporal categories
scored_consistent_overlaps <- function(gr1, gr2) {
  if (length(gr1)==0||length(gr2)==0) return(GRanges())
  if (!"conf_score"%in%names(mcols(gr1)))
    mcols(gr1)$conf_score <- mcols(gr1)$tier_score*(-log10(mcols(gr1)$pvalue+1e-300))
  if (!"conf_score"%in%names(mcols(gr2)))
    mcols(gr2)$conf_score <- mcols(gr2)$tier_score*(-log10(mcols(gr2)$pvalue+1e-300))
  proc <- function(g1s,g2s) {
    if (length(g1s)==0||length(g2s)==0) return(GRanges())
    ovl <- findOverlaps(g1s,g2s,minoverlap=1L)
    if (length(ovl)==0) return(GRanges())
    q_idx<-queryHits(ovl); s_idx<-subjectHits(ovl)
    cross<-tapply(mcols(g2s)$conf_score[s_idx],q_idx,max)
    hit  <- as.integer(names(cross)); res <- g1s[hit]
    mcols(res)$cross_conf_score <- as.numeric(cross)
    mcols(res)$conf_score <- sqrt(mcols(res)$conf_score * mcols(res)$cross_conf_score)
    res
  }
  c(proc(gr1[gr1$log2FC>0],gr2[gr2$log2FC>0]),
    proc(gr1[gr1$log2FC<0],gr2[gr2$log2FC<0]))
}

# [L3+L5] Combined limma + DESeq2 for n=2 per group contrasts
run_contrast_both_methods <- function(counts_sub, coldata_sub, peaks_gr,
                                       mode_label, contrast_name, grpA, grpB,
                                       ruv_ctrl_ids) {
  cat("  Running", contrast_name, "—", mode_label, "...\n")
  counts_sub <- round(counts_sub)
  # RUVg with pairwise controls
  valid_ctrl <- intersect(ruv_ctrl_ids, rownames(counts_sub))
  cat("    RUVg:", length(valid_ctrl), "controls\n")
  use_ruv <- length(valid_ctrl) >= 10
  if (use_ruv) {
    set     <- newSeqExpressionSet(counts=as.matrix(counts_sub),
                 phenoData=data.frame(condition=coldata_sub$group,
                                      row.names=colnames(counts_sub)))
    set_ruv <- RUVg(set, valid_ctrl, k=1)
    ruv_f   <- pData(set_ruv)[,grep("^W_",colnames(pData(set_ruv))),drop=FALSE]
    coldata_ruv <- cbind(coldata_sub, ruv_f)
    counts_use  <- counts(set_ruv)
    df_formula  <- as.formula(paste("~ 0 + group +",
                     paste(paste0("W_",seq_len(ncol(ruv_f))),collapse="+")))
  } else {
    coldata_ruv <- coldata_sub; counts_use <- counts_sub
    df_formula  <- ~ 0 + group
  }

  # DESeq2
  cat("    DESeq2...\n")
  dds    <- DESeqDataSetFromMatrix(countData=counts_use,
                                   colData=coldata_ruv, design=df_formula)
  dds    <- DESeq(dds, quiet=TRUE)
  grpA_f <- paste0("group",grpA); grpB_f <- paste0("group",grpB)
  res_d2 <- results(dds, contrast=c("group",grpA,grpB))
  d2_df  <- data.frame(
    peak_id=rownames(res_d2), chr=as.character(seqnames(peaks_gr)),
    start=start(peaks_gr)-1, end=end(peaks_gr),
    log2FC=res_d2$log2FoldChange, pvalue=res_d2$pvalue, padj=res_d2$padj,
    avg_expr=res_d2$baseMean,
    method=paste0(mode_label,"_DESeq2_",contrast_name),
    contrast=contrast_name, mode=mode_label, stringsAsFactors=FALSE)

  # limma-voom
  cat("    limma-voom...\n")
  dge  <- DGEList(counts=counts_use, group=coldata_sub$group)
  dge  <- calcNormFactors(dge, method="TMM")
  des  <- model.matrix(df_formula, data=coldata_ruv)
  cv   <- setNames(rep(0,ncol(des)),colnames(des))
  if (grpA_f%in%names(cv)) cv[grpA_f]<- 1
  if (grpB_f%in%names(cv)) cv[grpB_f]<--1
  v    <- voom(dge,des,plot=FALSE); fit <- lmFit(v,des)
  fit2 <- contrasts.fit(fit,contrasts=cv); fit2 <- eBayes(fit2)
  res_lm <- topTable(fit2,coef=1,number=Inf,sort.by="none")
  lm_df  <- data.frame(
    peak_id=rownames(res_lm), chr=as.character(seqnames(peaks_gr)),
    start=start(peaks_gr)-1, end=end(peaks_gr),
    log2FC=res_lm$logFC, pvalue=res_lm$P.Value, padj=res_lm$adj.P.Val,
    avg_expr=res_lm$AveExpr,
    method=paste0(mode_label,"_limma_",contrast_name),
    contrast=contrast_name, mode=mode_label, stringsAsFactors=FALSE)

  list(deseq2=d2_df, limma=lm_df)
}

filter_sig <- function(df) df %>%
  filter(!is.na(pvalue), pvalue<=pvalue_threshold, abs(log2FC)>=log2fc_threshold)

# [L7] Tier 3 diagnostic
diagnose_tier3_dars <- function(contrast_name, all_gr, counts_narrow,
                                 peaks_narrow, peaks_broad,
                                 grp_peaks_A, grp_peaks_B, ext_bp, out_dir) {
  tier3_gr <- all_gr[mcols(all_gr)$tier=="differential" &
                     grepl("EXTENDED",mcols(all_gr)$mode)]
  cat("  [L7] Tier3 diagnostic:", length(tier3_gr), "extended-only DARs\n")
  if (length(tier3_gr)==0) return(invisible(NULL))
  covered  <- countOverlaps(tier3_gr, peaks_narrow, minoverlap=1L) > 0
  in_A     <- if (!is.null(grp_peaks_A))
    countOverlaps(tier3_gr, grp_peaks_A, minoverlap=1L)>0 else rep(FALSE,length(tier3_gr))
  in_B     <- if (!is.null(grp_peaks_B))
    countOverlaps(tier3_gr, grp_peaks_B, minoverlap=1L)>0 else rep(FALSE,length(tier3_gr))
  flanked  <- GRanges(seqnames=seqnames(tier3_gr),
    ranges=IRanges(start=pmax(1L,start(tier3_gr)-ext_bp),end=end(tier3_gr)+ext_bp),
    strand="*")
  n_adj    <- countOverlaps(flanked, peaks_narrow, minoverlap=1L)
  cls      <- dplyr::case_when(
    covered & (in_A|in_B) ~ "tight_boundary",
    n_adj >= 2             ~ "adjacent_merged",
    covered                ~ "tight_boundary",
    TRUE                   ~ "diffuse")
  diag_df <- data.frame(
    chr=as.character(seqnames(tier3_gr)), start=start(tier3_gr)-1, end=end(tier3_gr),
    log2FC=mcols(tier3_gr)$log2FC, pvalue=mcols(tier3_gr)$pvalue,
    conf_score=mcols(tier3_gr)$conf_score,
    direction=ifelse(mcols(tier3_gr)$log2FC>0,"gained-open","gained-close"),
    has_narrow_peak=covered, in_grpA_peaks=in_A, in_grpB_peaks=in_B,
    n_adjacent_narrow=n_adj, classification=cls, stringsAsFactors=FALSE)
  out_f <- file.path(out_dir,paste0(contrast_name,"_Tier3_diagnostic.txt"))
  write.table(diag_df,out_f,sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
  tbl <- table(cls)
  cat("    Classification:",paste(names(tbl),tbl,sep="=",collapse=" | "),"\n")
  invisible(diag_df)
}

# Output writers
write_dar <- function(gr, tag, dir) {
  if (length(gr)==0) { cat("  Skip",tag,"(0)\n"); return(invisible(NULL)) }
  mc <- mcols(gr)
  df <- data.frame(
    chr=as.character(seqnames(gr)), start=start(gr)-1, end=end(gr),
    name=paste0(tag,"_",seq_along(gr)),
    score=round(-log10(mc$pvalue+1e-300)*10), strand=".",
    log2FC=mc$log2FC, pvalue=mc$pvalue,
    padj      =if("padj"      %in%names(mc)) mc$padj       else NA_real_,
    baseMean  =if("baseMean"  %in%names(mc)) mc$baseMean   else NA_real_,
    AveExpr   =if("AveExpr"   %in%names(mc)) mc$AveExpr    else NA_real_,
    conf_score=if("conf_score"%in%names(mc)) mc$conf_score else NA_real_,
    m_agree   =if("m_agree"   %in%names(mc)) mc$m_agree    else NA_real_,
    direction =ifelse(mc$log2FC>0,"gained-open","gained-close"),
    tier      =if("tier"      %in%names(mc)) mc$tier       else NA_character_,
    tier_score=if("tier_score"%in%names(mc)) mc$tier_score else NA_real_,
    mode      =if("mode"      %in%names(mc)) mc$mode       else NA_character_,
    method    =if("method"    %in%names(mc)) mc$method     else NA_character_,
    contrast  =if("contrast"  %in%names(mc)) mc$contrast   else NA_character_,
    stringsAsFactors=FALSE)
  write.table(df[,1:6],file.path(dir,paste0(tag,"_DARs.bed")),
              sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
  write.table(df,file.path(dir,paste0(tag,"_DARs_annotated.txt")),
              sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
  cat("  Saved",nrow(df),"→",paste0(tag,"_DARs.bed\n"))
}

# Unfiltered (pre-significance) results export -- matches the exact format
# plot_peak_fate_v11.r's read_full_results() already expects and has been
# silently degrading without (see its A3b/A5 comments): one row per tested
# peak, EVERY peak regardless of significance, columns chr/start/end/
# log2FoldChange/pvalue/padj. Without this, plot_peak_fate_v11.r can never
# rescue a locus that's real but just short of the contrast's own padj bar
# (exactly what happened to WntD at nc14late: real, HCR-validated, consistent
# direction across reps, but padj too high under n=2 to pass filter_sig()).
#
# Stacks NARROW + EXTENDED{ebp}, DESeq2 + limma raw results together (up to
# 4 rows per genomic region, at potentially different coordinates since
# EXTENDED windows are merged/wider). This is deliberate, not sloppy: the
# downstream fallback in plot_peak_fate_v11.r joins by genomic overlap, so
# stacking lets it see whichever mode/method had the most power for a given
# locus, rather than us pre-deciding which one "wins" at export time. We are
# NOT filtering on significance here -- that's the whole point of this file.
write_full_results <- function(df_list, tag, dir) {
  keep_cols <- c("chr","start","end","log2FC","pvalue","padj","method","mode")
  df_list <- lapply(df_list, function(d) {
    if (is.null(d) || nrow(d)==0) return(NULL)
    missing <- setdiff(keep_cols, names(d))
    for (mc in missing) d[[mc]] <- NA
    d[, keep_cols]
  })
  full <- do.call(rbind, df_list[!sapply(df_list, is.null)])
  if (is.null(full) || nrow(full)==0) {
    cat("  [write_full_results]", tag, "-- nothing to write (0 rows)\n")
    return(invisible(NULL))
  }
  # plot_peak_fate_v11.r's read_full_results() requires this exact column
  # name for the fold-change column: log2FoldChange (not log2FC).
  colnames(full)[colnames(full)=="log2FC"] <- "log2FoldChange"
  fn <- file.path(dir, paste0(tag, "_FULL_results.txt"))
  write.table(full, fn, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  cat("  [write_full_results]", tag, "->", nrow(full), "rows ->", basename(fn), "\n")
}

write_tier1_outputs <- function(cs_n, cs_e, contrast_name, t1dir) {
  dir.create(t1dir,recursive=TRUE,showWarnings=FALSE)
  write_one <- function(df, tag) {
    if (nrow(df)==0) return(invisible(NULL))
    df$direction<-ifelse(df$log2FC>0,"gained-open","gained-close")
    df$name<-paste0(tag,"_",seq_len(nrow(df))); df$strand<-"."
    write.table(df[,c("chr","start","end","name","pvalue","strand")],
      file.path(t1dir,paste0(tag,"_DARs.bed")),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
    write.table(df,file.path(t1dir,paste0(tag,"_DARs_annotated.txt")),
      sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
    cat("  T1 saved",nrow(df),"→",paste0(tag,"_DARs.bed\n"))
  }
  A <- rbind(cs_n$A_specific[cs_n$A_specific$log2FC<0,,drop=FALSE],
             cs_e$A_specific[cs_e$A_specific$log2FC<0,,drop=FALSE])
  A <- A[!duplicated(paste(A$chr,A$start,A$end)),]
  B <- rbind(cs_n$B_specific[cs_n$B_specific$log2FC>0,,drop=FALSE],
             cs_e$B_specific[cs_e$B_specific$log2FC>0,,drop=FALSE])
  B <- B[!duplicated(paste(B$chr,B$start,B$end)),]
  write_one(A,paste0("TIER1_",contrast_name,"_A_open_lost"))
  write_one(B,paste0("TIER1_",contrast_name,"_B_open_gained"))
  cat(sprintf("  TIER1 %s: A_lost=%d  B_gained=%d\n",contrast_name,nrow(A),nrow(B)))
}

################################################################################
# CONTRASTS — each may include extension_bp to override global default
################################################################################

comparisons_vent <- list(
  list(name="BOTv_vs_BOTCv",
       grpA="BOTv_nc14b",   grpB="BOTCv_nc14b",
       grpA_label="BOTv",   grpB_label="BOTCv",
       grpA_geno="BOTv",    grpB_geno="BOTCv"),

  list(name="BOTv_vs_BOTCv_nc14late",
       grpA="BOTv_nc14late", grpB="BOTCv_nc14late",
       grpA_label="BOTv",    grpB_label="BOTCv",
       grpA_geno="BOTv",     grpB_geno="BOTCv"),

  list(name="BOTv_temporal",
       grpA="BOTv_nc14late", grpB="BOTv_nc14b",
       grpA_label="BOTv_nc14late", grpB_label="BOTv_nc14b",
       grpA_geno="BOTv",           grpB_geno="BOTv",
       color_open="#1a6b1a", color_close="#2ca02c"),
  list(name="BOTv_gastr_vs_nc14b",
       grpA="BOTv_gastr",  grpB="BOTv_nc14b",
       grpA_label="BOTv_gastr",    grpB_label="BOTv_nc14b",
       grpA_geno="BOTv",           grpB_geno="BOTv",
       color_open="#0d3d0d", color_close="#2ca02c",
       extension_bp=250),
  list(name="BOTv_gastr_vs_nc14late",
       grpA="BOTv_gastr",  grpB="BOTv_nc14late",
       grpA_label="BOTv_gastr",    grpB_label="BOTv_nc14late",
       grpA_geno="BOTv",           grpB_geno="BOTv",
       color_open="#0d3d0d", color_close="#1a6b1a",
       extension_bp=250),

  list(name="BOTCv_temporal",
       grpA="BOTCv_nc14late", grpB="BOTCv_nc14b",
       grpA_label="BOTCv_nc14late", grpB_label="BOTCv_nc14b",
       grpA_geno="BOTCv",           grpB_geno="BOTCv",
       extension_bp=250)
)

################################################################################
# LOAD PEAKS AND COUNT READS
################################################################################

cat("\n=== VENT TEMPORAL DAR — v4 (L1-L5+L7, cache-accelerated) ===\n\n")
cat("Samples:\n"); print(coldata[,c("genotype","timepoint","group")]); cat("\n")

# ---------------------------------------------------------------------------
# LOAD UNION PEAK UNIVERSE
# ---------------------------------------------------------------------------
existing_peaks <- import_bed(full_universe_bed)
names(existing_peaks) <- paste0("peak_", seq_along(existing_peaks))

# ---------------------------------------------------------------------------
# NARROW COUNTS — load from cache or count on first run
# ---------------------------------------------------------------------------
narrow_cache_file <- file.path(CACHE_DIR, "narrow_counts.txt")
if (file.exists(narrow_cache_file)) {
  cat("Loading narrow counts from cache...\n")
  counts_all <- as.matrix(read.table(narrow_cache_file,
                                      header=TRUE, row.names=1, sep="\t"))
} else {
  cat("Cache not found — running featureCounts (first run).\n")
  cat("Run Precount_AllCounts_Matrices.r first to avoid this.\n")
  saf_main <- create_saf(existing_peaks, "peak")
  fc_main  <- featureCounts(files=bamfiles, annot.ext=saf_main, isPairedEnd=TRUE,
                            countMultiMappingReads=FALSE, primaryOnly=TRUE,
                            nthreads=FC_THREADS)
  counts_all <- fc_main$counts
  dir.create(CACHE_DIR, recursive=TRUE, showWarnings=FALSE)
  write.table(counts_all, narrow_cache_file,
              sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)
  cat("Saved for future runs.\n")
}
counts_all <- realign_cache_cols(counts_all, sample_names)
if (is.null(counts_all))
  stop("Could not align narrow counts cache columns to sample_names. ",
       "Check that the cache was built with the same BAMs as this script.")
rownames(counts_all) <- names(existing_peaks)[seq_len(nrow(counts_all))]
# Sync peaks object to rows present in counts
existing_peaks <- existing_peaks[rownames(counts_all)]

keep        <- rowSums(counts_all >= 5) >= 2
counts_filt <- counts_all[keep, ]
peaks_filt  <- existing_peaks[keep]
cat("Union peaks after filter:", length(peaks_filt), "\n\n")

# ---------------------------------------------------------------------------
# [L1] GROUP-SPECIFIC PEAK COUNTS — load from cache or count on first run
# ---------------------------------------------------------------------------
grp_cache_file <- file.path(CACHE_DIR, "grp_counts_cache.rds")
if (file.exists(grp_cache_file)) {
  cat("Loading group-specific peak counts from cache...\n")
  grp_cache <- readRDS(grp_cache_file)
  grp_cache <- realign_grp_cache(grp_cache, sample_names)
} else {
  cat("Cache not found — counting group peaks (first run).\n")
  grp_cache <- load_group_peak_counts(group_peak_beds, bamfiles, sample_names)
  dir.create(CACHE_DIR, recursive=TRUE, showWarnings=FALSE)
  saveRDS(grp_cache, grp_cache_file)
  cat("Saved:", grp_cache_file, "\n")
}
cat("\n")

# ---------------------------------------------------------------------------
# [L4] EXTENDED UNION PEAK COUNTS — load from cache or count on first run
# ---------------------------------------------------------------------------
cat("=== L4: Loading boundary-extended peak counts ===\n")
unique_ext_bps <- sort(unique(c(extension_bp_global,
  sapply(comparisons_vent, function(x) if (!is.null(x$extension_bp)) x$extension_bp else NA),
  NA))) %>% {.[!is.na(.)]}
ext_cache <- list()
for (ebp in unique_ext_bps) {
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
    if (nrow(cnts) < length(p_ext)) p_ext <- p_ext[names(p_ext) %in% cr$valid_peak_ids]
    rownames(cnts) <- names(p_ext)
    cnts <- realign_cache_cols(cnts, sample_names)
    kp <- rowSums(cnts >= 5) >= 2
    ext_cache[[key]] <- list(counts=cnts[kp, ], peaks=p_ext[kp])
    dir.create(CACHE_DIR, recursive=TRUE, showWarnings=FALSE)
    saveRDS(ext_cache[[key]], cache_file)
    cat("  Saved:", cache_file, "\n")
  }
}
cat("\n")

# ---------------------------------------------------------------------------
# [L1+L4] EXTENDED GROUP-SPECIFIC PEAK COUNTS — load from cache
# Eliminates 2 featureCounts calls per contrast (~12 calls saved across 6).
# ---------------------------------------------------------------------------
cat("=== L1+L4: Loading extended group-specific peak counts ===\n")
grp_ext_cache <- list()
for (ebp in unique_ext_bps) {
  key        <- as.character(ebp)
  cache_file <- file.path(CACHE_DIR, paste0("grp_ext", ebp, "_counts.rds"))
  if (file.exists(cache_file)) {
    cat("  Loading grp_ext", ebp, "bp from cache\n")
    grp_ext_cache[[key]] <- readRDS(cache_file)
    grp_ext_cache[[key]] <- realign_grp_cache(grp_ext_cache[[key]], sample_names)
  } else {
    cat("  Cache missing — counting extended group peaks (first run).\n")
    grp_counts_ext <- list(); grp_peaks_ext <- list()
    for (grp in names(group_peak_beds)) {
      bed <- group_peak_beds[[grp]]
      if (!file.exists(bed)) next
      gp     <- import_bed(bed)
      gp_ext <- build_extended_peaks_boundary(gp, ebp)
      cr     <- count_reads_in_peaks(gp_ext, bamfiles, paste0("grpext", ebp, "_", grp))
      cnts   <- cr$counts; colnames(cnts) <- sample_names
      gp_ext <- gp_ext[names(gp_ext) %in% cr$valid_peak_ids]
      rownames(cnts)        <- names(gp_ext)
      grp_counts_ext[[grp]] <- cnts
      grp_peaks_ext[[grp]]  <- gp_ext
    }
    grp_ext_cache[[key]] <- list(counts=grp_counts_ext, peaks=grp_peaks_ext)
    dir.create(CACHE_DIR, recursive=TRUE, showWarnings=FALSE)
    saveRDS(grp_ext_cache[[key]], cache_file)
    cat("  Saved:", cache_file, "\n")
  }
}
cat("\n")

idx <- function(grp) which(coldata$group == grp)

################################################################################
# RUN CONTRASTS
################################################################################

dar_list <- list()

for (cmp in comparisons_vent) {
  cat(strrep("-",60),"\n")
  cat("CONTRAST:",cmp$name,"\n\n")

  ebp    <- if (!is.null(cmp$extension_bp)) cmp$extension_bp else extension_bp_global
  cached <- ext_cache[[as.character(ebp)]]
  iA     <- idx(cmp$grpA); iB <- idx(cmp$grpB)
  cat("  extension_bp =",ebp,"(boundary-based, L4)\n")

  # [L2] CPM-scaled thresholds
  t1_thr <- compute_tier1_thresholds(counts_filt, iA, iB)

  # [L1] Group-specific Tier 1 peaks
  # For genotype contrasts: use each genotype's group peaks
  # For temporal contrasts within same genotype: use that genotype's group peaks for both
  gA <- cmp$grpA_geno; gB <- cmp$grpB_geno
  grp_A_peaks  <- grp_cache$peaks[[gA]]
  grp_A_counts <- grp_cache$counts[[gA]]
  grp_B_peaks  <- grp_cache$peaks[[gB]]
  grp_B_counts <- grp_cache$counts[[gB]]

  if (is.null(grp_A_counts) || is.null(grp_B_counts)) {
    cat("  WARNING: group peaks missing — skipping Tier 1 for this contrast\n")
    nc <- list(A_specific=data.frame(),B_specific=data.frame())
    bc <- list(A_specific=data.frame(),B_specific=data.frame())
  } else {
    nc <- identify_condition_specific_pair(
      grp_A_peaks, grp_A_counts, grp_B_peaks, grp_B_counts,
      counts_filt, iA, iB, cmp$grpA_label, cmp$grpB_label,
      t1_thr$min_reads, t1_thr$max_reads_other, "NARROW", cmp$name)

    # [L1+L4] Extended group-specific peaks — from grp_ext_cache (v4)
    ebp_key    <- as.character(ebp)
    grp_A_ext  <- grp_ext_cache[[ebp_key]]$peaks[[gA]]
    cnts_Ae    <- grp_ext_cache[[ebp_key]]$counts[[gA]]
    grp_B_ext  <- grp_ext_cache[[ebp_key]]$peaks[[gB]]
    cnts_Be    <- grp_ext_cache[[ebp_key]]$counts[[gB]]

    if (is.null(cnts_Ae) || is.null(cnts_Be)) {
      cat("  WARNING: grp_ext cache missing for", gA, "/", gB,
          "— falling back to featureCounts\n")
      grp_A_ext <- build_extended_peaks_boundary(grp_A_peaks, ebp)
      grp_B_ext <- build_extended_peaks_boundary(grp_B_peaks, ebp)
      cr_Ae <- count_reads_in_peaks(grp_A_ext, bamfiles, paste0("grpext_", gA))
      cr_Be <- count_reads_in_peaks(grp_B_ext, bamfiles, paste0("grpext_", gB))
      cnts_Ae <- cr_Ae$counts; colnames(cnts_Ae) <- sample_names
      cnts_Be <- cr_Be$counts; colnames(cnts_Be) <- sample_names
      grp_A_ext <- grp_A_ext[names(grp_A_ext) %in% cr_Ae$valid_peak_ids]
      grp_B_ext <- grp_B_ext[names(grp_B_ext) %in% cr_Be$valid_peak_ids]
    }

    bc <- identify_condition_specific_pair(
      grp_A_ext, cnts_Ae, grp_B_ext, cnts_Be,
      counts_filt, iA, iB, cmp$grpA_label, cmp$grpB_label,
      t1_thr$min_reads, t1_thr$max_reads_other,
      paste0("EXTENDED",ebp), cmp$name)
  }
  cat("  Tier1 narrow:", nrow(nc$A_specific),"+",nrow(nc$B_specific),"\n")
  cat("  Tier1 ext",ebp,"bp:", nrow(bc$A_specific),"+",nrow(bc$B_specific),"\n")
  write_tier1_outputs(nc, bc, cmp$name, tier1_dir)

  # [L3] Pairwise RUV controls
  ruv_ctrl_n <- select_ruv_controls_pairwise(counts_filt, iA, iB,
                                              cmp$grpA_label, cmp$grpB_label)
  ruv_ctrl_b <- select_ruv_controls_pairwise(cached$counts, iA, iB,
                                              cmp$grpA_label, cmp$grpB_label)

  # Tier 2/3 — both DESeq2 + limma with contrast-specific RUV (n=2 supports both)
  cdat_sub <- coldata[c(iA,iB),]
  nd <- run_contrast_both_methods(counts_filt[,c(iA,iB)], cdat_sub, peaks_filt,
          "NARROW", cmp$name, cmp$grpA, cmp$grpB, ruv_ctrl_n)
  bd <- run_contrast_both_methods(cached$counts[,c(iA,iB)], cdat_sub, cached$peaks,
          paste0("EXTENDED",ebp), cmp$name, cmp$grpA, cmp$grpB, ruv_ctrl_b)

  # Export UNFILTERED results before any significance filtering discards them
  # -- this is the file plot_peak_fate_v11.r has been missing for this exact
  # contrast (BOTv_vs_BOTCv_nc14late) the whole time.
  write_full_results(list(nd$deseq2, nd$limma, bd$deseq2, bd$limma),
                      cmp$name, output_dir)

  nd2 <- filter_sig(nd$deseq2); nlm <- filter_sig(nd$limma)
  bd2 <- filter_sig(bd$deseq2); blm <- filter_sig(bd$limma)
  cat("  Narrow  DESeq2:",nrow(nd2)," limma:",nrow(nlm),"\n")
  cat("  Extended DESeq2:",nrow(bd2)," limma:",nrow(blm),"\n")

  # [L5] Method agreement (here: 2 limma calls + 2 DESeq2 calls = up to 4)
  sig_n <- rbind(nd2,nlm); sig_n <- sig_n[!duplicated(sig_n$peak_id),]
  sig_b <- rbind(bd2,blm); sig_b <- sig_b[!duplicated(sig_b$peak_id),]
  all_sig <- unique(c(nd2$peak_id,nlm$peak_id,bd2$peak_id,blm$peak_id))
  in_nd2  <- all_sig%in%nd2$peak_id; in_nlm <- all_sig%in%nlm$peak_id
  in_bd2  <- all_sig%in%bd2$peak_id; in_blm <- all_sig%in%blm$peak_id
  n_meth  <- in_nd2+in_nlm+in_bd2+in_blm
  ma_map  <- setNames(ifelse(n_meth==4,2.0,ifelse(n_meth>=2,1.5,1.0)), all_sig)

  # Assemble
  all_gr <- c(
    df_to_gr(nc$A_specific,"condition_specific"),
    df_to_gr(nc$B_specific,"condition_specific"),
    df_to_gr(bc$A_specific,"condition_specific"),
    df_to_gr(bc$B_specific,"condition_specific"),
    df_to_gr(nd2,"differential",ma_map),
    df_to_gr(nlm,"differential",ma_map),
    df_to_gr(bd2,"differential",ma_map),
    df_to_gr(blm,"differential",ma_map))

  if (length(all_gr)>1) all_gr <- priority_deduplicate(all_gr)
  mcols(all_gr)$direction <- ifelse(mcols(all_gr)$log2FC>0,"gained-open","gained-close")
  dar_list[[cmp$name]] <- all_gr

  # [L7] Tier 3 diagnostic
  diagnose_tier3_dars(cmp$name, all_gr, counts_filt, peaks_filt,
                      cached$peaks, grp_cache$peaks[[gA]], grp_cache$peaks[[gB]],
                      ebp, diag_dir)

  cat("  Total DARs (after dedup):",length(all_gr),"\n\n")
}

################################################################################
# DERIVED TEMPORAL CATEGORIES — scored overlaps (L5)
################################################################################

cat("\nDeriving temporal categories (scored overlaps, L5)...\n")
gr_AvsC_b <- dar_list[["BOTv_vs_BOTCv"]]
gr_AvsC_l <- dar_list[["BOTv_vs_BOTCv_nc14late"]]
gr_A_temp <- dar_list[["BOTv_temporal"]]
gr_C_temp <- dar_list[["BOTCv_temporal"]]

stable_AvsC     <- scored_consistent_overlaps(gr_AvsC_b, gr_AvsC_l)
emerging_late   <- gr_AvsC_l[countOverlaps(gr_AvsC_l, gr_AvsC_b)==0]
resolving_early <- gr_AvsC_b[countOverlaps(gr_AvsC_b, gr_AvsC_l)==0]
temporal_div    <- if (length(gr_A_temp)>0&&length(gr_C_temp)>0) {
  c(subsetByOverlaps(gr_A_temp[gr_A_temp$log2FC>0],gr_C_temp[gr_C_temp$log2FC<0],minoverlap=1L),
    subsetByOverlaps(gr_A_temp[gr_A_temp$log2FC<0],gr_C_temp[gr_C_temp$log2FC>0],minoverlap=1L))
} else GRanges()

# Scored temporal consensus: stable signal confirmed in both timepoints with
# high confidence in both
temporal_consensus <- scored_consistent_overlaps(gr_AvsC_b, gr_AvsC_l)

cat(sprintf("  Stable BOTv/BOTCv (scored, both TPs): %d\n", length(stable_AvsC)))
cat(sprintf("  Emerging (late only)                 : %d\n", length(emerging_late)))
cat(sprintf("  Resolving (early only)               : %d\n", length(resolving_early)))
cat(sprintf("  Temporal divergent                   : %d\n", length(temporal_div)))
cat(sprintf("  High-confidence temporal consensus   : %d\n", length(temporal_consensus)))

################################################################################
# SAVE
################################################################################

cat("\nSaving outputs...\n")
for (nm in names(dar_list)) write_dar(dar_list[[nm]], nm, output_dir)
write_dar(stable_AvsC,        "BOTv_BOTCv_stable_both_timepoints",       output_dir)
write_dar(temporal_consensus, "BOTv_BOTCv_stable_scored_consensus",      output_dir)
write_dar(emerging_late,      "BOTv_BOTCv_emerging_late_only",           output_dir)
write_dar(resolving_early,    "BOTv_BOTCv_resolving_early_only",         output_dir)
write_dar(temporal_div,       "temporal_divergent_BOTv_vs_BOTCv",        output_dir)

# Summary
summary_df <- data.frame(
  Contrast    = names(dar_list),
  N_DARs      = sapply(dar_list, length),
  Tier1_n     = sapply(dar_list, function(gr) sum(mcols(gr)$tier=="condition_specific")),
  Tier2_n     = sapply(dar_list, function(gr)
                  sum(mcols(gr)$tier=="differential"&grepl("^NARROW",mcols(gr)$mode))),
  Tier3_n     = sapply(dar_list, function(gr)
                  sum(mcols(gr)$tier=="differential"&grepl("^EXTENDED",mcols(gr)$mode))),
  Consensus4  = sapply(dar_list, function(gr)
                  sum(!is.na(mcols(gr)$m_agree)&mcols(gr)$m_agree>=2.0)),
  Gained      = sapply(dar_list, function(gr) sum(mcols(gr)$direction=="gained-open")),
  Lost        = sapply(dar_list, function(gr) sum(mcols(gr)$direction=="gained-close"))
)
cat("\nSummary:\n"); print(summary_df)
write.table(summary_df, file.path(output_dir,"Split2_vent_v3_summary.txt"),
            sep="\t",quote=FALSE,row.names=FALSE)

# Tier 3 diagnostic summary
diag_files <- list.files(diag_dir,pattern="_Tier3_diagnostic.txt",full.names=TRUE)
if (length(diag_files)>0) {
  diag_summ <- do.call(rbind, lapply(diag_files, function(f) {
    d <- read.table(f,header=TRUE,sep="\t"); tbl <- table(d$classification)
    data.frame(contrast=gsub("_Tier3_diagnostic.txt","",basename(f)),
               total=nrow(d),
               tight_boundary  =sum(tbl["tight_boundary"],  na.rm=TRUE),
               adjacent_merged =sum(tbl["adjacent_merged"], na.rm=TRUE),
               diffuse         =sum(tbl["diffuse"],         na.rm=TRUE))
  }))
  cat("\nTier3 diagnostic summary:\n"); print(diag_summ)
  write.table(diag_summ,file.path(diag_dir,"Tier3_summary.txt"),
              sep="\t",quote=FALSE,row.names=FALSE)
}

cat("\n=== VENT TEMPORAL v3 COMPLETE ===\n")
cat("Output   :", output_dir, "\n")
cat("Tier 1   :", tier1_dir, "\n")
cat("Diagnostics:", diag_dir, "\n\n")
cat("New output columns: conf_score | m_agree | cross_conf_score\n\n")
