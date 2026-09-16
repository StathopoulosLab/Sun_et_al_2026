################################################################################
# 04_mixed_effects_logistic_regression_state_selective_accessibility.r
#
# NEW IN v8: writes lmm_dar_summary.txt per stratum (n_total, n_logFC_pos/neg,
# and BOTv-/BOTCv-preferential counts for the two stems where the sign
# convention is known -- see NEGATE_LFC_STEMS). This answers "total DARs /
# how many preferential each direction / anything excluded before fitting"
# directly and reproducibly, rather than by manual count.
#
# NEW IN v7: adds an optional, genuine bootstrap stability check (resampling
# WITH replacement -- see BOOTSTRAP STABILITY CHECK section below) alongside
# the primary single-fit result, so you can report both a model-based
# coefficient/SE/p-value AND an empirical measure of how sensitive that
# coefficient is to which peaks happen to be in the dataset. This is
# deliberately kept SEPARATE from the significance call (no Stouffer
# combination, no p-value inflation) -- see the header comment on
# bootstrap_stability_check() for details. Set RUN_BOOTSTRAP_STABILITY <-
# FALSE to skip it and behave exactly like v6.
#
# SANITIZATION NOTES (relative to v4) — read before trusting old outputs:
#
#  * BUG: with LMM_N_SUB=8000 and all reported DAR sets < 8000 peaks,
#    sample(n, n, replace=FALSE) returned a permutation of the FULL dataset
#    every time. All N_RUNS "independent" fits were the same model fit to
#    the same data. This is fixed below: the code now detects whether real
#    subsampling occurs (n_peaks_total > LMM_N_SUB) and takes a different,
#    honest path when it does not (REAL_SUBSAMPLE = FALSE, see below).
#
#  * BUG: combine_runs() summed z-scores across the (near-)duplicate runs
#    and divided by sqrt(N_RUNS) (Stouffer's method), which assumes
#    independence. With duplicate/highly-overlapping runs this inflates
#    z and deflates p far below what a single honest fit would give.
#    BH correction afterward cannot repair this — it controls FDR across
#    TFs given the p-values handed to it, not the inflation itself.
#    Fixed: when REAL_SUBSAMPLE is FALSE, significance now comes directly
#    from a single lmerTest fit (Satterthwaite df) on the full stratum.
#    When REAL_SUBSAMPLE is TRUE, combine_runs_bootstrap() reports the
#    EMPIRICAL between-run SD (not inverse-variance pooling) and is
#    labeled explicitly as a stability/bootstrap summary, not a
#    replicated-experiment meta-analysis.
#
#  * BUG: error bars used a pooled inverse-variance SE that shrinks with
#    more (non-independent) runs, while the plot caption claimed the bars
#    showed "variance across subsampling runs." Fixed: bars now show what
#    they say they show — single-fit model SE in the no-real-subsampling
#    path, empirical between-run CI in the real-subsampling path.
#
#  * ADDED: per-stratum diagnostics saved to lmm_model_stats.txt —
#    n_singular, n_nonconverged, and a max-pairwise-|r| collinearity flag
#    across the TF design matrix (several partner_peaks entries are
#    different FDR thresholds of the same ChIP dataset and are expected
#    to be correlated).
#
# Changes over v3 (carried forward, unchanged):
#
#  9. GLMMLASSO MODELS: Two additional regularized models are now fit for
#     every comparison × LFC-threshold combination, in addition to the
#     existing plain LMM and binary GLMM:
#
#       (a) Continuous regularized GLMM-LASSO  (glmmLasso, gaussian family)
#           Response  : logFC (continuous)
#           Penalty   : L1 on all TF fixed effects; chr as random intercept
#           λ selection: BIC-guided decreasing path (start.lambda → 0.01),
#                        then 5-fold CV on held-out MSE to pick best λ
#           Output    : glmmlasso_continuous_coefs.txt
#                       glmmlasso_continuous_lfc<thresh>.pdf
#
#       (b) Binary regularized GLMM-LASSO  (glmmLasso, binomial family)
#           Response  : Opening (1) vs Closing (0)
#           Same λ selection strategy as (a)
#           Output    : glmmlasso_binary_coefs.txt
#                       glmmlasso_binary_lfc<thresh>.pdf
#
#     Plot style matches the reference plot:
#       - TFs ordered by estimate on y-axis
#       - Colored by significance (p < 0.05 after BH correction)
#       - Gray zone annotation for |effect| < 0.020
#       - Footer: "Method: regularized_glmm_lasso"
#
#     All existing LMM / binary GLMM outputs are preserved unchanged.
#
# Changes over v2 (carried from v3):
#
#  6. MULTI-THRESHOLD LFC FILTERING
#  7. BUG FIX: %||% null-safe accessor
#  8. BUG FIX: Early return(NULL) when n_total == 0
#
# Run from: Enhanced_pseudopink_RuntD7integration/
# Usage:  source("RERUN_LMM_ONLY_v4.r")
#         # To force rerun specific stems:
#         # FORCE_RERUN_KEYS <- c("BOTv_vs_BOTCv","BOTv_vs_BOTCv_late",...)
#         # source("RERUN_LMM_ONLY_v4.r")
################################################################################

options(stringsAsFactors = FALSE, expressions = 5e5)
if (.Device == "null device" || identical(.Device, "RStudioGD")) {
  grDevices::pdf(NULL)
}

suppressPackageStartupMessages({
  library(GenomicRanges); library(rtracklayer)
  library(lme4); library(lmerTest)
  library(ggplot2); library(dplyr); library(tidyr)
})

# glmmLasso is loaded conditionally — warn clearly if missing so the rest
# of the pipeline (plain LMM) still runs without it
GLMMLASSO_AVAILABLE <- requireNamespace("glmmLasso", quietly = TRUE)
if (GLMMLASSO_AVAILABLE) {
  suppressPackageStartupMessages(library(glmmLasso))
  cat("  glmmLasso package found — regularized models will be fit.\n")
} else {
  warning(paste0(
    "glmmLasso package not found.\n",
    "  Install with: install.packages('glmmLasso')\n",
    "  Plain LMM / binary GLMM outputs will still be produced.\n",
    "  Regularized GLMM-LASSO outputs will be SKIPPED."
  ))
}

cat("\n", strrep("=",70), "\n")
cat("RERUN_LMM_ONLY v3 — multi-threshold LFC filtering + bug fixes\n")
cat(strrep("=",70), "\n\n")

################################################################################
# SETTINGS
################################################################################

dar_base_dir      <- "./ATAC_DARS"
temp_base_dir     <- "./ATAC_DARS"

# ── OUTPUT DIRECTORY — NEVER reuse a path you've written results to before ──
# This is the #1 way to lose or silently corrupt a prior run: a script that
# reads its own output paths from a previous version (v4 -> v5 -> v6 all
# pointed at the same "_newcic" folder here) will happily overwrite files
# that used a different schema, leaving a directory with a mix of old- and
# new-format files that no downstream script can parse consistently.
#
# RUN_TAG makes every run's output land in its own uniquely-named directory
# by default — set it once per run (or leave the Sys.Date() default, which
# changes every day) rather than hardcoding a bare folder name. If you
# genuinely want to overwrite a specific prior run's output on purpose,
# set RUN_TAG <- "" and point split1_output/split2_output at that exact
# path yourself — that should be a deliberate, explicit choice, never the
# default.
RUN_TAG <- format(Sys.Date(), "%Y%m%d_bootstrapCI")   # e.g. "20260903_bootstrapCI" -- change/remove deliberately, not by accident

split1_output <- if (nzchar(RUN_TAG)) sprintf("./Split1_NonVent_newcic_%s", RUN_TAG) else "./Split1_NonVent_newcic"
split2_output <- if (nzchar(RUN_TAG)) sprintf("./Split2_Vent_newcic_%s", RUN_TAG)    else "./Split2_Vent_newcic"

cat(sprintf("\n[OUTPUT] Split 1 results -> %s\n", normalizePath(split1_output, mustWork = FALSE)))
cat(sprintf("[OUTPUT] Split 2 results -> %s\n", normalizePath(split2_output, mustWork = FALSE)))
if (dir.exists(split1_output) || dir.exists(split2_output)) {
  cat("[OUTPUT] WARNING: one of the above directories already exists.\n")
  cat("         Files inside it may be overwritten by this run. If that's not\n")
  cat("         intended, stop now (Ctrl-C) and change RUN_TAG above.\n\n")
} else {
  cat("[OUTPUT] Both directories are new -- no existing results at risk.\n\n")
}

# Directory where Step1_SevenGeno wrote its output — used as fallback when
# annotated files are not yet copied to ATAC_DARS/
# Adjust if your output directory has a different version suffix
SEVENGENO_OUT_DIR <- "./Generate_fresh_counts/Output/SevenGeno_nc14b_v3"
VENT_TEMP_OUT_DIR <- "./Generate_fresh_counts/Output/Split2_Vent_temporal_v3"

N_RUNS          <- 20      # independent subsampling runs to average (20 recommended)
LFC_THRESHOLDS  <- c(0.50, 0.35, 0.25)  # |log2FC| cutoffs — one sub-dir per threshold
LMM_N_SUB       <- 8000   # peaks per run (larger = more power per run)
PADJ_METHOD     <- "BH"
RUN_BINARY      <- TRUE   # also fit Opening/Closing binary model (glmer)
RUN_ALL_PEAKS   <- TRUE   # also run on full DAR set with no secondary |logFC| filter
                           # outputs go to lmm/all_peaks/ alongside lfc*/ directories
                           # justification: DARs are already significance-filtered at
                           # DESeq2/limma level — no secondary threshold needed

# ── glmmLasso settings ────────────────────────────────────────────────────────
# Lambda grid: linear seq matching original 3_analysis_core_functions.r exactly.
# BIC is used to select the best lambda (not cross-validation).
LASSO_LAMBDA_START  <- 1.0   # upper end of lambda grid (most regularized)
LASSO_LAMBDA_END    <- 0.01  # lower end of lambda grid (least regularized)
LASSO_N_LAMBDA      <- 20    # number of steps — matches original seq(0.01,1,length.out=20)
LASSO_N_SUB         <- 5000  # subsample size (glmmLasso fits N_LAMBDA models per comparison)
# Effect threshold for gray zone in plot (matches reference: |effect| < 0.020)
LASSO_EFFECT_THRESH <- 0.020

# Keys to force-rerun even if lmm_coefficients.txt exists
# (use this to regenerate after sign fix)
FORCE_RERUN_KEYS <- if (exists("FORCE_RERUN_KEYS")) FORCE_RERUN_KEYS else
  c("BOTv_vs_BOTCv", "BOTv_vs_BOTCv_late")

# Stems whose logFC must be negated (Step1 v4 BOTv-numerator convention).
# ONLY the two direct BOTv_vs_BOTCv contrasts: grpA=BOTv, grpB=BOTCv → raw positive = BOTv open.
# temporal_divergent is derived from BOTv_temporal + BOTCv_temporal (both already correct); do NOT negate.
NEGATE_LFC_STEMS <- c("BOTv_vs_BOTCv", "BOTv_vs_BOTCv_nc14late")

################################################################################
# LOAD FUNCTIONS AND DATA  (same as v1)
################################################################################

cat("[1/5] Loading pipeline functions...\n")
if (file.exists("DAR_ULTIMATE_COMPLETE_LOADER.r")) source("DAR_ULTIMATE_COMPLETE_LOADER.r") else
  stop("DAR_ULTIMATE_COMPLETE_LOADER.r not found.")
if (file.exists("DAR_RUNTIME_PATCH.r")) source("DAR_RUNTIME_PATCH.r")
cat("  Done.\n\n")

cat("[2/5] Loading TF partner peaks...\n")
partner_peaks <- list(
  Bcd_FDR1                  = import("./data/tf_partners/bdtnp_Bcd1-2_FDR1_all_sites_dm6.bed"),
  Ftz_FDR1                  = import("./data/tf_partners/bdtnp_Ftz3_FDR1_dm6.bed"),
  Ftz_FDR25                 = import("./data/tf_partners/bdtnp_Ftz3_FDR25_dm6.bed"),
  Cic_Union                 = import("./data/tf_partners/Cic_ChIPseq_published_union_annotated.bed"),
  Cic_sfGFP                 = import("./data/tf_partners/CicsfGFP_GSE130584_MACS2_peaks.bed"), 
  Hb_ChIPseq                = import("./data/tf_partners/ChIP_Hb_GSE50771_dm6.bed"),
  Gt_ChIPseq                = import("./data/tf_partners/ChIP_Gt_GSE50771_dm6.bed"),
  Prd_FDR1                  = import("./data/tf_partners/bdtnp_Prd_FDR1_dm6.bed"),
  Tll_FDR1                  = import("./data/tf_partners/bdtnp_Tll_FDR1_dm6.bed"),
  Zld_2hr                   = import("./data/tf_partners/ChIP_Zld_2hr_GSM763061_peaks_dm6.bed"),
  Zld_3hr                   = import("./data/tf_partners/not using/ChIP_Zld_3hr_GSM763062_peaks_dm6.bed"),  
  Hkb_FDR1                  = import("./data/tf_partners/bdtnp_Hkb1-3_FDR1_all_sites_dm6.bed"),
  D_FDR1                    = import("./data/tf_partners/bdtnp_D_FDR1_dm6.bed"),
  Cad_FDR1                  = import("./data/tf_partners/bdtnp_Cad_FDR1_dm6.bed"),
  Hairy_FDR1                = import("./data/tf_partners/bdtnp_Hairy1-2_FDR1_all_sites_dm6.bed"),
  Opa_early                 = import("./data/tf_partners/ChIP_Opa_Early_rep1-2_q01_overlapping_sites.bed"),
  Opa_late                  = import("./data/tf_partners/ChIP_Opa_Late_rep1-2_q01_overlapping_sites.bed"),
  Run_ChIP_chip_FDR1        = import("./data/chipchip/bdtnp_Run1-2_FDR1_all_sites.bed"),
  Runt_ChIPseq_Reproducible = import("./data/chipseq/RuntAb_IgG_p0.05_reproducible.bed")
)
cat(sprintf("  %d TF datasets loaded.\n\n", length(partner_peaks)))

cat("[3/5] Loading peak universe...\n")
all_peaks <- GRanges()
for (ap_path in c("./ATAC_DARS/SevenGeno_nc14b_counts_matrix.txt",
                  "./SevenGeno_nc14b_all_peaks.bed",
                  "./ATAC_NarrowPeaks/FullUniverse_nc14b_AllGeno_union_peaks.bed")) {
  if (!file.exists(ap_path)) next
  if (grepl("\\.bed$", ap_path)) {
    all_peaks <- import(ap_path)
  } else {
    hdr <- readLines(ap_path, n=1)
    sep <- if (grepl("\t", hdr, fixed=TRUE)) "\t" else " "
    mat <- read.table(ap_path, header=TRUE, sep=sep, stringsAsFactors=FALSE, check.names=FALSE)
    cn  <- colnames(mat)
    if (all(c("chr","start","end") %in% cn)) {
      all_peaks <- GRanges(seqnames=mat$chr,
                           ranges=IRanges(start=mat$start+1, end=mat$end))
    } else if (grepl(":", mat[1,1])) {
      coords <- do.call(rbind, strsplit(gsub("-",":",mat[,1]),":"))
      all_peaks <- GRanges(seqnames=coords[,1],
                           ranges=IRanges(start=as.integer(coords[,2])+1,
                                          end=as.integer(coords[,3])))
    }
  }
  if (length(all_peaks) > 0) {
    cat(sprintf("  Loaded %d peaks from %s\n\n", length(all_peaks), ap_path)); break
  }
}
if (length(all_peaks) == 0)
  cat("  WARNING: peak universe not found — using DAR-only background.\n\n")

################################################################################
# HELPERS
################################################################################

load_dar_bed <- function(stem, base_dir) {
  path <- file.path(base_dir, paste0(stem, "_DARs.bed"))
  if (!file.exists(path)) { message("  [missing] ", path); return(GRanges()) }
  df <- read.table(path, sep="\t", header=FALSE, stringsAsFactors=FALSE)
  gr <- GRanges(seqnames=df[,1], ranges=IRanges(start=df[,2]+1, end=df[,3]), strand="*")
  names(gr) <- df[,4]
  negate_lfc <- stem %in% NEGATE_LFC_STEMS
  lfc_attached <- FALSE

  # Search for annotated file in multiple locations:
  # 1. ATAC_DARS/ (after cp from SevenGeno output)
  # 2. The SevenGeno output dir directly
  # 3. The temporal output dir
  ann_candidates <- c(
    file.path(base_dir, paste0(stem, "_DARs_annotated.txt")),
    file.path(base_dir, paste0(stem, "_DARs_annotated_peaks.csv")),
    file.path(SEVENGENO_OUT_DIR, paste0(stem, "_DARs_annotated.txt")),
    file.path(VENT_TEMP_OUT_DIR, paste0(stem, "_DARs_annotated.txt"))
  )

  for (ann_path in ann_candidates) {
    if (!file.exists(ann_path)) next
    fsize <- file.info(ann_path)$size
    if (is.na(fsize) || fsize < 10) next

    # Read first line to detect separator and BOM
    first_line <- tryCatch(
      readLines(ann_path, n = 1, warn = FALSE, encoding = "UTF-8"),
      error = function(e) ""
    )
    # Strip UTF-8 BOM if present
    first_line <- gsub("^\xef\xbb\xbf", "", first_line)
    sep <- if (grepl("\t", first_line, fixed = TRUE)) "\t" else ","

    ann <- tryCatch(
      read.table(ann_path, header = TRUE, sep = sep,
                 stringsAsFactors = FALSE, quote = "\"",
                 fill = TRUE, comment.char = "",
                 row.names = NULL, fileEncoding = "UTF-8"),
      error = function(e) {
        # Try without encoding specification as fallback
        tryCatch(
          read.table(ann_path, header = TRUE, sep = sep,
                     stringsAsFactors = FALSE, quote = "\"",
                     fill = TRUE, comment.char = "", row.names = NULL),
          error = function(e2) {
            message(sprintf("    [read error] %s: %s", basename(ann_path), e2$message))
            NULL
          }
        )
      }
    )
    if (is.null(ann) || nrow(ann) == 0) next
    colnames(ann)[1] <- gsub("^[^[:alnum:]]+", "", colnames(ann)[1])
    lfc_col  <- intersect(c("log2FoldChange","log2FC","LFC","log2fc"), colnames(ann))[1]
    if (is.na(lfc_col)) next
    name_col <- intersect(c("name","peakID","ID","peak"), colnames(ann))[1]

    # Build lfc vector aligned to gr
    # Annotated file may be a subset of BED rows — always try name-based first
    if (!is.na(name_col) && length(names(gr)) > 0 && !all(is.na(names(gr)))) {
      lfc_named <- setNames(as.numeric(ann[[lfc_col]]), ann[[name_col]])
      matched   <- lfc_named[names(gr)]   # NA for peaks not in annotated file
    } else {
      matched <- rep(NA_real_, length(gr))
    }

    # Positional fallback only when sizes match exactly AND name matching gave all NA
    if (all(is.na(matched)) && nrow(ann) == length(gr)) {
      matched <- as.numeric(ann[[lfc_col]])
      cat(sprintf("    [load_dar_bed] %s: positional match (%d rows)\n", stem, nrow(ann)))
    }

    # ── Fallback chain: coordinate match, then index match ──────────────────
    n_matched <- sum(!is.na(matched))

    # Fallback 2: coordinate-based (handles name format mismatches)
    if (n_matched == 0 && all(c("chr","start","end") %in% colnames(ann))) {
      ann_key <- paste(ann$chr, ann$start, ann$end, sep=":")
      gr_key  <- paste(as.character(seqnames(gr)), start(gr)-1, end(gr), sep=":")
      coord_m <- setNames(as.numeric(ann[[lfc_col]]), ann_key)[gr_key]
      if (sum(!is.na(coord_m)) > 0) {
        matched   <- coord_m
        n_matched <- sum(!is.na(matched))
        cat(sprintf("    [load_dar_bed] %s: coordinate match (%d/%d)\n",
                    stem, n_matched, length(gr)))
      }
    }

    # Fallback 3: trailing integer index (e.g. "_1" vs "_DAR_1")
    if (n_matched == 0 && !is.na(name_col)) {
      gr_idx  <- suppressWarnings(as.integer(gsub(".*_(\\d+)$", "\\1", names(gr))))
      ann_idx <- suppressWarnings(as.integer(gsub(".*_(\\d+)$", "\\1", ann[[name_col]])))
      if (!all(is.na(gr_idx)) && !all(is.na(ann_idx))) {
        lfc_by_idx <- setNames(as.numeric(ann[[lfc_col]]), ann_idx)
        matched    <- lfc_by_idx[as.character(gr_idx)]
        n_matched  <- sum(!is.na(matched))
        if (n_matched > 0)
          cat(sprintf("    [load_dar_bed] %s: index match (%d/%d)\n",
                      stem, n_matched, length(gr)))
      }
    }

    if (n_matched == 0) next   # nothing worked — try next candidate file
    cat(sprintf("    [load_dar_bed] %s: %d/%d peaks with log2FC\n",
                stem, n_matched, length(gr)))
    if (negate_lfc) {
      matched <- -matched
      cat(sprintf("    [sign fix] %s: log2FC negated\n", stem))
    }
    mcols(gr)$log2FoldChange <- matched
    mcols(gr)$direction <- ifelse(matched > 0, "gained-open", "gained-close")
    lfc_attached <- TRUE; break
  }
  if (!lfc_attached && ncol(df) >= 5) {
    scores <- as.numeric(df[,5])
    if (any(scores < 0, na.rm=TRUE)) {
      if (negate_lfc) scores <- -scores
      mcols(gr)$log2FoldChange <- scores
      mcols(gr)$direction <- ifelse(scores > 0, "gained-open", "gained-close")
      lfc_attached <- TRUE
    }
  }
  if (!lfc_attached) {
    message(sprintf("  [WARNING] No real log2FC for %s", stem))
    # Diagnose: check if any annotated file was found but matching failed
    for (ann_path in ann_candidates) {
      if (file.exists(ann_path) && file.info(ann_path)$size > 10) {
        hdr <- tryCatch(colnames(read.table(ann_path, nrows=1, header=TRUE,
                                            sep="\t", stringsAsFactors=FALSE,
                                            check.names=FALSE)),
                        error=function(e) NULL)
        if (!is.null(hdr)) {
          message(sprintf("    Found: %s", ann_path))
          message(sprintf("    Columns: %s", paste(hdr, collapse=", ")))
          message(sprintf("    log2FC col found: %s",
                          paste(intersect(c("log2FoldChange","log2FC","LFC","log2fc"), hdr), collapse=",")))
        }
        break
      }
    }
    if (!any(sapply(ann_candidates, function(p) file.exists(p) && !is.na(file.info(p)$size) && file.info(p)$size > 10))) {
      message(sprintf("    Annotated file searched in %d locations:", length(ann_candidates)))
      for (p in ann_candidates) {
        exists  <- file.exists(p)
        sz      <- if (exists) file.info(p)$size else NA
        message(sprintf("      [%s size=%s] %s",
                        if (exists) "EXISTS" else "missing",
                        if (!is.na(sz)) sz else "NA", p))
      }
    }
  }
  gr
}

# ── Bootstrap/stability summary — ONLY valid when runs used genuinely
#    different data (n_peaks_total > LMM_N_SUB, so sample() actually
#    subsampled). Reports empirical between-run spread, not a Stouffer
#    combination — real subsamples still overlap heavily peak-to-peak and
#    are NOT independent draws, so a naive z-sum is still not appropriate.
#    This is a stability diagnostic, not a replicated-experiment meta-analysis.
combine_runs_bootstrap <- function(run_list) {
  tfs <- Reduce(intersect, lapply(run_list, function(r) r$TF[r$TF != "(Intercept)"]))
  if (length(tfs) == 0) return(NULL)

  do.call(rbind, lapply(tfs, function(tf) {
    rows <- lapply(run_list, function(r) r[r$TF == tf, , drop=FALSE])
    ests <- sapply(rows, function(r) r$estimate[1])
    ests <- ests[!is.na(ests)]

    est       <- mean(ests)
    se_empir  <- sd(ests) / sqrt(length(ests))   # SE of the mean across bootstrap draws
    ci_lo     <- unname(quantile(ests, 0.025))
    ci_hi     <- unname(quantile(ests, 0.975))

    data.frame(TF = tf, estimate = est, se = se_empir,
               ci_lo_empirical = ci_lo, ci_hi_empirical = ci_hi,
               n_runs = length(ests), stringsAsFactors = FALSE)
  }))
}

# ── Collinearity flag across the TF design matrix used in a stratum.
#    partner_peaks includes multiple FDR thresholds of the same ChIP/BDTNP
#    dataset for some TFs (e.g. Ftz_FDR1 vs Ftz_FDR25) — these are expected
#    to correlate and can destabilize fixed-effect estimates if both are
#    kept in the same model.
tf_collinearity_flag <- function(model_data, keep_tf, thresh = 0.80) {
  if (length(keep_tf) < 2) return(list(max_abs_r = NA_real_, flagged_pairs = ""))
  cm <- suppressWarnings(cor(model_data[, keep_tf, drop = FALSE], use = "pairwise.complete.obs"))
  cm[lower.tri(cm, diag = TRUE)] <- NA
  idx <- which(abs(cm) > thresh, arr.ind = TRUE)
  pairs <- if (nrow(idx) > 0) {
    paste(sprintf("%s~%s(r=%.2f)", rownames(cm)[idx[,1]], colnames(cm)[idx[,2]],
                  cm[idx]), collapse = "; ")
  } else ""
  list(max_abs_r = if (all(is.na(cm))) NA_real_ else max(abs(cm), na.rm = TRUE),
       flagged_pairs = pairs)
}

################################################################################
# BOOTSTRAP STABILITY CHECK  (new in v7)
#
# This is deliberately SEPARATE from the primary single-fit result above.
# It answers a different question: "is this coefficient sensitive to which
# peaks happen to be in the dataset," not "is this coefficient significant."
# It resamples WITH REPLACEMENT (a real nonparametric case bootstrap), which
# is what actually produces a genuinely different dataset on every draw —
# unlike the original v4 script's without-replacement sampling, which
# returned the full dataset every time whenever n_peaks_total <= LMM_N_SUB.
#
# Output (per stratum, per model type): a bootstrap mean, bootstrap SD,
# percentile 95% CI, and a sign-consistency score (fraction of bootstrap
# fits whose coefficient sign matches the primary single-fit's sign) for
# every TF. These are written to a SEPARATE file
# (lmm_bootstrap_stability.txt / _binary.txt) — they are not combined into
# the primary p-value/significance call, so this cannot reintroduce the
# original Stouffer-style inflation.
#
# Limitation worth noting in the methods if you report this: `chr` has only
# ~5-6 levels in this dataset (Drosophila chromosome arms), too few for a
# meaningful cluster/block bootstrap at the chromosome level, so this
# resamples individual peaks rather than whole chromosomes. That is standard
# practice here, but means the bootstrap SD may modestly understate
# uncertainty if peaks within a chromosome are more correlated than the
# random-intercept term captures.
################################################################################

RUN_BOOTSTRAP_STABILITY <- TRUE   # set FALSE to skip this section entirely
N_BOOTSTRAP              <- 200   # resamples per stratum per model type.
                                   # 200 is a reasonable default for stable
                                   # SD/sign-consistency estimates; the
                                   # percentile CI tails are noisier at 200
                                   # than at, say, 1000, but 1000 x (every
                                   # stratum x both model types) is likely
                                   # too slow to run locally in one sitting.
                                   # Lower to ~50 for a first smoke-test of
                                   # runtime before committing to a full run.

# Optionally restrict the (slow) bootstrap to specific comparisons while you
# smoke-test runtime -- e.g. just your Figure 5 panels first. Leave both NULL
# to run on everything RUN_LMM_ONLY normally would.
BOOTSTRAP_ONLY_KEYS   <- NULL   # e.g. c("BOTv_vs_BOTCv", "BOTv_vs_BOTCv_late")
BOOTSTRAP_ONLY_THRESH <- NULL   # e.g. c("all_peaks")

bootstrap_stability_check <- function(formula, data, family_glmer = NULL,
                                       primary_coef, n_boot = N_BOOTSTRAP) {
  # family_glmer = NULL -> lmer (continuous); non-NULL -> glmer(family=...)
  is_glmm <- !is.null(family_glmer)
  n       <- nrow(data)
  boot_ests <- vector("list", n_boot)
  n_ok      <- 0L
  t0        <- Sys.time()

  for (b in seq_len(n_boot)) {
    d <- data[sample(n, n, replace = TRUE), ]   # WITH replacement -- genuinely differs every draw
    fit <- tryCatch({
      if (is_glmm) {
        glmer(formula, data = d, family = family_glmer,
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
      } else {
        lmer(formula, data = d, REML = FALSE,
             control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
      }
    }, error = function(e) NULL, warning = function(w) NULL)
    if (is.null(fit)) next

    ct <- tryCatch(as.data.frame(coef(summary(fit))), error = function(e) NULL)
    if (is.null(ct)) next
    ct$TF <- rownames(ct)
    names(ct)[names(ct) == "Estimate"] <- "estimate"
    boot_ests[[b]] <- ct[ct$TF != "(Intercept)", c("TF", "estimate")]
    n_ok <- n_ok + 1L
  }
  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  cat(sprintf("    [bootstrap] %d/%d resamples fit successfully (%.1fs)\n", n_ok, n_boot, elapsed))

  valid <- Filter(Negate(is.null), boot_ests)
  if (length(valid) == 0) return(NULL)
  long <- do.call(rbind, valid)

  do.call(rbind, lapply(split(long, long$TF), function(d) {
    tf <- d$TF[1]
    primary_sign <- sign(primary_coef$estimate[primary_coef$TF == tf][1])
    data.frame(
      TF                 = tf,
      boot_mean          = mean(d$estimate),
      boot_sd            = sd(d$estimate),
      boot_ci_lo         = unname(quantile(d$estimate, 0.025)),
      boot_ci_hi         = unname(quantile(d$estimate, 0.975)),
      pct_sign_consistent = round(100 * mean(sign(d$estimate) == primary_sign, na.rm = TRUE), 1),
      n_boot_used        = nrow(d),
      stringsAsFactors   = FALSE
    )
  }))
}

# ── Save plot ─────────────────────────────────────────────────────────────────
# is_binary=TRUE triggers binary-specific axis labels and caption.
# pos_color / neg_color: when both supplied, significant dots are colored by
#   direction (pos_color = positive estimate, neg_color = negative estimate)
#   and non-significant dots are grey. Used for BOTv_vs_BOTCv comparisons where
#   pink (BOTCv, Cic absent) = positive and green (BOTv, Cic present) = negative.
#   When NULL, falls back to single-color significance scheme (color arg).
save_coef_plot <- function(coef_plot, lmm_dir, filename, label, color,
                           n_runs, n_obs_mean, n_tfs,
                           is_binary = FALSE, n_peaks_total = NA,
                           pos_color = NULL, neg_color = NULL,
                           real_subsample = FALSE) {
  coef_plot <- coef_plot[coef_plot$TF != "(Intercept)", , drop=FALSE]
  if (nrow(coef_plot) == 0) return(invisible(NULL))

  has_padj <- "p_adj" %in% names(coef_plot) && !all(is.na(coef_plot$p_adj))
  coef_plot$sig_raw <- !is.na(coef_plot$p_value) & coef_plot$p_value < 0.05 &
                       abs(coef_plot$estimate) > 0.20
  coef_plot$sig_adj <- if (has_padj)
                         !is.na(coef_plot$p_adj) & coef_plot$p_adj < 0.05 &
                         abs(coef_plot$estimate) > 0.20
                       else FALSE

  n_sig_raw <- sum(coef_plot$sig_raw, na.rm = TRUE)
  n_sig_adj <- sum(coef_plot$sig_adj, na.rm = TRUE)

  size_note <- if (!is.na(n_peaks_total))
                 sprintf(" | n_peaks_total=%d", n_peaks_total)
               else ""
  sub <- sprintf("n_runs=%d | mean n_obs=%.0f | %d TFs | %d sig (raw) | %d sig (BH)%s",
                 n_runs, n_obs_mean, n_tfs, n_sig_raw, n_sig_adj, size_note)

  # Caption/label now depend on whether the coefficients came from a single
  # full-data fit (no real subsampling was possible: n_peaks_total <=
  # LMM_N_SUB) or a genuine bootstrap over subsamples. These are different
  # things and were previously described identically ("averaged over N
  # runs" / "Stouffer Z" / "variance across subsampling runs") even when
  # N==1 distinct dataset was fit N times.
  method_note <- if (real_subsample) {
    sprintf(paste0(
      "Estimate = mean across %d bootstrap subsamples (n=%d each, drawn without\n",
      "replacement from %s peaks — NOT independent draws, so treat as a stability\n",
      "diagnostic, not a replicated-experiment meta-analysis).\n",
      "p-value: one-sample t-test of the bootstrap mean vs 0, BH-corrected.\n",
      "Bars = empirical 95%% CI across the %d bootstrap runs (true run-to-run spread)."),
      n_runs, n_obs_mean, format(n_peaks_total, big.mark=","), n_runs)
  } else {
    sprintf(paste0(
      "Estimate = single fit on the full stratum (n=%s peaks; no subsampling was\n",
      "performed because n_peaks_total <= LMM_N_SUB).\n",
      "p-value: model-based (Satterthwaite df for LMM / Wald z for binary GLMM), BH-corrected.\n",
      "Bars = model-based 95%% CI on the coefficient (NOT run-to-run variance)."),
      format(n_peaks_total, big.mark=","))
  }

  if (is_binary) {
    x_label <- "Log-odds of Opening vs Closing"
    caption <- paste0(
      "Coefficients are log-odds from a logistic GLMM (binomial family, chr random intercept).\n",
      "Positive = TF binding enriched at Opening peaks; Negative = enriched at Closing peaks.\n",
      "Convert to odds ratio via exp(coef): e.g. coef=2 \u2192 OR\u22487.4\u00d7 more likely to open.\n",
      "Significance threshold: |log-odds| > 0.20 AND p_adj < 0.05.\n", method_note
    )
  } else {
    x_label <- "LMM coefficient on log2FC"
    caption <- paste0(
      "Coefficients are fixed effects from a linear mixed model (chr random intercept).\n",
      "Positive = TF binding associated with more open chromatin (higher log2FC).\n",
      "Significance threshold: |beta| > 0.20 AND p_adj < 0.05.\n", method_note
    )
  }

  use_dir_colors <- !is.null(pos_color) && !is.null(neg_color)

  if (use_dir_colors) {
    # Directional coloring: significant dots colored by direction, non-sig grey
    # pos_color = BOTCv/pink (Cic absent opens), neg_color = BOTv/green (Cic present opens)
    coef_plot$dot_color <- ifelse(
      !coef_plot$sig_raw, "ns",
      ifelse(coef_plot$estimate >= 0, "pos", "neg")
    )
    color_values <- c("ns" = "grey65", "pos" = pos_color, "neg" = neg_color)
    color_labels <- c("ns" = "Not significant",
                      "pos" = "Cic absent opens (BOTCv)",
                      "neg" = "Cic present opens (BOTv)")

    # Append direction key to caption
    caption <- paste0(caption, "\n",
      sprintf("Pink (%s) = opens when Cic absent (BOTCv); ",  pos_color),
      sprintf("Green (%s) = opens when Cic present (BOTv).", neg_color))

    p <- ggplot(coef_plot, aes(x=estimate, y=reorder(TF, estimate),
                                colour=dot_color)) +
      annotate("rect", xmin=-0.20, xmax=0.20, ymin=-Inf, ymax=Inf,
               fill="grey92", alpha=0.9) +
      geom_vline(xintercept=0, linewidth=0.4, colour="grey55") +
      geom_point(aes(size=sig_raw), shape=16) +
      {if (!all(is.na(coef_plot$se)))
         geom_errorbarh(aes(xmin=estimate-1.96*se, xmax=estimate+1.96*se),
                        height=0.3, linewidth=0.4)
       else list()} +
      scale_colour_manual(values=color_values, labels=color_labels,
                          name=NULL) +
      scale_size_manual(values=c("FALSE"=2.0, "TRUE"=3.2), guide="none") +
      guides(colour = guide_legend(override.aes = list(size=3))) +
      theme_minimal(base_size=10) +
      theme(
        axis.text.y      = element_text(size=8),
        panel.grid.minor = element_blank(),
        plot.subtitle    = element_text(size=8, colour="grey40"),
        plot.caption     = element_text(size=6.5, colour="grey50",
                                         hjust=0, lineheight=1.3),
        legend.position  = "bottom"
      ) +
      labs(title=label, subtitle=sub, caption=caption, x=x_label, y=NULL)

  } else {
    # Default: single-color significance scheme
    p <- ggplot(coef_plot, aes(x=estimate, y=reorder(TF, estimate),
                                colour=sig_raw)) +
      annotate("rect", xmin=-0.20, xmax=0.20, ymin=-Inf, ymax=Inf,
               fill="grey92", alpha=0.9) +
      geom_vline(xintercept=0, linewidth=0.4, colour="grey55") +
      geom_point(aes(size=sig_raw), shape=16) +
      {if (!all(is.na(coef_plot$se)))
         geom_errorbarh(aes(xmin=estimate-1.96*se, xmax=estimate+1.96*se),
                        height=0.3, linewidth=0.4)
       else list()} +
      scale_colour_manual(values=c("FALSE"="grey65", "TRUE"=color), guide="none") +
      scale_size_manual(values=c("FALSE"=2.0, "TRUE"=3.2), guide="none") +
      theme_minimal(base_size=10) +
      theme(
        axis.text.y      = element_text(size=8),
        panel.grid.minor = element_blank(),
        plot.subtitle    = element_text(size=8, colour="grey40"),
        plot.caption     = element_text(size=6.5, colour="grey50",
                                         hjust=0, lineheight=1.3)
      ) +
      labs(title=label, subtitle=sub, caption=caption, x=x_label, y=NULL)
  }

  tryCatch(ggsave(file.path(lmm_dir, filename), p,
                  width=7.5, height=max(4, nrow(coef_plot)*0.35+3.5)),
           error=function(e) message("    plot warning: ", conditionMessage(e)))
}


################################################################################
# GLMMLASSO HELPERS
# Aligned to original 3_analysis_core_functions.r logic:
#   - BIC-guided lambda selection over linear grid seq(0.01, 1, length.out=20)
#   - Manual standardize (scale) before fitting, back-transform after
#   - selected = abs(coef_STANDARDIZED) > 1e-6  (selection on std scale)
#   - Binary model: standardize TF predictors only, keep 0/1 response as-is
#   - glmmLasso called WITHOUT the family= argument for gaussian (default)
################################################################################

# ── Fit a single glmmLasso model safely ---------------------------------------
fit_glmmlasso_safe <- function(formula_obj, data, lambda, is_binomial = FALSE) {
  withCallingHandlers(
    tryCatch({
      if (is_binomial) {
        glmmLasso(
          fix     = formula_obj,
          rnd     = list(chr = ~1),
          data    = data,
          lambda  = lambda,
          family  = binomial(link = "logit"),
          control = list(print.iter = FALSE, standardize = FALSE)
        )
      } else {
        # Gaussian: omit family argument entirely (glmmLasso default = gaussian)
        glmmLasso(
          fix     = formula_obj,
          rnd     = list(chr = ~1),
          data    = data,
          lambda  = lambda,
          control = list(print.iter = FALSE, standardize = FALSE)
        )
      }
    }, error = function(e) NULL),
    warning = function(w) invokeRestart("muffleWarning")
  )
}

# ── BIC-guided lambda selection -----------------------------------------------
# Sweeps the full grid, reads fit$bic from each, returns min-BIC lambda.
# Diagnostic BIC path is printed so failures are visible in the log.
select_lambda_bic <- function(formula_obj, data, lambda_grid, is_binomial = FALSE) {
  bic_vals <- vapply(lambda_grid, function(lam) {
    fit <- fit_glmmlasso_safe(formula_obj, data, lam, is_binomial)
    if (is.null(fit)) return(Inf)
    bic <- tryCatch(as.numeric(fit$bic), error = function(e) Inf)
    if (length(bic) == 0 || is.na(bic) || !is.finite(bic)) Inf else bic
  }, numeric(1))

  # Diagnostic: print the BIC path so we can see what's happening
  finite_bics <- bic_vals[is.finite(bic_vals)]
  if (length(finite_bics) > 0) {
    cat(sprintf("    BIC path: %d/%d lambdas converged | BIC range [%.1f, %.1f]\n",
                length(finite_bics), length(lambda_grid),
                min(finite_bics), max(finite_bics)))
  } else {
    cat("    BIC path: WARNING — no lambdas converged (all BIC = Inf)\n")
  }

  best_idx    <- which.min(bic_vals)
  # If nothing converged, fall back to the smallest lambda (least regularized)
  best_lambda <- if (length(best_idx) > 0 && is.finite(bic_vals[best_idx]))
                   lambda_grid[[best_idx]]
                 else lambda_grid[1]  # seq goes low->high, so [1] = 0.01

  list(best_lambda = best_lambda,
       bic_df      = data.frame(lambda = lambda_grid, bic = bic_vals))
}

# ── Standardize TF predictors (and optionally response) -----------------------
# Returns scaled data + the scale params needed for back-transformation.
standardize_model_data <- function(model_data, tf_names, scale_response = TRUE,
                                   response = "logFC") {
  tf_means <- sapply(model_data[tf_names], mean, na.rm = TRUE)
  tf_sds   <- sapply(model_data[tf_names], sd,   na.rm = TRUE)
  resp_mean <- mean(model_data[[response]], na.rm = TRUE)
  resp_sd   <- sd(  model_data[[response]], na.rm = TRUE)

  scaled <- model_data
  for (tf in tf_names) {
    sd_tf <- tf_sds[tf]
    scaled[[tf]] <- if (!is.na(sd_tf) && sd_tf > 0)
                      (model_data[[tf]] - tf_means[tf]) / sd_tf
                    else 0
  }
  if (scale_response && !is.na(resp_sd) && resp_sd > 0) {
    scaled[[response]] <- (model_data[[response]] - resp_mean) / resp_sd
  }

  list(data      = scaled,
       tf_means  = tf_means,
       tf_sds    = tf_sds,
       resp_mean = resp_mean,
       resp_sd   = resp_sd)
}

# ── Back-transform standardized coefficients to original scale ----------------
backtransform_coefs <- function(coef_std, coef_names, scale_params) {
  resp_sd   <- scale_params$resp_sd
  tf_sds    <- scale_params$tf_sds
  resp_mean <- scale_params$resp_mean

  result        <- rep(NA_real_, length(coef_std))
  names(result) <- coef_names

  for (i in seq_along(coef_std)) {
    nm <- coef_names[i]
    if (nm == "(Intercept)") {
      result[i] <- coef_std[i] * resp_sd + resp_mean
    } else if (nm %in% names(tf_sds) && !is.na(tf_sds[nm]) && tf_sds[nm] > 0) {
      result[i] <- coef_std[i] * (resp_sd / tf_sds[nm])
    } else {
      result[i] <- coef_std[i] * resp_sd
    }
  }
  result
}

# ── Extract coefficient table -------------------------------------------------
# KEY FIX: selection criterion applied to STANDARDIZED coefficients.
# For continuous: back-transformed estimates stored in 'estimate';
#   standardized stored in 'estimate_std'. p_value = NA (matches original).
# For binary: coefficients are on log-odds scale (no back-transform).
#   SE and Wald z/p extracted if available.
extract_glmmlasso_coefs <- function(fit, lambda_used, scale_params = NULL,
                                     is_binary = FALSE) {
  if (is.null(fit)) return(NULL)

  coef_std   <- as.numeric(fit$coefficients)
  coef_names <- names(fit$coefficients)
  if (length(coef_std) == 0 || is.null(coef_names)) return(NULL)
  names(coef_std) <- coef_names

  tfs <- coef_names[coef_names != "(Intercept)"]
  if (length(tfs) == 0) return(NULL)

  # Selection is always on standardized scale
  selected_vec <- abs(coef_std[tfs]) > 1e-6

  if (!is_binary && !is.null(scale_params)) {
    coef_orig <- backtransform_coefs(coef_std, coef_names, scale_params)
    estimate_vec <- as.numeric(coef_orig[tfs])
  } else {
    estimate_vec <- as.numeric(coef_std[tfs])
  }

  df <- data.frame(
    TF           = tfs,
    estimate     = estimate_vec,
    estimate_std = as.numeric(coef_std[tfs]),
    selected     = selected_vec,
    lambda_used  = lambda_used,
    stringsAsFactors = FALSE
  )

  if (is_binary) {
    # Attempt Wald SE from StdErr matrix
    ses <- tryCatch({
      se_mat <- fit$StdErr
      if (!is.null(se_mat) && nrow(se_mat) == length(coef_std)) {
        se_vec        <- sqrt(diag(se_mat))
        names(se_vec) <- coef_names
        as.numeric(se_vec[tfs])
      } else {
        rep(NA_real_, length(tfs))
      }
    }, error = function(e) rep(NA_real_, length(tfs)))

    df$se      <- ses
    df$z_value <- ifelse(!is.na(ses) & ses > 0, df$estimate_std / ses, NA_real_)
    df$p_value <- ifelse(!is.na(df$z_value), 2 * pnorm(-abs(df$z_value)), NA_real_)
  } else {
    df$se      <- NA_real_
    df$z_value <- NA_real_
    df$p_value <- NA_real_
  }

  rownames(df) <- NULL
  df
}

# ── Plot (matches reference plot) ---------------------------------------------
# Coloring = "selected" (LASSO shrinkage), not p-values.
save_glmmlasso_plot <- function(coef_df, out_path, title_str, color,
                                 n_peaks, best_lambda, model_label,
                                 effect_thresh = LASSO_EFFECT_THRESH) {
  if (is.null(coef_df) || nrow(coef_df) == 0) return(invisible(NULL))

  df <- coef_df[coef_df$TF != "(Intercept)", , drop = FALSE]
  if (nrow(df) == 0) return(invisible(NULL))

  # Guard: if selected is all NA (no convergence), treat as all FALSE
  if (all(is.na(df$selected))) df$selected <- FALSE
  df$significant <- df$selected

  n_sig <- sum(df$significant, na.rm = TRUE)

  sub_str <- sprintf(
    "%d TFs analyzed | %d significant | Model: Linear Mixed Model",
    nrow(df), n_sig
  )
  footer_str <- sprintf(
    "Gray zone: |effect| < %.3f | Method: regularized_glmm_lasso | lambda=%.4f | n_peaks=%d",
    effect_thresh, best_lambda, n_peaks
  )

  has_se <- "se" %in% names(df) && !all(is.na(df$se))

  p <- ggplot(df, aes(x = estimate,
                       y = reorder(TF, estimate),
                       colour = significant)) +
    annotate("rect",
             xmin = -effect_thresh, xmax = effect_thresh,
             ymin = -Inf,           ymax  = Inf,
             fill = "grey92", alpha = 0.9) +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey55") +
    geom_point(aes(size = significant), shape = 16) +
    {if (has_se)
       geom_errorbarh(aes(xmin = estimate - 1.96 * se,
                           xmax = estimate + 1.96 * se),
                      height = 0.3, linewidth = 0.4)
     else list()} +
    scale_colour_manual(
      values = c("FALSE" = "#457B9D", "TRUE" = color),
      labels = c("FALSE" = "Not significant", "TRUE" = "Significant"),
      name   = NULL
    ) +
    scale_size_manual(values = c("FALSE" = 2.0, "TRUE" = 3.2), guide = "none") +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.y      = element_text(size = 8),
      panel.grid.minor = element_blank(),
      plot.subtitle    = element_text(size = 7.5, colour = "grey40"),
      plot.caption     = element_text(size = 7,   colour = "grey55"),
      legend.position  = "bottom"
    ) +
    labs(
      title    = title_str,
      subtitle = sub_str,
      caption  = footer_str,
      x        = "Effect Size on Chromatin Accessibility",
      y        = "Transcription Factor"
    )

  tryCatch(
    ggsave(out_path, p, width = 7,
           height = max(4, nrow(df) * 0.35 + 2.5)),
    error = function(e) message("    glmmLasso plot warning: ", conditionMessage(e))
  )
}

# ── Main glmmLasso runner (called once per comparison x lfc_thresh) -----------
run_glmmlasso <- function(key, model_data_cont, thresh_dir,
                           label, color, lfc_thresh, keep_tf,
                           force_rerun = FALSE) {

  if (!GLMMLASSO_AVAILABLE) return(invisible(NULL))
  if (length(keep_tf) == 0)  return(invisible(NULL))

  thresh_tag     <- sprintf("lfc%.2f", lfc_thresh)
  cont_coef_file <- file.path(thresh_dir, "glmmlasso_continuous_coefs.txt")
  bin_coef_file  <- file.path(thresh_dir, "glmmlasso_binary_coefs.txt")
  cont_pdf       <- file.path(thresh_dir,
                     sprintf("glmmlasso_continuous_%s.pdf", thresh_tag))
  bin_pdf        <- file.path(thresh_dir,
                     sprintf("glmmlasso_binary_%s.pdf", thresh_tag))

  run_cont <- force_rerun || !file.exists(cont_coef_file)
  run_bin  <- force_rerun || !file.exists(bin_coef_file)

  if (!run_cont && !run_bin) {
    cat(sprintf("  SKIP glmmLasso (done): %-40s [%s]\n", key, thresh_tag))
    return(invisible(NULL))
  }

  # Fixed subsample for glmmLasso — LASSO_N_SUB is a hard cap applied equally
  # to ALL comparisons regardless of peak set size. This ensures BIC-selected
  # lambda is comparable across comparisons and that large sets don't dominate.
  n_peaks_total <- nrow(model_data_cont)
  if (n_peaks_total < LASSO_N_SUB) {
    cat(sprintf("    NOTE: peak set (n=%d) smaller than LASSO_N_SUB=%d — using full set.\n",
                n_peaks_total, LASSO_N_SUB))
  }
  n_sub <- min(LASSO_N_SUB, n_peaks_total)
  mdat  <- model_data_cont[sample(nrow(model_data_cont), n_sub, replace = FALSE), ]
  mdat$chr <- as.factor(mdat$chr)

  # Linear lambda grid low->high: seq(0.01, 1, length.out=20)
  # BIC tends to decrease then increase; sweeping low->high is fine.
  lambda_grid <- seq(LASSO_LAMBDA_END, LASSO_LAMBDA_START, length.out = LASSO_N_LAMBDA)

  # ── (a) Continuous model ─────────────────────────────────────────────────────
  if (run_cont) {
    cat(sprintf("  [glmmLasso-cont] %-40s  |logFC|>%.2f  n=%d\n",
                key, lfc_thresh, nrow(mdat)))

    # Standardize both TF predictors and logFC response
    sc      <- standardize_model_data(mdat, keep_tf,
                                       scale_response = TRUE, response = "logFC")
    mdat_sc <- sc$data

    formula_cont <- as.formula(
      paste0("logFC ~ ", paste(keep_tf, collapse = " + ")))

    bic_res <- tryCatch(
      select_lambda_bic(formula_cont, mdat_sc, lambda_grid, is_binomial = FALSE),
      error = function(e) {
        cat(sprintf("    BIC sweep error (cont): %s\n", conditionMessage(e))); NULL
      }
    )

    if (!is.null(bic_res)) {
      best_lam <- bic_res$best_lambda
      cat(sprintf("    Best lambda (cont, BIC): %.4f\n", best_lam))

      fit_cont  <- fit_glmmlasso_safe(formula_cont, mdat_sc, best_lam,
                                       is_binomial = FALSE)
      coef_cont <- extract_glmmlasso_coefs(fit_cont, best_lam,
                                            scale_params = sc, is_binary = FALSE)

      if (!is.null(coef_cont) && nrow(coef_cont) > 0) {
        coef_cont$model_type    <- "glmmlasso_continuous"
        coef_cont$lfc_threshold <- lfc_thresh
        coef_cont$n_peaks_total <- n_peaks_total
        write.table(coef_cont, cont_coef_file,
                    sep = "\t", quote = FALSE, row.names = FALSE)
        cat(sprintf("    Saved: %s  (%d selected)\n",
                    cont_coef_file, sum(coef_cont$selected, na.rm = TRUE)))

        save_glmmlasso_plot(
          coef_df     = coef_cont,
          out_path    = cont_pdf,
          title_str   = sprintf("TF Effects on Chromatin Accessibility\n%s  [|logFC|>%.2f]",
                                 label, lfc_thresh),
          color       = color,
          n_peaks     = nrow(mdat),
          best_lambda = best_lam,
          model_label = "regularized_glmm_lasso"
        )
        cat(sprintf("    Saved: %s\n", cont_pdf))
      } else {
        cat("    glmmLasso-cont: no coefficients extracted.\n")
      }
    }
  }

  # ── (b) Binary Opening vs Closing ────────────────────────────────────────────
  if (run_bin && RUN_BINARY) {
    bin_data        <- mdat
    bin_data$binary <- as.integer(bin_data$logFC > 0)

    has_both <- length(unique(bin_data$binary)) == 2 &&
                sum( bin_data$binary) >= 10 &&
                sum(!bin_data$binary) >= 10

    if (!has_both) {
      cat(sprintf("  SKIP glmmLasso-bin [%s]: insufficient class balance\n", key))
    } else {
      cat(sprintf("  [glmmLasso-bin]  %-40s  |logFC|>%.2f  n=%d\n",
                  key, lfc_thresh, nrow(bin_data)))

      # Standardize TF predictors only — binary response stays as 0/1
      sc_bin  <- standardize_model_data(bin_data, keep_tf,
                                         scale_response = FALSE, response = "logFC")
      bin_sc  <- sc_bin$data
      bin_sc$binary <- bin_data$binary  # ensure 0/1 column is present and untouched

      formula_bin <- as.formula(
        paste0("binary ~ ", paste(keep_tf, collapse = " + ")))

      bic_res_bin <- tryCatch(
        select_lambda_bic(formula_bin, bin_sc, lambda_grid, is_binomial = TRUE),
        error = function(e) {
          cat(sprintf("    BIC sweep error (bin): %s\n", conditionMessage(e))); NULL
        }
      )

      if (!is.null(bic_res_bin)) {
        best_lam_bin <- bic_res_bin$best_lambda
        cat(sprintf("    Best lambda (bin, BIC): %.4f\n", best_lam_bin))

        fit_bin  <- fit_glmmlasso_safe(formula_bin, bin_sc, best_lam_bin,
                                        is_binomial = TRUE)
        coef_bin <- extract_glmmlasso_coefs(fit_bin, best_lam_bin,
                                             scale_params = NULL, is_binary = TRUE)

        if (!is.null(coef_bin) && nrow(coef_bin) > 0) {
          coef_bin$model_type    <- "glmmlasso_binary"
          coef_bin$lfc_threshold <- lfc_thresh
          coef_bin$n_peaks_total <- n_peaks_total
          write.table(coef_bin, bin_coef_file,
                      sep = "\t", quote = FALSE, row.names = FALSE)
          cat(sprintf("    Saved: %s  (%d selected)\n",
                      bin_coef_file, sum(coef_bin$selected, na.rm = TRUE)))

          save_glmmlasso_plot(
            coef_df     = coef_bin,
            out_path    = bin_pdf,
            title_str   = sprintf("TF Effects: Opening vs Closing\n%s  [|logFC|>%.2f]",
                                   label, lfc_thresh),
            color       = color,
            n_peaks     = nrow(bin_data),
            best_lambda = best_lam_bin,
            model_label = "regularized_glmm_lasso (binomial)"
          )
          cat(sprintf("    Saved: %s\n", bin_pdf))
        } else {
          cat("    glmmLasso-bin: no coefficients extracted.\n")
        }
      }
    }
  }

  invisible(NULL)
}

# ── Main multi-run LMM runner (loops over LFC_THRESHOLDS) ───────────────────
run_and_save_lmm <- function(key, gr, lmm_dir, label, color="#333333",
                              pos_color=NULL, neg_color=NULL) {

  if (length(gr) == 0) { cat("  SKIP (0 DARs):", key, "\n"); return(invisible(NULL)) }

  # ── Build full model data (once, outside threshold loop) ─────────────────
  bg <- if (length(all_peaks) > 0) all_peaks else gr
  tf_mat <- as.data.frame(lapply(partner_peaks, function(pp)
    as.integer(countOverlaps(bg, pp) > 0)))

  ov      <- findOverlaps(bg, gr, minoverlap=1)
  lfc_vec <- rep(NA_real_, length(bg))
  if ("log2FoldChange" %in% names(mcols(gr)))
    lfc_vec[queryHits(ov)] <- mcols(gr)$log2FoldChange[subjectHits(ov)]

  model_data_full <- cbind(
    data.frame(logFC  = lfc_vec,
               binary = as.integer(!is.na(lfc_vec) & lfc_vec > 0),
               chr    = as.character(seqnames(bg)),
               stringsAsFactors = FALSE),
    tf_mat)

  # ── Loop over LFC thresholds ─────────────────────────────────────────────
  for (lfc_thresh in LFC_THRESHOLDS) {

    thresh_tag  <- sprintf("lfc%.2f", lfc_thresh)
    thresh_dir  <- file.path(lmm_dir, thresh_tag)
    dir.create(thresh_dir, recursive = TRUE, showWarnings = FALSE)

    coef_file <- file.path(thresh_dir, "lmm_coefficients.txt")
    in_force  <- key %in% FORCE_RERUN_KEYS
    if (file.exists(coef_file) && !in_force) {
      cat(sprintf("  SKIP (done): %-40s [%s]\n", key, thresh_tag))
      next
    }
    if (in_force && file.exists(coef_file)) {
      cat(sprintf("  FORCE RERUN: %-40s [%s]\n", key, thresh_tag))
      file.remove(coef_file)
      rds <- file.path(thresh_dir, "lmm_result.rds")
      if (file.exists(rds)) file.remove(rds)
    }

    model_data_cont <- model_data_full[
      !is.na(model_data_full$logFC) & abs(model_data_full$logFC) > lfc_thresh, ,
      drop = FALSE]

    run_lmm_stratum(
      key             = key,
      model_data_cont = model_data_cont,
      n_peaks_total   = nrow(model_data_cont),
      thresh_dir      = thresh_dir,
      thresh_tag      = thresh_tag,
      lfc_thresh      = lfc_thresh,
      stratum_label   = sprintf("|logFC|>%.2f", lfc_thresh),
      label           = label,
      color           = color,
      pos_color       = pos_color,
      neg_color       = neg_color,
      force_rerun     = key %in% FORCE_RERUN_KEYS
    )

  }  # end threshold loop

  # ── All-peaks run (no |logFC| filter) ────────────────────────────────────
  # Runs on every DAR with a valid logFC — no secondary effect-size filter.
  # Outputs go to lmm/all_peaks/ alongside the threshold subdirectories.
  # Justification: the DAR universe is already significance-filtered at the
  # DESeq2/limma level, so no truly flat peaks are included. Removing the
  # secondary |logFC| filter maximises peak set size, equalises n across
  # comparisons, and avoids threshold-dependent result variation.
  if (RUN_ALL_PEAKS) {
    all_dir <- file.path(lmm_dir, "all_peaks")
    dir.create(all_dir, recursive = TRUE, showWarnings = FALSE)

    coef_file_all <- file.path(all_dir, "lmm_coefficients.txt")
    in_force_all  <- key %in% FORCE_RERUN_KEYS
    if (file.exists(coef_file_all) && !in_force_all) {
      cat(sprintf("  SKIP (done): %-40s [all_peaks]\n", key))
    } else {
      if (in_force_all && file.exists(coef_file_all)) {
        cat(sprintf("  FORCE RERUN: %-40s [all_peaks]\n", key))
        file.remove(coef_file_all)
        rds_all <- file.path(all_dir, "lmm_result.rds")
        if (file.exists(rds_all)) file.remove(rds_all)
      }

      model_data_all <- model_data_full[!is.na(model_data_full$logFC), , drop = FALSE]

      run_lmm_stratum(
        key             = key,
        model_data_cont = model_data_all,
        n_peaks_total   = nrow(model_data_all),
        thresh_dir      = all_dir,
        thresh_tag      = "all_peaks",
        lfc_thresh      = 0,
        stratum_label   = "all DARs (no |logFC| filter)",
        label           = label,
        color           = color,
        pos_color       = pos_color,
        neg_color       = neg_color,
        force_rerun     = in_force_all
      )
    }
  }

  invisible(NULL)
}

# ── run_lmm_stratum: inner worker shared by threshold loop + all-peaks run ────
# Fits LMM + binary GLMM + glmmLasso for one key x stratum combination.
run_lmm_stratum <- function(key, model_data_cont, n_peaks_total,
                             thresh_dir, thresh_tag, lfc_thresh, stratum_label,
                             label, color, pos_color=NULL, neg_color=NULL,
                             force_rerun = FALSE) {

  if (n_peaks_total < 20) {
    cat(sprintf("    SKIP [%s|%s] — fewer than 20 peaks\n", key, thresh_tag))
    return(invisible(NULL))
  }

  keep_tf <- names(partner_peaks)[sapply(names(partner_peaks), function(tf) {
    v <- model_data_cont[[tf]]
    !is.null(v) && var(v, na.rm = TRUE) > 0
  })]
  if (length(keep_tf) == 0) {
    cat(sprintf("    SKIP [%s|%s] — no variable TF columns\n", key, thresh_tag))
    return(invisible(NULL))
  }

  formula_cont     <- as.formula(
    paste0("logFC ~ ", paste(keep_tf, collapse = " + "), " + (1|chr)"))
  lmm_n_sub_actual <- min(LMM_N_SUB, n_peaks_total)

  # ── Is this stratum actually big enough for sample() to subsample? ──────
  # If n_peaks_total <= LMM_N_SUB, sample(n, n, replace=FALSE) returns a
  # permutation of the FULL dataset every time — repeated "runs" are the
  # same model fit to the same data, not independent replicates.
  REAL_SUBSAMPLE <- n_peaks_total > LMM_N_SUB

  coll <- tf_collinearity_flag(model_data_cont, keep_tf)
  if (!is.na(coll$max_abs_r) && coll$max_abs_r > 0.80) {
    cat(sprintf("    WARNING: TF design matrix has |r|>0.80 pairs: %s\n", coll$flagged_pairs))
  }

  cat(sprintf("  [LMM] %-40s  %s  %d peaks  %d TFs  (real_subsample=%s)\n",
              key, stratum_label, n_peaks_total, length(keep_tf), REAL_SUBSAMPLE))

  # ── DAR count summary (new in v8) ────────────────────────────────────
  # Answers "total DARs / #preferential each direction / anything excluded"
  # directly from the exact data the model sees, rather than a manual count
  # that could drift out of sync with what actually got fit.
  #   - n_peaks_total here is nrow(model_data_cont) -- i.e. AFTER whatever
  #     filtering happened upstream of run_lmm_stratum. If that doesn't match
  #     the raw DAR count from your DESeq2/DiffBind output for this
  #     comparison, the difference is exactly "excluded before fitting."
  #   - Direction labels (BOTv-preferential / BOTCv-preferential) are only
  #     meaningful for stems where logFC sign has a known convention. This
  #     script negates logFC for BOTv_vs_BOTCv and BOTv_vs_BOTCv_nc14late
  #     specifically (see NEGATE_LFC_STEMS above) so that, for those two
  #     stems ONLY, positive logFC = BOTCv-preferential (matches
  #     positive_means="TF enriched in BOTCv-open regions" in
  #     LMM_PlotGeneration). For any other stem, "pos"/"neg" are reported
  #     without a BOTv/BOTCv label -- don't assume the same sign convention
  #     applies elsewhere without checking NEGATE_LFC_STEMS for that stem.
  n_pos <- sum(model_data_cont$logFC > 0,  na.rm = TRUE)
  n_neg <- sum(model_data_cont$logFC <= 0, na.rm = TRUE)
  is_negated_stem <- key %in% NEGATE_LFC_STEMS
  dar_summary <- data.frame(
    comparison         = key,
    stratum             = thresh_tag,
    lfc_threshold       = lfc_thresh,
    n_total             = n_peaks_total,
    n_logFC_pos         = n_pos,
    n_logFC_neg         = n_neg,
    n_botcv_preferential = if (is_negated_stem) n_pos else NA_integer_,
    n_botv_preferential  = if (is_negated_stem) n_neg else NA_integer_,
    sign_convention_known = is_negated_stem,
    stringsAsFactors = FALSE)
  write.table(dar_summary, file.path(thresh_dir, "lmm_dar_summary.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  n_singular      <- 0L
  n_nonconverged  <- 0L

  if (!REAL_SUBSAMPLE) {
    # ── SINGLE FIT PATH ─────────────────────────────────────────────────
    # n_peaks_total <= LMM_N_SUB: there is exactly one distinct dataset
    # available, so fit it once and take inference directly from the model
    # (lmerTest / Satterthwaite df) instead of manufacturing 20 duplicate
    # "runs" and Stouffer-combining them (which inflates significance).
    cat(sprintf("    NOTE: n_peaks_total (%d) <= LMM_N_SUB (%d) — fitting once, no resampling.\n",
                n_peaks_total, LMM_N_SUB))

    fit <- tryCatch(
      lmerTest::lmer(formula_cont, data = model_data_cont, REML = FALSE,
                      control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
      error = function(e) NULL)

    if (is.null(fit)) { cat("    Fit failed.\n"); return(invisible(NULL)) }

    conv_msgs <- fit@optinfo$conv$lme4$messages
    if (isSingular(fit))                    n_singular     <- 1L
    if (!is.null(conv_msgs))                n_nonconverged <- 1L

    ct <- as.data.frame(coef(summary(fit)))
    ct$TF <- rownames(ct); rownames(ct) <- NULL
    names(ct) <- gsub("Estimate",          "estimate", names(ct))
    names(ct) <- gsub("Std\\. Error",      "se",       names(ct))
    names(ct) <- gsub("t value",           "t_value",  names(ct))
    names(ct) <- gsub("Pr\\(>\\|t\\|\\)",  "p_value",  names(ct))

    coef_tab <- ct[ct$TF != "(Intercept)", , drop = FALSE]
    coef_tab$n_runs <- 1L
    n_obs_mean <- nrow(model_data_cont)

  } else {
    # ── TRUE BOOTSTRAP PATH ─────────────────────────────────────────────
    # n_peaks_total > LMM_N_SUB: sample() genuinely draws different subsets
    # each run. Still not independent (draws overlap), so we report the
    # empirical between-run spread as a stability diagnostic rather than
    # Stouffer-combining z-scores as though the runs were independent
    # experiments.
    run_results_cont <- vector("list", N_RUNS)
    n_obs_all        <- numeric(N_RUNS)

    for (i in seq_len(N_RUNS)) {
      mdat <- model_data_cont[
        sample(nrow(model_data_cont), lmm_n_sub_actual, replace = FALSE), ]
      n_obs_all[i] <- nrow(mdat)

      fit <- tryCatch(
        lmer(formula_cont, data = mdat, REML = FALSE,
             control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
        error = function(e) NULL)
      if (is.null(fit)) next

      if (isSingular(fit)) n_singular <- n_singular + 1L
      if (!is.null(fit@optinfo$conv$lme4$messages)) n_nonconverged <- n_nonconverged + 1L

      ct <- as.data.frame(coef(summary(fit)))
      ct$TF <- rownames(ct); rownames(ct) <- NULL
      names(ct) <- gsub("Estimate",          "estimate", names(ct))
      names(ct) <- gsub("Std\\. Error",      "se",       names(ct))
      names(ct) <- gsub("t value",           "t_value",  names(ct))
      names(ct) <- gsub("Pr\\(>\\|t\\|\\)",  "p_value",  names(ct))
      run_results_cont[[i]] <- ct
      cat(sprintf("    run %2d/%d: n=%d singular=%s\n",
                  i, N_RUNS, nrow(mdat), isSingular(fit)))
    }

    valid_runs <- Filter(Negate(is.null), run_results_cont)
    if (length(valid_runs) == 0) { cat("    All runs failed.\n"); return(invisible(NULL)) }

    coef_tab <- combine_runs_bootstrap(valid_runs)
    if (is.null(coef_tab) || nrow(coef_tab) == 0) {
      cat("    Aggregation failed.\n"); return(invisible(NULL))
    }
    # Bootstrap path has no per-run p-value we can trust as independent;
    # p_value here is a one-sample t-test of the bootstrap mean against 0,
    # reported for reference/plotting only — treat n_sig_padj as a rough
    # stability screen, not a confirmatory significance test.
    coef_tab$p_value <- 2 * pt(-abs(coef_tab$estimate / coef_tab$se),
                                df = pmax(coef_tab$n_runs - 1, 1))
    n_obs_mean <- mean(n_obs_all[n_obs_all > 0])
  }

  coef_tab$p_adj         <- p.adjust(coef_tab$p_value, method = PADJ_METHOD)
  coef_tab$significant   <- !is.na(coef_tab$p_value) & coef_tab$p_value < 0.05
  coef_tab$sig_adj       <- !is.na(coef_tab$p_adj)   & coef_tab$p_adj   < 0.05
  coef_tab$model_type    <- "continuous_logFC"
  coef_tab$real_subsample<- REAL_SUBSAMPLE
  coef_tab$n_runs_used   <- coef_tab$n_runs[1]
  coef_tab$lfc_threshold <- lfc_thresh
  coef_tab$n_peaks_total <- n_peaks_total

  cat(sprintf("    %d sig (raw) | %d sig (BH) | n_singular=%d | n_nonconverged=%d\n",
              sum(coef_tab$significant, na.rm = TRUE),
              sum(coef_tab$sig_adj,     na.rm = TRUE),
              n_singular, n_nonconverged))

  write.table(coef_tab, file.path(thresh_dir, "lmm_coefficients.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("    Saved:", file.path(thresh_dir, "lmm_coefficients.txt"), "\n")

  # ── Bootstrap stability check (continuous model) ─────────────────────
  run_boot <- RUN_BOOTSTRAP_STABILITY &&
    (is.null(BOOTSTRAP_ONLY_KEYS)   || key       %in% BOOTSTRAP_ONLY_KEYS) &&
    (is.null(BOOTSTRAP_ONLY_THRESH) || thresh_tag %in% BOOTSTRAP_ONLY_THRESH)
  if (run_boot) {
    cat(sprintf("    [bootstrap] running %d resamples (continuous model)...\n", N_BOOTSTRAP))
    boot_tab <- bootstrap_stability_check(formula_cont, model_data_cont,
                                          family_glmer = NULL,
                                          primary_coef = coef_tab)
    if (!is.null(boot_tab)) {
      write.table(boot_tab, file.path(thresh_dir, "lmm_bootstrap_stability.txt"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      cat("    Saved:", file.path(thresh_dir, "lmm_bootstrap_stability.txt"), "\n")
    }
  }

  stats_df <- data.frame(
    comparison        = key,
    label              = label,
    stratum            = thresh_tag,
    lfc_threshold      = lfc_thresh,
    n_peaks_total      = n_peaks_total,
    real_subsample     = REAL_SUBSAMPLE,
    n_sub_per_run      = lmm_n_sub_actual,
    pct_used           = round(100 * lmm_n_sub_actual / n_peaks_total, 1),
    n_runs             = coef_tab$n_runs_used[1],
    n_obs_mean         = n_obs_mean,
    n_tfs              = length(keep_tf),
    n_singular         = n_singular,
    n_nonconverged     = n_nonconverged,
    max_abs_tf_cor     = coll$max_abs_r,
    flagged_tf_pairs   = coll$flagged_pairs,
    n_sig_rawp         = sum(coef_tab$significant, na.rm = TRUE),
    n_sig_padj         = sum(coef_tab$sig_adj,     na.rm = TRUE),
    stringsAsFactors = FALSE)
  write.table(stats_df, file.path(thresh_dir, "lmm_model_stats.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # ── Binary Opening vs Closing model ─────────────────────────────────
  # Same REAL_SUBSAMPLE branch as the continuous model above, for the same
  # reason: with n_peaks_total <= LMM_N_SUB there is only one distinct
  # dataset, so fit once rather than manufacturing duplicate "runs."
  if (RUN_BINARY && n_peaks_total >= 50) {
    bin_data        <- model_data_cont
    bin_data$binary <- as.integer(bin_data$logFC > 0)
    if (length(unique(bin_data$binary)) == 2 &&
        sum( bin_data$binary) >= 10 &&
        sum(!bin_data$binary) >= 10) {

      formula_bin <- as.formula(
        paste0("binary ~ ", paste(keep_tf, collapse = " + "), " + (1|chr)"))

      if (!REAL_SUBSAMPLE) {
        fit_b <- tryCatch(
          glmer(formula_bin, data = bin_data, family = binomial("logit"),
                control = glmerControl(optimizer = "bobyqa",
                                        optCtrl = list(maxfun = 2e5))),
          error = function(e) NULL)
        coef_bin <- NULL
        if (!is.null(fit_b)) {
          ct <- as.data.frame(coef(summary(fit_b)))
          ct$TF <- rownames(ct); rownames(ct) <- NULL
          names(ct) <- gsub("Estimate",          "estimate", names(ct))
          names(ct) <- gsub("Std\\. Error",      "se",       names(ct))
          names(ct) <- gsub("z value",           "z_value",  names(ct))
          names(ct) <- gsub("Pr\\(>\\|z\\|\\)",  "p_value",  names(ct))
          coef_bin <- ct[ct$TF != "(Intercept)", , drop = FALSE]
          coef_bin$n_runs <- 1L
        }
        n_bin_runs <- 1L
      } else {
        run_results_bin <- vector("list", N_RUNS)
        for (i in seq_len(N_RUNS)) {
          mdat_b <- bin_data[
            sample(nrow(bin_data), lmm_n_sub_actual, replace = FALSE), ]
          fit_b <- tryCatch(
            glmer(formula_bin, data = mdat_b, family = binomial("logit"),
                  control = glmerControl(optimizer = "bobyqa",
                                          optCtrl = list(maxfun = 2e5))),
            error = function(e) NULL)
          if (is.null(fit_b)) next
          ct <- as.data.frame(coef(summary(fit_b)))
          ct$TF <- rownames(ct); rownames(ct) <- NULL
          names(ct) <- gsub("Estimate",          "estimate", names(ct))
          names(ct) <- gsub("Std\\. Error",      "se",       names(ct))
          names(ct) <- gsub("z value",           "z_value",  names(ct))
          names(ct) <- gsub("Pr\\(>\\|z\\|\\)",  "p_value",  names(ct))
          run_results_bin[[i]] <- ct
        }
        valid_bin <- Filter(Negate(is.null), run_results_bin)
        coef_bin  <- if (length(valid_bin) > 0) combine_runs_bootstrap(valid_bin) else NULL
        if (!is.null(coef_bin) && nrow(coef_bin) > 0) {
          coef_bin$p_value <- 2 * pt(-abs(coef_bin$estimate / coef_bin$se),
                                      df = pmax(coef_bin$n_runs - 1, 1))
        }
        n_bin_runs <- length(valid_bin)
      }

      if (!is.null(coef_bin) && nrow(coef_bin) > 0) {
        coef_bin$p_adj         <- p.adjust(coef_bin$p_value, method = PADJ_METHOD)
        coef_bin$significant   <- !is.na(coef_bin$p_value) & coef_bin$p_value < 0.05
        coef_bin$sig_adj       <- !is.na(coef_bin$p_adj)   & coef_bin$p_adj   < 0.05
        coef_bin$model_type    <- "binary_Opening_vs_Closing"
        coef_bin$real_subsample<- REAL_SUBSAMPLE
        coef_bin$n_runs_used   <- coef_bin$n_runs[1]
        coef_bin$lfc_threshold <- lfc_thresh
        coef_bin$n_peaks_total <- n_peaks_total
        write.table(coef_bin,
                    file.path(thresh_dir, "lmm_coefficients_binary.txt"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
        cat(sprintf("    Binary model: %d sig (raw) | %d sig (BH) | real_subsample=%s\n",
                    sum(coef_bin$significant, na.rm = TRUE),
                    sum(coef_bin$sig_adj,     na.rm = TRUE), REAL_SUBSAMPLE))

        # ── Bootstrap stability check (binary model) ─────────────────────
        if (run_boot) {
          cat(sprintf("    [bootstrap] running %d resamples (binary model)...\n", N_BOOTSTRAP))
          boot_bin <- bootstrap_stability_check(formula_bin, bin_data,
                                                family_glmer = binomial("logit"),
                                                primary_coef = coef_bin)
          if (!is.null(boot_bin)) {
            write.table(boot_bin, file.path(thresh_dir, "lmm_bootstrap_stability_binary.txt"),
                        sep = "\t", quote = FALSE, row.names = FALSE)
            cat("    Saved:", file.path(thresh_dir, "lmm_bootstrap_stability_binary.txt"), "\n")
          }
        }

        save_coef_plot(coef_bin, thresh_dir, "lmm_binary_Opening_vs_Closing.pdf",
                       sprintf("%s [binary Opening/Closing]", label),
                       color, n_bin_runs, n_obs_mean, length(keep_tf),
                       is_binary = TRUE, n_peaks_total = n_peaks_total,
                       pos_color = pos_color, neg_color = neg_color,
                       real_subsample = REAL_SUBSAMPLE)
      }
    }
  }

  # ── Save RDS + continuous plot ───────────────────────────────────────
  saveRDS(list(coefficients   = coef_tab,
               stats          = stats_df,
               lfc_threshold  = lfc_thresh,
               stratum        = thresh_tag,
               real_subsample = REAL_SUBSAMPLE,
               n_runs         = coef_tab$n_runs_used[1],
               formula        = deparse(formula_cont)),
          file.path(thresh_dir, "lmm_result.rds"))

  save_coef_plot(coef_tab, thresh_dir,
                 sprintf("lmm_coefficients_rawp_%s.pdf", thresh_tag),
                 sprintf("%s  [%s]", label, stratum_label),
                 color, coef_tab$n_runs_used[1], n_obs_mean, length(keep_tf),
                 is_binary = FALSE, n_peaks_total = n_peaks_total,
                 pos_color = pos_color, neg_color = neg_color,
                 real_subsample = REAL_SUBSAMPLE)

  # ── glmmLasso ────────────────────────────────────────────────────────
  run_glmmlasso(
    key             = key,
    model_data_cont = model_data_cont,
    thresh_dir      = thresh_dir,
    label           = sprintf("%s  [%s]", label, stratum_label),
    color           = color,
    lfc_thresh      = lfc_thresh,
    keep_tf         = keep_tf,
    force_rerun     = force_rerun
  )

  cat(sprintf("    Done: %s\n\n", thresh_dir))
  invisible(NULL)
}


################################################################################
# COMPARISON CONFIGS  (identical to v1 — paths, stems, colors)
################################################################################

split1_configs <- list(
  BOT_vs_BOTR               = list(stem="BOT_vs_BOTR",                    label="Runt WT vs null (BOT vs BOTR)",             color="#1f77b4", lmm_dir=file.path(split1_output,"results_BOT_vs_BOTR","lmm")),
  BOT_vs_BOT_hR             = list(stem="BOT_vs_BOT_hR",                  label="Runt WT vs het (BOT vs BOT_hR)",            color="#17becf", lmm_dir=file.path(split1_output,"results_BOT_vs_BOT_hR","lmm")),
  BOT_hR_vs_BOTR            = list(stem="BOT_hR_vs_BOTR",                 label="Runt het vs null (BOT_hR vs BOTR)",         color="#17becf", lmm_dir=file.path(split1_output,"results_BOT_hR_vs_BOTR","lmm")),
  Runt_WT_dosage_gradient   = list(stem="Runt_WT_dosage_gradient_non_vent",label="Runt WT dosage gradient",                  color="#1f77b4", lmm_dir=file.path(split1_output,"results_Runt_WT_dosage_gradient","lmm")),
  BOT_vs_BOTC               = list(stem="BOT_vs_BOTC",                    label="Cic effect, Runt intact (BOT vs BOTC)",     color="#1f77b4", lmm_dir=file.path(split1_output,"results_BOT_vs_BOTC","lmm")),
  BOTC_vs_BOTR              = list(stem="BOTC_vs_BOTR",                   label="Cic-del vs Runt-null (BOTC vs BOTR)",       color="#ff7f0e", lmm_dir=file.path(split1_output,"results_BOTC_vs_BOTR","lmm")),
  BOTC_oR_vs_BOTC           = list(stem="BOTC_oR_vs_BOTC",               label="Runt rescue in Cic-del",                   color="#bcbd22", lmm_dir=file.path(split1_output,"results_BOTC_oR_vs_BOTC","lmm")),
  BOTC_oR_vs_BOTR           = list(stem="BOTC_oR_vs_BOTR",               label="oRunt vs Runt-null (BOTC_oR vs BOTR)",     color="#bcbd22", lmm_dir=file.path(split1_output,"results_BOTC_oR_vs_BOTR","lmm")),
  Runt_rescue_cross_context = list(stem="Runt_rescue_cross_Cic_context",  label="Runt rescue cross-Cic-context",            color="#bcbd22", lmm_dir=file.path(split1_output,"results_Runt_rescue_cross_context","lmm"))
)

split2_configs <- list(
  # pos_color = BOTCv/pink (Cic absent opens), neg_color = BOTv/green (Cic present opens)
  BOTv_vs_BOTCv        = list(stem="BOTv_vs_BOTCv",                    label="BOTv vs BOTCv (nc14b)",              color="#e377c2", pos_color="#e377c2", neg_color="#2ca02c", lmm_dir=file.path(split2_output,"results_BOTv_vs_BOTCv","lmm")),
  BOTv_vs_BOTCv_late   = list(stem="BOTv_vs_BOTCv_nc14late",           label="BOTv vs BOTCv (nc14late)",           color="#b5369a", pos_color="#b5369a", neg_color="#1a6b1a", lmm_dir=file.path(split2_output,"results_BOTv_vs_BOTCv_late","lmm")),
  tolrm9_inferred      = list(stem="tolrm9_inferred_bg",               label="tolrm9 inferred background",        color="#BDC3C7", lmm_dir=file.path(split2_output,"results_tolrm9_inferred","lmm")),
  tolrm9_empirical     = list(stem="tolrm9_empirical_bg",              label="tolrm9 empirical background",       color="#95A5A6", lmm_dir=file.path(split2_output,"results_tolrm9_empirical","lmm")),
  tolrm9_validated     = list(stem="tolrm9_validated_both_methods",    label="tolrm9 validated",                  color="#7F8C8D", lmm_dir=file.path(split2_output,"results_tolrm9_validated","lmm")),
  BOTCv_above_inf      = list(stem="BOTCv_above_inferred_bg",          label="BOTCv above inferred tolrm9",       color="#e377c2", lmm_dir=file.path(split2_output,"results_BOTCv_above_inf","lmm")),
  BOTCv_above_emp      = list(stem="BOTCv_above_empirical_bg",         label="BOTCv above empirical tolrm9",      color="#e377c2", lmm_dir=file.path(split2_output,"results_BOTCv_above_emp","lmm")),
  BOTv_above_inf       = list(stem="BOTv_above_inferred_bg",           label="BOTv above inferred tolrm9",        color="#2ca02c", lmm_dir=file.path(split2_output,"results_BOTv_above_inf","lmm")),
  BOTCv_only           = list(stem="BOTCv_only_crossvalidated",        label="BOTCv-specific (cross-validated)",  color="#e377c2", lmm_dir=file.path(split2_output,"results_BOTCv_only","lmm")),
  BOTv_only            = list(stem="BOTv_only_crossvalidated",         label="BOTv-specific (cross-validated)",   color="#2ca02c", lmm_dir=file.path(split2_output,"results_BOTv_only","lmm")),
  BOTv_BOTCv_div       = list(stem="BOTv_BOTCv_divergent_vs_BOTR",     label="BOTv vs BOTCv divergent",           color="#7F8C8D", lmm_dir=file.path(split2_output,"results_BOTv_BOTCv_div","lmm")),
  BOTv_temporal        = list(stem="BOTv_temporal",                    label="BOTv temporal (nc14b→nc14late)",    color="#1a6b1a", lmm_dir=file.path(split2_output,"results_BOTv_temporal","lmm")),
  BOTCv_temporal       = list(stem="BOTCv_temporal",                   label="BOTCv temporal (nc14b→nc14late)",   color="#b5369a", lmm_dir=file.path(split2_output,"results_BOTCv_temporal","lmm")),
  stable               = list(stem="BOTv_BOTCv_stable_both_timepoints",label="Stable (both timepoints)",          color="#117A65", lmm_dir=file.path(split2_output,"results_stable","lmm")),
  emerging             = list(stem="BOTv_BOTCv_emerging_late_only",    label="Emerging (late only)",              color="#F39C12", lmm_dir=file.path(split2_output,"results_emerging","lmm")),
  resolving            = list(stem="BOTv_BOTCv_resolving_early_only",  label="Resolving (early only)",            color="#8E44AD", lmm_dir=file.path(split2_output,"results_resolving","lmm")),
  temporal_divergent   = list(stem="temporal_divergent_BOTv_vs_BOTCv", label="Temporal divergent",                color="#E74C3C", lmm_dir=file.path(split2_output,"results_temporal_divergent","lmm"))
)

################################################################################
# RUN
################################################################################

run_all <- function(configs, split_label, base_dir) {
  cat("\n", strrep("=",70), "\n")
  cat(split_label, "\n")
  cat(strrep("=",70), "\n\n")
  for (key in names(configs)) {
    cfg <- configs[[key]]
    tryCatch({
      gr <- load_dar_bed(cfg$stem, base_dir)
      run_and_save_lmm(key=key, gr=gr, lmm_dir=cfg$lmm_dir,
                       label=cfg$label, color=cfg$color,
                       pos_color=cfg$pos_color %||% NULL,
                       neg_color=cfg$neg_color %||% NULL)
    }, error=function(e) {
      cat(sprintf("  ERROR in %s: %s\n", key, conditionMessage(e)))
    })
    gc(verbose=FALSE)
  }
}

cat("[4/5] Running Split 2 (vent — runs first)...\n")
run_all(split2_configs, "SPLIT 2 — VENTRALIZED", temp_base_dir)

cat("[5/5] Running Split 1 (non-vent)...\n")
run_all(split1_configs, "SPLIT 1 — NON-VENTRALIZED", dar_base_dir)

# Helper to rerun specific keys only (e.g. after sign fix)
rerun_missing <- function(keys_to_rerun) {
  all_configs <- c(split1_configs, split2_configs)
  base_dirs   <- c(setNames(rep(dar_base_dir,  length(split1_configs)), names(split1_configs)),
                   setNames(rep(temp_base_dir, length(split2_configs)), names(split2_configs)))
  for (key in keys_to_rerun) {
    cfg <- all_configs[[key]]
    if (is.null(cfg)) { cat(sprintf("  UNKNOWN: %s\n", key)); next }
    coef_file <- file.path(cfg$lmm_dir, "lmm_coefficients.txt")
    if (file.exists(coef_file)) file.remove(coef_file)
    rds <- file.path(cfg$lmm_dir, "lmm_result.rds")
    if (file.exists(rds)) file.remove(rds)
    tryCatch({
      gr <- load_dar_bed(cfg$stem, base_dirs[[key]])
      run_and_save_lmm(key=key, gr=gr, lmm_dir=cfg$lmm_dir,
                       label=cfg$label, color=cfg$color)
    }, error=function(e) cat(sprintf("  ERROR in %s: %s\n", key, e$message)))
    gc(verbose=FALSE)
  }
}

################################################################################
# SUMMARY
################################################################################

cat("\n", strrep("=",70), "\n")
cat("COMPLETE\n")
cat(strrep("=",70), "\n\n")
all_configs <- c(split1_configs, split2_configs)
n_ok <- 0; n_lasso_ok <- 0
for (key in names(all_configs)) {
  cfg <- all_configs[[key]]
  thresh_statuses <- sapply(LFC_THRESHOLDS, function(thr) {
    td          <- file.path(cfg$lmm_dir, sprintf("lfc%.2f", thr))
    has_coef    <- file.exists(file.path(td, "lmm_coefficients.txt"))
    has_rds     <- file.exists(file.path(td, "lmm_result.rds"))
    has_bin     <- file.exists(file.path(td, "lmm_coefficients_binary.txt"))
    has_lasso_c <- file.exists(file.path(td, "glmmlasso_continuous_coefs.txt"))
    has_lasso_b <- file.exists(file.path(td, "glmmlasso_binary_coefs.txt"))
    lasso_tag   <- if (has_lasso_c || has_lasso_b)
                     paste0("+lasso",
                            if (has_lasso_c) "c" else "",
                            if (has_lasso_b) "b" else "")
                   else ""
    if (has_coef && has_rds) {
      n_ok <<- n_ok + 1
      if (has_lasso_c || has_lasso_b) n_lasso_ok <<- n_lasso_ok + 1
      paste0("OK", if(has_bin) "+bin" else "", lasso_tag)
    } else if (has_coef) "no-rds"
    else "MISS"
  })
  cat(sprintf("  %-35s  %s\n", key,
              paste(sprintf("|lfc|>%.2f:%s", LFC_THRESHOLDS, thresh_statuses), collapse="  ")))
}
cat(sprintf("\n%d / %d threshold-comparison combinations complete\n", n_ok,
            length(all_configs) * length(LFC_THRESHOLDS)))
cat(sprintf("%d with glmmLasso output\n\n", n_lasso_ok))
cat("Output files per comparison x threshold (e.g. lfc0.50/):\n")
cat("  [Plain LMM / binary GLMM — unchanged from v3]\n")
cat("  lmm/lfc0.50/lmm_coefficients.txt                 — continuous LMM (averaged over N_RUNS subsampling runs)\n")
cat("  lmm/lfc0.50/lmm_coefficients_binary.txt          — binary GLMM Opening/Closing\n")
cat("  lmm/lfc0.50/lmm_result.rds\n")
cat("  lmm/lfc0.50/lmm_coefficients_rawp_lfc0.50.pdf\n")
cat("  lmm/lfc0.50/lmm_binary_Opening_vs_Closing.pdf\n")
cat("\n  [Regularized GLMM-LASSO — new in v4]\n")
cat("  lmm/lfc0.50/glmmlasso_continuous_coefs.txt       — LASSO continuous model (lambda by 5-fold CV)\n")
cat("  lmm/lfc0.50/glmmlasso_binary_coefs.txt           — LASSO binary Opening/Closing model\n")
cat("  lmm/lfc0.50/glmmlasso_continuous_lfc0.50.pdf\n")
cat("  lmm/lfc0.50/glmmlasso_binary_lfc0.50.pdf\n\n")
cat("To rerun specific keys after sign fix:\n")
cat("  rerun_missing(c('BOTv_vs_BOTCv', 'BOTv_vs_BOTCv_late'))\n\n")
cat("LASSO settings used:\n")
cat(sprintf("  Lambda grid: seq(%.2f, %.2f, length.out=%d)  [linear, matches original]\n",
            LASSO_LAMBDA_END, LASSO_LAMBDA_START, LASSO_N_LAMBDA))
cat(sprintf("  Lambda selection: BIC  |  Subsample n: %d  |  Effect threshold: %.3f\n\n",
            LASSO_N_SUB, LASSO_EFFECT_THRESH))

