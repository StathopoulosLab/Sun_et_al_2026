#!/usr/bin/env Rscript
# =============================================================================
# deg_limma_botv_botcv.R
#
# Limma-voom differential expression for FoxL1_BOTv vs HLH54F_BOTCv.
# Reads the same count matrix and sample annotations used by the Python
# pipeline, runs voom+eBayes across all five contrasts, and writes results
# into the same folder structure as qc_and_deg_botv_botcv.py so all DEG
# outputs live together.
#
# WHY LIMMA-VOOM INSTEAD OF DESEQ2 HERE
# --------------------------------------
# With n=2-3 replicates per arm, DESeq2's negative-binomial dispersion
# estimation hits a ceiling: it can't distinguish real gene-level variance
# from noise, so it assigns near-maximum dispersion to almost everything and
# almost nothing survives FDR correction (you'll see lfcSE ~1.2-1.3 for
# every gene and padj ~0.96 uniformly — that's the tell).
#
# Limma-voom handles this better because:
#   1. voom converts counts to log2-CPM with precision weights that
#      account for the mean-variance relationship, working in log space
#      where the linear model assumptions hold better at low n.
#   2. eBayes borrows variance information across ALL genes (empirical
#      Bayes shrinkage of gene-wise variances toward a common prior),
#      which is exactly what gives limma power when per-gene replication
#      is low — each gene "borrows" evidence from thousands of others.
#   3. RUVg-corrected counts are already on a scale where voom's
#      mean-variance modelling works well.
#
# This is the same approach used in your ATAC-seq DAR pipeline (limma-voom
# on RUVg-corrected data) so the statistical framework is consistent.
#
# CONTRASTS
# ---------
#  1. BOTCv_nc14b  vs  BOTv_nc14b        early genotype contrast
#  2. BOTCv_late   vs  BOTv_nc14d        late genotype contrast
#                   (BOTCv nc14d+gastr pooled  vs  BOTv nc14d)
#  3. BOTCv_late   vs  BOTCv_nc14b       BOTCv temporal
#  4. BOTv_nc14d   vs  BOTv_nc14b        BOTv temporal
#  5. BOTCv_late   vs  BOTv_nc14b        cross-genotype cross-time
#
# USAGE
# -----
#   Rscript deg_limma_botv_botcv.R
#   Rscript deg_limma_botv_botcv.R --fdr 0.10
#   Rscript deg_limma_botv_botcv.R --lfc 0.5 --fdr 0.10
#   Rscript deg_limma_botv_botcv.R --counts path/to/counts.tsv
#   Rscript deg_limma_botv_botcv.R --contrasts 1,3   # run only contrasts 1 and 3
#
# OUTPUT
# ------
#   results/combined_figures/qc_and_deg_botv_botcv/deg_limma/
#     <contrast_name>/results_full.tsv    all genes, sorted by adj.P.Val
#     <contrast_name>/results_sig.tsv     significant hits
#     <contrast_name>/volcano.pdf/png     volcano plot
#     summary_table.tsv                   one row per contrast
#
# REQUIRED R PACKAGES
# -------------------
#   limma, edgeR  — install with: BiocManager::install(c("limma","edgeR"))
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("limma", quietly=TRUE))
    stop("limma not found. Install with: BiocManager::install('limma')")
  if (!requireNamespace("edgeR", quietly=TRUE))
    stop("edgeR not found. Install with: BiocManager::install('edgeR')")
  library(limma)
  library(edgeR)
})

# =============================================================================
# CONFIG — change here or via CLI flags
# =============================================================================

COUNTS_PATH   <- "results/combined/all_genotypes_plus_set3_counts_matrix_qc.tsv"
SYMBOL_MAP    <- "gene_id_to_symbol.tsv"
OUT_DIR       <- "results/combined_figures/qc_and_deg_botv_botcv/deg_limma"

GROUP_A       <- "FoxL1_BOTv"
GROUP_B       <- "HLH54F_BOTCv"

# BOTCv "late" = nc14d + gastr pooled (low replicates individually)
# BOTv "late"  = nc14d only (has n>=2 there)
BOTCV_LATE_STAGES <- c("nc14d", "gastr")
BOTV_LATE_STAGE   <- "nc14d"

# Embryos to exclude (barcode substring matched against column names)
EXCLUDE_BCS   <- c("bc27")   # FoxL1_BOTv bc27 — same as Python pipeline

FDR_THRESH    <- 0.05
LFC_THRESH    <- 0.5   # log2(1.5) — same default as Python script
MIN_CPM       <- 1       # filter: gene must have CPM >= MIN_CPM in at
                         # least MIN_SAMPLES samples in the contrast.
TOP_LABEL_N   <- 15      # max labelled points on volcano

# =============================================================================
# CLI argument parsing
# =============================================================================

args <- commandArgs(trailingOnly=TRUE)
parse_arg <- function(flag, default) {
  i <- which(args == flag)
  if (length(i) && length(args) > i) return(args[i+1])
  return(default)
}
COUNTS_PATH <- parse_arg("--counts",   COUNTS_PATH)
OUT_DIR     <- parse_arg("--outdir",   OUT_DIR)
FDR_THRESH  <- as.numeric(parse_arg("--fdr", FDR_THRESH))
LFC_THRESH  <- as.numeric(parse_arg("--lfc", LFC_THRESH))
MIN_CPM     <- as.numeric(parse_arg("--min-cpm", MIN_CPM))

wanted_contrasts <- parse_arg("--contrasts", "1,2,3,4,5")
wanted_contrasts <- as.integer(strsplit(wanted_contrasts, ",")[[1]])

PRINT_COLS <- "--print-cols" %in% args   # diagnostic: dump all column names then exit

cat(sprintf("[limma] FDR=%.2f  |  log2FC threshold=%.3f\n", FDR_THRESH, LFC_THRESH))

# =============================================================================
# Helpers
# =============================================================================

get_stage <- function(colname) {
  # Extract stage from canonical column name, covering both naming conventions:
  #   Long form:  "FoxL1_BOTv_embryo_bc1"    -> nc14b (no stage tag = early set)
  #   Short form: "foxl1_Dm_nc14d_bc25"       -> nc14d
  #               "HLH54F_Dm_gastr_bc39"      -> gastr
  #               "HLH54F_BOTCv_Dm_nc14b_bc32"-> nc14b
  if (grepl("gastr",        colname, ignore.case=TRUE)) return("gastr")
  if (grepl("nc14d|nc14late|nc14_late", colname, ignore.case=TRUE)) return("nc14d")
  if (grepl("nc14b",        colname, ignore.case=TRUE)) return("nc14b")
  # "nc14" without b/d suffix (e.g. foxl1_Dm_nc14_bc20) — treat as nc14b
  if (grepl("nc14",         colname, ignore.case=TRUE)) return("nc14b")
  # No stage tag at all (set1/set2 embryos: FoxL1_BOTv_embryo_bc1) -> nc14b
  return("nc14b")
}

get_group <- function(colname) {
  # Match both the full data_loader canonical form AND the raw counts matrix form:
  #   GROUP_A (FoxL1_BOTv):   "FoxL1_BOTv_embryo_bc1"  OR  "foxl1_Dm_nc14d_bc25"
  #   GROUP_B (HLH54F_BOTCv): "HLH54F_BOTCv_Dm_nc14b_bc32"  OR  "HLH54F_Dm_nc14b_bc32"
  # Using case-insensitive matching and alternative patterns for each group.

  # Sublibrary 28273 used a single shared "foxl1_Dm_" demux prefix for BOTH
  # FoxL1_BOTv (bc25-30) AND B6_BOTC (bc19-24) — the prefix alone cannot
  # distinguish them (see Config_SampleMetadata.py's SET2_BARCODE_MAP / the note in
  # VALIDATED_FOXL1: "foxl1_Dm_nc14_bc20-23 are B6_BOTC, not FoxL1_BOTv").
  # Without this guard, bc19-24 would be wrongly matched into GROUP_A below
  # via the "foxl1_Dm" pattern, contaminating the FoxL1_BOTv arm with a
  # third, unrelated, non-ventralized genotype.
  if (grepl("^foxl1_Dm", colname, ignore.case=TRUE)) {
    m <- regmatches(colname, regexpr("bc([0-9]+)$", colname))
    if (length(m)) {
      bc_num <- suppressWarnings(as.integer(sub("bc", "", m)))
      if (!is.na(bc_num) && bc_num >= 19 && bc_num <= 24) {
        return(NA_character_)   # B6_BOTC — out of scope for this script
      }
    }
  }

  # "BOTv" matches set3's bare naming (BOTv_nc14b_bcN). Safe to add: it is
  # NOT a substring of "HLH54F_BOTCv" (that's "BOTCv" — note the extra C —
  # so no collision with GROUP_B).
  GROUP_A_PATTERNS <- c("FoxL1_BOTv", "foxl1_Dm", "FoxL1_Dm", "BOTv")
  GROUP_B_PATTERNS <- c("HLH54F_BOTCv", "HLH54F_Dm")

  for (pat in GROUP_A_PATTERNS)
    if (grepl(pat, colname, ignore.case=TRUE)) return(GROUP_A)
  for (pat in GROUP_B_PATTERNS)
    if (grepl(pat, colname, ignore.case=TRUE)) return(GROUP_B)
  return(NA_character_)
}

short_name <- function(colname) {
  # Return the barcode portion for readable labels
  m <- regmatches(colname, regexpr("bc[0-9]+", colname))
  if (length(m)) return(m)
  return(substr(colname, nchar(colname)-7, nchar(colname)))
}

save_plot <- function(stem, expr, env=parent.frame(), width=7, height=6) {
  # eval in the caller's environment so local variables (res, def, etc.) are visible
  pdf(paste0(stem, ".pdf"), width=width, height=height)
  tryCatch(eval(expr, envir=env), finally=dev.off())
  png(paste0(stem, ".png"), width=width*100, height=height*100, res=100)
  tryCatch(eval(expr, envir=env), finally=dev.off())
}

# =============================================================================
# Load data
# =============================================================================

cat("[limma] Loading counts matrix:", COUNTS_PATH, "\n")
if (!file.exists(COUNTS_PATH))
  stop(sprintf("Counts matrix not found: %s", COUNTS_PATH))

# quote="" is critical: gene names / row names sometimes contain quote
# characters that trip up R's default parser, causing truncated reads.
counts_raw <- read.table(COUNTS_PATH, sep="\t", header=TRUE,
                          row.names=1, check.names=FALSE, quote="",
                          comment.char="")
if ("GeneSymbol" %in% colnames(counts_raw))
  counts_raw <- counts_raw[, colnames(counts_raw) != "GeneSymbol"]

cat(sprintf("[limma] Raw counts: %d genes x %d embryos\n",
            nrow(counts_raw), ncol(counts_raw)))

# ── Column name normalisation ─────────────────────────────────────────────
# The counts TSV may have long pipeline column names like:
#   "28232_28233_28271_28272_28273_S1__HLH54F_BOTCv_Dm_nc14b_bc32"
# or short canonical names like:
#   "HLH54F_BOTCv_Dm_nc14b_bc32"
# Python data_loader strips everything up to and including "__".
# We do the same here, then print a sample so issues are visible.

strip_prefix <- function(nm) {
  if (grepl("__", nm, fixed=TRUE)) sub(".*__", "", nm) else nm
}

# Capture real batch identity BEFORE stripping the sample prefix, while the
# raw "sample__embryo" names (e.g. "28232_28233_..._S1__FoxL1_BOTv_embryo_bc1")
# are still intact — only the prefix carries the true sequencing-batch ID.
# Bare, already-stripped embryo names carry no batch information, so
# get_batch() must run on the raw names, not the post-strip ones used
# everywhere else in this script.
get_batch <- function(nm) {
  if (grepl("__", nm, fixed=TRUE)) return(sub("__.*", "", nm))
  return("unknown_batch")   # no sample prefix present — shouldn't normally occur
}
raw_cols_pre_strip <- colnames(counts_raw)
bare_names_pre_strip <- sapply(raw_cols_pre_strip, strip_prefix, USE.NAMES=FALSE)
batch_of_bare <- setNames(sapply(raw_cols_pre_strip, get_batch, USE.NAMES=FALSE),
                          bare_names_pre_strip)

colnames(counts_raw) <- sapply(colnames(counts_raw), strip_prefix,
                                USE.NAMES=FALSE)

cat("[limma] Column names after prefix-strip (first 6):\n")
cat(paste(" ", head(colnames(counts_raw), 6), collapse="\n"), "\n")

# Build sample annotation
all_cols  <- colnames(counts_raw)
group_vec <- sapply(all_cols, get_group, USE.NAMES=FALSE)
stage_vec <- sapply(all_cols, get_stage, USE.NAMES=FALSE)

cat(sprintf("[limma] GROUP_A ('%s') matches %d columns (via get_group(), not literal string match)\n",
            GROUP_A, sum(group_vec == GROUP_A, na.rm=TRUE)))
cat(sprintf("[limma] GROUP_B ('%s') matches %d columns (via get_group(), not literal string match)\n",
            GROUP_B, sum(group_vec == GROUP_B, na.rm=TRUE)))

# Keep only our two groups, drop excluded barcodes
keep_cols <- !is.na(group_vec) &
             !sapply(all_cols, function(c)
               any(sapply(EXCLUDE_BCS, function(bc) grepl(bc, c, fixed=TRUE))))
counts <- counts_raw[, keep_cols, drop=FALSE]
group  <- group_vec[keep_cols]
stage  <- stage_vec[keep_cols]
cols   <- colnames(counts)

cat(sprintf("[limma] %d embryos after exclusions\n", ncol(counts)))
for (g in c(GROUP_A, GROUP_B)) {
  idx <- group == g
  cat(sprintf("[limma]   %s: %d  (%s)\n", g, sum(idx),
              paste(sprintf("%s(%s)", short_name(cols[idx]), stage[idx]),
                    collapse=", ")))
}
if (ncol(counts) == 0 || PRINT_COLS) {
  cat("\n[limma] All column names in counts file (after prefix strip):\n")
  cat(paste(sprintf("  [%02d] %s  |  group=%s  stage=%s",
                    seq_along(all_cols), all_cols, group_vec, stage_vec),
            collapse="\n"), "\n")
  if (ncol(counts) == 0)
    stop("[limma] No embryos found for GROUP_A/GROUP_B. ",
         "Set GROUP_A / GROUP_B in the script to match the column names above,\n",
         "       or run:  Rscript deg_limma_botv_botcv.R --print-cols")
  else {
    cat("\n[limma] --print-cols done. Exiting.\n"); quit(save="no")
  }
}

# Symbol map (optional)
fbgn_to_sym <- NULL
if (file.exists(SYMBOL_MAP)) {
  sym_df <- read.table(SYMBOL_MAP, sep="\t", header=TRUE,
                       stringsAsFactors=FALSE, quote="", comment.char="")
  fbgn_to_sym <- setNames(sym_df$Symbol, sym_df$GeneID)
  cat(sprintf("[limma] Loaded %d gene symbol mappings\n", nrow(sym_df)))
}
sym <- function(gene_id) {
  if (is.null(fbgn_to_sym)) return(gene_id)
  s <- fbgn_to_sym[gene_id]
  ifelse(is.na(s) | s == "", gene_id, s)
}

# =============================================================================
# Contrast definitions
# =============================================================================
# Each: list(id, name, cols_A, cols_B, description)
# positive logFC = higher in A

make_cols <- function(grp, stg) {
  cols[group == grp & stage %in% stg]
}

BOTCV_LATE <- make_cols(GROUP_B, BOTCV_LATE_STAGES)
BOTV_LATE  <- make_cols(GROUP_A, BOTV_LATE_STAGE)
BOTCV_EARLY<- make_cols(GROUP_B, "nc14b")
BOTV_EARLY <- make_cols(GROUP_A, "nc14b")

contrast_defs <- list(
  list(id=1, name="BOTCv_nc14b_vs_BOTv_nc14b",
       A=BOTCV_EARLY, B=BOTV_EARLY,
       desc="BOTCv early vs BOTv early — same stage, different genotype"),
  list(id=2, name="BOTCv_late_vs_BOTv_nc14d",
       A=BOTCV_LATE,  B=BOTV_LATE,
       desc=sprintf("BOTCv late (%s pooled) vs BOTv %s — late genotype contrast",
                    paste(BOTCV_LATE_STAGES, collapse="+"), BOTV_LATE_STAGE)),
  list(id=3, name="BOTCv_late_vs_BOTCv_nc14b",
       A=BOTCV_LATE,  B=BOTCV_EARLY,
       desc="BOTCv late vs BOTCv early — temporal within BOTCv"),
  list(id=4, name="BOTv_nc14d_vs_BOTv_nc14b",
       A=BOTV_LATE,   B=BOTV_EARLY,
       desc=sprintf("BOTv %s vs BOTv early — temporal within BOTv", BOTV_LATE_STAGE)),
  list(id=5, name="BOTCv_late_vs_BOTv_nc14b",
       A=BOTCV_LATE,  B=BOTV_EARLY,
       desc="BOTCv late vs BOTv early — cross-genotype cross-time")
)

# =============================================================================
# Run one limma-voom contrast
# =============================================================================

run_contrast <- function(def, counts_all, out_base) {
  A <- def$A;  B <- def$B
  if (length(A) == 0 || length(B) == 0) {
    cat(sprintf("[limma]   SKIP %s: no samples for arm A (%d) or B (%d)\n",
                def$name, length(A), length(B)))
    return(NULL)
  }
  cat(sprintf("[limma]   Arm A (%d): %s\n", length(A),
              paste(short_name(A), collapse=", ")))
  cat(sprintf("[limma]   Arm B (%d): %s\n", length(B),
              paste(short_name(B), collapse=", ")))

  # Subset counts
  use_cols <- c(A, B)
  cnt <- counts_all[, use_cols]
  cnt <- round(cnt)   # ensure integer

  # Condition vector
  condition <- factor(c(rep("A", length(A)), rep("B", length(B))),
                      levels=c("B", "A"))   # B=reference so A-B is positive

  # ── Sublibrary / batch detection ─────────────────────────────────────────
  # Embryos from different sequencing sublibraries carry a batch effect that
  # can dominate over the genotype signal when both are present in the same
  # arm. batch_of_bare was captured from the RAW column names before prefix
  # stripping (see top of script) — a simple lookup here, not re-derivation,
  # since the bare names in `use_cols` no longer carry the "__" separator
  # real batch detection needs.
  batches <- unname(batch_of_bare[use_cols])
  n_batches <- length(unique(batches))

  if (n_batches > 1) {
    batch <- factor(batches)
    # Safety check: batch must not be perfectly confounded with condition
    tbl <- table(batch, condition)
    confounded <- any(apply(tbl, 1, function(r) sum(r > 0) == 1) &
                      apply(tbl, 2, function(r) sum(r > 0) == 1))
    if (confounded) {
      cat(sprintf("[limma]   NOTE: %d sublibraries detected but batch is confounded with condition — fitting without batch term\n", n_batches))
      batch <- NULL
    } else {
      cat(sprintf("[limma]   %d sublibraries detected — adding batch covariate: %s\n",
                  n_batches, paste(levels(batch), collapse=", ")))
    }
  } else {
    batch <- NULL
  }

  # Low-expression filter using edgeR's filterByExpr
  dge <- DGEList(counts=cnt)
  keep <- filterByExpr(dge, group=condition,
                       min.count=MIN_CPM, min.total.count=MIN_CPM*2,
                       min.prop=0)
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  cat(sprintf("[limma]   %d / %d genes pass expression filter\n",
              sum(keep), length(keep)))

  if (sum(keep) < 10) {
    cat("[limma]   SKIP: too few genes after filtering\n")
    return(NULL)
  }

  # TMM normalisation (edgeR)
  dge <- calcNormFactors(dge, method="TMM")

  # Design matrix — condition ± batch covariate
  if (!is.null(batch)) {
    design <- model.matrix(~ batch + condition)
  } else {
    design <- model.matrix(~ condition)
  }

  # voom + eBayes
  v    <- voom(dge, design, plot=FALSE)
  fit  <- lmFit(v, design)
  fit  <- eBayes(fit, trend=TRUE, robust=TRUE)

  # Extract results for the condition coefficient
  coef_name <- "conditionA"
  if (!coef_name %in% colnames(fit$coefficients)) {
    cat(sprintf("[limma]   ERROR: coefficient '%s' not found. Available: %s\n",
                coef_name, paste(colnames(fit$coefficients), collapse=", ")))
    return(NULL)
  }
  res <- topTable(fit, coef=coef_name,
                  number=Inf, sort.by="P", adjust.method="BH")

  # Rename columns for consistency with Python output
  res$gene_id  <- rownames(res)
  res$symbol   <- sapply(rownames(res), sym)
  res$contrast <- def$name
  names(res)[names(res) == "logFC"]     <- "log2FoldChange"
  names(res)[names(res) == "adj.P.Val"] <- "padj"
  names(res)[names(res) == "P.Value"]   <- "pvalue"
  names(res)[names(res) == "AveExpr"]   <- "baseMean"
  res <- res[, c("gene_id","symbol","baseMean","log2FoldChange",
                 "t","pvalue","padj","contrast")]

  # Write full results
  cdir <- file.path(out_base, def$name)
  dir.create(cdir, recursive=TRUE, showWarnings=FALSE)
  write.table(res, file.path(cdir, "results_full.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)

  # Significant hits
  sig <- res[!is.na(res$padj) & res$padj < FDR_THRESH &
               abs(res$log2FoldChange) > LFC_THRESH, ]
  sig <- sig[order(-sig$log2FoldChange), ]
  write.table(sig, file.path(cdir, "results_sig.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)

  n_up   <- sum(sig$log2FoldChange > 0, na.rm=TRUE)
  n_down <- sum(sig$log2FoldChange < 0, na.rm=TRUE)
  cat(sprintf("[limma]   → %d sig (padj<%.2f, |lFC|>%.2f): %d up-A / %d up-B\n",
              nrow(sig), FDR_THRESH, LFC_THRESH, n_up, n_down))

  # Contrast info
  batch_note <- if (!is.null(batch))
    sprintf("Batch covariate included (%d sublibraries: %s)",
            n_batches, paste(levels(batch), collapse=", "))
  else "No batch covariate (single sublibrary or confounded)"

  writeLines(c(
    sprintf("Contrast: %s", def$name),
    sprintf("Description: %s", def$desc),
    sprintf("Method: limma-voom + eBayes (trend=TRUE, robust=TRUE)"),
    batch_note,
    sprintf("Arm A (%d embryos): %s", length(A), paste(short_name(A), collapse=", ")),
    sprintf("Arm B (%d embryos): %s", length(B), paste(short_name(B), collapse=", ")),
    sprintf("Genes tested: %d", nrow(res)),
    sprintf("Significant (padj<%.2f, |lFC|>%.3f): %d (%d up-A, %d up-B)",
            FDR_THRESH, LFC_THRESH, nrow(sig), n_up, n_down)
  ), file.path(cdir, "contrast_info.txt"))

  # ── Volcano plot ────────────────────────────────────────────────────────
  # Uses ggplot2 + ggrepel for non-overlapping labels.
  # Falls back to base-R text() if ggrepel is not installed.
  use_ggrepel <- requireNamespace("ggplot2", quietly=TRUE) &&
                 requireNamespace("ggrepel", quietly=TRUE)

  draw_volcano <- function() {
    plot_res <- res[!is.na(res$padj), ]
    plot_res$neg_log10p <- -log10(pmax(plot_res$padj, 1e-300))
    plot_res$sig_class <- ifelse(
      plot_res$padj < FDR_THRESH & plot_res$log2FoldChange >  LFC_THRESH, "up",
      ifelse(plot_res$padj < FDR_THRESH & plot_res$log2FoldChange < -LFC_THRESH, "down", "ns"))

    # Label top genes by a combined rank score: -log10(padj) * |lFC|
    # Pick top N/2 from each direction separately so both sides are represented,
    # but fall back gracefully when one side has fewer than N/2 genes.
    sig_up   <- plot_res[plot_res$sig_class == "up",   ]
    sig_down <- plot_res[plot_res$sig_class == "down",  ]
    score <- function(d) d$neg_log10p * abs(d$log2FoldChange)
    n_each <- max(1L, TOP_LABEL_N %/% 2)
    top_up   <- if (nrow(sig_up)   > 0) sig_up[order(-score(sig_up)),  ][seq_len(min(n_each, nrow(sig_up))),   ] else sig_up[0,]
    top_down <- if (nrow(sig_down) > 0) sig_down[order(-score(sig_down)),][seq_len(min(n_each, nrow(sig_down))),] else sig_down[0,]
    to_label <- rbind(top_up, top_down)

    title_str <- gsub("_", " ", def$name)

    if (use_ggrepel) {
      library(ggplot2)
      library(ggrepel)
      col_map <- c("up"="#e31a1c", "down"="#1f78b4", "ns"="#cccccc")
      sz_map  <- c("up"=1.5,       "down"=1.5,        "ns"=0.8)
      alp_map <- c("up"=0.9,       "down"=0.9,        "ns"=0.4)

      plot_res$label <- ifelse(plot_res$gene_id %in% to_label$gene_id,
                               plot_res$symbol, "")
      p <- ggplot(plot_res, aes(x=log2FoldChange, y=neg_log10p,
                                colour=sig_class, size=sig_class,
                                alpha=sig_class)) +
        geom_point() +
        scale_colour_manual(values=col_map, guide="none") +
        scale_size_manual(values=sz_map, guide="none") +
        scale_alpha_manual(values=alp_map, guide="none") +
        geom_hline(yintercept=-log10(FDR_THRESH), lty=2, colour="#999999", linewidth=0.7) +
        geom_vline(xintercept=c(-LFC_THRESH, LFC_THRESH), lty=3,
                   colour=c("#1f78b4","#e31a1c"), linewidth=0.7) +
        geom_vline(xintercept=0, colour="#555555", linewidth=0.4) +
        geom_text_repel(
          aes(label=label),
          size=3, fontface="italic", colour="black",
          max.overlaps=30, box.padding=0.4, point.padding=0.3,
          segment.colour="#555555", segment.size=0.3,
          min.segment.length=0.2,
          seed=42
        ) +
        annotate("text", x=max(plot_res$log2FoldChange, na.rm=TRUE)*0.9,
                 y=max(plot_res$neg_log10p, na.rm=TRUE)*0.97,
                 label=sprintf("Up-A: %d", n_up), colour="#e31a1c", size=3.5, hjust=1) +
        annotate("text", x=min(plot_res$log2FoldChange, na.rm=TRUE)*0.9,
                 y=max(plot_res$neg_log10p, na.rm=TRUE)*0.97,
                 label=sprintf("Up-B: %d", n_down), colour="#1f78b4", size=3.5, hjust=0) +
        labs(x=expression(log[2]*" fold change  (A / B)"),
             y=expression(-log[10]*"(adjusted p-value)"),
             title=title_str,
             subtitle=def$desc) +
        theme_classic(base_size=12) +
        theme(plot.title=element_text(face="bold", size=11),
              plot.subtitle=element_text(size=9, colour="#444444"))
      print(p)
    } else {
      # Base-R fallback with staggered labels to reduce overlap
      is_up   <- plot_res$sig_class == "up"
      is_down <- plot_res$sig_class == "down"
      is_ns   <- plot_res$sig_class == "ns"
      xlim <- range(plot_res$log2FoldChange, na.rm=TRUE)
      ylim <- c(0, max(plot_res$neg_log10p, na.rm=TRUE) * 1.15)
      plot(plot_res$log2FoldChange[is_ns], plot_res$neg_log10p[is_ns],
           pch=16, cex=0.4, col="#cccccc", xlim=xlim, ylim=ylim,
           xlab=expression(log[2]*" fold change  (A / B)"),
           ylab=expression(-log[10]*"(adjusted p-value)"),
           main=title_str, cex.main=0.9)
      points(plot_res$log2FoldChange[is_up],  plot_res$neg_log10p[is_up],
             pch=16, cex=0.7, col="#e31a1c")
      points(plot_res$log2FoldChange[is_down], plot_res$neg_log10p[is_down],
             pch=16, cex=0.7, col="#1f78b4")
      abline(h=-log10(FDR_THRESH), lty=2, col="#999999", lwd=0.9)
      abline(v=c(-LFC_THRESH, LFC_THRESH), lty=3,
             col=c("#1f78b4","#e31a1c"), lwd=0.9)
      abline(v=0, lty=1, col="#555555", lwd=0.5)
      # Stagger labels vertically by alternating y offset to reduce overlap
      if (nrow(to_label) > 0) {
        y_offsets <- rep(c(0.05, 0.12, -0.04), length.out=nrow(to_label)) *
                     diff(ylim)
        text(to_label$log2FoldChange,
             to_label$neg_log10p + y_offsets,
             labels=to_label$symbol,
             cex=0.6, font=3,
             pos=ifelse(to_label$log2FoldChange > 0, 4, 2))
      }
      legend("topleft", bty="n", cex=0.85,
             legend=c(sprintf("Up-A: %d", n_up), sprintf("Up-B: %d", n_down)),
             text.col=c("#e31a1c","#1f78b4"))
    }
  }

  pdf(file.path(cdir, "volcano.pdf"), width=8, height=7)
  tryCatch(draw_volcano(), finally=dev.off())
  png(file.path(cdir, "volcano.png"), width=800, height=700, res=100)
  tryCatch(draw_volcano(), finally=dev.off())

  return(list(contrast=def$name, n_tested=nrow(res), batch_corrected=!is.null(batch),
              n_sig=nrow(sig), n_up_A=n_up, n_down_A=n_down,
              description=def$desc))
}

# =============================================================================
# Run all contrasts
# =============================================================================

dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)
summaries <- list()

for (def in contrast_defs) {
  if (!def$id %in% wanted_contrasts) next
  cat(sprintf("\n[limma] ── Contrast %d: %s ──\n", def$id, def$name))
  cat(sprintf("[limma]   %s\n", def$desc))
  result <- tryCatch(
    run_contrast(def, counts, OUT_DIR),
    error = function(e) { cat("[limma]   ERROR:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(result)) summaries[[length(summaries)+1]] <- result
}

# Summary table
if (length(summaries)) {
  sumdf <- do.call(rbind, lapply(summaries, as.data.frame))
  write.table(sumdf, file.path(OUT_DIR, "summary_table.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)
  cat("\n[limma] === Summary ===\n")
  print(sumdf[, c("contrast","n_tested","n_sig","n_up_A","n_down_A")],
        row.names=FALSE, right=FALSE)
}

cat(sprintf("\n[limma] Done. Results in %s/\n", OUT_DIR))
cat("[limma] To add results to volcano plots from the Python pipeline,\n")
cat("[limma] use results_full.tsv from deg_limma/ alongside deg/ outputs.\n")
