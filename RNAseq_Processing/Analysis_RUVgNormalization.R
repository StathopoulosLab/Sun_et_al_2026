#!/usr/bin/env Rscript
# =============================================================================
# ruvg_nextgen.R
#
# Fully corrected and robust RUVg pipeline
#
# Fixes:
#   ✓ No DESeq2 rank-deficiency errors (uses design = ~1)
#   ✓ Robust metadata/count reconciliation
#   ✓ Drops missing samples automatically
#   ✓ Handles sparse groups safely
#   ✓ Auto-selects empirical control genes
#   ✓ Optional TF-bound gene exclusion
#   ✓ Automatic k optimization
#   ✓ Overcorrection diagnostics
#   ✓ Safe PCA / silhouette handling
#
# =============================================================================

suppressPackageStartupMessages({
  library(RUVSeq)
  library(EDASeq)
  library(optparse)
  library(DESeq2)
  library(matrixStats)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(cluster)
  library(gridExtra)
})

# =============================================================================
# OPTIONS
# =============================================================================

option_list <- list(
  make_option("--counts", type="character",
              default="results/combined/all_genotypes_counts_matrix_qc.tsv"),
  
  make_option("--metadata", type="character",
              default="results/combined/sample_metadata.tsv"),
  
  make_option("--outdir", type="character",
              default="results/combined/ruvg_nextgen"),
  
  make_option("--tf-bound", type="character",
              default=NULL),
  
  make_option("--max-k", type="integer",
              default=5),
  
  make_option("--n-controls", type="integer",
              default=300),
  
  make_option("--min-count", type="integer",
              default=2),
  
  make_option("--min-samples", type="integer",
              default=3),

  make_option("--protect-genes", type="character",
              default=NULL,
              help="Comma-separated FBgn IDs to force past the min-count/min-samples filter regardless of expression. Mirrors PROTECTED_GROUPS in step5c_qc_filter.py and _SENTINEL_GENES in step6_tpm.py.")
)

opt <- parse_args(OptionParser(option_list=option_list))

dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(opt$outdir,"figures"),
           recursive=TRUE, showWarnings=FALSE)

cat("====================================\n")
cat("RUVg NextGen Corrected Pipeline\n")
cat("====================================\n")

# =============================================================================
# LOAD COUNTS
# =============================================================================

counts <- read.table(
  opt$counts,
  sep="\t",
  header=TRUE,
  row.names=1,
  check.names=FALSE,
  quote=""
)

if ("GeneSymbol" %in% colnames(counts)) {
  counts <- counts[, colnames(counts)!="GeneSymbol", drop=FALSE]
}

cat("Loaded counts:", nrow(counts), "genes x",
    ncol(counts), "samples\n")

# =============================================================================
# LOAD METADATA
# =============================================================================

meta <- read.table(
  opt$metadata,
  sep="\t",
  header=TRUE,
  stringsAsFactors=FALSE,
  check.names=FALSE,
  quote=""
)

required <- c("sample","group","batch")
missing_cols <- setdiff(required, colnames(meta))

if (length(missing_cols) > 0) {
  stop(
    paste("Metadata missing required columns:",
          paste(missing_cols, collapse=", "))
  )
}

# Keep shared samples only
shared <- intersect(colnames(counts), meta$sample)

if (length(shared) == 0) {
  stop("No shared sample names between counts and metadata.")
}

dropped_counts <- setdiff(colnames(counts), shared)
if (length(dropped_counts) > 0) {
  cat("Dropping count samples absent from metadata:\n")
  cat(paste(dropped_counts, collapse="\n"), "\n")
}

counts <- counts[, shared, drop=FALSE]

meta <- meta[match(shared, meta$sample), , drop=FALSE]
rownames(meta) <- meta$sample

meta$group <- factor(meta$group)
meta$batch <- factor(meta$batch)

cat("Samples retained:", ncol(counts), "\n")

# =============================================================================
# FILTER GENES
# =============================================================================

keep <- rowSums(counts >= opt$`min-count`) >= opt$`min-samples`

# Force-include protected genes regardless of whether they clear the
# generic expression filter above. Without this, a gene central to the
# study (e.g. one a genotype is named after) can be silently dropped by a
# blanket statistical threshold with no recourse -- exactly the gap
# PROTECTED_GROUPS (step5c_qc_filter.py) and _SENTINEL_GENES (step6_tpm.py)
# already close elsewhere in this pipeline.
if (!is.null(opt$`protect-genes`)) {
  protect_ids <- trimws(strsplit(opt$`protect-genes`, ",")[[1]])
  found   <- protect_ids[protect_ids %in% rownames(counts)]
  missing <- setdiff(protect_ids, rownames(counts))
  already_kept <- found[keep[found]]
  rescued      <- found[!keep[found]]
  if (length(rescued) > 0) {
    keep[rescued] <- TRUE
    cat("Protected genes rescued from filter (were below min-count/min-samples):\n")
    cat(paste(" ", rescued, collapse="\n"), "\n")
  }
  if (length(already_kept) > 0) {
    cat("Protected genes already passing filter:", paste(already_kept, collapse=", "), "\n")
  }
  if (length(missing) > 0) {
    cat(paste0("WARNING: protected gene(s) not found in counts matrix at all ",
               "(this is upstream of the RUVg filter -- check featureCounts/gene ",
               "annotation, not this script):\n"))
    cat(paste(" ", missing, collapse="\n"), "\n")
  }
}

counts <- counts[keep, , drop=FALSE]

cat("Genes retained:", nrow(counts), "\n")

# =============================================================================
# NORMALIZATION (SAFE DESIGN = ~1)
# =============================================================================

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData   = meta,
  design    = ~ 1
)

dds <- estimateSizeFactors(dds)
norm <- counts(dds, normalized=TRUE)

# =============================================================================
# CONTROL GENE SELECTION
# =============================================================================

cat("Selecting empirical controls...\n")

gene_mean <- rowMeans(norm)
gene_cv   <- rowSds(as.matrix(norm)) / pmax(gene_mean,1)

group_p <- apply(log2(norm+1),1,function(x){
  
  if (length(unique(meta$group)) < 2) return(1)
  
  fit <- tryCatch(
    lm(x ~ meta$group),
    error=function(e) NULL
  )
  
  if (is.null(fit)) return(1)
  
  pv <- tryCatch(
    anova(fit)$`Pr(>F)`[1],
    error=function(e) 1
  )
  
  ifelse(is.na(pv),1,pv)
})

controls <- data.frame(
  gene=rownames(norm),
  mean=gene_mean,
  cv=gene_cv,
  group_p=group_p,
  stringsAsFactors=FALSE
)

# Relaxed thresholds for small multi-genotype datasets where strict
# group-invariance filtering leaves too few candidates.
# Strategy: take moderately expressed genes (mean > 5) with the weakest
# group signal (top half by p-value), then rank by low CV + high mean.
controls <- controls %>%
  filter(mean > 5) %>%
  filter(group_p > 0.2)

# If still too few, relax further — take top-300 least variable regardless
if (nrow(controls) < 50) {
  controls <- data.frame(
    gene    = rownames(norm),
    mean    = gene_mean,
    cv      = gene_cv,
    group_p = group_p,
    stringsAsFactors = FALSE
  ) %>%
  filter(mean > 5) %>%
  arrange(cv)
  cat("WARNING: relaxed to CV-only control selection\n")
}

# Optional TF-bound removal
if (!is.null(opt$`tf-bound`)) {
  
  tf <- scan(opt$`tf-bound`,
             what="character",
             quiet=TRUE)
  
  controls <- controls[!(controls$gene %in% tf), ]
}

controls$score <- rank(controls$cv) +
  rank(-controls$mean)

controls <- controls[order(controls$score), ]

controls <- head(controls, opt$`n-controls`)

control_genes <- controls$gene

if (length(control_genes) < 10) {
  stop("Too few control genes found.")
}

write.table(
  controls,
  file.path(opt$outdir,"controls_selected.tsv"),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

cat("Controls selected:", length(control_genes), "\n")

# =============================================================================
# BUILD RUV OBJECT
# =============================================================================

x <- newSeqExpressionSet(
  counts=round(as.matrix(counts)),
  phenoData=AnnotatedDataFrame(meta)
)

# =============================================================================
# SAFE SCORING FUNCTION
# =============================================================================

score_model <- function(mat) {
  
  logmat <- t(log2(pmax(mat,0)+1))
  
  # remove zero variance genes
  vars <- apply(logmat,2,sd)
  logmat <- logmat[, vars > 0, drop=FALSE]
  
  pca <- prcomp(logmat)
  
  pc1 <- pca$x[,1]
  
  batch_r2 <- 0
  group_r2 <- 0
  sil_mean <- 0
  
  # batch association
  if (length(unique(meta$batch)) > 1) {
    batch_r2 <- summary(lm(pc1 ~ meta$batch))$r.squared
  }
  
  # group association
  if (length(unique(meta$group)) > 1) {
    group_r2 <- summary(lm(pc1 ~ meta$group))$r.squared
  }
  
  # silhouette if each group has >=2
  tab <- table(meta$group)
  
  if (length(tab) > 1 && all(tab >= 2)) {
    
    sil <- silhouette(
      as.numeric(meta$group),
      dist(logmat)
    )
    
    sil_mean <- mean(sil[,3])
  }
  
  score <- sil_mean + group_r2 - batch_r2
  
  data.frame(
    batch_r2=batch_r2,
    group_r2=group_r2,
    silhouette=sil_mean,
    score=score
  )
}

# =============================================================================
# OPTIMIZE K
# =============================================================================

cat("Testing k values...\n")

diag_list <- list()

for (k in 0:opt$`max-k`) {
  
  cat("k =", k, "\n")
  
  if (k == 0) {
    
    mat <- norm
    
  } else {
    
    if (length(control_genes) <= k) next
    
    fit <- tryCatch(
      RUVg(x, control_genes, k=k),
      error=function(e) NULL
    )
    
    if (is.null(fit)) next
    
    mat <- normCounts(fit)
  }
  
  s <- score_model(mat)
  s$k <- k
  
  diag_list[[length(diag_list)+1]] <- s
}

diag <- bind_rows(diag_list)

if (nrow(diag) == 0) {
  stop("No successful k fits.")
}

best_k <- diag$k[which.max(diag$score)]

write.table(
  diag,
  file.path(opt$outdir,"k_diagnostics.tsv"),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

cat("Best k =", best_k, "\n")

# =============================================================================
# FINAL MODEL
# =============================================================================

if (best_k == 0) {
  
  final_counts <- norm
  
} else {
  
  final <- RUVg(x, control_genes, k=best_k)
  
  final_counts <- normCounts(final)
  
  W <- pData(final)[,grep("^W_",colnames(pData(final))),drop=FALSE]
  
  write.table(
    W,
    file.path(opt$outdir,"W_factors.tsv"),
    sep="\t",
    quote=FALSE,
    col.names=NA
  )
}

write.table(
  as.data.frame(final_counts),
  file.path(opt$outdir,"corrected_counts.tsv"),
  sep="\t",
  quote=FALSE,
  col.names=NA
)

# =============================================================================
# PLOTS
# =============================================================================

p1 <- ggplot(diag, aes(k, score)) +
  geom_line() +
  geom_point(size=3) +
  geom_vline(xintercept=best_k,
             linetype=2,
             colour="red") +
  theme_bw(base_size=12) +
  labs(title="k optimization",
       x="k",
       y="score")

ggsave(
  file.path(opt$outdir,"figures","k_optimization.pdf"),
  p1,
  width=7,
  height=5
)

cat("====================================\n")
cat("DONE\n")
cat("Samples:", ncol(counts), "\n")
cat("Genes:", nrow(counts), "\n")
cat("Controls:", length(control_genes), "\n")
cat("Best k:", best_k, "\n")
cat("Output:", opt$outdir, "\n")
cat("====================================\n")