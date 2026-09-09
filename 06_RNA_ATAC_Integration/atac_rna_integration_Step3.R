#!/usr/bin/env Rscript
# =============================================================================
# atac_rna_integration_v3.R
#
# Expanded bistability analysis + temporal chromatin-expression integration.
#
# KEY FIXES vs v2
# ---------------
# 1. Sign convention: v5 of Step1_Vent already flips the DESeq2/limma contrast
#    so positive log2FC = BOTCv more open. BOTv_vs_BOTCv stems NO LONGER need
#    negation. Set USE_V5_ATAC=TRUE (default) to avoid double-negation.
#    Set USE_V5_ATAC=FALSE if your ATAC_DARS were produced by v4 or earlier.
#
# 2. Label clipping fixed throughout with expansion margins and clip="off".
#
# PARTS
# -----
#   PART 5b — Expanded bistability: annotated gene sets, multi-panel temporal
#             view, functional grouping of concordant bistable genes, and the
#             full complement of v5 temporal DARs (gastr comparisons included).
#
#   PART 6  — Temporal lead/lag analysis: for each gene with ATAC data at both
#             nc14b AND nc14late, classify whether chromatin difference precedes,
#             follows, or co-emerges with expression difference.
#             Four classes:
#               Chromatin-first : ATAC sig at nc14b, RNA not sig until late
#               Expression-first: RNA sig at nc14b, ATAC not sig until late
#               Co-emerging     : both become sig at nc14late
#               Stable          : both sig at nc14b already
#             Outputs a temporal concordance heatmap and per-class gene lists.
#
# USAGE
# -----
#   Rscript atac_rna_integration_v3.R
#   Rscript atac_rna_integration_v3.R --parts 5b,6
#   Rscript atac_rna_integration_v3.R --no-v5  # use v4 sign convention
#   Rscript atac_rna_integration_v3.R --bistable-thresh 0.4
# =============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(ggplot2); library(dplyr)
  library(patchwork);     library(scales);  library(tidyr)
})
has_ggrepel <- requireNamespace("ggrepel", quietly=TRUE)
if (has_ggrepel) suppressPackageStartupMessages(library(ggrepel))

# =============================================================================
# CONFIG
# =============================================================================

DAR_DIR     <- "./ATAC_DARS"
DEG_DIR     <- "../../seRNA-seq_03_19_26/celseq2_pipeline/results/combined_figures/qc_and_deg_botv_botcv/deg_limma"
CLUSTER_DIR <- "../../seRNA-seq_03_19_26/celseq2_pipeline/results/combined_figures/explore_clusters_botv_botcv"
OUT_DIR     <- "./Overview_Plots/atac_rna_integration"

# v5 of Step1_Vent already has positive=BOTCv open — no negate needed.
# Set FALSE only if using DARs produced by Step1 v4 or earlier.
USE_V5_ATAC     <- TRUE

# =============================================================================
# THRESHOLD PRESETS
# =============================================================================
# Single-embryo data is inherently sparse. Use these named presets rather than
# arbitrary numbers, and choose based on n-per-arm and analysis purpose.
#
# ── RNA-seq thresholds ───────────────────────────────────────────────────────
# With CEL-Seq2 single-embryo data the binding constraint is always power, not
# noise. With n=2–3 per arm, limma eBayes can only find large effects; with
# n=6–7 it finds moderate effects too.
#
#   RELAXED   FDR=0.20, |lFC|=0       Use when: n≤3 per arm, exploratory,
#                                       or as input to enrichment/integration.
#                                       Accepts ~20% FDR; lFC=0 means padj
#                                       is the only gate (less harsh multiple-
#                                       testing penalty). Best for gene set
#                                       membership and ranking, not final calls.
#
#   STANDARD  FDR=0.10, |lFC|=0.5     Use when: n=3–4 per arm, most analyses.
#                                       log2(1.5)× ≈ 50% change. Appropriate
#                                       for the BOTCv temporal contrast (n=2+3)
#                                       and cross-genotype late (n=2+2).
#
#   STRINGENT FDR=0.05, |lFC|=1.0     Use when: n≥5 per arm, figures or
#                                       validation prioritisation. 2-fold
#                                       change, 5% FDR. The BOTv nc14b arm
#                                       (n=7) can support this; most others
#                                       will give few or zero hits.
#
# ── ATAC-seq thresholds ──────────────────────────────────────────────────────
# n=2 per group throughout; r=0.96–0.97 between replicates (excellent for n=2).
# DESeq2+limma consensus + RUVg + conf_score weighting already accounts for the
# low n. The key parameter is the lFC threshold:
#
#   RELAXED   |lFC|=0.3, no padj gate  Use for bistability classification and
#                                       peak-set membership. Captures real but
#                                       modest accessibility differences. At n=2
#                                       the padj is conservative; lFC is more
#                                       informative than padj here.
#
#   STANDARD  |lFC|=0.5, padj<0.10    Use for concordance analysis and lead/lag
#                                       classification. Matches DESeq2+limma
#                                       consensus threshold in your pipeline.
#
#   STRINGENT |lFC|=1.0, padj<0.05    Use for figures and validation. A 2-fold
#                                       accessibility change at n=2 with padj<0.05
#                                       is very high confidence.
#
# Current active thresholds (change here or via CLI flags):
RNA_FDR_RELAXED   <- 0.20; RNA_LFC_RELAXED   <- 0.0
RNA_FDR_STANDARD  <- 0.10; RNA_LFC_STANDARD  <- 0.25
RNA_FDR_STRINGENT <- 0.05; RNA_LFC_STRINGENT <- 0.5

ATAC_LFC_RELAXED  <- 0.25
ATAC_LFC_STANDARD <- 0.5
ATAC_LFC_RELAXED  <- 0.25
ATAC_LFC_STANDARD <- 0.5
ATAC_LFC_STRINGENT<- 1.0

# Active thresholds used throughout this script
BISTABLE_THRESH <- ATAC_LFC_STANDARD   # 0.5 — bistable peak classification
STAB_THRESH     <- ATAC_LFC_RELAXED    # 0.3 — "stable" peak ceiling
FDR_RNA         <- RNA_FDR_STANDARD    # 0.10
LFC_RNA_THRESH  <- RNA_LFC_STANDARD    # 0.5 — for lead/lag classification
TSS_WINDOW      <- 5000
MIN_OVERLAP     <- 1L
TOP_N           <- 30

# v5 temporal stems (all produced by 02_dar_calling_ventralized_temporal.r)
# Sign: positive = BOTCv more open (for cross-genotype), positive = nc14b more open
# (for temporal within-genotype, because grpA=late, grpB=early → positive=grpB=early)
TEMPORAL_STEMS <- list(
  # Cross-genotype
  list(stem="BOTv_vs_BOTCv",           label="BOTCv vs BOTv (nc14b)",
       negate=FALSE, type="cross_geno", tp="nc14b",
       note="positive = BOTCv more open"),
  list(stem="BOTv_vs_BOTCv_nc14late",  label="BOTCv vs BOTv (nc14late)",
       negate=FALSE, type="cross_geno", tp="nc14late",
       note="positive = BOTCv more open"),
  # Temporal within-genotype
  # grpA=nc14late, grpB=nc14b → positive = nc14b more open = closing over time
  list(stem="BOTv_temporal",            label="BOTv temporal (nc14b→nc14late)",
       negate=FALSE, type="temporal",  tp="BOTv",
       note="positive = nc14b open = losing accessibility over time"),
  list(stem="BOTCv_temporal",           label="BOTCv temporal (nc14b→nc14late)",
       negate=FALSE, type="temporal",  tp="BOTCv",
       note="positive = nc14b open = losing accessibility over time"),
  # Gastr comparisons (BOTv only, from v5)
  list(stem="BOTv_gastr_vs_nc14b",      label="BOTv gastr vs nc14b",
       negate=FALSE, type="temporal",  tp="BOTv_gastr",
       note="positive = nc14b more open"),
  list(stem="BOTv_gastr_vs_nc14late",   label="BOTv gastr vs nc14late",
       negate=FALSE, type="temporal",  tp="BOTv_gastr_late",
       note="positive = nc14late more open")
)

# RNA contrasts — same as main script
RNA_CONTRASTS <- list(
  list(name="BOTCv_nc14b_vs_BOTv_nc14b",  label="BOTCv vs BOTv (nc14b)",
       tp_atac="nc14b",    positive_is="BOTCv"),
  list(name="BOTCv_late_vs_BOTv_nc14d",   label="BOTCv late vs BOTv nc14d",
       tp_atac="nc14late", positive_is="BOTCv"),
  list(name="BOTCv_late_vs_BOTCv_nc14b",  label="BOTCv temporal",
       tp_atac="BOTCv",    positive_is="BOTCv_late"),
  list(name="BOTv_nc14d_vs_BOTv_nc14b",   label="BOTv temporal",
       tp_atac="BOTv",     positive_is="BOTv_late"),
  list(name="BOTCv_late_vs_BOTv_nc14b",   label="BOTCv late vs BOTv early",
       tp_atac="nc14b",    positive_is="BOTCv")
)

# Colours
COL_BOTV <- "#2ca02c"; COL_BOTCV <- "#e377c2"
COL_OPEN <- "#e31a1c"; COL_CLOSE <- "#1f78b4"

# Temporal colour palette — light=early, dark=late within each genotype
COL_BOTCV_EARLY <- "#e377c2"   # BOTCv nc14b
COL_BOTCV_LATE  <- "#9e1f8e"   # BOTCv late (dark magenta)
COL_BOTV_EARLY  <- "#2ca02c"   # BOTv  nc14b
COL_BOTV_LATE   <- "#1a5c1a"   # BOTv  late (dark forest green)
COL_CROSS       <- "#7b4fb7"   # cross-genotype (purple)

# Lookup: contrast name → display colours for A-arm and B-arm
CONTRAST_ARM_COLS <- list(
  "BOTCv_late_vs_BOTCv_nc14b" = c(A=COL_BOTCV_LATE,  B=COL_BOTCV_EARLY),
  "BOTv_nc14d_vs_BOTv_nc14b"  = c(A=COL_BOTV_LATE,   B=COL_BOTV_EARLY),
  "BOTCv_late_vs_BOTv_nc14b"  = c(A=COL_BOTCV_LATE,  B=COL_BOTV_EARLY),
  "BOTCv_nc14b_vs_BOTv_nc14b" = c(A=COL_BOTCV_EARLY, B=COL_BOTV_EARLY),
  "BOTCv_late_vs_BOTv_nc14d"  = c(A=COL_BOTCV_LATE,  B=COL_BOTV_LATE)
)

# =============================================================================
# CLI
# =============================================================================
args <- commandArgs(trailingOnly=TRUE)
parse_arg <- function(flag, default) {
  i <- which(args==flag); if (length(i) && length(args)>i) args[i+1] else default
}
run_parts_raw   <- parse_arg("--parts","5b,6")
run_parts       <- strsplit(run_parts_raw,",")[[1]]
BISTABLE_THRESH <- as.numeric(parse_arg("--bistable-thresh", BISTABLE_THRESH))
TSS_WINDOW      <- as.integer(parse_arg("--tss-window",      TSS_WINDOW))
if ("--no-v5" %in% args) USE_V5_ATAC <- FALSE

cat(sprintf("[v3] Parts: %s  |  v5_atac=%s  |  bistable_thresh=%.1f\n",
            run_parts_raw, USE_V5_ATAC, BISTABLE_THRESH))
if (!USE_V5_ATAC) {
  cat("[v3] WARNING: --no-v5 set — applying negate to BOTv_vs_BOTCv stems\n")
  for (i in seq_along(TEMPORAL_STEMS))
    if (grepl("BOTv_vs_BOTCv", TEMPORAL_STEMS[[i]]$stem))
      TEMPORAL_STEMS[[i]]$negate <- TRUE
}

for (d in c("part5b_bistable_expanded","part6_temporal_lead_lag"))
  dir.create(file.path(OUT_DIR,d), recursive=TRUE, showWarnings=FALSE)

# =============================================================================
# HELPERS — label clipping fix throughout
# =============================================================================

save_plot <- function(stem, p, w=10, h=7, dpi=200) {
  ggplot2::ggsave(paste0(stem,".pdf"), p, width=w, height=h, device="pdf")
  ggplot2::ggsave(paste0(stem,".png"), p, width=w, height=h, dpi=dpi)
  cat(sprintf("[v3]   Saved %s\n", basename(stem)))
}

# Fix label clipping on coord_flip bars.
# expand_right: extra fraction on positive side (for labels extending right)
# expand_left:  extra fraction on negative side (for labels extending left,
#               e.g. NES bars or RNA lFC bars that go both +/-)
# plot.margin: left margin ensures y-axis labels aren't cut off by the device
bar_theme_fix <- function(expand_right=0.40, expand_left=0.30) {
  list(
    scale_y_continuous(expand=expansion(mult=c(expand_left, expand_right))),
    coord_flip(clip="off"),
    theme(plot.margin=margin(5, 20, 5, 20, "pt"))
  )
}

read_rna <- function(contrast_name) {
  path <- file.path(DEG_DIR, contrast_name, "results_full.tsv")
  if (!file.exists(path)) { message("[v3] Not found: ",path); return(NULL) }
  df <- tryCatch(read.table(path,header=TRUE,sep="\t",quote="",stringsAsFactors=FALSE),
                 error=function(e) NULL)
  if (!is.null(df)) df <- df[!is.na(df$log2FoldChange) & !is.na(df$t), ]
  df
}

read_bed_dar <- function(stem, negate=FALSE) {
  # IMPORTANT: try annotated files FIRST. The raw _DARs.bed col5 is conf_score
  # (always positive), NOT log2FC. Using it as lfc would make all peaks
  # appear same-signed and bistability classification would always return 0.
  for (ext in c("_DARs_annotated.txt","_DARs_annotated_peaks.csv","_DARs.bed")) {
    path <- file.path(DAR_DIR, paste0(stem, ext))
    if (!file.exists(path)) next
    sep_ch <- if (grepl("\\.bed$",path)) "\t"
              else if (grepl("\t",readLines(path,1,warn=FALSE),fixed=TRUE)) "\t" else ","
    df <- tryCatch(read.table(path,header=(ext!=".bed"),sep=sep_ch,stringsAsFactors=FALSE,
                              quote="",fill=TRUE,comment.char="",row.names=NULL),
                   error=function(e) NULL)
    if (is.null(df)||nrow(df)==0) next
    if (ext=="_DARs.bed") {
      # Only use BED col5 if it has both positive AND negative values
      # (meaning it really is log2FC, not conf_score)
      lfc_raw <- if (ncol(df)>=5) suppressWarnings(as.numeric(df[,5])) else NULL
      if (is.null(lfc_raw)||!any(lfc_raw<0,na.rm=TRUE)) {
        cat(sprintf("[v3]   BED col5 for %s is all-positive (likely conf_score not lFC) — skipping\n", stem))
        next
      }
      nm <- if (ncol(df)>=4) as.character(df[,4]) else paste0(stem,"_",seq_len(nrow(df)))
      return(data.frame(chr=as.character(df[,1]), start=as.integer(df[,2]),
                        end=as.integer(df[,3]), name=nm,
                        atac_lfc=lfc_raw*ifelse(negate,-1,1),
                        atac_padj=rep(NA_real_,nrow(df)),
                        stringsAsFactors=FALSE))
    }
    colnames(df)[1] <- gsub("^[^[:alnum:]]+","",colnames(df)[1])
    lfc_col  <- intersect(c("log2FoldChange","log2FC","LFC"),colnames(df))[1]
    padj_col <- intersect(c("padj","adj.P.Val","p_adj"),colnames(df))[1]
    chr_col  <- intersect(c("seqnames","chr","Chr"),colnames(df))[1]
    st_col   <- intersect(c("start","Start"),colnames(df))[1]
    en_col   <- intersect(c("end","End"),colnames(df))[1]
    nm_col   <- intersect(c("name","Name","peak_id"),colnames(df))[1]
    if (is.na(lfc_col)||is.na(chr_col)) next
    lfc_vals <- as.numeric(df[[lfc_col]])
    cat(sprintf("[v3]   %s: log2FC range [%.2f, %.2f]\n",
                stem, min(lfc_vals,na.rm=TRUE), max(lfc_vals,na.rm=TRUE)))

    # Completeness check: this branch trusts the annotated file's own row
    # count with no reference to the current bed file -- if the annotated
    # file has fallen out of sync (e.g. regenerated with more peaks after
    # the annotated file was last written), this silently returns fewer
    # peaks than actually exist. Compare against the bed file explicitly.
    bed_check_path <- file.path(DAR_DIR, paste0(stem, "_DARs.bed"))
    if (file.exists(bed_check_path)) {
      n_bed <- tryCatch(length(readLines(bed_check_path)), error=function(e) NA_integer_)
      if (!is.na(n_bed) && n_bed != nrow(df)) {
        cat(sprintf("[v3]   [WARNING] %s: annotated file has %d peaks but current bed file has %d -- this annotated file is likely OUT OF DATE. Regenerate it before trusting these results.\n",
                   stem, nrow(df), n_bed))
      }
    }

    return(data.frame(
      chr  = as.character(df[[chr_col]]),
      start= as.integer(df[[st_col]]),
      end  = as.integer(df[[en_col]]),
      name = if (!is.na(nm_col)) as.character(df[[nm_col]]) else paste0(stem,"_",seq_len(nrow(df))),
      atac_lfc  = lfc_vals*ifelse(negate,-1,1),
      atac_padj = if (!is.na(padj_col)) as.numeric(df[[padj_col]]) else NA_real_,
      stringsAsFactors=FALSE))
  }
  message(sprintf("[v3] WARNING: No usable DAR file found for '%s'", stem))
  NULL
}

annotate_peaks <- function(peak_df, tss_window=TSS_WINDOW, stem_hint=NULL) {
  # Annotate peaks to nearby genes.
  # Priority order:
  #   1. ChIPseeker output already produced by your ATAC pipeline
  #      (looks in ATAC_DARS/ and Overview_Plots/annotation/vent/ for
  #       *_DARs_annotated_ChIPseeker.csv and *_DARs_annotated.txt)
  #   2. TxDb + org.Dm.eg.db genomic overlap (TSS ± tss_window)
  #
  # Adds columns: symbol, gene_id, peak_annotation, feature_class
  # feature_class: "Gene-proximal" | "CRM-associated" | NA
  peak_df$symbol         <- NA_character_
  peak_df$gene_id        <- NA_character_
  peak_df$peak_annotation<- NA_character_
  peak_df$feature_class  <- NA_character_

  classify_feature <- function(ann_vec) {
    ann <- tolower(as.character(ann_vec))
    dplyr::case_when(
      is.na(ann_vec)                                   ~ NA_character_,
      grepl("unspecified|unknown|^na$", ann)           ~ NA_character_,
      grepl("promoter|5.utr|tss|1st exon|cds", ann)   ~ "Gene-proximal",
      grepl("intergenic|distal|intron|downstream", ann)~ "CRM-associated",
      TRUE                                             ~ NA_character_)
  }

  # ── 1. Try ChIPseeker / existing annotation files ────────────────────────
  chipseeker_paths <- character(0)
  if (!is.null(stem_hint)) {
    chipseeker_paths <- c(
      file.path("./Overview_Plots/annotation/vent",
                paste0(stem_hint,"_DARs_annotated_ChIPseeker.csv")),
      file.path(DAR_DIR, paste0(stem_hint,"_DARs_annotated_ChIPseeker.csv")),
      file.path(DAR_DIR, paste0(stem_hint,"_DARs_annotated.txt")),
      file.path(DAR_DIR, paste0(stem_hint,"_DARs_annotated_peaks.csv")))
  }

  for (p in chipseeker_paths) {
    if (!file.exists(p)) next
    sep_ch <- if (grepl("\t",readLines(p,1,warn=FALSE),fixed=TRUE)) "\t" else ","
    ann <- tryCatch(
      read.table(p,header=TRUE,sep=sep_ch,stringsAsFactors=FALSE,
                 quote="",fill=TRUE,comment.char="",row.names=NULL),
      error=function(e) NULL)
    if (is.null(ann)||nrow(ann)==0) next
    colnames(ann)[1] <- gsub("^[^[:alnum:]]+","",colnames(ann)[1])

    sym_col  <- intersect(c("SYMBOL","geneId","gene_name","nearest_gene",
                             "Gene.Name","symbol","geneName"), colnames(ann))[1]
    name_col <- intersect(c("name","peakID","peak_id","Name","V4"), colnames(ann))[1]
    ann_col  <- intersect(c("annotation","Annotation","peak_type","feature"),
                           colnames(ann))[1]
    gid_col  <- intersect(c("geneId","GeneID","gene_id","GENEID","ENTREZID"),
                           colnames(ann))[1]

    if (is.na(sym_col) || is.na(name_col)) {
      cat(sprintf("[v3]   ChIPseeker file %s has no sym/name cols — skipping\n",
                  basename(p)))
      next
    }
    m <- match(peak_df$name, ann[[name_col]])
    valid <- !is.na(m)
    if (sum(valid)==0) {
      # Try matching by chr:start-end if name doesn't match
      if (!is.na(name_col) && all(c("chr","start","end") %in% colnames(ann))) {
        ann_key <- paste0(ann$chr,":",ann$start,"-",ann$end)
        pk_key  <- paste0(peak_df$chr,":",peak_df$start,"-",peak_df$end)
        m2 <- match(pk_key, ann_key)
        valid2 <- !is.na(m2)
        if (sum(valid2)>0) { m <- m2; valid <- valid2 }
      }
    }
    if (sum(valid)==0) next
    peak_df$symbol[valid]          <- ann[[sym_col]][m[valid]]
    if (!is.na(ann_col))
      peak_df$peak_annotation[valid] <- ann[[ann_col]][m[valid]]
    if (!is.na(gid_col))
      peak_df$gene_id[valid]         <- as.character(ann[[gid_col]][m[valid]])
    peak_df$feature_class <- classify_feature(peak_df$peak_annotation)
    n_sym <- sum(!is.na(peak_df$symbol))
    cat(sprintf("[v3]   ChIPseeker annotation: %d/%d peaks → genes (%s)\n",
                n_sym, nrow(peak_df), basename(p)))
    if (n_sym > nrow(peak_df)*0.05) return(peak_df)   # accept if >5% annotated
  }

  # ── 2. TxDb fallback ────────────────────────────────────────────────────────
  valid <- !is.na(peak_df$chr) & !is.na(peak_df$start) & !is.na(peak_df$end)
  if (sum(valid)==0) return(peak_df)
  if (!requireNamespace("TxDb.Dmelanogaster.UCSC.dm6.ensGene",quietly=TRUE)) {
    cat("[v3]   TxDb not available for fallback annotation\n")
    cat("[v3]   Install: BiocManager::install('TxDb.Dmelanogaster.UCSC.dm6.ensGene')\n")
    return(peak_df)
  }
  suppressPackageStartupMessages(
    library(TxDb.Dmelanogaster.UCSC.dm6.ensGene))
  txdb     <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
  genes_gr <- suppressWarnings(genes(txdb))
  tss_gr   <- resize(genes_gr, width=1L, fix="start")
  tss_gr   <- suppressWarnings(
    trim(resize(tss_gr, width=tss_window*2+1L, fix="center")))
  pk_gr    <- GRanges(seqnames=peak_df$chr[valid],
                       ranges=IRanges(peak_df$start[valid]+1L, peak_df$end[valid]))
  mcols(pk_gr)$row_idx <- which(valid)
  hits <- suppressWarnings(findOverlaps(pk_gr, tss_gr))
  if (length(hits)==0) { cat("[v3]   TxDb: no peak-gene overlaps\n"); return(peak_df) }
  hit_df <- data.frame(row_idx=mcols(pk_gr)$row_idx[queryHits(hits)],
                        entrez=names(tss_gr)[subjectHits(hits)],
                        stringsAsFactors=FALSE)
  hit_df <- hit_df[!duplicated(hit_df$row_idx), ]

  if (requireNamespace("org.Dm.eg.db",quietly=TRUE)) {
    suppressPackageStartupMessages(library(org.Dm.eg.db))
    sym_map <- tryCatch(
      AnnotationDbi::select(org.Dm.eg.db, keys=unique(hit_df$entrez),
                            columns="SYMBOL", keytype="ENTREZID"),
      error=function(e) NULL)
    if (is.null(sym_map)||!any(!is.na(sym_map$SYMBOL)))
      sym_map <- tryCatch(
        AnnotationDbi::select(org.Dm.eg.db, keys=unique(hit_df$entrez),
                              columns="SYMBOL", keytype="FLYBASE"),
        error=function(e) NULL)
    if (!is.null(sym_map) && "SYMBOL" %in% colnames(sym_map)) {
      kc <- intersect(c("ENTREZID","FLYBASE"), colnames(sym_map))[1]
      hit_df <- merge(hit_df, sym_map[,c(kc,"SYMBOL"),drop=FALSE],
                      by.x="entrez", by.y=kc, all.x=TRUE)
      cat(sprintf("[v3]   TxDb annotation: %d/%d peaks → genes (TSS±%d bp)\n",
                  nrow(hit_df), nrow(peak_df), tss_window))
      cat(sprintf("[v3]     Symbol sample: %s\n",
                  paste(head(na.omit(hit_df$SYMBOL),5), collapse=", ")))
    } else hit_df$SYMBOL <- hit_df$entrez
  } else hit_df$SYMBOL <- hit_df$entrez

  peak_df$symbol[hit_df$row_idx]          <- hit_df$SYMBOL
  peak_df$gene_id[hit_df$row_idx]         <- hit_df$entrez
  peak_df$peak_annotation[hit_df$row_idx] <- "TSS_proximal"
  peak_df$feature_class <- classify_feature(peak_df$peak_annotation)
  peak_df
}

join_to_rna <- function(peak_df, rna_df) {
  if (is.null(peak_df)||is.null(rna_df)) return(NULL)
  rna_s <- rna_df[,intersect(c("symbol","gene_id","log2FoldChange","t","pvalue","padj"),
                               colnames(rna_df)),drop=FALSE]
  names(rna_s)[names(rna_s)=="log2FoldChange"] <- "rna_lfc"
  names(rna_s)[names(rna_s)=="t"]              <- "rna_t"
  names(rna_s)[names(rna_s)=="pvalue"]         <- "rna_pval"
  names(rna_s)[names(rna_s)=="padj"]           <- "rna_padj"
  m <- NULL
  if ("symbol"%in%colnames(peak_df) && "symbol"%in%colnames(rna_s)) {
    m <- merge(peak_df[!is.na(peak_df$symbol)&peak_df$symbol!="",],
               rna_s[!is.na(rna_s$symbol)&rna_s$symbol!="",],
               by="symbol",all=FALSE)
  }
  if ((is.null(m)||nrow(m)==0) &&
      "gene_id"%in%colnames(peak_df) && "gene_id"%in%colnames(rna_s)) {
    m2 <- merge(peak_df[!is.na(peak_df$gene_id)&peak_df$gene_id!="",],
                rna_s[!is.na(rna_s$gene_id)&rna_s$gene_id!="",],
                by="gene_id",all=FALSE)
    if (nrow(m2)>0) { if (!"symbol"%in%colnames(m2)) m2$symbol<-m2$gene_id; m<-m2 }
  }
  if (is.null(m)||nrow(m)==0) return(NULL)
  m[!is.na(m$atac_lfc)&!is.na(m$rna_lfc),]
}

sig_stars <- function(p) {
  dplyr::case_when(p<0.001~"***",p<0.01~"**",p<0.05~"*",p<0.10~".",TRUE~"")
}

# =============================================================================
# PART 5b — EXPANDED BISTABILITY ANALYSIS
# =============================================================================

if ("5b" %in% run_parts) {
  cat("\n",strrep("=",70),"\n",sep="")
  cat("PART 5b: EXPANDED BISTABILITY ANALYSIS\n\n")

  # Load all temporal DAR sets
  dar_cache <- list()
  for (cfg in TEMPORAL_STEMS) {
    df <- read_bed_dar(cfg$stem, negate=cfg$negate)
    if (!is.null(df)) {
      cat(sprintf("[v3]   %s: %d peaks (negate=%s)\n",
                  cfg$stem, nrow(df), cfg$negate))
      dar_cache[[cfg$stem]] <- df
    }
  }

  # ── 5b-1: Bistable peak classification ─────────────────────────────────────
  df_botv  <- dar_cache[["BOTv_temporal"]]
  df_botcv <- dar_cache[["BOTCv_temporal"]]
  if (is.null(df_botv)||is.null(df_botcv)) {
    cat("[v3]   Temporal DARs not found — skipping 5b\n")
  } else {
    gr_botv  <- GRanges(seqnames=df_botv$chr,
                         ranges=IRanges(df_botv$start+1L,df_botv$end),
                         lfc=df_botv$atac_lfc, name=df_botv$name)
    gr_botcv <- GRanges(seqnames=df_botcv$chr,
                         ranges=IRanges(df_botcv$start+1L,df_botcv$end),
                         lfc=df_botcv$atac_lfc, name=df_botcv$name)
    hits  <- findOverlaps(gr_botv,gr_botcv,minoverlap=MIN_OVERLAP)
    ov_w  <- width(pintersect(gr_botv[queryHits(hits)],gr_botcv[subjectHits(hits)]))
    recip <- pmin(ov_w/width(gr_botv[queryHits(hits)]),
                  ov_w/width(gr_botcv[subjectHits(hits)]))
    best  <- data.frame(i_bv=queryHits(hits),i_bc=subjectHits(hits),recip=recip) %>%
      group_by(i_bv) %>% slice_max(recip,n=1,with_ties=FALSE) %>% ungroup()

    shared <- data.frame(
      chr=as.character(seqnames(gr_botv[best$i_bv])),
      start=start(gr_botv[best$i_bv])-1L,
      end=end(gr_botv[best$i_bv]),
      name=mcols(gr_botv[best$i_bv])$name,
      lfc_botv=mcols(gr_botv[best$i_bv])$lfc,
      lfc_botcv=mcols(gr_botcv[best$i_bc])$lfc,
      stringsAsFactors=FALSE) %>%
      filter(!is.na(lfc_botv),!is.na(lfc_botcv)) %>%
      mutate(
        is_bistable = abs(lfc_botv)>BISTABLE_THRESH &
          abs(lfc_botcv)>BISTABLE_THRESH &
          sign(lfc_botv)!=sign(lfc_botcv),
        is_concordant_open  = lfc_botv< -BISTABLE_THRESH & lfc_botcv< -BISTABLE_THRESH,
        is_concordant_close = lfc_botv>  BISTABLE_THRESH & lfc_botcv>  BISTABLE_THRESH,
        is_stable = abs(lfc_botv)<STAB_THRESH & abs(lfc_botcv)<STAB_THRESH,
        divergence_score = abs(lfc_botv)+abs(lfc_botcv),
        chrom_bias = lfc_botcv - lfc_botv,   # + = BOTCv retains/gains more
        bistable_class = dplyr::case_when(
          is_bistable & lfc_botv>0 ~ "BOTv closes / BOTCv opens",
          is_bistable & lfc_botv<0 ~ "BOTv opens / BOTCv closes",
          is_concordant_open       ~ "Both open (late)",
          is_concordant_close      ~ "Both close (late)",
          is_stable                ~ "Stable",
          TRUE                     ~ "Intermediate"))

    cat(sprintf("[v3]   Shared peaks: %d  |  bistable: %d\n",
                nrow(shared), sum(shared$is_bistable)))

    # Annotate all shared peaks — try BOTv_temporal ChIPseeker first
    shared_ann <- annotate_peaks(shared, stem_hint="BOTv_temporal")

    # ── Feature class breakdown for bistable peaks ────────────────────────────
    bistable_ann <- filter(shared_ann, is_bistable)
    if (nrow(bistable_ann) > 0 && any(!is.na(shared_ann$feature_class))) {
      feat_sum <- bistable_ann %>%
        mutate(feature_class=ifelse(is.na(feature_class),"Unannotated",feature_class)) %>%
        count(bistable_class, feature_class) %>%
        group_by(bistable_class) %>%
        mutate(pct=100*n/sum(n)) %>% ungroup()
      p_feat <- ggplot(feat_sum,
                        aes(x=bistable_class, y=pct, fill=feature_class)) +
        geom_col(position="stack", width=0.7) +
        geom_text(aes(label=ifelse(pct>=8,paste0(round(pct),"%"),"")),
                  position=position_stack(vjust=0.5),
                  size=2.8, colour="white", fontface="bold") +
        scale_fill_manual(values=c("Gene-proximal"="#e31a1c",
                                    "CRM-associated"="#1f78b4",
                                    "Unannotated"="#cccccc"),
                          name="Feature type") +
        scale_y_continuous(expand=expansion(mult=c(0,0.05))) +
        labs(title="Bistable peaks: genomic feature composition",
             x=NULL, y="% of bistable peaks") +
        theme_classic(base_size=11) +
        theme(plot.title=element_text(face="bold"),
              axis.text.x=element_text(angle=15,hjust=1))
      save_plot(file.path(OUT_DIR,"part5b_bistable_expanded","feature_class_breakdown"),
                p_feat, w=8, h=5)
    }

    # Print top annotated bistable peaks for quick inspection
    top_bistable <- bistable_ann %>%
      filter(!is.na(symbol), symbol!="") %>%
      arrange(desc(divergence_score)) %>%
      head(30) %>%
      dplyr::select(chr, start, end, symbol, peak_annotation, feature_class,
                    bistable_class, lfc_botv, lfc_botcv, divergence_score)
    cat(sprintf("[v3]   Top %d annotated bistable peaks:\n", nrow(top_bistable)))
    print(as.data.frame(top_bistable), row.names=FALSE)

    # ── 5b-2: Class summary bar ──────────────────────────────────────────────
    class_sum <- shared_ann %>% count(bistable_class) %>%
      mutate(pct=100*n/sum(n))
    class_cols <- c(
      "BOTv opens / BOTCv closes"  = COL_BOTV,
      "BOTv closes / BOTCv opens"  = COL_BOTCV,
      "Both open (late)"           = "#F39C12",
      "Both close (late)"          = "#8E44AD",
      "Stable"                     = "#AAAAAA",
      "Intermediate"               = "#DDDDDD")
    p_class <- ggplot(class_sum,
                       aes(x=reorder(bistable_class,-n), y=n, fill=bistable_class)) +
      geom_col(width=0.72) +
      geom_text(aes(label=sprintf("%d\n(%.0f%%)",n,pct)),
                vjust=-0.25, size=3.2) +
      scale_fill_manual(values=class_cols, guide="none") +
      scale_y_continuous(expand=expansion(mult=c(0,0.18))) +
      labs(title=sprintf("Chromatin trajectory classes\n(shared BOTv/BOTCv temporal peaks, |lFC|>%.1f)",
                          BISTABLE_THRESH),
           x=NULL, y="Peaks") +
      theme_classic(base_size=11) +
      theme(plot.title=element_text(face="bold"),
            axis.text.x=element_text(angle=15,hjust=1))

    save_plot(file.path(OUT_DIR,"part5b_bistable_expanded","class_summary"),
              p_class, w=9, h=5)

    # ── 5b-3: For each RNA contrast, full bistable concordance ───────────────
    for (rna_cfg in RNA_CONTRASTS) {
      rna_df <- read_rna(rna_cfg$name)
      if (is.null(rna_df)) next

      merged <- join_to_rna(shared_ann, rna_df)
      if (is.null(merged)||nrow(merged)<2) {
        cat(sprintf("[v3]   Skipping %s: only %d merged peaks\n",
                    rna_cfg$name, if(is.null(merged)) 0L else nrow(merged)))
        next
      }
      cat(sprintf("[v3]   %s: %d peaks with RNA data (%d bistable)\n",
                  rna_cfg$name, nrow(merged),
                  sum(merged$is_bistable, na.rm=TRUE)))

      merged <- merged %>% mutate(
        rna_sig = !is.na(rna_padj) & rna_padj < FDR_RNA,
        # Concordance: chromatin bias and RNA bias in same direction
        # chrom_bias > 0 = BOTCv retains/gains more accessibility over time
        # rna_lfc > 0    = higher in A arm (BOTCv or BOTCv_late depending on contrast)
        conc_class = dplyr::case_when(
          is_bistable & chrom_bias>0 & rna_lfc>0 ~ "Concordant: BOTCv-biased",
          is_bistable & chrom_bias<0 & rna_lfc<0 ~ "Concordant: BOTv-biased",
          is_bistable                              ~ "Discordant bistable",
          is_concordant_open                       ~ "Both open + RNA",
          TRUE                                     ~ bistable_class),
        label_use = dplyr::case_when(
          is_bistable & rna_sig ~ symbol,
          TRUE                  ~ NA_character_))

      bistable <- filter(merged, is_bistable)

      # ── Scatter: BOTv lfc vs BOTCv lfc, coloured by RNA lfc ────────────────
      p_scatter <- ggplot(merged, aes(x=lfc_botv, y=lfc_botcv)) +
        # Background all peaks
        geom_point(data=filter(merged,!is_bistable),
                   colour="#e0e0e0", size=0.5, alpha=0.3) +
        # Bistable coloured by RNA
        geom_point(data=bistable,
                   aes(colour=rna_lfc, size=divergence_score,
                       shape=rna_sig), alpha=0.85) +
        scale_colour_gradient2(low=COL_BOTV, mid="#f5f5f5", high=COL_BOTCV,
                               midpoint=0, na.value="#cccccc",
                               name=sprintf("RNA lFC\n(%s+)",rna_cfg$positive_is)) +
        scale_size_continuous(range=c(1.5,5), name="ATAC\ndivergence",
                              guide="none") +
        scale_shape_manual(values=c("TRUE"=16,"FALSE"=1),
                           labels=c("TRUE"=sprintf("RNA padj<%.2f",FDR_RNA),
                                    "FALSE"="not sig"),
                           name="RNA") +
        geom_hline(yintercept=0,lty=2,colour="#cccccc",linewidth=0.4) +
        geom_vline(xintercept=0,lty=2,colour="#cccccc",linewidth=0.4) +
        geom_abline(slope=-1,intercept=0,lty=3,colour="#aaa",linewidth=0.4) +
        {
          if (has_ggrepel && nrow(bistable)>0) {
            lab_df <- bistable %>% filter(!is.na(symbol),symbol!="",rna_sig) %>%
              arrange(desc(abs(rna_lfc))) %>% head(TOP_N)
            if (nrow(lab_df)>0)
              geom_text_repel(data=lab_df,
                              aes(x=lfc_botv, y=lfc_botcv, label=symbol),
                              colour="black", size=2.5, fontface="italic",
                              max.overlaps=20, box.padding=0.4,
                              segment.size=0.3, seed=42, inherit.aes=FALSE)
          }
        } +
        annotate("text",
                 x=max(merged$lfc_botv,na.rm=TRUE)*0.75,
                 y=min(merged$lfc_botcv,na.rm=TRUE)*0.75,
                 label="BOTv opens\nBOTCv closes",
                 size=2.8, colour=COL_BOTV_LATE, fontface="bold") +
        annotate("text",
                 x=min(merged$lfc_botv,na.rm=TRUE)*0.75,
                 y=max(merged$lfc_botcv,na.rm=TRUE)*0.75,
                 label="BOTCv opens\nBOTv closes",
                 size=2.8, colour=COL_BOTCV_LATE, fontface="bold") +
        labs(title=sprintf("Chromatin bistability × Expression\n%s",rna_cfg$label),
             subtitle=sprintf("%d bistable peaks  |  Filled = RNA padj<%.2f  |  Dark pink=BOTCv  Dark green=BOTv",
                              nrow(bistable),FDR_RNA),
             x="BOTv temporal lFC (positive = nc14b open → closing over time)",
             y="BOTCv temporal lFC (positive = nc14b open → closing over time)") +
        theme_classic(base_size=11) +
        theme(plot.title=element_text(face="bold",size=10.5))

      # ── Chromatin-bias vs RNA-lfc scatter (bistable only) ──────────────────
      if (nrow(bistable)>3) {
        p_bias <- ggplot(bistable,
                          aes(x=chrom_bias, y=rna_lfc, colour=conc_class,
                              size=divergence_score)) +
          geom_hline(yintercept=0,lty=2,colour="#cccccc") +
          geom_vline(xintercept=0,lty=2,colour="#cccccc") +
          geom_point(alpha=0.8) +
          scale_colour_manual(values=c(
            "Concordant: BOTCv-biased" = COL_BOTCV_LATE,
            "Concordant: BOTv-biased"  = COL_BOTV_LATE,
            "Discordant bistable"      = "#aaaaaa",
            "Both open + RNA"          = "#F39C12"),
            name="Class") +
          scale_size_continuous(range=c(1.5,4), guide="none") +
          {
            if (has_ggrepel) {
              lab2 <- bistable %>%
                filter(!is.na(symbol),symbol!="",rna_sig,
                       grepl("Concordant",conc_class)) %>%
                arrange(desc(abs(rna_lfc)*divergence_score)) %>% head(15)
              if (nrow(lab2)>0)
                geom_text_repel(data=lab2,
                                aes(x=chrom_bias,y=rna_lfc,label=symbol),
                                colour="black",size=2.5,fontface="italic",
                                max.overlaps=15,box.padding=0.4,seed=42,
                                inherit.aes=FALSE)
            }
          } +
          geom_smooth(method="lm",se=TRUE,colour="#333",linewidth=0.6,
                      linetype="dashed",alpha=0.12) +
          labs(title=sprintf("Chromatin bias vs expression\n%s",rna_cfg$label),
               subtitle="chrom_bias = lfc_BOTCv − lfc_BOTv  |  + = BOTCv retains open longer",
               x="Chromatin bias (BOTCv − BOTv temporal lFC)",
               y=sprintf("RNA lFC (%s)",rna_cfg$positive_is)) +
          theme_classic(base_size=11) +
          theme(plot.title=element_text(face="bold",size=10.5))

        # Pearson correlation
        vr <- !is.na(bistable$chrom_bias) & !is.na(bistable$rna_lfc)
        if (sum(vr)>=5) {
          cr <- cor.test(bistable$chrom_bias[vr], bistable$rna_lfc[vr])
          cat(sprintf("[v3]   %s: chrom_bias×RNA r=%.3f, p=%.4f (n=%d)\n",
                      rna_cfg$name, cr$estimate, cr$p.value, sum(vr)))
        }
      } else p_bias <- NULL

      # ── Top concordant genes — lollipop with labels fixed ──────────────────
      top_conc <- bistable %>%
        filter(grepl("Concordant",conc_class),
               !is.na(symbol), symbol!="", !is.na(rna_lfc)) %>%
        mutate(combined_score = abs(rna_lfc)*divergence_score,
               atac_dir = ifelse(chrom_bias>0,"BOTCv-biased","BOTv-biased")) %>%
        arrange(desc(combined_score)) %>% head(TOP_N)

      if (nrow(top_conc)>0) {
        p_top <- ggplot(top_conc,
                         aes(x=reorder(symbol,combined_score),
                             y=combined_score, fill=atac_dir)) +
          geom_col(width=0.75) +
          geom_text(aes(label=sprintf("ATAC:%.2f/%.2f RNA:%.2f %s",
                                       lfc_botv,lfc_botcv,rna_lfc,
                                       sig_stars(rna_padj))),
                    hjust=-0.05, size=2.5, fontface="plain") +
          scale_fill_manual(values=c("BOTCv-biased"=COL_BOTCV_LATE,"BOTv-biased"=COL_BOTV_LATE),
                            name="Chromatin bias") +
          scale_y_continuous(expand=expansion(mult=c(0.05,0.55))) +
          coord_flip(clip="off") +
          labs(title=sprintf("Top concordant bistable genes\n%s",rna_cfg$label),
               subtitle="|RNA lFC| × chromatin divergence score",
               x=NULL, y="Combined score") +
          theme_classic(base_size=10) +
          theme(plot.title=element_text(face="bold"),
                plot.margin=margin(5,20,5,40,"pt"))
      } else p_top <- NULL

      # ── Concordance class bar ────────────────────────────────────────────────
      conc_sum <- bistable %>%
        group_by(conc_class) %>%
        summarise(n=n(), n_sig=sum(rna_sig,na.rm=TRUE), .groups="drop") %>%
        mutate(pct=100*n/sum(n))
      p_conc <- ggplot(conc_sum,
                        aes(x=reorder(conc_class,-n), y=n, fill=conc_class)) +
        geom_col(width=0.7) +
        geom_text(aes(label=sprintf("%d (%.0f%%)\n%d RNA sig",n,pct,n_sig)),
                  vjust=-0.2, size=2.8) +
        scale_fill_manual(values=c(
          "Concordant: BOTCv-biased"=COL_BOTCV_LATE,
          "Concordant: BOTv-biased" =COL_BOTV_LATE,
          "Discordant bistable"     ="#aaaaaa",
          "Both open + RNA"         ="#F39C12"),
          guide="none") +
        scale_y_continuous(expand=expansion(mult=c(0,0.2))) +
        labs(title="Concordance class", x=NULL, y="Bistable peaks") +
        theme_classic(base_size=10) +
        theme(axis.text.x=element_text(angle=15,hjust=1))

      # Full combined figure
      if (!is.null(p_bias) && !is.null(p_top)) {
        p_full <- (p_scatter | p_conc) / (p_bias | p_top) +
          plot_layout(heights=c(1,1))
      } else if (!is.null(p_bias)) {
        p_full <- (p_scatter | p_conc) / p_bias
      } else {
        p_full <- p_scatter | p_conc
      }

      stem_out <- file.path(OUT_DIR,"part5b_bistable_expanded",
                            sprintf("bistable_%s",rna_cfg$name))
      save_plot(stem_out, p_full, w=16, h=12)

      # Save full table
      write.table(merged,
                  paste0(stem_out,"_table.tsv"),
                  sep="\t",quote=FALSE,row.names=FALSE)
    }

    # ── 5b-4: Multi-stem comparison: which bistable peaks are consistent? ────
    # For each cross-genotype stem (nc14b and nc14late), flag whether each
    # bistable peak also shows differential accessibility at that timepoint
    if (!is.null(dar_cache[["BOTv_vs_BOTCv"]]) &&
        !is.null(dar_cache[["BOTv_vs_BOTCv_nc14late"]])) {
      cat("[v3]   Cross-referencing bistable peaks with cross-genotype DARs...\n")
      gr_nc14b   <- with(dar_cache[["BOTv_vs_BOTCv"]],
                          GRanges(chr,IRanges(start+1L,end),
                                  lfc=atac_lfc,name=name))
      gr_nc14late<- with(dar_cache[["BOTv_vs_BOTCv_nc14late"]],
                          GRanges(chr,IRanges(start+1L,end),
                                  lfc=atac_lfc,name=name))

      bistable_ann <- filter(shared_ann, is_bistable)
      if (nrow(bistable_ann)>0) {
        gr_bist <- GRanges(bistable_ann$chr,
                            IRanges(bistable_ann$start+1L, bistable_ann$end))
        bistable_ann$diff_at_nc14b    <- countOverlaps(gr_bist,gr_nc14b)>0
        bistable_ann$diff_at_nc14late <- countOverlaps(gr_bist,gr_nc14late)>0
        bistable_ann$temporal_profile <- dplyr::case_when(
          bistable_ann$diff_at_nc14b & bistable_ann$diff_at_nc14late ~
            "Bistable + cross-geno at both TPs",
          bistable_ann$diff_at_nc14b & !bistable_ann$diff_at_nc14late ~
            "Bistable + cross-geno at nc14b only",
          !bistable_ann$diff_at_nc14b & bistable_ann$diff_at_nc14late ~
            "Bistable + cross-geno at nc14late only",
          TRUE ~ "Bistable temporal only")
        tp_sum <- bistable_ann %>% count(temporal_profile) %>%
          mutate(pct=100*n/sum(n))
        cat("[v3]   Bistable temporal profile:\n")
        print(tp_sum)
        p_tp <- ggplot(tp_sum,
                        aes(x=reorder(temporal_profile,-n), y=n,
                            fill=temporal_profile)) +
          geom_col(width=0.7) +
          geom_text(aes(label=sprintf("%d (%.0f%%)",n,pct)),vjust=-0.2,size=3) +
          scale_fill_brewer(palette="Set2",guide="none") +
          scale_y_continuous(expand=expansion(mult=c(0,0.18))) +
          labs(title="Bistable peaks: cross-genotype DAR overlap at nc14b vs nc14late",
               x=NULL, y="Bistable peaks") +
          theme_classic(base_size=11) +
          theme(plot.title=element_text(face="bold"),
                axis.text.x=element_text(angle=15,hjust=1,size=9))
        save_plot(file.path(OUT_DIR,"part5b_bistable_expanded",
                            "bistable_temporal_profile"),
                  p_tp, w=10, h=5)
        write.table(bistable_ann,
                    file.path(OUT_DIR,"part5b_bistable_expanded",
                              "bistable_annotated_full.tsv"),
                    sep="\t",quote=FALSE,row.names=FALSE)
      }
    }
  }
}

# =============================================================================
# PART 6 — TEMPORAL LEAD/LAG ANALYSIS
# =============================================================================

if ("6" %in% run_parts) {
  cat("\n",strrep("=",70),"\n",sep="")
  cat("PART 6: TEMPORAL LEAD/LAG — does chromatin precede expression?\n\n")

  # Need: ATAC at nc14b AND nc14late cross-genotype, plus RNA at nc14b AND late
  df_atac_b    <- dar_cache[["BOTv_vs_BOTCv"]]
  df_atac_late <- dar_cache[["BOTv_vs_BOTCv_nc14late"]]
  rna_nc14b    <- read_rna("BOTCv_nc14b_vs_BOTv_nc14b")
  rna_late     <- read_rna("BOTCv_late_vs_BOTv_nc14d")

  if (is.null(df_atac_b)||is.null(df_atac_late)||
      is.null(rna_nc14b)||is.null(rna_late)) {
    cat("[v3]   Missing data for Part 6 — need ATAC nc14b+nc14late and RNA nc14b+late\n")
    if (!exists("dar_cache") || is.null(dar_cache[["BOTv_vs_BOTCv"]]))
      cat("[v3]   Note: run Part 5b first to load DARs, or check DAR_DIR\n")
  } else {
    # Annotate both ATAC sets
    cat("[v3]   Annotating nc14b ATAC peaks...\n")
    df_atac_b_ann    <- annotate_peaks(df_atac_b,    stem_hint="BOTv_vs_BOTCv")
    cat("[v3]   Annotating nc14late ATAC peaks...\n")
    df_atac_late_ann <- annotate_peaks(df_atac_late, stem_hint="BOTv_vs_BOTCv_nc14late")

    # Build gene-level ATAC summary: per gene, what is the ATAC lFC at each TP?
    # Use the peak with highest |lfc| if a gene has multiple peaks
    atac_gene_b <- df_atac_b_ann %>%
      filter(!is.na(symbol),symbol!="") %>%
      group_by(symbol) %>%
      slice_max(abs(atac_lfc),n=1,with_ties=FALSE) %>%
      ungroup() %>%
      dplyr::select(symbol,atac_lfc_nc14b=atac_lfc,
                    atac_padj_nc14b=atac_padj)

    atac_gene_late <- df_atac_late_ann %>%
      filter(!is.na(symbol),symbol!="") %>%
      group_by(symbol) %>%
      slice_max(abs(atac_lfc),n=1,with_ties=FALSE) %>%
      ungroup() %>%
      dplyr::select(symbol,atac_lfc_nc14late=atac_lfc,
                    atac_padj_nc14late=atac_padj)

    # RNA gene-level
    rna_gene_b <- rna_nc14b %>%
      dplyr::select(symbol,rna_lfc_nc14b=log2FoldChange,
                    rna_t_nc14b=t, rna_padj_nc14b=padj)
    rna_gene_late <- rna_late %>%
      dplyr::select(symbol,rna_lfc_nc14late=log2FoldChange,
                    rna_t_nc14late=t, rna_padj_nc14late=padj)

    # Merge all four on gene symbol
    combined <- Reduce(function(a,b) merge(a,b,by="symbol",all=FALSE),
                       list(atac_gene_b, atac_gene_late,
                            rna_gene_b,  rna_gene_late))
    cat(sprintf("[v3]   Lead/lag gene universe: %d genes with ATAC+RNA at both TPs\n",
                nrow(combined)))

    # Classify temporal relationship
    # "sig" for ATAC: padj < FDR_RNA (or if padj NA, use |lfc| > BISTABLE_THRESH)
    # "sig" for RNA:  padj < FDR_RNA
    combined <- combined %>% mutate(
      atac_sig_b    = (!is.na(atac_padj_nc14b)    & atac_padj_nc14b    < FDR_RNA) |
                      (is.na(atac_padj_nc14b)      & abs(atac_lfc_nc14b)   > BISTABLE_THRESH),
      atac_sig_late = (!is.na(atac_padj_nc14late)  & atac_padj_nc14late  < FDR_RNA) |
                      (is.na(atac_padj_nc14late)    & abs(atac_lfc_nc14late) > BISTABLE_THRESH),
      rna_sig_b     = !is.na(rna_padj_nc14b)    & rna_padj_nc14b    < FDR_RNA,
      rna_sig_late  = !is.na(rna_padj_nc14late)  & rna_padj_nc14late  < FDR_RNA,
      # [v3 fix] Directional concordance, checked against the SPECIFIC
      # timepoint pairing each lead/lag claim actually makes -- not one
      # blanket "late vs late" comparison for every class. "Chromatin-first"
      # claims that the EARLY ATAC signal predicts (same direction as) the
      # LATE RNA signal, so that pairing is atac_nc14b vs rna_nc14late, not
      # atac_nc14late vs rna_nc14late. The previous single concordant_direction
      # column used the wrong pairing for every class except Co-emerging, and
      # -- more importantly -- was never actually wired into the lead_lag
      # case_when below at all, so classification ran on significance
      # timing alone with no direction check whatsoever. That let genes
      # with flatly OPPOSITE-sign ATAC and RNA changes (chromatin going
      # BOTv-biased while expression simultaneously went MORE BOTCv-biased,
      # or vice versa) get labelled "Chromatin-first" purely because the
      # significance timing happened to match, which is what produced
      # "opposite calls to the ATAC" in the trajectory heatmap.
      concordant_stable      = sign(atac_lfc_nc14b)    == sign(rna_lfc_nc14b),
      concordant_chrom_first = sign(atac_lfc_nc14b)    == sign(rna_lfc_nc14late),
      concordant_expr_first  = sign(rna_lfc_nc14b)     == sign(atac_lfc_nc14late),
      concordant_co_emerge   = sign(atac_lfc_nc14late) == sign(rna_lfc_nc14late),
      # Lead/lag classification -- now gated on direction, not just timing.
      lead_lag = dplyr::case_when(
        # Both already different at nc14b, same direction
        atac_sig_b    & rna_sig_b    &  concordant_stable ~
          "Stable difference\n(both sig at nc14b)",
        # ATAC first at nc14b, RNA not sig until late, SAME direction --
        # this is the only pattern that actually supports "chromatin leads"
        atac_sig_b    & !rna_sig_b   & rna_sig_late  &  concordant_chrom_first ~
          "Chromatin-first\n(ATAC nc14b, RNA late)",
        # RNA first at nc14b, ATAC not sig until late, SAME direction
        !atac_sig_b   & rna_sig_b    & atac_sig_late &  concordant_expr_first ~
          "Expression-first\n(RNA nc14b, ATAC late)",
        # Both emerge at nc14late, SAME direction
        !atac_sig_b   & !rna_sig_b   & atac_sig_late & rna_sig_late & concordant_co_emerge ~
          "Co-emerging\n(both at nc14late)",
        # Same significance-timing pattern as above, but OPPOSITE direction:
        # a real, distinct biological pattern (chromatin and transcription
        # moving in conflicting directions -- e.g. a compensatory/buffering
        # response) rather than a lead/lag relationship. Kept as its own
        # set of classes instead of silently folding into the concordant
        # ones above, or dropping the genes entirely.
        atac_sig_b    & rna_sig_b    & !concordant_stable ~
          "Discordant\n(opposite direction, both sig nc14b)",
        atac_sig_b    & !rna_sig_b   & rna_sig_late  & !concordant_chrom_first ~
          "Discordant\n(chromatin leads, opposite direction)",
        !atac_sig_b   & rna_sig_b    & atac_sig_late & !concordant_expr_first ~
          "Discordant\n(expression leads, opposite direction)",
        !atac_sig_b   & !rna_sig_b   & atac_sig_late & rna_sig_late & !concordant_co_emerge ~
          "Discordant\n(co-emerging, opposite direction)",
        # Only ATAC changes, no RNA sig at either TP
        atac_sig_b    & !rna_sig_b   & !rna_sig_late ~ "ATAC-only\n(no RNA sig)",
        atac_sig_late & !rna_sig_b   & !rna_sig_late ~ "ATAC-only\n(no RNA sig)",
        # Only RNA changes
        rna_sig_b     & !atac_sig_b  & !atac_sig_late ~ "RNA-only\n(no ATAC sig)",
        rna_sig_late  & !atac_sig_b  & !atac_sig_late ~ "RNA-only\n(no ATAC sig)",
        # Neither sig
        TRUE ~ "Neither sig"),
      combined_effect = abs(atac_lfc_nc14late) + abs(rna_lfc_nc14late))

    # ── 6a: Lead/lag class bar ──────────────────────────────────────────────
    ll_sum <- combined %>%
      count(lead_lag) %>%
      mutate(pct=100*n/sum(n)) %>%
      filter(!grepl("Neither",lead_lag))   # focus on genes with some signal

    lead_lag_cols <- c(
      "Stable difference\n(both sig at nc14b)"         = "#c0392b",
      "Chromatin-first\n(ATAC nc14b, RNA late)"        = "#9b59b6",
      "Expression-first\n(RNA nc14b, ATAC late)"       = "#1abc9c",
      "Co-emerging\n(both at nc14late)"                = "#f39c12",
      "ATAC-only\n(no RNA sig)"                        = "#bdc3c7",
      "RNA-only\n(no ATAC sig)"                        = "#ecf0f1",
      # [v3 fix] same significance-timing pattern as the four classes above,
      # but opposite-direction -- muted/desaturated variants so they read
      # as "related but not the same claim" rather than competing for the
      # same visual weight as a real lead/lag call.
      "Discordant\n(opposite direction, both sig nc14b)"    = "#7f1d13",
      "Discordant\n(chromatin leads, opposite direction)"   = "#5b3568",
      "Discordant\n(expression leads, opposite direction)"  = "#0e6b5c",
      "Discordant\n(co-emerging, opposite direction)"       = "#a66a08")

    p_ll_bar <- ggplot(ll_sum,
                        aes(x=reorder(lead_lag,-n), y=n, fill=lead_lag)) +
      geom_col(width=0.7) +
      geom_text(aes(label=sprintf("%d\n(%.0f%%)",n,pct)),
                vjust=-0.2, size=3.2) +
      scale_fill_manual(values=lead_lag_cols, guide="none") +
      scale_y_continuous(expand=expansion(mult=c(0,0.18))) +
      labs(title="Temporal lead/lag: chromatin vs expression\n(BOTCv vs BOTv)",
           subtitle=sprintf("Genes with ATAC+RNA at both nc14b and nc14late  |  n=%d",
                            nrow(combined)),
           x=NULL, y="Genes") +
      theme_classic(base_size=11) +
      theme(plot.title=element_text(face="bold"),
            axis.text.x=element_text(angle=20,hjust=1,size=9))
    save_plot(file.path(OUT_DIR,"part6_temporal_lead_lag","lead_lag_summary"),
              p_ll_bar, w=11, h=6)

    # ── 6b: ATAC lFC (nc14b) vs RNA lFC (nc14late) scatter ───────────────────
    # Chromatin-first genes should have large ATAC lFC at nc14b
    # but small RNA lFC at nc14b, then RNA lFC catches up by nc14late
    focus <- filter(combined, !grepl("Neither",lead_lag))
    p_atac_early_rna_late <- ggplot(focus,
                                     aes(x=atac_lfc_nc14b, y=rna_lfc_nc14late,
                                         colour=lead_lag)) +
      geom_hline(yintercept=0,lty=2,colour="#cccccc") +
      geom_vline(xintercept=0,lty=2,colour="#cccccc") +
      geom_point(size=1.2, alpha=0.6) +
      scale_colour_manual(values=lead_lag_cols, name="Class") +
      {
        if (has_ggrepel) {
          top_cf <- focus %>%
            filter(grepl("Chromatin-first",lead_lag),!is.na(symbol)) %>%
            arrange(desc(abs(atac_lfc_nc14b)*abs(rna_lfc_nc14late))) %>% head(12)
          if (nrow(top_cf)>0)
            geom_text_repel(data=top_cf,
                            aes(x=atac_lfc_nc14b,y=rna_lfc_nc14late,label=symbol),
                            colour="black",size=2.5,fontface="italic",
                            max.overlaps=15,seed=42,inherit.aes=FALSE)
        }
      } +
      labs(title="Chromatin (nc14b) vs Expression (nc14late)\nBOTCv vs BOTv",
           subtitle="Chromatin-first genes: large |ATAC lFC| at nc14b, RNA difference emerges by nc14late",
           x="ATAC log2FC at nc14b (BOTCv vs BOTv)",
           y="RNA log2FC at nc14late (BOTCv late vs BOTv nc14d)") +
      theme_classic(base_size=11) +
      theme(plot.title=element_text(face="bold",size=10.5))
    save_plot(file.path(OUT_DIR,"part6_temporal_lead_lag","atac_nc14b_vs_rna_nc14late"),
              p_atac_early_rna_late, w=9, h=7)

    # ── 6c: Temporal trajectory heatmap ──────────────────────────────────────
    # For chromatin-first + co-emerging genes: show their ATAC and RNA lFC
    # at both timepoints as a 4-column heatmap
    traj_genes <- focus %>%
      filter(grepl("Chromatin-first|Co-emerging|Stable",lead_lag),
             !is.na(symbol),symbol!="") %>%
      arrange(desc(combined_effect)) %>%
      head(50) %>%
      dplyr::select(symbol, lead_lag,
                    `ATAC nc14b`=atac_lfc_nc14b,
                    `RNA nc14b`=rna_lfc_nc14b,
                    `ATAC late`=atac_lfc_nc14late,
                    `RNA late`=rna_lfc_nc14late) %>%
      pivot_longer(cols=c(`ATAC nc14b`,`RNA nc14b`,`ATAC late`,`RNA late`),
                   names_to="measurement", values_to="lfc") %>%
      mutate(measurement=factor(measurement,
                                 levels=c("ATAC nc14b","RNA nc14b",
                                          "ATAC late","RNA late")))

    if (nrow(traj_genes)>0) {
      p_traj <- ggplot(traj_genes,
                        aes(x=measurement,
                            y=reorder(symbol,lfc),
                            fill=lfc)) +
        geom_tile(colour="white",linewidth=0.4) +
        facet_grid(lead_lag~., scales="free_y", space="free_y") +
        scale_fill_gradient2(low=COL_BOTV, mid="white", high=COL_BOTCV,
                             midpoint=0, name="log2FC\n(BOTCv/+ arm)") +
        scale_x_discrete(expand=c(0,0)) +
        labs(title="Temporal trajectory: ATAC and RNA lFC at nc14b and nc14late",
             subtitle="Top 50 genes with chromatin-first, co-emerging, or stable signal",
             x=NULL, y="Gene") +
        theme_classic(base_size=9) +
        theme(strip.text.y=element_text(size=7,angle=0),
              axis.text.y=element_text(size=7, face="italic"),
              plot.title=element_text(face="bold",size=10))
      save_plot(file.path(OUT_DIR,"part6_temporal_lead_lag","trajectory_heatmap"),
                p_traj, w=8, h=max(5,nrow(traj_genes)/4*0.3+3))
    }

    # ── 6d: Per-class top gene lollipop ──────────────────────────────────────
    lead_lag_key_classes <- c("Chromatin-first\n(ATAC nc14b, RNA late)",
                               "Co-emerging\n(both at nc14late)",
                               "Stable difference\n(both sig at nc14b)")
    for (cls in lead_lag_key_classes) {
      cls_genes <- combined %>%
        # [v3 fix] no longer needs a concordant_direction filter here --
        # lead_lag itself is now concordance-gated at classification time
        # (see the case_when above), so every row already in one of these
        # classes is already direction-consistent by construction.
        filter(lead_lag==cls, !is.na(symbol), symbol!="") %>%
        arrange(desc(combined_effect)) %>% head(TOP_N)
      if (nrow(cls_genes)==0) next
      cls_name <- gsub("\n","_",gsub("[^A-Za-z_]","",cls))
      p_cls <- ggplot(cls_genes,
                       aes(x=reorder(symbol,combined_effect),
                           y=combined_effect,
                           colour=rna_lfc_nc14late>0)) +
        geom_segment(aes(xend=reorder(symbol,combined_effect), yend=0),
                     linewidth=1.4) +
        geom_point(size=3.5) +
        geom_text(aes(label=sprintf("A:%.2f/%.2f R:%.2f",
                                     atac_lfc_nc14b, atac_lfc_nc14late,
                                     rna_lfc_nc14late)),
                  hjust=-0.1, size=2.5) +
        scale_colour_manual(values=c("TRUE"=COL_BOTCV,"FALSE"=COL_BOTV),
                            labels=c("TRUE"="BOTCv-biased","FALSE"="BOTv-biased"),
                            name=NULL) +
        scale_y_continuous(expand=expansion(mult=c(0.05,0.55))) +
        coord_flip(clip="off") +
        labs(title=sprintf("Class: %s", gsub("\n"," ",cls)),
             subtitle="A=ATAC lFC (nc14b/late), R=RNA lFC (late)",
             x=NULL, y="|ATAC late lFC| + |RNA late lFC|") +
        theme_classic(base_size=10) +
        theme(plot.title=element_text(face="bold"),
              plot.margin=margin(5,20,5,40,"pt"))
      save_plot(file.path(OUT_DIR,"part6_temporal_lead_lag",
                          sprintf("top_genes_%s",cls_name)),
                p_cls, w=11, h=max(4,nrow(cls_genes)*0.42+2))
    }

    # Save full table
    write.table(combined,
                file.path(OUT_DIR,"part6_temporal_lead_lag","lead_lag_full_table.tsv"),
                sep="\t",quote=FALSE,row.names=FALSE)
    cat(sprintf("[v3]   Lead/lag table saved: %d genes\n",nrow(combined)))

    # ── Side-by-side nc14b vs nc14late ATAC scatter, temporally coloured ────
    # Two scatter panels: left = nc14b ATAC vs RNA nc14b, right = nc14late
    # ATAC vs RNA nc14late. Points are coloured by lead/lag class so you can
    # see whether the class identity shifts between timepoints.
    valid_both <- combined %>%
      filter(!is.na(atac_lfc_nc14b), !is.na(atac_lfc_nc14late),
             !is.na(rna_lfc_nc14b),  !is.na(rna_lfc_nc14late),
             !grepl("Neither", lead_lag))

    lead_lag_cols_v3 <- c(
      "Stable difference\n(both sig at nc14b)"        = COL_CROSS,
      "Chromatin-first\n(ATAC nc14b, RNA late)"       = COL_BOTCV_LATE,
      "Expression-first\n(RNA nc14b, ATAC late)"      = COL_BOTV_LATE,
      "Co-emerging\n(both at nc14late)"               = "#f39c12",
      "ATAC-only\n(no RNA sig)"                       = "#bdc3c7",
      "RNA-only\n(no ATAC sig)"                       = "#ecf0f1",
      "Discordant\n(opposite direction, both sig nc14b)"    = "#7f1d13",
      "Discordant\n(chromatin leads, opposite direction)"   = "#5b3568",
      "Discordant\n(expression leads, opposite direction)"  = "#0e6b5c",
      "Discordant\n(co-emerging, opposite direction)"       = "#a66a08")

    if (nrow(valid_both) >= 5) {
      p_nc14b_scatter <- ggplot(valid_both,
                                 aes(x=atac_lfc_nc14b, y=rna_lfc_nc14b,
                                     colour=lead_lag)) +
        geom_hline(yintercept=0,lty=2,colour="#dddddd") +
        geom_vline(xintercept=0,lty=2,colour="#dddddd") +
        geom_point(size=1.5, alpha=0.75) +
        scale_colour_manual(values=lead_lag_cols_v3, name="Lead/lag class") +
        labs(title="nc14b — ATAC vs RNA",
             subtitle="Colour = temporal class",
             x="ATAC log2FC (BOTCv vs BOTv, nc14b)",
             y="RNA log2FC (BOTCv vs BOTv, nc14b)") +
        theme_classic(base_size=10) +
        theme(plot.title=element_text(face="bold"))

      p_nc14late_scatter <- ggplot(valid_both,
                                    aes(x=atac_lfc_nc14late, y=rna_lfc_nc14late,
                                        colour=lead_lag)) +
        geom_hline(yintercept=0,lty=2,colour="#dddddd") +
        geom_vline(xintercept=0,lty=2,colour="#dddddd") +
        geom_point(size=1.5, alpha=0.75) +
        scale_colour_manual(values=lead_lag_cols_v3, name="Lead/lag class") +
        {
          if (has_ggrepel) {
            top_ll <- valid_both %>%
              filter(grepl("Chromatin-first",lead_lag),!is.na(symbol)) %>%
              arrange(desc(abs(atac_lfc_nc14b)*abs(rna_lfc_nc14late))) %>%
              head(12)
            if (nrow(top_ll)>0)
              geom_text_repel(data=top_ll,
                              aes(x=atac_lfc_nc14late, y=rna_lfc_nc14late,
                                  label=symbol),
                              colour="black",size=2.4,fontface="italic",
                              max.overlaps=12,seed=42,inherit.aes=FALSE)
          }
        } +
        labs(title="nc14late — ATAC vs RNA",
             subtitle="Chromatin-first genes labelled",
             x="ATAC log2FC (BOTCv vs BOTv, nc14late)",
             y="RNA log2FC (BOTCv late vs BOTv nc14d)") +
        theme_classic(base_size=10) +
        theme(plot.title=element_text(face="bold"))

      p_side <- (p_nc14b_scatter | p_nc14late_scatter) +
        plot_layout(guides="collect") &
        theme(legend.position="bottom")
      save_plot(file.path(OUT_DIR,"part6_temporal_lead_lag",
                          "nc14b_vs_nc14late_scatter"),
                p_side, w=14, h=6)
    }

    # ── ATAC change vs RNA change across timepoints — temporal arrows ────────
    # Arrow plot: each gene is drawn as an arrow from its nc14b position
    # to its nc14late position in ATAC lFC × RNA lFC space.
    # Genes that "travel" far in both axes are the lead/lag candidates.
    chrom_first <- combined %>%
      filter(grepl("Chromatin-first",lead_lag), !is.na(symbol), symbol!="") %>%
      arrange(desc(abs(atac_lfc_nc14b)*abs(rna_lfc_nc14late))) %>%
      head(20)

    if (nrow(chrom_first) >= 3) {
      arrow_df <- chrom_first %>%
        dplyr::select(symbol,
                       x0=atac_lfc_nc14b, y0=rna_lfc_nc14b,
                       x1=atac_lfc_nc14late, y1=rna_lfc_nc14late)
      p_arrow <- ggplot(arrow_df) +
        geom_hline(yintercept=0,lty=2,colour="#dddddd") +
        geom_vline(xintercept=0,lty=2,colour="#dddddd") +
        geom_segment(aes(x=x0,y=y0,xend=x1,yend=y1),
                     arrow=arrow(length=unit(0.12,"cm"),type="closed"),
                     colour=COL_BOTCV_LATE, linewidth=0.8, alpha=0.8) +
        geom_point(aes(x=x0,y=y0), colour=COL_BOTCV_EARLY, size=2) +
        geom_point(aes(x=x1,y=y1), colour=COL_BOTCV_LATE,  size=2) +
        {
          if (has_ggrepel)
            geom_text_repel(aes(x=x1,y=y1,label=symbol),
                            colour="black",size=2.5,fontface="italic",
                            max.overlaps=15,seed=42)
          else
            geom_text(aes(x=x1,y=y1,label=symbol),hjust=-0.1,size=2.5)
        } +
        labs(title="Chromatin-first genes: temporal trajectory",
             subtitle=sprintf("Arrow = nc14b (light, %s) to nc14late (dark, %s)\nn=%d genes",
                              COL_BOTCV_EARLY, COL_BOTCV_LATE, nrow(chrom_first)),
             x="ATAC log2FC (BOTCv vs BOTv)",
             y="RNA log2FC") +
        theme_classic(base_size=11) +
        theme(plot.title=element_text(face="bold"))
      save_plot(file.path(OUT_DIR,"part6_temporal_lead_lag",
                          "chromatin_first_arrows"),
                p_arrow, w=8, h=7)
    }
  }
}

cat(sprintf("\n[v3] Done. Outputs in %s/\n", OUT_DIR))
cat("[v3] Key outputs:\n")
cat("[v3]   part5b_bistable_expanded/ — bistable*.pdf, bistable_annotated_full.tsv\n")
cat("[v3]   part6_temporal_lead_lag/  — lead_lag_summary.pdf, trajectory_heatmap.pdf\n")
cat("[v3]                               lead_lag_full_table.tsv (all genes, all TPs)\n")
