#!/usr/bin/env Rscript
# =============================================================================
# atac_rna_integration_v1.R  v2
#
# Three-part ATAC × RNA integration for BOTv vs BOTCv.
# Run from the ATAC-seq working directory (same level as ATAC_DARS/ and
# Overview_Plots/). Point DEG_DIR and CLUSTER_DIR at your CEL-Seq2 outputs.
#
# USAGE
#   Rscript atac_rna_integration_v1.R
#   Rscript atac_rna_integration_v1.R --parts 1,2,3   (include Part 2 GSEA-style enrichment)
#   Rscript atac_rna_integration_v1.R --n-perm 1000
#   Rscript atac_rna_integration_v1.R --fdr-rna 0.10
# =============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(ggplot2); library(dplyr)
  library(patchwork);     library(scales)
})
has_ggrepel <- requireNamespace("ggrepel", quietly=TRUE)
has_tidyr   <- requireNamespace("tidyr",   quietly=TRUE)  # pivot_wider, used by 3f/3h

# =============================================================================
# CONFIG
# =============================================================================

DAR_DIR     <- "./ATAC_DARS"
FATE_DIR    <- "./Overview_Plots/peak_fate/BOTv_vs_BOTCv"
DEG_DIR     <- "../../seRNA-seq_03_19_26/celseq2_pipeline/results/combined_figures/qc_and_deg_botv_botcv/deg_limma"
CLUSTER_DIR <- "../../seRNA-seq_03_19_26/celseq2_pipeline/results/combined_figures/explore_clusters_botv_botcv"
OUT_DIR     <- "./Overview_Plots/atac_rna_integration"

FDR_RNA    <- 0.10
LFC_ATAC   <- 0.50
TSS_WINDOW <- 5000
N_PERM     <- 10000
TOP_LABEL  <- 30

# [v16 integration] 01_peak_fate_classification.r now writes a `nearby_genes` column
# alongside its single-nearest-gene SYMBOL call: every gene within a fixed
# annotation window, formatted "SYM(dist_bp);SYM(dist_bp);...", closest
# first (or the literal "none within Xkb"). GENE_WINDOW_BP must match that
# script's gene_window_bp so "all" below reflects its actual search radius,
# not an arbitrary cutoff. NEARBY_CONFIDENT_BP is this script's own choice
# of a stricter secondary radius for the "confident expansion" variant used
# in Part 3's robustness check (3m) -- tune it down if 10kb feels too
# permissive for your locus density.
GENE_WINDOW_BP      <- 10000
NEARBY_CONFIDENT_BP <- 5000

COL_BOTV  <- "#2ca02c"; COL_BOTCV <- "#e377c2"
COL_OPEN  <- "#e31a1c"; COL_CLOSE <- "#1f78b4"
COL_DISC  <- "#aaaaaa"; COL_SIG   <- "#ff7f00"

# Temporal colour palette — light=early, dark=late within each genotype
COL_BOTCV_EARLY <- "#e377c2"   # BOTCv nc14b  (light pink)
COL_BOTCV_LATE  <- "#9e1f8e"   # BOTCv late   (dark magenta)
COL_BOTV_EARLY  <- "#2ca02c"   # BOTv  nc14b  (light green)
COL_BOTV_LATE   <- "#1a5c1a"   # BOTv  late   (dark forest green)
COL_CROSS       <- "#7b4fb7"   # cross-genotype (purple)

# Named colour map for fate gene sets → display colour
# BOTCv fates → dark magenta; BOTv fates → dark green; Emerging → orange
#
# Determines which genotype a fate_short token is ULTIMATELY biased toward
# (i.e. which colour/label it should carry) -- NOT which genotype prefixes
# its raw key. For every category except Reversed these are the same thing
# ("BOTv_Deepened" stays BOTv-biased throughout, start to end). Reversed is
# the one place they diverge: plot_peak_fate_v11.r names "BOTv_Reversed" for
# a peak that STARTED BOTv-open and reversed INTO BOTCv-open, so its final
# state -- and therefore its correct colour -- is BOTC.v, not BOT.v.
# pi_fate_label() already gets this right for the text label (see its
# comment below); this is the same ending-genotype logic for colour, kept
# as ONE shared function so label and colour can never silently drift apart
# again the way fate_col_map/geno_group/fate_display_col previously did
# (each had its own independent, naive grepl("BOTCv"/"BOTv", fate_short)
# that only checked the raw key's prefix, so every Reversed fate rendered
# in its STARTING genotype's colour instead of its ending one).
fate_genotype_bias <- function(fate_short_vec) {
  vapply(fate_short_vec, function(x) {
    x <- sub("^genes_", "", x)
    if (grepl("^Emerg", x))         return("Emerged")
    if (grepl("^BOTv_Reversed",  x)) return("BOTCv")  # ended BOTCv-open
    if (grepl("^BOTCv_Reversed", x)) return("BOTv")   # ended BOTv-open
    if (grepl("^BOTv",  x))          return("BOTv")
    if (grepl("^BOTCv", x))          return("BOTCv")
    NA_character_
  }, character(1), USE.NAMES=FALSE)
}

fate_display_col <- function(fate_short_vec) {
  bias <- fate_genotype_bias(fate_short_vec)
  dplyr::case_when(
    bias == "BOTCv"   ~ COL_BOTCV_LATE,
    bias == "BOTv"    ~ COL_BOTV_LATE,
    bias == "Emerged" ~ "#F39C12",
    TRUE ~ "#888888")
}

FATE_COLS <- c(Deepened="#c9579a", Maintained="#e8a8cf",
               Reversed="#8B0057", Converged="#AAAAAA", Emerged="#F39C12")

# ── PI-requested display label translation ──────────────────────────────────
# fate_short (e.g. "BOTv_Deepened", "BOTCv_Converged_BOTv_opens") is the
# STABLE INTERNAL TOKEN derived from plot_peak_fate_v11.r's gene-list
# filenames, and is what all the grepl()-based matching/coloring in this
# script keys off of (fate_display_col above, genotype extraction, Converged
# filtering, etc.). Do NOT change fate_short itself -- that's exactly what
# broke once already when plot_peak_fate_v11.r's PI-facing display label
# changed format (see that script's A13 section for the fix that keeps
# fate_short stable regardless of display wording).
#
# This function ONLY translates fate_short to the PI-requested DISPLAY text,
# for use in scale_*(labels=pi_fate_label) or wrapped around fate_short when
# building a title/label string directly. Apply it at render time; never by
# overwriting the fate_short column, or every grepl() check downstream of
# that point silently breaks the same way the original bug did.
pi_fate_label <- function(fs) {
  vapply(fs, function(x) {
    x <- sub("^genes_", "", x)   # tolerate the raw gene-list-file key too
    # Exact matches for the genotype-collapsed significance-test categories
    # (all 6 Converged mechanisms unioned per genotype) -- must be checked
    # before the mechanism-level "^BOTv_Converged_" pattern below, since
    # these have no trailing "_mechanism" suffix to match against.
    if (x == "BOTv_Converged")  return("Converged (BOT.v mechanisms, combined)")
    if (x == "BOTCv_Converged") return("Converged (BOTC.v mechanisms, combined)")
    if (grepl("^BOTv_Converged_", x))
      return(paste0("Converged: ", gsub("_", " ", sub("^BOTv_Converged_", "", x))))
    if (grepl("^BOTCv_Converged_", x))
      return(paste0("Converged: ", gsub("_", " ", sub("^BOTCv_Converged_", "", x))))
    if (grepl("^Emerged", x)) return("Emerged")
    if (x == "BOTv_Deepened")    return("Deepened BOT.v bias")
    if (x == "BOTCv_Deepened")   return("Deepened BOTC.v bias")
    if (x == "BOTv_Maintained")  return("Maintained BOT.v bias")
    if (x == "BOTCv_Maintained") return("Maintained BOTC.v bias")
    # IMPORTANT: Reversed is keyed by ENDING genotype, matching
    # plot_peak_fate_v11.r's lbl_rev_A/lbl_rev_B convention exactly -- a
    # peak that STARTED BOTv-open and reversed INTO BOTCv-open is "Reversed
    # to BOTC.v", not "Reversed to BOT.v". Do not swap these.
    if (x == "BOTv_Reversed")    return("Reversed to BOTC.v")
    if (x == "BOTCv_Reversed")   return("Reversed to BOT.v")
    x   # unrecognized token -- show as-is rather than silently hide it
  }, character(1), USE.NAMES=FALSE)
}

ATAC_STEMS <- list(
  list(stem="BOTv_vs_BOTCv",          label="BOTCv_vs_BOTv_nc14b",
       negate=TRUE,  timepoint="nc14b"),
  list(stem="BOTv_vs_BOTCv_nc14late", label="BOTCv_vs_BOTv_nc14late",
       negate=TRUE,  timepoint="nc14late"),
  list(stem="BOTCv_temporal",          label="BOTCv_temporal",
       negate=FALSE, timepoint="BOTCv"),
  list(stem="BOTv_temporal",           label="BOTv_temporal",
       negate=FALSE, timepoint="BOTv")
)

RNA_CONTRASTS <- list(
  list(name="BOTCv_nc14b_vs_BOTv_nc14b",  label="BOTCv vs BOTv (nc14b)",
       atac_pairs=c("BOTCv_vs_BOTv_nc14b"), positive_is="BOTCv"),
  list(name="BOTCv_late_vs_BOTv_nc14d",   label="BOTCv late vs BOTv nc14d",
       atac_pairs=c("BOTCv_vs_BOTv_nc14late"), positive_is="BOTCv"),
  list(name="BOTCv_late_vs_BOTCv_nc14b",  label="BOTCv temporal",
       atac_pairs=c("BOTCv_temporal"), positive_is="BOTCv_late"),
  list(name="BOTv_nc14d_vs_BOTv_nc14b",   label="BOTv temporal",
       atac_pairs=c("BOTv_temporal"), positive_is="BOTv_late"),
  list(name="BOTCv_late_vs_BOTv_nc14b",   label="BOTCv late vs BOTv early",
       atac_pairs=c("BOTCv_vs_BOTv_nc14b","BOTCv_vs_BOTv_nc14late"),
       positive_is="BOTCv")
)

# =============================================================================
# CLI
# =============================================================================
args <- commandArgs(trailingOnly=TRUE)
parse_arg <- function(flag, default) {
  i <- which(args == flag)
  if (length(i) && length(args) > i) return(args[i+1])
  default
}
run_parts  <- as.integer(strsplit(parse_arg("--parts","1,3"),",")[[1]])
# [default changed] Part 2 (GSEA-style rank enrichment) is N_PERM-heavy and
# unused so far -- dropped from the default part list to skip its runtime
# cost automatically. It's fully self-contained (re-reads its own fate gene
# lists, nothing in Part 3 depends on anything it computes), so this is
# purely a time saving, not a behavior change to Parts 1/3. Add it back for
# a specific run with `--parts 1,2,3`.
FDR_RNA    <- as.numeric(parse_arg("--fdr-rna",  FDR_RNA))
LFC_ATAC   <- as.numeric(parse_arg("--lfc-atac", LFC_ATAC))
N_PERM     <- as.integer(parse_arg("--n-perm",   N_PERM))

cat(sprintf("[integrate] Parts: %s  |  RNA FDR=%.2f  |  N_PERM=%d\n",
            paste(run_parts,collapse=","), FDR_RNA, N_PERM))

for (d in c("part1_concordance","part2_enrichment","part3_cluster_fate"))
  dir.create(file.path(OUT_DIR,d), recursive=TRUE, showWarnings=FALSE)

# =============================================================================
# HELPERS
# =============================================================================

save_plot <- function(stem, p, w=9, h=6, dpi=200) {
  ggsave(paste0(stem,".pdf"), p, width=w, height=h, device="pdf")
  ggsave(paste0(stem,".png"), p, width=w, height=h, dpi=dpi)
  cat(sprintf("[integrate]   Saved %s\n", basename(stem)))
}

read_rna <- function(contrast_name) {
  path <- file.path(DEG_DIR, contrast_name, "results_full.tsv")
  if (!file.exists(path)) { message("[integrate] Not found: ", path); return(NULL) }
  df <- tryCatch(read.table(path, header=TRUE, sep="\t", quote="",
                             stringsAsFactors=FALSE),
                 error=function(e) NULL)
  if (is.null(df)) return(NULL)
  df <- df[!is.na(df$log2FoldChange) & !is.na(df$t), ]
  cat(sprintf("[integrate]   RNA %s: %d genes\n", contrast_name, nrow(df)))
  df
}

read_atac_annotated <- function(stem, negate=FALSE) {
  for (ext in c("_DARs_annotated.txt","_DARs_annotated_peaks.csv")) {
    path <- file.path(DAR_DIR, paste0(stem, ext))
    if (!file.exists(path)) next
    sep_ch <- if (grepl("\t", readLines(path,1,warn=FALSE), fixed=TRUE)) "\t" else ","
    df <- tryCatch(read.table(path, header=TRUE, sep=sep_ch,
                              stringsAsFactors=FALSE, quote="",
                              fill=TRUE, comment.char="", row.names=NULL),
                   error=function(e) NULL)
    if (is.null(df) || nrow(df)==0) next
    colnames(df)[1] <- gsub("^[^[:alnum:]]+","",colnames(df)[1])

    cat(sprintf("[integrate]   ATAC %s cols: %s\n", stem,
                paste(head(colnames(df),12), collapse=", ")))

    lfc_col  <- intersect(c("log2FoldChange","log2FC","LFC"),        colnames(df))[1]
    sym_col  <- intersect(c("SYMBOL","Symbol","symbol","geneSymbol"),colnames(df))[1]
    gene_col <- intersect(c("geneId","GeneID","gene_id","GENEID"),   colnames(df))[1]
    ann_col  <- intersect(c("annotation","Annotation"),              colnames(df))[1]
    chr_col  <- intersect(c("seqnames","chr","Chr","chrom"),         colnames(df))[1]
    start_col<- intersect(c("start","Start"),                        colnames(df))[1]
    end_col  <- intersect(c("end","End"),                            colnames(df))[1]

    if (is.na(lfc_col)) { message("[integrate]   No LFC col in ",path); next }

    out <- data.frame(
      peak_id    = if ("name" %in% colnames(df)) df$name
                   else paste0(stem,"_",seq_len(nrow(df))),
      chr        = if (!is.na(chr_col))   as.character(df[[chr_col]])   else NA_character_,
      start      = if (!is.na(start_col)) as.integer(df[[start_col]])   else NA_integer_,
      end        = if (!is.na(end_col))   as.integer(df[[end_col]])     else NA_integer_,
      atac_lfc   = as.numeric(df[[lfc_col]]) * ifelse(negate,-1,1),
      symbol     = if (!is.na(sym_col))   as.character(df[[sym_col]])   else NA_character_,
      gene_id    = if (!is.na(gene_col))  as.character(df[[gene_col]])  else NA_character_,
      annotation = if (!is.na(ann_col))   as.character(df[[ann_col]])   else NA_character_,
      stringsAsFactors=FALSE)

    # Completeness check: this function trusts the annotated file's own row
    # count with no reference to the current bed file -- if the annotated
    # file has fallen out of sync with the bed file (e.g. regenerated with
    # more peaks after the annotated file was last written), this silently
    # returns fewer peaks than actually exist, with no indication anything
    # is missing. Compare row counts against the bed file explicitly so
    # that gap is visible instead of invisible.
    bed_check_path <- file.path(DAR_DIR, paste0(stem, "_DARs.bed"))
    if (file.exists(bed_check_path)) {
      n_bed <- tryCatch(length(readLines(bed_check_path)), error=function(e) NA_integer_)
      if (!is.na(n_bed) && n_bed != nrow(out)) {
        cat(sprintf("[integrate]   [WARNING] %s: annotated file has %d peaks but current bed file has %d -- this annotated file is likely OUT OF DATE. Regenerate it before trusting these results.\n",
                   stem, nrow(out), n_bed))
      }
    }

    # If no chr/start/end from this file, try the BED file
    if (all(is.na(out$chr))) {
      bed_path <- file.path(DAR_DIR, paste0(stem, "_DARs.bed"))
      if (file.exists(bed_path)) {
        bed <- tryCatch(read.table(bed_path, sep="\t", header=FALSE,
                                   stringsAsFactors=FALSE), error=function(e) NULL)
        if (!is.null(bed) && nrow(bed) == nrow(out)) {
          out$chr   <- as.character(bed[[1]])
          out$start <- as.integer(bed[[2]])
          out$end   <- as.integer(bed[[3]])
          cat(sprintf("[integrate]     Coordinates loaded from BED file\n"))
        }
      }
    }

    n_sym <- sum(!is.na(out$symbol) & out$symbol != "")
    n_gid <- sum(!is.na(out$gene_id) & out$gene_id != "")
    cat(sprintf("[integrate]   ATAC %s: %d peaks  |  %d with symbol  |  %d with gene_id\n",
                stem, nrow(out), n_sym, n_gid))
    return(out)
  }
  # No annotated file — try just the BED
  bed_path <- file.path(DAR_DIR, paste0(stem, "_DARs.bed"))
  if (file.exists(bed_path)) {
    bed <- tryCatch(read.table(bed_path, sep="\t", header=FALSE,
                               stringsAsFactors=FALSE), error=function(e) NULL)
    if (!is.null(bed) && nrow(bed)>0) {
      lfc_raw <- if (ncol(bed)>=5) as.numeric(bed[[5]]) else rep(NA_real_, nrow(bed))
      out <- data.frame(
        peak_id=paste0(stem,"_",seq_len(nrow(bed))),
        chr=as.character(bed[[1]]), start=as.integer(bed[[2]]),
        end=as.integer(bed[[3]]),
        atac_lfc=lfc_raw * ifelse(negate,-1,1),
        symbol=NA_character_, gene_id=NA_character_, annotation=NA_character_,
        stringsAsFactors=FALSE)
      cat(sprintf("[integrate]   ATAC %s: %d peaks from BED (no annotation cols)\n",
                  stem, nrow(out)))
      return(out)
    }
  }
  message(sprintf("[integrate] WARNING: No ATAC file found for '%s'", stem))
  NULL
}

# Annotate peaks to genes via genomic overlap with TSS window
# Uses TxDb if available, falls back to returning unannotated df
annotate_peaks_to_genes <- function(atac_df, tss_window=TSS_WINDOW) {
  if (!requireNamespace("TxDb.Dmelanogaster.UCSC.dm6.ensGene", quietly=TRUE)) {
    cat("[integrate]   TxDb not available — skipping genomic annotation\n")
    return(atac_df)
  }
  if (all(is.na(atac_df$chr))) return(atac_df)

  suppressPackageStartupMessages({
    library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
  })
  txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene

  # Get TSS positions
  genes_gr <- suppressWarnings(genes(txdb))
  tss_gr   <- resize(genes_gr, width=1, fix="start")
  tss_gr   <- suppressWarnings(trim(resize(tss_gr,
                                            width=tss_window*2+1, fix="center")))

  # Build peak GRanges (use 1-based)
  valid <- !is.na(atac_df$chr) & !is.na(atac_df$start) & !is.na(atac_df$end)
  if (sum(valid)==0) return(atac_df)
  peak_gr <- GRanges(
    seqnames=atac_df$chr[valid],
    ranges=IRanges(start=atac_df$start[valid]+1L, end=atac_df$end[valid]))
  mcols(peak_gr)$peak_idx <- which(valid)

  # Overlap
  hits <- suppressWarnings(findOverlaps(peak_gr, tss_gr))
  if (length(hits)==0) {
    cat("[integrate]   No peak-gene overlaps found within TSS window\n")
    return(atac_df)
  }

  # Get gene_ids from hits
  hit_df <- data.frame(
    peak_idx = mcols(peak_gr)$peak_idx[queryHits(hits)],
    gene_id  = names(tss_gr)[subjectHits(hits)],
    stringsAsFactors=FALSE)

  # Add symbols if org.Dm.eg.db available
  if (requireNamespace("org.Dm.eg.db", quietly=TRUE)) {
    suppressPackageStartupMessages(library(org.Dm.eg.db))
    # TxDb.Dmelanogaster.UCSC.dm6.ensGene uses Entrez IDs as gene names.
    # Try ENTREZID keytype first; if that returns no symbols try FLYBASE.
    sym_map <- tryCatch(
      AnnotationDbi::select(org.Dm.eg.db,
                            keys=unique(hit_df$gene_id),
                            columns="SYMBOL", keytype="ENTREZID"),
      error=function(e) NULL)
    # Check if we actually got symbols
    got_syms <- !is.null(sym_map) && "SYMBOL" %in% colnames(sym_map) &&
                any(!is.na(sym_map$SYMBOL))
    if (!got_syms) {
      # Fallback: try FLYBASE keytype
      sym_map <- tryCatch(
        AnnotationDbi::select(org.Dm.eg.db,
                              keys=unique(hit_df$gene_id),
                              columns="SYMBOL", keytype="FLYBASE"),
        error=function(e) NULL)
      got_syms <- !is.null(sym_map) && "SYMBOL" %in% colnames(sym_map) &&
                  any(!is.na(sym_map$SYMBOL))
    }
    if (got_syms) {
      key_col <- intersect(c("ENTREZID","FLYBASE"),colnames(sym_map))[1]
      hit_df  <- merge(hit_df, sym_map[,c(key_col,"SYMBOL"),drop=FALSE],
                       by.x="gene_id", by.y=key_col, all.x=TRUE)
      cat(sprintf("[integrate]     Symbol sample (org.Db): %s\n",
                  paste(head(na.omit(hit_df$SYMBOL),5), collapse=", ")))
    } else {
      hit_df$SYMBOL <- hit_df$gene_id
      cat("[integrate]     org.Db returned no symbols — using gene IDs as labels\n")
    }
  } else {
    hit_df$SYMBOL <- hit_df$gene_id
  }

  # Take one gene per peak (closest TSS approximation: keep first hit)
  hit_df <- hit_df[!duplicated(hit_df$peak_idx), ]

  atac_df$symbol[hit_df$peak_idx]  <- hit_df$SYMBOL
  atac_df$gene_id[hit_df$peak_idx] <- hit_df$gene_id
  n_ann <- sum(!is.na(atac_df$symbol))
  cat(sprintf("[integrate]   Genomic overlap: %d / %d peaks annotated to genes\n",
              n_ann, nrow(atac_df)))
  atac_df
}

# Merge ATAC (with gene annotation) against RNA results
# Tries: symbol → gene_id → genomic annotation on the fly
merge_atac_rna <- function(atac_df, rna_df) {
  if (is.null(atac_df) || is.null(rna_df)) return(NULL)

  # If ATAC has no gene annotations, try to add them via genomic overlap
  if (all(is.na(atac_df$symbol)) && all(is.na(atac_df$gene_id))) {
    cat("[integrate]   No gene annotations in ATAC file — trying genomic overlap\n")
    atac_df <- annotate_peaks_to_genes(atac_df)
  }

  rna_sub <- rna_df[, intersect(c("symbol","gene_id","log2FoldChange","t","pvalue","padj"),
                                 colnames(rna_df)), drop=FALSE]
  # Standardise RNA col names
  names(rna_sub)[names(rna_sub)=="symbol"]         <- "rna_symbol"
  names(rna_sub)[names(rna_sub)=="gene_id"]        <- "rna_gene_id"
  names(rna_sub)[names(rna_sub)=="log2FoldChange"] <- "rna_lfc"
  names(rna_sub)[names(rna_sub)=="t"]              <- "rna_t"
  names(rna_sub)[names(rna_sub)=="pvalue"]         <- "rna_pval"
  names(rna_sub)[names(rna_sub)=="padj"]           <- "rna_padj"

  merged <- NULL

  # Try symbol merge
  n_atac_sym <- sum(!is.na(atac_df$symbol) & atac_df$symbol != "")
  n_rna_sym  <- if ("rna_symbol" %in% names(rna_sub))
    sum(!is.na(rna_sub$rna_symbol) & rna_sub$rna_symbol != "") else 0

  if (n_atac_sym > 0 && n_rna_sym > 0) {
    m1 <- merge(atac_df[!is.na(atac_df$symbol) & atac_df$symbol!="", ],
                rna_sub[!is.na(rna_sub$rna_symbol) & rna_sub$rna_symbol!="", ],
                by.x="symbol", by.y="rna_symbol", all=FALSE)
    cat(sprintf("[integrate]   Merge by symbol: %d rows\n", nrow(m1)))
    if (nrow(m1)>0) merged <- m1
  }

  # Try gene_id merge (Entrez or FBgn)
  if ((is.null(merged)||nrow(merged)==0) &&
      "rna_gene_id" %in% names(rna_sub)) {
    n_atac_gid <- sum(!is.na(atac_df$gene_id) & atac_df$gene_id!="")
    n_rna_gid  <- sum(!is.na(rna_sub$rna_gene_id) & rna_sub$rna_gene_id!="")
    if (n_atac_gid>0 && n_rna_gid>0) {
      m2 <- merge(atac_df[!is.na(atac_df$gene_id) & atac_df$gene_id!="", ],
                  rna_sub[!is.na(rna_sub$rna_gene_id) & rna_sub$rna_gene_id!="", ],
                  by.x="gene_id", by.y="rna_gene_id", all=FALSE)
      cat(sprintf("[integrate]   Merge by gene_id: %d rows\n", nrow(m2)))
      if (nrow(m2)>0) {
        if (!"symbol" %in% colnames(m2)) m2$symbol <- m2$gene_id
        merged <- m2
      }
    }
  }

  if (is.null(merged)||nrow(merged)==0) {
    cat("[integrate]   WARNING: 0 genes merged after all attempts.\n")
    cat("[integrate]   Install TxDb.Dmelanogaster.UCSC.dm6.ensGene + org.Dm.eg.db\n")
    cat("[integrate]   for automatic genomic overlap annotation:\n")
    cat("[integrate]   BiocManager::install(c('TxDb.Dmelanogaster.UCSC.dm6.ensGene',\n")
    cat("[integrate]                           'org.Dm.eg.db'))\n")
    # Show what we have for diagnosis
    cat(sprintf("[integrate]   ATAC: %d non-NA symbols, %d non-NA gene_ids\n",
                sum(!is.na(atac_df$symbol)), sum(!is.na(atac_df$gene_id))))
    if ("rna_symbol" %in% names(rna_sub))
      cat(sprintf("[integrate]   RNA:  first 5 symbols: %s\n",
                  paste(head(rna_sub$rna_symbol[!is.na(rna_sub$rna_symbol)],5),
                        collapse=", ")))
    return(NULL)
  }

  merged <- merged[!is.na(merged$atac_lfc) & !is.na(merged$rna_lfc), ]
  merged$quadrant <- with(merged, dplyr::case_when(
    atac_lfc>0 & rna_lfc>0 ~ "Concordant\nOpen+Up",
    atac_lfc<0 & rna_lfc<0 ~ "Concordant\nClose+Down",
    atac_lfc>0 & rna_lfc<0 ~ "Discordant\nOpen+Down",
    TRUE                    ~ "Discordant\nClose+Up"))
  merged$sig_rna  <- !is.na(merged$rna_padj) & merged$rna_padj < FDR_RNA
  merged$sig_atac <- abs(merged$atac_lfc) > LFC_ATAC

  tab <- table(merged$quadrant)
  cat(sprintf("[integrate]   Merged: %d genes\n", nrow(merged)))
  for (q in names(tab)) cat(sprintf("[integrate]     %s: %d\n",
                                     gsub("\n"," ",q), tab[[q]]))
  merged
}

read_fate_gene_lists <- function() {
  gene_dir <- file.path(FATE_DIR, "gene_lists")
  if (!dir.exists(gene_dir)) {
    message("[integrate] gene_lists/ not found — run plot_peak_fate_v11.r first")
    return(list())
  }
  files <- list.files(gene_dir, pattern="\\.csv$", full.names=TRUE)
  out   <- lapply(files, function(f) {
    df <- tryCatch(read.csv(f, stringsAsFactors=FALSE), error=function(e) NULL)
    if (is.null(df)) return(NULL)
    sc <- intersect(c("SYMBOL","symbol","Symbol"), colnames(df))[1]
    if (is.na(sc)) return(NULL)
    unique(df[[sc]][!is.na(df[[sc]]) & df[[sc]] != ""])
  })
  names(out) <- gsub("\\.csv$","", basename(files))
  out <- Filter(function(x) !is.null(x) && length(x) > 0, out)
  cat(sprintf("[integrate]   Fate gene lists: %s\n", paste(names(out), collapse=", ")))
  out
}

# Parses one 01_peak_fate_classification.r `nearby_genes` cell: "SYM(dist_bp);..."
# closest first, or the literal sentinel "none within Xkb" (or NA/"" for
# older exports). Returns data.frame(symbol, dist_bp) with 0 rows if there's
# nothing to parse -- never NULL, so callers can rbind without a NULL check.
parse_nearby_genes_cell <- function(x) {
  empty <- data.frame(symbol=character(0), dist_bp=numeric(0), stringsAsFactors=FALSE)
  if (is.na(x) || x == "" || grepl("^none within", x, ignore.case=TRUE)) return(empty)
  parts <- strsplit(x, ";", fixed=TRUE)[[1]]
  m <- regmatches(parts, regexec("^(.*)\\(([0-9]+(?:\\.[0-9]+)?)\\)$", parts))
  m <- Filter(function(mm) length(mm) == 3, m)
  if (length(m) == 0) return(empty)
  data.frame(symbol  = vapply(m, `[`, character(1), 2),
             dist_bp = as.numeric(vapply(m, `[`, character(1), 3)),
             stringsAsFactors=FALSE)
}

# Companion reader to read_fate_gene_lists() that ALSO uses the [v16]
# nearby_genes column to build two broader variants of each fate's gene
# set: "confident" (nearest-gene calls + any gene within NEARBY_CONFIDENT_BP)
# and "all" (nearest-gene calls + every gene within the full GENE_WINDOW_BP
# annotation window). Peak-to-gene assignment in ATAC data isn't always the
# single nearest gene -- ChIPseeker only ever reports one -- so Part 3 runs
# its enrichment tests on all three variants and reports which cluster-fate
# calls are robust to that choice, rather than silently picking one.
# Deliberately a SEPARATE function from read_fate_gene_lists() (not a
# replacement / not an added return field): every other part of this script,
# and 01_peak_fate_classification.r itself, treats the nearest-gene SYMBOL call as
# the primary annotation, and nothing about that should change silently.
read_fate_gene_lists_nearby_expanded <- function() {
  gene_dir <- file.path(FATE_DIR, "gene_lists")
  if (!dir.exists(gene_dir)) return(list(confident=list(), all=list()))
  files <- list.files(gene_dir, pattern="\\.csv$", full.names=TRUE)
  confident <- list(); all_exp <- list()
  had_nearby_col <- FALSE
  for (f in files) {
    df <- tryCatch(read.csv(f, stringsAsFactors=FALSE), error=function(e) NULL)
    if (is.null(df)) next
    sc <- intersect(c("SYMBOL","symbol","Symbol"), colnames(df))[1]
    if (is.na(sc)) next
    nearest <- unique(df[[sc]][!is.na(df[[sc]]) & df[[sc]] != ""])
    key <- gsub("\\.csv$","", basename(f))
    if (!("nearby_genes" %in% colnames(df))) {
      # No nearby_genes column in this export (older plot_peak_fate run) --
      # expanded == nearest, so this fate contributes no NEW evidence, but
      # it's still present (not silently dropped from the comparison).
      confident[[key]] <- nearest
      all_exp[[key]]   <- nearest
      next
    }
    had_nearby_col <- TRUE
    parsed <- do.call(rbind, lapply(df$nearby_genes, parse_nearby_genes_cell))
    conf_syms <- if (!is.null(parsed) && nrow(parsed) > 0)
      unique(parsed$symbol[parsed$dist_bp <= NEARBY_CONFIDENT_BP]) else character(0)
    all_syms_nb <- if (!is.null(parsed) && nrow(parsed) > 0)
      unique(parsed$symbol) else character(0)
    confident[[key]] <- unique(c(nearest, conf_syms))
    all_exp[[key]]   <- unique(c(nearest, all_syms_nb))
  }
  if (!had_nearby_col)
    message("[integrate]   nearby_genes column not found in any gene_lists CSV -- ",
            "3m robustness check will show no difference from nearest-only ",
            "(re-run 01_peak_fate_classification.r to populate it)")
  list(confident = Filter(function(x) length(x) > 0, confident),
       all       = Filter(function(x) length(x) > 0, all_exp))
}

# GSEA-style running enrichment (weighted KS)
running_enrichment <- function(gene_set, ranked_genes, t_scores) {
  n    <- length(ranked_genes)
  hits <- ranked_genes %in% gene_set
  nh   <- sum(hits)
  if (nh == 0 || nh == n)
    return(list(ES=0, running=rep(0,n), positions=which(hits)))
  wt   <- ifelse(hits, abs(t_scores), 0)
  wt_sum <- sum(wt)
  miss_p <- 1 / (n - nh)
  running <- cumsum(ifelse(hits, wt/wt_sum, -miss_p))
  ES_idx  <- which.max(abs(running))
  list(ES=running[ES_idx], running=running, positions=which(hits))
}

# NES: divide by mean of |permuted ES| (more stable than sd for small perm n)
compute_NES <- function(ES_obs, ES_perm) {
  pos_perm <- ES_perm[ES_perm >= 0]
  neg_perm <- ES_perm[ES_perm <  0]
  if (ES_obs >= 0) {
    mu <- if (length(pos_perm) > 0) mean(pos_perm) else mean(abs(ES_perm))
  } else {
    mu <- if (length(neg_perm) > 0) mean(abs(neg_perm)) else mean(abs(ES_perm))
  }
  if (mu < 1e-9) return(0)
  ES_obs / mu
}

# Simplify annotation categories
simplify_ann <- function(x) {
  dplyr::case_when(
    grepl("Promoter",     x, ignore.case=TRUE) ~ "Promoter",
    grepl("5.*UTR",       x, ignore.case=TRUE) ~ "5' UTR",
    grepl("3.*UTR",       x, ignore.case=TRUE) ~ "3' UTR",
    grepl("1st Exon",     x, ignore.case=TRUE) ~ "1st Exon",
    grepl("Exon",         x, ignore.case=TRUE) ~ "Other Exon",
    grepl("1st Intron",   x, ignore.case=TRUE) ~ "1st Intron",
    grepl("Intron",       x, ignore.case=TRUE) ~ "Other Intron",
    grepl("Downstream",   x, ignore.case=TRUE) ~ "Downstream",
    TRUE ~ "Distal Intergenic")
}

# Significance stars
sig_stars <- function(p) {
  case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*",
            p < 0.10  ~ ".",   TRUE ~ "")
}

# =============================================================================
# PART 1 — CONCORDANCE BINS
# =============================================================================

if (1 %in% run_parts) {
  cat("\n", strrep("=",70), "\n", sep="")
  cat("PART 1: CONCORDANCE BINS\n\n")

  atac_cache <- list()

  for (rna_cfg in RNA_CONTRASTS) {
    rna_df <- read_rna(rna_cfg$name)
    if (is.null(rna_df)) next

    for (atac_label in rna_cfg$atac_pairs) {
      atac_cfg <- Filter(function(x) x$label == atac_label, ATAC_STEMS)
      if (length(atac_cfg)==0) next
      atac_cfg <- atac_cfg[[1]]

      if (!atac_label %in% names(atac_cache))
        atac_cache[[atac_label]] <- read_atac_annotated(atac_cfg$stem, atac_cfg$negate)
      atac_df <- atac_cache[[atac_label]]
      if (is.null(atac_df)) next

      merged <- merge_atac_rna(atac_df, rna_df)
      if (is.null(merged) || nrow(merged) < 5) {
        cat(sprintf("[integrate]   Skipping plot for %s × %s (too few genes)\n",
                    rna_cfg$name, atac_label))
        next
      }

      # ── 1a: Concordance scatter ──────────────────────────────────────────
      quad_cols <- c("Concordant\nOpen+Up"    = COL_OPEN,
                     "Concordant\nClose+Down" = COL_CLOSE,
                     "Discordant\nOpen+Down"  = "#cccccc",
                     "Discordant\nClose+Up"   = "#dddddd")

      # Label top genes: most extreme on each axis among sig
      to_label <- merged %>%
        filter(sig_rna | sig_atac) %>%
        mutate(score = abs(rna_lfc) * abs(atac_lfc)) %>%
        arrange(desc(score)) %>%
        head(TOP_LABEL)

      p_scatter <- ggplot(merged, aes(x=atac_lfc, y=rna_lfc)) +
        geom_hline(yintercept=0, lty=2, colour="#cccccc") +
        geom_vline(xintercept=0, lty=2, colour="#cccccc") +
        geom_vline(xintercept=c(-LFC_ATAC,LFC_ATAC), lty=3, colour="#999", linewidth=0.4) +
        geom_point(aes(colour=quadrant,
                       alpha=ifelse(sig_atac | sig_rna, 0.85, 0.2)),
                   size=0.9) +
        scale_colour_manual(values=quad_cols, name="Quadrant") +
        scale_alpha_identity() +
        geom_point(data=filter(merged, sig_rna),
                   aes(x=atac_lfc, y=rna_lfc),
                   shape=21, colour=COL_SIG, fill=NA, size=2.5,
                   stroke=0.8, inherit.aes=FALSE) +
        {
          if (nrow(to_label)>0 && has_ggrepel) {
            library(ggrepel)
            ggrepel::geom_text_repel(data=to_label,
              aes(x=atac_lfc, y=rna_lfc, label=symbol),
              colour="black", size=2.8, fontface="italic",
              max.overlaps=20, box.padding=0.4,
              segment.colour="#888", segment.size=0.3, seed=42,
              inherit.aes=FALSE)
          }
        } +
        labs(title=sprintf("Accessibility × Expression concordance\n%s  |  ATAC: %s",
                           rna_cfg$label, atac_label),
             subtitle=sprintf("Circled = RNA padj<%.2f  |  Faded = |ATAC lFC|<%.1f  |  n=%d genes",
                              FDR_RNA, LFC_ATAC, nrow(merged)),
             x=sprintf("ATAC log2FC  (%s, positive=BOTCv/later open)", atac_label),
             y=sprintf("RNA log2FC  (%s, positive=%s)", rna_cfg$label, rna_cfg$positive_is)) +
        theme_classic(base_size=11) +
        theme(legend.position="right", plot.title=element_text(size=10.5, face="bold"),
              plot.subtitle=element_text(size=8.5, colour="#555"))

      # ── 1b: Quadrant summary bar ─────────────────────────────────────────
      quad_sum <- merged %>%
        group_by(quadrant) %>%
        summarise(n_total=n(), n_sig=sum(sig_rna), .groups="drop") %>%
        mutate(pct=100*n_total/sum(n_total),
               quad_short=gsub("\n"," ", quadrant))
      p_bar <- ggplot(quad_sum, aes(x=reorder(quad_short,-n_total), y=n_total,
                                     fill=quadrant)) +
        geom_col(width=0.65) +
        geom_text(aes(label=sprintf("%.0f%%\n(%d sig)", pct, n_sig)),
                  vjust=-0.2, size=3.2) +
        scale_fill_manual(values=quad_cols, guide="none") +
        labs(title="Quadrant summary", x=NULL, y="Genes") +
        theme_classic(base_size=10) +
        theme(axis.text.x=element_text(size=8))

      # ── 1c: Annotation breakdown for concordant genes ────────────────────
      conc <- merged %>% filter(grepl("Concordant", quadrant),
                                 !is.na(annotation))
      if (nrow(conc) > 5) {
        conc$ann_simple <- simplify_ann(conc$annotation)
        ann_sum <- conc %>% count(quadrant, ann_simple) %>%
          group_by(quadrant) %>% mutate(pct=100*n/sum(n)) %>% ungroup()
        ann_cols <- c(Promoter="#E74C3C","5' UTR"="#F39C12","1st Exon"="#2ECC71",
                      "Other Exon"="#27AE60","1st Intron"="#3498DB",
                      "Other Intron"="#2980B9",Downstream="#9B59B6",
                      "Distal Intergenic"="#BDC3C7","3' UTR"="#F1C40F")
        p_ann <- ggplot(ann_sum, aes(x=gsub("\n"," ",quadrant), y=pct,
                                      fill=ann_simple)) +
          geom_col(position="stack", width=0.7) +
          geom_text(aes(label=ifelse(pct>=8, paste0(round(pct),"%"), "")),
                    position=position_stack(vjust=0.5), size=2.8, colour="white",
                    fontface="bold") +
          scale_fill_manual(values=ann_cols, name="Feature") +
          coord_flip() +
          labs(title="Genomic features\n(concordant genes)",
               x=NULL, y="% peaks") +
          theme_minimal(base_size=9) +
          theme(panel.grid.major.y=element_blank())
        p_combined <- (p_scatter | (p_bar / p_ann)) +
          plot_layout(widths=c(3,1.4))
      } else {
        p_combined <- p_scatter | p_bar
      }

      stem_out <- file.path(OUT_DIR,"part1_concordance",
                            sprintf("concordance_%s_%s", rna_cfg$name, atac_label))
      save_plot(stem_out, p_combined, w=14, h=6)
      write.table(merged, paste0(stem_out,".tsv"), sep="\t", quote=FALSE, row.names=FALSE)
    }
  }
}

# =============================================================================
# PART 2 — GSEA-STYLE RANK ENRICHMENT
# =============================================================================

if (2 %in% run_parts) {
  cat("\n", strrep("=",70), "\n", sep="")
  cat("PART 2: RANK-BASED ENRICHMENT (GSEA-style)\n\n")

  fate_gene_lists <- read_fate_gene_lists()
  if (length(fate_gene_lists)==0) {
    cat("[integrate]   No fate gene lists — skipping Part 2\n")
  } else {
    enrich_rows <- list()

    for (rna_cfg in RNA_CONTRASTS) {
      rna_df <- read_rna(rna_cfg$name)
      if (is.null(rna_df) || nrow(rna_df)==0) next
      rna_ranked  <- rna_df[order(-rna_df$t), ]
      ranked_syms <- rna_ranked$symbol
      t_scores    <- rna_ranked$t

      for (fate_name in names(fate_gene_lists)) {
        gene_set <- intersect(fate_gene_lists[[fate_name]], ranked_syms)
        n_set    <- length(gene_set)
        if (n_set < 5) {
          cat(sprintf("[integrate]   SKIP %s: only %d genes in ranked list\n",
                      fate_name, n_set))
          next
        }

        # Observed enrichment
        res    <- running_enrichment(gene_set, ranked_syms, t_scores)
        ES_obs <- res$ES

        # Permutation null
        set.seed(42)
        ES_perm <- replicate(N_PERM, {
          perm_set <- sample(ranked_syms, n_set)
          running_enrichment(perm_set, ranked_syms, t_scores)$ES
        })
        NES  <- compute_NES(ES_obs, ES_perm)
        pval <- if (ES_obs >= 0) mean(ES_perm >= ES_obs) else mean(ES_perm <= ES_obs)
        pval <- max(pval, 1/N_PERM)   # floor at 1/N_PERM

        cat(sprintf("[integrate]   %s × %s: ES=%.3f, NES=%.2f, p=%.4f (n=%d)\n",
                    fate_name, rna_cfg$name, ES_obs, NES, pval, n_set))

        enrich_rows[[length(enrich_rows)+1]] <- data.frame(
          fate=fate_name, contrast=rna_cfg$name, contrast_label=rna_cfg$label,
          n_set=n_set, ES=round(ES_obs,4), NES=round(NES,3), pval=round(pval,4),
          stringsAsFactors=FALSE)

        # ── Mountain + barcode plot ─────────────────────────────────────────
        enrich_col <- if (ES_obs >= 0) COL_OPEN else COL_CLOSE
        n_total    <- length(ranked_syms)

        # Permutation null ribbon (5th–95th percentile)
        perm_mat <- replicate(min(100, N_PERM), {
          perm_set <- sample(ranked_syms, n_set)
          running_enrichment(perm_set, ranked_syms, t_scores)$running
        })
        perm_lo  <- apply(perm_mat, 1, quantile, 0.05)
        perm_hi  <- apply(perm_mat, 1, quantile, 0.95)

        plot_df <- data.frame(
          rank    = seq_len(n_total),
          running = res$running,
          perm_lo = perm_lo,
          perm_hi = perm_hi,
          t_score = t_scores,
          is_hit  = ranked_syms %in% gene_set)

        # Top hit genes to label
        hit_df <- data.frame(rank=res$positions,
                             symbol=ranked_syms[res$positions],
                             t_val=t_scores[res$positions]) %>%
          arrange(desc(abs(t_val))) %>% head(8)

        p_mountain <- ggplot(plot_df, aes(x=rank)) +
          geom_ribbon(aes(ymin=perm_lo, ymax=perm_hi), fill="#dddddd", alpha=0.7) +
          geom_hline(yintercept=0, colour="#aaaaaa", linewidth=0.5, lty=2) +
          geom_line(aes(y=running), colour=enrich_col, linewidth=0.9) +
          geom_hline(yintercept=ES_obs, colour=enrich_col,
                     linewidth=0.7, lty=3, alpha=0.7) +
          annotate("text", x=n_total*0.98,
                   y=ES_obs + diff(range(res$running))*0.06,
                   label=sprintf("ES = %.3f\nNES = %.2f\np = %.4f\nn = %d",
                                 ES_obs, NES, pval, n_set),
                   hjust=1, vjust=0, size=3.2, colour="#333333",
                   family="mono") +
          labs(title=sprintf("%s\nvs %s", pi_fate_label(fate_name), rna_cfg$label),
               subtitle="Grey ribbon = 5th–95th percentile of permuted null",
               x=NULL, y="Running enrichment score") +
          theme_classic(base_size=11) +
          theme(plot.title=element_text(size=10, face="bold"),
                axis.text.x=element_blank(), axis.ticks.x=element_blank())

        # Barcode strip with tick marks per hit
        p_barcode <- ggplot(data.frame(pos=res$positions)) +
          geom_segment(aes(x=pos, xend=pos, y=0, yend=1),
                       colour=enrich_col, linewidth=0.4, alpha=0.7) +
          scale_x_continuous(limits=c(1,n_total), expand=c(0,0),
                             breaks=c(1, n_total/2, n_total),
                             labels=c("High\n(up in A)","","Low\n(up in B)")) +
          labs(x="Gene rank", y=NULL) +
          theme_classic(base_size=9) +
          theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
                axis.line.y=element_blank(), panel.background=element_rect(fill="#f5f5f5"),
                plot.margin=margin(0,5,5,5))

        # t-statistic profile line
        p_tstat <- ggplot(plot_df, aes(x=rank, y=t_score)) +
          geom_area(fill="#e0e0e0", alpha=0.6) +
          geom_hline(yintercept=0, colour="#888", linewidth=0.4) +
          scale_x_continuous(limits=c(1,n_total), expand=c(0,0)) +
          labs(x=NULL, y="t-stat") +
          theme_classic(base_size=8) +
          theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

        p_full <- p_mountain / p_tstat / p_barcode +
          plot_layout(heights=c(5,1.2,0.8))

        # ── Per-test density plot: where do fate genes fall in ranking? ──────
        hit_ranks  <- which(ranked_syms %in% gene_set)
        miss_ranks <- which(!(ranked_syms %in% gene_set))
        dens_df    <- data.frame(
          rank  = c(hit_ranks, sample(miss_ranks, min(length(miss_ranks), 2000))),
          group = c(rep("ATAC fate genes",length(hit_ranks)),
                    rep("Background",    min(length(miss_ranks),2000))))

        p_density <- ggplot(dens_df, aes(x=rank, fill=group, colour=group)) +
          geom_density(alpha=0.4, adjust=1.2, linewidth=0.7) +
          scale_fill_manual(values=c("ATAC fate genes"=enrich_col,
                                      "Background"="#aaaaaa"), name=NULL) +
          scale_colour_manual(values=c("ATAC fate genes"=enrich_col,
                                        "Background"="#666666"), name=NULL) +
          geom_vline(xintercept=n_total/2, lty=2, colour="#888", linewidth=0.5) +
          labs(title=sprintf("Rank distribution: %s genes", pi_fate_label(fate_name)),
               subtitle=sprintf("Left = up in %s | Right = up in B | n=%d fate genes",
                                rna_cfg$positive_is, n_set),
               x="Gene rank", y="Density") +
          theme_classic(base_size=11) +
          theme(legend.position="top",
                plot.title=element_text(size=10, face="bold"))

        # Top-hit table panel
        if (nrow(hit_df) > 0) {
          hit_df$rank_pct <- round(100*hit_df$rank/n_total, 1)
          p_tbl <- ggplot(hit_df %>% arrange(rank),
                          aes(x=reorder(symbol,-t_val), y=t_val, fill=t_val>0)) +
            geom_col(width=0.7) +
            geom_text(aes(label=sprintf("rank:%d",rank)),
                      hjust=ifelse(hit_df$t_val>0,-0.1,1.1), size=2.8) +
            scale_fill_manual(values=c("TRUE"=COL_OPEN,"FALSE"=COL_CLOSE),
                              guide="none") +
            coord_flip() +
            labs(title="Top ranked fate genes", x=NULL, y="RNA t-statistic") +
            theme_classic(base_size=10)
          p_right <- p_density / p_tbl
        } else {
          p_right <- p_density
        }

        p_combo2 <- p_full | p_right

        stem_out <- file.path(OUT_DIR,"part2_enrichment",
                              sprintf("enrichment_%s__%s",
                                      gsub("[^A-Za-z0-9]","_",fate_name),
                                      rna_cfg$name))
        save_plot(stem_out, p_full,  w=9, h=7)
        save_plot(paste0(stem_out,"_detail"), p_combo2, w=16, h=7)
      }
    }

    # ── Summary bubble plot (all contrasts × fates) ─────────────────────────
    if (length(enrich_rows) > 0) {
      enrich_df <- do.call(rbind, enrich_rows)
      enrich_df$padj <- p.adjust(enrich_df$pval, method="BH")
      write.table(enrich_df,
                  file.path(OUT_DIR,"part2_enrichment","enrichment_summary.tsv"),
                  sep="\t", quote=FALSE, row.names=FALSE)

      # Clean up fate/contrast labels for plot
      enrich_df$fate_short     <- gsub("genes_","",enrich_df$fate)
      enrich_df$fate_display   <- pi_fate_label(enrich_df$fate_short)
      enrich_df$contrast_short <- gsub("_"," ",enrich_df$contrast_label)
      enrich_df$sig_label      <- sig_stars(enrich_df$padj)
      enrich_df$NES_capped     <- pmax(pmin(enrich_df$NES, 4), -4)

      p_bubble <- ggplot(enrich_df,
                         aes(x=contrast_short, y=fate_display,
                             fill=NES_capped,
                             size=-log10(pmax(enrich_df$pval,1/N_PERM)))) +
        geom_point(shape=21, colour="white") +
        geom_text(aes(label=sig_label), colour="black", size=4, vjust=0.4) +
        scale_fill_gradient2(low=COL_CLOSE, mid="white", high=COL_OPEN,
                             midpoint=0, name="NES\n(capped ±4)") +
        scale_size_continuous(range=c(2,12), name="-log10(p)",
                              breaks=c(0.5,1,2,3)) +
        scale_x_discrete(labels=function(x) gsub(" ","\n",x)) +
        labs(title="ATAC fate \u00d7 RNA contrast enrichment",
             subtitle=". p<0.10  * p<0.05  ** p<0.01  (permutation, uncorrected)",
             x=NULL, y="ATAC fate gene set") +
        theme_classic(base_size=11) +
        theme(axis.text.x=element_text(angle=30, hjust=1, size=9),
              axis.text.y=element_text(size=9),
              plot.title=element_text(face="bold", size=12),
              plot.subtitle=element_text(size=9, colour="#555"),
              legend.position="right")

      save_plot(file.path(OUT_DIR,"part2_enrichment","enrichment_overview"),
                p_bubble, w=12, h=7)

      # ── Ranked bar: top enrichments by NES × significance ────────────────
      top_enrich <- enrich_df %>%
        filter(pval < 0.20) %>%
        mutate(label=sprintf("%s\n\u00d7 %s", fate_display, contrast_short)) %>%
        arrange(NES)

      if (nrow(top_enrich) > 0) {
        p_top <- ggplot(top_enrich,
                        aes(x=reorder(label,NES), y=NES, fill=NES>0)) +
          geom_col(width=0.7) +
          geom_text(aes(label=sprintf("n=%d, p=%.3f", n_set, pval)),
                    hjust=ifelse(top_enrich$NES>0,-0.05,1.05), size=2.8) +
          scale_fill_manual(values=c("TRUE"=COL_OPEN,"FALSE"=COL_CLOSE),
                            guide="none") +
          coord_flip() +
          geom_hline(yintercept=0, colour="#555", linewidth=0.5) +
          labs(title="Top enrichments (p<0.20)", x=NULL, y="NES") +
          theme_classic(base_size=10) +
          theme(plot.title=element_text(face="bold"))
        save_plot(file.path(OUT_DIR,"part2_enrichment","enrichment_top_ranked"),
                  p_top, w=11, h=max(4, nrow(top_enrich)*0.45+2))
      }
    }
  }
}

# =============================================================================
# PART 3 — CLUSTER × FATE ENRICHMENT
# =============================================================================

if (3 %in% run_parts) {
  cat("\n", strrep("=",70), "\n", sep="")
  cat("PART 3: CLUSTER × FATE ENRICHMENT\n\n")

  clust_path <- file.path(CLUSTER_DIR, "cluster_assignments.tsv")
  if (!file.exists(clust_path)) {
    cat("[integrate]   cluster_assignments.tsv not found\n")
  } else {
    clust_df <- tryCatch(
      read.table(clust_path, header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE),
      error=function(e) NULL)
    fate_gene_lists <- read_fate_gene_lists()
    nearby_variants <- read_fate_gene_lists_nearby_expanded()  # [v16] for 3m

    if (!is.null(clust_df) && length(fate_gene_lists) > 0) {
      sym_col <- intersect(c("symbol","Symbol","SYMBOL"), colnames(clust_df))[1]
      if (is.na(sym_col)) {
        cat("[integrate]   No symbol column in cluster assignments\n")
      } else {
        all_syms    <- unique(clust_df[[sym_col]])
        cluster_ids <- sort(unique(clust_df$cluster))
        fate_names  <- names(fate_gene_lists)
        fate_short  <- gsub("^genes_","", fate_names)

        # Reusable Fisher grid builder. Called twice below: once on the
        # fine-grained fate taxonomy (for the composition plots, 3d/3h,
        # which only report % and don't need testing power) and once on a
        # genotype-collapsed taxonomy (for every plot that reports
        # significance -- 3a/3b/3c/3e -- plus cl_order). Splitting Converged
        # into 6 mechanisms x 2 genotypes fragments what were often modest
        # gene sets into slices too thin for Fisher's exact test to ever
        # call significant; collapsing per genotype before testing roughly
        # halves the test family (fewer BH comparisons) and thickens each
        # remaining cell (more power), without touching the fine labels
        # used for describing composition.
        run_fisher_grid <- function(gene_lists) {
          nm    <- names(gene_lists)
          short <- gsub("^genes_","", nm)
          rows  <- list()
          for (cl_id in cluster_ids) {
            cl_label   <- unique(clust_df$cluster_label[clust_df$cluster==cl_id])[1]
            in_cluster <- clust_df[[sym_col]][clust_df$cluster==cl_id]
            for (fi in seq_along(nm)) {
              fate <- nm[fi]
              in_fate <- intersect(gene_lists[[fate]], all_syms)
              a  <- length(intersect(in_cluster, in_fate))
              b  <- length(setdiff(in_cluster, in_fate))
              cc <- length(setdiff(in_fate, in_cluster))
              d  <- length(all_syms) - a - b - cc
              if (a == 0 || d < 0) next
              ft <- tryCatch(fisher.test(matrix(c(a,b,cc,d),2,2), alternative="greater"),
                             error=function(e) NULL)
              if (is.null(ft)) next
              rows[[length(rows)+1]] <- data.frame(
                cluster=cl_id, cluster_label=cl_label,
                fate=fate, fate_short=short[fi],
                n_cluster=length(in_cluster), n_fate=length(in_fate),
                n_overlap=a, pct_overlap=round(100*a/length(in_cluster),1),
                odds_ratio=round(ft$estimate,3), pval=ft$p.value,
                stringsAsFactors=FALSE)
            }
          }
          if (length(rows) == 0) return(NULL)
          df <- do.call(rbind, rows)
          # [diagnostic] geno_stratum groups fates by genotype_bias for the
          # per-genotype stratified-BH correction below. fate_short groups
          # by the individual fate category itself for an even finer,
          # per-fate correction: with padj_fate_stratified, growing one
          # fate category's gene set (more genes -> more tested clusters
          # for THAT fate) can never change another fate category's own
          # padj, because each fate's correction family is just its own
          # (up to 13) cluster tests, fully decoupled from every other
          # fate's size. This is the most literal answer to "more genes
          # should only ever ADD tests, never take significance away from
          # other tests" -- global and per-genotype correction both still
          # pool multiple fate categories into one shared budget, so a
          # category that grows can still dilute a DIFFERENT category's
          # padj through that shared budget; per-fate correction can't,
          # by construction. The real cost: with as few as 1-13 tests per
          # family, BH's actual FDR-control guarantee is weak (there's
          # barely a multiple-testing problem left to correct for), so
          # padj_fate_stratified trades rigor for sensitivity -- report it
          # as "recovered candidates," not as a stronger significance
          # claim than the global or per-genotype numbers.
          df$geno_stratum <- fate_genotype_bias(df$fate_short)
          df %>%
            group_by(geno_stratum) %>%
            mutate(padj_stratified = p.adjust(pval, "BH"),
                   sig_stratified  = sig_stars(padj_stratified)) %>%
            ungroup() %>%
            group_by(fate_short) %>%
            mutate(padj_fate_stratified = p.adjust(pval, "BH"),
                   sig_fate_stratified  = sig_stars(padj_fate_stratified)) %>%
            ungroup() %>%
            mutate(padj           = p.adjust(pval, "BH"),
                   neg_lp         = -log10(pmax(padj,1e-6)),
                   log2_OR        = log2(pmax(odds_ratio,0.01)),
                   sig            = sig_stars(padj),
                   cl_label_short = sprintf("C%d: %s", cluster, cluster_label))
        }

        # Union the 6 Converged mechanisms into one gene set per genotype.
        # Fine mechanism-level labels (e.g. "BOTv_Converged_both_shift")
        # stay untouched in fate_gene_lists itself / the fine table below;
        # this only builds a second, coarser list for the testing grid.
        #
        # [bugfix] read_fate_gene_lists() names its list entries from the
        # gene_lists/*.csv FILENAMES with only the ".csv" suffix stripped
        # (`names(out) <- gsub("\\.csv$","", basename(files))`) -- it does
        # NOT strip the "genes_" prefix that 01_peak_fate_classification.r's
        # save_gene_list() puts on every filename. So the raw keys reaching
        # this function look like "genes_BOTCv_Converged_BOTCv_closes", not
        # "BOTCv_Converged_BOTCv_closes". The original regex here
        # (^BOTv_Converged_ / ^BOTCv_Converged_) anchored to the true start
        # of the string and so never matched ANYTHING -- collapse_converged()
        # silently returned its input unchanged every time it was called,
        # meaning every "collapsed" significance test in this script (3a,
        # 3b, 3c, 3e, 3j, 3k, 3m, the stratified dotplot) has actually been
        # running on the full ~20-category fragmented taxonomy the whole
        # time, not the intended ~10-category collapsed one. Stripping the
        # optional "genes_" prefix before matching is the fix.
        collapse_converged <- function(gene_lists) {
          nm         <- sub("^genes_", "", names(gene_lists))
          conv_botv  <- names(gene_lists)[grepl("^BOTv_Converged_",  nm)]
          conv_botcv <- names(gene_lists)[grepl("^BOTCv_Converged_", nm)]
          out <- gene_lists[setdiff(names(gene_lists), c(conv_botv, conv_botcv))]
          if (length(conv_botv))
            out[["BOTv_Converged"]]  <- unique(unlist(gene_lists[conv_botv],  use.names=FALSE))
          if (length(conv_botcv))
            out[["BOTCv_Converged"]] <- unique(unlist(gene_lists[conv_botcv], use.names=FALSE))
          out
        }

        fish_df     <- run_fisher_grid(fate_gene_lists)                     # fine — composition only (3d/3h/3i)
        fish_sig_df <- run_fisher_grid(collapse_converged(fate_gene_lists)) # collapsed — all significance plots

        if (!is.null(fish_df) && !is.null(fish_sig_df)) {

          write.table(fish_sig_df,
                      file.path(OUT_DIR,"part3_cluster_fate",
                                "cluster_fate_enrichment_table.tsv"),
                      sep="\t", quote=FALSE, row.names=FALSE)
          write.table(fish_df,
                      file.path(OUT_DIR,"part3_cluster_fate",
                                "cluster_fate_composition_table_fine.tsv"),
                      sep="\t", quote=FALSE, row.names=FALSE)

          max_OR <- max(abs(fish_sig_df$log2_OR), na.rm=TRUE)
          cl_order <- fish_sig_df %>%
            group_by(cl_label_short) %>%
            summarise(max_sig=max(neg_lp)) %>%
            arrange(desc(max_sig)) %>%
            pull(cl_label_short)

          # ── 3a: Filled tile heatmap (OR) with significance stars ──────────
          # Uses fish_sig_df (Converged mechanisms collapsed per genotype)
          p_heat <- ggplot(fish_sig_df,
                           aes(x=fate_short,
                               y=factor(cl_label_short, levels=rev(cl_order)),
                               fill=log2_OR)) +
            geom_tile(colour="white", linewidth=0.6) +
            geom_text(aes(label=sig), size=5, colour="black", vjust=0.7) +
            scale_fill_gradient2(low=COL_CLOSE, mid="white", high=COL_OPEN,
                                 midpoint=0, name="log2(OR)",
                                 limits=c(-max_OR, max_OR)) +
            scale_x_discrete(labels=pi_fate_label) +
            labs(title="Expression cluster × ATAC fate: Fisher's odds ratio",
                 subtitle=". p<0.10  * p<0.05  ** p<0.01  *** p<0.001  (BH-corrected)",
                 x="ATAC fate gene set", y="Expression cluster") +
            theme_classic(base_size=11) +
            theme(axis.text.x=element_text(angle=30, hjust=1, size=9),
                  axis.text.y=element_text(size=9),
                  plot.title=element_text(face="bold", size=12),
                  plot.subtitle=element_text(size=9, colour="#555"),
                  panel.grid=element_blank())

          # ── 3b: Significance heatmap (−log10 padj) ───────────────────────
          p_sig <- ggplot(fish_sig_df,
                          aes(x=fate_short,
                              y=factor(cl_label_short, levels=rev(cl_order)),
                              fill=neg_lp)) +
            geom_tile(colour="white", linewidth=0.6) +
            geom_text(aes(label=sig), size=5, colour="white", vjust=0.7) +
            scale_fill_gradient(low="white", high="#c9579a",
                                name="-log10\n(adj.p)") +
            scale_x_discrete(labels=pi_fate_label) +
            labs(title="Significance (−log10 BH-padj)",
                 x="ATAC fate gene set", y=NULL) +
            theme_classic(base_size=11) +
            theme(axis.text.x=element_text(angle=30, hjust=1, size=9),
                  axis.text.y=element_blank(), axis.ticks.y=element_blank())

          p_dual <- p_heat | p_sig
          save_plot(file.path(OUT_DIR,"part3_cluster_fate",
                              "cluster_fate_heatmap"),
                    p_dual,
                    w=max(14, length(unique(fish_sig_df$fate_short))*1.4 + 8),
                    h=max(6,  length(cluster_ids)*0.55 + 2.5))

          # ── 3c: Top enriched pairs — lollipop chart ───────────────────────
          top_fish <- fish_sig_df %>%
            filter(padj < 0.20) %>%
            mutate(label=sprintf("%s\n+ %s", cl_label_short, pi_fate_label(fate_short))) %>%
            arrange(desc(log2_OR))

          if (nrow(top_fish) > 0) {
            p_lollipop <- ggplot(top_fish,
                                 aes(x=reorder(label, log2_OR), y=log2_OR,
                                     colour=log2_OR>0)) +
              geom_segment(aes(xend=reorder(label,log2_OR), yend=0),
                           linewidth=1.2) +
              geom_point(aes(size=neg_lp)) +
              geom_text(aes(label=sprintf("n=%d (%.0f%%)\n%s",
                                          n_overlap, pct_overlap, sig)),
                        hjust=ifelse(top_fish$log2_OR>0,-0.1,1.1),
                        size=2.8) +
              scale_colour_manual(values=c("TRUE"=COL_OPEN,"FALSE"=COL_CLOSE),
                                  guide="none") +
              scale_size_continuous(range=c(2,8), name="-log10\n(adj.p)") +
              coord_flip() +
              geom_hline(yintercept=0, colour="#555", linewidth=0.5) +
              labs(title="Top cluster-fate enrichments (padj<0.20)",
                   subtitle="Dot size = significance | Label = n overlap (%) and stars",
                   x=NULL, y="log2(Odds Ratio)") +
              theme_classic(base_size=10) +
              theme(plot.title=element_text(face="bold"),
                    plot.subtitle=element_text(size=8.5, colour="#555"))
            save_plot(file.path(OUT_DIR,"part3_cluster_fate",
                                "cluster_fate_top_enrichments"),
                      p_lollipop,
                      w=12, h=max(4, nrow(top_fish)*0.55+2.5))
          }

          # ── 3d: Fate composition per cluster (stacked bar) ────────────────
          # Build colour map: assign each fate_short a colour. Active fates
          # (Deepened/Maintained/Reversed) get bold, genotype-family colours
          # (BOTv=green, BOTCv=magenta, matching the rest of the pipeline).
          # Converged sub-mechanisms (6 per genotype: no_temporal_signal,
          # BOTv_opens, BOTv_closes, BOTCv_opens, BOTCv_closes, both_shift)
          # used to all collapse to one flat grey, making them impossible to
          # tell apart in the stacked bar. Instead, each mechanism gets a
          # distinct position along a light->dark ramp in its genotype's hue,
          # ordered from least informative (no_temporal_signal, lightest) to
          # most informative (both_shift, darkest) -- so ramp position
          # carries meaning and adjacent Converged segments stay separable.
          all_fate_shorts <- unique(fish_df$fate_short)

          botv_ramp  <- colorRampPalette(c("#d9ecd2", "#1a5c1a"))(8)
          botcv_ramp <- colorRampPalette(c("#f6d9ec", "#8B0057"))(8)
          conv_mechs <- c("no_temporal_signal","BOTv_opens","BOTv_closes",
                          "BOTCv_opens","BOTCv_closes","both_shift")

          fate_col_map <- vapply(all_fate_shorts, function(fs) {
            # fate_genotype_bias() correctly inverts Reversed categories to
            # their ENDING genotype (see its definition/comment near the top
            # of the file) -- using it here, instead of the old naive
            # grepl("^BOTv"/"^BOTCv", fs), is the fix for Reversed fates
            # previously being coloured by their starting genotype's ramp.
            genotype <- fate_genotype_bias(fs)
            ramp <- if (identical(genotype,"BOTv")) botv_ramp else botcv_ramp

            if (grepl("^Emerged", fs))    return("#F39C12")
            if (grepl("Deepened",   fs))  return(ramp[8])   # darkest = strongest active signal
            if (grepl("Reversed",   fs))  return(ramp[7])
            if (grepl("Maintained", fs))  return(ramp[5])
            if (grepl("Converged", fs)) {
              hit <- conv_mechs[vapply(conv_mechs, function(m) grepl(m, fs), logical(1))]
              pos <- if (length(hit)) match(hit[1], conv_mechs) else 1
              # spread the 6 mechanisms across ramp positions 1-6 (lightest->darker),
              # leaving 7-8 reserved for Reversed/Deepened so they stay visually "boldest"
              return(ramp[pos])
            }
            "#999999"
          }, character(1))
          names(fate_col_map) <- all_fate_shorts

          fate_comp <- fish_df %>%
            mutate(pct_of_cluster = 100*n_overlap/n_cluster)

          p_stacked <- ggplot(fate_comp,
                              aes(x=factor(cl_label_short, levels=cl_order),
                                  y=pct_of_cluster, fill=fate_short)) +
            geom_col(position="stack", width=0.75) +
            geom_text(aes(label=ifelse(pct_of_cluster>=2,
                                       paste0(round(pct_of_cluster,1),"%"),"")),
                      position=position_stack(vjust=0.5),
                      size=2.5, colour="white", fontface="bold") +
            scale_fill_manual(values=fate_col_map, labels=pi_fate_label, name="ATAC fate") +
            scale_y_continuous(expand=expansion(mult=c(0.02,0.40))) +
          coord_flip(clip="off") +
          theme(plot.margin=margin(5,160,5,5,"pt")) +
            labs(title="ATAC fate composition per expression cluster",
                 subtitle="% of cluster genes overlapping each fate gene set",
                 x="Expression cluster", y="% of cluster genes") +
            theme_classic(base_size=11) +
            theme(plot.title=element_text(face="bold"), legend.position="right")
          save_plot(file.path(OUT_DIR,"part3_cluster_fate",
                              "cluster_fate_composition"),
                    p_stacked,
                    w=12, h=max(5, length(cluster_ids)*0.5+2))

          cat("\n[integrate]   Top enrichments (padj<0.20):\n")
          top_print <- fish_df %>% filter(padj<0.20) %>%
            dplyr::select(cluster_label, fate_short, n_overlap, pct_overlap,
                          odds_ratio, pval, padj, sig) %>%
            arrange(padj)
          print(as.data.frame(top_print), row.names=FALSE)

          # ── 3e: Temporal-coloured dotplot with proper legend ─────────────
          # Significant enrichments (padj<0.20) coloured by genotype identity
          # of each ATAC fate: dark magenta=BOTCv, dark green=BOTv, orange=Emerging.
          # Uses scale_fill_manual with named values so the legend is built
          # automatically — no floating annotation text needed.
          # Pin to the collapsed significance-testing universe (fish_sig_df's
          # own fate_short values), NOT the fine 16-category composition
          # universe -- the fine per-mechanism Converged categories were
          # never tested here, so pinning to them would draw a dozen
          # permanently-empty columns that don't correspond to any test.
          all_fate_shorts_sig <- unique(fish_sig_df$fate_short)

          fish_sig <- fish_sig_df %>%
            filter(padj < 0.20) %>%
            mutate(
              fate_label_col = fate_display_col(fate_short),
              # ggplot2's discrete scales default to drop=TRUE, which
              # silently deletes any factor level absent from the plotted
              # rows even when declared -- so pin fate_short to the full
              # collapsed universe here, and pair with drop=FALSE on the
              # scales below, or the padj filter would remove whole axis
              # columns instead of just points.
              fate_short = factor(fate_short, levels = all_fate_shorts_sig),
              # fate_genotype_bias() correctly inverts Reversed categories
              # to their ENDING genotype -- see its definition near the top
              # of the file for why the old naive grepl("BOTCv"/"BOTv",
              # fate_short) here coloured every Reversed fate by its
              # STARTING genotype instead (dots were flipped: "Reversed to
              # BOTC.v" rendered dark green instead of dark pink, and vice
              # versa).
              geno_group = dplyr::case_when(
                fate_genotype_bias(as.character(fate_short)) == "BOTCv"   ~ "BOTCv fate",
                fate_genotype_bias(as.character(fate_short)) == "BOTv"    ~ "BOTv fate",
                fate_genotype_bias(as.character(fate_short)) == "Emerged" ~ "Emerging",
                TRUE ~ "Other"))

          if (nrow(fish_sig) > 0) {
            geno_cols <- c(
              "BOTCv fate" = COL_BOTCV_LATE,
              "BOTv fate"  = COL_BOTV_LATE,
              "Emerging"   = "#F39C12",
              "Other"      = "#888888")

            # Compute nice size breaks from actual data range
            max_nlp  <- max(fish_sig$neg_lp, na.rm=TRUE)
            sz_breaks <- unique(round(c(
              1,
              if (max_nlp >= 1.5) 1.5,
              if (max_nlp >= 2)   2,
              if (max_nlp >= 3)   3), 1))
            sz_breaks <- sz_breaks[sz_breaks <= max_nlp]

            p_tdot <- ggplot(fish_sig,
                              aes(x=fate_short,
                                  y=factor(cl_label_short, levels=rev(cl_order)),
                                  size=neg_lp,
                                  fill=geno_group)) +
              geom_point(shape=21, colour="white", alpha=0.92) +
              geom_text(aes(label=sig), colour="black", size=3.5, vjust=0.4) +
              scale_size_continuous(
                range=c(2,10),
                name="-log10\n(adj.p)",
                breaks=sz_breaks,
                labels=sprintf("%.1f", sz_breaks),
                guide=guide_legend(
                  title="-log10\n(adj.p)",
                  override.aes=list(fill="#888888", colour="white"),
                  order=2)) +
              scale_fill_manual(
                values=geno_cols,
                name="ATAC fate\ngenotype",
                labels=c(
                  "BOTCv fate" = "BOTC.v-biased",
                  "BOTv fate"  = "BOT.v-biased",
                  "Emerging"   = "Emerged",
                  "Other"      = "Other"),
                guide=guide_legend(
                  title="ATAC fate\ngenotype",
                  override.aes=list(size=4, colour="white"),
                  order=1)) +
              scale_x_discrete(labels=pi_fate_label, drop=FALSE) +
              scale_y_discrete(drop=FALSE) +
              labs(title="Cluster \u00d7 ATAC fate: significant enrichments",
                   subtitle=". p<0.10  * p<0.05  ** p<0.01  *** p<0.001 (BH-corrected)\nDark pink = BOTC.v-biased fates  |  Dark green = BOT.v-biased fates  |  Orange = Emerged",
                   x="ATAC fate", y="Expression cluster") +
              theme_classic(base_size=10) +
              theme(axis.text.x=element_text(angle=30,hjust=1,size=8),
                    axis.text.y=element_text(size=8.5),
                    plot.title=element_text(face="bold",size=11),
                    plot.subtitle=element_text(size=8.5,colour="#444"),
                    legend.title=element_text(size=9),
                    legend.text=element_text(size=8.5),
                    legend.spacing.y=unit(0.3,"cm"))
            # Sizing uses the full declared (collapsed) axis universe, not
            # just the significant subset, since drop=FALSE keeps every fate
            # column and every cluster row visible regardless of significance.
            save_plot(file.path(OUT_DIR,"part3_cluster_fate",
                                "cluster_fate_temporal_dotplot"),
                      p_tdot,
                      w=max(9, length(all_fate_shorts_sig)*1.4+3),
                      h=max(5, length(cluster_ids)*0.52+2.5))

            # [requested] Compact companion version: same data, but with
            # ggplot2's default drop=TRUE behaviour restored so unused
            # fate/cluster axis levels collapse away instead of rendering as
            # empty rows/columns -- the pre-fix look, kept as an option
            # alongside the full-grid version above (some audiences want
            # the complete grid for direct comparison against 3h; others
            # just want the compact view of the hits themselves).
            p_tdot_compact <- suppressMessages(
              p_tdot +
                scale_x_discrete(labels=pi_fate_label, drop=TRUE) +
                scale_y_discrete(drop=TRUE) +
                labs(subtitle=paste0(
                  ". p<0.10  * p<0.05  ** p<0.01  *** p<0.001 (BH-corrected)\n",
                  "Dark pink = BOTC.v-biased fates  |  Dark green = BOT.v-biased fates  |  Orange = Emerged\n",
                  "(compact view: empty rows/columns dropped -- see the full-grid version for the complete axis)")))
            n_fate_compact <- length(unique(as.character(fish_sig$fate_short)))
            n_cl_compact   <- length(unique(fish_sig$cl_label_short))
            save_plot(file.path(OUT_DIR,"part3_cluster_fate",
                                "cluster_fate_temporal_dotplot_compact"),
                      p_tdot_compact,
                      w=max(7, n_fate_compact*1.4+3),
                      h=max(4, n_cl_compact*0.52+2.5))

            # [requested] Stratified-significance companion versions: same
            # data, but padj computed within a smaller, more decoupled
            # family than the shared global budget (see run_fisher_grid's
            # padj_stratified / padj_fate_stratified, and 3n's diagnostic
            # further down). If one fate's growth is crowding out another's
            # significance under the shared budget, this is where that
            # recovered signal becomes visible -- a principled way to
            # "loosen" thresholding that's auditable (every recovered point
            # is still BH-corrected, just within a smaller family) rather
            # than just picking a bigger padj cutoff. Two variants, in
            # order of how decoupled the correction family is:
            #   - genotype-stratified: corrects within BOTv fates / BOTCv
            #     fates / Emerged separately (growth in one genotype arm
            #     can't dilute the other's budget, but fates WITHIN the
            #     same arm still share one)
            #   - fate-stratified: corrects within each INDIVIDUAL fate
            #     category's own cluster tests (growth in any other fate,
            #     same arm or not, can never dilute this one) -- the most
            #     literal "additive, not subtractive" reading, at the cost
            #     of very small per-family test counts where BH barely
            #     does any correcting at all; treat its hits as candidates,
            #     not confirmed findings.
            render_strat_dotplot <- function(padj_col, sig_col, neg_lp_label,
                                             out_stem, budget_desc) {
              fs <- fish_sig_df %>%
                filter(.data[[padj_col]] < 0.20) %>%
                mutate(
                  fate_label_col = fate_display_col(fate_short),
                  fate_short = factor(fate_short, levels = all_fate_shorts_sig),
                  geno_group = dplyr::case_when(
                    fate_genotype_bias(as.character(fate_short)) == "BOTCv"   ~ "BOTCv fate",
                    fate_genotype_bias(as.character(fate_short)) == "BOTv"    ~ "BOTv fate",
                    fate_genotype_bias(as.character(fate_short)) == "Emerged" ~ "Emerging",
                    TRUE ~ "Other"),
                  neg_lp_strat = -log10(pmax(.data[[padj_col]],1e-6)),
                  sig_strat    = .data[[sig_col]])
              if (nrow(fs) == 0) return(invisible(NULL))
              n_recovered <- length(setdiff(
                paste(fs$cluster, fs$fate_short),
                paste(fish_sig$cluster, fish_sig$fate_short)))
              p <- ggplot(fs,
                          aes(x=fate_short,
                              y=factor(cl_label_short, levels=rev(cl_order)),
                              size=neg_lp_strat, fill=geno_group)) +
                geom_point(shape=21, colour="white", alpha=0.92) +
                geom_text(aes(label=sig_strat), colour="black", size=3.5, vjust=0.4) +
                scale_size_continuous(range=c(2,10), name=neg_lp_label) +
                scale_fill_manual(
                  values=geno_cols, name="ATAC fate\ngenotype",
                  labels=c("BOTCv fate"="BOTC.v-biased",
                          "BOTv fate" ="BOT.v-biased",
                          "Emerging"  ="Emerged", "Other"="Other")) +
                scale_x_discrete(labels=pi_fate_label, drop=FALSE) +
                scale_y_discrete(drop=FALSE) +
                labs(title=sprintf("Cluster \u00d7 ATAC fate: significant enrichments (%s)", budget_desc),
                     subtitle=paste0(
                       ". p<0.10  * p<0.05  ** p<0.01  *** p<0.001 (BH-corrected ", budget_desc, ")\n",
                       sprintf("%d pair(s) recovered here that are NOT significant under the shared global-BH dotplot above",
                              n_recovered)),
                     x="ATAC fate", y="Expression cluster") +
                theme_classic(base_size=10) +
                theme(axis.text.x=element_text(angle=30,hjust=1,size=8),
                      axis.text.y=element_text(size=8.5),
                      plot.title=element_text(face="bold",size=11),
                      plot.subtitle=element_text(size=8,colour="#444"))
              save_plot(file.path(OUT_DIR,"part3_cluster_fate", out_stem), p,
                        w=max(9, length(all_fate_shorts_sig)*1.4+3),
                        h=max(5, length(cluster_ids)*0.52+2.5))
              cat(sprintf(paste0("[integrate]   %s dotplot: %d total significant pairs ",
                                 "(%d recovered vs the shared global-BH budget)\n"),
                          budget_desc, nrow(fs), n_recovered))
            }

            render_strat_dotplot("padj_stratified", "sig_stratified",
                                 "-log10\n(adj.p)\nper-genotype",
                                 "cluster_fate_temporal_dotplot_stratified",
                                 "WITHIN each genotype arm")
            render_strat_dotplot("padj_fate_stratified", "sig_fate_stratified",
                                 "-log10\n(adj.p)\nper-fate",
                                 "cluster_fate_temporal_dotplot_fate_stratified",
                                 "WITHIN each individual fate category")
          }

          # ── 3f: Net BOTv- vs BOTCv-driven chromatin dynamics per cluster ──
          # Collapses every BOTv-prefixed fate_short into one "BOTv" total and
          # every BOTCv-prefixed one into a "BOTCv" total (Emerged excluded --
          # it isn't genotype-specific), then plots the signed difference.
          # This answers a different question than 3d/3e: not "which exact
          # mechanism" but "does this expression cluster's chromatin story
          # skew toward BOTv or BOTCv overall?"
          if (!has_tidyr) {
            cat("[integrate]   Skipping 3f (net genotype plot): tidyr not installed\n")
          } else {
          geno_net_df <- fate_comp %>%
            mutate(genotype = dplyr::case_when(
              grepl("^BOTv",  fate_short) ~ "BOTv",
              grepl("^BOTCv", fate_short) ~ "BOTCv",
              TRUE ~ NA_character_)) %>%
            filter(!is.na(genotype)) %>%
            group_by(cl_label_short, genotype) %>%
            summarise(pct_total = sum(pct_of_cluster), .groups="drop") %>%
            tidyr::pivot_wider(names_from=genotype, values_from=pct_total,
                               values_fill=0) %>%
            mutate(net = BOTv - BOTCv)

          p_net <- ggplot(geno_net_df,
                          aes(x=reorder(cl_label_short, net), y=net,
                              fill=net>0)) +
            geom_col(width=0.7) +
            geom_text(aes(label=sprintf("%+.1f%%", net)),
                      hjust=ifelse(geno_net_df$net>0,-0.1,1.1), size=3) +
            geom_hline(yintercept=0, colour="#555", linewidth=0.5) +
            scale_fill_manual(values=c("TRUE"=COL_BOTV_LATE,"FALSE"=COL_BOTCV_LATE),
                              guide="none") +
            coord_flip(clip="off") +
            theme(plot.margin=margin(5,45,5,5,"pt")) +
            labs(title="Net genotype-driven chromatin dynamics per expression cluster",
                 subtitle="Positive (green) = cluster's ATAC fate genes skew BOT.v-biased | Negative (pink) = skew BOTC.v-biased",
                 x="Expression cluster",
                 y="Net % of cluster genes (\u03a3 BOTv-classified \u2212 \u03a3 BOTCv-classified)") +
            theme_classic(base_size=11) +
            theme(plot.title=element_text(face="bold"),
                  plot.subtitle=element_text(size=8.5,colour="#555"))
          save_plot(file.path(OUT_DIR,"part3_cluster_fate",
                              "cluster_fate_net_genotype"),
                    p_net,
                    w=10, h=max(5, length(unique(geno_net_df$cl_label_short))*0.5+2))
          }

          # ── 3g: Active vs Converged fraction per cluster ───────────────────
          # "Active" = Deepened/Maintained/Reversed/Emerged (locus still has a
          # live BOTv-vs-BOTCv difference at nc14late). "Converged" = the
          # difference disappeared. Answers: how much of each expression
          # cluster's chromatin signal is still "live" vs already resolved?
          active_conv_df <- fate_comp %>%
            mutate(status = ifelse(grepl("Converged", fate_short),
                                   "Converged", "Active")) %>%
            group_by(cl_label_short, status) %>%
            summarise(pct_total = sum(pct_of_cluster), .groups="drop") %>%
            mutate(cl_label_short = factor(cl_label_short, levels=cl_order))

          p_active <- ggplot(active_conv_df,
                             aes(x=cl_label_short, y=pct_total, fill=status)) +
            geom_col(position="stack", width=0.7) +
            geom_text(aes(label=paste0(round(pct_total,1),"%")),
                      position=position_stack(vjust=0.5),
                      size=3, colour="white", fontface="bold") +
            scale_fill_manual(values=c(Active="#c9579a", Converged="#AAAAAA"),
                              name="ATAC fate status") +
            coord_flip() +
            labs(title="Active vs Converged ATAC-fate signal per expression cluster",
                 subtitle="Active = Deepened/Maintained/Reversed/Emerged (still a BOTv-vs-BOTCv difference at nc14late)  |  Converged = difference resolved",
                 x="Expression cluster", y="% of cluster genes") +
            theme_classic(base_size=11) +
            theme(plot.title=element_text(face="bold"),
                  plot.subtitle=element_text(size=8.5,colour="#555"),
                  legend.position="right")
          save_plot(file.path(OUT_DIR,"part3_cluster_fate",
                              "cluster_fate_active_vs_converged"),
                    p_active,
                    w=10, h=max(5, length(unique(active_conv_df$cl_label_short))*0.5+2))

          # ── 3h: Hierarchically clustered fate-composition heatmap ─────────
          # 3a/3b order clusters by enrichment significance; this instead
          # reorders BOTH clusters and fate categories by similarity of their
          # composition profiles, surfacing groups of expression clusters that
          # share an ATAC-fate "fingerprint" even if no single cell is
          # individually significant.
          if (!has_tidyr) {
            cat("[integrate]   Skipping 3h (hierarchical heatmap): tidyr not installed\n")
          } else {
          comp_mat <- fate_comp %>%
            dplyr::select(cl_label_short, fate_short, pct_of_cluster) %>%
            tidyr::pivot_wider(names_from=fate_short, values_from=pct_of_cluster,
                               values_fill=0) %>%
            as.data.frame()
          rownames(comp_mat) <- comp_mat$cl_label_short
          comp_mat$cl_label_short <- NULL
          comp_mat <- as.matrix(comp_mat)

          if (nrow(comp_mat) >= 3 && ncol(comp_mat) >= 3) {
            row_hc <- hclust(dist(comp_mat), method="complete")
            col_hc <- hclust(dist(t(comp_mat)), method="complete")
            row_order_hc <- rownames(comp_mat)[row_hc$order]
            col_order_hc <- colnames(comp_mat)[col_hc$order]

            fate_comp_hc <- fate_comp %>%
              mutate(cl_label_short=factor(cl_label_short, levels=row_order_hc),
                     fate_short=factor(fate_short, levels=col_order_hc))

            p_hcheat <- ggplot(fate_comp_hc,
                               aes(x=fate_short, y=cl_label_short,
                                   fill=pct_of_cluster)) +
              geom_tile(colour="white", linewidth=0.6) +
              geom_text(aes(label=ifelse(pct_of_cluster>=2,
                                         round(pct_of_cluster,1),"")),
                        size=2.6, colour="#333") +
              scale_fill_gradient(low="white", high="#4a2545", name="% of\ncluster") +
              scale_x_discrete(labels=pi_fate_label) +
              labs(title="ATAC fate composition, hierarchically clustered",
                   subtitle="Rows and columns reordered by similarity (Euclidean distance, complete linkage) -- reveals clusters sharing a fate fingerprint",
                   x="ATAC fate", y="Expression cluster") +
              theme_classic(base_size=10) +
              theme(axis.text.x=element_text(angle=30,hjust=1,size=8),
                    axis.text.y=element_text(size=8.5),
                    plot.title=element_text(face="bold",size=12),
                    plot.subtitle=element_text(size=8.5,colour="#555"),
                    panel.grid=element_blank())
            save_plot(file.path(OUT_DIR,"part3_cluster_fate",
                                "cluster_fate_heatmap_clustered"),
                      p_hcheat,
                      w=max(11, length(col_order_hc)*1.1+4),
                      h=max(5, length(row_order_hc)*0.5+2.5))
          } else {
            cat("[integrate]   Skipping 3h (hierarchical heatmap): need >=3 clusters and >=3 fate categories\n")
          }
          }

          # ── 3i: Converged-mechanism-only breakdown per cluster ─────────────
          # Zooms into just the Converged rows -- the majority category in 3d
          # -- and asks *why* each cluster's converged loci converged. Uses
          # the same mechanism colour set as plot_peak_fate_v11.r's own
          # Converged-breakdown plots, so this stays visually consistent with
          # the ATAC-side figures. Faceted by which genotype's ATAC calls
          # produced the classification, since BOTv_Converged_* and
          # BOTCv_Converged_* are independent gene sets, not a single split.
          conv_mech_cols <- c(
            "BOTv_opens"          = "#7fb3a0",
            "BOTv_closes"         = "#a68a5b",
            "both_shift"          = "#8067b7",
            "BOTCv_opens"         = "#d99cc0",
            "BOTCv_closes"        = "#8a6a7a",
            "no_temporal_signal"  = "#d9d9d9")
          conv_mech_labels <- c(
            "BOTv_opens"="BOTv opens", "BOTv_closes"="BOTv closes",
            "both_shift"="Both shift", "BOTCv_opens"="BOTCv opens",
            "BOTCv_closes"="BOTCv closes", "no_temporal_signal"="No temporal signal")

          conv_only_df <- fate_comp %>%
            filter(grepl("Converged", fate_short)) %>%
            mutate(
              genotype_prefix = ifelse(grepl("^BOTv", fate_short),
                                       "BOT.v-biased loci", "BOTC.v-biased loci"),
              mechanism = vapply(fate_short, function(fs) {
                hit <- names(conv_mech_cols)[vapply(names(conv_mech_cols),
                                                     function(m) grepl(m, fs), logical(1))]
                if (length(hit)) hit[1] else "no_temporal_signal"
              }, character(1)),
              cl_label_short = factor(cl_label_short, levels=cl_order))

          if (nrow(conv_only_df) > 0) {
            p_conv <- ggplot(conv_only_df,
                             aes(x=cl_label_short, y=pct_of_cluster, fill=mechanism)) +
              geom_col(position="stack", width=0.7) +
              facet_wrap(~genotype_prefix) +
              coord_flip() +
              scale_fill_manual(values=conv_mech_cols, labels=conv_mech_labels,
                                name="Convergence\nmechanism") +
              labs(title="Converged-locus mechanism breakdown per expression cluster",
                   subtitle="Why the BOTv-vs-BOTCv chromatin difference disappeared, split by which genotype's ATAC calls produced the classification",
                   x="Expression cluster", y="% of cluster genes") +
              theme_classic(base_size=10) +
              theme(plot.title=element_text(face="bold",size=12),
                    plot.subtitle=element_text(size=8.5,colour="#555"),
                    axis.text.y=element_text(size=8.5),
                    legend.position="right")
            save_plot(file.path(OUT_DIR,"part3_cluster_fate",
                                "cluster_fate_converged_mechanism"),
                      p_conv,
                      w=13, h=max(5, length(unique(conv_only_df$cl_label_short))*0.5+2))
          }

          # ── 3j: Per-cluster combined-evidence test (Fisher's method) ──────
          # Rather than requiring one single cluster x fate cell to survive
          # correction across the whole test family, combine the evidence
          # across ALL collapsed fate categories for each cluster into one
          # omnibus p-value: X^2 = -2*sum(log(p_i)), df=2k (Fisher's method
          # for combining independent p-values). Categories with zero
          # overlap contribute p=1 (no evidence), so every cluster is
          # combined over the same fate universe and results are comparable
          # across clusters. This surfaces clusters whose chromatin-fate
          # profile is broadly, consistently shifted even when no single
          # fate individually clears a strict per-cell threshold -- the
          # pattern the hierarchical heatmap (3h) suggests visually but
          # that 3a-3e, by construction, can never confirm on their own.
          # Wrapped in a function so 3m can run the identical test again on
          # the nearby-gene-expanded fate lists for a direct comparison.
          run_combined_test <- function(fish_sig_x, gene_lists_collapsed_x,
                                        out_stem, title_suffix="") {
            universe_x <- names(gene_lists_collapsed_x)
            rows <- list()
            for (cl_id in cluster_ids) {
              cl_label <- unique(clust_df$cluster_label[clust_df$cluster==cl_id])[1]
              sub <- fish_sig_x[fish_sig_x$cluster==cl_id, ]
              pvals <- setNames(rep(1, length(universe_x)), universe_x)
              if (nrow(sub) > 0) pvals[sub$fate] <- sub$pval
              k  <- length(pvals)
              X2 <- -2*sum(log(pmax(pvals,1e-300)))
              p_combined <- pchisq(X2, df=2*k, lower.tail=FALSE)
              drv_n <- min(3, nrow(sub))
              driver_str <- if (drv_n > 0) {
                drv <- sub[order(sub$pval), ][seq_len(drv_n), ]
                paste(sprintf("%s (p=%.3g, OR=%.2f, n=%d)",
                              pi_fate_label(drv$fate_short), drv$pval,
                              drv$odds_ratio, drv$n_overlap),
                      collapse="; ")
              } else "none tested"
              rows[[length(rows)+1]] <- data.frame(
                cluster=cl_id, cluster_label=cl_label,
                cl_label_short=sprintf("C%d: %s", cl_id, cl_label),
                n_fates_tested=k, X2=round(X2,2), pval_combined=p_combined,
                top_drivers=driver_str, stringsAsFactors=FALSE)
            }
            out_df <- do.call(rbind, rows) %>%
              mutate(padj_combined = p.adjust(pval_combined, "BH"),
                     sig_combined  = sig_stars(padj_combined)) %>%
              arrange(padj_combined)

            write.table(out_df,
                        file.path(OUT_DIR,"part3_cluster_fate", paste0(out_stem,".tsv")),
                        sep="\t", quote=FALSE, row.names=FALSE)

            cat(sprintf("\n[integrate]   Per-cluster combined-evidence test (Fisher's method)%s:\n",
                        title_suffix))
            print(out_df %>% dplyr::select(cl_label_short, n_fates_tested, pval_combined,
                                           padj_combined, sig_combined),
                  row.names=FALSE)

            p_plot <- ggplot(out_df,
                             aes(x=reorder(cl_label_short,
                                          -log10(pmax(padj_combined,1e-6))),
                                 y=-log10(pmax(padj_combined,1e-6)))) +
              geom_col(aes(fill=padj_combined<0.20), width=0.65) +
              geom_text(aes(label=sig_combined), vjust=-0.3, size=4) +
              geom_hline(yintercept=-log10(0.20), linetype="dashed", colour="#888") +
              scale_fill_manual(values=c("TRUE"="#c9579a","FALSE"="#cccccc"), guide="none") +
              coord_flip() +
              labs(title=paste0("Per-cluster combined fate-enrichment evidence (Fisher's method)",
                                title_suffix),
                   subtitle="Combines evidence across every collapsed fate category into one p-value per cluster  |  dashed line = padj 0.20",
                   x="Expression cluster", y="-log10(BH-padj), combined") +
              theme_classic(base_size=10) +
              theme(plot.title=element_text(face="bold",size=11),
                    plot.subtitle=element_text(size=8,colour="#555"))
            save_plot(file.path(OUT_DIR,"part3_cluster_fate", out_stem),
                      p_plot, w=10, h=max(5, length(cluster_ids)*0.4+2))
            out_df
          }

          combined_df <- run_combined_test(fish_sig_df, collapse_converged(fate_gene_lists),
                                           "cluster_fate_combined_test")

          # ── 3k: Permutation test — are cluster fingerprints real? ─────────
          # 3j asks "does any single cluster show enrichment"; this asks a
          # complementary, fully global question: taken together, do
          # expression clusters differ from each other in fate composition
          # more than random gene-to-cluster shuffling would produce? This
          # is a PERMANOVA-style test (Anderson 2001): build a genes x
          # collapsed-fate-category membership matrix, compute the observed
          # between-cluster sum of squared deviations from the global mean
          # profile, then compare that statistic to its null distribution
          # under repeated random reassignment of genes to clusters (cluster
          # sizes held fixed). It never depends on any single cell being
          # individually significant -- it directly tests the whole
          # fingerprint pattern the hierarchical heatmap (3h) surfaces.
          # Wrapped in a function for the same reason as 3j above.
          run_permanova <- function(gene_lists_collapsed_x, out_stem, title_suffix="") {
            perm_syms <- intersect(all_syms, clust_df[[sym_col]])
            G <- vapply(gene_lists_collapsed_x, function(g) perm_syms %in% g,
                       logical(length(perm_syms)))
            rownames(G) <- perm_syms
            cl_vec <- clust_df$cluster[match(perm_syms, clust_df[[sym_col]])]

            between_ss_x <- function(cl_assignment) {
              global_mean <- colMeans(G)
              groups <- split(seq_along(cl_assignment), cl_assignment)
              sum(vapply(groups, function(idx) {
                n  <- length(idx)
                cm <- colMeans(G[idx, , drop=FALSE])
                n * sum((cm - global_mean)^2)
              }, numeric(1)))
            }

            T_obs      <- between_ss_x(cl_vec)
            n_perm_run <- min(N_PERM, 5000)
            set.seed(1)
            T_perm <- replicate(n_perm_run, between_ss_x(sample(cl_vec)))
            p_perm <- (1 + sum(T_perm >= T_obs)) / (n_perm_run + 1)

            cat(sprintf(paste0("\n[integrate]   PERMANOVA-style global test%s: observed ",
                               "between-cluster SS=%.2f, null mean=%.2f, p=%.4g ",
                               "(%d permutations)\n"),
                        title_suffix, T_obs, mean(T_perm), p_perm, n_perm_run))

            perm_plot_df <- data.frame(T=T_perm)
            p_plot <- ggplot(perm_plot_df, aes(x=T)) +
              geom_histogram(bins=40, fill="#cccccc", colour="white") +
              geom_vline(xintercept=T_obs, colour="#c9579a", linewidth=1) +
              annotate("text", x=T_obs, y=Inf,
                       label=sprintf("Observed\np=%.4g", p_perm),
                       colour="#c9579a", hjust=-0.05, vjust=1.3, size=3.5) +
              labs(title=paste0("Are cluster fate fingerprints more distinct than chance?",
                                title_suffix),
                   subtitle=sprintf(paste0("Permutation null (n=%d): shuffling gene-to-",
                                           "cluster labels, cluster sizes held fixed"),
                                    n_perm_run),
                   x="Between-cluster sum of squares (fate composition space)",
                   y="Permutations") +
              theme_classic(base_size=10) +
              theme(plot.title=element_text(face="bold",size=11),
                    plot.subtitle=element_text(size=8,colour="#555"))
            save_plot(file.path(OUT_DIR,"part3_cluster_fate", out_stem),
                      p_plot, w=9, h=5.5)
            list(T_obs=T_obs, T_perm=T_perm, p_perm=p_perm)
          }

          collapsed_lists   <- collapse_converged(fate_gene_lists)
          permanova_nearest <- run_permanova(collapsed_lists, "cluster_fate_permanova")

          # ── 3l: Effect-size ranked table (hypothesis generation) ──────────
          # With only 13 clusters, per-cell power is fundamentally limited --
          # some real, moderate effects will never clear a strict FDR bar no
          # matter how the test family is scoped. This ranks cluster x fate
          # pairs from the FINE-grained table by effect size (log2 OR),
          # filtered to a minimum overlap count so the ranking isn't
          # dominated by n=1/n=2 noise (log2 OR is unstable at tiny counts).
          # This is explicitly NOT a significance claim -- padj is reported
          # alongside for transparency -- but surfaces plausible candidates
          # (e.g. specific Converged mechanisms) worth targeted follow-up
          # even where the formal test family lacks power to confirm them.
          MIN_OVERLAP_FOR_RANKING <- 3
          rank_df <- fish_df %>%
            filter(n_overlap >= MIN_OVERLAP_FOR_RANKING) %>%
            mutate(fate_display = pi_fate_label(fate_short)) %>%
            arrange(desc(log2_OR)) %>%
            dplyr::select(cl_label_short, fate_display, n_overlap, pct_overlap,
                          odds_ratio, log2_OR, pval, padj, sig)

          write.table(rank_df,
                      file.path(OUT_DIR,"part3_cluster_fate",
                                "cluster_fate_ranked_by_effect_size.tsv"),
                      sep="\t", quote=FALSE, row.names=FALSE)

          cat(sprintf(paste0("\n[integrate]   Top effect sizes (n_overlap>=%d, ",
                             "NOT gated by significance -- candidates, not findings):\n"),
                      MIN_OVERLAP_FOR_RANKING))
          print(head(rank_df, 20), row.names=FALSE)

          top_rank <- head(rank_df, 20) %>%
            mutate(label = sprintf("%s\n+ %s", cl_label_short, fate_display))

          if (nrow(top_rank) > 0) {
            p_ranked <- ggplot(top_rank,
                               aes(x=reorder(label, log2_OR), y=log2_OR,
                                   colour=padj<0.20)) +
              geom_segment(aes(xend=reorder(label,log2_OR), yend=0), linewidth=1.1) +
              geom_point(aes(size=n_overlap)) +
              geom_text(aes(label=sprintf("n=%d", n_overlap)),
                        hjust=-0.2, size=2.6, colour="#333") +
              scale_colour_manual(values=c("TRUE"="#c9579a","FALSE"="#999999"),
                                  name="padj<0.20\n(fine test)") +
              scale_size_continuous(range=c(2,7), guide="none") +
              coord_flip() +
              geom_hline(yintercept=0, colour="#555", linewidth=0.5) +
              labs(title=sprintf("Top 20 cluster-fate pairs by effect size (n_overlap>=%d)",
                                 MIN_OVERLAP_FOR_RANKING),
                   subtitle="Ranked by log2(OR) regardless of significance -- candidates for targeted follow-up, not confirmed findings",
                   x=NULL, y="log2(Odds Ratio)") +
              theme_classic(base_size=10) +
              theme(plot.title=element_text(face="bold",size=11),
                    plot.subtitle=element_text(size=8,colour="#555"))
            save_plot(file.path(OUT_DIR,"part3_cluster_fate",
                                "cluster_fate_ranked_by_effect_size"),
                      p_ranked, w=12, h=max(4, nrow(top_rank)*0.5+2.5))
          }

          # ── 3m: Nearest-gene vs nearby-expanded robustness check ──────────
          # ATAC peak-to-gene assignment isn't always the single nearest
          # gene -- ChIPseeker only ever reports one, and 01_peak_fate_classification.r's
          # own comments give a real example (a locus whose true target,
          # HLH54F, was missed entirely by nearest-gene-only search because
          # a closer, likely-irrelevant gene sat in between). Rather than
          # picking one assignment rule and reporting only its results,
          # re-run the SAME pipeline (Fisher grid -> collapsed categories ->
          # combined test -> PERMANOVA) on a "confident" expanded gene set
          # (nearest + anything within NEARBY_CONFIDENT_BP) and compare
          # directly against the nearest-only results above. Enrichments
          # that hold up under BOTH assignment rules are the ones worth
          # trusting; enrichments that only appear under one are exactly
          # the kind of assignment-rule artifact this check exists to catch
          # -- in EITHER direction (nearest-only misses can matter just as
          # much as expansion-only additions).
          if (length(nearby_variants$confident) == 0) {
            cat("[integrate]   3m skipped: no nearby-gene-expanded fate lists available ",
                "(re-run 01_peak_fate_classification.r to populate the nearby_genes column)\n")
          } else {
            collapsed_confident <- collapse_converged(nearby_variants$confident)
            fish_sig_confident  <- run_fisher_grid(collapsed_confident)

            if (is.null(fish_sig_confident)) {
              cat("[integrate]   3m skipped: nearby-expanded Fisher grid produced no rows\n")
            } else {
              nb_suffix <- sprintf(" \u2014 nearby-expanded (\u2264%dkb)",
                                   as.integer(NEARBY_CONFIDENT_BP/1000))
              combined_df_confident <- run_combined_test(
                fish_sig_confident, collapsed_confident,
                "cluster_fate_combined_test_nearby_confident", nb_suffix)
              permanova_confident <- run_permanova(
                collapsed_confident,
                "cluster_fate_permanova_nearby_confident", nb_suffix)

              # Per-(cluster, collapsed fate) robustness: full outer join on
              # the UNFILTERED tables so a pair tested under only one
              # assignment rule (a==0, hence absent, under the other) is
              # still captured as "one-sided" rather than silently dropped
              # -- that's precisely the failure mode nearest-gene-only
              # annotation can produce.
              cl_lookup <- setNames(
                sprintf("C%d: %s", cluster_ids,
                        sapply(cluster_ids, function(c)
                          unique(clust_df$cluster_label[clust_df$cluster==c])[1])),
                cluster_ids)

              cmp_df <- dplyr::full_join(
                  fish_sig_df %>%
                    dplyr::select(cluster, fate_short,
                                  padj_nearest=padj, OR_nearest=odds_ratio,
                                  n_nearest=n_overlap),
                  fish_sig_confident %>%
                    dplyr::select(cluster, fate_short,
                                  padj_expanded=padj, OR_expanded=odds_ratio,
                                  n_expanded=n_overlap),
                  by=c("cluster","fate_short")) %>%
                mutate(cl_label_short = cl_lookup[as.character(cluster)],
                       fate_display   = pi_fate_label(fate_short),
                       in_nearest     = !is.na(padj_nearest)  & padj_nearest  < 0.20,
                       in_expanded    = !is.na(padj_expanded) & padj_expanded < 0.20,
                       category = dplyr::case_when(
                         in_nearest & in_expanded  ~ "Robust (both)",
                         in_nearest & !in_expanded ~ "Nearest-only",
                         !in_nearest & in_expanded ~ "Expanded-only",
                         TRUE ~ "Neither (tested, not significant either)")) %>%
                arrange(cluster, fate_short)

              write.table(cmp_df,
                          file.path(OUT_DIR,"part3_cluster_fate",
                                    "cluster_fate_robustness_nearest_vs_expanded.tsv"),
                          sep="\t", quote=FALSE, row.names=FALSE)

              cat(sprintf(paste0("\n[integrate]   Robustness (padj<0.20, nearest-gene-only vs ",
                                 "nearby-confident-expanded \u2264%dkb):\n"),
                          as.integer(NEARBY_CONFIDENT_BP/1000)))
              print(table(cmp_df$category))

              robust_df <- cmp_df %>%
                filter(category != "Neither (tested, not significant either)") %>%
                mutate(label = sprintf("%s\n+ %s", cl_label_short, fate_display))

              if (nrow(robust_df) > 0) {
                p_robust <- ggplot(robust_df,
                                   aes(x=reorder(label, as.integer(factor(category))),
                                       y=1, fill=category)) +
                  geom_tile(colour="white", width=0.9, height=0.8) +
                  scale_fill_manual(values=c(
                    "Robust (both)" ="#2ca02c",
                    "Nearest-only"  ="#1f78b4",
                    "Expanded-only" ="#e31a1c"),
                    name="padj<0.20 under") +
                  coord_flip() +
                  labs(title="Robustness to nearest-gene vs nearby-gene-expanded assignment",
                       subtitle=sprintf(paste0("Nearest-gene-only fate lists vs nearest + genes ",
                                               "within %dkb  |  padj<0.20 threshold"),
                                        as.integer(NEARBY_CONFIDENT_BP/1000)),
                       x=NULL, y=NULL) +
                  theme_classic(base_size=10) +
                  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                        plot.title=element_text(face="bold",size=11),
                        plot.subtitle=element_text(size=8,colour="#555"))
                save_plot(file.path(OUT_DIR,"part3_cluster_fate",
                                    "cluster_fate_robustness_nearest_vs_expanded"),
                          p_robust, w=10, h=max(4, nrow(robust_df)*0.35+2))
              } else {
                cat("[integrate]   3m: no pairs significant under either assignment rule\n")
              }
            }
          }

          # ── 3n: Genotype-arm diagnostic — why more genes can mean LESS ────
          # significant BOTv enrichment, and a principled way to loosen it
          # ------------------------------------------------------------------
          # Two independent mechanisms can each produce exactly this
          # pattern (more total genes, but enrichment collapses onto one
          # genotype arm):
          #
          # (1) ASYMMETRIC GROWTH. If BOTCv fate gene sets grew more (in
          #     count, or in how many clusters they now overlap at all)
          #     than BOTv ones did upstream in plot_peak_fate, that alone
          #     gives BOTCv categories more/stronger real hits without
          #     BOTv's own signal changing at all. The table below reports
          #     fate-set sizes per genotype so this is checkable directly.
          #
          # (2) SHARED FDR BUDGET. padj here is BH-corrected ACROSS THE
          #     WHOLE fish_sig_df FAMILY -- both genotypes and every
          #     cluster together. BH's adjustment for the i-th smallest
          #     p-value is p_(i) * m / i, where m is the TOTAL test count.
          #     If BOTCv's growth added many new tested (cluster, fate)
          #     cells -- even ones that don't reach significance -- m goes
          #     up for EVERYONE, which can push BOTv's own marginal p-values
          #     over the padj<0.20 line even though those p-values never
          #     changed. This is not a bug in the correction; BOTv and
          #     BOTCv enrichment are conceptually distinct hypothesis
          #     families ("is this cluster enriched for a BOTv-biased
          #     fate" vs "...a BOTCv-biased fate"), so pooling them into
          #     one FDR budget is a real, checkable modelling choice, not
          #     an inherent requirement of the method.
          #
          # run_fisher_grid() (above) already computes padj_stratified: BH
          # correction done WITHIN each genotype arm separately, restoring
          # an independent FDR budget per arm. This section reports both,
          # side by side, so you can see exactly how much of the BOTv
          # signal loss (if any) is explained by (2) vs actually absent.
          cat(paste0(
            "\n[integrate]   3n NOTE: collapse_converged() had a bug -- it never actually\n",
            "               collapsed the 6 Converged mechanisms per genotype (a \"genes_\"\n",
            "               filename prefix mismatch meant its regex matched nothing), so\n",
            "               every significance test in this script has been running on the\n",
            "               full ~20-category fragmented taxonomy, not the intended ~10-\n",
            "               category collapsed one, since that feature was first added. Now\n",
            "               fixed -- re-running should change the enrichment picture on its\n",
            "               own, separately from the stratified-BH options below.\n"))
          cat("\n[integrate]   3n: fate gene-set sizes by genotype arm (fine taxonomy):\n")
          geno_size_summary <- fish_df %>%
            distinct(fate_short, n_fate) %>%
            mutate(geno_stratum = fate_genotype_bias(fate_short)) %>%
            group_by(geno_stratum) %>%
            summarise(n_categories = n(),
                      total_genes_across_categories = sum(n_fate),
                      median_category_size = median(n_fate),
                      max_category_size = max(n_fate),
                      .groups="drop") %>%
            arrange(desc(total_genes_across_categories))
          print(as.data.frame(geno_size_summary), row.names=FALSE)
          write.table(geno_size_summary,
                      file.path(OUT_DIR,"part3_cluster_fate",
                                "diagnostic_genotype_fate_set_sizes.tsv"),
                      sep="\t", quote=FALSE, row.names=FALSE)

          cat("\n[integrate]   3n: significant (cluster, fate) pairs by genotype arm,\n",
              "               global BH (shared budget) vs stratified BH (per-arm budget):\n")
          sig_count_summary <- fish_sig_df %>%
            group_by(geno_stratum) %>%
            summarise(n_tests            = n(),
                      n_sig_global       = sum(padj < 0.20),
                      n_sig_stratified   = sum(padj_stratified < 0.20),
                      recovered_by_stratifying = sum(padj >= 0.20 & padj_stratified < 0.20),
                      .groups="drop") %>%
            arrange(desc(n_tests))
          print(as.data.frame(sig_count_summary), row.names=FALSE)
          write.table(sig_count_summary,
                      file.path(OUT_DIR,"part3_cluster_fate",
                                "diagnostic_genotype_significance_budget.tsv"),
                      sep="\t", quote=FALSE, row.names=FALSE)

          if (any(sig_count_summary$recovered_by_stratifying > 0)) {
            cat(sprintf(paste0("\n[integrate]   %d (cluster, fate) pair(s) are significant under ",
                               "per-genotype stratified correction but NOT under the shared global ",
                               "budget -- consistent with the growth of one genotype's fate sets ",
                               "crowding out the other's FDR budget rather than the other genotype's ",
                               "signal genuinely disappearing. See padj_stratified / sig_stratified in\n",
                               "               %s (already written above) for the recovered pairs.\n"),
                        sum(sig_count_summary$recovered_by_stratifying),
                        "cluster_fate_enrichment_table.tsv"))
          }

          # ── 3n continued: per-FATE-CATEGORY stratification ────────────────
          # Finer than the genotype-arm view above: correct within each
          # individual fate category's own cluster tests, so growing one
          # category's gene set can NEVER change a different category's
          # padj (see padj_fate_stratified's definition comment in
          # run_fisher_grid for the full reasoning and the power/FDR-rigor
          # trade-off that comes with it).
          cat("\n[integrate]   3n: fate gene-set sizes by INDIVIDUAL category (fine taxonomy):\n")
          fate_size_summary <- fish_df %>%
            distinct(fate_short, n_fate) %>%
            mutate(geno_stratum = fate_genotype_bias(fate_short)) %>%
            arrange(desc(n_fate))
          print(as.data.frame(fate_size_summary), row.names=FALSE)
          write.table(fate_size_summary,
                      file.path(OUT_DIR,"part3_cluster_fate",
                                "diagnostic_fate_category_sizes.tsv"),
                      sep="\t", quote=FALSE, row.names=FALSE)

          cat("\n[integrate]   3n: significant (cluster, fate) pairs, global vs per-genotype vs\n",
              "               per-fate-category stratified BH (increasingly decoupled budgets):\n")
          fate_sig_summary <- fish_sig_df %>%
            group_by(fate_short) %>%
            summarise(geno_stratum       = dplyr::first(geno_stratum),
                      n_tests            = n(),
                      n_sig_global       = sum(padj < 0.20),
                      n_sig_geno_strat   = sum(padj_stratified < 0.20),
                      n_sig_fate_strat   = sum(padj_fate_stratified < 0.20),
                      recovered_by_fate_stratifying =
                        sum(padj >= 0.20 & padj_fate_stratified < 0.20),
                      .groups="drop") %>%
            arrange(desc(recovered_by_fate_stratifying))
          print(as.data.frame(fate_sig_summary), row.names=FALSE)
          write.table(fate_sig_summary,
                      file.path(OUT_DIR,"part3_cluster_fate",
                                "diagnostic_fate_significance_budget.tsv"),
                      sep="\t", quote=FALSE, row.names=FALSE)

          n_recovered_fate <- sum(fate_sig_summary$recovered_by_fate_stratifying)
          if (n_recovered_fate > 0) {
            cat(sprintf(paste0("\n[integrate]   %d (cluster, fate) pair(s) recovered by per-fate ",
                               "stratification that weren't recovered by per-genotype stratification ",
                               "either -- these are candidates whose significance was being diluted by ",
                               "OTHER fate categories' growth specifically (not just their own genotype ",
                               "arm's). Treat these as leads for follow-up, not confirmed hits: with as ",
                               "few as %d tests in some of these per-fate families, BH has very little ",
                               "left to correct for.\n"),
                        n_recovered_fate, max(fate_sig_summary$n_tests, na.rm=TRUE)))
          }
        }
      }
    }
  }
}

cat(sprintf("\n[integrate] Done. All outputs in %s/\n", OUT_DIR))
