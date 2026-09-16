#!/usr/bin/env Rscript
# =============================================================================
# atac_rna_integration_v2.R
#
# Extended ATAC × RNA integration, adding:
#
#   PART 0 — Symbol fix for Part 1 concordance
#             Patches the Entrez→symbol mapping used in genomic overlap so
#             symbols match the RNA results_full.tsv format exactly.
#
#   PART 4 — CRM / Named-enhancer integration
#             Reads crm_dar_matrix.txt files from ./Split2_Vent/results_*/
#             Links named enhancers to the RNA data and asks:
#               (a) Do genes near named enhancers show concordant expression?
#               (b) Which named enhancers are both chromatinically different
#                   AND near differentially expressed genes?
#             Outputs a ranked enhancer × RNA concordance table and dotplot.
#
#   PART 5 — Bistability concordance
#             Reads bistable peaks (sign(lfc_botv) ≠ sign(lfc_botcv), both
#             |lfc| > threshold) from the temporal DARs, links them to genes
#             via the same TSS-window overlap used in Part 1, then asks whether
#             the expression of those genes is concordant or discordant with
#             the chromatin bistability.
#             A peak is chromatinally bistable if it opens in BOTv-late but
#             closes in BOTCv-late (or vice versa). Expression concordance
#             means the nearby gene also differs between genotypes in the same
#             direction. This is the "bistable enhancer → bistable gene" test.
#             Outputs:
#               bistable_peaks_rna.tsv           full table
#               bistable_concordance_scatter.pdf  chromatin vs expression LFC
#               bistable_concordance_bars.pdf     summary by class
#               bistable_top_genes.pdf            ranked concordant hits
#
# PREREQUISITES
#   Run atac_rna_integration_v1.R first (Parts 1-3).
#   The helpers defined there (read_rna, annotate_peaks_to_genes, etc.)
#   are re-declared here so this script can run standalone.
#
# USAGE
#   Rscript atac_rna_integration_v2.R
#   Rscript atac_rna_integration_v2.R --parts 4,5
#   Rscript atac_rna_integration_v2.R --bistable-thresh 0.4
#   Rscript atac_rna_integration_v2.R --tss-window 10000
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

DAR_DIR      <- "./ATAC_DARS"
SPLIT2_DIR   <- "./Split2_Vent"      # where results_{key}/all/02_crm/ live
CRM_BED      <- "./all_drosophila_melanogaster_crms_filt.bed"
DEG_DIR     <- "../../seRNA-seq_03_19_26/celseq2_pipeline/results/combined_figures/qc_and_deg_botv_botcv/deg_limma"
CLUSTER_DIR <- "../../seRNA-seq_03_19_26/celseq2_pipeline/results/combined_figures/explore_clusters_botv_botcv"
OUT_DIR      <- "./Overview_Plots/atac_rna_integration"

# Bistable thresholds (keep in sync with plot_bistable_lollipop_enhanced.r)
BISTABLE_THRESH <- 0.25    # |lfc| floor for a bistable peak
STAB_THRESH     <- 0.1    # |lfc| ceiling for a "stable" peak
MIN_OVERLAP     <- 1L     # min bp overlap for peak matching
TSS_WINDOW      <- 5000   # bp on each side of TSS for peak→gene linkage
FDR_RNA         <- 0.10
TOP_N           <- 30     # top gene labels on plots

# ATAC temporal stems (same as bistable scripts)
STEM_BOTV_TEMP  <- "BOTv_temporal"
STEM_BOTCV_TEMP <- "BOTCv_temporal"
STEM_NC14B      <- "BOTv_vs_BOTCv"         # negate for display
STEM_NC14LATE   <- "BOTv_vs_BOTCv_nc14late" # negate for display

# RNA contrasts to use for bistability concordance (temporal ones most relevant)
BISTABLE_RNA_CONTRASTS <- list(
  list(name="BOTCv_late_vs_BOTCv_nc14b",  label="BOTCv temporal",   positive_is="BOTCv_late"),
  list(name="BOTv_nc14d_vs_BOTv_nc14b",   label="BOTv temporal",    positive_is="BOTv_late"),
  list(name="BOTCv_late_vs_BOTv_nc14b",   label="BOTCv late vs BOTv early", positive_is="BOTCv")
)

# Colours
COL_BOTV  <- "#2ca02c"; COL_BOTCV <- "#e377c2"
COL_OPEN  <- "#e31a1c"; COL_CLOSE <- "#1f78b4"
COL_CONC  <- "#2ca02c"; COL_DISC  <- "#aaaaaa"

# Temporal palette — light=early, dark=late within each genotype
COL_BOTCV_EARLY <- "#e377c2"; COL_BOTCV_LATE <- "#9e1f8e"
COL_BOTV_EARLY  <- "#2ca02c"; COL_BOTV_LATE  <- "#1a5c1a"
COL_CROSS       <- "#7b4fb7"

# Derive CRM concordant colour from ATAC direction + DAR key
crm_conc_col <- function(atac_lfc, dar_key) {
  if (grepl("BOTCv", dar_key)) {
    ifelse(atac_lfc > 0, COL_BOTCV_LATE, COL_BOTCV_EARLY)
  } else if (grepl("BOTv", dar_key)) {
    ifelse(atac_lfc > 0, COL_BOTV_LATE, COL_BOTV_EARLY)
  } else {
    ifelse(atac_lfc > 0, COL_BOTCV_LATE, COL_BOTV_LATE)
  }
}

# =============================================================================
# CLI
# =============================================================================

args <- commandArgs(trailingOnly=TRUE)
parse_arg <- function(flag, default) {
  i <- which(args==flag); if (length(i) && length(args)>i) args[i+1] else default
}
run_parts       <- as.integer(strsplit(parse_arg("--parts","4,5"),",")[[1]])
BISTABLE_THRESH <- as.numeric(parse_arg("--bistable-thresh", BISTABLE_THRESH))
TSS_WINDOW      <- as.integer(parse_arg("--tss-window",      TSS_WINDOW))

cat(sprintf("[v2] Parts: %s  |  bistable_thresh=%.1f  |  tss_window=%d\n",
            paste(run_parts,collapse=","), BISTABLE_THRESH, TSS_WINDOW))

for (d in c("part4_crm","part5_bistable"))
  dir.create(file.path(OUT_DIR,d), recursive=TRUE, showWarnings=FALSE)

# =============================================================================
# SHARED HELPERS
# =============================================================================

save_plot <- function(stem, p, w=10, h=7, dpi=200) {
  ggsave(paste0(stem,".pdf"), p, width=w, height=h, device="pdf")
  ggsave(paste0(stem,".png"), p, width=w, height=h, dpi=dpi)
  cat(sprintf("[v2]   Saved %s\n", basename(stem)))
}

sig_stars <- function(p) {
  dplyr::case_when(p<0.001~"***",p<0.01~"**",p<0.05~"*",p<0.10~".",TRUE~"")
}

read_rna <- function(contrast_name) {
  path <- file.path(DEG_DIR, contrast_name, "results_full.tsv")
  if (!file.exists(path)) { message("[v2] Not found: ",path); return(NULL) }
  df <- tryCatch(read.table(path,header=TRUE,sep="\t",quote="",stringsAsFactors=FALSE),
                 error=function(e) NULL)
  if (is.null(df)) return(NULL)
  df[!is.na(df$log2FoldChange) & !is.na(df$t), ]
}

safe_join_annotated_col <- function(bed_chr, bed_start, bed_end, bed_name,
                                    ann, col_name, stem="") {
  n_bed <- length(bed_name)
  out <- rep(NA_real_, n_bed)
  if (is.null(ann) || is.null(col_name) || is.na(col_name)) return(out)
  matched_by <- "none"
  name_col <- intersect(c("name","peak_name","peakID","peak_id","peak","ID","Row.names"),
                        colnames(ann))[1]
  if (!is.na(name_col)) {
    idx <- match(bed_name, as.character(ann[[name_col]]))
    hit <- !is.na(idx)
    if (sum(hit) >= 0.5*n_bed) {
      out[hit] <- as.numeric(ann[[col_name]])[idx[hit]]
      matched_by <- paste0("name ('", name_col, "')")
    } else {
      norm <- function(x) gsub("_DAR(?=_|$)", "", x, perl=TRUE)
      idx_n <- match(norm(bed_name), norm(as.character(ann[[name_col]])))
      hit_n <- !is.na(idx_n)
      if (sum(hit_n) >= 0.5*n_bed) {
        out[hit_n] <- as.numeric(ann[[col_name]])[idx_n[hit_n]]
        matched_by <- paste0("normalized name ('", name_col, "', _DAR token stripped)")
      }
    }
  }
  if (matched_by == "none") {
    chr_col   <- intersect(c("chr","chrom","seqnames","Chr","Chrom"), colnames(ann))[1]
    start_col <- intersect(c("start","Start","chromStart"), colnames(ann))[1]
    end_col   <- intersect(c("end","End","chromEnd"), colnames(ann))[1]
    if (!is.na(chr_col) && !is.na(start_col) && !is.na(end_col)) {
      bed_key <- paste(bed_chr, bed_start, bed_end, sep=":")
      ann_key <- paste(ann[[chr_col]], ann[[start_col]], ann[[end_col]], sep=":")
      idx <- match(bed_key, ann_key)
      hit <- !is.na(idx)
      if (sum(hit) >= 0.5*n_bed) {
        out[hit] <- as.numeric(ann[[col_name]])[idx[hit]]
        matched_by <- "genomic coordinates"
      }
    }
  }
  if (matched_by == "none" && nrow(ann) == n_bed) {
    out <- as.numeric(ann[[col_name]])
    matched_by <- "row position (row counts matched exactly)"
  }
  n_matched <- sum(!is.na(out))
  if (matched_by == "none") {
    warning(sprintf(
      "[safe_join_annotated_col: %s] Could not safely align '%s' to %d bed peaks -- no usable name/coordinate match and row counts don't agree (%d annotated vs %d bed). Returning all-NA.",
      stem, col_name, n_bed, nrow(ann), n_bed))
  }
  cat(sprintf("  [%s] '%s' joined via %s (%d/%d matched)\n", stem, col_name, matched_by, n_matched, n_bed))
  out
}

read_bed_with_lfc <- function(stem, negate=FALSE) {
  # Reads BED + annotated file to get coordinates + real log2FC
  bed_path <- file.path(DAR_DIR, paste0(stem,"_DARs.bed"))
  if (!file.exists(bed_path)) { message("[v2] BED not found: ",bed_path); return(NULL) }
  bed <- tryCatch(read.table(bed_path,sep="\t",header=FALSE,stringsAsFactors=FALSE),
                  error=function(e) NULL)
  if (is.null(bed)||nrow(bed)==0) return(NULL)
  bed_name <- if(ncol(bed)>=4) as.character(bed[[4]]) else paste0(stem,"_",seq_len(nrow(bed)))
  lfc <- rep(NA_real_,nrow(bed))
  padj<- rep(NA_real_,nrow(bed))
  for (ext in c("_DARs_annotated.txt","_DARs_annotated_peaks.csv")) {
    ap <- file.path(DAR_DIR,paste0(stem,ext))
    if (!file.exists(ap)) next
    sep_ch <- if (grepl("\t",readLines(ap,1,warn=FALSE),fixed=TRUE)) "\t" else ","
    ann <- tryCatch(read.table(ap,header=TRUE,sep=sep_ch,stringsAsFactors=FALSE,
                               quote="",fill=TRUE,comment.char="",row.names=NULL),
                    error=function(e) NULL)
    if (is.null(ann)) next
    colnames(ann)[1] <- gsub("^[^[:alnum:]]+","",colnames(ann)[1])
    lc <- intersect(c("log2FoldChange","log2FC","LFC"),colnames(ann))[1]
    pc <- intersect(c("padj","adj.P.Val","FDR","p_adj"),colnames(ann))[1]
    if (!is.na(lc)) {
      lfc <- safe_join_annotated_col(bed[[1]], bed[[2]], bed[[3]], bed_name, ann, lc, stem)
      if (!is.na(pc)) padj <- safe_join_annotated_col(bed[[1]], bed[[2]], bed[[3]], bed_name, ann, pc, stem)
      break
    }
  }
  if (all(is.na(lfc)) && ncol(bed)>=5) {
    r5 <- suppressWarnings(as.numeric(bed[,5]))
    if (any(r5<0,na.rm=TRUE)) lfc <- r5
  }
  data.frame(
    chr  = as.character(bed[[1]]), start=as.integer(bed[[2]]),
    end  = as.integer(bed[[3]]),
    name = bed_name,
    atac_lfc = lfc * ifelse(negate,-1,1),
    atac_padj = padj,
    stringsAsFactors=FALSE)
}

# Link peaks to genes via TSS-window overlap, returning symbols
# Uses TxDb + org.Dm.eg.db when available; also tries the gene BED file
annotate_peaks <- function(peak_df, tss_window=TSS_WINDOW) {
  if (is.null(peak_df)||nrow(peak_df)==0) return(peak_df)
  peak_df$symbol  <- NA_character_
  peak_df$gene_id <- NA_character_

  has_txdb <- requireNamespace("TxDb.Dmelanogaster.UCSC.dm6.ensGene",quietly=TRUE)
  has_org  <- requireNamespace("org.Dm.eg.db", quietly=TRUE)

  valid <- !is.na(peak_df$chr) & !is.na(peak_df$start) & !is.na(peak_df$end)
  if (sum(valid)==0) return(peak_df)
  pk_gr <- GRanges(seqnames=peak_df$chr[valid],
                    ranges=IRanges(start=peak_df$start[valid]+1L,
                                   end=peak_df$end[valid]))
  mcols(pk_gr)$row_idx <- which(valid)

  if (has_txdb) {
    suppressPackageStartupMessages({
      library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
    })
    txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
    genes_gr <- suppressWarnings(genes(txdb))
    tss_gr   <- resize(genes_gr, width=1L, fix="start")
    tss_gr   <- suppressWarnings(trim(resize(tss_gr, width=tss_window*2+1L, fix="center")))

    hits <- suppressWarnings(findOverlaps(pk_gr, tss_gr))
    if (length(hits)>0) {
      hit_df <- data.frame(
        row_idx = mcols(pk_gr)$row_idx[queryHits(hits)],
        entrez  = names(tss_gr)[subjectHits(hits)],
        stringsAsFactors=FALSE)
      # One gene per peak: keep smallest TSS distance
      hit_df <- hit_df[!duplicated(hit_df$row_idx), ]

      if (has_org) {
        suppressPackageStartupMessages(library(org.Dm.eg.db))
        # Try ENTREZID → SYMBOL
        sym_map <- tryCatch(
          AnnotationDbi::select(org.Dm.eg.db, keys=unique(hit_df$entrez),
                                columns="SYMBOL", keytype="ENTREZID"),
          error=function(e) NULL)
        # Also try FLYBASE keytype in case TxDb uses FBgn IDs
        if (is.null(sym_map) || nrow(sym_map)==0 ||
            all(is.na(sym_map$SYMBOL))) {
          sym_map2 <- tryCatch(
            AnnotationDbi::select(org.Dm.eg.db, keys=unique(hit_df$entrez),
                                  columns="SYMBOL", keytype="FLYBASE"),
            error=function(e) NULL)
          if (!is.null(sym_map2) && any(!is.na(sym_map2$SYMBOL)))
            sym_map <- sym_map2
        }
        if (!is.null(sym_map) && "SYMBOL" %in% colnames(sym_map)) {
          key_col <- intersect(c("ENTREZID","FLYBASE","GENEID"),colnames(sym_map))[1]
          hit_df  <- merge(hit_df, sym_map[,c(key_col,"SYMBOL")],
                           by.x="entrez", by.y=key_col, all.x=TRUE)
        } else {
          hit_df$SYMBOL <- hit_df$entrez
        }
      } else {
        hit_df$SYMBOL <- hit_df$entrez
      }

      peak_df$symbol[hit_df$row_idx]  <- hit_df$SYMBOL
      peak_df$gene_id[hit_df$row_idx] <- hit_df$entrez
      cat(sprintf("[v2]   Annotated %d / %d peaks to genes (TSS±%d bp)\n",
                  sum(!is.na(peak_df$symbol)), nrow(peak_df), tss_window))
    }
  } else {
    # Fallback: use dm6_genes.bed if present
    gene_bed <- tryCatch(read.table("dm6_genes.bed",sep="\t",header=FALSE,
                                    stringsAsFactors=FALSE), error=function(e) NULL)
    if (!is.null(gene_bed) && ncol(gene_bed)>=6) {
      g_gr <- GRanges(seqnames=gene_bed[[1]],
                      ranges=IRanges(pmax(1,gene_bed[[2]]-tss_window),
                                     gene_bed[[3]]+tss_window))
      mcols(g_gr)$symbol <- gene_bed[[4]]
      hits <- suppressWarnings(findOverlaps(pk_gr, g_gr))
      if (length(hits)>0) {
        hit_df <- data.frame(
          row_idx=mcols(pk_gr)$row_idx[queryHits(hits)],
          symbol =mcols(g_gr)$symbol[subjectHits(hits)],
          stringsAsFactors=FALSE)
        hit_df <- hit_df[!duplicated(hit_df$row_idx), ]
        peak_df$symbol[hit_df$row_idx] <- hit_df$symbol
        cat(sprintf("[v2]   Gene-bed fallback: %d peaks annotated\n",
                    sum(!is.na(peak_df$symbol))))
      }
    } else {
      cat("[v2]   Install TxDb.Dmelanogaster.UCSC.dm6.ensGene + org.Dm.eg.db\n")
      cat("[v2]   for automatic peak→gene annotation:\n")
      cat("[v2]   BiocManager::install(c('TxDb.Dmelanogaster.UCSC.dm6.ensGene',\n")
      cat("[v2]                           'org.Dm.eg.db'))\n")
    }
  }
  peak_df
}

# Merge annotated peaks with RNA by symbol (preferred) or gene_id
join_peaks_rna <- function(peak_df, rna_df) {
  if (is.null(peak_df)||is.null(rna_df)) return(NULL)

  rna_cols <- intersect(c("symbol","gene_id","log2FoldChange","t","pvalue","padj"),
                         colnames(rna_df))
  rna_sub  <- rna_df[, rna_cols, drop=FALSE]
  # rename to avoid collision
  colnames(rna_sub)[colnames(rna_sub)=="log2FoldChange"] <- "rna_lfc"
  colnames(rna_sub)[colnames(rna_sub)=="t"]              <- "rna_t"
  colnames(rna_sub)[colnames(rna_sub)=="pvalue"]         <- "rna_pval"
  colnames(rna_sub)[colnames(rna_sub)=="padj"]           <- "rna_padj"
  colnames(rna_sub)[colnames(rna_sub)=="symbol"]         <- "rna_symbol"
  colnames(rna_sub)[colnames(rna_sub)=="gene_id"]        <- "rna_gene_id"

  merged <- NULL

  # Try symbol merge
  n_pk_sym <- sum(!is.na(peak_df$symbol) & peak_df$symbol!="")
  n_rn_sym <- if ("rna_symbol" %in% names(rna_sub))
    sum(!is.na(rna_sub$rna_symbol) & rna_sub$rna_symbol!="") else 0
  if (n_pk_sym>0 && n_rn_sym>0) {
    m1 <- merge(peak_df[!is.na(peak_df$symbol) & peak_df$symbol!="", ],
                rna_sub[!is.na(rna_sub$rna_symbol) & rna_sub$rna_symbol!="", ],
                by.x="symbol", by.y="rna_symbol", all=FALSE)
    if (nrow(m1)>0) { cat(sprintf("[v2]   Merged by symbol: %d\n",nrow(m1))); merged<-m1 }
  }

  # Fallback: gene_id
  if ((is.null(merged)||nrow(merged)==0) &&
      "gene_id" %in% colnames(peak_df) &&
      "rna_gene_id" %in% names(rna_sub)) {
    m2 <- merge(peak_df[!is.na(peak_df$gene_id)&peak_df$gene_id!="", ],
                rna_sub[!is.na(rna_sub$rna_gene_id)&rna_sub$rna_gene_id!="", ],
                by.x="gene_id", by.y="rna_gene_id", all=FALSE)
    if (nrow(m2)>0) {
      if (!"symbol" %in% colnames(m2)) m2$symbol <- m2$gene_id
      cat(sprintf("[v2]   Merged by gene_id: %d\n",nrow(m2))); merged<-m2
    }
  }

  if (is.null(merged)||nrow(merged)==0) {
    cat("[v2]   WARNING: 0 genes merged\n")
    # Show first 5 of each for diagnosis
    pk_s <- head(peak_df$symbol[!is.na(peak_df$symbol)],5)
    rn_s <- if ("rna_symbol" %in% names(rna_sub))
      head(rna_sub$rna_symbol[!is.na(rna_sub$rna_symbol)],5) else character(0)
    cat(sprintf("[v2]     Peak symbols: %s\n",paste(pk_s,collapse=", ")))
    cat(sprintf("[v2]     RNA  symbols: %s\n",paste(rn_s,collapse=", ")))
    return(NULL)
  }
  merged <- merged[!is.na(merged$atac_lfc) & !is.na(merged$rna_lfc), ]
  merged
}

# =============================================================================
# PART 4 — CRM / NAMED ENHANCER INTEGRATION
# =============================================================================

if (4 %in% run_parts) {
  cat("\n", strrep("=",70), "\n", sep="")
  cat("PART 4: CRM / NAMED ENHANCER INTEGRATION\n\n")

  # ── 4a: Read crm_dar_matrix files ──────────────────────────────────────────
  # Rows = named CRMs, columns include open/close peak counts and LFC summaries
  crm_keys <- c("BOTv_vs_BOTCv","BOTv_vs_BOTCv_late","BOTCv_temporal","BOTv_temporal")
  crm_data  <- list()

  for (key in crm_keys) {
    candidates <- c(
      file.path(SPLIT2_DIR, paste0("results_",key), "all","02_crm","crm_dar_matrix.txt"),
      file.path(SPLIT2_DIR, paste0("results_",key), "crm_dar_matrix.txt"),
      file.path(SPLIT2_DIR, paste0("results_",key), "crm","crm_dar_matrix.txt")
    )
    mat_path <- Filter(file.exists, candidates)[1]
    if (is.na(mat_path)) {
      cat(sprintf("[v2]   CRM matrix not found for %s (checked %d paths)\n",
                  key, length(candidates)))
      next
    }
    m <- tryCatch(read.table(mat_path, header=TRUE, sep="\t", quote="",
                              stringsAsFactors=FALSE, row.names=1,
                              check.names=FALSE),
                  error=function(e) { message("[v2] Error reading ",mat_path,": ",e$message); NULL })
    if (is.null(m)||nrow(m)==0) next
    cat(sprintf("[v2]   CRM matrix %s: %d CRMs × %d cols  (%s)\n",
                key, nrow(m), ncol(m), basename(mat_path)))
    crm_data[[key]] <- as.data.frame(m)
    crm_data[[key]]$crm_name <- rownames(m)
    crm_data[[key]]$dar_key  <- key
  }

  if (length(crm_data)==0) {
    cat("[v2]   No CRM matrices found. Check SPLIT2_DIR =", SPLIT2_DIR, "\n")
    cat("[v2]   Expected: Split2_Vent/results_{key}/all/02_crm/crm_dar_matrix.txt\n")
    cat("[v2]   Run MASTER_RUN_FILE_Split2_Vent.r with run_crm=TRUE first.\n")
  } else {
    # ── 4b: Load CRM BED to get genomic coordinates ─────────────────────────
    crm_gr <- NULL
    if (file.exists(CRM_BED)) {
      crm_bed <- tryCatch(read.table(CRM_BED,sep="\t",header=FALSE,
                                     stringsAsFactors=FALSE), error=function(e) NULL)
      if (!is.null(crm_bed) && ncol(crm_bed)>=4) {
        crm_gr <- GRanges(seqnames=crm_bed[[1]],
                           ranges=IRanges(as.integer(crm_bed[[2]])+1L,
                                          as.integer(crm_bed[[3]])),
                           name=as.character(crm_bed[[4]]))
        cat(sprintf("[v2]   CRM BED: %d named enhancers\n", length(crm_gr)))
      }
    } else {
      cat(sprintf("[v2]   CRM BED not found at %s\n", CRM_BED))
    }

    # ── 4c: For each RNA contrast, join CRM genes via gene column in matrix ──
    # The crm_dar_matrix usually has a column like "gene", "gene_symbol",
    # "nearest_gene" listing the gene each CRM is associated with in the
    # literature/database. If present, use it directly; otherwise fall back
    # to the CRM name (many named CRMs are named after their gene).

    for (rna_cfg in BISTABLE_RNA_CONTRASTS) {
      rna_df <- read_rna(rna_cfg$name)
      if (is.null(rna_df)) next

      for (key in names(crm_data)) {
        mat <- crm_data[[key]]

        # Find gene column
        gene_col_crm <- intersect(c("gene","gene_symbol","nearest_gene",
                                     "Gene","SYMBOL","associated_gene"),
                                   colnames(mat))[1]

        # If no gene column, parse gene name from CRM name
        # (most CRM names are like "eve_IAB2", "hkb_NE" — gene is prefix)
        if (is.na(gene_col_crm)) {
          mat$gene_from_name <- sub("[_-].*","", mat$crm_name)
          gene_col_crm <- "gene_from_name"
          cat(sprintf("[v2]   %s: no gene col — inferring from CRM name prefix\n", key))
        }

        # Find LFC column in CRM matrix
        lfc_col_crm <- intersect(c("mean_lfc","mean_log2FC","lfc","log2FC",
                                    "median_lfc","effect_size"), colnames(mat))[1]
        if (is.na(lfc_col_crm) && any(grepl("lfc|log2",colnames(mat),ignore.case=TRUE)))
          lfc_col_crm <- grep("lfc|log2",colnames(mat),ignore.case=TRUE,value=TRUE)[1]

        # Find n_open / n_close columns
        n_open_col  <- intersect(c("n_open","n_opening","n_gained","n_up"), colnames(mat))[1]
        n_close_col <- intersect(c("n_close","n_closing","n_lost","n_down"), colnames(mat))[1]

        cat(sprintf("[v2]   CRM %s × RNA %s: gene_col=%s, lfc_col=%s\n",
                    key, rna_cfg$name,
                    ifelse(is.na(gene_col_crm),"[none]",gene_col_crm),
                    ifelse(is.na(lfc_col_crm),"[none]",lfc_col_crm)))

        # Build CRM-gene link table
        crm_gene_df <- data.frame(
          crm_name  = mat$crm_name,
          crm_gene  = if (!is.na(gene_col_crm)) as.character(mat[[gene_col_crm]])
                      else mat$crm_name,
          crm_atac_lfc = if (!is.na(lfc_col_crm)) as.numeric(mat[[lfc_col_crm]])
                         else NA_real_,
          n_open    = if (!is.na(n_open_col))  as.integer(mat[[n_open_col]])  else NA_integer_,
          n_close   = if (!is.na(n_close_col)) as.integer(mat[[n_close_col]]) else NA_integer_,
          stringsAsFactors=FALSE)

        # ── Sign convention for concordance ───────────────────────────────────
        # The crm_atac_lfc sign depends on which DAR key produced it:
        #
        # Cross-genotype keys (BOTv_vs_BOTCv, BOTv_vs_BOTCv_late):
        #   v5 positive = BOTCv more open
        #   RNA positive = higher in BOTCv/late arm (A)
        #   → concordant when SAME sign
        #
        # Temporal keys (BOTCv_temporal, BOTv_temporal):
        #   positive = nc14b more open = LOSING accessibility over time
        #   RNA positive = higher in late arm (A)
        #   → concordant when opening chromatin (negative ATAC) AND up RNA
        #     i.e. concordant when OPPOSITE signs
        #   → flip ATAC sign before concordance test so same-sign = concordant
        #
        atac_sign_flip <- grepl("temporal", key, ignore.case=TRUE)
        crm_gene_df$crm_atac_lfc_for_conc <- crm_gene_df$crm_atac_lfc *
          ifelse(atac_sign_flip, -1, 1)
        if (atac_sign_flip)
          cat(sprintf("[v2]   %s: temporal key — flipping ATAC sign for concordance\n",key))

        # Merge with RNA
        rna_sub <- rna_df[, c("symbol","log2FoldChange","t","padj"), drop=FALSE]
        colnames(rna_sub) <- c("symbol","rna_lfc","rna_t","rna_padj")

        merged_crm <- merge(crm_gene_df,
                             rna_sub[!is.na(rna_sub$symbol)&rna_sub$symbol!="", ],
                             by.x="crm_gene", by.y="symbol", all=FALSE)
        if (nrow(merged_crm)==0) {
          cat(sprintf("[v2]   No CRM-RNA overlap for %s × %s\n",key,rna_cfg$name))
          next
        }
        cat(sprintf("[v2]   CRM-RNA merged: %d CRMs with expression data\n",
                    nrow(merged_crm)))

        # Concordant = same sign on the sign-corrected ATAC lFC and the RNA lFC
        merged_crm$concordant <- !is.na(merged_crm$crm_atac_lfc_for_conc) &
          !is.na(merged_crm$rna_lfc) &
          sign(merged_crm$crm_atac_lfc_for_conc) == sign(merged_crm$rna_lfc)
        merged_crm$rna_sig     <- !is.na(merged_crm$rna_padj) &
          merged_crm$rna_padj < FDR_RNA
        merged_crm$combined_score <- abs(merged_crm$crm_atac_lfc) *
          (-log10(pmax(merged_crm$rna_padj, 1e-10)))

        # Write full table
        write.table(merged_crm,
                    file.path(OUT_DIR,"part4_crm",
                              sprintf("crm_rna_%s_%s.tsv",key,rna_cfg$name)),
                    sep="\t", quote=FALSE, row.names=FALSE)

        # Plot: CRM ATAC lfc vs RNA lfc — labelled by CRM name
        if (is.na(lfc_col_crm)) next  # can't scatter without ATAC lfc

        to_label <- merged_crm %>%
          arrange(desc(combined_score)) %>%
          head(TOP_N)

        # Concordant colour reflects ATAC direction + genotype identity of the DAR key
        merged_crm <- merged_crm %>%
          mutate(conc_display_col = dplyr::case_when(
            !concordant ~ COL_DISC,
            crm_atac_lfc > 0 & grepl("BOTCv", key) ~ COL_BOTCV_LATE,
            crm_atac_lfc < 0 & grepl("BOTCv", key) ~ COL_BOTCV_EARLY,
            crm_atac_lfc > 0 & grepl("BOTv",  key) ~ COL_BOTV_LATE,
            crm_atac_lfc < 0 & grepl("BOTv",  key) ~ COL_BOTV_EARLY,
            crm_atac_lfc > 0 ~ COL_BOTCV_LATE,
            TRUE ~ COL_BOTV_LATE),
            conc_label = dplyr::case_when(
              !concordant ~ "Discordant",
              crm_atac_lfc > 0 ~ sprintf("Concordant: %s open",
                ifelse(grepl("BOTCv",key),"BOTCv/late","BOTCv")),
              TRUE ~ sprintf("Concordant: %s open",
                ifelse(grepl("BOTv",key),"BOTv/late","BOTv"))))

        conc_col_scale <- setNames(
          unique(merged_crm$conc_display_col),
          unique(merged_crm$conc_label))

        p_crm_scatter <- ggplot(merged_crm,
                                 aes(x=crm_atac_lfc_for_conc, y=rna_lfc,
                                     colour=conc_label, alpha=rna_sig)) +
          geom_hline(yintercept=0,lty=2,colour="#ccc") +
          geom_vline(xintercept=0,lty=2,colour="#ccc") +
          geom_point(size=1.8) +
          scale_colour_manual(values=conc_col_scale,
                              name="ATAC\u2194RNA") +
          scale_alpha_manual(values=c("TRUE"=0.9,"FALSE"=0.35),
                             labels=c("TRUE"=sprintf("RNA padj<%.2f",FDR_RNA),
                                      "FALSE"="Not sig"),
                             name="RNA") +
          {
            if (has_ggrepel && nrow(to_label)>0)
              geom_text_repel(data=to_label,
                              aes(x=crm_atac_lfc_for_conc, y=rna_lfc, label=crm_name),
                              colour="black", size=2.5, fontface="italic",
                              max.overlaps=20, box.padding=0.4,
                              segment.size=0.3, seed=42, inherit.aes=FALSE)
          } +
          labs(title=sprintf("Named enhancer × Expression concordance\n%s  |  RNA: %s",
                             key, rna_cfg$label),
               subtitle=sprintf("n=%d CRMs  |  %d concordant  |  %d RNA sig\n%s",
                                nrow(merged_crm),
                                sum(merged_crm$concordant, na.rm=TRUE),
                                sum(merged_crm$rna_sig, na.rm=TRUE),
                                ifelse(atac_sign_flip,
                                       "ATAC x-axis flipped: positive = more open in late/A arm",
                                       "ATAC x-axis: positive = BOTCv/late arm more open")),
               x=sprintf("CRM ATAC log2FC%s",
                          ifelse(atac_sign_flip," (sign-flipped: + = opening over time)","")),
               y=sprintf("RNA log2FC (%s)", rna_cfg$positive_is)) +
          theme_classic(base_size=11) +
          theme(plot.title=element_text(face="bold",size=10.5))

        # Top concordant CRMs — ranked lollipop
        conc_crms <- merged_crm %>%
          filter(concordant, !is.na(crm_atac_lfc)) %>%
          arrange(desc(combined_score)) %>%
          head(TOP_N) %>%
          mutate(label=sprintf("%s (%s)", crm_name, crm_gene))

        if (nrow(conc_crms)>0) {
          p_top_crm <- ggplot(conc_crms,
                               aes(x=reorder(label,combined_score),
                                   y=combined_score, fill=crm_atac_lfc>0)) +
            geom_col(width=0.7) +
            geom_text(aes(label=sprintf("lFC=%.2f / rna=%.2f",
                                        crm_atac_lfc, rna_lfc)),
                      hjust=-0.05, size=2.8) +
            scale_fill_manual(values=c("TRUE"=COL_BOTCV,"FALSE"=COL_BOTV),
                              labels=c("TRUE"="BOTCv more open",
                                       "FALSE"="BOTv more open"),
                              name="ATAC direction") +
            scale_y_continuous(expand=expansion(mult=c(0.02,0.35))) +
            coord_flip(clip="off") +
            labs(title="Top concordant CRMs",
                 subtitle="|ATAC lFC| × −log10(RNA padj)",
                 x=NULL, y="Combined score") +
            theme_classic(base_size=10) +
            theme(plot.title=element_text(face="bold"),
                  plot.margin=margin(5,120,5,5,"pt"))

          p_combined <- p_crm_scatter | p_top_crm
        } else {
          p_combined <- p_crm_scatter
        }

        stem_out <- file.path(OUT_DIR,"part4_crm",
                              sprintf("crm_concordance_%s_%s",key,rna_cfg$name))
        save_plot(stem_out, p_combined, w=14, h=6)
      }
    }

    # ── 4d: Summary across all CRM sets — which named enhancers appear ───────
    # consistently concordant across multiple contrasts?
    all_crm_rows <- list()
    for (key in names(crm_data)) {
      mat <- crm_data[[key]]
      gene_col_crm <- intersect(c("gene","gene_symbol","nearest_gene",
                                   "Gene","SYMBOL","associated_gene"),
                                 colnames(mat))[1]
      lfc_col_crm  <- intersect(c("mean_lfc","mean_log2FC","lfc","log2FC"),
                                 colnames(mat))[1]
      if (is.na(lfc_col_crm)) next
      all_crm_rows[[key]] <- data.frame(
        crm_name = mat$crm_name,
        dar_key  = key,
        atac_lfc = as.numeric(mat[[lfc_col_crm]]),
        stringsAsFactors=FALSE)
    }
    if (length(all_crm_rows)>0) {
      crm_summary <- do.call(rbind, all_crm_rows)
      # Which CRMs appear in >1 DAR set
      multi_crms <- crm_summary %>%
        group_by(crm_name) %>%
        summarise(n_dar_sets=n(),
                  mean_atac_lfc=mean(atac_lfc, na.rm=TRUE),
                  consistent_direction=all(sign(atac_lfc[!is.na(atac_lfc)])==
                                             sign(atac_lfc[!is.na(atac_lfc)][1])),
                  .groups="drop") %>%
        filter(n_dar_sets>1) %>%
        arrange(desc(abs(mean_atac_lfc)))

      write.table(multi_crms,
                  file.path(OUT_DIR,"part4_crm","crm_multi_dar_summary.tsv"),
                  sep="\t", quote=FALSE, row.names=FALSE)
      cat(sprintf("[v2]   %d CRMs appear in >1 DAR set; %d consistently directional\n",
                  nrow(multi_crms), sum(multi_crms$consistent_direction, na.rm=TRUE)))
    }
  }
}

# =============================================================================
# PART 5 — BISTABILITY CONCORDANCE
# =============================================================================

if (5 %in% run_parts) {
  cat("\n", strrep("=",70), "\n", sep="")
  cat("PART 5: BISTABILITY CONCORDANCE\n\n")

  # ── 5a: Load temporal DARs and find bistable peaks ───────────────────────
  cat("[v2]   Loading temporal DAR sets...\n")
  df_botv  <- read_bed_with_lfc(STEM_BOTV_TEMP,  negate=FALSE)
  df_botcv <- read_bed_with_lfc(STEM_BOTCV_TEMP, negate=FALSE)

  if (is.null(df_botv)||is.null(df_botcv)) {
    cat("[v2]   Cannot load temporal DARs — skipping Part 5\n")
  } else {
    cat(sprintf("[v2]   BOTv temporal: %d peaks\n", nrow(df_botv)))
    cat(sprintf("[v2]   BOTCv temporal: %d peaks\n", nrow(df_botcv)))

    # Build GRanges and find overlapping peaks
    gr_botv  <- GRanges(seqnames=df_botv$chr,
                         ranges=IRanges(df_botv$start+1L, df_botv$end),
                         name=df_botv$name, atac_lfc=df_botv$atac_lfc)
    gr_botcv <- GRanges(seqnames=df_botcv$chr,
                         ranges=IRanges(df_botcv$start+1L, df_botcv$end),
                         name=df_botcv$name, atac_lfc=df_botcv$atac_lfc)

    hits <- findOverlaps(gr_botv, gr_botcv, minoverlap=MIN_OVERLAP)
    ov_w <- width(pintersect(gr_botv[queryHits(hits)],
                              gr_botcv[subjectHits(hits)]))
    recip <- pmin(ov_w/width(gr_botv[queryHits(hits)]),
                  ov_w/width(gr_botcv[subjectHits(hits)]))

    # Best reciprocal overlap per BOTv peak
    best_hits <- data.frame(i_botv=queryHits(hits),
                             i_botcv=subjectHits(hits), recip=recip) %>%
      group_by(i_botv) %>% slice_max(recip, n=1, with_ties=FALSE) %>% ungroup()

    shared <- data.frame(
      chr       = as.character(seqnames(gr_botv[best_hits$i_botv])),
      start     = start(gr_botv[best_hits$i_botv]) - 1L,
      end       = end(gr_botv[best_hits$i_botv]),
      name      = mcols(gr_botv[best_hits$i_botv])$name,
      lfc_botv  = mcols(gr_botv[best_hits$i_botv])$atac_lfc,
      lfc_botcv = mcols(gr_botcv[best_hits$i_botcv])$atac_lfc,
      stringsAsFactors=FALSE) %>%
      filter(!is.na(lfc_botv), !is.na(lfc_botcv))

    cat(sprintf("[v2]   Shared (overlapping) peaks: %d\n", nrow(shared)))

    # Classify bistability
    shared <- shared %>% mutate(
      is_stable   = abs(lfc_botv) < STAB_THRESH & abs(lfc_botcv) < STAB_THRESH,
      is_bistable = !is_stable &
        abs(lfc_botv) > BISTABLE_THRESH &
        abs(lfc_botcv) > BISTABLE_THRESH &
        sign(lfc_botv) != sign(lfc_botcv),
      is_concordant_open  = !is_stable &
        lfc_botv > BISTABLE_THRESH &
        lfc_botcv > BISTABLE_THRESH,
      is_concordant_close = !is_stable &
        lfc_botv < -BISTABLE_THRESH &
        lfc_botcv < -BISTABLE_THRESH,
      divergence_score    = abs(lfc_botv) + abs(lfc_botcv),
      bistable_class = case_when(
        is_bistable & lfc_botv > 0 ~ "BOTv opens / BOTCv closes",
        is_bistable & lfc_botv < 0 ~ "BOTCv opens / BOTv closes",
        is_concordant_open          ~ "Both open",
        is_concordant_close         ~ "Both close",
        is_stable                   ~ "Stable",
        TRUE                        ~ "Intermediate"))

    cat(sprintf("[v2]   Bistable peaks: %d\n", sum(shared$is_bistable)))
    cat(sprintf("[v2]     BOTv-opens/BOTCv-closes: %d\n",
                sum(shared$bistable_class=="BOTv opens / BOTCv closes")))
    cat(sprintf("[v2]     BOTCv-opens/BOTv-closes: %d\n",
                sum(shared$bistable_class=="BOTCv opens / BOTv closes")))

    # ── 5b: Annotate bistable peaks to genes ─────────────────────────────────
    cat("[v2]   Annotating peaks to genes...\n")
    shared_ann <- annotate_peaks(shared, tss_window=TSS_WINDOW)

    # ── 5c: For each RNA contrast, test bistable concordance ─────────────────
    for (rna_cfg in BISTABLE_RNA_CONTRASTS) {
      rna_df <- read_rna(rna_cfg$name)
      if (is.null(rna_df)) next

      merged <- join_peaks_rna(shared_ann, rna_df)
      if (is.null(merged)||nrow(merged)<3) {
        cat(sprintf("[v2]   Insufficient merged data for %s\n", rna_cfg$name))
        next
      }
      cat(sprintf("[v2]   Bistable × %s: %d peaks with expression\n",
                  rna_cfg$name, nrow(merged)))

      # Chromatin-expression concordance for bistable peaks
      # A bistable peak is "expression-concordant" if:
      #   the gene it's near has higher expression in the genotype where the
      #   peak is more open. Since temporal RNA = comparing late vs early
      #   within a genotype, concordance means the RNA LFC direction matches
      #   the ATAC LFC direction in that genotype.
      merged <- merged %>% mutate(
        rna_sig       = !is.na(rna_padj) & rna_padj < FDR_RNA,
        # For cross-genotype contrasts: BOTCv late > BOTv early
        # rna_lfc > 0 = higher in BOTCv/late arm
        # ATAC sign: lfc_botcv > 0 = more open in BOTCv temporal
        chromatin_vs_rna = case_when(
          bistable_class=="BOTCv opens / BOTv closes" & rna_lfc>0 ~
            "Concordant (BOTCv-biased)",
          bistable_class=="BOTCv opens / BOTv closes" & rna_lfc<0 ~
            "Discordant",
          bistable_class=="BOTv opens / BOTCv closes" & rna_lfc<0 ~
            "Concordant (BOTv-biased)",
          bistable_class=="BOTv opens / BOTCv closes" & rna_lfc>0 ~
            "Discordant",
          bistable_class %in% c("Both open","Both close") ~
            "Co-directional chromatin",
          TRUE ~ bistable_class
        ))

      # Save full table
      write.table(merged,
                  file.path(OUT_DIR,"part5_bistable",
                            sprintf("bistable_rna_%s.tsv", rna_cfg$name)),
                  sep="\t", quote=FALSE, row.names=FALSE)

      # ── 5d: Scatter: BOTv lfc vs BOTCv lfc, coloured by RNA lfc ──────────
      bistable_only <- filter(merged, is_bistable)
      all_shared    <- merged  # use all shared for scatter background

      # Continuous RNA LFC as colour on chromatin bistability axes
      p_bist_scatter <- ggplot(all_shared, aes(x=lfc_botv, y=lfc_botcv)) +
        geom_hline(yintercept=0, lty=2, colour="#ccc") +
        geom_vline(xintercept=0, lty=2, colour="#ccc") +
        # Background: all shared peaks (grey)
        geom_point(data=filter(all_shared, !is_bistable),
                   colour="#dddddd", size=0.7, alpha=0.4) +
        # Bistable peaks coloured by RNA lfc
        geom_point(data=bistable_only,
                   aes(colour=rna_lfc, size=abs(rna_lfc)),
                   alpha=0.85) +
        scale_colour_gradient2(low=COL_BOTV, mid="#f0f0f0", high=COL_BOTCV,
                               midpoint=0,
                               name=sprintf("RNA lFC\n(%s)", rna_cfg$positive_is),
                               na.value="#cccccc") +
        scale_size_continuous(range=c(1,4), guide="none") +
        # Diagonal "bistable zone" lines
        geom_abline(slope=-1, intercept=0, lty=3, colour="#999", linewidth=0.5) +
        {
          if (has_ggrepel && nrow(bistable_only)>0) {
            top_b <- bistable_only %>%
              filter(!is.na(rna_lfc)) %>%
              arrange(desc(abs(rna_lfc))) %>%
              head(12)
            if (!is.null(top_b$symbol) && sum(!is.na(top_b$symbol))>0)
              geom_text_repel(data=top_b,
                              aes(x=lfc_botv, y=lfc_botcv, label=symbol),
                              colour="black", size=2.5, fontface="italic",
                              max.overlaps=15, box.padding=0.4,
                              segment.size=0.3, seed=42, inherit.aes=FALSE)
          }
        } +
        labs(title=sprintf("Chromatin bistability × Expression\n%s", rna_cfg$label),
             subtitle=sprintf("%d bistable peaks (|lFC|>%.1f, opposite signs)\nColour = RNA lFC near that peak's gene",
                              nrow(bistable_only), BISTABLE_THRESH),
             x="BOTv temporal log2FC (nc14b→late)",
             y="BOTCv temporal log2FC (nc14b→late)") +
        theme_classic(base_size=11) +
        theme(plot.title=element_text(face="bold",size=10.5))

      # Add quadrant annotation
      q_size <- diff(range(c(all_shared$lfc_botv, all_shared$lfc_botcv),
                            na.rm=TRUE)) * 0.04
      xl <- range(all_shared$lfc_botv, na.rm=TRUE)
      yl <- range(all_shared$lfc_botcv, na.rm=TRUE)
      p_bist_scatter <- p_bist_scatter +
        annotate("text", x=xl[2]*0.8, y=yl[1]*0.8,
                 label="BOTv opens\nBOTCv closes", size=3,
                 colour=COL_BOTV, fontface="bold") +
        annotate("text", x=xl[1]*0.8, y=yl[2]*0.8,
                 label="BOTCv opens\nBOTv closes", size=3,
                 colour=COL_BOTCV, fontface="bold")

      # ── 5e: Concordance class summary bar ────────────────────────────────
      conc_sum <- bistable_only %>%
        filter(!is.na(chromatin_vs_rna)) %>%
        group_by(bistable_class, chromatin_vs_rna) %>%
        summarise(n=n(), n_sig=sum(rna_sig, na.rm=TRUE), .groups="drop")

      p_conc_bar <- ggplot(conc_sum,
                            aes(x=bistable_class, y=n, fill=chromatin_vs_rna)) +
        geom_col(position="dodge", width=0.7) +
        geom_text(aes(label=sprintf("%d\n(%d sig)",n,n_sig)),
                  position=position_dodge(width=0.7),
                  vjust=-0.2, size=2.8) +
        scale_fill_manual(values=c(
          "Concordant (BOTCv-biased)" = COL_BOTCV,
          "Concordant (BOTv-biased)"  = COL_BOTV,
          "Discordant"                = COL_DISC,
          "Co-directional chromatin"  = "#F39C12"),
          name="Chromatin↔RNA") +
        labs(title="Bistable peak concordance with expression",
             subtitle=sprintf("(n sig) = RNA padj<%.2f", FDR_RNA),
             x="Bistable class", y="Peaks with expression data") +
        theme_classic(base_size=10) +
        theme(plot.title=element_text(face="bold"),
              axis.text.x=element_text(angle=15, hjust=1))

      # ── 5f: Top concordant bistable genes ────────────────────────────────
      top_conc <- bistable_only %>%
        filter(grepl("Concordant", chromatin_vs_rna),
               !is.na(symbol), symbol!="",
               !is.na(rna_lfc)) %>%
        mutate(combined_score = abs(rna_lfc) * divergence_score) %>%
        arrange(desc(combined_score)) %>%
        head(TOP_N) %>%
        mutate(direction_label = sprintf("ATAC: %.2f/%.2f | RNA: %.2f",
                                          lfc_botv, lfc_botcv, rna_lfc))

      if (nrow(top_conc)>0) {
        p_top_genes <- ggplot(top_conc,
                               aes(x=reorder(symbol, combined_score),
                                   y=combined_score,
                                   fill=bistable_class)) +
          geom_col(width=0.75) +
          geom_text(aes(label=direction_label), hjust=-0.03, size=2.5) +
          scale_fill_manual(values=c(
            "BOTCv opens / BOTv closes" = COL_BOTCV,
            "BOTv opens / BOTCv closes" = COL_BOTV),
            name="Bistable class") +
          scale_y_continuous(expand=expansion(mult=c(0.02,0.45))) +
          coord_flip(clip="off") +
          labs(title=sprintf("Top concordant bistable genes\n%s", rna_cfg$label),
               subtitle="Combined score = |RNA lFC| × chromatin divergence",
               x=NULL, y="Combined score") +
          theme_classic(base_size=10) +
          theme(plot.title=element_text(face="bold"),
                legend.position="bottom",
                plot.margin=margin(5,160,5,5,"pt"))

        p_full <- (p_bist_scatter | p_conc_bar) / p_top_genes +
          plot_layout(heights=c(1.4,1))
      } else {
        p_full <- p_bist_scatter | p_conc_bar
      }

      stem_out <- file.path(OUT_DIR,"part5_bistable",
                            sprintf("bistable_concordance_%s", rna_cfg$name))
      save_plot(stem_out, p_full, w=14, h=10)

      # ── 5g: Chromatin vs RNA LFC scatter for bistable peaks only ─────────
      if (nrow(bistable_only)>2) {
        # Express chromatin bistability as a single score:
        # positive = BOTCv-biased (BOTCv opens), negative = BOTv-biased
        bistable_only2 <- bistable_only %>%
          mutate(chrom_bias = lfc_botcv - lfc_botv)  # + = BOTCv-biased

        p_bias_scatter <- ggplot(bistable_only2,
                                  aes(x=chrom_bias, y=rna_lfc,
                                      colour=chromatin_vs_rna)) +
          geom_hline(yintercept=0,lty=2,colour="#ccc") +
          geom_vline(xintercept=0,lty=2,colour="#ccc") +
          geom_point(aes(size=divergence_score), alpha=0.8) +
          scale_colour_manual(values=c(
            "Concordant (BOTCv-biased)" = COL_BOTCV,
            "Concordant (BOTv-biased)"  = COL_BOTV,
            "Discordant"                = COL_DISC,
            "Co-directional chromatin"  = "#F39C12"),
            name="Class") +
          scale_size_continuous(range=c(1.5,5), name="Chrom.\ndivergence") +
          {
            if (has_ggrepel && !is.null(bistable_only2$symbol)) {
              top_b2 <- bistable_only2 %>%
                filter(!is.na(symbol), symbol!="") %>%
                arrange(desc(abs(rna_lfc)*divergence_score)) %>% head(12)
              if (nrow(top_b2)>0)
                geom_text_repel(data=top_b2,
                                aes(x=chrom_bias, y=rna_lfc, label=symbol),
                                colour="black", size=2.5, fontface="italic",
                                max.overlaps=15, seed=42, inherit.aes=FALSE)
            }
          } +
          geom_smooth(method="lm", se=TRUE, colour="#333", linewidth=0.6,
                      linetype="dashed", alpha=0.15) +
          labs(title=sprintf("Chromatin bias vs expression\n%s", rna_cfg$label),
               subtitle="Chrom. bias = lfc_BOTCv − lfc_BOTv (positive = BOTCv-biased chromatin)",
               x="Chromatin bias score",
               y=sprintf("RNA log2FC (%s)", rna_cfg$positive_is)) +
          theme_classic(base_size=11) +
          theme(plot.title=element_text(face="bold",size=10.5))

        save_plot(paste0(stem_out,"_bias_scatter"), p_bias_scatter, w=8, h=7)

        # Pearson correlation
        valid_r <- !is.na(bistable_only2$chrom_bias) & !is.na(bistable_only2$rna_lfc)
        if (sum(valid_r)>=5) {
          cr <- cor.test(bistable_only2$chrom_bias[valid_r],
                         bistable_only2$rna_lfc[valid_r])
          cat(sprintf("[v2]   Chromatin bias × RNA lFC correlation: r=%.3f, p=%.4f (n=%d)\n",
                      cr$estimate, cr$p.value, sum(valid_r)))
        }
      }
    }

    # ── 5h: Bistability class overview (all peaks, not just those with RNA) ──
    class_counts <- shared %>%
      count(bistable_class) %>%
      mutate(pct=100*n/sum(n))

    p_class_bar <- ggplot(class_counts,
                           aes(x=reorder(bistable_class,-n), y=n,
                               fill=bistable_class)) +
      geom_col(width=0.7) +
      geom_text(aes(label=sprintf("%d\n(%.0f%%)",n,pct)),
                vjust=-0.2, size=3) +
      scale_fill_manual(values=c(
        "BOTv opens / BOTCv closes"  = COL_BOTV,
        "BOTCv opens / BOTv closes"  = COL_BOTCV,
        "Both open"                  = "#F39C12",
        "Both close"                 = "#8E44AD",
        "Stable"                     = "#AAAAAA",
        "Intermediate"               = "#cccccc"),
        guide="none") +
      labs(title="Chromatin trajectory classes (shared BOTv/BOTCv temporal peaks)",
           subtitle=sprintf("Bistable threshold: |lFC|>%.1f, opposite signs | n=%d total",
                            BISTABLE_THRESH, nrow(shared)),
           x=NULL, y="Number of peaks") +
      theme_classic(base_size=11) +
      theme(plot.title=element_text(face="bold"),
            axis.text.x=element_text(angle=15, hjust=1))

    save_plot(file.path(OUT_DIR,"part5_bistable","bistable_class_summary"),
              p_class_bar, w=9, h=5)
  }
}

cat(sprintf("\n[v2] Done. Outputs in %s/\n", OUT_DIR))
cat("[v2] Key outputs:\n")
cat("[v2]   part4_crm/crm_concordance_*.pdf  — named enhancer × expression\n")
cat("[v2]   part4_crm/crm_multi_dar_summary.tsv — CRMs consistent across DAR sets\n")
cat("[v2]   part5_bistable/bistable_concordance_*.pdf — bistable chromatin × RNA\n")
cat("[v2]   part5_bistable/bistable_*.tsv — full tables for follow-up\n")
