################################################################################
# 00_tf_venn_shared_helpers.R
#
# Shared config + helpers for the Venn-diagram and TF-enrichment scripts
# (Split 1 Non-Ventralized and Split 2 Ventralized).
#
# WHAT'S IN HERE:
#   1. Master genotype color palette (the canonical 16-genotype list)
#   2. A "concept" color palette for non-genotype set names (decomposition
#      categories, background-correction methods, temporal categories, etc.)
#   3. get_color() — smart label -> color resolver (genotype substring match,
#      then concept keyword match, then a stable auto-assigned fallback)
#   4. make_venn2() / make_venn3() — built on VennDiagram::draw.pairwise.venn()/
#      draw.triple.venn() (NOT ggVennDiagram — its fill is a fixed count-based
#      gradient with no per-set color aesthetic, so genotype/concept colors
#      never actually rendered there). Category-name labels are wrapped and
#      pushed outward with margin/cat.dist so they don't clip at the plot
#      edge, with colors auto-resolved from genotype/concept names unless
#      explicitly overridden.
#   5. load_bed_generic() — GRanges loader. Uses a plain, proven-safe
#      read.table parser as the primary path (rtracklayer::import() does
#      C-level parsing that can hard-crash R on malformed input, so it's
#      only an opt-in last resort — see allow_rtracklayer_fallback).
#   6. PARTNER_PEAKS_PATHS — the canonical 17-TF partner peak file list.
#   7. TF "explanatory power" engine — for a DAR set, works out what fraction
#      of DARs are explained by binding of >=1 partner TF, builds pies/bars/
#      combinatorial summaries, and a cross-DAR-set comparison plot.
#
# Source this file AFTER NV_PATHS.r (or your own path config) in each script:
#   source("./NV_PATHS.r")
#   source("./00_tf_venn_shared_helpers.R")
################################################################################

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(scales)
})
# NOTE (v4): switched from ggVennDiagram to VennDiagram for the Venn plots.
# ggVennDiagram's circle fill is driven by a count-based continuous gradient
# and its circle borders are drawn with a fixed color internally — there is
# no per-set "identity" aesthetic for scale_color_manual()/scale_fill_manual()
# to actually attach to, so genotype/concept recoloring silently had no
# visible effect. VennDiagram::draw.pairwise.venn()/draw.triple.venn() take
# explicit fill=/col=/cat.col= arguments per circle, which is the only
# reliable way to get genotype-specific colors onto the diagram.
if (!requireNamespace("VennDiagram", quietly = TRUE)) install.packages("VennDiagram")
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(futile.logger))
flog.threshold(ERROR, name = "VennDiagramLogger")  # silence its internal logger
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)
have_rtracklayer <- requireNamespace("rtracklayer", quietly = TRUE)

################################################################################
# 1. MASTER GENOTYPE COLOR PALETTE
################################################################################

genotype_colors <- c(
  "yw_WT"          = "#555555",
  "BOTv"           = "#2ca02c",
  "BOTv_late"      = "#1a6b1a",
  "BOTv_gastr"     = "#0d3d0d",   # darkest green  -  BOTv gastrulation
  "BOTCv"          = "#e377c2",
  "BOTCv_late"     = "#b5369a",
  "BOT"            = "#1f77b4",
  "BOT_late"       = "#145380",   # singleton (rep1 removed)
  "BOT_hR"         = "#17becf",
  "BOT_hR_gastr"   = "#0e8a8a",
  "BOTR"           = "#9467bd",
  "BOTR_gastr"     = "#5c3585",
  "BOTC"           = "#ff7f0e",
  "BOTC_nc14late"  = "#d45e00",   # darker orange  -  BOTC late (rep2=interim proxy)
  "BOTC_oR"        = "#bcbd22",
  "BOTC_oR_late"   = "#8c8c00"
)

################################################################################
# 2. CONCEPT PALETTE — for set names that aren't single genotypes
#    (decomposition categories, background-correction methods, temporal
#    classes, rescue framing, etc). Chosen to stay visually distinct from
#    the genotype palette above and from each other.
################################################################################

concept_colors <- c(
  "Runt_specific"      = "#1f77b4",
  "Cic_specific"        = "#ff7f0e",
  "Convergent"          = "#2ca02c",
  "Opposing"            = "#d62728",
  "Rescue"              = "#bcbd22",
  "Rescue_cross_ctx"    = "#8B8000",
  "Dosage_gradient"     = "#0E6655",
  "tolrm9_inferred"     = "#95A5A6",
  "tolrm9_empirical"    = "#5D6D7E",
  "Validated_both"      = "#34495E",
  "Stable"              = "#117A65",
  "Emerging"            = "#F39C12",
  "Resolving"           = "#7D3C98",
  "Temporal_divergent"  = "#E74C3C",
  "Background_method"   = "#85929E",
  "WT_step"             = "#1f77b4",
  "Het_step"            = "#17becf",
  "Timepoint_early"     = "#2E86C1",
  "Timepoint_late"      = "#1B4F72",
  "Empirical_validation"= "#C0392B"
)

# Stable fallback palette for anything not matched above (Okabe-Ito-ish,
# colorblind-safe). A registry caches assignments so the same unmatched
# label always gets the same color within one R session.
.fallback_palette <- c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00",
                        "#56B4E9","#F0E442","#999999","#7570B3","#66A61E")
.fallback_registry <- new.env(parent = emptyenv())

#' Resolve a fill/outline color for an arbitrary set label.
#' Order of precedence: exact genotype match (longest substring wins) ->
#' concept keyword match -> stable auto-assigned fallback color.
get_color <- function(label) {
  if (is.null(label) || is.na(label) || !nzchar(label)) return("grey50")

  # --- genotype substring match, longest key first so "BOTC_oR" beats "BOTC"
  #     beats "BOT" ---
  geno_keys <- names(genotype_colors)[order(-nchar(names(genotype_colors)))]
  for (k in geno_keys) {
    if (str_detect(label, fixed(k))) return(unname(genotype_colors[k]))
  }

  # --- concept keyword match (case-insensitive, loose) ---
  concept_map <- list(
    Runt_specific     = c("runt.specific", "runt_specific"),
    Cic_specific      = c("cic.specific", "cic_specific"),
    Convergent        = c("convergent"),
    Opposing          = c("opposing"),
    Rescue_cross_ctx  = c("rescue.cross", "cross.ctx", "cross_ctx", "cross-cic"),
    Rescue            = c("rescue"),
    Dosage_gradient   = c("dosage", "gradient"),
    tolrm9_inferred   = c("inferred"),
    tolrm9_empirical  = c("empirical"),
    Validated_both    = c("validated"),
    Stable            = c("stable"),
    Emerging          = c("emerging"),
    Resolving         = c("resolving"),
    Temporal_divergent= c("temporal.divergent", "divergent"),
    Het_step          = c("het"),
    WT_step           = c("wt.?step", "wt.?stage")
  )
  label_l <- tolower(label)
  for (nm in names(concept_map)) {
    pats <- concept_map[[nm]]
    if (any(str_detect(label_l, pats))) return(unname(concept_colors[nm]))
  }

  # --- stable fallback ---
  if (!exists(label, envir = .fallback_registry, inherits = FALSE)) {
    n_used <- length(ls(envir = .fallback_registry))
    col <- .fallback_palette[(n_used %% length(.fallback_palette)) + 1]
    assign(label, col, envir = .fallback_registry)
  }
  get(label, envir = .fallback_registry, inherits = FALSE)
}

################################################################################
# 3. LABEL WRAPPING HELPER — prevents long category names from running off
#    the edge of Venn diagrams / axis text.
################################################################################

wrap_label <- function(x, width = 18) str_wrap(x, width = width)

# Null-coalesce helper (defined early since several plotting helpers below
# use it for default titles).
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && !nzchar(a))) b else a

################################################################################
# 4. RECOLORED, LABEL-SAFE VENN HELPERS
#
#    FIX FOR CLIPPED LABELS:
#      ggVennDiagram places category names just outside the circles using
#      plot-area coordinates with zero default expansion, so long labels
#      (or labels close to the panel edge) get clipped by the device. We
#      fix this by (a) wrapping long labels onto multiple lines, (b)
#      expanding the x/y scales so there's breathing room outside the
#      circles, (c) capping label/set text size so text doesn't grow wider
#      than the available margin, and (d) widening the saved PDF based on
#      how long the labels are so nothing gets clipped on export either.
################################################################################

#' 2-set Venn diagram with REAL per-circle colors + non-clipped labels.
#' Uses VennDiagram::draw.pairwise.venn(), which takes explicit fill=/col=
#' per circle — unlike ggVennDiagram, whose fill is a fixed count-based
#' gradient with no per-set color aesthetic to hook into.
#' @param g1,g2 GRanges sets
#' @param lab1,lab2 set labels (will be wrapped automatically)
#' @param title_str plot title
#' @param fname output filename (saved under output_dir)
#' @param col1,col2 optional explicit colors; if NULL, resolved via get_color()
#' @param label_type one of "count"|"percent"|"both"|"none"
#' @param wrap_width characters per line before wrapping a label
#' @param out_dir output directory (defaults to output_dir in caller env)
make_venn2 <- function(g1, g2, lab1, lab2, title_str, fname,
                        col1 = NULL, col2 = NULL,
                        label_type = "both", wrap_width = 16,
                        out_dir = output_dir) {
  if (length(g1) == 0 && length(g2) == 0) return(invisible(NULL))
  # FIX: use genomic overlap (any bp overlap), not exact chrom/start/end string
  # matching. peak_key()-based exact matching under-counts "shared" whenever
  # two independently-called DAR sets have slightly different peak boundaries
  # at the same locus -- which is the norm, not the exception, across
  # different DESeq2/limma runs. subsetByOverlaps(minoverlap=1L) is the same
  # definition already used by overlap_count() in the calling scripts and by
  # the Fisher enrichment test below, so Venn numbers now agree with the
  # overlap summary tables and enrichment output instead of silently using a
  # stricter, inconsistent definition of "shared".
  shared <- length(subsetByOverlaps(g1, g2, minoverlap = 1L))
  cat(sprintf("  %-35s  %d | shared: %d | %-35s  %d\n",
              lab1, length(g1), shared, lab2, length(g2)))

  col1 <- if (is.null(col1)) get_color(lab1) else col1
  col2 <- if (is.null(col2)) get_color(lab2) else col2

  lab1_w <- wrap_label(lab1, wrap_width)
  lab2_w <- wrap_label(lab2, wrap_width)

  print_mode <- switch(label_type,
                        "both"    = c("raw", "percent"),
                        "count"   = "raw",
                        "percent" = "percent",
                        "none"    = character(0),
                        c("raw", "percent"))

  # FIX: generous margin + cat.dist push the (possibly multi-line, wrapped)
  # category labels outward so they don't get clipped at the plot edge;
  # canvas width also scales with the longest wrapped label.
  max_chars <- max(nchar(lab1_w), nchar(lab2_w), na.rm = TRUE)
  out_w <- max(6.5, 6.5 + 0.05 * max_chars)
  out_h <- 5.5

  grDevices::pdf(file.path(out_dir, fname), width = out_w, height = out_h)
  grid::grid.newpage()
  invisible(draw.pairwise.venn(
    area1 = length(g1), area2 = length(g2), cross.area = shared,
    category   = c(lab1_w, lab2_w),
    fill       = c(col1, col2),
    col        = c(col1, col2),
    alpha      = 0.55, lwd = 2,
    cat.col    = c(col1, col2), cat.cex = 1.0, cat.fontface = "bold",
    cex        = 1.0, fontface = "plain",
    print.mode = print_mode,
    cat.dist   = c(0.045, 0.045), cat.pos = c(-30, 30),
    margin     = 0.12,
    ind        = TRUE
  ))
  grid::grid.text(title_str, x = 0.5, y = 0.97,
                   gp = grid::gpar(fontsize = 12, fontface = "bold"))
  grid::grid.text(sprintf("%s: %d  |  shared: %d  |  %s: %d",
                          lab1, length(g1), shared, lab2, length(g2)),
                   x = 0.5, y = 0.925,
                   gp = grid::gpar(fontsize = 9, col = "grey30"))
  grDevices::dev.off()
  invisible(NULL)
}

#' 3-set Venn diagram with REAL per-circle colors + non-clipped labels.
#' Uses VennDiagram::draw.triple.venn() for the same reason as make_venn2().
make_venn3 <- function(g1, g2, g3, lab1, lab2, lab3, title_str, fname,
                        col1 = NULL, col2 = NULL, col3 = NULL,
                        label_type = "both", wrap_width = 16,
                        out_dir = output_dir) {
  if (length(g1) == 0 || length(g2) == 0 || length(g3) == 0) return(invisible(NULL))

  col1 <- if (is.null(col1)) get_color(lab1) else col1
  col2 <- if (is.null(col2)) get_color(lab2) else col2
  col3 <- if (is.null(col3)) get_color(lab3) else col3

  lab1_w <- wrap_label(lab1, wrap_width)
  lab2_w <- wrap_label(lab2, wrap_width)
  lab3_w <- wrap_label(lab3, wrap_width)

  print_mode <- switch(label_type,
                        "both"    = c("raw", "percent"),
                        "count"   = "raw",
                        "percent" = "percent",
                        "none"    = character(0),
                        c("raw", "percent"))

  # FIX: genomic overlap (subsetByOverlaps, minoverlap=1L) instead of exact
  # chrom/start/end string matching -- see make_venn2() for rationale. n123
  # (all three) is computed directionally from g1: peaks in g1 that also
  # overlap g2 AND overlap g3. This matches the convention draw.triple.venn()
  # expects (n123 must be <= min(n12, n13)) and keeps the triple-overlap
  # count consistent with how n12/n13 are themselves defined from g1's side.
  n12  <- length(subsetByOverlaps(g1, g2, minoverlap = 1L))
  n23  <- length(subsetByOverlaps(g2, g3, minoverlap = 1L))
  n13  <- length(subsetByOverlaps(g1, g3, minoverlap = 1L))
  n123 <- length(subsetByOverlaps(subsetByOverlaps(g1, g2, minoverlap = 1L),
                                   g3, minoverlap = 1L))

  max_chars <- max(nchar(lab1_w), nchar(lab2_w), nchar(lab3_w), na.rm = TRUE)
  out_w <- max(7.5, 7.5 + 0.04 * max_chars)
  out_h <- 6.8

  grDevices::pdf(file.path(out_dir, fname), width = out_w, height = out_h)
  grid::grid.newpage()
  invisible(draw.triple.venn(
    area1 = length(g1), area2 = length(g2), area3 = length(g3),
    n12 = n12, n23 = n23, n13 = n13, n123 = n123,
    category   = c(lab1_w, lab2_w, lab3_w),
    fill       = c(col1, col2, col3),
    col        = c(col1, col2, col3),
    alpha      = 0.55, lwd = 2,
    cat.col    = c(col1, col2, col3), cat.cex = 1.0, cat.fontface = "bold",
    cex        = 1.0,
    print.mode = print_mode,
    cat.dist   = c(0.07, 0.07, 0.04),
    margin     = 0.12,
    ind        = TRUE
  ))
  grid::grid.text(title_str, x = 0.5, y = 0.97,
                   gp = grid::gpar(fontsize = 12, fontface = "bold"))
  grDevices::dev.off()
  invisible(NULL)
}

peak_key <- function(gr) paste(as.character(seqnames(gr)), start(gr), end(gr), sep = "_")

################################################################################
# 5. GENERIC GRANGES LOADER — prefers rtracklayer::import(), falls back to
#    read.table for plain 3-column BED if rtracklayer isn't installed.
################################################################################

#' Load a BED-like file into a GRanges, using only coordinates (chrom/start/
#' end[/name]) — nothing downstream (overlap counts, Fisher tests, Venn
#' membership) needs strand or score, so a plain, fast, proven-safe
#' read.table parser is the PRIMARY path. rtracklayer::import() does
#' compiled-code (C-level) parsing and can hard-crash the whole R session
#' (not just throw a catchable R error) on a file with any format quirk —
#' a header/track line, an odd strand value, a score outside 0-1000, etc.
#' tryCatch() cannot protect against that kind of crash, so rtracklayer is
#' now only attempted as a LAST-RESORT fallback if the plain parser fails
#' outright, and only when the user has explicitly opted in (see
#' `allow_rtracklayer_fallback` below).
allow_rtracklayer_fallback <- FALSE  # set TRUE only if you've confirmed your
                                      # BED files are clean and need import()'s
                                      # richer parsing (e.g. true narrowPeak)

load_bed_generic <- function(path) {
  if (!file.exists(path)) { message("  [missing] ", path); return(GRanges()) }

  df <- tryCatch({
    # Skip blank/comment/track-definition lines some BED exports include;
    # fill=TRUE tolerates a ragged trailing column without erroring.
    raw <- read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                       comment.char = "#", fill = TRUE,
                       blank.lines.skip = TRUE)
    raw <- raw[!grepl("^track", raw[, 1], ignore.case = TRUE), , drop = FALSE]
    raw
  }, error = function(e) NULL)

  if (!is.null(df) && nrow(df) > 0 && ncol(df) >= 3) {
    starts <- suppressWarnings(as.numeric(df[, 2]))
    ends   <- suppressWarnings(as.numeric(df[, 3]))
    keep   <- !is.na(starts) & !is.na(ends)
    if (any(!keep)) message("  [", path, "] dropped ", sum(!keep),
                            " row(s) with non-numeric coordinates")
    if (sum(keep) > 0) {
      gr <- GRanges(seqnames = df[keep, 1],
                    ranges   = IRanges(start = starts[keep] + 1, end = ends[keep]))
      if (ncol(df) >= 4) names(gr) <- df[keep, 4]
      return(gr)
    }
  }

  if (allow_rtracklayer_fallback && have_rtracklayer) {
    message("  [", path, "] plain parser failed  -  trying rtracklayer::import()",
            " (allow_rtracklayer_fallback=TRUE)")
    gr <- tryCatch(rtracklayer::import(path), error = function(e) NULL)
    if (!is.null(gr)) return(gr)
  }

  message("  [unreadable] ", path)
  GRanges()
}

################################################################################
# 6. PARTNER TF PEAK SET — the 17 TFs of interest for the explanatory
#    breakdown / enrichment analyses. Paths are relative to the working
#    directory each script is run from (same convention as the rest of the
#    pipeline). Loaded lazily via load_partner_peaks() so a script can run
#    even if some bed files aren't present on a given machine.
################################################################################

PARTNER_PEAKS_PATHS <- list(
  Bcd_FDR1                  = "./data/tf_partners/bdtnp_Bcd1-2_FDR1_all_sites_dm6.bed",
  Ftz_FDR1                  = "./data/tf_partners/bdtnp_Ftz3_FDR1_dm6.bed",
  Ftz_FDR25                 = "./data/tf_partners/bdtnp_Ftz3_FDR25_dm6.bed",
  Cic                       = "./data/tf_partners/Cic_ChIPseq_published_union_annotated.bed",
  Hb_ChIPseq                = "./data/tf_partners/ChIP_Hb_GSE50771_dm6.bed",
  Gt_ChIPseq                = "./data/tf_partners/ChIP_Gt_GSE50771_dm6.bed",
  Prd_FDR1                  = "./data/tf_partners/bdtnp_Prd_FDR1_dm6.bed",
  Tll_FDR1                  = "./data/tf_partners/bdtnp_Tll_FDR1_dm6.bed",
  Zld_2hr                   = "./data/tf_partners/ChIP_Zld_2hr_GSM763061_peaks_dm6.bed",
  Hkb_FDR1                  = "./data/tf_partners/bdtnp_Hkb1-3_FDR1_all_sites_dm6.bed",
  D_FDR1                    = "./data/tf_partners/bdtnp_D_FDR1_dm6.bed",
  Cad_FDR1                  = "./data/tf_partners/bdtnp_Cad_FDR1_dm6.bed",
  Hairy_FDR1                = "./data/tf_partners/bdtnp_Hairy1-2_FDR1_all_sites_dm6.bed",
  Opa_early                 = "./data/tf_partners/ChIP_Opa_Early_rep1-2_q01_overlapping_sites.bed",
  Opa_late                  = "./data/tf_partners/ChIP_Opa_Late_rep1-2_q01_overlapping_sites.bed",
  Run_ChIP_chip_FDR1        = "./data/chipchip/bdtnp_Run1-2_FDR1_all_sites.bed",
  Runt_ChIPseq_Reproducible = "./data/chipseq/RuntAb_IgG_p0.05_reproducible.bed"
)

#' Load every partner TF peak file into a named list of GRanges, skipping
#' (with a message) any that are missing on this machine. Output is flushed
#' after every file (flush.console()) so if R crashes mid-load, whatever
#' printed last in the console tells you exactly which file it died on.
load_partner_peaks <- function(paths = PARTNER_PEAKS_PATHS) {
  out <- list()
  for (nm in names(paths)) {
    cat(sprintf("  loading %-28s ...", nm)); flush.console()
    gr <- load_bed_generic(paths[[nm]])
    if (length(gr) > 0) {
      out[[nm]] <- gr
      cat(sprintf(" %d peaks\n", length(gr)))
    } else {
      cat(" skipped (0 peaks)\n")
    }
    flush.console()
  }
  out
}

# A dedicated, stable color per partner TF so the same TF is always the same
# color across every plot in every script (pies, bars, heatmaps).
.tf_palette_base <- c(
  "#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2",
  "#7f7f7f","#bcbd22","#17becf","#aec7e8","#ffbb78","#98df8a","#ff9896",
  "#c5b0d5","#c49c94","#f7b6d2"
)
tf_colors <- setNames(.tf_palette_base[seq_along(PARTNER_PEAKS_PATHS)],
                       names(PARTNER_PEAKS_PATHS))

################################################################################
# 7. TF "EXPLANATORY POWER" ENGINE
#
#    For a given DAR set, work out what fraction of DARs can be "explained"
#    by binding of at least one partner TF, how much of that explanatory
#    power is attributable to single vs. combinatorial TF binding, and which
#    individual TFs carry the most weight. Produces pies, bars, and a
#    cross-DAR-set comparison plot.
################################################################################

#' Build a logical DAR x TF overlap matrix.
tf_overlap_matrix <- function(dar_gr, tf_grs) {
  if (length(dar_gr) == 0 || length(tf_grs) == 0) {
    return(matrix(logical(0), nrow = 0, ncol = length(tf_grs),
                   dimnames = list(NULL, names(tf_grs))))
  }
  m <- vapply(tf_grs, function(g) {
    if (length(g) == 0) return(rep(FALSE, length(dar_gr)))
    countOverlaps(dar_gr, g, minoverlap = 1L) > 0
  }, FUN.VALUE = logical(length(dar_gr)))
  if (is.null(dim(m))) m <- matrix(m, nrow = length(dar_gr))
  colnames(m) <- names(tf_grs)
  m
}

#' Summarize explanatory power of a TF panel over one DAR set.
#' Returns a list with: overlap matrix, per-TF counts/percentages
#' (non-exclusive — a DAR bound by 2 TFs counts toward both), the
#' single/multi/unexplained breakdown, TF-bound-count distribution, and the
#' top combinatorial binding patterns.
tf_explanatory_summary <- function(dar_gr, tf_grs, label = "") {
  n_total <- length(dar_gr)
  m <- tf_overlap_matrix(dar_gr, tf_grs)
  if (n_total == 0 || ncol(m) == 0) {
    return(list(label = label, n_total = n_total, per_tf = NULL,
                n_explained = 0, n_unexplained = n_total,
                n_single = 0, n_multi = 0, count_dist = NULL,
                combo_table = NULL, matrix = m))
  }

  n_tf_per_dar <- rowSums(m)
  n_explained   <- sum(n_tf_per_dar > 0)
  n_unexplained <- n_total - n_explained
  n_single      <- sum(n_tf_per_dar == 1)
  n_multi       <- sum(n_tf_per_dar >= 2)

  per_tf <- data.frame(
    tf            = colnames(m),
    n_bound       = colSums(m),
    pct_of_total  = round(100 * colSums(m) / max(n_total, 1), 1),
    stringsAsFactors = FALSE
  ) %>% arrange(desc(n_bound))

  # which TF "wins" for DARs bound by exactly one TF
  single_tf_id <- rep(NA_character_, n_total)
  if (n_single > 0) {
    single_idx <- which(n_tf_per_dar == 1)
    single_tf_id[single_idx] <- colnames(m)[apply(m[single_idx, , drop = FALSE], 1, which.max)]
  }
  single_comp <- if (n_single > 0) {
    as.data.frame(table(single_tf_id), stringsAsFactors = FALSE) %>%
      rename(tf = single_tf_id, n = Freq) %>%
      mutate(pct_of_single = round(100 * n / n_single, 1)) %>%
      arrange(desc(n))
  } else NULL

  count_dist <- data.frame(n_tf_bound = n_tf_per_dar) %>%
    mutate(bucket = case_when(
      n_tf_bound == 0 ~ "0",
      n_tf_bound == 1 ~ "1",
      n_tf_bound == 2 ~ "2",
      n_tf_bound == 3 ~ "3",
      TRUE             ~ "4+"
    )) %>%
    count(bucket, name = "n") %>%
    mutate(bucket = factor(bucket, levels = c("0","1","2","3","4+")),
           pct = round(100 * n / n_total, 1)) %>%
    arrange(bucket)

  # top combinatorial binding patterns among explained DARs (collapse each
  # DAR's TRUE columns into a "TF_A+TF_B" string)
  combo_table <- NULL
  if (n_explained > 0) {
    explained_idx <- which(n_tf_per_dar > 0)
    combo_strs <- apply(m[explained_idx, , drop = FALSE], 1, function(r) {
      paste(sort(colnames(m)[r]), collapse = " + ")
    })
    combo_table <- as.data.frame(table(combo_strs), stringsAsFactors = FALSE) %>%
      rename(combo = combo_strs, n = Freq) %>%
      mutate(pct_of_explained = round(100 * n / n_explained, 1)) %>%
      arrange(desc(n))
  }

  list(label = label, n_total = n_total, per_tf = per_tf,
       n_explained = n_explained, n_unexplained = n_unexplained,
       n_single = n_single, n_multi = n_multi,
       single_composition = single_comp, count_dist = count_dist,
       combo_table = combo_table, matrix = m)
}

# --- shared donut/pie theme -------------------------------------------------
.pie_theme <- theme_void(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey30"),
        legend.title  = element_text(size = 9),
        legend.text   = element_text(size = 8.5),
        plot.margin   = margin(10, 20, 10, 20))

#' Donut: Unexplained vs Single-TF-explained vs Multi-TF-explained.
plot_explanatory_pie <- function(summ, title_str = NULL) {
  if (is.null(summ$per_tf) || summ$n_total == 0) return(NULL)
  df <- data.frame(
    category = factor(c("Unexplained","Explained\n(1 TF)","Explained\n(>=2 TFs)"),
                       levels = c("Unexplained","Explained\n(1 TF)","Explained\n(>=2 TFs)")),
    n = c(summ$n_unexplained, summ$n_single, summ$n_multi)
  ) %>% mutate(pct = round(100 * n / sum(n), 1),
               lab = sprintf("%s\n%d (%.1f%%)", category, n, pct))

  ggplot(df, aes(x = 2, y = n, fill = category)) +
    geom_col(width = 1, color = "white", linewidth = 0.8) +
    coord_polar(theta = "y") +
    xlim(0.4, 2.5) +
    geom_text(aes(label = ifelse(n > 0, lab, "")), position = position_stack(vjust = 0.5),
               size = 3.1, lineheight = 0.9) +
    scale_fill_manual(values = c("Unexplained" = "grey80",
                                  "Explained\n(1 TF)" = "#4C72B0",
                                  "Explained\n(>=2 TFs)" = "#DD8452"),
                       name = NULL) +
    labs(title = title_str %||% summ$label,
         subtitle = sprintf("n = %d DARs total", summ$n_total)) +
    .pie_theme +
    theme(legend.position = "none")
}

#' Pie: composition of single-TF-explained DARs, by which TF.
plot_single_tf_composition_pie <- function(summ, title_str = NULL, palette = tf_colors) {
  sc <- summ$single_composition
  if (is.null(sc) || nrow(sc) == 0) return(NULL)
  sc <- sc %>% mutate(tf = factor(tf, levels = tf[order(-n)]))
  cols <- palette[as.character(levels(sc$tf))]
  cols[is.na(cols)] <- "grey60"

  ggplot(sc, aes(x = "", y = n, fill = tf)) +
    geom_col(width = 1, color = "white", linewidth = 0.6) +
    coord_polar(theta = "y") +
    geom_text(aes(label = ifelse(pct_of_single >= 4,
                                  sprintf("%s\n%.0f%%", tf, pct_of_single), "")),
              position = position_stack(vjust = 0.5), size = 2.9, lineheight = 0.85) +
    scale_fill_manual(values = cols, name = "TF") +
    labs(title = title_str %||% sprintf("%s  -  single-TF composition", summ$label),
         subtitle = sprintf("Among %d DARs explained by exactly 1 TF", summ$n_single)) +
    .pie_theme
}

#' Horizontal bar: % of total DARs bound by each TF (non-exclusive — TFs can
#' co-occur, so bars need not sum to 100%).
plot_tf_binding_bar <- function(summ, title_str = NULL, palette = tf_colors) {
  pt <- summ$per_tf
  if (is.null(pt) || nrow(pt) == 0) return(NULL)
  pt <- pt %>% filter(n_bound > 0) %>% mutate(tf = factor(tf, levels = rev(tf)))
  cols <- palette[as.character(pt$tf)]
  cols[is.na(cols)] <- "grey60"

  ggplot(pt, aes(x = tf, y = pct_of_total, fill = tf)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%  (n=%d)", pct_of_total, n_bound)),
               hjust = -0.05, size = 3.0) +
    coord_flip(clip = "off") +
    scale_fill_manual(values = cols, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.30))) +
    labs(title = title_str %||% sprintf("%s  -  TF binding coverage", summ$label),
         subtitle = sprintf("%% of %d DARs overlapping each TF (TFs can co-occur)",
                             summ$n_total),
         x = NULL, y = "% of DARs bound") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          plot.title    = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, color = "grey30"),
          plot.margin   = margin(10, 40, 10, 10))
}

#' Donut: distribution of how many TFs bind each DAR (0,1,2,3,4+).
plot_tf_count_distribution <- function(summ, title_str = NULL) {
  cd <- summ$count_dist
  if (is.null(cd) || nrow(cd) == 0) return(NULL)
  pal <- c("0"="grey80","1"="#4C72B0","2"="#55A868","3"="#C44E52","4+"="#8172B2")

  ggplot(cd, aes(x = 2, y = n, fill = bucket)) +
    geom_col(width = 1, color = "white", linewidth = 0.8) +
    coord_polar(theta = "y") +
    xlim(0.4, 2.5) +
    geom_text(aes(label = ifelse(n > 0, sprintf("%s TF%s\n%.1f%%", bucket,
                                                  ifelse(bucket=="1","","s"), pct), "")),
              position = position_stack(vjust = 0.5), size = 3.0, lineheight = 0.9) +
    scale_fill_manual(values = pal, name = "# TFs bound") +
    labs(title = title_str %||% sprintf("%s  -  # TFs bound per DAR", summ$label),
         subtitle = sprintf("n = %d DARs total", summ$n_total)) +
    .pie_theme
}

#' Bar of the top-N combinatorial TF binding patterns among explained DARs.
plot_top_combinations <- function(summ, top_n = 10, title_str = NULL) {
  ct <- summ$combo_table
  if (is.null(ct) || nrow(ct) == 0) return(NULL)
  ct <- head(ct, top_n) %>% mutate(combo = factor(combo, levels = rev(combo)))

  ggplot(ct, aes(x = combo, y = n)) +
    geom_col(fill = "#4C72B0", width = 0.7) +
    geom_text(aes(label = sprintf("%d (%.1f%%)", n, pct_of_explained)),
               hjust = -0.05, size = 2.9) +
    coord_flip(clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.30))) +
    labs(title = title_str %||% sprintf("%s  -  top TF combinations", summ$label),
         subtitle = sprintf("Top %d binding patterns among %d explained DARs",
                             min(top_n, nrow(summ$combo_table)), summ$n_explained),
         x = NULL, y = "# DARs") +
    theme_minimal(base_size = 10.5) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          axis.text.y   = element_text(size = 8),
          plot.title    = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, color = "grey30"),
          plot.margin   = margin(10, 50, 10, 10))
}

#' Orchestrate the full explanatory breakdown for one DAR set: compute the
#' summary, draw all four plot types, save a combined PDF panel + per-plot
#' files, and return a one-row data.frame for cross-set comparison.
run_explanatory_breakdown <- function(dar_gr, tf_grs, label, out_dir,
                                       file_stub = NULL, top_n_combos = 10) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  file_stub <- file_stub %||% gsub("[^A-Za-z0-9]+", "_", label)

  summ <- tf_explanatory_summary(dar_gr, tf_grs, label = label)
  if (summ$n_total == 0) {
    message("  [skip explanatory breakdown  -  0 DARs] ", label)
    return(NULL)
  }

  p_pie    <- plot_explanatory_pie(summ)
  p_single <- plot_single_tf_composition_pie(summ)
  p_bar    <- plot_tf_binding_bar(summ)
  p_dist   <- plot_tf_count_distribution(summ)
  p_combo  <- plot_top_combinations(summ, top_n = top_n_combos)

  # Combined panel (skip any NULL plots gracefully)
  panel_plots <- Filter(Negate(is.null), list(p_pie, p_dist, p_single, p_bar))
  if (length(panel_plots) > 0) {
    combined <- wrap_plots(panel_plots, ncol = 2) +
      plot_annotation(title = sprintf("TF explanatory breakdown  -  %s", label),
                       theme = theme(plot.title = element_text(face = "bold", size = 14,
                                                                hjust = 0.5)))
    ggsave(file.path(out_dir, sprintf("Explanatory_%s_panel.pdf", file_stub)),
           combined, width = 12, height = 10)
  }
  if (!is.null(p_combo)) {
    ggsave(file.path(out_dir, sprintf("Explanatory_%s_top_combinations.pdf", file_stub)),
           p_combo, width = 8, height = max(3.5, 0.4 * top_n_combos + 1.5))
  }

  if (!is.null(summ$per_tf)) {
    write.table(summ$per_tf,
                file.path(out_dir, sprintf("Explanatory_%s_per_TF.txt", file_stub)),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  if (!is.null(summ$combo_table)) {
    write.table(summ$combo_table,
                file.path(out_dir, sprintf("Explanatory_%s_combinations.txt", file_stub)),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  data.frame(
    label         = label,
    n_total       = summ$n_total,
    n_explained   = summ$n_explained,
    n_unexplained = summ$n_unexplained,
    n_single      = summ$n_single,
    n_multi       = summ$n_multi,
    pct_explained = round(100 * summ$n_explained / summ$n_total, 1)
  )
}

#' Cross-DAR-set comparison: stacked bar of Unexplained / Single-TF /
#' Multi-TF percentages across every DAR set tested, sorted by % explained.
plot_explained_summary_across_sets <- function(summary_df, out_dir,
                                                fname = "Explanatory_summary_across_sets.pdf",
                                                title_str = "TF explanatory power across DAR sets") {
  if (is.null(summary_df) || nrow(summary_df) == 0) return(NULL)
  df <- summary_df %>%
    mutate(pct_single = round(100 * n_single / n_total, 1),
           pct_multi  = round(100 * n_multi  / n_total, 1),
           pct_none   = round(100 * n_unexplained / n_total, 1)) %>%
    select(label, pct_none, pct_single, pct_multi) %>%
    pivot_longer(cols = c(pct_none, pct_single, pct_multi),
                 names_to = "category", values_to = "pct") %>%
    mutate(category = case_when(
             category == "pct_none"   ~ "Unexplained",
             category == "pct_single" ~ "Explained (1 TF)",
             TRUE                      ~ "Explained (>=2 TFs)"
           ),
           category = factor(category, levels = c("Unexplained","Explained (1 TF)",
                                                    "Explained (>=2 TFs)")))
  order_lv <- summary_df %>% arrange(desc(pct_explained)) %>% pull(label)
  df$label <- factor(df$label, levels = order_lv)

  p <- ggplot(df, aes(x = label, y = pct, fill = category)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c("Unexplained" = "grey80",
                                  "Explained (1 TF)" = "#4C72B0",
                                  "Explained (>=2 TFs)" = "#DD8452"),
                       name = NULL) +
    scale_x_discrete(labels = setNames(wrap_label(levels(df$label), 18),
                                        levels(df$label))) +
    scale_y_continuous(labels = label_percent(scale = 1),
                        expand = expansion(mult = c(0, 0.03))) +
    labs(title = title_str,
         subtitle = "% of DARs in each set explained by >=1 partner TF",
         x = NULL, y = "% of DARs") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          panel.grid.major.x = element_blank(),
          plot.title    = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 9.5, color = "grey30"))

  out_w <- max(8, length(levels(df$label)) * 0.9)
  ggsave(file.path(out_dir, fname), p, width = out_w, height = 6)
  invisible(p)
}

cat("[00_tf_venn_shared_helpers.R] loaded: genotype + concept palettes, recolored Venn",
    "helpers, and the TF explanatory-breakdown engine are ready.\n")
